#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.reads        = params.reads        ?: 'data/*.fastq.gz'
params.out_dir      = params.out_dir      ?: 'results'

// Resolve bundled reference inside the workflow package
def BUNDLED_REF = file("${baseDir}/assets/reference.fasta")

// Choose reference: use user-provided --ref only if it exists; otherwise fall back to bundled
def ref_path = (params.ref && file(params.ref).exists()) ? params.ref : BUNDLED_REF.toString()
log.info "Reference FASTA: ${ref_path}"
// Resolve bundled reference inside the workflow package


params.threads      = params.threads      ?: 4
params.pseudocount  = params.pseudocount  ?: 1
params.water_id     = params.water_id     ?: 'water_control'
params.pos_cutoff   = params.pos_cutoff   ?: 10
params.inc_cutoff   = params.inc_cutoff   ?: 3

// micromamba runtime (advanced)
params.env_name     = params.env_name     ?: 'wf_viral_screen_env'
params.mamba_prefix = params.mamba_prefix ?: "${System.properties['user.home']}/.micromamba"

workflow {

  new File(params.out_dir).mkdirs()

  // Accept either a folder or a glob
  def reads_glob = params.reads
  def reads_path = file(params.reads)
  if( reads_path.exists() && reads_path.isDirectory() ) {
    reads_glob = "${params.reads}/*.{fastq,fastq.gz}"
  }

  samples_ch = Channel
    .fromPath(reads_glob, checkIfExists: true)
    .map { f ->
      def sample_id = f.getName().replaceFirst(/\.fastq(\.gz)?$/, '')
      tuple(sample_id, f)
    }

  ref_ch = Channel.fromPath(ref_path, checkIfExists: true).first()

  envfile_ch = setup_micromamba()

  bam_ch   = map_reads_minimap2(samples_ch, ref_ch, envfile_ch)
  count_ch = count_reads_per_virus(bam_ch, envfile_ch)

  water_counts_ch = count_ch
    .filter { sid, counts_tsv, bamstats_txt -> sid == params.water_id }
    .map    { sid, counts_tsv, bamstats_txt -> counts_tsv }
    .first()

  sample_countfiles_ch = count_ch
    .filter { sid, counts_tsv, bamstats_txt -> sid != params.water_id }
    .map    { sid, counts_tsv, bamstats_txt -> counts_tsv }
    .collect()

  compute_simple_ratios(sample_countfiles_ch, water_counts_ch, ref_ch)
}

process setup_micromamba {
  tag "micromamba-env"

  output:
    path "mm_env_prefix.txt"

  script:
  """
  set -euo pipefail

  MM="${workflow.projectDir}/bin/micromamba"

  PREFIX="${params.mamba_prefix}"
  if [[ "\$PREFIX" == "~"* ]]; then
    PREFIX="\$HOME\${PREFIX:1}"
  fi
  if [[ "\$PREFIX" != /* ]]; then
    PREFIX="\$PWD/\$PREFIX"
  fi
  mkdir -p "\$PREFIX"

  ENV_NAME="${params.env_name}"
  ENV_PREFIX="\$PREFIX/envs/\$ENV_NAME"

  if [[ ! -x "\$MM" ]]; then
    echo "ERROR: micromamba not found or not executable at: \$MM" >&2
    exit 1
  fi

  if [[ ! -d "\$ENV_PREFIX" ]]; then
    "\$MM" create -y -p "\$ENV_PREFIX" -c conda-forge -c bioconda minimap2=2.28 samtools=1.20
  fi

  echo "\$ENV_PREFIX" > mm_env_prefix.txt
  echo "Micromamba env prefix is: \$(cat mm_env_prefix.txt)" >&2
  """
}

process map_reads_minimap2 {
  tag "$sample_id"
  cpus params.threads

  publishDir "${params.out_dir}/bam", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(fastq_gz)
    path ref_fa
    path envfile

  output:
    tuple val(sample_id),
         path("${sample_id}.sorted.bam"),
         path("${sample_id}.sorted.bam.bai")

  script:
  """
  set -euo pipefail

  MM="${workflow.projectDir}/bin/micromamba"
  ENV_PREFIX=\$(cat "${envfile}")

  if [[ ! -d "\$ENV_PREFIX" ]]; then
    echo "ERROR: micromamba env prefix does not exist: \$ENV_PREFIX" >&2
    echo "envfile content:" >&2
    cat "${envfile}" >&2 || true
    exit 1
  fi

  "\$MM" run -p "\$ENV_PREFIX" minimap2 -ax map-ont --secondary=no -t ${task.cpus} ${ref_fa} ${fastq_gz} \\
    | "\$MM" run -p "\$ENV_PREFIX" samtools view -b - \\
    | "\$MM" run -p "\$ENV_PREFIX" samtools sort -o ${sample_id}.sorted.bam

  "\$MM" run -p "\$ENV_PREFIX" samtools index ${sample_id}.sorted.bam
  """
}

process count_reads_per_virus {
  tag "$sample_id"

  publishDir "${params.out_dir}/counts", mode: 'copy', overwrite: true

  input:
    tuple val(sample_id), path(bam), path(bai)
    path envfile

  output:
    tuple val(sample_id),
         path("${sample_id}.counts.tsv"),
         path("${sample_id}.bamstats.txt")

  script:
  """
  set -euo pipefail

  MM="${workflow.projectDir}/bin/micromamba"
  ENV_PREFIX=\$(cat "${envfile}")

  "\$MM" run -p "\$ENV_PREFIX" samtools stats ${bam} > ${sample_id}.bamstats.txt
  "\$MM" run -p "\$ENV_PREFIX" samtools idxstats ${bam} > ${sample_id}.counts.tsv
  """
}

process compute_simple_ratios {

  publishDir "${params.out_dir}/summary", mode: 'copy', overwrite: true

  input:
    path sample_count_files
    path water_counts
    path ref_fa

  output:
    path "positive_samples.txt"
    path "sample_summary.txt"
    path "wf-shrimpseq-report.html"

  script:
  """
  set -euo pipefail

  alpha=${params.pseudocount}
  pos_cut=${params.pos_cutoff}
  inc_cut=${params.inc_cutoff}

  mapfile -t viruses < <(grep '^>' ${ref_fa} | sed 's/^>//; s/ .*//')

  declare -A water
  while read -r c l m u; do
    [[ "\$c" == "*" ]] && continue
    water["\$c"]="\$m"
  done < ${water_counts}

  {
    printf "Sample"
    for v in "\${viruses[@]}"; do
      printf "\\t%s_Reads\\t%s_SimpleRatio\\t%s_Status" "\$v" "\$v" "\$v"
    done
    printf "\\n"
  } > positive_samples.txt

  printf "Sample\\tOverall_Status\\tPOS_Viruses\\tINC_Viruses\\tNum_POS\\tNum_INC\\tEHP_Interpretation\\n" > sample_summary.txt

  for f in ${sample_count_files}; do
    sample="\${f%.counts.tsv}"

    declare -A sample_counts
    while read -r c l m u; do
      [[ "\$c" == "*" ]] && continue
      sample_counts["\$c"]="\$m"
    done < "\$f"

    overall="NEG"; pos_list=""; inc_list=""
    num_pos=0; num_inc=0
    ehp_ssu_status="NEG"; ehp_swp_status="NEG"

    {
      printf "%s" "\$sample"
      for v in "\${viruses[@]}"; do
        s="\${sample_counts[\$v]:-0}"
        w="\${water[\$v]:-0}"

        simple=\$(awk -v s="\$s" -v w="\$w" -v a="\$alpha" 'BEGIN{printf "%.3f",(s+a)/(w+a)}')

        status="NEG"
        awk -v x="\$simple" -v c="\$pos_cut" 'BEGIN{exit !(x>=c)}' && status="POS"
        awk -v x="\$simple" -v c="\$inc_cut" 'BEGIN{exit !(x>=c)}' && [[ "\$status" != "POS" ]] && status="INC"

        if [[ "\$v" == "EHP_SSUrRNA" ]]; then
          ehp_ssu_status="\$status"
        elif [[ "\$v" == "EHP_spore-capsule" ]]; then
          ehp_swp_status="\$status"
        else
          [[ "\$status" == "POS" ]] && { overall="POS"; num_pos=\$((num_pos+1)); pos_list=\${pos_list:+\$pos_list,}\$v; }
          [[ "\$status" == "INC" ]] && { [[ "\$overall" != "POS" ]] && overall="INC"; num_inc=\$((num_inc+1)); inc_list=\${inc_list:+\$inc_list,}\$v; }
        fi

        printf "\\t%s\\t%s\\t%s" "\$s" "\$simple" "\$status"
      done
      printf "\\n"
    } >> positive_samples.txt

    # Final EHP call combines mapping evidence from both independent targets.
    if [[ "\$ehp_ssu_status" == "POS" && "\$ehp_swp_status" == "POS" ]]; then
      ehp_interpretation="EHP_POSITIVE"
      overall="POS"
      num_pos=\$((num_pos+1))
      pos_list=\${pos_list:+\$pos_list,}EHP
    elif [[ "\$ehp_ssu_status" == "POS" && "\$ehp_swp_status" != "POS" ]]; then
      ehp_interpretation="EHP_SUSPECT_SSU_ONLY"
      [[ "\$overall" != "POS" ]] && overall="INC"
      num_inc=\$((num_inc+1))
      inc_list=\${inc_list:+\$inc_list,}EHP
    elif [[ "\$ehp_ssu_status" != "POS" && "\$ehp_swp_status" == "POS" ]]; then
      ehp_interpretation="EHP_SUSPECT_SWP_ONLY"
      [[ "\$overall" != "POS" ]] && overall="INC"
      num_inc=\$((num_inc+1))
      inc_list=\${inc_list:+\$inc_list,}EHP
    elif [[ "\$ehp_ssu_status" == "NEG" && "\$ehp_swp_status" == "NEG" ]]; then
      ehp_interpretation="EHP_NEGATIVE"
    else
      ehp_interpretation="EHP_INCONCLUSIVE"
      [[ "\$overall" != "POS" ]] && overall="INC"
      num_inc=\$((num_inc+1))
      inc_list=\${inc_list:+\$inc_list,}EHP
    fi

    [[ -z "\$pos_list" ]] && pos_list="-"
    [[ -z "\$inc_list" ]] && inc_list="-"

    printf "%s\\t%s\\t%s\\t%s\\t%d\\t%d\\t%s\\n" \
      "\$sample" "\$overall" "\$pos_list" "\$inc_list" "\$num_pos" "\$num_inc" "\$ehp_interpretation" \
      >> sample_summary.txt
  done

  REPORT_HTML="wf-shrimpseq-report.html"

  cat > "\$REPORT_HTML" <<'HTML'
<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width,initial-scale=1" />
<title>wf-shrimpseq report</title>
<style>
  :root{
    --pos-bg:#e7f7ee; --pos-fg:#0b5d2a;
    --inc-bg:#fff7e6; --inc-fg:#7a4b00;
    --neg-bg:#fdecec; --neg-fg:#7a0b0b;
    --muted:#666;
    --border:#ddd;
  }
  body{ font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif; margin:18px; color:#111; }
  h1{ font-size:20px; margin:0 0 8px; }
  h2{ font-size:16px; margin:18px 0 8px; }
  p{ margin:6px 0; color:var(--muted); }
  .card{ border:1px solid var(--border); border-radius:10px; padding:12px; margin:12px 0; background:#fff; }
  table{ width:100%; border-collapse:collapse; font-size:13px; }
  th, td{ border:1px solid var(--border); padding:6px 8px; vertical-align:top; }
  th{ background:#f6f6f6; position:sticky; top:0; z-index:1; }
  td.num{ text-align:right; white-space:nowrap; }
  .status{ font-weight:700; text-align:center; white-space:nowrap; }
  .POS{ background:var(--pos-bg); color:var(--pos-fg); }
  .INC{ background:var(--inc-bg); color:var(--inc-fg); }
  .NEG{ background:var(--neg-bg); color:var(--neg-fg); }
  .EHP_POSITIVE{ background:var(--pos-bg); color:var(--pos-fg); }
  .EHP_INCONCLUSIVE, .EHP_SUSPECT_SSU_ONLY, .EHP_SUSPECT_SWP_ONLY{ background:var(--inc-bg); color:var(--inc-fg); }
  .EHP_NEGATIVE{ background:var(--neg-bg); color:var(--neg-fg); }
  .muted{ color:var(--muted); }
  .small{ font-size:12px; }
  .wrap{ overflow:auto; max-height:70vh; border-radius:10px; }
</style>
</head>
<body>
  <h1>wf-shrimpseq — Pathogen detection report</h1>
  <p class="small muted">Colored statuses:
    <span class="status POS">POS</span>
    <span class="status INC">INC</span>
    <span class="status NEG">NEG</span>
  </p>

  <div class="card">
    <h2>Sample summary</h2>
    <div class="wrap">
      <table>
HTML

  awk -F'\\t' '
    NR==1{
      print "<thead><tr>";
      for(j=1;j<=NF;j++) print "<th>" \$j "</th>";
      print "</tr></thead><tbody>";
      next
    }
    {
      print "<tr>";
      for(j=1;j<=NF;j++){
        v=\$j;
        if(j==2 || j==7){
          cls=v;
          print "<td class=\\"status " cls "\\">" v "</td>";
        } else {
          print "<td>" v "</td>";
        }
      }
      print "</tr>";
    }
    END{ print "</tbody>"; }
  ' sample_summary.txt >> "\$REPORT_HTML"

  cat >> "\$REPORT_HTML" <<'HTML'
      </table>
    </div>
  </div>

  <div class="card">
    <h2>Per-pathogen results</h2>
    <p class="small muted">Status columns are color-coded.</p>
    <div class="wrap">
      <table>
HTML

  awk -F'\\t' '
    NR==1{
      print "<thead><tr>";
      for(j=1;j<=NF;j++) print "<th>" \$j "</th>";
      print "</tr></thead><tbody>";
      next
    }
    {
      print "<tr>";
      for(j=1;j<=NF;j++){
        v=\$j;
        if(j==1){
          print "<td>" v "</td>";
        }
        else if( (j-1)%3==1 ){
          print "<td class=\\"num\\">" v "</td>";
        }
        else if( (j-1)%3==2 ){
          print "<td class=\\"num\\">" v "</td>";
        }
        else {
          cls=v;
          print "<td class=\\"status " cls "\\">" v "</td>";
        }
      }
      print "</tr>";
    }
    END{ print "</tbody>"; }
  ' positive_samples.txt >> "\$REPORT_HTML"

  cat >> "\$REPORT_HTML" <<'HTML'
      </table>
    </div>
  </div>

</body>
</html>
HTML
  """
}
