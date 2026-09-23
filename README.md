# wf-shrimpseq

Workflow that maps ONT reads to reference shrimp pathogens and reports per-pathogen SimpleRatio
using a water control.

## Inputs
- `--reads` : FASTQ/FASTQ.GZ files (glob), e.g. `data/*.fastq.gz`
- `--ref`   : reference FASTA, e.g. `reference/genomes.fasta`

## Outputs
Written to `--out_dir` (default: `./results`):
- `bam/` : sorted BAM + BAI per sample
- `counts/` : `*.counts.tsv` (samtools idxstats) and `*.bamstats.txt`
- `summary/` :
  - `positive_samples.txt`
  - `sample_summary.txt` (includes the final two-target EHP interpretation)

## EHP interpretation

EHP is evaluated from the two independent targets `EHP_SSUrRNA` and
`EHP_spore-capsule` (SWP). `EHP_combined` is not used as a reference target.

- Both POS: `EHP_POSITIVE`
- SSU POS and SWP NEG/INC: `EHP_SUSPECT_SSU_ONLY`
- SSU NEG/INC and SWP POS: `EHP_SUSPECT_SWP_ONLY`
- Both NEG: `EHP_NEGATIVE`
- INC+NEG, NEG+INC, or INC+INC: `EHP_INCONCLUSIVE`

## Example run (local)
```bash
nextflow run kuttichantran/wf-viral-screen -profile standard \
  --reads "data/*.fastq.gz" \
  --ref "reference/genomes.fasta" \
  --out_dir results


## Reference FASTA (bundled)

This workflow ships with a bundled reference at `assets/reference.fasta`, so you do **not** need to upload/provide `--ref` for normal runs.

- To use a different reference, either:
  - replace `assets/reference.fasta` in the workflow package, **or**
  - run with `--ref /path/to/your/custom.fasta` to override.
