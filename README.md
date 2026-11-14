# SalmonTElocus

Per-locus quantification of transposable elements (TEs) from RNA-seq with
[Salmon](https://github.com/COMBINE-lab/salmon). Every TE copy in the genome
gets its own count and its genomic coordinates.

SalmonTElocus is a rewrite of [SalmonTE](https://github.com/hyunhwan-jeong/SalmonTE)
(Jeong et al., 2018). SalmonTE quantifies one consensus sequence per TE
subfamily, so it cannot say which copy a read came from. SalmonTElocus cuts
each copy listed in a TE annotation out of the genome and quantifies the
copies individually. Subfamily totals are still reported as the sum of each
subfamily's copies.

Like SalmonTE, the index contains **only TEs**. No reads are assigned to
genes.

## Requirements

- Python ≥ 3.8 with pandas
- Salmon (tested with 2.3.0)
- For `test` only: R with DESeq2

```sh
conda env create -f environment.yml
conda activate salmontelocus
```

## Usage

### 1. Build the index (once per genome + annotation)

```sh
./SalmonTElocus.py index \
  --genome GRCh38.primary_assembly.genome.fa \
  --te-gtf hg38_rmsk_TE.gtf \
  --out hg38_TE_loci --threads 8
```

`--te-gtf` takes a TEtranscripts-format TE GTF: one line per copy, with
`gene_id` (subfamily), `transcript_id` (copy), `family_id` and `class_id`.
These are the files distributed with TEtranscripts and TElocal, such as
`hg38_rmsk_TE.gtf`. Both the genome and the GTF may be gzipped.

- Soft-masked genomes are handled: all sequence is uppercased, because
  soft-masking writes the TEs in lowercase.
- Copies on the minus strand are reverse-complemented.
- Identical copies are kept as separate targets (`--keepDuplicates`).
- Copies on sequences that are missing from the genome FASTA are skipped with
  a warning. The usual causes are `chr1` vs `1` naming, or UCSC alt-contig
  names against a GENCODE genome. Copies shorter than `--min-length` (default
  31 bp, Salmon's k-mer size) are also skipped with a warning.

### 2. Quantify

```sh
./SalmonTElocus.py quant --index hg38_TE_loci --out results --threads 8 fastq_dir/
```

Inputs are FASTQ files and/or directories (`.fastq`, `.fq`, optionally
`.gz`). Mates are paired by file name: `S_R1`/`S_R2`, `S_1`/`S_2`, and
Illumina `S_L001_R1_001`. A file without a mate is treated as single-end. For
anything else, such as several lanes per sample, use a sample sheet,
`--samples samples.csv`, with the columns `sample,r1,r2` (leave `r2` empty for
single-end).

Add `--gene-gtf gencode.v44.primary_assembly.annotation.gtf` to label each TE
copy by the genes it overlaps. This is annotation only: **reads are never
counted for genes**, and the counts are identical with or without it. Four
columns are added after `length`:

| Column | Meaning |
|---|---|
| `gene_context` | `exon` (overlaps an exon of any gene by ≥ 1 bp), `intron` (inside a gene but no exon), or `intergenic` |
| `overlapping_genes` | Gene names (`;`-separated): the genes with an overlapping exon for `exon`, the host genes for `intron` |
| `gene_types` | Matching `gene_type` values (protein_coding, lncRNA, …) |
| `gene_strand` | TE strand relative to the gene(s): `sense`, `antisense` or `both` |

The log reports how many expressed copies, and what share of TE reads, fall
in each context. To add the labels to an existing run, re-run `quant` with
`--gene-gtf`; Salmon is skipped for samples that are already done, so only the
tables are rebuilt.

Samples that already have a `quant.sf` are skipped, so an interrupted run can
be restarted. Use `--force` to redo them. Extra Salmon options can be passed
with `--salmon-args "--useEM --seqBias"`.

Output in `results/`:

| File | Content |
|---|---|
| `locus_counts.csv` | One row per TE copy that has reads (`--keep-zero` keeps all copies): `locus_id, chrom, start, end, strand, transcript_id, subfamily, family, class, length`, then one count column per sample |
| `locus_tpm.csv` | Same layout, with TPM values |
| `subfamily_counts.csv`, `subfamily_tpm.csv` | Copies summed per subfamily |
| `expressed_loci.bed` | BED6 file of copies with reads, for IGV or the UCSC browser (score = mean count, capped at 1000) |
| `mapping_info.csv` | Reads processed and mapped, and the detected library type, per sample |
| `conditions.csv` | Template to fill in for `test` |
| `salmon/<sample>/` | Raw Salmon output |

Locus IDs follow TElocal / TEfastR (`transcript_id:gene_id:family_id:class_id`,
for example `L1HS_dup1234:L1HS:L1:LINE`), so results can be joined directly
with TEfastR's `te_locus_counts()`. Coordinates are 1-based and inclusive, as
in the GTF; the BED file is 0-based, as BED requires.

### 3. Differential expression

Fill in `results/conditions.csv` (columns `sample,condition`; samples marked
`NA` are left out), then run:

```sh
# two groups, reference first
./SalmonTElocus.py test --in results --out de --level locus --conditions control,treatment
# numeric covariate (e.g. age): leave out --conditions
./SalmonTElocus.py test --in results --out de_age --level subfamily
```

DESeq2 runs on the rounded counts after dropping features with fewer than
`--min-count` reads (default 10). If there are too few features to fit a
dispersion trend, it falls back to gene-wise dispersions. Output:
`results_<level>.csv` (annotation + coordinates + DESeq2 statistics, sorted by
padj), `normalized_counts_<level>.csv`, `summary_by_class_<level>.csv` and
`ma_plot_<level>.pdf`. When `quant` was run with `--gene-gtf`, the locus
results also carry the gene-context columns, and
`summary_by_class_and_context_locus.csv` breaks the summary down by
exon / intron / intergenic.

## Interpreting per-locus counts

- **Young copies.** Near-identical copies (for example L1HS, AluY, SVA_E/F,
  HERVK) share most of their reads. Salmon's EM splits those reads between the
  copies, so per-copy counts for these subfamilies are estimates; their
  subfamily totals are robust. Older, divergent copies get mostly uniquely
  mapped reads, and their per-copy counts are reliable.
- **TEs inside genes.** Because no genes are in the index, any read that
  contains TE sequence is counted as TE. That includes Alus in 3′ UTRs and TE
  sequence in introns of pre-mRNA. Highly expressed genes that carry an
  embedded TE therefore show up as highly expressed loci. Use `--gene-gtf` and
  filter on `gene_context` before calling a hit autonomous TE expression. In a
  test on two testis samples and one liver sample (E-MTAB-2836), 44% of TE
  reads fell on copies overlapping exons and 42% on intronic copies.
  Testis/liver enrichment of HERVK-int was 6.8× across all copies but 38×
  across intergenic copies only.
- **Resources.** A test index of 82,467 hg38 copies (all chr21 TEs plus young
  L1, Alu, SVA and HERVK genome-wide) took 24 s and 1.8 GB of memory to build;
  quantifying 1M read pairs took 6 s. A full `hg38_rmsk_TE.gtf` index
  (about 4.7M copies) has not been benchmarked and will need far more memory,
  so build it on a cluster node.

## Tests

```sh
python3 -m unittest discover -s tests -v
```

The end-to-end tests build a synthetic soft-masked genome and check five
things: coordinates, minus-strand copies, an even split between identical
copies, subfamily sums, and a known 2× change recovered by DESeq2. They are
skipped if Salmon (or R/DESeq2) is not installed.

## Citation and licence

Built on SalmonTE. If you use this tool, please cite:

> Jeong, H.-H., Yalamanchili, H. K., Guo, C., Shulman, J. M., & Liu, Z. (2018).
> An ultra-fast and scalable quantification pipeline for transposable elements
> from next generation sequencing data. *Pacific Symposium on Biocomputing*, 23, 168–179.

and Salmon (Patro et al., *Nature Methods* 2017). Distributed under the GNU
GPL v3, as SalmonTE is (see `LICENSE`).
