# RNA-sequencing pipeline

This repository contains a reproducible RNA-seq analysis for the uploaded
MM1S/HUVEC raw count matrix.

## Design

The workflow derives sample metadata from the count-table column names and
excludes all `FNE-*` samples before any analysis. The retained groups are:

- `MM1S_monoculture`: `MM1S monoculture-1`, `MM1S monoculture-2`,
  `MM1S monoculture-3`
- `MM1S_coculture`: `MM1S cocultured-1`, `MM1S cocultured-2`,
  `MM1S cocultured-3`
- `HUVEC_MM1S`: `HUVEC-MM1S-1`, `HUVEC-MM1S-2`, `HUVEC-MM1S-3`

The primary contrast is:

- `MM1S_coculture_vs_MM1S_monoculture`

Additional pairwise contrasts support the requested three-group DEG overlap
and heatmap:

- `HUVEC_MM1S_vs_MM1S_monoculture`
- `HUVEC_MM1S_vs_MM1S_coculture`

## Methods

Because this cloud image does not include R/DESeq2, the executed workflow uses
a Python implementation:

1. Filter low-count genes.
2. Normalize counts with DESeq2-style median-of-ratios size factors.
3. Test pairwise differential expression using a negative-binomial Wald
   approximation with method-of-moments dispersion estimates.
4. Adjust p-values with Benjamini-Hochberg FDR.
5. Call DEGs at `FDR < 0.05` and `abs(log2 fold change) >= 1`.
6. Run preranked GSEA for `MM1S_coculture_vs_MM1S_monoculture` against
   `MSigDB_Hallmark_2020` through `gseapy`.

## Run

Install dependencies:

```bash
python3 -m pip install --user -r requirements.txt
```

Run the analysis:

```bash
python3 scripts/run_rnaseq_analysis.py \
  --counts /path/to/gene_count_Shabnam_s_data_RNAseq_HUVEC_d1c5.txt \
  --outdir results
```

## Outputs

- `results/tables/DE_MM1S_coculture_vs_MM1S_monoculture.tsv`
- `results/tables/DE_HUVEC_MM1S_vs_MM1S_monoculture.tsv`
- `results/tables/DE_HUVEC_MM1S_vs_MM1S_coculture.tsv`
- `results/plots/venn_DEG_three_pairwise_contrasts.png`
- `results/plots/heatmap_top_DEG_three_groups.png`
- `results/plots/volcano_MM1S_coculture_vs_monoculture.png`
- `results/gsea/MM1S_coculture_vs_monoculture_hallmark/gsea_results.tsv`
- `results/gsea/MM1S_coculture_vs_monoculture_hallmark/gsea_top20_NES.png`
- `results/analysis_summary.md`
