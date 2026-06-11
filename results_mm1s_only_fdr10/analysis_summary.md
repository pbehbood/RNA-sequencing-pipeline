# RNA-seq analysis summary

FNE samples were excluded before normalization, differential expression, and plotting.
For the MM1S-only scope, HUVEC-MM1S samples were also excluded before filtering and normalization.

## Samples used

| sample | group | replicate |
| --- | --- | --- |
| MM1S monoculture-1 | MM1S_monoculture | 1 |
| MM1S monoculture-2 | MM1S_monoculture | 2 |
| MM1S monoculture-3 | MM1S_monoculture | 3 |
| MM1S cocultured-1 | MM1S_coculture | 1 |
| MM1S cocultured-2 | MM1S_coculture | 2 |
| MM1S cocultured-3 | MM1S_coculture | 3 |

## Differential-expression settings

- Filter: raw count >= threshold in at least 3 retained samples.
- DEG call: FDR < 0.1 and abs(log2 fold change) >= 1.0.
- Normalization: DESeq2-style median-of-ratios size factors.
- Test: fast negative-binomial Wald approximation using method-of-moments dispersion.
- Filtered genes tested: 16,891.

## DEG counts

| contrast | up | down | total |
| --- | ---: | ---: | ---: |
| MM1S_coculture_vs_MM1S_monoculture | 482 | 1,238 | 1,720 |

## GSEA

- Library: MSigDB_Hallmark_2020.
- Pathways with GSEA FDR q-val < 0.25: 40.

Primary outputs are in `results_mm1s_only_fdr10/tables`, `results_mm1s_only_fdr10/plots`, and `results_mm1s_only_fdr10/gsea`.
