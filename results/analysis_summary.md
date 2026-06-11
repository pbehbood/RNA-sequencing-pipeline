# RNA-seq analysis summary

FNE samples were excluded before normalization, differential expression, and plotting.

## Samples used

| sample | group | replicate |
| --- | --- | --- |
| MM1S monoculture-1 | MM1S_monoculture | 1 |
| MM1S monoculture-2 | MM1S_monoculture | 2 |
| MM1S monoculture-3 | MM1S_monoculture | 3 |
| MM1S cocultured-1 | MM1S_coculture | 1 |
| MM1S cocultured-2 | MM1S_coculture | 2 |
| MM1S cocultured-3 | MM1S_coculture | 3 |
| HUVEC-MM1S-1 | HUVEC_MM1S | 1 |
| HUVEC-MM1S-2 | HUVEC_MM1S | 2 |
| HUVEC-MM1S-3 | HUVEC_MM1S | 3 |

## Differential-expression settings

- Filter: raw count >= threshold in at least 3 retained samples.
- DEG call: FDR < 0.05 and abs(log2 fold change) >= 1.0.
- Normalization: DESeq2-style median-of-ratios size factors.
- Test: fast negative-binomial Wald approximation using method-of-moments dispersion.
- Filtered genes tested: 17,573.

## DEG counts

| contrast | up | down | total |
| --- | ---: | ---: | ---: |
| MM1S_coculture_vs_MM1S_monoculture | 6 | 2 | 8 |
| HUVEC_MM1S_vs_MM1S_monoculture | 7 | 3 | 10 |
| HUVEC_MM1S_vs_MM1S_coculture | 3 | 2 | 5 |

## GSEA

- Library: MSigDB_Hallmark_2020.
- Pathways with GSEA FDR q-val < 0.25: 39.

Primary outputs are in `results/tables`, `results/plots`, and `results/gsea`.
