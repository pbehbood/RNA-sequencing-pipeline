#!/usr/bin/env python3
"""RNA-seq analysis for MM1S/HUVEC count data.

The workflow intentionally excludes FNE samples and derives sample metadata
from the count-table headers supplied with the project.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import gseapy as gp
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib_venn import venn3
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests


GROUP_ORDER = ["MM1S_monoculture", "MM1S_coculture", "HUVEC_MM1S"]
CONTRASTS = [
    ("MM1S_coculture", "MM1S_monoculture"),
    ("HUVEC_MM1S", "MM1S_monoculture"),
    ("HUVEC_MM1S", "MM1S_coculture"),
]
MM1S_ONLY_GROUPS = ["MM1S_monoculture", "MM1S_coculture"]
MM1S_ONLY_CONTRASTS = [("MM1S_coculture", "MM1S_monoculture")]
DEFAULT_INPUT = (
    "/home/ubuntu/.cursor/projects/workspace/uploads/"
    "gene_count_Shabnam_s_data_RNAseq_HUVEC_d1c5.txt"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run RNA-seq DE, visualization, and GSEA for MM1S/HUVEC samples."
    )
    parser.add_argument(
        "--counts",
        default=DEFAULT_INPUT,
        help="Tab-delimited count table with gene_id, gene_name, and sample columns.",
    )
    parser.add_argument(
        "--outdir",
        default="results",
        help="Directory for tables, plots, and GSEA outputs.",
    )
    parser.add_argument(
        "--analysis-scope",
        choices=["three_group", "mm1s_only"],
        default="three_group",
        help=(
            "Use three_group for MM1S/HUVEC pairwise contrasts or mm1s_only "
            "for only MM1S monoculture and coculture samples."
        ),
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=10,
        help="Keep genes with at least this count in at least --min-samples samples.",
    )
    parser.add_argument(
        "--min-samples",
        type=int,
        default=3,
        help="Minimum number of samples passing --min-count for gene filtering.",
    )
    parser.add_argument(
        "--lfc-threshold",
        type=float,
        default=1.0,
        help="Absolute log2 fold-change threshold for DEG calls.",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="Benjamini-Hochberg FDR threshold for DEG calls.",
    )
    parser.add_argument(
        "--gsea-library",
        default="MSigDB_Hallmark_2020",
        help="Enrichr/gseapy gene-set library for preranked GSEA.",
    )
    parser.add_argument(
        "--gsea-permutations",
        type=int,
        default=1000,
        help="Number of permutations for preranked GSEA.",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed.")
    return parser.parse_args()


def sample_group(sample: str) -> str | None:
    if sample.startswith("FNE"):
        return None
    if sample.startswith("MM1S cocultured"):
        return "MM1S_coculture"
    if sample.startswith("MM1S monoculture"):
        return "MM1S_monoculture"
    if sample.startswith("HUVEC-MM1S"):
        return "HUVEC_MM1S"
    raise ValueError(f"Unrecognized non-FNE sample column: {sample}")


def make_metadata(columns: Iterable[str]) -> pd.DataFrame:
    rows = []
    for sample in columns:
        group = sample_group(sample)
        if group is None:
            continue
        replicate = sample.rsplit("-", 1)[-1]
        rows.append({"sample": sample, "group": group, "replicate": replicate})
    metadata = pd.DataFrame(rows)
    metadata["group"] = pd.Categorical(metadata["group"], categories=GROUP_ORDER, ordered=True)
    return metadata.sort_values(["group", "replicate"]).reset_index(drop=True)


def read_counts(path: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t")
    required = {"gene_id", "gene_name"}
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f"Count table is missing required columns: {sorted(missing)}")

    metadata = make_metadata(raw.columns[2:])
    sample_cols = metadata["sample"].tolist()
    counts = raw[["gene_id", "gene_name", *sample_cols]].copy()
    counts["gene_name"] = counts["gene_name"].fillna(counts["gene_id"])
    counts[sample_cols] = counts[sample_cols].apply(pd.to_numeric, errors="raise").astype(float)
    return raw, counts, metadata


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    positive = counts.gt(0).all(axis=1)
    if positive.sum() < 100:
        library_sizes = counts.sum(axis=0)
        return library_sizes / np.median(library_sizes)

    log_counts = np.log(counts.loc[positive])
    geometric_means = np.exp(log_counts.mean(axis=1))
    ratios = counts.loc[positive].divide(geometric_means, axis=0)
    size_factors = ratios.median(axis=0)
    if (size_factors <= 0).any() or size_factors.isna().any():
        library_sizes = counts.sum(axis=0)
        size_factors = library_sizes / np.median(library_sizes)
    return size_factors


def bh_fdr(pvalues: pd.Series) -> np.ndarray:
    finite = np.isfinite(pvalues.to_numpy())
    padj = np.full(len(pvalues), np.nan)
    if finite.any():
        padj[finite] = multipletests(pvalues.to_numpy()[finite], method="fdr_bh")[1]
    return padj


def call_status(row: pd.Series, lfc_threshold: float, fdr_threshold: float) -> str:
    if not np.isfinite(row["padj"]) or row["padj"] >= fdr_threshold:
        return "not_significant"
    if row["log2FoldChange"] >= lfc_threshold:
        return "up"
    if row["log2FoldChange"] <= -lfc_threshold:
        return "down"
    return "fdr_only"


def differential_expression(
    gene_info: pd.DataFrame,
    norm_counts: pd.DataFrame,
    log_cpm: pd.DataFrame,
    metadata: pd.DataFrame,
    numerator_group: str,
    denominator_group: str,
    lfc_threshold: float,
    fdr_threshold: float,
) -> pd.DataFrame:
    num_samples = metadata.loc[metadata["group"] == numerator_group, "sample"].tolist()
    den_samples = metadata.loc[metadata["group"] == denominator_group, "sample"].tolist()
    if len(num_samples) < 2 or len(den_samples) < 2:
        raise ValueError(f"Need at least two replicates for {numerator_group} vs {denominator_group}")

    pseudo = 0.5
    num_mean = norm_counts[num_samples].mean(axis=1)
    den_mean = norm_counts[den_samples].mean(axis=1)
    base_mean = norm_counts[metadata["sample"].tolist()].mean(axis=1)
    log2fc = np.log2((num_mean + pseudo) / (den_mean + pseudo))

    global_mean = norm_counts.mean(axis=1)
    global_var = norm_counts.var(axis=1, ddof=1)
    dispersion = ((global_var - global_mean) / np.square(global_mean.replace(0, np.nan))).clip(
        lower=1e-4, upper=100.0
    )
    dispersion = dispersion.fillna(1.0)

    ln_fc = np.log((num_mean + pseudo) / (den_mean + pseudo))
    se = np.sqrt(
        ((1.0 / (num_mean + pseudo)) + dispersion) / len(num_samples)
        + ((1.0 / (den_mean + pseudo)) + dispersion) / len(den_samples)
    )
    wald_stat = ln_fc / se
    pvalue = pd.Series(2.0 * stats.norm.sf(np.abs(wald_stat)), index=norm_counts.index)

    welch = stats.ttest_ind(
        log_cpm[num_samples].to_numpy(),
        log_cpm[den_samples].to_numpy(),
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )

    result = gene_info.copy()
    result.insert(0, "gene_id", result.index)
    result["baseMean"] = base_mean
    result[f"mean_{numerator_group}"] = num_mean
    result[f"mean_{denominator_group}"] = den_mean
    result["log2FoldChange"] = log2fc
    result["dispersion_method_of_moments"] = dispersion
    result["stat"] = wald_stat
    result["pvalue"] = pvalue
    result["padj"] = bh_fdr(pvalue)
    result["welch_logCPM_pvalue"] = welch.pvalue
    result["comparison"] = f"{numerator_group}_vs_{denominator_group}"
    result["status"] = result.apply(
        call_status,
        axis=1,
        lfc_threshold=lfc_threshold,
        fdr_threshold=fdr_threshold,
    )
    return result.sort_values(["padj", "pvalue", "gene_name"], na_position="last")


def safe_neg_log10(values: pd.Series) -> pd.Series:
    clipped = values.clip(lower=np.nextafter(0, 1))
    return -np.log10(clipped)


def save_qc_plots(
    raw_counts: pd.DataFrame,
    norm_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    plots_dir: Path,
) -> None:
    sns.set_theme(style="whitegrid")
    sample_cols = metadata["sample"].tolist()
    library_sizes = raw_counts[sample_cols].sum(axis=0)
    lib_df = (
        pd.DataFrame({"sample": library_sizes.index, "library_size": library_sizes.values})
        .merge(metadata, on="sample")
        .sort_values(["group", "sample"])
    )
    plt.figure(figsize=(10, 5))
    ax = sns.barplot(data=lib_df, x="sample", y="library_size", hue="group", dodge=False)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_title("Library sizes after excluding FNE samples")
    ax.set_xlabel("")
    ax.set_ylabel("Raw reads assigned to genes")
    plt.tight_layout()
    plt.savefig(plots_dir / "qc_library_sizes.png", dpi=300)
    plt.savefig(plots_dir / "qc_library_sizes.pdf")
    plt.close()

    log_norm = np.log2(norm_counts[sample_cols] + 1).T
    scaled = StandardScaler().fit_transform(log_norm)
    pca = PCA(n_components=2, random_state=0)
    coords = pca.fit_transform(scaled)
    pca_df = metadata.copy()
    pca_df["PC1"] = coords[:, 0]
    pca_df["PC2"] = coords[:, 1]
    plt.figure(figsize=(7, 6))
    ax = sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="group", s=110)
    for _, row in pca_df.iterrows():
        ax.text(row["PC1"], row["PC2"], row["replicate"], fontsize=9, ha="left", va="bottom")
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    ax.set_title("PCA of log2 normalized counts")
    plt.tight_layout()
    plt.savefig(plots_dir / "qc_pca.png", dpi=300)
    plt.savefig(plots_dir / "qc_pca.pdf")
    plt.close()


def save_volcano(de: pd.DataFrame, plots_dir: Path, lfc_threshold: float, fdr_threshold: float) -> None:
    plot_df = de.copy()
    plot_df["neg_log10_padj"] = safe_neg_log10(plot_df["padj"].fillna(1.0))
    plot_df["volcano_class"] = np.select(
        [
            (plot_df["padj"] < fdr_threshold) & (plot_df["log2FoldChange"] >= lfc_threshold),
            (plot_df["padj"] < fdr_threshold) & (plot_df["log2FoldChange"] <= -lfc_threshold),
        ],
        ["Up in MM1S coculture", "Down in MM1S coculture"],
        default="Not significant",
    )
    palette = {
        "Up in MM1S coculture": "#d62728",
        "Down in MM1S coculture": "#1f77b4",
        "Not significant": "#bdbdbd",
    }
    plt.figure(figsize=(8, 7))
    ax = sns.scatterplot(
        data=plot_df,
        x="log2FoldChange",
        y="neg_log10_padj",
        hue="volcano_class",
        palette=palette,
        linewidth=0,
        s=12,
        alpha=0.75,
    )
    ax.axvline(lfc_threshold, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_threshold, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-math.log10(fdr_threshold), color="black", linestyle="--", linewidth=0.8)
    label_df = plot_df.loc[plot_df["padj"] < fdr_threshold].copy()
    label_df["rank_score"] = label_df["neg_log10_padj"] * label_df["log2FoldChange"].abs()
    for _, row in label_df.nlargest(15, "rank_score").iterrows():
        ax.text(row["log2FoldChange"], row["neg_log10_padj"], row["gene_name"], fontsize=8)
    ax.set_title("MM1S coculture vs MM1S monoculture")
    ax.set_xlabel("log2 fold change (coculture / monoculture)")
    ax.set_ylabel("-log10 FDR")
    plt.tight_layout()
    plt.savefig(plots_dir / "volcano_MM1S_coculture_vs_monoculture.png", dpi=300)
    plt.savefig(plots_dir / "volcano_MM1S_coculture_vs_monoculture.pdf")
    plt.close()


def save_venn(de_results: dict[str, pd.DataFrame], plots_dir: Path) -> pd.DataFrame:
    set_names = {
        "MM1S_coculture_vs_MM1S_monoculture": "MM1S coculture\nvs monoculture",
        "HUVEC_MM1S_vs_MM1S_monoculture": "HUVEC-MM1S\nvs MM1S mono",
        "HUVEC_MM1S_vs_MM1S_coculture": "HUVEC-MM1S\nvs MM1S coculture",
    }
    sig_sets = {
        name: set(df.loc[df["status"].isin(["up", "down"]), "gene_id"])
        for name, df in de_results.items()
    }
    ordered = list(set_names)
    plt.figure(figsize=(8, 7))
    venn3([sig_sets[name] for name in ordered], set_labels=[set_names[name] for name in ordered])
    plt.title("Overlap of differentially expressed genes")
    plt.tight_layout()
    plt.savefig(plots_dir / "venn_DEG_three_pairwise_contrasts.png", dpi=300)
    plt.savefig(plots_dir / "venn_DEG_three_pairwise_contrasts.pdf")
    plt.close()

    all_gene_ids = sorted(set().union(*sig_sets.values()))
    rows = []
    for gene_id in all_gene_ids:
        rows.append(
            {
                "gene_id": gene_id,
                **{name: gene_id in genes for name, genes in sig_sets.items()},
            }
        )
    return pd.DataFrame(rows)


def save_heatmap(
    norm_counts: pd.DataFrame,
    metadata: pd.DataFrame,
    de_results: dict[str, pd.DataFrame],
    plots_dir: Path,
    tables_dir: Path,
    all_degs: bool = False,
    filename_prefix: str = "heatmap_top_DEG_three_groups",
    title: str = "Top differentially expressed genes across MM1S and HUVEC-MM1S groups",
) -> None:
    sig_gene_ids = set()
    for df in de_results.values():
        sig_gene_ids.update(df.loc[df["status"].isin(["up", "down"]), "gene_id"])

    if sig_gene_ids:
        rank_rows = []
        for name, df in de_results.items():
            tmp = df.loc[df["gene_id"].isin(sig_gene_ids), ["gene_id", "gene_name", "padj", "log2FoldChange"]].copy()
            tmp["contrast"] = name
            rank_rows.append(tmp)
        rank_df = pd.concat(rank_rows, ignore_index=True)
        rank_df["abs_lfc"] = rank_df["log2FoldChange"].abs()
        ranked_gene_ids = (
            rank_df.sort_values(["padj", "abs_lfc"], ascending=[True, False])
            .drop_duplicates("gene_id")["gene_id"]
            .tolist()
        )
        top_gene_ids = ranked_gene_ids if all_degs else ranked_gene_ids[:100]
    else:
        main = de_results["MM1S_coculture_vs_MM1S_monoculture"]
        top_gene_ids = main.head(100)["gene_id"].tolist()

    sample_cols = metadata["sample"].tolist()
    matrix = np.log2(norm_counts.loc[top_gene_ids, sample_cols] + 1)
    z = matrix.sub(matrix.mean(axis=1), axis=0).div(matrix.std(axis=1).replace(0, np.nan), axis=0)
    z = z.fillna(0).clip(-2.5, 2.5)
    gene_labels = norm_counts.loc[top_gene_ids, ["gene_name"]].copy()
    labels = gene_labels["gene_name"].where(gene_labels["gene_name"].notna(), gene_labels.index.to_series())
    labels = labels + " (" + gene_labels.index.to_series().values + ")"
    z.index = labels

    ordered_samples = (
        metadata.sort_values(["group", "replicate"])["sample"].tolist()
    )
    z = z[ordered_samples]
    z.to_csv(tables_dir / f"{filename_prefix}_zscores.tsv", sep="\t")

    group_palette = {
        "MM1S_monoculture": "#4c78a8",
        "MM1S_coculture": "#f58518",
        "HUVEC_MM1S": "#54a24b",
    }
    col_colors = metadata.set_index("sample").loc[ordered_samples, "group"].map(group_palette)
    g = sns.clustermap(
        z,
        cmap="vlag",
        center=0,
        col_cluster=False,
        row_cluster=True,
        col_colors=col_colors,
        figsize=(10, max(6, min(24, 1.2 + 0.55 * len(top_gene_ids)))),
        xticklabels=True,
        yticklabels=True,
    )
    g.fig.suptitle(title, y=1.01)
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.savefig(plots_dir / f"{filename_prefix}.png", dpi=300, bbox_inches="tight")
    g.savefig(plots_dir / f"{filename_prefix}.pdf", bbox_inches="tight")
    plt.close(g.fig)


def save_deg_list(de: pd.DataFrame, tables_dir: Path) -> pd.DataFrame:
    degs = de.loc[de["status"].isin(["up", "down"])].copy()
    degs = degs.sort_values(["status", "padj", "gene_name"], ascending=[False, True, True])
    columns = [
        "gene_id",
        "gene_name",
        "status",
        "log2FoldChange",
        "padj",
        "pvalue",
        "baseMean",
        "mean_MM1S_coculture",
        "mean_MM1S_monoculture",
    ]
    degs[columns].to_csv(
        tables_dir / "DEG_list_MM1S_coculture_vs_monoculture.tsv",
        sep="\t",
        index=False,
    )
    return degs


def run_gsea(main_de: pd.DataFrame, outdir: Path, library: str, permutations: int, seed: int) -> pd.DataFrame:
    gsea_dir = outdir / "gsea" / "MM1S_coculture_vs_monoculture_hallmark"
    gsea_dir.mkdir(parents=True, exist_ok=True)

    ranking = main_de[["gene_name", "log2FoldChange", "pvalue"]].copy()
    ranking = ranking.replace([np.inf, -np.inf], np.nan).dropna()
    ranking = ranking[ranking["gene_name"].astype(str).str.len() > 0]
    ranking["rank_score"] = np.sign(ranking["log2FoldChange"]) * safe_neg_log10(ranking["pvalue"])
    ranking = (
        ranking.sort_values("rank_score", key=lambda s: s.abs(), ascending=False)
        .drop_duplicates("gene_name")
        .sort_values("rank_score", ascending=False)
    )
    rnk_path = gsea_dir / "MM1S_coculture_vs_monoculture_prerank.rnk"
    ranking[["gene_name", "rank_score"]].to_csv(rnk_path, sep="\t", index=False, header=False)

    prerank = gp.prerank(
        rnk=str(rnk_path),
        gene_sets=library,
        outdir=str(gsea_dir),
        min_size=15,
        max_size=500,
        permutation_num=permutations,
        seed=seed,
        threads=4,
        format="png",
        verbose=False,
    )
    result = prerank.res2d.copy()
    result.to_csv(gsea_dir / "gsea_results.tsv", sep="\t", index=False)

    q_col = "FDR q-val" if "FDR q-val" in result.columns else "FDR"
    nes_col = "NES"
    term_col = "Term"
    top = result.sort_values(q_col).head(20).copy()
    top[nes_col] = pd.to_numeric(top[nes_col], errors="coerce")
    top[q_col] = pd.to_numeric(top[q_col], errors="coerce")
    top = top.sort_values(nes_col)

    plt.figure(figsize=(10, 8))
    colors = np.where(top[nes_col] >= 0, "#d62728", "#1f77b4")
    plt.barh(top[term_col], top[nes_col], color=colors)
    plt.axvline(0, color="black", linewidth=0.8)
    plt.xlabel("Normalized enrichment score (NES)")
    plt.title(f"Top preranked GSEA pathways ({library})")
    plt.tight_layout()
    plt.savefig(gsea_dir / "gsea_top20_NES.png", dpi=300)
    plt.savefig(gsea_dir / "gsea_top20_NES.pdf")
    plt.close()
    return result


def write_summary(
    outdir: Path,
    metadata: pd.DataFrame,
    filtered_genes: int,
    de_results: dict[str, pd.DataFrame],
    gsea_result: pd.DataFrame,
    fdr_threshold: float,
    lfc_threshold: float,
    gsea_library: str,
    analysis_scope: str,
) -> None:
    output_root = outdir.as_posix()
    sample_rows = ["| sample | group | replicate |", "| --- | --- | --- |"]
    for _, row in metadata.iterrows():
        sample_rows.append(f"| {row['sample']} | {row['group']} | {row['replicate']} |")

    lines = [
        "# RNA-seq analysis summary",
        "",
        (
            "FNE samples were excluded before normalization, differential expression, "
            "and plotting."
        ),
        (
            "For the MM1S-only scope, HUVEC-MM1S samples were also excluded before "
            "filtering and normalization."
        )
        if analysis_scope == "mm1s_only"
        else "",
        "",
        "## Samples used",
        "",
        *sample_rows,
        "",
        "## Differential-expression settings",
        "",
        f"- Filter: raw count >= threshold in at least 3 retained samples.",
        f"- DEG call: FDR < {fdr_threshold} and abs(log2 fold change) >= {lfc_threshold}.",
        "- Normalization: DESeq2-style median-of-ratios size factors.",
        "- Test: fast negative-binomial Wald approximation using method-of-moments dispersion.",
        f"- Filtered genes tested: {filtered_genes:,}.",
        "",
        "## DEG counts",
        "",
        "| contrast | up | down | total |",
        "| --- | ---: | ---: | ---: |",
    ]
    for name, df in de_results.items():
        up = int((df["status"] == "up").sum())
        down = int((df["status"] == "down").sum())
        lines.append(f"| {name} | {up:,} | {down:,} | {up + down:,} |")

    q_col = "FDR q-val" if "FDR q-val" in gsea_result.columns else "FDR"
    sig_gsea = int((pd.to_numeric(gsea_result[q_col], errors="coerce") < 0.25).sum())
    lines.extend(
        [
            "",
            "## GSEA",
            "",
            f"- Library: {gsea_library}.",
            f"- Pathways with GSEA FDR q-val < 0.25: {sig_gsea:,}.",
            "",
            f"Primary outputs are in `{output_root}/tables`, `{output_root}/plots`, and `{output_root}/gsea`.",
            "",
        ]
    )
    (outdir / "analysis_summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    np.random.seed(args.seed)
    outdir = Path(args.outdir)
    tables_dir = outdir / "tables"
    plots_dir = outdir / "plots"
    gsea_dir = outdir / "gsea"
    for directory in (tables_dir, plots_dir, gsea_dir):
        directory.mkdir(parents=True, exist_ok=True)

    raw, counts_with_info, metadata = read_counts(args.counts)
    if args.analysis_scope == "mm1s_only":
        included_groups = MM1S_ONLY_GROUPS
        contrasts = MM1S_ONLY_CONTRASTS
    else:
        included_groups = GROUP_ORDER
        contrasts = CONTRASTS

    metadata = metadata.loc[metadata["group"].isin(included_groups)].copy()
    metadata["group"] = metadata["group"].cat.remove_unused_categories()
    counts_with_info = counts_with_info[["gene_id", "gene_name", *metadata["sample"].tolist()]].copy()
    metadata.to_csv(tables_dir / "sample_metadata.tsv", sep="\t", index=False)

    sample_cols = metadata["sample"].tolist()
    gene_info = counts_with_info[["gene_id", "gene_name"]].copy()
    raw_counts = counts_with_info.set_index("gene_id")[sample_cols]
    gene_info = gene_info.set_index("gene_id")

    keep = (raw_counts >= args.min_count).sum(axis=1) >= args.min_samples
    filtered_raw_counts = raw_counts.loc[keep].copy()
    filtered_gene_info = gene_info.loc[keep].copy()

    size_factors = median_ratio_size_factors(filtered_raw_counts)
    size_factors.to_frame("size_factor").to_csv(tables_dir / "size_factors.tsv", sep="\t")
    norm_counts = filtered_raw_counts.divide(size_factors, axis=1)
    norm_counts = filtered_gene_info.join(norm_counts)
    norm_counts.to_csv(tables_dir / "normalized_counts_filtered.tsv", sep="\t")

    norm_numeric = norm_counts[sample_cols]
    library_sizes = filtered_raw_counts.sum(axis=0)
    cpm = filtered_raw_counts.divide(library_sizes, axis=1) * 1_000_000
    log_cpm = np.log2(cpm + 0.5)

    save_qc_plots(filtered_raw_counts, norm_numeric, metadata, plots_dir)

    de_results: dict[str, pd.DataFrame] = {}
    for numerator, denominator in contrasts:
        name = f"{numerator}_vs_{denominator}"
        de = differential_expression(
            filtered_gene_info,
            norm_numeric,
            log_cpm,
            metadata,
            numerator,
            denominator,
            args.lfc_threshold,
            args.fdr_threshold,
        )
        de_results[name] = de
        de.to_csv(tables_dir / f"DE_{name}.tsv", sep="\t", index=False)

    main_name = "MM1S_coculture_vs_MM1S_monoculture"
    save_volcano(de_results[main_name], plots_dir, args.lfc_threshold, args.fdr_threshold)
    deg_list = save_deg_list(de_results[main_name], tables_dir)
    if args.analysis_scope == "three_group":
        venn_membership = save_venn(de_results, plots_dir)
        if not venn_membership.empty:
            venn_membership = venn_membership.merge(
                gene_info.reset_index(), on="gene_id", how="left"
            )
        venn_membership.to_csv(tables_dir / "venn_DEG_membership.tsv", sep="\t", index=False)
        save_heatmap(norm_counts, metadata, de_results, plots_dir, tables_dir)
    else:
        save_heatmap(
            norm_counts,
            metadata,
            de_results,
            plots_dir,
            tables_dir,
            all_degs=True,
            filename_prefix="heatmap_all_DEG_MM1S_coculture_vs_monoculture",
            title="All DEGs: MM1S coculture vs MM1S monoculture",
        )

    gsea_result = run_gsea(
        de_results[main_name],
        outdir,
        args.gsea_library,
        args.gsea_permutations,
        args.seed,
    )
    write_summary(
        outdir,
        metadata,
        int(keep.sum()),
        de_results,
        gsea_result,
        args.fdr_threshold,
        args.lfc_threshold,
        args.gsea_library,
        args.analysis_scope,
    )

    run_info = {
        "input": str(Path(args.counts).resolve()),
        "samples_excluded": [col for col in raw.columns[2:] if sample_group(col) is None],
        "samples_not_used": [col for col in raw.columns[2:] if col not in sample_cols],
        "analysis_scope": args.analysis_scope,
        "samples_used": sample_cols,
        "filtered_genes_tested": int(keep.sum()),
        "contrasts": list(de_results),
        "deg_count_main_contrast": int(len(deg_list)),
        "fdr_threshold": args.fdr_threshold,
        "lfc_threshold": args.lfc_threshold,
        "gsea_library": args.gsea_library,
        "gsea_permutations": args.gsea_permutations,
    }
    (outdir / "run_info.json").write_text(json.dumps(run_info, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
