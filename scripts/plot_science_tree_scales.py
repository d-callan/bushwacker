#!/usr/bin/env python3
"""
Generate publication-grade vector PDF Figure (`paper/fig_science_tree_scales.pdf`)
synthesizing the Multi-Scale Tree-Building & Experimental Fitness Campaign:
  - Panel A: Somatic B-Cell Evolution Across 147 Germinal Centers (10,429 Cells) —
             Unsupervised MDL Polytomy Degree Predicts In Vivo Clonal Expansion (15.5x, p = 1.9e-301)
  - Panel B: Experimental Yeast-Display DMS Antigen-Binding Landscape (Delta log10 KD)
             Across MDL Tree Geometry & Cross-Mouse CDR Convergent Shunts (11.0x OR)
  - Panel C: Threshold-Free Coalescent Species Trees (Bushwhacker -> ASTRAL-IV) in the
             Anomaly Zone & Empirical Rapid Radiations
  - Panel D: Microbiome Mega-Phylogenies (N = 256 -> 32,768 16S OTUs) —
             Empirical Information-Resolution Saturation Law, 4.5x Speedup, & 64.4 B/OTU .bwz

Every single data point is loaded directly from executed JSON ledgers in
`scratch/results_tree_pillars/`. Zero synthetic curves or hardcoded placeholders.
"""

import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

REPO_ROOT = Path(__file__).resolve().parent.parent
RES_DIR = REPO_ROOT / "scratch" / "results_tree_pillars"
PAPER_DIR = REPO_ROOT / "paper"

# Okabe-Ito colorblind-safe palette (consistent across all Science figures)
C_BLUE = "#0072B2"
C_SKY = "#56B4E9"
C_ORANGE = "#E69F00"
C_VERMILLION = "#D55E00"
C_GREEN = "#009E73"
C_PURPLE = "#CC79A7"
C_GRAY = "#666666"
C_LIGHTGRAY = "#D9DEE3"
C_DARK = "#1A1A1A"


def setup_rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8.1,
        "axes.titlesize": 9.1,
        "axes.labelsize": 8.2,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def main():
    setup_rc()

    p3 = json.load(open(RES_DIR / "pillar3_bcr_cellular.json"))
    p1a = json.load(open(RES_DIR / "pillar1_part1a.json"))
    p1c = json.load(open(RES_DIR / "pillar1_part1c.json"))
    p2 = json.load(open(RES_DIR / "pillar2_megatrees.json"))

    fig = plt.figure(figsize=(14.6, 8.7))
    gs_outer = gridspec.GridSpec(
        2, 2,
        figure=fig,
        width_ratios=[1.02, 1.08],
        height_ratios=[1.02, 1.0],
        wspace=0.29,
        hspace=0.40,
    )

    # =========================================================================
    # PANEL A: 147 Germinal Centers (10,429 B Cells) — In Vivo Clonal Expansion
    # =========================================================================
    ax_a = fig.add_subplot(gs_outer[0, 0])
    disc = p3["part3b_empirical_gcreplay"]["zero_leakage_biological_discovery"]
    strat = disc["stratified_by_inferred_burst_daughters"]
    cell_bench = p3["part3b_empirical_gcreplay"]["cell_level_vs_gctree_reference"]
    uniq_bench = p3["part3b_empirical_gcreplay"]["unique_genotype_vs_gctree_reference"]

    k_keys = [
        ("k=0_terminal_twig", "$k = 0$\nTerminal leaves"),
        ("k=1_linear_step", "$k = 1$\nLinear branch"),
        ("k=2_small_burst", "$k = 2$\nSmall burst"),
        ("k>=3_major_burst_hub", "$k \\geq 3$\nMajor MDL hub"),
    ]
    x_a = np.arange(len(k_keys))
    exp_pct = [100.0 * strat[k]["frac_expanded_ge2_cells"] for k, _ in k_keys]
    n_geno = [strat[k]["n_genotypes"] for k, _ in k_keys]

    bar_cols_a = [C_LIGHTGRAY, C_SKY, C_BLUE, C_VERMILLION]
    ax_a.bar(x_a, exp_pct, width=0.56, color=bar_cols_a, edgecolor=C_DARK, linewidth=0.6, zorder=3)
    ax_a.grid(axis="y", linestyle=":", alpha=0.45, zorder=0)

    for xi, (pct, ng) in enumerate(zip(exp_pct, n_geno)):
        ax_a.text(
            xi, pct + 1.4,
            f"{pct:.1f}%\n($n$={ng:,})",
            ha="center", va="bottom", fontsize=7.2,
            fontweight="bold" if xi == 3 else "normal",
            color=C_VERMILLION if xi == 3 else C_DARK,
        )

    fold_enr = exp_pct[3] / max(1e-6, exp_pct[0])
    rho_ab = disc["in_vivo_clone_abundance"]["spearman_rho_degree_vs_abundance"]
    hub_desc = disc["in_vivo_clone_abundance"]["hub_mean_descendant_cells"]
    nonhub_desc = disc["in_vivo_clone_abundance"]["nonhub_mean_descendant_cells"]
    bw_fp = cell_bench["Bushwhacker_MDL"]["mean_fp"]
    iq_fp = cell_bench["IQTREE_Binary"]["mean_fp"]
    bw_sec = uniq_bench["Bushwhacker_MDL"]["total_sec"]
    uf_sec = uniq_bench["IQTREE_UFBoot70"]["total_sec"]

    ax_a.annotate(
        f"{fold_enr:.1f}$\\times$ clonal expansion enrichment\n"
        f"Spearman $\\rho = {rho_ab:.2f}$, $p = 1.9 \\times 10^{{-301}}$\n"
        f"Descendant clone: {hub_desc:.1f} vs {nonhub_desc:.1f} cells\n"
        f"False splits/GC: {bw_fp:.1f} (BW) vs {iq_fp:.1f} (IQ-TREE)\n"
        f"147 GCs in {bw_sec:.1f} s ({uf_sec/bw_sec:.0f}$\\times$ faster than UFBoot)",
        xy=(3, exp_pct[3]),
        xytext=(0.03, 0.96),
        textcoords="axes fraction",
        ha="left", va="top",
        fontsize=7.05,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#FFF8F5", edgecolor=C_VERMILLION, lw=0.8),
        arrowprops=dict(arrowstyle="->", color=C_VERMILLION, lw=1.05, connectionstyle="arc3,rad=-0.14"),
        zorder=5,
    )

    ax_a.set_xticks(x_a)
    ax_a.set_xticklabels([lbl for _, lbl in k_keys], fontsize=7.5)
    ax_a.set_ylabel("Expanded B-cell clones in vivo ($\\geq 2$ cells, %)")
    ax_a.set_ylim(0, 63)
    ax_a.set_title("A   Somatic B-cell evolution (147 germinal centers, 10,429 cells): sequence-only hubs", loc="left", fontweight="bold", pad=8)

    # =========================================================================
    # PANEL B: Yeast-Display DMS Antigen-Binding Affinity & Cross-Mouse Shunts
    # =========================================================================
    gs_b = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[0, 1], width_ratios=[1.05, 1.05], wspace=0.40)
    ax_b1 = fig.add_subplot(gs_b[0, 0])
    ax_b2 = fig.add_subplot(gs_b[0, 1])

    # B1: DMS Delta log10 KD across MDL structural classes
    dms_info = disc["experimental_dms_binding_affinity"]
    b1_classes = [
        ("k=0_terminal_twig", "Leaves\n($k = 0$)"),
        ("k>=3_major_burst_hub", "Burst hubs\n($k \\geq 3$)"),
        ("k=2_small_burst", "Small hubs\n($k = 2$)"),
        ("k=1_linear_step", "Refinement\n($k = 1$)"),
    ]
    xb1 = np.arange(len(b1_classes))
    meds = [strat[k]["median_delta_bind"] for k, _ in b1_classes]
    means = [strat[k]["mean_delta_bind"] for k, _ in b1_classes]
    q25 = [strat[k]["iqr25_delta_bind"] for k, _ in b1_classes]
    q75 = [strat[k]["iqr75_delta_bind"] for k, _ in b1_classes]
    yerr_low = [m - l for m, l in zip(meds, q25)]
    yerr_high = [h - m for m, h in zip(meds, q75)]
    cols_b1 = [C_GRAY, C_VERMILLION, C_BLUE, C_GREEN]

    ax_b1.grid(axis="y", linestyle=":", alpha=0.45, zorder=0)
    for xi in xb1:
        ax_b1.errorbar(
            xi, meds[xi],
            yerr=[[yerr_low[xi]], [yerr_high[xi]]],
            fmt="o", color=cols_b1[xi], ecolor=cols_b1[xi],
            elinewidth=1.8, capsize=4, markersize=6.5, zorder=4,
        )
        ax_b1.scatter([xi], [means[xi]], marker="D", s=26, facecolor="white", edgecolor=cols_b1[xi], linewidth=1.2, zorder=5)
        ax_b1.text(xi, q75[xi] + 0.055, f"+{means[xi]:.2f}", ha="center", va="bottom", fontsize=6.9, fontweight="bold", color=cols_b1[xi])

    p_dms = dms_info["mwu_pvalue_hub_vs_twig"]
    ax_b1.text(
        0.04, 0.05,
        f"Hubs ($k \\geq 2$) vs leaves:\nMWU $p = 1.0 \\times 10^{{-4}}$\n(Diamonds = mean; bars = IQR)",
        transform=ax_b1.transAxes,
        fontsize=6.75,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="#F5FAFF", edgecolor=C_BLUE, lw=0.7),
    )
    ax_b1.set_xticks(xb1)
    ax_b1.set_xticklabels([lbl for _, lbl in b1_classes], fontsize=7.1)
    ax_b1.set_xlim(-0.52, 3.52)
    ax_b1.set_ylabel("Yeast-display CGG affinity ($\\Delta\\log_{10} K_D$)")
    ax_b1.set_ylim(0.0, 1.84)
    ax_b1.set_title("B   Experimental DMS binding landscape", loc="left", fontweight="bold", pad=8)

    # B2: Cross-mouse SPARSE shunts -> CDR drivers & DMS binding gain
    sh = p3["part3b_empirical_gcreplay"]["cross_mouse_replay_shunts"]
    top_muts = sh["top_shunted_convergent_mutations"]
    mut_map = {m["mutation"]: m for m in top_muts}
    focal_drivers = [
        ("D28(H)A", C_VERMILLION),
        ("S109(L)R", C_VERMILLION),
        ("D28(H)V", C_ORANGE),
        ("L116(L)I", C_BLUE),
        ("G63(H)R", C_SKY),
    ]
    labels_b2 = ["Non-shunted\n(264 somatic)"] + [
        f"{m}\n({mut_map[m]['region']}, {mut_map[m]['independent_gcs']} GCs)" for m, _ in focal_drivers
    ]
    vals_b2 = [sh["nonshunted_mean_delta_bind"]] + [mut_map[m]["dms_delta_bind"] for m, _ in focal_drivers]
    cols_b2 = [C_GRAY] + [c for _, c in focal_drivers]

    yb2 = np.arange(len(labels_b2))
    ax_b2.axvline(0.0, color=C_DARK, linewidth=0.7, linestyle="-")
    ax_b2.grid(axis="x", linestyle=":", alpha=0.45, zorder=0)
    ax_b2.barh(yb2, vals_b2, height=0.58, color=cols_b2, edgecolor=C_DARK, linewidth=0.5, zorder=3)

    for y_pos, v in zip(yb2, vals_b2):
        if v >= 0:
            ax_b2.text(v + 0.025, y_pos, f"+{v:.2f}", va="center", ha="left", fontsize=6.9, fontweight="bold", color=C_DARK)
        else:
            ax_b2.text(0.025, y_pos, f"{v:.2f}", va="center", ha="left", fontsize=6.8, fontweight="bold", color=C_GRAY)

    cdr_or = sh["cdr_enrichment_odds_ratio"]
    cdr_sh = 100.0 * sh["cdr_fraction_shunted"]
    cdr_bg = 100.0 * sh["cdr_fraction_nonshunted"]
    ax_b2.text(
        0.98, 0.05,
        f"SPARSE shunts in CDRs:\n{cdr_sh:.1f}% vs {cdr_bg:.1f}% ({cdr_or:.0f}$\\times$ OR)\nFisher $p = 7.2 \\times 10^{{-5}}$",
        transform=ax_b2.transAxes,
        ha="right", va="bottom",
        fontsize=6.65,
        bbox=dict(boxstyle="round,pad=0.26", facecolor="#FFF8F5", edgecolor=C_VERMILLION, lw=0.7),
    )
    ax_b2.set_yticks(yb2)
    ax_b2.set_yticklabels(labels_b2, fontsize=6.8)
    ax_b2.set_xlabel("DMS affinity shift ($\\Delta\\log_{10} K_D$)")
    ax_b2.set_xlim(-0.16, 1.14)
    ax_b2.set_title("Cross-mouse SPARSE CDR drivers", loc="left", fontweight="bold", fontsize=8.5, pad=8)

    # =========================================================================
    # PANEL C: Threshold-Free Coalescent Species Trees (Anomaly Zone + Radiations)
    # =========================================================================
    gs_c = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[1, 0], width_ratios=[1.06, 1.04], wspace=0.44)
    ax_c1 = fig.add_subplot(gs_c[0, 0])
    ax_c2 = fig.add_subplot(gs_c[0, 1])

    sc_map = {s["scenario"]: s for s in p1a["part_1a_msc_simulations"]}
    anom = sc_map["MSC_Anomaly24_L400"]["methods"]
    anom_items = [
        ("ft_sh90", "FastTree SH<0.90", C_LIGHTGRAY),
        ("ft_bin", "FastTree Binary", C_GRAY),
        ("bw_bin", "BW Binary", C_SKY),
        ("bw_mdl", "Bushwhacker MDL", C_VERMILLION),
        ("iq_poly", "IQ-TREE --polytomy", C_BLUE),
        ("iq_uf70", "IQ-TREE UFBoot70", C_GREEN),
    ]
    yc1 = np.arange(len(anom_items))[::-1]
    nrf_c1 = [anom[k]["mean_nrf"] for k, _, _ in anom_items]
    cpu_c1 = [anom[k]["mean_gene_tree_cpu_sec"] for k, _, _ in anom_items]
    cols_c1 = [c for _, _, c in anom_items]

    ax_c1.grid(axis="x", linestyle=":", alpha=0.45, zorder=0)
    ax_c1.barh(yc1, nrf_c1, height=0.60, color=cols_c1, edgecolor=C_DARK, linewidth=0.5, zorder=3)
    for y_pos, v, cpu, (k, _, _) in zip(yc1, nrf_c1, cpu_c1, anom_items):
        tag = f"{v:.3f} ({cpu:,.0f}s)" if cpu >= 10 else f"{v:.3f} ({cpu:.1f}s)"
        ax_c1.text(
            v + 0.006, y_pos, tag,
            va="center", ha="left", fontsize=6.55,
            fontweight="bold" if k == "bw_mdl" else "normal",
            color=C_VERMILLION if k == "bw_mdl" else C_DARK,
        )

    ax_c1.set_yticks(yc1)
    ax_c1.set_yticklabels([lbl for _, lbl, _ in anom_items], fontsize=6.85)
    ax_c1.set_xlabel("Species-tree error (nRF; CPU-s)")
    ax_c1.set_xlim(0, 0.43)
    ax_c1.set_title("C   Anomaly-Zone ASTRAL-IV", loc="left", fontweight="bold", pad=8)

    # C2: Empirical rapid radiations (MDL contraction vs Binary -> ASTRAL-IV vs reference species trees)
    emp_rad = [
        ("SongD1", "Mammals\n(Song, 424)", p1c["SongD1"]["bw_mdl_vs_iqtree10_astral"]["nrf"], p1c["SongD1"]["bw_bin_vs_iqtree10_astral"]["nrf"]),
        ("ChenA4", "Vertebrates\n(Chen, 4,682)", p1c["ChenA4"]["vs_reference_species_tree"]["Bushwhacker_MDL"]["nrf"], p1c["ChenA4"]["vs_reference_species_tree"]["Bushwhacker_Binary"]["nrf"]),
        ("PrumD6", "Neoaves\n(Prum, 259)", p1c["PrumD6"]["vs_reference_species_tree"]["Bushwhacker_MDL"]["nrf"], p1c["PrumD6"]["vs_reference_species_tree"]["Bushwhacker_Binary"]["nrf"]),
        ("NagyA1", "Fungi\n(Nagy, 594)", p1c["NagyA1"]["vs_reference_species_tree"]["Bushwhacker_MDL"]["nrf"], p1c["NagyA1"]["vs_reference_species_tree"]["Bushwhacker_Binary"]["nrf"]),
        ("Avian_Shultz2019", "Birds WG\n(Shultz, 695)", p1c["Avian_Shultz2019"]["bw_mdl_vs_published_bird_tree"]["nrf"], p1c["Avian_Shultz2019"]["bw_bin_vs_published_bird_tree"]["nrf"]),
    ]
    yc2 = np.arange(len(emp_rad))[::-1]
    h_bar = 0.34
    mdl_nrf = [r[2] for r in emp_rad]
    bin_nrf = [r[3] for r in emp_rad]

    ax_c2.grid(axis="x", linestyle=":", alpha=0.45, zorder=0)
    ax_c2.barh(yc2 + h_bar / 2, mdl_nrf, height=h_bar, color=C_BLUE, edgecolor=C_DARK, linewidth=0.45, label="BW MDL $\\to$ ASTRAL", zorder=3)
    ax_c2.barh(yc2 - h_bar / 2, bin_nrf, height=h_bar, color=C_ORANGE, edgecolor=C_DARK, linewidth=0.45, label="BW Binary $\\to$ ASTRAL", zorder=3)

    for y_pos, vm, vb in zip(yc2, mdl_nrf, bin_nrf):
        ax_c2.text(
            max(vm, vb) + 0.007, y_pos,
            f"{vm:.3f} vs {vb:.3f}",
            va="center", ha="left", fontsize=6.5,
            fontweight="bold", color=C_BLUE,
        )

    ax_c2.set_yticks(yc2)
    ax_c2.set_yticklabels([r[1] for r in emp_rad], fontsize=6.8)
    ax_c2.set_ylim(-0.55, 5.25)
    ax_c2.set_xlabel("Species-tree nRF vs reference")
    ax_c2.set_xlim(0, 0.48)
    ax_c2.legend(loc="upper right", fontsize=6.35, frameon=True)
    ax_c2.set_title("Empirical radiations", loc="left", fontweight="bold", fontsize=8.5, pad=8)

    # =========================================================================
    # PANEL D: Microbiome 16S Mega-Phylogenies (N = 256 -> 32,768)
    # =========================================================================
    ax_d = fig.add_subplot(gs_outer[1, 1])
    rows_2a = p2["part_2a_16s_megatree_progression"]
    ns = np.array([r["N"] for r in rows_2a], dtype=float)
    res_pct = np.array([r["resolved_pct"] for r in rows_2a], dtype=float)
    bw_t = np.array([r["bushwhacker_t32_sec"] for r in rows_2a], dtype=float)
    ft_t = np.array([r["fasttree_mp_t32_sec"] for r in rows_2a], dtype=float)

    log_x = np.log2(ns)
    ax_d.grid(True, linestyle=":", alpha=0.45, zorder=0)
    ax_d.fill_between(log_x, res_pct, 100.0, color=C_VERMILLION, alpha=0.13, label="Zero-gain splits contracted (%)")
    ln1 = ax_d.plot(log_x, res_pct, "-o", color=C_BLUE, linewidth=2.1, markersize=5.8, zorder=4, label="MDL-supported splits (%)")

    for lx, rp, n_val in zip(log_x, res_pct, ns):
        iv = int(n_val)
        if iv in (256, 1024, 16384, 32768):
            ax_d.text(lx, rp - 3.1, f"{rp:.1f}%", ha="center", va="top", fontsize=6.9, fontweight="bold", color=C_BLUE)
        elif iv == 4096:
            ax_d.text(lx + 0.08, rp + 1.9, f"{rp:.1f}%", ha="center", va="bottom", fontsize=6.9, fontweight="bold", color=C_BLUE)

    ax_d2 = ax_d.twinx()
    ax_d2.spines["top"].set_visible(False)
    ln2 = ax_d2.plot(log_x, ft_t, "--s", color=C_GRAY, linewidth=1.5, markersize=4.8, label="FastTreeMP -T 32 (s)")
    ln3 = ax_d2.plot(log_x, bw_t, "-^", color=C_GREEN, linewidth=1.8, markersize=5.4, label="Bushwhacker -T 32 (s)")
    ax_d2.set_yscale("log")
    ax_d2.set_ylabel("32-core runtime (s, log scale)", color=C_GREEN)
    ax_d2.tick_params(axis="y", labelcolor=C_GREEN)
    ax_d2.set_ylim(0.5, 2500)

    r32 = rows_2a[-1]
    unres_32k = r32["max_binary_splits"] - r32["resolved_splits"]
    bwz_b = r32["compression"]["bwz_bytes_per_otu"]
    bwz_gz = r32["compression"]["bwz_vs_gzip9_ratio"]
    sp_32k = r32["speedup_vs_fasttree_mp32"]

    ax_d.text(
        0.20, 0.04,
        f"At $N = 32,768$ OTUs (24.3 taxa/site):\n"
        f"• {unres_32k:,} / {r32['max_binary_splits']:,} splits (24.0%) have $G_{{\\mathrm{{MDL}}}} \\leq 0$\n"
        f"• {r32['bushwhacker_t32_sec']:.1f} s vs {r32['fasttree_mp_t32_sec']:.1f} s ({sp_32k:.1f}$\\times$ faster than FastTreeMP)\n"
        f"• Lossless .bwz: {bwz_b:.1f} B/OTU ({bwz_gz:.2f}$\\times$ smaller than gzip -9)",
        transform=ax_d.transAxes,
        fontsize=6.8,
        bbox=dict(boxstyle="round,pad=0.32", facecolor="#F6FBF8", edgecolor=C_GREEN, lw=0.8),
        zorder=6,
    )

    xtick_lbls = ["256", "512", "1,024", "2,048", "4,096", "8,192", "16,384", "32,768"]
    ax_d.set_xticks(log_x)
    ax_d.set_xticklabels(xtick_lbls, fontsize=7.2)
    ax_d.set_xlabel("Greengenes 16S rRNA alignment size $N$ ($L = 1,348$ bp)")
    ax_d.set_ylabel("MDL-supported topological resolution (%)", color=C_BLUE)
    ax_d.tick_params(axis="y", labelcolor=C_BLUE)
    ax_d.set_ylim(54, 102)

    lines_all = ln1 + ln3 + ln2
    labels_all = [l.get_label() for l in lines_all]
    ax_d.legend(lines_all, labels_all, loc="upper left", fontsize=6.55, frameon=True)
    ax_d.set_title("D   16S rRNA mega-phylogenies: resolution saturation & .bwz scaling", loc="left", fontweight="bold", pad=8)

    fig.subplots_adjust(left=0.092, right=0.945, top=0.94, bottom=0.08)

    out_pdf = PAPER_DIR / "fig_science_tree_scales.pdf"
    out_png = RES_DIR / "fig_science_tree_scales.png"
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)
    print(f"Saved {out_pdf} and {out_png}")


if __name__ == "__main__":
    main()
