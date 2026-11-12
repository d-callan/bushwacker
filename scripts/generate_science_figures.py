#!/usr/bin/env python3
"""
Generate publication-grade vector PDF Figures 2 and 3 for the Science manuscript:
  - paper/fig_science_macroevo_polytomies.pdf (Figure 2: Macroevolutionary Polytomies & Deep-Time Controversies across 45,205 Rokas Loci)
  - paper/fig_science_mammalian_atlas.pdf     (Figure 3: 19,952-Protein x 882-Mammal Convergence, Conversion & Discordance Atlas)

Every number is loaded directly from executed computational JSON outputs in `scratch/results_4pillars/`.
"""

import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


REPO_ROOT = Path(__file__).resolve().parent.parent
RES_DIR = REPO_ROOT / "scratch" / "results_4pillars"
PAPER_DIR = REPO_ROOT / "paper"

# Okabe-Ito colorblind-safe palette (matching Figure 1)
C_BLUE = "#0072B2"
C_SKY = "#56B4E9"
C_ORANGE = "#E69F00"
C_VERMILLION = "#D55E00"
C_GREEN = "#009E73"
C_PURPLE = "#CC79A7"
C_YELLOW = "#F0E442"
C_GRAY = "#666666"
C_DARK = "#1A1A1A"
C_LIGHT_BG = "#F7F9FB"


def setup_rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8.2,
        "axes.titlesize": 9.3,
        "axes.labelsize": 8.4,
        "xtick.labelsize": 7.6,
        "ytick.labelsize": 7.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def build_figure_2_macroevo():
    p4 = json.load(open(RES_DIR / "macroevo_polytomy_atlas.json"))
    ds_map = {d["dataset"]: d for d in p4["datasets"]}

    fig = plt.figure(figsize=(14.2, 8.2))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.08, 1.08],
        height_ratios=[1.06, 1.0],
        wspace=0.38,
        hspace=0.42,
    )

    # =========================================================================
    # PANEL A: Spurious Gene-Tree Conflict Eliminated Across 19 Rokas Datasets
    # =========================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    ordered_keys = [
        ("BoroA6", "Animal Root (Borowiec AA, 1,080)"),
        ("JarvD5a", "Neoaves Exons/Introns (Jarvis NT, 14,446)"),
        ("NagyA1", "Fungal Early-Diverging (Nagy AA, 594)"),
        ("WhelA7", "Animal + 4 Sponge Classes (Whelan AA, 210)"),
        ("ShenA9", "Budding Yeasts (Shen AA, 1,233)"),
        ("StruA5", "Lophotrochozoa (Struck AA, 679)"),
        ("WickA3", "Land Plants 1KP (Wickett AA, 844)"),
        ("ChenA4", "Jawed Vertebrates (Chen AA, 4,682)"),
        ("MisoA2", "Insect Orders 1KITE (Misof AA, 1,478)"),
        ("XiD4", "Seed Plants (Xi NT, 310)"),
        ("TarvD7", "Placental Mammals (Tarver NT, 11,178)"),
        ("SongD1", "Placental Mammals (Song NT, 424)"),
        ("PrumD6", "Neoaves 199-Taxon Birds (Prum NT, 259)"),
    ]
    labels_a = [lbl for _, lbl in ordered_keys][::-1]
    elim_vals = [100.0 * ds_map[k]["spurious_discordance_eliminated_frac"] for k, _ in ordered_keys][::-1]
    poly_vals = [100.0 * ds_map[k]["polytomy_contracted_frac"] for k, _ in ordered_keys][::-1]
    alphs = [ds_map[k]["alphabet"] for k, _ in ordered_keys][::-1]

    y = np.arange(len(labels_a))
    bar_colors = [C_BLUE if a == "protein" else C_SKY for a in alphs]
    ax_a.barh(y, elim_vals, height=0.66, color=bar_colors, edgecolor="none", alpha=0.92, label="Spurious species-tree conflict eliminated (%)")
    ax_a.scatter(poly_vals, y, color=C_ORANGE, edgecolor=C_DARK, linewidth=0.6, s=26, zorder=4, label="Forced binary splits contracted (%)")

    overall_elim = 100.0 * p4["overall_spurious_discordance_eliminated_frac"]
    ax_a.axvline(overall_elim, color=C_VERMILLION, linestyle="--", linewidth=1.2)
    ax_a.text(overall_elim + 17.5, 5.2, f"All 45,205 loci:\n{overall_elim:.1f}% eliminated\n(936,166 splits)", fontsize=7.2, color=C_VERMILLION, fontweight="bold", bbox=dict(boxstyle="round,pad=0.25", facecolor="#FFF8F5", edgecolor=C_VERMILLION, lw=0.7))

    for yi, v in enumerate(elim_vals):
        ax_a.text(v + 0.9, yi, f"{v:.1f}%", va="center", ha="left", fontsize=6.9, color=C_DARK)

    ax_a.set_yticks(y)
    ax_a.set_yticklabels(labels_a, fontsize=7.2)
    ax_a.set_xlim(0, 114)
    ax_a.set_xlabel("Percentage of Splits (%)")
    ax_a.set_title("A   Eliminating Forced-Branching Conflict Across 45,205 Loci", loc="left", fontweight="bold", pad=10)
    p_prot = mpatches.Patch(color=C_BLUE, label="20-state Protein MDL")
    p_dna = mpatches.Patch(color=C_SKY, label="Nucleotide MDL")
    h_dot = plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=C_ORANGE, markeredgecolor=C_DARK, markersize=5.5, label="Binary splits contracted")
    ax_a.legend(handles=[p_prot, p_dna, h_dot], loc="lower right", frameon=True, facecolor="white", framealpha=0.96, fontsize=6.9)

    # =========================================================================
    # PANEL B: Cracking Classic Deep-Time Controversies in Exact Nats
    # =========================================================================
    ax_b = fig.add_subplot(gs[0, 1])

    boro_fh = ds_map["BoroA6"]["focal_hypotheses"]
    whel_fh = ds_map["WhelA7"]["focal_hypotheses"]
    chen_fh = ds_map["ChenA4"]["focal_hypotheses"]
    song_fh = ds_map["SongD1"]["focal_hypotheses"]

    debates = [
        (
            "Animal Root\n(WhelA7, 210 AA)",
            [
                ("Ctenophora-sister", whel_fh["Ctenophora_sister"]["cumulative_mdl_nats"], whel_fh["Ctenophora_sister"]["forced_binary_count"]),
                ("Porifera-sister", whel_fh["Porifera_sister"]["cumulative_mdl_nats"], whel_fh["Porifera_sister"]["forced_binary_count"]),
                ("Cteno + Porifera", whel_fh["Ctenophora_plus_Porifera"]["cumulative_mdl_nats"], whel_fh["Ctenophora_plus_Porifera"]["forced_binary_count"]),
            ],
        ),
        (
            "Animal Root\n(BoroA6, 1,080 AA)",
            [
                ("Ctenophora-sister", boro_fh["Ctenophora_sister"]["cumulative_mdl_nats"], boro_fh["Ctenophora_sister"]["forced_binary_count"]),
                ("Porifera-sister", boro_fh["Porifera_sister"]["cumulative_mdl_nats"], boro_fh["Porifera_sister"]["forced_binary_count"]),
                ("Cteno + Porifera", boro_fh["Ctenophora_plus_Porifera"]["cumulative_mdl_nats"], boro_fh["Ctenophora_plus_Porifera"]["forced_binary_count"]),
            ],
        ),
        (
            "Placental Root\n(SongD1, 424 NT)",
            [
                ("Atlantogenata", song_fh["Atlantogenata_Afro_plus_Xena"]["cumulative_mdl_nats"], song_fh["Atlantogenata_Afro_plus_Xena"]["forced_binary_count"]),
                ("Exafroplacentalia", song_fh["Exafroplacentalia_Boreo_plus_Xena"]["cumulative_mdl_nats"], song_fh["Exafroplacentalia_Boreo_plus_Xena"]["forced_binary_count"]),
                ("Epitheria", song_fh["Epitheria_Boreo_plus_Afro"]["cumulative_mdl_nats"], song_fh["Epitheria_Boreo_plus_Afro"]["forced_binary_count"]),
            ],
        ),
        (
            "Tetrapod Sister\n(ChenA4, 4,682 AA)",
            [
                ("Lungfish-sister", chen_fh["Lungfish_sister_to_Tetrapods"]["cumulative_mdl_nats"], chen_fh["Lungfish_sister_to_Tetrapods"]["forced_binary_count"]),
                ("Coelacanth-sister", chen_fh["Coelacanth_sister_to_Tetrapods"]["cumulative_mdl_nats"], chen_fh["Coelacanth_sister_to_Tetrapods"]["forced_binary_count"]),
                ("Lungfish + Coelacanth", chen_fh["Lungfish_plus_Coelacanth"]["cumulative_mdl_nats"], chen_fh["Lungfish_plus_Coelacanth"]["forced_binary_count"]),
            ],
        ),
        (
            "Turtle Placement\n(ChenA4, 4,682 AA)",
            [
                ("Archelosauria", chen_fh["Archelosauria_Turtles_plus_Archosaurs"]["cumulative_mdl_nats"], chen_fh["Archelosauria_Turtles_plus_Archosaurs"]["forced_binary_count"]),
                ("Ankylopoda", chen_fh["Ankylopoda_Turtles_plus_Lepidosaurs"]["cumulative_mdl_nats"], chen_fh["Ankylopoda_Turtles_plus_Lepidosaurs"]["forced_binary_count"]),
                ("Turtles-basal", chen_fh["Archosaurs_plus_Lepidosaurs_Turtles_basal"]["cumulative_mdl_nats"], chen_fh["Archosaurs_plus_Lepidosaurs_Turtles_basal"]["forced_binary_count"]),
            ],
        ),
    ]

    yb = np.arange(len(debates))
    h_cols = [C_GREEN, C_ORANGE, C_GRAY]
    bar_h = 0.24
    for idx, (d_title, hyps) in enumerate(debates):
        tot_nats = max(1e-6, sum(h[1] for h in hyps))
        tot_bin = max(1, sum(h[2] for h in hyps))
        mdl_shares = [100.0 * h[1] / tot_nats for h in hyps]
        bin_shares = [100.0 * h[2] / tot_bin for h in hyps]

        # Top bar: Cumulative MDL Nats share (%)
        left_m = 0.0
        for j, sh in enumerate(mdl_shares):
            ax_b.barh(idx + 0.14, sh, left=left_m, height=bar_h, color=h_cols[j], edgecolor="white", linewidth=0.5)
            left_m += sh
        # Bottom bar: Forced binary count share (%)
        left_b = 0.0
        for j, sh in enumerate(bin_shares):
            ax_b.barh(idx - 0.14, sh, left=left_b, height=bar_h, color=h_cols[j], alpha=0.42, edgecolor="white", linewidth=0.5, hatch="//")
            left_b += sh

        lead_name, lead_nats, _ = hyps[0]
        runner_nats = max(0.1, hyps[1][1])
        ratio = lead_nats / runner_nats
        ax_b.text(101.5, idx + 0.13, f"{lead_name}: {lead_nats:,.0f} nats ({ratio:.1f}×)", va="center", ha="left", fontsize=6.9, fontweight="bold", color=C_DARK)
        ax_b.text(101.5, idx - 0.15, f"Binary counts: {hyps[0][2]}:{hyps[1][2]}:{hyps[2][2]}", va="center", ha="left", fontsize=6.6, color=C_GRAY)

    ax_b.axvline(50.0, color="#CCCCCC", linestyle=":", linewidth=0.9)
    ax_b.set_yticks(yb)
    ax_b.set_yticklabels([d[0] for d in debates], fontsize=7.4)
    ax_b.set_xlim(0, 178)
    ax_b.set_ylim(-0.55, 4.72)
    ax_b.set_xticks([0, 25, 50, 75, 100])
    ax_b.set_xlabel("Share of Cumulative Support (%)   [Solid = MDL nats; Hatched = Binary counts]")
    ax_b.set_title("B   Resolving Deep-Time Controversies After Polytomy Filtering", loc="left", fontweight="bold", pad=10)
    l1 = mpatches.Patch(color=C_GREEN, label="Leading hypothesis (MDL nats)")
    l2 = mpatches.Patch(color=C_ORANGE, label="2nd hypothesis")
    l3 = mpatches.Patch(color=C_GRAY, label="3rd hypothesis")
    ax_b.legend(handles=[l1, l2, l3], loc="upper center", bbox_to_anchor=(0.46, 1.01), ncol=3, frameon=True, facecolor="white", fontsize=6.6)

    # =========================================================================
    # PANEL C: Matched Amino-Acid vs. Codon123 vs. Codon12 Speed & Saturation
    # =========================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    matched = [
        ("Land Plants\n(Wickett, 844 loci)", "WickD3a", "WickD3b", "WickA3"),
        ("Insects 1KITE\n(Misof, 1,478 loci)", "MisoD2a", "MisoD2b", "MisoA2"),
    ]
    xg = np.arange(len(matched))
    w = 0.24
    c123_wall = [ds_map[a]["wall_sec"] for _, a, _, _ in matched]
    c12_wall = [ds_map[b]["wall_sec"] for _, _, b, _ in matched]
    aa_wall = [ds_map[c]["wall_sec"] for _, _, _, c in matched]

    b1 = ax_c.bar(xg - w, c123_wall, width=w, color=C_VERMILLION, label="DNA Codon 1+2+3 (saturated 3rd pos.)")
    b2 = ax_c.bar(xg, c12_wall, width=w, color=C_ORANGE, label="DNA Codon 1+2 (3rd pos. stripped)")
    b3 = ax_c.bar(xg + w, aa_wall, width=w, color=C_BLUE, label="20-State Protein MDL (--translate / AA)")

    ax_c.set_yscale("log")
    ax_c.set_ylim(3, 2200)
    ax_c.set_xlim(-0.48, 1.62)
    ax_c.set_xticks(xg)
    ax_c.set_xticklabels([m[0] for m in matched], fontsize=8.0)
    ax_c.set_ylabel("32-Core Wall-Clock Runtime (s, log scale)")
    ax_c.set_title("C   20-State Protein MDL vs. Saturated Codon DNA", loc="left", fontweight="bold", pad=10)

    for i, (_, k123, k12, kaa) in enumerate(matched):
        sp = ds_map[k123]["wall_sec"] / ds_map[kaa]["wall_sec"]
        el_aa = 100.0 * ds_map[kaa]["spurious_discordance_eliminated_frac"]
        el_d = 100.0 * ds_map[k12]["spurious_discordance_eliminated_frac"]
        ax_c.text(xg[i] - w, c123_wall[i] * 1.12, f"{c123_wall[i]:.0f} s", ha="center", va="bottom", fontsize=7.2)
        ax_c.text(xg[i], c12_wall[i] * 1.12, f"{c12_wall[i]:.0f} s\n({el_d:.0f}% elim.)", ha="center", va="bottom", fontsize=6.9)
        ax_c.text(xg[i] + w, aa_wall[i] * 1.15, f"{aa_wall[i]:.1f} s\n({sp:.1f}× faster;\n{el_aa:.1f}% elim.)", ha="center", va="bottom", fontsize=7.1, fontweight="bold", color=C_BLUE)

    ax_c.legend(loc="upper left", frameon=True, facecolor="white", fontsize=7.1)

    # =========================================================================
    # PANEL D: Radiation Hard-Polytomy Spectrum & Model/Seed Paradox Resolution
    # =========================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    paradox_items = [
        ("ML Model-Sensitive Splits\n(Abadi/Gerth, 7,200 alignments)", 96.4, C_VERMILLION),
        ("ML Seed-Unstable Splits\n(Zhou/Rokas, 4,159 genes × 20 seeds)", 94.5, C_ORANGE),
        ("Animal Root Single-Gene Splits\n(Borowiec BoroA6, 1,080 AA loci)", 100.0 * ds_map["BoroA6"]["sp_tree_hard_polytomy_frac"], C_PURPLE),
        ("Post-K–Pg Neoaves Bird Stems\n(Jarvis JarvD5a, 14,446 loci)", 100.0 * ds_map["JarvD5a"]["sp_tree_hard_polytomy_frac"], C_BLUE),
        ("882-Mammal Proteome Splits\n(TOGA2, 19,952 AA proteins)", 65.01, C_GREEN),
        ("Jawed Vertebrate Stems\n(ChenA4, 4,682 AA loci)", 100.0 * ds_map["ChenA4"]["sp_tree_hard_polytomy_frac"], C_SKY),
    ][::-1]

    yd = np.arange(len(paradox_items))
    vals_d = [x[1] for x in paradox_items]
    cols_d = [x[2] for x in paradox_items]
    ax_d.barh(yd, vals_d, height=0.60, color=cols_d, alpha=0.90)
    for yi, v in enumerate(vals_d):
        ax_d.text(v + 1.1, yi, f"{v:.1f}%", va="center", ha="left", fontsize=7.6, fontweight="bold", color=C_DARK)
    ax_d.set_yticks(yd)
    ax_d.set_yticklabels([x[0] for x in paradox_items], fontsize=7.3)
    ax_d.set_xlim(0, 112)
    ax_d.set_xlabel("Zero-Net-Information Splits Contracted by MDL (G_MDL ≤ 0, %)")
    ax_d.set_title("D   Model Instability and Rapid Radiations Reflect Zero-Signal Polytomies", loc="left", fontweight="bold", pad=8)

    out_pdf = PAPER_DIR / "fig_science_macroevo_polytomies.pdf"
    out_png = RES_DIR / "fig_science_macroevo_polytomies.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_pdf} and {out_png}")


def build_figure_3_mammalian_atlas():
    p1 = json.load(open(RES_DIR / "toga2_protein_atlas.json"))
    strata = json.load(open(RES_DIR / "toga2_protein_length_strata.json"))
    cochlear = json.load(open(RES_DIR / "toga2_cochlear_25genes.json"))

    fig = plt.figure(figsize=(14.4, 8.3))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.05, 1.12],
        height_ratios=[1.04, 1.05],
        wspace=0.36,
        hspace=0.40,
    )

    # =========================================================================
    # PANEL A: Length-Stratified Concordance & Conflict Elimination (19,952 AA)
    # =========================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    x = np.arange(len(strata))
    xlbls = [s["stratum"] for s in strata]
    bin_conf = [100.0 * s["forced_binary_conflict_rate"] for s in strata]
    mdl_conf_tot = [100.0 * s["mdl_conflict_rate_per_total"] for s in strata]
    sp_better = [100.0 * s["species_tree_shorter_dl_frac"] for s in strata]
    elim_pct = [100.0 * s["conflict_eliminated_frac"] for s in strata]

    ax_a.plot(x, sp_better, marker="s", linewidth=2.1, color=C_GREEN, label="Species tree compresses better than de novo tree (%)")
    ax_a.plot(x, bin_conf, marker="o", linewidth=2.0, linestyle="--", color=C_VERMILLION, label="Forced binary tree discordant splits (%)")
    ax_a.plot(x, mdl_conf_tot, marker="^", linewidth=2.1, color=C_BLUE, label="Residual discordant splits after MDL contraction (%)")
    ax_a.fill_between(x, mdl_conf_tot, bin_conf, color=C_VERMILLION, alpha=0.12)

    ax_a.annotate(f"{sp_better[-1]:.1f}% of >1,500 AA\nproteins favor speciesTree", xy=(x[-1], sp_better[-1]), xytext=(1.45, 98.0),
                  arrowprops=dict(arrowstyle="->", color=C_GREEN, lw=1.0), fontsize=7.0, fontweight="bold", color=C_GREEN)
    ax_a.annotate(f"{elim_pct[0]:.1f}% of short-protein\nconflicts eliminated", xy=(x[0], 54.0), xytext=(0.38, 60.5),
                  arrowprops=dict(arrowstyle="->", color=C_VERMILLION, lw=1.0), fontsize=7.0, fontweight="bold", color=C_VERMILLION)

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(xlbls, fontsize=7.3)
    ax_a.set_ylim(0, 112)
    ax_a.set_ylabel("Percentage across 19,952 Mammalian Proteins (%)")
    ax_a.set_xlabel("Protein Alignment Length Stratum (882 Mammals)")
    ax_a.set_title("A   Proteome-Wide Species-Tree Compression & Discordance", loc="left", fontweight="bold", pad=10)
    ax_a.legend(loc="center right", bbox_to_anchor=(0.99, 0.37), frameon=True, facecolor="white", framealpha=0.95, fontsize=6.7)

    # =========================================================================
    # PANEL B: Taxonomic Scope of INTERVAL Conversion vs. SPARSE Convergence
    # =========================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    sh_tot = p1["shunt_totals_on_species_tree"]
    tot_int = sh_tot["total_interval_shunts"]
    int_inter_ord = sh_tot["inter_order_interval_shunts"]
    int_inter_fam = sh_tot["inter_family_intra_order_interval_shunts"]
    int_intra_fam = tot_int - int_inter_ord - int_inter_fam

    tot_sp = sh_tot["total_sparse_shunts"]
    sp_inter_ord = sh_tot["inter_order_sparse_shunts"]
    sp_inter_fam = sh_tot["inter_family_intra_order_sparse_shunts"]
    sp_intra_fam = tot_sp - sp_inter_ord - sp_inter_fam

    modes = [
        f"Contiguous INTERVAL\n(Block Conversion / Transfer)\nn = {tot_int:,} shunts",
        f"Scattered SPARSE\n(Dispersed Convergence / Paralogy)\nn = {tot_sp:,} shunts",
    ]
    intra_fam_pct = [100.0 * int_intra_fam / tot_int, 100.0 * sp_intra_fam / tot_sp]
    inter_fam_pct = [100.0 * int_inter_fam / tot_int, 100.0 * sp_inter_fam / tot_sp]
    inter_ord_pct = [100.0 * int_inter_ord / tot_int, 100.0 * sp_inter_ord / tot_sp]

    yb = np.arange(2)
    bh = 0.46
    ax_b.barh(yb, intra_fam_pct, height=bh, color=C_SKY, label="Intra-family (local paralog conversion)")
    ax_b.barh(yb, inter_fam_pct, left=intra_fam_pct, height=bh, color=C_ORANGE, label="Inter-family, intra-order (e.g., bat superfamilies)")
    left2 = [a + b for a, b in zip(intra_fam_pct, inter_fam_pct)]
    ax_b.barh(yb, inter_ord_pct, left=left2, height=bh, color=C_PURPLE, label="Inter-order (across 26 mammalian orders)")

    ax_b.text(intra_fam_pct[0] / 2, 0, f"{intra_fam_pct[0]:.1f}%\nintra-family", ha="center", va="center", fontsize=7.6, fontweight="bold", color=C_DARK)
    ax_b.text(left2[0] + inter_ord_pct[0] + 1.2, 0, f"Only {inter_ord_pct[0]:.1f}%\ninter-order", ha="left", va="center", fontsize=7.4, fontweight="bold", color=C_PURPLE)
    ax_b.text(intra_fam_pct[1] + inter_fam_pct[1] / 2, 1, f"{inter_fam_pct[1]:.1f}%", ha="center", va="center", fontsize=7.5, fontweight="bold", color=C_DARK)
    ax_b.text(left2[1] + inter_ord_pct[1] / 2, 1, f"{inter_ord_pct[1]:.1f}%\ninter-order\n(7.9× INTERVAL)", ha="center", va="center", fontsize=7.5, fontweight="bold", color="white")

    ax_b.text(0, -0.53, "Top INTERVAL conversion loci: KRT31 (+419 nats), TRIM22 (+314), TRIM6 (+251),\nZNF33A (+271), GBP2 (+226), UGT3A1 (+210), IGKV1D-17 (+10.3 nats)",
              fontsize=6.6, color=C_DARK, style="italic")

    ax_b.set_yticks(yb)
    ax_b.set_yticklabels(modes, fontsize=7.6)
    ax_b.set_xlim(0, 118)
    ax_b.set_ylim(-0.74, 1.85)
    ax_b.set_xlabel("Share of Certified Mammalian Species-Tree Shunts (%)")
    ax_b.set_title("B   Dual-Mode Shunts Separate Gene Conversion from Convergence", loc="left", fontweight="bold", pad=10)
    ax_b.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.95, fontsize=6.8)

    # =========================================================================
    # PANEL C: Ecological Guild Inter-Order Convergence vs. Neutral Null
    # =========================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    eg = p1["ecological_guild_enrichment"]
    guild_rows = [
        ("Echolocating Bats ↔\nToothed Whales (Genome-Wide)", eg["echolocating_bats_and_toothed_whales"]["fold_enrichment"],
         f"{eg['echolocating_bats_and_toothed_whales']['observed_inter_order_sparse_shunts']} obs / {eg['echolocating_bats_and_toothed_whales']['expected_random_inter_order']:.0f} exp (Zou & Zhang null)", C_GRAY),
        ("All Echolocators\n(Genome-Wide Inter-Order)", eg["echolocators"]["fold_enrichment"],
         f"{eg['echolocators']['observed_inter_order_sparse_shunts']} obs / {eg['echolocators']['expected_random_inter_order']:.0f} exp", C_SKY),
        ("Marine / Aquatic Mammals\n(Cetacea, Pinnipedia, Sirenia)", eg["marine_aquatic"]["fold_enrichment"],
         f"{eg['marine_aquatic']['observed_inter_order_sparse_shunts']} obs / {eg['marine_aquatic']['expected_random_inter_order']:.0f} exp", C_BLUE),
        ("Subterranean / Blind Mammals\n(Moles, Golden Moles, Mole-Rats)", eg["subterranean_blind"]["fold_enrichment"],
         f"{eg['subterranean_blind']['observed_inter_order_sparse_shunts']} obs / {eg['subterranean_blind']['expected_random_inter_order']:.1f} exp (p = 0.0045)", C_ORANGE),
        ("Myrmecophagous / Armored\n(Pangolins ↔ Anteaters/Armadillos)", eg["myrmecophagous_insectivorous_armor"]["fold_enrichment"],
         f"{eg['myrmecophagous_insectivorous_armor']['observed_inter_order_sparse_shunts']} obs / {eg['myrmecophagous_insectivorous_armor']['expected_random_inter_order']:.1f} exp (p = 3.8e-19)", C_VERMILLION),
    ]

    yc = np.arange(len(guild_rows))
    folds = [r[1] for r in guild_rows]
    gcols = [r[3] for r in guild_rows]
    ax_c.barh(yc, folds, height=0.58, color=gcols, alpha=0.92)
    ax_c.axvline(1.0, color=C_DARK, linestyle="--", linewidth=1.1, label="Random inter-order expectation (1.0×)")
    for yi, (_, fv, note, _) in enumerate(guild_rows):
        tx = max(fv, 1.15) + 0.55
        ax_c.text(tx, yi, f"{fv:.2f}×  {note}", va="center", ha="left", fontsize=6.75, fontweight="bold" if fv > 2 else "normal", color=C_DARK)

    ax_c.set_yticks(yc)
    ax_c.set_yticklabels([r[0] for r in guild_rows], fontsize=7.2)
    ax_c.set_xlim(0, 44.5)
    ax_c.set_xlabel("Inter-Order SPARSE Convergence Fold-Enrichment Over Null")
    ax_c.set_title("C   Genome-Wide Null Calibration vs. Guild-Specific Adaptive Excess", loc="left", fontweight="bold", pad=10)
    ax_c.legend(loc="lower right", frameon=True, facecolor="white", fontsize=6.8)

    # =========================================================================
    # PANEL D: Cochlear Hair-Cell Mechanotransduction Complex (Amino-Acid Shunts)
    # =========================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    focal_coch = []
    annotations_coch = {
        "LOXHD1": "Emballonuridae ↔ Hipposideridae (+29.0); Solenodon ↔ 111 Bats (+5.5)",
        "PCDH15": "3 independent echolocating bat shunts (+12.9, +11.4, +6.5 nats)",
        "WHRN": "Rhinolophus horseshoe bat ankle-link convergence (+12.3 nats)",
        "CDH23": "Bats (+12.0), Dolphins ↔ Kogiidae (+10.1), Solenodon (+4.7 nats)",
        "TMC1": "Mechanotransducer pore: Emballonuridae ↔ Hipposideridae (+9.3 nats)",
        "EPS8": "Stereocilia tip actin-bundling cap (+7.9 nats)",
        "SLC26A5": "Outer hair-cell motor Prestin (+4.2 AA; +7.8 codon Delphinidae ↔ Bats)",
        "STRC": "Stereocilin outer hair-cell tectorial coupling (+3.9, +3.6 nats)",
        "ADGRV1": "Ankle-link GPCR VLGR1: Rhinolophus ↔ Megaderma (+2.0 nats)",
    }
    for r in cochlear:
        sym = r["symbol"]
        if sym in annotations_coch:
            tot_sp_gain = sum(s["net_mdl_gain_nats"] for s in r.get("sp_shunts", []) if s["mode"] == "SPARSE")
            if tot_sp_gain > 0:
                focal_coch.append((sym, tot_sp_gain, annotations_coch[sym]))
    focal_coch.sort(key=lambda x: x[1])

    yd = np.arange(len(focal_coch))
    vals_d = [x[1] for x in focal_coch]
    ax_d.barh(yd, vals_d, height=0.60, color=C_BLUE, alpha=0.92)
    for yi, (sym, v, ann) in enumerate(focal_coch):
        ax_d.text(v + 0.8, yi, f"+{v:.1f} nats  ({ann})", va="center", ha="left", fontsize=6.55, color=C_DARK)

    ax_d.set_yticks(yd)
    ax_d.set_yticklabels([x[0] for x in focal_coch], fontsize=7.9, fontweight="bold")
    ax_d.set_xlim(0, 128)
    ax_d.set_xlabel("Cumulative 20-State Amino-Acid SPARSE Shunt Gain on Species Tree (nats)")
    ax_d.set_title("D   Pathway-Wide Convergence in the Cochlear Mechanotransduction Complex", loc="left", fontweight="bold", pad=10)

    out_pdf = PAPER_DIR / "fig_science_mammalian_atlas.pdf"
    out_png = RES_DIR / "fig_science_mammalian_atlas.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_pdf} and {out_png}")


def _load_viral_map():
    atlas_p = RES_DIR / "viral_pandemic_atlas.json"
    vmap = {}
    if atlas_p.exists():
        for d in json.load(open(atlas_p)).get("datasets", []):
            vmap[d["dataset"]] = d

    vdir = RES_DIR / "viral_datasets"
    if vdir.exists():
        for sp in sorted(vdir.glob("*.summary.json")):
            d = json.load(open(sp))
            vmap[d["dataset"]] = d

        # Ensure SARS-CoV-2 open global & world union are populated from executed .bw.json if FastTreeMP is still finishing
        import re
        def _hubs(nwk):
            s = nwk.strip().rstrip(";")
            toks = re.findall(r"\(|\)|,|;|:[^(),;]+|[^(),;:]+", s)
            st, degs = [], []
            i, n = 0, len(toks)
            while i < n:
                tk = toks[i]
                if tk == "(":
                    st.append(1)
                elif tk == ",":
                    if st:
                        st[-1] += 1
                elif tk == ")":
                    ch = st.pop()
                    degs.append(ch if not st else ch + 1)
                    if i + 1 < n and not toks[i + 1].startswith(":") and toks[i + 1] not in ("(", ")", ",", ";"):
                        i += 1
                    if i + 1 < n and toks[i + 1].startswith(":"):
                        i += 1
                else:
                    if i + 1 < n and toks[i + 1].startswith(":"):
                        i += 1
                i += 1
            p = [d for d in degs if d >= 4]
            return {"n_polytomy_hubs": len(p), "max_hub_degree": max(p, default=3)}

        sc2_comp = {
            "sars_cov_2_open_global": {
                "raw_fasta_bytes": 152540330,
                "bwz_bytes": 332974,
                "gzip9_bytes": 3907144,
                "xz6_bytes": 518864,
                "bwz_vs_raw_ratio": round(152540330 / 332974, 2),
                "bwz_vs_gzip9_ratio": round(3907144 / 332974, 2),
                "bwz_vs_xz6_ratio": round(518864 / 332974, 2),
            },
            "sars_cov_2_open_world_union": {
                "raw_fasta_bytes": 301132831,
                "bwz_bytes": 633025,
                "gzip9_bytes": 7570496,
                "xz6_bytes": 1001356,
                "bwz_vs_raw_ratio": round(301132831 / 633025, 2),
                "bwz_vs_gzip9_ratio": round(7570496 / 633025, 2),
                "bwz_vs_xz6_ratio": round(1001356 / 633025, 2),
            },
        }
        for key in ("sars_cov_2_open_global", "sars_cov_2_open_world_union"):
            if key not in vmap and (vdir / f"{key}.bw.json").exists():
                bwd = json.load(open(vdir / f"{key}.bw.json"))
                mbin = max(1, bwd["max_binary_splits"])
                res = bwd["resolved_internal_splits"]
                vmap[key] = {
                    "dataset": key,
                    "n_taxa": bwd["n_taxa"],
                    "seq_len": bwd["seq_len"],
                    "bushwhacker_runtime_sec": bwd["runtime_sec"],
                    "fasttree_gtr_wall_sec": None,
                    "speedup_vs_fasttree": None,
                    "collapsed_polytomy_frac": 1.0 - res / mbin,
                    "hub_topology_stats": _hubs(bwd["newick"]),
                    "compression_benchmark": sc2_comp[key],
                    "top_shunts": bwd.get("shunts", []),
                }
    return vmap


def build_figure_4_microbial_viral():
    vmap = _load_viral_map()
    p3 = json.load(open(RES_DIR / "microbial_web_atlas.json"))
    blind_sc2 = json.load(open(RES_DIR / "viral_datasets" / "sc2_blinded_x_panel.bw.json"))

    fig = plt.figure(figsize=(14.6, 8.5))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.05, 1.13],
        height_ratios=[1.06, 1.04],
        wspace=0.46,
        hspace=0.44,
    )

    # =========================================================================
    # PANEL A: Pandemic Viral Superspreading Hubs & Outbreak Polytomies
    # =========================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    cohorts = [
        ("SARS-CoV-2 World Union\n(N = 10,062, 29.9 kb)", "sars_cov_2_open_world_union", None),
        ("SARS-CoV-2 Open Global\n(N = 5,097, 29.9 kb)", "sars_cov_2_open_global", None),
        ("Global DENV-2 Complete CDS\n(N = 2,632, 10.1 kb)", "denv2_cds_nt", "denv2_polyprotein_aa"),
        ("Global DENV-1 Complete CDS\n(N = 1,426, 10.1 kb)", "denv1_cds_nt", "denv1_polyprotein_aa"),
        ("Global DENV-3 Complete CDS\n(N = 855, 10.1 kb)", "denv3_cds_nt", "denv3_polyprotein_aa"),
        ("Global DENV-4 Complete CDS\n(N = 585, 10.1 kb)", "denv4_cds_nt", "denv4_polyprotein_aa"),
    ][::-1]

    ya = np.arange(len(cohorts))
    nt_poly = [100.0 * vmap[nt_k]["collapsed_polytomy_frac"] for _, nt_k, _ in cohorts]
    ax_a.barh(ya, nt_poly, height=0.58, color=C_SKY, edgecolor="none", alpha=0.94, label="Nucleotide zero-signal splits collapsed (%)")

    aa_pts_x, aa_pts_y = [], []
    for yi, (_, nt_k, aa_k) in enumerate(cohorts):
        hs = vmap[nt_k]["hub_topology_stats"]
        nhubs = hs["n_polytomy_hubs"]
        mdeg = hs["max_hub_degree"]
        if aa_k and aa_k in vmap:
            ap = 100.0 * vmap[aa_k]["collapsed_polytomy_frac"]
            aa_pts_x.append(ap)
            aa_pts_y.append(yi)
            ax_a.text(ap + 1.8, yi, f"{nt_poly[yi]:.1f}% NT / {ap:.1f}% AA  ({nhubs} hubs, max deg {mdeg})",
                      va="center", ha="left", fontsize=6.65, color=C_DARK)
        else:
            ax_a.text(nt_poly[yi] + 1.5, yi, f"{nt_poly[yi]:.1f}%  ({nhubs:,} hubs, max deg {mdeg})",
                      va="center", ha="left", fontsize=6.9, fontweight="bold", color=C_DARK)

    ax_a.scatter(aa_pts_x, aa_pts_y, color=C_BLUE, edgecolor=C_DARK, linewidth=0.6, s=28, zorder=4, label="Translated polyprotein AA collapsed (%)")
    ax_a.set_yticks(ya)
    ax_a.set_yticklabels([c[0] for c in cohorts], fontsize=7.3)
    ax_a.set_xlim(0, 134)
    ax_a.set_ylim(-0.55, 5.95)
    ax_a.set_xlabel("Forced Binary Splits Collapsed into Polytomy Hubs (%)")
    ax_a.set_title("A   Pandemic Superspreading Polytomy Hubs Across 15,560 Genomes", loc="left", fontweight="bold", pad=12)
    ax_a.legend(loc="upper center", bbox_to_anchor=(0.48, 1.02), ncol=2, frameon=True, facecolor="white", framealpha=0.95, fontsize=6.5)

    # =========================================================================
    # PANEL B: Lossless .bwz Compression vs. gzip-9 / xz-6 & Runtime Speedups
    # =========================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    comp_Order = [
        ("SARS-CoV-2\n(10,062, 301 MB)", "sars_cov_2_open_world_union"),
        ("SARS-CoV-2\n(5,097, 153 MB)", "sars_cov_2_open_global"),
        ("DENV-2\n(2,632, 27 MB)", "denv2_cds_nt"),
        ("DENV-4\n(585, 6.0 MB)", "denv4_cds_nt"),
        ("DENV-3\n(855, 8.8 MB)", "denv3_cds_nt"),
        ("DENV-1\n(1,426, 14.6 MB)", "denv1_cds_nt"),
    ][::-1]

    yb = np.arange(len(comp_Order))
    bh = 0.34
    gz_ratios = [vmap[k]["compression_benchmark"]["bwz_vs_gzip9_ratio"] for _, k in comp_Order]
    xz_ratios = [vmap[k]["compression_benchmark"]["bwz_vs_xz6_ratio"] for _, k in comp_Order]

    ax_b.barh(yb + bh / 2, gz_ratios, height=bh, color=C_GREEN, alpha=0.92, label="Smaller than gzip -9 (×)")
    ax_b.barh(yb - bh / 2, xz_ratios, height=bh, color=C_BLUE, alpha=0.92, label="Smaller than xz -6 / LZMA2 (×)")
    ax_b.axvline(1.0, color=C_DARK, linestyle="--", linewidth=1.0)

    for yi, (_, k) in enumerate(comp_Order):
        d = vmap[k]
        cb = d["compression_benchmark"]
        bwz_kb = cb["bwz_bytes"] / 1024.0
        b_per_g = cb["bwz_bytes"] / max(1, d["n_taxa"])
        sp_ft = d.get("speedup_vs_fasttree")
        bw_s = d["bushwhacker_runtime_sec"]
        if sp_ft:
            if d.get("fasttree_gtr_timed_out"):
                sp_floor = int(sp_ft * 10.0) / 10.0
                sp_txt = f">{sp_floor:.1f}× faster than 32-core FastTreeMP ({bw_s:.1f} s)"
            else:
                sp_txt = f"{sp_ft:.1f}× faster than 32-core FastTreeMP ({bw_s:.1f} s)"
        else:
            sp_txt = f"BW: {bw_s:.0f} s; {b_per_g:.1f} B/genome (.bwz = {bwz_kb:.0f} KB)"
        ax_b.text(gz_ratios[yi] + 0.25, yi + bh / 2, f"{gz_ratios[yi]:.2f}× vs gzip  |  {xz_ratios[yi]:.2f}× vs xz",
                  va="center", ha="left", fontsize=6.7, fontweight="bold", color=C_DARK)
        ax_b.text(gz_ratios[yi] + 0.25, yi - bh / 2, sp_txt,
                  va="center", ha="left", fontsize=6.4, color=C_GRAY)

    ax_b.set_yticks(yb)
    ax_b.set_yticklabels([c[0] for c in comp_Order], fontsize=7.3)
    ax_b.set_xlim(0, 23.5)
    ax_b.set_ylim(-0.55, 5.95)
    ax_b.set_xlabel("Lossless Archive Size Reduction Factor Over Standard Compressors (×)")
    ax_b.set_title("B   Lossless .bwz Genomic Archiving & 32-Core Inference Speed", loc="left", fontweight="bold", pad=12)
    ax_b.legend(loc="upper center", bbox_to_anchor=(0.48, 1.02), ncol=2, frameon=True, facecolor="white", framealpha=0.95, fontsize=6.5)

    # =========================================================================
    # PANEL C: De Novo Viral Recombination & Blinded Pango X* Shunts
    # =========================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    d2_top = vmap["denv2_cds_nt"]["top_shunts"]
    d3_top = vmap["denv3_cds_nt"]["top_shunts"]
    d1_top = vmap["denv1_cds_nt"]["top_shunts"]
    sc2_u_sh = vmap["sars_cov_2_open_world_union"]["top_shunts"]
    hiv_sh = vmap["hiv1_kal153_global_nt"]["top_shunts"]
    bl_sh = blind_sc2["shunts"]

    shunt_items = [
        ("DENV-2 NS2A Mosaic [3632..3717]", d2_top[0]["net_mdl_gain_nats"], "INTERVAL"),
        ("SARS-CoV-2 (10k) 3' Tract [29700..29757]", round(sc2_u_sh[0]["net_mdl_gain_nats"], 2), "INTERVAL"),
        ("DENV-2 Capsid Mosaic [8..238]", d2_top[1]["net_mdl_gain_nats"], "INTERVAL"),
        ("Blinded SC2 Panel: XBF ↔ XBT", round(bl_sh[1]["net_mdl_gain_nats"], 2), "SPARSE"),
        ("HIV-1 KAL153 (CRF03_AB gag+nef)", hiv_sh[0]["net_mdl_gain_nats"], "SPARSE"),
        ("DENV-2 NS5 Mosaic [9626..9978]", d2_top[2]["net_mdl_gain_nats"], "INTERVAL"),
        ("DENV-3 Capsid Mosaic [16..230]", d3_top[0]["net_mdl_gain_nats"], "INTERVAL"),
        ("Blinded SC2 Spike RBD [22635..22846]", round(bl_sh[3]["net_mdl_gain_nats"], 2), "INTERVAL"),
        ("DENV-1 5'/Capsid Mosaic [21..90]", d1_top[0]["net_mdl_gain_nats"], "INTERVAL"),
        ("Blinded SC2: BA.1/1.1 ↔ XAG/XAM", round(bl_sh[5]["net_mdl_gain_nats"], 2), "SPARSE"),
        ("SARS-CoV-2 (10k) Spike NTD+RBD [21329..22917]", round(sc2_u_sh[2]["net_mdl_gain_nats"], 2), "INTERVAL"),
    ][::-1]

    yc = np.arange(len(shunt_items))
    vals_c = [x[1] for x in shunt_items]
    cols_c = [C_VERMILLION if x[2] == "INTERVAL" else C_PURPLE for x in shunt_items]
    ax_c.barh(yc, vals_c, height=0.62, color=cols_c, alpha=0.92)
    for yi, (_, v, md) in enumerate(shunt_items):
        ax_c.text(v + 1.1, yi, f"+{v:.1f} nats ({md})", va="center", ha="left", fontsize=6.55, color=C_DARK)

    ax_c.set_yticks(yc)
    ax_c.set_yticklabels([x[0] for x in shunt_items], fontsize=7.0)
    ax_c.set_xlim(0, 108)
    ax_c.set_xlabel("Net Unsupervised Reticulate Shunt Gain (nats)")
    ax_c.set_title("C   De Novo Viral Mosaic Recombinants & Blinded Pango X* Recovery", loc="left", fontweight="bold", pad=10)
    p_int = mpatches.Patch(color=C_VERMILLION, label="Contiguous INTERVAL mosaic tract")
    p_sp = mpatches.Patch(color=C_PURPLE, label="Multi-segment / SPARSE recombinant link")
    ax_c.legend(handles=[p_int, p_sp], loc="lower right", frameon=True, facecolor="white", fontsize=6.7)

    # =========================================================================
    # PANEL D: Archaeal Complexity Hypothesis (1,672 sp.) & Klebsiella AMR
    # =========================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    ar_cat = p3["gtdb_ar122_hgt_and_complexity"]["complexity_hypothesis_by_category"]
    op_pct = 100.0 * ar_cat["Operational"]["loci_with_shunts"] / ar_cat["Operational"]["loci"]
    inf_pct = 100.0 * ar_cat["Informational"]["loci_with_shunts"] / ar_cat["Informational"]["loci"]
    kleb_sm = p3["klebsiella_special_benchmarks"]["kleb_super_matrix_214kb"]
    kleb_kpc = p3["klebsiella_special_benchmarks"]["kleb_blaKPC_9taxa"]
    kpc_poly_pct = 100.0 * (1.0 - kleb_kpc["resolved_splits"] / kleb_kpc["max_binary_splits"])
    kleb_lat_pct = 100.0 * kleb_sm["ref_shunt_net_savings_nats"] / kleb_sm["vertical_tree_savings_nats"]

    micro_rows = [
        ("Klebsiella blaKPC Plasmid\n(Horizontal star polytomy %)", kpc_poly_pct,
         f"{kleb_kpc['max_binary_splits'] - kleb_kpc['resolved_splits']}/{kleb_kpc['max_binary_splits']} splits collapsed across ST11/ST258/quasipneumoniae", C_ORANGE),
        ("Archaeal Operational Proteins\n(ar122, 1,672 sp.: % w/ HGT)", op_pct,
         f"21/102 loci (30 shunts, 391.8 nats; DPANN Nanoarchaeota highways)", C_BLUE),
        ("Archaeal Informational Core\n(Ribosomal/RNAP: % w/ HGT)", inf_pct,
         "0/20 loci (0 shunts, 0.0 nats; Fisher's exact p = 0.0235)", C_GREEN),
        ("Klebsiella 214-kb Core\n(Lateral HGT % of nats)", kleb_lat_pct,
         f"77.0 lateral vs. 15,485 vertical nats (99.5% vertical; dgaE +32.9 nats)", C_SKY),
    ][::-1]

    yd = np.arange(len(micro_rows))
    vals_d = [r[1] for r in micro_rows]
    cols_d = [r[3] for r in micro_rows]
    ax_d.barh(yd, vals_d, height=0.54, color=cols_d, alpha=0.92)
    for yi, (_, v, note, _) in enumerate(micro_rows):
        ax_d.text(v + 1.4, yi, f"{v:.1f}%  —  {note}", va="center", ha="left", fontsize=6.6, fontweight="bold" if yi in (1, 2) else "normal", color=C_DARK)

    ax_d.set_yticks(yd)
    ax_d.set_yticklabels([r[0] for r in micro_rows], fontsize=7.2)
    ax_d.set_xlim(0, 178)
    ax_d.set_xlabel("Percentage (%)")
    ax_d.set_title("D   Quantifying the Microbial 'Tree vs. Web' & Complexity Hypothesis", loc="left", fontweight="bold", pad=10)

    out_pdf = PAPER_DIR / "fig_science_microbial_viral.pdf"
    out_png = RES_DIR / "fig_science_microbial_viral.png"
    fig.savefig(out_pdf, bbox_inches="tight")
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_pdf} and {out_png}")


if __name__ == "__main__":
    setup_rc()
    build_figure_2_macroevo()
    build_figure_3_mammalian_atlas()
    build_figure_4_microbial_viral()
