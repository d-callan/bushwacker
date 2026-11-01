# BUSHWHACKER

**B**ipartition **U**nification, **S**hunt **H**unting, & **W**eighted **H**ub **A**ggregation via **C**ombinadic **K**olmogorov **E**ntropy **R**eduction

[![Language: Rust](https://img.shields.io/badge/Language-Rust%20%28Edition%202021%29-DEA584.svg)](#installation--rust-cli)
[![Technical Report](https://img.shields.io/badge/Technical%20Report-PDF%20%2817%20pp.%29-blue.svg)](paper/bushwhacker_report.pdf)

> *"Clear the statistical underbrush."*

**BUSHWHACKER** is a first-principles **Minimum Description Length (MDL) / algorithmic Kolmogorov complexity** engine for phylogenetic graph induction, implemented natively in zero-dependency **Rust** ([`src/lib.rs`](src/lib.rs), [`src/main.rs`](src/main.rs)) with a reference **Python/NumPy** implementation ([`python/mdl_phylo.py`](python/mdl_phylo.py), [`python/bushwhacker.py`](python/bushwhacker.py)).

Instead of forcing every multiple sequence alignment onto a strictly bifurcating tree and evaluating it by in-sample point-estimate maximum likelihood `ln L(X | T, θ)`, **BUSHWHACKER** treats an alignment `X` (`N` taxa × `L` sites) as a discrete data table to be transmitted losslessly over a noiseless two-part Rissanen–Wallace sender–receiver channel:

$$
G^* = \arg\min_{G = (V, E_{\text{tree}}, E_{\text{shunt}})} \left( L_{\text{model}}(G) + L_{\text{data}}(X \mid G) \right)
$$

Every proposed evolutionary structure—an internal tree split `e = (u, v)`, a multifurcating polytomy hub, a recombinant mosaic tract (`INTERVAL` shunt), or a sparse convergent replacement event (`SPARSE` shunt)—must pay for its own exact combinatorial address cost in **nats** before it is admitted into the induced evolutionary graph.

---

## Key Capabilities

1. **Zero-Parameter, Bootstrap-Free Polytomy Calibration (`0.0% FPR`):**
   Every internal split dividing a `(d_u + d_v)`-way node into subtrees of degrees `d_u ≥ 2` and `d_v ≥ 2` pays an exact combinatorial topology penalty `L_topo(d_u, d_v) = ln((d_u + d_v)/2 - 1) + ln C(d_u + d_v, d_u)` plus a combinadic branch header `L*(m_e) + ln C(L_eff, m_e)`. Unsupported branches contract deterministically (`MDL_Gain ≤ 0 nats`) into explicit multifurcations without needing 1,000 bootstrap replicates or arbitrary support cutoffs.
2. **Intrinsic Immunity to Long-Branch Attraction (`96.0%` True Split in the Felsenstein Zone):**
   In the Felsenstein Zone (`b_long = 0.40`, `b_short = 0.03`, `b_int = 0.035`), a random mutation on a long branch with `m_c` mutations carries marginal combinadic positional cost `ln((L_eff - m_c) / (m_c + 1)) ≈ 0.6 nats`, whereas a mutation on a short branch carries `≈ 3.8 nats`. Uniting the two short branches yields `+5.8 nats` of compression per synapomorphy, while uniting the two long branches fails to overcome the quartet topology penalty `L_topo(2, 2) = ln 3 nats`.
3. **Dual-Mode Reticulate Shunt Peeling (`INTERVAL` vs. `SPARSE`):**
   When localized recombination or multi-site adaptive convergence occurs across distant lineages, BUSHWHACKER preserves the vertical species backbone `T_0` and peels off the horizontal/convergent signal as a certified reticulate shunt `e_s = (u, v, M_s)`, automatically classifying contiguous mosaic blocks (`INTERVAL` mode, `L_mask = ln 2 + 2 ln L`) versus dispersed functional convergence (`SPARSE` mode, `L_mask(m) = ln 2 + ln C(M_pool, m)`).

---

## Installation & Rust CLI

BUSHWHACKER is written in pure, zero-dependency Rust (`Edition 2021`) for instant compilation and native speed (`0.023s` on 17-taxon β-globin; `0.358s` on 116-taxon mammalian *Prestin*).

```bash
cargo build --release
```

### Usage

```text
./target/release/bushwhacker <input.fasta> [options]

Options:
  -o, --out-newick <PATH>      Output Newick file path (with MDL branch supports & lengths)
  -j, --out-json <PATH>        Output JSON report path (full MDL ledger, splits, & shunts)
  --no-coding                  Disable automatic coding-frame 3rd-position entropy weighting
  --no-shunts                  Disable Stage 4 reticulate shunt peeling
  --min-shunt-sites <N>        Minimum parallel mutations to evaluate a shunt (default: 2)
  -h, --help                   Print help message
```

### Example: Mammalian *Prestin* (`prestin.fasta`, 116 Mammals, 2,235 bp)

```bash
./target/release/bushwhacker prestin.fasta -o prestin.nwk -j prestin.json
```

```text
==============================================================================
BUSHWHACKER Summary: prestin.fasta  (N=116 taxa, L=2235 sites)
------------------------------------------------------------------------------
  Raw Uncompressed DL   :     351367.4 nats
  Star Tree MDL (Null)  :      80036.0 nats
  Induced Tree MDL      :      43812.9 nats  (54.74% of Star)
  Induced Network MDL   :      43808.0 nats  (54.74% of Star)
  Resolved Splits       :   95 /  113 binary splits (12 polytomy hubs)
  Certified Shunts      :    1
    [1] SPARSE  N140 <-> N167  (+4.94 nats, P=99.3%, m=4 sites, codons=[7, 384, 392, 601])
  Runtime               : 0.3577 s
==============================================================================
```

Without any species-tree supervision or empirical amino acid matrices, **BUSHWHACKER**:
- Resolves the 116-mammal species backbone at **88.0% precision** against the curated NCBI Taxonomy (`RF = 37`, outperforming raw and support-collapsed IQ-TREE 2, FastTree 2, BioNJ, and Maximum Parsimony), and
- Automatically discovers and certifies the classic **echolocation convergence shunt** (`SPARSE N140 Dolphins <-> N167 Echolocating Bats`, `+4.94 nats`, `P = 99.3%`) on the **exact four classic echolocation codons `7, 384, 392, 601` (`N7T, I384T, S392A, G601A`)** of Liu et al. (*PNAS* 2010) and Li et al. (*Curr. Biol.* 2010).

---

## Summary of Benchmark & Stress-Test Results

Full derivations, proofs, and multi-panel figures are available in [`paper/bushwhacker_report.pdf`](paper/bushwhacker_report.pdf).

### 1. Polytomy Calibration & Support-Threshold Comparison (`N = 20` Taxa, 16 Zero-Length Polytomy Branches)

| Method / Support Cutoff | Mean False Splits (FP / 16) | False Positive Rate (FPR) | Precision (%) | Recall (%) | RF Error |
| :--- | :---: | :---: | :---: | :---: | :---: |
| **BUSHWHACKER (Zero Bootstrap)** | **0.00** | **0.0%** | **100.0%** | **96.9%** | **0.40** |
| IQ-TREE 2 (Raw ML) | 16.00 | 100.0% | 19.6% | 97.5% | 16.10 |
| IQ-TREE 2 (`SH-aLRT ≥ 80%`) | 0.90 | 5.6% | 83.2% | 95.0% | 1.10 |
| IQ-TREE 2 (`aBayes ≥ 0.95`) | 0.70 | 4.4% | 86.2% | 95.0% | 0.90 |
| IQ-TREE 2 (`UFBoot ≥ 95%`) | 0.30 | 1.9% | 94.0% | 87.5% | 0.80 |
| FastTree 2 (Raw) | 16.00 | 100.0% | 19.6% | 97.5% | 16.10 |
| FastTree 2 (`SH ≥ 0.90`) | 0.50 | 3.1% | 88.9% | 95.0% | 0.70 |
| BioNJ (`JC69`) | 16.00 | 100.0% | 20.0% | 100.0% | 16.00 |
| Maximum Parsimony (`NNI`) | 16.00 | 100.0% | 20.0% | 100.0% | 16.00 |

### 2. Classical Felsenstein Zone (`L = 1,000 bp`, 50 Replicates)

| Method | True Split `(A,B) \| (C,D)` (%) | LBA Trap `(A,C) \| (B,D)` (%) | Polytomy (%) |
| :--- | :---: | :---: | :---: |
| **BUSHWHACKER (Ours)** | **96.0%** | **0.0%** | **2.0%** |
| IQ-TREE 2 (`GTR+G`) | 94.0% | 4.0% | 0.0% |
| FastTree 2 (`GTR+CAT`) | 90.0% | 4.0% | 0.0% |

### 3. Multi-Clade Parasite LBA (32 Taxa, 4 Long-Branch Parasites Nested in 4 Clades)

| Method | True Parasite–Host Sister Recall (%) | False Parasite–Parasite LBA Trap (%) |
| :--- | :---: | :---: |
| **BUSHWHACKER (Ours)** | **91.7%** | **0.0%** |
| IQ-TREE 2 (`GTR+G`) | 79.2% | 8.3% |
| FastTree 2 (`GTR+CAT`) | 45.8% | 4.2% |
| BioNJ (`JC69`) | 29.2% | 58.3% |
| Maximum Parsimony (`NNI`) | 8.3% | 75.0% |

### 4. Empirical Mammalian *Prestin* (116 Mammals vs. Curated NCBI Species Tree)

| Method / Support Rule | Resolved Splits | True Positives (TP) | False Positives (FP) | Precision (%) | Recall (%) | F1 Score (%) | RF Distance |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **BUSHWHACKER (+ Echoloc. Shunt)** | **83** | **73** | **10** | **88.0%** | **73.0%** | **79.8%** | **37** |
| IQ-TREE 2 (`GTR+G` Raw) | 100 | 76 | 24 | 76.0% | 76.0% | 76.0% | 48 |
| IQ-TREE 2 (`SH-aLRT ≥ 80%`) | 63 | 55 | 8 | 87.3% | 55.0% | 67.5% | 53 |
| IQ-TREE 2 (`aBayes ≥ 0.95`) | 60 | 52 | 8 | 86.7% | 52.0% | 65.0% | 56 |
| FastTree 2 (`GTR+CAT` Raw) | 100 | 74 | 26 | 74.0% | 74.0% | 74.0% | 52 |
| FastTree 2 (`SH ≥ 0.90`) | 66 | 56 | 10 | 84.8% | 56.0% | 67.5% | 54 |
| BioNJ (`JC69`) | 100 | 77 | 23 | 77.0% | 77.0% | 77.0% | 46 |
| Maximum Parsimony (`NNI`) | 100 | 70 | 30 | 70.0% | 70.0% | 70.0% | 60 |

---

## Scaling Performance

Asymptotically, **BUSHWHACKER** scales as `O(N² L + N³ + R_NNI · N · L log L)`, which is **linear in alignment length `L`** (up to a `log L` sort) and **sub-quadratic to quadratic in taxa `N`** for typical datasets:

| Taxa (`N`) | Sites (`L`) | Rust Engine Runtime (s) | Scaling Mode |
| :---: | :---: | :---: | :--- |
| **16** | 1,000 | **0.0069 s** | Taxa scaling (`L = 1,000 bp`) |
| **32** | 1,000 | **0.0216 s** | Taxa scaling (`L = 1,000 bp`) |
| **64** | 1,000 | **0.0677 s** | Taxa scaling (`L = 1,000 bp`) |
| **128** | 1,000 | **0.2285 s** | Taxa scaling (`L = 1,000 bp`) |
| **256** | 1,000 | **0.7511 s** | Taxa scaling (`L = 1,000 bp`) |
| **512** | 1,000 | **3.1132 s** | Taxa scaling (`L = 1,000 bp`) |
| 64 | **500** | **0.0308 s** | Site scaling (`N = 64 taxa`) |
| 64 | **2,000** | **0.1692 s** | Site scaling (`N = 64 taxa`) |
| 64 | **5,000** | **0.4821 s** | Site scaling (`N = 64 taxa`) |
| 64 | **10,000** | **0.9813 s** | Site scaling (`N = 64 taxa`) |

---

## Repository Structure

```text
bushwacker/
├── Cargo.toml                             # Zero-dependency Rust crate manifest
├── README.md                              # Project documentation
├── src/
│   ├── lib.rs                             # Native Rust BUSHWHACKER core engine
│   └── main.rs                            # Native Rust CLI binary
├── python/
│   ├── bushwhacker.py                     # Reference Python CLI
│   └── mdl_phylo.py                       # Reference Python/NumPy MDL engine
├── benchmarks/
│   ├── benchmark_first_principles.py      # 5-suite synthetic & empirical benchmark harness
│   ├── benchmark_results.json             # Executed 5-suite benchmark outputs
│   ├── stress_test_vs_paradigms.py        # 5-axis adversarial stress-test harness vs. ML/Bayesian/NJ/MP
│   ├── stress_test_results.json           # Executed stress-test outputs
│   └── generate_report_figures.py         # Publication vector PDF figure generator
└── paper/
    ├── bushwhacker_report.pdf             # Compiled 17-page technical report
    ├── bushwhacker_report.tex             # LaTeX source manuscript
    ├── fig1_synthetic_benchmarks.pdf      # Figure 1 (4-regime synthetic benchmarks)
    ├── fig2_empirical_compression.pdf     # Figure 2 (Empirical compression & polytomy ratios)
    └── fig3_stress_tests.pdf              # Figure 3 (5-axis ruthless stress tests)
```
