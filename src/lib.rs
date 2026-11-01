//! BUSHWHACKER: Bipartition Unification, Shunt Hunting, & Weighted Hub Aggregation
//! via Combinadic Kolmogorov Entropy Reduction.
//!
//! Native Rust implementation of the first-principles Minimum Description Length (MDL)
//! phylogenetic graph induction engine.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

// ---------------------------------------------------------------------------
// Combinatorial Information-Theoretic Primitives (All in Nats)
// ---------------------------------------------------------------------------

/// High-precision Lanczos approximation of ln(Gamma(z)) for z > 0.
#[inline]
pub fn ln_gamma(z: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    if z < 0.5 {
        let pi = std::f64::consts::PI;
        return pi.ln() - (pi * z).sin().ln() - ln_gamma(1.0 - z);
    }
    let z_m1 = z - 1.0;
    let mut x = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        x += c / (z_m1 + i as f64);
    }
    let t = z_m1 + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z_m1 + 0.5) * t.ln() - t + x.ln()
}

/// Continuous exact log-binomial coefficient ln C(n, k) via ln_gamma.
#[inline]
pub fn log_binom(n: f64, k: f64) -> f64 {
    if k <= 0.0 || k >= n || n <= 0.0 {
        return 0.0;
    }
    ln_gamma(n + 1.0) - ln_gamma(k + 1.0) - ln_gamma(n - k + 1.0)
}

/// Rissanen universal prefix code length (in nats) for non-negative integer m.
#[inline]
pub fn universal_int_nats(m: f64) -> f64 {
    if m <= 0.0 {
        return 2.0_f64.ln();
    }
    2.0_f64.ln() + (1.0 + m).ln() + 0.5 * (1.0 + (1.0 + m).ln()).ln()
}

/// Exact enumerative code length (in nats) to transmit m_e mutations along a
/// single branch e over L_eff effective sites.
#[inline]
pub fn branch_mutation_code_nats(
    m_e: f64,
    l_eff: f64,
    subst_nats: f64,
    alphabet_size: usize,
) -> f64 {
    if m_e <= 1e-9 {
        return 0.0;
    }
    let m_pos = m_e.min(0.49 * l_eff);
    let header = universal_int_nats(m_e);
    let mut site_positions = log_binom(l_eff, m_pos);
    if m_e > 0.49 * l_eff {
        site_positions += (m_e - 0.49 * l_eff) * 2.0_f64.ln();
    }
    let default_state_cost = m_e * ((alphabet_size.saturating_sub(1).max(1)) as f64).ln();
    let state_cost = if subst_nats > 0.0 {
        subst_nats
    } else {
        default_state_cost
    };
    header + site_positions + state_cost
}

/// Exact combinatorial cost (in nats) of resolving a k = (deg_a + deg_b)-way
/// polytomy into two sub-nodes of outward degrees deg_a >= 2 and deg_b >= 2.
#[inline]
pub fn split_topology_cost_nats(deg_a: usize, deg_b: usize) -> f64 {
    let k = deg_a + deg_b;
    if k < 4 || deg_a < 2 || deg_b < 2 {
        return 0.0;
    }
    let num_sizes = ((k / 2).saturating_sub(1)).max(1) as f64;
    let mut subset_cost = log_binom(k as f64, deg_a as f64);
    if deg_a == deg_b {
        subset_cost -= 2.0_f64.ln();
    }
    num_sizes.ln() + subset_cost
}

// ---------------------------------------------------------------------------
// Data Structures for the MDL Phylogenetic Compression Graph
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
pub struct ShuntEdge {
    pub node_u: usize,
    pub node_v: usize,
    pub label_u: String,
    pub label_v: String,
    pub mode: String, // "INTERVAL" or "SPARSE"
    pub sites: Vec<usize>,
    pub interval_span: Option<(usize, usize)>,
    pub shared_mutations: usize,
    pub data_savings_nats: f64,
    pub mask_cost_nats: f64,
    pub endpoint_cost_nats: f64,
    pub net_mdl_gain_nats: f64,
    pub posterior_prob: f64,
}

#[derive(Clone, Debug)]
pub struct BranchRecord {
    pub u: usize,
    pub v: usize,
    pub is_internal: bool,
    pub mutation_count: f64,
    pub branch_length: f64,
    pub mdl_gain_nats: f64,
    pub mdl_support: f64,
    pub synapomorphy_sites: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct BushwhackerResult {
    pub names: Vec<String>,
    pub n_taxa: usize,
    pub n_sites: usize,
    pub l_eff: f64,
    pub newick: String,
    pub adjacency: BTreeMap<usize, Vec<usize>>,
    pub node_labels: BTreeMap<usize, String>,
    pub branches: Vec<BranchRecord>,
    pub shunts: Vec<ShuntEdge>,
    pub n_internal_branches: usize,
    pub max_binary_internal_branches: usize,
    pub polytomy_nodes: usize,
    pub max_polytomy_degree: usize,
    pub raw_uncompressed_nats: f64,
    pub star_tree_nats: f64,
    pub backbone_tree_nats: f64,
    pub full_graph_nats: f64,
    pub runtime_sec: f64,
}

// ---------------------------------------------------------------------------
// Genetic Code & Alignment Encoding Utilities
// ---------------------------------------------------------------------------

pub const IS_PURINE: [u8; 4] = [1, 0, 1, 0];

/// Translate a 6-bit encoded codon `(b0 << 4) | (b1 << 2) | b2` (where A=0, C=1, G=2, T=3)
/// into its single-letter amino acid (or b'*' for stop).
#[inline]
pub fn translate_codon_idx(idx: u8) -> u8 {
    // 64 codons in order AAA(0), AAC(1), AAG(2), AAT(3), ACA(4), ...
    const TABLE: &[u8; 64] = b"KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
    if (idx as usize) < 64 {
        TABLE[idx as usize]
    } else {
        b'X'
    }
}

#[inline]
pub fn translate_bases(b0: i8, b1: i8, b2: i8) -> u8 {
    if b0 < 0 || b1 < 0 || b2 < 0 {
        return b'X';
    }
    let idx = ((b0 as u8) << 4) | ((b1 as u8) << 2) | (b2 as u8);
    translate_codon_idx(idx)
}

pub fn parse_fasta_file<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, Vec<String>), String> {
    let file = File::open(path.as_ref()).map_err(|e| format!("Failed to open FASTA: {e}"))?;
    let reader = BufReader::new(file);
    let mut names = Vec::new();
    let mut seqs = Vec::new();
    let mut curr_name: Option<String> = None;
    let mut curr_seq = String::new();

    for line_res in reader.lines() {
        let line = line_res.map_err(|e| format!("Read error: {e}"))?;
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(name) = curr_name.take() {
                names.push(name);
                seqs.push(std::mem::take(&mut curr_seq));
            }
            let first_tok = rest.split_whitespace().next().unwrap_or("").to_string();
            curr_name = Some(first_tok);
        } else {
            curr_seq.push_str(trimmed);
        }
    }
    if let Some(name) = curr_name {
        names.push(name);
        seqs.push(curr_seq);
    }
    Ok((names, seqs))
}

/// Encode N aligned nucleotide sequences of length L into flat row-major arrays:
/// - `mat`: `(N * L)` `i8` in `{0, 1, 2, 3, -1}`
/// - `bitmasks`: `(N * L)` `u8` in `{1, 2, 4, 8, 15}`
/// - `majority_state`: `(L)` `i8`
pub fn encode_nucleotide_alignment(seqs: &[String]) -> (Vec<i8>, Vec<u8>, Vec<i8>, usize, usize) {
    let n = seqs.len();
    let l = seqs.first().map(|s| s.len()).unwrap_or(0);
    let mut mat = vec![-1_i8; n * l];
    let mut bitmasks = vec![15_u8; n * l];
    let mut counts = vec![[0_u32; 4]; l];

    for (i, s) in seqs.iter().enumerate() {
        let row_off = i * l;
        for (j, &byte) in s.as_bytes().iter().take(l).enumerate() {
            let v: i8 = match byte {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' | b'U' | b'u' => 3,
                _ => -1,
            };
            if v >= 0 {
                mat[row_off + j] = v;
                bitmasks[row_off + j] = 1_u8 << (v as u8);
                counts[j][v as usize] += 1;
            }
        }
    }

    let mut majority_state = vec![0_i8; l];
    for j in 0..l {
        let mut best_b = 0_i8;
        let mut best_c = counts[j][0];
        for b in 1..4 {
            if counts[j][b] > best_c {
                best_c = counts[j][b];
                best_b = b as i8;
            }
        }
        majority_state[j] = best_b;
    }

    (mat, bitmasks, majority_state, n, l)
}

pub fn detect_coding_frame_and_weights(
    mat: &[i8],
    n: usize,
    l: usize,
) -> (Vec<f64>, Vec<f64>, bool) {
    let mut dist_weights = vec![1.0_f64; l];
    let mut tree_weights = vec![1.0_f64; l];
    let mut is_coding = false;

    if l >= 90 && l % 3 == 0 {
        let n_codons = l / 3;
        let mut stop_count = 0usize;
        let mut valid_codons = 0usize;
        let check_n = n.min(32);
        for i in 0..check_n {
            let row = &mat[i * l..(i + 1) * l];
            for c in 0..(n_codons.saturating_sub(1)) {
                let b0 = row[3 * c];
                let b1 = row[3 * c + 1];
                let b2 = row[3 * c + 2];
                if b0 >= 0 && b1 >= 0 && b2 >= 0 {
                    valid_codons += 1;
                    if translate_bases(b0, b1, b2) == b'*' {
                        stop_count += 1;
                    }
                }
            }
        }
        if valid_codons > 0 && (stop_count as f64 / valid_codons as f64) < 0.02 {
            is_coding = true;
            for c in 0..n_codons {
                let c0 = 3 * c;
                let c1 = 3 * c + 1;
                let c2 = 3 * c + 2;
                let mut cod_mask = 0_u64;
                let mut aa_mask = 0_u32;
                for i in 0..n {
                    let row = &mat[i * l..(i + 1) * l];
                    let b0 = row[c0];
                    let b1 = row[c1];
                    let b2 = row[c2];
                    if b0 >= 0 && b1 >= 0 && b2 >= 0 {
                        let idx = ((b0 as u8) << 4) | ((b1 as u8) << 2) | (b2 as u8);
                        cod_mask |= 1_u64 << idx;
                        let aa = translate_codon_idx(idx);
                        let aa_bit = if aa == b'*' {
                            26
                        } else {
                            (aa.wrapping_sub(b'A')).min(25)
                        };
                        aa_mask |= 1_u32 << aa_bit;
                    }
                }
                if cod_mask.count_ones() > 1 {
                    if aa_mask.count_ones() == 1 {
                        for col in [c0, c1, c2] {
                            dist_weights[col] = 2.0;
                            tree_weights[col] = 1.25;
                        }
                    } else {
                        dist_weights[c0] = 0.12;
                        dist_weights[c1] = 0.08;
                        dist_weights[c2] = 1.00;
                        tree_weights[c0] = 0.45;
                        tree_weights[c1] = 0.40;
                        tree_weights[c2] = 1.15;
                    }
                }
            }
        }
    }
    (dist_weights, tree_weights, is_coding)
}

#[inline]
fn median_f64(vals: &mut [f64]) -> f64 {
    if vals.is_empty() {
        return 0.0;
    }
    vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mid = vals.len() / 2;
    if vals.len() % 2 == 1 {
        vals[mid]
    } else {
        0.5 * (vals[mid - 1] + vals[mid])
    }
}

// ---------------------------------------------------------------------------
// Core Engine: BushwhackerEngine
// ---------------------------------------------------------------------------

pub struct BushwhackerEngine {
    pub coding_aware: bool,
    pub detect_shunts: bool,
    pub min_shunt_sites: usize,
}

impl Default for BushwhackerEngine {
    fn default() -> Self {
        Self {
            coding_aware: true,
            detect_shunts: true,
            min_shunt_sites: 2,
        }
    }
}

impl BushwhackerEngine {
    pub fn new(coding_aware: bool, detect_shunts: bool, min_shunt_sites: usize) -> Self {
        Self {
            coding_aware,
            detect_shunts,
            min_shunt_sites,
        }
    }

    pub fn infer_from_fasta<P: AsRef<Path>>(&self, path: P) -> Result<BushwhackerResult, String> {
        let (names, seqs) = parse_fasta_file(path)?;
        self.infer(&names, &seqs)
    }

    pub fn infer(&self, names: &[String], seqs: &[String]) -> Result<BushwhackerResult, String> {
        let t0 = Instant::now();
        let n = names.len();
        if n == 0 {
            return Err("Empty alignment provided.".to_string());
        }
        let (mat, bitmasks, majority_state, _, l) = encode_nucleotide_alignment(seqs);

        let total_valid: usize = mat.iter().filter(|&&x| x >= 0).count();
        let l_eff = ((total_valid as f64) / (n as f64)).clamp(1.0, l as f64);

        let (ti_cost, tv_cost) = self.estimate_subst_surprise(&mat, n, l);
        let (mut dist_weights, mut tree_weights, is_coding) = if self.coding_aware {
            detect_coding_frame_and_weights(&mat, n, l)
        } else {
            (vec![1.0; l], vec![1.0; l], false)
        };

        let (dw_mod, tw_mod) =
            self.compute_dual_site_coherence_weights(&mat, n, l, &dist_weights, &tree_weights);
        dist_weights = dw_mod;
        tree_weights = tw_mod;

        // Step 1 & 2: Sub-quadratic O(N log N * L) 3-point terminal loads & candidate backbone
        let (terminal_loads, mut adj, node_labels) = self.build_subquadratic_loads_and_backbone(
            names,
            &mat,
            &bitmasks,
            &dist_weights,
            &tree_weights,
            n,
            l,
            l_eff,
            ti_cost,
            tv_cost,
        );

        // Step 3: MDL Edge Contraction + Polytomy Hub Clique Expansion
        let (adj_c, mut branches, _, edge_mutations) = self.contract_unsupported_edges_mdl(
            adj,
            &bitmasks,
            &majority_state,
            &terminal_loads,
            n,
            l,
            l_eff,
            ti_cost,
            tv_cost,
            &tree_weights,
            &mut node_labels.clone(),
            true,
        );
        adj = adj_c;
        // Refresh node labels for any hub-expanded nodes
        let mut full_labels = node_labels;
        for &u in adj.keys() {
            full_labels.entry(u).or_insert_with(|| format!("N{u}"));
        }

        // Step 4: Dual-Mode Residual Shunt Peeling + Coordinate-Descent Re-Verification
        let mut shunts = Vec::new();
        if self.detect_shunts && n >= 4 {
            shunts = self.peel_residual_shunts_mdl(
                &adj,
                &full_labels,
                &bitmasks,
                &majority_state,
                &terminal_loads,
                &edge_mutations,
                n,
                l,
                l_eff,
                ti_cost,
                tv_cost,
                &mat,
                is_coding,
            );
            let int_shunts: Vec<&ShuntEdge> = shunts
                .iter()
                .filter(|s| {
                    s.mode == "INTERVAL" && s.interval_span.is_some() && s.net_mdl_gain_nats >= 15.0
                })
                .collect();
            if !int_shunts.is_empty() {
                let mut masked_bitmasks = bitmasks.clone();
                let mut masked_any = false;
                for sh in int_shunts {
                    if let Some((a_sp, b_sp)) = sh.interval_span {
                        let rec_node = if sh.node_u < n {
                            Some(sh.node_u)
                        } else if sh.node_v < n {
                            Some(sh.node_v)
                        } else {
                            None
                        };
                        if let Some(r_idx) = rec_node {
                            let row = &mut masked_bitmasks[r_idx * l..(r_idx + 1) * l];
                            for s_i in a_sp..=b_sp.min(l - 1) {
                                row[s_i] = 15;
                            }
                            masked_any = true;
                        }
                    }
                }
                if masked_any {
                    let (adj_re, branches_re, _, _) = self.contract_unsupported_edges_mdl(
                        adj,
                        &masked_bitmasks,
                        &majority_state,
                        &terminal_loads,
                        n,
                        l,
                        l_eff,
                        ti_cost,
                        tv_cost,
                        &tree_weights,
                        &mut full_labels,
                        false,
                    );
                    adj = adj_re;
                    branches = branches_re;
                }
            }
        }

        // Step 5: MDL bit/nat budgets & Newick serialization
        let raw_nats = (n as f64) * l_eff * 4.0_f64.ln();
        let star_nats = self.compute_star_tree_mdl_nats(&bitmasks, n, l, l_eff);
        let backbone_nats = self.compute_tree_total_mdl_nats(&adj, &branches, l_eff);
        let shunt_gain_total: f64 = shunts.iter().map(|s| s.net_mdl_gain_nats).sum();
        let full_graph_nats = (l_eff * 4.0_f64.ln()).max(backbone_nats - shunt_gain_total);

        let mut polytomy_degrees = Vec::new();
        for (&u, nbs) in &adj {
            if u >= n && nbs.len() > 3 {
                polytomy_degrees.push(nbs.len());
            }
        }
        let n_internal_branches = branches.iter().filter(|b| b.is_internal).count();
        let max_binary_internal_branches = n.saturating_sub(3);
        let max_polytomy_degree = polytomy_degrees
            .iter()
            .copied()
            .max()
            .unwrap_or(if n >= 3 { 3 } else { n });

        let newick = self.serialize_newick(&adj, names, &branches, n);
        let runtime_sec = t0.elapsed().as_secs_f64();

        Ok(BushwhackerResult {
            names: names.to_vec(),
            n_taxa: n,
            n_sites: l,
            l_eff,
            newick,
            adjacency: adj,
            node_labels: full_labels,
            branches,
            shunts,
            n_internal_branches,
            max_binary_internal_branches,
            polytomy_nodes: polytomy_degrees.len(),
            max_polytomy_degree,
            raw_uncompressed_nats: raw_nats,
            star_tree_nats: star_nats,
            backbone_tree_nats: backbone_nats,
            full_graph_nats,
            runtime_sec,
        })
    }

    fn estimate_subst_surprise(&self, mat: &[i8], n: usize, l: usize) -> (f64, f64) {
        let k_sub = n.min(64);
        let indices: Vec<usize> = if n <= 64 {
            (0..n).collect()
        } else {
            (0..k_sub)
                .map(|i| (i * (n - 1)) / (k_sub - 1))
                .collect()
        };
        let mut ti = 1.0_f64;
        let mut tv = 2.0_f64;
        for i in 0..indices.len() {
            let r_i = &mat[indices[i] * l..(indices[i] + 1) * l];
            for j in (i + 1)..indices.len() {
                let r_j = &mat[indices[j] * l..(indices[j] + 1) * l];
                for s in 0..l {
                    let bi = r_i[s];
                    let bj = r_j[s];
                    if bi >= 0 && bj >= 0 && bi != bj {
                        if IS_PURINE[bi as usize] == IS_PURINE[bj as usize] {
                            ti += 1.0;
                        } else {
                            tv += 1.0;
                        }
                    }
                }
            }
        }
        let tot = ti + tv;
        let p_ti = (ti / tot).clamp(0.20, 0.85);
        let p_tv_each = ((1.0 - p_ti) / 2.0).max(0.05);
        (-p_ti.ln(), -p_tv_each.ln())
    }

    fn compute_dual_site_coherence_weights(
        &self,
        mat: &[i8],
        n: usize,
        l: usize,
        dist_weights: &[f64],
        tree_weights: &[f64],
    ) -> (Vec<f64>, Vec<f64>) {
        if n < 6 || l < 10 {
            return (dist_weights.to_vec(), tree_weights.to_vec());
        }

        let mut n_states = vec![0_u8; l];
        let mut w_mult = vec![1.0_f64; l];
        for s in 0..l {
            let mut mask = 0_u8;
            for i in 0..n {
                let v = mat[i * l + s];
                if v >= 0 {
                    mask |= 1 << (v as u8);
                }
            }
            let st = mask.count_ones() as u8;
            n_states[s] = st;
            if st == 3 {
                w_mult[s] = 0.90;
            } else if st >= 4 {
                w_mult[s] = 0.75;
            }
        }

        // Extract parsimony-informative column-state vectors packed into 128-bit [u64; 2] words
        let k_coh = n.min(128);
        let coh_taxa: Vec<usize> = if n <= 128 {
            (0..n).collect()
        } else {
            (0..k_coh).map(|i| (i * (n - 1)) / (k_coh - 1)).collect()
        };

        let mut char_vecs: Vec<[u64; 2]> = Vec::new();
        let mut char_bar: Vec<[u64; 2]> = Vec::new();
        let mut char_sites: Vec<usize> = Vec::new();
        let mut char_counts: Vec<usize> = Vec::new();

        for s in 0..l {
            if n_states[s] < 2 {
                continue;
            }
            let mut n_val = 0usize;
            for &ti in &coh_taxa {
                if mat[ti * l + s] >= 0 {
                    n_val += 1;
                }
            }
            if n_val < 4 {
                continue;
            }
            for b in 0..4_i8 {
                let mut cnt = 0usize;
                for &ti in &coh_taxa {
                    if mat[ti * l + s] == b {
                        cnt += 1;
                    }
                }
                if cnt >= 2 && cnt <= n_val / 2 {
                    let mut c_bits = [0_u64; 2];
                    let mut cb_bits = [0_u64; 2];
                    for (idx_t, &ti) in coh_taxa.iter().enumerate() {
                        let v = mat[ti * l + s];
                        let word = idx_t >> 6;
                        let bit = 1_u64 << (idx_t & 63);
                        if v == b {
                            c_bits[word] |= bit;
                        } else if v >= 0 {
                            cb_bits[word] |= bit;
                        }
                    }
                    char_vecs.push(c_bits);
                    char_bar.push(cb_bits);
                    char_sites.push(s);
                    char_counts.push(cnt);
                }
            }
        }

        let m_chars = char_vecs.len();
        if m_chars < 8 {
            let md = dist_weights
                .iter()
                .zip(&w_mult)
                .map(|(d, w)| d * w)
                .collect();
            let mt = tree_weights
                .iter()
                .zip(&w_mult)
                .map(|(t, w)| t * w)
                .collect();
            return (md, mt);
        }

        let mut ref_idx: Vec<usize> = (0..m_chars).collect();
        if m_chars > 350 {
            ref_idx.sort_by_key(|&idx| {
                let c = char_counts[idx];
                let bal = c.min(k_coh.saturating_sub(c));
                (std::cmp::Reverse(bal), idx)
            });
            ref_idx.truncate(350);
        }

        let delta_local = (l / 6).clamp(40, 180);
        let mut coh_scores = vec![1.0_f64; m_chars];

        for k in 0..m_chars {
            let ck = char_vecs[k];
            let cbk = char_bar[k];
            let sk = char_sites[k];
            let mut sum_overlap = 0.0_f64;
            let mut cnt_overlap = 0usize;
            let mut sum_long = 0.0_f64;
            let mut cnt_long = 0usize;

            for &r in &ref_idx {
                let sr = char_sites[r];
                if sk.abs_diff(sr) <= delta_local {
                    continue;
                }
                let cr = char_vecs[r];
                let cbr = char_bar[r];
                let q11 = (ck[0] & cr[0]).count_ones() + (ck[1] & cr[1]).count_ones();
                let q10 = (ck[0] & cbr[0]).count_ones() + (ck[1] & cbr[1]).count_ones();
                let q01 = (cbk[0] & cr[0]).count_ones() + (cbk[1] & cr[1]).count_ones();
                let q00 = (cbk[0] & cbr[0]).count_ones() + (cbk[1] & cbr[1]).count_ones();

                let viol = q11.min(q10).min(q01.min(q00)) as f64;
                let compat = (-1.1 * viol).exp();
                sum_long += compat;
                cnt_long += 1;
                if q11 >= 1 {
                    sum_overlap += compat;
                    cnt_overlap += 1;
                }
            }
            coh_scores[k] = if cnt_overlap >= 3 {
                sum_overlap / (cnt_overlap as f64)
            } else if cnt_long >= 1 {
                sum_long / (cnt_long as f64)
            } else {
                1.0
            };
        }

        let mut site_coh = vec![1.0_f64; l];
        let mut site_cnt = vec![0_u32; l];
        for (k, &s) in char_sites.iter().enumerate() {
            if site_cnt[s] == 0 {
                site_coh[s] = coh_scores[k];
            } else {
                site_coh[s] = site_coh[s].max(coh_scores[k]);
            }
            site_cnt[s] += 1;
        }

        let mut inf_vals: Vec<f64> = (0..l)
            .filter(|&s| site_cnt[s] > 0)
            .map(|s| site_coh[s])
            .collect();
        if !inf_vals.is_empty() {
            let med_coh = median_f64(&mut inf_vals);
            if med_coh > 1e-6 {
                for s in 0..l {
                    if site_cnt[s] > 0 {
                        site_coh[s] = (site_coh[s] / med_coh).clamp(0.65, 1.30);
                    }
                }
            }
        }

        let var_count = n_states.iter().filter(|&&st| st >= 2).count();
        let multi_count = n_states.iter().filter(|&&st| st >= 3).count();
        let frac_var = (var_count as f64) / (l as f64);
        let frac_multi = if var_count > 0 {
            (multi_count as f64) / (var_count as f64)
        } else {
            0.0
        };

        let mut mod_dist = vec![0.0_f64; l];
        let mut mod_tree = vec![0.0_f64; l];
        let deep_coding = frac_var >= 0.25 && frac_multi > 0.40;

        for s in 0..l {
            let w_mt = if frac_multi > 0.50 { 1.0 } else { w_mult[s] };
            let coh_tree_factor = (site_coh[s] * w_mt).clamp(0.75, 1.25);
            let (eff_d, eff_t) = if deep_coding && coh_tree_factor >= 0.90 {
                (1.0, 1.0)
            } else {
                (dist_weights[s], tree_weights[s])
            };
            mod_dist[s] = eff_d * site_coh[s].powf(1.25) * w_mult[s];
            mod_tree[s] = eff_t * coh_tree_factor;
        }
        (mod_dist, mod_tree)
    }

    #[inline]
    fn compute_pair_jc_dist(
        row_i: &[i8],
        row_j: &[i8],
        site_weights: &[f64],
        l_eff: f64,
    ) -> f64 {
        let mut w_valid = 0.0_f64;
        let mut w_diffs = 0.0_f64;
        for s in 0..row_i.len() {
            let bi = row_i[s];
            let bj = row_j[s];
            if bi >= 0 && bj >= 0 {
                let w = site_weights[s];
                w_valid += w;
                if bi != bj {
                    w_diffs += w;
                }
            }
        }
        let p_hat = (w_diffs / w_valid.max(1.0)).clamp(0.0, 0.745);
        -0.75 * (1.0 - (4.0 / 3.0) * p_hat).ln() * l_eff
    }

    #[inline]
    fn compute_pair_unweighted_jc_dist(row_i: &[i8], row_j: &[i8], l_eff: f64) -> f64 {
        let mut n_valid = 0usize;
        let mut n_diffs = 0usize;
        for s in 0..row_i.len() {
            let bi = row_i[s];
            let bj = row_j[s];
            if bi >= 0 && bj >= 0 {
                n_valid += 1;
                if bi != bj {
                    n_diffs += 1;
                }
            }
        }
        let p_hat = ((n_diffs as f64) / (n_valid.max(1) as f64)).clamp(0.0, 0.745);
        -0.75 * (1.0 - (4.0 / 3.0) * p_hat).ln() * l_eff
    }

    fn partition_items_recursive(
        items: Vec<usize>,
        seqs: &HashMap<usize, Vec<i8>>,
        site_weights: &[f64],
        l_eff: f64,
        max_cluster_size: usize,
        out_clusters: &mut Vec<Vec<usize>>,
    ) {
        if items.len() <= max_cluster_size {
            out_clusters.push(items);
            return;
        }
        let p0 = items[0];
        let seq_p0 = &seqs[&p0];
        let mut p1 = p0;
        let mut max_d1 = -1.0_f64;
        for &x in &items {
            if x != p0 {
                let d = Self::compute_pair_jc_dist(seq_p0, &seqs[&x], site_weights, l_eff);
                if d > max_d1 {
                    max_d1 = d;
                    p1 = x;
                }
            }
        }
        let seq_p1 = &seqs[&p1];
        let mut p2 = p0;
        let mut max_d2 = -1.0_f64;
        for &x in &items {
            if x != p1 {
                let d = Self::compute_pair_jc_dist(seq_p1, &seqs[&x], site_weights, l_eff);
                if d > max_d2 {
                    max_d2 = d;
                    p2 = x;
                }
            }
        }
        let seq_p2 = &seqs[&p2];
        let mut projections: Vec<(f64, usize)> = items
            .into_iter()
            .map(|x| {
                let sx = &seqs[&x];
                let d1 = Self::compute_pair_jc_dist(sx, seq_p1, site_weights, l_eff);
                let d2 = Self::compute_pair_jc_dist(sx, seq_p2, site_weights, l_eff);
                (d1 - d2, x)
            })
            .collect();
        projections.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.1.cmp(&b.1))
        });
        let mid = projections.len() / 2;
        let left: Vec<usize> = projections[..mid].iter().map(|x| x.1).collect();
        let right: Vec<usize> = projections[mid..].iter().map(|x| x.1).collect();
        Self::partition_items_recursive(left, seqs, site_weights, l_eff, max_cluster_size, out_clusters);
        Self::partition_items_recursive(right, seqs, site_weights, l_eff, max_cluster_size, out_clusters);
    }

    #[allow(clippy::too_many_arguments)]
    fn agglomerate_bounded_cluster(
        &self,
        cluster_nodes: &[usize],
        seqs: &mut HashMap<usize, Vec<i8>>,
        site_weights: &[f64],
        l: usize,
        l_eff: f64,
        adj: &mut BTreeMap<usize, Vec<usize>>,
        node_labels: &mut BTreeMap<usize, String>,
        next_node: &mut usize,
        stop_at_one: bool,
    ) -> usize {
        let m = cluster_nodes.len();
        if m == 1 {
            return cluster_nodes[0];
        }
        if m == 2 {
            let u = cluster_nodes[0];
            let v = cluster_nodes[1];
            if stop_at_one {
                let w = *next_node;
                *next_node += 1;
                adj.entry(w).or_default().extend([u, v]);
                adj.entry(u).or_default().push(w);
                adj.entry(v).or_default().push(w);
                node_labels.insert(w, format!("N{w}"));
                let su = &seqs[&u];
                let sv = &seqs[&v];
                let sw: Vec<i8> = (0..l)
                    .map(|s| if su[s] >= 0 { su[s] } else { sv[s] })
                    .collect();
                seqs.insert(w, sw);
                return w;
            } else {
                adj.entry(u).or_default().push(v);
                adj.entry(v).or_default().push(u);
                return u;
            }
        }

        let cap = 2 * m + 2;
        let mut d_mat = vec![0.0_f64; cap * cap];
        let mut v_mat = vec![0.0_f64; cap * cap];
        let mut local_to_global: Vec<usize> = cluster_nodes.to_vec();

        for i in 0..m {
            let si = &seqs[&local_to_global[i]];
            for j in (i + 1)..m {
                let sj = &seqs[&local_to_global[j]];
                let d_val = Self::compute_pair_jc_dist(si, sj, site_weights, l_eff);
                let v_val = d_val.max(0.5) * (1.0 + (d_val / l_eff.max(1.0)).powi(2));
                d_mat[i * cap + j] = d_val;
                d_mat[j * cap + i] = d_val;
                v_mat[i * cap + j] = v_val;
                v_mat[j * cap + i] = v_val;
            }
        }

        let mut active: Vec<usize> = (0..m).collect();
        let mut row_sums = vec![0.0_f64; cap];
        for i in 0..m {
            let mut s = 0.0_f64;
            for j in 0..m {
                s += d_mat[i * cap + j];
            }
            row_sums[i] = s;
        }

        let min_active = if stop_at_one { 2 } else { 3 };
        while active.len() > min_active {
            let k = active.len();
            let mut best_q = f64::INFINITY;
            let mut best_pair = (0usize, 1usize);
            let k_m2 = ((k as f64) - 2.0).max(1.0);

            for a in 0..k {
                let la = active[a];
                let ra = row_sums[la];
                for b in (a + 1)..k {
                    let lb = active[b];
                    let q_val = k_m2 * d_mat[la * cap + lb] - ra - row_sums[lb];
                    if q_val < best_q {
                        best_q = q_val;
                        best_pair = (a, b);
                    }
                }
            }

            let (a_pos, b_pos) = best_pair;
            let lu = active[a_pos];
            let lv = active[b_pos];
            let gu = local_to_global[lu];
            let gv = local_to_global[lv];

            let gw = *next_node;
            *next_node += 1;
            let lw = local_to_global.len();
            local_to_global.push(gw);

            adj.entry(gw).or_default().extend([gu, gv]);
            adj.entry(gu).or_default().push(gw);
            adj.entry(gv).or_default().push(gw);
            node_labels.insert(gw, format!("N{gw}"));

            if stop_at_one {
                let su = &seqs[&gu];
                let sv = &seqs[&gv];
                let sw: Vec<i8> = (0..l)
                    .map(|s| if su[s] >= 0 { su[s] } else { sv[s] })
                    .collect();
                seqs.insert(gw, sw);
            }

            let d_uv = d_mat[lu * cap + lv];
            let v_uv = v_mat[lu * cap + lv].max(1e-6);
            let others: Vec<usize> = active
                .iter()
                .copied()
                .filter(|&x| x != lu && x != lv)
                .collect();
            let lam = if !others.is_empty() {
                let lam_num: f64 = others
                    .iter()
                    .map(|&x| v_mat[lv * cap + x] - v_mat[lu * cap + x])
                    .sum();
                (0.5 + lam_num / (2.0 * (others.len() as f64) * v_uv)).clamp(0.05, 0.95)
            } else {
                0.5
            };

            let limb_u = (0.5 * d_uv + (row_sums[lu] - row_sums[lv]) / (2.0 * k_m2)).max(0.0);
            let limb_v = (d_uv - limb_u).max(0.0);

            let mut rw_sum = 0.0_f64;
            for &lx in &others {
                let d_ux = d_mat[lu * cap + lx];
                let d_vx = d_mat[lv * cap + lx];
                let d_wx = (lam * (d_ux - limb_u) + (1.0 - lam) * (d_vx - limb_v)).max(0.0);
                let v_wx = (lam * v_mat[lu * cap + lx]
                    + (1.0 - lam) * v_mat[lv * cap + lx]
                    - lam * (1.0 - lam) * v_uv)
                    .max(0.25);
                d_mat[lw * cap + lx] = d_wx;
                d_mat[lx * cap + lw] = d_wx;
                v_mat[lw * cap + lx] = v_wx;
                v_mat[lx * cap + lw] = v_wx;
                row_sums[lx] += d_wx - d_ux - d_vx;
                rw_sum += d_wx;
            }
            row_sums[lw] = rw_sum;

            active.retain(|&x| x != lu && x != lv);
            active.push(lw);
        }

        if stop_at_one {
            // Merge the final 2 active nodes into the cluster root
            let gu = local_to_global[active[0]];
            let gv = local_to_global[active[1]];
            let gw = *next_node;
            *next_node += 1;
            adj.entry(gw).or_default().extend([gu, gv]);
            adj.entry(gu).or_default().push(gw);
            adj.entry(gv).or_default().push(gw);
            node_labels.insert(gw, format!("N{gw}"));
            let su = &seqs[&gu];
            let sv = &seqs[&gv];
            let sw: Vec<i8> = (0..l)
                .map(|s| if su[s] >= 0 { su[s] } else { sv[s] })
                .collect();
            seqs.insert(gw, sw);
            gw
        } else {
            // Top-level unrooted join of final 3 active nodes
            let gw = *next_node;
            *next_node += 1;
            node_labels.insert(gw, format!("N{gw}"));
            for &lx in &active {
                let gx = local_to_global[lx];
                adj.entry(gw).or_default().push(gx);
                adj.entry(gx).or_default().push(gw);
            }
            gw
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn build_subquadratic_loads_and_backbone(
        &self,
        names: &[String],
        mat: &[i8],
        bitmasks: &[u8],
        dist_weights: &[f64],
        tree_weights: &[f64],
        n: usize,
        l: usize,
        l_eff: f64,
        ti_cost: f64,
        tv_cost: f64,
    ) -> (Vec<f64>, BTreeMap<usize, Vec<usize>>, BTreeMap<usize, String>) {
        let mut node_labels: BTreeMap<usize, String> =
            (0..n).map(|i| (i, names[i].clone())).collect();
        let mut adj: BTreeMap<usize, Vec<usize>> = (0..n).map(|i| (i, Vec::new())).collect();

        if n == 1 {
            return (vec![0.5; 1], adj, node_labels);
        }
        if n == 2 {
            let d01 = Self::compute_pair_unweighted_jc_dist(&mat[0..l], &mat[l..2 * l], l_eff);
            let ld = (d01 * 0.5).max(0.5);
            adj.get_mut(&0).unwrap().push(1);
            adj.get_mut(&1).unwrap().push(0);
            return (vec![ld, ld], adj, node_labels);
        }
        if n <= 128 {
            // Exact small-N path (N <= 128) with incremental O(N^2) agglomeration
            let mut uw_mat = vec![0.0_f64; n * n];
            for i in 0..n {
                let ri = &mat[i * l..(i + 1) * l];
                for j in (i + 1)..n {
                    let rj = &mat[j * l..(j + 1) * l];
                    let d = Self::compute_pair_unweighted_jc_dist(ri, rj, l_eff);
                    uw_mat[i * n + j] = d;
                    uw_mat[j * n + i] = d;
                }
            }
            let mut loads = vec![0.0_f64; n];
            for i in 0..n {
                let mut j_near = if i == 0 { 1 } else { 0 };
                let mut d_min = uw_mat[i * n + j_near];
                for j in 0..n {
                    if j != i && uw_mat[i * n + j] < d_min {
                        d_min = uw_mat[i * n + j];
                        j_near = j;
                    }
                }
                let mut deltas: Vec<f64> = (0..n)
                    .filter(|&k| k != i && k != j_near)
                    .map(|k| uw_mat[i * n + k] - uw_mat[j_near * n + k])
                    .collect();
                let med = median_f64(&mut deltas);
                loads[i] = (0.5 * (d_min + med)).clamp(0.4, d_min.max(0.8));
            }

            if n == 3 {
                let center = n;
                adj.insert(center, vec![0, 1, 2]);
                node_labels.insert(center, format!("N{center}"));
                for i in 0..3 {
                    adj.get_mut(&i).unwrap().push(center);
                }
                return (loads, adj, node_labels);
            }
            if n == 4 {
                let t: Vec<(Vec<u8>, f64)> = (0..4)
                    .map(|i| (bitmasks[i * l..(i + 1) * l].to_vec(), loads[i]))
                    .collect();
                let splits = [([0, 1], [2, 3]), ([0, 2], [1, 3]), ([0, 3], [1, 2])];
                let mut best_idx = 0;
                let mut best_gain = f64::NEG_INFINITY;
                for (idx, (p1, p2)) in splits.iter().enumerate() {
                    let (g, _, _) = self.compute_split_mdl_gain_from_children(
                        &[t[p1[0]].clone(), t[p1[1]].clone()],
                        &[t[p2[0]].clone(), t[p2[1]].clone()],
                        l,
                        l_eff,
                        ti_cost,
                        tv_cost,
                        tree_weights,
                    );
                    if g > best_gain {
                        best_gain = g;
                        best_idx = idx;
                    }
                }
                if best_gain > 0.0 {
                    let (p1, p2) = splits[best_idx];
                    let (u, v) = (4, 5);
                    adj.insert(u, vec![p1[0], p1[1], v]);
                    adj.insert(v, vec![p2[0], p2[1], u]);
                    node_labels.insert(u, format!("N{u}"));
                    node_labels.insert(v, format!("N{v}"));
                    for &x in &p1 {
                        adj.get_mut(&x).unwrap().push(u);
                    }
                    for &x in &p2 {
                        adj.get_mut(&x).unwrap().push(v);
                    }
                } else {
                    let center = 4;
                    adj.insert(center, vec![0, 1, 2, 3]);
                    node_labels.insert(center, format!("N{center}"));
                    for i in 0..4 {
                        adj.get_mut(&i).unwrap().push(center);
                    }
                }
                return (loads, adj, node_labels);
            }

            let mut seqs_map: HashMap<usize, Vec<i8>> = (0..n)
                .map(|i| (i, mat[i * l..(i + 1) * l].to_vec()))
                .collect();
            let all_nodes: Vec<usize> = (0..n).collect();
            let mut next_node = n;
            self.agglomerate_bounded_cluster(
                &all_nodes,
                &mut seqs_map,
                dist_weights,
                l,
                l_eff,
                &mut adj,
                &mut node_labels,
                &mut next_node,
                false,
            );
            self.refine_backbone_nni_mdl(
                &mut adj,
                bitmasks,
                &loads,
                n,
                l,
                l_eff,
                ti_cost,
                tv_cost,
                tree_weights,
                5,
            );
            return (loads, adj, node_labels);
        }

        // Sub-Quadratic O(N log N * L) Hierarchical Tournament Path for N > 128
        let mut seqs_map: HashMap<usize, Vec<i8>> = (0..n)
            .map(|i| (i, mat[i * l..(i + 1) * l].to_vec()))
            .collect();
        let mut leaf_clusters: Vec<Vec<usize>> = Vec::new();
        Self::partition_items_recursive(
            (0..n).collect(),
            &seqs_map,
            dist_weights,
            l_eff,
            64,
            &mut leaf_clusters,
        );

        // Compute 3-point terminal loads within each leaf cluster + 32 global anchors in O(N * L)
        let n_anc = 32.min(n);
        let global_anchors: Vec<usize> = (0..n_anc).map(|i| (i * (n - 1)) / (n_anc - 1)).collect();
        let mut loads = vec![0.5_f64; n];
        for cluster in &leaf_clusters {
            let mut ref_set: Vec<usize> = cluster.clone();
            for &ga in &global_anchors {
                if !ref_set.contains(&ga) {
                    ref_set.push(ga);
                }
            }
            let m_ref = ref_set.len();
            let mut dist_to_ref = vec![0.0_f64; m_ref * m_ref];
            for a in 0..m_ref {
                let sa = &seqs_map[&ref_set[a]];
                for b in (a + 1)..m_ref {
                    let sb = &seqs_map[&ref_set[b]];
                    let d = Self::compute_pair_unweighted_jc_dist(sa, sb, l_eff);
                    dist_to_ref[a * m_ref + b] = d;
                    dist_to_ref[b * m_ref + a] = d;
                }
            }
            for (loc_i, &gi) in cluster.iter().enumerate() {
                let mut best_b = if loc_i == 0 { 1 } else { 0 };
                let mut d_min = dist_to_ref[loc_i * m_ref + best_b];
                for b in 0..m_ref {
                    if b != loc_i && dist_to_ref[loc_i * m_ref + b] < d_min {
                        d_min = dist_to_ref[loc_i * m_ref + b];
                        best_b = b;
                    }
                }
                let mut deltas: Vec<f64> = (0..m_ref)
                    .filter(|&k| k != loc_i && k != best_b)
                    .map(|k| dist_to_ref[loc_i * m_ref + k] - dist_to_ref[best_b * m_ref + k])
                    .collect();
                let med = median_f64(&mut deltas);
                loads[gi] = (0.5 * (d_min + med)).clamp(0.4, d_min.max(0.8));
            }
        }

        // Hierarchical bottom-up cluster agglomeration into a single degree-3 binary backbone
        let mut next_node = n;
        let mut current_roots: Vec<usize> = Vec::with_capacity(leaf_clusters.len());
        for cluster in &leaf_clusters {
            let r = self.agglomerate_bounded_cluster(
                cluster,
                &mut seqs_map,
                dist_weights,
                l,
                l_eff,
                &mut adj,
                &mut node_labels,
                &mut next_node,
                true,
            );
            current_roots.push(r);
        }

        while current_roots.len() > 64 {
            let mut super_clusters: Vec<Vec<usize>> = Vec::new();
            Self::partition_items_recursive(
                current_roots,
                &seqs_map,
                dist_weights,
                l_eff,
                64,
                &mut super_clusters,
            );
            let mut next_roots = Vec::with_capacity(super_clusters.len());
            for sc in &super_clusters {
                let r = self.agglomerate_bounded_cluster(
                    sc,
                    &mut seqs_map,
                    dist_weights,
                    l,
                    l_eff,
                    &mut adj,
                    &mut node_labels,
                    &mut next_node,
                    true,
                );
                next_roots.push(r);
            }
            current_roots = next_roots;
        }

        self.agglomerate_bounded_cluster(
            &current_roots,
            &mut seqs_map,
            dist_weights,
            l,
            l_eff,
            &mut adj,
            &mut node_labels,
            &mut next_node,
            false,
        );

        let nni_rounds = if n <= 256 { 5 } else { 3 };
        self.refine_backbone_nni_mdl(
            &mut adj,
            bitmasks,
            &loads,
            n,
            l,
            l_eff,
            ti_cost,
            tv_cost,
            tree_weights,
            nni_rounds,
        );
        (loads, adj, node_labels)
    }

    fn propagate_directed_fitch_outward(
        adj: &BTreeMap<usize, Vec<usize>>,
        states: &mut HashMap<(usize, usize), (Vec<u8>, f64)>,
        mut queue: VecDeque<(usize, usize, u8)>,
        n: usize,
    ) {
        while let Some((x, y, depth)) = queue.pop_front() {
            if x < n {
                continue;
            }
            let children: Vec<(Vec<u8>, f64)> = adj[&x]
                .iter()
                .copied()
                .filter(|&nb| nb != y)
                .map(|nb| states[&(nb, x)].clone())
                .collect();
            let new_state = Self::combine_fitch_children(&children);
            let changed = match states.get(&(x, y)) {
                Some((old_mask, old_load)) => {
                    new_state.0 != *old_mask || (new_state.1 - *old_load).abs() > 1e-4
                }
                None => true,
            };
            states.insert((x, y), new_state);
            let max_d = if n <= 128 { 12 } else { 3 };
            if changed && depth < max_d && y >= n {
                for &z in &adj[&y] {
                    if z != x {
                        queue.push_back((y, z, depth + 1));
                    }
                }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn refine_backbone_nni_mdl(
        &self,
        adj: &mut BTreeMap<usize, Vec<usize>>,
        bitmasks: &[u8],
        terminal_loads: &[f64],
        n: usize,
        l: usize,
        l_eff: f64,
        ti_cost: f64,
        tv_cost: f64,
        site_weights: &[f64],
        max_rounds: usize,
    ) {
        let mut states =
            self.compute_all_directed_fitch_states(adj, bitmasks, n, l, terminal_loads, l_eff);

        for round_idx in 0..max_rounds {
            let mut improved = false;

            let mut internal_edges = Vec::new();
            for (&u, nbs) in adj.iter() {
                if u >= n && nbs.len() == 3 {
                    for &v in nbs {
                        if u < v && v >= n && adj.get(&v).map(|x| x.len()) == Some(3) {
                            internal_edges.push((u, v));
                        }
                    }
                }
            }

            for (u, v) in internal_edges {
                if !adj.get(&u).map(|x| x.contains(&v)).unwrap_or(false)
                    || adj[&u].len() != 3
                    || adj[&v].len() != 3
                {
                    continue;
                }
                let u_nbs: Vec<usize> = adj[&u].iter().copied().filter(|&x| x != v).collect();
                let v_nbs: Vec<usize> = adj[&v].iter().copied().filter(|&y| y != u).collect();
                let (a, b) = (u_nbs[0], u_nbs[1]);
                let (c, d) = (v_nbs[0], v_nbs[1]);

                let ta = states[&(a, u)].clone();
                let tb = states[&(b, u)].clone();
                let tc = states[&(c, v)].clone();
                let td = states[&(d, v)].clone();

                let (g0, _, _) = self.compute_split_mdl_gain_from_children(
                    &[ta.clone(), tb.clone()],
                    &[tc.clone(), td.clone()],
                    l,
                    l_eff,
                    ti_cost,
                    tv_cost,
                    site_weights,
                );
                let (g1, _, _) = self.compute_split_mdl_gain_from_children(
                    &[ta.clone(), tc.clone()],
                    &[tb.clone(), td.clone()],
                    l,
                    l_eff,
                    ti_cost,
                    tv_cost,
                    site_weights,
                );
                let (g2, _, _) = self.compute_split_mdl_gain_from_children(
                    &[ta, td],
                    &[tb.clone(), tc],
                    l,
                    l_eff,
                    ti_cost,
                    tv_cost,
                    site_weights,
                );

                let best_g = g0.max(g1).max(g2);
                if best_g > g0 + 0.10 {
                    let (swap_v, keep_v) = if g1 >= g2 { (c, d) } else { (d, c) };
                    let t_swap_v = states[&(swap_v, v)].clone();

                    adj.get_mut(&u).unwrap().retain(|&x| x != b);
                    adj.get_mut(&b).unwrap().retain(|&x| x != u);
                    adj.get_mut(&v).unwrap().retain(|&x| x != swap_v);
                    adj.get_mut(&swap_v).unwrap().retain(|&x| x != v);

                    adj.get_mut(&u).unwrap().push(swap_v);
                    adj.get_mut(&swap_v).unwrap().push(u);
                    adj.get_mut(&v).unwrap().push(b);
                    adj.get_mut(&b).unwrap().push(v);

                    states.remove(&(b, u));
                    states.remove(&(u, b));
                    states.remove(&(swap_v, v));
                    states.remove(&(v, swap_v));
                    states.insert((b, v), tb);
                    states.insert((swap_v, u), t_swap_v);

                    let st_uv = Self::combine_fitch_children(&[
                        states[&(a, u)].clone(),
                        states[&(swap_v, u)].clone(),
                    ]);
                    let st_vu = Self::combine_fitch_children(&[
                        states[&(keep_v, v)].clone(),
                        states[&(b, v)].clone(),
                    ]);
                    states.insert((u, v), st_uv);
                    states.insert((v, u), st_vu);

                    let q = VecDeque::from([(u, a, 0), (u, swap_v, 0), (v, keep_v, 0), (v, b, 0)]);
                    Self::propagate_directed_fitch_outward(adj, &mut states, q, n);
                    improved = true;
                }
            }

            // 2-Hop Cross-Node Attractor Swap (Radius-2 SPR across u - w - v)
            if round_idx < 2 && n <= 250 {
                let int_nodes: Vec<usize> = adj
                    .iter()
                    .filter(|(&w, nbs)| w >= n && nbs.len() == 3)
                    .map(|(&w, _)| w)
                    .collect();

                for w in int_nodes {
                    if adj.get(&w).map(|x| x.len()) != Some(3) {
                        continue;
                    }
                    let nbs = adj[&w].clone();
                    let int_nbs: Vec<usize> = nbs
                        .iter()
                        .copied()
                        .filter(|&x| x >= n && adj.get(&x).map(|v| v.len()) == Some(3))
                        .collect();
                    if int_nbs.len() < 2 {
                        continue;
                    }
                    for idx_i in 0..int_nbs.len() {
                        for idx_j in (idx_i + 1)..int_nbs.len() {
                            let u = int_nbs[idx_i];
                            let v = int_nbs[idx_j];
                            if !adj[&u].contains(&w)
                                || !adj[&v].contains(&w)
                                || adj[&u].len() != 3
                                || adj[&v].len() != 3
                            {
                                continue;
                            }
                            let m_list: Vec<usize> =
                                adj[&w].iter().copied().filter(|&x| x != u && x != v).collect();
                            if m_list.len() != 1 {
                                continue;
                            }
                            let m_node = m_list[0];
                            let u_ch: Vec<usize> =
                                adj[&u].iter().copied().filter(|&x| x != w).collect();
                            let v_ch: Vec<usize> =
                                adj[&v].iter().copied().filter(|&y| y != w).collect();
                            if u_ch.len() != 2 || v_ch.len() != 2 {
                                continue;
                            }
                            let (a, b) = (u_ch[0], u_ch[1]);
                            let (c, d) = (v_ch[0], v_ch[1]);

                            let ta = states[&(a, u)].clone();
                            let tb = states[&(b, u)].clone();
                            let tc = states[&(c, v)].clone();
                            let td = states[&(d, v)].clone();
                            let tm = states[&(m_node, w)].clone();

                            let t_u_curr = Self::combine_fitch_children(&[ta.clone(), tb.clone()]);
                            let t_v_curr = Self::combine_fitch_children(&[tc.clone(), td.clone()]);
                            let (g_uw_0, _, _) = self.compute_split_mdl_gain_from_children(
                                &[ta.clone(), tb.clone()],
                                &[tm.clone(), t_v_curr],
                                l,
                                l_eff,
                                ti_cost,
                                tv_cost,
                                site_weights,
                            );
                            let (g_wv_0, _, _) = self.compute_split_mdl_gain_from_children(
                                &[tc.clone(), td.clone()],
                                &[tm.clone(), t_u_curr],
                                l,
                                l_eff,
                                ti_cost,
                                tv_cost,
                                site_weights,
                            );
                            let base_2hop = g_uw_0 + g_wv_0;
                            let mut best_swap_gain = base_2hop + 0.50;
                            let mut best_pair: Option<(usize, usize)> = None;

                            for (u_swap, t_us, t_uk) in
                                [(a, ta.clone(), tb.clone()), (b, tb.clone(), ta.clone())]
                            {
                                for (v_swap, t_vs, t_vk) in
                                    [(c, tc.clone(), td.clone()), (d, td.clone(), tc.clone())]
                                {
                                    let t_u_new = Self::combine_fitch_children(&[
                                        t_uk.clone(),
                                        t_vs.clone(),
                                    ]);
                                    let t_v_new = Self::combine_fitch_children(&[
                                        t_vk.clone(),
                                        t_us.clone(),
                                    ]);
                                    let (g_uw_n, _, _) = self.compute_split_mdl_gain_from_children(
                                        &[t_uk.clone(), t_vs],
                                        &[tm.clone(), t_v_new],
                                        l,
                                        l_eff,
                                        ti_cost,
                                        tv_cost,
                                        site_weights,
                                    );
                                    let (g_wv_n, _, _) = self.compute_split_mdl_gain_from_children(
                                        &[t_vk, t_us.clone()],
                                        &[tm.clone(), t_u_new],
                                        l,
                                        l_eff,
                                        ti_cost,
                                        tv_cost,
                                        site_weights,
                                    );
                                    if (g_uw_n + g_wv_n) > best_swap_gain {
                                        best_swap_gain = g_uw_n + g_wv_n;
                                        best_pair = Some((u_swap, v_swap));
                                    }
                                }
                            }

                            if let Some((u_swap, v_swap)) = best_pair {
                                let u_keep = if u_swap == a { b } else { a };
                                let v_keep = if v_swap == c { d } else { c };
                                let t_us = states[&(u_swap, u)].clone();
                                let t_vs = states[&(v_swap, v)].clone();

                                adj.get_mut(&u).unwrap().retain(|&x| x != u_swap);
                                adj.get_mut(&u_swap).unwrap().retain(|&x| x != u);
                                adj.get_mut(&v).unwrap().retain(|&x| x != v_swap);
                                adj.get_mut(&v_swap).unwrap().retain(|&x| x != v);

                                adj.get_mut(&u).unwrap().push(v_swap);
                                adj.get_mut(&v_swap).unwrap().push(u);
                                adj.get_mut(&v).unwrap().push(u_swap);
                                adj.get_mut(&u_swap).unwrap().push(v);

                                states.remove(&(u_swap, u));
                                states.remove(&(u, u_swap));
                                states.remove(&(v_swap, v));
                                states.remove(&(v, v_swap));
                                states.insert((u_swap, v), t_us);
                                states.insert((v_swap, u), t_vs);

                                let st_uw = Self::combine_fitch_children(&[
                                    states[&(u_keep, u)].clone(),
                                    states[&(v_swap, u)].clone(),
                                ]);
                                let st_vw = Self::combine_fitch_children(&[
                                    states[&(v_keep, v)].clone(),
                                    states[&(u_swap, v)].clone(),
                                ]);
                                states.insert((u, w), st_uw.clone());
                                states.insert((v, w), st_vw.clone());

                                let st_wu = Self::combine_fitch_children(&[
                                    states[&(m_node, w)].clone(),
                                    st_vw,
                                ]);
                                let st_wv = Self::combine_fitch_children(&[
                                    states[&(m_node, w)].clone(),
                                    st_uw,
                                ]);
                                states.insert((w, u), st_wu);
                                states.insert((w, v), st_wv);

                                let q = VecDeque::from([
                                    (u, u_keep, 0),
                                    (u, v_swap, 0),
                                    (v, v_keep, 0),
                                    (v, u_swap, 0),
                                    (w, m_node, 0),
                                ]);
                                Self::propagate_directed_fitch_outward(adj, &mut states, q, n);
                                improved = true;
                            }
                        }
                    }
                }
            }

            if !improved {
                break;
            }
        }
    }

    pub fn compute_all_directed_fitch_states(
        &self,
        adj: &BTreeMap<usize, Vec<usize>>,
        bitmasks: &[u8],
        n: usize,
        l: usize,
        terminal_loads: &[f64],
        _l_eff: f64,
    ) -> HashMap<(usize, usize), (Vec<u8>, f64)> {
        let root = adj
            .iter()
            .find(|(&u, nbs)| u >= n && !nbs.is_empty())
            .or_else(|| adj.iter().find(|(_, nbs)| !nbs.is_empty()))
            .map(|(&u, _)| u)
            .unwrap_or(0);

        let mut parent: HashMap<usize, Option<usize>> = HashMap::new();
        parent.insert(root, None);
        let mut order = vec![root];
        let mut q = VecDeque::new();
        q.push_back(root);

        while let Some(curr) = q.pop_front() {
            if let Some(nbs) = adj.get(&curr) {
                for &nb in nbs {
                    if !parent.contains_key(&nb) {
                        parent.insert(nb, Some(curr));
                        order.push(nb);
                        q.push_back(nb);
                    }
                }
            }
        }

        let mut states: HashMap<(usize, usize), (Vec<u8>, f64)> =
            HashMap::with_capacity(order.len() * 2);

        // Bottom-up pass
        for &u in order.iter().rev() {
            let Some(Some(p)) = parent.get(&u) else {
                continue;
            };
            let p = *p;
            if u < n {
                states.insert(
                    (u, p),
                    (bitmasks[u * l..(u + 1) * l].to_vec(), terminal_loads[u]),
                );
            } else {
                let children: Vec<(Vec<u8>, f64)> = adj[&u]
                    .iter()
                    .copied()
                    .filter(|&nb| nb != p)
                    .map(|c| states[&(c, u)].clone())
                    .collect();
                states.insert((u, p), Self::combine_fitch_children(&children));
            }
        }

        // Top-down pass (with O(d * L) leave-one-out hub acceleration for degree d >= 4)
        for &u in &order {
            let p_opt = parent[&u];
            let children: Vec<usize> = adj[&u]
                .iter()
                .copied()
                .filter(|&nb| Some(nb) != p_opt)
                .collect();
            if u < n {
                for c in children {
                    states.insert(
                        (u, c),
                        (bitmasks[u * l..(u + 1) * l].to_vec(), terminal_loads[u]),
                    );
                }
            } else if adj[&u].len() >= 4 {
                let all_nbs = &adj[&u];
                let mut tot_inv = 0.0_f64;
                let mut tot_counts = vec![[0_u16; 4]; l];
                for &nb in all_nbs {
                    let (ref s_mask, m_c) = states[&(nb, u)];
                    tot_inv += 1.0 / m_c.max(0.4);
                    for s in 0..l {
                        let m = s_mask[s];
                        if m != 15 {
                            tot_counts[s][0] += (m & 1) as u16;
                            tot_counts[s][1] += ((m >> 1) & 1) as u16;
                            tot_counts[s][2] += ((m >> 2) & 1) as u16;
                            tot_counts[s][3] += ((m >> 3) & 1) as u16;
                        }
                    }
                }
                for c in children {
                    let (s_c, m_c) = states[&(c, u)].clone();
                    let inv_rem = (tot_inv - 1.0 / m_c.max(0.4)).max(1e-6);
                    let par_load = (1.0 / inv_rem).max(0.4);
                    let mut out_mask = vec![15_u8; l];
                    for s in 0..l {
                        let mc = s_c[s];
                        let mut c0 = tot_counts[s][0];
                        let mut c1 = tot_counts[s][1];
                        let mut c2 = tot_counts[s][2];
                        let mut c3 = tot_counts[s][3];
                        if mc != 15 {
                            c0 -= (mc & 1) as u16;
                            c1 -= ((mc >> 1) & 1) as u16;
                            c2 -= ((mc >> 2) & 1) as u16;
                            c3 -= ((mc >> 3) & 1) as u16;
                        }
                        let max_c = c0.max(c1).max(c2.max(c3));
                        if max_c > 0 {
                            let mut mask = 0_u8;
                            if c0 == max_c {
                                mask |= 1;
                            }
                            if c1 == max_c {
                                mask |= 2;
                            }
                            if c2 == max_c {
                                mask |= 4;
                            }
                            if c3 == max_c {
                                mask |= 8;
                            }
                            out_mask[s] = mask;
                        }
                    }
                    states.insert((u, c), (out_mask, par_load));
                }
            } else {
                for c in children {
                    let other_nbs: Vec<(Vec<u8>, f64)> = adj[&u]
                        .iter()
                        .copied()
                        .filter(|&nb| nb != c)
                        .map(|nb| states[&(nb, u)].clone())
                        .collect();
                    if other_nbs.is_empty() {
                        states.insert((u, c), (vec![15_u8; l], 1.0));
                    } else {
                        states.insert((u, c), Self::combine_fitch_children(&other_nbs));
                    }
                }
            }
        }

        states
    }

    #[inline]
    pub fn combine_fitch_children(child_list: &[(Vec<u8>, f64)]) -> (Vec<u8>, f64) {
        if child_list.len() == 1 {
            return child_list[0].clone();
        }
        let l = child_list[0].0.len();
        if child_list.len() == 2 {
            let (s1, m1) = &child_list[0];
            let (s2, m2) = &child_list[1];
            let mut out = vec![0_u8; l];
            let mut n_mut = 0usize;
            for i in 0..l {
                let b1 = s1[i];
                let b2 = s2[i];
                let inter = b1 & b2;
                if inter == 0 {
                    out[i] = b1 | b2;
                    if b1 != 15 && b2 != 15 {
                        n_mut += 1;
                    }
                } else {
                    out[i] = inter;
                }
            }
            let par_load =
                1.0 / (1.0 / m1.max(0.4) + 1.0 / m2.max(0.4)) + 0.25 * (n_mut as f64);
            return (out, par_load.max(0.4));
        }

        let mut out = vec![15_u8; l];
        let mut inv_sum = 0.0_f64;
        for (_, m_c) in child_list {
            inv_sum += 1.0 / m_c.max(0.4);
        }
        for s in 0..l {
            let mut counts = [0_u16; 4];
            for (s_mask, _) in child_list {
                let m = s_mask[s];
                if m != 15 {
                    counts[0] += (m & 1) as u16;
                    counts[1] += ((m >> 1) & 1) as u16;
                    counts[2] += ((m >> 2) & 1) as u16;
                    counts[3] += ((m >> 3) & 1) as u16;
                }
            }
            let max_c = counts[0].max(counts[1]).max(counts[2].max(counts[3]));
            if max_c > 0 {
                let mut mask = 0_u8;
                for b in 0..4 {
                    if counts[b] == max_c {
                        mask |= 1 << b;
                    }
                }
                out[s] = mask;
            }
        }
        let par_load = (1.0 / inv_sum.max(1e-6)).max(0.4);
        (out, par_load)
    }

    #[allow(clippy::too_many_arguments)]
    fn expand_polytomy_hubs_mdl(
        &self,
        adj: &mut BTreeMap<usize, Vec<usize>>,
        node_labels: &mut BTreeMap<usize, String>,
        bitmasks: &[u8],
        terminal_loads: &[f64],
        n: usize,
        l: usize,
        l_eff: f64,
        ti_cost: f64,
        tv_cost: f64,
        site_weights: &[f64],
    ) {
        let mut next_id = adj.keys().copied().max().unwrap_or(n) + 1;
        let mut states =
            self.compute_all_directed_fitch_states(adj, bitmasks, n, l, terminal_loads, l_eff);
        let hubs: Vec<usize> = adj
            .iter()
            .filter(|(&h, nbs)| h >= n && (4..=40).contains(&nbs.len()))
            .map(|(&h, _)| h)
            .collect();

        for h in hubs {
            let Some(nbs_ref) = adj.get(&h) else {
                continue;
            };
            if nbs_ref.len() < 4 {
                continue;
            }
            let nbs = nbs_ref.clone();
            let d = nbs.len();
            let nb_states: Vec<(Vec<u8>, f64)> =
                nbs.iter().map(|&c| states[&(c, h)].clone()).collect();

            let mut seed_weights: HashMap<Vec<usize>, f64> = HashMap::new();
            for s in 0..l {
                let mut n_val = 0usize;
                for st in &nb_states {
                    if st.0[s] != 15 {
                        n_val += 1;
                    }
                }
                if n_val < 4 {
                    continue;
                }
                for b in 0..4_u8 {
                    let mut idx_vec: Vec<usize> = (0..d)
                        .filter(|&i| {
                            let m = nb_states[i].0[s];
                            m != 15 && ((m >> b) & 1) != 0
                        })
                        .collect();
                    let cnt = idx_vec.len();
                    if cnt >= 2 && cnt <= d - 2 {
                        if cnt > d / 2 || (cnt == d / 2 && !idx_vec.contains(&0)) {
                            idx_vec = (0..d).filter(|i| !idx_vec.contains(i)).collect();
                        }
                        if idx_vec.len() >= 2 && idx_vec.len() <= d - 2 {
                            *seed_weights.entry(idx_vec).or_insert(0.0) += site_weights[s];
                        }
                    }
                }
            }

            if seed_weights.is_empty() {
                continue;
            }

            let mut top_seeds: Vec<(Vec<usize>, f64)> = seed_weights.into_iter().collect();
            top_seeds.sort_by(|a, b| {
                b.1.partial_cmp(&a.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.0.cmp(&b.0))
            });
            top_seeds.truncate(25);

            let mut positive_attractors: Vec<(f64, BTreeSet<usize>)> = Vec::new();
            let mut seen_sets: HashSet<BTreeSet<usize>> = HashSet::new();

            for (idx_vec, _) in top_seeds {
                let mut u_set: BTreeSet<usize> = idx_vec.into_iter().collect();
                let u_ch: Vec<(Vec<u8>, f64)> =
                    u_set.iter().map(|&i| nb_states[i].clone()).collect();
                let v_ch: Vec<(Vec<u8>, f64)> = (0..d)
                    .filter(|j| !u_set.contains(j))
                    .map(|j| nb_states[j].clone())
                    .collect();
                let (mut g_cur, _, sites_cur) = self.compute_split_mdl_gain_from_children(
                    &u_ch,
                    &v_ch,
                    l,
                    l_eff,
                    ti_cost,
                    tv_cost,
                    site_weights,
                );
                if !sites_cur.is_empty() && d <= 20 {
                    for k_try in 0..d {
                        let mut cand_set = u_set.clone();
                        if !cand_set.remove(&k_try) {
                            cand_set.insert(k_try);
                        }
                        if cand_set.len() < 2 || cand_set.len() > d - 2 {
                            continue;
                        }
                        let u_try: Vec<(Vec<u8>, f64)> =
                            cand_set.iter().map(|&i| nb_states[i].clone()).collect();
                        let v_try: Vec<(Vec<u8>, f64)> = (0..d)
                            .filter(|j| !cand_set.contains(j))
                            .map(|j| nb_states[j].clone())
                            .collect();
                        let (g_try, _, _) = self.compute_split_mdl_gain_from_children(
                            &u_try,
                            &v_try,
                            l,
                            l_eff,
                            ti_cost,
                            tv_cost,
                            site_weights,
                        );
                        if g_try > g_cur + 0.25 {
                            u_set = cand_set;
                            g_cur = g_try;
                        }
                    }
                }
                if g_cur > 0.0 && seen_sets.insert(u_set.clone()) {
                    positive_attractors.push((g_cur, u_set));
                }
            }

            if positive_attractors.is_empty() {
                continue;
            }

            positive_attractors
                .sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
            let mut chosen_splits: Vec<BTreeSet<usize>> = Vec::new();
            let mut used_indices: BTreeSet<usize> = BTreeSet::new();
            for (_, u_set) in positive_attractors {
                if u_set.is_disjoint(&used_indices)
                    && (d - used_indices.len() - u_set.len()) >= 2
                {
                    used_indices.extend(&u_set);
                    chosen_splits.push(u_set);
                }
            }

            let had_splits = !chosen_splits.is_empty();
            for u_set in chosen_splits {
                let new_node = next_id;
                next_id += 1;
                adj.insert(new_node, vec![h]);
                node_labels.insert(new_node, format!("N{new_node}"));
                adj.get_mut(&h).unwrap().push(new_node);
                for idx_c in u_set {
                    let c_node = nbs[idx_c];
                    if adj[&h].contains(&c_node) {
                        adj.get_mut(&h).unwrap().retain(|&x| x != c_node);
                        adj.get_mut(&c_node).unwrap().retain(|&x| x != h);
                        adj.get_mut(&new_node).unwrap().push(c_node);
                        adj.get_mut(&c_node).unwrap().push(new_node);
                    }
                }
            }

            if had_splits {
                states = self.compute_all_directed_fitch_states(
                    adj,
                    bitmasks,
                    n,
                    l,
                    terminal_loads,
                    l_eff,
                );
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn contract_unsupported_edges_mdl(
        &self,
        mut adj: BTreeMap<usize, Vec<usize>>,
        bitmasks: &[u8],
        majority_state: &[i8],
        terminal_loads: &[f64],
        n: usize,
        l: usize,
        l_eff: f64,
        ti_cost: f64,
        tv_cost: f64,
        site_weights: &[f64],
        node_labels: &mut BTreeMap<usize, String>,
        expand_hubs: bool,
    ) -> (
        BTreeMap<usize, Vec<usize>>,
        Vec<BranchRecord>,
        HashMap<usize, Vec<i8>>,
        HashMap<(usize, usize), Vec<(usize, i8, i8)>>,
    ) {
        let passes = if expand_hubs { 2 } else { 1 };
        for outer_pass in 0..passes {
            loop {
                let states = self.compute_all_directed_fitch_states(
                    &adj,
                    bitmasks,
                    n,
                    l,
                    terminal_loads,
                    l_eff,
                );
                let mut internal_edges = Vec::new();
                for (&u, nbs) in &adj {
                    if u >= n && nbs.len() >= 3 {
                        for &v in nbs {
                            if u < v && v >= n && adj.get(&v).map(|x| x.len()).unwrap_or(0) >= 3 {
                                internal_edges.push((u, v));
                            }
                        }
                    }
                }
                if internal_edges.is_empty() {
                    break;
                }

                let mut bad_edges: Vec<(f64, usize, usize)> = Vec::new();
                for (u, v) in internal_edges {
                    let (gain, _, _) = self.evaluate_internal_edge_mdl_gain(
                        u,
                        v,
                        &adj,
                        &states,
                        l,
                        l_eff,
                        ti_cost,
                        tv_cost,
                        site_weights,
                    );
                    if gain <= 0.0 {
                        bad_edges.push((gain, u, v));
                    }
                }
                if bad_edges.is_empty() {
                    break;
                }

                bad_edges
                    .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
                let mut locked_nodes: HashSet<usize> = HashSet::new();
                let mut alias: HashMap<usize, usize> = HashMap::new();
                let mut contracted_any = false;

                for (gain, raw_u, raw_v) in bad_edges {
                    let mut u = raw_u;
                    while let Some(&au) = alias.get(&u) {
                        u = au;
                    }
                    let mut v = raw_v;
                    while let Some(&av) = alias.get(&v) {
                        v = av;
                    }
                    if u == v {
                        continue;
                    }
                    if (n <= 128 || gain > -1.0)
                        && (locked_nodes.contains(&u) || locked_nodes.contains(&v))
                    {
                        continue;
                    }
                    if !adj.get(&u).map(|x| x.contains(&v)).unwrap_or(false) {
                        continue;
                    }
                    adj.get_mut(&u).unwrap().retain(|&x| x != v);
                    adj.get_mut(&v).unwrap().retain(|&x| x != u);
                    let v_nbs = adj[&v].clone();
                    for nb in v_nbs {
                        adj.get_mut(&nb).unwrap().retain(|&x| x != v);
                        if nb != u && !adj[&u].contains(&nb) {
                            adj.get_mut(&u).unwrap().push(nb);
                            adj.get_mut(&nb).unwrap().push(u);
                        }
                    }
                    adj.get_mut(&v).unwrap().clear();
                    alias.insert(v, u);
                    locked_nodes.insert(u);
                    locked_nodes.insert(v);
                    contracted_any = true;
                }
                if !contracted_any {
                    break;
                }
            }

            if expand_hubs && outer_pass == 0 && n >= 5 {
                self.expand_polytomy_hubs_mdl(
                    &mut adj,
                    node_labels,
                    bitmasks,
                    terminal_loads,
                    n,
                    l,
                    l_eff,
                    ti_cost,
                    tv_cost,
                    site_weights,
                );
            }
        }

        let states =
            self.compute_all_directed_fitch_states(&adj, bitmasks, n, l, terminal_loads, l_eff);
        let (node_states, edge_mutations) = self.reconstruct_nodal_states_and_mutations(
            &adj,
            bitmasks,
            majority_state,
            &states,
            n,
            l,
        );

        let mut branches = Vec::new();
        for (&u, nbs) in &adj {
            let mut sorted_nbs = nbs.clone();
            sorted_nbs.sort_unstable();
            for v in sorted_nbs {
                if u < v {
                    let is_int = u >= n && v >= n;
                    let m_count = edge_mutations
                        .get(&(u, v))
                        .map(|m| m.len() as f64)
                        .unwrap_or(0.0);
                    let (gain, eff_mut, syn_sites) = if is_int {
                        let (g, m_syn, sites) = self.evaluate_internal_edge_mdl_gain(
                            u,
                            v,
                            &adj,
                            &states,
                            l,
                            l_eff,
                            ti_cost,
                            tv_cost,
                            site_weights,
                        );
                        (g, m_count.max(m_syn), sites)
                    } else {
                        (m_count * 2.0, m_count, Vec::new())
                    };
                    let support = 1.0 / (1.0 + (-gain.clamp(-50.0, 50.0)).exp());
                    branches.push(BranchRecord {
                        u,
                        v,
                        is_internal: is_int,
                        mutation_count: eff_mut,
                        branch_length: eff_mut / l_eff,
                        mdl_gain_nats: gain,
                        mdl_support: support,
                        synapomorphy_sites: syn_sites,
                    });
                }
            }
        }

        (adj, branches, node_states, edge_mutations)
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_internal_edge_mdl_gain(
        &self,
        u: usize,
        v: usize,
        adj: &BTreeMap<usize, Vec<usize>>,
        states: &HashMap<(usize, usize), (Vec<u8>, f64)>,
        l: usize,
        l_eff: f64,
        ti_cost: f64,
        tv_cost: f64,
        site_weights: &[f64],
    ) -> (f64, f64, Vec<usize>) {
        let mut u_sub: Vec<usize> = adj[&u].iter().copied().filter(|&nb| nb != v).collect();
        let mut v_sub: Vec<usize> = adj[&v].iter().copied().filter(|&nb| nb != u).collect();
        let deg_u = u_sub.len();
        let deg_v = v_sub.len();
        if deg_u < 2 || deg_v < 2 {
            return (-10.0, 0.0, Vec::new());
        }
        if u_sub.len() > 16 {
            u_sub.sort_by(|&a, &b| {
                states[&(a, u)]
                    .1
                    .partial_cmp(&states[&(b, u)].1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            u_sub.truncate(12);
        }
        if v_sub.len() > 16 {
            v_sub.sort_by(|&a, &b| {
                states[&(a, v)]
                    .1
                    .partial_cmp(&states[&(b, v)].1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });
            v_sub.truncate(12);
        }
        let u_children: Vec<(Vec<u8>, f64)> =
            u_sub.iter().map(|&nb| states[&(nb, u)].clone()).collect();
        let v_children: Vec<(Vec<u8>, f64)> =
            v_sub.iter().map(|&nb| states[&(nb, v)].clone()).collect();
        let (gain_sub, m_opt, sites) = self.compute_split_mdl_gain_from_children(
            &u_children,
            &v_children,
            l,
            l_eff,
            ti_cost,
            tv_cost,
            site_weights,
        );
        let extra_topo = split_topology_cost_nats(deg_u, deg_v)
            - split_topology_cost_nats(u_children.len(), v_children.len());
        (gain_sub - extra_topo, m_opt, sites)
    }

    #[allow(clippy::too_many_arguments)]
    pub fn compute_split_mdl_gain_from_children(
        &self,
        u_children: &[(Vec<u8>, f64)],
        v_children: &[(Vec<u8>, f64)],
        l: usize,
        l_eff: f64,
        ti_cost: f64,
        tv_cost: f64,
        site_weights: &[f64],
    ) -> (f64, f64, Vec<usize>) {
        let deg_a = u_children.len();
        let deg_b = v_children.len();
        let topo_cost = split_topology_cost_nats(deg_a, deg_b);

        let all_children: Vec<&(Vec<u8>, f64)> =
            u_children.iter().chain(v_children.iter()).collect();
        let k_ch = all_children.len();
        let mut refined_loads = Vec::with_capacity(k_ch);

        let max_load = all_children.iter().map(|c| c.1).fold(0.0_f64, f64::max);
        if k_ch == 4 && max_load >= 40.0 {
            for idx_c in 0..4 {
                let others: Vec<&[u8]> = (0..4)
                    .filter(|&j| j != idx_c)
                    .map(|j| all_children[j].0.as_slice())
                    .collect();
                let mask_c = &all_children[idx_c].0;
                let mut n_anc = 0usize;
                let mut n_diff = 0usize;
                for s in 0..l {
                    let co = others[0][s] & others[1][s] & others[2][s];
                    let is_single = co == 1 || co == 2 || co == 4 || co == 8;
                    if is_single && mask_c[s] != 15 {
                        n_anc += 1;
                        if (mask_c[s] & co) == 0 {
                            n_diff += 1;
                        }
                    }
                }
                if n_anc >= 50 {
                    let p_c = ((n_diff as f64) / (n_anc as f64)).min(0.73);
                    let m_cond = -0.75 * (1.0 - (4.0 / 3.0) * p_c).ln() * l_eff * 1.05;
                    let m_prior = all_children[idx_c].1;
                    let r_load = if m_cond < 0.75 * m_prior {
                        m_prior.min(m_cond)
                    } else {
                        0.5 * (m_prior + m_cond)
                    };
                    refined_loads.push(r_load.max(1.5));
                } else {
                    refined_loads.push(all_children[idx_c].1);
                }
            }
        } else {
            for c in &all_children {
                refined_loads.push(c.1);
            }
        }

        let u_loads = &refined_loads[..deg_a];
        let v_loads = &refined_loads[deg_a..];

        let mut u_weights = Vec::with_capacity(deg_a);
        let mut v_weights = Vec::with_capacity(deg_b);
        for &m_c in u_loads {
            let mc = m_c.clamp(0.5, 0.48 * l_eff);
            u_weights.push(((l_eff - mc).max(1.5) / (mc + 1.0)).ln());
        }
        for &m_c in v_loads {
            let mc = m_c.clamp(0.5, 0.48 * l_eff);
            v_weights.push(((l_eff - mc).max(1.5) / (mc + 1.0)).ln());
        }

        let u_asym = if u_weights.len() >= 2 {
            u_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - u_weights.iter().copied().fold(f64::INFINITY, f64::min)
        } else {
            0.0
        };
        let v_asym = if v_weights.len() >= 2 {
            v_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - v_weights.iter().copied().fold(f64::INFINITY, f64::min)
        } else {
            0.0
        };

        let has_long = u_asym >= 1.5 || v_asym >= 1.5;
        let p_long_u = if u_asym >= 1.5 {
            (u_loads.iter().copied().fold(0.0_f64, f64::max) / l_eff).min(0.65)
        } else {
            0.0
        };
        let p_long_v = if v_asym >= 1.5 {
            (v_loads.iter().copied().fold(0.0_f64, f64::max) / l_eff).min(0.65)
        } else {
            0.0
        };
        let min_uw = u_weights.iter().copied().fold(f64::INFINITY, f64::min);
        let min_vw = v_weights.iter().copied().fold(f64::INFINITY, f64::min);

        let mut pos_sites: Vec<(f64, usize, f64, f64)> = Vec::new();
        let mut clean_sites: Vec<usize> = Vec::new();

        for s in 0..l {
            let mut u_valid = 0i16;
            let mut v_valid = 0i16;
            let mut u_miss = [0i16; 4];
            let mut v_miss = [0i16; 4];
            let mut u_shared = [0i16; 4];
            let mut v_shared = [0i16; 4];
            let mut cu = [0.0_f64; 4];
            let mut cv = [0.0_f64; 4];

            for ((s_mask, _), &w_c) in u_children.iter().zip(&u_weights) {
                let m = s_mask[s];
                if m != 15 {
                    u_valid += 1;
                    for b in 0..4 {
                        if ((m >> b) & 1) != 0 {
                            u_shared[b] += 1;
                        } else {
                            u_miss[b] += 1;
                            cu[b] += w_c;
                        }
                    }
                }
            }
            if u_valid < 2 {
                continue;
            }
            for ((s_mask, _), &w_c) in v_children.iter().zip(&v_weights) {
                let m = s_mask[s];
                if m != 15 {
                    v_valid += 1;
                    for b in 0..4 {
                        if ((m >> b) & 1) != 0 {
                            v_shared[b] += 1;
                        } else {
                            v_miss[b] += 1;
                            cv[b] += w_c;
                        }
                    }
                }
            }
            if v_valid < 2 {
                continue;
            }

            let mut steps_star = i16::MAX;
            let mut c_star = f64::INFINITY;
            let mut b_u_opt = 0usize;
            let mut b_v_opt = 0usize;
            for b in 0..4 {
                let st = u_miss[b] + v_miss[b];
                if st < steps_star {
                    steps_star = st;
                }
                let cs = cu[b] + cv[b];
                if cs < c_star {
                    c_star = cs;
                }
                if cu[b] < cu[b_u_opt] {
                    b_u_opt = b;
                }
                if cv[b] < cv[b_v_opt] {
                    b_v_opt = b;
                }
            }

            if steps_star < 2 || b_u_opt == b_v_opt {
                continue;
            }

            let min_cu = cu[b_u_opt];
            let min_cv = cv[b_v_opt];
            let steps_res = u_miss[b_u_opt] + v_miss[b_v_opt];

            let u_clean = u_shared[b_u_opt] >= 2 && u_miss[b_u_opt] == 0;
            let v_clean = v_shared[b_v_opt] >= 2 && v_miss[b_v_opt] == 0;
            let both_parsimony = (steps_res + 1) < steps_star;
            let raw_delta = c_star - (min_cu + min_cv);

            if both_parsimony && u_clean && v_clean {
                clean_sites.push(s);
            }

            let masked_u = v_clean && u_asym >= 1.5 && raw_delta >= 2.8;
            let masked_v = u_clean && v_asym >= 1.5 && raw_delta >= 2.8;
            let both_asym = u_asym >= 1.8 && v_asym >= 1.8 && raw_delta >= 3.0;

            if !(both_parsimony || masked_u || masked_v || both_asym) {
                continue;
            }

            let same_purine = IS_PURINE[b_u_opt] == IS_PURINE[b_v_opt];
            let subst_s = (if same_purine { ti_cost } else { tv_cost }) * site_weights[s];
            let delta_steps = if has_long {
                ((steps_star - steps_res) as f64).max(1.0)
            } else {
                1.0
            };

            let mut jc_mult_hit = 0.0_f64;
            if p_long_u > 0.05 && (masked_u || both_asym) {
                jc_mult_hit += 0.55 * p_long_u * (min_uw + subst_s);
            }
            if p_long_v > 0.05 && (masked_v || both_asym) {
                jc_mult_hit += 0.55 * p_long_v * (min_vw + subst_s);
            }

            let delta_pos = raw_delta * site_weights[s] + jc_mult_hit;
            let site_sav = delta_pos + delta_steps * subst_s;
            if site_sav > 0.0 {
                pos_sites.push((site_sav, s, subst_s, site_weights[s]));
            }
        }

        if pos_sites.is_empty() {
            return (-topo_cost, 0.0, Vec::new());
        }

        pos_sites.sort_by(|a, b| {
            b.0.partial_cmp(&a.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.1.cmp(&b.1))
        });

        let scale_pos_cost = if clean_sites.len() >= 7
            && (clean_sites.last().unwrap() - clean_sites[0]) >= 500
        {
            let uniq_codons: HashSet<usize> = clean_sites.iter().map(|&s| s / 3).collect();
            uniq_codons.len() >= 7
        } else {
            false
        };

        let mut cum_sav = 0.0_f64;
        let mut cum_subst = 0.0_f64;
        let mut cum_pos_cost = 0.0_f64;
        let mut prev_raw_pos = 0.0_f64;

        let mut best_data_gain = f64::NEG_INFINITY;
        let mut best_idx = 0usize;

        for (idx, &(sav, _, subst_s, w_s)) in pos_sites.iter().enumerate() {
            let m = (idx + 1) as f64;
            cum_sav += sav;
            cum_subst += subst_s;
            let header = 2.0_f64.ln() + (1.0 + m).ln() + 0.5 * (1.0 + (1.0 + m).ln()).ln();
            let raw_pos = log_binom(l_eff, m.min(0.49 * l_eff));
            if scale_pos_cost {
                let marg = raw_pos - prev_raw_pos;
                cum_pos_cost += marg * w_s.min(1.0);
                prev_raw_pos = raw_pos;
            } else {
                cum_pos_cost = raw_pos;
            }
            let branch_cost = header + cum_pos_cost + cum_subst;
            let net_gain = cum_sav - branch_cost;
            if net_gain > best_data_gain {
                best_data_gain = net_gain;
                best_idx = idx;
            }
        }

        let m_opt = (best_idx + 1) as f64;
        let chosen_sites: Vec<usize> = pos_sites[..=best_idx].iter().map(|x| x.1).collect();
        (best_data_gain - topo_cost, m_opt, chosen_sites)
    }

    pub fn reconstruct_nodal_states_and_mutations(
        &self,
        adj: &BTreeMap<usize, Vec<usize>>,
        bitmasks: &[u8],
        majority_state: &[i8],
        states: &HashMap<(usize, usize), (Vec<u8>, f64)>,
        n: usize,
        l: usize,
    ) -> (
        HashMap<usize, Vec<i8>>,
        HashMap<(usize, usize), Vec<(usize, i8, i8)>>,
    ) {
        let root = adj
            .iter()
            .find(|(&u, nbs)| u >= n && !nbs.is_empty())
            .or_else(|| adj.iter().find(|(_, nbs)| !nbs.is_empty()))
            .map(|(&u, _)| u)
            .unwrap_or(0);

        let mut parent: HashMap<usize, Option<usize>> = HashMap::new();
        parent.insert(root, None);
        let mut order = vec![root];
        let mut q = VecDeque::new();
        q.push_back(root);
        while let Some(curr) = q.pop_front() {
            if let Some(nbs) = adj.get(&curr) {
                for &nb in nbs {
                    if !parent.contains_key(&nb) {
                        parent.insert(nb, Some(curr));
                        order.push(nb);
                        q.push_back(nb);
                    }
                }
            }
        }

        let mut node_state: HashMap<usize, Vec<i8>> = HashMap::new();
        let mut node_mask_dict: HashMap<usize, Vec<u8>> = HashMap::new();
        let mut edge_mutations: HashMap<(usize, usize), Vec<(usize, i8, i8)>> = HashMap::new();

        let maj_bit: Vec<u8> = majority_state.iter().map(|&b| 1_u8 << (b as u8)).collect();

        for &u in &order {
            let p_opt = parent[&u];
            let c_mask = if u < n {
                bitmasks[u * l..(u + 1) * l].to_vec()
            } else {
                let children: Vec<(Vec<u8>, f64)> = adj[&u]
                    .iter()
                    .copied()
                    .filter(|&nb| Some(nb) != p_opt)
                    .map(|c| states[&(c, u)].clone())
                    .collect();
                Self::combine_fitch_children(&children).0
            };
            node_mask_dict.insert(u, c_mask.clone());

            if let Some(p) = p_opt {
                let p_st = &node_state[&p];
                let p_bit: Vec<u8> = p_st.iter().map(|&b| 1_u8 << (b as u8)).collect();
                let mpr_mask = if u < n {
                    c_mask
                } else {
                    let mut star_list: Vec<(Vec<u8>, f64)> = adj[&u]
                        .iter()
                        .copied()
                        .filter(|&nb| nb != p)
                        .map(|c| states[&(c, u)].clone())
                        .collect();
                    star_list.push((p_bit.clone(), 1.0));
                    Self::combine_fitch_children(&star_list).0
                };
                let mut st = vec![0_i8; l];
                for s in 0..l {
                    let m = mpr_mask[s];
                    if (m & p_bit[s]) != 0 {
                        st[s] = p_st[s];
                    } else if (m & maj_bit[s]) != 0 {
                        st[s] = majority_state[s];
                    } else {
                        for b in 0..4_i8 {
                            if ((m >> b) & 1) != 0 {
                                st[s] = b;
                                break;
                            }
                        }
                    }
                }
                node_state.insert(u, st);

                let st_p = &node_state[&p];
                let st_u = &node_state[&u];
                let mask_p = &node_mask_dict[&p];
                let mask_u = &node_mask_dict[&u];
                let mut muts = Vec::new();
                for s in 0..l {
                    if mask_p[s] != 15 && mask_u[s] != 15 && st_p[s] != st_u[s] {
                        if p < u {
                            muts.push((s, st_p[s], st_u[s]));
                        } else {
                            muts.push((s, st_u[s], st_p[s]));
                        }
                    }
                }
                edge_mutations.insert((p.min(u), p.max(u)), muts);
            } else {
                let mut st = vec![0_i8; l];
                for s in 0..l {
                    let m = c_mask[s];
                    if (m & maj_bit[s]) != 0 {
                        st[s] = majority_state[s];
                    } else {
                        for b in 0..4_i8 {
                            if ((m >> b) & 1) != 0 {
                                st[s] = b;
                                break;
                            }
                        }
                    }
                }
                node_state.insert(u, st);
            }
        }

        (node_state, edge_mutations)
    }

    #[allow(clippy::too_many_arguments)]
    fn peel_residual_shunts_mdl(
        &self,
        adj: &BTreeMap<usize, Vec<usize>>,
        node_labels: &BTreeMap<usize, String>,
        bitmasks: &[u8],
        majority_state: &[i8],
        terminal_loads: &[f64],
        edge_mutations: &HashMap<(usize, usize), Vec<(usize, i8, i8)>>,
        n: usize,
        l: usize,
        l_eff: f64,
        ti_cost: f64,
        tv_cost: f64,
        mat: &[i8],
        is_coding: bool,
    ) -> Vec<ShuntEdge> {
        let active_nodes: Vec<usize> = adj
            .iter()
            .filter(|(_, nbs)| !nbs.is_empty())
            .map(|(&u, _)| u)
            .collect();
        let v_count = active_nodes.len().max(4) as f64;
        let endpoint_cost = 2.0 * v_count.ln();

        let mut site_mut_counts = vec![0_i32; l];
        for mut_list in edge_mutations.values() {
            for &(s, _, _) in mut_list {
                site_mut_counts[s] += 1;
            }
        }
        let variable_sites = site_mut_counts.iter().filter(|&&c| c > 0).count().max(10);

        let root = adj
            .iter()
            .find(|(&u, nbs)| u >= n && !nbs.is_empty())
            .or_else(|| adj.iter().find(|(_, nbs)| !nbs.is_empty()))
            .map(|(&u, _)| u)
            .unwrap_or(0);

        let mut parent: HashMap<usize, Option<usize>> = HashMap::new();
        parent.insert(root, None);
        let mut q = VecDeque::new();
        q.push_back(root);
        while let Some(curr) = q.pop_front() {
            if let Some(nbs) = adj.get(&curr) {
                for &nb in nbs {
                    if !parent.contains_key(&nb) {
                        parent.insert(nb, Some(curr));
                        q.push_back(nb);
                    }
                }
            }
        }

        let mut branch_by_child: HashMap<usize, Vec<(usize, i8, i8)>> = HashMap::new();
        let mut branch_mut_load: HashMap<usize, f64> = HashMap::new();
        for (&u, &p_opt) in &parent {
            let Some(p) = p_opt else { continue };
            let key = (p.min(u), p.max(u));
            let raw_list = edge_mutations.get(&key).cloned().unwrap_or_default();
            let directed: Vec<(usize, i8, i8)> = if p < u {
                raw_list
            } else {
                raw_list.into_iter().map(|(s, bu, bp)| (s, bp, bu)).collect()
            };
            branch_mut_load.insert(u, (directed.len() as f64).max(0.5));
            branch_by_child.insert(u, directed);
        }

        let states_dir =
            self.compute_all_directed_fitch_states(adj, bitmasks, n, l, terminal_loads, l_eff);
        let (node_state, _) = self.reconstruct_nodal_states_and_mutations(
            adj,
            bitmasks,
            majority_state,
            &states_dir,
            n,
            l,
        );

        let mut node_valid: HashMap<usize, Vec<bool>> = HashMap::new();
        for &u in &active_nodes {
            if u < n {
                let row = &mat[u * l..(u + 1) * l];
                node_valid.insert(u, row.iter().map(|&x| x >= 0).collect());
            } else {
                let nb_check = parent
                    .get(&u)
                    .and_then(|&p| p)
                    .unwrap_or_else(|| adj[&u][0]);
                let m = &states_dir[&(u, nb_check)].0;
                node_valid.insert(u, m.iter().map(|&x| x != 15).collect());
            }
        }

        let orig_branch_by_child = branch_by_child.clone();
        let orig_branch_load = branch_mut_load.clone();
        let mut origin_node: HashMap<(usize, usize), usize> = HashMap::new();
        for (&u, muts) in &orig_branch_by_child {
            for &(s, _, _) in muts {
                origin_node.insert((u, s), u);
            }
        }

        let mut eff_anc_node: HashMap<usize, usize> = HashMap::new();
        for &u in &active_nodes {
            let Some(Some(p)) = parent.get(&u) else {
                continue;
            };
            let mut eff_anc = *p;
            if is_coding && u >= n {
                let mut cum_load = orig_branch_load.get(&u).copied().unwrap_or(1.0);
                let mut curr_w = *p;
                let mut steps_up = 0usize;
                let mut existing_sites: HashSet<usize> = branch_by_child
                    .get(&u)
                    .map(|v| v.iter().map(|x| x.0).collect())
                    .unwrap_or_default();
                let st_u_arr = &node_state[&u];
                let mut added_any = false;
                while curr_w >= n && steps_up < 3 {
                    let Some(Some(next_p)) = parent.get(&curr_w).copied() else {
                        break;
                    };
                    let add_load = orig_branch_load.get(&curr_w).copied().unwrap_or(1.0);
                    if cum_load + add_load > 48.0 {
                        break;
                    }
                    cum_load += add_load;
                    let st_anc_arr = &node_state[&next_p];
                    if let Some(w_muts) = orig_branch_by_child.get(&curr_w) {
                        for &(s_int, a_st, b_st) in w_muts {
                            if !existing_sites.contains(&s_int)
                                && node_valid[&u][s_int]
                                && st_u_arr[s_int] == b_st
                            {
                                let c_start = 3 * (s_int / 3);
                                if c_start + 2 < l {
                                    let aa_u = translate_bases(
                                        st_u_arr[c_start],
                                        st_u_arr[c_start + 1],
                                        st_u_arr[c_start + 2],
                                    );
                                    let aa_a = translate_bases(
                                        st_anc_arr[c_start],
                                        st_anc_arr[c_start + 1],
                                        st_anc_arr[c_start + 2],
                                    );
                                    if aa_u != aa_a {
                                        branch_by_child
                                            .entry(u)
                                            .or_default()
                                            .push((s_int, a_st, b_st));
                                        origin_node.insert((u, s_int), curr_w);
                                        existing_sites.insert(s_int);
                                        added_any = true;
                                    }
                                }
                            }
                        }
                    }
                    eff_anc = next_p;
                    curr_w = next_p;
                    steps_up += 1;
                }
                if added_any {
                    let entry = branch_mut_load.entry(u).or_insert(1.0);
                    *entry = entry.max(cum_load * 0.65);
                }
            }
            eff_anc_node.insert(u, eff_anc);
        }

        let mut site_to_branches: Vec<Vec<(usize, i8, i8)>> = vec![Vec::new(); l];
        for (&u, muts) in &branch_by_child {
            for &(s, a, b) in muts {
                site_to_branches[s].push((u, a, b));
            }
        }

        let is_ancestor = |a_node: usize, b_node: usize| -> bool {
            let mut curr = Some(b_node);
            while let Some(c) = curr {
                if c == a_node {
                    return true;
                }
                curr = parent.get(&c).copied().flatten();
            }
            false
        };

        let mut pair_shared_sites: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        let mut pair_priority_score: HashMap<(usize, usize), f64> = HashMap::new();

        for (s, entries) in site_to_branches.iter().enumerate() {
            if entries.len() < 2 || entries.len() > 36 {
                continue;
            }
            let bg_c = site_mut_counts[s].max(1) as f64;
            for i in 0..entries.len() {
                let (u, _, b_u) = entries[i];
                let w_u = origin_node.get(&(u, s)).copied().unwrap_or(u);
                for j in (i + 1)..entries.len() {
                    let (v, _, b_v) = entries[j];
                    if b_u == b_v {
                        let w_v = origin_node.get(&(v, s)).copied().unwrap_or(v);
                        if w_u == w_v {
                            continue;
                        }
                        if is_ancestor(u, v) || is_ancestor(v, u) {
                            continue;
                        }
                        if is_ancestor(w_u, v) || is_ancestor(w_v, u) {
                            continue;
                        }
                        if w_u != u || w_v != v {
                            let pw_u = parent.get(&w_u).copied().flatten();
                            let pw_v = parent.get(&w_v).copied().flatten();
                            let p_pw_u = pw_u.and_then(|x| parent.get(&x).copied().flatten());
                            let p_pw_v = pw_v.and_then(|x| parent.get(&x).copied().flatten());
                            if pw_u == pw_v || p_pw_u == pw_v || pw_u == p_pw_v {
                                continue;
                            }
                        }
                        let pair = (u.min(v), u.max(v));
                        pair_shared_sites.entry(pair).or_default().push(s);
                        let mut w_s = 1.0 / bg_c;
                        if is_coding {
                            let c_start = 3 * (s / 3);
                            if c_start + 2 < l {
                                let p_u_eff = parent
                                    .get(&w_u)
                                    .copied()
                                    .flatten()
                                    .or_else(|| eff_anc_node.get(&u).copied());
                                let p_v_eff = parent
                                    .get(&w_v)
                                    .copied()
                                    .flatten()
                                    .or_else(|| eff_anc_node.get(&v).copied());
                                if let (Some(pu), Some(pv)) = (p_u_eff, p_v_eff) {
                                    let st_u = &node_state[&u];
                                    let st_v = &node_state[&v];
                                    let st_pu = &node_state[&pu];
                                    let st_pv = &node_state[&pv];
                                    let aa_u = translate_bases(
                                        st_u[c_start],
                                        st_u[c_start + 1],
                                        st_u[c_start + 2],
                                    );
                                    let aa_v = translate_bases(
                                        st_v[c_start],
                                        st_v[c_start + 1],
                                        st_v[c_start + 2],
                                    );
                                    let aa_pu = translate_bases(
                                        st_pu[c_start],
                                        st_pu[c_start + 1],
                                        st_pu[c_start + 2],
                                    );
                                    let aa_pv = translate_bases(
                                        st_pv[c_start],
                                        st_pv[c_start + 1],
                                        st_pv[c_start + 2],
                                    );
                                    if aa_u == aa_v && aa_u != aa_pu && aa_v != aa_pv {
                                        w_s *= 4.0;
                                    }
                                }
                            }
                        }
                        *pair_priority_score.entry(pair).or_insert(0.0) += w_s;
                    }
                }
            }
        }

        let mut sorted_pairs: Vec<((usize, usize), Vec<usize>)> =
            pair_shared_sites.into_iter().collect();
        sorted_pairs.sort_by(|a, b| {
            let sa = pair_priority_score.get(&a.0).copied().unwrap_or(0.0);
            let sb = pair_priority_score.get(&b.0).copied().unwrap_or(0.0);
            sb.partial_cmp(&sa)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.1.len().cmp(&a.1.len()))
                .then_with(|| a.0.cmp(&b.0))
        });
        sorted_pairs.truncate(500);

        let mut candidate_shunts: Vec<ShuntEdge> = Vec::new();
        for ((u, v), mut shared_sites) in sorted_pairs {
            shared_sites.sort_unstable();
            shared_sites.dedup();
            let avail_sites: Vec<usize> = shared_sites
                .into_iter()
                .filter(|&s| node_valid[&u][s] && node_valid[&v][s])
                .collect();
            if avail_sites.len() < self.min_shunt_sites {
                continue;
            }

            let (best_int_gain, best_int_sites, best_span, int_data_sav, int_mask_cost) = self
                .evaluate_interval_shunt_mdl(
                    u,
                    v,
                    &avail_sites,
                    &site_mut_counts,
                    &node_state,
                    &node_valid,
                    &eff_anc_node,
                    &parent,
                    adj,
                    &branch_mut_load,
                    l,
                    l_eff,
                    variable_sites,
                    endpoint_cost,
                    ti_cost,
                    tv_cost,
                    is_coding,
                );

            let pu_direct = parent.get(&u).copied().flatten();
            let pv_direct = parent.get(&v).copied().flatten();
            let local_triad = pu_direct == pv_direct
                || (pu_direct.is_some()
                    && parent.get(&pu_direct.unwrap()).copied().flatten() == pv_direct)
                || (pv_direct.is_some()
                    && pu_direct == parent.get(&pv_direct.unwrap()).copied().flatten());

            let (best_sp_gain, best_sp_sites, sp_data_sav, sp_mask_cost) = if local_triad {
                (-10.0, Vec::new(), 0.0, 0.0)
            } else {
                self.evaluate_sparse_shunt_mdl(
                    u,
                    v,
                    &avail_sites,
                    &site_mut_counts,
                    &node_state,
                    &eff_anc_node,
                    &branch_mut_load,
                    l_eff,
                    variable_sites,
                    endpoint_cost,
                    ti_cost,
                    tv_cost,
                    is_coding,
                )
            };

            if best_int_gain <= 0.0 && best_sp_gain <= 0.0 {
                continue;
            }

            let (chosen_mode, chosen_sites, chosen_span, data_sav, mask_cost, net_gain) =
                if best_int_gain >= best_sp_gain && best_int_gain > 0.0 {
                    (
                        "INTERVAL".to_string(),
                        best_int_sites,
                        best_span,
                        int_data_sav,
                        int_mask_cost,
                        best_int_gain,
                    )
                } else {
                    (
                        "SPARSE".to_string(),
                        best_sp_sites,
                        None,
                        sp_data_sav,
                        sp_mask_cost,
                        best_sp_gain,
                    )
                };

            let post_prob = 1.0 / (1.0 + (-net_gain.clamp(-50.0, 50.0)).exp());
            candidate_shunts.push(ShuntEdge {
                node_u: u,
                node_v: v,
                label_u: node_labels
                    .get(&u)
                    .cloned()
                    .unwrap_or_else(|| format!("N{u}")),
                label_v: node_labels
                    .get(&v)
                    .cloned()
                    .unwrap_or_else(|| format!("N{v}")),
                mode: chosen_mode,
                shared_mutations: chosen_sites.len(),
                sites: chosen_sites,
                interval_span: chosen_span,
                data_savings_nats: data_sav,
                mask_cost_nats: mask_cost,
                endpoint_cost_nats: endpoint_cost,
                net_mdl_gain_nats: net_gain,
                posterior_prob: post_prob,
            });
        }

        candidate_shunts.sort_by(|a, b| {
            b.net_mdl_gain_nats
                .partial_cmp(&a.net_mdl_gain_nats)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let mut certified_shunts: Vec<ShuntEdge> = Vec::new();
        let mut used_sites_per_node: HashMap<usize, HashSet<usize>> = HashMap::new();

        for mut sh in candidate_shunts {
            let (u, v) = (sh.node_u, sh.node_v);
            let overlap: HashSet<usize> = sh
                .sites
                .iter()
                .copied()
                .filter(|&s| {
                    let wu = origin_node.get(&(u, s)).copied().unwrap_or(u);
                    let wv = origin_node.get(&(v, s)).copied().unwrap_or(v);
                    used_sites_per_node
                        .get(&u)
                        .map(|set| set.contains(&s))
                        .unwrap_or(false)
                        || used_sites_per_node
                            .get(&v)
                            .map(|set| set.contains(&s))
                            .unwrap_or(false)
                        || used_sites_per_node
                            .get(&wu)
                            .map(|set| set.contains(&s))
                            .unwrap_or(false)
                        || used_sites_per_node
                            .get(&wv)
                            .map(|set| set.contains(&s))
                            .unwrap_or(false)
                })
                .collect();

            if sh.sites.len().saturating_sub(overlap.len()) < self.min_shunt_sites {
                continue;
            }
            if !overlap.is_empty() {
                let rem_sites: Vec<usize> = sh
                    .sites
                    .iter()
                    .copied()
                    .filter(|s| !overlap.contains(s))
                    .collect();
                let (g_r, sites_r, span_r, sav_r, mc_r) = if sh.mode == "INTERVAL" {
                    self.evaluate_interval_shunt_mdl(
                        u,
                        v,
                        &rem_sites,
                        &site_mut_counts,
                        &node_state,
                        &node_valid,
                        &eff_anc_node,
                        &parent,
                        adj,
                        &branch_mut_load,
                        l,
                        l_eff,
                        variable_sites,
                        endpoint_cost,
                        ti_cost,
                        tv_cost,
                        is_coding,
                    )
                } else {
                    let (g, st, sv, mc) = self.evaluate_sparse_shunt_mdl(
                        u,
                        v,
                        &rem_sites,
                        &site_mut_counts,
                        &node_state,
                        &eff_anc_node,
                        &branch_mut_load,
                        l_eff,
                        variable_sites,
                        endpoint_cost,
                        ti_cost,
                        tv_cost,
                        is_coding,
                    );
                    (g, st, None, sv, mc)
                };
                if g_r <= 0.0 {
                    continue;
                }
                sh.shared_mutations = sites_r.len();
                sh.sites = sites_r;
                sh.interval_span = span_r;
                sh.data_savings_nats = sav_r;
                sh.mask_cost_nats = mc_r;
                sh.net_mdl_gain_nats = g_r;
                sh.posterior_prob = 1.0 / (1.0 + (-g_r.clamp(-50.0, 50.0)).exp());
            }

            for &s in &sh.sites {
                let wu = origin_node.get(&(u, s)).copied().unwrap_or(u);
                let wv = origin_node.get(&(v, s)).copied().unwrap_or(v);
                used_sites_per_node.entry(u).or_default().insert(s);
                used_sites_per_node.entry(v).or_default().insert(s);
                used_sites_per_node.entry(wu).or_default().insert(s);
                used_sites_per_node.entry(wv).or_default().insert(s);
            }
            certified_shunts.push(sh);
        }

        certified_shunts.sort_by(|a, b| {
            b.net_mdl_gain_nats
                .partial_cmp(&a.net_mdl_gain_nats)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        certified_shunts
    }

    #[allow(clippy::too_many_arguments)]
    fn compute_parallel_site_savings(
        &self,
        u: usize,
        v: usize,
        shared_sites: &[usize],
        site_mut_counts: &[i32],
        node_state: &HashMap<usize, Vec<i8>>,
        eff_anc_node: &HashMap<usize, usize>,
        branch_mut_load: &HashMap<usize, f64>,
        l_eff: f64,
        variable_sites: usize,
        ti_cost: f64,
        tv_cost: f64,
        is_coding: bool,
    ) -> Vec<(usize, f64)> {
        let (Some(&p_u), Some(&p_v)) = (eff_anc_node.get(&u), eff_anc_node.get(&v)) else {
            return Vec::new();
        };
        let st_u = &node_state[&u];
        let st_pu = &node_state[&p_u];
        let st_v = &node_state[&v];
        let st_pv = &node_state[&p_v];
        let l = st_u.len();
        let m_u = branch_mut_load.get(&u).copied().unwrap_or(1.0);
        let m_v = branch_mut_load.get(&v).copied().unwrap_or(1.0);

        let mut results = Vec::new();
        for &s in shared_sites {
            if st_u[s] != st_pu[s] && st_v[s] != st_pv[s] && st_u[s] == st_v[s] {
                let bg_muts = (site_mut_counts[s] - 2).max(0) as f64;
                let b_anc = st_pu[s] as usize;
                let b_der = st_u[s] as usize;
                let subst_s = if IS_PURINE[b_anc] == IS_PURINE[b_der] {
                    ti_cost
                } else {
                    tv_cost
                };

                let mut codon_bonus = 0.0_f64;
                let mut is_synonymous = false;
                if is_coding {
                    let c_start = 3 * (s / 3);
                    if c_start + 2 < l {
                        let aa_u =
                            translate_bases(st_u[c_start], st_u[c_start + 1], st_u[c_start + 2]);
                        let aa_pu =
                            translate_bases(st_pu[c_start], st_pu[c_start + 1], st_pu[c_start + 2]);
                        let aa_v =
                            translate_bases(st_v[c_start], st_v[c_start + 1], st_v[c_start + 2]);
                        let aa_pv =
                            translate_bases(st_pv[c_start], st_pv[c_start + 1], st_pv[c_start + 2]);
                        if aa_u != aa_pu && aa_v != aa_pv && aa_u == aa_v {
                            codon_bonus = 19.0_f64.ln() - subst_s;
                        } else if aa_u == aa_pu && aa_v == aa_pv {
                            is_synonymous = true;
                        }
                    }
                }

                let expected_collision = (m_u * m_v) / (3.0 * l_eff).max(1.0);
                let (bp_term, conservation_surprise) = if is_synonymous {
                    let bp = ((variable_sites as f64) / (0.5 * (m_u + m_v)).max(0.5))
                        .max(1.0)
                        .ln()
                        - expected_collision.max(1.0).ln().max(0.0);
                    (bp, 0.0)
                } else {
                    let bp = (l_eff / (0.5 * (m_u + m_v)).max(0.5)).max(1.0).ln()
                        - expected_collision.max(1.0).ln().max(0.0);
                    let eff_var_pool = ((variable_sites as f64)
                        * if codon_bonus > 0.0 { 0.45 } else { 1.0 })
                    .max(1.0);
                    let cs = ((l_eff / eff_var_pool) * (3.0 / (bg_muts + 0.35)))
                        .ln()
                        .max(0.0);
                    (bp, cs)
                };

                let sav = bp_term + subst_s + codon_bonus + conservation_surprise;
                if sav > 0.0 {
                    results.push((s, sav));
                }
            }
        }
        results
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_interval_shunt_mdl(
        &self,
        u: usize,
        v: usize,
        shared_sites: &[usize],
        site_mut_counts: &[i32],
        node_state: &HashMap<usize, Vec<i8>>,
        node_valid: &HashMap<usize, Vec<bool>>,
        eff_anc_node: &HashMap<usize, usize>,
        parent: &HashMap<usize, Option<usize>>,
        adj: &BTreeMap<usize, Vec<usize>>,
        branch_mut_load: &HashMap<usize, f64>,
        l: usize,
        l_eff: f64,
        variable_sites: usize,
        endpoint_cost: f64,
        ti_cost: f64,
        tv_cost: f64,
        is_coding: bool,
    ) -> (f64, Vec<usize>, Option<(usize, usize)>, f64, f64) {
        let site_sav_pairs = self.compute_parallel_site_savings(
            u,
            v,
            shared_sites,
            site_mut_counts,
            node_state,
            eff_anc_node,
            branch_mut_load,
            l_eff,
            variable_sites,
            ti_cost,
            tv_cost,
            is_coding,
        );
        let k = site_sav_pairs.len();
        let min_k = self.min_shunt_sites.max(3);
        if k < min_k {
            return (-1e9, Vec::new(), None, 0.0, 0.0);
        }

        let cand_sites: Vec<usize> = site_sav_pairs.iter().map(|x| x.0).collect();
        let cand_savs: Vec<f64> = site_sav_pairs.iter().map(|x| x.1).collect();
        let base_mask_cost = 2.0_f64.ln() + 2.0 * (l.max(2) as f64).ln();

        let val_u = &node_valid[&u];
        let val_v = &node_valid[&v];
        let both_val: Vec<bool> = (0..l).map(|s| val_u[s] && val_v[s]).collect();

        let mut u_reps = vec![u];
        if let Some(nbs) = adj.get(&u) {
            for &nb in nbs {
                if parent.get(&nb).copied().flatten() == Some(u) {
                    u_reps.push(nb);
                }
            }
        }
        let mut v_reps = vec![v];
        if let Some(nbs) = adj.get(&v) {
            for &nb in nbs {
                if parent.get(&nb).copied().flatten() == Some(v) {
                    v_reps.push(nb);
                }
            }
        }

        let mut disc_prefs: Vec<Vec<i32>> = Vec::new();
        let mut disc_pos_list: Vec<Vec<usize>> = Vec::new();
        for &ur in &u_reps {
            let st_ur = &node_state[&ur];
            for &vr in &v_reps {
                let st_vr = &node_state[&vr];
                let mut pref = Vec::with_capacity(l + 1);
                pref.push(0);
                let mut pos = Vec::new();
                let mut acc = 0i32;
                for s in 0..l {
                    if both_val[s] && st_ur[s] != st_vr[s] {
                        acc += 1;
                        pos.push(s);
                    }
                    pref.push(acc);
                }
                disc_prefs.push(pref);
                disc_pos_list.push(pos);
            }
        }

        let mut gap_pref = Vec::with_capacity(l + 1);
        gap_pref.push(0i32);
        let mut g_acc = 0i32;
        for &bv in &both_val {
            if !bv {
                g_acc += 1;
            }
            gap_pref.push(g_acc);
        }

        let mut sav_pref = Vec::with_capacity(k + 1);
        sav_pref.push(0.0_f64);
        let mut s_acc = 0.0_f64;
        for &sv in &cand_savs {
            s_acc += sv;
            sav_pref.push(s_acc);
        }

        let mut best_gain = -1e9_f64;
        let mut best_sites = Vec::new();
        let mut best_span = None;
        let mut best_data_sav = 0.0_f64;
        let mut best_mask_cost = base_mask_cost;

        for i in 0..k {
            let a = cand_sites[i];
            for j in (i + min_k - 1)..k {
                let b = cand_sites[j];
                let span_len = b - a + 1;
                let m_in = j - i + 1;
                if span_len < 15 {
                    continue;
                }
                let gaps_in = (gap_pref[b + 1] - gap_pref[a]) as usize;
                if (gaps_in as f64) > 0.12 * (span_len as f64) {
                    continue;
                }
                let mut best_rep_idx = 0usize;
                let mut discordant = i32::MAX;
                for (r_idx, dp) in disc_prefs.iter().enumerate() {
                    let d_val = dp[b + 1] - dp[a];
                    if d_val < discordant {
                        discordant = d_val;
                        best_rep_idx = r_idx;
                    }
                }
                let disc_f = discordant as f64;
                if disc_f > 0.25 * (m_in as f64) || (disc_f / (span_len as f64)) > 0.03 {
                    continue;
                }
                let disc_cost = if discordant > 0 {
                    branch_mutation_code_nats(disc_f, span_len as f64, 0.0, 4)
                } else {
                    0.0
                };
                let total_mask_cost = base_mask_cost + disc_cost;
                let data_sav = sav_pref[j + 1] - sav_pref[i];
                let net_gain = data_sav - endpoint_cost - total_mask_cost;
                if net_gain > best_gain {
                    best_gain = net_gain;
                    best_sites = cand_sites[i..=j].to_vec();
                    let disc_positions = &disc_pos_list[best_rep_idx];
                    let s_left = disc_positions
                        .iter()
                        .copied()
                        .rfind(|&p| p < a)
                        .unwrap_or_else(|| a.saturating_sub(20));
                    let s_right = disc_positions
                        .iter()
                        .copied()
                        .find(|&p| p > b)
                        .unwrap_or_else(|| (b + 20).min(l - 1));
                    best_span = Some(((s_left + a) / 2, (b + s_right) / 2));
                    best_data_sav = data_sav;
                    best_mask_cost = total_mask_cost;
                }
            }
        }

        (
            best_gain,
            best_sites,
            best_span,
            best_data_sav,
            best_mask_cost,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_sparse_shunt_mdl(
        &self,
        u: usize,
        v: usize,
        shared_sites: &[usize],
        site_mut_counts: &[i32],
        node_state: &HashMap<usize, Vec<i8>>,
        eff_anc_node: &HashMap<usize, usize>,
        branch_mut_load: &HashMap<usize, f64>,
        l_eff: f64,
        variable_sites: usize,
        endpoint_cost: f64,
        ti_cost: f64,
        tv_cost: f64,
        is_coding: bool,
    ) -> (f64, Vec<usize>, f64, f64) {
        let site_sav_pairs = self.compute_parallel_site_savings(
            u,
            v,
            shared_sites,
            site_mut_counts,
            node_state,
            eff_anc_node,
            branch_mut_load,
            l_eff,
            variable_sites,
            ti_cost,
            tv_cost,
            is_coding,
        );
        if site_sav_pairs.len() < self.min_shunt_sites {
            return (-1e9, Vec::new(), 0.0, 0.0);
        }

        let mut strict_parallel: Vec<(f64, usize)> =
            site_sav_pairs.into_iter().map(|(s, sav)| (sav, s)).collect();
        strict_parallel.sort_by(|a, b| {
            b.0.partial_cmp(&a.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| b.1.cmp(&a.1))
        });

        let mut m_pool = (variable_sites as f64)
            .max((strict_parallel.len() as f64) * 2.0)
            .max(20.0);
        if is_coding {
            m_pool = ((variable_sites as f64) * 0.35)
                .min(l_eff / 3.0)
                .max(20.0);
        }

        let mut best_gain = -1e9_f64;
        let mut best_sites = Vec::new();
        let mut best_data_sav = 0.0_f64;
        let mut best_mask_cost = 0.0_f64;
        let mut cum_sav = 0.0_f64;
        let mut cur_sites = Vec::new();

        for (idx, &(sav, s)) in strict_parallel.iter().enumerate() {
            let m_idx = idx + 1;
            cum_sav += sav;
            cur_sites.push(s);
            if m_idx >= self.min_shunt_sites {
                let mask_cost = 2.0_f64.ln() + log_binom(m_pool, m_idx as f64);
                let net_gain = cum_sav - endpoint_cost - mask_cost;
                if net_gain > best_gain {
                    best_gain = net_gain;
                    let mut sorted_s = cur_sites.clone();
                    sorted_s.sort_unstable();
                    best_sites = sorted_s;
                    best_data_sav = cum_sav;
                    best_mask_cost = mask_cost;
                }
            }
        }

        (best_gain, best_sites, best_data_sav, best_mask_cost)
    }

    fn compute_star_tree_mdl_nats(
        &self,
        bitmasks: &[u8],
        n: usize,
        l: usize,
        l_eff: f64,
    ) -> f64 {
        let child_list: Vec<(Vec<u8>, f64)> = (0..n)
            .map(|i| (bitmasks[i * l..(i + 1) * l].to_vec(), 1.0))
            .collect();
        let (root_mask, _) = Self::combine_fitch_children(&child_list);
        let mut total_nats = l_eff * 4.0_f64.ln();
        for i in 0..n {
            let row = &bitmasks[i * l..(i + 1) * l];
            let mut muts = 0usize;
            for s in 0..l {
                if row[s] != 15 && (row[s] & root_mask[s]) == 0 {
                    muts += 1;
                }
            }
            total_nats += branch_mutation_code_nats(muts as f64, l_eff, 0.0, 4);
        }
        total_nats
    }

    fn compute_tree_total_mdl_nats(
        &self,
        adj: &BTreeMap<usize, Vec<usize>>,
        branches: &[BranchRecord],
        l_eff: f64,
    ) -> f64 {
        let mut total_nats = l_eff * 4.0_f64.ln();
        for b in branches {
            total_nats += branch_mutation_code_nats(b.mutation_count, l_eff, 0.0, 4);
            if b.is_internal {
                let deg_u = adj[&b.u].len().saturating_sub(1);
                let deg_v = adj[&b.v].len().saturating_sub(1);
                total_nats += split_topology_cost_nats(deg_u, deg_v);
            }
        }
        total_nats
    }

    fn serialize_newick(
        &self,
        adj: &BTreeMap<usize, Vec<usize>>,
        names: &[String],
        branches: &[BranchRecord],
        n: usize,
    ) -> String {
        if n == 1 {
            return format!("{}:0.0;", names[0]);
        }
        if n == 2 {
            let b_len = branches.first().map(|b| b.branch_length).unwrap_or(0.0);
            return format!("({}:{:.6},{}:{:.6});", names[0], b_len / 2.0, names[1], b_len / 2.0);
        }

        let mut edge_map: HashMap<(usize, usize), &BranchRecord> = HashMap::new();
        for b in branches {
            edge_map.insert((b.u.min(b.v), b.u.max(b.v)), b);
        }

        let root = adj
            .iter()
            .filter(|(&u, nbs)| u >= n && !nbs.is_empty())
            .max_by_key(|(_, nbs)| nbs.len())
            .map(|(&u, _)| u)
            .unwrap_or(n);

        let mut visited: HashSet<usize> = HashSet::new();

        fn build_sub(
            u: usize,
            p: Option<usize>,
            adj: &BTreeMap<usize, Vec<usize>>,
            names: &[String],
            edge_map: &HashMap<(usize, usize), &BranchRecord>,
            visited: &mut HashSet<usize>,
            n: usize,
        ) -> String {
            visited.insert(u);
            let label = if u < n {
                names[u].clone()
            } else {
                let children: Vec<usize> = adj[&u]
                    .iter()
                    .copied()
                    .filter(|&nb| Some(nb) != p && !visited.contains(&nb))
                    .collect();
                let sub_strs: Vec<String> = children
                    .into_iter()
                    .map(|c| build_sub(c, Some(u), adj, names, edge_map, visited, n))
                    .collect();
                format!("({})", sub_strs.join(","))
            };
            if let Some(parent_node) = p {
                let rec = edge_map.get(&(u.min(parent_node), u.max(parent_node)));
                let bl = rec.map(|r| r.branch_length).unwrap_or(0.0);
                if u >= n {
                    if let Some(r) = rec {
                        return format!("{label}{:.3}:{bl:.6}", r.mdl_support);
                    }
                }
                format!("{label}:{bl:.6}")
            } else {
                format!("{label};")
            }
        }

        build_sub(root, None, adj, names, &edge_map, &mut visited, n)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_star_radiation_contracts_to_zero_splits() {
        // 8 taxa radiating from an ancestral AAA... sequence with independent terminal mutations
        let l = 300;
        let names: Vec<String> = (0..8).map(|i| format!("T{i}")).collect();
        let mut seqs = Vec::new();
        for i in 0..8 {
            let mut s = vec![b'A'; l];
            for k in 0..5 {
                s[i * 10 + k] = b'G';
            }
            seqs.push(String::from_utf8(s).unwrap());
        }
        let engine = BushwhackerEngine::default();
        let res = engine.infer(&names, &seqs).unwrap();
        assert_eq!(res.n_internal_branches, 0);
        assert!(res.polytomy_nodes >= 1);
    }

    #[test]
    fn test_clean_clade_resolves_with_high_support() {
        let l = 300;
        let names: Vec<String> = vec!["A1", "A2", "B1", "B2", "B3"]
            .into_iter()
            .map(String::from)
            .collect();
        let mut s_a1 = vec![b'A'; l];
        let mut s_a2 = vec![b'A'; l];
        let mut s_b1 = vec![b'A'; l];
        let mut s_b2 = vec![b'A'; l];
        let mut s_b3 = vec![b'A'; l];
        // 5 clean synapomorphies uniting (A1, A2) vs (B1, B2, B3)
        for pos in [11, 41, 89, 149, 221] {
            s_a1[pos] = b'G';
            s_a2[pos] = b'G';
        }
        // Terminal mutations
        s_a1[1] = b'C';
        s_a2[2] = b'T';
        s_b1[3] = b'C';
        s_b2[4] = b'T';
        s_b3[5] = b'G';
        let seqs = vec![
            String::from_utf8(s_a1).unwrap(),
            String::from_utf8(s_a2).unwrap(),
            String::from_utf8(s_b1).unwrap(),
            String::from_utf8(s_b2).unwrap(),
            String::from_utf8(s_b3).unwrap(),
        ];
        let engine = BushwhackerEngine::default();
        let res = engine.infer(&names, &seqs).unwrap();
        assert_eq!(res.n_internal_branches, 1);
        assert!(res.backbone_tree_nats < res.star_tree_nats);
    }
}

