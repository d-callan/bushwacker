use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::env;
use std::fs;
use std::path::Path;
use std::process;

use bushwhacker::{
    compress_alignment_to_bwz, decompress_alignment_from_bwz, parse_alignment_file,
    BushwhackerEngine, BushwhackerResult,
};

fn print_usage(prog: &str) {
    eprintln!(
        "BUSHWHACKER: First-Principles Kolmogorov / MDL Phylogenetic Graph Induction Engine\n\
         (Bipartition Unification, Shunt Hunting, & Weighted Hub Aggregation via\n\
          Combinadic Kolmogorov Entropy Reduction)\n\n\
         Usage:\n  \
           {prog} <input.fasta|phylip|nexus|bwz> [options]\n\n\
         Supported Input Formats (auto-detected):\n  \
           - FASTA   (.fa, .fasta, .fna)\n  \
           - PHYLIP  (.phy, .phylip; strict 10-char or relaxed, sequential or interleaved)\n  \
           - NEXUS   (.nex, .nexus, .nxs; DATA or CHARACTERS blocks, sequential or interleaved)\n  \
           - BWZ     (.bwz; Bushwhacker lossless phylogenetic archive for --decompress)\n\n\
         Options:\n  \
           -o, --out-newick <PATH>      Output Newick file path (with MDL branch supports & lengths)\n  \
           -j, --out-json <PATH>        Output JSON report path (full MDL ledger, splits, & shunts)\n  \
           -c, --compress <PATH>        Write lossless .bwz phylogenetic archive to <PATH>\n  \
           -d, --decompress <PATH>      Decompress input .bwz archive into FASTA at <PATH>\n  \
           --no-coding                  Disable automatic coding-frame 3rd-position entropy weighting\n  \
           --no-shunts                  Disable Stage 4 reticulate shunt peeling\n  \
           --min-shunt-sites <N>        Minimum parallel mutations to evaluate a shunt (default: 2)\n  \
           -h, --help                   Print this help message"
    );
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn clade_leaves_for_shunt_endpoint(
    adj: &BTreeMap<usize, Vec<usize>>,
    names: &[String],
    n: usize,
    u: usize,
    other: usize,
) -> Vec<String> {
    if u < n {
        return vec![names[u].clone()];
    }
    // Find neighbor of `u` on the path toward `other`
    let mut parent: HashMap<usize, usize> = HashMap::new();
    let mut q = VecDeque::from([u]);
    let mut visited: HashSet<usize> = HashSet::from([u]);
    while let Some(curr) = q.pop_front() {
        if curr == other {
            break;
        }
        if let Some(nbs) = adj.get(&curr) {
            for &nb in nbs {
                if visited.insert(nb) {
                    parent.insert(nb, curr);
                    q.push_back(nb);
                }
            }
        }
    }
    let mut step_to_other = None;
    let mut cur = other;
    while let Some(&p) = parent.get(&cur) {
        if p == u {
            step_to_other = Some(cur);
            break;
        }
        cur = p;
    }

    // Collect leaves in each branch attached to `u` except `step_to_other`
    let mut branch_leaf_sets: Vec<Vec<String>> = Vec::new();
    if let Some(nbs) = adj.get(&u) {
        for &nb in nbs {
            if Some(nb) == step_to_other {
                continue;
            }
            let mut leaves = Vec::new();
            let mut bq = VecDeque::from([nb]);
            let mut bvis = HashSet::from([u, nb]);
            while let Some(x) = bq.pop_front() {
                if x < n {
                    leaves.push(names[x].clone());
                } else if let Some(xnbs) = adj.get(&x) {
                    for &y in xnbs {
                        if bvis.insert(y) {
                            bq.push_back(y);
                        }
                    }
                }
            }
            if !leaves.is_empty() {
                branch_leaf_sets.push(leaves);
            }
        }
    }
    let total_away: usize = branch_leaf_sets.iter().map(|v| v.len()).sum();
    if total_away > 0 && total_away <= n / 2 {
        let mut all_away: Vec<String> = branch_leaf_sets.into_iter().flatten().collect();
        all_away.sort();
        all_away.truncate(16);
        return all_away;
    }
    if let Some(mut smallest) = branch_leaf_sets.into_iter().min_by_key(|v| v.len()) {
        smallest.sort();
        smallest.truncate(16);
        return smallest;
    }
    Vec::new()
}

fn serialize_json(res: &BushwhackerResult, fasta_path: &str) -> String {
    let internal_branches: Vec<_> = res
        .branches
        .iter()
        .filter(|b| b.is_internal)
        .map(|b| {
            format!(
                "    {{\"u\": {}, \"v\": {}, \"mutations\": {:.2}, \"branch_length\": {:.6}, \"mdl_gain_nats\": {:.3}, \"mdl_support\": {:.4}, \"num_synapomorphies\": {}}}",
                b.u,
                b.v,
                b.mutation_count,
                b.branch_length,
                b.mdl_gain_nats,
                b.mdl_support,
                b.synapomorphy_sites.len()
            )
        })
        .collect();

    let shunt_entries: Vec<_> = res
        .shunts
        .iter()
        .map(|s| {
            let sites_str = s
                .sites
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(", ");
            let codons_str = s
                .sites
                .iter()
                .map(|x| (x / 3 + 1).to_string())
                .collect::<Vec<_>>()
                .join(", ");
            let span_str = match s.interval_span {
                Some((a, b)) => format!("[{a}, {b}]"),
                None => "null".to_string(),
            };
            let leaves_u = clade_leaves_for_shunt_endpoint(
                &res.adjacency,
                &res.names,
                res.n_taxa,
                s.node_u,
                s.node_v,
            )
            .into_iter()
            .map(|x| json_escape(&x))
            .collect::<Vec<_>>()
            .join(", ");
            let leaves_v = clade_leaves_for_shunt_endpoint(
                &res.adjacency,
                &res.names,
                res.n_taxa,
                s.node_v,
                s.node_u,
            )
            .into_iter()
            .map(|x| json_escape(&x))
            .collect::<Vec<_>>()
            .join(", ");
            format!(
                "    {{\n      \"node_u\": {},\n      \"node_v\": {},\n      \"label_u\": {},\n      \"label_v\": {},\n      \"leaves_u\": [{}],\n      \"leaves_v\": [{}],\n      \"mode\": {},\n      \"shared_mutations\": {},\n      \"sites_0based\": [{}],\n      \"codons_1based\": [{}],\n      \"interval_span\": {},\n      \"data_savings_nats\": {:.3},\n      \"mask_cost_nats\": {:.3},\n      \"endpoint_cost_nats\": {:.3},\n      \"net_mdl_gain_nats\": {:.3},\n      \"posterior_prob\": {:.4}\n    }}",
                s.node_u,
                s.node_v,
                json_escape(&s.label_u),
                json_escape(&s.label_v),
                leaves_u,
                leaves_v,
                json_escape(&s.mode),
                s.shared_mutations,
                sites_str,
                codons_str,
                span_str,
                s.data_savings_nats,
                s.mask_cost_nats,
                s.endpoint_cost_nats,
                s.net_mdl_gain_nats,
                s.posterior_prob
            )
        })
        .collect();

    format!(
        "{{\n  \"tool\": \"BUSHWHACKER\",\n  \"input\": {},\n  \"n_taxa\": {},\n  \"seq_len\": {},\n  \"raw_Description_Length_nats\": {:.2},\n  \"star_tree_MDL_nats\": {:.2},\n  \"tree_MDL_nats\": {:.2},\n  \"network_MDL_nats\": {:.2},\n  \"compression_ratio_vs_star\": {:.4},\n  \"resolved_internal_splits\": {},\n  \"max_binary_splits\": {},\n  \"polytomy_nodes\": {},\n  \"certified_shunts\": {},\n  \"runtime_sec\": {:.4},\n  \"newick\": {},\n  \"binary_newick\": {},\n  \"internal_branches\": [\n{}\n  ],\n  \"shunts\": [\n{}\n  ]\n}}\n",
        json_escape(fasta_path),
        res.n_taxa,
        res.n_sites,
        res.raw_uncompressed_nats,
        res.star_tree_nats,
        res.backbone_tree_nats,
        res.full_graph_nats,
        res.full_graph_nats / res.star_tree_nats.max(1e-9),
        res.n_internal_branches,
        res.max_binary_internal_branches,
        res.polytomy_nodes,
        res.shunts.len(),
        res.runtime_sec,
        json_escape(&res.newick),
        json_escape(&res.binary_newick),
        internal_branches.join(",\n"),
        shunt_entries.join(",\n")
    )
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let prog = args.first().map(|s| s.as_str()).unwrap_or("bushwhacker");

    if args.len() < 2 {
        print_usage(prog);
        process::exit(1);
    }

    let mut fasta_path: Option<String> = None;
    let mut out_newick: Option<String> = None;
    let mut out_json: Option<String> = None;
    let mut out_compress: Option<String> = None;
    let mut out_decompress: Option<String> = None;
    let mut no_coding = false;
    let mut no_shunts = false;
    let mut min_shunt_sites: usize = 2;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                print_usage(prog);
                process::exit(0);
            }
            "-o" | "--out-newick" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: missing argument for {}", args[i - 1]);
                    process::exit(1);
                }
                out_newick = Some(args[i].clone());
            }
            "-j" | "--out-json" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: missing argument for {}", args[i - 1]);
                    process::exit(1);
                }
                out_json = Some(args[i].clone());
            }
            "-c" | "--compress" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: missing argument for {}", args[i - 1]);
                    process::exit(1);
                }
                out_compress = Some(args[i].clone());
            }
            "-d" | "--decompress" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: missing argument for {}", args[i - 1]);
                    process::exit(1);
                }
                out_decompress = Some(args[i].clone());
            }
            "--no-coding" => {
                no_coding = true;
            }
            "--no-shunts" => {
                no_shunts = true;
            }
            "--min-shunt-sites" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("Error: missing argument for --min-shunt-sites");
                    process::exit(1);
                }
                min_shunt_sites = args[i].parse().unwrap_or_else(|_| {
                    eprintln!("Error: invalid integer for --min-shunt-sites: {}", args[i]);
                    process::exit(1);
                });
            }
            arg if arg.starts_with('-') => {
                eprintln!("Error: unknown option '{arg}'");
                print_usage(prog);
                process::exit(1);
            }
            arg => {
                if fasta_path.is_none() {
                    fasta_path = Some(arg.to_string());
                } else {
                    eprintln!("Error: unexpected positional argument '{arg}'");
                    process::exit(1);
                }
            }
        }
        i += 1;
    }

    let Some(fasta) = fasta_path else {
        print_usage(prog);
        process::exit(1);
    };

    if let Some(dec_fasta_out) = out_decompress {
        let bwz_bytes = fs::read(&fasta).unwrap_or_else(|e| {
            eprintln!("Error reading .bwz archive '{fasta}': {e}");
            process::exit(1);
        });
        let (names, seqs, newick) = decompress_alignment_from_bwz(&bwz_bytes).unwrap_or_else(|e| {
            eprintln!("Error decompressing '{fasta}': {e}");
            process::exit(1);
        });
        let mut out_str = String::with_capacity(names.len() * (seqs[0].len() + 32));
        for (nm, sq) in names.iter().zip(&seqs) {
            out_str.push('>');
            out_str.push_str(nm);
            out_str.push('\n');
            out_str.push_str(sq);
            out_str.push('\n');
        }
        fs::write(&dec_fasta_out, out_str).unwrap_or_else(|e| {
            eprintln!("Error writing decompressed FASTA to '{dec_fasta_out}': {e}");
            process::exit(1);
        });
        println!(
            "Decompressed {} taxa x {} sites from '{}' -> '{}'",
            names.len(),
            seqs.first().map(|s| s.len()).unwrap_or(0),
            fasta,
            dec_fasta_out
        );
        if let Some(nwk_path) = out_newick {
            fs::write(&nwk_path, format!("{newick}\n")).unwrap_or_else(|e| {
                eprintln!("Error writing Newick to '{nwk_path}': {e}");
                process::exit(1);
            });
            println!("Extracted embedded Newick tree to: {nwk_path}");
        }
        return;
    }

    let (names, seqs) = parse_alignment_file(Path::new(&fasta)).unwrap_or_else(|e| {
        eprintln!("Error reading alignment '{fasta}': {e}");
        process::exit(1);
    });

    let engine = BushwhackerEngine::new(!no_coding, !no_shunts, min_shunt_sites);
    let res = engine.infer(&names, &seqs).unwrap_or_else(|e| {
        eprintln!("Error running BUSHWHACKER on '{fasta}': {e}");
        process::exit(1);
    });

    println!("{}", "=".repeat(78));
    println!(
        "BUSHWHACKER Summary: {}  (N={} taxa, L={} sites)",
        fasta, res.n_taxa, res.n_sites
    );
    println!("{}", "-".repeat(78));
    println!("  Raw Uncompressed DL   : {:12.1} nats", res.raw_uncompressed_nats);
    println!("  Star Tree MDL (Null)  : {:12.1} nats", res.star_tree_nats);
    println!(
        "  Induced Tree MDL      : {:12.1} nats  ({:.2}% of Star)",
        res.backbone_tree_nats,
        100.0 * res.backbone_tree_nats / res.star_tree_nats.max(1e-9)
    );
    println!(
        "  Induced Network MDL   : {:12.1} nats  ({:.2}% of Star)",
        res.full_graph_nats,
        100.0 * res.full_graph_nats / res.star_tree_nats.max(1e-9)
    );
    println!(
        "  Resolved Splits       : {:4} / {:4} binary splits ({} polytomy hubs)",
        res.n_internal_branches, res.max_binary_internal_branches, res.polytomy_nodes
    );
    println!("  Certified Shunts      : {:4}", res.shunts.len());
    for (idx, s) in res.shunts.iter().enumerate() {
        let codons: Vec<usize> = s.sites.iter().map(|&x| x / 3 + 1).collect();
        println!(
            "    [{}] {}  {} <-> {}  (+{:.2} nats, P={:.1}%, m={} sites, codons={:?})",
            idx + 1,
            s.mode,
            s.label_u,
            s.label_v,
            s.net_mdl_gain_nats,
            100.0 * s.posterior_prob,
            s.shared_mutations,
            codons
        );
    }
    println!("  Runtime               : {:.4} s", res.runtime_sec);
    println!("{}", "=".repeat(78));

    if let Some(nwk_path) = out_newick {
        fs::write(&nwk_path, format!("{}\n", res.newick)).unwrap_or_else(|e| {
            eprintln!("Error writing Newick to '{nwk_path}': {e}");
            process::exit(1);
        });
        println!("Wrote Newick tree to: {nwk_path}");
    } else if out_compress.is_none() {
        println!("\nNewick Tree:\n{}", res.newick);
    }

    if let Some(json_path) = out_json {
        let json_str = serialize_json(&res, &fasta);
        fs::write(&json_path, json_str).unwrap_or_else(|e| {
            eprintln!("Error writing JSON to '{json_path}': {e}");
            process::exit(1);
        });
        println!("Wrote JSON report to: {json_path}");
    }

    if let Some(bwz_path) = out_compress {
        let (bwz_bytes, stats) =
            compress_alignment_to_bwz(&names, &seqs, &res).unwrap_or_else(|e| {
                eprintln!("Error compressing alignment to BWZ: {e}");
                process::exit(1);
            });
        fs::write(&bwz_path, &bwz_bytes).unwrap_or_else(|e| {
            eprintln!("Error writing BWZ archive to '{bwz_path}': {e}");
            process::exit(1);
        });
        println!(
            "Wrote lossless BWZ archive to: {bwz_path}\n  Sequence Payload : {:8} B ({:9.1} nats, {:5.1}x vs raw ASCII {:8} B)\n  Taxon Headers    : {:8} B\n  Total Archive    : {:8} B",
            stats.compressed_seq_bytes,
            stats.seq_mdl_nats,
            stats.raw_seq_bytes as f64 / stats.compressed_seq_bytes.max(1) as f64,
            stats.raw_seq_bytes,
            stats.compressed_header_bytes,
            stats.total_bwz_bytes
        );
    }
}
