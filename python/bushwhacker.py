#!/usr/bin/env python3
"""
BUSHWHACKER CLI & Library Entrypoint
====================================
Bipartition Unification, Shunt Hunting, & Weighted Hub Aggregation via
Combinadic Kolmogorov Entropy Reduction.
"""

from mdl_phylo import (
    BranchRecord,
    BushwhackerEngine,
    BushwhackerResult,
    MDLGraphResult,
    MDLPhyloEngine,
    ShuntEdge,
    branch_mutation_code_nats,
    detect_coding_frame_and_weights,
    encode_nucleotide_alignment,
    log_binom,
    log_binom_vec,
    main,
    parse_fasta,
    split_topology_cost_nats,
    universal_int_nats,
)

__all__ = [
    "BushwhackerEngine",
    "BushwhackerResult",
    "MDLPhyloEngine",
    "MDLGraphResult",
    "BranchRecord",
    "ShuntEdge",
    "parse_fasta",
    "encode_nucleotide_alignment",
    "detect_coding_frame_and_weights",
    "log_binom",
    "log_binom_vec",
    "universal_int_nats",
    "split_topology_cost_nats",
    "branch_mutation_code_nats",
    "main",
]

if __name__ == "__main__":
    main()
