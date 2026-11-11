/* tslint:disable */
/* eslint-disable */

export class WasmBushwhacker {
    free(): void;
    [Symbol.dispose](): void;
    /**
     * Compress an alignment or character matrix (`fasta_text`) and its JSON-serialized
     * `BushwhackerResult` into a self-contained, bit-for-bit lossless `.bwz` binary archive.
     */
    compress_bwz(fasta_text: string, result_json: string): Uint8Array;
    /**
     * Decompress a `.bwz` binary archive back into a JSON object `{ fasta, newick, n_taxa, n_sites }`.
     */
    decompress_bwz(bwz_bytes: Uint8Array | number[]): string;
    /**
     * Induce an MDL phylogenetic graph with explicit or auto-detected data modality
     * (`auto`, `nucleotide`, `protein`, `translate`, `binary`, `discrete`, `ordered`).
     */
    infer_alignment(input_text: string, data_mode: string, coding_aware: boolean, detect_shunts: boolean, min_shunt_sites: number): string;
    /**
     * Backwards-compatible entrypoint that auto-detects Nucleotide, Protein, Binary (`0/1`),
     * or Discrete/Ordered CSV/TSV matrices.
     */
    infer_fasta(fasta_text: string, coding_aware: boolean, detect_shunts: boolean, min_shunt_sites: number): string;
    constructor();
}

export type InitInput = RequestInfo | URL | Response | BufferSource | WebAssembly.Module;

export interface InitOutput {
    readonly memory: WebAssembly.Memory;
    readonly __wbg_wasmbushwhacker_free: (a: number, b: number) => void;
    readonly wasmbushwhacker_compress_bwz: (a: number, b: number, c: number, d: number, e: number) => [number, number, number, number];
    readonly wasmbushwhacker_decompress_bwz: (a: number, b: number, c: number) => [number, number, number, number];
    readonly wasmbushwhacker_infer_alignment: (a: number, b: number, c: number, d: number, e: number, f: number, g: number, h: number) => [number, number, number, number];
    readonly wasmbushwhacker_infer_fasta: (a: number, b: number, c: number, d: number, e: number, f: number) => [number, number, number, number];
    readonly wasmbushwhacker_new: () => number;
    readonly __wbindgen_externrefs: WebAssembly.Table;
    readonly __wbindgen_malloc: (a: number, b: number) => number;
    readonly __wbindgen_realloc: (a: number, b: number, c: number, d: number) => number;
    readonly __wbindgen_free: (a: number, b: number, c: number) => void;
    readonly __externref_table_dealloc: (a: number) => void;
    readonly __wbindgen_start: () => void;
}

export type SyncInitInput = BufferSource | WebAssembly.Module;

/**
 * Instantiates the given `module`, which can either be bytes or
 * a precompiled `WebAssembly.Module`.
 *
 * @param {{ module: SyncInitInput }} module - Passing `SyncInitInput` directly is deprecated.
 *
 * @returns {InitOutput}
 */
export function initSync(module: { module: SyncInitInput } | SyncInitInput): InitOutput;

/**
 * If `module_or_path` is {RequestInfo} or {URL}, makes a request and
 * for everything else, calls `WebAssembly.instantiate` directly.
 *
 * @param {{ module_or_path: InitInput | Promise<InitInput> }} module_or_path - Passing `InitInput` directly is deprecated.
 *
 * @returns {Promise<InitOutput>}
 */
export default function __wbg_init (module_or_path?: { module_or_path: InitInput | Promise<InitInput> } | InitInput | Promise<InitInput>): Promise<InitOutput>;
