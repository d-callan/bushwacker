// BUSHWHACKER // Dedicated WebAssembly Worker (Off-Main-Thread MDL Engine)
import init, { WasmBushwhacker } from './pkg/bushwhacker.js';

let engine = null;
let ready = false;

async function ensureInitialized() {
  if (ready && engine) return;
  await init();
  engine = new WasmBushwhacker();
  ready = true;
}

self.onmessage = async (ev) => {
  const msg = ev.data || {};
  const { id, type } = msg;
  try {
    if (type === 'init') {
      await ensureInitialized();
      self.postMessage({ id, type: 'ready' });
      return;
    }

    await ensureInitialized();

    if (type === 'infer_alignment' || type === 'infer_fasta') {
      const t0 = performance.now();
      const mode = String(msg.dataMode || 'auto');
      const jsonStr = typeof engine.infer_alignment === 'function'
        ? engine.infer_alignment(
            String(msg.fastaText || ''),
            mode,
            Boolean(msg.codingAware),
            Boolean(msg.detectShunts),
            Number(msg.minShuntSites || 4)
          )
        : engine.infer_fasta(
            String(msg.fastaText || ''),
            Boolean(msg.codingAware),
            Boolean(msg.detectShunts),
            Number(msg.minShuntSites || 4)
          );
      const wallMs = performance.now() - t0;
      self.postMessage({ id, type: 'infer_result', jsonStr, wallMs });
      return;
    }

    if (type === 'compress_bwz') {
      const bytes = engine.compress_bwz(
        String(msg.fastaText || ''),
        String(msg.resultJson || '')
      );
      self.postMessage({ id, type: 'compress_bwz_result', bytes }, [bytes.buffer]);
      return;
    }

    if (type === 'decompress_bwz') {
      const u8 = msg.bytes instanceof Uint8Array ? msg.bytes : new Uint8Array(msg.bytes);
      const jsonStr = engine.decompress_bwz(u8);
      self.postMessage({ id, type: 'decompress_bwz_result', jsonStr });
      return;
    }

    throw new Error(`Unknown worker command: ${type}`);
  } catch (err) {
    const errMsg = (err && err.message) ? err.message : String(err);
    self.postMessage({ id, type: 'error', error: errMsg });
  }
};
