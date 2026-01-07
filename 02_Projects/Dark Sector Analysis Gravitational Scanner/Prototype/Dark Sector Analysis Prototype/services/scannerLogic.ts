
import { ScanResult, ScannerParams } from '../types';

/**
 * Standard MT19937 Implementation (Mersenne Twister)
 * This matches the initialization used in NumPy and other scientific libraries.
 * Internal state is 32-bit as per specification, but results are handled as 64-bit floats.
 */
class MersenneTwister {
  private mt = new Uint32Array(624);
  private index = 0;

  constructor(seed: number) {
    this.mt[0] = seed >>> 0;
    for (this.index = 1; this.index < 624; this.index++) {
      const s = this.mt[this.index - 1] ^ (this.mt[this.index - 1] >>> 30);
      this.mt[this.index] = (Math.imul(1812433253, s) + this.index) >>> 0;
    }
  }

  private twist() {
    for (let i = 0; i < 624; i++) {
      const y = (this.mt[i] & 0x80000000) | (this.mt[(i + 1) % 624] & 0x7fffffff);
      this.mt[i] = this.mt[(i + 397) % 624] ^ (y >>> 1);
      if (y % 2 !== 0) this.mt[i] ^= 0x9908b0df;
    }
    this.index = 0;
  }

  next(): number {
    if (this.index >= 624) this.twist();
    let y = this.mt[this.index++];
    y ^= y >>> 11;
    y ^= (y << 7) & 0x9d2c5680;
    y ^= (y << 15) & 0xefc60000;
    y ^= y >>> 18;
    return (y >>> 0) / 4294967296.0;
  }

  nextNormal(): number {
    let u = 0, v = 0;
    while (u === 0) u = this.next();
    while (v === 0) v = this.next();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
  }
}

/**
 * Converts RA/Dec to HMS/DMS format matching astropy's hmsdms precision.
 */
function formatCoordinates(ra: number, dec: number): string {
  const raVal = ((ra % 24) + 24) % 24;
  const h = Math.floor(raVal);
  const m_exact = (raVal - h) * 60;
  const m = Math.floor(m_exact);
  const s = Math.round((m_exact - m) * 60);
  
  const sign = dec >= 0 ? '+' : '-';
  const absDec = Math.abs(dec);
  const d = Math.floor(absDec);
  const dm_exact = (absDec - d) * 60;
  const dm = Math.floor(dm_exact);
  const ds = Math.round((dm_exact - dm) * 60);

  return `${h.toString().padStart(2, '0')}h${m.toString().padStart(2, '0')}m${s.toString().padStart(2, '0')}s ${sign}${d.toString().padStart(2, '0')}d${dm.toString().padStart(2, '0')}m${ds.toString().padStart(2, '0')}s`;
}

export const runN12Scanner = (params: ScannerParams): ScanResult[] => {
  const { start_xyz, search_radius, base_ra, base_dec } = params;
  const N_ensemble = 25;
  const LY_per_unit = 100;
  const steps = 100;
  const results: ScanResult[] = [];

  for (let i = 0; i < steps; i++) {
    // Force 64-bit float precision for p calculation
    const p: number = start_xyz + (search_radius * i) / (steps - 1);
    
    // Procedural Seed Locking (Numpy style)
    const seed_val = (Math.floor(p * 10000) % 4294967295) + 1;
    const rng = new MersenneTwister(seed_val);
    
    // np.random.rand() generates the scale for the normal distribution
    const rand_scale = rng.next();
    
    let sumReal = 0;
    let sumImag = 0;

    for (let j = 0; j < N_ensemble; j++) {
      // phases = np.random.normal(0, np.random.rand(), N_ensemble)
      const phase = rng.nextNormal() * rand_scale;
      sumReal += Math.cos(phase);
      sumImag += Math.sin(phase);
    }
    
    // complex_mean = np.mean(np.exp(1j * phases))
    const meanReal = sumReal / N_ensemble;
    const meanImag = sumImag / N_ensemble;
    const sync_val = Math.sqrt(meanReal * meanReal + meanImag * meanImag);
    
    // Hard-Sync Gain Filter: np.clip(sync_val / (1.0 - (sync_val * 0.5)), sync_val, 0.9999)
    const hard_sync = Math.min(0.9999, Math.max(sync_val, sync_val / (1.0 - (sync_val * 0.5))));

    if (hard_sync > 0.0) {
      const phase_angle = Math.atan2(meanImag, meanReal);
      const normalized_phase = (phase_angle + Math.PI) / (2 * Math.PI);
      const floor = normalized_phase * N_ensemble;
      
      const floor_delta = Math.abs(floor - 12.00) + 0.01;
      const pull = (hard_sync ** 2) / Math.sqrt(floor_delta);

      // --- COORDINATE PROGRESS LOGIC ---
      const progress = search_radius > 0 ? (p - start_xyz) / search_radius : 0;
      
      // ra_linear = (base_ra + (p * 0.0005)) 
      // ra_wobble = np.sin(p * 0.2) * (1.0 + progress * 5.0)
      const ra_linear = base_ra + (p * 0.0005);
      const ra_wobble = Math.sin(p * 0.2) * (1.0 + progress * 5.0);
      const ra = ((ra_linear + ra_wobble) % 24 + 24) % 24;
      
      // dec_swing = np.cos(p * 0.08) * (2.0 + progress * 15.0)
      const dec_swing = Math.cos(p * 0.08) * (2.0 + progress * 15.0);
      let dec = base_dec + dec_swing;
      dec = Math.max(-89.0, Math.min(89.0, dec));

      const dist_ly = p * LY_per_unit;
      
      // Simbad Identification Logic
      let simbad_res = "GRAVITATIONAL INVERSION | NO OBJECT IN SIMBAD";
      
      // Mocked simulation of granular search result for high-tension points
      if (pull > 2.5) {
        simbad_res = `DEEP FIELD ANOMALY [TX-${Math.floor(p % 1000).toString().padStart(3, '0')}]`;
      } else if (start_xyz > 50000000) {
        simbad_res = "DEEP SECTOR VOID TARGET";
      }

      results.push({
        xyz: p,
        hms_dms: formatCoordinates(ra, dec),
        ra,
        dec,
        sync: hard_sync,
        pull: pull,
        floor: floor,
        dist: dist_ly,
        simbad_id: simbad_res
      });
    }
  }

  // Explicit sort by XYZ Depth
  return results.sort((a, b) => a.xyz - b.xyz);
};
