
export interface ScanResult {
  xyz: number;
  hms_dms: string;
  ra: number;
  dec: number;
  sync: number;
  pull: number;
  floor: number;
  dist: number;
  simbad_id: string;
}

export interface ScannerParams {
  start_xyz: number;
  search_radius: number;
  base_ra: number;
  base_dec: number;
}

export interface AIAnalysis {
  summary: string;
  hypotheses: string[];
  anomalies: string;
}
