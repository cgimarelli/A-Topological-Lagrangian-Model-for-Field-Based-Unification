import numpy as np
import warnings
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u

Simbad.TIMEOUT = 300
Simbad.add_votable_fields('main_id', 'otype', 'flux(V)')
warnings.filterwarnings('ignore')
def dark_sector_scanner(start_xyz, search_radius, base_ra, base_dec):
    N_ensemble = 25
    LY_per_unit = 100 
    granular_steps = [1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3]    
    scan_points = np.linspace(start_xyz, start_xyz + search_radius, 100)
    results = []
    print(f"\n--- INITIATING DARK SECTOR SCAN ---")
    print(f"SECTOR: RA {base_ra} | DEC {base_dec}")
    print(f"DEPTH:  {start_xyz} to {start_xyz + search_radius} XYZ Units")
    print("DATA CAP: OFF (Raw Metric Tension)\n")
    for p in scan_points:        
        seed_val = (int(p * 10000) % 4294967295) + 1
        np.random.seed(seed_val)        
        phases = np.random.normal(0, np.random.rand(), N_ensemble)
        complex_mean = np.mean(np.exp(1j * phases))
        sync_val = np.abs(complex_mean)      
        hard_sync = np.clip(sync_val / (1.0 - (sync_val * 0.5)), sync_val, 0.9999)
        if hard_sync > 0.0: 
            phase_angle = np.angle(complex_mean)
            normalized_phase = (phase_angle + np.pi) / (2 * np.pi)
            floor = normalized_phase * N_ensemble            
            floor_delta = abs(floor - 12.00) + 0.01
            pull = (hard_sync**2) / np.sqrt(floor_delta)            
            dist_ly = p * LY_per_unit            
            if start_xyz > 5000:
                ra = base_ra
                dec = base_dec
            else:
                ra = (base_ra + (p * 0.05)) % 24
                dec = np.clip(base_dec + (p * 0.02), -90, 90)            
            peak_coord = SkyCoord(ra=ra*u.hourangle, dec=dec*u.deg)            
            simbad_res = "GRAVITATIONAL INVERSION | NO OBJECT IN SIMBAD"
            is_dark_sector = 11.0 <= floor <= 13.5            
            simbad_res = "GRAVITATIONAL INVERSION | NO OBJECT IN SIMBAD"
            is_dark_sector = 11.0 <= floor <= 13.5            
            for r in granular_steps:
                try:
                    local_map = Simbad.query_region(peak_coord, radius=r * u.arcmin)
                    if local_map:
                        colnames = local_map.colnames
                        id_col = next((c for c in colnames if c.lower() in ['main_id', 'id']), colnames[0])
                        obj_name = local_map[0][id_col]
                        if hasattr(obj_name, 'decode'): obj_name = obj_name.decode()                        
                        if is_dark_sector:
                            simbad_res = f"ARC MINS {r}': {obj_name}"
                        else:
                            simbad_res = f"MULTI-DIMENSIONAL OBJECT DETECTED @ {r}': {obj_name}"
                        break
                except Exception:
                    simbad_res = "SCAN ERROR (QUERY TIMEOUT)"
                    break
            results.append({
                'xyz': p, 
                'hms_dms': peak_coord.to_string('hmsdms'),
                'sync': hard_sync, 
                'pull': pull, 
                'floor': floor, 
                'dist': dist_ly, 
                'simbad_id': simbad_res
            })

    return results

##########################################################################
###### Adjust base_ra, base_dec, start_xyz, and search_radius below ######
##########################################################################

final_hits = dark_sector_scanner(
    base_ra=14.83, 
    base_dec=46.00,
    start_xyz=7000000, 
    search_radius=2500000,     
)
header = f"{'XYZ':<12} | {'Coordinate Address':<25} | {'PULL':<8} | {'Sync':<7} | {'FLOOR':<7} | {'SIMBAD RESOLUTION'}"
print("\n" + header)
print("-" * 150)
for h in final_hits:
    print(f"{h['xyz']:<12.0f} | {h['hms_dms']:<25} | {h['pull']:<8.4f} | {h['sync']:<7.4f} | {h['floor']:<7.2f} | {h['simbad_id']}")
