import numpy as np
import warnings
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u

Simbad.TIMEOUT = 300
# Updated 'flux(V)' to 'V' to remove DeprecationWarning
Simbad.add_votable_fields('main_id', 'otype', 'V')
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
            # Calculate progress (0 to 1) of the current search radius
            progress = (p - start_xyz) / search_radius if search_radius > 0 else 0
            ra_linear = (base_ra + (p * 0.0005)) 
            ra_wobble = np.sin(p * 0.2) * (1.0 + progress * 5.0)
            ra = (ra_linear + ra_wobble) % 24
            dec_swing = np.cos(p * 0.08) * (2.0 + progress * 15.0)
            dec = base_dec + dec_swing
            dec = np.clip(dec, -89.0, 89.0)
            peak_coord = SkyCoord(ra=ra*u.hourangle, dec=dec*u.deg)
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
                        simbad_res = f"{obj_name}"
                        break
                except Exception:
                    simbad_res = "SCAN ERROR (QUERY TIMEOUT)"
                    break                    
            results.append({
                'xyz': p,
                'ra': ra,
                'dec': dec,
                'hms_dms': peak_coord.to_string('hmsdms'),
                'sync': hard_sync,
                'pull': pull,
                'floor': floor,
                'simbad_id': simbad_res
            })
    return results
def generate_visuals(results):
    plt.style.use('dark_background')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    fig.patch.set_facecolor('#0b0e14')
    xyz = [r['xyz'] for r in results]
    sync = [r['sync'] for r in results]
    pull = [r['pull'] for r in results]
    ra = [r['ra'] for r in results]
    dec = [r['dec'] for r in results]    
    sc = ax1.scatter(ra, dec, s=[p*50 for p in pull], c=pull, cmap='plasma', alpha=0.7, edgecolors='white', linewidth=0.5)
    ax1.set_title("METRIC TENSION MAP (RA VS DEC)", color='#d1d5db', loc='left', fontsize=10)
    ax1.set_xlabel("RA (h)", color='#9ca3af')
    ax1.set_ylabel("DEC (°)", color='#9ca3af')
    ax1.grid(True, linestyle='--', alpha=0.2)    
    ax2.plot(xyz, sync, color='#8b5cf6', marker='o', markersize=4, linewidth=1.5, alpha=0.8)
    ax2.set_title("SYNC STABILITY VS. XYZ DISTANCE", color='#d1d5db', loc='left', fontsize=10)
    ax2.set_xlabel("XYZ Units", color='#9ca3af')
    ax2.set_ylabel("Synchronization", color='#9ca3af')
    ax2.set_ylim(0, 1.1)
    ax2.grid(True, linestyle='--', alpha=0.2)    
    plt.tight_layout()
    plt.show()

##########################################################################
###### Adjust base_ra, base_dec, start_xyz, and search_radius below ######
##########################################################################

final_results = dark_sector_scanner(
    base_ra=0,
    base_dec=0,
    start_xyz=0,
    search_radius=100000,
)
generate_visuals(final_results)
header = f"{'XYZ':<12} | {'Coordinate Address':<25} | {'PULL':<8} | {'Sync':<7} | {'FLOOR':<7} | {'SIMBAD RESOLUTION'}"
print("\n" + header)
print("-" * 150)
for h in final_results:
    print(f"{h['xyz']:<12.0f} | {h['hms_dms']:<25} | {h['pull']:<8.4f} | {h['sync']:<7.4f} | {h['floor']:<7.2f} | {h['simbad_id']}")
