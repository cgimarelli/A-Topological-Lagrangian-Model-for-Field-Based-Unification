import sys
import threading
import tkinter as tk
import customtkinter as ctk
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import hashlib
import warnings
import logging


try:
    from astroquery.simbad import Simbad
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.exceptions import NoResultsWarning

    warnings.filterwarnings("ignore")
    logging.getLogger('astroquery').setLevel(logging.ERROR)

    Simbad.add_votable_fields('main_id', 'otype', 'ra(d)', 'dec(d)')
    Simbad.TIMEOUT = 120
    SIMBAD_AVAILABLE = True
except ImportError:
    SIMBAD_AVAILABLE = False
    print("!! WARNING: 'astroquery' or 'astropy' not installed. Simbad lookups will be disabled.")

PLANCK_LENGTH = 1.616255e-35
HUBBLE_CONST = 2.27e-18
C_LIGHT = 2.99792458e8
LAMBDA_COSMO = 1.7e-52
N_DIRAC = ( (C_LIGHT / (HUBBLE_CONST * PLANCK_LENGTH)) )**(2/3)
BETA_STIFFNESS = LAMBDA_COSMO * N_DIRAC
GAMMA_MODULUS = PLANCK_LENGTH**2

def perform_simbad_lookup(coord_str):
    if SIMBAD_AVAILABLE == False: return "SIMBAD N/A", "N/A", None, None
    try:
        c = SkyCoord(coord_str, frame='icrs')
        result_table = Simbad.query_region(c, radius=5 * u.arcmin)
        if result_table:
            raw_cols = result_table.colnames
            def find_col(possible_names):
                for name in possible_names:
                    for raw in raw_cols:
                        if name.lower() in raw.lower(): return raw
                return None
            id_col = find_col(['main_id', 'id', 'name']) or raw_cols[0]
            ra_col = find_col(['ra(d)', 'ra_d', 'ra'])
            dec_col = find_col(['dec(d)', 'dec_d', 'dec'])
            otype_col = find_col(['otype', 'type'])
            ras = np.array(result_table[ra_col]).astype(float)
            decs = np.array(result_table[dec_col]).astype(float)
            res_coords = SkyCoord(ra=ras*u.deg, dec=decs*u.deg, frame='icrs')
            seps = c.separation(res_coords)
            idx = np.argmin(seps)
            row = result_table[idx]
            name = str(row[id_col])
            obj_type = str(row[otype_col]) if otype_col else "Object"
            return f"MATCH: {name} ({obj_type})", f"{seps[idx].arcsecond:.2f}\"", ras[idx], decs[idx]
        return "NO MATCH", "N/A", None, None
    except Exception:
        return "SIMBAD FAIL", "N/A", None, None

class VolumetricScanner:
    def __init__(self, N_ensemble=25, shear_modulus=1.25):
        self.N = N_ensemble
        self.shear = shear_modulus
        self.BARYONIC_FLOOR = 12.0

    def hash_coordinate(self, x, y, z):
        coord_str = f"{x:.4f}|{y:.4f}|{z:.4f}"
        hash_obj = hashlib.md5(coord_str.encode())
        return int(hash_obj.hexdigest(), 16) % 4294967295

    def get_vacuum_state(self, x, y, z):
        seed_val = self.hash_coordinate(x, y, z)
        np.random.seed(seed_val)
        limit = (np.pi / 4.0) * self.shear
        phases = np.random.uniform(-limit, limit, self.N)
        complex_mean = np.mean(np.exp(1j * phases))
        sync_val = np.abs(complex_mean)
        hard_sync = np.clip(sync_val / (1.0 - (sync_val * 0.5)), sync_val, 0.9999)
        return hard_sync, complex_mean

    def get_topological_invariants(self, x, y, z):
        eps = 0.02
        _, z0 = self.get_vacuum_state(x, y, z)
        theta0 = np.angle(z0)
        _, zx = self.get_vacuum_state(x + eps, y, z)
        _, zy = self.get_vacuum_state(x, y + eps, z)
        _, zz = self.get_vacuum_state(x, y, z + eps)

        I_vec = np.array([(np.angle(zx) - theta0) / eps,
                          (np.angle(zy) - theta0) / eps,
                          (np.angle(zz) - theta0) / eps])
        shear_vel = np.linalg.norm(I_vec)

        _, zyx = self.get_vacuum_state(x + eps, y + eps, z)
        vorticity_z = ((np.angle(zyx) - np.angle(zx)) / eps) - ((np.angle(zyx) - np.angle(zy)) / eps)
        helicity = shear_vel * abs(vorticity_z)
        return shear_vel, helicity

    def analyze_topology(self, hard_sync, complex_mean):
        if hard_sync <= 0.0: return 0, 0, 0
        phase_angle = np.angle(complex_mean)
        normalized_phase = (phase_angle + np.pi) / (2 * np.pi)
        floor = normalized_phase * self.N
        floor_delta = abs(floor - self.BARYONIC_FLOOR) + 0.01
        pull = (hard_sync**2) / np.sqrt(floor_delta)
        m_raw = pull / (self.N * (hard_sync**2))
        m = m_raw * 5.5
        d_confinement = (GAMMA_MODULUS / BETA_STIFFNESS)**0.25 * 1e15
        rms_radius_fm = 0.775 * (d_confinement / 2.0)
        return pull, m, rms_radius_fm

    def classify_matter(self, m, hard_sync, kinetic, helicity, pull):
        if pull < 1.0 and hard_sync > 0.99: return "Minkowskian Baseline"
        if pull > 2.0 and m < 0.15: return "Non-Local Metric Deformation"
        if helicity > 50.0:
            if m > 1.5: return "Spacetime Curvature Limit (Black Hole)"
            elif kinetic > 20.0: return "Localized Torsional Inflection"
        if hard_sync > 0.95:
            if 0.28 <= m < 0.32: return "Solenoidal Metric Schema (Iron)"
            if 0.32 <= m < 0.36: return "Helical Metric Schema"
            if 0.36 <= m < 0.45: return "Trefoil Metric Schema"
            if 0.80 <= m < 1.20: return "Baryonic Metric Substrate"
        if 0.15 <= m < 0.28: return "Gasious Substrate Schema"
        return "Metric Spacetime Anomaly"

    def cartesian_to_astronomical(self, x, y, z):
        r = np.sqrt(x**2 + y**2 + z**2) or 1e-9
        dec_rad = np.arcsin(z / r)
        ra_rad = np.arctan2(y, x)
        ra_deg = (np.degrees(ra_rad) + 360) % 360
        dec_deg = np.degrees(dec_rad)
        alpha_inv = 137.035999084
        drag_shift = (1.5 / self.N) * (1.0 / alpha_inv)
        ra_deg = (ra_deg + drag_shift) % 360
        ra_h = int(ra_deg / 15.0)
        ra_m = int((ra_deg / 15.0 - ra_h) * 60)
        ra_s = ((ra_deg / 15.0 - ra_h) * 60 - ra_m) * 60
        sign = "+" if dec_deg >= 0 else "-"
        abs_dec = abs(dec_deg)
        dec_d = int(abs_dec)
        dec_m = int((abs_dec - dec_d) * 60)
        dec_s = ((abs_dec - dec_d) * 60 - dec_m) * 60
        return f"{ra_h:02d}h{ra_m:02d}m{ra_s:05.2f}s {sign}{dec_d:02d}d{dec_m:02d}m{dec_s:05.2f}s"

    def scan_sector(self, center_coord, radius):
        volume = (4/3) * np.pi * (radius**3)
        num_points = int(volume * 0.1)
        if num_points > 8000: num_points = 8000

        print(f"Scanning Volume... Sampling {num_points} points.")
        rng = np.random.default_rng(int(abs(sum(center_coord) * radius)) % 4294967295)

        phi = rng.uniform(0, 2*np.pi, num_points)
        costheta = rng.uniform(-1, 1, num_points)
        u_rand = rng.uniform(0, 1, num_points)
        theta = np.arccos(costheta)
        r_dist = radius * np.cbrt(u_rand)

        xs = r_dist * np.sin(theta) * np.cos(phi)
        ys = r_dist * np.sin(theta) * np.sin(phi)
        zs = r_dist * np.cos(theta)

        xs = center_coord[0] + xs
        ys = center_coord[1] + ys
        zs = center_coord[2] + zs

        plot_data = []
        for i in range(num_points):
            k, z_complex = self.get_vacuum_state(xs[i], ys[i], zs[i])
            kin, hel = self.get_topological_invariants(xs[i], ys[i], zs[i])
            if k > 0.6 or hel > 4.0:
                pull, m, rad = self.analyze_topology(k, z_complex)
                identity = self.classify_matter(m, k, kin, hel, pull)
                coord_str = self.cartesian_to_astronomical(xs[i], ys[i], zs[i])

                # Re-calculate simple floats for the 2D plot
                r = np.sqrt(xs[i]**2 + ys[i]**2 + zs[i]**2) or 1e-9
                dec_rad = np.arcsin(zs[i] / r)
                ra_rad = np.arctan2(ys[i], xs[i])
                ra_deg = (np.degrees(ra_rad) + 360) % 360
                alpha_inv = 137.035999084
                drag_shift = (1.5 / self.N) * (1.0 / alpha_inv)
                ra_deg = (ra_deg + drag_shift) % 360
                ra_h = ra_deg / 15.0
                dec_d = np.degrees(dec_rad)

                plot_data.append({
                    'ra': ra_h,
                    'dec': dec_d,
                    'mode': m,
                    'radius': rad,
                    'matter_id': identity,
                    'helicity': hel,
                    'coord_str': coord_str
                })
        return plot_data

# --- PLOTTING FUNCTION ---
def plot_spectral_map(results):
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(12, 10))
    fig.patch.set_facecolor('#0b0e14')
    ax.set_facecolor('#0b0e14')

    ra = [r['ra'] for r in results]
    dec = [r['dec'] for r in results]
    modes = [r['mode'] for r in results]
    radii = [r['radius'] * 60 for r in results]

    sc = ax.scatter(ra, dec, s=radii, c=modes, cmap='nipy_spectral', alpha=0.8, edgecolors='none')

    iron_ra = [r['ra'] for r in results if "Solenoidal" in r['matter_id'] or "Iron" in r['matter_id']]
    iron_dec = [r['dec'] for r in results if "Solenoidal" in r['matter_id'] or "Iron" in r['matter_id']]
    if iron_ra:
        ax.scatter(iron_ra, iron_dec, s=250, facecolors='none', edgecolors='lime', linewidth=2, label='Stable Cores (Iron)')

    bh_ra = [r['ra'] for r in results if "Black Hole" in r['matter_id']]
    bh_dec = [r['dec'] for r in results if "Black Hole" in r['matter_id']]
    if bh_ra:
        ax.scatter(bh_ra, bh_dec, s=150, marker='x', color='red', linewidth=3, label='Metric Fractures (BH)')

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Topological Mode (m) - [Spectral Stress]')
    ax.set_title("VOLUMETRIC TOPOLOGY MAP: Projected Density", loc='left', color='cyan', fontsize=16, fontweight='bold')
    ax.set_xlabel("Right Ascension (Hours)", fontsize=12, color='white')
    ax.set_ylabel("Declination (Degrees)", fontsize=12, color='white')
    ax.invert_xaxis()
    ax.legend(loc='upper right', frameon=True, facecolor='#1a1a1a', edgecolor='cyan')
    ax.grid(True, color='#2a2a35', alpha=0.3)
    plt.tight_layout()
    plt.show()

# --- EXECUTION WRAPPER ---
def execute_volumetric_scan(base_ra, base_dec, start_xyz, search_radius):
    print(f"--- INITIATING VOLUMETRIC METRIC SCAN ---")
    print(f"Target: RA {base_ra}, DEC {base_dec}")
    print(f"Depth: {start_xyz} XYZ | Radius: {search_radius}")

    # Convert RA/Dec to Cartesian Center for the Scanner
    ra_rad = np.radians(base_ra * 15)
    dec_rad = np.radians(base_dec)
    x_c = start_xyz * np.cos(dec_rad) * np.cos(ra_rad)
    y_c = start_xyz * np.cos(dec_rad) * np.sin(ra_rad)
    z_c = start_xyz * np.sin(dec_rad)

    # Run Scan
    scanner = VolumetricScanner()
    scan_results = scanner.scan_sector((x_c, y_c, z_c), search_radius)

    if not scan_results:
        print("No structure detected in this volume.")
    else:
        print(f"Structure Detected: {len(scan_results)} metric nodes.")

        # --- 1. PRINT EVERYTHING FIRST ---
        # AUDIT HEADERS
        print("\n" + "="*140)
        print("FULL SIMBAD REALITY CHECK - ALL NODES")
        print("="*140)
        print(f"{'PREDICTED COORDINATES':<32} | {'MODE (m)':<8} | {'HELICITY':<8} | {'PREDICTED IDENTITY':<35} | {'OFFSET':<10} | {'REALITY CHECK'}")
        print("-" * 140)

        # AUDIT LOOP (Print the rows)
        for node in scan_results:
            coord = node['coord_str']
            mode = node['mode']
            hel = node['helicity']
            ident = node['matter_id']
            sim_res, off_arc, real_ra, real_dec = perform_simbad_lookup(coord)

            # Formatted String for GUI box
            row = f"{coord:<32} | {mode:<8.4f} | {hel:<8.2f} | {ident:<35} | {off_arc:<10} | {sim_res}"
            print(row)

        # --- 2. OPEN THE PLOT LAST ---
        # Now that text is printed, we open the plot.
        # Ensure your plot_spectral_map function uses plt.show() (blocking) at the end.
        print("\n--- Opening Visualization ---")
        plot_spectral_map(scan_results)


########################################################################################
############################# BOTTOM WRAPPER DO NOT MODIFY #############################
########################################################################################

class PrintLogger:
    def __init__(self, text_widget):
        self.text_widget = text_widget
        self.terminal = sys.stdout

    def write(self, message):
        # 1. Safely write to the "real" console only if it exists
        if self.terminal is not None:
            try:
                self.terminal.write(message)
                self.terminal.flush()
            except Exception:
                pass # If the console is broken/missing, just ignore it

        # 2. Write to the GUI Textbox (This always works)
        try:
            self.text_widget.configure(state="normal") # Unlock
            self.text_widget.insert("end", message)    # Write
            self.text_widget.see("end")                # Scroll to bottom
            self.text_widget.configure(state="disabled") # Lock
        except Exception:
            pass

    def flush(self):
        # Safely flush only if terminal exists
        if self.terminal is not None:
            try:
                self.terminal.flush()
            except Exception:
                pass
        self.terminal.flush()

ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class ScannerApp(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Space Scanner | Desktop v3")
        self.geometry("1100x700")

        # Layout: 2 Columns. Left (Inputs), Right (Console Output)
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # --- LEFT SIDEBAR (INPUTS) ---
        self.sidebar_frame = ctk.CTkFrame(self, width=200, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")

        self.logo_label = ctk.CTkLabel(self.sidebar_frame, text="SCANNER INPUTS", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_label.grid(row=0, column=0, padx=20, pady=(20, 10))

        # RA Input
        self.ra_label = ctk.CTkLabel(self.sidebar_frame, text="Right Ascension (Hrs):", anchor="w")
        self.ra_label.grid(row=1, column=0, padx=20, pady=(10, 0))
        self.ra_entry = ctk.CTkEntry(self.sidebar_frame, placeholder_text="e.g. 5")
        self.ra_entry.grid(row=2, column=0, padx=20, pady=(0, 10))
        self.ra_entry.insert(0, "5")

        # Dec Input
        self.dec_label = ctk.CTkLabel(self.sidebar_frame, text="Declination (Deg):", anchor="w")
        self.dec_label.grid(row=3, column=0, padx=20, pady=(10, 0))
        self.dec_entry = ctk.CTkEntry(self.sidebar_frame, placeholder_text="e.g. -5")
        self.dec_entry.grid(row=4, column=0, padx=20, pady=(0, 10))
        self.dec_entry.insert(0, "-5")

        # Distance Input
        self.dist_label = ctk.CTkLabel(self.sidebar_frame, text="Distance (XYZ units):", anchor="w")
        self.dist_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.dist_entry = ctk.CTkEntry(self.sidebar_frame, placeholder_text="e.g. 200")
        self.dist_entry.grid(row=6, column=0, padx=20, pady=(0, 10))
        self.dist_entry.insert(0, "200")

        # Radius Input
        self.rad_label = ctk.CTkLabel(self.sidebar_frame, text="Scan Radius:", anchor="w")
        self.rad_label.grid(row=7, column=0, padx=20, pady=(10, 0))
        self.rad_entry = ctk.CTkEntry(self.sidebar_frame, placeholder_text="e.g. 8")
        self.rad_entry.grid(row=8, column=0, padx=20, pady=(0, 10))
        self.rad_entry.insert(0, "8")

        # RUN BUTTON
        self.run_button = ctk.CTkButton(self.sidebar_frame, text="INITIATE SCAN", command=self.start_scan_thread)
        self.run_button.grid(row=9, column=0, padx=20, pady=30)

        # --- RIGHT SIDE (CONSOLE OUTPUT) ---
        self.console_label = ctk.CTkLabel(self, text="DATA LOG & SIMBAD CHECK", font=ctk.CTkFont(size=14, weight="bold"))
        self.console_label.grid(row=0, column=1, padx=20, pady=(10,0), sticky="nw")

        # The Text Box (Monospaced font is crucial for your table alignment)
        self.textbox = ctk.CTkTextbox(self, width=800, font=("Consolas", 12))
        self.textbox.grid(row=0, column=1, padx=(20, 20), pady=(40, 20), sticky="nsew")
        self.textbox.configure(state="disabled") # Read only initially

        # Redirect sys.stdout to this textbox
        sys.stdout = PrintLogger(self.textbox)

    def start_scan_thread(self):
        # Disable button so user doesn't double click
        self.run_button.configure(state="disabled", text="SCANNING...")

        # Run the heavy math in a separate thread so GUI doesn't freeze
        threading.Thread(target=self.run_scan, daemon=True).start()

    def run_scan(self):
        try:
            # Get values
            base_ra = float(self.ra_entry.get())
            base_dec = float(self.dec_entry.get())
            start_xyz = float(self.dist_entry.get())
            search_radius = float(self.rad_entry.get())

            # Execute your existing logic
            execute_volumetric_scan(base_ra, base_dec, start_xyz, search_radius)

        except ValueError:
            print("ERROR: Please enter valid numbers for coordinates.")
        except Exception as e:
            print(f"CRITICAL ERROR: {e}")
        finally:
            # Re-enable button
            self.run_button.configure(state="normal", text="INITIATE SCAN")

if __name__ == "__main__":
    app = ScannerApp()
    app.mainloop()
