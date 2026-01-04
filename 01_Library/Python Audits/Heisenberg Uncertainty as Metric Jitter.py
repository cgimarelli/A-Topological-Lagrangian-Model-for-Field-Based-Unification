import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, PillowWriter

# Use 'Agg' backend to ensure reliable saving without display errors
import matplotlib
matplotlib.use('Agg')

# High-fidelity dark theme
plt.style.use('dark_background')

# --- 1. PHYSICS ENGINE SETUP (The Solver) ---
N_phys = 60  # Resolution for physics grid
L = 5.0
dx = L / N_phys
dt = 0.05
c = 1.0
gamma = 0.15
lambda_pot = 0.0 # Starts at 0, ramps up (The Quench)

# Physics Grid
x_p = np.linspace(-2.5, 2.5, N_phys)
y_p = np.linspace(-2.5, 2.5, N_phys)
X_P, Y_P = np.meshgrid(x_p, y_p)

# Initial State: Random Vacuum Noise
np.random.seed(137)
psi = np.random.uniform(-np.pi, np.pi, (N_phys, N_phys))
velocity = np.zeros((N_phys, N_phys))

# --- 2. VISUALIZATION SETUP (Your Geometry) ---
fig = plt.figure(figsize=(12, 16), facecolor='black') # Adjusted slightly for video
ax = fig.add_subplot(111, projection='3d', facecolor='black')

# Grid resolution for visual metric warping (Your 'res' variable)
res = 80 # Slightly lower than 120 to make animation render feasible
x_v = np.linspace(-2.5, 2.5, res)
y_v = np.linspace(-2.5, 2.5, res)
X_V, Y_V = np.meshgrid(x_v, y_v)

# --- HELPER FUNCTIONS FROM YOUR CODE ---

def apply_torsion_warp(X_g, Y_g, centers_x, centers_y, strength):
    """Your exact warping logic, adapted for dynamic centers"""
    X_out, Y_out = X_g.copy(), Y_g.copy()
    if len(centers_x) == 0: return X_out, Y_out

    for wx, wy in zip(centers_x, centers_y):
        dx, dy = X_g - wx, Y_g - wy
        dist = np.sqrt(dx**2 + dy**2)
        angle = strength * np.exp(-dist / 0.45)
        X_rot = wx + dx * np.cos(angle) - dy * np.sin(angle)
        Y_rot = wy + dx * np.sin(angle) + dy * np.cos(angle)
        mask = dist < 1.4
        X_out[mask] = X_rot[mask]
        Y_out[mask] = Y_rot[mask]
    return X_out, Y_out

def get_z_geometry(base_h, depth, width, centers_x, centers_y):
    """Your exact Z-geometry logic"""
    well_sum = 0
    if len(centers_x) == 0: return base_h

    for wx, wy in zip(centers_x, centers_y):
        dist = np.sqrt((X_V - wx)**2 + (Y_V - wy)**2)
        well_sum += depth * np.exp(-dist / width)
    return base_h - well_sum

def laplacian(Z):
    """Physics solver helper"""
    Z_top = Z[0:-2, 1:-1]
    Z_left = Z[1:-1, 0:-2]
    Z_bottom = Z[2:, 1:-1]
    Z_right = Z[1:-1, 2:]
    Z_center = Z[1:-1, 1:-1]
    return (Z_top + Z_left + Z_bottom + Z_right - 4 * Z_center) / dx**2

# --- THE MAIN LOOP ---

def update(frame):
    global psi, velocity, lambda_pot

    # === A. PHYSICS STEP ===
    # Ramp coupling strength (The Mega-Snap)
    if lambda_pot < 4.0:
        lambda_pot += 0.05

    # Sine-Gordon Solver
    lap = np.zeros_like(psi)
    lap[1:-1, 1:-1] = laplacian(psi)
    force_topo = -lambda_pot * np.sin(psi)
    force_friction = -gamma * velocity
    accel = (c**2 * lap + force_friction + force_topo)
    velocity += accel * dt
    psi += velocity * dt

    # Boundaries
    psi[0, :] = psi[1, :]
    psi[-1, :] = psi[-2, :]
    psi[:, 0] = psi[:, 1]
    psi[:, -1] = psi[:, -2]

    # === B. DEFECT DETECTION ===
    # Find where the physics engine has created a knot
    grad_x, grad_y = np.gradient(psi)
    stress = np.sqrt(grad_x**2 + grad_y**2)

    # Threshold for a "Well"
    defect_mask = stress > 2.0
    dx_idx, dy_idx = np.where(defect_mask)

    # Extract coordinates for visualizer (limit to 8 strongest to keep visual clean)
    # We take a strided sample to avoid clumping
    step = max(1, len(dx_idx) // 8)
    wells_x = []
    wells_y = []

    if len(dx_idx) > 0:
        # Convert grid indices to physical coordinates
        for i in range(0, len(dx_idx), step):
            ix, iy = dx_idx[i], dy_idx[i]
            # Map physics grid to physical space (-2.5 to 2.5)
            wx = x_p[ix] # Note using physics grid coordinate array
            wy = y_p[iy]
            wells_x.append(wx)
            wells_y.append(wy)

    # === C. RENDERING (Your Exact Style) ===
    ax.clear()
    ax.set_axis_off()
    ax.view_init(elev=22, azim=-45 + frame*0.2) # Slight rotation
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_zlim(0, 10)

    # Title
    ax.text2D(0.5, 0.95, f"Torsional Lagrangian Solver\nCoupling: {lambda_pot:.2f}",
              transform=ax.transAxes, color='white', ha='center', fontsize=18)

    # --- BASE CALCULATIONS ---
    Z1_base = 0.15 * np.sin(X_V) * np.cos(Y_V)

    # Layer 4: Source (Purple/Magenta)
    Z4 = get_z_geometry(7.5 + Z1_base, 1.0, 0.5, wells_x, wells_y)
    X4, Y4 = apply_torsion_warp(X_V, Y_V, wells_x, wells_y, 1.5)
    n4 = plt.Normalize(Z4.min(), Z4.max())
    ax.plot_surface(X4, Y4, Z4, facecolors=cm.cool(n4(Z4)), alpha=0.3, shade=False)
    ax.plot_wireframe(X4, Y4, Z4, color='#ff00ff', alpha=0.1, lw=0.3, rstride=4, cstride=4)

    # Layer 3: Ensemble (Gold)
    Z3 = get_z_geometry(5.0 + Z1_base, 1.4, 0.35, wells_x, wells_y)
    X3, Y3 = apply_torsion_warp(X_V, Y_V, wells_x, wells_y, 3.5)
    n3 = plt.Normalize(Z3.min(), Z3.max())
    ax.plot_surface(X3, Y3, Z3, facecolors=cm.Wistia(n3(Z3)), alpha=0.15, shade=False)
    ax.plot_wireframe(X3, Y3, Z3, color='#fffb00', alpha=0.1, lw=0.3, rstride=4, cstride=4)

    # Layer 2: Gravity (Red/Yellow) - Shear Flow Locus
    Z2 = get_z_geometry(2.5 + Z1_base, 2.0, 0.22, wells_x, wells_y)
    X2, Y2 = apply_torsion_warp(X_V, Y_V, wells_x, wells_y, 9.0)
    n2 = plt.Normalize(Z2.min(), Z2.max())
    ax.plot_surface(X2, Y2, Z2, facecolors=cm.gist_heat(1.0 - n2(Z2)), alpha=0.55, shade=False)
    ax.plot_wireframe(X2, Y2, Z2, color='orange', alpha=0.35, lw=0.5, rstride=3, cstride=3)

    # Layer 1: Substrate (Cyan)
    ax.plot_wireframe(X_V, Y_V, Z1_base, color='#00f2ff', alpha=0.35, lw=0.7, rstride=4, cstride=4)

    # --- CASCADING DNA DRILLS (Dynamic) ---
    for wx, wy in zip(wells_x, wells_y):
        z_target = 0.15 * np.sin(wx) * np.cos(wy)

        # Reduced points slightly for animation speed (1500 -> 300)
        t = np.linspace(0, 1, 300)
        z_path = z_target + (t * 8.5)
        r_drill = 0.12 * (0.35 + 0.65 * t)

        for phase in [0, np.pi]:
            theta = - (20 + 50 * (1-t)) * t + phase + (frame * 0.1) # Added spin rotation
            micro_freq = 250
            micro_r = 0.02

            xf = wx + r_drill * np.cos(theta) + micro_r * np.cos(micro_freq * t)
            yf = wy + r_drill * np.sin(theta) + micro_r * np.sin(micro_freq * t)
            zf = z_path + micro_r * np.sin(micro_freq * t)

            # Simplified color logic for speed
            ax.plot(xf, yf, zf, color='#00f2ff', lw=1.0, alpha=0.5)

        # Core Singularity (White Needle)
        ax.plot([wx, wx], [wy, wy], [z_target, z_target + 8.8], color='white', lw=0.8, alpha=0.5)
        ax.scatter(wx, wy, z_target, color='#00f2ff', s=50, edgecolors='white', alpha=1.0, zorder=30)

    if frame % 5 == 0:
        print(f"Rendering frame {frame}/60...")

# --- EXECUTION ---
print("Initializing Physics-Driven Geometry Solver...")
# 60 frames @ 15fps = 4 seconds of animation
anim = FuncAnimation(fig, update, frames=60, interval=100)

output_file = "torsional_stack_solver.gif"
print(f"Saving to {output_file} (This is computationally heavy, please wait)...")
anim.save(output_file, writer=PillowWriter(fps=15))
print("Done.")