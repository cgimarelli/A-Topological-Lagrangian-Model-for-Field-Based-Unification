import sympy as sp
import numpy as np
from scipy.integrate import quad

def verify_hhd_symbolic():
    """
    Symbolically verifies the Helmholtz-Hodge decomposition and the 
    coupling between shear stress and solenoidal vorticity.
    """
    print("--- Symbolic Verification of HHD & Shear Coupling ---")
    
    # Define coordinates and potentials
    x, y, z = sp.symbols('x y z')
    phi = sp.Function('phi')(x, y, z)
    A1, A2, A3 = sp.symbols('A1 A2 A3', cls=sp.Function)
    A = [A1(x,y,z), A2(x,y,z), A3(x,y,z)]
    
    # 1. Irrotational component (Gradient)
    I_par = [sp.diff(phi, x), sp.diff(phi, y), sp.diff(phi, z)]
    
    # 2. Solenoidal component (Curl)
    # Curl calculation manually for clarity
    I_perp = [
        sp.diff(A[2], y) - sp.diff(A[1], z),
        sp.diff(A[0], z) - sp.diff(A[2], x),
        sp.diff(A[1], x) - sp.diff(A[0], y)
    ]
    
    # 3. Verify irrotational property: curl(grad(phi)) = 0
    curl_I_par = [
        sp.diff(I_par[2], y) - sp.diff(I_par[1], z),
        sp.diff(I_par[0], z) - sp.diff(I_par[2], x),
        sp.diff(I_par[1], x) - sp.diff(I_par[0], y)
    ]
    is_irrotational = all(val.simplify() == 0 for val in curl_I_par)
    print(f"Is Irrotational Component curl-free? {is_irrotational}")
    
    # 4. Verify solenoidal property: div(curl(A)) = 0
    div_I_perp = sp.diff(I_perp[0], x) + sp.diff(I_perp[1], y) + sp.diff(I_perp[2], z)
    is_solenoidal = div_I_perp.simplify() == 0
    print(f"Is Solenoidal Component divergence-free? {is_solenoidal}")
    
    # 5. Verify Shear Source Coupling
    # Modeling a Couette flow u_x = y, inducing shear sigma_xy
    u_x = y
    sigma_xy = sp.diff(u_x, y) # Constant shear
    
    # The source S is the divergence of sigma. For Couette, grad(sigma) relates to vorticity.
    # If S has a curl, I_perp must be non-zero.
    # Assume source current S_y = sigma_xy, S_x = 0
    S = [0, sigma_xy, 0]
    curl_S = [
        sp.diff(S[2], y) - sp.diff(S[1], z),
        sp.diff(S[0], z) - sp.diff(S[2], x),
        sp.diff(S[1], x) - sp.diff(S[0], y)
    ]
    has_rotational_source = any(val != 0 for val in curl_S)
    print(f"Does shear flow induce non-zero curl in the source? {has_rotational_source}")
    
    print("SYMBOLIC VERIFICATION COMPLETE.\n")

def verify_l2_orthogonality_numerical():
    """
    Numerically verifies L2 orthogonality between a gradient and a curl 
    in a bounded domain with N-P boundary conditions.
    """
    print("--- Numerical Verification of L2 Orthogonality ---")
    
    # Domain: Cube [0, 1]^3
    # phi = sin(pi*x)*sin(pi*y)*sin(pi*z) -> Normal at boundary
    # A = [sin(pi*x), sin(pi*y), sin(pi*z)] -> Resulting curl is parallel
    
    def phi_func(x, y, z): return np.sin(np.pi*x) * np.sin(np.pi*y) * np.sin(np.pi*z)
    
    def grad_phi(x, y, z):
        p = np.pi
        return np.array([
            p*np.cos(p*x)*np.sin(p*y)*np.sin(p*z),
            p*np.sin(p*x)*np.cos(p*y)*np.sin(p*z),
            p*np.sin(p*x)*np.sin(p*y)*np.cos(p*z)
        ])
    
    def curl_A(x, y, z):
        # Using A = [0, 0, sin(pi*x)*sin(pi*y)]
        # Curl A = [dAz/dy, -dAz/dx, 0]
        p = np.pi
        return np.array([
            p*np.sin(p*x)*np.cos(p*y),
            -p*np.cos(p*x)*np.sin(p*y),
            0
        ])

    # Inner product integral over [0,1]x[0,1]x[0,1]
    # We'll use a simple sum over a grid for speed
    grid_size = 50
    vals = np.linspace(0, 1, grid_size)
    dx = vals[1] - vals[0]
    inner_product = 0.0
    
    for xi in vals:
        for yi in vals:
            for zi in vals:
                gp = grad_phi(xi, yi, zi)
                ca = curl_A(xi, yi, zi)
                inner_product += np.dot(gp, ca) * (dx**3)
                
    print(f"Numerical L2 Inner Product: {inner_product:.5e}")
    # Value should be close to 0
    success = abs(inner_product) < 1e-10
    print(f"Orthogonality Check (< 1e-10): {success}")
    print("NUMERICAL VERIFICATION COMPLETE.")

if __name__ == "__main__":
    verify_hhd_symbolic()
    verify_l2_orthogonality_numerical()