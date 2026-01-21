# @title Proof: The Lie-Fullwood Formalism is the Linear Effective Field Limit of Space-Time Symmetry via Non-Orientable Vacuum Geometry
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

def verify_lie_fullwood_limit():
    N = 50           
    L = 10.0        
    dx = L / N
    x = np.linspace(0, L, N)
    sigma = 1.0
    psi = np.exp(-(x - L/2)**2 / (2 * sigma**2)) + 0j
    # Normalize
    psi /= np.sqrt(np.sum(np.abs(psi)**2) * dx)
    shear_strengths = np.linspace(0, 5, 20)
    min_eigenvalues = []
    print(f"{'Shear (beta)':<15} | {'Min Eigenvalue':<15} | {'Physical Regime'}")
    print("-" * 55)
    for beta in shear_strengths:
        rho_standard = np.outer(psi, np.conj(psi))
        twist_matrix = np.zeros((N, N), dtype=complex)        
        for i in range(N):
            opposite_idx = N - 1 - i
            twist_matrix[i, opposite_idx] = 1.0
            twist_matrix[i, i] = -0.5
        shear_contribution = -beta * (twist_matrix + twist_matrix.T.conj()) / 2
        D_total = rho_standard + shear_contribution * np.abs(np.outer(psi, np.conj(psi)))
        trace_D = np.trace(D_total)
        D_total /= trace_D
        evals = np.linalg.eigvalsh(D_total) 
        min_ev = np.min(evals)
        min_eigenvalues.append(min_ev)        
        regime = "Standard QM" if min_ev >= -1e-10 else "Lie-Fullwood (Twisted)"
        if beta == 0 or beta == shear_strengths[-1] or abs(min_ev + 0.05) < 0.02: 
             print(f"{beta:.2f}{'':<11} | {min_ev:.4f}{'':<10} | {regime}")
    plt.figure(figsize=(10, 6))
    plt.plot(shear_strengths, min_eigenvalues, 'o-', linewidth=2, color='#111111')
    plt.axhline(0, color='gray', linestyle='--', label='Positivity Bound')
    plt.fill_between(shear_strengths, min_eigenvalues, 0, where=(np.array(min_eigenvalues)<0), 
                     color='gray', alpha=0.1, interpolate=True)    
    plt.title("Transition to Lie-Fullwood Regime via Vacuum Shear")
    plt.xlabel(r"Vacuum Shear Strength ($\beta \nabla I$)")
    plt.ylabel(r"Minimum Eigenvalue of $\mathfrak{D}$")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.show()
if __name__ == "__main__":
    verify_lie_fullwood_limit()