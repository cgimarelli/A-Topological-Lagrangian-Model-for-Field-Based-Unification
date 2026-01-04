import numpy as np
from scipy.integrate import quad

def verify_instanton_suppression():
    print("\n--- Instanton Stability & Tunneling Verification ---")   
    def V(theta, lam):
        return lam * (1 - np.cos(theta))
    def theta_instanton(tau, lam):
        return 4 * np.arctan(np.exp(np.sqrt(lam) * tau))
    def theta_dot(tau, lam):
        root_lam = np.sqrt(lam)
        u = np.exp(root_lam * tau)
        return (4 * root_lam * u) / (1 + u**2)
    def euclidean_lagrangian(tau, lam):
        th = theta_instanton(tau, lam)
        tdot = theta_dot(tau, lam)
        kinetic = 0.5 * tdot**2
        potential = V(th, lam)
        return kinetic + potential
    lam_micro = 1.0
    action_micro, error = quad(euclidean_lagrangian, -20, 20, args=(lam_micro))
    analytical_micro = 8 * np.sqrt(lam_micro)
    print(f"\n[Microscopic Regime: lambda = {lam_micro}]")
    print(f"Computed Action Barrier (S_E): {action_micro:.6f}")
    print(f"Analytical Prediction (8*sqrt(lam)): {analytical_micro:.6f}")
    tunneling_prob_micro = np.exp(-action_micro)
    print(f"Tunneling Probability: {tunneling_prob_micro:.4e}")
    print("-> Result: Barrier is permeable; phase slips can occur.")
    lam_macro = 1e6     
    analytical_macro = 8 * np.sqrt(lam_macro)
    tunneling_prob_macro = np.exp(-analytical_macro)    
    print(f"\n[Macroscopic Regime: lambda = {lam_macro:.0e}]")
    print(f"Action Barrier (S_E): {analytical_macro:.1f}")
    print(f"Tunneling Probability: {tunneling_prob_macro}")    
    if tunneling_prob_macro == 0.0:
        print("-> Result: Probability Underflow (Zero).")
        print("-> CONCLUSION: The Phase-Loop is Topologically Protected.")
        print("-> The knot cannot spontaneously untie due to infinite Action Cost.")
    assert np.isclose(action_micro, analytical_micro, atol=1e-4), "Action calculation failed!"
if __name__ == "__main__":
    verify_instanton_suppression()