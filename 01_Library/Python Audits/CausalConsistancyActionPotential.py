import numpy as np
def verify_causal_filtering_action():
    print("\n--- Causal Consistency Test (Action-Potential Model) ---")
    tau_0 = 1.0 
    num_histories = 100000 
    vacuum_stiffness = 50.0 
    np.random.seed(42)
    omegas = np.random.uniform(0.5, 10.5, num_histories) * (2 * np.pi)
    delta_thetas = omegas * tau_0
    potential = 1.0 - np.cos(delta_thetas)
    weights = np.exp(-vacuum_stiffness * potential)
    winding_numbers = delta_thetas / (2 * np.pi)
    nearest_integers = np.round(winding_numbers)
    deviations = np.abs(winding_numbers - nearest_integers)    
    weighted_mean_deviation = np.average(deviations, weights=weights)    
    print(f"Total Histories: {num_histories}")
    print(f"Vacuum Stiffness (Lambda): {vacuum_stiffness}")
    print(f"Random Noise Deviation: 0.250000")
    print(f"Action-Filtered Deviation: {weighted_mean_deviation:.6f}")    
    if weighted_mean_deviation < 0.05:
        print("-> SUCCESS: The topological potential forces strict quantization.")
        print("-> Causal paradoxes are physically suppressed by the action barrier.")
    else:
        print("-> FAIL: Vacuum too soft.")
if __name__ == "__main__":
    verify_causal_filtering_action()