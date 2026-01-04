import unittest; 
import numpy as np; 
import matplotlib.pyplot as plt
class TestCausalityConsistency(unittest.TestCase):
    def test_ctc_selection_mechanism(self):
        print("\n--- Testing CTC Causal Selection Mechanism ---"); tau_0 = 1.0; np.random.seed(55)
        num_histories = 1000; omegas = np.random.uniform(0.5, 10.5, num_histories) * (2 * np.pi)
        consistent_count = 0; paradox_count = 0; consistency_errors = []; threshold = 0.05
        for omega in omegas:
            delta_theta = omega * tau_0; consistency_error = 1.0 - np.cos(delta_theta)
            consistency_errors.append(consistency_error)
            if consistency_error < threshold: consistent_count += 1
            else: paradox_count += 1
        survival_rate = consistent_count / num_histories
        print(f"Total: {num_histories} | Paradoxical: {paradox_count} | Consistent: {consistent_count}")
        print(f"Causal Selection Rate: {survival_rate*100:.2f}%")
        self.assertTrue(survival_rate < 0.2); self.assertTrue(survival_rate > 0.0)
        print("SUCCESS: Cyclic topology filters paradoxical histories.")
        survivor_omegas = [w for w, err in zip(omegas, consistency_errors) if err < threshold]
        for w in survivor_omegas[:5]:
            n = (w * tau_0) / (2 * np.pi); print(f" Survivor Mode n = {n:.4f}"); self.assertAlmostEqual(n, round(n), delta=0.1)
if __name__ == '__main__': unittest.main(argv=['first-arg-is-ignored'], exit=False)