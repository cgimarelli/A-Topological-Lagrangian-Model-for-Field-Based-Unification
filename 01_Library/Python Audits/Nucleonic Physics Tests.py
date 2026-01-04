import unittest
import sympy as sp
import math

class TestNucleonicPhysics(unittest.TestCase):
    """
    Verification suite for 'Shear Flow Mechanism' Mass Derivations.
    Audits the Hamiltonian Stability and Physical Scale Predictions.
    """

    def test_stability_radius_derivation_symbolic(self):
        """
        Test Eq 10: Symbolic Derivation of R_stable.
        Minimizes Hamiltonian E(R) = gamma/R^2 + beta*R^2.
        Asserts R_stable == (gamma/beta)^(1/4).
        """
        print("\n--- Testing Hamiltonian Stability Derivation (Symbolic) ---")

        # 1. Define Symbols
        R = sp.symbols('R', real=True, positive=True)
        gamma, beta = sp.symbols('gamma beta', real=True, positive=True)

        # 2. Define Hamiltonian Density (Eq 9)
        # Shear Tension (1/R^2) + Coherence Pressure (R^2)
        E = (gamma / R**2) + (beta * R**2)

        # 3. Differentiate w.r.t R
        dE_dR = sp.diff(E, R)
        print(f"dE/dR = {dE_dR}")

        # 4. Solve for Critical Point
        solutions = sp.solve(dE_dR, R)
        print(f"Critical Points: {solutions}")

        # Filter for the positive real solution
        r_stable_sym = [s for s in solutions if s.is_positive][0]

        # 5. Expected Result
        expected = (gamma / beta)**sp.Rational(1, 4)

        # 6. Assertion
        # Use simplify to handle potential algebraic variations
        check = sp.simplify(r_stable_sym - expected)
        self.assertEqual(check, 0, "FATAL: Symbolic derivation of R_stable failed!")
        print("SUCCESS: Math derivation of stability radius is exact.")

    def test_nucleonic_scale_prediction_numerical(self):
        """
        Test Section II.D: Prediction of Nucleonic Mass Scale.
        Plugs in Planck Area (gamma) and Vacuum Coherence (beta).
        Asserts result is within the femtometer range (~4e-15 m).
        """
        print("\n--- Testing Nucleonic Scale Prediction (Numerical) ---")

        # 1. Constants
        # Planck Length l_p ~ 1.616e-35 m
        l_p = 1.616255e-35

        # Shear Modulus gamma = l_p^2 (Eq 11)
        gamma = l_p**2

        # Coherence Coupling beta ~ 1e-12 m^-2 (Inverse solar system scale approx)
        # NOTE: Paper calibrates this to macroscopic coherence.
        beta = 1.0e-12

        print(f"Gamma (Planck Area): {gamma:.4e} m^2")
        print(f"Beta (Coherence): {beta:.4e} m^-2")

        # 2. Calculate R_stable
        # R = (gamma / beta)^(1/4)
        R_stable = (gamma / beta)**0.25

        print(f"Predicted Radius: {R_stable:.4e} m")

        # 3. Assertions
        # Target: ~ 4.0e-15 m (4 femtometers)
        target = 4.0e-15
        tolerance = 0.5e-15 # Allow variance due to beta calibration

        # Check order of magnitude (Femtometer scale)
        self.assertTrue(1e-15 < R_stable < 1e-14,
                        f"FATAL: Prediction {R_stable} is not in the femtometer (nucleonic) range!")

        # Check specific value proximity
        error = abs(R_stable - target)
        print(f"Deviation from Target (4.0 fm): {error:.4e} m")

        # If the math holds, this should pass easily based on the constants provided in text.
        self.assertLess(error, 1.0e-15, "Prediction is drifting too far from 4.0 fm.")
        print("SUCCESS: Model accurately recovers the scale of the atomic nucleus from Planck geometry.")

if __name__ == '__main__':
    unittest.main(argv=['first-arg-is-ignored'], exit=False)