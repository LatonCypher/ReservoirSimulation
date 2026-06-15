using static System.Math;
namespace ReservoirSimulator
{
    public class GlobalNewtonSolver
    {
        // A delegate representing the system's physics/equations.
        // Given the current state vector x, it fills and returns the Residual vector and Jacobian matrix.
        public delegate void SystemEvaluator(double[] x, double[] residual, double[,] jacobian);

        private const int MAXITS = 200;
        private const double ALF = 1.0e-4;
        private const double EPS = 1.11e-16; // Machine precision (equivalent to double epsilon)

        /// <summary>
        /// Globally Convergent Multi-dimensional Newton-Raphson Solver from Numerical Recipes.
        /// </summary>
        /// <param name="x">Initial guess vector (updated in-place with the solution).</param>
        /// <param name="evaluator">User function providing residuals and Jacobian matrix coefficients.</param>
        /// <param name="tolx">Convergence criterion on structural state changes.</param>
        /// <param name="tolf">Convergence criterion on residual function values.</param>
        /// <returns>True if convergence succeeded; False if it trapped at a local minimum of f.</returns>
        public static bool Newt(double[] x, SystemEvaluator evaluator, double tolx = 1e-7, double tolf = 1e-6)
        {
            int n = x.Length;

            // Memory Allocations
            double[] g = new double[n];          // Gradient of f = 0.5 * F * F
            double[] p = new double[n];          // Update/search direction step vector (dx)
            double[] xold = new double[n];       // Retained previous state vector
            double[] residual = new double[n];   // Active equation residual vector F(x)
            double[,] jacobian = new double[n, n]; // Active System Jacobian Matrix J(x)

            // Evaluate initial system state
            evaluator(x, residual, jacobian);

            // Compute initial objective function value: f = 0.5 * sum(F_i^2)
            double f = 0.0;
            for (int i = 0; i < n; i++) f += residual[i] * residual[i];
            f *= 0.5;

            // Calculate stpmax: Maximum allowable step length for line search
            double sum = 0.0;
            for (int i = 0; i < n; i++) sum += x[i] * x[i];
            double stpmax = 100.0 * Max(Sqrt(sum), (double)n);

            for (int its = 1; its <= MAXITS; its++)
            {
                // 1. Compute Gradient vector g = J^T * F
                // And check if the current state vector is already fully converged
                double test = 0.0;
                for (int i = 0; i < n; i++)
                {
                    double gradSum = 0.0;
                    for (int j = 0; j < n; j++)
                    {
                        gradSum += jacobian[j, i] * residual[j];
                    }
                    g[i] = gradSum;

                    double den = Max(Abs(x[i]), 1.0);
                    test = Max(test, Abs(residual[i]) * den); // Normalised residual check
                }

                if (test < tolf) return true; // Converged on residual bounds!

                // 2. Set up the Newton direction right-hand side vector: p = -Residual
                for (int i = 0; i < n; i++) p[i] = -residual[i];

                // 3. Solve the linear system: J * p = -R via your preferred solver
                // (e.g., Gaussian Elimination, LU, or your native SepalSolver engine)
                SolveLinearSystem(jacobian, p);

                // Copy current state to memory buffers before launching line search
                double fold = f;
                Array.Copy(x, xold, n);

                // 4. Execute the Backtracking Line Search routine
                bool check = LineSearch(xold, fold, g, p, x, out f, stpmax, evaluator, residual);

                // 5. Test state-change convergence bounds (TolX)
                test = 0.0;
                for (int i = 0; i < n; i++)
                {
                    double temp = Abs(x[i] - xold[i]) / Max(Abs(x[i]), 1.0);
                    if (temp > test) test = temp;
                }

                if (test < tolx)
                {
                    // If line search flagged a boundary dead-end AND state variables aren't moving,
                    // the solver has stalled on a local minimum of f instead of a true root.
                    if (check)
                    {
                        Console.WriteLine("[Warning] Stalled at a local minimum of the residual function.");
                        return false;
                    }
                    return true; // Clean convergence on state changes
                }

                // 6. Re-evaluate the full system Jacobian for the next step cycle
                evaluator(x, residual, jacobian);
            }

            Console.WriteLine("[Warning] Maximum iterations exceeded in globally convergent Newton loop.");
            return false;
        }

        /// <summary>
        /// Backtracking line search routine from Section 9.7.
        /// </summary>
        private static bool LineSearch(double[] xold, double fold, double[] g, double[] p,
                                       double[] x, out double f, double stpmax,
                                       SystemEvaluator evaluator, double[] trialResidual)
        {
            int n = xold.Length;
            bool check = false;

            // Verify if step length exceeds maximum constraint, scale back if necessary
            double sum = 0.0;
            for (int i = 0; i < n; i++) sum += p[i] * p[i];
            sum = Math.Sqrt(sum);

            if (sum > stpmax)
            {
                double scale = stpmax / sum;
                for (int i = 0; i < n; i++) p[i] *= scale;
            }

            // Compute the initial directional derivative slope: g^T * p
            double slope = 0.0;
            for (int i = 0; i < n; i++) slope += g[i] * p[i];
            if (slope >= 0.0)
                throw new InvalidOperationException("Roundoff error or non-descent direction in LineSearch.");

            // Compute minimum step allocation scale (alamin)
            double test = 0.0;
            for (int i = 0; i < n; i++)
            {
                double temp = Abs(p[i]) / Max(Abs(xold[i]), 1.0);
                if (temp > test) test = temp;
            }
            double alamin = tolx_internal(test);

            // Line Search Backtracking Loop
            double alam = 1.0; // Start with full Newton step length
            double alam2 = 0.0;
            double f2 = 0.0;
            double[,] dummyJ = new double[n, n]; // Pre-allocated throwaway array for evaluation signatures

            while (true)
            {
                // Advance the trial state vector: x = xold + alam * p
                for (int i = 0; i < n; i++) x[i] = xold[i] + alam * p[i];

                // Evaluate the updated residual state (Jacobian values are ignored during line search)
                evaluator(x, trialResidual, dummyJ);

                // Compute candidate objective function value
                f = 0.0;
                for (int i = 0; i < n; i++) f += trialResidual[i] * trialResidual[i];
                f *= 0.5;

                if (alam < alamin)
                {
                    // Step size shrunk past resolution limits without finding a decrease
                    Array.Copy(xold, x, n);
                    check = true;
                    return check;
                }
                else if (f <= fold + ALF * alam * slope)
                {
                    // Sufficient decrease verified (Armijo / Alpha Condition met)!
                    return check;
                }
                else
                {
                    // Backtrack: Model f(alam) using a polynomial curve fit to find the optimal next step
                    double tmplam;
                    if (alam == 1.0)
                    {
                        // First backtrack utilizes a simple quadratic fit
                        tmplam = -slope / (2.0 * (f - fold - slope));
                    }
                    else
                    {
                        // Subsequent backtracks utilize a cubic fit using historical step parameters
                        double rhs1 = f - fold - alam * slope;
                        double rhs2 = f2 - fold - alam2 * slope;
                        double a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
                        double b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);

                        if (a == 0.0)
                        {
                            tmplam = -slope / (2.0 * b);
                        }
                        else
                        {
                            double disc = b * b - 3.0 * a * slope;
                            if (disc < 0.0) tmplam = 0.5 * alam;
                            else if (b <= 0.0) tmplam = (-b + Sqrt(disc)) / (3.0 * a);
                            else tmplam = -slope / (b +     Sqrt(disc));
                        }
                        tmplam =    Min(tmplam, 0.5 * alam); // Bound step length reduction to at most 0.5
                    }

                    alam2 = alam;
                    f2 = f;
                    // Enforce a minimum step reduction bound (at least 0.1 * alam) to prevent infinitely small crawls
                    alam = Max(tmplam, 0.1 * alam);
                }
            }
        }

        private static double tolx_internal(double test) => test > 0.0 ? EPS / test : 0.0;

        /// <summary>
        /// Standard Gaussian Elimination Solver for J * p = -R.
        /// </summary>
        private static void SolveLinearSystem(double[,] J, double[] p)
        {
            int n = p.Length;
            for (int i = 0; i < n; i++)
            {
                int pivot = i;
                for (int k = i + 1; k < n; k++)
                    if (Abs(J[k, i]) > Abs(J[pivot, i])) pivot = k;

                for (int j = i; j < n; j++)
                {
                    (J[pivot, j], J[i, j])=(J[i, j], J[pivot, j]);
                }
                (p[pivot], p[i])=(p[i], p[pivot]);
                if (Abs(J[i, i]) < 1e-15)
                    throw new InvalidOperationException("Jacobian is structurally singular.");

                for (int k = i + 1; k < n; k++)
                {
                    double factor = J[k, i] / J[i, i];
                    p[k] -= factor * p[i];
                    for (int j = i; j < n; j++) J[k, j] -= factor * J[i, j];
                }
            }

            for (int i = n - 1; i >= 0; i--)
            {
                for (int j = i + 1; j < n; j++) p[i] -= J[i, j] * p[j];
                p[i] /= J[i, i];
            }
        }
    }
}