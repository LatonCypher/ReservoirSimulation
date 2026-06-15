using MKLNET;
using System.Collections.Generic;
using static System.Math;

namespace ReservoirSimulator
{
    public class WolfeNewtonSolver
    {
        // Delegate matching your forward-mode AutoDiff signature
        internal delegate void SystemEvaluator(double[] x, double[] residual, Jacobian jacobian);

        // High-performance CSR Matrix Vector multiplication tool from your engine
        // Custom signature: public double[] Mult(List<double> a_value, List<int> a_index, List<int> a_start, double[] x)

        internal static (bool, int) Solve(
            SystemEvaluator evaluator,
            double[] xInitial,
            Func<Jacobian, double[], double[]> linearSolver, // Function wrapper for J * p = -R
            double tolx = 1e-7,
            double tolf = 1e-6,
            int maxNewtonIterations = 10,
            int maxLineSearchAttempts = 10)
        {
            int n = xInitial.Length;
            double[] x = (double[])xInitial.Clone();

            // Memory allocations
            double[] xOld = new double[n];
            double[] xTrial = new double[n];
            double[] residual = new double[n];
            Jacobian jacobian = new();

            double[] trialResidual = new double[n];
            Jacobian trialJacobian = new();

            // Wolfe Condition Hyperparameters
            double c1 = 1e-4; // Sufficient decrease parameter (Armijo)
            double c2 = 0.9;  // Curvature condition parameter (for Newton methods, 0.9 is standard)

            // Initial Evaluation
            evaluator(x, residual, jacobian);


            Console.WriteLine($"""
                    iter  |   Residual Norm  
                ----------+----------------
                """);
            for (int iter = 1; iter <= maxNewtonIterations; iter++)
            {
                // 1. Check convergence based on the residual infinity norm
                double maxRes = residual.Max(Abs);
                Console.WriteLine($"  {iter,4}    |    {maxRes:F4}");
                if (maxRes < tolf)
                {
                    Array.Copy(x, xInitial, n); // Output final answer in-place
                    return (true, iter);
                }

                // 2. Compute the Newton update direction p by solving: J * p = -R
                double[] p = linearSolver(jacobian, residual);

                // 3. Compute baseline objective function and initial slope
                double fOld = 0.5 * ComputeNormSquared(residual);

                // Initial slope shortcut: grad(f)^T * p = -||R||^2 when using exact Newton step
                double slopeStart = -2*fOld;
                if (slopeStart >= 0.0)
                {
                    Console.WriteLine("[Error] Non-descent direction generated. Check Jacobian properties.");
                    return (false, iter);
                }

                // Save previous state layout
                Array.Copy(x, xOld, n);

                // 4. Line Search with Derivative Tracking via Strong Wolfe Conditions
                double alpha = 1.0; // Always attempt full Newton step first
                double alphaOld = 0.0;
                double fOldLs = fOld;
                double slopeOldLs = slopeStart;

                bool stepAccepted = false;

                for (int lsAttempt = 0; lsAttempt < maxLineSearchAttempts; lsAttempt++)
                {
                    // Advance trial state: xTrial = xOld + alpha * p
                    for (int i = 0; i < n; i++) xTrial[i] = xOld[i] + alpha * p[i];

                    // Evaluate BOTH residual and Jacobian (forced by Forward AutoDiff)
                    evaluator(xTrial, trialResidual, trialJacobian);

                    double fTrial = 0.5 * ComputeNormSquared(trialResidual);

                    // Compute the fresh slope at the trial position using the new Jacobian:
                    // slopeTrial = trialResidual^T * (trialJacobian * p)
                    double[] Jp = MultiplySparseMatrixByVector(trialJacobian, p);
                    double slopeTrial = DotProduct(trialResidual, Jp);

                    // Condition A: Sufficient Decrease (Armijo Condition)
                    if (fTrial > fOld + c1 * alpha * slopeStart || (lsAttempt > 0 && fTrial >= fOldLs))
                    {
                        // Overshot or error increased: Interpolate using cubic curve between alphaOld and alpha
                        alpha = CubicHermiteInterpolate(alphaOld, fOldLs, slopeOldLs, alpha, fTrial, slopeTrial);
                        // Update bounds tracking
                        fOldLs = fTrial; slopeOldLs = slopeTrial; alphaOld = alpha;
                        continue;
                    }

                    // Condition B: Strong Curvature Condition
                    if (Abs(slopeTrial) <= c2 * Abs(slopeStart))
                    {
                        // Step matches safe gradient flatness bounds! Accept it.
                        stepAccepted = true;
                        fOld = fTrial;
                        Array.Copy(trialResidual, residual, n);
                        jacobian = trialJacobian.Duplicate();
                        Array.Copy(xTrial, x, n);
                        break;
                    }

                    // If it satisfies decrease but the slope is still too steeply negative,
                    // we haven't stepped far enough into the valley. 
                    if (slopeTrial >= 0.0)
                    {
                        // Passed the local minimum, backtrack using interpolation
                        alpha = CubicHermiteInterpolate(alpha, fTrial, slopeTrial, alphaOld, fOldLs, slopeOldLs);
                        stepAccepted = true; // Force exit to evaluate next global step layout
                        break;
                    }

                    // Otherwise, safely extrapolate/step forward
                    alphaOld = alpha;
                    fOldLs = fTrial;
                    slopeOldLs = slopeTrial;
                    alpha = Min(2.0 * alpha, 1.0); // Restrict forward bounds to 1.0max
                }

                if(double.IsNaN(alpha) || alpha <= 1e-12)
                {

                    Console.WriteLine($"  {iter,4}    |    {maxRes:F4}");
                    Console.WriteLine("[Warning] Line search produced non-positive or NaN alpha.");
                    return (false, iter);
                }

                if (!stepAccepted)
                {
                    Console.WriteLine($"""
                                  |    - [Warning] Line search reached maximum limits. 
                                  |                Accepting best found alpha = {alpha:F4}
                        """);
                        
                    Array.Copy(xTrial, x, n);
                    Array.Copy(trialResidual, residual, n);
                    jacobian = trialJacobian.Duplicate();
                }

                // 5. Test state convergence bound (TolX)
                //double deltaXNorm = 0.0;
                //for (int i = 0; i < n; i++)
                //{
                //    double temp = Abs(x[i] - xOld[i]) / Max(Abs(x[i]), 1.0);
                //    deltaXNorm = Max(deltaXNorm, temp);
                //}
                //if (deltaXNorm < tolx)
                //{
                //    maxRes = residual.Max(Abs);
                //    Console.WriteLine($"  {iter,4}*   |    {maxRes:F4}");
                //    Array.Copy(x, xInitial, n);
                //    return (true, iter);
                //}
            }
            Console.WriteLine("[Warning] Solver exited without meeting convergence criteria.");
            return (false, maxNewtonIterations);
        }

        /// <summary>
        /// Performs Cubic Hermite Interpolation to minimize the step length dynamically
        /// using values and slopes evaluated at both boundary parameters.
        /// </summary>
        static double CubicHermiteInterpolate(double a0, double f0, double d0, double a1, double f1, double d1)
        {
            double da = a1 - a0;
            if (Abs(da) < 1e-12) return 0.5 * (a0 + a1);

            double d1_val = d0 + d1 - 3.0 * (f0 - f1) / da;
            double d2_val = Sign(da) * Sqrt(d1_val * d1_val - d0 * d1);

            double numerator = d1 - d0 + 2.0 * d2_val;
            if (Abs(numerator) < 1e-12) return 0.5 * (a0 + a1);

            double alphaNew = a1 - da * ((d1 + d2_val - d1_val) / numerator);

            // Safeguard safeguarding boundaries to prevent extreme compressions/extrapolations
            double minBound = a0 + 0.1 * da;
            double maxBound = a0 + 0.9 * da;
            return Clamp(alphaNew, Min(minBound, maxBound), Max(minBound, maxBound));
        }

        static double ComputeNormSquared(double[] vec)
        {
            double sum = 0.0;
            for (int i = 0; i < vec.Length; i++) sum += vec[i] * vec[i];
            return sum;
        }

        static double DotProduct(double[] v1, double[] v2)
        {
            double sum = 0.0;
            for (int i = 0; i < v1.Length; i++) sum += v1[i] * v2[i];
            return sum;
        }

        static double[] MultiplySparseMatrixByVector(Jacobian jacobian, double[] vec)=>
            Mult(jacobian.Value, jacobian.Index, jacobian.Start, vec);

        static double[] Mult(List<double> a_value, List<int> a_index, List<int> a_start, double[] x)
        {
            // The number of rows is determined by the size of the a_start array minus 1
            int numRows = a_start.Count - 1;

            // Allocate the output column vector y
            double[] y = new double[numRows];

            // Iterate through each row of the sparse matrix
            for (int i = 0; i < numRows; i++)
            {
                // Identify the bounds of the current row within the values list
                int rowStart = a_start[i];
                int rowEnd = a_start[i + 1];

                double sum = 0.0;

                // Perform the dot product of the sparse row elements and the vector x
                for (int k = rowStart; k < rowEnd; k++)
                {
                    int col = a_index[k];
                    sum += a_value[k] * x[col];
                }

                // Assign the accumulated result to the output vector
                y[i] = sum;
            }

            return y;
        }
    }
}