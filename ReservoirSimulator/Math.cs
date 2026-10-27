using System.Runtime.CompilerServices;
using SepalSolver;

namespace ReservoirSimulator
{
    internal static class Math
    {
        private const double Epsilon = 1e-14;
        internal const double pi = SepalSolver.Math.pi;
        internal static int DerivativeCapacity = 31;

        // --- Scalar Wrappers ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Floor(double x) => SepalSolver.Math.Floor(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double[] Linspace(double x, double y, int N) => SepalSolver.Math.Linspace(x, y, N);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Sqr(double x) => SepalSolver.Math.Sqr(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Sin(double x) => SepalSolver.Math.Sin(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Cos(double x) => SepalSolver.Math.Cos(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Tan(double x) => SepalSolver.Math.Tan(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Asin(double x) => SepalSolver.Math.Asin(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Acos(double x) => SepalSolver.Math.Acos(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Atan(double x) => SepalSolver.Math.Atan(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Sinh(double x) => SepalSolver.Math.Sinh(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Cosh(double x) => SepalSolver.Math.Cosh(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Tanh(double x) => SepalSolver.Math.Tanh(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Asinh(double x) => SepalSolver.Math.Asinh(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Acosh(double x) => SepalSolver.Math.Acosh(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Atanh(double x) => SepalSolver.Math.Atanh(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Sqrt(double x) => SepalSolver.Math.Sqrt(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Cbrt(double x) => SepalSolver.Math.Pow(x, 1.0 / 3.0);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Log(double x) => SepalSolver.Math.Log(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Exp(double x) => SepalSolver.Math.Exp(x);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Hypot(double x, double y) => SepalSolver.Math.Hypot(x, y);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Min(double x, double y) => x <= y ? x : y;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Max(double x, double y) => x >= y ? x : y;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static int Min(int x, int y) => x <= y ? x : y;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static int Max(int x, int y) => x >= y ? x : y;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static int Sign(double x) => (x > 0) ? 1 : ((x < 0) ? -1 : 0);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Pow(double x, double y) => SepalSolver.Math.Pow(x, y);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Abs(double x) => x >= 0 ? x : -x;

        // --- ADiff<T> Overloads ---

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Abs<T>(ADiff<T> a) where T : struct, IVectorND<T>
            => (a.Value >= 0.0) ? a : -a;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Max<T>(ADiff<T> a, ADiff<T> b) where T : struct, IVectorND<T>
            => (a.Value >= b.Value) ? a : b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Min<T>(ADiff<T> a, ADiff<T> b) where T : struct, IVectorND<T>
            => (a.Value <= b.Value) ? a : b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Pow<T>(ADiff<T> a, double b) where T : struct, IVectorND<T>
        {
            double powerVal = Pow(a.Value, b);
            // Protect against 0 * NaN or division by zero when a.Value == 0 and b > 1 (e.g. relative permeability curves)
            double derivativeScale;
            if (System.Math.Abs(a.Value) < 1e-30)
                derivativeScale = (b > 1.0) ? 0.0 : ((b == 1.0) ? 1.0 : double.PositiveInfinity);
            else
                derivativeScale = b * Pow(a.Value, b - 1.0);
            return new ADiff<T>(powerVal, a.Derivative * derivativeScale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Exp<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double expVal = Exp(a.Value);
            return new ADiff<T>(expVal, a.Derivative * expVal);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Log<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double logVal = Log(a.Value);
            double invVal = 1.0 / a.Value;
            return new ADiff<T>(logVal, a.Derivative * invVal);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Sqrt<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double sqrtVal = Sqrt(a.Value);
            double inv2Sqrt = 0.5 / sqrtVal;
            return new ADiff<T>(sqrtVal, a.Derivative * inv2Sqrt);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Cbrt<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double cbrtVal = Cbrt(a.Value);
            double derivativeScale = 1.0 / (3.0 * cbrtVal * cbrtVal);
            return new ADiff<T>(cbrtVal, a.Derivative * derivativeScale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Sin<T>(ADiff<T> a) where T : struct, IVectorND<T>
            => new ADiff<T>(Sin(a.Value), a.Derivative * Cos(a.Value));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Cos<T>(ADiff<T> a) where T : struct, IVectorND<T>
            => new ADiff<T>(Cos(a.Value), -a.Derivative * Sin(a.Value));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Tan<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double tanVal = Tan(a.Value);
            double cosVal = Cos(a.Value);
            double derivativeScale = 1.0 / (cosVal * cosVal);
            return new ADiff<T>(tanVal, a.Derivative * derivativeScale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Asin<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double scale = 1.0 / Sqrt(1.0 - a.Value * a.Value);
            return new ADiff<T>(Asin(a.Value), a.Derivative * scale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Acos<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double scale = -1.0 / Sqrt(1.0 - a.Value * a.Value);
            return new ADiff<T>(Acos(a.Value), a.Derivative * scale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Atan<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double scale = 1.0 / (1.0 + a.Value * a.Value);
            return new ADiff<T>(Atan(a.Value), a.Derivative * scale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Sinh<T>(ADiff<T> a) where T : struct, IVectorND<T>
            => new ADiff<T>(Sinh(a.Value), a.Derivative * Cosh(a.Value));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Cosh<T>(ADiff<T> a) where T : struct, IVectorND<T>
            => new ADiff<T>(Cosh(a.Value), a.Derivative * Sinh(a.Value));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Tanh<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double coshVal = Cosh(a.Value);
            double derivativeScale = 1.0 / (coshVal * coshVal);
            return new ADiff<T>(Tanh(a.Value), a.Derivative * derivativeScale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Asinh<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double scale = 1.0 / Sqrt(a.Value * a.Value + 1.0);
            return new ADiff<T>(Asinh(a.Value), a.Derivative * scale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Acosh<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double scale = 1.0 / Sqrt(a.Value * a.Value - 1.0);
            return new ADiff<T>(Acosh(a.Value), a.Derivative * scale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Atanh<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double scale = 1.0 / (1.0 - a.Value * a.Value);
            return new ADiff<T>(Atanh(a.Value), a.Derivative * scale);
        }

        // --- Linear Algebra / SpMV ---

        /// <summary>
        /// Compressed Sparse Row (CSR) matrix-vector product: y = A * x.
        /// Uses ReadOnlySpan for zero heap allocation across linear solver iterations.
        /// </summary>
        public static void Mult(
            ReadOnlySpan<double> a_value,
            ReadOnlySpan<int> a_index,
            ReadOnlySpan<int> a_start,
            ReadOnlySpan<double> x,
            Span<double> y)
        {
            int numRows = a_start.Length - 1;

            for (int i = 0; i < numRows; i++)
            {
                int rowStart = a_start[i];
                int rowEnd = a_start[i + 1];

                double sum = 0.0;
                for (int k = rowStart; k < rowEnd; k++)
                {
                    sum += a_value[k] * x[a_index[k]];
                }

                y[i] = sum;
            }
        }
    }
}