using SepalSolver;
using System.Runtime.CompilerServices;
namespace ReservoirSimulator
{
    internal static class Math
    {
        private const double Epsilon = 1e-14;
        internal const double pi = SepalSolver.Math.pi;
        internal static int DerivativeCapacity = 31;

        public static double Floor(double x) =>
            SepalSolver.Math.Floor(x);

        public static double[] Linspace(double x, double y, int N) =>
            SepalSolver.Math.Linspace(x, y, N);

        public static double Sqr(double x) =>
            SepalSolver.Math.Sqr(x);

        public static double Sin(double x) =>
            SepalSolver.Math.Sin(x);

        public static double Cos(double x) =>
            SepalSolver.Math.Cos(x);

        public static double Tan(double x) =>
           SepalSolver.Math.Tan(x);

        public static double Asin(double x) =>
           SepalSolver.Math.Asin(x);

        public static double Acos(double x) =>
            SepalSolver.Math.Acos(x);

        public static double Atan(double x) =>
            SepalSolver.Math.Atan(x);

        public static double Sinh(double x) =>
            SepalSolver.Math.Sinh(x);

        public static double Cosh(double x) =>
            SepalSolver.Math.Cosh(x);

        public static double Tanh(double x) =>
            SepalSolver.Math.Tanh(x);

        public static double Asinh(double x) =>
           SepalSolver.Math.Asinh(x);

        public static double Acosh(double x) =>
            SepalSolver.Math.Acosh(x);

        public static double Atanh(double x) =>
            SepalSolver.Math.Atanh(x);
        public static double Sqrt(double x) =>
            SepalSolver.Math.Sqrt(x);

        public static double Cbrt(double x) =>
            SepalSolver.Math.Pow(x, 1.0/3.0);

        public static double Log(double x) =>
            SepalSolver.Math.Log(x);

        public static double Exp(double x) =>
        SepalSolver.Math.Exp(x);

        public static double Hypot(double x, double y) =>
           SepalSolver.Math.Hypot(x, y);

        public static double Min(double x, double y) => x <= y ? x : y;

        public static double Max(double x, double y) => x >= y ? x : y;
        public static int Min(int x, int y) => x <= y ? x : y;

        public static int Max(int x, int y) => x >= y ? x : y;

        public static int Sign(double x) => (x > 0) ? 1 : ((x < 0) ? -1 : 0);

        public static double Pow(double x, double y) =>
           SepalSolver.Math.Pow(x, y);
        public static double Abs(double x) => x >= 0 ? x : -x;


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
            // Optimizing the derivative evaluation path to minimize redundant Math.Pow calls
            double derivativeScale = b * Pow(a.Value, b - 1.0);
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
            return new ADiff<T>(logVal, a.Derivative / a.Value);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Sqrt<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double sqrtVal = Sqrt(a.Value);
            return new ADiff<T>(sqrtVal, a.Derivative / (2 * sqrtVal));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Cbrt<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double cbrtVal = Cbrt(a.Value);
            // d(∛u) = du / (3 * ∛u^2)
            double derivativeScale = 1.0 / (3.0 * (cbrtVal * cbrtVal));
            return new ADiff<T>(cbrtVal, a.Derivative * derivativeScale);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Sin<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            // Computing both together if needed, but standard Math calls are heavily optimized by the JIT
            return new ADiff<T>(Sin(a.Value), a.Derivative * Cos(a.Value));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Cos<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            return new ADiff<T>(Cos(a.Value), -a.Derivative * Sin(a.Value));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Tan<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double tanVal = Tan(a.Value);
            double cosVal = Cos(a.Value);
            // d(tan u) = du / cos^2(u)
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
        {
            return new ADiff<T>(Sinh(a.Value), a.Derivative * Cosh(a.Value));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Cosh<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            return new ADiff<T>(Cosh(a.Value), a.Derivative * Sinh(a.Value));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff<T> Tanh<T>(ADiff<T> a) where T : struct, IVectorND<T>
        {
            double coshVal = Cosh(a.Value);
            // d(tanh u) = du / cosh^2(u)
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

        public static double[] Mult(List<double> a_value, List<int> a_index, List<int> a_start, double[] x)
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
