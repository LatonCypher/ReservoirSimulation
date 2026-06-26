using MKLNET;
using SepalSolver;
using System.Drawing;
using System.Runtime.CompilerServices;
using static SepalSolver.Ode;
using static SepalSolver.Statistics;
using static System.Runtime.InteropServices.JavaScript.JSType;
namespace ReservoirSimulator
{
    internal static class Math
    {
        private const double Epsilon = 1e-14;
        internal const double pi = SepalSolver.Math.pi;
        internal static int DerivativeCapacity = 31;

        public static double[] Linspace(double x, double y, int N) =>
            SepalSolver.Math.Linspace(x, y, N);

        public static double Sqrt(double x) =>
            SepalSolver.Math.Sqrt(x);

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


        public static double Sqr(double x) => x*x;

        public static double Abs(double x) => x >= 0 ? x : -x;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Abs(ADiff a) => (a.Value >= 0.0) ? a : -a;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Max(ADiff a, ADiff b) => (a.Value >= b.Value) ? a : b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Min(ADiff a, ADiff b) => (a.Value <= b.Value) ? a : b;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Pow(ADiff a, double b) => new ADiff().AddInPlace(a).PowInPlace(b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Exp(ADiff a) => new ADiff().AddInPlace(a).ExpInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Log(ADiff a) => new ADiff().AddInPlace(a).LogInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Sqrt(ADiff a) => new ADiff().AddInPlace(a).SqrtInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Cbrt(ADiff a) => new ADiff().AddInPlace(a).CbrtInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Sin(ADiff a) => new ADiff().AddInPlace(a).SinInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Cos(ADiff a) => new ADiff().AddInPlace(a).CosInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Tan(ADiff a) => new ADiff().AddInPlace(a).TanInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Asin(ADiff a) => new ADiff().AddInPlace(a).AsinInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Acos(ADiff a) => new ADiff().AddInPlace(a).AcosInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Atan(ADiff a) => new ADiff().AddInPlace(a).AtanInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Sinh(ADiff a) => new ADiff().AddInPlace(a).SinhInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Cosh(ADiff a) => new ADiff().AddInPlace(a).CoshInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Tanh(ADiff a) => new ADiff().AddInPlace(a).TanhInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Asinh(ADiff a) => new ADiff().AddInPlace(a).AsinhInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Acosh(ADiff a) => new ADiff().AddInPlace(a).AcoshInPlace();

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static ADiff Atanh(ADiff a) => new ADiff().AddInPlace(a).AtanhInPlace();

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
