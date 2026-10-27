using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct SwofResult(double swi, double sor, double krwEndpoint, double kroEndpoint,
                          double nw, double no, double pd, double lambda,
                          double pcRSquared, double krwRSquared, double kroRSquared)
    {
        // Endpoints
        public readonly double Swi = swi;
        public readonly double Sor = sor;
        public readonly double KrwEndpoint = krwEndpoint;
        public readonly double KroEndpoint = kroEndpoint;

        // Corey Exponents
        public readonly double Nw = nw;
        public readonly double No = no;

        // Brooks-Corey Capillary parameters
        public readonly double Pd = pd;
        public readonly double InvLambda = 1.0 / lambda; // Pre-inverted (1.0 / lambda) for register speed

        // Diagnostics
        public double PcRSquared { get; } = pcRSquared;
        public double KrwRSquared { get; } = krwRSquared;
        public double KroRSquared { get; } = kroRSquared;

        // High-Performance Inlined AutoDiff Math Models
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrw<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            // Normalization: swNorm = (sw - swi) / (1.0 - swi - sor)
            ADiff<T> swNorm = (sw - Swi) / (1.0 - Swi - Sor);
            return KrwEndpoint * Pow(swNorm, Nw);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrow<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            ADiff<T> swNorm = (sw - Swi) / (1.0 - Swi - Sor);
            return KroEndpoint * Pow(1.0 - swNorm, No);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluatePcow<T>(ADiff<T> sw) where T : struct, IVectorND<T>
        {
            // Brooks-Corey Sw* = (sw - swi) / (1.0 - swi)
            ADiff<T> swStar = (sw - Swi) / (1.0 - Swi);
            return Pd * Pow(swStar, -InvLambda);
        }
    }

    public class SwofRow(double sw, double krw, double krow, double pc)
    {
        public double Sw { get; set; } = sw;
        public double Krw { get; set; } = krw;
        public double Krow { get; set; } = krow;
        public double Pc { get; set; } = pc;
    }

    public class SWOFAnalysis(List<double> _swof)
    {
        public readonly List<SwofRow> swofTable = [..Enumerable.Range(0, _swof.Count/4).
            Select(i => new SwofRow(_swof[4*i+0], _swof[4*i+1], _swof[4*i+2], _swof[4*i+3]))];

        public SwofResult Result { get; private set; }

        public void FitProperties()
        {
            if (swofTable == null || swofTable.Count < 3)
                throw new ArgumentException("The SWOF table must contain at least 3 rows.");

            var sorted = swofTable.OrderBy(r => r.Sw).ToList();

            // 1. Core Endpoints
            double swi = sorted.First().Sw;
            double kroEndpoint = sorted.First().Krow;

            var sorRow = sorted.FirstOrDefault(r => r.Krow <= 0.0) ?? sorted.Last();
            double swMax = sorRow.Sw;
            double sor = 1.0 - swMax;
            double krwEndpoint = sorRow.Krw;

            // 2. Prepare Log Arrays
            List<double> pcX = []; List<double> pcY = [];
            List<double> krwX = []; List<double> krwY = [];
            List<double> kroX = []; List<double> kroY = [];

            foreach (var row in sorted)
            {
                if (row.Sw > swi && row.Pc > 0)
                {
                    double swStar = (row.Sw - swi) / (1.0 - swi);
                    pcX.Add(Log(swStar));
                    pcY.Add(Log(row.Pc));
                }

                if (row.Sw > swi && row.Sw < swMax)
                {
                    double swNorm = (row.Sw - swi) / (1.0 - swi - sor);
                    if (row.Krw > 0)
                    {
                        krwX.Add(Log(swNorm));
                        krwY.Add(Log(row.Krw));
                    }
                    if (row.Krow > 0)
                    {
                        kroX.Add(Log(1.0 - swNorm));
                        kroY.Add(Log(row.Krow));
                    }
                }
            }

            // 3. Execute Least-Squares Regressions
            var (mPc, bPc, rSqPc) = LinearRegression(pcX, pcY);
            var (nw, _, rSqW) = LinearRegression(krwX, krwY);
            var (no, _, rSqO) = LinearRegression(kroX, kroY);

            double lambda = -1.0 / mPc;
            double pd = Exp(bPc);

            Result = new(swi, sor, krwEndpoint, kroEndpoint, nw, no, pd, lambda, rSqPc, rSqW, rSqO);
        }

        private static (double slope, double intercept, double rSquared) LinearRegression(List<double> x, List<double> y)
        {
            int n = x.Count;
            if (n < 2) return (1.0, 0.0, 0.0);

            double sumX = x.Sum(); double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumXSq = x.Sum(xi => Sqr(xi));

            double denominator = (n * sumXSq) - (sumX * sumX);
            if (Abs(denominator) < 1e-11) return (1.0, 0.0, 0.0);

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            double yMean = sumY / n;
            double ssTot = y.Sum(yi => Sqr(yi - yMean));
            double ssRes = y.Zip(x, (yi, xi) => Sqr(yi - (m * xi + b))).Sum();
            double rSq = ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0;

            return (m, b, rSq);
        }
    }
}