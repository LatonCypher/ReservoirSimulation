using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvdgResult(double pRef, double cg, double bgAtPref, double cMug,
                                      double viscGAtPref, double bgRSquared, double viscGRSquared)
    {
        public readonly double PRef = pRef;
        public readonly double Cg = cg;                // Gas compressibility factor
        public readonly double BgAtPref = bgAtPref;
        public readonly double CMug = cMug;            // Gas viscosity exponent
        public readonly double ViscGAtPref = viscGAtPref;

        public double BgRSquared { get; } = bgRSquared;
        public double ViscGRSquared { get; } = viscGRSquared;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBg<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return BgAtPref * Exp(-Cg * (p - PRef));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMug<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return ViscGAtPref * Exp(CMug * (p - PRef));
        }
    }

    public class PvdgRow(double pressure, double bg, double viscosity)
    {
        public double Pressure { get; set; } = pressure;
        public double Bg { get; set; } = bg;
        public double Viscosity { get; set; } = viscosity;
    }

    public class PVDGAnalysis(List<double> _pvdg)
    {
        private readonly List<PvdgRow> pvdgTable = [..Enumerable.Range(0, _pvdg.Count/3).
            Select(i => new PvdgRow(_pvdg[3*i+0], _pvdg[3*i+1], _pvdg[3*i+2]))];

        public PvdgResult Result { get; private set; }

        public void FitProperties()
        {
            if (pvdgTable == null || pvdgTable.Count < 2)
                throw new ArgumentException("PVDG table requires at least 2 rows.");

            var sorted = pvdgTable.OrderBy(r => r.Pressure).ToList();
            double pRef = sorted[0].Pressure;

            List<double> deltaP = [];
            List<double> lnBg = [];
            List<double> lnViscG = [];

            foreach (var row in sorted)
            {
                if (row.Viscosity <= 0 || row.Bg <= 0)
                    throw new InvalidOperationException("Properties must be strictly positive.");

                deltaP.Add(row.Pressure - pRef);
                lnBg.Add(Log(row.Bg));
                lnViscG.Add(Log(row.Viscosity));
            }

            var (mBg, bBg, rSqBg) = LinearRegression(deltaP, lnBg);
            var (mViscG, bViscG, rSqViscG) = LinearRegression(deltaP, lnViscG);

            Result =  new(pRef, -mBg, Exp(bBg), mViscG, Exp(bViscG), rSqBg, rSqViscG);
        }

        private static (double slope, double intercept, double rSquared) LinearRegression(List<double> x, List<double> y)
        {
            int n = x.Count;
            double sumX = x.Sum(); double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumXSq = x.Sum(xi => Sqr(xi));

            double denominator = (n * sumXSq) - (sumX * sumX);
            if (Abs(denominator) < 1e-11) return (0.0, 0.0, 1.0);

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            double yMean = sumY / n;
            double ssTot = y.Sum(yi => Sqr(yi - yMean));
            double ssRes = y.Zip(x, (yi, xi) => Sqr(yi - (m * xi + b))).Sum();
            return (m, b, ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0);
        }
    }
}