using SepalSolver;
using System.Drawing.Drawing2D;
using System.Runtime.CompilerServices;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvtgResult(double rvBase, double rvExponent, double rvRSquared,
                      double bgBase, double bgExponent, double bgRSquared,
                      double viscGBase, double viscGExponent, double viscGRSquared)
    {
        // Rv vs P power law coefficients
        public readonly double RvBase = rvBase;
        public readonly double RvExponent = rvExponent;
        public readonly double RvRSquared = rvRSquared;

        // Bg vs Rv exponential coefficients
        public readonly double BgBase = bgBase;
        public readonly double BgExponent = bgExponent;
        public readonly double BgRSquared = bgRSquared;

        // Viscosity vs Rv exponential coefficients
        public readonly double ViscGBase = viscGBase;
        public readonly double ViscGExponent = viscGExponent;
        public readonly double ViscGRSquared = viscGRSquared;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateRv<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return RvBase * Pow(p, RvExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBg<T>(ADiff<T> rv) where T : struct, IVectorND<T>
        {
            return BgBase * Exp(BgExponent * rv);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMug<T>(ADiff<T> rv) where T : struct, IVectorND<T>
        {
            return ViscGBase * Exp(ViscGExponent * rv);
        }
    }

    public class PvtgRow(double rv, double pressure, double bg, double viscosity)
    {
        public double Rv { get; set; } = rv;
        public double Pressure { get; set; } = pressure;
        public double Bg { get; set; } = bg;
        public double Viscosity { get; set; } = viscosity;
    }
    public class PVTGAnalysis
    {
        private readonly List<PvtgRow> dewPoints = [];

        public PvtgResult Result { get; private set; }

        public PVTGAnalysis(List<double> _pvtg)
        {
            double currentRv = -1.0;
            PvtgRow potentialDewPoint = null;

            for (int i = 0; i < _pvtg.Count/4; i++)
            {
                double rv = _pvtg[4*i+0];
                double p = _pvtg[4*i+1];
                double bg = _pvtg[4*i+2];
                double mu = _pvtg[4*i+3];

                if (Math.Abs(rv - currentRv) > 1e-7)
                {
                    if (potentialDewPoint != null) dewPoints.Add(potentialDewPoint);
                    currentRv = rv;
                    potentialDewPoint = new PvtgRow(rv, p, bg, mu);
                }
                else
                {
                    // For gas, Bg decreases with pressure up to the dew point inflection
                    if (bg < potentialDewPoint.Bg)
                    {
                        potentialDewPoint = new PvtgRow(rv, p, bg, mu);
                    }
                }
            }
            if (potentialDewPoint != null) dewPoints.Add(potentialDewPoint);
        }

        public PvtgResult result;

        public void FitProperties()
        {
            if (dewPoints.Count < 2)
                throw new ArgumentException("PVTG analysis requires at least 2 distinct entries.");

            List<double> rvValues = [];
            List<double> lnBg = [];
            List<double> lnViscG = [];
            List<double> lnPressure = [];
            List<double> lnRv = [];

            foreach (var row in dewPoints)
            {
                if (row.Viscosity <= 0 || row.Bg <= 0 || row.Pressure <= 0 || row.Rv <= 0)
                    continue;

                rvValues.Add(row.Rv);
                lnBg.Add(Log(row.Bg));
                lnViscG.Add(Log(row.Viscosity));

                lnPressure.Add(Log(row.Pressure));
                lnRv.Add(Log(row.Rv));
            }

            var (mBg, bBg, rSqBg) = LinearRegression(rvValues, lnBg);
            var (mViscG, bViscG, rSqViscG) = LinearRegression(rvValues, lnViscG);
            var (mRv, bRv, rSqRv) = LinearRegression(lnPressure, lnRv);

            result = new (Exp(bRv), mRv, rSqRv, Exp(bBg), mBg, rSqBg, Exp(bViscG), mViscG, rSqViscG);
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