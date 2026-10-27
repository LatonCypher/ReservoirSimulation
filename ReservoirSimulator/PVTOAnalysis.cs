using SepalSolver;
using System.Drawing.Drawing2D;
using System.Runtime.CompilerServices;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct PvtoResult(double rsBase, double rsExponent, double rsRSquared,
                      double boBase, double boExponent, double boRSquared,
                      double viscBase, double viscExponent, double viscRSquared)
    {
        // Raw Coefficients held as raw fields (Blittable, ultra-lean)
        public readonly double RsBase = rsBase;
        public readonly double RsExponent = rsExponent;
        public readonly double RsRSquared = rsRSquared;

        public readonly double BoBase = boBase;
        public readonly double BoExponent = boExponent;
        public readonly double BoRSquared = boRSquared;

        public readonly double ViscBase = viscBase;
        public readonly double ViscExponent = viscExponent;
        public readonly double ViscRSquared = viscRSquared;

        // High-Performance Inlined Mathematical Models
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateRs<T>(ADiff<T> p) where T : struct, IVectorND<T>
        {
            return RsBase * Pow(p, RsExponent);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateBo<T>(ADiff<T> rs) where T : struct, IVectorND<T>
        {
            return BoBase * Exp(BoExponent * rs);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateMuo<T>(ADiff<T> rs) where T : struct, IVectorND<T>
        {
            return ViscBase * Exp(-ViscExponent * rs);
        }
    }

    public class PvtoRow(double rs, double pressure, double bo, double viscosity)
    {
        public double Rs { get; set; } = rs;
        public double Pressure { get; set; } = pressure;
        public double Bo { get; set; } = bo;
        public double Viscosity { get; set; } = viscosity;
    }

    public class PVTOAnalysis
    {
        readonly List<PvtoRow> saturatedPoints = [];

        public PVTOAnalysis(List<double> _pvto)
        {
            double currentRs = -1.0;
            PvtoRow potentialBubblePoint = null;

            for (int i = 0; i < _pvto.Count / 4; i++)
            {
                double rs = _pvto[4*i+0];
                double p = _pvto[4*i+1];
                double bo = _pvto[4*i+2];
                double mu = _pvto[4*i+3];

                // When a new Rs block is encountered
                if (Abs(rs - currentRs) > 1e-5)
                {
                    if (potentialBubblePoint != null)
                        saturatedPoints.Add(potentialBubblePoint);
                    currentRs = rs;
                    potentialBubblePoint = new PvtoRow(rs, p, bo, mu);
                }
                else
                {
                    // Within an Rs block, Bo increases up to the bubble point.
                    // If Bo begins to drop, the previous row was the saturated inflection point.
                    if (bo > potentialBubblePoint.Bo)
                        potentialBubblePoint = new PvtoRow(rs, p, bo, mu);

                }
            }

            // Catch the final record block
            if (potentialBubblePoint != null)
                saturatedPoints.Add(potentialBubblePoint);
        }

        public PvtoResult Result { get; private set; }  

        public void FitProperties()
        {
            if (saturatedPoints.Count < 2)
            {
                throw new ArgumentException("PVTO analysis requires at least 2 distinct Rs saturation entries to regress slopes.");
            }

            // Arrays for Rs, Bo, Viscosity regressions (Log-linear spaces)
            List<double> rsValues = [];
            List<double> lnBo = [];
            List<double> lnVisc = [];

            // Arrays for Rs vs Pressure regression (Log-Log power law space)
            List<double> lnPressure = [];
            List<double> lnRs = [];

            foreach (var row in saturatedPoints)
            {
                if (row.Viscosity <= 0 || row.Bo <= 0 || 
                    row.Pressure <= 0 || row.Rs <= 0)
                    continue;

                rsValues.Add(row.Rs);
                lnBo.Add(Log(row.Bo));
                lnVisc.Add(Log(row.Viscosity));

                // Log transformations for Rs = alpha * P^beta
                lnPressure.Add(Log(row.Pressure));
                lnRs.Add(Log(row.Rs));
            }

            // 1. Run regressions for fluid properties driven by Rs
            var (mBo, bBo, rSqBo) = LinearRegression(rsValues, lnBo);
            var (mVisc, bVisc, rSqVisc) = LinearRegression(rsValues, lnVisc);

            // 2. Run the log-log regression for Rs driven by Pressure along the saturation boundary
            var (mRs, bRs, rSqRs) = LinearRegression(lnPressure, lnRs);

            Result = new (Exp(bRs), mRs, rSqRs, Exp(bBo), mBo, rSqBo, Exp(bVisc), -mVisc, rSqVisc );
        }

        private static (double slope, double intercept, double rSquared) LinearRegression(List<double> x, List<double> y)
        {
            int n = x.Count;
            double sumX = x.Sum();
            double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumXSq = x.Sum(Sqr);

            double denominator = (n * sumXSq) - (sumX * sumX);
            if (Abs(denominator) < 1e-11) return (0.0, 0.0, 1.0);

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            double yMean = sumY / n;
            double ssTot = y.Select(yi => Sqr(yi - yMean)).Sum();
            double ssRes = y.Zip(x, (yi, xi) => Sqr(yi - (m * xi + b))).Sum();
            double rSq = ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0;

            return (m, b, rSq);
        }
    }
}