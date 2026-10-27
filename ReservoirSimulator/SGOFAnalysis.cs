using System.Runtime.CompilerServices;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public readonly struct SgofResult(double sgc, double sorg, double swi, double krgEndpoint, double krogEndpoint,
                                          double ng, double nog, double pdgo, double lambdaGo,
                                          double pcRSquared, double krgRSquared, double krogRSquared)
    {
        public readonly double Sgc = sgc;                      // Critical Gas Saturation
        public readonly double Sorg = sorg;                    // Residual Oil to Gas
        public readonly double Swi = swi;                      // Connate Water background parameter
        public readonly double KrgEndpoint = krgEndpoint;
        public readonly double KrogEndpoint = krogEndpoint;

        public readonly double Ng = ng;                        // Gas Corey Exponent
        public readonly double Nog = nog;                      // Oil-to-Gas Corey Exponent

        public readonly double Pdgo = pdgo;                    // Gas-oil displacement capillary entry
        public readonly double InvLambdaGo = 1.0 / lambdaGo;

        public double PcRSquared { get; } = pcRSquared;
        public double KrgRSquared { get; } = krgRSquared;
        public double KrogRSquared { get; } = krogRSquared;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrg<T>(ADiff<T> sg) where T : struct, IVectorND<T>
        {
            // sgNorm = (sg - sgc) / (1.0 - swi - sorg - sgc)
            ADiff<T> sgNorm = (sg - Sgc) / (1.0 - Swi - Sorg - Sgc);
            return KrgEndpoint * Pow(sgNorm, Ng);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluateKrog<T>(ADiff<T> sg) where T : struct, IVectorND<T>
        {
            ADiff<T> sgNorm = (sg - Sgc) / (1.0 - Swi - Sorg - Sgc);
            return KrogEndpoint * Pow(1.0 - sgNorm, Nog);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public ADiff<T> EvaluatePcgo<T>(ADiff<T> sg) where T : struct, IVectorND<T>
        {
            // Gas-oil capillary expansion operates via scaled gas track projection
            ADiff<T> sgStar = sg / (1.0 - Swi);
            return Pdgo * Pow(1.0 - sgStar, -InvLambdaGo);
        }
    }

    public class SgofRow(double sg, double krg, double krog, double pcg)
    {
        public double Sg { get; set; } = sg;
        public double Krg { get; set; } = krg;
        public double Krog { get; set; } = krog;
        public double Pcg { get; set; } = pcg;
    }

    public class SGOFAnalysis(List<double> _sgof)
    {
        public readonly List<SgofRow> sgofTable = [..Enumerable.Range(0, _sgof.Count/4).
            Select(i => new SgofRow(_sgof[4*i+0], _sgof[4*i+1], _sgof[4*i+2], _sgof[4*i+3]))];

        public SgofResult Result { get; private set; }

        public void FitProperties(double fixedSwi)
        {
            if (sgofTable == null || sgofTable.Count < 3)
                throw new ArgumentException("The SGOF table must contain at least 3 rows.");

            var sorted = sgofTable.OrderBy(r => r.Sg).ToList();

            // 1. Resolve Gas Endpoints
            double sgc = sorted.First().Sg; // Critical gas initialization entry
            double krogEndpoint = sorted.First().Krog; // Krog at Sg = 0

            var activeGasLimit = sorted.FirstOrDefault(r => r.Krog <= 0.0) ?? sorted.Last();
            double sgMax = activeGasLimit.Sg;
            double sorg = 1.0 - fixedSwi - sgMax;
            double krgEndpoint = activeGasLimit.Krg;

            List<double> pcX = []; List<double> pcY = [];
            List<double> krgX = []; List<double> krgY = [];
            List<double> krogX = []; List<double> krogY = [];

            foreach (var row in sorted)
            {
                if (row.Sg > sgc && row.Pcg > 0)
                {
                    double sgStar = row.Sg / (1.0 - fixedSwi);
                    pcX.Add(Log(1.0 - sgStar));
                    pcY.Add(Log(row.Pcg));
                }

                if (row.Sg > sgc && row.Sg < sgMax)
                {
                    double sgNorm = (row.Sg - sgc) / (1.0 - fixedSwi - sorg - sgc);
                    if (row.Krg > 0)
                    {
                        krgX.Add(Log(sgNorm));
                        krgY.Add(Log(row.Krg));
                    }
                    if (row.Krog > 0)
                    {
                        krogX.Add(Log(1.0 - sgNorm));
                        krogY.Add(Log(row.Krog));
                    }
                }
            }

            var (mPc, bPc, rSqPc) = LinearRegression(pcX, pcY);
            var (ng, _, rSqG) = LinearRegression(krgX, krgY);
            var (nog, _, rSqOg) = LinearRegression(krogX, krogY);

            double lambdaGo = -1.0 / mPc;
            double pdgo = Exp(bPc);

            Result = new(sgc, sorg, fixedSwi, krgEndpoint, krogEndpoint, ng, nog, pdgo, lambdaGo, rSqPc, rSqG, rSqOg);
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