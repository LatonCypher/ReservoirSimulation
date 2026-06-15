using FFmpeg.AutoGen;
using System;
using SepalSolver;

namespace ReservoirSimulator
{
    public class SWOFAnalysis
    {
        List<SwofRow> swofTable;
        public SWOFAnalysis(Matrix _swof)
        {
            swofTable = [..Enumerable.Range(0, _swof.Rows).
                Select(i => new SwofRow(_swof[i, 0], _swof[i, 1], _swof[i, 2], _swof[i, 3]))];
        }
        public class PcResult
        {
            public double Swi { get; set; }      // Irreducible Water Saturation
            public double Pd { get; set; }       // Entry/Displacement Pressure
            public double Lambda { get; set; }   // Pore-size distribution index
            public double RSquared { get; set; } // Quality of fit metric
        }

        public class PermResult
        {
            public double Swi { get; set; }
            public double Sor { get; set; }
            public double KrwEndpoint { get; set; } // k_rw^o
            public double KroEndpoint { get; set; } // k_ro^o
            public double Nw { get; set; }          // Water Corey exponent
            public double No { get; set; }          // Oil Corey exponent
            public double KrwRSquared { get; set; }
            public double KroRSquared { get; set; }
        }

        public class SwofRow
        {
            public double Sw { get; set; }
            public double Krw { get; set; }
            public double Krow { get; set; }
            public double Pc { get; set; }

            public SwofRow(double sw, double krw, double krow, double pc)
            {
                Sw = sw;
                Krw = krw;
                Krow = krow;
                Pc = pc;
            }
        }

        public PcResult FitBrooksCoreyPc()
        {
            if (swofTable == null || swofTable.Count < 3)
            {
                throw new ArgumentException("The SWOF table must contain at least 3 rows to perform regression.");
            }

            // 1. Identify Swi from the very first row (assumed sorted by Sw ascending)
            var sortedTable = swofTable.OrderBy(r => r.Sw).ToList();
            double swi = sortedTable.First().Sw;

            // 2. Filter data: Logarithms cannot handle Pc <= 0 or Sw* <= 0
            // We need rows where Sw > Swi and Pc > 0
            List<double> xLogValues = new List<double>();
            List<double> yLogValues = new List<double>();

            foreach (var row in sortedTable)
            {
                if (row.Sw > swi && row.Pc > 0)
                {
                    double swStar = (row.Sw - swi) / (1.0 - swi);

                    xLogValues.Add(SepalSolver.Math.Log(swStar)); // ln(Sw*)
                    yLogValues.Add(SepalSolver.Math.Log(row.Pc));   // ln(Pc)
                }
            }

            if (xLogValues.Count < 2)
            {
                return new PcResult
                {
                    Swi = swi,
                    Pd = 0,
                    Lambda = 1,
                    RSquared = 1
                };
            }

            // 3. Perform Least-Squares Linear Regression (y = mx + b)
            int n = xLogValues.Count;
            double sumX = xLogValues.Sum();
            double sumY = yLogValues.Sum();
            double sumXY = xLogValues.Zip(yLogValues, (x, y) => x * y).Sum();
            double sumXSq = xLogValues.Sum(x => x * x);
            double sumYSq = yLogValues.Sum(y => y * y);

            // Calculate Slope (m) and Intercept (b)
            double denominator = (n * sumXSq) - (sumX * sumX);
            if (SepalSolver.Math.Abs(denominator) < 1e-9)
            {
                throw new InvalidOperationException("Linear regression failed due to a zero denominator (possible vertical line data).");
            }

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            // 4. Extract Brooks-Corey Parameters
            double lambda = -1.0 / m;
            double pd = SepalSolver.Math.Exp(b);

            // 5. Compute R-Squared (Coefficient of Determination) for validation
            double yMean = sumY / n;
            double ssTot = yLogValues.Select(y => SepalSolver.Math.Pow(y - yMean, 2)).Sum();
            double ssRes = yLogValues.Zip(xLogValues, (y, x) => SepalSolver.Math.Pow(y - (m * x + b), 2)).Sum();
            double rSquared = ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0;

            return new PcResult
            {
                Swi = swi,
                Pd = pd,
                Lambda = lambda,
                RSquared = rSquared
            };
        }

        public PermResult FitCoreyPermeability()
        {
            if (swofTable == null || swofTable.Count < 3)
            {
                throw new ArgumentException("Table requires at least 3 rows for regression.");
            }

            var sorted = swofTable.OrderBy(r => r.Sw).ToList();

            // 1. Resolve critical saturation endpoints directly from the table boundaries
            double swi = sorted.First().Sw;
            double kroEndpoint = sorted.First().Krow; // Kro at Swi

            // Find where oil stops flowing (Krow == 0) to get Sor
            var sorRow = sorted.FirstOrDefault(r => r.Krow <= 0.0);
            if (sorRow == null)
                throw new InvalidOperationException("Could not find a row where Krow cleanly reaches 0 to extract Sor.");

            double swMax = sorRow.Sw; // Sw at residual oil
            double sor = 1.0 - swMax;
            double krwEndpoint = sorRow.Krw; // Krw at Sor

            // 2. Prepare Data Arrays for Regressions
            var krwDataX = new List<double>();
            var krwDataY = new List<double>();

            var kroDataX = new List<double>();
            var kroDataY = new List<double>();

            foreach (var row in sorted)
            {
                // Domain checking to avoid Log(0) or Log(Negative)
                if (row.Sw > swi && row.Sw < swMax)
                {
                    double swStar = (row.Sw - swi) / (1.0 - swi - sor);

                    if (row.Krw > 0)
                    {
                        krwDataX.Add(SepalSolver.Math.Log(swStar));
                        krwDataY.Add(SepalSolver.Math.Log(row.Krw));
                    }

                    if (row.Krow > 0)
                    {
                        kroDataX.Add(SepalSolver.Math.Log(1.0 - swStar));
                        kroDataY.Add(SepalSolver.Math.Log(row.Krow));
                    }
                }
            }

            // 3. Compute Linear Regressions
            var (nw, ln_krwO, rSqW) = LinearRegression(krwDataX, krwDataY);
            var (no, ln_kroO, rSqO) = LinearRegression(kroDataX, kroDataY);

            return new PermResult
            {
                Swi = swi,
                Sor = sor,
                // We use the computed endpoints if the fit is perfect, 
                // but taking them directly from data matches physical table constraints better.
                KrwEndpoint = krwEndpoint,
                KroEndpoint = kroEndpoint,
                Nw = nw,
                No = no,
                KrwRSquared = rSqW,
                KroRSquared = rSqO
            };
        }

        private static (double slope, double intercept, double rSquared) LinearRegression(List<double> x, List<double> y)
        {
            int n = x.Count;
            if (n < 2) return (1.0, 0.0, 0.0);

            double sumX = x.Sum();
            double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumXSq = x.Select(xi => xi * xi).Sum();

            double denominator = (n * sumXSq) - (sumX * sumX);
            if (SepalSolver.Math.Abs(denominator) < 1e-9) return (1.0, 0.0, 0.0);

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            // R-Squared validation
            double yMean = sumY / n;
            double ssTot = y.Select(yi => SepalSolver.Math.Pow(yi - yMean, 2)).Sum();
            double ssRes = y.Zip(x, (yi, xi) => SepalSolver.Math.Pow(yi - (m * xi + b), 2)).Sum();
            double rSq = ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0;

            return (m, b, rSq);
        }
    }
}