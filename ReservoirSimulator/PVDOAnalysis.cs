using SepalSolver;

namespace ReservoirSimulator
{
    public class PVDOAnalysis
    {
        List<PvdoRow> pvdoTable;
        public class PvdoResult
        {
            public double PRef { get; set; }                    // The specified anchor pressure
            public double ConstantCompressibility { get; set; } // c_o derived from Bo slope
            public double BoAtPRef { get; set; }               // Fitted Bo exactly at PRef
            public double ViscosityAtPRef { get; set; }        // Fitted Viscosity exactly at PRef
            public double ViscosityExponent { get; set; }       // c_mu derived from mu slope
            public double BoRSquared { get; set; }              // Quality of Bo model fit
            public double ViscRSquared { get; set; }            // Quality of Viscosity model fit
        }

        public class PvdoRow
        {
            public double Pressure { get; set; }
            public double Bo { get; set; }
            public double Viscosity { get; set; }

            public PvdoRow(double pressure, double bo, double viscosity)
            {
                Pressure = pressure;
                Bo = bo;
                Viscosity = viscosity;
            }
        }

        public PVDOAnalysis(Matrix _pvdo)
        {
            pvdoTable = [..Enumerable.Range(0, _pvdo.Rows).
                Select(i => new PvdoRow(_pvdo[i, 0], _pvdo[i, 1], _pvdo[i, 2]))];
        }
        public PvdoResult FitPvdoProperties(double pRef)
        {
            if (pvdoTable == null || pvdoTable.Count < 2)
            {
                throw new ArgumentException("PVDO table requires at least 2 rows to establish exponential parameters.");
            }

            var sorted = pvdoTable.OrderBy(r => r.Pressure).ToList();

            // Coordinate transformations
            List<double> deltaP = []; // x = (P - PRef)
            List<double> lnBo = [];   // y1 = ln(Bo)
            List<double> lnVisc = []; // y2 = ln(mu_o)

            foreach (var row in sorted)
            {
                if (row.Viscosity <= 0 || row.Bo <= 0)
                    throw new InvalidOperationException("PVT properties must be strictly greater than zero for log conversion.");

                deltaP.Add(row.Pressure - pRef);
                lnBo.Add(Math.Log(row.Bo));
                lnVisc.Add(Math.Log(row.Viscosity));
            }

            // Run regressions in the transformed log-linear spaces
            var (mBo, bBo, rSqBo) = LinearRegression(deltaP, lnBo);
            var (mVisc, bVisc, rSqVisc) = LinearRegression(deltaP, lnVisc);

            // Back-transform intercepts to resolve physical units
            double co = -mBo;                      // Bo decreases with P -> co is -slope
            double boAtPref = Math.Exp(bBo);       // e^intercept_Bo
            double cMu = mVisc;                     // Viscosity increases with P -> c_mu is slope
            double viscAtPref = Math.Exp(bVisc);   // e^intercept_Visc

            return new PvdoResult
            {
                PRef = pRef,
                ConstantCompressibility = co,
                BoAtPRef = boAtPref,
                ViscosityAtPRef = viscAtPref,
                ViscosityExponent = cMu,
                BoRSquared = rSqBo,
                ViscRSquared = rSqVisc
            };
        }

        private static (double slope, double intercept, double rSquared) LinearRegression(List<double> x, List<double> y)
        {
            int n = x.Count;
            double sumX = x.Sum();
            double sumY = y.Sum();
            double sumXY = x.Zip(y, (xi, yi) => xi * yi).Sum();
            double sumXSq = x.Select(xi => xi * xi).Sum();

            double denominator = (n * sumXSq) - (sumX * sumX);
            if (Math.Abs(denominator) < 1e-11) return (0.0, 0.0, 1.0);

            double m = ((n * sumXY) - (sumX * sumY)) / denominator;
            double b = (sumY - (m * sumX)) / n;

            double yMean = sumY / n;
            double ssTot = y.Select(yi => Math.Pow(yi - yMean, 2)).Sum();
            double ssRes = y.Zip(x, (yi, xi) => Math.Pow(yi - (m * xi + b), 2)).Sum();
            double rSq = ssTot > 0 ? 1.0 - (ssRes / ssTot) : 1.0;

            return (m, b, rSq);
        }
    }
}