using EclipseDeckReader;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class Well<T> where T : struct, IVectorND<T>
    {
        int lastindex = 0, nextindex = 1;
        public string Name { get; }
        public Phase PreferedPhase { get; set; }
        public List<int> PerfNatIndex { get; set; } = [];
        public List<double> Perforation_WI { get; set; } = [];
        public double MinPressure { get; }
        public double MaxPressure { get; }

        public double Target_Qo { get; set; }
        public double Qo { get; set; }
        public double Target_Qw { get; set; }
        public double Qw { get; set; }
        public double Target_Ql { get; set; }
        public double Ql { get { return Qo + Qw; } }
        public double Target_Qg { get; set; }
        public double Qg { get; set; }
        public double Fw { get{  return Qw*100 / (Qo + Qw); } }
        public double Zref { get; set; }
        public double ZControl { get; set; }
        public int I { get; }
        public int J { get; }
        public int ref_Index { get; set; } 
        public ConstraintType ConstraintType { get; set; } = ConstraintType.LRATE;
        public WellType WellType { get; set; }

        public int ActiveIndex { get; set; }
        public double HighestZref { get; set; }
        public double LowestZref { get; set; }

        public int[] Index;
        public double Pwf_n, Qflow_n;
        public ADiff<T> Pwf_np1, Qflow_np1;
        double ResWell, ResCon;
        IVectorND<T>[] JacResWell, JacResCon;
        IVectorND<T> JacResWellm, JacResConm;
        WellData WellData;

        public Well(WellData welldata)
        {
            WellData = welldata;
            Name = WellData.Specification.WellName;
            Zref = WellData.Specification.ReferenceDepth;
            I = WellData.Specification.HeadI;
            J = WellData.Specification.HeadJ;
            ref_Index = WellType == WellType.PRODUCER ? int.MaxValue : 0;
        }

        internal void ComputeProductivityIndex<T>(int Nx, int Ny, double[] Kx, double[] Ky, Grid<T>[] Grids, DateTime Time) where T : struct, IVectorND<T>
        {
            if (WellData.CompletionHistory[0].EffectiveDate > Time) return;

            PerfNatIndex.Clear(); Perforation_WI.Clear();
            int m = 0; double alpha_well = 1.127e-3 * 2 * pi;
            while (Time > WellData.ControlHistory[m].EffectiveDate) m++;
            var Control = WellData.ControlHistory[m];
            switch (Control.WellType)
            {
                case "PRODUCER":
                    WellType = WellType.PRODUCER;
                    ZControl = double.MaxValue;
                    break;
                case "INJECTOR":
                    WellType = WellType.INJECTOR;
                    ZControl = double.MinValue;
                    break;
            }
            m = 0;
            while (Time > WellData.CompletionHistory[m].EffectiveDate) m++;
            var comp = WellData.CompletionHistory[m];
            double radius = comp.WellBoreDiameter / 2, skin = comp.EffectiveSkin;
            foreach (int k in Enumerable.Range(comp.KStart, comp.KEnd - comp.KStart + 1))
            {
                if (comp.Status != "OPEN") continue;
                PerfNatIndex.Add(m = I + J * Nx + k * Nx * Ny);
                var grid = Grids[m]; grid.numWells++;

                switch (Control.WellType)
                {
                    case "PRODUCER":
                        WellType = WellType.PRODUCER;
                        if (grid.Center.Z > ZControl)
                        {
                            ZControl = grid.Center.Z;
                            Pwf_n = grid.Po_n;
                        }
                        break;
                    case "INJECTOR":
                        WellType = WellType.INJECTOR;
                        if (grid.Center.Z < ZControl)
                        {
                            ZControl = grid.Center.Z;
                            Pwf_n = grid.Po_n;
                        }
                        break;
                }
                double DX = 0, DY = 0, DZ = 0;
                for (int c = 0; c < 8; c++)
                {
                    // Determine local binary topological offsets
                    double fx = (c == 1 || c == 2 || c == 5 || c == 6) ? 1 : -1;
                    double fy = (c == 2 || c == 3 || c == 6 || c == 7) ? 1 : -1;
                    double fz = (c >= 4) ? 1 : -1;
                    DX += fx * grid.Corners[c].X;
                    DY += fy * grid.Corners[c].Y;
                    DZ += fz * grid.Corners[c].Z;
                }
                DX /= 4; DY /= 4; DZ /= 4;
                double re = 0.28 * Hypot(Pow(Ky[m] / Kx[m], 0.25) * DX, Pow(Kx[m] / Ky[m], 0.25) * DY) /
                                        (Pow(Ky[m] / Kx[m], 0.25)   +   Pow(Kx[m] / Ky[m], 0.25));
                Perforation_WI.Add(alpha_well * Sqrt(Kx[m] * Ky[m]) * DZ / (Log(re / radius) + skin));
            }
        }

        internal ADiff<T> Constraint<T>(DateTime Time, ADiff<T> Pressure, ADiff<T> Rate) where T : struct, IVectorND<T>
        {
            ADiff<T> Cons = 0;
            int m = 0;
            Target_Qo = 0; Target_Qw = 0; Target_Qg = 0; Target_Ql = 0;
            while (WellData.ControlHistory[m].EffectiveDate > Time) m++;
            var Control = WellData.ControlHistory[m];
            switch (Control.ControlMode)
            {
                case "ORAT":
                    Target_Qo = -Control.OilRateTarget;
                    Cons = Rate - Target_Qo;
                    break;
                case "WRAT":
                    Target_Qw = -Control.WaterRateTarget;
                    Cons = Rate - Target_Qw;
                    break;
                case "GRAT":
                    Target_Qg = -Control.GasRateTarget;
                    Cons = Rate - Target_Qg;
                    break;
                case "LRAT":
                    Target_Ql = -Control.LiquidRateTarget;
                    Cons = Rate - Target_Ql;
                    break;
                case "RATE":
                    switch (WellData.Specification.PreferredPhase)
                    {
                        case "WATER":
                            PreferedPhase = Phase.WATER;
                            Target_Qw = Control.TargetRate;
                            Cons = Rate - Target_Qw;
                            break;
                        case "GAS":
                            PreferedPhase = Phase.GAS;
                            Target_Qg = Control.TargetRate;
                            Cons = Rate - Target_Qg;
                            break;
                    }
                    break;
                case "BHP":
                    switch (WellData.Specification.PreferredPhase)
                    {
                        case "WATER":
                            PreferedPhase = Phase.WATER;
                            Target_Qw = Control.TargetRate;
                            Cons = Rate - Target_Qw;
                            break;
                        case "GAS":
                            PreferedPhase = Phase.GAS;
                            Target_Qg = Control.TargetRate;
                            Cons = Rate - Target_Qg;
                            break;
                    }
                    break;
            }
            return Cons;
        }

        public void ReconcileBHP()
        {

        }

        internal double[][] GetTrajectoryCoordinates<T>(int Nx, int Ny, int Nz, Grid<T>[] Grids) where T : struct, IVectorND<T>
        {
            var gridTop = Grids[I + J * Nx]; 
            var center1 = gridTop.Faces.First().Center;  
            double x1 = center1.X, y1 = center1.Y, z1 = center1.Z - 100;
            var gridBottom = Grids[I + J * Nx + (Nz - 1) * Nx * Ny];
            var center2 = gridBottom.Faces.Last().Center;
            double x2 = center2.X, y2 = center2.Y, z2 = center2.Z;
            double[][] coordinates = [[x1, y1, z1], [x2, y2, z2]];
            return coordinates;
        }
    }
}
