using EclipseDeckReader;
using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class Well<T> where T : struct, IVectorND<T>
    {
        public string Name { get; set; }
        public Phase PreferedPhase { get; set; }
        public List<int> PerfNatIndex { get; set; } = [];
        public List<double> Perforation_WI { get; set; } = [];
        public double MinPressure { get; }
        public double MaxPressure { get; }

        public double Target_Qo { get; set; }
        public ADiff<T> Qo { get; set; }
        public double Target_Qw { get; set; }
        public ADiff<T> Qw { get; set; }
        public double Target_Ql { get; set; }
        public ADiff<T> Ql { get { return Qo + Qw; } }
        public double Target_Qg { get; set; }
        public ADiff<T> Qg { get; set; }
        public double Fw { get{  return Qw.Value*100 / (Qo + Qw).Value; } }
        public double ZControl { get; set; }
        public int I { get; }
        public int J { get; }

        public int Nx { get; set; }
        public int Ny { get; set; }
        public int Nz { get; set; }
        public Grid<T>[] Grids { get; set; }
        public int ref_Index { get; set; } 
        public ConstraintType ConstraintType { get; set; } = ConstraintType.LRATE;
        public WellType WellType { get; set; }

        public WellControlRecord CurrentControl
        {
            get { return currentControl; }
            set
            {
                if(currentControl != value)
                {
                    currentControl = value;
                    var Control = currentControl;
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


                    ConstraintType = currentControl.ControlMode switch
                    {
                        "ORAT" => ConstraintType.ORATE,
                        "WRAT" => ConstraintType.WRATE,
                        "GRAT" => ConstraintType.GRATE,
                        "LRAT" => ConstraintType.LRATE,
                        "BHP" => WellType == WellType.INJECTOR ? ConstraintType.MAXP : ConstraintType.MINP,
                        _ => ConstraintType.ORATE // '_' acts as the default case
                    };

                }
            }
        }
        public WellCompletionRecord CurrentCompletionRecord
        {
            get { return currentcompletionRecord; }
            set
            {
                if(currentcompletionRecord != value)
                {
                    currentcompletionRecord = value;

                    var comp = currentcompletionRecord;
                    PerfNatIndex = [];
                    double radius = comp.WellBoreDiameter / 2, skin = comp.EffectiveSkin;
                    foreach (int k in Enumerable.Range(comp.KStart, comp.KEnd - comp.KStart + 1))
                    {
                        if (comp.Status != "OPEN") continue;
                        int m = I + J * Nx + k * Nx * Ny;
                        PerfNatIndex.Add(m);
                        var grid = Grids[m];

                        switch (CurrentControl.WellType)
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
                        double ky = grid.Faces[1].Permeability, kx = grid.Faces[2].Permeability;
                        double re = 0.28 * Hypot(Pow(ky / kx, 0.25) * DX, Pow(kx / ky, 0.25) * DY) /
                                                (Pow(ky / kx, 0.25) + Pow(kx / ky, 0.25));
                        grid.WI = alpha_well * Sqrt(kx * ky) * DZ / (Log(re / radius) + skin);
                    }
                }
            }
        }
        WellCompletionRecord currentcompletionRecord = null;
        WellControlRecord currentControl = null;
        public int ActiveIndex { get; set; }
        public double HighestZref { get; set; }
        public double LowestZref { get; set; }
        public double Pwf_n, Qflow_n;
        WellData WellData;
        double alpha_well = 1.127e-3 * 2 * pi;


        public Well(WellData welldata)
        {
            WellData = welldata;
            Name = WellData.Specification.WellName;
            Zref = WellData.Specification.ReferenceDepth;
            I = WellData.Specification.HeadI;
            J = WellData.Specification.HeadJ;
            ref_Index = WellType == WellType.PRODUCER ? int.MaxValue : 0;
        }

        internal List<DateTime> GetEventTimes()
        {
            return [.. WellData.ControlHistory.Select(history => history.EffectiveDate),
                    .. WellData.CompletionHistory.Select(history => history.EffectiveDate)];
        }

        internal void UpdateControl(DateTime time)
        {
            if (WellData.ControlHistory.Count == 0) return;
            if (WellData.ControlHistory[0].EffectiveDate > time) return;
            int m = 0;
            for (m = 0; m < WellData.ControlHistory.Count; m++)
            {
                if (time < WellData.ControlHistory[m].EffectiveDate)
                    break;
            }
            CurrentControl = WellData.ControlHistory[m - 1];
        }

        internal void UpdateCompletion(DateTime time)
        {
            if (WellData.CompletionHistory.Count == 0) return;
            if (WellData.CompletionHistory[0].EffectiveDate > time) return;
            int m = 0;
            for (m = 0; m < WellData.CompletionHistory.Count; m++)
            {
                if (time < WellData.CompletionHistory[m].EffectiveDate)
                    break;
            }
            CurrentCompletionRecord = WellData.CompletionHistory[m - 1];
        }

        internal ADiff<T> Constraint(ADiff<T> pressure,  ADiff<T> rate)
        {
            
            Target_Qo = 0; Target_Qw = 0; Target_Qg = 0; Target_Ql = 0;
            ADiff<T> cons = 0;
            switch (ConstraintType)
            {
                case ConstraintType.MINP:
                    {
                        // Producers draw down to a minimum BHP; 
                        double targetBhp = currentControl.MinBhpLimit;
                        cons = pressure - targetBhp;
                        break;
                    }

                case ConstraintType.MAXP:
                    {
                        // Injectors pump up to a maximum BHP
                        double targetBhp = currentControl.MaxBhpLimit;
                        cons = pressure - targetBhp;
                        break;
                    }

                case ConstraintType.ORATE:
                    {
                        Target_Qo = currentControl.OilRateTarget;
                        // Target is phase rate, so enforce against accumulated oil rate:
                        cons = Qo - Target_Qo;
                        break;
                    }

                case ConstraintType.WRATE:
                    {
                        Target_Qw = currentControl.WaterRateTarget;
                        cons = Qw - Target_Qw;
                        break;
                    }

                case ConstraintType.GRATE:
                    {
                        Target_Qg = currentControl.GasRateTarget;
                        cons = Qg - Target_Qg;
                        break;
                    }

                case ConstraintType.LRATE:
                    {
                        Target_Ql = currentControl.LiquidRateTarget;
                        cons = (Qo + Qw) - Target_Ql;
                        break;
                    }

                case ConstraintType.RATE:
                    {
                        // General rate target based on specification
                        double target = currentControl.TargetRate;
                        switch (WellData.Specification.PreferredPhase.ToUpperInvariant())
                        {
                            case "OIL":
                                PreferedPhase = Phase.OIL;
                                Target_Qo = target;
                                cons = Qo - Target_Qo;
                                break;

                            case "WATER":
                                PreferedPhase = Phase.WATER;
                                Target_Qw = target;
                                cons = Qw - Target_Qw;
                                break;

                            case "GAS":
                                PreferedPhase = Phase.GAS;
                                Target_Qg = target;
                                cons = Qg - Target_Qg;
                                break;

                            default: // TOTAL FLUID
                                cons = rate - target;
                                break;
                        }
                        break;
                    }

                default:
                    throw new NotSupportedException($"Unsupported well control mode: {currentControl.ControlMode}");
            }

            return cons;
        }

        internal double[][] GetTrajectoryCoordinates<T>(Grid<T>[] Grids) where T : struct, IVectorND<T>
        {
            var gridTop = Grids[PerfNatIndex.First()]; 
            var center1 = gridTop.Faces.First().Center;  
            double x1 = center1.X, y1 = center1.Y, z1 = center1.Z - 100;
            var gridBottom = Grids[PerfNatIndex.Last()];
            var center2 = gridBottom.Faces.Last().Center;
            double x2 = center2.X, y2 = center2.Y, z2 = center2.Z;
            double[][] coordinates = [[x1, y1, z1], [x2, y2, z2]];
            return coordinates;
        }
    }
}
