using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class Well
    {
        int lastindex = 0, nextindex = 1;
        public string Name { get; }
        public int[] Perforation_NatIndex { get { return perforation_NatIndex; } }
        public double MinPressure { get; }
        public double MaxPressure { get; }
        public double ProdRate(double time)
        {
            while (nextindex < ProductionProfile.Time.Count && ProductionProfile.Time[nextindex] < time)
                (lastindex, nextindex) = (nextindex, nextindex + 1);
            double rate = ProductionProfile.Rate[lastindex];
            switch (WellType)
            {
                case WellType.Producer:
                    return -Abs(rate);
                case WellType.Injector:
                    return Abs(rate);
                default:
                    return 0;
            }
        }
        public double OilRate { get; set; }
        public double WaterRate { get; set; }
        public double WaterCut
        {
            get
            {
                double totalRate = OilRate + WaterRate;
                return WaterRate*100 / totalRate;
            }
        }
        public double Radius { get; }
        public double Skin { get; }
        public double Zref { get; set; }
        public int I { get; }
        public int J { get; }
        public int[] PerfInterval { get; }
        public ConstraintType ConstraintType { get; set; } = ConstraintType.FlowRate;
        public WellType WellType { get; set; }
        public (List<double> Time, List<double> Rate) ProductionProfile { get; set; }

        int[] perforation_NatIndex;

        public Well(WellType welltype, string name, double radius, double skin, double minPressure, double maxPressure,
            int i, int j, int[] perfInterval, List<double> time, List<double> rate)
        {
            WellType = welltype;
            Name = name;
            Radius = radius;
            Skin = skin;
            MinPressure = minPressure;
            MaxPressure = maxPressure;
            I = i; J = j;
            PerfInterval = perfInterval;
            ProductionProfile = (Time: time, Rate: rate);
        }

        internal void ComputeNaturalIndex(int Nx, int Ny)
        {
            var rng = Enumerable.Range(PerfInterval[0], PerfInterval[1] - PerfInterval[0] + 1);
            perforation_NatIndex = [.. rng.Select(k => I + J*Nx + k*Nx*Ny)];
        }
        internal ADiff Constraint(double time, ADiff Pressure, ADiff Rate )
        {
            return ConstraintType switch
            {
                ConstraintType.FlowRate => Rate - ProdRate(time),
                ConstraintType.MaxPressure => Pressure - MaxPressure,
                ConstraintType.MinPressure => Pressure - MinPressure,
                _ => Rate - ProdRate(time),
            };
        }

        internal ADiff Constraint(double time, ADiff Pressure, ADiff Rate, ADiff ConstraintValue)
        {
            ConstraintValue.Clear();
            return ConstraintType switch
            {
                ConstraintType.FlowRate => ConstraintValue.AddInPlace(Rate).SubtractInPlace(ProdRate(time)),
                ConstraintType.MaxPressure => ConstraintValue.AddInPlace(Pressure).SubtractInPlace(MaxPressure),
                ConstraintType.MinPressure => ConstraintValue.AddInPlace(Pressure).SubtractInPlace(MinPressure),
                _ => Rate - ProdRate(time),
            };
        }
    }
}
