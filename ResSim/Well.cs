namespace ResSim
{
    public class Well
    {
        public string Name { get; }
        public int[] Perforation_NatIndex { get { return perforation_NatIndex; } }
        public double MinPressure { get; }
        public double MaxPressure { get; }
        public double ProdRate(double time)
        {
            int index = ProductionProfile.Time.FindLastIndex(t => t <= time);
            return index >= 0 ? ProductionProfile.Rate[index] : 0;
        }
        public double OilRate { get; set; }
        public double WaterRate { get; set; }
        public double Radius { get; }
        public double Skin { get; }
        public double Zref { get; set; }
        public int I { get; }
        public int J { get; }
        public int[] PerfInterval { get; }
        public ConstraintType ConstraintType { get; set; } = ConstraintType.FlowRate;
        public WellType WellType { get; set; } = WellType.Producer;

        public (List<double> Time, List<double> Rate) ProductionProfile { get; set; }

        int[] perforation_NatIndex;



        public Well(string name, double radius, double skin, double minPressure, double maxPressure, 
            double oilRate, double waterRate, int i, int j, int[] perfInterval)
        {
            Name = name;
            Radius = radius;
            MinPressure = minPressure;
            MaxPressure = maxPressure;
            OilRate = oilRate;
            WaterRate = waterRate;
            I = i; J = j;
            PerfInterval = perfInterval;
        }
        

        public void ComputeNaturalIndex(int Nx, int Ny)
        {
            perforation_NatIndex = [..PerfInterval.Select(k => I + J*Nx + k*Nx*Ny)];
        }

        public double Constraint(double time, double Pressure, double Rate)
        {
            return ConstraintType switch
            {
                ConstraintType.FlowRate => Rate - ProdRate(time),
                ConstraintType.MaxPressure => Pressure - MaxPressure,
                ConstraintType.MinPressure => Pressure - MinPressure,
                _ => Rate - ProdRate(time),
            };
        }

        public AutoDiff Constraint(double time, AutoDiff Pressure, AutoDiff Rate)
        {
            return ConstraintType switch
            {
                ConstraintType.FlowRate => Rate - ProdRate(time),
                ConstraintType.MaxPressure => Pressure - MaxPressure,
                ConstraintType.MinPressure => Pressure - MinPressure,
                _ => Rate - ProdRate(time),
            };
        }
    }
}
