using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class Well
    {
        int lastindex = 0, nextindex = 1;
        public string Name { get; }
        public int[] Perforation_NatIndex { get { return perforation_NatIndex; } }
        public double[] Perforation_WI { get { return productivity_Index; } }
        public double MinPressure { get; }
        public double MaxPressure { get; }
        public double ComputeRate(double time)
        {
            var (Time, Rate) = ConstraintType switch
            {
                ConstraintType.LiqRate => LiqRateProfile,
                ConstraintType.OilRate => OilRateProfile,
                ConstraintType.WaterRate => WaterRateProfile,
                ConstraintType.GasRate => GasRateProfile,
                ConstraintType.Rate => RateProfile,
                _ => throw new NotImplementedException()
            };

            while (nextindex < Time.Count && Time[nextindex] < time)
                (lastindex, nextindex) = (nextindex, nextindex + 1);
            double qlast = lastindex > 0 ? Rate[lastindex-1] : 0;
            double rate = Rate[lastindex], ramp = 1 - Exp(Time[lastindex] - time);
            return WellType switch
            {
                WellType.Producer => -Abs(qlast + ramp*(rate - qlast)),
                WellType.Injector => Abs(qlast + ramp*(rate - qlast)),
                _ => 0,
            };
        }
        public double OilRate { get; set; }
        public double WaterRate { get; set; }
        public double LiqRate { get { return OilRate + WaterRate; } }
        public double GasRate { get; set; }
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
        public List<double> ProfileTime { get => profileTime; }
        public ConstraintType ConstraintType { get; set; } = ConstraintType.LiqRate;
        public WellType WellType { get; set; }
        public (List<double> Time, List<double> Rate) LiqRateProfile { get; set; }
        public (List<double> Time, List<double> Rate) RateProfile { get; set; }
        public (List<double> Time, List<double> Rate) OilRateProfile { get; set; }
        public (List<double> Time, List<double> Rate) WaterRateProfile { get; set; }
        public (List<double> Time, List<double> Rate) GasRateProfile { get; set; }


        int[] perforation_NatIndex;
        double[] productivity_Index;
        List<double> profileTime;


        public Well(WellType welltype, string name, double radius, double skin, double minPressure, double maxPressure,
            int i, int j, int[] perfInterval, List<double> time, List<double> orate = null, List<double> wrate = null, 
            List<double> grate = null, List<double> lrate = null, List<double> rate = null)
        {
            WellType = welltype;
            Name = name;
            Radius = radius;
            Skin = skin;
            MinPressure = minPressure;
            MaxPressure = maxPressure;
            I = i; J = j;
            PerfInterval = perfInterval;
            profileTime = time;
            if (lrate is not null)
            { LiqRateProfile = (Time: time, Rate: lrate); ConstraintType = ConstraintType.LiqRate; }
            if (orate is not null)
            { OilRateProfile = (Time: time, Rate: orate); ConstraintType = ConstraintType.OilRate; }
            if (wrate is not null)
            { WaterRateProfile = (Time: time, Rate: wrate); ConstraintType = ConstraintType.WaterRate; }
            if (grate is not null)
            { GasRateProfile = (Time: time, Rate: grate); ConstraintType = ConstraintType.GasRate; }
            if (rate is not null)
            { RateProfile = (Time: time, Rate: rate); ConstraintType = ConstraintType.Rate; }
        }

        internal void ComputeNaturalIndex(int Nx, int Ny)
        {
            var rng = Enumerable.Range(PerfInterval[0], PerfInterval[1] - PerfInterval[0] + 1);
            perforation_NatIndex = [.. rng.Select(k => I + J*Nx + k*Nx*Ny)];
        }

        internal void ComputeProductivityIndex(double[] Kx, double[] Ky, double[] Dx, double[] Dy, double[] Dz)
        {
            productivity_Index = new double[PerfInterval[1] - PerfInterval[0] + 1];
            int i = 0; double alpha_well = 1.127e-3*2*pi;
            foreach (int m in perforation_NatIndex)
            {
                double re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                productivity_Index[i++] = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/Radius) + Skin);
            }
        }
        internal ADiff Constraint(double time, ADiff Pressure, ADiff Rate)
        {
            return ConstraintType switch
            {
                ConstraintType.MaxPressure => Pressure - MaxPressure,
                ConstraintType.MinPressure => Pressure - MinPressure,
                _ => Rate - ComputeRate(time)
            };
        }

        internal ADiff Constraint(double time, ADiff Pressure, ADiff Rate, ADiff ConstraintValue)
        {
            return ConstraintType switch
            {
                ConstraintType.MaxPressure => ConstraintValue.CopyFrom(Pressure).SubtractInPlace(MaxPressure),
                ConstraintType.MinPressure => ConstraintValue.CopyFrom(Pressure).SubtractInPlace(MinPressure),
                _ => ConstraintValue.CopyFrom(Rate).SubtractInPlace(ComputeRate(time)),
            };
        }

        internal double[][] GetTrajectoryCoordinates(double[] xCoords, double[] yCoords, double[] zCoords)
        {
            double x1 = xCoords[I] + 0.5*(xCoords[I+1] - xCoords[I]), x2 = x1;
            double y1 = yCoords[J] + 0.5*(yCoords[J+1] - yCoords[J]), y2 = y1;
            double z1 = zCoords[0] - 100, z2 = zCoords[PerfInterval[1]+1];
            double[][] coordinates = [[x1, y1, z1], [x2, y2, z2]];
            return coordinates;
        }
    }
}
