using System.Numerics;

namespace ReservoirSimulator
{
    public class Grid
    {
        public int[] Indices {get; set;} = [];
        public double NTG { get; set; }
        public double Porosity { get; set; } = 0;
        public Face[] Faces { get; set; } = [];
        public double Volume { get; set; } = 0;
        public double PoreVolume { get; set; } = 0;
        internal AutoDiff Po { get; set; } = new() { Value = 0, Derivative = new Vector4d(1, 0, 0, 0) };
        internal AutoDiff Sw { get; set; } = new() { Value = 0, Derivative = new Vector4d(0, 1, 0, 0) };
        internal AutoDiff Pb { get; set; } = new() { Value = 0, Derivative = new Vector4d(0, 0, 1, 0) };
        internal AutoDiff Sg { get; set; } = new() { Value = 0, Derivative = new Vector4d(0, 0, 0, 1) };
        internal AutoDiff Pw, So, Pg;

        internal Func<AutoDiff, AutoDiff> Swe, Pcow, Pcog, Bo, Bw, Bg, μo, μw, γo, γw, Kro, Krw, Krg, Er;
        internal Vector3[] Corners { get; set; } = [];

        internal void UpdateGridParameters()
        {

        }
    }
}
