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
        internal ADiff Po { get; set; } = new();
        internal ADiff Sw { get; set; } = new();
        internal ADiff Pb { get; set; } = new();
        internal ADiff Sg { get; set; } = new();

        internal Vector3[] Corners { get; set; } = [];

    }
}
