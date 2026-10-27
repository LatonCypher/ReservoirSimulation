using SepalSolver;

namespace ReservoirSimulator
{
    public class Face
    {
        public int ActiveGridIndex { get; set; }
        public List<int> GridIndices { get; set; } = [];
        public List<double> Transmissibility { get; set; } = [];
        public double Permeability { get; set; }
        public FlowDirection Direction { get; set; }
        public Vector3d[] Corners { get; set; }

        public Vector3d Center { get; set; }

        public double PerpendicularDistance2Center { get; set; }
        public AABB BoundingBox  { get; set; }

        public void MakeBoundBox()
        {
            BoundingBox = new(Corners.Min(c => c.X), Corners.Min(c => c.Y), Corners.Min(c => c.Z),
                              Corners.Max(c => c.X), Corners.Max(c => c.Y), Corners.Max(c => c.Z));
        }
    }
}
