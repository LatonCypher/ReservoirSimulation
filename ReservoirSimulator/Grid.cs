using System.Numerics;
using SepalSolver;

namespace ReservoirSimulator
{
    public class Grid<T> where T : struct, IVectorND<T>
    {
        // Define conversion constants
        const double alpha = 1.127e-3, beta = 5.615;
        public List<double> PoList, SwList, SgList, PbList, PdList;
        public bool LogResult { get; set; }
        public bool isActive { get; set; } = true;
        public bool isSaturated { get; set; } = true;
        public double NTG { get; set; }
        public double Porosity { get; set; } = 0;
        public Face[] Faces { get; set; } = [.. Enumerable.Range(0, 6).Select(i => new Face() { Direction = (FlowDirection)i })];
        public double PoreVolume { get; set; } = 0;
        public double Po_n, Pw_n, Pg_n, Sw_n, So_n, Sg_n, Pb_n, Pd_n, ErSo_Bo_n, ErSw_Bw_n, ErSg_Bg_n;
        public ADiff<T> Po_np1, Pw_np1, Pg_np1, Sw_np1, So_np1, Sg_np1, Pb_np1, Pd_np1, ErSo_Bo_np1, ErSw_Bw_np1, ErSg_Bg_np1;
        public Func<ADiff<T>, ADiff<T>> Krw, Krow, Krg, Krog, Pcow, Pcgo, Pcgw, Bo, Bw, Bg, μo, μw, μg, γo, γw, γg, Er, Rs, Rv; 
        public int ActiveIndex { get; set; }
        public (int i, int j, int k) IJKIndex { get; set; }

        internal Vector3d[] Corners { get; set; } = [];
        internal Vector3d Center { get; set; }
        double ResO, ResW, ResG, ResB, ResD;
        // If declaring as class fields or local variables
        IVectorND<T>[] JacResO, JacResW, JacResG, JacResB,  JacResD;
        IVectorND<T> JacResOm, JacResWm, JacResGm, JacResBm, JacResDm;
        internal List<int> indexOfConnectingGrids = [];
        internal int numWells = 0;

        internal void MakeFaces()
        {
            Center = new();
            for (int i = 0; i < 8; i++) Center += Corners[i];
            Center /= 8.0;
            Faces[0].Corners = [Corners[0], Corners[3], Corners[2], Corners[1]]; Faces[0].MakeBoundBox();
            Faces[1].Corners = [Corners[0], Corners[1], Corners[5], Corners[4]]; Faces[1].MakeBoundBox();
            Faces[2].Corners = [Corners[0], Corners[4], Corners[7], Corners[3]]; Faces[2].MakeBoundBox();
            Faces[3].Corners = [Corners[1], Corners[2], Corners[6], Corners[5]]; Faces[3].MakeBoundBox();
            Faces[4].Corners = [Corners[3], Corners[7], Corners[6], Corners[2]]; Faces[4].MakeBoundBox();
            Faces[5].Corners = [Corners[4], Corners[5], Corners[7], Corners[6]]; Faces[5].MakeBoundBox();
            PoreVolume = Utility.ComputeHexahedronVolume(this) / beta * Porosity;

            foreach (var face in Faces)
            {
                face.Center = new();
                foreach (var corner in face.Corners) face.Center += corner;
                face.Center /= 4.0;
                face.PerpendicularDistance2Center =
                    Utility.ComputePerpendicularDistance(Center, face.Corners[0], face.Corners[1], face.Corners[2]);
            }
        }

        internal void FinalizeGrid()
        {
            if (LogResult)
            { PoList = []; SwList = []; SgList = []; PbList = []; }
        }

        public void PrepareGrid<T>(double[] Xnp1)
        {
           
        }
    }
}