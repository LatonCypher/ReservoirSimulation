using System;
using System.Collections.Generic;
using System.Text;

namespace ReservoirSimulator
{
    public class Face
    {
        public List<int> CellIndices { get; set; } = [];
        public List<double> Transmissibility { get; set; } = [];
        public FlowDirection Direction { get; set; }
        public double[] Coordinates { get; set; }
    }
}
