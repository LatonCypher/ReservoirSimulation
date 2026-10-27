using System;
using System.Collections.Generic;

namespace EclipseDeckReader
{
    public sealed class EclipseDataDeck
    {
        // 1. Grid Dimensions
        public int NX { get; set; }
        public int NY { get; set; }
        public int NZ { get; set; }

        // 2. Run Control Flags
        public List<string> PhaseFlags { get; } = [];

        // 3. Flat Spatial Parameter Arrays
        public List<double> PorosityArray { get; set; } = [];
        public List<double> PermXArray { get; set; } = [];
        public List<double> PermYArray { get; set; } = [];
        public List<double> PermZArray { get; set; } = [];
        public List<double> NtgArray { get; set; } = [];
        public List<double> MultXArray { get; set; } = [];
        public List<double> MultYArray { get; set; } = [];
        public List<double> MultZArray { get; set; } = [];
        public List<int> ActnumArray { get; set; } = [];
        public List<int> SatnumArray { get; set; } = [];
        public List<int> PvtnumArray { get; set; } = [];
        public List<int> EquinumArray { get; set; } = []; 

        // 4. Block-Centered Grid Dimension Arrays
        public List<double> DxArray { get; set; } = [];
        public List<double> DyArray { get; set; } = [];
        public List<double> DzArray { get; set; } = [];
        public List<double> ZtopArray { get; set; } = [];

        // 5. Corner-Point Geometry Coordinates
        public List<double> CoordArray { get; set; } = [];
        public List<double> ZcornArray { get; set; } = [];

        // 6. Surface Fluid Densities
        public double SurfaceOilDensity { get; set; }
        public double SurfaceWaterDensity { get; set; }
        public double SurfaceGasDensity { get; set; }

        // 7. Rock Compressibility Properties
        public double RockReferencePressure { get; set; }
        public double RockCompressibility { get; set; }

        // 8. Tabular Fluid, Rock, and Equilibration Data Blocks
        public List<List<double>> RawPvdoTables { get; } = [];
        public List<List<double>> RawPvtoTables { get; } = [];
        public List<List<double>> RawPvdgTables { get; } = [];
        public List<List<double>> RawPvtgTables { get; } = [];
        public List<List<double>> RawPvtwTables { get; } = [];
        public List<List<double>> RawSwofTables { get; } = [];
        public List<List<double>> RawSgofTables { get; } = [];
        public List<List<double>> RawEquilTables { get; } = []; 

        // 9. Fault Structures
        public Dictionary<string, FaultRecord> Faults { get; } = new(StringComparer.OrdinalIgnoreCase);

        // 10. Dynamic Well Networks Block
        public Dictionary<string, WellData> Wells { get; } = new(StringComparer.OrdinalIgnoreCase);

        // 11. Dynamic Timeline Synchronization State
        public DateTime StartDate { get; set; }                       // Master baseline anchor point (Day 0)
        public DateTime CurrentSimDate { get; set; } = DateTime.Now;  // Tracks active day as scheduler advances
    }
}