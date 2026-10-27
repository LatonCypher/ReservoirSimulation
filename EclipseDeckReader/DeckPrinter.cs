using System;
using System.Collections.Generic;
using System.Linq;

namespace EclipseDeckReader
{
    public static class DeckPrinter
    {
        /// <summary>
        /// Analyzes a populated EclipseDataDeck instance and outputs a structured summary to the console.
        /// </summary>
        public static void PrintSummary(EclipseDataDeck deck)
        {
            Console.WriteLine("==================================================================");
            Console.WriteLine("                  ECLIPSE DATA DECK SUMMARY                       ");
            Console.WriteLine("==================================================================");

            // 1. Grid Structural Block
            int totalCells = deck.NX * deck.NY * deck.NZ;
            Console.WriteLine("## GRID GEOMETRY");
            Console.WriteLine($"* Dimensions (NX x NY x NZ) : {deck.NX} x {deck.NY} x {deck.NZ}");
            Console.WriteLine($"* Total Grid Block Cells    : {totalCells:N0}");
            Console.WriteLine();

            // 2. Active Fluid Domain Phases
            Console.WriteLine("## FLUID SYSTEMS AND TRACKED PHASES");
            Console.WriteLine($"* Active Keyword Flags      : {string.Join(", ", deck.PhaseFlags)}");
            Console.WriteLine($"* Surface Oil Baseline Dens : {deck.SurfaceOilDensity:F2}");
            Console.WriteLine($"* Surface Water Baseline Dens: {deck.SurfaceWaterDensity:F2}");
            Console.WriteLine($"* Surface Gas Baseline Dens  : {deck.SurfaceGasDensity:F4}");
            Console.WriteLine();

            // 3. Spatial Parameter Arrays Profiler
            Console.WriteLine("## CONTINUOUS RECURRENT SPATIAL ARRAYS");
            PrintArrayStats("Porosity (PORO)", deck.PorosityArray);
            PrintArrayStats("Net-To-Gross (NTG)", deck.NtgArray);
            PrintArrayStats("Permeability X (PERMX)", deck.PermXArray);
            PrintArrayStats("Permeability Y (PERMY)", deck.PermYArray);
            PrintArrayStats("Permeability Z (PERMZ)", deck.PermZArray);
            PrintArrayStats("Transmissibility Mult X (MULTX)", deck.MultXArray);
            PrintArrayStats("Transmissibility Mult Y (MULTY)", deck.MultYArray);
            PrintArrayStats("Transmissibility Mult Z (MULTZ)", deck.MultZArray);
            PrintArrayStats("Coordinate Lines (COORD)", deck.CoordArray);
            PrintArrayStats("Corner Heights (ZCORN)", deck.ZcornArray);
            Console.WriteLine();

            // 4. Discrete Assignment Regions
            Console.WriteLine("## DISCRETE ASSIGNMENT REGIONS");
            PrintRegionStats("Active Cell Flags (ACTNUM)", deck.ActnumArray);
            PrintRegionStats("Saturation Regions (SATNUM)", deck.SatnumArray);
            PrintRegionStats("PVT Properties Regions (PVTNUM)", deck.PvtnumArray);
            PrintRegionStats("Equilibrium Regions (EQUINUM)", deck.EquinumArray); 
            Console.WriteLine();

            // 5. Multi-Region Hydrostatic Equilibration Summary block
            Console.WriteLine("## RESERVOIR INITIAL HYDROSTATIC EQUILIBRIUM");
            Console.WriteLine($"* Registered EQUIL Tables Count: {deck.RawEquilTables.Count}");
            for (int e = 0; e < deck.RawEquilTables.Count; e++)
            {
                var eq = deck.RawEquilTables[e];
                if (eq.Count >= 5)
                {
                    Console.WriteLine($"  └─ Equilibrium Region {e + 1} -> Datum Z: {eq[0]:F1} | Datum P: {eq[1]:F2} | WOC Z: {eq[2]:F1} | GOC Z: {eq[4]:F1}");
                }
            }
            Console.WriteLine();


            // 5. Rock Mechanical Properties
            Console.WriteLine("## ROCK MECHANICS AND PORE EXPANSION");
            Console.WriteLine($"* Reference Pressure ($P_{{ref}}$)  : {deck.RockReferencePressure:F2} psi");
            Console.WriteLine($"* Formation Compressibility ($c_r$): {deck.RockCompressibility:E4} 1/psi");
            Console.WriteLine();

            // 6. Fault Multipliers
            Console.WriteLine("## FAULT MULTIPLIERS");
            Console.WriteLine($"* Registered Faults Count  : {deck.Faults.Count}");
            foreach (var kvp in deck.Faults)
            {
                Console.WriteLine($"  └─ Fault: {kvp.Key,-12} | Face: {kvp.Value.Face} | Mult: {kvp.Value.TransmissibilityMultiplier:F4}");
            }
            Console.WriteLine();

            // 7. Dynamic Well Management Block (Added Well Network Summary Engine)
            Console.WriteLine("## WELL GEOMETRY AND OPERATIONS NETWORK");
            Console.WriteLine($"* Total Wells Registered    : {deck.Wells.Count}");
            if (deck.Wells.Count == 0)
            {
                Console.WriteLine("* [NO WELL NETWORKS REGISTED OR PARSED]");
            }
            else
            {
                foreach (var wellKvp in deck.Wells)
                {
                    string name = wellKvp.Key;
                    WellData data = wellKvp.Value;

                    Console.WriteLine($"\n* Well Entity: {name}");

                    // A. Print Master Well Specifications
                    if (data.Specification != null)
                    {
                        var spec = data.Specification;
                        Console.WriteLine($"  ├─ Group Reference : {spec.GroupName,-10} | Head (I,J): ({spec.HeadI}, {spec.HeadJ}) | Ref Depth: {spec.ReferenceDepth:F1} | Primary Phase: {spec.PreferredPhase}");
                    }
                    else
                    {
                        Console.WriteLine("  ├─ Group Reference : [NO INITIAL SPECIFICATION RECORD LOADED]");
                    }

                    // B. Print Interval Completions (COMPDAT)
                    Console.WriteLine($"  ├─ Interval Completions Count: {data.CompletionHistory.Count}");
                    for (int c = 0; c < data.CompletionHistory.Count; c++)
                    {
                        var comp = data.CompletionHistory[c];
                        Console.WriteLine($"  │  └─ Comp [{c + 1}]: Local Cell Coordinates (I, J): ({comp.I}, {comp.J}) | K-Layers: {comp.KStart}-{comp.KEnd} | Status: {comp.Status,-4} | Diameter: {comp.WellBoreDiameter:F3} | Skin: {comp.EffectiveSkin:F2}");
                    }

                    // C. Print Operational Constraints Over Time Timeline (WCONPROD / WCONINJE)
                    Console.WriteLine($"  └─ Control Constraint Schedules Count: {data.ControlHistory.Count}");

                    // Sort historical profiles chronologically to show clear step actions
                    var sortedControls = data.ControlHistory.OrderBy(ch => ch.EffectiveDate).ToList();
                    for (int h = 0; h < sortedControls.Count; h++)
                    {
                        var ctrl = sortedControls[h];
                        Console.Write($"     └─ Date [{ctrl.EffectiveDate:yyyy-MM-dd}] Mode: {ctrl.ControlMode,-5} ");

                        if (ctrl.ControlMode.Contains("RAT")) // Rate-driven constraints
                        {
                            Console.WriteLine($"| Targets -> Oil: {ctrl.OilRateTarget:F1} | Water: {ctrl.WaterRateTarget:F1} | Gas: {ctrl.GasRateTarget:F1} | Liq: {ctrl.LiquidRateTarget:F1} | Min BHP Limit: {ctrl.MinBhpLimit:F1}");
                        }
                        else // Pressure or Injection control constraints
                        {
                            Console.WriteLine($"| Limits  -> Max BHP Limit: {ctrl.MaxBhpLimit:F1} | Fixed Water Inj Rate: {ctrl.WaterRateTarget:F1} | Fixed Gas Inj Rate: {ctrl.GasRateTarget:F1}");
                        }
                    }
                }
            }
            Console.WriteLine();

            // 8. Fluid Properties Summary Block
            Console.WriteLine("## FLUID PVT REGIONS COUNT");
            PrintTableGroupOverview("Oil Live PVT Properties (PVTO)", deck.RawPvtoTables);
            PrintTableGroupOverview("Oil Dead PVT Properties (PVDO)", deck.RawPvdoTables);
            PrintTableGroupOverview("Water PVT Properties (PVTW)", deck.RawPvtwTables);
            PrintTableGroupOverview("Gas Live PVT Properties (PVTG)", deck.RawPvtgTables);
            PrintTableGroupOverview("Gas Dead PVT Properties (PVDG)", deck.RawPvdgTables);
            Console.WriteLine();

            // 9. Rock-Fluid Properties Block
            Console.WriteLine("## ROCK-FLUID INTERFACIAL RELATIONS");
            PrintTableGroupOverview("Water-Oil Saturation Functions (SWOF)", deck.RawSwofTables);
            PrintTableGroupOverview("Gas-Oil Saturation Functions (SGOF)", deck.RawSgofTables);
            Console.WriteLine("==================================================================");
        }

        private static void PrintArrayStats(string arrayName, List<double> array)
        {
            if (array == null || array.Count == 0)
            {
                Console.WriteLine($"* {arrayName,-35}: [NOT PRESENT OR EMPTY]");
                return;
            }
            double min = array.Min();
            double max = array.Max();
            double avg = array.Average();
            Console.WriteLine($"* {arrayName,-35}: Loaded Count = {array.Count:N0} | Min = {min,7:F4} | Max = {max,7:F4} | Avg = {avg,7:F4}");
        }

        private static void PrintRegionStats(string arrayName, List<int> array)
        {
            if (array == null || array.Count == 0)
            {
                Console.WriteLine($"* {arrayName,-35}: [NOT PRESENT OR EMPTY]");
                return;
            }
            int min = array.Min();
            int max = array.Max();
            int uniqueRegions = array.Distinct().Count();
            Console.WriteLine($"* {arrayName,-35}: Loaded Count = {array.Count:N0} | Distinct Regions = {uniqueRegions} (IDs: {min} to {max})");
        }

        private static void PrintTableGroupOverview(string tableName, List<List<double>> tables)
        {
            if (tables == null || tables.Count == 0)
            {
                Console.WriteLine($"* {tableName,-35}: [0 Tables Found]");
                return;
            }

            Console.WriteLine($"* {tableName,-35}: [{tables.Count} Table Region(s) Registered]");
            for (int r = 0; r < tables.Count; r++)
            {
                Console.WriteLine($"  └─ Region {r + 1}: Parsed Numeric Values Flat Elements Count = {tables[r].Count}");
            }
        }
    }
}