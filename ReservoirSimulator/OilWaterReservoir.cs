using EclipseDeckReader;
using SepalSolver;
using System.Text;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class OilWaterReservoir
    {
        public int funcall;
        public List<double[]> P, S, Rate, OilRate, WaterRate, GasRate, Pwf, WaterCut;
        public List<double> Time, SweepEff;
        readonly Func<ADiff<Vector2d>, ADiff<Vector2d>> Sws, Swe, Pc_D, Pc_I, Bo, Bw, μo, μw, γo, γw, Kro, Krw, Er;
        // Define conversion constants
        const double alpha = 1.127e-3,        // Darcy to Field units factor
                beta = 5.615;                 // ft3 to bbl conversion factor


        double [] xs;

        List<double> a_value = [], b = [];
        List<int> a_index = [], a_start = [0];

        static double betweenab(double a, double b, double f) => a + f*(b-a);
        static double interps(List<double> X, List<double> Y, double x)
        {
            int i = 0, j = 1;
            while (j < X.Count && X[j] < x) (i, j) = (j, j + 1);
            if (i == 0 || j == X.Count) return Y[i];
            double f = (x - X[i-1])/(X[i] - X[i-1]);
            return betweenab(Y[i-1], Y[i], f);
        }
        static double[] interpa(List<double> X, List<double[]> Y, double x)
        {
            int i = 0, j = 1;
            while (j < X.Count && X[j] < x) (i, j) = (j, j + 1);
            if (i == 0 || j == X.Count) return Y[i];
            double f = (x - X[i-1])/(X[i] - X[i-1]);
            return [.. Y[i-1].Zip(Y[i], (a, b) => betweenab(a, b, f))];
        }
        (double[] Po, double[] Sw, double[] Pwell, double[] Qwell, double maxDP, double maxDS, int i_p, int i_s) ExtractSolution()
        {
            int indx = 0, i_p = 0, i_s = 0; double maxDp = 0, maxDs = 0;
            double[] Po = [.. Enumerable.Repeat(double.NaN, Ngrids)],
                     Sw = [.. Enumerable.Repeat(double.NaN, Ngrids)],
                     Pwell = new double[Nwells], Qwell = new double[Nwells];
            for (int i = 0; i < Ngrids; i++)
            {
                if (Actnum[i] == 0) { indx += 2; continue; }
                Po[i] = xs[indx++]; // Matches Po index
                Sw[i] = xs[indx++]; // Matches Sw index
                if (Abs(Po[i] - Po_n[i]) > maxDp)
                { i_p = i; maxDp = Abs(Po[i] - Po_n[i]); }
                if (Abs(Sw[i] - Sw_n[i]) > maxDs)
                { i_s = i; maxDs = Abs(Sw[i] - Sw_n[i]); }
            }
            for (int i = 0; i < Nwells; i++)
            {
                Pwell[i] = xs[indx++]; // Matches Pwf index
                Qwell[i] = xs[indx++]; // Matches Q index
            }
            return (Po, Sw, Pwell, Qwell, maxDp, maxDs, i_p, i_s);
        }
        double dt, t = 0;
        void InitialGuess()
        {
            int indx = 0, n = Time.Count-1;
            double omega = n < 2 ? 0 : dt/(Time[n]-Time[n-1]);
            for (int i = 0; i < Ngrids; i++)
            {
                if (Actnum[i] == 0) { indx += 2; continue; }
                xs[indx++] = Po_n[i];                                   // Matches Po index
                if (n >= 2) xs[indx] += omega*(Po_n[i] - P[n-1][i]);    // Extrapolation
                xs[indx++] = Sw_n[i];                                   // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                xs[indx++] = Pwells_n[i];                               // Matches Pwf index
                xs[indx++] = Qwells_n[i];                               // Matches Q index
            }
        }
        double[] LinSolve(List<double> a_value, List<int> a_index, List<int> a_start, List<double> b)
        {
            return MklSparseSolver.Solve([.. a_value], [.. a_index], [.. a_start], [.. b]);
        }

        readonly int Nx, Ny, Nz, NxNy, Ngrids, Nwells, varNum;
        readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z, zTop;
        double[] Po_n, Sw_n, Pw_n, So_n, Qwells_n, Pwells_n, ErSw_Bw_n, ErSo_Bo_n, Coord, Zcorn, PoreVolume; 
        readonly double krw0, kro0, Pb, Prefw, Prefr, Pe, Pw_woc, Sw_r, So_r,
            Bo0, Bw0, Po_woc, Z_woc, co, cw, cr, bo, bw, nw, no, np,
            μo0, μw0, γo0, γw0, P_datum, Z_datum, Pc_max;
        Grid<Vector2d>[] Grids;
        List<Well<Vector2d>> Wells;
        Aquifer Aquifer;
        double[,] Trans;
        bool[,] Water_Upstream, Oil_Upstream;
        int[] Actnum, Satnum, Pvtnum;
        List<int> Columns2Delete = [];
        Dictionary<(int, int), EQUILPressureDepth> EquilPVTCombinations = [];
        EclipseDataDeck DATA;
        PVTWAnalysis[] PVTWs;
        PVDOAnalysis[] PVDOs;
        SWOFAnalysis[] SWOFs;

        public OilWaterReservoir(EclipseDataDeck Data)
        {
            DATA = Data;
            // RUNSPEC
            Nx = Data.NX; Ny = Data.NY; Nz = Data.NZ; NxNy = Nx * Ny; Ngrids = Nx * Ny * Nz;

            // GRIDS
            if (Data.DxArray.Count > 0 && Data.DyArray.Count > 0 && Data.DzArray.Count > 0)
            {
                Dx = [.. Data.DxArray]; Dy = [.. Data.DyArray]; Dz = [.. Data.DzArray];
                zTop = Data.ZtopArray.Count > 0 ? [.. Data.ZtopArray] : new double[NxNy];
            }
            else { Coord = [.. Data.CoordArray]; Zcorn = [.. Data.ZcornArray]; }
            Φ = [..Data.PorosityArray]; Kx = [..Data.PermXArray]; Ky = [..Data.PermYArray]; Kz = [..Data.PermZArray];
            Grids = Zcorn is not null ? Utility.ProcessGrid<Vector2d>(Nx, Ny, Nz, Coord, Zcorn, Φ, Kx, Ky, Kz) :
                Utility.ProcessGrid<Vector2d>(Nx, Ny, Nz, Dx, Dy, Dz, zTop, Φ, Kx, Ky, Kz);

            // PVTW
            PVTWs = [.. Data.RawPvtwTables.Select(pvtw => new PVTWAnalysis(pvtw))];
            foreach (var pvtw in PVTWs) pvtw.FitProperties();

            // PVDO
            PVDOs = [..Data.RawPvdoTables.Select(pvdo=>new PVDOAnalysis(pvdo))];
            foreach (var pvdo in PVDOs) pvdo.FitProperties();

            // SWOF
            SWOFs = [.. Data.RawSwofTables.Select(swof => new SWOFAnalysis(swof))];
            foreach (var swof in SWOFs) swof.FitProperties();

            // ROCK
            cr = Data.RockCompressibility; Prefr = Data.RockReferencePressure;
            Er = P => Exp(cr * (P - Prefr));

            // DENSITIES
            γo0 = Data.SurfaceOilDensity/144; γw0 = Data.SurfaceWaterDensity/144;

            // REGION

            Utility.SetAcitveGrids(Grids, Data.ActnumArray, Er);
            Utility.SetSWOF(Grids, Data.SatnumArray, SWOFs);
            Utility.SetPVDO(Grids, Data.PvtnumArray, PVDOs, γo0); 
            Utility.SetPVTW(Grids, Data.PvtnumArray, PVTWs, γw0);

            // TRANSMISSIBILITY
            Utility.ComputeFaceTransmissibility(Grids);

            // EQUIL
            EQUILAnalysis[] EQUILs = [.. Data.RawEquilTables.Select(equil => new EQUILAnalysis(equil))];
            double MaxZ = Grids.Max(grid => grid.Center.Z), MinZ = Grids.Min(grid => grid.Center.Z);
            for(int e = 0; e < EQUILs.Length; e++)
            {
                for (int p = 0; p < PVDOs.Length; p++)
                {
                    var equil = EQUILs[e]; var pvdo = PVDOs[p]; var pvtw = PVTWs[p];
                    double BoAtPref = pvdo.Result.BoAtPref, BwAtPref = pvtw.Result.BwAtPref;
                    double Co = pvdo.Result.Co, Cw = pvtw.Result.Cw;
                    double PRefo = pvdo.Result.PRef, PRefw = pvtw.Result.PRef;
                    Func<double, double> γo = P => γo0 / (BoAtPref * Exp(-Co * (P - PRefo))), 
                                         γw = P => γw0 / (BwAtPref * Exp(-Cw * (P - PRefw)));
                    var equildata = new EQUILPressureDepth(equil.Result);
                    int N = (int)(Abs(equildata.ZDatum - equildata.ZWOC) / (MaxZ - MinZ) * Nz * 10);
                    if (equildata.ZDatum < equildata.ZWOC)
                    {
                        var oil2WOC = Utility.IntegrateLegRK4(equildata.ZDatum, equildata.PDatum, equildata.ZWOC, N, γo);
                        equildata.PoWOC = oil2WOC.Last().P; equildata.PwWOC = equildata.PoWOC - equildata.PcWOC;
                    }
                    else
                    {
                        var water2WOC = Utility.IntegrateLegRK4(equil.ZDatum, equil.PDatum, equildata.ZWOC, N, γw);
                        equildata.PwWOC = water2WOC.Last().P; equildata.PoWOC = equildata.PwWOC + equildata.PcWOC;
                    }
                    N = (int)(Abs(MinZ - equildata.ZWOC) / (MaxZ - MinZ) * Nz * 10);
                    var OilPressureHeight1 = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PoWOC, MinZ, N, γo);
                    var WaterPressureHeight1 = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, MinZ, N, γw);
                    OilPressureHeight1.Reverse(); equildata.OilZonePressureDepthPairs.AddRange(OilPressureHeight1);
                    WaterPressureHeight1.Reverse(); equildata.WaterZonePressureDepthPairs.AddRange(WaterPressureHeight1);

                    N = (int)(Abs(MaxZ - equildata.ZWOC) / (MaxZ - MinZ) * Nz * 10);
                    var WaterPressureHeight2 = Utility.IntegrateLegRK4(equildata.ZWOC, equildata.PwWOC, MinZ, N, γw);
                    equildata.WaterZonePressureDepthPairs.AddRange(WaterPressureHeight2[1..]);
                    equildata.OilZonePressureDepthPairs.AddRange(WaterPressureHeight2[1..].
                        Select(pair => new PressureDepthPair(pair.Z, pair.P + equildata.PcWOC)));

                    EquilPVTCombinations[(e, p)] = equildata;
                }
            }
        }

        public void Initialize()
        {
            funcall = 0;
            for (int m = 0; m < Ngrids; m++)
            {
                var grid = Grids[m];
                if(!grid.isActive) continue;
                int s = DATA.SatnumArray[m], e = DATA.EquinumArray[m], p = DATA.PvtnumArray[m];

                var swoftable = SWOFs[s].swofTable;
                var equildata = EquilPVTCombinations[(e, p)];
                var oiltable = equildata.OilZonePressureDepthPairs;
                var watertable = equildata.WaterZonePressureDepthPairs;
                grid.Po_n = Utility.Interpolate(oiltable, grid.Center.Z, pair => pair.Z, pair => pair.P);
                grid.Pw_n = Utility.Interpolate(watertable, grid.Center.Z, pair => pair.Z, pair => pair.P);
                double pc = grid.Po_n - grid.Pw_n;
                if (pc == 0)
                    grid.Sw_n = grid.Center.Z < equildata.ZWOC ? swoftable.Min(swof => swof.Sw) : swoftable.Max(swof => swof.Sw);
                else grid.Sw_n = Utility.Interpolate(swoftable, pc, swof => swof.Pc, swof => swof.Sw);
                grid.So_n = 1 - grid.Sw_n;
            }

            foreach (var well in Wells)
            {
                well.ComputeProductivityIndex<Vector2d>(Nx, Ny, Kx, Ky, Grids, DATA.StartDate);
                well.Zref = 0;
            }
        }

        public static (List<int> start, List<int> index, List<double> value) DeleteCol(
    List<int> a_start, List<int> a_index, List<double> a_value, List<int> index2delete)
        {
            if (index2delete == null || index2delete.Count == 0)
                return (a_start, a_index, a_value);

            // 1. Convert deletion targets to a HashSet for O(1) lookups
            HashSet<int> deleteSet = [..index2delete];

            // 2. Build a global column remapping array safely
            int maxCol = a_index.Count > 0 ? a_index.Max() + 1 : 0;
            int[] colMap = new int[maxCol];
            int currentNewCol = 0;

            for (int c = 0; c < maxCol; c++)
            {
                colMap[c] = deleteSet.Contains(c) ? -1 : currentNewCol++;
            }

            // 3. Allocate clean tracking structures with estimated capacity bounds
            List<int> new_start = new(a_start.Count) {0};
            List<int> new_index = new(a_index.Count);
            List<double> new_value = new(a_value.Count);

            // 4. Process the entire matrix in a single linear sweep
            int rowCount = a_start.Count - 1;
            for (int i = 0; i < rowCount; i++)
            {
                int rowStart = a_start[i];
                int rowEnd = a_start[i + 1];

                for (int j = rowStart; j < rowEnd; j++)
                {
                    int globalCol = a_index[j];

                    // Check boundary to guard against tracking arrays larger than structural coordinates
                    if (globalCol >= 0 && globalCol < maxCol && colMap[globalCol] != -1)
                    {
                        new_index.Add(colMap[globalCol]);
                        new_value.Add(a_value[j]);
                    }
                }
                // The next row boundary is exactly the accumulated count of valid items
                new_start.Add(new_index.Count);
            }

            // FIX: Return the newly calculated compressed collections
            return (new_start, new_index, new_value);
        }

        public void Simulate2Phase(double[] ResultTime, List<Well<Vector2d>> Wells)
        {
            int Lx = Nx - 1, Ly = Ny - 1, Lz = Nz - 1;
            Direction xDir = Direction.X, yDir = Direction.Y, zDir = Direction.Z;
            bool Upstreamlock = false;
            int totalItems = Ngrids;
            int cores = Environment.ProcessorCount;
            int chunkSize = totalItems / cores;
            Trans = new double[Ngrids, 6];
            Water_Upstream = new bool[Ngrids, 6];
            Oil_Upstream = new bool[Ngrids, 6];

            void Residual(double[] xnp1, double time)
            {
                funcall++;

                //  HIGH CPU USAGE: Batched work gives cores long, uninterrupted calculations
                Parallel.For(0, cores, coreId =>
                {
                    int start = coreId * chunkSize;
                    int end = (coreId == cores - 1) ? totalItems : start + chunkSize;

                    // A single core executes this inner sequential loop at absolute maximum hardware speed
                    for (int m = start; m < end; m++)
                    {
                        if (Actnum[m] == 0) continue;


                        int indx1 = 2*m, indx2 = 2*m + 1;
                        ADiff<Vector2d> po_m = xnp1[2*m], sw_m = xnp1[2*m+1],
                        pw_m = po_m - Pc_I(sw_m), so_m = 1 - sw_m;
                        var meanP = so_m*po_m + sw_m*pw_m;
                        ErSo_Bo_np1[m] = Er(meanP)*so_m/Bo(po_m);
                        ErSw_Bw_np1[m] = Er(meanP)*sw_m/Bw(pw_m);
                        Res[indx1] = -PoreVolume[m]*(ErSo_Bo_np1[m] - ErSo_Bo_n[m]);
                        Res[indx2] = PoreVolume[m]*(ErSw_Bw_np1[m] - ErSw_Bw_n[m]);

                        int k = m / NxNy, rem = m % NxNy,
                        j = rem / Nx, i = rem % Nx;

                        // Add fluxes to the residual for this grid block from each of its 6 neighbors
                        if (i > 0)
                            AddFluxes(Xminus, m, m-1);
                        if (i < Lx)
                            AddFluxes(Xplus, m, m+1);
                        if (j > 0)
                            AddFluxes(Yminus, m, m-Nx);
                        if (j < Ly)
                            AddFluxes(Yplus, m, m+Nx);
                        if (k > 0)
                            AddFluxes(Zminus, m, m-NxNy);
                        if (k < Lz)
                            AddFluxes(Zplus, m, m+NxNy);
                        // Add aquifer fluxes if an aquifer is present and connected to this grid block
                        if (Aquifer != null && Aquifer.IsthereAquiferFlow(i, j, k))
                            AddAquiferFluxes(Aquifer.FlowDirection, m, Aquifer.Pa, Aquifer.Connectivity_Efficiency);
                    }
                });

                for (int nwell = 0; nwell < Nwells; nwell++)
                {
                    var well = Wells[nwell];
                    pwell = xnp1[2*Ngrids + 2*nwell];
                    qwell = xnp1[2*Ngrids + 2*nwell + 1];
                    Res[2*Ngrids + 2*nwell] = qwell*dt;
                    Res[2 * Ngrids + 2 * nwell + 1] = well.Constraint(time, pwell, qwell);
                    well.Qw = 0; well.Qo = 0; Zref = well.Zref;
                    ADiff water_rate = 0, oil_rate = 0;
                    switch (well.WellType)
                    {
                        case WellType.Producer:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                if (Actnum[m] == 0) continue;
                                int indx1 = 2*m, indx2 = 2*m+1;
                                ADiff po_m = xnp1[indx1], sw_m = xnp1[indx2],
                                pw_m = po_m - Pc_I(sw_m), so_m = 1 - sw_m,
                                go_avg = 0.5*(γo(po_m) + γo(pwell)),
                                gw_avg = 0.5*(γw(pw_m) + γw(pwell)),
                                DΦo = pwell - po_m - go_avg*(Zref - Z[m]),
                                DΦw = pwell - pw_m - gw_avg*(Zref - Z[m]);
                                WI = well.Perforation_WI[well.Perforation_NatIndex.IndexOf(m)];

                                WIo = WI*Kro(so_m)/(μo(po_m)*Bo(po_m));
                                oil_rate = DΦo*WIo;
                                well.Qo += oil_rate.Value;
                                Res[indx1] += oil_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= oil_rate*dt;

                                WIw = WI*Krw(sw_m)/(μw(pw_m)*Bw(pw_m));
                                water_rate = DΦw*WIw;
                                well.Qw += water_rate.Value;
                                Res[indx2] += water_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= water_rate*dt;
                            }
                            break;

                        case WellType.INJECTOR:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                if (Actnum[m] == 0) continue;
                                int indx1 = 2*m, indx2 = 2*m+1;
                                ADiff po_m = xnp1[indx1], sw_m = xnp1[indx2],
                                pw_m = po_m - Pc_I(sw_m), so_m = 1 - sw_m,
                                gw_avg = 0.5*(γw(pw_m) + γw(pwell)),
                                DΦw = pwell - pw_m - gw_avg*(Zref - Z[m]);
                                WI = well.Perforation_WI[well.Perforation_NatIndex.IndexOf(m)];

                                WIw = WI*krw0/(μw(pwell)*Bw(pwell));
                                water_rate = DΦw*WIw;
                                well.Qw += water_rate.Value;
                                Res[indx2] += water_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= water_rate*dt;
                            }
                            break;
                    }
                }
                b.Clear(); a_value.Clear();
                a_index.Clear(); a_start.Clear();
                a_start.Add(0);
                foreach (var res in Res.Where(r=> r is not null))
                {
                    b.Add(res.Value);
                    var sdic = res.Derivatives.OrderBy(kvp => kvp.Key);
                    a_value.AddRange(sdic.Select(kvp => kvp.Value));
                    a_index.AddRange(sdic.Select(kvp => kvp.Key));
                    a_start.Add(a_value.Count);
                }
                if(Columns2Delete.Count > 0)
                    (a_start, a_index, a_value) = DeleteCol(a_start, a_index, a_value, Columns2Delete);
            }

            dt = 0.1;
            // Initialize historical data tracking containers for plotting and reporting
            P = [Po_n]; S = [Sw_n]; Rate = [Qwells_n];
            Pwf = [Pwells_n]; WaterCut = [new double[Nwells]];
            Time = [0.0]; SweepEff = [0.0];
            List<double> Interval = [0, .. ResultTime];
            foreach (var well in Wells)
                Interval.AddRange(well.ProfileTime);
            Interval = [.. Interval.Distinct().OrderBy(x => x)];
            Console.WriteLine($"""
                    Total Grids = {Ngrids}, 
                    Active Grids = {Ngrids - Columns2Delete.Count/2}
                    ======================================================================
                                            Starting simulation
                    """);
            double L = 0, tnp1, history_dtmax = 0, rnorm = 0;
            bool staterejected, isComplete;
            for (int i = 1; i < Interval.Count; i++)
            {
                L = Interval[i];
                dt = Min(dt, Abs(L - t));
                isComplete = false;
                while (!isComplete)
                {
                    staterejected = false;
                    tnp1 = t + dt;
                    InitialGuess();

                    // Call the Newton-Raphson nonlinear solver to
                    // find the next state solution
                    Console.WriteLine($"""



                    Time: 
                    {tnp1:E4} days
                        iter  |   Residual Norm  
                    ----------+------------------
                    """);
                    int iter; bool isConverged = false;
                    for (iter = 1; iter < 10; iter++)
                    {
                        //Upstreamlock = iter > 4;
                        // Solve the nonlinear system using Newton-Raphson method
                        Residual(xs, tnp1); rnorm = b.Max(Abs);
                        double[] dx = LinSolve(a_value, a_index, a_start, b);

                        int activecount = 0;
                        for (int v = 0; v < 2*Ngrids; v++)
                            if (Actnum[v/2] == 1) 
                                xs[v].Value -= dx[activecount++];

                        for (int v = 2*Ngrids; v < varNum; v++)
                            xs[v].Value -= dx[activecount++];
                        
                        Console.WriteLine($"  {iter,4}    |    {rnorm:E4}");
                        isConverged = rnorm < 1e-6;
                        if (isConverged) break;
                    }

                    // Check convergence. If non-converged,
                    // chop the time step (time-step cuts) and retry.
                    if (!isConverged)
                    {
                        dt = 0.01*dt;
                        Console.WriteLine("""
                                ============================
                                 Rejected (Non-Convergence)
                                ============================
                                """);
                        continue;
                    }

                    // Unpack solution values to evaluate operational constraint validations
                    var (Po, Sw, Pwell, Qwell, maxDP, maxDS, i_p, i_s) = ExtractSolution();
                    if (maxDP > 400 || maxDS > 0.15)
                    {
                        int j = maxDP > 400 ? i_p : i_s;
                        double limit = maxDP > 400 ? 400 : 0.15;
                        double change = maxDP > 400 ? maxDP : maxDS;
                        Console.ForegroundColor = ConsoleColor.Yellow;
                        Console.WriteLine($"\n[Warning] Time Step Rejected: State Change Limit Exceeded.");
                        Console.ResetColor();
                        Console.WriteLine($"    ↳ Location       : Cell Index {j}");
                        Console.WriteLine($"    ↳ Violation      : {(maxDP > 400 ? "ΔP" : "ΔS")} = {change:F4} (Max Allowed: {limit:F4})");
                        Console.WriteLine($"    ↳ Current State  : P = {Po_n[j]:F1} psi | Sw = {Sw_n[j]:F4}");
                        Console.WriteLine($"    ↳ Attempted State: P = {Po[j]:F1} psi | Sw = {Sw[j]:F4}");
                        Console.WriteLine($"    ↳ Solver Action  : Rolling back to t = {t:F3} days.");
                        Console.WriteLine($"    ↳ Time-Step Chop : Reducing Δt from {dt:F2} to {0.5*dt:F2} days (Factor: 0.5).\n");
                        dt = 0.5*dt;
                        continue;
                    }

                    // Track the maximum time step that still yields convergence for potential future use in adaptive time-stepping
                    history_dtmax = Max(dt, history_dtmax);
                    for (int n = 0; n < Nwells; n++)
                    {
                        // Validate Producer constraints:
                        if (Wells[n].WellType == WellType.Producer)
                        {
                            // switch to BHP control if pressure falls below minimum limits
                            if (Wells[n].ConstraintType != ConstraintType.MinPressure &&
                                Pwells_n[n] < Wells[n].MinPressure)
                            {
                                Wells[n].ConstraintType = ConstraintType.MinPressure;
                                Console.WriteLine($"""
                                ===============================================================
                                      Rejected (Minimum Pressure Violated) @ {Wells[n].Name}
                                ===============================================================
                                """);
                                staterejected = true;
                            }
                        }

                        // Validate Injector constraints:
                        if (Wells[n].WellType == WellType.INJECTOR)
                        {

                            // switch to BHP control if pressure exceeds fracturing limits
                            if (Wells[n].ConstraintType != ConstraintType.MaxPressure &&
                                Pwells_n[n] > Wells[n].MaxPressure)
                            {
                                Wells[n].ConstraintType = ConstraintType.MaxPressure;
                                Console.WriteLine("""
                                ================================================
                                      Rejected (Maximum Pressure Violated) 
                                ================================================
                                """);
                                staterejected = true;
                            }
                        }
                    }
                    if (staterejected) continue;
                    

                    // Log verified parameters to performance history arrays
                    Po_n = Po; Sw_n = Sw; Pwells_n = Pwell; Qwells_n = Qwell;
                    Time.Add(t = tnp1); isComplete = Abs(t - L) < 1e-13;
                    P.Add(Po_n); S.Add(Sw_n); Rate.Add(Qwells_n); Pwf.Add(Pwells_n);
                    WaterCut.Add([.. Wells.Select(w => w.Fw)]);
                    SweepEff.Add((Sw_n.Sum() - S[0].Sum())*100/(Ngrids - S[0].Sum()));
                    for (int m = 0; m < Ngrids; m++)
                    {
                        if (Actnum[m] == 1)
                        {
                            ErSo_Bo_n[m] = ErSo_Bo_np1[m].Value;
                            ErSw_Bw_n[m] = ErSw_Bw_np1[m].Value;
                        }
                    }

                    // Adaptive Time-Stepping Logic:
                    // scale dt up if convergence is fast, scale down if slow
                    if (iter < 4) dt = 1.25*dt;
                    if (iter > 8) dt = 0.5*dt;
                    if (dt < 1e-12) throw new Exception("time step is too small");

                    // Prevent Overshoot
                    if (!isComplete) dt = Min(dt, L - t);
                    dt = Min(dt, 100);
                }
                if (history_dtmax > 0) dt = 0.1*history_dtmax;
            }
        }

        private double[] xCoords, yCoords, zCoords;
        public void ExportParaViewVTR(string folderPath)
        {
            xCoords = CalculateCoordinates([.. Enumerable.Range(0, Nx).Select(i => Dx[i])]);
            yCoords = CalculateCoordinates([.. Enumerable.Range(0, Ny).Select(i => Dy[i*Nx])]);
            zCoords = CalculateCoordinates([.. Enumerable.Range(0, Nz).Select(i => Dz[i*Nx*Ny])]);
            ExportSimulation(folderPath, "Test1", Nx, Ny, Nz,
                xCoords, yCoords, zCoords, Time, P, S);
        }

        /// <summary>
        /// Exports a single-block simulation time-series tracking Saturation at node corners.
        /// Perfectly mirrors ParaView's native parallel XML layout (PVD + PVTR + VTR).
        /// </summary>
        public static void ExportSimulation(
                string outputDirectory,
                string caseName,
                int Nx, int Ny, int Nz,
                double[] xCoord,        // Length: Nx + 1 (Nodes)
                double[] yCoord,        // Length: Ny + 1 (Nodes)
                double[] zCoord,        // Length: Nz + 1 (Nodes)
                List<double> Time,      // Length: Total timesteps
                List<double[]> P,       // List of flat cell-centered arrays (length: Nx * Ny * Nz)
                List<double[]> S)       // List of flat cell-centered arrays (length: Nx * Ny * Nz)
        {
            Directory.CreateDirectory(outputDirectory);
            DirectoryInfo directory = new DirectoryInfo(outputDirectory);
            foreach (FileInfo file in directory.EnumerateFiles())
                file.Delete();

            string fFormat = "G17"; // High-precision float string format
            var pvdSb = new StringBuilder();

            // 1. GENERATE THE TIMELINE COLLECTION MASTER (.pvd)
            pvdSb.AppendLine("<?xml version=\"1.0\"?>");
            pvdSb.AppendLine("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">");
            pvdSb.AppendLine("  <Collection>");

            int Nt = (int)(Time.Last()/20);
            double delt = Time.Last()/Nt;
            for (int framenum = 0; framenum <= 1000; framenum++)
            {
                double t = framenum*delt;
                double currentTime = t;
                double[] cellS = interpa(Time, S, t);
                double[] cellP = interpa(Time, P, t);

                // Unique filename for each timestep's grid data in the same directory
                string vtrFileName = $"{caseName}Time{framenum}.vtr";

                // PVD now links directly to the standalone VTR file
                pvdSb.AppendLine($"    <DataSet timestep=\"{currentTime.ToString(fFormat)}\" group=\"\" part=\"0\" file=\"{vtrFileName}\"/>");

                // -----------------------------------------------------------------
                // 2. GENERATE THE DYNAMIC DATA SNAPSHOT (.vtr)
                // -----------------------------------------------------------------
                var vtrSb = new StringBuilder();
                vtrSb.AppendLine("<VTKFile type=\"RectilinearGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">");
                vtrSb.AppendLine($"  <RectilinearGrid WholeExtent=\"0 {Nx} 0 {Ny} 0 {Nz}\">");

                // Optional time value metadata
                vtrSb.AppendLine("    <FieldData>");
                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"ascii\">");
                vtrSb.AppendLine($"        {currentTime.ToString(fFormat)}");
                vtrSb.AppendLine("      </DataArray>");
                vtrSb.AppendLine("    </FieldData>");

                vtrSb.AppendLine($"    <Piece Extent=\"0 {Nx} 0 {Ny} 0 {Nz}\">");

                // CellData Block maps simulation variables directly to the centers of your grid blocks
                vtrSb.AppendLine("    <CellData Scalars=\"Saturation,Pressure\">");

                // Saturation Array
                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Saturation\" format=\"ascii\" RangeMin=\"{cellS.Min()}\" RangeMax=\"{cellS.Max()}\">");
                vtrSb.Append("        ");
                int cellIndx = 0;
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            vtrSb.Append(cellS[cellIndx++].ToString(fFormat) + " ");
                        }
                        vtrSb.Append("\n        ");
                    }
                    vtrSb.Append("\n\n        ");
                }
                vtrSb.AppendLine("\n      </DataArray>\n");

                // Pressure Array
                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\" RangeMin=\"{cellP.Min()}\" RangeMax=\"{cellP.Max()}\">");
                vtrSb.Append("        ");
                cellIndx = 0;
                for (int k = 0; k < Nz; k++)
                {
                    for (int j = 0; j < Ny; j++)
                    {
                        for (int i = 0; i < Nx; i++)
                        {
                            vtrSb.Append(cellP[cellIndx++].ToString(fFormat) + " ");
                        }
                        vtrSb.Append("\n        ");
                    }
                    vtrSb.Append("\n\n        ");
                }
                vtrSb.AppendLine("\n      </DataArray>");
                vtrSb.AppendLine("    </CellData>");

                // Coordinates Block defines the physical bounding nodes of your grid blocks
                vtrSb.AppendLine("    <Coordinates>");

                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"X\" format=\"ascii\" RangeMin=\"{xCoord.Min()}\" RangeMax=\"{xCoord.Max()}\">");
                vtrSb.Append("        ");
                for (int i = 0; i <= Nx; i++) vtrSb.Append(xCoord[i].ToString(fFormat) + " ");
                vtrSb.AppendLine("\n      </DataArray>");

                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Y\" format=\"ascii\" RangeMin=\"{yCoord.Min()}\" RangeMax=\"{yCoord.Max()}\">");
                vtrSb.Append("        ");
                for (int j = 0; j <= Ny; j++) vtrSb.Append(yCoord[j].ToString(fFormat) + " ");
                vtrSb.AppendLine("\n      </DataArray>");

                vtrSb.AppendLine($"      <DataArray type=\"Float64\" Name=\"Z\" format=\"ascii\" RangeMin=\"{zCoord.Min()}\" RangeMax=\"{zCoord.Max()}\">");
                vtrSb.Append("        ");
                for (int k = 0; k <= Nz; k++) vtrSb.Append(zCoord[k].ToString(fFormat) + " ");
                vtrSb.AppendLine("\n      </DataArray>");

                vtrSb.AppendLine("    </Coordinates>");

                vtrSb.AppendLine("    </Piece>");
                vtrSb.AppendLine("  </RectilinearGrid>");
                vtrSb.AppendLine("</VTKFile>");

                File.WriteAllText(Path.Combine(outputDirectory, vtrFileName), vtrSb.ToString());
            }

            pvdSb.AppendLine("  </Collection>");
            pvdSb.AppendLine("</VTKFile>");
            File.WriteAllText(Path.Combine(outputDirectory, $"{caseName}.pvd"), pvdSb.ToString());
        }

        public void ExportWells(string outputDirectory)
        {
            Directory.CreateDirectory(outputDirectory);
            foreach (var well in Wells)
            {
                double[][] Trajectory = well.GetTrajectoryCoordinates(xCoords, yCoords, zCoords);
                ExportWellAsVtk(Path.Combine(outputDirectory, $"{well.Name}.vtk"), well.Name, Trajectory);
            }
        }

        /// <summary>
        /// Computes cumulative edge locations from cell delta sizes.
        /// Length of returning array is cell count + 1.
        /// </summary>
        private double[] CalculateCoordinates(double[] deltas)
        {
            double[] coords = new double[deltas.Length + 1];
            coords[0] = 0.0;
            for (int i = 0; i < deltas.Length; i++)
            {
                coords[i + 1] = coords[i] + deltas[i];
            }
            return coords;
        }

        public void ExportWellAsVtk(string outputPath, string wellName, double[][] trajectories)
        {
            // trajectories array expected shape: [numPoints][3] -> X, Y, Z coordinates
            int numPoints = trajectories.Length;

            using StreamWriter writer = new(outputPath, false, Encoding.ASCII);

            // 1. Write standard VTK Polydata Header
            writer.WriteLine("# vtk DataFile Version 3.0");
            writer.WriteLine($"Well Path: {wellName}");
            writer.WriteLine("ASCII");
            writer.WriteLine("DATASET POLYDATA");

            // 2. Write Geometry Nodes (Points)
            writer.WriteLine($"POINTS {numPoints} double");
            for (int i = 0; i < numPoints; i++)
            {
                writer.WriteLine($"{trajectories[i][0]:F2} {trajectories[i][1]:F2} {trajectories[i][2]:F2}");
            }
            writer.WriteLine();

            // 3. Topology: Link the points into a single continuous wireframe string
            // Format: LINES [NumberOfLines] [TotalIntegerCountInBlock]
            // Total count = NumberOfLines + total number of indices mapped
            writer.WriteLine($"LINES 1 {numPoints + 1}");
            writer.Write($"{numPoints} ");
            for (int i = 0; i < numPoints; i++)
            {
                writer.Write($"{i} ");
            }
            writer.WriteLine();
            writer.WriteLine();

            // 4. Optional Metadata Attributions (Coloring / Labeling in ParaView)
            writer.WriteLine($"CELL_DATA 1");
            writer.WriteLine("FIELD FieldData 1");
            writer.WriteLine($"WellID 1 {wellName.Length} string");
            writer.WriteLine(wellName);
        }


        // 1. Generates the individual active mesh snapshot for a specific time step
        public void ExportStepToVtu(
            string directoryPath, string baseName, int stepIndex,
            int nx, int ny, int nz, double[] dx, double[] dy, double[] dz, double[] zTop,
            int[] Actnum, double[] pressure, double[] saturation)
        {
            // Build the specific filename for this step (e.g., RunTest1_0.vtu)
            string fileName = $"{baseName}_{stepIndex}.vtu";
            string fullPath = Path.Combine(directoryPath, fileName);

            // 2. Compute explicit structural coordinates for vertices
            // A regular mesh has (nx+1)*(ny+1)*(nz+1) absolute grid vertex coordinates
            double[] xCoords = new double[nx + 1], yCoords = new double[ny + 1];

            for (int i = 1; i <= nx; i++) xCoords[i] = xCoords[i - 1] + dx[i - 1];
            for (int j = 1; j <= ny; j++) yCoords[j] = yCoords[j - 1] + dy[j - 1];

            StringBuilder xmlPoints = new(), xmlConnectivity = new(), xmlOffsets = new();
            StringBuilder xmlTypes = new(), xmlPressure = new(), xmlSaturation = new();

            int pointOffset = 0;
            int Ngrids = pressure.Length;
            int numActiveCells = Actnum.Where(a => a > 0).Count();
            int numPoints = numActiveCells*8; 
            for (int m = 0; m < Ngrids; m++)
            {
                if (Actnum[m] == 0) continue;
                Grid grid = Grids[m];

                for (int c = 0; c < 8; c++)
                    xmlPoints.AppendLine($"{grid.Corners[c].X} {grid.Corners[c].Y} {grid.Corners[c].Z}");     

                // Hexahedron connectivity mapping string index pointers
                xmlConnectivity.AppendLine($"{pointOffset} {pointOffset+1} {pointOffset+2} {pointOffset+3} {pointOffset+4} {pointOffset+5} {pointOffset+6} {pointOffset+7}");

                pointOffset += 8;
                xmlOffsets.AppendLine($"{pointOffset}");
                xmlTypes.AppendLine("12"); // VTK_HEXAHEDRON structural code identifier

                // Map data properties directly using the corresponding active global identifier
                xmlPressure.Append($"{pressure[m]} ");
                xmlSaturation.Append($"{saturation[m]} ");
            }

            // 3. Assemble complete XML output file buffer
            StringBuilder vtuFile = new();
            vtuFile.AppendLine("<?xml version=\"1.0\"?>");
            vtuFile.AppendLine("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">");
            vtuFile.AppendLine($"  <UnstructuredGrid>");
            vtuFile.AppendLine($"    <Piece NumberOfPoints=\"{numPoints}\" NumberOfCells=\"{numActiveCells}\">");

            // Point Data Section
            vtuFile.AppendLine("      <Points>");
            vtuFile.AppendLine("        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">");
            vtuFile.Append(xmlPoints.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("      </Points>");

            // Cell Connectivity Specification Section
            vtuFile.AppendLine("      <Cells>");
            vtuFile.AppendLine("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">");
            vtuFile.Append(xmlConnectivity.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">");
            vtuFile.Append(xmlOffsets.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">");
            vtuFile.Append(xmlTypes.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("      </Cells>");

            // Simulation Physical Property Field Section
            vtuFile.AppendLine("      <CellData Scalars=\"Pressure\">");
            vtuFile.AppendLine("        <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">");
            vtuFile.AppendLine(xmlPressure.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("        <DataArray type=\"Float64\" Name=\"Water_Saturation\" format=\"ascii\">");
            vtuFile.AppendLine(xmlSaturation.ToString());
            vtuFile.AppendLine("        </DataArray>");
            vtuFile.AppendLine("      </CellData>");

            vtuFile.AppendLine("    </Piece>");
            vtuFile.AppendLine("  </UnstructuredGrid>");
            vtuFile.AppendLine("</VTKFile>");

            File.WriteAllText(fullPath, vtuFile.ToString());
        }

        // 2. Generates the master pointer file that ParaView opens
        public void ExportParaViewVTU(string outputDirectory)
        {
            Directory.CreateDirectory(outputDirectory);
            DirectoryInfo directory = new DirectoryInfo(outputDirectory);
            foreach (FileInfo file in directory.EnumerateFiles())
                file.Delete();

            string baseName = "Test1";
            StringBuilder pvd = new();
            pvd.AppendLine("<?xml version=\"1.0\"?>");
            pvd.AppendLine("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">");
            pvd.AppendLine("  <Collection>");

            
            int Nt = (int)(Time.Last()/20);
            double delt = Time.Last()/Nt;
            for (int framenum = 0; framenum <= Nt; framenum++)
            {
                double t = framenum*delt;
                double currentTime = t;
                double[] cellS = interpa(Time, S, t);
                double[] cellP = interpa(Time, P, t);
                ExportStepToVtu(outputDirectory, baseName, framenum, Nx, Ny, Nz, Dx, Dy, Dz, zTop, Actnum, cellP, cellS);
                // Direct reference matching your exact naming convention
                string vtuFileName = $"{baseName}_{framenum}.vtu";
                pvd.AppendLine($"    <DataSet timestep=\"{currentTime}\" group=\"\" part=\"0\" file=\"{vtuFileName}\"/>");
            }

            pvd.AppendLine("  </Collection>");
            pvd.AppendLine("</VTKFile>");

            string pvdPath = Path.Combine(outputDirectory, $"{baseName}.pvd");
            File.WriteAllText(pvdPath, pvd.ToString());
        }

        public static double GetZcornValue(double[] zcorn, int nx, int ny, int nz,
                                   int i, int j, int k, int corner)
        {
            // Deconstruct the 0-7 corner index into local binary offsets
            int localX = (corner == 1 || corner == 3 || corner == 5 || corner == 6) ? 1 : 0;
            int localY = (corner == 2 || corner == 3 || corner == 6 || corner == 7) ? 1 : 0;
            int localZ = (corner >= 4) ? 1 : 0;

            // Standard ECLIPSE/OPM flat array offset formula
            int index = (k * 2 * ny * 2 * nx) +
                        (localZ * ny * 2 * nx) +
                        (j * 2 * nx) +
                        (localY * nx * 2) +
                        (i * 2) +
                        localX;

            return zcorn[index];
        }

    }
}
