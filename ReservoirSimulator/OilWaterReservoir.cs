using SepalSolver;
using System.Text;
using System.Timers;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class OilWaterReservoir
    {
        public int funcall;
        public List<double[]> P, S, Rate, OilRate, WaterRate, GasRate, Pwf, WaterCut;
        public List<double> Time, SweepEff;
        readonly Func<ADiff, ADiff> Sws, Swe, Pc_D, Pc_I, Bo, Bw, μo, μw, γo, γw, Kro, Krw, Er;
        // Define conversion constants
        const double alpha = 1.127e-3,        // Darcy to Field units factor
                beta = 5.615;                 // ft3 to bbl conversion factor

        double Transmissibility(Direction d, int m, int n)
        {
            return d switch
            {
                Direction.X => alpha * Dy[m] * Dz[m] * Harmmean(Kx[m] / Dx[m], Kx[n] / Dx[n]),
                Direction.Y => alpha * Dx[m] * Dz[m] * Harmmean(Ky[m] / Dy[m], Ky[n] / Dy[n]),
                Direction.Z => alpha * Dx[m] * Dy[m] * Harmmean(Kz[m] / Dz[m], Kz[n] / Dz[n]),
                _ => throw new ArgumentException("Invalid direction"),
            };
        }

        ADiff[] ErSo_Bo_np1, ErSw_Bw_np1, Res, xs;

        List<double> a_value = [], b = [];
        List<int> a_index = [], a_start = [0];

        static double betweenab(double a, double b, double f) => a + f*(b-a);
        static double interps(List<double> X, List<double> Y, double x)
        {
            int i = 0, j = 1;
            while (j < X.Count && X[j] < x) (i, j) = (j, j + 1);
            double f = (x - X[i-1])/(X[i] - X[i-1]);
            return betweenab(Y[i-1], Y[i], f);
        }
        static double[] interpa(List<double> X, List<double[]> Y, double x)
        {
            int i = 0, j = 1;
            while (j < X.Count && X[j] < x) (i, j) = (j, j + 1);
            double f = (x - X[i-1])/(X[i] - X[i-1]);
            return [.. Y[i-1].Zip(Y[i], (a, b) => betweenab(a, b, f))];
        }
        double Harmmean(double x1, double x2) => 2/(1/x1 + 1/x2);
        (double[] Po, double[] Sw, double[] Pwell, double[] Qwell, double maxDP, double maxDS, int i_p, int i_s) ExtractSolution()
        {
            int indx = 0, i_p = 0, i_s = 0; double maxDp = 0, maxDs = 0;
            double[] Po = [.. Enumerable.Repeat(double.NaN, Ngrids)],
                     Sw = [.. Enumerable.Repeat(double.NaN, Ngrids)],
                     Pwell = new double[Nwells], Qwell = new double[Nwells];
            for (int i = 0; i < Ngrids; i++)
            {
                if (Actnum[i] == 0) { indx += 2; continue; }
                Po[i] = xs[indx++].Value; // Matches Po index
                Sw[i] = xs[indx++].Value; // Matches Sw index
                if (Abs(Po[i] - Po_n[i]) > maxDp)
                { i_p = i; maxDp = Abs(Po[i] - Po_n[i]); }
                if (Abs(Sw[i] - Sw_n[i]) > maxDs)
                { i_s = i; maxDs = Abs(Sw[i] - Sw_n[i]); }
            }
            for (int i = 0; i < Nwells; i++)
            {
                Pwell[i] = xs[indx++].Value; // Matches Pwf index
                Qwell[i] = xs[indx++].Value; // Matches Q index
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
                xs[indx++].Value = Po_n[i];                          // Matches Po index
                if (n >= 2) xs[indx] += omega*(Po_n[i] - P[n-1][i]); // Extrapolation
                xs[indx++].Value = Sw_n[i];                          // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                xs[indx++].Value = Pwells_n[i];                        // Matches Pwf index
                xs[indx++].Value = Wells[i].ConstraintType switch
                {
                    ConstraintType.MaxPressure => Qwells_n[i],
                    ConstraintType.MinPressure => Qwells_n[i],
                    _ => Wells[i].ComputeRate(dt + Time[n])
                }; // Matches Q index

            }
        }
        double[] LinSolve(List<double> a_value, List<int> a_index, List<int> a_start, List<double> b)
        {
            return MklSparseSolver.Solve([.. a_value], [.. a_index], [.. a_start], [.. b]);
        }

        readonly int Nx, Ny, Nz, NxNy, Ngrids, Nwells, varNum;
        readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z;
        double[] Po_n, Sw_n, Pw_n, So_n, Qwells_n, Pwells_n, ErSw_Bw_n, ErSo_Bo_n;
        readonly double krw0, kro0, Pb, Prefw, Prefr, Pe, Pw_woc, Sw_r, So_r,
            Bo0, Bw0, Po_woc, Z_woc, co, cw, cr, bo, bw, nw, no, np,
            μo0, μw0, γo0, γw0, P_datum, Z_datum, Pc_max;
        List<Well> Wells;
        Aquifer Aquifer;
        double[,] Trans;
        bool[,] Water_Upstream, Oil_Upstream;
        int[] Actnum;
        List<int> Columns2Delete = [];
        public OilWaterReservoir(
            // DIMENS
            int _nx, int _ny, int _nz,

            // GRID
            double[] _dx, double[] _dy, double[] _dz, double[] _zTop,
            double[] _perm, double[] _phi, double _mult_z,

            // PVTW
            double _pref_w, double _bw0, double _cw, double _μw0, double _bw,

            // PVDO (Pressure, FVF, Muo) 
            Matrix _pvdo,

            // ROCK
            double _pref_r, double _cr,

            // SWOF (Sw, Krw, Kro, Pcwo)
            Matrix _swof,

            // DENSITY (Oil density, water density) 
            double _ρo0, double _ρw0,

            // EQUIL     
            double _datum, double _pdatun, double _z_woc, double _pcwoc,

            // WELL
            List<Well> _wells,

            // AQUIFER
            Aquifer _aquifer = null,

            // ACTNUM
            int[] _actnum = null)
        {
            ADiff.capacity = 16;
            Nx = _nx; Ny = _ny; Nz = _nz; NxNy = Nx*Ny; Ngrids = Nx*Ny*Nz;
            Dx = _dx; Dy = _dy; Dz = _dz; Kx = _perm; Ky = _perm;
            Kz = [.. _perm.Select(k => k*_mult_z)]; Φ = _phi;
            Wells = _wells; Aquifer = _aquifer;

            // Extract the Z
            Z = new double[Ngrids];
            if (_zTop.Length != NxNy)
                throw new Exception($"Number of elements in TOPS must be equal to {NxNy}");
            for (int i = 0; i < NxNy; i++)
                Z[i] = _zTop[i] + 0.5*Dz[i];
            for (int i = NxNy; i < Ngrids; i++)
                Z[i] = Z[i-NxNy] + 0.5*(Dz[i-NxNy] + Dz[i]);

            // Extract PVTW
            Bw0 = _bw0; μw0 = _μw0; γw0 = _ρw0/144; bw = _bw; cw = _cw;
            Bw = Pw => Bw0*Exp(-cw*(Pw - Prefw));
            μw = Pw => μw0*Exp(bw*(Pw - Prefw));
            γw = Pw => γw0*Exp(cw*(Pw - Prefw));

            // Extract PVDO
            Pb = double.PositiveInfinity;
            foreach (var well in Wells)
                Pb  = Min(Pb, well.MinPressure);

            PVDOAnalysis pvdo = new(_pvdo);
            var pvdoresult = pvdo.FitPvdoProperties(Pb);
            Bo0 = pvdoresult.BoAtPRef;
            μo0 = pvdoresult.ViscosityAtPRef;
            γo0 = _ρo0/144;
            co = pvdoresult.ConstantCompressibility;
            bo = pvdoresult.ViscosityExponent;
            Bo = Po => Bo0*Exp(-co*(Po - Pb));
            μo = Po => μo0*Exp(bo*(Po - Pb));
            γo = Po => γo0*Exp(co*(Po - Pb));


            // Extract SWOF
            SWOFAnalysis swof = new(_swof);
            var pcresult = swof.FitBrooksCoreyPc();
            var krresult = swof.FitCoreyPermeability();
            Pe = pcresult.Pd; np = pcresult.Lambda; Sw_r = pcresult.Swi;
            So_r = krresult.Sor; kro0 = krresult.KroEndpoint;
            krw0 = krresult.KrwEndpoint; no = krresult.No; nw = krresult.Nw;
            Pc_max = 150; 

            Sws = Sw => (Sw - Sw_r)/(1 - Sw_r);
            Swe = Sw => (Sw - Sw_r)/(1 - Sw_r - So_r);
            Pc_D = Sw => Pe <= 1e-12 ? 0.0 :
                    Sw <= Sw_r ? Pc_max :
                    Sw >= (1.0 - So_r) ? Pe :
                    Pe * Pow(Sws(Sw), -1.0 / np);
            Pc_I = Sw => Pe <= 1e-12 ? 0.0 :
                    Sw <= Sw_r ? Pc_max :
                    Sw >= 1.0 ? 0.0 :
                    Pe * (Pow(Swe(Sw), -1.0 / np) - 1.0);
            Kro = So => So <= So_r ? 0 : kro0 * Pow(1 - Swe(1 - So), no);
            Krw = Sw => Sw <= Sw_r ? 0 : krw0 * Pow(Swe(Sw), nw);

            //Extract ROCK
            cr = _cr; Prefr = _pref_r;
            Er = P => Exp(cr*(P - Prefr));

            // Extract datum
            Z_woc = _z_woc; P_datum = _pdatun; Z_datum = _datum;
            if (_datum < _z_woc)
            {
                Po_woc = P_datum + γo(P_datum).Value*(Z_woc - Z_datum);
                Pw_woc = Po_woc - Pe;
            }
            else
            {
                Pw_woc = P_datum + γw(P_datum).Value*(Z_woc - Z_datum);
                Po_woc = Pw_woc + Pe;
            }
            Wells = _wells; Nwells = Wells.Count; varNum = 2*Ngrids + 2*Nwells;
            Aquifer = _aquifer;
            Actnum = _actnum is not null? _actnum:[.. Enumerable.Repeat(1, Ngrids)];
        }
        public OilWaterReservoir(int _nx, int _ny, int _nz,
            double[] _perm, double[] _phi, double[] _dx, double[] _dy,
            double[] _dz, double[] _z, double _peow, double _pw_woc,
            double _z_woc, double _mult_z, double _sw_r, double _so_r,
            double _bo0, double _bw0, double _μo0, double _μw0, double _γo0,
            double _γw0, double _krw0, double _kro0, double _co, double _cw,
            double _cr, double _bo, double _bw, double _nw, double _no,
            double _pb, double _pref, List<Well> _wells, Aquifer _aquifer = null, int[] _actnum = null)
        {
            ADiff.capacity = 16;
            Kx = _perm; Ky = _perm; Kz = [.. _perm.Select(k => k*_mult_z)];
            Nx = _nx; Ny = _ny; Nz = _nz; NxNy = Nx*Ny; Ngrids = Nx*Ny*Nz;
            Dx = _dx; Dy = _dy; Dz = _dz; Z = _z; Φ = _phi; kro0 = _kro0;
            krw0 = _krw0; Bo0 = _bo0; Bw0 = _bw0; Pb = _pb; Prefw = _pref;
            Prefr = _pref; Pe = _peow; So_r = _so_r; Sw_r = _sw_r; co = _co;
            cw = _cw; cr = _cr; bo = _bo; bw = _bw; no = _no; nw = _nw;
            np = 0.7; μo0 = _μo0; μw0 = _μw0; γo0 = _γo0; γw0 = _γw0;
            Pw_woc = _pw_woc; Po_woc = Pw_woc + Pe; Z_woc = _z_woc;
            Wells = _wells; Nwells = Wells.Count; varNum = 2*Ngrids + 2*Nwells;
            Aquifer = _aquifer;
            Actnum = _actnum is not null ? _actnum : [.. Enumerable.Repeat(1, Ngrids)];


            Sws = Sw => (Sw - Sw_r)/(1 - Sw_r);
            Swe = Sw => (Sw - Sw_r)/(1 - Sw_r - So_r);
            Pc_D = Sw => Pe <= 1e-12 ? 0.0 :
                Sw <= Sw_r ? Pc_max :
                Sw >= (1.0 - So_r) ? Pe :
                Pe * Pow(Sws(Sw), -1.0 / np);
            Pc_I = Sw => Pe <= 1e-12 ? 0.0 :
                Sw <= Sw_r ? Pc_max :
                Sw >= 1.0 ? 0.0 :
                Pe * (Pow(Swe(Sw), -1.0 / np) - 1.0);
            Bo = Po => Bo0*Exp(-co*(Po - Pb));
            Bw = Pw => Bw0*Exp(-cw*(Pw - Prefw));
            μo = Po => μo0*Exp(bo*(Po - Pb));
            μw = Pw => μw0*Exp(bw*(Pw - Prefw));
            γo = Po => γo0*Exp(co*(Po - Pb));
            γw = Pw => γw0*Exp(cw*(Pw - Prefw));
            Er = P => Exp(cr*(P - Prefr));
            Krw = Sw => Sw <= Sw_r ? 0 : Sw > 1 - So_r ? krw0 : krw0 * Pow(Swe(Sw), nw);
            Kro = So => So <= So_r ? 0 : So > 1 - Sw_r ? kro0 : kro0 * Pow(1 - Swe(1 - So), no);
        }

        public void Initialize()
        {
            funcall = 0;

            // 2. Initialize the spatial grid blocks
            Pw_n = [..Enumerable.Repeat(double.NaN, Ngrids)]; 
            Sw_n = [.. Enumerable.Repeat(double.NaN, Ngrids)];
            Po_n = [.. Enumerable.Repeat(double.NaN, Ngrids)]; 
            So_n = [.. Enumerable.Repeat(double.NaN, Ngrids)];
            ErSo_Bo_n = [.. Enumerable.Repeat(double.NaN, Ngrids)]; 
            ErSw_Bw_n = [.. Enumerable.Repeat(double.NaN, Ngrids)];
            ErSo_Bo_np1 = new ADiff[Ngrids]; ErSw_Bw_np1 = new ADiff[Ngrids];
            Res = new ADiff[varNum]; xs = new ADiff[varNum];

            if (Pe == 0)
            {
                for (int i = 0; i < Ngrids; i++)
                {
                    if (Actnum[i] == 0)
                    {
                        Columns2Delete.Add(2*i);
                        Columns2Delete.Add(2*i + 1);
                        continue;
                    }

                    if (Z[i] < Z_woc)
                    {
                        double p = Po_woc + γo(Po_woc).Value * (Z[i] - Z_woc);
                        Po_n[i] = Po_woc + 0.5*(γo(p) + γo(Po_woc)).Value * (Z[i] - Z_woc);
                        Pw_n[i] = Po_n[i]; Sw_n[i] = Sw_r;
                    }
                    else
                    {
                        double p = Pw_woc + γw(Pw_woc).Value * (Z[i] - Z_woc);
                        Pw_n[i] = Pw_woc + 0.5*(γw(p) + γw(Pw_woc)).Value * (Z[i] - Z_woc);
                        Po_n[i] = Pw_n[i]; Sw_n[i] = 1;
                    }

                    So_n[i] = 1.0 - Sw_n[i];
                    xs[2*i] = new(Po_n[i], 2*i);
                    xs[2*i+1] = new(Sw_n[i], 2*i+1);
                    Res[2*i] = new(); Res[2*i+1] = new();
                    double meanP = So_n[i]*Po_n[i] + Sw_n[i]*Pw_n[i];
                    ErSo_Bo_n[i] = (Er(meanP) * So_n[i] / Bo(Po_n[i])).Value;
                    ErSw_Bw_n[i] = (Er(meanP) * Sw_n[i] / Bw(Pw_n[i])).Value;
                    ErSo_Bo_np1[i] = new(); ErSw_Bw_np1[i] = new();
                }
            }
            else
            {
                // 1. Pre-generate a fine lookup table for the inverse Capillary Pressure relationship
                List<double> Sw_Table = [.. Linspace(1.0-So_r-1e-5, Sw_r+1e-5, 50)];
                // Calculate Pc for each Sw point in our table
                List<double> Pc_Table = [.. Sw_Table.Select(s => Pc_D(s).Value)];

                for (int i = 0; i < Ngrids; i++)
                {
                    if (Actnum[i] == 0)
                    {
                        Columns2Delete.Add(2*i);
                        Columns2Delete.Add(2*i + 1);
                        continue;
                    }

                    Pw_n[i] = Pw_woc + γw(Pw_woc).Value * (Z[i] - Z_woc);
                    Po_n[i] = Po_woc + γo(Po_woc).Value * (Z[i] - Z_woc);
                    double pc = Po_n[i] - Pw_n[i];

                    // 3. Directly interpolate saturation instead of using an iterative solver
                    if (pc > Pc_Table.First() && pc <= Pc_Table.Last())
                        // interps(List<double> X, List<double> Y, double x)
                        Sw_n[i] = interps(Pc_Table, Sw_Table, pc);
                    else if (pc > Pc_Table.Last())
                        // Capillary pressure exceeds our table limit; clamp to residual water
                        Sw_n[i] = Sw_r + 1e-5;
                    else
                        // Below or at the entry boundary threshold
                        Sw_n[i] = 1.0;

                    So_n[i] = 1.0 - Sw_n[i];
                    xs[2*i] = new(Po_n[i], 2*i);
                    xs[2*i+1] = new(Sw_n[i], 2*i+1);
                    Res[2*i] = new(); Res[2*i+1] = new();
                    double meanP = So_n[i]*Po_n[i] + Sw_n[i]*Pw_n[i];
                    ErSo_Bo_n[i] = (Er(meanP) * So_n[i] / Bo(Po_n[i])).Value;
                    ErSw_Bw_n[i] = (Er(meanP) * Sw_n[i] / Bw(Pw_n[i])).Value;
                    ErSo_Bo_np1[i] = new(); ErSw_Bw_np1[i] = new();
                }
            }

            Pwells_n = new double[Nwells]; Qwells_n = new double[Nwells];
            for (int i = 0; i < Nwells; i++)
            {
                Qwells_n[i] = 0; // Start with no flow
                Wells[i].ComputeNaturalIndex(Nx, Ny);
                Wells[i].ComputeProductivityIndex(Kx, Ky, Dx, Dy, Dz);
                switch (Wells[i].WellType)
                {
                    case WellType.Producer:
                        Pwells_n[i] = Po_n[Wells[i].Perforation_NatIndex.First()];
                        Wells[i].Zref = Z[Wells[i].Perforation_NatIndex.First()];
                        break;
                    case WellType.Injector:
                        Pwells_n[i] = Pw_n[Wells[i].Perforation_NatIndex.Last()];
                        Wells[i].Zref = Z[Wells[i].Perforation_NatIndex.Last()];
                        break;
                }
                xs[2*Ngrids + 2*i] = new(Pwells_n[i], 2*Ngrids + 2*i);
                xs[2*Ngrids + 2*i + 1] = new(Qwells_n[i], 2*Ngrids + 2*i + 1);
                Res[2*Ngrids + 2*i] = new(); Res[2*Ngrids + 2*i+1] = new();
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

        public void Simulate2Phase(double[] ResultTime, List<Well> Wells)
        {
            int Lx = Nx - 1, Ly = Ny - 1, Lz = Nz - 1;
            Direction xDir = Direction.X, yDir = Direction.Y, zDir = Direction.Z;
            FlowDirection Xminus = FlowDirection.Iminus, Xplus = FlowDirection.Iplus,
                          Yminus = FlowDirection.Jminus, Yplus = FlowDirection.Jplus,
                          Zminus = FlowDirection.Kminus, Zplus = FlowDirection.Kplus;
            bool Upstreamlock = false;
            int totalItems = Ngrids;
            int cores = Environment.ProcessorCount;
            int chunkSize = totalItems / cores;
            Trans = new double[Ngrids, 6];
            Water_Upstream = new bool[Ngrids, 6];
            Oil_Upstream = new bool[Ngrids, 6];

            Parallel.For(0, cores, coreId =>
            {
                int start = coreId * chunkSize;
                int end = (coreId == cores - 1) ? totalItems : start + chunkSize;

                // A single core executes this inner sequential loop at absolute maximum hardware speed
                for (int m = start; m < end; m++)
                {
                    int k = m / NxNy, rem = m % NxNy,
                    j = rem / Nx, i = rem % Nx;

                    // Add fluxes to the residual for this grid block from each of its 6 neighbors
                    if (i > 0) Trans[m, (int)Xminus] = Transmissibility(xDir, m, m - 1);
                    if (i < Lx) Trans[m, (int)Xplus] = Transmissibility(xDir, m, m + 1);
                    if (j > 0) Trans[m, (int)Yminus] = Transmissibility(yDir, m, m - Nx);
                    if (j < Ly) Trans[m, (int)Yplus] = Transmissibility(yDir, m, m + Nx);
                    if (k > 0) Trans[m, (int)Zminus] = Transmissibility(zDir, m, m - NxNy);
                    if (k < Lz) Trans[m, (int)Zplus] = Transmissibility(zDir, m, m + NxNy);
                }
            });

            void CheckUpstream(FlowDirection Flowdir, int m, ADiff OilPot, ADiff WaterPot)
            {
                Oil_Upstream[m, (int)Flowdir] = OilPot > 0;
                Water_Upstream[m, (int)Flowdir] = WaterPot > 0;
            }

            void Residual(ADiff[] xnp1, double time)
            {
                funcall++;
                double WI, Zref; ADiff WIw, WIo, pwell, qwell;

                void AddFluxes(FlowDirection Flowdir, int m, int n)
                {
                    if (Actnum[n] == 0) return;
                    ADiff Po_up, Pw_up, So_up, Sw_up, Tw, To;
                    double Tr = Trans[m, (int)Flowdir];
                    int indx1, indx2;
                    ADiff po_n = xnp1[indx1 = 2*n],
                        sw_n = xnp1[indx2 = 2*n+1],
                        pw_n = po_n - Pc_I(sw_n),
                        so_n = 1 - sw_n,
                        po_m = xnp1[indx1 = 2*m],
                        sw_m = xnp1[indx2 = 2*m+1],
                        pw_m = po_m - Pc_I(sw_m),
                        so_m = 1 - sw_m,
                        go_avg = 0.5*(γo(po_m) + γo(po_n)),
                        gw_avg = 0.5*(γw(pw_m) + γw(pw_n)),
                        DΦo = po_n - po_m - go_avg * (Z[n] - Z[m]),
                        DΦw = pw_n - pw_m - gw_avg * (Z[n] - Z[m]);

                    if (!Upstreamlock) CheckUpstream(Flowdir, m, DΦo, DΦw);

                    (Po_up, So_up) = Oil_Upstream[m, (int)Flowdir] ? (po_n, so_n) : (po_m, so_m);
                    To = Tr*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                    Res[indx1] += To* DΦo * dt;

                    (Pw_up, Sw_up) = Water_Upstream[m, (int)Flowdir] ? (pw_n, sw_n) : (pw_m, sw_m);
                    Tw = Tr*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                    Res[indx2] += Tw* DΦw * dt;
                }

                void AddAquiferFluxes(FlowDirection Flowdir, int m, double Pa, double Ceff)
                {
                    ADiff Tw;
                    double Tr = Flowdir switch
                    {
                        FlowDirection.Iminus => Transmissibility(xDir, m, m),
                        FlowDirection.Iplus => Transmissibility(xDir, m, m),
                        FlowDirection.Jminus => Transmissibility(yDir, m, m),
                        FlowDirection.Jplus => Transmissibility(yDir, m, m),
                        FlowDirection.Kminus => Transmissibility(zDir, m, m),
                        FlowDirection.Kplus => Transmissibility(zDir, m, m),
                        _ => 0
                    };
                    int indx1, indx2;
                    ADiff po_m = xnp1[indx1 = 2*m],
                        sw_m = xnp1[indx2 = 2*m+1],
                        pw_m = po_m - Pc_I(sw_m), DP = Pa - pw_m;
                    Tw = Tr*krw0/(μw(pw_m)*Bw(pw_m));
                    Res[indx2] += Ceff * Tw * dt * Flowdir switch
                    {
                        FlowDirection.Kminus => DP + 0.5*γw(pw_m)*Dz[m],
                        FlowDirection.Kplus => DP - 0.5*γw(pw_m)*Dz[m],
                        _ => DP
                    };
                }

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
                        double PV = -Dx[m]*Dy[m]*Dz[m]*Φ[m]/beta;
                        ADiff po_m = xnp1[2*m], sw_m = xnp1[2*m+1],
                            pw_m = po_m - Pc_I(sw_m),
                            so_m = 1 - sw_m;
                        var meanP = so_m*po_m + sw_m*pw_m;
                        ErSo_Bo_np1[m] = Er(meanP)*so_m/Bo(po_m);
                        ErSw_Bw_np1[m] = Er(meanP)*sw_m/Bw(pw_m);
                        Res[indx1] = PV*(ErSo_Bo_np1[m] - ErSo_Bo_n[m]);
                        Res[indx2] = PV*(ErSw_Bw_np1[m] - ErSw_Bw_n[m]);

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
                    well.WaterRate = 0; well.OilRate = 0; Zref = well.Zref;
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
                                well.OilRate += oil_rate.Value;
                                Res[indx1] += oil_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= oil_rate*dt;

                                WIw = WI*Krw(sw_m)/(μw(pw_m)*Bw(pw_m));
                                water_rate = DΦw*WIw;
                                well.WaterRate += water_rate.Value;
                                Res[indx2] += water_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= water_rate*dt;
                            }
                            break;

                        case WellType.Injector:
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
                                well.WaterRate += water_rate.Value;
                                Res[indx2] += water_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= water_rate*dt;
                            }
                            break;
                    }
                    Res[2 * Ngrids + 2 * nwell + 1] = well.Constraint(time, pwell, qwell);
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

            dt = 0.01;
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
                        if (Wells[n].WellType == WellType.Injector)
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
                    WaterCut.Add([.. Wells.Select(w => w.WaterCut)]);
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
        public void ExportParaView(string folderPath)
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
                List<double> time,      // Length: Total timesteps
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

            double delt = time.Last()/1000;

            for (int framenum = 0; framenum <= 1000; framenum++)
            {
                double t = framenum*delt;
                double currentTime = t;
                double[] cellS = interpa(time, S, t);
                double[] cellP = interpa(time, P, t);

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

    }
}
