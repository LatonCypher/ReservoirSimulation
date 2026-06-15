using SepalSolver;
using System.Text;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class OilWaterReservoir
    {
        public int funcall;
        public List<double[]> P, S, Rate, Pwf, WaterCut;
        public List<double> Time, SweepEff;
        readonly Func<ADiff, ADiff> Sws, Swe, Pc_D, Pc_I, Bo, Bw, μo, μw, γo, γw, Kro, Krw, Er;
        // Define conversion constants
        const double alpha = 1.127e-3,        // Darcy to Field units factor
          alpha_well = 1.127e-3*2*pi,         // Darcy to Field units factor for wells
                beta = 5.615;                 // ft3 to bbl conversion factor

        double[,] GridCoordinates;

        double OverlapArea(int m, int n)
        {
            return 0;
        }
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
        double betweenab(double a, double b, double f) => a + f*(b-a);
        double interps(List<double> X, List<double> Y, double x)
        {
            int i = X.FindIndex(xi => xi > x);
            double f = (x - X[i-1])/(X[i] - X[i-1]);
            return betweenab(Y[i-1], Y[i], f);
        }
        double[] interpa(List<double> X, List<double[]> Y, double x)
        {
            int i = X.FindIndex(xi => xi > x);
            double f = (x - X[i-1])/(X[i] - X[i-1]);
            return [.. Y[i-1].Zip(Y[i], (a, b) => betweenab(a, b, f))];
        }
        double Harmmean(double x1, double x2) => 2/(1/x1 + 1/x2);
        (double[], double[], double[], double[]) Unpack(ADiff[] x)
        {
            int indx = 0;
            double[] Po = new double[Ngrids], Sw = new double[Ngrids],
                     Pwell = new double[Nwells], Qwell = new double[Nwells];
            for (int i = 0; i < Ngrids; i++)
            {
                Po[i] = x[indx++].Value; // Matches Po index
                Sw[i] = x[indx++].Value; // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                Pwell[i] = x[indx++].Value; // Matches Pwf index
                Qwell[i] = x[indx++].Value; // Matches Q index
            }
            return (Po, Sw, Pwell, Qwell);
        }
        void Pack(double[] Po, double[] Sw, double[] Pwell, double[] Qwell, ADiff[] x)
        {
            int indx = 0;
            for (int i = 0; i < Ngrids; i++)
            {
                x[indx++].Value = Po[i];  // Matches Po index
                x[indx++].Value = Sw[i];  // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                x[indx++].Value = Pwell[i]; // Matches Pwf index
                x[indx++].Value = Qwell[i]; // Matches Q index
            }
        }

        public readonly int Nx, Ny, Nz, NxNy, Ngrids, Nwells, varNum;
        public readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z;
        double[] Po_n, Sw_n, Pw_n, So_n, Qwells_n, Pwells_n, ErSw_Bw_n, ErSo_Bo_n;
        readonly double krw0, kro0, Pb, Prefw, Prefr, Pe, Pw_woc, Sw_r, So_r,
            Bo0, Bw0, Po_woc, Z_woc, co, cw, cr, bo, bw, nw, no, np,
            μo0, μw0, γo0, γw0, P_datum, Z_datum, Pc_max;
        List<Well> Wells;

        public OilWaterReservoir(
            // DIMENS
            int _nx, int _ny, int _nz,

            // GRID
            double[] _dx, double[] _dy,  double[] _dz, double[] _zTop, 
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
            List<Well> _wells)
        {
            ADiff.capacity = 16;
            Nx = _nx; Ny = _ny; Nz = _nz; NxNy = Nx*Ny; Ngrids = Nx*Ny*Nz;
            Dx = _dx; Dy = _dy; Dz = _dz; Kx = _perm; Ky = _perm; 
            Kz = [.. _perm.Select(k => k*_mult_z)]; Φ = _phi;
            Wells = _wells;

            // Extract the Z
            Z = new double[Ngrids];
            if (_zTop.Length != NxNy)
                throw new Exception($"Number of elements in TOPS must be equal to {NxNy}");
            for (int i = 0; i < NxNy; i++)
                Z[i]= _zTop[i] + 0.5*Dz[i];
            for (int i = NxNy; i < Ngrids; i++)
                Z[i]= Z[i-NxNy] + 0.5*(Dz[i-NxNy] + Dz[i]);

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
        }
        public OilWaterReservoir(int _nx, int _ny, int _nz,
            double[] _perm, double[] _phi, double[] _dx, double[] _dy,
            double[] _dz, double[] _z, double _peow, double _pw_woc,
            double _z_woc, double _mult_z, double _sw_r, double _so_r,
            double _bo0, double _bw0, double _μo0, double _μw0, double _γo0,
            double _γw0, double _krw0, double _kro0, double _co, double _cw,
            double _cr, double _bo, double _bw, double _nw, double _no,
            double _pb, double _pref, List<Well> _wells)
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
            GridCoordinates = new double[Ngrids, 24];


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
            Pw_n = new double[Ngrids]; Sw_n = new double[Ngrids];
            Po_n = new double[Ngrids]; So_n = new double[Ngrids];

            if (Pe == 0)
            {
                for (int i = 0; i < Ngrids; i++)
                {
                    Pw_n[i] = Pw_woc + (Z[i] > Z_woc ? γw(Pw_woc) : γo(Pw_woc)).Value * (Z[i] - Z_woc);
                    Po_n[i] = Pw_n[i]; Sw_n[i] = Z[i] > Z_woc ? 1 : Sw_r; So_n[i] = 1.0 - Sw_n[i];
                }
            }
            else
            {
                // 1. Pre-generate a fine lookup table for the inverse Capillary Pressure relationship
                List<double> Sw_Table = [.. Linspace(1.0-So_r-1e-5, Sw_r+1e-5, 50)];
                // Calculate Pc for each Sw point in our table
                List<double> Pc_Table = [.. Sw_Table.Select(sw => Pc_D(sw).Value)];

                for (int i = 0; i < Ngrids; i++)
                {
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
                }
            }

            Qwells_n = new double[Nwells]; Pwells_n = new double[Nwells];
            for (int i = 0; i < Nwells; i++)
            {
                Qwells_n[i] = 0; // Start with no flow
                Wells[i].ComputeNaturalIndex(Nx, Ny);
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
            }
            (ErSo_Bo_n, ErSw_Bw_n) = PoreVolumeFraction(Po_n, Pw_n, So_n, Sw_n);
        }

        (double[], double[]) PoreVolumeFraction(double[] Po, double[] Pw, double[] So, double[] Sw)
        {
            double[] Fo = new double[Ngrids], Fw = new double[Ngrids];
            for (int i = 0; i < Ngrids; i++)
            {
                double meanP = So[i]*Po[i] + Sw[i]*Pw[i];
                Fw[i] = (Er(meanP) * Sw[i] / Bw(Pw[i])).Value;
                Fo[i] = (Er(meanP) * So[i] / Bo(Po[i])).Value;
            }
            return (Fo, Fw);
        }

        public static (List<int>, List<int>, List<double>) DeleteCol(List<int> a_start, 
            List<int> a_index, List<double> a_value, List<int> index2delete)
        {
            if(index2delete.Count == 0)
                return (a_start, a_index, a_value);

            // 1. Build a global column remapping array (O(ColCount))
            int maxCol = a_index.Count > 0 ? a_index.Max() + 1 : 0;
            int[] colMap = new int[maxCol];
            int currentNewCol = 0;
            for (int c = 0; c<maxCol; c++)
                colMap[c] = index2delete.Contains(c) ? -1 : currentNewCol++;

            // 2. Allocate clean, high-performance tracking structures
            List<int> new_start = [0], new_index = [];
            List<double> new_value = [];

            // 3. Process the entire matrix in a single linear sweep (O(NNZ))
            int rowCount = a_start.Count - 1;
            for (int i = 0; i<rowCount; i++)
            {
                int rowStart = a_start[i];
                int rowEnd = a_start[i + 1];

                for (int j = rowStart; j<rowEnd; j++)
                {
                    int globalCol = a_index[j];
                    if (colMap[globalCol] != -1) // If it's a keeper, copy directly
                    {
                        new_index.Add(colMap[globalCol]);
                        new_value.Add(a_value[j]);
                    }
                }
                // The next row boundary is exactly the accumulated count of valid items
                new_start.Add(new_index.Count);
            }

            // 4. Update your master collections via O(1) reference swapping
            return (a_start, a_index, a_value);
        }


        public void Simulate2Phase(double[] ResultTime, List<Well> Wells)
        {
            int Lx = Nx - 1, Ly = Ny - 1, Lz = Nz - 1;
            double dt, t = 0;

            Direction xDir = Direction.X, yDir = Direction.Y, zDir = Direction.Z;

            List<double> a_value = [], b = [];
            List<int> a_index = [], a_start = [0];
            ADiff[] Res = new ADiff[varNum], xs;
            void Residual(ADiff[] xnp1, double time)
            {
                funcall++;
                double re, WI, Zref; ADiff WIw, WIo, pwell, qwell;

                void AddFluxes(Direction dir, int m, int n)
                {
                    ADiff Po_up, Pw_up, So_up, Sw_up, Tw, To;
                    double Tr = Transmissibility(dir, m, n); 
                    int indx1, indx2;
                    ADiff po_n = xnp1[indx1 = 2*n],
                        sw_n = xnp1[indx2 = 2*n+1],
                        pw_n = po_n - Pc_I(sw_n),
                        so_n = 1 - sw_n,
                        pmgoz_n = po_n - γo(po_n)*Z[n],
                        pmgwz_n = pw_n - γw(pw_n)*Z[n],
                        po_m = xnp1[indx1 = 2*m],
                        sw_m = xnp1[indx2 = 2*m+1],
                        pw_m = po_m - Pc_I(sw_m),
                        so_m = 1 - sw_m,
                        pmgoz_m = po_m - γo(po_m)*Z[m],
                        pmgwz_m = pw_m - γw(pw_m)*Z[m];

                    (Po_up, So_up) = pmgoz_m >= pmgoz_n ? (po_m, so_m) : (po_n, so_n);
                    To = Tr*Kro(So_up)/(μo(0.5*(po_m + po_n))*Bo(0.5*(po_m + po_n)));
                    //To = Tr*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                    Res[indx1] += To*(po_n - po_m - γo(Po_up)*(Z[n] - Z[m]))*dt;

                    (Pw_up, Sw_up) = pmgwz_m >= pmgwz_n ? (pw_m, sw_m) : (pw_n, sw_n);
                    Tw = Tr*Krw(Sw_up)/(μw(0.5*(pw_m + pw_n))*Bw(0.5*(pw_m + pw_n)));
                    //Tw = Tr*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                    Res[indx2] += Tw*(pw_n - pw_m - γw(Pw_up)*(Z[n] - Z[m]))*dt;
                }

                //  HIGH CPU USAGE: Batched work gives cores long, uninterrupted calculations
                int totalItems = Ngrids;
                int cores = Environment.ProcessorCount;
                int chunkSize = totalItems / cores;
                Parallel.For(0, cores, coreId =>
                {
                    int start = coreId * chunkSize;
                    int end = (coreId == cores - 1) ? totalItems : start + chunkSize;

                    // A single core executes this inner sequential loop at absolute maximum hardware speed
                    for (int m = start; m < end; m++)
                    {
                        int indx1 = 2*m, indx2 = 2*m + 1;
                        double V = Dx[m]*Dy[m]*Dz[m]/beta;
                        ADiff po_m = xnp1[2*m], sw_m = xnp1[2*m+1],
                            pw_m = po_m - Pc_I(sw_m),
                            so_m = 1 - sw_m;
                        var meanP = so_m*po_m + sw_m*pw_m;
                        var erso_bo_np1 = Er(meanP)*so_m/Bo(po_m);
                        var ersw_bw_np1 = Er(meanP)*sw_m/Bw(pw_m);
                        Res[indx1] = -V*Φ[m]*(erso_bo_np1 - ErSo_Bo_n[m]);
                        Res[indx2] = -V*Φ[m]*(ersw_bw_np1 - ErSw_Bw_n[m]);

                        int k = m / NxNy, rem = m % NxNy,
                        j = rem / Nx, i = rem % Nx;

                        if (i > 0) AddFluxes(xDir, m, m-1);
                        if (i < Lx) AddFluxes(xDir, m, m+1);
                        if (j > 0) AddFluxes(yDir, m, m-Nx);
                        if (j < Ly) AddFluxes(yDir, m, m+Nx);
                        if (k > 0) AddFluxes(zDir, m, m-NxNy);
                        if (k < Lz) AddFluxes(zDir, m, m+NxNy);
                    }
                });

                for (int nwell = 0; nwell < Nwells; nwell++)
                {
                    var well = Wells[nwell];
                    pwell = xnp1[2*Ngrids + 2*nwell];
                    qwell = xnp1[2*Ngrids + 2*nwell + 1];
                    Res[2*Ngrids + 2*nwell] = qwell*dt;
                    Res[2*Ngrids + 2*nwell + 1] = well.Constraint(time, pwell, qwell);
                    well.WaterRate = 0; well.OilRate = 0; Zref = well.Zref;
                    ADiff water_rate = 0, oil_rate = 0;

                    switch (well.WellType)
                    {
                        case WellType.Producer:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                int indx1 = 2*m, indx2 = 2*m+1;
                                ADiff po_m = xnp1[indx1], sw_m = xnp1[indx2],
                                pw_m = po_m - Pc_I(sw_m), so_m = 1 - sw_m;
                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);

                                WIo = WI*Kro(so_m)/(μo(po_m)*Bo(po_m));
                                oil_rate = (pwell - po_m - γo(po_m)*(Zref - Z[m]))*WIo;

                                well.OilRate += oil_rate.Value;
                                Res[indx1] += oil_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= oil_rate*dt;

                                WIw = WI*Krw(sw_m)/(μw(pw_m)*Bw(pw_m));
                                water_rate = (pwell - pw_m - γw(pw_m)*(Zref - Z[m]))*WIw;

                                well.WaterRate += water_rate.Value;
                                Res[indx2] += water_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= water_rate*dt;
                            }
                            break;

                        case WellType.Injector:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                int indx1 = 2*m, indx2 = 2*m+1;
                                ADiff po_m = xnp1[indx1], sw_m = xnp1[indx2],
                                pw_m = po_m - Pc_I(sw_m), so_m = 1 - sw_m;
                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);

                                WIw = WI*krw0/(μw(pwell)*Bw(pwell));
                                water_rate = (pwell - pw_m - γw(pwell)*(Zref - Z[m]))*WIw;

                                well.WaterRate += water_rate.Value;
                                Res[indx2] += water_rate*dt;
                                Res[2*Ngrids + 2*nwell] -= water_rate*dt;
                            }
                            break;
                    }
                }
                b.Clear(); a_value.Clear();
                a_index.Clear(); a_start.Clear();
                a_start.Add(0);
                foreach (var res in Res)
                {
                    b.Add(res.Value);
                    var sdic = res.Derivatives.OrderBy(kvp => kvp.Key);
                    a_value.AddRange(sdic.Select(kvp => kvp.Value));
                    a_index.AddRange(sdic.Select(kvp => kvp.Key));
                    a_start.Add(a_value.Count);
                }
                int[] index2delete = [];
                if (index2delete.Length!=0)
                {
                    for (int i = 1; i < a_start.Count; i++)
                    {
                        for (int j = a_start[i]; j --> a_start[i-1]; )
                        {
                            if (index2delete.Contains(a_index[j]))
                            {
                                a_index.RemoveAt(j);
                                a_value.RemoveAt(j);
                                for (int k = j; k < a_start[i]; k++)
                                    a_index[k] -= 1;
                                for (int k = i; k < a_start.Count; k++)
                                    a_start[k] -= 1;
                            }
                        }
                    }
                }
            }

            dt = 0.001;
            xs = [.. Enumerable.Range(0, varNum).Select(i => new ADiff(0, i))];
            Pack(Po_n, Sw_n, Pwells_n, Qwells_n, xs);
            // Initialize historical data tracking containers for plotting and reporting
            P = [Po_n]; S = [Sw_n]; Rate = [Qwells_n]; 
            Pwf = [Pwells_n]; WaterCut = [new double[Nwells]];
            Time = [0.0]; SweepEff = [0.0];
            List<double>  Interval = [0, .. ResultTime];
            foreach (var well in Wells)
                Interval.AddRange(well.ProductionProfile.Time);
            Interval = [.. Interval.Distinct().OrderBy(x => x)];
            Console.WriteLine($"""
                    ======================================================================
                                            Starting simulation
                    """);
            double L = 0, tnp1, history_dtmax = 0, rnorm = 0;
            bool staterejected, isComplete;
            // Track the history of residual norms for this specific timestep
            //List<(int, double)> residualHistory = new List<(int, double)>(15);
            //double steplength = 1;
            for (int i = 1; i < Interval.Count; i++)
            {
                L = Interval[i];
                dt = Min(dt, Abs(L - t));
                isComplete = false;
                while (!isComplete)
                {
                    staterejected = false;
                    tnp1 = t + dt;

                    // Call the Newton-Raphson nonlinear solver to
                    // find the next state solution
                    Console.WriteLine($"""
                    Time: 
                    {tnp1:F3} days
                        iter  |   Residual Norm  
                    ----------+----------------
                    """);
                    int iter; bool isConverged = false;
                    //residualHistory.Clear();
                    //steplength = 1;
                    for (iter = 1; iter < 10; iter++)
                    {
                        // Solve the nonlinear system using Newton-Raphson method
                        Residual(xs, tnp1); rnorm = b.Max(Abs);
                        double[] dx = MklSparseSolver.Solve([.. a_value], [.. a_index], [.. a_start], [.. b]);
                        for (int v = 0; v < varNum; v++) xs[v].Value -= dx[v];
                        Console.WriteLine($"  {iter,4}    |    {rnorm:F4}");
                        isConverged = rnorm < 1e-6;
                        if (isConverged) break;
                    }
                    Console.WriteLine("\n\n\n\n");

                    // Check convergence. If non-converged,
                    // chop the time step (time-step cuts) and retry.
                    if (!isConverged)
                    {
                        dt = 0.25*dt;
                        Console.WriteLine("""
                                ================================================
                                           Rejected (Non-Convergence)
                                ================================================
                                """);
                        Pack(Po_n, Sw_n, Pwells_n, Qwells_n, xs);
                        continue;
                    }

                    // Unpack solution values to evaluate operational constraint validations
                    double maxDP = 0, maxDS = 0;
                    for (int m = 0; m < Ngrids; m++)
                    {
                        maxDP = Max(maxDP, Abs(xs[2*m].Value - Po_n[m]));
                        maxDS = Max(maxDS, Abs(xs[2*m+1].Value - Sw_n[m]));
                    }
                    if (maxDP > 100 || maxDS > 0.15)
                    {
                        dt = 0.5*dt; 
                        Pack(Po_n, Sw_n, Pwells_n, Qwells_n, xs);
                        continue;
                    }
                    (Po_n, Sw_n, Pwells_n, Qwells_n) = Unpack(xs);
                    history_dtmax = Max(dt, history_dtmax);

                    for (int n = 0; n < Nwells; n++)
                    {
                        // Validate Producer constraints:
                        if (Wells[n].WellType == WellType.Producer)
                        {
                            // switch to BHP control if pressure falls below minimum limits
                            if (Wells[n].ConstraintType == ConstraintType.FlowRate &&
                                Pwells_n[n] < Wells[n].MinPressure)
                            {
                                Wells[n].ConstraintType = ConstraintType.MinPressure;
                                Console.WriteLine("""
                                ================================================
                                      Rejected (Minimum Pressure Violated) 
                                ================================================
                                """);
                                staterejected = true;
                                break;
                            }
                        }

                        // Validate Injector constraints:
                        if (Wells[n].WellType == WellType.Injector)
                        {

                            // switch to BHP control if pressure exceeds fracturing limits
                            if (Wells[n].ConstraintType == ConstraintType.FlowRate &&
                                Pwells_n[n] > Wells[n].MaxPressure)
                            {
                                Wells[n].ConstraintType = ConstraintType.MaxPressure;
                                Console.WriteLine("""
                                ================================================
                                      Rejected (Maximum Pressure Violated) 
                                ================================================
                                """);
                                staterejected = true;
                                break;
                            }
                        }
                    }
                    if (staterejected)
                    {
                        Pack(Po_n, Sw_n, Pwells_n, Qwells_n, xs);
                        continue;
                    }

                    // Log verified parameters to performance history arrays
                    Time.Add(t = tnp1); isComplete = Abs(t - L) < 1e-13;
                    P.Add(Po_n); S.Add(Sw_n); Rate.Add(Qwells_n); Pwf.Add(Pwells_n);
                    WaterCut.Add([.. Wells.Select(w => w.WaterCut)]);
                    SweepEff.Add((Sw_n.Sum() - S[0].Sum())*100/(Ngrids - S[0].Sum()));
                    (ErSo_Bo_n, ErSw_Bw_n) = PoreVolumeFraction(Po_n, Pw_n, So_n, Sw_n);

                    // Adaptive Time-Stepping Logic:
                    // scale dt up if convergence is fast, scale down if slow
                    if (iter < 4) dt = 1.25*dt;
                    if (iter > 8) dt = 0.5*dt;
                    if (dt < 1e-5) throw new Exception("time step is too small");

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

            for (int t = 0; t < time.Count; t++)
            {
                double currentTime = time[t];
                double[] cellS = S[t];
                double[] cellP = P[t];

                // Unique filename for each timestep's grid data in the same directory
                string vtrFileName = $"{caseName}Time{t}.vtr";

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

    }
}
