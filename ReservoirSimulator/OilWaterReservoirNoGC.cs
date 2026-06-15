using Microsoft.CodeAnalysis.CSharp.Syntax;
using SepalSolver;
using System.Text;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class OilWaterReservoirNoGC
    {
        public List<double[]> P, S, Rate, Pwf, WaterCut;
        public List<double> Time, SweepEff;
        public int funcall;

        // Define conversion constants
        const double alpha = 1.127e-3,        // Darcy to Field units factor
          alpha_well = 1.127e-3*2*pi,         // Darcy to Field units factor for wells
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
        List<int> a_index = [], a_start = [];
        
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
     

        /// <summary>
        /// Sws = Sw => (Sw - Sw_r)/(1 - Sw_r);
        /// </summary>
        /// <param name="Sw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Sws(ADiff Sw, ADiff result)
        {
            result.CopyFrom(Sw).
                SubtractInPlace(Sw_r).
                DivideInPlace(1 - Sw_r);  
            return result;
        }
        double Sws(double Sw)=>
            (Sw - Sw_r)/(1 - Sw_r);

        /// <summary>
        /// Swe = Sw => (Sw - Sw_r)/(1 - Sw_r - So_r);
        /// </summary>
        /// <param name="Sw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Swe(ADiff Sw, ADiff result)
        {
            result.CopyFrom(Sw).
                SubtractInPlace(Sw_r).
                DivideInPlace(1 - Sw_r - So_r);
            return result;
        }
        double Swe(double Sw) =>
            (Sw - Sw_r)/(1 - Sw_r - So_r);

        /// <summary>
        /// Pc_D = Sw => Pe * Pow(Sws(Sw), -1.5);
        /// </summary>
        /// <param name="Sw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Pc_D(ADiff Sw, ADiff result)
        {
            if (Pe < 1e-12)
                return result.Clear();
            if (Sw.Value >= 1.0)
                return result.Clear();
            Sws(Sw, result).PowInPlace(-1.5).
                MultiplyInPlace(Pe);
            return result;
        }
        double Pc_D(double Sw)=>
            Pe <= 1e-12 ? 0.0 :
            Sw >= (1.0 - So_r) ? Pe :
            Pe * Pow(Sws(Sw), -np);

        /// <summary>
        /// Pc_I = Sw => Pe * (Pow(Swe(Sw), -1.5) - 1);
        /// </summary>
        /// <param name="Sw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Pc_I(ADiff Sw, ADiff result) 
        {
            if (Pe < 1e-12)
                return result.Clear();
            if(Sw.Value >= 1.0)
                return result.Clear();
            Swe(Sw, result).PowInPlace(-1.5).
                SubtractInPlace(1).
                MultiplyInPlace(Pe);
            return result;
        }
        double Pc_I(double Sw)=>
            Pe <= 1e-12 ? 0.0 :
            Sw >= 1.0 ? 0.0 :
            Pe * (Pow(Swe(Sw), -np) - 1.0);

        /// <summary>
        /// Bo = Po => Bo0*Exp(-co*(Po - Pb));
        /// </summary>
        /// <param name="Po"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Bo(ADiff Po, ADiff result)
        {
            result.CopyFrom(Po).
                SubtractInPlace(Pb).
                MultiplyInPlace(-co).
                ExpInPlace().
                MultiplyInPlace(Bo0);
            return result;
        }
        double Bo(double Po) =>
            Bo0*Exp(-co*(Po - Pb));

        /// <summary>
        /// Bw = Pw => Bw0*Exp(-cw*(Pw - Pref));
        /// </summary>
        /// <param name="Pw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Bw(ADiff Pw, ADiff result)
        {
            result.CopyFrom(Pw).
                SubtractInPlace(Prefw).
                MultiplyInPlace(-cw).
                ExpInPlace().
                MultiplyInPlace(Bw0);
            return result;
        }
        double Bw(double Pw) =>
            Bw0*Exp(-cw*(Pw - Prefw));

        /// <summary>
        /// μo = Po => μo0*Exp(bo*(Po - Pb));
        /// </summary>
        /// <param name="Po"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff μo(ADiff Po, ADiff result)
        {
            result.CopyFrom(Po).
                SubtractInPlace(Pb).
                MultiplyInPlace(bo).
                ExpInPlace().
                MultiplyInPlace(μo0);
            return result;
        }
        double μo(double Po) =>
            μo0*Exp(bo*(Po - Pb));

        /// <summary>
        /// μw = Pw => μw0*Exp(bw*(Pw - Pref));
        /// </summary>
        /// <param name="Pw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff μw(ADiff Pw, ADiff result)
        {
            result.CopyFrom(Pw).
                SubtractInPlace(Prefw).
                MultiplyInPlace(bw).
                ExpInPlace().
                MultiplyInPlace(μw0);
            return result;
        }
        double μw(double Pw)=>
            μw0*Exp(bw*(Pw - Prefw));

        /// <summary>
        /// γo = Po => γo0*Exp(co*(Po - Pb));
        /// </summary>
        /// <param name="Po"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff γo(ADiff Po, ADiff result)
        {
            result.CopyFrom(Po).
                SubtractInPlace(Pb).
                MultiplyInPlace(co).
                ExpInPlace().
                MultiplyInPlace(γo0);
            return result;
        }
        double γo(double Po)=>
            γo0*Exp(co*(Po - Pb));

        /// <summary>
        /// γw = Pw => γw0*Exp(cw*(Pw - Pref));
        /// </summary>
        /// <param name="Pw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff γw(ADiff Pw, ADiff result)
        {
            result.CopyFrom(Pw).
                SubtractInPlace(Prefw).
                MultiplyInPlace(cw).
                ExpInPlace().
                MultiplyInPlace(γw0);
            return result;
        }
        double γw( double Pw) => 
            γw0*Exp(cw*(Pw - Prefw));

        /// <summary>
        /// Er = P => Exp(cr*(P - Pref));
        /// </summary>
        /// <param name="P"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Er(ADiff P, ADiff result)
        {
            result.CopyFrom(P).
                SubtractInPlace(Prefr).
                MultiplyInPlace(cr).
                ExpInPlace();
            return result;
        }
        double Er(double P) => 
            Exp(cr*(P - Prefr));

        /// <summary>
        /// Kro = So => kro0 * Pow(1 - Swe(1 - So), no);
        /// </summary>
        /// <param name="So"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Kro(ADiff So, ADiff result)
        {
            double so = So.Value;

            // Physical boundary check: If oil is at or below residual saturation, Kro is zero
            if (so <= So_r)
            {
                result.Value = 0.0;
                result.Derivatives.Clear();
                return result;
            }

            // Step 1: Compute effective oil saturation
            double denom = 1.0 - Sw_r - So_r;
            double so_e = (so - So_r) / denom;

            // Step 2: Assign the primal value directly to your destination scratchpad
            result.Value = kro0 * Pow(so_e, no);

            // Step 3: Compute the exact analytical chain-rule scale factor
            double scale = (no * result.Value) / (so - So_r);

            // Step 4: Clear destination scratchpad and map the scaled derivatives from the input 'So'
            result.Derivatives.Clear();
            foreach (var kvp in So.Derivatives)
                result.Derivatives[kvp.Key] = kvp.Value * scale;

            return result;
        }
        double Kro(double So) => 
            So <= So_r ? 0 : kro0 * Pow(1 - Swe(1 - So), no);

        /// <summary>
        /// Krw = Sw => krw0 * Pow(Swe(Sw), nw);
        /// </summary>
        /// <param name="Sw"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        ADiff Krw(ADiff Sw, ADiff result)
        {
            double sw = Sw.Value;

            // Physical boundary check: If water is at or below connate saturation, Krw is zero
            if (sw <= Sw_r)
            {
                result.Value = 0.0;
                result.Derivatives.Clear();
                return result;
            }

            // Step 1: Compute effective water saturation
            double denom = 1.0 - Sw_r - So_r;
            double sw_e = (sw - Sw_r) / denom;

            // Step 2: Assign the primal value directly to your destination scratchpad
            result.Value = krw0 * Pow(sw_e, nw);

            // Step 3: Compute the exact analytical chain-rule scale factor
            double scale = (nw * result.Value) / (sw - Sw_r);

            // Step 4: Clear destination scratchpad and map the scaled derivatives from the input 'Sw'
            result.Derivatives.Clear();
            foreach (var kvp in Sw.Derivatives)
                result.Derivatives[kvp.Key] = kvp.Value * scale;

            return result;
        }
        double Krw(double Sw) => 
            Sw <= Sw_r ? 0 : krw0 * Pow(Swe(Sw), nw);

        readonly int Nx, Ny, Nz, NxNy, Ngrids, Nwells, varNum;
        readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z;
        double[] Po_n, Sw_n, Pw_n, So_n, Qwells_n, Pwells_n, ErSw_Bw_n, ErSo_Bo_n;
        readonly double krw0, kro0, Pb, Prefw, Prefr, Pe, Pw_woc, Sw_r, So_r,
            Bo0, Bw0, Po_woc, Z_woc, co, cw, cr, bo, bw, nw, no, np,
            μo0, μw0, γo0, γw0, P_datum, Z_datum;
        List<Well> Wells;
        Matrix SWOF, SGOF, PVTO, PVDO, PVDG;

        public OilWaterReservoirNoGC(
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

            // Extract SWOF
            SWOFAnalysis swof = new(_swof);
            var pcresult = swof.FitBrooksCoreyPc();
            var krresult = swof.FitCoreyPermeability();
            Pe = pcresult.Pd; np = pcresult.Lambda; Sw_r = pcresult.Swi;
            So_r = krresult.Sor; kro0 = krresult.KroEndpoint;
            krw0 = krresult.KrwEndpoint; no = krresult.No; nw = krresult.Nw;

            //Extract ROCK
            cr = _cr; Prefr = _pref_r;

            // Extract datum
            Z_woc = _z_woc; P_datum = _pdatun; Z_datum = _datum;
            if (_datum < _z_woc)
            {
                Po_woc = P_datum + γo(P_datum)*(Z_woc - Z_datum);
                Pw_woc = Po_woc - Pe;
            }
            else
            {
                Pw_woc = P_datum + γw(P_datum)*(Z_woc - Z_datum);
                Po_woc = Pw_woc + Pe;
            }
            Wells = _wells; Nwells = Wells.Count; varNum = 2*Ngrids + 2*Nwells;
        }
        
        public OilWaterReservoirNoGC(int _nx, int _ny, int _nz,
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
            np = -1.5; μo0 = _μo0; μw0 = _μw0; γo0 = _γo0; γw0 = _γw0;
            Pw_woc = _pw_woc; Po_woc = Pw_woc + Pe; Z_woc = _z_woc;
            Wells = _wells; Nwells = Wells.Count; varNum = 2*Ngrids + 2*Nwells;
        }
        
        public void Initialize()
        {
            funcall = 0;
            // 2. Initialize the spatial grid blocks
            Pw_n = new double[Ngrids]; Sw_n = new double[Ngrids];
            Po_n = new double[Ngrids]; So_n = new double[Ngrids];
            xs = new ADiff[varNum];   Res = new ADiff[varNum];  
            ErSo_Bo_np1 = new ADiff[Ngrids]; ErSw_Bw_np1 = new ADiff[Ngrids];

            if (Pe == 0)
            {
                for (int i = 0; i < Ngrids; i++)
                {
                    Pw_n[i] = Pw_woc + (Z[i] > Z_woc ? γw(Pw_woc) : γo(Pw_woc)) * (Z[i] - Z_woc);
                    Po_n[i] = Pw_n[i]; Sw_n[i] = Z[i] > Z_woc ? 1 : Sw_r; So_n[i] = 1.0 - Sw_n[i];
                    xs[2*i] = new(Po_n[i], 2*i); xs[2*i+1] = new(Sw_n[i], 2*i+1);
                    ErSo_Bo_np1[i] = new(); ErSw_Bw_np1[i] = new(); 
                    Res[2*i] = new(); Res[2*i+1] = new();
                }
            }
            else
            {
                // 1. Pre-generate a fine lookup table for the inverse Capillary Pressure relationship
                List<double> Sw_Table = [.. Linspace(1.0-So_r-1e-5, Sw_r+1e-5, 50)];
                // Calculate Pc for each Sw point in our table
                List<double> Pc_Table = [.. Sw_Table.Select(Pc_D)];

                // Initialize the ADiff scratchpad arrays for the next time step
                ErSo_Bo_np1 = new ADiff[Ngrids]; ErSw_Bw_np1 = new ADiff[Ngrids];

                Res = new ADiff[varNum]; xs = new ADiff[varNum];

                for (int i = 0; i < Ngrids; i++)
                {
                    Pw_n[i] = Pw_woc + γw(Pw_woc) * (Z[i] - Z_woc);
                    Po_n[i] = Po_woc + γo(Po_woc) * (Z[i] - Z_woc);
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
                    xs[2*i] = new(Po_n[i], 2*i); xs[2*i+1] = new(Sw_n[i], 2*i+1);
                    Res[2*i] = new(); Res[2*i+1] = new();
                    ErSo_Bo_np1[i] = new(); ErSw_Bw_np1[i] = new();
                }
            }

            Qwells_n = new double[Nwells];  Pwells_n = new double[Nwells];
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
                xs[2*Ngrids + 2*i] = new(Pwells_n[i], 2*Ngrids + 2*i);
                xs[2*Ngrids + 2*i + 1] = new(Qwells_n[i], 2*Ngrids + 2*i + 1); 
                Res[2*Ngrids + 2*i] = new(); Res[2*Ngrids + 2*i+1] = new();
            }
            (ErSo_Bo_n, ErSw_Bw_n) = PoreVolumeFraction(Po_n, Pw_n, So_n, Sw_n);
        }

        void GuessReset()
        {
            int indx = 0;
            for (int i = 0; i < Ngrids; i++)
            {
                xs[indx++].Value = Po_n[i];  // Matches Po index
                xs[indx++].Value = Sw_n[i];  // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                xs[indx++].Value = Pwells_n[i]; // Matches Pwf index
                xs[indx++].Value = Qwells_n[i]; // Matches Q index
            }
        }

        void ExtractSolution()
        {
            int indx = 0;
            for (int i = 0; i < Ngrids; i++)
            {
                Po_n[i] = xs[indx++].Value;  // Matches Po index
                Sw_n[i] = xs[indx++].Value;  // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                Pwells_n[i] = xs[indx++].Value; // Matches Pwf index
                Qwells_n[i] = xs[indx++].Value; // Matches Q index
            }
        }

        (double[], double[]) PoreVolumeFraction(double[] Po, double[] Pw, double[] So, double[] Sw)
        {
            double[] Fo = new double[Ngrids], Fw = new double[Ngrids];
            for (int i = 0; i < Ngrids; i++)
            {
                double meanP = So[i]*Po[i] + Sw[i]*Pw[i];
                Fw[i] = Er(meanP) * Sw[i] / Bw(Pw[i]);
                Fo[i] = Er(meanP) * So[i] / Bo(Po[i]);
            }
            return (Fo, Fw);
        }

        void PoreVolumeUpdate()
        {
            for (int i = 0; i < Ngrids; i++)
            {
                ErSo_Bo_n[i] = ErSo_Bo_np1[i].Value;
                ErSw_Bw_n[i] = ErSw_Bw_np1[i].Value;
            }
        }

        static ADiff GetPo(ADiff[] x, int i, ADiff Po) =>
            Po.AddInPlace(x[2*i]);

        static ADiff GetSw(ADiff[] x, int i, ADiff Sw) =>
            Sw.AddInPlace(x[2*i+1]);

        static ADiff GetSo(ADiff Sw, ADiff So) =>
            So.AddInPlace(Sw).NegateInPlace().AddInPlace(1);

        ADiff GetPw(ADiff Po, ADiff Sw, ADiff Pw) =>
            Pc_I(Sw, Pw).NegateInPlace().AddInPlace(Po);

        public static (List<int>, List<int>, List<double>) DeleteCol(List<int> a_start,
            List<int> a_index, List<double> a_value, List<int> index2delete)
        {
            if (index2delete.Count == 0)
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
            ADiff pwell = new(), qwell = new(), tempvar = new();

            void Residual(ADiff[] xnp1, double time)
            {
                funcall++;
                double re, WI, Zref; 
                ScratchPad Wellscratchpad = new();

                void ComputeBlockPVT(ScratchPad scratchpad, int m)
                {
                    scratchpad.ClearMVariables();
                    GetPo(xnp1, m, scratchpad.Po_m);
                    GetSw(xnp1, m, scratchpad.Sw_m);
                    GetSo(scratchpad.Sw_m, scratchpad.So_m);
                    GetPw(scratchpad.Po_m, scratchpad.Sw_m, scratchpad.Pw_m);
                    Bo(scratchpad.Po_m, scratchpad.Bo_m);
                    Bw(scratchpad.Pw_m, scratchpad.Bw_m);
                    γo(scratchpad.Po_m, scratchpad.Go_m);
                    γw(scratchpad.Pw_m, scratchpad.Gw_m);
                    μo(scratchpad.Po_m, scratchpad.Uo_m);
                    μw(scratchpad.Pw_m, scratchpad.Uw_m);
                    Kro(scratchpad.So_m, scratchpad.Kro_m);
                    Krw(scratchpad.Sw_m, scratchpad.Krw_m);
                    scratchpad.PmGZo_m.
                        SubtractInPlace(scratchpad.Go_m).
                        MultiplyInPlace(Z[m]).
                        AddInPlace(scratchpad.Po_m);
                    scratchpad.PmGZw_m.
                        SubtractInPlace(scratchpad.Gw_m).
                        MultiplyInPlace(Z[m]).
                        AddInPlace(scratchpad.Pw_m);
                }

                void AddFluxes(Direction dir, int m, int n, ScratchPad scratchpad)
                {
                    scratchpad.ClearNVariables();
                    GetPo(xnp1, n, scratchpad.Po_n);
                    GetSw(xnp1, n, scratchpad.Sw_n);
                    GetSo(scratchpad.Sw_n, scratchpad.So_n);
                    GetPw(scratchpad.Po_n, scratchpad.Sw_n, scratchpad.Pw_n);
                    Bo(scratchpad.Po_n, scratchpad.Bo_n);
                    Bw(scratchpad.Pw_n, scratchpad.Bw_n);
                    γo(scratchpad.Po_n, scratchpad.Go_n);
                    γw(scratchpad.Pw_n, scratchpad.Gw_n);
                    μo(scratchpad.Po_n, scratchpad.Uo_n);
                    μw(scratchpad.Pw_n, scratchpad.Uw_n);
                    Kro(scratchpad.So_n, scratchpad.Kro_n);
                    Krw(scratchpad.Sw_n, scratchpad.Krw_n);
                    scratchpad.PmGZo_n.
                        SubtractInPlace(scratchpad.Go_n).
                        MultiplyInPlace(Z[n]).
                        AddInPlace(scratchpad.Po_n);
                    scratchpad.PmGZw_n.
                        SubtractInPlace(scratchpad.Gw_n).
                        MultiplyInPlace(Z[n]).
                        AddInPlace(scratchpad.Pw_n);

                    ADiff Tr = scratchpad.Tr.AddInPlace(Transmissibility(dir, m, n)),
                          po_m = scratchpad.Po_m, bo_m = scratchpad.Bo_m, go_m = scratchpad.Go_m,
                          μo_m = scratchpad.Uo_m, kro_m = scratchpad.Kro_m, pmgzo_m = scratchpad.PmGZo_m,

                          pw_m = scratchpad.Pw_m, bw_m = scratchpad.Bw_m, gw_m = scratchpad.Gw_m,
                          μw_m = scratchpad.Uw_m, krw_m = scratchpad.Krw_m, pmgzw_m = scratchpad.PmGZw_m,

                          po_n = scratchpad.Po_n, bo_n = scratchpad.Bo_n, go_n = scratchpad.Go_n,
                          μo_n = scratchpad.Uo_n, kro_n = scratchpad.Kro_n, pmgzo_n = scratchpad.PmGZo_n,

                          pw_n = scratchpad.Pw_n, bw_n = scratchpad.Bw_n, gw_n = scratchpad.Gw_n,
                          μw_n = scratchpad.Uw_n, krw_n = scratchpad.Krw_n, pmgzw_n = scratchpad.PmGZw_n;


                    if (pmgzo_m >= pmgzo_n)
                    {
                        scratchpad.Oil_Rate.
                            CopyFrom(po_n).
                            SubtractInPlace(po_m).
                            SubtractProductInPlace(go_m, Z[n] - Z[m]).
                            MultiplyInPlace(Tr).
                            MultiplyInPlace(kro_m).
                            DivideInPlace(bo_m).
                            DivideInPlace(μo_m);
                    }
                    else
                    {
                        scratchpad.Oil_Rate.
                            CopyFrom(po_n).
                            SubtractInPlace(po_m).
                            SubtractProductInPlace(go_n, Z[n] - Z[m]).
                            MultiplyInPlace(Tr).
                            MultiplyInPlace(kro_n).
                            DivideInPlace(bo_n).
                            DivideInPlace(μo_n);
                    }
                    Res[2*m].AddProductInPlace(scratchpad.Oil_Rate, dt);


                    if (pmgzw_m >= pmgzw_n)
                    {
                        scratchpad.Water_Rate.
                            CopyFrom(pw_n).
                            SubtractInPlace(pw_m).
                            SubtractProductInPlace(gw_m, Z[n] - Z[m]).
                            MultiplyInPlace(Tr).
                            MultiplyInPlace(krw_m).
                            DivideInPlace(bw_m).
                            DivideInPlace(μw_m);
                    }
                    else
                    {
                        scratchpad.Water_Rate.
                            CopyFrom(pw_n).
                            SubtractInPlace(pw_m).
                            SubtractProductInPlace(gw_n, Z[n] - Z[m]).
                            MultiplyInPlace(Tr).
                            MultiplyInPlace(krw_n).
                            DivideInPlace(bw_n).
                            DivideInPlace(μw_n);
                    }
                    Res[2*m+1].AddProductInPlace(scratchpad.Water_Rate, dt);
                }

                Parallel.For(0, Ngrids,
                    // Thread-local initializer: This creates ONE scratchpad instance 
                    // per CPU thread, rather than one per loop iteration!
                    () => new ScratchPad(),

                    // The actual loop work
                    (m, loopState, scratchpad) =>
                    {
                        int indx1 = 2*m, indx2 = 2*m + 1;
                        scratchpad.ClearMVariables();
                        ComputeBlockPVT(scratchpad, m);

                        double PV = -Dx[m]*Dy[m]*Dz[m]*Φ[m]/beta;

                        ADiff SoPo = scratchpad.SoPo_m, SwPw = scratchpad.SwPw_m;

                        SoPo.AddInPlace(scratchpad.So_m).
                        MultiplyInPlace(scratchpad.Po_m);

                        SwPw.AddInPlace(scratchpad.Sw_m).
                        MultiplyInPlace(scratchpad.Pw_m);

                        scratchpad.meanP.CopyFrom(SoPo).AddInPlace(SwPw);
                        Er(scratchpad.meanP, scratchpad.Er);

                        ErSo_Bo_np1[m].
                        CopyFrom(scratchpad.Er).
                        MultiplyInPlace(scratchpad.So_m).
                        DivideInPlace(scratchpad.Bo_m);

                        ErSw_Bw_np1[m].
                        CopyFrom(scratchpad.Er).
                        MultiplyInPlace(scratchpad.Sw_m).
                        DivideInPlace(scratchpad.Bw_m);

                        Res[indx1].CopyFrom(ErSo_Bo_np1[m]).
                        SubtractInPlace(ErSo_Bo_n[m]).
                        MultiplyInPlace(PV);

                        Res[indx2].CopyFrom(ErSw_Bw_np1[m]).
                        SubtractInPlace(ErSw_Bw_n[m]).
                        MultiplyInPlace(PV);

                        int k = m / NxNy, rem = m % NxNy,
                        j = rem / Nx, i = rem % Nx;

                        // Add fluxes to the residual for this grid block from each of its 6 neighbors
                        if (i > 0) AddFluxes(xDir, m, m-1, scratchpad);
                        if (i < Lx) AddFluxes(xDir, m, m+1, scratchpad);
                        if (j > 0) AddFluxes(yDir, m, m-Nx, scratchpad);
                        if (j < Ly) AddFluxes(yDir, m, m+Nx, scratchpad);
                        if (k > 0) AddFluxes(zDir, m, m-NxNy, scratchpad);
                        if (k < Lz) AddFluxes(zDir, m, m+NxNy, scratchpad);

                        // Pass the scratchpad forward to be reused by the next iteration on this thread
                        return scratchpad;
                    },

                    // Thread-local cleanup (nothing to do here, but required by syntax)
                    doubleResult => { }
                );


                for (int nwell = 0; nwell < Nwells; nwell++)
                {
                    var well = Wells[nwell];
                    pwell.CopyFrom(xnp1[2*Ngrids + 2*nwell]);
                    qwell.CopyFrom(xnp1[2*Ngrids + 2*nwell + 1]);
                    Res[2*Ngrids + 2*nwell].Clear().AddProductInPlace(qwell, dt);
                    well.Constraint(time, pwell, qwell, Res[2*Ngrids + 2*nwell + 1]);
                    well.WaterRate = 0; well.OilRate = 0; Zref = well.Zref;

                    switch (well.WellType)
                    {
                        case WellType.Producer:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                int indx1 = 2*m, indx2 = 2*m+1;
                                Wellscratchpad.ClearMVariables();
                                ComputeBlockPVT(Wellscratchpad, m);

                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);

                                Wellscratchpad.Oil_Rate.CopyFrom(pwell).
                                    SubtractInPlace(Wellscratchpad.Po_m).
                                    SubtractProductInPlace(Wellscratchpad.Go_m, Zref - Z[m]).
                                    MultiplyInPlace(Wellscratchpad.Kro_m).
                                    DivideInPlace(Wellscratchpad.Uo_m).
                                    DivideInPlace(Wellscratchpad.Bo_m).
                                    MultiplyInPlace(WI);

                                well.OilRate += Wellscratchpad.Oil_Rate.Value;
                                Res[indx1].AddProductInPlace(Wellscratchpad.Oil_Rate, dt);
                                Res[2*Ngrids + 2*nwell].SubtractProductInPlace(Wellscratchpad.Oil_Rate, dt);

                                Wellscratchpad.Water_Rate.CopyFrom(pwell).
                                    SubtractInPlace(Wellscratchpad.Pw_m).
                                    SubtractProductInPlace(Wellscratchpad.Gw_m, Zref - Z[m]).
                                    MultiplyInPlace(Wellscratchpad.Krw_m).
                                    DivideInPlace(Wellscratchpad.Uw_m).
                                    DivideInPlace(Wellscratchpad.Bw_m).
                                    MultiplyInPlace(WI);

                                well.WaterRate += Wellscratchpad.Water_Rate.Value;
                                Res[indx2].AddProductInPlace(Wellscratchpad.Water_Rate, dt);
                                Res[2*Ngrids + 2*nwell].SubtractProductInPlace(Wellscratchpad.Water_Rate, dt);
                            }
                            break;

                        case WellType.Injector:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                int indx1 = 2*m, indx2 = 2*m+1;
                                Wellscratchpad.ClearMVariables();
                                ComputeBlockPVT(Wellscratchpad, m);

                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);

                                Wellscratchpad.Water_Rate.CopyFrom(pwell).
                                    SubtractInPlace(Wellscratchpad.Pw_m).
                                    SubtractProductInPlace(γw(pwell, tempvar), Zref - Z[m]).
                                    MultiplyInPlace(krw0).
                                    DivideInPlace(μw(pwell, tempvar)).
                                    DivideInPlace(Bw(pwell, tempvar)).
                                    MultiplyInPlace(WI);

                                well.WaterRate += Wellscratchpad.Water_Rate.Value;
                                Res[indx2].AddProductInPlace(Wellscratchpad.Water_Rate, dt);
                                Res[2*Ngrids + 2*nwell].SubtractProductInPlace(Wellscratchpad.Water_Rate, dt);
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
            }

            dt = 0.001;
            // Initialize historical data tracking containers for plotting and reporting
            P = [Po_n]; S = [Sw_n]; Rate = [Qwells_n];
            Pwf = [Pwells_n]; WaterCut = [new double[Nwells]];
            Time = [0.0]; SweepEff = [0.0];
            List<double> Interval = [0, .. ResultTime];
            foreach (var well in Wells)
                Interval.AddRange(well.ProductionProfile.Time);
            Interval = [.. Interval.Distinct().OrderBy(x => x)];
            Console.WriteLine($"""
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

                    // Call the Newton-Raphson nonlinear solver to
                    // find the next state solution
                    Console.WriteLine($"""
                    Time: 
                    {tnp1:F3} days
                        iter  |   Residual Norm  
                    ----------+----------------
                    """);
                    int iter; bool isConverged = false;
                    for (iter = 1; iter < 10; iter++)
                    {
                        // Solve the nonlinear system using Newton-Raphson method
                        Residual(xs, tnp1);  rnorm = b.Max(Abs);
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
                        GuessReset();
                        continue;
                    }

                    // Unpack solution values to evaluate operational constraint validations
                    double maxDP = 0, maxDS = 0;
                    for (int m = 0; m < Ngrids; m++)
                    {
                        maxDP = Max(maxDP, Abs(xs[2*m].Value - Po_n[m]));
                        maxDS = Max(maxDS, Abs(xs[2*m+1].Value - Sw_n[m]));
                    }
                    if(maxDP > 100 || maxDS > 0.15)
                    {
                        dt = 0.5*dt;
                        GuessReset();
                        continue;
                    }

                    ExtractSolution();
                    // Quality Check:
                    
                    // Track the maximum time step that still yields convergence for potential future use in adaptive time-stepping
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
                        GuessReset();
                        continue;
                    }

                    // Log verified parameters to performance history arrays

                    Time.Add(t = tnp1); isComplete = Abs(t - L) < 1e-13;
                    P.Add([..Po_n]); S.Add([..Sw_n]); Rate.Add([..Qwells_n]); Pwf.Add([..Pwells_n]);
                    WaterCut.Add([.. Wells.Select(w => w.WaterCut)]);
                    SweepEff.Add((Sw_n.Sum() - S[0].Sum())*100/(Ngrids - S[0].Sum()));
                    PoreVolumeUpdate();
                    // Adaptive Time-Stepping Logic:
                    // scale dt up if convergence is fast, scale down if slow
                    if (iter < 4) dt = 1.25*dt;
                    if (iter > 8) dt = 0.5*dt;
                    if (dt < 1e-5) throw new Exception("time step is too small");

                    // Prevent Overshoot
                    if (!isComplete) dt = Min(dt, L - t);
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
