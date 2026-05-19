using static System.Runtime.InteropServices.JavaScript.JSType;

namespace ResSim
{
    public class Reservoir2PhasesAD
    {

        Func<AutoDiff, AutoDiff> Sws, Swe, Pc_D, Pc_I, Bo, Bw, μo, μw, γo, γw, Er, Krw, Kro;
        // Define conversion constants
        double alpha = 1.127e-3,        // Darcy to Field units factor
          alpha_well = 1.127e-3*2*pi,   // Darcy to Field units factor for wells
                beta = 5.615;           // ft3 to bbl conversion factor
        double Transmissibility(Direction d, int m, int n)
        {
            return d switch
            {
                Direction.X => alpha * Harmmean(Dy[m] * Dz[m] * Kx[m] / Dx[m], Dy[n] * Dz[n] * Kx[n] / Dx[n]),
                Direction.Y => alpha * Harmmean(Dx[m] * Dz[m] * Ky[m] / Dy[m], Dx[n] * Dz[n] * Ky[n] / Dy[n]),
                Direction.Z => alpha * Harmmean(Dx[m] * Dy[m] * Kz[m] / Dz[m], Dx[n] * Dy[n] * Kz[n] / Dz[n]),
                _ => throw new ArgumentException("Invalid direction"),
            };
        }
        double betweenab(double a, double b, double f) => a + f*(b-a);
        double interps(List<double> X, List<double> Y, double x)
        {
            int i = X.FindIndex(xi => xi>x);
            double f = (x-X[i-1])/(X[i]-X[i-1]);
            return betweenab(Y[i-1], Y[i], f);
        }
        double[] interpa(List<double> X, List<double[]> Y, double x)
        {
            int i = X.FindIndex(xi => xi>x);
            double f = (x-X[i-1])/(X[i]-X[i-1]);
            return [.. Y[i-1].Zip(Y[i], (a, b) => betweenab(a, b, f))];
        }
        void WriteArray(double[] V)
        {
            for (int i = 0; i < Ny; i++)
                Console.WriteLine(string.Join(", ", V[(i*Nx)..((i+1)*Nx)].Select(x => x.ToString("F3"))));
        }
        double Harmmean(double x1, double x2) => 2/(1/x1 + 1/x2);
        (AutoDiff[], AutoDiff[], AutoDiff[], AutoDiff[]) Unpack(AutoDiff[] x)
        {
            int indx = 0;
            AutoDiff[] Po = new AutoDiff[Ngrids], Sw = new AutoDiff[Ngrids],
                     Pwells = new AutoDiff[Nwells], Qwells = new AutoDiff[Nwells];
            for (int i = 0; i < Ngrids; i++)
            {
                Po[i] = x[indx++]; // Matches Po index
                Sw[i] = x[indx++]; // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                Pwells[i] = x[indx++]; // Matches Pwf index
                Qwells[i] = x[indx++]; // Matches Q index
            }
            return (Po, Sw, Pwells, Qwells);
        }

        (double[], double[], double[], double[]) AutoDiff2Double(AutoDiff[] x)
        {
            int indx = 0;
            double[] Po = new double[Ngrids], Sw = new double[Ngrids],
                     Pwells = new double[Nwells], Qwells = new double[Nwells];
            for (int i = 0; i < Ngrids; i++)
            {
                Po[i] = x[indx++].Value; // Matches Po index
                Sw[i] = x[indx++].Value; // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                Pwells[i] = x[indx++].Value; // Matches Pwf index
                Qwells[i] = x[indx++].Value; // Matches Q index
            }
            return (Po, Sw, Pwells, Qwells);
        }
        AutoDiff[] Pack(AutoDiff[] Ro, AutoDiff[] Rw, AutoDiff[] Rwells, AutoDiff[] Rdecision)
        {
            int indx = 0;
            AutoDiff[] R_total = new AutoDiff[2*Ngrids + 2*Nwells];
            for (int i = 0; i < Ngrids; i++)
            {
                R_total[indx++] = Ro[i]; // Matches Po index
                R_total[indx++] = Rw[i]; // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                R_total[indx++] = Rwells[i];    // Matches Pwf index
                R_total[indx++] = Rdecision[i]; // Matches Q index
            }
            return R_total;
        }

        AutoDiff[] Double2AutoDiff(double[] Ro, double[] Rw, double[] Rwells, double[] Rcontrol)
        {
            int indx = 0, varNum = 2*Ngrids + 2*Nwells;
            AutoDiff[] R_total = new AutoDiff[varNum];
            for (int i = 0; i < Ngrids; i++)
            {
                R_total[indx] = new(Ro[i], varNum, indx); indx++; // Matches Po index
                R_total[indx] = new(Rw[i], varNum, indx); indx++; // Matches Sw index
            }
            for (int i = 0; i < Nwells; i++)
            {
                R_total[indx] = new(Rwells[i], varNum, indx); indx++;   // Matches Pwf index
                R_total[indx] = new(Rcontrol[i], varNum, indx); indx++; // Matches Q index
            }
            return R_total;
        }

        int Nx, Ny, Nz, NxNy, Ngrids, Nwells;
        double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z;
        double[] Po_n, Sw_n, Pw_n, So_n, Qwells_n, Pwells_n, ErSw_Bw_n, ErSo_Bo_n;
        double krw0, kro0, Pb, Pref, Pe, Pw_woc, Sw_r, So_r, Bo0, Bw0, Po_woc, Z_woc;
        public Reservoir2PhasesAD(int _nx, int _ny, int _nz, double[] _perm, double[] _phi,
            double[] dx, double[] dy, double[] dz, double[] z, double peow, double pw_woc,
            double z_woc, double mult_y, double sw_r, double so_r, double bo0, double bw0,
            double μo0, double μw0, double γo0, double γw0, double _krw0, double _kro0,
            double co, double cw, double cr, double bo, double bw, double nw, double no,
            double pb, double pref)
        {
            Kx = _perm; Ky = _perm; Kz = [.. _perm.Select(k => k*mult_y)];
            Nx = _nx; Ny = _ny; Nz = _nz; NxNy = Nx*Ny; Ngrids = Nx*Ny*Nz;
            Dx = dx; Dy = dy; Dz = dz; Z = z; Φ = _phi;
            kro0 = _kro0; krw0 = _krw0;
            Pb = pb; Pref = pref; Pe = peow; 
            Pw_woc = pw_woc; Po_woc = Pw_woc + Pe; Z_woc = z_woc;
            Sws = Sw => (Sw - Sw_r)/(1 - Sw_r);
            Swe = Sw => (Sw - Sw_r)/(1 - Sw_r - So_r);
            Pc_D = Sw => Pe * Pow(Sws(Sw), -0.5);
            Pc_I = Sw => Pe * (Pow(Swe(Sw), -0.5) - 1);
            Bo = Po => Bo0*Exp(co*(Pb - Po));
            Bw = Pw => Bw0*Exp(cw*(Pref - Pw));
            μo = Po => μo0*Exp(bo*(Po - Pb));
            μw = Pw => μw0*Exp(bw*(Pw - Pref));
            γo = Po => γo0*Exp(bo*(Po - Pb));
            γw = Pw => γw0*Exp(bw*(Pw - Pref));
            Er = P => Exp(cr*(P-1500));
            Kro = So => kro0 * Pow(1 - Swe(1 - So), no);
            Krw = Sw => krw0 * Pow(Swe(Sw), nw);
        }

        public void Initialize(List<Well> Wells)
        {
            // Hydrostatic pressure calculation (assuming Z increases downward)
            // P = P_woc + rho * g * deltaZ
            double CalculateWaterPressure(int i) => (Pw_woc + γw(Pw_woc) * (Z[i] - Z_woc)).Value;

            // Oil pressure accounts for Capillary Entry Pressure (Pe) at the Water-Oil Contact (WOC)
            double CalculateOilPressure(int i) => (Po_woc + γo(Po_woc) * (Z[i] - Z_woc)).Value;

            // 1. Pre-generate a fine lookup table for the inverse Capillary Pressure relationship
            int tablePoints = 50;
            List<double> Sw_Table = [.. Linspace(Sw_r + 1e-5, 1.0, tablePoints)];

            // Calculate Pc for each Sw point in our table
            List<double> Pc_Table = [.. Sw_Table.Select(sw => Pc_D(sw).Value)];

            // CRITICAL STEP: Because interps expects the independent variable array (X) to be 
            // strictly increasing, and Pc drops monotonically as Sw increases, we must reverse 
            // both collections so that Pc moves from lowest to highest.
            Pc_Table.Reverse();
            Sw_Table.Reverse();

            // 2. Initialize the spatial grid blocks
            for (int i = 0; i < Ngrids; i++)
            {
                Pw_n[i] = CalculateWaterPressure(i);
                Po_n[i] = CalculateOilPressure(i);
                double pc = Po_n[i] - Pw_n[i];

                // 3. Directly interpolate saturation instead of using an iterative solver
                if (pc > Pe && pc <= Pc_Table.Last())
                    // interps(List<double> X, List<double> Y, double x)
                    Sw_n[i] = interps(Pc_Table, Sw_Table, pc);
                else if (pc > Pc_Table.Last())
                    // Capillary pressure exceeds our table limit; clamp to residual water
                    Sw_n[i] = Sw_r;
                else
                    // Below or at the entry boundary threshold
                    Sw_n[i] = 1.0;

                So_n[i] = 1.0 - Sw_n[i];
            }

            // Corrected accumulation updates mapping block-specific localized pressures
            ErSw_Bw_n = [.. Sw_n.Zip(Pw_n, (sw, pw) => (Er(pw) * sw / Bw(pw)).Value)];
            ErSo_Bo_n = [.. So_n.Zip(Po_n, (so, po) => (Er(po) * so / Bo(po)).Value)];

            for (int i = 0; i < Nwells; i++)
            {
                Qwells_n[i] = 0; // Start with no flow
                switch (Wells[i].WellType)
                {
                    case WellType.Producer:
                        Pwells_n[i] = Po_n[Wells[i].Perforation_NatIndex.First()];
                        break;
                    case WellType.Injector:
                        Pwells_n[i] = Pw_n[Wells[i].Perforation_NatIndex.Last()];
                        break;
                }
            }
        }

        public void Simulate2Phase(double[] ResultTime, List<Well> Wells)
        {
            int Lx = Nx - 1, Ly = Ny - 1, Lz = Nz - 1; Nwells = Wells.Count;
            double dt, t = 0;

            Phase wPhase = Phase.Water, oPhase = Phase.Oil;
            Direction xDir = Direction.X, yDir = Direction.Y, zDir = Direction.Z;

            (ColVec Res, Matrix Jac) Residual(AutoDiff[] xnp1, double time)
            {
                double re, WI; AutoDiff WIw, WIo;
                var (Po_np1, Sw_np1, Pwells_np1, Qwells_np1) = Unpack(xnp1);
                AutoDiff[] Pw_np1 = [.. Po_np1.Zip(Sw_np1, (po, sw) => po - Pc_I(sw))],
                    So_np1 = [.. Sw_np1.Select(sw => 1 - sw)];

                bool ComparePotential(Phase phase, int m, int n)
                {
                    return phase switch
                    {
                        Phase.Water => Pw_np1[m] - γw(Pw_np1[m])*Z[m] > Pw_np1[n] - γw(Pw_np1[n])*Z[n],
                        Phase.Oil => Po_np1[m] - γo(Po_np1[m])*Z[m] > Po_np1[n] - γo(Po_np1[n])*Z[n],
                        _ => throw new ArgumentException("Invalid phase"),
                    };
                }

                AutoDiff Flux(Phase phase, Direction dir, int m, int n)
                {
                    AutoDiff Po_up, Pw_up, So_up, Sw_up, Tw, To;
                    double Tr = Transmissibility(dir, m, n);
                    switch (phase)
                    {
                        case Phase.Water:
                            (Pw_up, Sw_up) = ComparePotential(phase, m, n) ?
                                (Pw_np1[m], Sw_np1[m]) : (Pw_np1[n], Sw_np1[n]);
                            Tw = Tr*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                            return Tw*(Pw_np1[m] - Pw_np1[n] - γw(Pw_up)*(Z[m] - Z[n]));
                        case Phase.Oil:
                            (Po_up, So_up) = ComparePotential(phase, m, n) ?
                                (Po_np1[m], So_np1[m]) : (Po_np1[n], So_np1[n]);
                            To = Tr*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                            return To*(Po_np1[m] - Po_np1[n] - γo(Po_up)*(Z[m] - Z[n]));
                        default:
                            throw new ArgumentException("Invalid phase");
                    }
                }

                AutoDiff[] Rw = new AutoDiff[Ngrids], Ro = new AutoDiff[Ngrids],
                    Rwells = new AutoDiff[Nwells], Rcontrol = new AutoDiff[Nwells];

                for (int m = 0; m < Ngrids; m++)
                {
                    double V = Dx[m]*Dy[m]*Dz[m]/beta;
                    Rw[m] = -V*Φ[m]*(Er(Pw_np1[m])*Sw_np1[m]/Bw(Pw_np1[m]) - ErSw_Bw_n[m])/dt;
                    Ro[m] = -V*Φ[m]*(Er(Po_np1[m])*So_np1[m]/Bo(Po_np1[m]) - ErSo_Bo_n[m])/dt;

                    var (k, rem) = DivRem(m, NxNy);
                    var (j, i) = DivRem(rem, Nx);

                    if (i > 0)
                    {
                        Rw[m] += Flux(wPhase, xDir, m-1, m);
                        Ro[m] += Flux(oPhase, xDir, m-1, m);
                    }
                    if (i < Lx)
                    {
                        Rw[m] += Flux(wPhase, xDir, m+1, m);
                        Ro[m] += Flux(oPhase, xDir, m+1, m);
                    }

                    if (j > 0)
                    {
                        Rw[m] += Flux(wPhase, yDir, m-Nx, m);
                        Ro[m] += Flux(oPhase, yDir, m-Nx, m);
                    }
                    if (j < Ly)
                    {
                        Rw[m] += Flux(wPhase, yDir, m+Nx, m);
                        Ro[m] += Flux(oPhase, yDir, m+Nx, m);
                    }

                    if (k > 0)
                    {
                        Rw[m] += Flux(wPhase, zDir, m-NxNy, m);
                        Ro[m] += Flux(oPhase, zDir, m-NxNy, m);
                    }
                    if (k < Lz)
                    {
                        Rw[m] += Flux(wPhase, zDir, m+NxNy, m);
                        Ro[m] += Flux(oPhase, zDir, m+NxNy, m);
                    }
                }

                for (int count = 0; count < Nwells; count++)
                {
                    var well = Wells[count];
                    Rwells[count] += Qwells_np1[count];
                    well.WaterRate = 0; well.OilRate = 0;
                    AutoDiff water_rate, oil_rate; double Zref = well.Zref;
                    switch (well.WellType)
                    {
                        case WellType.Producer:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);
                                WIw = WI*Krw(Sw_np1[m])/(μw(Pw_np1[m])*Bw(Pw_np1[m]));
                                WIo = WI*Kro(So_np1[m])/(μo(Po_np1[m])*Bo(Po_np1[m]));
                                water_rate = (Pwells_np1[count] - Pw_np1[m] - γw(Pwells_np1[count])*(Zref - Z[m]))*WIw;
                                oil_rate = (Pwells_np1[count] - Po_np1[m] - γo(Pwells_np1[count])*(Zref - Z[m]))*WIo;
                                well.WaterRate += water_rate.Value; well.OilRate += oil_rate.Value;
                                Rw[m] += water_rate; Ro[m] += oil_rate;
                                Rwells[count] -= oil_rate + water_rate;
                            }
                            break;

                        case WellType.Injector:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);
                                WIw = WI*krw0/(μw(Pw_np1[m])*Bw(Pw_np1[m]));
                                water_rate = (Pwells_np1[count] - Pw_np1[m] - γw(Pwells_np1[count])*(Zref - Z[m]))*WIw;
                                well.WaterRate += water_rate.Value;
                                Rw[m] += water_rate;
                                Rwells[count] -= water_rate;
                            }
                            break;
                    }
                    Rcontrol[count] = well.Constraint(time, Pwells_np1[count], Qwells_np1[count]);
                }
                var rtotal = Pack(Ro, Rw, Rwells, Rcontrol);
                ColVec res = rtotal.Select(r => r.Value).ToArray();
                Matrix jac = rtotal.Select(r => (RowVec)r.Derivatives).ToArray();
                return (res, jac);
            }

            dt = 0.01;
            Initialize(Wells);
            AutoDiff[] xs = Double2AutoDiff(Po_n, Sw_n, Pwells_n, Qwells_n);
            // Initialize historical data tracking containers for plotting and reporting
            List<double[]> P = [Po_n], S = [Sw_n];
            List<double> Time = [0.0], WaterCut = [0.0], SweepEff = [0.0],
                ProdRate = [Qwells_n[0]], InjRate = [Qwells_n[1]],
                ProdPwf = [Pwells_n[0]], InjPwf = [Pwells_n[1]];
            Console.WriteLine($"""
                    ======================================================================
                                            Starting simulation

                    Time: 
                    {0:F2} days



                    """);
            while (t < Time.Last())
            {
                int iter = 0; bool isConverged = false;
                for (iter = 0; iter < 10; iter++)
                {
                    // Solve the nonlinear system using Newton-Raphson method
                    var (res, jac) = Residual(xs, t + dt);
                    var dx = LinSolve(jac, res);
                    foreach(var (xi, dxi) in xs.Zip(dx, (xi, dxi) => (xi, dxi)))
                        xi.Value -= dxi;
                    isConverged = res.Select(Abs).Max() < 1e-6;
                    if(isConverged) break;
                }

                // Check convergence. If non-converged,
                // chop the time step (time-step cuts) and retry.
                if (!isConverged)
                {
                    dt = 0.25*dt;
                    xs = Double2AutoDiff(Po_n, Sw_n, Pwells_n, Qwells_n);
                    Console.WriteLine("""
                                ================================================
                                           Rejected (Non-Convergence)
                                ================================================
                                """);
                    continue;
                }

                // Unpack solution values to evaluate operational constraint validations
                var (Po_s, Sw_s, Pwells_s, Qwells_s) = AutoDiff2Double(xs);

                // Validate Producer constraints:
                // switch to BHP control if pressure falls below minimum limits
                for (var i = 0; i < Nwells; i++)
                {
                    if (Pwells_s[i] < Wells[i].MinPressure)
                    {
                        Wells[i].ConstraintType = ConstraintType.MinPressure;
                        Console.WriteLine("""
                                ================================================
                                      Rejected (Minimum Pressure Violated) 
                                ================================================
                                """);
                        continue;
                    }
                    // Validate Injector constraints:
                    // switch to BHP control if pressure exceeds fracturing limits
                    if (Pwells_s[i] > Wells[i].MaxPressure)
                    {
                        Wells[i].ConstraintType = ConstraintType.MaxPressure;
                        Console.WriteLine("""
                                ================================================
                                      Rejected (Maximum Pressure Violated) 
                                ================================================
                                """);
                        continue;
                    }
                }

                // Accept the validated time step and update internal tracking states
                (Po_n, Sw_n, Pwells_n, Qwells_n) = (Po_s, Sw_s, Pwells_s, Qwells_s);
                Pw_n = [.. Po_n.Zip(Sw_n, (po, sw) => (po - Pc_I(sw)).Value)];
                So_n = [.. Sw_n.Select(sw => 1 - sw)];
                ErSw_Bw_n = [.. Sw_n.Zip(Pw_n, (sw, pw) => (Er(pw) * sw / Bw(pw)).Value)];
                ErSo_Bo_n = [.. So_n.Zip(Po_n, (so, po) => (Er(po) * so / Bo(po)).Value)];

                // Log verified parameters to performance history arrays
                P.Add(Po_n); S.Add(Sw_n);
                //ProdPwf.Add(Pwells_n[0]); InjPwf.Add(Pwells_n[1]);
                //ProdRate.Add(Qwells_n[0]); InjRate.Add(Qwells_n[1]);
                //WaterCut.Add(Producer.WaterRate*100/(Producer.WaterRate + Producer.OilRate));
                SweepEff.Add((Sw_n.Sum() - S[0].Sum())*100/(Ngrids - S[0].Sum()));
                Time.Add(Time.Last() + dt);

                // Adaptive Time-Stepping Logic:
                // scale dt up if convergence is fast, scale down if slow
                if (iter < 4) dt = 1.25*dt;
                if (iter > 8) dt = 0.5*dt;
                if (dt < 1e-5) throw new Exception("time step is too small");
            }
        }

    }
}
