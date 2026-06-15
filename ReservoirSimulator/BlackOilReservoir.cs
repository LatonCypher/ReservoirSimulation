using SepalSolver;
using static ReservoirSimulator.Math;

namespace ReservoirSimulator
{
    public class BlackOilReservoir
    {
        readonly Func<ADiff, ADiff> Sws, Swe, Pc_D, Pc_I, Bo, Bw, μo, μw, γo, γw, Kro, Krw, Er;
        // Define conversion constants
        const double alpha = 1.127e-3,        // Darcy to Field units factor
          alpha_well = 1.127e-3*2*pi,         // Darcy to Field units factor for wells
                beta = 5.615;                 // ft3 to bbl conversion factor
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

        readonly int Nx, Ny, Nz, NxNy, Ngrids, Nwells, varNum;
        readonly double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z;
        double[] Po_n, Sw_n, Pw_n, So_n, Qwells_n, Pwells_n, ErSw_Bw_n, ErSo_Bo_n;
        readonly double krw0, kro0, Pb, Pref, Pe, Pw_woc, Sw_r, So_r,
            Bo0, Bw0, Po_woc, Z_woc, co, cw, cr, bo, bw, nw, no,
            μo0, μw0, γo0, γw0;
        List<Well> Wells;
        public BlackOilReservoir(int _nx, int _ny, int _nz,
            double[] _perm, double[] _phi, double[] _dx, double[] _dy,
            double[] _dz, double[] _z, double _peow, double _pw_woc,
            double _z_woc, double _mult_z, double _sw_r, double _so_r,
            double _bo0, double _bw0, double _μo0, double _μw0, double _γo0,
            double _γw0, double _krw0, double _kro0, double _co, double _cw,
            double _cr, double _bo, double _bw, double _nw, double _no,
            double _pb, double _pref, List<Well> _wells)
        {
            ADiff.capacity = 23;
            Kx = _perm; Ky = _perm; Kz = [.. _perm.Select(k => k*_mult_z)];
            Nx = _nx; Ny = _ny; Nz = _nz; NxNy = Nx*Ny; Ngrids = Nx*Ny*Nz;
            Dx = _dx; Dy = _dy; Dz = _dz; Z = _z; Φ = _phi; kro0 = _kro0;
            krw0 = _krw0; Bo0 = _bo0; Bw0 = _bw0; Pb = _pb; Pref = _pref;
            Pe = _peow; So_r = _so_r; Sw_r = _sw_r; co = _co; cw = _cw;
            cr = _cr; bo = _bo; bw = _bw; no = _no; nw = _nw; μo0 = _μo0;
            μw0 = _μw0; γo0 = _γo0; γw0 = _γw0; Pw_woc = _pw_woc;
            Po_woc = Pw_woc + Pe; Z_woc = _z_woc;
            Wells = _wells; Nwells = Wells.Count; varNum = 2*Ngrids + 2*Nwells;

            Sws = Sw => (Sw - Sw_r)/(1 - Sw_r);
            Swe = Sw => (Sw - Sw_r)/(1 - Sw_r - So_r);
            Pc_D = Sw => Pe * Pow(Sws(Sw), -1.5);
            Pc_I = Sw => Pe * (Pow(Swe(Sw), -1.5) - 1);
            Bo = Po => Bo0*Exp(co*(Pb - Po));
            Bw = Pw => Bw0*Exp(cw*(Pref - Pw));
            μo = Po => μo0*Exp(bo*(Po - Pb));
            μw = Pw => μw0*Exp(bw*(Pw - Pref));
            γo = Po => γo0*Exp(bo*(Po - Pb));
            γw = Pw => γw0*Exp(bw*(Pw - Pref));
            Er = P => Exp(cr*(P - Pref));
            Kro = So => kro0 * Pow(1 - Swe(1 - So), no);
            Krw = Sw => krw0 * Pow(Swe(Sw), nw);

            // 1. Pre-generate a fine lookup table for the inverse Capillary Pressure relationship
            List<double> Sw_Table = [.. Linspace(1.0-So_r-1e-5, Sw_r+1e-5, 50)];
            // Calculate Pc for each Sw point in our table
            List<double> Pc_Table = [.. Sw_Table.Select(sw => Pc_D(sw).Value)];

            // 2. Initialize the spatial grid blocks
            Pw_n = new double[Ngrids]; Sw_n = new double[Ngrids];
            Po_n = new double[Ngrids]; So_n = new double[Ngrids];

            for (int i = 0; i < Ngrids; i++)
            {
                Pw_n[i] = Pw_woc + γw(Pw_woc).Value * (Z[i] - Z_woc);
                Po_n[i] = Po_woc + γo(Po_woc).Value * (Z[i] - Z_woc);
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

            Qwells_n = new double[Nwells];
            Pwells_n = new double[Nwells];
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
            (ErSo_Bo_n, ErSw_Bw_n) = FluidInGrid(Po_n, Pw_n, So_n, Sw_n);
        }

        (double[], double[]) FluidInGrid(double[] Po, double[] Pw, double[] So, double[] Sw)
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

        public void Simulate2Phase(double[] ResultTime, List<Well> Wells)
        {
            int Lx = Nx - 1, Ly = Ny - 1, Lz = Nz - 1, indx1, indx2;
            double dt, t = 0;

            Phase wPhase = Phase.Water, oPhase = Phase.Oil;
            Direction xDir = Direction.X, yDir = Direction.Y, zDir = Direction.Z;

            List<double> a_value = [], b = [];
            List<int> a_index = [], a_start = [0];
            ADiff[] Res = new ADiff[varNum], xs;
            void Residual(ADiff[] xnp1, double time)
            {
                double re, WI, V, Zref; ADiff WIw, WIo, pwell, qwell;

                ADiff Flux(Phase phase, Direction dir, int m, int n)
                {
                    ADiff Po_up, Pw_up, So_up, Sw_up, Tw, To;
                    double Tr = Transmissibility(dir, m, n);
                    ADiff po_m = xnp1[2*m], sw_m = xnp1[2*m+1],
                        pw_m = po_m - Pc_I(sw_m),
                        so_m = 1 - sw_m,
                        po_n = xnp1[2*n], sw_n = xnp1[2*n+1],
                        pw_n = po_n - Pc_I(sw_n),
                        so_n = 1 - sw_n;
                    switch (phase)
                    {
                        case Phase.Water:
                            (Pw_up, Sw_up) = pw_m - γw(pw_m)*Z[m] > pw_n - γw(pw_n)*Z[n] ?
                                (pw_m, sw_m) : (pw_n, sw_n);
                            Tw = Tr*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                            return Tw*(pw_n - pw_m - γw(Pw_up)*(Z[n] - Z[m]));
                        case Phase.Oil:
                            (Po_up, So_up) = po_m - γo(po_m)*Z[m] > po_n - γo(po_n)*Z[n] ?
                                (po_m, so_m) : (po_n, so_n);
                            To = Tr*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                            return To*(po_n - po_m - γo(Po_up)*(Z[n] - Z[m]));
                        default:
                            throw new ArgumentException("Invalid phase");
                    }
                }

                //for (int m = 0; m < Ngrids; m++)
                Parallel.For(0, Ngrids, m =>
                {
                    int indx1 = 2*m, indx2 = 2*m + 1;
                    V = Dx[m]*Dy[m]*Dz[m]/beta;
                    ADiff po_m = xnp1[2*m], sw_m = xnp1[2*m+1],
                        pw_m = po_m - Pc_I(sw_m),
                        so_m = 1 - sw_m;
                    var meanP = so_m*po_m + sw_m*pw_m;
                    var erso_bo_np1 = Er(meanP)*so_m/Bo(po_m);
                    var ersw_bw_np1 = Er(meanP)*sw_m/Bw(pw_m);
                    Res[indx1] = -V*Φ[m]*(erso_bo_np1 - ErSo_Bo_n[m])/dt;
                    Res[indx2] = -V*Φ[m]*(ersw_bw_np1 - ErSw_Bw_n[m])/dt;

                    int k = m / NxNy, rem = m % NxNy,
                    j = rem / Nx, i = rem % Nx;

                    if (i > 0)
                    {
                        Res[indx1] += Flux(oPhase, xDir, m, m-1);
                        Res[indx2] += Flux(wPhase, xDir, m, m-1);
                    }
                    if (i < Lx)
                    {
                        Res[indx1] += Flux(oPhase, xDir, m, m+1);
                        Res[indx2] += Flux(wPhase, xDir, m, m+1);
                    }

                    if (j > 0)
                    {
                        Res[indx1] += Flux(oPhase, yDir, m, m-Nx);
                        Res[indx2] += Flux(wPhase, yDir, m, m-Nx);
                    }
                    if (j < Ly)
                    {
                        Res[indx1] += Flux(oPhase, yDir, m, m+Nx);
                        Res[indx2] += Flux(wPhase, yDir, m, m+Nx);
                    }

                    if (k > 0)
                    {
                        Res[indx1] += Flux(oPhase, zDir, m, m-NxNy);
                        Res[indx2] += Flux(wPhase, zDir, m, m-NxNy);
                    }
                    if (k < Lz)
                    {
                        Res[indx1] += Flux(oPhase, zDir, m, m+NxNy);
                        Res[indx2] += Flux(wPhase, zDir, m, m+NxNy);
                    }
                });

                for (int nwell = 0; nwell < Nwells; nwell++)
                {
                    var well = Wells[nwell];
                    pwell = xnp1[2*Ngrids + 2*nwell];
                    qwell = xnp1[2*Ngrids + 2*nwell + 1];
                    Res[2*Ngrids + 2*nwell] = qwell;
                    Res[2*Ngrids + 2*nwell + 1] = well.Constraint(time, pwell, qwell);
                    well.WaterRate = 0; well.OilRate = 0; Zref = well.Zref;
                    ADiff water_rate = 0, oil_rate = 0;
                    switch (well.WellType)
                    {
                        case WellType.Producer:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                indx1 = 2*m; indx2 = 2*m+1;
                                ADiff po_m = xnp1[indx1], sw_m = xnp1[indx2],
                                pw_m = po_m - Pc_I(sw_m), so_m = 1 - sw_m;
                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);
                                WIw = WI*Krw(sw_m)/(μw(pw_m)*Bw(pw_m));
                                WIo = WI*Kro(so_m)/(μo(po_m)*Bo(po_m));
                                water_rate = (pwell - pw_m - γw(pw_m)*(Zref - Z[m]))*WIw;
                                oil_rate = (pwell - po_m - γo(po_m)*(Zref - Z[m]))*WIo;
                                well.WaterRate += water_rate.Value;
                                well.OilRate += oil_rate.Value;
                                Res[indx1] += oil_rate; Res[indx2] += water_rate;
                                Res[2*Ngrids + 2*nwell] -= oil_rate + water_rate;
                            }
                            break;

                        case WellType.Injector:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                indx1 = 2*m; indx2 = 2*m+1;
                                ADiff po_m = xnp1[indx1], sw_m = xnp1[indx2],
                                pw_m = po_m - Pc_I(sw_m), so_m = 1 - sw_m;
                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);
                                WIw = WI*krw0/(μw(pwell)*Bw(pwell));
                                water_rate = (pwell - pw_m - γw(pwell)*(Zref - Z[m]))*WIw;
                                well.WaterRate += water_rate.Value;
                                Res[indx2] += water_rate;
                                Res[2*Ngrids + 2*nwell] -= water_rate;
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

            dt = 0.01;
            xs = [.. Enumerable.Range(0, varNum).Select(i => new ADiff(0, i))];
            Pack(Po_n, Sw_n, Pwells_n, Qwells_n, xs);
            // Initialize historical data tracking containers for plotting and reporting
            List<double[]> P = [Po_n], S = [Sw_n], Rate = [Qwells_n],
                Pwf = [Pwells_n], WaterCut = [new double[Nwells]];
            List<double> Time = [0.0], SweepEff = [0.0], Interval = [0, .. ResultTime];
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
                        Residual(xs, tnp1);
                        double[] dx = MklSparseSolver.Solve([.. a_value], [.. a_index], [.. a_start], [.. b]);
                        for (int v = 0; v < varNum; v++) xs[v].Value -= dx[v];
                        rnorm = Sqrt(b.Select(x => x*x).Sum());
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

                            //// switch to rate control if production rate is above target
                            //if (Wells[n].ConstraintType == ConstraintType.MinPressure
                            //    && Qwells_n[n] < Wells[n].ProdRate(t + dt))
                            //{
                            //    Wells[n].ConstraintType = ConstraintType.FlowRate;
                            //    Console.WriteLine("""
                            //    ================================================
                            //          Rejected (Target Production Rate Exceeded) 
                            //    ================================================
                            //    """);
                            //    staterejected = true;
                            //    break;
                            //}
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

                            //// switch to rate control if injection rate exceeds rate target
                            //if (Wells[n].ConstraintType == ConstraintType.MaxPressure &&
                            //    Qwells_n[n] > Wells[n].ProdRate(t + dt))
                            //{
                            //    Wells[n].ConstraintType = ConstraintType.FlowRate;
                            //    Console.WriteLine("""
                            //    ================================================
                            //          Rejected (Target Injection Rate Exceeded) 
                            //    ================================================
                            //    """);
                            //    staterejected = true;
                            //    break;
                            //}
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
                    (ErSo_Bo_n, ErSw_Bw_n) = FluidInGrid(Po_n, Pw_n, So_n, Sw_n);
                    // Adaptive Time-Stepping Logic:
                    // scale dt up if convergence is fast, scale down if slow
                    if (iter < 4) dt = 1.25*dt;
                    if (iter > 8) dt = 0.5*dt;
                    if (dt < 1e-5) throw new Exception("time step is too small");

                    // Prevent Overshoot
                    if (!isComplete) dt = Min(dt, L - t);
                }

                if (history_dtmax > 0)
                    dt = history_dtmax;
            }
        }
    }
}
