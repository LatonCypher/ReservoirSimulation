using ScottPlot;
using System;
using System.Collections.Generic;
using System.Text;

namespace ResSim
{
    public class Reservoir2Phases
    {

        Func<double, double> Sws, Swe, Pc_D, Pc_I, Bo, Bw, μo, μw, γo, γw, Er, Krw, Kro;
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
        (double[], double[], double[], double[]) Unpack(double[] x)
        {
            int indx = 0;
            double[] Po = Zeros(M), Sw = Zeros(M),
                     Pwells = Zeros(Nwells), Qwells = Zeros(Nwells);
            for (int i = 0; i < M; i++)
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
        double[] Pack(double[] Ro, double[] Rw, double[] Rwells, double[] Rdecision)
        {
            int indx = 0;
            double[] R_total = Zeros(2*M + 2*Nwells);
            for (int i = 0; i < M; i++)
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

        int Nx, Ny, Nz, NxNy, M, Nwells;
        double[] Kx, Ky, Kz, Φ, Dx, Dy, Dz, Z;
        double[] Po_n, Sw_n, Pw_n, So_n, ErSw_Bw_n, ErSo_Bo_n;
        double krw0, kro0, Pb, Pref, Pe, Pw_woc, Sw_r, So_r, Bo0, Bw0, Po_woc, Z_woc;
        public Reservoir2Phases(int _nx, int _ny, int _nz, double[] _perm, double[] _phi,
            double[] dx, double[] dy, double[] dz, double[] z, double peow, double pw_woc,
            double z_woc, double mult_y, double sw_r, double so_r, double bo0, double bw0,
            double μo0, double μw0, double γo0, double γw0, double _krw0, double _kro0,
            double co, double cw, double cr, double bo, double bw, double nw, double no,
            double pb, double pref)
        {
            Kx = _perm; Ky = _perm; Kz = [.. _perm.Select(k => k*mult_y)];
            Nx = _nx; Ny = _ny; Nz = _nz; NxNy = Nx*Ny; M = Nx*Ny*Nz;
            Dx = dx; Dy = dy; Dz = dz; Z = z; Φ = _phi;
            kro0 = _kro0; krw0 = _krw0; Bo0 = bo0; Bw0 = bw0;
            Pb = pb; Pref = pref; Pe = peow; Sw_r = sw_r; So_r = so_r;
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

        public void Simulate2Phase(double[] Time, List<Well> Wells)
        {
            int Lx = Nx - 1, Ly = Ny - 1, Lz = Nz - 1;  Nwells = Wells.Count;
            double dt;
            double[] Po_n, Sw_n, Pw_n, So_n;
            bool[] RateControl = [true, true];

            Phase wPhase = Phase.Water, oPhase = Phase.Oil;
            Direction xDir = Direction.X, yDir = Direction.Y, zDir = Direction.Z;

            double[] Residual(double[] xnp1)
            {
                double re, WI, WIw, WIo;
                var (Po_np1, Sw_np1, Pwells_np1, Qwells_np1) = Unpack(xnp1);
                double[] Pw_np1 = [.. Po_np1.Zip(Sw_np1, (po, sw) => po - Pc_I(sw))],
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

                double Flux(Phase phase, Direction dir, int m, int n)
                {
                    double Po_up, Pw_up, So_up, Sw_up, Tr, Tw, To;
                    Tr = Transmissibility(dir, m, n);
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

                double[] Rw = Zeros(M), Ro = Zeros(M),
                    Rwells = Zeros(Nwells), Rcontrol = Zeros(Nwells);

                for (int m = 0; m < M; m++)
                {
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

                    double V = Dx[m]*Dy[m]*Dz[m]/beta;
                    Rw[m] -= V*Φ[m]*(Er(Pw_np1[m])*Sw_np1[m]/Bw(Pw_np1[m]) - Er(Pw_n[m])*Sw_n[m]/Bw(Pw_n[m]))/dt;
                    Ro[m] -= V*Φ[m]*(Er(Po_np1[m])*So_np1[m]/Bo(Po_np1[m]) - Er(Po_n[m])*So_n[m]/Bo(Po_n[m]))/dt;
                }

                for (int count = 0; count < Nwells; count++)
                {
                    var well = Wells[count];
                    Rwells[count] += Qwells_np1[count];
                    well.WaterRate = 0; well.OilRate = 0;
                    double water = 0, oil = 0, Zref = well.Zref;
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
                                water = (Pwells_np1[count] - Pw_np1[m] - γw(Pwells_np1[count])*(Zref - Z[m]))*WIw;
                                oil = (Pwells_np1[count] - Po_np1[m] - γo(Pwells_np1[count])*(Zref - Z[m]))*WIo;
                                well.WaterRate += water; well.OilRate += oil;
                                Rw[m] += water; Ro[m] += oil;
                                Rwells[count] -= oil + water;
                            }
                            break;

                        case WellType.Injector:
                            foreach (int m in well.Perforation_NatIndex)
                            {
                                re = 0.28*Hypot(Pow(Ky[m]/Kx[m], 0.25)*Dx[m], Pow(Kx[m]/Ky[m], 0.25)*Dy[m])/
                                    (Pow(Ky[m]/Kx[m], 0.25) + Pow(Kx[m]/Ky[m], 0.25));
                                WI = alpha_well*Sqrt(Kx[m]*Ky[m])*Dz[m]/(Log(re/well.Radius) + well.Skin);
                                WIw = WI*krw0/(μw(Pw_np1[m])*Bw(Pw_np1[m]));
                                water = (Pwells_np1[count] - Pw_np1[m] - γw(Pwells_np1[count])*(Zref - Z[m]))*WIw;
                                well.WaterRate += water;
                                Rw[m] += water;
                                Rwells[count] -= water;
                            }
                            break;
                    }
                    Rcontrol[count] = well.Constraint(Time.Last() + dt, Pwells_np1[count], Qwells_np1[count]);
                }
                return Pack(Ro, Rw, Rwells, Rcontrol);
            }
        }

        public void Initialize()
        {
            // Hydrostatic pressure calculation (assuming Z increases downward)
            // P = P_woc + rho * g * deltaZ
            double CalculateWaterPressure(int i) => Pw_woc + γw(Pw_woc) * (Z[i] - Z_woc);

            // Oil pressure accounts for Capillary Entry Pressure (Pe) at the Water-Oil Contact (WOC)
            double CalculateOilPressure(int i) => Po_woc + γo(Po_woc) * (Z[i] - Z_woc);

            // 1. Pre-generate a fine lookup table for the inverse Capillary Pressure relationship
            int tablePoints = 500;
            List<double> Sw_Table = [.. Linspace(Sw_r + 1e-5, 1.0, tablePoints)];

            // Calculate Pc for each Sw point in our table
            List<double> Pc_Table = [.. Sw_Table.Select(sw => Pc_D(sw))];

            // CRITICAL STEP: Because interps expects the independent variable array (X) to be 
            // strictly increasing, and Pc drops monotonically as Sw increases, we must reverse 
            // both collections so that Pc moves from lowest to highest.
            Pc_Table.Reverse();
            Sw_Table.Reverse();

            // 2. Initialize the spatial grid blocks
            for (int i = 0; i < M; i++)
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
            ErSw_Bw_n = [.. Sw_n.Zip(Pw_n, (sw, pw) => Er(pw) * sw / Bw(pw))];
            ErSo_Bo_n = [.. So_n.Zip(Po_n, (so, po) => Er(po) * so / Bo(po))];
        }
    }
}
