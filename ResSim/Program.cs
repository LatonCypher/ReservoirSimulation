using ScottPlot.Colormaps;
using System.Text;
{
    // =========================================================================
    // STAGE 1: GLOBAL SIMULATION PARAMETERS & FIELD CHARACTERISTICS
    // =========================================================================

    // Set root project repository path for exporting tracking diagnostics/GIFs
    folderpath = "C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\";

    // Nblocks: Spatial discretization grid cells
    // L: Maximum block boundary index for flow limits
    // M: Total number of discrete primary variables tracking
    // the porous media blocks (2 variables * 25 blocks)
    // Nwells: Count of source/sink boundaries active in the system
    int Nblocks = 25, L = Nblocks - 1, M = 2*Nblocks, Nwells = 2;

    // Thermodynamic properties, initial pressures (psi),
    // initial fluid saturation baselines,
    // residual saturation limits (Sw_r, So_r),
    // and baseline phase viscosities (cp)
    double Pinit = 3000, Sinit = 0.2, Sw_r = 0.10, So_r = 0.15,
           μw0 = 5.005, μo0 = 2, kro0 = 1.0, krw0 = 0.30, Pe = 2,
           co = 2e-5, cw = 4e-6, cr = 1e-5, bo = 2e-5, bw = 4e-10,
           Bw0 = 1.005, Bo0 = 1.4, no = 2.5, nw = 3;

    // Well operational specifications using named tuples:
    // constraints, metrics, and grid block locations
    var Producer =
        (MinPressure: 1500.0, ProdRate: 0.0,
        OilRate: 0.0, WaterRate: 0.0, Index: 0);
    var Injector =
        (MaxPressure: 4500.0, InjRate: 0.0,
        OilRate: 0.0, WaterRate: 0.0, Index: 24);

    // =========================================================================
    // STAGE 2: DATA STRUCTURE PACKING & UNPACKING (VECTOR TO PHYSICAL FIELDS)
    // =========================================================================

    // Unpacks a monolithic 1D solver array back into distinct,
    // physically meaningful array fields
    (double[], double[], double[], double[]) Unpack(double[] x)
    {
        int indx = 0;
        double[] Po = Zeros(Nblocks), Sw = Zeros(Nblocks),
                 Pwells = Zeros(Nwells), Qwells = Zeros(Nwells);

        // Extract interleaved cell variables: [Po_0, Sw_0, Po_1, Sw_1, ...]
        for (int i = 0; i < Nblocks; i++)
        {
            // Extract Oil Pressure
            Po[i] = x[indx++];
            // Extract Water Saturation
            Sw[i] = x[indx++];
        }

        // Extract interleaved well variables: [Pwf_0, Q_0, Pwf_1, Q_1, ...]
        for (int i = 0; i < Nwells; i++)
        {
            // Extract Well Bottomhole Pressure (BHP)
            Pwells[i] = x[indx++];
            // Extract Total Well Volumetric Flow Rate
            Qwells[i] = x[indx++];
        }
        return (Po, Sw, Pwells, Qwells);
    }

    // Packs separate grid-block and well residuals into a single
    // consolidated vector for Fsolve
    double[] Pack(double[] Ro, double[] Rw, double[] Rwells, double[] Rcontrol)
    {
        int indx = 0;
        double[] R_total = Zeros(M + 2*Nwells);

        // Map conservation equations sequentially to match variable indexing
        for (int i = 0; i < Nblocks; i++)
        {
            // Oil mass conservation residual
            R_total[indx++] = Ro[i];
            // Water mass conservation residual
            R_total[indx++] = Rw[i];
        }
        for (int i = 0; i < Nwells; i++)
        {
            // Wellbore mass balance residual
            R_total[indx++] = Rwells[i];
            // Well operating constraint/control residual
            R_total[indx++] = Rcontrol[i];
        }
        return R_total;
    }

    // =========================================================================
    // STAGE 3: CONSTITUTIVE EQUATIONS & PETROPHYSICAL RELATIONS
    // =========================================================================

    // Normalized and effective water saturations used to
    // calculate relative permeabilities and capillary pressures
    double Sws(double Sw) => (Sw - Sw_r)/(1 - Sw_r);
    double Swe(double Sw) => (Sw - Sw_r)/(1 - Sw_r - So_r);

    // Capillary Pressure relationships (Drainage vs Imbibition curves)
    double Pc_D(double Sw) => Pe * Pow(Sws(Sw), -0.5);
    double Pc_I(double Sw) => Pe * (Pow(Swe(Sw), -0.5) - 1);

    // Formation Volume Factors (B-factors) modeling fluid compressibility
    double Bo(double Po) => Bo0 * Exp(co*(1500 - Po));
    double Bw(double Pw) => Bw0 * Exp(cw*(2500 - Pw));

    // Pressure-dependent dynamic fluid viscosities
    double μo(double Po) => μo0 * Exp(bo*(Po - 1500));
    double μw(double Pw) => μw0 * Exp(bw*(Pw - 2500));

    // Modified Corey Brooks models evaluating relative permeability curves
    double Krw(double Sw) => krw0 * Pow(Swe(Sw), nw);
    double Kro(double So) => kro0 * Pow(1 - Swe(1 - So), no);

    // Inter-block transmissibility calculation using the
    // harmonic mean of raw cell permeabilities
    double Harmmean(double x1, double x2) => 2/(1/x1 + 1/x2);

    // Unit conversion coefficients for oilfield standard parameters
    // =================================================================
    // Transmissibility conversion factor (Darcy to Field Units)
    double alpha = 1.127e-3;
    // Geometric productivity factor conversion for radial well inflows
    double alpha_well = alpha*2*pi;
    // Reservoir volume factor (Cubic feet to Stock Tank Barrels)
    double beta = 5.615;

    // Spatial dimensions (ft), block pore volume (STB),
    // and well radius parameters
    double dt, Dx = 200, Dy = 1000, Dz = 20, Ax = Dy*Dz,
           V = Dx*Dy*Dz/beta, rw = 0.5, re, WI, WIw, WIo;

    // Heterogeneous field instantiation using a
    // normal distribution for porosity and permeability
    // Localized Porosity array
    double[] Phi = Randn(Nblocks, 0.2, 0.01);
    // Localized Permeability array (md)
    double[] K = Randn(Nblocks, 900.0, 300.0);

    // Arrays storing baseline values at historical time level (n)
    double[] Po_n, Sw_n, Pwells_n, Qwells_n, Pw_n, So_n;
    // Tracks operational control mode for each well
    // (True = Rate, False = Pressure)
    bool[] RateControl = [true, true];

    // =========================================================================
    // STAGE 4: NONLINEAR RESIDUAL FORMULATION (MASS CONSERVATION & WELL EQUATIONS)
    // =========================================================================
    double[] Residual(double[] xnp1)
    {
        double Po_up, Pw_up, So_up, Sw_up, Tr, Tw, To;

        // Unpack current iteration guess variables for evaluation
        var (Po_np1, Sw_np1, Pwells_np1, Qwells_np1) = Unpack(xnp1);

        // Calculate dependent properties using
        // capillary pressure and saturation definitions
        double[] Pw_np1 = [.. Po_np1.Zip(Sw_np1, (po, sw) => po - Pc_I(sw))],
                 So_np1 = [.. Sw_np1.Select(sw => 1 - sw)];

        double[] Rw = Zeros(Nblocks), Ro = Zeros(Nblocks),
                 Rwells = Zeros(Nwells), Rcontrol = Zeros(Nwells);

        // Grid block mass balance accumulation loop
        for (int i = 0; i < Nblocks; i++)
        {
            // Inter-block inter-flux flux logic (Left neighbor interaction)
            if (i > 0)
            {
                Tr = alpha*Ax*Harmmean(K[i-1], K[i]);

                // Single-point upstream weighting for water phase stability
                (Pw_up, Sw_up) = Pw_np1[i-1] > Pw_np1[i] ?
                    (Pw_np1[i-1], Sw_np1[i-1]) : (Pw_np1[i], Sw_np1[i]);
                Tw = Tr*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                Rw[i] += Tw*(Pw_np1[i-1] - Pw_np1[i])/Dx;

                // Single-point upstream weighting for oil phase stability
                (Po_up, So_up) = Po_np1[i-1] > Po_np1[i] ?
                    (Po_np1[i-1], So_np1[i-1]) : (Po_np1[i], So_np1[i]);
                To = Tr*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                Ro[i] += To*(Po_np1[i-1] - Po_np1[i])/Dx;
            }

            // Inter-block inter-flux flux logic (Right neighbor interaction)
            if (i < L)
            {
                Tr = alpha*Ax*Harmmean(K[i], K[i+1]);

                // Single-point upstream weighting for water phase stability
                (Pw_up, Sw_up) = Pw_np1[i+1] > Pw_np1[i] ?
                    (Pw_np1[i+1], Sw_np1[i+1]) : (Pw_np1[i], Sw_np1[i]);
                Tw = Tr*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                Rw[i] += Tw*(Pw_np1[i+1] - Pw_np1[i])/Dx;

                // Single-point upstream weighting for oil phase stability
                (Po_up, So_up) = Po_np1[i+1] > Po_np1[i] ?
                    (Po_np1[i+1], So_np1[i+1]) : (Po_np1[i], So_np1[i]);
                To = Tr*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                Ro[i] += To*(Po_np1[i+1] - Po_np1[i])/Dx;
            }

            // Time accumulation terms (implicit backward Euler method)
            Rw[i] -= V*Phi[i]*(Sw_np1[i]/Bw(Pw_np1[i]) - Sw_n[i]/Bw(Pw_n[i]))/dt;
            Ro[i] -= V*Phi[i]*(So_np1[i]/Bo(Po_np1[i]) - So_n[i]/Bo(Po_n[i]))/dt;
        }

        // Peaceman radius calculation for an isolated wellbore
        re = 0.14*Hypot(Dx, Dy);
        int idx;

        // --- Well Formulation: Producer Section ---
        Rwells[0] += Qwells_np1[0];
        idx = Producer.Index;
        // Base Well Index calculation
        WI = alpha_well*K[idx]*Dz/Log(re/rw);
        // Water mobility at well
        WIw = WI*Krw(Sw_np1[idx])/(μw(Pw_np1[idx])*Bw(Pw_np1[idx]));
        // Oil mobility at well
        WIo = WI*Kro(So_np1[idx])/(μo(Po_np1[idx])*Bo(Po_np1[idx]));
        // Calculate phase flow rates using the wellbore
        // pressure and connected grid cell properties
        Producer.WaterRate = (Pwells_np1[0] - Pw_np1[Producer.Index])*WIw;
        Producer.OilRate = (Pwells_np1[0] - Po_np1[Producer.Index])*WIo;
        // Inject sink terms directly back into the connected grid cell
        Rw[Producer.Index] += Producer.WaterRate;
        Ro[Producer.Index] += Producer.OilRate;
        // Well balance check
        Rwells[0] -= Producer.WaterRate + Producer.OilRate;
        // Evaluate dynamic constraint swapping logic
        // (Switch between Target Rate vs Minimum Allowed BHP)
        Rcontrol[0] = RateControl[0] ?
            Qwells_np1[0] - Producer.ProdRate : Pwells_np1[0] - Producer.MinPressure;

        // --- Well Formulation: Injector Section ---
        Rwells[1] += Qwells_np1[1];
        idx = Injector.Index;
        // Base Well Index calculation
        WI = alpha_well*K[idx]*Dz/Log(re/rw);
        // Injecting pure water phase
        WIw = WI*krw0/(μw(Pw_np1[idx])*Bw(Pw_np1[idx]));
        // Calculate injection rate using the
        // wellbore pressure and connected grid cell properties
        Injector.WaterRate = (Pwells_np1[1] - Pw_np1[idx])*WIw;
        // Inject source terms directly back into the connected grid cell
        Rw[idx] += Injector.WaterRate;
        // Well balance check
        Rwells[1] -= Injector.WaterRate;
        // Evaluate dynamic constraint swapping logic
        // (Switch between Target Rate vs Maximum Allowed BHP)
        Rcontrol[1] = RateControl[1] ?
            Qwells_np1[1] - Injector.InjRate : Pwells_np1[1] - Injector.MaxPressure;

        // Consolidate values back into a single vector for solver optimization loops
        return Pack(Ro, Rw, Rwells, Rcontrol);
    }

    // Mathematical mapping and linear interpolation helper utilities
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

    // chunk writer utility for formatted console output of array values
    string Write(double[] data)
    {
        var chunks = data.Chunk(8);
        List<string> sb = [];
        foreach (var chunk in chunks)
            sb.Add(string.Join(", ", chunk.Select(x => x.ToString("F2"))));
        return string.Join(", \n", sb);
    }

    // =========================================================================
    // STAGE 5: RUNTIME EXECUTION MANAGEMENT LOOP & SENSITIVITY DESIGN
    // =========================================================================
    double EndTime = 10000; // Target simulation duration limit (days)
    double delt = EndTime/300; // Time steps for reporting intervals

    // Sensitivity loop testing performance metrics across varying injection-production rates
    for (int rate = 200; rate <= 500; rate += 100)
    {
        dt = 0.01; // Reset to a safe initial time step value

        // Initialize state vectors for a new sensitivity run
        Po_n = Repmat(Pinit, Nblocks); Sw_n = Repmat(Sinit, Nblocks);
        Pwells_n = Repmat(Pinit, Nwells); Qwells_n = Zeros(Nwells);
        Pw_n = [.. Po_n.Zip(Sw_n, (po, sw) => po - Pc_I(sw))];
        So_n = [.. Sw_n.Select(sw => 1 - sw)];

        // Initialize historical data tracking containers for plotting and reporting
        List<double[]> P = [Po_n], S = [Sw_n];
        List<double> Time = [0.0], WaterCut = [0.0], SweepEff = [0.0],
            ProdRate = [Qwells_n[0]], InjRate = [Qwells_n[1]],
            ProdPwf = [Pwells_n[0]], InjPwf = [Pwells_n[1]];

        // Define specific production/injection targets for this scenario run
        Producer.ProdRate = -rate;
        Injector.InjRate = rate;

        // Initialize Plot Subplots for Real-Time Visual Feedback Diagnostics
        Subplot(7, 4, [0, 1, 4, 5]);
        var Pbhp = Plot([0], [0], "r", 2);
        Axis([0, EndTime, 0, Injector.MaxPressure*1.1]);
        Title("Producer BHP");

        Subplot(7, 4, [8, 9, 12, 13]);
        var Prate = Plot([0], [0], "r", 2);
        Axis([0, EndTime, 0, Producer.ProdRate*1.1]);
        Title("Producer Rate");

        Subplot(7, 4, [16, 17, 20, 21]);
        var Pbsw = Plot([0], [0], "r", 2);
        Axis([0, EndTime, 0, 105]);
        Title("Producer WaterCut");

        Subplot(7, 4, [2, 3, 6, 7]);
        var Ibhp = Plot([0], [0], "b", 2);
        Axis([0, EndTime, 0, Injector.MaxPressure*1.1]);
        Title("Injector BHP");

        Subplot(7, 4, [10, 11, 14, 15]);
        var Irate = Plot([0], [0], "b", 2);
        Axis([0, EndTime, 0, Injector.InjRate*1.1]);
        Title("Injector Rate");

        Subplot(7, 4, [18, 19, 22, 23]);
        var Iswp = Plot([0], [0], "b", 2);
        Axis([0, EndTime, 0, 105]);
        Title("Injector Sweep Efficiency");

        // Spatial Profile Schematic setup for the moving Water Saturation Front
        Subplot(7, 4, [24, 25, 26, 27]);
        double[] xpf = [0.3, 0.7], ypf = Linspace(0.3, 0.7, 5);
        var (xperf, yperf) = Meshgrid(xpf, ypf);
        Rectangle([Producer.Index + 0.45, 0.2, 0.1, 1.6]); HoldOn();
        Plot(Producer.Index + xperf, yperf, "k", 2);
        Rectangle([Injector.Index + 0.45, 0.2, 0.1, 1.6]);
        Plot(Injector.Index + xperf, yperf, "k", 2);
        double[] xres = [0, 0, Nblocks, Nblocks], yres = [0, 1, 1, 0];
        Fill(xres, yres, [1, 0, 0], 0.5);
        xpf = Linspace(0.5, Nblocks - 0.5, Nblocks);
        ypf = Repmat(Sinit, Nblocks + 2);
        double[] xplot = [0, 0, .. xpf, Nblocks, Nblocks];
        double[] yplot = [0, .. ypf, 0];
        var Water = Fill(xplot, yplot, [0, 0, 1], 0.5, 0);
        Axis([0, Nblocks, 0, 2.5]);
        Title("Water Saturation Front");
        HoldOff();

        // Package the initial state vector and configure the nonlinear solver options
        double[] xs = Pack(Po_n, Sw_n, Pwells_n, Qwells_n), xn;
        var opts = SolverSet(MaxIter: 10, AbsTol: 1e-6, UseParallel: true);

        Console.WriteLine($"""
                    ======================================================================
                               Starting simulation for Rate = {rate} STB/day

                    Time: 
                    {0:F2} days

                    Producer BHP: 
                    {Pinit:F2} psi

                    Injector BHP: 
                    {Pinit:F2} psi
                    
                    Pressure: 
                    {Write(Po_n)}
                    
                    Saturation:
                    {Write(Sw_n)}



                    """);

        // =========================================================================
        // STAGE 6: CORE TIME-STEPPING INTERACTION LOOP (NEWTON RAPHSON + CONTROLS)
        // =========================================================================
        double[] displaytimes = Linspace(0, EndTime, 21);
        double printtime = 0;
        while (Time.Last() < EndTime)
        {
            if (displaytimes.Any(t =>
            (Time.Last() - t)*(Time.Last() + dt - t) < 0))
            {
                opts.Display = true;
                printtime = Array.Find(displaytimes,
                    t => (Time.Last() - t)*(Time.Last() + dt - t) < 0);
            }
            else
            {
                opts.Display = false;
            }
            if(opts.Display)
                Console.WriteLine($"""


                    Time: 
                    {printtime:F2} days
                    """);

            // Call the Newton-Raphson nonlinear solver to
            // find the next state solution
            xn = [.. Fsolve(Residual, xs, opts)];

            // Check convergence. If non-converged,
            // chop the time step (time-step cuts) and retry.
            if (!opts.ans.IsConverged)
            {
                dt = 0.25*dt;
                Console.WriteLine("""
                                ================================================
                                           Rejected (Non-Convergence)
                                ================================================
                                """);
                continue;
            }

            // Unpack solution values to evaluate operational constraint validations
            var (Po_s, Sw_s, Pwells_s, Qwells_s) = Unpack(xn);

            // Validate Producer constraints:
            // switch to BHP control if pressure falls below minimum limits
            if (Pwells_s[0] < Producer.MinPressure)
            {
                RateControl[0] = Pwells_s[0] > Producer.MinPressure;
                Console.WriteLine("""
                                ================================================
                                      Rejected (Minimum Pressure Violated) 
                                ================================================
                                """);
                continue;
            }
            // Validate Injector constraints:
            // switch to BHP control if pressure exceeds fracturing limits
            if (Pwells_s[1] > Injector.MaxPressure)
            {
                RateControl[1] = Pwells_s[1] < Injector.MaxPressure;
                Console.WriteLine("""
                                ================================================
                                      Rejected (Maximum Pressure Violated) 
                                ================================================
                                """);
                continue;
            }

            // Accept the validated time step and update internal tracking states
            (Po_n, Sw_n, Pwells_n, Qwells_n) = (Po_s, Sw_s, Pwells_s, Qwells_s);
            Pw_n = [.. Po_n.Zip(Sw_n, (po, sw) => po - Pc_I(sw))];
            So_n = [.. Sw_n.Select(sw => 1 - sw)];

            // Log verified parameters to performance history arrays
            P.Add(Po_n); S.Add(Sw_n);
            ProdPwf.Add(Pwells_n[0]); InjPwf.Add(Pwells_n[1]);
            ProdRate.Add(Qwells_n[0]); InjRate.Add(Qwells_n[1]);
            WaterCut.Add(Producer.WaterRate*100/(Producer.WaterRate + Producer.OilRate));
            SweepEff.Add((Sw_n.Sum()/Nblocks - Sinit)*100/(1 - Sinit));
            Time.Add(Time.Last() + dt);

            // Adaptive Time-Stepping Logic:
            // scale dt up if convergence is fast, scale down if slow
            if (opts.ans.Iter < 4) dt = 1.25*dt;
            if (opts.ans.Iter > 8) dt = 0.5*dt;
            if (dt < 1e-5) throw new Exception("time step is too small");

            // Set current solution vector as the next time step initialization guess (xs)
            xs = [.. xn];


            if (opts.Display)
            {
                Console.WriteLine($"""
                    Producer BHP: 
                    {interps(Time, ProdPwf, printtime):F2} psi

                    Injector BHP: 
                    {interps(Time, InjPwf, printtime):F2} psi
                    
                    Pressure: 
                    {Write(interpa(Time, P, printtime))}
                    
                    Saturation:
                    {Write(interpa(Time, S, printtime))}



                    """);
            }
        }

        // =========================================================================
        // STAGE 7: POST-PROCESSING VISUALIZATION & ANIMATION EXPORT
        // =========================================================================

        // Dynamic closure map function mapping raw state histories
        // directly onto figure frames
        byte[] Animfun(int i)
        {
            // Sync time dimensions to interpolation index targets
            Pbhp.Xdata = Vcart(Pbhp.Xdata, i*delt);
            Pbhp.Ydata = Vcart(Pbhp.Ydata, interps(Time, ProdPwf, i*delt));

            Prate.Xdata = Pbhp.Xdata;
            Prate.Ydata = Vcart(Prate.Ydata, interps(Time, ProdRate, i*delt));

            Pbsw.Xdata = Pbhp.Xdata;
            Pbsw.Ydata = Vcart(Pbsw.Ydata, interps(Time, WaterCut, i*delt));

            Ibhp.Xdata = Pbhp.Xdata;
            Ibhp.Ydata = Vcart(Ibhp.Ydata, interps(Time, InjPwf, i*delt));

            Irate.Xdata = Pbhp.Xdata;
            Irate.Ydata = Vcart(Irate.Ydata, interps(Time, InjRate, i*delt));

            Iswp.Xdata = Pbhp.Xdata;
            Iswp.Ydata = Vcart(Iswp.Ydata, interps(Time, SweepEff, i*delt));

            // Generate spatial mapping visualization showing fluid front changes over time
            Sw_n = interpa(Time, S, i*delt);
            yplot = [0, Sw_n.First(), .. Sw_n, Sw_n.Last(), 0];
            Water.Ydata = yplot;

            // Capture rendered visual asset frame dimensions
            return GetFrame(800, 1000);
        }

        // Compile and output the processed timeline visualization to an animated GIF
        Console.WriteLine(DateTime.Now);
        AnimationMaker(Animfun, $"1D_WaterFlooding_{rate}.gif", 30, 300);
        Console.WriteLine(DateTime.Now);

        // Clean up graphics devices to free memory hooks for the next loop run
        CloseFig();
    }
}

{
    // 2D - 2Phase
    folderpath = "C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\";
    int Nx = 25, Ny = 5, Lx = Nx - 1, Ly = Ny - 1, M = 2*Nx*Ny, Nwells = 2, WellIndex;
    double Pinit = 3000, Sinit = 0.2, Sw_r = 0.10, So_r = 0.15,
           μw0 = 5.005, μo0 = 2, kro0 = 1.0, krw0 = 0.30, Pe = 2,
           co = 2e-5, cw = 4e-6, cr = 1e-5, bo = 2e-5, bw = 4e-10,
           Bw0 = 1.005, Bo0 = 1.4, no = 2.5, nw = 3;
    var Producer = (MinPressure: 1500.0, ProdRate: 0.0, OilRate: 0.0, WaterRate: 0.0, I: 0, J: 2);
    var Injector = (MaxPressure: 4500.0, InjRate: 0.0, OilRate: 0.0, WaterRate: 0.0, I: 24, J: 2);

    (double[], double[], double[], double[]) Unpack(double[] x)
    {
        int indx = 0;
        double[] Po = Zeros(Nx*Ny), Sw = Zeros(Nx*Ny),
                 Pwells = Zeros(Nwells), Qwells = Zeros(Nwells);
        for (int i = 0; i < Nx*Ny; i++)
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
        double[] R_total = Zeros(M + 2*Nwells);
        for (int i = 0; i < Nx*Ny; i++)
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
    double Sws(double Sw) => (Sw - Sw_r)/(1 - Sw_r);
    double Swe(double Sw) => (Sw - Sw_r)/(1 - Sw_r - So_r);
    double Pc_D(double Sw) => Pe * Pow(Sws(Sw), -0.5);
    double Pc_I(double Sw) => Pe * (Pow(Swe(Sw), -0.5) - 1);
    double Bo(double Po) => Bo0*Exp(co*(2000 - Po));
    double Bw(double Pw) => Bw0*Exp(cw*(2500 - Pw));
    double μo(double Po) => μo0*Exp(bo*(Po - 2000));
    double μw(double Pw) => μw0*Exp(bw*(Pw - 2500));
    double Krw(double Sw) => krw0 * Pow(Swe(Sw), nw);
    double Kro(double So) => kro0 * Pow(1 - Swe(1 - So), no);
    double Harmmean(double x1, double x2) => 2/(1/x1 + 1/x2);

    // Define conversion constants
    double alpha = 1.127e-3;               // Darcy to Field units factor
    double alpha_well = alpha*2*pi;        // Darcy to Field units factor for wells
    double beta = 5.615;                   // ft3 to bbl conversion factor

    double dt, Dx = 200, Dy = 200, Dz = 20, Ax = Dy*Dz, Ay = Dx*Dz,
           V = Dx*Dy*Dz/beta, rw = 0.5, re, WI, WIw, WIo;

    double[] Phi = Randn(Nx*Ny, 0.2, 0.01);      // Porosity
    double[] K = Randn(Nx*Ny, 900.0, 300.0);     // Permeability

    double[] Po_n, Sw_n, Pwells_n, Qwells_n, Pw_n, So_n;
    bool[] RateControl = [true, true];

    double[] Residual(double[] xnp1)
    {
        double Po_up, Pw_up, So_up, Sw_up, Tw, To;
        var (Po_np1, Sw_np1, Pwells_np1, Qwells_np1) = Unpack(xnp1);
        double[] Pw_np1 = [.. Po_np1.Zip(Sw_np1, (po, sw) => po - Pc_I(sw))],
            So_np1 = [.. Sw_np1.Select(sw => 1 - sw)];

        double[] Rw = Zeros(Nx*Ny), Ro = Zeros(Nx*Ny),
            Rwells = Zeros(Nwells), Rcontrol = Zeros(Nwells);

        for (int m = 0; m < Nx*Ny; m++)
        {
            var (j, i) = DivRem(m, Nx);
            Rw[m] -= V*Phi[m]*(Sw_np1[m]/Bw(Pw_np1[m]) - Sw_n[m]/Bw(Pw_n[m]))/dt;
            Ro[m] -= V*Phi[m]*(So_np1[m]/Bo(Po_np1[m]) - So_n[m]/Bo(Po_n[m]))/dt;

            if (i > 0)
            {
                (Pw_up, Sw_up) = Pw_np1[m-1] > Pw_np1[m] ? (Pw_np1[m-1], Sw_np1[m-1]) : (Pw_np1[m], Sw_np1[m]);
                Tw = alpha*Ax*Harmmean(K[m-1], K[m])*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                Rw[m] += Tw*(Pw_np1[m-1] - Pw_np1[m])/Dx;

                (Po_up, So_up) = Po_np1[m-1] > Po_np1[m] ? (Po_np1[m-1], So_np1[m-1]) : (Po_np1[m], So_np1[m]);
                To = alpha*Ax*Harmmean(K[m-1], K[m])*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                Ro[m] += To*(Po_np1[m-1] - Po_np1[m])/Dx;
            }
            if (i < Lx)
            {
                (Pw_up, Sw_up) = Pw_np1[m+1] > Pw_np1[m] ? (Pw_np1[m+1], Sw_np1[m+1]) : (Pw_np1[m], Sw_np1[m]);
                Tw = alpha*Ax*Harmmean(K[m], K[m+1])*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                Rw[m] += Tw*(Pw_np1[m+1] - Pw_np1[m])/Dx;

                (Po_up, So_up) = Po_np1[m+1] > Po_np1[m] ? (Po_np1[m+1], So_np1[m+1]) : (Po_np1[m], So_np1[m]);
                To = alpha*Ax*Harmmean(K[m], K[m+1])*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                Ro[m] += To*(Po_np1[m+1] - Po_np1[m])/Dx;
            }

            if (j > 0)
            {
                (Pw_up, Sw_up) = Pw_np1[m-Nx] > Pw_np1[m] ? (Pw_np1[m-Nx], Sw_np1[m-Nx]) : (Pw_np1[m], Sw_np1[m]);
                Tw = alpha*Ay*Harmmean(K[m-Nx], K[m])*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                Rw[m] += Tw*(Pw_np1[m-Nx] - Pw_np1[m])/Dy;

                (Po_up, So_up) = Po_np1[m-Nx] > Po_np1[m] ? (Po_np1[m-Nx], So_np1[m-Nx]) : (Po_np1[m], So_np1[m]);
                To = alpha*Ay*Harmmean(K[m-Nx], K[m])*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                Ro[m] += To*(Po_np1[m-Nx] - Po_np1[m])/Dy;
            }
            if (j < Ly)
            {
                (Pw_up, Sw_up) = Pw_np1[m+Nx] > Pw_np1[m] ? (Pw_np1[m+Nx], Sw_np1[m+Nx]) : (Pw_np1[m], Sw_np1[m]);
                Tw = alpha*Ay*Harmmean(K[m], K[m+Nx])*Krw(Sw_up)/(μw(Pw_up)*Bw(Pw_up));
                Rw[m] += Tw*(Pw_np1[m+Nx] - Pw_np1[m])/Dy;
                (Po_up, So_up) = Po_np1[m+Nx] > Po_np1[m] ? (Po_np1[m+Nx], So_np1[m+Nx]) : (Po_np1[m], So_np1[m]);
                To = alpha*Ay*Harmmean(K[m], K[m+Nx])*Kro(So_up)/(μo(Po_up)*Bo(Po_up));
                Ro[m] += To*(Po_np1[m+Nx] - Po_np1[m])/Dy;
            }
        }

        // Producer
        Rwells[0] += Qwells_np1[0];
        WellIndex = Producer.I + Producer.J*Nx;
        re = 0.14*Hypot(Dx, Dy); WI = alpha_well*K[WellIndex]*Dz/Log(re/rw);
        WIw = WI*Krw(Sw_np1[WellIndex])/(μw(Pw_np1[WellIndex])*Bw(Pw_np1[WellIndex]));
        WIo = WI*Kro(So_np1[WellIndex])/(μo(Po_np1[WellIndex])*Bo(Po_np1[WellIndex]));
        Producer.WaterRate = (Pwells_np1[0] - Pw_np1[WellIndex])*WIw;
        Producer.OilRate = (Pwells_np1[0] - Po_np1[WellIndex])*WIo;
        Rw[WellIndex] += Producer.WaterRate;
        Ro[WellIndex] += Producer.OilRate;
        Rwells[0] -= Producer.WaterRate + Producer.OilRate;
        Rcontrol[0] = RateControl[0] ? Qwells_np1[0] - Producer.ProdRate : Pwells_np1[0] - Producer.MinPressure;

        // Injector
        Rwells[1] += Qwells_np1[1];
        WellIndex = Injector.I + Injector.J*Nx;
        re = 0.14*Hypot(Dx, Dy); WI = alpha_well*K[WellIndex]*Dz/Log(re/rw);
        WIw = WI*krw0/(μw(Pw_np1[WellIndex])*Bw(Pw_np1[WellIndex]));
        Rw[WellIndex] += (Pwells_np1[1] - Pw_np1[WellIndex])*WIw;
        Rwells[1] -= (Pwells_np1[1] - Pw_np1[WellIndex])*WIw;
        Rcontrol[1] = RateControl[1] ? Qwells_np1[1] - Injector.InjRate : Pwells_np1[1] - Injector.MaxPressure;

        return Pack(Ro, Rw, Rwells, Rcontrol);
    }

    double EndTime = 10000;
    double delt = EndTime/300;
    for (int rate = 200; rate <= 500; rate += 100)
    {
        dt = 0.01;
        Po_n = Repmat(Pinit, Nx*Ny); Sw_n = Repmat(Sinit, Nx*Ny);
        Pwells_n = Repmat(Pinit, Nwells); Qwells_n = Zeros(Nwells);
        Pw_n = [.. Po_n.Zip(Sw_n, (po, sw) => po - Pc_I(sw))];
        So_n = [.. Sw_n.Select(sw => 1 - sw)];

        List<double[]> P = [Po_n], S = [Sw_n];
        List<double> Time = [0.0], WaterCut = [0.0], SweepEff = [0.0],
                     ProdRate = [Qwells_n[0]], InjRate = [Qwells_n[1]],
                     ProdPwf = [Pwells_n[0]], InjPwf = [Pwells_n[1]];

        Producer.ProdRate = -rate;
        Injector.InjRate = rate;

        // Plot of Initial State
        Subplot(8, 4, [0, 1, 4, 5]);
        var Pbhp = Plot([0], [0], "r", 2);
        Axis([0, EndTime, 0, Injector.MaxPressure*1.1]);
        Title("Producer BHP");

        Subplot(8, 4, [8, 9, 12, 13]);
        var Prate = Plot([0], [0], "r", 2);
        Axis([0, EndTime, 0, Producer.ProdRate*1.1]);
        Title("Producer Rate");

        Subplot(8, 4, [16, 17, 20, 21]);
        var Pbsw = Plot([0], [0], "r", 2);
        Axis([0, EndTime, 0, 105]);
        Title("Producer WaterCut");

        Subplot(8, 4, [2, 3, 6, 7]);
        var Ibhp = Plot([0], [0], "b", 2);
        Axis([0, EndTime, 0, Injector.MaxPressure*1.1]);
        Title("Injector BHP");

        Subplot(8, 4, [10, 11, 14, 15]);
        var Irate = Plot([0], [0], "b", 2);
        Axis([0, EndTime, 0, Injector.InjRate*1.1]);
        Title("Injector Rate");

        Subplot(8, 4, [18, 19, 22, 23]);
        var Iswp = Plot([0], [0], "b", 2);
        Axis([0, EndTime, 0, 105]);
        Title("Injector Sweep Efficiency");


        Subplot(8, 4, [24, 25, 26, 27, 28, 29, 30, 31]);
        RectHandle[,] Water = new RectHandle[Nx, Ny];
        HoldOn();
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                Water[i, j] = Rectangle([i, j, 1, 1]);
                Water[i, j].FillAlpha = 0.5;
                Water[i, j].FillColor = [1-Sinit, 0, Sinit];
                Water[i, j].LineAlpha = 0.3;
            }
        }
        Axis([0, Nx, 0, Ny]);
        Title("Water Saturation Front");
        HoldOff();

        double[] xs = Pack(Po_n, Sw_n, Pwells_n, Qwells_n), xn;
        var opts = SolverSet(Display: true, MaxIter: 10, AbsTol: 1e-6, UseParallel: true);

        while (Time.Last() < EndTime)
        {
            xn = [.. Fsolve(Residual, xs, opts)];
            if (!opts.ans.IsConverged)
            {
                // if not converged, reduce time step and repeat
                dt = 0.25*dt; continue;
            }
            var (Po_s, Sw_s, Pwells_s, Qwells_s) = Unpack(xn);
            if (Pwells_s[0] < Producer.MinPressure)
            {
                // if min pressure violated, change to pressure control and repeat
                RateControl[0] = Pwells_s[0] > Producer.MinPressure; continue;
            }
            if (Pwells_s[1] > Injector.MaxPressure)
            {
                // if max pressure violated, change to pressure control and repeat
                RateControl[1] = Pwells_s[1] < Injector.MaxPressure; continue;
            }
            (Po_n, Sw_n, Pwells_n, Qwells_n) = (Po_s, Sw_s, Pwells_s, Qwells_s);
            Pw_n = [.. Po_n.Zip(Sw_n, (po, sw) => po - Pc_I(sw))];
            So_n = [.. Sw_n.Select(sw => 1 - sw)];
            P.Add(Po_n); S.Add(Sw_n);
            ProdPwf.Add(Pwells_n[0]);
            InjPwf.Add(Pwells_n[1]);
            ProdRate.Add(Qwells_n[0]);
            InjRate.Add(Qwells_n[1]);
            WaterCut.Add(Producer.WaterRate*100/(Producer.WaterRate + Producer.OilRate));
            SweepEff.Add((Sw_n.Sum()/(Nx*Ny) - Sinit)*100/(1 - Sinit));
            Time.Add(Time.Last() + dt);

            if (opts.ans.Iter < 6) dt = 1.25*dt;
            if (opts.ans.Iter > 9) dt = 0.5*dt;
            if (dt < 1e-5) throw new Exception("time step is too small");

            // call guess function
            xs = [.. xn];

            // write to console
            Console.WriteLine($"Solution at time t = {Time.Last():F3}");
            //Console.WriteLine("\n Pressure = ");
            //WriteArray(Po_n);

            //Console.WriteLine("\n Saturation = ");
            //WriteArray(Sw_n);
            //Console.WriteLine("\n\n\n\n");
        }

        // Post Processing
        byte[] Animfun(int framenum)
        {
            double t = framenum*delt;
            Pbhp.Xdata = Vcart(Pbhp.Xdata, t);
            Pbhp.Ydata = Vcart(Pbhp.Ydata, interps(Time, ProdPwf, t));

            Prate.Xdata = Pbhp.Xdata;
            Prate.Ydata = Vcart(Prate.Ydata, interps(Time, ProdRate, t));

            Pbsw.Xdata = Pbhp.Xdata;
            Pbsw.Ydata = Vcart(Pbsw.Ydata, interps(Time, WaterCut, t));

            Ibhp.Xdata = Pbhp.Xdata;
            Ibhp.Ydata = Vcart(Ibhp.Ydata, interps(Time, InjPwf, t));

            Irate.Xdata = Pbhp.Xdata;
            Irate.Ydata = Vcart(Irate.Ydata, interps(Time, InjRate, t));

            Iswp.Xdata = Pbhp.Xdata;
            Iswp.Ydata = Vcart(Iswp.Ydata, interps(Time, SweepEff, t));

            Sw_n = interpa(Time, S, t);
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                    Water[i, j].FillColor = [1-Sw_n[i + Nx*j], 0, Sw_n[i + Nx*j]];
            }
            return GetFrame(800, 1000);
        }
        Console.WriteLine(DateTime.Now);
        AnimationMaker(Animfun, $"2D_WaterFlooding_{rate}.gif", 30, 300);
        Console.WriteLine(DateTime.Now);
        CloseFig();
        Console.WriteLine("=======================================================");
    }
}