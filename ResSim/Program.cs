using EclipseDeckReader;
using System.Diagnostics;
{
    //// 1. Grid Specifications (10x10x3)
    //int nx = 50, ny = 1, nz = 1;
    //int ngrids = nx * ny * nz;

    //double[] dx = Repmat(10, ngrids);
    //double[] dy = Repmat(10, ngrids);
    //double[] dz = Repmat(20, ngrids);
    //double[] z = [..Enumerable.Range(0, ngrids).
    //                Select(i=>10 + i/(nx*ny)*20)];          
    //double[] perm = Repmat(500, ngrids);                   // 100 mD
    //double[] phi = Repmat(0.20, ngrids);                   // 20% porosity

    //// 2. Instantiate Wells
    //List<Well> wells =
    //[
    //    // Injector at Block 0
    //    new Well(WellType.Producer, "WP3051X", 0.5, 0, 1500, 6000, 49, 0, [0,0], [0], [100]),
    //    // Producer at Block 9
    //    new Well(WellType.Injector, "WP5032X", 0.5, 0, 1500, 6000, 0, 0, [0,0], [0], [100])
    //];

    //// 3. Physical Parameters mapped precisely to your signature
    //double peow = 0;         // Disabled capillary entry pressure for the baseline check
    //double pw_woc = 2000.0;  // Reference pressure at Water-Oil Contact
    //double z_woc = 25.0;    // Water-Oil Contact
    //double mult_z = 1.0;     // Vertical permeability multiplier anisotropy
    //double sw_r = 0.2;       // Residual water saturation
    //double so_r = 0.2;       // Residual oil saturation
    //double bo0 = 1.1;        // Initial Oil FVF (Formation Volume Factor)
    //double bw0 = 1.038;      // Initial Water FVF
    //double muo0 = 2.0;       // Oil Viscosity (2.0 cP)
    //double muw0 = 0.118;     // Water Viscosity (1.0 cP)
    //double gamo0 = 0.356;    // Oil gravity gradient (psi/ft)
    //double gamw0 = 0.433;    // Water gravity gradient (psi/ft)
    //double krw0 = 0.4;       // Max water relative permeability endpoint
    //double kro0 = 0.8;       // Max oil relative permeability endpoint

    //// Setting compressibilities to zero isolates purely incompressible flow dynamics
    //double co = 8e-6;
    //double cw = 6e-6;
    //double cr = 4e-6;

    //double bo = 3e-7;        // Viscosity exponential factors (Oil)
    //double bw = 0.0;         // Viscosity exponential factors (Water)
    //double nw = 1.0;         // Corey relative permeability exponent (Water)
    //double no = 1.0;         // Corey relative permeability exponent (Oil)
    //double pb = 1500.0;      // Bubble point pressure
    //double pref = 6000;      // Reference system pressure


    //// 4. Run for 30 days broken into intervals
    //double[] resultTime = Linspace(0, 2000, 2);

    //// 5. Initialization via exact signature mapping
    //var reservoir1 = new OilWaterReservoir(
    //    nx, ny, nz, perm, phi, dx, dy, dz, z, peow, pw_woc,
    //    z_woc, mult_z, sw_r, so_r, bo0, bw0, muo0, muw0, gamo0, gamw0,
    //    krw0, kro0, co, cw, cr, bo, bw, nw, no, pb, pref, wells
    //);
    //reservoir1.Initialize();
    //tic();
    //Console.WriteLine("Executing Core Simulation Test Loop...");
    //reservoir1.Simulate2Phase(resultTime, wells);
    //Console.WriteLine($"Simulation completed successfully without crashing in {toc():F2} seconds!");
    //reservoir1.ExportParaView("C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\RunTest1");

    //// 5. Initialization via exact signature mapping
    //var reservoir2 = new OilWaterReservoirNoGC(
    //    nx, ny, nz, perm, phi, dx, dy, dz, z, peow, pw_woc,
    //    z_woc, mult_z, sw_r, so_r, bo0, bw0, muo0, muw0, gamo0, gamw0,
    //    krw0, kro0, co, cw, cr, bo, bw, nw, no, pb, pref, wells
    //);
    //reservoir2.Initialize();
    //tic();
    //Console.WriteLine("Executing Core Simulation Test Loop...");
    //reservoir2.Simulate2Phase(resultTime, wells);
    //Console.WriteLine($"Simulation completed successfully without crashing in {toc():F2} seconds!");
    //reservoir2.ExportParaView("C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\RunTest2");


    //Console.WriteLine($"Number of function calls = {reservoir1.funcall:F0}!");
    //Console.WriteLine($"Number of function calls = {reservoir2.funcall:F0}!");
}

{
    //// 1. Define the path to your master ECLIPSE file (.DATA or .GRDECL)
    //// If testing with the Norne deck, point it to your local file path
    //string masterDeckPath = @"C:\Users\lateef.a.kareem\Documents\GitHub\opm-data\norne\NORNE_ATW2013.DATA";

    //Console.WriteLine("==================================================");
    //Console.WriteLine("        LAUNCHING ECLIPSE INGESTION PIPELINE     ");
    //Console.WriteLine("==================================================");

    //if (!File.Exists(masterDeckPath))
    //{
    //    Console.ForegroundColor = ConsoleColor.Red;
    //    Console.WriteLine($"[Error] Could not locate file: {masterDeckPath}");
    //    Console.ResetColor();
    //    return;
    //}

    //// 2. Initialize your reader and track performance
    //Stopwatch timer = Stopwatch.StartNew();

    //try
    //{
    //    DeckReader reader = new();

    //    Console.WriteLine($"Reading and parsing deck stream...");
    //    EclipseDataDeck loadedDeck = reader.LoadDeck(masterDeckPath);

    //    timer.Stop();

    //    // 3. Print out the success metrics and structural logs
    //    loadedDeck.PrintDeckSummary();
    //    Console.WriteLine($"\n✓ Core parsing engine finished in: {timer.ElapsedMilliseconds} ms");

    //    // 4. Verification Check: Confirm arrays aren't empty before handing to math loops
    //    loadedDeck.VerifyArrayIntegrity();
    //}
    //catch (Exception ex)
    //{
    //    Console.ForegroundColor = ConsoleColor.Red;
    //    Console.WriteLine($"\n[Critical Ingestion Failure]: {ex.Message}");
    //    if (ex.InnerException != null)
    //        Console.WriteLine($" ↳ Detail: {ex.InnerException.Message}");
    //    Console.ResetColor();
    //}

    //Console.WriteLine("\nPipeline idle. Press any key to close context window...");
    //Console.ReadKey();
}

static int[] ReadActnumFile(string filePath)
{
    var intList = new List<int>(8405); // Presize with your expected 41x41x5 grid dimension

    // Read line-by-line to avoid loading giant buffers all at once
    foreach (string line in File.ReadLines(filePath))
    {
        if (string.IsNullOrWhiteSpace(line)) continue;

        ReadOnlySpan<char> span = line.AsSpan().Trim();

        // 1. Strip out the header metadata if present
        if (span.StartsWith("["))
        {
            int closeBracketIdx = span.IndexOf(']');
            if (closeBracketIdx != -1)
            {
                span = span.Slice(closeBracketIdx + 1).Trim();
            }
        }

        // 2. Extract elements separated by spaces
        while (span.Length > 0)
        {
            int spaceIdx = span.IndexOf(' ');

            if (spaceIdx == -1)
            {
                // Last element on the line
                if (int.TryParse(span, out int lastVal))
                {
                    intList.Add(lastVal);
                }
                break;
            }

            ReadOnlySpan<char> token = span.Slice(0, spaceIdx);
            if (int.TryParse(token, out int val))
            {
                intList.Add(val);
            }

            // Advance past the processed token and continuous spaces
            span = span.Slice(spaceIdx + 1).TrimStart();
        }
    }

    return intList.ToArray();
}
{   // Eclipse Style Input

    // DIMENS
    int nx = 46, ny = 112, nz = 5;

    // GRID
    double[] dx = Repmat(25, nx*ny*nz), dy = Repmat(25, nx*ny*nz),
        dz = Repmat(10, nx*ny*nz), zTop = Repmat(0, nx*ny);
    double[] phi = Repmat(0.2, nx * ny * nz),
        perm = Repmat(1000, nx * ny * nz);
    double mult_z = 0.2;
    int[] Actnum = ReadActnumFile("C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\ResSim\\actnum.txt");// read(int)

    // PVTW
    double pref_w = 6000, bw0 = 1, cw = 8e-7, μw0 = 0.3, bw = 1e-10;

    // PVDO (Pressure, FVF, Muo) 
    Matrix pvdo = new double[,]
    {
        { 14,       1.0001,   0.9999 },
        { 100,      1.0000,   1.0000 },
        { 1000,     0.9999,   1.0001 },
        { 2000,     0.9998,   1.0002 },
        { 3000,     0.9997,   1.0003 },
        { 4000,     0.9996,   1.0004 },
        { 5000,     0.9995,   1.0005 },
        { 6000,     0.9994,   1.0006 },
        { 7000,     0.9993,   1.0007 },
        { 8000,     0.9992,   1.0008 },
        { 9000,     0.9991,   1.0009 },
        { 10000,    0.9990,   1.0010 }
    };

    // ROCK
    double pref_r = 6000, cr = 5e-6;

    // SWOF (Sw, Krw, Kro, Pcwo)
    Matrix swof = new double[,]
    {
        { 0.200,  0.0000,  1.0000,  0.0 },
        { 0.250,  0.0069,  0.8403,  0.0 },
        { 0.300,  0.0278,  0.6944,  0.0 },
        { 0.350,  0.0625,  0.5625,  0.0 },
        { 0.400,  0.1111,  0.4444,  0.0 },
        { 0.450,  0.1736,  0.3403,  0.0 },
        { 0.500,  0.2500,  0.2500,  0.0 },
        { 0.550,  0.3403,  0.1736,  0.0 },
        { 0.600,  0.4444,  0.1111,  0.0 },
        { 0.650,  0.5625,  0.0625,  0.0 },
        { 0.700,  0.6944,  0.0278,  0.0 },
        { 0.750,  0.8403,  0.0069,  0.0 },
        { 0.800,  1.0000,  0.0000,  0.0 }
    };

    //

    // DENSITY (Oil density, water density) 
    double ρo0 = 43.68, ρw0 = 62.43;

    // EQUIL     
    double datum = 0, pdatun = 3000, z_woc = 50, pcwoc = 0;

    // WELL
    List<Well> wells =
    [
        // Injector at Block 0
        new Well(WellType.Producer, "WP1", 0.5, 0, 1500, 9000, 36, 87, [0,0], [0, 100, 200, 300, 400], [0, 300, 600, 900, 1200]),
        // Producer at Block 9
        new Well(WellType.Producer, "WP2", 0.5, 0, 1500, 9000, 12, 87, [0,0], [0, 100, 200, 300, 400], [0, 300, 600, 900, 1200]),
        // Producer at Block 9
        new Well(WellType.Injector, "WI1", 0.5, 0, 1500, 9000, 17, 15, [4,4], [0, 100, 200, 300, 400], [0, 700, 1400, 2100, 2800]),
    ];

    //--AQUID              I1   I2    J1   J2     K1    K2                   FACE          CNCT_EFF
    Aquifer Aquifer = new([0, nx-1], [0, ny-1], [nz-1, nz-1], FlowDirection.Kplus, 3050, 0);

    var reservoir = new OilWaterReservoir(
            // DIMENS
            nx, ny, nz,

            // GRID
            dx, dy, dz, zTop, perm, phi, mult_z,

            // PVTW
            pref_w, bw0, cw, μw0, bw,

            // PVDO (Pressure, FVF, Muo) 
            pvdo,

            // ROCK
            pref_r, cr,

            // SWOF (Sw, Krw, Kro, Pcwo)
            swof,

            // DENSITY (Oil density, water density) 
            ρo0, ρw0,

            // EQUIL     
            datum, pdatun, z_woc, pcwoc,

            // WELL
            wells,

            // AQUIFER
            Aquifer,

            // ACTNUM
            Actnum
    );
    reservoir.Initialize();
    double[] resultTime = Linspace(0, 10000, 2);


    tic();
    Console.WriteLine("Executing Core Simulation Test Loop...");
    reservoir.Simulate2Phase(resultTime, wells);
    Console.WriteLine($"Simulation completed successfully without crashing in {toc():F2} seconds!");


    reservoir.ExportParaView("C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\RunTest1");
    reservoir.ExportWells("C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\RunTest1");
}




{
    // 2D - 2Phase
    folderpath = "C:\\Users\\lateef.a.kareem\\Documents\\GitHub\\ReservoirSimulation\\";
    int Nx = 25, Ny = 7, Lx = Nx - 1, Ly = Ny - 1, M = 2*Nx*Ny, Nwells = 2, WellIndex;
    double Pinit = 3000, Sinit = 0.2, Sw_r = 0.10, So_r = 0.15,
           μw0 = 5.005, μo0 = 2, kro0 = 1.0, krw0 = 0.30, Pe = 2,
           co = 2e-5, cw = 4e-6, cr = 1e-5, bo = 2e-5, bw = 4e-10,
           Bw0 = 1.005, Bo0 = 1.4, no = 2.5, nw = 3;
    var Producer = (MinPressure: 1500.0, ProdRate: 0.0, OilRate: 0.0, WaterRate: 0.0, I: 0, J: 3);
    var Injector = (MaxPressure: 4500.0, InjRate: 0.0, OilRate: 0.0, WaterRate: 0.0, I: 24, J: 0);

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
    for (int i = 0; i < Ny-3; i++) K[12 + i*25] *= 0.00001;

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

    double EndTime = 5000, delt = EndTime/300;
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
        AnimationMaker(Animfun, $"2D_WaterFlooding_{rate}_with lower barrier.gif", 30, 300);
        Console.WriteLine(DateTime.Now);
        CloseFig();
        Console.WriteLine("=======================================================");
    }
}