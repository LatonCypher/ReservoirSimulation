using EclipseDeckReader;
using System.Diagnostics;
using SepalSolver;
using System.Text;


{
    // 1. Define the path to your master ECLIPSE file (.DATA or .GRDECL)
    // If testing with the Norne deck, point it to your local file path
    string masterDeckPath = @"C:\Users\lateef.a.kareem\Documents\GitHub\opm-data\norne\NORNE_ATW2013.DATA";

    Console.WriteLine("==================================================");
    Console.WriteLine("        LAUNCHING ECLIPSE INGESTION PIPELINE     ");
    Console.WriteLine("==================================================");

    if (!File.Exists(masterDeckPath))
    {
        Console.ForegroundColor = ConsoleColor.Red;
        Console.WriteLine($"[Error] Could not locate file: {masterDeckPath}");
        Console.ResetColor();
        return;
    }

    // 2. Initialize your reader and track performance
    Stopwatch timer = Stopwatch.StartNew();
    EclipseDataDeck parsedDeck = null;
    try
    {
        // 1. Initialize stateless reader tool
        var reader = new DeckReader();

        // 2. Load the input deck into memory
        parsedDeck = reader.LoadDeckFile(masterDeckPath);

        // 3. Print out the formatted analysis report
        DeckPrinter.PrintSummary(parsedDeck);
    }
    catch (Exception ex)
    {
        Console.ForegroundColor = ConsoleColor.Red;
        Console.WriteLine($"\n[Critical Ingestion Failure]: {ex.Message}");
        if (ex.InnerException != null)
            Console.WriteLine($" ↳ Detail: {ex.InnerException.Message}");
        Console.ResetColor();
    }

    Console.WriteLine("\nPipeline idle. Press any key to close context window...");

    OilWaterReservoir Reservoir = new(parsedDeck);

    Reservoir.Simulate2Phase([0, 3650]);
}