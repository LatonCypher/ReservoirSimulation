using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EclipseDeckReader
{
    public sealed class DeckReader
    {
        public EclipseDataDeck LoadDeckFile(string filepath)
        {
            var deck = new EclipseDataDeck();
            if (!File.Exists(filepath))
                throw new FileNotFoundException("Master deck file not found.", filepath);

            Console.WriteLine($"\n[Pipeline] MASTER DECK OPENED: {Path.GetFileName(filepath)}");
            string rawText = File.ReadAllText(filepath);
            List<Token> tokens = EclipseTokenizer.Tokenize(rawText.AsMemory());
            ParseTokenStream(tokens, deck, filepath);
            return deck;
        }

        private void ParseSubFile(string subFilePath, EclipseDataDeck deck)
        {
            Console.WriteLine($"   └─ [Pipeline] Entering Include File: {Path.GetFileName(subFilePath)}");
            string rawText = File.ReadAllText(subFilePath);
            List<Token> tokens = EclipseTokenizer.Tokenize(rawText.AsMemory());
            ParseTokenStream(tokens, deck, subFilePath);
            Console.WriteLine($"   └─ [Pipeline] Exiting Include File: {Path.GetFileName(subFilePath)}");
        }

        private void ParseTokenStream(List<Token> tokens, EclipseDataDeck deck, string currentFilePath)
        {
            int i = 0;
            while (i < tokens.Count)
            {
                if (tokens[i].Type == TokenType.Keyword)
                {
                    string keyword = tokens[i].Content.ToString().ToUpper();
                    i++;

                    while (i < tokens.Count &&
                          (tokens[i].Type == TokenType.NewLine ||
                           string.IsNullOrWhiteSpace(tokens[i].Content.ToString()) ||
                           tokens[i].Content.Span.StartsWith("--")))
                    {
                        i++;
                    }

                    switch (keyword)
                    {
                        case "DIMENS":
                        case "SPECGRID":
                            var dim = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            if (dim.Count >= 3)
                            {
                                deck.NX = dim[0]; deck.NY = dim[1]; deck.NZ = dim[2];
                                InitializeMultiplierArrays(deck);
                            }
                            break;

                        case "OIL":
                        case "WATER":
                        case "GAS":
                        case "DISGAS":
                        case "VAPOIL":
                            deck.PhaseFlags.Add(keyword);
                            break;

                        case "PORO":
                            deck.PorosityArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current PORO array status: {deck.PorosityArray.Count} entries.");
                            break;

                        case "PERMX":
                            deck.PermXArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current PERMX array status: {deck.PermXArray.Count} entries.");
                            break;

                        case "PERMY":
                            deck.PermYArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current PERMY array status: {deck.PermYArray.Count} entries.");
                            break;

                        case "PERMZ":
                            deck.PermZArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current PERMZ array status: {deck.PermZArray.Count} entries.");
                            break;

                        case "NTG":
                            deck.NtgArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current NTG array status: {deck.NtgArray.Count} entries.");
                            break;

                        case "MULTX":
                            deck.MultXArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current MULTX array status: {deck.MultXArray.Count} entries.");
                            break;

                        case "MULTY":
                            deck.MultYArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current MULTY array status: {deck.MultYArray.Count} entries.");
                            break;

                        case "MULTZ":
                            deck.MultZArray = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            Console.WriteLine($"      [Data] Current MULTZ array status: {deck.MultZArray.Count} entries.");
                            break;

                        case "ACTNUM":
                            var actData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            if (deck.ActnumArray.Count < actData.Count) { deck.ActnumArray = actData; }
                            Console.WriteLine($"      [Data] Current ACTNUM array status: {deck.ActnumArray.Count} entries.");
                            break;

                        case "COORD":
                            var coordData = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (deck.CoordArray.Count < coordData.Count) { deck.CoordArray = coordData; }
                            Console.WriteLine($"      [Data] Current COORD array status: {deck.CoordArray.Count} entries.");
                            break;

                        case "ZCORN":
                            var zcornData = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (deck.ZcornArray.Count < zcornData.Count) { deck.ZcornArray = zcornData; }
                            Console.WriteLine($"      [Data] Current ZCORN array status: {deck.ZcornArray.Count} entries.");
                            break;

                        case "SATNUM":
                            var satData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            // Shift 1-based ECLIPSE pointers to 0-based C# array indices
                            for (int idx = 0; idx < satData.Count; idx++) { satData[idx] -= 1; }
                            if (deck.SatnumArray.Count < satData.Count) { deck.SatnumArray = satData; }
                            Console.WriteLine($"      [Data] SATNUM: Loaded and converted {deck.SatnumArray.Count} entries to 0-based indexing.");
                            break;

                        case "PVTNUM":
                            var pvtData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            // Shift 1-based ECLIPSE pointers to 0-based C# array indices
                            for (int idx = 0; idx < pvtData.Count; idx++) { pvtData[idx] -= 1; }
                            if (deck.PvtnumArray.Count < pvtData.Count) { deck.PvtnumArray = pvtData; }
                            Console.WriteLine($"      [Data] PVTNUM: Loaded and converted {deck.PvtnumArray.Count} entries to 0-based indexing.");
                            break;

                        case "EQUALS":
                            ParseEqualsModifierBlock(tokens, ref i, deck);
                            break;

                        case "COPY":
                            ParseCopyModifierBlock(tokens, ref i, deck);
                            break;

                        case "MULTIPLY":
                            ParseMultiplyModifierBlock(tokens, ref i, deck);
                            break;

                        case "START": // Initial Day 0 baseline timestamp card
                            ParseStartBlock(tokens, ref i, deck);
                            break;

                        case "DATES": // Parsed timeline card block synchronization
                            ParseDatesBlock(tokens, ref i, deck);
                            break;

                        case "WELSPECS":
                            ParseWellSpecsBlock(tokens, ref i, deck);
                            break;

                        case "COMPDAT":
                            ParseCompdatBlock(tokens, ref i, deck);
                            break;

                        case "WCONPROD":
                            ParseWconprodBlock(tokens, ref i, deck);
                            break;

                        case "WCONINJE":
                            ParseWconinjeBlock(tokens, ref i, deck);
                            break;

                        case "PVDO":
                            deck.RawPvdoTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            break;
                        case "PVDG":
                            deck.RawPvdgTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            break;

                        case "PVTO":
                            deck.RawPvtoTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: false));
                            break;

                        case "PVTW":
                            deck.RawPvtwTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            Console.WriteLine($"      [Data] Current PVTW table status: {deck.RawPvtwTables.Count} region table(s) registered.");
                            break;

                        case "EQUIL": 
                            deck.RawEquilTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            Console.WriteLine($"      [Data] EQUIL: Successfully loaded {deck.RawEquilTables.Count} hydrostatic initialization record region(s).");
                            break;

                        case "EQUINUM":
                            var equilData = KeywordExtractor.ExtractIntArray(tokens, ref i);
                            // Shift 1-based ECLIPSE pointers to 0-based C# array indices
                            for (int idx = 0; idx < equilData.Count; idx++) { equilData[idx] -= 1; }
                            if (deck.EquinumArray.Count < equilData.Count) { deck.EquinumArray = equilData; }
                            Console.WriteLine($"      [Data] EQUINUM: Loaded and converted {deck.EquinumArray.Count} entries to 0-based indexing.");
                            break;

                        case "DENSITY": 
                            var densities = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (densities.Count >= 3)
                            {
                                deck.SurfaceOilDensity = densities[0];
                                deck.SurfaceWaterDensity = densities[1];
                                deck.SurfaceGasDensity = densities[2];
                                Console.WriteLine($"      [Data] DENSITY: Loaded Surface Baselines -> Oil: {deck.SurfaceOilDensity:F2} | Water: {deck.SurfaceWaterDensity:F2} | Gas: {deck.SurfaceGasDensity:F4}");
                            }
                            break;

                        case "ROCK": 
                            var rockProps = KeywordExtractor.ExtractDoubleArray(tokens, ref i);
                            if (rockProps.Count >= 2)
                            {
                                deck.RockReferencePressure = rockProps[0];
                                deck.RockCompressibility = rockProps[1];
                                Console.WriteLine($"      [Data] ROCK: Mechanical Properties -> Ref Pressure: {deck.RockReferencePressure:F2} | Compressibility: {deck.RockCompressibility:E4}");
                            }
                            break;

                        case "PVTG":
                            deck.RawPvtgTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: false));
                            break;

                        case "SWOF":
                            deck.RawSwofTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            break;
                        case "SGOF":
                            deck.RawSgofTables.AddRange(KeywordExtractor.ExtractMultiTableGroup(tokens, ref i, splitOnSingleSlash: true));
                            break;

                        case "FAULTS":
                            ParseFaultsSection(tokens, ref i, deck.Faults);
                            break;

                        case "MULTFLT":
                            ParseMultfltSection(tokens, ref i, deck.Faults);
                            break;

                        case "INCLUDE":
                            HandleIncludeKeyword(tokens, ref i, currentFilePath, deck);
                            break;

                        default:
                            while (i < tokens.Count && tokens[i].Type != TokenType.Slash && tokens[i].Type != TokenType.Keyword) i++;
                            if (i < tokens.Count && tokens[i].Type == TokenType.Slash) i++;
                            break;
                    }
                }
                else
                {
                    i++;
                }
            }
        }

        private static void InitializeMultiplierArrays(EclipseDataDeck deck)
        {
            int totalCells = deck.NX * deck.NY * deck.NZ;
            if (totalCells <= 0) return;

            InitializeUniformArray(deck.MultXArray, totalCells, 1.0);
            InitializeUniformArray(deck.MultYArray, totalCells, 1.0);
            InitializeUniformArray(deck.MultZArray, totalCells, 1.0);
            Console.WriteLine($"      [Grid Initialization] MULTX, MULTY, and MULTZ baseline profiles initialized with {totalCells} entries (Default factor = 1.0).");
        }

        private void ParseStartBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            ClearLeadingLayoutSpacers(tokens, ref index);
            if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; return; }

            if (int.TryParse(tokens[index].Content.ToString(), out int day))
            {
                string monthStr = tokens[index + 1].Content.ToString().Trim('\'').Trim('"').ToUpper();
                int.TryParse(tokens[index + 2].Content.ToString(), out int year);

                int month = monthStr switch
                {
                    "JAN" => 1,
                    "FEB" => 2,
                    "MAR" => 3,
                    "APR" => 4,
                    "MAY" => 5,
                    "JUN" => 6,
                    "JUL" => 7,
                    "AUG" => 8,
                    "SEP" => 9,
                    "OCT" => 10,
                    "NOV" => 11,
                    "DEC" => 12,
                    _ => 1
                };

                // Initialize both the baseline Start Date and the tracking Current Date
                deck.StartDate = new DateTime(year, month, day);
                deck.CurrentSimDate = deck.StartDate;

                Console.WriteLine($"      [Schedule Init] Master Simulation START date initialized to: {deck.StartDate:yyyy-MM-dd}");
            }

            index += 3;
            AdvanceAndClearTrailingJunk(tokens, ref index);
        }

        private void ParseDatesBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; break; }

                if (int.TryParse(tokens[index].Content.ToString(), out int day))
                {
                    string monthStr = tokens[index + 1].Content.ToString().Trim('\'').Trim('"').ToUpper();
                    int.TryParse(tokens[index + 2].Content.ToString(), out int year);

                    int month = monthStr switch
                    {
                        "JAN" => 1,
                        "FEB" => 2,
                        "MAR" => 3,
                        "APR" => 4,
                        "MAY" => 5,
                        "JUN" => 6,
                        "JUL" => 7,
                        "AUG" => 8,
                        "SEP" => 9,
                        "OCT" => 10,
                        "NOV" => 11,
                        "DEC" => 12,
                        _ => 1
                    };

                    deck.CurrentSimDate = new DateTime(year, month, day);
                    Console.WriteLine($"      [Schedule Sync] Simulator internal baseline timeline stepped to: {deck.CurrentSimDate:yyyy-MM-dd}");
                }
                index += 3;
                AdvanceAndClearTrailingJunk(tokens, ref index);
            }
        }

        private void ParseWellSpecsBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; break; }

                string wellName = tokens[index].Content.ToString().Trim('\'').Trim('"');
                string groupName = tokens[index + 1].Content.ToString().Trim('\'').Trim('"');
                int.TryParse(tokens[index + 2].Content.ToString(), out int headI);
                int.TryParse(tokens[index + 3].Content.ToString(), out int headJ);
                double.TryParse(tokens[index + 4].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double refDepth);
                string prefPhase = (index + 5 < tokens.Count) ? tokens[index + 5].Content.ToString().Trim('\'').Trim('"') : "OIL";

                index += 6;
                AdvanceAndClearTrailingJunk(tokens, ref index);

                if (!deck.Wells.TryGetValue(wellName, out var wellData))
                {
                    wellData = new WellData();
                    deck.Wells[wellName] = wellData;
                }

                wellData.Specification = new WellSpecRecord
                {
                    WellName = wellName,
                    GroupName = groupName,
                    HeadI = headI,
                    HeadJ = headJ,
                    ReferenceDepth = refDepth,
                    PreferredPhase = prefPhase
                };
            }
        }
        private void ParseCompdatBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; break; }

                string wellName = tokens[index].Content.ToString().Trim('\'').Trim('"');
                int.TryParse(tokens[index + 1].Content.ToString(), out int iRaw);
                int.TryParse(tokens[index + 2].Content.ToString(), out int jRaw);
                int.TryParse(tokens[index + 3].Content.ToString(), out int kStart);
                int.TryParse(tokens[index + 4].Content.ToString(), out int kEnd);
                string status = tokens[index + 5].Content.ToString().ToUpper();

                double.TryParse(tokens[index + 7].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double dia);
                double.TryParse(tokens[index + 9].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double skin);

                index += 10;
                AdvanceAndClearTrailingJunk(tokens, ref index);

                if (deck.Wells.TryGetValue(wellName, out var wellData))
                {
                    // Convert 1-based ECLIPSE block locations to 0-based internal matrix indexing
                    int iZeroBased = iRaw - 1;
                    int jZeroBased = jRaw - 1;

                    // Expand multi-layer completions cell-by-cell into the timeline record history array
                    for (int k = kStart; k <= kEnd; k++)
                    {
                        wellData.CompletionHistory.Add(new WellCompletionRecord
                        {
                            EffectiveDate = deck.CurrentSimDate, // Stamped timeline tracking
                            I = iZeroBased,
                            J = jZeroBased,
                            KStart = k,
                            KEnd = k,
                            Status = status,
                            WellBoreDiameter = dia,
                            EffectiveSkin = skin
                        });
                    }
                }
            }
        }
        private void ParseWconprodBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; break; }

                string wellName = tokens[index].Content.ToString().Trim('\'').Trim('"');
                string wellStatus = tokens[index + 1].Content.ToString().ToUpper(); // Capture master switch status (OPEN/SHUT)
                string mode = tokens[index + 2].Content.ToString().ToUpper();

                double.TryParse(tokens[index + 3].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double oRat);
                double.TryParse(tokens[index + 4].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double wRat);
                double.TryParse(tokens[index + 5].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double gRat);
                double.TryParse(tokens[index + 6].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double lRat);
                double.TryParse(tokens[index + 8].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double minBhp);

                index += 9;
                AdvanceAndClearTrailingJunk(tokens, ref index);

                if (deck.Wells.TryGetValue(wellName, out var wellData))
                {
                    wellData.ControlHistory.Add(new WellControlRecord
                    {
                        EffectiveDate = deck.CurrentSimDate, // Bound to active calendar step
                        WellType = "PRODUCER",
                        InjectionPhase = "NONE",
                        WellStatus = wellStatus,
                        ControlMode = mode,
                        OilRateTarget = oRat,
                        WaterRateTarget = wRat,
                        GasRateTarget = gRat,
                        LiquidRateTarget = lRat,
                        MinBhpLimit = minBhp
                    });
                }
            }
        }
        private void ParseWconinjeBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                ClearLeadingLayoutSpacers(tokens, ref index);
                if (index >= tokens.Count || tokens[index].Type == TokenType.Slash) { if (index < tokens.Count) index++; break; }

                string wellName = tokens[index].Content.ToString().Trim('\'').Trim('"');
                string phaseType = tokens[index + 1].Content.ToString().ToUpper(); // WATER, GAS, etc.
                string wellStatus = tokens[index + 2].Content.ToString().ToUpper(); // Capture master switch status (OPEN/SHUT)
                string mode = tokens[index + 3].Content.ToString().ToUpper();
                double.TryParse(tokens[index + 4].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double injRate);
                double.TryParse(tokens[index + 5].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double maxBhp);

                index += 6;
                AdvanceAndClearTrailingJunk(tokens, ref index);

                if (deck.Wells.TryGetValue(wellName, out var wellData))
                {
                    var control = new WellControlRecord
                    {
                        EffectiveDate = deck.CurrentSimDate, // Bound to active calendar step
                        WellType = "INJECTOR",
                        InjectionPhase = phaseType,
                        WellStatus = wellStatus,
                        ControlMode = mode,
                        MaxBhpLimit = maxBhp
                    };

                    if (phaseType == "WATER") control.WaterRateTarget = injRate;
                    else if (phaseType == "GAS") control.GasRateTarget = injRate;
                    control.TargetRate = injRate; // Sync shared shorthand indicator field

                    wellData.ControlHistory.Add(control);
                }
            }
        }
        private static void ClearLeadingLayoutSpacers(List<Token> tokens, ref int index)
        {
            while (index < tokens.Count &&
                  (tokens[index].Type == TokenType.NewLine ||
                   string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                   tokens[index].Content.Span.StartsWith("--")))
            {
                index++;
            }
        }

        private static void AdvanceAndClearTrailingJunk(List<Token> tokens, ref int index)
        {
            while (index < tokens.Count && tokens[index].Type != TokenType.Slash)
            {
                index++;
            }
            if (index < tokens.Count && tokens[index].Type == TokenType.Slash)
            {
                index++;
            }
            while (index < tokens.Count)
            {
                if (tokens[index].Type == TokenType.NewLine)
                {
                    index++;
                    break;
                }
                index++;
            }
        }

        private void ParseMultiplyModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                while (index < tokens.Count &&
                      (tokens[index].Type == TokenType.NewLine ||
                       string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                       tokens[index].Content.Span.StartsWith("--")))
                {
                    index++;
                }

                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash) { index++; break; }

                string targetProp = tokens[index].Content.ToString().ToUpper().Trim('\'').Trim('"');

                if (tokens[index].Type == TokenType.Keyword &&
                    targetProp != "PERMX" && targetProp != "PERMY" && targetProp != "PERMZ" && targetProp != "PORO" &&
                    targetProp != "NTG" && targetProp != "ACTNUM" && targetProp != "MULTX" && targetProp != "MULTY" && targetProp != "MULTZ" &&
                    targetProp != "DX" && targetProp != "DY" && targetProp != "DZ" && targetProp != "ZTOP")
                {
                    break;
                }

                if (index + 7 >= tokens.Count) break;

                double.TryParse(tokens[index + 1].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double multiplierValue);

                int.TryParse(tokens[index + 2].Content.ToString(), out int iStart);
                int.TryParse(tokens[index + 3].Content.ToString(), out int iEnd);
                int.TryParse(tokens[index + 4].Content.ToString(), out int jStart);
                int.TryParse(tokens[index + 5].Content.ToString(), out int jEnd);
                int.TryParse(tokens[index + 6].Content.ToString(), out int kStart);
                int.TryParse(tokens[index + 7].Content.ToString(), out int kEnd);

                index += 8;

                while (index < tokens.Count && tokens[index].Type != TokenType.Slash) index++;
                if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;

                while (index < tokens.Count)
                {
                    if (tokens[index].Type == TokenType.NewLine) { index++; break; }
                    index++;
                }

                int totalCells = deck.NX * deck.NY * deck.NZ;
                if (targetProp == "ACTNUM" && deck.ActnumArray.Count == totalCells)
                {
                    int modificationCount = 0;
                    for (int k = kStart; k <= kEnd; k++)
                    {
                        for (int j = jStart; j <= jEnd; j++)
                        {
                            for (int i = iStart; i <= iEnd; i++)
                            {
                                int flatIndex = (k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1);
                                if (flatIndex >= 0 && flatIndex < totalCells)
                                {
                                    deck.ActnumArray[flatIndex] = (int)(deck.ActnumArray[flatIndex] * multiplierValue);
                                    modificationCount++;
                                }
                            }
                        }
                    }
                    Console.WriteLine($"      [ARRAY MODIFICATION - MULTIPLY] Modified '{targetProp}': Multiplied {modificationCount} cells by {multiplierValue} [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
                    continue;
                }

                List<double> targetArray = targetProp switch
                {
                    "PERMX" => deck.PermXArray,
                    "PERMY" => deck.PermYArray,
                    "PERMZ" => deck.PermZArray,
                    "PORO" => deck.PorosityArray,
                    "NTG" => deck.NtgArray,
                    "MULTX" => deck.MultXArray,
                    "MULTY" => deck.MultYArray,
                    "MULTZ" => deck.MultZArray,
                    "DX" => deck.DxArray,
                    "DY" => deck.DyArray,
                    "DZ" => deck.DzArray,
                    "ZTOP" => deck.ZtopArray,
                    _ => null
                };

                if (targetArray != null && targetArray.Count == totalCells)
                {
                    int modificationCount = 0;
                    for (int k = kStart; k <= kEnd; k++)
                    {
                        for (int j = jStart; j <= jEnd; j++)
                        {
                            for (int i = iStart; i <= iEnd; i++)
                            {
                                int flatIndex = (k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1);
                                if (flatIndex >= 0 && flatIndex < totalCells)
                                {
                                    targetArray[flatIndex] *= multiplierValue;
                                    modificationCount++;
                                }
                            }
                        }
                    }
                    Console.WriteLine($"      [ARRAY MODIFICATION - MULTIPLY] Modified '{targetProp}': Multiplied {modificationCount} cells by {multiplierValue} [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
                }
            }
        }

        private void ParseCopyModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count && tokens[index].Type != TokenType.Slash)
            {
                if (tokens[index].Type == TokenType.NewLine ||
                    string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                    tokens[index].Content.Span.StartsWith("--"))
                {
                    index++;
                    continue;
                }

                if (index + 1 >= tokens.Count) break;

                string sourceProp = tokens[index].Content.ToString().ToUpper();
                string targetProp = tokens[index + 1].Content.ToString().ToUpper();
                index += 2;

                while (index < tokens.Count && tokens[index].Type != TokenType.Slash && tokens[index].Type != TokenType.Keyword) index++;
                if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;

                if (sourceProp == "ACTNUM" && deck.ActnumArray.Count > 0)
                {
                    var clonedArray = new List<int>(deck.ActnumArray);
                    if (targetProp == "ACTNUM") deck.ActnumArray = clonedArray;
                    Console.WriteLine($"      [ARRAY MODIFICATION - COPY] Modified '{targetProp}': Overwrote array values by cloning all {clonedArray.Count} elements from '{sourceProp}'");
                    continue;
                }

                List<double> sourceArray = sourceProp switch
                {
                    "PERMX" => deck.PermXArray,
                    "PERMY" => deck.PermYArray,
                    "PERMZ" => deck.PermZArray,
                    "PORO" => deck.PorosityArray,
                    "NTG" => deck.NtgArray,
                    "MULTX" => deck.MultXArray,
                    "MULTY" => deck.MultYArray,
                    "MULTZ" => deck.MultZArray,
                    "DX" => deck.DxArray,
                    "DY" => deck.DyArray,
                    "DZ" => deck.DzArray,
                    "ZTOP" => deck.ZtopArray,
                    _ => null
                };

                if (sourceArray != null && sourceArray.Count > 0)
                {
                    var clonedArray = new List<double>(sourceArray);
                    switch (targetProp)
                    {
                        case "PERMX": deck.PermXArray = clonedArray; break;
                        case "PERMY": deck.PermYArray = clonedArray; break;
                        case "PERMZ": deck.PermZArray = clonedArray; break;
                        case "PORO": deck.PorosityArray = clonedArray; break;
                        case "NTG": deck.NtgArray = clonedArray; break;
                        case "MULTX": deck.MultXArray = clonedArray; break;
                        case "MULTY": deck.MultYArray = clonedArray; break;
                        case "MULTZ": deck.MultZArray = clonedArray; break;
                        case "DX": deck.DxArray = clonedArray; break;
                        case "DY": deck.DyArray = clonedArray; break;
                        case "DZ": deck.DzArray = clonedArray; break;
                        case "ZTOP": deck.ZtopArray = clonedArray; break;
                    }
                    Console.WriteLine($"      [ARRAY MODIFICATION - COPY] Modified '{targetProp}': Overwrote array values by cloning all {clonedArray.Count} elements from '{sourceProp}'");
                }
            }
            if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;
        }

        private void ParseEqualsModifierBlock(List<Token> tokens, ref int index, EclipseDataDeck deck)
        {
            while (index < tokens.Count)
            {
                while (index < tokens.Count &&
                      (tokens[index].Type == TokenType.NewLine ||
                       string.IsNullOrWhiteSpace(tokens[index].Content.ToString()) ||
                       tokens[index].Content.Span.StartsWith("--")))
                {
                    index++;
                }

                if (index >= tokens.Count) break;
                if (tokens[index].Type == TokenType.Slash) { index++; break; }

                string targetProp = tokens[index].Content.ToString().ToUpper().Trim('\'').Trim('"');

                if (tokens[index].Type == TokenType.Keyword &&
                    targetProp != "PERMX" && targetProp != "PERMY" && targetProp != "PERMZ" && targetProp != "PORO" &&
                    targetProp != "NTG" && targetProp != "ACTNUM" && targetProp != "MULTX" && targetProp != "MULTY" && targetProp != "MULTZ" &&
                    targetProp != "DX" && targetProp != "DY" && targetProp != "DZ" && targetProp != "ZTOP")
                {
                    break;
                }
                index++;

                if (index >= tokens.Count) break;

                if (!double.TryParse(tokens[index].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double value))
                {
                    value = 0.0;
                }
                index++;

                int iStart = 1, iEnd = deck.NX;
                int jStart = 1, jEnd = deck.NY;
                int kStart = 1, kEnd = deck.NZ;

                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int iS)) { iStart = iS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int iE)) { iEnd = iE; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int jS)) { jStart = jS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int jE)) { jEnd = jE; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int kS)) { kStart = kS; index++; }
                if (index < tokens.Count && int.TryParse(tokens[index].Content.ToString(), out int kE)) { kEnd = kE; index++; }

                while (index < tokens.Count && tokens[index].Type != TokenType.Slash) index++;
                if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;

                while (index < tokens.Count)
                {
                    if (tokens[index].Type == TokenType.NewLine) { index++; break; }
                    index++;
                }

                int totalCells = deck.NX * deck.NY * deck.NZ;
                if (totalCells <= 0) continue;

                if (targetProp == "MULTX" || targetProp == "MULTY" || targetProp == "MULTZ" ||
                    targetProp == "PERMX" || targetProp == "PERMY" || targetProp == "PERMZ" ||
                    targetProp == "PORO" || targetProp == "NTG" || targetProp == "DX" ||
                    targetProp == "DY" || targetProp == "DZ" || targetProp == "ZTOP")
                {
                    List<double> targetArray = targetProp switch
                    {
                        "PERMX" => deck.PermXArray,
                        "PERMY" => deck.PermYArray,
                        "PERMZ" => deck.PermZArray,
                        "PORO" => deck.PorosityArray,
                        "NTG" => deck.NtgArray,
                        "MULTX" => deck.MultXArray,
                        "MULTY" => deck.MultYArray,
                        "MULTZ" => deck.MultZArray,
                        "DX" => deck.DxArray,
                        "DY" => deck.DyArray,
                        "DZ" => deck.DzArray,
                        "ZTOP" => deck.ZtopArray,
                        _ => null
                    };

                    if (targetArray != null)
                    {
                        // Ensure it's correctly allocated if it was completely empty prior to this EQUALS run
                        if (targetArray.Count == 0) InitializeUniformArray(targetArray, totalCells, 1.0);

                        int modificationCount = 0;
                        for (int k = kStart; k <= kEnd; k++)
                        {
                            for (int j = jStart; j <= jEnd; j++)
                            {
                                for (int i = iStart; i <= iEnd; i++)
                                {
                                    int flatIndex = (k - 1) * deck.NX * deck.NY + (j - 1) * deck.NX + (i - 1);
                                    if (flatIndex >= 0 && flatIndex < totalCells)
                                    {
                                        targetArray[flatIndex] = value;
                                        modificationCount++;
                                    }
                                }
                            }
                        }
                        Console.WriteLine($"      [ARRAY MODIFICATION - EQUALS] Modified '{targetProp}': Assigned value {value} across {modificationCount} cells [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
                    }
                }
                else
                {
                    // The value read from the EQUALS card
                    int adjustedValue = (int)value;

                    switch (targetProp)
                    {
                        case "ACTNUM":
                            if (deck.ActnumArray.Count == 0) InitializeUniformIntArray(deck.ActnumArray, totalCells, 1);
                            ApplyIntEqualsOverride(deck.ActnumArray, adjustedValue, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "SATNUM":
                            if (deck.SatnumArray.Count == 0) InitializeUniformIntArray(deck.SatnumArray, totalCells, 0); // Default to index 0
                            // Convert user input value (1-based) to 0-based internal storage
                            ApplyIntEqualsOverride(deck.SatnumArray, adjustedValue - 1, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "PVTNUM":
                            if (deck.PvtnumArray.Count == 0) InitializeUniformIntArray(deck.PvtnumArray, totalCells, 0); // Default to index 0
                            // Convert user input value (1-based) to 0-based internal storage
                            ApplyIntEqualsOverride(deck.PvtnumArray, adjustedValue - 1, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;

                        case "EQUINUM":
                            if (deck.EquinumArray.Count == 0) InitializeUniformIntArray(deck.EquinumArray, totalCells, 0); // Default to index 0
                            // Convert user input value (1-based) to 0-based internal storage
                            ApplyIntEqualsOverride(deck.EquinumArray, adjustedValue - 1, iStart, iEnd, jStart, jEnd, kStart, kEnd, deck.NX, deck.NY, totalCells, targetProp);
                            break;
                    }
                }
            }
        }

        private static void ApplyIntEqualsOverride(List<int> array, double value, int iStart, int iEnd, int jStart, int jEnd, int kStart, int kEnd, int nx, int ny, int totalCells, string propName)
        {
            int modificationCount = 0;
            for (int k = kStart; k <= kEnd; k++)
            {
                for (int j = jStart; j <= jEnd; j++)
                {
                    for (int i = iStart; i <= iEnd; i++)
                    {
                        int flatIndex = (k - 1) * nx * ny + (j - 1) * nx + (i - 1);
                        if (flatIndex >= 0 && flatIndex < totalCells)
                        {
                            array[flatIndex] = (int)value;
                            modificationCount++;
                        }
                    }
                }
            }
            Console.WriteLine($"      [ARRAY MODIFICATION - EQUALS] Modified '{propName}': Assigned region index {(int)value} across {modificationCount} cells [Box I:{iStart}-{iEnd}, J:{jStart}-{jEnd}, K:{kStart}-{kEnd}]");
        }

        private static void InitializeUniformArray(List<double> array, int size, double value)
        {
            array.Clear();
            array.Capacity = size;
            for (int c = 0; c < size; c++) array.Add(value);
        }

        private static void InitializeUniformIntArray(List<int> array, int size, int value)
        {
            array.Clear();
            array.Capacity = size;
            for (int c = 0; c < size; c++) array.Add(value); 
        }

        private void HandleIncludeKeyword(List<Token> tokens, ref int index, string baseFilePath, EclipseDataDeck deck)
        {
            if (index < tokens.Count && (tokens[index].Type == TokenType.StringLiteral || tokens[index].Type == TokenType.DataValue))
            {
                string relativePath = tokens[index].Content.ToString().Trim('\'').Trim('"');
                index++;

                string baseDir = Path.GetDirectoryName(baseFilePath) ?? "";
                string normalizedPath = relativePath.Replace('/', Path.DirectorySeparatorChar).Replace('\\', Path.DirectorySeparatorChar);
                string fullSubFilePath = Path.Combine(baseDir, normalizedPath);

                if (!File.Exists(fullSubFilePath))
                {
                    string fileName = Path.GetFileName(normalizedPath);
                    string subDirRel = Path.GetDirectoryName(normalizedPath) ?? "";
                    string targetDir = Path.GetFullPath(Path.Combine(baseDir, subDirRel));

                    if (Directory.Exists(targetDir))
                    {
                        string matchedFile = Directory.GetFiles(targetDir)
                            .FirstOrDefault(f => Path.GetFileName(f).Equals(fileName, StringComparison.OrdinalIgnoreCase));

                        if (matchedFile != null) fullSubFilePath = matchedFile;
                    }
                }

                if (File.Exists(fullSubFilePath))
                {
                    ParseSubFile(fullSubFilePath, deck);
                }
                else
                {
                    Console.BackgroundColor = ConsoleColor.DarkRed;
                    Console.ForegroundColor = ConsoleColor.White;
                    Console.WriteLine($"\n[CRITICAL WARNING] INCLUDE FILE NOT FOUND: \"{relativePath}\"");
                    Console.ResetColor();
                }
            }
            if (index < tokens.Count && tokens[index].Type == TokenType.Slash) index++;
        }

        private static void ParseFaultsSection(List<Token> tokens, ref int index, Dictionary<string, FaultRecord> faults)
        {
            int linesRead = 0;
            while (index < tokens.Count && tokens[index].Type != TokenType.Keyword)
            {
                if (tokens[index].Type == TokenType.Slash || tokens[index].Type == TokenType.NewLine) { index++; continue; }
                if (index + 7 >= tokens.Count) break;

                string name = tokens[index].Content.ToString().Trim('\'').Trim('"');

                int.TryParse(tokens[index + 1].Content.ToString(), out int i1);
                int.TryParse(tokens[index + 2].Content.ToString(), out int i2);
                int.TryParse(tokens[index + 3].Content.ToString(), out int j1);
                int.TryParse(tokens[index + 4].Content.ToString(), out int j2);
                int.TryParse(tokens[index + 5].Content.ToString(), out int k1);
                int.TryParse(tokens[index + 6].Content.ToString(), out int k2);

                string rawFace = tokens[index + 7].Content.ToString().Trim('\'').Trim('"');
                char face = rawFace.Length > 0 ? rawFace[0] : 'X';

                faults[name] = new FaultRecord { FaultName = name, IStart = i1, IEnd = i2, JStart = j1, JEnd = j2, KStart = k1, KEnd = k2, Face = face };
                linesRead++;
                index += 8;
            }
            Console.WriteLine($"      [Data] FAULTS section processed: {linesRead} structural boundary records mapped.");
        }

        private static void ParseMultfltSection(List<Token> tokens, ref int index, Dictionary<string, FaultRecord> faults)
        {
            int linesUpdated = 0;
            while (index < tokens.Count && tokens[index].Type != TokenType.Keyword)
            {
                if (tokens[index].Type == TokenType.Slash || tokens[index].Type == TokenType.NewLine) { index++; continue; }
                if (index + 1 >= tokens.Count) break;

                string faultName = tokens[index].Content.ToString().Trim('\'').Trim('"');

                if (double.TryParse(tokens[index + 1].Content.ToString(), System.Globalization.CultureInfo.InvariantCulture, out double multiplier))
                {
                    if (faults.TryGetValue(faultName, out var record))
                    {
                        record.TransmissibilityMultiplier = multiplier;
                        linesUpdated++;
                    }
                }
                index += 2;
            }
            Console.WriteLine($"      [Data] MULTFLT section processed: {linesUpdated} transmissibility overrides applied.");
        }
    }
}