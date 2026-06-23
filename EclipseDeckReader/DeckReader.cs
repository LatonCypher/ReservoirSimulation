using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Text.RegularExpressions;

namespace EclipseDeckReader;

public sealed class DeckReader
{
    public EclipseDataDeck LoadDeck(string masterFilePath)
    {
        if (!File.Exists(masterFilePath))
            throw new FileNotFoundException($"Master ECLIPSE deck path not found: {masterFilePath}");

        // 1. Resolve recursive INCLUDE streams into one in-memory block
        string combinedText = ReadAndResolveIncludes(masterFilePath);

        // 2. Isolate keyword sections using the updated robust lookahead parser
        var rawSections = ExtractRawKeywordBlocks(combinedText);

        // 3. Construct empty container and map dimensions
        EclipseDataDeck deck = new();
        ParseGridDimensions(rawSections, deck);

        // 4. Populate clean primitive arrays
        if (rawSections.TryGetValue("COORD", out var coordVal)) deck.Coord = UnpackDoubleArray(coordVal);
        if (rawSections.TryGetValue("ZCORN", out var zcornVal)) deck.Zcorn = UnpackDoubleArray(zcornVal);
        if (rawSections.TryGetValue("ACTNUM", out var actnumVal)) deck.Actnum = UnpackIntArray(actnumVal);
        if (rawSections.TryGetValue("PORO", out var poroVal)) deck.Porosity = UnpackDoubleArray(poroVal);
        if (rawSections.TryGetValue("NTG", out var ntgVal)) deck.Ntg = UnpackDoubleArray(ntgVal);

        // 5. Parse out structural FAULTS planes (Norne Style)
        if (rawSections.TryGetValue("FAULTS", out var faultsVal))
        {
            ParseFaultRecords(faultsVal, deck);
        }

        // 6. Overwrite with explicit MULTFLT trans multipliers if present
        if (rawSections.TryGetValue("MULTFLT", out var multfltVal))
        {
            ApplyFaultMultipliers(multfltVal, deck);
        }

        return deck;
    }

    private string ReadAndResolveIncludes(string filePath)
    {
        string baseDir = Path.GetDirectoryName(filePath) ?? Directory.GetCurrentDirectory();
        string fileText = File.ReadAllText(filePath);

        string includePattern = @"(?i)\bINCLUDE\b\s+['""]([^'""]+)['""]\s*/";

        return Regex.Replace(fileText, includePattern, match =>
        {
            string relativePath = match.Groups[1].Value.Trim()
                .Replace('\\', Path.DirectorySeparatorChar)
                .Replace('/', Path.DirectorySeparatorChar);

            string fullIncludePath = Path.Combine(baseDir, relativePath);

            if (File.Exists(fullIncludePath))
            {
                return ReadAndResolveIncludes(fullIncludePath);
            }
            throw new FileNotFoundException($"ECLIPSE Include statement broken. File missing at: {fullIncludePath}");
        }, RegexOptions.Singleline);
    }

    private Dictionary<string, string> ExtractRawKeywordBlocks(string fullText)
    {
        var sections = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);

        var targetKeywords = new HashSet<string>(StringComparer.OrdinalIgnoreCase)
    {
        "SPECGRID", "DIMENS", "COORD", "ZCORN", "ACTNUM", "PORO", "NTG", "FAULTS", "MULTFLT"
    };

        string[] lines = fullText.Split(new[] { "\r\n", "\r", "\n" }, StringSplitOptions.None);

        string? currentKeyword = null;
        StringBuilder currentSectionData = new();

        foreach (string rawLine in lines)
        {
            string line = rawLine.Trim();
            if (string.IsNullOrEmpty(line)) continue;

            // Strip comments instantly so they never interfere with keyword boundaries
            int commentIdx = line.IndexOf("--");
            if (commentIdx == 0) continue; // Full comment line
            if (commentIdx > 0) line = line.Substring(0, commentIdx).Trim();

            // Check if the clean line is a single standalone keyword word
            if (targetKeywords.Contains(line))
            {
                // If we were already tracking a keyword, save its accumulated text block
                if (currentKeyword != null)
                {
                    sections[currentKeyword] = currentSectionData.ToString().Trim();
                }

                // Pivot to the new keyword section and flush the buffer
                currentKeyword = line.ToUpper();
                currentSectionData.Clear();
                continue;
            }

            // If we are currently inside a targeted section, accumulate the data line
            if (currentKeyword != null)
            {
                currentSectionData.AppendLine(line);
            }
        }

        // Capture the final section remaining at EOF
        if (currentKeyword != null && currentSectionData.Length > 0)
        {
            sections[currentKeyword] = currentSectionData.ToString().Trim();
        }

        return sections;
    }
    private void ParseGridDimensions(Dictionary<string, string> sections, EclipseDataDeck deck)
    {
        string? dimContent = null;
        if (sections.TryGetValue("SPECGRID", out var specgrid)) dimContent = specgrid;
        else if (sections.TryGetValue("DIMENS", out var dimens)) dimContent = dimens;

        if (dimContent == null)
            throw new InvalidDataException("Incomplete specification deck: Neither SPECGRID nor DIMENS parameters found.");

        string[] tokens = EclipseTokenizer.TokenizeSection(dimContent);
        deck.Nx = int.Parse(tokens[0]);
        deck.Ny = int.Parse(tokens[1]);
        deck.Nz = int.Parse(tokens[2]);
    }

    private double[] UnpackDoubleArray(string content)
    {
        string[] tokens = EclipseTokenizer.TokenizeSection(content);
        List<double> values = new(tokens.Length);

        foreach (string token in tokens)
        {
            if (token.Contains('*'))
            {
                var parts = token.Split('*');
                int count = int.Parse(parts[0]);
                double val = double.Parse(parts[1]);

                if (count > 1000) values.Capacity += count;
                for (int i = 0; i < count; i++) values.Add(val);
            }
            else
            {
                values.Add(double.Parse(token));
            }
        }
        return values.ToArray();
    }

    private int[] UnpackIntArray(string content)
    {
        string[] tokens = EclipseTokenizer.TokenizeSection(content);
        List<int> values = new(tokens.Length);

        foreach (string token in tokens)
        {
            if (token.Contains('*'))
            {
                var parts = token.Split('*');
                int count = int.Parse(parts[0]);
                int val = int.Parse(parts[1]);
                for (int i = 0; i < count; i++) values.Add(val);
            }
            else
            {
                values.Add(int.Parse(token));
            }
        }
        return values.ToArray();
    }

    private void ParseFaultRecords(string content, EclipseDataDeck deck)
    {
        string[] tokens = EclipseTokenizer.TokenizeSection(content);
        int idx = 0;

        // Reads the standard 8-column structural plane layout
        while (idx < tokens.Length)
        {
            if (tokens[idx] == "/") { idx++; continue; }
            if (idx + 8 > tokens.Length) break;

            try
            {
                FaultRecord fault = new()
                {
                    FaultName = tokens[idx++].Replace("'", "").Replace("\"", ""),
                    IMin = int.Parse(tokens[idx++]),
                    IMax = int.Parse(tokens[idx++]),
                    JMin = int.Parse(tokens[idx++]),
                    JMax = int.Parse(tokens[idx++]),
                    KMin = int.Parse(tokens[idx++]),
                    KMax = int.Parse(tokens[idx++]),
                    Direction = tokens[idx++].ToUpper(),
                    TransmissibilityMultiplier = 1.0 // Default baseline until MULTFLT modifies it
                };
                deck.Faults.Add(fault);
            }
            catch
            {
                idx++; // Advance index window safely if minor layout variations pop up
            }
        }
    }

    private void ApplyFaultMultipliers(string content, EclipseDataDeck deck)
    {
        string[] tokens = EclipseTokenizer.TokenizeSection(content);
        int idx = 0;

        // Process potential MULTFLT overrides and match them by name to existing structural records
        while (idx < tokens.Length)
        {
            if (tokens[idx] == "/") { idx++; continue; }
            if (idx + 9 > tokens.Length) break;

            try
            {
                string targetName = tokens[idx++].Replace("'", "").Replace("\"", "");
                string dir = tokens[idx++];
                int iMin = int.Parse(tokens[idx++]);
                int iMax = int.Parse(tokens[idx++]);
                int jMin = int.Parse(tokens[idx++]);
                int jMax = int.Parse(tokens[idx++]);
                int kMin = int.Parse(tokens[idx++]);
                int kMax = int.Parse(tokens[idx++]);
                double multiplier = double.Parse(tokens[idx++]);

                // Query and find the matching fault tracking profile to assign the explicit historical multiplier
                var matchedFault = deck.Faults.Find(f =>
                    string.Equals(f.FaultName, targetName, StringComparison.OrdinalIgnoreCase) &&
                    f.IMin == iMin && f.JMin == jMin && f.KMin == kMin);

                if (matchedFault != null)
                {
                    // Update property value context on the reference entity
                    // (Requires setting the property to private set or modifying it if needed, or re-adding)
                    // For direct data ingestion transfer, we modify the existing tracked list element
                    typeof(FaultRecord).GetProperty(nameof(FaultRecord.TransmissibilityMultiplier))?
                                      .SetValue(matchedFault, multiplier);
                }
            }
            catch
            {
                idx++;
            }
        }
    }
}