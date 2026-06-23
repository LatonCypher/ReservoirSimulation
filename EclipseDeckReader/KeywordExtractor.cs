using System.Text.RegularExpressions;

namespace EclipseDeckReader;

internal static partial class KeywordExtractor
{
    // High-performance .NET Source Generator Regex for matching ECLIPSE keyword blocks
    // Captures everything between the keyword boundaries up to a trailing slash '/'
    [GeneratedRegex(@"(?i)\b(SPECGRID|DIMENS|COORD|ZCORN|ACTNUM|PORO|NTG|MULTFLT)\b\s+(.*?)\s*/", RegexOptions.Singleline)]
    private static partial Regex EclipseKeywordRegex();

    public static Dictionary<string, string> ParseRawBlocks(string fullResolvedText)
    {
        var sections = new Dictionary<string, string>(StringComparer.OrdinalIgnoreCase);
        var matches = EclipseKeywordRegex().Matches(fullResolvedText);

        foreach (Match match in matches)
        {
            string keyword = match.Groups[1].Value.ToUpper();
            string payload = match.Groups[2].Value.Trim();

            // Handle cases where a file might overwrite or define a keyword sequentially
            sections[keyword] = payload;
        }

        return sections;
    }
}