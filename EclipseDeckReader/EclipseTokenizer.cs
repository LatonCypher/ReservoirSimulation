using System.Text;
using System.Text.RegularExpressions;

namespace EclipseDeckReader;

internal static partial class EclipseTokenizer
{
    [GeneratedRegex(@"\s+", RegexOptions.Compiled)]
    private static partial Regex WhitespaceRegex();

    public static string[] TokenizeSection(string content)
    {
        string[] lines = content.Split(new[] { "\r\n", "\r", "\n" }, StringSplitOptions.None);
        StringBuilder cleanTextBuilder = new(content.Length); // Pre-size buffer to prevent structural re-allocations

        foreach (string line in lines)
        {
            string trimmed = line.Trim();
            if (string.IsNullOrEmpty(trimmed)) continue;

            // Strip inline comments
            int commentIdx = trimmed.IndexOf("--");
            if (commentIdx == 0) continue;
            if (commentIdx > 0) trimmed = trimmed.Substring(0, commentIdx).Trim();

            // Catch explicit data line breaks or stops
            if (trimmed == "/") break;

            cleanTextBuilder.Append(trimmed).Append(' ');
        }

        string[] rawTokens = WhitespaceRegex().Split(cleanTextBuilder.ToString().Trim());
        List<string> filteredTokens = new(rawTokens.Length);

        foreach (string token in rawTokens)
        {
            if (!string.IsNullOrWhiteSpace(token))
            {
                filteredTokens.Add(token);
            }
        }

        return filteredTokens.ToArray();
    }
}