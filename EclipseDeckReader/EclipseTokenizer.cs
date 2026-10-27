using System;
using System.Collections.Generic;

namespace EclipseDeckReader
{
    public enum TokenType
    {
        Keyword,
        DataValue,
        StringLiteral,
        Slash,
        NewLine,
        EndOfFile
    }

    public readonly struct Token
    {
        public TokenType Type { get; }
        public ReadOnlyMemory<char> Content { get; }

        public Token(TokenType type, ReadOnlyMemory<char> content)
        {
            Type = type;
            Content = content;
        }

        public override string ToString() => Content.ToString();
    }

    public static class EclipseTokenizer
    {
        public static List<Token> Tokenize(ReadOnlyMemory<char> fileContent)
        {
            var tokens = new List<Token>(32768);
            int index = 0;
            int length = fileContent.Length;
            ReadOnlySpan<char> span = fileContent.Span;

            while (index < length)
            {
                // 1. Skip Whitespace but PRESERVE Newlines (\r and \n)
                while (index < length && char.IsWhiteSpace(span[index]))
                {
                    if (span[index] == '\r' || span[index] == '\n')
                    {
                        break;
                    }
                    index++;
                }

                if (index >= length) break;

                // 2. Explicitly Handle Newline Tokens
                if (span[index] == '\r' || span[index] == '\n')
                {
                    tokens.Add(new Token(TokenType.NewLine, fileContent.Slice(index, 1)));
                    index++;
                    continue;
                }

                // 3. Handle Comments (-- line entries)
                if (index < length - 1 && span[index] == '-' && span[index + 1] == '-')
                {
                    while (index < length && span[index] != '\n' && span[index] != '\r')
                    {
                        index++;
                    }
                    continue;
                }

                // 4. Handle Keyword / Data Block Termination Marker
                if (span[index] == '/')
                {
                    tokens.Add(new Token(TokenType.Slash, fileContent.Slice(index, 1)));
                    index++;
                    continue;
                }

                // 5. Robust String/Quote Literal Extraction Block ('YES', 'W1', etc.)
                if (span[index] == '\'' || span[index] == '"')
                {
                    char quoteChar = span[index];
                    int stringStart = index;
                    index++;

                    while (index < length && span[index] != quoteChar)
                    {
                        index++;
                    }

                    if (index < length) index++;

                    ReadOnlyMemory<char> stringContent = fileContent.Slice(stringStart, index - stringStart);
                    tokens.Add(new Token(TokenType.StringLiteral, stringContent));
                    continue;
                }

                // 6. Standard Word / Number / Repetition string (e.g. 1500*0.33)
                int start = index;
                while (index < length && !char.IsWhiteSpace(span[index]) && span[index] != '/' && span[index] != '\'' && span[index] != '"')
                {
                    if (span[index] == '-' && index < length - 1 && span[index + 1] == '-')
                    {
                        break;
                    }
                    index++;
                }

                if (index == start)
                {
                    index++;
                    continue;
                }

                ReadOnlyMemory<char> tokenContent = fileContent.Slice(start, index - start);

                bool isKeyword = IsLikelyKeyword(tokenContent.Span);
                TokenType type = isKeyword ? TokenType.Keyword : TokenType.DataValue;

                tokens.Add(new Token(type, tokenContent));
            }

            return tokens;
        }

        private static bool IsLikelyKeyword(ReadOnlySpan<char> span)
        {
            if (span.Length > 8 || span.Length == 0) return false;

            foreach (char c in span)
            {
                if (!char.IsAsciiLetterUpper(c) && !char.IsDigit(c)) return false;
            }

            return char.IsAsciiLetterUpper(span[0]);
        }
    }
}