//using SkiaSharp;
using static System.Math;
using System.Collections.Generic;

namespace ReservoirSimulator
{
    public class ContourDataExtractor
    {
    //    public readonly struct ColorMapEntry
    //    {
    //        public SKColor Color { get; }
    //        public double ZValue { get; }

    //        public ColorMapEntry(byte r, byte g, byte b, double zValue)
    //        {
    //            Color = new SKColor(r, g, b);
    //            ZValue = zValue;
    //        }
    //    }

    //    public readonly struct MapBoundingBox
    //    {
    //        public int XStart { get; }
    //        public int XEnd { get; }
    //        public int YStart { get; }
    //        public int YEnd { get; }

    //        public MapBoundingBox(int xStart, int xEnd, int yStart, int yEnd)
    //        {
    //            if (xStart >= xEnd || yStart >= yEnd)
    //                throw new ArgumentException("Start coordinates must be strictly less than End coordinates.");

    //            XStart = xStart;
    //            XEnd = xEnd;
    //            YStart = yStart;
    //            YEnd = yEnd;
    //        }
    //    }

    //    private readonly List<ColorMapEntry> _colorMap;
    //    private readonly double _colorToleranceSq;

    //    public ContourDataExtractor(List<ColorMapEntry> colorMap, double colorTolerance = 15.0)
    //    {
    //        _colorMap = colorMap ?? throw new ArgumentNullException(nameof(colorMap));
    //        _colorToleranceSq = colorTolerance * colorTolerance;
    //    }

    //    public double[,] ExtractZValues(string imagePath, int Nx, int Ny, MapBoundingBox bbox)
    //    {
    //        using (SKBitmap bitmap = SKBitmap.Decode(imagePath))
    //        {
    //            return ExtractZValues(bitmap, Nx, Ny, bbox);
    //        }
    //    }

    //    public double[,] ExtractZValues(SKBitmap bitmap, int Nx, int Ny, MapBoundingBox bbox)
    //    {
    //        // Safety checks to ensure the provided bounding box fits inside the loaded image
    //        int xStart = Clamp(bbox.XStart, 0, bitmap.Width - 1);
    //        int xEnd = Clamp(bbox.XEnd, 0, bitmap.Width - 1);
    //        int yStart = Clamp(bbox.YStart, 0, bitmap.Height - 1);
    //        int yEnd = Clamp(bbox.YEnd, 0, bitmap.Height - 1);

    //        double[,] zValues = new double[Ny, Nx];

    //        // Calculate steps based strictly on the bounding box dimensions
    //        double activeWidth = xEnd - xStart;
    //        double activeHeight = yEnd - yStart;

    //        double xStep = activeWidth / Nx;
    //        double yStep = activeHeight / Ny;

    //        IntPtr pixelsPtr = bitmap.GetPixels();

    //        unsafe
    //        {
    //            uint* ptr = (uint*)pixelsPtr.ToPointer();
    //            int rowWords = bitmap.RowBytes / 4;

    //            for (int j = 0; j < Ny; j++)
    //            {
    //                // Correct for Y-axis flip: 
    //                // We map simulation row 0 to the BOTTOM of the bounding box (yEnd)
    //                double centerY = yEnd - ((j + 0.5) * yStep);
    //                int pixelY = Clamp((int)centerY, yStart, yEnd);

    //                uint* rowPtr = ptr + (pixelY * rowWords);

    //                for (int i = 0; i < Nx; i++)
    //                {
    //                    // Map simulation column to the bounding box X region
    //                    double centerX = xStart + ((i + 0.5) * xStep);
    //                    int pixelX = Clamp((int)centerX, xStart, xEnd);

    //                    uint rawColor = rowPtr[pixelX];
    //                    SKColor pixelColor = (SKColor)rawColor;

    //                    byte r = pixelColor.Red;
    //                    byte g = pixelColor.Green;
    //                    byte b = pixelColor.Blue;

    //                    // Instantly drop background white space and grid lines
    //                    if ((r >= 250 && g >= 250 && b >= 250) || (r <= 5 && g <= 5 && b <= 5))
    //                    {
    //                        zValues[j, i] = double.NaN;
    //                        continue;
    //                    }

    //                    zValues[j, i] = GetZValueFromColor(pixelColor);
    //                }
    //            }
    //        }

    //        return zValues;
    //    }

    //    private double GetZValueFromColor(SKColor targetColor)
    //    {
    //        double minDistanceSq = double.MaxValue;
    //        double bestZValue = double.NaN;

    //        for (int k = 0; k < _colorMap.Count; k++)
    //        {
    //            ColorMapEntry entry = _colorMap[k];
    //            double distSq = CalculateColorDistanceSquared(targetColor, entry.Color);

    //            if (distSq < minDistanceSq)
    //            {
    //                minDistanceSq = distSq;
    //                bestZValue = entry.ZValue;
    //            }
    //        }

    //        if (minDistanceSq > _colorToleranceSq)
    //        {
    //            return double.NaN;
    //        }

    //        return bestZValue;
    //    }

    //    private static double CalculateColorDistanceSquared(SKColor c1, SKColor c2)
    //    {
    //        double deltaR = c1.Red - c2.Red;
    //        double deltaG = c1.Green - c2.Green;
    //        double deltaB = c1.Blue - c2.Blue;
    //        return (0.299 * deltaR * deltaR) + (0.587 * deltaG * deltaG) + (0.114 * deltaB * deltaB);
    //    }
    }
}
