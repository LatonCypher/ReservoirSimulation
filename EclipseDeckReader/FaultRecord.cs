namespace EclipseDeckReader;

/// <summary>
/// Represents a fault transmissibility multiplier record extracted from an ECLIPSE model.
/// </summary>
public sealed class FaultRecord
{
    public string FaultName { get; set; } = string.Empty;

    // Face boundaries: I1-I2, J1-J2, K1-K2
    public int IStart { get; set; }
    public int IEnd { get; set; }
    public int JStart { get; set; }
    public int JEnd { get; set; }
    public int KStart { get; set; }
    public int KEnd { get; set; }

    // Orientation Face: 'X', 'Y', or 'Z'
    public char Face { get; set; }

    // Transmissibility multiplier applied to connection faces (from MULTFLT)
    public double TransmissibilityMultiplier { get; set; } = 1.0;
}