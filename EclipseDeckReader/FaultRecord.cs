namespace EclipseDeckReader;

/// <summary>
/// Represents a fault transmissibility multiplier record extracted from an ECLIPSE model.
/// </summary>
public sealed class FaultRecord
{
    public required string FaultName { get; init; }
    public required string Direction { get; init; } // I, J, or K vertical/horizontal alignment
    public required double TransmissibilityMultiplier { get; init; }

    // Bounding grid coordinate box where the fault is active
    public required int IMin { get; init; }
    public required int IMax { get; init; }
    public required int JMin { get; init; }
    public required int JMax { get; init; }
    public required int KMin { get; init; }
    public required int KMax { get; init; }
}