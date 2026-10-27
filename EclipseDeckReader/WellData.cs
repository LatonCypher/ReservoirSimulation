using System;
using System.Collections.Generic;

namespace EclipseDeckReader
{
    public sealed class WellSpecRecord
    {
        public string WellName { get; set; } = string.Empty;
        public string GroupName { get; set; } = string.Empty;
        public int HeadI { get; set; }
        public int HeadJ { get; set; }
        public double ReferenceDepth { get; set; }
        public string PreferredPhase { get; set; } = string.Empty;
    }

    public sealed class WellCompletionRecord
    {
        // Timestamp when this specific mechanical completion configuration/status becomes active
        public DateTime EffectiveDate { get; set; }
        public int I { get; set; }
        public int J { get; set; }
        public int KStart { get; set; }
        public int KEnd { get; set; }
        public string Status { get; set; } = "OPEN"; // Connection-level switch ("OPEN" / "SHUT")
        public double WellBoreDiameter { get; set; }
        public double EffectiveSkin { get; set; }
    }

    public sealed class WellControlRecord
    {
        // Timestamp when these operational limits/controls become active
        public DateTime EffectiveDate { get; set; }
        public string WellType { get; set; } = "PRODUCER";       // "PRODUCER" or "INJECTOR"
        public string InjectionPhase { get; set; } = "NONE";      // "WATER", "GAS", or "NONE"
        public string WellStatus { get; set; } = "OPEN";          // Wellhead master switch ("OPEN" / "SHUT")
        public string ControlMode { get; set; } = string.Empty;   // ORAT, GRAT, LRAT, BHP, etc.
        public double TargetRate { get; set; }
        public double OilRateTarget { get; set; }
        public double WaterRateTarget { get; set; }
        public double GasRateTarget { get; set; }
        public double LiquidRateTarget { get; set; }
        public double MinBhpLimit { get; set; }
        public double MaxBhpLimit { get; set; }
    }

    public sealed class WellData
    {
        public WellSpecRecord Specification { get; set; } = new();
        public List<WellCompletionRecord> CompletionHistory { get; } = []; // Track sequential mechanical events
        public List<WellControlRecord> ControlHistory { get; } = [];       // Track sequential operational events
    }
}