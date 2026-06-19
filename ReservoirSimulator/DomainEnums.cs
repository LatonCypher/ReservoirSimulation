namespace ReservoirSimulator
{
    internal class Jacobian
    {
        public List<double> Value = [];
        public List<int> Index = [];
        public List<int> Start = [0];
        public Jacobian Duplicate() =>
            new() { Value = [.. Value], Index = [.. Index], Start = [.. Start] };
    }

    public enum ConstraintType
    {
        MaxPressure,
        MinPressure,
        FlowRate
    }

    public enum WellType
    {
        Injector,
        Producer
    }

    public enum Phase
    {
        Oil,
        Water
    }

    public enum Direction
    {
        X,
        Y,
        Z
    }

    public enum AquiferFlowDirection
    {
        Iplus,
        Iminus,
        Jplus,
        Jminus,
        Kplus,
        Kminus
    }
}
