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

    public enum WellStatus
    {
        OPEN,
        CLOSED,
    }

    public enum ConstraintType
    {
        MAXP,
        MINP,
        LRATE,
        ORATE,
        WRATE,
        GRATE,
        RATE
    }

    public enum WellType
    {
        INJECTOR,
        PRODUCER
    }

    public enum Phase
    {
        OIL,
        WATER,
        GAS,
        DISGAS,
        VAPOIL
    }

    public enum Direction
    {
        X,
        Y,
        Z
    }

    public enum FlowDirection
    {
        KM,
        JM,
        IM,
        IP,
        JP,
        KP
    }
}
