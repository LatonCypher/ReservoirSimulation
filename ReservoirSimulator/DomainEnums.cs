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
        LiqRate,
        OilRate,
        WaterRate,
        GasRate
    }

    public enum WellType
    {
        Injector,
        Producer
    }

    public enum Phase
    {
        Oil,
        Water,
        Gas,
        DisGas,
        VapOil
    }

    public enum Direction
    {
        X,
        Y,
        Z
    }

    public enum FlowDirection
    {
        Iminus,
        Iplus,
        Jminus,
        Jplus,
        Kminus,
        Kplus
    }
}
