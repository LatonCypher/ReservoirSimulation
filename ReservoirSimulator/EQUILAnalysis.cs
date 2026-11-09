using SepalSolver;
using static SepalSolver.Math;

namespace ReservoirSimulator
{
    public readonly struct EquilResult(double zDatum, double pDatum, double zWoc, double pcWoc, 
         double zGoc, double pcGoc, double poWoc, double pwWoc, double poGoc, double pgGoc) 
    {
        public double ZDatum { get; } = zDatum;
        public double PDatum { get; } = pDatum;
        public double ZWOC { get; } = zWoc;
        public double PcWOC { get; } = pcWoc;
        public double ZGOC { get; } = zGoc;
        public double PcGOC { get; } = pcGoc;
        public double PoWOC { get; } = poWoc;
        public double PwWOC { get; } = pwWoc;
        public double PoGOC { get; } = poGoc;
        public double PgGOC { get; } = pgGoc;
    }

    

    public class EQUILAnalysis
    {
        public double ZDatum { get; set; }
        public double PDatum { get; set; }
        public double zWOC { get; set; }
        public double pcWOC { get; set; }
        public double zGOC { get; set; }
        public double pcGOC { get; set; }
        public double PoWOC { get; set; }
        public double PoGOC { get; set; }
        public double PwWOC { get; set; }
        public double PgGOC { get; set; }

        public EquilResult Result { get; private set; }

        public EQUILAnalysis(List<double> equil)
        {
            ZDatum = equil[0]; 
            PDatum = equil[1]; 
            zWOC = equil[2];
            pcWOC = equil.Count >= 4 ? equil[3] : 0;
            zGOC = equil.Count >= 5 ? equil[4] : 0;
            pcGOC = equil.Count >= 6 ? equil[5] : 0;
        }

        public void CalculateEquilProperties<T>(Func<ADiff<T>, ADiff<T>> γo, Func<ADiff<T>, ADiff<T>> γw, Func<ADiff<T>, ADiff<T>> γg) where T : struct, IVectorND<T>
        {
            if (ZDatum <= zGOC)
            {
                double ResidualGOC(double p) =>
                    p - (PDatum + 0.5 * (γg(PDatum) + γg(p)) * (zGOC - ZDatum)).Value;
                PgGOC = Fsolve(ResidualGOC, PDatum);  PoGOC = PgGOC - pcGOC;

                double ResidualWOC(double p) =>
                    p - (PoGOC + 0.5 * (γo(PoGOC) + γo(p)) * (zWOC - zGOC)).Value;
                PoWOC = Fsolve(ResidualWOC, PoGOC);  PwWOC = PoWOC - pcWOC;
            }
            else if (ZDatum > zGOC && ZDatum < zWOC)
            {
                double ResidualWOC(double p) =>
                    p - (PDatum + 0.5 * (γo(PDatum) + γo(p)) * (zWOC - ZDatum)).Value;
                PoWOC = Fsolve(ResidualWOC, PDatum); PwWOC = PoWOC - pcWOC;

                double ResidualGOC(double p) =>
                    p - (PDatum + 0.5 * (γo(PDatum) + γo(p)) * (zGOC - ZDatum)).Value;
                PoGOC = Fsolve(ResidualGOC, PDatum); PgGOC = PoGOC + pcGOC;
            }
            else
            {
                double ResidualWOC(double p) =>
                    p - (PDatum + 0.5 * (γw(PDatum) + γw(p)) * (zWOC - ZDatum)).Value;
                PwWOC = Fsolve(ResidualWOC, PDatum); PoWOC = PwWOC + pcWOC;

                double ResidualGOC(double p) =>
                    p - (PoWOC + 0.5 * (γo(PoWOC) + γo(p)) * (zGOC - zWOC)).Value;
                PoGOC = Fsolve(ResidualGOC, PoWOC); PgGOC = PoGOC + pcGOC;
            }
            Result = new(ZDatum, PDatum, zWOC, pcWOC, zGOC, pcGOC, PoWOC, PwWOC, PoGOC, PgGOC);
        }
    }
}
