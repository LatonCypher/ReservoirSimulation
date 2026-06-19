using System;
using System.Collections.Generic;
using System.Text;

namespace ReservoirSimulator
{
    public class Aquifer
    {
        public int[] Iinterval = [];
        public int[] Jinterval = [];
        public int[] Kinterval = [];
        public AquiferFlowDirection FlowDirection;
        public double Connectivity_Efficiency, Pa;

        public bool IsthereAquiferFlow(int i, int j, int k)
        {
            return (Iinterval[0] <= i && i <= Iinterval[1]) &&
                   (Jinterval[0] <= j && j <= Jinterval[1]) &&
                   (Kinterval[0] <= k && k <= Kinterval[1]);
        }

        public Aquifer(int[] iinterval, int[] jinterval, int[] kinterval, AquiferFlowDirection flowDirection, double Pinit, double connectivity_Efficiency)
        {
            Iinterval = iinterval;
            Jinterval = jinterval;
            Kinterval = kinterval;
            FlowDirection = flowDirection;
            Connectivity_Efficiency = connectivity_Efficiency;
            Pa = Pinit;

            if (FlowDirection == AquiferFlowDirection.Iminus)
            {
                if (Iinterval[0] != Iinterval[1] && Iinterval[0] != 0)
                    throw new ArgumentException("Invalid flow direction for the given I interval.");
            }
            else if (FlowDirection == AquiferFlowDirection.Iplus)
            {
                if (Iinterval[0] != Iinterval[1] && Iinterval[1] >= 0)
                    throw new ArgumentException("Invalid flow direction for the given I interval.");
            }

            if (FlowDirection == AquiferFlowDirection.Jminus)
            {
                if (Jinterval[0] != Jinterval[1] && Jinterval[0] != 0)
                    throw new ArgumentException("Invalid flow direction for the given J interval.");
            }
            else if (FlowDirection == AquiferFlowDirection.Jplus)
            {
                if (Jinterval[0] != Jinterval[1] && Jinterval[1] >= 0)
                    throw new ArgumentException("Invalid flow direction for the given J interval.");
            }

            if(FlowDirection == AquiferFlowDirection.Kminus)
            {
                if (Kinterval[0] != Kinterval[1] && Kinterval[0] != 0)
                    throw new ArgumentException("Invalid flow direction for the given K interval.");
            }
            else if (FlowDirection == AquiferFlowDirection.Kplus)
            {
                if (Kinterval[0] != Kinterval[1] && Kinterval[1] >= 0)
                    throw new ArgumentException("Invalid flow direction for the given K interval.");
            }
        }
    }
}
