namespace ReservoirSimulator
{
    internal class ScratchPad
    {
        public ADiff Tr = new(), Er = new(), Oil_Rate = new(), Oil_Vol = new(), Water_Rate = new(), Water_Vol = new();

        public ADiff Po_m = new(), So_m = new(), Go_m = new(), Bo_m = new(), Uo_m = new(), Kro_m = new(), PmGZo_m = new(),
                     Pw_m = new(), Sw_m = new(), Gw_m = new(), Bw_m = new(), Uw_m = new(), Krw_m = new(), PmGZw_m = new(),
                     SoPo_m = new(), SwPw_m = new(), meanP = new();

        public ADiff Po_n = new(), So_n = new(), Go_n = new(), Bo_n = new(), Uo_n = new(), Kro_n = new(), PmGZo_n = new(),
                     Pw_n = new(), Sw_n = new(), Gw_n = new(), Bw_n = new(), Uw_n = new(), Krw_n = new(), PmGZw_n = new();

        public void ClearMVariables()
        {
            Po_m.Clear(); So_m.Clear(); Go_m.Clear(); Bo_m.Clear(); Uo_m.Clear(); Kro_m.Clear(); PmGZo_m.Clear();
            Pw_m.Clear(); Sw_m.Clear(); Gw_m.Clear(); Bw_m.Clear(); Uw_m.Clear(); Krw_m.Clear(); PmGZw_m.Clear();
            Oil_Rate.Clear(); Oil_Vol.Clear(); Water_Rate.Clear(); Water_Vol.Clear(); Tr.Clear();
            SoPo_m.Clear(); SwPw_m.Clear(); meanP.Clear(); Er.Clear();
        }
        public void ClearNVariables()
        {
            Po_n.Clear(); So_n.Clear(); Go_n.Clear(); Bo_n.Clear(); Uo_n.Clear(); Kro_n.Clear(); PmGZo_n.Clear();
            Pw_n.Clear(); Sw_n.Clear(); Gw_n.Clear(); Bw_n.Clear(); Uw_n.Clear(); Krw_n.Clear(); PmGZw_n.Clear();
            Oil_Rate.Clear(); Oil_Vol.Clear(); Water_Rate.Clear(); Water_Vol.Clear(); Tr.Clear(); 
        }
    }
}
