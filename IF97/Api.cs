using System;

namespace IF97
{
    public sealed class Api
    {
        //Table 1. Numerical values of the coefficients of the B23-equation, Eqs. (5) and (6), for
        //defining the boundary between regions 2 and 3
        static double[] Region23data = {
            0.34805185628969e3,
            -0.11671859879975e1,
            0.10192970039326e-2,
            0.57254459862746e3,
            0.13918839778870e2
        };

        static double Region23_T(double T)
        {
            double p_star = 1, T_star = 1, theta = T / T_star;
            double PI = Region23data[0] + Region23data[1] * theta + Region23data[2] * theta * theta;
            return PI * p_star;
        }
        static double Region23_p(double p)
        {
            double p_star = 1, T_star = 1, PI = p / p_star;
            double THETA = Region23data[3] + Math.Sqrt((PI - Region23data[4]) / Region23data[2]);
            return THETA * T_star;
        }
        static IF97REGIONS RegionDetermination_TP(double T, double p)
        {
            if (T > Constants.Text)
            {
                throw new ArgumentOutOfRangeException("Temperature out of range");
            }
            else if (T > Constants.Tmax && T <= Constants.Text)
            {
                if (p <= Constants.Pext)
                {
                    return IF97REGIONS.REGION_5;
                }
                else
                {
                    throw new ArgumentOutOfRangeException("Pressure out of range");
                }
            }
            else if (T > Constants.T23min && T <= Constants.Tmax)
            {
                if (p > Constants.Pmax)
                {
                    throw new ArgumentOutOfRangeException("Pressure out of range");
                }
                else if (p < 16.5292)
                { // Check this one first to avoid the call to 2-3 boundary curve (a little bit faster)
                    return IF97REGIONS.REGION_2;
                }
                else if (p > Region23_T(T))
                {
                    return IF97REGIONS.REGION_3;
                }
                else
                {
                    return IF97REGIONS.REGION_2;
                }
            }
            else if (T >= Constants.Tmin && T <= Constants.T23min)
            {
                if (p > Constants.Pmax)
                    throw new ArgumentOutOfRangeException("Pressure out of range");
                else if (p > Region4.p_T(T))
                    return IF97REGIONS.REGION_1;
                else if (p < Region4.p_T(T))
                    return IF97REGIONS.REGION_2;
                else
                    return IF97REGIONS.REGION_4;
            }
            else
            {
                throw new ArgumentOutOfRangeException("Temperature out of range");
            }
        }
        static readonly Region1 R1 = new Region1();
        static readonly Region2 R2 = new Region2();
        static readonly Region5 R5 = new Region5();
        public static double RegionOutput(IF97Parameters outkey, double T, double p, SatState State)
        {

            IF97REGIONS region = RegionDetermination_TP(T, p);

            switch (region)
            {
                case IF97REGIONS.REGION_1:
                    if (State == SatState.VAPOR)
                        return R2.output(outkey, T, p);  // On saturation curve and need the Vapor phase
                    else
                        return R1.output(outkey, T, p);  // otherwise, use Liquid Region 1
                case IF97REGIONS.REGION_2:
                    if (State == SatState.LIQUID)
                        return R1.output(outkey, T, p);  // On saturation curve and need the Liquid phase
                    else
                        return R2.output(outkey, T, p);  // otherwise, use Vapor Region 2
                case IF97REGIONS.REGION_3:
                    return Region3.output(outkey, T, p, State);
                case IF97REGIONS.REGION_4:
                    if (State == SatState.VAPOR)
                    {
                        return R2.output(outkey, T, p);
                    }
                    else if (State == SatState.LIQUID)
                    {
                        return R1.output(outkey, T, p);
                    }
                    else
                    {
                        throw new ArgumentOutOfRangeException("Cannot use Region 4 with T and p as inputs");
                    }
                case IF97REGIONS.REGION_5: return R5.output(outkey, T, p);
            }
            throw new ArgumentOutOfRangeException("Unable to match region");
        }
        static IF97REGIONS RegionDetermination_pX(double p, double X, IF97Parameters inkey)
        {
            // Setup needed Region Equations for region determination

            // Saturation Region Limit Variables
            double Tsat = 0;
            double Xliq = 0;
            double Xvap = 0;

            // Check overall boundary limits
            if ((p < Constants.Pmin) || (p > Constants.Pmax))
                throw new ArgumentOutOfRangeException("Pressure out of range");
            double Xmin = R1.output(inkey, Constants.Tmin, p);
            double Xmax = R2.output(inkey, Constants.Tmax, p);
            if ((X < Xmin) || (X > (Xmax + 1.0E-10)))
            {
                if (inkey == IF97Parameters.h)
                {
                    throw new ArgumentOutOfRangeException("Enthalpy out of range");
                }
                else
                {
                    throw new ArgumentOutOfRangeException("Entropy out of range");
                }
            }

            // Check saturation Dome first
            if (p <= Constants.Pcrit)
            {
                Tsat = Region4.Tsat97(p);
                Xliq = R1.output(inkey, Tsat, p);
                Xvap = R2.output(inkey, Tsat, p);
                if ((Xliq <= X) && (X <= Xvap))
                {    // Within Saturation Dome
                    return IF97REGIONS.REGION_4;               //    Region 4
                }
            }
            // End Check saturation Dome

            // Check values below 16.529 MPa
            if (p <= Constants.P23min)
            {                        // p <= P23min (saturation dome)
                if (X <= Xliq) return IF97REGIONS.REGION_1;
                else if (X >= Xvap) return IF97REGIONS.REGION_2;
                else return IF97REGIONS.REGION_4;
            }
            // Check values above 16.529 MPa
            else if (X <= R1.output(inkey, Constants.T23min, p))
                return IF97REGIONS.REGION_1;
            else if (X >= R2.output(inkey, Region23_p(p), p))
                return IF97REGIONS.REGION_2;
            else
                return IF97REGIONS.REGION_3;
        } // Region Output backward

        static int BackwardRegion(double p, double X, IF97Parameters inkey)
        {
            // This routine is for testing purposes only.  It returns the
            // Region as an integer based on the backward evaluation of either
            // (p,h) or (p,s)

            // Make sure input and output keys are valid for Backward formulas
            if ((inkey != IF97Parameters.h) && (inkey != IF97Parameters.s))
                throw new ArgumentException("Backward Formulas take variable inputs of Enthalpy or Entropy only.");

            IF97REGIONS region = RegionDetermination_pX(p, X, inkey);

            switch (region)
            {
                case IF97REGIONS.REGION_1: return 1;
                case IF97REGIONS.REGION_2: return 2;
                case IF97REGIONS.REGION_3: return 3;
                case IF97REGIONS.REGION_4: return 4;
                default: return 0;
            }
        }

        static readonly Region1H B1H = new Region1H();
        static readonly Region1S B1S = new Region1S();
        static readonly Region2aH B2aH = new Region2aH();
        static readonly Region2bH B2bH = new Region2bH();
        static readonly Region2cH B2cH = new Region2cH();
        static readonly Region2aS B2aS = new Region2aS();
        static readonly Region2bS B2bS = new Region2bS();
        static readonly Region2cS B2cS = new Region2cS();
        static readonly Region3aH B3aH = new Region3aH();
        static readonly Region3bH B3bH = new Region3bH();
        static readonly Region3aS B3aS = new Region3aS();
        static readonly Region3bS B3bS = new Region3bS();

        static double RegionOutputBackward(double p, double X, IF97Parameters inkey)
        {
            // Note that this routine returns only temperature (IF97_T).  All other values should be
            // calculated from this temperature and the known pressure using forward equations.
            // Setup Backward Regions for output

            // Make sure input and output keys are valid for Backward formulas
            if ((inkey != IF97Parameters.h) && (inkey != IF97Parameters.s))
                throw new ArgumentException("Backward Formulas take variable inputs of Enthalpy or Entropy only.");

            // Get Saturation Parameters

            IF97REGIONS region = RegionDetermination_pX(p, X, inkey);

            switch (region)
            {
                case IF97REGIONS.REGION_1:
                    if (inkey == IF97Parameters.h)
                        return B1H.T_pX(p, X);
                    else
                        return B1S.T_pX(p, X);
                case IF97REGIONS.REGION_2:
                    if (inkey == IF97Parameters.h)
                    {
                        if (p <= 4.0)
                            return B2aH.T_pX(p, X);
                        else if (X >= BackwardsRegion.H2b2c_p(p))
                            return B2bH.T_pX(p, X);
                        else
                            return B2cH.T_pX(p, X);
                    }
                    else
                    {
                        if (p <= 4.0)
                            return B2aS.T_pX(p, X);
                        else if (X >= Constants.S2bc)
                            return B2bS.T_pX(p, X);
                        else
                            return B2cS.T_pX(p, X);
                    };
                case IF97REGIONS.REGION_3:
                    if (inkey == IF97Parameters.h)
                    {
                        if (X <= BackwardsRegion.H3ab_p(p))
                            return B3aH.T_pX(p, X);
                        else
                            return B3bH.T_pX(p, X);
                    }
                    else
                    {
                        if (X <= Constants.Scrit)
                            return B3aS.T_pX(p, X);
                        else
                            return B3bS.T_pX(p, X);
                    };
                case IF97REGIONS.REGION_4: return Tsat97(p);
                default: throw new ArgumentOutOfRangeException("Unable to match region");
            }
        }  // Region Output backward

        static double rho_pX(double p, double X, IF97Parameters inkey)
        {
            // NOTE: This implementation works, and with the 2016 Supplementary Release
            //       for v(p,T) for Region 3 implemented, it is no longer iterative.  
            //       However, the 2014 Supplementary Release for v(p,h) and v(p,s) are 
            //       more direct and may be slightly faster, since only one algebraic 
            //       equation is needed instead of two in Region 3.
            double T = RegionOutputBackward(p, X, inkey);
            if (RegionDetermination_pX(p, X, inkey) == IF97REGIONS.REGION_4)
            {      // If in saturation dome
                double Tsat = Tsat97(p);
                double Xliq = R1.output(inkey, Tsat, p);
                double Xvap = R2.output(inkey, Tsat, p);
                double vliq = 1.0 / R1.output(IF97Parameters.d, Tsat, p);
                double vvap = 1.0 / R2.output(IF97Parameters.d, Tsat, p);
                return 1.0 / (vliq + (X - Xliq) * (vvap - vliq) / (Xvap - Xliq));  //    Return Mixture Density
            }
            else
            {                                                   // else
                return RegionOutput(IF97Parameters.d, T, p, SatState.NONE);
            }
        }

        static double Q_pX(double p, double X, IF97Parameters inkey)
        {
            double Xliq, Xvap;
            if ((p < Constants.Pmin) || (p > Constants.Pmax))
            {
                throw new ArgumentOutOfRangeException("Pressure out of range");
            }
            else if (p < Constants.Ptrip)
            {
                return 0;  //Liquid, at all temperatures
            }
            else if (p > Constants.Pcrit)
            {
                double t;
                switch (inkey)
                {
                    case IF97Parameters.h:
                    case IF97Parameters.s:
                        t = RegionOutputBackward(p, X, inkey); break;
                    case IF97Parameters.u:
                    case IF97Parameters.d:
                    default:
                        // There are no reverse functions for t(p,U) or t(p,rho)
                        throw new ArgumentException("Quality cannot be determined for these inputs.");
                }
                if (t < Constants.Tcrit)
                    return 1.0;  // Vapor, at all pressures above critical point
                else
                    // Supercritical Region (p>Pcrit) && (t>Tcrit)
                    throw new ArgumentException("Quality not defined in supercritical region.");
            }
            else
            {
                switch (inkey)
                {
                    case IF97Parameters.h:
                    case IF97Parameters.s:
                    case IF97Parameters.u:
                        Xliq = RegionOutput(inkey, Tsat97(p), p, SatState.LIQUID);
                        Xvap = RegionOutput(inkey, Tsat97(p), p, SatState.VAPOR);
                        return Math.Min(1.0, Math.Max(0.0, (X - Xliq) / (Xvap - Xliq)));
                    case IF97Parameters.d:
                        Xliq = 1.0 / RegionOutput(IF97Parameters.d, Tsat97(p), p, SatState.LIQUID);
                        Xvap = 1.0 / RegionOutput(IF97Parameters.d, Tsat97(p), p, SatState.VAPOR);
                        X = 1.0 / X;
                        return Math.Min(1.0, Math.Max(0.0, (X - Xliq) / (Xvap - Xliq)));
                    default:
                        throw new ArgumentException("Quality cannot be determined for these inputs.");
                }
            }
            // If all else fails, which it shouldn't...
            throw new ArgumentException("Quality cannot be determined for these inputs.");
        }
        static double X_pQ(IF97Parameters inkey, double p, double Q)
        {
            double Xliq, Xvap;
            if ((p < Constants.Ptrip) || (p > Constants.Pcrit))
                throw new ArgumentOutOfRangeException("Pressure out of range");
            if ((Q < 0.0) || (Q > 1.0))
                throw new ArgumentOutOfRangeException("Quality out of range");
            switch (inkey)
            {
                case IF97Parameters.h:
                case IF97Parameters.s:
                case IF97Parameters.u:
                    Xliq = RegionOutput(inkey, Tsat97(p), p, SatState.LIQUID);
                    Xvap = RegionOutput(inkey, Tsat97(p), p, SatState.VAPOR);
                    return Q * Xvap + (1 - Q) * Xliq;
                case IF97Parameters.d:
                    Xliq = 1.0 / RegionOutput(IF97Parameters.d, Tsat97(p), p, SatState.LIQUID);
                    Xvap = 1.0 / RegionOutput(IF97Parameters.d, Tsat97(p), p, SatState.VAPOR);
                    return 1.0 / (Q * Xvap + (1 - Q) * Xliq);
                default:
                    throw new ArgumentException("Mixture property undefined");
            }
        }


        static double[] HTmaxdata = {
            1.00645619394616E4,
            1.94706669580164E5,
            -4.67105212810782E5,
            -3.38175262587035E4
        };

        static double Hmax(double s)
        {
            // This function covers the top and right domain boundaries of constant Pmax and Tmax
            double s_star = 1, h_star = 1, sigma = s / s_star;
            if (s < Constants.STPmax)  // Use forward equation along Pmax using T(Pmax,s) as Temperature
                return RegionOutput(IF97Parameters.h, RegionOutputBackward(Constants.Pmax, s, IF97Parameters.s), Constants.Pmax, SatState.NONE);
            else
            {
                // Determining H(s) along Tmax is difficult because there is no direct p(T,s) formulation.
                // This linear combination fit h(s)=a*ln(s)+b/s+c/s²+d is not perfect, but it's close
                // and can serve as a limit along that Tmax boundary. Coefficients in HTmaxdata above.
                // There is a better way to do this using Newton-Raphson on Tmax = T(p,s), but it is iterative and slow.
                double ETA = HTmaxdata[0] * Math.Log(sigma) + HTmaxdata[1] / sigma + HTmaxdata[2] / Math.Pow(sigma, 2) + HTmaxdata[3];
                return ETA * h_star;
            }
        }

        static double Hmin(double s)
        {
            if (s < Constants.Sgtrip)  // Interpolate through Region 4
                return (s - Constants.Sftrip) * (Constants.Hgtrip - Constants.Hftrip) / (Constants.Sgtrip - Constants.Sftrip) + Constants.Hftrip;
            else             // Use forward equation along Pmin using T(Pmin,s) as Temperature
                return RegionOutput(IF97Parameters.h, RegionOutputBackward(Constants.Pmin, s, IF97Parameters.s), Constants.Pmin, SatState.NONE);
        }

        static Boundary13HS b13 = new Boundary13HS();
        static Boundary23HS b23hs = new Boundary23HS();
        static Region2cHS R2c = new Region2cHS();
        static IF97BACKREGIONS RegionDetermination_HS(double h, double s)
        {

            // Check Overall Boundaries
            if ((s < Constants.Smin) || (s > Constants.Smax))
                throw new ArgumentOutOfRangeException("Entropy out of range");
            if ((h > Hmax(s)) || (h < Hmin(s)))
                throw new ArgumentOutOfRangeException("Enthalpy out of range");

            // ============================================================================
            // Start at the low entropy curves and work our way up.
            // =================================== Region 1 Check =========================
            if (s <= Constants.SfT23)
            {
                if (h < BackwardsRegion.Hsat_s(s))     //   If below Saturated Liquid Curve
                    return IF97BACKREGIONS.BACK_4;                //       REGION 4
                else if (s < Constants.S13min)              //   If below H13 Curve
                    return IF97BACKREGIONS.BACK_1;                //       REGION 1
                else
                {                            //   IF within H13 Curve (S13min < s < SfT23)
                    if (h < b13.h_s(s))           //       below curve
                        return IF97BACKREGIONS.BACK_1;            //           REGION 1
                    else                          //       above curve
                        return IF97BACKREGIONS.BACK_3A;           //           REGION 3 
                }                            //
            }
            else if (s <= Constants.Scrit)
            {  //========== Region 3a Check (S < Scrit) ==============
                if (h < BackwardsRegion.Hsat_s(s))     //  If below Saturated Liquid Curve
                    return IF97BACKREGIONS.BACK_4;                //      REGION 4
                else                              //  If above curve
                    return IF97BACKREGIONS.BACK_3A;               //      REGION 3(a)
            }
            else if (s <= Constants.S23min)
            {  //========== Region 3b Check            ==============
                if (h < BackwardsRegion.Hsat_s(s))     //  If below Saturated Liquid Curve
                    return IF97BACKREGIONS.BACK_4;                //      REGION 4
                else                              //  If above curve
                    return IF97BACKREGIONS.BACK_3B;               //      REGION 3(a)
            }
            else if (s <= Constants.S23max)
            {  //========== Region 3b/2c Check Along B23 Curve ======
                if (h < BackwardsRegion.Hsat_s(s))     //  if below Saturated Vapor Curve
                    return IF97BACKREGIONS.BACK_4;                //      REGION 4
                else if (h < Constants.H23min)              //  if below bounding box
                    return IF97BACKREGIONS.BACK_3B;               //      REGION 3(b)
                else if (h > Constants.H23max)              //  if above bounding box
                    return IF97BACKREGIONS.BACK_2C;               //      REGION 2(c)
                else
                {                            //  Need to check TB23 Curve
                    double TB23 = b23hs.t_hs(h, s);       //  Calc TB23(h,s)
                    double PB23 = Region23_T(TB23);      //  Calc Corresponding PB23
                    double P = R2c.p_hs(h, s);         //  Calc P(h,s) using Region 2c
                    if (P > PB23)                        //  Above B23 Curve
                        return IF97BACKREGIONS.BACK_3B;                  //      REGION 3(b)
                    else                                 //  Below B23 Curve
                        return IF97BACKREGIONS.BACK_2C;                  //      REGION 2(c)
                }
            }
            else if (s <= Constants.S2bc)
            {   //========== Region 3b Check            ==============
                if (h < BackwardsRegion.Hsat_s(s))     //  If below Saturated Liquid Curve
                    return IF97BACKREGIONS.BACK_4;                //      REGION 4
                else                              //  If above curve
                    return IF97BACKREGIONS.BACK_2C;               //      REGION 2(c)
            }
            else if (s < Constants.Sgtrip)
            { //========== Region 2a/2b (s > S2bc) above Sat. Curve ==
                if (h < BackwardsRegion.Hsat_s(s))     //  If below Saturated Vapor Curve
                    return IF97BACKREGIONS.BACK_4;                //      REGION_4
                else
                {                            //  If above Curve then
                    if (h > BackwardsRegion.H2ab_s(s))   //      if h > h2ab(s) Curve (P=4 MPa)
                        return IF97BACKREGIONS.BACK_2B;             //          REGION 2(b)
                    else                            //      if h < h2ab(s) Curve (P=4 MPa)
                        return IF97BACKREGIONS.BACK_2A;             //          REGION 2(a)
                }
            }
            else
                return IF97BACKREGIONS.BACK_2A;      //========== Region 2a fall thru ========================
        }// Region Determination HS

        static Region1HS B1HS = new Region1HS();
        static Region2aHS B2aHS = new Region2aHS();
        static Region2bHS B2bHS = new Region2bHS();
        static Region2cHS B2cHS = new Region2cHS();
        static Region3aHS B3aHS = new Region3aHS();
        static Region3bHS B3bHS = new Region3bHS();
        static Region4HS B4HS = new Region4HS();
        static double BackwardOutputHS(IF97Parameters outkey, double h, double s)
        {
            // Note that this routine returns only temperature (IF97_T).  All other values should be
            // calculated from this temperature and the known pressure using forward equations.
            // Setup Backward Regions for output
            //
            double Pval = 0, Tval = 0;

            // Make sure output keys are valid for Backward_HS formulas
            if ((outkey != IF97Parameters.p) && (outkey != IF97Parameters.T))
                throw new ArgumentException("Backward HS Formulas output Temperature or Pressure only.");

            // Get Saturation Parameters

            IF97BACKREGIONS region = RegionDetermination_HS(h, s);

            switch (region)
            {
                case IF97BACKREGIONS.BACK_1: Pval = B1HS.p_hs(h, s); break;
                case IF97BACKREGIONS.BACK_2A: Pval = B2aHS.p_hs(h, s); break;
                case IF97BACKREGIONS.BACK_2B: Pval = B2bHS.p_hs(h, s); break;
                case IF97BACKREGIONS.BACK_2C: Pval = B2cHS.p_hs(h, s); break;
                case IF97BACKREGIONS.BACK_3A: Pval = B3aHS.p_hs(h, s); break;
                case IF97BACKREGIONS.BACK_3B: Pval = B3bHS.p_hs(h, s); break;
                case IF97BACKREGIONS.BACK_4:
                    if (s >= Constants.SgT23)   // T(h,s) only defined over part of the 2-phase region
                        Tval = B4HS.t_hs(h, s);
                    else
                        throw new ArgumentOutOfRangeException("Entropy out of range");
                    break;
                default: throw new ArgumentOutOfRangeException("Unable to match region");
            }
            if (outkey == IF97Parameters.p)            // Returning Pressure (IF97_P)
                if (region == IF97BACKREGIONS.BACK_4)        //
                    return psat97(Tval);     //       Not REGION 4, already have pressure
                else                         // 
                    return Pval;             //       REGION 4, Calculate Psat from Tsat
                                             //
            else                             // ELSE Returning Temperature
                if (region == IF97BACKREGIONS.BACK_4)        //
                return Tval;                                     // REGION 4, already have Temperature
            else                                                 //
                return RegionOutputBackward(Pval, h, IF97Parameters.h);  // Not REGION 4 Calc from Backward T(p,h)
        }  // Region Output backward

        // ******************************************************************************** //
        //                                     API                                          //
        // ******************************************************************************** //

        /// Get the mass density [kg/m^3] as a function of T [K] and p [Pa]
        public static double rhomass_Tp(double T, double p) { return RegionOutput(IF97Parameters.d, T, p, SatState.NONE); }
        /// Get the mass enthalpy [J/kg] as a function of T [K] and p [Pa]
        public static double hmass_Tp(double T, double p) { return RegionOutput(IF97Parameters.h, T, p, SatState.NONE); }
        /// Get the mass entropy [J/kg/K] as a function of T [K] and p [Pa]
        public static double smass_Tp(double T, double p) { return RegionOutput(IF97Parameters.s, T, p, SatState.NONE); }
        /// Get the mass internal energy [J/kg] as a function of T [K] and p [Pa]
        public static double umass_Tp(double T, double p) { return RegionOutput(IF97Parameters.u, T, p, SatState.NONE); }
        /// Get the mass constant-pressure specific heat [J/kg/K] as a function of T [K] and p [Pa]
        public static double cpmass_Tp(double T, double p) { return RegionOutput(IF97Parameters.cp, T, p, SatState.NONE); }
        /// Get the mass constant-volume specific heat [J/kg/K] as a function of T [K] and p [Pa]
        public static double cvmass_Tp(double T, double p) { return RegionOutput(IF97Parameters.cv, T, p, SatState.NONE); }
        /// Get the speed of sound [m/s] as a function of T [K] and p [Pa]
        public static double speed_sound_Tp(double T, double p) { return RegionOutput(IF97Parameters.w, T, p, SatState.NONE); }
        /// Get the [d(rho)/d(p)]T [kg/m³/Pa] as a function of T [K] and p [Pa]
        public static double drhodp_Tp(double T, double p) { return RegionOutput(IF97Parameters.drhodp, T, p, SatState.NONE); }

        // ******************************************************************************** //
        //                            Transport Properties                                  //
        // ******************************************************************************** //

        /// Get the viscosity [Pa-s] as a function of T [K] and Rho [kg/m³]
        static double visc_TRho(double T, double rho)
        {
            // Since we have density, we don't need to determine the region for viscosity.
            // All regions use base region equations for visc(T,rho).
            return R1.visc(T, rho);
        }
        /// Get the viscosity [Pa-s] as a function of T [K] and p [Pa]
        static double visc_Tp(double T, double p) { return RegionOutput(IF97Parameters.mu, T, p, SatState.NONE); }
        /// Get the thermal conductivity [W/m-K] as a function of T [K] and p [Pa]
        static double tcond_Tp(double T, double p) { return RegionOutput(IF97Parameters.k, T, p, SatState.NONE); }
        /// Calculate the Prandtl number [dimensionless] as a function of T [K] and p [Pa]
        static double prandtl_Tp(double T, double p)
        {
            return visc_Tp(T, p) * cpmass_Tp(T, p) * 1000 / tcond_Tp(T, p);
        }

        // ******************************************************************************** //
        //                             Saturated Vapor/Liquid Functions                     //
        // ******************************************************************************** //
        /// Get the saturated liquid mass density [kg/m^3] as a function of p [Pa]
        public static double rholiq_p(double p) { return RegionOutput(IF97Parameters.d, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor mass density [kg/m^3] as a function of p [Pa]
        public static double rhovap_p(double p) { return RegionOutput(IF97Parameters.d, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid mass enthalpy [J/kg] as a function of p [Pa]
        public static double hliq_p(double p) { return RegionOutput(IF97Parameters.h, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor mass enthalpy [J/kg] as a function of p [Pa]
        public static double hvap_p(double p) { return RegionOutput(IF97Parameters.h, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid mass entropy [J/kg/K] as a function of p [Pa]
        public static double sliq_p(double p) { return RegionOutput(IF97Parameters.s, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor mass entropy [J/kg/K] as a function of p [Pa]
        public static double svap_p(double p) { return RegionOutput(IF97Parameters.s, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid mass internal energy [J/kg] as a function of p [Pa]
        public static double uliq_p(double p) { return RegionOutput(IF97Parameters.u, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor mass internal energy [J/kg] as a function of p [Pa]
        public static double uvap_p(double p) { return RegionOutput(IF97Parameters.u, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid mass isobaric specific heat [J/kg/K] as a function of T [K] and p [Pa]
        public static double cpliq_p(double p) { return RegionOutput(IF97Parameters.cp, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor mass isobaric specific heat [J/kg/K] as a function of T [K] and p [Pa]
        public static double cpvap_p(double p) { return RegionOutput(IF97Parameters.cp, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid mass isochoric specific heat [J/kg/K] as a function of T [K] and p [Pa]
        public static double cvliq_p(double p) { return RegionOutput(IF97Parameters.cv, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor mass isochoric specific heat [J/kg/K] as a function of T [K] and p [Pa]
        public static double cvvap_p(double p) { return RegionOutput(IF97Parameters.cv, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid speed of sound [m/s] as a function of T [K] and p [Pa]
        public static double speed_soundliq_p(double p) { return RegionOutput(IF97Parameters.w, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor speed of sound [m/s] as a function of T [K] and p [Pa]
        public static double speed_soundvap_p(double p) { return RegionOutput(IF97Parameters.w, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid viscosity [Pa-s] as a function of p [Pa]
        public static double viscliq_p(double p) { return RegionOutput(IF97Parameters.mu, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor viscosity [Pa-s] as a function of p [Pa]
        public static double viscvap_p(double p) { return RegionOutput(IF97Parameters.mu, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Get the saturated liquid thermal conductivity [W/m-K] as a function of p [Pa]
        public static double tcondliq_p(double p) { return RegionOutput(IF97Parameters.k, Tsat97(p), p, SatState.LIQUID); }
        /// Get the saturated vapor thermal conductivity [W/m-K] as a function of p [Pa]
        public static double tcondvap_p(double p) { return RegionOutput(IF97Parameters.k, Tsat97(p), p, SatState.VAPOR); }
        // ******************************************************************************** //
        /// Calculate the saturated liquid Prandtl number [dimensionless] as a function of p [Pa]
        public static double prandtlliq_p(double p) { return viscliq_p(p) * cpliq_p(p) * 1000 / tcondliq_p(p); }
        /// Calculate the saturated vapor Prandtl number [dimensionless] as a function of p [Pa]
        public static double prandtlvap_p(double p) { return viscvap_p(p) * cpvap_p(p) * 1000 / tcondvap_p(p); }


        // ******************************************************************************** //
        //                               2-Phase Functions                                  //
        // ******************************************************************************** //
        /// Get the saturation temperature [K] as a function of p [Pa]
        public static double Tsat97(double p)
        {
            return Region4.T_p(p);
        }
        /// Get the saturation pressure [Pa] as a function of T [K]
        public static double psat97(double T)
        {
            return Region4.p_T(T);
        }
        /// Get surface tension [N/m] as a function of T [K]
        public static double sigma97(double T)
        {
            return Region4.sigma_t(T);
        }
        // ******************************************************************************** //
        //                              Backward Functions                                  //
        // ******************************************************************************** //
        public static double T_phmass(double p, double h)
        {
            return RegionOutputBackward(p, h, IF97Parameters.h);
        }
        public static double rhomass_phmass(double p, double h)
        {
            return rho_pX(p, h, IF97Parameters.h);
        }
        public static double T_psmass(double p, double s)
        {
            return RegionOutputBackward(p, s, IF97Parameters.s);
        }
        public static double rhomass_psmass(double p, double s)
        {
            return rho_pX(p, s, IF97Parameters.s);
        }
        public static double p_hsmass(double h, double s)
        {
            return BackwardOutputHS(IF97Parameters.p, h, s);
        }
        public static double T_hsmass(double h, double s)
        {
            return BackwardOutputHS(IF97Parameters.T, h, s);
        }
        public static int Region_ph(double p, double h)
        {
            return BackwardRegion(p, h, IF97Parameters.h);
        }
        public static int Region_ps(double p, double s)
        {
            return BackwardRegion(p, s, IF97Parameters.s);
        }
        // ******************************************************************************** //
        //                              Trivial Functions                                   //
        // ******************************************************************************** //
        /// Get the Triple Point Temperature and Pressure
        public static double Ttrip { get { return Constants.Ttrip; } }
        public static double ptrip { get { return Constants.Ptrip; } }
        /// Get the Critical Point Temperature and Pressure and Density
        public static double Tcrit { get { return Constants.Tcrit; } }
        public static double pcrit { get { return Constants.Pcrit; } }
        public static double rhocrit { get { return Constants.Rhocrit; } }
        /// Get the Max and Min Temperatures and Pressures
        public static double  Tmin { get { return Constants.Tmin; } }
        public static double Pmin { get { return Constants.Pmin; } }
        public static double Tmax { get { return Constants.Tmax; } }
        public static double Pmax { get { return Constants.Pmax; } }
        /// Get physical constants
        public static double MW { get { return Constants.MW; } }
        public static double Rgas { get { return Constants.Rgas; } }
        public static double Acentric { get { return -Math.Log(psat97(0.7 * Constants.Tcrit) / Constants.Pcrit) - 1; } }
        // ******************************************************************************** //
        //                              Utility Functions                                   //
        // ******************************************************************************** //
        const string IF97VERSION = "v2.1.2";
        public static string get_if97_version()
        {
#if IAPWS_UNITS
            return IF97VERSION + " (IAPWS Units)";
#else
            return IF97VERSION;
#endif
        }
        public static double hmass_pQ(double p, double Q)
        {
            return X_pQ(IF97Parameters.h, p, Q);
        }
        public static double umass_pQ(double p, double Q)
        {
            return X_pQ(IF97Parameters.u, p, Q);
        }
        public static double smass_pQ(double p, double Q)
        {
            return X_pQ(IF97Parameters.s, p, Q);
        }
        public static double v_pQ(double p, double Q)
        {
            return 1.0 / X_pQ(IF97Parameters.d, p, Q);
        }
        public static double rhomass_pQ(double p, double Q)
        {
            return X_pQ(IF97Parameters.d, p, Q);
        }
        public static double Q_phmass(double p, double h)
        {
            return Q_pX(p, h, IF97Parameters.h);
        }
        public static double Q_pumass(double p, double u)
        {
            return Q_pX(p, u, IF97Parameters.u);
        }
        public static double Q_psmass(double p, double s)
        {
            return Q_pX(p, s, IF97Parameters.s);
        }
        public static double Q_prhomass(double p, double rho)
        {
            return Q_pX(p, rho, IF97Parameters.d);
        }
        public static double Q_pv(double p, double v)
        {
            return Q_pX(p, 1.0 / v, IF97Parameters.d);
        }
        /*************************************************************************/
        /* Vapor Quality as a function of H and S cannot be implemented at this  */
        /* time as there IAPWS has not released a backward formula that covers   */
        /* the entire vapor dome.  T(H,S) is only currently defined for the      */
        /* range of s > s"(T23sat).  Leaving this code here in case IAPWS        */
        /* releases a fit for the entire vapor dome.                             */
        /*************************************************************************/
        /*  double Q_hsmass(double h, double s){
                double hliq, hvap;
                if ((s < Global.Smin) || (s > Global.Smax)) {
                    throw new ArgumentOutOfRangeException("Entropy out of range");
                } else if ((h < Hmin(s)) || (h > Hmax(s))) {
                    throw new ArgumentOutOfRangeException("Enthalpy out of range");
                }
                double p = BackwardOutputHS(Parameters.p, h, s);
                double t = RegionOutputBackward( p, h, Parameters.h);
                if ((p > Global.Pcrit) && (t > Global.Tcrit)) {
                    throw new ArgumentOutOfRangeException("Temperature out of range");
                } else if (BackwardRegion(p, h, Parameters.h) == 4) {
                    hliq = RegionOutput(Parameters.h, Tsat97(p), p, SatState.LIQUID);
                    hvap = RegionOutput(Parameters.h, Tsat97(p), p, SatState.VAPOR);
                    return Math.Min(1.0,Math.Max(0.0,(h-hliq)/(hvap-hliq)));
                } else if (p > Global.Pcrit) {
                    return 0.0;
                } else {
                    return 1.0;
                }
            }                                                                 */
        /*************************************************************************/

        /* namespace IF97 */

#if ENABLE_CATCH

        static ValueTuple<char,double,double,double>[] Table5 = {
            // Table 5 data 
            ('A', 630, 50e6, 1.470853100e-3),
            ('A', 670, 80e6, 1.503831359e-3),
            ('B', 710, 50e6, 2.204728587e-3),
            ('B', 750, 80e6, 1.973692940e-3),
            ('C', 630, 20e6, 1.761696406e-3),
            ('C', 650, 30e6, 1.819560617e-3),
            ('D', 656, 26e6, 2.245587720e-3),
            ('D', 670, 30e6, 2.506897702e-3),
            ('E', 661, 26e6, 2.970225962e-3),
            ('E', 675, 30e6, 3.004627086e-3),
            ('F', 671, 26e6, 5.019029401e-3),
            ('F', 690, 30e6, 4.656470142e-3),
            ('G', 649, 23.6e6, 2.163198378e-3),
            ('G', 650, 24e6, 2.166044161e-3),
            ('H', 652, 23.6e6, 2.651081407e-3),
            ('H', 654, 24e6, 2.967802335e-3),
            ('I', 653, 23.6e6, 3.273916816e-3),
            ('I', 655, 24e6, 3.550329864e-3),
            ('J', 655, 23.5e6, 4.545001142e-3),
            ('J', 660, 24e6, 5.100267704e-3),
            ('K', 660, 23e6, 6.109525997e-3),
            ('K', 670, 24e6, 6.427325645e-3),
            ('L', 646, 22.6e6, 2.117860851e-3),
            ('L', 646, 23e6, 2.062374674e-3),
            ('M', 648.6, 22.6e6, 2.533063780e-3),
            ('M', 649.3, 22.8e6, 2.572971781e-3),
            ('N', 649.0, 22.6e6, 2.923432711e-3),
            ('N', 649.7, 22.8e6, 2.913311494e-3),
            ('O', 649.1, 22.6e6, 3.131208996e-3),
            ('O', 649.9, 22.8e6, 3.221160278e-3),
            ('P', 649.4, 22.6e6, 3.715596186e-3),
            ('P', 650.2, 22.8e6, 3.664754790e-3),
            ('Q', 640, 21.1e6, 1.970999272e-3),
            ('Q', 643, 21.8e6, 2.043919161e-3),
            ('R', 644, 21.1e6, 5.251009921e-3),
            ('R', 648, 21.8e6, 5.256844741e-3),
            ('S', 635, 19.1e6, 1.932829079e-3),
            ('S', 638, 20e6, 1.985387227e-3),
            ('T', 626, 17e6, 8.483262001e-3),
            ('T', 640, 20e6, 6.227528101e-3),
            // Table 11
            ('U', 644.6, 21.5e6, 2.268366647e-3),
            ('U', 646.1, 22e6, 2.296350553e-3),
            ('V', 648.6, 22.5e6, 2.832373260e-3),
            ('V', 647.9, 22.3e6, 2.811424405e-3),
            ('W', 647.5, 22.15e6, 3.694032281e-3),
            ('W', 648.1, 22.3e6, 3.622226305e-3),
            ('X', 648, 22.11e6, 4.528072649e-3),
            ('X', 649, 22.3e6, 4.556905799e-3),
            ('Y', 646.84, 22e6, 2.698354719e-3),
            ('Y', 647.05, 22.064e6, 2.717655648e-3),
            ('Z', 646.89, 22e6, 3.798732962e-3),
            ('Z', 647.15, 22.064e6, 3.701940010e-3)
        };


        void print_IF97_Table5()
        {
            for (int i = 0; i < Table5.Length; ++i){
                double v = Region3Backwards.Region3_v_TP(Table5[i].Item1, Table5[i].Item2, Table5[i].Item3);
                if (Math.Abs(v - Table5[i].Item4) > 1e-12){
                        Console.WriteLine(String.Format("{0} {1} {2} {3} {4}", Table5[i].Item1, Table5[i].Item2, Table5[i].Item3, v, Table5[i].Item4);
                }
                char region = Region3Backwards.RegionDetermination(Table5[i].Item2, Table5[i].Item3);
                if (region != Table5[i].Item1)
                {
                        Region3Backwards.RegionDetermination(Table5[i].Item2, Table5[i].Item3);
                        Console.WriteLine(String.Format("{0} {1} {2} {3} {4}", Table5[i].Item1, Table5[i].Item2, Table5[i].Item3, v, Table5[i].Item4);
                }
            }
        }


        static ValueTuple<LineEnum,double,double>[] Table3 = {    
            // Table 3
            (LineEnum.AB, 40e6, 6.930341408e2),
            (LineEnum.CD, 25e6, 6.493659208e2),
            (LineEnum.EF, 40e6, 7.139593992e2),
            (LineEnum.GH, 23e6, 6.498873759e2),
            (LineEnum.IJ, 23e6, 6.515778091e2),
            (LineEnum.JK, 23e6, 6.558338344e2),
            (LineEnum.MN, 22.8e6, 6.496054133e2),
            (LineEnum.OP, 22.8e6, 6.500106943e2),
            (LineEnum.QU, 22e6, 6.456355027e2),
            (LineEnum.RX, 22e6, 6.482622754e2),
            // Table 11
            (LineEnum.UV, 22.3e6, 6.477996121e2),
            (LineEnum.WX, 22.3e6, 6.482049480e2)
        };

        void print_boundary_line_Table3()
        {
            for (int i = 0; i < Table3.Length; ++i){
                double T = Region3Backwards.DividingLine(Table3[i].Item1, Table3[i].Item2);
                if (Math.Abs(T - Table3[i].Item3) > 1e-7){
                    Console.WriteLine("p: {0} errT: {1}", Table3[i].Item2, T - Table3[i].Item3);
                }
            }
        }

#endif
    }
    enum IF97REGIONS { REGION_1, REGION_2, REGION_3, REGION_4, REGION_5 }
    enum IF97BACKREGIONS { BACK_1, BACK_2A, BACK_2B, BACK_2C, BACK_3A, BACK_3B, BACK_4 }
}
