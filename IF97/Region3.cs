using System;

namespace IF97
{
    public class Region3
    {
        static readonly ValueTuple<int, int, double>[] Region3residdata = {
            (0, 0,    0.10658070028513e1),
            (0, 0,   -0.15732845290239e2),
            (0, 1,    0.20944396974307e2),
            (0, 2,   -0.76867707878716e1),
            (0, 7,    0.26185947787954e1),
            (0, 10,  -0.28080781148620e1),
            (0, 12,   0.12053369696517e1),
            (0, 23,  -0.84566812812502e-2),
            (1, 2,   -0.12654315477714e1),
            (1, 6,   -0.11524407806681e1),
            (1, 15,   0.88521043984318),
            (1, 17,  -0.64207765181607),
            (2, 0,    0.38493460186671),
            (2, 2,   -0.85214708824206),
            (2, 6,    0.48972281541877e1),
            (2, 7,   -0.30502617256965e1),
            (2, 22,   0.39420536879154e-1),
            (2, 26,   0.12558408424308),
            (3, 0,   -0.27999329698710),
            (3, 2,    0.13899799569460e1),
            (3, 4,   -0.20189915023570e1),
            (3, 16,  -0.82147637173963e-2),
            (3, 26,  -0.47596035734923),
            (4, 0,    0.43984074473500e-1),
            (4, 2,   -0.44476435428739),
            (4, 4,    0.90572070719733),
            (4, 26,   0.70522450087967),
            (5, 1,    0.10770512626332),
            (5, 3,   -0.32913623258954),
            (5, 26,  -0.50871062041158),
            (6, 0,   -0.22175400873096e-1),
            (6, 2,    0.94260751665092e-1),
            (6, 26,   0.16436278447961),
            (7, 2,   -0.13503372241348e-1),
            (8, 26,  -0.14834345352472e-1),
            (9, 2,    0.57922953628084e-3),
            (9, 26,   0.32308904703711e-2),
            (10, 0,   0.80964802996215e-4),
            (10, 1,  -0.16557679795037e-3),
            (11, 26, -0.44923899061815e-4)
        };
        static int[] Ir, Jr;
        static double[] nr;
        /// For Viscosity Calculations
        static int[] muJ0;
        static double[] mun0;
        static int[] muIr, muJr;
        static double[] munr;
        /// For Thermal Conductivity Calculations
        static int[] lamJ0;
        static double[] lamn0;
        static int[] lamIr, lamJr;
        static double[] lamnr;
        static readonly double R;

        static Region3()
        {
            int len = Region3residdata.Length;
            if (len > 0)
            {
                Ir = new int[len];
                Jr = new int[len];
                nr = new double[len];
                for (int i = 0; i < len; i++)
                {
                    var tuple = Region3residdata[i];
                    Ir[i] = tuple.Item1;
                    Jr[i] = tuple.Item2;
                    nr[i] = tuple.Item3;
                }
            }

            len = Constants.Hresiddata.Length;
            muIr = new int[len];
            muJr = new int[len];
            munr = new double[len];
            for (int i = 0; i < len; i++)
            {
                var tuple = Constants.Hresiddata[i];
                muIr[i] = tuple.Item1;
                muJr[i] = tuple.Item2;
                munr[i] = tuple.Item3;
            }
            len = Constants.Hidealdata.Length;
            mun0 = new double[len];
            muJ0 = new int[len];
            for (int i = 0; i < len; i++)
            {
                var tuple = Constants.Hidealdata[i];
                muJ0[i] = tuple.Item1;
                mun0[i] = tuple.Item2;
            }
            len = Constants.Lresiddata.Length;
            lamIr = new int[len];
            lamJr = new int[len];
            lamnr = new double[len];
            for (int i = 0; i < len; i++)
            {
                var tuple = Constants.Lresiddata[i];
                lamIr[i] = tuple.Item1;
                lamJr[i] = tuple.Item2;
                lamnr[i] = tuple.Item3;
            }
            len = Constants.Lidealdata.Length;
            lamn0 = new double[len];
            lamJ0 = new int[len];
            for (int i = 0; i < len; i++)
            {
                var tuple = Constants.Lidealdata[i];
                lamJ0[i] = tuple.Item1;
                lamn0[i] = tuple.Item2;
            }
            R = Constants.Rgas;
        }

        static double phi(double T, double rho)
        {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
            double summer = nr[0] * Math.Log(delta);
            for (int i = 1; i < 40; ++i)
            {
                summer += nr[i] * FastPow.Pow(delta, Ir[i]) * FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
#if REGION3_ITERATE
        //
        // These two extra terms Needed to evaluate Newton-Raphson
        // ****************************************************************************
        double dphi_ddelta(double T, double rho) {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
        double summer = nr[0] / delta;
            for (int i = 1; i< 40; ++i){
                summer += nr[i]* Ir[i]* FastPow.Pow(delta, Ir[i]-1)* FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
        double d2phi_ddelta2(double T, double rho) {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
        double summer = -nr[0] / (delta * delta);
            for (int i = 1; i< 40; ++i){
                summer += nr[i]* Ir[i]* (Ir[i]-1.0)* FastPow.Pow(delta, Ir[i]-2)* FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
        // ****************************************************************************
#endif
        static double delta_dphi_ddelta(double T, double rho)
        {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
            double summer = nr[0];
            for (int i = 1; i < 40; ++i)
            {
                summer += nr[i] * Ir[i] * FastPow.Pow(delta, Ir[i]) * FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
        static double tau_dphi_dtau(double T, double rho)
        {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
            double summer = 0;
            for (int i = 1; i < 40; ++i)
            {
                summer += nr[i] * Jr[i] * FastPow.Pow(delta, Ir[i]) * FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
        static double delta2_d2phi_ddelta2(double T, double rho)
        {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
            double summer = -nr[0];
            for (int i = 1; i < 40; ++i)
            {
                summer += nr[i] * Ir[i] * (Ir[i] - 1) * FastPow.Pow(delta, Ir[i]) * FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
        static double tau2_d2phi_dtau2(double T, double rho)
        {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
            double summer = 0;
            for (int i = 1; i < 40; ++i)
            {
                summer += nr[i] * Jr[i] * (Jr[i] - 1) * FastPow.Pow(delta, Ir[i]) * FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
        static double deltatau_d2phi_ddelta_dtau(double T, double rho)
        {
            double delta = rho / Constants.Rhocrit, tau = Constants.Tcrit / T;
            double summer = 0;
            for (int i = 1; i < 40; ++i)
            {
                summer += nr[i] * Jr[i] * Ir[i] * FastPow.Pow(delta, Ir[i]) * FastPow.Pow(tau, Jr[i]);
            }
            return summer;
        }
        static double mu0(double T)
        {
            double T_bar = T / Constants.Tcrit;
            double summer = 0.0;
            for (int i = 0; i < muJ0.Length; ++i)
            {
                summer += mun0[i] / FastPow.Pow(T_bar, muJ0[i]);
            }
            return 100.0 * Math.Sqrt(T_bar) / summer;
        }
        static double mu1(double T, double rho)
        {
            double rho_bar = rho / Constants.Rhocrit;
            double summer = 0.0;
            for (int i = 0; i < muJr.Length; ++i)
            {
                summer += rho_bar * FastPow.Pow(Trterm(T), muIr[i]) * munr[i] * FastPow.Pow(Rhorterm(rho), muJr[i]);
            }
            return Math.Exp(summer);
        }
        static double lambda0(double T)
        {
            double T_bar = T / Constants.Tcrit;
            double summer = 0.0;
            for (int i = 0; i < lamJ0.Length; ++i)
            {
                summer += lamn0[i] / FastPow.Pow(T_bar, lamJ0[i]);
            }
            return Math.Sqrt(T_bar) / summer;
        }
        static double lambda1(double T, double rho)
        {
            double rho_bar = rho / Constants.Rhocrit;
            double summer = 0.0;
            for (int i = 0; i < lamJr.Length; ++i)
            {
                summer += rho_bar * FastPow.Pow(Trterm(T), lamIr[i]) * lamnr[i] * FastPow.Pow(Rhorterm(rho), lamJr[i]);
            }
            return Math.Exp(summer);
        }
        static double lambda2(double T, double _, double rho)
        {
            double y, Cpbar, mubar, k, Z, zeta, delChi, Cpcalc;
            double rhobar = rho / Constants.Rhocrit;   /// Dimensionless
            const double LAMBDA = 177.8514;      /// Dimensionless
            const double qD = 1.0 / 0.40;      /// 1/nm
            const double Tr = 1.5 * Constants.Tcrit;     /// Dimensionless
            const double xi0 = 0.13;          /// nm
            const double nu = 0.630;         /// Dimensionless
            const double gam = 1.239;         /// Dimensionless
            const double GAMMA0 = 0.06;          /// Dimensionless
            double PI = 2 * Math.Acos(0.0);   /// Have to define this in C++
            const double Cpstar = 0.46151805;  /// Note: Slightly lower than IF97 Rgas  {J/kg-K}
            Cpcalc = cpmass(T, rho);                                  /// J/kg-K
            Cpbar = Cpcalc / Cpstar;                                   /// Unit-less
            if ((Cpbar < 0) || (Cpbar > 1.0E13)) Cpbar = 1.0E13;     /// Unit-less
            k = Cpcalc / cvmass(T, rho);                                /// Unit-less
            mubar = visc(T, rho) / 1.0E-6;                              /// Unit-less
            zeta = Constants.Pcrit / Constants.Rhocrit * drhodp(T, rho);                      /// 
            if ((zeta < 0) || (zeta > 1.0E13)) zeta = 1.0E13;
            delChi = rhobar * (zeta - delTr(rho) * Tr / T);
            y = qD * xi0 * Math.Pow(delChi / GAMMA0, nu / gam);
            if (y < 1.2E-7)
                Z = 0.0;
            else
                Z = 2.0 / (PI * y) * (((1.0 - 1.0 / k) * Math.Atan(y) + y / k) - (1.0 - Math.Exp(-1.0 / (1.0 / y + y * y / (3.0 * rhobar * rhobar)))));
            return LAMBDA * rhobar * Cpbar * T / (Constants.Tcrit * mubar) * Z;
        }
        static double Trterm(double T)
        {
            return Constants.Tcrit / T - 1.0;
        }
        static double Rhorterm(double rho)
        {
            return rho / Constants.Rhocrit - 1.0;
        }
        static double p(double T, double rho)
        {
            return rho * R * T * delta_dphi_ddelta(T, rho) * 1000;
        }

#if REGION3_ITERATE
        //
        // Newton-Raphson Technique for solving p(T,rho) for rho
        //    Solves to find root of p - rho*R*T*delta*dphi_ddelta = 0
        //    The equation is rearranged to solve for rho and turned
        //    into functions f(T,P,rho0) and f'(T,P,rho0) for the
        //    Newton-Raphson technique.  Functions for
        //    dphi/ddelta and d²phi/ddelta² were also required.  These
        //    additional Taylor functions are defined above.
        //
        double f(double T, double p, double rho0) {
            return 1.0/(rho0* rho0) - R* T* dphi_ddelta(T, rho0)/(p* Constants.Rhocrit) * 1000;
        }
        double df(double T, double p, double rho0) {
            const double rho_c = 322.0, rho_c2 = rho_c * rho_c;
            return -2.0/(rho0* rho0* rho0) - R* T* d2phi_ddelta2(T, rho0)/(p* rho_c2) * 1000;
        }
        double rhomass(double T, double p, double rho0) {
            int iter = 100;
        double f_T_p_rho0 = f(T, p, rho0);
            while (Math.Abs(f_T_p_rho0) > 1.0e-14 )
            {
                rho0 -= (f_T_p_rho0/df(T, p, rho0) );
                // don't go more than 100 iterations or throw an exception
                if (--iter == 0) throw new Exception("Failed to converge!"); 
                f_T_p_rho0 = f(T, p, rho0);
            }
            return rho0;
        }
        // END Newton-Raphson
#endif

        static double umass(double T, double rho)
        {
            return R * T * tau_dphi_dtau(T, rho);
        }
        static double smass(double T, double rho)
        {
            return R * (tau_dphi_dtau(T, rho) - phi(T, rho));
        }
        static double hmass(double T, double rho)
        {
            return R * T * (tau_dphi_dtau(T, rho) + delta_dphi_ddelta(T, rho));
        }
        static double cpmass(double T, double rho)
        {
            return R * (-tau2_d2phi_dtau2(T, rho) + Math.Pow(delta_dphi_ddelta(T, rho) - deltatau_d2phi_ddelta_dtau(T, rho), 2) / (2 * delta_dphi_ddelta(T, rho) + delta2_d2phi_ddelta2(T, rho)));
        }
        static double cvmass(double T, double rho)
        {
            return R * (-tau2_d2phi_dtau2(T, rho));
        }
        static double speed_sound(double T, double rho)
        {
            double RHS = 2 * delta_dphi_ddelta(T, rho) + delta2_d2phi_ddelta2(T, rho) - Math.Pow(delta_dphi_ddelta(T, rho) - deltatau_d2phi_ddelta_dtau(T, rho), 2) / tau2_d2phi_dtau2(T, rho);
            return Math.Sqrt(R * 1000 * T * RHS);
        }
        static double visc(double T, double rho)
        {
            /// This base region function was not inherited 
            const double mu_star = 1.0E-6; // Reference viscosity [Pa-s]
            const double mu2 = 1.0;        // For Industrial Formulation (IF97), mu2 = 1.0
            return mu_star * mu0(T) * mu1(T, rho) * mu2;
        }
        static double tcond(double T, double p, double rho)
        {
            /// This base region function was not inherited in Region3
            const double lambda_star = 0.001;
            double lambda_bar = lambda0(T) * lambda1(T, rho) + lambda2(T, p, rho);
            return lambda_star * lambda_bar;
        }
        static double drhodp(double T, double rho)
        /// Derived from IAPWS Revised Advisory Note No. 3 (See Table 2, Section 3.1 & 3.3)
        /// NOTE: rho is passed in here, not p as it is in Regions 1, 2, & 5.  This is done
        ///       because p(T,rho) is a simple algebraic in Region 3, the Helmholz functions
        ///       are functions of (T,rho), and the work has already been done by the output()
        ///       function to convert p to rho @ T.
        {
            return (rho / p(T, rho)) / (2.0 + delta2_d2phi_ddelta2(T, rho) / delta_dphi_ddelta(T, rho));
        }
        static double delTr(double rho)
        {
            /// This is the IF97 correlation for drhodp at the reducing temperature, Tr
            double rhobar = rho / Constants.Rhocrit;
            double summer = 0;
            int j;
            //
            if (rhobar <= 0.310559006) j = 0;
            else if (rhobar <= 0.776397516) j = 1;
            else if (rhobar <= 1.242236025) j = 2;
            else if (rhobar <= 1.863354037) j = 3;
            else j = 4;
            //
            double pow_rhobar = 1.0;
            for (int i = 0; i < 6; ++i)
            {
                summer += Constants.A[i, j] * pow_rhobar;
                pow_rhobar *= rhobar;
            }
            return 1.0 / summer;
        }
        static char SatSubRegionAdjust(SatState State, double p, char subregion)
        {
            switch (State)      // See if saturated state is requested
            {
                // If looking for Saturated Vapor...
                case SatState.VAPOR:
                    {     // ...force below saturation curve
                        if (subregion == 'C') return 'T';
                        else if (subregion == 'S')
                        {
                            if (p < 20.5)
                                return 'T';
                            else
                                return 'R';
                        }
                        else if (subregion == 'U')
                        {
                            if (p < 21.90096265)
                                return 'X';
                            else
                                return 'Z';
                        }
                        else if (subregion == 'Y') return 'Z';
                        break;
                    }

                // If looking for Saturated Liquid...
                case SatState.LIQUID:
                    {     // ...force above saturation curve
                        if (subregion == 'Z')
                        {
                            if (p > 21.93161551)
                                return 'Y';
                            else
                                return 'U';
                        }
                        else if (subregion == 'X') return 'U';
                        else if ((subregion == 'R') || (subregion == 'K')) return 'S';
                        else if (subregion == 'T')
                        {
                            if (p > 19.00881189173929)
                                return 'S';
                            else
                                return 'C';
                        }
                        break;
                    }
                case SatState.NONE:
                default: return subregion;
            }
            return subregion;  // in case no adjustment needs to be made
        }  // SatSubRegionAdjust

        public static double output(IF97Parameters key, double T, double p, SatState State)
        {
            double rho;
            char region = BackwardsRegion3.RegionDetermination(T, p);

            // if this is a saturated vapor or liquid function, make sure we're on
            // the correct side of the saturation curve and adjust region before
            // calculating density.
            region = SatSubRegionAdjust(State, p, region);

            rho = 1 / BackwardsRegion3.Region3_v_TP(region, T, p);

#if REGION3_ITERATE
            // Use previous rho value from algebraic equations 
            //      as an initial guess to solve rhomass iteratively 
            //      with Newton-Raphson
            rho = rhomass(T, p, rho);
#endif
            switch (key)                 // return all properties using the new rho value
            {
                case IF97Parameters.d: return rho;
                case IF97Parameters.h: return hmass(T, rho);
                case IF97Parameters.s: return smass(T, rho);
                case IF97Parameters.u: return umass(T, rho);
                case IF97Parameters.cp: return cpmass(T, rho);
                case IF97Parameters.cv: return cvmass(T, rho);
                case IF97Parameters.w: return speed_sound(T, rho);
                case IF97Parameters.mu: return visc(T, rho);
                case IF97Parameters.k: return tcond(T, p, rho);
                case IF97Parameters.drhodp: return drhodp(T, rho);

                default:
                    throw new ArgumentException("Bad key to output");  // JPH: changed this to invalid_argument exception
            }
        }
    }
}
