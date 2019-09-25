using System;

namespace IF97
{
    public abstract class RegionBase
    {
        int[] Ir, Jr;
        double[] nr;
        int[] J0;
        double[] n0;
        protected double T_star, p_star;
        protected const double R = Constants.Rgas;
        /// For Viscosity Calculations
        int[] muJ0;
        double[] mun0;
        int[] muIr, muJr;
        double[] munr;
        /// For Thermal Conductivity Calculations
        int[] lamJ0;
        double[] lamn0;
        int[] lamIr, lamJr;
        double[] lamnr;

        public RegionBase(ValueTuple<int, int, double>[] resid, ValueTuple<int, double>[] ideal)
        {
            int len = resid.Length;
            if (len > 0)
            {
                Ir = new int[len];
                Jr = new int[len];
                nr = new double[len];
                for (int i = 0; i < len; i++)
                {
                    var tuple = resid[i];
                    Ir[i] = tuple.Item1;
                    Jr[i] = tuple.Item2;
                    nr[i] = tuple.Item3;
                }
            }
            len = ideal.Length;
            n0 = new double[len];
            J0 = new int[len];
            for (int i = 0; i < len; i++)
            {
                var tuple = ideal[i];
                J0[i] = tuple.Item1;
                n0[i] = tuple.Item2;
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
        }

        double rhomass(double T, double p)
        {
            return p_star / (R * T) * 1000.0 / (dgammar0_dPI(T, p) + dgammar_dPI(T, p));
        }
        double hmass(double T, double p)
        {
            return R * T_star * (dgammar0_dTAU(T, p) + dgammar_dTAU(T, p));
        }
        double smass(double T, double p)
        {
            double tau = T_star / T;
            return R * (tau * (dgammar0_dTAU(T, p) + dgammar_dTAU(T, p)) - (gammar(T, p) + gammar0(T, p)));
        }
        double umass(double T, double p)
        {
            double tau = T_star / T, PI = p / p_star;
            return R * T * (tau * (dgammar0_dTAU(T, p) + dgammar_dTAU(T, p)) - PI * (dgammar0_dPI(T, p) + dgammar_dPI(T, p)));
        }
        double cpmass(double T, double p)
        {
            double tau = T_star / T;
            return -R * tau * tau * (d2gammar_dTAU2(T, p) + d2gammar0_dTAU2(T, p));
        }
        protected virtual double cvmass(double T, double p)
        {
            double tau = T_star / T, PI = p / p_star;
            return cpmass(T, p) - R * FastPow.Pow(1 + PI * dgammar_dPI(T, p) - tau * PI * d2gammar_dPIdTAU(T, p), 2) / (1 - PI * PI * d2gammar_dPI2(T, p));
        }
        protected virtual double speed_sound(double T, double p)
        {
            double tau = T_star / T, PI = p / p_star;
            double RHS = (1 + 2 * PI * dgammar_dPI(T, p) + PI * PI * FastPow.Pow(dgammar_dPI(T, p), 2)) / ((1 - PI * PI * d2gammar_dPI2(T, p)) + FastPow.Pow(1 + PI * dgammar_dPI(T, p) - tau * PI * d2gammar_dPIdTAU(T, p), 2) / (tau * tau * (d2gammar0_dTAU2(T, p) + d2gammar_dTAU2(T, p))));
            return Math.Sqrt(R * 1000 * T * RHS);
        }
        public double visc(double T, double rho)
        {
            /// This base region function is valid for all IF97 regions since it is a function
            /// of density, not pressure, and can be called from any region instance.
            const double mu_star = 1.0E-6; // Reference viscosity [Pa-s]
            const double mu2 = 1.0;        // For Industrial Formulation (IF97), mu2 = 1.0
            return mu_star * mu0(T) * mu1(T, rho) * mu2;
        }
        double tcond(double T, double p, double rho)
        {
            /// This base region function is valid for all IF97 regions 
            const double lambda_star = 0.001;  // Reference conductivity [W/m-K]
            double lambda_bar = lambda0(T) * lambda1(T, rho) + lambda2(T, p, rho);
            return lambda_star * lambda_bar;
        }
        protected virtual double drhodp(double T, double p)
        {
            /// Only valid for regions 2 and 5.  Will be overridden in Regions 1 and 3.
            /// Derived from IAPWS Revised Advisory Note No. 3 (See Table 2, Section 4.1 & 4.3)
            double PI = p / p_star;
            return (rhomass(T, p) / p) * ((1.0 - PI * PI * d2gammar_dPI2(T, p)) / (1.0 + PI * dgammar_dPI(T, p)));
        }
        double delTr(double rho)
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
            for (int i = 0; i < 6; i++)
                summer += Constants.A[i, j] * FastPow.Pow(rhobar, i);
            return 1.0 / summer;
        }
        protected abstract double PIrterm(double p);
        protected abstract double TAUrterm(double T);
        protected abstract double TAU0term(double T);
        public double output(IF97Parameters key, double T, double p)
        {
            switch (key)
            {
                case IF97Parameters.T: return T;
                case IF97Parameters.p: return p;
                case IF97Parameters.d: return rhomass(T, p);
                case IF97Parameters.h: return hmass(T, p);
                case IF97Parameters.s: return smass(T, p);
                case IF97Parameters.u: return umass(T, p);
                case IF97Parameters.cp: return cpmass(T, p);
                case IF97Parameters.cv: return cvmass(T, p);
                case IF97Parameters.w: return speed_sound(T, p);
                case IF97Parameters.mu: return visc(T, rhomass(T, p));   // Viscosity is a function of rho.
                case IF97Parameters.tc: return tcond(T, p, rhomass(T, p)); // Conductivity needs p and rho.
                case IF97Parameters.drhodp: return drhodp(T, p);       // For verification testing.
            }
            throw new ArgumentOutOfRangeException("Unable to match input parameters");
        }

        protected double gammar(double T, double p)
        {
            double _PI = PIrterm(p), _TAU = TAUrterm(T);
            double summer = 0;
            for (int i = 0; i < Jr.Length; ++i)
            {
                summer += nr[i] * FastPow.Pow(_PI, Ir[i]) * FastPow.Pow(_TAU, Jr[i]);
            }
            return summer;
        }
        protected double dgammar_dPI(double T, double p)
        {
            double _PI = PIrterm(p), _TAU = TAUrterm(T);
            double summer = 0;
            for (int i = 0; i < Jr.Length; ++i)
            {
                summer += nr[i] * Ir[i] * FastPow.Pow(_PI, Ir[i] - 1) * FastPow.Pow(_TAU, Jr[i]);
            }
            return summer;
        }
        protected double d2gammar_dPI2(double T, double p)
        {
            double _PI = PIrterm(p), _TAU = TAUrterm(T);
            double summer = 0;
            for (int i = 0; i < Jr.Length; ++i)
            {
                summer += nr[i] * Ir[i] * (Ir[i] - 1) * FastPow.Pow(_PI, Ir[i] - 2) * FastPow.Pow(_TAU, Jr[i]);
            }
            return summer;
        }
        protected double dgammar_dTAU(double T, double p)
        {
            double _PI = PIrterm(p), _TAU = TAUrterm(T);
            double summer = 0;
            for (int i = 0; i < Jr.Length; ++i)
            {
                summer += nr[i] * Jr[i] * FastPow.Pow(_PI, Ir[i]) * FastPow.Pow(_TAU, Jr[i] - 1);
            }
            return summer;
        }
        protected double d2gammar_dPIdTAU(double T, double p)
        {
            double _PI = PIrterm(p), _TAU = TAUrterm(T);
            double summer = 0;
            for (int i = 0; i < Jr.Length; ++i)
            {
                summer += nr[i] * Jr[i] * Ir[i] * FastPow.Pow(_PI, Ir[i] - 1) * FastPow.Pow(_TAU, Jr[i] - 1);
            }
            return summer;
        }
        protected double d2gammar_dTAU2(double T, double p)
        {
            double _PI = PIrterm(p), _TAU = TAUrterm(T);
            double summer = 0;
            for (int i = 0; i < Jr.Length; ++i)
            {
                summer += nr[i] * Jr[i] * (Jr[i] - 1) * FastPow.Pow(_PI, Ir[i]) * FastPow.Pow(_TAU, Jr[i] - 2);
            }
            return summer;
        }
        double gammar0(double T, double p)
        {
            if (J0.Length == 0) { return 0; } // Region 1 has no term
            double PI = p / p_star, _TAU = TAU0term(T);
            double summer = Math.Log(PI);
            for (int i = 0; i < n0.Length; ++i)
            {
                summer += n0[i] * FastPow.Pow(_TAU, J0[i]);
            }
            return summer;
        }
        double dgammar0_dPI(double _, double p)
        {
            if (J0.Length == 0) { return 0; } // Region 1 has no term
            double PI = p / p_star;
            return 1.0 / PI;
        }
        double d2gammar0_dPI2(double _, double p)
        {
            if (J0.Length == 0) { return 0; } // Region 1 has no term
            double PI = p / p_star;
            return -1.0 / (PI * PI);
        }
        double dgammar0_dTAU(double T, double _)
        {
            double _TAU = TAU0term(T);
            double summer = 0;
            for (int i = 0; i < J0.Length; ++i)
            {
                summer += n0[i] * J0[i] * FastPow.Pow(_TAU, J0[i] - 1);
            }
            return summer;
        }
        double d2gammar0_dTAU2(double T, double _)
        {
            double _TAU = TAU0term(T);
            double summer = 0;
            for (int i = 0; i < J0.Length; ++i)
            {
                summer += n0[i] * J0[i] * (J0[i] - 1) * FastPow.Pow(_TAU, J0[i] - 2);
            }
            return summer;
        }
        double mu0(double T)
        {
            double T_bar = T / Constants.Tcrit;
            double summer = 0.0;
            for (int i = 0; i < muJ0.Length; ++i)
            {
                summer += mun0[i] / FastPow.Pow(T_bar, muJ0[i]);
            }
            return 100.0 * Math.Sqrt(T_bar) / summer;
        }
        double mu1(double T, double rho)
        {
            double rho_bar = rho / Constants.Rhocrit;
            double summer = 0.0;
            for (int i = 0; i < muJr.Length; ++i)
            {
                summer += rho_bar * FastPow.Pow(Trterm(T), muIr[i]) * munr[i] * FastPow.Pow(Rhorterm(rho), muJr[i]);
            }
            return Math.Exp(summer);
        }
        double lambda0(double T)
        {
            double T_bar = T / Constants.Tcrit;
            double summer = 0.0;
            for (int i = 0; i < lamJ0.Length; ++i)
            {
                summer += lamn0[i] / FastPow.Pow(T_bar, lamJ0[i]);
            }
            return Math.Sqrt(T_bar) / summer;
        }
        double lambda1(double T, double rho)
        {
            double rho_bar = rho / Constants.Rhocrit;
            double summer = 0.0;
            for (int i = 0; i < lamJr.Length; ++i)
            {
                summer += rho_bar * FastPow.Pow(Trterm(T), lamIr[i]) * lamnr[i] * FastPow.Pow(Rhorterm(rho), lamJr[i]);
            }
            return Math.Exp(summer);
        }
        protected virtual double lambda2(double T, double p, double rho)
        {
            double y, Cpbar, mubar, k, Z, delChi;
            double rhobar = rho / Constants.Rhocrit;
            const double LAMBDA = 177.8514;
            const double qD = 1.0 / 0.40;
            double Tr = 1.5 * Constants.Tcrit;
            const double xi0 = 0.13;
            const double nu = 0.630;
            const double gam = 1.239;
            const double GAMMA0 = 0.06;
            const double PI = 3.141592654;
            double Cpstar = 0.46151805;  /// Note: Slightly lower than IF97 Rgas
            Cpbar = cpmass(T, p) / Cpstar;
            if ((Cpbar < 0) || (Cpbar > 1.0E13)) Cpbar = 1.0E13;     /// Unit-less
            k = cpmass(T, p) / cvmass(T, p);
            mubar = visc(T, rho) / 1.0E-6;
            delChi = rhobar * (Constants.Pcrit / Constants.Rhocrit * drhodp(T, p) - delTr(rho) * Tr / T);
            if (delChi > 0)                            /// At low (T,p), delChi can go negative, causing
                y = qD * xi0 * Math.Pow(delChi / GAMMA0, nu / gam);  ///   y to be imaginary from this nth-root equation.
            else                                       ///   
                y = 0.0;                               ///   Limit delChi to > 0, values.
            if (y < 1.2E-7)                            /// Z is not calculated if y < 1.2E-7 since the
                Z = 0.0;                               ///   critical enhancement becomes insignificant.
            else
                Z = 2.0 / PI / y * (((1.0 - 1.0 / k) * Math.Atan(y) + y / k) - (1.0 - Math.Exp(-1.0 / (1.0 / y + y * y / (3.0 * rhobar * rhobar)))));
            return LAMBDA * rhobar * Cpbar * T / (Constants.Tcrit * mubar) * Z;
        }
        double Trterm(double T)
        {
            return Constants.Tcrit / T - 1.0;
        }
        double Rhorterm(double rho)
        {
            return rho / Constants.Rhocrit - 1.0;
        }
    }
}
