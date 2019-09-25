using System;

namespace IF97
{
    /// This "region" is the saturation curve
    public class Region4
    {
        static double[] n;
        static double p_star, T_star;
        static readonly ValueTuple<int, double>[] sat = {
            (1,  0.11670521452767e4),
            (2, -0.72421316703206e6),
            (3, -0.17073846940092e2),
            (4,  0.12020824702470e5),
            (5, -0.32325550322333e7),
            (6,  0.14915108613530e2),
            (7, -0.48232657361591e4),
            (8,  0.40511340542057e6),
            (9, -0.23855557567849),
            (10, 0.65017534844798e3)
        };
        static Region4()
        {
            p_star = 1;
            T_star = 1;
            n = new double[sat.Length + 1];
            n[0] = 0;
            for (int i = 1; i < n.Length; ++i)
            {
                n[i] = sat[i - 1].Item2;
            }
        }
        public static double p_T(double T)
        {
            // Allow extrapolation down to Pmin = P(Tmin=273.15K) = 611.213 Pa
            if ((T < Constants.Tmin) || (T > Constants.Tcrit))
            {
                throw new ArgumentOutOfRangeException("Temperature out of range");
            }
            double theta = T / T_star + n[9] / (T / T_star - n[10]);
            double A = theta * theta + n[1] * theta + n[2];
            double B = n[3] * theta * theta + n[4] * theta + n[5];
            double C = n[6] * theta * theta + n[7] * theta + n[8];
            return p_star * FastPow.Pow(2 * C / (-B + Math.Sqrt(B * B - 4 * A * C)), 4);
        }
        public static double T_p(double p)
        {
            // Allow extrapolation down to Pmin = P(Tmin=273.15K) = 611.213 Pa
            if ((p < Constants.Pmin) || (p > Constants.Pcrit))
            {
                throw new ArgumentOutOfRangeException("Pressure out of range");
            }


            double beta2 = Math.Sqrt(p / p_star);
            double beta = Math.Sqrt(beta2);
            double[] EFG = new double[3];
            // Each cycle can be vectorized
            EFG[0] = 1.0; EFG[1] = n[1]; EFG[2] = n[2];
            for (int i = 0; i < 3; ++i)
            {
                EFG[i] *= beta;
            }
            for (int i = 0; i < 3; ++i)
            {
                EFG[i] += n[i + 3];
            }
            for (int i = 0; i < 3; ++i)
            {
                EFG[i] *= beta;
            }
            for (int i = 0; i < 3; ++i)
            {
                EFG[i] += n[i + 6];
            }
            double E = EFG[0], F = EFG[1], G = EFG[2];
            double D = 2 * G / (-F - Math.Sqrt(F * F - 4 * E * G));
            double n10pD = n[10] + D;
            return T_star * 0.5 * (n10pD - Math.Sqrt(n10pD * n10pD - 4 * (n[9] + n[10] * D)));
        }

        public static double sigma_t(double T)
        {
            // Surface Tension [mN/m] in two-phase region as a function of temperature [K]
            // Implemented from IAPWS R1-76(2014).
            // May be extrapolated down to -25C in the super-cooled region.
            if ((T < (Constants.Ttrip - 25.0)) || (T > Constants.Tcrit))
            {
                throw new ArgumentOutOfRangeException("Temperature out of range");
            }
            double Tau = 1.0 - T / Constants.Tcrit;
            const double B = 235.8 / 1000;  // Published value in [mN/m]; Convert to SI [N/m] in all cases 
            const double b = -0.625;
            const double mu = 1.256;
            return B * Math.Pow(Tau, mu) * (1.0 + b * Tau);
        }

        // ******************************************************************************** //
        //                               2-Phase Functions                                  //
        // ******************************************************************************** //
        /// Get the saturation temperature [K] as a function of p [Pa]
        public static double Tsat97(double p)
        {
            return T_p(p);
        }
        /// Get the saturation pressure [Pa] as a function of T [K]
        public static double psat97(double T)
        {
            return p_T(T);
        }
        /// Get surface tension [N/m] as a function of T [K]
        public static double sigma97(double T)
        {
            return sigma_t(T);
        }
    }
}
