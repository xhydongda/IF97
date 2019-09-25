using System;

namespace IF97
{
    public class Region1 : RegionBase
    {
        // Note: the coefficients of n_i have been multiplied by -1**I_i such that all Gibbs terms are of the form (PI-7.1)**(I_i) rather than (7.1-PI)**(I_i)
        static readonly ValueTuple<int, int, double>[] Region1residdata = {
            (0, -2, 0.14632971213167),
            (0, -1, -0.84548187169114),
            (0, 0, -3.756360367204),
            (0, 1, 3.3855169168385),
            (0, 2, -0.95791963387872),
            (0, 3, 0.15772038513228),
            (0, 4, -0.016616417199501),
            (0, 5, 0.00081214629983568),
            (1, -9, -0.00028319080123804),
            (1, -7, 0.00060706301565874),
            (1, -1, 0.018990068218419),
            (1, 0, 0.032529748770505),
            (1, 1, 0.021841717175414),
            (1, 3, 0.00005283835796993),
            (2, -3, -0.00047184321073267),
            (2, 0, -0.00030001780793026),
            (2, 1, 0.000047661393906987),
            (2, 3, -4.4141845330846E-06),
            (2, 17, -7.2694996297594E-16),
            (3, -4, 0.000031679644845054),
            (3, 0, 2.8270797985312E-06),
            (3, 6, 8.5205128120103E-10),
            (4, -5, -0.0000022425281908),
            (4, -2, -6.5171222895601E-07),
            (4, 10, -1.4341729937924E-13),
            (5, -8, 4.0516996860117E-07),
            (8, -11, -1.2734301741641E-09),
            (8, -6, -1.7424871230634E-10),
            (21, -29, 6.8762131295531E-19),
            (23, -31, -1.4478307828521E-20),
            (29, -38, -2.6335781662795E-23),
            (30, -39, -1.1947622640071E-23),
            (31, -40, -1.8228094581404E-24),
            (32, -41, -9.3537087292458E-26)
        };
        static readonly ValueTuple<int, double>[] Region1idealdata = { };
        public Region1() : base(Region1residdata, Region1idealdata)
        {
            T_star = 1386;
            p_star = 16.53;
        }

        protected override double speed_sound(double T, double p)
        {
            // Evidently this formulation is special for some reason, and cannot be implemented using the base class formulation
            // see Table 3
            double tau = T_star / T;
            double RHS = Math.Pow(dgammar_dPI(T, p), 2) / (Math.Pow(dgammar_dPI(T, p) - tau * d2gammar_dPIdTAU(T, p), 2) / (tau * tau * d2gammar_dTAU2(T, p)) - d2gammar_dPI2(T, p));
            return Math.Sqrt(R * 1000 * T * RHS);
        }

        protected override double cvmass(double T, double p)
        {
            // Evidently this formulation is special for some reason, and cannot be implemented using the base class formulation
            // see Table 3
            double tau = T_star / T;
            return R * (-tau * tau * d2gammar_dTAU2(T, p) + Math.Pow(dgammar_dPI(T, p) - tau * d2gammar_dPIdTAU(T, p), 2) / d2gammar_dPI2(T, p));
        }
        protected override double drhodp(double T, double p)
        {
            //double PI = p/p_star;
            /// This one is different as well...
            /// Derived from IAPWS Revised Advisory Note No. 3 (See Table 2, Section 4.1 & 4.2)
            return -d2gammar_dPI2(T, p) / (Math.Pow(dgammar_dPI(T, p), 2) * R * T) * 1000;
        }
        protected override double TAUrterm(double T)
        {
            return T_star / T - 1.222;
        }
        protected override double PIrterm(double p)
        {
            return p / p_star - 7.1;
        }
        protected override double TAU0term(double _)
        {
            return 0.0;
        }
    }
}
