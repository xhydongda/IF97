using System;

namespace IF97
{
    public class Region5 : RegionBase
    {
        static readonly ValueTuple<int, int, double>[] Region5residdata = {
        (1, 1,  0.15736404855259e-2),
        (1, 2,  0.90153761673944e-3),
        (1, 3, -0.50270077677648e-2),
        (2, 3,  0.22440037409485e-5),
        (2, 9, -0.41163275453471e-5),
        (3, 7,  0.37919454822955e-7)
        };
        static readonly ValueTuple<int, double>[] Region5idealdata = {
        ( 0, -0.13179983674201e2),
        ( 1,  0.68540841634434e1),
        (-3, -0.24805148933466e-1),
        (-2,  0.36901534980333),
        (-1, -0.31161318213925e1),
        ( 2, -0.32961626538917)
        };
        public Region5() : base(Region5residdata, Region5idealdata)
        {
            T_star = 1000; p_star = 1;
        }
        protected override double TAUrterm(double T)
        {
            return T_star / T;
        }
        protected override double PIrterm(double p)
        {
            return p / p_star;
        }
        protected override double TAU0term(double T)
        {
            return T_star / T;
        }
        protected override double lambda2(double T, double p, double rho)
        {
            return 0;
        }
    }
}
