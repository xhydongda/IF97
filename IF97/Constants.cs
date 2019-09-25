using System;

namespace IF97
{
    public sealed class Constants
    {
        private Constants() { }

        public const double Tcrit = 647.096;               // K
        public const double Pcrit = 22.064;                // Pa
        public const double Rhocrit = 322.0;               // kg/m³
        public const double Scrit = 4.41202148223476;      // J/kg-K (needed for backward eqn. in Region 3(a)(b)
        public const double Ttrip = 273.16;                // K
        public const double Ptrip = 0.000611656;           // Pa
        public const double Tmin = 273.15;                 // K
        public const double Tmax = 1073.15;                // K
        public const double Pmin = 0.000611213;            // Pa
        public const double Pmax = 100.0;                  // Pa
        public const double Rgas = 0.461526;               // J/kg-K : mass based!
        public const double MW = 0.018015268;              // kg/mol
                                                           // Bounds for Region Determination
        public const double Text = 2273.15;                // Extended (Region 5) Temperature Limit (Region 5) [K]
        public const double Pext = 50.0;                   // Extended (Region 5) Pressure Limit (Region 5) [Pa]
        public const double P23min = 16.529164252605;      // Min Pressure on Region23 boundary curve; Max is Pmax
        public const double T23min = 623.15;               // Min Temperature on Region23 boundary curve
        public const double T23max = 863.15;               // Max Temperature on Region23 boundary curve
        public const double P2bcmin = 6.54670;             // Min Pressure [MPa] on H2b2c boundary curve; Max is Pmax
        public const double S2bc = 5.85;                   // Min Pressure [MPa] on H2b2c boundary curve; Max is Pmax
                                                           // Bounds for Backward p(h,s), t(h,s) Determination
        public const double Smin = 0.0;                    // Min Entropy [kJ/kg-K] for Backward p(h,s)
        public const double Smax = 11.921054825051103;     // Max Entropy [kJ/kg-K] for Backward p(h,s)
        public const double STPmax = 6.04048367171238;     // S(Tmax,Pmax)
        public const double Sgtrip = 9.155492076509681;    // Sat. Vapor  Entropy [kJ/kg-K] at Triple Point
        public const double Sftrip = -4.09187776773977E-7; // Sat. Liquid Entropy [kJ/kg-K] at Triple Point
        public const double Hgtrip = 2500.9109532932;      // Sat. Vapor  Enthalpy [kJ/kg] at Triple Point
        public const double Hftrip = 5.16837786577998E-4;  // Sat. Liquid Enthalpy [kJ/kg] at Triple Point
        public const double SfT23 = 3.778281340;           // Sat. Liquid Entropy [KJ/kg-K] at T23min
        public const double SgT23 = 5.210887825;           // Sat. Vapor  Entropy [KJ/kg-K] at T23min
        public const double S13min = 3.397782955;          // Entropy at (T13,Pmax)
        public const double S23min = 5.048096828;          // B23 Bounding Box
        public const double S23max = 5.260578707;          // B23 Bounding Box
        public const double H23min = 2.563592004E3;        // B23 Bounding Box
        public const double H23max = 2.812942061E3;        // B23 Bounding Box
        public static readonly ValueTuple<int, int, double>[] Hresiddata = {        // Residual H for viscosity
            (0, 0,  5.20094e-1),
            (1, 0,  8.50895e-2),
            (2, 0, -1.08374   ),
            (3, 0, -2.89555e-1),
            (0, 1,  2.22531e-1),
            (1, 1,  9.99115e-1),
            (2, 1,  1.88797   ),
            (3, 1,  1.26613   ),
            (5, 1,  1.20573e-1),
            (0, 2, -2.81378e-1),
            (1, 2, -9.06851e-1),
            (2, 2, -7.72479e-1),
            (3, 2, -4.89837e-1),
            (4, 2, -2.57040e-1),
            (0, 3,  1.61913e-1),
            (1, 3,  2.57399e-1),
            (0, 4, -3.25372e-2),
            (3, 4,  6.98452e-2),
            (4, 5,  8.72102e-3),
            (3, 6, -4.35673e-3),
            (5, 6, -5.93264e-4)
        };
        public static readonly ValueTuple<int, double>[] Hidealdata = {         // Ideal H for viscosity
            (0,  1.67752   ),
            (1,  2.20462   ),
            (2,  0.6366564 ),
            (3, -0.241605  )
        };
        public static readonly ValueTuple<int, int, double>[] Lresiddata = {       // Residual L for Thermal Conductivity
            ( 0, 0,  1.60397357000 ),
            ( 1, 0,  2.33771842000 ),
            ( 2, 0,  2.19650529000 ),
            ( 3, 0, -1.21051378000 ),
            ( 4, 0, -2.72033700000 ),
            ( 0, 1, -0.64601352300 ),
            ( 1, 1, -2.78843778000 ),
            ( 2, 1, -4.54580785000 ),
            ( 3, 1,  1.60812989000 ),
            ( 4, 1,  4.57586331000 ),
            ( 0, 2,  0.11144390600 ),
            ( 1, 2,  1.53616167000 ),
            ( 2, 2,  3.55777244000 ),
            ( 3, 2, -0.62117814100 ),
            ( 4, 2, -3.18369245000 ),
            ( 0, 3,  0.10299735700 ),
            ( 1, 3, -0.46304551200 ),
            ( 2, 3, -1.40944978000 ),
            ( 3, 3,  0.07163732240 ),
            ( 4, 3,  1.11683480000 ),
            ( 0, 4, -0.05041236340 ),
            ( 1, 4,  0.08328270190 ),
            ( 2, 4,  0.27541827800 ),
            ( 3, 4,  0.00000000000 ),
            ( 4, 4, -0.19268305000 ),
            ( 0, 5,  0.00609859258 ),
            ( 1, 5, -0.00719201245 ),
            ( 2, 5, -0.02059388160 ),
            ( 3, 5,  0.00000000000 ),
            ( 4, 5,  0.01291384200 )
        };
        public static readonly ValueTuple<int, double>[] Lidealdata = {          // Ideal L for thermal conductivity
            (0,  2.443221E-3),
            (1,  1.323095E-2),
            (2,  6.770357E-3),
            (3, -3.454586E-3),
            (4,  4.096266E-4)
        };

        public static readonly double[,] A = {
            {  6.53786807199516,  6.52717759281799,   5.35500529896124,   1.55225959906681,   1.11999926419994  },
            { -5.61149954923348, -6.30816983387575,  -3.96415689925446,   0.464621290821181,  0.595748562571649 },
            {  3.39624167361325,  8.08379285492595,   8.91990208918795,   8.93237374861479,   9.88952565078920  },
            { -2.27492629730878, -9.82240510197603, -12.03387295057900, -11.03219600611260, -10.32550511470400  },
            { 10.26318546627090, 12.13584137913950,   9.19494865194302,   6.16780999933360,   4.66861294457414  },
            {  1.97815050331519, -5.54349664571295,  -2.16866274479712,  -0.965458722086812, -0.503243546373828 },
        };
    }
}
