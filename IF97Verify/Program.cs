﻿using System;
using System.Diagnostics;

namespace IF97
{
    class Program
    {
        static void Main(string[] args)
        {
            verifyIF97();
        }

        static void testpow(int start, int end)
        {
            double x = 2;
            for (int j = start; j < end; j++)
            {
                double y1 = FastPow.Pow(x, j);
                double y2 = Math.Pow(x, j);
                Console.WriteLine(String.Format("{0}:{1}", j, y1 - y2));
            }
            Console.ReadLine();
        }

        static void testpowperformance(int start, int end, int times = 10000000)
        {
            Random r = new Random();
            Stopwatch sw = new Stopwatch();
            sw.Start();
            for (int i = 0; i < times; i++)
            {
                for (int j = start; j < end; j++)
                {
                    double y1 = FastPow.Pow(r.NextDouble(), j);
                }

            }
            sw.Stop();
            Console.WriteLine(sw.ElapsedTicks);
            sw.Reset();
            sw.Start();
            for (int i = 0; i < times; i++)
            {
                for (int j = start; j < end; j++)
                {
                    double y2 = Math.Pow(r.NextDouble(), j);
                }
            }
            sw.Stop();
            Console.WriteLine(sw.ElapsedTicks);
            Console.ReadLine();
        }

        static void verifyIF97()
        {
            Console.Write("*****************************************************************\n");
            Console.Write("******************** Verification Tables ************************\n");
            Console.Write("* Tables below are printed for verification.  Unless otherwise  *\n");
            Console.Write("* noted, tables are reproduced from the                         *\n");
            Console.Write("* \"Revised Release on the IAPWS Industrial Formulation 1997\"    *\n");
            Console.Write("* IAPWS R7-97(2012).                                            *\n");
            Console.Write("*****************************************************************\n\n\n");


            double T1 = 300, T2 = 300, T3 = 500, p1 = 3, p2 = 80, p3 = 3;
            Console.Write("*****************************************************************\n");
            Console.Write("******************** Table 5 - Region 1 *************************\n");
            Console.Write("*****************************************************************\n");
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "v ", Api.vmass_Tp(T1, p1), Api.vmass_Tp(T2, p2), Api.vmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "h ", Api.hmass_Tp(T1, p1), Api.hmass_Tp(T2, p2), Api.hmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "u ", Api.umass_Tp(T1, p1), Api.umass_Tp(T2, p2), Api.umass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "s ", Api.smass_Tp(T1, p1), Api.smass_Tp(T2, p2), Api.smass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "cp", Api.cpmass_Tp(T1, p1), Api.cpmass_Tp(T2, p2), Api.cpmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "w ", Api.speed_sound_Tp(T1, p1), Api.speed_sound_Tp(T2, p2), Api.speed_sound_Tp(T3, p3)));
            Console.Write("***************************************************************\n\n");

            T1 = 300; T2 = 700; T3 = 700; p1 = 0.0035; p2 = 0.0035; p3 = 30;
            Console.Write("***************************************************************\n");
            Console.Write("******************* Table 15 - Region 2 ***********************\n");
            Console.Write("***************************************************************\n");
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "v ", Api.vmass_Tp(T1, p1), Api.vmass_Tp(T2, p2), Api.vmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "h ", Api.hmass_Tp(T1, p1), Api.hmass_Tp(T2, p2), Api.hmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "u ", Api.umass_Tp(T1, p1), Api.umass_Tp(T2, p2), Api.umass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "s ", Api.smass_Tp(T1, p1), Api.smass_Tp(T2, p2), Api.smass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "cp", Api.cpmass_Tp(T1, p1), Api.cpmass_Tp(T2, p2), Api.cpmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "w ", Api.speed_sound_Tp(T1, p1), Api.speed_sound_Tp(T2, p2), Api.speed_sound_Tp(T3, p3)));
            Console.Write("***************************************************************\n\n");

            T1 = 650; T2 = 650; T3 = 750; p1 = 25.5837018; p2 = 22.2930643; p3 = 78.3095639;
            Console.Write("***************************************************************\n");
            Console.Write("******************* Table 33* - Region 3 **********************\n");
            Console.Write("***************************************************************\n");
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "rho", Api.rhomass_Tp(T1, p1), Api.rhomass_Tp(T2, p2), Api.rhomass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "h  ", Api.hmass_Tp(T1, p1), Api.hmass_Tp(T2, p2), Api.hmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "s  ", Api.smass_Tp(T1, p1), Api.smass_Tp(T2, p2), Api.smass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "cp ", Api.cpmass_Tp(T1, p1), Api.cpmass_Tp(T2, p2), Api.cpmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "w  ", Api.speed_sound_Tp(T1, p1), Api.speed_sound_Tp(T2, p2), Api.speed_sound_Tp(T3, p3)));
            Console.Write("***************************************************************\n");
            Console.Write("* NOTE: This table is evaluated by first evaluating density    \n");
            Console.Write("        from the three pressure values of Table 33 in the      \n");
            Console.Write("        2007 Revised Release document, using the reverse       \n");
            Console.Write("        formulas from the (2014) \"Revised Supplementary       \n");
            Console.Write("        Release on Backward Equations for Specific Volume\".   \n");
            Console.Write("        As a result, the values below will not be exactly      \n");
            Console.Write("        the Table 33 values, but should be within +/-1.0E-6.  \n");
            Console.Write("        of the published values.                               \n\n");

            T1 = 1500; T2 = 1500; T3 = 2000; p1 = 0.5; p2 = 30; p3 = 30;
            Console.Write("***************************************************************\n");
            Console.Write("******************* Table 42 - Region 5 ***********************\n");
            Console.Write("***************************************************************\n");
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "v ", Api.vmass_Tp(T1, p1), Api.vmass_Tp(T2, p2), Api.vmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "h ", Api.hmass_Tp(T1, p1), Api.hmass_Tp(T2, p2), Api.hmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "s ", Api.smass_Tp(T1, p1), Api.smass_Tp(T2, p2), Api.smass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "cp", Api.cpmass_Tp(T1, p1), Api.cpmass_Tp(T2, p2), Api.cpmass_Tp(T3, p3)));
            Console.Write(String.Format("{0} {1:E8} {2:E8} {3:E8}\n", "w ", Api.speed_sound_Tp(T1, p1), Api.speed_sound_Tp(T2, p2), Api.speed_sound_Tp(T3, p3)));
            Console.Write("***************************************************************\n\n\n");

            Console.Write("**************** Reverse Functions T(p,h) & T(p,s) ********************\n\n");

            p1 = 3.0; p2 = 80; p3 = 80; double h1 = 500, h2 = 500, h3 = 1500;
            double s1 = 0.5, s2 = 0.5, s3 = 3;
            Console.Write("_______________________________________________________________________\n");
            Console.Write("                                Region 1                               \n");
            Console.Write("_______________________________________________________________________\n");
            Console.Write("              Table 7              |               Table 9             \n");
            Console.Write("___________________________________|___________________________________\n");
            Console.Write(" p/MPa  h/(kJ/kg)        T/K       | p/MPa s/(kJ/kg-K)       T/K       \n");
            Console.Write("___________________________________|___________________________________\n");
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F}  {4:F1}        {5:E8}\n", p1, h1, Api.T_phmass(p1, h1),
                                                                                    p1, s1, Api.T_psmass(p1, s1)));
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F}  {4:F1}        {5:E8}\n", p2, h2, Api.T_phmass(p2, h2),
                                                                                    p2, s2, Api.T_psmass(p2, s2)));
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F}  {4:F1}        {5:E8}\n", p3, h3, Api.T_phmass(p3, h3),
                                                                                    p3, s3, Api.T_psmass(p3, s3)));
            Console.Write("_______________________________________________________________________\n\n\n");

            p1 = 0.001; p2 = 3; p3 = 3; h1 = 3000; h2 = 3000; h3 = 4000;
            double p1s = 0.1, p2s = 0.1, p3s = 2.5; s1 = 7.5; s2 = 8; s3 = 8;
            Console.Write("_______________________________________________________________________\n");
            Console.Write("                           Region 2a, 2b, 2c                           \n");
            Console.Write("_______________________________________________________________________\n");
            Console.Write("              Table 24             |               Table 29            \n");
            Console.Write("___________________________________|___________________________________\n");
            Console.Write(" p/MPa  h/(kJ/kg)        T/K       | p/MPa s/(kJ/kg-K)       T/K       \n");
            Console.Write("___________________________________|___________________________________\n");
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p1, h1, Api.T_phmass(p1, h1),
                                                                                  p1s, s1, Api.T_psmass(p1s, s1)));
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p2, h2, Api.T_phmass(p2, h2),
                                                                                    p2s, s2, Api.T_psmass(p2s, s2)));
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p3, h3, Api.T_phmass(p3, h3),
                                                                                    p3s, s3, Api.T_psmass(p3s, s3)));
            Console.Write("___________________________________|___________________________________\n");
            p1 = 5; p2 = 5; p3 = 25; h1 = 3500; h2 = 4000; h3 = 3500;
            p1s = 8; p2s = 8; p3s = 90; s1 = 6; s2 = 7.5; s3 = 6;
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p1, h1, Api.T_phmass(p1, h1),
                                                                                    p1s, s1, Api.T_psmass(p1s, s1)));
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p2, h2, Api.T_phmass(p2, h2),
                                                                                    p2s, s2, Api.T_psmass(p2s, s2)));
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p3, h3, Api.T_phmass(p3, h3),
                                                                                    p3s, s3, Api.T_psmass(p3s, s3)));
            Console.Write("_______________________________________________________________________\n\n");
            p1 = 40; p2 = 60; p3 = 60; h1 = 2700; h2 = 2700; h3 = 3200;
            p1s = 20; p2s = 80; p3s = 80; s1 = 5.75; s2 = 5.25; s3 = 5.75;
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p1, h1, Api.T_phmass(p1, h1),
                                                                                  p1s, s1, Api.T_psmass(p1s, s1)));
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p2, h2, Api.T_phmass(p2, h2),
                                                                                  p2s, s2, Api.T_psmass(p2s, s2)));
            Console.Write(string.Format("{0:F3}  {1:F}       {2:E8} {3:F1} {4:F1}       {5:E8}\n", p3, h3, Api.T_phmass(p3, h3),
                                                                                  p3s, s3, Api.T_psmass(p3s, s3)));
            Console.Write("_______________________________________________________________________\n\n\n");

            Console.Write("_______________________________________________________________________\n");
            Console.Write("                             Region 3a, 3b                             \n");
            Console.Write("_______________________________________________________________________\n");
            Console.Write("              Table 5*             |               Table 12*           \n");
            Console.Write("___________________________________|___________________________________\n");
            Console.Write(" p/MPa  h/(kJ/kg)        T/K       | p/MPa s/(kJ/kg-K)       T/K       \n");
            Console.Write("___________________________________|___________________________________\n");
            p1 = 20; p2 = 50; p3 = 100; h1 = 1700; h2 = 2000; h3 = 2100;
            s1 = 3.8; s2 = 3.6; s3 = 4.0;
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F} {4:F1}       {5:E8}\n", p1, h1, Api.T_phmass(p1, h1),
                                                                                  p1, s1, Api.T_psmass(p1, s1)));
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F} {4:F1}       {5:E8}\n", p2, h2, Api.T_phmass(p2, h2),
                                                                                  p2, s2, Api.T_psmass(p2, s2)));
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F} {4:F1}       {5:E8}\n", p3, h3, Api.T_phmass(p3, h3),
                                                                                  p3, s3, Api.T_psmass(p3, s3)));
            Console.Write("_______________________________________________________________________\n\n");
            h1 = 2500; h2 = 2400; h3 = 2700;
            s1 = 5.0; s2 = 4.5; s3 = 5.0;
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F} {4:F1}       {5:E8}\n", p1, h1, Api.T_phmass(p1, h1),
                                                                                  p1, s1, Api.T_psmass(p1, s1)));
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F} {4:F1}       {5:E8}\n", p2, h2, Api.T_phmass(p2, h2),
                                                                                  p2, s2, Api.T_psmass(p2, s2)));
            Console.Write(string.Format("{0:F}  {1:F}        {2:E8} {3:F} {4:F1}       {5:E8}\n", p3, h3, Api.T_phmass(p3, h3),
                                                                                  p3, s3, Api.T_psmass(p3, s3)));
            Console.Write("_______________________________________________________________________\n");
            Console.Write("* The Region 3a, 3b formulation comes from the 2014 \"Revised \n");
            Console.Write("  Supplementary Release on Backward Equations for the Functions\n");
            Console.Write("  T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 [IAPWS SR3-03(2014)]\" \n\n\n");

            Console.Write("**************** Reverse Functions p(h,s) & T(h,s) ********************\n");
            Console.Write(" Verification tables taken from the \"Revised Supplementary Release on \n");
            Console.Write(" Backward Equations for Pressure as a Function of Enthalpy and Entropy \n");
            Console.Write(" p(h,s) for Regions 1 and 2 of the IAPWS IF97\" [IAPWS SR2-01(2014)]   \n");
            Console.Write("***********************************************************************\n\n");

            h1 = 0.001; h2 = 90; h3 = 1500;
            s1 = 0; s2 = 0; s3 = 3.4;
            Console.Write("_________________________________________________\n");
            Console.Write("                     Region 1                    \n");
            Console.Write("_________________________________________________\n");
            Console.Write("                      Table 3                    \n");
            Console.Write("_________________________________________________\n");
            Console.Write("     h/(kJ/kg)     s/(kJ/kg-K)       P/MPa       \n");
            Console.Write("_________________________________________________\n");
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h1, s1, Api.p_hsmass(h1, s1)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h2, s2, Api.p_hsmass(h2, s2)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h3, s3, Api.p_hsmass(h3, s3)));
            Console.Write("_________________________________________________\n\n\n");


            Console.Write("_________________________________________________\n");
            Console.Write("                 Region 2a, 2b, 2c               \n");
            Console.Write("_________________________________________________\n");
            Console.Write("                      Table 9                    \n");
            Console.Write("_________________________________________________\n");
            Console.Write("     h/(kJ/kg)     s/(kJ/kg-K)       P/MPa       \n");
            Console.Write("_________________________________________________\n");
            h1 = 2800; h2 = 2800; h3 = 4100;
            s1 = 6.5; s2 = 9.5; s3 = 9.5;
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h1, s1, Api.p_hsmass(h1, s1)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h2, s2, Api.p_hsmass(h2, s2)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h3, s3, Api.p_hsmass(h3, s3)));
            Console.Write("_________________________________________________\n");
            h1 = 2800; h2 = 3600; h3 = 3600;
            s1 = 6.0; s2 = 6.0; s3 = 7.0;
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h1, s1, Api.p_hsmass(h1, s1)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h2, s2, Api.p_hsmass(h2, s2)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h3, s3, Api.p_hsmass(h3, s3)));
            Console.Write("_________________________________________________\n");
            h1 = 2800; h2 = 2800; h3 = 3400;
            s1 = 5.1; s2 = 5.8; s3 = 5.8;
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h1, s1, Api.p_hsmass(h1, s1)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h2, s2, Api.p_hsmass(h2, s2)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h3, s3, Api.p_hsmass(h3, s3)));
            Console.Write("_________________________________________________\n\n\n");


            Console.Write("_________________________________________________\n");
            Console.Write("                   Region 3a, 3b                 \n");
            Console.Write("_________________________________________________\n");
            Console.Write("                      Table 5*                   \n");
            Console.Write("_________________________________________________\n");
            Console.Write("     h/(kJ/kg)     s/(kJ/kg-K)       P/MPa       \n");
            Console.Write("_________________________________________________\n");
            h1 = 1700; h2 = 2000; h3 = 2100;
            s1 = 3.8; s2 = 4.2; s3 = 4.3;
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h1, s1, Api.p_hsmass(h1, s1)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h2, s2, Api.p_hsmass(h2, s2)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h3, s3, Api.p_hsmass(h3, s3)));
            Console.Write("_________________________________________________\n");
            h1 = 2600; h2 = 2400; h3 = 2700;
            s1 = 5.1; s2 = 4.7; s3 = 5.0;
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h1, s1, Api.p_hsmass(h1, s1)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h2, s2, Api.p_hsmass(h2, s2)));
            Console.Write(string.Format("{0:F3}     {1:F1}       {2:E9}\n", h3, s3, Api.p_hsmass(h3, s3)));
            Console.Write("_________________________________________________\n");
            Console.Write("* The Region 3a, 3b formulation comes from the 2014 \"Revised \n");
            Console.Write("  Supplementary Release on Backward Equations p(h,s) for Region 3,\n");
            Console.Write("  Equations as a Function of h and s for the Region Boundaries,\n");
            Console.Write("  and and Equation Tsat(h,s) for Region 4 of the IAPWS IF97...\" \n");
            Console.Write("  [IAPWS SR4-04(2014)]\" \n\n\n");

            Console.Write("_________________________________________________\n");
            Console.Write("    Region 4 (2-phase Saturation Temperature)    \n");
            Console.Write("_________________________________________________\n");
            Console.Write("                      Table 29*                  \n");
            Console.Write("_________________________________________________\n");
            Console.Write("     h/(kJ/kg)     s/(kJ/kg-K)        T/K        \n");
            Console.Write("_________________________________________________\n");
            h1 = 1800; h2 = 2400; h3 = 2500;
            s1 = 5.3; s2 = 6.0; s3 = 5.5;
            Console.Write(string.Format("{0:F3}     {1:F1}        {2:E9}\n", h1, s1, Api.T_hsmass(h1, s1)));
            Console.Write(string.Format("{0:F3}     {1:F1}        {2:E9}\n", h2, s2, Api.T_hsmass(h2, s2)));
            Console.Write(string.Format("{0:F3}     {1:F1}        {2:E9}\n", h3, s3, Api.T_hsmass(h3, s3)));
            Console.Write("_________________________________________________\n");
            Console.Write("* The Region 3a, 3b formulation comes from the 2014 \"Revised \n");
            Console.Write("  Supplementary Release on Backward Equations p(h,s) for Region 3,\n");
            Console.Write("  Equations as a Function of h and s for the Region Boundaries,\n");
            Console.Write("  and and Equation Tsat(h,s) for Region 4 of the IAPWS IF97...\" \n");
            Console.Write("  [IAPWS SR4-04(2014)]\" \n\n\n");
        }
    }
}
