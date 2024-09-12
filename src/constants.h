#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1

#include <cmath>

const bool Verb = false; //turns on/off verbose output

const double m0 = 931.478;
const double c = 30.;

const double vfact = c/sqrt(m0);  //velocity (cm/ns) = vfact *sqrt(2.*E(MeV)/A(amu))
const double pi = acos(-1.);
const double rad_to_deg = acos(-1.)/180.;
//masses excesses from AME2016 compilation
const double mass_n = 8.07132;
const double mass_p = 7.28897;
const double mass_d = 13.13572;
const double mass_t = 14.9498;
const double mass_3He = 14.93121;
const double mass_alpha = 2.42491;
const double mass_6He = 17.5921;
const double mass_8He = 31.6096;
const double mass_5Li = 11.678886;
const double mass_6Li = 14.0868;
const double mass_7Li = 14.9071;
const double mass_8Li = 20.9458;
const double mass_9Li = 24.9549;
const double mass_6Be = 18.375033;
const double mass_7Be = 15.768999;
const double mass_8Be = 4.9416;
const double mass_9Be = 11.3484;
const double mass_10Be = 12.6074;
const double mass_11Be = 20.1771;
const double mass_8B = 22.9215;
const double mass_9B = 12.416488;
const double mass_10B = 12.0506;
const double mass_11B = 8.6677;
const double mass_9C = 28.910972;
const double mass_10C = 15.698672;
const double mass_11C = 10.649396;
const double mass_12C = 0.;
const double mass_13C = 3.12500888;
const double mass_14C = 3.019892;
const double mass_15C = 9.8731;
const double mass_11N = 24.303559;
const double mass_12N = 17.338068;
const double mass_13N = 5.345481;
const double mass_14N = 2.863416;
const double mass_15N = 0.101438;
const double mass_16N = 5.6839;
const double mass_17N = 7.8700;
const double mass_13O = 23.115432;
const double mass_14O = 8.007781;
const double mass_15O = 2.855605;
const double mass_16O = -4.737001;
const double mass_17O = -0.808763;
const double mass_18O = -0.7828;
const double mass_19O = 3.333;
const double mass_14F = 31.964402;
const double mass_15F = 16.566751;
const double mass_16F = 10.675;
const double mass_17F = 1.951702;
const double mass_18F = .873113;
const double mass_19F = -1.4874;
const double mass_20F = -0.01746;
const double mass_21F = -0.0476;
const double mass_15Ne = 40.215;
const double mass_16Ne = 23.987;
const double mass_17Ne = 16.500447;
const double mass_18Ne = 5.317614;
const double mass_19Ne = 1.7520;
const double mass_20Ne = -7.0419;
const double mass_21Ne = -5.7317;
const double mass_22Ne = -8.0247;
const double mass_19Na = 12.9293;

//total masses
const double Mass_n = m0+mass_n;
const double Mass_p = m0+mass_p;
const double Mass_d = 2.*m0+mass_d;
const double Mass_t = 3.*m0+mass_t;
const double Mass_3He = 3.*m0+mass_3He;
const double Mass_alpha = 4.*m0+mass_alpha;
const double Mass_6He = 6.*m0+mass_6He;
const double Mass_8He = 8.*m0+mass_8He;
const double Mass_5Li = 5.*m0+mass_5Li;
const double Mass_6Li = 6.*m0+mass_6Li;
const double Mass_7Li = 7.*m0+mass_7Li;
const double Mass_8Li = 8.*m0+mass_8Li;
const double Mass_9Li = 9.*m0+mass_6Li;
const double Mass_6Be = 6.*m0+mass_6Be;
const double Mass_7Be = 7.*m0+mass_7Be;
const double Mass_8Be = 8.*m0+mass_8Be;
const double Mass_9Be = 9.*m0+mass_9Be;
const double Mass_10Be = 10.*m0+mass_10Be;
const double Mass_11Be = 11.*m0+mass_11Be;
const double Mass_8B = 8.*m0+mass_8B;
const double Mass_9B = 9.*m0+mass_9B;
const double Mass_10B = 10.*m0+mass_10B;
const double Mass_11B = 11.*m0+mass_11B;
const double Mass_9C = 9.*m0+mass_9C;
const double Mass_10C = 10.*m0+mass_10C;
const double Mass_11C = 11.*m0+mass_11C;
const double Mass_12C = 12.*m0+mass_12C;
const double Mass_13C = 13.*m0+mass_13C;
const double Mass_14C = 14.*m0+mass_14C;
const double Mass_15C = 15.*m0+mass_15C;
const double Mass_11N = 11.*m0+mass_11N;
const double Mass_12N = 12.*m0+mass_12N;
const double Mass_13N = 13.*m0+mass_13N;
const double Mass_14N = 14.*m0+mass_14N;
const double Mass_15N = 15.*m0+mass_15N;
const double Mass_16N = 16.*m0+mass_16N;
const double Mass_17N = 17.*m0+mass_17N;
const double Mass_13O = 13.*m0+mass_13O;
const double Mass_14O = 14.*m0+mass_14O;
const double Mass_15O = 15.*m0+mass_15O;
const double Mass_16O = 16.*m0+mass_16O;
const double Mass_17O = 17.*m0+mass_17O;
const double Mass_18O = 18.*m0+mass_18O;
const double Mass_19O = 19.*m0+mass_19O;
const double Mass_14F = 14.*m0+mass_14F;
const double Mass_15F = 15.*m0+mass_15F;
const double Mass_16F = 16.*m0+mass_16F;
const double Mass_17F = 17.*m0+mass_17F;
const double Mass_18F = 18.*m0+mass_18F;
const double Mass_19F = 19.*m0+mass_19F;
const double Mass_20F = 20.*m0+mass_20F;
const double Mass_21F = 21.*m0+mass_21F;
const double Mass_15Ne = 15.*m0+mass_15Ne;
const double Mass_16Ne = 16.*m0+mass_16Ne;
const double Mass_17Ne = 17.*m0+mass_17Ne;
const double Mass_18Ne = 18.*m0 + mass_18Ne;
const double Mass_19Ne = 19.*m0 + mass_19Ne;
const double Mass_20Ne = 20.*m0 + mass_20Ne;
const double Mass_21Ne = 21.*m0 + mass_21Ne;
const double Mass_22Ne = 22.*m0 + mass_22Ne;
const double Mass_19Na = 19.*m0 + mass_19Na;


#endif
