/* The copyright in this software is being made available under the BSD
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2010-2016, ITU/ISO/IEC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  * Neither the name of the ITU/ISO/IEC nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/** \file     TEncSearch.cpp
 \brief    encoder search class
 */

#include "TLibCommon/CommonDef.h"
#include "TLibCommon/TComRom.h"
#include "TLibCommon/TComMotionInfo.h"
#include "TEncSearch.h"
#include "TLibCommon/TComTU.h"
#include "TLibCommon/Debug.h"
#include <math.h>
#include <limits>
#include <fstream>
using namespace std;
// variables decleration

int counter_i ;
long int array[100000];

extern int counter_ME=0;
extern int counter_FME=0;
//extern int CNT;
int CTUH, CTUW;
int CTUH1, CTUW1;
int CTUH2, CTUW2;
int CTUH3, CTUW3;
int index_ref=0;
double A, B, C2, D, E, F;
double  C, H1, H2, V1, V2, U1, U2, U3, U4;
double A00, A10, A20, AN10, AN20, A01, A02,A0N1,A0N2, A1N1, A2N1, AN1N1, AN2N1, A1N2, A2N2, AN1N2, AN2N2, A11, A21, AN11, AN21, A12, A22, AN12, AN22;
double  R31, R32, R33, R34, R35, R36, R41, R42, R43, R44, R45, R46, R61, R62, R63, R64, R65, R66, R71, R72, R73, R74, R75, R76;
//double C1, C2, C3, C4, C5, C6, C7, C8, C9;
double F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11, F12, F13, F14, F15, F16, F17, F18, F19, F20, F21, F22, F23, F24, F25;
double  F34YLD, F34YLU, F34RYD, F34RYU, F34YN, F34YP, F34XN, F34XP;
double SADD1U, SADD1D, SADD2U, SADD2D, SADD3U, SADD3D, SADD4U, SADD4D, SADD5U, SADD5D;
double FN2, FN1, F0;
long int SAD_Best;
signed short MVX_Best;
signed short MVY_Best;
double CH, CW;
double XMIN;
double YMIN;
double MT_VCX=0;
double MT_VCY=0;
signed short MVX_HALF=0;
signed short MVX_QRTER=0;
signed short MVY_HALF=0;
signed short MVY_QRTER=0;
double CNYL1=0;
double  CNYL2=0;

double FQRTERYL=0;
double FHALFYL=0;
double CNYR1=0;
double   CNYR2=0;
double FQRTERYR=0;
double FHALFYR=0;
double CNROW11D=0;
double  CNROW12D=0;
double FQRTERD1=0;
double FHALFYD1=0;
double CNROW21D=0;
double  CNROW22D=0;
double FQRTERD2=0;
double FHALFYD2=0;
double	CNX1=0;
double   CNX2=0;
double FQRTERX=0;
double FHALFX=0;
double FQRTERXN=0;
double FHALFXN=0;
double FQRTERY=0;
double FHALFY=0;
double FQRTERYN=0;
double FHALFYN=0;
double	CNY1=0;
double   CNY2=0;
int flag_start;
int flag_2point;
int flag_star;

// NEW


double FQRTERYLD;
double FHALFYLD;
double FQRTERYLU;
double FHALFYLU;

double FQRTERYRD;
double FHALFYRD;
double FQRTERYRU;
double FHALFYRU;
double CNROW11;
double CNROW12;

double CNROW21;
double CNROW22;
double CNROW41;
double CNROW42;
double CNROW51;
double CNROW52;




double SADD1, SADD2, SADD3, SADD4, SADD5, SADD6, SADD7, SADD8, SADD9, SADD10, SADD11, SADD12, SADD13, SADD14, SADD15, SADD16, SADD17, SADD18, SADD19, SADD20, SADD21, SADD22, SADD23, SADD24, SADD25;
double SAD1, SAD2,SAD3, SAD4, SAD5, SAD6, SAD7, SAD8,SAD9, SAD10,SAD11,SAD12,SAD13,SAD14,SAD15, SAD16, SAD17,SAD18, SAD19, SAD20, SAD21, SAD22,SAD23, SAD24, SAD25;
double SAD26, SAD27, SAD28, SAD29, SAD30, SAD31, SAD32, SAD33, SAD34, SAD35, SAD36, SAD37, SAD38, SAD39, SAD40, SAD41, SAD42, SAD43, SAD44, SAD45, SAD46, SAD47, SAD48, SAD49;
//END


// ehab modification
double IN1, IN2, IN3, IN4, IN5, IN6, IN7, IN8, IN9;



double input1 = (IN1 - 20243.63083) / 140021.1865;
double input2 = (IN2 - 15269.39007) / 115404.6936;
double input3 = (IN3 - 21755.52693) / 150640.1796;
double input4 = (IN4 - 16539.1768) / 128876.155;
double input5 = (IN5 - 8755.97938) / 98305.47864;
double input6 = (IN6 - 16416.20724) / 128638.1432;
double input7 = (IN7 - 22069.84058) / 151819.7337;
double input8 = (IN8 - 15592.27331) / 116696.186;
double input9 = (IN9 - 20410.4815) / 141399.6972;


double IN[9] = { input1, input2, input3, input4, input5, input6, input7, input8, input9 };
double X1[22];
double X2[20];
double OUT[49];

double in_h1[22][9] = {
	{ 0.4188, -0.8390, 0.3634, 1.3826, -0.7641, -1.4384, -0.4837, 1.3704, 0.2488 },
	{ -3.1622, -2.7183, 2.5697, 0.1008, -1.3657, 0.9375, 1.1973, -0.9195, -2.4992 },
	{ -0.8567, -0.5960, -0.1544, -0.5812, -0.9051, 0.0221, -1.1282, -0.3519, -1.5814 },
	{ -0.0985, 2.0418, -0.9661, -0.1207, 2.2722, 0.2192, -0.0041, -0.4663, -0.2389 },
	{ 0.1930, -0.5540, -0.1735, -0.1261, -0.9707, 0.5252, 0.1230, 1.0118, -0.6391 },
	{ -0.2402, 0.2126, 1.3759, 0.5494, -3.7594, -0.1021, 0.6856, -0.4631, 0.4901 },
	{ 2.1232, 1.7249, 0.0937, -0.9380, -6.6722, 0.1826, 0.7575, -2.5180, 0.4793 },
	{ 0.6028, -1.5483, 0.5561, -0.0741, -0.3986, 0.2925, -1.3962, 0.7853, 0.3631 },
	{ 1.0222, 0.7314, 1.1179, 0.8878, 1.5671, 0.0654, 0.9228, 0.2818, 0.6221 },
	{ -1.4192, 1.0719, 0.0106, -2.4433, 0.3615, 0.7044, 0.1426, 0.9880, -0.4487 },
	{ 0.2109, 0.6300, 1.6213, 0.0885, -2.2724, -0.5426, -0.0174, 0.7288, 0.2117 },
	{ -0.4780, 0.4952, -0.7796, 0.3969, -0.3429, -0.5198, 0.0804, -1.0327, 0.1319 },
	{ -1.6205, -1.6765, -0.5206, -0.9994, -0.2993, -0.5279, -0.6249, -1.0725, -2.0734 },
	{ 1.3951, 0.4778, -0.5137, 0.7893, 0.1896, 1.4127, 0.0945, -1.1629, 0.9523 },
	{ -0.8336, 0.2422, 1.1381, -1.1056, -0.3907, 0.4724, 1.4863, 1.0001, -0.1881 },
	{ -0.5792, -0.2834, 0.3922, 0.7466, 0.5168, 0.5310, 1.1926, -0.8168, -0.3209 },
	{ 1.1406, 1.5332, -2.6224, -0.7757, 0.9054, -1.3309, -1.6413, -0.6074, 0.9761 },
	{ 0.7994, 0.0626, 0.4927, -0.1799, -8.6309, -0.1736, 0.2583, 0.0840, 0.1358 },
	{ 0.6231, -0.3739, -0.6663, -2.0499, -2.2619, 1.4795, 0.5982, -0.3513, 0.5603 },
	{ -1.0471, -1.0738, -0.8969, -1.1787, -1.6139, -1.0951, -1.3595, -1.4937, -1.4949 },
	{ -0.3897, 1.5351, 0.7348, 0.6107, 2.1876, -0.3406, -0.4328, -1.1603, -0.1508 },
	{ -0.3646, 0.5300, 1.2964, 0.4213, -1.1555, -0.4844, -0.0847, 0.2718, -0.8641 }
};

double b1[22] = { 0.0968, -0.0250, -0.2454, 0.1010, 0.1318, 0.0542, 0.1058, 0.2305, -0.1033
, 0.0374, 0.0646, -0.4426, -0.4149, 0.2440, 0.0675, 0.1511, -0.0579, -0.0196, 0.0649, -0.4930
, 0.0456, 0.0823 };

double h1_h2[20][22] = {
	{ -7.7249e-01, 2.6431e-01, -9.6478e+00, 5.7659e-01, -6.8524e-01, -8.5642e-01
	, 5.7903e-01, -2.4538e-01, -5.5378e-01, -4.0176e-01, -9.7965e-01, -4.1804e-01
	, -1.6703e+00, 5.5812e-01, -3.5718e-01, 6.3046e-01, -6.5166e-02, -2.3084e+00
	, 2.0373e-01, -1.4436e+01, 1.0545e+00, -8.0646e-01 },
	{ -9.2299e-01, -1.3398e-01, 1.2166e+01, -8.9669e-01, 3.9931e-01, 1.4148e+00
	, 3.4808e-01, -4.6797e-01, -8.4117e-01, 1.3640e-01, 5.2777e-01, -6.4318e-02
	, 1.9923e+00, -6.9894e-03, -6.9813e-01, -1.2943e-01, -4.8058e-01, 3.6424e+00
	, 3.3312e-01, 2.1540e+01, -6.3291e-01, 6.1915e-01 },
	{ -4.4339e-01, -7.3786e-01, 7.4142e-02, -2.5513e-02, -6.1128e-02, 6.2453e-02
	, -3.7724e-01, 2.3511e-02, 1.6668e-01, 1.0700e-01, -3.0010e-01, 2.3017e-01
	, -4.2833e-01, -6.1653e-01, -3.6665e-01, -7.4501e-01, 6.3364e-01, 1.3370e-01
	, -6.3892e-01, -1.6140e-01, -2.5458e-01, -6.1352e-01 },
	{ -2.5342e-01, -1.1557e+00, -8.8972e+00, 1.1200e+00, -3.8263e-01, -1.2180e+00
	, 3.7304e-01, -1.5087e-01, -1.7179e-01, 1.4643e+00, 6.5448e-01, 1.3033e-01
	, -1.8577e+00, -1.6357e-01, 3.6182e-01, -9.3574e-01, 1.4875e+00, -1.1326e+00
	, -2.2030e-01, -1.3803e+01, 5.7300e-01, 7.9046e-02 },
	{ 1.9505e-01, -9.3900e-02, 8.7619e+00, -2.5106e-01, 5.5512e-01, -1.3546e+00
	, -4.2487e-01, 5.8492e-02, -8.7926e-01, -1.1287e-01, -2.5815e-01, -9.6214e-02
	, 2.3680e+00, -1.1635e+00, 3.3251e-01, -2.1654e-01, 4.4577e-02, 9.1294e-02
	, 6.0819e-01, 1.7054e+01, 9.0263e-02, -3.8808e-01 },
	{ 1.3376e+00, 3.6215e-01, -3.2805e+00, -8.0277e-01, -2.2490e-01, 9.0664e-01
	, -4.5395e-01, -2.8347e-01, -5.1327e-01, -4.8614e-01, 1.1241e+00, 4.0925e-01
	, -2.4976e+00, -7.5775e-02, -8.3905e-01, -7.7445e-03, 1.8012e-01, 2.3571e-01
	, -2.9843e+00, -5.3056e+00, -3.1086e-01, 1.0124e+00 },
	{ 1.0320e+00, -3.8213e-03, 3.9394e+00, -3.6281e-01, 9.9283e-01, -5.9571e-01
	, -2.6135e+00, 9.8830e-01, 2.3075e-02, 1.1617e-01, -7.3448e-02, 2.7877e-01
	, 9.0161e-02, -1.5880e-01, 1.5812e-01, 3.6147e-01, -8.7079e-01, -2.1118e+00
	, 2.1737e-01, 5.7140e+00, 7.2454e-02, 1.2460e-01 },
	{ 1.4191e+00, -1.3795e+00, 1.4226e+01, 3.5731e-01, -4.3278e-01, -1.3036e-01
	, 1.3185e+00, 6.3229e-02, -1.5569e-01, -2.5354e+00, 2.1023e-01, -3.3875e-01
	, 4.5020e+00, 2.8151e-01, -1.1676e+00, 9.1117e-02, 1.6466e+00, 8.3898e-01
	, -7.5561e-01, 2.2389e+01, 2.4598e-01, -2.6032e-01 },
	{ -1.0307e+00, -9.5127e-01, -1.8096e+00, -7.4152e-01, 5.4390e-01, -6.8462e-01
	, 5.9298e-01, 1.6526e-01, -2.5950e-01, 1.2674e+00, -2.8528e-02, 4.4995e-01
	, -9.6574e-01, 1.5723e-02, 6.8176e-01, -8.3339e-01, 6.3247e-01, 2.5262e-01
	, 1.6069e+00, -1.9496e+00, -2.3618e+00, -6.2261e-01 },
	{ -8.5784e-01, 2.5519e+00, 5.8740e+00, -1.4510e+00, 2.3619e-01, 2.1438e+00
	, -6.2898e-01, -3.9396e-01, -5.9588e-01, -1.3421e+00, 8.4696e-01, 1.7535e-01
	, -3.1885e+00, 8.2163e-02, 2.1861e-01, 7.5717e-01, -2.9919e+00, 1.9338e-01
	, -1.5780e-01, 1.3698e+01, 8.2885e-01, 6.8846e-01 },
	{ -2.0094e+00, 1.8651e-01, 1.3487e+01, 1.0498e+00, 4.6990e-01, 6.5322e-01
	, -4.4782e-01, -2.7279e-01, -6.0518e-01, 1.5358e+00, 7.8206e-01, -1.0685e-01
	, -7.9487e-01, -1.9577e-01, 4.2578e-01, 8.5561e-01, -2.0111e+00, 1.3913e+00
	, 5.7941e-01, 2.1091e+01, 1.4987e-01, 2.1623e-01 },
	{ -1.3040e+00, 7.1162e-01, 1.2711e+01, 2.0955e-01, 3.8926e-01, 6.4854e-01
	, 6.3198e-01, -9.3707e-02, -3.3682e-01, 9.3046e-01, 1.1383e+00, -1.2226e-01
	, -6.1993e-02, 1.6249e-02, 5.8634e-01, -6.4534e-01, -9.0463e-01, 2.0703e+00
	, 3.2393e-01, 2.1376e+01, 9.8415e-03, 1.0165e-01 },
	{ -2.6833e+00, 8.3828e-01, 6.7023e+00, 5.3934e-01, 4.6061e-02, -5.3399e-01
	, 7.2944e-02, -7.2018e-01, -2.4472e-01, 1.6687e+00, -9.6605e-01, -7.9306e-02
	, 3.9227e+00, -4.6172e-01, 1.7079e+00, 8.2681e-02, -2.8339e-01, -2.2906e+00
	, 1.1397e+00, 1.0299e+01, 8.7729e-01, -2.6431e-01 },
	{ 1.2610e+00, -3.7321e-01, 1.1007e+01, 1.0322e+00, 4.6022e-01, -1.4049e-01
	, -1.7034e-01, 4.0909e-01, -1.3915e+00, 8.8395e-01, 1.9081e+00, 5.8712e-02
	, 3.8753e+00, -2.3236e+00, -3.2894e-01, -9.1990e-01, 1.0913e+00, 2.0320e+00
	, -1.5536e+00, 1.7559e+01, -5.3088e-01, 9.8249e-01 },
	{ -4.0896e-01, 7.1712e-01, 1.1668e+01, 1.8688e-01, -2.7097e-01, 8.5820e-01
	, 8.7188e-01, -4.7675e-01, -5.9425e-01, -1.2257e+00, 1.1923e+00, 7.9077e-02
	, -1.7815e+00, 7.9465e-02, -1.5918e-01, 7.2800e-01, -1.5230e+00, 2.1438e+00
	, -2.1096e+00, 2.0425e+01, 1.1798e+00, 9.6891e-01 },
	{ 5.9876e-01, 2.0323e-01, 1.3961e+01, 4.4481e-01, -4.7372e-02, 1.0660e-01
	, -1.0577e+00, 4.4640e-01, -1.4008e-01, -1.1253e+00, 7.7748e-01, -2.0840e-01
	, 5.7582e-01, 1.7009e-01, -5.8938e-01, -6.3158e-02, 2.9879e-01, 2.0206e+00
	, -5.3090e-01, 2.2437e+01, -1.1929e+00, 1.3144e-01 },
	{ 4.0286e-01, 4.9295e-01, 1.2400e+01, -2.3767e+00, 2.8774e-01, 9.6609e-01
	, -1.3715e-01, 6.9674e-01, -5.9024e-01, -1.4304e+00, 6.0725e-01, -1.8678e-04
	, 3.2892e+00, 1.5649e-01, -4.0359e-01, -1.7953e-01, -7.6833e-01, 3.3000e+00
	, -6.6463e-01, 1.9821e+01, -2.7221e+00, 1.5754e-01 },
	{ -1.0893e+00, 1.4672e+00, 1.2701e+01, -3.3829e-01, 8.1630e-01, -8.5677e-01
	, -9.5359e-02, 2.9946e-01, -1.0282e-01, 2.3705e-01, -6.6309e-01, -2.1834e-01
	, 1.8922e+00, 4.5198e-01, -2.2855e-01, 4.1517e-01, -3.3573e-01, 1.2768e+00
	, 3.7947e-01, 2.1938e+01, -4.4816e-01, -5.3077e-01 },
	{ -1.7448e-01, -1.3368e+00, 1.2834e+01, -1.0565e+00, 2.0146e-01, 9.7325e-01
	, 1.6952e+00, 4.3453e-01, -5.4840e-01, -1.3475e+00, 5.4966e-01, -3.9675e-02
	, 1.7955e+00, 1.0098e+00, -5.2014e-01, -2.4874e-01, 9.3802e-01, 2.3873e+00
	, 1.6493e+00, 2.1541e+01, -1.8942e+00, -7.7259e-01 },
	{ 6.1627e-01, -1.2645e+00, 1.2434e+01, 8.5836e-01, 3.1937e-01, 5.5226e-01
	, 5.0396e-01, -4.4247e-01, -5.3835e-01, 1.0235e+00, 1.0570e+00, -3.0606e-01
	, -5.4422e-01, 1.4594e-01, 3.3781e-01, -6.0054e-01, 6.3995e-01, 2.5336e+00
	, -2.8540e-01, 2.1549e+01, -9.7358e-01, 4.3888e-01 }
};

double b2[20] = { 0.2465, 0.1810, 0.1084, 0.2927, -0.0100, 0.2559, 0.2297, 0.4111
, 0.2209, 0.0900, 0.4601, 0.3575, 0.2452, 0.1910, 0.4810, 0.7003, 0.3756, 0.5572
, 0.2701, 0.3773 };

double varin[9] = { 0.1662, 0.3118, 0.1513, 0.5052, 0.0907, 0.6853, 0.2147, 0.4113, 0.4587 };

double var1[22] = { 0.5674, 0.8600, 0.1275, 0.4265, 0.0188, 0.1379, 0.0459, 0.0968
, 0.3672, 0.1667, 0.7585, 0.9622, 0.4536, 0.5257, 0.1639, 0.0152, 0.1569, 0.6782
, 0.6824, 0.4533, 0.4969, 0.4567 };

double var2[20] = { 0.5645, 0.1914, 0.3639, 0.3568, 0.9555, 0.7649, 0.8503, 0.5039
, 0.0313, 0.4370, 0.3266, 0.9795, 0.1352, 0.9763, 0.5507, 0.5641, 0.2268, 0.3005
, 0.4263, 0.4836 };

double h2_out[49][20] = {
	{ -8.8796e-03, -9.3443e-01, -2.2276e-01, -9.4226e-01, -5.1887e+00, 4.2630e-01
	, 6.5824e-01, -2.5433e+00, -1.5217e+00, -7.2961e-02, -8.3379e-02, -3.6639e-01
	, 1.8563e-01, -7.7006e-01, -4.3896e-01, -8.2330e-01, -1.7236e+00, -2.1389e-01
	- 2.0685e+00, -5.4638e-01
	},
	{ 4.7074e-01, -1.1674e+00, 2.1995e-01, -2.2405e+00, -3.6187e+00, -2.8966e-01
	, 1.1311e+00, -2.1008e+00, -2.9138e-01, 2.2138e-01, -2.6763e-01, -3.8986e-01
	, 4.3320e-01, -1.3554e+00, -5.5308e-01, -1.0014e+00, -6.6871e-01, -2.9217e-01
	- 7.1842e-01, -1.5300e+00
	},
	{ 6.3480e-01, -2.0556e+00, 5.6842e-01, -1.1938e+00, -4.2157e+00, -1.8049e-01
	, 1.0510e+00, -1.3859e+00, -3.4426e-01, 4.2034e-01, -4.2664e-01, -7.1688e-01
	, -8.3181e-01, -1.1719e+00, -9.2771e-01, -9.5589e-01, 5.8814e-01, 7.7269e-02
	- 3.2058e-01, -1.6074e+00
	},
	{ 1.0519e+00, -1.1068e+00, -7.0303e-02, -9.5274e-01, -8.3602e+00, -5.4129e-01
	, -5.1777e-01, -6.6817e-01, -1.1922e+00, -3.4785e-01, -8.8002e-01, -7.9787e-01
	, -4.4483e-01, -2.2763e-01, -8.5153e-01, -4.5265e-01, 5.1680e-01, 5.1663e-02
	, -1.6592e-01, -1.3746e+00
	},
	{ 8.3152e-01, -2.2336e+00, -1.3866e-02, -1.1343e+00, -3.3852e+00, -4.6231e-01
	, 1.9046e-02, 1.7627e-02, 2.0916e-01, -1.1375e+00, -1.1419e+00, -1.4611e+00
	, -1.3891e+00, -1.6739e+00, -1.0819e+00, -6.7560e-01, 6.9206e-01, -3.9025e-02
	, 5.3254e-02, -1.3794e+00
	},
	{ 7.1629e-01, -2.4959e+00, -1.9558e-01, -1.2967e+00, -2.6469e+00, 7.5783e-02
	, -1.6102e-01, 1.8379e-01, -2.4085e-01, -2.3470e+00, -1.2259e+00, -1.4597e+00
	, -1.0507e+00, -1.4831e+00, -1.1225e+00, -6.1025e-01, 1.0113e-04, -3.3680e-01
	, 1.4678e-01, -1.2923e+00
	},
	{ 3.5613e-01, -2.2676e+00, -9.7936e-02, -7.7185e-02, -5.4554e+00, -4.6930e-01
	, -1.0870e-01, 9.0707e-02, 2.7009e-01, -3.1112e+00, -8.9237e-01, -1.0544e+00
	, -1.5368e+00, -9.8739e-01, -1.1964e+00, -5.7208e-01, -1.0080e+00, -1.6637e-01
	, -4.3348e-01, -1.0231e+00
	},
	{ -1.2834e+00, -6.8776e-01, 3.7490e-01, -4.1912e-01, -2.5641e+00, 4.6594e-01
	, 6.7419e-01, -2.2200e+00, -1.2352e-01, -6.5136e-01, 2.7605e-01, -2.4444e-01
	, 3.5218e-01, -4.1196e-01, -5.7811e-01, -9.5387e-01, -2.1090e+00, -7.4350e-01
	, -1.8140e+00, -4.5688e-01
	},
	{ -8.2038e-01, -6.6839e-01, -6.4986e-02, -2.1472e+00, -1.3892e+00, 2.1598e-01
	, 1.1966e+00, -2.1047e+00, -7.9914e-02, 5.3176e-03, -2.4808e-03, -3.0155e-02
	, 6.3824e-01, 1.3777e-01, -1.7073e-01, -1.0345e+00, -1.6662e+00, -5.5943e-01
	, -9.6971e-01, -7.0671e-01
	},
	{ -7.3914e-01, -2.4241e+00, 3.9511e-01, -1.9270e+00, -3.8525e+00, 5.7844e-01
	, 1.4038e+00, -1.6548e+00, 8.3917e-02, 1.3292e-01, -2.2919e-01, -4.2001e-01
	, -5.0191e-01, -2.5788e-01, -9.5524e-01, -1.1329e+00, -2.1894e-01, -3.9584e-01
	, -2.4279e-01, -1.1363e+00
	},
	{ -1.6371e-02, -5.0738e-01, 8.1411e-01, -1.0688e+00, -2.9087e+00, -9.7867e-02
	, 1.0418e+00, -4.9857e-01, 1.6778e-01, 3.8505e-02, -8.6645e-01, -1.0691e+00
	, -4.7543e-01, -3.6935e-01, -6.6374e-01, -6.6293e-01, 1.5834e+00, 2.1931e-01
	, -1.3344e-01, -1.6509e+00
	},
	{ -1.5738e-01, -2.4066e+00, -3.5091e-01, -8.8250e-01, -3.7278e+00, 4.3627e-01
	, 6.2754e-01, -1.2864e-01, 2.5365e-01, -1.4086e+00, -1.1656e+00, -1.3009e+00
	, -7.3654e-01, -1.1468e+00, -1.4129e+00, -5.4417e-01, 4.9812e-01, -4.3431e-01
	, 1.0757e-01, -1.1817e+00
	},
	{ -2.0766e-01, -2.6832e+00, -2.1598e-01, -1.2289e+00, -8.3120e-01, 4.4710e-01
	, -7.3113e-02, 3.8168e-01, -8.4527e-03, -2.9213e+00, -1.0906e+00, -1.0889e+00
	, -3.3709e-01, -8.6031e-01, -1.3418e+00, -4.3204e-01, -2.1498e-01, -4.3794e-01
	, 6.6690e-02, -8.7411e-01
	},
	{ -8.0709e-01, -2.6568e+00, -3.0165e-02, 3.5911e-01, -2.9098e+00, 2.4462e-01
	, -2.3430e-01, 1.1504e-01, 6.1064e-01, -2.9499e+00, -1.1606e+00, -1.3722e+00
	, -1.7896e+00, -7.4530e-01, -1.5791e+00, -3.9661e-01, -1.1269e+00, -5.1628e-01
	, -3.1986e-01, -6.9310e-01
	},
	{ -3.5598e-01, -1.0476e+00, 1.9286e-01, 5.3692e-01, -2.7188e+00, -4.1726e-01
	, 2.3446e-01, -1.7329e+00, 5.0654e-01, -1.8978e+00, 8.5674e-02, -5.7106e-01
	, 2.9948e-01, -3.5226e-01, -7.7723e-01, -9.2975e-01, -2.3965e+00, -6.4689e-01
	, -1.3375e+00, -3.8681e-01
	},
	{ -6.5232e-01, -1.3299e+00, -5.1636e-02, -5.5785e-01, -2.0279e+00, -1.4866e+00
	, 2.0022e-01, -1.8307e+00, 4.5635e-01, -6.6030e-01, 2.2247e-02, -3.7371e-01
	, 8.7689e-01, -7.1572e-01, -8.1632e-01, -1.0429e+00, -2.0299e+00, -5.2471e-01
	, -4.8315e-01, -6.4764e-01
	},
	{ -2.8892e-01, -1.2936e-01, 7.4722e-02, -1.4027e+00, -4.3107e+00, -2.1870e-01
	, 5.3461e-01, -2.2106e+00, 4.5546e-01, 2.7846e-01, 2.6829e-01, 2.6110e-02
	, -4.5874e-01, -4.8876e-01, -6.9529e-01, -7.3001e-01, -6.9270e-01, -5.1580e-01
	, -1.3372e-01, -2.7748e-01
	},
	{ -8.4779e-01, -4.1042e-01, 4.8052e-01, -1.6485e+00, -6.7977e+00, -2.7869e-01
	, 3.0455e-01, -1.2556e+00, 5.1304e-01, 1.8535e-01, -3.0734e-01, -2.1726e-01
	, -2.5600e+00, -7.6607e-01, -5.1537e-01, -2.4669e-01, 2.2424e+00, -9.0553e-02
	, -1.2450e-01, -6.1147e-01
	},
	{ -8.2141e-01, -3.6031e+00, -4.5201e-01, -6.8642e-01, -5.5348e+00, 5.0287e-01
	, -4.7636e-02, -2.7401e-02, 3.3190e-01, -2.7540e-01, -1.3013e+00, -9.2780e-01
	, -1.8868e+00, -3.8816e-01, -7.0779e-01, -3.0638e-01, 1.2783e+00, -7.3617e-01
	, 1.0844e-01, -2.5776e-01
	},
	{ -5.8284e-01, -2.9290e+00, -5.8127e-03, -7.1343e-01, -2.5672e+00, 1.4451e+00
	, -5.2150e-01, 3.2066e-01, -8.2624e-01, -2.3322e+00, -1.0415e+00, -1.2489e+00
	, -9.3233e-01, -3.6429e-01, -9.7437e-01, -3.9675e-01, -4.7318e-01, -8.4142e-01
	, -2.0060e-01, -4.8331e-01
	},
	{ -1.0055e+00, -3.2412e+00, 2.8195e-01, 1.0445e+00, -2.5354e+00, 1.4572e+00
	, -6.9751e-01, 6.6011e-02, -2.6899e-01, -3.2484e+00, -1.1111e+00, -1.2310e+00
	, -2.1882e+00, 1.8537e-01, -1.2614e+00, -3.3854e-01, -1.4302e+00, -6.4285e-01
	, -5.2961e-01, -2.9892e-01
	},
	{ -3.7267e-01, -1.8196e-01, 8.0399e-02, 1.2162e+00, -7.2333e+00, -1.1294e+00
	, -4.6020e-01, -1.1028e+00, 1.0718e-01, -2.7367e+00, -5.5830e-01, -1.4790e-01
	, -2.4331e-01, 7.0998e-01, -7.2557e-01, -4.9420e-01, -1.1283e+00, -4.8874e-01
	, -1.0482e+00, -2.2012e-01
	},
	{ -2.3679e-01, 8.2772e-02, -8.0701e-02, 8.5082e-01, -1.9936e+00, -3.5844e+00
	, -3.1974e-01, -1.2630e+00, 6.7901e-01, -2.3278e+00, -2.9890e-01, -1.8305e-01
	, 3.0530e-01, 4.9596e-01, -7.6488e-01, -5.6843e-01, -1.3721e+00, -3.4460e-01
	, -1.2868e-01, -1.2279e-01
	},
	{ -1.0138e+00, 7.1795e-01, 4.0708e-02, 4.8599e-01, -6.8950e+00, -4.0967e+00
	, -1.0551e+00, -1.1680e+00, 5.9334e-01, -1.6236e-02, 1.8102e-01, -1.9450e-01
	, -5.4009e-01, -3.3661e-02, -8.4119e-01, -5.0327e-01, -1.3604e+00, -4.9975e-01
	, 1.2562e-01, 1.3608e-01
	},
	{ -1.7417e+00, 2.7904e+00, -3.4670e-02, -8.5813e-01, 1.2576e+01, -1.6387e+00
	, -1.3490e+00, -1.7936e-01, -9.8649e-01, 6.7953e-02, -8.5740e-02, 9.2728e-02
	, -5.4293e-01, 3.5053e-01, -1.2419e-01, 2.7367e-01, 4.5119e-01, 5.7138e-02
	, -4.4356e-02, 2.2363e-01
	},
	{ -2.9222e+00, -1.2837e+00, 9.1539e-02, 5.2343e-02, -6.2186e+00, 1.0593e+00
	, -1.1400e+00, -2.8466e-01, -2.1597e+00, -6.4877e-02, -7.4170e-01, -3.8495e-01
	, -2.4046e+00, 7.1718e-01, 2.2751e-01, -2.7010e-01, -9.2537e-02, -7.8861e-01
	, 1.1558e-02, 1.6387e-01
	},
	{ -1.8723e+00, -1.1885e+00, 2.3537e-01, 6.6637e-01, -2.8713e+00, 1.4326e+00
	, -3.7156e-01, 2.1280e-01, -1.8041e+00, -1.4957e+00, -9.7072e-01, -5.3555e-01
	, -3.4676e-01, 1.1226e+00, -1.6467e-01, -3.8561e-01, -6.3428e-01, -7.4461e-01
	, -9.8832e-01, -2.1176e-03
	},
	{ -8.3268e-01, -9.1609e-01, 1.1393e-01, 1.0801e+00, -7.0975e+00, 6.3484e-01
	, -3.9767e-01, -6.5457e-01, -8.3609e-01, -2.9016e+00, -8.8460e-01, -2.4927e-01
	, -4.8758e-01, 8.2601e-01, -5.5559e-01, -4.4935e-01, -1.0486e+00, -5.3399e-01
	, -1.1190e+00, -1.7410e-01
	},
	{ -4.8193e-01, -2.2935e+00, 4.5017e-02, 1.2933e+00, -3.8685e+00, -2.3498e+00
	, -9.1214e-01, -5.8650e-01, 1.2468e+00, -3.1203e+00, -6.7203e-01, -7.0250e-01
	, -1.1814e+00, -6.6742e-01, -1.0255e+00, -6.9717e-01, -2.4347e+00, -3.8902e-01
	, -2.4751e-01, -4.7109e-01
	},
	{ 3.0185e-01, -2.1535e+00, -1.5563e-01, 2.2590e-01, -1.4744e+00, -2.2101e+00
	, -1.0128e+00, -6.3995e-01, 5.9887e-01, -3.3811e+00, -3.6816e-01, -9.3643e-02
	, -1.8883e-01, -1.0856e+00, -9.4051e-01, -6.6552e-01, -1.8883e+00, -3.4401e-01
	, 1.0668e-01, -6.9795e-01
	},
	{ 5.1771e-01, -1.4710e+00, -1.3550e-01, 1.2037e-01, -3.4345e+00, -1.7844e+00
	, -2.1069e+00, -1.4560e-01, -4.3460e-01, -8.0962e-01, -8.1033e-02, -8.6725e-02
	, 1.3662e-01, -1.4411e+00, -8.6814e-01, -8.1351e-01, -1.6989e+00, -5.5062e-01
	, 3.4111e-01, -7.8890e-02
	},
	{ 3.3904e-01, 8.2847e-01, -4.3646e-01, -7.8717e-01, -5.7526e+00, -8.0716e-01
	, -2.9486e+00, -7.5468e-02, -1.4636e+00, 5.6076e-02, -3.2006e-01, 6.2459e-02
	, 2.2884e-01, -1.8045e+00, 2.1517e-01, -6.8700e-01, -1.3534e+00, -1.8158e-01
	, 3.7168e-02, -4.6347e-02
	},
	{ -1.6737e+00, -1.1322e+00, 7.7339e-01, -3.4195e-01, -4.7508e+00, 1.0272e+00
	, -1.5029e+00, -2.5778e-01, -2.1104e+00, 2.0086e-01, -3.1178e-01, -2.9652e-01
	, 3.0371e-03, -1.2250e+00, 6.9642e-01, -7.0325e-01, -1.1950e+00, -7.1733e-01
	, -3.3197e-01, -1.6507e-01
	},
	{ -1.6072e+00, -2.0603e+00, 4.3131e-01, -4.4671e-01, -2.7138e+00, 1.5546e+00
	, -2.9501e-01, -4.2162e-01, -2.0792e+00, -3.7016e-01, -6.8997e-01, -4.9821e-01
	, -1.0747e+00, 3.2590e-01, 2.8504e-01, -8.8750e-01, -1.3113e+00, -9.3401e-01
	, -1.7274e+00, -5.3122e-01
	},
	{ -1.0511e+00, -1.9247e+00, -1.4242e-01, 5.5262e-01, -4.3468e+00, 2.0159e+00
	, 1.6395e-01, -1.1034e+00, -1.8419e+00, -1.6754e+00, -3.4201e-01, -6.6965e-01
	, -1.2002e+00, 6.8580e-01, 4.5353e-02, -8.4223e-01, -1.6190e+00, -8.4221e-01
	, -1.7356e+00, -4.6853e-01
	},
	{ -3.9692e-02, -2.4203e+00, 1.5644e-02, 1.0260e+00, -1.9879e+00, -1.9299e+00
	, -1.1852e+00, -3.0705e-02, 2.8723e-01, -3.3955e+00, -8.8404e-01, -8.3973e-01
	, -1.1982e+00, -1.1028e+00, -1.0992e+00, -7.6294e-01, -2.3387e+00, -4.1025e-01
	, -8.2719e-03, -4.3880e-01
	},
	{ 6.2091e-01, -2.1534e+00, 2.8370e-01, -4.6866e-01, -1.6580e-01, -2.2832e+00
	, -1.2162e+00, 4.1357e-02, -1.3659e-01, -3.4686e+00, -6.7962e-01, -4.9695e-01
	, -2.2147e-01, -1.0530e+00, -9.5756e-01, -7.4174e-01, -1.5271e+00, -3.9755e-01
	, 2.6641e-01, -7.7928e-01
	},
	{ 1.2574e+00, -2.3225e+00, -5.7411e-02, -5.3822e-01, -1.7285e+00, -1.5017e+00
	, -1.8813e+00, 7.2878e-02, -6.2272e-01, -1.9075e+00, -4.1389e-01, -1.9220e-01
	, -1.7733e-01, -1.3043e+00, -5.3729e-01, -7.8339e-01, -1.7860e+00, -7.5231e-01
	, 2.6428e-01, -7.4150e-01
	},
	{ 1.6889e+00, -3.1283e-01, 9.7593e-04, -7.2068e-01, -2.6172e+00, -2.4634e+00
	, -1.7930e+00, 1.4391e-01, -1.9280e+00, 1.0069e-01, -5.0860e-01, -5.3526e-01
	, 9.1847e-01, -6.3177e-01, -2.6638e-01, -7.4912e-01, -7.6928e-01, -2.1407e-01
	, -1.2314e-01, -7.4792e-01
	},
	{ 4.1100e-01, -2.0204e+00, -5.6237e-01, -6.3982e-01, -1.5651e+00, -3.0832e-01
	, -1.8866e+00, -3.6622e-01, -1.9920e+00, -1.4897e-01, -2.7793e-01, -5.7470e-01
	, 9.6123e-01, -1.3398e+00, 4.1709e-01, -1.1115e+00, -1.9587e+00, -7.7792e-01
	, -7.5892e-01, -4.9877e-01
	},
	{ -1.6433e+00, -1.1931e+00, 9.2423e-02, -5.4025e-01, -2.3120e+00, 9.6693e-01
	, -1.0290e-01, -8.0921e-01, -2.3514e+00, -2.2863e-02, -2.3164e-02, -2.0438e-01
	, 6.4025e-01, -5.1558e-01, 1.6085e-01, -9.1764e-01, -1.0057e+00, -6.2366e-01
	, -2.2346e+00, -6.1000e-01
	},
	{ -1.6619e+00, -1.1436e+00, -1.9180e-01, 2.3507e-01, -2.9320e+00, 1.6000e+00
	, 1.5664e-01, -1.5466e+00, -2.4510e+00, -5.4767e-01, -3.2955e-02, -3.8746e-01
	, -6.1878e-01, 1.6079e-01, 2.4252e-01, -1.1011e+00, -2.0636e+00, -9.1237e-01
	, -2.1326e+00, -5.3097e-01
	},
	{ 7.9552e-01, -2.1182e+00, 3.9360e-01, 3.5163e-01, -4.5251e+00, -1.4848e+00
	, -5.6838e-01, -1.8277e-01, -8.9536e-02, -3.4784e+00, -8.1300e-01, -9.4205e-01
	, -1.3861e+00, -1.0434e+00, -1.1470e+00, -7.0490e-01, -1.4145e+00, -2.1478e-01
	, -1.4518e-01, -9.8913e-01
	},
	{ 1.4701e+00, -2.0194e+00, 1.7165e-02, -1.3912e+00, -2.7195e+00, -2.0561e+00
	, -1.1073e+00, -1.4142e-01, 4.2658e-01, -2.7501e+00, -8.1945e-01, -8.8771e-01
	, -9.3392e-01, -1.6390e+00, -1.1105e+00, -6.3433e-01, -1.3193e+00, -1.0124e-01
	, 7.3663e-02, -8.5452e-01
	},
	{ 1.5880e+00, -1.6855e+00, 1.0324e-01, -1.1680e+00, -3.1638e+00, -1.1791e+00
	, -1.4786e+00, 1.3161e-01, -9.1269e-02, -1.4099e+00, -7.9533e-01, -9.7736e-01
	, -1.1950e+00, -1.8103e+00, -7.4375e-01, -8.4785e-01, -9.9407e-01, -1.7630e-01
	, 6.3859e-02, -9.9320e-01
	},
	{ 1.3675e+00, -1.0945e+00, 1.8701e-01, -9.8098e-01, -8.5314e+00, -1.3896e+00
	, -1.1424e+00, -6.2893e-01, -1.1345e+00, -3.0533e-01, -8.3957e-01, -6.7144e-01
	, 1.1683e-01, -4.2777e-01, -6.5817e-01, -5.3660e-01, -2.5386e-01, 6.0994e-03
	, -2.1694e-01, -1.2484e+00
	},
	{ 1.5294e+00, -1.4965e+00, -3.1763e-01, -1.1929e+00, -3.7355e+00, -8.0087e-01
	, -1.6059e+00, -8.2569e-01, -1.9619e+00, 3.7679e-01, -4.8385e-01, -7.8596e-01
	, 3.5319e-01, -1.8522e+00, 1.5575e-01, -5.8826e-01, -1.4066e+00, -4.1145e-01
	, -3.9778e-01, -1.2670e+00
	},
	{ 4.9471e-01, -1.4294e+00, 5.3741e-02, -1.7904e+00, -2.6849e+00, 4.5577e-01
	, -3.5644e-01, -1.0962e+00, -2.4250e+00, 3.2019e-01, -4.2371e-01, -7.4663e-01
	, 1.6397e-01, -1.4172e+00, 1.3252e-01, -8.7773e-01, -1.0928e+00, -4.3745e-01
	, -1.8633e+00, -1.1042e+00
	},
	{ -2.2778e-01, -1.5190e+00, 2.2858e-01, -5.8312e-01, -5.0983e+00, 1.0695e+00
	, 2.5056e-01, -1.7277e+00, -2.7653e+00, -1.3574e-01, 2.4409e-02, -4.4029e-01
	, -1.9820e-03, -5.1547e-01, -3.5138e-01, -7.6257e-01, -1.5551e+00, -3.9364e-01
	, -2.4294e+00, -8.8016e-01
	}
};

double bout[49] = { -0.5331, -0.9636, -1.2114, -0.7465, -0.7748, -0.7380, -0.7105
, -1.1860, -1.2326, -1.1127, -1.4849, -1.0036, -0.9175, -0.7242, -0.9370, -0.8835, -1.2616
, -1.1517, -0.8409, -0.8018, -0.6119, -1.0316, -1.4345, -0.7738, -0.1712, -0.5656, -1.3207
, -0.9472, -0.6523, -1.2741, -0.8502, -0.7972, -0.9676, -0.7192, -1.0251, -0.9046, -0.8703
, -1.1919, -1.2313, -0.8956, -1.0853, -1.0595, -0.5924, -0.9088, -0.6481, -0.6158, -0.7292
, -0.5779, -0.5256 };

double relu(double x)
{
	//return std::max(x, 0);
	if (x<0)
	{
		return 0;
	}
	else
	{
		return x;
	}
}

double sigmoid(double x){
	return (1 / (1 + exp(-x)));
}
//end of modification
//! \ingroup TLibEncoder
//! \{

static const TComMv s_acMvRefineH[9] =
{
  TComMv(  0,  0 ), // 0
  TComMv(  0, -1 ), // 1
  TComMv(  0,  1 ), // 2
  TComMv( -1,  0 ), // 3
  TComMv(  1,  0 ), // 4
  TComMv( -1, -1 ), // 5
  TComMv(  1, -1 ), // 6
  TComMv( -1,  1 ), // 7
  TComMv(  1,  1 )  // 8
};

static const TComMv s_acMvRefineQ[9] =
{
  TComMv(  0,  0 ), // 0
  TComMv(  0, -1 ), // 1
  TComMv(  0,  1 ), // 2
  TComMv( -1, -1 ), // 5
  TComMv(  1, -1 ), // 6
  TComMv( -1,  0 ), // 3
  TComMv(  1,  0 ), // 4
  TComMv( -1,  1 ), // 7
  TComMv(  1,  1 )  // 8
};

static Void offsetSubTUCBFs(TComTU &rTu, const ComponentID compID)
{
        TComDataCU *pcCU              = rTu.getCU();
  const UInt        uiTrDepth         = rTu.GetTransformDepthRel();
  const UInt        uiAbsPartIdx      = rTu.GetAbsPartIdxTU(compID);
  const UInt        partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> 1;

  //move the CBFs down a level and set the parent CBF

  UChar subTUCBF[2];
  UChar combinedSubTUCBF = 0;

  for (UInt subTU = 0; subTU < 2; subTU++)
  {
    const UInt subTUAbsPartIdx = uiAbsPartIdx + (subTU * partIdxesPerSubTU);

    subTUCBF[subTU]   = pcCU->getCbf(subTUAbsPartIdx, compID, uiTrDepth);
    combinedSubTUCBF |= subTUCBF[subTU];
  }

  for (UInt subTU = 0; subTU < 2; subTU++)
  {
    const UInt subTUAbsPartIdx = uiAbsPartIdx + (subTU * partIdxesPerSubTU);
    const UChar compositeCBF = (subTUCBF[subTU] << 1) | combinedSubTUCBF;

    pcCU->setCbfPartRange((compositeCBF << uiTrDepth), compID, subTUAbsPartIdx, partIdxesPerSubTU);
  }
}


TEncSearch::TEncSearch()
: m_puhQTTempTrIdx(NULL)
, m_pcQTTempTComYuv(NULL)
, m_pcEncCfg (NULL)
, m_pcTrQuant (NULL)
, m_pcRdCost (NULL)
, m_pcEntropyCoder (NULL)
, m_iSearchRange (0)
, m_bipredSearchRange (0)
, m_motionEstimationSearchMethod (MESEARCH_FULL)
, m_pppcRDSbacCoder (NULL)
, m_pcRDGoOnSbacCoder (NULL)
, m_pTempPel (NULL)
, m_isInitialized (false)
{
  for (UInt ch=0; ch<MAX_NUM_COMPONENT; ch++)
  {
    m_ppcQTTempCoeff[ch]                           = NULL;
#if ADAPTIVE_QP_SELECTION
    m_ppcQTTempArlCoeff[ch]                        = NULL;
#endif
    m_puhQTTempCbf[ch]                             = NULL;
    m_phQTTempCrossComponentPredictionAlpha[ch]    = NULL;
    m_pSharedPredTransformSkip[ch]                 = NULL;
    m_pcQTTempTUCoeff[ch]                          = NULL;
#if ADAPTIVE_QP_SELECTION
    m_ppcQTTempTUArlCoeff[ch]                      = NULL;
#endif
    m_puhQTTempTransformSkipFlag[ch]               = NULL;
  }

  for (Int i=0; i<MAX_NUM_REF_LIST_ADAPT_SR; i++)
  {
    memset (m_aaiAdaptSR[i], 0, MAX_IDX_ADAPT_SR * sizeof (Int));
  }
  for (Int i=0; i<AMVP_MAX_NUM_CANDS+1; i++)
  {
    memset (m_auiMVPIdxCost[i], 0, (AMVP_MAX_NUM_CANDS+1) * sizeof (UInt) );
  }

  setWpScalingDistParam( NULL, -1, REF_PIC_LIST_X );
}


Void TEncSearch::destroy()
{
  assert (m_isInitialized);
  if ( m_pTempPel )
  {
    delete [] m_pTempPel;
    m_pTempPel = NULL;
  }

  if ( m_pcEncCfg )
  {
    const UInt uiNumLayersAllocated = m_pcEncCfg->getQuadtreeTULog2MaxSize()-m_pcEncCfg->getQuadtreeTULog2MinSize()+1;

    for (UInt ch=0; ch<MAX_NUM_COMPONENT; ch++)
    {
      for (UInt layer = 0; layer < uiNumLayersAllocated; layer++)
      {
        delete[] m_ppcQTTempCoeff[ch][layer];
#if ADAPTIVE_QP_SELECTION
        delete[] m_ppcQTTempArlCoeff[ch][layer];
#endif
      }
      delete[] m_ppcQTTempCoeff[ch];
      delete[] m_puhQTTempCbf[ch];
#if ADAPTIVE_QP_SELECTION
      delete[] m_ppcQTTempArlCoeff[ch];
#endif
    }

    for( UInt layer = 0; layer < uiNumLayersAllocated; layer++ )
    {
      m_pcQTTempTComYuv[layer].destroy();
    }
  }

  delete[] m_puhQTTempTrIdx;
  delete[] m_pcQTTempTComYuv;

  for (UInt ch=0; ch<MAX_NUM_COMPONENT; ch++)
  {
    delete[] m_pSharedPredTransformSkip[ch];
    delete[] m_pcQTTempTUCoeff[ch];
#if ADAPTIVE_QP_SELECTION
    delete[] m_ppcQTTempTUArlCoeff[ch];
#endif
    delete[] m_phQTTempCrossComponentPredictionAlpha[ch];
    delete[] m_puhQTTempTransformSkipFlag[ch];
  }
  m_pcQTTempTransformSkipTComYuv.destroy();

  m_tmpYuvPred.destroy();
  m_isInitialized = false;
}

TEncSearch::~TEncSearch()
{
  if (m_isInitialized)
  {
    destroy();
  }
}




Void TEncSearch::init(TEncCfg*       pcEncCfg,
                      TComTrQuant*   pcTrQuant,
                      Int            iSearchRange,
                      Int            bipredSearchRange,
                      MESearchMethod motionEstimationSearchMethod,
                      const UInt     maxCUWidth,
                      const UInt     maxCUHeight,
                      const UInt     maxTotalCUDepth,
                      TEncEntropy*   pcEntropyCoder,
                      TComRdCost*    pcRdCost,
                      TEncSbac***    pppcRDSbacCoder,
                      TEncSbac*      pcRDGoOnSbacCoder
                      )
{
  assert (!m_isInitialized);
  m_pcEncCfg                     = pcEncCfg;
  m_pcTrQuant                    = pcTrQuant;
  m_iSearchRange                 = iSearchRange;
  m_bipredSearchRange            = bipredSearchRange;
  m_motionEstimationSearchMethod = motionEstimationSearchMethod;
  m_pcEntropyCoder               = pcEntropyCoder;
  m_pcRdCost                     = pcRdCost;

  m_pppcRDSbacCoder              = pppcRDSbacCoder;
  m_pcRDGoOnSbacCoder            = pcRDGoOnSbacCoder;
  
  for (UInt iDir = 0; iDir < MAX_NUM_REF_LIST_ADAPT_SR; iDir++)
  {
    for (UInt iRefIdx = 0; iRefIdx < MAX_IDX_ADAPT_SR; iRefIdx++)
    {
      m_aaiAdaptSR[iDir][iRefIdx] = iSearchRange;
    }
  }

  // initialize motion cost
  for( Int iNum = 0; iNum < AMVP_MAX_NUM_CANDS+1; iNum++)
  {
    for( Int iIdx = 0; iIdx < AMVP_MAX_NUM_CANDS; iIdx++)
    {
      if (iIdx < iNum)
      {
        m_auiMVPIdxCost[iIdx][iNum] = xGetMvpIdxBits(iIdx, iNum);
      }
      else
      {
        m_auiMVPIdxCost[iIdx][iNum] = MAX_INT;
      }
    }
  }

  const ChromaFormat cform=pcEncCfg->getChromaFormatIdc();
  initTempBuff(cform);

  m_pTempPel = new Pel[maxCUWidth*maxCUHeight];

  const UInt uiNumLayersToAllocate = pcEncCfg->getQuadtreeTULog2MaxSize()-pcEncCfg->getQuadtreeTULog2MinSize()+1;
  const UInt uiNumPartitions = 1<<(maxTotalCUDepth<<1);
  for (UInt ch=0; ch<MAX_NUM_COMPONENT; ch++)
  {
    const UInt csx=::getComponentScaleX(ComponentID(ch), cform);
    const UInt csy=::getComponentScaleY(ComponentID(ch), cform);
    m_ppcQTTempCoeff[ch] = new TCoeff* [uiNumLayersToAllocate];
#if ADAPTIVE_QP_SELECTION
    m_ppcQTTempArlCoeff[ch]  = new TCoeff*[uiNumLayersToAllocate];
#endif
    m_puhQTTempCbf[ch] = new UChar  [uiNumPartitions];

    for (UInt layer = 0; layer < uiNumLayersToAllocate; layer++)
    {
      m_ppcQTTempCoeff[ch][layer] = new TCoeff[(maxCUWidth*maxCUHeight)>>(csx+csy)];
#if ADAPTIVE_QP_SELECTION
      m_ppcQTTempArlCoeff[ch][layer]  = new TCoeff[(maxCUWidth*maxCUHeight)>>(csx+csy) ];
#endif
    }

    m_phQTTempCrossComponentPredictionAlpha[ch]    = new SChar  [uiNumPartitions];
    m_pSharedPredTransformSkip[ch]                 = new Pel   [MAX_CU_SIZE*MAX_CU_SIZE];
    m_pcQTTempTUCoeff[ch]                          = new TCoeff[MAX_CU_SIZE*MAX_CU_SIZE];
#if ADAPTIVE_QP_SELECTION
    m_ppcQTTempTUArlCoeff[ch]                      = new TCoeff[MAX_CU_SIZE*MAX_CU_SIZE];
#endif
    m_puhQTTempTransformSkipFlag[ch]               = new UChar [uiNumPartitions];
  }
  m_puhQTTempTrIdx   = new UChar  [uiNumPartitions];
  m_pcQTTempTComYuv  = new TComYuv[uiNumLayersToAllocate];
  for( UInt ui = 0; ui < uiNumLayersToAllocate; ++ui )
  {
    m_pcQTTempTComYuv[ui].create( maxCUWidth, maxCUHeight, pcEncCfg->getChromaFormatIdc() );
  }
  m_pcQTTempTransformSkipTComYuv.create( maxCUWidth, maxCUHeight, pcEncCfg->getChromaFormatIdc() );
  m_tmpYuvPred.create(MAX_CU_SIZE, MAX_CU_SIZE, pcEncCfg->getChromaFormatIdc());
  m_isInitialized = true;
}


__inline Void TEncSearch::xTZSearchHelp( const TComPattern* const pcPatternKey, IntTZSearchStruct& rcStruct, const Int iSearchX, const Int iSearchY, const UChar ucPointNr, const UInt uiDistance )
{
  Distortion  uiSad = 0;

  const Pel* const  piRefSrch = rcStruct.piRefY + iSearchY * rcStruct.iYStride + iSearchX;

  //-- jclee for using the SAD function pointer
  m_pcRdCost->setDistParam( pcPatternKey, piRefSrch, rcStruct.iYStride,  m_cDistParam );

  setDistParamComp(COMPONENT_Y);

  // distortion
  m_cDistParam.bitDepth = pcPatternKey->getBitDepthY();
  m_cDistParam.m_maximumDistortionForEarlyExit = rcStruct.uiBestSad;

  if((m_pcEncCfg->getRestrictMESampling() == false) && m_pcEncCfg->getMotionEstimationSearchMethod() == MESEARCH_SELECTIVE)
  {
    Int isubShift = 0;
    // motion cost
    Distortion uiBitCost = m_pcRdCost->getCostOfVectorWithPredictor( iSearchX, iSearchY );

    // Skip search if bit cost is already larger than best SAD
    if (uiBitCost < rcStruct.uiBestSad)
    {
      if ( m_cDistParam.iRows > 32 )
      {
        m_cDistParam.iSubShift = 4;
      }
      else if ( m_cDistParam.iRows > 16 )
      {
        m_cDistParam.iSubShift = 3;
      }
      else if ( m_cDistParam.iRows > 8 )
      {
        m_cDistParam.iSubShift = 2;
      }
      else
      {
        m_cDistParam.iSubShift = 1;
      }

      Distortion uiTempSad = m_cDistParam.DistFunc( &m_cDistParam );
      if((uiTempSad + uiBitCost) < rcStruct.uiBestSad)
      {
        uiSad += uiTempSad >>  m_cDistParam.iSubShift;
        while(m_cDistParam.iSubShift > 0)
        {
          isubShift         = m_cDistParam.iSubShift -1;
          m_cDistParam.pOrg = pcPatternKey->getROIY() + (pcPatternKey->getPatternLStride() << isubShift);
          m_cDistParam.pCur = piRefSrch + (rcStruct.iYStride << isubShift);
          uiTempSad = m_cDistParam.DistFunc( &m_cDistParam );
          uiSad += uiTempSad >>  m_cDistParam.iSubShift;
          if(((uiSad << isubShift) + uiBitCost) > rcStruct.uiBestSad)
          {
            break;
          }

          m_cDistParam.iSubShift--;
        }

        if(m_cDistParam.iSubShift == 0)
        {
          uiSad += uiBitCost;
          if( uiSad < rcStruct.uiBestSad )
          {
            rcStruct.uiBestSad      = uiSad;
            rcStruct.iBestX         = iSearchX;
            rcStruct.iBestY         = iSearchY;
            rcStruct.uiBestDistance = uiDistance;
            rcStruct.uiBestRound    = 0;
            rcStruct.ucPointNr      = ucPointNr;
            m_cDistParam.m_maximumDistortionForEarlyExit = uiSad;
          }
        }
      }
    }
  }
  else
  {
    // fast encoder decision: use subsampled SAD when rows > 8 for integer ME
    if ( m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE3 )
    {
      if ( m_cDistParam.iRows > 8 )
      {
        m_cDistParam.iSubShift = 1;
      }
    }

    uiSad = m_cDistParam.DistFunc( &m_cDistParam );

	// modification here

	array[counter_i] = uiSad;
	//cout << "\narray    " << array[counter_i];
	//cout << "   number   " << counter_i;
	// end of modification
    // only add motion cost if uiSad is smaller than best. Otherwise pointless
    // to add motion cost.
    if( uiSad < rcStruct.uiBestSad )
    {
      // motion cost
      uiSad += m_pcRdCost->getCostOfVectorWithPredictor( iSearchX, iSearchY );

      if( uiSad < rcStruct.uiBestSad )
      {
        rcStruct.uiBestSad      = uiSad;
        rcStruct.iBestX         = iSearchX;
        rcStruct.iBestY         = iSearchY;
        rcStruct.uiBestDistance = uiDistance;
        rcStruct.uiBestRound    = 0;
        rcStruct.ucPointNr      = ucPointNr;
        m_cDistParam.m_maximumDistortionForEarlyExit = uiSad;
      }
    }
  }
  counter_i = counter_i + 1;
}

__inline Void TEncSearch::xTZ2PointSearch( const TComPattern* const pcPatternKey, IntTZSearchStruct& rcStruct, const TComMv* const pcMvSrchRngLT, const TComMv* const pcMvSrchRngRB )
{
  Int   iSrchRngHorLeft   = pcMvSrchRngLT->getHor();
  Int   iSrchRngHorRight  = pcMvSrchRngRB->getHor();
  Int   iSrchRngVerTop    = pcMvSrchRngLT->getVer();
  Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

  // 2 point search,                   //   1 2 3
  // check only the 2 untested points  //   4 0 5
  // around the start point            //   6 7 8
  Int iStartX = rcStruct.iBestX;
  Int iStartY = rcStruct.iBestY;
  switch( rcStruct.ucPointNr )
  {
    case 1:
    {
      if ( (iStartX - 1) >= iSrchRngHorLeft )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX - 1, iStartY, 0, 2 );
      }
      if ( (iStartY - 1) >= iSrchRngVerTop )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iStartY - 1, 0, 2 );
      }
    }
      break;
    case 2:
    {
      if ( (iStartY - 1) >= iSrchRngVerTop )
      {
        if ( (iStartX - 1) >= iSrchRngHorLeft )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX - 1, iStartY - 1, 0, 2 );
        }
        if ( (iStartX + 1) <= iSrchRngHorRight )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX + 1, iStartY - 1, 0, 2 );
        }
      }
    }
      break;
    case 3:
    {
      if ( (iStartY - 1) >= iSrchRngVerTop )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iStartY - 1, 0, 2 );
      }
      if ( (iStartX + 1) <= iSrchRngHorRight )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX + 1, iStartY, 0, 2 );
      }
    }
      break;
    case 4:
    {
      if ( (iStartX - 1) >= iSrchRngHorLeft )
      {
        if ( (iStartY + 1) <= iSrchRngVerBottom )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX - 1, iStartY + 1, 0, 2 );
        }
        if ( (iStartY - 1) >= iSrchRngVerTop )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX - 1, iStartY - 1, 0, 2 );
        }
      }
    }
      break;
    case 5:
    {
      if ( (iStartX + 1) <= iSrchRngHorRight )
      {
        if ( (iStartY - 1) >= iSrchRngVerTop )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX + 1, iStartY - 1, 0, 2 );
        }
        if ( (iStartY + 1) <= iSrchRngVerBottom )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX + 1, iStartY + 1, 0, 2 );
        }
      }
    }
      break;
    case 6:
    {
      if ( (iStartX - 1) >= iSrchRngHorLeft )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX - 1, iStartY , 0, 2 );
      }
      if ( (iStartY + 1) <= iSrchRngVerBottom )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iStartY + 1, 0, 2 );
      }
    }
      break;
    case 7:
    {
      if ( (iStartY + 1) <= iSrchRngVerBottom )
      {
        if ( (iStartX - 1) >= iSrchRngHorLeft )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX - 1, iStartY + 1, 0, 2 );
        }
        if ( (iStartX + 1) <= iSrchRngHorRight )
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX + 1, iStartY + 1, 0, 2 );
        }
      }
    }
      break;
    case 8:
    {
      if ( (iStartX + 1) <= iSrchRngHorRight )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX + 1, iStartY, 0, 2 );
      }
      if ( (iStartY + 1) <= iSrchRngVerBottom )
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iStartY + 1, 0, 2 );
      }
    }
      break;
    default:
    {
      assert( false );
    }
      break;
  } // switch( rcStruct.ucPointNr )
}




__inline Void TEncSearch::xTZ8PointSquareSearch( const TComPattern* const pcPatternKey, IntTZSearchStruct& rcStruct, const TComMv* const pcMvSrchRngLT, const TComMv* const pcMvSrchRngRB, const Int iStartX, const Int iStartY, const Int iDist )
{
  const Int   iSrchRngHorLeft   = pcMvSrchRngLT->getHor();
  const Int   iSrchRngHorRight  = pcMvSrchRngRB->getHor();
  const Int   iSrchRngVerTop    = pcMvSrchRngLT->getVer();
  const Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

  // 8 point search,                   //   1 2 3
  // search around the start point     //   4 0 5
  // with the required  distance       //   6 7 8
  assert( iDist != 0 );
  const Int iTop        = iStartY - iDist;
  const Int iBottom     = iStartY + iDist;
  const Int iLeft       = iStartX - iDist;
  const Int iRight      = iStartX + iDist;
  rcStruct.uiBestRound += 1;

  if ( iTop >= iSrchRngVerTop ) // check top
  {
    if ( iLeft >= iSrchRngHorLeft ) // check top left
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iTop, 1, iDist );
    }
    // top middle
    xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iTop, 2, iDist );

    if ( iRight <= iSrchRngHorRight ) // check top right
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iRight, iTop, 3, iDist );
    }
  } // check top
  if ( iLeft >= iSrchRngHorLeft ) // check middle left
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iStartY, 4, iDist );
  }
  if ( iRight <= iSrchRngHorRight ) // check middle right
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iRight, iStartY, 5, iDist );
  }
  if ( iBottom <= iSrchRngVerBottom ) // check bottom
  {
    if ( iLeft >= iSrchRngHorLeft ) // check bottom left
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iBottom, 6, iDist );
    }
    // check bottom middle
    xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist );

    if ( iRight <= iSrchRngHorRight ) // check bottom right
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iRight, iBottom, 8, iDist );
    }
  } // check bottom
}


//additing other square search

__inline Void TEncSearch::xTZ8PointSquareSearch2( const TComPattern* const pcPatternKey, IntTZSearchStruct& rcStruct, const TComMv* const pcMvSrchRngLT, const TComMv* const pcMvSrchRngRB, const Int iStartX, const Int iStartY, const Int iDist )
{
  const Int   iSrchRngHorLeft   = pcMvSrchRngLT->getHor();
  const Int   iSrchRngHorRight  = pcMvSrchRngRB->getHor();
  const Int   iSrchRngVerTop    = pcMvSrchRngLT->getVer();
  const Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

  // 8 point search,                   //   1 2 3
  // search around the start point     //   4 0 5
  // with the required  distance       //   6 7 8
  assert( iDist != 0 );
  const Int iTop        = iStartY - iDist;
  const Int iBottom     = iStartY + iDist;
  const Int iLeft       = iStartX - iDist;
  const Int iRight      = iStartX + iDist;
  rcStruct.uiBestRound += 1;
// check top
  if ( iTop >= iSrchRngVerTop ) // check top
  {
	 if ( iLeft >= iSrchRngHorLeft ) // check top left
    {
		xTZSearchHelp(pcPatternKey, rcStruct, iLeft, iTop, 9, iDist);
    }
	  
	 if ( iLeft >= iSrchRngHorLeft ) // check top left
    {
		xTZSearchHelp(pcPatternKey, rcStruct, iStartX - 1, iTop, 10, iDist);
    }
    xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iTop, 11, iDist );
	
	if (iRight <= iSrchRngHorRight) // check top left
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iStartX +1, iTop, 12, iDist );
    }
	
    if ( iRight <= iSrchRngHorRight ) // check top right
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iRight, iTop, 13, iDist );
    }
  }

  if ( iLeft >= iSrchRngHorLeft ) // check middle left
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iStartY-1, 14, iDist );
  }

  if (iRight <= iSrchRngHorRight) // check middle left
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iRight, iStartY-1, 15, iDist );
  }
  
  
  if ( iLeft >= iSrchRngHorLeft ) // check middle left
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iStartY, 16, iDist );
  }
  
  
  if ( iRight <= iSrchRngHorRight ) // check middle right
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iRight, iStartY, 17, iDist );
  }
  
  if ( iLeft >= iSrchRngHorLeft ) // check middle left
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iStartY+1, 18, iDist );
  }
  
  if (iRight <= iSrchRngHorRight) // check middle left
  {
    xTZSearchHelp( pcPatternKey, rcStruct, iRight, iStartY+1, 19, iDist );
  }
  
  
  
  if ( iBottom <= iSrchRngVerBottom ) // check bottom
  {
	  
	if ( iLeft >= iSrchRngHorLeft ) // check bottom left
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iBottom, 20, iDist );
    }  
	  
	if ( iLeft >= iSrchRngHorLeft ) // check bottom left
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iStartX - 1, iBottom, 21, iDist );
    }   
	  
	  
    
    // check bottom middle
    xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iBottom, 22, iDist );

	if ( iRight <= iSrchRngHorRight ) // check bottom right
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iStartX + 1, iBottom, 23, iDist );
    }
	
    if ( iRight <= iSrchRngHorRight ) // check bottom right
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iRight, iBottom, 24, iDist );
    }
  } 
  
  // check bottom
}











__inline Void TEncSearch::xTZ8PointDiamondSearch( const TComPattern*const  pcPatternKey,
                                                  IntTZSearchStruct& rcStruct,
                                                  const TComMv*const  pcMvSrchRngLT,
                                                  const TComMv*const  pcMvSrchRngRB,
                                                  const Int iStartX,
                                                  const Int iStartY,
                                                  const Int iDist,
                                                  const Bool bCheckCornersAtDist1 )
{
  const Int   iSrchRngHorLeft   = pcMvSrchRngLT->getHor();
  const Int   iSrchRngHorRight  = pcMvSrchRngRB->getHor();
  const Int   iSrchRngVerTop    = pcMvSrchRngLT->getVer();
  const Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

  // 8 point search,                   //   1 2 3
  // search around the start point     //   4 0 5
  // with the required  distance       //   6 7 8
  assert ( iDist != 0 );
  const Int iTop        = iStartY - iDist;
  const Int iBottom     = iStartY + iDist;
  const Int iLeft       = iStartX - iDist;
  const Int iRight      = iStartX + iDist;
  rcStruct.uiBestRound += 1;

  if ( iDist == 1 )
  {
    if ( iTop >= iSrchRngVerTop ) // check top
    {
      if (bCheckCornersAtDist1)
      {
        if ( iLeft >= iSrchRngHorLeft) // check top-left
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iTop, 1, iDist );
        }
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iTop, 2, iDist );
        if ( iRight <= iSrchRngHorRight ) // check middle right
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iRight, iTop, 3, iDist );
        }
      }
      else
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iTop, 2, iDist );
      }
    }
    if ( iLeft >= iSrchRngHorLeft ) // check middle left
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iStartY, 4, iDist );
    }
    if ( iRight <= iSrchRngHorRight ) // check middle right
    {
      xTZSearchHelp( pcPatternKey, rcStruct, iRight, iStartY, 5, iDist );
    }
    if ( iBottom <= iSrchRngVerBottom ) // check bottom
    {
      if (bCheckCornersAtDist1)
      {
        if ( iLeft >= iSrchRngHorLeft) // check top-left
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iBottom, 6, iDist );
        }
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist );
        if ( iRight <= iSrchRngHorRight ) // check middle right
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iRight, iBottom, 8, iDist );
        }
      }
      else
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist );
      }
    }
  }
  else
  {
    if ( iDist <= 8 )
    {
      const Int iTop_2      = iStartY - (iDist>>1);
      const Int iBottom_2   = iStartY + (iDist>>1);
      const Int iLeft_2     = iStartX - (iDist>>1);
      const Int iRight_2    = iStartX + (iDist>>1);

      if (  iTop >= iSrchRngVerTop && iLeft >= iSrchRngHorLeft &&
          iRight <= iSrchRngHorRight && iBottom <= iSrchRngVerBottom ) // check border
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX,  iTop,      2, iDist    );
        xTZSearchHelp( pcPatternKey, rcStruct, iLeft_2,  iTop_2,    1, iDist>>1 );
        xTZSearchHelp( pcPatternKey, rcStruct, iRight_2, iTop_2,    3, iDist>>1 );
        xTZSearchHelp( pcPatternKey, rcStruct, iLeft,    iStartY,   4, iDist    );
        xTZSearchHelp( pcPatternKey, rcStruct, iRight,   iStartY,   5, iDist    );
        xTZSearchHelp( pcPatternKey, rcStruct, iLeft_2,  iBottom_2, 6, iDist>>1 );
        xTZSearchHelp( pcPatternKey, rcStruct, iRight_2, iBottom_2, 8, iDist>>1 );
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX,  iBottom,   7, iDist    );
      }
      else // check border
      {
        if ( iTop >= iSrchRngVerTop ) // check top
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iTop, 2, iDist );
        }
        if ( iTop_2 >= iSrchRngVerTop ) // check half top
        {
          if ( iLeft_2 >= iSrchRngHorLeft ) // check half left
          {
            xTZSearchHelp( pcPatternKey, rcStruct, iLeft_2, iTop_2, 1, (iDist>>1) );
          }
          if ( iRight_2 <= iSrchRngHorRight ) // check half right
          {
            xTZSearchHelp( pcPatternKey, rcStruct, iRight_2, iTop_2, 3, (iDist>>1) );
          }
        } // check half top
        if ( iLeft >= iSrchRngHorLeft ) // check left
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iStartY, 4, iDist );
        }
        if ( iRight <= iSrchRngHorRight ) // check right
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iRight, iStartY, 5, iDist );
        }
        if ( iBottom_2 <= iSrchRngVerBottom ) // check half bottom
        {
          if ( iLeft_2 >= iSrchRngHorLeft ) // check half left
          {
            xTZSearchHelp( pcPatternKey, rcStruct, iLeft_2, iBottom_2, 6, (iDist>>1) );
          }
          if ( iRight_2 <= iSrchRngHorRight ) // check half right
          {
            xTZSearchHelp( pcPatternKey, rcStruct, iRight_2, iBottom_2, 8, (iDist>>1) );
          }
        } // check half bottom
        if ( iBottom <= iSrchRngVerBottom ) // check bottom
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iBottom, 7, iDist );
        }
      } // check border
    }
    else // iDist > 8
    {
      if ( iTop >= iSrchRngVerTop && iLeft >= iSrchRngHorLeft &&
          iRight <= iSrchRngHorRight && iBottom <= iSrchRngVerBottom ) // check border
      {
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iTop,    0, iDist );
        xTZSearchHelp( pcPatternKey, rcStruct, iLeft,   iStartY, 0, iDist );
        xTZSearchHelp( pcPatternKey, rcStruct, iRight,  iStartY, 0, iDist );
        xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iBottom, 0, iDist );
        for ( Int index = 1; index < 4; index++ )
        {
          const Int iPosYT = iTop    + ((iDist>>2) * index);
          const Int iPosYB = iBottom - ((iDist>>2) * index);
          const Int iPosXL = iStartX - ((iDist>>2) * index);
          const Int iPosXR = iStartX + ((iDist>>2) * index);
          xTZSearchHelp( pcPatternKey, rcStruct, iPosXL, iPosYT, 0, iDist );
          xTZSearchHelp( pcPatternKey, rcStruct, iPosXR, iPosYT, 0, iDist );
          xTZSearchHelp( pcPatternKey, rcStruct, iPosXL, iPosYB, 0, iDist );
          xTZSearchHelp( pcPatternKey, rcStruct, iPosXR, iPosYB, 0, iDist );
        }
      }
      else // check border
      {
        if ( iTop >= iSrchRngVerTop ) // check top
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iTop, 0, iDist );
        }
        if ( iLeft >= iSrchRngHorLeft ) // check left
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iLeft, iStartY, 0, iDist );
        }
        if ( iRight <= iSrchRngHorRight ) // check right
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iRight, iStartY, 0, iDist );
        }
        if ( iBottom <= iSrchRngVerBottom ) // check bottom
        {
          xTZSearchHelp( pcPatternKey, rcStruct, iStartX, iBottom, 0, iDist );
        }
        for ( Int index = 1; index < 4; index++ )
        {
          const Int iPosYT = iTop    + ((iDist>>2) * index);
          const Int iPosYB = iBottom - ((iDist>>2) * index);
          const Int iPosXL = iStartX - ((iDist>>2) * index);
          const Int iPosXR = iStartX + ((iDist>>2) * index);

          if ( iPosYT >= iSrchRngVerTop ) // check top
          {
            if ( iPosXL >= iSrchRngHorLeft ) // check left
            {
              xTZSearchHelp( pcPatternKey, rcStruct, iPosXL, iPosYT, 0, iDist );
            }
            if ( iPosXR <= iSrchRngHorRight ) // check right
            {
              xTZSearchHelp( pcPatternKey, rcStruct, iPosXR, iPosYT, 0, iDist );
            }
          } // check top
          if ( iPosYB <= iSrchRngVerBottom ) // check bottom
          {
            if ( iPosXL >= iSrchRngHorLeft ) // check left
            {
              xTZSearchHelp( pcPatternKey, rcStruct, iPosXL, iPosYB, 0, iDist );
            }
            if ( iPosXR <= iSrchRngHorRight ) // check right
            {
              xTZSearchHelp( pcPatternKey, rcStruct, iPosXR, iPosYB, 0, iDist );
            }
          } // check bottom
        } // for ...
      } // check border
    } // iDist <= 8
  } // iDist == 1
}

Distortion TEncSearch::xPatternRefinement( TComPattern* pcPatternKey,
                                           TComMv baseRefMv,
                                           Int iFrac, TComMv& rcMvFrac,
                                           Bool bAllowUseOfHadamard
                                         )
{
  Distortion  uiDist;
  Distortion  uiDistBest  = std::numeric_limits<Distortion>::max();
  UInt        uiDirecBest = 0;

  Pel*  piRefPos;
  Int iRefStride = m_filteredBlock[0][0].getStride(COMPONENT_Y);

  m_pcRdCost->setDistParam( pcPatternKey, m_filteredBlock[0][0].getAddr(COMPONENT_Y), iRefStride, 1, m_cDistParam, m_pcEncCfg->getUseHADME() && bAllowUseOfHadamard );

  const TComMv* pcMvRefine = (iFrac == 2 ? s_acMvRefineH : s_acMvRefineQ);

  for (UInt i = 0; i < 9; i++)
  {
    TComMv cMvTest = pcMvRefine[i];
    cMvTest += baseRefMv;

    Int horVal = cMvTest.getHor() * iFrac;
    Int verVal = cMvTest.getVer() * iFrac;
    piRefPos = m_filteredBlock[ verVal & 3 ][ horVal & 3 ].getAddr(COMPONENT_Y);
    if ( horVal == 2 && ( verVal & 1 ) == 0 )
    {
      piRefPos += 1;
    }
    if ( ( horVal & 1 ) == 0 && verVal == 2 )
    {
      piRefPos += iRefStride;
    }
    cMvTest = pcMvRefine[i];
    cMvTest += rcMvFrac;

    setDistParamComp(COMPONENT_Y);

    m_cDistParam.pCur = piRefPos;
    m_cDistParam.bitDepth = pcPatternKey->getBitDepthY();
    uiDist = m_cDistParam.DistFunc( &m_cDistParam );
    uiDist += m_pcRdCost->getCostOfVectorWithPredictor( cMvTest.getHor(), cMvTest.getVer() );

    if ( uiDist < uiDistBest )
    {
      uiDistBest  = uiDist;
      uiDirecBest = i;
      m_cDistParam.m_maximumDistortionForEarlyExit = uiDist;
    }
  }

  rcMvFrac = pcMvRefine[uiDirecBest];

  return uiDistBest;
}



Void
TEncSearch::xEncSubdivCbfQT(TComTU      &rTu,
                            Bool         bLuma,
                            Bool         bChroma )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx         = rTu.GetAbsPartIdxTU();
  const UInt uiTrDepth            = rTu.GetTransformDepthRel();
  const UInt uiTrMode             = pcCU->getTransformIdx( uiAbsPartIdx );
  const UInt uiSubdiv             = ( uiTrMode > uiTrDepth ? 1 : 0 );
  const UInt uiLog2LumaTrafoSize  = rTu.GetLog2LumaTrSize();

  if( pcCU->isIntra(0) && pcCU->getPartitionSize(0) == SIZE_NxN && uiTrDepth == 0 )
  {
    assert( uiSubdiv );
  }
  else if( uiLog2LumaTrafoSize > pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() )
  {
    assert( uiSubdiv );
  }
  else if( uiLog2LumaTrafoSize == pcCU->getSlice()->getSPS()->getQuadtreeTULog2MinSize() )
  {
    assert( !uiSubdiv );
  }
  else if( uiLog2LumaTrafoSize == pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) )
  {
    assert( !uiSubdiv );
  }
  else
  {
    assert( uiLog2LumaTrafoSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) );
    if( bLuma )
    {
      m_pcEntropyCoder->encodeTransformSubdivFlag( uiSubdiv, 5 - uiLog2LumaTrafoSize );
    }
  }

  if ( bChroma )
  {
    const UInt numberValidComponents = getNumberValidComponents(rTu.GetChromaFormat());
    for (UInt ch=COMPONENT_Cb; ch<numberValidComponents; ch++)
    {
      const ComponentID compID=ComponentID(ch);
      if( rTu.ProcessingAllQuadrants(compID) && (uiTrDepth==0 || pcCU->getCbf( uiAbsPartIdx, compID, uiTrDepth-1 ) ))
      {
        m_pcEntropyCoder->encodeQtCbf(rTu, compID, (uiSubdiv == 0));
      }
    }
  }

  if( uiSubdiv )
  {
    TComTURecurse tuRecurse(rTu, false);
    do
    {
      xEncSubdivCbfQT( tuRecurse, bLuma, bChroma );
    } while (tuRecurse.nextSection(rTu));
  }
  else
  {
    //===== Cbfs =====
    if( bLuma )
    {
      m_pcEntropyCoder->encodeQtCbf( rTu, COMPONENT_Y, true );
    }
  }
}




Void
TEncSearch::xEncCoeffQT(TComTU &rTu,
                        const ComponentID  component,
                        Bool         bRealCoeff )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  const UInt uiTrDepth=rTu.GetTransformDepthRel();

  const UInt  uiTrMode        = pcCU->getTransformIdx( uiAbsPartIdx );
  const UInt  uiSubdiv        = ( uiTrMode > uiTrDepth ? 1 : 0 );

  if( uiSubdiv )
  {
    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
      xEncCoeffQT( tuRecurseChild, component, bRealCoeff );
    } while (tuRecurseChild.nextSection(rTu) );
  }
  else if (rTu.ProcessComponentSection(component))
  {
    //===== coefficients =====
    const UInt  uiLog2TrafoSize = rTu.GetLog2LumaTrSize();
    UInt    uiCoeffOffset   = rTu.getCoefficientOffset(component);
    UInt    uiQTLayer       = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrafoSize;
    TCoeff* pcCoeff         = bRealCoeff ? pcCU->getCoeff(component) : m_ppcQTTempCoeff[component][uiQTLayer];

    if (isChroma(component) && (pcCU->getCbf( rTu.GetAbsPartIdxTU(), COMPONENT_Y, uiTrMode ) != 0) && pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag() )
    {
      m_pcEntropyCoder->encodeCrossComponentPrediction( rTu, component );
    }

    m_pcEntropyCoder->encodeCoeffNxN( rTu, pcCoeff+uiCoeffOffset, component );
  }
}




Void
TEncSearch::xEncIntraHeader( TComDataCU*  pcCU,
                            UInt         uiTrDepth,
                            UInt         uiAbsPartIdx,
                            Bool         bLuma,
                            Bool         bChroma )
{
  if( bLuma )
  {
    // CU header
    if( uiAbsPartIdx == 0 )
    {
      if( !pcCU->getSlice()->isIntra() )
      {
        if (pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
        {
          m_pcEntropyCoder->encodeCUTransquantBypassFlag( pcCU, 0, true );
        }
        m_pcEntropyCoder->encodeSkipFlag( pcCU, 0, true );
        m_pcEntropyCoder->encodePredMode( pcCU, 0, true );
      }
      m_pcEntropyCoder  ->encodePartSize( pcCU, 0, pcCU->getDepth(0), true );

      if (pcCU->isIntra(0) && pcCU->getPartitionSize(0) == SIZE_2Nx2N )
      {
        m_pcEntropyCoder->encodeIPCMInfo( pcCU, 0, true );

        if ( pcCU->getIPCMFlag (0))
        {
          return;
        }
      }
    }
    // luma prediction mode
    if( pcCU->getPartitionSize(0) == SIZE_2Nx2N )
    {
      if (uiAbsPartIdx==0)
      {
        m_pcEntropyCoder->encodeIntraDirModeLuma ( pcCU, 0 );
      }
    }
    else
    {
      UInt uiQNumParts = pcCU->getTotalNumPart() >> 2;
      if (uiTrDepth>0 && (uiAbsPartIdx%uiQNumParts)==0)
      {
        m_pcEntropyCoder->encodeIntraDirModeLuma ( pcCU, uiAbsPartIdx );
      }
    }
  }

  if( bChroma )
  {
    if( pcCU->getPartitionSize(0) == SIZE_2Nx2N || !enable4ChromaPUsInIntraNxNCU(pcCU->getPic()->getChromaFormat()))
    {
      if(uiAbsPartIdx==0)
      {
         m_pcEntropyCoder->encodeIntraDirModeChroma ( pcCU, uiAbsPartIdx );
      }
    }
    else
    {
      UInt uiQNumParts = pcCU->getTotalNumPart() >> 2;
      assert(uiTrDepth>0);
      if ((uiAbsPartIdx%uiQNumParts)==0)
      {
        m_pcEntropyCoder->encodeIntraDirModeChroma ( pcCU, uiAbsPartIdx );
      }
    }
  }
}




UInt
TEncSearch::xGetIntraBitsQT(TComTU &rTu,
                            Bool         bLuma,
                            Bool         bChroma,
                            Bool         bRealCoeff /* just for test */ )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  const UInt uiTrDepth=rTu.GetTransformDepthRel();
  m_pcEntropyCoder->resetBits();
  xEncIntraHeader ( pcCU, uiTrDepth, uiAbsPartIdx, bLuma, bChroma );
  xEncSubdivCbfQT ( rTu, bLuma, bChroma );

  if( bLuma )
  {
    xEncCoeffQT   ( rTu, COMPONENT_Y,      bRealCoeff );
  }
  if( bChroma )
  {
    xEncCoeffQT   ( rTu, COMPONENT_Cb,  bRealCoeff );
    xEncCoeffQT   ( rTu, COMPONENT_Cr,  bRealCoeff );
  }
  UInt   uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();

  return uiBits;
}

UInt TEncSearch::xGetIntraBitsQTChroma(TComTU &rTu,
                                       ComponentID compID,
                                       Bool         bRealCoeff /* just for test */ )
{
  m_pcEntropyCoder->resetBits();
  xEncCoeffQT   ( rTu, compID,  bRealCoeff );
  UInt   uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();
  return uiBits;
}

Void TEncSearch::xIntraCodingTUBlock(       TComYuv*    pcOrgYuv,
                                            TComYuv*    pcPredYuv,
                                            TComYuv*    pcResiYuv,
                                            Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
                                      const Bool        checkCrossCPrediction,
                                            Distortion& ruiDist,
                                      const ComponentID compID,
                                            TComTU&     rTu
                                      DEBUG_STRING_FN_DECLARE(sDebug)
                                           ,Int         default0Save1Load2
                                     )
{
  if (!rTu.ProcessComponentSection(compID))
  {
    return;
  }
  const Bool           bIsLuma          = isLuma(compID);
  const TComRectangle &rect             = rTu.getRect(compID);
        TComDataCU    *pcCU             = rTu.getCU();
  const UInt           uiAbsPartIdx     = rTu.GetAbsPartIdxTU();
  const TComSPS       &sps              = *(pcCU->getSlice()->getSPS());

  const UInt           uiTrDepth        = rTu.GetTransformDepthRelAdj(compID);
  const UInt           uiFullDepth      = rTu.GetTransformDepthTotal();
  const UInt           uiLog2TrSize     = rTu.GetLog2LumaTrSize();
  const ChromaFormat   chFmt            = pcOrgYuv->getChromaFormat();
  const ChannelType    chType           = toChannelType(compID);
  const Int            bitDepth         = sps.getBitDepth(chType);

  const UInt           uiWidth          = rect.width;
  const UInt           uiHeight         = rect.height;
  const UInt           uiStride         = pcOrgYuv ->getStride (compID);
        Pel           *piOrg            = pcOrgYuv ->getAddr( compID, uiAbsPartIdx );
        Pel           *piPred           = pcPredYuv->getAddr( compID, uiAbsPartIdx );
        Pel           *piResi           = pcResiYuv->getAddr( compID, uiAbsPartIdx );
        Pel           *piReco           = pcPredYuv->getAddr( compID, uiAbsPartIdx );
  const UInt           uiQTLayer        = sps.getQuadtreeTULog2MaxSize() - uiLog2TrSize;
        Pel           *piRecQt          = m_pcQTTempTComYuv[ uiQTLayer ].getAddr( compID, uiAbsPartIdx );
  const UInt           uiRecQtStride    = m_pcQTTempTComYuv[ uiQTLayer ].getStride(compID);
  const UInt           uiZOrder         = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
        Pel           *piRecIPred       = pcCU->getPic()->getPicYuvRec()->getAddr( compID, pcCU->getCtuRsAddr(), uiZOrder );
        UInt           uiRecIPredStride = pcCU->getPic()->getPicYuvRec()->getStride  ( compID );
        TCoeff        *pcCoeff          = m_ppcQTTempCoeff[compID][uiQTLayer] + rTu.getCoefficientOffset(compID);
        Bool           useTransformSkip = pcCU->getTransformSkip(uiAbsPartIdx, compID);

#if ADAPTIVE_QP_SELECTION
        TCoeff        *pcArlCoeff       = m_ppcQTTempArlCoeff[compID][ uiQTLayer ] + rTu.getCoefficientOffset(compID);
#endif

  const UInt           uiChPredMode     = pcCU->getIntraDir( chType, uiAbsPartIdx );
  const UInt           partsPerMinCU    = 1<<(2*(sps.getMaxTotalCUDepth() - sps.getLog2DiffMaxMinCodingBlockSize()));
  const UInt           uiChCodedMode    = (uiChPredMode==DM_CHROMA_IDX && !bIsLuma) ? pcCU->getIntraDir(CHANNEL_TYPE_LUMA, getChromasCorrespondingPULumaIdx(uiAbsPartIdx, chFmt, partsPerMinCU)) : uiChPredMode;
  const UInt           uiChFinalMode    = ((chFmt == CHROMA_422)       && !bIsLuma) ? g_chroma422IntraAngleMappingTable[uiChCodedMode] : uiChCodedMode;

  const Int            blkX                                 = g_auiRasterToPelX[ g_auiZscanToRaster[ uiAbsPartIdx ] ];
  const Int            blkY                                 = g_auiRasterToPelY[ g_auiZscanToRaster[ uiAbsPartIdx ] ];
  const Int            bufferOffset                         = blkX + (blkY * MAX_CU_SIZE);
        Pel  *const    encoderLumaResidual                  = resiLuma[RESIDUAL_ENCODER_SIDE ] + bufferOffset;
        Pel  *const    reconstructedLumaResidual            = resiLuma[RESIDUAL_RECONSTRUCTED] + bufferOffset;
  const Bool           bUseCrossCPrediction                 = isChroma(compID) && (uiChPredMode == DM_CHROMA_IDX) && checkCrossCPrediction;
  const Bool           bUseReconstructedResidualForEstimate = m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate();
        Pel *const     lumaResidualForEstimate              = bUseReconstructedResidualForEstimate ? reconstructedLumaResidual : encoderLumaResidual;

#if DEBUG_STRING
  const Int debugPredModeMask=DebugStringGetPredModeMask(MODE_INTRA);
#endif

  //===== init availability pattern =====
  DEBUG_STRING_NEW(sTemp)

#if !DEBUG_STRING
  if( default0Save1Load2 != 2 )
#endif
  {
    const Bool bUseFilteredPredictions=TComPrediction::filteringIntraReferenceSamples(compID, uiChFinalMode, uiWidth, uiHeight, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag());

    initIntraPatternChType( rTu, compID, bUseFilteredPredictions DEBUG_STRING_PASS_INTO(sDebug) );

    //===== get prediction signal =====
    predIntraAng( compID, uiChFinalMode, piOrg, uiStride, piPred, uiStride, rTu, bUseFilteredPredictions );

    // save prediction
    if( default0Save1Load2 == 1 )
    {
      Pel*  pPred   = piPred;
      Pel*  pPredBuf = m_pSharedPredTransformSkip[compID];
      Int k = 0;
      for( UInt uiY = 0; uiY < uiHeight; uiY++ )
      {
        for( UInt uiX = 0; uiX < uiWidth; uiX++ )
        {
          pPredBuf[ k ++ ] = pPred[ uiX ];
        }
        pPred += uiStride;
      }
    }
  }
#if !DEBUG_STRING
  else
  {
    // load prediction
    Pel*  pPred   = piPred;
    Pel*  pPredBuf = m_pSharedPredTransformSkip[compID];
    Int k = 0;
    for( UInt uiY = 0; uiY < uiHeight; uiY++ )
    {
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        pPred[ uiX ] = pPredBuf[ k ++ ];
      }
      pPred += uiStride;
    }
  }
#endif

  //===== get residual signal =====
  {
    // get residual
    Pel*  pOrg    = piOrg;
    Pel*  pPred   = piPred;
    Pel*  pResi   = piResi;

    for( UInt uiY = 0; uiY < uiHeight; uiY++ )
    {
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        pResi[ uiX ] = pOrg[ uiX ] - pPred[ uiX ];
      }

      pOrg  += uiStride;
      pResi += uiStride;
      pPred += uiStride;
    }
  }

  if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
  {
    if (bUseCrossCPrediction)
    {
      if (xCalcCrossComponentPredictionAlpha( rTu, compID, lumaResidualForEstimate, piResi, uiWidth, uiHeight, MAX_CU_SIZE, uiStride ) == 0)
      {
        return;
      }
      TComTrQuant::crossComponentPrediction ( rTu, compID, reconstructedLumaResidual, piResi, piResi, uiWidth, uiHeight, MAX_CU_SIZE, uiStride, uiStride, false );
    }
    else if (isLuma(compID) && !bUseReconstructedResidualForEstimate)
    {
      xStoreCrossComponentPredictionResult( encoderLumaResidual, piResi, rTu, 0, 0, MAX_CU_SIZE, uiStride );
    }
  }

  //===== transform and quantization =====
  //--- init rate estimation arrays for RDOQ ---
  if( useTransformSkip ? m_pcEncCfg->getUseRDOQTS() : m_pcEncCfg->getUseRDOQ() )
  {
    m_pcEntropyCoder->estimateBit( m_pcTrQuant->m_pcEstBitsSbac, uiWidth, uiHeight, chType );
  }

  //--- transform and quantization ---
  TCoeff uiAbsSum = 0;
  if (bIsLuma)
  {
    pcCU       ->setTrIdxSubParts ( uiTrDepth, uiAbsPartIdx, uiFullDepth );
  }

  const QpParam cQP(*pcCU, compID);

#if RDOQ_CHROMA_LAMBDA
  m_pcTrQuant->selectLambda     (compID);
#endif

  m_pcTrQuant->transformNxN     ( rTu, compID, piResi, uiStride, pcCoeff,
#if ADAPTIVE_QP_SELECTION
    pcArlCoeff,
#endif
    uiAbsSum, cQP
    );

  //--- inverse transform ---

#if DEBUG_STRING
  if ( (uiAbsSum > 0) || (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask) )
#else
  if ( uiAbsSum > 0 )
#endif
  {
    m_pcTrQuant->invTransformNxN ( rTu, compID, piResi, uiStride, pcCoeff, cQP DEBUG_STRING_PASS_INTO_OPTIONAL(&sDebug, (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask)) );
  }
  else
  {
    Pel* pResi = piResi;
    memset( pcCoeff, 0, sizeof( TCoeff ) * uiWidth * uiHeight );
    for( UInt uiY = 0; uiY < uiHeight; uiY++ )
    {
      memset( pResi, 0, sizeof( Pel ) * uiWidth );
      pResi += uiStride;
    }
  }


  //===== reconstruction =====
  {
    Pel* pPred      = piPred;
    Pel* pResi      = piResi;
    Pel* pReco      = piReco;
    Pel* pRecQt     = piRecQt;
    Pel* pRecIPred  = piRecIPred;

    if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
    {
      if (bUseCrossCPrediction)
      {
        TComTrQuant::crossComponentPrediction( rTu, compID, reconstructedLumaResidual, piResi, piResi, uiWidth, uiHeight, MAX_CU_SIZE, uiStride, uiStride, true );
      }
      else if (isLuma(compID))
      {
        xStoreCrossComponentPredictionResult( reconstructedLumaResidual, piResi, rTu, 0, 0, MAX_CU_SIZE, uiStride );
      }
    }

 #if DEBUG_STRING
    std::stringstream ss(stringstream::out);
    const Bool bDebugPred=((DebugOptionList::DebugString_Pred.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
    const Bool bDebugResi=((DebugOptionList::DebugString_Resi.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));
    const Bool bDebugReco=((DebugOptionList::DebugString_Reco.getInt()&debugPredModeMask) && DEBUG_STRING_CHANNEL_CONDITION(compID));

    if (bDebugPred || bDebugResi || bDebugReco)
    {
      ss << "###: " << "CompID: " << compID << " pred mode (ch/fin): " << uiChPredMode << "/" << uiChFinalMode << " absPartIdx: " << rTu.GetAbsPartIdxTU() << "\n";
      for( UInt uiY = 0; uiY < uiHeight; uiY++ )
      {
        ss << "###: ";
        if (bDebugPred)
        {
          ss << " - pred: ";
          for( UInt uiX = 0; uiX < uiWidth; uiX++ )
          {
            ss << pPred[ uiX ] << ", ";
          }
        }
        if (bDebugResi)
        {
          ss << " - resi: ";
        }
        for( UInt uiX = 0; uiX < uiWidth; uiX++ )
        {
          if (bDebugResi)
          {
            ss << pResi[ uiX ] << ", ";
          }
          pReco    [ uiX ] = Pel(ClipBD<Int>( Int(pPred[uiX]) + Int(pResi[uiX]), bitDepth ));
          pRecQt   [ uiX ] = pReco[ uiX ];
          pRecIPred[ uiX ] = pReco[ uiX ];
        }
        if (bDebugReco)
        {
          ss << " - reco: ";
          for( UInt uiX = 0; uiX < uiWidth; uiX++ )
          {
            ss << pReco[ uiX ] << ", ";
          }
        }
        pPred     += uiStride;
        pResi     += uiStride;
        pReco     += uiStride;
        pRecQt    += uiRecQtStride;
        pRecIPred += uiRecIPredStride;
        ss << "\n";
      }
      DEBUG_STRING_APPEND(sDebug, ss.str())
    }
    else
#endif
    {

      for( UInt uiY = 0; uiY < uiHeight; uiY++ )
      {
        for( UInt uiX = 0; uiX < uiWidth; uiX++ )
        {
          pReco    [ uiX ] = Pel(ClipBD<Int>( Int(pPred[uiX]) + Int(pResi[uiX]), bitDepth ));
          pRecQt   [ uiX ] = pReco[ uiX ];
          pRecIPred[ uiX ] = pReco[ uiX ];
        }
        pPred     += uiStride;
        pResi     += uiStride;
        pReco     += uiStride;
        pRecQt    += uiRecQtStride;
        pRecIPred += uiRecIPredStride;
      }
    }
  }

  //===== update distortion =====
  ruiDist += m_pcRdCost->getDistPart( bitDepth, piReco, uiStride, piOrg, uiStride, uiWidth, uiHeight, compID );
}




Void
TEncSearch::xRecurIntraCodingLumaQT(TComYuv*    pcOrgYuv,
                                    TComYuv*    pcPredYuv,
                                    TComYuv*    pcResiYuv,
                                    Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
                                    Distortion& ruiDistY,
#if HHI_RQT_INTRA_SPEEDUP
                                    Bool        bCheckFirst,
#endif
                                    Double&     dRDCost,
                                    TComTU&     rTu
                                    DEBUG_STRING_FN_DECLARE(sDebug))
{
  TComDataCU   *pcCU          = rTu.getCU();
  const UInt    uiAbsPartIdx  = rTu.GetAbsPartIdxTU();
  const UInt    uiFullDepth   = rTu.GetTransformDepthTotal();
  const UInt    uiTrDepth     = rTu.GetTransformDepthRel();
  const UInt    uiLog2TrSize  = rTu.GetLog2LumaTrSize();
        Bool    bCheckFull    = ( uiLog2TrSize  <= pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() );
        Bool    bCheckSplit   = ( uiLog2TrSize  >  pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) );

        Pel     resiLumaSplit [NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];
        Pel     resiLumaSingle[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];

        Bool    bMaintainResidual[NUMBER_OF_STORED_RESIDUAL_TYPES];
        for (UInt residualTypeIndex = 0; residualTypeIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; residualTypeIndex++)
        {
          bMaintainResidual[residualTypeIndex] = true; //assume true unless specified otherwise
        }

        bMaintainResidual[RESIDUAL_ENCODER_SIDE] = !(m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate());

#if HHI_RQT_INTRA_SPEEDUP
  Int maxTuSize = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize();
  Int isIntraSlice = (pcCU->getSlice()->getSliceType() == I_SLICE);
  // don't check split if TU size is less or equal to max TU size
  Bool noSplitIntraMaxTuSize = bCheckFull;
  if(m_pcEncCfg->getRDpenalty() && ! isIntraSlice)
  {
    // in addition don't check split if TU size is less or equal to 16x16 TU size for non-intra slice
    noSplitIntraMaxTuSize = ( uiLog2TrSize  <= min(maxTuSize,4) );

    // if maximum RD-penalty don't check TU size 32x32
    if(m_pcEncCfg->getRDpenalty()==2)
    {
      bCheckFull    = ( uiLog2TrSize  <= min(maxTuSize,4));
    }
  }
  if( bCheckFirst && noSplitIntraMaxTuSize )

  {
    bCheckSplit = false;
  }
#else
  Int maxTuSize = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize();
  Int isIntraSlice = (pcCU->getSlice()->getSliceType() == I_SLICE);
  // if maximum RD-penalty don't check TU size 32x32
  if((m_pcEncCfg->getRDpenalty()==2)  && !isIntraSlice)
  {
    bCheckFull    = ( uiLog2TrSize  <= min(maxTuSize,4));
  }
#endif
  Double     dSingleCost                        = MAX_DOUBLE;
  Distortion uiSingleDistLuma                   = 0;
  UInt       uiSingleCbfLuma                    = 0;
  Bool       checkTransformSkip  = pcCU->getSlice()->getPPS()->getUseTransformSkip();
  Int        bestModeId[MAX_NUM_COMPONENT] = { 0, 0, 0};
  checkTransformSkip           &= TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(COMPONENT_Y), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());
  checkTransformSkip           &= (!pcCU->getCUTransquantBypass(0));

  assert (rTu.ProcessComponentSection(COMPONENT_Y));
  const UInt totalAdjustedDepthChan   = rTu.GetTransformDepthTotalAdj(COMPONENT_Y);

  if ( m_pcEncCfg->getUseTransformSkipFast() )
  {
    checkTransformSkip       &= (pcCU->getPartitionSize(uiAbsPartIdx)==SIZE_NxN);
  }

  if( bCheckFull )
  {
    if(checkTransformSkip == true)
    {
      //----- store original entropy coding status -----
      m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_ROOT ] );

      Distortion singleDistTmpLuma                    = 0;
      UInt       singleCbfTmpLuma                     = 0;
      Double     singleCostTmp                        = 0;
      Int        firstCheckId                         = 0;

      for(Int modeId = firstCheckId; modeId < 2; modeId ++)
      {
        DEBUG_STRING_NEW(sModeString)
        Int  default0Save1Load2 = 0;
        singleDistTmpLuma=0;
        if(modeId == firstCheckId)
        {
          default0Save1Load2 = 1;
        }
        else
        {
          default0Save1Load2 = 2;
        }


        pcCU->setTransformSkipSubParts ( modeId, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan );
        xIntraCodingTUBlock( pcOrgYuv, pcPredYuv, pcResiYuv, resiLumaSingle, false, singleDistTmpLuma, COMPONENT_Y, rTu DEBUG_STRING_PASS_INTO(sModeString), default0Save1Load2 );

        singleCbfTmpLuma = pcCU->getCbf( uiAbsPartIdx, COMPONENT_Y, uiTrDepth );

        //----- determine rate and r-d cost -----
        if(modeId == 1 && singleCbfTmpLuma == 0)
        {
          //In order not to code TS flag when cbf is zero, the case for TS with cbf being zero is forbidden.
          singleCostTmp = MAX_DOUBLE;
        }
        else
        {
          UInt uiSingleBits = xGetIntraBitsQT( rTu, true, false, false );
          singleCostTmp     = m_pcRdCost->calcRdCost( uiSingleBits, singleDistTmpLuma );
        }
        if(singleCostTmp < dSingleCost)
        {
          DEBUG_STRING_SWAP(sDebug, sModeString)
          dSingleCost   = singleCostTmp;
          uiSingleDistLuma = singleDistTmpLuma;
          uiSingleCbfLuma = singleCbfTmpLuma;

          bestModeId[COMPONENT_Y] = modeId;
          if(bestModeId[COMPONENT_Y] == firstCheckId)
          {
            xStoreIntraResultQT(COMPONENT_Y, rTu );
            m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiFullDepth ][ CI_TEMP_BEST ] );
          }

          if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
          {
            const Int xOffset = rTu.getRect( COMPONENT_Y ).x0;
            const Int yOffset = rTu.getRect( COMPONENT_Y ).y0;
            for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
            {
              if (bMaintainResidual[storedResidualIndex])
              {
                xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex], resiLumaSingle[storedResidualIndex], rTu, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
              }
            }
          }
        }
        if (modeId == firstCheckId)
        {
          m_pcRDGoOnSbacCoder->load ( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_ROOT ] );
        }
      }

      pcCU ->setTransformSkipSubParts ( bestModeId[COMPONENT_Y], COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan );

      if(bestModeId[COMPONENT_Y] == firstCheckId)
      {
        xLoadIntraResultQT(COMPONENT_Y, rTu );
        pcCU->setCbfSubParts  ( uiSingleCbfLuma << uiTrDepth, COMPONENT_Y, uiAbsPartIdx, rTu.GetTransformDepthTotalAdj(COMPONENT_Y) );

        m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ uiFullDepth ][ CI_TEMP_BEST ] );
      }
    }
    else
    {
      //----- store original entropy coding status -----
      if( bCheckSplit )
      {
        m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_ROOT ] );
      }
      //----- code luma/chroma block with given intra prediction mode and store Cbf-----
      dSingleCost   = 0.0;

      pcCU ->setTransformSkipSubParts ( 0, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan );
      xIntraCodingTUBlock( pcOrgYuv, pcPredYuv, pcResiYuv, resiLumaSingle, false, uiSingleDistLuma, COMPONENT_Y, rTu DEBUG_STRING_PASS_INTO(sDebug));

      if( bCheckSplit )
      {
        uiSingleCbfLuma = pcCU->getCbf( uiAbsPartIdx, COMPONENT_Y, uiTrDepth );
      }
      //----- determine rate and r-d cost -----
      UInt uiSingleBits = xGetIntraBitsQT( rTu, true, false, false );

      if(m_pcEncCfg->getRDpenalty() && (uiLog2TrSize==5) && !isIntraSlice)
      {
        uiSingleBits=uiSingleBits*4;
      }

      dSingleCost       = m_pcRdCost->calcRdCost( uiSingleBits, uiSingleDistLuma );

      if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
      {
        const Int xOffset = rTu.getRect( COMPONENT_Y ).x0;
        const Int yOffset = rTu.getRect( COMPONENT_Y ).y0;
        for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
        {
          if (bMaintainResidual[storedResidualIndex])
          {
            xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex], resiLumaSingle[storedResidualIndex], rTu, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
          }
        }
      }
    }
  }

  if( bCheckSplit )
  {
    //----- store full entropy coding status, load original entropy coding status -----
    if( bCheckFull )
    {
      m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_TEST ] );
      m_pcRDGoOnSbacCoder->load ( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_ROOT ] );
    }
    else
    {
      m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_ROOT ] );
    }
    //----- code splitted block -----
    Double     dSplitCost      = 0.0;
    Distortion uiSplitDistLuma = 0;
    UInt       uiSplitCbfLuma  = 0;

    TComTURecurse tuRecurseChild(rTu, false);
    DEBUG_STRING_NEW(sSplit)
    do
    {
      DEBUG_STRING_NEW(sChild)
#if HHI_RQT_INTRA_SPEEDUP
      xRecurIntraCodingLumaQT( pcOrgYuv, pcPredYuv, pcResiYuv, resiLumaSplit, uiSplitDistLuma, bCheckFirst, dSplitCost, tuRecurseChild DEBUG_STRING_PASS_INTO(sChild) );
#else
      xRecurIntraCodingLumaQT( pcOrgYuv, pcPredYuv, pcResiYuv, resiLumaSplit, uiSplitDistLuma, dSplitCost, tuRecurseChild DEBUG_STRING_PASS_INTO(sChild) );
#endif
      DEBUG_STRING_APPEND(sSplit, sChild)
      uiSplitCbfLuma |= pcCU->getCbf( tuRecurseChild.GetAbsPartIdxTU(), COMPONENT_Y, tuRecurseChild.GetTransformDepthRel() );
    } while (tuRecurseChild.nextSection(rTu) );

    UInt    uiPartsDiv     = rTu.GetAbsPartIdxNumParts();
    {
      if (uiSplitCbfLuma)
      {
        const UInt flag=1<<uiTrDepth;
        UChar *pBase=pcCU->getCbf( COMPONENT_Y );
        for( UInt uiOffs = 0; uiOffs < uiPartsDiv; uiOffs++ )
        {
          pBase[ uiAbsPartIdx + uiOffs ] |= flag;
        }
      }
    }
    //----- restore context states -----
    m_pcRDGoOnSbacCoder->load ( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_ROOT ] );
    
    //----- determine rate and r-d cost -----
    UInt uiSplitBits = xGetIntraBitsQT( rTu, true, false, false );
    dSplitCost       = m_pcRdCost->calcRdCost( uiSplitBits, uiSplitDistLuma );

    //===== compare and set best =====
    if( dSplitCost < dSingleCost )
    {
      //--- update cost ---
      DEBUG_STRING_SWAP(sSplit, sDebug)
      ruiDistY += uiSplitDistLuma;
      dRDCost  += dSplitCost;

      if (pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
      {
        const Int xOffset = rTu.getRect( COMPONENT_Y ).x0;
        const Int yOffset = rTu.getRect( COMPONENT_Y ).y0;
        for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
        {
          if (bMaintainResidual[storedResidualIndex])
          {
            xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex], resiLumaSplit[storedResidualIndex], rTu, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE);
          }
        }
      }

      return;
    }

    //----- set entropy coding status -----
    m_pcRDGoOnSbacCoder->load ( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_TEST ] );

    //--- set transform index and Cbf values ---
    pcCU->setTrIdxSubParts( uiTrDepth, uiAbsPartIdx, uiFullDepth );
    const TComRectangle &tuRect=rTu.getRect(COMPONENT_Y);
    pcCU->setCbfSubParts  ( uiSingleCbfLuma << uiTrDepth, COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan );
    pcCU ->setTransformSkipSubParts  ( bestModeId[COMPONENT_Y], COMPONENT_Y, uiAbsPartIdx, totalAdjustedDepthChan );

    //--- set reconstruction for next intra prediction blocks ---
    const UInt  uiQTLayer   = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
    const UInt  uiZOrder    = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;
    const UInt  uiWidth     = tuRect.width;
    const UInt  uiHeight    = tuRect.height;
    Pel*  piSrc       = m_pcQTTempTComYuv[ uiQTLayer ].getAddr( COMPONENT_Y, uiAbsPartIdx );
    UInt  uiSrcStride = m_pcQTTempTComYuv[ uiQTLayer ].getStride  ( COMPONENT_Y );
    Pel*  piDes       = pcCU->getPic()->getPicYuvRec()->getAddr( COMPONENT_Y, pcCU->getCtuRsAddr(), uiZOrder );
    UInt  uiDesStride = pcCU->getPic()->getPicYuvRec()->getStride  ( COMPONENT_Y );

    for( UInt uiY = 0; uiY < uiHeight; uiY++, piSrc += uiSrcStride, piDes += uiDesStride )
    {
      for( UInt uiX = 0; uiX < uiWidth; uiX++ )
      {
        piDes[ uiX ] = piSrc[ uiX ];
      }
    }
  }
  ruiDistY += uiSingleDistLuma;
  dRDCost  += dSingleCost;
}


Void
TEncSearch::xSetIntraResultLumaQT(TComYuv* pcRecoYuv, TComTU &rTu)
{
  TComDataCU *pcCU        = rTu.getCU();
  const UInt uiTrDepth    = rTu.GetTransformDepthRel();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  UInt uiTrMode     = pcCU->getTransformIdx( uiAbsPartIdx );
  if(  uiTrMode == uiTrDepth )
  {
    UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
    UInt uiQTLayer    = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;

    //===== copy transform coefficients =====

    const TComRectangle &tuRect=rTu.getRect(COMPONENT_Y);
    const UInt coeffOffset = rTu.getCoefficientOffset(COMPONENT_Y);
    const UInt numCoeffInBlock = tuRect.width * tuRect.height;

    if (numCoeffInBlock!=0)
    {
      const TCoeff* srcCoeff = m_ppcQTTempCoeff[COMPONENT_Y][uiQTLayer] + coeffOffset;
      TCoeff* destCoeff      = pcCU->getCoeff(COMPONENT_Y) + coeffOffset;
      ::memcpy( destCoeff, srcCoeff, sizeof(TCoeff)*numCoeffInBlock );
#if ADAPTIVE_QP_SELECTION
      const TCoeff* srcArlCoeff = m_ppcQTTempArlCoeff[COMPONENT_Y][ uiQTLayer ] + coeffOffset;
      TCoeff* destArlCoeff      = pcCU->getArlCoeff (COMPONENT_Y)               + coeffOffset;
      ::memcpy( destArlCoeff, srcArlCoeff, sizeof( TCoeff ) * numCoeffInBlock );
#endif
      m_pcQTTempTComYuv[ uiQTLayer ].copyPartToPartComponent( COMPONENT_Y, pcRecoYuv, uiAbsPartIdx, tuRect.width, tuRect.height );
    }

  }
  else
  {
    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
      xSetIntraResultLumaQT( pcRecoYuv, tuRecurseChild );
    } while (tuRecurseChild.nextSection(rTu));
  }
}


Void
TEncSearch::xStoreIntraResultQT(const ComponentID compID, TComTU &rTu )
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiTrDepth = rTu.GetTransformDepthRel();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  const UInt uiTrMode     = pcCU->getTransformIdx( uiAbsPartIdx );
  if ( compID==COMPONENT_Y || uiTrMode == uiTrDepth )
  {
    assert(uiTrMode == uiTrDepth);
    const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
    const UInt uiQTLayer    = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;

    if (rTu.ProcessComponentSection(compID))
    {
      const TComRectangle &tuRect=rTu.getRect(compID);

      //===== copy transform coefficients =====
      const UInt uiNumCoeff    = tuRect.width * tuRect.height;
      TCoeff* pcCoeffSrc = m_ppcQTTempCoeff[compID] [ uiQTLayer ] + rTu.getCoefficientOffset(compID);
      TCoeff* pcCoeffDst = m_pcQTTempTUCoeff[compID];

      ::memcpy( pcCoeffDst, pcCoeffSrc, sizeof( TCoeff ) * uiNumCoeff );
#if ADAPTIVE_QP_SELECTION
      TCoeff* pcArlCoeffSrc = m_ppcQTTempArlCoeff[compID] [ uiQTLayer ] + rTu.getCoefficientOffset(compID);
      TCoeff* pcArlCoeffDst = m_ppcQTTempTUArlCoeff[compID];
      ::memcpy( pcArlCoeffDst, pcArlCoeffSrc, sizeof( TCoeff ) * uiNumCoeff );
#endif
      //===== copy reconstruction =====
      m_pcQTTempTComYuv[ uiQTLayer ].copyPartToPartComponent( compID, &m_pcQTTempTransformSkipTComYuv, uiAbsPartIdx, tuRect.width, tuRect.height );
    }
  }
}


Void
TEncSearch::xLoadIntraResultQT(const ComponentID compID, TComTU &rTu)
{
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiTrDepth = rTu.GetTransformDepthRel();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  const UInt uiTrMode     = pcCU->getTransformIdx( uiAbsPartIdx );
  if ( compID==COMPONENT_Y || uiTrMode == uiTrDepth )
  {
    assert(uiTrMode == uiTrDepth);
    const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
    const UInt uiQTLayer    = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
    const UInt uiZOrder     = pcCU->getZorderIdxInCtu() + uiAbsPartIdx;

    if (rTu.ProcessComponentSection(compID))
    {
      const TComRectangle &tuRect=rTu.getRect(compID);

      //===== copy transform coefficients =====
      const UInt uiNumCoeff = tuRect.width * tuRect.height;
      TCoeff* pcCoeffDst = m_ppcQTTempCoeff[compID] [ uiQTLayer ] + rTu.getCoefficientOffset(compID);
      TCoeff* pcCoeffSrc = m_pcQTTempTUCoeff[compID];

      ::memcpy( pcCoeffDst, pcCoeffSrc, sizeof( TCoeff ) * uiNumCoeff );
#if ADAPTIVE_QP_SELECTION
      TCoeff* pcArlCoeffDst = m_ppcQTTempArlCoeff[compID] [ uiQTLayer ] + rTu.getCoefficientOffset(compID);
      TCoeff* pcArlCoeffSrc = m_ppcQTTempTUArlCoeff[compID];
      ::memcpy( pcArlCoeffDst, pcArlCoeffSrc, sizeof( TCoeff ) * uiNumCoeff );
#endif
      //===== copy reconstruction =====
      m_pcQTTempTransformSkipTComYuv.copyPartToPartComponent( compID, &m_pcQTTempTComYuv[ uiQTLayer ], uiAbsPartIdx, tuRect.width, tuRect.height );

      Pel*    piRecIPred        = pcCU->getPic()->getPicYuvRec()->getAddr( compID, pcCU->getCtuRsAddr(), uiZOrder );
      UInt    uiRecIPredStride  = pcCU->getPic()->getPicYuvRec()->getStride (compID);
      Pel*    piRecQt           = m_pcQTTempTComYuv[ uiQTLayer ].getAddr( compID, uiAbsPartIdx );
      UInt    uiRecQtStride     = m_pcQTTempTComYuv[ uiQTLayer ].getStride  (compID);
      UInt    uiWidth           = tuRect.width;
      UInt    uiHeight          = tuRect.height;
      Pel* pRecQt               = piRecQt;
      Pel* pRecIPred            = piRecIPred;
      for( UInt uiY = 0; uiY < uiHeight; uiY++ )
      {
        for( UInt uiX = 0; uiX < uiWidth; uiX++ )
        {
          pRecIPred[ uiX ] = pRecQt   [ uiX ];
        }
        pRecQt    += uiRecQtStride;
        pRecIPred += uiRecIPredStride;
      }
    }
  }
}

Void
TEncSearch::xStoreCrossComponentPredictionResult(       Pel    *pResiDst,
                                                  const Pel    *pResiSrc,
                                                        TComTU &rTu,
                                                  const Int     xOffset,
                                                  const Int     yOffset,
                                                  const Int     strideDst,
                                                  const Int     strideSrc )
{
  const Pel *pSrc = pResiSrc + yOffset * strideSrc + xOffset;
        Pel *pDst = pResiDst + yOffset * strideDst + xOffset;

  for( Int y = 0; y < rTu.getRect( COMPONENT_Y ).height; y++ )
  {
    ::memcpy( pDst, pSrc, sizeof(Pel) * rTu.getRect( COMPONENT_Y ).width );
    pDst += strideDst;
    pSrc += strideSrc;
  }
}

SChar
TEncSearch::xCalcCrossComponentPredictionAlpha(       TComTU &rTu,
                                                const ComponentID compID,
                                                const Pel*        piResiL,
                                                const Pel*        piResiC,
                                                const Int         width,
                                                const Int         height,
                                                const Int         strideL,
                                                const Int         strideC )
{
  const Pel *pResiL = piResiL;
  const Pel *pResiC = piResiC;

        TComDataCU *pCU = rTu.getCU();
  const Int  absPartIdx = rTu.GetAbsPartIdxTU( compID );
  const Int diffBitDepth = pCU->getSlice()->getSPS()->getDifferentialLumaChromaBitDepth();

  SChar alpha = 0;
  Int SSxy  = 0;
  Int SSxx  = 0;

  for( UInt uiY = 0; uiY < height; uiY++ )
  {
    for( UInt uiX = 0; uiX < width; uiX++ )
    {
      const Pel scaledResiL = rightShift( pResiL[ uiX ], diffBitDepth );
      SSxy += ( scaledResiL * pResiC[ uiX ] );
      SSxx += ( scaledResiL * scaledResiL   );
    }

    pResiL += strideL;
    pResiC += strideC;
  }

  if( SSxx != 0 )
  {
    Double dAlpha = SSxy / Double( SSxx );
    alpha = SChar(Clip3<Int>(-16, 16, (Int)(dAlpha * 16)));

    static const SChar alphaQuant[17] = {0, 1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8};

    alpha = (alpha < 0) ? -alphaQuant[Int(-alpha)] : alphaQuant[Int(alpha)];
  }
  pCU->setCrossComponentPredictionAlphaPartRange( alpha, compID, absPartIdx, rTu.GetAbsPartIdxNumParts( compID ) );

  return alpha;
}

Void
TEncSearch::xRecurIntraChromaCodingQT(TComYuv*    pcOrgYuv,
                                      TComYuv*    pcPredYuv,
                                      TComYuv*    pcResiYuv,
                                      Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE],
                                      Distortion& ruiDist,
                                      TComTU&     rTu
                                      DEBUG_STRING_FN_DECLARE(sDebug))
{
  TComDataCU         *pcCU                  = rTu.getCU();
  const UInt          uiTrDepth             = rTu.GetTransformDepthRel();
  const UInt          uiAbsPartIdx          = rTu.GetAbsPartIdxTU();
  const ChromaFormat  format                = rTu.GetChromaFormat();
  UInt                uiTrMode              = pcCU->getTransformIdx( uiAbsPartIdx );
  const UInt          numberValidComponents = getNumberValidComponents(format);

  if(  uiTrMode == uiTrDepth )
  {
    if (!rTu.ProcessChannelSection(CHANNEL_TYPE_CHROMA))
    {
      return;
    }

    const UInt uiFullDepth = rTu.GetTransformDepthTotal();

    Bool checkTransformSkip = pcCU->getSlice()->getPPS()->getUseTransformSkip();
    checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(COMPONENT_Cb), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());

    if ( m_pcEncCfg->getUseTransformSkipFast() )
    {
      checkTransformSkip &= TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(COMPONENT_Y), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize());

      if (checkTransformSkip)
      {
        Int nbLumaSkip = 0;
        const UInt maxAbsPartIdxSub=uiAbsPartIdx + (rTu.ProcessingAllQuadrants(COMPONENT_Cb)?1:4);
        for(UInt absPartIdxSub = uiAbsPartIdx; absPartIdxSub < maxAbsPartIdxSub; absPartIdxSub ++)
        {
          nbLumaSkip += pcCU->getTransformSkip(absPartIdxSub, COMPONENT_Y);
        }
        checkTransformSkip &= (nbLumaSkip > 0);
      }
    }


    for (UInt ch=COMPONENT_Cb; ch<numberValidComponents; ch++)
    {
      const ComponentID compID = ComponentID(ch);
      DEBUG_STRING_NEW(sDebugBestMode)

      //use RDO to decide whether Cr/Cb takes TS
      m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[uiFullDepth][CI_QT_TRAFO_ROOT] );

      const Bool splitIntoSubTUs = rTu.getRect(compID).width != rTu.getRect(compID).height;

      TComTURecurse TUIterator(rTu, false, (splitIntoSubTUs ? TComTU::VERTICAL_SPLIT : TComTU::DONT_SPLIT), true, compID);

      const UInt partIdxesPerSubTU = TUIterator.GetAbsPartIdxNumParts(compID);

      do
      {
        const UInt subTUAbsPartIdx   = TUIterator.GetAbsPartIdxTU(compID);

        Double     dSingleCost               = MAX_DOUBLE;
        Int        bestModeId                = 0;
        Distortion singleDistC               = 0;
        UInt       singleCbfC                = 0;
        Distortion singleDistCTmp            = 0;
        Double     singleCostTmp             = 0;
        UInt       singleCbfCTmp             = 0;
        SChar      bestCrossCPredictionAlpha = 0;
        Int        bestTransformSkipMode     = 0;

        const Bool checkCrossComponentPrediction =    (pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, subTUAbsPartIdx) == DM_CHROMA_IDX)
                                                   &&  pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag()
                                                   && (pcCU->getCbf(subTUAbsPartIdx,  COMPONENT_Y, uiTrDepth) != 0);

        const Int  crossCPredictionModesToTest = checkCrossComponentPrediction ? 2 : 1;
        const Int  transformSkipModesToTest    = checkTransformSkip            ? 2 : 1;
        const Int  totalModesToTest            = crossCPredictionModesToTest * transformSkipModesToTest;
              Int  currModeId                  = 0;
              Int  default0Save1Load2          = 0;

        for(Int transformSkipModeId = 0; transformSkipModeId < transformSkipModesToTest; transformSkipModeId++)
        {
          for(Int crossCPredictionModeId = 0; crossCPredictionModeId < crossCPredictionModesToTest; crossCPredictionModeId++)
          {
            pcCU->setCrossComponentPredictionAlphaPartRange(0, compID, subTUAbsPartIdx, partIdxesPerSubTU);
            DEBUG_STRING_NEW(sDebugMode)
            pcCU->setTransformSkipPartRange( transformSkipModeId, compID, subTUAbsPartIdx, partIdxesPerSubTU );
            currModeId++;

            const Bool isOneMode  = (totalModesToTest == 1);
            const Bool isLastMode = (currModeId == totalModesToTest); // currModeId is indexed from 1

            if (isOneMode)
            {
              default0Save1Load2 = 0;
            }
            else if (!isOneMode && (transformSkipModeId == 0) && (crossCPredictionModeId == 0))
            {
              default0Save1Load2 = 1; //save prediction on first mode
            }
            else
            {
              default0Save1Load2 = 2; //load it on subsequent modes
            }

            singleDistCTmp = 0;

            xIntraCodingTUBlock( pcOrgYuv, pcPredYuv, pcResiYuv, resiLuma, (crossCPredictionModeId != 0), singleDistCTmp, compID, TUIterator DEBUG_STRING_PASS_INTO(sDebugMode), default0Save1Load2);
            singleCbfCTmp = pcCU->getCbf( subTUAbsPartIdx, compID, uiTrDepth);

            if (  ((crossCPredictionModeId == 1) && (pcCU->getCrossComponentPredictionAlpha(subTUAbsPartIdx, compID) == 0))
               || ((transformSkipModeId    == 1) && (singleCbfCTmp == 0))) //In order not to code TS flag when cbf is zero, the case for TS with cbf being zero is forbidden.
            {
              singleCostTmp = MAX_DOUBLE;
            }
            else if (!isOneMode)
            {
              UInt bitsTmp = xGetIntraBitsQTChroma( TUIterator, compID, false );
              singleCostTmp  = m_pcRdCost->calcRdCost( bitsTmp, singleDistCTmp);
            }

            if(singleCostTmp < dSingleCost)
            {
              DEBUG_STRING_SWAP(sDebugBestMode, sDebugMode)
              dSingleCost               = singleCostTmp;
              singleDistC               = singleDistCTmp;
              bestCrossCPredictionAlpha = (crossCPredictionModeId != 0) ? pcCU->getCrossComponentPredictionAlpha(subTUAbsPartIdx, compID) : 0;
              bestTransformSkipMode     = transformSkipModeId;
              bestModeId                = currModeId;
              singleCbfC                = singleCbfCTmp;

              if (!isOneMode && !isLastMode)
              {
                xStoreIntraResultQT(compID, TUIterator);
                m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiFullDepth ][ CI_TEMP_BEST ] );
              }
            }

            if (!isOneMode && !isLastMode)
            {
              m_pcRDGoOnSbacCoder->load ( m_pppcRDSbacCoder[ uiFullDepth ][ CI_QT_TRAFO_ROOT ] );
            }
          }
        }

        if(bestModeId < totalModesToTest)
        {
          xLoadIntraResultQT(compID, TUIterator);
          pcCU->setCbfPartRange( singleCbfC << uiTrDepth, compID, subTUAbsPartIdx, partIdxesPerSubTU );

          m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ uiFullDepth ][ CI_TEMP_BEST ] );
        }

        DEBUG_STRING_APPEND(sDebug, sDebugBestMode)
        pcCU ->setTransformSkipPartRange                ( bestTransformSkipMode,     compID, subTUAbsPartIdx, partIdxesPerSubTU );
        pcCU ->setCrossComponentPredictionAlphaPartRange( bestCrossCPredictionAlpha, compID, subTUAbsPartIdx, partIdxesPerSubTU );
        ruiDist += singleDistC;
      } while (TUIterator.nextSection(rTu));

      if (splitIntoSubTUs)
      {
        offsetSubTUCBFs(rTu, compID);
      }
    }
  }
  else
  {
    UInt    uiSplitCbf[MAX_NUM_COMPONENT] = {0,0,0};

    TComTURecurse tuRecurseChild(rTu, false);
    const UInt uiTrDepthChild   = tuRecurseChild.GetTransformDepthRel();
    do
    {
      DEBUG_STRING_NEW(sChild)

      xRecurIntraChromaCodingQT( pcOrgYuv, pcPredYuv, pcResiYuv, resiLuma, ruiDist, tuRecurseChild DEBUG_STRING_PASS_INTO(sChild) );

      DEBUG_STRING_APPEND(sDebug, sChild)
      const UInt uiAbsPartIdxSub=tuRecurseChild.GetAbsPartIdxTU();

      for(UInt ch=COMPONENT_Cb; ch<numberValidComponents; ch++)
      {
        uiSplitCbf[ch] |= pcCU->getCbf( uiAbsPartIdxSub, ComponentID(ch), uiTrDepthChild );
      }
    } while ( tuRecurseChild.nextSection(rTu) );


    UInt uiPartsDiv = rTu.GetAbsPartIdxNumParts();
    for(UInt ch=COMPONENT_Cb; ch<numberValidComponents; ch++)
    {
      if (uiSplitCbf[ch])
      {
        const UInt flag=1<<uiTrDepth;
        ComponentID compID=ComponentID(ch);
        UChar *pBase=pcCU->getCbf( compID );
        for( UInt uiOffs = 0; uiOffs < uiPartsDiv; uiOffs++ )
        {
          pBase[ uiAbsPartIdx + uiOffs ] |= flag;
        }
      }
    }
  }
}




Void
TEncSearch::xSetIntraResultChromaQT(TComYuv*    pcRecoYuv, TComTU &rTu)
{
  if (!rTu.ProcessChannelSection(CHANNEL_TYPE_CHROMA))
  {
    return;
  }
  TComDataCU *pcCU=rTu.getCU();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  const UInt uiTrDepth   = rTu.GetTransformDepthRel();
  UInt uiTrMode     = pcCU->getTransformIdx( uiAbsPartIdx );
  if(  uiTrMode == uiTrDepth )
  {
    UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
    UInt uiQTLayer    = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;

    //===== copy transform coefficients =====
    const TComRectangle &tuRectCb=rTu.getRect(COMPONENT_Cb);
    UInt uiNumCoeffC    = tuRectCb.width*tuRectCb.height;//( pcCU->getSlice()->getSPS()->getMaxCUWidth() * pcCU->getSlice()->getSPS()->getMaxCUHeight() ) >> ( uiFullDepth << 1 );
    const UInt offset = rTu.getCoefficientOffset(COMPONENT_Cb);

    const UInt numberValidComponents = getNumberValidComponents(rTu.GetChromaFormat());
    for (UInt ch=COMPONENT_Cb; ch<numberValidComponents; ch++)
    {
      const ComponentID component = ComponentID(ch);
      const TCoeff* src           = m_ppcQTTempCoeff[component][uiQTLayer] + offset;//(uiNumCoeffIncC*uiAbsPartIdx);
      TCoeff* dest                = pcCU->getCoeff(component) + offset;//(uiNumCoeffIncC*uiAbsPartIdx);
      ::memcpy( dest, src, sizeof(TCoeff)*uiNumCoeffC );
#if ADAPTIVE_QP_SELECTION
      TCoeff* pcArlCoeffSrc = m_ppcQTTempArlCoeff[component][ uiQTLayer ] + offset;//( uiNumCoeffIncC * uiAbsPartIdx );
      TCoeff* pcArlCoeffDst = pcCU->getArlCoeff(component)                + offset;//( uiNumCoeffIncC * uiAbsPartIdx );
      ::memcpy( pcArlCoeffDst, pcArlCoeffSrc, sizeof( TCoeff ) * uiNumCoeffC );
#endif
    }

    //===== copy reconstruction =====

    m_pcQTTempTComYuv[ uiQTLayer ].copyPartToPartComponent( COMPONENT_Cb, pcRecoYuv, uiAbsPartIdx, tuRectCb.width, tuRectCb.height );
    m_pcQTTempTComYuv[ uiQTLayer ].copyPartToPartComponent( COMPONENT_Cr, pcRecoYuv, uiAbsPartIdx, tuRectCb.width, tuRectCb.height );
  }
  else
  {
    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
      xSetIntraResultChromaQT( pcRecoYuv, tuRecurseChild );
    } while (tuRecurseChild.nextSection(rTu));
  }
}



Void
TEncSearch::estIntraPredLumaQT(TComDataCU* pcCU,
                               TComYuv*    pcOrgYuv,
                               TComYuv*    pcPredYuv,
                               TComYuv*    pcResiYuv,
                               TComYuv*    pcRecoYuv,
                               Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
                               DEBUG_STRING_FN_DECLARE(sDebug))
{
  const UInt         uiDepth               = pcCU->getDepth(0);
  const UInt         uiInitTrDepth         = pcCU->getPartitionSize(0) == SIZE_2Nx2N ? 0 : 1;
  const UInt         uiNumPU               = 1<<(2*uiInitTrDepth);
  const UInt         uiQNumParts           = pcCU->getTotalNumPart() >> 2;
  const UInt         uiWidthBit            = pcCU->getIntraSizeIdx(0);
  const ChromaFormat chFmt                 = pcCU->getPic()->getChromaFormat();
  const UInt         numberValidComponents = getNumberValidComponents(chFmt);
  const TComSPS     &sps                   = *(pcCU->getSlice()->getSPS());
  const TComPPS     &pps                   = *(pcCU->getSlice()->getPPS());
        Distortion   uiOverallDistY        = 0;
        UInt         CandNum;
        Double       CandCostList[ FAST_UDI_MAX_RDMODE_NUM ];
        Pel          resiLumaPU[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE];

        Bool    bMaintainResidual[NUMBER_OF_STORED_RESIDUAL_TYPES];
        for (UInt residualTypeIndex = 0; residualTypeIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; residualTypeIndex++)
        {
          bMaintainResidual[residualTypeIndex] = true; //assume true unless specified otherwise
        }

        bMaintainResidual[RESIDUAL_ENCODER_SIDE] = !(m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate());

  // Lambda calculation at equivalent Qp of 4 is recommended because at that Qp, the quantisation divisor is 1.
#if FULL_NBIT
  const Double sqrtLambdaForFirstPass= (m_pcEncCfg->getCostMode()==COST_MIXED_LOSSLESS_LOSSY_CODING && pcCU->getCUTransquantBypass(0)) ?
                sqrt(0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12) / 3.0)))
              : m_pcRdCost->getSqrtLambda();
#else
  const Double sqrtLambdaForFirstPass= (m_pcEncCfg->getCostMode()==COST_MIXED_LOSSLESS_LOSSY_CODING && pcCU->getCUTransquantBypass(0)) ?
                sqrt(0.57 * pow(2.0, ((LOSSLESS_AND_MIXED_LOSSLESS_RD_COST_TEST_QP_PRIME - 12 - 6 * (sps.getBitDepth(CHANNEL_TYPE_LUMA) - 8)) / 3.0)))
              : m_pcRdCost->getSqrtLambda();
#endif

  //===== set QP and clear Cbf =====
  if ( pps.getUseDQP() == true)
  {
    pcCU->setQPSubParts( pcCU->getQP(0), 0, uiDepth );
  }
  else
  {
    pcCU->setQPSubParts( pcCU->getSlice()->getSliceQp(), 0, uiDepth );
  }

  //===== loop over partitions =====
  TComTURecurse tuRecurseCU(pcCU, 0);
  TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth==0)?TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);

  do
  {
    const UInt uiPartOffset=tuRecurseWithPU.GetAbsPartIdxTU();
//  for( UInt uiPU = 0, uiPartOffset=0; uiPU < uiNumPU; uiPU++, uiPartOffset += uiQNumParts )
  //{
    //===== init pattern for luma prediction =====
    DEBUG_STRING_NEW(sTemp2)

    //===== determine set of modes to be tested (using prediction signal only) =====
    Int numModesAvailable     = 35; //total number of Intra modes
    UInt uiRdModeList[FAST_UDI_MAX_RDMODE_NUM];
    Int numModesForFullRD = m_pcEncCfg->getFastUDIUseMPMEnabled()?g_aucIntraModeNumFast_UseMPM[ uiWidthBit ] : g_aucIntraModeNumFast_NotUseMPM[ uiWidthBit ];

    // this should always be true
    assert (tuRecurseWithPU.ProcessComponentSection(COMPONENT_Y));
    initIntraPatternChType( tuRecurseWithPU, COMPONENT_Y, true DEBUG_STRING_PASS_INTO(sTemp2) );

    Bool doFastSearch = (numModesForFullRD != numModesAvailable);
    if (doFastSearch)
    {
      assert(numModesForFullRD < numModesAvailable);

      for( Int i=0; i < numModesForFullRD; i++ )
      {
        CandCostList[ i ] = MAX_DOUBLE;
      }
      CandNum = 0;

      const TComRectangle &puRect=tuRecurseWithPU.getRect(COMPONENT_Y);
      const UInt uiAbsPartIdx=tuRecurseWithPU.GetAbsPartIdxTU();

      Pel* piOrg         = pcOrgYuv ->getAddr( COMPONENT_Y, uiAbsPartIdx );
      Pel* piPred        = pcPredYuv->getAddr( COMPONENT_Y, uiAbsPartIdx );
      UInt uiStride      = pcPredYuv->getStride( COMPONENT_Y );
      DistParam distParam;
      const Bool bUseHadamard=pcCU->getCUTransquantBypass(0) == 0;
      m_pcRdCost->setDistParam(distParam, sps.getBitDepth(CHANNEL_TYPE_LUMA), piOrg, uiStride, piPred, uiStride, puRect.width, puRect.height, bUseHadamard);
      distParam.bApplyWeight = false;
      for( Int modeIdx = 0; modeIdx < numModesAvailable; modeIdx++ )
      {
        UInt       uiMode = modeIdx;
        Distortion uiSad  = 0;

        const Bool bUseFilter=TComPrediction::filteringIntraReferenceSamples(COMPONENT_Y, uiMode, puRect.width, puRect.height, chFmt, sps.getSpsRangeExtension().getIntraSmoothingDisabledFlag());

        predIntraAng( COMPONENT_Y, uiMode, piOrg, uiStride, piPred, uiStride, tuRecurseWithPU, bUseFilter, TComPrediction::UseDPCMForFirstPassIntraEstimation(tuRecurseWithPU, uiMode) );

        // use hadamard transform here
        uiSad+=distParam.DistFunc(&distParam);

        UInt   iModeBits = 0;

        // NB xModeBitsIntra will not affect the mode for chroma that may have already been pre-estimated.
        iModeBits+=xModeBitsIntra( pcCU, uiMode, uiPartOffset, uiDepth, CHANNEL_TYPE_LUMA );

        Double cost      = (Double)uiSad + (Double)iModeBits * sqrtLambdaForFirstPass;

#if DEBUG_INTRA_SEARCH_COSTS
        std::cout << "1st pass mode " << uiMode << " SAD = " << uiSad << ", mode bits = " << iModeBits << ", cost = " << cost << "\n";
#endif

        CandNum += xUpdateCandList( uiMode, cost, numModesForFullRD, uiRdModeList, CandCostList );
      }

      if (m_pcEncCfg->getFastUDIUseMPMEnabled())
      {
        Int uiPreds[NUM_MOST_PROBABLE_MODES] = {-1, -1, -1};

        Int iMode = -1;
        pcCU->getIntraDirPredictor( uiPartOffset, uiPreds, COMPONENT_Y, &iMode );

        const Int numCand = ( iMode >= 0 ) ? iMode : Int(NUM_MOST_PROBABLE_MODES);

        for( Int j=0; j < numCand; j++)
        {
          Bool mostProbableModeIncluded = false;
          Int mostProbableMode = uiPreds[j];

          for( Int i=0; i < numModesForFullRD; i++)
          {
            mostProbableModeIncluded |= (mostProbableMode == uiRdModeList[i]);
          }
          if (!mostProbableModeIncluded)
          {
            uiRdModeList[numModesForFullRD++] = mostProbableMode;
          }
        }
      }
    }
    else
    {
      for( Int i=0; i < numModesForFullRD; i++)
      {
        uiRdModeList[i] = i;
      }
    }

    //===== check modes (using r-d costs) =====
#if HHI_RQT_INTRA_SPEEDUP_MOD
    UInt   uiSecondBestMode  = MAX_UINT;
    Double dSecondBestPUCost = MAX_DOUBLE;
#endif
    DEBUG_STRING_NEW(sPU)
    UInt       uiBestPUMode  = 0;
    Distortion uiBestPUDistY = 0;
    Double     dBestPUCost   = MAX_DOUBLE;

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
    UInt max=numModesForFullRD;

    if (DebugOptionList::ForceLumaMode.isSet())
    {
      max=0;  // we are forcing a direction, so don't bother with mode check
    }
    for ( UInt uiMode = 0; uiMode < max; uiMode++)
#else
    for( UInt uiMode = 0; uiMode < numModesForFullRD; uiMode++ )
#endif
    {
      // set luma prediction mode
      UInt uiOrgMode = uiRdModeList[uiMode];

      pcCU->setIntraDirSubParts ( CHANNEL_TYPE_LUMA, uiOrgMode, uiPartOffset, uiDepth + uiInitTrDepth );

      DEBUG_STRING_NEW(sMode)
      // set context models
      m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST] );

      // determine residual for partition
      Distortion uiPUDistY = 0;
      Double     dPUCost   = 0.0;
#if HHI_RQT_INTRA_SPEEDUP
      xRecurIntraCodingLumaQT( pcOrgYuv, pcPredYuv, pcResiYuv, resiLumaPU, uiPUDistY, true, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode) );
#else
      xRecurIntraCodingLumaQT( pcOrgYuv, pcPredYuv, pcResiYuv, resiLumaPU, uiPUDistY, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode) );
#endif

#if DEBUG_INTRA_SEARCH_COSTS
      std::cout << "2nd pass [luma,chroma] mode [" << Int(pcCU->getIntraDir(CHANNEL_TYPE_LUMA, uiPartOffset)) << "," << Int(pcCU->getIntraDir(CHANNEL_TYPE_CHROMA, uiPartOffset)) << "] cost = " << dPUCost << "\n";
#endif

      // check r-d cost
      if( dPUCost < dBestPUCost )
      {
        DEBUG_STRING_SWAP(sPU, sMode)
#if HHI_RQT_INTRA_SPEEDUP_MOD
        uiSecondBestMode  = uiBestPUMode;
        dSecondBestPUCost = dBestPUCost;
#endif
        uiBestPUMode  = uiOrgMode;
        uiBestPUDistY = uiPUDistY;
        dBestPUCost   = dPUCost;

        xSetIntraResultLumaQT( pcRecoYuv, tuRecurseWithPU );

        if (pps.getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
        {
          const Int xOffset = tuRecurseWithPU.getRect( COMPONENT_Y ).x0;
          const Int yOffset = tuRecurseWithPU.getRect( COMPONENT_Y ).y0;
          for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
          {
            if (bMaintainResidual[storedResidualIndex])
            {
              xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex], resiLumaPU[storedResidualIndex], tuRecurseWithPU, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE );
            }
          }
        }

        UInt uiQPartNum = tuRecurseWithPU.GetAbsPartIdxNumParts();

        ::memcpy( m_puhQTTempTrIdx,  pcCU->getTransformIdx()       + uiPartOffset, uiQPartNum * sizeof( UChar ) );
        for (UInt component = 0; component < numberValidComponents; component++)
        {
          const ComponentID compID = ComponentID(component);
          ::memcpy( m_puhQTTempCbf[compID], pcCU->getCbf( compID  ) + uiPartOffset, uiQPartNum * sizeof( UChar ) );
          ::memcpy( m_puhQTTempTransformSkipFlag[compID],  pcCU->getTransformSkip(compID)  + uiPartOffset, uiQPartNum * sizeof( UChar ) );
        }
      }
#if HHI_RQT_INTRA_SPEEDUP_MOD
      else if( dPUCost < dSecondBestPUCost )
      {
        uiSecondBestMode  = uiOrgMode;
        dSecondBestPUCost = dPUCost;
      }
#endif
    } // Mode loop

#if HHI_RQT_INTRA_SPEEDUP
#if HHI_RQT_INTRA_SPEEDUP_MOD
    for( UInt ui =0; ui < 2; ++ui )
#endif
    {
#if HHI_RQT_INTRA_SPEEDUP_MOD
      UInt uiOrgMode   = ui ? uiSecondBestMode  : uiBestPUMode;
      if( uiOrgMode == MAX_UINT )
      {
        break;
      }
#else
      UInt uiOrgMode = uiBestPUMode;
#endif

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
      if (DebugOptionList::ForceLumaMode.isSet())
      {
        uiOrgMode = DebugOptionList::ForceLumaMode.getInt();
      }
#endif

      pcCU->setIntraDirSubParts ( CHANNEL_TYPE_LUMA, uiOrgMode, uiPartOffset, uiDepth + uiInitTrDepth );
      DEBUG_STRING_NEW(sModeTree)

      // set context models
      m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST] );

      // determine residual for partition
      Distortion uiPUDistY = 0;
      Double     dPUCost   = 0.0;

      xRecurIntraCodingLumaQT( pcOrgYuv, pcPredYuv, pcResiYuv, resiLumaPU, uiPUDistY, false, dPUCost, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sModeTree));

      // check r-d cost
      if( dPUCost < dBestPUCost )
      {
        DEBUG_STRING_SWAP(sPU, sModeTree)
        uiBestPUMode  = uiOrgMode;
        uiBestPUDistY = uiPUDistY;
        dBestPUCost   = dPUCost;

        xSetIntraResultLumaQT( pcRecoYuv, tuRecurseWithPU );

        if (pps.getPpsRangeExtension().getCrossComponentPredictionEnabledFlag())
        {
          const Int xOffset = tuRecurseWithPU.getRect( COMPONENT_Y ).x0;
          const Int yOffset = tuRecurseWithPU.getRect( COMPONENT_Y ).y0;
          for (UInt storedResidualIndex = 0; storedResidualIndex < NUMBER_OF_STORED_RESIDUAL_TYPES; storedResidualIndex++)
          {
            if (bMaintainResidual[storedResidualIndex])
            {
              xStoreCrossComponentPredictionResult(resiLuma[storedResidualIndex], resiLumaPU[storedResidualIndex], tuRecurseWithPU, xOffset, yOffset, MAX_CU_SIZE, MAX_CU_SIZE );
            }
          }
        }

        const UInt uiQPartNum = tuRecurseWithPU.GetAbsPartIdxNumParts();
        ::memcpy( m_puhQTTempTrIdx,  pcCU->getTransformIdx()       + uiPartOffset, uiQPartNum * sizeof( UChar ) );

        for (UInt component = 0; component < numberValidComponents; component++)
        {
          const ComponentID compID = ComponentID(component);
          ::memcpy( m_puhQTTempCbf[compID], pcCU->getCbf( compID  ) + uiPartOffset, uiQPartNum * sizeof( UChar ) );
          ::memcpy( m_puhQTTempTransformSkipFlag[compID],  pcCU->getTransformSkip(compID)  + uiPartOffset, uiQPartNum * sizeof( UChar ) );
        }
      }
    } // Mode loop
#endif

    DEBUG_STRING_APPEND(sDebug, sPU)

    //--- update overall distortion ---
    uiOverallDistY += uiBestPUDistY;

    //--- update transform index and cbf ---
    const UInt uiQPartNum = tuRecurseWithPU.GetAbsPartIdxNumParts();
    ::memcpy( pcCU->getTransformIdx()       + uiPartOffset, m_puhQTTempTrIdx,  uiQPartNum * sizeof( UChar ) );
    for (UInt component = 0; component < numberValidComponents; component++)
    {
      const ComponentID compID = ComponentID(component);
      ::memcpy( pcCU->getCbf( compID  ) + uiPartOffset, m_puhQTTempCbf[compID], uiQPartNum * sizeof( UChar ) );
      ::memcpy( pcCU->getTransformSkip( compID  ) + uiPartOffset, m_puhQTTempTransformSkipFlag[compID ], uiQPartNum * sizeof( UChar ) );
    }

    //--- set reconstruction for next intra prediction blocks ---
    if( !tuRecurseWithPU.IsLastSection() )
    {
      const TComRectangle &puRect=tuRecurseWithPU.getRect(COMPONENT_Y);
      const UInt  uiCompWidth   = puRect.width;
      const UInt  uiCompHeight  = puRect.height;

      const UInt  uiZOrder      = pcCU->getZorderIdxInCtu() + uiPartOffset;
            Pel*  piDes         = pcCU->getPic()->getPicYuvRec()->getAddr( COMPONENT_Y, pcCU->getCtuRsAddr(), uiZOrder );
      const UInt  uiDesStride   = pcCU->getPic()->getPicYuvRec()->getStride( COMPONENT_Y);
      const Pel*  piSrc         = pcRecoYuv->getAddr( COMPONENT_Y, uiPartOffset );
      const UInt  uiSrcStride   = pcRecoYuv->getStride( COMPONENT_Y);

      for( UInt uiY = 0; uiY < uiCompHeight; uiY++, piSrc += uiSrcStride, piDes += uiDesStride )
      {
        for( UInt uiX = 0; uiX < uiCompWidth; uiX++ )
        {
          piDes[ uiX ] = piSrc[ uiX ];
        }
      }
    }

    //=== update PU data ====
    pcCU->setIntraDirSubParts     ( CHANNEL_TYPE_LUMA, uiBestPUMode, uiPartOffset, uiDepth + uiInitTrDepth );
	
  } while (tuRecurseWithPU.nextSection(tuRecurseCU));


  if( uiNumPU > 1 )
  { // set Cbf for all blocks
    UInt uiCombCbfY = 0;
    UInt uiCombCbfU = 0;
    UInt uiCombCbfV = 0;
    UInt uiPartIdx  = 0;
    for( UInt uiPart = 0; uiPart < 4; uiPart++, uiPartIdx += uiQNumParts )
    {
      uiCombCbfY |= pcCU->getCbf( uiPartIdx, COMPONENT_Y,  1 );
      uiCombCbfU |= pcCU->getCbf( uiPartIdx, COMPONENT_Cb, 1 );
      uiCombCbfV |= pcCU->getCbf( uiPartIdx, COMPONENT_Cr, 1 );
    }
    for( UInt uiOffs = 0; uiOffs < 4 * uiQNumParts; uiOffs++ )
    {
      pcCU->getCbf( COMPONENT_Y  )[ uiOffs ] |= uiCombCbfY;
      pcCU->getCbf( COMPONENT_Cb )[ uiOffs ] |= uiCombCbfU;
      pcCU->getCbf( COMPONENT_Cr )[ uiOffs ] |= uiCombCbfV;
    }
  }

  //===== reset context models =====
  m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);

  //===== set distortion (rate and r-d costs are determined later) =====
  pcCU->getTotalDistortion() = uiOverallDistY;
}




Void
TEncSearch::estIntraPredChromaQT(TComDataCU* pcCU,
                                 TComYuv*    pcOrgYuv,
                                 TComYuv*    pcPredYuv,
                                 TComYuv*    pcResiYuv,
                                 TComYuv*    pcRecoYuv,
                                 Pel         resiLuma[NUMBER_OF_STORED_RESIDUAL_TYPES][MAX_CU_SIZE * MAX_CU_SIZE]
                                 DEBUG_STRING_FN_DECLARE(sDebug))
{
  const UInt    uiInitTrDepth  = pcCU->getPartitionSize(0) != SIZE_2Nx2N && enable4ChromaPUsInIntraNxNCU(pcOrgYuv->getChromaFormat()) ? 1 : 0;

  TComTURecurse tuRecurseCU(pcCU, 0);
  TComTURecurse tuRecurseWithPU(tuRecurseCU, false, (uiInitTrDepth==0)?TComTU::DONT_SPLIT : TComTU::QUAD_SPLIT);
  const UInt    uiQNumParts    = tuRecurseWithPU.GetAbsPartIdxNumParts();
  const UInt    uiDepthCU=tuRecurseWithPU.getCUDepth();
  const UInt    numberValidComponents = pcCU->getPic()->getNumberValidComponents();

  do
  {
    UInt       uiBestMode  = 0;
    Distortion uiBestDist  = 0;
    Double     dBestCost   = MAX_DOUBLE;

    //----- init mode list -----
    if (tuRecurseWithPU.ProcessChannelSection(CHANNEL_TYPE_CHROMA))
    {
      UInt uiModeList[FAST_UDI_MAX_RDMODE_NUM];
      const UInt  uiQPartNum     = uiQNumParts;
      const UInt  uiPartOffset   = tuRecurseWithPU.GetAbsPartIdxTU();
      {
        UInt  uiMinMode = 0;
        UInt  uiMaxMode = NUM_CHROMA_MODE;

        //----- check chroma modes -----
        pcCU->getAllowedChromaDir( uiPartOffset, uiModeList );

#if ENVIRONMENT_VARIABLE_DEBUG_AND_TEST
        if (DebugOptionList::ForceChromaMode.isSet())
        {
          uiMinMode=DebugOptionList::ForceChromaMode.getInt();
          if (uiModeList[uiMinMode]==34)
          {
            uiMinMode=4; // if the fixed mode has been renumbered because DM_CHROMA covers it, use DM_CHROMA.
          }
          uiMaxMode=uiMinMode+1;
        }
#endif

        DEBUG_STRING_NEW(sPU)

        for( UInt uiMode = uiMinMode; uiMode < uiMaxMode; uiMode++ )
        {
          //----- restore context models -----
          m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[uiDepthCU][CI_CURR_BEST] );
          
          DEBUG_STRING_NEW(sMode)
          //----- chroma coding -----
          Distortion uiDist = 0;
          pcCU->setIntraDirSubParts  ( CHANNEL_TYPE_CHROMA, uiModeList[uiMode], uiPartOffset, uiDepthCU+uiInitTrDepth );
          xRecurIntraChromaCodingQT       ( pcOrgYuv, pcPredYuv, pcResiYuv, resiLuma, uiDist, tuRecurseWithPU DEBUG_STRING_PASS_INTO(sMode) );

          if( pcCU->getSlice()->getPPS()->getUseTransformSkip() )
          {
            m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[uiDepthCU][CI_CURR_BEST] );
          }

          UInt    uiBits = xGetIntraBitsQT( tuRecurseWithPU, false, true, false );
          Double  dCost  = m_pcRdCost->calcRdCost( uiBits, uiDist );

          //----- compare -----
          if( dCost < dBestCost )
          {
            DEBUG_STRING_SWAP(sPU, sMode);
            dBestCost   = dCost;
            uiBestDist  = uiDist;
            uiBestMode  = uiModeList[uiMode];

            xSetIntraResultChromaQT( pcRecoYuv, tuRecurseWithPU );
            for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberValidComponents; componentIndex++)
            {
              const ComponentID compID = ComponentID(componentIndex);
              ::memcpy( m_puhQTTempCbf[compID], pcCU->getCbf( compID )+uiPartOffset, uiQPartNum * sizeof( UChar ) );
              ::memcpy( m_puhQTTempTransformSkipFlag[compID], pcCU->getTransformSkip( compID )+uiPartOffset, uiQPartNum * sizeof( UChar ) );
              ::memcpy( m_phQTTempCrossComponentPredictionAlpha[compID], pcCU->getCrossComponentPredictionAlpha(compID)+uiPartOffset, uiQPartNum * sizeof( SChar ) );
            }
          }
        }

        DEBUG_STRING_APPEND(sDebug, sPU)

        //----- set data -----
        for (UInt componentIndex = COMPONENT_Cb; componentIndex < numberValidComponents; componentIndex++)
        {
          const ComponentID compID = ComponentID(componentIndex);
          ::memcpy( pcCU->getCbf( compID )+uiPartOffset, m_puhQTTempCbf[compID], uiQPartNum * sizeof( UChar ) );
          ::memcpy( pcCU->getTransformSkip( compID )+uiPartOffset, m_puhQTTempTransformSkipFlag[compID], uiQPartNum * sizeof( UChar ) );
          ::memcpy( pcCU->getCrossComponentPredictionAlpha(compID)+uiPartOffset, m_phQTTempCrossComponentPredictionAlpha[compID], uiQPartNum * sizeof( SChar ) );
        }
      }

      if( ! tuRecurseWithPU.IsLastSection() )
      {
        for (UInt ch=COMPONENT_Cb; ch<numberValidComponents; ch++)
        {
          const ComponentID compID    = ComponentID(ch);
          const TComRectangle &tuRect = tuRecurseWithPU.getRect(compID);
          const UInt  uiCompWidth     = tuRect.width;
          const UInt  uiCompHeight    = tuRect.height;
          const UInt  uiZOrder        = pcCU->getZorderIdxInCtu() + tuRecurseWithPU.GetAbsPartIdxTU();
                Pel*  piDes           = pcCU->getPic()->getPicYuvRec()->getAddr( compID, pcCU->getCtuRsAddr(), uiZOrder );
          const UInt  uiDesStride     = pcCU->getPic()->getPicYuvRec()->getStride( compID);
          const Pel*  piSrc           = pcRecoYuv->getAddr( compID, uiPartOffset );
          const UInt  uiSrcStride     = pcRecoYuv->getStride( compID);

          for( UInt uiY = 0; uiY < uiCompHeight; uiY++, piSrc += uiSrcStride, piDes += uiDesStride )
          {
            for( UInt uiX = 0; uiX < uiCompWidth; uiX++ )
            {
              piDes[ uiX ] = piSrc[ uiX ];
            }
          }
        }
      }

      pcCU->setIntraDirSubParts( CHANNEL_TYPE_CHROMA, uiBestMode, uiPartOffset, uiDepthCU+uiInitTrDepth );
      pcCU->getTotalDistortion      () += uiBestDist;
    }

  } while (tuRecurseWithPU.nextSection(tuRecurseCU));

  //----- restore context models -----

  if( uiInitTrDepth != 0 )
  { // set Cbf for all blocks
    UInt uiCombCbfU = 0;
    UInt uiCombCbfV = 0;
    UInt uiPartIdx  = 0;
    for( UInt uiPart = 0; uiPart < 4; uiPart++, uiPartIdx += uiQNumParts )
    {
      uiCombCbfU |= pcCU->getCbf( uiPartIdx, COMPONENT_Cb, 1 );
      uiCombCbfV |= pcCU->getCbf( uiPartIdx, COMPONENT_Cr, 1 );
    }
    for( UInt uiOffs = 0; uiOffs < 4 * uiQNumParts; uiOffs++ )
    {
      pcCU->getCbf( COMPONENT_Cb )[ uiOffs ] |= uiCombCbfU;
      pcCU->getCbf( COMPONENT_Cr )[ uiOffs ] |= uiCombCbfV;
    }
  }

  m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[uiDepthCU][CI_CURR_BEST] );
}




/** Function for encoding and reconstructing luma/chroma samples of a PCM mode CU.
 * \param pcCU pointer to current CU
 * \param uiAbsPartIdx part index
 * \param pOrg pointer to original sample arrays
 * \param pPCM pointer to PCM code arrays
 * \param pPred pointer to prediction signal arrays
 * \param pResi pointer to residual signal arrays
 * \param pReco pointer to reconstructed sample arrays
 * \param uiStride stride of the original/prediction/residual sample arrays
 * \param uiWidth block width
 * \param uiHeight block height
 * \param compID texture component type
 */
Void TEncSearch::xEncPCM (TComDataCU* pcCU, UInt uiAbsPartIdx, Pel* pOrg, Pel* pPCM, Pel* pPred, Pel* pResi, Pel* pReco, UInt uiStride, UInt uiWidth, UInt uiHeight, const ComponentID compID )
{
  const UInt uiReconStride   = pcCU->getPic()->getPicYuvRec()->getStride(compID);
  const UInt uiPCMBitDepth   = pcCU->getSlice()->getSPS()->getPCMBitDepth(toChannelType(compID));
  const Int  channelBitDepth = pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
  Pel* pRecoPic = pcCU->getPic()->getPicYuvRec()->getAddr(compID, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu()+uiAbsPartIdx);

  const Int pcmShiftRight=(channelBitDepth - Int(uiPCMBitDepth));

  assert(pcmShiftRight >= 0);

  for( UInt uiY = 0; uiY < uiHeight; uiY++ )
  {
    for( UInt uiX = 0; uiX < uiWidth; uiX++ )
    {
      // Reset pred and residual
      pPred[uiX] = 0;
      pResi[uiX] = 0;
      // Encode
      pPCM[uiX] = (pOrg[uiX]>>pcmShiftRight);
      // Reconstruction
      pReco   [uiX] = (pPCM[uiX]<<(pcmShiftRight));
      pRecoPic[uiX] = pReco[uiX];
    }
    pPred += uiStride;
    pResi += uiStride;
    pPCM += uiWidth;
    pOrg += uiStride;
    pReco += uiStride;
    pRecoPic += uiReconStride;
  }
}


//!  Function for PCM mode estimation.
Void TEncSearch::IPCMSearch( TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv )
{
  UInt              uiDepth      = pcCU->getDepth(0);
  const Distortion  uiDistortion = 0;
  UInt              uiBits;

  Double dCost;

  for (UInt ch=0; ch < pcCU->getPic()->getNumberValidComponents(); ch++)
  {
    const ComponentID compID  = ComponentID(ch);
    const UInt width  = pcCU->getWidth(0)  >> pcCU->getPic()->getComponentScaleX(compID);
    const UInt height = pcCU->getHeight(0) >> pcCU->getPic()->getComponentScaleY(compID);
    const UInt stride = pcPredYuv->getStride(compID);

    Pel * pOrig    = pcOrgYuv->getAddr  (compID, 0, width);
    Pel * pResi    = pcResiYuv->getAddr(compID, 0, width);
    Pel * pPred    = pcPredYuv->getAddr(compID, 0, width);
    Pel * pReco    = pcRecoYuv->getAddr(compID, 0, width);
    Pel * pPCM     = pcCU->getPCMSample (compID);

    xEncPCM ( pcCU, 0, pOrig, pPCM, pPred, pResi, pReco, stride, width, height, compID );

  }

  m_pcEntropyCoder->resetBits();
  xEncIntraHeader ( pcCU, uiDepth, 0, true, false);
  uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();

  dCost = m_pcRdCost->calcRdCost( uiBits, uiDistortion );

  m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST]);

  pcCU->getTotalBits()       = uiBits;
  pcCU->getTotalCost()       = dCost;
  pcCU->getTotalDistortion() = uiDistortion;

  pcCU->copyToPic(uiDepth);
}




Void TEncSearch::xGetInterPredictionError( TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPartIdx, Distortion& ruiErr, Bool /*bHadamard*/ )
{
  motionCompensation( pcCU, &m_tmpYuvPred, REF_PIC_LIST_X, iPartIdx );

  UInt uiAbsPartIdx = 0;
  Int iWidth = 0;
  Int iHeight = 0;
  pcCU->getPartIndexAndSize( iPartIdx, uiAbsPartIdx, iWidth, iHeight );

  DistParam cDistParam;

  cDistParam.bApplyWeight = false;


  m_pcRdCost->setDistParam( cDistParam, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA),
                            pcYuvOrg->getAddr( COMPONENT_Y, uiAbsPartIdx ), pcYuvOrg->getStride(COMPONENT_Y),
                            m_tmpYuvPred .getAddr( COMPONENT_Y, uiAbsPartIdx ), m_tmpYuvPred.getStride(COMPONENT_Y),
                            iWidth, iHeight, m_pcEncCfg->getUseHADME() && (pcCU->getCUTransquantBypass(iPartIdx) == 0) );

  ruiErr = cDistParam.DistFunc( &cDistParam );
}

//! estimation of best merge coding
Void TEncSearch::xMergeEstimation( TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPUIdx, UInt& uiInterDir, TComMvField* pacMvField, UInt& uiMergeIndex, Distortion& ruiCost, TComMvField* cMvFieldNeighbours, UChar* uhInterDirNeighbours, Int& numValidMergeCand )
{
	//ofstream myfile;

	//myfile.open("C:\\RACEHOECES2_CTU.csv", ios::app);

  UInt uiAbsPartIdx = 0;
  Int iWidth = 0;
  Int iHeight = 0;

 
  pcCU->getPartIndexAndSize( iPUIdx, uiAbsPartIdx, iWidth, iHeight );
  UInt uiDepth = pcCU->getDepth( uiAbsPartIdx );
 // CTUH1 = iHeight;
 // CTUW1 = iWidth;
 // myfile << CTUH1 << ',' << CTUW1 << endl;
  PartSize partSize = pcCU->getPartitionSize( 0 );
  if ( pcCU->getSlice()->getPPS()->getLog2ParallelMergeLevelMinus2() && partSize != SIZE_2Nx2N && pcCU->getWidth( 0 ) <= 8 )
  {
    if ( iPUIdx == 0 )
    {
      pcCU->setPartSizeSubParts( SIZE_2Nx2N, 0, uiDepth ); // temporarily set
      pcCU->getInterMergeCandidates( 0, 0, cMvFieldNeighbours,uhInterDirNeighbours, numValidMergeCand );
      pcCU->setPartSizeSubParts( partSize, 0, uiDepth ); // restore
    }
  }
  else
  {
    pcCU->getInterMergeCandidates( uiAbsPartIdx, iPUIdx, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand );
  }

  xRestrictBipredMergeCand( pcCU, iPUIdx, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand );

  ruiCost = std::numeric_limits<Distortion>::max();
  for( UInt uiMergeCand = 0; uiMergeCand < numValidMergeCand; ++uiMergeCand )
  {
    Distortion uiCostCand = std::numeric_limits<Distortion>::max();
    UInt       uiBitsCand = 0;

    PartSize ePartSize = pcCU->getPartitionSize( 0 );

    pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField( cMvFieldNeighbours[0 + 2*uiMergeCand], ePartSize, uiAbsPartIdx, 0, iPUIdx );
    pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField( cMvFieldNeighbours[1 + 2*uiMergeCand], ePartSize, uiAbsPartIdx, 0, iPUIdx );

    xGetInterPredictionError( pcCU, pcYuvOrg, iPUIdx, uiCostCand, m_pcEncCfg->getUseHADME() );
    uiBitsCand = uiMergeCand + 1;
    if (uiMergeCand == m_pcEncCfg->getMaxNumMergeCand() -1)
    {
        uiBitsCand--;
    }
    uiCostCand = uiCostCand + m_pcRdCost->getCost( uiBitsCand );
    if ( uiCostCand < ruiCost )
    {
      ruiCost = uiCostCand;
      pacMvField[0] = cMvFieldNeighbours[0 + 2*uiMergeCand];
      pacMvField[1] = cMvFieldNeighbours[1 + 2*uiMergeCand];
      uiInterDir = uhInterDirNeighbours[uiMergeCand];
      uiMergeIndex = uiMergeCand;
    }
  }
 
}

/** convert bi-pred merge candidates to uni-pred
 * \param pcCU
 * \param puIdx
 * \param mvFieldNeighbours
 * \param interDirNeighbours
 * \param numValidMergeCand
 * \returns Void
 */
Void TEncSearch::xRestrictBipredMergeCand( TComDataCU* pcCU, UInt puIdx, TComMvField* mvFieldNeighbours, UChar* interDirNeighbours, Int numValidMergeCand )
{
	
  if ( pcCU->isBipredRestriction(puIdx) )
  {
    for( UInt mergeCand = 0; mergeCand < numValidMergeCand; ++mergeCand )
    {
      if ( interDirNeighbours[mergeCand] == 3 )
      {
        interDirNeighbours[mergeCand] = 1;
        mvFieldNeighbours[(mergeCand << 1) + 1].setMvField(TComMv(0,0), -1);
      }
    }
  }
}

//! search of the best candidate for inter prediction
#if AMP_MRG
Void TEncSearch::predInterSearch( TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv DEBUG_STRING_FN_DECLARE(sDebug), Bool bUseRes, Bool bUseMRG )
#else
Void TEncSearch::predInterSearch( TComDataCU* pcCU, TComYuv* pcOrgYuv, TComYuv* pcPredYuv, TComYuv* pcResiYuv, TComYuv* pcRecoYuv, Bool bUseRes )
#endif
{
  for(UInt i=0; i<NUM_REF_PIC_LIST_01; i++)
  {
    m_acYuvPred[i].clear();
  }
  m_cYuvPredTemp.clear();
  pcPredYuv->clear();

  if ( !bUseRes )
  {
    pcResiYuv->clear();
  }

  pcRecoYuv->clear();
  
 // ofstream myfile;

//  myfile.open("C:\\RACEHOECES1_CTU.csv", ios::app);
  TComMv       cMvSrchRngLT;
  TComMv       cMvSrchRngRB;

  TComMv       cMvZero;
  TComMv       TempMv; //kolya

  TComMv       cMv[2];
  TComMv       cMvBi[2];
  TComMv       cMvTemp[2][33];

  Int          iNumPart    = pcCU->getNumPartitions();
  Int          iNumPredDir = pcCU->getSlice()->isInterP() ? 1 : 2;

  TComMv       cMvPred[2][33];

  TComMv       cMvPredBi[2][33];
  Int          aaiMvpIdxBi[2][33];

  Int          aaiMvpIdx[2][33];
  Int          aaiMvpNum[2][33];

  AMVPInfo     aacAMVPInfo[2][33];

  Int          iRefIdx[2]={0,0}; //If un-initialized, may cause SEGV in bi-directional prediction iterative stage.
  Int          iRefIdxBi[2];

  UInt         uiPartAddr;
  Int          iRoiWidth, iRoiHeight;

  UInt         uiMbBits[3] = {1, 1, 0};

  UInt         uiLastMode = 0;
  Int          iRefStart, iRefEnd;

  PartSize     ePartSize = pcCU->getPartitionSize( 0 );

  Int          bestBiPRefIdxL1 = 0;
  Int          bestBiPMvpL1 = 0;
  Distortion   biPDistTemp = std::numeric_limits<Distortion>::max();
  counter_ME = counter_ME + 1;

  TComMvField cMvFieldNeighbours[MRG_MAX_NUM_CANDS << 1]; // double length for mv of both lists
  UChar uhInterDirNeighbours[MRG_MAX_NUM_CANDS];
  Int numValidMergeCand = 0 ;

  for ( Int iPartIdx = 0; iPartIdx < iNumPart; iPartIdx++ )
  {
    Distortion   uiCost[2] = { std::numeric_limits<Distortion>::max(), std::numeric_limits<Distortion>::max() };
    Distortion   uiCostBi  =   std::numeric_limits<Distortion>::max();
    Distortion   uiCostTemp;

    UInt         uiBits[3];
    UInt         uiBitsTemp;
    Distortion   bestBiPDist = std::numeric_limits<Distortion>::max();

    Distortion   uiCostTempL0[MAX_NUM_REF];
    for (Int iNumRef=0; iNumRef < MAX_NUM_REF; iNumRef++)
    {
      uiCostTempL0[iNumRef] = std::numeric_limits<Distortion>::max();
    }
    UInt         uiBitsTempL0[MAX_NUM_REF];

    TComMv       mvValidList1;
    Int          refIdxValidList1 = 0;
    UInt         bitsValidList1 = MAX_UINT;
    Distortion   costValidList1 = std::numeric_limits<Distortion>::max();

    xGetBlkBits( ePartSize, pcCU->getSlice()->isInterP(), iPartIdx, uiLastMode, uiMbBits);

    pcCU->getPartIndexAndSize( iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight );
	
	
#if AMP_MRG
    Bool bTestNormalMC = true;

    if ( bUseMRG && pcCU->getWidth( 0 ) > 8 && iNumPart == 2 )
    {
      bTestNormalMC = false;
    }

    if (bTestNormalMC)
    {
#endif

    //  Uni-directional prediction
    for ( Int iRefList = 0; iRefList < iNumPredDir; iRefList++ )
    {
      RefPicList  eRefPicList = ( iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0 );

      for ( Int iRefIdxTemp = 0; iRefIdxTemp < pcCU->getSlice()->getNumRefIdx(eRefPicList); iRefIdxTemp++ )
      {
        uiBitsTemp = uiMbBits[iRefList];
        if ( pcCU->getSlice()->getNumRefIdx(eRefPicList) > 1 )
        {
          uiBitsTemp += iRefIdxTemp+1;
          if ( iRefIdxTemp == pcCU->getSlice()->getNumRefIdx(eRefPicList)-1 )
          {
            uiBitsTemp--;
          }
        }
        xEstimateMvPredAMVP( pcCU, pcOrgYuv, iPartIdx, eRefPicList, iRefIdxTemp, cMvPred[iRefList][iRefIdxTemp], false, &biPDistTemp);
        aaiMvpIdx[iRefList][iRefIdxTemp] = pcCU->getMVPIdx(eRefPicList, uiPartAddr);
        aaiMvpNum[iRefList][iRefIdxTemp] = pcCU->getMVPNum(eRefPicList, uiPartAddr);

        if(pcCU->getSlice()->getMvdL1ZeroFlag() && iRefList==1 && biPDistTemp < bestBiPDist)
        {
          bestBiPDist = biPDistTemp;
          bestBiPMvpL1 = aaiMvpIdx[iRefList][iRefIdxTemp];
          bestBiPRefIdxL1 = iRefIdxTemp;
        }

        uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdx[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];

        if ( m_pcEncCfg->getFastMEForGenBLowDelayEnabled() && iRefList == 1 )    // list 1
        {
          if ( pcCU->getSlice()->getList1IdxToList0Idx( iRefIdxTemp ) >= 0 )
          {
            cMvTemp[1][iRefIdxTemp] = cMvTemp[0][pcCU->getSlice()->getList1IdxToList0Idx( iRefIdxTemp )];
            uiCostTemp = uiCostTempL0[pcCU->getSlice()->getList1IdxToList0Idx( iRefIdxTemp )];
            /*first subtract the bit-rate part of the cost of the other list*/
            uiCostTemp -= m_pcRdCost->getCost( uiBitsTempL0[pcCU->getSlice()->getList1IdxToList0Idx( iRefIdxTemp )] );
            /*correct the bit-rate part of the current ref*/
            m_pcRdCost->setPredictor  ( cMvPred[iRefList][iRefIdxTemp] );
            uiBitsTemp += m_pcRdCost->getBitsOfVectorWithPredictor( cMvTemp[1][iRefIdxTemp].getHor(), cMvTemp[1][iRefIdxTemp].getVer() );
            /*calculate the correct cost*/
            uiCostTemp += m_pcRdCost->getCost( uiBitsTemp );
          }
          else
          {
            xMotionEstimation ( pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp );
          }
        }
        else
        {
          xMotionEstimation ( pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPred[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp );
        }
        xCopyAMVPInfo(pcCU->getCUMvField(eRefPicList)->getAMVPInfo(), &aacAMVPInfo[iRefList][iRefIdxTemp]); // must always be done ( also when AMVP_MODE = AM_NONE )
        xCheckBestMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPred[iRefList][iRefIdxTemp], aaiMvpIdx[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp);

        if ( iRefList == 0 )
        {
          uiCostTempL0[iRefIdxTemp] = uiCostTemp;
          uiBitsTempL0[iRefIdxTemp] = uiBitsTemp;
        }
        if ( uiCostTemp < uiCost[iRefList] )
        {
          uiCost[iRefList] = uiCostTemp;
          uiBits[iRefList] = uiBitsTemp; // storing for bi-prediction

          // set motion
          cMv[iRefList]     = cMvTemp[iRefList][iRefIdxTemp];
          iRefIdx[iRefList] = iRefIdxTemp;
        }

        if ( iRefList == 1 && uiCostTemp < costValidList1 && pcCU->getSlice()->getList1IdxToList0Idx( iRefIdxTemp ) < 0 )
        {
          costValidList1 = uiCostTemp;
          bitsValidList1 = uiBitsTemp;

          // set motion
          mvValidList1     = cMvTemp[iRefList][iRefIdxTemp];
          refIdxValidList1 = iRefIdxTemp;
        }
      }
    }

    //  Bi-predictive Motion estimation
    if ( (pcCU->getSlice()->isInterB()) && (pcCU->isBipredRestriction(iPartIdx) == false) )
    {

      cMvBi[0] = cMv[0];            cMvBi[1] = cMv[1];
      iRefIdxBi[0] = iRefIdx[0];    iRefIdxBi[1] = iRefIdx[1];

      ::memcpy(cMvPredBi, cMvPred, sizeof(cMvPred));
      ::memcpy(aaiMvpIdxBi, aaiMvpIdx, sizeof(aaiMvpIdx));

      UInt uiMotBits[2];

      if(pcCU->getSlice()->getMvdL1ZeroFlag())
      {
        xCopyAMVPInfo(&aacAMVPInfo[1][bestBiPRefIdxL1], pcCU->getCUMvField(REF_PIC_LIST_1)->getAMVPInfo());
        pcCU->setMVPIdxSubParts( bestBiPMvpL1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
        aaiMvpIdxBi[1][bestBiPRefIdxL1] = bestBiPMvpL1;
        cMvPredBi[1][bestBiPRefIdxL1]   = pcCU->getCUMvField(REF_PIC_LIST_1)->getAMVPInfo()->m_acMvCand[bestBiPMvpL1];

        cMvBi[1] = cMvPredBi[1][bestBiPRefIdxL1];
        iRefIdxBi[1] = bestBiPRefIdxL1;
        pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllMv( cMvBi[1], ePartSize, uiPartAddr, 0, iPartIdx );
        pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllRefIdx( iRefIdxBi[1], ePartSize, uiPartAddr, 0, iPartIdx );
        TComYuv* pcYuvPred = &m_acYuvPred[REF_PIC_LIST_1];
        motionCompensation( pcCU, pcYuvPred, REF_PIC_LIST_1, iPartIdx );

        uiMotBits[0] = uiBits[0] - uiMbBits[0];
        uiMotBits[1] = uiMbBits[1];

        if ( pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1) > 1 )
        {
          uiMotBits[1] += bestBiPRefIdxL1+1;
          if ( bestBiPRefIdxL1 == pcCU->getSlice()->getNumRefIdx(REF_PIC_LIST_1)-1 )
          {
            uiMotBits[1]--;
          }
        }

        uiMotBits[1] += m_auiMVPIdxCost[aaiMvpIdxBi[1][bestBiPRefIdxL1]][AMVP_MAX_NUM_CANDS];

        uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];

        cMvTemp[1][bestBiPRefIdxL1] = cMvBi[1];
      }
      else
      {
        uiMotBits[0] = uiBits[0] - uiMbBits[0];
        uiMotBits[1] = uiBits[1] - uiMbBits[1];
        uiBits[2] = uiMbBits[2] + uiMotBits[0] + uiMotBits[1];
      }

      // 4-times iteration (default)
      Int iNumIter = 4;

      // fast encoder setting: only one iteration
      if ( m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE2 || pcCU->getSlice()->getMvdL1ZeroFlag() )
      {
        iNumIter = 1;
      }

      for ( Int iIter = 0; iIter < iNumIter; iIter++ )
      {
        Int         iRefList    = iIter % 2;

        if ( m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode()==FASTINTERSEARCH_MODE2 )
        {
          if( uiCost[0] <= uiCost[1] )
          {
            iRefList = 1;
          }
          else
          {
            iRefList = 0;
          }
        }
        else if ( iIter == 0 )
        {
          iRefList = 0;
        }
        if ( iIter == 0 && !pcCU->getSlice()->getMvdL1ZeroFlag())
        {
          pcCU->getCUMvField(RefPicList(1-iRefList))->setAllMv( cMv[1-iRefList], ePartSize, uiPartAddr, 0, iPartIdx );
          pcCU->getCUMvField(RefPicList(1-iRefList))->setAllRefIdx( iRefIdx[1-iRefList], ePartSize, uiPartAddr, 0, iPartIdx );
          TComYuv*  pcYuvPred = &m_acYuvPred[1-iRefList];
          motionCompensation ( pcCU, pcYuvPred, RefPicList(1-iRefList), iPartIdx );
        }

        RefPicList  eRefPicList = ( iRefList ? REF_PIC_LIST_1 : REF_PIC_LIST_0 );

        if(pcCU->getSlice()->getMvdL1ZeroFlag())
        {
          iRefList = 0;
          eRefPicList = REF_PIC_LIST_0;
        }

        Bool bChanged = false;

        iRefStart = 0;
        iRefEnd   = pcCU->getSlice()->getNumRefIdx(eRefPicList)-1;

        for ( Int iRefIdxTemp = iRefStart; iRefIdxTemp <= iRefEnd; iRefIdxTemp++ )
        {
          uiBitsTemp = uiMbBits[2] + uiMotBits[1-iRefList];
          if ( pcCU->getSlice()->getNumRefIdx(eRefPicList) > 1 )
          {
            uiBitsTemp += iRefIdxTemp+1;
            if ( iRefIdxTemp == pcCU->getSlice()->getNumRefIdx(eRefPicList)-1 )
            {
              uiBitsTemp--;
            }
          }
          uiBitsTemp += m_auiMVPIdxCost[aaiMvpIdxBi[iRefList][iRefIdxTemp]][AMVP_MAX_NUM_CANDS];
          // call ME
          xMotionEstimation ( pcCU, pcOrgYuv, iPartIdx, eRefPicList, &cMvPredBi[iRefList][iRefIdxTemp], iRefIdxTemp, cMvTemp[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp, true );

          xCopyAMVPInfo(&aacAMVPInfo[iRefList][iRefIdxTemp], pcCU->getCUMvField(eRefPicList)->getAMVPInfo());
          xCheckBestMVP(pcCU, eRefPicList, cMvTemp[iRefList][iRefIdxTemp], cMvPredBi[iRefList][iRefIdxTemp], aaiMvpIdxBi[iRefList][iRefIdxTemp], uiBitsTemp, uiCostTemp);

          if ( uiCostTemp < uiCostBi )
          {
            bChanged = true;

            cMvBi[iRefList]     = cMvTemp[iRefList][iRefIdxTemp];
            iRefIdxBi[iRefList] = iRefIdxTemp;

            uiCostBi            = uiCostTemp;
            uiMotBits[iRefList] = uiBitsTemp - uiMbBits[2] - uiMotBits[1-iRefList];
            uiBits[2]           = uiBitsTemp;

            if(iNumIter!=1)
            {
              //  Set motion
              pcCU->getCUMvField( eRefPicList )->setAllMv( cMvBi[iRefList], ePartSize, uiPartAddr, 0, iPartIdx );
              pcCU->getCUMvField( eRefPicList )->setAllRefIdx( iRefIdxBi[iRefList], ePartSize, uiPartAddr, 0, iPartIdx );

              TComYuv* pcYuvPred = &m_acYuvPred[iRefList];
              motionCompensation( pcCU, pcYuvPred, eRefPicList, iPartIdx );
            }
          }
        } // for loop-iRefIdxTemp

        if ( !bChanged )
        {
          if ( uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1] )
          {
            xCopyAMVPInfo(&aacAMVPInfo[0][iRefIdxBi[0]], pcCU->getCUMvField(REF_PIC_LIST_0)->getAMVPInfo());
            xCheckBestMVP(pcCU, REF_PIC_LIST_0, cMvBi[0], cMvPredBi[0][iRefIdxBi[0]], aaiMvpIdxBi[0][iRefIdxBi[0]], uiBits[2], uiCostBi);
            if(!pcCU->getSlice()->getMvdL1ZeroFlag())
            {
              xCopyAMVPInfo(&aacAMVPInfo[1][iRefIdxBi[1]], pcCU->getCUMvField(REF_PIC_LIST_1)->getAMVPInfo());
              xCheckBestMVP(pcCU, REF_PIC_LIST_1, cMvBi[1], cMvPredBi[1][iRefIdxBi[1]], aaiMvpIdxBi[1][iRefIdxBi[1]], uiBits[2], uiCostBi);
            }
          }
          break;
        }
      } // for loop-iter
    } // if (B_SLICE)

#if AMP_MRG
    } //end if bTestNormalMC
#endif
    //  Clear Motion Field
    pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvField( TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx );
    pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvField( TComMvField(), ePartSize, uiPartAddr, 0, iPartIdx );
    pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd    ( cMvZero,       ePartSize, uiPartAddr, 0, iPartIdx );
    pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd    ( cMvZero,       ePartSize, uiPartAddr, 0, iPartIdx );

    pcCU->setMVPIdxSubParts( -1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
    pcCU->setMVPNumSubParts( -1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
    pcCU->setMVPIdxSubParts( -1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
    pcCU->setMVPNumSubParts( -1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

    UInt uiMEBits = 0;
    // Set Motion Field_
    cMv[1] = mvValidList1;
	
    iRefIdx[1] = refIdxValidList1;
    uiBits[1] = bitsValidList1;
    uiCost[1] = costValidList1;

#if AMP_MRG
    if (bTestNormalMC)
    {
#endif
    if ( uiCostBi <= uiCost[0] && uiCostBi <= uiCost[1])
    {
      uiLastMode = 2;
      pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMv( cMvBi[0], ePartSize, uiPartAddr, 0, iPartIdx );
      pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx( iRefIdxBi[0], ePartSize, uiPartAddr, 0, iPartIdx );
      pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMv( cMvBi[1], ePartSize, uiPartAddr, 0, iPartIdx );
      pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx( iRefIdxBi[1], ePartSize, uiPartAddr, 0, iPartIdx );

      TempMv = cMvBi[0] - cMvPredBi[0][iRefIdxBi[0]];
      pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd    ( TempMv,                 ePartSize, uiPartAddr, 0, iPartIdx );

      TempMv = cMvBi[1] - cMvPredBi[1][iRefIdxBi[1]];
      pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd    ( TempMv,                 ePartSize, uiPartAddr, 0, iPartIdx );

      pcCU->setInterDirSubParts( 3, uiPartAddr, iPartIdx, pcCU->getDepth(0) );

      pcCU->setMVPIdxSubParts( aaiMvpIdxBi[0][iRefIdxBi[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
      pcCU->setMVPNumSubParts( aaiMvpNum[0][iRefIdxBi[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
      pcCU->setMVPIdxSubParts( aaiMvpIdxBi[1][iRefIdxBi[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
      pcCU->setMVPNumSubParts( aaiMvpNum[1][iRefIdxBi[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

      uiMEBits = uiBits[2];
    }
    else if ( uiCost[0] <= uiCost[1] )
    {
      uiLastMode = 0;
      pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMv( cMv[0], ePartSize, uiPartAddr, 0, iPartIdx );
      pcCU->getCUMvField(REF_PIC_LIST_0)->setAllRefIdx( iRefIdx[0], ePartSize, uiPartAddr, 0, iPartIdx );

      TempMv = cMv[0] - cMvPred[0][iRefIdx[0]];
      pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd    ( TempMv,                 ePartSize, uiPartAddr, 0, iPartIdx );

      pcCU->setInterDirSubParts( 1, uiPartAddr, iPartIdx, pcCU->getDepth(0) );

      pcCU->setMVPIdxSubParts( aaiMvpIdx[0][iRefIdx[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
      pcCU->setMVPNumSubParts( aaiMvpNum[0][iRefIdx[0]], REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

      uiMEBits = uiBits[0];
    }
    else
    {
      uiLastMode = 1;
      pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMv( cMv[1], ePartSize, uiPartAddr, 0, iPartIdx );
      pcCU->getCUMvField(REF_PIC_LIST_1)->setAllRefIdx( iRefIdx[1], ePartSize, uiPartAddr, 0, iPartIdx );

      TempMv = cMv[1] - cMvPred[1][iRefIdx[1]];
      pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd    ( TempMv,                 ePartSize, uiPartAddr, 0, iPartIdx );

      pcCU->setInterDirSubParts( 2, uiPartAddr, iPartIdx, pcCU->getDepth(0) );

      pcCU->setMVPIdxSubParts( aaiMvpIdx[1][iRefIdx[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
      pcCU->setMVPNumSubParts( aaiMvpNum[1][iRefIdx[1]], REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));

      uiMEBits = uiBits[1];
    }
#if AMP_MRG
    } // end if bTestNormalMC
#endif

    if ( pcCU->getPartitionSize( uiPartAddr ) != SIZE_2Nx2N )
    {
      UInt uiMRGInterDir = 0;
      TComMvField cMRGMvField[2];
      UInt uiMRGIndex = 0;

      UInt uiMEInterDir = 0;
      TComMvField cMEMvField[2];

      m_pcRdCost->selectMotionLambda( true, 0, pcCU->getCUTransquantBypass(uiPartAddr) );

#if AMP_MRG
      // calculate ME cost
      Distortion uiMEError = std::numeric_limits<Distortion>::max();
      Distortion uiMECost  = std::numeric_limits<Distortion>::max();

      if (bTestNormalMC)
      {
        xGetInterPredictionError( pcCU, pcOrgYuv, iPartIdx, uiMEError, m_pcEncCfg->getUseHADME() );
        uiMECost = uiMEError + m_pcRdCost->getCost( uiMEBits );
      }
#else
      // calculate ME cost
      Distortion uiMEError = std::numeric_limits<Distortion>::max();
      xGetInterPredictionError( pcCU, pcOrgYuv, iPartIdx, uiMEError, m_pcEncCfg->getUseHADME() );
      Distortion uiMECost = uiMEError + m_pcRdCost->getCost( uiMEBits );
#endif
      // save ME result.
      uiMEInterDir = pcCU->getInterDir( uiPartAddr );
      TComDataCU::getMvField( pcCU, uiPartAddr, REF_PIC_LIST_0, cMEMvField[0] );
      TComDataCU::getMvField( pcCU, uiPartAddr, REF_PIC_LIST_1, cMEMvField[1] );

      // find Merge result
      Distortion uiMRGCost = std::numeric_limits<Distortion>::max();

      xMergeEstimation( pcCU, pcOrgYuv, iPartIdx, uiMRGInterDir, cMRGMvField, uiMRGIndex, uiMRGCost, cMvFieldNeighbours, uhInterDirNeighbours, numValidMergeCand);

      if ( uiMRGCost < uiMECost )
      {
        // set Merge result
        pcCU->setMergeFlagSubParts ( true,          uiPartAddr, iPartIdx, pcCU->getDepth( uiPartAddr ) );
        pcCU->setMergeIndexSubParts( uiMRGIndex,    uiPartAddr, iPartIdx, pcCU->getDepth( uiPartAddr ) );
        pcCU->setInterDirSubParts  ( uiMRGInterDir, uiPartAddr, iPartIdx, pcCU->getDepth( uiPartAddr ) );
        pcCU->getCUMvField( REF_PIC_LIST_0 )->setAllMvField( cMRGMvField[0], ePartSize, uiPartAddr, 0, iPartIdx );
        pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllMvField( cMRGMvField[1], ePartSize, uiPartAddr, 0, iPartIdx );

        pcCU->getCUMvField(REF_PIC_LIST_0)->setAllMvd    ( cMvZero,            ePartSize, uiPartAddr, 0, iPartIdx );
        pcCU->getCUMvField(REF_PIC_LIST_1)->setAllMvd    ( cMvZero,            ePartSize, uiPartAddr, 0, iPartIdx );

        pcCU->setMVPIdxSubParts( -1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
        pcCU->setMVPNumSubParts( -1, REF_PIC_LIST_0, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
        pcCU->setMVPIdxSubParts( -1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
        pcCU->setMVPNumSubParts( -1, REF_PIC_LIST_1, uiPartAddr, iPartIdx, pcCU->getDepth(uiPartAddr));
      }
      else
      {
        // set ME result
        pcCU->setMergeFlagSubParts( false,        uiPartAddr, iPartIdx, pcCU->getDepth( uiPartAddr ) );
        pcCU->setInterDirSubParts ( uiMEInterDir, uiPartAddr, iPartIdx, pcCU->getDepth( uiPartAddr ) );
        pcCU->getCUMvField( REF_PIC_LIST_0 )->setAllMvField( cMEMvField[0], ePartSize, uiPartAddr, 0, iPartIdx );
        pcCU->getCUMvField( REF_PIC_LIST_1 )->setAllMvField( cMEMvField[1], ePartSize, uiPartAddr, 0, iPartIdx );
      }
    }

    //  MC
    motionCompensation ( pcCU, pcPredYuv, REF_PIC_LIST_X, iPartIdx );

  } //  end of for ( Int iPartIdx = 0; iPartIdx < iNumPart; iPartIdx++ )

  setWpScalingDistParam( pcCU, -1, REF_PIC_LIST_X );
 // CTUW = iRoiWidth;
 // CTUH = iRoiHeight;
 // myfile << CTUH << ',' << CTUW << endl;
  return;
}


// AMVP
Void TEncSearch::xEstimateMvPredAMVP( TComDataCU* pcCU, TComYuv* pcOrgYuv, UInt uiPartIdx, RefPicList eRefPicList, Int iRefIdx, TComMv& rcMvPred, Bool bFilled, Distortion* puiDistBiP )
{
	
	//ofstream myfile;

	//myfile.open("C:\\RACEHOECES3_CTU.csv", ios::app);
  AMVPInfo*  pcAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAMVPInfo();

  TComMv     cBestMv;
  Int        iBestIdx   = 0;
  TComMv     cZeroMv;
  TComMv     cMvPred;
  Distortion uiBestCost = std::numeric_limits<Distortion>::max();
  UInt       uiPartAddr = 0;
  Int        iRoiWidth, iRoiHeight;
  Int        i;
 
  pcCU->getPartIndexAndSize( uiPartIdx, uiPartAddr, iRoiWidth, iRoiHeight );
 

  // Fill the MV Candidates
  if (!bFilled)
  {
    pcCU->fillMvpCand( uiPartIdx, uiPartAddr, eRefPicList, iRefIdx, pcAMVPInfo );
  }

  // initialize Mvp index & Mvp
  iBestIdx = 0;
  cBestMv  = pcAMVPInfo->m_acMvCand[0];
  if (pcAMVPInfo->iN <= 1)
  {
    rcMvPred = cBestMv;

    pcCU->setMVPIdxSubParts( iBestIdx, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
    pcCU->setMVPNumSubParts( pcAMVPInfo->iN, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));

    if(pcCU->getSlice()->getMvdL1ZeroFlag() && eRefPicList==REF_PIC_LIST_1)
    {
      (*puiDistBiP) = xGetTemplateCost( pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, rcMvPred, 0, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight);
    }
    return;
  }

  if (bFilled)
  {
    assert(pcCU->getMVPIdx(eRefPicList,uiPartAddr) >= 0);
    rcMvPred = pcAMVPInfo->m_acMvCand[pcCU->getMVPIdx(eRefPicList,uiPartAddr)];
    return;
  }

  m_cYuvPredTemp.clear();
  //-- Check Minimum Cost.
  for ( i = 0 ; i < pcAMVPInfo->iN; i++)
  {
    Distortion uiTmpCost;
    uiTmpCost = xGetTemplateCost( pcCU, uiPartAddr, pcOrgYuv, &m_cYuvPredTemp, pcAMVPInfo->m_acMvCand[i], i, AMVP_MAX_NUM_CANDS, eRefPicList, iRefIdx, iRoiWidth, iRoiHeight);
    if ( uiBestCost > uiTmpCost )
    {
      uiBestCost = uiTmpCost;
      cBestMv   = pcAMVPInfo->m_acMvCand[i];
      iBestIdx  = i;
      (*puiDistBiP) = uiTmpCost;
    }
  }

  m_cYuvPredTemp.clear();

  // Setting Best MVP
  rcMvPred = cBestMv;
  pcCU->setMVPIdxSubParts( iBestIdx, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
  pcCU->setMVPNumSubParts( pcAMVPInfo->iN, eRefPicList, uiPartAddr, uiPartIdx, pcCU->getDepth(uiPartAddr));
  CTUW2 = iRoiWidth;
  CTUH2 = iRoiHeight;

//  myfile << CTUH2 << ',' << CTUW2 << endl;
  return;
  
}

UInt TEncSearch::xGetMvpIdxBits(Int iIdx, Int iNum)
{
  assert(iIdx >= 0 && iNum >= 0 && iIdx < iNum);

  if (iNum == 1)
  {
    return 0;
  }

  UInt uiLength = 1;
  Int iTemp = iIdx;
  if ( iTemp == 0 )
  {
    return uiLength;
  }

  Bool bCodeLast = ( iNum-1 > iTemp );

  uiLength += (iTemp-1);

  if( bCodeLast )
  {
    uiLength++;
  }

  return uiLength;
}

Void TEncSearch::xGetBlkBits( PartSize eCUMode, Bool bPSlice, Int iPartIdx, UInt uiLastMode, UInt uiBlkBit[3])
{
  if ( eCUMode == SIZE_2Nx2N )
  {
    uiBlkBit[0] = (! bPSlice) ? 3 : 1;
    uiBlkBit[1] = 3;
    uiBlkBit[2] = 5;
  }
  else if ( (eCUMode == SIZE_2NxN || eCUMode == SIZE_2NxnU) || eCUMode == SIZE_2NxnD )
  {
    UInt aauiMbBits[2][3][3] = { { {0,0,3}, {0,0,0}, {0,0,0} } , { {5,7,7}, {7,5,7}, {9-3,9-3,9-3} } };
    if ( bPSlice )
    {
      uiBlkBit[0] = 3;
      uiBlkBit[1] = 0;
      uiBlkBit[2] = 0;
    }
    else
    {
      ::memcpy( uiBlkBit, aauiMbBits[iPartIdx][uiLastMode], 3*sizeof(UInt) );
    }
  }
  else if ( (eCUMode == SIZE_Nx2N || eCUMode == SIZE_nLx2N) || eCUMode == SIZE_nRx2N )
  {
    UInt aauiMbBits[2][3][3] = { { {0,2,3}, {0,0,0}, {0,0,0} } , { {5,7,7}, {7-2,7-2,9-2}, {9-3,9-3,9-3} } };
    if ( bPSlice )
    {
      uiBlkBit[0] = 3;
      uiBlkBit[1] = 0;
      uiBlkBit[2] = 0;
    }
    else
    {
      ::memcpy( uiBlkBit, aauiMbBits[iPartIdx][uiLastMode], 3*sizeof(UInt) );
    }
  }
  else if ( eCUMode == SIZE_NxN )
  {
    uiBlkBit[0] = (! bPSlice) ? 3 : 1;
    uiBlkBit[1] = 3;
    uiBlkBit[2] = 5;
  }
  else
  {
    printf("Wrong!\n");
    assert( 0 );
  }
}

Void TEncSearch::xCopyAMVPInfo (AMVPInfo* pSrc, AMVPInfo* pDst)
{
  pDst->iN = pSrc->iN;
  for (Int i = 0; i < pSrc->iN; i++)
  {
    pDst->m_acMvCand[i] = pSrc->m_acMvCand[i];
  }
}

Void TEncSearch::xCheckBestMVP ( TComDataCU* pcCU, RefPicList eRefPicList, TComMv cMv, TComMv& rcMvPred, Int& riMVPIdx, UInt& ruiBits, Distortion& ruiCost )
{
  AMVPInfo* pcAMVPInfo = pcCU->getCUMvField(eRefPicList)->getAMVPInfo();
  
  assert(pcAMVPInfo->m_acMvCand[riMVPIdx] == rcMvPred);

  if (pcAMVPInfo->iN < 2)
  {
    return;
  }

  m_pcRdCost->selectMotionLambda( true, 0, pcCU->getCUTransquantBypass(0) );
  m_pcRdCost->setCostScale ( 0    );

  Int iBestMVPIdx = riMVPIdx;

  m_pcRdCost->setPredictor( rcMvPred );
  Int iOrgMvBits  = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer());
  iOrgMvBits += m_auiMVPIdxCost[riMVPIdx][AMVP_MAX_NUM_CANDS];
  Int iBestMvBits = iOrgMvBits;

  for (Int iMVPIdx = 0; iMVPIdx < pcAMVPInfo->iN; iMVPIdx++)
  {
    if (iMVPIdx == riMVPIdx)
    {
      continue;
    }

    m_pcRdCost->setPredictor( pcAMVPInfo->m_acMvCand[iMVPIdx] );

    Int iMvBits = m_pcRdCost->getBitsOfVectorWithPredictor(cMv.getHor(), cMv.getVer());
    iMvBits += m_auiMVPIdxCost[iMVPIdx][AMVP_MAX_NUM_CANDS];

    if (iMvBits < iBestMvBits)
    {
      iBestMvBits = iMvBits;
      iBestMVPIdx = iMVPIdx;
    }
  }

  if (iBestMVPIdx != riMVPIdx)  //if changed
  {
    rcMvPred = pcAMVPInfo->m_acMvCand[iBestMVPIdx];

    riMVPIdx = iBestMVPIdx;
    UInt uiOrgBits = ruiBits;
    ruiBits = uiOrgBits - iOrgMvBits + iBestMvBits;
    ruiCost = (ruiCost - m_pcRdCost->getCost( uiOrgBits ))  + m_pcRdCost->getCost( ruiBits );
  }
  
}


Distortion TEncSearch::xGetTemplateCost( TComDataCU* pcCU,
                                         UInt        uiPartAddr,
                                         TComYuv*    pcOrgYuv,
                                         TComYuv*    pcTemplateCand,
                                         TComMv      cMvCand,
                                         Int         iMVPIdx,
                                         Int         iMVPNum,
                                         RefPicList  eRefPicList,
                                         Int         iRefIdx,
                                         Int         iSizeX,
                                         Int         iSizeY
                                         )
{
  Distortion uiCost = std::numeric_limits<Distortion>::max();

  TComPicYuv* pcPicYuvRef = pcCU->getSlice()->getRefPic( eRefPicList, iRefIdx )->getPicYuvRec();

  pcCU->clipMv( cMvCand );

  // prediction pattern
  if ( pcCU->getSlice()->testWeightPred() && pcCU->getSlice()->getSliceType()==P_SLICE )
  {
    xPredInterBlk( COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, &cMvCand, iSizeX, iSizeY, pcTemplateCand, true, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) );
  }
  else
  {
    xPredInterBlk( COMPONENT_Y, pcCU, pcPicYuvRef, uiPartAddr, &cMvCand, iSizeX, iSizeY, pcTemplateCand, false, pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) );
  }

  if ( pcCU->getSlice()->testWeightPred() && pcCU->getSlice()->getSliceType()==P_SLICE )
  {
    xWeightedPredictionUni( pcCU, pcTemplateCand, uiPartAddr, iSizeX, iSizeY, eRefPicList, pcTemplateCand, iRefIdx );
  }

  // calc distortion

  uiCost = m_pcRdCost->getDistPart( pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA), pcTemplateCand->getAddr(COMPONENT_Y, uiPartAddr), pcTemplateCand->getStride(COMPONENT_Y), pcOrgYuv->getAddr(COMPONENT_Y, uiPartAddr), pcOrgYuv->getStride(COMPONENT_Y), iSizeX, iSizeY, COMPONENT_Y, DF_SAD );
  uiCost = (UInt) m_pcRdCost->calcRdCost( m_auiMVPIdxCost[iMVPIdx][iMVPNum], uiCost, DF_SAD );
  return uiCost;
}


Void TEncSearch::xMotionEstimation( TComDataCU* pcCU, TComYuv* pcYuvOrg, Int iPartIdx, RefPicList eRefPicList, TComMv* pcMvPred, Int iRefIdxPred, TComMv& rcMv, UInt& ruiBits, Distortion& ruiCost, Bool bBi  )
{
  UInt          uiPartAddr;
  Int           iRoiWidth;
  Int           iRoiHeight;

  TComMv        cMvHalf, cMvQter;
  TComMv        cMvSrchRngLT;
  TComMv        cMvSrchRngRB;

  TComYuv*      pcYuv = pcYuvOrg;
 
 //ofstream myfile;
 // ofstream myfile2;
  //ofstream myfile3;
 // myfile.open("C:\\FB_CTU.csv", ios::app);
//  myfile2.open("C:\\BLOWINGBUBBLES_MV.csv", ios::app);
  
  
  assert(eRefPicList < MAX_NUM_REF_LIST_ADAPT_SR && iRefIdxPred<Int(MAX_IDX_ADAPT_SR));
  m_iSearchRange = m_aaiAdaptSR[eRefPicList][iRefIdxPred];

  Int           iSrchRng      = ( bBi ? m_bipredSearchRange : m_iSearchRange );
  TComPattern   tmpPattern;
  TComPattern*  pcPatternKey  = &tmpPattern;

  Double        fWeight       = 1.0;

  pcCU->getPartIndexAndSize( iPartIdx, uiPartAddr, iRoiWidth, iRoiHeight );
  CH = iRoiHeight;
  CW = iRoiWidth;
 //myfile << CH << ',' << CW << endl;

  if ( bBi ) // Bipredictive ME
  {
    TComYuv*  pcYuvOther = &m_acYuvPred[1-(Int)eRefPicList];
    pcYuv                = &m_cYuvPredTemp;

    pcYuvOrg->copyPartToPartYuv( pcYuv, uiPartAddr, iRoiWidth, iRoiHeight );

    pcYuv->removeHighFreq( pcYuvOther, uiPartAddr, iRoiWidth, iRoiHeight, pcCU->getSlice()->getSPS()->getBitDepths().recon, m_pcEncCfg->getClipForBiPredMeEnabled() );

    fWeight = 0.5;
  }
  m_cDistParam.bIsBiPred = bBi;

  //  Search key pattern initialization
  pcPatternKey->initPattern( pcYuv->getAddr  ( COMPONENT_Y, uiPartAddr ),
                             iRoiWidth,
                             iRoiHeight,
                             pcYuv->getStride(COMPONENT_Y),
                             pcCU->getSlice()->getSPS()->getBitDepth(CHANNEL_TYPE_LUMA) );

  Pel*        piRefY      = pcCU->getSlice()->getRefPic( eRefPicList, iRefIdxPred )->getPicYuvRec()->getAddr( COMPONENT_Y, pcCU->getCtuRsAddr(), pcCU->getZorderIdxInCtu() + uiPartAddr );
  Int         iRefStride  = pcCU->getSlice()->getRefPic( eRefPicList, iRefIdxPred )->getPicYuvRec()->getStride(COMPONENT_Y);

  TComMv      cMvPred = *pcMvPred;

  if ( bBi )
  {
	  
    xSetSearchRange   ( pcCU, rcMv   , iSrchRng, cMvSrchRngLT, cMvSrchRngRB );
  }
  else
  {
	  
    xSetSearchRange   ( pcCU, cMvPred, iSrchRng, cMvSrchRngLT, cMvSrchRngRB );
  }

  m_pcRdCost->selectMotionLambda(true, 0, pcCU->getCUTransquantBypass(uiPartAddr) );

  m_pcRdCost->setPredictor  ( *pcMvPred );
  m_pcRdCost->setCostScale  ( 2 );

  setWpScalingDistParam( pcCU, iRefIdxPred, eRefPicList );
  //  Do integer search
  if ( (m_motionEstimationSearchMethod==MESEARCH_FULL) || bBi )
  {
    xPatternSearch      ( pcPatternKey, piRefY, iRefStride, &cMvSrchRngLT, &cMvSrchRngRB, rcMv, ruiCost );
  }
  else
  {
    rcMv = *pcMvPred;
    const TComMv *pIntegerMv2Nx2NPred=0;
    if (pcCU->getPartitionSize(0) != SIZE_2Nx2N || pcCU->getDepth(0) != 0)
    {
      pIntegerMv2Nx2NPred = &(m_integerMv2Nx2N[eRefPicList][iRefIdxPred]);
    }
    xPatternSearchFast  ( pcCU, pcPatternKey, piRefY, iRefStride, &cMvSrchRngLT, &cMvSrchRngRB, rcMv, ruiCost, pIntegerMv2Nx2NPred );
    if (pcCU->getPartitionSize(0) == SIZE_2Nx2N)
    {
      m_integerMv2Nx2N[eRefPicList][iRefIdxPred] = rcMv;
    }
  }

  m_pcRdCost->selectMotionLambda( true, 0, pcCU->getCUTransquantBypass(uiPartAddr) );
  m_pcRdCost->setCostScale ( 1 );

  const Bool bIsLosslessCoded = pcCU->getCUTransquantBypass(uiPartAddr) != 0;
  xPatternSearchFracDIF( bIsLosslessCoded, pcPatternKey, piRefY, iRefStride, &rcMv, cMvHalf, cMvQter, ruiCost );

  m_pcRdCost->setCostScale( 0 );

  

 //myfile2 << rcMv.getHor() << ',' << rcMv.getVer() << ',' << cMvHalf.getHor() << ',' << cMvQter.getHor() << ',' << cMvHalf.getVer()<<','<<cMvQter.getVer() << endl;

  TComMv MV_HALF;

  TComMv MV_QRTER;

   
 MV_HALF.setHor(MVX_HALF);
 MV_HALF.setVer(MVY_HALF);
 MV_QRTER.setHor(MVX_QRTER);
 MV_QRTER.setVer(MVY_QRTER);
 
// ofstream myfile2;

  //myfile2.open("C:\\FB_MV.csv", ios::app);
// myfile2 << MV_HALF.getHor() << ',' << MV_QRTER.getHor() << ',' << MV_HALF.getVer() << ',' << MV_QRTER.getVer() << ',' << cMvHalf.getHor() << ',' << cMvQter.getHor() << ',' << cMvHalf.getVer() << ',' << cMvQter.getVer() << ','<<rcMv.getHor()<<','<<rcMv.getVer() <<endl;

  rcMv <<= 2;

  //rcMv += (MV_HALF <<= 1);
 //rcMv += MV_QRTER;
  rcMv += (cMvHalf <<= 1);
 rcMv += cMvQter;

  UInt uiMvBits = m_pcRdCost->getBitsOfVectorWithPredictor( rcMv.getHor(), rcMv.getVer() );

  ruiBits      += uiMvBits;
  ruiCost       = (Distortion)( floor( fWeight * ( (Double)ruiCost - (Double)m_pcRdCost->getCost( uiMvBits ) ) ) + (Double)m_pcRdCost->getCost( ruiBits ) );
}


Void TEncSearch::xSetSearchRange ( const TComDataCU* const pcCU, const TComMv& cMvPred, const Int iSrchRng,
                                   TComMv& rcMvSrchRngLT, TComMv& rcMvSrchRngRB )
{
  Int  iMvShift = 2;
  TComMv cTmpMvPred = cMvPred;
  pcCU->clipMv( cTmpMvPred );

  rcMvSrchRngLT.setHor( cTmpMvPred.getHor() - (iSrchRng << iMvShift) );
  rcMvSrchRngLT.setVer( cTmpMvPred.getVer() - (iSrchRng << iMvShift) );

  rcMvSrchRngRB.setHor( cTmpMvPred.getHor() + (iSrchRng << iMvShift) );
  rcMvSrchRngRB.setVer( cTmpMvPred.getVer() + (iSrchRng << iMvShift) );
  pcCU->clipMv        ( rcMvSrchRngLT );
  pcCU->clipMv        ( rcMvSrchRngRB );

#if ME_ENABLE_ROUNDING_OF_MVS
  rcMvSrchRngLT.divideByPowerOf2(iMvShift);
  rcMvSrchRngRB.divideByPowerOf2(iMvShift);
#else
  rcMvSrchRngLT >>= iMvShift;
  rcMvSrchRngRB >>= iMvShift;
#endif
}


Void TEncSearch::xPatternSearch(const TComPattern* const pcPatternKey,
	const Pel*               piRefY,
	const Int                iRefStride,
	const TComMv* const      pcMvSrchRngLT,
	const TComMv* const      pcMvSrchRngRB,
	TComMv&      rcMv,
	Distortion&  ruiSAD)
{
	Int   iSrchRngHorLeft = pcMvSrchRngLT->getHor();
	Int   iSrchRngHorRight = pcMvSrchRngRB->getHor();
	Int   iSrchRngVerTop = pcMvSrchRngLT->getVer();
	Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

	Distortion  uiSad;
	Distortion  uiSadBest = std::numeric_limits<Distortion>::max();
	Int         iBestX = 0;
	Int         iBestY = 0;

	



	

	//////////modified variables///////////
	int dimHor = 0;
	int dimVer = 0;
	//-- jclee for using the SAD function pointer
	m_pcRdCost->setDistParam(pcPatternKey, piRefY, iRefStride, m_cDistParam);

	// fast encoder decision: use subsampled SAD for integer ME
	if (m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE1 || m_pcEncCfg->getFastInterSearchMode() == FASTINTERSEARCH_MODE3)
	{
		if (m_cDistParam.iRows > 8)
		{
			m_cDistParam.iSubShift = 1;
		}
	}

	piRefY += (iSrchRngVerTop * iRefStride);

	



	int Distance_x;
	int Distance_y;


	if (iSrchRngHorLeft <= 0 && iSrchRngHorRight>0)
	{
		Distance_x = ((abs(iSrchRngHorLeft) + abs(iSrchRngHorRight)) / 2);
		
	}
	else if (iSrchRngHorLeft <= 0 && iSrchRngHorRight <= 0)
	{
		Distance_x = ((abs(iSrchRngHorLeft) - abs(iSrchRngHorRight)) / 2);
		
	}

	else if (iSrchRngHorLeft >= 0 && iSrchRngHorRight >= 0)
	{

		Distance_x = ((abs(iSrchRngHorRight) - abs(iSrchRngHorLeft)) / 2);
		

	}





	if (iSrchRngVerTop <= 0 && iSrchRngVerBottom > 0)
	{
		Distance_y = ((abs(iSrchRngVerTop) + abs(iSrchRngVerBottom)) / 2);
		
	}
	else if (iSrchRngVerTop <= 0 && iSrchRngVerBottom <= 0)
	{
		Distance_y = ((abs(iSrchRngVerTop) - abs(iSrchRngVerBottom)) / 2);
		
	}

	else if (iSrchRngVerTop >= 0 && iSrchRngVerBottom >= 0)
	{

		Distance_y = ((abs(iSrchRngVerBottom) - abs(iSrchRngVerTop)) / 2);
	

	}



	/////////////modification_Here////////////

	if (((iSrchRngVerBottom - iSrchRngVerTop) % 2 == 0) && ((iSrchRngHorRight - iSrchRngHorLeft) % 2 == 0))
	{

		dimHor = ((2 * Distance_x) + 1);
		dimVer = ((2 * Distance_y) + 1);

	}

	else if (((iSrchRngVerBottom - iSrchRngVerTop) % 2 == 0) && ((iSrchRngHorRight - iSrchRngHorLeft) % 2 != 0))
	{
		dimHor = ((2 * Distance_x) + 2);
		dimVer = ((2 * Distance_y) + 1);

	}

	else if (((iSrchRngVerBottom - iSrchRngVerTop) % 2 != 0) && ((iSrchRngHorRight - iSrchRngHorLeft) % 2 == 0))
	{
		dimHor = ((2 * Distance_x) + 1);
		dimVer = ((2 * Distance_y) + 2);

	}

	else if (((iSrchRngVerBottom - iSrchRngVerTop) % 2 != 0) && ((iSrchRngHorRight - iSrchRngHorLeft) % 2 != 0))
	{
		dimHor = ((2 * Distance_x) + 2);
		dimVer = ((2 * Distance_y) + 2);

	}

	/////////end of modification////////////////





	/////////////////modifying//////////////

	int shift_y;
	int shift_x;

	if (abs(iSrchRngVerTop) == abs(iSrchRngVerBottom) && abs(iSrchRngHorLeft) == abs(iSrchRngHorRight))
	{
		shift_y = Distance_y;
		shift_x = Distance_x;
	}

	else

	{

		if (iSrchRngVerTop <= 0 && abs(iSrchRngHorLeft) == abs(iSrchRngHorRight))


		{
			shift_y = abs(iSrchRngVerTop);
			shift_x = Distance_x;
		}

		else if (iSrchRngVerTop > 0 && abs(iSrchRngHorLeft) == abs(iSrchRngHorRight))

		{
			shift_y = -abs(iSrchRngVerTop);
			shift_x = Distance_x;
		}

		else if (iSrchRngHorLeft <= 0 && abs(iSrchRngVerTop) == abs(iSrchRngVerBottom))
		{
			shift_x = abs(iSrchRngHorLeft);
			shift_y = Distance_y;
		}

		else if (iSrchRngHorLeft > 0 && abs(iSrchRngVerTop) == abs(iSrchRngVerBottom))
		{

			shift_x = -abs(iSrchRngHorLeft);
			shift_y = Distance_y;
		}

		else if (iSrchRngHorLeft > 0 && iSrchRngVerTop > 0)
		{

			shift_x = -abs(iSrchRngHorLeft);
			shift_y = -abs(iSrchRngVerTop);
		}

		else if (iSrchRngHorLeft <= 0 && iSrchRngVerTop <= 0)

		{
			shift_x = abs(iSrchRngHorLeft);
			shift_y = abs(iSrchRngVerTop);

		}

		else if (iSrchRngHorLeft <= 0 && iSrchRngVerTop > 0)

		{
			shift_x = abs(iSrchRngHorLeft);
			shift_y = -abs(iSrchRngVerTop);

		}

		else if (iSrchRngHorLeft > 0 && iSrchRngVerTop <= 0)

		{
			shift_x = -abs(iSrchRngHorLeft);
			shift_y = abs(iSrchRngVerTop);

		}




	}




	///////////////end of modifying////////////////






	//	iSrchRngVerTop = -Distance_y;
	//  iSrchRngVerBottom = Distance_y;

	//	iSrchRngHorLeft = -Distance_x;
	//iSrchRngHorRight = Distance_x;


	int **dynamicArray = 0;
	dynamicArray = new int *[dimVer];
	for (int i = 0; i < dimVer; i++)
		dynamicArray[i] = new int[dimHor];

	//cout << "\nshift_y" << shift_y;
	//cout << "      shift_x " << shift_x;



	for (Int y = iSrchRngVerTop; y <= iSrchRngVerBottom; y++)
	{
		for (Int x = iSrchRngHorLeft; x <= iSrchRngHorRight; x++)
		{



			// end of decleration

			//  find min. distortion position
			m_cDistParam.pCur = piRefY + x;

			setDistParamComp(COMPONENT_Y);

			m_cDistParam.bitDepth = pcPatternKey->getBitDepthY();
			uiSad = m_cDistParam.DistFunc(&m_cDistParam);

			//cout << "\nX=" << x;
			//cout << "    y=" << y;


			// motion cost
			uiSad += m_pcRdCost->getCostOfVectorWithPredictor(x, y);






			dynamicArray[y + shift_y][x + shift_x] = uiSad;





			if (uiSad < uiSadBest)
			{
				uiSadBest = uiSad;
				iBestX = x;
				iBestY = y;
				m_cDistParam.m_maximumDistortionForEarlyExit = uiSad;
			}
		}
		piRefY += iRefStride;
	}

	//cout << "\n\nbestsad     " << uiSadBest;
	//cout << "\nXbest= " << iBestX;
	//cout << " Ybest  " << iBestY;


	// initiate the SAD window

	///// modifying to get the error window




	if ((((iBestX - 4) >= iSrchRngHorLeft) && ((iBestX + 4) <= iSrchRngHorRight)) && (((iBestY - 4) >= iSrchRngVerTop) && ((iBestY + 4) <= iSrchRngVerBottom)))
		 

	
	{





		// horizontal axis
		float A00 = dynamicArray[iBestY + shift_y][iBestX + shift_x];
		float A10 = dynamicArray[iBestY + shift_y][iBestX + shift_x + 1];
		float A20 = dynamicArray[iBestY + shift_y][iBestX + shift_x + 2];
		
		float AN10 = dynamicArray[iBestY + shift_y][iBestX + shift_x - 1];
		float AN20 = dynamicArray[iBestY + shift_y][iBestX + shift_x - 2];
		

		// vertical axis

		float A01 = dynamicArray[iBestY + shift_y + 1][iBestX + shift_x];
		float A02 = dynamicArray[iBestY + shift_y + 2][iBestX + shift_x];
		
		float A0N1 = dynamicArray[iBestY + shift_y - 1][iBestX + shift_x];
		float A0N2 = dynamicArray[iBestY + shift_y - 2][iBestX + shift_x];
		

		// u1,u2 corner

		float A1N1 = dynamicArray[iBestY + shift_y - 1][iBestX + shift_x + 1];
		float A2N1 = dynamicArray[iBestY + shift_y - 1][iBestX + shift_x + 2];
		
		float AN1N1 = dynamicArray[iBestY + shift_y - 1][iBestX + shift_x - 1];
		float AN2N1 = dynamicArray[iBestY + shift_y - 1][iBestX + shift_x - 2];
		


		float A1N2 = dynamicArray[iBestY + shift_y - 2][iBestX + shift_x + 1];
		float A2N2 = dynamicArray[iBestY + shift_y - 2][iBestX + shift_x + 2];
		
		float AN1N2 = dynamicArray[iBestY + shift_y - 2][iBestX + shift_x - 1];
		float AN2N2 = dynamicArray[iBestY + shift_y - 2][iBestX + shift_x - 2];
		


		


		// U3,U4 Corner


		float A11 = dynamicArray[iBestY + shift_y + 1][iBestX + shift_x + 1];
		float A21 = dynamicArray[iBestY + shift_y + 1][iBestX + shift_x + 2];
		
		float AN11 = dynamicArray[iBestY + shift_y + 1][iBestX + shift_x - 1];
		float AN21 = dynamicArray[iBestY + shift_y + 1][iBestX + shift_x - 2];
		


		float A12 = dynamicArray[iBestY + shift_y + 2][iBestX + shift_x + 1];
		float A22 = dynamicArray[iBestY + shift_y + 2][iBestX + shift_x + 2];
		
		float AN12 = dynamicArray[iBestY + shift_y + 2][iBestX + shift_x - 1];
		float AN22 = dynamicArray[iBestY + shift_y + 2][iBestX + shift_x - 2];
		


	

		

		// final result 

		if (MT_VCX == 0)

		{

			MVX_HALF = 0;
			MVX_QRTER = 0;

		}


		else if (MT_VCX == 0.25)

		{



			MVX_HALF = 0;
			MVX_QRTER = 1;

		}

		else if (MT_VCX == 0.5)

		{


			MVX_HALF = 1;
			MVX_QRTER = 0;

		}
		else if (MT_VCX == 0.75)

		{


			MVX_HALF = 1;
			MVX_QRTER = 1;

		}

		else if (MT_VCX == -0.5)

		{


			MVX_HALF = -1;
			MVX_QRTER = 0;

		}

		else if (MT_VCX == -0.25)

		{


			MVX_HALF = 0;
			MVX_QRTER = -1;

		}

		else if (MT_VCX == -0.75)

		{


			MVX_HALF = -1;
			MVX_QRTER = -1;

		}

		//MVY
		if (MT_VCY == 0)

		{

			MVY_HALF = 0;
			MVY_QRTER = 0;

		}


		else if (MT_VCY == 0.25)

		{



			MVY_HALF = 0;
			MVY_QRTER = 1;

		}

		else if (MT_VCY == 0.5)

		{

			MVY_HALF = 1;
			MVY_QRTER = 0;

		}
		else if (MT_VCY == 0.75)

		{

			MVY_HALF = 1;
			MVY_QRTER = 1;

		}


		else if (MT_VCY == -0.5)

		{


			MVY_HALF = -1;
			MVY_QRTER = 0;

		}

		else if (MT_VCY == -0.25)

		{


			MVY_HALF = 0;
			MVY_QRTER = -1;

		}
		else if (MT_VCY == -0.75)

		{


			MVY_HALF = -1;
			MVY_QRTER = -1;

		}


		


	}


	





	  
		  
			 



		 






	  


	  


		// end of modifying the SAD Window


		for (int i = 0; i < dimVer; i++)
			delete[] dynamicArray[i];
		delete[] dynamicArray;


		rcMv.set(iBestX, iBestY);


		ruiSAD = uiSadBest - m_pcRdCost->getCostOfVectorWithPredictor(iBestX, iBestY);

		//getchar();
		return;
	}


Void TEncSearch::xPatternSearchFast( const TComDataCU* const  pcCU,
                                     const TComPattern* const pcPatternKey,
                                     const Pel* const         piRefY,
                                     const Int                iRefStride,
                                     const TComMv* const      pcMvSrchRngLT,
                                     const TComMv* const      pcMvSrchRngRB,
                                     TComMv&                  rcMv,
                                     Distortion&              ruiSAD,
                                     const TComMv* const      pIntegerMv2Nx2NPred )
{
	//ofstream myfile;
	//ofstream myfile2;
	//myfile.open("C:\\FB_SSE.csv", ios::app);
	//myfile2.open("C:\\FB_Flag.csv", ios::app);
  assert (MD_LEFT < NUM_MV_PREDICTORS);
  pcCU->getMvPredLeft       ( m_acMvPredictors[MD_LEFT] );
  assert (MD_ABOVE < NUM_MV_PREDICTORS);
  pcCU->getMvPredAbove      ( m_acMvPredictors[MD_ABOVE] );
  assert (MD_ABOVE_RIGHT < NUM_MV_PREDICTORS);
  pcCU->getMvPredAboveRight ( m_acMvPredictors[MD_ABOVE_RIGHT] );

  switch ( m_motionEstimationSearchMethod )
  {
    case MESEARCH_DIAMOND:
      xTZSearch( pcCU, pcPatternKey, piRefY, iRefStride, pcMvSrchRngLT, pcMvSrchRngRB, rcMv, ruiSAD, pIntegerMv2Nx2NPred, false );
	  
	  C = array[0];
	  for (int i = 1; i <=index_ref - 1; i++)
	  {
		  if (array[i] < C)
			  C = array[i];

	  }
	 
	 // index_ref = index_ref + 1;
	  U1 = array[index_ref];
	  V1 = array[index_ref + 1];
	  U2 = array[index_ref + 2];
	  H1 = array[index_ref + 3];
	  
	  H2 = array[index_ref + 4];
	  U3 = array[index_ref + 5];
	  V2 = array[index_ref + 6];
	  U4 = array[index_ref + 7];
	 
	 
	  
	   A00 = C;
	   A10 = H2;
	   A20 = array[index_ref + 16];

	  AN10 = H1;
	  AN20 = array[index_ref + 15];


	  // vertical axis

	   A01 = V2;
	   A02 = array[index_ref + 21];

	   A0N1 = V1;
	   A0N2 = array[index_ref + 10];


	  // u1,u2 corner

	   A1N1 = U2;
	   A2N1 = array[index_ref + 14];

	   AN1N1 = U1;
	   AN2N1 = array[index_ref + 13];



	   A1N2 = array[index_ref + 11];
	   A2N2 = array[index_ref + 12];

	   AN1N2 = array[index_ref + 9];
	   AN2N2 = array[index_ref + 8];






	  // U3,U4 Corner


	   A11 = U4;
	   A21 = array[index_ref + 18];

	   AN11 = U3;
	   AN21 = array[index_ref + 17];



	   A12 = array[index_ref + 22];
	   A22 = array[index_ref + 23];

	  AN12 = array[index_ref + 20];
	  AN22 = array[index_ref + 19];
	  
	  


	  // neural network implementation
	  
	  
	  IN1 = U1;
	  IN2 = V1;
	  IN3 = U2;
	  IN4 = H1;
	  IN5 = C;
	  IN6 = H2;
	  IN7 = U3;
	  IN8 = V2;
	  IN9 = U4;

	  
	  // hidden layer

	  for(int i=0;i<9;i++){
		IN[i] = IN[i] * varin[i];	  
	  }
	  
	  for(int i=0;i<22;i++){
		for(int j=0;j<9;j++){
			X1[i] += in_h1[i][j] * IN[j];
		}
		X1[i] += b1[i];
		X1[i] = relu(X1[i]) * var1[i];
	  }

	  for(int i=0;i<20;i++){
		for(int j=0;j<22;j++){
			X2[i] += h1_h2[i][j] * X1[j];
		}
		X2[i] += b2[i];
		X2[i] = relu(X2[i]) * var2[i];
	  }	  
	  
	  
	  // OUTPUT LAYER
	  
	  for(int i=0;i<49;i++){
		for(int j=0;j<20;j++){
			OUT[i] += h2_out[i][j] * X2[j];
		}
		OUT[i] += bout[i];
		OUT[i] = sigmoid(OUT[i]);
	  }	  
	

	 
//1
	 
if((OUT[0]>OUT[1])&&(OUT[0]>OUT[2])&&(OUT[0]>OUT[3])&&(OUT[0]>OUT[4])&&(OUT[0]>OUT[5])&&(OUT[0]>OUT[6])&&(OUT[0]>OUT[7])&&(OUT[0]>OUT[8])&&(OUT[0]>OUT[9])&&(OUT[0]>OUT[10])
&&(OUT[0]>OUT[11])&&(OUT[0]>OUT[12])&&(OUT[0]>OUT[13])&&(OUT[0]>OUT[14])&&(OUT[0]>OUT[15])&&(OUT[0]>OUT[16])&&(OUT[0]>OUT[17])&&(OUT[0]>OUT[18])&&(OUT[0]>OUT[19])
&&(OUT[0]>OUT[20])&&(OUT[0]>OUT[21])&&(OUT[0]>OUT[22])&&(OUT[0]>OUT[23])&&(OUT[0]>OUT[24])&&(OUT[0]>OUT[25])&&(OUT[0]>OUT[26])&&(OUT[0]>OUT[27])&&(OUT[0]>OUT[28])
&&(OUT[0]>OUT[29])&&(OUT[0]>OUT[30])&&(OUT[0]>OUT[31])&&(OUT[0]>OUT[32])&&(OUT[0]>OUT[33])&&(OUT[0]>OUT[34])&&(OUT[0]>OUT[35])&&(OUT[0]>OUT[36])&&(OUT[0]>OUT[37])
&&(OUT[0]>OUT[38])&&(OUT[0]>OUT[39])&&(OUT[0]>OUT[40])&&(OUT[0]>OUT[41])&&(OUT[0]>OUT[42])&&(OUT[0]>OUT[43])&&(OUT[0]>OUT[44])&&(OUT[0]>OUT[45])&&(OUT[0]>OUT[46])
&&(OUT[0]>OUT[47])&&(OUT[0]>OUT[48]) ) 
	 
{
MT_VCX=-0.75;
MT_VCY=-0.75;		 
		 
}
		 
//2
	 
	else if((OUT[1]>OUT[0])&&(OUT[1]>OUT[2])&&(OUT[1]>OUT[3])&&(OUT[1]>OUT[4])&&(OUT[1]>OUT[5])&&(OUT[1]>OUT[6])&&(OUT[1]>OUT[7])&&(OUT[1]>OUT[8])&&(OUT[1]>OUT[9])&&(OUT[1]>OUT[10])
&&(OUT[1]>OUT[11])&&(OUT[1]>OUT[12])&&(OUT[1]>OUT[13])&&(OUT[1]>OUT[14])&&(OUT[1]>OUT[15])&&(OUT[1]>OUT[16])&&(OUT[1]>OUT[17])&&(OUT[1]>OUT[18])&&(OUT[1]>OUT[19])
&&(OUT[1]>OUT[20])&&(OUT[1]>OUT[21])&&(OUT[1]>OUT[22])&&(OUT[1]>OUT[23])&&(OUT[1]>OUT[24])&&(OUT[1]>OUT[25])&&(OUT[1]>OUT[26])&&(OUT[1]>OUT[27])&&(OUT[1]>OUT[28])
&&(OUT[1]>OUT[29])&&(OUT[1]>OUT[30])&&(OUT[1]>OUT[31])&&(OUT[1]>OUT[32])&&(OUT[1]>OUT[33])&&(OUT[1]>OUT[34])&&(OUT[1]>OUT[35])&&(OUT[1]>OUT[36])&&(OUT[1]>OUT[37])
&&(OUT[1]>OUT[38])&&(OUT[1]>OUT[39])&&(OUT[1]>OUT[40])&&(OUT[1]>OUT[41])&&(OUT[1]>OUT[42])&&(OUT[1]>OUT[43])&&(OUT[1]>OUT[44])&&(OUT[1]>OUT[45])&&(OUT[1]>OUT[46])
&&(OUT[1]>OUT[47])&&(OUT[1]>OUT[48]) ) 
{
MT_VCX=-0.5;
MT_VCY=-0.75;		
	
}

//3

else if((OUT[2]>OUT[0])&&(OUT[2]>OUT[1])&&(OUT[2]>OUT[3])&&(OUT[2]>OUT[4])&&(OUT[2]>OUT[5])&&(OUT[2]>OUT[6])&&(OUT[2]>OUT[7])&&(OUT[2]>OUT[8])&&(OUT[2]>OUT[9])&&(OUT[2]>OUT[10])
&&(OUT[2]>OUT[11])&&(OUT[2]>OUT[12])&&(OUT[2]>OUT[13])&&(OUT[2]>OUT[14])&&(OUT[2]>OUT[15])&&(OUT[2]>OUT[16])&&(OUT[2]>OUT[17])&&(OUT[2]>OUT[18])&&(OUT[2]>OUT[19])
&&(OUT[2]>OUT[20])&&(OUT[2]>OUT[21])&&(OUT[2]>OUT[22])&&(OUT[2]>OUT[23])&&(OUT[2]>OUT[24])&&(OUT[2]>OUT[25])&&(OUT[2]>OUT[26])&&(OUT[2]>OUT[27])&&(OUT[2]>OUT[28])
&&(OUT[2]>OUT[29])&&(OUT[2]>OUT[30])&&(OUT[2]>OUT[31])&&(OUT[2]>OUT[32])&&(OUT[2]>OUT[33])&&(OUT[2]>OUT[34])&&(OUT[2]>OUT[35])&&(OUT[2]>OUT[36])&&(OUT[2]>OUT[37])
&&(OUT[2]>OUT[38])&&(OUT[2]>OUT[39])&&(OUT[2]>OUT[40])&&(OUT[2]>OUT[41])&&(OUT[2]>OUT[42])&&(OUT[2]>OUT[43])&&(OUT[2]>OUT[44])&&(OUT[2]>OUT[45])&&(OUT[2]>OUT[46])
&&(OUT[2]>OUT[47])&&(OUT[2]>OUT[48]) ) 
{
MT_VCX=-0.25;
MT_VCY=-0.75;		
	
}
	 
//4
else if((OUT[3]>OUT[0])&&(OUT[3]>OUT[1])&&(OUT[3]>OUT[2])&&(OUT[3]>OUT[4])&&(OUT[3]>OUT[5])&&(OUT[3]>OUT[6])&&(OUT[3]>OUT[7])&&(OUT[3]>OUT[8])&&(OUT[3]>OUT[9])&&(OUT[3]>OUT[10])
&&(OUT[3]>OUT[11])&&(OUT[3]>OUT[12])&&(OUT[3]>OUT[13])&&(OUT[3]>OUT[14])&&(OUT[3]>OUT[15])&&(OUT[3]>OUT[16])&&(OUT[3]>OUT[17])&&(OUT[3]>OUT[18])&&(OUT[3]>OUT[19])
&&(OUT[3]>OUT[20])&&(OUT[3]>OUT[21])&&(OUT[3]>OUT[22])&&(OUT[3]>OUT[23])&&(OUT[3]>OUT[24])&&(OUT[3]>OUT[25])&&(OUT[3]>OUT[26])&&(OUT[3]>OUT[27])&&(OUT[3]>OUT[28])
&&(OUT[3]>OUT[29])&&(OUT[3]>OUT[30])&&(OUT[3]>OUT[31])&&(OUT[3]>OUT[32])&&(OUT[3]>OUT[33])&&(OUT[3]>OUT[34])&&(OUT[3]>OUT[35])&&(OUT[3]>OUT[36])&&(OUT[3]>OUT[37])
&&(OUT[3]>OUT[38])&&(OUT[3]>OUT[39])&&(OUT[3]>OUT[40])&&(OUT[3]>OUT[41])&&(OUT[3]>OUT[42])&&(OUT[3]>OUT[43])&&(OUT[3]>OUT[44])&&(OUT[3]>OUT[45])&&(OUT[3]>OUT[46])
&&(OUT[3]>OUT[47])&&(OUT[3]>OUT[48]) ) 
{
MT_VCX=0;
MT_VCY=-0.75;		
	
}

//5

else if((OUT[4]>OUT[0])&&(OUT[4]>OUT[1])&&(OUT[4]>OUT[2])&&(OUT[4]>OUT[3])&&(OUT[4]>OUT[5])&&(OUT[4]>OUT[6])&&(OUT[4]>OUT[7])&&(OUT[4]>OUT[8])&&(OUT[4]>OUT[9])&&(OUT[4]>OUT[10])
&&(OUT[4]>OUT[11])&&(OUT[4]>OUT[12])&&(OUT[4]>OUT[13])&&(OUT[4]>OUT[14])&&(OUT[4]>OUT[15])&&(OUT[4]>OUT[16])&&(OUT[4]>OUT[17])&&(OUT[4]>OUT[18])&&(OUT[4]>OUT[19])
&&(OUT[4]>OUT[20])&&(OUT[4]>OUT[21])&&(OUT[4]>OUT[22])&&(OUT[4]>OUT[23])&&(OUT[4]>OUT[24])&&(OUT[4]>OUT[25])&&(OUT[4]>OUT[26])&&(OUT[4]>OUT[27])&&(OUT[4]>OUT[28])
&&(OUT[4]>OUT[29])&&(OUT[4]>OUT[30])&&(OUT[4]>OUT[31])&&(OUT[4]>OUT[32])&&(OUT[4]>OUT[33])&&(OUT[4]>OUT[34])&&(OUT[4]>OUT[35])&&(OUT[4]>OUT[36])&&(OUT[4]>OUT[37])
&&(OUT[4]>OUT[38])&&(OUT[4]>OUT[39])&&(OUT[4]>OUT[40])&&(OUT[4]>OUT[41])&&(OUT[4]>OUT[42])&&(OUT[4]>OUT[43])&&(OUT[4]>OUT[44])&&(OUT[4]>OUT[45])&&(OUT[4]>OUT[46])
&&(OUT[4]>OUT[47])&&(OUT[4]>OUT[48]) ) 
{
MT_VCX=0.25;
MT_VCY=-0.75;		
	
}

//6


else if((OUT[5]>OUT[0])&&(OUT[5]>OUT[1])&&(OUT[5]>OUT[2])&&(OUT[5]>OUT[3])&&(OUT[5]>OUT[4])&&(OUT[5]>OUT[6])&&(OUT[5]>OUT[7])&&(OUT[5]>OUT[8])&&(OUT[5]>OUT[9])&&(OUT[5]>OUT[10])
&&(OUT[5]>OUT[11])&&(OUT[5]>OUT[12])&&(OUT[5]>OUT[13])&&(OUT[5]>OUT[14])&&(OUT[5]>OUT[15])&&(OUT[5]>OUT[16])&&(OUT[5]>OUT[17])&&(OUT[5]>OUT[18])&&(OUT[5]>OUT[19])
&&(OUT[5]>OUT[20])&&(OUT[5]>OUT[21])&&(OUT[5]>OUT[22])&&(OUT[5]>OUT[23])&&(OUT[5]>OUT[24])&&(OUT[5]>OUT[25])&&(OUT[5]>OUT[26])&&(OUT[5]>OUT[27])&&(OUT[5]>OUT[28])
&&(OUT[5]>OUT[29])&&(OUT[5]>OUT[30])&&(OUT[5]>OUT[31])&&(OUT[5]>OUT[32])&&(OUT[5]>OUT[33])&&(OUT[5]>OUT[34])&&(OUT[5]>OUT[35])&&(OUT[5]>OUT[36])&&(OUT[5]>OUT[37])
&&(OUT[5]>OUT[38])&&(OUT[5]>OUT[39])&&(OUT[5]>OUT[40])&&(OUT[5]>OUT[41])&&(OUT[5]>OUT[42])&&(OUT[5]>OUT[43])&&(OUT[5]>OUT[44])&&(OUT[5]>OUT[45])&&(OUT[5]>OUT[46])
&&(OUT[5]>OUT[47])&&(OUT[5]>OUT[48]) ) 
{
MT_VCX=0.5;
MT_VCY=-0.75;		
	
}

//7

else if((OUT[6]>OUT[0])&&(OUT[6]>OUT[1])&&(OUT[6]>OUT[2])&&(OUT[6]>OUT[3])&&(OUT[6]>OUT[4])&&(OUT[6]>OUT[5])&&(OUT[6]>OUT[7])&&(OUT[6]>OUT[8])&&(OUT[6]>OUT[9])&&(OUT[6]>OUT[10])
&&(OUT[6]>OUT[11])&&(OUT[6]>OUT[12])&&(OUT[6]>OUT[13])&&(OUT[6]>OUT[14])&&(OUT[6]>OUT[15])&&(OUT[6]>OUT[16])&&(OUT[6]>OUT[17])&&(OUT[6]>OUT[18])&&(OUT[6]>OUT[19])
&&(OUT[6]>OUT[20])&&(OUT[6]>OUT[21])&&(OUT[6]>OUT[22])&&(OUT[6]>OUT[23])&&(OUT[6]>OUT[24])&&(OUT[6]>OUT[25])&&(OUT[6]>OUT[26])&&(OUT[6]>OUT[27])&&(OUT[6]>OUT[28])
&&(OUT[6]>OUT[29])&&(OUT[6]>OUT[30])&&(OUT[6]>OUT[31])&&(OUT[6]>OUT[32])&&(OUT[6]>OUT[33])&&(OUT[6]>OUT[34])&&(OUT[6]>OUT[35])&&(OUT[6]>OUT[36])&&(OUT[6]>OUT[37])
&&(OUT[6]>OUT[38])&&(OUT[6]>OUT[39])&&(OUT[6]>OUT[40])&&(OUT[6]>OUT[41])&&(OUT[6]>OUT[42])&&(OUT[6]>OUT[43])&&(OUT[6]>OUT[44])&&(OUT[6]>OUT[45])&&(OUT[6]>OUT[46])
&&(OUT[6]>OUT[47])&&(OUT[6]>OUT[48]) ) 
{
MT_VCX=0.75;
MT_VCY=-0.75;		
	
}

//8

else if((OUT[7]>OUT[0])&&(OUT[7]>OUT[1])&&(OUT[7]>OUT[2])&&(OUT[7]>OUT[3])&&(OUT[7]>OUT[4])&&(OUT[7]>OUT[5])&&(OUT[7]>OUT[6])&&(OUT[7]>OUT[8])&&(OUT[7]>OUT[9])&&(OUT[7]>OUT[10])
&&(OUT[7]>OUT[11])&&(OUT[7]>OUT[12])&&(OUT[7]>OUT[13])&&(OUT[7]>OUT[14])&&(OUT[7]>OUT[15])&&(OUT[7]>OUT[16])&&(OUT[7]>OUT[17])&&(OUT[7]>OUT[18])&&(OUT[7]>OUT[19])
&&(OUT[7]>OUT[20])&&(OUT[7]>OUT[21])&&(OUT[7]>OUT[22])&&(OUT[7]>OUT[23])&&(OUT[7]>OUT[24])&&(OUT[7]>OUT[25])&&(OUT[7]>OUT[26])&&(OUT[7]>OUT[27])&&(OUT[7]>OUT[28])
&&(OUT[7]>OUT[29])&&(OUT[7]>OUT[30])&&(OUT[7]>OUT[31])&&(OUT[7]>OUT[32])&&(OUT[7]>OUT[33])&&(OUT[7]>OUT[34])&&(OUT[7]>OUT[35])&&(OUT[7]>OUT[36])&&(OUT[7]>OUT[37])
&&(OUT[7]>OUT[38])&&(OUT[7]>OUT[39])&&(OUT[7]>OUT[40])&&(OUT[7]>OUT[41])&&(OUT[7]>OUT[42])&&(OUT[7]>OUT[43])&&(OUT[7]>OUT[44])&&(OUT[7]>OUT[45])&&(OUT[7]>OUT[46])
&&(OUT[7]>OUT[47])&&(OUT[7]>OUT[48]) ) 
{
MT_VCX=-0.75;
MT_VCY=-0.5;		
	
}

//9


else if((OUT[8]>OUT[0])&&(OUT[8]>OUT[1])&&(OUT[8]>OUT[2])&&(OUT[8]>OUT[3])&&(OUT[8]>OUT[4])&&(OUT[8]>OUT[5])&&(OUT[8]>OUT[6])&&(OUT[8]>OUT[7])&&(OUT[8]>OUT[9])&&(OUT[8]>OUT[10])
&&(OUT[8]>OUT[11])&&(OUT[8]>OUT[12])&&(OUT[8]>OUT[13])&&(OUT[8]>OUT[14])&&(OUT[8]>OUT[15])&&(OUT[8]>OUT[16])&&(OUT[8]>OUT[17])&&(OUT[8]>OUT[18])&&(OUT[8]>OUT[19])
&&(OUT[8]>OUT[20])&&(OUT[8]>OUT[21])&&(OUT[8]>OUT[22])&&(OUT[8]>OUT[23])&&(OUT[8]>OUT[24])&&(OUT[8]>OUT[25])&&(OUT[8]>OUT[26])&&(OUT[8]>OUT[27])&&(OUT[8]>OUT[28])
&&(OUT[8]>OUT[29])&&(OUT[8]>OUT[30])&&(OUT[8]>OUT[31])&&(OUT[8]>OUT[32])&&(OUT[8]>OUT[33])&&(OUT[8]>OUT[34])&&(OUT[8]>OUT[35])&&(OUT[8]>OUT[36])&&(OUT[8]>OUT[37])
&&(OUT[8]>OUT[38])&&(OUT[8]>OUT[39])&&(OUT[8]>OUT[40])&&(OUT[8]>OUT[41])&&(OUT[8]>OUT[42])&&(OUT[8]>OUT[43])&&(OUT[8]>OUT[44])&&(OUT[8]>OUT[45])&&(OUT[8]>OUT[46])
&&(OUT[8]>OUT[47])&&(OUT[8]>OUT[48]) ) 
{
MT_VCX=-0.5;
MT_VCY=-0.5;	
	
}


//10


else if((OUT[9]>OUT[0])&&(OUT[9]>OUT[1])&&(OUT[9]>OUT[2])&&(OUT[9]>OUT[3])&&(OUT[9]>OUT[4])&&(OUT[9]>OUT[5])&&(OUT[9]>OUT[6])&&(OUT[9]>OUT[7])&&(OUT[9]>OUT[8])&&(OUT[9]>OUT[10])
&&(OUT[9]>OUT[11])&&(OUT[9]>OUT[12])&&(OUT[9]>OUT[13])&&(OUT[9]>OUT[14])&&(OUT[9]>OUT[15])&&(OUT[9]>OUT[16])&&(OUT[9]>OUT[17])&&(OUT[9]>OUT[18])&&(OUT[9]>OUT[19])
&&(OUT[9]>OUT[20])&&(OUT[9]>OUT[21])&&(OUT[9]>OUT[22])&&(OUT[9]>OUT[23])&&(OUT[9]>OUT[24])&&(OUT[9]>OUT[25])&&(OUT[9]>OUT[26])&&(OUT[9]>OUT[27])&&(OUT[9]>OUT[28])
&&(OUT[9]>OUT[29])&&(OUT[9]>OUT[30])&&(OUT[9]>OUT[31])&&(OUT[9]>OUT[32])&&(OUT[9]>OUT[33])&&(OUT[9]>OUT[34])&&(OUT[9]>OUT[35])&&(OUT[9]>OUT[36])&&(OUT[9]>OUT[37])
&&(OUT[9]>OUT[38])&&(OUT[9]>OUT[39])&&(OUT[9]>OUT[40])&&(OUT[9]>OUT[41])&&(OUT[9]>OUT[42])&&(OUT[9]>OUT[43])&&(OUT[9]>OUT[44])&&(OUT[9]>OUT[45])&&(OUT[9]>OUT[46])
&&(OUT[9]>OUT[47])&&(OUT[9]>OUT[48]) ) 
{
MT_VCX=-0.25;
MT_VCY=-0.5;	
	
}

//11

else if((OUT[10]>OUT[0])&&(OUT[10]>OUT[1])&&(OUT[10]>OUT[2])&&(OUT[10]>OUT[3])&&(OUT[10]>OUT[4])&&(OUT[10]>OUT[5])&&(OUT[10]>OUT[6])&&(OUT[10]>OUT[7])&&(OUT[10]>OUT[8])&&(OUT[10]>OUT[9])
&&(OUT[10]>OUT[11])&&(OUT[10]>OUT[12])&&(OUT[10]>OUT[13])&&(OUT[10]>OUT[14])&&(OUT[10]>OUT[15])&&(OUT[10]>OUT[16])&&(OUT[10]>OUT[17])&&(OUT[10]>OUT[18])&&(OUT[10]>OUT[19])
&&(OUT[10]>OUT[20])&&(OUT[10]>OUT[21])&&(OUT[10]>OUT[22])&&(OUT[10]>OUT[23])&&(OUT[10]>OUT[24])&&(OUT[10]>OUT[25])&&(OUT[10]>OUT[26])&&(OUT[10]>OUT[27])&&(OUT[10]>OUT[28])
&&(OUT[10]>OUT[29])&&(OUT[10]>OUT[30])&&(OUT[10]>OUT[31])&&(OUT[10]>OUT[32])&&(OUT[10]>OUT[33])&&(OUT[10]>OUT[34])&&(OUT[10]>OUT[35])&&(OUT[10]>OUT[36])&&(OUT[10]>OUT[37])
&&(OUT[10]>OUT[38])&&(OUT[10]>OUT[39])&&(OUT[10]>OUT[40])&&(OUT[10]>OUT[41])&&(OUT[10]>OUT[42])&&(OUT[10]>OUT[43])&&(OUT[10]>OUT[44])&&(OUT[10]>OUT[45])&&(OUT[10]>OUT[46])
&&(OUT[10]>OUT[47])&&(OUT[10]>OUT[48]) ) 
{
MT_VCX=0;
MT_VCY=-0.5;	
		
	
}

//12

else if((OUT[11]>OUT[0])&&(OUT[11]>OUT[1])&&(OUT[11]>OUT[2])&&(OUT[11]>OUT[3])&&(OUT[11]>OUT[4])&&(OUT[11]>OUT[5])&&(OUT[11]>OUT[6])&&(OUT[11]>OUT[7])&&(OUT[11]>OUT[8])&&(OUT[11]>OUT[9])
&&(OUT[11]>OUT[10])&&(OUT[11]>OUT[12])&&(OUT[11]>OUT[13])&&(OUT[11]>OUT[14])&&(OUT[11]>OUT[15])&&(OUT[11]>OUT[16])&&(OUT[11]>OUT[17])&&(OUT[11]>OUT[18])&&(OUT[11]>OUT[19])
&&(OUT[11]>OUT[20])&&(OUT[11]>OUT[21])&&(OUT[11]>OUT[22])&&(OUT[11]>OUT[23])&&(OUT[11]>OUT[24])&&(OUT[11]>OUT[25])&&(OUT[11]>OUT[26])&&(OUT[11]>OUT[27])&&(OUT[11]>OUT[28])
&&(OUT[11]>OUT[29])&&(OUT[11]>OUT[30])&&(OUT[11]>OUT[31])&&(OUT[11]>OUT[32])&&(OUT[11]>OUT[33])&&(OUT[11]>OUT[34])&&(OUT[11]>OUT[35])&&(OUT[11]>OUT[36])&&(OUT[11]>OUT[37])
&&(OUT[11]>OUT[38])&&(OUT[11]>OUT[39])&&(OUT[11]>OUT[40])&&(OUT[11]>OUT[41])&&(OUT[11]>OUT[42])&&(OUT[11]>OUT[43])&&(OUT[11]>OUT[44])&&(OUT[11]>OUT[45])&&(OUT[11]>OUT[46])
&&(OUT[11]>OUT[47])&&(OUT[11]>OUT[48]) ) 
{
MT_VCX=0.25;
MT_VCY=-0.5;	
	
	
}


//13

else if((OUT[12]>OUT[0])&&(OUT[12]>OUT[1])&&(OUT[12]>OUT[2])&&(OUT[12]>OUT[3])&&(OUT[12]>OUT[4])&&(OUT[12]>OUT[5])&&(OUT[12]>OUT[6])&&(OUT[12]>OUT[7])&&(OUT[12]>OUT[8])&&(OUT[12]>OUT[9])
&&(OUT[12]>OUT[10])&&(OUT[12]>OUT[11])&&(OUT[12]>OUT[13])&&(OUT[12]>OUT[14])&&(OUT[12]>OUT[15])&&(OUT[12]>OUT[16])&&(OUT[12]>OUT[17])&&(OUT[12]>OUT[18])&&(OUT[12]>OUT[19])
&&(OUT[12]>OUT[20])&&(OUT[12]>OUT[21])&&(OUT[12]>OUT[22])&&(OUT[12]>OUT[23])&&(OUT[12]>OUT[24])&&(OUT[12]>OUT[25])&&(OUT[12]>OUT[26])&&(OUT[12]>OUT[27])&&(OUT[12]>OUT[28])
&&(OUT[12]>OUT[29])&&(OUT[12]>OUT[30])&&(OUT[12]>OUT[31])&&(OUT[12]>OUT[32])&&(OUT[12]>OUT[33])&&(OUT[12]>OUT[34])&&(OUT[12]>OUT[35])&&(OUT[12]>OUT[36])&&(OUT[12]>OUT[37])
&&(OUT[12]>OUT[38])&&(OUT[12]>OUT[39])&&(OUT[12]>OUT[40])&&(OUT[12]>OUT[41])&&(OUT[12]>OUT[42])&&(OUT[12]>OUT[43])&&(OUT[12]>OUT[44])&&(OUT[12]>OUT[45])&&(OUT[12]>OUT[46])
&&(OUT[12]>OUT[47])&&(OUT[12]>OUT[48]) ) 
{
MT_VCX=0.5;
MT_VCY=-0.5;		
	
}

//14

else if((OUT[13]>OUT[0])&&(OUT[13]>OUT[1])&&(OUT[13]>OUT[2])&&(OUT[13]>OUT[3])&&(OUT[13]>OUT[4])&&(OUT[13]>OUT[5])&&(OUT[13]>OUT[6])&&(OUT[13]>OUT[7])&&(OUT[13]>OUT[8])&&(OUT[13]>OUT[9])
&&(OUT[13]>OUT[10])&&(OUT[13]>OUT[11])&&(OUT[13]>OUT[12])&&(OUT[13]>OUT[14])&&(OUT[13]>OUT[15])&&(OUT[13]>OUT[16])&&(OUT[13]>OUT[17])&&(OUT[13]>OUT[18])&&(OUT[13]>OUT[19])
&&(OUT[13]>OUT[20])&&(OUT[13]>OUT[21])&&(OUT[13]>OUT[22])&&(OUT[13]>OUT[23])&&(OUT[13]>OUT[24])&&(OUT[13]>OUT[25])&&(OUT[13]>OUT[26])&&(OUT[13]>OUT[27])&&(OUT[13]>OUT[28])
&&(OUT[13]>OUT[29])&&(OUT[13]>OUT[30])&&(OUT[13]>OUT[31])&&(OUT[13]>OUT[32])&&(OUT[13]>OUT[33])&&(OUT[13]>OUT[34])&&(OUT[13]>OUT[35])&&(OUT[13]>OUT[36])&&(OUT[13]>OUT[37])
&&(OUT[13]>OUT[38])&&(OUT[13]>OUT[39])&&(OUT[13]>OUT[40])&&(OUT[13]>OUT[41])&&(OUT[13]>OUT[42])&&(OUT[13]>OUT[43])&&(OUT[13]>OUT[44])&&(OUT[13]>OUT[45])&&(OUT[13]>OUT[46])
&&(OUT[13]>OUT[47])&&(OUT[13]>OUT[48]) ) 
{
MT_VCX=0.75;
MT_VCY=-0.5;		
	
}

//15

else if((OUT[14]>OUT[0])&&(OUT[14]>OUT[1])&&(OUT[14]>OUT[2])&&(OUT[14]>OUT[3])&&(OUT[14]>OUT[4])&&(OUT[14]>OUT[5])&&(OUT[14]>OUT[6])&&(OUT[14]>OUT[7])&&(OUT[14]>OUT[8])&&(OUT[14]>OUT[9])
&&(OUT[14]>OUT[10])&&(OUT[14]>OUT[11])&&(OUT[14]>OUT[12])&&(OUT[14]>OUT[13])&&(OUT[14]>OUT[15])&&(OUT[14]>OUT[16])&&(OUT[14]>OUT[17])&&(OUT[14]>OUT[18])&&(OUT[14]>OUT[19])
&&(OUT[14]>OUT[20])&&(OUT[14]>OUT[21])&&(OUT[14]>OUT[22])&&(OUT[14]>OUT[23])&&(OUT[14]>OUT[24])&&(OUT[14]>OUT[25])&&(OUT[14]>OUT[26])&&(OUT[14]>OUT[27])&&(OUT[14]>OUT[28])
&&(OUT[14]>OUT[29])&&(OUT[14]>OUT[30])&&(OUT[14]>OUT[31])&&(OUT[14]>OUT[32])&&(OUT[14]>OUT[33])&&(OUT[14]>OUT[34])&&(OUT[14]>OUT[35])&&(OUT[14]>OUT[36])&&(OUT[14]>OUT[37])
&&(OUT[14]>OUT[38])&&(OUT[14]>OUT[39])&&(OUT[14]>OUT[40])&&(OUT[14]>OUT[41])&&(OUT[14]>OUT[42])&&(OUT[14]>OUT[43])&&(OUT[14]>OUT[44])&&(OUT[14]>OUT[45])&&(OUT[14]>OUT[46])
&&(OUT[14]>OUT[47])&&(OUT[14]>OUT[48]) ) 
{
MT_VCX=-0.75;
MT_VCY=-0.25;			
	
}

//16

else if((OUT[15]>OUT[0])&&(OUT[15]>OUT[1])&&(OUT[15]>OUT[2])&&(OUT[15]>OUT[3])&&(OUT[15]>OUT[4])&&(OUT[15]>OUT[5])&&(OUT[15]>OUT[6])&&(OUT[15]>OUT[7])&&(OUT[15]>OUT[8])&&(OUT[15]>OUT[9])
&&(OUT[15]>OUT[10])&&(OUT[15]>OUT[11])&&(OUT[15]>OUT[12])&&(OUT[15]>OUT[13])&&(OUT[15]>OUT[14])&&(OUT[15]>OUT[16])&&(OUT[15]>OUT[17])&&(OUT[15]>OUT[18])&&(OUT[15]>OUT[19])
&&(OUT[15]>OUT[20])&&(OUT[15]>OUT[21])&&(OUT[15]>OUT[22])&&(OUT[15]>OUT[23])&&(OUT[15]>OUT[24])&&(OUT[15]>OUT[25])&&(OUT[15]>OUT[26])&&(OUT[15]>OUT[27])&&(OUT[15]>OUT[28])
&&(OUT[15]>OUT[29])&&(OUT[15]>OUT[30])&&(OUT[15]>OUT[31])&&(OUT[15]>OUT[32])&&(OUT[15]>OUT[33])&&(OUT[15]>OUT[34])&&(OUT[15]>OUT[35])&&(OUT[15]>OUT[36])&&(OUT[15]>OUT[37])
&&(OUT[15]>OUT[38])&&(OUT[15]>OUT[39])&&(OUT[15]>OUT[40])&&(OUT[15]>OUT[41])&&(OUT[15]>OUT[42])&&(OUT[15]>OUT[43])&&(OUT[15]>OUT[44])&&(OUT[15]>OUT[45])&&(OUT[15]>OUT[46])
&&(OUT[15]>OUT[47])&&(OUT[15]>OUT[48]) ) 
{
MT_VCX=-0.5;
MT_VCY=-0.25;		
	
}

//17

else if((OUT[16]>OUT[0])&&(OUT[16]>OUT[1])&&(OUT[16]>OUT[2])&&(OUT[16]>OUT[3])&&(OUT[16]>OUT[4])&&(OUT[16]>OUT[5])&&(OUT[16]>OUT[6])&&(OUT[16]>OUT[7])&&(OUT[16]>OUT[8])&&(OUT[16]>OUT[9])
&&(OUT[16]>OUT[10])&&(OUT[16]>OUT[11])&&(OUT[16]>OUT[12])&&(OUT[16]>OUT[13])&&(OUT[16]>OUT[14])&&(OUT[16]>OUT[15])&&(OUT[16]>OUT[17])&&(OUT[16]>OUT[18])&&(OUT[16]>OUT[19])
&&(OUT[16]>OUT[20])&&(OUT[16]>OUT[21])&&(OUT[16]>OUT[22])&&(OUT[16]>OUT[23])&&(OUT[16]>OUT[24])&&(OUT[16]>OUT[25])&&(OUT[16]>OUT[26])&&(OUT[16]>OUT[27])&&(OUT[16]>OUT[28])
&&(OUT[16]>OUT[29])&&(OUT[16]>OUT[30])&&(OUT[16]>OUT[31])&&(OUT[16]>OUT[32])&&(OUT[16]>OUT[33])&&(OUT[16]>OUT[34])&&(OUT[16]>OUT[35])&&(OUT[16]>OUT[36])&&(OUT[16]>OUT[37])
&&(OUT[16]>OUT[38])&&(OUT[16]>OUT[39])&&(OUT[16]>OUT[40])&&(OUT[16]>OUT[41])&&(OUT[16]>OUT[42])&&(OUT[16]>OUT[43])&&(OUT[16]>OUT[44])&&(OUT[16]>OUT[45])&&(OUT[16]>OUT[46])
&&(OUT[16]>OUT[47])&&(OUT[16]>OUT[48]) ) 
{
MT_VCX=-0.25;
MT_VCY=-0.25;		
	
}

//18
	
else if((OUT[17]>OUT[0])&&(OUT[17]>OUT[1])&&(OUT[17]>OUT[2])&&(OUT[17]>OUT[3])&&(OUT[17]>OUT[4])&&(OUT[17]>OUT[5])&&(OUT[17]>OUT[6])&&(OUT[17]>OUT[7])&&(OUT[17]>OUT[8])&&(OUT[17]>OUT[9])
&&(OUT[17]>OUT[10])&&(OUT[17]>OUT[11])&&(OUT[17]>OUT[12])&&(OUT[17]>OUT[13])&&(OUT[17]>OUT[14])&&(OUT[17]>OUT[15])&&(OUT[17]>OUT[16])&&(OUT[17]>OUT[18])&&(OUT[17]>OUT[19])
&&(OUT[17]>OUT[20])&&(OUT[17]>OUT[21])&&(OUT[17]>OUT[22])&&(OUT[17]>OUT[23])&&(OUT[17]>OUT[24])&&(OUT[17]>OUT[25])&&(OUT[17]>OUT[26])&&(OUT[17]>OUT[27])&&(OUT[17]>OUT[28])
&&(OUT[17]>OUT[29])&&(OUT[17]>OUT[30])&&(OUT[17]>OUT[31])&&(OUT[17]>OUT[32])&&(OUT[17]>OUT[33])&&(OUT[17]>OUT[34])&&(OUT[17]>OUT[35])&&(OUT[17]>OUT[36])&&(OUT[17]>OUT[37])
&&(OUT[17]>OUT[38])&&(OUT[17]>OUT[39])&&(OUT[17]>OUT[40])&&(OUT[17]>OUT[41])&&(OUT[17]>OUT[42])&&(OUT[17]>OUT[43])&&(OUT[17]>OUT[44])&&(OUT[17]>OUT[45])&&(OUT[17]>OUT[46])
&&(OUT[17]>OUT[47])&&(OUT[17]>OUT[48]) ) 
{
MT_VCX=0;
MT_VCY=-0.25;		
	
}

//19

else if((OUT[18]>OUT[0])&&(OUT[18]>OUT[1])&&(OUT[18]>OUT[2])&&(OUT[18]>OUT[3])&&(OUT[18]>OUT[4])&&(OUT[18]>OUT[5])&&(OUT[18]>OUT[6])&&(OUT[18]>OUT[7])&&(OUT[18]>OUT[8])&&(OUT[18]>OUT[9])
&&(OUT[18]>OUT[10])&&(OUT[18]>OUT[11])&&(OUT[18]>OUT[12])&&(OUT[18]>OUT[13])&&(OUT[18]>OUT[14])&&(OUT[18]>OUT[15])&&(OUT[18]>OUT[16])&&(OUT[18]>OUT[17])&&(OUT[18]>OUT[19])
&&(OUT[18]>OUT[20])&&(OUT[18]>OUT[21])&&(OUT[18]>OUT[22])&&(OUT[18]>OUT[23])&&(OUT[18]>OUT[24])&&(OUT[18]>OUT[25])&&(OUT[18]>OUT[26])&&(OUT[18]>OUT[27])&&(OUT[18]>OUT[28])
&&(OUT[18]>OUT[29])&&(OUT[18]>OUT[30])&&(OUT[18]>OUT[31])&&(OUT[18]>OUT[32])&&(OUT[18]>OUT[33])&&(OUT[18]>OUT[34])&&(OUT[18]>OUT[35])&&(OUT[18]>OUT[36])&&(OUT[18]>OUT[37])
&&(OUT[18]>OUT[38])&&(OUT[18]>OUT[39])&&(OUT[18]>OUT[40])&&(OUT[18]>OUT[41])&&(OUT[18]>OUT[42])&&(OUT[18]>OUT[43])&&(OUT[18]>OUT[44])&&(OUT[18]>OUT[45])&&(OUT[18]>OUT[46])
&&(OUT[18]>OUT[47])&&(OUT[18]>OUT[48]) ) 
{
MT_VCX=0.25;
MT_VCY=-0.25;		
	
}

//20

else if((OUT[19]>OUT[0])&&(OUT[19]>OUT[1])&&(OUT[19]>OUT[2])&&(OUT[19]>OUT[3])&&(OUT[19]>OUT[4])&&(OUT[19]>OUT[5])&&(OUT[19]>OUT[6])&&(OUT[19]>OUT[7])&&(OUT[19]>OUT[8])&&(OUT[19]>OUT[9])
&&(OUT[19]>OUT[10])&&(OUT[19]>OUT[11])&&(OUT[19]>OUT[12])&&(OUT[19]>OUT[13])&&(OUT[19]>OUT[14])&&(OUT[19]>OUT[15])&&(OUT[19]>OUT[16])&&(OUT[19]>OUT[17])&&(OUT[19]>OUT[18])
&&(OUT[19]>OUT[20])&&(OUT[19]>OUT[21])&&(OUT[19]>OUT[22])&&(OUT[19]>OUT[23])&&(OUT[19]>OUT[24])&&(OUT[19]>OUT[25])&&(OUT[19]>OUT[26])&&(OUT[19]>OUT[27])&&(OUT[19]>OUT[28])
&&(OUT[19]>OUT[29])&&(OUT[19]>OUT[30])&&(OUT[19]>OUT[31])&&(OUT[19]>OUT[32])&&(OUT[19]>OUT[33])&&(OUT[19]>OUT[34])&&(OUT[19]>OUT[35])&&(OUT[19]>OUT[36])&&(OUT[19]>OUT[37])
&&(OUT[19]>OUT[38])&&(OUT[19]>OUT[39])&&(OUT[19]>OUT[40])&&(OUT[19]>OUT[41])&&(OUT[19]>OUT[42])&&(OUT[19]>OUT[43])&&(OUT[19]>OUT[44])&&(OUT[19]>OUT[45])&&(OUT[19]>OUT[46])
&&(OUT[19]>OUT[47])&&(OUT[19]>OUT[48]) ) 
{
MT_VCX=0.5;
MT_VCY=-0.25;		
	
}

//21

else if((OUT[20]>OUT[0])&&(OUT[20]>OUT[1])&&(OUT[20]>OUT[2])&&(OUT[20]>OUT[3])&&(OUT[20]>OUT[4])&&(OUT[20]>OUT[5])&&(OUT[20]>OUT[6])&&(OUT[20]>OUT[7])&&(OUT[20]>OUT[8])&&(OUT[20]>OUT[9])
&&(OUT[20]>OUT[10])&&(OUT[20]>OUT[11])&&(OUT[20]>OUT[12])&&(OUT[20]>OUT[13])&&(OUT[20]>OUT[14])&&(OUT[20]>OUT[15])&&(OUT[20]>OUT[16])&&(OUT[20]>OUT[17])&&(OUT[20]>OUT[18])
&&(OUT[20]>OUT[19])&&(OUT[20]>OUT[21])&&(OUT[20]>OUT[22])&&(OUT[20]>OUT[23])&&(OUT[20]>OUT[24])&&(OUT[20]>OUT[25])&&(OUT[20]>OUT[26])&&(OUT[20]>OUT[27])&&(OUT[20]>OUT[28])
&&(OUT[20]>OUT[29])&&(OUT[20]>OUT[30])&&(OUT[20]>OUT[31])&&(OUT[20]>OUT[32])&&(OUT[20]>OUT[33])&&(OUT[20]>OUT[34])&&(OUT[20]>OUT[35])&&(OUT[20]>OUT[36])&&(OUT[20]>OUT[37])
&&(OUT[20]>OUT[38])&&(OUT[20]>OUT[39])&&(OUT[20]>OUT[40])&&(OUT[20]>OUT[41])&&(OUT[20]>OUT[42])&&(OUT[20]>OUT[43])&&(OUT[20]>OUT[44])&&(OUT[20]>OUT[45])&&(OUT[20]>OUT[46])
&&(OUT[20]>OUT[47])&&(OUT[20]>OUT[48]) ) 
{
MT_VCX=0.75;
MT_VCY=-0.25;		
	
}

//22

else if((OUT[21]>OUT[0])&&(OUT[21]>OUT[1])&&(OUT[21]>OUT[2])&&(OUT[21]>OUT[3])&&(OUT[21]>OUT[4])&&(OUT[21]>OUT[5])&&(OUT[21]>OUT[6])&&(OUT[21]>OUT[7])&&(OUT[21]>OUT[8])&&(OUT[21]>OUT[9])
&&(OUT[21]>OUT[10])&&(OUT[21]>OUT[11])&&(OUT[21]>OUT[12])&&(OUT[21]>OUT[13])&&(OUT[21]>OUT[14])&&(OUT[21]>OUT[15])&&(OUT[21]>OUT[16])&&(OUT[21]>OUT[17])&&(OUT[21]>OUT[18])
&&(OUT[21]>OUT[19])&&(OUT[21]>OUT[20])&&(OUT[21]>OUT[22])&&(OUT[21]>OUT[23])&&(OUT[21]>OUT[24])&&(OUT[21]>OUT[25])&&(OUT[21]>OUT[26])&&(OUT[21]>OUT[27])&&(OUT[21]>OUT[28])
&&(OUT[21]>OUT[29])&&(OUT[21]>OUT[30])&&(OUT[21]>OUT[31])&&(OUT[21]>OUT[32])&&(OUT[21]>OUT[33])&&(OUT[21]>OUT[34])&&(OUT[21]>OUT[35])&&(OUT[21]>OUT[36])&&(OUT[21]>OUT[37])
&&(OUT[21]>OUT[38])&&(OUT[21]>OUT[39])&&(OUT[21]>OUT[40])&&(OUT[21]>OUT[41])&&(OUT[21]>OUT[42])&&(OUT[21]>OUT[43])&&(OUT[21]>OUT[44])&&(OUT[21]>OUT[45])&&(OUT[21]>OUT[46])
&&(OUT[21]>OUT[47])&&(OUT[21]>OUT[48]) ) 
{
MT_VCX=-0.75;
MT_VCY=0;	
	
}

//23

else if((OUT[22]>OUT[0])&&(OUT[22]>OUT[1])&&(OUT[22]>OUT[2])&&(OUT[22]>OUT[3])&&(OUT[22]>OUT[4])&&(OUT[22]>OUT[5])&&(OUT[22]>OUT[6])&&(OUT[22]>OUT[7])&&(OUT[22]>OUT[8])&&(OUT[22]>OUT[9])
&&(OUT[22]>OUT[10])&&(OUT[22]>OUT[11])&&(OUT[22]>OUT[12])&&(OUT[22]>OUT[13])&&(OUT[22]>OUT[14])&&(OUT[22]>OUT[15])&&(OUT[22]>OUT[16])&&(OUT[22]>OUT[17])&&(OUT[22]>OUT[18])
&&(OUT[22]>OUT[19])&&(OUT[22]>OUT[20])&&(OUT[22]>OUT[21])&&(OUT[22]>OUT[23])&&(OUT[22]>OUT[24])&&(OUT[22]>OUT[25])&&(OUT[22]>OUT[26])&&(OUT[22]>OUT[27])&&(OUT[22]>OUT[28])
&&(OUT[22]>OUT[29])&&(OUT[22]>OUT[30])&&(OUT[22]>OUT[31])&&(OUT[22]>OUT[32])&&(OUT[22]>OUT[33])&&(OUT[22]>OUT[34])&&(OUT[22]>OUT[35])&&(OUT[22]>OUT[36])&&(OUT[22]>OUT[37])
&&(OUT[22]>OUT[38])&&(OUT[22]>OUT[39])&&(OUT[22]>OUT[40])&&(OUT[22]>OUT[41])&&(OUT[22]>OUT[42])&&(OUT[22]>OUT[43])&&(OUT[22]>OUT[44])&&(OUT[22]>OUT[45])&&(OUT[22]>OUT[46])
&&(OUT[22]>OUT[47])&&(OUT[22]>OUT[48]) ) 
{
MT_VCX=-0.5;
MT_VCY=0;		
	
}

//24

else if((OUT[23]>OUT[0])&&(OUT[23]>OUT[1])&&(OUT[23]>OUT[2])&&(OUT[23]>OUT[3])&&(OUT[23]>OUT[4])&&(OUT[23]>OUT[5])&&(OUT[23]>OUT[6])&&(OUT[23]>OUT[7])&&(OUT[23]>OUT[8])&&(OUT[23]>OUT[9])
&&(OUT[23]>OUT[10])&&(OUT[23]>OUT[11])&&(OUT[23]>OUT[12])&&(OUT[23]>OUT[13])&&(OUT[23]>OUT[14])&&(OUT[23]>OUT[15])&&(OUT[23]>OUT[16])&&(OUT[23]>OUT[17])&&(OUT[23]>OUT[18])
&&(OUT[23]>OUT[19])&&(OUT[23]>OUT[20])&&(OUT[23]>OUT[21])&&(OUT[23]>OUT[22])&&(OUT[23]>OUT[24])&&(OUT[23]>OUT[25])&&(OUT[23]>OUT[26])&&(OUT[23]>OUT[27])&&(OUT[23]>OUT[28])
&&(OUT[23]>OUT[29])&&(OUT[23]>OUT[30])&&(OUT[23]>OUT[31])&&(OUT[23]>OUT[32])&&(OUT[23]>OUT[33])&&(OUT[23]>OUT[34])&&(OUT[23]>OUT[35])&&(OUT[23]>OUT[36])&&(OUT[23]>OUT[37])
&&(OUT[23]>OUT[38])&&(OUT[23]>OUT[39])&&(OUT[23]>OUT[40])&&(OUT[23]>OUT[41])&&(OUT[23]>OUT[42])&&(OUT[23]>OUT[43])&&(OUT[23]>OUT[44])&&(OUT[23]>OUT[45])&&(OUT[23]>OUT[46])
&&(OUT[23]>OUT[47])&&(OUT[23]>OUT[48]) ) 
{
MT_VCX=-0.25;
MT_VCY=0;		
	
}

//25

else if((OUT[24]>OUT[0])&&(OUT[24]>OUT[1])&&(OUT[24]>OUT[2])&&(OUT[24]>OUT[3])&&(OUT[24]>OUT[4])&&(OUT[24]>OUT[5])&&(OUT[24]>OUT[6])&&(OUT[24]>OUT[7])&&(OUT[24]>OUT[8])&&(OUT[24]>OUT[9])
&&(OUT[24]>OUT[10])&&(OUT[24]>OUT[11])&&(OUT[24]>OUT[12])&&(OUT[24]>OUT[13])&&(OUT[24]>OUT[14])&&(OUT[24]>OUT[15])&&(OUT[24]>OUT[16])&&(OUT[24]>OUT[17])&&(OUT[24]>OUT[18])
&&(OUT[24]>OUT[19])&&(OUT[24]>OUT[20])&&(OUT[24]>OUT[21])&&(OUT[24]>OUT[22])&&(OUT[24]>OUT[23])&&(OUT[24]>OUT[25])&&(OUT[24]>OUT[26])&&(OUT[24]>OUT[27])&&(OUT[24]>OUT[28])
&&(OUT[24]>OUT[29])&&(OUT[24]>OUT[30])&&(OUT[24]>OUT[31])&&(OUT[24]>OUT[32])&&(OUT[24]>OUT[33])&&(OUT[24]>OUT[34])&&(OUT[24]>OUT[35])&&(OUT[24]>OUT[36])&&(OUT[24]>OUT[37])
&&(OUT[24]>OUT[38])&&(OUT[24]>OUT[39])&&(OUT[24]>OUT[40])&&(OUT[24]>OUT[41])&&(OUT[24]>OUT[42])&&(OUT[24]>OUT[43])&&(OUT[24]>OUT[44])&&(OUT[24]>OUT[45])&&(OUT[24]>OUT[46])
&&(OUT[24]>OUT[47])&&(OUT[24]>OUT[48]) ) 
{
MT_VCX=0;
MT_VCY=0;	
	
}

//26

else if((OUT[25]>OUT[0])&&(OUT[25]>OUT[1])&&(OUT[25]>OUT[2])&&(OUT[25]>OUT[3])&&(OUT[25]>OUT[4])&&(OUT[25]>OUT[5])&&(OUT[25]>OUT[6])&&(OUT[25]>OUT[7])&&(OUT[25]>OUT[8])&&(OUT[25]>OUT[9])
&&(OUT[25]>OUT[10])&&(OUT[25]>OUT[11])&&(OUT[25]>OUT[12])&&(OUT[25]>OUT[13])&&(OUT[25]>OUT[14])&&(OUT[25]>OUT[15])&&(OUT[25]>OUT[16])&&(OUT[25]>OUT[17])&&(OUT[25]>OUT[18])
&&(OUT[25]>OUT[19])&&(OUT[25]>OUT[20])&&(OUT[25]>OUT[21])&&(OUT[25]>OUT[22])&&(OUT[25]>OUT[23])&&(OUT[25]>OUT[24])&&(OUT[25]>OUT[26])&&(OUT[25]>OUT[27])&&(OUT[25]>OUT[28])
&&(OUT[25]>OUT[29])&&(OUT[25]>OUT[30])&&(OUT[25]>OUT[31])&&(OUT[25]>OUT[32])&&(OUT[25]>OUT[33])&&(OUT[25]>OUT[34])&&(OUT[25]>OUT[35])&&(OUT[25]>OUT[36])&&(OUT[25]>OUT[37])
&&(OUT[25]>OUT[38])&&(OUT[25]>OUT[39])&&(OUT[25]>OUT[40])&&(OUT[25]>OUT[41])&&(OUT[25]>OUT[42])&&(OUT[25]>OUT[43])&&(OUT[25]>OUT[44])&&(OUT[25]>OUT[45])&&(OUT[25]>OUT[46])
&&(OUT[25]>OUT[47])&&(OUT[25]>OUT[48]) ) 
{
MT_VCX=0.25;
MT_VCY=0;		
	
}

//27

else if((OUT[26]>OUT[0])&&(OUT[26]>OUT[1])&&(OUT[26]>OUT[2])&&(OUT[26]>OUT[3])&&(OUT[26]>OUT[4])&&(OUT[26]>OUT[5])&&(OUT[26]>OUT[6])&&(OUT[26]>OUT[7])&&(OUT[26]>OUT[8])&&(OUT[26]>OUT[9])
&&(OUT[26]>OUT[10])&&(OUT[26]>OUT[11])&&(OUT[26]>OUT[12])&&(OUT[26]>OUT[13])&&(OUT[26]>OUT[14])&&(OUT[26]>OUT[15])&&(OUT[26]>OUT[16])&&(OUT[26]>OUT[17])&&(OUT[26]>OUT[18])
&&(OUT[26]>OUT[19])&&(OUT[26]>OUT[20])&&(OUT[26]>OUT[21])&&(OUT[26]>OUT[22])&&(OUT[26]>OUT[23])&&(OUT[26]>OUT[24])&&(OUT[26]>OUT[25])&&(OUT[26]>OUT[27])&&(OUT[26]>OUT[28])
&&(OUT[26]>OUT[29])&&(OUT[26]>OUT[30])&&(OUT[26]>OUT[31])&&(OUT[26]>OUT[32])&&(OUT[26]>OUT[33])&&(OUT[26]>OUT[34])&&(OUT[26]>OUT[35])&&(OUT[26]>OUT[36])&&(OUT[26]>OUT[37])
&&(OUT[26]>OUT[38])&&(OUT[26]>OUT[39])&&(OUT[26]>OUT[40])&&(OUT[26]>OUT[41])&&(OUT[26]>OUT[42])&&(OUT[26]>OUT[43])&&(OUT[26]>OUT[44])&&(OUT[26]>OUT[45])&&(OUT[26]>OUT[46])
&&(OUT[26]>OUT[47])&&(OUT[26]>OUT[48]) ) 
{
MT_VCX=0.5;
MT_VCY=0;		
	
}


//28

else if((OUT[27]>OUT[0])&&(OUT[27]>OUT[1])&&(OUT[27]>OUT[2])&&(OUT[27]>OUT[3])&&(OUT[27]>OUT[4])&&(OUT[27]>OUT[5])&&(OUT[27]>OUT[6])&&(OUT[27]>OUT[7])&&(OUT[27]>OUT[8])&&(OUT[27]>OUT[9])
&&(OUT[27]>OUT[10])&&(OUT[27]>OUT[11])&&(OUT[27]>OUT[12])&&(OUT[27]>OUT[13])&&(OUT[27]>OUT[14])&&(OUT[27]>OUT[15])&&(OUT[27]>OUT[16])&&(OUT[27]>OUT[17])&&(OUT[27]>OUT[18])
&&(OUT[27]>OUT[19])&&(OUT[27]>OUT[20])&&(OUT[27]>OUT[21])&&(OUT[27]>OUT[22])&&(OUT[27]>OUT[23])&&(OUT[27]>OUT[24])&&(OUT[27]>OUT[25])&&(OUT[27]>OUT[26])&&(OUT[27]>OUT[28])
&&(OUT[27]>OUT[29])&&(OUT[27]>OUT[30])&&(OUT[27]>OUT[31])&&(OUT[27]>OUT[32])&&(OUT[27]>OUT[33])&&(OUT[27]>OUT[34])&&(OUT[27]>OUT[35])&&(OUT[27]>OUT[36])&&(OUT[27]>OUT[37])
&&(OUT[27]>OUT[38])&&(OUT[27]>OUT[39])&&(OUT[27]>OUT[40])&&(OUT[27]>OUT[41])&&(OUT[27]>OUT[42])&&(OUT[27]>OUT[43])&&(OUT[27]>OUT[44])&&(OUT[27]>OUT[45])&&(OUT[27]>OUT[46])
&&(OUT[27]>OUT[47])&&(OUT[27]>OUT[48]) ) 
{
MT_VCX=0.75;
MT_VCY=0;		
	
}

//29

else if((OUT[28]>OUT[0])&&(OUT[28]>OUT[1])&&(OUT[28]>OUT[2])&&(OUT[28]>OUT[3])&&(OUT[28]>OUT[4])&&(OUT[28]>OUT[5])&&(OUT[28]>OUT[6])&&(OUT[28]>OUT[7])&&(OUT[28]>OUT[8])&&(OUT[28]>OUT[9])
&&(OUT[28]>OUT[10])&&(OUT[28]>OUT[11])&&(OUT[28]>OUT[12])&&(OUT[28]>OUT[13])&&(OUT[28]>OUT[14])&&(OUT[28]>OUT[15])&&(OUT[28]>OUT[16])&&(OUT[28]>OUT[17])&&(OUT[28]>OUT[18])
&&(OUT[28]>OUT[19])&&(OUT[28]>OUT[20])&&(OUT[28]>OUT[21])&&(OUT[28]>OUT[22])&&(OUT[28]>OUT[23])&&(OUT[28]>OUT[24])&&(OUT[28]>OUT[25])&&(OUT[28]>OUT[26])&&(OUT[28]>OUT[27])
&&(OUT[28]>OUT[29])&&(OUT[28]>OUT[30])&&(OUT[28]>OUT[31])&&(OUT[28]>OUT[32])&&(OUT[28]>OUT[33])&&(OUT[28]>OUT[34])&&(OUT[28]>OUT[35])&&(OUT[28]>OUT[36])&&(OUT[28]>OUT[37])
&&(OUT[28]>OUT[38])&&(OUT[28]>OUT[39])&&(OUT[28]>OUT[40])&&(OUT[28]>OUT[41])&&(OUT[28]>OUT[42])&&(OUT[28]>OUT[43])&&(OUT[28]>OUT[44])&&(OUT[28]>OUT[45])&&(OUT[28]>OUT[46])
&&(OUT[28]>OUT[47])&&(OUT[28]>OUT[48]) ) 
{
MT_VCX=-0.75;
MT_VCY=0.25;		
	
}

//30

else if((OUT[29]>OUT[0])&&(OUT[29]>OUT[1])&&(OUT[29]>OUT[2])&&(OUT[29]>OUT[3])&&(OUT[29]>OUT[4])&&(OUT[29]>OUT[5])&&(OUT[29]>OUT[6])&&(OUT[29]>OUT[7])&&(OUT[29]>OUT[8])&&(OUT[29]>OUT[9])
&&(OUT[29]>OUT[10])&&(OUT[29]>OUT[11])&&(OUT[29]>OUT[12])&&(OUT[29]>OUT[13])&&(OUT[29]>OUT[14])&&(OUT[29]>OUT[15])&&(OUT[29]>OUT[16])&&(OUT[29]>OUT[17])&&(OUT[29]>OUT[18])
&&(OUT[29]>OUT[19])&&(OUT[29]>OUT[20])&&(OUT[29]>OUT[21])&&(OUT[29]>OUT[22])&&(OUT[29]>OUT[23])&&(OUT[29]>OUT[24])&&(OUT[29]>OUT[25])&&(OUT[29]>OUT[26])&&(OUT[29]>OUT[27])
&&(OUT[29]>OUT[28])&&(OUT[29]>OUT[30])&&(OUT[29]>OUT[31])&&(OUT[29]>OUT[32])&&(OUT[29]>OUT[33])&&(OUT[29]>OUT[34])&&(OUT[29]>OUT[35])&&(OUT[29]>OUT[36])&&(OUT[29]>OUT[37])
&&(OUT[29]>OUT[38])&&(OUT[29]>OUT[39])&&(OUT[29]>OUT[40])&&(OUT[29]>OUT[41])&&(OUT[29]>OUT[42])&&(OUT[29]>OUT[43])&&(OUT[29]>OUT[44])&&(OUT[29]>OUT[45])&&(OUT[29]>OUT[46])
&&(OUT[29]>OUT[47])&&(OUT[29]>OUT[48]) ) 
{
MT_VCX=-0.5;
MT_VCY=0.25;		
	
}

//31

else if((OUT[30]>OUT[0])&&(OUT[30]>OUT[1])&&(OUT[30]>OUT[2])&&(OUT[30]>OUT[3])&&(OUT[30]>OUT[4])&&(OUT[30]>OUT[5])&&(OUT[30]>OUT[6])&&(OUT[30]>OUT[7])&&(OUT[30]>OUT[8])&&(OUT[30]>OUT[9])
&&(OUT[30]>OUT[10])&&(OUT[30]>OUT[11])&&(OUT[30]>OUT[12])&&(OUT[30]>OUT[13])&&(OUT[30]>OUT[14])&&(OUT[30]>OUT[15])&&(OUT[30]>OUT[16])&&(OUT[30]>OUT[17])&&(OUT[30]>OUT[18])
&&(OUT[30]>OUT[19])&&(OUT[30]>OUT[20])&&(OUT[30]>OUT[21])&&(OUT[30]>OUT[22])&&(OUT[30]>OUT[23])&&(OUT[30]>OUT[24])&&(OUT[30]>OUT[25])&&(OUT[30]>OUT[26])&&(OUT[30]>OUT[27])
&&(OUT[30]>OUT[28])&&(OUT[30]>OUT[29])&&(OUT[30]>OUT[31])&&(OUT[30]>OUT[32])&&(OUT[30]>OUT[33])&&(OUT[30]>OUT[34])&&(OUT[30]>OUT[35])&&(OUT[30]>OUT[36])&&(OUT[30]>OUT[37])
&&(OUT[30]>OUT[38])&&(OUT[30]>OUT[39])&&(OUT[30]>OUT[40])&&(OUT[30]>OUT[41])&&(OUT[30]>OUT[42])&&(OUT[30]>OUT[43])&&(OUT[30]>OUT[44])&&(OUT[30]>OUT[45])&&(OUT[30]>OUT[46])
&&(OUT[30]>OUT[47])&&(OUT[30]>OUT[48]) ) 
{
MT_VCX=-0.25;
MT_VCY=0.25;		
	
}

//32

else if((OUT[31]>OUT[0])&&(OUT[31]>OUT[1])&&(OUT[31]>OUT[2])&&(OUT[31]>OUT[3])&&(OUT[31]>OUT[4])&&(OUT[31]>OUT[5])&&(OUT[31]>OUT[6])&&(OUT[31]>OUT[7])&&(OUT[31]>OUT[8])&&(OUT[31]>OUT[9])
&&(OUT[31]>OUT[10])&&(OUT[31]>OUT[11])&&(OUT[31]>OUT[12])&&(OUT[31]>OUT[13])&&(OUT[31]>OUT[14])&&(OUT[31]>OUT[15])&&(OUT[31]>OUT[16])&&(OUT[31]>OUT[17])&&(OUT[31]>OUT[18])
&&(OUT[31]>OUT[19])&&(OUT[31]>OUT[20])&&(OUT[31]>OUT[21])&&(OUT[31]>OUT[22])&&(OUT[31]>OUT[23])&&(OUT[31]>OUT[24])&&(OUT[31]>OUT[25])&&(OUT[31]>OUT[26])&&(OUT[31]>OUT[27])
&&(OUT[31]>OUT[28])&&(OUT[31]>OUT[29])&&(OUT[31]>OUT[30])&&(OUT[31]>OUT[32])&&(OUT[31]>OUT[33])&&(OUT[31]>OUT[34])&&(OUT[31]>OUT[35])&&(OUT[31]>OUT[36])&&(OUT[31]>OUT[37])
&&(OUT[31]>OUT[38])&&(OUT[31]>OUT[39])&&(OUT[31]>OUT[40])&&(OUT[31]>OUT[41])&&(OUT[31]>OUT[42])&&(OUT[31]>OUT[43])&&(OUT[31]>OUT[44])&&(OUT[31]>OUT[45])&&(OUT[31]>OUT[46])
&&(OUT[31]>OUT[47])&&(OUT[31]>OUT[48]) ) 
{
	
MT_VCX=0;
MT_VCY=0.25;		
}

//33

else if((OUT[32]>OUT[0])&&(OUT[32]>OUT[1])&&(OUT[32]>OUT[2])&&(OUT[32]>OUT[3])&&(OUT[32]>OUT[4])&&(OUT[32]>OUT[5])&&(OUT[32]>OUT[6])&&(OUT[32]>OUT[7])&&(OUT[32]>OUT[8])&&(OUT[32]>OUT[9])
&&(OUT[32]>OUT[10])&&(OUT[32]>OUT[11])&&(OUT[32]>OUT[12])&&(OUT[32]>OUT[13])&&(OUT[32]>OUT[14])&&(OUT[32]>OUT[15])&&(OUT[32]>OUT[16])&&(OUT[32]>OUT[17])&&(OUT[32]>OUT[18])
&&(OUT[32]>OUT[19])&&(OUT[32]>OUT[20])&&(OUT[32]>OUT[21])&&(OUT[32]>OUT[22])&&(OUT[32]>OUT[23])&&(OUT[32]>OUT[24])&&(OUT[32]>OUT[25])&&(OUT[32]>OUT[26])&&(OUT[32]>OUT[27])
&&(OUT[32]>OUT[28])&&(OUT[32]>OUT[29])&&(OUT[32]>OUT[30])&&(OUT[32]>OUT[31])&&(OUT[32]>OUT[33])&&(OUT[32]>OUT[34])&&(OUT[32]>OUT[35])&&(OUT[32]>OUT[36])&&(OUT[32]>OUT[37])
&&(OUT[32]>OUT[38])&&(OUT[32]>OUT[39])&&(OUT[32]>OUT[40])&&(OUT[32]>OUT[41])&&(OUT[32]>OUT[42])&&(OUT[32]>OUT[43])&&(OUT[32]>OUT[44])&&(OUT[32]>OUT[45])&&(OUT[32]>OUT[46])
&&(OUT[32]>OUT[47])&&(OUT[32]>OUT[48]) ) 
{
MT_VCX=0.25;
MT_VCY=0.25;		
	
}

//34

else if((OUT[33]>OUT[0])&&(OUT[33]>OUT[1])&&(OUT[33]>OUT[2])&&(OUT[33]>OUT[3])&&(OUT[33]>OUT[4])&&(OUT[33]>OUT[5])&&(OUT[33]>OUT[6])&&(OUT[33]>OUT[7])&&(OUT[33]>OUT[8])&&(OUT[33]>OUT[9])
&&(OUT[33]>OUT[10])&&(OUT[33]>OUT[11])&&(OUT[33]>OUT[12])&&(OUT[33]>OUT[13])&&(OUT[33]>OUT[14])&&(OUT[33]>OUT[15])&&(OUT[33]>OUT[16])&&(OUT[33]>OUT[17])&&(OUT[33]>OUT[18])
&&(OUT[33]>OUT[19])&&(OUT[33]>OUT[20])&&(OUT[33]>OUT[21])&&(OUT[33]>OUT[22])&&(OUT[33]>OUT[23])&&(OUT[33]>OUT[24])&&(OUT[33]>OUT[25])&&(OUT[33]>OUT[26])&&(OUT[33]>OUT[27])
&&(OUT[33]>OUT[28])&&(OUT[33]>OUT[29])&&(OUT[33]>OUT[30])&&(OUT[33]>OUT[31])&&(OUT[33]>OUT[32])&&(OUT[33]>OUT[34])&&(OUT[33]>OUT[35])&&(OUT[33]>OUT[36])&&(OUT[33]>OUT[37])
&&(OUT[33]>OUT[38])&&(OUT[33]>OUT[39])&&(OUT[33]>OUT[40])&&(OUT[33]>OUT[41])&&(OUT[33]>OUT[42])&&(OUT[33]>OUT[43])&&(OUT[33]>OUT[44])&&(OUT[33]>OUT[45])&&(OUT[33]>OUT[46])
&&(OUT[33]>OUT[47])&&(OUT[33]>OUT[48]) ) 
{
MT_VCX=0.5;
MT_VCY=0.25;		
	
}


//35

else if((OUT[34]>OUT[0])&&(OUT[34]>OUT[1])&&(OUT[34]>OUT[2])&&(OUT[34]>OUT[3])&&(OUT[34]>OUT[4])&&(OUT[34]>OUT[5])&&(OUT[34]>OUT[6])&&(OUT[34]>OUT[7])&&(OUT[34]>OUT[8])&&(OUT[34]>OUT[9])
&&(OUT[34]>OUT[10])&&(OUT[34]>OUT[11])&&(OUT[34]>OUT[12])&&(OUT[34]>OUT[13])&&(OUT[34]>OUT[14])&&(OUT[34]>OUT[15])&&(OUT[34]>OUT[16])&&(OUT[34]>OUT[17])&&(OUT[34]>OUT[18])
&&(OUT[34]>OUT[19])&&(OUT[34]>OUT[20])&&(OUT[34]>OUT[21])&&(OUT[34]>OUT[22])&&(OUT[34]>OUT[23])&&(OUT[34]>OUT[24])&&(OUT[34]>OUT[25])&&(OUT[34]>OUT[26])&&(OUT[34]>OUT[27])
&&(OUT[34]>OUT[28])&&(OUT[34]>OUT[29])&&(OUT[34]>OUT[30])&&(OUT[34]>OUT[31])&&(OUT[34]>OUT[32])&&(OUT[34]>OUT[33])&&(OUT[34]>OUT[35])&&(OUT[34]>OUT[36])&&(OUT[34]>OUT[37])
&&(OUT[34]>OUT[38])&&(OUT[34]>OUT[39])&&(OUT[34]>OUT[40])&&(OUT[34]>OUT[41])&&(OUT[34]>OUT[42])&&(OUT[34]>OUT[43])&&(OUT[34]>OUT[44])&&(OUT[34]>OUT[45])&&(OUT[34]>OUT[46])
&&(OUT[34]>OUT[47])&&(OUT[34]>OUT[48]) ) 
{
MT_VCX=0.75;
MT_VCY=0.25;		
	
}


//36

else if((OUT[35]>OUT[0])&&(OUT[35]>OUT[1])&&(OUT[35]>OUT[2])&&(OUT[35]>OUT[3])&&(OUT[35]>OUT[4])&&(OUT[35]>OUT[5])&&(OUT[35]>OUT[6])&&(OUT[35]>OUT[7])&&(OUT[35]>OUT[8])&&(OUT[35]>OUT[9])
&&(OUT[35]>OUT[10])&&(OUT[35]>OUT[11])&&(OUT[35]>OUT[12])&&(OUT[35]>OUT[13])&&(OUT[35]>OUT[14])&&(OUT[35]>OUT[15])&&(OUT[35]>OUT[16])&&(OUT[35]>OUT[17])&&(OUT[35]>OUT[18])
&&(OUT[35]>OUT[19])&&(OUT[35]>OUT[20])&&(OUT[35]>OUT[21])&&(OUT[35]>OUT[22])&&(OUT[35]>OUT[23])&&(OUT[35]>OUT[24])&&(OUT[35]>OUT[25])&&(OUT[35]>OUT[26])&&(OUT[35]>OUT[27])
&&(OUT[35]>OUT[28])&&(OUT[35]>OUT[29])&&(OUT[35]>OUT[30])&&(OUT[35]>OUT[31])&&(OUT[35]>OUT[32])&&(OUT[35]>OUT[33])&&(OUT[35]>OUT[34])&&(OUT[35]>OUT[36])&&(OUT[35]>OUT[37])
&&(OUT[35]>OUT[38])&&(OUT[35]>OUT[39])&&(OUT[35]>OUT[40])&&(OUT[35]>OUT[41])&&(OUT[35]>OUT[42])&&(OUT[35]>OUT[43])&&(OUT[35]>OUT[44])&&(OUT[35]>OUT[45])&&(OUT[35]>OUT[46])
&&(OUT[35]>OUT[47])&&(OUT[35]>OUT[48]) ) 
{
MT_VCX=-0.75;
MT_VCY=0.5;		
	
}

//37

else if((OUT[36]>OUT[0])&&(OUT[36]>OUT[1])&&(OUT[36]>OUT[2])&&(OUT[36]>OUT[3])&&(OUT[36]>OUT[4])&&(OUT[36]>OUT[5])&&(OUT[36]>OUT[6])&&(OUT[36]>OUT[7])&&(OUT[36]>OUT[8])&&(OUT[36]>OUT[9])
&&(OUT[36]>OUT[10])&&(OUT[36]>OUT[11])&&(OUT[36]>OUT[12])&&(OUT[36]>OUT[13])&&(OUT[36]>OUT[14])&&(OUT[36]>OUT[15])&&(OUT[36]>OUT[16])&&(OUT[36]>OUT[17])&&(OUT[36]>OUT[18])
&&(OUT[36]>OUT[19])&&(OUT[36]>OUT[20])&&(OUT[36]>OUT[21])&&(OUT[36]>OUT[22])&&(OUT[36]>OUT[23])&&(OUT[36]>OUT[24])&&(OUT[36]>OUT[25])&&(OUT[36]>OUT[26])&&(OUT[36]>OUT[27])
&&(OUT[36]>OUT[28])&&(OUT[36]>OUT[29])&&(OUT[36]>OUT[30])&&(OUT[36]>OUT[31])&&(OUT[36]>OUT[32])&&(OUT[36]>OUT[33])&&(OUT[36]>OUT[34])&&(OUT[36]>OUT[35])&&(OUT[36]>OUT[37])
&&(OUT[36]>OUT[38])&&(OUT[36]>OUT[39])&&(OUT[36]>OUT[40])&&(OUT[36]>OUT[41])&&(OUT[36]>OUT[42])&&(OUT[36]>OUT[43])&&(OUT[36]>OUT[44])&&(OUT[36]>OUT[45])&&(OUT[36]>OUT[46])
&&(OUT[36]>OUT[47])&&(OUT[36]>OUT[48]) ) 
{
MT_VCX=-0.5;
MT_VCY=0.5;		
	
}

//38

else if((OUT[37]>OUT[0])&&(OUT[37]>OUT[1])&&(OUT[37]>OUT[2])&&(OUT[37]>OUT[3])&&(OUT[37]>OUT[4])&&(OUT[37]>OUT[5])&&(OUT[37]>OUT[6])&&(OUT[37]>OUT[7])&&(OUT[37]>OUT[8])&&(OUT[37]>OUT[9])
&&(OUT[37]>OUT[10])&&(OUT[37]>OUT[11])&&(OUT[37]>OUT[12])&&(OUT[37]>OUT[13])&&(OUT[37]>OUT[14])&&(OUT[37]>OUT[15])&&(OUT[37]>OUT[16])&&(OUT[37]>OUT[17])&&(OUT[37]>OUT[18])
&&(OUT[37]>OUT[19])&&(OUT[37]>OUT[20])&&(OUT[37]>OUT[21])&&(OUT[37]>OUT[22])&&(OUT[37]>OUT[23])&&(OUT[37]>OUT[24])&&(OUT[37]>OUT[25])&&(OUT[37]>OUT[26])&&(OUT[37]>OUT[27])
&&(OUT[37]>OUT[28])&&(OUT[37]>OUT[29])&&(OUT[37]>OUT[30])&&(OUT[37]>OUT[31])&&(OUT[37]>OUT[32])&&(OUT[37]>OUT[33])&&(OUT[37]>OUT[34])&&(OUT[37]>OUT[35])&&(OUT[37]>OUT[36])
&&(OUT[37]>OUT[38])&&(OUT[37]>OUT[39])&&(OUT[37]>OUT[40])&&(OUT[37]>OUT[41])&&(OUT[37]>OUT[42])&&(OUT[37]>OUT[43])&&(OUT[37]>OUT[44])&&(OUT[37]>OUT[45])&&(OUT[37]>OUT[46])
&&(OUT[37]>OUT[47])&&(OUT[37]>OUT[48]) ) 
{
MT_VCX=-0.25;
MT_VCY=0.5;		
	
}
//39

else if((OUT[38]>OUT[0])&&(OUT[38]>OUT[1])&&(OUT[38]>OUT[2])&&(OUT[38]>OUT[3])&&(OUT[38]>OUT[4])&&(OUT[38]>OUT[5])&&(OUT[38]>OUT[6])&&(OUT[38]>OUT[7])&&(OUT[38]>OUT[8])&&(OUT[38]>OUT[9])
&&(OUT[38]>OUT[10])&&(OUT[38]>OUT[11])&&(OUT[38]>OUT[12])&&(OUT[38]>OUT[13])&&(OUT[38]>OUT[14])&&(OUT[38]>OUT[15])&&(OUT[38]>OUT[16])&&(OUT[38]>OUT[17])&&(OUT[38]>OUT[18])
&&(OUT[38]>OUT[19])&&(OUT[38]>OUT[20])&&(OUT[38]>OUT[21])&&(OUT[38]>OUT[22])&&(OUT[38]>OUT[23])&&(OUT[38]>OUT[24])&&(OUT[38]>OUT[25])&&(OUT[38]>OUT[26])&&(OUT[38]>OUT[27])
&&(OUT[38]>OUT[28])&&(OUT[38]>OUT[29])&&(OUT[38]>OUT[30])&&(OUT[38]>OUT[31])&&(OUT[38]>OUT[32])&&(OUT[38]>OUT[33])&&(OUT[38]>OUT[34])&&(OUT[38]>OUT[35])&&(OUT[38]>OUT[36])
&&(OUT[38]>OUT[37])&&(OUT[38]>OUT[39])&&(OUT[38]>OUT[40])&&(OUT[38]>OUT[41])&&(OUT[38]>OUT[42])&&(OUT[38]>OUT[43])&&(OUT[38]>OUT[44])&&(OUT[38]>OUT[45])&&(OUT[38]>OUT[46])
&&(OUT[38]>OUT[47])&&(OUT[38]>OUT[48]) ) 
{
MT_VCX=0;
MT_VCY=0.5;		
	
}

//40
else if((OUT[39]>OUT[0])&&(OUT[39]>OUT[1])&&(OUT[39]>OUT[2])&&(OUT[39]>OUT[3])&&(OUT[39]>OUT[4])&&(OUT[39]>OUT[5])&&(OUT[39]>OUT[6])&&(OUT[39]>OUT[7])&&(OUT[39]>OUT[8])&&(OUT[39]>OUT[9])
&&(OUT[39]>OUT[10])&&(OUT[39]>OUT[11])&&(OUT[39]>OUT[12])&&(OUT[39]>OUT[13])&&(OUT[39]>OUT[14])&&(OUT[39]>OUT[15])&&(OUT[39]>OUT[16])&&(OUT[39]>OUT[17])&&(OUT[39]>OUT[18])
&&(OUT[39]>OUT[19])&&(OUT[39]>OUT[20])&&(OUT[39]>OUT[21])&&(OUT[39]>OUT[22])&&(OUT[39]>OUT[23])&&(OUT[39]>OUT[24])&&(OUT[39]>OUT[25])&&(OUT[39]>OUT[26])&&(OUT[39]>OUT[27])
&&(OUT[39]>OUT[28])&&(OUT[39]>OUT[29])&&(OUT[39]>OUT[30])&&(OUT[39]>OUT[31])&&(OUT[39]>OUT[32])&&(OUT[39]>OUT[33])&&(OUT[39]>OUT[34])&&(OUT[39]>OUT[35])&&(OUT[39]>OUT[36])
&&(OUT[39]>OUT[37])&&(OUT[39]>OUT[38])&&(OUT[39]>OUT[40])&&(OUT[39]>OUT[41])&&(OUT[39]>OUT[42])&&(OUT[39]>OUT[43])&&(OUT[39]>OUT[44])&&(OUT[39]>OUT[45])&&(OUT[39]>OUT[46])
&&(OUT[39]>OUT[47])&&(OUT[39]>OUT[48]) ) 
{
MT_VCX=0.25;
MT_VCY=0.5;		
	
}

//41

else if((OUT[40]>OUT[0])&&(OUT[40]>OUT[1])&&(OUT[40]>OUT[2])&&(OUT[40]>OUT[3])&&(OUT[40]>OUT[4])&&(OUT[40]>OUT[5])&&(OUT[40]>OUT[6])&&(OUT[40]>OUT[7])&&(OUT[40]>OUT[8])&&(OUT[40]>OUT[9])
&&(OUT[40]>OUT[10])&&(OUT[40]>OUT[11])&&(OUT[40]>OUT[12])&&(OUT[40]>OUT[13])&&(OUT[40]>OUT[14])&&(OUT[40]>OUT[15])&&(OUT[40]>OUT[16])&&(OUT[40]>OUT[17])&&(OUT[40]>OUT[18])
&&(OUT[40]>OUT[19])&&(OUT[40]>OUT[20])&&(OUT[40]>OUT[21])&&(OUT[40]>OUT[22])&&(OUT[40]>OUT[23])&&(OUT[40]>OUT[24])&&(OUT[40]>OUT[25])&&(OUT[40]>OUT[26])&&(OUT[40]>OUT[27])
&&(OUT[40]>OUT[28])&&(OUT[40]>OUT[29])&&(OUT[40]>OUT[30])&&(OUT[40]>OUT[31])&&(OUT[40]>OUT[32])&&(OUT[40]>OUT[33])&&(OUT[40]>OUT[34])&&(OUT[40]>OUT[35])&&(OUT[40]>OUT[36])
&&(OUT[40]>OUT[37])&&(OUT[40]>OUT[38])&&(OUT[40]>OUT[39])&&(OUT[40]>OUT[41])&&(OUT[40]>OUT[42])&&(OUT[40]>OUT[43])&&(OUT[40]>OUT[44])&&(OUT[40]>OUT[45])&&(OUT[40]>OUT[46])
&&(OUT[40]>OUT[47])&&(OUT[40]>OUT[48]) ) 
{
MT_VCX=0.5;
MT_VCY=0.5;		
	
}

//42

else if((OUT[41]>OUT[0])&&(OUT[41]>OUT[1])&&(OUT[41]>OUT[2])&&(OUT[41]>OUT[3])&&(OUT[41]>OUT[4])&&(OUT[41]>OUT[5])&&(OUT[41]>OUT[6])&&(OUT[41]>OUT[7])&&(OUT[41]>OUT[8])&&(OUT[41]>OUT[9])
&&(OUT[41]>OUT[10])&&(OUT[41]>OUT[11])&&(OUT[41]>OUT[12])&&(OUT[41]>OUT[13])&&(OUT[41]>OUT[14])&&(OUT[41]>OUT[15])&&(OUT[41]>OUT[16])&&(OUT[41]>OUT[17])&&(OUT[41]>OUT[18])
&&(OUT[41]>OUT[19])&&(OUT[41]>OUT[20])&&(OUT[41]>OUT[21])&&(OUT[41]>OUT[22])&&(OUT[41]>OUT[23])&&(OUT[41]>OUT[24])&&(OUT[41]>OUT[25])&&(OUT[41]>OUT[26])&&(OUT[41]>OUT[27])
&&(OUT[41]>OUT[28])&&(OUT[41]>OUT[29])&&(OUT[41]>OUT[30])&&(OUT[41]>OUT[31])&&(OUT[41]>OUT[32])&&(OUT[41]>OUT[33])&&(OUT[41]>OUT[34])&&(OUT[41]>OUT[35])&&(OUT[41]>OUT[36])
&&(OUT[41]>OUT[37])&&(OUT[41]>OUT[38])&&(OUT[41]>OUT[39])&&(OUT[41]>OUT[40])&&(OUT[41]>OUT[42])&&(OUT[41]>OUT[43])&&(OUT[41]>OUT[44])&&(OUT[41]>OUT[45])&&(OUT[41]>OUT[46])
&&(OUT[41]>OUT[47])&&(OUT[41]>OUT[48]) ) 
{
MT_VCX=0.75;
MT_VCY=0.5;		
	
}
//43

else if((OUT[42]>OUT[0])&&(OUT[42]>OUT[1])&&(OUT[42]>OUT[2])&&(OUT[42]>OUT[3])&&(OUT[42]>OUT[4])&&(OUT[42]>OUT[5])&&(OUT[42]>OUT[6])&&(OUT[42]>OUT[7])&&(OUT[42]>OUT[8])&&(OUT[42]>OUT[9])
&&(OUT[42]>OUT[10])&&(OUT[42]>OUT[11])&&(OUT[42]>OUT[12])&&(OUT[42]>OUT[13])&&(OUT[42]>OUT[14])&&(OUT[42]>OUT[15])&&(OUT[42]>OUT[16])&&(OUT[42]>OUT[17])&&(OUT[42]>OUT[18])
&&(OUT[42]>OUT[19])&&(OUT[42]>OUT[20])&&(OUT[42]>OUT[21])&&(OUT[42]>OUT[22])&&(OUT[42]>OUT[23])&&(OUT[42]>OUT[24])&&(OUT[42]>OUT[25])&&(OUT[42]>OUT[26])&&(OUT[42]>OUT[27])
&&(OUT[42]>OUT[28])&&(OUT[42]>OUT[29])&&(OUT[42]>OUT[30])&&(OUT[42]>OUT[31])&&(OUT[42]>OUT[32])&&(OUT[42]>OUT[33])&&(OUT[42]>OUT[34])&&(OUT[42]>OUT[35])&&(OUT[42]>OUT[36])
&&(OUT[42]>OUT[37])&&(OUT[42]>OUT[38])&&(OUT[42]>OUT[39])&&(OUT[42]>OUT[40])&&(OUT[42]>OUT[41])&&(OUT[42]>OUT[43])&&(OUT[42]>OUT[44])&&(OUT[42]>OUT[45])&&(OUT[42]>OUT[46])
&&(OUT[42]>OUT[47])&&(OUT[42]>OUT[48]) ) 
{
MT_VCX=-0.75;
MT_VCY=0.75;		
	
}

//44

else if((OUT[43]>OUT[0])&&(OUT[43]>OUT[1])&&(OUT[43]>OUT[2])&&(OUT[43]>OUT[3])&&(OUT[43]>OUT[4])&&(OUT[43]>OUT[5])&&(OUT[43]>OUT[6])&&(OUT[43]>OUT[7])&&(OUT[43]>OUT[8])&&(OUT[43]>OUT[9])
&&(OUT[43]>OUT[10])&&(OUT[43]>OUT[11])&&(OUT[43]>OUT[12])&&(OUT[43]>OUT[13])&&(OUT[43]>OUT[14])&&(OUT[43]>OUT[15])&&(OUT[43]>OUT[16])&&(OUT[43]>OUT[17])&&(OUT[43]>OUT[18])
&&(OUT[43]>OUT[19])&&(OUT[43]>OUT[20])&&(OUT[43]>OUT[21])&&(OUT[43]>OUT[22])&&(OUT[43]>OUT[23])&&(OUT[43]>OUT[24])&&(OUT[43]>OUT[25])&&(OUT[43]>OUT[26])&&(OUT[43]>OUT[27])
&&(OUT[43]>OUT[28])&&(OUT[43]>OUT[29])&&(OUT[43]>OUT[30])&&(OUT[43]>OUT[31])&&(OUT[43]>OUT[32])&&(OUT[43]>OUT[33])&&(OUT[43]>OUT[34])&&(OUT[43]>OUT[35])&&(OUT[43]>OUT[36])
&&(OUT[43]>OUT[37])&&(OUT[43]>OUT[38])&&(OUT[43]>OUT[39])&&(OUT[43]>OUT[40])&&(OUT[43]>OUT[41])&&(OUT[43]>OUT[42])&&(OUT[43]>OUT[44])&&(OUT[43]>OUT[45])&&(OUT[43]>OUT[46])
&&(OUT[43]>OUT[47])&&(OUT[43]>OUT[48]) ) 
{
MT_VCX=-0.5;
MT_VCY=0.75;	
	
}

//45
else if((OUT[44]>OUT[0])&&(OUT[44]>OUT[1])&&(OUT[44]>OUT[2])&&(OUT[44]>OUT[3])&&(OUT[44]>OUT[4])&&(OUT[44]>OUT[5])&&(OUT[44]>OUT[6])&&(OUT[44]>OUT[7])&&(OUT[44]>OUT[8])&&(OUT[44]>OUT[9])
&&(OUT[44]>OUT[10])&&(OUT[44]>OUT[11])&&(OUT[44]>OUT[12])&&(OUT[44]>OUT[13])&&(OUT[44]>OUT[14])&&(OUT[44]>OUT[15])&&(OUT[44]>OUT[16])&&(OUT[44]>OUT[17])&&(OUT[44]>OUT[18])
&&(OUT[44]>OUT[19])&&(OUT[44]>OUT[20])&&(OUT[44]>OUT[21])&&(OUT[44]>OUT[22])&&(OUT[44]>OUT[23])&&(OUT[44]>OUT[24])&&(OUT[44]>OUT[25])&&(OUT[44]>OUT[26])&&(OUT[44]>OUT[27])
&&(OUT[44]>OUT[28])&&(OUT[44]>OUT[29])&&(OUT[44]>OUT[30])&&(OUT[44]>OUT[31])&&(OUT[44]>OUT[32])&&(OUT[44]>OUT[33])&&(OUT[44]>OUT[34])&&(OUT[44]>OUT[35])&&(OUT[44]>OUT[36])
&&(OUT[44]>OUT[37])&&(OUT[44]>OUT[38])&&(OUT[44]>OUT[39])&&(OUT[44]>OUT[40])&&(OUT[44]>OUT[41])&&(OUT[44]>OUT[42])&&(OUT[44]>OUT[43])&&(OUT[44]>OUT[45])&&(OUT[44]>OUT[46])
&&(OUT[44]>OUT[47])&&(OUT[44]>OUT[48]) ) 
{
MT_VCX=-0.25;
MT_VCY=0.75;	
	
}

//46
else if((OUT[45]>OUT[0])&&(OUT[45]>OUT[1])&&(OUT[45]>OUT[2])&&(OUT[45]>OUT[3])&&(OUT[45]>OUT[4])&&(OUT[45]>OUT[5])&&(OUT[45]>OUT[6])&&(OUT[45]>OUT[7])&&(OUT[45]>OUT[8])&&(OUT[45]>OUT[9])
&&(OUT[45]>OUT[10])&&(OUT[45]>OUT[11])&&(OUT[45]>OUT[12])&&(OUT[45]>OUT[13])&&(OUT[45]>OUT[14])&&(OUT[45]>OUT[15])&&(OUT[45]>OUT[16])&&(OUT[45]>OUT[17])&&(OUT[45]>OUT[18])
&&(OUT[45]>OUT[19])&&(OUT[45]>OUT[20])&&(OUT[45]>OUT[21])&&(OUT[45]>OUT[22])&&(OUT[45]>OUT[23])&&(OUT[45]>OUT[24])&&(OUT[45]>OUT[25])&&(OUT[45]>OUT[26])&&(OUT[45]>OUT[27])
&&(OUT[45]>OUT[28])&&(OUT[45]>OUT[29])&&(OUT[45]>OUT[30])&&(OUT[45]>OUT[31])&&(OUT[45]>OUT[32])&&(OUT[45]>OUT[33])&&(OUT[45]>OUT[34])&&(OUT[45]>OUT[35])&&(OUT[45]>OUT[36])
&&(OUT[45]>OUT[37])&&(OUT[45]>OUT[38])&&(OUT[45]>OUT[39])&&(OUT[45]>OUT[40])&&(OUT[45]>OUT[41])&&(OUT[45]>OUT[42])&&(OUT[45]>OUT[43])&&(OUT[45]>OUT[44])&&(OUT[45]>OUT[46])
&&(OUT[45]>OUT[47])&&(OUT[45]>OUT[48]) ) 
{
MT_VCX=0;
MT_VCY=0.75;
	
}
//47

else if((OUT[46]>OUT[0])&&(OUT[46]>OUT[1])&&(OUT[46]>OUT[2])&&(OUT[46]>OUT[3])&&(OUT[46]>OUT[4])&&(OUT[46]>OUT[5])&&(OUT[46]>OUT[6])&&(OUT[46]>OUT[7])&&(OUT[46]>OUT[8])&&(OUT[46]>OUT[9])
&&(OUT[46]>OUT[10])&&(OUT[46]>OUT[11])&&(OUT[46]>OUT[12])&&(OUT[46]>OUT[13])&&(OUT[46]>OUT[14])&&(OUT[46]>OUT[15])&&(OUT[46]>OUT[16])&&(OUT[46]>OUT[17])&&(OUT[46]>OUT[18])
&&(OUT[46]>OUT[19])&&(OUT[46]>OUT[20])&&(OUT[46]>OUT[21])&&(OUT[46]>OUT[22])&&(OUT[46]>OUT[23])&&(OUT[46]>OUT[24])&&(OUT[46]>OUT[25])&&(OUT[46]>OUT[26])&&(OUT[46]>OUT[27])
&&(OUT[46]>OUT[28])&&(OUT[46]>OUT[29])&&(OUT[46]>OUT[30])&&(OUT[46]>OUT[31])&&(OUT[46]>OUT[32])&&(OUT[46]>OUT[33])&&(OUT[46]>OUT[34])&&(OUT[46]>OUT[35])&&(OUT[46]>OUT[36])
&&(OUT[46]>OUT[37])&&(OUT[46]>OUT[38])&&(OUT[46]>OUT[39])&&(OUT[46]>OUT[40])&&(OUT[46]>OUT[41])&&(OUT[46]>OUT[42])&&(OUT[46]>OUT[43])&&(OUT[46]>OUT[44])&&(OUT[46]>OUT[45])
&&(OUT[46]>OUT[47])&&(OUT[46]>OUT[48]) ) 
{
MT_VCX=0.25;
MT_VCY=0.75;	
	
}

//48
else if((OUT[47]>OUT[0])&&(OUT[47]>OUT[1])&&(OUT[47]>OUT[2])&&(OUT[47]>OUT[3])&&(OUT[47]>OUT[4])&&(OUT[47]>OUT[5])&&(OUT[47]>OUT[6])&&(OUT[47]>OUT[7])&&(OUT[47]>OUT[8])&&(OUT[47]>OUT[9])
&&(OUT[47]>OUT[10])&&(OUT[47]>OUT[11])&&(OUT[47]>OUT[12])&&(OUT[47]>OUT[13])&&(OUT[47]>OUT[14])&&(OUT[47]>OUT[15])&&(OUT[47]>OUT[16])&&(OUT[47]>OUT[17])&&(OUT[47]>OUT[18])
&&(OUT[47]>OUT[19])&&(OUT[47]>OUT[20])&&(OUT[47]>OUT[21])&&(OUT[47]>OUT[22])&&(OUT[47]>OUT[23])&&(OUT[47]>OUT[24])&&(OUT[47]>OUT[25])&&(OUT[47]>OUT[26])&&(OUT[47]>OUT[27])
&&(OUT[47]>OUT[28])&&(OUT[47]>OUT[29])&&(OUT[47]>OUT[30])&&(OUT[47]>OUT[31])&&(OUT[47]>OUT[32])&&(OUT[47]>OUT[33])&&(OUT[47]>OUT[34])&&(OUT[47]>OUT[35])&&(OUT[47]>OUT[36])
&&(OUT[47]>OUT[37])&&(OUT[47]>OUT[38])&&(OUT[47]>OUT[39])&&(OUT[47]>OUT[40])&&(OUT[47]>OUT[41])&&(OUT[47]>OUT[42])&&(OUT[47]>OUT[43])&&(OUT[47]>OUT[44])&&(OUT[47]>OUT[45])
&&(OUT[47]>OUT[46])&&(OUT[47]>OUT[48]) ) 
{
MT_VCX=0.5;
MT_VCY=0.75;	
	
}

//49
else if((OUT[48]>OUT[0])&&(OUT[48]>OUT[1])&&(OUT[48]>OUT[2])&&(OUT[48]>OUT[3])&&(OUT[48]>OUT[4])&&(OUT[48]>OUT[5])&&(OUT[48]>OUT[6])&&(OUT[48]>OUT[7])&&(OUT[48]>OUT[8])&&(OUT[48]>OUT[9])
&&(OUT[48]>OUT[10])&&(OUT[48]>OUT[11])&&(OUT[48]>OUT[12])&&(OUT[48]>OUT[13])&&(OUT[48]>OUT[14])&&(OUT[48]>OUT[15])&&(OUT[48]>OUT[16])&&(OUT[48]>OUT[17])&&(OUT[48]>OUT[18])
&&(OUT[48]>OUT[19])&&(OUT[48]>OUT[20])&&(OUT[48]>OUT[21])&&(OUT[48]>OUT[22])&&(OUT[48]>OUT[23])&&(OUT[48]>OUT[24])&&(OUT[48]>OUT[25])&&(OUT[48]>OUT[26])&&(OUT[48]>OUT[27])
&&(OUT[48]>OUT[28])&&(OUT[48]>OUT[29])&&(OUT[48]>OUT[30])&&(OUT[48]>OUT[31])&&(OUT[48]>OUT[32])&&(OUT[48]>OUT[33])&&(OUT[48]>OUT[34])&&(OUT[48]>OUT[35])&&(OUT[48]>OUT[36])
&&(OUT[48]>OUT[37])&&(OUT[48]>OUT[38])&&(OUT[48]>OUT[39])&&(OUT[48]>OUT[40])&&(OUT[48]>OUT[41])&&(OUT[48]>OUT[42])&&(OUT[48]>OUT[43])&&(OUT[48]>OUT[44])&&(OUT[48]>OUT[45])
&&(OUT[48]>OUT[46])&&(OUT[48]>OUT[47]) ) 
{
MT_VCX=0.75;
MT_VCY=0.75;	
	
}


// final result 

if (MT_VCX == 0)

{

	MVX_HALF = 0;
	MVX_QRTER = 0;

}


else if (MT_VCX == 0.25)

{



	MVX_HALF = 0;
	MVX_QRTER = 1;

}

else if (MT_VCX == 0.5)

{


	MVX_HALF = 1;
	MVX_QRTER = 0;

}
else if (MT_VCX == 0.75)

{


	MVX_HALF = 1;
	MVX_QRTER = 1;

}

else if (MT_VCX == -0.5)

{


	MVX_HALF = -1;
	MVX_QRTER = 0;

}

else if (MT_VCX == -0.25)

{


	MVX_HALF = 0;
	MVX_QRTER = -1;

}

else if (MT_VCX == -0.75)

{


	MVX_HALF = -1;
	MVX_QRTER = -1;

}

//MVY
if (MT_VCY == 0)

{

	MVY_HALF = 0;
	MVY_QRTER = 0;

}


else if (MT_VCY == 0.25)

{



	MVY_HALF = 0;
	MVY_QRTER = 1;

}

else if (MT_VCY == 0.5)

{

	MVY_HALF = 1;
	MVY_QRTER = 0;

}
else if (MT_VCY == 0.75)

{

	MVY_HALF = 1;
	MVY_QRTER = 1;

}


else if (MT_VCY == -0.5)

{


	MVY_HALF = -1;
	MVY_QRTER = 0;

}

else if (MT_VCY == -0.25)

{


	MVY_HALF = 0;
	MVY_QRTER = -1;

}
else if (MT_VCY == -0.75)

{


	MVY_HALF = -1;
	MVY_QRTER = -1;

}







	
	  //end of neural network code


	  for (int k = 0; k <= counter_i; k++)
	  {
		  array[k] = 0;
	  }

	  counter_i = 0;
	  index_ref = 0;
	// myfile2 << flag_start << ',' << flag_2point << ',' << flag_star <<  endl;

	  flag_start = 0;
	  flag_2point = 0;
	  flag_star = 0;

      break;



    case MESEARCH_SELECTIVE:
      xTZSearchSelective( pcCU, pcPatternKey, piRefY, iRefStride, pcMvSrchRngLT, pcMvSrchRngRB, rcMv, ruiSAD, pIntegerMv2Nx2NPred );
      break;

    case MESEARCH_DIAMOND_ENHANCED:
      xTZSearch( pcCU, pcPatternKey, piRefY, iRefStride, pcMvSrchRngLT, pcMvSrchRngRB, rcMv, ruiSAD, pIntegerMv2Nx2NPred, true );
      break;

    case MESEARCH_FULL: // shouldn't get here.
    default:
      break;
  }
}


Void TEncSearch::xTZSearch( const TComDataCU* const pcCU,
                            const TComPattern* const pcPatternKey,
                            const Pel* const         piRefY,
                            const Int                iRefStride,
                            const TComMv* const      pcMvSrchRngLT,
                            const TComMv* const      pcMvSrchRngRB,
                            TComMv&                  rcMv,
                            Distortion&              ruiSAD,
                            const TComMv* const      pIntegerMv2Nx2NPred,
                            const Bool               bExtendedSettings)
{
  const Bool bUseAdaptiveRaster                      = bExtendedSettings;
  const Int  iRaster                                 = 5;
  const Bool bTestOtherPredictedMV                   = bExtendedSettings;
  const Bool bTestZeroVector                         = true;
  const Bool bTestZeroVectorStart                    = bExtendedSettings;
  const Bool bTestZeroVectorStop                     = false;
  const Bool bFirstSearchDiamond                     = true;  // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const Bool bFirstCornersForDiamondDist1            = bExtendedSettings;
  const Bool bFirstSearchStop                        = m_pcEncCfg->getFastMEAssumingSmootherMVEnabled();
  const UInt uiFirstSearchRounds                     = 3;     // first search stop X rounds after best match (must be >=1)
  const Bool bEnableRasterSearch                     = true;
  const Bool bAlwaysRasterSearch                     = bExtendedSettings;  // true: BETTER but factor 2 slower
  const Bool bRasterRefinementEnable                 = false; // enable either raster refinement or star refinement
  const Bool bRasterRefinementDiamond                = false; // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const Bool bRasterRefinementCornersForDiamondDist1 = bExtendedSettings;
  const Bool bStarRefinementEnable                   = true;  // enable either star refinement or raster refinement
  const Bool bStarRefinementDiamond                  = true;  // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const Bool bStarRefinementCornersForDiamondDist1   = bExtendedSettings;
  const Bool bStarRefinementStop                     = false;
  const UInt uiStarRefinementRounds                  = 2;  // star refinement stop X rounds after best match (must be >=1)
  const Bool bNewZeroNeighbourhoodTest               = bExtendedSettings;

  UInt uiSearchRange = m_iSearchRange;
  pcCU->clipMv( rcMv );
#if ME_ENABLE_ROUNDING_OF_MVS
  rcMv.divideByPowerOf2(2);
#else
  rcMv >>= 2;
#endif
  // init TZSearchStruct
  IntTZSearchStruct cStruct;
  cStruct.iYStride    = iRefStride;
  cStruct.piRefY      = piRefY;
  cStruct.uiBestSad   = MAX_UINT;

  // set rcMv (Median predictor) as start point and as best point
  xTZSearchHelp( pcPatternKey, cStruct, rcMv.getHor(), rcMv.getVer(), 0, 0 );

  // test whether one of PRED_A, PRED_B, PRED_C MV is better start point than Median predictor
  if ( bTestOtherPredictedMV )
  {
    for ( UInt index = 0; index < NUM_MV_PREDICTORS; index++ )
    {
      TComMv cMv = m_acMvPredictors[index];
      pcCU->clipMv( cMv );
#if ME_ENABLE_ROUNDING_OF_MVS
      cMv.divideByPowerOf2(2);
#else
      cMv >>= 2;
#endif
      if (cMv != rcMv && (cMv.getHor() != cStruct.iBestX && cMv.getVer() != cStruct.iBestY))
      {
        // only test cMV if not obviously previously tested.
        xTZSearchHelp( pcPatternKey, cStruct, cMv.getHor(), cMv.getVer(), 0, 0 );
      }
    }
  }

  // test whether zero Mv is better start point than Median predictor
  if ( bTestZeroVector )
  {
    if ((rcMv.getHor() != 0 || rcMv.getVer() != 0) &&
        (0 != cStruct.iBestX || 0 != cStruct.iBestY))
    {
      // only test 0-vector if not obviously previously tested.
      xTZSearchHelp( pcPatternKey, cStruct, 0, 0, 0, 0 );
    }
  }

  Int   iSrchRngHorLeft   = pcMvSrchRngLT->getHor();
  Int   iSrchRngHorRight  = pcMvSrchRngRB->getHor();
  Int   iSrchRngVerTop    = pcMvSrchRngLT->getVer();
  Int   iSrchRngVerBottom = pcMvSrchRngRB->getVer();

  if (pIntegerMv2Nx2NPred != 0)
  {
    TComMv integerMv2Nx2NPred = *pIntegerMv2Nx2NPred;
    integerMv2Nx2NPred <<= 2;
    pcCU->clipMv( integerMv2Nx2NPred );
#if ME_ENABLE_ROUNDING_OF_MVS
    integerMv2Nx2NPred.divideByPowerOf2(2);
#else
    integerMv2Nx2NPred >>= 2;
#endif
    if ((rcMv != integerMv2Nx2NPred) &&
        (integerMv2Nx2NPred.getHor() != cStruct.iBestX || integerMv2Nx2NPred.getVer() != cStruct.iBestY))
    {
      // only test integerMv2Nx2NPred if not obviously previously tested.
      xTZSearchHelp(pcPatternKey, cStruct, integerMv2Nx2NPred.getHor(), integerMv2Nx2NPred.getVer(), 0, 0);
    }

    // reset search range
    TComMv cMvSrchRngLT;
    TComMv cMvSrchRngRB;
    Int iSrchRng = m_iSearchRange;
    TComMv currBestMv(cStruct.iBestX, cStruct.iBestY );
    currBestMv <<= 2;
    xSetSearchRange( pcCU, currBestMv, iSrchRng, cMvSrchRngLT, cMvSrchRngRB );
    iSrchRngHorLeft   = cMvSrchRngLT.getHor();
    iSrchRngHorRight  = cMvSrchRngRB.getHor();
    iSrchRngVerTop    = cMvSrchRngLT.getVer();
    iSrchRngVerBottom = cMvSrchRngRB.getVer();
  }

  // start search
  Int  iDist = 0;
  Int  iStartX = cStruct.iBestX;
  Int  iStartY = cStruct.iBestY;

  const Bool bBestCandidateZero = (cStruct.iBestX == 0) && (cStruct.iBestY == 0);

  // first search around best position up to now.
  // The following works as a "subsampled/log" window search around the best candidate
  for (iDist = 1; iDist <= (Int)uiSearchRange; iDist *= 2)
	  
  {
	  flag_start = 1;
	  flag_2point = 0;
	  flag_star = 0;
    if ( bFirstSearchDiamond == 1 )
    {
      xTZ8PointDiamondSearch ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist, bFirstCornersForDiamondDist1 );
    }
    else
    {
      xTZ8PointSquareSearch  ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist );
    }

    if ( bFirstSearchStop && ( cStruct.uiBestRound >= uiFirstSearchRounds ) ) // stop criterion
    {
      break;
    }
  }

  if (!bNewZeroNeighbourhoodTest)
  {
    // test whether zero Mv is a better start point than Median predictor
    if ( bTestZeroVectorStart && ((cStruct.iBestX != 0) || (cStruct.iBestY != 0)) )
    {
      xTZSearchHelp( pcPatternKey, cStruct, 0, 0, 0, 0 );
      if ( (cStruct.iBestX == 0) && (cStruct.iBestY == 0) )
      {
        // test its neighborhood
        for ( iDist = 1; iDist <= (Int)uiSearchRange; iDist*=2 )
        {
          xTZ8PointDiamondSearch( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, 0, 0, iDist, false );
          if ( bTestZeroVectorStop && (cStruct.uiBestRound > 0) ) // stop criterion
          {
            break;
          }
        }
      }
    }
  }
  else
  {
    // Test also zero neighbourhood but with half the range
    // It was reported that the original (above) search scheme using bTestZeroVectorStart did not
    // make sense since one would have already checked the zero candidate earlier
    // and thus the conditions for that test would have not been satisfied
    if (bTestZeroVectorStart == true && bBestCandidateZero != true)
    {
      for ( iDist = 1; iDist <= ((Int)uiSearchRange >> 1); iDist*=2 )
      {
        xTZ8PointDiamondSearch( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, 0, 0, iDist, false );
        if ( bTestZeroVectorStop && (cStruct.uiBestRound > 2) ) // stop criterion
        {
          break;
        }
      }
    }
  }

  // calculate only 2 missing points instead 8 points if cStruct.uiBestDistance == 1
  if ( cStruct.uiBestDistance == 1 )
  {
    cStruct.uiBestDistance = 0;
	flag_2point = 1;
    xTZ2PointSearch( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB );
  }

  // raster search if distance is too big
  if (bUseAdaptiveRaster)
  {
    int iWindowSize = iRaster;
    Int   iSrchRngRasterLeft   = iSrchRngHorLeft;
    Int   iSrchRngRasterRight  = iSrchRngHorRight;
    Int   iSrchRngRasterTop    = iSrchRngVerTop;
    Int   iSrchRngRasterBottom = iSrchRngVerBottom;

    if (!(bEnableRasterSearch && ( ((Int)(cStruct.uiBestDistance) > iRaster))))
    {
      iWindowSize ++;
      iSrchRngRasterLeft /= 2;
      iSrchRngRasterRight /= 2;
      iSrchRngRasterTop /= 2;
      iSrchRngRasterBottom /= 2;
    }
    cStruct.uiBestDistance = iWindowSize;
    for ( iStartY = iSrchRngRasterTop; iStartY <= iSrchRngRasterBottom; iStartY += iWindowSize )
    {
      for ( iStartX = iSrchRngRasterLeft; iStartX <= iSrchRngRasterRight; iStartX += iWindowSize )
      {
        xTZSearchHelp( pcPatternKey, cStruct, iStartX, iStartY, 0, iWindowSize );
      }
    }
  }
  else
  {
    if ( bEnableRasterSearch && ( ((Int)(cStruct.uiBestDistance) > iRaster) || bAlwaysRasterSearch ) )
    {
      cStruct.uiBestDistance = iRaster;
      for ( iStartY = iSrchRngVerTop; iStartY <= iSrchRngVerBottom; iStartY += iRaster )
      {
        for ( iStartX = iSrchRngHorLeft; iStartX <= iSrchRngHorRight; iStartX += iRaster )
        {
          xTZSearchHelp( pcPatternKey, cStruct, iStartX, iStartY, 0, iRaster );
        }
      }
    }
  }

  // raster refinement

  if ( bRasterRefinementEnable && cStruct.uiBestDistance > 0 )
  {
    while ( cStruct.uiBestDistance > 0 )
    {
      iStartX = cStruct.iBestX;
      iStartY = cStruct.iBestY;
      if ( cStruct.uiBestDistance > 1 )
      {
        iDist = cStruct.uiBestDistance >>= 1;
        if ( bRasterRefinementDiamond == 1 )
        {
          xTZ8PointDiamondSearch ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist, bRasterRefinementCornersForDiamondDist1 );
        }
        else
        {
          xTZ8PointSquareSearch  ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist );
        }
      }

      // calculate only 2 missing points instead 8 points if cStruct.uiBestDistance == 1
      if ( cStruct.uiBestDistance == 1 )
      {
        cStruct.uiBestDistance = 0;
        if ( cStruct.ucPointNr != 0 )
        {
          xTZ2PointSearch( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB );
        }
      }
    }
  }

  // star refinement
  if ( bStarRefinementEnable && cStruct.uiBestDistance > 0 )
  {
	  flag_start = 0;
	  flag_2point = 0;
	  flag_star = 1;
    while ( cStruct.uiBestDistance > 0 )
    {
      iStartX = cStruct.iBestX;
      iStartY = cStruct.iBestY;
      cStruct.uiBestDistance = 0;
      cStruct.ucPointNr = 0;
      for ( iDist = 1; iDist < (Int)uiSearchRange + 1; iDist*=2 )
      {
        if ( bStarRefinementDiamond == 1 )
        {
          xTZ8PointDiamondSearch ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist, bStarRefinementCornersForDiamondDist1 );
        }
        else
        {
          xTZ8PointSquareSearch  ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist );
        }
        if ( bStarRefinementStop && (cStruct.uiBestRound >= uiStarRefinementRounds) ) // stop criterion
        {
          break;
        }
      }

      // calculate only 2 missing points instead 8 points if cStrukt.uiBestDistance == 1
      if ( cStruct.uiBestDistance == 1 )
      {
        cStruct.uiBestDistance = 0;
        if ( cStruct.ucPointNr != 0 )
        {
			flag_2point =1;
          xTZ2PointSearch( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB );
        }
      }
    }
  }

  // getting the 8 SAD points
  iDist = 1;
  iStartX = cStruct.iBestX;
  iStartY = cStruct.iBestY;
  index_ref = counter_i;
  SAD_Best=cStruct.uiBestSad;
  MVX_Best=cStruct.iBestX;
  MVY_Best=cStruct.iBestY;
  
  xTZ8PointSquareSearch(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);

  iDist = 2;
  xTZ8PointSquareSearch2(pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist);
  
  cStruct.iBestX=MVX_Best;
  cStruct.iBestY=MVY_Best;
  cStruct.uiBestSad=SAD_Best;

  // write out best match
  rcMv.set( cStruct.iBestX, cStruct.iBestY );
  ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCostOfVectorWithPredictor( cStruct.iBestX, cStruct.iBestY );
}


Void TEncSearch::xTZSearchSelective( const TComDataCU* const   pcCU,
                                     const TComPattern* const  pcPatternKey,
                                     const Pel* const          piRefY,
                                     const Int                 iRefStride,
                                     const TComMv* const       pcMvSrchRngLT,
                                     const TComMv* const       pcMvSrchRngRB,
                                     TComMv                   &rcMv,
                                     Distortion               &ruiSAD,
                                     const TComMv* const       pIntegerMv2Nx2NPred )
{
  const Bool bTestOtherPredictedMV    = true;
  const Bool bTestZeroVector          = true;
  const Bool bEnableRasterSearch      = true;
  const Bool bAlwaysRasterSearch      = false;  // 1: BETTER but factor 15x slower
  const Bool bStarRefinementEnable    = true;   // enable either star refinement or raster refinement
  const Bool bStarRefinementDiamond   = true;   // 1 = xTZ8PointDiamondSearch   0 = xTZ8PointSquareSearch
  const Bool bStarRefinementStop      = false;
  const UInt uiStarRefinementRounds   = 2;  // star refinement stop X rounds after best match (must be >=1)
  const UInt uiSearchRange            = m_iSearchRange;
  const Int  uiSearchRangeInitial     = m_iSearchRange >> 2;
  const Int  uiSearchStep             = 4;
  const Int  iMVDistThresh            = 8;

  Int   iSrchRngHorLeft         = pcMvSrchRngLT->getHor();
  Int   iSrchRngHorRight        = pcMvSrchRngRB->getHor();
  Int   iSrchRngVerTop          = pcMvSrchRngLT->getVer();
  Int   iSrchRngVerBottom       = pcMvSrchRngRB->getVer();
  Int   iFirstSrchRngHorLeft    = 0;
  Int   iFirstSrchRngHorRight   = 0;
  Int   iFirstSrchRngVerTop     = 0;
  Int   iFirstSrchRngVerBottom  = 0;
  Int   iStartX                 = 0;
  Int   iStartY                 = 0;
  Int   iBestX                  = 0;
  Int   iBestY                  = 0;
  Int   iDist                   = 0;

  pcCU->clipMv( rcMv );
#if ME_ENABLE_ROUNDING_OF_MVS
  rcMv.divideByPowerOf2(2);
#else
  rcMv >>= 2;
#endif
  // init TZSearchStruct
  IntTZSearchStruct cStruct;
  cStruct.iYStride    = iRefStride;
  cStruct.piRefY      = piRefY;
  cStruct.uiBestSad   = MAX_UINT;
  cStruct.iBestX = 0;
  cStruct.iBestY = 0;


  // set rcMv (Median predictor) as start point and as best point
  xTZSearchHelp( pcPatternKey, cStruct, rcMv.getHor(), rcMv.getVer(), 0, 0 );

  // test whether one of PRED_A, PRED_B, PRED_C MV is better start point than Median predictor
  if ( bTestOtherPredictedMV )
  {
    for ( UInt index = 0; index < NUM_MV_PREDICTORS; index++ )
    {
      TComMv cMv = m_acMvPredictors[index];
      pcCU->clipMv( cMv );
#if ME_ENABLE_ROUNDING_OF_MVS
      cMv.divideByPowerOf2(2);
#else
      cMv >>= 2;
#endif
      xTZSearchHelp( pcPatternKey, cStruct, cMv.getHor(), cMv.getVer(), 0, 0 );
    }
  }

  // test whether zero Mv is better start point than Median predictor
  if ( bTestZeroVector )
  {
    xTZSearchHelp( pcPatternKey, cStruct, 0, 0, 0, 0 );
  }

  if ( pIntegerMv2Nx2NPred != 0 )
  {
    TComMv integerMv2Nx2NPred = *pIntegerMv2Nx2NPred;
    integerMv2Nx2NPred <<= 2;
    pcCU->clipMv( integerMv2Nx2NPred );
#if ME_ENABLE_ROUNDING_OF_MVS
    integerMv2Nx2NPred.divideByPowerOf2(2);
#else
    integerMv2Nx2NPred >>= 2;
#endif
    xTZSearchHelp(pcPatternKey, cStruct, integerMv2Nx2NPred.getHor(), integerMv2Nx2NPred.getVer(), 0, 0);

    // reset search range
    TComMv cMvSrchRngLT;
    TComMv cMvSrchRngRB;
    Int iSrchRng = m_iSearchRange;
    TComMv currBestMv(cStruct.iBestX, cStruct.iBestY );
    currBestMv <<= 2;
    xSetSearchRange( pcCU, currBestMv, iSrchRng, cMvSrchRngLT, cMvSrchRngRB );
    iSrchRngHorLeft   = cMvSrchRngLT.getHor();
    iSrchRngHorRight  = cMvSrchRngRB.getHor();
    iSrchRngVerTop    = cMvSrchRngLT.getVer();
    iSrchRngVerBottom = cMvSrchRngRB.getVer();
  }

  // Initial search
  iBestX = cStruct.iBestX;
  iBestY = cStruct.iBestY; 
  iFirstSrchRngHorLeft    = ((iBestX - uiSearchRangeInitial) > iSrchRngHorLeft)   ? (iBestX - uiSearchRangeInitial) : iSrchRngHorLeft;
  iFirstSrchRngVerTop     = ((iBestY - uiSearchRangeInitial) > iSrchRngVerTop)    ? (iBestY - uiSearchRangeInitial) : iSrchRngVerTop;
  iFirstSrchRngHorRight   = ((iBestX + uiSearchRangeInitial) < iSrchRngHorRight)  ? (iBestX + uiSearchRangeInitial) : iSrchRngHorRight;  
  iFirstSrchRngVerBottom  = ((iBestY + uiSearchRangeInitial) < iSrchRngVerBottom) ? (iBestY + uiSearchRangeInitial) : iSrchRngVerBottom;    

  for ( iStartY = iFirstSrchRngVerTop; iStartY <= iFirstSrchRngVerBottom; iStartY += uiSearchStep )
  {
    for ( iStartX = iFirstSrchRngHorLeft; iStartX <= iFirstSrchRngHorRight; iStartX += uiSearchStep )
    {
      xTZSearchHelp( pcPatternKey, cStruct, iStartX, iStartY, 0, 0 );
      xTZ8PointDiamondSearch ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, 1, false );
      xTZ8PointDiamondSearch ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, 2, false );
    }
  }

  Int iMaxMVDistToPred = (abs(cStruct.iBestX - iBestX) > iMVDistThresh || abs(cStruct.iBestY - iBestY) > iMVDistThresh);

  //full search with early exit if MV is distant from predictors
  if ( bEnableRasterSearch && (iMaxMVDistToPred || bAlwaysRasterSearch) )
  {
    for ( iStartY = iSrchRngVerTop; iStartY <= iSrchRngVerBottom; iStartY += 1 )
    {
      for ( iStartX = iSrchRngHorLeft; iStartX <= iSrchRngHorRight; iStartX += 1 )
      {
        xTZSearchHelp( pcPatternKey, cStruct, iStartX, iStartY, 0, 1 );
      }
    }
  }
  //Smaller MV, refine around predictor
  else if ( bStarRefinementEnable && cStruct.uiBestDistance > 0 )
  {
    // start refinement
    while ( cStruct.uiBestDistance > 0 )
    {
      iStartX = cStruct.iBestX;
      iStartY = cStruct.iBestY;
      cStruct.uiBestDistance = 0;
      cStruct.ucPointNr = 0;
      for ( iDist = 1; iDist < (Int)uiSearchRange + 1; iDist*=2 )
      {
        if ( bStarRefinementDiamond == 1 )
        {
          xTZ8PointDiamondSearch ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist, false );
        }
        else
        {
          xTZ8PointSquareSearch  ( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB, iStartX, iStartY, iDist );
        }
        if ( bStarRefinementStop && (cStruct.uiBestRound >= uiStarRefinementRounds) ) // stop criterion
        {
          break;
        }
      }

      // calculate only 2 missing points instead 8 points if cStrukt.uiBestDistance == 1
      if ( cStruct.uiBestDistance == 1 )
      {
        cStruct.uiBestDistance = 0;
        if ( cStruct.ucPointNr != 0 )
        {
          xTZ2PointSearch( pcPatternKey, cStruct, pcMvSrchRngLT, pcMvSrchRngRB );
        }
      }
    }
  }

  // write out best match
  rcMv.set( cStruct.iBestX, cStruct.iBestY );
  ruiSAD = cStruct.uiBestSad - m_pcRdCost->getCostOfVectorWithPredictor( cStruct.iBestX, cStruct.iBestY );

}


Void TEncSearch::xPatternSearchFracDIF(
                                       Bool         bIsLosslessCoded,
                                       TComPattern* pcPatternKey,
                                       Pel*         piRefY,
                                       Int          iRefStride,
                                       TComMv*      pcMvInt,
                                       TComMv&      rcMvHalf,
                                       TComMv&      rcMvQter,
                                       Distortion&  ruiCost
                                      )
{
  //  Reference pattern initialization (integer scale)
	
  TComPattern cPatternRoi;
  Int         iOffset    = pcMvInt->getHor() + pcMvInt->getVer() * iRefStride;
  cPatternRoi.initPattern(piRefY + iOffset,
                          pcPatternKey->getROIYWidth(),
                          pcPatternKey->getROIYHeight(),
                          iRefStride,
                          pcPatternKey->getBitDepthY());

  //  Half-pel refinement
  xExtDIFUpSamplingH ( &cPatternRoi );

  rcMvHalf = *pcMvInt;   rcMvHalf <<= 1;    // for mv-cost
  TComMv baseRefMv(0, 0);
  ruiCost = xPatternRefinement( pcPatternKey, baseRefMv, 2, rcMvHalf, !bIsLosslessCoded );

  m_pcRdCost->setCostScale( 0 );

  xExtDIFUpSamplingQ ( &cPatternRoi, rcMvHalf );
  baseRefMv = rcMvHalf;
  baseRefMv <<= 1;

  rcMvQter = *pcMvInt;   rcMvQter <<= 1;    // for mv-cost
  rcMvQter += rcMvHalf;  rcMvQter <<= 1;
  ruiCost = xPatternRefinement( pcPatternKey, baseRefMv, 1, rcMvQter, !bIsLosslessCoded );
}


//! encode residual and calculate rate-distortion for a CU block
Void TEncSearch::encodeResAndCalcRdInterCU( TComDataCU* pcCU, TComYuv* pcYuvOrg, TComYuv* pcYuvPred,
                                            TComYuv* pcYuvResi, TComYuv* pcYuvResiBest, TComYuv* pcYuvRec,
                                            Bool bSkipResidual DEBUG_STRING_FN_DECLARE(sDebug) )
{
  assert ( !pcCU->isIntra(0) );

  const UInt cuWidthPixels      = pcCU->getWidth ( 0 );
  const UInt cuHeightPixels     = pcCU->getHeight( 0 );
  const Int  numValidComponents = pcCU->getPic()->getNumberValidComponents();
  const TComSPS &sps=*(pcCU->getSlice()->getSPS());

  // The pcCU is not marked as skip-mode at this point, and its m_pcTrCoeff, m_pcArlCoeff, m_puhCbf, m_puhTrIdx will all be 0.
  // due to prior calls to TComDataCU::initEstData(  );

  if ( bSkipResidual ) //  No residual coding : SKIP mode
  {
    pcCU->setSkipFlagSubParts( true, 0, pcCU->getDepth(0) );

    pcYuvResi->clear();

    pcYuvPred->copyToPartYuv( pcYuvRec, 0 );
    Distortion distortion = 0;

    for (Int comp=0; comp < numValidComponents; comp++)
    {
      const ComponentID compID=ComponentID(comp);
      const UInt csx=pcYuvOrg->getComponentScaleX(compID);
      const UInt csy=pcYuvOrg->getComponentScaleY(compID);
      distortion += m_pcRdCost->getDistPart( sps.getBitDepth(toChannelType(compID)), pcYuvRec->getAddr(compID), pcYuvRec->getStride(compID), pcYuvOrg->getAddr(compID),
                                               pcYuvOrg->getStride(compID), cuWidthPixels >> csx, cuHeightPixels >> csy, compID);
    }

    m_pcRDGoOnSbacCoder->load(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_CURR_BEST]);
    m_pcEntropyCoder->resetBits();

    if (pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
    {
      m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, 0, true);
    }

    m_pcEntropyCoder->encodeSkipFlag(pcCU, 0, true);
    m_pcEntropyCoder->encodeMergeIndex( pcCU, 0, true );

    UInt uiBits = m_pcEntropyCoder->getNumberOfWrittenBits();
    pcCU->getTotalBits()       = uiBits;
    pcCU->getTotalDistortion() = distortion;
    pcCU->getTotalCost()       = m_pcRdCost->calcRdCost( uiBits, distortion );

    m_pcRDGoOnSbacCoder->store(m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_TEMP_BEST]);

#if DEBUG_STRING
    pcYuvResiBest->clear(); // Clear the residual image, if we didn't code it.
    for(UInt i=0; i<MAX_NUM_COMPONENT+1; i++)
    {
      sDebug+=debug_reorder_data_inter_token[i];
    }
#endif

    return;
  }

  //  Residual coding.

   pcYuvResi->subtract( pcYuvOrg, pcYuvPred, 0, cuWidthPixels );

  TComTURecurse tuLevel0(pcCU, 0);

  Double     nonZeroCost       = 0;
  UInt       nonZeroBits       = 0;
  Distortion nonZeroDistortion = 0;
  Distortion zeroDistortion    = 0;

  m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ pcCU->getDepth( 0 ) ][ CI_CURR_BEST ] );

  xEstimateInterResidualQT( pcYuvResi,  nonZeroCost, nonZeroBits, nonZeroDistortion, &zeroDistortion, tuLevel0 DEBUG_STRING_PASS_INTO(sDebug) );

  // -------------------------------------------------------
  // set the coefficients in the pcCU, and also calculates the residual data.
  // If a block full of 0's is efficient, then just use 0's.
  // The costs at this point do not include header bits.

  m_pcEntropyCoder->resetBits();
  m_pcEntropyCoder->encodeQtRootCbfZero( );
  const UInt   zeroResiBits = m_pcEntropyCoder->getNumberOfWrittenBits();
  const Double zeroCost     = (pcCU->isLosslessCoded( 0 )) ? (nonZeroCost+1) : (m_pcRdCost->calcRdCost( zeroResiBits, zeroDistortion ));

  if ( zeroCost < nonZeroCost || !pcCU->getQtRootCbf(0) )
  {
    const UInt uiQPartNum = tuLevel0.GetAbsPartIdxNumParts();
    ::memset( pcCU->getTransformIdx()     , 0, uiQPartNum * sizeof(UChar) );
    for (Int comp=0; comp < numValidComponents; comp++)
    {
      const ComponentID component = ComponentID(comp);
      ::memset( pcCU->getCbf( component ) , 0, uiQPartNum * sizeof(UChar) );
      ::memset( pcCU->getCrossComponentPredictionAlpha(component), 0, ( uiQPartNum * sizeof(SChar) ) );
    }
    static const UInt useTS[MAX_NUM_COMPONENT]={0,0,0};
    pcCU->setTransformSkipSubParts ( useTS, 0, pcCU->getDepth(0) );
#if DEBUG_STRING
    sDebug.clear();
    for(UInt i=0; i<MAX_NUM_COMPONENT+1; i++)
    {
      sDebug+=debug_reorder_data_inter_token[i];
    }
#endif
  }
  else
  {
    xSetInterResidualQTData( NULL, false, tuLevel0); // Call first time to set coefficients.
  }

  // all decisions now made. Fully encode the CU, including the headers:
  m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[pcCU->getDepth(0)][CI_CURR_BEST] );

  UInt finalBits = 0;
  xAddSymbolBitsInter( pcCU, finalBits );
  // we've now encoded the pcCU, and so have a valid bit cost

  if ( !pcCU->getQtRootCbf( 0 ) )
  {
    pcYuvResiBest->clear(); // Clear the residual image, if we didn't code it.
  }
  else
  {
    xSetInterResidualQTData( pcYuvResiBest, true, tuLevel0 ); // else set the residual image data pcYUVResiBest from the various temp images.
  }
  m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ pcCU->getDepth( 0 ) ][ CI_TEMP_BEST ] );

  pcYuvRec->addClip ( pcYuvPred, pcYuvResiBest, 0, cuWidthPixels, sps.getBitDepths() );

  // update with clipped distortion and cost (previously unclipped reconstruction values were used)

  Distortion finalDistortion = 0;
  for(Int comp=0; comp<numValidComponents; comp++)
  {
    const ComponentID compID=ComponentID(comp);
    finalDistortion += m_pcRdCost->getDistPart( sps.getBitDepth(toChannelType(compID)), pcYuvRec->getAddr(compID ), pcYuvRec->getStride(compID ), pcYuvOrg->getAddr(compID ), pcYuvOrg->getStride(compID), cuWidthPixels >> pcYuvOrg->getComponentScaleX(compID), cuHeightPixels >> pcYuvOrg->getComponentScaleY(compID), compID);
  }

  pcCU->getTotalBits()       = finalBits;
  pcCU->getTotalDistortion() = finalDistortion;
  pcCU->getTotalCost()       = m_pcRdCost->calcRdCost( finalBits, finalDistortion );
}



Void TEncSearch::xEstimateInterResidualQT( TComYuv    *pcResi,
                                           Double     &rdCost,
                                           UInt       &ruiBits,
                                           Distortion &ruiDist,
                                           Distortion *puiZeroDist,
                                           TComTU     &rTu
                                           DEBUG_STRING_FN_DECLARE(sDebug) )
{
  TComDataCU *pcCU        = rTu.getCU();
  const UInt uiAbsPartIdx = rTu.GetAbsPartIdxTU();
  const UInt uiDepth      = rTu.GetTransformDepthTotal();
  const UInt uiTrMode     = rTu.GetTransformDepthRel();
  const UInt subTUDepth   = uiTrMode + 1;
  const UInt numValidComp = pcCU->getPic()->getNumberValidComponents();
  DEBUG_STRING_NEW(sSingleStringComp[MAX_NUM_COMPONENT])

  assert( pcCU->getDepth( 0 ) == pcCU->getDepth( uiAbsPartIdx ) );
  const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();

  UInt SplitFlag = ((pcCU->getSlice()->getSPS()->getQuadtreeTUMaxDepthInter() == 1) && pcCU->isInter(uiAbsPartIdx) && ( pcCU->getPartitionSize(uiAbsPartIdx) != SIZE_2Nx2N ));
#if DEBUG_STRING
  const Int debugPredModeMask = DebugStringGetPredModeMask(pcCU->getPredictionMode(uiAbsPartIdx));
#endif

  Bool bCheckFull;

  if ( SplitFlag && uiDepth == pcCU->getDepth(uiAbsPartIdx) && ( uiLog2TrSize >  pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) ) )
  {
    bCheckFull = false;
  }
  else
  {
    bCheckFull =  ( uiLog2TrSize <= pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() );
  }

  const Bool bCheckSplit  = ( uiLog2TrSize >  pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) );

  assert( bCheckFull || bCheckSplit );

  // code full block
  Double     dSingleCost = MAX_DOUBLE;
  UInt       uiSingleBits                                                                                                        = 0;
  Distortion uiSingleDistComp            [MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = {{0,0},{0,0},{0,0}};
  Distortion uiSingleDist                                                                                                        = 0;
  TCoeff     uiAbsSum                    [MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = {{0,0},{0,0},{0,0}};
  UInt       uiBestTransformMode         [MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = {{0,0},{0,0},{0,0}};
  //  Stores the best explicit RDPCM mode for a TU encoded without split
  UInt       bestExplicitRdpcmModeUnSplit[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = {{3,3}, {3,3}, {3,3}};
  SChar      bestCrossCPredictionAlpha   [MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/] = {{0,0},{0,0},{0,0}};

  m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_ROOT ] );

  if( bCheckFull )
  {
    Double minCost[MAX_NUM_COMPONENT][2/*0 = top (or whole TU for non-4:2:2) sub-TU, 1 = bottom sub-TU*/];
    Bool checkTransformSkip[MAX_NUM_COMPONENT];
    pcCU->setTrIdxSubParts( uiTrMode, uiAbsPartIdx, uiDepth );

    m_pcEntropyCoder->resetBits();

    memset( m_pTempPel, 0, sizeof( Pel ) * rTu.getRect(COMPONENT_Y).width * rTu.getRect(COMPONENT_Y).height ); // not necessary needed for inside of recursion (only at the beginning)

    const UInt uiQTTempAccessLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
    TCoeff *pcCoeffCurr[MAX_NUM_COMPONENT];
#if ADAPTIVE_QP_SELECTION
    TCoeff *pcArlCoeffCurr[MAX_NUM_COMPONENT];
#endif

    for(UInt i=0; i<numValidComp; i++)
    {
      minCost[i][0] = MAX_DOUBLE;
      minCost[i][1] = MAX_DOUBLE;
    }

    Pel crossCPredictedResidualBuffer[ MAX_TU_SIZE * MAX_TU_SIZE ];

    for(UInt i=0; i<numValidComp; i++)
    {
      checkTransformSkip[i]=false;
      const ComponentID compID=ComponentID(i);
      const Int channelBitDepth=pcCU->getSlice()->getSPS()->getBitDepth(toChannelType(compID));
      pcCoeffCurr[compID]    = m_ppcQTTempCoeff[compID][uiQTTempAccessLayer] + rTu.getCoefficientOffset(compID);
#if ADAPTIVE_QP_SELECTION
      pcArlCoeffCurr[compID] = m_ppcQTTempArlCoeff[compID ][uiQTTempAccessLayer] +  rTu.getCoefficientOffset(compID);
#endif

      if(rTu.ProcessComponentSection(compID))
      {
        const QpParam cQP(*pcCU, compID);

        checkTransformSkip[compID] = pcCU->getSlice()->getPPS()->getUseTransformSkip() &&
                                     TUCompRectHasAssociatedTransformSkipFlag(rTu.getRect(compID), pcCU->getSlice()->getPPS()->getPpsRangeExtension().getLog2MaxTransformSkipBlockSize()) &&
                                     (!pcCU->isLosslessCoded(0));

        const Bool splitIntoSubTUs = rTu.getRect(compID).width != rTu.getRect(compID).height;

        TComTURecurse TUIterator(rTu, false, (splitIntoSubTUs ? TComTU::VERTICAL_SPLIT : TComTU::DONT_SPLIT), true, compID);

        const UInt partIdxesPerSubTU = TUIterator.GetAbsPartIdxNumParts(compID);

        do
        {
          const UInt           subTUIndex             = TUIterator.GetSectionNumber();
          const UInt           subTUAbsPartIdx        = TUIterator.GetAbsPartIdxTU(compID);
          const TComRectangle &tuCompRect             = TUIterator.getRect(compID);
          const UInt           subTUBufferOffset      = tuCompRect.width * tuCompRect.height * subTUIndex;

                TCoeff        *currentCoefficients    = pcCoeffCurr[compID] + subTUBufferOffset;
#if ADAPTIVE_QP_SELECTION
                TCoeff        *currentARLCoefficients = pcArlCoeffCurr[compID] + subTUBufferOffset;
#endif
          const Bool isCrossCPredictionAvailable      =    isChroma(compID)
                                                         && pcCU->getSlice()->getPPS()->getPpsRangeExtension().getCrossComponentPredictionEnabledFlag()
                                                         && (pcCU->getCbf(subTUAbsPartIdx, COMPONENT_Y, uiTrMode) != 0);

          SChar preCalcAlpha = 0;
          const Pel *pLumaResi = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix( COMPONENT_Y, rTu.getRect( COMPONENT_Y ).x0, rTu.getRect( COMPONENT_Y ).y0 );

          if (isCrossCPredictionAvailable)
          {
            const Bool bUseReconstructedResidualForEstimate = m_pcEncCfg->getUseReconBasedCrossCPredictionEstimate();
            const Pel  *const lumaResidualForEstimate       = bUseReconstructedResidualForEstimate ? pLumaResi                                                     : pcResi->getAddrPix(COMPONENT_Y, tuCompRect.x0, tuCompRect.y0);
            const UInt        lumaResidualStrideForEstimate = bUseReconstructedResidualForEstimate ? m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y) : pcResi->getStride(COMPONENT_Y);

            preCalcAlpha = xCalcCrossComponentPredictionAlpha(TUIterator,
                                                              compID,
                                                              lumaResidualForEstimate,
                                                              pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
                                                              tuCompRect.width,
                                                              tuCompRect.height,
                                                              lumaResidualStrideForEstimate,
                                                              pcResi->getStride(compID));
          }

          const Int transformSkipModesToTest    = checkTransformSkip[compID] ? 2 : 1;
          const Int crossCPredictionModesToTest = (preCalcAlpha != 0)        ? 2 : 1; // preCalcAlpha cannot be anything other than 0 if isCrossCPredictionAvailable is false

          const Bool isOneMode                  = (crossCPredictionModesToTest == 1) && (transformSkipModesToTest == 1);

          for (Int transformSkipModeId = 0; transformSkipModeId < transformSkipModesToTest; transformSkipModeId++)
          {
            pcCU->setTransformSkipPartRange(transformSkipModeId, compID, subTUAbsPartIdx, partIdxesPerSubTU);

            for (Int crossCPredictionModeId = 0; crossCPredictionModeId < crossCPredictionModesToTest; crossCPredictionModeId++)
            {
              const Bool isFirstMode          = (transformSkipModeId == 0) && (crossCPredictionModeId == 0);
              const Bool bUseCrossCPrediction = crossCPredictionModeId != 0;

              m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_ROOT ] );
              m_pcEntropyCoder->resetBits();

              pcCU->setTransformSkipPartRange(transformSkipModeId, compID, subTUAbsPartIdx, partIdxesPerSubTU);
              pcCU->setCrossComponentPredictionAlphaPartRange((bUseCrossCPrediction ? preCalcAlpha : 0), compID, subTUAbsPartIdx, partIdxesPerSubTU );

              if ((compID != COMPONENT_Cr) && ((transformSkipModeId == 1) ? m_pcEncCfg->getUseRDOQTS() : m_pcEncCfg->getUseRDOQ()))
              {
                m_pcEntropyCoder->estimateBit(m_pcTrQuant->m_pcEstBitsSbac, tuCompRect.width, tuCompRect.height, toChannelType(compID));
              }

#if RDOQ_CHROMA_LAMBDA
              m_pcTrQuant->selectLambda(compID);
#endif

              Pel *pcResiCurrComp = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
              UInt resiStride     = m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID);

              TCoeff bestCoeffComp   [MAX_TU_SIZE*MAX_TU_SIZE];
              Pel    bestResiComp    [MAX_TU_SIZE*MAX_TU_SIZE];

#if ADAPTIVE_QP_SELECTION
              TCoeff bestArlCoeffComp[MAX_TU_SIZE*MAX_TU_SIZE];
#endif
              TCoeff     currAbsSum   = 0;
              UInt       currCompBits = 0;
              Distortion currCompDist = 0;
              Double     currCompCost = 0;
              UInt       nonCoeffBits = 0;
              Distortion nonCoeffDist = 0;
              Double     nonCoeffCost = 0;

              if(!isOneMode && !isFirstMode)
              {
                memcpy(bestCoeffComp,    currentCoefficients,    (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#if ADAPTIVE_QP_SELECTION
                memcpy(bestArlCoeffComp, currentARLCoefficients, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#endif
                for(Int y = 0; y < tuCompRect.height; y++)
                {
                  memcpy(&bestResiComp[y * tuCompRect.width], (pcResiCurrComp + (y * resiStride)), (sizeof(Pel) * tuCompRect.width));
                }
              }

              if (bUseCrossCPrediction)
              {
                TComTrQuant::crossComponentPrediction(TUIterator,
                                                      compID,
                                                      pLumaResi,
                                                      pcResi->getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
                                                      crossCPredictedResidualBuffer,
                                                      tuCompRect.width,
                                                      tuCompRect.height,
                                                      m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
                                                      pcResi->getStride(compID),
                                                      tuCompRect.width,
                                                      false);

                m_pcTrQuant->transformNxN(TUIterator, compID, crossCPredictedResidualBuffer, tuCompRect.width, currentCoefficients,
#if ADAPTIVE_QP_SELECTION
                                          currentARLCoefficients,
#endif
                                          currAbsSum, cQP);
              }
              else
              {
                m_pcTrQuant->transformNxN(TUIterator, compID, pcResi->getAddrPix( compID, tuCompRect.x0, tuCompRect.y0 ), pcResi->getStride(compID), currentCoefficients,
#if ADAPTIVE_QP_SELECTION
                                          currentARLCoefficients,
#endif
                                          currAbsSum, cQP);
              }

              if(isFirstMode || (currAbsSum == 0))
              {
                if (bUseCrossCPrediction)
                {
                  TComTrQuant::crossComponentPrediction(TUIterator,
                                                        compID,
                                                        pLumaResi,
                                                        m_pTempPel,
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
                                                        tuCompRect.width,
                                                        tuCompRect.height,
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
                                                        tuCompRect.width,
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
                                                        true);

                  nonCoeffDist = m_pcRdCost->getDistPart( channelBitDepth, m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix( compID, tuCompRect.x0, tuCompRect.y0 ),
                                                          m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride( compID ), pcResi->getAddrPix( compID, tuCompRect.x0, tuCompRect.y0 ),
                                                          pcResi->getStride(compID), tuCompRect.width, tuCompRect.height, compID); // initialized with zero residual distortion
                }
                else
                {
                  nonCoeffDist = m_pcRdCost->getDistPart( channelBitDepth, m_pTempPel, tuCompRect.width, pcResi->getAddrPix( compID, tuCompRect.x0, tuCompRect.y0 ),
                                                          pcResi->getStride(compID), tuCompRect.width, tuCompRect.height, compID); // initialized with zero residual distortion
                }

                m_pcEntropyCoder->encodeQtCbfZero( TUIterator, toChannelType(compID) );

                if ( isCrossCPredictionAvailable )
                {
                  m_pcEntropyCoder->encodeCrossComponentPrediction( TUIterator, compID );
                }

                nonCoeffBits = m_pcEntropyCoder->getNumberOfWrittenBits();
                nonCoeffCost = m_pcRdCost->calcRdCost( nonCoeffBits, nonCoeffDist );
              }

              if((puiZeroDist != NULL) && isFirstMode)
              {
                *puiZeroDist += nonCoeffDist; // initialized with zero residual distortion
              }

              DEBUG_STRING_NEW(sSingleStringTest)

              if( currAbsSum > 0 ) //if non-zero coefficients are present, a residual needs to be derived for further prediction
              {
                if (isFirstMode)
                {
                  m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_ROOT ] );
                  m_pcEntropyCoder->resetBits();
                }

                m_pcEntropyCoder->encodeQtCbf( TUIterator, compID, true );

                if (isCrossCPredictionAvailable)
                {
                  m_pcEntropyCoder->encodeCrossComponentPrediction( TUIterator, compID );
                }

                m_pcEntropyCoder->encodeCoeffNxN( TUIterator, currentCoefficients, compID );
                currCompBits = m_pcEntropyCoder->getNumberOfWrittenBits();

                pcResiCurrComp = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix( compID, tuCompRect.x0, tuCompRect.y0 );

                m_pcTrQuant->invTransformNxN( TUIterator, compID, pcResiCurrComp, m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID), currentCoefficients, cQP DEBUG_STRING_PASS_INTO_OPTIONAL(&sSingleStringTest, (DebugOptionList::DebugString_InvTran.getInt()&debugPredModeMask)) );

                if (bUseCrossCPrediction)
                {
                  TComTrQuant::crossComponentPrediction(TUIterator,
                                                        compID,
                                                        pLumaResi,
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
                                                        tuCompRect.width,
                                                        tuCompRect.height,
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID     ),
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID     ),
                                                        true);
                }

                currCompDist = m_pcRdCost->getDistPart( channelBitDepth, m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix( compID, tuCompRect.x0, tuCompRect.y0 ),
                                                        m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
                                                        pcResi->getAddrPix( compID, tuCompRect.x0, tuCompRect.y0 ),
                                                        pcResi->getStride(compID),
                                                        tuCompRect.width, tuCompRect.height, compID);

                currCompCost = m_pcRdCost->calcRdCost(currCompBits, currCompDist);
                  
                if (pcCU->isLosslessCoded(0))
                {
                  nonCoeffCost = MAX_DOUBLE;
                }
              }
              else if ((transformSkipModeId == 1) && !bUseCrossCPrediction)
              {
                currCompCost = MAX_DOUBLE;
              }
              else
              {
                currCompBits = nonCoeffBits;
                currCompDist = nonCoeffDist;
                currCompCost = nonCoeffCost;
              }

              // evaluate
              if ((currCompCost < minCost[compID][subTUIndex]) || ((transformSkipModeId == 1) && (currCompCost == minCost[compID][subTUIndex])))
              {
                bestExplicitRdpcmModeUnSplit[compID][subTUIndex] = pcCU->getExplicitRdpcmMode(compID, subTUAbsPartIdx);

                if(isFirstMode) //check for forced null
                {
                  if((nonCoeffCost < currCompCost) || (currAbsSum == 0))
                  {
                    memset(currentCoefficients, 0, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));

                    currAbsSum   = 0;
                    currCompBits = nonCoeffBits;
                    currCompDist = nonCoeffDist;
                    currCompCost = nonCoeffCost;
                  }
                }

#if DEBUG_STRING
                if (currAbsSum > 0)
                {
                  DEBUG_STRING_SWAP(sSingleStringComp[compID], sSingleStringTest)
                }
                else
                {
                  sSingleStringComp[compID].clear();
                }
#endif

                uiAbsSum                 [compID][subTUIndex] = currAbsSum;
                uiSingleDistComp         [compID][subTUIndex] = currCompDist;
                minCost                  [compID][subTUIndex] = currCompCost;
                uiBestTransformMode      [compID][subTUIndex] = transformSkipModeId;
                bestCrossCPredictionAlpha[compID][subTUIndex] = (crossCPredictionModeId == 1) ? pcCU->getCrossComponentPredictionAlpha(subTUAbsPartIdx, compID) : 0;

                if (uiAbsSum[compID][subTUIndex] == 0)
                {
                  if (bUseCrossCPrediction)
                  {
                    TComTrQuant::crossComponentPrediction(TUIterator,
                                                          compID,
                                                          pLumaResi,
                                                          m_pTempPel,
                                                          m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0),
                                                          tuCompRect.width,
                                                          tuCompRect.height,
                                                          m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(COMPONENT_Y),
                                                          tuCompRect.width,
                                                          m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID),
                                                          true);
                  }
                  else
                  {
                    pcResiCurrComp = m_pcQTTempTComYuv[uiQTTempAccessLayer].getAddrPix(compID, tuCompRect.x0, tuCompRect.y0);
                    const UInt uiStride = m_pcQTTempTComYuv[uiQTTempAccessLayer].getStride(compID);
                    for(UInt uiY = 0; uiY < tuCompRect.height; uiY++)
                    {
                      memset(pcResiCurrComp, 0, (sizeof(Pel) * tuCompRect.width));
                      pcResiCurrComp += uiStride;
                    }
                  }
                }
              }
              else
              {
                // reset
                memcpy(currentCoefficients,    bestCoeffComp,    (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#if ADAPTIVE_QP_SELECTION
                memcpy(currentARLCoefficients, bestArlCoeffComp, (sizeof(TCoeff) * tuCompRect.width * tuCompRect.height));
#endif
                for (Int y = 0; y < tuCompRect.height; y++)
                {
                  memcpy((pcResiCurrComp + (y * resiStride)), &bestResiComp[y * tuCompRect.width], (sizeof(Pel) * tuCompRect.width));
                }
              }
            }
          }

          pcCU->setExplicitRdpcmModePartRange            (   bestExplicitRdpcmModeUnSplit[compID][subTUIndex],                            compID, subTUAbsPartIdx, partIdxesPerSubTU);
          pcCU->setTransformSkipPartRange                (   uiBestTransformMode         [compID][subTUIndex],                            compID, subTUAbsPartIdx, partIdxesPerSubTU );
          pcCU->setCbfPartRange                          ((((uiAbsSum                    [compID][subTUIndex] > 0) ? 1 : 0) << uiTrMode), compID, subTUAbsPartIdx, partIdxesPerSubTU );
          pcCU->setCrossComponentPredictionAlphaPartRange(   bestCrossCPredictionAlpha   [compID][subTUIndex],                            compID, subTUAbsPartIdx, partIdxesPerSubTU );
        } while (TUIterator.nextSection(rTu)); //end of sub-TU loop
      } // processing section
    } // component loop

    for(UInt ch = 0; ch < numValidComp; ch++)
    {
      const ComponentID compID = ComponentID(ch);
      if (rTu.ProcessComponentSection(compID) && (rTu.getRect(compID).width != rTu.getRect(compID).height))
      {
        offsetSubTUCBFs(rTu, compID); //the CBFs up to now have been defined for two sub-TUs - shift them down a level and replace with the parent level CBF
      }
    }

    m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_ROOT ] );
    m_pcEntropyCoder->resetBits();

    if( uiLog2TrSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) )
    {
      m_pcEntropyCoder->encodeTransformSubdivFlag( 0, 5 - uiLog2TrSize );
    }

    for(UInt ch = 0; ch < numValidComp; ch++)
    {
      const UInt chOrderChange = ((ch + 1) == numValidComp) ? 0 : (ch + 1);
      const ComponentID compID=ComponentID(chOrderChange);
      if( rTu.ProcessComponentSection(compID) )
      {
        m_pcEntropyCoder->encodeQtCbf( rTu, compID, true );
      }
    }

    for(UInt ch = 0; ch < numValidComp; ch++)
    {
      const ComponentID compID=ComponentID(ch);
      if (rTu.ProcessComponentSection(compID))
      {
        if(isChroma(compID) && (uiAbsSum[COMPONENT_Y][0] != 0))
        {
          m_pcEntropyCoder->encodeCrossComponentPrediction( rTu, compID );
        }

        m_pcEntropyCoder->encodeCoeffNxN( rTu, pcCoeffCurr[compID], compID );
        for (UInt subTUIndex = 0; subTUIndex < 2; subTUIndex++)
        {
          uiSingleDist += uiSingleDistComp[compID][subTUIndex];
        }
      }
    }

    uiSingleBits = m_pcEntropyCoder->getNumberOfWrittenBits();

    dSingleCost = m_pcRdCost->calcRdCost( uiSingleBits, uiSingleDist );
  } // check full

  // code sub-blocks
  if( bCheckSplit )
  {
    if( bCheckFull )
    {
      m_pcRDGoOnSbacCoder->store( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_TEST ] );
      m_pcRDGoOnSbacCoder->load ( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_ROOT ] );
    }
    Distortion uiSubdivDist = 0;
    UInt       uiSubdivBits = 0;
    Double     dSubdivCost = 0.0;

    //save the non-split CBFs in case we need to restore them later

    UInt bestCBF     [MAX_NUM_COMPONENT];
    UInt bestsubTUCBF[MAX_NUM_COMPONENT][2];
    for(UInt ch = 0; ch < numValidComp; ch++)
    {
      const ComponentID compID=ComponentID(ch);

      if (rTu.ProcessComponentSection(compID))
      {
        bestCBF[compID] = pcCU->getCbf(uiAbsPartIdx, compID, uiTrMode);

        const TComRectangle &tuCompRect = rTu.getRect(compID);
        if (tuCompRect.width != tuCompRect.height)
        {
          const UInt partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> 1;

          for (UInt subTU = 0; subTU < 2; subTU++)
          {
            bestsubTUCBF[compID][subTU] = pcCU->getCbf ((uiAbsPartIdx + (subTU * partIdxesPerSubTU)), compID, subTUDepth);
          }
        }
      }
    }


    TComTURecurse tuRecurseChild(rTu, false);
    const UInt uiQPartNumSubdiv = tuRecurseChild.GetAbsPartIdxNumParts();

    DEBUG_STRING_NEW(sSplitString[MAX_NUM_COMPONENT])

    do
    {
      DEBUG_STRING_NEW(childString)
      xEstimateInterResidualQT( pcResi, dSubdivCost, uiSubdivBits, uiSubdivDist, bCheckFull ? NULL : puiZeroDist,  tuRecurseChild DEBUG_STRING_PASS_INTO(childString));
#if DEBUG_STRING
      // split the string by component and append to the relevant output (because decoder decodes in channel order, whereas this search searches by TU-order)
      std::size_t lastPos=0;
      const std::size_t endStrng=childString.find(debug_reorder_data_inter_token[MAX_NUM_COMPONENT], lastPos);
      for(UInt ch = 0; ch < numValidComp; ch++)
      {
        if (lastPos!=std::string::npos && childString.find(debug_reorder_data_inter_token[ch], lastPos)==lastPos)
        {
          lastPos+=strlen(debug_reorder_data_inter_token[ch]); // skip leading string
        }
        std::size_t pos=childString.find(debug_reorder_data_inter_token[ch+1], lastPos);
        if (pos!=std::string::npos && pos>endStrng)
        {
          lastPos=endStrng;
        }
        sSplitString[ch]+=childString.substr(lastPos, (pos==std::string::npos)? std::string::npos : (pos-lastPos) );
        lastPos=pos;
      }
#endif
    } while ( tuRecurseChild.nextSection(rTu) ) ;

    UInt uiCbfAny=0;
    for(UInt ch = 0; ch < numValidComp; ch++)
    {
      UInt uiYUVCbf = 0;
      for( UInt ui = 0; ui < 4; ++ui )
      {
        uiYUVCbf |= pcCU->getCbf( uiAbsPartIdx + ui * uiQPartNumSubdiv, ComponentID(ch),  uiTrMode + 1 );
      }
      UChar *pBase=pcCU->getCbf( ComponentID(ch) );
      const UInt flags=uiYUVCbf << uiTrMode;
      for( UInt ui = 0; ui < 4 * uiQPartNumSubdiv; ++ui )
      {
        pBase[uiAbsPartIdx + ui] |= flags;
      }
      uiCbfAny|=uiYUVCbf;
    }

    m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_ROOT ] );
    m_pcEntropyCoder->resetBits();

    // when compID isn't a channel, code Cbfs:
    xEncodeInterResidualQT( MAX_NUM_COMPONENT, rTu );
    for(UInt ch = 0; ch < numValidComp; ch++)
    {
      xEncodeInterResidualQT( ComponentID(ch), rTu );
    }

    uiSubdivBits = m_pcEntropyCoder->getNumberOfWrittenBits();
    dSubdivCost  = m_pcRdCost->calcRdCost( uiSubdivBits, uiSubdivDist );

    if (!bCheckFull || (uiCbfAny && (dSubdivCost < dSingleCost)))
    {
      rdCost += dSubdivCost;
      ruiBits += uiSubdivBits;
      ruiDist += uiSubdivDist;
#if DEBUG_STRING
      for(UInt ch = 0; ch < numValidComp; ch++)
      {
        DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[ch])
        DEBUG_STRING_APPEND(sDebug, sSplitString[ch])
      }
#endif
    }
    else
    {
      rdCost  += dSingleCost;
      ruiBits += uiSingleBits;
      ruiDist += uiSingleDist;

      //restore state to unsplit

      pcCU->setTrIdxSubParts( uiTrMode, uiAbsPartIdx, uiDepth );

      for(UInt ch = 0; ch < numValidComp; ch++)
      {
        const ComponentID compID=ComponentID(ch);

        DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[ch])
        if (rTu.ProcessComponentSection(compID))
        {
          DEBUG_STRING_APPEND(sDebug, sSingleStringComp[compID])

          const Bool splitIntoSubTUs   = rTu.getRect(compID).width != rTu.getRect(compID).height;
          const UInt numberOfSections  = splitIntoSubTUs ? 2 : 1;
          const UInt partIdxesPerSubTU = rTu.GetAbsPartIdxNumParts(compID) >> (splitIntoSubTUs ? 1 : 0);

          for (UInt subTUIndex = 0; subTUIndex < numberOfSections; subTUIndex++)
          {
            const UInt  uisubTUPartIdx = uiAbsPartIdx + (subTUIndex * partIdxesPerSubTU);

            if (splitIntoSubTUs)
            {
              const UChar combinedCBF = (bestsubTUCBF[compID][subTUIndex] << subTUDepth) | (bestCBF[compID] << uiTrMode);
              pcCU->setCbfPartRange(combinedCBF, compID, uisubTUPartIdx, partIdxesPerSubTU);
            }
            else
            {
              pcCU->setCbfPartRange((bestCBF[compID] << uiTrMode), compID, uisubTUPartIdx, partIdxesPerSubTU);
            }

            pcCU->setCrossComponentPredictionAlphaPartRange(bestCrossCPredictionAlpha[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
            pcCU->setTransformSkipPartRange(uiBestTransformMode[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
            pcCU->setExplicitRdpcmModePartRange(bestExplicitRdpcmModeUnSplit[compID][subTUIndex], compID, uisubTUPartIdx, partIdxesPerSubTU);
          }
        }
      }

      m_pcRDGoOnSbacCoder->load( m_pppcRDSbacCoder[ uiDepth ][ CI_QT_TRAFO_TEST ] );
    }
  }
  else
  {
    rdCost  += dSingleCost;
    ruiBits += uiSingleBits;
    ruiDist += uiSingleDist;
#if DEBUG_STRING
    for(UInt ch = 0; ch < numValidComp; ch++)
    {
      const ComponentID compID=ComponentID(ch);
      DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[compID])

      if (rTu.ProcessComponentSection(compID))
      {
        DEBUG_STRING_APPEND(sDebug, sSingleStringComp[compID])
      }
    }
#endif
  }
  DEBUG_STRING_APPEND(sDebug, debug_reorder_data_inter_token[MAX_NUM_COMPONENT])
}



Void TEncSearch::xEncodeInterResidualQT( const ComponentID compID, TComTU &rTu )
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();
  const UInt uiCurrTrMode = rTu.GetTransformDepthRel();
  assert( pcCU->getDepth( 0 ) == pcCU->getDepth( uiAbsPartIdx ) );
  const UInt uiTrMode = pcCU->getTransformIdx( uiAbsPartIdx );

  const Bool bSubdiv = uiCurrTrMode != uiTrMode;

  const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();

  if (compID==MAX_NUM_COMPONENT)  // we are not processing a channel, instead we always recurse and code the CBFs
  {
    if( uiLog2TrSize <= pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() && uiLog2TrSize > pcCU->getQuadtreeTULog2MinSizeInCU(uiAbsPartIdx) )
    {
      if((pcCU->getSlice()->getSPS()->getQuadtreeTUMaxDepthInter() == 1) && (pcCU->getPartitionSize(uiAbsPartIdx) != SIZE_2Nx2N))
      {
        assert(bSubdiv); // Inferred splitting rule - see derivation and use of interSplitFlag in the specification.
      }
      else
      {
        m_pcEntropyCoder->encodeTransformSubdivFlag( bSubdiv, 5 - uiLog2TrSize );
      }
    }

    assert( !pcCU->isIntra(uiAbsPartIdx) );

    const Bool bFirstCbfOfCU = uiCurrTrMode == 0;

    for (UInt ch=COMPONENT_Cb; ch<pcCU->getPic()->getNumberValidComponents(); ch++)
    {
      const ComponentID compIdInner=ComponentID(ch);
      if( bFirstCbfOfCU || rTu.ProcessingAllQuadrants(compIdInner) )
      {
        if( bFirstCbfOfCU || pcCU->getCbf( uiAbsPartIdx, compIdInner, uiCurrTrMode - 1 ) )
        {
          m_pcEntropyCoder->encodeQtCbf( rTu, compIdInner, !bSubdiv );
        }
      }
      else
      {
        assert( pcCU->getCbf( uiAbsPartIdx, compIdInner, uiCurrTrMode ) == pcCU->getCbf( uiAbsPartIdx, compIdInner, uiCurrTrMode - 1 ) );
      }
    }

    if (!bSubdiv)
    {
      m_pcEntropyCoder->encodeQtCbf( rTu, COMPONENT_Y, true );
    }
  }

  if( !bSubdiv )
  {
    if (compID != MAX_NUM_COMPONENT) // we have already coded the CBFs, so now we code coefficients
    {
      if (rTu.ProcessComponentSection(compID))
      {
        if (isChroma(compID) && (pcCU->getCbf(uiAbsPartIdx, COMPONENT_Y, uiTrMode) != 0))
        {
          m_pcEntropyCoder->encodeCrossComponentPrediction(rTu, compID);
        }

        if (pcCU->getCbf(uiAbsPartIdx, compID, uiTrMode) != 0)
        {
          const UInt uiQTTempAccessLayer = pcCU->getSlice()->getSPS()->getQuadtreeTULog2MaxSize() - uiLog2TrSize;
          TCoeff *pcCoeffCurr = m_ppcQTTempCoeff[compID][uiQTTempAccessLayer] + rTu.getCoefficientOffset(compID);
          m_pcEntropyCoder->encodeCoeffNxN( rTu, pcCoeffCurr, compID );
        }
      }
    }
  }
  else
  {
    if( compID==MAX_NUM_COMPONENT || pcCU->getCbf( uiAbsPartIdx, compID, uiCurrTrMode ) )
    {
      TComTURecurse tuRecurseChild(rTu, false);
      do
      {
        xEncodeInterResidualQT( compID, tuRecurseChild );
      } while (tuRecurseChild.nextSection(rTu));
    }
  }
}




Void TEncSearch::xSetInterResidualQTData( TComYuv* pcResi, Bool bSpatial, TComTU &rTu ) // TODO: turn this into two functions for bSpatial=true and false.
{
  TComDataCU* pcCU=rTu.getCU();
  const UInt uiCurrTrMode=rTu.GetTransformDepthRel();
  const UInt uiAbsPartIdx=rTu.GetAbsPartIdxTU();
  assert( pcCU->getDepth( 0 ) == pcCU->getDepth( uiAbsPartIdx ) );
  const UInt uiTrMode = pcCU->getTransformIdx( uiAbsPartIdx );
  const TComSPS *sps=pcCU->getSlice()->getSPS();

  if( uiCurrTrMode == uiTrMode )
  {
    const UInt uiLog2TrSize = rTu.GetLog2LumaTrSize();
    const UInt uiQTTempAccessLayer = sps->getQuadtreeTULog2MaxSize() - uiLog2TrSize;

    if( bSpatial )
    {
      // Data to be copied is in the spatial domain, i.e., inverse-transformed.

      for(UInt i=0; i<pcResi->getNumberValidComponents(); i++)
      {
        const ComponentID compID=ComponentID(i);
        if (rTu.ProcessComponentSection(compID))
        {
          const TComRectangle &rectCompTU(rTu.getRect(compID));
          m_pcQTTempTComYuv[uiQTTempAccessLayer].copyPartToPartComponentMxN    ( compID, pcResi, rectCompTU );
        }
      }
    }
    else
    {
      for (UInt ch=0; ch < getNumberValidComponents(sps->getChromaFormatIdc()); ch++)
      {
        const ComponentID compID   = ComponentID(ch);
        if (rTu.ProcessComponentSection(compID))
        {
          const TComRectangle &rectCompTU(rTu.getRect(compID));
          const UInt numCoeffInBlock    = rectCompTU.width * rectCompTU.height;
          const UInt offset             = rTu.getCoefficientOffset(compID);
          TCoeff* dest                  = pcCU->getCoeff(compID)                        + offset;
          const TCoeff* src             = m_ppcQTTempCoeff[compID][uiQTTempAccessLayer] + offset;
          ::memcpy( dest, src, sizeof(TCoeff)*numCoeffInBlock );

#if ADAPTIVE_QP_SELECTION
          TCoeff* pcArlCoeffSrc            = m_ppcQTTempArlCoeff[compID][uiQTTempAccessLayer] + offset;
          TCoeff* pcArlCoeffDst            = pcCU->getArlCoeff(compID)                        + offset;
          ::memcpy( pcArlCoeffDst, pcArlCoeffSrc, sizeof( TCoeff ) * numCoeffInBlock );
#endif
        }
      }
    }
  }
  else
  {

    TComTURecurse tuRecurseChild(rTu, false);
    do
    {
      xSetInterResidualQTData( pcResi, bSpatial, tuRecurseChild );
    } while (tuRecurseChild.nextSection(rTu));
  }
}




UInt TEncSearch::xModeBitsIntra( TComDataCU* pcCU, UInt uiMode, UInt uiPartOffset, UInt uiDepth, const ChannelType chType )
{
  // Reload only contexts required for coding intra mode information
  m_pcRDGoOnSbacCoder->loadIntraDirMode( m_pppcRDSbacCoder[uiDepth][CI_CURR_BEST], chType );

  // Temporarily set the intra dir being tested, and only
  // for absPartIdx, since encodeIntraDirModeLuma/Chroma only use
  // the entry at absPartIdx.

  UChar &rIntraDirVal=pcCU->getIntraDir( chType )[uiPartOffset];
  UChar origVal=rIntraDirVal;
  rIntraDirVal = uiMode;
  //pcCU->setIntraDirSubParts ( chType, uiMode, uiPartOffset, uiDepth + uiInitTrDepth );

  m_pcEntropyCoder->resetBits();
  if (isLuma(chType))
  {
    m_pcEntropyCoder->encodeIntraDirModeLuma ( pcCU, uiPartOffset);
  }
  else
  {
    m_pcEntropyCoder->encodeIntraDirModeChroma ( pcCU, uiPartOffset);
  }

  rIntraDirVal = origVal; // restore

  return m_pcEntropyCoder->getNumberOfWrittenBits();
}




UInt TEncSearch::xUpdateCandList( UInt uiMode, Double uiCost, UInt uiFastCandNum, UInt * CandModeList, Double * CandCostList )
{
  UInt i;
  UInt shift=0;

  while ( shift<uiFastCandNum && uiCost<CandCostList[ uiFastCandNum-1-shift ] )
  {
    shift++;
  }

  if( shift!=0 )
  {
    for(i=1; i<shift; i++)
    {
      CandModeList[ uiFastCandNum-i ] = CandModeList[ uiFastCandNum-1-i ];
      CandCostList[ uiFastCandNum-i ] = CandCostList[ uiFastCandNum-1-i ];
    }
    CandModeList[ uiFastCandNum-shift ] = uiMode;
    CandCostList[ uiFastCandNum-shift ] = uiCost;
    return 1;
  }

  return 0;
}





/** add inter-prediction syntax elements for a CU block
 * \param pcCU
 * \param uiQp
 * \param uiTrMode
 * \param ruiBits
 * \returns Void
 */
Void  TEncSearch::xAddSymbolBitsInter( TComDataCU* pcCU, UInt& ruiBits )
{
  if(pcCU->getMergeFlag( 0 ) && pcCU->getPartitionSize( 0 ) == SIZE_2Nx2N && !pcCU->getQtRootCbf( 0 ))
  {
    pcCU->setSkipFlagSubParts( true, 0, pcCU->getDepth(0) );

    m_pcEntropyCoder->resetBits();
    if(pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
    {
      m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, 0, true);
    }
    m_pcEntropyCoder->encodeSkipFlag(pcCU, 0, true);
    m_pcEntropyCoder->encodeMergeIndex(pcCU, 0, true);

    ruiBits += m_pcEntropyCoder->getNumberOfWrittenBits();
  }
  else
  {
    m_pcEntropyCoder->resetBits();

    if(pcCU->getSlice()->getPPS()->getTransquantBypassEnableFlag())
    {
      m_pcEntropyCoder->encodeCUTransquantBypassFlag(pcCU, 0, true);
    }

    m_pcEntropyCoder->encodeSkipFlag ( pcCU, 0, true );
    m_pcEntropyCoder->encodePredMode( pcCU, 0, true );
    m_pcEntropyCoder->encodePartSize( pcCU, 0, pcCU->getDepth(0), true );
    m_pcEntropyCoder->encodePredInfo( pcCU, 0 );

    Bool codeDeltaQp = false;
    Bool codeChromaQpAdj = false;
    m_pcEntropyCoder->encodeCoeff   ( pcCU, 0, pcCU->getDepth(0), codeDeltaQp, codeChromaQpAdj );

    ruiBits += m_pcEntropyCoder->getNumberOfWrittenBits();
  }
}





/**
 * \brief Generate half-sample interpolated block
 *
 * \param pattern Reference picture ROI
 * \param biPred    Flag indicating whether block is for biprediction
 */
Void TEncSearch::xExtDIFUpSamplingH( TComPattern* pattern )
{
  Int width      = pattern->getROIYWidth();
  Int height     = pattern->getROIYHeight();
  Int srcStride  = pattern->getPatternLStride();

  Int intStride = m_filteredBlockTmp[0].getStride(COMPONENT_Y);
  Int dstStride = m_filteredBlock[0][0].getStride(COMPONENT_Y);
  Pel *intPtr;
  Pel *dstPtr;
  Int filterSize = NTAPS_LUMA;
  Int halfFilterSize = (filterSize>>1);
  Pel *srcPtr = pattern->getROIY() - halfFilterSize*srcStride - 1;

  const ChromaFormat chFmt = m_filteredBlock[0][0].getChromaFormat();

  m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, m_filteredBlockTmp[0].getAddr(COMPONENT_Y), intStride, width+1, height+filterSize, 0, false, chFmt, pattern->getBitDepthY());
  m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, m_filteredBlockTmp[2].getAddr(COMPONENT_Y), intStride, width+1, height+filterSize, 2, false, chFmt, pattern->getBitDepthY());

  intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + halfFilterSize * intStride + 1;
  dstPtr = m_filteredBlock[0][0].getAddr(COMPONENT_Y);
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width+0, height+0, 0, false, true, chFmt, pattern->getBitDepthY());

  intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride + 1;
  dstPtr = m_filteredBlock[2][0].getAddr(COMPONENT_Y);
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width+0, height+1, 2, false, true, chFmt, pattern->getBitDepthY());

  intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + halfFilterSize * intStride;
  dstPtr = m_filteredBlock[0][2].getAddr(COMPONENT_Y);
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width+1, height+0, 0, false, true, chFmt, pattern->getBitDepthY());

  intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
  dstPtr = m_filteredBlock[2][2].getAddr(COMPONENT_Y);
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width+1, height+1, 2, false, true, chFmt, pattern->getBitDepthY());
}





/**
 * \brief Generate quarter-sample interpolated blocks
 *
 * \param pattern    Reference picture ROI
 * \param halfPelRef Half-pel mv
 * \param biPred     Flag indicating whether block is for biprediction
 */
Void TEncSearch::xExtDIFUpSamplingQ( TComPattern* pattern, TComMv halfPelRef )
{
  Int width      = pattern->getROIYWidth();
  Int height     = pattern->getROIYHeight();
  Int srcStride  = pattern->getPatternLStride();

  Pel *srcPtr;
  Int intStride = m_filteredBlockTmp[0].getStride(COMPONENT_Y);
  Int dstStride = m_filteredBlock[0][0].getStride(COMPONENT_Y);
  Pel *intPtr;
  Pel *dstPtr;
  Int filterSize = NTAPS_LUMA;

  Int halfFilterSize = (filterSize>>1);

  Int extHeight = (halfPelRef.getVer() == 0) ? height + filterSize : height + filterSize-1;

  const ChromaFormat chFmt = m_filteredBlock[0][0].getChromaFormat();

  // Horizontal filter 1/4
  srcPtr = pattern->getROIY() - halfFilterSize * srcStride - 1;
  intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y);
  if (halfPelRef.getVer() > 0)
  {
    srcPtr += srcStride;
  }
  if (halfPelRef.getHor() >= 0)
  {
    srcPtr += 1;
  }
  m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, intPtr, intStride, width, extHeight, 1, false, chFmt, pattern->getBitDepthY());

  // Horizontal filter 3/4
  srcPtr = pattern->getROIY() - halfFilterSize*srcStride - 1;
  intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y);
  if (halfPelRef.getVer() > 0)
  {
    srcPtr += srcStride;
  }
  if (halfPelRef.getHor() > 0)
  {
    srcPtr += 1;
  }
  m_if.filterHor(COMPONENT_Y, srcPtr, srcStride, intPtr, intStride, width, extHeight, 3, false, chFmt, pattern->getBitDepthY());

  // Generate @ 1,1
  intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
  dstPtr = m_filteredBlock[1][1].getAddr(COMPONENT_Y);
  if (halfPelRef.getVer() == 0)
  {
    intPtr += intStride;
  }
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1, false, true, chFmt, pattern->getBitDepthY());

  // Generate @ 3,1
  intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
  dstPtr = m_filteredBlock[3][1].getAddr(COMPONENT_Y);
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3, false, true, chFmt, pattern->getBitDepthY());

  if (halfPelRef.getVer() != 0)
  {
    // Generate @ 2,1
    intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
    dstPtr = m_filteredBlock[2][1].getAddr(COMPONENT_Y);
    if (halfPelRef.getVer() == 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 2, false, true, chFmt, pattern->getBitDepthY());

    // Generate @ 2,3
    intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
    dstPtr = m_filteredBlock[2][3].getAddr(COMPONENT_Y);
    if (halfPelRef.getVer() == 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 2, false, true, chFmt, pattern->getBitDepthY());
  }
  else
  {
    // Generate @ 0,1
    intPtr = m_filteredBlockTmp[1].getAddr(COMPONENT_Y) + halfFilterSize * intStride;
    dstPtr = m_filteredBlock[0][1].getAddr(COMPONENT_Y);
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 0, false, true, chFmt, pattern->getBitDepthY());

    // Generate @ 0,3
    intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + halfFilterSize * intStride;
    dstPtr = m_filteredBlock[0][3].getAddr(COMPONENT_Y);
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 0, false, true, chFmt, pattern->getBitDepthY());
  }

  if (halfPelRef.getHor() != 0)
  {
    // Generate @ 1,2
    intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
    dstPtr = m_filteredBlock[1][2].getAddr(COMPONENT_Y);
    if (halfPelRef.getHor() > 0)
    {
      intPtr += 1;
    }
    if (halfPelRef.getVer() >= 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1, false, true, chFmt, pattern->getBitDepthY());

    // Generate @ 3,2
    intPtr = m_filteredBlockTmp[2].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
    dstPtr = m_filteredBlock[3][2].getAddr(COMPONENT_Y);
    if (halfPelRef.getHor() > 0)
    {
      intPtr += 1;
    }
    if (halfPelRef.getVer() > 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3, false, true, chFmt, pattern->getBitDepthY());
  }
  else
  {
    // Generate @ 1,0
    intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride + 1;
    dstPtr = m_filteredBlock[1][0].getAddr(COMPONENT_Y);
    if (halfPelRef.getVer() >= 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1, false, true, chFmt, pattern->getBitDepthY());

    // Generate @ 3,0
    intPtr = m_filteredBlockTmp[0].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride + 1;
    dstPtr = m_filteredBlock[3][0].getAddr(COMPONENT_Y);
    if (halfPelRef.getVer() > 0)
    {
      intPtr += intStride;
    }
    m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3, false, true, chFmt, pattern->getBitDepthY());
  }

  // Generate @ 1,3
  intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
  dstPtr = m_filteredBlock[1][3].getAddr(COMPONENT_Y);
  if (halfPelRef.getVer() == 0)
  {
    intPtr += intStride;
  }
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 1, false, true, chFmt, pattern->getBitDepthY());

  // Generate @ 3,3
  intPtr = m_filteredBlockTmp[3].getAddr(COMPONENT_Y) + (halfFilterSize-1) * intStride;
  dstPtr = m_filteredBlock[3][3].getAddr(COMPONENT_Y);
  m_if.filterVer(COMPONENT_Y, intPtr, intStride, dstPtr, dstStride, width, height, 3, false, true, chFmt, pattern->getBitDepthY());
}





//! set wp tables
Void  TEncSearch::setWpScalingDistParam( TComDataCU* pcCU, Int iRefIdx, RefPicList eRefPicListCur )
{
  if ( iRefIdx<0 )
  {
    m_cDistParam.bApplyWeight = false;
    return;
  }

  TComSlice       *pcSlice  = pcCU->getSlice();
  WPScalingParam  *wp0 , *wp1;

  m_cDistParam.bApplyWeight = ( pcSlice->getSliceType()==P_SLICE && pcSlice->testWeightPred() ) || ( pcSlice->getSliceType()==B_SLICE && pcSlice->testWeightBiPred() ) ;

  if ( !m_cDistParam.bApplyWeight )
  {
    return;
  }

  Int iRefIdx0 = ( eRefPicListCur == REF_PIC_LIST_0 ) ? iRefIdx : (-1);
  Int iRefIdx1 = ( eRefPicListCur == REF_PIC_LIST_1 ) ? iRefIdx : (-1);

  getWpScaling( pcCU, iRefIdx0, iRefIdx1, wp0 , wp1 );

  if ( iRefIdx0 < 0 )
  {
    wp0 = NULL;
  }
  if ( iRefIdx1 < 0 )
  {
    wp1 = NULL;
  }

  m_cDistParam.wpCur  = NULL;

  if ( eRefPicListCur == REF_PIC_LIST_0 )
  {
    m_cDistParam.wpCur = wp0;
  }
  else
  {
    m_cDistParam.wpCur = wp1;
  }
}



//! \}
