#ifndef lint
static const char	RCSid[] = "$Id: spec_rgb.c,v 2.31 2023/12/08 18:48:09 greg Exp $";
#endif
/*
 * Convert colors and spectral ranges.
 * Added von Kries white-balance calculations 10/01 (GW).
 *
 * Externals declared in color.h
 */

#include "copyright.h"

#include <stdio.h>
#include <string.h>
#include "color.h"

#define CEPS	1e-4			/* color epsilon */

#define CEQ(v1,v2)	(((v1) <= (v2)+CEPS) & ((v2) <= (v1)+CEPS))

#define XYEQ(c1,c2)	(CEQ((c1)[CIEX],(c2)[CIEX]) & CEQ((c1)[CIEY],(c2)[CIEY]))


RGBPRIMS  stdprims = STDPRIMS;	/* standard primary chromaticities */
RGBPRIMS  xyzprims = {1,0, 0,1, 0,0, 1.f/3.f,1.f/3.f};

COLOR  cblack = BLKCOLOR;		/* global black color */
COLOR  cwhite = WHTCOLOR;		/* global white color */

float  xyneu[2] = {1./3., 1./3.};	/* neutral xy chromaticities */

/*
 *	The tables below encode the CIE 1931 2° standard observer
 * functions for X, Y, and Z.  These tables are cumulative, as are
 * ones that follow.
 */

#define	CIE_X_WLMIN	362
#define CIE_X_WLMAX	774

static const unsigned short	cie_x_cumul[CIE_X_WLMAX-CIE_X_WLMIN+1] = {
0,1,1,1,1,1,1,2,2,2,3,3,3,4,5,5,6,7,8,9,10,11,12,14,16,18,20,23,26,29,
33,37,41,47,53,60,68,77,86,97,108,121,135,151,170,190,214,241,271,304,
342,385,432,486,545,612,686,768,860,961,1073,1195,1326,1467,1618,
1776,1943,2117,2298,2485,2677,2875,3076,3281,3489,3700,3912,4125,4340,4555,
4769,4983,5197,5409,5620,5830,6038,6244,6448,6651,6851,7049,7245,7437,
7627,7813,7995,8173,8347,8517,8682,8842,8996,9143,9284,9418,9545,9665,
9778,9884,9984,10077,10164,10245,10321,10390,10454,10513,10566,10615,
10659,10698,10734,10766,10794,10819,10842,10861,10879,10893,10906,10917,
10926,10933,10939,10944,10948,10951,10953,10955,10957,10958,10960,10961,
10964,10967,10971,10977,10984,10994,11006,11020,11038,11060,11085,11114,
11148,11187,11231,11280,11335,11396,11464,11537,11618,11705,11799,11901,
12010,12126,12249,12380,12518,12664,12818,12980,13150,13328,13515,13709,
13913,14124,14345,14574,14813,15060,15317,15582,15858,16142,16437,16741,
17055,17379,17713,18057,18412,18776,19151,19536,19931,20337,20753,21180,
21616,22063,22520,22988,23465,23953,24450,24957,25474,26000,26535,27080,
27633,28195,28765,29343,29929,30522,31122,31729,32342,32961,33585,34215,
34849,35487,36129,36775,37423,38073,38724,39375,40027,40679,41329,41978,
42625,43270,43911,44548,45181,45808,46429,47044,47652,48253,48845,49430,
50005,50571,51128,51674,52209,52733,53245,53745,54232,54706,55167,55614,
56048,56468,56876,57269,57651,58019,58376,58720,59052,59373,59681,59979,
60265,60539,60803,61056,61298,61529,61750,61962,62163,62355,62538,62712,
62877,63034,63183,63325,63459,63586,63706,63820,63927,64028,64123,64213,
64297,64376,64451,64520,64586,64647,64704,64758,64808,64855,64899,64941,
64980,65016,65051,65083,65114,65143,65169,65194,65218,65240,65260,65278,
65296,65312,65327,65341,65354,65366,65377,65387,65397,65406,65415,65423,
65430,65437,65444,65450,65455,65461,65466,65470,65475,65479,65483,65486,
65489,65493,65495,65498,65501,65503,65505,65507,65509,65511,65513,65514,
65516,65517,65518,65519,65520,65521,65522,65523,65524,65525,65526,65526,
65527,65527,65528,65528,65529,65529,65530,65530,65530,65531,65531,65531,
65532,65532,65532,65532,65532,65533,65533,65533,65533,65533,65533,65533,
65533,65534,65534,65534,65534,65534,65534,65534,65534,65534,65534,65534,
65534,65534,65534,65534,65535
};

#define CIE_Y_WLMIN	386
#define CIE_Y_WLMAX	760

static const unsigned short	cie_y_cumul[CIE_Y_WLMAX-CIE_Y_WLMIN+1] = {
0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,5,5,6,7,8,8,10,11,12,14,
15,17,19,22,25,28,31,35,40,45,50,56,63,70,78,86,95,105,115,126,
138,150,164,178,193,208,225,242,260,280,300,321,343,367,391,
417,443,472,501,532,564,598,633,670,708,748,790,833,879,926,
975,1027,1080,1136,1194,1255,1318,1384,1453,1525,1601,1679,1761,
1846,1935,2027,2123,2223,2327,2435,2548,2665,2787,2915,3048,3187,
3332,3484,3643,3808,3981,4162,4352,4550,4757,4975,5203,5442,5691,
5952,6225,6509,6805,7114,7435,7769,8116,8476,8849,9235,9634,10045,
10469,10904,11351,11808,12275,12752,13239,13734,14239,14752,15273,
15801,16337,16880,17429,17984,18546,19112,19684,20260,20841,21426,
22015,22608,23203,23802,24403,25007,25612,26220,26828,27439,28050,
28662,29275,29888,30501,31114,31727,32340,32951,33561,34170,34777,
35382,35985,36585,37182,37777,38368,38955,39539,40119,40695,41266,
41832,42393,42950,43501,44046,44586,45119,45646,46167,46682,47189,
47690,48183,48670,49149,49621,50085,50542,50991,51432,51866,52293,
52711,53122,53524,53919,54306,54685,55057,55420,55775,56123,56463,
56795,57119,57435,57743,58044,58337,58623,58901,59172,59435,59691,
59939,60180,60414,60640,60859,61070,61274,61471,61661,61844,62019,
62188,62351,62507,62657,62802,62940,63073,63201,63323,63441,63553,
63660,63763,63861,63954,64043,64128,64208,64285,64358,64427,64493,
64555,64614,64670,64723,64773,64820,64865,64907,64947,64984,65019,
65052,65083,65113,65140,65166,65190,65212,65233,65253,65271,65288,
65304,65319,65334,65347,65360,65371,65383,65393,65403,65412,65420,
65428,65435,65442,65449,65454,65460,65465,65470,65474,65478,65482,
65485,65488,65492,65494,65497,65500,65502,65504,65506,65508,65510,
65512,65513,65515,65516,65517,65519,65520,65521,65522,65523,65523,
65524,65525,65526,65526,65527,65527,65528,65528,65529,65529,65530,
65530,65530,65531,65531,65531,65532,65532,65532,65532,65532,65533,
65533,65533,65533,65533,65533,65533,65534,65534,65534,65534,65534,
65534,65534,65534,65534,65534,65534,65534,65534,65534,65534,65534,65535
};

#define CIE_Z_WLMIN	359
#define CIE_Z_WLMAX	624

static const unsigned short	cie_z_cumul[CIE_Z_WLMAX-CIE_Z_WLMIN+1] = {
0,1,1,2,2,3,4,5,6,7,8,9,11,12,14,16,19,22,25,28,32,37,41,47,52,59,66,75,84,
95,107,121,136,154,173,196,221,250,283,321,362,408,458,512,573,640,717,804,
903,1015,1142,1285,1446,1627,1830,2058,2313,2598,2917,3272,3668,4108,4597,
5135,5723,6360,7044,7773,8544,9356,10205,11090,12006,12952,13923,14918,
15934,16967,18016,19076,20148,21227,22312,23401,24492,25585,26678,27771,
28861,29950,31037,32121,33202,34281,35355,36424,37487,38542,39588,40624,
41647,42657,43652,44631,45590,46527,47438,48321,49173,49993,50783,51542,
52270,52968,53636,54275,54885,55465,56018,56543,57042,57514,57961,58384,
58784,59162,59519,59856,60175,60477,60762,61032,61287,61529,61757,61974,
62179,62374,62559,62734,62901,63059,63211,63355,63492,63622,63745,63862,
63971,64075,64172,64263,64347,64427,64500,64569,64632,64692,64747,64798,
64846,64891,64933,64973,65010,65045,65078,65109,65139,65166,65192,65216,
65239,65260,65280,65298,65315,65331,65345,65359,65371,65383,65393,65403,
65412,65420,65428,65435,65441,65447,65452,65457,65462,65466,65470,65473,
65476,65479,65482,65485,65487,65489,65491,65493,65495,65497,65498,65500,
65501,65503,65504,65505,65506,65508,65509,65510,65511,65512,65513,65514,
65515,65516,65517,65518,65519,65520,65520,65521,65522,65523,65523,65524,
65525,65525,65526,65527,65527,65528,65528,65529,65529,65530,65530,65531,
65531,65531,65532,65532,65532,65532,65533,65533,65533,65533,65533,65533,
65534,65534,65534,65534,65534,65534,65534,65534,65534,65535
};

#define SCOTOPIC_WLMIN	382
#define SCOTOPIC_WLMAX	684

static const unsigned short	scotopic_cumul[SCOTOPIC_WLMAX-SCOTOPIC_WLMIN+1] = {
0, 1, 1, 2, 3, 4, 5, 6, 7, 9, 10, 12, 15, 17, 20, 24, 28, 33, 38, 44,
52, 60, 69, 80, 93, 107, 123, 142, 163, 186, 213, 242, 275, 312, 353,
398, 448, 502, 562, 627, 698, 775, 859, 949, 1046, 1150, 1261, 1380,
1507, 1642, 1785, 1936, 2096, 2265, 2442, 2628, 2823, 3027, 3239, 3461,
3691, 3930, 4178, 4435, 4700, 4974, 5257, 5548, 5847, 6154, 6469, 6793,
7123, 7462, 7809, 8162, 8524, 8892, 9268, 9651, 10041, 10438, 10843, 11254,
11673, 12099, 12532, 12973, 13422, 13878, 14342, 14814, 15293, 15780, 16276,
16779, 17290, 17809, 18336, 18872, 19415, 19967, 20526, 21093, 21667, 22249,
22839, 23436, 24039, 24649, 25266, 25890, 26519, 27154, 27795, 28441, 29092,
29747, 30405, 31068, 31734, 32402, 33074, 33747, 34420, 35095, 35771, 36446,
37119, 37793, 38464, 39132, 39798, 40460, 41118, 41772, 42421, 43064, 43701,
44332, 44957, 45575, 46185, 46787, 47381, 47967, 48543, 49110, 49668, 50215,
50753, 51280, 51797, 52302, 52797, 53281, 53754, 54215, 54665, 55104, 55531,
55947, 56352, 56744, 57125, 57495, 57853, 58200, 58536, 58860, 59174, 59477,
59769, 60051, 60322, 60583, 60834, 61075, 61306, 61528, 61741, 61944, 62139,
62326, 62504, 62674, 62836, 62991, 63138, 63278, 63412, 63538, 63659, 63773,
63881, 63983, 64080, 64172, 64259, 64340, 64418, 64490, 64559, 64623, 64684,
64742, 64795, 64846, 64893, 64937, 64978, 65017, 65053, 65087, 65119, 65149,
65176, 65202, 65226, 65248, 65269, 65289, 65307, 65323, 65339, 65353, 65367,
65379, 65391, 65402, 65412, 65421, 65430, 65438, 65445, 65452, 65458, 65464,
65469, 65474, 65479, 65483, 65487, 65491, 65494, 65497, 65500, 65503, 65505,
65507, 65509, 65511, 65513, 65514, 65516, 65517, 65519, 65520, 65521, 65522,
65523, 65524, 65525, 65525, 65526, 65527, 65527, 65528, 65528, 65529, 65529,
65529, 65530, 65530, 65530, 65531, 65531, 65531, 65531, 65532, 65532, 65532,
65532, 65532, 65533, 65533, 65533, 65533, 65533, 65533, 65533, 65533, 65533,
65533, 65534, 65534, 65534, 65534, 65534, 65534, 65534, 65534, 65535
};

#define MELANOPIC_WLMIN	380
#define MELANOPIC_WLMAX	651

static const unsigned short	melanopic_cumul[MELANOPIC_WLMAX-MELANOPIC_WLMIN+1] = {
0, 1, 2, 3, 4, 5, 7, 8, 10, 12, 14, 17, 20, 23, 27, 32, 37, 42, 49, 56,
65, 75, 86, 99, 114, 131, 150, 173, 200, 230, 265, 303, 347, 395, 448,
508, 574, 650, 734, 828, 931, 1041, 1158, 1283, 1415, 1554, 1703, 1862,
2031, 2211, 2400, 2600, 2809, 3028, 3257, 3497, 3748, 4012, 4288, 4576,
4876, 5187, 5509, 5842, 6185, 6540, 6905, 7283, 7673, 8075, 8489, 8915,
9351, 9799, 10259, 10729, 11211, 11705, 12211, 12729, 13258, 13799,
14351, 14915, 15491, 16078, 16676, 17286, 17908, 18541, 19184, 19837,
20498, 21168, 21846, 22532, 23226, 23928, 24637, 25353, 26075, 26802,
27532, 28267, 29005, 29746, 30488, 31233, 31979, 32726, 33473, 34220,
34966, 35711, 36455, 37196, 37935, 38671, 39403, 40129, 40851, 41568,
42279, 42983, 43680, 44369, 45050, 45724, 46388, 47043, 47688, 48322,
48945, 49557, 50156, 50743, 51317, 51879, 52428, 52964, 53487, 53997,
54493, 54976, 55445, 55901, 56343, 56771, 57186, 57587, 57976, 58351,
58712, 59061, 59397, 59720, 60031, 60330, 60616, 60891, 61154, 61405,
61646, 61875, 62094, 62303, 62501, 62690, 62870, 63040, 63201, 63354,
63498, 63634, 63763, 63884, 63998, 64105, 64206, 64300, 64389, 64472,
64549, 64622, 64689, 64753, 64811, 64866, 64917, 64964, 65008, 65049,
65086, 65121, 65154, 65184, 65211, 65237, 65260, 65282, 65302, 65321,
65338, 65354, 65368, 65381, 65394, 65405, 65415, 65425, 65434, 65442,
65449, 65456, 65463, 65468, 65474, 65479, 65483, 65487, 65491, 65494,
65498, 65501, 65503, 65506, 65508, 65510, 65512, 65514, 65515, 65517,
65518, 65520, 65521, 65522, 65523, 65524, 65525, 65525, 65526, 65527,
65527, 65528, 65528, 65529, 65529, 65530, 65530, 65530, 65531, 65531,
65531, 65531, 65532, 65532, 65532, 65532, 65532, 65533, 65533, 65533,
65533, 65533, 65533, 65533, 65533, 65533, 65534, 65534, 65534, 65535
};

COLORMAT  xyz2rgbmat = {		/* XYZ to RGB (no white balance) */
	{(CIE_y_g - CIE_y_b - CIE_x_b*CIE_y_g + CIE_y_b*CIE_x_g)/CIE_C_rD,
	 (CIE_x_b - CIE_x_g - CIE_x_b*CIE_y_g + CIE_x_g*CIE_y_b)/CIE_C_rD,
	 (CIE_x_g*CIE_y_b - CIE_x_b*CIE_y_g)/CIE_C_rD},
	{(CIE_y_b - CIE_y_r - CIE_y_b*CIE_x_r + CIE_y_r*CIE_x_b)/CIE_C_gD,
	 (CIE_x_r - CIE_x_b - CIE_x_r*CIE_y_b + CIE_x_b*CIE_y_r)/CIE_C_gD,
	 (CIE_x_b*CIE_y_r - CIE_x_r*CIE_y_b)/CIE_C_gD},
	{(CIE_y_r - CIE_y_g - CIE_y_r*CIE_x_g + CIE_y_g*CIE_x_r)/CIE_C_bD,
	 (CIE_x_g - CIE_x_r - CIE_x_g*CIE_y_r + CIE_x_r*CIE_y_g)/CIE_C_bD,
	 (CIE_x_r*CIE_y_g - CIE_x_g*CIE_y_r)/CIE_C_bD}
};

COLORMAT  rgb2xyzmat = {		/* RGB to XYZ (no white balance) */
	{CIE_x_r*CIE_C_rD/CIE_D,CIE_x_g*CIE_C_gD/CIE_D,CIE_x_b*CIE_C_bD/CIE_D},
	{CIE_y_r*CIE_C_rD/CIE_D,CIE_y_g*CIE_C_gD/CIE_D,CIE_y_b*CIE_C_bD/CIE_D},
	{(1.-CIE_x_r-CIE_y_r)*CIE_C_rD/CIE_D,
	 (1.-CIE_x_g-CIE_y_g)*CIE_C_gD/CIE_D,
	 (1.-CIE_x_b-CIE_y_b)*CIE_C_bD/CIE_D}
};

COLORMAT  vkmat = {		/* Sharp primary matrix */
	{ 1.2694, -0.0988, -0.1706},
	{-0.8364,  1.8006,  0.0357},
	{ 0.0297, -0.0315,  1.0018}
};

COLORMAT  ivkmat = {		/* inverse Sharp primary matrix */
	{ 0.8156,  0.0472,  0.1372},
	{ 0.3791,  0.5769,  0.0440},
	{-0.0123,  0.0167,  0.9955}
};


void
spec_cie(			/* compute XYZ color from a spectral range */
COLOR  col,		/* returned color */
int  s,			/* starting and ending wavelengths */
int  e
)
{
	if (s > e) { int i = s; s = e; e = i; }

	if ((s >= CIE_X_WLMAX) | (e <= CIE_X_WLMIN))
		col[CIEX] = 0;
	else
		col[CIEX] = (cie_x_cumul[e>=CIE_X_WLMAX ? CIE_X_WLMAX-CIE_X_WLMIN : e-CIE_X_WLMIN] -
			cie_x_cumul[(s>CIE_X_WLMIN)*(s-CIE_X_WLMIN)])*(1./65535.);

	if ((s >= CIE_Y_WLMAX) | (e <= CIE_Y_WLMIN))
		col[CIEY] = 0;
	else
		col[CIEY] = (cie_y_cumul[e>=CIE_Y_WLMAX ? CIE_Y_WLMAX-CIE_Y_WLMIN : e-CIE_Y_WLMIN] -
			cie_y_cumul[(s>CIE_Y_WLMIN)*(s-CIE_Y_WLMIN)])*(1./65535.);

	if ((s >= CIE_Z_WLMAX) | (e <= CIE_Z_WLMIN))
		col[CIEZ] = 0;
	else
		col[CIEZ] = (cie_z_cumul[e>=CIE_Z_WLMAX ? CIE_Z_WLMAX-CIE_Z_WLMIN : e-CIE_Z_WLMIN] -
			cie_z_cumul[(s>CIE_Z_WLMIN)*(s-CIE_Z_WLMIN)])*(1./65535.);
}


void
spec_rgb(			/* compute RGB color from spectral range */
COLOR  col,
int  s,
int  e
)
{
	COLOR  ciecolor;

	spec_cie(ciecolor, s, e);
	colortrans(col, xyz2rgbmat, ciecolor);
}


static double
spec_dot(			/* spectrum dot-product with cumulative observer */
	SCOLOR scol,
	int  ncs,
	const float wlpt[4],
	const unsigned short cumul[],
	int  wlmin,
	int  wlmax
)
{
	const double	wlstp = (wlpt[0] - wlpt[3])/(double)ncs;
	int		wl1 = (int)wlpt[3];
	int		n = 0;
	double		sum = 0;

	if (wl1 < wlmin) {
		n += (int)((wlmin - wl1)/wlstp);
		wl1 = wlmin;
	} else if (wl1 >= wlmax)
		return(.0);

	while (n < ncs) {
		int	wl0 = wl1;
		wl1 = (int)(ncs==3 ? wlpt[2-n] : wlpt[3] + (n+1)*wlstp);
		if (wl1 >= wlmax) {
			sum += (65535 - cumul[wl0-wlmin]) * scol[ncs-1-n];
			break;
		}
		sum += (cumul[wl1-wlmin] - cumul[wl0-wlmin]) * scol[ncs-1-n++];
	}
	return(sum * (1./65535.));
}


void
scolor2cie(			/* accurate conversion from spectrum to XYZ */
	COLOR col,
	SCOLOR scol,
	int ncs,
	const float wlpt[4]
)
{
	if (ncs == 3) {		/* not a spectrum! */
		rgb_cie(col, scol);
		return;
	}
	col[CIEX] = spec_dot(scol, ncs, wlpt, cie_x_cumul, CIE_X_WLMIN, CIE_X_WLMAX);
	col[CIEY] = spec_dot(scol, ncs, wlpt, cie_y_cumul, CIE_Y_WLMIN, CIE_Y_WLMAX);
	col[CIEZ] = spec_dot(scol, ncs, wlpt, cie_z_cumul, CIE_Z_WLMIN, CIE_Z_WLMAX);
}


void
scolor2rgb(			/* accurate conversion from spectrum to RGB */
	COLOR col,
	SCOLOR scol,
	int ncs,
	const float wlpt[4]
)
{
	COLOR	ciecolor;

	if (ncs == 3) {		/* not really a spectrum, so we punt... */
		copycolor(col, scol);
		return;
	}
	scolor2cie(ciecolor, scol, ncs, wlpt);
	cie_rgb(col, ciecolor);
}


double
scolor2photopic(		/* compute scotopic integral for spectral color */
	SCOLOR  scol,
	int ncs,
	const float wlpt[4]
)
{
	if (ncs == 3)
		return bright(scol);

	return(spec_dot(scol, ncs, wlpt, cie_y_cumul, CIE_Y_WLMIN, CIE_Y_WLMAX));
}


double
scolor2scotopic(		/* compute Y channel for spectral color */
	SCOLOR  scol,
	int ncs,
	const float wlpt[4]
)
{
	return(spec_dot(scol, ncs, wlpt, scotopic_cumul, SCOTOPIC_WLMIN, SCOTOPIC_WLMAX));
}


double
scolor2melanopic(		/* compute melanopic integral for spectral color */
	SCOLOR  scol,
	int ncs,
	const float wlpt[4]
)
{
	return(spec_dot(scol, ncs, wlpt, melanopic_cumul, MELANOPIC_WLMIN, MELANOPIC_WLMAX));
}


double
scolor_photopic(		/* compute scotopic integral for spectral color */
	SCOLOR  scol
)
{
	return(scolor2photopic(scol, NCSAMP, WLPART));
}


double
scolor_scotopic(		/* compute Y channel for spectral color */
	SCOLOR  scol
)
{
	return(scolor2scotopic(scol, NCSAMP, WLPART));
}


double
scolor_melanopic(		/* compute melanopic integral for spectral color */
	SCOLOR  scol
)
{
	return(scolor2melanopic(scol, NCSAMP, WLPART));
}


void
cie_rgb(			/* convert CIE color to standard RGB */
COLOR	rgb,
COLOR  xyz
)
{
	colortrans(rgb, xyz2rgbmat, xyz);
	clipgamut(rgb, xyz[CIEY], CGAMUT_LOWER, cblack, cwhite);
}


int
clipgamut(			/* clip to gamut cube */
COLOR  col,
double  brt,
int  gamut,
COLOR  lower,
COLOR  upper
)
{
	int  rflags = 0;
	double  brtmin, brtmax, v, vv;
	COLOR  cgry;
	int  i;
					/* check for no check */
	if (gamut == 0) return(0);
					/* check brightness limits */
	brtmin = 1./3.*(lower[0]+lower[1]+lower[2]);
	if (gamut & CGAMUT_LOWER && brt < brtmin) {
		copycolor(col, lower);
		return(CGAMUT_LOWER);
	}
	brtmax = 1./3.*(upper[0]+upper[1]+upper[2]);
	if (gamut & CGAMUT_UPPER && brt > brtmax) {
		copycolor(col, upper);
		return(CGAMUT_UPPER);
	}
					/* compute equivalent grey */
	v = (brt - brtmin)/(brtmax - brtmin);
	for (i = 0; i < 3; i++)
		cgry[i] = v*upper[i] + (1.-v)*lower[i];
	vv = 1.;			/* check each limit */
	for (i = 0; i < 3; i++)
		if (gamut & CGAMUT_LOWER && col[i] < lower[i]) {
			v = (lower[i] - cgry[i])/(col[i] - cgry[i]);
			if (v < vv) vv = v;
			rflags |= CGAMUT_LOWER;
		} else if (gamut & CGAMUT_UPPER && col[i] > upper[i]) {
			v = (upper[i] - cgry[i])/(col[i] - cgry[i]);
			if (v < vv) vv = v;
			rflags |= CGAMUT_UPPER;
		}
	if (rflags)			/* desaturate to cube face */
		for (i = 0; i < 3; i++)
			col[i] = vv*col[i] + (1.-vv)*cgry[i];
	return(rflags);
}


void
colortrans(			/* convert c1 by mat and put into c2 */
COLOR  c2,
COLORMAT  mat,
COLOR  c1
)
{
	COLOR	cout;

	cout[0] = mat[0][0]*c1[0] + mat[0][1]*c1[1] + mat[0][2]*c1[2];
	cout[1] = mat[1][0]*c1[0] + mat[1][1]*c1[1] + mat[1][2]*c1[2];
	cout[2] = mat[2][0]*c1[0] + mat[2][1]*c1[1] + mat[2][2]*c1[2];

	copycolor(c2, cout);
}


void
multcolormat(			/* multiply m1 by m2 and put into m3 */
COLORMAT  m3,			/* m3 can be either m1 or m2 w/o harm */
COLORMAT  m2,
COLORMAT  m1
)
{
	COLORMAT  mt;
	int  i, j;

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			mt[i][j] =	m1[i][0]*m2[0][j] +
					m1[i][1]*m2[1][j] +
					m1[i][2]*m2[2][j] ;
	cpcolormat(m3, mt);
}


int
colorprimsOK(			/* are color primaries reasonable? */
RGBPRIMS  pr
)
{
	int	i, j;
				/* check white point */
	if ((pr[3][CIEX] <= CEPS) | (pr[3][CIEX] >= 1.-CEPS) |
			(pr[3][CIEY] <= CEPS) | (pr[3][CIEY] >= 1.-CEPS))
		return(0);
	for (i = 3; i--; )	/* check for XYZ color primaries */
		if (!CEQ(pr[i][CIEX],(i==0)) | !CEQ(pr[i][CIEY],(i==1)))
			break;
	if (i < 0)
		return(-1);	/* flag as XYZ color space */
				/* check color primaries */
	for (i = 0; i < 3; i++) {
		if ((pr[i][CIEX] <= -2.) | (pr[i][CIEY] <= -2.))
			return(0);
		if ((pr[i][CIEX] >= 3.) | (pr[i][CIEY] >= 3.))
			return(0);
		if (pr[i][CIEX] + pr[i][CIEY] <= -2.)
			return(0);
		if (pr[i][CIEX] + pr[i][CIEY] >= 3.)
			return(0);
	}
	for (i = 0; i < 4; i++)	/* make sure space is 3-dimensional */
		for (j = i+1; j < 4; j++)
			if (CEQ(pr[i][CIEX],pr[j][CIEX]) &
					CEQ(pr[i][CIEY],pr[j][CIEY]))
				return(0);
	return(1);
}



int
compxyz2rgbmat(			/* compute conversion from CIE to RGB space */
COLORMAT  mat,
RGBPRIMS  pr
)
{
	double  C_rD, C_gD, C_bD;

	if (pr == stdprims) {	/* can use xyz2rgbmat */
		cpcolormat(mat, xyz2rgbmat);
		return(1);
	}
	if (pr == xyzprims) {	/* identity */
		memset(mat, 0, sizeof(COLORMAT));
		mat[0][0] = mat[1][1] = mat[2][2] = 1.f;
		return(1);
	}
	if (CEQ(pr[WHT][CIEX],0.) | CEQ(pr[WHT][CIEY],0.))
		return(0);
	C_rD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[GRN][CIEY] - pr[BLU][CIEY]) -
			  pr[WHT][CIEY]*(pr[GRN][CIEX] - pr[BLU][CIEX]) +
		  pr[GRN][CIEX]*pr[BLU][CIEY] - pr[BLU][CIEX]*pr[GRN][CIEY] ) ;
	if (CEQ(C_rD,0.))
		return(0);
	C_gD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[BLU][CIEY] - pr[RED][CIEY]) -
			  pr[WHT][CIEY]*(pr[BLU][CIEX] - pr[RED][CIEX]) -
		  pr[RED][CIEX]*pr[BLU][CIEY] + pr[BLU][CIEX]*pr[RED][CIEY] ) ;
	if (CEQ(C_gD,0.))
		return(0);
	C_bD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[RED][CIEY] - pr[GRN][CIEY]) -
			  pr[WHT][CIEY]*(pr[RED][CIEX] - pr[GRN][CIEX]) +
		  pr[RED][CIEX]*pr[GRN][CIEY] - pr[GRN][CIEX]*pr[RED][CIEY] ) ;
	if (CEQ(C_bD,0.))
		return(0);
	mat[0][0] = (pr[GRN][CIEY] - pr[BLU][CIEY] -
			pr[BLU][CIEX]*pr[GRN][CIEY] +
			pr[BLU][CIEY]*pr[GRN][CIEX])/C_rD ;
	mat[0][1] = (pr[BLU][CIEX] - pr[GRN][CIEX] -
			pr[BLU][CIEX]*pr[GRN][CIEY] +
			pr[GRN][CIEX]*pr[BLU][CIEY])/C_rD ;
	mat[0][2] = (pr[GRN][CIEX]*pr[BLU][CIEY] -
			pr[BLU][CIEX]*pr[GRN][CIEY])/C_rD ;
	mat[1][0] = (pr[BLU][CIEY] - pr[RED][CIEY] -
			pr[BLU][CIEY]*pr[RED][CIEX] +
			pr[RED][CIEY]*pr[BLU][CIEX])/C_gD ;
	mat[1][1] = (pr[RED][CIEX] - pr[BLU][CIEX] -
			pr[RED][CIEX]*pr[BLU][CIEY] +
			pr[BLU][CIEX]*pr[RED][CIEY])/C_gD ;
	mat[1][2] = (pr[BLU][CIEX]*pr[RED][CIEY] -
			pr[RED][CIEX]*pr[BLU][CIEY])/C_gD ;
	mat[2][0] = (pr[RED][CIEY] - pr[GRN][CIEY] -
			pr[RED][CIEY]*pr[GRN][CIEX] +
			pr[GRN][CIEY]*pr[RED][CIEX])/C_bD ;
	mat[2][1] = (pr[GRN][CIEX] - pr[RED][CIEX] -
			pr[GRN][CIEX]*pr[RED][CIEY] +
			pr[RED][CIEX]*pr[GRN][CIEY])/C_bD ;
	mat[2][2] = (pr[RED][CIEX]*pr[GRN][CIEY] -
			pr[GRN][CIEX]*pr[RED][CIEY])/C_bD ;
	return(1);
}


int
comprgb2xyzmat(			/* compute conversion from RGB to CIE space */
COLORMAT  mat,
RGBPRIMS  pr
)
{
	double  C_rD, C_gD, C_bD, D;

	if (pr == stdprims) {	/* can use rgb2xyzmat */
		cpcolormat(mat, rgb2xyzmat);
		return(1);
	}
	if (pr == xyzprims) {	/* identity */
		memset(mat, 0, sizeof(COLORMAT));
		mat[0][0] = mat[1][1] = mat[2][2] = 1.f;
		return(1);
	}
	if (CEQ(pr[WHT][CIEX],0.) | CEQ(pr[WHT][CIEY],0.))
		return(0);
	C_rD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[GRN][CIEY] - pr[BLU][CIEY]) -
			  pr[WHT][CIEY]*(pr[GRN][CIEX] - pr[BLU][CIEX]) +
		  pr[GRN][CIEX]*pr[BLU][CIEY] - pr[BLU][CIEX]*pr[GRN][CIEY] ) ;
	C_gD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[BLU][CIEY] - pr[RED][CIEY]) -
			  pr[WHT][CIEY]*(pr[BLU][CIEX] - pr[RED][CIEX]) -
		  pr[RED][CIEX]*pr[BLU][CIEY] + pr[BLU][CIEX]*pr[RED][CIEY] ) ;
	C_bD = (1./pr[WHT][CIEY]) *
			( pr[WHT][CIEX]*(pr[RED][CIEY] - pr[GRN][CIEY]) -
			  pr[WHT][CIEY]*(pr[RED][CIEX] - pr[GRN][CIEX]) +
		  pr[RED][CIEX]*pr[GRN][CIEY] - pr[GRN][CIEX]*pr[RED][CIEY] ) ;
	D = pr[RED][CIEX]*(pr[GRN][CIEY] - pr[BLU][CIEY]) +
			pr[GRN][CIEX]*(pr[BLU][CIEY] - pr[RED][CIEY]) +
			pr[BLU][CIEX]*(pr[RED][CIEY] - pr[GRN][CIEY]) ;
	if (CEQ(D,0.))
		return(0);
	mat[0][0] = pr[RED][CIEX]*C_rD/D;
	mat[0][1] = pr[GRN][CIEX]*C_gD/D;
	mat[0][2] = pr[BLU][CIEX]*C_bD/D;
	mat[1][0] = pr[RED][CIEY]*C_rD/D;
	mat[1][1] = pr[GRN][CIEY]*C_gD/D;
	mat[1][2] = pr[BLU][CIEY]*C_bD/D;
	mat[2][0] = (1.-pr[RED][CIEX]-pr[RED][CIEY])*C_rD/D;
	mat[2][1] = (1.-pr[GRN][CIEX]-pr[GRN][CIEY])*C_gD/D;
	mat[2][2] = (1.-pr[BLU][CIEX]-pr[BLU][CIEY])*C_bD/D;
	return(1);
}


int
comprgb2rgbmat(			/* compute conversion from RGB1 to RGB2 */
COLORMAT  mat,
RGBPRIMS  pr1,
RGBPRIMS  pr2
)
{
	COLORMAT  pr1toxyz, xyztopr2;

	if (pr1 == pr2) {
		memset(mat, 0, sizeof(COLORMAT));
		mat[0][0] = mat[1][1] = mat[2][2] = 1.f;
		return(1);
	}
	if (!comprgb2xyzmat(pr1toxyz, pr1))
		return(0);
	if (!compxyz2rgbmat(xyztopr2, pr2))
		return(0);
				/* combine transforms */
	multcolormat(mat, pr1toxyz, xyztopr2);
	return(1);
}


int
compxyzWBmat(			/* CIE von Kries transform from wht1 to wht2 */
COLORMAT  mat,
float  wht1[2],
float  wht2[2]
)
{
	COLOR	cw1, cw2;
	if (XYEQ(wht1,wht2)) {
		memset(mat, 0, sizeof(COLORMAT));
		mat[0][0] = mat[1][1] = mat[2][2] = 1.f;
		return(1);
	}
	if (CEQ(wht1[CIEX],0.) | CEQ(wht1[CIEY],0.))
		return(0);
	cw1[RED] = wht1[CIEX]/wht1[CIEY];
	cw1[GRN] = 1.;
	cw1[BLU] = (1. - wht1[CIEX] - wht1[CIEY])/wht1[CIEY];
	colortrans(cw1, vkmat, cw1);
	if (CEQ(wht2[CIEX],0.) | CEQ(wht2[CIEY],0.))
		return(0);
	cw2[RED] = wht2[CIEX]/wht2[CIEY];
	cw2[GRN] = 1.;
	cw2[BLU] = (1. - wht2[CIEX] - wht2[CIEY])/wht2[CIEY];
	colortrans(cw2, vkmat, cw2);
	if (CEQ(cw1[RED],0.) | CEQ(cw1[GRN],0.) | CEQ(cw1[BLU],0.))
		return(0);
	mat[0][0] = cw2[RED]/cw1[RED];
	mat[1][1] = cw2[GRN]/cw1[GRN];
	mat[2][2] = cw2[BLU]/cw1[BLU];
	mat[0][1] = mat[0][2] = mat[1][0] =
	mat[1][2] = mat[2][0] = mat[2][1] = 0.0;
	multcolormat(mat, vkmat, mat);
	multcolormat(mat, mat, ivkmat);
	return(1);
}


int
compxyz2rgbWBmat(		/* von Kries conversion from CIE to RGB space */
COLORMAT  mat,
RGBPRIMS  pr
)
{
	COLORMAT	wbmat;

	if (!compxyz2rgbmat(mat, pr))
		return(0);
	if (XYEQ(pr[WHT],xyneu))
		return(1);
	if (!compxyzWBmat(wbmat, xyneu, pr[WHT]))
		return(0);
	multcolormat(mat, wbmat, mat);
	return(1);
}

int
comprgb2xyzWBmat(		/* von Kries conversion from RGB to CIE space */
COLORMAT  mat,
RGBPRIMS  pr
)
{
	COLORMAT	wbmat;
	
	if (!comprgb2xyzmat(mat, pr))
		return(0);
	if (XYEQ(pr[WHT],xyneu))
		return(1);
	if (!compxyzWBmat(wbmat, pr[WHT], xyneu))
		return(0);
	multcolormat(mat, mat, wbmat);
	return(1);
}

int
comprgb2rgbWBmat(		/* von Kries conversion from RGB1 to RGB2 */
COLORMAT  mat,
RGBPRIMS  pr1,
RGBPRIMS  pr2
)
{
	COLORMAT  pr1toxyz, xyztopr2, wbmat;

	if (pr1 == pr2) {
		memset(mat, 0, sizeof(COLORMAT));
		mat[0][0] = mat[1][1] = mat[2][2] = 1.f;
		return(1);
	}
	if (!comprgb2xyzmat(pr1toxyz, pr1))
		return(0);
	if (!compxyzWBmat(wbmat, pr1[WHT], pr2[WHT]))
		return(0);
	if (!compxyz2rgbmat(xyztopr2, pr2))
		return(0);
				/* combine transforms */
	multcolormat(mat, pr1toxyz, wbmat);
	multcolormat(mat, mat, xyztopr2);
	return(1);
}
