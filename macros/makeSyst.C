//--------------------------------------------------------------------------------------------------
// Compute the systematic uncertainties from the reduced trees
//
// Authors: L. Di Matteo                                                                  (Sep 2013)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TFile.h"
#include "TSystem.h"
#include "TSystem.h"
#include <TTree.h>
#include <THStack.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TVector2.h>
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitMonoPhoton/Core/MitGPTreeReduced.h"
#endif

// First we list the percentage errors as a function of jetEta and jetPt, which are given in Summer13_V1_MC_Uncertainty_AK5PF.txt. I have formatted them here
// into an array. The divisions by eta and pt for the rows and columns are given in the get_jet_uncertainty function.
// C. Ferko
const float JECErrors[40][44]=
{
	{0.0634,0.0602,0.0577,0.0554,0.0536,0.0521,0.0507,0.0495,0.0484,0.0473,0.0465,0.0457,0.045,0.0444,0.0439,0.0434,0.043,0.0426,0.0423,0.042,0.0418,0.0417,0.0415,0.0414,0.0414,0.0412,0.0412,0.0411,0.0412,0.0413,0.0414,0.0415,0.0416,0.0417,0.0419,0.042,0.0421,0.0422,0.0423,0.0424,0.0425,0.0426,0.0427,0.0429},
	{0.0629,0.0598,0.0574,0.0552,0.0534,0.052,0.0506,0.0494,0.0483,0.0473,0.0464,0.045,0.0423,0.0402,0.0392,0.0387,0.0381,0.0377,0.0374,0.0371,0.0369,0.0367,0.0365,0.0364,0.0363,0.0362,0.0361,0.036,0.0361,0.0362,0.0364,0.0365,0.0366,0.0367,0.0369,0.037,0.0371,0.0373,0.0374,0.0375,0.0376,0.0377,0.0379,0.038},
	{0.0635,0.0603,0.0577,0.0554,0.0536,0.0521,0.0507,0.0495,0.0484,0.0473,0.0464,0.045,0.0424,0.0402,0.0387,0.0381,0.0384,0.0394,0.0391,0.0388,0.0386,0.0384,0.0383,0.0381,0.0381,0.0379,0.0378,0.0378,0.0379,0.038,0.0381,0.0382,0.0384,0.0385,0.0386,0.0387,0.0389,0.039,0.0391,0.0392,0.0393,0.0394,0.0396,0.0397},
	{0.0716,0.0662,0.0618,0.0582,0.0556,0.0537,0.0519,0.0504,0.0491,0.0479,0.0468,0.0453,0.0426,0.0404,0.0389,0.0382,0.0385,0.0399,0.0399,0.0396,0.0394,0.0392,0.0391,0.039,0.0389,0.0388,0.0387,0.0386,0.0387,0.0388,0.0389,0.039,0.0392,0.0393,0.0394,0.0395,0.0396,0.0398,0.0399,0.04,0.0401,0.0402,0.0403,0.0405},
	{0.071,0.0635,0.0571,0.0518,0.048,0.0451,0.0425,0.0402,0.0383,0.0364,0.0349,0.0336,0.0323,0.0312,0.0302,0.0294,0.0286,0.0279,0.0274,0.027,0.0267,0.0265,0.0263,0.0261,0.026,0.0258,0.0257,0.0256,0.0257,0.0259,0.026,0.0262,0.0264,0.0266,0.0267,0.0269,0.0271,0.0273,0.0274,0.0276,0.0277,0.0279,0.0281,0.0283},
	{0.0542,0.0504,0.0473,0.0445,0.0422,0.0403,0.0385,0.0369,0.0354,0.0339,0.0327,0.0316,0.0304,0.0293,0.0283,0.0274,0.0266,0.0259,0.0253,0.0249,0.0246,0.0244,0.0243,0.0242,0.0241,0.0241,0.0239,0.0239,0.024,0.0242,0.0244,0.0246,0.0248,0.0249,0.0251,0.0253,0.0255,0.0257,0.0259,0.026,0.0262,0.0264,0.0266,0.0268},
	{0.0627,0.0563,0.0509,0.0464,0.043,0.0405,0.0381,0.036,0.0341,0.0323,0.0308,0.0293,0.0274,0.0255,0.0238,0.0223,0.0207,0.0193,0.0181,0.0172,0.0166,0.0162,0.0159,0.0158,0.0159,0.0161,0.0168,0.0172,0.0174,0.0176,0.0179,0.0181,0.0184,0.0186,0.0189,0.0192,0.0194,0.0197,0.0198,0.02,0.0203,0.0205,0.0208,0.021},
	{0.0774,0.0683,0.0602,0.0535,0.0488,0.0453,0.0422,0.0395,0.0371,0.0349,0.0331,0.0313,0.0291,0.0269,0.025,0.0232,0.0213,0.0196,0.0183,0.0173,0.0165,0.0159,0.0156,0.0154,0.0154,0.0156,0.0163,0.017,0.0176,0.018,0.0183,0.0185,0.0188,0.019,0.0193,0.0195,0.0197,0.02,0.0202,0.0204,0.0206,0.0209,0.0211,0.0213},
	{0.0613,0.0551,0.0495,0.0448,0.0412,0.0384,0.0357,0.0333,0.0311,0.029,0.0272,0.0256,0.0239,0.0223,0.0208,0.0195,0.018,0.0167,0.0156,0.0147,0.014,0.0134,0.0129,0.0125,0.0122,0.0119,0.0117,0.0117,0.0123,0.0128,0.0131,0.0135,0.0138,0.0141,0.0145,0.0148,0.0151,0.0154,0.0157,0.0159,0.0162,0.0165,0.0169,0.0172},
	{0.0554,0.0503,0.0459,0.042,0.039,0.0366,0.0342,0.032,0.0301,0.0282,0.0265,0.025,0.0233,0.0218,0.0204,0.0191,0.0177,0.0164,0.0154,0.0145,0.0138,0.0132,0.0127,0.0123,0.0119,0.0116,0.0114,0.0114,0.0119,0.0126,0.0132,0.0137,0.0141,0.0144,0.0147,0.0151,0.0154,0.0157,0.0159,0.0162,0.0165,0.0168,0.0171,0.0174},
	{0.0501,0.0458,0.0423,0.039,0.0364,0.0342,0.032,0.03,0.0281,0.0262,0.0246,0.0231,0.0217,0.0203,0.019,0.0178,0.0166,0.0155,0.0145,0.0137,0.013,0.0125,0.012,0.0116,0.0113,0.0109,0.0105,0.0104,0.0108,0.0112,0.0117,0.0122,0.0125,0.0129,0.0133,0.0136,0.014,0.0143,0.0146,0.0149,0.0152,0.0155,0.0158,0.0162},
	{0.0496,0.0455,0.0422,0.039,0.0365,0.0343,0.0322,0.0302,0.0284,0.0266,0.025,0.0235,0.022,0.0206,0.0193,0.0181,0.0169,0.0157,0.0147,0.0138,0.0131,0.0126,0.0121,0.0116,0.0113,0.0109,0.0106,0.0104,0.0109,0.0114,0.0119,0.0124,0.013,0.0135,0.0138,0.0142,0.0145,0.0149,0.0151,0.0153,0.0157,0.016,0.0163,0.0166},
	{0.0492,0.0451,0.0417,0.0385,0.0359,0.0338,0.0316,0.0296,0.0277,0.0259,0.0242,0.0228,0.0214,0.02,0.0188,0.0176,0.0164,0.0153,0.0144,0.0136,0.013,0.0124,0.012,0.0116,0.0112,0.0108,0.0105,0.0103,0.0107,0.0111,0.0115,0.012,0.0124,0.0128,0.0132,0.0136,0.0139,0.0143,0.0145,0.0148,0.0151,0.0154,0.0158,0.0161},
	{0.0488,0.0447,0.0413,0.0382,0.0355,0.0333,0.0312,0.0291,0.0272,0.0253,0.0236,0.0221,0.0207,0.0194,0.0182,0.017,0.0158,0.0147,0.0138,0.013,0.0124,0.0118,0.0114,0.011,0.0106,0.0102,0.0098,0.0096,0.01,0.0104,0.0108,0.0112,0.0116,0.012,0.0124,0.0128,0.0132,0.0136,0.0138,0.0141,0.0145,0.0148,0.0152,0.0155},
	{0.0491,0.045,0.0416,0.0384,0.0358,0.0336,0.0315,0.0295,0.0276,0.0257,0.0241,0.0226,0.0212,0.0198,0.0185,0.0173,0.0161,0.015,0.014,0.0132,0.0125,0.0119,0.0114,0.011,0.0106,0.0102,0.0098,0.0096,0.01,0.0105,0.0109,0.0114,0.0118,0.0123,0.0127,0.0132,0.0135,0.0139,0.0142,0.0144,0.0148,0.0151,0.0154,0.0158},
	{0.0496,0.0455,0.0422,0.0391,0.0365,0.0344,0.0323,0.0303,0.0285,0.0267,0.0251,0.0237,0.0222,0.0207,0.0193,0.0181,0.0167,0.0155,0.0144,0.0135,0.0128,0.0121,0.0116,0.0111,0.0107,0.0102,0.0098,0.0096,0.0101,0.0106,0.0111,0.0117,0.0122,0.0128,0.0133,0.0138,0.0143,0.0146,0.0148,0.0151,0.0154,0.0158,0.0161,0.0164},
	{0.05,0.046,0.0427,0.0396,0.037,0.0349,0.0329,0.0309,0.0291,0.0274,0.0258,0.0244,0.0228,0.0213,0.0199,0.0186,0.0172,0.0159,0.0147,0.0138,0.013,0.0123,0.0117,0.0112,0.0107,0.0102,0.0098,0.0096,0.01,0.0106,0.0112,0.0118,0.0123,0.0129,0.0135,0.0141,0.0146,0.0149,0.0151,0.0154,0.0157,0.016,0.0163,0.0167},
	{0.0493,0.0452,0.0419,0.0388,0.0362,0.034,0.0319,0.0299,0.0281,0.0263,0.0246,0.0232,0.0217,0.0203,0.019,0.0178,0.0165,0.0153,0.0142,0.0134,0.0127,0.012,0.0115,0.0111,0.0107,0.0102,0.0098,0.0096,0.01,0.0105,0.0109,0.0114,0.0119,0.0124,0.0129,0.0134,0.0138,0.0141,0.0144,0.0147,0.015,0.0153,0.0157,0.016},
	{0.0496,0.0456,0.0423,0.0392,0.0367,0.0346,0.0325,0.0306,0.0288,0.027,0.0254,0.024,0.0224,0.021,0.0196,0.0183,0.017,0.0157,0.0146,0.0137,0.0129,0.0122,0.0116,0.0112,0.0107,0.0102,0.0098,0.0096,0.01,0.0105,0.0111,0.0116,0.0122,0.0127,0.0133,0.0138,0.0142,0.0146,0.0148,0.0151,0.0154,0.0157,0.0161,0.0164},
	{0.0488,0.0447,0.0414,0.0382,0.0357,0.0335,0.0313,0.0293,0.0274,0.0256,0.0239,0.0224,0.021,0.0196,0.0184,0.0172,0.016,0.0149,0.0139,0.0131,0.0125,0.0119,0.0114,0.011,0.0106,0.0102,0.0098,0.0096,0.01,0.0104,0.0108,0.0113,0.0117,0.0121,0.0126,0.013,0.0134,0.0137,0.014,0.0143,0.0146,0.015,0.0153,0.0156},
	{0.0488,0.0447,0.0414,0.0382,0.0357,0.0335,0.0313,0.0293,0.0274,0.0256,0.0239,0.0224,0.021,0.0196,0.0184,0.0172,0.016,0.0149,0.0139,0.0131,0.0125,0.0119,0.0114,0.011,0.0106,0.0102,0.0098,0.0096,0.01,0.0104,0.0108,0.0113,0.0117,0.0121,0.0126,0.013,0.0134,0.0137,0.014,0.0143,0.0146,0.015,0.0153,0.0156},
	{0.0496,0.0456,0.0423,0.0392,0.0367,0.0346,0.0325,0.0306,0.0288,0.027,0.0254,0.024,0.0224,0.021,0.0196,0.0183,0.017,0.0157,0.0146,0.0137,0.0129,0.0122,0.0116,0.0112,0.0107,0.0102,0.0098,0.0096,0.01,0.0105,0.0111,0.0116,0.0122,0.0127,0.0133,0.0138,0.0142,0.0146,0.0148,0.0151,0.0154,0.0157,0.0161,0.0164},
	{0.0493,0.0452,0.0419,0.0388,0.0362,0.034,0.0319,0.0299,0.0281,0.0263,0.0246,0.0232,0.0217,0.0203,0.019,0.0178,0.0165,0.0153,0.0142,0.0134,0.0127,0.012,0.0115,0.0111,0.0107,0.0102,0.0098,0.0096,0.01,0.0105,0.0109,0.0114,0.0119,0.0124,0.0129,0.0134,0.0138,0.0141,0.0144,0.0147,0.015,0.0153,0.0157,0.016},
	{0.05,0.046,0.0427,0.0396,0.037,0.0349,0.0329,0.0309,0.0291,0.0274,0.0258,0.0244,0.0228,0.0213,0.0199,0.0186,0.0172,0.0159,0.0147,0.0138,0.013,0.0123,0.0117,0.0112,0.0107,0.0102,0.0098,0.0096,0.01,0.0106,0.0112,0.0118,0.0123,0.0129,0.0135,0.0141,0.0146,0.0149,0.0151,0.0154,0.0157,0.016,0.0163,0.0167},
	{0.0496,0.0455,0.0422,0.0391,0.0365,0.0344,0.0323,0.0303,0.0285,0.0267,0.0251,0.0237,0.0222,0.0207,0.0193,0.0181,0.0167,0.0155,0.0144,0.0135,0.0128,0.0121,0.0116,0.0111,0.0107,0.0102,0.0098,0.0096,0.0101,0.0106,0.0111,0.0117,0.0122,0.0128,0.0133,0.0138,0.0143,0.0146,0.0148,0.0151,0.0154,0.0158,0.0161,0.0164},
	{0.0491,0.045,0.0416,0.0384,0.0358,0.0336,0.0315,0.0295,0.0276,0.0257,0.0241,0.0226,0.0212,0.0198,0.0185,0.0173,0.0161,0.015,0.014,0.0132,0.0125,0.0119,0.0114,0.011,0.0106,0.0102,0.0098,0.0096,0.01,0.0105,0.0109,0.0114,0.0118,0.0123,0.0127,0.0132,0.0135,0.0139,0.0142,0.0144,0.0148,0.0151,0.0154,0.0158},
	{0.0488,0.0447,0.0413,0.0382,0.0355,0.0333,0.0312,0.0291,0.0272,0.0253,0.0236,0.0221,0.0207,0.0194,0.0182,0.017,0.0158,0.0147,0.0138,0.013,0.0124,0.0118,0.0114,0.011,0.0106,0.0102,0.0098,0.0096,0.01,0.0104,0.0108,0.0112,0.0116,0.012,0.0124,0.0128,0.0132,0.0136,0.0138,0.0141,0.0145,0.0148,0.0152,0.0155},
	{0.0492,0.0451,0.0417,0.0385,0.0359,0.0338,0.0316,0.0296,0.0277,0.0259,0.0242,0.0228,0.0214,0.02,0.0188,0.0176,0.0164,0.0153,0.0144,0.0136,0.013,0.0124,0.012,0.0116,0.0112,0.0108,0.0105,0.0103,0.0107,0.0111,0.0115,0.012,0.0124,0.0128,0.0132,0.0136,0.0139,0.0143,0.0145,0.0148,0.0151,0.0154,0.0158,0.0161},
	{0.0496,0.0455,0.0422,0.039,0.0365,0.0343,0.0322,0.0302,0.0284,0.0266,0.025,0.0235,0.022,0.0206,0.0193,0.0181,0.0169,0.0157,0.0147,0.0138,0.0131,0.0126,0.0121,0.0116,0.0113,0.0109,0.0106,0.0104,0.0109,0.0114,0.0119,0.0124,0.013,0.0135,0.0138,0.0142,0.0145,0.0149,0.0151,0.0153,0.0157,0.016,0.0163,0.0166},
	{0.0501,0.0458,0.0423,0.039,0.0364,0.0342,0.032,0.03,0.0281,0.0262,0.0246,0.0231,0.0217,0.0203,0.019,0.0178,0.0166,0.0155,0.0145,0.0137,0.013,0.0125,0.012,0.0116,0.0113,0.0109,0.0105,0.0104,0.0108,0.0112,0.0117,0.0122,0.0125,0.0129,0.0133,0.0136,0.014,0.0143,0.0146,0.0149,0.0152,0.0155,0.0158,0.0162},
	{0.0554,0.0503,0.0459,0.042,0.039,0.0366,0.0342,0.032,0.0301,0.0282,0.0265,0.025,0.0233,0.0218,0.0204,0.0191,0.0177,0.0164,0.0154,0.0145,0.0138,0.0132,0.0127,0.0123,0.0119,0.0116,0.0114,0.0114,0.0119,0.0126,0.0132,0.0137,0.0141,0.0144,0.0147,0.0151,0.0154,0.0157,0.0159,0.0162,0.0165,0.0168,0.0171,0.0174},
	{0.0613,0.0551,0.0495,0.0448,0.0412,0.0384,0.0357,0.0333,0.0311,0.029,0.0272,0.0256,0.0239,0.0223,0.0208,0.0195,0.018,0.0167,0.0156,0.0147,0.014,0.0134,0.0129,0.0125,0.0122,0.0119,0.0117,0.0117,0.0123,0.0128,0.0131,0.0135,0.0138,0.0141,0.0145,0.0148,0.0151,0.0154,0.0157,0.0159,0.0162,0.0165,0.0169,0.0172},
	{0.0774,0.0683,0.0602,0.0535,0.0488,0.0453,0.0422,0.0395,0.0371,0.0349,0.0331,0.0313,0.0291,0.0269,0.025,0.0232,0.0213,0.0196,0.0183,0.0173,0.0165,0.0159,0.0156,0.0154,0.0154,0.0156,0.0163,0.017,0.0176,0.018,0.0183,0.0185,0.0188,0.019,0.0193,0.0195,0.0197,0.02,0.0202,0.0204,0.0206,0.0209,0.0211,0.0213},
	{0.0627,0.0563,0.0509,0.0464,0.043,0.0405,0.0381,0.036,0.0341,0.0323,0.0308,0.0293,0.0274,0.0255,0.0238,0.0223,0.0207,0.0193,0.0181,0.0172,0.0166,0.0162,0.0159,0.0158,0.0159,0.0161,0.0168,0.0172,0.0174,0.0176,0.0179,0.0181,0.0184,0.0186,0.0189,0.0192,0.0194,0.0197,0.0198,0.02,0.0203,0.0205,0.0208,0.021},
	{0.0542,0.0504,0.0473,0.0445,0.0422,0.0403,0.0385,0.0369,0.0354,0.0339,0.0327,0.0316,0.0304,0.0293,0.0283,0.0274,0.0266,0.0259,0.0253,0.0249,0.0246,0.0244,0.0243,0.0242,0.0241,0.0241,0.0239,0.0239,0.024,0.0242,0.0244,0.0246,0.0248,0.0249,0.0251,0.0253,0.0255,0.0257,0.0259,0.026,0.0262,0.0264,0.0266,0.0268},
	{0.071,0.0635,0.0571,0.0518,0.048,0.0451,0.0425,0.0402,0.0383,0.0364,0.0349,0.0336,0.0323,0.0312,0.0302,0.0294,0.0286,0.0279,0.0274,0.027,0.0267,0.0265,0.0263,0.0261,0.026,0.0258,0.0257,0.0256,0.0257,0.0259,0.026,0.0262,0.0264,0.0266,0.0267,0.0269,0.0271,0.0273,0.0274,0.0276,0.0277,0.0279,0.0281,0.0283},
	{0.0716,0.0662,0.0618,0.0582,0.0556,0.0537,0.0519,0.0504,0.0491,0.0479,0.0468,0.0453,0.0426,0.0404,0.0389,0.0382,0.0385,0.0399,0.0399,0.0396,0.0394,0.0392,0.0391,0.039,0.0389,0.0388,0.0387,0.0386,0.0387,0.0388,0.0389,0.039,0.0392,0.0393,0.0394,0.0395,0.0396,0.0398,0.0399,0.04,0.0401,0.0402,0.0403,0.0405},
	{0.0635,0.0603,0.0577,0.0554,0.0536,0.0521,0.0507,0.0495,0.0484,0.0473,0.0464,0.045,0.0424,0.0402,0.0387,0.0381,0.0384,0.0394,0.0391,0.0388,0.0386,0.0384,0.0383,0.0381,0.0381,0.0379,0.0378,0.0378,0.0379,0.038,0.0381,0.0382,0.0384,0.0385,0.0386,0.0387,0.0389,0.039,0.0391,0.0392,0.0393,0.0394,0.0396,0.0397},
	{0.0629,0.0598,0.0574,0.0552,0.0534,0.052,0.0506,0.0494,0.0483,0.0473,0.0464,0.045,0.0423,0.0402,0.0392,0.0387,0.0381,0.0377,0.0374,0.0371,0.0369,0.0367,0.0365,0.0364,0.0363,0.0362,0.0361,0.036,0.0361,0.0362,0.0364,0.0365,0.0366,0.0367,0.0369,0.037,0.0371,0.0373,0.0374,0.0375,0.0376,0.0377,0.0379,0.038},
	{0.0634,0.0602,0.0577,0.0554,0.0536,0.0521,0.0507,0.0495,0.0484,0.0473,0.0465,0.0457,0.045,0.0444,0.0439,0.0434,0.043,0.0426,0.0423,0.042,0.0418,0.0417,0.0415,0.0414,0.0414,0.0412,0.0412,0.0411,0.0412,0.0413,0.0414,0.0415,0.0416,0.0417,0.0419,0.042,0.0421,0.0422,0.0423,0.0424,0.0425,0.0426,0.0427,0.0429}
};


using namespace std;
using namespace mithep;

TString getEnv(const char* name);

float getJECError(float eta, float pt);

void  getPhotonSyst(float et, TRandom3& randGen, float& etUp, float& etDown);

void  getMetSyst(TVector2& met, TVector2& metUp, TVector2& metDown, TString mode, 
                 TVector2& singleVec, TVector2& singleVecUp, TVector2& singleVecDown,
                 MitGPTreeReduced& myTree);

//==================================================================================================
void makeSyst()
{
  // read all environment variables
  TString here   = getEnv("PWD");
  TString home   = getEnv("HOME");
  TString mitHgg = getEnv("MIT_MONOPH_DIR");
  TString hstDir = getEnv("MIT_ANA_HIST");
  TString anaCfg = getEnv("MIT_PLOT_CFG");
  TString prdCfg = getEnv("MIT_PROD_CFG");

  //DEBUG
  cout << JECErrors[0][0] << endl;
  
  // define samples
  TaskSamples* samples = new TaskSamples(prdCfg.Data(),hstDir.Data());
  samples->SetNameTxt(anaCfg.Data());
  samples->ReadFile((mitHgg + TString("/config")).Data());
  vector<const Sample*> listOfSamples;
  for (UInt_t iSample=0; iSample < *samples->NSamples(); iSample++) listOfSamples.push_back(samples->GetSample(iSample));  
  vector<const Sample*> listOfDatasets;
  for (UInt_t iSample=0; iSample < *samples->NDataSamples(); iSample++) listOfDatasets.push_back(samples->GetDataSample(iSample));
  // define infolder
  TString inFileName = here + "/monoph-2013-July9_reduced.root";
  std::cout << "inFileName " << inFileName << std::endl;
  
  // prepare the random number generator used for the smearings
  TRandom3 theRandGenerator;

  // prepare the systematic names : keep an eye to the ordering !
  vector<TString>  v_systNames;
  v_systNames.push_back("Nominal");
  v_systNames.push_back("PuUp");
  v_systNames.push_back("PuDown");
  v_systNames.push_back("PhotonRecoUp");
  v_systNames.push_back("PhotonRecoDown");
  v_systNames.push_back("LeptonRecoUp");
  v_systNames.push_back("LeptonRecoDown");
  v_systNames.push_back("JetRecoUp");
  v_systNames.push_back("JetRecoDown");
  const int nSyst = v_systNames.size();

  // prepare the systematic matrix to store all the yields
  vector<vector<float> > m_systYields; 

  // prepare the shortened sample list to organize the output in an easy way
  vector<TString>  v_groupSampleName;
  
  // loop through the samples and compute the yields for different conditions
  int   theGroupSampleCounter = 0;
  vector<float> v_thisYield;
  for (UInt_t iSample=0; iSample < listOfSamples.size(); iSample++) {
  //for (UInt_t iSample=0; iSample < 10; iSample++) {

    //Say which sample we are processing
    cout << "Processing sample " << *(listOfSamples.at(iSample)->Name()) << endl;

    //determine if the sample is first of the series
    bool isFirstOfSerie = (*listOfSamples.at(iSample)->Legend()).CompareTo(" ");
    bool isLastOfSerie = false;
    if (iSample == listOfSamples.size() - 1) isLastOfSerie = true;
    if (iSample < listOfSamples.size() - 1 && (*listOfSamples.at(iSample+1)->Legend()).CompareTo(" ") != 0) isLastOfSerie = true;
    
    //if sample first of the list insert the shortname in the appropriate vector and refresh the yield
    if (isFirstOfSerie) {
      v_groupSampleName.push_back(*(listOfSamples.at(iSample)->Legend()));
      for (int iSyst=0; iSyst < nSyst; iSyst++) 
        v_thisYield.push_back(0.);
    }

    //get the sample scale 
    float thisScale = *(listOfSamples.at(iSample)->Scale());
 
    //get a nice GPTree for this sample 
    MitGPTreeReduced thistree;
    thistree.LoadTree(inFileName,*(listOfSamples.at(iSample)->Name()));
    thistree.InitTree();

    // Loop on the tree and make all the computations
    Int_t nEntries = thistree.tree_->GetEntries();
    for ( Int_t iEntry = 0; iEntry < nEntries; iEntry++ ) {
      thistree.tree_-> GetEntry(iEntry);
      
      //Step1: get the relevant observables for the selection
      float nominalPhotonEt = thistree.phoEt_;
      float nominalMet = thistree.met_;
      float nominalJetPt = thistree.jet1Pt_;
      float nominalLepPt = thistree.lep1Pt_;
      TVector2 nominalMetVec;
      nominalMetVec.SetMagPhi(thistree.met_,thistree.metPhi_);

      //Step2: get the nominal yield for the standard selection and PU systematics
      if ( nominalPhotonEt >= 160 && nominalLepPt < 10. && nominalJetPt < 100. && nominalMet >= 140 && thistree.ncosmics_ == 0 ) {
        //Nominal
        v_thisYield[0] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_;
        //PuUp
        v_thisYield[1] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weightup_;
        //PuDown
        v_thisYield[2] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weightdo_;
      }
      
      //Step3: propagate the photon uncertanties
      float upPhotonEt, downPhotonEt;
      getPhotonSyst(nominalPhotonEt, theRandGenerator, upPhotonEt, downPhotonEt);
      TVector2 nominalPhotonVec, photonSystUpPhotonVec, photonSystDownPhotonVec;
      nominalPhotonVec.SetMagPhi(nominalPhotonEt, thistree.phoPhi_);
      photonSystUpPhotonVec.SetMagPhi(upPhotonEt, thistree.phoPhi_);
      photonSystDownPhotonVec.SetMagPhi(downPhotonEt, thistree.phoPhi_);
      TVector2 photonSystUpMetVec,photonSystDownMetVec;
      
      getMetSyst(nominalMetVec, photonSystUpMetVec, photonSystDownMetVec, "single", 
                 nominalPhotonVec, photonSystUpPhotonVec, photonSystDownPhotonVec, 
                 thistree);
                
      //Compute the photon eff systematics: multiplicative efficiencies error propagation
      float photonEffSystUp = (1 + 0.008*sqrt(thistree.nphotons_));
      float photonEffSystDown = (1 - 0.008*sqrt(thistree.nphotons_));

      if ( upPhotonEt >= 160 && nominalLepPt < 10. && nominalJetPt < 100. && photonSystUpMetVec.Mod() >= 140 && thistree.ncosmics_ == 0 )
        v_thisYield[3] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_*photonEffSystUp;
      if ( downPhotonEt >= 160 && nominalLepPt < 10. && nominalJetPt < 100. && photonSystDownMetVec.Mod() >= 140 && thistree.ncosmics_ == 0 )
        v_thisYield[4] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_*photonEffSystDown;
      
      //Step4: propagate the lepton uncertainties
      if (nominalLepPt < 10.) {
        v_thisYield[5] = v_thisYield[0];
        v_thisYield[6] = v_thisYield[0];
      }
      else {
        TVector2 dummy, leptonSystUpLeptonVec, leptonSystDownLeptonVec;
        TVector2 leptonSystUpMetVec,leptonSystDownMetVec;
        getMetSyst(nominalMetVec, leptonSystUpMetVec, leptonSystDownMetVec, "leptons", 
                   dummy, leptonSystUpLeptonVec, leptonSystDownLeptonVec,
                   thistree);
        //Compute the lepton eff systematics: multiplicative efficiencies error propagation
        float leptonEffSystUp = (1 + 0.01*sqrt(thistree.nlep_));
        float leptonEffSystDown = (1 - 0.01*sqrt(thistree.nlep_));

        if ( nominalPhotonEt >= 160 && leptonSystUpLeptonVec.Mod() < 10. && nominalJetPt < 100. && leptonSystUpMetVec.Mod() >= 140 && thistree.ncosmics_ == 0 )
          v_thisYield[5] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_*leptonEffSystUp;
        if ( nominalPhotonEt >= 160 && leptonSystDownLeptonVec.Mod() < 10. && nominalJetPt < 100. && leptonSystDownMetVec.Mod() >= 140 && thistree.ncosmics_ == 0 )
          v_thisYield[6] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_*leptonEffSystDown;
      }
      
      //Step5: propagate the jets uncertainties
      TVector2 dummy, jetSystUpJetVec, jetSystDownJetVec;
      TVector2 jetSystUpMetVec,jetSystDownMetVec;
      jetSystUpMetVec.SetMagPhi(thistree.jet1Pt_,thistree.jet1Phi_);
      jetSystDownMetVec.SetMagPhi(thistree.jet1Pt_,thistree.jet1Phi_);
      getMetSyst(nominalMetVec, jetSystUpMetVec, jetSystDownMetVec, "jets", 
                 dummy, jetSystUpJetVec, jetSystDownJetVec,
                 thistree);
      //DEBUG
      //cout << nominalJetPt << " " << jetSystUpJetVec.Mod() << " " << jetSystDownJetVec.Mod() << endl;
      //cout << nominalMet << " " << jetSystUpMetVec.Mod() << " " << jetSystDownMetVec.Mod() << endl;
      if ( nominalPhotonEt >= 160 && nominalLepPt < 10. && jetSystUpJetVec.Mod() < 100. && jetSystUpMetVec.Mod() >= 140 && thistree.ncosmics_ == 0 )
        v_thisYield[7] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_;
      if ( nominalPhotonEt >= 160 && nominalLepPt < 10. && jetSystDownJetVec.Mod() < 100. && jetSystDownMetVec.Mod() >= 140 && thistree.ncosmics_ == 0 )
        v_thisYield[8] += thisScale*thistree.evt_weight_*thistree.kf_weight_*thistree.pu_weight_;
            
      
    }
    
    //DEBUG
    //cout << "this yield is " <<  v_thisYield[0] << " " << v_thisYield[1] << " " << v_thisYield[2] << endl;
    
    //add the yields to the syst matrix if this sample is last of the series:
    //either last sample or ~ sample followed by non ~ sample
    if (isLastOfSerie) {
       //Add me to the matrix!
       //clear the yield vector
       m_systYields.push_back(v_thisYield);
       v_thisYield.clear();
       theGroupSampleCounter++;
    }
    
  }//end loop on samples

  //test the matrix
  for (int iGroupSample = 0; iGroupSample < m_systYields.size(); iGroupSample++) {
    cout << v_groupSampleName.at(iGroupSample);
    for (int iSyst = 0; iSyst < nSyst; iSyst++) cout << " " << m_systYields.at(iGroupSample)[iSyst];
    cout << endl;
  }
        
  return;
}

TString getEnv(const char* name)
{
  if (! gSystem->Getenv(name)) {
    printf(" Environment variable: %s  not defined. EXIT!\n",name);
    return TString("");
  } 
  return TString(gSystem->Getenv(name));  
}

// This function shows which values of jetEta and jetPt the rows and columns of the above double refer to. 
// For example, the first row is for -5.4<jetPt<-5 and the first column is for 9<jetPt<11.
// C.Ferko
float getJECError(float eta, float pt)
{
	int index1;
	int index2;
	
	if (-5.4<=eta && eta<=-5) index1=0;
	if (-5<eta && eta<=-4.4) index1=1;
	if (-4.4<eta && eta<=-4) index1=2;
	if (-4<eta && eta<=-3.5) index1=3;
	if (-3.5<eta && eta<=-3) index1=4;
	if (-3<eta && eta<=-2.8) index1=5;
	if (-2.8<eta && eta<=-2.6) index1=6;
	if (-2.6<eta && eta<=-2.4) index1=7;
	if (-2.4<eta && eta<=-2.2) index1=8;
	if (-2.2<eta && eta<=-2) index1=9;
	if (-2<eta && eta<=-1.8) index1=10;
	if (-1.8<eta && eta<=-1.6) index1=11;
	if (-1.6<eta && eta<=-1.4) index1=12;
	if (-1.4<eta && eta<=-1.2) index1=13;
	if (-1.2<eta && eta<=-1) index1=14;
	if (-1<eta && eta<=-0.8) index1=15;
	if (-0.8<eta && eta<=-0.6) index1=16;
	if (-0.6<eta && eta<=-0.4) index1=17;
	if (-0.4<eta && eta<=-0.2) index1=18;
	if (-0.2<eta && eta<=0) index1=19;
	if (0<eta && eta<=0.2) index1=20;
	if (0.2<eta && eta<=0.4) index1=21;
	if (0.4<eta && eta<=0.6) index1=22;
	if (0.6<eta && eta<=0.8) index1=23;
	if (0.8<eta && eta<=1) index1=24;
	if (1<eta && eta<=1.2) index1=25;
	if (1.2<eta && eta<=1.4) index1=26;
	if (1.4<eta && eta<=1.6) index1=27;
	if (1.6<eta && eta<=1.8) index1=28;
	if (1.8<eta && eta<=2) index1=29;
	if (2<eta && eta<=2.2) index1=30;
	if (2.2<eta && eta<=2.4) index1=31;
	if (2.4<eta && eta<=2.6) index1=32;
	if (2.6<eta && eta<=2.8) index1=33;
	if (2.8<eta && eta<=3) index1=34;
	if (3<eta && eta<=3.5) index1=35;
	if (3.5<eta && eta<=4) index1=36;
	if (4<eta && eta<=4.4) index1=37;
	if (4.4<eta && eta<=5) index1=38;
	if (5<eta && eta<=5.4) index1=39;
	
	if (9<pt && pt<=11) index2=0;	
	if (11<pt && pt<=13.5) index2=1;
	if (13.5<pt && pt<=16.5) index2=2;
	if (16.5<pt && pt<=19.5) index2=3;
	if (19.5<pt && pt<=22.5) index2=4;
	if (22.5<pt && pt<=26) index2=5;
	if (26<pt && pt<=30) index2=6;
	if (30<pt && pt<=34.5) index2=7;
	if (34.5<pt && pt<=40) index2=8;
	if (40<pt && pt<=46) index2=9;
	if (46<pt && pt<=52.5) index2=10;
	if (52.5<pt && pt<=60) index2=11;
	if (60<pt && pt<=69) index2=12;
	if (69<pt && pt<=79) index2=13;
	if (79<pt && pt<=90.5) index2=14;
	if (90.5<pt && pt<=105.5) index2=15;
	if (105.5<pt && pt<=123.5) index2=16;
	if (123.5<pt && pt<=143) index2=17;

	if (143<pt && pt<=163.5) index2=18;
	if (163.5<pt && pt<=185) index2=19;
	if (185<pt && pt<=208) index2=20;
	if (208<pt && pt<=232.5) index2=21;
	if (232.5<pt && pt<=258.5) index2=22;
	if (258.5<pt && pt<=286) index2=23;
	if (286<pt && pt<=331) index2=24;
	if (331<pt && pt<=396) index2=25;
	if (396<pt && pt<=468.5) index2=26;
	if (468.5<pt && pt<=549.5) index2=27;
	if (549.5<pt && pt<=639) index2=28;
	if (639<pt && pt<=738) index2=29;
	if (738<pt && pt<=847.5) index2=30;
	if (847.5<pt && pt<=968.5) index2=31;
	if (968.5<pt && pt<=1102) index2=32;
	if (1102<pt && pt<=1249.5) index2=33;
	if (1249.5<pt && pt<=1412) index2=34;
	if (1412<pt && pt<=1590.5) index2=35;
	if (1590.5<pt && pt<=1787) index2=36;
	if (1787<pt && pt<=1945) index2=37;
	if (1945<pt && pt<=2119) index2=38;
	if (2119<pt && pt<=2369) index2=39;
	if (2369<pt && pt<=2643.5) index2=40;
	if (2643.5<pt && pt<=2945) index2=41;
	if (2945<pt && pt<=3276.5) index2=42;
	if (pt>3276.5) index2=43;
	
	return JECErrors[index1][index2];
}

void getPhotonSyst(float et, TRandom3& randGen, float& etUp, float& etDown) {
  
  //take the worst syst for each category
  float smear = 0.0025;
  float scale = 0.002;
  float etSmear = randGen.Gaus(et,smear*et);
  etUp = et*(1.+scale);
  etDown = et*(1.-scale);
  
  //convolve the two uncertainties
  if (etSmear > etUp) etUp = etSmear;
  if (etSmear < etDown) etDown = etSmear;
    
  return;
}

void  getMetSyst(TVector2& met, TVector2& metUp, TVector2& metDown, TString mode, 
                 TVector2& singleVec, TVector2& singleVecUp, TVector2& singleVecDown,
                 MitGPTreeReduced& myTree) {
  //Main function used to propagate different uncertainties to met computation
  
  //case 1: consider just a single object with some variation
  if (mode == "single") {
    metUp = met + singleVec - singleVecUp;
    metDown = met + singleVec - singleVecDown;    
  }
  //case 2: consider the uncertainties on leptons
  else if (mode == "leptons") {
    float eleScaleUnc = 0.006;
    float muScaleUnc = 0.002;
    for (int iLep = 0; iLep < myTree.nlep_; iLep++) {
      TVector2 myLep;
      //Assume by default electron
      float thisScaleUnc = eleScaleUnc;
      if (iLep == 0) {
        //DEBUG phi set to zero
        myLep.SetMagPhi(myTree.lep1Pt_, 0);
        if (myTree.lep1Id_ > 12) thisScaleUnc = muScaleUnc;
      }
      else if (iLep == 1) {
        //DEBUG
        myLep.SetMagPhi(myTree.lep2Pt_, 0);
        if (myTree.lep2Id_ > 12) thisScaleUnc = muScaleUnc;
      }
      else if (iLep == 2) {
        //DEBUG
        myLep.SetMagPhi(myTree.lep3Pt_, 0);
        if (myTree.lep3Id_ > 12) thisScaleUnc = muScaleUnc;
      }
      else 
        continue;
      metUp = met + myLep - myLep*(1.+thisScaleUnc);
      metDown = met + myLep - myLep*(1.-thisScaleUnc);

      //Store the variated leading lepton vector
      if (iLep == 0) {
        singleVecUp.SetMagPhi(myLep.Mod()*(1.+thisScaleUnc), myLep.Phi());
        singleVecDown.SetMagPhi(myLep.Mod()*(1.-thisScaleUnc), myLep.Phi());
      }
    }    
  }
  //case 3: consider the uncertainties on jets
  else if (mode == "jets") {
    TVector2 myUncEnergy(met); //unclustered energy
    metUp = met;
    metDown = met;
    for (int iJet = 0; iJet < myTree.nalljets_; iJet++) {
      TVector2 myJet;
      float thisScaleUnc;
      if (iJet == 0) {
        myJet.SetMagPhi(myTree.jet1Pt_, myTree.jet1Phi_);
        thisScaleUnc = getJECError(myTree.jet1Eta_, myTree.jet1Pt_);
      }
      else if (iJet == 1) {
        myJet.SetMagPhi(myTree.jet2Pt_, myTree.jet2Phi_);
        thisScaleUnc = getJECError(myTree.jet2Eta_, myTree.jet2Pt_);
      }
      else if (iJet == 2) {
        myJet.SetMagPhi(myTree.jet3Pt_, myTree.jet3Phi_);
        thisScaleUnc = getJECError(myTree.jet3Eta_, myTree.jet3Pt_);
      }
      else if (iJet == 3) {
        myJet.SetMagPhi(myTree.jet4Pt_, myTree.jet4Phi_);
        thisScaleUnc = getJECError(myTree.jet4Eta_, myTree.jet4Pt_);
      }
      else 
        continue;
      myUncEnergy = myUncEnergy + myJet;
      metUp = metUp + myJet - myJet*(1.+thisScaleUnc);
      metDown = metDown + myJet - myJet*(1.-thisScaleUnc);

      //Store the variated leading jet vector
      if (iJet == 0) {
        singleVecUp.SetMagPhi(myJet.Mod()*(1.+thisScaleUnc), myJet.Phi());
        singleVecDown.SetMagPhi(myJet.Mod()*(1.-thisScaleUnc), myJet.Phi());
      }
    }
    //Finalize the uncertainty by including the effect on the unclustered energy
    TVector2 thePhoton;
    thePhoton.SetMagPhi(myTree.phoEt_,myTree.phoPhi_);
    myUncEnergy = -1.*(myUncEnergy + thePhoton);
    metUp = metUp + myUncEnergy - myUncEnergy*(1.+0.1);
    metDown = metDown + myUncEnergy - myUncEnergy*(1.-0.1);    
  }
  
  return;
}


