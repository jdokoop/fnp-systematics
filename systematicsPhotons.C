//----------------------------------------------------------
// Take fit to the direct photon spectra from the following
// publications
// * PPG060 [PRL 98 012002]
// * PPG086 [PRL 104 132301]
// * PPG 136 [PRD 86 072008]
//
// J. Orjuela-Koop
// May. 2017
//----------------------------------------------------------

//----------------------
// Variables
//----------------------

const int NPOINTSNLO136 = 19;
const int NPOINTSNLO060 = 27;
const int NPOINTS060 = 17;
const int NPOINTS086 = 3;//5;
const int NPOINTS136 = 18;

TGraph *g_nlo_136;
TGraph *g_nlo_136_npf;

TGraph *g_nlo_060;
TGraph *g_nlo_060_npf;

TGraphErrors *g_060_spectrum;
TGraphErrors *g_060_spectrum_npf;

TGraphErrors *g_086_spectrum;
TGraphErrors *g_086_spectrum_npf;

TGraphErrors *g_136_spectrum;
TGraphErrors *g_136_spectrum_npf;

//Combining all published values
TGraphErrors *g_combined;
TGraphErrors *g_combined_npf;

TF1 *f_published_060_spectrum_fit;
TF1 *f_published_060_spectrum_fit_extrapolated;
TF1 *f_published_060_spectrum_fit_extrapolated_npf;

TF1 *f_published_086_spectrum_fit;
TF1 *f_published_086_spectrum_fit_extrapolated;
TF1 *f_published_086_spectrum_fit_extrapolated_npf;

TF1 *f_published_136_spectrum_fit;
TF1 *f_published_136_spectrum_fit_extrapolated;
TF1 *f_published_136_spectrum_fit_extrapolated_npf;

TF1 *f_combined_fit;
TF1 *f_combined_fit_npf;

TBox *systematicErrors060[NPOINTS060];
TBox *systematicErrors_npf060[NPOINTS060];

TBox *systematicErrors086[NPOINTS086];
TBox *systematicErrors_npf086[NPOINTS086];

TBox *systematicErrors136[NPOINTS136];
TBox *systematicErrors_npf136[NPOINTS136];

TBox *systematicErrorsCombined[NPOINTS136 + NPOINTS086 + NPOINTS060];
TBox *systematicErrors_npfCombined[NPOINTS136 + NPOINTS086 + NPOINTS060];

///////////////////////////////////////////////////////////
// PPG060 - Dd^3\sigma/dp^3 = (1/2pi pT) dN/dpT
///////////////////////////////////////////////////////////

float data_x_ppg060[NPOINTS060] = {3.25,
                                   3.75,
                                   4.25,
                                   4.75,
                                   5.25,
                                   5.75,
                                   6.25,
                                   6.75,
                                   7.25,
                                   7.75,
                                   8.25,
                                   8.75,
                                   9.25,
                                   9.75,
                                   11.0,
                                   13.0,
                                   15.0
                                  };

float data_y_ppg060[NPOINTS060] = {2.22E-05,
                                   9.84E-06,
                                   4.38E-06,
                                   1.64E-06,
                                   9.82E-07,
                                   5.24E-07,
                                   3.03E-07,
                                   2.42E-07,
                                   8.90E-08,
                                   1.13E-07,
                                   6.68E-08,
                                   3.30E-08,
                                   2.57E-08,
                                   2.16E-08,
                                   1.17E-08,
                                   1.83E-09,
                                   2.04E-09
                                  };

float err_y_stat_ppg060[NPOINTS060] = {1.10E-06,
                                       5.41E-07,
                                       2.94E-07,
                                       1.74E-07,
                                       1.10E-07,
                                       7.15E-08,
                                       4.85E-08,
                                       3.51E-08,
                                       2.39E-08,
                                       2.03E-08,
                                       1.49E-08,
                                       1.10E-08,
                                       8.73E-09,
                                       6.76E-09,
                                       2.02E-09,
                                       9.05E-10,
                                       6.55E-10
                                      };

float err_y_syst_ppg060[NPOINTS060] = {2.54E-05,
                                       8.74E-06,
                                       3.19E-06,
                                       1.17E-06,
                                       5.34E-07,
                                       2.54E-07,
                                       1.27E-07,
                                       7.85E-08,
                                       3.62E-08,
                                       2.84E-08,
                                       1.71E-08,
                                       9.51E-09,
                                       6.47E-09,
                                       4.97E-09,
                                       2.31E-09,
                                       4.61E-10,
                                       3.42E-10
                                      };

///////////////////////////////////////////////////////////
// PPG086 - Dd^3\sigma/dp^3 = (1/2pi pT) dN/dpT
///////////////////////////////////////////////////////////

float data_x_ppg086[NPOINTS086] = {1.69126, 2.19773, 2.70364};//, 3.3488, 4.37494};
float data_y_ppg086[NPOINTS086] = {0.000609745, 0.000163268, 4.65505e-05};//, 7.78764e-06, 1.2448e-06};
float err_y_stat_ppg086[NPOINTS086] = {0.000156004, 4.78465e-05, 1.89656e-05};//, 5.06527e-06, 1.57247e-06};
float err_y_syst_ppg086[NPOINTS086] = {0.000339391, 6.91015e-05, 1.9658e-05};//, 3.99712e-06, 7.87164e-07};

///////////////////////////////////////////////////////////
// PPG136 - Dd^3\sigma/dp^3 = (1/2pi pT) dN/dpT
///////////////////////////////////////////////////////////

float data_x_ppg136[NPOINTS136] = {5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75, 11.00, 13.00, 15.00, 17.00, 19.00, 21.00, 23.00, 25.00};
float data_y_ppg136[NPOINTS136] = {1.14e+03, 6.13e+02, 3.48e+02, 2.31e+02, 1.36e+02, 9.29e+01, 6.70e+01, 4.83e+01, 3.21e+01, 2.04e+01, 9.81e+00, 2.97e+00, 1.06e+00, 3.38e-01, 1.73e-01, 8.82e-02, 4.22e-02, 2.87e-02};
float err_y_stat_ppg136[NPOINTS136] = {3.04e+01, 1.92e+01, 1.27e+01, 8.50e+00, 6.12e+00, 4.41e+00, 3.22e+00, 2.45e+00, 1.89e+00, 1.46e+00, 4.23e-01, 1.89e-01, 9.85e-02, 5.51e-02, 3.37e-02, 2.06e-02, 1.52e-02, 1.41e-02};
float err_y_syst_ppg136[NPOINTS136] = {4.78e+02, 2.21e+02, 1.01e+02, 6.24e+01, 3.13e+01, 1.95e+01, 1.34e+01, 9.18e+00, 6.10e+00, 3.68e+00, 1.67e+00, 4.75e-01, 1.69e-01, 5.42e-02, 2.77e-02, 1.50e-02, 7.18e-03, 4.30e-03};

///////////////////////////////////////////////////////////
// COMBINED SPECTRA
///////////////////////////////////////////////////////////

float data_x[NPOINTS060 + NPOINTS086] = {0.0};
float data_y[NPOINTS060 + NPOINTS086] = {0.0};
float err_y_syst[NPOINTS060 + NPOINTS086] = {0.0};
float err_y_stat[NPOINTS060 + NPOINTS086] = {0.0};

float data_y_npf[NPOINTS060 + NPOINTS086] = {0.0};
float err_y_syst_npf[NPOINTS060 + NPOINTS086] = {0.0};
float err_y_stat_npf[NPOINTS060 + NPOINTS086] = {0.0};

///////////////////////////////////////////////////////////
// NLO
///////////////////////////////////////////////////////////

float data_x_nlo_136[NPOINTSNLO136] = {1.0186335403726707, 1.204968944099379, 1.391304347826087, 1.5900621118012421, 1.8136645962732918, 2.062111801242236, 2.31055900621118,  2.62111801242236,  2.9192546583850927, 3.2546583850931676, 3.6894409937888195, 4.049689440993789, 4.447204968944099, 4.84472049689441,  5.24223602484472,  5.577639751552795, 6.049689440993788, 6.534161490683229, 6.981366459627328};
float data_y_nlo_136[NPOINTSNLO136] = {0.003815385822083397, 0.001701254279852592, 0.000861756594301430, 0.000436515832240166, 0.000221113563918577, 0.000116867110851035, 0.000064451164598286, 0.000032647216092701, 0.000018786511982616, 0.000010360593046256, 0.000005248074602497, 0.000003019951720402, 0.000001813266084165, 0.000001088737237013, 7.117175782354078e-7, 5.065419111692708e-7, 3.173499699587764e-7, 1.988206725093311e-7, 1.356149887611501e-7};

float data_x_nlo_060[NPOINTSNLO060] = {2.33587786259542,
                                       2.549618320610687,
                                       2.793893129770992,
                                       3.099236641221374,
                                       3.435114503816794,
                                       3.770992366412213,
                                       4.076335877862595,
                                       4.595419847328245,
                                       4.961832061068701,
                                       5.358778625954198,
                                       5.786259541984732,
                                       6.305343511450381,
                                       6.793893129770993,
                                       7.221374045801525,
                                       7.740458015267176,
                                       8.259541984732824,
                                       8.80916030534351,
                                       9.572519083969466,
                                       10.09160305343511,
                                       10.73282442748091,
                                       11.34351145038167,
                                       12.04580152671755,
                                       12.68702290076335,
                                       13.35877862595419,
                                       14.27480916030534,
                                       14.97709923664122,
                                       15.52671755725190
                                      };

float data_y_nlo_060[NPOINTSNLO060] = {72997.2063589229,
                                       42843.8015982454,
                                       25140.6813539916,
                                       13944.0572906539,
                                       7519.04715836717,
                                       4534.39609789014,
                                       2812.64269954585,
                                       1472.94378746564,
                                       913.262978030331,
                                       582.181720713649,
                                       381.570078673493,
                                       229.813037996709,
                                       150.558395845290,
                                       101.476732339781,
                                       66.4666933043638,
                                       46.0397726713958,
                                       31.0044336918035,
                                       18.1274485595098,
                                       12.9125354739595,
                                       8.69008477856125,
                                       6.01554477495446,
                                       4.04670939914393,
                                       2.88008961896363,
                                       2.04935534991032,
                                       1.30167875532326,
                                       0.92602427003898,
                                       0.71720242076655
                                      };

float err_x[NPOINTS060] = {0};

//Variations on points
TGraphErrors *g_060_spectrum_var1;
TGraphErrors *g_060_spectrum_var1_npf;
TF1 *f_spectrum_fit_var1;
TF1 *f_spectrum_fit_var1_extrapolated;
TF1 *f_spectrum_fit_var1_extrapolated_npf;

TGraphErrors *g_060_spectrum_var2;
TGraphErrors *g_060_spectrum_var2_npf;
TF1 *f_spectrum_fit_var2;
TF1 *f_spectrum_fit_var2_extrapolated;
TF1 *f_spectrum_fit_var2_extrapolated_npf;

//Variations on low-pT behavior
TF1 *f_spectrum_fit_var3;
TF1 *f_spectrum_fit_var3_extrapolated;
TF1 *f_spectrum_fit_var3_extrapolated_npf;

TF1 *f_spectrum_fit_var4;
TF1 *f_spectrum_fit_var4_extrapolated;
TF1 *f_spectrum_fit_var4_extrapolated_npf;

//----------------------
// Functions
//----------------------

void makeNLO()
{
  g_nlo_136 = new TGraph(NPOINTSNLO136, data_x_nlo_136, data_y_nlo_136);

  //Multiply NLO calculation by the phase space factor 2pi x pT
  for (int i = 0; i < NPOINTSNLO136; i++)
  {
    data_y_nlo_136[i] = 2 * TMath::Pi() * data_x_nlo_136[i] * data_y_nlo_136[i];
  }

  g_nlo_136_npf = new TGraph(NPOINTSNLO136, data_x_nlo_136, data_y_nlo_136);

  g_nlo_060 = new TGraph(NPOINTSNLO060, data_x_nlo_060, data_y_nlo_060);

  //Multiply NLO calculation by the phase space factor 2pi x pT
  //Also multiply by 1E-9 to make the units mb instead of pb
  for (int i = 0; i < NPOINTSNLO060; i++)
  {
    data_y_nlo_060[i] = 1E-9 * 2 * TMath::Pi() * data_x_nlo_060[i] * data_y_nlo_060[i];
  }

  g_nlo_060_npf = new TGraph(NPOINTSNLO060, data_x_nlo_060, data_y_nlo_060);
}


void makePublishedSpectrum060()
{
  //Make spectrum with no phase space factor
  float data_y_ppg060_npf[NPOINTS060];
  float err_y_stat_ppg060_npf[NPOINTS060];
  for (int i = 0; i < NPOINTS060; i++)
  {
    data_y_ppg060_npf[i] = 2 * TMath::Pi() * data_x_ppg060[i] * data_y_ppg060[i];
    err_y_stat_ppg060_npf[i] = 2 * TMath::Pi() * data_x_ppg060[i] * err_y_stat_ppg060[i];
  }

  for (int i = 0; i < NPOINTS060; i++)
  {
    systematicErrors_npf060[i] = new TBox(data_x_ppg060[i] - 0.15, data_y_ppg060_npf[i] - (2 * TMath::Pi() * data_x_ppg060[i]*err_y_syst_ppg060[i] / 2.0), data_x_ppg060[i] + 0.15, data_y_ppg060_npf[i] + (2 * TMath::Pi() * data_x_ppg060[i]*err_y_syst_ppg060[i] / 2.0));
    systematicErrors_npf060[i]->SetLineColor(kBlack);
    systematicErrors_npf060[i]->SetFillStyle(0);
  }

  //Make normal spectrum
  for (int i = 0; i < NPOINTS060; i++)
  {
    systematicErrors060[i] = new TBox(data_x_ppg060[i] - 0.15, data_y_ppg060[i] - (err_y_syst_ppg060[i] / 2.0), data_x_ppg060[i] + 0.15, data_y_ppg060[i] + (err_y_syst_ppg060[i] / 2.0));
    systematicErrors060[i]->SetLineColor(kBlack);
    systematicErrors060[i]->SetFillStyle(0);
  }

  g_060_spectrum = new TGraphErrors(NPOINTS060, data_x_ppg060, data_y_ppg060, err_x, err_y_stat_ppg060);
  g_060_spectrum_npf = new TGraphErrors(NPOINTS060, data_x_ppg060, data_y_ppg060_npf, err_x, err_y_stat_ppg060_npf);

  f_published_060_spectrum_fit = new TF1("f_published_060_spectrum_fit", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 3.25, 15.0);
  f_published_060_spectrum_fit->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);
  f_published_060_spectrum_fit->SetLineColor(kBlack);

  f_published_060_spectrum_fit_extrapolated = new TF1("f_published_060_spectrum_fit_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_published_060_spectrum_fit_extrapolated->SetParameters(f_published_060_spectrum_fit->GetParameter(0), f_published_060_spectrum_fit->GetParameter(1), f_published_060_spectrum_fit->GetParameter(2), f_published_060_spectrum_fit->GetParameter(3), f_published_060_spectrum_fit->GetParameter(4));
  f_published_060_spectrum_fit_extrapolated->SetLineColor(kBlack);

  f_published_060_spectrum_fit_extrapolated_npf = new TF1("f_published_060_spectrum_fit_extrapolated_npf", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_published_060_spectrum_fit_extrapolated_npf->SetParameters(f_published_060_spectrum_fit->GetParameter(0), f_published_060_spectrum_fit->GetParameter(1), f_published_060_spectrum_fit->GetParameter(2), f_published_060_spectrum_fit->GetParameter(3), f_published_060_spectrum_fit->GetParameter(4));
  f_published_060_spectrum_fit_extrapolated_npf->SetLineColor(kBlack);
}


void makePublishedSpectrum086()
{
  //Make spectrum with no phase space factor
  float data_y_ppg086_npf[NPOINTS086];
  float err_y_stat_ppg086_npf[NPOINTS086];
  for (int i = 0; i < NPOINTS086; i++)
  {
    data_y_ppg086_npf[i] = 2 * TMath::Pi() * data_x_ppg086[i] * data_y_ppg086[i];
    err_y_stat_ppg086_npf[i] = 2 * TMath::Pi() * data_x_ppg086[i] * err_y_stat_ppg086[i];
  }

  for (int i = 0; i < NPOINTS086; i++)
  {
    systematicErrors_npf086[i] = new TBox(data_x_ppg086[i] - 0.15, data_y_ppg086_npf[i] - (2 * TMath::Pi() * data_x_ppg086[i]*err_y_syst_ppg086[i] / 2.0), data_x_ppg086[i] + 0.15, data_y_ppg086_npf[i] + (2 * TMath::Pi() * data_x_ppg086[i]*err_y_syst_ppg086[i] / 2.0));
    systematicErrors_npf086[i]->SetLineColor(kBlack);
    systematicErrors_npf086[i]->SetFillStyle(0);
  }

  //Make normal spectrum
  for (int i = 0; i < NPOINTS086; i++)
  {
    systematicErrors086[i] = new TBox(data_x_ppg086[i] - 0.15, data_y_ppg086[i] - (err_y_syst_ppg086[i] / 2.0), data_x_ppg086[i] + 0.15, data_y_ppg086[i] + (err_y_syst_ppg086[i] / 2.0));
    systematicErrors086[i]->SetLineColor(kBlack);
    systematicErrors086[i]->SetFillStyle(0);
  }

  g_086_spectrum = new TGraphErrors(NPOINTS086, data_x_ppg086, data_y_ppg086, err_x, err_y_stat_ppg086);
  g_086_spectrum_npf = new TGraphErrors(NPOINTS086, data_x_ppg086, data_y_ppg086_npf, err_x, err_y_stat_ppg086_npf);

  f_published_086_spectrum_fit = new TF1("f_published_086_spectrum_fit", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 3.25, 15.0);
  f_published_086_spectrum_fit->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);
  f_published_086_spectrum_fit->SetLineColor(kBlack);

  f_published_086_spectrum_fit_extrapolated = new TF1("f_published_086_spectrum_fit_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_published_086_spectrum_fit_extrapolated->SetParameters(f_published_086_spectrum_fit->GetParameter(0), f_published_086_spectrum_fit->GetParameter(1), f_published_086_spectrum_fit->GetParameter(2), f_published_086_spectrum_fit->GetParameter(3), f_published_086_spectrum_fit->GetParameter(4));
  f_published_086_spectrum_fit_extrapolated->SetLineColor(kBlack);

  f_published_086_spectrum_fit_extrapolated_npf = new TF1("f_published_086_spectrum_fit_extrapolated_npf", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_published_086_spectrum_fit_extrapolated_npf->SetParameters(f_published_086_spectrum_fit->GetParameter(0), f_published_086_spectrum_fit->GetParameter(1), f_published_086_spectrum_fit->GetParameter(2), f_published_086_spectrum_fit->GetParameter(3), f_published_086_spectrum_fit->GetParameter(4));
  f_published_086_spectrum_fit_extrapolated_npf->SetLineColor(kBlack);
}


void makePublishedSpectrum136()
{
  //Multiply each point by 1E-9 to convert the units from pb to mb, to match the other measurements
  for (int i = 0; i < NPOINTS136; i++)
  {
    data_y_ppg136[i] = 1E-9 * data_y_ppg136[i];
    err_y_syst_ppg136[i] = 1E-9 * err_y_syst_ppg136[i];
    err_y_stat_ppg136[i] = 1E-9 * err_y_stat_ppg136[i];
  }

  //Make spectrum with no phase space factor
  float data_y_ppg136_npf[NPOINTS136];
  float err_y_stat_ppg136_npf[NPOINTS136];
  for (int i = 0; i < NPOINTS136; i++)
  {
    data_y_ppg136_npf[i] = 2 * TMath::Pi() * data_x_ppg136[i] * data_y_ppg136[i];
    err_y_stat_ppg136_npf[i] = 2 * TMath::Pi() * data_x_ppg136[i] * err_y_stat_ppg136[i];
  }

  for (int i = 0; i < NPOINTS136; i++)
  {
    systematicErrors_npf136[i] = new TBox(data_x_ppg136[i] - 0.15, data_y_ppg136_npf[i] - (2 * TMath::Pi() * data_x_ppg136[i]*err_y_syst_ppg136[i] / 2.0), data_x_ppg136[i] + 0.15, data_y_ppg136_npf[i] + (2 * TMath::Pi() * data_x_ppg136[i]*err_y_syst_ppg136[i] / 2.0));
    systematicErrors_npf136[i]->SetLineColor(kBlack);
    systematicErrors_npf136[i]->SetFillStyle(0);
  }

  //Make normal spectrum
  for (int i = 0; i < NPOINTS136; i++)
  {
    systematicErrors136[i] = new TBox(data_x_ppg136[i] - 0.15, data_y_ppg136[i] - (err_y_syst_ppg136[i] / 2.0), data_x_ppg136[i] + 0.15, data_y_ppg136[i] + (err_y_syst_ppg136[i] / 2.0));
    systematicErrors136[i]->SetLineColor(kBlack);
    systematicErrors136[i]->SetFillStyle(0);
  }

  g_136_spectrum = new TGraphErrors(NPOINTS136, data_x_ppg136, data_y_ppg136, err_x, err_y_stat_ppg136);
  g_136_spectrum_npf = new TGraphErrors(NPOINTS136, data_x_ppg136, data_y_ppg136_npf, err_x, err_y_stat_ppg136_npf);

  f_published_136_spectrum_fit = new TF1("f_published_136_spectrum_fit", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 3.25, 15.0);
  f_published_136_spectrum_fit->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);
  f_published_136_spectrum_fit->SetLineColor(kBlack);

  f_published_136_spectrum_fit_extrapolated = new TF1("f_published_136_spectrum_fit_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_published_136_spectrum_fit_extrapolated->SetParameters(f_published_136_spectrum_fit->GetParameter(0), f_published_136_spectrum_fit->GetParameter(1), f_published_136_spectrum_fit->GetParameter(2), f_published_136_spectrum_fit->GetParameter(3), f_published_136_spectrum_fit->GetParameter(4));
  f_published_136_spectrum_fit_extrapolated->SetLineColor(kBlack);

  f_published_136_spectrum_fit_extrapolated_npf = new TF1("f_published_136_spectrum_fit_extrapolated_npf", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_published_136_spectrum_fit_extrapolated_npf->SetParameters(f_published_136_spectrum_fit->GetParameter(0), f_published_136_spectrum_fit->GetParameter(1), f_published_136_spectrum_fit->GetParameter(2), f_published_136_spectrum_fit->GetParameter(3), f_published_136_spectrum_fit->GetParameter(4));
  f_published_136_spectrum_fit_extrapolated_npf->SetLineColor(kBlack);
}


void combinePublications()
{
  for (int i = 0; i < NPOINTS060; i++)
  {
    data_x[i] = data_x_ppg060[i];
    data_y[i] = data_y_ppg060[i];
    err_y_stat[i] = err_y_stat_ppg060[i];
    err_y_syst[i] = err_y_syst_ppg060[i];

    data_y_npf[i] = data_y_ppg060[i] * 2 * TMath::Pi() * data_x_ppg060[i];
    err_y_stat_npf[i] = err_y_stat_ppg060[i] * 2 * TMath::Pi() * data_x_ppg060[i];
    err_y_syst_npf[i] = err_y_syst_ppg060[i] * 2 * TMath::Pi() * data_x_ppg060[i];
  }

  for (int i = 0; i < NPOINTS086; i++)
  {
    data_x[i + NPOINTS060] = data_x_ppg086[i];
    data_y[i + NPOINTS060] = data_y_ppg086[i];
    err_y_stat[i + NPOINTS060] = err_y_stat_ppg086[i];
    err_y_syst[i + NPOINTS060] = err_y_syst_ppg086[i];

    data_y_npf[i + NPOINTS060] = data_y_ppg086[i] * 2 * TMath::Pi() * data_x_ppg086[i];
    err_y_stat_npf[i + NPOINTS060] = err_y_stat_ppg086[i] * 2 * TMath::Pi() * data_x_ppg086[i];
    err_y_syst_npf[i + NPOINTS060] = err_y_syst_ppg086[i] * 2 * TMath::Pi() * data_x_ppg086[i];
  }

  /*
    for (int i = 0; i < NPOINTS136; i++)
    {
      data_x[i + NPOINTS060 + NPOINTS086] = data_x_ppg136[i];
      data_y[i + NPOINTS060 + NPOINTS086] = data_y_ppg136[i];
      err_y_stat[i + NPOINTS060 + NPOINTS086] = err_y_stat_ppg136[i];
      err_y_syst[i + NPOINTS060 + NPOINTS086] = err_y_syst_ppg136[i];

      data_y_npf[i + NPOINTS060 + NPOINTS086] = data_y_ppg136[i] * 2 * TMath::Pi() * data_x_ppg136[i];
      err_y_stat_npf[i + NPOINTS060 + NPOINTS086] = err_y_stat_ppg136[i] * 2 * TMath::Pi() * data_x_ppg136[i];
      err_y_syst_npf[i + NPOINTS060 + NPOINTS086] = err_y_syst_ppg136[i] * 2 * TMath::Pi() * data_x_ppg136[i];
    }
    */

  g_combined = new TGraphErrors(NPOINTS060 + NPOINTS086, data_x, data_y, err_x, err_y_stat);
  g_combined_npf = new TGraphErrors(NPOINTS060 + NPOINTS086, data_x, data_y_npf, err_x, err_y_stat_npf);

  //Now, make systematic boxes
  for (int i = 0; i < NPOINTS060 + NPOINTS086; i++)
  {
    systematicErrorsCombined[i] = new TBox(data_x[i] - 0.15, data_y[i] - (err_y_syst[i] / 2.0), data_x[i] + 0.15, data_y[i] + (err_y_syst[i] / 2.0));
    systematicErrorsCombined[i]->SetLineColor(kBlack);
    systematicErrorsCombined[i]->SetFillStyle(0);

    systematicErrors_npfCombined[i] = new TBox(data_x[i] - 0.15, data_y_npf[i] - (err_y_syst_npf[i] / 2.0), data_x[i] + 0.15, data_y_npf[i] + (err_y_syst_npf[i] / 2.0));
    systematicErrors_npfCombined[i]->SetLineColor(kBlack);
    systematicErrors_npfCombined[i]->SetFillStyle(0);
  }

  //Fit the resulting combined spectra with a modified Hagedorn function
  f_combined_fit = new TF1("f_combined_fit", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0, 18.0);
  f_combined_fit->SetParameters(2.42345e-01, -8.27585e-02, 9.18447e-03, 4.13943e+00, 1.36974e+01);
  //f_combined_fit->SetParameters(0.00488828, 0.8888921, 0.209916, 1.43194, 6.48914); //Looks most like PPG162
  g_combined->Fit(f_combined_fit, "0R");

  f_combined_fit_npf = new TF1("f_combined_fit_npf", "2*TMath::Pi()*x*[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0, 18.0);
  f_combined_fit_npf->SetParameters(f_combined_fit->GetParameter(0), f_combined_fit->GetParameter(1), f_combined_fit->GetParameter(2), f_combined_fit->GetParameter(3), f_combined_fit->GetParameter(4));

  //Result of above fit
  //f_combined_fit_npf->SetParameters(2.42345e-01, -8.27585e-02, 9.18447e-03, 4.13943e+00, 1.36974e+01);
}


void defineVariation1()
{
  //Define the seventh point in the spectrum as a tipping point, and tilt clockwise about that point by an amount proportional to the point's pT
  float data_y_var1[NPOINTS060];
  int tippingPointIndex = 6;
  float pTextreme = data_x_ppg060[0];
  float pT0 = data_x[tippingPointIndex];

  for (int i = 0; i < NPOINTS060; i++)
  {
    float scaling = (err_y_syst[i] / 2.0) * ((data_x[i] - pT0) / (pTextreme - pT0));
    data_y_var1[i] = data_y[i] + scaling;
  }

  g_060_spectrum_var1 = new TGraphErrors(NPOINTS060, data_x, data_y_var1, err_x, err_y_stat);
  g_060_spectrum_var1->SetMarkerColor(kRed);
  g_060_spectrum_var1->SetLineColor(kRed);
  g_060_spectrum_var1->SetMarkerStyle(20);
  g_060_spectrum_var1->SetMarkerSize(0.7);

  //Fit the data with a modified Hagedorn function
  f_spectrum_fit_var1 = new TF1("f_spectrum_fit_var1", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 3.25, 15.0);
  f_spectrum_fit_var1->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);
  g_060_spectrum_var1->Fit(f_spectrum_fit_var1, "Q0R");

  f_spectrum_fit_var1_extrapolated = new TF1("f_spectrum_fit_var1_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_spectrum_fit_var1_extrapolated->SetParameters(f_spectrum_fit_var1->GetParameter(0), f_spectrum_fit_var1->GetParameter(1), f_spectrum_fit_var1->GetParameter(2), f_spectrum_fit_var1->GetParameter(3), f_spectrum_fit_var1->GetParameter(4));

  f_spectrum_fit_var1_extrapolated_npf = new TF1("f_spectrum_fit_var1_extrapolated_npf", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
  f_spectrum_fit_var1_extrapolated_npf->SetParameters(f_spectrum_fit_var1->GetParameter(0), f_spectrum_fit_var1->GetParameter(1), f_spectrum_fit_var1->GetParameter(2), f_spectrum_fit_var1->GetParameter(3), f_spectrum_fit_var1->GetParameter(4));

  //Spectrum without phase space factor
  float data_y_npf[NPOINTS060];

  for (int i = 0; i < NPOINTS060; i++)
  {
    data_y_npf[i] = 2 * TMath::Pi() * data_x[i] * data_y_var1[i];
  }

  g_060_spectrum_var1_npf = new TGraphErrors(NPOINTS060, data_x, data_y_npf, err_x, err_y_stat);
  g_060_spectrum_var1_npf->SetMarkerColor(kRed);
  g_060_spectrum_var1_npf->SetLineColor(kRed);
  g_060_spectrum_var1_npf->SetMarkerStyle(20);
  g_060_spectrum_var1_npf->SetMarkerSize(0.7);
}


void defineVariation2()
{
  //Define the seventh point in the spectrum as a tipping point, and tilt clockwise about that point by an amount proportional to the point's pT
  float data_y_ppg060_var2[NPOINTS060];
  int tippingPointIndex = 6;
  float pTextreme = data_x_ppg060[0];
  float pT0 = data_x_ppg060[tippingPointIndex];

  for (int i = 0; i < NPOINTS060; i++)
  {
    float scaling = -1 * (err_y_syst_ppg060[i] / 2.0) * ((data_x_ppg060[i] - pT0) / (pTextreme - pT0));
    data_y_ppg060_var2[i] = data_y_ppg060[i] + scaling;
  }

  g_060_spectrum_var2 = new TGraphErrors(NPOINTS060, data_x_ppg060, data_y_ppg060_var2, err_x, err_y_stat_ppg060);
  g_060_spectrum_var2->SetMarkerColor(kBlue);
  g_060_spectrum_var2->SetLineColor(kBlue);
  g_060_spectrum_var2->SetMarkerStyle(20);
  g_060_spectrum_var2->SetMarkerSize(0.7);

  //Fit the data with a modified Hagedorn function
  f_spectrum_fit_var2 = new TF1("f_spectrum_fit_var2", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0, 15.0);
  //f_spectrum_fit_var2->SetParameters(26.1, 2.0, 2.5, 0.28, 5.89);
  f_spectrum_fit_var2->SetParameters(5.03, -0.45, 0.014, 0.27, 4.72);
  g_060_spectrum_var2->Fit(f_spectrum_fit_var2, "Q0R");
  f_spectrum_fit_var2->SetLineColor(kBlue);

  f_spectrum_fit_var2_extrapolated = new TF1("f_spectrum_fit_var2_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_spectrum_fit_var2_extrapolated->SetParameters(f_spectrum_fit_var2->GetParameter(0), f_spectrum_fit_var2->GetParameter(1), f_spectrum_fit_var2->GetParameter(2), f_spectrum_fit_var2->GetParameter(3), f_spectrum_fit_var2->GetParameter(4));
  f_spectrum_fit_var2_extrapolated->SetLineColor(kBlue);

  f_spectrum_fit_var2_extrapolated_npf = new TF1("f_spectrum_fit_var2_extrapolated_npf", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
  f_spectrum_fit_var2_extrapolated_npf->SetParameters(f_spectrum_fit_var2->GetParameter(0), f_spectrum_fit_var2->GetParameter(1), f_spectrum_fit_var2->GetParameter(2), f_spectrum_fit_var2->GetParameter(3), f_spectrum_fit_var2->GetParameter(4));
  f_spectrum_fit_var2_extrapolated_npf->SetLineColor(kBlue);

  //Spectrum without phase space factor
  float data_y_ppg060_npf[NPOINTS060];

  for (int i = 0; i < NPOINTS060; i++)
  {
    data_y_ppg060_npf[i] = 2 * TMath::Pi() * data_x_ppg060[i] * data_y_ppg060_var2[i];
  }

  g_060_spectrum_var2_npf = new TGraphErrors(NPOINTS060, data_x_ppg060, data_y_ppg060_npf, err_x, err_y_stat_ppg060);
  g_060_spectrum_var2_npf->SetMarkerColor(kBlue);
  g_060_spectrum_var2_npf->SetLineColor(kBlue);
  g_060_spectrum_var2_npf->SetMarkerStyle(20);
  g_060_spectrum_var2_npf->SetMarkerSize(0.7);
}


void defineVariation3()
{
  //Fit published spectrum with a Tsallis functional form
  f_spectrum_fit_var3 = new TF1("f_spectrum_fit_var3", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 10.0*([1] - 1))*([1]*[2] + 10.0)) * pow(([1]*[2] + TMath::Sqrt(10.0*10.0 + x*x))/([1]*[2]+10.0),-1*[1])", 0, 18.0);
  f_spectrum_fit_var3->SetParameter(0, 1.3);
  f_spectrum_fit_var3->SetParameter(1, 5.5);
  f_spectrum_fit_var3->SetParameter(2, 0.001);
  g_combined->Fit(f_spectrum_fit_var3, "Q0R");
  f_spectrum_fit_var3->SetLineColor(kOrange - 3);

  f_spectrum_fit_var3_extrapolated = new TF1("f_spectrum_fit_var3_extrapolated", "([0]*(([1]-1)*([1]-1))/(([1]*[2] + 10.0*([1] - 1))*([1]*[2] + 10.0)) * pow(([1]*[2] + TMath::Sqrt(10.0*10.0 + x*x))/([1]*[2]+10.0),-1*[1]))", 0.0, 18.0);
  f_spectrum_fit_var3_extrapolated->SetParameters(f_spectrum_fit_var3->GetParameter(0), f_spectrum_fit_var3->GetParameter(1), f_spectrum_fit_var3->GetParameter(2));
  f_spectrum_fit_var3_extrapolated->SetLineColor(kOrange - 3);

  f_spectrum_fit_var3_extrapolated_npf = new TF1("f_spectrum_fit_var3_extrapolated_npf", "2*TMath::Pi()*x*(([0]*(([1]-1)*([1]-1))/(([1]*[2] + 10.0*([1] - 1))*([1]*[2] + 10.0)) * pow(([1]*[2] + TMath::Sqrt(10.0*10.0 + x*x))/([1]*[2]+10.0),-1*[1])))", 0.0, 18.0);
  f_spectrum_fit_var3_extrapolated_npf->SetParameters(f_spectrum_fit_var3->GetParameter(0), f_spectrum_fit_var3->GetParameter(1), f_spectrum_fit_var3->GetParameter(2));
  f_spectrum_fit_var3_extrapolated_npf->SetLineColor(kOrange - 3);
}


void defineVariation4()
{
  //Fit published spectrum with a Tsallis functional form
  f_spectrum_fit_var4 = new TF1("f_spectrum_fit_var4", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1])", 0, 18.0);
  f_spectrum_fit_var4->SetParameter(0, 1.3);
  f_spectrum_fit_var4->SetParameter(1, 5.5);
  f_spectrum_fit_var4->SetParameter(2, 0.001);
  g_060_spectrum->Fit(f_spectrum_fit_var4, "Q0R");
  f_spectrum_fit_var4->SetLineColor(kSpring - 6);

  f_spectrum_fit_var4_extrapolated = new TF1("f_spectrum_fit_var4_extrapolated", "([0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1]))", 0.0, 18.0);
  f_spectrum_fit_var4_extrapolated->SetParameters(f_spectrum_fit_var4->GetParameter(0), f_spectrum_fit_var4->GetParameter(1), f_spectrum_fit_var4->GetParameter(2));
  f_spectrum_fit_var4_extrapolated->SetLineColor(kSpring - 6);

  f_spectrum_fit_var4_extrapolated_npf = new TF1("f_spectrum_fit_var4_extrapolated_npf", "2*TMath::Pi()*x*(([0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1])))", 0.0, 18.0);
  f_spectrum_fit_var4_extrapolated_npf->SetParameters(f_spectrum_fit_var4->GetParameter(0), f_spectrum_fit_var4->GetParameter(1), f_spectrum_fit_var4->GetParameter(2));
  f_spectrum_fit_var4_extrapolated_npf->SetLineColor(kSpring - 6);
}


void plotDataPublishedFit()
{
  TCanvas *c = new TCanvas("c", "PPG036 Direct Photon Spectrum + Hagedorn Fit", 700, 700);
  c->SetLogy();

  TH1F *hTemplate = new TH1F("hTemplate", "hTemplate", 100, 0, 20);
  hTemplate->SetTitle("");
  hTemplate->GetXaxis()->SetTitleFont(62);
  hTemplate->GetXaxis()->SetLabelFont(62);
  hTemplate->GetXaxis()->SetRangeUser(0, 20);
  hTemplate->GetYaxis()->SetTitleFont(62);
  hTemplate->GetYaxis()->SetLabelFont(62);
  hTemplate->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hTemplate->GetYaxis()->SetTitle("(1/2#pi p_{T}) dN/dp_{T}");
  hTemplate->GetYaxis()->SetRangeUser(1e-10, 10);
  hTemplate->Draw();


  g_combined->SetTitle("");
  g_combined->GetXaxis()->SetTitleFont(62);
  g_combined->GetXaxis()->SetLabelFont(62);
  g_combined->GetXaxis()->SetRangeUser(0, 18);
  g_combined->GetYaxis()->SetTitleFont(62);
  g_combined->GetYaxis()->SetLabelFont(62);
  g_combined->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  g_combined->GetYaxis()->SetTitle("(1/2#pi p_{T}) dN/dp_{T}");
  g_combined->GetYaxis()->SetRangeUser(1e-10, 5e-5);
  g_combined->SetMarkerStyle(20);
  g_combined->SetMarkerSize(0.8);
  g_combined->SetMarkerColor(kBlack);
  g_combined->Draw("P,same");
  //g_060_spectrum_var1->Draw("P,same");
  //g_060_spectrum_var2->Draw("P,same");
  f_published_060_spectrum_fit_extrapolated->Draw("same");
  f_combined_fit->Draw("same");
  //f_spectrum_fit_var1_extrapolated->Draw("same");
  //f_spectrum_fit_var2->Draw("same");

  //f_spectrum_fit_var3_extrapolated->Draw("same");
  //f_spectrum_fit_var4->Draw("same");

  for (int i = 0; i < NPOINTS060; i++)
  {
    systematicErrorsCombined[i]->Draw("same");
  }
}


void plotDataNoPhaseFactor()
{
  TCanvas *cNP = new TCanvas("cNP", "No Phase Space Factor", 700, 700);
  cNP->SetLogy();

  TH1F *hTemplate2 = new TH1F("hTemplate2", "hTemplate2", 100, 0, 20);
  hTemplate2->SetTitle("");
  hTemplate2->GetXaxis()->SetTitleFont(62);
  hTemplate2->GetXaxis()->SetLabelFont(62);
  hTemplate2->GetXaxis()->SetRangeUser(0, 20);
  hTemplate2->GetYaxis()->SetTitleFont(62);
  hTemplate2->GetYaxis()->SetLabelFont(62);
  hTemplate2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hTemplate2->GetYaxis()->SetTitle("dN/dp_{T}");
  hTemplate2->GetYaxis()->SetTitleOffset(1.3);
  hTemplate2->GetYaxis()->SetRangeUser(1e-8, 150);
  hTemplate2->Draw();

  g_combined_npf->SetTitle("");
  g_combined_npf->GetXaxis()->SetTitleFont(62);
  g_combined_npf->GetXaxis()->SetLabelFont(62);
  g_combined_npf->GetXaxis()->SetRangeUser(0, 18);
  g_combined_npf->GetYaxis()->SetTitleFont(62);
  g_combined_npf->GetYaxis()->SetLabelFont(62);
  g_combined_npf->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  g_combined_npf->GetYaxis()->SetTitle("dN/dp_{T}");
  g_combined_npf->GetYaxis()->SetRangeUser(1e-8, 150);
  g_combined_npf->SetMarkerStyle(20);
  g_combined_npf->SetMarkerSize(0.8);
  g_combined_npf->SetMarkerColor(kBlack);

  g_060_spectrum_npf->SetTitle("");
  g_060_spectrum_npf->GetXaxis()->SetTitleFont(62);
  g_060_spectrum_npf->GetXaxis()->SetLabelFont(62);
  g_060_spectrum_npf->GetXaxis()->SetRangeUser(0, 18);
  g_060_spectrum_npf->GetYaxis()->SetTitleFont(62);
  g_060_spectrum_npf->GetYaxis()->SetLabelFont(62);
  g_060_spectrum_npf->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  g_060_spectrum_npf->GetYaxis()->SetTitle("dN/dp_{T}");
  g_060_spectrum_npf->GetYaxis()->SetRangeUser(1e-8, 150);
  g_060_spectrum_npf->SetMarkerStyle(20);
  g_060_spectrum_npf->SetMarkerSize(0.8);
  g_060_spectrum_npf->SetMarkerColor(kBlack);

  g_086_spectrum_npf->SetTitle("");
  g_086_spectrum_npf->GetXaxis()->SetTitleFont(62);
  g_086_spectrum_npf->GetXaxis()->SetLabelFont(62);
  g_086_spectrum_npf->GetXaxis()->SetRangeUser(0, 18);
  g_086_spectrum_npf->GetYaxis()->SetTitleFont(62);
  g_086_spectrum_npf->GetYaxis()->SetLabelFont(62);
  g_086_spectrum_npf->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  g_086_spectrum_npf->GetYaxis()->SetTitle("dN/dp_{T}");
  g_086_spectrum_npf->GetYaxis()->SetRangeUser(1e-8, 150);
  g_086_spectrum_npf->SetMarkerStyle(20);
  g_086_spectrum_npf->SetMarkerSize(0.8);
  g_086_spectrum_npf->SetMarkerColor(kBlack);

  g_136_spectrum_npf->SetTitle("");
  g_136_spectrum_npf->GetXaxis()->SetTitleFont(62);
  g_136_spectrum_npf->GetXaxis()->SetLabelFont(62);
  g_136_spectrum_npf->GetXaxis()->SetRangeUser(0, 18);
  g_136_spectrum_npf->GetYaxis()->SetTitleFont(62);
  g_136_spectrum_npf->GetYaxis()->SetLabelFont(62);
  g_136_spectrum_npf->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  g_136_spectrum_npf->GetYaxis()->SetTitle("dN/dp_{T}");
  g_136_spectrum_npf->GetYaxis()->SetRangeUser(1e-8, 150);
  g_136_spectrum_npf->SetMarkerStyle(20);
  g_136_spectrum_npf->SetMarkerSize(0.8);
  g_136_spectrum_npf->SetMarkerColor(kBlack);

  g_nlo_136_npf->Draw("L,same");
  g_nlo_060_npf->Draw("L, same");

  g_nlo_136_npf->SetLineStyle(7);
  g_nlo_060_npf->SetLineStyle(7);

  g_combined_npf->Draw("P,same");
  f_combined_fit_npf->Draw("same");
  //g_060_spectrum_npf->Draw("P,same");
  //g_136_spectrum_npf->Draw("P,same");
  //g_086_spectrum_npf->Draw("P,same");
  //g_060_spectrum_var1_npf->Draw("P,same");
  //g_060_spectrum_var2_npf->Draw("P,same");
  f_published_060_spectrum_fit_extrapolated_npf->Draw("same");
  //f_spectrum_fit_var1_extrapolated_npf->Draw("same");
  //f_spectrum_fit_var2_extrapolated_npf->Draw("same");
  //f_spectrum_fit_var3_extrapolated_npf->Draw("same");
  //f_spectrum_fit_var4_extrapolated_npf->Draw("same");

  for (int i = 0; i < NPOINTS060; i++)
  {
    //systematicErrors_npf060[i]->Draw("same");
  }

  for (int i = 0; i < NPOINTS086; i++)
  {
    //systematicErrors_npf086[i]->Draw("same");
  }

  for (int i = 0; i < NPOINTS086 + NPOINTS060; i++)
  {
    systematicErrors_npfCombined[i]->Draw("same");
  }

  TLatex latex;
  latex.SetNDC();
  latex.SetTextSize(0.025);
  latex.DrawLatex(.15, .85, "PHENIX Direct Photon Spectrum");
  latex.DrawLatex(.15, .82, "PRL 98, 012002 [PPG060]");
}


void systematicsPhotons()
{
  gStyle->SetOptStat(0);

  makeNLO();
  makePublishedSpectrum060();
  makePublishedSpectrum086();
  makePublishedSpectrum136();
  combinePublications();
  //defineVariation1();
  //defineVariation2();
  //defineVariation3();
  //defineVariation4();
  plotDataPublishedFit();
  plotDataNoPhaseFactor();
}

