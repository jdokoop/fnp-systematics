//-------------------------------------------------------
// Take fit to the direct photon spectrum from PPG060
// and tweak it to account for systematic uncertainties
//
// J. Orjuela-Koop
// May. 2017
//-------------------------------------------------------

//----------------------
// Variables
//----------------------

const int NPOINTS = 25;

TGraphErrors *g_spectrum;
TGraphErrors *g_spectrum_npf;
TF1 *f_published_spectrum_fit;
TF1 *f_published_spectrum_fit_extrapolated;
TF1 *f_published_spectrum_fit_extrapolated_npf;
TBox *systematicErrors[NPOINTS];
TBox *systematicErrors_npf[NPOINTS];

//Points from published spectrum Dd^3\sigma/dp^3 = (1/2pi pT) dN/dpT
float data_x[NPOINTS] = {0.616,
                         0.866,
                         1.215,
                         1.719,
                         2.223,
                         2.725,
                         3.228,
                         3.730,
                         4.231,
                         4.732,
                         5.234,
                         5.735,
                         6.237,
                         6.738,
                         7.238,
                         7.739,
                         8.240,
                         8.740,
                         9.241,
                         9.741,
                         10.880,
                         12.900,
                         14.910,
                         16.920,
                         18.930
                        };

float data_y[NPOINTS] = {5.953e+00,
                         1.783e+00,
                         3.960e-01,
                         6.148e-02,
                         1.300e-02,
                         3.285e-03,
                         9.786e-04,
                         3.280e-04,
                         1.190e-04,
                         4.887e-05,
                         2.164e-05,
                         1.038e-05,
                         5.191e-06,
                         2.748e-06,
                         1.564e-06,
                         9.132e-07,
                         5.453e-07,
                         3.474e-07,
                         2.135e-07,
                         1.332e-07,
                         5.916e-08,
                         1.494e-08,
                         3.699e-09,
                         1.248e-09,
                         5.470e-10
                        };

float err_y_stat[NPOINTS] = {3.075e-03,
                             9.725e-04,
                             2.780e-04,
                             6.855e-05,
                             2.453e-05,
                             1.028e-05,
                             4.934e-06,
                             2.547e-06,
                             2.236e-07,
                             1.329e-07,
                             8.269e-08,
                             5.426e-08,
                             3.646e-08,
                             2.542e-08,
                             1.826e-08,
                             1.340e-08,
                             9.973e-09,
                             7.586e-09,
                             5.801e-09,
                             4.511e-09,
                             1.433e-09,
                             6.880e-10,
                             3.018e-10,
                             1.652e-10,
                             1.142e-10
                            };

float err_y_syst[NPOINTS] = {7.516e-01,
                             1.299e-01,
                             2.713e-02,
                             4.562e-03,
                             1.026e-03,
                             2.689e-04,
                             8.193e-05,
                             2.786e-05,
                             1.094e-05,
                             4.518e-06,
                             2.011e-06,
                             9.687e-07,
                             4.863e-07,
                             2.585e-07,
                             1.479e-07,
                             8.683e-08,
                             5.216e-08,
                             3.348e-08,
                             2.067e-08,
                             1.296e-08,
                             6.023e-09,
                             1.661e-09,
                             4.493e-10,
                             1.736e-10,
                             8.630e-11
                            };

float err_x[NPOINTS] = {0};

//Variations on points
TGraphErrors *g_spectrum_var1;
TGraphErrors *g_spectrum_var1_npf;
TF1 *f_spectrum_fit_var1;
TF1 *f_spectrum_fit_var1_extrapolated;
TF1 *f_spectrum_fit_var1_extrapolated_npf;

TGraphErrors *g_spectrum_var2;
TGraphErrors *g_spectrum_var2_npf;
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

TH1F *h_ratio_var1;
TH1F *h_ratio_var2;

//----------------------
// Functions
//----------------------

void makePublishedSpectrum()
{
	//Make spectrum with no phase space factor
	float data_y_npf[NPOINTS];
	float err_y_stat_npf[NPOINTS];
	for (int i = 0; i < NPOINTS; i++)
	{
		data_y_npf[i] = 2 * TMath::Pi() * data_x[i] * data_y[i];
		err_y_stat_npf[i] = 2 * TMath::Pi() * data_x[i] * err_y_stat[i];
	}

	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors_npf[i] = new TBox(data_x[i] - 0.15, data_y_npf[i] - (2 * TMath::Pi() * data_x[i]*err_y_syst[i] / 2.0), data_x[i] + 0.15, data_y_npf[i] + (2 * TMath::Pi() * data_x[i]*err_y_syst[i] / 2.0));
		systematicErrors_npf[i]->SetLineColor(kBlack);
		systematicErrors_npf[i]->SetFillStyle(0);
	}

	//Make normal spectrum
	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors[i] = new TBox(data_x[i] - 0.15, data_y[i] - (err_y_syst[i] / 2.0), data_x[i] + 0.15, data_y[i] + (err_y_syst[i] / 2.0));
		systematicErrors[i]->SetLineColor(kBlack);
		systematicErrors[i]->SetFillStyle(0);
	}

	g_spectrum = new TGraphErrors(NPOINTS, data_x, data_y, err_x, err_y_stat);
	g_spectrum_npf = new TGraphErrors(NPOINTS, data_x, data_y_npf, err_x, err_y_stat_npf);

	f_published_spectrum_fit = new TF1("f_published_spectrum_fit", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 3.25, 15.0);
	f_published_spectrum_fit->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_published_spectrum_fit->SetLineColor(kBlack);

	f_published_spectrum_fit_extrapolated = new TF1("f_published_spectrum_fit_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_spectrum_fit_extrapolated->SetParameters(f_published_spectrum_fit->GetParameter(0), f_published_spectrum_fit->GetParameter(1), f_published_spectrum_fit->GetParameter(2), f_published_spectrum_fit->GetParameter(3), f_published_spectrum_fit->GetParameter(4));
	f_published_spectrum_fit_extrapolated->SetLineColor(kBlack);

	f_published_spectrum_fit_extrapolated_npf = new TF1("f_published_spectrum_fit_extrapolated_npf", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_spectrum_fit_extrapolated_npf->SetParameters(f_published_spectrum_fit->GetParameter(0), f_published_spectrum_fit->GetParameter(1), f_published_spectrum_fit->GetParameter(2), f_published_spectrum_fit->GetParameter(3), f_published_spectrum_fit->GetParameter(4));
	f_published_spectrum_fit_extrapolated_npf->SetLineColor(kBlack);
}


void defineVariation1()
{
	//Define the seventh point in the spectrum as a tipping point, and tilt clockwise about that point by an amount proportional to the point's pT
	float data_y_var1[NPOINTS];
	int tippingPointIndex = 8;
	float pTextreme = data_x[0];
	float pT0 = data_x[tippingPointIndex];

	for (int i = 0; i < NPOINTS; i++)
	{
		float scaling = (err_y_syst[i] / 2.0) * ((data_x[i] - pT0) / (pTextreme - pT0));
		data_y_var1[i] = data_y[i] + scaling;
	}

	g_spectrum_var1 = new TGraphErrors(NPOINTS, data_x, data_y_var1, err_x, err_y_stat);
	g_spectrum_var1->SetMarkerColor(kRed);
	g_spectrum_var1->SetLineColor(kRed);
	g_spectrum_var1->SetMarkerStyle(20);
	g_spectrum_var1->SetMarkerSize(0.7);

	//Fit the data with a modified Hagedorn function
	f_spectrum_fit_var1 = new TF1("f_spectrum_fit_var1", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0.0, 18.0);
	f_spectrum_fit_var1->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);
	g_spectrum_var1->Fit(f_spectrum_fit_var1, "Q0R");

	f_spectrum_fit_var1_extrapolated = new TF1("f_spectrum_fit_var1_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_spectrum_fit_var1_extrapolated->SetParameters(f_spectrum_fit_var1->GetParameter(0), f_spectrum_fit_var1->GetParameter(1), f_spectrum_fit_var1->GetParameter(2), f_spectrum_fit_var1->GetParameter(3), f_spectrum_fit_var1->GetParameter(4));

	f_spectrum_fit_var1_extrapolated_npf = new TF1("f_spectrum_fit_var1_extrapolated_npf", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_spectrum_fit_var1_extrapolated_npf->SetParameters(f_spectrum_fit_var1->GetParameter(0), f_spectrum_fit_var1->GetParameter(1), f_spectrum_fit_var1->GetParameter(2), f_spectrum_fit_var1->GetParameter(3), f_spectrum_fit_var1->GetParameter(4));

	//Spectrum without phase space factor
	float data_y_npf[NPOINTS];

	for (int i = 0; i < NPOINTS; i++)
	{
		data_y_npf[i] = 2 * TMath::Pi() * data_x[i] * data_y_var1[i];
	}

	g_spectrum_var1_npf = new TGraphErrors(NPOINTS, data_x, data_y_npf, err_x, err_y_stat);
	g_spectrum_var1_npf->SetMarkerColor(kRed);
	g_spectrum_var1_npf->SetLineColor(kRed);
	g_spectrum_var1_npf->SetMarkerStyle(20);
	g_spectrum_var1_npf->SetMarkerSize(0.7);
}


void defineVariation2()
{
	//Define the seventh point in the spectrum as a tipping point, and tilt clockwise about that point by an amount proportional to the point's pT
	float data_y_var2[NPOINTS];
	int tippingPointIndex = 8;
	float pTextreme = data_x[0];
	float pT0 = data_x[tippingPointIndex];

	for (int i = 0; i < NPOINTS; i++)
	{
		float scaling = -1*(err_y_syst[i] / 2.0) * ((data_x[i] - pT0) / (pTextreme - pT0));
		data_y_var2[i] = data_y[i] + scaling;
	}

	g_spectrum_var2 = new TGraphErrors(NPOINTS, data_x, data_y_var2, err_x, err_y_stat);
	g_spectrum_var2->SetMarkerColor(kBlue);
	g_spectrum_var2->SetLineColor(kBlue);
	g_spectrum_var2->SetMarkerStyle(20);
	g_spectrum_var2->SetMarkerSize(0.7);

	//Fit the data with a modified Hagedorn function
	f_spectrum_fit_var2 = new TF1("f_spectrum_fit_var2", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0.0, 18.0);
	f_spectrum_fit_var2->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);
	g_spectrum_var2->Fit(f_spectrum_fit_var2, "Q0R");

	f_spectrum_fit_var2_extrapolated = new TF1("f_spectrum_fit_var2_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_spectrum_fit_var2_extrapolated->SetParameters(f_spectrum_fit_var2->GetParameter(0), f_spectrum_fit_var2->GetParameter(1), f_spectrum_fit_var2->GetParameter(2), f_spectrum_fit_var2->GetParameter(3), f_spectrum_fit_var2->GetParameter(4));

	f_spectrum_fit_var2_extrapolated_npf = new TF1("f_spectrum_fit_var2_extrapolated_npf", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_spectrum_fit_var2_extrapolated_npf->SetParameters(f_spectrum_fit_var2->GetParameter(0), f_spectrum_fit_var2->GetParameter(1), f_spectrum_fit_var2->GetParameter(2), f_spectrum_fit_var2->GetParameter(3), f_spectrum_fit_var2->GetParameter(4));

	//Spectrum without phase space factor
	float data_y_npf[NPOINTS];

	for (int i = 0; i < NPOINTS; i++)
	{
		data_y_npf[i] = 2 * TMath::Pi() * data_x[i] * data_y_var2[i];
	}

	g_spectrum_var2_npf = new TGraphErrors(NPOINTS, data_x, data_y_npf, err_x, err_y_stat);
	g_spectrum_var2_npf->SetMarkerColor(kBlue);
	g_spectrum_var2_npf->SetLineColor(kBlue);
	g_spectrum_var2_npf->SetMarkerStyle(20);
	g_spectrum_var2_npf->SetMarkerSize(0.7);

	f_spectrum_fit_var2_extrapolated_npf->SetLineColor(kBlue);
}


void defineVariation3()
{
	//Fit published spectrum with a Tsallis functional form
	f_spectrum_fit_var3 = new TF1("f_spectrum_fit_var3", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 10.0*([1] - 1))*([1]*[2] + 10.0)) * pow(([1]*[2] + TMath::Sqrt(10.0*10.0 + x*x))/([1]*[2]+10.0),-1*[1])", 0, 18.0);
	f_spectrum_fit_var3->SetParameter(0, 1.3);
	f_spectrum_fit_var3->SetParameter(1, 5.5);
	f_spectrum_fit_var3->SetParameter(2, 0.001);
	g_spectrum->Fit(f_spectrum_fit_var3, "Q0R");
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
	g_spectrum->Fit(f_spectrum_fit_var4, "Q0R");
	f_spectrum_fit_var4->SetLineColor(kSpring - 6);

	f_spectrum_fit_var4_extrapolated = new TF1("f_spectrum_fit_var4_extrapolated", "([0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1]))", 0.0, 18.0);
	f_spectrum_fit_var4_extrapolated->SetParameters(f_spectrum_fit_var4->GetParameter(0), f_spectrum_fit_var4->GetParameter(1), f_spectrum_fit_var4->GetParameter(2));
	f_spectrum_fit_var4_extrapolated->SetLineColor(kSpring - 6);

	f_spectrum_fit_var4_extrapolated_npf = new TF1("f_spectrum_fit_var4_extrapolated_npf", "2*TMath::Pi()*x*(([0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1])))", 0.0, 18.0);
	f_spectrum_fit_var4_extrapolated_npf->SetParameters(f_spectrum_fit_var4->GetParameter(0), f_spectrum_fit_var4->GetParameter(1), f_spectrum_fit_var4->GetParameter(2));
	f_spectrum_fit_var4_extrapolated_npf->SetLineColor(kSpring - 6);
}


void getRatios()
{
	h_ratio_var1 = new TH1F("h_ratio_var1", "h_ratio_var1", 100, 0, 20.0);
	h_ratio_var2 = new TH1F("h_ratio_var2", "h_ratio_var2", 100, 0, 20.0);

	for (int i = 0; i < 95; i++)
	{
		float ratio1 = f_spectrum_fit_var1_extrapolated_npf->Eval(h_ratio_var1->GetBinCenter(i + 1)) / f_published_spectrum_fit_extrapolated_npf->Eval(h_ratio_var1->GetBinCenter(i + 1));
		float ratio2 = f_spectrum_fit_var2_extrapolated_npf->Eval(h_ratio_var2->GetBinCenter(i + 1)) / f_published_spectrum_fit_extrapolated_npf->Eval(h_ratio_var2->GetBinCenter(i + 1));

		h_ratio_var1->SetBinContent(i, ratio1);
		h_ratio_var2->SetBinContent(i, ratio2);
	
		h_ratio_var1->SetLineColor(kRed);
		h_ratio_var2->SetLineColor(kBlue);
	
		h_ratio_var1->SetLineWidth(2);
		h_ratio_var2->SetLineWidth(2);
	}

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


	g_spectrum->SetTitle("");
	g_spectrum->GetXaxis()->SetTitleFont(62);
	g_spectrum->GetXaxis()->SetLabelFont(62);
	g_spectrum->GetXaxis()->SetRangeUser(0, 18);
	g_spectrum->GetYaxis()->SetTitleFont(62);
	g_spectrum->GetYaxis()->SetLabelFont(62);
	g_spectrum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_spectrum->GetYaxis()->SetTitle("(1/2#pi p_{T}) dN/dp_{T}");
	g_spectrum->GetYaxis()->SetRangeUser(1e-10, 5e-5);
	g_spectrum->SetMarkerStyle(20);
	g_spectrum->SetMarkerSize(0.8);
	g_spectrum->SetMarkerColor(kBlack);
	g_spectrum->Draw("P,same");
	//g_spectrum_var1->Draw("P,same");
	//g_spectrum_var2->Draw("P,same");
	f_published_spectrum_fit_extrapolated->Draw("same");
	f_spectrum_fit_var1->Draw("same");
	//f_spectrum_fit_var2_extrapolated->Draw("same");

	//f_spectrum_fit_var3_extrapolated->Draw("same");
	//f_spectrum_fit_var4->Draw("same");

	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors[i]->Draw("same");
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

	g_spectrum_npf->SetTitle("");
	g_spectrum_npf->GetXaxis()->SetTitleFont(62);
	g_spectrum_npf->GetXaxis()->SetLabelFont(62);
	g_spectrum_npf->GetXaxis()->SetRangeUser(0, 18);
	g_spectrum_npf->GetYaxis()->SetTitleFont(62);
	g_spectrum_npf->GetYaxis()->SetLabelFont(62);
	g_spectrum_npf->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_spectrum_npf->GetYaxis()->SetTitle("dN/dp_{T}");
	g_spectrum_npf->GetYaxis()->SetRangeUser(1e-8, 150);
	g_spectrum_npf->SetMarkerStyle(20);
	g_spectrum_npf->SetMarkerSize(0.8);
	g_spectrum_npf->SetMarkerColor(kBlack);
	g_spectrum_npf->Draw("P,same");
	//g_spectrum_var1_npf->Draw("P,same");
	//g_spectrum_var2_npf->Draw("P,same");
	f_published_spectrum_fit_extrapolated_npf->Draw("same");
	f_spectrum_fit_var1_extrapolated_npf->Draw("same");
	f_spectrum_fit_var2_extrapolated_npf->Draw("same");
	//f_spectrum_fit_var3_extrapolated_npf->Draw("same");
	//f_spectrum_fit_var4_extrapolated_npf->Draw("same");

	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors_npf[i]->Draw("same");
	}

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.025);
	latex.DrawLatex(.15, .85, "PHENIX Pizero Spectrum");
	latex.DrawLatex(.15, .82, "PRD 76, 051106(R) [PPG063]");
}


void systematicsPizeros()
{
	gStyle->SetOptStat(0);

	makePublishedSpectrum();
	defineVariation1();
	defineVariation2();
	defineVariation3();
	defineVariation4();
	getRatios();
	//plotDataPublishedFit();
	plotDataNoPhaseFactor();
}

