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

const int NPOINTS = 12;

int varColors[4] = {kGreen + 2, kRed, kBlue, kOrange-3};

TGraphErrors *g_spectrum;
TGraphErrors *g_spectrum_npf;
TF1 *f_published_spectrum_fit;
TF1 *f_published_spectrum_fit_extrapolated;
TF1 *f_published_spectrum_fit_extrapolated_npf;
TBox *systematicErrors[NPOINTS];
TBox *systematicErrors_npf[NPOINTS];

//mT-scaled spectrum to guide extrapolation at low pT
TF1 *f_mtscaled_spectrum;
TF1 *f_mtscaled_spectrum_npf;

//Points from published spectrum Dd^3\sigma/dp^3 = (1/2pi pT) dN/dpT
float data_x[NPOINTS] = {2.75,
                         3.25,
                         3.75,
                         4.25,
                         4.75,
                         5.25,
                         5.75,
                         6.5,
                         7.5,
                         8.5,
                         9.5,
                         11
                        };

float data_y[NPOINTS] = {1.30E-3,
                         3.78E-4,
                         1.37E-4,
                         5.49E-5,
                         2.22E-5,
                         1.08E-5,
                         5.66E-6,
                         2.02E-6,
                         6.99E-7,
                         1.81E-7,
                         1.02E-7,
                         2.21E-8
                        };

float err_y_stat[NPOINTS] = {2.25E-4,
                             8.09E-6,
                             3.53E-6,
                             1.81E-6,
                             1.03E-6,
                             6.82E-7,
                             4.39E-7,
                             1.75E-7,
                             9.14E-8,
                             4.04E-8,
                             2.80E-8,
                             8.52E-9
                            };

float err_y_syst[NPOINTS] = {7.91E-6,
                             1.37E-6,
                             4.03E-7,
                             1.47E-7,
                             5.76E-8,
                             2.83E-8,
                             1.52E-8,
                             5.58E-9,
                             2.00E-9,
                             5.35E-10,
                             3.08E-10,
                             1.35E-10
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

//Ratio of different variations relative to central fit
TH1F *h_ratio_var1;
TH1F *h_ratio_var2;
TH1F *h_ratio_var3;
TH1F *h_ratio_var4;

//----------------------
// Functions
//----------------------

void makePublishedSpectrum()
{
	//Define spectrum obtained from mT scaling
	f_mtscaled_spectrum = new TF1("f_mtscaled_spectrum", "TMath::Power(TMath::Sqrt(1 + (0.135/x)*(0.135/x)),-1.0)*0.5*TMath::Sqrt(1 + (0.135/x)*(0.135/x))*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_mtscaled_spectrum->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_mtscaled_spectrum->SetLineColor(kBlack);
	f_mtscaled_spectrum->SetLineStyle(7);

	f_mtscaled_spectrum_npf = new TF1("f_mtscaled_spectrum_npf", "TMath::Power(TMath::Sqrt(1 + (0.135/x)*(0.135/x)),-1.0)*0.5*TMath::Sqrt(1 + (0.135/x)*(0.135/x))*(2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_mtscaled_spectrum_npf->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_mtscaled_spectrum_npf->SetLineColor(kBlack);
	f_mtscaled_spectrum_npf->SetLineStyle(7);

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
	f_published_spectrum_fit->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);
	f_published_spectrum_fit->SetLineColor(varColors[2]);

	f_published_spectrum_fit_extrapolated = new TF1("f_published_spectrum_fit_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_spectrum_fit_extrapolated->SetParameters(f_published_spectrum_fit->GetParameter(0), f_published_spectrum_fit->GetParameter(1), f_published_spectrum_fit->GetParameter(2), f_published_spectrum_fit->GetParameter(3), f_published_spectrum_fit->GetParameter(4));
	f_published_spectrum_fit_extrapolated->SetLineColor(varColors[2]);

	f_published_spectrum_fit_extrapolated_npf = new TF1("f_published_spectrum_fit_extrapolated_npf", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_spectrum_fit_extrapolated_npf->SetParameters(f_published_spectrum_fit->GetParameter(0), f_published_spectrum_fit->GetParameter(1), f_published_spectrum_fit->GetParameter(2), f_published_spectrum_fit->GetParameter(3), f_published_spectrum_fit->GetParameter(4));
	f_published_spectrum_fit_extrapolated_npf->SetLineColor(varColors[2]);
}


void defineVariation1()
{
	//Define the seventh point in the spectrum as a tipping point, and tilt clockwise about that point by an amount proportional to the point's pT
	float data_y_var1[NPOINTS];
	int tippingPointIndex = 9;
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
	g_spectrum_var1_npf->SetMarkerColor(varColors[0]);
	g_spectrum_var1_npf->SetLineColor(varColors[0]);
	g_spectrum_var1_npf->SetMarkerStyle(20);
	g_spectrum_var1_npf->SetMarkerSize(0.7);

	f_spectrum_fit_var1_extrapolated_npf->SetLineColor(varColors[0]);
}


void defineVariation2()
{
	//Define the seventh point in the spectrum as a tipping point, and tilt clockwise about that point by an amount proportional to the point's pT
	float data_y_var2[NPOINTS];
	int tippingPointIndex = 9;
	float pTextreme = data_x[0];
	float pT0 = data_x[tippingPointIndex];

	for (int i = 0; i < NPOINTS; i++)
	{
		float scaling = -1 * (err_y_syst[i] / 2.0) * ((data_x[i] - pT0) / (pTextreme - pT0));
		data_y_var2[i] = data_y[i] + scaling;
	}

	g_spectrum_var2 = new TGraphErrors(NPOINTS, data_x, data_y_var2, err_x, err_y_stat);
	g_spectrum_var2->SetMarkerColor(kBlue);
	g_spectrum_var2->SetLineColor(kBlue);
	g_spectrum_var2->SetMarkerStyle(20);
	g_spectrum_var2->SetMarkerSize(0.7);

	//Fit the data with a modified Hagedorn function
	f_spectrum_fit_var2 = new TF1("f_spectrum_fit_var2", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 3.25, 15.0);
	f_spectrum_fit_var2->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);
	g_spectrum_var2->Fit(f_spectrum_fit_var2, "Q0R");
	f_spectrum_fit_var2->SetLineColor(varColors[1]);

	f_spectrum_fit_var2_extrapolated = new TF1("f_spectrum_fit_var2_extrapolated", "([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_spectrum_fit_var2_extrapolated->SetParameters(f_spectrum_fit_var2->GetParameter(0), f_spectrum_fit_var2->GetParameter(1), f_spectrum_fit_var2->GetParameter(2), f_spectrum_fit_var2->GetParameter(3), f_spectrum_fit_var2->GetParameter(4));
	f_spectrum_fit_var2_extrapolated->SetLineColor(varColors[1]);

	f_spectrum_fit_var2_extrapolated_npf = new TF1("f_spectrum_fit_var2_extrapolated_npf", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_spectrum_fit_var2_extrapolated_npf->SetParameters(f_spectrum_fit_var2->GetParameter(0), f_spectrum_fit_var2->GetParameter(1), f_spectrum_fit_var2->GetParameter(2), f_spectrum_fit_var2->GetParameter(3), f_spectrum_fit_var2->GetParameter(4));
	f_spectrum_fit_var2_extrapolated_npf->SetLineColor(varColors[1]);

	//Spectrum without phase space factor
	float data_y_npf[NPOINTS];

	for (int i = 0; i < NPOINTS; i++)
	{
		data_y_npf[i] = 2 * TMath::Pi() * data_x[i] * data_y_var2[i];
	}

	g_spectrum_var2_npf = new TGraphErrors(NPOINTS, data_x, data_y_npf, err_x, err_y_stat);
	g_spectrum_var2_npf->SetMarkerColor(varColors[1]);
	g_spectrum_var2_npf->SetLineColor(varColors[1]);
	g_spectrum_var2_npf->SetMarkerStyle(20);
	g_spectrum_var2_npf->SetMarkerSize(0.7);
}


void defineVariation3()
{
	//Fit published spectrum with a Tsallis functional form
	f_spectrum_fit_var3 = new TF1("f_spectrum_fit_var3", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 10.0*([1] - 1))*([1]*[2] + 10.0)) * pow(([1]*[2] + TMath::Sqrt(10.0*10.0 + x*x))/([1]*[2]+10.0),-1*[1])", 0, 18.0);
	//f_spectrum_fit_var3->SetParameter(0, 1.3);
	//f_spectrum_fit_var3->SetParameter(1, 5.5);
	//f_spectrum_fit_var3->SetParameter(2, 0.001);
	f_spectrum_fit_var3->SetParameter(0, 50.0);
	f_spectrum_fit_var3->SetParameter(1, 20.5);
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
	f_spectrum_fit_var4->SetParameter(0, 7.86561e-01);
	f_spectrum_fit_var4->SetParameter(1, 9.24328e+00);
	f_spectrum_fit_var4->SetParameter(2, 1.00848e-01);
	//f_spectrum_fit_var3->SetParameter(0, 1.3);
	//f_spectrum_fit_var3->SetParameter(1, 5.5);
	//f_spectrum_fit_var3->SetParameter(2, 0.001);
	g_spectrum->Fit(f_spectrum_fit_var4, "Q0R");

	f_spectrum_fit_var4->SetLineColor(varColors[3]);

	f_spectrum_fit_var4_extrapolated = new TF1("f_spectrum_fit_var4_extrapolated", "([0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1]))", 0.0, 18.0);
	f_spectrum_fit_var4_extrapolated->SetParameters(f_spectrum_fit_var4->GetParameter(0), f_spectrum_fit_var4->GetParameter(1), f_spectrum_fit_var4->GetParameter(2));
	f_spectrum_fit_var4_extrapolated->SetLineColor(varColors[3]);

	f_spectrum_fit_var4_extrapolated_npf = new TF1("f_spectrum_fit_var4_extrapolated_npf", "2*TMath::Pi()*x*(([0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1])))", 0.0, 18.0);
	f_spectrum_fit_var4_extrapolated_npf->SetParameters(f_spectrum_fit_var4->GetParameter(0), f_spectrum_fit_var4->GetParameter(1), f_spectrum_fit_var4->GetParameter(2));
	f_spectrum_fit_var4_extrapolated_npf->SetLineColor(varColors[3]);
}


void getRatios()
{
	h_ratio_var1 = new TH1F("h_ratio_var1", "h_ratio_var1", 100, 0, 20.0);
	h_ratio_var2 = new TH1F("h_ratio_var2", "h_ratio_var2", 100, 0, 20.0);
	h_ratio_var3 = new TH1F("h_ratio_var3", "h_ratio_var3", 100, 0, 20.0);
	h_ratio_var4 = new TH1F("h_ratio_var4", "h_ratio_var4", 100, 0, 20.0);

	for (int i = 0; i < 95; i++)
	{
		float ratio1 = f_spectrum_fit_var1_extrapolated_npf->Eval(h_ratio_var1->GetBinCenter(i + 1)) / f_mtscaled_spectrum_npf->Eval(h_ratio_var1->GetBinCenter(i + 1));
		float ratio2 = f_spectrum_fit_var2_extrapolated_npf->Eval(h_ratio_var2->GetBinCenter(i + 1)) / f_mtscaled_spectrum_npf->Eval(h_ratio_var2->GetBinCenter(i + 1));
		float ratio3 = f_published_spectrum_fit_extrapolated_npf->Eval(h_ratio_var3->GetBinCenter(i + 1)) / f_mtscaled_spectrum_npf->Eval(h_ratio_var3->GetBinCenter(i + 1));
		float ratio4 = f_spectrum_fit_var4_extrapolated_npf->Eval(h_ratio_var4->GetBinCenter(i + 1)) / f_mtscaled_spectrum_npf->Eval(h_ratio_var4->GetBinCenter(i + 1));

		h_ratio_var1->SetBinContent(i, ratio1);
		h_ratio_var2->SetBinContent(i, ratio2);
		h_ratio_var3->SetBinContent(i, ratio3);
		h_ratio_var4->SetBinContent(i, ratio4);

		h_ratio_var1->SetLineColor(varColors[0]);
		h_ratio_var2->SetLineColor(varColors[1]);
		h_ratio_var3->SetLineColor(varColors[2]);
		h_ratio_var4->SetLineColor(varColors[3]);

		h_ratio_var1->SetLineWidth(2);
		h_ratio_var2->SetLineWidth(2);
		h_ratio_var3->SetLineWidth(2);
		h_ratio_var4->SetLineWidth(2);
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
	//f_spectrum_fit_var1->Draw("same");
	//f_spectrum_fit_var2_extrapolated->Draw("same");

	f_spectrum_fit_var3_extrapolated->Draw("same");
	//f_spectrum_fit_var4->Draw("same");

	f_mtscaled_spectrum->Draw("same");

	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors[i]->Draw("same");
	}
}


void plotDataNoPhaseFactor()
{
	TCanvas *cNP = new TCanvas("cNP", "Stacked Representation of Fit", 700, 900);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
	pad1->SetLogy();
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();

	TH1F *hTemplate2 = new TH1F("hTemplate2", "hTemplate2", 100, 0, 20);
	hTemplate2->SetTitle("");
	hTemplate2->GetXaxis()->SetTitleFont(62);
	hTemplate2->GetXaxis()->SetLabelFont(62);
	hTemplate2->GetXaxis()->SetRangeUser(0, 18);
	hTemplate2->GetYaxis()->SetTitleFont(62);
	hTemplate2->GetYaxis()->SetLabelFont(62);
	hTemplate2->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate2->GetYaxis()->SetTitle("dN/dp_{T}");
	hTemplate2->GetYaxis()->SetTitleOffset(1.3);
	hTemplate2->GetYaxis()->SetRangeUser(2e-8, 200);
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
	f_spectrum_fit_var4_extrapolated_npf->Draw("same");

	f_mtscaled_spectrum_npf->Draw("same");

	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors_npf[i]->Draw("same");
	}

	g_spectrum_npf->Draw("P,same");

	TLatex latex;
	latex.SetNDC();
	latex.SetTextSize(0.025);
	latex.DrawLatex(.15, .85, "PHENIX Eta Spectrum");
	latex.DrawLatex(.15, .82, "PRD 76, 051106(R) [PPG063]");

	TLegend *legend = new TLegend(0.45, 0.45, 0.88, 0.65);
	legend->AddEntry(f_mtscaled_spectrum_npf, "Fit to Spectrum (based on m_{T} scaling)", "l");
	legend->AddEntry(f_published_spectrum_fit_extrapolated_npf, "Variation 1: Modified Hagedorn", "l");
	legend->AddEntry(f_spectrum_fit_var1_extrapolated_npf, "Variation 2: Counterclockwise Tilt", "l");
	legend->AddEntry(f_spectrum_fit_var2_extrapolated_npf, "Variation 3: Clockwise Tilt", "l");
	legend->AddEntry(f_spectrum_fit_var4_extrapolated_npf, "Variation 4: Tsallis", "l");
	legend->SetFillStyle(0.0);
	legend->SetLineColor(kWhite);
	legend->Draw("same");

	cNP->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();
	h_ratio_var1->SetMarkerStyle(7);
	h_ratio_var1->GetXaxis()->SetTitleFont(62);
	h_ratio_var1->GetXaxis()->SetLabelFont(62);
	h_ratio_var1->GetXaxis()->SetRangeUser(0, 18);
	h_ratio_var1->GetYaxis()->SetTitleFont(62);
	h_ratio_var1->GetYaxis()->SetLabelFont(62);
	h_ratio_var1->SetTitle("");
	h_ratio_var1->GetYaxis()->CenterTitle();
	h_ratio_var1->GetYaxis()->SetRangeUser(0.0, 4.8);
	h_ratio_var1->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_ratio_var1->GetYaxis()->SetTitle("Variation / Fit");
	h_ratio_var1->GetYaxis()->SetTitleSize(0.085);
	h_ratio_var1->GetYaxis()->SetTitleOffset(0.4);
	h_ratio_var1->GetYaxis()->SetLabelSize(0.085);
	h_ratio_var1->GetXaxis()->SetTitleSize(0.085);
	h_ratio_var1->GetXaxis()->SetTitleOffset(1.2);
	h_ratio_var1->GetXaxis()->SetLabelSize(0.085);
	h_ratio_var1->Draw("L");
	h_ratio_var2->Draw("L,same");
	h_ratio_var3->Draw("L,same");
	h_ratio_var4->Draw("L,same");

	cNP->cd();
}


void systematicsEtas()
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

