//-------------------------------------------------------
// Fit the direct photon spectrum from PPG036
//
// J. Orjuela-Koop
// Feb. 2017
//-------------------------------------------------------

#include <iostream>

using namespace std;

//-------------------------
// Variables
//-------------------------

const int NPOINTS = 17;

TGraphErrors *g_spectrum;

//-------------------------
// Functions
//-------------------------

void fitPPG060PhotonSpectrum()
{
	//Initialize spectrum with data from PPG060 [data table from published paper]
	float data_x[NPOINTS] = {3.25,
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

	float data_y[NPOINTS] = {2.22E-05,
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

	float err_y_stat[NPOINTS] = {1.10E-06,
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

	float err_y_syst[NPOINTS] = {2.54E-05,
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

	float err_x[NPOINTS] = {0};

	g_spectrum = new TGraphErrors(NPOINTS, data_x, data_y, err_x, err_y_stat);

	//Fit the data with a Hagedorn function
	TF1 *fit = new TF1("fit", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 3.25, 15.0);
	fit->SetParameters(41.0, -0.89, 0.35, 0.67, 8.15);
	g_spectrum->Fit(fit, "Q0R");

	//Fit the data with a Tsallis function
	TF1 *tsallisFit = new TF1("tsallisFit", "[0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.140*([1] - 1))*([1]*[2] + 0.140)) * pow(([1]*[2] + TMath::Sqrt(0.140*0.140 + x*x))/([1]*[2]+0.140),-1*[1])", 3.25, 15.0);
	tsallisFit->SetParameter(0, 1.3);
	tsallisFit->SetParameter(1, 5.5);
	tsallisFit->SetParameter(2, 0.001);
	g_spectrum->Fit(tsallisFit, "Q0R");

	//Prepare array of TBoxes with systematic errors
	TBox *systematicErrors[NPOINTS];
	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors[i] = new TBox(data_x[i] - 0.15, data_y[i] - (err_y_syst[i] / 2.0), data_x[i] + 0.15, data_y[i] + (err_y_syst[i] / 2.0));
		systematicErrors[i]->SetLineColor(kBlack);
		systematicErrors[i]->SetFillStyle(0);
	}

	//Draw data+fit
	TCanvas *c = new TCanvas("c", "PPG036 Direct Photon Spectrum + Hagedorn Fit", 700, 700);
	c->SetLogy();
	g_spectrum->SetTitle("");
	g_spectrum->GetXaxis()->SetTitleFont(62);
	g_spectrum->GetXaxis()->SetLabelFont(62);
	g_spectrum->GetYaxis()->SetTitleFont(62);
	g_spectrum->GetYaxis()->SetLabelFont(62);
	g_spectrum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_spectrum->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}}");
	g_spectrum->GetYaxis()->SetRangeUser(1e-10, 5e-5);
	g_spectrum->SetMarkerStyle(20);
	g_spectrum->SetMarkerSize(0.8);
	g_spectrum->SetMarkerColor(kBlack);
	g_spectrum->Draw("AP");
	g_spectrum->GetFunction("tsallisFit")->Draw("same");

	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors[i]->Draw("same");
	}

	TLatex *tParA = new TLatex(.6, .8, Form("A = %g", fit->GetParameter(0)));
	tParA->SetNDC(true);
	tParA->SetTextAlign(22);
	tParA->SetTextColor(kRed + 2);
	tParA->Draw("same");

	TLatex *tParB = new TLatex(.6, .75, Form("B = %g", fit->GetParameter(1)));
	tParB->SetNDC(true);
	tParB->SetTextAlign(22);
	tParB->SetTextColor(kRed + 2);
	tParB->Draw("same");

	TLatex *tParC = new TLatex(.6, .7, Form("C = %g", fit->GetParameter(2)));
	tParC->SetNDC(true);
	tParC->SetTextAlign(22);
	tParC->SetTextColor(kRed + 2);
	tParC->Draw("same");

	TLatex *tParD = new TLatex(.6, .65, Form("D = %g", fit->GetParameter(3)));
	tParD->SetNDC(true);
	tParD->SetTextAlign(22);
	tParD->SetTextColor(kRed + 2);
	tParD->Draw("same");

	TLatex *tParE = new TLatex(.6, .6, Form("E = %g", fit->GetParameter(4)));
	tParE->SetNDC(true);
	tParE->SetTextAlign(22);
	tParE->SetTextColor(kRed + 2);
	tParE->Draw("same");

	//c->SaveAs("PhotonFitPPG060.pdf");
}