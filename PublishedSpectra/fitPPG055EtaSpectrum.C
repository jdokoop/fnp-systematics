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

const int NPOINTS = 12;

TGraphErrors *g_spectrum;

//-------------------------
// Functions
//-------------------------

void fitPPG055EtaSpectrum()
{
	//Initialize spectrum with data-thieved points
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

	float error_x[NPOINTS] = {0};

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

	g_spectrum = new TGraphErrors(NPOINTS, data_x, data_y);

	//Fit the data with a Hagedorn function
	TF1 *fit = new TF1("fit", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 2.735455021, 11.0158413);
	fit->SetParameters(545, -0.109245, 0.0, 0.587309, 8.05);
	g_spectrum->Fit(fit, "0R");

	//Prepare array of TBoxes with systematic errors
	TBox *systematicErrors[NPOINTS];
	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors[i] = new TBox(data_x[i] - 0.15, data_y[i] - (err_y_syst[i] / 2.0), data_x[i] + 0.15, data_y[i] + (err_y_syst[i] / 2.0));
		systematicErrors[i]->SetLineColor(kBlack);
		systematicErrors[i]->SetFillStyle(0);
	}

	//Draw data+fit
	TCanvas *c = new TCanvas("c", "PPG055 Eta Spectrum + Hagedorn Fit", 700, 700);
	c->SetLogy();
	g_spectrum->SetTitle("");
	g_spectrum->GetXaxis()->SetTitleFont(62);
	g_spectrum->GetXaxis()->SetLabelFont(62);
	g_spectrum->GetYaxis()->SetTitleFont(62);
	g_spectrum->GetYaxis()->SetLabelFont(62);
	g_spectrum->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_spectrum->GetYaxis()->SetTitle("E#frac{d^{3}#sigma}{dp^{3}}");
	g_spectrum->GetYaxis()->SetRangeUser(1e-10, 1e2);
	g_spectrum->SetMarkerStyle(20);
	g_spectrum->SetMarkerColor(kBlack);
	g_spectrum->Draw("AP");
	g_spectrum->GetFunction("fit")->Draw("same");

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

	//c->SaveAs("EtaFitPPG055.pdf");
}