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

const int NPOINTS = 25;

TGraphErrors *g_spectrum;

//-------------------------
// Functions
//-------------------------

void fitPPG063PizeroSpectrum()
{
	//Initialize spectrum with data-thieved points
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

	g_spectrum = new TGraphErrors(NPOINTS, data_x, data_y, err_x, err_y_stat);

	//Fit the data with a Hagedorn function
	TF1 *fit = new TF1("fit", "[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0.5, 19.0);
	fit->SetParameters(545.1, -0.109245, 0.0797, 0.5951, 7.8939);
	g_spectrum->Fit(fit, "0R");

	//Get ratio data/fit
	for (int i = 0; i < NPOINTS; i++)
	{
		double x, y;
		g_spectrum->GetPoint(i, x, y);
		cout << y / fit->Eval(x) << endl;
	}

	//Prepare array of TBoxes with systematic errors
	TBox *systematicErrors[NPOINTS];
	for (int i = 0; i < NPOINTS; i++)
	{
		systematicErrors[i] = new TBox(data_x[i] - 0.15, data_y[i] - (err_y_syst[i] / 2.0), data_x[i] + 0.15, data_y[i] + (err_y_syst[i] / 2.0));
		systematicErrors[i]->SetLineColor(kBlack);
		systematicErrors[i]->SetFillStyle(0);
	}

	//Draw data+fit
	TCanvas *c = new TCanvas("c", "PPG063 Pizero Spectrum + Hagedorn Fit", 700, 700);
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

	//c->SaveAs("PizeroFitPPG063.pdf");
}