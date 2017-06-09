#include <iostream>

using namespace std;

void systematicsRatio()
{
	TF1* f_published_spectrum_fit_extrapolated_npf_eta = new TF1("f_published_spectrum_fit_extrapolated_npf_eta", "1.3*2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_spectrum_fit_extrapolated_npf_eta->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);
	f_published_spectrum_fit_extrapolated_npf_eta->SetLineColor(kBlue);

	TF1 *f_published_spectrum_fit_extrapolated_npf_pizero = new TF1("f_published_spectrum_fit_extrapolated_npf_pizero", "1.3*2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_spectrum_fit_extrapolated_npf_pizero->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_published_spectrum_fit_extrapolated_npf_pizero->SetLineColor(kRed);

	TF1 *f_published_spectrum_fit_extrapolated_npf_photons = new TF1("f_published_spectrum_fit_extrapolated_npf_photons", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_published_spectrum_fit_extrapolated_npf_photons->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);
	f_published_spectrum_fit_extrapolated_npf_photons->SetLineColor(kGreen+2);

	TCanvas *c = new TCanvas("c", "PPG036 Direct Photon Spectrum + Hagedorn Fit", 700, 700);
	c->SetLogy();
	f_published_spectrum_fit_extrapolated_npf_photons->Draw();
	f_published_spectrum_fit_extrapolated_npf_eta->Draw("same");
	f_published_spectrum_fit_extrapolated_npf_pizero->Draw("same");
}