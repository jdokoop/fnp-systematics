//------------------------------------------------------------------
// - Take fit to pizero dN/dpT
// - Use chain rule to get dN/dmT for pizeros
// - Multiply by ~ 0.5 to get dN/dmT for etas, according to AN
// - Convert back to dN/dpT for etas
//
// Published dN/dmT and ratios available in PHENIX AN:
// https://www.phenix.bnl.gov/phenix/WWW/p/info/an/1041/an-07.pdf
//------------------------------------------------------------------

//--------------------------
// Variables
//--------------------------

//Fit to pizero dN/dpT
TF1 *f_pizero_dNdpT;

//Pizero dN/dmT
TF1 *f_pizero_dNdmT;

//Eta dN/dmT
TF1 *f_eta_dNdmT;

//Eta dNdpT
TF1 *f_eta_dNdpT;

//--------------------------
// Functions
//--------------------------
void mTscalingPizeros()
{
	//Define fit to published pizero dNdpT
	f_pizero_dNdpT = new TF1("f_pizero_dNdpT", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_pizero_dNdpT->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_pizero_dNdpT->SetLineColor(kRed);

	//Construct pizero dNdmT
	// dN/dmT = \sqrt{1 + (m/p_{T})^2} \times dN/dpT
	f_pizero_dNdmT = new TF1("f_pizero_dNdmT", "TMath::Sqrt(1 + (0.135/x)*(0.135/x))*(2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_pizero_dNdmT->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_pizero_dNdmT->SetLineColor(kBlue);

	//Construct eta dNdmT = 0.5 * pizero dNdmT
	f_eta_dNdmT = new TF1("f_eta_dNdmT", "0.5*TMath::Sqrt(1 + (0.135/x)*(0.135/x))*(2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_eta_dNdmT->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_eta_dNdmT->SetLineColor(kBlue);

	//Convert back to eta dNdpT
	// dN/dpT = [\sqrt{1 + (m/p_{T})^2}]^{-1} \times dN/dmT
	f_eta_dNdpT = new TF1("f_eta_dNdpT", "TMath::Power(TMath::Sqrt(1 + (0.135/x)*(0.135/x)),-1.0)*0.5*TMath::Sqrt(1 + (0.135/x)*(0.135/x))*(2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_eta_dNdpT->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_eta_dNdpT->SetLineColor(kGreen + 2);

	TCanvas *c = new TCanvas("c", "c", 600, 600);
	f_pizero_dNdpT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	f_pizero_dNdpT->GetYaxis()->SetTitle("dN/dp_{T}");
	f_pizero_dNdpT->Draw();
	f_pizero_dNdmT->Draw("same");
	f_eta_dNdpT->Draw("same");
}