//--------------------------------------------------------------------
// Catalog of fits to each published spectrum, plus all corresponding
// variations
//--------------------------------------------------------------------

//------------------------------
// Variables
//------------------------------

//Etas
TF1 *f_eta_spectrum[5];

//Pizeros
TF1 *f_pizero_spectrum[3];

//Photons
TF1 *f_photon_spectrum[4];

//Colors
int plotColors[5] = {kBlack, kBlue, kRed, kGreen + 2, kOrange - 3};

//------------------------------
// Functions
//------------------------------

void plotPizero()
{
	TCanvas *c1 = new TCanvas("c1", "Pizeros", 600, 600);
	c1->SetLogy();
	f_pizero_spectrum[0]->Draw();
	f_pizero_spectrum[1]->Draw("same");
	f_pizero_spectrum[2]->Draw("same");
}

void plotPhoton()
{
	TCanvas *c2 = new TCanvas("c2", "Photons", 600, 600);
	c2->SetLogy();
	f_photon_spectrum[0]->Draw();
	f_photon_spectrum[1]->Draw("same");
	f_photon_spectrum[2]->Draw("same");
	f_photon_spectrum[3]->Draw("same");
}

void plotEta()
{
	TCanvas *c3 = new TCanvas("c3", "Etas", 600, 600);
	c3->SetLogy();
	f_eta_spectrum[0]->Draw();
	f_eta_spectrum[1]->Draw("same");
	f_eta_spectrum[2]->Draw("same");
	f_eta_spectrum[3]->Draw("same");
	f_eta_spectrum[4]->Draw("same");
}

void spectraCatalog()
{
	f_photon_spectrum[0] = new TF1("f_photon_spectrum_0", "2*TMath::Pi()*x*[p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4])", 0.0, 18.0);
	f_photon_spectrum[0]->SetParameters(0.242345, -0.0827585, 0.00918447, 4.13943, 13.6974);
	f_photon_spectrum[0]->SetLineColor(plotColors[0]);

	f_photon_spectrum[1] = new TF1("f_photon_spectrum_0", "2*TMath::Pi()*x*([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4]))", 0.0, 18.0);
	f_photon_spectrum[1]->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);
	f_photon_spectrum[1]->SetLineColor(plotColors[1]);

	f_photon_spectrum[2] = new TF1("f_photon_spectrum_0", "2*TMath::Pi()*x*(([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4])))", 0.0, 18.0);
	f_photon_spectrum[2]->SetParameters(0.171216, -0.0312528, 0.00442746, 11.6472, 29.538);
	f_photon_spectrum[2]->SetLineColor(plotColors[2]);

	f_photon_spectrum[3] = new TF1("f_photon_spectrum_0", "2*TMath::Pi()*x*(([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4])))", 0.0, 18.0);
	f_photon_spectrum[3]->SetParameters(0.012215, -0.0539112, 0.00527051, 4.69066, 11.7111);
	f_photon_spectrum[3]->SetLineColor(plotColors[3]);


	f_pizero_spectrum[0] = new TF1("f_pizero_spectrum_0", "2*TMath::Pi()*x*([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4]))", 0.0, 18.0);
	f_pizero_spectrum[0]->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_pizero_spectrum[0]->SetLineColor(plotColors[0]);

	f_pizero_spectrum[1] = new TF1("f_pizero_spectrum_0", "2*TMath::Pi()*x*(([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4])))", 0.0, 18.0);
	f_pizero_spectrum[1]->SetParameters(339.852, 0.417738, 0.0599019, 0.719301, 8.33941);
	f_pizero_spectrum[1]->SetLineColor(plotColors[1]);

	f_pizero_spectrum[2] = new TF1("f_pizero_spectrum_0", "2*TMath::Pi()*x*(([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4])))", 0.0, 18.0);
	f_pizero_spectrum[2]->SetParameters(184.594, 0.524929, 0.0111792, 0.762548, 8.24103);
	f_pizero_spectrum[2]->SetLineColor(plotColors[2]);


	f_eta_spectrum[0] = new TF1("f_eta_spectrum_0", "TMath::Power(TMath::Sqrt(1+(0.135/x)*(0.135/x)),-1.0)*0.5*TMath::Sqrt(1+(0.135/x)*(0.135/x))*(2*TMath::Pi()*x*([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4])))", 0.0, 18.0);
	f_eta_spectrum[0]->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);
	f_eta_spectrum[0]->SetLineColor(plotColors[0]);

	f_eta_spectrum[1] = new TF1("f_eta_spectrum_0", "2*TMath::Pi()*x*(([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4])))", 0.0, 18.0);
	f_eta_spectrum[1]->SetParameters(1201.68, -0.0135258, 0.0215599, 0.782018, 9.31206);
	f_eta_spectrum[1]->SetLineColor(plotColors[1]);

	f_eta_spectrum[2] = new TF1("f_eta_spectrum_0", "2*TMath::Pi()*x*([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4]))", 0.0, 18.0);
	f_eta_spectrum[2]->SetParameters(1372.01, 0.0103744, 0.0225786, 0.71981, 9.07938);
	f_eta_spectrum[2]->SetLineColor(plotColors[2]);

	f_eta_spectrum[3] = new TF1("f_eta_spectrum_0", "2*TMath::Pi()*x*([p0]/TMath::Power((TMath::Exp(-[p1]*x-[p2]*x*x)+(x/[p3])),[p4]))", 0.0, 18.0);
	f_eta_spectrum[3]->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);
	f_eta_spectrum[3]->SetLineColor(plotColors[3]);

	f_eta_spectrum[4] = new TF1("f_eta_spectrum_0", "2*TMath::Pi()*x*(([p0]*(([p1]-1)*([p1]-1))/(([p1]*[p2]+0.53*([p1]-1))*([p1]*[p2]+0.53))*pow(([p1]*[p2]+TMath::Sqrt(0.53*0.53+x*x))/([p1]*[p2]+0.53),-1*[p1])))", 0.0, 18.0);
	f_eta_spectrum[4]->SetParameters(0.786558, 9.24328, 0.100848);
	f_eta_spectrum[4]->SetLineColor(plotColors[4]);

	plotPhoton();
	plotPizero();
	plotEta();
}