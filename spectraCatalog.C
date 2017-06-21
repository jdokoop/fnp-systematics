//--------------------------------------------------------------------
// Catalog of fits to each published spectrum, plus all corresponding
// variations
//--------------------------------------------------------------------

//------------------------------
// Variables
//------------------------------

//Etas
TF1 *f_eta_fit;
TF1 *f_eta_var1;
TF1 *f_eta_var2;
TF1 *f_eta_var3;
TF1 *f_eta_var4;

//Pizeros
TF1 *f_pizero_fit;
TF1 *f_pizero_var1;
TF1 *f_pizero_var2;

//Photons
TF1 *f_photon_fit;
TF1 *f_photon_var1;
TF1 *f_photon_var2;
TF1 *f_photon_var3;

//------------------------------
// Functions
//------------------------------

void spectraCatalog()
{
	//////////////////////////
	// ETAS
	//////////////////////////
	f_eta_fit = new TF1("f_eta_fit", "TMath::Power(TMath::Sqrt(1 + (0.135/x)*(0.135/x)),-1.0)*0.5*TMath::Sqrt(1 + (0.135/x)*(0.135/x))*(2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_eta_fit->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);

	f_eta_var1 = new TF1("f_eta_var1", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_eta_var1->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);

	f_eta_var2 = new TF1("f_eta_var2", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_eta_var2->SetParameters(1201.68, -0.0135258, 0.0215599, 0.782018, 9.31206);

	f_eta_var3 = new TF1("f_eta_var3", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_eta_var3->SetParameters(1372.01, 0.0103744, 0.0225786, 0.71981, 9.07938);

	f_eta_var4 = new TF1("f_eta_var4", "2*TMath::Pi()*x*(([0]*(([1]-1)*([1]-1))/(([1]*[2] + 0.53*([1] - 1))*([1]*[2] + 0.53)) * pow(([1]*[2] + TMath::Sqrt(0.53*0.53 + x*x))/([1]*[2]+0.53),-1*[1])))", 0.0, 18.0);
	f_eta_var4->SetParameters(0.786558, 9.24328, 0.100848);

	//////////////////////////
	// PIZEROS
	//////////////////////////
	f_pizero_fit = new TF1("f_pizero_fit", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_pizero_fit->SetParameters(254.066, 0.470186, 0.0380066, 0.737713, 8.28442);

	f_pizero_var1 = new TF1("f_pizero_var1", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_pizero_var1->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);

	f_pizero_var2 = new TF1("f_pizero_var2", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_pizero_var2->SetParameters(545.206, 0.421156, -0.00156003, 0.58791, 8.04994);

	//////////////////////////
	// PIZEROS
	//////////////////////////
	f_photon_fit = new TF1("f_photon_fit", "2*TMath::Pi()*x*[0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])", 0, 18.0);
	f_photon_fit->SetParameters(41.0, -0.89, 0.35, 0.67, 8.15);

	f_photon_var1 = new TF1("f_photon_var1", "2*TMath::Pi()*x*([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4]))", 0.0, 18.0);
	f_photon_var1->SetParameters(13.0488, -0.191588, 0.0164036, 0.999159, 8.42105);

	f_photon_var2 = new TF1("f_photon_var2", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_photon_var2->SetParameters(2.42345e-01, -8.27585e-02, 9.18447e-03, 4.13943e+00, 1.36974e+01);

	f_photon_var3 = new TF1("f_photon_var3", "2*TMath::Pi()*x*(([0]/TMath::Power((TMath::Exp(-[1]*x-[2]*x*x)+(x/[3])),[4])))", 0.0, 18.0);
	f_photon_var3->SetParameters(2.42345e-01, -8.27585e-02, 9.18447e-03, 4.13943e+00, 1.36974e+01);
}