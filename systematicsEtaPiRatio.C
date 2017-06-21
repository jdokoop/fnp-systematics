#include <iostream>

using namespace std;

//-------------------------------------------
// Variables
//-------------------------------------------

const int NPT = 12;

float x[NPT] = {2.75,
                3.25,
                3.75,
                4.25,
                4.75,
                5.25,
                5.75,
                6.50,
                7.50,
                8.50,
                9.50,
                11.0
               };

float etapiratio[NPT] = {0.439555,
                         0.42086 ,
                         0.446485,
                         0.473496,
                         0.468078,
                         0.51012 ,
                         0.561144,
                         0.540097,
                         0.596452,
                         0.426385,
                         0.588294,
                         0.419133
                        };

float err_y_stat[NPT] = {0.0756827,
                         0.00899723,
                         0.0115065,
                         0.0155712,
                         0.0218055,
                         0.0323237,
                         0.043522,
                         0.0468465,
                         0.0779586,
                         0.0951103,
                         0.161896,
                         0.16185
                        };

float err_y_syst[NPT] = {0.107157,
                         0.040213,
                         0.0298148,
                         0.0312415,
                         0.0285021,
                         0.0307196,
                         0.0336704,
                         0.0325939,
                         0.0367457,
                         0.0270475,
                         0.0385867,
                         0.0290273
                        };

float err_x[NPT] = {0.0};

//////////////////////////////////////////////////////

float pythia_x[24] = {0.13161101632951522,
                      0.3234218864245675,
                      0.4962222763831343,
                      0.7250792103339021,
                      0.9900073117231292,
                      1.2354374847672436,
                      1.5732390933463318,
                      1.9844016573239094,
                      2.4128686327077746,
                      3.026322203265903,
                      3.6388008774067764,
                      4.158176943699732,
                      4.713867901535462,
                      5.2698025834755065,
                      5.9553984889105545,
                      6.492322690714112,
                      7.0480136485498415,
                      7.677796734097002,
                      8.36339263953205,
                      8.863514501584206,
                      9.77114306604923,
                      10.437972215452106,
                      11.271508652205702,
                      11.864245673897146
                     };

float pythia_y[24] = {0.02130148671703,
                      0.09285888374360,
                      0.15910309529612,
                      0.23075798196441,
                      0.29198147696807,
                      0.34262734584450,
                      0.39088471849865,
                      0.43144040945649,
                      0.45888374360224,
                      0.48418230563002,
                      0.49895198635145,
                      0.50821350231537,
                      0.50967584694126,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098,
                      0.5098
                     };

TGraphErrors *g_etapi_ratio;
TF1 *f_etapi_ratio_const_fit;
TF1 *f_etapi_ratio_const_published;
TF1 *f_etapi_ratio_const_published_var1;
TF1 *f_etapi_ratio_const_published_var2;

TGraph *g_pythia;

//-------------------------------------------
// Functions
//-------------------------------------------

void makePublishedRatio()
{
	g_etapi_ratio = new TGraphErrors(NPT, x, etapiratio, err_x, err_y_stat);
	g_etapi_ratio->SetMarkerColor(kBlack);
	g_etapi_ratio->SetMarkerColor(kBlack);
	g_etapi_ratio->SetMarkerStyle(20);
	g_etapi_ratio->SetMarkerSize(0.7);

	g_pythia = new TGraph(24, pythia_x, pythia_y);
	g_pythia->SetMarkerColor(kRed);
	g_pythia->SetMarkerColor(kRed);
	g_pythia->SetMarkerStyle(20);
	g_pythia->SetMarkerSize(0.7);
	g_pythia->SetLineStyle(7);
	g_pythia->SetLineColor(kRed);

	///////////////////////////////////////

	f_etapi_ratio_const_published = new TF1("f_etapi_ratio_const_published", "0.48", 0, 18);
	f_etapi_ratio_const_published_var1 = new TF1("f_etapi_ratio_const_published_var1", "0.52", 0, 18);
	f_etapi_ratio_const_published_var2 = new TF1("f_etapi_ratio_const_published_var2", "0.44", 0, 18);

	f_etapi_ratio_const_published->SetLineWidth(2);
	f_etapi_ratio_const_published_var1->SetLineWidth(2);
	f_etapi_ratio_const_published_var2->SetLineWidth(2);

	f_etapi_ratio_const_published->SetLineColor(kBlack);
	f_etapi_ratio_const_published_var1->SetLineColor(kBlue);
	f_etapi_ratio_const_published_var2->SetLineColor(kGreen+2);

	f_etapi_ratio_const_published->SetLineStyle(7);
}


void fitConstant()
{
	f_etapi_ratio_const_fit = new TF1("f_etapi_ratio_const_fit", "[0]", 0, 18);
	g_etapi_ratio->Fit(f_etapi_ratio_const_fit, "0R");
}


void plotEtaPiRatio()
{
	TCanvas * c = new TCanvas("c", "c", 600, 600);

	TH1F *hTemplate = new TH1F("hTemplate", "hTemplate", 100, 0, 18);
	hTemplate->SetTitle("");
	hTemplate->GetXaxis()->SetTitleFont(62);
	hTemplate->GetXaxis()->SetLabelFont(62);
	hTemplate->GetXaxis()->SetRangeUser(0, 20);
	hTemplate->GetYaxis()->SetTitleFont(62);
	hTemplate->GetYaxis()->SetLabelFont(62);
	hTemplate->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	hTemplate->GetYaxis()->SetTitle("#eta / #pi^{0} Ratio");
	hTemplate->GetYaxis()->SetTitleOffset(1.3);
	hTemplate->GetXaxis()->SetTitleOffset(1.3);
	hTemplate->GetYaxis()->SetRangeUser(0, 1);
	hTemplate->Draw();

	g_etapi_ratio->Draw("P,same");
	g_pythia->Draw("C,same");
	//f_etapi_ratio_const_fit->Draw("same");
	f_etapi_ratio_const_published->Draw("same");
	f_etapi_ratio_const_published_var1->Draw("same");
	f_etapi_ratio_const_published_var2->Draw("same");

	TLegend *legend = new TLegend(0.45, 0.45, 0.88, 0.65);
	legend->AddEntry(g_pythia, "PYTHIA", "l");
	legend->AddEntry(f_etapi_ratio_const_published, "Variation 1: Published fit [PPG065]", "l");
	legend->AddEntry(f_etapi_ratio_const_published_var1, "Variation 2: + Uncertainty on published fit", "l");
	legend->AddEntry(f_etapi_ratio_const_published_var2, "Variation 3: - Uncertainty on published fit", "l");
	legend->SetFillStyle(0.0);
	legend->SetLineColor(kWhite);
	legend->Draw("same");

}


void systematicsEtaPiRatio()
{
	gStyle->SetOptStat(0);

	makePublishedRatio();
	fitConstant();
	plotEtaPiRatio();
}