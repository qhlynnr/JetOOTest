
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "HistVar1D.h"
#include "GetVar.h"

using namespace std;

vector<float> GetDPhi(const char* inFileName, string inTreeName, const char* branchName = "jtphi"){

    float phi1;
    float phi2;
    float dphi;
    float pi = TMath::Pi();
    int nref;
    Float_t jtphi[100];
    vector<float> jtdphi = {};

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* Tree = (TTree*)inFile->Get(inTreeName.c_str());
    Tree->SetBranchAddress(branchName,jtphi);
    Tree->SetBranchAddress("nref",&nref);

    for (int i = 0; i < Tree->GetEntries(); i++){
        Tree->GetEntry(i);
        phi1 = jtphi[0];
        phi2 = jtphi[1];
        dphi = phi1 - phi2;
        if (nref <= 2) {
            jtdphi.push_back(-999);  // or continue
            continue;
        }
        if (pi <= dphi && dphi <= 2 * pi){
            dphi = (dphi - 2*pi);
        }
        else if (-2*pi <= dphi && dphi <= -pi){
            dphi = dphi + 2*pi;
        }
        dphi = TMath::Abs(dphi);
        jtdphi.push_back(dphi);
    }
    inFile->Close();
    delete inFile;
    return jtdphi;
}

void DrawSelectionLatex(double x = 0.6, double y = 0.2, double dy = 0.05) {
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetNDC();
    latex.DrawLatex(x, y + 7*dy, "HighPurity = 1");
    latex.DrawLatex(x, y + 6*dy, "PVFilter = 1");
    latex.DrawLatex(x, y + 5*dy, "CCFilter = 1");
    latex.DrawLatex(x, y + 4*dy, "nVtx > 0");
    latex.DrawLatex(x, y + 3*dy, "|V_{Z}| < 15 cm");
    latex.DrawLatex(x, y + 2*dy, "jtpt1 #geq 50");
    latex.DrawLatex(x, y + 1*dy, "jtpt2 #geq 30");
    latex.DrawLatex(x, y + 0*dy, "|jt#eta| < 2");
    latex.DrawLatex(x, y - 1*dy, "#Delta#phi > 5#pi/6");
}

int run() {
    string ForestFolder = "/eos/cms/store/group/phys_heavyions/jdlang/Run3_OO_2025Data_QuickForest/";
    vector<string> ForestSubfolder = {
        "OO_394153_PhysicsIonPhysics0_Prompt_v3/crab_OO_394153_PhysicsIonPhysics0_Prompt_v3/250706_213337/0000/",
        "OO_394153_PhysicsIonPhysics2_Prompt_v3/crab_OO_394153_PhysicsIonPhysics2_Prompt_v3/250706_213442/0000/",
        "OO_394153_PhysicsIonPhysics4_Prompt_v3/crab_OO_394153_PhysicsIonPhysics4_Prompt_v3/250706_213505/0000/",
        "OO_394154_PhysicsIonPhysics0/crab_OO_394154_PhysicsIonPhysics0/250705_130926/0000/",
        "OO_394154_PhysicsIonPhysics5/crab_OO_394154_PhysicsIonPhysics5/250706_211713/0000/",
        "OO_394169_PhysicsIonPhysics0/crab_OO_394169_PhysicsIonPhysics0/250706_210830/0000/",
        "OO_394169_PhysicsIonPhysics2/crab_OO_394169_PhysicsIonPhysics2/250706_210847/0000/",
        "OO_394169_PhysicsIonPhysics4/crab_OO_394169_PhysicsIonPhysics4/250706_232638/0000/"};

    string JetAnalyserTreeString = "ak0PFJetAnalyzer/t";
    string HFAdcanalyserTreeString = "HFAdcana/adc";
    string SkimTreeString = "skimanalysis/HltTree";

    float jtpt1 = -1;
    float jtpt2 = -1;
    float A_J = 0;
    vector<float>* zVtx = nullptr;
    vector<float>* ptSumVtx = nullptr;
    vector<bool>* highPurity = nullptr;
    vector<float>* trackEta = nullptr;

    Float_t jtpt[50];
    Float_t jteta[50];
    int nref = 0;
    int BestVertex = -1;
    int PVFilter, CCFilter,nVtx,nTrk;
    int mMaxL1HFAdcPlus, mMaxL1HFAdcMinus;
    int jetcount = 0;
    int jetcountTotal = 0;
    int testcontrolledjetcut = 1;
    int numberOfEventsWithAJ = 0;
    HistVar1D hvar = {"HiForestMiniAOD Jet A_{J} Distribution",
        "A_{J}",
        "Number of Jets",
        "Plots/",
        "HiForestMiniAOD_AJ_lynntest.png",
        30, 0, 1, 0, 100};

    TCanvas* canvas = new TCanvas("c", "canvas", 800, 600);
    TH1F hjtpt1("hjtpt1", "jtpt1", 100, 0, 500);
    TH1F hNtrk("hNtrk", "Number of Tracks (highPurity)", 100, 0, 1500);
    TH1F hNtrkOriginal("hNtrkOriginal", "Number of Tracks (original)", 100, 0, 1500);
    TH1F h_highntrk("h_highntrk", "Jet A_{J} of OO Sample (ntrk > 400)", hvar.nbin, hvar.xmin, hvar.xmax);
    TH1F h_lowntrk("h_lowntrk", "Jet A_{J} of OO Sample (ntrk < 200)", hvar.nbin, hvar.xmin, hvar.xmax);
    TH1F h_midntrk("h_midntrk", "Jet A_{J} of OO Sample (200 < ntrk < 00)", hvar.nbin, hvar.xmin, hvar.xmax);

    hjtpt1.SetDirectory(0);
    hNtrk.SetDirectory(0);
    hNtrkOriginal.SetDirectory(0);
    h_highntrk.SetDirectory(0);
    h_lowntrk.SetDirectory(0);
    h_midntrk.SetDirectory(0);

    for(string subfolder : ForestSubfolder) {
        string inputFileFolder = ForestFolder + subfolder;
        cout << "Processing subfolder: " << subfolder << endl;
        for(int i = 1; i < 120; i++){
            bool etaCut = false;
            bool highPurityBool = true;

            string inputFileName = inputFileFolder + Form("HiForestMiniAOD_%i.root",i);
            TFile* inFile = new TFile(inputFileName.c_str(),"READ");
            if (!inFile || inFile->IsZombie()) {
                cout << "Error: Could not open the file!" << endl;
                break;
            }

            cout << "Processing file: " << i << endl;
            
            TTree* JetAnalyserTree = (TTree*)inFile->Get(JetAnalyserTreeString.c_str());
            TTree* HFAdcanaTree = (TTree*)inFile->Get(HFAdcanalyserTreeString.c_str());
            TTree* SkimTree = (TTree*)inFile->Get(SkimTreeString.c_str());
            TTree* PPTracksTree = (TTree*)inFile->Get("ppTracks/trackTree");

            SkimTree->SetBranchAddress("pprimaryVertexFilter",&PVFilter);
            SkimTree->SetBranchAddress("pclusterCompatibilityFilter",&CCFilter);
            PPTracksTree->SetBranchAddress("nVtx",&nVtx);
            PPTracksTree->SetBranchAddress("zVtx",&zVtx);
            PPTracksTree->SetBranchAddress("ptSumVtx",&ptSumVtx);
            PPTracksTree->SetBranchAddress("nTrk",&nTrk);
            PPTracksTree->SetBranchAddress("trkEta",&trackEta);
            PPTracksTree->SetBranchAddress("highPurity",&highPurity);
            HFAdcanaTree->SetBranchAddress("mMaxL1HFAdcPlus",&mMaxL1HFAdcPlus);
            HFAdcanaTree->SetBranchAddress("mMaxL1HFAdcMinus",&mMaxL1HFAdcMinus);
            JetAnalyserTree->SetBranchAddress("nref",&nref);
            JetAnalyserTree->SetBranchAddress("jtpt",jtpt);
            JetAnalyserTree->SetBranchAddress("jteta",jteta);

            Long64_t nEntries = JetAnalyserTree->GetEntries();
            vector<float> DphiVector = GetDPhi(inputFileName.c_str(),JetAnalyserTreeString);
            float dPhi = 0;
            int nTrkNew = 0;
            for (Long64_t entrynum = 0; entrynum < nEntries; entrynum++){
                JetAnalyserTree->GetEntry(entrynum);
                SkimTree->GetEntry(entrynum);
                PPTracksTree->GetEntry(entrynum);
                HFAdcanaTree->GetEntry(entrynum);

                if (nref <= 2) continue;

                jtpt1 = jtpt[0];
                jtpt2 = jtpt[1];
                A_J = (jtpt1-jtpt2)/(jtpt1+jtpt2);
                dPhi = DphiVector.at(entrynum);

                //h1.Fill(A_J);
                hjtpt1.Fill(jtpt1);

                BestVertex = -1;
                for (int vertexnum = 0; vertexnum < nVtx; vertexnum++) {
                    if (ptSumVtx->at(vertexnum) > 0 || BestVertex == -1 || ptSumVtx->at(vertexnum) > ptSumVtx->at(BestVertex)) {
                        BestVertex = vertexnum;
                    }
                }
                nTrkNew = 0;
                for(int trk = 0; trk < trackEta->size(); trk++) {
                    if (highPurity->at(trk) && fabs(trackEta->at(trk)) < 2.4) {
                        nTrkNew++;
                    }
                }

                float VZ = zVtx->at(BestVertex);
                etaCut = true;
                for (int jtnum = 0; jtnum < nref; jtnum++) {
                    if (abs(jteta[jtnum]) >= 2.4) {
                        etaCut = false;
                        break;
                    }
                }
                for (int n=0; n < nref; n++){
                    jetcountTotal++;
                    if (jtpt[n] > 45){
                        jetcount++;
                    }
                }
                hNtrk.Fill(nTrkNew);
                hNtrkOriginal.Fill(nTrk); // Fill the original number of tracks
                if (!etaCut) continue; // Apply eta cut
                if (PVFilter != 1 || CCFilter != 1 || fabs(VZ) >= 15 || nVtx <= 0) continue;
                if (mMaxL1HFAdcMinus < 14 && mMaxL1HFAdcPlus < 14) continue;
                if (jtpt1 < 50) continue;
                if (jtpt2 < 30) continue; // Apply jet pt cuts
                if (dPhi < (5.0/6) * TMath::Pi()) continue; // Apply dPhi cut
                numberOfEventsWithAJ++;
                cout << "Event num: " << entrynum 
                    << ", dPhi: " << dPhi 
                    << ", A_J: " << A_J 
                    << ", jtpt1: " << jtpt1 
                    << ", jtpt2: " << jtpt2 << endl;
                    if (nTrkNew > 400) {
                        h_highntrk.Fill(A_J); // Fill the histogram for events with nTrk > 600
                    } else if (nTrkNew < 400 && nTrkNew > 200) {
                        h_midntrk.Fill(A_J); // Fill the histogram for events with nTrk < 600 and > 200
                    } else {
                        h_lowntrk.Fill(A_J); // Fill the histogram for events with nTrk
                    }
            }
            cout << "Closing file: " << i << endl;
            inFile->Close();
            delete inFile;
            /*if (jetcount > testcontrolledjetcut){
                break;
            }*/
            cout << "Total number of jets with pt > 45 GeV After File " << i << ": " << jetcount << endl;
            cout << "Total number of jets with After File " << i << ": " << jetcountTotal << endl;
            cout << "Number of events drawn: " << numberOfEventsWithAJ << endl;

        }
        cout << "Processed subfolder: " << subfolder << endl;
        cout << "Total number of jets with pt > 45 GeV After folder" << subfolder << ": " << jetcount << endl;
        cout << "Total number of jets After folder" << subfolder << ": " << jetcountTotal << endl;
        cout << "Number of events drawn: " << numberOfEventsWithAJ << endl;

       /* if (jetcount > testcontrolledjetcut){
            break;
        }*/

    }
    cout << "Total number of jets with pt > 45 GeV FINAL: " << jetcount << endl;
    gStyle->SetOptStat(0); // Disable statistics box

    TCanvas* c_jtpt = new TCanvas("c_jtpt", "c_jtpt", 800, 600);
    hjtpt1.SetLineColor(kBlue+1);
    hjtpt1.SetLineWidth(2);
    hjtpt1.Draw("HIST");
    hjtpt1.SetXTitle("Leading jtpt (GeV/c)");

    canvas->cd();
    DrawSelectionLatex();

        TCanvas* c_ntrk = new TCanvas("c_ntrk", "c_ntrk", 800, 600);
        hNtrk.Draw("HIST");
        hNtrkOriginal.Draw("HIST SAME"); 
        hNtrk.SetLineColor(kRed);
        hNtrkOriginal.SetLineColor(kBlack);
        hNtrk.SetLineWidth(2);
        hNtrkOriginal.SetLineWidth(2);
        gStyle->SetOptStat(0); // Disable statistics box

        TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->AddEntry(&hNtrk, "highPurity", "l");
        leg->AddEntry(&hNtrkOriginal, "original", "l");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();


    TCanvas* canvas_highntrk = new TCanvas("c_highntrk", "chighntrk", 800, 800);
    h_highntrk.SetLineColor(kBlack);
    h_highntrk.SetLineWidth(2);
    h_highntrk.Draw("HIST");
    h_highntrk.SetXTitle(hvar.xLabel.c_str()); // Set the x
    h_highntrk.SetYTitle(hvar.yLabel.c_str());
    DrawSelectionLatex();

    TCanvas* canvas_midntrk = new TCanvas("c_midntrk", "chmidntrk", 800, 800);
    h_midntrk.SetLineColor(kBlack);
    h_midntrk.SetLineWidth(2);
    h_midntrk.Draw("HIST");
    h_midntrk.SetXTitle(hvar.xLabel.c_str()); // Set the x
    h_midntrk.SetYTitle(hvar.yLabel.c_str()); DrawSelectionLatex();

    TCanvas* canvas_lowntrk = new TCanvas("c_lowntrk", "clowntrk", 800, 800);
    h_lowntrk.SetLineColor(kBlack);
    h_lowntrk.SetLineWidth(2);
    h_lowntrk.Draw("HIST");
    h_lowntrk.SetXTitle(hvar.xLabel.c_str()); // Set the x
    h_lowntrk.SetYTitle(hvar.yLabel.c_str());
    DrawSelectionLatex();


   // canvas->SetLogy(); // Set log scale if specified
    canvas->SaveAs(hvar.outFileName.c_str());
    c_jtpt->SetLogy(); // Set log scale if specified
    c_jtpt->SaveAs("Plots/LeadingJetPt.png");
    c_ntrk->SetLogy(); // Set log scale if specified
    c_ntrk->SaveAs("Plots/NumberOfTracksCDF.png");
    canvas_highntrk->SaveAs("Plots/JetAJ_HighNtrk.png");
    canvas_lowntrk->SaveAs("Plots/JetAJ_LowNtrk.png");
    canvas_midntrk->SaveAs("Plots/JetAJ_MidNtrk.png");

    return 0;
}