#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TMath.h>

using namespace std;

TTree* ReturnTree(const char* inputFileName,string treeName){

    TFile* inFile = new TFile(inputFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }
    TTree* tree = (TTree*)inFile->Get(treeName.c_str());
    return tree;
}

struct HistVar1D{
	const char* canvasName = "CanvasName";
	const char* histTitle = "Hist Title";
    const char* outFileName = "Plots/Unknown.png";
	const char* xLabel = "xLabel";
    int nbin = 50;
    int xmin = 0;
    int xmax = TMath::Pi()+0.0001;
    bool isLog = false;
};

void DrawHistVectorNoTree(vector<float>& dataVec,
		HistVar1D HistParam){

    TCanvas* canvas = new TCanvas(HistParam.canvasName, HistParam.canvasName, 800, 600);
    TH1F h1("h1", HistParam.histTitle, HistParam.nbin, HistParam.xmin, HistParam.xmax);
    if (HistParam.isLog == true){    
        canvas->SetLogy();
    }

    for (float value : dataVec){
	h1.Fill(value);
    }
   

    h1.SetXTitle(HistParam.xLabel); // Set the x-axis label
    h1.Draw("HIST");
    canvas->SaveAs(HistParam.outFileName);

    delete canvas;
}
void DrawHistFVector(const char* inFileName,
		     string inTreeName,
		const char* branchName, //must be a float_t array
		const char* indexName = "nref",
		const char* canvasName = "canvasName",
                const char* histTitle = "default title",
                const char* outFileName = "default.png",
                int nbin = 50,
                int xmin = 0,
                int xmax = 200,
		int pfidset = 0,
		bool isLog = false){

    TFile* inFile_p = new TFile(inFileName,"READ");
    if (!inFile_p || inFile_p->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
        return;
}

    TTree* tree = (TTree*)inFile_p->Get(inTreeName.c_str());

    int n;

    vector<float>* branch = nullptr;
    vector<int>* pfid = nullptr;

    tree->SetBranchAddress(indexName,&n);
    tree->SetBranchAddress(branchName,&branch);
    tree->SetBranchAddress("pfId",&pfid);
    //    tree->Scan("pfEta:pfPhi");

    cout << "In Function Branch Name:" << branchName << endl;

    TCanvas* canvas = new TCanvas(canvasName, canvasName, 800, 600);
    TH1F h1("h1", histTitle, nbin, xmin, xmax);

    if (isLog == true){    
        canvas->SetLogy();
    }
   
    for (int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        for (int j = 0; j < n; j++){
	    if (!branch->empty()){
		if (pfid->at(j) == pfidset){
		    //cout << "pfid: " << pfid << endl;
		    //cout << branch->at(j) << endl;
                    h1.Fill(branch->at(j));
		}
            }
        }
    }

    h1.SetXTitle(branchName); // Set the x-axis label
    h1.Draw("HIST");
    canvas->SaveAs(outFileName);

    inFile_p->Close();
    delete inFile_p;
}

void DrawHistFArray(TTree* tree,
		const char* branchName, //must be a float_t array
		const char* indexName = "nref",
		const char* canvasName = "canvasName",
                const char* histTitle = "default title",
                const char* outFileName = "default.png",
                int nbin = 50,
                int xmin = 0,
                int xmax = 200,
		bool isLog = false){
    int n;
    Float_t branch[100];

    tree->SetBranchAddress(indexName,&n);
    tree->SetBranchAddress(branchName,branch);

    if (!n||!branch[0]){
	    cerr << "BranchAintFound" << endl;
	    cerr << branchName << endl;
	    cerr << indexName << endl;
	    return;
    }

    TCanvas* canvas = new TCanvas(canvasName, canvasName, 800, 600);
    TH1F h1("h1", histTitle, nbin, xmin, xmax);

    if (isLog == true){    
        canvas->SetLogy();
    }
   

    for (int i = 0; i < tree->GetEntries(); i++){
        tree->GetEntry(i);
        if (tree->GetEntry(i) <= 0) {  // Check if GetEntry fails
            cerr << "Error: Failed to read entry " << i << endl;
            break;
        }
        for (int j = 0; j < n; j++){
            if (!branch[j]){
                cerr << "Error here in Drawhist: " << endl;
                break;
            }
            else{
                    h1.Fill(branch[j]);
		}
        }
    }

    h1.SetXTitle(branchName); // Set the x-axis label
    h1.Draw("HIST");
    canvas->SaveAs(outFileName);

    delete canvas;
}

tuple<vector<float>,vector<float>,vector<float>>GetPt1Pt2(TTree* JetAnalyserTree){

    float jtpt1 = -1;
    float jtpt2 = -1;
    float A_J = 0;
    Float_t jtpt[100];
    int nref;
    
    vector<float> jtpt1vector;
    vector<float> jtpt2vector;
    vector<float> ajvector;

    JetAnalyserTree->SetBranchAddress("nref",&nref);
    JetAnalyserTree->SetBranchAddress("jtpt",jtpt);

    for (int i = 0; i < JetAnalyserTree->GetEntries(); i++){
        JetAnalyserTree->GetEntry(i);
	vector<float> orderedJtpt = {};
        for (int j = 0; j < nref; j++){
            if (!jtpt[j]){
		cerr << "Error in entry number: "<< i <<j << endl;
                break;
            }
            else{
                orderedJtpt.push_back(jtpt[j]);
	    }
	}

	sort(orderedJtpt.begin(),orderedJtpt.end(),greater<float>());
	if (nref >= 2){
		jtpt1 = orderedJtpt[0];
		jtpt2 = orderedJtpt[1];
		A_J = (jtpt1-jtpt2)/(jtpt1+jtpt2);
		jtpt1vector.push_back(jtpt1);
		jtpt2vector.push_back(jtpt2);
		ajvector.push_back(A_J);
	    }
    }
    return make_tuple(jtpt1vector,jtpt2vector, ajvector);
}

void GetPt1Pt2_New(TTree* Tree){

    float jtpt1 = -1;
    float jtpt2 = -1;
    float A_J = 0;
    Float_t jtpt[50];
    int nref;

    Tree->SetBranchAddress("nref",&nref);
    Tree->SetBranchAddress("jtpt",jtpt);

    for (int i = 0; i < Tree->GetEntries(); i++){
        Tree->GetEntry(i);
	if (jtpt[0] < jtpt[1]){
	    cout << "i: " << i << endl;
            for (int j = 0; j < nref; j++){
		cout << j << ":" << jtpt[j] << endl;
	    }
	}  
    }
}

//Fixlater
vector<float> GetEmiss(TTree* tree){
    float EMissjpt = 0;
    float jptxtot, jptytot;

    int nref;
    Float_t jtpt[100];
    Float_t jtphi[100];

    vector<float> EMiss;

    tree->SetBranchAddress("nref",&nref);
    tree->SetBranchAddress("jtpt",jtpt);
    tree->SetBranchAddress("jtphi",jtphi);

    for (int i = 0; i < tree->GetEntries(); i++){
	EMissjpt = 0;
	jptxtot = 0;
	jptytot = 0;
        tree->GetEntry(i);
        for (int j = 0; j < nref; j++){
            if (!jtpt[j]){
		cerr << "Error in entry number: "<< i <<j << endl;
                break;
            }
	    else{
		jptxtot += TMath::Cos(jtphi[j])*jtpt[j];
		jptytot += TMath::Sin(jtphi[j])*jtpt[j];
	    }
        //SUbleading jet with pt>50
	}
	EMissjpt = TMath::Sqrt(jptxtot*jptxtot + jptytot*jptytot);
	EMiss.push_back(EMissjpt);
    }
    return EMiss;
}

void draw2D(TTree* tree){
	int n;
	int selectionE = 2000;
	vector<float>* pfpt = nullptr;
	vector<float>* pfE = nullptr;
	vector<float>* pfEta = nullptr;
	vector<float>* pfPhi = nullptr;

	vector<int> event_array = {};

	tree->SetBranchAddress("nPF",&n);
	tree->SetBranchAddress("pfPt",&pfpt);
	tree->SetBranchAddress("pfE",&pfE);
	tree->SetBranchAddress("pfPhi",&pfPhi);
	tree->SetBranchAddress("pfEta",&pfEta);

        auto c1 = new TCanvas("c","canvas",0,0,600,400);
        TH2D* graph = new TH2D ("h2","Distribution of 1 Event;Phi;Eta", 50, -6, 6, 50,-4,4);


	for(int i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		for(int j = 0;j < n;j++){
			if (pfE->at(j) > 2000){
				cout << "event: " << i << " part index: " << j << endl;
				cout << "pfpt: " << pfpt->at(j) << endl;
				cout << "pfE: " << pfE->at(j) << endl;
				event_array.push_back(i);
				cout << endl;
			}
		}
	}
//Edit this later
	tree->GetEntry(1);
	for(int i = 0; i < n; i++){
		graph->Fill(pfPhi->at(i),pfEta->at(i),pfpt->at(i));
	}
	graph->Draw("LEGO");
	c1->SetTheta(45.);
    	c1->SetPhi(60.);
	c1->SetLogz();
	c1->SaveAs("Plots_partialskim_Event5139_ptE2801_Logz.jpg");

}
struct CutsVector{
    string name;
	vector<float> cutVector;
	float cutmin;
	float cutmax;
};

struct CutsBranches{
    string treename;
    string name;
    float cutmin;
    float cutmax;
};

//with dphi and eta cuts:
//In the future make vectors and trees all together
//automate this so I can feed in arbitrary numbers of vectors and tree branches and it would change the cut and caption accordingly
void PlotBranchWithCuts(const char* inFileName,
			const char* branchName,
			HistVar1D h,
            CutsVector dPhiCut,
            CutsBranches etaCut){

    int nref, hiBin;
    float dPhi;
    Float_t branch[100];
    Float_t jteta[100];

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* JetAnalyzerTree = (TTree*)inFile->Get("akCs4FlowPFJetAnalyzer/t");
    TTree* HiTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");

    TCanvas* canvas = new TCanvas(h.canvasName, h.canvasName, 900, 900);
    TH1F h1("h1", h.histTitle, h.nbin, h.xmin, h.xmax);


    if (h.isLog == true){
        canvas->SetLogy();
    }

    JetAnalyzerTree->SetBranchAddress(etaCut.name.c_str(),jteta);
    JetAnalyzerTree->SetBranchAddress(branchName,branch);
    JetAnalyzerTree->SetBranchAddress("nref",&nref);

    for (int i = 0; i < JetAnalyzerTree->GetEntries(); i++){
	    JetAnalyzerTree->GetEntry(i);
	    HiTree->GetEntry(i);
        dPhi = dPhiCut.cutVector[i];
        if (dPhi > dPhiCut.cutmin && dPhi < dPhiCut.cutmax){
            for (int j = 0; j < nref; j++){
                if (jteta[j] > etaCut.cutmin && jteta[j] < etaCut.cutmax){
                    cout << "branch: " << branch[j] << endl;
                    cout << "eta: " << jteta[j] << endl;
                    cout << "dPhi: " << dPhi << endl;
                    h1.Fill(branch[j]);
                }
            }
        }
    }

    gStyle->SetOptStat(0);
    gPad->SetTicks();
    h1.SetXTitle(h.xLabel); // Set the x-axis label
    h1.SetYTitle("# of events");
    h1.Draw("HIST E1");

    TLatex cuts1;
    TLatex cuts2;

    cuts1.SetNDC(); // use normalized device coordinates (0–1)
    cuts1.SetTextSize(0.03);
    cuts2.SetNDC();
    cuts2.SetTextSize(0.03);

    if (dPhiCut.cutmin > 0 || dPhiCut.cutmax < TMath::Pi()){
        cuts1.DrawLatex(0.3, 0.75, Form("%.2f < %s < %.2f",dPhiCut.cutmin,dPhiCut.name.c_str(),dPhiCut.cutmax));
    }
    if (etaCut.cutmin > 0 || etaCut.cutmax < TMath::Pi()){
    	cuts2.DrawLatex(0.3,0.70, Form("%.2f < |%s| < %.2f",etaCut.cutmin,etaCut.name.c_str(),etaCut.cutmax));
    }
    canvas->SaveAs(h.outFileName);
}


void PlotAJ(const char* inFileName,
	     vector<float> dphiVector,
	     const char* canvasName = "canvasName",
         const char* histTitle = "Jet A_J (cuts dPHi > X)",
         const char* outFileName = "Plots/PartialSkim_dPhi.png",
	     int centmin = 0,
	     int centmax = 200,
	     float dphimin = 0,
	     float dphimax = 7,
         int nbin = 50,
         int xmin = 0,
         int xmax = TMath::Pi()+0.0001){
    float dphi, A_J;
    int HiBin;
    Float_t jtphi[100];
    Float_t jtpt[100];

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* JetAnalyzerTree = (TTree*)inFile->Get("akCs4FlowPFJetAnalyzer/t");
    TTree* HiTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");

    TCanvas* canvas = new TCanvas(canvasName, canvasName, 900, 900);
    TH1F h1("h1", histTitle, nbin, xmin, xmax);

    TLatex dphicuts;
    dphicuts.SetNDC(); // use normalized device coordinates (0–1)
    dphicuts.SetTextSize(0.04);

    TLatex latex;
    latex.SetNDC(); // use normalized device coordinates (0–1)
    latex.SetTextSize(0.04);

    JetAnalyzerTree->SetBranchAddress("jtpt",jtpt);
    HiTree->SetBranchAddress("hiBin",&HiBin);

    int eventcounter = 0;
    for (int i = 0; i < dphiVector.size(); i++){
	dphi = dphiVector[i];
	HiTree->GetEntry(i);
	JetAnalyzerTree->GetEntry(i);
	A_J = (jtpt[0]-jtpt[1])/(jtpt[0]+jtpt[1]);
	if (dphi > dphimin && dphi < dphimax){
	    if (HiBin > centmin && HiBin < centmax){
		eventcounter += 1;
		h1.Fill(A_J);
	    }
	}
    }  
    cout << dphimin << "< dphi <" << dphimax << endl;
    cout << centmin << "< cent <" << centmax << endl;
    cout << eventcounter << " Events with the cut" << endl;

    gStyle->SetOptStat(0);
    gPad->SetTicks();
    h1.SetXTitle("A_{J}"); // Set the x-axis label
    h1.SetYTitle("# of events");
    h1.Draw("HIST E1");
    

    latex.DrawLatex(0.3, 0.85, Form("cent: %d-%d%%",TMath::Nint(centmin/2),TMath::Nint(centmax/2)));
    if (dphimin > 0 || dphimax < TMath::Pi()){
        dphicuts.DrawLatex(0.3, 0.75, Form("dphi cuts: %.2f-%.2f",dphimin,dphimax));
    }

    canvas->SaveAs(outFileName);
}
//under construction
void PlotDPhi(const char* inFileName,
  	     string inTreeName,
	     const char* branchName = "jtphi",
	     const char* canvasName = "canvasName",
         const char* histTitle = "dPhi",
         const char* outFileName = "Plots/PartialSkim_dPhi.png",
	     int centmin = 0,
	     int centmax = 200,
         int nbin = 50,
         int xmin = 0,
         int xmax = TMath::Pi()+0.0001){

    float phi1;
    float phi2;
    float dphi;
    float pi = TMath::Pi();
    int HiBin;
    Float_t jtphi[100];
    vector<float> jtdphi = {};

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* Tree = (TTree*)inFile->Get(inTreeName.c_str());
    TTree* HiTree = (TTree*)inFile->Get("hiEvtAnalyzer/HiTree");

    TCanvas* canvas = new TCanvas(canvasName, canvasName, 800, 600);
    TH1F h1("h1", histTitle, nbin, xmin, xmax);

    TLatex latex;
    latex.SetNDC(); // use normalized device coordinates (0–1)
    latex.SetTextSize(0.04);
   
    Tree->SetBranchAddress(branchName,jtphi);
    HiTree->SetBranchAddress("hiBin",&HiBin);

    for (int i = 0; i < Tree->GetEntries(); i++){
        Tree->GetEntry(i);
	HiTree->GetEntry(i);

	phi1 = jtphi[0];
	phi2 = jtphi[1];
	dphi = phi1 - phi2;
	if (pi <= dphi && dphi <= 2 * pi){
	    dphi = (dphi - 2*pi);
	}
	else if (-2*pi <= dphi && dphi <= -pi){
	    dphi = dphi + 2*pi;
	}
	dphi = TMath::Abs(dphi);

	//cout << "Centrality: " << HiBin << " " << endl;
	if(HiBin > centmin && HiBin < centmax){
	    h1.Fill(dphi);
	    jtdphi.push_back(dphi);
	}
    }  

    h1.SetXTitle("dPhi"); // Set the x-axis label
    h1.Draw("HIST");

    latex.DrawLatex(0.5, 0.85, Form("cent: %d-%d%%",TMath::Nint(centmin/2),TMath::Nint(centmax/2)));
    canvas->SaveAs(outFileName);
}

vector<float> GetDPhi(const char* inFileName,
  	     string inTreeName,
	     const char* branchName = "jtphi"){

    float phi1;
    float phi2;
    float dphi;
    float pi = TMath::Pi();
    Float_t jtphi[100];
    vector<float> jtdphi = {};

    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* Tree = (TTree*)inFile->Get(inTreeName.c_str());
    Tree->SetBranchAddress(branchName,jtphi);

    for (int i = 0; i < Tree->GetEntries(); i++){
        Tree->GetEntry(i);
	phi1 = jtphi[0];
	phi2 = jtphi[1];
	dphi = phi1 - phi2;

	if (pi <= dphi && dphi <= 2 * pi){
	    dphi = (dphi - 2*pi);
	    //cout << "New dphi:" << dphi << endl;
	}
	else if (-2*pi <= dphi && dphi <= -pi){
	    dphi = dphi + 2*pi;
	}
	dphi = TMath::Abs(dphi);
	//cout << dphi << endl;
	jtdphi.push_back(dphi);
    }
    return jtdphi;
}

/*
void printWithDphiCuts(){
//Missing pt, phi, eta, etc etc etc
    vector<float> Dphi = GetDPhi(inputFileName,JetAnalyserTreeString);

    vector<int> events_with_smalldphi = {};
    vector<int> events_with_bigdphi = {};

    for (int i = 0; i < Dphi.size(); i++) {
        //cout << Dphi[i] << endl;
	if (Dphi[i] < 0.8){
		events_with_smalldphi.push_back(i);
	}
	if (Dphi[i] > 3.14){
		events_with_bigdphi.push_back(i);
	}
    }

    cout << "———————— Events with Dphi < 0.8 –––––––––" << endl;
    for (int i = 0; i < 10; i++) {
        int eventIndex = events_with_smalldphi[i];
        JetAnalyserTree->GetEntry(eventIndex);
	HiTree->GetEntry(eventIndex);

        cout << "Event num: " << eventIndex << endl;
        cout << "jet index | pt       phi      eta |" << endl;

        for (int j = 0; j < nref; j++) {
            cout << j << "         | "
                 << jtpt[j] << "   "
                 << jtphi[j] << "   "
                 << jteta[j] << "   "<<endl;
	}
	
	cout << "HiBin: " << HiBin << endl;
        cout << endl;
    }

   
    cout << "———————— Events with Dphi > 3.14 –––––––––" << endl;
    for (int i = 0; i < 10; i++) {
        int eventIndex = events_with_bigdphi[i];
        JetAnalyserTree->GetEntry(eventIndex);
	HiTree->GetEntry(eventIndex);

        cout << "Event num: " << eventIndex << endl;
        cout << "jet index | pt       phi      eta|" << endl;

        for (int j = 0; j < nref; j++) {
            cout << j << "         | "
                 << jtpt[j] << "   "
                 << jtphi[j] << "   "
                 << jteta[j] << endl;
	}

	cout << "HiBin: " << HiBin << endl;
        cout << endl;
	
    }
}
*/

void run(const char* inputFileName) {
    string JetAnalyserTreeString = "akCs4FlowPFJetAnalyzer/t";
    string PFTreeString = "particleFlowAnalyser/pftree";
    string HiTreeString = "hiEvtAnalyzer/HiTree";

    TTree* JetAnalyserTree = ReturnTree(inputFileName,JetAnalyserTreeString);
    TTree* PFTree = ReturnTree(inputFileName, PFTreeString);
    TTree* HiTree = ReturnTree(inputFileName, HiTreeString);

    vector<float> Dphi = GetDPhi(inputFileName,JetAnalyserTreeString);

    int nref;
    int HiBin;

    Float_t jtpt[100];
    Float_t jtphi[100];
    Float_t jteta[100];
    Float_t jtPfCHF[100];

    JetAnalyserTree->SetBranchAddress("nref",&nref);
    JetAnalyserTree->SetBranchAddress("jtpt",jtpt);
    JetAnalyserTree->SetBranchAddress("jtphi",jtphi);
    JetAnalyserTree->SetBranchAddress("jteta",jteta);

    HiTree->SetBranchAddress("hiBin",&HiBin);

    HistVar1D DefaultHistParams{"StructCvsName",
        "Struct Title",
        "Plots/TestStruct.png",
        "xLabel",
        50,0,1,true};
    
    CutsVector dPhiCut_low{"dPhi",Dphi, 0, 0.8};
    CutsVector dPhiCut_high{"dPhi",Dphi, TMath::Pi()*(7.0/8.0), TMath::Pi()};
    CutsBranches jtEtaCut{"JetAnalyzer", "jteta", 0, 2};        
    
   /* for (int i = 0; i < 2; i++){
	    for (int j = 0;j < branchnamearr.size() ;j++){
		    string branchname = branchnamearr[j];
		    string title = "";
		    string filename = "";

		    title += branchname;
		    title += " of jets ";
		    title += dphititle;

		    filename += "Plots/PartialSkim_";
		    filename += branchname;
		    filename += "_";
		    filename += dphicut;
		    filename += "_log";
		    filename += ".png";		    

		    cout << "BranchName: " << branchname << endl;
		    cout << "Title: " << title << endl;
		    cout << filename << endl;
		    cout << dphimin << " < dphi < " << dphimax << endl;
		    cout << endl;
            PlotBranchWithCuts(inputFileName, "jtPfCHF",
                DefaultHistParams, dPhiCut_low, jtEtaCut);        
        }
    }  */  
}
