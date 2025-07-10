#pragma once
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>
#include <TMath.h>

inline void GetPt1Pt2(TTree* Tree){

    if (!Tree) {
        std::cerr << "Error: Tree pointer is null." << std::endl;
        return;
    }

    float jtpt1 = -1;
    float jtpt2 = -1;
    float A_J = 0;

    Float_t jtpt[50];
    int nref = 0;

    Tree->SetBranchAddress("nref",&nref);
    Tree->SetBranchAddress("jtpt",jtpt);

    Long64_t nEntries = Tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++){
        Tree->GetEntry(i);
        if (nref >= 2){
            cout << "Nref: " << nref << std::endl;
            for (int j = 0; j < nref; j++){
                std::cout << "trkpt " << j << ":" << jtpt[j] << std::endl;
            }
        }  
    }
}