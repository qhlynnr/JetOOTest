#pragma once
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>
#include <vector>
#include <TMath.h>

inline TTree* ReturnTree(std::string inputFileName,
    std::string treeName){

    TFile* inFile = new TFile(inputFileName.c_str(),"READ");
    if (!inFile || inFile->IsZombie()) {
        cout << "Error: Could not open the file!" << endl;
    }

    TTree* tree = (TTree*)inFile->Get(treeName.c_str());
    return tree;
}