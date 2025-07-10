#pragma once
#include <vector>
#include <TMath.h>

struct HistVar1D{
	std::string histTitle = "Hist Title";
	std::string xLabel = "xLabel";
    std::string yLabel = "yLabel";
    std::string outFolderName = "Plots/";
    std::string outFileName = "Unknown.png";
    int nbin = 50;
    float xmin = 0;
    float xmax = TMath::Pi()+0.0001;
    float ymin = 0;
    float ymax = 100;
};

