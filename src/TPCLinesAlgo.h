////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesAlgo
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_ALGO_H
#define TPC_LINES_ALGO_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>


#include "TPCLinesParameters.cpp"
#include "TPCSimpleHits.h"
#include "TPCLinesHough.h"
#include "TPCSimpleClusters.h"
#include "TPCLinesTrackFinder.h"
#include "TPCLinesDisplay.cpp"


class TPCLinesAlgo {
    private:
        // Algorithm parameters
        TPCLinesAlgoPsetType fTPCLinesPset;

        // Channel boundaries
        std::map<std::string, std::vector<int>> fChB =  {
            {"U0", {0, 1983}},
            {"V0", {1984, 3967}},
            {"C0", {3968, 5631}},
            {"U1", {5632, 7631}},
            {"V1", {7632, 9599}},
            {"C1", {9600, 11263}}
        };

        // Input variables
        size_t fNTotalHits;
        std::vector<SHit> fHitList;
        SVertex fVertex;

        // Hough algorithm
        TPCLinesHough fHoughAlgo;

        // Track finder and clustering algorithm
        TPCLinesTrackFinder fTrackFinder;

        // Display app directory
        TPCLinesDisplay fDisplay;
        std::string fDisplayAppPath;

    public:
        // constructor
        TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset, std::string displayPath="");
        
        // Function to set the input variables
        void SetHitList(std::string view,
                        std::vector<int>& vertex,
                        std::vector<int> *_X,
                        std::vector<double> *_Y,
                        std::vector<double> *_Int,
                        std::vector<double> *_Wi,
                        std::vector<double> *_ST,
                        std::vector<double> *_ET,
                        std::string eventLabel="");
        
        // Function to analyze the view
        std::map<std::string, double> AnaView();        

        // Display
        void Display();                            
};

#endif // TPC_SIMPLE_LINES_H