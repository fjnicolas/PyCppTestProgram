////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesTrackFinder.h
//
// \brief Definition of TPCLinesTrackFinder
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_TRACKFINDER_H
#define TPC_LINES_TRACKFINDER_H


#include "TPCLinesParameters.cpp"
#include "TPCSimpleHits.h"
#include "TPCSimpleLines.h"
#include "TPCSimpleClusters.h"
#include "TPCLinesDisplay.cpp"


class TPCLinesTrackFinder {
    private:
        TrackFinderAlgorithmPsetType fTPCLinesTrackFinderPset;
        double GetHitLinearDensity(std::vector<SHit> data);
        LineEquation PerformPCA2D(std::vector<SHit>& data);
        void GetHitsInTube(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD);
        void GetHitsInTubeSingleWire(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD);
        void Get2DClusters(std::vector<SHit> & hits, double epsilon, int minPts);

        // Display app directory
        TPCLinesDisplay fDisplay;


    public:
        // constructor
        TPCLinesTrackFinder(TrackFinderAlgorithmPsetType tpcLinesTrackFinderPset);

        // main function
        std::vector<SCluster> ReconstructTracksFromHoughDirection(std::vector<SHit> &hitList, LineEquation houghLine, int trkIter);
        
};

#endif // TPC_SIMPLE_LINES_H