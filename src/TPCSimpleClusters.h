////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleClusters.h
//
// \brief Definition of SimpleTPCClusters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_CLUSTERS_H
#define TPC_SIMPLE_CLUSTERS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

#include "TPCSimpleHits.h"
#include "TPCSimpleLines.h"

class SCluster {
    private:
        int fNHits = 0;

        std::vector<SHit> fHitList;
        std::vector<double> fConnectednessV;
        std::vector<double> fCompactnessV;
        std::vector<double> fConnectedness1DV;
        std::vector<double> fWidths;
        
        double fCompactness;
        double fCompactnessRMS;
        double fConnectedness;
        double fConnectednessRMS;
        double fConnectedness1D;
        double fConnectedness1DRMS;
        double fAverageWidth;

        static constexpr double DefaultMaxDistance = std::numeric_limits<double>::max();
    
    public:
        SCluster(std::vector<SHit> hitList);

        int NHits(){return fNHits;}
        void AddHit(SHit& hit);
        std::vector<SHit> GetHits(){return fHitList;};

        double GetMinDistanceToCluster(SHit h);
        double GetMinDistanceToClusterW(SHit h);
        double GetMinDistanceToClusterOverlap(SHit h);
        double GetMinDistanceToCluster1D(SHit h);
        std::pair<SHit, double> GetClosestHitToPoint(SHit h);

        friend std::ostream& operator<<(std::ostream& out, SCluster & c);

};


#endif // TPC_SIMPLE_CLUSTERS_H