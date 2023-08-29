////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleClusters.cpp
//
// \brief Definition of SimpleTPCClusters
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCSimpleClusters.h"
#include "DistanceUtils.h"

SCluster::SCluster(std::vector<SHit> hitList){
    fNHits = hitList.size();

    // sort hits by X coordinate
    std::sort(hitList.begin(), hitList.end(), [](SHit& a, SHit& b) {
        return a.X() < b.X();
    });
    if (fNHits > 0) {

        for(auto &hit: hitList){
            double dminComp = DefaultMaxDistance;
            double dminConn = DefaultMaxDistance;
            double dminX = DefaultMaxDistance;
            
            for(auto &hit2: hitList){
                if (&hit == &hit2) continue;

                // distance between centers               
                double d = GetHitDistance(hit, hit2);
                if(d<dminComp) dminComp=d;

                // calculate connectedeness
                d = GetHitDistanceW(hit, hit2);
                //d = GetHitDistanceOverlap(hit, hit2)
                if(d<dminConn) dminConn=d;

                // calculate connectedness 1D
                d = GetHitDistance1D(hit, hit2);
                if(d<dminX) dminX=d;
            }

            SHit newHit = hit;
            newHit.SetHitConnectivity(dminComp, dminConn, dminX);

            fCompactnessV.push_back(dminComp);
            fConnectednessV.push_back(dminConn);
            fConnectedness1DV.push_back(dminX);
            fWidths.push_back(newHit.Width());
            fHitList.push_back(newHit);
        }
            
        fCompactness = CalculateMean(fCompactnessV);
        fCompactnessRMS = CalculateStdDev(fCompactnessV);
        fConnectedness = CalculateMean(fConnectednessV);
        fConnectednessRMS = CalculateStdDev(fConnectednessV);
        fConnectedness1D = CalculateMean(fConnectedness1DV);
        fConnectedness1DRMS = CalculateStdDev(fConnectedness1DV);
        fAverageWidth = CalculateMean(fWidths);
    }
    else {
        fCompactness = 0;
        fCompactnessRMS = 0;
        fConnectedness = 0;
        fConnectednessRMS = 0;
        fConnectedness1D = 0;
        fConnectedness1DRMS = 0;
        fAverageWidth = 0;
    }
}

void SCluster::AddHit(SHit& hit) {
    fHitList.push_back(hit);
    fNHits++;
}

std::ostream& operator<<(std::ostream& out, SCluster & c)
{
    out <<  " --- Created Simple Cluster:\n"
    << "  Comp=" << c.fCompactness << " pm " << c.fCompactnessRMS
    << "  Conn=" << c.fConnectedness << " pm " << c.fConnectednessRMS
    << "  Conn1D=" <<c.fConnectedness1D << " pm " << c.fConnectedness1DRMS;

    return out;
}


double SCluster::GetMinDistanceToCluster(SHit h) {
    double minD = 1e6;
    for (auto& hit : fHitList) {
        double d = GetHitDistance(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}


double SCluster::GetMinDistanceToClusterW(SHit h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = GetHitDistanceW(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

double SCluster::GetMinDistanceToClusterOverlap(SHit h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = GetHitDistanceOverlap(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}

double SCluster::GetMinDistanceToCluster1D(SHit h) {
    double minD = 1e4;
    for (auto& hit : fHitList) {
        double d = GetHitDistance1D(hit, h);
        if (d < minD) minD = d;
    }
    return minD;
}


std::pair<SHit, double> SCluster::GetClosestHitToPoint(SHit h) {
    double minD = 1e4;
    size_t closestHitIx = 0;

    for (size_t ix =0; ix<fHitList.size(); ix++){
        SHit hit = fHitList[ix];
        double d = GetHitDistance(hit, h);
        if (d < minD){
            minD = d;
            closestHitIx = ix;
        }
    }

    std::pair<SHit, double> result(fHitList[closestHitIx], minD);
    return result;
}