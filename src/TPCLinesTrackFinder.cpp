////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesTrackFinder.h
//
// \brief Definition of TPCLinesTrackFinder
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>


#include "TPCLinesTrackFinder.h"
#include "DistanceUtils.h"
#include "TPCLinesDBSCAN.cpp"

// Constructor
TPCLinesTrackFinder::TPCLinesTrackFinder(TrackFinderAlgorithmPsetType tpcLinesTrackFinderPset):
    fTPCLinesTrackFinderPset(tpcLinesTrackFinderPset),
    fDisplay(TPCLinesDisplay())
{}


double TPCLinesTrackFinder::GetHitLinearDensity(std::vector<SHit> data){
    int minX = 1e6;
    int maxX = -1e6;
    for (SHit& h : data) {
        if(h.X()<minX) minX = h.X();
        if(h.X()>maxX) maxX = h.X();
    }
    
    //std::cout<<"MinX "<<minX<<" MaxX "<<maxX<<" "<<data.size()<<" "<<(maxX+1-minX)<<std::endl;
    return 1.*data.size()/(maxX+1-minX);
}
    

LineEquation TPCLinesTrackFinder::PerformPCA2D(std::vector<SHit>& data) {

    // first center the data
    double meanX = 0.0, meanY = 0.0;
    for (SHit& h : data) {
        meanX += h.X();
        meanY += h.Y();
    }
    meanX /= data.size();
    meanY /= data.size();

    // get the covariance and fill centered vectors
    double covXX = 0.0, covYY = 0.0, covXY = 0.0;
    for (SHit& h : data) {
        double x = h.X()-meanX;
        double y = h.Y()-meanY;
        covXX += x * x;
        covYY += y * y;
        covXY += x * y;
    }
    covXX /= data.size();
    covYY /= data.size();
    covXY /= data.size();

    double eigenvalue1 = (covXX + covYY + std::sqrt((covXX - covYY) * (covXX - covYY) + 4 * covXY * covXY)) / 2.0;
    double eigenvalue2 = (covXX + covYY - std::sqrt((covXX - covYY) * (covXX - covYY) + 4 * covXY * covXY)) / 2.0;

    double eigenvectorX1 = eigenvalue1 - covYY;
    double eigenvectorY1 = covXY;
    double eigenvectorX2 = eigenvalue2 - covYY;
    double eigenvectorY2 = covXY;

    std::cout << "Eigenvalues: " << eigenvalue1 << ", " << eigenvalue2 << std::endl;
    std::cout << "Eigenvectors: (" << eigenvectorX1 << ", " << eigenvectorY1 << "), (" << eigenvectorX2 << ", " << eigenvectorY2 << ")" << std::endl;

    // Calculate slopes (m) and intercepts (b) for lines
    double slope1 = eigenvectorY1 / eigenvectorX1;
    double yIntercept1 = meanY - slope1 * meanX; // Since the line passes through the mean
    /*double slope2 = eigenvectorY2 / eigenvectorX2;
    double yIntercept2 = 0.0; // Since the line passes through the origin (0, 0)*/

    std::cout << "Slope-Intercept Form of Lines:" << std::endl;
    std::cout << "Line 1: y = " << slope1 << "x + " << yIntercept1 << std::endl;

    return LineEquation(slope1, yIntercept1);
}


// Get hits in the tube
void TPCLinesTrackFinder::GetHitsInTube(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD) {
    hitsIn.clear();
    hitsOut.clear();

    for (auto& hit : hitList) {
        SPoint pProj = line.GetLineClosestPoint( SPoint(hit.X(), hit.Y()) );
        double d = std::sqrt(std::pow(hit.X() - pProj.X(), 2) + std::pow(hit.Y() - pProj.Y(), 2));
        if (d < maxD) {
            hitsIn.push_back(hit);
        } else {
            hitsOut.push_back(hit);
        }
    }

    return;
}


void TPCLinesTrackFinder::GetHitsInTubeSingleWire(std::vector<SHit> & hitsIn, std::vector<SHit> & hitsOut, LineEquation& line, std::vector<SHit>& hitList, double maxD){
    hitsIn.clear();
    hitsOut.clear();

    std::unordered_map<int, int> hitIndexDict;
    std::unordered_map<int, double> hitDistDict;
    std::unordered_map<int, std::vector<SHit>> discardedhitsDict;

    for (int hitix = 0; hitix < hitList.size(); hitix++) {
        SHit hit = hitList[hitix];
        SPoint pProj = line.GetLineClosestPoint( SPoint(hit.X(), hit.Y()) );
        double d = std::sqrt(std::pow(hit.X() - pProj.X(), 2) + std::pow(hit.Y() - pProj.Y(), 2));
        if (d < maxD) {
            int xIndex = static_cast<int>(hit.X());
            //std::cout<<hit<<std::endl;
            if (hitIndexDict.count(xIndex) > 0) {
                if (d < hitDistDict[xIndex]) {
                    int refHitIndex = hitIndexDict[static_cast<int>(hit.X())];
                    SHit oldHit = hitsIn[refHitIndex];
                    hitsIn[refHitIndex] = hit;
                    hitDistDict[xIndex] = d;
                    discardedhitsDict[xIndex].push_back(oldHit);
                    //std::cout<<" Found index, Replacing\n";
                } else {
                    hitsOut.push_back(hit);
                    //std::cout<<"Found index, Not replacing\n";
                }
            } else {
                //std::cout<<" Not found index, adding\n";
                hitsIn.push_back(hit);
                hitIndexDict[xIndex] = hitsIn.size() - 1;
                hitDistDict[xIndex] = d;
            }
        }
        else {
            hitsOut.push_back(hit);
        }
    }

    for (auto& entry : discardedhitsDict) {
        int wireix = entry.first;
        auto& discardedHitList = entry.second;
        int refHitIndex = hitIndexDict[wireix];
        std::vector<SHit> refHitV = { hitsIn[refHitIndex] };
        std::sort(discardedHitList.begin(), discardedHitList.end(), [&](SHit& a, SHit& b) {
            return std::abs(a.Y() - refHitV[0].Y()) < std::abs(b.Y() - refHitV[0].Y());
        });

        for (auto& hit : discardedHitList) {
            bool doOverlap = false;
            for (auto& refHit : refHitV) {
                if (HitWidthOverlap(hit, refHit)) {
                    doOverlap = true;
                }
            }
            if (doOverlap) {
                refHitV.push_back(hit);
                hitsIn.push_back(hit);
            } else {
                hitsOut.push_back(hit);
            }
        }
    }

    return;
}

// Function to get frequency of elements in a vector
map<int, int> GetListFrequency(vector<int> vList) {
    map<int, int> frequency_dict;
    for (int element : vList) {
        if (element == -1) continue;
        frequency_dict[element]++;
    }
    return frequency_dict;
}


// 2DClusters
void TPCLinesTrackFinder::Get2DClusters(std::vector<SHit> & hits, double epsilon, int minPts){


    DBSCAN dbscan(epsilon, minPts);
    //dbscan.setDistanceFunction(manhattanDistance);  // You can change this to euclideanDistance or any other custom function
    dbscan.setDistanceFunction(GetpClusterDistanceW);
    dbscan.fit(hits);
    
    std::vector<int>& clusterAssignment = dbscan.getClusterAssignment();
    for (int i = 0; i < hits.size(); ++i) {
        if (clusterAssignment[i] == -1) {
            std::cout << "SPoint (" << hits[i].X() << ", " << hits[i].Y() << ") is noise" << std::endl;
        } else {
            std::cout << "SPoint (" << hits[i].X() << ", " << hits[i].Y() << ") belongs to cluster " << clusterAssignment[i] << std::endl;
        }
    }

    return 0;
}



// Main function
std::vector<SCluster> TPCLinesTrackFinder::ReconstructTracksFromHoughDirection(std::vector<SHit> & hitList, LineEquation houghLine, int trkIter) {
    std::vector<SCluster> recoTracks;

    // Get hits in the Hough tube
    std::vector<SHit> hitHoughTubeList, hitOutHoughTubeList;
    GetHitsInTube(hitHoughTubeList, hitOutHoughTubeList, houghLine, hitList,  fTPCLinesTrackFinderPset.MaxDTube);
    double hitDensity = GetHitLinearDensity(hitHoughTubeList);
    std::cout << "   Hit density: " << hitDensity << std::endl;

    if(fTPCLinesTrackFinderPset.Verbose>=2){
        fDisplay.Show("Hough direction", hitList, houghLine, hitHoughTubeList);
    }
    

    // Get hits in the refined PCA direction
    LineEquation pcaLine = PerformPCA2D(hitHoughTubeList);
    std::vector<SHit> hitPCATubeList, hitPCAOutTubeList;
    //hitPCATubeList = GetHitsInTube(pcaLine, hitList,  fTPCLinesTrackFinderPset.MaxDCluster);
    if (!fTPCLinesTrackFinderPset.SingleWireMode && hitDensity > 2) {
        GetHitsInTube(hitPCATubeList, hitPCAOutTubeList, pcaLine, hitList,  fTPCLinesTrackFinderPset.MaxDCluster);
    }
    else {
        std::cout<<"Single wire mode\n";
        GetHitsInTubeSingleWire(hitPCATubeList, hitPCAOutTubeList, pcaLine, hitList,  fTPCLinesTrackFinderPset.MaxDCluster);
    }

    if(fTPCLinesTrackFinderPset.Verbose>=2){
        fDisplay.Show("PCA direction", hitList, pcaLine, hitPCATubeList);
    }

    // Create cluster with PCA hits
    SCluster pcaCluster(hitPCATubeList);
    if(pcaCluster.NHits()< fTPCLinesTrackFinderPset.MinTrackHits) return recoTracks;

    // CLUSTERING: CONNECTEDNESS
    std::cout<<"\n******** Making connectedness clusters\n";
    std::cout<<pcaCluster<<std::endl;

    double epsilon = 1.8268;
    std::vector<SHit> pcaClusterHits = pcaCluster.GetHits();
    Get2DClusters(pcaClusterHits, epsilon, fTPCLinesTrackFinderPset.MinTrackHits);

    
    recoTracks.push_back(pcaCluster);

    // return the unmatched hits
    std::cout<<" Remaining hits: "<<hitPCAOutTubeList.size()<<std::endl;
    hitList.clear();
    for (SHit& h : hitPCAOutTubeList) {
        hitList.push_back(h);
    }
    
    return recoTracks;
}
