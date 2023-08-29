////////////////////////////////////////////////////////////////////////////
//
// \file DistanceUtils.h
//
// \brief Definition of fistsance functions
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef DISTANCE_UTILS_H
#define DISTANCE_UTILS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

#include "TPCSimpleHits.h"

// Mean and Std Dev function for std::vector
double CalculateMean( std::vector<double>& values) {

    double sum = 0;
    for ( double& value : values) {
        sum += value;
    }
    return sum / values.size();
}

double CalculateStdDev( std::vector<double>& values) {

    double mean = CalculateMean(values);
    double sumSqDiff = 0;
    for ( double& value : values) {
        double diff = value - mean;
        sumSqDiff += diff * diff;
    }
    return std::sqrt(sumSqDiff / values.size());
}


// Hit distance funcitons
double GetHitDistance( SHit hit1,  SHit hit2) {

    double xdif = hit1.X() - hit2.X();
    double ydif = hit1.Y() - hit2.Y();

    return std::sqrt( xdif*xdif + ydif*ydif);
}

double GetHitDistance1D( SHit hit1,  SHit hit2) {
    return std::abs(hit1.X() - hit2.X());
}

double GetHitDistanceW( SHit hit1,  SHit hit2) {
    
    // distance in X
    double dX = std::pow(hit1.X() - hit2.X(), 2);
    
    // distance between centers
    double d0 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y(), 2));
    
    // dist of hit 1 to width of hit 2
    double dp1 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() + hit2.Width(), 2));
    double dm1 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() - hit2.Width(), 2));
    double d1 = std::min(dp1, dm1);
    
    // dist of hit 2  to width of hit 1
    double dp2 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() + hit1.Width(), 2));
    double dm2 = std::sqrt(dX + std::pow(hit1.Y() - hit2.Y() - hit1.Width(), 2));
    double d2 = std::min(dp2, dm2);
    
    // return min
    return std::min(d0, (d1 + d2) / 2.);
}

double GetHitDistanceOverlap( SHit& hit1,  SHit& hit2) {
    double y_range1_min = hit1.Y() - hit1.Width();
    double y_range1_max = hit1.Y() + hit1.Width();
    double y_range2_min = hit2.Y() - hit2.Width();
    double y_range2_max = hit2.Y() + hit2.Width();
    bool overlap = (y_range1_min <= y_range2_max) && (y_range1_max >= y_range2_min);
    double dX = std::pow(hit1.X() - hit2.X(), 2);
    double dY = 1;
    if (!overlap) {
        dY = std::pow(hit1.Y() - hit2.Y(), 2);
    }
    return std::sqrt(dX + dY);
}


bool HitWidthOverlap(SHit& hit1, SHit& hit2) {
    std::pair<double, double> y_range1 = {hit1.Y() - hit1.Width(), hit1.Y() + hit1.Width()};
    std::pair<double, double> y_range2 = {hit2.Y() - hit2.Width(), hit2.Y() + hit2.Width()};
    
    return y_range1.first <= y_range2.second && y_range1.second >= y_range2.first;
}







/*

double GetpClusterDistanceW(SPoint p1, SPoint p2) {
    double d0 = std::sqrt(std::pow(p1[0] - p2[0], 2) + std::pow(p1[1] - p2[1], 2));
    double dp = std::sqrt(std::pow(p1[0] - p2[0], 2) + std::pow(p1[1] - p2[1] + p2[2], 2));
    double dm = std::sqrt(std::pow(p1[0] - p2[0], 2) + std::pow(p1[1] - p2[1] - p2[2], 2));
    return std::min(d0, std::min(dp, dm));
}

double GetpClusterDistanceOverlap(SPoint p1, SPoint p2) {
    double y_range1_min = p1[1] - p1[2];
    double y_range1_max = p1[1] + p1[2];
    double y_range2_min = p2[1] - p2[2];
    double y_range2_max = p2[1] + p2[2];
    bool overlap = (y_range1_min <= y_range2_max) && (y_range1_max >= y_range2_min);
    double dX = p1[0] - p2[0];
    double dY = 1;
    if (!overlap) {
        dY = std::sqrt(std::pow(p1[0] - p2[0], 2) + std::pow(p1[1] - p2[1], 2));
    }
    return std::sqrt(std::pow(dX, 2) + std::pow(dY, 2));
}

double GetClusterMinDistance(const Cluster& cluster1, const Cluster& cluster2) {
    double minDistance = 1e4;
    for (const auto& hit1 : cluster1.HitList) {
        for (const auto& hit2 : cluster2.HitList) {
            double distance = std::sqrt(std::pow(hit1.X - hit2.X, 2) + std::pow(hit1.Y - hit2.Y, 2));
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
    }
    return minDistance;
}

double GetClusterConnectedness(const Cluster& cluster1, const Cluster& cluster2) {
    double minDistance = 1e4;
    for (const auto& hit1 : cluster1.HitList) {
        for (const auto& hit2 : cluster2.HitList) {
            double distance = GetHitDistanceW(hit1, hit2);
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
    }
    return minDistance;
}

double GetClusterConnectednessOverlap(const Cluster& cluster1, const Cluster& cluster2) {
    double minDistance = 1e4;
    for (const auto& hit1 : cluster1.HitList) {
        for (const auto& hit2 : cluster2.HitList) {
            double distance = GetHitDistanceOverlap(hit1, hit2);
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
    }
    return minDistance;
}

void MergeIsolatedHits2(std::vector<LinearCluster>& recoTrackList, std::vector<Hit>& hitList, double dCleaning) {
    std::vector<std::vector<Hit>> trackHitDict(recoTrackList.size());
    for (size_t ix = 0; ix < recoTrackList.size(); ++ix) {
        trackHitDict[ix] = {};
    }

    for (const auto& hit : hitList) {
        std::vector<int> candidateTracksIx;
        for (size_t ix = 0; ix < recoTrackList.size(); ++ix) {
            double d = recoTrackList[ix].GetDistance(hit);
            if (d < dCleaning) {
                candidateTracksIx.push_back(ix);
            }
        }
        std::cout << "\n" << hit.X << " " << hit.Y << " CandidateTracks ";
        for (const auto& ix : candidateTracksIx) {
            std::cout << ix << " ";
        }
        std::cout << std::endl;

        if (!candidateTracksIx.empty()) {
            double minDCluster = 1e4;
            int selectedTrackIx = -1;
            for (const auto& ix : candidateTracksIx) {
                double minD = 1e4;
                for (const auto& hitT : recoTrackList[ix].HitCluster.HitList) {
                    double d = GetHitDistanceW(hit, hitT);
                    if (d < minD) {
                        minD = d;
                    }
                }
                std::cout << "Track: " << ix << " dTrack " << minD << " DCluster: " << minD << std::endl;
                std::cout << "Track conn " << recoTrackList[ix].HitCluster.Connectedness << " Comp " << recoTrackList[ix].HitCluster.Compactness << std::endl;

                if (minD < minDCluster) {
                    selectedTrackIx = ix;
                    minDCluster = minD;
                }

                std::cout << selectedTrackIx << " " << minDCluster << std::endl;
            }
            double trackConn = recoTrackList[selectedTrackIx].HitCluster.Connectedness;
            double trackConnRMS = recoTrackList[selectedTrackIx].HitCluster.ConnectednessRMS;
            if (minDCluster < 2.5 * trackConn) {
                std::cout << "Adding to track " << selectedTrackIx << std::endl;
                recoTrackList[selectedTrackIx].HitCluster.HitList.push_back(hit);
                recoTrackList[selectedTrackIx].HitCluster.Connectedness = trackConn;
                recoTrackList[selectedTrackIx].HitCluster.ConnectednessRMS = trackConnRMS;
            }
        }
    }
}

void MergeIsolatedHits(std::vector<LinearCluster>& recoTrackList, std::vector<Hit>& hitList, double dCleaning1D, double dTh = 3) {
    std::cout << "\n\n +-+-+-+-+-+-+- Isolated hit merger +-+-+-+-+-+-+-" << std::endl;
    std::vector<LinearCluster> newRecoTrackList;

    std::vector<std::vector<int>> hitToTrackDict(hitList.size());
    std::vector<bool> mergedHitsCounter(hitList.size(), false);
    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        hitToTrackDict[hitix] = {};
    }

    for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
        const auto& hit = hitList[hitix];
        double minD = 1e3;
        for (size_t trkix = 0; trkix < recoTrackList.size(); ++trkix) {
            const auto& track = recoTrackList[trkix];
            double trackComp = track.HitCluster.Compactness;
            double hitTrackDist = track.HitCluster.GetMinDistanceToCluster1D(hit);
            if (hitTrackDist < dCleaning1D) {
                double d = track.TrackEquation.GetDistance(hit);
                double hypoY = track.TrackEquation.Slope * hit.X + track.TrackEquation.Intercept;
                if (d < dTh * trackComp && d < minD) {
                    minD = d;
                    hitToTrackDict[hitix].push_back(trkix);
                }
            }
        }
    }

    for (size_t trkix = 0; trkix < recoTrackList.size(); ++trkix) {
        const auto& track = recoTrackList[trkix];
        double trackComp = track.HitCluster.Compactness;
        double trackConn = track.HitCluster.Connectedness;
        std::cout << "\n Merging track " << trkix << " COMP " << trackComp << " CONN " << trackConn << std::endl;

        std::vector<Hit> candidateHits;
        std::vector<int> candidateHitsIx;
        for (size_t hitix = 0; hitix < hitList.size(); ++hitix) {
            if (!mergedHitsCounter[hitix]) {
                if (std::find(hitToTrackDict[hitix].begin(), hitToTrackDict[hitix].end(), trkix) != hitToTrackDict[hitix].end()) {
                    candidateHits.push_back(hitList[hitix]);
                    candidateHitsIx.push_back(hitix);
                }
            }
        }

        std::vector<double> candidateHitsDist;
        for (const auto& hit : candidateHits) {
            double distToCluster = track.HitCluster.GetMinDistanceToCluster(hit);
            candidateHitsDist.push_back(distToCluster);
        }

        std::vector<Hit> sortedCandidateHits(candidateHits.size());
        std::vector<int> sortedCandidateHitsIx(candidateHitsIx.size());
        std::iota(sortedCandidateHits.begin(), sortedCandidateHits.end(), 0);
        std::sort(sortedCandidateHits.begin(), sortedCandidateHits.end(), [&candidateHitsDist](int i, int j) {
            return candidateHitsDist[i] < candidateHitsDist[j];
        });
        std::sort(sortedCandidateHitsIx.begin(), sortedCandidateHitsIx.end(), [&candidateHitsDist](int i, int j) {
            return candidateHitsDist[i] < candidateHitsDist[j];
        });

        std::vector<Hit> newHitList = track.HitCluster.HitList;
        Cluster currentCluster = track.HitCluster;

        std::cout << "  candidates hits: " << std::endl;
        for (size_t ix = 0; ix < candidateHits.size(); ++ix) {
            const auto& hit = candidateHits[sortedCandidateHits[ix]];
            double minD = currentCluster.GetMinDistanceToCluster(hit);
            double minDConn = currentCluster.GetMinDistanceToClusterOverlap(hit);
            std::cout << hit.X << " " << hit.Y << " " << minD << " " << minDConn << std::endl;
            std::cout << "       " << hit.Id << " " << hit.X << " " << hit.Y << " d " << minD << " dConn " << minDConn << std::endl;
            if (minDConn < dTh * trackConn) {
                newHitList.push_back(hit);
                mergedHitsCounter[candidateHitsIx[sortedCandidateHits[ix]]] = true;
                currentCluster.HitList = newHitList;
            }
        }

        newRecoTrackList.push_back({currentCluster, track.Slope, track.Intercept});
    }

    recoTrackList = newRecoTrackList;
}

int main() {
    std::vector<LinearCluster> recoTrackList;
    std::vector<Hit> hitList;
    double dCleaning = 0.0;
    MergeIsolatedHits2(recoTrackList, hitList, dCleaning);
    MergeIsolatedHits(recoTrackList, hitList, dCleaning);
    return 0;
}

Copy*/

#endif // TPC_SIMPLE_CLUSTERS_H