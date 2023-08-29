////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesAlgo
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////

#include "TPCLinesAlgo.h"

// Constructor
TPCLinesAlgo::TPCLinesAlgo(TPCLinesAlgoPsetType tpcLinesAlgoPset, std::string displayPath):
    fTPCLinesPset(tpcLinesAlgoPset),
    fNTotalHits(0),
    fHitList({}),
    fVertex(SPoint(-1, -1), ""),
    fHoughAlgo(tpcLinesAlgoPset.HoughAlgorithmPset),
    fTrackFinder(tpcLinesAlgoPset.TrackFinderAlgorithmPset),
    fDisplayAppPath(displayPath),
    fDisplay(TPCLinesDisplay())
{}


// Ana function
void TPCLinesAlgo::SetHitList(std::string view,
                            std::vector<int>& vertex,
                            std::vector<int> *_X,
                            std::vector<double> *_Y,
                            std::vector<double> *_Int,
                            std::vector<double> *_Wi,
                            std::vector<double> *_ST,
                            std::vector<double> *_ET,
                            std::string eventLabel)
{   
    // reset variables
    fHitList.clear();
    fNTotalHits=0;

    size_t nTotalHits = _X->size();

    // Define the vertex
    double vertexX = vertex[2];
    if (view == "U0" || view == "U1") vertexX = vertex[0];
    if (view == "V0" || view == "V1") vertexX = vertex[1];
    double vertexY = vertex[3] / fTPCLinesPset.DriftConversion;
    std::cout << "  ** Vertex XY: " << vertexX << " " << vertexY << std::endl;
    if (vertexX != -1){

        // Get hits in the selected plane
        std::vector<double> filteredX, filteredY, filteredInt, filteredST, filteredET, filteredWi;
        for (int i = 0; i < nTotalHits; i++) {
            int x = _X->at(i);
            
            // filter channels for the view
            if ( x > fChB[view][0] && x <= fChB[view][1]) {
                double y = _Y->at(i)/fTPCLinesPset.DriftConversion;
                
                // filter distance to vertex
                double d = std::sqrt(std::pow(x - vertexX, 2) + std::pow(y - vertexY, 2));
                if (d < fTPCLinesPset.MaxRadius) {
                    filteredX.push_back( x );
                    filteredY.push_back( y );
                    filteredInt.push_back( _Int->at(i) );
                    filteredWi.push_back( _Wi->at(i) / fTPCLinesPset.DriftConversion );
                    filteredST.push_back( _ST->at(i) / fTPCLinesPset.DriftConversion );
                    filteredET.push_back( _ET->at(i) / fTPCLinesPset.DriftConversion );
                }
            }
        }

        if (filteredX.size() > 3) {
            // Shift hits to have origin in (0,0)
            double minX = *std::min_element(filteredX.begin(), filteredX.end()) - 3;
            double minY = *std::min_element(filteredY.begin(), filteredY.end()) - 3;

            
            for (int i = 0; i < filteredX.size(); i++) {
                SHit hit(i, filteredX[i] - minX, filteredY[i] - minY, filteredWi[i], filteredInt[i], filteredST[i] - minY, filteredET[i] - minY);
                fHitList.push_back(hit);
            }
            fNTotalHits = fHitList.size();

            // create vertex with common origin
            vertexX = vertexX - minX;
            vertexY = vertexY - minY;

            fVertex = SVertex(SPoint(vertexX, vertexY), view);
            std::cout << "  ** Origin vertex: " << fVertex;
            std::cout << "  ** NInputHits:"<<fNTotalHits<<std::endl;

        }
        else {
            std::cout << "   SKIPPED NHits selected near the vertex " << filteredX.size() << std::endl;
        }

    }
}


std::map<std::string, double> TPCLinesAlgo::AnaView()
{
    // Map for the final results
    std::map<std::string, double> anaResults;

    // loop over the hough tracks
    int trkIter = 0;
    std::vector<SHit> hitListForHough = fHitList;

    while(trkIter<fTPCLinesPset.MaxHoughTracks){

        std::cout<<" **** Track finder iteration: "<<trkIter<< " # hough hits:"<<hitListForHough.size()<<std::endl;

        // Get the best Hough line       
        HoughLine houghLine = fHoughAlgo.GetBestHoughLine(hitListForHough, fVertex);
        std::cout<<"    Hough line results: "<<houghLine.NHits()<<" Score: "<<houghLine.Score()<<std::endl;

        // Skip if not enough hits
        if(houghLine.NHits() < fTPCLinesPset.MinTrackHits){
            std::cout<<"   Hough lines has <"<<fTPCLinesPset.MinTrackHits<<", ending the search\n";
            trkIter = fTPCLinesPset.MaxHoughTracks;
        }

        // Call the track/cluster algorithm 
        std::vector<SCluster> linearClusterV = fTrackFinder.ReconstructTracksFromHoughDirection(hitListForHough, houghLine.GetLineEquation(), trkIter);

        fDisplay.Show("Final Algo", fHitList, houghLine.GetLineEquation(), hitListForHough);

        // Skip if not enough hits
        if(hitListForHough.size() < fTPCLinesPset.MinTrackHits){
            std::cout<<"   Remaining hits are "<<hitListForHough.size()<<", ending the search\n";
            trkIter = fTPCLinesPset.MaxHoughTracks;
        }

                                   
        //if(houghLine.NHits < houghPset.MinHoughHits): continue

        trkIter++;
    }
    Display();
    return anaResults;
}


void TPCLinesAlgo::Display(){


    fDisplay.Show(
        "aa", 
        fHitList,
        LineEquation(0, 0),
        fHitList);

    if(fDisplayAppPath!="acxds"){
        std::string command = "python3 "+fDisplayAppPath+"/DisplayTPCLines.py";
        int status = system(command.c_str());
    }
    

    return;
}
