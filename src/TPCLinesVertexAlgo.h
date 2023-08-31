
////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesVertexAlgo
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPCLINES_VERTEXALGO_H
#define TPCLINES_VERTEXALGO_H

#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "TPCSimpleHits.h"
#include "TPCSimpleLines.h"
#include "TPCSimpleTriangles.h"




class TPCLinesVertexAlgo {
    private:



    public:
        // constructor
        TPCLinesVertexAlgo(){};

        // main function
        void GetOrigins(std::vector<SLinearCluster> trackList, std::vector<STriangle>& vertexList, std::vector<SPoint> &originList, SLinearCluster &mainDirection);

        
};

#endif // TPC_SIMPLE_LINES_H
