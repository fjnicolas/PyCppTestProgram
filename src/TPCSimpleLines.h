////////////////////////////////////////////////////////////////////////////
//
// \file TPCLineObjects.h
//
// \brief Definition of TPCLinebjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_LINES_H
#define TPC_SIMPLE_LINES_H

#include <cmath>
#include "TPCSimpleHits.h" 

class LineEquation {
    private:
        float fM;
        float fN;

    public:
        LineEquation(float slope, float intercept);
        
        float Slope() {return fM;}
        float Intercept() {return fN;}

        SPoint GetLineClosestPoint(double a, double b, double c, SPoint P);
        SPoint GetLineClosestPoint(SPoint P);
        float GetDistance(SPoint P);
        float Evaluate(SPoint p);
        float EvaluateX(double x);
};


class HoughLine {
    private:
        LineEquation fEquation;
        float fScore;
        int fNHits;
    public:
        HoughLine(LineEquation line = LineEquation(0, 0), float score = -1, int nHits = 0);
        
        void SetLineEquation(LineEquation eq){ fEquation = eq;}
        void SetScore(float sc){ fScore = sc;}
        void SetNHits(float nhits){ fNHits = nhits;}

        float Score(){return fScore;}
        float NHits(){return fNHits;}
        LineEquation GetLineEquation(){return fEquation;}
};


#endif // TPC_SIMPLE_LINES_H