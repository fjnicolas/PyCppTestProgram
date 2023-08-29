////////////////////////////////////////////////////////////////////////////
//
// \file TPCSimpleObjects.h
//
// \brief Definition of SimpleTPCObjects
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_SIMPLE_HITS_H
#define TPC_SIMPLE_HITS_H

#include <iostream>
#include <string>
#include <iosfwd>

class SPoint {
    private:
        double fX;
        double fY;

    public:
        SPoint(int x, int y);
        SPoint(double x, double y);
        double X(){return fX;}
        double Y(){return fY;}

        friend std::ostream& operator<<(std::ostream& out, SPoint const& p);
};

class SVertex {
    private:
        SPoint fP;
        std::string fView;
        bool fActive;

    public:
        SVertex();
        SVertex(SPoint p, std::string view);
        
        SPoint Point(){return fP;}
        std::string View(){return fView;}
        bool IsActive(){return fActive;}

        double X(){return fP.X();}
        double Y(){return fP.Y();}

        friend std::ostream& operator<<(std::ostream& out, SVertex const& v);
};

class SHit {
    private:
        int fId;
        float fX;
        float fY;
        float fWidth;
        float fStartT;
        float fEndT;
        float fIntegral;
        float fXProj;
        float fYProj;
        float fCompactness;
        float fConnectedness;
        float fConnectednes1D;

    public:
        // Constructor
        SHit(int id, float x, float y, float w, float integral, float st, float et);

        int Id(){return fId;}
        float X(){return fX;}
        float Y(){return fY;}
        float Width(){return fWidth;}
        float StartT(){return fStartT;}
        float EndT(){return fEndT;}
        float Integral(){return fIntegral;}
        float XProj(){return fXProj;}
        float YProj(){return fYProj;}
        float Compactness(){return fCompactness;}
        float Connectednes(){return fConnectedness;}
        float Connectednes1D(){return fConnectednes1D;}

        // Set the hit connectivity
        void SetHitConnectivity(float comp, float conn, float conn1d);

        // Overload cout
        friend std::ostream& operator<<(std::ostream& out, SHit const& h);
};

#endif // TPC_SIMPLE_HITS_H