////////////////////////////////////////////////////////////////////////////
//
// \file TPCLinesAlgo.h
//
// \brief Definition of TPCLinesDisplay
//
// \author fjnicolas@ugr.es
//
////////////////////////////////////////////////////////////////////////////
#ifndef TPC_LINES_DISPLAY_H
#define TPC_LINES_DISPLAY_H

#include "TGraphErrors.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "TPCSimpleHits.h"

class TPCLinesDisplay {
    private:
        void DrawHitScatter(std::vector<SHit> hitV, int color, int style, double size, double errorAlpha);
        void DrawLine(LineEquation line, double xmin, double xmax, int color, int style);
        TH2F GetFrame(std::vector<SHit> hitsV);
    public:
        TPCLinesDisplay(){};
        void Show(
            std::string eventLabel,
            std::vector<SHit> allHitsV,
            LineEquation houghLine,
            std::vector<SHit> selectedHitsV);
};


TH2F TPCLinesDisplay::GetFrame(std::vector<SHit> hitsV){

    std::vector<double> x, y, err;
    for(size_t ix=0; ix<hitsV.size(); ix++){
        x.push_back(hitsV[ix].X());
        y.push_back(hitsV[ix].Y());
        err.push_back(hitsV[ix].Width());
    }
    double minX = *std::min_element(x.begin(), x.end());
    double maxX = *std::max_element(x.begin(), x.end());
    double minY = *std::min_element(y.begin(), y.end());
    double maxY = *std::max_element(y.begin(), y.end());
    double overX = 0.1*(maxX-minX);
    double overY = 0.1*(maxY-minY);

    TH2F hFrame("frame", ";WireIndex;TimeIndex", 200, minX-overX, maxX+overX, 200, minY-overY, maxY+overY);
    hFrame.SetStats(0);
    return hFrame;
}

void TPCLinesDisplay::DrawLine(LineEquation line, double xmin, double xmax, int color, int style){
    // Draw a horizontal line from x = 1 to x = 4 at y = 2
    double y1 = line.EvaluateX(xmin);
    double y2 = line.EvaluateX(xmax);
    TLine* horizontalLine = new TLine(xmin, y1, xmax, y2);
    std::cout<<" Drawing line "<<xmin<<" "<<y1<<" "<<xmax<<" "<<y2<<std::endl;
    horizontalLine->SetLineColor(color); // Set line color
    horizontalLine->SetLineStyle(style); // Set line color
    horizontalLine->SetLineWidth(2);     // Set line width
    horizontalLine->Draw();
}

void TPCLinesDisplay::DrawHitScatter(std::vector<SHit> hitsV, int color, int style, double size, double errorAlpha){


    std::vector<double> x, y, err;
    for(size_t ix=0; ix<hitsV.size(); ix++){
        //std::cout<<ix<<" "<<hitsV[ix].X()<<" "<<hitsV[ix].Y()<<std::endl;
        x.push_back(hitsV[ix].X());
        y.push_back(hitsV[ix].Y());
        err.push_back(hitsV[ix].Width());
    }

    TGraph *g = new TGraphErrors(x.size(),&x[0],&y[0], 0, &err[0]); 
    g->SetMarkerColorAlpha(color, 0.5);
    g->SetMarkerStyle(style);
    g->SetMarkerSize(size);
    g->SetLineColorAlpha(kGray, errorAlpha);
    g->Draw("p");
    
    return;
}



void TPCLinesDisplay::Show(
    std::string eventLabel,
    std::vector<SHit> allHitsV,
    LineEquation houghLine,
    std::vector<SHit> selectedHitsV)
{


    TCanvas c("c", eventLabel.c_str(), 0, 0, 800, 800);
    c.Divide(1,1);
    c.cd(1);

    // general frame
    TH2F hFrame=GetFrame(allHitsV);
    hFrame.Draw();

    // all hits scatter
    DrawHitScatter(allHitsV, 65, 8, 1, 0.6);

    // selected hits scatter
    DrawHitScatter(selectedHitsV, kRed, 24, 1.1, 0);

    // hough line
    DrawLine(houghLine, hFrame.GetXaxis()->GetXmin(), hFrame.GetXaxis()->GetXmax(), kBlack, kDashed);

    c.cd();
    c.Update();c.WaitPrimitive(); 

    return;

}

#endif // TPC_SIMPLE_LINES_H