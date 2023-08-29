#!/usr/bin/env python3mport numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import matplotlib.patches as patches
import os

plt_params = {'legend.fontsize': 'medium',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'large',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'large'}


def GetHitDictFromList(hitList):
    X = []; Y=[]; Z=[]; Wi=[]; ST=[]; ET=[]
    for ix in range(len(hitList)):
        X.append(hitList[ix].X)
        Y.append(hitList[ix].Y)
        Z.append(hitList[ix].Integral)
        Wi.append(hitList[ix].Width)
        ST.append(hitList[ix].StartT)
        ET.append(hitList[ix].EndT)

    return {"X":X,"Y":Y,"Z":Z,"Wi":Wi, "ST":ST, "ET":ET}

def Display(eventLabel, hitList, lineEq=HoughLine(), hitSelectedList=[], clusterTrackV=[], vertex=SVertex(), distV=[], recoVertexList=[], trackEdges=[], parJuntions=[], mainTrack=None, angleList=None, invisibleTracks=[], show=False, outputDir="plots/"):

    fPlotStartEndPoints = False
    plt.rcParams.update(plt_params)
    fig = plt.figure(eventLabel, figsize=(8, 8))
    labelStr=""
    if(isinstance(lineEq, HoughLine)):
        labelStr+="All hits"
    gs = GridSpec(4, 4, hspace=0, wspace=0)
    ax1 = fig.add_subplot(gs[0:4, 0:4])
    
    #### convert to dicitonaries
    hitDict = GetHitDictFromList(hitList)
    hitSelectedDict = GetHitDictFromList(hitSelectedList)
    
    height = max(hitDict["Y"]) - min(hitDict["Y"])
    width = max(hitDict["X"]) - min(hitDict["X"])
    L=math.sqrt(width*width+height*height)

    alpha=1.
    if(vertex.Active==True):
        alpha=0.6
    #### plot all the active hits
    ax1.scatter(hitDict["X"], hitDict["Y"], label=labelStr, clim=[-1, 2], alpha=alpha)
    for ix in range(len(hitDict["X"])):
        
        # plot start/end points
        if(fPlotStartEndPoints==True):
            ax1.plot( [hitDict["X"][ix], hitDict["X"][ix]], [hitDict["ST"][ix], hitDict["ET"][ix]], c="lightgrey", alpha=alpha/2)
        
        # plot width
        ax1.plot( [hitDict["X"][ix], hitDict["X"][ix]], [hitDict["Y"][ix]-hitDict["Wi"][ix], hitDict["Y"][ix]+hitDict["Wi"][ix]], c="grey", alpha=alpha)
    #### plot the selected hits
    if(len(hitSelectedDict)>0):
        ax1.scatter(hitSelectedDict["X"], hitSelectedDict["Y"], marker="x", c="red", label="Selected hits")
    

    #### Plot hough tracks
    if( isinstance(lineEq, HoughLine)):
        if(lineEq.Score>=-1):
            m = lineEq.Equation.m
            n = lineEq.Equation.n
            labelHough="Hough dir (Score="+str(lineEq.Score)+", NHits="+str(lineEq.NHits)+")"
            ax1.plot( [-L, L], [-m*L+n, m*L+n], c="black", label=labelHough, ls="-." )
    else:
        m = lineEq.m
        n = lineEq.n
        ax1.plot( [-L, L], [-m*L+n, m*L+n], c="black", ls="-." )

    #### Plot the clusters
    for cIx in range(len(clusterTrackV)):
        clusterHitDict = GetHitDictFromList( clusterTrackV[cIx].HitCluster.HitList )
        clusterLabel = "Cluster:"+str(cIx)
        if(clusterTrackV[cIx].GetId()!=-1): clusterLabel = "Cluster:"+str(clusterTrackV[cIx].GetId())
        ax1.scatter(clusterHitDict["X"], clusterHitDict["Y"], marker="o", facecolors='none', edgecolor="C"+str(cIx+1), linewidth=1.5, label=clusterLabel)

        if(clusterTrackV[cIx].HasResidualHits):
            clusterResHitDict = GetHitDictFromList( clusterTrackV[cIx].ResidualHitCluster.HitList )
            if(len(clusterResHitDict["X"])>0):
                ax1.scatter(clusterResHitDict["X"], clusterResHitDict["Y"], color="grey", marker="x", s=50)

        
        fPlotStartEndSlopes=True
        if(fPlotStartEndSlopes):
            Xpoints = [clusterTrackV[cIx].minX, clusterTrackV[cIx].maxX]
            
            m1=clusterTrackV[cIx].TrackEquationStart.m
            n1=clusterTrackV[cIx].TrackEquationStart.n
            Ypoints = m1*np.array(Xpoints)+n1
            ax1.plot( Xpoints, Ypoints, c="C"+str(cIx+1) )

            m1=clusterTrackV[cIx].TrackEquationEnd.m
            n1=clusterTrackV[cIx].TrackEquationEnd.n
            Ypoints = m1*np.array(Xpoints)+n1
            ax1.plot( Xpoints, Ypoints, c="C"+str(cIx+1) )
            
            m1=clusterTrackV[cIx].TrackEquation.m
            n1=clusterTrackV[cIx].TrackEquation.n
            Ypoints = m1*np.array(Xpoints)+n1
            ax1.plot( [-L, L], [-m1*L+n1, m1*L+n1], c="C"+str(cIx+1), alpha=0.7, ls=":" )

        # Plot fit
        xV = np.arange(clusterTrackV[cIx].minX-6, clusterTrackV[cIx].maxX+6, 1)
        if hasattr(clusterTrackV[cIx], "TrackFit"):
            yV = clusterTrackV[cIx].TrackFit.evaluate(xV)
            ax1.plot( xV, yV, c="C"+str(cIx+1), alpha=0.7, ls=":" )

        m1=clusterTrackV[cIx].TrackEquation.m
        n1=clusterTrackV[cIx].TrackEquation.n
        yV1 = m1*np.array(xV)+n1
        ax1.plot(xV , yV1, c="C"+str(cIx+1), alpha=0.7, ls="-." )
        
        track = clusterTrackV[cIx]
        if(track.HasStartEndPoints==True):
            vertexLabel = ""
            if(cIx==0): 
                vertexLabel = "Edge hits"
            ax1.scatter([track.StartPoint.X], [track.StartPoint.Y], edgecolor="fuchsia", facecolors='none', s=60, label=vertexLabel)
            ax1.scatter([track.EndPoint.X], [track.EndPoint.Y], edgecolor="fuchsia", facecolors='none', s=60)


    #### plot the vertex
    if(vertex.Active == True):
        ax1.scatter([vertex.P.X], [vertex.P.Y], c="red", label="PANDORA vertex", marker="*")

    #### plot the reco vertexes
    for vix, v in enumerate(recoVertexList):
        vertexLabel = ""
        if(vix==0): vertexLabel = "Intersections ("+str(len(recoVertexList))+")"
        ax1.scatter([v.X], [v.Y], c="cyan", label=vertexLabel, marker="o")

    #### plot the parallel junctions
    for jix, junc in enumerate(parJuntions):
        from matplotlib.patches import FancyArrow
        juncLabel = ""
        if(jix==0): juncLabel = "ParJunctions"
        ax1.arrow(junc[0][0], junc[0][1], junc[1][0]-junc[0][0], junc[1][1]-junc[0][1], 
        head_width=0.2, head_length=0.2, fc='black', ec='black', label=juncLabel)

    #### plot the main track
    if(mainTrack!=None):
        clusterHitDict = GetHitDictFromList(mainTrack.HitCluster.HitList )
        ax1.scatter(clusterHitDict["X"], clusterHitDict["Y"], edgecolor="black", facecolors='none', s=60, label="MainTrack")

    #### plot triangles
    if(angleList!=None):
        for triangle in angleList:
            print("Drawing triangle", triangle.MainVertex.X, triangle.MainVertex.Y, triangle.VertexB.X, triangle.VertexB.Y)
            trianglePatch = patches.Polygon([(triangle.MainVertex.X, triangle.MainVertex.Y),
                                             (triangle.VertexB.X, triangle.VertexB.Y),
                                             (triangle.VertexC.X, triangle.VertexC.Y)], closed=True, edgecolor='b',color="grey", alpha=0.6, fill=True)
            
            overlayL=40
            XArrow = np.array([triangle.MidPoint.X, triangle.MainVertex.X])
            XArrow = np.array(XArrow) 
            XArrowOver = [XArrow[0], XArrow[1]+overlayL] if(triangle.MidPoint.X<triangle.MainVertex.X) else [XArrow[0]-overlayL, XArrow[1]]
            XArrowOver = np.array(XArrowOver)
            
            YArrow = triangle.Direction.Slope()*XArrow + triangle.Direction.Intercept()
            YArrowOver = triangle.Direction.Slope()*XArrowOver + triangle.Direction.Intercept()

            arrowDir = (XArrow[1]-XArrow[0], YArrow[1]-YArrow[0])

            print(arrowDir)

            arrow = patches.FancyArrowPatch((XArrow[0], YArrow[0]), (XArrow[1], YArrow[1]),  mutation_scale=5, linewidth=1)
            ax1.add_patch(arrow)
            ax1.plot((XArrowOver[0], XArrowOver[1]), (YArrowOver[0], YArrowOver[1]),  c="black", linewidth=1)


            XArrowOver = [triangle.MainVertex.X-overlayL, triangle.MainVertex.X+overlayL]
            XArrowOver = np.array(XArrowOver)
            ax1.plot(XArrowOver, triangle.MomentumHypo1.Slope()*XArrowOver+triangle.MomentumHypo1.Intercept() ,  c="gold", linewidth=1)
            ax1.plot(XArrowOver, triangle.MomentumHypo2.Slope()*XArrowOver+triangle.MomentumHypo2.Intercept() ,  c="gold", linewidth=1)

            ## draw the vertex hit

            ax1.scatter([triangle.MainVertexHit.X], [triangle.MainVertexHit.Y], marker="*", c="midnightblue",s=20)


            ax1.add_patch(trianglePatch)

    #### plot invisible tracks
    for ix, invTrack in enumerate(invisibleTracks):

            m=invTrack.m
            n=invTrack.n

            x0 = recoVertexList[ix].X

            Xpoints = np.arange(x0-10, x0+10, 1)
            Ypoints = m*np.array(Xpoints)+n

            print(Xpoints, Ypoints)
            
            labelName = ""
            if(ix==0):
                labelName = "InvisibleTrack"

            ax1.plot( Xpoints, Ypoints, c="midnightblue", label=labelName)



    ax1.set_xlabel("WireIx")
    ax1.set_ylabel("TimeIx")
    ax1.set_xlim(-2, max(hitDict["X"])); ax1.set_ylim(-2, max(hitDict["Y"]))
    ax1.legend()

    if(len(distV)>0):
        ax2 = fig.add_subplot(gs[1, 3])
        ax2.hist(distV)
        ax2.set_xlabel("d")


    directory_path = outputDir
    # Check if the directory exists
    if not os.path.exists(directory_path):
        # Create the directory if it doesn't exist
        os.makedirs(directory_path)
    
    fig.savefig(directory_path+eventLabel+".pdf")
    fig.savefig("/Users/franciscojaviernicolas/Work/scratch_plots/DisplayTPCLines.pdf")

    if(show==True):
        plt.show()


def HelloWorldDisplay(x):
    print("Hello World TPCLINES!", x)

    print("PLOTTING")


    fig = plt.figure("AAA", figsize=(8, 8))
    plt.plot([1, 2], [1, 2])
    print("PLOTTING")
    plt.show()


    return 

HelloWorldDisplay(2)

