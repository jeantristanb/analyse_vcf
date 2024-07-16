# -*- coding: latin-1 -*-
import dadi
import numpy
import pylab
import glob
import os
import sys
from numpy import array
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerLine2D

def GetPatternPopHW(Data1D) :
    ns = (Data1D.sample_sizes[0],)
    ### model le plus simple et neutre
    func = dadi.Demographics1D.snm
    ### dans le cas d'un spectre de frequence 1d
    pts = max(ns)
    ### quel params and why? depend du model
    params = None
    ## estimation du spectre de frequence
    Estimsnm=func(params,ns,pts)
    theta = dadi.Inference.optimal_sfs_scaling(Estimsnm, Data1D)
    return (Estimsnm, theta)


def GetMat2DSF(FileFreq2D, Pop1, Pop2) :
    Lire=open(FileFreq2D)
    FsI=[]
    Sum=0
    CmtCol=0
    Mat1DPop1=[]
    Mat1DPop2=[]
    for Ligne in Lire :
       FsI.append([int(x) for x in Ligne.split()])
       Sum+=sum(FsI[CmtCol])
       if CmtCol==0 :
          Mat1DPop2=FsI[CmtCol]
       else :
          Mat1DPop2=array(Mat1DPop2) + array(FsI[CmtCol])
       Mat1DPop1.append(sum(FsI[CmtCol]))
       CmtCol+=1
    Lire.close()
    FsIInversed=[]
    #for ligne in reversed(FsI):
    #    FsIInversed.append(ligne[::-1])
    print Mat1DPop1, Mat1DPop2
    return [dadi.Spectrum(FsI, pop_ids=[Pop1,Pop2]),dadi.Spectrum(list(Mat1DPop1),pop_ids=[Pop1]), dadi.Spectrum(list(Mat1DPop2),pop_ids=[Pop2]),Sum]

def Plot2D(Mat,File, SaveFig=True):
    dadi.Plotting.plot_single_2d_sfs(Mat, vmin=1)
    if SaveFig :
       pylab.savefig(File, dpi=100)
       pylab.close()


def WriteMatSFS2D(Mat,File) :
    Chaine=""
    for x in Mat :
       Chaine+="\t".join([str(y)  for y in x])+"\n"
    Chaine=Chaine.replace("--","0")
    Ecrire=open(File, 'w')
    Ecrire.write(Chaine)
    Ecrire.close()

def WriteRaportMat1D(Mat, Entete=False) :
    if Entete==True :
       return["D","S","Pi","Wat","n"]
    if Mat!=None :
       return [Mat.Tajima_D(),Mat.S(), Mat.pi() ,Mat.Watterson_theta(), len(Mat)-1]
    
def PlotSFS1(Mat, Col, Legend):
    Tmp=[x for x in Mat]
    Tmp=Tmp[1:(len(Tmp)-1)]
    TmpS=sum(Tmp)
    Tmp=[x/float(TmpS) for x in Tmp]
    return plt.plot(range(1,(len(Tmp)+1)), Tmp, Col, label=Legend)


def PlotSFSWithEst(Mat, File) :
    (MatTheo,ns)=GetPatternPopHW(Mat)
    dadi.Plotting.plot_1d_comp_Poisson(MatTheo*ns, Mat)
    plt.savefig(File , dpi=100)
    pylab.close()

def PlotSFSCmp2(Mat1DPop1, Mat1DPop2,  NomPop1, NomPop2, Entete) :
   ## #Comparaison des SFS
   MinNbInd=min(len(Mat1DPop2), len(Mat1DPop1))-1
   MatProj1=Mat1DPop1.project([MinNbInd])
   MatProj2=Mat1DPop2.project([MinNbInd])
   (Mat1DTheo,ns)=GetPatternPopHW(MatProj2)
   line1,=PlotSFS1(MatProj1, 'b', NomPop1)
   line1,=PlotSFS1(MatProj2 ,'r', NomPop2)
   line1,=PlotSFS1(Mat1DTheo, 'k', "Attendu")
   plt.legend(handler_map={line1: HandlerLine2D(numpoints=4)})
   pylab.savefig(Entete+"_CmpSFS1D.jpeg", dpi=100)
   pylab.close()
   Char=NomPop1+"\t"+"\t".join([str(x) for x in MatProj1])+"\n"
   Char+="\t"+"\t".join([str(x) for x in MatProj2])+"\n"
   Char+="\t"+"\t".join([str(x) for x in Mat1DTheo])+"\n"
   Ecrire=open(Entete+"_CmpSFS1D.txt" ,'w')
   Ecrire.write(Char)




#FileFreq2D="/home/jeantristan/Travail/ExtractDataForABC/SFS/AllTrip/MexCBDNFCaraibe/SFS/SFS2d.CBD.Car.out"
FileFreq2D=sys.argv[1]
DirSortie="Figure/"
SplitFreqName=FileFreq2D.split(".")
NomPop1=SplitFreqName[-3]
NomPop2=SplitFreqName[-2]

[Mat2D,Mat1DPop1,Mat1DPop2, NbSnpAll]=GetMat2DSF(FileFreq2D, NomPop1, NomPop2)

Dir=DirSortie+"/"
try :
  os.makedirs(Dir)
except :
  print "Dir "+ Dir + " Exist "

MatTheo=Mat2D.scramble_pop_ids()
Entete=DirSortie+"/Fig_"+NomPop1+"_"+NomPop2

plt.subplot(211)
Plot2D(MatTheo, Entete+'_scramble.jpeg', False)
plt.subplot(212)
Plot2D(Mat2D,Entete+'_2D.jpeg')
WriteMatSFS2D(MatTheo,Entete+'_scramble.out')
WriteMatSFS2D(Mat2D,Entete+'_i.out')
PlotSFSCmp2(Mat1DPop1,Mat1DPop2,  NomPop1, NomPop2, Entete)

Chaine="NomPop\t"+"\t".join([str(x) for x in WriteRaportMat1D(Mat1DPop1, Entete=True)])+"\tNbSnpAll\n"
Chaine+=NomPop1+"\t"+"\t".join([str(x) for x in WriteRaportMat1D(Mat1DPop1)])+"\t"+str(NbSnpAll)+"\n"
Chaine+=NomPop2+"\t"+"\t".join([str(x) for x in WriteRaportMat1D(Mat1DPop2)])+"\t"+str(NbSnpAll)+"\n"
Ecrire=open(Entete+"_statres.txt" , 'w')
Ecrire.write(Chaine)
