# -*- coding: utf-8 -*-
import sys
import os
import gzip
from fct_compdiv import *
from fct_parsevcf import *

def initialise_div(PosVcfIndPop):
    ListeCmtPop=range(len(PosVcfIndPop))
    NbPop=len(PosVcfIndPop)
    MaxNbInd=max([len(x) for x in PosVcfIndPop])
    PreComputea1=ComputePiHaploide.PrecomputingDa1(MaxNbInd)
    PreComputea2=ComputePiHaploide.PrecomputingDa2(MaxNbInd)
    AllDtaj=[[0 for x in range(7)] for x in range(NbPop)]
    AllFst=[[[0 for x in range(4)] for x in range(NbPop)] for x in range(NbPop)]
    return(ListeCmtPop,NbPop,MaxNbInd,PreComputea1, PreComputea2, AllDtaj, AllFst)

def intialisefile(headerout, DicInfoPop):
    PopName=DicInfoPop['PopIndic']
    NbPop=len(PopName)
    ChaineHeader=""
    NameHeaderDtaj=['Dtj', 'S', 'pi', 'TetWat', 'n', 'NbPos', 'NbPoly']
    for Pop in range(0,len(PopName)):
      for Header in NameHeaderDtaj :
       ChaineHeader+= "\t"+Header+"_"+PopName[Pop]
    ChaineEnteteFst=""
    CmtPop=0
    while CmtPop<NbPop-1 :
      CmtPop2=CmtPop+1
      while CmtPop2<NbPop :
         for Entete in NameEnteteFst :
             ChaineEnteteFst+="\t"+Entete+"_"+NamePop[CmtPop]+"_"+NamePop[CmtPop2]
         CmtPop2+=1
      CmtPop+=1
 
    WriteDtajFen=open(headerout+".wind.dtaj",'w')
    WriteDtajFen.write("Chro\tPosDeb\tPosFin"+NomHeaderDtaj+"\n")
    WriteFstFen=open(headerout+".wind.fst",'w')
    WriteFstFen.write("Chro\tPosDeb\tPosFin"+ChaineHeaderFst+"\n")
    return (WriteDtajFen, WriteFstFen)
    


def computedstat_allvcf(VcfFileList,PosVcfIndPop,freqminna,IsGZ)
    (ListeCmtPop,NbPop,MaxNbInd,PreComputea1, PreComputea2, AllDtaj, AllFst)=initialise_div(PosVcfIndPop)
    ## for each vcf file
    for FileVCF in ListeVCF :



