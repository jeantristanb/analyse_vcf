import ComputePiHaploide
import sys 
import os 
import random
import numpy

def GetNa(Mat, NbErreur) :
    if NbErreur== 0 :
       return Mat
    CmtNbErreur=0
    NbCol=len(Mat[0])-1
    NbLigne=len(Mat)-1
    while CmtNbErreur<NbErreur:
         Col=random.randint(0,NbCol)
         Ligne=random.randint(0,NbLigne)
         if Mat[Ligne][Col]!= -1:
            Mat[Ligne][Col] = -1
            CmtNbErreur+=1
    return Mat
    

NMax=10
NSnp=1000
NbSim=1000
MinFreqNa=0.5
NbSnpI=NMax*NSnp

PreComputea1=ComputePiHaploide.PrecomputingDa1(NMax)
PreComputea2=ComputePiHaploide.PrecomputingDa2(NMax)

Ecrire=open("TestD/TestD_"+str(NMax)+"_"+str(NSnp)+"_"+str(MinFreqNa)+".out",'w')
for CmtSim in range(1,NbSim):
    NbErreur=0
    Lire=os.popen("/home/jeantristan/Travail/bin/ms "+ str(NMax)+ " 1 -s "+str(NSnp))
    MatData=numpy.mat([[int(y)+1 for y in x.strip() ] for x in Lire.readlines()[6::]]).T
    MatData=MatData.tolist()
    Nb1=[x.count(1) for x in MatData]
    Nb2=[x.count(2) for x in MatData]
    #ComputePiHaploide.ComputeDivDtajima([Nb1.count(x) for x in range(0,NMax+1)])
    Lire.close()
    DtajI=ComputePiHaploide.TajimaDWithMissingData(Nb1,Nb2,PreComputea1,PreComputea2,NMax, MinFreqNa)
    for FreqCmt in [0.05, 0.1, 0.2, 0.3,0.49] :
        NbErreur=int(NbSnpI*FreqCmt) - NbErreur
        MatData=GetNa(MatData, NbErreur)
        Nb1=[x.count(1) for x in MatData]
        Nb2=[x.count(2) for x in MatData]
        DtajI.append(ComputePiHaploide.TajimaDWithMissingData(Nb1,Nb2,PreComputea1,PreComputea2,NMax, MinFreqNa)[0])
    Ecrire.write("\t".join([str(x) for x in DtajI])+"\n")
    Lire.close()


Ecrire.close()
    
