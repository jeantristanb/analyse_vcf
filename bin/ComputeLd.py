#!/usr/bin/python
# coding=utf-8
from utils import *
from utils_div import *
import argparse
import re
import os
import sys
import glob
def parseArguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--vcf',type=str,required=True)
    #parser.add_argument('--rs',type=str,required=False)
    parser.add_argument('--bp',type=str,required=False)
    parser.add_argument('--out',type=str,required=False, default='out')
    parser.add_argument('--chr',type=str,required=False)
    args = parser.parse_args()
    return args


def getvcf(File):
    headvcf=GetHeaderVcf(File) 
    readvcf=OpenFile(File, Type='r')
    InfoDic={}
    InfoPos={}
    for line in readvcf:
      if line[0]!="#":
        (Chro, Pos, Ref, ListAlt,Geno)=GetInfoVcf(line, None)
        if  Chro not in InfoDic :
           InfoDic[Chro]=[]
           InfoPos[Chro]=[]
        InfoDic[Chro].append([Chro, Pos, Ref, ListAlt,Geno])
        InfoPos[Chro].append(Pos)
    return (InfoDic,InfoPos)

def ComputedLDPosChro(AllInfo,InfoPos, Chro=None, Pos=None):
   if (Chro not in InfoPos) and (Pos not in InfoPos):
       sys.exit('error Chro '+Chro+' Pos '+str(Pos)+'not found')       
   #compute_geno_bin(Snp1,Snp2, MinFreqBon=0)
   # (Chro, Pos, Ref, ListeAlt,Geno)
   ChroInfo=AllInfo[Chro]
   Geno1=AllInfo[Chro]
   #for Pos2 in InfoPos[Chro] :
   #compute_geno_bin(,Snp2, MinFreqBon=0)

#def ComputedLD (AllInfo,InfoPos ):
#  for Chro in InfoPos:
#   NbPos=len(InfoPos[Chro])
#   AllInfChro=AllInfo[Chro]
#   for CmtPos1 in range(0,NbPos-1) :
#     for CmtPos2 in range(CmtPos1+1,NbPos-1) :
#         Snp1=AllInfChro[CmtPos1]
#         Snp2=AllInfChro[CmtPos2]
#         if len(Snp1[3])==1 and len(Snp2[3])==1:
#            print(compute_geno_bin(Snp1[4],Snp2[4])) 
#mcld

def LaunchMCLD( AllInfo,InfoPos ,out,tmpfile='.akfmlksdvnfdk'):
  AllW=open(out, 'w')
  #Loci   Aver. |Delta-prime|   Aver. |Corr.|     Approx. p-value   Permut. p-value
  Head=["Chro","Pos1", "Pos2", "Alt1", "Alt2", "Dp_mcld", "R2_mcld", "AppPval_mcld", "Dp", "Rp", "p1", "p2", "n1"]
  AllHead="\t".join(Head)
  AllW.write(AllHead+'\n')
  for Chro in InfoPos: 
    ChroInfo=AllInfo[Chro]
    ListChro=InfoPos[Chro]
    NInd=len(ChroInfo[0][4])
    NLocus=len(ChroInfo)
    AllChaine=""
    for CmtInd in range(NInd) :
       InfoInd=""
       for CmtPos in range(NLocus):
          InfoInd+=" "+str(ChroInfo[CmtPos][4][CmtInd][2])+" "+str(ChroInfo[CmtPos][4][CmtInd][3])
       AllChaine+=InfoInd+"\n" 
    writeall=open(tmpfile,'w')
    writeall.write(AllChaine)
    writeall.close()
    os.system(dirpyth+"/mcld -file="+tmpfile+ ' -perm=1 > '+tmpfile+".out")
    readfil=open(tmpfile+".out")
    Head=readfil.readline() 
    for line in readfil :
      spllin=line.split()
      numpos=spllin[0].split('/')
      CmtPos1=int(numpos[0])-1
      CmtPos2=int(numpos[1])-1
      Pos1=ListChro[CmtPos1]
      Pos2=ListChro[CmtPos2]
      Snp1=ChroInfo[CmtPos1]
      Snp2=ChroInfo[CmtPos2]
      # (Dprime,Rcar,pa1,pb1, n)
      if len(Snp1[3])==1 and  len(Snp2[3])==1:
        res=compute_geno_bin(Snp1[4],Snp2[4]) 
      else :
        res=["NA"]*5
      Alt1=",".join(Snp1[3])
      Alt2=",".join(Snp2[3])
      ResumResPos=[Chro,Pos1,Pos2,Alt1, Alt2,spllin[1],spllin[2],spllin[3], res[0], res[1], res[2], res[3],res[4]]
      AllW.write("\t".join([str(x) for x in ResumResPos])+"\n")
    
      

dirpyth=os.path. dirname(sys.argv[0])
args=parseArguments()
(AllInfo,PosInfo)=getvcf(args.vcf)
LaunchMCLD( AllInfo,PosInfo,args.out,tmpfile='.akfmlksdvnfdk')
#ComputedLD (AllInfo,PosInfo)
#if args.bp and args.chr :
#    ComputedLDPosChro(AllInfo,PosInfo, args.chr, args.bp)

