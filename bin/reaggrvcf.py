#!/usr/bin/env python3
from utils import *
from utils_div import *
import sys
import argparse

def parseArguments():
    parser = argparse.ArgumentParser(description='extract annotation for specific position')
    parser.add_argument('--vcfi',type=str,required=True, help="file contains chro and list of files with annotation")
    parser.add_argument('--vcfout', type=str,help="output")
    args = parser.parse_args()
    return args



def getvcf(File):
    headvcf=GetHeaderVcf(File)
    readvcf=OpenFile(File, Type='r')
    InfoDic={}
    InfoPos={}
    Head=""
    for line in readvcf:
      if line[0]!="#":
        (Chro, Pos, Ref, ListAlt,Geno)=GetInfoVcf(line, None)
        if  Chro not in InfoDic :
           InfoDic[Chro]={}
           InfoPos[Chro]={}
        if Pos not in InfoDic[Chro].keys():
          InfoDic[Chro][Pos]=[[Chro, Pos, Ref, ListAlt,Geno, line]]
        else :
          InfoDic[Chro][Pos].append([Chro, Pos, Ref, ListAlt,Geno, line])
      else :
        Head+=line 
    return (InfoDic ,Head)

def MergePosMulti(Chro, Pos, PosMultiInfo):
    ## check that ref is ok
    Ref=PosMultiInfo[0][2]
    Pos=PosMultiInfo[0][1]
    ListAll=[]
    for x in PosMultiInfo :
      if x[2]!=Ref :
         print('error ref for pos '+str(Pos)+' refI '+Ref+' ref2 '+ x[2])
         sys.exit()
      ListAll+=x[3]
    ListNewGeno=[[Ref,Ref] for x in range(len(PosMultiInfo[0][4]))]
    ListNewGeno2=[[0,0] for x in range(len(PosMultiInfo[0][4]))]
    ListNewGenoF=["" for x in range(len(PosMultiInfo[0][4]))]
    for Cmtind in range(len(ListNewGeno)):
        balisefound=True
        for CmtPos in range(len(PosMultiInfo)):
          InfoInf=PosMultiInfo[CmtPos][4][Cmtind]
          if InfoInf[0]==None :
            balisefound=False   
          else :
            if InfoInf[0]!=Ref :
               ListNewGeno[Cmtind][0]=InfoInf[0]     
               ListNewGeno2[Cmtind][0]=ListAll.index(InfoInf[0])
            if InfoInf[1]!=Ref :
               ListNewGeno[Cmtind][1]=InfoInf[1]
               ListNewGeno2[Cmtind][1]=ListAll.index(InfoInf[1])
        if InfoInf[4]==0 :
           sep="/"
        else :
           sep="|"
        if balisefound :
          ListNewGenoF[Cmtind]=str(ListNewGeno2[Cmtind][0])+sep+str(ListNewGeno2[Cmtind][1])
        else :
          ListNewGenoF[Cmtind]="."
    #19      45353086        .       C       T       .       PASS    .       GT      
    return "\t".join([Chro,Pos, ".", Ref, ",".join(ListAll), ".", "PASS", ".", "GT"])+"\t"+"\t".join(ListNewGenoF)+"\n"

    
    

def mergepos(Dic, Head, FileOut) :
    writevcf=open(FileOut, 'w')
    writevcf.write(Head)
    for chro in Dic.keys():
      for pos in Dic[chro] :
         if len(Dic[chro][pos])>2:
            writevcf.write(MergePosMulti(chro, pos,Dic[chro][pos]))
         else :
            writevcf.write(InfoDic[chro][pos][0][5])
    

args = parseArguments()
## first 
(InfoDic, HeadVcf)=getvcf(args.vcfi)
mergepos(InfoDic, HeadVcf, args.vcfout)



