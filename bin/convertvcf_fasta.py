from utils import *
import argparse



def parseArguments():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--vcf',type=str,required=True)
    parser.add_argument('--out',type=str,required=True)
    parser.add_argument('--indel',type=str,required=False)
    parser.add_argument('--chr',type=str,required=False)
    parser.add_argument('--begin', type=int,required=False)
    parser.add_argument('--end', type=int,required=False)
    args = parser.parse_args()
    return args

args=parseArguments()
UseIndel=False
if args.indel and args.indel[0]=="T" :
 UseIndel=True


Header=GetHeaderVcf(args.vcf)
Header=Header[len(Header)-1]
listhead=Header.split()
listhead=[x.upper() for x in listhead[9::]]
counthead=range(0,len(listhead))
fastseq=[['',''] for x in counthead]

ReadVcf=OpenFile(args.vcf)
if args.chr and args.begin and args.end :
    CheckChro=True
    chroch=args.chr
    beginch=args.begin
    endch=args.end

        
cmtpos=0
for line in ReadVcf :
 if line[0]!="#" :
   spll=line.split()
   if CheckChro :
      pos=int(spll[1])
      if chroch!=spll[0] or pos<beginch or endch>endch:
         continue
   infovcf=GetInfoVcf(line)
   for cmti in counthead:
       A1=infovcf[4][cmti][0]
       if not A1 or (not UseIndel and len(A1)>0):
          A1='N'
       A2=infovcf[4][cmti][1]
       if not A2 or (not UseIndel and len(A2)>0):
          A2='N'
       fastseq[cmti][0]+=A1
       fastseq[cmti][1]+=A2
      
Fasta=""
for cmtind in counthead :
   Fasta+=">"+listhead[cmtind]+"_A1\n"+fastseq[cmti][0]+"\n" 
   Fasta+=">"+listhead[cmtind]+"_A2\n"+fastseq[cmti][1]+"\n" 

open(args.out)
