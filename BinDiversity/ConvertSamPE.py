#!/usr/bin/python
import fileinput
import sys

def PrintFormatIllumina(EcrireSortie,ligneSam):
   splitline=line.split()
   EcrireSortie.write("@"+splitline[0]+"\n")
   EcrireSortie.write(splitline[9]+"\n")
   EcrireSortie.write("+\n")
   EcrireSortie.write(splitline[10]+"\n")


if len(sys.argv)!=3 :
   print "Exe SortieR1 SortieR2"
   sys.exit()


try :
    EcrireR1=open(sys.argv[1], 'w')
except :
   print "file "+ sys.argv[1] + " Can't open"


try :
    EcrireR2=open(sys.argv[2], 'w')
except :
   print "file "+ sys.argv[2] + " Can't open"

print sys.argv[2] + " " + sys.argv[1]

cmt=1
BaliseEcrire=False
for line in sys.stdin:
    if line[0]!='@' or BaliseEcrire:
       if cmt%2==0: 
          PrintFormatIllumina(EcrireR2,line)
       else :
          PrintFormatIllumina(EcrireR1,line)
       BaliseEcrire=True
       cmt+=1

EcrireR1.close() 
EcrireR2.close() 

    
