#!/usr/bin/python
import sys

if(len(sys.argv)!=3):
    print "exe sortiesamtoolstat SortieBilan"
    sys.exit()


try :
    Fichier=sys.argv[1]
    Lire=open(Fichier)
except:
    print "impossible d'ouvrir "+Fichier 
    sys.exit()

try :
    FichierSortie=sys.argv[2]
    ecrire=open(FichierSortie, "w")
except:
    print "impossible d'ouvrir "+FichierSortie
    sys.exit()

#Fichier="ResumeSamTool-Ak3.txt"
#Lire=open(Fichier)

DonneeTemp=Lire.readlines()
Donnee=[DonneeTemp[0]]
Donnee+=DonneeTemp[3::]
Lire.close()
Parser="\t"

Tot=float(Donnee[0].split()[0])
NbDuplique=float(Donnee[1].split()[0])
Mapped=float(Donnee[2].split()[0])
#Paired2=Donnee[3].split()[0]
BienPairedBonInsert=float(Donnee[6].split()[0])
BienPairedAll=float(Donnee[7].split()[0])
Singleton=float(Donnee[8].split()[0])
DiffChro=float(Donnee[9].split()[0])
DiffChroQ5=float(Donnee[10].split()[0])

ecrire.write("Total_Seq"+Parser+str(Tot)+Parser+str(Tot/Tot)+Parser+str(Tot/2.0)+"\n")
ecrire.write("Nb_Seq_Mapp"+Parser+str(Mapped)+Parser+str((Mapped)/Tot)+Parser+str(Singleton+BienPairedAll/2.0) +"\n")
ecrire.write("Nb_Seq_Non_Mappe"+Parser+str(Tot-Mapped)+Parser+str((Tot-Mapped)/Tot)+ Parser+str((Tot-Mapped)/2.0)+"\n")
ecrire.write("Nb_Seq_Dup"+Parser+str(NbDuplique)+Parser+str((NbDuplique)/Tot)+Parser+str(NbDuplique/2.0)+"\n")
ecrire.write("Nb_SeqMap_BonInsert"+Parser+ str(BienPairedBonInsert)+Parser+str((BienPairedBonInsert)/Tot)+Parser+str(BienPairedBonInsert/2.0) +"\n")
ecrire.write("Nb_Seq_Bien_aligne_Tout_Insert"+Parser+ str(BienPairedAll)+Parser+str(BienPairedAll/Tot)+ Parser+ str(BienPairedAll/2.0)+ "\n") 
ecrire.write("Nb_Singleton"+Parser +str(Singleton)+Parser+str((Singleton)/Tot)+Parser+str(Singleton)+"\n")
ecrire.write("Nb_Map_Diff_Chro"+Parser+str(DiffChro)+Parser+str((DiffChro)/Tot)+Parser+str((DiffChro)/2.0)+"\n")
ecrire.write("Nb_Map_Diff_Chro_Q_Sup_5"+Parser+str(DiffChroQ5) +Parser+str((DiffChroQ5)/Tot)+Parser+str((DiffChroQ5)/2.0)+"\n")
PEBumapped=Tot/2.0-(Singleton+BienPairedAll/2.0)
ecrire.write("Nb_Seq_Both_Unm"+Parser+str(PEBumapped*2) +Parser+str((PEBumapped*2)/Tot)+Parser+str(PEBumapped)+"\n")
ecrire.close()
