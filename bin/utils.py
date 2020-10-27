# coding=utf-8
import gzip

def Str(Number, Lim=2):
    if type(Number)==float :
       return str(round(Number,Lim))
    return str(Number)

def Concatener(Vecteur, sep="\t"):
    Chaine=""
    for El in Vecteur :
       Chaine+=str(El)+sep
    return Chaine


def GetMatrix(Col, Row):
    return [[0 for x in xrange(Col)] for x in xrange(Row)]
def intersect(liste1, liste2):
    return list(set(liste1) & set(liste2))

def diff(liste1, liste2):
    return list(set(liste1) - set(liste2))

def GetInfoVcf(lines, postosave=None):
    def getgeno(x,l):
      if x[0]=='.' :
        return [None, None,0,0, None]
      x=x.split(':')[0]#re.split([':|/'])
      if '/' in x :
         x=x.split('/')
         phase=0 
      elif '|' in x:
         x=x.split('|')
         phase=1 
      return [l[int(x[0])], l[int(x[1])], int(x[0])+1, int(x[0])+1,phase]
    TmpSplit=lines.split()
    Chro=TmpSplit[0]
    Pos=TmpSplit[1]
    NomPos=Chro+"_"+Pos
    Ref=TmpSplit[3]
    ListeAlt=TmpSplit[4].split(',')
    TmpPos=[Ref]
    TmpPos+=ListeAlt
    if postosave :
     return (Chro, Pos, Ref, ListeAlt,[getgeno(TmpSplit[x], TmpPos)  for x in postosave ])
    else :
      return (Chro, Pos, Ref, ListeAlt,[getgeno(x, TmpPos)  for x in TmpSplit[9:]])

def GetHeaderVcf(File):
    VcfRead=OpenFile(File)
    Entete=[]
    for line in VcfRead:
       if line[0]=="#" :
         Entete.append(line.lower().replace("\n",""))
       else :
        VcfRead.close()
        return Entete
    VcfRead.close()
    return Entete


def OpenFile(Fichier, Type='r'):
    try :
       if Fichier.endswith('.gz') or Fichier.endswith('.gzip'):
         LireFich=gzip.open(Fichier) 
       else :
         LireFich=open(Fichier, Type)
    except IOError:
       sys.exit("File "+ Fichier + " OuvrirFichier in mode " + Type + "can't OuvrirFichier")
    return LireFich


