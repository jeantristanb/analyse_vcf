from scipy import stats

class Binom :
    def GetMatricePvalue(self,Max,p, Loi) :
        tab=[]
        for NbN in range(0,Max+1) :
           TabNb=[]
           for Nb in range(0,NbN+1) :
               TabNb.append(Loi(Nb,NbN, p))
           tab.append(TabNb)
        return tab
    def Test(self, N,NTot):
       if NTot>self.Max:
          return self.Loi(N,NTot, self.p)
       else :
          return self.tab[NTot][N]
    def __init__(self, Max, p, Type='pvalue'):
      self.p=p
      self.Max=Max
      if Type=='pvalue':
           self.Loi=stats.binom_test
      elif Type=='probLog':
           self.Loi=stats.binom.logpmf
      else :
           raise
      self.tab=self.GetMatricePvalue(Max,p,self.Loi)


