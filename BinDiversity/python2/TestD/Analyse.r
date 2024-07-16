
a=read.table("TestD_10_1000_0.5.out")


a2<-a[,-c(2,3,4,5)]
apply(a2, 2, mean, na.rm=T)
apply(a2, 2, mean, na.rm=T)
boxplot(a2)
