GetGStat<-function(VarSet){
#Calculates G statistic for each SNP.     
    
    VarSet$GStat<-apply(VarSet,1,function(x){
        obs<-c(as.numeric(x[18]),as.numeric(x[9]),as.numeric(x[19]),as.numeric(x[10]))
        exp<-c(0.5*as.numeric(x[15]),0.5*as.numeric(x[6]),0.5*as.numeric(x[15]),0.5*as.numeric(x[6]))
        2 * sum(obs * log(obs / exp),na.rm = T)
    })
    VarSet
    
}