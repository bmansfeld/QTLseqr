GetGStat <- function(VarSet) {
  VarSet$GStat <- apply(VarSet[,5:ncol(VarSet)], 1, function(x) {
    obs <-
      c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
    exp <-
      c(0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]), 0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]))
    2 * sum(obs * log(obs / exp), na.rm = T)
  })
  VarSet
  
}

GetGStat2 <- function(VarSet) {
  VarSet$GStat2 <- apply(VarSet[,5:ncol(VarSet)], 1, function(x) {
    obs <-
      c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
    exp <- rep((((x["AD_REF.LOW"]) + (x["AD_ALT.LOW"])) * ((x["AD_REF.HIGH"]) + (x["AD_ALT.HIGH"]))) / sum(obs), 4)
    2 * sum(obs * log(obs / exp), na.rm = T)
  })
  VarSet
  
}