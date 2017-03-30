GetGStat <- function(SNPset) {
  SNPset$GStat <- apply(SNPset[,5:ncol(SNPset)], 1, function(x) {
    obs <-
      c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
    exp <-
      c(0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]), 0.5 * (x["DP.LOW"]), 0.5 * (x["DP.HIGH"]))
    2 * sum(obs * log(obs / exp), na.rm = T)
  })
  SNPset

}

GetGStat2 <- function(SNPset) {
  SNPset$GStat2 <- apply(SNPset[,5:ncol(SNPset)], 1, function(x) {
    obs <-
      c((x["AD_REF.LOW"]), (x["AD_REF.HIGH"]), (x["AD_ALT.LOW"]), (x["AD_ALT.HIGH"]))
    exp <- rep((((x["AD_REF.LOW"]) + (x["AD_ALT.LOW"])) * ((x["AD_REF.HIGH"]) + (x["AD_ALT.HIGH"]))) / sum(obs), 4)
    2 * sum(obs * log(obs / exp), na.rm = T)
  })
  SNPset

}
