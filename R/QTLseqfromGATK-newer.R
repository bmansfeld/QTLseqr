library(ggplot2)
library(plotly)
library(plyr)
library(reshape2)
source("format_genomic.R")

#get and plot confidence intervals 
#Bulks of 15 individuals
confint<-read.table(file="15individuals.txt", header = T)

melt_confint<-melt(confint,measure.vars = 2:7,variable.name = "interval")

ggplot(melt_confint, aes(y=value, x=DEPTH,color=interval)) +
      geom_line() +
      scale_colour_manual(values = rep(rev(c("#F8766D","#7CAE00","#00BFC4")), each=2))

#Bulks of 30 individuals
confint30<-read.table(file="30individuals.txt", header = T)

melt_confint30<-melt(confint30,measure.vars = 2:7,variable.name = "interval")

ggplot(melt_confint30, aes(y=value, x=DEPTH,color=interval)) +
      geom_line() +
      scale_colour_manual(values = rep(rev(c("#F8766D","#7CAE00","#00BFC4")), each=2))

#Function that imports VariantsToTable results and returns a SNP set ready for Sliding window
ImportVarToTable<-function(VarTable=character(),ChromList=NULL){
      #import VariantsToTable results
      BULK<-read.table(file=VarTable, header = T, stringsAsFactors = F)
      #Remove any unwanted contigs or chromosomes
      if (!is.null(ChromList)) 
            {BULK<-BULK[BULK$CHROM %in% ChromList,]}
      #remove variant 2 caused by selecting variants after filtering = indels. might be removed eventually
      #BULK<-BULK[,1:4]
      #rename columns
      names(BULK)[5:7]<-c("AD","DP","GQ")
      #extract depth of ALT allele from AD field
      BULK$AD_REF<-as.numeric(gsub(",.*$","",x=BULK[,"AD"]))
      BULK$AD_ALT<-BULK$DP-BULK$AD_REF
      BULK$POS<-as.numeric(BULK$POS)
      BULK$SNPindex<-BULK$AD_ALT/BULK$DP
      BULK<-BULK[complete.cases(BULK),]
}



#Function that merges the BULKs based on CHROM and POS,
#fills in 0s on the SNP index for missing SNPs in the low bulk,
#then calculates Delta SNP index and minimum DEPTH
GetDeltaSNP<-function(HighBulk, LowBulk){
      tmp<-merge(x=HighBulk,y=LowBulk, all.x = T, by= c("CHROM","POS"), suffixes = c(".HIGH",".LOW"))
      #discard SNPs with less than 0.3 SNP index *AND* 7 DEPTH in HIGH bulk
      tmp<-subset(tmp, !(SNPindex.HIGH < 0.3 & DP.HIGH < 7))
      #discard SNPs with less than 0.3 SNP index *AND* 7 DEPTH in Low bulk 
      #but keep them if they have higher than 0.3 in high bulk or if they are NAs (ie. missing low bulk SNPs)
      tmp<-subset(tmp, (!(SNPindex.LOW < 0.3 & DP.LOW < 7 & SNPindex.HIGH < 0.3) | is.na(DP.LOW)))
      #set SNP index in missing SNPs in Low bulk as zero
      tmp$SNPindex.LOW[is.na(tmp$SNPindex.LOW)]<-0
      #Get delta SNP index by subtracting Low Bulk from High BULK snp index 
      tmp$deltaSNP<-tmp$SNPindex.HIGH-tmp$SNPindex.LOW
      #get minimum depth for each snp. If no snp exists in low bulk use high bulk depth
      tmp$minDepth<-apply(tmp,1,function(X) as.numeric(min(X["DP.HIGH"],X["DP.LOW"],na.rm = T)))
      tmp$maxDepth<-apply(tmp,1,function(X) as.numeric(max(X["DP.HIGH"],X["DP.LOW"],na.rm = T)))
      tmp
}


#Function that merges the BULKs based on CHROM and POS,
#fills in 0s on the SNP index for missing SNPs in the low bulk,
#then calculates Delta SNP index and minimum DEPTH
GetDeltaSNP_GQfilt<-function(HighBulk, LowBulk, minGQ=20){
      tmp<-merge(x=HighBulk,y=LowBulk, all.x = T, by= c("CHROM","POS"), suffixes = c(".HIGH",".LOW"))
      
      #discard SNPs with GQ < minGQ  in Low bulk 
      #but keep them if are NAs (ie. missing low bulk SNPs)
      tmp<-subset(tmp, (!(GQ.LOW < minGQ) | is.na(GQ.LOW)))
      
      #discard SNPs with less than minGQ GQ high bulk
      tmp<-subset(tmp, GQ.HIGH >= minGQ)
      
      #set SNP index in missing SNPs in Low bulk as zero
      tmp$SNPindex.LOW[is.na(tmp$SNPindex.LOW)]<-0
      #Get delta SNP index by subtracting Low Bulk from High BULK snp index 
      tmp$deltaSNP<-tmp$SNPindex.HIGH-tmp$SNPindex.LOW
      #get minimum depth for each snp. If no snp exists in low bulk use high bulk depth
      tmp$minDepth<-apply(tmp,1,function(X) as.numeric(min(X["DP.HIGH"],X["DP.LOW"],na.rm = T)))
      tmp$maxDepth<-apply(tmp,1,function(X) as.numeric(max(X["DP.HIGH"],X["DP.LOW"],na.rm = T)))
tmp

      }


#Function that merges the BULKs based on CHROM and POS,
#fills in SNP index for missing SNPs in the low bulk by assuming the median coverage at each pos
#then calculates Delta SNP index and minimum DEPTH
GetDeltaSNP_GQfilt_impute<-function(HighBulk, LowBulk, minGQ=20){
      tmp<-merge(x=HighBulk,y=LowBulk, all.x = T, by= c("CHROM","POS"), suffixes = c(".HIGH",".LOW"))
      
      #discard SNPs with GQ < minGQ  in Low bulk 
      #but keep them if are NAs (ie. missing low bulk SNPs)
      tmp<-subset(tmp, (!(GQ.LOW < minGQ) | is.na(GQ.LOW)))
      
      #discard SNPs with less than minGQ GQ high bulk
      tmp<-subset(tmp, GQ.HIGH >= minGQ)
      
      #set DP.LOW as median of DP in the low bulk
      tmp$DP.LOW[is.na(tmp$SNPindex.LOW)]<-median(tmp$DP.LOW, na.rm = T)
      tmp$AD_REF.LOW[is.na(tmp$SNPindex.LOW)]<-median(tmp$DP.LOW, na.rm = T)
      tmp$AD_ALT.LOW[is.na(tmp$SNPindex.LOW)]<-0
      
      #set SNP index in missing SNPs in Low bulk as zero
      tmp$SNPindex.LOW<- tmp$AD_ALT.LOW/tmp$DP.LOW
      #Get delta SNP index by subtracting Low Bulk from High BULK snp index 
      tmp$deltaSNP<-tmp$SNPindex.HIGH-tmp$SNPindex.LOW
      #get minimum depth for each snp. If no snp exists in low bulk use high bulk depth
      tmp$minDepth<-apply(tmp,1,function(X) as.numeric(min(X["DP.HIGH"],X["DP.LOW"],na.rm = T)))
      tmp$maxDepth<-apply(tmp,1,function(X) as.numeric(max(X["DP.HIGH"],X["DP.LOW"],na.rm = T)))
      tmp
      
}



GetGStat<-function(VarSet){
      VarSet$GStat<-apply(VarSet,1,function(x){
            obs<-c(as.numeric(x[18]),as.numeric(x[9]),as.numeric(x[19]),as.numeric(x[10]))
            exp<-c(0.5*as.numeric(x[15]),0.5*as.numeric(x[6]),0.5*as.numeric(x[15]),0.5*as.numeric(x[6]))
            2 * sum(obs * log(obs / exp),na.rm = T)
      })
      VarSet
      
}

#
GetGprime<-function(VarSet, WinSize = 1e4)
#     For each SNP calculates the G' statistic, a weighted average of G across neighboring SNPs.
#     G' is calculated over a window of WinSize SNPs and weighted with a tricube kernel      
{
#     Create empty dataframe to rebuild VarSet
      SW<-data.frame()
#     Calculte tricube kernel weights for WinSize
      Dvector <- abs(seq(from =-1,to=1,length=WinSize+1))
      KNum <- (1- Dvector^3)^3
      KDen <- sum(KNum)
      K<-KNum/KDen
      
#     Calculate G' for each SNP within each chromosome
      for (x in levels(as.factor(VarSet$CHROM))) {
            
            chr<-as.data.frame(subset(VarSet, CHROM == x))
            
            chr$Gprime<-rollapply(chr$GStat, width = WinSize+1, align = "center",fill=NA, 
                                  partial=T,FUN=function(x) {
                                        sum(x*K[length(K)-length(x)+1:length(K)], na.rm = T)
                                  })
            
            chr$pval<- 1-pnorm((log(chr$Gprime)-mean(chr$Gprime)/sd(chr$Gprime)))
            
            SW<-rbind(SW,chr)
      }
      SW
      
}








      
##Function that performs sliding window analysis, 
#returns df of positions steps and their respective mean SNPindex, mean depth, and count of SNPs
#also gets the confidence intervals from the supplied confintTable 
#and filters out windows with less than minSNPs SNPs.
GetVARSlidingWindow<-function(VarSet, WinSize = 1e6, StepSize=1e4, confintTable, minSNPs=10){
      SW<-data.frame(POS=numeric(0),HighSNPindex=numeric(0), LowSNPindex=numeric(0), deltaSNPindex=numeric(0),CHROM=integer(0),DEPTH=integer(0),nSNPs=integer(0))
      for (x in paste0(rep("Chr",7),1:7)) {
            
            chr<-as.data.frame(subset(VarSet, CHROM == x))
            
            #create sliding window step bins
            bin <- seq(1,max(chr[,"POS"]),StepSize)
            
            SWHighindex<-sapply(1:(length(bin)),
                                 function(i)
                                       mean(chr[chr$POS >= bin[i]-(0.5*WinSize) &
                                                      chr$POS < bin[i]+(0.5*WinSize), "SNPindex.HIGH"], na.rm=T))
            SWLowindex<-sapply(1:(length(bin)),
                                function(i)
                                      mean(chr[chr$POS >= bin[i]-(0.5*WinSize) &
                                                     chr$POS < bin[i]+(0.5*WinSize), "SNPindex.LOW"], na.rm=T))
            
            SWDeltaindex<-sapply(1:(length(bin)),
                            function(i)
                                  mean(chr[chr$POS >= bin[i]-(0.5*WinSize) &
                                                 chr$POS < bin[i]+(0.5*WinSize), "deltaSNP"], na.rm=T))
            #get mean SNPindex per step
            # SWindex2<-sapply(1:(length(bin)),
            #                 function(i)
            #                       mean(chr[chr$POS>=bin[i] &
            #                                      chr$POS < bin[i]+ WinSize, "SNPindex"], na.rm=T))
            #get mean depth of window step

            SWdepth<-sapply(1:(length(bin)),
                            function(i)
                                  mean(chr[chr$POS>=bin[i]-(0.5*WinSize) &
                                                 chr$POS < bin[i]+(0.5*WinSize),"minDepth"], na.rm=T))
            # SWdepth2<-sapply(1:(length(bin)),
            #                 function(i)
            #                       mean(chr[chr$POS>=bin[i] &
            #                                      chr$POS < bin[i]+WinSize, "DP"], na.rm=T))
            
            #get count of SNPs in window step
            count<-sapply(1:(length(bin)),
                          function(i)
                                nrow(chr[chr$POS>=bin[i]-(0.5*WinSize)  &
                                               chr$POS < bin[i]+(0.5*WinSize),]))
            
            chrom <- data.frame(POS=bin, HighSNPindex=SWHighindex, LowSNPindex=SWLowindex, deltaSNPindex=SWDeltaindex, CHROM=x, DEPTH=SWdepth, nSNPs=count)  
            SW<-rbind(SW,chrom)
      }
      SW$conf99H<-ifelse(SW$DEPTH > 300 ,confintTable[300,"P_H_99"],confintTable[floor(SW$DEPTH),"P_H_99"])
      SW$conf99L<-ifelse(SW$DEPTH > 300 ,confintTable[300,"P_L_99"],confintTable[floor(SW$DEPTH),"P_L_99"])
      SW$conf95H<-ifelse(SW$DEPTH > 300 ,confintTable[300,"P_H_95"],confintTable[floor(SW$DEPTH),"P_H_95"])
      SW$conf95L<-ifelse(SW$DEPTH > 300 ,confintTable[300,"P_L_95"],confintTable[floor(SW$DEPTH),"P_L_95"])
      SW$conf90H<-ifelse(SW$DEPTH > 300 ,confintTable[300,"P_H_90"],confintTable[floor(SW$DEPTH),"P_H_90"])
      SW$conf90L<-ifelse(SW$DEPTH > 300 ,confintTable[300,"P_L_90"],confintTable[floor(SW$DEPTH),"P_L_90"])
      
      #skip windows with less than minSNPs SNPs
      SW<-subset(SW, nSNPs>= minSNPs)
      SW

}

#plot SNP index graphs for high and low bulks

PlotSNPindex<-function(SW_deltaSNPTable_with_CI,deltaSNPTable, PlotName, PlotPoints=F, HighBulkLab = "High Bulk SNPindex", LowBulkLab = "Low Bulk SNPindex"){
      #Set Labels for faceting
      bulklabels=c(HighSNPindex = HighBulkLab, LowSNPindex = LowBulkLab)
      #melt SNP index df for ploting 
      SW_SNPindex<-melt(SW_deltaSNPTable_with_CI, measure.vars = c("HighSNPindex","LowSNPindex"),id.vars = c("CHROM","POS"), 
                        value.name = "SNPindex",variable.name = "BULK")
      
      SNPindex<-melt(deltaSNPTable, measure.vars = c("SNPindex.HIGH","SNPindex.LOW"),id.vars = c("CHROM","POS"), 
                     value.name = "SNPindex",variable.name = "BULK")
      levels(SNPindex$BULK)<-c("HighSNPindex","LowSNPindex")
      
      p<-ggplot(SW_SNPindex) +
      {if (PlotPoints) geom_point(data=subset(SNPindex, CHROM %in% paste0(rep("Chr",7),1:7)), aes(x=POS,y=SNPindex),
                                  color="lightblue", alpha=0.3)} + 
            geom_line(aes(x=POS, y=SNPindex), size=1) + 
            facet_grid(BULK ~ CHROM, scales = "free_x",labeller=labeller(BULK = bulklabels)) +
            labs(title=PlotName) +
            scale_x_continuous(labels=format_genomic()) +
            theme_bw()
      p
}

PlotSNPindex2<-function(SW_deltaSNPTable_with_CI,deltaSNPTable, PlotName, PlotPoints=F){
      #melt SNP index df for ploting 
      SW_SNPindex<-melt(SW_deltaSNPTable_with_CI, measure.vars = c("HighSNPindex","LowSNPindex"),id.vars = c("CHROM","POS"), 
                        value.name = "SNPindex",variable.name = "BULK")
      
      SNPindex<-melt(deltaSNPTable, measure.vars = c("SNPindex.HIGH","SNPindex.LOW"),id.vars = c("CHROM","POS","minDepth"), 
                     value.name = "SNPindex",variable.name = "BULK")
      levels(SNPindex$BULK)<-c("HighSNPindex","LowSNPindex","minDepth")
      
      p<-ggplot(SW_SNPindex) +
      {if (PlotPoints) geom_point(data=subset(SNPindex, CHROM %in% paste0(rep("Chr",7),1:7)), aes(x=POS,y=SNPindex, color=minDepth),
                                  alpha=0.3)} + 
            geom_line(aes(x=POS, y=SNPindex), size=1) + 
            facet_grid(BULK ~ CHROM, scales = "free_x") +
            labs(title=PlotName) +
            scale_x_continuous(labels=format_genomic()) +
            theme_bw()
      p
}

# #plot delta SNP index plots with confidence intervals
# PlotDeltaSNP<-function(deltaSNPTable, PlotName){
#       deltaSNP_CI<-melt(deltaSNPTable,id.vars = 1:2,measure.vars = 12:17,variable.name = "interval")
#       p<-ggplot(deltaSNPTable) + 
#             geom_line(data =deltaSNP_CI, aes(x=POS, y=value, color=interval)) +
#             geom_point(aes(x=POS, y=deltaSNPindex),size=1,alpha=0.8) +
#             facet_grid(~ CHROM, scales = "free_x") +
#             scale_colour_manual(values = rep(c("#F8766D","#7CAE00","#00BFC4"), each=2)) +
#             labs(x="Position", y=expression(paste(Delta,"(SNP-index)")),
#                   title=PlotName) +
#             theme_bw() +
#             theme(panel.spacing.x =unit(0, "lines"))
#       p
# }

PlotDeltaSNP<-function(SW_deltaSNPTable_with_CI,deltaSNPTable, PlotName, PlotPoints=F){
      deltaSNP_CI<-melt(SW_deltaSNPTable_with_CI,id.vars = c("CHROM","POS"),
                        measure.vars = c("conf99H","conf99L","conf95H","conf95L","conf90H","conf90L"),
                        variable.name = "interval")
      ggplot(SW_deltaSNPTable_with_CI) + 
      {if (PlotPoints) geom_point(data=subset(deltaSNPTable, CHROM %in% paste0(rep("Chr",7),1:7)), aes(x=POS,y=deltaSNP),
                                  color="lightblue", alpha=0.3)} + 
            geom_line(data =deltaSNP_CI, aes(x=POS, y=value, color=interval)) +
            scale_colour_manual(values = rep(c("#F8766D","#7CAE00","#00BFC4"), each=2)) +
            geom_line(aes(x=POS, y=deltaSNPindex), size=1) + 
            facet_grid( ~ CHROM, scales = "free_x") +
            labs(x="Position", y=expression(paste(Delta,"(SNP-index)")),
                 title=PlotName) +
            ylim(-0.6,0.6) +
            scale_x_continuous(labels=format_genomic()) +
            geom_abline(aes(intercept = 0, slope = 0))+
            theme_bw()
}

PlotDeltaSNPChrom<-function(SW_deltaSNPTable_with_CI,deltaSNPTable, PlotName, Chr="Chr1",PlotPoints=F){
      Chrom<-subset(SW_deltaSNPTable_with_CI, CHROM==Chr)
      deltaSNP_CI<-melt(Chrom,id.vars = c("CHROM","POS"),
                        measure.vars = c("conf99H","conf99L","conf95H","conf95L","conf90H","conf90L"),
                        variable.name = "interval")
      ggplot(Chrom) + 
      {if (PlotPoints) geom_point(data=subset(deltaSNPTable, CHROM %in% Chr), aes(x=POS,y=deltaSNP),
                                  color="lightblue", alpha=0.3)} + 
            geom_smooth(data =deltaSNP_CI, aes(x=POS, y=value, color=interval)) +
            scale_colour_manual(values = rep(c("#F8766D","#7CAE00","#00BFC4"), each=2)) +
            geom_line(aes(x=POS, y=deltaSNPindex), size=1) + 
            facet_grid( ~ CHROM, scales = "free_x") +
            labs(x="Position", y=expression(paste(Delta,"(SNP-index)")),
                 title=PlotName) +
            scale_x_continuous(labels=format_genomic()) +
            theme_bw()
}


#Function that returns positions of sliding windows that pass a user 
#defined confidence interval threshold
sigWindows<-function(SW_deltaSNPTable_with_CI, CIthreshold="90") {
      col<-if(CIthreshold == "99") c(8:9) else
                  if(CIthreshold == "95") c(10:11) else
                         if(CIthreshold == "90") c(12:13) else
                                stop("Confidence interval thresholds should be strings with values of \"90\", \"95\" or \"99\"")
      whichWind<-which(SW_deltaSNPTable_with_CI$deltaSNPindex > SW_deltaSNPTable_with_CI[col[1]] | SW_deltaSNPTable_with_CI$deltaSNPindex < SW_deltaSNPTable_with_CI[col[2]])
      SW_deltaSNPTable_with_CI[whichWind,]
}


###############################################################################
#########################   VLASPIK VS. GY   ##################################
###############################################################################

#Get Data
R_BULK_filt<-ImportVarToTable(VarTable="VGH1_recal_1_SNPS.table")
S_BULK_filt<-ImportVarToTable(VarTable="VGL1_recal_1_SNPS.table")

#get Delta SNP indices
VGH1vsVGL1<-GetDeltaSNP(R_BULK_filt,S_BULK_filt)

#get sliding window values for both bulks
SW_VGH1vsVGL1<-GetVARSlidingWindow(VGH1vsVGL1, WinSize = 1e6, StepSize = 1e4, confint)

#plots
png(filename = "VGH1vsVGL1_SNPindex.png", width=1200, height=600, res = 100)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_VGH1vsVGL1, deltaSNPTable = VGH1vsVGL1, PlotName = "VGH1vsVGL1", PlotPoints = F)
dev.off()

png(filename = "VGH1vsVGL1_DELTASNP2.png", width=1200, height=600, res = 100)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_VGH1vsVGL1, deltaSNPTable = VGH1vsVGL1, PlotName = "VGH1vsVGL1", PlotPoints = T)
dev.off()
# #SW_Parent<-GetVARSlidingWindow(Parent, WinSize = 1e6, StepSize = 1e4)

PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_VGH1vsVGL1, deltaSNPTable = VGH1vsVGL1, PlotName = "VGH1vsVGL1", Chr = "Chr5", PlotPoints = T)


#############################################################################
R_BULK_filt<-ImportVarToTable(VarTable="VGH2_recal_1_SNPS.table")
S_BULK_filt<-ImportVarToTable(VarTable="VGL2_recal_1_SNPS.table")

#get Delta SNP indices
VGH2vsVGL2<-GetDeltaSNP(R_BULK_filt,S_BULK_filt)

#get sliding window values for both bulks
SW_VGH2vsVGL2<-GetVARSlidingWindow(VGH2vsVGL2, WinSize = 1e6, StepSize = 1e4, confint)

#plots
png(filename = "VGH2vsVGL2_SNPindex.png", width=1200, height=600, res = 100)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_VGH2vsVGL2, deltaSNPTable = VGH2vsVGL2, PlotName = "VGH2vsVGL2", PlotPoints = F)
dev.off()
png(filename = "VGH2vsVGL2_DELTASNP.png", width=1200, height=600, res = 100)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_VGH2vsVGL2, deltaSNPTable = VGH2vsVGL2, PlotName = "VGH2vsVGL2", PlotPoints = F)
dev.off()

PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_VGH2vsVGL2, deltaSNPTable = VGH2vsVGL2, PlotName = "VGH2vsVGL2", Chr = "Chr5", PlotPoints = T)


###############################################################################
#########################   POINSETTE VS. GY   ################################
###############################################################################

R_BULK_filt<-ImportVarToTable(VarTable="PGH1_recal_1_SNPS.table")
S_BULK_filt<-ImportVarToTable(VarTable="PGL1_recal_1_SNPS.table")

#get Delta SNP indices
PGH1vsPGL1<-GetDeltaSNP(R_BULK_filt,S_BULK_filt)

#get sliding window values for both bulks
SW_PGH1vsPGL1<-GetVARSlidingWindow(PGH1vsPGL1, WinSize = 1e6, StepSize = 1e4, confint)

#plots
png(filename = "PGH1vsPGL1_SNPindex.png", width=1200, height=600, res = 100)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_PGH1vsPGL1, deltaSNPTable = PGH1vsPGL1, PlotName = "PGH1vsPGL1", PlotPoints = T)
dev.off()
png(filename = "PGH1vsPGL1_DELTASNP.png", width=1200, height=600, res = 100)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_PGH1vsPGL1, deltaSNPTable = PGH1vsPGL1, PlotName = "PGH1vsPGL1", PlotPoints = F)
dev.off()

PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_PGH1vsPGL1, deltaSNPTable = PGH1vsPGL1, PlotName = "PGH1vsPGL1", Chr = "Chr5", PlotPoints = T)




# #get sliding window values for both bulks
# SW_H_BULK<-GetVARSlidingWindow(H_BULK_filt, WinSize = 1e6, StepSize = 1e4)
# SW_L_BULK<-GetVARSlidingWindow(L_BULK_filt, WinSize = 1e6, StepSize = 1e4)
# 
# #get Delta SNP indices
# PGH1vsPGL1<-GetDeltaSNP(SW_H_BULK,SW_L_BULK,confint)
# #plot
# PlotSNPindex(PGH1vsPGL1, "PGH1vsPGL1")
# PlotDeltaSNP(PGH1vsPGL1, "PGH1vsPGL1")

ggplot(subset(tmp, CHROM %in% paste0(rep("Chr",7),1:7)), aes(x=POS,y=deltaSNP)) + geom_point() +
      facet_grid(~ CHROM, scales = "free_x") 



###############################################################################
###############################################################################
###############################   OTHERs   ####################################
###############################################################################
###############################################################################

############Test with GQ filtering
#Get Data
R_BULK_filt<-ImportVarToTable(VarTable="VGH1_recal_1_SNPS.table")
S_BULK_filt<-ImportVarToTable(VarTable="VGL1_recal_1_SNPS.table")

#get Delta SNP indices
VGH1vsVGL1_GQfilt<-GetDeltaSNP_GQfilt(R_BULK_filt,S_BULK_filt, minGQ = 98)

#get sliding window values for both bulks
SW_VGH1vsVGL1_GQfilt<-GetVARSlidingWindow(VGH1vsVGL1_GQfilt, WinSize = 2e6, StepSize = 1e4, confint)

#plots
png(filename = "VGH1vsVGL1_SNPindex_GQfilt.png", width=2400, height=1200, res = 200)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_VGH1vsVGL1_GQfilt, deltaSNPTable = VGH1vsVGL1_GQfilt, PlotName = "VGH1vsVGL1_GQfilt", PlotPoints = T)
dev.off()

png(filename = "VGH1vsVGL1_DELTASNP_GQfilt.png", width=2400, height=1200, res = 200)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_VGH1vsVGL1_GQfilt, deltaSNPTable = VGH1vsVGL1_GQfilt, PlotName = "VGH1vsVGL1_GQfilt", PlotPoints = F)
dev.off()
# #SW_Parent<-GetVARSlidingWindow(Parent, WinSize = 1e6, StepSize = 1e4)

#p<-PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_VGH1vsVGL1_GQfilt, deltaSNPTable = VGH1vsVGL1_GQfilt, PlotName = "VGH1vsVGL1", Chr = "Chr6", PlotPoints = T)
#ggplotly(p)

################################  SEASON 2 ######################################

R_BULK_filt<-ImportVarToTable(VarTable="VGH2_recal_1_SNPS.table")
S_BULK_filt<-ImportVarToTable(VarTable="VGL2_recal_1_SNPS.table")

#get Delta SNP indices
VGH2vsVGL2_GQfilt<-GetDeltaSNP_GQfilt(R_BULK_filt,S_BULK_filt, minGQ = 98)

#get sliding window values for both bulks
SW_VGH2vsVGL2_GQfilt<-GetVARSlidingWindow(VGH2vsVGL2_GQfilt, WinSize = 2e6, StepSize = 1e4, confint)

#plots
png(filename = "VGH2vsVGL2_SNPindex_GQfilt.png", width=2400, height=1200, res = 200)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_VGH2vsVGL2_GQfilt, deltaSNPTable = VGH2vsVGL2_GQfilt, PlotName = "VGH2vsVGL2_GQfilt", PlotPoints = T)
dev.off()
png(filename = "VGH2vsVGL2_DELTASNP_GQfilt.png", width=2400, height=1200, res = 200)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_VGH2vsVGL2_GQfilt, deltaSNPTable = VGH2vsVGL2_GQfilt, PlotName = "VGH2vsVGL2_GQfilt", PlotPoints = F)
dev.off()

#PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_VGH2vsVGL2, deltaSNPTable = VGH2vsVGL2, PlotName = "VGH2vsVGL2", Chr = "Chr6", PlotPoints = T)


#############################################################################
#########################   POINSETTE VS. GY   ##############################
#############################################################################

R_BULK_filt<-ImportVarToTable(VarTable="PGH1_recal_1_SNPS.table")
S_BULK_filt<-ImportVarToTable(VarTable="PGL1_recal_1_SNPS.table")

#get Delta SNP indices
PGH1vsPGL1_GQfilt<-GetDeltaSNP_GQfilt(R_BULK_filt,S_BULK_filt,minGQ = 98)

#get sliding window values for both bulks
SW_PGH1vsPGL1_GQfilt<-GetVARSlidingWindow(PGH1vsPGL1_GQfilt, WinSize = 2e6, StepSize = 1e4, confint)

#plots
png(filename = "PGH1vsPGL1_SNPindex_GQfilt.png", width=2400, height=1200, res = 200)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_PGH1vsPGL1_GQfilt, deltaSNPTable = PGH1vsPGL1_GQfilt, PlotName = "PGH1vsPGL1_GQfilt", PlotPoints = T)
dev.off()
png(filename = "PGH1vsPGL1_DELTASNP_GQfilt.png", width=2400, height=1200, res = 200)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_PGH1vsPGL1_GQfilt, deltaSNPTable = PGH1vsPGL1_GQfilt, PlotName = "PGH1vsPGL1_GQfilt", PlotPoints = F)
dev.off()

PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_PGH1vsPGL1, deltaSNPTable = PGH1vsPGL1, PlotName = "PGH1vsPGL1", Chr = "Chr5", PlotPoints = T)

################################  SEASON 2 ######################################


R_BULK_filt<-ImportVarToTable(VarTable="PGH2_recal_1_SNPS.table")
S_BULK_filt<-ImportVarToTable(VarTable="PGL2_recal_1_SNPS.table")

#get Delta SNP indices
PGH2vsPGL2_GQfilt<-GetDeltaSNP_GQfilt(R_BULK_filt,S_BULK_filt,minGQ = 40)

#get sliding window values for both bulks
SW_PGH2vsPGL2_GQfilt<-GetVARSlidingWindow(PGH2vsPGL2_GQfilt, WinSize = 2e6, StepSize = 1e4, confint)

#plots
png(filename = "PGH2vsPGL2_SNPindex_GQfilt.png", width=2400, height=1200, res = 200)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_PGH2vsPGL2_GQfilt, deltaSNPTable = PGH2vsPGL2_GQfilt, PlotName = "PGH2vsPGL2_GQfilt", PlotPoints = T)
dev.off()
png(filename = "PGH2vsPGL2_DELTASNP_GQfilt.png", width=2400, height=1200, res = 200)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_PGH2vsPGL2_GQfilt, deltaSNPTable = PGH2vsPGL2_GQfilt, PlotName = "PGH2vsPGL2_GQfilt", PlotPoints = F)
dev.off()

PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_PGH2vsPGL2, deltaSNPTable = PGH2vsPGL2, PlotName = "PGH2vsPGL2", Chr = "Chr5", PlotPoints = T)




################################  MERGED SEASONS ######################################

Chromosomes<-paste0(rep("Chr",7),1:7)

R_BULK_filt<-ImportVarToTable(VarTable="PGH_mergedSE_SNPs.table", ChromList = Chromosomes)
S_BULK_filt<-ImportVarToTable(VarTable="PGL_mergedSE_SNPs.table", ChromList = Chromosomes)

#get Delta SNP indices
PGHvsPGL_mergedSE_GQfilt<-GetDeltaSNP_GQfilt(R_BULK_filt,S_BULK_filt,minGQ = 30)

#subset by min depth <= n
PGHvsPGL_mergedSE_GQfilt <- subset(PGHvsPGL_mergedSE_GQfilt, minDepth <= 2*mean(PGHvsPGL_mergedSE_GQfilt$minDepth, na.rm = T))

#get sliding window values for both bulks
SW_PGHvsPGL_mergedSE_GQfilt<-GetVARSlidingWindow(PGHvsPGL_mergedSE_GQfilt, WinSize = 2e6, StepSize = 1e4, confint30,minSNPs = 10)

#plots
png(filename = "PGHvsPGL_mergedSE_SNPindex_GQfilt.png", width=3000, height=1000, res = 200)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_PGHvsPGL_mergedSE_GQfilt, deltaSNPTable = PGHvsPGL_mergedSE_GQfilt, PlotName = "PGHvsPGL_mergedSE_GQfilt", PlotPoints = T,HighBulkLab = "Resistant Bulk", LowBulkLab = "Susceptible Bulk")
dev.off()

# png(filename = "PGHvsPGL_mergedSE_SNPindex_GQfilt2.png", width=2400, height=1200, res = 200)
# PlotSNPindex2(SW_deltaSNPTable_with_CI = SW_PGHvsPGL_mergedSE_GQfilt, deltaSNPTable = PGHvsPGL_mergedSE_GQfilt, PlotName = "PGHvsPGL_mergedSE_GQfilt", PlotPoints = F)
# dev.off()

png(filename = "PGHvsPGL_mergedSE_DELTASNP_GQfilt.png", width=3000, height=1000, res = 200)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_PGHvsPGL_mergedSE_GQfilt, deltaSNPTable = PGHvsPGL_mergedSE_GQfilt, PlotName = "PGHvsPGL_mergedSE_GQfilt", PlotPoints = T)
dev.off()

PlotDeltaSNPChrom(SW_deltaSNPTable_with_CI = SW_PGHvsPGL_mergedSE_GQfilt, deltaSNPTable = PGHvsPGL_mergedSE_GQfilt, PlotName = "PGHvsPGL", Chr = "Chr3", PlotPoints = T)

sigW<-sigWindows(SW_deltaSNPTable_with_CI = SW_PGHvsPGL_mergedSE_GQfilt, CIthreshold = "90")

####Vlaspik vs GY merged SE ##############
R_BULK_filt<-ImportVarToTable(VarTable="VGH_mergedSE_SNP.table")
S_BULK_filt<-ImportVarToTable(VarTable="VGL_mergedSE_SNP.table")

#get Delta SNP indices
VGHvsVGL_mergedSE_GQfilt<-GetDeltaSNP_GQfilt(R_BULK_filt,S_BULK_filt,minGQ = 20)

#subset by min depth <= n
VGHvsVGL_mergedSE_GQfilt<-subset(VGHvsVGL_mergedSE_GQfilt, minDepth <= 60)

#get sliding window values for both bulks
SW_VGHvsVGL_mergedSE_GQfilt<-GetVARSlidingWindow(VGHvsVGL_mergedSE_GQfilt, WinSize = 2e6, StepSize = 1e4, confint30)

#plots
png(filename = "VGHvsVGL_mergedSE_SNPindex_GQfilt.png", width=2400, height=1200, res = 200)
PlotSNPindex(SW_deltaSNPTable_with_CI = SW_VGHvsVGL_mergedSE_GQfilt, deltaSNPTable = VGHvsVGL_mergedSE_GQfilt, PlotName = "VGHvsVGL_mergedSE_GQfilt", PlotPoints = F)
dev.off()
png(filename = "VGHvsVGL_mergedSE_DELTASNP_GQfilt.png", width=2400, height=1200, res = 200)
PlotDeltaSNP(SW_deltaSNPTable_with_CI = SW_VGHvsVGL_mergedSE_GQfilt, deltaSNPTable = VGHvsVGL_mergedSE_GQfilt, PlotName = "VGHvsVGL_mergedSE_GQfilt", PlotPoints = F)
      
dev.off()



###READ Parents SNPS###
poinsett <- ImportVarToTable(VarTable = "PG_recal_1_SNPS.table")
vlaspik <- ImportVarToTable(VarTable = "V_recal_1_SNPS.table")
gy <- ImportVarToTable(VarTable = "Gy_recal_2_SNPS.table")



# ################OTHER OPTIONS IN CASE OF MIS-LABELING#########################
# ##############################################################################
# H_BULK_filt<-ImportVarToTable(VarTable="VGH1_SNPs_results.table")
# L_BULK_filt<-ImportVarToTable(VarTable="VGL2_results.table")
# 
# #get sliding window values for both bulks
# SW_H_BULK<-GetVARSlidingWindow(H_BULK_filt, WinSize = 1e6, StepSize = 1e4)
# SW_L_BULK<-GetVARSlidingWindow(L_BULK_filt, WinSize = 1e6, StepSize = 1e4)
# 
# #get Delta SNP indices
# VGH1vsVGL2<-GetDeltaSNP(SW_H_BULK,SW_L_BULK,confint)
# 
# #############################################################################
# 
# ggplot(SW_Parent) + 
#       geom_point(aes(x=POS, y=SNPindex)) + 
#       +
#       theme_bw()
# 
# ggplot(subset(SNP, CHROM=="Chr7")) + 
#       geom_point(aes(x=POS, y=SNPindex, clor=CHROM)) + 
#       facet_grid(BULK ~ CHROM, scales = "free_x") +
#       theme(panel.spacing.x =unit(0, "lines"))
# 
# facet_grid(~ CHROM, scales = "free_x") 
# #plot SetChr
# SetChr<-"Chr6"
# p<-ggplot(subset(deltaSNP, CHROM==SetChr)) + 
#       geom_line(data =subset(deltaSNP_CI, CHROM==SetChr), aes(x=POS, y=value, color=interval)) +
#       scale_colour_manual(values = rep(c("#F8766D","#7CAE00","#00BFC4"), each=2)) +
#       geom_point(aes(x=POS, y=deltaSNPindex),size=1,alpha=0.8) +
#       facet_grid(~ CHROM, scales = "free_x") +
#       labs(x="Position", y=expression(paste(Delta,"(SNP-index)"))) +
#       theme_bw() +
#       theme(panel.spacing.x =unit(0, "lines"))
# p
# #generate interactive plot
# ggplotly(p)
# 
# #plot Chr
# SetChr<-"Chr6"
# p<-ggplot(subset(deltaSNP, CHROM==SetChr)) + 
#       geom_point(aes(x=POS, y=deltaSNPindex,color=CHROM)) +
#       geom_line(aes(x=POS, y=conf99H),color="red", ) +
#       geom_line(aes(x=POS, y=conf95H), color="blue") +
#       geom_line(aes(x=POS, y=conf99L),color="red") +
#       geom_line(aes(x=POS, y=conf95L), color="blue") +
#       facet_grid(~ CHROM, scales = "free_x") +
#       theme(panel.spacing.x =unit(0, "lines"))
# p
# #generate interactive plot
# ggplotly(p)
