GetDeltaSNP_GQfilt_impute<-function(HighBulk, LowBulk, minGQ=20){

#Takes SNPs from two bulks and merges them to one Variant Set.
#Filtering is done by Genotype Quality using the minGQ argument.
#For SNPs that appear in the high bulk but are missing in the low bulk,
#the function assumes that all reads at that site had REF alleles and
#uses the median read depth for calculating SNP index

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
