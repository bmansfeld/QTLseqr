ImportVarToTable<-function(VarTable=character(),ChromList=NULL){

#Function that imports VariantsToTable (GATK) results and returns a SNP set ready for Sliding window
#Use the ChromList paramater to subset your data to include only chromosomes in this list


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
