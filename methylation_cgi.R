#Analyzing a CGI sums file from beginning to end 
#REQUIRED LIBRARIES
require(optparse)
#Input example:
  #chr1 500 600 78  98
  #assuming the file does not have a header
cgi_sums.df<-read.table("h1_rrbs_cgi_sum.bed", header=FALSE,sep="\t",row.names=NULL)
  #Add column names 
colnames(cgi_sums.df)<-c("Chromosome","Start","End",
                         "Methylated","Unmethylated")
  #calculate the columns Length, Total,Weighted
  
cgi_sums.df$Length<-cgi_sums.df$End-cgi_sums.df$Start
cgi_sums.df$Total<-cgi_sums.df$Methylated+cgi_sums.df$Unmethylated
cgi_sums.df$Weighted<-cgi_sums.df$Methylated/cgi_sums.df$Total

#Apply the methylation caller function
  #THE METHYLATION CALLER FUNCTION
  
  #--------------------------
  #METHYLATION CALLER FUNCTION
  # -input: weighted value
  # -methylation and unmethylation threshold
  #all function names are camelBack
  
  methylationCaller<-function(weighted, meth.tsh, unmeth.tsh){
    if(weighted>=meth.tsh){
      return("methylated")
    }
    else if(weighted<=unmeth.tsh){
      return("unmethylated")
    }
    else{
      return("partial")
    }
  }
  
methylation_call.vector<-sapply(cgi_sums.df$Weighted,
                                  methylationCaller, meth.tsh=0.8, unmeth.tsh=0.2)
  #append methylation_call.vector to the main dataframe
cgi_sums.df$Status<-methylation_call.vector

cgi_sums_bed.df<-cgi_sums.df[which(cgi_sums.df$Status=='methylated' | cgi_sums.df$Status=='unmethylated'), 
                             c("Chromosome","Start","End", "Status")]
                                   
                                   
                                
  
#write the BED like file 
write.table(cgi_sums_bed.df,"cgi_bed.bed", quote=FALSE, sep="\t", row.names=FALSE,
            col.names=FALSE)
  
#for kmerSVM, produce two separate files 
cgi_sums_methylated<-cgi_sums.df[cgi_sums.df$Status %in% c('methylated'),c(1,2,3)]
cgi_sums_unmethylated<-cgi_sums.df[cgi_sums.df$Status %in% c('unmethylated','unmethylated'),c(1,2,3)]

#write to BED like file
write.table(cgi_sums_methylated,"cgi_methylated.bed", quote=FALSE, sep="\t", row.names=FALSE,
            col.names=FALSE)
write.table(cgi_sums_unmethylated,"cgi_unmethylated.bed", quote=FALSE, sep="\t", row.names=FALSE,
            col.names=FALSE)

  