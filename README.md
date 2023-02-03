Insplico Simulation
================
Luis P Iniguez

### **Generating GTF (bash):**

``` bash
wget https://ftp.ensembl.org/pub/release-85/gtf/homo_sapiens/Homo_sapiens.GRCh38.85.gtf.gz  && gunzip Homo_sapiens.GRCh38.85.gtf.gz

perl -e 'open($fh,"Homo_sapiens.GRCh38.85.gtf");
while(<$fh>){
    chomp;
    if(substr($_,0,1) eq "#"){print $_."\n";next;}
    @fs=split("\t",$_);
  if($fs[0] =~ /^\d+/ || $fs[0] eq "X" || $fs[0] eq "Y" || $fs[0] eq "MT"){$fs[0]="chr".$fs[0];}
    print join("\t",@fs)."\n";
}' > Hsa38.gtf
rm Homo_sapiens.GRCh38.85.gtf
```

### **File parsing (bash):**

``` bash

cat Hsa38.gtf |\
perl -F'\t' -ane 'if(@F<5 || /transcript_support_level \"1\"/){print $_}' >\ Hsa38_only_transcript_support_level_1.gtf


grep '^chr.+\texon\t' -P Hsa38_only_transcript_support_level_1.gtf |\
awk '{ print $10,$14,$20;}'| sed 's/[\",\;]//g' | sort | uniq -c |\
awk '{for (i=1;i<=NF;i++){if($i ~ /.+/){printf "\t%s", $i;}}printf "\n"}' |\
sed 's/^\t//g' | awk '{print $3,$4,$2,$1;}' OFS="\t" > Transcript_Exon_Count.tab


mkdir -p chroms
FAPATH="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.22_GRCh38.p7/GCA_000001405.22_GRCh38.p7_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"
wget ${FAPATH}chr1.fna.gz && gunzip chr1.fna.gz && mv chr1.fna chroms/chr1.fa
wget ${FAPATH}chr2.fna.gz && gunzip chr2.fna.gz && mv chr2.fna chroms/chr2.fa
wget ${FAPATH}chr3.fna.gz && gunzip chr3.fna.gz && mv chr3.fna chroms/chr3.fa
wget ${FAPATH}chr4.fna.gz && gunzip chr4.fna.gz && mv chr4.fna chroms/chr4.fa
wget ${FAPATH}chr5.fna.gz && gunzip chr5.fna.gz && mv chr5.fna chroms/chr5.fa
wget ${FAPATH}chr6.fna.gz && gunzip chr6.fna.gz && mv chr6.fna chroms/chr6.fa
wget ${FAPATH}chr7.fna.gz && gunzip chr7.fna.gz && mv chr7.fna chroms/chr7.fa
wget ${FAPATH}chr8.fna.gz && gunzip chr8.fna.gz && mv chr8.fna chroms/chr8.fa
wget ${FAPATH}chr9.fna.gz && gunzip chr9.fna.gz && mv chr9.fna chroms/chr9.fa
wget ${FAPATH}chr10.fna.gz && gunzip chr10.fna.gz && mv chr10.fna chroms/chr10.fa
wget ${FAPATH}chr11.fna.gz && gunzip chr11.fna.gz && mv chr11.fna chroms/chr11.fa
wget ${FAPATH}chr12.fna.gz && gunzip chr12.fna.gz && mv chr12.fna chroms/chr12.fa
wget ${FAPATH}chr13.fna.gz && gunzip chr13.fna.gz && mv chr13.fna chroms/chr13.fa
wget ${FAPATH}chr14.fna.gz && gunzip chr14.fna.gz && mv chr14.fna chroms/chr14.fa
wget ${FAPATH}chr15.fna.gz && gunzip chr15.fna.gz && mv chr15.fna chroms/chr15.fa
wget ${FAPATH}chr16.fna.gz && gunzip chr16.fna.gz && mv chr16.fna chroms/chr16.fa
wget ${FAPATH}chr17.fna.gz && gunzip chr17.fna.gz && mv chr17.fna chroms/chr17.fa
wget ${FAPATH}chr18.fna.gz && gunzip chr18.fna.gz && mv chr18.fna chroms/chr18.fa
wget ${FAPATH}chr19.fna.gz && gunzip chr19.fna.gz && mv chr19.fna chroms/chr19.fa
wget ${FAPATH}chr20.fna.gz && gunzip chr20.fna.gz && mv chr20.fna chroms/chr20.fa
wget ${FAPATH}chr21.fna.gz && gunzip chr21.fna.gz && mv chr21.fna chroms/chr21.fa
wget ${FAPATH}chr22.fna.gz && gunzip chr22.fna.gz && mv chr22.fna chroms/chr22.fa
wget ${FAPATH}chrX.fna.gz && gunzip chrX.fna.gz && mv chrX.fna chroms/chrX.fa
wget ${FAPATH}chrY.fna.gz && gunzip chrY.fna.gz && mv chrY.fna chroms/chrY.fa 
```

### **R session setup:**

``` r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("pacman", quietly = TRUE))
    BiocManager::install("pacman")
library(pacman)
p_load(dplyr)
p_load(Biostrings)
p_load(polyester)
p_load(ggplot2)
p_load(tibble)
p_load(rtracklayer)
p_load(plyranges)
p_load(tidyr)
```

### **Initial Parameters:**

``` r
readlengths<-c(75,150)
fragsizes<-c(250,500)
nsim<-1000
nsamples<-20
seed<-20221107

gtf_file<-"Hsa38_only_transcript_support_level_1.gtf"
fastas_folder<-"chroms"
my_obj <- import(gtf_file)
```

### **RNA-seq reads simulation:**

``` r
for (i in 1:nsamples){
  
    ########################################
    # Transcript_Exon_Count.tab is a tab delimited file containing gene info and number 
    # of exon per transcript.
    # Select one random middle exon from transcripts with more than 2 exons.
    ########################################
  
    set.seed(seed+i)
    numexo_rand<-read.table("Transcript_Exon_Count.tab", col.names = c("Transcript","GeneName","GeneID","Nexo")) %>%
      as_tibble()  %>% filter(Nexo>2) %>%
      rowwise() %>% mutate(RandomExo=round(runif(1, min = 1.5, max=Nexo-0.5),0))
    
    ########################################
    # At the end only one transcript per gene should be selected per simulated sample, 
    # therefore genes are weighted based on the median number of possible selected 
    # exons in their transcripts. This means that genes with transcripts with many 
    # exons have a grater possibility o being selected.
    ########################################
    
    set.seed(seed+i)
    geneweights<-numexo_rand %>%
      group_by(GeneID) %>%
      summarize(Ntranscripts=dplyr::n(),medianExo=median(Nexo-2)) %>%
      ungroup() %>% mutate( prob=(medianExo)/sum(medianExo)) %>%
      slice_sample(n=nsim,weight_by = prob)
    
    ########################################
    # Once gene is selected, then a transcript out of it is selected, and again transcripts
    # with more exons are more likely to be selected. Then a random value between 0 and 1 
    # is assinged as UPFI, DOFI=1-UPFI for an uniform distribution simulation. 
    # Additionally a random number between 10 and 20 is selected simulating the 
    # expression/number of total copies of the transcript (UPFI+DOFI).
    ########################################
    
    set.seed(seed+i)
    numexo_rand_subset<- numexo_rand %>% right_join(geneweights, by="GeneID") %>%
                      group_by(GeneID) %>% mutate(prob2=((Nexo-2)/sum(Nexo-2))) %>%
                      slice_sample(n=1,weight_by = prob2) %>% ungroup()%>%
                      dplyr::select(Transcript,GeneName,GeneID,Nexo,RandomExo)%>%
                      column_to_rownames("Transcript")%>%
                      mutate(UPFI=round(runif(nsim),1), DOFI=1-UPFI,
                             Ncopies=round(runif(nsim,min = 10.5, max=20.5),0),
                             Ncopies_UPFI=round(UPFI*Ncopies,0),
                             Ncopies_DOFI=round(DOFI*Ncopies,0))
    ########################################
    # Selected Transcripts are filtred from the original GTF. 
    ########################################
    
    transcript_ofI<-my_obj %>%
      filter(transcript_id %in% rownames(numexo_rand_subset) & type %in% c("exon")) %>%
      split(f = .$transcript_id)
    
    ########################################
    # Random and up/down-stream exons  are filtered from the transcripts.
    ########################################
    
    exons_ofI<-lapply(names(transcript_ofI),function(x){
      transcript_ofI[[x]] %>% 
        plyranges::filter(exon_number == numexo_rand_subset[x,"RandomExo"])})
    names(exons_ofI)<-names(transcript_ofI)
    exons_ofI_upstream<-lapply(names(transcript_ofI),function(x){
      transcript_ofI[[x]] %>% 
        plyranges::filter(exon_number == (numexo_rand_subset[x,"RandomExo"]-1))})
    names(exons_ofI_upstream)<-names(transcript_ofI)
    exons_ofI_downstream<-lapply(names(transcript_ofI),function(x){
      transcript_ofI[[x]] %>%
        plyranges::filter(exon_number == (numexo_rand_subset[x,"RandomExo"]+1))})
    names(exons_ofI_downstream)<-names(transcript_ofI)
    
    ########################################
    # Extracting additional information and bulidng the simulation table output.
    ########################################
    
    exonlengt<-sapply(names(transcript_ofI),function(x){lengths(exons_ofI[[x]])})
    exonstart<-sapply(names(transcript_ofI),function(x){start(exons_ofI[[x]])})
    exonend<-sapply(names(transcript_ofI),function(x){end(exons_ofI[[x]])})
    intronuplengt<-sapply(names(transcript_ofI),function(x){
                          distance(exons_ofI[[x]], exons_ofI_upstream[[x]])})
    introndownlengt<-sapply(names(transcript_ofI),function(x){
                          distance(exons_ofI[[x]], exons_ofI_downstream[[x]])})
    numexo_rand_subset <- numexo_rand_subset %>% as_tibble(rownames = "Transcript") %>%
                  left_join(tibble("Transcript"=names(exonlengt),
                                   "ExoStart"=exonstart,"ExoEnd"=exonend,
                                   "UpIntronLength"=intronuplengt,
                                   "ExonLength"=exonlengt,"DownIntronLength"=introndownlengt),
                            by="Transcript")
    
    ########################################
    # Random selected exons are substituted by new exons, which are defined by merging
    # the up and down stream intron + exon of the exon of interest, emulating UPFI and DOFI 
    ########################################
    
    dofiexons<-data.frame(start=sapply(names(transcript_ofI),function(x){
                                        min(start(exons_ofI_upstream[[x]]),start(exons_ofI[[x]]))}),
                          end=sapply(names(transcript_ofI),function(x){
                                        max(end(exons_ofI_upstream[[x]]),end(exons_ofI[[x]]) )}),
                          seqnames=sapply(exons_ofI,chrom) %>% sapply(function(x){
                                        as.character(x@values)}),
                          strand=sapply(exons_ofI,strand)%>% sapply(function(x){
                                        as.character(x@values)}),
                          transcript_id= sapply(names(transcript_ofI),function(x){
                                        exons_ofI[[x]]$transcript_id}),
                          gene_id= sapply(names(transcript_ofI),function(x){
                                        exons_ofI[[x]]$gene_id})) %>%
              as_granges()%>% split(f = .$transcript_id)
    upfiexons<-data.frame(start=sapply(names(transcript_ofI),function(x){
                                        min(start(exons_ofI_downstream[[x]]),
                                            start(exons_ofI[[x]]))}),
                          end=sapply(names(transcript_ofI),function(x){
                                        max(end(exons_ofI_downstream[[x]]),
                                            end(exons_ofI[[x]]))}),
                          seqnames=sapply(exons_ofI,chrom) %>% sapply(function(x){
                                        as.character(x@values)}),
                          strand=sapply(exons_ofI,strand)%>% sapply(function(x){
                                        as.character(x@values)}),
                          transcript_id= sapply(names(transcript_ofI),function(x){
                                        exons_ofI[[x]]$transcript_id}),
                          gene_id= sapply(names(transcript_ofI),function(x){
                                        exons_ofI[[x]]$gene_id})) %>%
              as_granges()%>% split(f = .$transcript_id)
    
    ########################################
    # For building the UPFI and DOFI transcripts, firstly the matched exons are removed from
    # the original transcript and the newly merged exon (exon-intron-exon) is inserted. It
    # is important to sort the new exons in ascending order in the + strand and in 
    # descending order in the negative strand.
    ########################################
    
    transcript_upfi<-lapply(names(transcript_ofI),function(x){
        res<-filter_by_non_overlaps(transcript_ofI[[x]],upfiexons[[x]]) %>%
                c(upfiexons[[x]]) %>% select(transcript_id,gene_id) %>%
                mutate(transcript_id=paste0(transcript_id,"_upfi"),type="exon")
      if(length(unique(strand(res)))!=1){break}
      if(unique(strand(res))=="+"){
        res<- res%>%arrange(start(.))}else{
        res<-res[order(start(res),decreasing = T),]}
      return(res)
    })  %>% as(., "GRangesList") %>% unlist
    transcript_dofi<-lapply(names(transcript_ofI),function(x){
      res<-filter_by_non_overlaps(transcript_ofI[[x]],dofiexons[[x]]) %>%
              c(dofiexons[[x]]) %>% select(transcript_id,gene_id) %>%
              mutate(transcript_id=paste0(transcript_id,"_dofi"),type="exon")
      if(length(unique(strand(res)))!=1){break}
      if(unique(strand(res))=="+"){
        res<- res%>%arrange(start(.))}else{
        res<-res[order(start(res),decreasing = T),]}
      return(res)
    }) %>% as(., "GRangesList") %>% unlist
    
    ########################################
    # Join UPFI and DOFI transcripts and write a gtf file. 
    ########################################
    
    all_transcripts<-c(transcript_upfi,transcript_dofi)
    write_gff(all_transcripts,file = paste0("simulations_",i,".gtf"))
    
    ########################################
    # Get the sequence from the GTF and write it in disk
    ########################################
    
    sequen<-seq_gtf(gtf=paste0("simulations_",i,".gtf"), seqs=fastas_folder)
    writeXStringSet(x=sequen, filepath= paste0("simulations_",i,".fa"))

    ########################################
    # Parse the simulated table to obtain UPFI and DOFI expression. 
    ########################################
    
    temp4libsize<-numexo_rand_subset %>%
      dplyr::select(Transcript,Ncopies_UPFI,Ncopies_DOFI) %>%
      pivot_longer(-c(Transcript), values_to = "Copies", names_to = "Type")%>%
      mutate(Transcript=ifelse(grepl(pattern = "UPFI",Type),
                                        paste0(Transcript,"_upfi"),
                                        paste0(Transcript,"_dofi")))
   
    ########################################
    # Based on the number of copies, the transcript and read lengths the number of reads needs
    # to be calculated. 
    ########################################
    
    temp4libsize<-tibble(Transcript= names(sequen), Lengths=lengths(sequen)) %>%
      left_join(temp4libsize, by="Transcript",multiple="all") %>%
      mutate(R_75=round((Lengths/75)*Copies,0),
             R_150=round((Lengths/150)*Copies,0)) %>%
      mutate(Type2="uniform") 
    
    ########################################
    # Reads for UPFI and DOFI transcripts are simulated, for the selected readlengths
    ########################################
    
    for(readl in readlengths ){
      
        ########################################
        # if readlength > 102 then a uniform error model is applied rather than the illumina5
        # that comes with polyester package. 
        ########################################
      
        errormodel<-ifelse(readl>102,"uniform","illumina5") 
        a<-temp4libsize %>% select(Transcript,contains(as.character(readl))) 
        countmat<-a%>% column_to_rownames("Transcript")
        colnames(countmat)<-"expression"
        tibblecounts<- a %>% mutate(x=ifelse(grepl(Transcript,pattern = "upfi"),
                                             "UPFI_counts","DOFI_counts"),
                                    Transcript = gsub("_upfi|_dofi","",Transcript)) %>%
                      pivot_wider(Transcript,
                                  names_from = x, values_from = paste0("R_",readl))
        outdir_temp<-paste0("simulations/sim",i,"/R_",readl,"_SE/")
        
        ########################################
        # First SE reads are simulated.
        ########################################
        
        simulate_experiment_countmat(readmat = as.matrix(countmat),
                                     fasta=paste0("simulations_",i,".fa"),
                                     outdir=outdir_temp,  readlen = as.numeric(readl),
                                     error_model = errormodel,paired = FALSE,
                                     strand_specific=F,gzip=TRUE,
                                     error_rate = 0.005,seed = seed+i, bias="none")
        
        ########################################
        # Thene PE for the different insert sizes. 
        ########################################
        
        for(fragsize in fragsizes){
          outdir_temp<-paste0("simulations/sim",i,"/R_",readl,"_Fsize_",fragsize,"_PE/")
          simulate_experiment_countmat(readmat = as.matrix(countmat),
                                      fasta=paste0("simulations_",i,".fa"),
                                      outdir=outdir_temp, fraglen = fragsize,
                                      fragsd = 25, readlen = as.numeric(readl),
                                      error_model = errormodel,paired = TRUE,
                                      strand_specific=F,gzip=TRUE,
                                      error_rate = 0.005,seed = seed+i,bias="none")
        }
        file_count<-paste0("simulations/sim",i,"/R_",readl,"_countmat.txt")
        numexo_rand_subset %>% 
          left_join(tibblecounts,by = "Transcript") %>%
          write.table(file=file_count,quote=F,sep="\t", row.names = F)
      }
}
```

### **Alignment of simulated reads (bash):**

``` bash
for rlen in 75 150
do
for sim in {1..20}
do
STAR --readFilesIn fastas/sim${sim}/R_${rlen}_SE/*fasta.gz \
 --runThreadN 4 --chimSegmentMin 20 \
 --alignSJoverhangMin 8 --twopassMode Basic \
 --alignIntronMin 20 --alignIntronMax 200000 \
 --alignSJDBoverhangMin 5 --outFilterMultimapNmax 100 \
 --outFileNamePrefix results/sim${sim}/R_${rlen}_SE/ \
 --quantMode GeneCounts --outSAMattributes All \
 --runMode alignReads --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverLmax 0.1 \
 --seedSearchStartLmax 25 --outSAMtype BAM SortedByCoordinate \
 --readFilesCommand zcat \
 --outSAMstrandField intronMotif --genomeDir hg38_reference

for insize in 250 500
do
STAR --readFilesIn fastas/sim${sim}/R_${rlen}_Fsize_${insize}_PE/*fasta.gz \
 --runThreadN 4 --chimSegmentMin 20 \
 --alignSJoverhangMin 8 --twopassMode Basic \
 --alignIntronMin 20 --alignIntronMax 200000 \
 --alignSJDBoverhangMin 5 --outFilterMultimapNmax 100 \
 --outFileNamePrefix results/sim${sim}/R_${rlen}_Fsize_${insize}_PE/ \
 --quantMode GeneCounts --outSAMattributes All \
 --runMode alignReads --outFilterMismatchNmax 999 \
 --outFilterMismatchNoverLmax 0.1 \
 --seedSearchStartLmax 25 --outSAMtype BAM SortedByCoordinate \
 --readFilesCommand zcat \
 --outSAMstrandField intronMotif --genomeDir hg38_reference
done
done
done
```
