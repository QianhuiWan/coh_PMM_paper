

# prepare SAF file for featureCount: 
## TE up&down 100bp SAF files + 5683690 random regions up&down 100bp SAF files
## strand not considered for simplify this step
## strand will be considered in later steps
## + start 100bp and end 100bp within TE SAF files

# load pkgs ####################################################################
library(tidyverse)
library(data.table)
library(magrittr)
library(plyranges)
library(BSgenome.Hsapiens.UCSC.hg38)

# load data ####################################################################
## load rmsk
rmsk_hg38_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_repeatMasker_hg38"
rmsk_hg38_te_ids_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_hg38_rmsk_teCoord_ids.txt"

rmsk_hg38 <- fread(rmsk_hg38_path)
rmsk_hg38_te_ids <- fread(rmsk_hg38_te_ids_path,
                          col.names = c("genoName", "genoStart", "genoEnd", 
                                        "strand", "te_id")) %>%
  mutate(genoEnd = genoEnd-1)

rmsk_hg38_addAnno <- rmsk_hg38 %>%
  dplyr::mutate(div = milliDiv/1000) %>% ## Convert percent divergence of RepeatMasker output to divergence
  dplyr::mutate(div_jc = -(0.75*log(1 - (4*div)/3))) %>% ## Apply Juke's Cantor correction to divergence distance
  dplyr::mutate(age=(div_jc*1000) / 2.2) %>%
  # dplyr::filter(age<=40) %>%
  dplyr::left_join(rmsk_hg38_te_ids,
                   by = c("genoName" ="genoName", "genoStart"="genoStart",
                          "genoEnd"="genoEnd", "strand" = "strand")) %>%
  as.data.frame() %>%
  dplyr::select(genoName, genoStart, genoEnd, strand,
                te_id, repName, repClass, repFamily, age)
write_tsv(rmsk_hg38_addAnno,
          file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_repeatMasker_hg38_addTEid.tsv")
## rmsk is 0-based, change to 1-based for later SAF
rmsk_hg38_addAnno_GR <- rmsk_hg38_addAnno %>% 
  `colnames<-`(c("seqnames", "start", "end", "strand", 
                 "te_id", "repName", "repClass", "repFamily", "age")) %>% 
  dplyr::mutate(start = start+1) %>% # 0 to 1 based
  as_granges()

## get hg38 seqinfo
new_seqinfo <- Seqinfo(seqnames = seqnames(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:25],
                       seqlengths = seqlengths(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:25],
                       isCircular = isCircular(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:25],
                       genome = genome(seqinfo(BSgenome.Hsapiens.UCSC.hg38))[1:25])
## Replace Seqinfo with the new Seqinfo object
seqlevels(rmsk_hg38_addAnno_GR) <- seqlevels(new_seqinfo)
seqinfo(rmsk_hg38_addAnno_GR) <- new_seqinfo

###################### get TE up&down 100bp SAF files ##########################
# get TE upstream 100bp regions SAF ("GeneID", "Chr", "Start","End", "Strand")
## start-100 to start
rmsk_hg38_start_regions <- rmsk_hg38_addAnno_GR %>%
  mutate(end = start, 
         # ensure no negative start 
         start = pmax(start - 100, 1)) %>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()

## save to SAF
rmsk_hg38_start_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = te_id, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_rmsk_hg38_start100bp.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")

# get TE downstream 100bp regions SAF
## end to end+100
rmsk_hg38_end_regions <- rmsk_hg38_addAnno_GR %>%
  mutate(start = end, 
         # end = end + 100) %>% 
         # Ensure end <= chromosome length
         end = pmin(end + 100, seqlengths(.)[as.character(seqnames(.))])) %>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()

## save to SAF
rmsk_hg38_end_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = te_id, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_rmsk_hg38_end100bp.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")

###################### get within TE start & end 100bp SAF files ###############
# within TE, not up/down stream 100bp
# get TE start 100bp within TE regions SAF ("GeneID", "Chr", "Start","End", "Strand")
## start to start+100 regions
rmsk_hg38_start100inTE_regions <- rmsk_hg38_addAnno_GR %>%
  mutate(end = start + 100) %>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()

## save to SAF
rmsk_hg38_start100inTE_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = te_id, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_rmsk_hg38_start100bpInTE.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")

# get TE end 100bp within TE regions SAF
## end-100 to end
rmsk_hg38_end100inTE_regions <- rmsk_hg38_addAnno_GR %>%
  mutate(start = pmax(end - 100, 1)) %>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()

## save to SAF
rmsk_hg38_end100inTE_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = te_id, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_rmsk_hg38_end100bpInTE.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")


###################### get Random up&down 100bp SAF files ######################
## load random regions
rendomRegion_path <- "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_p14/random_regions_matchTElengths.bed"
random_hg38 <- fread(rendomRegion_path)
colnames(random_hg38) <- c("seqnames", "start", "end", "strand", "regionID")

random_hg38_GR <- random_hg38 %>% 
  dplyr::mutate(start = start+1) %>% # 0 to 1 based
  as_granges()

## Replace Seqinfo with the new Seqinfo object
seqlevels(random_hg38_GR) <- seqlevels(new_seqinfo)
seqinfo(random_hg38_GR) <- new_seqinfo

## Randomly select 5683690 GRanges
set.seed(5566)  # Set seed for reproducibility
random_hg38_GR <- 
  random_hg38_GR[sample(seq_along(random_hg38_GR), 
                        size = 5683690, replace = FALSE),] 

random_hg38_GR %>% GenomicRanges::sort() %>% as.data.frame() %>%
  write_tsv("/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_randomRegions5683690_GR_hg38.tsv")

# get random region upstream 100bp regions SAF
## start-100 to start
random_hg38_start_regions <- random_hg38_GR %>%
  mutate(end = start, start = pmax(start - 100, 1)) %>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()

## save to SAF
random_hg38_start_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = regionID, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_random5683690_hg38_start100bp.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")

# get random region downstream 100bp regions SAF
## end to end+100
random_hg38_end_regions <- random_hg38_GR %>%
  mutate(start = end, 
         # Ensure end <= chromosome length
         end = pmin(end + 100, seqlengths(.)[as.character(seqnames(.))]))%>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()


## save to SAF
random_hg38_end_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = regionID, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_random5683690_hg38_end100bp.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")


###################### get within random start & end 100bp SAF files ###########
# within random, not up/down stream 100bp
# get random start 100bp within random regions SAF 
# ("GeneID", "Chr", "Start","End", "Strand")
## start to start+100 regions
random_hg38_start100inRandom_regions <- random_hg38_GR %>%
  mutate(end = start + 100) %>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()

## save to SAF
random_hg38_start100inRandom_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = regionID, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_random5683690_hg38_start100bpInRandom.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")

# get TE end 100bp within TE regions SAF
## end-100 to end
random_hg38_end100inRandom_regions <- random_hg38_GR %>%
  mutate(start = pmax(end - 100, 1)) %>% 
  GenomicRanges::trim() %>% 
  GenomicRanges::sort()

## save to SAF
random_hg38_end100inRandom_regions %>% as.data.frame() %>% 
  dplyr::select(GeneID = regionID, 
                Chr = seqnames, Start = start, End = end, Strand = strand) %>% 
  write.table(file = "/net/nfs-irwrsrchnas01/labs/dschones/bioresearch/qianhui/hg38_2024/hg38_rmskInfo/filtered_random5683690_hg38_end100bpInRandom.saf", 
              row.names = F, col.names = T, 
              quote = F, sep = "\t")

# check regioin ids
table(random_hg38_start_regions$regionID%in%random_hg38_end_regions$regionID)
table(random_hg38_start100inRandom_regions$regionID%in%random_hg38_end100inRandom_regions$regionID)


