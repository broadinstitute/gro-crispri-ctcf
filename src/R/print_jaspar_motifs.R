#########################################################
## Print the JASPAR motifs for the selected TFs


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("JASPAR2020")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("TFBSTools")

## Functions are defined at the end of the file

##The inputs of this file come from the outputs of filter_motif_by_pam.py

library(JASPAR2020)
library(TFBSTools)
library(ggseqlogo)

tf.strong.candidates <- read.table('tf_candidates.txt')

logo.list <- list()
no.logo.list <- list()

tfs.type <- c()

# Extract the positive matches
tfs.table <- matrix(ncol=2)
no.tfs.table <- matrix(ncol=2)
for (tf in tf.strong.candidates$V1){
  pfm <- getMatrixByName(JASPAR2020, name=tf)
  icm.pfm <- toICM(pfm)
  ## Filter by ICM matrix: need at least 5 strong nucleotides in the guide
  is.logo = filter.logo.list(icm.pfm@profileMatrix)
  if (is.logo){
    logo.list[[tf]] = icm.pfm@profileMatrix
    tfs.type <- c(tfs.type, icm.pfm@matrixClass)
    tfs.table <- rbind(tfs.table, c(tf, icm.pfm@matrixClass))
  } else{
    no.logo.list[[tf]] = icm.pfm@profileMatrix
    no.tfs.table <- rbind(no.tfs.table, c(tf, icm.pfm@matrixClass))
  }
}

## Plot
pdf('tf_candidates.pdf', height = 8, width = 7)
ggseqlogo(logo.list[1:16], ncol=2)
ggseqlogo(logo.list[17:32], ncol=2)
ggseqlogo(logo.list[33:48], ncol=2)
ggseqlogo(logo.list[49:65], ncol=2)
ggseqlogo(logo.list[66:82], ncol=2)
ggseqlogo(logo.list[83:88], ncol=2)
dev.off()


tf.negatives.strong.candidates <- read.table('tf_negative_candidates.txt')

for (tf in tf.negatives.strong.candidates$V1){
  pfm <- getMatrixByName(JASPAR2020, name=tf)
  icm.pfm <- toICM(pfm)
  no.logo.list[[tf]] = icm.pfm@profileMatrix
  no.tfs.table <- rbind(no.tfs.table, c(tf, icm.pfm@matrixClass))
}

pdf('tf_negative_candidates.pdf', height = 8, width = 7)
ggseqlogo(no.logo.list[1:16], ncol=2)
ggseqlogo(no.logo.list[17:32], ncol=2)
ggseqlogo(no.logo.list[33:48], ncol=2)
ggseqlogo(no.logo.list[49:64], ncol=2)
ggseqlogo(no.logo.list[65:80], ncol=2)
ggseqlogo(no.logo.list[81:96], ncol=2)
ggseqlogo(no.logo.list[97:112], ncol=2)
ggseqlogo(no.logo.list[113:128], ncol=2)
ggseqlogo(no.logo.list[129:144], ncol=2)
ggseqlogo(no.logo.list[145:160], ncol=2)
ggseqlogo(no.logo.list[161:176], ncol=2)
ggseqlogo(no.logo.list[177:192], ncol=2)
ggseqlogo(no.logo.list[193:208], ncol=2)
ggseqlogo(no.logo.list[209:224], ncol=2)
ggseqlogo(no.logo.list[225:240], ncol=2)
ggseqlogo(no.logo.list[241:256], ncol=2)
ggseqlogo(no.logo.list[257:272], ncol=2)
ggseqlogo(no.logo.list[273:288], ncol=2)
ggseqlogo(no.logo.list[289:304], ncol=2)
ggseqlogo(no.logo.list[305:320], ncol=2)
ggseqlogo(no.logo.list[321:336], ncol=2)
ggseqlogo(no.logo.list[337:352], ncol=2)
ggseqlogo(no.logo.list[353:368], ncol=2)
ggseqlogo(no.logo.list[369:384], ncol=2)
ggseqlogo(no.logo.list[385:400], ncol=2)
ggseqlogo(no.logo.list[401:416], ncol=2)
ggseqlogo(no.logo.list[417:432], ncol=2)
ggseqlogo(no.logo.list[433:448], ncol=2)
ggseqlogo(no.logo.list[449:464], ncol=2)
ggseqlogo(no.logo.list[465:480], ncol=2)
ggseqlogo(no.logo.list[481:496], ncol=2)
ggseqlogo(no.logo.list[497:512], ncol=2)
ggseqlogo(no.logo.list[513:528], ncol=2)
ggseqlogo(no.logo.list[529:544], ncol=2)
ggseqlogo(no.logo.list[545:561], ncol=2)
ggseqlogo(no.logo.list[562:578], ncol=2)
ggseqlogo(no.logo.list[579:595], ncol=2)
ggseqlogo(no.logo.list[596:612], ncol=2)
ggseqlogo(no.logo.list[613:629], ncol=2)
ggseqlogo(no.logo.list[630:646], ncol=2)
ggseqlogo(no.logo.list[647:658], ncol=2)
dev.off()



table(tfs.type)



###########

write.table(no.tfs.table, file = 'no_tfs_table.txt',
            sep='\t',col.names = F, row.names = F, quote = F)
write.table(tfs.table, file = '.tfs_table.txt',
            sep='\t',col.names = F, row.names = F, quote = F)
write.table(table(tfs.type), file = 'tfs_type.txt',
            sep='\t',col.names = F, row.names = F, quote = F)


########## Functions #######

filter.logo.list <- function(profile.matrix){
  
  is.logo = FALSE
   
  mid = ncol(profile.matrix)/2
  upper.mid.matrix <- profile.matrix[,c(1:(as.integer(mid)))]
  lower.mid.matrix <- profile.matrix[,c((as.integer(mid)+1):ncol(profile.matrix))]
  
  ## Upper matrix
  top.letter <- apply(upper.mid.matrix, 2, function(x){if(max(x)>1){return(which.max(x))}else{return(NA)}})
  letter = 2 #C
  letter.index <- which(rle(top.letter)$values==letter&rle(top.letter)$lengths>=2)
  if (length(letter.index)>0){
    letter.index <- min(letter.index)
    repeated.letter <- rle(top.letter)$values[letter.index]
    if(length(repeated.letter)>0){
      sub.matrix <- profile.matrix[,letter.index:ncol(profile.matrix)]
      if (sum(apply(sub.matrix, 2, max)>0.9)>5){
        is.logo = TRUE
      }
    }
  }
  
  ## Lower matrix
  top.letter <- apply(lower.mid.matrix, 2, function(x){if(max(x)>1){return(which.max(x))}else{return(NA)}})
  letter = 3 #G
  letter.index <- which(rle(top.letter)$values==letter&rle(top.letter)$lengths>=2)
  if (length(letter.index)>0){
    letter.index <- max(which(rle(top.letter)$values==letter&rle(top.letter)$lengths>=2))
    repeated.letter <- rle(top.letter)$values[letter.index]
    if(length(repeated.letter)>0){
      index.end <- sum(rle(top.letter)$lengths[1:letter.index]) + as.integer(mid)
      sub.matrix <- profile.matrix[,1:index.end]
      if (sum(apply(sub.matrix, 2, max)>0.9)>5){
        is.logo = TRUE
      }
    }
  }
  return(is.logo)
}





