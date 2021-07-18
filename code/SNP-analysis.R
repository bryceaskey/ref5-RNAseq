library(stringr)

# Load gene annotations
araport11 <- read.csv(file="C:/Users/Bryce/Research/ref5-RNAseq/data/Araport11.csv", header=FALSE)
araport11 <- araport11[, c(1, 3, 4, 6, 7, 8, 13)]
colnames(araport11) <- c("locus", "short_name", "name", "chromosome", "start", "end", "aliases")

# Load variant data - SNPs and indels
H10 <- read.table(file="C:/Users/Bryce/Research/ref5-RNAseq/data/variants/H10.vcf")[ ,c(1,2,4,5,6,8)]
colnames(H10) <- c("chromosome", "position", "reference", "alternative", "quality", "info")

H55 <- read.table(file="C:/Users/Bryce/Research/ref5-RNAseq/data/variants/H55.vcf")[ ,c(1,2,4,5,6,8)]
colnames(H55) <- c("chromosome", "position", "reference", "alternative", "quality", "info")

R5 <- read.table(file="C:/Users/Bryce/Research/ref5-RNAseq/data/variants/R5.vcf")[ ,c(1,2,4,5,6,8)]
colnames(R5) <- c("chromosome", "position", "reference", "alternative", "quality", "info")

WT <- read.table(file="C:/Users/Bryce/Research/ref5-RNAseq/data/variants/WT.vcf")[ ,c(1,2,4,5,6,8)]
colnames(WT) <- c("chromosome", "position", "reference", "alternative", "quality", "info")

Y6 <- read.table(file="C:/Users/Bryce/Research/ref5-RNAseq/data/variants/Y6.vcf")[ ,c(1,2,4,5,6,8)]
colnames(Y6) <- c("chromosome", "position", "reference", "alternative", "quality", "info")

# Associate genes with variants
add_variants <- function(genome, genotype_name, variants){
  genome[[paste(genotype_name, "_variants", sep="")]] <- NA
  for(gene in 1:nrow(genome)){
    variant_subset <- variants[variants$chromosome == genome$chromosome[gene], ]
    variant_positions <- variant_subset$position[variant_subset$position >= genome$start[gene] & variant_subset$position <= genome$end[gene]]
    variant_refs <- variant_subset$reference[variant_subset$position >= genome$start[gene] & variant_subset$position <= genome$end[gene]]
    variant_alts <- variant_subset$alternative[variant_subset$position >= genome$start[gene] & variant_subset$position <= genome$end[gene]]
    if(length(variant_positions) == 0){
      variant_data <- NA
    }else if(length(variant_positions) == 1){
      variant_data <- paste(variant_positions, "(", variant_refs, ":", variant_alts, ")", sep="")
    }else{
      variant_array <- array(c(variant_positions, variant_refs, variant_alts), dim=c(length(variant_positions), 3))
      variant_data <- toString(paste(variant_array[, 1], "(", variant_array[, 2], ":", variant_array[, 3], ")", sep=""))
    }
    #print(variant_data)
    genome[[paste(genotype_name, "_variants", sep="")]][gene] <- variant_data
  }
  return(genome)
}

araport11 <- add_variants(araport11, "H10", H10)
araport11 <- add_variants(araport11, "H55", H55)
araport11 <- add_variants(araport11, "R5", R5)
araport11 <- add_variants(araport11, "WT", WT)
araport11 <- add_variants(araport11, "Y6", Y6)

# Determine variant intersections between genotypes
variant_intersections <- function(genome, intersect, exclude, mutation_type){
  # Reorganize variant columns so that genotypes to intersect from are first, and genotypes to exclude from are last
  genotypes <- unlist(strsplit(colnames(genome)[8:length(colnames(genome))], split="_variants"))
  genotypes_to_drop <- setdiff(genotypes, intersect)
  genome_subset <- genome[ , !colnames(genome) %in% paste(genotypes_to_drop, "_variants", sep="")]
  for(genotype in exclude){
    genome_subset[[paste(genotype, "_variants", sep="")]] <- genome[[paste(genotype, "_variants", sep="")]]
  }
  
  genome_subset <- genome_subset[rowSums(!is.na(genome_subset[, 8:(8 + length(intersect) - 1)])) == length(intersect), ]

  if(length(exclude) > 1){
    genome_subset <- genome_subset[rowSums(!is.na(genome_subset[, (8 + length(intersect)):ncol(genome_subset)])) == 0, ]
  }else{
    genome_subset <- genome_subset[!is.na(genome_subset[, (8 + length(intersect)):ncol(genome_subset)]) == 0, ]
  }
  
  # Filter for mutation type
  mutation_match <- data.frame(array(data=NA, dim=c(nrow(genome_subset), length(intersect))))
  colnames(mutation_match) <- intersect
  for(genotype in intersect){
    for(gene in 1:nrow(genome_subset)){
      gene_variants <- strsplit(genome_subset[[paste(genotype, "_variants", sep="")]][gene], split=", ")[[1]]
      if(length(gene_variants) == 1){
        variant_type <- str_extract(gene_variants, pattern="(?<=\\().*(?=\\))")
      }else{
        variant_type <- unlist(lapply(gene_variants, str_extract, pattern="(?<=\\().*(?=\\))"))
      }
      mutation_match[[genotype]][gene] <- sum(mutation_type %in% variant_type)
    }
  }
  
  genome_subset <- genome_subset[apply(mutation_match, 1, all), ]
  return(genome_subset)
}

H10H55intersect_R5exclude <- variant_intersections(araport11, intersect=c("H10", "H55"), exclude="R5", mutation_type=c("G:T", "A:C"))
