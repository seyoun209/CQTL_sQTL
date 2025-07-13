#encode data only  (eclip)
library(utils)

setwd("/work/users/s/e/seyoun/CQTL_sQTL")

chromosome_names <- c(
  paste0("chr", 1:22),
  "chrX",
  "chrY"
)

meta <- fread("/work/users/s/e/seyoun/CQTL_sQTL/external_data/encode_rbp/metadata.tsv")
targetInfo <-fread("/work/users/s/e/seyoun/CQTL_sQTL/external_data/encode_rbp/targetInfo_encode.txt")
#target <- fread("/work/users/s/e/seyoun/CQTL_sQTL/external_data/encode_rbp/targetInfo.txt")

directory <- "/work/users/s/e/seyoun/CQTL_sQTL/external_data/encode_rbp"
file_list <- list.files(directory, pattern = "bed\\.gz$", full.names = TRUE)


#Finding target and gene
accession_code <- sub(".*?/(ENCF\\d+)\\.bed\\.gz", "\\1", basename(file_list[i]))
accession_code <- sub("\\.bed\\.gz$", "", accession_code)

exp_accession <- meta[which(meta$`File accession` == accession_code),`Experiment accession`]
gene_cellLine <- targetInfo |> dplyr::filter(Accession == exp_accession) |> 
  dplyr::select("Target gene symbol","Biosample term name") |> 
  dplyr::rename("rbp_nm" ="Target gene symbol","cell_line" = "Biosample term name")


list_rbp <- list()
rbp_all.df <- data.frame()
for (i in 1:length(file_list)) {
  print(i)
  accession_code <- sub(".*?/(ENCF\\d+)\\.bed\\.gz", "\\1", basename(file_list[i]))
  accession_code <- sub("\\.bed\\.gz$", "", accession_code)
  exp_accession <- meta[which(meta$`File accession` == accession_code),`Experiment accession`]
  bio_rep <- meta[which(meta$`File accession` == accession_code),`Biological replicate(s)`]
  bio_rep_n <-ncol(str_split(bio_rep,",",simplify = TRUE))
  gene_cellLine <- targetInfo |> dplyr::filter(Accession == exp_accession) |> 
    dplyr::select("Target gene symbol","Biosample term name") |> 
    dplyr::rename("rbp_nm" ="Target gene symbol","cell_line" = "Biosample term name")
  
  bed_data_temp <- fread(file_list[i]) |> dplyr::select(V1,V2,V3,V4,V6,V7,V8)
  if(bed_data_temp$V4[1] == "."){
    bed_data <- bed_data_temp |>
      dplyr::rename("chr" = "V1",
                    "start" ="V2",
                    "end" ="V3",
                    "strand" = "V6" ,
                    "-log2FC" = "V7",
                    "-log10pval" = "V8") |> dplyr::select(-"V4") |> 
      mutate(rbp_nm = rep(gene_cellLine$rbp_nm),
             cell_line= rep(gene_cellLine$cell_line),
             rep_n=rep(bio_rep_n) )
  }else{
  bed_data <- bed_data_temp |> mutate(V1 = str_split(V1,"_",simplify = TRUE)[,1],
           rbp_nm = str_split(V4,"_",simplify = TRUE)[,1],
  cell_line = str_split(V4,"_",simplify = TRUE)[,2],
  rep_n = gsub("rep","",str_split(V4,"_",simplify = TRUE)[,3]) |> as.double() ) |>
    dplyr::rename("chr" = "V1",
                  "start" ="V2",
                  "end" ="V3",
                  "strand" = "V6" ,
                  "-log2FC" = "V7",
                  "-log10pval" = "V8") |> dplyr::select(-"V4")
}
  list_rbp[[i]] <- bed_data
  names(list_rbp)[i] <- paste0(bed_data[1,]$cell_line,"_", bed_data[1,]$rbp_nm,"_",accession_code)
  
  rbp_all.df <-  bind_rows(rbp_all.df,bed_data)
}


write.table(rbp_all.df,file="external_data/encode_rbp/rbp_all_dataframe.txt",sep="\t",quote=FALSE, col.names = TRUE, row.names = FALSE)
save(list_rbp,file="external_data/encode_rbp/rbplist")

rbp_all.df <- fread("external_data/encode_rbp/rbp_all_dataframe.txt")

# I am going to seperate by rbp and cell lines 
top_rbps <- rbp_all.df %>%
  group_by(rbp_nm) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

cell_line_counts <- rbp_all.df %>%
  group_by(cell_line) %>%
  summarise(count = n()) %>%
  arrange(desc(count))



# Step 2: Function to create bed files for a given tissue
create_bed_files <- function(data, tissue_name) {
  data %>%
    dplyr::filter(cell_line == tissue_name) %>%
    group_by(rbp_nm) %>%
    do({
      bed_data <- select(., "chr","start","end","strand","rbp_nm","cell_line",
                         "rep_n","-log2FC","-log10pval")
      file_name <- paste0("external_data/encode_rbp/rbp_prep/", tissue_name, "_", .$rbp_nm[1], ".bed")
      write_tsv(bed_data, file_name, col_names = FALSE, quote = "none")
      tibble()
    })
}

# Step 3: Create bed files for each tissue
tissues_to_process <- c("K562", "HepG2", "SM-9MVZL")



for (tissue in tissues_to_process) {
  create_bed_files(rbp_all.df, tissue)
  cat("Completed processing for", tissue, "\n")
}


# Function to create bed files for each RBP
create_rbp_bed_files <- function(data) {
  data %>%
    group_by(rbp_nm) %>%
    do({
      bed_data <- select(., "chr","start","end","strand","rbp_nm","cell_line",
                         "rep_n","-log2FC","-log10pval")
      file_name <- paste0("external_data/encode_rbp/rbp_prep_rbpOnly/", unique(.$rbp_nm), ".bed")
      write_tsv(bed_data, file_name, col_names = FALSE, quote = "none")
      tibble()
    })
}

# Create RBP bed files
create_rbp_bed_files(rbp_all.df)






