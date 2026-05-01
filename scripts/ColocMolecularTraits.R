library(tidyverse)
library(data.table)
library(coloc)
library(optparse)

########## FUNCTIONS ############

extract_credible_sets <- function(finemapping_res) {
valid_credible_sets <- finemapping_res %>% 
    distinct(cs_index) %>% 
    mutate(cs_index = as.numeric(str_remove(cs_index,'L'))) %>% 
    pull(cs_index)  
valid_credible_sets  
}

run_coloc <- function(SusielbfParquet_protein,
                      SusielbfParquet_expression,
                      SusieParquet_protein,
                      SusieParquet_expression,
                      gene,
                      data_dir = 'coloc_temp/',
                      coloc_out = 'coloc_out/'
                     ) {
    message('Loading fine-mapping res')
    pqtl_cs <- arrow::read_parquet(SusieParquet_protein) %>% 
                    extract_credible_sets()
    eqtl_cs <- arrow::read_parquet(SusieParquet_expression) %>% 
                    extract_credible_sets()
    number_credible_sets_eqtl <- eqtl_cs %>% length()
    number_credible_sets_pqtl <- pqtl_cs %>% length()
    
    message(paste0('number credible sets eqtl:',number_credible_sets_pqtl ))
    message(paste0('number credible sets pqtl:',number_credible_sets_eqtl ))
    
    
    if (number_credible_sets_eqtl > 0 & number_credible_sets_pqtl > 0) {
    message('Loading lbf parquets')
    lbf_protein <- arrow::read_parquet(paste0(data_dir,basename(SusielbfParquet_protein)))
    lbf_expression <- arrow::read_parquet(paste0(data_dir,basename(SusielbfParquet_expression)))
    
    lbf_matrix_eqtl <- lbf_expression %>% 
                    select(variant,contains('lbf')) %>% 
                    column_to_rownames('variant') %>%
                    dplyr::select(eqtl_cs) %>% 
                    t()  
    lbf_matrix_pqtl <- lbf_protein %>% 
                    select(variant,contains('lbf')) %>% 
                    column_to_rownames('variant') %>% 
                    dplyr::select(eqtl_cs) %>% 
                    t() 
    message('Performing coloc')
    res <- coloc.bf_bf(lbf_matrix_eqtl,lbf_matrix_pqtl)
    res_out <- res$summary %>% mutate(gene = gene)
    res_out
        }
    }


################ PARSE COMMAND LINE #############
option_list <- list(
  make_option("--gene_id", type = "character"),
  make_option("--eqtl_susie", type = "character", default = NA),
  make_option("--sqtl_susie", type = "character", default = NA),
  make_option("--pqtl_susie", type = "character", default = NA),
  make_option("--output", type = "character"),
)

opt <- parse_args(OptionParser(option_list = option_list))Splice

load_susie_safe <- function(path) {
  if (is.na(path) || path == "") return(NULL)
  readRDS(path)
}

eqtl <- load_susie_safe(opt$eqtl_susie)
sqtl <- load_susie_safe(opt$sqtl_susie)
pqtl <- load_susie_safe(opt$pqtl_susie)


susie_list <- list(
  eQTL = eqtl,
  sQTL = sqtl,
  pQTL = pqtl
)

# drop missing
susie_list <- susie_list[!sapply(susie_list, is.null)]

