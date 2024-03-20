library(coloc)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(locuscomparer)

# ## 函数的输入
# gwas_path <-
#   "data/analysis/DNAJC16/coloc/sumstats_coloc.feather"  # 改成直接传入df
# eQTLs_path <-
#   "data/analysis/DNAJC16/coloc/eQTL_coloc.feather"  # 改成直接传入df
# output_folder <- "data/analysis/DNAJC16/coloc/output" # 保持不变
# n_gwas <- 4869  # 用户传入
# n_eQTL <- 73  # 用户传入
# gwas_type <- "quant"  # 用户传入
# gwas_sdY <- 1  # 用户传入
# data_eQTLs <-
#   feather::read_feather(eQTLs_path)
# data_gwas <-
#   feather::read_feather(gwas_path)

runColocAnalysis <-
  function(data_gwas,
           data_eQTLs,
           output_folder,
           n_gwas,
           n_eQTL,
           gwas_sdY,
           gwas_type
           ) {

    "
    data_gwas: GWAS summary statistics, must have columns: snp, beta, varbeta, position, pval
    data_eQTLs: eQTL summary statistics, must have columns: snp, beta, varbeta, maf, position, pval, rsid, gene_id
    
    "
    
    
    if (!dir.exists(output_folder)) {
      dir.create(output_folder)
    }
    
    data_gwas$pvalue <- data_gwas$pval  #TODO: pvalue cols => pval
    gene_id_list <- unique(data_eQTLs$gene_id)
    
    ### get a gene eQTLs
    coloc_res <- list()
    coloc_snp_res <- list()
    
    for (idx in 1:length(gene_id_list)) {
      current_gene_id <- gene_id_list[idx]
      current_output_dir <- file.path(output_folder, current_gene_id)
      if (!dir.exists(current_output_dir)) {
        dir.create(current_output_dir)
      }
      
      print(paste("Gene ID:", current_gene_id))
      data_eQTLs_gene <-
        filter(data_eQTLs, gene_id == current_gene_id)[, c("snp",
                                                           "beta",
                                                           "varbeta",
                                                           "maf",
                                                           "position",
                                                           "pvalue",
                                                           "rsid")]
      
      merged_data <-
        inner_join(
          data_eQTLs_gene,
          data_gwas[, c("snp",
                        "beta",
                        "varbeta",
                        "position",
                        "pvalue")],
          by = c("snp", "position"),
          suffix = c("_eQTL", "_gwas")
        )
      
      print(paste("Number of SNPs in common:", dim(merged_data)[1]))
      
      ## this is the most important part of code !
      re <- coloc.abf(
        dataset1 = list(
          beta = merged_data$beta_gwas,
          varbeta = merged_data$varbeta_gwas,
          snp = merged_data$snp,
          position = merged_data$position,
          sdY = gwas_sdY,
          N = n_gwas,
          type = gwas_type
        ),
        dataset2 = list(
          beta = merged_data$beta_eQTL,
          varbeta = merged_data$varbeta_eQTL,
          snp = merged_data$snp,
          position = merged_data$position,
          N = n_eQTL,
          type = "quant",
          MAF = merged_data$maf
        )
      )
      
      ## locuscompare
      gwas_pltdata <-
        merged_data[, c("rsid", "pvalue_gwas")]
      
      
      eqtl_pltdata <-
        merged_data[, c("rsid", "pvalue_eQTL")]
      
      tmp_gwas_pltdata_path <-
        file.path(current_output_dir, "gwas.tsv")
      tmp_eqtl_pltdata_path <-
        file.path(current_output_dir, "eqtl.tsv")
      write.table(
        gwas_pltdata,
        tmp_gwas_pltdata_path,
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
      )
      write.table(
        eqtl_pltdata,
        tmp_eqtl_pltdata_path,
        col.names = T,
        row.names = F,
        sep = "\t",
        quote = F
      )
      locuscompare(
        tmp_gwas_pltdata_path,
        tmp_eqtl_pltdata_path,
        legend = F,
        title1 = "GWAS",
        title2 = "eQTL",
        genome = "hg38",
        population = "EUR",
        marker_col1 = "rsid",
        pval_col1 = "pvalue_gwas",
        marker_col2 = "rsid",
        pval_col2 = "pvalue_eQTL"
      )
      ggsave(
        file = file.path(current_output_dir, "locuscompare.png"),
        width = 10,
        height = 5
      )
      
      ## save result
      
      
      
      
      current_coloc_res <- re$summary
      

      gene_symbol <- current_gene_id
      tryCatch({
        name <-
          bitr(
            current_gene_id,
            fromType = 'ENSEMBL',
            toType = 'SYMBOL',
            OrgDb = 'org.Hs.eg.db'
          )
        gene_symbol <- name$SYMBOL[1]
      }, warning = function(w) {
        gene_symbol <- "Gene symbol not found"
      }, error = function(e) {
        gene_symbol <- "Gene symbol not found"
      })
      
      current_coloc_res$gene_id <- current_gene_id
      # print(gene_symbol)
      current_coloc_res$gene_name <- gene_symbol
      
      
      
      
      if (re$summary['PP.H4.abf'] > 0.95) {
        current_coloc_res$is_significant <- TRUE
      }
      else {
        current_coloc_res$is_significant <- FALSE
      }
      
      re_df <- re$results[, c("snp", "SNP.PP.H4")]
      colnames(re_df) <- c("snp", paste0(gene_symbol, "_PP.H4"))
      res_df <- inner_join(merged_data, re_df)
      coloc_res[[current_gene_id]] <- current_coloc_res
      coloc_snp_res[[current_gene_id]] <- res_df
      print(
        paste(
          "Gene:",
          gene_symbol,
          "PP.H4.abf",
          re$summary['PP.H4.abf']
        )
      )
      
    }
    coloc_res_df <- bind_rows(coloc_res)
    return (list(coloc_res_df = coloc_res_df, coloc_snp_res = coloc_snp_res))
  }


# runColocAnalysis(data_gwas, data_eQTLs, output_folder, n_gwas, n_eQTL, gwas_sdY, gwas_type)
