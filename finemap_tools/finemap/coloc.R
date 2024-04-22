library(coloc)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(locuscomparer)
coloc_analysis <- function(sumstats1,
                           sumstats2,
                           output_folder,
                           n1,
                           n2,
                           sdY1,
                           sdY2 = NULL,
                           type1 = "quant",
                           type2 = "quant",
                           num_causal_snps = 1,
                           ld_df = NULL,
                           sumstats1_suffix = "gwas",
                           sumstats2_suffix = "eQTL")   {
  if (!dir.exists(output_folder)) {
    dir.create(output_folder)
  }
  
  if (num_causal_snps > 1) {
    if (is.null(ld_df)) {
      stop("ld_df must be provided when num_causal_snps > 1")
    }
  }
  if (is.null(sdY2) && type2 == "quant") {
    if (!("maf" %in% colnames(sumstats2))) {
      stop("sdY2 must be provided when type2 is quant and maf is provided in sumstats2")
    }
  }
  
  
  if ("gene_id" %in% colnames(sumstats2)) {
    print("gene_id column found in sumstats2, using gene_id to split the data")
  }
  
  ## 由于ld_d 的index和columns 被R会修改成1.xx.a.t的形式，因此sumstats1 和2的index 默认也是这个样子，所以这里的snp，直接用rownames
  sumstats1$snpidx <- rownames(sumstats1)
  sumstats2$snpidx <- rownames(sumstats2)
  
  merged_data <-
    inner_join(
      sumstats2,
      sumstats1[, c("snpidx",
                    "beta",
                    "varbeta",
                    "position",
                    "pvalue")],
      by = c("snpidx", "position"),
      
      suffix = c(
        paste0("_", sumstats1_suffix),
        paste0("_", sumstats2_suffix)
      )
    )
  
  print(paste0("Number of SNPs in common:", dim(merged_data)[1]))
  
  ## this is the most important part of code !
  dataset1 = list(
    beta = merged_data[, paste0("beta_", sumstats1_suffix)],
    varbeta = merged_data[, paste0("varbeta_", sumstats1_suffix)],
    snp = merged_data[, "snpidx"],
    # index is the snp name with 1.xx.a.t like ld_df do
    position = merged_data$position,
    sdY = sdY1,
    N = n1,
    type = type1
  )
  
  dataset2 = list(
    beta = merged_data[, paste0("beta_", sumstats2_suffix)],
    varbeta = merged_data[, paste0("varbeta_", sumstats2_suffix)],
    snp = merged_data[, "snpidx"],
    # index is the snp name with 1.xx.a.t like ld_df do
    position = merged_data$position,
    N = n2,
    type = type2
  )


  if ("maf" %in% colnames(merged_data)) {
    dataset2$MAF = merged_data$maf
  }
  if (!(is.null(sdY2))) {
    dataset2$sdY = sdY2
  }

  if (num_causal_snps > 1) {
    print(dim(ld_df))
    print(dim(merged_data))
    dataset1$LD = as.matrix(ld_df)
    dataset2$LD = as.matrix(ld_df)
    
    save(sumstats1,
         sumstats2,
         merged_data,
         dataset1,
         dataset2,
         ld_df,
         file = "/home/xutingfeng/GIFT/data/analysis/DNAJC16/coloc/coloc_test.RData")
    
    if (!(is.null(check_dataset(dataset1)))) {
      stop("dataset1 is not valid")
    }
    if (!(is.null(check_dataset(dataset2)))) {
      stop("dataset2 is not valid")
    }
    
    ds1 <- runsusie(dataset1)
    ds2 <- runsusie(dataset2)
    if (requireNamespace("susieR", quietly = TRUE)) {
      susie.res = coloc.susie(ds1, ds2)
      print(susie.res$summary)
    }
  }
    
  else {
      print("ld_df is not provided, using coloc.abf")
      re <- coloc.abf(dataset1 = dataset1,
                      dataset2 = dataset2)
      print("Finish coloc.abf")
      print(re$summary)
      # save sensitivity plot
      png(
        filename =  file.path(output_folder, "sensitivity.png"),
        # 文件名称
        width = 10,
        # 宽
        height = 10,
        # 高
        units = "in",
        # 单位
        bg = "white",
        # 背景颜色
        res = 400
      )             # 分辨率
      sensitivity(re, "H4 > 0.95")
      dev.off()
    }
    
    
    ## locuscompare
    pltdata1 <-
      merged_data[, c("rsid", paste0("pvalue_", sumstats1_suffix))]
    
    
    pltdata2 <-
      merged_data[, c("rsid", paste0("pvalue_", sumstats2_suffix))]
    
    tmp_pltdata1_path <-
      file.path(output_folder, paste0(sumstats1_suffix, ".tsv"))
    
    tmp_pltdata2_path <-
      file.path(output_folder, paste0(sumstats2_suffix, ".tsv"))
    
    write.table(
      pltdata1,
      tmp_pltdata1_path,
      col.names = T,
      row.names = F,
      sep = "\t",
      quote = F
    )
    write.table(
      pltdata2,
      tmp_pltdata2_path,
      col.names = T,
      row.names = F,
      sep = "\t",
      quote = F
    )
    locuscompare(
      tmp_pltdata1_path,
      tmp_pltdata2_path,
      legend = F,
      title1 = sumstats1_suffix,
      title2 = sumstats2_suffix,
      genome = "hg38",
      population = "EUR",
      marker_col1 = "rsid",
      pval_col1 = paste0("pvalue_", sumstats1_suffix),
      marker_col2 = "rsid",
      pval_col2 = paste0("pvalue_", sumstats2_suffix)
    )
    ggsave(
      file = file.path(output_folder, "locuscompare.png"),
      width = 10,
      height = 5
    )
    
    ## save result
    current_coloc_res <- re$summary
    
    if (re$summary['PP.H4.abf'] > 0.95) {
      current_coloc_res$is_significant <- TRUE
    }
    else {
      current_coloc_res$is_significant <- FALSE
    }
    
    re_df <- re$results[, c("snp", "SNP.PP.H4")]
    # print(colnames(re_df))
    # print(colnames(merged_data))
    res_df <- inner_join(merged_data, re_df)
    print(colnames(res_df))
    # return (list(coloc_res = current_coloc_res, coloc_snp_res = res_df))
    return (list(coloc_res = current_coloc_res, coloc_snp_res = res_df))
    
  }







# runColocAnalysis <-
#   function(sumstats1,
#            sumstats2,
#            output_folder,
#            n1,
#            n2,
#            sdY1,
#            sdY2=NULL,
#            type1="quant",
#            type2="quant",
#            num_causal_snps=1,
#            ld_df=NULL,
#            sumstats1_suffix="gwas",
#             sumstats2_suffix="eQTL"
#            ) {

#     "
#     sumstats1: GWAS summary statistics, must have columns: snp, beta, varbeta, position, pval
#     sumstats2: eQTL summary statistics, must have columns: snp, beta, varbeta, maf, position, pval, rsid, gene_id
    
#     "
    
    

#     # sumstats1$pvalue <- sumstats1$pval  #TODO: pvalue cols => pval
#     gene_id_list <- unique(sumstats2$gene_id)
    
#     ### get a gene eQTLs
#     coloc_res <- list()
#     coloc_snp_res <- list()
    
#     for (idx in 1:length(gene_id_list)) {
#       current_gene_id <- gene_id_list[idx]
#       current_output_dir <- file.path(output_folder, current_gene_id)
#       if (!dir.exists(current_output_dir)) {
#         dir.create(current_output_dir)
#       }
      
#       print(paste0("Gene ID:", current_gene_id))
#       sumstats2_gene <-
#         filter(sumstats2, gene_id == current_gene_id)[, c("snp",
#                                                            "beta",
#                                                            "varbeta",
#                                                            "maf",
#                                                            "position",
#                                                            "pvalue",
#                                                            "rsid")]
      
#       merged_data <-
#         inner_join(
#           sumstats2_gene,
#           sumstats1[, c("snp",
#                         "beta",
#                         "varbeta",
#                         "position",
#                         "pvalue")],
#           by = c("snp", "position"),
#           suffix = c("_eQTL", "_gwas")
#         )
      
#       print(paste0("Number of SNPs in common:", dim(merged_data)[1]))
      
#       ## this is the most important part of code !
#       re <- coloc.abf(
#         dataset1 = list(
#           beta = merged_data$beta_gwas,
#           varbeta = merged_data$varbeta_gwas,
#           snp = merged_data$snp,
#           position = merged_data$position,
#           sdY = sdY1,
#           N = n1,
#           type = type1,
#         ),
#         dataset2 = list(
#           beta = merged_data$beta_eQTL,
#           varbeta = merged_data$varbeta_eQTL,
#           snp = merged_data$snp,
#           position = merged_data$position,
#           N = n2,
#           type = "quant",
#           MAF = merged_data$maf
#         )
#       )
#       # save sensitivity plot
#       png( 
#         filename =  file.path(current_output_dir, "sensitivity.png"), # 文件名称
#         width = 10,            # 宽
#         height = 10,           # 高
#         units = "in",          # 单位
#         bg = "white",          # 背景颜色
#         res = 400)             # 分辨率
#       sensitivity(my.res,"H4 > 0.95")
#       dev.off()
      
#       ## locuscompare
#       gwas_pltdata <-
#         merged_data[, c("rsid", "pvalue_gwas")]
      
      
#       eqtl_pltdata <-
#         merged_data[, c("rsid", "pvalue_eQTL")]
      
#       tmp_gwas_pltdata_path <-
#         file.path(current_output_dir, "gwas.tsv")
#       tmp_eqtl_pltdata_path <-
#         file.path(current_output_dir, "eqtl.tsv")
#       write.table(
#         gwas_pltdata,
#         tmp_gwas_pltdata_path,
#         col.names = T,
#         row.names = F,
#         sep = "\t",
#         quote = F
#       )
#       write.table(
#         eqtl_pltdata,
#         tmp_eqtl_pltdata_path,
#         col.names = T,
#         row.names = F,
#         sep = "\t",
#         quote = F
#       )
#       locuscompare(
#         tmp_gwas_pltdata_path,
#         tmp_eqtl_pltdata_path,
#         legend = F,
#         title1 = "GWAS",
#         title2 = "eQTL",
#         genome = "hg38",
#         population = "EUR",
#         marker_col1 = "rsid",
#         pval_col1 = "pvalue_gwas",
#         marker_col2 = "rsid",
#         pval_col2 = "pvalue_eQTL"
#       )
#       ggsave(
#         file = file.path(current_output_dir, "locuscompare.png"),
#         width = 10,
#         height = 5
#       )
      
#       ## save result
      
      
      
      
#       current_coloc_res <- re$summary
      

#       gene_symbol <- current_gene_id
#       tryCatch({
#         name <-
#           bitr(
#             current_gene_id,
#             fromType = 'ENSEMBL',
#             toType = 'SYMBOL',
#             OrgDb = 'org.Hs.eg.db'
#           )
#         gene_symbol <- name$SYMBOL[1]
#       }, warning = function(w) {
#         gene_symbol <- "Gene symbol not found"
#       }, error = function(e) {
#         gene_symbol <- "Gene symbol not found"
#       })
      
#       current_coloc_res$gene_id <- current_gene_id
#       # print(gene_symbol)
#       current_coloc_res$gene_name <- gene_symbol
      
      
      
      
#       if (re$summary['PP.H4.abf'] > 0.95) {
#         current_coloc_res$is_significant <- TRUE
#       }
#       else {
#         current_coloc_res$is_significant <- FALSE
#       }
      
#       re_df <- re$results[, c("snp", "SNP.PP.H4")]
#       colnames(re_df) <- c("snp", paste00(gene_symbol, "_PP.H4"))
#       res_df <- inner_join(merged_data, re_df)
#       coloc_res[[current_gene_id]] <- current_coloc_res
#       coloc_snp_res[[current_gene_id]] <- res_df
#       print(
#         paste0(
#           "Gene:",
#           gene_symbol,
#           "PP.H4.abf",
#           re$summary['PP.H4.abf']
#         )
#       )
      
#     }
#     coloc_res_df <- bind_rows(coloc_res)
#     return (list(coloc_res_df = coloc_res_df, coloc_snp_res = coloc_snp_res))
#   }
# }
