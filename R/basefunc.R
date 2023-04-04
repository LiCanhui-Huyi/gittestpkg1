append_extra_info <- function(fileOfRNA,fileOfATAC,fileOfpeakMap){
  #读取数据文件
  rna <- data.table::fread(fileOfRNA,nrows = 1,header = FALSE)  #只读取第一行，包含样本名，用于筛选相同的样本
  atac <-data.table::fread(fileOfATAC,nrows = 1,header = FALSE)
  atacmap <<- data.table::fread(fileOfpeakMap,sep = "\t")   # 全局，用于peak注释
  
  #获取两个矩阵中相同的列名（样本）
  rnacol <- rna %in% atac; rnacol[1] <- TRUE; rnacol <-which(rnacol==TRUE)
  RNA_in_ATAC <- as.data.frame(data.table::fread(fileOfRNA,select = rnacol)) #用fread读取文件，数据格式不是简单的dataframe
  RNA_in_ATAC <- RNA_in_ATAC[stringr::str_sort(colnames(RNA_in_ATAC))] #进行排序
  geneid <- RNA_in_ATAC[,1];RNA_in_ATAC[,1]<- gsub("\\.\\d*","",geneid) #提取基因id
  
  #使用biomaRt获取基因ID对应的基因名
  while(! exists("my_newid")){
    try(my_mart <-biomaRt::useMart("ensembl"),silent = TRUE)#创建mart对象
    if(! exists("my_mart"))next #如果mymart对象不存在，进入下一次循环
    try(my_dataset <- biomaRt::useDataset("hsapiens_gene_ensembl",mart = my_mart) )#选择数据库    
    try(my_newid <- biomaRt::getBM(attributes = c("hgnc_symbol","ensembl_gene_id","chromosome_name",'start_position','end_position'), #想要获得的数据类型
                                   filters = "ensembl_gene_id", #提供的数据类型
                                   values = RNA_in_ATAC[,1], #提供的数据类型对应数据
                                   mart = my_dataset))
  }
  #合并处理biomaRt获得的数据框
  RNA_in_ATAC <- merge(my_newid,RNA_in_ATAC,by.x = "ensembl_gene_id",by.y = "Ensembl_ID")
  RNA_in_ATAC[,3] <- paste0("chr",RNA_in_ATAC[,3])  #染色体名加上chr，与ATAC数据一致
  RNA <<-  RNA_in_ATAC  #创建全局变量RNA存储矩阵
  
  #获取两个矩阵中相同的列名（样本）
  ataccol <- atac %in% rna;ataccol[1] <- TRUE;ataccol <- which(ataccol==TRUE)
  ATAC_in_RNA <-as.data.frame(data.table::fread(fileOfATAC,select = ataccol))
  ATAC_in_RNA <- ATAC_in_RNA[stringr::str_sort(colnames(ATAC_in_RNA))]
  ATAC <<- as.data.frame(cbind(atacmap[,3:5],ATAC_in_RNA[,-1]))
  
  
}

peak_anno <- function(atacmap){
  atacmapbed <- atacmap[,3:6] #map文件3到6列包含bed文件的内容
  write.table(atacmapbed,file = "atac.bed",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
  
  #所需人类基因注释
  txdb  <-  TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  #注释peak
  x <- ChIPseeker::annotatePeak("atac.bed", tssRegion=c(-1000, 1000), TxDb=txdb,annoDb = "org.Hs.eg.db")
  PeakAnno <<- as.data.frame(x)  #保存为数据框，便于查看
  ChIPseeker::upsetplot(x, vennpie=TRUE)    #注释可视化
}


my_cor_test <- function(pair){ #相关性检测       辅助函数    传入变量为基因和peak的位置对
  rna <- RNA[pair[1],]
  atac <- ATAC[pair[2],]
  r <- cor.test(as.numeric(rna[,-1:-5]),as.numeric(atac[,-1:-3]),method = "spearman",exact=FALSE)
  return(c(r$estimate,r$p.value))
}
 

linkage <- function(dataOfRNA,dataOfATAC,geneset,rho=0.3,p=0.1){  #该函数目前遇到基因集全部来自同一染色体时会报错
  
  
  # 创建一个空的结果数据框
  result <- data.frame(RNA_index = integer(), ATAC_index = integer())
  
  # 对于每个基因进行操作
  for (i in seq_along(geneset)) {
    gene <- geneset[i]
    # print(gene)
    # 找到该基因在RNA数据框中的行数
    gene_row <- match(gene, dataOfRNA$hgnc_symbol)
    # print(gene_row)
    
    # 如果该基因不在RNA数据框中，则跳过
    if (is.na(gene_row)) {
      next
    }
    
    # 找到该基因的位置信息
    chr <- dataOfRNA$chromosome_name[gene_row]
    start_pos <- dataOfRNA$start_position[gene_row]
    end_pos <- dataOfRNA$end_position[gene_row]
    
    # 找到所有满足条件的peak的行数
    peak_rows <- which(dataOfATAC$chrom == chr & 
                         dataOfATAC$chromStart >= start_pos - 500000 &
                         dataOfATAC$chromEnd <= end_pos + 500000)
    
    # 将结果添加到结果数据框中
    for (j in peak_rows) {
      result <- rbind(result, data.frame(RNA_index = gene_row, ATAC_index = j))
    }
  }
  
  
  
  
  cor_p <- apply(result, 1, my_cor_test) %>% t()
  loc_cor <- (abs(cor_p[,1]) > rho & cor_p[,2] < p) %>% {cbind(result[.,],cor_p[.,])} %>% as.data.frame()
  loc_cor <- cbind(RNA[loc_cor[,1],2],ATAC[loc_cor[,2],1:3],loc_cor)
  colnames(loc_cor) <- c("hgnc_symbol", "chrom", "chromStart",  "chromEnd" ,"idx_gene","idx_peak","rho","p_value")
  loc_cor <<- loc_cor
  
}


#可视化相关性
plot_gene_peak_correlation <- function(loc_cor,RNA,ATAC){
  # merged_data <- rbind(RNA[location_pair[1,1],6:ncol(RNA)], ATAC[location_pair[1,2],4:ncol(ATAC)])
  # return(merged_data)
  
  gene_expr <- RNA[loc_cor[1,5],6:ncol(RNA)] %>% as.numeric()
  peak_data <- ATAC[loc_cor[1,6],4:ncol(ATAC)] %>% as.numeric()
  
  
  p <- ggplot2::ggplot(data.frame(gene_expr, peak_data), ggplot2::aes(x = gene_expr, y = peak_data)) +
    ggplot2::geom_point() + # 添加散点图
    ggplot2::stat_smooth(method = "lm", se = FALSE) + # 添加回归线
    ggplot2::theme_classic() + # 使用白色背景并去除网格线
    ggplot2::xlab("RNA-seq") + ggplot2::ylab("ATAC-seq") + # 添加X轴和Y轴标签
    ggplot2::ggtitle(loc_cor[1,1], paste(loc_cor[1,2],":",loc_cor[1,3],"-",loc_cor[1,4])) + # 添加标题和副标题
    ggplot2::annotate("text", x = min(gene_expr), y = max(peak_data), # 在左上角添加相关系数和FDR值
                      label = paste0("Corr =",round(as.numeric(loc_cor[1,7]), 4), "\np_value =", format.pval(as.numeric(loc_cor[1,8]))),
                      hjust = 0, vjust = 1, size = 3)
  
  # 显示图形
  print(p)
}


##
library(magrittr)
rna <- "TCGA-BRCA.htseq_fpkm-uq.tsv"
atac <- "brca_brca_peak_Log2Counts_dedup.brca_brca_peak_log2counts_dedup"
peakmap <- "brca_brca_peak.probeMap"
append_extra_info(rna,atac,peakmap)
# peak_anno(atacmap)

timestart <- Sys.time()
linkage(RNA,ATAC,RNA$hgnc_symbol[1:5])
timeend <- Sys.time()
print(timeend-timestart)


plot_gene_peak_correlation(loc_cor[1,],RNA,ATAC)


# 尚未完成（顺序无关）
# 
# 6封装成包
# 7函数文档
