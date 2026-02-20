library(sleuth)
sleuth_table <- data.frame(
  sample = c("SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"),
  condition = c("2dpi", "6dpi", "2dpi", "6dpi"),
  path = c(
    "quant_reads/SRR5660030",
    "quant_reads/SRR5660033",
    "quant_reads/SRR5660044",
    "quant_reads/SRR5660045"
  ),
  stringsAsFactors = FALSE
)
so = sleuth_prep(sleuth_table)
so = sleuth_fit(so, ~condition, 'full')
so = sleuth_fit(so, ~1, 'reduced')
so = sleuth_lrt(so, 'reduced', 'full')

library(dplyr)
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 
#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 
columns <- c('target_id','test_stat','pval','qval')
subset_sig <- sleuth_significant[,columns]
#write FDR < 0.05 transcripts to file
write.table(subset_sig, file="sleuth_out.txt",quote = FALSE,row.names = FALSE, sep = "\t")
#head(dplyr::select(sleuth_significant, target_id, pval, qval), n=10)