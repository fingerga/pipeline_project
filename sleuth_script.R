library(sleuth)

#creating necessary input table for sleuth with samples, conditions, and paths to reads"
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
#initialize sleuth object using sleuth_prep function from sleuth library
so = sleuth_prep(sleuth_table)
#fit a model comparing the two conditions 
so = sleuth_fit(so, ~condition, 'full')
#fit the reduced model to compare in the likelihood ratio test 
so = sleuth_fit(so, ~1, 'reduced')
#perform the likelihood ratio test for differential expression between conditions 
so = sleuth_lrt(so, 'reduced', 'full')

library(dplyr)
#extract the test results from the sleuth object 
sleuth_table = sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) 
#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant = dplyr::filter(sleuth_table, qval <= 0.05) |> dplyr::arrange(pval) 
#selecting columns to output
columns <- c('target_id','test_stat','pval','qval')
subset_sig <- sleuth_significant[,columns]
#write FDR < 0.05 transcripts to file
write.table(subset_sig, file="sleuth_out.txt",quote = FALSE,row.names = FALSE, sep = "\t")