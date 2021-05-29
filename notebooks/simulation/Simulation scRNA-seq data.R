"Установка необходимых библиотек"

library("devtools")
devtools::install_github('YosefLab/SymSim')
library("SymSim")

"Набор необходимых генов"

data(gene_len_pool)
gene_len <- sample(gene_len_pool, 20000, replace = FALSE)

"Датасет 1 - 500 клеток, 2 батча и 3 типа клеток"

"Генерация данных и изображение кластеров"

true_counts_res <- SimulateTrueCounts(ncells_total=500, min_popsize=50, i_minpop=2, ngenes=20000, nevf=10, evf_type="discrete", n_de_evf=9, vary="s", Sigma=0.4, phyla=Phyla3(), randseed=0)
true_counts_res_dis <- true_counts_res
tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], data=log2(true_counts_res[[1]]+1), evf_type="discrete", n_pc=20, label='pop', saving = F, plotname="discrete populations (true counts)")
tsne_true_counts[[2]]

"Перерасчет транскриптов в UMI(unicular molecular identifiers"

observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], meta_cell=true_counts_res[[3]], protocol="UMI", alpha_mean=0.05, alpha_sd=0.02, gene_len=gene_len, depth_mean=5e4, depth_sd=3e3)
tsne_UMI_counts <- PlotTsne(meta=observed_counts[[2]], data=log2(observed_counts[[1]]+1), evf_type="discrete", n_pc=20, label='pop', saving = F, plotname="observed counts UMI")
tsne_UMI_counts[[2]]

"Добавление двух батч-эффектов"

observed_counts_2batches <- DivideBatches(observed_counts_res = observed_counts, nbatch = 2, batch_effect_size = 1)
tsne_batches <- PlotTsne(meta=observed_counts_2batches[[2]], data=log2(observed_counts_2batches[[1]]+1), evf_type="discrete", n_pc=20, label='batch', saving = F, plotname="observed counts in batches")
tsne_batches[[2]]

"Сохранение файла"

install.packages("anndata")
library(anndata)
install_anndata(method = "auto", conda = "auto")
expression_matrix <- observed_counts_2batches[[1]]
var <- data.frame(cell_id = observed_counts_2batches[[2]]$cellid, cell_type = observed_counts_2batches[[2]]$pop, batch = observed_counts_2batches[[2]]$batch, row.names = observed_counts_2batches[[2]]$cellid)
obs <- data.frame(gene_name = seq(1, 20000), row.names = seq(1, 20000))
ad <- AnnData(T, X = expression_matrix, var = var, obs = obs)
ad$write_h5ad("c500t3b2.h5ad")

"Датасет 2 - 15к клеток, 4 батча и 5 типов клеток"

"График с отображением кластеров"

true_counts_res_2 <- SimulateTrueCounts(ncells_total=15000, min_popsize=50, i_minpop=2, ngenes=20000, nevf=10, evf_type="discrete", n_de_evf=9, vary="s", Sigma=0.4, phyla=Phyla5(), randseed=0)
true_counts_res_dis_2 <- true_counts_res_2
tsne_true_counts_2 <- PlotTsne(meta=true_counts_res_2[[3]], data=log2(true_counts_res_2[[1]]+1), evf_type="discrete", n_pc=20, label='pop', saving = F, plotname="discrete populations (true counts)")
tsne_true_counts_2[[2]]  

"Генерация и сохранение датасета"

observed_counts_2 <- True2ObservedCounts(true_counts=true_counts_res_2[[1]], meta_cell=true_counts_res_2[[3]], protocol="UMI", alpha_mean=0.05, alpha_sd=0.02, gene_len=gene_len, depth_mean=5e4, depth_sd=3e3)
observed_counts_2batches_2 <- DivideBatches(observed_counts_res = observed_counts_2, nbatch = 4, batch_effect_size = 1)
expression_matrix_2 <- observed_counts_2batches_2[[1]]
var_2 <- data.frame(cell_id = observed_counts_2batches_2[[2]]$cellid, cell_type = observed_counts_2batches_2[[2]]$pop, batch = observed_counts_2batches_2[[2]]$batch, row.names = observed_counts_2batches_2[[2]]$cellid)
ad <- AnnData(T, X = expression_matrix_2, var = var_2, obs = obs)
ad$write_h5ad("c1500t5b4.h5ad")

"Датасет 2 - 100к клеток, 2 батча и 5 типов клеток"

"График с отображением кластеров"

true_counts_res_3 <- SimulateTrueCounts(ncells_total=100000, min_popsize=50, i_minpop=2, ngenes=20000, nevf=10, evf_type="discrete", n_de_evf=9, vary="s", Sigma=0.4, phyla=Phyla5(), randseed=0)
true_counts_res_dis_3 <- true_counts_res_3
tsne_true_counts_3 <- PlotTsne(meta=true_counts_res_3[[3]], data=log2(true_counts_res_3[[1]]+1), evf_type="discrete", n_pc=20, label='pop', saving = F, plotname="discrete populations (true counts)")
tsne_true_counts_3[[2]]

"Генерация и сохранение датасета"

observed_counts_3 <- True2ObservedCounts(true_counts=true_counts_res_3[[1]], meta_cell=true_counts_res_3[[3]], protocol="UMI", alpha_mean=0.05, alpha_sd=0.02, gene_len=gene_len, depth_mean=5e4, depth_sd=3e3)
observed_counts_2batches_3 <- DivideBatches(observed_counts_res = observed_counts_3, nbatch = 2, batch_effect_size = 1)
expression_matrix_3 <- observed_counts_2batches_3[[1]]
var_3 <- data.frame(cell_id = observed_counts_2batches_3[[2]]$cellid, cell_type = observed_counts_2batches_3[[2]]$pop, batch = observed_counts_2batches_3[[2]]$batch, row.names = observed_counts_2batches_3[[2]]$cellid)
ad <- AnnData(T, X = expression_matrix_3, var = var_3, obs = obs)
ad$write_h5ad("c10000t5b2.h5ad")