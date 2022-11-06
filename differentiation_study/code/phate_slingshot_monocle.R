
setwd("./")
source("./seurat_fx.thesis.R")

## source ~/anaconda3/bin/activate r-4.0

library(Seurat)
library(dplyr)
library(reshape2)
library(tidyverse)
library(gtools)
library(furrr)
library(readr)
library(tibble)
library(tidyr)
library(stringr)
library(monocle) # version 2.20
library(slingshot)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(tibble)
library(psych)
library(stats)
library(cartography)
library(colorRamps)
library(scales)
library(colorspace)
library(rcartocolor)
library(destiny)
library(phateR)
library(plotly)
library(unikn)
args <- commandArgs(trailingOnly=TRUE)
set.seed(42)



make_colors_II = function(N) {
        # https://hneth.github.io/unikn/articles/colors.html
        pal = seecol(pal = usecol(c(pal_unikn_pref), n = N ))
        return (pal)
        }




# setting I 
# mk dir and define filename
# ---------------------------------- # 

outdir <- paste0("~/project/slingshot_out/revision/")
if(!dir.exists(outdir)) dir.create(outdir)
rds_path <- "~/project/data/"
rds <- "project_data.RDS"
cell_type <- "Annotation_Level_1"



# Read
# ---------------------------------- #
R <- readRDS(paste0( rds_path, rds ))
include_celltype <- c("Fib1", "Fib2", "Fib3")
R$Annotation_Level_1 %>% gsub("/", ".", .) %>% as.character() -> R$Annotation_Level_1
Idents(R) <- "Annotation_Level_1"
R$cell_type <- Idents(R) %>% as.character()
R@assays$ECMRegulators=NULL



# Colors
# ---------------------------------- #
palette_celltype = brewer.pal(n = length(unique(R$cell_type)), name = "Set3")
real_colors <- palette_celltype
names(real_colors) <- unique(R$cell_type) %>% mixedsort()
umap_colors <- real_colors





# Subset
# ---------------------------------- #
R <- subset(R, cell_type %in% include_celltype )
root_cell <- "Fib1"
Idents(R) <- R$cell_type %>% as.character()


# colors used for slingshot
# ---------------------------------- #
palette = plasma(100)
palette_celltype = brewer.pal(n = length(unique(R$cell_type)), name = "Set3")

C <- palette_celltype[as.factor(R$cell_type)]
names(C) <- as.factor(R$cell_type)

color <- unique(C)
names(color) <- unique(names(C))


# run Diffusion map by destiny
# ---------------------------------- #
norm <- R[["RNA"]]@data
dm <- DiffusionMap(as.matrix(t(norm)), n_pcs = 30)
tmp <- data.frame(DC1 = eigenvectors(dm)[, 1], DC2 = eigenvectors(dm)[, 2])
tmp %>% as.data.frame() %>% write.table(paste0(outdir, "tmp.embedding.txt"), sep="\t", row.names = TRUE, col.names = TRUE)
colnames(tmp) <- c("dm_1", "dm_2")
R[["diffusion_map"]] <- CreateDimReducObject(embeddings = as.matrix(tmp), key = "dm_", assay = DefaultAssay(R))



## run Phate
## ---------------------------------- #
seurat_data <- as.data.frame(norm)
phate_data_input <- t(seurat_data)

for ( n in c(2,3)) {
	for ( k in c(5,10,15,20)) {
			d = 40
			p = 10
			phate_output <- phate(phate_data_input, ndim = n, decay=d, knn = k, npca=p,  mds.solver = "smacof") 
			R[[paste0("phateR.", as.character(d), ".", as.character(n), ".", as.character(k))]] <- CreateDimReducObject(embeddings = phate_output$embedding, key = "PHATE_", assay = DefaultAssay(R))
			}}

R_celltype <- R[[]] %>% dplyr::select(cell_type) %>%
                        dplyr::mutate("cell_id" = rownames(.))

R_DFMAP <- Embeddings(object = R, reduction = "diffusion_map") %>% data.frame() %>%
                                dplyr::mutate("cell_id" = rownames(.)) %>%
                                left_join(., R_celltype, by="cell_id")

R_UMAP <- Embeddings(object = R, reduction = "umap") %>% data.frame() %>%
                                dplyr::mutate("cell_id" = rownames(.)) %>%
                                left_join(., R_celltype, by="cell_id")


##
make_umap_p <- function(df, color_column, color_list, identifier) {

	colnames(df)[1:2] <- c("dimension_1", "dimension_2")
	pdf(file = paste0(outdir, identifier, '.colored.by.', color_column, '.rev.knn.pdf'), height=5, width=6.5)
        print(ggplot(df, aes(x=dimension_1, y=dimension_2, group = eval(parse(text=color_column)), colour = eval(parse(text=color_column)))) +
                  geom_point(size=0.8, alpha=0.8) +
                  theme_light() + theme_classic() + theme(legend.title=element_blank())+
                  scale_colour_manual(values=color_list))
        dev.off()
        }



##
make_umap_p(R_DFMAP, "cell_type", umap_colors[include_celltype], "dfmap")
make_umap_p(R_UMAP, "cell_type", umap_colors[include_celltype], "umap")
saveRDS(R, file = paste0(rds_path, rds, ".", args[2], '.revision.knn.RDS'))

## 2d 
for ( d in c(40)) {
	for ( k in c(5,10,15,20) ) {
	R_PHATE <- Embeddings(object = R, reduction = paste0("phateR.", as.character(d), ".2.", as.character(k) )) %>% data.frame() %>%
					dplyr::mutate("cell_id" = rownames(.)) %>%
					left_join(., R_celltype, by="cell_id")
	print(make_umap_p(R_PHATE, "cell_type", umap_colors[include_celltype], paste0("phateR.", as.character(d), ".2.", as.character(k))))

	R_PHATE <- Embeddings(object = R, reduction = paste0("phateR.", as.character(d), ".3.", as.character(k) )) %>% data.frame() %>%
                                        dplyr::mutate("cell_id" = rownames(.)) %>%
                                        left_join(., R_celltype, by="cell_id")
	print(make_umap_p(R_PHATE, "cell_type", umap_colors[include_celltype], paste0("phateR.", as.character(d), ".3.", as.character(k))))
	}
	}


##
R_ALL <- R



# run Slingshots 
# ---------------------------------- #

run_pseudotime <- function(your_R, reduction_type, root_cell ) {
	
	R <- your_R
	start.clus <- root_cell
	reduction = reduction_type
	sds = slingshot(Embeddings(R, reduction), clusterLabels = Idents(R), start.clus = start.clus )
	R@tools[['slingshot']] = SlingshotDataSet(sds)
	pseudotime = slingPseudotime(sds)
	sds_all = slingshot(Embeddings(R_ALL, reduction), clusterLabels = Idents(R_ALL), start.clus = start.clus )
	
	
	
	# Plot slingshot curves
	# ---------------------------------- # 
	pseudotime %>% as.data.frame() %>%
                        write.table(paste0(outdir, "slingshot.pseudotime.rev.knn.txt"), sep="\t", row.names = TRUE, col.names = TRUE)
	curves = colnames(pseudotime)

	# colors
	C1 = umap_colors[mixedsort(include_celltype)]
	C = C1[R$cell_type]

	pdf(file = paste0(outdir, 'slingshot_curves.', reduction, '.rev.knn.pdf'), height=5, width = 5.2)
	print(plot(sds$reducedDim, col = C, pch = 16, cex = 0.5) + ## 1 
                lines(SlingshotDataSet(sds), lwd = 2, col = 'black') 
		)
	dev.off()
	
	
	
	# Plot slingshot curves
	# ---------------------------------- #
	sds_all$reducedDim %>% rownames() -> cell_id
	umap_colors[R$cell_type] -> cell_color
	names(cell_color) <- colnames(R)



	# Plot slingshot curve II by cell type
	# ---------------------------------- #
	pseudotime_orig <- pseudotime
	sds_orig <- sds
	pdf(file = paste0(outdir, 'slingshot_curves.separate.', reduction, '.rev.knn.pdf'), width = 7.2)
	par(mfrow = c(2, 3))
	for ( c_num in seq(1, length(curves))) {
                        
			sds <- sds[,curves[c_num]]
                        pseudotime = slingPseudotime(sds)
                        print(plot(sds$reducedDim, col = C, pch = 16, cex = 0.5, main = curves[c_num] ) +
                                lines(SlingshotDataSet(sds), linInd = c_num, lwd = 2, col = 'black'))
                     
			sds <- sds_orig
			pseudotime <- pseudotime_orig
			}

	dev.off()




	R_TIME <- pseudotime_orig %>% as.data.frame() %>% dplyr::mutate("cell_id" = rownames(.))
	R_META <- R[["cell_type"]] %>% dplyr::mutate("cell_id" = rownames(.))
	##R_META <- R_ALL[["cell_type"]] %>% dplyr::mutate("cell_id" = rownames(.))
	##R_UMAP <- Embeddings(object = R_ALL, reduction = reduction) %>% data.frame() %>% 
	R_UMAP <- Embeddings(object = R, reduction = reduction) %>% data.frame() %>% 
			dplyr::mutate("cell_id" = rownames(.)) %>% 
				left_join(., R_META, by="cell_id") %>% 
				left_join(., R_TIME, by="cell_id") 


	
	colnames(R_UMAP)[1:2] <- c("dim_1", "dim_2")
	
	
	make_umap <- function(df, color_column, color_list) {
		pdf(file = paste0(outdir, 'slingshot_curves.separate.', reduction, '.colored.by.', color_column, '.NEW.rev.knn.pdf'), width = 7)
		print(ggplot(df, aes(x=dim_1, y=dim_2, group = eval(parse(text=color_column)), colour = eval(parse(text=color_column)))) +
	          	geom_point(size=0.8, alpha=0.8) + 
			theme_light() + theme_classic() +
			scale_colour_manual(values=color_list) + 
			theme(legend.position = "none", panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
			scale_alpha(guide = 'none'))# to remove extra legend 
		dev.off()
		}

	make_umap(R_UMAP, "cell_type", umap_colors)





	make_umap_c <- function(df, color_column) {
		pdf(file = paste0(outdir, 'slingshot_curves.separate.', reduction, '.colored.by.', color_column, '.NEW.rev.knn.pdf'), width = 7)
        	print(ggplot(R_UMAP, aes(x=dim_1, y=dim_2, colour = eval(parse(text=color_column)) )) +
                	  geom_point(size=0.8, alpha=0.8) +
                	  theme_light() + theme_classic() +
			  scale_colour_viridis_c(na.value="#D3D3D3", option = "C") +
                	  theme(legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
                	  scale_alpha(guide = 'none'))# to remove extra legend 
        	dev.off()
		}  

	make_umap_c(R_UMAP, "Lineage1")
	for (i in c("Lineage2", "Lineage3", "Lineage4", "Lineage5", "Lineage6")) {
		if ( i %in% colnames(R_UMAP)) {
			make_umap_c(R_UMAP, i)
		}
	}




	# Plot slingshot curve II by cell type
	# ---------------------------------- #
	pseudotime_orig <- pseudotime
	sds_orig <- sds
	pdf(file = paste0(outdir, 'slingshot_curves.separate.', reduction, '.colored.by.time.rev.knn.pdf'), width = 7.2)
	par(mfrow = c(2, 3))
	for ( c_num in seq(1, length(curves))) {

			##sds_all <- sds
                        sds <- sds[,curves[c_num]]
			##pseudotime_all = slingPseudotime(sds_all)
                        pseudotime = slingPseudotime(sds)
                        colors = palette[cut(pseudotime[,1], breaks = 100)]
			colors[ is.na(colors) ] = "#D3D3D3"
                        print(plot(sds$reducedDim, col = colors, pch = 16, cex = 0.5, main = curves[c_num] ) +
                                lines(SlingshotDataSet(sds), linInd = c_num, lwd = 2, col = 'black'))

                        sds <- sds_orig
                        pseudotime <- pseudotime_orig
                        }

	dev.off()
	
	
	
	
	# Add pseudotimes (arclength) to meta data for visualisation
	# ---------------------------------- #
	# > R[[c("Lineage1","Lineage2")]] %>% as.data.frame() -> b
	# > pseudotime %>% as.data.frame() -> a
	# > all.equal(a,b)
	# [1] TRUE
	# ---------------------------------- #
	for ( curve in curves ) {
		pseudotime_sub <- pseudotime[colnames(R),curve]
		R <- AddMetaData(object = R,
                         metadata = pseudotime_sub,
                         col.name = curve
                         )
         	}
	
	
	
	
	# Condition density along pseudotime
	# ---------------------------------- #
	make_density <- function(your_obj, curve, color_list) {
		pdf(file = paste0(outdir, 'density_condition.', reduction, ".", curve, '.rev.knn.pdf'), width = 7, height = 5)
		df <- data.frame(your_obj[["cell_type"]], your_obj[[curve]]) 
		colnames(df) <- c("cell_type", "Lineage")
		na.omit(df) -> df
		your_obj$cell_type <- factor(your_obj$cell_type, levels=mixedsort(unique(your_obj$cell_type)))
		p <- ggplot(df, aes(x=Lineage, fill=cell_type)) +
			geom_density(alpha=0.8) + theme_classic()+
			#scale_fill_manual(values=color_list)	
			scale_fill_manual(values=C) 
		print(p)
		dev.off()
		}

	R$cell_type <- factor(R$cell_type, levels=mixedsort(unique(R$cell_type)))
	make_density(R, "Lineage1", umap_colors[mixedsort(include_celltype)])
	for (i in c("Lineage2", "Lineage3", "Lineage4", "Lineage5", "Lineage6")) {
		if ( i %in% colnames(R_UMAP)) {
			make_density(R, i, umap_colors[mixedsort(include_celltype)])
		}
	}


	return (R)
	}




pairs <- list( c(40,2), c(40,3))
ks <- c(5,10,15,20)

phateR_slingshot_set <- c()
for ( n in c(1:length(pairs)) ) {
	for (k in ks) {
		phateR_slingshot_set[[paste0(n, ".", k)]] <- run_pseudotime(R, paste0("phateR.", as.character(unlist(pairs[n])[1]), ".", as.character(unlist(pairs[n])[2]), ".", as.character(k)) , root_cell)
	}}




run_monocle_DE <- function(your_obj, your_column, VGAM_opt, VGAM_opt_name, identifier) {


	# subset based on lineage 2 
	# ---------------------------------- #
	L <- your_obj[[your_column]] %>% deframe()
	names(L) <- rownames(your_obj[[your_column]])
	L[!is.na(L)] %>% names() -> L2_cell
	#R$Lineage2[!is.na(R$Lineage2)] %>% names() -> L2_cell
	your_obj[, L2_cell] -> new_R 
	Idents(new_R) <- new_R$cell_type %>% as.character()

	
	
		
	# Transfer into cds
	# ---------------------------------- #
	cds <- as.CellDataSet(new_R)
	# Estimate size factor
	cds <- estimateSizeFactors(cds)
	cds <- estimateDispersions(cds)
	
	
	
	# call Monocle2
	# ---------------------------------- #
	# install https://www.bioconductor.org/packages/3.14/bioc/src/contrib/Archive/monocle/
	# http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories
	# https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/monocle2.html#monocle2-process
	# ---------------------------------- #
	# select superset of feature genes as genes expressed in at least 5% of all the cells.
	cds <- detectGenes(cds, min_expr = 0.1)
	fData(cds)$use_for_ordering <- fData(cds)$num_cells_expressed > 0.05 * ncol(cds)
	cds_genes <-  row.names(subset(fData(cds),num_cells_expressed >= 10))
	
	
	
	# get genes used for ordering cells 
	# ---------------------------------- #
	# while removing batch by using fullModelFormulaStr 
	# https://www.biostars.org/p/316204/
	# ---------------------------------- #
	# colnames(pData(cds)
	# df = 1 : a linear fit 
	# df = 2 : affords a little nonlinearity
	# df = 3 : VGAM
	# http://www2.uaem.mx/r-mirror/web/packages/VGAM/vignettes/categoricalVGAM.rev.knn.pdf
	# ---------------------------------- #
	clustering_DEG_genes <- differentialGeneTest(cds[cds_genes,],
						fullModelFormulaStr = '~cell_type',
					     	cores = 10)



	cds_ordering_df <- clustering_DEG_genes %>% filter(qval < 0.01 & use_for_ordering == TRUE) %>% arrange(qval)
	cds_ordering_df[1:1000, ] %>% select(gene_short_name) %>% deframe() %>% as.character() -> cds_ordering_genes
	cds_ordering_genes[ is.na(cds_ordering_genes)==FALSE ] -> cds_ordering_genes 

	# Clustering Genes by Pseudotemporal Expression Pattern by Monocle2
	# ---------------------------------- #
	pData(cds)[[your_column]] -> pData(cds)$Pseudotime
	diff_test_res <- differentialGeneTest(cds[cds_ordering_genes,],
	                                fullModelFormulaStr = VGAM_opt, 
	                                cores = 10)
	
	diff_test_res %>% filter(qval <0.01) %>% arrange(qval) %>% write.table(., paste0(outdir, identifier, "slingshot.DEgenes.", your_column, ".", VGAM_opt_name, ".rev.txt"), quote = FALSE, sep="\t", row.names = TRUE, col.names = TRUE)


	

	# make plots
	# ---------------------------------- #
	DE_N_Set <- c(30) # c(20, 30, 50)

	for ( DE_N in DE_N_Set ) { 
	
		diff_test_res %>% filter(qval <0.01) %>% arrange(qval) %>% top_n(., DE_N, wt=-qval) %>% rownames() -> sig_gene_names


	        hm <- get_pseudotime_matrix(cds[sig_gene_names,],  
	                                    cluster_rows = TRUE,
	                                    hclust_method = "ward.D",
	                                    num_clusters = 6,
	                                    hmcols = NULL,
	                                    add_annotation_row = NULL,
        	                            add_annotation_col = NULL,
        	                            show_rownames = FALSE,
        	                            use_gene_short_name = TRUE,
        	                            norm_method = "log",
        	                            scale_max=3,
        	                            scale_min=-3,
        	                            trend_formula = VGAM_opt, 
        	                            return_heatmap=TRUE,
        	                            cores=1)


		bks = c(seq(min(hm), 0, length.out=ceiling(200/2) + 1),
	                seq(max(hm)/200, max(hm),length.out=floor(200/2)))
		
	        #my_color = colorRampPalette(c("#00008B", "#F5F5F5", "#cd5c5c"))(length(bks))
		#my_color3 = viridis(length(bks))
		my_color4 = plasma(length(bks))
		my_color5 = colorRampPalette(rev(rcartocolor::carto_pal(7, "Sunset")))(length(bks))
		my_color6 = colorRampPalette(rcartocolor::carto_pal(7, "ag_Sunset"))(length(bks))
		my_color7 = colorRampPalette(rev(rcartocolor::carto_pal(7, "SunsetDark")))(length(bks))
			
		my_color_set <- list(my_color7)
		#my_color_set <- list(my_color4, my_color5, my_color6, my_color7)
		my_color_name <- c("SunsetDark")
		#my_color_name <- c("plasma", "Sunset", "ag_Sunset", "SunsetDark")
		


		# cluster and re-order rows
		#ALL_HCS <- c( "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
		ALL_HCS <- c("ward.D")

		for ( sub_color in seq(1,length(my_color_set)))  {
        		for ( ALL_HC in c(ALL_HCS) ) {	
		
				pdf(file = paste0(outdir, identifier, 'slingshot.top', DE_N , '.DEgenes.', your_column, '.', ALL_HC, '.', my_color_name[sub_color], '.', VGAM_opt_name, '.rev.knn.pdf'), width = 4, height = 6)
				print(monocle::plot_pseudotime_heatmap(cds[sig_gene_names,],
	        	        	        ##add_annotation_col = your_column,
						cluster_rows = TRUE,
						trend_formula = VGAM_opt,
						hclust_method = ALL_HC, 
						num_clusters = 1,
						hmcols = my_color_set[sub_color][[1]],
						scale_max = 3, 
						scale_min = -3,
	        	        	        cores = 1,
	        	        	        show_rownames = T,
						return_heatmap = FALSE))
				dev.off()
				}

			}

		}
	
	

        colors = palette[cut(pData(cds)$Pseudotime, breaks = 100)]
	phenoData(cds)[["color"]] <- colors
	GENE_OF_INTEREST <- c("Gene1", "Gene2", "Gene3", "Gene4")
	pdf(file = paste0(outdir, identifier, 'slingshot.gene_of_interest.pseudotime.', your_column,  '.plasma.rev.knn.pdf'), width = 4, height = 10)
	print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = your_column ) +         
	      				scale_color_viridis(option = "C") 
	      				#scale_color_viridis()
	      				)
	dev.off()

	pdf(file = paste0(outdir, identifier, 'slingshot.gene_of_interest.pseudotime.', your_column,  '.inferno.rev.knn.pdf'), width = 4, height = 10)
        print(plot_genes_in_pseudotime(cds[GENE_OF_INTEREST,], color_by = your_column ) +
                                        scale_color_viridis(option = "B")
                                        #scale_color_viridis()
                                        )

	dev.off()


	}



run_monocle_DE(phateR_slingshot_set[["1.10"]], "Lineage1", "~sm.ns(Pseudotime, df=1)", "linear", "phateR.40.2.10.Lineage1")
run_monocle_DE(phateR_slingshot_set[["1.10"]], "Lineage2", "~sm.ns(Pseudotime, df=1)", "linear", "phateR.40.2.10.Lineage2")
run_monocle_DE(phateR_slingshot_set[["1.10"]], "Lineage3", "~sm.ns(Pseudotime, df=1)", "linear", "phateR.40.2.10.Lineage3")



