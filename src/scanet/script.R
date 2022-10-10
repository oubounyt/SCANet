suppressPackageStartupMessages(library(DESeq2, quietly = T))
suppressPackageStartupMessages(library(BioNERO, quietly = T))
suppressPackageStartupMessages(library(ddpcr, quietly = T))
suppressPackageStartupMessages(library(tibble, quietly = T))
suppressPackageStartupMessages(library(R.utils, quietly = T))
suppressPackageStartupMessages(library(testthat, quietly = T))


convertor_r <- function(data_r, annotation_r){
    
    data_r = t(data_r)
    exp_filt <- PC_correction(data_r)
    final_exp <- SummarizedExperiment(data_r, colData = annotation_r)
    print(dim(final_exp))
    print("Removing outlier cells ...")
    quiet(exp_filt <- ZKfiltering(final_exp, cor_method = "pearson"))
    print(paste("Removed",dim(final_exp)[2]-dim(exp_filt)[2] ,"cells ..."))
    print(dim(exp_filt))

    return(exp_filt)
}

plot_power_r <- function(final_exp, path, network_type="unsigned", cor_method="pearson", w=400, h=400){
    

    path_ = paste(path, "/power_plot.png", sep = "")
    png(file = path_, width = w, height = h) 
    quiet(sft <- SFT_fit(final_exp, net_type=network_type, cor_method=cor_method))
    print(sft$plot)
    dev.off()

    return(sft$power)
}

co_expression_r <- function(final_exp, power, module_merging_threshold, network_type, cor_method, w=400,h=400){

    if (power == 0) {power = NULL}
    
    net <- exp2gcn(final_exp, module_merging_threshold=module_merging_threshold, net_type=network_type, 
                   SFTpower=power, cor_method=cor_method)
    names(net)

    return(net)
}

plot_dendrogram_r <- function(net, path, w, h){
    
    path_ = paste(path, "/dendrogram_plot.png", sep = "")
    png(file = path_, width = w, height = h) 
    print(plot_dendro_and_colors(net))
    dev.off()

}

plot_eigengene_network_r <- function(net, path, w, h){
    
    path_ = paste(path, "/eigengene_network_plot.png", sep = "")
    png(file = path_, width = w, height = h) 
    print(plot_eigengene_network(net))
    dev.off()

}

plot_modules_r <- function(net, path, w, h){
    
    path_ = paste(path, "/modules_plot.png", sep = "")
    png(file = path_, width = w, height = h) 
    print(plot_ngenes_per_module(net))
    dev.off()
    genes_and_modules <- net$genes_and_modules
    return(genes_and_modules)

}

module_to_annotation_cor_r <- function(final_exp, net, cor_method){
    
    pdf(file = NULL)
    invisible(capture.output(MEtrait <- module_trait_cor(exp=final_exp, MEs=net$MEs, cor_method=cor_method,
                                                         cex.lab.x=0.01, cex.lab.y = 0.01, cex.text = 0.01)))
    dev.off()
    return(MEtrait)

}

plot_module_membership_r <- function(net){

    invisible(capture.output(MEs = net$MEs))
    return(net$MEs)

}

hub_genes_r <- function(final_exp, net){

    invisible(capture.output(hubs <- get_hubs_gcn(final_exp, net)))
    return(hubs)

}

module_to_network_r <- function(net, module, co_cutoff){

    
    edges_filtered <- get_edge_list(net, module=module, filter=TRUE, method="min_cor", rcutoff = co_cutoff)
    if (dim(edges_filtered)[1] != 0){message_fit = capture_message(check_SFT(edges_filtered))
                                    m = message_fit[1]$message
                                    cat(m)}
    
    return(edges_filtered)

}

network_statistics_r <- function(edges_filtered){

    adjacency_matrix_from_data_frame <- function(dat, symmetric = TRUE,node_columns = c(1, 2)) {
        if (length(node_columns) != 2) {
            stop("length of `node_columns` must be exactly 2")
        }
            
        # remove self-interactions
        col1 <- node_columns[1]
        col2 <- node_columns[2]
        nodes1 <- as.character(dat[[col1]])
        nodes2 <- as.character(dat[[col2]])
        dat <- dat[nodes1 != nodes2, ]
        
        # create empty matrix
        proteins <- unique(c(nodes1, nodes2))
        n_proteins <- length(proteins)
        adjacency <- matrix(0, nrow = n_proteins, ncol = n_proteins,
                            dimnames = list(proteins, proteins))
        
        # identify pairwise interactions
        adjacency[as.matrix(dat[, node_columns])] <- 1
        # add symmetric interactions
        if (symmetric) {
            adjacency[as.matrix(dat[, rev(node_columns)])] <- 1
        }
        return(adjacency)
        }

    adj <- adjacency_matrix_from_data_frame(edges_filtered, node_columns = c(1, 2))
    network_statistics <- net_stats(adj_matrix = adj, net_type = "gcn", calculate_additional = TRUE)
    return(network_statistics)

}

module_consensus_powers_r <- function(final_exp, samples, nPerm){

print("print the samples")
print(samples)

assy = assay(final_exp)
assy = t(assy)
assy = as.data.frame(assy)
print(dim(assy))

type_ = colData(final_exp)
type_ = as.data.frame(type_)
print(dim(type_))

# Make index as a colmun
assy = rownames_to_column(assy)
print(dim(assy))
type_ = rownames_to_column(type_)
print(dim(type_))

# Data to DF
DF = merge(assy, type_, by='rowname' )
print(dim(DF))
#print(DF)

# Split based on cell type
data_list <-split(DF, f = DF$X__clusters__)
print(length(data_list))
print(names(data_list))

# Getting data from samples
print("")
print("getting data from samples")

samples = samples
list_SummarizedExperiment <- c()
list_names <- c()
for (i in samples) {
    a = SummarizedExperiment(data_list[c(i)])
    a = assay(a)
    a = as.data.frame(a)
    rownames(a) <- a$rowname
    annota_ = a$X__clusters__
    print(class(annota_))
    annota_ <- data.frame(annota_)
    print(class(annota_))
    a = subset(a, select = -c(rowname,X__clusters__) )
    a= t(a)
    print(dim(a))
    #list_df <- c(list_df, a)
    tmp <- SummarizedExperiment(a, colData = annota_)
    colnames(tmp)
    n = annota_$annota_
    n = n[1]
    n = as.character(n)
    print(n)
    list_SummarizedExperiment <- c(list_SummarizedExperiment, tmp)
    list_names <- c(list_names, n)
 }

names(list_SummarizedExperiment) <- list_names
print("here")
print(names(list_SummarizedExperiment[c(1)]))
print(names(list_SummarizedExperiment[c(2)]))
print("here")
#new
print("new")
cons_sft <- consensus_SFT_fit(list_SummarizedExperiment, setLabels = c(names(list_SummarizedExperiment[c(1)]), names(list_SummarizedExperiment[c(2)])), cor_method = "pearson")
powers <- cons_sft$power
print(powers)
png(file = "outs/power_plot_consensus_modules.png", width = 800, height = 400) 
print(cons_sft$plot)
dev.off()
print("end new")
#end new

return(powers)
}

module_consensus_r <- function(final_exp, samples, nPerm, powers){

print("print the samples")
print(samples)

assy = assay(final_exp)
assy = t(assy)
assy = as.data.frame(assy)
print(dim(assy))

type_ = colData(final_exp)
type_ = as.data.frame(type_)
print(dim(type_))

# Make index as a colmun
assy = rownames_to_column(assy)
print(dim(assy))
type_ = rownames_to_column(type_)
print(dim(type_))

# Data to DF
DF = merge(assy, type_, by='rowname' )
print(dim(DF))

# Split based on cell type
data_list <-split(DF, f = DF$X__clusters__)
print(length(data_list))
print(names(data_list))

# getting data from samples
print("")
print("getting data from samples")

samples = samples
list_SummarizedExperiment <- c()
list_names <- c()
for (i in samples) {
    print(i)
    a = SummarizedExperiment(data_list[c(i)])
    a = assay(a)
    a = as.data.frame(a)
    rownames(a) <- a$rowname
    annota_ = a$X__clusters__
    print(class(annota_))
    annota_ <- data.frame(annota_)
    print(class(annota_))
    a = subset(a, select = -c(rowname,X__clusters__) )
    a= t(a)
    print(dim(a))
    tmp <- SummarizedExperiment(a, colData = annota_)
    n = annota_$annota_
    n = n[1]
    n = as.character(n)
    print(n)
    list_SummarizedExperiment <- c(list_SummarizedExperiment, tmp)
    list_names <- c(list_names, n)
 }

#names(list_SummarizedExperiment) <- list_names
exp_ <- lapply(list_SummarizedExperiment, function(x) ZKfiltering(x))
exp_ <- lapply(exp_, function(x) PC_correction(x))        

consensus <- consensus_modules(list_SummarizedExperiment, power = powers, cor_method = "spearman")
print(head(consensus$genes_cmodules))
consensus_trait <- consensus_trait_cor(consensus, cor_method = "spearman")

return(consensus_trait)

}


module_conservation_r <- function(final_exp, samples, nPerm, powers){

print("print the samples")
print(samples)

assy = assay(final_exp)
assy = t(assy)
assy = as.data.frame(assy)
print(dim(assy))

type_ = colData(final_exp)
type_ = as.data.frame(type_)
print(dim(type_))

# Make index as a colmun
assy = rownames_to_column(assy)
print(dim(assy))
type_ = rownames_to_column(type_)
print(dim(type_))

# Data to DF
DF = merge(assy, type_, by='rowname' )
print(dim(DF))

# Split based on cell type
data_list <-split(DF, f = DF$X__clusters__)
print(length(data_list))
print(names(data_list))

# getting data from samples
print("")
print("getting data from samples")

samples = samples
list_SummarizedExperiment <- c()
list_names <- c()
for (i in samples) {
    print(i)
    a = SummarizedExperiment(data_list[c(i)])
    a = assay(a)
    a = as.data.frame(a)
    rownames(a) <- a$rowname
    annota_ = a$X__X__clusters____
    print(class(annota_))
    annota_ <- data.frame(annota_)
    print(class(annota_))
    a = subset(a, select = -c(rowname,X__X__clusters____) )
    a= t(a)
    print(dim(a))
    #list_df <- c(list_df, a)
    tmp <- SummarizedExperiment(a, colData = annota_)
    n = annota_$annota_
    n = n[1]
    n = as.character(n)
    print(n)
    list_SummarizedExperiment <- c(list_SummarizedExperiment, tmp)
    list_names <- c(list_names, n)
 }

names(list_SummarizedExperiment) <- list_names

#exp_ <- lapply(list_SummarizedExperiment, function(x) filter_by_variance(x))
exp_ <- lapply(list_SummarizedExperiment, function(x) ZKfiltering(x))
exp_ <- lapply(exp_, function(x) PC_correction(x))  
               
power_ortho <- lapply(list_SummarizedExperiment, SFT_fit, cor_method="pearson")
png(file = "outs/power_plot_conservation_1.png", width = 800, height = 400) 
print(power_ortho[[1]]$plot)
dev.off()
png(file = "outs/power_plot_conservation_2.png", width = 800, height = 400) 
print(power_ortho[[2]]$plot)
dev.off()

gcns <- lapply(seq_along(power_ortho), function(n) 
  exp2gcn(list_SummarizedExperiment[[n]], SFTpower = powers[[n]], 
          cor_method = "pearson"))
# Using rice as reference and maize as test
outputs = capture_message(pres <- module_preservation(list_SummarizedExperiment, ref_net = gcns[[1]], 
                                                      test_net = gcns[[2]], nPerm=nPerm, nThreads = 20), 
                                                      entrace = FALSE)
return(list(outputs, gcns[[1]]))

}

module_conservation_plot_r <- function(data){

    net = data[[2]]
    genes_and_modules <- net$genes_and_modules
    return(genes_and_modules)
}

consrv_module_to_network_r <- function(data, module, co_cutoff){

    
    png(file = "outs/consrv__scale_free_topology.png", width = 400, height = 400) 
    quiet(edges <- get_edge_list(data[[2]], module=module, filter=TRUE), all = TRUE)
    edges_filtered <- get_edge_list(data[[2]], module=module, filter=TRUE, method="min_cor", rcutoff = co_cutoff)
    dev.off()
    return(edges_filtered)

}