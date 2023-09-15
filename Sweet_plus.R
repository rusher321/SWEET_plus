Sweet_raw <- function(dat, s_method = "pearson", n_method = "pearson", balance = 0.1, PID = "all"){
  
  # calculate the weight 
  print("row is sample and colmun is feature!\n")
  corM <- cor(t(dat), method = s_method)
  pat_N <- nrow(dat)
  value <- (apply(corM, 1, sum)-1)/(pat_N-1)
  rmax <- max(value)
  rmin <- min(value)
  diff <- rmax - rmin +0.01
  value <- (value - rmin+0.01)/diff
  value <- value * balance * pat_N
  value_d <- data.frame(weight = value, row.names = rownames(dat))
  
  # calculate the coefficient of raw edge for each sample
  print("n_method need to be one of pearson, spearman, kendall and sparcc!\n")
  # remove feature with sd =0 
  dat <- dat[, apply(dat, 2, sd)!=0]
  
  if(n_method != "sparcc"){
    corN <- cor(dat, method = n_method)
  }else{
    print("if you use sparcc, the input need is count data!\n")
    corN <- SpiecEasi::sparcc(dat)
  }
  
  if(PID != "all"){
  # remove the id not in the dat
    PID <- rownames(dat)[which(rownames(dat) %in% PID)] 
  }else{
    PID <- rownames(dat)
  }
  
  if(n_method != "sparcc"){
    netlist <- lapply(PID, function(x){tmp_dat <- rbind(dat, dat[x,]); 
                           p_cor <- cor(tmp_dat, method = n_method);
                           p_cor <- value_d[x,1]*(p_cor-corN)+corN;
                           diag(p_cor) <- 0; p_cor})
  }else{
    netlist <- lapply(PID, function(x){tmp_dat <- rbind(dat, dat[x,]); 
                           p_cor <- spiec.easi::sparcc(tmp_dat);
                           p_cor <- value_d[x,1]*(p_cor-corN)+corN;
                           diag(p_cor) <- 0; p_cor})
  }
  
  # calculate the z-score of edge for each sample 
  all_value <- lapply(netlist, function(x){x[lower.tri(x)]})
  all_value <- unlist(all_value)
  mean_n <- mean(all_value)
  sd_n <- sd(all_value)
  print(mean_n)
  print(sd_n)
  outlist <- lapply(netlist, function(x){p_adj <- (x-mean_n)/sd_n;
                                        diag(p_adj) <- 0; p_adj})
  
  names(outlist) <- PID
  
  return(outlist)
    
}

net_stat <- function(cor_matrix, z_cut = 1.96){
  
  diag(cor_matrix) <- 0
  feature_name <- rownames(cor_matrix)
  # translate the correlation matrix to adjacency matrix
  net_adj <- ifelse(cor_matrix < -z_cut & cor_matrix < 0, 
                    -1, ifelse(cor_matrix > z_cut & cor_matrix > 0, 1, 0))
  # compute the degree between closeness 
  edge_pos <-   sum(net_adj[lower.tri(cor_matrix)] == 1) 
  edge_neg <-   sum(net_adj[lower.tri(cor_matrix)] == -1)
  
  net_adj <- adj2igraph(Matrix::Matrix(net_adj, sparse = TRUE))
  degree_n <-    data.frame(degree = igraph::degree(net_adj), row.names = feature_name)
  E(net_adj)$weight <- abs(E(net_adj)$weight) # translate weight to postive 
  betwenness_n <- data.frame(betweeness = igraph::betweenness(net_adj), row.names= feature_name)
  closeness_n <- data.frame(closeness = igraph::closeness(net_adj), row.names = feature_name)
  out <- list(edge_pos, edge_neg, degree_n, betwenness_n, 
              closeness_n)
  names(out) <- c("edge_pos", "edge_net", "degree", "betwenness", "closeness")
  
  return(out)
}

net_stat_com <- function(netlist, z_cut = 1.96){
  
  all_res <- lapply(netlist, net_stat, z_cut = z_cut)
  pID <- names(netlist)
  edge_num <- data.frame(pos_edge = sapply(all_res, function(x){x[[1]]}),
                         neg_edge = sapply(all_res, function(x){x[[2]]}), row.names = pID)
  
  degree_stat <- do.call("cbind", lapply(all_res, function(x){x[[3]]}))
  betw_stat <- do.call("cbind", lapply(all_res, function(x){x[[4]]}))
  clos_stat <- do.call("cbind", lapply(all_res, function(x){x[[5]]}))
  
  colnames(degree_stat) <- colnames(betw_stat) <- colnames(clos_stat) <- pID
  
  out <- list(edge_num, degree_stat, betw_stat, clos_stat)
  names(out) <- c("edge_number", "degree", "betweeness", "closeness")
  return(out)
  
}


Sweet_plus <- function(dat1, dat2, weight = "mean", s_method1 = "pearson", 
                        s_method2 = "pearson", n_method = "pearson", balance = 0.1, 
                        PID = "all"){
  # make sure consistence of the id 
  id <- intersect(rownames(dat1), rownames(dat2))
  dat1 <- dat1[id, ]
  dat2 <- dat2[id, ]
  
  # calculate the weight 
  print("row is sample and colmun is feature!\n")
  corM <- cor(t(dat1), method = s_method1)
  pat_N <- nrow(dat1)
  value <- (apply(corM, 1, sum)-1)/(pat_N-1)
  rmax <- max(value)
  rmin <- min(value)
  diff <- rmax - rmin +0.01
  value <- (value - rmin+0.01)/diff
  value <- value * balance * pat_N
  value_d1 <- data.frame(weight = value, row.names = id)
  
  corM <- cor(t(dat2), method = s_method2)
  pat_N <- nrow(dat2)
  value <- (apply(corM, 1, sum)-1)/(pat_N-1)
  rmax <- max(value)
  rmin <- min(value)
  diff <- rmax - rmin +0.01
  value <- (value - rmin+0.01)/diff
  value <- value * balance * pat_N
  value_d2 <- data.frame(weight = value, row.names = id)
  
  value_d <- data.frame(weight = (value_d1$weight+value_d2$weight)/2,
                         row.names = id)
  # calculate the coefficient of raw edge for each sample
  print("n_method need to be one of pearson, spearman, kendall and sparcc!\n")
  # remove feature with sd =0 
  dat1 <- dat1[, apply(dat1, 2, sd)!=0]
  dat2 <- dat2[, apply(dat2, 2, sd)!=0]
  
  corN <- cor(dat1, dat2, method = n_method)
  
  if(PID != "all"){
    # remove the id not in the data
    PID <- id[which(id %in% PID)] 
  }else{
    PID <- id
  }
  
  if(weight == "mean"){
    netlist <- lapply(PID, function(x){tmp_dat1 <- rbind(dat1, dat1[x,]); 
                                     tmp_dat2 <- rbind(dat2, dat2[x,]); 
                                     p_cor <- cor(tmp_dat1, tmp_dat2, method = n_method);
                                     p_cor <- value_d[x,1]*(p_cor-corN)+corN;
                                     p_cor})
  }else if(weight == "matrix1"){
    netlist <- lapply(PID, function(x){tmp_dat1 <- rbind(dat1, dat1[x,]); 
    tmp_dat2 <- rbind(dat2, dat2[x,]); 
    p_cor <- cor(tmp_dat1, tmp_dat2, method = n_method);
    p_cor <- value_d1[x,1]*(p_cor-corN)+corN;
    p_cor})    
  }else{
    netlist <- lapply(PID, function(x){tmp_dat1 <- rbind(dat1, dat1[x,]); 
    tmp_dat2 <- rbind(dat2, dat2[x,]); 
    p_cor <- cor(tmp_dat1, tmp_dat2, method = n_method);
    p_cor <- value_d2[x,1]*(p_cor-corN)+corN;
    p_cor})
    
  }
  # calculate the z-score of edge for each sample 
  all_value <- unlist(netlist)
  mean_n <- mean(all_value)
  sd_n <- sd(all_value)
  outlist <- lapply(netlist, function(x){p_adj <- (x-mean_n)/sd_n; p_adj})
  names(outlist) <- PID
  
  return(outlist) 

}


net_stat_plus <- function(cor_matrix, row = T, z_cut = 1.96){
  
  feature_name <- rownames(cor_matrix)
  feature_name2 <- colnames(cor_matrix)
  # translate the correlation matrix to adjacency matrix
  net_adj <- ifelse(cor_matrix < -z_cut & cor_matrix < 0, 
                    -1, ifelse(cor_matrix > z_cut & cor_matrix > 0, 1, 0))
  # compute the degree 
  edge_pos <-   sum(net_adj == 1) 
  edge_neg <-   sum(net_adj == -1)
  
  # get degree of feature in matrix1 
  var1_num <-  apply(net_adj, 1, function(x){sum(abs(x))})
  var1_num <- data.frame(var_num = var1_num, row.names = feature_name)
  # get degree of feature in matrix2 
  var2_num <-  apply(net_adj, 2, function(x){sum(abs(x))})
  var2_num <- data.frame(var_num = var2_num, row.names = feature_name2)
  
  out <- list(edge_pos, edge_neg, var1_num, var2_num)
  names(out) <- c("edge_pos", "edge_net", "degree_d1", "degree_d2")
  
  return(out)
  
}

net_stat_com_plus <- function(netlist, z_cut = 1.96){
  
  all_res <- lapply(netlist, net_stat_plus, z_cut = z_cut)
  pID <- names(netlist)
  edge_num <- data.frame(pos_edge = sapply(all_res, function(x){x[[1]]}),
                         neg_edge = sapply(all_res, function(x){x[[2]]}), row.names = pID)
  
  degree_stat1 <- do.call("cbind", lapply(all_res, function(x){x[[3]]}))
  degree_stat2 <- do.call("cbind", lapply(all_res, function(x){x[[4]]}))
  
  colnames(degree_stat1) <- colnames(degree_stat2) <-  pID
  
  out <- list(edge_num, degree_stat1, degree_stat2)
  names(out) <- c("edge_number", "degree_m1", "degree_m2")
  
  return(out)
  
}



