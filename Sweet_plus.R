Sweet_sin <- function(dat, s_method = "pearson", n_method = "pearson", balance = 0.1, PID = "all"){
  
  # calculate the weight 
  print("row is sample and colmun is feature!\n")
  corM <- cor(t(dat), s_method)
  pat_N <- nrow(dat)
  value <- (apply(corM, 1, sum)-1)/(pat_N-1)
  value <- value * balance * pat_N
  value_d <- data.frame(weight = value, row.names = rownames(dat))
  
  # calculate the coefficient of raw edge for each sample
  print("n_method need to be one of pearson, spearman, kendall and sparcc!\n")
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
                           diag(p_cor) <- 0; p_cor})
  }else{
    netlist <- lapply(PID, function(x){tmp_dat <- rbind(dat, dat[x,]); 
                           p_cor <- spiec.easi::sparcc(tmp_dat);
                           diag(p_cor) <- 0; p_cor})
  }
  
  # calculate the z-score of edge for each sample 
  all_value <- lapply(netlist, function(x){x[lower.tri(x)]})
  all_value <- unlist(all_value)
  mean_n <- mean(all_value)
  sd_n <- sd(all_value)
  outlist <- lapply(netlist, function(x){p_adj <- (x-mean_n)/sd_n;
                                        diag(p_adj) <- 0; p_adj})
  
  names(outlist) <- PID
  
  return(outlist)
    
}

Sweet_plus <- function(dat1, dat2, method1, method2, method3, balance = 0.1){
  
  
  
  
}