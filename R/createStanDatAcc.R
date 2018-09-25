createStanDatAcc<-function(d=NULL, acc=NULL,form=NULL){
  
  subj <- as.integer(factor(d$subj))
  N_subj <- length(unique(subj))
  item <- as.integer(factor(d$item))
  N_items <- length(unique(item))
  X <- unname(model.matrix(form, d))  
  attr(X, which="assign") <- NULL
  
  stanDatAcc <- list(N = nrow(X),           
                  P = ncol(X),              
                  n_u = ncol(X),             
                  n_w = ncol(X),            
                  X = X,                     
                  Z_u = X,                 
                  Z_w = X,                   
                  J = N_subj,                
                  K = N_items,
                  acc = acc,                    
                  subj = subj,
                  item = item)
  stanDatAcc
}
