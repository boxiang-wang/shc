q.shc <- function(x, y, qvcv=F, qv=1, qv.list=c(0.5, 1, 2, 5, 100), silence = T, 
                  lambda2.list = c(1e-4, 1e-3, 1e-2, 1e-1, 1, 5, 10), ...) {
  lambda2.list <- rev(sort(lambda2.list))
  print(lambda2.list)
  if(qvcv == F){
    ## cv error for l1 Holder classifier
    cv_l1_out <- cv.shc(x, y, qv=qv, lambda2=0, ...)
    l1_l1 <- cv_l1_out$lambda.min
    l1_cvm <- cv_l1_out$cvm.min
    ## cv error for enet Holder classifier
    mm.cvm <- Inf
    for(i in seq.int(length(lambda2.list))){
      cv_enet_out <- cv.shc(x, y, qv=qv, lambda2=lambda2.list[i], ...)
      if(mm.cvm > cv_enet_out$cvm.min){
        mm.cvm <- cv_enet_out$cvm.min
        mm.lambda <- cv_enet_out$lambda.min
        lambda2.id <- i
      }
    enet_l2 <- lambda2.list[lambda2.id]
    enet_l1 <- mm.lambda
    enet_cvm <- mm.cvm
    return(c(l1_cvm=l1_cvm, l1_l1=l1_l1, enet_cvm=enet_cvm, 
             enet_l1=enet_l1, enet_l2=enet_l2))
    }
  } else {
    res <- matrix(NA, length(qv.list), 6)
    res[,1] <- qv.list
    colnames(res) <- c('qv', 'l1.cvm', 'l1.l1', 'enet.cvm', 'l2.l1', 'l2.l2')
    for(j in seq.int(length(qv.list))){
      cv_l1_out <- cv.shc(x, y, qv=qv.list[j], lambda2=0, ...)
      res[j, 3] <- cv_l1_out$lambda.min
      res[j, 2] <- cv_l1_out$cvm.min
      ## cv error for enet Holder classifier
      mm.cvm <- Inf
      for(i in seq.int(length(lambda2.list))){
        cv_enet_out <- cv.shc(x, y, qv=qv.list[j], lambda2=lambda2.list[i], ...)
        if(mm.cvm > cv_enet_out$cvm.min){
          mm.cvm <- cv_enet_out$cvm.min
          mm.lambda <- cv_enet_out$lambda.min
          lambda2.id <- i
        }
        res[j, 6] <- lambda2.list[lambda2.id]
        res[j, 5] <- mm.lambda
        res[j, 4] <- mm.cvm
      }
      if(silence == F) {
        print(paste('q = ', qv.list[j], ' is completed.', sep=""))
      }
    }
    return(res)
  }
}  


