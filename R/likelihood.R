#################################################################################
Error_LLH<-function(TRFUM)   # compute error-likelihood, this function depends on the TRansition FUnction Matrix
   {                                             # error_prior is a vector not a matrix
        InOutPair<-list(); ii<-0; root_node<-numeric(); root<-0
        for ( i in 1:nrow(TRFUM))          #  i row   find each input output combination
         {
          if (sum(TRFUM[i,])>0)
            {
             ii<-ii+1
             p<-numeric(); k<-1; jj<-numeric()
             for (j in 1:nrow(TRFUM))  # j column
               if (TRFUM[i,j]!=0)
                 {
                  p[k]<-TRFUM[i,j];jj[k]<-j
                  InOutPair[[ii]]<-c(i,jj,p[k])      # note that InOutPair is a list with vector as its entries.
                  k<-k+1                             # this determines the directionality, i.e. g_k=f(g_i,g_j), i<j
                 }
            }
          if (sum(TRFUM[i,])==0)
            {
             root<-root+1
             root_node[root]<-i
            }
         }
          mismatch<-0      # count the mismatches
          if (length(InOutPair)>0)
          for (i in 1:length(InOutPair))
            {
             if (InOutPair[[i]][length(InOutPair[[i]])]<=10)
               {
                in_node<-c(InOutPair[[i]][2], InOutPair[[i]][3])
                out_node<-InOutPair[[i]][1]
                if (InOutPair[[i]][length(InOutPair[[i]])]==1)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitAnd(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==2)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-bitAnd(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==3)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitOr(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==4)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-bitOr(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==5)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitOr(1-GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==6)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitOr(GeneData[in_node[1],1:(SampleSize-1)], 1-GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==7)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitAnd(1-GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==8)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitAnd(GeneData[in_node[1],1:(SampleSize-1)], 1-GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==9)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], bitXor(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
                if (InOutPair[[i]][length(InOutPair[[i]])]==10)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-bitXor(GeneData[in_node[1],1:(SampleSize-1)], GeneData[in_node[2],1:(SampleSize-1)])))
               }
             if (InOutPair[[i]][length(InOutPair[[i]])]>10)
               {
                in_node<-InOutPair[[i]][2]
                out_node<-InOutPair[[i]][1]
                if (InOutPair[[i]][length(InOutPair[[i]])]==11)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], GeneData[in_node,1:(SampleSize-1)]))
                if (InOutPair[[i]][length(InOutPair[[i]])]==12)
                  mismatch<-mismatch+sum(bitXor(GeneData[out_node,2:SampleSize], 1-GeneData[in_node,1:(SampleSize-1)]))
               }
            }
         pseudo_count<-0.0001
         mismatch<-mismatch+pseudo_count   # mismatch may be 0
         Perror<-mismatch/(ii*SampleSize+pseudo_count); #print("Perror="); print(Perror)
         ErrorFactor<-mismatch*log(Perror)+(ii*SampleSize-mismatch+pseudo_count)*log(1-Perror)
         ErrorPrior<-(prior_para[num.node+1,1]-1)*log(Perror)+(prior_para[num.node+1,2]-1)*log(1-Perror)
         RootFactor<-numeric(); RootPrior<-numeric()
         succ_count<-numeric(); succ_prob<-numeric()
         for (i in 1:length(root_node))
           {
            succ_count[i]<-sum(GeneData[root_node[i],])+pseudo_count  # succ_count may be 0
            succ_prob[i]<-succ_count[i]/(SampleSize+pseudo_count)
            RootFactor[i]<-succ_count[i]*log(succ_prob[i])+(SampleSize+pseudo_count-succ_count[i])*log(1-succ_prob[i])

            RootPrior[i]<-(prior_para[root_node[i],1]-1)*log(succ_prob[i])+(prior_para[root_node[i],2]-1)*log(1-succ_prob[i])
           }
         if (length(InOutPair)==0)        # all nodes are root nodes
           {
            likelihood<-sum(RootFactor)
            post_para<-sum(RootFactor) + sum(RootPrior)
            log_post_model<-0
            for (i in 1:nrow(TRFUM))
              {
               nume<-lbeta(prior_para[i,1]+sum(GeneData[i,]), prior_para[i,2]+SampleSize-sum(GeneData[i,])) # lbeta=log(beta)
               deno<-lbeta(prior_para[i,1], prior_para[i,2])
               log_post_model<-log_post_model+nume-deno
              }
            log_post_model=log_post_model+length(TRFUM[TRFUM>0])*log(penalty)
            ErrorFactor<-NA; Perror<-NA; mismatch<-NA
           }
         if (length(InOutPair)>0)        # exist non root nodes
           {
            likelihood<-ErrorFactor+sum(RootFactor)      # this is the likelihood
            post_para<-ErrorFactor+sum(RootFactor) + sum(RootPrior)+ErrorPrior    # this is the posterior of (T, F, theta)
            log_post_model<-0                             # this is the posterior of (T, F)
            for (i in 1:length(root_node))
              {
               index<-root_node[i]
               nume<-lbeta(prior_para[index,1]+sum(GeneData[index,]), prior_para[index,2]+SampleSize-sum(GeneData[index,]))
               deno<-lbeta(prior_para[index,1], prior_para[index,2])
               log_post_model<-log_post_model+nume-deno
              }
            noise_nume<-lbeta(mismatch+prior_para[num.node+1,1], length(InOutPair)*SampleSize-mismatch+prior_para[num.node+1,2])
            noise_deno<-lbeta(prior_para[num.node+1,1], prior_para[num.node+1,2])
            log_post_model<-log_post_model+noise_nume-noise_deno
            log_post_model=log_post_model+length(TRFUM[TRFUM>0])*log(penalty)    # add penalty to the posterior of (T, F)
           }
         result<-list()
         result[[1]]<-c(ErrorFactor, RootFactor, likelihood, post_para, log_post_model)
         para_sample<-rep(NA,num.node)
         for (i in 1:num.node)
          if (i %in% root_node==T)
            for (j in 1:length(root_node))
              if (i==root_node[j])
                para_sample[i]<-succ_prob[j]
         result[[2]]<-c(para_sample, mismatch, Perror)
     return(result)
   }
#############################################################################################################################
#############################################################################################################################
BF1<-function(test.stat)    # model g_k=g_i and g_j
  { #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[4],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
  }
BF2<-function(test.stat)    # model g_k=complement(g_i and g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[3],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
    return(BIC.value)
 }
BF3<-function(test.stat)    # model g_k=(g_i or g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[5],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF4<-function(test.stat)    # model g_k=complement(g_i or g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[6],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF5<-function(test.stat)    # model g_k=complement(g_i) or g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[3],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF6<-function(test.stat)    # model g_k=g_i or complement( g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[5],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF7<-function(test.stat)    # model g_k=complement(g_i) and g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[6],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF8<-function(test.stat)    # model g_k=g_i and complement(g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[4],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF9<-function(test.stat)    # model g_k=g_i xor g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
    return(BIC.value)
 }
 BF10<-function(test.stat)    # model g_k=complement(g_i xor g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF11<-function(test.stat)    # model g_k=g_i
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[4],test.stat[5],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF12<-function(test.stat)    # model g_k=g_j
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[2],test.stat[3],test.stat[6],test.stat[7])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
BF13<-function(test.stat)    # model g_k=complement(g_i)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[3],test.stat[6],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
 BF14<-function(test.stat)    # model g_k=complement( g_j)
  {  #test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)
   test.stat<-test.stat+pseudo.count  # prevent come cells from being 0
   false.count<-sum(test.stat[1],test.stat[4],test.stat[5],test.stat[8])
   error.estimate=false.count/(SampleSize+pseudo.count*8)
   BIC.value=-2*(false.count*log(error.estimate)+(SampleSize+pseudo.count*8-false.count)*log(1-error.estimate))+2*log(SampleSize)
   post.data=exp(-0.5*BIC.value)
   if (false.count<=threshold)
   # return(false.count)
     return(BIC.value)
 }
#######################################################################################