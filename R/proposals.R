 ProposalConstruction<-function(GeneData, SampleSize)
 {
gene.data=GeneData
num.node<-nrow(gene.data)
error.prop<-0.4; pseudo.count<-0.01
SampleSize<-SampleSize+pseudo.count*8; threshold<-SampleSize*error.prop
candidate.prior<-list()
kk<-1
for (i in 1: nrow(gene.data))
  for (j in 1: nrow(gene.data))
   if (j!=i)
   for (k in 1: nrow(gene.data))
    if (k!=i && k!=j)
      {  # print(i)
       test.result<-list()  # save test results for all possible relations
       gene.triplet<-rbind(rbind(gene.data[i,1:(SampleSize-1)], gene.data[j,1:(SampleSize-1)]),gene.data[k,2:SampleSize])
        c000<-0; c001<-0; c010<-0; c011<-0; c100<-0; c101<-0; c110<-0; c111<-0
        for (ii in 1: (ncol(gene.triplet)-1))   # counts in each cell
     {
      if (gene.triplet[1,ii]==0 && gene.triplet[2,ii]==0 && gene.triplet[3,(ii+1)]==0)
         c000<-c000+1
      if (gene.triplet[1,ii]==0 && gene.triplet[2,ii]==0 && gene.triplet[3,(ii+1)]==1)
         c001<-c001+1
      if (gene.triplet[1,ii]==0 && gene.triplet[2,ii]==1 && gene.triplet[3,(ii+1)]==0)
         c010<-c010+1
      if (gene.triplet[1,ii]==0 && gene.triplet[2,ii]==1 && gene.triplet[3,(ii+1)]==1)
         c011<-c011+1
      if (gene.triplet[1,ii]==1 && gene.triplet[2,ii]==0 && gene.triplet[3,(ii+1)]==0)
         c100<-c100+1
      if (gene.triplet[1,ii]==1 && gene.triplet[2,ii]==0 && gene.triplet[3,(ii+1)]==1)
         c101<-c101+1
      if (gene.triplet[1,ii]==1 && gene.triplet[2,ii]==1 && gene.triplet[3,(ii+1)]==0)
         c110<-c110+1
      if (gene.triplet[1,ii]==1 && gene.triplet[2,ii]==1 && gene.triplet[3,(ii+1)]==1)
         c111<-c111+1
      }
       test.stat<-c(c000, c001, c010, c011, c100, c101, c110, c111)   #  generate random sample
       test.result[[1]]<-c(i,j, k, 1, BF1(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[2]]<-c(i,j, k, 2, BF2(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[3]]<-c(i,j, k, 3, BF3(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[4]]<-c(i,j, k, 4, BF4(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[5]]<-c(i,j, k, 5, BF5(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[6]]<-c(i,j, k, 6, BF6(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[7]]<-c(i,j, k, 7, BF7(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[8]]<-c(i,j, k, 8, BF8(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[9]]<-c(i,j, k, 9, BF9(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[10]]<-c(i,j, k, 10, BF10(test.stat, pseudo.count, SampleSize, threshold))
       test.result[[11]]<-c(i,j, k, 11, BF11(test.stat, pseudo.count, SampleSize, threshold)) # model g_k=g_i
       test.result[[12]]<-c(i,j, k, 12, BF12(test.stat, pseudo.count, SampleSize, threshold)) # model g_k=g_j
       test.result[[13]]<-c(i,j, k, 13, BF13(test.stat, pseudo.count, SampleSize, threshold)) # model g_k=complement(g_i)
       test.result[[14]]<-c(i,j, k, 14, BF14(test.stat, pseudo.count, SampleSize, threshold))  # model g_k=complement(g_j)
       # save the most likely of all 14 relations by their false counts.
       # *******************here it may filter out some true ones due to noise and model uncertainty ****************
       miscount<-numeric(); jj<-1
       for (ii in 1:length(test.result))
         if (length(test.result[[ii]])==5)
          {
            miscount[jj]<-test.result[[ii]][5]
            jj<-jj+1
          }

         if (length(miscount)>0)
         {
            min.miscount<-min(miscount)
           for (ii in 1:length(test.result))
           if ( length(test.result[[ii]])==5 && test.result[[ii]][5]==min.miscount)
             {
              candidate.prior[[kk]]<-test.result[[ii]]
              kk<-kk+1
             }
         }   # for each pariwise genes or gene triplet, only the most likely one is output
    }
candidate<-matrix(nrow=length(candidate.prior), ncol=length(candidate.prior[[1]]))
for ( i in 1:length(candidate.prior))
  candidate[i,]<-candidate.prior[[i]]
order.candidate<-candidate[order(candidate[,3]),]   #order by output variables
##################################   # transform BIC into proportion
 num.node<-max(order.candidate[,3])
 triplet<-list();j1<-1; pairwise<-list(); j2<-1
  for (j in 1:nrow(order.candidate))
   {
   if (order.candidate[j,4]<=10)
     {
       triplet[[j1]]<-order.candidate[j,]
       j1<-j1+1
     }
   if (order.candidate[j,4]>10)
     {
      if (order.candidate[j,4]==11)
      {
        pairwise[[j2]]<-c(order.candidate[j,1], order.candidate[j,3], 11, order.candidate[j,5])
        j2<-j2+1
      }
       if (order.candidate[j,4]==12)
      {
        pairwise[[j2]]<-c(order.candidate[j,2], order.candidate[j,3], 11, order.candidate[j,5])  # model 11, g_k=g_i
        j2<-j2+1
      }
       if (order.candidate[j,4]==13)
      {
        pairwise[[j2]]<-c(order.candidate[j,1], order.candidate[j,3], 12, order.candidate[j,5])
        j2<-j2+1
      }
      if (order.candidate[j,4]==14)
      {
        pairwise[[j2]]<-c(order.candidate[j,2], order.candidate[j,3], 12, order.candidate[j,5])# model 12: complement relation
        j2<-j2+1
      }
     }
   }
candidate.triplet<-matrix(nrow=length(triplet),ncol=5); weighted.triplet<-matrix(nrow=length(triplet),ncol=6)
candidate.pairwise<-matrix(nrow=length(pairwise),ncol=4)
for (i in 1:length(triplet))
  {
   candidate.triplet[i,]<-triplet[[i]]
   for (j in 1:5)
    weighted.triplet[i,j]=candidate.triplet[i,j]
   weighted.triplet[i,6]=candidate.triplet[i,5]
  }
for (i in 1:length(pairwise))
  candidate.pairwise[i,]<-pairwise[[i]]
 unique.pairwise<-unique(candidate.pairwise)
 weighted.pairwise<-matrix(nrow=length(unique.pairwise),ncol=5)
 weighted.pairwise=cbind(unique.pairwise,unique.pairwise[,4])
constant<-0.001
for (i in 1:num.node)
 if (i %in% candidate.triplet[,3]==T)
 for (j in 1:nrow(candidate.triplet))
   if (candidate.triplet[j,3]==i)
       {
        trip.matrix<-matrix(ncol=5)
        trip.matrix<-candidate.triplet[candidate.triplet[,3]==i,1:5]
        trip.matrix<-data.matrix(trip.matrix)  # avoid pair.matrix to be "numeric"
        if (ncol(trip.matrix)==1)
         trip.matrix<-t(trip.matrix)
        score<-1/trip.matrix[,5]          # use reciprocal of miscunt as weight
        prop<-(score+constant)/(sum(score+constant))
        trip.matrix[,5]<-prop
        weighted.triplet[weighted.triplet[,3]==i,]<-cbind(trip.matrix, weighted.triplet[weighted.triplet[,3]==i,6])
        }

for (i in 1:num.node)
  if (i %in% unique.pairwise[,2]==T)
  for (j in 1:nrow(unique.pairwise))
   if (unique.pairwise[j,2]==i)
       {
        pair.matrix<-matrix(ncol=4)
        pair.matrix<-unique.pairwise[unique.pairwise[,2]==i,1:4]
        pair.matrix<-data.matrix(pair.matrix)  # avoid pair.matrix to be "numeric"
        if (ncol(pair.matrix)==1)
        pair.matrix<-t(pair.matrix)
        score<-1/(pair.matrix[,4]+constant)      # note here
        prop<-score/(sum(score))
        pair.matrix[,4]<-prop
        weighted.pairwise[weighted.pairwise[,2]==i,]<-cbind(pair.matrix, weighted.pairwise[weighted.pairwise[,2]==i,5])
        }
CandidateTriplet=weighted.triplet
CandidatePairwise=weighted.pairwise
Candidate=list()
Candidate[[1]]=CandidateTriplet
Candidate[[2]]=CandidatePairwise
return(Candidate)
 }
############
#' @importFrom stats runif
ConstructInitial <- function(Candidate, num.node)
{
prior.triplet=Candidate[[1]]
prior.pairwise=Candidate[[2]]
trans_func_matrix<-matrix(0,nrow=num.node,ncol=num.node)
prop_trans_func_matrix<-matrix(0,nrow=num.node, ncol=num.node)
prop_incid_matrix<-matrix(0, nrow=num.node, ncol=num.node)
prop_ances_matrix<-matrix(0, nrow=num.node, ncol=num.node)
incid_matrix<-matrix(0,nrow=num.node,ncol=num.node)
ratio<-0.5
for (i in 1:num.node)
  {
      prop_trans_func_matrix<-trans_func_matrix
      prop_incid_matrix<-incid_matrix
   #################################
      pairwise.prior.set<-matrix()        # clear variables     # determine candidate parents for  node i
      pairwise.prior.pare<-matrix()
      pairwise.prior.set<-prior.pairwise[prior.pairwise[,2]==i,]    # priors for node update_order[k]
      pairwise.prior.set<-data.matrix(pairwise.prior.set)
      if (ncol(pairwise.prior.set)==1)
       {
         pairwise.prior.set<-t(pairwise.prior.set)
         pairwise.prior.pare<-t(data.matrix(pairwise.prior.set[,1]))
       }
      if (ncol(pairwise.prior.set)>1 && nrow(pairwise.prior.set)>0)
        pairwise.prior.pare<-pairwise.prior.set[,1]

      triplet.prior.set<-matrix()
      triplet.prior.pare<-matrix()
      triplet.prior.set<-prior.triplet[prior.triplet[,3]==i,]
      triplet.prior.set<-data.matrix(triplet.prior.set)
      if (ncol(triplet.prior.set)>1 && nrow(triplet.prior.set)>0)
        triplet.prior.pare<-triplet.prior.set[,1:2]
      if (ncol(triplet.prior.set)==1)
        {
          triplet.prior.set<-t(triplet.prior.set)
          triplet.prior.pare<-t(data.matrix(triplet.prior.set[,1:2]))
         }
      prop_prob<-runif(1); aa<-0
      ##############################################
      if (prop_prob>=ratio)   # consider two parents
       {
        if (nrow(triplet.prior.set)>=1)
                 {
                    candidate.pare.set<-list()
                    for (jj in 1:nrow(triplet.prior.set))
                      if (length(triplet.prior.set[jj,1:2])==2)
                        {
                         aa<-aa+1
                         candidate.pare.set[[aa]]<-triplet.prior.set[jj,]
                        }
                    if (length(candidate.pare.set)>0)
                        {
                         score<-numeric()
                         for (jj in 1:length(candidate.pare.set))
                           score[jj]<-candidate.pare.set[[jj]][5]
                         bbb<-matrix(nrow=length(candidate.pare.set),ncol=6)
                         for (jjj in 1:length(candidate.pare.set))
                           bbb[jjj,]<-candidate.pare.set[[jjj]]

                         max.score.candidate<-bbb[bbb[,5]==max(score),]
                         if (is.numeric(max.score.candidate)==T)
                           add_parent=max.score.candidate
                         if (is.matrix(max.score.candidate)==T)
                         if (nrow(max.score.candidate)>1)
                           add_parent<-max.score.candidate[sample.int(nrow(max.score.candidate),1),]

                         add_two_parent<-c(add_parent[1],add_parent[2])
                      }
                 }
           if (aa>0)
             {
              func_order<-add_parent[4]
              prop_trans_func_matrix[i, add_two_parent[1]]<-func_order
              prop_trans_func_matrix[i, add_two_parent[2]]<-func_order
              for (ii in 1:nrow(prop_trans_func_matrix))
                for (jj in 1:ncol(prop_trans_func_matrix))
                 if (prop_trans_func_matrix[ii,jj]>0)
                   prop_incid_matrix[ii,jj]<-1
              prop_ances_matrix<-update_ancestor_matrix(prop_incid_matrix)
              if (check_ances_matrix(prop_ances_matrix)==0)
               {
                trans_func_matrix[i,]<-prop_trans_func_matrix[i,]
                incid_matrix[i,]<-prop_incid_matrix[i,]
               }
             }
          }
      ##########################################
      if (prop_prob<ratio) # consider one parent
       if (nrow(pairwise.prior.set)>0)
         {
            candidate.pare.set<-pairwise.prior.set
            if (is.numeric(candidate.pare.set)==T)
              add_parent<-candidate.pare.set
            if (is.matrix(candidate.pare.set)==T)
              add_parent<-candidate.pare.set[candidate.pare.set[,4]==max(candidate.pare.set[,4]),]
            if (is.matrix(add_parent)==T)
              if ( nrow(add_parent)>1)
                add_parent<-add_parent[sample.int(nrow(add_parent),1),]

            add_one_parent<-add_parent[1]

           func_order<-add_parent[3]
           prop_trans_func_matrix[i, add_one_parent]<-func_order
           for (ii in 1:nrow(prop_trans_func_matrix))
             for (jj in 1:ncol(prop_trans_func_matrix))
               if (prop_trans_func_matrix[ii,jj]>0)
                 prop_incid_matrix[ii,jj]<-1
           prop_ances_matrix<-update_ancestor_matrix(prop_incid_matrix)
           if (check_ances_matrix(prop_ances_matrix)==0)
            {
             trans_func_matrix[i,]<-prop_trans_func_matrix[i,]
             incid_matrix[i,]<-prop_incid_matrix[i,]
            }
          }
  }
  return(trans_func_matrix)
}
#########################################################################################
Prop_Trans_Func_Matrix<-function(prop_incid_matrix)      # based on incidence matrix, define transition function matrix
  {
       prop_trans_func_matrix<-prop_incid_matrix; jj<-numeric()
             for (i in 1:nrow(prop_incid_matrix))
               {
                 if (sum(prop_incid_matrix[i,])==1)    #pairwise relation
                   for (j in 1: nrow(prop_incid_matrix))
                     if (prop_incid_matrix[i,j]==1)
                       prop_trans_func_matrix[i,j]<-10+sample.int(2,1,replace=FALSE)

                 if (sum(prop_incid_matrix[i,])==2)    # triplet relation
                  {     j1<-1
                   for (j in 1:ncol(prop_incid_matrix))
                     if (prop_incid_matrix[i,j]==1)
                       {
                        jj[j1]<-j
                        j1<-j1+1
                       }
                    prop_trans_func_matrix[i,jj[1]]<-sample.int(10,1,replace=FALSE)
                    prop_trans_func_matrix[i,jj[2]]<-prop_trans_func_matrix[i,jj[1]]
                   }
                 }
      return(prop_trans_func_matrix)
  }
#################################################################################