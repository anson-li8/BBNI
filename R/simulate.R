#############################################################################################################################
GenerateNetwork<-function(num.node)
{
 loop=1
 while (loop!=0)
  {
   tent_incid_matrix<-matrix(0,nrow=num.node, ncol=num.node)   # define random incidence matrix, ancestor matrix, transition matrix
   tent_trans_matrix<-matrix(0,nrow=num.node, ncol=num.node)
   for (i in 1:num.node)
     {
      u1=runif(1); uu=4       # uu determines the network complexity
      if (u1>1/10 && u1<=uu/10)
       {
        position=sample(seq(num.node)[-i],1); tent_incid_matrix[i,position]=1
        u2=runif(1)
        if (u2>0.5)
          tent_trans_matrix[i,position]=11
        if (u2<0.5)
          tent_trans_matrix[i,position]=12
       }
      if (u1>uu/10)
       {
        position=sample(seq(num.node)[-i],2); tent_incid_matrix[i,position[1]]=1; tent_incid_matrix[i,position[2]]=1
        func=sample(10,1); tent_trans_matrix[i,position[1]]=func; tent_trans_matrix[i,position[2]]=func
       }
     }
   tent_ances_matrix=update.ancestor_matrix(tent_incid_matrix)
   loop=check.ances.matrix(tent_ances_matrix)
  }
  return(tent_trans_matrix)
}
#############################################################################################################################
#############################################################################################################################
GenerateSample<-function(trans_matrix)
{
 node_ances=matrix(nrow=num.node, ncol=2)
 GeneData=matrix(0, nrow=num.node, ncol=SampleSize)
 incid_matrix<-trans_matrix
    for (i in 1:nrow(trans_matrix))
     for (j in 1:ncol(trans_matrix))
      if (trans_matrix[i,j]>0)
       incid_matrix[i,j]<-1
 ances_matrix=update.ancestor_matrix(incid_matrix)
 for (i in 1:num.node)
  {
   node_ances[i,1]=i
   node_ances[i,2]=sum(ances_matrix[i,])
  }
 node_ances=node_ances[order(node_ances[,2]),]
 for (i in 1:nrow(node_ances))
   {
    if (node_ances[i,2]==0)
      GeneData[node_ances[i,1],]=rbinom(SampleSize,1,prob=para[node_ances[i,1]])
    if (node_ances[i,2]!=0)
      {
       parent<-numeric(); ii=1
       for (j in 1:ncol(incid_matrix))
         if (incid_matrix[node_ances[i,1], j]!=0)
           {
            parent[ii]=j; ii=ii+1
           }
       func=trans_matrix[node_ances[i,1],parent[1]]
       if (func==1)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitAnd(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==2)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-bitAnd(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==3)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitOr(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==4)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-bitOr(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==5)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitOr(1-GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==6)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitOr(GeneData[parent[1],1:(SampleSize-1)], 1-GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==7)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitAnd(1-GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==8)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitAnd(GeneData[parent[1],1:(SampleSize-1)], 1-GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==9)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(bitXor(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==10)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-bitXor(GeneData[parent[1],1:(SampleSize-1)], GeneData[parent[2],1:(SampleSize-1)]), error[node_ances[i,1],1:(SampleSize-1)])
       if (func==11)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(GeneData[parent[1],1:(SampleSize-1)], error[node_ances[i,1],1:(SampleSize-1)])
       if (func==12)
         GeneData[node_ances[i,1],2:SampleSize]=bitXor(1-GeneData[parent[1],1:(SampleSize-1)], error[node_ances[i,1],1:(SampleSize-1)])
      }
   }
 return(GeneData)
}
#############################################################################################################################