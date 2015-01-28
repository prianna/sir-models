N = length(WT06P4V1$V1)
KidType1Run4V1 <- vector()
KidType2Run4V1 <- vector()
KidBothRun4V1 <- vector()
PType1Run4V1 <- vector()
PType2Run4V1 <- vector()
PBothRun4V1 <- vector()
NumV1Run4V1 <- vector()
NumV2Run4V1 <- vector()
for(i in 1:N) 
{
  if(WT06P4V1$V6[i] < 3)
  {
    if( WT06P4V1$V1[i] > 0 && WT06P4V1$V2[i] > 0 )
    {
      PBothRun4V1 <- c(PBothRun4V1, WT06P4V1$V7[i])
    }
    else if(WT06P4V1$V1[i] == 0 && WT06P4V1$V2[i] > 0)
    {
      PType2Run4V1 <- c(PType2Run4V1, WT06P4V1$V7[i])
    }
    else ##( WT06P4V1$V1[i] > 0 && WT06P4V1$V2[i] == 0 )
    {
      PType1Run4V1 <- c(PType1Run4V1, WT06P4V1$V7[i])
    }
  }
  if(WT06P4V1$V6[i] > 0)
  {
    if( WT06P4V1$V1[i] > 0)
    {
      NumV1Run4V1 <- c(NumV1Run4V1, WT06P4V1$V1[i])
    }
    if( WT06P4V1$V2[i] > 0)
    {
      NumV2Run4V1 <- c(NumV2Run4V1, WT06P4V1$V2[i])
    }
    if( WT06P4V1$V1[i] > 0 && WT06P4V1$V2[i] > 0 )
    {
      KidBothRun4V1 <- c(KidBothRun4V1, WT06P4V1$V1[i]+WT06P4V1$V2[i])
    }
    else if(WT06P4V1$V1[i] == 0 && WT06P4V1$V2[i] > 0)
    {
      KidType2Run4V1 <-c(KidType2Run4V1, WT06P4V1$V2[i])
    }
    else if( WT06P4V1$V1[i] > 0 && WT06P4V1$V2[i] == 0 )
    {
      KidType1Run4V1 <-c(KidType1Run4V1, WT06P4V1$V1[i])
    }
    else
    {
      ##
    }
  }
}
sumPBothRun4V1 <- sum(PBothRun4V1)
sumPType1Run4V1 <- sum(PType1Run4V1)
sumPType2Run4V1 <- sum(PType2Run4V1)