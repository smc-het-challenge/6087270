library(RWeka)
args<-commandArgs(TRUE)
vcfdat = read.table(args[1],sep='\t',comment.char='#')
batdat = read.table(args[2],sep='\t',comment.char='#')
datacol=as.integer(args[3]) + 9
namecol=9
tumour_stat = data.frame(do.call(rbind, strsplit(as.vector(vcfdat[,datacol]), split = ":", fixed = TRUE)))
colnames(tumour_stat) = strsplit(as.vector(unique(vcfdat[,namecol])),':')[[1]]
TDP = as.integer(as.vector(tumour_stat[,'DP']))
TAD1 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',1)))
TAD2 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',2)))
BQ <- as.integer(as.vector(tumour_stat[,'BQ']))

if(datacol==11){datacol=10}else{datacol=11}
tumour_stat = data.frame(do.call(rbind, strsplit(as.vector(vcfdat[,datacol]), split = ":", fixed = TRUE)))
colnames(tumour_stat) = strsplit(as.vector(unique(vcfdat[,namecol])),':')[[1]]
NDP <- as.integer(as.vector(tumour_stat[,'DP']))
NAD1 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',1)))
NAD2 <-  as.integer(unlist(lapply(strsplit(as.vector(tumour_stat[,'AD']),','),'[[',2)))

chrvcf <-as.vector(vcfdat[,1])
posvcf<-as.integer(as.vector(vcfdat[,2]))
vcfref<- as.vector(vcfdat[,4])
vcfalt<- as.vector(vcfdat[,5])

chrbat <- as.vector(batdat[2:length(batdat[,1]),1])
strbat <- as.integer(as.vector(batdat[2:length(batdat[,1]),2]))
endbat <- as.integer(as.vector(batdat[2:length(batdat[,1]),3]))
BAFbat <- as.double(as.vector(batdat[2:length(batdat[,1]),4]))
LRbat <- as.double(as.vector(batdat[2:length(batdat[,1]),6]))
BAF=NULL
LR=NULL
A=NULL
T=NULL
C=NULL
G=NULL
tA=NULL
tT=NULL
tC=NULL
tG=NULL
BAF[length(TDP)]=NA
LR[length(TDP)]=NA
newx=NULL
X=TAD2/TDP
tempchr=NULL
temppos=NULL
l = 1
for(i in 1 : length(TDP)){
  for(j in 1:(length(batdat[,1])-1)){
    
    if(vcfref[i]=="A"){
      A[i]=1
    }else{
      A[i]=0
    }
    if(vcfref[i]=="T"){
      T[i]=1
    }else{
      T[i]=0
    }
    if(vcfref[i]=="C"){
      C[i]=1
    }else{
      C[i]=0
    }
    if(vcfref[i]=="G"){
      G[i]=1
    }else{
      G[i]=0
    }
    if(vcfalt[i]=="A"){
      tA[i]=1
    }else{
      tA[i]=0
    }
    if(vcfalt[i]=="T"){
      tT[i]=1
    }else{
      tT[i]=0
    }
    if(vcfalt[i]=="C"){
      tC[i]=1
    }else{
      tC[i]=0
    }
    if(vcfalt[i]=="G"){
      tG[i]=1
    }else{
      tG[i]=0
    }
    if(chrvcf[i] == chrbat[j] && posvcf[i]>=strbat[j] && posvcf[i]<=endbat[j]){
      BAF[i]=BAFbat[j]
      LR[i]=LRbat[j]
    }
    
  }
}
TFA= TAD2/TDP
NFA=NAD2/NDP
CNV=TDP/NDP

BBAF=BAF
if(length(TDP)>5000){
  ttt= data.frame(TFA,NAD2,TDP,NDP,CNV,BQ,A,T,G,C,tA,tT,tC,tG)
  ttt=as.data.frame(scale(ttt))
  LR=scale(LR)
  BAF=scale(BAF)
  LR=as.factor(LR)
  BAF = as.factor(BAF)
  ttt=data.frame(ttt,BAF,LR)
}else{
  ttt= data.frame(TFA,NAD2,TDP,NDP,CNV,BQ,C,tC)
  ttt=as.data.frame(scale(ttt))
  LR=scale(LR)
  BAF=scale(BAF)
  LR=as.factor(LR)
  BAF = as.factor(BAF)
  ttt=data.frame(ttt,BAF,LR)
}
write.arff(ttt,file="tempfile.arff")

out = system("java -cp /usr/share/java/weka.jar -mx1024m weka.clusterers.EM -t tempfile.arff -I 1000 -N -1 -X 10 -max -1 -ll-cv 5.0E-7 -ll-iter 5.0E-7 -M 5.0E-7 -K 1000 -num-slots 4 -S 100",TRUE)
clusterNum = 0
 
 for(c in 1 : length(out)) {
   if(out[c] == "Clustered Instances") {
     c = c + 1
     for(p in c : length(out)) {
       if(out[p] != "") {
         clusterNum = clusterNum + 1;
       }else {
         break
       }
    }
     break
   }
}
runSKM= data.frame(TFA,BQ,NAD2)
AAAA=SimpleKMeans(runSKM,Weka_control(N=clusterNum,S=100))

clust=AAAA$class_ids
fuckkkkk=as.data.frame(AAAA$class_ids)
kk0=1
mk0=1
m0 = NULL
B0 = NULL
cellularity=NULL
cell=1
resultNum=NULL
for(j in 1 : (clusterNum)){
  m0 = NULL
  B0 = NULL
  mk0=1
  kk0=1
  for(i in 1 : length(TDP)){
    if(clust[i]==(j-1)){
      m0[kk0]=TAD2[i]/(TAD1[i]+TAD2[i])
      
      if(is.na(BBAF[i])){
        
      }else{
        B0[mk0]=BBAF[i]
        mk0=mk0+1
      }
      kk0=kk0+1
    }
  }
  cellularity[j]=(mean(m0)/(1.5-mean(B0)))
  resultNum[j]=kk0-1
}
mcell=min(cellularity)
for(i in 1 : clusterNum ){
  if(cellularity[i] == mcell){
    cellularity[i]=0
    tempi = i
  }
}
tcell=cellularity[clusterNum]
tresu=resultNum[length(resultNum)]
cellularity[length(cellularity)]=cellularity[tempi]
resultNum[length(resultNum)]=resultNum[tempi]
cellularity[tempi]=tcell
resultNum[tempi]=tresu
clust = clust + 1
#resultNum[tempi]=resultNum[1]
#cellularity[tempi]=cellularity[1]
#cellularity[1]=0
#resultNum[1]=0
#anscel=cellularity[2:length(cellularity)]
#ansnum=resultNum[2:length(resultNum)]
#anscluster=NULL
#clustercount=1
#for(i in 1:length(AAAA$class_ids)){
#  if(AAAA$class_ids[i]!=tempi-1){
#    anscluster[clustercount]=AAAA$class_ids[i]
#    clustercount=clustercount+1
#  }
#}
if(tempi!=clusterNum){
  for(i in 1:length(clust)){
    if(clust[i]==tempi){
      clust[i]=clusterNum
    }else if (clust[i]==clusterNum){
      clust[i]=tempi
    }
  }
}
  
write.table(max(cellularity),"subchallenge1A.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(clusterNum-1,"subchallenge1B.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(cbind(1:clusterNum,resultNum,cellularity),"subchallenge1C.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(clust,"subchallenge2A.txt",row.names=F,col.names=F,quote=F,sep="\t")
