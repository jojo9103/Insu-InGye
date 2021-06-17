argv1<-commandArgs(trailingOnly = T)

total_result<-NULL
for(f in argv1){

lm_result<-read.table(f,sep='\t',header=T,stringsAsFactors = F)
print(dim(lm_result))
drug_name<-strsplit(f,'/')[[1]]
drug_name<-drug_name[length(drug_name)]
drug_name<-strsplit(drug_name,'_Selected_CL.txt')[[1]][1]
drug_name<-strsplit(drug_name,'ssgsea_')[[1]][2]
print(drug_name)

rownames(lm_result)=lm_result[,'TF']
lm_result<-lm_result[,-1]
if(ncol(lm_result)>2){
lm_p<-lm_result[,grep('lm_p',colnames(lm_result))]
lm<-lm_result[,-grep('lm_p',colnames(lm_result))]

colnames(lm_p)<-gsub('lm',drug_name,colnames(lm_p))
colnames(lm)<-gsub('lm',drug_name,colnames(lm))

sub_result<-NULL
for(n in 1:ncol(lm)){
  n1=colnames(lm)[n]
  n1=strsplit(n1,'_')[[1]]
  n1=n1[length(n1)]
  sub_summary=cbind.data.frame(TF=rownames(lm),lm=lm[,n],lm_p=lm_p[,n],drug_name=drug_name,dose=n1,drug_id=paste0(drug_name,'_',n1))
  sub_result<-rbind.data.frame(sub_result,sub_summary)
}

}else{
  n1=colnames(lm_result)[1]
  n1=strsplit(n1,'_')[[1]]
  n1=n1[length(n1)]
  sub_result=cbind.data.frame(TF=rownames(lm_result),lm_result,drug_name=drug_name,dose=n1,drug_id=paste0(drug_name,'_',n1))
  colnames(sub_result)=c('TF','lm','lm_p','drug_name','dose','drug_id')
}


total_result<-rbind.data.frame(total_result,sub_result)

}
write.table(total_result,'/home/jojo9103/LINCS/LINCS/20210415_New_GSVA/lm_result_Summary.txt',sep='\t',row.names = F,quote = F)
