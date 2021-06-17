library(GSVA)
library(openxlsx)
library(stringr)
library(pheatmap)
# argv1[1] : pathway file, argv1[2]  input dataset path, argv1[3] output path

argv1<-commandArgs(trailingOnly = T)

geneset<-read.xlsx(argv1[1],sheet=1) # sheet 2 == symbol gene

geneset_list<-list()
for(l in 1:nrow(geneset)){
  l1=geneset[l,8:ncol(geneset)]
  l1=l1[!is.na(l1)]
  geneset_list[[geneset$TF[l]]]=as.character(l1)
}

filename1=list.files(argv1[2])

DMSO=read.table(gzfile(paste0(argv1[2],'DMSO_Selected_CL.txt.gz')),sep='\t',header=T,stringsAsFactors = F)
rownames(DMSO)<-DMSO[,3]
ford_mat=DMSO[,1:3]
DMSO=DMSO[,4:ncol(DMSO)]

for(drug_f in filename1){
  i_n=strsplit(drug_f,'_')[[1]][1]
  if(i_n!='DMSO'){
    drug_f1=read.table(gzfile(paste0(argv1[2],drug_f)),sep='\t',header = T,stringsAsFactors = F,quote = '')
    rownames(drug_f1)=drug_f1[,3]
    if(ncol(drug_f1)>3){
      drug_f2=cbind.data.frame(DMSO,drug_f1[,4:ncol(drug_f1)])
      colnames(drug_f2)=c(colnames(DMSO),colnames(drug_f1)[4:ncol(drug_f1)])
      drug_f1<-as.matrix(drug_f2)
      ssgsea_r<-gsva(drug_f1,geneset_list,method='ssgsea',kcdf="Gaussian",abs.ranking=FALSE,
                     min.sz=2,max.sz=Inf,parallel.sz=10,mx.diff=TRUE,
                     ssgsea.norm=TRUE,
                     verbose=TRUE)
      # t-test
      ssgsea_r_ctl<-ssgsea_r[,grep('DMSO',colnames(ssgsea_r))]
      ssgsea_r_case<-as.data.frame(ssgsea_r[,-grep('DMSO',colnames(ssgsea_r))])
      if(ncol(ssgsea_r_case)>4){
        pvalue<-c()
        tf_n<-c()
        for(i in 1:nrow(ssgsea_r_ctl)){
          t_test<-t.test(ssgsea_r_ctl[i,],ssgsea_r_case[i,])
          pvalue<-c(pvalue,t_test$p.value)
          tf_n<-c(tf_n,rownames(ssgsea_r_ctl)[i])
        }
        t_test_mat<-cbind.data.frame(TF_name=tf_n,pvalue=pvalue)
      }
      drug_f<-strsplit(drug_f,'[.]')[[1]]
      drug_f<-paste0(drug_f[1],'.txt')
      print(drug_f)
      write.table(ssgsea_r,paste0(argv1[3],'ssgsea/','ssgsea_',drug_f),sep='\t',quote = F)
      drug_f1<-cbind.data.frame(ford_mat,drug_f1)
      write.table(drug_f1,paste0(argv1[3],'input_data/','DMSO_',drug_f),sep='\t',quote = F)
      write.table(t_test_mat,paste0(argv1[3],'ssgsea_ttest/',drug_f,'_','t_test_mat.txt'),sep='\t',quote=F,row.names=F)
      }
  }
  
}

