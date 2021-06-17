library(GSVA)
library(stringr)
library(openxlsx)

argv1<-commandArgs(trailingOnly = T)

geneset<-read.xlsx(argv1[1])
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


o_drug<-c() # outlier drug

TF_diff=NULL
drug_name=c()

for(drug_f in filename1){
  i_n=strsplit(drug_f,'_')[[1]][1]
  i_n1=strsplit(drug_f,'.gz')[[1]][1]
  i_n2=strsplit(drug_f,'_Selected_CL.txt.gz')[[1]][1]
  print(i_n2)
  if(i_n!='DMSO'){
    drug_mat=read.table(gzfile(paste0(argv1[2],drug_f)),sep='\t',header = T,stringsAsFactors = F,quote = '')
    if(ncol(drug_mat)>3){
      drug_mat1<-as.data.frame(drug_mat[,4:ncol(drug_mat)])
      rownames(drug_mat1)=drug_mat[,3]
      colnames(drug_mat1)=colnames(drug_mat)[4:ncol(drug_mat)]
      drug_f1_coln<-str_split_fixed(colnames(drug_mat1),'_',4)
      if(i_n2=='BRD-K68548958'){
	      drug_f1_coln[drug_f1_coln[,2]=='6'&drug_f1_coln[,3]=='25',2]='0'
	      drug_f1_coln[drug_f1_coln[,2]=='0'&drug_f1_coln[,3]=='25',3]='0'
      }
      if(i_n2=='pazopanib'){
	      drug_f1_coln[drug_f1_coln[,2]=='6'&drug_f1_coln[,3]=='40',2]='0'
	      drug_f1_coln[drug_f1_coln[,2]=='0'&drug_f1_coln[,3]=='40',3]='0'
      }
      if(i_n2=='pyrazolanthrone'){
	      drug_f1_coln[drug_f1_coln[,2]=='6'&drug_f1_coln[,3]=='25',2]='0'
	      drug_f1_coln[drug_f1_coln[,2]=='0'&drug_f1_coln[,3]=='25',3]='0'	  
      }


      dose_max<-max(as.numeric(drug_f1_coln[,3]))
      time_max<-max(as.numeric(drug_f1_coln[,2]))
      drug_mat2<-as.data.frame(drug_mat1[,drug_f1_coln[,2]==time_max&drug_f1_coln[,3]==dose_max])
      rownames(drug_mat2)=drug_mat[,3]
      colnames(drug_mat2)=colnames(drug_mat1)[drug_f1_coln[,2]==time_max&drug_f1_coln[,3]==dose_max]
      DMSO_sub<-DMSO[,grep(paste0('_',time_max,'_'),colnames(DMSO))]
      if(ncol(drug_mat)==0){o_drug<-c(o_drug,drug_f)
      }else{
        input_data=cbind.data.frame(DMSO_sub,drug_mat2)
        colnames(input_data)<-c(colnames(DMSO_sub),colnames(drug_mat2))
        input_data<-as.matrix(input_data)
        ssgsea_r<-gsva(input_data,geneset_list,method='ssgsea',kcdf="Gaussian",abs.ranking=FALSE,
                       min.sz=2,max.sz=Inf,parallel.sz=10,mx.diff=TRUE,
                       ssgsea.norm=TRUE,
                       verbose=F) # verbose =TRUE
        ssgsea_r_zscore<-t(apply(ssgsea_r,1,scale))
        colnames(ssgsea_r_zscore)=colnames(ssgsea_r)
        rownames(ssgsea_r_zscore)=rownames(ssgsea_r)
        write.table(ssgsea_r,paste0(argv1[3],'/ssgsea/','ssgsea1_',i_n1,'.txt'),sep='\t',quote = F)
        write.table(ssgsea_r_zscore,paste0(argv1[3],'/ssgsea_zscore/','ssgsea_zscore_',i_n1,'.txt'),sep='\t',quote = F)
        # ssgsea_control
        ssgsea_r_zscore_ctl<-ssgsea_r_zscore[,colnames(ssgsea_r_zscore)%in%colnames(DMSO_sub)]
        ssgsea_r_zscore_case<-as.data.frame(ssgsea_r_zscore[,colnames(ssgsea_r_zscore)%in%colnames(drug_mat2)])
	colnames(ssgsea_r_zscore_case)<-colnames(drug_mat2)
        TF<-c()
        TF_mean<-c()
        for(i in 1:nrow(ssgsea_r_zscore_ctl)){
          mean_diff=mean(as.numeric(ssgsea_r_zscore_case[i,]))-mean(as.numeric(ssgsea_r_zscore_ctl[i,]))
          TF=append(TF,rownames(ssgsea_r_zscore_ctl)[i])
          TF_mean=append(TF_mean,mean_diff)
        }
	if(!is.null(TF_diff)){TF_diff=cbind.data.frame(TF_diff,TF_mean)}
        if(is.null(TF_diff)){TF_diff=cbind.data.frame(TF,TF_mean)}
      }
      drug_name<-append(drug_name,i_n2)
    }
  }
  
}
colnames(TF_diff)=c('TF',drug_name)
write.table(TF_diff,paste0(argv1[3],'Combine_anal_ssgsea_diff.txt'),sep='\t',quote = F,row.names=F)
write.table(o_drug,paste0(argv1[3],'outlier.txt'),sep='\t',quote = F)


