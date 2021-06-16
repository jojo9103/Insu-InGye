# setwd('~/data1/2020Year/drug_repositioning/')
# argv1<-c('./Cancer_Exp/', 'BRCA', '01,06', './BRCA_all.txt')
# argv[1]= help, All, ACC,BRCA ..  argv1[3]= All, 01, 02
# argv1<-c('./Cancer_Exp/', 'BRCA', '01,06', './BRCA_all.txt', 'TCGA_mapping.txt')
library(stringr)
argv1<-commandArgs(trailingOnly = T)
dir_path=argv1[1] # TCGA dir 

#dir_path='./data1/2020Year/drug_repositioning/Cancer_Exp/'
print(sprintf('%s-%s',argv1[2],argv1[3]))

# TN_list = 샘플에서 Tumor인지, Normal, Metastasis인지 확인
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes 참
TN_list<-list()
for(i in 1:40){
  if(i %in%c(1,3,5,9)){
    i1=paste0('0',i)
  TN_list[[i]]='T'
  }
  if(i %in%c(2,4,40)){
    if(i <10){
    i1=paste0('0',i)
    }
    TN_list[[i]]='R'
  }
  if(i %in%c(6,7)){
    if(i <10){
      i1=paste0('0',i)
    }
    TN_list[[i]]='M'
  }
  if(i %in%c(10,11,12,14)){
    TN_list[[i]]='M'
  }
}

# 설명충 어떻게 진행하는지 설명
if(argv1[1]=='help'|argv1[1]=='h'){
  print('Rscript TCGA_Sample_filter.R TCGA_Sample_folder_path Cancer_Type Sample_Type Output')
  print('For Example:')
  print('TCGA_Sample_folder_path (Cancer_Exp)-- ACC')
  print('                                     | BLCA')
  print('                                     | BRCA')
  print('                                     | CESC')
  print('                                     | CHOL')
  print('                                     L COAD')
  print('I want BRCA Cancer type and Sample type 01,06 output file name is BRCA_all.txt.gz')
  print('Rscript TCGA_Sample_filter.R ./Cancer_Exp/ BRCA 01,06 ./BRCA_all.txt TCGA_mapping.txt # BRCA Sample, Sample type : 01,06 , output : BRCA_all.txt.gz')
  print('Another')
  print('I want All TCGA Cancer type and Sample type All(Tumor,Normal) output file name is all.txt.gz')
  print('Rscript TCGA_Sample_filter.R ./Cancer_Exp/ All All ./all.txt TCGA_mapping.txt # All TCGA, Sample type : all , output : all.txt.gz')
  
}else{

######### data load
# 전체인경우 All or all, ALL
  print('error')
  }
  all_data=data_mat
}
######### data load
## output all_data

if(argv1[3]=='All'|argv1[3]=='ALL'|argv1[3]=='all'){
  col_n<-colnames(all_data)
  col_n<-str_split_fixed(col_n,'-',4)[,4]
  TN<-unlist(apply(as.matrix(col_n),1,function(x)TN_list[[as.numeric(x)]]))
  col_n<-paste0(TN,'_',col_n)
  col_n[1]='Sample_type'
  all_data_test<-rbind.data.frame(Type=col_n,all_data)
}else{
  re_list<-list()
  ty=strsplit(argv1[3],',')[[1]]
  for(i in ty){
    TN<-TN_list[[as.numeric(i)]]
    col_n<-colnames(all_data)
    col_n<-str_split_fixed(col_n,'-',4)[,4]
    col_n<-paste0(TN,'_',col_n)
    col_n[1]='Sample_type'
    all_data1<-rbind.data.frame(Type=col_n,all_data)
    row_n<-all_data1[,1]
    if(length(grep(i,x = all_data1[1,]))!=0){
	c_name=colnames(all_data1)[grep(i,x=all_data1[1,])]
	all_data1<-all_data1[,grep(i,x = all_data1[1,])]
    	re_list[[i]]= as.data.frame(all_data1)
	colnames(re_list[[i]])=c_name
    }
  }
  all_data_test<-row_n
  for(i in re_list){
  all_data_test<-cbind.data.frame(all_data_test,i)
  }
}
# Gene = TCGA 에서 duplicated, withdrawn 
Gene<-read.table(argv1[5],sep='\t',header=T,stringsAsFactors = F)
Gene<-Gene[order(Gene$Entrez),]
colnames(Gene)[2]='old_symbol'
colnames(all_data_test)[1]='old_symbol'
all_data_test1<-merge(Gene,all_data_test,key='old_symbol',all=T)
print(dim(all_data_test1))

ti_type=all_data_test1[grep('Sample_type',all_data_test1[,1]),]
val_type=all_data_test1[-grep('Sample_type',all_data_test1[,1]),]
val_type=val_type[-grep('Withdrawn',val_type$Entrez),]
print(dim(val_type))
val_type<-val_type[order(as.numeric(val_type$Entrez)),]
all_data_test<-rbind.data.frame(ti_type,val_type)

all_data_test<-all_data_test[,c(1,4:ncol(all_data_test))]

dup_gene<-all_data_test$Symbol[duplicated(all_data_test$Symbol)]
dup_gene<-names(table(dup_gene))
dup_all_data_test<-all_data_test[all_data_test$Symbol%in%dup_gene,]
n_dup_all_data_test<-all_data_test[!all_data_test$Symbol%in%dup_gene,]
print(dim(n_dup_all_data_test))
for(g in dup_gene){
  mat<-all_data_test[all_data_test$Symbol%in%g,]
  mat1<-mat[,4:ncol(mat)]
  mat1<-apply(as.matrix(mat1),2,as.numeric)
  mat1<-apply(mat1,1,sum)
  n_dup_all_data_test<-rbind.data.frame(n_dup_all_data_test,mat[c(1:length(mat1))[mat1==max(mat1)][1],])
}
		   
file1<-trimws(argv1[4])
file2<-paste0(dir_path,argv1[2],'/',file1)
write.table(n_dup_all_data_test,file2,sep='\t',quote = F,row.names = F)
}













