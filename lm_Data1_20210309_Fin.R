library(stringr)

argv1<-commandArgs(trailingOnly = T)

for(i in argv1){
  i1=strsplit(i,'/')[[1]]
  i1=i1[length(i1)]
  i1=strsplit(i1,'[.]')[[1]]
  i1=i1[1]
  
 print(i1)
  ssgsea<-read.table(i,sep='\t',header=T,stringsAsFactors = F,quote="",fill=T)
  ssgsea_1<-as.matrix(ssgsea[,-grep('DMSO',colnames(ssgsea))])
  
  if(ncol(ssgsea_1)>1){
    coln_mat<-str_split_fixed(colnames(ssgsea_1),'_',4)
#   if(i1=='ssgsea_BRD-K68548958_Selected_CL'){
#	   coln_mat=coln_mat[coln_mat[,2]!='6'&coln_mat[,3]!='25',]
#   }
#    if(i1=='ssgsea_pazopanib_Selected_CL'){
#	    coln_mat=coln_mat[coln_mat[,2]!='6'&coln_mat[,3]!='40',]
#    }
#    if(i1=='ssgsea_pyrazolanthrone_Selected_CL'){
#	    coln_mat=coln_mat[coln_mat[,2]!='6'&coln_mat[,3]!='25',]    
#    }
    time1<-as.numeric(coln_mat[,2])
    dose1<-as.numeric(coln_mat[,3])
    ti_do<-cbind.data.frame(time1,dose1)
    dose_result<-data.frame(TF=rownames(ssgsea_1))
    time_result<-data.frame(TF=rownames(ssgsea_1))
    # dose1
    for(dos in names(table(ti_do[,'dose1']))){
      dos_mat=ssgsea_1[,ti_do[,'dose1']==dos]
      dos_time=time1[dose1==dos]
      if(dim(as.data.frame(dos_mat))[2]>1){
      time_m<-apply(dos_mat,1,function(x)lm(x~dos_time)$coefficients['dos_time'])
      if (length(table(dos_time))!=1){
        time_mp<-apply(dos_mat,1,function(x)summary(lm(x~dos_time))$coefficients[2,'Pr(>|t|)'])
      }else{
        time_mp<-rep(NA,nrow(dos_mat))
      }
      result<-cbind.data.frame(TF=names(time_m),time_m,time_mp)
      colnames(result)<-c('TF',paste0('lm_d',dos),paste0('lm_p_d',dos))
      time_result<-merge(time_result,result,all=T,by = 'TF')
      # print(dim(time_result))
    }
    }
    rownames(time_result)=rownames(result)
    # time1
    for(ti in names(table(ti_do[,'time1']))){
      ti_mat=ssgsea_1[,ti_do[,'time1']==ti]
      ti_dose=dose1[time1==ti]
      if(dim(as.data.frame(ti_mat))[2]>1){
      dose_m<-apply(ti_mat,1,function(x)lm(x~ti_dose)$coefficients['ti_dose'])
      if(length(table(ti_dose))!=1){
        dose_mp<-apply(ti_mat,1,function(x)summary(lm(x~ti_dose))$coefficients[2,'Pr(>|t|)'])
      }else{
        dose_mp<-rep(NA,nrow(ti_mat))
      }
      result<-cbind.data.frame(TF=names(time_m),dose_m,dose_mp)
      colnames(result)<-c('TF',paste0('lm_t',ti),paste0('lm_p_t',ti))
      dose_result<-merge(dose_result,result,all=T,by = 'TF')
    }
    }
    rownames(dose_result)=rownames(result)

    Fin_result<-merge(time_result,dose_result)
    Fin_result<-Fin_result[,!is.na(Fin_result[1,])]
    if(dim(as.data.frame(Fin_result))[2]>1){
    write.table(Fin_result,paste0('./lm_result/',i1,'.txt'),sep='\t',row.names = F,quote=F)
    }
  }
}

