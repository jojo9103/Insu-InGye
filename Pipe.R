library(stringr)

# 이 코드의 목적
# L-score를 계산하기 위함
# 그러기 위해서는 LINCS에서 각 TF별로 dose,time기준으로 같은 경향성을 갖고 있어야함.


LINCS<-read.table('~/data1/2021Year/LINCS/20210415/Combine_anal_ssgsea_diff.txt',sep='\t',header=T,stringsAsFactors = F,check.names = F)
LINCS_linear<-read.table('~/data1/2021Year/LINCS/20210415/lm_result_Summary.txt',sep='\t',header=T,stringsAsFactors = F)
LINCS_linear1<-LINCS_linear[LINCS_linear$lm_p<0.05,]
# TCGA 결과파일
BRCA<-read.table('~/data1/2021Year/LINCS/20210415/BRCA_Summary.txt',sep='\t',header=T,stringsAsFactors = F)
BRCA<-BRCA[BRCA$TF_name%in%LINCS$TF,]
print(table(LINCS$TF==BRCA$TF_name))
BRCA$diff=BRCA$Tumor_mean-BRCA$Normal_mean
# fdr이용해서 p-value 보정해줌.
# BRCA$adj=p.adjust(BRCA$pvalue,method='fdr') (ttest_Tumor_normal.R에서 fdr로 변경함)
BRCA$logp=-log10(BRCA$adj)
# Up regul in BRCA +, Up regul in Normal -
BRCA[BRCA$diff<0,"logp"]=BRCA[BRCA$diff<0,"logp"]*-1

# combine LINCS data, TCGA BRCA data 

#LINCS_data =  TCGA, LINCS데이터를 합친 파일 (Scatter plot만들때 필요함)
LINCS_data<-NULL
for(i in 2:ncol(LINCS)){
  i1=cbind.data.frame(TF=BRCA$TF_name,BRCA_TF=BRCA$logp,LINCS_TF=LINCS[,i],Drug_name=colnames(LINCS)[i])
  LINCS_data<-rbind.data.frame(LINCS_data,i1)
}
# write.table(LINCS_data,'~/data1/2021Year/LINCS/20210415/TCGA_LINCS_TF_score.txt',sep='\t',row.names = F,quote = F)


## LINCS filter ##
# LINCS_data1<-LINCS_data
LINCS_linear_1<-LINCS_linear
LINCS_linear_1$pos=''
# slope이 양수 음수 혹은  0인 경우를 확인함.
LINCS_linear_1$pos[LINCS_linear_1$lm>0]='pos'
LINCS_linear_1$pos[LINCS_linear_1$lm<0]='neg'
LINCS_linear_1$pos[LINCS_linear_1$lm==0]='zero'

drug_name<-names(table(LINCS_linear_1$drug_name))
# 약물마다 dose기준, time기준 linear regression slope에 대해 같은 경향성을 확인함.
for(dn in drug_name){
  # print(dn)
  drug_data=LINCS_linear_1[LINCS_linear_1$drug_name==dn,]
  drug_status=str_split_fixed(drug_data$dose,'',2)
  if(length(table(drug_status[,1]))>1){print(dn)}
  time_status<-drug_data[drug_status[,1]=='t',] # time
  dose_status<-drug_data[drug_status[,1]=='d',] # dosez
  # time
  drug_time<-drug_status[drug_status[,1]=='t',]
  if(!is.character(drug_time)){
    time_status<-time_status[as.numeric(drug_time[,2])==max(as.numeric(drug_time[,2])),]
    combine_status<-cbind.data.frame(TF=time_status$TF,pos_t=time_status$pos)
  }else{combine_status<-cbind.data.frame(TF=time_status$TF,pos_t=time_status$pos)}
  
  #dose
  if(nrow(dose_status)!=0){
    drug_dose<-drug_status[drug_status[,1]=='d',]
    dose_status<-dose_status[as.numeric(drug_dose[,2])==max(as.numeric(drug_dose[,2])),]
    combine_status<-merge(combine_status,cbind.data.frame(TF=dose_status$TF,pos_d=dose_status$pos),all=T)
    combine_status<-na.omit(combine_status)
    if(nrow(combine_status)>0){
      print('dd')
    combine_status$diff=''
    combine_status[combine_status$pos_t==combine_status$pos_d,'diff']=combine_status$pos_d[combine_status$pos_t==combine_status$pos_d]
    combine_status[combine_status$pos_t!=combine_status$pos_d,'diff']='zero'
    combine_status<-combine_status[,c(1,4)]
    }
  }
  if(nrow(combine_status)>0){
  colnames(combine_status)<-c('TF','diff')
  # LINCS filter
  LINCS_sub=LINCS[LINCS$TF%in%combine_status$TF,c('TF',dn)]
  LINCS_sub$dr_pos=''
  LINCS_sub$dr_pos[LINCS_sub[,2]>0]='pos'
  LINCS_sub$dr_pos[LINCS_sub[,2]<0]='neg'
  LINCS_sub$dr_pos[LINCS_sub[,2]==0]='zero'
  Filter_result<-merge(LINCS_sub,combine_status)
  Filter_result$filter=''
  Filter_result[Filter_result$dr_pos==Filter_result$diff,"filter"]='not'
  Filter_result[Filter_result$dr_pos!=Filter_result$diff,"filter"]='Filter'
  Filter_result_Gene<-Filter_result$TF[Filter_result$filter=='Filter']
  if(length(Filter_result_Gene)>0){
    LINCS[LINCS$TF%in%Filter_result_Gene,dn]=0
  }
  }
}


LINCS_score<-LINCS[,2:ncol(LINCS)]
dim(LINCS_score)
LINCS_score<-apply(LINCS_score,2,function(x)x*BRCA$logp)
dim(LINCS_score)
LINCS_score<-apply(LINCS_score,2,sum)
LINCS_score<-LINCS_score/519
LINCS_score<-cbind.data.frame(Drug_n=names(LINCS_score),L_score=LINCS_score)
LINCS_score<-LINCS_score[order(LINCS_score$L_score,decreasing = F),]

# 800 * 500 L_score_hist
hist(LINCS_score$L_score,main='L-score distribution',xlab='L-score')

table(LINCS$TF==BRCA$TF_name)
LINCS_data1<-NULL
for(i in 2:ncol(LINCS)){
  i1=cbind.data.frame(TF=BRCA$TF_name,BRCA_TF=BRCA$logp,LINCS_TF=LINCS[,i],Drug_name=colnames(LINCS)[i])
  LINCS_data1<-rbind.data.frame(LINCS_data1,i1)
}

# write.table(LINCS_data1,'~/data1/2021Year/LINCS/20210415/TCGA_LINCS_TF_score_Lscore.txt',sep='\t',row.names = F,quote = F)

# 800*500 L_score (linear filter)
library(ggplot2)
ggplot(data=LINCS_data1)+geom_point(aes(x=LINCS_TF,y=BRCA_TF),size=.5)+theme_bw()+xlab('LINCS enrichment score')+ylab('TCGA enrichment score')+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=10))+geom_vline(xintercept = 0,linetype='dashed',col='brown')+geom_hline(yintercept = 0,linetype='dashed',col='brown')

# 800*500 L_score1 (raw points)
ggplot(data=LINCS_data)+geom_point(aes(x=LINCS_TF,y=BRCA_TF),size=.5)+theme_bw()+xlab('LINCS enrichment score')+ylab('TCGA enrichment score')+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=10))+geom_vline(xintercept = 0,linetype='dashed',col='brown')+geom_hline(yintercept = 0,linetype='dashed',col='brown')

LINCS_score2<-LINCS_score[order(LINCS_score$L_score,decreasing = F),]
# top candidate
for(i in 1:10){
  LINCS_data2=LINCS_data1
  LINCS_data2$select=''
  print(LINCS_score2$Drug_n[i])
  LINCS_data2$select[LINCS_data2$Drug_name==LINCS_score2$Drug_n[i]]='Candidate drug'
  p<-ggplot(data=LINCS_data2)+geom_point(aes(x=LINCS_TF,y=BRCA_TF,color=select),size=0.5)+
    geom_point(data=subset(LINCS_data1,Drug_name==LINCS_score2$Drug_n[i]),aes(x=LINCS_TF,y=BRCA_TF),col="#E69F00",size=0.5)+
    theme_bw()+xlab('LINCS score')+ylab('TCGA score')+
    scale_color_manual(values=c("#999999","#E69F00"))+ggtitle(LINCS_score2$Drug_n[i])+
    theme(axis.title = element_text(size=15),axis.text = element_text(size=10))+
    geom_vline(xintercept = 0,linetype='dashed',col='brown')+
    geom_hline(yintercept = 0,linetype='dashed',col='brown')+
    geom_smooth(data=subset(LINCS_data1,Drug_name==LINCS_score2$Drug_n[i]),aes(x=LINCS_TF,y=BRCA_TF),col="#E69F00",size=1.5,method='auto',se=F)
  ggsave(filename = paste0('/home/jojo9103/data1/2021Year/LINCS/20210415/Candidate_drug/top/',LINCS_score2$Drug_n[i],'.tiff'),plot = p,width = 200,height = 100,units = 'mm')
}

# bottom
LINCS_score2<-LINCS_score[order(LINCS_score$L_score,decreasing = T),]
for(i in 1:10){
  LINCS_data2=LINCS_data1
  LINCS_data2$select=''
  print(LINCS_score2$Drug_n[i])
  LINCS_data2$select[LINCS_data2$Drug_name==LINCS_score2$Drug_n[i]]='Candidate drug'
  p<-ggplot(data=LINCS_data2)+geom_point(aes(x=LINCS_TF,y=BRCA_TF,color=select),size=0.5)+
    geom_point(data=subset(LINCS_data1,Drug_name==LINCS_score2$Drug_n[i]),aes(x=LINCS_TF,y=BRCA_TF),col="#E69F00",size=0.5)+
    theme_bw()+xlab('LINCS score')+ylab('TCGA score')+
    scale_color_manual(values=c("#999999","#E69F00"))+ggtitle(LINCS_score2$Drug_n[i])+
    theme(axis.title = element_text(size=15),axis.text = element_text(size=10))+
    geom_vline(xintercept = 0,linetype='dashed',col='brown')+
    geom_hline(yintercept = 0,linetype='dashed',col='brown')+
    geom_smooth(data=subset(LINCS_data1,Drug_name==LINCS_score2$Drug_n[i]),aes(x=LINCS_TF,y=BRCA_TF),col="#E69F00",size=1.5,method='auto',se=F)
  ggsave(filename = paste0('/home/jojo9103/data1/2021Year/LINCS/20210415/Candidate_drug/bottom/',LINCS_score2$Drug_n[i],'.tiff'),plot = p,width = 200,height = 100,units = 'mm')
}



## Combination score

Drug_name<-c()
Drug_lscore<-c()
# Combine_score_mat<-LINCS$TF
for(i in 2:ncol(LINCS)){
  print(sprintf('%s/%s',i,ncol(LINCS)))
  for(l in 2:ncol(LINCS)){
    if(i!=l){
      combine_mat<-cbind.data.frame(LINCS[,c(i,l)])
      combine_mat1<-apply(combine_mat,2,abs)
      combine_mat$TF_score=0
      combine_mat$TF_score[combine_mat1[,1]>combine_mat1[,2]]=combine_mat[combine_mat1[,1]>combine_mat1[,2],1]
      combine_mat$TF_score[combine_mat1[,1]<combine_mat1[,2]]=combine_mat[combine_mat1[,1]<combine_mat1[,2],2]
      combine_mat$TF_score[combine_mat1[,1]==combine_mat1[,2]]=combine_mat[combine_mat1[,1]==combine_mat1[,2],2]
      Drug_name<-c(Drug_name,paste0(colnames(LINCS)[i],'_',colnames(LINCS)[l]))
      # Combine_score_mat<-cbind.data.frame(Combine_score_mat,combine_mat)
      Drug_lscore<-c(Drug_lscore,sum(combine_mat$TF_score))
    }else{
      # Combine_score_mat<-cbind.data.frame(Combine_score_mat,LINCS[,i])
      combine_mat<-cbind.data.frame(LINCS[,c(i,l)])
      Drug_lscore<-c(Drug_lscore,sum(combine_mat[,1]))
      Drug_name<-c(Drug_name,colnames(LINCS)[i])
    }
  }
}


# write.table(LINCS,'~/data1/2021Year/LINCS/20210415/Combine_score/LINCS.txt',sep='\t',row.names = F,quote = F)















