# BRCA = 각 TF에 대해서 Tumor mean, Normal mean, t-test p-value를 갖는 파일 
BRCA<-read.table('~/data1/2021Year/LINCS/20210415/BRCA_Summary.txt',sep='\t',header=T,stringsAsFactors = F)
LINCS<-read.table('~/data1/2021Year/LINCS/20210415/Combine_anal_ssgsea_diff.txt',sep='\t',header=T,stringsAsFactors = F,check.names = F)
BRCA<-BRCA[BRCA$TF_name%in%LINCS$TF,]
print(table(LINCS$TF==BRCA$TF_name))
BRCA$diff=BRCA$Tumor_mean-BRCA$Normal_mean
BRCA$adj=p.adjust(BRCA$pvalue,method='fdr')
# TCGA ssgsea score
BRCA_raw<-read.table('~/data1/2021Year/LINCS/20210415/plot_summary/ssgsea_symbol.txt',sep='\t',header=T,stringsAsFactors = F,check.names = F)

library(pheatmap)
library(RColorBrewer)
case<-BRCA_raw[,grep('-01',colnames(BRCA_raw))]
control<-BRCA_raw[,grep('-11',colnames(BRCA_raw))]

hc_case<-hclust(dist(t(case)))
hc_control<-hclust(dist(t(control)))

plot(hc_case,hang=-1)

hc_case_cut<-cutree(hc_case,4)
hc_case_cut<-hc_case_cut[order(hc_case_cut)]
case<-case[,names(hc_case_cut)]

hc_control_cut<-cutree(hc_control,4)
hc_control_cut<-hc_control_cut[order(hc_control_cut)]
control<-control[,names(hc_control_cut)]
bar_info<-data.frame(sample=c(rep('control',ncol(control)),rep('case',ncol(case))))
rownames(bar_info)<-c(names(control),names(case))
# 800*600
pheatmap(cbind.data.frame(control,case),cluster_cols = F,annotation_col = bar_info,scale = 'row',fontsize_row = 3,fontsize_col = 3,#show_colnames = F,show_rownames = F,treeheight_row = 0,
         color=colorRampPalette(rev(brewer.pal(n=7,name ='RdYlBu')))(100),breaks = seq(from=-2,to=2,length.out = 100) )


combine_mat<-cbind.data.frame(control,case)
Summary_mat<-read.table('~/data1/2021Year/LINCS/20210415/BRCA_Summary.txt',sep='\t',header=T,stringsAsFactors = F)
Summary_mat$adj<-p.adjust(Summary_mat$pvalue,method = 'fdr')
Summary_mat<-Summary_mat[Summary_mat$adj<0.05,]
combine_mat<-combine_mat[rownames(combine_mat)%in%Summary_mat$TF_name,]
# 800*600
pheatmap(combine_mat,cluster_cols = F,annotation_col = bar_info,scale = 'row',fontsize_row = 3,fontsize_col = 3,#show_colnames = F,show_rownames = F,treeheight_row = 0,
         color=colorRampPalette(rev(brewer.pal(n=7,name ='RdYlBu')))(100),breaks = seq(from=-2,to=2,length.out = 100))

# Summary_mat<-Summary_mat[Summary_mat$adj<0.0005,]
# combine_mat<-combine_mat[rownames(combine_mat)%in%Summary_mat$TF_name,]
# pheatmap(combine_mat,cluster_cols = F,annotation_col = bar_info,scale = 'row',fontsize_row = 3,fontsize_col = 3,
#          color=colorRampPalette(rev(brewer.pal(n=7,name ='RdYlBu')))(100),breaks = seq(from=-2,to=2,length.out = 100))

TF<-c()
p<-c()
for(i in 1:nrow(case)){
  ttest<-t.test(case[i,],control[i,])
  TF<-c(TF,rownames(case)[i])
  p<-c(p,ttest$p.value)
}
sub_result<-cbind.data.frame(TF,p)
sub_result<-sub_result[sub_result$TF%in%LINCS$TF,]



# Tumor<-read.table('~/data1/2021Year/LINCS//20210415/plot_summary/TCGA_RNA-Seq_BRCA_T.txt',sep='\t',header=T,stringsAsFactors = F,quote = '',check.names = F)
# Tumor<-Tumor[-1,]
# rownames(Tumor)=Tumor[,3]
# Tumor<-Tumor[,4:ncol(Tumor)]
# normal<-read.table('~/data1/2021Year/LINCS//20210415/plot_summary/TCGA_RNA-Seq_BRCA_N.txt',sep='\t',header=T,stringsAsFactors = F,quote = '',check.names = F)
# normal<-normal[-1,]
# rownames(normal)=normal[,3]
# normal<-normal[,4:ncol(normal)]
# 
# Tumor1<-apply(Tumor,2,as.numeric)
# rownames(Tumor1)=rownames(Tumor)
# normal1<-apply(normal,2,as.numeric)
# rownames(normal1)=rownames(normal)
# table(rownames(Tumor1)==rownames(normal1))
# com<-cbind.data.frame(normal1,Tumor1)
# 
# com1_mean<-apply(com,1,mean)
# com<-com[com1_mean!=0,]

#
com<-read.table('/home/jojo9103/data1/2021Year/LINCS/20210428/ssgsea_symbol.txt',sep='\t',header=T,stringsAsFactors = F,check.names = F)


Tumor1<-Tumor1[rownames(Tumor1)%in%rownames(com),]
normal1<-normal1[rownames(normal1)%in%rownames(com),]

exp_hc_case<-hclust(dist(t(Tumor1)))
exp_hc_control<-hclust(dist(t(normal1)))


exp_hc_case_cut<-cutree(exp_hc_case,4)
exp_hc_case_cut<-exp_hc_case_cut[order(exp_hc_case_cut)]
Tumor<-Tumor[,names(exp_hc_case_cut)]

exp_hc_control_cut<-cutree(exp_hc_control,4)
exp_hc_control_cut<-exp_hc_control_cut[order(exp_hc_control_cut)]
normal<-normal[,names(exp_hc_control_cut)]


bar_info<-data.frame(sample=c(rep('Normal',ncol(normal1)),rep('Tumor',ncol(Tumor1))))
rownames(bar_info)<-c(names(normal1),names(Tumor1))
# 800*600
# pheatmap(com,cluster_cols = F,annotation_col = bar_info,scale = 'row',fontsize_row = 3,fontsize_col = 3)



