import pandas as pd
import numpy as np
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.concat as cg
from cmapPy.pandasGEXpress.parse import parse
import cmapPy.pandasGEXpress.write_gctx as wg
import cmapPy.pandasGEXpress.subset_gctoo as sg
LINCS='/home/jojo9103/LINCS/LINCS/'
LINCS_G=pd.read_csv(f'{LINCS}/LINCS_Gene.txt',sep='\t')
Cell_info=pd.read_csv(f'{LINCS}/LINCS_Cell.txt',sep='\t')
inst_info=pd.read_csv(f'{LINCS}/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz',sep='\t',compression='gzip')

inst_info1=inst_info[inst_info['cell_id'].isin(Cell_info['cell_id'])]
#inst_info1
LINCS_data=parse(f'{LINCS}GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx')
LINCS_data=LINCS_data.data_df
#LINCS_data.data_df.apply(np.mean,axis=1)
dup_g=LINCS_G[LINCS_G['Symbol'].duplicated()]['Symbol']
for g in dup_g:
    dup_g1=LINCS_G[LINCS_G['Symbol'].isin([g])]['pr_gene_id']
    dup_mat=list(map(str,dup_g1.to_list()))
    dup_m=LINCS_data.loc[dup_mat,].apply(np.mean,axis=1)
    selet_dup_g=dup_m.index[dup_m!=max(dup_m)]
    LINCS_data=LINCS_data.loc[~LINCS_data.index.isin(selet_dup_g)]

# remove WithDrawn 
WithDr=LINCS_G[LINCS_G['Entrez']=='Withdrawn']['pr_gene_id']
for g in WithDr:
     LINCS_data=LINCS_data.loc[~LINCS_data.index.isin([str(g)])]

LINCS_G.index=LINCS_G['pr_gene_id']
LINCS_G1=LINCS_G.loc[map(int,LINCS_data.index)]
LINCS_G1=LINCS_G1.loc[:,['pr_gene_id','Entrez','Symbol']]

new_inst_id=pd.read_csv(f'{LINCS}/inst_id_new_inst_id_20210125.txt',sep='\t')

for dg in set(new_inst_id['pert_iname']):
    print(dg)
    s_inst=new_inst_id[new_inst_id['pert_iname']==dg]
    result_mat={}
    for dose in set(s_inst['pert_dose']):
        dose_inst=s_inst[s_inst['pert_dose']==dose]
        for cl in set(dose_inst['cell_id']):
                cl_inst=dose_inst[dose_inst['cell_id']==cl]
                for hr in set(cl_inst['pert_time']):
                    hr_inst1=cl_inst[cl_inst['pert_time']==hr]
                    key_id=hr_inst1['id_1'].values[0]
                    if dg=='DMSO':
                        key_id1=key_id.split('-666_')
                        key_id=''.join(key_id1)
                    Ldata_sub=LINCS_data.loc[:,hr_inst1['inst_id']]
                    Ldata_sub=Ldata_sub.apply(np.mean,1)
                    result_mat[hr_inst1['id_1'].values[0]]=Ldata_sub
    LINCS_data_mat=pd.DataFrame(result_mat)
    dg=dg.split(' ')
    dg='_'.join(dg)
    dg=dg.split('(')
    dg='_'.join(dg)
    dg=dg.split(')')
    dg='_'.join(dg)
    dg=dg.split('/')
    dg='_'.join(dg)
    dg=dg.split("'")
    dg='_'.join(dg)
    dg=dg.split('&')
    dg='_'.join(dg)
    LINCS_data_mat.index=Ldata_sub.index
    LINCS_data_mat=pd.concat([LINCS_G1.reset_index(drop=True),LINCS_data_mat.reset_index(drop=True)],axis=1)
    if dg=='DMSO':
        LINCS_data_mat.columns=[''.join(i.split('-666_')) for i in LINCS_data_mat.columns]
    LINCS_data_mat.to_csv(f'./Data/{dg}_subset.txt.gz',compression='gzip',sep='\t',index=False)
    LINCS_data_mat=LINCS_data.loc[:,s_inst['inst_id']]
    LINCS_data_mat.columns=s_inst['id_2']
#    if dg=='DMSO':
#        LINCS_data_mat.columns=[''.join(i.split('-666_')) for i in test.columns]
#    LINCS_data_mat.to_csv(f'./Data_1/{dg}_subset.txt.gz',compression='gzip',sep='\t')

