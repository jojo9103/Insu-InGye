import pandas as pd
import numpy as np
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.concat as cg
from cmapPy.pandasGEXpress.parse import parse
import cmapPy.pandasGEXpress.write_gctx as wg
import cmapPy.pandasGEXpress.subset_gctoo as sg
# cmapPy = LINCS data 다루기위한 module
# https://clue.io/cmapPy/pandasGEXpress.html 참고

#파일 경로
LINCS='/home/jojo9103/LINCS/LINCS/'
# LINCS_Cell = Cell line 정보
# LINCS_G = Gene들의 최신화 (duplication, withdrawn등)
# GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz = LINCS id 정보
LINCS_G=pd.read_csv(f'{LINCS}/LINCS_Gene.txt',sep='\t')
Cell_info=pd.read_csv(f'{LINCS}/LINCS_Cell.txt',sep='\t')
inst_info=pd.read_csv(f'{LINCS}/GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz',sep='\t',compression='gzip')

inst_info1=inst_info[inst_info['cell_id'].isin(Cell_info['cell_id'])]

# LINCS_data = gctx 핵심 데이터! 정리해야할 데이터 읽어오기
LINCS_data=parse(f'{LINCS}GSE70138_Broad_LINCS_Level3_INF_mlr12k_n345976x12328_2017-03-06.gctx')
LINCS_data=LINCS_data.data_df

# dup_g = duplicated 유전자 확인
dup_g=LINCS_G[LINCS_G['Symbol'].duplicated()]['Symbol']
# for문 이용해서 duplicated gene 제거 (mean으로 합침)
for g in dup_g:
    dup_g1=LINCS_G[LINCS_G['Symbol'].isin([g])]['pr_gene_id']
    dup_mat=list(map(str,dup_g1.to_list()))
    dup_m=LINCS_data.loc[dup_mat,].apply(np.mean,axis=1)
    selet_dup_g=dup_m.index[dup_m!=max(dup_m)]
    LINCS_data=LINCS_data.loc[~LINCS_data.index.isin(selet_dup_g)]

# WithDrawn 제거
WithDr=LINCS_G[LINCS_G['Entrez']=='Withdrawn']['pr_gene_id']
for g in WithDr:
     LINCS_data=LINCS_data.loc[~LINCS_data.index.isin([str(g)])]

LINCS_G.index=LINCS_G['pr_gene_id']
LINCS_G1=LINCS_G.loc[map(int,LINCS_data.index)]
LINCS_G1=LINCS_G1.loc[:,['pr_gene_id','Entrez','Symbol']]

# inst_id_new_inst_id_20210125 = 34만개정도 되는 id에서 중복이 안되게 새로운 id를 만듬 (추후 과정에서 필요)
new_inst_id=pd.read_csv(f'{LINCS}/inst_id_new_inst_id_20210125.txt',sep='\t')

# 농도나, Cell line, 적용시간이 중복되는 경우 평균값으로 합침
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
                    # DMSO = 농도가 없기때문에 column에 -666인 경우 제외하고 진행
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
    # 각 약물에 대해서 파일로 저장하기
    LINCS_data_mat.to_csv(f'./Data/{dg}_subset.txt.gz',compression='gzip',sep='\t',index=False)
    LINCS_data_mat=LINCS_data.loc[:,s_inst['inst_id']]
    LINCS_data_mat.columns=s_inst['id_2']

