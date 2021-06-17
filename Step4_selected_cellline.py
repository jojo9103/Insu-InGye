import pandas as pd

LINCS='/home/jojo9103/LINCS/LINCS/'
# selected cell line info
s_LINCS=pd.read_csv(f'{LINCS}/LINCS_Cell.txt',sep='\t')

new_inst_id=pd.read_csv(f'{LINCS}/inst_id_new_inst_id_20210125.txt',sep='\t')

for dg in set(new_inst_id['pert_iname']):
        ori_dg=dg
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
        print(dg)
        mat=pd.read_csv(f'{LINCS}/Data/{dg}_subset.txt.gz',sep='\t',compression='gzip',dtype='str')
        mat['Entrez']=pd.to_numeric(mat['Entrez'])
        mat=mat.sort_values(by=['Entrez'],axis=0)
        selected=['pr_gene_id','Entrez','Symbol']
        # 원하는 Cell line들만 가져오기
        selected1=[i for i in mat.columns if i.split('_')[0] in s_LINCS['cell_id'].values]
        selected_1=[i.split('_') for i in selected1]
        # DMSO가 아닌경우 농도와 시간이 있는데, 농도 기준 정렬후에 시간기준 정렬함.
        if dg!='DMSO':
            selected1=sorted(selected_1,key=lambda x : (float(x[2]),float(x[1])))
        # DMSO가 아닌경우 농도와 시간이 있는데, 시간기준 정렬함.
        if dg=='DMSO':
            selected1=sorted(selected_1,key=lambda x : float(x[1]))
        selected1=['_'.join(i) for i in selected1]
        selected+=selected1
        mat1=mat[selected]
        if dg=='DMSO':
            mat1.columns=[''.join(i.split('_-666')) for i in mat1.columns]
        if mat1.shape[1]>3:
            mat1.to_csv(f'{LINCS}/Data1/{dg}_Selected_CL.txt.gz',sep='\t',compression='gzip' ,index=False) #index= False

