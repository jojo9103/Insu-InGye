import pandas as pd
import numpy as np
LINCS='/home/jojo9103/LINCS/LINCS/'

LINCS_inst=pd.read_csv(f'{LINCS}/inst_id_new_inst_id_20210125.txt',sep='\t')

for dg in set(LINCS_inst['pert_iname']):
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
  mat=pd.read_csv(f'./Data_1/{dg}_subset.txt.gz',sep='\t',compression='gzip')
  L_inst=LINCS_inst[LINCS_inst['pert_iname']==ori_dg]
  m_dir={}
  for cl in set(L_inst['cell_id']):
    L_inst1=L_inst[L_inst['cell_id']==cl]
    for ds in set(L_inst['pert_dose']):
      ds_inst1=L_inst1[L_inst1['pert_dose']==ds]
      for hr in set(ds_inst1['pert_time']):
          hr_inst1=ds_inst1[ds_inst1['pert_time']==hr]
          key_id=hr_inst1['id_1'].values[0]
          if dg =='DMSO':
              key_id1=key_id.split('-666_')
              key_id=''.join(key_id1)
          mat1=mat.loc[:,hr_inst1['id_2']]
          m_dir[key_id]=mat1.apply(np.mean,1)
  result=pd.DataFrame(m_dir)
  result.index=mat['rid']
  print(dg)
  result.to_csv(f'./Data_2/{dg}_average_subset.txt.gz',sep='\t',compression='gzip')

