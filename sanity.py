import pandas as pd
d = pd.read_csv('g205_z4_host_halo.csv')
d.loc[d.Sub_SubhaloRankInGr != 0,
      ['snap','Redshift','Group_GroupID','Sub_SubhaloID','Sub_SubhaloRankInGr']]

gmc[gmc.Run == 'g205z4'].groupby('Snapshot').HaloID.unique()

d.plot(x='Redshift', y='phys_Group_Group_M_Crit200_Msun', logy=True)
