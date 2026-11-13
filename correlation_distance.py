import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
fname = "clumps/gmc_catalog/gmc_output.hdf5"
group_name = "GMC_alpha20_n100"
fname = "clumps/gmc_catalog/gmc_output.hdf5"
group_name = "GMC_alpha20_n100"
def load_catalog(filename, group):
    with h5py.File(filename, "r") as f:
        g = f[group]
        data = {k: g[k][()] for k in g.keys() if len(g[k].shape) <= 1}
    return pd.DataFrame(data)
df = load_catalog(fname, group_name)
if df is not None:
    metal_cols = ['C_Fraction', 'N_Fraction', 'O_Fraction', 'Ne_Fraction', 'Mg_Fraction', 'Si_Fraction', 'Fe_Fraction']
    df['Total_Z'] = df[metal_cols].sum(axis=1)
    z_sun = 0.0134
    df['log_Z_Zsun'] = np.log10(df['Total_Z'] / z_sun + 1e-10)
    df_filtered = df[df['Number_of_Particles'] > 30].copy()
    df_filtered['Run'] = df_filtered['Run'].str.decode('utf-8')
    host_metallicity = df_filtered.groupby('Run')['log_Z_Zsun'].transform('mean')
    df_filtered['Host_Disk_Z'] = host_metallicity
    print(df_filtered.columns)
    plt.figure(figsize=(8, 6))
    plt.scatter(df_filtered['EffectiveRadius'], df_filtered['log_Z_Zsun'], 
                alpha=0.5, edgecolors='k', c='orange')
    plt.xscale('log') 
    plt.xlabel('Effective Radius (kpc)')
    plt.ylabel(r'$\log(Z_{clump} / Z_{\odot})$')
    plt.title('Clump Metallicity vs. Effective Radius (z=4)')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.savefig('radius_metallicity_correlation.png', dpi=300)
    plt.show()
