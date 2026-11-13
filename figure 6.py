import h5py
import pandas as pd
import numpy as np
fname = "clumps/gmc_catalog/gmc_output.hdf5"
group_name = "GMC_alpha20_n100"
def load_catalog(filename, group):
    try:
        with h5py.File(filename, "r") as f:
            if group not in f:
                print(f"Error: Group '{group}' not found.")
                return None
            g = f[group]
            data = {} 
            for key in g.keys():
                if len(g[key].shape) <= 1:
                    data[key] = g[key][()]
                else:
                    continue 
            return pd.DataFrame(data)
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        return None
df = load_catalog(fname, group_name)
if df is not None:
    metal_columns = [
        'C_Fraction', 'N_Fraction', 'O_Fraction', 
        'Ne_Fraction', 'Mg_Fraction', 'Si_Fraction', 'Fe_Fraction'
    ]
    df['Total_Z'] = df[metal_columns].sum(axis=1)
    z_sun = 0.0134
    df['log_Z_Zsun'] = np.log10(df['Total_Z'] / z_sun + 1e-10)
    df_filtered = df[df['Number_of_Particles'] > 30].copy()
    print(f"Successfully loaded {len(df_filtered)} resolved clumps.")
    print("\n--- Figure 6 Data Preview ---")
    cols_to_show = ['Run', 'StellarMass', 'SFR_inst', 'Total_Z', 'log_Z_Zsun']
    print(df_filtered[cols_to_show].head())

else:
    print("Failed to load data. Please check your folder structure.")
df_filtered['Run'] = df_filtered['Run'].str.decode('utf-8')
pd.options.display.float_format = '{:.4f}'.format
print(df_filtered[cols_to_show].head())
