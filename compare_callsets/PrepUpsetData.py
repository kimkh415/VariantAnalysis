import joblib
import pandas as pd
import numpy as np


def resolve_type(row, type_col, id_col):
    """
    Extracts SV type from the column if available; otherwise parses the ID string.
    """
    if pd.notna(row.get(type_col)): return row[type_col]
    val = str(row.get(id_col, ''))
    if 'DEL' in val: return 'DEL'
    if 'INS' in val: return 'INS'
    if 'DUP' in val: return 'DUP'
    if 'INV' in val: return 'INV'
    return np.nan


def process_base_run(jl_file, base_name, comp_name):
    """
    Anchors on the base
    - Match 'tpbase' to 'tp' via MatchId.
    """
    data = joblib.load(jl_file)
    # Base (Matches + Misses) and Comp (Matches only, for ID lookup)
    base_df = data[data['state'].isin(['tpbase', 'fn'])].copy()
    comp_df = data[data['state'] == 'tp'][['MatchId', 'id', 'svtype']].copy()
    # Create Link Keys (MatchId tuple)
    base_df['join_key'] = base_df['MatchId'].apply(lambda x: tuple(x) if x is not None else None)
    comp_df['join_key'] = comp_df['MatchId'].apply(lambda x: tuple(x) if x is not None else None)
    # Merge Comp info onto Base
    merged = pd.merge(base_df, comp_df, on='join_key', how='left', suffixes=('_base', '_comp'))
    # Rename Columns dynamically
    merged = merged.rename(columns={
        'id_base': f'ID_{base_name}',
        'svtype_base': f'Type_{base_name}',
        'id_comp': f'ID_{comp_name}',
        'svtype_comp': f'Type_{comp_name}'
    })
    # Resolve Types
    merged[f'Type_{base_name}'] = merged.apply(lambda r: resolve_type(r, f'Type_{base_name}', f'ID_{base_name}'), axis=1)
    merged[f'Type_{comp_name}'] = merged.apply(lambda r: resolve_type(r, f'Type_{comp_name}', f'ID_{comp_name}'), axis=1)
    # Prioritize Match (tpbase=0) over FN (fn=1)
    merged['sort_rank'] = merged['state'].map({'tpbase': 0, 'fn': 1})
    df_unique = merged.sort_values([f'ID_{base_name}', 'sort_rank']).drop_duplicates(subset=[f'ID_{base_name}'], keep='first')
    return df_unique[[f'ID_{base_name}', f'Type_{base_name}', f'ID_{comp_name}', f'Type_{comp_name}']]


def process_residuals(jl_file, base_name, comp_name):
    """
    Handles variants missed by C.
    """
    data = joblib.load(jl_file)
    # S&P + S-only
    base_centric = process_base_run(jl_file, base_name, comp_name)
    # P-only
    comp_only = data[data['state'] == 'fp'].copy()
    comp_only = comp_only.rename(columns={'id': f'ID_{comp_name}', 'svtype': f'Type_{comp_name}'})
    comp_only[f'ID_{base_name}'] = np.nan
    comp_only[f'Type_{base_name}'] = np.nan
    comp_only[f'Type_{comp_name}'] = comp_only.apply(lambda r: resolve_type(r, f'Type_{comp_name}', f'ID_{comp_name}'), axis=1)
    comp_only = comp_only.drop_duplicates(subset=[f'ID_{comp_name}'])
    return pd.concat([base_centric, comp_only], ignore_index=True)


# load Truvari bench dir
# C-P
df_cp = process_base_run("result_cp.jl", "C", "P")
# C-S
df_cs = process_base_run("result_cs.jl", "C", "S")
master_c = pd.merge(df_cp, df_cs[[f'ID_C', f'ID_S', f'Type_S']], on='ID_C', how='outer')

# S-P
df_sp = process_residuals("result_sp.jl", "S", "P")
df_sp['ID_C'] = np.nan
df_sp['Type_C'] = np.nan

final_cols = ['ID_C', 'Type_C', 'ID_P', 'Type_P', 'ID_S', 'Type_S']
final_df = pd.concat([master_c[final_cols], df_sp[final_cols]], ignore_index=True)

summary = final_df.groupby(['Type_C', 'Type_P', 'Type_S'], dropna=False).size().reset_index(name='Count')
print(summary.to_string(index=False))

output_file = "data_ids.csv"
final_df.to_csv(output_file, index=False)

summary.to_csv("data_upset.csv", index=False)

