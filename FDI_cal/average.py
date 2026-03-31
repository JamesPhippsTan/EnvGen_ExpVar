import pandas as pd


# Read DataFrame from a CSV file (adjust path as needed)
df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/veqtl_trans_v2/MAF_map_veqtl/HS_FDD_summary.csv")

# Calculate average of the 'FDD' column
fdd_mean = df["FDD"].mean()
print(f"Average FDD veqtl HS: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/MAF_map_eqtl/HS_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FDD eqtl HS: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/MAF_map_eqtl/Ctrl_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FDD eqtl Ctrl: {fdd_mean}")
'''
df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/veqtl_trans_v2/MAF_map_nonveqtl/HS_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FDD non veqtl HS: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/veqtl_trans_v2/MAF_map_nonveqtl/Ctrl_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FDD non veqtl Ctr: {fdd_mean}")
'''
df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/MAF_map_noneqtl/HS_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FDD non eqtl HS: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/MAF_map_noneqtl/Ctrl_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FDD non eqtl Ctrl: {fdd_mean}")

#Minor
df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/Minor_map_noneqtl/Ctrl_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD non eqtl Ctrl: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/Minor_map_noneqtl/HS_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD non eqtl HS: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/Minor_map_eqtl/Ctrl_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD eqtl Ctrl: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/Minor_map_eqtl/HS_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD eqtl HS: {fdd_mean}")


df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/veqtl_trans_v2/Minor_map_nonveqtl/HS_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD non veqtl HS: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/veqtl_trans_v2/Minor_map_nonveqtl/Ctrl_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD non veqtl Ctr: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/veqtl_trans_v2/Minor_map_veqtl/HS_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD veqtl HS: {fdd_mean}")

df = pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/veqtl_trans_v2/Minor_map_veqtl/Ctrl_FDD_summary.csv")
fdd_mean = df["FDD"].mean()
print(f"Average FMD veqtl Ctr: {fdd_mean}")