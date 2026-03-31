import pandas as pd
df=pd.read_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/FDI_FDD_percent.csv", index_col=0)
df["FDD_withoutbetween"]=df["FDD_all"]/(df["FDD_all"]+df["FDI_all"])
df["FDI_withoutbetween"]=df["FDI_all"]/(df["FDD_all"]+df["FDI_all"])
df.to_csv("/tmp/global2/jtanshengyi/hxu_test/sampling/FDD_FDI_withoutbetween.csv", index=False)

