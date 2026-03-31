# Author: Huiting
import multiprocessing
import random
import pandas as pd
import numpy as np
import subprocess
import sys
import time
import pandas as pd
import random

"""
def get_subsample(non_QTL_df,size,df,seed,name):
    subsample=non_QTL_df.groupby('Quintile').sample(n=size, random_state=seed)
    subsample.to_csv(f'subsample_{group}_{seed}.csv', index=False)
    header_cmd = f"head -n 1 {df} > mapped_{group}_{seed}.csv"
    subprocess.run(["bash", "-c", header_cmd])
    i=0
    for snp in subsample['SNP']:
        i+=1
        print(i/len(subsample['SNP']))
        cmd = f"(grep {snp} {df} | shuf --random-source=<(yes {seed}) -n 1) >> mapped_{group}_{seed}.csv"
        subprocess.run(["bash", "-c", cmd])
    # Combine all filtered SNP dataframes into one
    #mapped_df = pd.concat(dfs, ignore_index=True)
    #map_list=[df[df['SNP'] == i].sample(n=1, random_state=seed) for i in subsample['SNP']]
    
    mapped_df = pd.read_csv(f"mapped_{group}_{seed}.csv")
    assert len(mapped_df) == len(subsample), f"Length mismatch: mapped_df={len(mapped_df)}, non_QTL_df={len(non_QTL_df)}"
    return mapped_df

def get_subsample(non_QTL_df,size,df_path,seed,name):
    subsample=non_QTL_df.groupby('Quintile').sample(n=size, random_state=seed)
    subsample.to_csv(f'MAF_map/subsample_{group}_{seed}.csv', index=False)    
    i=0
    chunk_size=2000000
    for chunk in pd.read_csv(df_path,chunksize=chunk_size):
        # Merge
        merged_chunk = pd.merge(
        chunk,
        subsample.loc[:, ["SNP"]],
        on="SNP",
        how='inner')
        output_file=f'/fml/ag-pallares_projects/NexDNAseq_G30/hxu/MAF_map/mapped_{group}_{seed}.csv'
        if i ==0:
            merged_chunk.to_csv(output_file, 
                             header=True,  # write header only for first chunk
                             index=False)
        else:
            # Write the merged chunk to file
            merged_chunk.to_csv(output_file, 
                               mode='a',  # append mode
                               header=False,  # write header only for first chunk
                               index=False)
          
        i+=1
        print(f"{chunk_size*i}/1066308130","done")
    merged_df = pd.read_csv(f"/fml/ag-pallares_projects/NexDNAseq_G30/hxu/MAF_map/mapped_{group}_{seed}.csv")
    mapped_df = merged_df.sample(frac=1, random_state=seed).drop_duplicates(subset=["SNP"])
    assert len(mapped_df) == len(subsample), f"Length mismatch: mapped_df={len(mapped_df)}, non_QTL_df={len(non_QTL_df)}"
    return mapped_df
"""

# get random SNPs info from subsample with Criteria
def get_subsample(non_QTL_df,size,df_path,seed,group,output_dir):
    # get subsample with Criteria
    subsample=non_QTL_df.groupby('Quintile').sample(n=size, random_state=seed)
    
    snp_info={x:(y,z) for (x,y,z) in zip(subsample["SNP"],subsample["MAF"],subsample["Quintile"])}
    output={x:[] for x in subsample["SNP"]}
    random.seed(seed)
    #map SNPs in the subsampke in big dataset
    with open(df_path) as f:
        header=f.readline()
        header_list = header.strip().split(",") 
        if "SNP" in header:
            snp_index = header_list.index("SNP")
            print(f'"SNP" is in column {snp_index}')
        else:
            print('"SNP" not found in header')
        
        for idx,line in enumerate(f):
            snp=line.split(",")[snp_index]
            # process
            if idx%100000==0:
                print("   ",idx,"     ",end="\r")
            if snp in output.keys():
                output[snp].append(line)
    for snp,snp_tbl_info in output.items():
        # get random sample
        output[snp]=random.choice(snp_tbl_info)
        
    output_file=f'{output_dir}/mapped_{group}_{seed}.csv'
    with open(output_file,"w") as f:
        f.write(f"SNP,MAF,Quintile,{header}")
        for snp in output:
            f.write(f"{snp},{snp_info[snp][0]},{snp_info[snp][1]},{output[snp]}")
    return pd.read_table(output_file,sep=',' )
    
   
# calculate FDI from a dataframe    
def Fraction_Derived_Increased(df, alt_col, slope):
    # Keep necessary columns only
    df = df[["SNP", "ANC", alt_col, slope]]
    # Direction of effect with respect to derived allele
    df['DER_Direction'] = np.where(df["ANC"] != df[alt_col], 
                                     df[slope], 
                                     -df[slope])

    # Determine if derived allele increases trait
    df["Derived_Increased"] = (df["DER_Direction"] > 0).astype(int)
    percent=df["Derived_Increased"].mean()
    return percent


def Fraction_Minor_Increased(df,minor_df_path, alt_col, slope):
    minor_df=pd.read_table(minor_df_path,sep= r'\s+')
    merged_df = df.merge(minor_df, on='SNP', how='left')
    # Keep necessary columns only
    df = merged_df[["SNP", "Minor", alt_col, slope]]
    # Direction of effect with respect to derived allele
    df['Direction'] = np.where(df["Minor"] == df[alt_col], 
                                     df[slope], 
                                     -df[slope])
    # Determine if derived allele increases trait
    df["Minor_Increased"] = (df["Direction"] > 0).astype(int)
    percent=df["Minor_Increased"].mean()
    return percent
    
def main ():
    seed=int(sys.argv[1])
    group=str(sys.argv[2])
    df_path=str(sys.argv[3])
    non_QTL_df=pd.read_table(str(sys.argv[4]),sep= r'\s+')
    criteria_df=pd.read_table(str(sys.argv[5]),sep= r'\s+')
    minor_df_path=str(sys.argv[6])
    output=str(sys.argv[7])
    
    
    labels = criteria_df['Quintile']
    bins = list(criteria_df['Min'])+ [criteria_df['Max'].iloc[-1]]
    bins = [float(b) for b in bins]
    
    non_QTL_df['Quintile'] = pd.cut(non_QTL_df['MAF'], bins=bins, labels=labels, right=True)

    # size of snps based on criteria
    size=list(criteria_df['NumSNPstoSample'])[0]
    
    mapped_df=get_subsample(non_QTL_df,size,df_path,seed,group,output)
    
    #fdi=Fraction_Derived_Increased(mapped_df,'ALT','COR')
    #with open(f'{output}/{group}_fdi_SEED_{seed}.txt', "w") as f:
    #   f.write(f'mapped_{group}_{seed} {fdi}')
    
    fmi=Fraction_Minor_Increased(mapped_df,minor_df_path,'ALT','COR')
    
    with open(f'{output}/{group}_fmi_SEED_{seed}.txt', "w") as f:
        f.write(f'mapped_{group}_{seed} {fmi}')
    
  

    
if __name__ == '__main__':
    main()

   


'''
    df_path=f"/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/{group}.trans_qtl_pairs_filtered.csv"
    non_QTL_df=pd.read_table(f"/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/{group}_non-trans-eQTL_downsubsample_SNPs_v2.txt",sep=r'\s+' )
    criteria=pd.read_table(f'/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/{group}_trans-eQTL_MAF_Criteria_v2.txt',sep=r'\s+' )
    df_path=f"/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/{group}_trans_eqtl_sig_filtered.csv"
    non_QTL_df=pd.read_table(f"/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/{group}_trans-eQTL_downsubsample_SNPs_v2.txt",sep=r'\s+' )
    criteria=pd.read_table(f'/tmp/global2/jtanshengyi/hxu_test/sampling/eqtl_trans/new/{group}_trans-eQTL_MAF_Criteria_v2.txt',sep=r'\s+' )
    
    # Build list of argument tuples for each call
    seeds = list(range(10)) 
    args_list = [(non_QTL_df, size, df_path, s, group) for s in seeds]
    with multiprocessing.Pool(processes=len(seeds)) as pool:
        results = pool.starmap(FDI_cal_subsample, args_list)
    df = pd.DataFrame(results, columns=['subsample', 'FDI'])
    df.to_csv(f'MAF_map/{group}_fdi_result.csv', index=False)
'''