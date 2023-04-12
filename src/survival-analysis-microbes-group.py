#!/usr/bin/env python
# coding: utf-8

# In[4]:


import pandas as pd
import lifelines
from lifelines import KaplanMeierFitter
       
OS_data = pd.read_table("survival.data.tsv")
OS_result = open("OS_result.txt","w")
PFS_result = open("PFS_result.txt","w")
with open("projects_list.txt","r") as cancer_list:
    for cancer_name in cancer_list:
        cancer_name = cancer_name.rstrip("\n")
        cancer_data = OS_data.loc[OS_data["cancer"]== cancer_name]
    
        with open("taxa_list.txt","r") as taxa_list:
            for taxa_name in taxa_list:
                taxa_name = taxa_name.rstrip("\n")
                cancer_median = cancer_data[taxa_name].median()
                
                ### classify group
                def get_group(x):
                     if x <= cancer_median:
                        return "low"
                     else: return "high"
                cancer_data[taxa_name +"_group"] = cancer_data[taxa_name].apply(get_group)
                
                ## survival analysis
                kmf = KaplanMeierFitter()
                OS_T = cancer_data["OS_time"]
                OS_E = cancer_data["OS"]
                
                PFS_T = cancer_data["PFS_time"]
                PFS_E = cancer_data["PFS"]
                
                dem = (cancer_data[taxa_name + "_group"] == "high")
                
                OS_sur = logrank_test(OS_T[dem], OS_T[~dem],OS_E[dem], OS_E[~dem], alpha=.99)
                PFS_sur = logrank_test(PFS_T[dem], PFS_T[~dem],PFS_E[dem], PFS_E[~dem], alpha=.99)
                
                ## find p value
                OS_p = str(OS_sur.p_value)
                PFS_p = str(PFS_sur.p_value)
                
                print(cancer_name + "：" + taxa_name + ":" , OS_p, file=OS_result)
                print(cancer_name + "：" + taxa_name + ":" , PFS_p, file=PFS_result)

                cancer_data.to_csv( cancer_name + "_group"+".csv")
OS_result.close()
PFS_result.close()
               