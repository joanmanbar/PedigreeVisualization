#!/usr/bin/env python
# coding: utf-8




# # Visualizing Pedigree

# Created by Joan Barreto Ortiz (March 2022)
# email: jbarreto@umn.edu




# # Imports/Settings

import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from collections import Counter
import seaborn as sns

# import session_info
# session_info.show()



# # Make sure Jupyter Notebook shows all outputs from the same cell
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = "all"

# # Set max number of rows to visualize
# pd.set_option('display.max_rows', 10)

# # plt.rcParams['figure.figsize'] = [12, 8]
# plt.rcParams['figure.dpi'] = 80 


# List of 100 custom parent names
P100 = [ "P0" + str(int(s)+1) for s in range(9)]
P100 = P100 + [ "P" + str(int(s)+1) for s in range(9,99)]
# P100




# ## Custom Functions

# ### Remove Backcrosses (`removeBC`)
# `parents = True` returns the pedigree without backcrosses. \
# `parents = False` returns the backcrosses
# 
# **NOTE** = This will not work well if there has been more than 9 backcrosses to a recurrent parent, i.e. _A/B*10_. \
# Additional work is necessary for double-digit backcrosses.

def removeBC(s, parents = True):
    s = s.replace('//', '/2/')
    BC_pattern = re.compile(r"/{1}\d{1}\*{1}|\*{1}\d{1}")
    # \d{1} matches only one ({1}) digit
    # \* matches only one asterisk
    # | means OR
    
    if parents == True:
        s2 = re.sub(BC_pattern, "/", s)
        s2 = s2.replace('//', '/')
        s2 = s2.rstrip('/')
    else:
        s2 = re.findall(BC_pattern, s)
        s2 = [x.replace('/', '') for x in s2]
    
    return(s2)




# ### Combine Two Dictionaries (`MergeDD`)

def MergeDD(D1,D2):
    Merged_D = {k: D1.get(k, 0) + D2.get(k, 0) for k in set(D1) | set(D2)}
    return(Merged_D)

# D1 = {'B': 0.39, 'A': 0.39}
# D2 = {'B': 0.21, 'C': 0.39, 'F': 0.10}
# MergeDD(D1,D2)




# ### Simplest Cross (`SimpleCross`)

# Function for simplest cross (A/B)
# Input is a string, residual R, parent dict
# Output is a dictionary Parents:Values

def SimpleCross(s, R):
    ThisD = dict()
    # List of primary crosses
    List_PC = s.split("/")
    List_PC = sorted(List_PC, key=len)
    # Recurrent is a longer string
    NonRecurrent = List_PC[0]
    Recurrent = List_PC[1]
    # Number of backcrosses
    n = len(Recurrent)/3
    
    # Recurrent's contribution
    Rec_Cont =  ((2**n)-1) / (2**n)
    NonRec_Cont = 1 - Rec_Cont

    # Rename keys
    NonRecurrent = NonRecurrent[:3]
    Recurrent = Recurrent[:3]
    # Add values
    ThisD[Recurrent] = R * Rec_Cont
    ThisD[NonRecurrent] = R * NonRec_Cont
        
    return(ThisD)

# # Examples
# s = 'P01/P02P02P02'
# SimpleCross(s, R=.5)




# ### Double Cross (`DoubleCross`)
# This is used when a pedigree has more than one cross. In other words, it's not 'A', or 'A/B'.

def DoubleCross(s, R):
       
    D2 = dict()
    CH = s.split('/2/')
    CH.sort(key=lambda x: x.count('/'))
    
    # count can be either 3 or 4
    SlashCount = s.count('/')
    if SlashCount==4:
        SubR = R/len(CH)
    
    for s2 in CH:
        
        try:
            SubR
        except:
            SubR=R-sum(D2.values())
            
        if s2.count('/')==0:
            n = len(s2)/3
            Rec_Cont =  ((2**n)-1) / (2**n)
            CH_R = SubR * Rec_Cont
            D1 = dict()
            D1[s2[:3]] = CH_R
            D2 = MergeDD(D1,D2)
        elif s2.count('/')==1: 
            D1 = SimpleCross(s2, SubR)
            D2 = MergeDD(D1,D2)
        else:
            print('Error! \n Neither single nor simple cross.')
            break  
        SubR=R-sum(D2.values())
        
    return(D2)

# # Examples:
# s = 'P01/P02/2/P03'
# s = 'P06/P07P07/2/P08P08P08'
# DoubleCross(s, R=1)




# ### Preprocess Pedigree (`PreProcess`)
# This creates transformed string (`s`), and three dictionaries: values, keys, and parents

def PreProcess(pedigree):
    
    pedigree = pedigree.replace('//','/2/')
    s = pedigree
#     Simple_CP = re.compile('P\d{2}/{1}P\d{2}')
    
    # Clean string from backcrosses
    s = removeBC(s, parents=True)    
    
    # Create dict with parents
    List_parents = [p for p in s.split("/") if p.isdigit()==False]
    List_parents = list(filter(None, List_parents))
    D_parents = dict((p,0) for p in List_parents)
    # Change original names
    D_key = zip(List_parents,P100[:len(List_parents)])
    D_key = dict(D_key)
    
    D_values = dict((D_key[key], value) for (key, value) in D_parents.items())
    
    pattern = re.compile('|'.join(re.escape(key) for key in D_key.keys()))
    result = pattern.sub(lambda x: D_key[x.group()], s)
    
    s = result
            
    # Find backcrosses
    All_L = [p for p in pedigree.split("/")]
    Backcrosses = [x for x in All_L if "*" in x]
    Backcrosses = '/'.join(Backcrosses)
    Backcrosses = '/' + Backcrosses + '/'
    Backcrosses = removeBC(Backcrosses,False)
    Backcrosses = [x.replace('*', '') for x in Backcrosses]

    BC_index = [ind for ind,x in enumerate(All_L) if "*" in x]
    
    # Add backcross to new string
    Transformed_L = [p for p in s.split("/")]
    
    for ind,val in enumerate(BC_index):
        Transformed_L[val] = Transformed_L[val]*int(Backcrosses[ind])
    
    # Concatenate    
    s = '/'.join(Transformed_L)
    s
    
    return(s, D_values, D_key, D_parents)




# ### Ancestral Contributions (`AnCon`)

def AnCon(pedigree, R=1, SimplifyPed=True, ShowPedigree=True):
    
    R_updated = R
    
    # Preprocessing
    s, D_values, D_key, D_parents = PreProcess(pedigree)
    Transformed_s = s
    
    # Get number of crosses (after second crosses)
    List_CrossesASC = [x for x in s.split("/") if x.isdigit()]
    
    if len(List_CrossesASC)==0 and s.count('/')==0:
        D_values[Transformed_s] = R_updated
        s = ''
    
    while len(s) > 0:

        # Get max cross number, and the previous cross
        if len(List_CrossesASC) > 0:
            MAX = max(List_CrossesASC)
            NextCross = "/" + str(int(MAX)-1) + "/"
            MaxCross = "/" + MAX + "/"
        else:
            MAX = '1'
            MaxCross = "/" + MAX + "/"
            List_CrossesASC.append(MAX)
        
        # Separate female and male
        Fem_Male = s.split(MaxCross)
        Fem_Male.sort(key=lambda x: x.count('/'))
        
        # Current half
        CH = Fem_Male[0]
        SlashCount = CH.count('/')
        
        if len(Fem_Male) == 1:
            s = ''
        else:
            s = s
        
        # No cross
        if SlashCount == 0:
            n = len(CH)/3
            Rec_Cont =  ((2**n)-1) / (2**n)
            CH_R = R_updated * Rec_Cont
            D1 = dict()
            D1[CH[:3]] = CH_R
            D_values = MergeDD(D1,D_values)
            List_CrossesASC.remove(MAX)
            R_updated = R - sum(D_values.values())
            
            try: 
                s = Fem_Male[1]
            except:
                break
    
        # Simple cross        
        elif SlashCount == 1:
            CH_R = (R_updated/len(Fem_Male))
            D1 = SimpleCross(CH, CH_R)
            D_values = MergeDD(D1,D_values)
            List_CrossesASC.remove(MAX)
            R_updated = R - sum(D_values.values())
            
            try: 
                s = Fem_Male[1]
            except:
                break
            
        # Check for double cross 
        elif SlashCount==3 or SlashCount==4:
            D1 = DoubleCross(CH,R_updated/len(Fem_Male))
            D_values = MergeDD(D1,D_values)
            R_updated = R - sum(D_values.values())
            List_CrossesASC.remove(MAX)
            try : 
                s = Fem_Male[1]
            except:
                break    
        else:
            print("Not currently available for higher number of repeated crosses")
            break
             
                
    if SimplifyPed == True:
        OUTPUT = D_values
        if ShowPedigree == True:
            print(Transformed_s)
        
    else:
        for key, value in D_key.items():
            D_parents[key] = D_values[value]
        OUTPUT = D_parents
        if ShowPedigree == True:
            print(pedigree)
    
    return(OUTPUT)



# # Exmaples:
    
# pedigree = 'A/*3B//C'    
# AnCon(pedigree, R=1, SimplifyPed=True)

# pedigree = 'II53.388/ANDES//YAKTANA 54/NORIN 10 BREVOR/3/LERMA ROJO/4/B4946.A.4.18.2.1Y/YAQUI 53//3*YAQUI 50'
# pedigree = 'P01/P02/2/P03/P04'
# pedigree='A/5/2*B/C5//D/3/E//F5*3/G/4/H'
# pedigree = 'AGATHA/3*YECORA_F_70'
# pedigree = 'A/B*4//3*C/3/D/C/4/F*3'
# pedigree = 'A/5/B/C//D/3/E//F/G/4/H'
# pedigree = 'A/5/B/C//D/3/F/G/4/H'
# pedigree = 'A/B/2/C/3/D/4/E/F/2/G/H'
# pedigree = 'A'
# pedigree = 'PENJAMO T 62/GABO 55'
# Dtest = AnCon(pedigree, R=1, SimplifyPed=False)
# Dtest
# sum(Dtest.values())




# ### Visualize Ancestral Contributions (`PedigreeVis`)
# `Log = False` plots the original contributions. \
# `Log = True` plots the Log10 contributions. This is useful to visualize which parents have been used more often in the current pedigree.

# Creating autocpt arguments
def func(pct, allvalues):
    absolute = int(pct / 100.*np.sum(allvalues))
    return "{:.1f}%".format(pct, absolute)

def PedigreeVis(DD, Log = False):

    # Data to plot
    labels = list(DD.keys())
    
    if Log == True:
        DD.update((x, y+1) for x, y in DD.items())
        sizes = np.log10(list(DD.values()))
        Min = abs(min(sizes))
        sizes = [s+Min for s in sizes]
    else:
        sizes = list(DD.values())
    
    # Creating plot
    fig, ax = plt.subplots(figsize =(10, 7))
    wedges, texts, autotexts = ax.pie(sizes,
                                      autopct = lambda pct: func(pct, sizes),
                                      textprops = dict(color ="white"))

    # Adding legend
    ax.legend(wedges, labels,
              title ="Parents",
              loc ="center left",
              bbox_to_anchor =(1, 0, 0.5, 1))

    plt.setp(autotexts, size = 10, weight ="bold")
    ax.set_title("Expected Genetic Contribution", color = 'white')

    # show plot
    plt.show()

# pedigree = P['Pedigree'][8]
# pedigree = 'A/B/2/C/3/D*2/4/E/F/2/G/H*5'
# DD = AnCon(pedigree, R=1, SimplifyPed=False)
# PedigreeVis(DD, Log = True)




# ### Obtain Unique Parents  (`UniqueParents`)
# This creates a list with unique parents in pedigrees. 
# It requires a dataframe with Pedigree and GID as columns

def UniqueParents(df):
    
    # Expand all parents
    Mdf = df.copy()
    Mdf["Pedigree"] = Mdf["Pedigree"].apply(lambda x: removeBC(str(x)))
    Mdf = Mdf["Pedigree"].str.split('/', expand = True)
    Mdf = Mdf.fillna(value=np.nan)
    
    # Get all potential parents
    Parents = []
    # Iterate over the index range from
    # 0 to max number of columns in dataframe
    for Col in range(Mdf.shape[1]):
        # Select column by index position using iloc[]
        L = list(Mdf.iloc[:,Col])
        Parents.extend(L)
        
    # Count occurrences
    D = Counter(Parents)
    # D.most_common()
    D.pop(np.nan)
    
    # create list with elements to remove
    L2R = ['','2','3','4','5','6','7','8','9','10',
           '11','12','13','14']

    # collection is not exactly a dict type. That's why dict(D)
    for key in dict(D).keys():
        if key in L2R:
            D.pop(key)
    
    Dkeys = list(D.keys())
    UP = list(np.unique(Dkeys))
    
    # Return a list of unique parents (UP)
    return(UP)

# UP = UniqueParents(df=P)



# ### Get ancestral contribution for all GID in dataframe  (`df_AnCon`)
# This creates a df with GID as rows and unique parents as columns. 
# It requires a dataframe with Pedigree and GID as columns.

def df_AnCon(df):
    
    P2 = df.copy()
    P2 = P2.reset_index()
    Max = max(P2['index']) + 1 
    
    # Get colnames
    UP = UniqueParents(df)
    
    # Create empty df with UP
    df_all = pd.DataFrame(columns=UP)
    
    for i in range(Max):
        # Get ancestral contribution for current row
        D = AnCon(P2['Pedigree'][i], R=1, SimplifyPed=False, ShowPedigree=False)
        # Add data to dictionary using GID as index
        df_D = pd.DataFrame(D, index=[P2['GID'][i],])
        # Update df and fill na with 0
        df_all = pd.concat([df_all, df_D])
        df_all = df_all.fillna(0)

    
    return(df_all)


# df1 = df_AnCon(df=P)
# # Verify that sum of contributions is exactly 1
# df1['SUM'] = df1.sum(axis=1)
# df1.loc[df1['SUM'] != 1]




# ### Plot Relationships between GIDs/Parents (`PlotRelationship`)
# `Log = False` plots the original contributions. \
# `Log = True` plots the Log10 contributions. This is useful to visualize which parents have been used more often in the current pedigree.

def PlotRelationship(df, Cov=False, PlotSize=(12,8), LW=0.2, LabSize=6):
        
    if Cov==True:
        df = df.T.cov()
    else:
        df = df
    
    # sns.set(font_scale=0.2)
    fig, ax = plt.subplots(figsize=PlotSize)
    ax = sns.heatmap(df, linewidths=LW, 
                     xticklabels=True, yticklabels=True,
                    cmap='YlOrBr')
    ax.xaxis.tick_top() # x axis on top
    ax.xaxis.set_label_position('top')
    plt.tick_params(labelsize=LabSize)
    plt.xticks(rotation=90)
    plt.show()
    
    
    
