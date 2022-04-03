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

D1 = {'B': 0.39, 'A': 0.39}
D2 = {'B': 0.21, 'C': 0.39, 'F': 0.10}
MergeDD(D1,D2)




# ### Simplest Cross (`SimpleCross`)

# Function for simplest cross (A/B)
# Input is a string, residual R, parent dict
# Output is a dictionary Parents:Values

def SimpleCross(s, R):
    ThisD = dict()
    # List of primary crosses
    List_PC = s.split("/")
    LastMale = List_PC[1]
    LastFemale = List_PC[0]
    # Contributions
    Ratio_M = len(LastMale)/3
    Ratio_F = len(LastFemale)/3
    Total = Ratio_M + Ratio_F
    # Rename keys
    LastMale = LastMale[:3]
    LastFemale = LastFemale[:3]
    # Add values
    ThisD[LastMale] = R * Ratio_M/Total
    ThisD[LastFemale] = R * Ratio_F/Total
    
#     print(ThisD[LastFemale] + ThisD[LastMale])
    
    return(ThisD)




# ### Single or No Cross (`Simp_Doub`)

def Simp_Doub(s,R):
    if s.count('/') == 0:
        Ratio_s = len(s)/3
        s = s[:3]
        D1 = dict()
        D1[s] = R - (R/2)**Ratio_s
    elif s.count('/') == 1:
        D1 = SimpleCross(s,R/2)
    else: 
        print('Error!')
        
    return(D1)




# ### Double Cross (`DoubleCross`)
# This is used when a pedigree has more than one cross. In other words, it's not 'A', or 'A/B'.

def DoubleCross(s, R):

    D2 = dict()
    CH = s.split('/2/')
    CH.sort(key=lambda x: x.count('/'))
    SubR = R/len(CH)
    
    for s2 in CH:
        D1 = Simp_Doub(s2,SubR)
        D2 = MergeDD(D1,D2)
        
    return(D2)




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

def AnCon(pedigree, R=1, SimplifyPed=True):
    
    # Preprocessing
    s, D_values, D_key, D_parents = PreProcess(pedigree)
    Transformed_s = s
    
    # Get number of crosses (after second crosses)
    List_CrossesASC = [x for x in s.split("/") if x.isdigit()]    
    
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
        
        # Check for single or no cross
        if SlashCount < 2:
            D1 = Simp_Doub(CH,R)
            # Parent Ratio
            PR = sum(D1.values())
            if PR>1:
                # Contribution
                C = R * (PR/(PR+1))
                D1 = D1.update((x, C) for x, y in D1.items())
            else:
                R = R - PR
                
            D_values = MergeDD(D1,D_values)
            List_CrossesASC.remove(MAX)
            
            try : 
                s = Fem_Male[1]
            except:
                break
            
        # Check for double cross 
        elif SlashCount == 3:
            D1 = DoubleCross(CH,R)
            PR = sum(D1.values())
            D_values = MergeDD(D1,D_values)
            R = R - PR
            List_CrossesASC.remove(MAX)
            try : 
                s = Fem_Male[1]
            except:
                break
            
        # Break otherwise    
        else:
            break  
    # If there is no other parent, add D1 again
    if len(s)==0:
        D_values = MergeDD(D1,D_values)
    
    print('Sum of contributions = ', sum(D_values.values()))
    
    # Return original names if required        
    if SimplifyPed == True:
        print(Transformed_s)
        return(D_values)
    else:
        for key, value in D_key.items():
            D_parents[key] = D_values[value]
        print(pedigree)
        return(D_parents)




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
        sizes = np.log10(list(DD.values()))
        Min = abs(min(sizes))
        sizes = [s+Min for s in sizes]
#         sizes
    else:
        sizes = list(DD.values())
    
    # Creating plot
    fig, ax = plt.subplots(figsize =(10, 7))
    wedges, texts, autotexts = ax.pie(sizes,
                                      autopct = lambda pct: func(pct, sizes),
                                      textprops = dict(color ="white"),
                                      normalize=False)

    # Adding legend
    ax.legend(wedges, labels,
              title ="Parents",
              loc ="center left",
              bbox_to_anchor =(1, 0, 0.5, 1))

    plt.setp(autotexts, size = 10, weight ="bold")
    ax.set_title("Expected Genetic Contribution", color = 'white')

    # show plot
    plt.show()



