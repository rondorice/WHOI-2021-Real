#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Get pandas and cartopy
import pandas as pd
import numpy as np

get_ipython().system('pip install git+https://github.com/SciTools/cartopy.git')
get_ipython().system('pip uninstall shapely -y')
get_ipython().system('pip install shapely --no-binary shapely --force')

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature


# In[62]:


#Read in data and trim down taxonomy
ev_params = pd.read_csv('environmental_parameters.csv',sep='\t')
ab_matrix = pd.read_csv('abundance_matrix.csv',sep='\t')
ab_matrix['taxonomy'] = ab_matrix['taxonomy'].str.rsplit(';').str[-1] 
ab_matrix = ab_matrix.groupby(['taxonomy']).sum()
#Transpose the data and rename columns
df_tr = ab_matrix.T
df_tr = df_tr.reset_index()
df_tr = df_tr.rename(columns={"index":'Sample_ID'})
#Merge
df_welp = pd.merge(df_tr, ev_params, on=['Sample_ID'], how='left',indicator='Exist')
#Removing unnecessary rows of DF (probably a better way to do this, but alas)
df_welp.drop(df_welp.iloc[:, 1:21], inplace = True, axis = 1)
df_welp.drop(df_welp.iloc[:, 8:19],inplace = True, axis = 1)
df_welp.drop(df_welp.iloc[:, 9:10],inplace = True, axis = 1)
df_welp.drop(df_welp.iloc[:, 12:84],inplace = True, axis = 1)
df_welp.drop(df_welp.iloc[:, 0:1],inplace = True, axis = 1)
df_welp.drop(df_welp.iloc[:, 8:9],inplace = True, axis = 1)

#Aggregate and sum
aggregation = {' Micromonas':'sum',' Micromonas+Clade-A.ABC.12':'sum',' Micromonas+Clade-B.4':'sum',
            ' Micromonas+Clade-B.E.3':'sum',' Micromonas+Clade-B.X':'sum',' Micromonas+Clade-B_arctic':'sum',
            ' Micromonas+Clade-C.D.5':'sum','Latitude':'first','Longitude':'first'}
df_agg = df_welp.groupby(df_welp['Station']).aggregate(aggregation)
df_agg = df_agg.reset_index()
#Rename the columns
df_agg = df_agg.rename(columns={" Micromonas+Clade-B.E.3":'Micromonas Bravo'," Micromonas+Clade-A.ABC.12":'Micromonas commoda',
                               ' Micromonas+Clade-B_arctic':'Micromonas Polaris',' Micromonas+Clade-C.D.5':'Micrmonas pusilla'})
df_agg


# In[5]:


# Make a pie chart to get a legend
transposed = df_agg.T
transposed = transposed.drop(['Latitude','Longitude',' Micromonas'])
new_header = transposed.iloc[0] #Grab the first row for the header
transposed = transposed[1:] #Take the data less the header row
transposed.columns = new_header #Set the header row as the df header
transposed = transposed.reset_index()
transposed
plt.figure(figsize=(16,8))
ax1 = plt.subplot(121, aspect='equal')
transposed.plot(kind='pie', y = 'TARA_006', ax=ax1, autopct='%1.1f%%', 
startangle=90, shadow=False, labels=transposed['index'], legend = True, fontsize=14)


# In[6]:


#All size fractions- mapping 18s as a test
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

ax = plt.axes(projection=ccrs.Robinson())
ax.coastlines(resolution='110m')
ax.stock_img() 


def plot_pie_inset(data,ilon,ilat,ax,width):
    ax_sub= inset_axes(ax, width=width, height=width, loc=10, 
                       bbox_to_anchor=(ilon, ilat),
                       bbox_transform=ax.transData, 
                       borderpad=0)
    wedges,texts= ax_sub.pie(data)
    
    ax_sub.set_aspect("equal")
#for column in merged_df.columns[0:1]:
for index, row in df_agg.iterrows():
    lon,lat = df_agg.iloc[index]['Longitude'],df_agg.iloc[index]['Latitude']
    lonr,latr =  ccrs.Robinson().transform_point(lon,lat, ccrs.PlateCarree())
    plot_pie_inset([df_agg.iloc[index][' Micromonas+Clade-A.ABC.12'],
                       df_agg.iloc[index][' Micromonas+Clade-B.4'],df_agg.iloc[index][' Micromonas+Clade-B.E.3'],
                       df_agg.iloc[index][' Micromonas+Clade-B.X'],df_agg.iloc[index][' Micromonas+Clade-B_arctic'],
                       df_agg.iloc[index][' Micromonas+Clade-C.D.5']],lonr,latr,ax,0.2)


# In[29]:


## MetaG Mapping for micromonas
metaG_df = pd.read_csv('mappingrates.csv')
fullMG_df = pd.read_excel('PRJEB4352_metaG_wenv_PE.xlsx')

fullMG_df.drop(fullMG_df.columns.difference(['run_accession','Latitude','Longitude','OS_region']), 1, inplace = True)
metaG_df = metaG_df.rename(columns={'directorylist.txt':'run_accession','[2021-07-01 23:32:59.160] [jointLog] [info] Mapping rate = 0.350353%':'New_MR'})
metaG_df['New_MR'] = metaG_df['New_MR'].str.split('=').str[1]
metaG_df['New_MR'] = metaG_df['New_MR'].str.split('%').str[0]

cond = ~fullMG_df['run_accession'].isin(metaG_df['run_accession'])
fullMG_df.drop(fullMG_df[cond].index, inplace = True)
fullMG_df = fullMG_df.reset_index(drop = True)

metaG_df['Longitude'] = pd.Series(fullMG_df['Longitude'])
metaG_df['Latitude'] = pd.Series(fullMG_df['Latitude'])

metaG_df = metaG_df.dropna()
metaG_df['New_MR'] = pd.to_numeric(metaG_df['New_MR'])
#metaG_df
#pd.set_option('display.max_rows', 500)
metaG_df


# In[86]:


df_micro_full = pd.merge(metaG_df, df_agg, on=['Longitude','Latitude'], how='left',indicator='Exist')
df_micro_full = df_micro_full.dropna()
df_micro_full['Total'] = df_micro_full['Total']
df_micro_full['18s/mG'] = df_micro_full['Total'] / df_micro_full['New_MR']
df_micro_full.replace([np.inf, -np.inf], np.nan, inplace=True)
df_micro_full = df_micro_full.dropna()
df_micro_full = df_micro_full.reset_index()
df_micro_full


# In[85]:


# Read in metaT df for micromonas and lat, long
metaT_df = pd.read_csv('mappingrates_t.csv')
fullMT_df = pd.read_excel('PRJEB6609_metaT_wenv_PE.xlsx')

# Get only needed columns and rename them, split names as needed
fullMT_df.drop(fullMT_df.columns.difference(['run_accession','Latitude','Longitude']), 1, inplace = True)
metaT_df = metaT_df.rename(columns={'directorylist.txt':'run_accession','[2021-07-01 22:37:53.234] [jointLog] [info] Mapping rate = 0.442336%':'New_MR'})
metaT_df['New_MR'] = metaT_df['New_MR'].str.split('=').str[1]
metaT_df['New_MR'] = metaT_df['New_MR'].str.split('%').str[0]

# Find only the needed rows
cond = ~fullMT_df['run_accession'].isin(metaT_df['run_accession'])
fullMT_df.drop(fullMT_df[cond].index, inplace = True)
fullMT_df = fullMT_df.reset_index(drop = True)

metaT_df['Longitude'] = pd.Series(fullMT_df['Longitude'])
metaT_df['Latitude'] = pd.Series(fullMT_df['Latitude'])

metaT_df = metaT_df.dropna()
metaT_df['New_MR'] = pd.to_numeric(metaT_df['New_MR'])
metaT_df


# In[93]:


# Merge into metaT df with 18s and metaT data
df_micro_full_mT = pd.merge(metaT_df, df_agg, on=['Longitude','Latitude'], how='left',indicator='Exist')
df_micro_full_mT = df_micro_full_mT.dropna()
df_micro_full_mT['Total'] = df_micro_full_mT['Total']
df_micro_full_mT['18s/mT'] = df_micro_full_mT['Total'] / df_micro_full_mT['New_MR']
df_micro_full_mT.replace([np.inf, -np.inf], np.nan, inplace=True)
df_micro_full_mT = df_micro_full_mT.dropna()
df_micro_full_mT = df_micro_full_mT.reset_index()
df_micro_full_mT


# In[107]:


#Plot metaT vs 18s data
import seaborn
from  matplotlib import pyplot

seaborn.set(style='ticks')


fg = seaborn.FacetGrid(data=df_micro_full_mT, hue='Station')
fg.map(pyplot.scatter, 'Total','New_MR').add_legend()
#fg.set(xscale='log')
#fg.set(yscale='log')
fg.fig.suptitle('MetaT vs 18s')
fg.set(xlabel='18s', ylabel='MetaT')


# In[106]:


#Plot metaG vs 18s data
import seaborn
from  matplotlib import pyplot

seaborn.set(style='ticks')


fg = seaborn.FacetGrid(data=df_micro_full, hue='Station')
fg.map(pyplot.scatter, 'Total','New_MR').add_legend()
#fg.set(xscale='log')
#fg.set(yscale='log')
fg.fig.suptitle('MetaG vs 18s')
fg.set(xlabel='18s', ylabel='MetaG')


# In[ ]:




