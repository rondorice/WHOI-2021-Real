#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install git+https://github.com/SciTools/cartopy.git')
get_ipython().system('pip uninstall shapely -y')
get_ipython().system('pip install shapely --no-binary shapely --force')

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature


# In[16]:


import matplotlib.pyplot as plt
import pandas as pd

#Reading in data and grouping
micro_init = pd.read_csv('Desktop/rename_m.merged.numreads.csv', sep='\t')
micro_init['Name'] = micro_init['Name'].str.split('_').str[:1].str.join('_')
micro_init = micro_init.groupby(['Name'], axis=0, as_index=False).sum()
#Transpose
micro_init_t = micro_init.transpose()
micro_init_t.reset_index(level=0, inplace=True)
new_header = micro_init_t.iloc[0]
micro_init_t = micro_init_t[1:]
micro_init_t.columns = new_header
#Rename to improve legend
micro_init_t.columns = ['run_accession','Micromonas Polaris','Micromonas Polaris','Micromonas Polaris','Micromonas Brave','Micromonas Puscilla','Micromonas Commoda']
micro_init['Name'] = micro_init['Name'].str.replace('MMETSP0802', 'M. polaris [MMETSP0802]')
micro_init['Name'] = micro_init['Name'].str.replace('MMETSP1327', 'M. polaris [MMETSP1327]')
micro_init['Name'] = micro_init['Name'].str.replace('MMETSP1390', 'M. polaris [MMETSP1390]')
micro_init['Name'] = micro_init['Name'].str.replace('MMETSP1401', 'M. bravo [MMETSP1401]')
micro_init['Name'] = micro_init['Name'].str.replace('MicpuC3v2', 'M. puscilla [MicpuC3v2]')
micro_init['Name'] = micro_init['Name'].str.replace('MicpuN3v2', 'M. commoda [MicpuN3v2]')

#Pie chart legend for mapping data
plt.figure(figsize=(16,8))
ax1 = plt.subplot(121, aspect='equal')
micro_init.plot(kind='pie', y = 'ERR868508', ax=ax1, autopct='%1.1f%%', 
startangle=90, shadow=False, labels=micro_init['Name'], legend = True, fontsize=14)


# In[67]:


# Read in other data and merge df
other_data = pd.read_excel('Desktop/PRJEB4352_metaG_wenv_PE.xlsx')
other_data
good = other_data[['run_accession','Latitude','Longitude']]
merged_df = good.merge(micro_init_t, how = 'inner', on =['run_accession'])
#merged_df.set_index(['run_accession','Latitude','Longitude'], inplace= True)
merged_df


# In[68]:


#Micromonas mapping
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
for index, row in merged_df.iterrows():
    lon,lat = merged_df.iloc[index]['Longitude'],merged_df.iloc[index]['Latitude']
    lonr,latr =  ccrs.Robinson().transform_point(lon,lat, ccrs.PlateCarree())
    plot_pie_inset([merged_df.iloc[index]['MMETSP0802'],merged_df.iloc[index]['MMETSP1327'],
                       merged_df.iloc[index]['MMETSP1390'],merged_df.iloc[index]['MMETSP1401'],
                       merged_df.iloc[index]['MicpuC3v2'],merged_df.iloc[index]['MicpuN3v2']],lonr,latr,ax,0.2)


plt.show()


# In[ ]:




