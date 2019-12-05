#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import numpy as np
import hypertools as hyp
import seaborn as sns
from timeit import default_timer as timer
import os
import altair as alt


get_ipython().run_line_magic('matplotlib', 'inline')


# In[3]:


import warnings
warnings.filterwarnings('ignore')


# In[4]:


cdr = pd.read_csv('/home/msouza/mcs/study/code/bioinfo/data/analysis/attila_cdr/Final_wHAIPI009594-32_1aafreq.csv').sort_values(by='quantity', ascending=False)


# In[5]:


cdr.head()


# In[6]:


c = cdr.iloc[:1000]


# In[7]:


c.head()


# In[8]:


alt.Chart(c).mark_point().encode(
    x=alt.X('length'),
    y=alt.Y('quantity')
)


# In[9]:


sequence = c.pop('cdr3')


# In[10]:


cluster_labels = hyp.cluster(c, n_clusters=10)


# In[11]:


hyp.plot(c, '.', reduce='TSNE', ndims=3, group=cluster_labels, legend=None, title='FinalRound_R4_S4_L001_R1_001aafreq')


# In[12]:


hyp.plot(c, '.', reduce='TSNE', ndims=2, group=cluster_labels, legend=None, title='FinalRound_R4_S4_L001_R1_001aafreq');


# In[ ]:


hyp.plot(cdr, '.', reduce='TSNE', ndims=3, group=cluster_labels, legend=None, title='FinalRound_R4_S4_L001_R1_001aafreq')


# In[ ]:




