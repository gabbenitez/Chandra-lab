#!/usr/bin/env python
# coding: utf-8

# 
# # Itzkovitz Population Mixture Model
# 
# 
# #### http://shalevlab.weizmann.ac.il/wp-content/uploads/2016/02/telomeres.pdf
# 

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy
from scipy.optimize import curve_fit
import math
import pandas as pd
import os


# 
# ## Average telomere length
# 

# In[2]:


#we set up the equation
def avg_tel(t, M, delta, alpha):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# In[3]:


# set up time
t = np.linspace (18, 100, 1000)

# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

GS_age_DNAmTL = df[["age", "DNAmTL"]]
x_GS = GS_age_DNAmTL.age
y_GS = GS_age_DNAmTL.DNAmTL


# In[4]:


# for playing around with

#parameters
L0 = 8.5
M = 1000
delta = 0.00002
alpha = 0.02

# plotting scatter of GS data vs model
plt.plot(t, avg_tel(t, M, delta, alpha,), color = 'red')
plt.scatter(x_GS, y_GS, s = 2)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# ---

# 
# ### Find LO based on dataset
# 

# No scientific basis for this.... just estimating L0 based on being higher than initial values of earliest data??? bc with L0 as a paramter, results are funky...
# 
# 
# How best to infer L0 from our data? also given that it starts from a non-0 age...

# In[5]:


# function to find L0, as an average of y values corresponding to first x values 
def find_L0(df):
    sorted_age_DNAmTL = df.sort_values('age')
    temp = sorted_age_DNAmTL.head(300)
    L0 = temp['DNAmTL'].median()
    L0 = L0 + 0.4
    return L0


# 
# ## Curve Fitting
# 

# Should constrain bounds for parameters based on scientific evidence/intuition but just for this sake will let python do its thing
# 
# Also weird that L0 if it gets too low, with these parameters, the curve of the function flips

# #### Fit to GS data

# In[6]:


#get L0 for dataset
L0 = find_L0(GS_age_DNAmTL)
L0


# In[7]:


GS_pars, GS_testcov = curve_fit(avg_tel, x_GS, y_GS)
print(GS_pars)


# In[8]:


plt.plot(t, avg_tel(t, *GS_pars), color = 'red')
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# #### GS fit to LBC data

# In[9]:


df = pd.read_csv('../Data/sheets/LBC.csv', index_col=0)

LBC_age_DNAmTL = df[["age", "DNAmTL"]]
x_LBC = LBC_age_DNAmTL.age
y_LBC = LBC_age_DNAmTL.DNAmTL


# In[10]:


#to play with L0
#L0 = 8

plt.plot(t, avg_tel(t, *GS_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# #### Doing it the other way round just for continuity

# #### Fit to LBC

# In[11]:


#get L0 for dataset
L0 = find_L0(LBC_age_DNAmTL)

#inflate LBC a lot more by arbitrary bc dataset begins at such late ages 
L0 = L0 + 1.2


# In[12]:


LBC_pars, LBC_testcov = curve_fit(avg_tel, x_LBC, y_LBC)
print(LBC_pars)


# In[13]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# #### LBC fit to GS data

# In[14]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


