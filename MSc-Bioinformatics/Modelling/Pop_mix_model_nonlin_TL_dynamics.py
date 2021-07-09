#!/usr/bin/env python
# coding: utf-8

# # Itzkovitz Population Mixture Model
# 
# #### http://shalevlab.weizmann.ac.il/wp-content/uploads/2016/02/telomeres.pdf

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy
from scipy.optimize import curve_fit
import math
import pandas as pd


# ## Average telomere length
# 

# In[2]:


#we set up the equation
def avg_tel(t, L0, M, delta, alpha):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# In[3]:


# set up time
t = np.linspace (0, 100, 1000)

# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_values_GS = age_DNAmTL.age
y_values_GS = age_DNAmTL.DNAmTL


# In[4]:


# for playing around with

#parameters
L0 = 8
M = 1000
delta = 0.00001
alpha = 0.01

# plotting scatter of GS data
# vs our model

plt.plot(t, avg_tel(t, L0, M, delta, alpha,), color = 'red')
plt.scatter(x_values_GS, y_values_GS, s = 2)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# ## Curve Fitting

# Infers parameters based on data

# In[5]:


x_GS=np.array(x_values_GS)
y_GS=np.array(y_values_GS)


# #### Fit to GS data

# In[6]:


#give it initial values
p_init = [8., 1000, 0.00001, 0.01]

GS_pars, GS_testcov = curve_fit(avg_tel, x_GS, y_GS, p0 = p_init)
print(GS_pars)


# In[7]:


plt.plot(t, avg_tel(t, *GS_pars), color = 'red')
plt.scatter(x_values_GS, y_values_GS, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# #### GS fit to LBC data

# See how model fitted to GS describes LBC

# In[8]:


df = pd.read_csv('../Data/sheets/LBC.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_LBC = age_DNAmTL.age
y_LBC = age_DNAmTL.DNAmTL


# In[9]:


plt.plot(t, avg_tel(t, *GS_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

#to play with L0 and shift it down to account for batch effects
#z = (avg_tel(t, *GS_pars)) - 0.35
#plt.plot(t, z, color = 'green')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# 
# ### Doing it the other way round just for continuity
# 
# Not expecting to see anything amazing given not much of an incline with LBC

# #### Fit to LBC

# In[10]:


#give it initial values also make runtime longer
p_init= [7.8, 1000, 0.00001, 0.01]


LBC_pars, LBC_testcov = curve_fit(avg_tel, x_LBC, y_LBC, p0 = p_init, maxfev=10000)
print(LBC_pars)


# In[11]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# #### LBC fit to GS data

# In[12]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5)

#to play with L0 and shift it down to account for batch effects
#z = (avg_tel(t, *LBC_pars)) - 0.35
#plt.plot(t, z, color = 'green')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[ ]:
