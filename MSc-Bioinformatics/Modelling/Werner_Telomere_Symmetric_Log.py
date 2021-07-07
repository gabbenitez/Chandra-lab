#!/usr/bin/env python
# coding: utf-8

# # Werner Telomere Shortening and Distributions

# https://elifesciences.org/articles/08687#s3

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy
from scipy.optimize import curve_fit
import math
import pandas as pd
import os


# In[4]:


## average telomere length of a cell shortens by a constant factor during each division
## this underestimates number of senescent cells once telomeres become critically short
## telomere shortening can only persist if they occur on stem cell population
## uses telomere lengths in mature cells as proxy for dist. of TL in SCs


# Asymmetric division - 1 HSc -> 1 HSc' + 1 progenitor (leaving Sc pool)
# 
#     1 HSc' keeps HSc properties, enters state i + 1, and telomeres shorten by delta.c
#     
# Symmetric division  - 1 HSc -> 1 HSc' + 1 HSc'
# 
#     Both keep HSc properties, enter state i + 1, telomeres shorten by delta.c

# ###### List of Parameters
# 
# 1 + c = accessible telomere states of stem cells
# 
# i = a state, with cells of equal average TL
# 
# NO = number of cells in state 0
# 
# p = probability of symmetric division
# 
# 1 - p = probability of asymmetric division
# 
# c = initial length
# 
# r = proliferation rate of cell

# In[18]:


# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_values = age_DNAmTL.age
y_values = age_DNAmTL.DNAmTL


# ## Expected Telomere Length

# #### Asymmetric Division
# Expected Telomere Length at time T = (c * NO - delta.c * r * t) / NO

# In[17]:


# we setup the equation
def exp_asy(t, c, NO, delta, r):
    y = (c*NO-delta*r*t)/NO
    return y


# In[45]:


#set up parameters (to play with)

t = np.linspace (18, 100, 1000)
c = 7.9
NO = 10
delta = 0.1
r = 1.2

#plotting scatter with model on top

plt.scatter(x_values, y_values, s = 2)
plt.plot(t, exp_asy(t, c, NO, delta , r), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

#plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# #### Symmetric division

# In[105]:


# we setup the equation
def exp_sym(t, c, delta, p, r, NO):
    
    y = c - delta * ((1+p)/p)*np.log(r*p/NO*t+1)
    
    return y


# In[149]:


#set up parameters (to play with)
t = np.linspace (18, 100, 1000)
c = 8
delta = 0.06
p = 0.08
r = 30
NO = 100


#plotting scatter with model on top
plt.scatter(x_values, y_values, s = 2)
plt.plot(t, exp_sym(t, c, delta, p, r, NO,), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

#plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# ## Curve Fitting

# In[ ]:


# to find initial CO specific to our data

#function to find LO, as an average of y values corresponding to head x values 
def find_LO(df):
    sorted_age_DNAmTL = df.sort_values('age')
    temp = sorted_age_DNAmTL.head(300)
    LO = temp['DNAmTL'].median()
    return LO


#get LO for dataset
LO = find_LO(age_DNAmTL)
LO


# #### Asymmetric Divison

# In[151]:


# optimise parameters
pars, testcov = curve_fit(exp_asy, x_values, y_values)
print(pars)
print(testcov)


# In[152]:


#plotting scatter with model on top
plt.scatter(x_values, y_values, s = 2)
plt.plot(t, exp_asy(t, *pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

#plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# #### Asymmetric Division

# In[154]:


# optimise parameters
pars, testcov = curve_fit(exp_sym, x_values, y_values)
print(pars)
print(testcov)


# In[155]:


#plotting scatter with model on top
plt.scatter(x_values, y_values, s = 2)
plt.plot(t, exp_sym(t, *pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

#plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# In[ ]:




