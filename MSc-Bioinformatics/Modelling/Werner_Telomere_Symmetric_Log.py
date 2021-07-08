#!/usr/bin/env python
# coding: utf-8

# # Werner Telomere Shortening and Distributions

# https://elifesciences.org/articles/08687#s3

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy
from scipy.optimize import curve_fit
import math
import pandas as pd
import os


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

# In[2]:


# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_GS = age_DNAmTL.age
y_GS = age_DNAmTL.DNAmTL


# ## Expected Telomere Length

# #### Asymmetric Division
# Expected Telomere Length at time T = (c * NO - delta.c * r * t) / NO

# In[3]:


# we setup the equation
def exp_asy(t, c, NO, delta, r):
    y = (c*NO-delta*r*t)/NO
    return y


# In[4]:


#set up parameters (to play with)

t = np.linspace (18, 100, 1000)
c = 7.9
NO = 10
delta = 0.1
r = 1.2

#plotting scatter with model on top

plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_asy(t, c, NO, delta , r), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((18, 100))

plt.show()


# #### Symmetric division

# In[5]:


# we setup the equation
def exp_sym(t, c, delta, p, r, NO):
    
    y = c - delta * ((1+p)/p)*np.log(r*p/NO*t+1)
    
    return y


# In[6]:


#set up parameters (to play with)
t = np.linspace (18, 100, 1000)
c = 8
delta = 0.06
p = 0.08
r = 30
NO = 100


#plotting scatter with model on top
plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_sym(t, c, delta, p, r, NO,), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((18, 100))

plt.show()


# ## Curve Fitting

# #### Asymmetric Divison

# In[7]:


# optimise parameters
pars, testcov = curve_fit(exp_asy, x_GS, y_GS)
print(pars)
print(testcov)


# In[8]:


#plotting scatter with model on top
plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_asy(t, *pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((18, 100))

plt.show()


# ---

# # Symmetric Division - Fitting Parameters

# Using the _symmetric_ model (Model 2 in the paper, the one w/ fit better), fit parameters using GS

# In[9]:


# optimise parameters
GS_pars, GS_testcov = curve_fit(exp_sym, x_GS, y_GS)
print(GS_pars)
print(GS_testcov)


# In[10]:


#plotting scatter with model on top
plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_sym(t, *GS_pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((18, 100))

plt.show()


# With the new parameters, fit to GS, we see how the trajectory looks with LBC data

# In[11]:


#get LBC data and extract what we need
df = pd.read_csv('../Data/sheets/LBC.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_LBC = age_DNAmTL.age
y_LBC = age_DNAmTL.DNAmTL


# In[12]:


#plotting scatter with model on top
plt.scatter(x_LBC, y_LBC, s = 2)
plt.plot(t, exp_sym(t, *GS_pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((18, 100))

plt.show()


# To tick all the boxes, we do it the other way i.e. fit our parameters to LBC and see how it maps to GS

# In[13]:


pars, testcov = curve_fit(exp_sym, x_LBC, y_LBC)
print(pars)
print(testcov)


# In[14]:


#plotting scatter with model on top
plt.scatter(x_LBC, y_LBC, s = 2)
plt.plot(t, exp_sym(t, *pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((18, 100))

plt.show()


# In[15]:


#plotting scatter with model on top
plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_sym(t, *pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((18, 100))

plt.show()
