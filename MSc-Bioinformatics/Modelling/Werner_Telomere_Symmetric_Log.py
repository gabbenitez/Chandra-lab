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


# ###### List of Parameters
# 
# c = telomere length
# 
# delta = rate of telomere attrition
# 
# p = probabikity of symmetric division 
# 
# r = proliferation rate of cell
#  
# N0 = number of cells in state 0 

# In[2]:


# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_GS = age_DNAmTL.age
y_GS = age_DNAmTL.DNAmTL


# ## Expected Telomere Length

# #### Asymmetric Division

# In[3]:


# we setup the equation
def exp_asy(t, c, NO, delta, r):
    y = (c*NO-delta*c*r*t)/NO
    return y


# In[4]:


#set up parameters (to play with)

t = np.linspace (0, 100, 1000)

c = 8
NO = 2
delta = 0.06
r = 0.06

#plotting scatter with model on top

plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_asy(t, c, NO, delta , r), color = 'red')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# #### Symmetric division

# Similar to the paper, this fits our data better and also has the flattening out at later ages

# In[5]:


# we setup the equation
def exp_sym(t, c, delta, p, r, NO):
    y = c - delta * c * ((1+p)/p)*np.log(r*p/NO*t+1)
    return y


# In[6]:


#set up parameters (to play with)
c = 8.5
delta = 0.05
p = 0.02
r = 30
NO = 500

#plotting scatter with model on top
plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_sym(t, c, delta, p, r, NO,), color = 'red')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# ## Curve Fitting

# #### Asymmetric Divison
# 
# Just done for comparisons sake

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
plt.xlim((0, 100))

plt.show()


# ---

# # Symmetric Division - Fitting Parameters

# Using the _symmetric_ model (Model 2 in the paper, the one w/ fit better), fit parameters using GS
# 
# Interesting the curve is concave, slow down at later ages, whcih is opposite to human survivorship curve

# #### Fit to GS Data

# - c = telomere length
# 
# - delta = rate of telomere attrition
# 
# - p = probabikity of symmetric division 
# 
# - r = proliferation rate of cell
#  
# - N0 = number of cells in state 0 

# In[9]:


# optimise parameters, and give realistic bounds

p_init = [8, 0.05, 0.01, 30, 100]

bounds=[[7,0, 0, 0, 0],[10, 1, 1, 1000, 1000]]

GS_pars, GS_testcov = curve_fit(exp_sym, x_GS, y_GS, p0 = p_init, bounds = bounds)
print(GS_pars)


# In[10]:


#plotting scatter with model on top
plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_sym(t, *GS_pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# #### GS fit to LBC data

# looks quite nice! RE batch effects, can just shift the values down a bit but not sure how scientifically valid this is

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

#to play with
#z = (exp_sym(t, *GS_pars))-0.3
#plt.plot(t, z, color = 'green')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# ### Fit to LBC

# To tick all the boxes, we do it the other way i.e. fit our parameters to LBC and see how it maps to GS

# In[13]:


# optimise parameters, and give realistic bounds

p_init = [8, 0.05, 0.01, 30, 100]

bounds=[[7,0, 0, 0, 0],[10, 1, 1, 1000, 1000]]

pars, testcov = curve_fit(exp_sym, x_LBC, y_LBC, bounds = bounds)
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
plt.xlim((0, 100))

plt.show()


# #### LBC fit to GS

# In[15]:


#plotting scatter with model on top
plt.scatter(x_GS, y_GS, s = 2)
plt.plot(t, exp_sym(t, *pars), color = 'red')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[ ]:

