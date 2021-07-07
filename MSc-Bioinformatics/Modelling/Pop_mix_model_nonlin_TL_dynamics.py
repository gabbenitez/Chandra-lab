#!/usr/bin/env python
# coding: utf-8

# # Itzkovitz Population Mixture Model
# 
# #### http://shalevlab.weizmann.ac.il/wp-content/uploads/2016/02/telomeres.pdf

# In[261]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy
from scipy.optimize import curve_fit
import math
import pandas as pd
import os


# ## Average telomere length
# 

# In[208]:


# we setup the equation
def avg_tel(t, M, delta, alpha):
    y = LO - ((2*M*delta)/alpha)*(1-math.e**(-alpha*t))
    return y


# In[235]:


#parameters and set up
t = np.linspace (18, 100, 1000)

#parameters [initial length, division rate, fraction of attrition, death - division)
#LO is just from the GS, everything else has no scientific basis rn haha

LO = 8
M = 1000
delta = 0.00001
alpha = 0.01


# In[236]:


# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_values = age_DNAmTL.age
y_values = age_DNAmTL.DNAmTL


# In[237]:


# plotting scatter of GS data
# vs our model

plt.plot(t, avg_tel(t, M, delta, alpha,), color = 'red')
plt.scatter(x_values, y_values, s = 2)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# ## Curve Fitting

# #### scipy.optimize.curve_fit

# In[238]:


# Using LO based on our specific df
# and fitting
# M = rate of division
# delta = fraction of attrition
# alpha = rate of death - rate of division


# In[239]:


# function to find LO, as an average of y values corresponding to head x values 
def find_LO(df):
    sorted_age_DNAmTL = df.sort_values('age')
    temp = sorted_age_DNAmTL.head(300)
    LO = temp['DNAmTL'].median()
    return LO


#get LO for dataset
LO = find_LO(age_DNAmTL)
LO


# In[240]:


# optimise other parameters
pars, testcov = curve_fit(avg_tel, x_values, y_values)
print(pars)
print(testcov)


# In[241]:


r = len(x_values)

t = np.linspace (0, 100, r)

#this bit is just for playing around with 
#pars = [0.05917522, 0.06062484, -0.01415034]

plt.plot(t, avg_tel_LO(t, *pars), color = 'red')
plt.scatter(x_values, y_values, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# #### lmfit.Model

# In[242]:


# try it again with another function - Model, from lmfit
from lmfit import Model


# In[243]:


gmodel = Model(avg_tel_LO)
print('parameter names: {}'.format(gmodel.param_names))
print('independent variables: {}'.format(gmodel.independent_vars))


# In[244]:


#make params
params = gmodel.make_params(M = 1000, delta = 0.00001, alpha = 0.01)
params


# In[245]:


#evaluate 
x_eval = np.linspace (0, 100, 1000)
y_eval = gmodel.eval(params, x = x_eval)


# In[246]:


result = gmodel.fit(y_values, x=x_values, M = 1000, delta = 0.00001, alpha = 0.01) 
print(result.fit_report())


# In[260]:


r = len(x_values)

t = np.linspace (0, 100, r)

plt.scatter(x_values, y_values, s=10, alpha = 0.5)

plt.plot(x_values, result.init_fit, 'k--', label='initial fit')
plt.plot(x_values, result.best_fit, color = 'green', label='gmodel fit')
plt.plot(t, avg_tel(t, *pars), color = 'red', label = 'curve fit' )
plt.legend(loc='best')


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# ### Using Fit LO

# In[31]:


# here we also try and fit the LO value... which we get weird results with
def avg_tel_no_LO(t, LO, M, delta, alpha):
    y = LO - ((2*M*delta)/alpha)*(1-math.e**(-alpha*t))
    return y


# In[32]:


pars, testcov = curve_fit(avg_tel_no_LO, x_values, y_values)
print(pars)


# In[33]:


r = len(x_values)

t = np.linspace (0, 100, r)

plt.plot(t, avg_tel_no_LO(t, *pars), color = 'red')
plt.scatter(x_values, y_values, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((10, 100))

plt.show()


# In[ ]:


# mathematically, idgi rn

