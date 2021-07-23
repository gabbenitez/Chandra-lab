#!/usr/bin/env python
# coding: utf-8

# # Itzkovitz Population Mixture Model
# 
# #### http://shalevlab.weizmann.ac.il/wp-content/uploads/2016/02/telomeres.pdf

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import scipy
from scipy.optimize import curve_fit
import math
import pandas as pd


# Exponential decrease in average TL with time
# 
# Nice model, has a flattening at later ages but assumes TL in SC population are of constant length - biggest case of bias in data from individuals with abnormal haematopoiesis
# 
# 
# Validated model against longitudinal leukocyte TL measurements in baboons at different ages, fits better than a linear delcine in TL with age
# 
# 
# Linear TL dynamics at single-cell level give rise to exponential dynamics at pop. level, from repop. pool of SC w/ constant length, and derived pool of cells with TL decreasing linearly w/ each division

# ---

# - Two pools of cells, one of repopulating SC and a derviced pool of progenitor/differentiated cells
# 
# 
# - TL in SC remains constant
# 
# 
# - TL in derived cells reduces by constant fraction $\Delta$
# 
# 
# - Derived cell has constant _M_ rate of dividing
# 
# 
# - Derived cell has constant _D_ rate of dying (independent of TL length)
# 
# 
# - _M_ in SC << _M_ in derived cells
# 
# 
# - _D_ - _M_ = $\alpha$ '
# 
# 
# - _L_<sub>0</sub> as the intial TL, identical for repopulating cells

# ##### Average Telomere Length equation

# L(t) = _L_<sub>0</sub> - $\frac{2 M \Delta}{\alpha}$ (1 - $e^{ \alpha' t})$.

# ## Average telomere length
# 

# In[3]:


#we set up the equation
def avg_tel(t, L0, M, delta, alpha):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# In[74]:


# set up time
t = np.linspace (0, 100, 1000)

# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_values_GS = age_DNAmTL.age
y_values_GS = age_DNAmTL.DNAmTL


# In[75]:


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

# In[6]:


x_GS=np.array(x_values_GS)
y_GS=np.array(y_values_GS)


# #### Fit to GS data

# In[55]:


#give it initial values
p_init = [8, 1000, 0.00001, 0.01]

GS_pars, GS_testcov = curve_fit(avg_tel, x_GS, y_GS, p0 = p_init)
print(GS_pars)

#calculate variance
GS_variance = np.sqrt(np.diagonal(GS_testcov))

#this is one standard deviation 
print(GS_variance)

#?????
GS_varup = GS_pars + GS_variance
GS_vardown = GS_pars - GS_variance


# In[58]:


from uncertainties import ufloat
a = ufloat(GS_pars[0], GS_variance[0])
b = ufloat(GS_pars[1], GS_variance[1])
c = ufloat(GS_pars[2], GS_variance[2])
d = ufloat(GS_pars[3], GS_variance[3])

text_res = "Best fit parameters:\na = {}\nb = {}\nc = {}\nd = {}".format(a, b, c, d)
print(text_res)


# In[68]:


plt.plot(t, avg_tel(t, *GS_pars), color = 'red')
plt.scatter(x_values_GS, y_values_GS, s=10, alpha = 0.5)

#plt.plot(t, avg_tel(t, *GS_varup), color = 'black')
#plt.plot(t, avg_tel(t, *GS_vardown), color = 'black')


#bound_upper = avg_tel(t, *(GS_pars+ GS_variance))
#bound_lower = avg_tel(t, *(GS_pars - GS_variance))
# plotting the confidence intervals
#plt.fill_between(t, bound_lower, bound_upper,
                 #color = 'black', alpha = 0.15)


plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

#plt.ylim((5.5, 9.0))
#plt.xlim((0, 100))

plt.show()


# #### GS fit to LBC data

# See how model fitted to GS describes LBC

# In[9]:


df = pd.read_csv('../Data/sheets/LBC.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_LBC = age_DNAmTL.age
y_LBC = age_DNAmTL.DNAmTL


# In[10]:


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

# In[11]:


#give it initial values also make runtime longer
p_init= [7.8, 1000, 0.00001, 0.01]

# bounds = [[7, 900, 0, 0], [9, 2000, 1, 1]], 

LBC_pars, LBC_testcov = curve_fit(avg_tel, x_LBC, y_LBC, p0 = p_init, maxfev=10000)
print(LBC_pars)


# In[12]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# #### LBC fit to GS data

# In[13]:


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


# ### Playing with bounds here

# In[14]:


#### Fit to LBC

#give it initial values also make runtime longer
p_init= [7.8, 1000, 0.00001, 0.01]

bounds = [[7, 900, 0, 0], [9, 2000, 1, 1]], 

LBC_pars, LBC_testcov = curve_fit(avg_tel, x_LBC, y_LBC, p0 = p_init, bounds = [[7, 900, 0, 0], [9, 2000, 1, 1]], maxfev=10000)
print(LBC_pars)


# In[15]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[16]:


#### LBC fit to GS data

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


#  ---

# # Fit to Individual LBC Trajectories
# 
# As a test, we're going to keep the initial telomere length static, so keeping the same offset. 
# Then:
# - fit to all LBC data
# - fit to individuals data (will take a while!)
# - get spread of other parameters
# - plot the functions using this full range of parameters to GS
# 
# #### L0 will be = 8.28
# 
#  

# In[17]:


L0 = 8.28581491


# In[18]:


#we set up the new equation
def tel(t, M, delta, alpha):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# ### Fit to all LBC Data (INCORRECT!)
# 
# this is actually NOT good, it infers the wrong parameters based on the wrong offset. Keeping this here to humble myself.
# 
# What should be done, and is below this, is fitting to the LBC values, but then changing the offset to LBC after.

# In[19]:


p_init = [1000, 0.00001, 0.01]

LBC_lo_pars, LBC_lo_testcov = curve_fit(tel, x_LBC, y_LBC, p0 = p_init, maxfev = 100000)

print(LBC_lo_pars)

plt.plot(t, tel(t, *LBC_lo_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# ### Fit to LBC, then change to GS L0
# 
# So it looks like the issue is that the previous LBC fit massively overestimates the L0, meaning that the rate of attrition is incredibly low...

# In[20]:


# we use the same LBC params saved from earlier
LBC_pars_noL0 = LBC_pars[1:]
LBC_pars_noL0 = np.insert(LBC_pars_noL0, 0, L0)


# In[21]:


plt.plot(t, avg_tel(t, *LBC_pars_noL0), color = 'red')
#plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)
#plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

#plt.ylim((5.5, 9.0))
#plt.xlim((0, 100))

plt.show()


# ## Some cool stuff 
# 
# Try fitting with specific parameter bounds of GS using LBC data...
# 
# Is it because the dynamics at later ages are so different that we can't infer good fits???
# 
# So - use GS params to infer LBC L0...
# 
# Do: 
# - fit to GS (cross sectional over much time gives inference to the average population dynamics)
# - keep all params except GS L0 (we assume that dynamics between LBC and GS are similar, the offset however is not bc of batch effects)
# - get LBC L0 fit (so given similar dynamics, what is the average L0)
# - From this L0 fit, keep it static, and infer the other parameters...
# - New LBC L0 vs GS L0 is the offset
# 
# Assumes dynamics/ratio of parameters are similar between the datasets....

# In[22]:


#we set up the new equation to get the fit of L0 to
def tel_l0(t, L0):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# In[23]:


#Get GS fits and keep all params except GS L0
GS_pars_nol0 = GS_pars[1:]

#set up the GS parameters for the function
M = 1047.7183331711515
delta = 1.0437599240704384e-05
alpha = 0.006803590426130318

#sanity test, we get the right value for GS
test, test_1 =curve_fit(tel_l0, x_GS, y_GS, p0 = 8.28, maxfev = 10000)
print(test)

#now we try to get the L0 for LBC
#keep all params (dynamics) and infer offset
LBC_l0, LBC_lo_cov =curve_fit(tel_l0, x_LBC, y_LBC, p0 = 8.28, maxfev = 10000)
print(LBC_l0)

offset = test - LBC_l0
print("the offset is", offset)


# In[24]:


#now using this inferred L0 for LBC, we fit the other parameters
L0 = LBC_l0

# we use the function tel, with params are M, delta, alpha
# and help it out with the init

p_init = [1047.7183331711515, 1.0437599240704384e-05, 0.006803590426130318]

LBC_l0_pars, LBC_l0_cov = curve_fit(tel, x_LBC, y_LBC, p0 = p_init,  maxfev = 10000)
LBC_l0_pars


# In[25]:


plt.plot(t, tel(t, *LBC_l0_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)
#plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# Now how do the LBC cell pars fit with the GS l0?

# In[26]:


L0 = 8.28581381

plt.plot(t, tel(t, *LBC_l0_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[27]:


L0 = 8.28581381

#just for bants, adding the offset to LBC...

plt.plot(t, tel(t, *LBC_l0_pars), color = 'red')
plt.scatter(x_LBC, y_LBC + offset, s=10, alpha = 0.5)
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[28]:


offset


# #### Get variance of L0
# 
# Can use variance to look at bounds of parameters

# In[69]:


print(LBC_l0)
print(LBC_lo_cov)


# In[70]:


#calculate one standard deviation of the parameters
perr = np.sqrt(np.diag(LBC_lo_cov))
print("variance of L0", perr)


# In[71]:


np.sqrt(np.diag(LBC_l0_cov))


# # Fit to individual trajectories

# In[72]:


#get the data in the correct format (i.e. per Sample_Name)
#so an individuals DNAmTL measure at each wave with their age


#maybe loop to go through list of unique LBC Sample_Names
#each time make a subset of age and DNAmTL for each person 
#fit parameters to the persons trajectory
#add parameters to empty array, new row each time

#then, plot (with good alpha vibes obvs) inferred full function...


# In[76]:


#get data
# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/LBC_with_LTL_update.csv', index_col=0)
df = df[["Sample_Name", "age", "DNAmTL"]]

#get list of all sample names
Sample_Names = df['Sample_Name'].unique()

#x_LBC = age_DNAmTL.age
#y_LBC = age_DNAmTL.DNAmTL


# In[111]:


Sample_Names = df['Sample_Name'].unique()

LBC_l0


# In[ ]:


blankpars = []
blankvars = []

for i in Sample_Names:

    subset = df.loc[df['Sample_Name'] == i]
    x_sub = df.age
    y_sub = df.DNAmTL
    
    
    #fit params to trajectory of subset with inferred LBC_lo
    #but give bounds for L0 within 1 SDEV of inferred

    p_init = [7.93586611, 1000, 0.00001, 0.01]
    bounds = [[7.93205815, 900, 0, 0], [7.93967407, 2000, 1, 1]]

    sub_pars, sub_cov = curve_fit(avg_tel, x_sub, y_sub, p0 = p_init, bounds = [[7.93205815, 900, 0, 0], [7.93967407, 2000, 1, 1]], maxfev = 10000)
    
    #calculate variance
    sub_var = np.sqrt(np.diagonal(sub_cov))

    #add to blank arrays for parameters and variance 
    blankpars.append(sub_pars)
    blankvars.append(sub_var)

    plt.plot(t, avg_tel(t, *sub_pars), alpha = 0.2)
    plt.scatter(x_LBC, y_LBC + offset, s=10, alpha = 0.5)
    #plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

    plt.ylabel('DNAmTL')
    plt.yscale('linear')
    plt.xlabel('Age')

    plt.ylim((5.5, 9.0))
    plt.xlim((0, 100))

    plt.show()


# In[ ]:





# In[ ]:



