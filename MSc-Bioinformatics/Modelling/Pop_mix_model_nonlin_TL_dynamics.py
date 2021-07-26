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

# In[2]:#!/usr/bin/env python
# coding: utf-8

# # Itzkovitz Population Mixture Model
# 
# #### http://shalevlab.weizmann.ac.il/wp-content/uploads/2016/02/telomeres.pdf

# In[3]:


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

# In[4]:


#we set up the equation
def avg_tel(t, L0, M, delta, alpha):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# In[5]:


# set up time
t = np.linspace (0, 100, 1000)

# get GS data and extract what we need
df = pd.read_csv('../Data/sheets/Generation_Scotland.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_values_GS = age_DNAmTL.age
y_values_GS = age_DNAmTL.DNAmTL


# In[6]:


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

# In[7]:


x_GS=np.array(x_values_GS)
y_GS=np.array(y_values_GS)


# #### Fit to GS data

# In[8]:


#give it initial values
p_init = [8, 1000, 0.00001, 0.01]

GS_pars, GS_testcov = curve_fit(avg_tel, x_GS, y_GS, p0 = p_init)
print(GS_pars)

#calculate variance
GS_variance = np.sqrt(np.diagonal(GS_testcov))

#this is one standard deviation 
print(GS_variance)


# In[9]:


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


# In[10]:


#to play with

GS_pars1, GS_testcov1 = curve_fit(avg_tel, x_GS, y_GS, p0 = p_init)
print(GS_pars1)

#calculate variance
GS_variance1 = np.sqrt(np.diagonal(GS_testcov1))

#this is one standard deviation 
print(GS_variance1)


# #### GS fit to LBC data

# See how model fitted to GS describes LBC

# In[11]:


df = pd.read_csv('../Data/sheets/LBC.csv', index_col=0)

age_DNAmTL = df[["age", "DNAmTL"]]
x_LBC = age_DNAmTL.age
y_LBC = age_DNAmTL.DNAmTL


# In[12]:


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

# In[13]:


#give it initial values also make runtime longer
p_init= [7.8, 1000, 0.00001, 0.01]

#bounds = [[8, 900, 0, 0], [9, 2000, 1, 1]] 

LBC_pars, LBC_testcov = curve_fit(avg_tel, x_LBC, y_LBC, p0 = p_init, maxfev=100000)
print(LBC_pars)


# In[161]:


np.sqrt(np.diagonal(LBC_testcov))


# In[14]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 15.0))
plt.xlim((0, 100))

plt.show()


# #### LBC fit to GS data

# In[15]:


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

# In[16]:


#### Fit to LBC

#give it initial values also make runtime longer
p_init= [7.8, 1000, 0.00001, 0.01]

bounds = [[8, 900, 0, 0], [12, 2000, 1, 1]], 

LBC_pars, LBC_testcov = curve_fit(avg_tel, x_LBC, y_LBC, p0 = p_init, bounds = [[7, 900, 0, 0], [12, 2000, 1, 1]], maxfev=10000)
print(LBC_pars)


# In[17]:


np.sqrt(np.diag(LBC_testcov))


# In[18]:


plt.plot(t, avg_tel(t, *LBC_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

#plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[19]:


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

# In[20]:


#testing things
#get the L0 from GS
L0 = 8.28581491


# In[21]:


#we set up the new equation
def tel(t, M, delta, alpha):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# ### Fit to all LBC Data (INCORRECT!)
# 
# this is actually NOT good, it infers the wrong parameters based on the wrong offset. Keeping this here to humble myself.
# 
# What should be done, and is below this, is fitting to the LBC values, but then changing the offset to LBC after.

# In[22]:


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
# So it looks like the issue is that the previous LBC fit massively overestimates the L0, so what if we just use GS L0 value.
# 
# Infers LBC dynamics from its data, then changes the L0

# In[23]:


# we use the same LBC params saved from earlier
LBC_pars_GSL0 = LBC_pars[1:]

#and insert the GS L0 value
LBC_pars_GSL0 = np.insert(LBC_pars_GSL0, 0, L0)
print (L0)


# In[24]:


plt.plot(t, avg_tel(t, *LBC_pars_GSL0), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


#  ---

#  ---

# # Start Reading Here

# # Fitting to Individual LBC Trajectories
# 
# Things to try out
# - Using GS parameters, fit the L0 for LBC. (This assumes similar dynamics between populations in the datasets)

# ## Using GS dynamics to find L0 fit for LBC
# 
# If telomere dynamics between individuals in the two datasets behave roughly the same, can we use GS parameters to infer the L0 of LBC (or at least get it in biologically relevant bounds)?
# 
# What is preventing us from getting good fits with LBC: lack of data at earlier times, and different dynamics in later ages?
# 
# Do: 
# - Get GS params
# - keep all params except GS L0 (we assume that dynamics between LBC and GS are similar, and offset bc of batch effects)
# - get LBC L0 fit (so given similar dynamics, what is the average L0 in LBC)
# - Keep the new LBC L0 set, and infer the other parameters
# - New LBC L0 vs GS L0 is the offset

# In[25]:


#we set up the new equation to get the fit of L0 to
def tel_l0(t, L0):
    y = L0 - ((2*M*delta)/alpha) * (1-np.exp(-alpha*t))
    return y


# In[26]:


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


#i dont know how mathematically true this could be
offset = test - LBC_l0
print("the offset is", offset)


# In[162]:


np.sqrt(np.diag(LBC_lo_cov))


# In[27]:


plt.plot(t, tel_l0(t, LBC_l0), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)
#plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[28]:


#now using this inferred L0 for LBC, we fit the other parameters
L0 = LBC_l0

# we use the function tel, with params are M, delta, alpha
# and help it out with the init

p_init = [1047.7183331711515, 1.0437599240704384e-05, 0.006803590426130318]

LBC_l0_pars, LBC_l0_cov = curve_fit(tel, x_LBC, y_LBC, p0 = p_init,  maxfev = 10000)
LBC_l0_pars


# In[163]:


np.sqrt(np.diag(LBC_l0_cov))


# In[29]:


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

# In[30]:


#L0 = 8.28581381

plt.plot(t, tel(t, *LBC_l0_pars), color = 'red')
plt.scatter(x_LBC, y_LBC, s=10, alpha = 0.5)
plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

plt.ylabel('DNAmTL')
plt.yscale('linear')
plt.xlabel('Age')

plt.ylim((5.5, 9.0))
plt.xlim((0, 100))

plt.show()


# In[31]:


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


# In[32]:


print (LBC_l0_pars)
print (GS_pars)


# #### Get variance of L0

# In[33]:


print(LBC_l0)
print(LBC_lo_cov)
#calculate one standard deviation of the parameters
perr = np.sqrt(np.diag(LBC_lo_cov))
print("variance of L0", perr)
np.sqrt(np.diag(LBC_l0_cov))


# # Fit to individual trajectories

# In[34]:


#get the data in the correct format (i.e. per Sample_Name)
#so an individuals DNAmTL measure at each wave with their age


#maybe loop to go through list of unique LBC Sample_Names
#each time make a subset of age and DNAmTL for each person 
#fit parameters to the persons trajectory
#add parameters to empty array, new row each time

#then, plot (with good alpha vibes obvs) inferred full function...


# In[145]:


#get data
# get LBC data and extract what we need
df = pd.read_csv('../Data/sheets/LBC_with_LTL_update.csv', index_col=0)
df = df[["Sample_Name", "age", "DNAmTL"]]

#get list of all sample names
Sample_Names = df['Sample_Name'].unique()

#x_LBC = age_DNAmTL.age
#y_LBC = age_DNAmTL.DNAmTL


# In[146]:


Sample_Names = df['Sample_Name'].unique()


# In[147]:


def makepars(g,):
    for i in g:

        subset = df.loc[df['Sample_Name'] == i]
        x_sub = subset.age
        y_sub = subset.DNAmTL
    
    
        #fit params to trajectory of subset with inferred LBC_lo
        #but give bounds for L0 within 1 SDEV of inferred
        #7.93586611
        p_init = [7.93586611, 1000, 0.00001, 0.01]
        #bounds = [[7.93205815, 900, 0, 0], [7.93967407, 2000, 1, 1]]

        #we give L0 bounds from 7-9
        sub_pars, sub_cov = curve_fit(avg_tel, x_sub, y_sub, p0 = p_init, bounds = [[7.5, 900, 0, 0], [13, 2000, 1, 1]], maxfev = 1000000)
    
        #calculate variance
        sub_var = np.sqrt(np.diagonal(sub_cov))

        #add to blank arrays for parameters and variance 
        blankpars.append(sub_pars)
        blankvars.append(sub_var)
        
        print(len(blankpars))

    return


# In[148]:


blankpars = []
blankvars = []
makepars(Sample_Names)


# In[151]:


LBC_pars = blankpars 
LBC_vars = blankvars


# In[40]:


for a in blankpars:
    
    plt.plot(t, avg_tel(t, *a), linewidth = 0.5, alpha = 0.1)
    #plt.scatter(x_GS, y_GS, s=10, alpha = 0.5, color = 'grey')

    plt.ylabel('DNAmTL (kB)')
    plt.yscale('linear')
    plt.xlabel('Age')

    plt.ylim((0, 14))
    plt.xlim((0, 100))

print("done")
plt.scatter(x_LBC, y_LBC, s=1, alpha = 0.5, color = 'black')
plt.show()


# There's 3 "groups" of different dynamics here, one which goees to the p_init start value of 8, one which goes to the max of 9, and others which are more spread out

# In[41]:


#what is the collection at the minimum value???? is this because there is no spread???
#we try and do this again for LBC_21 and LBC_36 separately...


# In[113]:


#get LBC 21 LTL data and extract what we need
df = pd.read_csv('../Data/sheets/LBC_21.csv', index_col=0)

temp = df[["age", "DNAmTL", "Sample_Name"]]
x_21 = temp.age
y_21 = temp.DNAmTL

Sample_Names_21 = df['Sample_Name'].unique()


# In[114]:





# In[115]:


blankpars = []
blankvars = []
makepars(Sample_Names_21)


# In[116]:


LBC_21_pars = blankpars
LBC_21_vars = blankvars


# In[164]:


for a in blankpars_21:
    
    plt.plot(t, avg_tel(t, *a), alpha = 0.1)
    #plt.scatter(x_GS, y_GS, s=10, alpha = 0.1, color = 'grey')

    plt.ylabel('DNAmTL')
    plt.yscale('linear')
    plt.xlabel('Age')

    plt.ylim((0, 14.0))
    plt.xlim((0, 100))

print("done")
plt.scatter(x_21, y_21, s=10, alpha = 0.8, color = 'black')
plt.show()


# In[117]:


#get LBC 36 LTL data and extract what we need
df = pd.read_csv('../Data/sheets/LBC_36.csv', index_col=0)

temp = df[["age", "DNAmTL"]]
x_36 = temp.age
y_36 = temp.DNAmTL

Sample_Names_36 = df['Sample_Name'].unique()


# In[118]:


blankpars = []
blankvars = []
makepars(Sample_Names_36)


# In[119]:


LBC_36_pars = blankpars
LBC_36_vars = blankvars


# In[52]:


for a in blankpars_36:
    
    plt.plot(t, avg_tel(t, *a), alpha = 0.1)
    plt.scatter(x_GS, y_GS, s=10, alpha = 0.1, color = 'grey')

    plt.ylabel('DNAmTL')
    plt.yscale('linear')
    plt.xlabel('Age')

    plt.ylim((0, 14))
    plt.xlim((0, 100))

print("done")
plt.scatter(x_36, y_36, s=10, alpha = 0.8, color = 'black')
plt.show()


# In[155]:


len(LBC_pars)


# In[72]:


#get variance


# In[156]:


pars_36 = pd.DataFrame(list(map(np.ravel, LBC_36_pars)))
pars_21 = pd.DataFrame(list(map(np.ravel, LBC_21_pars)))
pars_LBC = pd.DataFrame(list(map(np.ravel, LBC_pars)))


# In[159]:


pars_LBC


# In[160]:


print(pars_LBC.var())

print(pars_36.var())

print(pars_21.var())


# In[165]:


pars_36


# In[150]:





# In[ ]:







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
