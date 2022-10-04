#!/usr/bin/env python
# coding: utf-8

# # Plotting with matplotlib

# ## Learning Goals
# 
# After going through this notebook, you should be able to:
# 
# 1. import matplotlib in a Jupyter Notebook 
# 2. Make a basic plot of a 1D function with matplotlib
# 3. Change the size of a maplotlib plot
# 4. Change font sizes in matplotlib
# 5. Change tickmarks
# 6. Add a legend
# 7. Add labels to axes

# ## Importing libraries

# matplotlib is a python library used for making plots.  It is not native to python and thus we must tell our python kernel to import it.  In addition, and something specific to Jupyter Notebooks, we will have to tell the notebook that we want to display the plot inline

# In[1]:


import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# while we are at it we will import numpy

# In[2]:


import numpy as np


# ## Basic plotting

# Let's plot $f(x) = x^2$.  To have a computer do this, we must explicitly define our $x$ values and compute the corresponding $f(x)$ value. 

# In[3]:


# create a list of x values
xValues = []   # initialize list
for i in range(100):
    xValues.append(-5 + i*0.1)
print(xValues)


# In[4]:


# determine f(x) values
yValues = []   # initialize list
for i in range(100):
    yValues.append(xValues[i]*xValues[i])   # f(x) = x^2 = x*x
print(yValues)


# In[5]:


# plot
plt.plot(xValues,yValues)


# The above lists and plots could also be all done in one code snippet.  I will do so below and also use numpy routines and operations to make the code look cleaner.

# In[6]:


x_vals = np.arange(-5,5,0.1)
y_vals = x_vals**2
plt.plot(x_vals,y_vals)


# ## Multiple curves on the same plot

# In[7]:


x_vals = np.arange(-5,5,0.1)
y_vals_1 = x_vals**2
y_vals_2 = x_vals**3
plt.plot(x_vals,y_vals_1)
plt.plot(x_vals,y_vals_2)


# ## Changing linewidth

# In[8]:


x_vals = np.arange(-5,5,0.1)
y_vals_1 = x_vals**2
y_vals_2 = x_vals**3
plt.plot(x_vals,y_vals_1,lw=2)
plt.plot(x_vals,y_vals_2,lw=5)


# ## Changing the size of the plot

# What I mean by size, here, is the image size.  If you want to plot a larger domain/range you must simply change the number of x and y values you compute.
# 
# There are likely an almost infinite number of ways of changing the image size of a matplotlib plot.  Here we will use the approach of declaring a figure object and then modifying that object.

# In[9]:


# declare a figure object and set the figure size
fig = plt.figure(figsize=(6,6),dpi=80)
plt.plot(x_vals,y_vals)


# ## Changing the font size in matplotlib plots

# Changing the fontsizes in matplotlib can be tricky because one can change the font size of different parts of the plot in different spots.  For the current plot, we only have text for along our axes so we will modify those.  We will, however, modify the default fontsize which will affect our subsequent plots as well.

# In[10]:


# update default font size to 16
plt.rcParams.update({'font.size': 16})


# In[11]:


# declare a figure object and set the figure size
fig = plt.figure(figsize=(6,6),dpi=80)
plt.plot(x_vals,y_vals)


# In[ ]:




