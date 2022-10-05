#!/usr/bin/env python
# coding: utf-8

# # Functions
# 
# Functions are pieces of code that we will want to use numerous times.  This could be as simple as $f(x) = sin(x)$, or something much more complicated.

# ## Learning goals
# 
# After this notebook, you should be able to
# 
# 1. Define a python function
# 2. Use a python function

# ## Syntax

# If we want to compute $f(x) = sin(x)$ numerous times throughout our code or our notebook, it can be useful to define a function that performs this manipulation and returns the result.  The syntax for this is 
# 
# `def function_name(arguments):`
# 
# Tabbed in below the definition, you write the lines of code that you want to perform.  Importantly, the last line should read 
# 
# `return some_value`
# 
# where your code needs to have defined variable `some_value` prior to the return line. This line returns this value just like $f(x)$ returns the value of $f$ at position $x$.
# 
# I will define a couple of example functions below

# In[1]:


import numpy as np
def function_sin(x):
    return np.sin(x)


# In[2]:


x = np.pi/4
print("sin(pi/4) = ", function_sin(x))


# Now this function is a bit silly for two reasons.  First is that the function is pretty simple so it is not really necessary to define a function to do this manipulation.  Second, our function actually just calls a `numpy` function for sin.  We didn't really write our own.  But that is ok for now.

# In[3]:


def polynomial_function(x):
    return x**5 + 4*x**2 - 2*x + 6


# In[4]:


x = 3
print("3^5 + 4*3^2 - 2*3 +6 = ", polynomial_function(x))


# ## Functions with more than one argument

# You can feed multiple arguments to fuction.  This could be because your function has more than one variable, e.g. $f(x,y)$, or because there are parameters/constants that need to be fed to the function at different times. 
# 
# We will use the ideal gas law as an example: $PV = nRT$.  This can be written as $P$ as a function of $n, V, T$
# \begin{equation}
# P(n,V,T) = \frac{nRT}{V}
# \end{equation}
# 
# but we see that we also have $R$, the universal gas constant to deal with.  The value of this constant will depend on the units that you are using for volume and the desired units of pressure.

# In[5]:


def p_ideal(n, V, T, R):
    return n*R*T/V


# In[6]:


R = 0.08206 # units of L atm mol^{-1} K^{-1}
n = 1       # units of mol
T = 300     # K
V = 22.4    # L
print("P(1mol, 22.4 L, 300K) = ", p_ideal(n,V,T,R), "atm")


# In[ ]:




