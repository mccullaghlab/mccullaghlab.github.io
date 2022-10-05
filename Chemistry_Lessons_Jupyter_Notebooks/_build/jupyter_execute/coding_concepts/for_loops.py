#!/usr/bin/env python
# coding: utf-8

# # for loops

# ## Learning goals

# ### For loops in python

# Loops are a common structure in coding/scripting.  

# In[1]:


for i in range(10):
    print(i)


# Loops are often used to compute sums.  For example if we want to know $\sum_{i=0}^{10} i$

# In[2]:


partialSum = 0
for i in range(11):
    partialSum += i
    print(i, partialSum)
print("Final sum is:", partialSum)


# Loops are often used to compute sums.  For example if we want to know $\sum_{i=1}^{10} i$

# In[3]:


partialSum = 0
for i in range(1,11):
    partialSum += i
    print(i, partialSum)
print("Final sum is:", partialSum)

