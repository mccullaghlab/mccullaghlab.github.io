#!/usr/bin/env python
# coding: utf-8

# # Arrays and Matrices

# ## Learning Goals
# 
# After this notebook, you will be able to:
# 
# 1. Generate an manipulate python lists 
# 2. Import numpy package
# 3. Use numpy arrays
# 4. Use numpy matrices

# ## Lists

# The python native version of an array is called a list.  A list contains multiple elements and can be of various data types.  There are multiple ways to initialize a list.  We will start by investigating how one initializes a list.

# ### Initializing lists

# In[1]:


# initialize list
a = [1,2,3]
print("a= ", a,type(a))
print("a[0] = ", a[0],type(a[0]))
print("a[1] = ", a[1],type(a[1]))
print("a[2] = ", a[2],type(a[2]))


# In the above code we created a list called `a` with three elements.  The first element is the integer value `1`, the second element is the integer value `2`, the third element is integer value `3`.  
# 
# Python is a C-based language so indexing of lists starts at 0.  So the first element the list is `a[0]`
# 
# We will now look at ways of adding to lists

# In[2]:


# initialize an empty list
a = []
print("a = ", a)
# append a value to a
a.append(3)
print("a = ", a)
# append another value to a
a.append(6)
print("a = ",  a)


# ### Data types within lists

# Lists can have various data types

# In[3]:


# initialize a list of varying data types
a = [1,2.0,3]
print(a,type(a))
print("a[0] = ", a[0],type(a[0]))
print("a[1] = ", a[1],type(a[1]))
print("a[2] = ", a[2],type(a[2]))


# In[4]:


a = [1,2.0,str(3)]
print(a,type(a))
print("a[0] = ", a[0],type(a[0]))
print("a[1] = ", a[1],type(a[1]))
print("a[2] = ", a[2],type(a[2]))


# ### Useful functions on lists

# There are various attributes of list objects and functions applied to lists that might be useful.

# In[5]:


a = [1,2,3]
print(a)
print("Length of a:", len(a))


# ## Numpy

# Lists can be useful but when we want to perform vector and matrix math they can be quite cumberson.  Python is a very common coding language and is actively developed in a variety of areas.  As such it has a number of ''libraries'' or bits of code to do common things that are already written.  If we import those libraries we can utilize these functions.  The most common example it the numpy library (http://www.numpy.org/).  This library contains a lot of functions to perform vector, matrix and other common mathematical manipulations. The first step is to 
# 
# `import numpy as np`. 

# In[6]:


import numpy as np


# ### Vectors

# In[7]:


# defining an array/vector
a = np.array([2,3,1])
print("a = ", a, type(a))


# In[8]:


# vectors have some attributes that might be of interest
print("Size of array a:", a.size)
print("Data type of array a:", a.dtype)


# For a list of all attributes of arrays see: https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.ndarray.html

# In[9]:


# vector manipulations
# for this I start by defining another vector in R3:
b = np.array([-1,0,7])
print("a=",a)
print("b=",b)
print("a+b=",a+b)
print("2a+3b=",2*a+3*b)
print("a*b=",a*b)
print("a*b=",np.dot(a,b))
# what are the differences between the last two lines?
print("a x b =",np.cross(a,b))


# ### Matrices

# We can define a matrix in a number of ways.  The first is to cast a list of lists as a numpy matrix.

# In[10]:


# defining a matrix - these are just 2D arrays
a = np.matrix([[1,2],[3,4]])
b = np.matrix([[1,2],[3,4]],dtype=float)
c = np.matrix([[5,6],[-1,0]],dtype=float)
print(a)
print(b)
print(c)


# We can also declare the matrices and populate them element by element.

# In[11]:


a = np.empty((3,3),dtype=float)   # declare an empty matrix of size 3x3 and type float
for i in range(3):
    for j in range(3):
        a[i,j] = i+j              # place the value of i+j into the i,jth element of the matrix
print(a)


# Elements in a matrix can be accessed by square brackets and the integer indeces

# In[12]:


print("a[0,0]=",a[0,0])
print("a[1,2]=",a[1,2])


# ### Manipulating Matrices

# Numpy has a lot of built in functions to manipulate or assess a matrix.  

# In[13]:


# trace
print("Tr(a)=",np.trace(a))
# determinant
print("|a|=det(a)=",np.linalg.det(a))
# transpose
print("a^T=",a.T)


# More complicated matrix manipulations such as inverses and diagonalizations can also be done

# In[14]:


# inverse
np.linalg.inv(a)


# <div class="alert alert-success">
# Note that not all matrices are invertible!
# </div>

# In[8]:


# eigenvalues and eigenvectors
e,v = np.linalg.eig(a)
print("Eigenvalues:",e)
print("Eigenvectors:",v)

