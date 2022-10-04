#!/usr/bin/env python
# coding: utf-8

# # Jupyter Notebook Tutorial

# Jupyter notebook is a collection of cells.  These cells can be one of three main types:
# 1. Code
# 
#     This includes executable code (usually python).  You can have different snippets of code in different cells.
#     
# 2. Markdown
# 
#     This is text with some formatting including numbered lists etc.  The current cell is an example.  You can google the text formatting syntax if interested.  Markdown also allows you to do headers by starting the cell with # (main header), ## (sub header), or ### (sub-sub header).  I use these in this and other notebooks.
#        
# 3. Raw
# 
#     Literally prints what you type.  Looks like code.  Is often used to convey code without it being executable.
#     
#     
# The cell type must be selected from dropdown menu in toolbar above.

# You can edit a cell by clicking or double clicking on it and simply deleting and/or typing.

# You execute markdown and code cells by hitting shift+enter.  Once you have edited a cell it must be run again for edits to take place
This is an example of a Raw cell.  Raw cells do not get executed.
# Today we will go through basics of coding in code cells as well as using matplotlib for plotting.

# ## Code Cells (python)

# ### Basics

# In[1]:


# assigning a variable
a = 3
print(a)


# In[2]:


# addition
# first we will just print the value 6 by adding 3 to a and printing.  
# Note that the value of three is still stored in variable 'a' from above
print(a+3)
# we did not save the value 6 to any variable.  We will now save that value into b and print (should get same result)
b = a+3
print(b)


# In[3]:


# note that variables previously defined in cells above will save their values in subsequent cells 
#(as long as cells above have been run)
print(a,b)


# In[4]:


# Subtraction and multiplication are similar to addition
print(3*4)
a = 3*4
print(a)
print(3-4)
a = 3-4
print(a)


# In[5]:


# division is only slightly more challenging.  The division performed depends on the type of number (i.e. integer or real)
print(3//4)
print(3.0/4.0)


# ### Data types

# As was pointed out in the division code cell above, there are different data types in python (and other) programming languages.  There are five main data types: integers, floats, complex, strings, and boolean.  Of these, the three number types (integer, float and complex) will be the most relevant for us.

# Integers are exactly what they sound like.  These are used in a variety of ways and will come up frequently.  A couple of things to keep in mind: (1) integer division is different than float (/real number) division and (2) array indeces must be integers.  Integer division gives the integer number of the quotient (truncation).  In the case of 3/4, the answer should be 0.75 but integer division yields 0 because it truncates 0.75 to 0.  What do you think 5/4 would be?  ... try it.

# In[6]:


# assign an integer to a variable - the variable will now have int data type
a = 5
print(type(a))


# Floats are real numbers.  Float division is what you would expect.

# In[7]:


# assign a float to a variable - the variable will now have float data type
a = float(5)
print(type(a))


# Complex numbers are exactly what they sound like.

# In[8]:


# assign a complext to a variable (use j instead of i and don't use * separator) - the variable will now have complex data type
j=2
a = 5.0+3j
print(type(a))


# ### Operations

# In addition to addition, subtraction, multiplication and division there are some other operations that are commonly used in coding.  These include the remainder (or modulus) operation, %, and the add (subtract, divide, or multiply) by operation.  I will demonstrate these in code snippets.

# In[9]:


# remainder operation is %
# I will do 1 to 8 divided (integer division) by four and then compute the remained
print("1 divided by 4 is", 1//4, " remainder ", 1%4)
print("2 divided by 4 is", 2//4, " remainder ", 2%4)
print("3 divided by 4 is", 3//4, " remainder ", 3%4)
print("4 divided by 4 is", 4//4, " remainder ", 4%4)
print("5 divided by 4 is", 5//4, " remainder ", 5%4)
print("6 divided by 4 is", 6//4, " remainder ", 6%4)
print("7 divided by 4 is", 7//4, " remainder ", 7%4)
print("8 divided by 4 is", 8//4, " remainder ", 8%4)


# In[10]:


# the add to or multiply by operation
# start by assigning a variable
a = 4
print("a=",a)
# now add two to a:
a = a + 2
print("a=",a)


# In[11]:


# This operation can be achieved using a shorthand of +=
a = 4
print("a=",a)
# now add two to a:
a += 2
print("a=",a)


# ### Vectors and matrices

# The next math concepts to understand how to perform in code are vectors and matrices.  Python is a very common coding language and is actively developed in a variety of areas.  As such it has a number of ''libraries'' or bits of code to do common things that are already written.  If we import those libraries we can utilize these functions.  The most common example it the numpy library (http://www.numpy.org/).  This library contains a lot of functions to perform vector, matrix and other common mathematical manipulations. The first step is to import numpy. 

# In[12]:


# in this code cell we will import numpy
import numpy as np


# In[13]:


# defining an array/vector
a = np.array([2,3,1])
print(a)


# In[14]:


# vectors have some attributes that might be of interest
print("Size of array a:", a.size)
print("Data type of array a:", a.dtype)


# For a list of all attributes of arrays see: https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.ndarray.html

# In[15]:


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


# In[16]:


# defining a matrix - these are just 2D arrays
a = np.matrix([[1,2],[3,4]])
b = np.matrix([[1,2],[3,4]],dtype=float)
print(a)
print("shape of a:",a.shape)
print("size of a:",a.size)
print("data type of a:", a.dtype, "data type of b:", b.dtype)
print("a.T=",a.T)


# In[17]:


# manipulating matrices
c = np.matrix([[5,6],[-1,0]],dtype=float)
print(c)
print("b*c=",b*c)
print("b*c=",np.dot(b,c))
print("b*cT=",np.dot(b,c.T))


# In[18]:


# elements in a matrix


# ### For loops in python

# Loops are a common structure in coding/scripting.  

# In[19]:


for i in range(10):
    print(i)


# Loops are often used to compute sums.  For example if we want to know $\sum_{i=0}^{10} i$

# In[20]:


partialSum = 0
for i in range(11):
    partialSum += i
    print(i, partialSum)
print("Final sum is:", partialSum)


# Loops are often used to compute sums.  For example if we want to know $\sum_{i=1}^{10} i$

# In[21]:


partialSum = 0
for i in range(1,11):
    partialSum += i
    print(i, partialSum)
print("Final sum is:", partialSum)


# ## Plotting with matplotlib

# Let's plot $f(x) = x^2$.  To have a computer do this, we must explicitly define our $x$ values and compute the corresponding $f(x)$ value. 

# In[22]:


# create a list of x values
xValues = []   # initialize list
for i in range(100):
    xValues.append(i*0.1)
print(xValues)


# In[23]:


# determine f(x) values
yValues = []   # initialize list
for i in range(100):
    yValues.append(xValues[i]*xValues[i])   # f(x) = x^2 = x*x
print(yValues)


# In[24]:


# plot
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
plt.plot(xValues,yValues)


# In[ ]:




