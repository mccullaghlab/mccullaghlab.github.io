#!/usr/bin/env python
# coding: utf-8

# # Variables in Python

# ## Learning Goals
# 
# After going through this notebook, you should be able to:
# 1. Identify a variable in a python code snippet of a jupyter notebook
# 2. Assign values to a variable in python
# 3. Perform mathematical operations to variables in python
# 4. Assign variables of different data types

# ## Assigning a variable

# A variable in code is much like a variable in math.  Unlike in math, variables must be assigned specific values before the code operates.  We will start by assigning a value to a variable called `a` and then print that value.

# In[1]:


# assigning a variable
a = 3
print(a)


# Operations such as addition can be performed on variables with numbers (integers, floats, doubles, etc).  New variables can be assigned with the resulting value or the value can be simply printed or used on the fly

# In[2]:


# addition
# first we will just print the value 6 by adding 3 to a and printing.  
# Note that the value of three is still stored in variable 'a' from above
print(a+3)
# we did not save the value 6 to any variable.  We will now save that value into b and print (should get same result)
b = a+3
print(b)


# In a Jupyter Notebook, variables that have been assigned in cells that have already run (not necessarily in one only above the current one) can be recalled in cells that are run subsequently.

# In[3]:


# note that variables previously defined in cells above will save their values in subsequent cells 
#(as long as cells above have been run)
print(a,b)


# Other operations such as multiplication and addition are available.  The specific result (and in the case of python specifically) and/or operators available depend on the particular data types.  More on data types below.

# In[4]:


# Subtraction and multiplication are similar to addition
print(3*4)
a = 3*4
print(a)
print(3-4)
a = 3-4
print(a)


# Variables can also be updated in value.  For example, we will assign `a` a value and then assign `a` that same values multiplied by 4.

# In[5]:


a = 3
print(a)
a = a*4
print(a)


# In python, there is a shorthand for this type of operation which you will see in the code snippet below

# In[6]:


a = 3
print(a)
a *= 4
print(a)


# ## Data types

# As was pointed out in the division code cell above, there are different data types in python (and other) programming languages.  There are five main data types: integers, floats, complex, strings, and boolean.  Of these, the three number types (integer, float and complex) will be the most relevant for us.

# Integers are exactly what they sound like.  These are used in a variety of ways and will come up frequently.  A couple of things to keep in mind: (1) integer division is different than float (/real number) division and (2) array indeces must be integers.  Integer division gives the integer number of the quotient (truncation).  In the case of 3/4, the answer should be 0.75 but integer division yields 0 because it truncates 0.75 to 0.  What do you think 5/4 would be?  ... try it.

# In[7]:


# assign an integer to a variable - the variable will now have int data type
a = 5
print(a, type(a))


# Floats are real numbers.  Float division is what you would expect.

# In[8]:


# assign a float to a variable - the variable will now have float data type
a = float(5)
print(a, type(a))


# Complex numbers are exactly what they sound like.

# In[9]:


# assign a complex to a variable (use j instead of i and don't use * separator) - the variable will now have complex data type
j=2
a = 5.0+3j
print(a, type(a))


# Strings

# In[10]:


a = str(5)
print(a, type(a))


# Boolean

# In[11]:


a = bool(5)
print(a,type(a))


# ## Mathematical Operators

# In[12]:


# division is only slightly more challenging.  The division performed depends on the type of number (i.e. integer or real)
print(3//4)    # integer devision
print(3/4)     # float division


# In addition to addition, subtraction, multiplication and division there are some other operations that are commonly used in coding.  These include the remainder (or modulus) operation, %, and the add (subtract, divide, or multiply) by operation.  I will demonstrate these in code snippets.

# In[13]:


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


# In[14]:


# the add to or multiply by operation
# start by assigning a variable
a = 4
print("a=",a)
# now add two to a:
a = a + 2
print("a=",a)


# In[15]:


# This operation can be achieved using a shorthand of +=
a = 4
print("a=",a)
# now add two to a:
a += 2
print("a=",a)


# In[ ]:




