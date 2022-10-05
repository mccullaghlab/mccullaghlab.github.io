#!/usr/bin/env python
# coding: utf-8

# # Jupyter Notebook Tutorial

# ## What is a Jupyter Notebook?
# 
# A Jupyter Notebook, such as the one your are currently looking at, is a collection of rich text elements (text, images, equations, etc) and computer code (typically python).  These notebooks are stored as files analagous to html or xml that can be executed/intepreted by a Jupyter Notebook App (comes with Anaconda python).  The code can be executed in an active kernel either through the App or using a online server/hub, if setup properly.  

# ## A Notebook is a Collection of Cells

# A Jupyter Notebook is a collection of cells.  Cells can be added and deleted using the buttons on the top or by keyboard shortcuts.

# Each cell can be one of three types.  These types can be changed using the dropdown menu above (when in an active kernel)
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

# In an active kernel, you can modify a markdown cell by double clicking on the cell.

# You execute markdown and code cells by hitting shift+enter or using the buttons on the top.  Once you have edited a cell it must be run again for edits to take place
This is an example of a Raw cell.  Raw cells do not get executed.
# ## Markdown Cells

# Markdown cells are the richt text cells, such as this one.  Raw text is typed in and then it is passed through a markdown interpreter to be formatted.  There is much that can be done with markdown but here we highlight three important aspects:
# 
# 1. Equations
# 2. Embedding images
# 3. Tables

# ### Equations in Markdown Cells

# Equations can be embedded in markdown cells using LaTeX syntax.  The following latex math functions are understood
# 
# `$ \phi = \pi $`  
# $\phi = \pi$
# 
# `\begin{equation} 
# x = 2 \
# \end{equation}` \begin{equation}x=2\end{equation}
# 
# `\begin{eqnarray} 
# x &=& 2\\
# x + y&=& 5
# \end{eqnarray}` 
# 
# \begin{eqnarray}
# x&=&2 \\
# x + y&=&5
# \end{eqnarray}
# 
# `\begin{align} 
# x &=& 2\\
# x + y&=& 5
# \end{align}` 
# 
# \begin{align}
# x&=2 \\
# x + y&=5
# \end{align}

# ### Embedding Images

# Images can be embedded in markdown cells.  This is similar to how you would embed an image in html.  The image file must be able to be found, locally.  Here I use the logo of my book as an example.
# 
# `<img src="logo.png" width="300" align="center">`
# 
# <img src="logo.png" width="300" align="center">

# ### Tables

# Tables can be created but, somewhat interestingly, it does not take latex format.  Markdown as its own format
# 
# `| Energy Function                  | Definition                    | Differential Form                | Natural Variables |         
# | :---------------------------- | :-------------------------------------- | :------------------------------ | :------------------------------ | 
# | U - Internal Energy              | $\Delta U = q + w $           | $dU=TdS-PdV$                 | $S$ and $V$|
# | H - Enthalpy                     | $H = U + PV$                  | $dH = TdS + VdP$                 | $S$ and $P$|
# | A - Helmholtz Free Energy        | $A = U - TS$                  | $dA = -SdT - PdV$                | $T$ and $V$|
# | G - Gibbs Free Energy            | $G = H - TS$                  | $dG = -SdT + VdP$                | $T$ and $P$|`
# 
# | Energy Function           | Definition                 | Differential Form          | Natural Variables |         
# | :------------------------ | :------------------------- | :------------------------- | :------------------------ | 
# | U - Internal Energy                | $\Delta U = q + w $           | $dU=TdS-PdV$                 | $S$ and $V$|
# | H - Enthalpy                       | $H = U + PV$                  | $dH = TdS + VdP$                 | $S$ and $P$|
# | A - Helmholtz Free Energy          | $A = U - TS$                  | $dA = -SdT - PdV$                | $T$ and $V$|
# | G - Gibbs Free Energy              | $G = H - TS$                  | $dG = -SdT + VdP$                | $T$ and $P$|

# ## Code Cells (python)

# ### Basics

# In[1]:


# assigning a variable
a = 3
print(a)


# ### Libraries

# Libraries can be imported but they must be available to the python interpreter/kernel.  This can sometimes be challenging if launching remotely.

# In[2]:


import numpy as np

