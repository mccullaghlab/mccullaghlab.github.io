#!/usr/bin/env python
# coding: utf-8

# # First Derivative and Maxwell Relations

# ## Learning goals:
# 
# After this section, student's should be able to:
# 
# 1. Identify the two first partial derivative relationships of each Thermodynamic energy function
# 2. Define a Maxwell relation
# 3. Identify the Maxwell relation in each Thermodynamic energy function (of two variables)

# ## First Partial Derivatives of Thermodynamic Energy Functions

# We have already seen the differential forms of each of the four Thermodynamic energy function and how these can be transformed from one to the other using Legendre transformations.  As a reminder, here are the four differential forms:
# 
# \begin{align}
# dU &= TdS - PdV \\
# dH &= TdS + VdP \\
# dA &= -SdT - PdV \\
# dG &= -SdT + VdP
# \end{align}
# 
# Recall that for a differntiable function of two variables, e.g. $f(x,y)$, the *exact* or *total* differential is given as
# 
# \begin{align}
# df &= \left( \frac{\partial f}{\partial x}\right)_ydx + \left( \frac{\partial f}{\partial y}\right)_xdy
# \end{align}
# 
# Let's consider each Thermodynamic energy function as a function of its natural variables and compare the two forms of its differential.

# ### Internal Energy
# 
# For internal energy, $U(S,V)$, we have the following two forms for its differential:
# 
# \begin{align}
# dU &= TdS - PdV \\
#  &= \left( \frac{\partial U}{\partial S}\right)_VdS + \left( \frac{\partial U}{\partial V}\right)_SdV
# \end{align}
# 
# where the first one comes from the standard form of the internal energy and the second is just the mathematical relationship between a differential and the partials of the independent variables.  Notice what these equalities imply:
# 
# \begin{align}
# T & = \left( \frac{\partial U}{\partial S}\right)_V \\
# -P &=\left( \frac{\partial U}{\partial V}\right)_S
# \end{align}

# ### Enthalpy
# 
# For enthalpy, $H(S,P)$, we have the following two forms for its differential:
# 
# \begin{align}
# dH &= TdS + VdP \\
#  &= \left( \frac{\partial H}{\partial S}\right)_PdS + \left( \frac{\partial H}{\partial P}\right)_SdP
# \end{align}
# 
# where the first one comes from the standard form of the enthalpy and the second is just the mathematical relationship between a differential and the partials of the independent variables.  Notice what these equalities imply:
# 
# \begin{align}
# T & = \left( \frac{\partial H}{\partial S}\right)_P \\
# V &=\left( \frac{\partial H}{\partial P}\right)_S
# \end{align}

# ### Helmholtz Free Energy
# 
# For Helmholtz free energy, $A(T,V)$, we have the following two forms for its differential:
# 
# \begin{align}
# dA &= -SdT - PdV \\
#  &= \left( \frac{\partial A}{\partial T}\right)_VdT + \left( \frac{\partial A}{\partial V}\right)_TdV
# \end{align}
# 
# where the first one comes from the standard form of the Helmholtz free energy and the second is just the mathematical relationship between a differential and the partials of the independent variables.  Notice what these equalities imply:
# 
# \begin{align}
# -S & = \left( \frac{\partial A}{\partial T}\right)_V \\
# -P &=\left( \frac{\partial A}{\partial V}\right)_T
# \end{align}

# ### Gibbs Free Energy
# 
# For Gibbs free energy, $G(T,P)$, we have the following two forms for its differential:
# 
# \begin{align}
# dG &= -SdT + VdP \\
#  &= \left( \frac{\partial G}{\partial T}\right)_PdT + \left( \frac{\partial G}{\partial P}\right)_TdP
# \end{align}
# 
# where the first one comes from the standard form of the Helmholtz free energy and the second is just the mathematical relationship between a differential and the partials of the independent variables.  Notice what these equalities imply:
# 
# \begin{align}
# -S & = \left( \frac{\partial G}{\partial T}\right)_P \\
# V &=\left( \frac{\partial G}{\partial P}\right)_T
# \end{align}

# ### Summary of first derivative relationships
# 
# | Energy Function                    | Differential Form   $ $     $ $ $ $      $ $ | First Derivative Relationships |         
# | :--------------------------------- | :--------------------------------------------------- | :--------------------------------------------------- | 
# | U - Internal Energy                | $ dU = TdS - PdV      $          | $T  = \left( \frac{\partial U}{\partial S}\right)_V$ and $-P =\left( \frac{\partial U}{\partial V}\right)_S$ |
# | H - Enthalpy                       | $dH = TdS + VdP$                 | $T  = \left( \frac{\partial H}{\partial S}\right)_P$ and $V =\left( \frac{\partial H}{\partial P}\right)_S$|
# | A - Helmholtz Free Energy          | $dA = -SdT - PdV$                | $-S  = \left( \frac{\partial A}{\partial T}\right)_V$ and $-P =\left( \frac{\partial A}{\partial V}\right)_T$|
# | G - Gibbs Free Energy              | $dG = -SdT + VdP$                | $-S  = \left( \frac{\partial G}{\partial T}\right)_P$ and $V =\left( \frac{\partial G}{\partial P}\right)_T$|

# ## Maxwell Relations
# 
# The Maxwell relations are the equivalencies of the second cross partials of the Thermodynamic energy functions

# The Maxwell relations rely on an additional mathematical property of *exact* differentials. Namely, that for an exact differential of (at least) two independent variables, the ***second cross paritial derivatives*** are equivalent. 
# 
# What is a ***second cross partial derivative***?  Let's again consider a function, $f(x,y)$, of independent variables $x$ and $y$.  The differential, $df$, is expressed in terms of first partial derivatives of each independent variable:
# 
# \begin{align}
# df &= \left( \frac{\partial f}{\partial x}\right)_ydx + \left( \frac{\partial f}{\partial y}\right)_xdy
# \end{align}
# 
# where $\left( \frac{\partial f}{\partial x}\right)$ is the first partial derivative of $f$ with respect to (wrt) $x$ and $\left( \frac{\partial f}{\partial y}\right)$ is the first partial derivative of $f$ wrt $y$.  A ***second cross paritial derivative*** is the derivative of a first partial derivative wrt the other independent variable.  In this case, there are two cross partials:
# 
# \begin{align}
# \left(\frac{\partial}{\partial y} \left( \frac{\partial f}{\partial x}\right)_y\right)_x
# &= \left( \frac{\partial^2 f}{\partial y\partial x}\right) \\
# \left(\frac{\partial}{\partial x} \left( \frac{\partial f}{\partial y}\right)_x\right)_y &= \left( \frac{\partial^2 f}{\partial x\partial y}\right)
# \end{align}
# 
# This is as opposed to the second partial derivatives wrt to the same variables
# 
# \begin{align}
# \frac{\partial}{\partial x} \left( \frac{\partial f}{\partial x}\right) 
# &= \left( \frac{\partial^2 f}{\partial x^2}\right) \\
# \frac{\partial}{\partial y} \left( \frac{\partial f}{\partial y}\right) &= \left( \frac{\partial^2 f}{\partial y^2}\right)
# \end{align}

# The ***second cross paritial derivatives*** of an *exact* differential are equivalent, meaning:
# 
# \begin{align}
#  \left( \frac{\partial^2 f}{\partial y\partial x}\right) &= \left( \frac{\partial^2 f}{\partial x\partial y}\right)
# \end{align}

# ### Example: Second Cross Partials are Equivalent
# 
# Consider the function $f(x,y) = 2xy + x^2 + \frac{y^3}{x}$.  Show that the second cross partials are equivalent.
# 
# Strategy: Determine each cross partial by first determining the first partial derivatives and subsequently differentiating those wrt the other variable.
# 
# \begin{align}
# \left( \frac{\partial f}{\partial x}\right)_y &= 2y + 2x - \frac{y^3}{x^2} \\
# \left( \frac{\partial f}{\partial y}\right)_y &= 2x + \frac{3y^2}{x} \\
# \end{align}
# 
# Now compute the second cross partials
# 
# \begin{align}
# \left( \frac{\partial ^2f}{\partial y\partial x}\right) &= 2  - \frac{3y^2}{x^2} \\
# \left( \frac{\partial ^2f}{\partial x\partial y}\right)_y &= 2 - \frac{3y^2}{x^2} \\
# \end{align}
# 
# We see that these two are equivalent.

# ### Thermodynamic Energy Functions
# 
# So what does this mean for the Thermodynamic energy functions?  We will consider the internal energy and extrapolate to the rest
# 
# \begin{align}
# dU =&  \left( \frac{\partial U}{\partial S}\right)_VdS + \left( \frac{\partial U}{\partial V}\right)_SdV \\
# \Rightarrow \left(\frac{\partial}{\partial V}\left( \frac{\partial U}{\partial S}\right)_V\right)_S =& \left(\frac{\partial}{\partial S}\left( \frac{\partial U}{\partial V}\right)_S\right)_V
# \end{align}
# 
# The real power of this comes from plugging in for the first derivatives of the energy functions (namely $T  = \left( \frac{\partial U}{\partial S}\right)_V$ and $-P =\left( \frac{\partial U}{\partial V}\right)_S$ in this case)
# 
# \begin{align}
#  \left(\frac{\partial}{\partial V}\left(T\right)\right)_S =& \left(\frac{\partial}{\partial S}\left( -P\right)\right)_V\\
#  \left(\frac{\partial T}{\partial V}\right)_S =& -\left(\frac{\partial P}{\partial S}\right)_V
# \end{align}
# 
# Similar relationship can be derived from the differential forms of $H$, $A$, $G$.  These are summmarized in following table.

# ### Summary of Maxwell Relations
# 
# | Energy Function                    | Differential Form   $$ $$             | Maxwell   Relations $$ $$|         
# | :--------------------------------- | :--------------------------------------------------- | :--------------------------------------------------- | 
# | U - Internal Energy                | $ dU = TdS - PdV      $          | $\left(\frac{\partial T}{\partial V}\right)_S = -\left(\frac{\partial P}{\partial S}\right)_V$ |
# | H - Enthalpy                       | $dH = TdS + VdP$                 | $\left(\frac{\partial T}{\partial P}\right)_S = \left(\frac{\partial V}{\partial S}\right)_P$ |
# | A - Helmholtz Free Energy          | $dA = -SdT - PdV$                | $\left(\frac{\partial S}{\partial V}\right)_T = \left(\frac{\partial P}{\partial T}\right)_V$ |
# | G - Gibbs Free Energy              | $dG = -SdT + VdP$                | $\left(\frac{\partial S}{\partial P}\right)_T = -\left(\frac{\partial V}{\partial T}\right)_P$ |
