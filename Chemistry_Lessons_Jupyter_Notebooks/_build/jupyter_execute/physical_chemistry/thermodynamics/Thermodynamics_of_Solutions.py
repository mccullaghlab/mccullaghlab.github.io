#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# # Thermodynamics of Solutions

# ## Learning goals
# 
# 1. Express each of the four Thermodynamic energy functions for multicomponent systems
# 2. Define chemical potential both mathematically and qualitatively
# 3. Recognize that two component systems tend toward the component with the smaller chemical potential
# 4. Demonstrate that the equilibrium condition of a two component system requires that the chemical potentials of each component be equivalent
# 5. Recognize the Gibbs-Duhem equation and state its implication.

# ## What makes solutions different from gasses?
# 
# Any number of things.  But here are a few:
# 
# 1. Density of a solution is larger than that of a gas.
# 2. Solutions contain multiple species (by definition).
# 3. PV work doesn't seem so relevant for solutions as compared to gasses.
# 4. Solutions are typically held under constant $P$ and $T$ conditions.

# ## Thermodynamic Energy Functions for Multicomponent Systems

# Though it was never (or at least not frequently) stated, the Thermodynamic energy functions we have encountered so far are only for single component systems.  e.g. a single type of gas.  For a multicomponent system, we must consider the energy associated with each component to, for example, compute the energy difference between two different concentrations of a solution.  
# 
# We start by looking at the differential form of the internal energy
# 
# \begin{equation}
# dU = TdS - PdV
# \end{equation}
# 
# We now say that the infintesimal change in internal energy, $dU$, can be changed by changes in the number of moles/molecules of each component.  We will write this for a general $M$-component system
# 
# \begin{equation}
# dU = TdS - PdV + \sum_{i=1}^M \mu_idN_i
# \end{equation}
# 
# where $\mu_i$ is call the ***chemical potential*** of component $i$ because it is the energy per mole/molecule of component $i$ and $dN_i$ is the infintesimal change in the amount (number of molecules/moles) of component $i$.

# Notice that the formula
# 
# \begin{equation}
# dU = TdS - PdV + \sum_{i=1}^M \mu_idN_i
# \end{equation}
# 
# dictates that 
# 
# \begin{equation}
# \mu_i = \left(\frac{\partial U}{\partial N_i}\right)_{S,V,N_{j\neq i}}
# \end{equation}
# 
# or $\mu_i$ is the slope in internal energy along the dimension $N_i$.  i.e. $\mu_i$ is the internal energy per molecule of component $i$.

# Also notice that conversion to the other three Thermodynamic energy functions ($H$, $A$, and $G$) will not affect the chemical potential terms.  More to the point, given these chemical potential terms in $dU$, we can write the three other differential forms:
# 
# \begin{align}
# dH &= TdS + VdP + \sum_{i=1}^M \mu_idN_i \\
# dA &= -SdT - PdV + \sum_{i=1}^M \mu_idN_i \\
# dG &= -SdT + VdP + \sum_{i=1}^M \mu_idN_i
# \end{align}
# 
# Due to many solutions being held under constant $P$ and $T$, the Gibbs free energy function is the most useful Thermoydnamic energy function.  As such, the most common definition of the chemical potential is
# 
# \begin{equation}
# \mu_i = \left(\frac{\partial G}{\partial N_i}\right)_{T,P,N_{j\neq i}}
# \end{equation}

# ## Example: Chemical potential of a two component system
# 
# Consider a two component system with components $A$ (e.g. H$_2$O(l)) and $B$ (e.g. H$_2$O(g)) that can exchange between one another.  Show that, under constant temperature, pressure, and total number of moles, the following two statements hold. 
# 
# (a) The system will proceed towards the component with the lower chemical potential.
# 
# (b) The chemical potentials of each component are equal at equilibrium.

# We start by thinking about the physical system.  We can consider a phase equilibrium between water vapor and liquid water.  The liquid water molecules can become vapor and the vapor molecules can condense to become liquid.  The total number of molecules or moles ($N$) is fixed, however.  We can write that
# 
# \begin{align}
# N = N_A + N_B
# \end{align}
# 
# We also know that total number of moles is fixed meaning that
# 
# \begin{align}
# dN =& 0 \\
#  =& dN_A + dN_B \\
#  \Rightarrow dN_A =& -dN_B
# \end{align}
# 
# where the last equality simply tells us that the change in amount of A is negative the change in amount of B (which is true if there are only these two components in the system).

# Now onto part (a).  The key to this statement is that the "system will proceed".  What does that mean?  Under constant $P$ and $T$ it means that $dG < 0$.  So let's write out our expression for dG
# 
# \begin{align}
# dG =& -SdT + VdP + \mu_AdN_A + \mu_BdN_B \\
#  =& \mu_AdN_A + \mu_BdN_B \\
#  <& 0
# \end{align}
# 
# where I have used that $dT = dP = 0$.  We can now use the relationship $dN_A = -dN_B$ to yield
# 
# \begin{align}
# dG =& \mu_AdN_A - \mu_BdN_A < 0 \\
#  = & (\mu_A-\mu_B)dN_A < 0
# \end{align}
# 
# The last inequality can be achieved only if $dN_A < 0$ and $(\mu_A-\mu_B)>0$ or $dN_A > 0$ and $(\mu_A-\mu_B)<0$.  Mathematically, this is equivalent to saying that if the number of moles of component $A$ decreases ($dN_A < 0$), then the chemical potential of component $B$ must be less than that of $A$ ($(\mu_A-\mu_B)>0$ or, equivalently, $\mu_A>\mu_B$).  Similarly, if the number of moles of $A$ increases ($dN_A > 0$), then the chemical potential of $A$ must be less than that of $B$.

# No for part (b).  The solution to this is quite similar but we must recognize that, under these conditions, the condition of equilibirium dictates that $dG=0$.  Again, we start with the differential form of $G$:
# 
# \begin{align}
# dG = \mu_AdN_A + \mu_BdN_B = 0
# \end{align}
# 
# again, we plug in the relationship $dN_A = -dN_B$ to yield
# 
# \begin{align}
# dG &= \mu_AdN_A - \mu_BdN_A = 0 \\
# &= (\mu_A-\mu_B)dN_A = 0
# \end{align}
# 
# For the last equality to hold $dN_A$ could be zero and/or $(\mu_A-\mu_B)$ could be zero.  While $dN_A$ could be zero, it need not be zero under dynamic equilibrium.  At equilibrium, moleules of liquid water will be converting to vapor and vis versa.  Thus, $(\mu_A-\mu_B)=0$ or $\mu_A=\mu_B$ at equilibrium.

# ## The Gibbs-Duhem Equation
# 
# The Gibbs-Duhem equation demonstrates the relationship between all intensive variables of a system.  The implication of this is that if there are $M$ intesive variables of a system, you need only measure $M-1$ of them to get all of them.  Seems unimportant but it isn't.

# To derive the Gibbs-Duhem equation, we start with the fact that the chemical potential of component is the Gibbs free energy per molecule (or mole) of that component.  Thus, we can write the total Gibbs free energy as
# 
# \begin{align}
# G = \sum_i \mu_iN_i
# \end{align}
# 
# We now take the differential of both sides to get
# 
# \begin{align}
# dG = \sum_i \left[ \mu_idN_i + N_id\mu_i\right]
# \end{align}
# 
# where we use the product rule on each term of the sum.  We now split this into two sums and relate it to the other equation for $dG$ that we have
# 
# \begin{align}
# dG =& \sum_i  \mu_idN_i +\sum_i N_id\mu_i \\
#  =& -SdT + VdP + \sum_i  \mu_idN_i
# \end{align}
# 
# now substract $\sum_i  \mu_idN_i$ from both sides to get the Gibbs-Duhem equation
# 
# \begin{align}
# \sum_i N_id\mu_i =& -SdT + VdP \\
# \Rightarrow 0 =& -SdT + VdP - \sum_i N_id\mu_i
# \end{align}
# 
# where $T$, $P$ and all of the $\mu_i$s are the intensive variables of the system.  The Gibbs-Duhem equation suggest that state functions exist.
