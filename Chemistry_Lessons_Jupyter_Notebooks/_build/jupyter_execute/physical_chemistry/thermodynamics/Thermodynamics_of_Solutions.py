#!/usr/bin/env python
# coding: utf-8

# # Thermodynamics of Solutions

# ## Learning goals
# 
# After working through this notebook, you will be able to:
# 
# 1. Express each of the four Thermodynamic energy functions for multicomponent systems
# 2. Define chemical potential both mathematically and qualitatively
# 3. Derive the Gibbs-Duhem equation.
# 4. Use the Gibbs-Duhem equation to relate the chemical potentials of two component solutions.

# ## Coding Concepts:
# 
# The following coding concepts are used in this notebook:
# 
# 1. [Variables](../../coding_concepts/variables.ipynb)

# ## What makes solutions different from gasses?
# 
# Any number of things.  But here are a few:
# 
# 1. Density of a solution is larger than that of a gas.
# 2. Solutions contain multiple species (by definition).
# 3. PV work doesn't seem so relevant for solutions as compared to gasses.
# 4. Solutions are typically held under constant $P$ and $T$ conditions.

# ## A Two Component Solution
# 
# Consider a two component solution with components $A$ (e.g. H$_2$O(l)) and $B$ (e.g. CH$_3$OH(l)).  How does the chemical potential of one component compare to the other under constant $T$ and $P$ conditions?
# 
# Unlike the case of a two phase equilibrium, the components of this solution cannot interconvert.  Thus, we do not necessarily have that the chemical potentials are equal at equilibrium.

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

# ## The Gibbs-Duhem Equation for a Two Component Solution
# 
# For a two component system the Gibbs-Duhem equation is
# 
# \begin{align}
# 0 =& -SdT + VdP - N_Ad\mu_A - N_Bd\mu_B
# \end{align}
# 
# Under constant $T$ and $P$ conditions this yields
# 
# \begin{align}
# N_Ad\mu_A &=- N_Bd\mu_B 
# \end{align}
# 
# If we now divide by $N = N_A + N_B$ (constant total number of moles) we get
# 
# \begin{align}
# x_Ad\mu_A &=- x_Bd\mu_B 
# \end{align}
# 
# where $x_A = \frac{N_A}{N}$ is the mole fraction of component $A$ and $x_B = \frac{N_B}{N}$ is the mole fraction of component $B$.  This equation can be rearranged to solve for $d\mu_A$ or $d\mu_B$ which is one of the more common forms of the Gibbs-Duhem equation.
# 
# \begin{align}
# d\mu_A &=- \frac{x_B}{x_A}d\mu_B 
# \end{align}
# 
# This equation can be used to solve for the chemical potential of one species if you know the chemical potential of the other species as a function of composition ($x_A$ and/or $x_B$).

# ## Example: The Chemical Potentials of a Two Component Solution
# 
# Consider a solution of two components held under constant $T$ and $P$.  The chemical potential of component $A$ is determined to obey the functional form $\mu_A = \mu_A^* + RT\ln x_A$, where $\mu_A^*$ is the chemical potential of pure A, and $x_A$ is the mole fraction of species A.  Show that the chemical potential of species $B$ must be given as $\mu_B = \mu_B^* + RT\ln x_B$, where $\mu_B^*$ is the chemical potential of pure B, and $x_B$ is the mole fraction of species B. 
# 

# We start with the Gibbs-Duhem relation at constant $P$ and $T$ for a binary system:
# 	\begin{eqnarray}
# 	d\mu_B = -\frac{x_A}{x_B}d\mu_A
# 	\end{eqnarray}
# We now determine the differential of $\mu_A$ from the equation given:
# 	\begin{eqnarray}
# 	d\mu_A &=& d(\mu_A^* + RT\ln x_A) \\
# 	&=& \frac{RT}{x_A}dx_A \\
# 	&=& -\frac{RT}{x_A}dx_B,
# 	\end{eqnarray}
# where the last equality holds because $x_A + x_B = 1$, or $dx_A = -dx_B$.  Now plug into Gibbs-Duhem and integrate:
# \begin{eqnarray}
# 	d\mu_B &=& -\frac{x_A}{x_B}\left( -\frac{RT}{x_A}dx_B \right) \\
# 	&=& \frac{RT}{x_B}dx_B \\
# 	\Rightarrow \int_1^{x_B}d\mu_B &=& \int_1^{x_B} \frac{RT}{x_B}dx_B \\
# 	\mu_B - \mu^*_B &=& RT\ln(x_B) \\
# 	\Rightarrow \mu_B  &=& \mu^*_B + RT\ln(x_B)
# \end{eqnarray}

# ## Solutions in Equilibria with their Vapors

# The discussion thus far has been focused on single phase solutions.  But solutions can also be in equilibrium with their vapors.  Take the example of a solution of water and methanol.  The vapor above this solution will contain both water and methanol.  But, again, clearly water and methanol are not inter converting so what can we say about this situation?
# 
# In this case, at constant $T$ and $P$ and at equilibrium, the chemical potentials of the vapor and the liquid will be equal for each component.  

# In[ ]:




