#!/usr/bin/env python
# coding: utf-8

# # Thermodynamic Properties from the Partition Function

# ## Motivation

# Thermoydnamic properties of a system, such as total internal energy, pressure and heat capacity, can be determined from the partition function for the system.  These notes will discuss how to do that.

# ## Learning Goals

# After reading these notes and attending lecture, students should be able to:
# 1. Write out the expression for internal energy as a function of the Canonical partition function,
# 2. Write out the expression for pressure in terms of the Canonical partition function,
# 3. Write out the expression for heat capacity at constant volumne in terms of the Canonical partition function.

# ## Overview

# For a system with fixed number of particles, $N$, volume, $V$, and temperature $T$, the partition function is a sum of the Boltzmann factors for all possible energy levels:
# 
# \begin{equation}
# Q(N,V,T) = \sum_i e^{-\beta E_i}.
# \end{equation}
# 
# where $\beta = \frac{1}{k_BT}$.
# 
# Various macroscopic properties of the system can be determined from partition function, $Q(N,V,T)$.  In these notes we will show that 
# 
# \begin{eqnarray}
# \langle E \rangle &=& -\left( \frac{\partial \ln Q}{\partial \beta} \right) \\
# P &=& k_BT \left( \frac{\partial \ln Q}{\partial V} \right) \\
# C_V &=& k_B\beta^2\left( \frac{\partial^2 \ln Q}{\partial \beta^2} \right)
# \end{eqnarray}

# ## Internal Energy

# The internal energy for a system, often referred to by the symbol $U$ in classical Thermodynamics, can be determined from the partition function for a system.  The first step is to accept that the internal energy, $U$, is equivalent to the average energy of the system, $\langle E \rangle$.  
# 
# The average energy can be expressed as a weighted average over a given probabilty distribution in the standard way,
# \begin{equation}
# \langle E \rangle = \sum_i E_i P_i.
# \end{equation}
# Recall that, in the context of a system with constant $N$, $V$, and $T$, the probability of state $i$ is given as a ratio of the Boltzmann factor for that state and the partition function,
# \begin{equation}
# P_i = \frac{e^{-\beta E_i}}{Q}.
# \end{equation}
# Substituting this in to the equation for $\langle E \rangle$ gives
# \begin{equation}
# \langle E \rangle = \frac{\sum_i E_ie^{-\beta E_i}}{Q}.
# \end{equation}

# The key step in simplifying this equation is to recognize that the term $E_ie^{-\beta E_i}$ looks like $xe^{\alpha x}$ where the prefactor of the exponent happens to be the derivative of the argument of the exponent with respect to $\alpha$.  In terms of the energy, this means
# \begin{equation}
# E_ie^{-\beta E_i} = -\frac{\partial}{\partial \beta} \left(e^{-\beta E_i}\right).
# \end{equation}
# In this equation, the symbol $\partial$ refers to the partial derivative, instead of the total derivative denoted by a straight-backed $d$.  Partial derivatives show-up when we consider derivatives of functions of more than one variable.  $Q$, for instance, we refer to as $Q(N,V,T)$ meaning that it is a function of $N$, $V$ and $T$.  This is analagous to writing a function, $f(x,y,z)$, as a function of $x$, $y$, and $z$.  If we consider taking the derivative of $f$ we must consider which variable ($x$, $y$, or $z$) to differentiate by.  
# 
# 

# We now substitute the key relationshop into the average energy equation to get
# \begin{eqnarray}
# \langle E \rangle &=& \frac{\sum_i E_ie^{-\beta E_i}}{Q} \\
# &=& \frac{\sum_i -\frac{\partial}{\partial \beta} \left(e^{-\beta E_i}\right)}{Q} \\
# &=& \frac{-\frac{\partial}{\partial \beta} \sum_i e^{-\beta E_i}}{Q},
# \end{eqnarray}
# where in the last step we recognize that the derivative operator can be distributed over the sum.  Now we recognized that $Q = \sum_i e^{-\beta E_i}$ to get
# \begin{eqnarray}
# \langle E \rangle &=& \frac{-\frac{\partial}{\partial \beta} Q}{Q} \\
# &=& -\frac{1}{Q} \frac{\partial Q}{\partial \beta}.
# \end{eqnarray}
# 
# The final step in the derivation is to recognize that $\frac{1}{Q} \frac{\partial Q}{\partial \beta} = \frac{\partial \ln Q}{\partial \beta}$.  This can be hard to recognize but it should follow if you do it in the reverse.  That is, ask 
# \begin{equation}
# \frac{\partial \ln f(x)}{\partial x} = ?
# \end{equation}
# Using the chain rule for differentiation, this becomes
# \begin{equation}
# \frac{\partial \ln f(x)}{\partial x} = \frac{1}{f(x)} \frac{\partial f(x)}{\partial x}
# \end{equation}

# So, finally, we get that
# \begin{eqnarray}
# \langle E \rangle &=&  -\frac{1}{Q} \frac{\partial Q}{\partial \beta} \\
# &=& - \frac{\partial \ln Q}{\partial \beta}.
# \end{eqnarray}

# ## Pressure

# To get the pressure from the partition function we will use the average energy expression derived above.  Addtionally, we will need to employ a fact from classical Thermodynamics that we haven't yet seen.  Namely that
# \begin{equation}
# p = -\left(\frac{\partial U}{\partial V}\right)_{N}.
# \end{equation}
# We will show why this is true in a few weeks when we discuss the interal energy from a classical Thermodynamics perspective.
# 
# We start in an analogous manner as for $\langle E \rangle$
# \begin{eqnarray}
# p &=& \frac{\sum_i p_i e^{-\beta E_i}}{Q} \\
# &=&-\frac{\sum_i \frac{\partial E_i}{\partial V} e^{-\beta E_i}}{Q},
# \end{eqnarray}
# where I have plugged in that $p_i = -\frac{\partial E_i}{\partial V}$ in the last step.  
# 
# We now investigate the volume derivate of the parition function:
# \begin{eqnarray}
# \frac{\partial Q}{\partial V} &=& \frac{\partial }{\partial V} \sum_i e^{-\beta E_i} \\
# &=& \sum_i \frac{\partial e^{-\beta E_i} }{\partial V}  \\
# &=& \sum_i \frac{\partial (-\beta E_i)}{\partial V} e^{-\beta E_i}  \\
# &=& -\beta \sum_i \frac{\partial  E_i}{\partial V} e^{-\beta E_i},
# \end{eqnarray}
# or
# \begin{eqnarray}
# \sum_i \frac{\partial  E_i}{\partial V} e^{-\beta E_i} = -\frac{1}{\beta}\frac{\partial Q}{\partial V}.
# \end{eqnarray}
# 
# We can plug this result into the previous expression for $p$ to get
# \begin{eqnarray}
# p &=& -\frac{\sum_i \frac{\partial E_i}{\partial V} e^{-\beta E_i}}{Q} \\
# &=& -\frac{-\frac{1}{\beta}\frac{\partial Q}{\partial V}}{Q} \\
# &=& \frac{1}{\beta}\frac{\frac{\partial Q}{\partial V}}{Q} \\
# &=& \frac{1}{\beta}\frac{\partial \ln Q}{\partial V}
# \end{eqnarray}

# ## Heat Capacity

# The heat capacity for a system measures the amount of energy it takes to raise the system by one Kelvin (or, equivalently, one degree Celsius).  Heat capacity can be determined under constant volume condidtions, denoted $C_V$, or constant pressure conditions, denoted $C_P$.  Here, we will concern ourselves with $C_V$.  From classical Thermodynamics we know that the heat capacity at constant volume is
# \begin{equation}
# C_V = \left(\frac{\partial U}{\partial T}\right)_V.
# \end{equation}
# 
# In order to determine this quantity from a partition function, we recognize that $U = \langle E \rangle$ yielding
# \begin{eqnarray}
# C_V &=& \left(\frac{\partial \langle E \rangle}{\partial T}\right)_V \\
# &=& \left(\frac{\partial}{\partial T} \left(- \frac{\partial \ln Q}{\partial \beta}\right)\right)_V
# \end{eqnarray}
# 
# The difficultly here is resolving the issue of $\beta$ vs $T$.  We do this by starting with the equation for $\beta$ in terms of $T$ and differentiating both sides
# \begin{eqnarray}
# \beta &=& \frac{1}{k_BT} \\
# \Rightarrow d\beta &=& - \frac{1}{k_BT^2}dT \\
# &=& -k_B \beta^2 dT \\
# \Rightarrow dT &=& -\frac{d\beta}{k_B \beta^2}
# \end{eqnarray}
# 
# Substituting this into the above equation for $C_V$ yields
# \begin{eqnarray}
# C_V &=& -k_B \beta^2 \left(\frac{\partial}{\partial \beta} \left(- \frac{\partial \ln Q}{\partial \beta}\right)\right)_V \\
# &=&k_B \beta^2 \left( \frac{\partial^2 \ln Q}{\partial \beta^2}\right)_V
# \end{eqnarray}

# ## Example 1: Average Energy of a Monatomic Ideal Gas

# The partition function for an ideal gas is with fixed number of particles, $N$, volume, $V$, and temperature, $T$, is
# 
# \begin{equation}
# Q=\frac{q^N}{N!}.
# \end{equation}
# 
# For a monatomic gas, we have that
# 
# \begin{equation}
# q=\left(\frac{2\pi mk_BT}{h^2}\right)^{\frac{3}{2}}V
# \end{equation}
# 
# What is the average energy, $\langle E \rangle$, of a monatomic ideal gas?

# From the notes above, we know that
# \begin{equation}
# \langle E \rangle = -\left( \frac{\partial \ln Q}{\partial \beta}\right)
# \end{equation}
# So, the strategy to solve this problem will be:
# 
# 1. get an equation for $\ln Q$ in terms of $\beta$ 
# 2. differentiate that equation with respect to $\beta$.  
# 
# 
# Doing this problem stepwise, in particular taking the $\ln Q$ outside of the derivative can really help simplify things as we will have quite a few terms that do not dependend on $\beta$

# ### 1. Equation for $\ln Q$

# We start by getting an overall equation for $Q$ and then taking the natural log
# 
# \begin{eqnarray}
# Q &=& \frac{q^N}{N!} \\
# &=&\frac{\left[ \left(\frac{2\pi mk_BT}{h^2}\right)^{\frac{3}{2}}V \right]^N}{N!} \\
# &=& \frac{\left[ \left(\frac{2\pi m}{\beta h^2}\right)^{\frac{3}{2}}V \right]^N}{N!} \\
# &=& \frac{ \left(\frac{2\pi m}{\beta h^2}\right)^{\frac{3N}{2}}V^N }{N!},
# \end{eqnarray}
# where I substituted that $k_BT = \frac{1}{\beta}$ to get Q in terms of $\beta$ and then distributed the power of $N$.  Now we take the natural log and isolate term(s) that depend on $\beta$ (we will make heavy use of the rules of logarithms $\ln(A\cdot B) = \ln A + \ln B$ and $\ln(A^N) = N\ln A$)
# \begin{eqnarray}
# \ln Q &=& \ln \left[\frac{\left(\frac{2\pi m}{\beta h^2}\right)^{\frac{3N}{2}}V^N }{N!}\right] \\
# &=& \ln\left[ \left(\frac{2\pi m}{\beta h^2}\right)^{\frac{3N}{2}}\right] + \ln V^N - \ln N! \\
# &=& \frac{3N}{2}\ln \left(\frac{2\pi m}{\beta h^2}\right) + \ln V^N - \ln N! \\
# &=& \frac{3N}{2}\ln\frac{1}{\beta} + \frac{3N}{2}\ln \left(\frac{2\pi m}{h^2}\right)+ \ln V^N - \ln N! \\
# &=& -\frac{3N}{2}\ln\beta + \frac{3N}{2}\ln \left(\frac{2\pi m}{h^2}\right)+ \ln V^N - \ln N!
# \end{eqnarray}
# 
# Notice that only the first term of the last equation for $\ln Q$ depends on $\beta$.  Thus, when we differentiate with respect to (w.r.t) $\beta$ the last three terms will all yield zeros.  

# ### 2. Differentiate with respect to $\beta$

# \begin{eqnarray}
# \langle E \rangle &=& - \left( \frac{\partial \ln Q}{\partial \beta}\right) \\
# &=& -\frac{\partial}{\partial \beta} \left[ -\frac{3N}{2}\ln\beta + \frac{3N}{2}\ln \left(\frac{2\pi m}{h^2}\right)+ \ln V^N - \ln N! \right] \\
# &=& \frac{3N}{2}\frac{\partial}{\partial \beta} \ln\beta \\
# &=& \frac{3N}{2}\frac{1}{\beta} \\
# &=& \frac{3Nk_BT}{2}
# \end{eqnarray}

# ## Example 2: Heat Capacity of a Monatomic Ideal Gas
# <div class="alert alert-block alert-info">
# What is the heat capacity at constant volume, $C_V$, of a monoatomic ideal gas?
# </div>

# ## Example 3: Pressure of a Monatomic Ideal Gas
# <div class="alert alert-block alert-info">
# What is the pressure, $p$, of a monoatomic ideal gas?
# </div>
