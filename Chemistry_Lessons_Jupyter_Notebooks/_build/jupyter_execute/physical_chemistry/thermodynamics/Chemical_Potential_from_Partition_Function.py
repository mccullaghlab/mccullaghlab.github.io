#!/usr/bin/env python
# coding: utf-8

# # Chemical Potential from the Partition Function

# ## Motivation
# 
# The chemical potential, or molar Gibbs free energy of a component, is an important quantity for multicomponent systems.  We saw that two phase equilibrium at constant temperature and pressure dictates that the chemical potentials of each phase be equivalent.  Here, we will show that the chemical potential can be determined from the partition function of the system.  This provides a tie between this quantity and the molecular level.

# ## Learning Goals:
# 
# After working through these notes, you will be able to:
# 
# 1. Write out the chemical potential in terms of derivatives of all four energy functions
# 2. Express the chemical potential in terms of the canonical partition function
# 3. Calculate the standard state chemical potential for an ideal gas

# ## Coding Concepts: 
# 
# None, yet.

# ## Review of the Partition Function

# Recall that the partition function, $Q$, of a system under constant $V$ and $T$ is defined as a sum over the Boltzmann factors of all energy levels of the system
# \begin{equation}
# Q = \sum_i e^{-\beta E_i},
# \end{equation}
# where $E_i$ is the energy of the system in energy level $i$.
# 
# For an ideal gas system the particles do not interact and thus we can write
# \begin{equation}
# Q = \frac{q^N}{N!},
# \end{equation}
# where $q$ is the molecular partition function defined as a sum over molecular energy levels
# \begin{equation}
# q = \sum_j e^{-\beta \epsilon_j},
# \end{equation}
# where $\epsilon_j$ is the energy of a molecule is molecular state $j$.
# 
# From the partition function we can compute important thermodynamic quantities such as
# \begin{eqnarray}
# U &=& k_BT^2\left( \frac{\partial lnQ}{\partial T}\right)_{N,V} \\
# P &=&  \frac{1}{\beta}\frac{\partial \ln Q}{\partial V} \\
# C_V &=& k_B \beta^2 \left( \frac{\partial^2 \ln Q}{\partial \beta^2}\right)_{N,V} \\
# S &=& k_BT\left( \frac{\partial lnQ}{\partial T}\right)_{N,V} + k_B\ln Q
# \end{eqnarray}

# ## Entropy from the Partition Function
# 
# Actually, we have not yet shown the last equation above.  I will do so here a bit of an aside.  We start with the Gibbs entropy formula
# \begin{equation}
# S = -k_B\sum_ip_i\ln p_i
# \end{equation}
# where $p_i$ is the probability of state $i$ and in the Canonical ensemble is given as
# \begin{equation}
# p_i = \frac{e^{-\beta E_i}}{Q}
# \end{equation}
# Plugging this formula into the Gibbs entropy equation yields
# \begin{eqnarray}
# S &=& -k_B\sum_i \frac{e^{-\beta E_i}}{Q} \ln \left( \frac{e^{-\beta E_i}}{Q} \right)  \\
# &=& -k_B\sum_i \frac{e^{-\beta E_i}}{Q} \left( -\beta E_i - \ln Q\right) \\
# &=& \beta k_B\sum_i E_i\frac{e^{-\beta E_i}}{Q} + k_B\sum_i\ln Q\frac{e^{-\beta E_i}}{Q} \\
# &=& \frac{1}{T}\langle E \rangle  + \frac{k_B\ln Q}{Q} \sum_i e^{-\beta E_i} \\
# &=& \frac{U}{T} + k_B\ln Q
# \end{eqnarray}

# ## Helmholtz Free Energy from Partition Function
# 
# Given the above equalities and the Thermodynamic definition of $A$, it is straightforward to determine an equation for $A$ in terms of the partition function.
# 
# \begin{eqnarray}
# A &=& U - TS \\
# &=& k_BT^2\left( \frac{\partial lnQ}{\partial T}\right)_{N,V} - T\left( k_BT\left( \frac{\partial lnQ}{\partial T}\right)_{N,V} + k_B\ln Q \right) \\
# &=& -k_BT\ln Q 
# \end{eqnarray}

# ## Chemical Potential

# To derive an equation for the chemical potential in terms of the partition function we start be writing the differential form of the Helmholtz free energy for a single component system
# \begin{eqnarray}
# dA &=& \left(\frac{\partial A}{\partial T} \right)_{N,V}dT + \left(\frac{\partial A}{\partial V} \right)_{N,T}dV + \left(\frac{\partial A}{\partial N} \right)_{T,V}dN \\
# &=& -SdT - PdV + \mu dN 
# \end{eqnarray}
# where $\mu$ is the chemical potential of the species and $N$ is the number of molecules in the system.
# 
# The chemical potential can also be expressed as partial derivatives of other energy functions,
# \begin{equation}
# \mu = \left(\frac{\partial A}{\partial N} \right)_{T,V} = \left(\frac{\partial G}{\partial N} \right)_{T,P} = \left(\frac{\partial H}{\partial N} \right)_{S,P} = \left(\frac{\partial U}{\partial N} \right)_{S,V},
# \end{equation}
# as long as the appropriate natural variables are kept constant.
# 
# It is more convenient to work in terms of moles in which case we convert $N$ to $n$ in the equation for $dA$ and express the molar chemical potential as
# \begin{equation}
# \mu = \left(\frac{\partial A}{\partial n} \right)_{T,V} = \left(\frac{\partial G}{\partial n} \right)_{T,P} = \left(\frac{\partial H}{\partial n} \right)_{S,P} = \left(\frac{\partial U}{\partial n} \right)_{S,V}.
# \end{equation}
# The units of chemical potential written in this way are energy per mole.

# ## Chemical Potential from the Partition Function

# Since the chemical potential can be determined from any energy function and we have a simple expression for $A$ in terms of $Q$, we will derive the formula for chemical potential from $A$.
# 
# We start with the definition of the molar chemical potential in terms of $A$
# \begin{eqnarray}
# \mu &=& \left(\frac{\partial A}{\partial n} \right)_{T,V} \\
# &=& \frac{\partial}{\partial n} \left(-k_BT\ln Q \right)_{T,V} \\
# &=& -k_BT\left(\frac{\partial \ln Q}{\partial n} \right)_{T,V} \\
# &=& -RT\left(\frac{\partial \ln Q}{\partial N} \right)_{T,V}
# \end{eqnarray}

# ## Chemical Potential of an Ideal Gas

# Recall that, for an ideal gas, we can write 
# \begin{equation}
# Q = \frac{q^N}{N!}
# \end{equation}
# To compute the chemical potential, we must differentiate $\ln Q$ with respect to $N$.  We express $\ln Q$ as
# \begin{equation}
# \ln Q = \ln\left(\frac{q^N}{N!} \right) = N\ln q - \ln N! = N\ln q - N\ln N + N
# \end{equation}
# where the last equality employed the use of the Stirling approximation.
# 
# Thus, the chemical potential of an ideal gas is
# \begin{equation}
# \mu = -RT(\ln q - \ln N -1 + 1) = -RT\ln\frac{q}{N}
# \end{equation}
# 
# Recall that $q$ is a function of $V$ and $T$ and that $q \propto V$.  Thus $q/V$ is independent of $V$ and only a function of $T$. This allows us to rewrite $\mu$ as
# \begin{eqnarray}
# \mu = -RT\ln\frac{q}{N} &=& -RT\ln\left[ \frac{q}{V}\frac{V}{N} \right] \\
# &=& -RT\ln\left[\frac{q}{V}\frac{k_BT}{P} \right] \\
# &=& -RT\ln\frac{qk_BT}{V} + RT\ln P
# \end{eqnarray}
# where we used that $\frac{V}{N} = \frac{k_BT}{P}$ from the ideal gas law.
# 
# From the above equation, we call the first term the standard state chemical potential, $\mu^\circ$, 
# \begin{equation}
# \mu^\circ(T) = -RT\ln\left[ \frac{qk_BT}{VP^\circ} \right],
# \end{equation}
# where $P^\circ$ is the standard state pressure of 1 bar or $10^5$ Pa.

# ## Example: Chemical Potential of a Monatomic Ideal Gas

# Compute the standard state chemical potential of Ar(g) at 298.15 K.  Assume ideal gas behavior.

# For a monatomic ideal gas we have that
# \begin{equation}
# q = \left(\frac{2\pi m k_B T}{h^2}\right)^{3/2}V
# \end{equation}
# or 
# \begin{equation}
# \frac{q}{V} = \left(\frac{2\pi m k_B T}{h^2}\right)^{3/2}
# \end{equation}
# 
# Given the temperature and mass of the substance, we can compute $q/V$:
# \begin{eqnarray}
# \frac{q}{V} &=& \left(\frac{2\pi\cdot0.03995\cdot 1.3806\times10^{-23}\cdot298.15}{6.022\times10^{23}(6.626\times10^{-34})^2}\right)^{3/2} \\
# &=& 2.444\times10^{32} \quad\text{m}^{-3}
# \end{eqnarray}
# and $\frac{k_BT}{P^\circ}$:
# \begin{eqnarray}
# \frac{k_BT}{P^\circ} &=& \frac{1.3806\times10^{-23}\cdot 298.15}{1.00\times10^5}  \\
# &=& 4.116\times10^{-26}\quad\text{m}^3
# \end{eqnarray}
# 
# Plugging this into equation above yields
# \begin{eqnarray}
# \mu^\circ(298.15) &=& -8.314\cdot298.15 \cdot \ln(2.444\times10^{32}\cdot 4.116\times10^{-26}) \\
# &=& -39.97 \quad\text{kJ}\cdot\text{mol}^{-1}
# \end{eqnarray}
