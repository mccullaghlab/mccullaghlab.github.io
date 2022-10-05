#!/usr/bin/env python
# coding: utf-8

# # Boltzmann Factor

# ## Learning Goals

# After going through these notes, you should be able to:
# 
# 1. Write out the Boltzmann factor for a given state of a system,
# 2. Define all of the terms/variables in a Boltzmann factor,
# 3. Determine the relative probability of two or more states given distributions and energies,
# 4. Determine the partition function for a system given complete set of states.

# ## Coding Concepts
# 
# The following coding concepts are used in this notebook:
# 1. [Variables](../../coding_concepts/variables.ipynb)
# 2. [Functions](../../coding_concepts/functions.ipynb)
# 3. [Plotting with matplotlib](../../coding_concepts/plotting_with_matplotlib.ipynb)

# ## Introduction and Motivation

# Our goal for learning Thermodynamics is to be able to predict if a process (e.g. chemical reaction) occurs under certain conditions.  The tools/laws of Thermodynamics will help us do this but classical Thermodynamics deals with macroscopic quantities (e.g. temperature, heat, or work for a mole of a substance).  As chemists, we are more used to thinking about chemicals and chemical reactions on the molecular level so first we must come up with a way to going from molecular properties to macroscopic properties.  This is the field of statistical Thermodynamics and the Boltzmann factor is one of the most fundamental concepts in this discipline.  

# ## The Energy of a Macroscopic System

# Consider a macroscopic system such a liter of gas with a fixed number of particles, $N$, and volume, $V$.  Even though there may be almost a mole of molecules of this gas, we can still consider determining the energy of the system.  This can be done either using the Schrodinger equation if we want to model the system quantum mechanically,
# \begin{equation}
# \hat{H}_N \Psi_j = E_j \Psi_j
# \end{equation}
# where $j$ denotes the energy level and the energy is a function of $N$ and $V$ ($E(N,V)$), or we can define a classic $N$-body Hamiltoninan/energy function
# \begin{equation}
# H(R^N, p^N) = \sum_{i}^N \frac{p_i^2}{2m_i} + \sum_{i}^N U_1(R_i) + \sum_{i,j>i}^N U_2(R_{i},R_j) + ...
# \end{equation}
# where $R^N$ are the 3D coordinates of the $N$ particles and $p^N$ are the momenta of the $N$ particles.
# 
# In the special case that the gas in an ideal gas, the energy of the system will simply be a sum of individual particle energies.  Quantum mechanically that would be
# \begin{equation}
# E_j(N,V) = \sum_{i}^N \epsilon_i,
# \end{equation}
# where $\epsilon_i$ is the energy of particle $i$.  If the gas is monatomic (no internal potential energy terms) then $\epsilon_i$ is solution to the particle in the 3D box
# \begin{equation}
# \epsilon_{n_x,n_y,n_z} = \frac{h^2}{8ma^2}\left(n_x^2+n_y^2+n_z^2\right),
# \end{equation}
# where $n_x, n_y, n_z = 1, 2, ...$ are the quantum numbers in the three directions, $h$ is Planck's constant, $m$ is the mass of the particle, and $a = V^{1/3}$ is related to the volume of the box.
# 
# 
# In classical mechanics, the energy of a monatomic ideal gas is simply the kinetic energy.  Each particle has a kinetic energy and that is it, thus the classical Hamiltonian is
# \begin{equation}
# H(R^N,p^N) = \sum_{i}^N \frac{p_i^2}{2m_i}
# \end{equation}

# ## An Ensemble of Systems

# We now consider the question, what is the probability that my system has a particular energy?  Quantum mechanically, this would be what is the probability of observing a particular value for $E_j(N,V)$.  Classically, it would be what the probability of observing a particular value of $H(R^N)$.  To address this question, we consider an ensemble of identical systems.  Each of these systems is in contact with the same thermal bath and thus it can be considered that each system maintains a constant $N$, $V$, and $T$.  
# 
# We start by asking a related but slighly simpler question: what is the probability that a system takes on an energy $E_1$ compared to probability the system takes on an energy $E_2$?  We denote the total number of systems with energy $E_1$ as $a_1$ and the total number of systems with energy $E_2$ as $a_2$.  The relative number of systems in these two states is
# \begin{equation}
# \frac{a_2}{a_1} = f(E_1,E_2) = f(E_1-E_2),
# \end{equation}
# where $f$ is some function of both energies and the last equality must hold because only differences in energies matter.
# 
# We can also compare a state 3 to states 1 and 2 to get:
# \begin{eqnarray}
# \frac{a_3}{a_1} &=& f(E_1-E_3) \\
# \frac{a_3}{a_2} &=& f(E_2-E_3).
# \end{eqnarray}
# We can combine these results with above to demonstrate
# \begin{equation}
# \frac{a_3}{a_1} = \frac{a_2}{a_1}\cdot\frac{a_3}{a_2},
# \end{equation}
# or 
# \begin{equation}
# f(E_1-E_3) = f(E_1-E_2)f(E_2-E_3).
# \end{equation}
# 
# We see that multiplication of $f$ functions leads to addition of their arguments ($E_1-E_2 + E_2-E_3 = E_1-E_3$).  A function for which this is true is the exponential function. Thus, we can see that
# \begin{equation}
# f(E) = e^{\beta E},
# \end{equation}
# where $\beta$ is an arbitrary constant.  And
# \begin{equation}
# a_j = Ce^{-\beta E_j},
# \end{equation}
# where $C$ is an arbitrary constant.
# 
# In this form, $a_j$ is known as the Boltzmann factor for state $j$.

# ## Partition Function

# The Boltzmann factor given above has two components that need to be defined: $\beta$ and $C$.  Presented without derivation,
# \begin{equation}
# \beta = \frac{1}{k_BT},
# \end{equation}
# where $k_B$ is the Boltzmann constant and $T$ is the temperature in Kelvin.  
# 
# $C$ is determined by summing (integrating in classical mechanics) over all possible states
# \begin{equation}
# \sum_i a_i = C\sum_i e^{-\beta E_i}.
# \end{equation}
# The sum of all states/systems must be equal to the total number of systems, $A$, thus we have
# \begin{eqnarray}
# A &=& C\sum_i e^{-\beta E_i} \\
# \Rightarrow C &=& \frac{A}{\sum_i e^{-\beta E_i}}
# \end{eqnarray}
# 
# Plugging this back into the equation for $a_j$ we get
# \begin{equation}
# \frac{a_j}{A} = \frac{e^{-\beta E_j}}{\sum_i e^{-\beta E_i}}
# \end{equation}
# 
# $\frac{a_j}{A}$ is the fraction of systems with energy $E_j$ and thus we recognize this as a probability, or written as 
# \begin{equation}
# P_j = \frac{e^{-\beta E_j}}{\sum_i e^{-\beta E_i}}
# \end{equation}
# where $P_j$ is now the probability of observing a system with energy $E_j$.  

# The denominator of the probability equation, $\sum_i e^{-\beta E_i}$, is referred to as the Canoncial partition function and is denoted $Q$.  This function is depdent on $N$, $V$ and $T$ and is often written as 
# \begin{equation}
# Q(N,V,T) = \sum_i e^{-\frac{E_i}{k_B T}},
# \end{equation}
# or 
# \begin{equation}
# Q(N,V,\beta) = \sum_i e^{-\beta E_i},
# \end{equation}

# ## Example: Quantized States

# Consider a single monoatomic ideal gas particle with the mass of helium.  The system is kept at a temperature that allows it to populate two states: $(n_x, n_y, n_z) = (1,1,1)$ and $(n_x, n_y, n_z) = (2,1,1)$. Compute the Boltzmann factors and relative probability of these two states.  Leave your answers in terms of $\beta$, $a$, and $m$.

# We start by writing out the energy of the system from above:
# 
# \begin{equation}
# E_j(N,V) = \epsilon_i,
# \end{equation}
# 
# where 
# 
# \begin{equation}
# \epsilon_{n_x,n_y,n_z} = \frac{h^2}{8ma^2}\left(n_x^2+n_y^2+n_z^2\right).
# \end{equation}
# 
# Plugging in the two allowed sets of quantum numbers we get two possible total energies of the system:
# 
# \begin{eqnarray}
# E_1 &=& \frac{3h^2}{8ma^2} \\
# E_2 &=& \frac{6h^2}{8ma^2}
# \end{eqnarray}
# 
# The partition function is:
# \begin{eqnarray}
# Q &=& e^{-\beta E_1} + e^{-\beta E_2} \\
# &=& e^{-\beta \frac{3h^2}{8ma^2}} + e^{-\beta\frac{6h^2}{8ma^2}} \\
# &=& e^{-\beta \frac{3h^2}{8ma^2}} \left(e + e^{2} \right)
# \end{eqnarray}
# 
# The Boltzmann factors are:
# \begin{eqnarray}
# e^{-\beta E_1} &=& e^{-\beta \frac{3h^2}{8ma^2}} \\
# e^{-\beta E_2} &=& e^{-\beta \frac{6h^2}{8ma^2}}
# \end{eqnarray}
# 
# Relative probability is:
# \begin{eqnarray}
# \frac{a_1}{a_2} &=& \frac{e^{-\beta E_1}}{e^{-\beta E_2}} \\
# &=& \frac{e^{-\beta \frac{3h^2}{8ma^2}}}{e^{-\beta \frac{6h^2}{8ma^2}}} \\
# &=& e^{\beta \frac{3h^2}{8ma^2}}
# \end{eqnarray}

# We see that the ratio of the first two translational states is an exponential function that depends on temperature, mass of the particle, and size of the box.  If we consider the particle to have the mass of a Helium atom (4.002 amu) and to be constrained in a micrometer length box, we get the following ratio of these as a function of temperature. 

# In[1]:


# import numpy library for use later
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
# define a function of a1/a2
k = 1.38e-23
h = 6.626e-34
def a1_a2(T,m,a):
    return np.exp(-3*h**2/(k*T*8*m*a**2))

# define variables from the problem
a = 1e-6             # in m
m = 4.002*1.66e-27   # in kg
# define temperature domain
T = np.arange(1.0,100,0.1)
# make a plot
plt.plot(T,a1_a2(T,m,a))


# In[ ]:




