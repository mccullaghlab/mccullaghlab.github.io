#!/usr/bin/env python
# coding: utf-8

# # Classical Ideal Gas

# ## Motivation
# 
# So far we have considered monatomic and diatomic ideal gasses composed of quantum particles.  The energy levels of the various aspects of the particles is quantized (though we introduced some continuous approximations).  Here we ask what happens when we consider classical particles?  How does this affect the partition function for a system?

# ## Learning Goals
# After these notes, students should be able to:
# 1. Describe the difference between the variables that dictate the state of a quantum system and a classical system
# 2. Write out the classical energy for a system of a monatomic ideal gas
# 3. Discuss the units of a p.f.

# ## Coding Concepts
# 
# 1. Variables
# 2. Functions
# 3. Numpy

# ## Introduction to Classical Particles
# 
# If we consider a gas to be made up of classical particles, we can no longer describe its energy using the quantized energy functions.  Instead, we must consider what it means for a particle to be "classical" and how we describe the energy of such particles.  
# 
# A classical particle will have a position, $\vec{r} = (r_x, r_y, r_z)$, and a momentum, $\vec{p} = (p_x, p_y, p_z) = m\vec{v}$.  If the particle is in a box with linear dimension $L$ (cubic box such that $L^3=V$) then each component of $\vec{r}$ can take on values of $0 < r_x,r_y,r_z < L$.  The components of the momentum of the particle are not limited in the same manner.  Each component of the momentum can take on values of $-\infty < p_x,p_y,p_z < \infty$.
# 
# 
# In a system of $N$ classical particles, we can denote the collection of particle positions as $\vec{R}^N$ and particle momenta as $\vec{P}^N$.  Each of these collections is made up of $N$ different three dimensional vectors denoting the position or momentum of each particle.  
# 
# The energy of a system of classical particles may depend on the particular configuration of the particles, $\vec{R}^N$, and values of the momenta, $\vec{P}^N$.  We might write that the energy, denoted $H$, is
# \begin{equation}
# H(\vec{R}^N,\vec{P}^N) = K(\vec{P}^N) + U(\vec{R}^N),
# \end{equation}
# where $K(\vec{P}^N)$ is the kinetic energy and $U(\vec{R}^N)$ is the potential energy.

# ## Example 1: Position and Potential Energy with Classical Particles
# 
# Consider a system of two classical particles that interact with a Lennard-Jones potential of the form
# \begin{equation}
# U(\vec{r}^1,\vec{r}^2) = 4\left[ \left( \frac{2}{|\vec{r}^{12}|}\right)^{12} -  \left( \frac{2}{|\vec{r}^{12}|}\right)^{6} \right]
# \end{equation}
# where $|\vec{r}^{12}| = \sqrt{(r^1_x - r^2_x)^2 + (r^1_y - r^2_y)^2 + (r^1_z - r^2_z)^2}$ is the distance between particles 1 and 2.  Consider this potential to be in units of $\epsilon$.  Compute the potential energy when the particle have postitions $\vec{r}^1 = (0.25, 0.6, 1.0)$ and $\vec{r}^2 = (3.0, 1.2, 0.8)$

# In[5]:


import numpy as np
def lj(r1,r2,eps,sigma):
    r12 = np.sqrt((r1[0]-r2[0])**2+(r1[1]-r2[1])**2+(r1[2]-r2[2])**2)
    return 4*eps* ((sigma/r12)**12 - (sigma/r12)**6)
r1 = np.array([0.25,0.6,1.0])
r2 = np.array([3.0,1.2,0.8])
print("Potential energy (units of epsilon):", np.round(lj(r1,r2,2,2),2))


# ## Partition Function for a Classical Monatomic Ideal Gas
# 
# If we consider a collection of ideal particles to be an ideal gas then there is no potential energy between the particles.  That is, $U(\vec{R}^N) = 0$ for all configurations of the particles.  This allows us to write that 
# \begin{equation}
# H(\vec{R}^N,\vec{P}^N) = K(\vec{P}^N).
# \end{equation}
# Furthermore, the kinetic energy of a classical particle is simply $\frac{\vec{p}^2}{2m} = \frac{1}{2}m\vec{v}^2$.  The kinetic energy of a collection of classical particles is just the sum over individual particle kinetic energies
# \begin{equation}
# H(\vec{R}^N,\vec{P}^N) = K(\vec{P}^N) = \sum_{i=1}^N \frac{\vec{p}_i^2}{2m_i}.
# \end{equation}
# 
# Note that $\vec{p}^2 = p_x^2 + p_y^2 + p_z^2$ is the squared norm of the momentum.  Also note that if the mass is the ame for all particles it can be factored out of the sum.
# 
# If we wish to estimate the partition function for this ideal gas, we must estimate the Boltzmann factor for all possible configurations for the system.  Since phase space, positions and momenta of the particles, is continuous,  we must consider an integral (rather than a sum) over all of these variables.  That is
# \begin{eqnarray}
# Q = \int e^{-\beta H(\vec{R}^N,\vec{P}^N)} d\vec{R}^Nd\vec{P}^N,
# \end{eqnarray}
# where $d\vec{R}^N = dr_{x1}dr_{y1}dr_{z1}...dr_{xN}dr_{yN}dr_{zN}$ and $d\vec{P}^N = dp_{x1}dp_{y1}dp_{z1}...dp_{xN}dp_{yN}dp_{zN}$.  So this is a $6N$ dimensional integral.  
# 
# Now we plug in the equation above for $H(\vec{R}^N,\vec{P}^N)$ to get
# \begin{eqnarray}
# Q &=& \int e^{-\beta \sum_{i=1}^N \frac{\vec{p}_i^2}{2m_i} } d\vec{R}^Nd\vec{P}^N,\\
#  &=& \int d\vec{R}^N\int e^{-\beta \sum_{i=1}^N \frac{\vec{p}_i^2}{2m_i} } d\vec{P}^N
# \end{eqnarray}
# where the last equality holds because the integrand, $e^{-\beta \sum_{i=1}^N \frac{\vec{p}_i^2}{2m_i} }$, does not depend on any of the particle positions and thus can be pulled out of those integrals.  We now must talk about bounds of integration for the $3N$ positional integrals.  Since we consider this to be a cubic box, all positions are resricted to be in the domain $0 < r_x < L$ thus yielding
# 
# \begin{eqnarray}
# \int d\vec{R}^N &=& \int_0^L dr_{x1} \int_0^L dr_{y1} \int_0^Ldr_{z1}...\int_0^Ldr_{xN}\int_0^Ldr_{yN}\int_0^Ldr_{zN} \\
# &=& L\cdot L\cdot L\cdot ... L\cdot L\cdot L\cdot \\
# &=& V\cdot ... \cdot V \\
# &=& V^{N}
# \end{eqnarray}
# 
# Each one of these integrals is simply equal to $L = \int_0^L dr_{x1}$ and $L^3 = V$.  The integral for $Q$ is now
# \begin{eqnarray}
# Q &=& V^N\int e^{-\beta \sum_{i=1}^N \frac{\vec{p}_i^2}{2m_i} } d\vec{P}^N \\
# &=& V^N \int \prod_{i=1}^N e^{-\beta \frac{\vec{p}_i^2}{2m_i} } d\vec{P}^N \\
# &=& V^N \int\int\int e^{-\beta \frac{\vec{p}_1^2}{2m_1} } dp_{x1}dp_{y1}dp_{z1}\int\int\int e^{-\beta \frac{\vec{p}_2^2}{2m_2} } dp_{x2}dp_{y2}dp_{z2}...\int\int\int e^{-\beta \frac{\vec{p}_N^2}{2m_N} } dp_{xN}dp_{yN}dp_{zN} \\
# &=& V^N \left( \int\int\int e^{-\beta \frac{\vec{p}^2}{2m} } dp_{x}dp_{y}dp_{z}\right)^N\cdot\frac{1}{N!}
# \end{eqnarray}
# where the last equality holds because we consider all particles to be the same mass and indistinguishable.  
# 
# We can simplify this even further if we consider that $\vec{p}^2 = p_x^2 + p_y^2 + p_z^2$ to get that
# \begin{eqnarray}
# \int\int\int e^{-\beta \frac{\vec{p}^2}{2m} } dp_{x}dp_{y}dp_{z} &=& \int\int\int e^{-\beta \frac{p_x^2}{2m} }e^{-\beta \frac{p_y^2}{2m} }e^{-\beta \frac{p_z^2}{2m} }dp_{x}dp_{y}dp_{z} \\
# &=& \int_{-\infty}^{\infty} e^{-\beta \frac{p_x^2}{2m}}dp_x\int_{-\infty}^{\infty} e^{-\beta \frac{p_y^2}{2m} }dp_y\int_{-\infty}^{\infty} e^{-\beta \frac{p_z^2}{2m} }dp_{z} \\
# &=& \left( \int_{-\infty}^{\infty} e^{-\beta \frac{p^2}{2m}}dp \right)^3,
# \end{eqnarray}
# where the last equality holds because there is no reason to consider the momentum in the $x$ direction to be different than that of the $y$ or $z$ directions.
# 
# Finally, we plug this back into the equation for $Q$ to get
# \begin{eqnarray}
# Q &=& V^N \left( \int_{-\infty}^{\infty} e^{-\beta \frac{p^2}{2m}}dp \right)^{3N}\cdot\frac{1}{N!}
# \end{eqnarray}
# 
# The solution to $\int_{-\infty}^{\infty} e^{-\beta \frac{p^2}{2m}}dp$ is left for homework (problem 1a).  Ultimately you will get that
# \begin{eqnarray}
# Q &=& \frac{V^N \left( 2\pi m k_BT \right)^{3N/2}}{N!}
# \end{eqnarray}

# ## A Note on Units of $Q$
# 
# You will notice above that the classical p.f. has units of $V$^N $p$^{3N}.  The quantum p.f. is unitless.  It is typical to introduce $h^{-1}$ in each direction for classical p.f.s so that it is unitless and there is greater correspondance between classical and quantum.
