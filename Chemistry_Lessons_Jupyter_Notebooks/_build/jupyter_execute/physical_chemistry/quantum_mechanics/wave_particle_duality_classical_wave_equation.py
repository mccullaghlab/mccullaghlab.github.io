#!/usr/bin/env python
# coding: utf-8

# # Wave Particle Duality and the Classical Wave Equation

# ## Motivation
# 
# We have already used the concept that light and matter both exhibit aspects of wave-like and particle-like behavior.  de Broglie formalized this hypothesis in his famous wavelength equation.  Additionally, we will discuss how waves are treated classically.

# ## Learning Goals
# 
# After working through these notes, you will be able to:
# 
# 1. Use the de Broglie wavelength equation to convert between wavelength and particle mass
# 2. Solve the classical wave equation using the methods of separation of variables

# ## Coding Concepts
# 
# The following coding concepts are used in this notebook:
# 
# 1. [Variables](../../coding_concepts/variables.ipynb)
# 2. [Plotting with matplotliub](../../coding_concepts/plotting_with_matplotlib.ipynb)

# ## Matter has wavelike properties

# Both light and matter have wave-like and particle-like properties.  This is more readily observed in light.  White light splitting into colors of the visible spectrum when shone through a prism is an example of light displaying wave-like character.  The photoelectric effect, where a photon of light with energy $h\nu$ collides with an electron causing it to leave the metal if the incident energy is larger than the work function, is an example of light displaying matter-like character.  In 1925, de Broglie reasoned that both light and matter display this behavior.  Whats more, the realatonship between wavelength and momentum derived for light by Einstein also applies to matter
# \begin{equation}
# \lambda = \frac{h}{p} = \frac{h}{mv}    
# \end{equation}

# ### Example: The de Broglie Wavelength of Matter
# 
# Compute the de Broglie Wavelength of the following two objects:
# 
# 1. An electron traveling at $3\times10^{5}$ m$\cdot$s$^{-1}$
# 2. A baseball (5.0 oz) traveling at 90 mph

# ### Solution: The de Broglie Wavelength of Matter
# 
# 1. An electron traveling at $3\times10^5$ m$\cdot$s$^{-1}$.  For this we compute the momentum of the electron:
# 
# \begin{align}
# p = m_ev &= (9.109\times10^{-31}\text{ kg})(3\times10^5 \text{ m}\cdot\text{s}^{-1}) \\
# &= 2.73 \times 10^{-25} \text{ kg}\cdot\text{m}\cdot\text{s}^{-1}
# \end{align}
# 
# The de Broglie wavelength is then
# \begin{align}
# \lambda &= \frac{6.626\times10^{-34} \text{ J}\cdot\text{s}}{2.73 \times 10^{-25} \text{ kg}\cdot\text{m}\cdot\text{s}^{-1}} \\
# &= 2.43 \text{ nm}
# \end{align}
# 
# 2. For a baseball traveling at 90 mph
# 
# \begin{align}
# m = 5\text{ oz}\frac{1\text{ lb}}{16 \text{ oz}} \frac{0.454 \text{kg}}{1 \text{ lb}} = 0.14 \text {kg}
# \end{align}
# 
# \begin{align}
# p = m_ev &= (0.14\text{ kg})\left(90\text{ mph} \frac{1610\text{ m}}{1 \text{mi}} \frac{1 \text{hr}}{3600 \text{s}}\right) \\
# &= 5.6 \text{ kg}\cdot\text{m}\cdot\text{s}^{-1}
# \end{align}
# 
# \begin{align}
# \lambda = \frac{6.626\times10^{-34}}{5.6} = 1.2\times10^{-34} \text{ m}
# \end{align}

# In[5]:


print(9.109e-31*3e5)
print(6.626e-34/(9.109e-31*3e5))


# From this example we see that, at this velocity, the wavelength of electron motion is of a similar magnitude to the radius of hydrogen atom orbits.  This demonstrates that the wavelike character of the electron is important for its behavior.  
# 
# The wavelength of the baseball, however, is miniscule, even for atomic distances.  Thus, the wavelike behavior of a baseball is irrelevant for its behavior.  
# 
# The wavelike behavior of electrons is not only important for determining the structure and behavior of atoms and molecules.  It is utilized in electron microscopes and can display similar wavelike defraction patterns to X-rays.  

# ## The 1-D Wave Equation

# The behavior and physics of waves might be foreign to us.  Waves are, however, readily treated using classical mechanics.  Consider, for example, a string of tied between two objects a length $l$ apart.  We consider the motion of the string (i.e. strumming a guitar string). The amplitude of the motion normal to the direction of the wave is described by a time dependent function $u(x,t)$ where $x$ is the position along the end separation vector and $t$ is time.  The classical wave equation relates the second derivative of the function $u$ in position and in time
# \begin{equation}
# \frac{\partial^2 u(x,t)}{\partial x^2} = \frac{1}{v^2}\frac{\partial^2u(x,t)}{\partial t^2},
# \end{equation}
# where $v$ is the linear velocity of the wave.  This is a linear partial differential equation: linear because the second derivatives only appear to the first power and there are not cross partial derivatives and partial because derivatives of the same function w.r.t different variables show up in the same equation.  
# 
# The above equation is subject to the following boundary conditions:
# \begin{align}
# u(0,t) =& 0\\
# u(l,t) =& 0
# \end{align}
# Since the string is anchored at both ends.

# In[20]:


import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
# make an array containing domain of wavelengths to consider
x = np.arange(0,np.pi,0.0001)
# setup plot parameters
fig = plt.figure(figsize=(10,5), dpi= 80, facecolor='w', edgecolor='k')
ax = plt.subplot(111)
ax.set_xlabel(r'$x$',size=20)
plt.tick_params(axis='x',labelsize=20)
plt.ylim((0.0,1.05))
# plot quantum result
ax.plot(x,np.sin(x),lw=3)
plt.xlim(0,np.pi)
# ticks
plt.xticks([0,np.pi],[0,'l'])
plt.yticks(ticks=None);


# ## Solving the 1-D Wave Equation with Separation of Variables

# The first major step to solve the above differential equation is to assume that the dependence of the two different variables can be separated.  Namely that
# \begin{equation}
# u(x,t) = X(x)T(t)
# \end{equation}
# where $X(x)$ is some function strictly of $x$ and $T(t)$ is another function strictly of $t$.  Substituting this for for $u$ into the differential equation above yields
# \begin{equation}
# T(t)\frac{d^2X(x)}{dx^2} = \frac{1}{v^2}X(x)\frac{d^2T(t)}{dt}
# \end{equation}
# Dividing both sides of the above equation by $u(x,t) = X(x)T(t)$ will complete the separation of variables
# \begin{equation}
# \frac{1}{X(x)}\frac{d^2X(x)}{dx^2} = \frac{1}{v^2T(t)}\frac{d^2T(t)}{dt}
# \end{equation}
# 
# Since the left-hand side depends only on $x$ and the right-hand side depends only on $t$, and these two coordinates are independent, the only way for this to be achieved is for both sides to be equal to the same constant.  That is
# \begin{align}
# \frac{1}{X(x)}\frac{d^2X(x)}{dx^2} &= K \\
# \frac{1}{v^2T(t)}\frac{d^2T(t)}{dt} &= K
# \end{align}
# 
# Each of these equations can be rewritten as
# \begin{align}
# \frac{d^2X(x)}{dx^2} - KX(x) = 0 \\
# \frac{d^2T(t)}{dt} - v^2KT(t) = 0
# \end{align}
# 
# These are ordinary differential equations that can be solved using standard ODE procedures.  We now consider the solutions to these equations for three scenarios
# \begin{align}
# K &= 0\\
# K &>0\\
# K &< 0
# \end{align}

# ### K =0 

# $K = 0$ is what is called the *trivial* solution.  Namely this solution is such that
# \begin{equation}
# u(x,t) = 0
# \end{equation}
# 
# This is a mathematically valid solution to our problem but of no physical interest.

# ### K > 0

# If $K>0$, we get the following general form for $X(x)$ and $T(t)$:
# \begin{align}
# \frac{d^2y(x)}{dx^2} - k^2y(x) = 0 
# \end{align}
# 
# Or you can write this as 
# \begin{align}
# \frac{d^2y(x)}{dx^2} = k^2y(x) 
# \end{align}
# which implies that $y$ is a function that, when you take two derivatives, returns itself times a constant squared.  This type of behavior is exhibited by functions of the form $e^{\alpha x}$ the value of $\alpha$ needs to be determined.  So, we propose that $y(x) = e^{\alpha x}$ and plug this into the above equation to determine the expression for $\alpha$.
# 
# \begin{align}
# &\frac{d^2e^{\alpha x}}{dx^2} - k^2e^{\alpha x} = 0 \\
# & \alpha^2e^{\alpha x} - k^2e^{\alpha x} = 0 \\
# & (\alpha^2 - k^2)e^{\alpha x} = 0 \\
# \Rightarrow & (\alpha^2 - k^2) = 0 \\
# \Rightarrow & \alpha = \pm k
# \end{align}
# 
# This means that there are two possible values for $\alpha$ and thus two possible solutions, $y(x)$.  Namely
# \begin{align}
# y(x) = e^{kx} \\
# y(x) = e^{-kx}
# \end{align}
# 
# In general, when there are two possible solutions for a given differential equation, then any linear combination of these two solutions is also a solution.  Thus, it is more complete to write that the solution is
# \begin{equation}
# y(x) = c_1e^{kx} + c_2e^{-kx},
# \end{equation}
# where $c_1$ and $c_2$ are constants.

# ### K < 0

# If $K<0$, we can set $K = (ki)^2$ for arbitrary constant $k>0$.  This implies the general form for $X(x)$ and $T(t)$:
# \begin{align}
# \frac{d^2y(x)}{dx^2} - (ki)^2y(x) = 0 
# \end{align}
# 
# Or you can write this as 
# \begin{align}
# \frac{d^2y(x)}{dx^2} = (ki)^2y(x) 
# \end{align}
# which implies that $y$ is a function that, when you take two derivatives, returns itself times a constant squared.  This type of behavior is exhibited by functions of the form $e^{\alpha x}$ the value of $\alpha$ needs to be determined.  So, we propose that $y(x) = e^{\alpha x}$ and plug this into the above equation to determine the expression for $\alpha$.
# 
# \begin{align}
# &\frac{d^2e^{\alpha x}}{dx^2} - (ik)^2e^{\alpha x} = 0 \\
# & \alpha^2e^{\alpha x} - (ki)^2e^{\alpha x} = 0 \\
# & (\alpha^2 - (ki)^2)e^{\alpha x} = 0 \\
# \Rightarrow & (\alpha^2 - (ki)^2) = 0 \\
# \Rightarrow & \alpha = \pm ki
# \end{align}
# 
# This means that there are two possible values for $\alpha$ and thus two possible solutions, $y(x)$.  Namely
# \begin{align}
# y(x) = e^{ikx} \\
# y(x) = e^{-ikx}
# \end{align}
# 
# In general, when there are two possible solutions for a given differential equation, then any linear combination of these two solutions is also a solution.  Thus, it is more complete to write that the solution is
# \begin{equation}
# y(x) = c_1e^{ikx} + c_2e^{-ikx},
# \end{equation}
# where $c_1$ and $c_2$ are constants.
# 
# In this case, we can also recall/recognize that 
# \begin{equation}
# e^{ix} = \cos(x) + i \sin(x)
# \end{equation}
# 
# Using this relationship is can be readily shown that
# \begin{equation}
# y(x) = (c_1 + c_2)\cos(kx) + i(c_1-c_2)\sin(kx)
# \end{equation}
