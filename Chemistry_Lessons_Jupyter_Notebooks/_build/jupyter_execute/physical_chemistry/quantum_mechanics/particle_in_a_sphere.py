#!/usr/bin/env python
# coding: utf-8

# # Particle in a Sphere

# ## Motivation:
# 
# Gearing up to work on the hydrogen atom, we recognize that the hydrogen atom is a spherically symmetric system.  That is, we care about the distance of the electron from the proton/nucleus not the specific $x$, $y$, and $z$ coordinates of the electron.  Here, we consider what constraining a particle to a sphere does as opposed to a cube.

# ## Learning Goals:
# 
# After working through these notes, you will be able to:
# 
# 1. Convert between Cartesian and spherical polar coordinates
# 2. Write out the Laplacian in spherical polar coordinates
# 3. Write out the Hamiltonian for the particle in a the sphere
# 4. Perform separation of variables on the particle in a sphere Hamiltonian

# ## Coding Concepts:
# 
# The following coding concepts are used in this notebook
# 
# 1. [Variables](../../coding_concepts/variables.ipynb)
# 2. [Functions](../../coding_concepts/functions.ipynb)
# 3. [Plotting with matplotlib](../../coding_concepts/plotting_with_matplotlib.ipynb)

# ## Statement of the Problem

# Consider a particle restricted to be ***in a sphere***.  This is equivalent to a free particle (no potential) inside a fixed radius and an infinte potential at $r_0$.  
# 
# We would like to know the allowed energies and probability of positions for such a system.  To get these we will setup and solve the Schrodinger equation.
# 
# We start by defining the Hamiltonian for such a system.

# ## The Hamiltonian

# We start constructing the Hamiltonian in the same way for every problem.  The Hamiltonian is always (for the sake of this class/these notes) a sum of the Kinetic Energy operator and the Potential Energy operator:
# \begin{equation}
# \hat{H} = \hat{K} + \hat{V}.
# \end{equation}
# 
# The potential energy for this problem is zero since the particle is *free* inside the sphere and cannot go outside of the sphere.  Thus
# \begin{equation}
# \hat{H} = \hat{K}
# \end{equation}
# 
# The Kinetic Energy operator in three dimensions for a single particle is
# \begin{equation}
# \hat{K} = -\frac{\hbar^2}{2m}\left( \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2} + \frac{\partial^2}{\partial z^2}\right)
# \end{equation}
# 
# We could attempt to solve this problem in Cartesian coordinates.  Solving the differential equation will yield the same initial solutions as the 3D particle in a box:
# \begin{align}
# X(x) = A\sin k_x x + B\cos k_x x \\
# Y(y) = A\sin k_y y + B\cos k_y y \\
# Z(z) = A\sin k_z z + B\cos k_z z 
# \end{align}
# where the total wavefunction is $\psi(x,y,z) = X(x)Y(y)Z(z)$
# 
# The issue is that applying the boundary condition that
# \begin{equation}
# \psi(r=r_0) = 0
# \end{equation}
# is difficult in Cartesian coordinates in part because it requires a coupling of $x$, $y$, and $z$.  Instead, we will solve the problem in spherical polar coordinates.

# ## Spherical Polar Coordinates

# In spherical polar coordinates, we define a point in $\mathbb{R}^3$ by a distane from the origin, $r$, and two angles, $\theta$ and $\phi$.  
# 
# ![title](img/SphericalCoordinates_1201.svg)
# 
# These can be expressed in terms of Cartesian coordinates as
# \begin{align}
# r &= \sqrt{x^2+ y^2 + z^2} \\
# \phi &= \tan^{-1}\frac{y}{x} \\
# \theta &= \cos^{-1}\frac{z}{r}
# \end{align}
# Note that the following domain limits are placed on these coordinates
# \begin{align}
# 0 \leq &r < \infty \\
# 0 \leq &\theta \leq \pi \\
# 0 \leq &\phi \leq 2\pi
# \end{align}
# The Cartesian coordinated in terms of their spherical counterparts are:
# \begin{align}
# x &= r\sin\theta\cos\phi \\
# y &= r\sin\theta\sin\phi \\
# z &= r\cos\theta
# \end{align}
# 
# In order to describe the particle in a sphere problem in spherical coordinates, we need to define the Hamiltonian in spherical coordinates.  This can be done with the equations above and their derivatives and second derivatives.  The derivation is long and tedious so here we present the result for the Laplacian in spherical polar coordinates:
# \begin{equation}
# \nabla^2_{r\theta\phi} = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{r^2\sin^2\theta}\frac{\partial^2}{\partial^2\phi}
# \end{equation}

# ## Separation of Variables in Spherical Coordinates

# The Schrodinger equation in spherical polar coordinates is
# 
# \begin{align}
# \hat{H}\psi(r,\theta,\phi) &= E \psi(r,\theta,\phi) \\
# -\frac{\hbar^2}{2m}\nabla^2_{r\theta\phi} &= E \psi(r,\theta,\phi) \\
# -\frac{\hbar^2}{2m}\left[\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{r^2\sin^2\theta}\frac{\partial^2}{\partial^2\phi}\right] \psi(r,\theta,\phi) &= E \psi(r,\theta,\phi)
# \end{align}
# 
# This equation does not look immediately separable.  Indeed it looks quite intimidating.  It turns out, however, that the equation is separable into $r$, $\theta$, and $\phi$ components. We will start by separating $r$ from $\theta$ and $\phi$.  
# 
# Multiply both sides of the above equation by $2mr^2$ to get
# 
# \begin{equation}
# -\hbar^2\left[\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right) - \frac{\hbar^2}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2}{\partial^2\phi}\right]\psi(r,\theta,\phi)  = 2mr^2E\psi(r,\theta,\phi)
# \end{equation}
# 
# Rearrange (combine terms dependent on $r$) to get
# 
# \begin{equation}
# -\hbar^2\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2}{\partial^2\phi}\right]\psi(r,\theta,\phi) -\hbar^2\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)\psi(r,\theta,\phi) = 2mr^2E\psi(r,\theta,\phi)
# \end{equation}
# 
# Notice that we have a sum of a term that depends on $\theta$ and $\phi$ and then two terms that depend only on $r$.  Thus we can posit a separation of variables such that
# \begin{equation}
# \psi(r,\theta,\phi) = R(r)Y(\theta,\phi)
# \end{equation}
# 
# Additionally, we will define the angular momentum operator (that only operates on $\theta$ and $\phi$) as
# \begin{equation}
# \hat{L}^2 = -\hbar^2\left(\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2}{\partial^2\phi}\right)
# \end{equation}
# 
# Substite the above equation for $\psi$ into the Schrodinger equation to get
# \begin{align}
# \hat{L}^2R(r)Y(\theta,\phi) -\hbar^2 \frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)R(r)Y(\theta,\phi) &= 2mr^2ER(r)Y(\theta,\phi) \\
# \Rightarrow \frac{-\hbar^2}{Y(\theta,\phi)}\hat{L}^2Y(\theta,\phi) -\frac{\hbar^2}{R(r)}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)R(r) &= 2mr^2E
# \end{align}
# 
# Notice that the left-hand term depends on $\theta$ and $\phi$ but the other terms do not.  Thus, we must have that
# \begin{equation}
# \frac{1}{Y(\theta,\phi)}\hat{L}^2Y(\theta,\phi) = E_{\theta,\phi},
# \end{equation}
# where $E_{\theta,\phi}$ is a constant.
# 
# Additionally, we must have that
# \begin{equation}
# \frac{\hbar^2}{R(r)}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)R(r) + 2mr^2E = E_{\theta,\phi}
# \end{equation}
# 
# 
# We will solve these two equations separately.

# ## Solutions to $\theta$ and $\phi$ equation

# It is the goal of this section of the notes to find a solution to the differential equation
# \begin{equation}
# \frac{1}{Y(\theta,\phi)}\hat{L}^2Y(\theta,\phi) = E_{\theta,\phi},
# \end{equation}
# where $E_{\theta,\phi}$ is a constant.
# 
# This equation can be rearranged to the eigenvalue equation
# \begin{equation}
# \hat{L}^2Y(\theta,\phi) = E_{\theta,\phi}Y(\theta,\phi)
# \end{equation}
# Now plugging back in the expanded form of $\hat{L}^2$ we get
# \begin{align}
# -\hbar^2\left(\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2}{\partial^2\phi}\right)Y(\theta,\phi) &= E_{\theta,\phi}Y(\theta,\phi) \\
# \Rightarrow \left(\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2}{\partial^2\phi}\right)Y(\theta,\phi) &= \frac{-E_{\theta,\phi}}{\hbar^2}Y(\theta,\phi)
# \end{align}
# 
# Multiplying the above equation by $\sin^2\theta$ yields
# \begin{align}
# \left(\sin\theta\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\frac{\partial^2}{\partial^2\phi}\right)Y(\theta,\phi) &= \frac{-E_{\theta,\phi}\sin^2\theta}{\hbar^2}Y(\theta,\phi)
# \end{align}
# We notice that this motivates a separation of variables into 
# \begin{equation}
# Y(\theta,\phi) = \Theta(\theta)\Phi(\phi)
# \end{equation}
# 
# Substituting this and $\beta = \frac{E_{\theta,\phi}}{\hbar^2}$ yields
# 
# \begin{equation}
# \left[\sin\theta\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)+\beta\sin^2\theta \right] \Theta(\theta)\Phi(\phi) = -\frac{\partial^2}{\partial^2\phi} \Theta(\theta)\Phi(\phi).
# \end{equation}
# 
# Divide both sides of the equation above by $\Theta(\theta)\Phi(\phi)$ and note that the operator on the right-hand side is independent of $\theta$ to yield
# 
# \begin{equation}
# \frac{\sin\theta}{\Theta(\theta)}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)\Theta(\theta)+\beta\sin^2\theta = -\frac{1}{\Phi(\phi)}\frac{\partial^2}{\partial^2\phi} \Phi(\phi).
# \end{equation}
# 
# Now the left-hand side is independent of $\phi$ and the right-hand side is independent of $\theta$.  Since these two things are equal but independent of the other's variable they must be constant.  We will define this constant as $m^2$ (for reasons that will become clear later) and solve the following two equations independently
# 
# $m^2 = -\frac{1}{\Phi(\phi)}\frac{\partial^2}{\partial^2\phi} \Phi(\phi) \tag{1}$
# 
# $m^2 = \frac{\sin\theta}{\Theta(\theta)}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)\Theta(\theta)+\beta\sin^2\theta \tag{2}$

# ### Solutions to $\phi$ equation (1)

# Equation (1) above can me simply rearranged to give
# 
# $-m^2 \Phi(\phi)= \frac{\partial^2}{\partial^2\phi} \Phi(\phi)$
# 
# which is straightforward eigenvalue-eigenvector problem with solutions
# 
# $\Phi(\phi) = A_me^{im\phi}\quad \mathrm{and}\quad A_{-m}e^{-im\phi}$.
# 
# Applying boundary conditions, $\Phi(\phi+2\pi) = \Phi(\phi)$ yields the quantization
# 
# $m=0,\pm 1, \pm 2, ...$
# 
# Thus we can write
# 
# $\Phi(\phi) = Ae^{im\phi} \quad m=0,\pm 1, \pm 2, ...$.
# 
# Normalization yields $A=\frac{1}{\sqrt{2\pi}}$.

# ### Solutions to $\theta$ equation (2)

# The solutions to equation (2) above are not as straightforward as those for equation (1).  We start by rewriting the original equation here:
# 
# $m^2 = \frac{\sin\theta}{\Theta(\theta)}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right)\Theta(\theta)+\beta\sin^2\theta$.
# 
# Now make a change of variable $x = \cos\theta$ which yields $\frac{dx}{-\sin\theta}=d\theta$ and define $P(x) = \Theta(\theta)$.  Plugging these in an performing some rearrangements yields the Legendre equation
# 
# $(1-x^2)\frac{d^2}{dx^2}P(x)-2x\frac{d}{dx}P(x)+\left[\beta-\frac{m^2}{1-x^2}\right]P(x) = 0$.
# 
# For $\Theta(\theta)$ to be continous $\beta=J(J+1)$ where $J=0,1,2...$.  Note that this also puts a limit on $m$ with $m=0,\pm 1, \pm 2, ... , \pm J$.   The quantization of $\beta$ leads to the quantization of energy
# 
# $E_J = \frac{\hbar^2}{2I}J(J+1)$.
# 
# The solutions, $P(x)$, to the Legendre equation are known as the Associated Legendre polynomials.  
# 
# $P_\nu^m = (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_\nu(x)$
# 
# where
# 
# $P_\nu(x) = \sum_{k=0}^{\infty}\frac{(-\nu)_k(\nu+1)_k}{k!^2}\left(\frac{1-x}{2}\right)^k$
# 
# and $(\nu)_k = \frac{(\nu+k-1)!}{(\nu-1)!}$.

# In[8]:


# plot of some of the Legendre polynomials
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from scipy.special import lpmv
x = np.arange(-1,1,0.001)
plt.figure(figsize=(12,6),dpi= 80, facecolor='w', edgecolor='k')
plt.tick_params(axis='both',labelsize=20)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
for l in range(4):
    for m in range(l):
        label = "l=" + str(l) + ", m=" + str(m)
        plt.plot(x,lpmv(m,l,x),lw=4,label=label)
plt.legend(fontsize=16);


# ### Combining Solutions to $\theta$ and $\phi$

# The total wavefunctions are the product of $\Phi(\phi)$ and $\Theta(\theta)$.  It is easy to see
# 
# $Y_l^m(\theta,\phi)\propto P_l^{|m|}(\cos\theta)e^{im\phi}$.
# 
# These are the spherical harmonics. We will now look at some of these.

# In[10]:


# make two plots of the same spherical harmonic
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm
get_ipython().run_line_magic('matplotlib', 'inline')
from scipy.special import sph_harm
def plot_spherical_harmonic(m,l,theta=np.linspace(0,np.pi,100),phi=np.linspace(0,2*np.pi,100)):
    THETA, PHI = np.meshgrid(theta, phi)
    X = np.sin(THETA) * np.cos(PHI)
    Y = np.sin(THETA) * np.sin(PHI)
    Z = np.cos(THETA)
    # Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
    fcolors = sph_harm(m, l, PHI, THETA).real
    s = sph_harm(m, l, PHI, THETA).real
    s /= s.max()
    fmax, fmin = fcolors.max(), fcolors.min()
    fcolors = (fcolors - fmin)/(fmax - fmin)
    

    # Set the aspect ratio to 1 so our sphere looks spherical
    fig = plt.figure(figsize=(24,12),dpi= 80, facecolor='w', edgecolor='k')
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.plot_surface(X, Y, Z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
    ax.set_axis_off()
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    ax.plot_surface(X*s, Y*s, Z*s,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
    # Turn off the axis planes
    ax.set_axis_off()
    plt.show();
plot_spherical_harmonic(0,1)
plot_spherical_harmonic(0,2);


# ## Solution to the $r$ differential equation

# It is the goal of this section to solve the equation
# \begin{equation}
# \frac{\hbar^2}{R(r)}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right)R(r) + 2mr^2E = E_{\theta,\phi}
# \end{equation}
# 
# The above equation is a form of the Helmholtz equation and has the solution of the spherical Bessel functions.

# In[12]:


import matplotlib.pyplot as plt
from scipy.special import spherical_jn
x = np.arange(0.0, 10.0, 0.01)
fig, ax = plt.subplots()
ax.set_ylim(-0.5, 1.5)
ax.set_title(r'Spherical Bessel functions $j_n$')
for n in np.arange(0, 4):
    ax.plot(x, spherical_jn(n, x), label=rf'$j_{n}$')
plt.legend(loc='best')
plt.show()


# In[16]:


import matplotlib.pyplot as plt
from scipy.special import spherical_yn
x = np.arange(0.0, 10.0, 0.01)
fig, ax = plt.subplots()
ax.set_ylim(-1.5, 0.5)
ax.set_title(r'Spherical Bessel functions $y_n$')
for n in np.arange(0, 4):
    ax.plot(x, spherical_yn(n, x), label=rf'$y_{n}$')
plt.legend(loc='best')
plt.show()

