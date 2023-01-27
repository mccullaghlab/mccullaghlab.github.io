#!/usr/bin/env python
# coding: utf-8

# # Schrodinger Equation

# ## Motivation
# 
# The Schrodinger equation is the fundamental equation of quantum mechanics.  In these notes, we introduce the Schrodinger equation and its relationship to the classical wave equation.  Motivated by the Hamiltonian operator in the Schrodinger equation, we then introduce the concept of operators in quantum mechanics.  

# ## Learning Goals
# 
# After working through these notes, you will be able to:
# 
# 1. Write out and identify the Schrodinger equation
# 2. Demonstrate the equivalence of the Schrodinger equation and the classical wave equation
# 3. Define an opertor
# 4. Apply operators to functions 
# 5. Define a linear operator
# 6. Identify linear operators

# ## Coding Concepts
# 
# No coding concepts are used.  Yet.

# ## Schrodinger Equation from the Classical Wave Equation

# The Schrodinger equation is the basis for all of the quantum mechanics we will discuss in this class.  Specifically, the time independent Schrodinger equation is all we will discuss in this class.  This equation has the well known form
# \begin{equation}
# \hat{H}\Psi = E\Psi
# \end{equation}
# where $\bar{H}$ is the Hamiltonian operator, $\Psi$ is the wave function, and $E$ is the energy.  This equation cannot be derived from any more basic equations, rather we have to accept the Schrodinger equations as one of the postulates of quantum mechanics.  
# 
# The classical wave equation can be manipulated into a Schrodinger equation.  To see this, recall that the classical wave equation can be written as 
# \begin{equation}
# \frac{d^2\psi(x)}{dx^2} + \frac{\omega^2}{v^2}\psi(x) = 0
# \end{equation}
# for solution
# \begin{equation}
# u(x,t) = \psi(x)\cos\omega t
# \end{equation}
# 
# Now recall that $\omega = 2\pi\nu$ and $v = \nu\lambda$ meaning we can write the differential equation as
# \begin{align}
# &\frac{d^2\psi(x)}{dx^2} + \frac{4\pi^2\nu^2}{\nu^2\lambda^2}\psi(x) = 0 \\
# &= \frac{d^2\psi(x)}{dx^2} + \frac{4\pi^2}{\lambda^2}\psi(x) = 0
# \end{align}
# 
# Now recall from the de Broglie wave equation that we can write
# \begin{align}
# \lambda &= \frac{h}{mv} \\
# &= \frac{h}{\sqrt{2m[E-V(x)]}}
# \end{align}
# where the last equality holds because $E = K + V = \frac{1}{2}mv^2 + V$.  Plugging this relationship in for $\lambda$ in the classical wave equation we get
# \begin{align}
# &\frac{d^2\psi(x)}{dx^2} + \frac{4\pi^2\cdot2m[E-V(x)]}{2mh^2}\psi(x) = 0 \\
# =& \frac{d^2\psi(x)}{dx^2} + \frac{2m[E-V(x)]}{\hbar^2}\psi(x) = 0 \\
# \Rightarrow &  \frac{-\hbar^2}{2m}\frac{d^2\psi(x)}{dx^2}  + V(x)\psi(x) =  E\psi(x) \\
# \end{align}
# 
# We see that this looks like the Schrodinger equation with 
# \begin{equation}
# \hat{H} =  \frac{-\hbar^2}{2m}\frac{d^2}{dx^2}  + V(x)
# \end{equation}

# ## Operators

# We see that the Hamiltonian in the Schrodinger equation is an *operator*.  Operators are symbols that describe performing a mathematical operation.  Examples of operators include derivative with respect to $x$, $\frac{d}{dx}$, multiply by $x$, simply denoted $x$, and integrate, denoted $\int$.  
# 
# We denote generic operators by characters (often capital letters) with hats over them.  For example below we define the operator $\hat{A}$ as the first derivative with respect to $x$:
# \begin{equation}
# \hat{A} = \frac{d}{dx}
# \end{equation}
# Operators by themselves are pretty useless.  The must be applied to a function. We might apply the operator $\hat{A}$, above, a generic function of one dimension:
# \begin{align}
# \hat{A}f(x) = \frac{d}{dx}f(x) = \frac{df}{dx}
# \end{align}
# 
# If we know the specific function, $f(x) = 2x$ for example, we can determine the result of operating $\hat{A}$ on $f(x)$:
# \begin{align}
# \hat{A}f(x) = \frac{d2x}{dx} = 2
# \end{align}
# 
# Operators are very important in quantum mechanics.  Specifically, applying operators to the wavefunction for a system can be used to determine physical obervables for that system.  Below is a table of classical observables and the corresponding quantum mechanical operator. 

# Classical and quantum correspondence. 
# ![operators](img/operators.png)

# Operators can be divided into linear operators and non-linear operators.  We will focus on linear operators.

# ### Example: Operators
# 
# Perform the following operations
# 
# 1. $\hat{A}(2x)$, $\hat{A} = \frac{d^2}{dx^2}$
# 2. $\hat{A}(x^2)$, $\hat{A} = \frac{d^2}{dx^2} + 2\frac{d}{dx} + 3$
# 3. $\hat{A}(xy^3)$, $\hat{A} = \frac{\partial}{\partial y}$

# ### Linear Operators

# A linear operator, $\hat{A}$, has the following property:
# 
# $\hat{A} \left[ c_1f_1(x) + c_2f_2(x)\right] = c_1\hat{A}f_1(x) + c_2\hat{A}f_2(x)$
# 
# There are plenty of examples of linear operators including differentiation and integration.  Taking the square root, however, is an example of an operator that is *not* linear.  
# 
# In order to determine if an operator is linear, we must apply it to a sum of arbitrary functions and demonstrate that the above equality holds.  

# ### Example: Linear Operators
# 
# Which of the following are linear operators?
# 
# 1. $\hat{A} = \frac{d}{dx}$
# 2. $\hat{A} = \sqrt{}$
# 3. $\hat{A} = \frac{\partial}{\partial y} + \frac{\partial}{\partial x}$
# 
