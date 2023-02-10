#!/usr/bin/env python
# coding: utf-8

# # Operator Algebra

# ## Motivation
# 
# Operators are an important component of quantum mechanics.  Postulate 2 states that all physical observables have an associated (linear and Hermitian) operator.  Each of these operators can be built from the basic position and momentum operators.  To do so, we need to know the rules of operator algebra.

# ## Learning Goals
# 
# After working through these notes, you will be able to:
# 
# 1. Determine if two operators are equal
# 2. Determine the single operator equivalent to the addition of other operators
# 3. Determine the single operator equivalent to the multiplication of other operators
# 4. Determine the commutator of two operators.

# ## Coding Concepts
# 
# The following coding concepts are used in this notebook:
# 
# 1. None. Yet.

# ## Operator Equality

# Two operators $\hat{A}$ and $\hat{B}$ are said to be equivalent if they achieve the same result when applied to ***any*** function.  Assuming that $\hat{A}$ and $\hat{B}$ depend only on a single variable $x$ then if
# \begin{align}
# \hat{A}f(x) &= \hat{B}f(x)
# \end{align}
# for all functions $f(x)$ then $\hat{A} = \hat{B}$.
# 
# It is important to note that it is not assumed that $f(x)$ is an eigenfunction of $\hat{A}$.  It could be that applying $\hat{A}$ to $f(x)$ creates a new function $g(x)$, $\hat{A}f(x) = g(x)$.  In order for $\hat{A}$ and $\hat{B}$ to be equal $\hat{B}f(x)=g(x)$ as well.
# 
# It is also important to note that this behavior must be consistent for any function $f(x)$.  It is possible to find a particular $f(x)$ for which the result of two different operators is the same.  For example
# \begin{align}
# f(x) &= e^{-x^2} \\
# \hat{A} &= \frac{d}{dx} \\
# \hat{B} &= -2x
# \end{align}
# 
# If we now apply to the two operators to $f(x)$ we get
# \begin{align}
# \hat{A}f(x) &= \frac{d}{dx} e^{-x^2} = -2xe^{-x^2} \\
# \hat{B}f(x) &= -2x e^{-x^2} 
# \end{align}
# We see that the result is the same for this particular function.  We cannot say that $\hat{A}$ and $\hat{B}$ are equivalent, however, because this behavior does not hold true for other functions (e.g. $f(x) = x^2$).

# ## Operator Addition

# Operator addition (and substraction) follows what we might consider to be *standard* rules.  That is if
# \begin{align}
# \hat{A}f &= g \\
# \hat{B}f &= h
# \end{align}
# then
# \begin{align}
# \left(\hat{A} + \hat{B}\right)f = \hat{A}f + \hat{B}f = g + h
# \end{align}
# 
# This can also be thought of as the distributive law for functions over operators.  

# ### Example: Operator Addition
# 
# Given $\hat{A} = \frac{d}{dx}$ and $\hat{B} = 3x^2 -x$, determine the result of $(\hat{A} + \hat{B})f(x)$ for $f(x) = 3e^{x^2}.
# 
# \begin{align}
# (\hat{A} + \hat{B})f(x) &= \left( \frac{d}{dx}+  3x^2 -x \right) 3e^{x^2}  \\
# &= \frac{d}{dx}\left(3e^{x^2}\right) +  3x^2\left(3e^{x^2}\right) -x\left(3e^{x^2}\right) \\
# &= 6xe^{x^2} +  9x^2e^{x^2} -3xe^{x^2} \\
# &= 3xe^{x^2} +  9x^2e^{x^2} \\
# &= \left(3x +  9x^2\right)e^{x^2}
# \end{align}

# ## Operator Multiplication

# Operator multiplication does not follow *standard* rules.  Specifically, some operators do not commute with others.  Operator multiplication, however, is extremely important for building up the set of operators pertaining to observables.  Thus, it is important to define what we mean by
# \begin{equation}
# \hat{A}\hat{B}
# \end{equation}
# for generic operators $\hat{A}$ and $\hat{B}$.  
# 
# The operation $\hat{A}\hat{B}$ means to first apply operator $\hat{B}$ to a function and then apply operator $\hat{A}$.  This is perhaps more concisely written as
# \begin{equation}
# \hat{A}\hat{B}f(x) = \hat{A}\left[\hat{B}f(x)\right]
# \end{equation}
# 
# Note that since operators do not commute we should typically assume that
# \begin{equation}
# \hat{A}\hat{B} \neq \hat{B}\hat{A}
# \end{equation}
# 
# This implies that the order in which you apply the operator is hugely important.  

# ### Example: Operator Multiplication
# 
# Given $\hat{A} = \frac{d}{dx}$ and $\hat{B} = 3x^2 -x$, determine the result of $\hat{A}\hat{B}f(x)$ and $\hat{B}\hat{A}f(x)$ for $f(x) = 3e^{x^2}$.
# 
# \begin{align}
# \hat{A}\hat{B}f(x) &= \frac{d}{dx}\left[\left( 3x^2 -x \right) 3e^{x^2}\right]  \\
# &= \frac{d}{dx}\left[9x^2e^{x^2} -3xe^{x^2}\right] \\
# &= \frac{d}{dx}\left(9x^2e^{x^2}\right) -\frac{d}{dx}\left(3xe^{x^2}\right) \\
# &= 18xe^{x^2} + 18x^3e^{x^2}  -3e^{x^2} - 6x^2e^{x^2} \\
# &= \left( 6x^3 - 2x^2 + 6x -1 \right)3e^{x^2}
# \end{align}
# 
# \begin{align}
# \hat{B}\hat{A}f(x) &= \left( 3x^2 -x \right) \left[\frac{d}{dx}3e^{x^2}\right]  \\
# &= \left( 3x^2 -x \right)\left( 6xe^{x^2}\right) \\
# &= \left( 3x^3 -2x^2 \right) 3e^{x^2}
# \end{align}

# ## Operator Commutation

# Whether operators commute or not turns out to be quite important in quantum mechanics.  We define the ***commutator*** of two operators as
# \begin{equation}
# \left[\hat{A},\hat{B}\right] = \hat{A}\hat{B} - \hat{B}\hat{A}
# \end{equation}
# 
# If this value is zero then the operators commute.  If this value is non-zero, the operators do not commute.  We determine a commutator for arbitrary function.

# ### Example: Commutators
# 
# Estimate the commutator of $\hat{A} = \frac{d}{dx}$ and $\hat{B} = 3x^2 -x$.
# 
# We know that these two operators don't commute. Let's compute the commutator, however.
# 
# \begin{align}
# [\hat{A},\hat{B}]f(x) &= [\hat{A}\hat{B} - \hat{B}\hat{A}]f(x) \\
# &= \hat{A}\hat{B}f(x) - \hat{B}\hat{A}]f(x) \\
# &= \frac{d}{dx}\left[\left( 3x^2 -x \right) f(x)\right] - \left( 3x^2 -x \right) \left[\frac{d}{dx}f(x)\right] \\
# &= \frac{d}{dx}\left( 3x^2f -xf \right)- \left( 3x^2 -x \right)\frac{df}{dx} \\
# &=  6xf + 3x^2\frac{df}{dx} - f - x\frac{df}{dx} - 3x^2\frac{df}{dx} +x\frac{df}{dx} \\
# &=  (6x-1)f 
# \end{align}
# Thus,
# \begin{align}
# [\hat{A},\hat{B}] = 6x-1
# \end{align}
# 
# We see that if we plug in that $f(x) = 3e^{x^2}$ and substract the commutator will be the same as the difference between the two answers in the multiplication example.
