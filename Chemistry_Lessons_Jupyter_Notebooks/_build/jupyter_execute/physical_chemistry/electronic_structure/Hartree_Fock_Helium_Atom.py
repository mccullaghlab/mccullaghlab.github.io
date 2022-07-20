#!/usr/bin/env python
# coding: utf-8

# # Hartree-Fock

# For the sake of these notes we will consider a two electron Hamiltonian of the form
# 
# $ \hat{H} = -\frac{1}{2}\nabla_1^2 + \frac{Z}{r_1} - \frac{1}{2}\nabla_2^2 + \frac{Z}{r_2} + \frac{1}{r_{12}} $,
# 
# where $r_1$ is the distance of electron $1$ from the nucleus, $Z$ is the nuclear charge and $r_{12}$ is the distance between the two electrons.  The goal, as always, is to solve the Schrodinger equation with this Hamiltonian.  This will yield the energies and wavefunctions of the system. 
# 
# As we have seen, the Schrodinger equation with this Hamiltonian cannot be solved analytically.  Instead, we have previously utilized both perturbation theory and the variational approach to solve for the approximate energies and wavefunctions.  The Hartree-Fock approach is based heavily on the variational approach so will focus on that.
# 
# For the variational approach, we posit a trial wavefunction with variational parameters.  We previously used
# 
# $\psi_t(r_1,r_2) = \frac{\zeta^3}{\pi}e^{-\zeta(r_1+r_2)}$,
# 
# due to this being a product of two Hydrogen 1s like orbitals and $\zeta$ is our variational parameter.  The first step in the variational approach is to solve for the expectation value of the Hamiltonian in this trial wavefunction.  We showed that
# 
# $\langle E \rangle_t = \frac{ \langle \psi_t | \hat{H} | \psi_t \rangle}{\langle \psi_t | \psi_t \rangle} = \zeta^2 - \frac{27}{8}\zeta$
# 
# We then minimize $\langle E \rangle_t$ with respect to the variational parameter
# 
# $\frac{d \langle E \rangle_t}{d\zeta}|_{\zeta_{min}} = 0$
# 
# $\Rightarrow \zeta_{min} = \frac{27}{16}$
# 
# $\Rightarrow E_{min} = -2.8477$ Hartree

# In following this procedure for the variational approach we have completely ignored the spin of the electron.  Recall from postulate 6 of quantum mechanics, an electronic wavefunction must be antisymmetric with respect to exhange of an electron.

# ## Antisymmetric Wavefunctions

# We need to ensure that our wavefunctions our antisymmetric with respect to exchange of electrons.  This will be achieved using slater determinants of spin-orbits.  For the minimal basis helium atom we will have two spin orbits
# 
# $\chi_1(1) = \psi_1(1) \alpha(1)$
# 
# $\chi_2(1) = \psi_1(1) \beta(1)$
# 
# where $\psi_1$ is the spatial function and $\alpha$ and $\beta$ are the spin functions. The $1$ in parantheses, $(1)$, denotes that these are functions of the coordinates of electron $1$.   We construct an antisymmetric wavefunction using Slater determinants:
# 
# $|\Psi(1,2)\rangle = \frac{1}{\sqrt{2}} \begin{vmatrix}\chi_1(1) & \chi_1(2)\\ \chi_2(1) & \chi_2(2)  \end{vmatrix}  = \frac{1}{\sqrt{2}} (\chi_1(1)\chi_2(2) - \chi_2(1) \chi_1(2)) $

# ## Coulomb and Exhange Integrals

# We now want to determine the expectation value of the Hamiltonian in the antisymmetric trial wavefunction.  Let's start by rewriting the Hamiltonian as
# 
# \begin{align}
# \hat{H} &=& -\frac{1}{2}\nabla_1^2 + \frac{Z}{r_1} - \frac{1}{2}\nabla_2^2 + \frac{Z}{r_2} + \frac{1}{r_{12}} \\
#  &=& h(1) + h(2) + \frac{1}{r_{12}},
# \end{align}
# 
# where $h(i) = -\frac{1}{2}\nabla_i^2 + \frac{Z}{r_i}$ is the one electron operator for electron $i$.  Using this form of Hamiltonian, we can write the expectation value
# 
# \begin{align}
# \langle\Psi(1,2)| \hat{H} |\Psi(1,2)\rangle &=& \langle\Psi(1,2)| h(1) + h(2) + \frac{1}{r_{12}} |\Psi(1,2)\rangle\\
# &=& \langle\Psi(1,2)| h(1)|\Psi(1,2)\rangle  +  \langle\Psi(1,2)|h(2)|\Psi(1,2)\rangle +  \langle\Psi(1,2)|\frac{1}{r_{12}} |\Psi(1,2)\rangle
# \end{align}
# 
# We can now evaluate these terms separately.  
# 
# \begin{align}
# \langle\Psi(1,2)| h(1)|\Psi(1,2)\rangle &=& \langle \frac{1}{\sqrt{2}} (\chi_1(1)\chi_2(2) - \chi_2(1) \chi_1(2)) | h(1) | \frac{1}{\sqrt{2}} (\chi_1(1)\chi_2(2) - \chi_2(1) \chi_1(2)) \rangle\\
# &=& \frac{1}{2} \langle(\chi_1(1)\chi_2(2) - \chi_2(1) \chi_1(2)) | h(1) | (\chi_1(1)\chi_2(2) - \chi_2(1) \chi_1(2)) \rangle \\
# &=& \frac{1}{2} \left( \langle\chi_1(1)\chi_2(2) | h(1) | \chi_1(1)\chi_2(2) \rangle + \langle\chi_2(1)\chi_1(2) | h(1) | \chi_2(1)\chi_1(2) \rangle - \langle\chi_1(1)\chi_2(2) | h(1) | \chi_2(1)\chi_1(2) \rangle - \langle\chi_2(1)\chi_1(2) | h(1) | \chi_1(1)\chi_2(2) \rangle \right)\\
# &=& \frac{1}{2}\left(\langle\chi_1(1) | h(1) | \chi_1(1) \rangle + \langle\chi_2(1) | h(1) | \chi_2(1) \rangle\right)
# \end{align}
# 
# where the last equality holds if $\chi_1$ is orthogonal to $\chi_2$.  It can be readily seen that you will get a similar result for the $h(2)$ term.
# 
# The two electron integral is more interesting.  
# 
# \begin{align}
# \langle\Psi(1,2)| \frac{1}{r_{12}}|\Psi(1,2)\rangle &=& \langle \frac{1}{\sqrt{2}} (\chi_1(1)\chi_2(2) - \chi_2(1) \chi_1(2)) | \frac{1}{r_{12}} | \frac{1}{\sqrt{2}} (\chi_1(1)\chi_2(2) - \chi_2(1) \chi_1(2)) \rangle\\
# &=&\frac{1}{2} \left( \langle\chi_1(1)\chi_2(2) | \frac{1}{r_{12}} | \chi_1(1)\chi_2(2) \rangle + \langle\chi_2(1)\chi_1(2) | \frac{1}{r_{12}} | \chi_2(1)\chi_1(2) \rangle - \langle\chi_1(1)\chi_2(2) |\frac{1}{r_{12}} | \chi_2(1)\chi_1(2) \rangle - \langle\chi_2(1)\chi_1(2) | \frac{1}{r_{12}} | \chi_1(1)\chi_2(2) \rangle \right)\\
# &=& \langle\chi_1(1)\chi_2(2) | \frac{1}{r_{12}} | \chi_1(1)\chi_2(2) \rangle - \langle\chi_1(1)\chi_2(2) |\frac{1}{r_{12}} | \chi_2(1)\chi_1(2) \rangle
# \end{align}
# 
# where the final equality holds due to commutative multiplication in this case.  The first term in the last line is called a Coulomb integral and the second term is called an Exchange integral.

# ## Self-consistent Field (SCF)

# 
