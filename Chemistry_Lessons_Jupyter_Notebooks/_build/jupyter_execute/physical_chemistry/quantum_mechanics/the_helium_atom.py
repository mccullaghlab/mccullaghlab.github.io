#!/usr/bin/env python
# coding: utf-8

# # The Helium Atom

# The helium atom is composed of a two electrons, two protons and two neutrons.  If these particles were classical, the interaction energy would be the Coulombic attraction between electron and the nuclei as well as the Coulombic repulsion between the two electrons (assuming the nuclei are stationary):
# 
# $V(\vec{r_1},\vec{r_2}) = -\frac{2e^2}{4\pi\epsilon_0|\vec{r_1}|}-\frac{2e^2}{4\pi\epsilon_0|\vec{r_2}|} + \frac{e^2}{4\pi\epsilon_0|\vec{r_1}-\vec{r_2}|}$,
# 
# where $\vec{r_1}$ is the separation vector between the electron 1 and the nucleus, $r_2$ is the separation vector between electron 2 and the nucleus, $e$ is the charge of an electron and $\epsilon_0$ is the permittivity of free space. In atomic units the potential simply becomes
# 
# $V(\vec{r_1},\vec{r_2}) = -\frac{2}{|\vec{r_1}|}-\frac{2}{|\vec{r_2}|} + \frac{1}{|\vec{r_1}-\vec{r_2}|}$.
# 
# In the quantum mechanical picture, we will treat this same function as the potential energy operator.  The Hamiltonian in atomic units for a helium atom is then, including the kinetic energy,
# 
# $\hat{H}(\vec{r_1},\vec{r_2}) = -\frac{1}{2}\nabla_1^2-\frac{1}{2}\nabla_2^2-\frac{2}{|\vec{r_1}|}-\frac{2}{|\vec{r_2}|} + \frac{1}{|\vec{r_1}-\vec{r_2}|}$
# 
# where $\nabla_i^2$ is the Laplacian for particle $i$ in three dimensions.  We can rearrange the right hand side of the above equation to yield
# 
# $\hat{H}(\vec{r_1},\vec{r_2}) = -\frac{1}{2}\nabla_1^2-\frac{2}{|\vec{r_1}|}-\frac{1}{2}\nabla_2^2-\frac{2}{|\vec{r_2}|} + \frac{1}{|\vec{r_1}-\vec{r_2}|}$
# 
# $ = \hat{H}^{Z=2}_H(\vec{r_1}) + \hat{H}^{Z=2}_H(\vec{r_2}) + \frac{1}{|\vec{r_1}-\vec{r_2}|}$,
# 
# where $\hat{H}^{Z=2}_H(\vec{r_i})$ is the Hydrogen-like Hamiltonian for particle $i$ with a nuclear charge of 2$e$.

# The Schrodinger equation for the Helium atom cannot be solved analytically.  We will show how to approximate the ground-state energy using first-order perturbation theory and the variational approach.

# ## First-order Perturbation Theory

# To apply first-order perturbation theory we must write the Hamiltonian asa zeroth-order Hamiltonian plus a perturbation.  Here we treat the electron-electron repulsion as the perturbation (not a great approximation as we will see)
# 
# $\hat{H}(\vec{r_1},\vec{r_2}) = \hat{H}^0 + \hat{H}^1$,
# 
# where
# 
# $\hat{H}^0 = \hat{H}^{Z=2}_H(\vec{r_1}) + \hat{H}^{Z=2}_H(\vec{r_2})$
# 
# and 
# 
# $\hat{H}^1 = \frac{1}{|\vec{r_1}-\vec{r_2}|}$.
# 
# In first-order perturbation theory, the total energy of the system is approximated as the eigenvalue of the zeroth-order Schrodinger equation plus the expectation value of the perturbation in the zeroth-order wavefunctions.  The zeroth-order wavefunction is the solution to the Schrodingder equation
# 
# $H^0(\vec{r_1},\vec{r_2}) \Psi^0(\vec{r_1},\vec{r_2})= E^0\Psi^0(\vec{r_1},\vec{r_2})$.
# 
# 
# Since we can write the Hamiltonian here as a sum of independent hydrogen-like Hamiltonians, the wavefunctions (eigenfunctions of this Hamiltonian) are products of hydrogen-like wavefunctions
# 
# $\Psi^0(\vec{r_1},\vec{r_2}) = \phi_{1s}^{Z=2}(\vec{r_1})\phi_{1s}^{Z=2}(\vec{r_2})$,
# 
# where 
# 
# $\phi_{1s}^{Z}(\vec{r_1}) = \frac{Z^{3/2}}{\sqrt{\pi}}\exp(-Z|\vec{r_1}|)$
# 
# can be found in Table 7.5 (page 340) of McQuarrie Quantum Chemistry.  The zeroth order energy is then a sum of hydrogen-like energies
# 
# $E^0 = -\frac{1}{2} - \frac{1}{2} = -1$ Hartree.
# 
# First order perturbation theory dictates that 
# 
# $E^1 = \langle \Psi^0 (\vec{r_1},\vec{r_2})| \frac{1}{|\vec{r_1}-\vec{r_2}|} | \Psi^0 (\vec{r_1},\vec{r_2})\rangle $.
# 
# 

# ## Variational Approach

# In[ ]:




