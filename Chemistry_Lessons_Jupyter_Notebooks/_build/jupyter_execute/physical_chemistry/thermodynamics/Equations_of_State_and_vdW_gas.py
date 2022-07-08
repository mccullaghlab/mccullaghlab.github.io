#!/usr/bin/env python
# coding: utf-8

# ## Equations of State

# An equation of state is a function that relates all intensive state functions of a system.

# ### Ideal Gas as an example

# The ideal gas law is the most common example of an equation of state.  The equation is typically written as
# 
# $PV = nRT$,
# 
# where $P$ is the pressure, $V$ is the volume, $n$ is the number of moles, $R$ is the universal gas constant, and $T$ is the temperature of the system.  Note that $P$ and $T$ are intensive variables while $n$ and $V$ are extensive variables.

# ### An aside: extensive and intensive variables

# An ***extensive*** variable is one that is proportional to the system size (e.g. volume, number of moles, mass)
# 
# An ***intensive*** variable is one that is *not* proportional so the system size (e.g. pressure, temperature, density)
# 
# An extensive variable can be made into an intensive variable by dividing by number of particles/moles.

# An equation of state relates all of the ***intensive*** properties of a system and yet the ideal gas law contains bouth $V$ and $n$ which are ***extensive***.  We can divide both sides of the ideal gas law by $V$ and get
# 
# $ P = \frac{nRT}{V}$
# 
# $\Rightarrow P = \rho RT$
# 
# where $\rho=\frac{n}{V}$ is the volume density (e.g. moles per liter) of particles and is an ***intensive*** variable.  This form of the ideal gas law relates only intensive variables of the system.

# ### State functions are only known for a select few systems

# Knowing these functions is necessary to compute the thermodynamic energy functions (for example what we did for the ideal gas and $\Delta U$).  These are only known for select systems.  Here are some examples:
# 
# $P(V-nb) = nRT$  - hard spheres
# 
# $(P+\frac{n^2a}{V^2})(V-nb) = nRT$  - Van der Waals gas
# 
# $(P + \frac{n^2a}{TV^2})(V-nb) = nRT$ - Berthdot gas
# 
# where $a$ and $b$ are constants for specific gasses. $a$ represents the attraction between particles and $b$ represents the finite volume of particles.

# ## VdW Gas

# How does changing the equation of state affect the work done during a process?  For this, we will investigate a VdW gas in expansion/contraction processes and compare to our previous results for an ideal gas.

# ### The problem

# Compute the work required to compress $n$ moles of a VdW gas from $V_1$ to $V_2$ in the following processes:
# 
# (i) Reversibly and isothermally
# 
# (ii) Against a constant external pressure, $P_{ext}$.

# ### (i) Reversibly and isothermally

# Compute the work to compress $n$ moles of a VdW gas from $V_1$ to $V_2$ during a reversible isothermal process.  Start with 
# 
# $w = -\int_{V_1}^{V_2} P_{ext}dV$.
# 
# We now utilize that the process is done reversibly and thus $P_{ext}=P_{sys}$ every step along the path.  Thus
# 
# $w = -\int_{V_1}^{V_2} P_{sys}dV$.

# We can now use the equation of state of the system and plug in for $P_{sys}$.  First we will need to rearrange the VdW equation of state to have pressure on one side:
# 
# $(P_{sys}+\frac{n^2a}{V^2})(V-nb) = nRT$ 
# 
# $\Rightarrow P_{sys}+\frac{n^2a}{V^2} = \frac{nRT}{V-nb}$
# 
# $\Rightarrow P_{sys} = \frac{nRT}{V-nb} - \frac{n^2a}{V^2}$

# Now we will plug in $P$ into the equation for work:
# 
# $w = -\int_{V_1}^{V_2} P_{sys}dV$
# 
# $ = -\int_{V_1}^{V_2} \left[ \frac{nRT}{V-nb} - \frac{n^2a}{V^2} \right]dV$
# 
# $ = -\int_{V_1}^{V_2} \frac{nRT}{V-nb}dV + \int_{V_1}^{V_2} \frac{n^2a}{V^2} dV$
# 
# Since the process is isothermal, $T$ is constant.  $n$, $a$ and $b$ are also constant.  Thus, 
# 
# $w = -nRT\int_{V_1}^{V_2} \frac{1}{V-nb}dV + n^2a\int_{V_1}^{V_2} \frac{1}{V^2} dV$

# We now solve the integrals from $V_1$ to $V_2$:
# 
# $w = -nRT\ln\left( \frac{V_2-nb}{V_1-nb}\right) - n^2a\left[\frac{1}{V_2}-\frac{1}{V_2}\right]$
# 
# where I used that 
# 
# $\int\frac{1}{V-nb}dV = \ln(V-nb) + C$
# 
# and
# 
# $\int\frac{1}{V^2} dV = -\frac{1}{V} + C$

# ### (ii) Against a constant external pressure

# Compute the work to compress $n$ moles of a VdW gas from $V_1$ to $V_2$ by a constant external pressure, $P_{ext}$.  Start with 
# 
# $w = -\int_{V_1}^{V_2} P_{ext}dV$.
# 
# Since $P_{ext}=constant$, we can pull it out of the integral and get
# 
# $w = -P_{ext}\int_{V_1}^{V_2} dV = -P_{ext}\Delta V$.

# Since $P_{ext}$ will have to be equal to $P_{sys}$ at the end of the process, we can write:
# 
# $w  = -P_{2}\Delta V$
# 
# $ = -\left(\frac{nRT_2}{V_2-nb} - \frac{n^2a}{V_2^2}\right)\Delta V$
