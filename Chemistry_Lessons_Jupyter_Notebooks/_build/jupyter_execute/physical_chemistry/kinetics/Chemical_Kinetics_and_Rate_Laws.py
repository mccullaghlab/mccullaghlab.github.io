#!/usr/bin/env python
# coding: utf-8

# # Chemical Kinetics and Rate Laws

# ## Motivation
# 
# The rate of a chemical reaction is a measure of how fast the reaction occurs.  This rate will not be the same for all reactions and will not stay constant as a reaction progresses.  How fast a reaction happens is an extremely important aspect of a reaction for all sorts of things.  For example, how fast food rots/decays is important in our daily lives.  All food will rot/decay but we often choose to slow this down by putting it in the refrigerator.  

# ## Learning Goals
# 
# After working throught these notes, you should be able to:
# 
# 1. Write out the differential rate law (rate of reaction in terms of rate of change of reactants and products) for an arbitrary reaction
# 2. Differentiate between a *differential rate law* and a *rate law*
# 3. Describe what is meant by overall order of reaction and order of reaction with respect to a specific reactant
# 4. Identify the correct units for the rate constant based on overall rate of reaction
# 5. Use the mehtod of initial rates to estimate the orders of reaction with respect to reactants and the rate constant

# ## Coding Concepts
# 
# The following coding concepts are used in this notebook:
# 
# 1. numpy

# ## The Differential Rate Law

# The rate of a chemical reaction, denoted $v(t)$, is related the rate of dissapearance of the reactants and appearance of the products.  The equation relating all of these things is called the differential rate law.  
# 
# The rate of disappearance and appearance of products in a certain reaction will depend on their stoichiometric coefficients.  In the reaction
# \begin{equation}
# 2NO(g) + O_2(g) \rightarrow 2NO_2(g),
# \end{equation}
# for example, $NO(g)$ will dissapear twice as fast as $O_2(g)$ because two moles of $NO$ are required to react with one mole of $O_2$ to form $2NO_2$.  Thus, when we define a rate of reaction, we must do so in a consistent manner with regards to the stoichiometric coefficients.   
# 
# Consider the general reaction
# \begin{equation}
# aA + bB \rightarrow cC + dD
# \end{equation}
# We define the extent of reaction, $\xi$ using the moles of each reactant and product
# \begin{eqnarray}
# n_A(t) &=& n_A(0) - a\xi(t) \\
# n_B(t) &=& n_B(0) - b\xi(t) \\
# n_C(t) &=& n_C(0) + c\xi(t) \\
# n_D(t) &=& n_D(0) + d\xi(t) 
# \end{eqnarray}
# This is similar to what we did in our discussion of Chemical Equilibria but here we have noted that the number of moles of reactants/products and the extent of reaction are time dependent quantities.
# 
# The rate of disappearance of reactants $A$ and $B$ as well as the rate of appearance of $C$ and $D$ can be related to the single rate of change of the extent of reaction.  To see this, we differentiate the four equations above with respect to time, $t$,:
# \begin{eqnarray}
# \frac{dn_A(t)}{t} &=& - a\frac{d\xi(t)}{dt} \\
# \frac{dn_B(t)}{t} &=& - b\frac{d\xi(t)}{dt} \\
# \frac{dn_C(t)}{t} &=&  c\frac{d\xi(t)}{dt} \\
# \frac{dn_D(t)}{t} &=&  d\frac{d\xi(t)}{dt} 
# \end{eqnarray}
# 
# Each equation can be rearranged to solve for $\frac{d\xi(t)}{dt}$ to yield:
# \begin{eqnarray}
# -\frac{1}{a}\frac{dn_A(t)}{t} &=& \frac{d\xi(t)}{dt} \\
# -\frac{1}{b}\frac{dn_B(t)}{t} &=& \frac{d\xi(t)}{dt} \\
# \frac{1}{c}\frac{dn_C(t)}{t} &=&  \frac{d\xi(t)}{dt} \\
# \frac{1}{d}\frac{dn_D(t)}{t} &=&  \frac{d\xi(t)}{dt} 
# \end{eqnarray}
# The rate of reaction is defined as the rate of change of extent of the reaction yielding
# \begin{equation}
# v(t) = \frac{d\xi(t)}{dt} = -\frac{1}{a}\frac{dn_A(t)}{t} = -\frac{1}{b}\frac{dn_B(t)}{t} = \frac{1}{c}\frac{dn_C(t)}{t} = \frac{1}{d}\frac{dn_D(t)}{t}
# \end{equation}
# 
# We typically measure amount of reactants and/or products in concentration, so we typically write the differential rate law in terms of these by dividing the above equation by volumne, $V$:
# \begin{equation}
# v(t) = \frac{1}{V}\frac{d\xi(t)}{dt} = -\frac{1}{a}\frac{d[A]}{t} = -\frac{1}{b}\frac{d[B]}{t} = \frac{1}{c}\frac{d[C]}{t} = \frac{1}{d}\frac{d[D]}{t}
# \end{equation}
# This equation is what is called the ***differential rate law*** as it relates the rate of reaction to the time derivatives of the concentrations of products and reactants.

# ### Example: Determine the Differential Rate Law 

# Write out the differential rate law for the following balanced chemical reaction
# 
# \begin{equation}
# 2NO(g) + 2H_2(g) \rightarrow N_2(g) + 2H_2O(g)
# \end{equation}

# The differential rate law equates the rate of reaction to the rate of dissappearance of *all reactants* and the rate of appearance of *all products*.  For this reaction it looks like:
# \begin{equation}
# v(t) = -\frac{1}{2}\frac{d[NO]}{dt} = -\frac{1}{2}\frac{d[H_2]}{dt} = \frac{d[N_2]}{dt} = \frac{1}{2}\frac{d[H_2O]}{dt}
# \end{equation}

# ## The Rate Law

# The rate law for a chemical reaction is an equation that relates the rate of a reaction to the concentrations of the reactants.  For example, our reaction
# \begin{equation}
# 2NO(g) + O_2(g) \rightarrow 2NO_2(g) 
# \end{equation}
# Has the ***differential rate law***
# \begin{equation}
# v(t) = -\frac{1}{2}\frac{d[NO]}{dt} = -\frac{d[O_2]}{dt} = \frac{1}{2}\frac{d[NO_2]}{dt}
# \end{equation}
# but the ***rate law*** for this reaction, one that relates the rate, $v(t)$, to the concentrations of reactants cannot be inferred from the balanced chemical reaction.  The rate law for this reaction has been experimentally determined to be
# \begin{equation}
# v(t) = k[NO]^2[O_2]
# \end{equation}
# 
# The rate law introduces a few new quantities and concepts to our discussion of rate of reaction.  One of these, $k$, is called the rate constant.  This important quantity determines the proportionality or scaling between the concentrations of reactants to their appropriate power and the rate of a reaction.  
# 
# The rate law also introduces the concept of *order of reaction*.  This reaction is said to be *second order with resepect to $NO_2$*, *first order with respect to $O_2$*, and *third order overall*.  Second order with resepect to $NO_2$ can be inferred from the rate law because $[NO]$ is raised to the second power.  First order with respect to $O_2$ can be inferred from the rate law because $[O_2]$ is raised to the first power.  Third order overall simply stems from adding up all powers of the reactants in the rate law to get a value of three.  
# 
# The order of a reaction for any reactant ***cannot*** be inferred from the overall balanced chemical reaction. Notice, however, that, in this case, the order of reaction follows the stoichiometry of the balanced chemical reaction.  
# 
# The overall order of the reaction dictate the units of the force constant.  The rate of a reaction is always in units of amount per time (typically concentration per time).  The rate law equates concentration per time to some overall power of concentration times $k$ and thus the units of $k$ must allow the two sides to be equal.  Here is a table to help with this concept
# 
# | Rate Law    | Order | Units of $k$    |
# |:------------|:------| :------------------- |
# |$v = k$      | 0     |  M$\cdot$s$^{-1}$      |
# | $v =k[A]$   | 1     |  s$^{-1}$              |
# | $v =k[A]^2$ | 2     |  M$^{-1}\cdot$s$^{-1}$ |
# | $v =k[A]^3$ | 3     |  M$^{-2}\cdot$s$^{-1}$ |

# ## Rate Laws Must be Determined Experimentally

# Rate laws must be determined experimentally.  The order of a reaction cannot be inferred by the stoichiometry of the balanced chemical reaction or any other *a priori* information about the reaction.  
# 
# To see how we might go about determining an order of reaction consider the generic reaction
# \begin{equation}
# aA + bB \rightarrow cC + dD
# \end{equation}
# The generic rate law can be written as
# \begin{equation}
# v(t) = k[A]^{m_A}[B]^{m_B}
# \end{equation}
# It is our goal to determine the powers $m_A$ and $m_B$.  
# 
# Here we consider the method of measuring initial rate (with some experimental technique) and varying the initial concentrations of $A$ and $B$ in different trials.  It is easiest to see how this method works in an example given below.

# ### Example: Rate Law from Initial Rates

# Consider the following initial rate data for the reaction
# \begin{equation}
# 2NO_2(g) + F_2(g) \rightarrow 2NO_2F(g)
# \end{equation}
# 
# | Run/Trial   | $[NO_2]_0$/mol$\cdot$dm$^3$ | $[F_2]_0$/mol$\cdot$dm$^{-3}$ | $v_0$/mol$\cdot$dm$^{-3}\cdot$s$^{-1}$ |
# | :-----------|:-------------------------- | :---------------------------- | :------------------------------- |
# | 1   | 1.15 | 1.15 | $6.12\times10^{-4}$ |
# | 2   | 1.72 | 1.15 | $1.36\times10^{-3}$ |
# | 3   | 1.15 | 2.30 | $1.22\times10^{-3}$ |
# 
# Determine the reaction rate law and value of the rate constant.

# We start by writing out the generic form of the rate law
# \begin{equation}
# v = k[NO_2]^{m_{NO_2}}[F_2]^{m_{F_2}}
# \end{equation}
# Here, we are using initial rates which can be written as
# \begin{equation}
# v_0 = k[NO_2]_0^{m_{NO_2}}[F_2]_0^{m_{F_2}}
# \end{equation}
# 
# We now look at the trials/runs provided and recognize that the concentration of $NO_2$ is held constant for trials 1 and 3.  Thus, the change in initial rate between these two runs must be due to the change in initial concentration of $F_2$.  We start by writing out the initial rate law for each of these trials
# \begin{eqnarray}
# \text{Run 1: } 6.12\times10^{-4} &=& k (1.15)^{m_{NO_2}}(1.15)^{m_{F_2}} \\
# \text{Run 3: } 1.22\times10^{-3} &=& k (1.15)^{m_{NO_2}}(2.30)^{m_{F_2}}
# \end{eqnarray}
# 
# We now divide the initial rate law from Run 3 by that of Run 1 to get:
# \begin{eqnarray}
# \frac{1.22\times10^{-3}}{6.12\times10^{-4}} &=& \frac{k (1.15)^{m_{NO_2}}(2.30)^{m_{F_2}}}{k (1.15)^{m_{NO_2}}(1.15)^{m_{F_2}}} \\
# &=& \left( \frac{2.30}{1.15}\right)^{m_{F_2}} \\
# \Rightarrow \ln \frac{1.22\times10^{-3}}{6.12\times10^{-4}} &=& m_{F_2}\ln\frac{2.30}{1.15} \\
# \Rightarrow m_{F_2} &=& \frac{\ln \frac{1.22\times10^{-3}}{6.12\times10^{-4}}}{\ln\frac{2.30}{1.15}} \\
# &=& 0.9952776 \approx 1
# \end{eqnarray}

# In[1]:


import numpy as np
np.log(1.22e-3/6.12e-4)/np.log(2.30/1.15)


# To now determine $m_{NO_2}$ we use runs 1 and 2
# \begin{eqnarray}
# \text{Run 1: } 6.12\times10^{-4} &=& k (1.15)^{m_{NO_2}}(1.15)^{m_{F_2}} \\
# \text{Run 2: } 1.36\times10^{-3} &=& k (1.72)^{m_{NO_2}}(1.15)^{m_{F_2}}
# \end{eqnarray}
# 
# We now divide the initial rate law from Run 2 by that of Run 1 to get:
# \begin{eqnarray}
# \frac{1.36\times10^{-3}}{6.12\times10^{-4}} &=& \frac{k (1.72)^{m_{NO_2}}(1.15)^{m_{F_2}}}{k (1.15)^{m_{NO_2}}(1.15)^{m_{F_2}}} \\
# &=& \left( \frac{1.72}{1.15}\right)^{m_{NO_2}} \\
# \Rightarrow \ln \frac{1.36\times10^{-3}}{6.12\times10^{-4}} &=& m_{NO_2}\ln\frac{1.72}{1.15} \\
# \Rightarrow m_{NO_2} &=& \frac{\ln \frac{1.46\times10^{-3}}{6.12\times10^{-4}}}{\ln\frac{1.72}{1.15}} \\
# &=& 1.98 \approx 2
# \end{eqnarray}

# In[2]:


import numpy as np
np.log(1.36e-3/6.12e-4)/np.log(1.72/1.15)


# With the orders of reaction determined with respect to all reactants we can now determine the value of the rate constant.  This can be done for with any run.  Here I choose to use Run 1:
# \begin{eqnarray}
# \text{Run 1: } 6.12\times10^{-4} &=& k (1.15)^{2}(1.15)^{1} \\
# \Rightarrow k &=& \frac{6.12\times10^{-4}}{1.15^3} \\
# &=& 4.02\times10^{-4} \text{ dm}^6\cdot\text{mol}^{-2}\text{s}^{-1}
# \end{eqnarray}

# In[3]:


6.12e-4/1.15**3

