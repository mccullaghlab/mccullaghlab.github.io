#!/usr/bin/env python
# coding: utf-8

# In[1]:


# plot of some of the Legendre polynomials
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
from scipy.special import lpmv
x = np.arange(-1,1,0.0001)
plt.figure(figsize=(12,6),dpi= 80, facecolor='w', edgecolor='k')
plt.tick_params(axis='both',labelsize=20)
plt.grid(b=True, which='major', axis='both', color='#808080', linestyle='--')
for l in range(2):
    for m in range(-l,l+1):
        label = "l=" + str(l) + ", m=" + str(m)
        plt.plot(x,lpmv(m,l,x),lw=4,label=label)
plt.title("Associated Lengedre Polynomials",fontsize=16)
plt.legend(fontsize=16);


# In[2]:


from scipy import integrate


# In[3]:


def aleg(x,l,m):
    return lpmv(m,l,x)
def aleg2(x,l,m):
    return lpmv(m,l,x)**2


# In[4]:


integrate.quad(aleg,-1.0,1.0,args=(1,0))[0]


# In[5]:


h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
print(hbar**2*np.pi**2/(2*me*r0**2))
print(hbar**2*np.pi**2/(2*me*r0**2)*6.022e23/1000)


# In[6]:


import numpy as np
h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
n2 = 2
n1 = 1
print(hbar**2/(2*me*r0**2)*6.022e23/1000*(spherical_jn_zero(1,1)**2-spherical_jn_zero(0,1)**2))


# In[8]:


import numpy as np
h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
print(hbar**2/(2*me*r0**2)*6.022e23/1000*(spherical_jn_zero(1,1)**2))


# In[2]:


import numpy as np
h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
n2 = 7.725
n1 = 4.493
print(hbar**2/(2*me*r0**2)*6.022e23/1000*(n2**2-n1**2))


# In[1]:


# make two plots of the same spherical harmonic
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm
from scipy.special import spherical_jn
from scipy.special import lpmv
get_ipython().run_line_magic('matplotlib', 'inline')
from scipy.optimize import root
    

def spherical_jn_zero(l, n, ngrid=100):
    """Returns nth zero of spherical bessel function of order l
    """
    if l > 0:
        # calculate on a sensible grid
        x = np.linspace(l, l + 2*n*(np.pi * (np.log(l)+1)), ngrid)
        y = spherical_jn(l, x)
    
        # Find m good initial guesses from where y switches sign
        diffs = np.sign(y)[1:] - np.sign(y)[:-1]
        ind0s = np.where(diffs)[0][:n]  # first m times sign of y changes
        x0s = x[ind0s]
    
        def fn(x):
            return spherical_jn(l, x)

        return [root(fn, x0).x[0] for x0 in x0s][-1]
    else:
        return n*np.pi
    
def particle_in_sphere_wf(r,theta,phi,n,l,m):
    denom = spherical_jn_zero(l, n)
    return sph_harm(m, l, phi, theta).real*spherical_jn(l, r*denom)


# In[2]:


import numpy as np
h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
for l in range(4):
    for n in range(1,5):
        E = spherical_jn_zero(l, n)
        print(l,n,E**2,hbar**2/(2*me*r0**2)*6.022e23/1000*E**2)


# In[3]:


import numpy as np
h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
for l in range(4):
    n = 1
    E1 = spherical_jn_zero(l, n)
    E2 = spherical_jn_zero(l+1, n)
    print(l,l+1,E2**2-E1**2,hbar**2/(2*me*r0**2)*6.022e23/1000*(E2**2-E1**2))


# In[5]:


import numpy as np
h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
for n in range(1,5):
    l = 0
    E1 = spherical_jn_zero(l, n)
    E2 = spherical_jn_zero(l, n+1)
    print(n,n+1,E2**2-E1**2,hbar**2/(2*me*r0**2)*6.022e23/1000*(E2**2-E1**2))


# In[6]:


import numpy as np
h = 6.626e-34
hbar = h/(2*np.pi)
me = 9.1093e-31 # kg
r0 = 10e-9
for n in range(1,5):
    l = 1
    E1 = spherical_jn_zero(l, n)
    E2 = spherical_jn_zero(l, n+1)
    print(n,n+1,E2**2-E1**2,hbar**2/(2*me*r0**2)*6.022e23/1000*(E2**2-E1**2))


# In[ ]:




