#!/usr/bin/env python
"""
Code for Kreuzer's surface science class to plot the coverage
for different interaction strengths
"""
import Numeric, pylab
e = Numeric.e

v_list = [0,-3,3,10,20,50,100]
colors = ['k','b','r','g','m','c','y']
x = Numeric.arange(0.01,50,0.1,typecode=Numeric.Float)
theta = []

for v in v_list:
    theta.append( ( e**(x-v) + (e**(2*(x-v))-e**(x-v)+2*e**(x))*((e**(x-v)-1)**2+4*e**(x))**-0.5 ) / \
                  ( 1 + e**(x-v) + ( (e**(x-v)-1)**2 + 4*e**(x))**0.5 ) )

for i in range(len(v_list)):
    pylab.plot(x,theta[i],colors[i]+'-',linewidth=2,label='v = '+str(v_list[i]))

pylab.title('Coverage for various interaction strengths')
pylab.xlabel(r'$\tt{(\mu - E_{s}) \ / \ k_{B}T}$')
pylab.ylabel(r'$\tt{\theta(T,\mu)}$')
pylab.legend(loc=0)
pylab.savefig('coverage.png')
