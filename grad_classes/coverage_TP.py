#!/usr/bin/env python
"""
Code for Kreuzer's surface science class to plot the coverage
for different interaction strengths
"""
import Numeric, pylab
e = Numeric.e

P_list = [0.01,0.1,1,10,100,1000,10000,100000]
colors = ['k','b','r','g','m','c','y','k']
v = Numeric.arange(-5.99,12,0.1,typecode=Numeric.Float)
theta = []
B = 1
A = B*e**-v

for P in P_list:
    theta.append( (A*P + ((A*P)**2 - A*P + 2*B*P)*((A*P-1)**2+4*B*P)**-0.5) / \
                  (1 + A*P + ((A*P-1)**2+4*B*P)**0.5) )

last = []
pylab.plot(v,[0.5 for i in v],'r--',linewidth=2,label=r'$\tt{\theta(T,P)\ =\ 0.5}$')
for i in range(len(P_list)):
    pylab.plot(v,theta[i],colors[i]+'-',linewidth=2,label='P = '+str(P_list[i]))
    last.append(theta[i][-1])

pylab.title('Coverage vs. Interaction strength for various pressures')
pylab.xlabel(r'$\tt{V_{nn}}$')
pylab.ylabel(r'$\tt{\theta(T,P)}$')
pylab.legend(loc=0)
pylab.savefig('coverage_TP.png')

pylab.clf()
"""
V_list = [-5,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]
V_list = Numeric.arange(-5.99,10,1,typecode=Numeric.Float)*1.0
V_list = v
T = []
print v
print V_list
for V in V_list:
    T.append([])
    print V
    print theta[list(V_list).index(V)][list(v).index(V)]
    T[list(V_list).index(V)].append(theta[list(V_list).index(V)][list(v).index(V)])
for t in T:
    pylab.plot(P_list,t,'bo-',linewidth=2)

pylab.title('P vs. Theta for various interaction strengths')
pylab.semilogx()
pylab.xlabel('log(P)')
pylab.ylabel(r'$\tt{\theta(T,P)}$')
pylab.legend(loc=0)
pylab.savefig('PvsTheta.png')
"""
