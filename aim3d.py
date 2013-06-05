#!/usr/bin/env python
"""
Make contour plots of aim.surf files
"""
import os,sys,commands,glob
import pylab,Numeric

def main():

    try:
        f = open(sys.argv[1],'r')
        header1 = f.readline()
        header2 = f.readline()
        header3 = f.readline()
        lines = f.readlines()
        tail = lines.pop(-1)
        f.close()
    except IndexError:
        print '\nusage: '+sys.argv[0]+'  aim.surf\n'
        sys.exit(0)

    x,y,z,c = [],[],[],[]
    for line in lines:
        row = line.split()
        x.append(float(row[0]))
        y.append(float(row[1]))
        z.append(float(row[2]))
        c.append(float(row[3]))

    q = int(header3.split()[0])
    x = x[::q]
    y = y[0:q]
    v,u = [],[]
    for i in range(len(c)/q):
        v.append([])
        u.append([])
        for j in range(q):
            v[i].append(c.pop(0))
            u[i].append(z.pop(0))

#    print Numeric.shape(x),Numeric.shape(y),Numeric.shape(z),Numeric.shape(v),Numeric.shape(c),Numeric.shape(u),q

    fig = pylab.figure(figsize=(8,6),facecolor='w',edgecolor='k')
    pylab.title('Bader Surface',fontsize=18)
    pylab.contourf(x,y,Numeric.transpose(v),10,linewidth=0.05)
    pylab.show()

if __name__ == '__main__':
    main()
