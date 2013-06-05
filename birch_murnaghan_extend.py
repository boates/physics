#!/usr/bin/env python

import Numeric,sys,math

def main():

    try:
        BM, rs = [], []
        BM.append(float(sys.argv[1]))
        BM.append(float(sys.argv[2]))
        BM.append(float(sys.argv[3]))
        rs.append(float(sys.argv[4]))
        rs.append(float(sys.argv[5]))
    except:
        print '\nusage: '+sys.argv[0]+' BM_param[0], BM_param[1], BM_param[2], min_rs, max_rs\n'
        sys.exit(0)

    RS = Numeric.arange(rs[0],rs[1],0.001,typecode=Numeric.Float)

    V = RS**3*4./3.*math.pi*320.0

    P = (3*BM[0]/2.0)*( (BM[1]/V)**(7./3.) - (BM[1]/V)**(5./3.) )*( 1+3/4.*(BM[2]-4.)*((BM[1]/V)**(2./3.)-1) )

#    print P, V

    out = open('BM_range.dat','w')
    for i in range(len(P)):
        out.write(str(RS[i])+' '+str(P[i])+'\n')
    out.close()


if __name__ == '__main__':
    main()
