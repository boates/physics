#!/usr/bin/env python

def testsort(a,b):
    """
    a is sorted increasingly
    b is sorted decreasingly
    c is the combined sorted list
    """
    c = []
    n = len(a) + len(b)
    a.append(1000)
    b = [1000]+b
    i, j = 0, len(b)-1
    while len(c) < n:
        if a[i] <= b[j]:
            c.append(a[i])
            i+=1
        else:
            c.append(b[j])
            j-=1

    return c


def main():
    a = [1,2,3,4,5,6,7]
    b = [9,7,5,4,3,1]
    s = a+b
    s.sort()
    print s
    
    c = testsort(a,b)
    print c


if __name__ == '__main__':
    main()
