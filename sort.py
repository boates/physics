#!/usr/bin/env python
"""
sort.py
Author: Brian Boates
"""
import random

def bubblesort(x):
    N = len(x)
    if N <= 1: return x
    swapped = True
    while swapped == True:
        swapped = False
        for i in range(N-1):
            if x[i+1] < x[i]:
                x[i], x[i+1] = x[i+1], x[i]
                swapped = True
    return x

def quicksort(x):
    N = len(x)
    if N <= 1: return x
    pivot = x.pop( N/2 )
    less, greater = [], []
    for v in x:
        if v <= pivot: less.append(v)
        else: greater.append(v)
    return quicksort(less) + [pivot] + quicksort(greater)

def mergesort(x):
    # broken
    N = len(x)
    if N <= 1: return x
    left = [x[i] for i in range(N/2)]
    right = [x[i] for i in range(N/2,N)]
    left = mergesort(left)
    right = mergesort(right)
    result = []
    nL, nR = len(left), len(right)
    while nL > 0 or nR > 0:
        if nL > 0 and nR > 0:
            if left[0] <= right[0]: result.append(left.pop(0))
            else: result.append(right.pop(0))
        elif nL > 0: result.append(left.pop(0))
        elif nR > 0: result.append(right.pop(0))
    return result

x = range(10)
random.shuffle(x)
print x
#x = [4,2,11,5,23,8,0,3,1,15,6,5]
print bubblesort(x)
print quicksort(x)
#print mergesort(x)
