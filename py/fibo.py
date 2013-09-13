# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:42:12 2013

@author: adamwolf
"""

def fib(n):    # write Fibonacci series up to n
     "Print a Fibonacci series up to n"
     a, b = 0, 1
     while b < n:
        print b,
        a, b = b, a+b
        
if __name__ == "__main__":
    import sys
    fib(int(sys.argv[1]))        