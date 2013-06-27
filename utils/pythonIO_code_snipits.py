'''
Created on 6 Jun 2013

@author: musselle
'''
import os 
import sys 

import numpy as np 
import matplotlib.pyplot as plt

import argparse


parser = argparse.ArgumentParser(description='Procedure to write CDHIT.clstr file to FastQ')







# raw_input() asks the user for a string of data (ended with a newline), and simply returns the string. 

# input() uses raw_input to read a string of data, and then attempts to evaluate it as if it were a Python program
x = None
while not x:
    try:
        x = int(raw_input())
    except ValueError:
        print 'Invalid Number'

# $ python program.py hello there programmer!
sys.argv[0] # --> program.py
sys.argv[1] # --> hello
sys.argv[2] # --> there







if __name__ == '__main__':
    pass