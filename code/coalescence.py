#!/usr/bin/env python3

__appname__ = '[coalescence.py]'
__author__ = 'Pablo Lechon (plechon@ucm.es)'
__version__ = '0.0.1'

## IMPORTS ##

import sys
import pandas as pd


## CONSTANTS ##


## FUNCTIONS ##

def main(argv):
    '''Main function'''
    #Load data from assembly simulations
    assembly_data = pd.read_csv('../data/simulation_results.csv')
    #Calculate cohesion of each community
    assembly_data['Cohesion'] = assembly_data['F'] - assembly_data['C']
    #Start coalescing communities from group l=0.5


    

    return 0

## CODE ##

if (__name__ == '__main__'):
    status = main(sys.argv)
    sys.exit(status)
     

