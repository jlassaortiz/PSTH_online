#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 13:51:54 2020

@author: javi_lassaortiz
"""

from load_intan_rhd_format.load_intan_rhd_format import read_data

# directorio = input('directorio: ')
directorio = '/Users/javi_lassaortiz/Documents/LSD/Finch dormido/Test procesamiento online'

datos = read_data(directorio + '/info.rhd')
