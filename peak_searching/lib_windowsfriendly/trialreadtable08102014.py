# -*- coding: utf-8 -*-
"""
Created on Wed Oct 08 15:15:23 2014

@author: lina2532
"""
import numpy as np

damagefiletrial = open("distancefile_trial.txt", "r")

trialdamage = np.genfromtxt("distancefile_trial.txt", names=True, delimiter='\t', dtype=None)

print trialdamage['basetype']

