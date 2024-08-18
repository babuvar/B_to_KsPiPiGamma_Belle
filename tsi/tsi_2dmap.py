#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np

#clean
os.system('rm -f output.txt')

for Sin in np.arange(-0.6, 0.7, 0.1):
	for Ain in np.arange(-0.3, 0.4, 0.1):

		os.system('./calc_observed2true init_observed2true.txt %s %s'%(Sin, Ain))
