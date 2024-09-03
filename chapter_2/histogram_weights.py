# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:41:58 2023

@author: finkt
"""

import json
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import lognorm

f = open('congress_network_data.json')
data = json.load(f)

inList = data[0]['inList']
inWeight = data[0]['inWeight']
outList = data[0]['outList']
outWeight = data[0]['outWeight']
usernameList = data[0]['usernameList']

all_weights = []
for single_node_weights in outWeight:
    all_weights.extend(single_node_weights)

n, bins, patches = plt.hist(all_weights, bins=100, density=True)

# best fit of data for lognorm distribution
s, loc, scale=lognorm.fit(all_weights,floc=0.0)

x=np.linspace(0,0.14,10000)
y=lognorm.pdf(x,s, loc=loc, scale=scale)
plt.plot(x,y,label='lognorm',linewidth=4)

plt.legend()