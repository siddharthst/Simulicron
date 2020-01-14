#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
module_path = os.path.abspath(os.path.join('../src/simulicronalpha/'))
if module_path not in sys.path:
    sys.path.append(module_path)


# In[ ]:


# Imports
import inspect
import random
import tempfile
import uuid
from functools import wraps
import collections
import seaborn as sns
import matplotlib.pyplot as plt
import coloredlogs
import logging
import numpy as np
import pandas as pd
from coolname import generate_slug
import warnings
import cProfile

warnings.filterwarnings("ignore")

# Logger
logger = logging.getLogger(__name__)
coloredlogs.install(level="DEBUG")


# In[ ]:


import RunSim as init


# In[ ]:


k = init.RunSim(tcount=2, insize=100)


# In[ ]:


z = k.createSim()


# In[ ]:


k.resetFitness()


# In[ ]:


k.runSimulation()

