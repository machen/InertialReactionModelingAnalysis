import pandas as pd
import numpy as np
import re, os
import matplotlib.pyplot as plt
import seaborn as sns

"""Script is meant to do slicing operations on data sets
and calculate average values in the various slices entered. It should take in
the chemdata format or flow data format as defined by COMSOL and then given a list
of slice locations and IDs, calculate average values that can then be used for follow up analysis"""

