import os
import glob
import sys
import numpy as np
import matplotlib.pyplot as plt

filenames = glob.glob(os.getcwd() + '/*.dat')

all_data = []
# Load data
for filename in filenames:
    f = open(filename)
    data = []
    for line in f:
        data.append(line.split()[3])
    all_data.append(data)

labels = [0.4, 0.5, 0.6, 0.7, 0.8]

# Plot
for i,data in enumerate(all_data):
    label = "$d = " + str(labels[i]) + "$"
    plt.plot(data, label=label)

plt.xlabel('Timesteps')
plt.ylabel('Tumor Radius')
plt.title('Tumor Size Over Time, Varying Death Rates')
plt.legend(loc=2)
plt.show()