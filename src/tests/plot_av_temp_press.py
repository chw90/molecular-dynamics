#!/usr/bin/python3

import yaml
from matplotlib import pyplot as plt
import numpy as np

dt = 0.002

file = open('av_temp_press.yaml', 'rb')
data = np.asarray(yaml.safe_load(file)['data'])
data_time = data[:,0] *dt
data_temp = data[:,1]
data_press = data[:,2]

fig, axs = plt.subplots(2)
fig.suptitle('Convergence Plot')

axs[0].plot(data_time, data_temp)
axs[0].set_ylabel('Temperature in K')
axs[1].plot(data_time, data_press)
axs[1].set_ylabel('Pressure in Pa')
axs[1].set_xlabel('Time in s')

plt.savefig("av_temp_press.pdf")
plt.show()
