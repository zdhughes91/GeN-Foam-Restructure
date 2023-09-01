
### IMPORTS

import sys
import matplotlib.pyplot as plt

### MAIN

fig, (ax1) = plt.subplots(1, sharex=True)

ax1.set_ylabel("power (W)")
ax1.set_xlabel("time(s)")

file=open('transient_diff/GeN-Foam.dat', "r")


list=[];
time=[];
power=[];
file.readline()
file.readline()

for line in file:
	newline=line.strip()
	list+=[newline.split(';')]

for i in range(0,len(list)):
	time.append(float(list[i][0]))
	power.append(float(list[i][2]))

t0 = time[0]
for i in range(len(time)) :
	time[i] = time[i] - t0  

ax1.plot(time, power, label='transient_diff/GeN-Foam.dat')

ax1.legend()
ax1.grid(True)
plt.show()

