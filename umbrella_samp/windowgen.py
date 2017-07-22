""" This code generates a set of windows for site distances """

import re
import numpy as np

g = globals()
inps = open('wginput.txt')
inpn = np.genfromtxt('wginput.txt', delimiter='=')
line = inps.readlines()
w = np.array([])

for i in range(len(line)):
	if line[i][0] == "C":
		clen = inpn[i, 1]
	elif line[i][0] == "W":
		wgap = inpn[i, 1]
	elif line[i][0] == "N":
		n = inpn[i, 1]

ovr = (n * wgap - clen) / (n - 1)

output = open('wout.txt', 'w')

for i in range(int(n)):
	ww = [i * (wgap - ovr), (i + 1) * wgap - i * ovr]
	print(ww)
	print(i)
	output.write(str(ww)+'\n')
inps.close()
output.close()

"""for i in inps:
	if i[0]=="b'CLEN'":
		print('doonneee')

for i in range(int(n)):
		g["window"+str(i)] = [i*(wgap-ovr),(i+1)*wgap-i*ovr]
		print(g["window"+str(i)])
		w = np.append(w, ww, axis=0)

"""
