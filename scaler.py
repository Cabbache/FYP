#!/usr/bin/python3

import sys

if len(sys.argv) != 3:
	print("Usage: " + sys.argv[0] + " [obj file] [multiplier]")
	sys.exit(1)

filename = sys.argv[1]

try:
	scale = float(sys.argv[2])
except ValueError:
	print("multiplier not a float");
	sys.exit(1)

with open(filename) as objfile:
	lines = [line.rstrip() for line in objfile]

for line in lines:
	if 'f' in line:
		print(line)
		continue
	print('v ' + ' '.join([str(round(scale*float(x), 7)) for x in line.split(' ')[1:]]))
