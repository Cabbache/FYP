#!/usr/bin/python3

import sys

if len(sys.argv) != 2:
	print("Usage: " + sys.argv[0] + " [uniform multiplier]")
	sys.exit(1)

try:
	scale = float(sys.argv[1])
except ValueError:
	print("multiplier not a float");
	sys.exit(1)

while True:
	try:
		line = input()
	except EOFError:
		break

	if 'f' in line:
		print(line)
		continue
	print('v ' + ' '.join([str(round(scale*float(x), 7)) for x in line.split(' ')[1:]]))
