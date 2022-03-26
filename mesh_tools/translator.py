#!/usr/bin/python3

import sys

if len(sys.argv) != 4:
	print("Usage: " + sys.argv[0] + " [x] [y] [z]")
	sys.exit(1)

try:
	t_x = float(sys.argv[1])
	t_y = float(sys.argv[2])
	t_z = float(sys.argv[3])
except ValueError:
	print("x,y, or z are not a float");
	sys.exit(1)

while True:
	try:
		line = input()
	except EOFError:
		break

	if 'f' in line:
		print(line)
		continue
	bits = line.split(' ')[1:]
	print("v " + str(round(float(bits[0])+t_x, 7)) + " " + str(round(float(bits[1])+t_y, 7)) + " " + str(round(float(bits[2])+t_z, 7)))
