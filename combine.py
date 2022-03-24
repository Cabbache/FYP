#!/usr/bin/python3

#This program combines two obj files trivially by
#simply appending the vertices and faces
#Duplicate faces are not removed, use trimmer.py for this

import sys

if len(sys.argv) != 3:
	print("Usage: " + sys.argv[0] + " [obj1] [obj2]")
	sys.exit(1)

filename1 = sys.argv[1]
filename2 = sys.argv[2]

with open(filename1, "r") as f1:
	f1lines = [x.strip() for x in f1.readlines()]
with open(filename2, "r") as f2:
	f2lines = [x.strip() for x in f2.readlines()]

vertices1 = [x for x in f1lines if "v" in x]
vertices2 = [x for x in f2lines if "v" in x]
faces1 = [x for x in f1lines if "f" in x]
faces2 = []
for x in range(len(f2lines)):
	line = f2lines[x]
	if "v" in line:
		continue
	faces2.append("f " + " ".join([str(int(part) + len(vertices1)) for part in line.split(" ")[1:]]))

result = vertices1 + vertices2 + faces1 + faces2
for line in result:
	print(line)
