#!/usr/bin/python3

import sys

lines = []

while True:
	try:
		line = input()
		lines.append(line)
	except EOFError:
		break

def get_indexes(lines, line):
	indexes = []
	for x in range(len(lines)):
		if line == lines[x]:
			indexes.append(x)
	return indexes

verts = set([x for x in lines if "v" in x])
for vert in verts:
	indexes = get_indexes(lines, vert)
	if len(indexes) == 1:
		continue
	if len(indexes) > 2:
		print("Panic")
		sys.exit(1)

	#remove last one
	remove = indexes[-1]
	del lines[remove]

	#update faces
	for a in range(len(lines)):
		ln = lines[a]
		if "f" not in ln:
			continue
		bits = ln.split(" ")[1:]
		for b in range(len(bits)):
			bit = bits[b]
			if int(bit) > remove+1:
				bits[b] = str(int(bit)-1)
			elif int(bit) == remove+1:
				bits[b] = str(indexes[0]+1)
		lines[a] = "f " + " ".join(bits)

#remove identical faces
faces = "\n".join([
			("f " + " ".join(face)) for face in set(
				[
					frozenset(line.split(" ")[1:])
					for line in
					[x for x in lines if "f" in x]
				]
		)
	])

print("\n".join([x for x in lines if "v" in x]))
print(faces)
