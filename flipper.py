while True:
	try:
		line = input()
	except EOFError:
		break
	
	count = line.count(' ')
	if line[0] != 'v':
		print(line)
		continue
	sp = line.split(' ')
	print(' '.join([sp[0], sp[2], sp[1], sp[3]]))
