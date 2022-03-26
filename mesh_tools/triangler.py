while True:
	try:
		line = input()
	except EOFError:
		break
	
	count = line.count(' ')
	if line[0] != 'f' or count == 3:
		print(line)
		continue
	elif count == 4:
		sp = line.split(' ')
		print(' '.join(sp[:4]))
		print(' '.join([sp[0], sp[1], sp[2], sp[4]]))
	else:
		print("bad:")
		print(line)
		break
