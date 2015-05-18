import random

def cross_types(geom_opts, lat_opts):
	while True:
		method = [random.choice(geom_opts),random.choice(lat_opts)]
		if method == [[1,1]]:
			break
		elif method  == [[2,2]]:
			break
		else:
			return method

