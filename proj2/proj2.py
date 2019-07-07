import os
for t in [ 1, 2, 4 ]:
	for s in [ 2, 4, 8, 16, 32, 64, 128, 256 ]:
		cmd = "g++ -DNUMNODES=%d -DNUMT=%d proj2.cpp -lm -fopenmp" % ( s, t )
		os.system( cmd )
		cmd = "a.out"
		os.system( cmd )