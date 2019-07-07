import os
for t in [ 1, 2, 4, 8, 16 ]:
	for s in [ 100, 1000, 10000, 40000, 80000, 200000, 400000, 600000, 800000, 1000000 ]:
		cmd = "g++ -DNUMTRIALS=%d -DNUMT=%d proj_1.cpp -o prog -lm -fopenmp" % ( s, t )
		os.system( cmd )
		cmd = "./prog"
		os.system( cmd )