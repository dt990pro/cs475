import os

arr = [ 1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000, 50000000 ]

for i in arr:
	print(i)

for s in arr:
	cmd = "g++ -DARR_SIZE=%d simd.p4.cpp proj_4.cpp -o prog -lm -fopenmp" % ( s )
	os.system( cmd )
	cmd = "./prog"
	os.system( cmd )
	cmd = "rm -f prog"
	os.system( cmd )