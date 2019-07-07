import os

k = 1024
BLOCKSIZEs = [16, 32, 64]
SIZE  = [512*k]

for arr in SIZE:
	for BLOCKSIZE in BLOCKSIZEs:
		cmd = '/usr/local/apps/cuda/cuda-9.2/bin/nvcc -DBLOCKSIZE=%d -DSIZE=%d -o arrayMul  arrayMul.cu' % (BLOCKSIZE, arr)
		os.system( cmd )
		cmd = "./arrayMul"
		os.system( cmd )
		cmd = "rm -f arrayMul"
		os.system( cmd )