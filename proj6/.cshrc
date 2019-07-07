setenv INTEL_LICENSE_FILE 28518@linlic.engr.oregonstate.edu
setenv ICCPATH /nfs/guille/a2/rh80apps/intel/studio.2013-sp1/composer_xe_2015/bin/
set path=( $path $ICCPATH )
source /nfs/guille/a2/rh80apps/intel/studio.2013-sp1/bin/iccvars.csh intel64
setenv CUDA_PATH /usr/local/apps/cuda/cuda-9.2
setenv LD_LIBRARY_PATH $CUDA_PATH/lib64:$LD_LIBRARY_PATH
set path = ( $path $CUDA_PATH/bin )