# issue installing volcanofinder on lisc - fails in the make step

#-----------------------------------------------------------------------------------------------------------------------
tar -xzvf volcanofinder_v1.0.tar.gz
cd volcanofinder_v1.0
make


#attempt 1 error messages
gcc -c VolcanoFinder.c -O3 -Wall
gcc -c freq.c -O3 -Wall
gcc -c factorials.c -O3 -Wall
gcc -c bfgs.c -O3 -Wall
bfgs.c: In function ‘setulb_’:
bfgs.c:252:46: warning: variable ‘l3’ set but not used [-Wunused-but-set-variable]
  252 |     static integer lsnd, lsgo, lygo, l1, l2, l3, ld, lr, lt;
      |                                              ^~
bfgs.c:252:42: warning: variable ‘l2’ set but not used [-Wunused-but-set-variable]
  252 |     static integer lsnd, lsgo, lygo, l1, l2, l3, ld, lr, lt;
      |                                          ^~
bfgs.c:252:38: warning: variable ‘l1’ set but not used [-Wunused-but-set-variable]
  252 |     static integer lsnd, lsgo, lygo, l1, l2, l3, ld, lr, lt;
      |                                      ^~
gcc -c sort.c -O3 -Wall
gcc -c my_rand.c -O3 -Wall
gcc -o VolcanoFinder VolcanoFinder.o -O3 -Wall freq.o factorials.o bfgs.o sort.o my_rand.o -lm
/usr/bin/ld: freq.o:(.bss+0x10): multiple definition of `data'; VolcanoFinder.o:(.bss+0xd38): first defined here
/usr/bin/ld: freq.o:(.bss+0x0): multiple definition of `data_bvalue'; VolcanoFinder.o:(.bss+0xd28): first defined here
/usr/bin/ld: freq.o:(.bss+0x8): multiple definition of `data_rec'; VolcanoFinder.o:(.bss+0xd30): first defined here
collect2: error: ld returned 1 exit status
make: *** [Makefile:6: VolcanoFinder] Error 1

#-----------------------------------------------------------------------------------------------------------------------



#cluster support: make an issue on the programme github or you stepwise go back with the GCC versions, maybe back to GCC 4.8, which was the default in Redhat Linux 8.

#try to go back to gcc 4.8

# check alternative versions of gcc on cluster
module avail g* # does not contain 'gcc'
module avail gcc # does not give output
applocate gcc # does not have older versions of gcc
# module not avail - check within conda

#-----------------------------------------------------------------------------------------------------------------------

# search within conda - is giving newer versions of the package
(base) [pawar@login01 hg19_1000g]$ conda search gcc
Loading channels: done
# Name                       Version           Build  Channel             
gcc                            8.5.0      h143be6b_1  conda-forge         

  (base) [pawar@login01 hg19_1000g]$ conda search "gcc=4.8"
Loading channels: done
No match found for: gcc=4.8. Search: *gcc*=4.8
# Name                       Version           Build  Channel             
libgcc                         4.8.4               1  conda-forge         
libgcc                         4.8.5               2  conda-forge 
#try to install older version

# create temp env - where install lower version of gcc
conda create -n gcc_4.8.5
conda activate gcc_4.8.5
conda install -c conda-forge/label/cf201901 gcc=4.8.5

Solving environment: failed

PackagesNotFoundError: The following packages are not available from current channels:

  - gcc=4.8.5*

conda install -c conda-forge gcc=4.8.5

Solving environment: failed

PackagesNotFoundError: The following packages are not available from current channels:

  - gcc=4.8.5*

#-----------------------------------------------------------------------------------------------------------------------
#try to install  gcc=4.8.5 using the tar version from the website

(gcc_4.8.5) [pawar@login01 harvi]$ conda install https://anaconda.org/conda-forge/gcc/4.8.5/download/osx-64/gcc-4.8.5-8.tar.bz2

Downloading and Extracting Packages:
                                                                                                                                                
[Errno 122] Disk quota exceeded: '/lisc/scratch/admixlab/harvi/miniforge3/pkgs/gcc-4.8.5-8/info/repodata_record.json'

#-----------------------------------------------------------------------------------------------------------------------
#try to install using a docker image for this version

(gcc_4.8.5) [pawar@login01 ~]$ docker pull gcc:4.8.5
Emulate Docker CLI using podman. Create /etc/containers/nodocker to quiet msg.
Error: creating runtime static files directory "/.local/share/containers/storage/libpod": mkdir /.local: permission denied

gcc_4.8.5) [pawar@login01 ~]$ docker info
Emulate Docker CLI using podman. Create /etc/containers/nodocker to quiet msg.
Error: creating runtime static files directory "/.local/share/containers/storage/libpod": mkdir /.local: permission denied
#do not have permission to use docker
