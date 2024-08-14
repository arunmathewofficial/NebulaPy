#!/bin/bash
#
# Authors: Jonathan Mackey, Harpreet Dhanoa, others?
#
# - 2012-2020: ongoing development to add options for new machines.
# - 2023.10.03 JM: updates for new MacOS
# 
# MacOS: use homebrew, and install:
# brew install readline python3 numpy scipy gcc git 
# 

mkdir include
mkdir bin
mkdir lib
mkdir lib64
mkdir share
mkdir cmake

NCORES=4
export CC=gcc
export CXX=g++
export FC=gfortran
SHARED=YES
HDF5=YES
HDF5_LIBS="/usr/include/hdf5/serial,/usr/lib/x86_64-linux-gnu/hdf5/serial"
export PYTHON=/usr/bin/python3 

COMPILE_SILO=yes

export WGET='wget'

case $HOSTNAME in
  login[0-9].kay.ichec.ie)
    echo "Compiling on KAY/ICHEC"
    source /usr/share/Modules/init/bash
    module purge
    module load cmake3
    module load gcc
#module load cmake3/3.12.3
#module load python py/intel
#module load python numpy
    HDF5=NO
    SHARED=YES
    module list
    ;;
esac

DDD=`uname -a | grep "Darwin"`
if [ ! -z "$DDD" ]; then
  #export CXX=/usr/local/bin/g++-13
  #export CC=/usr/local/bin/gcc-13
  #export FC=/usr/local/bin/gfortran-13
  echo "***** COMPILING WITH OS-X: host ${HOST}: COMPILERS ARE $CC $CXX "  
  MAKE_UNAME=OSX
  #NCORES=1
  path=`pwd`
  export PYTHON=`which python3`
fi


export NCORES
CURDIR=`pwd`

##################################
##########     SILO     ##########
##################################
if [ "$COMPILE_SILO" == "yes" ]
then

  echo "********************************"
  echo "*** INSTALLING SILO LIBRARY ****"
  echo "********************************"
#################################
# Change these for new versions:
#  FILE=silo-4.10.2-bsd.tgz
#  SRC_DIR=silo-4.10.2-bsd
  VER=4.11.1
  REMOTE_URL=https://github.com/LLNL/Silo/releases/download/${VER}/silo-${VER}-bsd-smalltest.tar.xz
  FILE=silo-${VER}-bsd-smalltest.tar.xz
  SRC_DIR=silo-${VER}-bsd
#################################


  if [ -e $FILE ]; then
          echo "*** File exists, no need to download ***"
  else 
          echo "***** File does not exist ******"
          echo "********************************"
          echo "*** DOWNLOADING SILO LIBRARY ***"
          echo "********************************"
          $WGET --no-check-certificate $REMOTE_URL -O $FILE
          if [ $? != 0 ]; then
            echo "** failed to download with wget, trying curl instead **"
            rm $FILE
            curl --location --remote-name $REMOTE_URL -o $FILE
          fi
          # check it downloaded.
          if [ ! -f $FILE ]; then
            echo "File not found! : $FILE"
            echo "Download of Silo Library Failed... quitting"
            exit
          fi
  fi 
#echo "********************************"
#echo "*** EXTRACTING SILO LIBRARY ***"
#echo "********************************"
  tar --xz -xf $FILE
#echo "********************************"
#echo "*** RUNNING CONFIGURE ***"
#echo "********************************"
  BASE_PATH=`pwd`
  mkdir -p ${HOME}/.local/
  mkdir -p ${HOME}/.local/silo
  echo "***Path = $BASE_PATH ***"
  cd $SRC_DIR
  make clean

  if [[ MAKE_UNAME == "OSX" ]]
    then
    ./configure --prefix=$HOME/.local/silo \
      --enable-browser --with-readline=no \
      --disable-fortran \
      --disable-silex --disable-fpzip \
      --enable-pythonmodule --enable-shared
  else
    ./configure --prefix=$HOME/.local/silo \
      --enable-browser \
      --disable-fortran \
      --disable-silex \
      --enable-pythonmodule --enable-shared
  fi


  echo "********************************"
  echo "*** RUNNING MAKE ***"
  echo "********************************"
  make -j $NCORES
  echo "********************************"
  echo "*** INSTALLING SILO LIBRARY ***"
  echo "********************************"
  make install
  cd $CURDIR
  echo "********************************"
  echo "*** FINISHED! ***"
  echo "********************************"
fi


cd $CURDIR
