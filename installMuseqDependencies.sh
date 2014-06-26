#!/bin/bash 
## This is meant for linux
if [ $# -ne 1 ]; then
    echo "usage ./installMe [/full/path/to/dest/directory]"
    exit 1
fi
 
## initialize
cd $1 
#mkdir Python
mkdir PythonDep
pPath=$1
bPath=${pPath}/bin
lPath=${pPath}/lib

## helper functions
function errexit(){
    echo "Oops! couldn't download $1 :("
    exit 1
}

function pause(){
   read -p "Press enter to continue ..."
}

function label(){
    echo "=================================" >> ${pPath}/install.log
    echo "log of $1 " >> ${pPath}/install.log
    echo "---------------------------------" >> ${pPath}/install.log
}
function gettar(){
    wget $1
    modname=$(basename $1)
    if [ -f $modname ]; then
	echo "extracting ..."
	tar xvf $modname 1> /dev/null 2>> ${pPath}/install.log
	rm $modname
	dirname=$(basename $(basename $(basename $1 .tar.gz) .tar.bz2) .tgz)
	cd $dirname
    else 
	errexit $modname
    fi
}

function getinst(){
    echo "installing $(basename `pwd`) ..."
    ${bPath}/python2.7 setup.py install 1> /dev/null 2>> ${pPath}/install.log   
}

function getmake(){
    echo "configuring ..."
    ./configure --prefix=$1 1>/dev/null 2>> ${pPath}/install.log
    echo "making ..."
    make 1> /dev/null 2>> ${pPath}/install.log
    echo "installing ..."
    make install 1> /dev/null 2>> ${pPath}/install.log
    echo "cleaning ..."
    make clean 1> /dev/null 2>> ${pPath}/install.log
}

label "install python"
cd ${pPath}
gettar http://www.python.org/ftp/python/2.7.5/Python-2.7.5.tar.bz2
getmake ${pPath}

pause
label "install atlas+lapack"
cd ${pPath}
cd PythonDep
gettar http://www.netlib.org/lapack/lapack-3.4.2.tgz
gettar http://users.wfu.edu/cottrell/lapack/lapack-3.4.0-autoconf.tar.gz
getmake ${lPath}
export ATLAS=${lPath}/lib


pause
label "install nympy"
cd ${pPath}
cd PythonDep
gettar http://sourceforge.net/projects/numpy/files/NumPy/1.7.1/numpy-1.7.1.tar.gz 
getinst

pause
label "edit site.cfg in numpy module in order for scipy to install"
cd ${pPath}
echo $(pwd)
pause
cp PythonDep/numpy-1.7.1/site.cfg.example ${lPath}/python2.7/site-packages/numpy/distutils/site.cfg
echo "
[DEFAULT]
library_dirs = ${lPath}:${lPath}/lib

[blas_opt]
libraries = ${lPath}/lib

[lapack_opt]
libraries = ${lPath}/lib

[blas_opt]
libraries = ${lPath}/lib

[lapack_opt]
libraries = ${lPath}/lib

[lapack]
libraries = ${lPath}/lib
" >> ${lPath}/python2.7/site-packages/numpy/distutils/site.cfg 2>> ${pPath}/install.log

pause
label "install scipy"
cd ${pPath}
cd PythonDep
gettar http://sourceforge.net/projects/scipy/files/scipy/0.12.0/scipy-0.12.0.tar.gz
getinst

pause
label "install freetype2"
cd ${pPath}
cd PythonDep
wget http://download.savannah.gnu.org/releases/freetype/freetype-2.4.12.tar.gz
if [ -f freetype-2.4.12.tar.gz ]; then
    tar xvfz freetype-2.4.12.tar.gz 1> /dev/null 2>>${pPath}/install.log
    rm freetype-2.4.12.tar.gz
    cd freetype-2.4.12/
    echo "installing freetype2 ..."
    ./configure --prefix=${lPath} 1> /dev/null 2>>${pPath}/install.log
    make 1> /dev/null 2>>${pPath}/install.log
    mkdir ${lPath}/include/freetype2/freetype/internal
    pause
    make install 1> /dev/null 2>>${pPath}/install.log
    PKG_CONFIG_PATH=${lPath}/lib/pkgconfig:$PKG_CONFIG_PATH 
    export PKG_CONFIG_PATH
else 
    errexit freetype2
fi

pause
label "install matplotlib"
cd ${pPath}
cd PythonDep
gettar https://downloads.sourceforge.net/project/matplotlib/matplotlib/matplotlib-1.2.1/matplotlib-1.2.1.tar.gz
getinst

pause
label "isntall scikit-learn"
cd ${pPath}
cd PythonDep
gettar http://sourceforge.net/projects/scikit-learn/files/scikit-learn-0.13.1.tar.gz
getinst

pause
label "install ruffus"
export PATH=${bPath}:$PATH
cd ${pPath}
cd PythonDep
wget --no-check-certificate https://pypi.python.org/packages/2.7/s/setuptools/setuptools-0.6c11-py2.7.egg 
chmod +x setuptools-0.6c11-py2.7.egg
sh setuptools-0.6c11-py2.7.egg
cd ${lPath}/python2.7/site-packages/
gettar https://ruffus.googlecode.com/files/ruffus-2.2.tar.gz 
getinst 


pause
label "remove the garbage"
# cd ${pPath}
# echo "clean the mess?"
# pause
# rm -rf ../PythonDep 2>> ${pPath}/install.log

##=====================================
## you may want to export the following
## env variables or add them to .bashrc
##-------------------------------------
#PYTHONPATH=$PYTHONPATH:${lPath}:${lPath}/lib
#PATH=$PATH:${bPath}
#LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${lPath}/lib
#export PYTHONPAH 
#export PATH
#export LD_LIBRARY_PATH



