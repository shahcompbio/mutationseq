python setup.py build &&
cp build/lib.linux-x86_64-2.7/pybam.so test &&
cd test &&
python nmer.py
cd ..
