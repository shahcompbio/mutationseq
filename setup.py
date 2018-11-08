import os
import sys
from setuptools import setup, find_packages
import versioneer
import museq
import museq.pickle_models
import museq.setup

try:
    boost_source_arg = filter(lambda a: a.startswith('--boost_source='), sys.argv)[0]
    boost_source = os.path.expanduser(boost_source_arg.split('=')[1])
    sys.argv.remove(boost_source_arg)
except IndexError:
    sys.stderr.write('Please specify the boost source directory using --boost_source=dir\n')
    sys.exit(1)

currentdir = os.getcwd()
os.chdir("museq")
#museq.setup.compile(boost_source)
os.chdir(currentdir)


#museq.pickle_models.setup_museq_models()

setup(
    name='mutationseq',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='mutationseq',
    author='Jafar Taghiyar Renani',
    author_email='dgrewal@bccrc.ca',
    entry_points={'console_scripts': ['museq = museq.classify:main', 'museq_het = museq.preprocess:main']},
    package_data={'':['*.pickle', '*.so', "*.config"]}
)


