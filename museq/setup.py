import fnmatch
import os
import sys
import platform

from distutils.core import setup
from distutils.extension import Extension

try:
    boost_source_arg = filter(lambda a: a.startswith('--boost_source='), sys.argv)[0]
    boost_source = os.path.expanduser(boost_source_arg.split('=')[1])
    sys.argv.remove(boost_source_arg)
except IndexError:
    sys.stderr.write('Please specify the boost source directory using --boost_source=dir\n')
    sys.exit(1)

def get_filenames(directory, filter):
    for root, dirnames, filenames in os.walk(directory):
        for filename in fnmatch.filter(filenames, filter):
            yield os.path.join(root, filename)

boost_python_source = list(get_filenames(boost_source + '/libs/python/src', '*.cpp'))

bamtools_source = [
 'bamtools/src/api/BamAlignment.cpp',
 'bamtools/src/api/BamMultiReader.cpp',
 'bamtools/src/api/BamReader.cpp',
 'bamtools/src/api/BamWriter.cpp',
 'bamtools/src/api/SamHeader.cpp',
 'bamtools/src/api/SamProgram.cpp',
 'bamtools/src/api/SamProgramChain.cpp',
 'bamtools/src/api/SamReadGroup.cpp',
 'bamtools/src/api/SamReadGroupDictionary.cpp',
 'bamtools/src/api/SamSequence.cpp',
 'bamtools/src/api/SamSequenceDictionary.cpp',
 'bamtools/src/api/internal/bam/BamHeader_p.cpp',
 'bamtools/src/api/internal/bam/BamMultiReader_p.cpp',
 'bamtools/src/api/internal/bam/BamRandomAccessController_p.cpp',
 'bamtools/src/api/internal/bam/BamReader_p.cpp',
 'bamtools/src/api/internal/bam/BamWriter_p.cpp',
 'bamtools/src/api/internal/index/BamIndexFactory_p.cpp',
 'bamtools/src/api/internal/index/BamStandardIndex_p.cpp',
 'bamtools/src/api/internal/index/BamToolsIndex_p.cpp',
 'bamtools/src/api/internal/io/BamDeviceFactory_p.cpp',
 'bamtools/src/api/internal/io/BamFile_p.cpp',
 'bamtools/src/api/internal/io/BamFtp_p.cpp',
 'bamtools/src/api/internal/io/BamHttp_p.cpp',
 'bamtools/src/api/internal/io/BamPipe_p.cpp',
 'bamtools/src/api/internal/io/BgzfStream_p.cpp',
 'bamtools/src/api/internal/io/ByteArray_p.cpp',
 'bamtools/src/api/internal/io/HostAddress_p.cpp',
 'bamtools/src/api/internal/io/HostInfo_p.cpp',
 'bamtools/src/api/internal/io/HttpHeader_p.cpp',
 'bamtools/src/api/internal/io/ILocalIODevice_p.cpp',
 'bamtools/src/api/internal/io/RollingBuffer_p.cpp',
 'bamtools/src/api/internal/io/TcpSocket_p.cpp',
 'bamtools/src/api/internal/io/TcpSocketEngine_p.cpp',
 'bamtools/src/api/internal/io/TcpSocketEngine_unix_p.cpp',
 'bamtools/src/api/internal/sam/SamFormatParser_p.cpp',
 'bamtools/src/api/internal/sam/SamFormatPrinter_p.cpp',
 'bamtools/src/api/internal/sam/SamHeaderValidator_p.cpp',
 'bamtools/src/api/internal/utils/BamException_p.cpp',
 'bamtools/src/utils/bamtools_pileup_engine.cpp',
 'bamtools/src/utils/bamtools_fasta.cpp',
 'bamtools/src/utils/bamtools_utilities.cpp',
 ]

extra_compile_args = ['-Wno-unused-variable']
if platform.platform().startswith('Darwin'):
    extra_compile_args.append('-Wno-unneeded-internal-declaration')
    extra_compile_args.append('-Wno-unused-private-field')
    extra_compile_args.append('-Wno-mismatched-tags')

extra_link_args = []
if platform.platform().startswith('Darwin'):
    extra_link_args.append('-Wl,-no_compact_unwind')

ext_modules = [Extension('pybam',
                     ['src/pybam.cpp'] + boost_python_source + bamtools_source,
                     language='c++',
                     include_dirs=['./', './src', boost_source, './bamtools/src/'],
                     extra_compile_args=extra_compile_args,
                     libraries=['z'],
                     extra_link_args=extra_link_args
                     )]

setup(name='pybam', ext_modules=ext_modules)

