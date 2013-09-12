from distutils.core import setup, Extension

samtools = ['samtools/kstring.c','samtools/sam_header.c','samtools/razf.c',
            'samtools/bgzf.c','samtools/kaln.c','samtools/kprobaln.c',
            'samtools/bam.c','samtools/faidx.c','samtools/sam.c',
            'samtools/bam_pileup.c','samtools/bam_index.c',
            'samtools/bam_aux.c','samtools/bam_import.c',
            'samtools/bam_md.c','samtools/errmod.c']

old_pybam = Extension('old_pybam',
        sources = samtools + ['src/base.c','src/fasta.c','src/old_pybam.c','src/old_pybam_init.c'],
        library_dirs = [],
        include_dirs = ['samtools/','samtools/bcftools','include/samtools'],
        libraries = ['m', 'z'])

setup (name = 'old_pybam',
        version = '0.1',
        description = 'old_pybam',
        ext_modules = [old_pybam])
