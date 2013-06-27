import pybam

f = pybam.Fasta("/share/lustre/jknaggs/Homo_sapiens_assembly19.fasta")
f.load("1")
print f.nmer(1000234, 5)
print f.nmer(1000243, 5)
