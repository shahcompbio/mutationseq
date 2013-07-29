import sys

data = {}
for line in open(sys.argv[1]):
    line = line.strip().split()
   
    t = (line[2], line[3], line[4])
    pair = (line[0], line[1])
    if pair not in data:
        data[pair] = []
    
    data[pair].append(t)

for pair in data:
    case = pair[1].split('.')[0]
    case = case.split('/')[-1]
    case = case.split('_')[0]
    if not case:
        continue
    fname = sys.argv[2] + '_' + case + '.pos'
    fd = open(fname, 'w')
    print >> fd, "# normal", pair[0]
    print >> fd, "# tumour", pair[1]
    print >> fd, "# reference /share/lustre/reference/genomes/human_all.fasta"
    for t in data[pair]:
        print >> fd, t[0], t[1], t[2]
    fd.close()
