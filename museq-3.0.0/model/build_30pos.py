import sys
sqval=str(30)
label_set=set(['SOMATIC', 'GERMLINE', 'WILDTYPE'])
label=None
name=sys.argv[1].split('.')[0]

for line in open(sys.argv[1]):
    l=line.split()
    f=l[0]
    #print f
    #print f[0]
    if ((f[0]=="#")and (f[1]=="#")and (f[2]=="/")):
        reference=l[3].split(':')[1]
        normal=l[1].split(':')[1]
        tumour=l[2].split(':')[1]
        with open(name+".pos", 'w') as myfile:
            myfile.write('#'+' '+'normal'+' '+normal+'\n')
            myfile.write('#'+' '+'tumour'+' '+tumour+'\n')
            myfile.write('#'+' '+'reference'+' '+reference+'\n')
    if (f[0]!="#"):
        chromosome=l[0]
        position=l[1]
        label=l[3]
    print label    
    if label in label_set:
        with open(name+".pos", 'a+') as myfile:
            myfile.write(str(chromosome)+" "+str(position)+" "+str(label)+'\n')
           
