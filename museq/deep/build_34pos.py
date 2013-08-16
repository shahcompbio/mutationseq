import sys

data=set()
label_set=set(['SOMATIC', 'GERMLINE', 'WILDTYPE'])

for line in open(sys.argv[1]):
    l=line.split()
    label=l[0]
    if label in label_set:
    	ll=l[1].split('_')
    	d=ll[1]
    	name="34_"+str(d)
    	chromosome=ll[2]
    	position=ll[3]
    	if d in data:
            with open(name+".pos", 'a+') as myfile:
            	myfile.write(str(chromosome)+" "+str(position)+" "+str(label)+'\n')
    	else:
            data.add(d)
            with open(name+".pos", "w") as myfile:
            	myfile.write("# normal "+"/share/lustre/archive/"+"?/"+"illumina_exoncapture/"+'\n')
                myfile.write("# tumour "+"/share/lustre/archive/"+d+"/illumina_exoncapture/"+'\n')
                myfile.write("# reference "+"/share/data/genomes/human_all.fasta"+'\n')
            	myfile.write(str(chromosome)+" "+str(position)+" "+str(label)+'\n')

        
            
    