from string import *
fp = open("TRANSFAC_and_JASPAR_PWMs.txt","r")
lines = fp.readlines()
fout = open("TRANSFAC_and_JASPAR_edgelist.txt","w")
for line in lines:
    data = line.strip("\r\n").split("\t")
    tf = data[0]
    targets = data[2:]
    for gene in targets:
        new_line = tf+"\t"+gene+"\n"
        fout.write(new_line)
fout.close()
fp.close()
