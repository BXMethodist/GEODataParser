import os

def runFastqdump(list):
    f = open(list, 'r')
    ids = [x.strip() for x in f.readlines()]
    f.close()

    for link in ids:
        id = link[link.rfind('/')+1:-9]
        print id
        pbs = open(id+".pbs", "w")
        pbs.write("#!/bin/bash\n")
        pbs.write("#PBS -r n\n")
        pbs.write("#PBS -N "+id+'\n')
        pbs.write("#PBS -q mediummem\n")
        pbs.write("#PBS -m e\n")
        pbs.write("#PBS -M bxia@houstonmethodist.org\n")
        pbs.write("#PBS -l walltime=96:00:00\n")
        pbs.write("#PBS -l nodes=1:ppn=8\n")
        pbs.write("#PBS -l pmem=16000mb\n")
        pbs.write("cd /home/tmhbxx3/archive/H3K4me3/ENCODE_with_input/FASTQ\n")
        pbs.write("module load python/2.7.11\n")
        pbs.write('wget '+link +'\n')
        pbs.write('gunzip ' + link[link.rfind('/')+1:])

        #pbs.write('cat ' + bowtie_path + name + '_1.fastq ' + bowtie_path + name + '_2.fastq >' + bowtie_path + name + '.fastq'+'\n')
        pbs.write('bowtie -p 8 -m 1 --chunkmbs 512 --best /archive/tmhkxc48/ref_data/hg19/bowtie/hg19 '
                  + id + ".fastq " + id + ".bowtie\n")
        pbs.close()
        os.system('qsub '+id+".pbs")
        break
    return

list = '/home/tmhbxx3/archive/H3K4me3/ENCODE_with_input/FASTQ/single.txt'


runFastqdump(list)