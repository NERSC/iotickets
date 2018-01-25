#!/usr/bin/env python

import argparse
import re
import json
import collections
import itertools
import sys
import os
import yaml
import stat

def main():
    parser = argparse.ArgumentParser(description='generate commands')
    parser.add_argument("-p", "--param",required=False, help="fasta style paramfile")
    parser.add_argument("-q", "--qsub", required=False,default=False,action='store_true',help="file")
    args = parser.parse_args()

    bindir = os.path.abspath(os.path.dirname(sys.argv[0]))
    wd = os.getcwd()
    array = []
    counter = 0

    param = ""
    f = None
    if (args.param):
        with open(args.param,'r') as f:
            param = f.read()
        f.close()
    else:
        param = get_param_example()

    header = ""
    for line in param.split("\n"):
        print line
        line = line.rstrip()
        result = re.search('^>(.*)$',line)
        if hasattr(result,'group'):
            if result.group(1):
                header = result.group(1)
                array.append([])
                counter = (counter + 1)
        else:
            array[counter-1].append(line)

#    sys.stderr.write(yaml.dump(array,default_flow_style=False))

    cmds = list(itertools.product(*array))

    
    dir = os.getcwd()
    mode = 0755
    #cmds bash
    counter = 1
    maq = int()
    trimq = int()
    for i in (range(1,len(cmds))):
        iter = "_%04d" % counter
        iterdir = dir + '/' + iter
        cmdbash = iterdir +  '/' + iter + '.bash'

        commandline = " ".join(cmds[i])
        result = re.search('^.*maq=(\d+).*$',commandline)
        if hasattr(result,'group'):
            if result.group(1):
                maq = result.group(1)
        result = re.search('^.*trimq=(\d+).*$',commandline)
        if hasattr(result,'group'):
            if result.group(1):
                trimq = result.group(1)
        
        maq = int(maq)
        trimq = int(trimq)
        if maq < trimq: #skip
            print str(maq) + " *********skipped******** " + str(trimq)
            continue

        
        counter=counter + 1
        os.mkdir(iterdir)        

#        ($maq=$_)=~s/^.*maq=(\d+).*$/$1/;($trimq=$_)=~s/^.*trimq=(\d+).*$/$1/;print if ($maq < $trimq)' | perl -pne 's/^(_\d+).*$/$1/' | xargs -i mv {} EXCLUDE/
        bash  = open(cmdbash,"w")
        cmd = "#!/usr/bin/env bash\nset -e\ncd " + iterdir + ";\n\n"
       # cmd = cmd + "touch stdout.log;touch stderr.log;\n"
        lines = commandline.split(';')
        for line in lines:
            cmd = cmd + line + " 1>> stdout.log 2>> stderr.log ;\n"
        bash.write(cmd + "\n")
        bash.close()        
        os.chmod(cmdbash, mode)        

    arraybash = dir + '/array.bash' 
    #array job bash
    bash  = open(arraybash,"w")
    cmd = "#!/usr/bin/env bash\nset -e;\n\n" 
    cmd = cmd + 'id=$[$SGE_TASK_ID - 0]' +  ";\n"
    cmd = cmd + 'printf -v id "_%04d" $id' + ";\n"
    cmd = cmd + 'cd ' + dir + '/${id}' + ";\n"
    cmd = cmd + 'bash ${id}.bash' + ";\n"
    bash.write(cmd + "\n")
    bash.close()
    os.chmod(arraybash, mode)        

    qsubbash = dir + '/qsub.bash'
    #qsub bash
    bash  = open(qsubbash,"w")
    cmd = "#!/usr/bin/env bash\n"
    cmd = cmd + '#SBATCH -N 1' + "\n"
    cmd = cmd + '#SBATCH -J cori_bigmem' + "\n"
    cmd = cmd + '#SBATCH -t 12:00:00' + "\n"
    cmd = cmd + '#SBATCH --mem=122G'
    cmd = cmd + '#SBATCH -M escori' + "\n"
    cmd = cmd + '#SBATCH --array=1-' + str(int(len(cmds)-1)) + "\n"
    cmd = cmd + "set -e;\n\n"
    cmd = cmd + "cd " + dir + ";\n"
    cmd = cmd + "bash " + arraybash + ";\n"
    bash.write( cmd + "\n")
    bash.close()
    os.chmod(qsubbash, mode)        

#    print json.dumps(array,indent=4)

'''
module load bbtools;rqcfilter.sh in=/global/dna/dm_archive/sdm/illumina/00/85/55/8555.2.104187.GGTAGC.fastq path=(output path) rna=f minlength=45 filterk=25 phix=t 
 qtrim=r maq=5 trimq=6 trimfragadapter=t hdist=1 filterk=25 
'''
 
def get_param_example():
    param_example = """>assemblers
/global/cscratch1/sd/aclum/spades_3.10.1/SPAdes-3.10.1-Linux/bin/spades.py -m 700 -o spades --only-assembler -k 33,55,77,99,127 --meta -t 32 --12 
>inputs
/test/data
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC04.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC06.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC13.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/low.biofuel.150.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/med.bog.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/mock.150.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/antarctic.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC04.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC06.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC13.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/low.biofuel.150.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/med.bog.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/mock.150.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/antarctic.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC04.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC06.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/MC13.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/low.biofuel.150.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/med.bog.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/mock.150.bfc21.fastq
/global/cscratch1/sd/aclum/metagenome_benchmarks/antarctic.bfc21.fastq"""

    return param_example

if __name__ == "__main__":
    main()

'''
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/mock.150.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/mock.250.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/mock.250.bbqc.sampled.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/nat.mock.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/low.biofuel.150.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/low.biofuel.250.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/low.biofuel.250.bbqc.sampled.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/150.wet.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/250.wet.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/250.wet.bbqc.sampled.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/hot.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/med.bog.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/mock2pct.150.bbqc.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/mock2pct.250.bbqc.sampled.fastq
/global/dna/shared/data/functests/assembly/Meta/Bench.20150730/bbqc_data/cami.bbqc.fastq
'''
