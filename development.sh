####
# Computerome
####

rm -rf metafunk2
git clone https://github.com/anttonalberdi/metafunk2.git
python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 --includesteps 3

workdir="/home/projects/ku-cbd/people/antalb/metafun2_test"
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/metafun2_test.err -o ${workdir}/metafun2_test.out -l nodes=1:ppn=8,mem=8gb,walltime=1:00:00:00 -N metafun2_test -de python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 --includesteps 3

cat metafunk2_test/testsamp.log
li metafunk2_test/genomes

mkdir metafunk2_test
mkdir metafunk2_test/quality_filtering

module load AdapterRemoval/2.1.3
module load pigz/2.3.4
module load seqkit/0.7.1
module load jre/1.8.0
module load bbmap/36.49
module load samtools/1.8
module load bwa/0.7.15
python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r /home/projects/ku-cbd/people/antalb/mMyoMyo_m19_AffsNnoesSC.p1.fa.gz -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 --includesteps 3




####
# SCI-PY tutorial
####

#Go to project path
cd /Users/jpl786/github/metafunk2/

#Build in developer mode
python setup.py build
python setup.py install



#https://python-packaging-tutorial.readthedocs.io/en/latest/setup_py.html
#https://www.youtube.com/watch?v=xiI1i525ljE
