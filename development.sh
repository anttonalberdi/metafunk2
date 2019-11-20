####
# Computerome
####

rm -rf metafunk2
git clone https://github.com/anttonalberdi/metafunk2.git
python metafunk2/metafunk2/run_metafunk2.py -n blank -1 metafunk2_test/Blank06092017_1.fastq.gz -2 metafunk2_test/Blank06092017_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 8 --includesteps 6

#3
python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 250 --includesteps 4 5 6

#4
python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 100 --includesteps 4

workdir="/home/projects/ku-cbd/people/antalb/metafunk2_test"
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/metafunk2_test.err -o ${workdir}/metafunk2_test.out -l nodes=1:ppn=8,mem=250gb,walltime=0:01:00:00 -N metafunk2_test -de python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 250 --includesteps 4 5 6


cat metafunk2_test/testsamp.log
cat metafunk2_test/testsamp.stats

li metafunk2_test/blank.binning/metabat

mkdir metafunk2_test
mkdir metafunk2_test/quality_filtering

module unload gcc/5.1.0
module load AdapterRemoval/2.1.3
module load pigz/2.3.4
module load seqkit/0.7.1
module load jre/1.8.0
module load bbmap/36.49
module load samtools/1.8
module load bwa/0.7.15
module load anaconda3/4.0.0
module load spades/3.13.1
module load perl/5.24.0
module load metabat/2.12.1
module load maxbin/2.2.7
python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r /home/projects/ku-cbd/people/antalb/mMyoMyo_m19_AffsNnoesSC.p1.fa.gz -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 --includesteps 3

Blank06092017*

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -l nodes=1:ppn=1,mem=8gb,walltime=1:00:00:00 -N remove_dogwolf -de 'rm -r /home/projects/ku-cbd/people/antalb/dog_wolf2'

## Metafunk2_merged
module load ncbi-blast/2.6.0+
module load cd-hit/4.8.1
module load MUMmer/3.23
module load kentUtils/350
module load amos/20121115
