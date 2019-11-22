####
# Computerome
####

rm -rf metafunk2
git clone https://github.com/anttonalberdi/metafunk2.git
python metafunk2/metafunk2/run_metafunk2.py -n blank -1 metafunk2_test/Blank06092017_1.fastq.gz -2 metafunk2_test/Blank06092017_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 8 --includesteps 3,4

python metafunk2/metafunk2/run_metafunk2.py -n blank2 -1 metafunk2_test/Blank06092017_1.fastq.gz -2 metafunk2_test/Blank06092017_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 8 -i 1,2,3,4

li metafunk2_test

#MERGED
rm -rf metafunk2
git clone https://github.com/anttonalberdi/metafunk2.git
#python metafunk2/metafunk2/run_metafunk2_merged.py -n tests -p /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 8

workdir="/home/projects/ku-cbd/people/antalb/metafunk2_test"
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/merged.err -o ${workdir}/merged.out -l nodes=1:ppn=8,mem=50gb,walltime=0:01:00:00 -N merged -de python metafunk2/metafunk2/run_metafunk2_merged.py -n tests -p /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 -m 8

cat metafunk2_test/testsamp.log
cat metafunk2_test/testsamp.stats

li metafunk2_test/blank.binning/metabat

mkdir metafunk2_test
mkdir metafunk2_test/quality_filtering

## Metafunk2
module load AdapterRemoval/2.1.3
module load pigz/2.3.4
module load seqkit/0.7.1
module load jre/1.8.0
module load bbmap/36.49
module load samtools/1.9
module load bwa/0.7.15
module load anaconda3/2.2.0
module load spades/3.13.1
module load perl/5.20.2

module load metabat/2.12.1 #conflicts with gcc
module load maxbin/2.2.7

## Metafunk2_merged
module load perl/5.20.2
module load ncbi-blast/2.6.0+
module load cd-hit/4.8.1
module load MUMmer/3.23
module load kentUtils/350
module load amos/20121115

#####################
### Real sample test
#####################

cp fish_metagenomics/metafunk_june19/RawData/GH1_7b_*.fastq.gz metafunk2_test2
cp fish_metagenomics/metafunk_june19/RawData/GH3_6b_*.fastq.gz metafunk2_test2
cp fish_metagenomics/metafunk_june19/RawData/AI0_11b_*.fastq.gz metafunk2_test2


workdir="/home/projects/ku-cbd/people/antalb/metafunk2_test2"
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/GH1_7b.err -o ${workdir}/GH1_7b.out -l nodes=1:ppn=8,mem=50gb,walltime=1:00:00:00 -N GH1_7b -de python metafunk2/metafunk2/run_metafunk2.py -n GH1_7b -1 metafunk2_test2/GH1_7b_1.fastq.gz -2 metafunk2_test2/GH1_7b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test2 -t 8 -m 50 -i 1,2,3,4 -k
xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/GH3_6b.err -o ${workdir}/GH3_6b.out -l nodes=1:ppn=8,mem=50gb,walltime=1:00:00:00 -N GH3_6b -de python metafunk2/metafunk2/run_metafunk2.py -n GH3_6b -1 metafunk2_test2/GH3_6b_1.fastq.gz -2 metafunk2_test2/GH3_6b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test2 -t 8 -m 50 -i 1,2,3,4 -k

xqsub -V -A ku-cbd -W group_list=ku-cbd -d `pwd` -e ${workdir}/AI0_11b.err -o ${workdir}/AI0_11b.out -l nodes=1:ppn=8,mem=50gb,walltime=1:00:00:00 -N AI0_11b -de python metafunk2/metafunk2/run_metafunk2.py -n AI0_11b -1 metafunk2_test2/AI0_11b_1.fastq.gz -2 metafunk2_test2/AI0_11b_2.fastq.gz -r 'gambusia=databases/GCA_003097735.1_ASM309773v1_genomic.fna.gz,human=databases/GCF_000001405.39_GRCh38.p13_genomic.fna.gz' -o /home/projects/ku-cbd/people/antalb/metafunk2_test2 -t 8 -m 50 -i 1,2,3,4 -k
