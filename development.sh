####
# Computerome
####

rm -rf metafunk2
git clone https://github.com/anttonalberdi/metafunk2.git
python metafunk2/metafunk2/run_metafunk2.py -n testsamp -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r /home/projects/ku-cbd/people/antalb/mMyoMyo_m19_AffsNnoesSC.p1.fa.gz -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8 --includesteps 3

cat metafunk2_test/testsamp.log

mkdir metafunk2_test
mkdir metafunk2_test/quality_filtering

module load AdapterRemoval/2.1.3
module load pigz/2.3.4
module load seqkit/0.7.1
module load jre/1.8.0
module load bbmap/36.49
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
