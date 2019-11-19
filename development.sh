####
# Computerome
####

rm -rf metafunk2
git clone https://github.com/anttonalberdi/metafunk2.git
mkdir metafunk2_test
mkdir metafunk2_test/quality_filtering

module load AdapterRemoval/2.1.3
python metafunk2/metafunk2/run_metafunk2.py -n sample1 -1 metafunk2_test/GH2_3b_1.fastq.gz -2 metafunk2_test/GH2_3b_2.fastq.gz -r refgenome -o /home/projects/ku-cbd/people/antalb/metafunk2_test -t 8

####
# SCI-PY tutorial
####

#Go to project path
cd /Users/jpl786/github/metafunk2/

#Build in developer mode
python setup.py develop



#https://python-packaging-tutorial.readthedocs.io/en/latest/setup_py.html
#https://www.youtube.com/watch?v=xiI1i525ljE
