# metafunk2
Taxonomic and functional metagenomics pipeline

Metafunk2 is a pipeline optimised to use in Computerome.


metafunk2 pipeline
1) Quality filtering using AdapterRemoval
2) Duplicate read removal using seqkit rmdup
3) Mapping reads against reference genome(s) using bwa
4) Metagenomic assembly using metaSpades

metafunk2_merged pipeline (starts from individual assemblies generated with the metafunk2 pipeline)
1) Reassembly (merging) of metagenomic assemblies using Minimus2
2) Mapping reads from individual samples against reassembly using hisat2
3) Binning using Maxbin and Metabat
4) Bin refinement using Binning_refiner and DAS_Tool
5) Bin quality check using checkm
6) Bin annotation using CAT bin
7) Mapping reads from individual samples to bins using bwa
