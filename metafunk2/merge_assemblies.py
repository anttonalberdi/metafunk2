#!/usr/bin/env python
# -*- coding: utf-8

"""The script for merging assemblies"""

import os
import sys
import random
import argparse
import subprocess
import time
import gzip

def merge_assemblies(projectname,projectpath,threads,memory,logfilepath):
    #Create merged and merged/assembly directories if do not exist
    merged_dir = "merged"
    merged_abs = os.path.join(projectpath, merged_dir)
    if not os.path.exists(merged_abs):
        os.makedirs(merged_abs)
    assembly_dir = "assembly"
    assembly_abs = os.path.join(merged_abs, assembly_dir)
    if not os.path.exists(assembly_abs):
        os.makedirs(assembly_abs)

    assembliespath = os.path.join(projectpath,'*.assembly', 'contigs.fasta')
    mergedassembliespath = os.path.join(assembly_abs, 'assemblies.fna')
    mergedassembliesbase = os.path.join(assembly_abs, 'assemblies')
    nrassembliespath = os.path.join(assembly_abs, 'assemblies.nr.fna')
    afgassembliespath = os.path.join(assembly_abs, 'assemblies.afg')

    #Concatenate assemblies
    logfile=open(logfilepath,"a+")
    current_time = time.strftime("%m.%d.%y %H:%M", time.localtime())
    logfile.write("{0} |    Concatenating assemblies \r\n".format(current_time))
    logfile.close()
    concCmd = 'cat '+assembliespath+' > '+mergedassembliespath+''
    subprocess.check_call(concCmd, shell=True)

    #Removing redundant contigs
    #cdhitCmd = 'cd-hit -i '+mergedassembliespath+' -o '+nrassembliespath+' -T '+threads+' -M 0 -c 0.99 -d 100 -aS 0.9'
    #subprocess.check_call(cdhitCmd, shell=True)

    #Modify merged assembly to afg format
    toamosCmd = 'toAmos -s '+mergedassembliespath+' -o '+afgassembliespath+''
    subprocess.check_call(toamosCmd, shell=True)

    #Decomposed minimus2 pipeline (#https://github.com/nathanhaigh/amos/blob/master/src/Pipeline/minimus2.acf)
    mergedassemblies_bnk = os.path.join(mergedassembliesbase + '.bnk')
    mergedassemblies_afg = os.path.join(mergedassembliesbase + '.afg')
    mergedassemblies_refseq = os.path.join(mergedassembliesbase + '.ref.seq')
    mergedassemblies_qryseq = os.path.join(mergedassembliesbase + '.qry.seq')
    mergedassemblies_delta = os.path.join(mergedassembliesbase + '.delta')
    mergedassemblies_coords = os.path.join(mergedassembliesbase + '.coords')
    mergedassemblies_ovl = os.path.join(mergedassembliesbase + '.ovl')
    mergedassemblies_OVL = os.path.join(mergedassembliesbase + '.OVL')
    mergedassemblies_contig = os.path.join(mergedassembliesbase + '.contig')
    mergedassemblies_reassembly  = os.path.join(merged_abs, 'reassembly.fna')

    #Remove path if does not exist
    rmbankCmd = 'rm -fr '+mergedassemblies_bnk+''
    subprocess.check_call(rmbankCmd, shell=True)

    #Create bank
    bankCmd = 'bank-transact -c -z -b '+mergedassemblies_bnk+' -m '+mergedassemblies_afg+''
    subprocess.check_call(bankCmd, shell=True)

    #Dump1
    dump1Cmd = 'dumpreads '+mergedassemblies_bnk+' -M 0 > '+mergedassemblies_refseq+''
    subprocess.check_call(dump1Cmd, shell=True)

    #Dump2
    dump2Cmd = 'dumpreads '+mergedassemblies_bnk+' -m 0 > '+mergedassemblies_qryseq+''
    subprocess.check_call(dump2Cmd, shell=True)

    #Nucmer
    nucmerCmd = 'nucmer -maxmatch -c 100 '+mergedassemblies_refseq+' '+mergedassemblies_qryseq+' -p '+mergedassembliesbase+''
    subprocess.check_call(nucmerCmd, shell=True)

    #Coords
    coordsCmd = 'show-coords -H -c -l -o -r -I 95 '+mergedassemblies_delta+' | nucmerAnnotate | egrep "BEGIN|END|CONTAIN|IDENTITY" > '+mergedassemblies_coords+''
    subprocess.check_call(coordsCmd, shell=True)

    #ovl
    ovlCmd = 'nucmer2ovl -ignore 20 -tab '+mergedassemblies_coords+' | sort2 > '+mergedassemblies_ovl+''
    subprocess.check_call(ovlCmd, shell=True)

    #OVL
    OVLCmd = 'ovl2OVL '+mergedassemblies_ovl+' > '+mergedassemblies_OVL+''
    subprocess.check_call(OVLCmd, shell=True)

    #Transact
    transactCmd = 'bank-transact -z -b '+mergedassemblies_bnk+' -m '+mergedassemblies_OVL+''
    subprocess.check_call(transactCmd, shell=True)

    #Tigger
    tiggerCmd = 'tigger -b '+mergedassemblies_bnk+''
    subprocess.check_call(tiggerCmd, shell=True)

    #Consensus
    consensusCmd = 'make-consensus -B -e 0.06 -b '+mergedassemblies_bnk+' -w 15'
    subprocess.check_call(consensusCmd, shell=True)

    #Contig
    contigCmd = 'bank2contig '+mergedassemblies_bnk+' > '+mergedassemblies_contig+''
    subprocess.check_call(contigCmd, shell=True)

    #Fasta
    fastaCmd = 'bank2fasta -b '+mergedassemblies_bnk+' > '+mergedassemblies_reassembly+''
    subprocess.check_call(fastaCmd, shell=True)
