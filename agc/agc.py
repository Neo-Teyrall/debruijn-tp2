#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "William Margerit"
__copyright__ = "Universite Paris"
__credits__ = ["William Margerit"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "William Margerit"
__email__ = "mailliw.marg@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

#==============================================================
# Main program
#==============================================================
def read_fasta(amplicon_file: str,minseqlen):
    if ".zip" in amplicon_file:
        filin = gzip.open(amplicon_file)
    else:
        filin = open(amplicon_file)

    seq = ""
    for line in filin:
        if ">" in line:
            if len(seq) >= minseqlen: 
                yield(seq)
            seq = ""
        else :
            seq += line[:-1]
    if len(seq) >= minseqlen: 
                yield(seq)

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    out ={}
    for i in read_fasta(amplicon_file,minseqlen):
        if i in out:
            out[i]+= 1
        else :
            out.setdefault(i,1)
    out2 =[]
    for i in out:
        if out[i] >= mincount:
            out2.append([out[i], i])
    out2.sort()
    out2 = out2[::-1]
    for i in range(len(out2)):
        print([out2[i][1],out2[i][0]])
        yield [out2[i][1],out2[i][0]]
        pass

def get_chunks(sequence, chunk_size):
    i = 0
    out = []
    while i < len(sequence)and i/chunk_size < 4:
        if i + chunk_size > len(sequence):
            raise(ValueError)
        out.append(sequence[i:i+chunk_size])
        i += chunk_size
    return out

def cut_kmer(sequence, kmer_size):
    for i in range(len(sequence)-kmer_size+1):
        yield(sequence[i:i+kmer_size])


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    for i in cut_kmer(sequence,kmer_size):
        if i not in kmer_dict :
            kmer_dict.setdefault(i,[id_seq])
        else : 
            kmer_dict[i].append(id_seq)
    return kmer_dict

def search_mates(kmer_dict, sequence, kmer_size):
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
    cpt = 0 
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            cpt +=1
    return (cpt/len(alignment_list[0])*100)

 

def detect_chimera(perc_identity_matrix):
    exrat= []
    diff1 = False
    diff2 = False
    for i in perc_identity_matrix:
        sdt = statistics.stdev(i)
        exrat.append(sdt)
        if i[0] > i[1] and sdt>5:
            diff1 = True
        if i[0] < i[1] and sdt>5:
            diff2 = True
    mean = statistics.mean(exrat)
    if mean > 5 and diff1 and diff2 :
        return True
    return False

def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))

def get_unique(ids):
    return {}.fromkeys(ids).keys()



def chimera_removal(amplicon_file, minseqlen, mincount,
                    chunk_size, kmer_size):
    l = list dereplication_fulllength(amplicon_file,minseqlen, mincount):
        
        
        pass

def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    amplicon_file, minseqlen, mincount = "../tests/test_sequences.fasta", 200, 3
    
    a= list(dereplication_fulllength("../tests/test_sequences.fasta", 200, 3))

    a = cut_kmer("ATGGTCGT",2)
if __name__ == '__main__':
    main()
