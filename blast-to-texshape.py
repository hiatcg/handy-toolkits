#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright @ ZHANG Yang
# Last Update: 2020-04-27

# Things to do before running this script
# 1. Blast, 
# Blast Result File Name: PlasmidID_PlasmidName_PrimerID
# Down

# 1. Mapping the Sample ID to NCBI Reference Sequence ID

import os
import re
import sys
import pprint

from Bio import SeqIO
from Bio import Entrez
from Bio import Geo

Entrez.email = "zy0781@connect.hku.hk"

sample_label = dict()

def get_gene_decription():
    title = ""

    with open("./sequence.gb", "r") as f:
        for record in SeqIO.parse(f, format="genbank"):
            title += record.description + ", " + (record.id).replace("_", "\_") + ": "
            url = "https://www.ncbi.nlm.nih.gov/nuccore/" + record.id
            if record.features:
                for feature in record.features:
                    if feature.type == "CDS":
                        # start point + 1 to match BLAST Result
                        title += str(feature.location.start + 1) + ".." + str(feature.location.end)

    return ("\\href{" + url + "}{" + title + "}\n\n")

def map_sample_to_template():
    # Mapping the Sample ID to NCBI Reference Sequence ID
    for file in os.listdir("."):
        if file.startswith("PL"):
            print(file)
            plasmid_id = os.path.splitext(file)[0].split('_')[0].strip()
            break
    
    with open("./sequence.gb", "r") as f:
        for record in SeqIO.parse(f, format="genbank"):
            sample_label[plasmid_id] =  record.id

    pprint.pprint(sample_label)
    


def reverse_complement(seq):
    result = ""
    for nucleotide in seq[::-1]:
        if nucleotide == "A":
            result += "T"
        elif nucleotide == "T":
            result += "A"
        elif nucleotide == "C":
            result += "G"
        elif nucleotide == "G":
            result += "C"
        else:
            result += nucleotide
    return result


def transform_blast_to_texshade():
    latex_string = ""
    latex_string_main = ""
    plasmid_id = ""

    for file in os.listdir("."):
        if file.endswith(".blast"):
            # plasmid_id: SAMPLE ID
            # subject_name: Matched NCBI Reference Sequence ID
            print(os.path.splitext(file)[0].split('_'))
            plasmid_id, plasmid_name, primer_name = os.path.splitext(file)[0].split('_')[0:3]
            subject_name = sample_label[plasmid_id].strip()

            # Calculate Label Length, for the alignment of the result
            plasmid_id_length = len(plasmid_id)
            subject_name_length = len(subject_name)

            if plasmid_id_length < subject_name_length:
                plasmid_id += (subject_name_length - plasmid_id_length) * " " 
            else:
                subject_name += (plasmid_id_length - subject_name_length) * " "
            aln_result_prefix = (max(plasmid_id_length, subject_name_length) + 1) * " "

            result_aln = ""
            full_query_string = ""
            full_subject_string = ""
            full_result_string = ""
            
            with open(file, "r") as f:
                query_start_num = 0
                subject_start_num = 0
                query_line = f.readline()

                # Find Result Starting Number
                # query_line like:# Query  3     AAGAGAACCCCAGTGGAGACCAAGATGAGAAGACCCCTGTTCTAGGGGATGCGAAATCTG  62
                query_start_num = re.search("\s*(\d{1,})\D{1,}(\d{1,})", query_line).groups()[0]
                print(query_start_num)
                
                result_line = f.readline()
                subject_line = f.readline()

                # Figure out the orentation of the Matched Sequence
                m = re.search("\s*(\d{1,})\D{1,}(\d{1,})", subject_line)
                if int(m.groups()[0]) > int(m.groups()[1]):
                    reverse_flag = True
                    # Following line is not essential
                    # subject_start_num = m.groups()[1]
                else:
                    reverse_flag = False
                    subject_start_num = m.groups()[0]
                    
                while query_line:
                    m = re.search("\s*(\d{1,})\D{1,}(\d{1,})", subject_line)
                    if reverse_flag is True:
                        subject_start_num = int(m.groups()[1])

                    # Extract PURE Sequences for Sample and Template
                    # 1. Dicard Numbers in Sample Sequence
                    clean_query_line = re.sub("\d{1}", " ", query_line)
                    clean_subject_line = re.sub("\d{1}", " ", subject_line)
                    
                    # 2. Find the start position of the sequences of Sample and Template
                    query_string_start = re.search(r"\w*\s*", clean_query_line).end()
                    subject_string_start = re.search(r"\w*\s*", clean_subject_line).end()
                    string_start = min(query_string_start, subject_string_start)
                    query_string = clean_query_line[string_start:].rstrip()
                    subject_string = clean_subject_line[string_start:].rstrip()

                    # Make sure no trailing gap (blank) in blast result.
                    query_string_end = len(query_string)
                    subject_string_end = len(subject_string)
                    if query_string_end < subject_string_end:
                        query_string += (subject_string_end - query_string_end) * " "
                    else:
                        subject_string += (query_string_end - subject_string_end) * " "

                    result_string = result_line[string_start: (string_start + max(query_string_end, subject_string_end))]

                    full_query_string += query_string
                    full_subject_string += subject_string
                    for index,  c in enumerate(result_string):
                        if c == "|":
                            full_result_string += "*"
                        elif query_string[index] == "-" or subject_string[index] == "-":
                            full_result_string += " "
                        else:
                            full_result_string += "."

                    blank_line = f.readline()
                    query_line = f.readline()
                    result_line = f.readline()
                    subject_line = f.readline()

                if reverse_flag:
                    full_query_string = reverse_complement(full_query_string.upper()) # Some BLAST Results contain small case characters 
                    full_subject_string = reverse_complement(full_subject_string.upper()) # Some BLAST Results contain small case characters 

                result_aln += plasmid_id + " " + full_query_string + "\n"
                result_aln += subject_name + " " + full_subject_string + "\n"
                result_aln += aln_result_prefix + full_result_string + "\n"
                print(result_aln)

            with open(file[:-5] + "texshape", "w", encoding="utf-8") as f:
                f.write(result_aln)

            latex_string += r"\begin{blastResult}{./" + f"{plasmid_id.strip()}_{plasmid_name}/" + f"{file[:-5] + 'texshape'}" + "}" + \
                "\n\t" + r"\startnumber{1}{" + f"{query_start_num}" + "}" + \
                "\n\t" + r"\startnumber{2}{" + f"{subject_start_num}" + "}\n\t" + \
                r"\showcaption[top]{\refsample{" + f"{plasmid_id.strip()}" + "} " + \
                f"{plasmid_name} sequencing result using " + r"\refsample{" + \
                f"{primer_name}}}" + "}\n\end{blastResult}\n\n"
            
            latex_string_main += r"\begin{blastResult}{./" + f"{file[:-5] + 'texshape'}" + "}" + \
                "\n\t" + r"\startnumber{1}{" + f"{query_start_num}" + "}" + \
                "\n\t" + r"\startnumber{2}{" + f"{subject_start_num}" + "}\n\t" + \
                r"\showcaption[top]{\refsample{" + f"{plasmid_id.strip()}" + "} " + \
                f"{plasmid_name} sequencing result using " + r"\refsample{" + \
                f"{primer_name}}}" + "}\n\end{blastResult}\n\n"
                

    with open(f"./{plasmid_id.strip()}_sequencing_results.tex", "w", encoding="utf-8") as f:
        f.write(get_gene_decription())
        f.write(latex_string)
    
    with open(f"./{plasmid_id.strip()}_sequencing_results_main.tex", "w", encoding="utf-8") as f:
        f.write(r"""
% !TEX TS-program = xelatex
% !TEX encoding = UTF-8 Unicode

\documentclass{article}

\usepackage[margin=0.5in, top=1in, bottom=1in]{geometry}


\usepackage{texshade}

\newenvironment{blastResult}[1]
{
	\begin{texshade}{#1}
		\shadingcolors{greens}
		\shownumbering{leftright}
		\showruler{bottom}{2}
		\rulersteps{5}
		\setsize{residues}{footnotesize}
		\setsize{names}{footnotesize}
		\setsize{numbering}{footnotesize}
		\hideconsensus
		\noblockskip
}
{ \end{texshade}}

\usepackage{hyperref}
\usepackage{hyperxmp}
\hypersetup{
	unicode=true, pdftoolbar=true,
	pdfmenubar=true, pdffitwindow=true,
	pdfstartview={FitH}, pdfnewwindow=true,
	colorlinks=true, citecolor=green,
	filecolor=magenta, urlcolor=cyan,
	pdfcreator={ZHANG Yang, zy0781@connect.hku.hk},	
	pdfauthor={ZHANG Yang, zy0781@connect.hku.hk},
	pdfproducer={ZHANG Yang, zy0781@connect.hku.hk},
	pdfkeywords={Sequencing Results, Proteomics},
	pdfcopyright={ZHANG Yang, zy0781@connect.hku.hk}
}

\newcommand{\refsample}[1]{#1}

\begin{document}

        """)
        f.write(get_gene_decription())
        f.write(latex_string_main)
        f.write(r"\end{document}")




try:
    map_sample_to_template()
    transform_blast_to_texshade()
except Exception as ex:
    print(ex)
    raise(ex)