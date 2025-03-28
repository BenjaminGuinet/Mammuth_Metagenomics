#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import pysam
import re,os
from sys import argv


dna_damage_dict = {"CpG": [0] * 31, "CT" : [0] * 31 ,"other":[0] * 31}
cpg_sites = [0] * 31
real_path = os.path.realpath(argv[1])
directory_file=os.path.dirname(real_path)+"/"

directory_file=argv[2]+"/"

focal_bamfile = (pysam.AlignmentFile(argv[1], "rb"))
readcounter = 0
for read in focal_bamfile.fetch():
    readcounter += 1
    read_mapping_quality = read.mapping_quality
    if read_mapping_quality > 0:
        read_length = len(read.query_sequence)
        read_sequence = read.query_sequence
        edit_distance = read.get_tag("NM")
        MD_string = read.get_tag("MD")
        cigar = read.cigarstring
        if "N" in read_sequence or "N" in MD_string or "I" in cigar or "D" in cigar or "N" in cigar or "S" in cigar or "H" in cigar or 'P' in cigar:
            continue
        else:
            """
            recreate reference sequence
            """
            refseq = ""
            newread_seq = ""
            MD = MD_string.split()[0].rstrip('\n')
            MDlist = re.findall('(\d+|\D+)',MD)
            MDcounter = 0
            for base in MDlist:
                if base.isdigit():
                    base = int(base)
                    refseq += read_sequence[MDcounter:MDcounter+base]
                    newread_seq += read_sequence[MDcounter:MDcounter+base]
                    MDcounter += int(base)
                elif '^' in base:
                    continue
                elif base.isalpha():
                    refseq += base
                    newread_seq += read_sequence[MDcounter]
                    MDcounter += len(base)
            if "^" not in MD_string:
                readpos = -1
                for refseq_allel,read_allel in zip(refseq,newread_seq):
                    readpos += 1
                    if readpos < 31 and readpos < len(newread_seq)-1: # and len(newread_seq) >= 35:
                        if refseq_allel == "C" and refseq[readpos+1] == "G":
                            cpg_sites[readpos] += 1
                            if read_allel == "T":
                                dna_damage_dict["CpG"][readpos] += 1
                        if refseq_allel == read_allel:
                            continue
                        elif "N" in refseq_allel or "N" in read_allel:
                            continue
                        else:
                            mutation_combination = refseq_allel + read_allel
                            if mutation_combination in dna_damage_dict:
                                dna_damage_dict[mutation_combination][readpos] += 1
                            else:
                                dna_damage_dict["other"][readpos] += 1
title = re.sub(".bam","",os.path.basename(argv[1]))
fig, ax = plt.subplots(figsize=(10,10))
plt.rcParams['figure.figsize'] = [20, 10]
plt.rcParams["figure.autolayout"] = True
plt.rcParams['axes.linewidth'] = 1.25
plt.rcParams['font.size'] = 12
ax.set_xticks(list(range(0, 31, 2)))
ax.set_ylim(0,0.5)
ax.set_xlim(0 - 0.5, 30 + 0.5)
ax.set_xlabel("distance from read end (bp)",weight='bold',fontsize=14)
ax.set_ylabel("mismatch frequency",weight='bold',fontsize=14)
for key,value in dna_damage_dict.items():
    x_axis_damage = []
    y_axis_damage = []
    for i in range(len(value)):
        x_axis_damage += [i]
        if key == "CT":
            y_axis_damage += [(value[i]*4)/readcounter]
        elif key == "CpG":
            if cpg_sites[i] > 0:
                y_axis_damage += [value[i]/cpg_sites[i]]
        elif key == "other":
            y_axis_damage += [value[i]/(readcounter)]
    if key == "other":
        ax.plot(x_axis_damage,y_axis_damage,linewidth=2,linestyle='dotted',color='black',label="other")
    elif key == "CT":
        ax.plot(x_axis_damage,y_axis_damage,label="C to T",linewidth=2,color="grey")
    elif key == "CpG":
      try:
        #print(x_axis_damage)
        #print(y_axis_damage)
        ax.plot(x_axis_damage,y_axis_damage,label="CpG to TpG",linestyle='dashed',linewidth=2,color="red")
      except:
        continue 
#
ax.legend(loc="upper right",edgecolor="black")
ax.set_title(title, fontsize=16, weight='bold') 
plt.savefig(directory_file+os.path.basename(argv[1])[0:-4] + "_dna_damage_plot.pdf",dpi=300)
print("Plot saved to : ",directory_file+os.path.basename(argv[1])[0:-4] + "_dna_damage_plot.pdf" )
plt.close()
outputfile = open(directory_file+os.path.basename(argv[1])[0:-4] + "_dna_damage.txt","w")
outputfile.write("distance_from_read_end\tC_to_T\tCpG_to_TpG\tother\n")
for i in range(31):
    try:
     outputfile.write(str(i+1) + "\t" + str((dna_damage_dict["CT"][i]*4)/readcounter) + "\t" + str(dna_damage_dict["CpG"][i]/cpg_sites[i]) + "\t" + str(dna_damage_dict["other"][i]/readcounter) + "\n")
    except:
     outputfile.write(str(i+1) + "\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\n")
outputfile.close()

"""
except:
      print("BAM file is wrong")
      outputfile = open(os.path.basename(argv[1])[0:-4] + "_dna_damage.txt","w")
      outputfile.write("distance_from_read_end\tC_to_T\tCpG_to_TpG\tother\n")
      from reportlab.lib.pagesizes import letter
      from reportlab.pdfgen import canvas
      c = canvas.Canvas(os.path.basename(argv[1])[0:-4] + "_dna_damage_plot.pdf", pagesize=letter)
      c.save()
"""
