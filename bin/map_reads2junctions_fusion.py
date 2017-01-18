#!/usr/bin/python

import sys
import os
import commands
import re

if len(sys.argv) >= 8:
    run_aligner = int(sys.argv[1])
    fusion_ref_filename = sys.argv[2]
    unmapped_reads_filename = sys.argv[3]
    num_threads = sys.argv[4]
    fusion_candidates_filename = sys.argv[5]
    fusion_output_filename = sys.argv[6]
    margin = int(sys.argv[7])
    min_junction_overlap_len = int(sys.argv[8])
    temp_foldername = sys.argv[9]
    python_path =  sys.argv[10]
    bin_foldername = sys.argv[11]
    splice_mapper_path = sys.argv[12]
    # - indicates its installed, so if it is we will just use STAR
    if splice_mapper_path == '-': splice_mapper_path = ""
    else: splice_mapper_path += '/'
    splice_mapper_type = sys.argv[13]
    splice_mapper_options = sys.argv[14:]
    
else:
    sys.exit(1)
################################################################################

# Each ref sequence should be in seperate file (mapsplice_restriction)
if (splice_mapper_type == "STAR"):
    splice_mapper_ref_foldername = temp_foldername + 'star_ref/'
    splice_mapper_out_foldername = temp_foldername + "star_out/"
    junctions_filename = splice_mapper_out_foldername + 'SJ.out.tab_update'
    num_uniq_reads_idx = 6
    num_multi_reads_idx = 7
    num_nonredun_reads_idx = 9
    start_pos_offset = -1
    end_pos_offset = 0
elif (splice_mapper_type == "MapSplice"):
    splice_mapper_ref_foldername = temp_foldername + 'mapsplice_ref/'
    splice_mapper_out_foldername = temp_foldername + "mapsplice_out/"
    junctions_filename = splice_mapper_out_foldername + 'junctions.txt'
    num_uniq_reads_idx = 20
    num_multi_reads_idx = 21
    num_nonredun_reads_idx = 29
    start_pos_offset = 0
    end_pos_offset = -1
os.system("mkdir " + splice_mapper_ref_foldername)
os.system("mkdir " + splice_mapper_out_foldername)

num_lines = int(commands.getstatusoutput('wc -l ' + fusion_ref_filename)[1].split()[0])
# To make sure they are sorted in junctions.txt (otherwise ref_10 comes before ref_2)
ref_str_len = len(str(num_lines / 2))

junction_ref_file = open(fusion_ref_filename, 'r')
i = 1

# Lets not store these references in separate files.
file = open(splice_mapper_ref_foldername + 'ref_1.fa', 'w')
while (True):
    line = junction_ref_file.readline()
    if not line:
        break
    if (line[0] != '>'):
        print("Err: Invalid fasta format")
        exit(1)
    line = junction_ref_file.readline()
    if (line[0] == '>'):
        print("Err: Invalid fasta format")
        exit(1)
    filename = 'ref_' + ''.join(['0'] * (ref_str_len - len(str(i)))) + str(i)
    file.write('>' + filename + '\n')
    file.write(line.rstrip()+"\n")
    i += 1
file.close()

junction_ref_file.close()

if (run_aligner):
    if (splice_mapper_type == "STAR"):
        os.system("mkdir " + splice_mapper_out_foldername + "/genome")
        star_cmnd = (splice_mapper_path + "STAR --genomeChrBinNbits 16 --runMode genomeGenerate --genomeDir " + splice_mapper_out_foldername + 
                     "/genome --genomeFastaFiles " + splice_mapper_ref_foldername + "/*   --runThreadN " + str(num_threads))
        print(star_cmnd)
        os.system(star_cmnd)

        begin_dir = os.getcwd()
        os.chdir(splice_mapper_out_foldername)
        star_cmnd = (splice_mapper_path + "STAR --genomeDir " + splice_mapper_out_foldername + 
                     "/genome --readFilesIn " + unmapped_reads_filename + " --runThreadN " + str(num_threads) + 
                      " --outSJfilterCountUniqueMin 0 0 0 0 --outSJfilterCountTotalMin 1 1 1 1 " + 
                      " --outSJfilterOverhangMin " + ' '.join([str(min_junction_overlap_len)] * 4) + " " + 
                      " ".join(splice_mapper_options) +
                      " --outFilterMultimapNmax " + str(num_lines / 2))
        
        print(star_cmnd)
        os.system(star_cmnd)
        #print "Warning: STAR is commented out "
        junction_update_cmnd = (python_path + " " + bin_foldername + "add_nonredundant_to_star.py Aligned.out.sam SJ.out.tab " +
                str(min_junction_overlap_len) + " | sort -k1 > SJ.out.tab_update")
        print(junction_update_cmnd)
        os.system(junction_update_cmnd)
        os.chdir(begin_dir)
        
        
    elif (splice_mapper_type == "MapSplice"):
        # Run mapsplice
        # Build bowtie index, for large number of reference files, bowtie-build seems to have problem with building indices
        os.system("cat " + splice_mapper_ref_foldername + "/* > " + splice_mapper_out_foldername + "/ref.fa")
        os.system("bowtie-build " + splice_mapper_out_foldername + "/ref.fa " + splice_mapper_out_foldername + "/built_bowtie_index ")

        mapsplice_cmnd = ("python " + splice_mapper_path + "/mapsplice.py -1 " + unmapped_reads_filename +
                          " -x " + splice_mapper_out_foldername + "/built_bowtie_index " + " -c " + splice_mapper_ref_foldername + " -p " + num_threads + 
                          " -o " + splice_mapper_out_foldername  + " " + " ".join(splice_mapper_options))
        print(">>>Running Mapsplice:")
        print(mapsplice_cmnd)
        os.system(mapsplice_cmnd)
        #print "Warning: mapsplice is commented out "
else:
    print(splice_mapper_type + " command is skipped.")

# Generate a dictionary, mapping reads to ref sequences
ref_dict = dict()
fusion_ref_file_info = open(fusion_ref_filename + '.info', 'r')
for line in fusion_ref_file_info:
    fields = line.split()
    if (ref_dict.has_key(fields[1])):
        ref_dict[fields[1]].append([fields[0], int(fields[2]), int(fields[3]), fields[4]])
    else:
        ref_dict[fields[1]] = [[fields[0], int(fields[2]), int(fields[3]), fields[4]]]
    

# Parse junctions file
fusion_fa_file = open(fusion_ref_filename, 'r')
junction_out_file = open(temp_foldername + 'junctions.txt.notused', 'w')

ref_num = 1
fusion_junction = ''
ref_name = ''
fusion_score = 0
junctions_file = open(junctions_filename, 'r')  # Note:  1-based pos format and intron neighbour positions 
fusion_candidate_dict = dict()
skip_flag = True   # To pass first "fields[0] != ref_name" condition
for line in junctions_file:
    fields = line.split()
    fields[1] = int(fields[1]) + start_pos_offset
    fields[2] = int(fields[2]) + end_pos_offset
    while (fields[0] != ref_name):
        if (not skip_flag):
            for read_idx in range(num_reads):
                if (fusion_junction[read_idx] != ''):
                    junction_out_file.write('\t'.join([ref_dict[ref_line][read_idx][0], ref_name, ref_line,
                                                       str(fusion_junction[read_idx][0]), str(fusion_junction[read_idx][1]),
                                                       str(num_uniq_reads_ls[read_idx]), str(num_multi_reads_ls[read_idx]),
                                                       str(num_nonredun_reads_ls[read_idx])]) + '\n')
                    fusion_candidate_dict[ref_dict[ref_line][read_idx][0]] = fusion_junction[read_idx] + [ref_dict[ref_line][read_idx][3]]
         
        skip_flag = False
        ref_name = 'ref_' + ''.join(['0'] * (ref_str_len - len(str(ref_num)))) + str(ref_num)
        ref_line = fusion_fa_file.readline().strip()[1:]
        if not ref_line:
            break
        fusion_fa_file.readline()  # skip the seq line
        #line_f1 = fusion_psl_temp_file.readline()
        #line_f2 = fusion_psl_temp_file.readline()
        gene_pos = re.split(r'/|_', ref_line)
        num_reads = len(ref_dict[ref_line])
        fusion_junction = [''] * num_reads
        fusion_score = [0] * num_reads
        num_uniq_reads_ls = [0] * num_reads
        num_multi_reads_ls = [0] * num_reads
        num_nonredun_reads_ls = [0] * num_reads
        
        ref_num += 1
    
    num_uniq_reads = int(fields[num_uniq_reads_idx])
    for read_idx in range(num_reads):
        if ((abs(fields[1] - ref_dict[ref_line][read_idx][1]) < margin) and
            (abs(fields[2] - ref_dict[ref_line][read_idx][2]) < margin)):
            if (num_uniq_reads > fusion_score[read_idx]):
                fusion_score[read_idx] = num_uniq_reads
                num_uniq_reads_ls[read_idx] = int(fields[num_uniq_reads_idx])
                num_multi_reads_ls[read_idx] = int(fields[num_multi_reads_idx])
                num_nonredun_reads_ls[read_idx] = int(fields[num_nonredun_reads_idx])
                
                junc_pos1 = int(gene_pos[1]) 
                if (gene_pos[3] == '+'):
                    junc_pos1 += fields[1]  # includes start of gap
                else:
                    junc_pos1 += (int(gene_pos[2]) - int(gene_pos[1])) - fields[1]  # does not includes end of gap
                junc_pos2 = int(gene_pos[5]) # doesnt include end of gap
                if (gene_pos[7] == '+'):
                    junc_pos2 += fields[2] - (int(gene_pos[2]) - int(gene_pos[1])) # doesnt include end of gap
                else:
                    junc_pos2 += (int(gene_pos[6]) - int(gene_pos[5])) - (fields[2] - (int(gene_pos[2]) - int(gene_pos[1])))
                fusion_junction[read_idx] = [junc_pos1, junc_pos2]
if (not skip_flag):
    for read_idx in range(num_reads):
        if (fusion_junction[read_idx] != ''):
            junction_out_file.write('\t'.join([ref_dict[ref_line][read_idx][0], ref_name, ref_line,
                                                       str(fusion_junction[read_idx][0]), str(fusion_junction[read_idx][1]),
                                                       str(num_uniq_reads_ls[read_idx]), str(num_multi_reads_ls[read_idx]),
                                                       str(num_nonredun_reads_ls[read_idx])]) + '\n')
            fusion_candidate_dict[ref_dict[ref_line][read_idx][0]] = fusion_junction[read_idx] + [ref_dict[ref_line][read_idx][3]]
                                                       
fusion_fa_file.close()
junction_out_file.close()

fusion_candidates_file = open(fusion_candidates_filename, 'r')
fusion_gpd_file = open(fusion_output_filename, 'w')
fusion_candidate_ls = fusion_candidate_dict.keys()
for line in fusion_candidates_file:
    readname = re.split(r'\+|\-', (line.split()[0].split('|')[1]))[0]
    if (readname in fusion_candidate_ls):
        # Update junction pos in readname
        fields = line.strip().split()
        update_readname = fields[0].split('|')
        update_readname[1] = re.split(r'(\+|\-)',update_readname[1])
        if (fusion_candidate_dict[readname][2] == '+'):
            update_readname[1][-1] = str(fusion_candidate_dict[readname][0])
            fusion_candidate_dict[readname][2] = '-'
        else:
            update_readname[1][-1] = str(fusion_candidate_dict[readname][1])
            fusion_candidate_dict[readname][2] = '+'
        update_readname[1] = ''.join(update_readname[1])
        fusion_gpd_file.write('|'.join(update_readname) + '\t' + '\t'.join(fields[1:]) + '\n')
    
fusion_candidates_file.close()
fusion_gpd_file.close()


        
