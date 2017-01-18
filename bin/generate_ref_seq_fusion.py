#!/usr/bin/python

import sys
import os
import re

if len(sys.argv) >= 6:
    gap_len = int(sys.argv[1])
    segment_len = int(sys.argv[2])
    ref_filename = sys.argv[3]
    fusion_gpd_filename = sys.argv[4]
    output_filename = sys.argv[5]
else:
    sys.exit(1)


########
# Fasta read
chr_seq = dict()
file = open(ref_filename, 'r')
p = re.compile('>([^\n]+)\n([^>]+)')
for m in p.finditer(file.read()):
  m1 = re.match('(\S+)',m.group(1))
  seq = m.group(2).replace("\n",'')
  seq = seq.replace("\r",'')
  seq = seq.replace(" ",'')
  chr_seq[m1.group(1)] = seq
file.close()

########
base_dict = dict()
base_dict['A'] = 'T'
base_dict['a'] = 't'
base_dict['C'] = 'G'
base_dict['c'] = 'g'
base_dict['G'] = 'C'
base_dict['g'] = 'c'
base_dict['T'] = 'A'
base_dict['t'] = 'a'

def reverse_complement(str_in, strand):
    
    if (strand == '+'):
        return str_in

    str_out = list(str_in) # To be more mem efficient
    len_str = len(str_in)
    for i in range(len_str):
        if (base_dict.has_key(str_in[len_str - i - 1])):
            str_out[i] = base_dict[str_in[len_str - i - 1]]
        else:
            str_out[i] = str_in[len_str - i - 1]
    return ''.join(str_out)

def reverse_seg(fields):
    num_blocks = int(fields[8])
    chr = fields[2]
    if (num_blocks == 1):
        end_ls = [int(fields[5])]
        start_pos = max(end_ls[0] - segment_len, int(fields[4]))
    else:
        # Include both start and end position (0-based)
        start_ls = [int(i) for i in fields[9].split(',')[:-1]]
        end_ls = [int(i) for i in fields[10].split(',')[:-1]]
        seg_len = 0
        for block_idx in range(num_blocks-1, -1, -1):
            block_len =  end_ls[block_idx] - start_ls[block_idx]
            if ((seg_len + block_len) <= segment_len):
                start_pos = start_ls[block_idx]
            else:
                start_pos = end_ls[block_idx] - (segment_len - seg_len)
        
            seg_len += block_len
            if (seg_len >= segment_len):
                break
        
    return [start_pos, end_ls[-1] + gap_len]
    

def forward_seg(fields):
    
    num_blocks = int(fields[8])
    chr = fields[2]
    if (num_blocks == 1):
        start_ls = [int(fields[4])]
        end_pos = min(start_ls[0] + segment_len, int(fields[5]))
    else:
        start_ls = [int(i) for i in fields[9].split(',')[:-1]]
        end_ls = [int(i) for i in fields[10].split(',')[:-1]]
        seg_len = 0
        for block_idx in range(num_blocks):
            block_len = end_ls[block_idx] - start_ls[block_idx] 
            if ((seg_len + block_len) <= segment_len):
                end_pos = end_ls[block_idx]
            else:
                end_pos = start_ls[block_idx] + (segment_len - seg_len)
        
            seg_len += block_len
            if (seg_len >= segment_len):
                break
            
    return [start_ls[0] - gap_len, end_pos]


# Note: we need to check both forward and reverse order for compatibility
def add_seq_info(ref_seq_dict, common_item_ls, readname, seq_1_info, seq_2_info):
    chr_1 = seq_1_info[-1]
    chr_2 = seq_2_info[-1]

    update = False

    # Format of storing each ref information [start_1 end_1 dir_1 start_2 end_2 dir_2 [readnames] [direction]]
    updated_item_1 = -1  # This is used for merging if current read overlap with two or more disjoint references
    updated_item_2 = -1  
    if (ref_seq_dict.has_key(chr_1)):
        if (ref_seq_dict[chr_1].has_key(chr_2)):
            item_idx_1 = 0
            while (item_idx_1 < len(ref_seq_dict[chr_1][chr_2])):
                item = ref_seq_dict[chr_1][chr_2][item_idx_1]
                # Check if there is any overlap, then extend the ref sequence
                if ((item[2] == seq_1_info[2]) and
                    (min(item[1], seq_1_info[1]) - max(item[0], seq_1_info[0])) > 0):
                    if ((item[5] == seq_2_info[2]) and
                        (min(item[4], seq_2_info[1]) - max(item[3], seq_2_info[0])) > 0):
                        if (updated_item_1 == -1):
                            ref_seq_dict[chr_1][chr_2][item_idx_1] = [min(item[0], seq_1_info[0]), max(item[1], seq_1_info[1]), item[2],
                                                                      min(item[3], seq_2_info[0]), max(item[4], seq_2_info[1]), item[5],
                                                                      item[6] + [readname], item[7] + ['+']]
                            updated_item_1 = item_idx_1
                        else:
                            ref_seq_dict[chr_1][chr_2][updated_item_1] = [min(item[0], ref_seq_dict[chr_1][chr_2][updated_item_1][0]), max(item[1], ref_seq_dict[chr_1][chr_2][updated_item_1][1]), item[2],
                                                                          min(item[3], ref_seq_dict[chr_1][chr_2][updated_item_1][3]), max(item[4], ref_seq_dict[chr_1][chr_2][updated_item_1][4]), item[5],
                                                                          item[6] + ref_seq_dict[chr_1][chr_2][updated_item_1][6], item[7] + ref_seq_dict[chr_1][chr_2][updated_item_1][7]]
                            print "Warning: following condition is not tested throughly on real data"
                            print "Two reference region merged on " + readname
                            print "\t" + str(ref_seq_dict[chr_1][chr_2][updated_item_1]) + "\t" + str(ref_seq_dict[chr_1][chr_2][item_idx_1])
                            del ref_seq_dict[chr_1][chr_2][item_idx_1]
                            item_idx_1 -= 1  # It will be post incremented
                        update = True
                        
                item_idx_1 += 1
        else:
            ref_seq_dict[chr_1][chr_2] = list()
    else:
        ref_seq_dict[chr_1] = dict()
        ref_seq_dict[chr_1][chr_2] = list()
    
    if (ref_seq_dict.has_key(chr_2)):
        if (ref_seq_dict[chr_2].has_key(chr_1)):
            item_idx_2 = 0
            while (item_idx_2 < len(ref_seq_dict[chr_2][chr_1])):
                item = ref_seq_dict[chr_2][chr_1][item_idx_2]
                # Check if there is any overlap, then extend the ref sequence
                if ((item[2] != seq_2_info[2]) and
                    (min(item[1], seq_2_info[1]) - max(item[0], seq_2_info[0])) > 0):
                    if ((item[5] != seq_1_info[2]) and
                        (min(item[4], seq_1_info[1]) - max(item[3], seq_1_info[0])) > 0):
                        if ((updated_item_1 == -1) and
                            (updated_item_2 == -1)):
                            ref_seq_dict[chr_2][chr_1][item_idx_2] = [min(item[0], seq_2_info[0]), max(item[1], seq_2_info[1]), item[2],
                                                                      min(item[3], seq_1_info[0]), max(item[4], seq_1_info[1]), item[5],
                                                                      item[6] + [readname], item[7] + ['-']]
                            updated_item_2 = item_idx_2
                        elif (updated_item_1 == -1):
                            ref_seq_dict[chr_2][chr_1][updated_item_2] = [min(item[0], ref_seq_dict[chr_2][chr_1][updated_item_2][0]), max(item[1], ref_seq_dict[chr_2][chr_1][updated_item_2][1]), item[2],
                                                                          min(item[3], ref_seq_dict[chr_2][chr_1][updated_item_2][3]), max(item[4], ref_seq_dict[chr_2][chr_1][updated_item_2][4]), item[5],
                                                                          item[6] + ref_seq_dict[chr_2][chr_1][updated_item_2][6], item[7] + ref_seq_dict[chr_2][chr_1][updated_item_2][7]]
                            print "Warning: following condition is not tested throughly on real data"
                            print "Two reference region merged on " + readname
                            print "\t" + str(ref_seq_dict[chr_2][chr_1][updated_item_2]) + "\t" + str(ref_seq_dict[chr_2][chr_1][item_idx_2])                                                        
                            del ref_seq_dict[chr_2][chr_1][item_idx_2]
                            item_idx_2 -= 1  # It will be post incremented
                        else:
                            # As this is the reverse orientation need to change dir "+" -> "-" and "-" to "+"
                            dir_item = []
                            for dir_elem in item[7]:
                                if dir_elem == "+":
                                    dir_item.append("-")
                                else:
                                    dir_item.append("+")
            
                            ref_seq_dict[chr_1][chr_2][updated_item_1] = [min(item[3], ref_seq_dict[chr_1][chr_2][updated_item_1][0]), max(item[4], ref_seq_dict[chr_1][chr_2][updated_item_1][1]), ref_seq_dict[chr_1][chr_2][updated_item_1][2],
                                                                          min(item[0], ref_seq_dict[chr_1][chr_2][updated_item_1][3]), max(item[1], ref_seq_dict[chr_1][chr_2][updated_item_1][4]), ref_seq_dict[chr_1][chr_2][updated_item_1][5],
                                                                          item[6] + ref_seq_dict[chr_1][chr_2][updated_item_1][6], dir_item + ref_seq_dict[chr_1][chr_2][updated_item_1][7]]
                            print "Warning: following condition is not tested throughly on real data"
                            print "Two reference region merged on " + readname
                            print "\t" + str(ref_seq_dict[chr_1][chr_2][updated_item_1]) + "\t" + str(ref_seq_dict[chr_2][chr_1][item_idx_2])                            
                            del ref_seq_dict[chr_2][chr_1][item_idx_2]                        

                            if ((chr_1 == chr_2) and
                                (item_idx_2 <  updated_item_1)):
                                updated_item_1 -= 1
                            item_idx_2 -= 1  # It will be post incremented
                        update = True
                item_idx_2 += 1
        else:
            ref_seq_dict[chr_2][chr_1] = list()
    else:
        ref_seq_dict[chr_2] = dict()
        ref_seq_dict[chr_2][chr_1] = list()

    if (not update):
        ref_list = seq_1_info[0:3] + seq_2_info[0:3] + [[readname], ['+']]

        ref_seq_dict[chr_1][chr_2].append(ref_list)


"""
readname_compatible_ls = []
gpd_compatible_file = open(fusion_compatible_gpd_filename, 'r')
for line in gpd_compatible_file:
    readname = re.split(r'\+|\-', (line.split()[0].split('|')[1]))[0]
    readname_compatible_ls.append(readname)
gpd_compatible_file.close()
"""

ref_seq_dict = dict()
readname_info = dict()   # Keep reads information to be save in info file
common_item_ls = list()
gpd_file = open(fusion_gpd_filename, 'r')
output_file = open(output_filename, 'w')
output_info_file = open(output_filename + '.info', 'w') # Keep the info about ref sequences to be used in MLE_input file
while (True):
    line1 = gpd_file.readline()
    if (not line1):
        break
    line2 = gpd_file.readline()

    fields_1 = line1.split('\t')
    readname = re.split(r'\+|\-', (fields_1[0].split('|')[1]))[0]
    #if (readname not in readname_compatible_ls):
    #   continue

    chr = fields_1[2]
    if (not chr_seq.has_key(chr)):
        continue
       
    fields_2 = line2.split('\t')
    chr = fields_2[2]
    if (not chr_seq.has_key(chr)):
        continue
    "Note: Sequences include both exon and intron segments"
    readname_info[readname] = [fields_1[2]]
    if (fields_1[3] == '+'):
        seq_1_info = reverse_seg(fields_1) + ['+', fields_1[2]]
        readname_info[readname] += ['+', seq_1_info[1] - gap_len]
    else:
        seq_1_info = forward_seg(fields_1) + ['-', fields_1[2]]
        readname_info[readname] += ['-', seq_1_info[0] + gap_len]

    readname_info[readname] += [fields_2[2]]
    if (fields_2[3] == '-'):
        seq_2_info = reverse_seg(fields_2) + ['-', fields_2[2]]
        readname_info[readname] += ['-', seq_2_info[1] - gap_len]
    else:
        seq_2_info = forward_seg(fields_2) + ['+', fields_2[2]]
        readname_info[readname] += ['+', seq_2_info[0] + gap_len]

    add_seq_info(ref_seq_dict, common_item_ls, readname, seq_1_info, seq_2_info)
# TODO: Join overlapping items
if (len(common_item_ls) > 0):
    print "Warning: overlapping reference seqeuces: " + str(common_item_ls)

for chr1 in ref_seq_dict.keys():
    for chr2 in ref_seq_dict[chr1].keys():
        for item in ref_seq_dict[chr1][chr2]:
            refname = ('_'.join([chr1, str(item[0]), str(item[1]), item[2]]) + '/' +
                       '_'.join([chr2, str(item[3]), str(item[4]), item[5]]))
            output_file.write('>' + refname + '\n')
            output_file.write(reverse_complement(chr_seq[chr1][item[0]:item[1]], item[2]) +
                              reverse_complement(chr_seq[chr2][item[3]:item[4]], item[5]) + '\n')
            readname_idx = -1 # Pre-increment
            for readname in item[6]:
                readname_idx += 1
                output_info_file.write(readname + '\t' + refname)
                seg_1_len = item[1] - item[0]
                seg_2_len = item[4] - item[3]
                # Check if it is in the same order
                ref_order = ''
                if (item[7][readname_idx] == '+'):
                       seg_1_junction_point = -item[0] + readname_info[readname][2]
                       seg_2_junction_point = -item[3] + readname_info[readname][5]

                else:
                    seg_1_junction_point = -item[0] + readname_info[readname][5]
                    seg_2_junction_point = -item[3] + readname_info[readname][2]

                # To find junction location in ref, if seq is reversed complement we need to adjust properly
                if (item[2] == '-'):
                    seg_1_junction_point = seg_1_len - seg_1_junction_point
                if (item[5] == '-'):
                    seg_2_junction_point = seg_2_len - seg_2_junction_point
                seg_2_junction_point += seg_1_len
                output_info_file.write('\t' + '\t'.join([str(seg_1_junction_point), str(seg_2_junction_point), item[7][readname_idx]]) +
                                       '\n')

output_file.close()
gpd_file.close()
output_info_file.close()
