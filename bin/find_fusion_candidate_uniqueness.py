#!/usr/bin/env python
import sys, re, os
import psl_basics, genepred_basics

#  Pre: <psl file> - of the format newname4_LR.bestpsl_fusion_pair
#       <bedGraph of uniqueness> - only needs the "1" entries
#       <junction.txt-like output>
#       <LR_fusion_candidates.gpd output>
#       <L_min_intron>
#  post: Writes out a "junctions.txt" type file to the filename you select, 
#        a tab delimited file with entries as follows
#    <fusion id> - derived from a piece of long read evidence
#    <artifical fusion name> - used for split read mapping short reads
#    <fasta header> - another name used for the artifical fusion name
#                     of the format chromosome, 0-indexed base, 1-indexed base, direction of sequence in artifical fusion, forward slash then the same
#                     This is mostly meaningless here since it is not an artifical sequence that gets generated and is mearly in the same format
#    <left junction> - if junction is on the left side of the read the value indicates
#                        the start, or leftmost position of the read before the breakpoit
#                        and is 0-indexed
#                      if junction is on the right side of the read the value indicated is
#                        the end, or rightmost base coordinate of the read before the breakpoint
#                        and is 1-indexed
#    <right junction> - if junction is on the left side of the read the value indicates
#                        the start, or leftmost position of the read before the breakpoit
#                        and is 0-indexed
#                      if junction is on the right side of the read the value indicated is
#                        the end, or rightmost base coordinate of the read before the breakpoint
#                        and is 1-indexed
#   <Type 1 nonredundant> - Not used
#   <Type 1 unique> - Not used
#   <Type 2> - Not used
#   <Type 3> - Not used
#   <Type 4> - Not used
#   <long read direction> - direction the long read evidence was oriented with respect to the artifical fusion
#                           this may be useful if putting together some polished reporting
#                           its the last field in the LR_fusion_candidates.fa.info file.
#   <junction count> - Number of different junctions within our long read fusion point margin of error
#   <junction string> - a string that starts with the index-1 base thats the last base before and the last base after the
#                       the long read fusion breakpoint.  Then it is followed by a "|" then 
#                       junctions (short read) and their read support counts for each junction 
#                       within our margine or long read fusion point error.  
#                       These coordiantes are 1-indexed base positions for the coordiante just before the breakpoint on either side.
#   <max of min length>
#   <max of min unique count>
#   <max of min flanking uniqueness> - Not used?

# Modifies:  Populates the  temp folder above with files necessary for 
#        counting read contributions fusion events
#        2.) temp/uniqueness/artificial_fusion_names.txt
#            This file is for converting the index of an artifical fusion to its name used in its sam file
#            <artifical fusion index> <artifical fusion name>
#        3.) temp/uniqueness/long_read_locations.txt
#            This file has the junction evidence observed in the long reads 
#            <artifical fusion index> <fusion id> <fasta header> <left junction> <right junction>

gap_or_overlap_search_multiplier = 2.5
# search for long read and short read where distance between the two
# is less than min search distance
# or is less than the larger of either gap or overlap value times 
#    the above variable, but is less than the max search distance

def main():
  if len(sys.argv) != 6:
    print sys.argv[0] + ' <psl file of fusion candidates> <uniqueness bedGraph> <junction output file> <fusion genepred out file> <L_min_intron>'
    sys.exit()

  psl_filename = sys.argv[1]
  uniqueness_bedgraph_filename = sys.argv[2]
  output_file = sys.argv[3]
  output_genepred_filename = sys.argv[4]
  min_gap_in_block = int(sys.argv[5]) # this value is used for smoothing indels when doing the psl to genepred conversion so there not so many blocks
                       # L_min_intron is used for this

  #1. For fusion ID's we have evidence for, convert the psl file into a genepred file
  #    we need our accurate junction point information so we can make observations about the base unqiueness
  print "make the fusion genepred file"
  fusion_gpds = make_fusion_genepred(psl_filename,min_gap_in_block,output_genepred_filename)

  #2. For fusion ID's we have evidence for get uniqueness information around the short read junction point
  print "set all fusion coordiantes on any fusion long read to zero"
  fusion_coordinates = initialize_fusion_coordinates(fusion_gpds)
  print "get the uniqueness value for all fusion coordinates"
  # modifies fusion_coordinates
  get_coordinate_uniqueness(uniqueness_bedgraph_filename, fusion_coordinates)
  print "calculate the uniqueness of each fusion"
  fid_least_unique = get_fusion_uniqueness(fusion_gpds, fusion_coordinates)  

  #3. Now we can enrich our per-long read file with information about the
  #    shortest long read support around around the junction
  print "add annotations to junctions.txt file"
  add_annotations_per_fusion(fusion_gpds,fid_least_unique, output_file)


def add_annotations_per_fusion(fusion_gpds, fid_least_unique, output_file):
  ofile = open(output_file,'w')
  for fid in fusion_gpds:
    ofile.write(fid + "\t" + '.' + "\t")
    #print fusion_gpds[fid]
    left = fusion_gpds[fid][0][0]
    right = fusion_gpds[fid][1][0]
    ml = re.search('\|F[\d]+([+-])',left['gene_name'])
    mr = re.search('\|F[\d]+([+-])',right['gene_name'])
    ldir = ml.group(1)
    rdir = mr.group(1)

    lname = left['chrom'] + '_' + str(left['txStart']) + '_' + str(left['txEnd'])+'_'+left['strand']
    rname = right['chrom'] + '_' + str(right['txStart']) + '_' + str(right['txEnd'])+'_'+right['strand']
    name = lname + '/' + rname
    ofile.write(name + "\t")

    lcoord = left['txStart']
    ltrue = lcoord + 1
    if ldir == '-':
      lcoord = left['txEnd']
      ltrue = lcoord
    rcoord = right['txStart']
    rtrue = rcoord +1
    if rdir == '-':
      rcoord = right['txEnd']
      rtrue = rcoord + 1
    ofile.write(str(lcoord)+"\t")
    ofile.write(str(rcoord)+"\t")
    ofile.write('.' + "\t")
    ofile.write('.' + "\t")
    ofile.write('.' + "\t")
    ofile.write('.' + "\t")
    ofile.write('.' + "\t")
    ofile.write('+' + "\t")
    ofile.write('.' + "\t")
    entry = left['chrom']+":"+str(ltrue) + '/' + right['chrom']+":"+str(rtrue)
    ofile.write(entry+"\t")
    osets = []
    ofile.write(str(fid_least_unique[fid]['shortest_total_length'])+"\t")
    ofile.write(str(fid_least_unique[fid]['shortest_unique_count'])+"\t")
    proxdic = fid_least_unique[fid]['lowest_proximal_unique_counts']
    proxdists = map(int,proxdic.keys())
    proxdists.sort()
    for dist in proxdists:
      osets.append(str(dist)+':'+str(proxdic[dist]))
    ostring = ",".join(osets)
    ofile.write(ostring+"\n")
  ofile.close()

# pre: genepred entries dictionary - keyed by fid with a list of two entries per fid
#      fusion_coordinates - for every fusion coordinate it holds the uniqueness of zero or 1
# post: for each fusion, check the left and right sides in the list, and get the total length of the long read support, and also get the uniqueness in the sequences near by the fusion events.
#       this will be the shortest value or lowest unique count for either side of the fusion
#       we get: shortest total length, lowest total unique base count, and lowest flanking unique count for various distances proximal to the fusion point
def get_fusion_uniqueness(fusion_gpds, fusion_coordinates):
  output = {}
  for fid in fusion_gpds:
    #get uniqueness of proximal 10, 50, 100
    proxcheck = [10, 50, 100] # these are the sizes to check
    proxuniqcount = {}
    shortest_length = 9999999999
    shortest_unique = 9999999999
    for line in fusion_gpds[fid]:
      [g, side] = line
      total_base_count = 0
      unique_base_count = 0
      unique_array = []
      for i in range(0,g['exonCount']):
        for j in range(g['exonStarts'][i],g['exonEnds'][i]):
          total_base_count += 1
          if g['chrom'] in fusion_coordinates:
            coord = j+1
            if coord in fusion_coordinates[g['chrom']]:
              unique_array.append(fusion_coordinates[g['chrom']][coord])
              if fusion_coordinates[g['chrom']][coord] == 1:
                unique_base_count+=1
            else:
              unique_array.append(0)
          else:
            unique_array.append(0)
      #now for this line we have our unique counts
      # save the shortest
      if total_base_count < shortest_length: shortest_length = total_base_count
      if unique_base_count < shortest_unique: shortest_unique = unique_base_count
      if side == '-': #fusion is on the right side
        unique_array.reverse() # reverse in place so the first elements are proximal to the junction
      for dist in proxcheck:
        sum = 0
        maxdist = dist
        if maxdist > len(unique_array): maxdist = len(unique_array)
        for k in range(0,maxdist): 
          sum += unique_array[k]
        if dist not in proxuniqcount:
          proxuniqcount[dist] = int(sum)
        elif sum < proxuniqcount[dist]:
          proxuniqcount[dist] = int(sum)  # set it to the lower value
    #Now we have our shortest results to report back
    o = {}
    o['shortest_total_length'] = shortest_length
    o['shortest_unique_count'] = shortest_unique
    o['lowest_proximal_unique_counts'] = proxuniqcount
    output[fid] = o  
  return output    

#pre: uniqueness_bedgraph is a bed graph format file with unique entries labled 1
#     fusion_coordiantes is initialized to all values equal to zero it is a dictionary keyed by chromosome then 1 indexed coordinate
#post: returns nothing, but it modifies fusion_coordinates to have the uniqueness score
#
def get_coordinate_uniqueness(uniqueness_bedgraph_filename, fusion_coordinates):
  with open(uniqueness_bedgraph_filename) as infile:
    for line in infile:
      f = line.rstrip().split("\t")
      chrom = f[0]
      start = int(f[1])
      finish = int(f[2])
      uniqueness = float(f[3])
      if uniqueness != 1: continue  # we are only interested in cataloging unique entries right now
      for i in range(start, finish):
        if chrom in fusion_coordinates:
          coord = i+1
          if coord in fusion_coordinates[chrom]:
            fusion_coordinates[chrom][coord] = uniqueness
  return

# pre: the per fusion dictionary of genepreds 
# post: a dictionary keyed by chromosome then 1-indexed coordinate initalized to zero
def initialize_fusion_coordinates(fusion_gpds):
  d = {}
  for fid in fusion_gpds:
    for line in fusion_gpds[fid]:
      entry = line[0]
      if entry['chrom'] not in d:
        d[entry['chrom']] = {}  
      for i in range(0,entry['exonCount']):
        for j in range(entry['exonStarts'][i],entry['exonEnds'][i]):
          d[entry['chrom']][j+1] = 0
  return d

# We are going to make an attempt at making LR_fusion.gpd
# We will try to keep the name thats being used the same except for
# changing the length and the fusion point. 
# It looks like there is some smooothing done to ignore indels
# where blocks are combined.
# pre: a <psl file>
#      query name is of the format:
#        A_B|C[+/-]D|ccs[+/-]E
#        A: identity
#        B: length
#        C: F# for fusion id
#        +/-: Indicates whether the fusion occurs at the start of the alignment or end
#        D: pos where the junction occurs
#           NOTE:  This is what we really need to know what this format actually is so we can recode position correctly.  We know where the true fusion point is.  
#                  And in this function we pass in the 1-indexed base of the coordinate of the last base of the fragment proximal to the fusion event
#        +/-: optional
#        E: optional
#      <min gap in block size>  this i used to smooth out the genepred file a bit so blocks are more like exons
#        gaps smaller than this get joined together.
#      <output fusion genepred filename>
# post: write a genepred file
#       the gene_name is taken from qName in the psl file.
#       the name is updated to have the length of genepred entry, and
#       the new fusion position.
#       save the genepred entries in a dictionary keyed by fusion
#       that has first and second entries in an array, and each of those entries has a genepred entry and a +(left) -(right) side the fusion is on
def make_fusion_genepred(psl_filename,min_gap_in_block_size,outfile):
  i = 0
  ofile = open(outfile,'w')
  fusion_gpd = {}
  with open(psl_filename) as infile:
    for line in infile:
      m = i % 2
      i += 1
      entry = psl_basics.read_psl_entry(line)
      gline = psl_basics.convert_entry_to_genepred_line(entry)
      e = genepred_basics.genepred_line_to_dictionary(gline)
      # m == 0 is left, m == 1 is right
      look = re.search('\|(F[\d]+)([+-])',e['gene_name'])
      fid = look.group(1)
      fsidesign = look.group(2) #which side is the fusion on

      e['cdsStart'] = e['txStart']        
      e['cdsEnd'] = e['txEnd']        
      
      f = genepred_basics.smooth_gaps(e,min_gap_in_block_size)
 
      # get new transcript length
      newlen = 0
      for j in range(0,f['exonCount']): newlen += f['exonEnds'][j] - f['exonStarts'][j]

      # make a new name based on the new fusion site and new length
      look = re.match('^([\d\.]+_)\d+(\|F\d+[+-].*)$',f['gene_name'])
      v1 = look.group(1)
      v2 = str(newlen)
      v3 = look.group(2)
      f['gene_name'] = v1 + v2 + v3

      gpline = genepred_basics.genepred_entry_to_genepred_line(f)
      ofile.write(gpline+"\n")
      if fid not in fusion_gpd: 
        n = [0,0]
      fusion_gpd[fid] = n
      fusion_gpd[fid][m] = [f, fsidesign] #send back the genepred and the side the fusion is on
  ofile.close()
  return fusion_gpd

def opposite(sign):
  if sign == '+': return '-'
  return '+'

  








main()
