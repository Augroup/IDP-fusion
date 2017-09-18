#!/usr/bin/env python
import sys, re
import genepred_basics, sequence_basics

# Pre: <reference genepred>
#      <genome fasta>
#      <junctions.txt.accurate>
#      <distance to combine> if zero or empty then don't do any combining
# Post: output report to standard output
#      Fields:
#        <junction id> - unqiue for each junction
#        <left chromosome>
#        <left coordinate> - 1-indexed base coordinate proximal to the fusion
#        <left sign> - "+" means fusion is on the right side "-" means left for fusion
#        <End of reference exon> - "Y" means it is the last base of a reference exon
#        <right chromosome>
#        <right coordinate> - 1-indexed base coordinate proximal to the fusion
#        <right sign> - "-" means fusion is on the right side "+" means left for fusion
#        <Start of reference exon> - "Y" means it is the last base of a reference exon
#        <Fusion description> - the two fusion genes separated by a dash.  these may or may not be ordered properly.  see next three fields
#        <Transcript strands> - for those genes in fusion description, which strand are they on
#        <Evidence for correct 5-3prime ordering -transcript> - the transcript orientations were
#          used to arrange the left and right ordering of the sides because they appeared informative "Y"
#          "N" indicates transcript directions were either uninformative or in disagreement with the orientations of the fusion
#        <Evidece for correct 5-3prime ordering -SpliceJunctionBases> - similar to the above.  it states whether the flanking intronic bases support
#          these as being potential canoncical splice junctions.  In absence of exonic evidence, cononical splice site junctions will be used
#          to order the left and right if they are "Y" present.
#        <Gene 1 name> -left (hopefully 5' if ordered properly)
#        <Gene 2 name> - right
#        <LR count> - long reads supporting this
#        <nr Type1> - non redundant unique short read fusion support.  
#          This unqiue means that each counted read maps to only one place among fusions
#          or reference sequences, and that is in support of this junction.  non redudnant
#          means that if multiple differet reads map to this exact same location it is only counted once here
#        <Type1> - count of short reads that are uniquely mapped to this fusion junction.  It means
#          They map to one place and one place only in the whole of the reference or fusion sequences and it is here.
#        <Type2> - Multimapped fusion reads.  These reads map to mutliple locations, but those multiple locations are all in support of different fusion events, not reference
#        <Type3> - Reference Unique - these reads map to a fusion, but they also map to a unique location in the reference sequence
#        <Type4> - Reference multimapped - these reads map to a fusion, but they also map to multiple locations in the reference genome
#        <Max-min LR side length> - For each long read supporting the fusion, take the shorter of the left or right side.  Then from all those shortest sides return the longest.
#        <Max-min unique base count> - For each long read supporting the fusion, 
#          take the number of unqiuely mappable bases on the left and right sides and return the count for the side with the least
#          then return the largest of those numbers for all long reads
#        <Max-min flanking unqiue base counts> - Similar to the above count but this is over several windows of distances
#          starting at the fusion point and going out.  so rather than looking at the entire long read, we are only considering the nearest bases.
#          if it is 10:0,50:20 it means that no bases within ten base pairs were
#          unqiue but 20 base pairs within 50 were unique.

def main():
  if len(sys.argv) < 4:
    print sys.argv[0] + ' <reference genepred> <genome fasta> <junctions.txt.accurate> <distance to combine>'
    sys.exit()
  reference_genepred_filename = sys.argv[1]
  genome_filename = sys.argv[2]
  junctions_filename = sys.argv[3]
  combine_distance = 0
  if len(sys.argv) == 5:
    combine_distance = int(sys.argv[4])
  juncs = read_junctions_and_group_long_reads(junctions_filename)
  splicebases = get_splice_bases(juncs,genome_filename)
  annot = annotate_junctions(juncs,reference_genepred_filename,splicebases)
  if combine_distance > 0:
    newannot = combine_nearby(annot,combine_distance)
  else:
    for id in annot:
      print annot[id]

def combine_nearby(annot,combine_distance):
  ids = annot.keys()
  poolofsets = set()
  for i in range(0,len(ids)):
    for j in range(i+1,len(ids)):
      d1 = parse_id(ids[i])
      d2 = parse_id(ids[j])
      dist = max(distance(d1,d2),distance(d2,d1)) # do they match
      if dist < combine_distance and dist != -1:
        newset = set()
        newset.add(ids[i])        
        newset.add(ids[j])
        poolofsets.add("\t".join(newset))
      else:
        poolofsets.add(ids[i])
        poolofsets.add(ids[j])
  #now combine sets until we have a minimal sized set
  numsets = len(poolofsets)
  prevnumsets = 0
  while(prevnumsets != len(poolofsets)):
    #print len(poolofsets)
    prevnumsets = len(poolofsets)
    try_to_combine(poolofsets)
  #print "done combining"
  newannot = {}
  for subset in poolofsets:
    entrylist = subset.split("\t")
    entry = entrylist[0]
    updated_line = annot[entry]
    for entry2 in range(1,len(entrylist)):
      #print "one " + entrylist[entry2] + " " + updated_line 
      updated_line = add_two_lines(updated_line,annot[entrylist[entry2]])
      #print "two " + entrylist[entry2] + " " + updated_line
    print updated_line

def add_two_lines(l1,l2):
  f1 = l1.split("\t")
  f2 = l2.split("\t")
  for i in range(15,21):
    if f1[i] != '.' and f2[i] != '.':
      f1[i] = str(int(f1[i])+int(f2[i]))
    elif f1[i] == '.':
      f1[i] = f2[i]
  for i in range(21,23):
    f1[i] = str(max(int(f1[i]),int(f2[i])))
  val1 = f1[23].split(',')
  val2 = f2[23].split(',')
  nval = []
  for i in range(0,len(val1)):
    res1 = val1[i].split(":")
    res2 = val2[i].split(":")
    nval.append(res1[0]+":"+str(max(int(res1[1]),int(res2[1]))))
  f1[23] = ",".join(nval)
  return "\t".join(f1)

def try_to_combine(poolofsets):
  setlist = list(poolofsets)
  for i in range(0,len(setlist)):
    for j in range(i+1,len(setlist)):
      if setlist[i] in poolofsets and setlist[j] in poolofsets:
        if should_combine(setlist[i],setlist[j]):
          newset = set()
          for e1 in setlist[i].split("\t"):
            newset.add(e1)
          for e2 in setlist[j].split("\t"):
            newset.add(e2)
          poolofsets.remove(setlist[i])
          poolofsets.remove(setlist[j])
          poolofsets.add("\t".join(newset))
  return        

def should_combine(s1,s2):
  for e1 in s1.split("\t"):
    if e1 in s2.split("\t"):
      return True
  return False

def distance(d1,d2):
  if d1['chr1'] != d2['chr1'] or d1['chr2'] != d2['chr2']:  return -1
  if d1['dir1'] != d2['dir1'] or d1['dir2'] != d2['dir2']:  return -1
  dist = max(abs(d1['coo1']-d2['coo1']),abs(d1['coo2']-d2['coo2']))
  return dist

def parse_id(id):
  m = re.match('^([^:]+):([\d]+)([+-])/([^:]+):(\d+)([+-])$',id)
  if not m: 
    print "problem parsing id "+id
    sys.exit()
  d = {}
  d['chr1'] = m.group(1)
  d['coo1'] = int(m.group(2))
  d['dir1'] = m.group(3)
  d['chr2'] = m.group(4)
  d['coo2'] = int(m.group(5))
  d['dir2'] = m.group(6)
  return d

# pre: get the junction dictionary by id
#      get the genome
# post: return the flanking intronic bases on the left and right side
#       if strand on the left is on the negative side, reverse complement the bases
#       if strand on the right is on the negative side, reverse complement the bases
#       This way 'GT AG' indicates a normal splice site junction
#       'CT AC' indicates the orientation should be reversed

def get_splice_bases(junctions,genome_filename):
  genome = sequence_basics.read_fasta_into_hash(genome_filename)
  bases = {}
  for id in junctions:
    chr1 = junctions[id]['chr1']
    coo1 = int(junctions[id]['coo1'])
    dir1 = junctions[id]['dir1']
    chr2 = junctions[id]['chr2']
    coo2 = int(junctions[id]['coo2'])
    dir2 = junctions[id]['dir2']
    bases1 = '??'
    if dir1 == '+' and len(genome[chr1]) > coo1+1:
      bases1 = genome[chr1][coo1] + genome[chr1][coo1+1]
    elif dir1 == '-' and len(genome[chr1]) > coo1-2:
      bases1 = sequence_basics.rc(genome[chr1][coo1-3] + genome[chr1][coo1-2])
    bases2 = '??'
    if dir2 == '+' and len(genome[chr2]) > coo2-2:
      bases2 = genome[chr2][coo2-3] + genome[chr2][coo2-2]
    elif dir2 == '-' and len(genome[chr2]) > coo2+1:
      bases2 = sequence_basics.rc(genome[chr2][coo2] + genome[chr2][coo2+1])
    bases[id] = bases1.upper() + " " + bases2.upper()
  return bases

# pre: junction dictionary by id
#      reference genepred file
#      splicebases dictionary by id with flanking, possibly intronic, bases
# post: print the following fields
#       <systematic id of junction>
#       <left chromosome>
#       <left coordinate>
#       <left direction>
#       <right chromosome>
#       <right coordinate>
#       <right direction>
#       <fusion name>
#       <strands of anntotated fusion genes>
#       <evidence for orientation of fusion based on reference transcripts>
#       <evidence for orientation of fusion based on connical splice site "GU-AG">
#       <left fusion gene>
#       <right fusion gene>
#       <long read count>
#       <non-redundant Type 1>  Unique + nonredundant fusion
#       <Type 1> Unique fusion
#       <Type 2> Multimapped fusion
#       <Type 3> Supports fusion but also uniquely maps to reference
#       <Type 4> Supports fusion but also multimaps to reference

def annotate_junctions(junctions,reference_genepred_filename, splicebases):
  annot_struct = genepred_basics.get_gene_annotation_data_structure(reference_genepred_filename)
  genepred = genepred_basics.get_per_chromosome_array(reference_genepred_filename)
  result = {}
  for junc_name in junctions:
    splice = splicebases[junc_name]
    junc = junctions[junc_name]
    genepred_entries = []
    if junc['chr1'] in genepred:  genepred_entries = genepred[junc['chr1']]
    leftend = 'N'  # are we on the end of an exon on the left side
    if junc['dir1'] == '+' and genepred_basics.is_exon_end(junc['chr1'],int(junc['coo1']),genepred_entries):
       leftend = 'Y'
    elif junc['dir1'] == '-' and genepred_basics.is_exon_start(junc['chr1'],int(junc['coo1']),genepred_entries):
       leftend = 'Y'

    genepred_entries = []
    if junc['chr2'] in genepred:  genepred_entries = genepred[junc['chr2']]
    rightend = 'N'
    if junc['dir2'] == '+' and genepred_basics.is_exon_start(junc['chr2'],int(junc['coo2']),genepred_entries):
      rightend = 'Y'
    elif junc['dir2'] == '-' and genepred_basics.is_exon_end(junc['chr2'],int(junc['coo2']),genepred_entries):
      rightend = 'Y'
    annot1 = genepred_basics.gene_annotate_by_coordinates(junc['chr1'],int(junc['coo1'])-1,int(junc['coo1']),annot_struct)
    annot2 = genepred_basics.gene_annotate_by_coordinates(junc['chr2'],int(junc['coo2'])-1,int(junc['coo2']),annot_struct)
    best1 = junc['chr1']+":"+str(junc['coo1'])
    best1dir = ''
    best2 = junc['chr2']+":"+str(junc['coo2'])
    best2dir = ''
    name1 = ''
    name2 = ''
    if annot1:
      best1 = annot1[0][0]
      name1 = annot1[0][0]
      best1dir = annot1[0][1]
    if annot2:
      best2 = annot2[0][0]
      name2 = annot2[0][0]
      best2dir = annot2[0][1]
    evidence_transcript = 'N'
    evidence_splice = 'N'
    action = 'unknown'
    if best1dir == '+' and best2dir == '+':
      if junc['dir1'] == '-' and junc['dir2'] == '-':
        action= "reverse"
      elif junc['dir1'] == '+' and junc['dir2'] == '+':
        action = "no reverse"
      else:
        action = "unknown"
    elif best1dir == '-' and best2dir == '-':
      if junc['dir1'] == '+' and junc['dir2'] == '+':
        action = "reverse"
      elif junc['dir1'] == '-' and junc['dir2'] == '-':
        action = "no reverse"
      else:
        action = "unknown"
    if best1dir == '+' and best2dir == '-':
      if junc['dir1'] == '-' and junc['dir2'] == '+':
        action = "reverse"
      elif junc['dir1'] == '+' and junc['dir2'] == '-':
        action = "no reverse"
      else:
        action = "unknown"
    if best1dir == '-' and best2dir == '+':
      if junc['dir1'] == '+' and junc['dir2'] == '-':
        action = "reverse"
      elif junc['dir1'] == '-' and junc['dir2'] == '+':
        action = "no reverse"
      else:
        action = "unknown"
    if best1dir == '-' and best2dir == '':
      if junc['dir1'] == '+':
        action = "reverse"
      elif junc['dir1'] == '-':
        action = "no reverse"
      else:
        action = "unknown"
    if best1dir == '+' and best2dir == '':
      if junc['dir1'] == '-':
        action = "reverse"
      elif junc['dir1'] == '+':
        action = "no reverse"
      else:
        action = "unknown"
    if best1dir == '' and best2dir == '+':
      if junc['dir2'] == '-':
        action = "reverse"
      elif junc['dir2'] == '+':
        action = "no reverse"
      else:
        action = "unknown"
    if best1dir == '' and best2dir == '-':
      if junc['dir2'] == '+':
        action = "reverse"
      elif junc['dir2'] == '-':
        action = "no reverse"
      else:
        action = "unknown"

    #set if we have transcript evidence
    if action != "unknown":
      evidence_transcript = 'Y'

    # set if we have splice evidence
    if action == "no reverse" and splice == 'GT AG':
      evidence_splice = 'Y'
    if action == "reverse" and splice == 'CT AC':
      evidence_splice = 'Y'
    if action == "unknown" and (splice == 'GT AG' or splice == 'CT AC'):
      evidence_splice = 'Y'

    # if action is still unknown, use splice site to make choice to reverse
    if action == "unknown":
      if splice == 'GT AG':
        action = "no reverse"
      elif splice == 'CT AC':
        action = "reverse"

    # switch our left and right end assignments if reverse
    if action == "reverse":
      temp1 = leftend
      leftend = rightend
      rightend = temp1

    names = ''
    if action == "no reverse":
      names = junc['chr1'] + "\t" + junc['coo1'] + "\t" + junc['dir1'] + "\t" + leftend + "\t" + junc['chr2'] + "\t" + junc['coo2'] + "\t" + junc['dir2'] + "\t" + rightend + "\t" + best1 + "-" + best2 + "\t" + best1dir + "/" + best2dir + "\t"+evidence_transcript+"\t"+evidence_splice+"\t"+name1 + "\t" + name2
    elif action == "reverse":
      names = junc['chr2'] + "\t" + junc['coo2'] + "\t" + opposite(junc['dir2']) + "\t" + leftend + "\t" + junc['chr1'] + "\t" + junc['coo1'] + "\t" + opposite(junc['dir1']) + "\t" + rightend+ "\t" +best2 + "-" + best1 + "\t" + best2dir + "/" + best1dir + "\t"+evidence_transcript+"\t"+evidence_splice+"\t" + name2 + "\t" + name1
    else:
      names = junc['chr1'] + "\t" + junc['coo1'] + "\t" + junc['dir1'] + "\t" + leftend + "\t" +junc['chr2'] + "\t" + junc['coo2'] + "\t" + junc['dir2'] + "\t" + rightend + "\t" + best1 + "-" + best2 + "\t" + best1dir + "/" + best2dir + "\t"+evidence_transcript+"\t"+evidence_splice+"\t" + name1 + "\t" + name2
    result[junc_name] = junc_name + "\t" + names + "\t" + str(junc['lr_count']) + "\t" + junc['type1nr'] + "\t" + junc['type1'] + "\t" + junc['type2'] + "\t" + junc['type3'] + "\t" + junc['type4'] + "\t" + str(junc['max_min_side_lengths']) + "\t" + str(junc['max_min_unique_side_counts']) + "\t" + junc['max_min_flanking_distance']
  return result

# Pre: a "junctions.txt.accurate" type filename
# Post: a dictionary of junction events with the number of long reads
#       supporting that junction event and the various counts of short read 
#       evidence.
#       The chromosomes, coordinates, and directions of the sides of the
#       fusion are included too
def read_junctions_and_group_long_reads(junctions_filename):
  juncs = {}
  with open(junctions_filename) as infile:
    for line in infile:
      f = line.rstrip().split("\t")
      m = re.match('^([\S]+)_[\d]+_[\d]+_([+-])/([\S]+)_[\d]+_[\d]+_([+-])',f[2])
      if not m:
        print 'malformed entry in junctions file '+"\n"+f[2]
        sys.exit()
      v = {}
      v['chr1'] = m.group(1)
      v['dir1'] = m.group(2)
      v['coo1'] = int(f[3])
      v['coo1true'] = v['coo1']
      if v['dir1'] == '-': v['coo1true'] += 1
      v['chr2'] = m.group(3)
      v['dir2'] = m.group(4)
      v['coo2'] = int(f[4])
      v['coo2true'] = v['coo2']
      if v['dir2'] == '+': v['coo2true'] += 1
      #coo1true and coo2true contain the 'true' last base before the
      #breakpoints (index-1).  These can be used to interrogate the genepred
      #file because they will tell us what gene the breakpoint is in.
      # dir1 + means the break point of the first sequence is on the right side
      # dir2 + means the break point of the second sequence is on the left side
      uniq = v['chr1'] + ':' + str(v['coo1true'])+v['dir1'] +"/"+ v['chr2'] + ':' + str(v['coo2true'])+v['dir2']
      uniqrc = v['chr2'] + ':' + str(v['coo2true'])+opposite(v['dir2'])+"/"+v['chr1']+':'+str(v['coo1true'])+opposite(v['dir1'])
      vals = [uniq, uniqrc]
      vals.sort()
      mine = vals[0]  # use the lower alphabetical order version for directionless fusion
      s = {} # choose the sorted type of this
      m = re.match('^([\S]+):([\d]+)([+-])/([\S]+):([\d]+)([+-])$',mine)
      s['chr1'] = m.group(1)
      s['coo1'] = m.group(2)
      s['dir1'] = m.group(3)
      s['chr2'] = m.group(4)
      s['coo2'] = m.group(5)
      s['dir2'] = m.group(6)
      s['longreads'] = set()
      s['type1nr'] = f[5]
      s['type1'] = f[6]
      s['type2'] = f[7]
      s['type3'] = f[8]
      s['type4'] = f[9]
      s['min_side_lengths'] = set()
      s['min_unique_side_counts'] = set()
      s['flanking_distance_keys'] = []
      s['flanking_distance'] = {}
      fpair = f[15].split(",")
      for v1 in fpair:
        [keynum, keyval] = v1.split(':')
        s['flanking_distance_keys'].append(int(keynum))
        s['flanking_distance'][int(keynum)] = set()
      s['flanking_distance_unique_counts'] = {}
      # if we want to find the best long read support this would be a good place to do it
      if mine not in juncs:
        juncs[mine] = s
      juncs[mine]['longreads'].add(f[0])
      juncs[mine]['min_side_lengths'].add(int(f[13]))
      juncs[mine]['min_unique_side_counts'].add(int(f[14]))
      for v1 in fpair:              
        [keynum, keyval] = v1.split(':')
        s['flanking_distance'][int(keynum)].add(int(keyval))
  for mine in juncs:
    juncs[mine]['lr_count'] = len(juncs[mine]['longreads'])
    juncs[mine]['max_min_side_lengths'] = max(juncs[mine]['min_side_lengths'])
    juncs[mine]['max_min_unique_side_counts'] = max(juncs[mine]['min_unique_side_counts'])
    opairs = []
    for num in juncs[mine]['flanking_distance_keys']:
      opairs.append(str(num)+":"+str(min(juncs[mine]['flanking_distance'][num])))
    juncs[mine]['max_min_flanking_distance'] = ",".join(opairs)
  return juncs

def opposite(dir):
  if dir != '-' and dir != '+': return dir
  if dir == '+':
    return '-'
  return '+'

main()
