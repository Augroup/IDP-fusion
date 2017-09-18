#!/usr/bin/env python
import sys, re
import genepred_basics, sequence_basics

# Pre: <a genepred file>
#      <output of make_fusion_report.py>
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
# Post: 
#   The same format report as described above, but with annotations added for
#   fusions based on another genepred file (here designed to be an ensembl based)

def main():
  if len(sys.argv) != 3:
    print sys.argv[0] + ' <reference genepred> <report file>'
    sys.exit()
  reference_genepred_filename = sys.argv[1]
  report_filename = sys.argv[2]

  juncs = read_report(report_filename)
  annotate_junctions(juncs,reference_genepred_filename)

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

# pre: report array
#      reference genepred file
# post: print the new report

def annotate_junctions(junctions,reference_genepred_filename):
  annot_struct = genepred_basics.get_gene_annotation_data_structure(reference_genepred_filename)
  genepred = genepred_basics.get_per_chromosome_array(reference_genepred_filename)
  for entry_set in junctions:
    f = entry_set['entries']
    d = entry_set['dictionary']
    if d['Gene_1'] != '' and d['Gene_2'] != '':
      print "\t".join(f)  #already annotated.  move on
      continue

    [best1dir,best2dir] = d['Transcript_strands'].split("/")

    leftend = d['End_of_reference_exon_1']  # are we on the end of an exon on the left side
    name1 = d['Gene_1']
    best1 = name1
    name1change = 0
    if d['Gene_1'] == '':
      best1 = d['Chromosome_1']+":"+str(d['Coordinate_1'])
      genepred_entries = []    
      if d['Chromosome_1'] in genepred:  
        genepred_entries = genepred[d['Chromosome_1']]
        if d['Direction_1'] == '+' and genepred_basics.is_exon_end(d['Chromosome_1'],int(d['Coordinate_1']),genepred_entries):
          leftend = 'Y'
        elif d['Direction_1'] == '-' and genepred_basics.is_exon_start(d['Chromosome_1'],int(d['Coordinate_1']),genepred_entries):
          leftend = 'Y'
        annot1 = genepred_basics.gene_annotate_by_coordinates(d['Chromosome_1'],int(d['Coordinate_1'])-1,int(d['Coordinate_1']),annot_struct)
        if annot1:
          name1change = 1
          best1 = annot1[0][0]
          name1 = annot1[0][0]
          best1dir = annot1[0][1]

    rightend = d['Start_of_reference_exon_2']
    name2 = d['Gene_2']
    best2 = name2
    name2change = 0
    if d['Gene_2'] == '':
      best2 = d['Chromosome_2']+":"+str(d['Coordinate_2'])
      genepred_entries = []
      if d['Chromosome_2'] in genepred:  
        genepred_entries = genepred[d['Chromosome_2']]
        if d['Direction_2'] == '+' and genepred_basics.is_exon_start(d['Chromosome_2'],int(d['Coordinate_2']),genepred_entries):
          rightend = 'Y'
        elif d['Direction_2'] == '-' and genepred_basics.is_exon_end(d['Chromosome_2'],int(d['Coordinate_2']),genepred_entries):
          rightend = 'Y'
        annot2 = genepred_basics.gene_annotate_by_coordinates(d['Chromosome_2'],int(d['Coordinate_2'])-1,int(d['Coordinate_2']),annot_struct)
        if annot2:
          name2change = 1
          best2 = annot2[0][0]
          name2 = annot2[0][0]
          best2dir = annot2[0][1]

    if name1change == 0 and name2change == 0:
      print "\t".join(f)
      continue # didn't find any additional name information

    evidence_transcript = 'N'
    evidence_splice = 'N'
    action = 'unknown'
    junc = {}
    junc['dir1'] = d['Direction_1']
    junc['dir2'] = d['Direction_2']
    junc['coo1'] = d['Coordinate_1']
    junc['coo2'] = d['Coordinate_2']
    junc['chr1'] = d['Chromosome_1']
    junc['chr2'] = d['Chromosome_2']

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

    # switch our left and right end assignments if reverse
    if action == "reverse":
      temp1 = leftend
      leftend = rightend
      rightend = temp1

    names = ''
    if action == "no reverse":
      names = junc['chr1'] + "\t" + str(junc['coo1']) + "\t" + junc['dir1'] + "\t" + leftend + "\t" + junc['chr2'] + "\t" + str(junc['coo2']) + "\t" + junc['dir2'] + "\t" + rightend + "\t" + best1 + "-" + best2 + "\t" + best1dir + "/" + best2dir + "\t"+evidence_transcript+"\t"+evidence_splice+"\t"+name1 + "\t" + name2
    elif action == "reverse":
      names = junc['chr2'] + "\t" + str(junc['coo2']) + "\t" + opposite(junc['dir2']) + "\t" + leftend + "\t" + junc['chr1'] + "\t" + str(junc['coo1']) + "\t" + opposite(junc['dir1']) + "\t" + rightend+ "\t" +best2 + "-" + best1 + "\t" + best2dir + "/" + best1dir + "\t"+evidence_transcript+"\t"+evidence_splice+"\t" + name2 + "\t" + name1
    else:
      names = junc['chr1'] + "\t" + str(junc['coo1']) + "\t" + junc['dir1'] + "\t" + leftend + "\t" +junc['chr2'] + "\t" + str(junc['coo2']) + "\t" + junc['dir2'] + "\t" + rightend + "\t" + best1 + "-" + best2 + "\t" + best1dir + "/" + best2dir + "\t"+evidence_transcript+"\t"+evidence_splice+"\t" + name1 + "\t" + name2
    print d['Alphabetized_junction_name'] + "\t" + names + "\t" + str(d['LR_count']) + "\t" + str(d['nr_Type1_count']) + "\t" + str(d['Type1_count']) + "\t" + str(d['Type2_count']) + "\t" + str(d['Type3_count']) + "\t" + str(d['Type4_count']) + "\t" + str(d['Max_min_side_length']) + "\t" + str(d['Max_min_side_unique']) + "\t" + d['Max_min_side_flank']

# Pre: a report type filename
# Post: An array of entries
def read_report(report_filename):
  juncs = []
  with open(report_filename) as infile:
    for line in infile:
      if re.match('^#',line): continue
      d = idpfusion_line_to_dictionary(line)
      f = line.rstrip().split("\t")
      val = {}
      val['dictionary'] = d
      val['entries'] = f
      juncs.append(val)
  return juncs

def opposite(dir):
  if dir != '-' and dir != '+': return dir
  if dir == '+':
    return '-'
  return '+'


# Pre:  IDP-fusion output file line
# Post:  A dictionary with the still 1-indexed coordinates and evidence
#        preceeding the fusion point.
def idpfusion_line_to_dictionary(line):
  f = line.rstrip().split("\t")
  for i in range(16,21):
    if f[i] == '.': f[i] = 0
  d = {}
  d['Alphabetized_junction_name'] = f[0]
  d['Chromosome_1'] = f[1]
  d['Coordinate_1'] = f[2]
  d['Direction_1'] = f[3]
  d['End_of_reference_exon_1'] = f[4]
  d['Chromosome_2'] = f[5]
  d['Coordinate_2'] = f[6]
  d['Direction_2'] = f[7]
  d['Start_of_reference_exon_2'] = f[8]
  d['Fusion_description'] = f[9]
  d['Transcript_strands'] = f[10]
  d['Evidence_for_correct_5-3prime_ordering-Transcript'] = f[11]
  d['Evidence_for_correct_5-3prime_ordering-SpliceJunctionBases'] =f[12]
  d['Gene_1'] = f[13]
  d['Gene_2'] = f[14]
  d['LR_count'] = int(f[15])
  d['nr_Type1_count'] = int(f[16])
  d['Type1_count'] = int(f[17])
  d['Type2_count'] = int(f[18])
  d['Type3_count'] = int(f[19])
  d['Type4_count'] = int(f[20])
  d['Max_min_side_length'] = int(f[21])
  d['Max_min_side_unique'] = int(f[22])
  d['Max_min_side_flank'] = f[23]
  return d

main()
