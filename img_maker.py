import HTSeq
import sys
import numpy as np
import draw_script

# output formats for images:
OUTPUT_FORMATS = ['png', 'pdf']


def get_gtf_reader(gtf_path):
  """
  Takes a path to a GTF (a string) and returns a HTSeq.GFF_Reader object
  """
  return HTSeq.GFF_Reader(gtf_path)


def find_features(gtf_path, transcript_id):
  """
  Takes a species (String) and a transcript_id (e.g. ENST0000...) and returns the Set of HTSeq.GenomicFeature's corresponding
  to this transcript that are in the GTF for this species.  Identifies the transcript by the 'transcript_id' attribute in the 9th field of the GTF   
  """
  gtf = get_gtf_reader(gtf_path)
  feature_set = set()
  try:
    for feature in gtf:
      if feature.attr['transcript_id'] == transcript_id:
        feature_set.add(feature)
    return feature_set
  except IOError:
    sys.exit("Exception!  Could not find a GTF file at: "+str(gtf_path))


def check_chromosome_naming(bam_file):
  try:
    bam = HTSeq.BAM_Reader(bam_file)
    return 'chr' == bam.get_header_dict()['SQ'][0]['SN'][:3]
  except IOError:
    sys.exit("Exception!  Could not find the BAM file located at: "+str(bam_file))


def create_window(features, buffer_left, buffer_right, bam_has_chr_prefix):
  """
  Creates and returns a HTSeq.GenomicInterval in which to search for aligning reads.
  Does this by searching through the set of genomic features for the desired transcript and locating the 'transcript' type (which gives the length of this transcript)
  Adds on a buffer on either side depending on the buffer_* args.  
  """
  for f in features:
    if f.type == "transcript":
      chr = str(f.iv.chrom)
      if bam_has_chr_prefix and chr[:3] != 'chr':
        chr = 'chr'+str(chr)
      window = HTSeq.GenomicInterval(chr, f.iv.start-buffer_left, f.iv.end+buffer_right, '.')
      return window
  print "Could not locate a GenomicFeature of type 'transcript', so could not create a GenomicInterval."
  sys.exit(1)


def get_coverage_profile(bam_file, window):
  """
  Reads through the BAM file and returns a numpy array giving the coverage for the transcript.  
  """
  try:
    bam = HTSeq.BAM_Reader(bam_file)
    coverage = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    try:
      for alnmt in bam[window]:
        if alnmt.aligned:
          coverage[alnmt.iv]+=1
    except ValueError:
      sys.exit("""Exception when reading the BAM file.
                  This is common for two situations:
                  1: There is no .bai file for your BAM file  
                  2: It is possible that your BAM file and GTF file do not have the same chromosome names.  
                     Check the chromosome names in the sam or index file and those in the GTF for agreement (and fix as necessary).""")  
    #we now have coverage, which is a generator for tuples.  
    #Each tuple has a GenomicInterval and an integer for the read-depth.  To eventually plot, we need to make this into a numpy array
    cvg_list = []
    it = coverage.steps() #an iterator
    try:
      step = it.next() #get the first object from the iterator so we can enter the while loop
      while step:
        if step[0].start<=window.end and step[0].end>=window.start:  #if step overlaps with window
          if step[0].start<=window.start and step[0].end>=window.end: #if the step is longer than the window (unlikely, but possible)
            cvg_list.extend(window.length*[step[1]])
          elif step[0].start<=window.start:
            cvg_list.extend(abs(step[0].end-window.start)*[step[1]])            
          elif step[0].end>=window.end:
            cvg_list.extend(abs(window.end-step[0].start)*[step[1]])
          else:
            cvg_list.extend(step[0].length*[step[1]])         
        step = it.next()
    except StopIteration: #when the generator is done, it throws an exception, which we catch and move on
      pass
    if len(cvg_list) == 0:
      sys.exit("Could not find any coverage data for the genomic region: "+str(window)+" in BAM file: "+str(bam_file))
    cvg = np.array(cvg_list)
    return cvg
  except IOError:
    sys.exit("Could not locate the BAM file at "+str(bam_file))


def main():
  """
  Script can be invoked with 5 args:
    1: the BAM file (NEEDS TO BE SORTED/INDEXED!)
    2: path to a GTF file
    3: transcript_id (ENST id or other that is in the GTF file)
    4: buffer_left (how many bp to the left of the transcript should we plot)
    5: buffer_right (how many bp to the right of the transcript should we plot)
    6: output format for the figure (png, pdf)
  """

  if len(sys.argv) == 7:
    bam_file = sys.argv[1]
    gtf_path = sys.argv[2]
    transcript_id = sys.argv[3]
    try:
      buffer_left = int(sys.argv[4])
      buffer_right = int(sys.argv[5])
    except ValueError:
      sys.exit("Could not parse the left and/or right buffers as integers.")
    
    save_format = sys.argv[6]
    if save_format.lower() not in OUTPUT_FORMATS:
      sys.exit("Please specify a save format from: "+ str(OUTPUT_FORMATS))

    features = find_features(gtf_path, transcript_id)
    has_chr_prefix = check_chromosome_naming(bam_file) #whether the BAM file has chromosomes named like 'chr12' or just '12'
    window = create_window(features, buffer_left, buffer_right, has_chr_prefix)
    cvg = get_coverage_profile(bam_file, window)

    draw_script.make_gene_plot(cvg, features, buffer_left, save_format)

  else:
    print "Please look at the input arguments to launch the script correctly."
    sys.exit(1)


if __name__ == "__main__":
  main()
  
