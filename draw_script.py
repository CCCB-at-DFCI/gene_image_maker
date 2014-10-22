import sys

try:
  import HTSeq
  import matplotlib.pyplot as plt
  from matplotlib.patches import Rectangle, Polygon
  import numpy as np
except ImportError as ex:
  sys.exit("Failed to import a required dependency: "+str(ex.message))

# features that we are including in the gene model picture:
# need to match those in the GTF (the third column)
EXON='exon'
UTR='UTR'
INTRON='intron'

# features that are retained from parsing the GTF
target_features = [EXON, UTR]

class Config():
  """
  This class acts as a container for some parameters controlling the appearance of the plot (e.g. colors, relative sizes)
  """
  def __init__(self):
    self.GENE_MODEL_COLOR = 'cadetblue'
    self.COVERAGE_COLOR = 'darkgrey'
    self.GENE_MODEL_HEIGHT = 0.3 #the fractional height of the gene model image in relation to the coverage plot
    self.GENE_MODEL_FEATURE_HEIGHTS = {EXON: 1.0, UTR: 0.5, INTRON: 0.2} #relative heights of the features w.r.t GENE_MODEL_HEIGHT
    self.VERTICAL_OFFSET = 0.55 #distance between the coverage plot x-axis and gene model, relative to the height of the coverage plot
    self.LEFT_START = 0 #pixels shown to the left of the plot's origin
    self.ASPECT_RATIO = 0.1 #controls how wide vs. tall the coverage plot is.  Higher numbers compress the x-axis and make the y-axis taller
    self.DIRECTION_MARKER_SIZE=0.5 #what percentage of the exon height the width of the direction triangle should be
    self.DIRECTION_MARKER_RATIO=1.4 #controls the aspect ratio of the directional marker.  Higher number means a longer triangle
    self.DIRECTION_MARKER_COLOR='white'


  def setup(self, max_depth, region_size, gene):
    self.ASPECT_RATIO = self.ASPECT_RATIO/float(max_depth/float(region_size))
    self.GENE_MODEL_HEIGHT = self.GENE_MODEL_HEIGHT*max_depth
    self.GENE_MODEL_FEATURE_HEIGHTS={k:v*self.GENE_MODEL_HEIGHT for k,v in self.GENE_MODEL_FEATURE_HEIGHTS.iteritems()}
    self.VERTICAL_OFFSET = self.VERTICAL_OFFSET*max_depth
    self.DIRECTION_MARKER_SIZE=(0.5*self.DIRECTION_MARKER_SIZE)/(np.sqrt(3)/2.0)
    self.DIRECTION_MARKER_VERTICES = self.generate_direction_marker(gene.get_strand() == "-")
    self.DIRECTION_MARKER_HEIGHT = (self.DIRECTION_MARKER_VERTICES[:,0]-self.DIRECTION_MARKER_VERTICES[:,1])[0]
    self.DIRECTION_MARKER_VERTICES[0,:] = self.DIRECTION_MARKER_VERTICES[0,:]+(self.DIRECTION_MARKER_HEIGHT/2.0-self.DIRECTION_MARKER_VERTICES[0,0]) #adjust so that the midpoint of the triangle's height falls at the origin

  def generate_direction_marker(self, is_on_negative_strand):
    directions = np.array([0, 120, 240]) #angles for a right-pointing triangle

    # rotate if gene is on negative strand:
    if is_on_negative_strand:
      directions = directions + 60

    # setup the directions to the vertices, then add vertices appropriately scaled to produce a reasonable looking triangle
    directions = np.radians(directions)
    unit_vectors = np.array([np.cos(directions),np.sin(directions)])
    vertices =  (self.DIRECTION_MARKER_SIZE*self.GENE_MODEL_FEATURE_HEIGHTS[EXON])*unit_vectors
    scaling_mtx = np.diag([self.DIRECTION_MARKER_RATIO*self.ASPECT_RATIO, 1]) #for stretching lengthwise so that it appears as a triangle on the stretched figure space
    return scaling_mtx.dot(vertices)
    


class GeneModel():
  """
  This class holds the data and some logic methods for the gene model
  """

  def __init__(self):
    self.features=[]


  def add_features(self, all_features):
    #takes a set of GenomicFeature's
    self.features = list(all_features)

  @staticmethod
  def break_intervals(fA, fB):
    """
    In the case where an exon and UTR annotation overlap their intervals (which is true in general), 
    we break this overlap into disjoint segments.  
    Each resulting segment is assigned an identity (e.g. is it a UTR region or not) so that
    we can plot the UTR and exons with different heights if desired.

    fA and fB are HTSeq.GenomicFeature objects and are sorted by start and end positions
    """
    locs = []
    locs.append(fA.iv.start)
    locs.append(fB.iv.start)
    if fA.iv.end < fB.iv.end:
      locs.append(fA.iv.end)
      locs.append(fB.iv.end)
    else:
      locs.append(fB.iv.end)
      locs.append(fA.iv.end)
    
    potential_intervals = [HTSeq.GenomicInterval(fA.iv.chrom, locs[i], locs[i+1], fA.iv.strand) for i in range(3)]
    tuples = [(i, (fA.iv.overlaps(i), fB.iv.overlaps(i))) for i in potential_intervals if i.length>0]
    new_intervals = []
    for gi, state in tuples:
      if state == (1,1):
        type = UTR
      elif state == (1,0):
        type = fA.type
      elif state == (0,1):          
        type = fB.type
      new_intervals.append(HTSeq.GenomicFeature('',type, gi))
    return sorted(new_intervals, key=lambda f: (f.iv.start, f.iv.end))


  @staticmethod
  def edit_for_overlaps(features):
    """
    Takes a sequence of HTSeq.GenomicFeature's and edits so that the features have disjoint intervals.  
    This allows plotting of UTRs, which would otherwise by covered by the overlapping exon feature
    """
    edited_intervals = []
    last_feature = None
    for i in range(len(features)-1):
      if not last_feature:
        fA = features[i]
      else:
        fA = last_feature
      fB = features[i+1]
      if fA.iv.overlaps(fB.iv) and fA.type != fB.type:
        new_intervals = GeneModel.break_intervals(fA,fB)
        edited_intervals.extend(new_intervals[:-1])
        last_feature = new_intervals[-1]
      else:
        edited_intervals.append(fA)
        last_feature = None
    edited_intervals.append(fB)
    return sorted(edited_intervals, key=lambda f: (f.iv.start, f.iv.end))
            

  def get_structure(self):
    """
    Returns a sorted sequence of HTSeq.GenomicFeature's so the gene model can be plotted with the coverage plot
    Filters for features contained in the target_feature list
    The intervals for the resulting GenomicFeatures are disjoint (non-overlapping) so that UTR's and exons can be plotted
    """
    filtered_features = filter(lambda f: f.type in target_features, self.features)
    sf = sorted(filtered_features, key=lambda f: (f.iv.start, f.iv.end)) #sort by start position, then by end. Thus if we get two intervals with the same start, the shorter one comes first
    sf = GeneModel.edit_for_overlaps(sf)
    return sf


  def get_gene_length(self):
    for f in self.features:
      if f.type == 'transcript':
        return abs(f.iv.start-f.iv.end)


  def get_strand(self):
    return self.features[0].iv.strand


  def get_gene_name(self):
    try:
      return self.features[0].attr['gene_name']
    except KeyError:
      print "Could not locate gene name.  Using gene_id field instead..."
      return self.features[0].attr['gene_id']



def add_gene_feature_to_plot(ax, start_x, width, type, config, is_on_negative_strand):
  """
  Draws gene features (exon, intron, etc) onto the plot.
  Args are: 
    -ax: axes object (the axes in which we're drawing), 
    -start_x: a double (in axes/genomic coordinates) indicating where the feature starts
    -width: a double giving the length of the feature (in axes/genomic coords)
    -type: a string indicating the type of feature (exon, intron, etc).  
    -config: Config object which has the plot parameters
    -is_on_negative_strand: boolean indicating the strand
  """

  y_start = -config.VERTICAL_OFFSET
  try:
    height = config.GENE_MODEL_FEATURE_HEIGHTS[type]
  except KeyError:
    sys.exit('Could not draw feature of type <'+str(type)+'> since it has not been assigned a size.  Needs to be in the set: '+str(config.GENE_MODEL_FEATURE_HEIGHTS.keys()))

  # add the feature:  
  start_pt = (start_x, y_start-height/2.0)
  ax.add_patch(Rectangle(start_pt, width, height, edgecolor = config.GENE_MODEL_COLOR, facecolor=config.GENE_MODEL_COLOR, linewidth=0))
  
  # if we have a exon, draw the transcription direction:
  if type == EXON and config.DIRECTION_MARKER_HEIGHT < width*0.75:
    center_pt = np.array([[start_x+width/2.0], [y_start]]) # (2x1) vector
    vertices = config.DIRECTION_MARKER_VERTICES + center_pt
    ax.add_patch(Polygon([x for x in vertices.transpose()], facecolor=config.DIRECTION_MARKER_COLOR, edgecolor=config.DIRECTION_MARKER_COLOR, linewidth=0))
  return ax


def create_plot(cvg_depth, config):
  """
  Creates the figure, plots the coverage depth.
  Input args are:
    - cvg_depth: numpy array of the coverage
    - config: a Config object for setting some of the plot details/paramters
  """

  fig = plt.figure()
  ax1 = fig.add_subplot(1,1,1)
  max_depth = np.max(cvg_depth)
  width = cvg_depth.shape[0]

  # customize the axes:
  ax1.set_yticks([max_depth-1]) #ticks need to be less than the actual data max for them to show up
  ax1.get_xaxis().tick_bottom()
  ax1.get_yaxis().tick_left()
  ax1.spines['top'].set_visible(False)
  ax1.spines['bottom'].set_position('zero')
  ax1.spines['left'].set_bounds(0,max_depth*1.1)
  ax1.spines['right'].set_visible(False)
  ax1.set_aspect(config.ASPECT_RATIO)

  # plot the coverage profile:
  x_coord = np.arange(0, cvg_depth.shape[0])  
  line = ax1.fill_between(x_coord, cvg_depth, 0)
  line.set_facecolor(config.COVERAGE_COLOR)
  line.set_edgecolor(config.COVERAGE_COLOR)
  
  ax1.set_ylim(bottom=-(config.VERTICAL_OFFSET+config.GENE_MODEL_HEIGHT))
  ax1.set_xlim(right=width)

  font = {'family':'serif', 'size':8, 'weight':'light'}
  plt.rc("font", **font)

  return fig  


def plot_gene_model(fig, gene_model, buffer_left, config):
  """
  Adds the gene model to the plot, aligned appropriately with the coverage profile
  """
  ax = fig.gca()

  strand = gene_model.get_strand()
  marker = None
  gene_structure = gene_model.get_structure()
  abs_start = gene_structure[0].iv.start-buffer_left  # the start position for the window we're looking at (in genomic coords).
  for feature in gene_model.get_structure():
    start = feature.iv.start-abs_start
    end = feature.iv.end-abs_start
    type = feature.type
    if type == EXON:
      ax = add_gene_feature_to_plot(ax, start, abs(start-end), type, config, gene_model.get_strand() == "-")   
    elif type == UTR:
      ax = add_gene_feature_to_plot(ax, start, abs(start-end), type, config, gene_model.get_strand() == "-")   
    else:
      sys.exit('Unrecognized feature')

    # have to also plot any introns (the gene model does not track these)
    if marker and marker<start:
      ax = add_gene_feature_to_plot(ax, marker, abs(start-marker), INTRON, config, gene_model.get_strand() == "-")   
    marker = end


def save_figure(fig, name, format):
  plt.savefig(str(name)+"."+str(format), format = format,  bbox_inches='tight')
  

def create_gene_model(all_features):
  gene = GeneModel()
  gene.add_features(all_features)
  return gene


def add_gene_name(fig, gene_name):
  ax = fig.gca()
  ax.set_title(gene_name)


def make_gene_plot(cvg, feature_array, buffer_left, format):
  """
  Function that is called from external script
  Given a coverage array (numpy array) and a set of HTSeq.GenomicFeature's.
  """
  config = Config()
  region_size = len(cvg)
  max_depth = np.max(cvg)
  gene = create_gene_model(feature_array)
  config.setup(max_depth, region_size, gene)
  fig = create_plot(cvg, config)
  plot_gene_model(fig, gene, buffer_left, config)
  gene_name = gene.get_gene_name()
  add_gene_name(fig, gene_name)
  save_figure(fig, gene_name, format)

