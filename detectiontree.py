#parse gene info and store as dic to build NCBI taxonomy labelled tree
from ete3 import Tree, PhyloTree
from ete3 import TreeStyle , NodeStyle , RectFace , AttrFace, faces , TextFace, CircleFace 
import pickle
from colour import Color
from Bio import AlignIO, SeqIO
import taxa
from Bio import Entrez
import ujson as json
from csb.bio.io.hhpred import HHOutputParser
import uniprot as uni
import glob
import numpy as np
from sklearn.manifold import MDS
from multiprocessing import Pool
import unicodedata


def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'r') as f:
        return pickle.load(f)




seqfiles = '../mergeLineages/treefiles/'
Entrez.email = "dmoi@iibintech.com.ar"
#download iformation
dl = False
# use ubiquitin data to make a complete speicies tree
useBGTree = True
#cut uninteresting clades
cut_clades = False
root = 'Eukaryota'
dead_branches = ['Fungi' , 'Chordata' ]
includenames = [ 'Eukaryota','Eutheria' , 'Chordata' , 'Methateria' , 'Arthropoda' , 'Nematoda' , 'Viridiplantae', 'Alveolata', 'Eumetazoa', 'Fungi', 'Amoebozoa', 'Heterolobosea' , 'Opisthokonta' ,  'Rhizaria' , 'Kinetoplastida' , 'Rhodophyta' , 'Stramenopiles' , 'Choanoflagellida']
outputspecies = [ 'Eukaryota', 'Fungi' , 'Metazoa'] 

#restrict results displayed
restrictfiles = False
allowedFiles = ['']
tidyOrphans = False
filename = 'hap2tree'
trim = 40
minLevel = 4

subtrees = [ 'Arthropoda' , 'Viridiplantae' ]

if dl == True:
	org_dict={}
	genedict = {}
	print 'downloading taxonomy data'
	if useBGTree == True:
		#all clades
		taxlist = []
		namelist = []
		with open('../phylofiles/species/refs/ubiquitin.txt','r') as info:
			for i,line in enumerate(info):
				if i>0:
					taxid = line.split('	')[-2]
					taxlist.append(taxid)
			taxlist = set(taxlist)
		#save all organisms lineage and genome info
		save_obj(taxa.grabGenomes(taxlist), 'bgtaxa')

	
	#all possible sequences to be detected
	fastas = glob.glob(seqfiles + '*.fa*')
	print fastas
	print 'loading species info'
	for fasta in fastas:
		names , Lineages, genes = taxa.get_taxinfo(fasta)
		genedict[fasta] = genes

	#3rd dataset
	#all sequences detected
	detection = load_obj('detection')
	for results in detection.keys():
		names, lineages, genes = taxa.uniprotIDlist_Lineages( detection[results] ) 
		#eliminate redundant entries
		taxlist = set(taxlist)
		genedict[results] = genes

	print 'download finished'
	save_obj(genedict, 'genedict')


def grab_IDs(seqfile):
    record_iterator = SeqIO.parse(seqfile, "fasta")
    idlist = []
    for prot in record_iterator:
    	code = prot.id
    	if '_' in prot.id:
    		code = prot.id.split('_')[0]
        if '|' in prot.id:
        	code = prot.id.split('|')[1]
        idlist.append(code)
    return idlist



##### use ncbi lineages to generate a speices tree 
def create_tree(filename, genedict,tax_dict , root = 'Eukaryota' , completeGenomes = True):
	print 'making tree'
	t = Tree()
	for ref in genedict:
		subdict = genedict[ref]
		for protcode in subdict:
			try:
				genomeLink = subdict[protcode][4]
				if genomeLink != 'noGenome' or completeGenomes == False:
					lineage = subdict[protcode][3]
					words = lineage.split(';')
					node = t
					if root in lineage:
						for level,taxa in enumerate(words):
							createNew = True
							for c in node.children:
								if c.name == taxa:
									createNew = False
									node = c
									array = json.loads(node.refs)
									if ref not in array:
										array.append(ref)
									node.refs= json.dumps(array)
									array = json.loads(node.codes)
									if protcode not in array:
										array.append(protcode)
									node.codes = json.dumps(array)
							if createNew == True:
								newnode = node.add_child(name=taxa)
								node = newnode
								node.add_features( codes = json.dumps([protcode]) )
								node.add_features( refs = json.dumps([ref]) )
			except:
				print protcode
				print subdict[protcode]
	
	#load the rest without protcodes if you want a global taxonomic tree
	print 'making background tree'
	for ref in tax_dict:
		genomeLink = tax_dict[ref][2]
		if completeGenomes == False or  genomeLink != 'noGenome':
			lineage = tax_dict[ref][1]
			if root in lineage:	
				words = lineage.split(';')
				node = t
				for level,taxa in enumerate(words):
					createNew = True
					for c in node.children:
						if c.name == taxa:
							createNew = False
							node = c
							break
					if createNew == True:
						print 'adding new node'
						print taxa
						newnode = node.add_child(name=taxa)
						node = newnode
						node.add_features( codes = json.dumps([]) )
						node.add_features( refs = json.dumps([]) )
		else :
			print ref
			print tax_dict[ref]
	print t
	print 'tree loaded '
	t.write(format=1, outfile= filename)
	return t


###use isomap to embed species tree in 2d color space
def tree_to_speciescolors(t,distmat = None, columndict = None):
	if distmat == None or columndict == None:
		columndict, distmat = create_distmat(t)
		save_obj(columndict, 'columndict')
		save_obj(distmat, 'distmat')
	proj = create_2dprojection(distmat)
	colors = proj_tocolordict(proj,columndict,t)
	save_obj(colors , 'colors')
	return colors

def get_dist(args):
	n,m = args
	return n.get_distance(m, topology_only = True)

def create_distmat(t):
	#creates a rough distance matrix between species to be used for coloring purposes in other applications
	#since this is only a topological tree the distances aren't actual evolutionary distances...
	distmat = np.zeros((len(t.get_leaves()), len(t.get_leaves())))
	column_dict = {}
	print 'creating distance matrix'
	print len(t.get_leaves())
	print 'species'
	jobs = []
	coords = []
	for i,n in enumerate(t.get_leaves()):
		column_dict[n.name] = i
		for j, m in enumerate(t.get_leaves()):
			if i < j :
				coords.append((i,j))
				jobs.append((n,m))

	print len(jobs)
	print 'getting distances'
	pool = Pool()
	results = pool.map_async(get_dist,jobs).get()
	print 'DONE'

	for k, coords in enumerate(coords):
		i,j = coords
		distmat[i,j] = results[k]
	distmat = distmat + distmat.T
	print 'DONE'
	return column_dict, distmat

def create_2dprojection(distmat):
	#uses isomap to return a species distance map in 2d based on the topological distmat of all species in tree
	print 'map to 3d space'
	mapper=MDS(n_components=3, metric=True, n_init=4, max_iter=300, verbose=0, eps=0.001, n_jobs=-1, random_state=0, dissimilarity='precomputed')
	projmat =mapper.fit_transform(distmat)
	print 'DONE'
	return projmat

def proj_tocolordict(projmat, columndict , t):
	#map columns to 0-1
	print 'creating color dictionary'
	colors = {}
	for i in range(projmat.shape[1]):
		projmat[:,i] -= np.amin(projmat[:,i])
		projmat[:,i] /= np.amax(projmat[:,i])
		
	for name in columndict: 
		i=columndict[name]
		c = Color(rgb = (projmat[i,0], projmat[i,1], projmat[i,2]) ).hex
		#define hue and saturation based on 2d mapping
		colors[name] = c

	# assign avg colors to all upstream clades
	for j,n in enumerate(t.traverse()):
		if n.name not in colors:
			rgb = np.zeros(3)
			print n.name
			rgb = 0
			i = 0
			for l in n.get_leaves():
				if l.name in colors:
					rgb += np.asarray(Color(colors[l.name]).rgb)
					i +=1
				else:
					for key in  colors:
						genusSpecies = key.split()
						for name in genusSpecies:
							if name in l.name:
								rgb += np.asarray(Color(colors[key]).rgb)
								i += 1
								break
			if i !=0:
				print i
				rgb /= i
				colors[n.name] = Color(rgb = rgb).hex
	print 'DONE'
	return colors

##final formating and output
def format( tree , genedict , detection , includenames , filename , speciescolors = None):
	print 'final output...'
	red = Color('red')
	blue = Color('blue')
	colorvec = list(red.range_to(blue, len(genedict.keys())))
	colormap = {}
	columnmap = {}
	for i,ref in enumerate(genedict):
		if ref not in detection:	
			columnmap[ref] = 3
			if 'hybrid' in ref.lower():
				columnmap[ref] = 0
			if 'eff' in ref.lower():
				columnmap[ref] = 1
			if 'hap' in ref.lower():
				columnmap[ref] = 2
			colormap[ref] = colorvec[i].hex
	for i,ref in enumerate(detection):
		columnmap[ref] = 3 + i
		colormap[ref] = colorvec[i].hex 

	print columnmap
	print colormap

	circledict = {}
	for n in t.traverse():
		
		nst = NodeStyle()
		nst["size"] = 0
		nst["fgcolor"] = 'black'
		nst["hz_line_width"] = 4
		nst["vt_line_width"]= 4
		nst.show_name = False
		if n.is_leaf():
			if speciescolors != None and n.name in speciescolors:
				nst["bgcolor"] = colors[n.name]
			nst.show_name = True			
			n.add_face( AttrFace(attr = 'name', ftype='Helvetica', fgcolor='black', fsize =18 ,fstyle = 'normal'   ), column =0 )
			refs = json.loads(n.refs)
			for ref in genedict:
				if ref in refs and ref in detection:
					n.add_face( CircleFace ( 10 , colormap[ref]), column =  2 + columnmap[ref] )
					n.img_style = nst
				if ref in refs and ref not in detection:
					n.add_face( RectFace ( 20 , 20 , colormap[ref], colormap[ref] ), column =  2 + columnmap[ref] )
					n.img_style = nst
				if ref not in refs and ref not in detection:
					n.add_face( RectFace ( 20 , 20 , colormap[ref], 'white' ), column =  2 + columnmap[ref] )
					n.img_style = nst
			###color by species
			if n.name in speciescolors:
				nst['bgcolor'] = speciescolors[n.name]
		else:
			if n.name.strip() in includenames:
				n.add_face( AttrFace(attr = 'name', ftype='Helvetica', fgcolor='black', fsize =20 ,fstyle = 'normal'   ), column =0 )
				nst.size = 2
				n.img_style = nst
			else:
				nst.size = 0
				n.img_style = nst 
	ts = TreeStyle()
	for i,ref in enumerate(colormap.keys()):
		if 'ubi' not in ref:
			ts.title.add_face(TextFace(ref, fsize=12), column=0)
			ts.title.add_face( RectFace(10 , 10 , colormap[ref] , colormap[ref]), column = 1)
	ts.show_leaf_name=False 

	"""ts.mode = "c"
				ts.arc_start = 270
				ts.arc_span = 359
				ts.root_opening_factor = 1
			"""
	ts.scale =  190
	
	t.show(tree_style = ts)
	t.render(filename + ".png", tree_style = ts)
	t.render( filename +".svg" , tree_style = ts)

def prune_tree(cut_clades, dead_branches ,trim , minLevel, tidyOrphans):
	#topological tree
	print 'pruning'
	#cut off uninteresting clades
		
	prunevec = []
	for n in t.traverse():
		prunevec.append(n)

	# remove layers
	descendantsup = [t]
	descendants = []
	prunevec = []
	for i in range(trim):
		print 'trimming ' + str(i)
		for node in descendantsup:
			for d in node.children:
				descendants.append(d)
				prunevec.append(d)
		descendantsup = descendants
		descendants = []	
	t.prune(prunevec)
	print 'DONE'
	

	print dead_branches
	if cut_clades == True:
		for n in t.traverse():
			if n.name.strip() in dead_branches:
				print 'cutting ' + n.name
				for d in n.get_descendants():
					prunevec.remove(d)
		t.prune(prunevec)

	#tidy up monoclade branches...
	if tidyOrphans == True:
		keepPruning = True
		while(keepPruning):
			keepPruning = False
			for n in t.get_leaves():
				if len(n.up.children) == 1 and n.name.strip() not in dead_branches :
					keepPruning = True
					prunevec.remove(n)
			t.prune(prunevec)
	print 'done pruning'

def output_list(species,genedict):
	#output csv of taxonomic information for each sequence in specific clades
	print species
	for specie in species:
		for filename in genedict:
			if 'fasta' in filename:
				output = []
				retrievelist = []
				genelist = genedict[filename]
				print len(genelist)
				for gene in genelist:
					if len(genelist[gene])>3:
						if specie in genelist[gene][3]:
							if 'noGenome'  not in genelist[gene]:
								if gene.strip() not in retrievelist:
									retrievelist.append(gene.strip())
				outfile = filename+specie+'.fasta'			
				record_iterator = SeqIO.parse(filename, "fasta")
				for prot in record_iterator:
					protcode = prot.id
					if '|' in protcode:
						protcode = protcode.split('|')[1]
					if '_' in protcode:
						protcode = prot.id.split('_')[0]
					if ' ' in protcode:
						protcode = prot.id.split(' ')[0]
					#if prot.id in categories[category]:
					if protcode in retrievelist:
						if protcode not in output:
							output.append(prot)
				print 'filtered for ' + specie + 'found nseqs='
				print len(output)
				handle = open(outfile, 'w')
				SeqIO.write(output, handle, "fasta")
				handle.close()

##### run the functions #####	
tax_dict = load_obj('bgtaxa')
genedict = load_obj('genedict')
detection = {}
try:
	detection = load_obj('detection')
except:
	pass

print len(tax_dict)
print tax_dict.values()[0]
#copute colors from precalculated distmat
#it's prob best to calculate the distmat on a cluster since the size is nspecies**2 
t = create_tree(filename,genedict,tax_dict)
colors = tree_to_speciescolors(t , distmat = load_obj('distmat') , columndict = load_obj('columndict'))
output_list(outputspecies, genedict )
prune_tree(cut_clades, dead_branches, trim , minLevel , tidyOrphans)
format( t , genedict ,detection , includenames , filename , colors)

