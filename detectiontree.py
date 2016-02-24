#parse gene info and store as dic to build NCBI taxonomy labelled tree
from ete3 import Tree, PhyloTree, TreeStyle , NodeStyle , RectFace , AttrFace, faces , TextFace, CircleFace
from Bio import SeqIO
import pickle
from colour import Color
from Bio import AlignIO
import taxa
from Bio import Entrez
import ujson as json

import glob

def save_obj(obj, name ):
    with open( name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name ):
    with open( name + '.pkl', 'r') as f:
        return pickle.load(f)

treefile = 'haptree.phy'
seqfiles = '../database_cluster/hap2ectoin/'
Entrez.email = "dmoi@iibintech.com.ar"
org_dict = {}
lineage = {}
orglookup = {}
#download iformation
dl = False 
#cut uninteresting clades
cut_clades = False
dead_branches = ['Fungi' , 'Chordata']
includenames = [ 'Eutheria' , 'Chordata' , 'Methateria' , 'Arthropoda' , 'Nematoda' , 'Viridiplantae', 'Alveolata', 'Eumetazoa', 'Fungi', 'Amoebozoa', 'Rhizaria' ]
#restrict results displayed
restrictfiles = False
allowedFiles = ['']
tidyOrphans = True
filename = 'hap2tree'
trim = 5

if dl == True:
	print 'downloading taxonomy data'
	with open('detection_taxlist.txt', 'w') as taxout:
		#all clades
		taxlist = []
		namelist = []
		with open('../phylofiles/species/refs/ubiquitin.txt','r') as info:
			for i,line in enumerate(info):
				if i>0:
					taxid = line.split('	')[-2]
					taxlist.append(taxid)
					taxout.write(taxid + '\n')
			taxlist = set(taxlist)
			record = Entrez.read(Entrez.efetch(db="taxonomy", id=list(taxlist)))
			for rec in record:
				namelist.append(rec['Lineage'])
		org_dict['../phylofiles/species/refs/ubiquitin.txt'] = namelist
		protlist = []
		

		#all possible sequences to be detected
		seqfiles = glob.glob(seqfiles + '*.fasta')
		for fasta in seqfiles:
			protlist = []
			namelist = []
			print fasta
			record_iterator = SeqIO.parse(fasta, "fasta")
			for prot in record_iterator:
				protcode = prot.id
				if '|' in protcode:
					protcode = protcode.split('|')[1]
				if '_' in protcode:
					protcode = prot.id.split('_')[0]
				protlist.append(protcode)
			names , Lineages, genes = taxa.uniprotIDlist_Lineages(protlist)

			org_dict[fasta] = Lineages
		
		#3rd dataset
		#all sequences detected
		detection = load_obj('detection')
		for results in detection.keys():
			names, lineages, genes = taxa.uniprotIDlist_Lineages( detection[results] ) 
			#eliminate redundant entries
			taxlist = set(taxlist)
			if 'model' in results:
				if 'LOMETS' in org_dict.keys():
					org_dict['LOMETS'] += lineages
				else :
					org_dict['LOMETS'] = lineages
			else :
				org_dict[results] = lineages

		#add manually detected stuff without uniprot codes here
		taxlist = ['505711','33657','461836','342808','221722']
		for taxid in taxlist:
			taxout.write(taxid + '\n')	
		print taxlist
		record = Entrez.read(Entrez.efetch(db="taxonomy", id=list(taxlist)))
		for rec in record:
			namelist.append(rec['Lineage'])
		print namelist
		org_dict['LOMETS'] += namelist
	print 'download finished'
	save_obj(orglookup, 'org_dict_positive_protcode')
	save_obj(org_dict, 'org_dict_detection')


def create_tree(org_dict,filename):
	print 'making tree'
	t = Tree()
	for ref in org_dict:
		for lineage in org_dict[ref]:
			words = lineage.split(';')
			node = t
			for level,taxa in enumerate(words):
				createNew = True
				for c in node.children:
					if c.name == taxa:
						createNew = False
						node = c
						string = node.refs
						array = json.loads(string)
						array.append(ref)
						node.refs= json.dumps(array)
				if createNew == True:
					newnode = node.add_child(name=taxa)
					node = newnode
					node.add_features( refs = json.dumps([refs]))
	print t
	print 'tree loaded '
	t.write(format=1, outfile= filename)
	return t


def format( tree , org_dict ,includenames , filename ):
	print 'final output...'
	red = Color('red')
	blue = Color('blue')

	colorvec = list(red.range_to(blue, len(refs)))
	colormap = {}
	columnmap = {}

	for i,ref in enumerate(refs):
		print ref
		if 'LOMETS' in ref or 'pdb' in ref:
			columnmap[ref] = 2
		else:
			columnmap[ref] = 0
		colormap[ref] = colorvec[i].hex
	circledict = {}
	for n in t.traverse():
		nst = NodeStyle()
		nst["size"] = 0
		nst["fgcolor"] = 'black'
		nst["hz_line_width"] = 2
		nst["vt_line_width"]= 2
		nst.show_name = False
		#n.add_face( AttrFace(attr = 'name', ftype='Helvetica', fgcolor='black', fsize =10 ,fstyle = 'normal'   ), column =0 )
		if n.is_leaf():
			nst.show_name = True			
			for ref in org_dict.keys():
				alreadynoted=[]
				for lineage in org_dict[ref]:
					if 'ubi' not in ref and n not in alreadynoted and  ref in json.loads(n.refs):	
						alreadynoted.append(n)
						n.add_face( RectFace ( 10 , 10 , colormap[ref], colormap[ref] ), column =  2 + columnmap[ref] )
						n.img_style = nst
				else:
					n.img_style = nst
		else:
				if n.name.strip() in includenames:
					n.add_face( AttrFace(attr = 'name', ftype='Helvetica', fgcolor='black', fsize =10 ,fstyle = 'normal'   ), column =0 )
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
	#ts.show_leaf_name=False 
	t.show(tree_style = ts)
	t.render(filename + ".png")
	t.render(filename +".svg")

org_dict = load_obj('org_dict_detection')
orglookup = load_obj( 'org_dict_positive_protcode')
refs = org_dict.keys()

def prune_tree(cut_clades, trim , tidyOrphans):
	#topological tree
	print 'pruning'
	#cut off uninteresting clades
		
	prunevec = []
	for n in t.traverse():
		prunevec.append(n)

	# remove layers
	for i in range(trim):
		for leaf in t.get_leaves():		
			prunevec.remove(leaf)		
		t.prune(prunevec)
	
	if cut_clades ==True:
		for n in t.traverse():
			if n.name in dead_branches:
				for d in n.get_descendants():
					try:
						prunevec.remove(d)
					except:
						pass
		t.prune(prunevec)

	#tidy up the single children...
	if tidyOrphans == True:
		keepPruning = True
		while(keepPruning):
			keepPruning = False
			for n in t.get_leaves():
				if len(n.up.children) == 1 and n.name.strip() not in dead_branches :
					print 'single child'
					print n.name
					keepPruning = True
					prunevec.remove(n)
			t.prune(prunevec)
	print 'done pruning'

t = create_tree(org_dict,filename)
prune_tree(cut_clades, trim , tidyOrphans)
format( t , org_dict ,includenames , filename)

