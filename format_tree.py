#parse gene info and store as dic to build NCBI taxonomy labelled tree
from ete3 import Tree, PhyloTree, TreeStyle , NodeStyle , RectFace , AttrFace, faces , TextFace, CircleFace
from Bio import Entrez
import pickle
from colour import Color
import uniprot as uni
from Bio import AlignIO
import taxa
import glob
import pickle

#this program loads a phylogenetic tree based on an alignment of homologous proteins
#it shows the the evolution of a set of proteins or genes in the form of a binary tree
#I use phyml to calculate phylogenies but it works with any tree
# the formating and color scheme for species is identical to the detection tree 
# divergence from the taxonomic tree indicates important evolutionary events like duplications or losses.



#load a tree and associated alignment
#treefile = '/home/cactuskid/Dropbox/IIB/mergeLineages/yuyo/restricted/hapAndeff/strcutres_and_tcoffeeset_aln_struct.phy_phyml_tree.txtlabels.txt'

folder = '/home/cactuskid/Dropbox/IIB/mergeLineages/yuyo/*/*labels.txt'
#folder = '/home/cactuskid/Dropbox/IIB/mergeLineages/yuyo/*/*/*labels.txt'

treefiles = glob.glob(folder)
for treefile in treefiles:
	print treefile
	colorSepcies = False
	#alg = '/home/cactuskid/Dropbox/IIB/mergeLineages/phylogeny/hybrid/merged_curate_aln.fasta'
	t = PhyloTree( treefile, sp_naming_function=None) #, alignment=alg, alg_format="fasta")
	# Calculate the midpoint node
	R = t.get_midpoint_outgroup()
	# and set it as tree outgroup
	t.set_outgroup(R)

	def save_obj(obj, name ):
	    with open( name + '.pkl', 'wb') as f:
	        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

	def load_obj(name ):
	    with open( name + '.pkl', 'r') as f:
	        return pickle.load(f)

	genedict = load_obj('genedict')
	speciescolors = load_obj('colors')
	red = Color('red')
	blue = Color('blue')
	colorvec = list(red.range_to(blue, len(genedict)))
	colormap = {}
	columnmap = {}
	for i,fasta in enumerate(genedict):
		columnmap[fasta] = i
		colormap[fasta] = colorvec[i].hex
	annotated = [] 
	print speciescolors
	for fasta in genedict:
		for leaf in t.get_leaves():

			nst = NodeStyle()
			nst["size"] = 0
			nst["fgcolor"] = 'black'
			nst["hz_line_width"] = 2
			nst["vt_line_width"]= 2
			nst.show_name = True
			if leaf.name.split('/')[0] in genedict[fasta]:
				if 'HH' not in fasta and 'LOMETS' not in fasta: 
					leaf.add_face( RectFace ( 10 , 10 , colormap[fasta], colormap[fasta] ), column = columnmap[fasta] )
					if leaf not in annotated:
						try:
							face = leaf.add_face( TextFace ( text = genedict[fasta][leaf.name.split('/')[0]][2]) , column = 10  )
							annotated.append(leaf)
						except:
							print genedict[fasta][leaf.name.split('/')[0]]	
					
					if colorSepcies == True:
						lineage = genedict[fasta][leaf.name.split('/')[0]][3].split(';') + [genedict[fasta][leaf.name.split('/')[0]][2]] + genedict[fasta][leaf.name.split('/')[0]][2].split()
						
						for taxa in lineage[-1:]:
							if taxa in speciescolors:
								print taxa
								print speciescolors[taxa]
								nst['bgcolor'] = speciescolors[taxa]	
								break
						else:
							#fall back on approximate coloring...
							print leaf.name
							for key in  speciescolors:
								genusSpecies = key.split()
								for name in genusSpecies:
									if name in lineage:
										nst['bgcolor'] = speciescolors[key]	
										break
					leaf.img_style = nst 

	nst = NodeStyle()
	nst["size"] = 0
	nst["fgcolor"] = 'black'
	nst["hz_line_width"] = 2
	nst["vt_line_width"]= 2
	for n in t.traverse():
		if n.is_leaf() == False:
			n.img_style = nst 

	ts = TreeStyle()
	for i,ref in enumerate(colormap.keys()):
		if 'ubi' not in ref and 'HH' not in ref and 'LOMETS' not in ref:
			ts.title.add_face(TextFace(ref, fsize=12), column=0)
			ts.title.add_face( RectFace(10 , 10 , colormap[ref] , colormap[ref]), column = 1)
	#ts.mode = "c"
	#ts.arc_start = -180 #0 degrees = 3 o'clock
	ts.scale =  150 
	ts.show_leaf_name = False
	ts.show_branch_support = True
	#t.show(tree_style = ts)
	outpng = treefile+'outnocolor.png'
	t.render(outpng , tree_style = ts)
