#utility functions for getting taxonomic data

from Bio import Entrez
from Bio import SeqIO
import uniprot as uni
from itertools import islice , chain
Entrez.email = 'dmoi@iibintech.com.ar'

def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]


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

def grab_seqs(seqfile,codes):
    record_iterator = SeqIO.parse(seqfile, "fasta")
    prots = []
    for prot in record_iterator:
    	code = prot.id
    	if '_' in prot.id:
    		code = prot.id.split('_')[0]
        if '|' in prot.id:
        	code = prot.id.split('|')[1]
        if code in codes:
        	prots.append(prot)
    return prots

def filter_lineages(Fasta, ancestor):
	# ancestor must be an NCBI taxa
	protcodes = grab_IDs(Fasta)
	for i, codes in enumerate(chunks(protcodes , 100) ):
		print codes
		newnames, newlineages, newgenes = uniprotIDlist_Lineages(codes)
		if i == 0:
			names, Lineages, genes = (newnames, newlineages, newgenes)
		names += newnames
		Lineages += newlineages
		for gene in newgenes:
			genes[gene] = newgenes[gene]
	finalcodes = []

	for gene in genes:
		if len(genes[gene]) > 3 and ancestor in genes[gene][3]:
			finalcodes.append(gene)
	print len(finalcodes)

	outfile =  Fasta + '_'+ ancestor + '.fasta' 
	seqs = grab_seqs(Fasta, finalcodes)
	handle = open( outfile, 'w')
	SeqIO.write(seqs, handle, 'fasta')
	handle.close()

def get_taxinfo(Fasta):
	# ancestor must be an NCBI taxa
	protcodes = grab_IDs(Fasta)
	print protcodes
	for i, codes in enumerate(chunks(protcodes ,99) ):
		print codes
		if i == 0:
			names, Lineages, genes = uniprotIDlist_Lineages(codes)
		else:
			newnames, newlineages, newgenes = uniprotIDlist_Lineages(codes)
			names += newnames
			Lineages += newlineages
			for gene in newgenes:
				genes[gene] = newgenes[gene]
	return names, Lineages, genes

def uniprotIDlist_Lineages(IdList):
	lastset = []
	taxlist = []
	names = []
	Lineages =[]
	genes = {}
	try:
		genes  = uni.map(list(IdList), f='ACC', t='P_ENTREZGENEID')
		# handle = Entrez.efetch(db="nuccore", id=list(lastset), retmode = 'xml')
		if len(genes)>0:
			#values = genes.values()
			#handle = Entrez.elink(dbfrom="nucleotide", id=values, linkname="gene_taxonomy")
			#taxa = Entrez.read(handle, validate = False)
			#for i,entry in enumerate(taxa):
			#	taxlist.append(entry['LinkSetDb'][0]['Link'][0]['Id'])
			for Uniprot in genes:	
				genes[Uniprot] = list(genes[Uniprot])
				#if entry['IdList'][0] in genes[Uniprot]:
				#genes[Uniprot].append(entry['LinkSetDb'][0]['Link'][0]['Id'])
	except: 
		pass
	
	codes = getTax(IdList)
	print len(codes)
	for code in codes.keys():
		if code not in genes:
			genes[code] = [ 'noGeneID' , codes[code] ]
		else :
			genes[code].append(codes[code])	
		taxlist.append(codes[code])
	print len(taxlist)

	taxonomydata = Entrez.read(Entrez.efetch(db="taxonomy", id=taxlist))
	for i,record in enumerate(taxonomydata):
		name = record['ScientificName']
		Lineage = record['Lineage']
		print record
		for gene in genes:
			if len(genes[gene])>1 and genes[gene][1] == taxlist[i]:
				genes[gene].append(name)
				genes[gene].append(Lineage)
		names.append(name)
		Lineages.append(Lineage)

	#annotate the availability of genome and the reference genomes available
	handle = Entrez.elink(dbfrom="taxonomy", id=taxlist, linkname="taxonomy_genome")
	genomes = Entrez.read(handle, validate = False)
	for i,entry in enumerate(genomes):
		for protcode in genes:
			genomefound = False
			if entry['IdList'][0] in genes[protcode]:
				try:
					genes[protcode].append(entry['LinkSetDb'][0]['Link'][0]['Id'])
					genomefound = True
					break
				except:
					print entry
		if genomefound == False:
			genes[protcode].append('noGenome')
	return names, Lineages, genes 

def grabGenomes(taxlist):
	genomedict = {}
	
	
	
	
	taxlist = list(taxlist)

	for i, chunk in enumerate(chunks(taxlist , 400) ):
		handle = Entrez.efetch(db="taxonomy", id=chunk)
		if i == 0:
			taxrecords = Entrez.read(handle, validate = False)
		else :
			taxrecords+=(Entrez.read(handle, validate = False))
		print i 
	for i,taxa in enumerate(taxlist):
		rec = taxrecords[i]
		genomedict[taxa] = []
		genomedict[taxa].append(rec['ScientificName'])
		genomedict[taxa].append(rec['Lineage'])

	for i,  chunk in enumerate(chunks(taxlist ,400) ):
		handle = Entrez.elink(dbfrom="taxonomy", id=chunk, linkname="taxonomy_genome")
		if i == 0:
			genomes = Entrez.read(handle, validate = False)
		else:
			genomes += Entrez.read(handle, validate = False)
		print i 
	for taxa in taxlist:
		genomefound = False
		for i,entry in enumerate(genomes):
			if entry['IdList'][0] in taxa:
				try:
					genomedict[taxa].append(entry['LinkSetDb'][0]['Link'][0]['Id'])
					genomefound = True 
					break
				except:
					print entry
				
		if genomefound == False:
			genomedict[taxa].append('noGenome')
	return genomedict
	
def test_tax():
	names , genes = uniprotIDlist_Lineages(['F4JP36'])
	print names
	print genes

def test_genome():
	genomedict = grabGenomes(['664439' , '944289' , '8010' , '483514' ])
	print genomedict


def getTax_species_name(taxid):
    #fetch taxonomy info for each version of hap2
    Entrez.email = "dmoi@iibintech.com.ar"
    try:
        record = Entrez.read(Entrez.efetch(db="taxonomy", id=taxid))
    except:
        print 'error species name'
        return None
    return record[0]['ScientificName'], record[0]['Lineage']


def getTax(prots):
	#fetch taxonomy info for each version of hap2
	taxdict= {}
	uniprotEntry = uni.retrieve(prots)
	lines = uniprotEntry.split('\n')
	#print lines
	code = ''
	ready  = False
	for line in lines:
		if 'ID' in line[0:2]:
			for prot in prots:
				if prot in line:
					code = prot
		if 'OX' in line[0:2]:
			taxid =line.split(' ')[3].replace('NCBI_TaxID=','').replace(';','')
			ready = True
		if '//' in line and ready == True and code != '':
			taxdict[code] = taxid
	return taxdict
#filter_lineages( '../mergeLineages/hmmerscan/all_hmmercut.fasta' , 'Eukaryota')