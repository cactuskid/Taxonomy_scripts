#utility functions for getting taxonomic data

from Bio import Entrez
import uniprot as uni

def uniprotIDlist_Lineages(IdList):
#fetch a batch of uniprot IDs, some may now be mapped by entrez link. Do not use if you need 100% mapping
	genes  = uni.map(list(IdList), f='ACC', t='P_ENTREZGENEID')
	lastset = []
	taxlist = []
	names = []
	Lineages =[]
	Entrez.email = 'dmoi@iibintech.com.ar'
	# handle = Entrez.efetch(db="nuccore", id=list(lastset), retmode = 'xml')
	if len(genes)>0:
		values = genes.values()
		handle = Entrez.elink(dbfrom="nucleotide", id=values, linkname="gene_taxonomy")
		taxa = Entrez.read(handle, validate = False)
		for i,entry in enumerate(taxa):
			taxlist.append(entry['LinkSetDb'][0]['Link'][0]['Id'])
			for Uniprot in genes:	
				genes[Uniprot] = list(genes[Uniprot])
				if entry['IdList'][0] in genes[Uniprot]:
					genes[Uniprot].append(entry['LinkSetDb'][0]['Link'][0]['Id'])
	missing = []
	
	print len(taxlist)

	for prot in IdList:
		if prot not in genes.keys():
			missing.append(prot)
	missingcodes = getTax(missing)

	for code in missingcodes.keys():
		genes[code] = [ 'noGeneID' , missingcodes[code] ]
		taxlist.append(missingcodes[code])
	
	print len(taxlist)

	taxonomydata = Entrez.read(Entrez.efetch(db="taxonomy", id=taxlist))
	print taxonomydata
	for i,record in enumerate(taxonomydata):
		name = record['ScientificName']
		Lineage = record['Lineage']
		for gene in genes:
			if len(genes[gene])>1 and genes[gene][1] == taxlist[i]:
				genes[gene].append(name)
				genes[gene].append(Lineage)
		names.append(name)
		Lineages.append(Lineage)
	return names, Lineages, genes 

def test_tax():
	names , genes = uniprotIDlist_Lineages(['F4JP36'])
	print names
	print genes


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
	i =0
	for line in lines:
		if 'AC' in line[0:2] and prots[i] in line:
			code = prots[i]
			i += 1
		else:
			for prot in prots:
				if prot in line:
					code = prot
		if 'OX' in line[0:2]:
			taxid =line.split(' ')[3].replace('NCBI_TaxID=','').replace(';','')
		if '//' in line:
			taxdict[code] = taxid
	return taxdict

