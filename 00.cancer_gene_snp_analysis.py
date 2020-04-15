import string, os, sys

cwd= os.getcwd()
pa= cwd.split('/')[-1]

files= os.listdir(cwd)

cgs= {}
for line in open("../Vogelstein_2013_Science_genes.txt").readlines():
	f= line.rstrip().split("\t")
	cgs[f[0]]= f[1]		# ALK= Oncogene, APC= TSG

docm= {}
for line in open("../DoCM_variants_3.2.tsv").readlines():
	f= line.rstrip().split("\t")
	if f[0] == "hgvs": continue
	docm[ f[1]+':'+f[2]+'-'+f[3]+'_'+f[4]+'>'+f[5] ]= f[-2]		#1:1169361-1169361_C>G= renal carcinoma


muts_interest= {}

for file in files:
	if not file.endswith(".tsv"): continue	# Mutect files

	lines= open(file).readlines()
	header= []
	for line in lines:

		tmp_f= line.rstrip().split("\t")
		f= []
		for item in tmp_f:
			f.append(item.strip())

		if f[0] == "type":
			header= f
			continue

		if f[header.index("judgement")] == "REJECT": continue

		mut= f[1]+':'+f[2]+'-'+f[2]+'_'+f[4]+'>'+f[5]

		genes= f[header.index("gene")].replace(';', ',').split(',')
		for gene in genes:
			if not gene in cgs: continue
			funcs= f[header.index("func")].split(';')
			exonicfunc= f[header.index("exonicfunc")]
			aachange= f[header.index("aachange")]
			cadd_phred= f[header.index("cadd_phred")]

			if mut in docm: muts_interest[mut]= [gene, cgs[gene], ";".join(funcs), exonicfunc, aachange, cadd_phred]
			elif "splicing" in funcs:
				if cgs[gene] == "TSG": muts_interest[f[1]+':'+f[2]+'-'+f[2]+'_'+f[4]+'>'+f[5]]= [gene, cgs[gene], ";".join(funcs), exonicfunc, aachange, cadd_phred]
			elif exonicfunc == "stopgain":
				if cgs[gene] == "TSG": muts_interest[f[1]+':'+f[2]+'-'+f[2]+'_'+f[4]+'>'+f[5]]= [gene, cgs[gene], ";".join(funcs), exonicfunc, aachange, cadd_phred]


fo= open(pa+"_cancer_gene_snps.txt", 'w')

for mut in muts_interest:
	if mut in docm:
		fo.write(mut+"\t"+"\t".join(muts_interest[mut])+"\t"+docm[mut]+"\n")
	else:
		fo.write(mut+"\t"+"\t".join(muts_interest[mut])+"\tNA\n")

fo.close()

