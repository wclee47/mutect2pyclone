import string, os, sys

cwd = os.getcwd()
pa= cwd.split('/')[-1]

depth_cut1= 50
depth_cut2= 30
normal_af_cut= 0.01
tumor_af_cut= 0.05


###########################
# public database filtering
###########################

# snp138NonFlagged
dbsnp= {}
for row in open("../ex1.hg19_snp138NonFlagged_dropped").readlines():
	tmp= row.rstrip().split("\t")
	dbsnp[tmp[2]+"\t"+tmp[3]+"\t"+tmp[4]+"\t"+tmp[5]+"\t"+tmp[6]]= ''


# 1000G, etc.
public_db= {}	# mutations >= 1% in public database (mutations to be removed)
for row in open("../muts_in_public_database.txt").readlines():
	tmp= row.rstrip().split("\t")
	public_db[tmp[0]+"\t"+tmp[1]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+tmp[3]]= ''

#####
cancer_gene_muts= {}
for row in open(pa+"_cancer_gene_snps.txt").readlines():
	tmp= row.rstrip().split("\t")
	cancer_gene_muts[tmp[0]]= tmp[-1]	# 3:178936082-178936082_G>A => DoCM


files= os.listdir('.')

union_muts= []

for file in files:
	if not file.endswith(".tsv"): continue
	if file.endswith("neoantigens_terrence.tsv"): continue

	id= file.split('.')[0]

	print(id)

	lines= open(file).readlines()
	header= []
	for line in lines:
		tmp= line.rstrip().split("\t")
		f= []
		for col in tmp:
			f.append(col.strip())
		if f[0] == "type":
			header= f
			continue

		if f[header.index("judgement")] == "REJECT": continue

		mut= f[header.index("chr")] + "\t" + f[header.index("start")] + "\t" + f[header.index("end")] + "\t" + f[header.index("ref_allele")] + "\t" + f[header.index("alt_allele")]

		mut_alt_key= f[header.index("chr")] + ":" + f[header.index("start")] + "-" + f[header.index("end")] + "_" + f[header.index("ref_allele")] + ">" + f[header.index("alt_allele")]
		if mut_alt_key in cancer_gene_muts:
			if cancer_gene_muts[mut_alt_key] != "NA":	# DoCM reported
				if not mut in union_muts: union_muts.append(mut)
				continue

		t_ref_count= int(f[header.index("t_ref_count")])
		t_alt_count= int(f[header.index("t_alt_count")])
		t_dp= t_ref_count + t_alt_count
		if t_dp < depth_cut1: continue
		t_f= t_alt_count / float(t_dp)

		n_ref_count= int(f[header.index("n_ref_count")])
		n_alt_count= int(f[header.index("n_alt_count")])
		n_dp= n_ref_count + n_alt_count
		if n_dp < depth_cut2: continue
		n_f= n_alt_count / float(n_dp)

		if t_f < tumor_af_cut: continue
		if n_f > normal_af_cut: continue

		if mut in dbsnp or mut in public_db: continue

		if not mut in union_muts: union_muts.append(mut)


for file in files:
	if not file.endswith(".tsv"): continue
	if file.endswith("neoantigens_terrence.tsv"): continue

	id= file.split('.')[0]

	print(id)

	fo= open(id+"_mutect_filtered.bed", 'w')

	lines= open(file).readlines()
	header= []
	for line in lines:
		tmp= line.rstrip().split("\t")
		f= []
		for col in tmp:
			f.append(col.strip())
		if f[0] == "type":
			header= f
			continue

		mut= f[header.index("chr")] + "\t" + f[header.index("start")] + "\t" + f[header.index("end")] + "\t" + f[header.index("ref_allele")] + "\t" + f[header.index("alt_allele")]

		if mut in union_muts:
			mut_f= mut.split("\t")
			fo.write(mut_f[0]+"\t"+str(int(mut_f[1])-1)+"\t"+mut_f[2]+"\t"+mut_f[3]+"\t"+mut_f[4]+"\t"+f[header.index("gene")]+"\t"+f[header.index("aachange")]+"\t"+f[header.index("func")]+"\t"+f[header.index("exonicfunc")]+"\n")

	fo.close()

	os.system("bedtools intersect -a %s_mutect_filtered.bed -b /rsrch2/iacs/ngs/captureBeds/E1_sureSelectV4_hg19.bed -wa -u > %s_mutect_filtered_targeted.bed" %(id, id))
	os.system("rm %s_mutect_filtered.bed" %id)

