import string, os

ref_peudo_count= "200"

files= os.listdir('.')

union_muts= []

for file in files:
	if not file.endswith("mutect_filtered_targeted.bed"): continue

	lines= open(file).readlines()

	for line in lines:
		f= line.rstrip().split("\t")
		mut= f[0]+"\t"+str(int(f[1])+1)+"\t"+f[2]+"\t"+f[3]+"\t"+f[4]
		if not mut in union_muts: union_muts.append(mut)

for file in files:
	if not file.endswith(".tsv"): continue

	id= file.split('.')[0]	# should customized

	dict= {}

	lines= open(file).readlines()
	header= []
	for line in lines:
		tmp= line.rstrip().split("\t")
		f= []
		for col in tmp:
			f.append(col.strip())
		if f[0] == "type":	# should customized
			header= f
			continue

		mut= f[header.index("chr")]+"\t"+f[header.index("start")]+"\t"+f[header.index("end")]+"\t"+f[header.index("ref_allele")]+"\t"+f[header.index("alt_allele")]

		ref_allele_lowered= f[header.index("ref_allele")].lower()
		alt_allele_lowered= f[header.index("alt_allele")].lower()
		t_ref_count, t_alt_count= 0, 0

		if f[header.index(ref_allele_lowered+"_sample")] == "NA": t_ref_count= 0
		else: t_ref_count= int(f[header.index(ref_allele_lowered+"_sample")])
		if f[header.index(alt_allele_lowered+"_sample")] == "NA": t_alt_count= 0
		else: t_alt_count= int(f[header.index(alt_allele_lowered+"_sample")])

		dict[mut]= str(t_ref_count)+"\t"+str(t_alt_count)

	fo= open(id+"_pyclone_snp_intermediate.txt", 'w')	# should customized
	fo.write("mutation_id\tref_counts\tvar_counts\n")

	for mut in union_muts:
		mut_f= mut.split("\t")
		mut_id= mut_f[0]+':'+mut_f[1]+'-'+mut_f[2]+'_'+mut_f[3]+'>'+mut_f[4]
		if mut in dict:
			fo.write(mut_id+"\t"+dict[mut]+"\n")
		else:
			fo.write(mut_id+"\t"+ref_peudo_count+"\t0\n")

	fo.close()


