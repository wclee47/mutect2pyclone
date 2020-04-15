import string, os

local_files= os.listdir('.')

dict= {}
ids= []
mut_ids= []
mut_ids_to_be_removed= []

for local_file in local_files:
	if not local_file.endswith("_pyclone_snp_intermediate.txt"): continue

	lines= open(local_file).readlines()

	id= local_file.split("_pyclone")[0]
	ids.append(id)

	cn_file= id+"_segments.txt"

	header= []
	for line in lines:
		f= line.rstrip().split("\t")
		if f[0] == "mutation_id":
			header= f
			continue

		mut_id= f[0]
		if not mut_id in mut_ids: mut_ids.append(mut_id)

		mut_chr= mut_id.split(':')[0]
		mut_start= int( mut_id.split(':')[1].split('_')[0].split('-')[0] )
		ref_count= f[1]
		var_count= f[2]
		normal_cn= '2'
		minor_cn= "na"
		major_cn= "na"
		for cn_line in open(cn_file).readlines():
			cn_f= cn_line.rstrip().split("\t")
			cn_chr= cn_f[0][1:-1]
			if cn_chr == "chromosome": continue
			cn_start= int(cn_f[1])
			cn_end= int(cn_f[2])
			if mut_chr != cn_chr: continue
			if cn_start < mut_start < cn_end:
				if cn_f[10] == "NA" or cn_f[11] == "NA": break
				cnA= int(cn_f[10])
				cnB= int(cn_f[11])
				minor_cn= cnA
				major_cn= cnB
				if minor_cn > major_cn:
					minor_cn= cnB
					major_cn= cnA
				break

		dict[id+':'+mut_id]= [mut_id, ref_count, var_count, normal_cn, str(minor_cn), str(major_cn)]

		if str(major_cn) == '0': mut_ids_to_be_removed.append(mut_id)   # newly added because of error

		if minor_cn == "na" or major_cn == "na": mut_ids_to_be_removed.append(mut_id)

for id in ids:
	fo= open(id+"_pyclone_input_only_mutect.txt", 'w')
	fo.write("mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\n")
	for mut_id in mut_ids:
		if mut_id in mut_ids_to_be_removed: continue

		fo.write("\t".join(dict[id+':'+mut_id])+"\n")

	fo.close()

