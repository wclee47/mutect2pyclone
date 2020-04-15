import string, os

database_cut= 0.01

lines= open("all_mutect_catted.tsv").readlines()
fo= open("muts_in_public_database.txt", 'w')

header= []

muts= []

for line in lines:
	tmp= line.rstrip().split("\t")
	f= []
	for col in tmp:
		f.append(col.strip())

	if f[0] == "type":
		header= f
		continue

	mut= f[header.index("chr")] + "\t" + f[header.index("start")] + "\t" + f[header.index("ref_allele")] + "\t" + f[header.index("alt_allele")]

	esp6500siv2_all= f[header.index("esp6500siv2_all")]
	exac_all= f[header.index("exac_all")]
	x1kg2015aug_max= f[header.index("x1kg2015aug_max")]

	will_skip= False
	if esp6500siv2_all != '.':
		if float(esp6500siv2_all) >= database_cut: will_skip= True
	if exac_all != '.':
		if float(exac_all) >= database_cut: will_skip= True
	if x1kg2015aug_max != '.':
		if float(x1kg2015aug_max) >= database_cut: will_skip= True

	if will_skip:
		if not mut in muts:
			muts.append(mut)
			fo.write(mut+"\n")

fo.close()
