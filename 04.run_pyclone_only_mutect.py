import string, os

local_files= os.listdir('.')

cwd = os.getcwd()
pa= cwd.split('/')[-1]

ids= []
input_file_dict= {}
purity_dict= {}

for local_file in local_files:
	if not local_file.endswith("_pyclone_input_only_mutect.txt"): continue
	id= local_file[:local_file.find('_')]
	ids.append(id)
	input_file_dict[id]= local_file

	lines= open(id+"_alternative_solutions.txt").readlines()
	purity_dict[id]= lines[1].rstrip().split("\t")[0]

ids.sort()

fo= open("pyclone_only_mutect.lsf", 'w')
fo.write("#BSUB -J "+pa+"\n")
fo.write("#BSUB -W 240:00\n")
fo.write("#BSUB -o "+os.path.join(cwd, "pyclone_only_mutect_stdout.txt")+"\n")
fo.write("#BSUB -e "+os.path.join(cwd, "pyclone_only_mutect_stderr.txt")+"\n")
fo.write("#BSUB -cwd "+cwd+"\n")
fo.write("#BSUB -q long\n")
fo.write("#BSUB -n 1\n")
fo.write("#BSUB -M 16384\n")
fo.write("#BSUB -R rusage[mem=16384]\n\n")

input_file_thread= []
tumor_content_thread= []
for id in ids:
	input_file_thread.append( input_file_dict[id] )
	tumor_content_thread.append( purity_dict[id] )

fo.write("PyClone run_analysis_pipeline --in_files "+ ' '.join(input_file_thread) +" --working_dir pyclone_analysis_only_mutect_w_adjusted_counts --density pyclone_beta_binomial --num_iters 10000 --burnin 1000 --tumour_contents "+ ' '.join(tumor_content_thread) +" --samples "+ ' '.join(ids) +"\n")

fo.close()

os.system("bsub < pyclone_only_mutect.lsf")
