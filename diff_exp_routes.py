from flask import Blueprint, render_template, abort, request
import sqlite3
from sqlitedict import SqliteDict
import ast
import os
import time
import numpy as np
from bisect import bisect_left 
from random import shuffle
from math import log
import config
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code
import riboflask_diff
import collections
from flask_login import current_user

# Differential expression page, used to find diffentially expressed genes via z-score
diff_plotpage_blueprint = Blueprint("diffpage", __name__, template_folder="templates")
@diff_plotpage_blueprint.route('/<organism>/<transcriptome>/differential/')
def diffpage(organism,transcriptome):
	global user_short_passed
	user_short_passed = True
	global local
	try:
		print local
	except:
		local = False
	try:
		user = current_user.name
	except:
		user = None
	organism = str(organism)
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}';".format(organism))
	result = (cursor.fetchone())
	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]
	studyinfo_dict = fetch_study_info(organism)
	# holds all values the user could possibly pass in the url (keywords are after request.args.get), anything not passed by user will be a string: "None"
	html_args = {"user_short":str(request.args.get('short')),
				 "minreads":str(request.args.get('minread')),
				 "minzscore":str(request.args.get('minzscore')),
				 "region":str(request.args.get('region')),
				 "ambig":str(request.args.get('ambig')),
				 "plottype":str(request.args.get('plottype')),
				 "gene_list":str(request.args.get('gene_list')),
				 "transcriptome":str(transcriptome),
				 "min_cov":str(request.args.get('min_cov'))}
	html_args["riboseq_files_1"] = "None"
	html_args["riboseq_files_2"] = "None"
	html_args["rnaseq_files_1"] = "None"
	html_args["rnaseq_files_2"] = "None"
	
	user_files = request.args.get('riboseq_files_1')
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["riboseq_files_1"] = [str(x) for x in user_files]

	user_files = request.args.get('riboseq_files_2')
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["riboseq_files_2"] = [str(x) for x in user_files]

	user_files = request.args.get('rnaseq_files_1')
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["rnaseq_files_1"] = [str(x) for x in user_files]

	user_files = request.args.get('rnaseq_files_2')
	if user_files != None:
		if len(user_files) != 0:
			user_files = user_files.split(",")
			html_args["rnaseq_files_2"] = [str(x) for x in user_files]

	user_labels = request.args.get('riboseq_labels_1')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["riboseq_labels_1"] = [str(x) for x in user_labels]
	else:
		html_args["riboseq_labels_1"] = []

	user_labels = request.args.get('riboseq_labels_2')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["riboseq_labels_2"] = [str(x) for x in user_labels]
	else:
		html_args["riboseq_labels_2"] = []

	user_labels = request.args.get('rnaseq_labels_1')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["rnaseq_labels_1"] = [str(x) for x in user_labels]
	else:
		html_args["rnaseq_labels_1"] = []

	user_labels = request.args.get('rnaseq_labels_2')
	if user_labels != None:
		user_labels = user_labels.split(",")
		html_args["rnaseq_labels_2"] = [str(x) for x in user_labels]
	else:
		html_args["rnaseq_labels_2"] = []
		
	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict,accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)
	connection.close()
	return render_template('index_diff.html', studies_dict=accepted_studies, accepted_files=accepted_files,organism=organism, default_tran=default_tran,local=local,transcriptome=transcriptome,html_args=html_args,studyinfo_dict=studyinfo_dict,seq_types=seq_types)

# Creates/serves the z-score plot for differential expression
diffquery_blueprint = Blueprint("diffquery", __name__, template_folder="templates")
@diffquery_blueprint.route('/diffquery', methods=['POST'])
def diffquery():
	global user_short_passed
	data = ast.literal_eval(request.data)
	plottype = data["plottype"]
	html_args = data["html_args"]
	organism = data["organism"]
	transcriptome = data["transcriptome"]
	store_de_results = data["store_de_results"]
	cond_desc = data["cond_desc"]
	gene_list = data["gene_list"]
	if data["min_cov"] != "undefined":
		min_cov = float(data["min_cov"])
	else:
		min_cov = 0
	if min_cov > 1:
		return "Minimum coverage should be a value between 0 and 1"

	filename = organism+"_differential_translation_"+str(time.time())+".csv"
	csv_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
	master_file_dict = data["master_file_dict"]
	master_transcript_dict = {}
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]
		
	if owner == 1:
		transhelve = sqlite3.connect("{0}{1}/{1}.v2.sqlite".format(config.ANNOTATION_DIR,organism))
	else:
		transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.v2.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
	trancursor = transhelve.cursor()
	trancursor.execute("SELECT * from transcripts WHERE principal = 1 AND tran_type = 'coding'")
	transcript_list =  trancursor.fetchall()
	traninfo_dict = {}
	for result in transcript_list:
		if result[1] == "DDX3Y":
			print "RESULT", result
		traninfo_dict[result[0]] = {"transcript":result[0] , "gene":result[1], "length":result[2] , "cds_start":result[3] , "cds_stop":result[4] , "seq":result[5] ,
				"strand":result[6], "stop_list":result[7].split(","),"start_list":result[8].split(","), "exon_junctions":result[9].split(","),
				"tran_type":result[10], "principal":result[11]}
	transhelve.close()
	longest_tran_list = traninfo_dict.keys()
	
	#This will hold one or more z-scores for every transcript
	master_dict = {}
	
	#Exclusively used for the ribo_vs_rna plot, keys are genes, values are ribo fc and rna fc
	ribo_vs_rna_dict = {}
	#List of tuples, (gene, z_Score) used to find the top 500 most DE genes
	absolute_zcores = []

	riboseq1_filepaths = {}
	file_paths_dict = fetch_file_paths(master_file_dict["riboseq1"]["file_ids"], organism)
	for seq_type in file_paths_dict:
		for file_id in file_paths_dict[seq_type]:
			riboseq1_filepaths[file_id] = file_paths_dict[seq_type][file_id]
	riboseq2_filepaths = {}
	file_paths_dict = fetch_file_paths(master_file_dict["riboseq2"]["file_ids"], organism)
	for seq_type in file_paths_dict:
		for file_id in file_paths_dict[seq_type]:
			riboseq2_filepaths[file_id] = file_paths_dict[seq_type][file_id]
	try:
		rnaseq1_filepaths = {}
		file_paths_dict = fetch_file_paths(master_file_dict["rnaseq1"]["file_ids"], organism)
		for seq_type in file_paths_dict:
			for file_id in file_paths_dict[seq_type]:
				rnaseq1_filepaths[file_id] = file_paths_dict[seq_type][file_id]
		rnaseq2_filepaths = {}
		file_paths_dict = fetch_file_paths(master_file_dict["rnaseq2"]["file_ids"], organism)
		for seq_type in file_paths_dict:
			for file_id in file_paths_dict[seq_type]:
				rnaseq2_filepaths[file_id] = file_paths_dict[seq_type][file_id]
	except:
		rnaseq1_filepaths = {}
		rnaseq2_filepaths = {}
		
	minreads = float(data["minreads"])
	minzscore = float(data["minzscore"])
	
	if "mapped_reads_norm" in data:
		mapped_reads_norm = True
	else:
		mapped_reads_norm = False
		
	if "ambiguous" in data:
		ambiguous = True
	else:
		ambiguous = False
	# User can decide to look at just riboseq fold-change, rnaseq fold-change or TE fold-change
	#if len(master_file_dict["riboseq1"]["file_ids"]) != len(master_file_dict["riboseq2"]["file_ids"]):
	#	return "Error: Both Riboseq Condition boxes need to have an equal number of files."
	#if len(master_file_dict["rnaseq1"]["file_ids"]) != len(master_file_dict["rnaseq2"]["file_ids"]):
	#	return "Error: Both RNA-Seq Condition boxes need to have an equal number of files."
	if len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and (len(master_file_dict["riboseq1"]["file_ids"]) != len(master_file_dict["rnaseq1"]["file_ids"])):
		return "Error: If RNA-Seq boxes are not empty they need to have an equal number of files as the riboseq boxes"
	
	if len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq2"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq2"]["file_ids"]) >= 1:
		label = "TE"
	elif len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq2"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq1"]["file_ids"]) == 0 and len(master_file_dict["rnaseq2"]["file_ids"]) == 0:
		label = "Riboseq"
	elif len(master_file_dict["riboseq1"]["file_ids"]) == 0 and len(master_file_dict["riboseq2"]["file_ids"]) == 0 and len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq2"]["file_ids"]) >= 1:
		label = "Rnaseq"
	else:
		label = "TE"
	#	return "ERROR IMBALANCED OR NO FILES: Either all 4 boxes (RIBO-seq files 1, RIBO-seq files 2, mRNA-seq files 1, mRNA-seq files 2) must have a file associated with it OR both riboseq boxes OR both rnaseq boxes, you currently have {} files in RIBO-seq condition 1 files, {} in RIBO-seq condition 2 files, {} in mRNA-seq condition 1 files and {} in mRNA-seq condition 2 files".format(len(riboseq1_filepaths),len(riboseq2_filepaths),len(rnaseq1_filepaths),len(rnaseq2_filepaths))
	
	no_groups = max(len(master_file_dict["riboseq1"]["file_ids"]),len(master_file_dict["rnaseq1"]["file_ids"]))
	region = data["region"] #can be cds,fiveprime, or threeprime

	riboseq1_tot_reads = 0.001
	riboseq2_tot_reads = 0.001
	rnaseq1_tot_reads = 0.001
	rnaseq2_tot_reads = 0.001
	csv_file.write("Transcript,Gene,")
	# Group together files according to the order they were given in, if a mismatch in group number raise an error. 
	for i in range(0,no_groups):
		if label == "TE":
			csv_file.write("Group,Fold change (log2),Geometric Mean (log2), Bin mean, Bin standard deviation, Z-score,")
			csv_file.write("RiboSeq_Cond1_count,RiboSeq_Cond2_count,mRNA-Seq_Cond1_count,mRNA-Seq_Cond2_count,")
			if mapped_reads_norm:
				csv_file.write("Normalised_RiboSeq_Cond1_count,Normalised_RiboSeq_Cond2_count,Normalised_mRNA-Seq_Cond1_count,Normalised_mRNA-Seq_Cond2_count,")
		elif label == "Riboseq":
			csv_file.write("Group,Fold change (log2),Geometric Mean (log2), Bin mean, Bin standard deviation, Z-score,")
			csv_file.write("RiboSeq_Cond1_count,RiboSeq_Cond2_count,")
			if mapped_reads_norm:
				csv_file.write("Normalised_RiboSeq_Cond1_count ,Normalised_RiboSeq_Cond2_count,")
		else:
			csv_file.write("Group,Fold change (log2),Geometric Mean (log2), Bin mean, Bin standard deviation, Z-score,")
			csv_file.write("mRNA-Seq_Cond1_count,mRNA-Seq_Cond2_count,")
			if mapped_reads_norm:
				csv_file.write("Normalised_mRNA-Seq_Cond1_count ,Normalised_mRNA-Seq_Cond2_count,")
		if plottype == "z_score":
			if len(riboseq1_filepaths) != 0:
				riboseq1_filepath = riboseq1_filepaths[int(master_file_dict["riboseq1"]["file_ids"][i])]
				riboseq2_filepath = riboseq2_filepaths[int(master_file_dict["riboseq2"]["file_ids"][i])]
			else:
				riboseq1_filepath = ""
				riboseq2_filepath = ""
				
			if len(rnaseq1_filepaths) != 0:
				rnaseq1_filepath = rnaseq1_filepaths[int(master_file_dict["rnaseq1"]["file_ids"][i])]
				rnaseq2_filepath = rnaseq2_filepaths[int(master_file_dict["rnaseq2"]["file_ids"][i])]
			else:
				rnaseq1_filepath = ""
				rnaseq2_filepath = ""

			transcript_dict,groupname = calculate_zscore(riboseq1_filepath, riboseq2_filepath, rnaseq1_filepath, rnaseq2_filepath, master_dict, longest_tran_list, mapped_reads_norm,label,region,traninfo_dict,minreads, minzscore,ambiguous,min_cov)
			if transcript_dict == "error":
				return groupname
			master_transcript_dict[groupname] = transcript_dict

	del_list = []
	for tran in master_dict:
		#If the transcript is missing from one of the groups or if the coverage of the transcript is less than the minimum allowed, discard it here
		if (len(master_dict[tran].keys())-1) != no_groups:
			del_list.append(tran)
		elif min_cov > 0:
			max_group_cov = 0
			for group in master_transcript_dict:
				try:
					if master_transcript_dict[group][tran]["max_cov"] > max_group_cov:
						max_group_cov = master_transcript_dict[group][tran]["max_cov"]
				except Exception as e:
					print e
					pass
			if max_group_cov < min_cov:
				del_list.append(tran)
			
	for tran in del_list:
		del master_dict[tran]

	csv_file.write("Average fold change (log2),Average Z-score")
	if no_groups == 2:
		csv_file.write(",FDR")
	csv_file.write("\n")
	aggregated_values = []
	total_cases = 0

	#If there are 2 replicates in each group, emprically calculate the fdr for each gene
	if no_groups == 2:
		all_genuine_z_scores = []
		false_positives_z_dict = {}
		for i in range(0,no_groups):
			false_positives_z_dict[i] = []
		for tran in master_dict:
			geo_mean_list = []
			fc_list = []
			z_scores = []
			for groupname in master_dict[tran]:
				if groupname == "skip":
					continue
				z_score = (master_dict[tran][groupname]["fold_change"]-master_dict[tran][groupname]["mean"])/master_dict[tran][groupname]["standard_dev"]
				z_scores.append(z_score)
			average_z_score = sum(z_scores)/len(z_scores)
			master_dict[tran]["z_score"] = average_z_score
			all_genuine_z_scores.append(average_z_score)
			#If z_scores not in same direction for all groups, this counts as a false positive. 
			if not all(item >= 0 for item in z_scores) and not all(item <= 0 for item in z_scores):
				for i in range(0,len(z_scores)):
					false_positives_z_dict[i].append(abs(float(z_scores[i])))
		p_val_dict = {}
		all_shuffled_z_scores = []
		for i in range(0,10):
			for x in range(0,no_groups):
				shuffle(false_positives_z_dict[x])
			for y in range(0,len(false_positives_z_dict[0])):
				shuffled_z_scores = []
				for t in range(0,no_groups):
					shuffled_z_scores.append(false_positives_z_dict[t][y])
				average_z_score = sum(shuffled_z_scores)/len(shuffled_z_scores)
				all_shuffled_z_scores.append(average_z_score)
		all_shuffled_z_scores.sort()
		genuine_te_set = list(set(all_genuine_z_scores))
		temp_count = 0
		p_val_dict[0] = len(all_shuffled_z_scores)
		if all_shuffled_z_scores != []:
			for gen_te in genuine_te_set:
				if gen_te  > all_shuffled_z_scores[-1]:
					p_val_dict[gen_te] = 0
					continue
				success1s = 0
				success2s = 0
				pos = bisect_left(all_shuffled_z_scores, abs(gen_te))
				reversed_pos = len(all_shuffled_z_scores)-pos
				p_val_dict[abs(gen_te)] = reversed_pos
			for tran in master_dict:
				if "z_score" in master_dict[tran]:
					z_score = master_dict[tran]["z_score"]
					vala = float(len([1 for bb2 in all_genuine_z_scores if bb2 >= abs(z_score)]))
					try:
						valb = float(p_val_dict[abs(z_score)])/10
					except:
						fdr = 1
					try:
						fdr = round(valb/(vala),2)
					except:
						fdr = 0.0
					if fdr > 1:
						fdr = 1
					master_dict[tran]["fdr"] = fdr
	for tran in master_dict:
		total_cases += 1
		geo_mean_list = []
		fc_list = []
		z_scores = []
		gene = (traninfo_dict[tran]["gene"]).replace(",","_").replace(";","_")
		csv_file.write("{},{},".format(tran, gene))
		for groupname in master_dict[tran]:
			if groupname in ["skip","fdr","z_score"]:
				continue
			transcript_dict = master_transcript_dict[groupname]
			geo_mean_list.append(master_dict[tran][groupname]["geometric_mean"])
			fc_list.append(master_dict[tran][groupname]["fold_change"])
			fc = master_dict[tran][groupname]["fold_change"]
			z_score = (fc-master_dict[tran][groupname]["mean"])/master_dict[tran][groupname]["standard_dev"]
			z_scores.append(z_score)
			geo_mean = master_dict[tran][groupname]["geometric_mean"]
			if label == "TE":
				csv_file.write("{},{},{},{},{},{},".format(groupname,fc,geo_mean,master_dict[tran][groupname]["mean"],master_dict[tran][groupname]["standard_dev"],z_score))
				riboseq1_count = (2**transcript_dict[tran]["riboseq1"])*transcript_dict["ribo1_modifier"]-0.0001
				riboseq2_count = (2**transcript_dict[tran]["riboseq2"])*transcript_dict["ribo2_modifier"]-0.0001
				rnaseq1_count = (2**transcript_dict[tran]["rnaseq1"])*transcript_dict["rna1_modifier"]-0.0001
				rnaseq2_count = (2**transcript_dict[tran]["rnaseq2"])*transcript_dict["rna2_modifier"]-0.0001
				if riboseq1_count < 0.1:
					riboseq1_count = 0
				if riboseq2_count < 0.1:
					riboseq2_count = 0
				if rnaseq1_count < 0.1:
					rnaseq1_count = 0
				if rnaseq2_count < 0.1:
					rnaseq2_count = 0
				csv_file.write("{},{},{},{},".format(riboseq1_count,riboseq2_count,rnaseq1_count,rnaseq2_count))
				if gene not in ribo_vs_rna_dict:
					ribo_vs_rna_dict[gene] = {"tran":tran,
											"ribo1":2**transcript_dict[tran]["riboseq1"],
											"ribo2":2**transcript_dict[tran]["riboseq2"],
											"rna1":2**transcript_dict[tran]["rnaseq1"],
											"rna2":2**transcript_dict[tran]["rnaseq2"]}
				else:
					ribo_vs_rna_dict[gene]["ribo1"] += 2**transcript_dict[tran]["riboseq1"]
					ribo_vs_rna_dict[gene]["ribo2"] += 2**transcript_dict[tran]["riboseq2"]
					ribo_vs_rna_dict[gene]["rna1"] += 2**transcript_dict[tran]["rnaseq1"]
					ribo_vs_rna_dict[gene]["rna2"] += 2**transcript_dict[tran]["rnaseq2"]
				if mapped_reads_norm:
					csv_file.write("{},{},{},{},".format((2**transcript_dict[tran]["riboseq1"])-0.0001,(2**transcript_dict[tran]["riboseq2"])-0.0001,(2**transcript_dict[tran]["rnaseq1"])-0.0001,(2**transcript_dict[tran]["rnaseq2"])-0.0001))
					
			elif label == "Riboseq":
				csv_file.write("{},{},{},{},{},{},".format(groupname,fc,geo_mean,master_dict[tran][groupname]["mean"],master_dict[tran][groupname]["standard_dev"],z_score))
				riboseq1_count = (2**transcript_dict[tran]["riboseq1"])*transcript_dict["ribo1_modifier"]-0.0001
				riboseq2_count = (2**transcript_dict[tran]["riboseq2"])*transcript_dict["ribo2_modifier"]-0.0001
				if riboseq1_count < 0.1:
					riboseq1_count = 0
				if riboseq2_count < 0.1:
					riboseq2_count = 0
				csv_file.write("{},{},".format(riboseq1_count, riboseq2_count))
				if gene not in ribo_vs_rna_dict:
					ribo_vs_rna_dict[gene] = {"tran":tran,
											"ribo1":2**transcript_dict[tran]["riboseq1"],
											"ribo2":2**transcript_dict[tran]["riboseq2"]}
				else:
					ribo_vs_rna_dict[gene]["ribo1"] += 2**transcript_dict[tran]["riboseq1"]
					ribo_vs_rna_dict[gene]["ribo2"] += 2**transcript_dict[tran]["riboseq2"]
				if mapped_reads_norm:
					csv_file.write("{},{},".format((2**transcript_dict[tran]["riboseq1"])-0.0001,(2**transcript_dict[tran]["riboseq2"])-0.0001))
			elif label == "Rnaseq":
				csv_file.write("{},{},{},{},{},{},".format(groupname,fc,geo_mean,master_dict[tran][groupname]["mean"],master_dict[tran][groupname]["standard_dev"],z_score))
				rnaseq1_count = (2**transcript_dict[tran]["rnaseq1"])*transcript_dict["rna1_modifier"]-0.0001
				rnaseq2_count = (2**transcript_dict[tran]["rnaseq2"])*transcript_dict["rna2_modifier"]-0.0001
				if rnaseq1_count < 0.1:
					rnaseq1_count = 0
				if rnaseq2_count < 0.1:
					rnaseq2_count = 0
				csv_file.write("{},{},".format(rnaseq1_count, rnaseq2_count))
				if gene not in ribo_vs_rna_dict:
					ribo_vs_rna_dict[gene] = {"tran":tran,
											"rna1":2**transcript_dict[tran]["rnaseq1"],
											"rna2":2**transcript_dict[tran]["rnaseq2"]}
				else:
					ribo_vs_rna_dict[gene]["rna1"] += 2**transcript_dict[tran]["rnaseq1"]
					ribo_vs_rna_dict[gene]["rna2"] += 2**transcript_dict[tran]["rnaseq2"]
				if mapped_reads_norm:
					csv_file.write("{},{},".format((2**transcript_dict[tran]["rnaseq1"])-0.0001,(2**transcript_dict[tran]["rnaseq2"])-0.0001))	
		try:
			average_geo_mean = (sum(geo_mean_list)/len(geo_mean_list))
		except:
			average_geo_mean = 0
		try:
			average_fc_mean = (sum(fc_list)/len(fc_list))
		except:
			average_fc_mean = 0
		average_z_score = sum(z_scores)/len(z_scores)
		master_dict[tran]["z_score"] = average_z_score
		absolute_zcores.append((str(gene),abs(average_z_score),average_z_score))
		
		fc = average_fc_mean
		aggregated_values.append([tran,average_geo_mean,fc,gene])
		csv_file.write("{},{}".format(fc, average_z_score))
		if no_groups == 2:
			if "fdr" in master_dict[tran]:
				csv_file.write(",{}".format(master_dict[tran]["fdr"]))
		csv_file.write("\n")
	if store_de_results:
		sorted_z_scores = sorted(absolute_zcores,key=lambda x: x[1])
		top500 =  sorted_z_scores[-500:]
		cursor.execute("SELECT organism_id from organisms WHERE organism_name = '{}';".format(organism))
		result = cursor.fetchone()
		organism_id = result[0]
		cursor.execute("SELECT study_id from files WHERE file_id = {};".format(file_id))
		result = cursor.fetchone()
		study_id = result[0]
		
		#Insert the condition description to the gene_reg table (don't pass a condition id, will auto increment)
		cursor.execute("INSERT INTO gene_regulation VALUES(null,'{}',{},{});".format(cond_desc,study_id, organism_id))
		#Grab the condition id that this description was just assigned
		cursor.execute("SELECT MAX(condition_id) from gene_regulation;")
		result = cursor.fetchone()
		condition_id = result[0]
		for tup in top500:
			gene = tup[0]
			z_score = tup[2]
			cursor.execute("INSERT INTO gene_reg_z_scores VALUES('{}',{},{})".format(gene,condition_id,z_score))
	
	connection.commit()
	sorted_aggregated_values = sorted(aggregated_values, key=lambda x: x[1])
	bin_list = []
	for i in range(0, len(sorted_aggregated_values), 300):
		# if we are not near the end of the list calculate mean and std dev
		if i < (len(sorted_aggregated_values)-300):
			#for every transcript in this bin add the min exp to bin_values list
			bin_values = []
			for x in range(i,i+300):
				bin_values.append(sorted_aggregated_values[x][2])
			#work out the mean and standard deviation of bin_values
			mean = np.mean(bin_values)
			standard_dev = np.std(bin_values)
			# Append mean and std dev to sorted_aggregated_values so we can work out z-score for each gene later
			for x in range(i,i+300):
				sorted_aggregated_values[x].append(mean)
				sorted_aggregated_values[x].append(standard_dev)
			y = minzscore*(standard_dev)
			upper_threshold = y+mean
			lower_threshold = (-1*y)+mean
			bin_list.append([mean, standard_dev, upper_threshold,lower_threshold])
		else:
			bin_list.append(bin_list[-1])
			# Append mean and std dev to sorted_aggregated_values so we can work out z-score for each gene later
			for x in range(i,len(sorted_aggregated_values)):
				sorted_aggregated_values[x].append(bin_list[-1][0])
				sorted_aggregated_values[x].append(bin_list[-1][1])
	try:
		user = current_user.name
	except:
		user = None
	connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
	connection.text_factory = str
	cursor = connection.cursor()
	background_col = config.BACKGROUND_COL
	uga_col = config.UGA_COL
	uag_col = config.UAG_COL
	uaa_col = config.UAA_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE
	if user != None:
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		#get a list of organism id's this user can access
		cursor.execute("SELECT background_col,title_size,subheading_size,axis_label_size,marker_size from user_settings WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchone())
		background_col = result[0]
		title_size = result[1]
		subheading_size = result[2]
		axis_label_size = result[3]
		marker_size = result[4]
		connection.close()
	if html_args["user_short"] == "None" or user_short_passed == True:
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"differential")
	else:
		short_code = html_args["user_short"]
		user_short_passed = True
	if plottype == "z_score":
		return riboflask_diff.generate_plot(sorted_aggregated_values,
										bin_list,
										organism,
										label,
										transcriptome,
										master_file_dict["riboseq1"]["file_ids"],
										master_file_dict["riboseq2"]["file_ids"],
										master_file_dict["rnaseq1"]["file_ids"],
										master_file_dict["rnaseq2"]["file_ids"],
										background_col,
										short_code,
										mapped_reads_norm,
										filename,
										no_groups,
										str(title_size)+"pt",
										str(axis_label_size)+"pt", 
										str(subheading_size)+"pt",
										str(marker_size)+"pt",
										ambiguous,
										gene_list)
	else:
		return riboflask_diff.ribo_vs_rna(ribo_vs_rna_dict,
										organism,
										transcriptome,
										master_file_dict["riboseq1"]["file_ids"],
										master_file_dict["riboseq2"]["file_ids"],
										master_file_dict["rnaseq1"]["file_ids"],
										master_file_dict["rnaseq2"]["file_ids"],
										background_col,
										short_code,
										mapped_reads_norm,
										filename,
										no_groups,
										str(title_size)+"pt",
										str(axis_label_size)+"pt", 
										str(subheading_size)+"pt",
										str(marker_size)+"pt",
										ambiguous,
										gene_list,
										label)


# Given either two or four filepaths, calculates a z-score, places the z-scores in a master dict 
def calculate_zscore(riboseq1_filepath, riboseq2_filepath, rnaseq1_filepath, rnaseq2_filepath,master_dict, longest_tran_list, mapped_reads_norm,label,region,traninfo_dict, minreads, minzscore,ambiguous,min_cov):
	print "CALCULATE Z SCORE CALLED"
	#if minreads != 0:
	#	minreads = log(minreads,2)
	#else:
	#	minreads = -1000
	riboseq1_tot_reads = 0.001
	riboseq2_tot_reads = 0.001
	rnaseq1_tot_reads = 0.001
	rnaseq2_tot_reads = 0.001
	transcript_dict ={}
	if ambiguous == False:
		ambig_type = "unambiguous"
	elif ambiguous == True:
		ambig_type = "ambiguous"
		
		
	groupname = ""
	covdict = {}
	if riboseq1_filepath:
		groupname += (riboseq1_filepath.split("/")[-1]).replace(".sqlite","")+"_"
		if os.path.isfile(riboseq1_filepath):
			sqlite_db = SqliteDict(riboseq1_filepath, autocommit=False)
		else:
			return "error","File not found: {}".format(groupname)
		if region == "fiveprime":
			opendict =sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict =sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict =sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict =sqlite_db["{}_all_totals".format(ambig_type)]

		if mapped_reads_norm == True:
			riboseq1_tot_reads += float(sqlite_db["noncoding_counts"])
			riboseq1_tot_reads += float(sqlite_db["coding_counts"])
		if min_cov > 0:
			if "{}_all_coverage".format(ambig_type) not in sqlite_db:
				calculate_coverages(sqlite_db,longest_tran_list,ambig_type, region,traninfo_dict)
			if region == "fiveprime":
				covdict =sqlite_db["{}_fiveprime_coverage".format(ambig_type)]
			elif region == "cds":
				covdict =sqlite_db["{}_cds_coverage".format(ambig_type)]
			elif region == "threeprime":
				covdict =sqlite_db["{}_threeprime_coverage".format(ambig_type)]
			elif region == "all":
				covdict =sqlite_db["{}_all_coverage".format(ambig_type)]
		sqlite_db.close()
		for transcript in longest_tran_list:
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001,
												"riboseq1_cov":0.0001,"riboseq2_cov":0.0001, "rnaseq1_cov":0.0001, "rnaseq2_cov":0.0001,
												"max_cov":0}
			if transcript in opendict:
				transcript_dict[transcript]["riboseq1"] += opendict[transcript]
			if transcript in covdict:
				transcript_dict[transcript]["riboseq1_cov"] += covdict[transcript]
				
	if riboseq2_filepath:
		groupname += (riboseq2_filepath.split("/")[-1]).replace(".sqlite","")+"_"
		if os.path.isfile(riboseq2_filepath):
			sqlite_db = SqliteDict(riboseq2_filepath, autocommit=False)
		else:
			return "error","File not found, please report this to tripsvizsite@gmail.com or via the contact page."
		if region == "fiveprime":
			opendict =sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict =sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict =sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict =sqlite_db["{}_all_totals".format(ambig_type)]
		if mapped_reads_norm == True:
			riboseq2_tot_reads += float(sqlite_db["noncoding_counts"])
			riboseq2_tot_reads += float(sqlite_db["coding_counts"])
		if min_cov > 0:
			if "{}_all_coverage".format(ambig_type) not in sqlite_db:
				calculate_coverages(sqlite_db,longest_tran_list,ambig_type, region,traninfo_dict)
			if region == "fiveprime":
				covdict =sqlite_db["{}_fiveprime_coverage".format(ambig_type)]
			elif region == "cds":
				covdict =sqlite_db["{}_cds_coverage".format(ambig_type)]
			elif region == "threeprime":
				covdict =sqlite_db["{}_threeprime_coverage".format(ambig_type)]
			elif region == "all":
				covdict =sqlite_db["{}_all_coverage".format(ambig_type)]
		sqlite_db.close()
		for transcript in longest_tran_list:
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001,
												"riboseq1_cov":0.0001,"riboseq2_cov":0.0001, "rnaseq1_cov":0.0001, "rnaseq2_cov":0.0001,
												"max_cov":0}
			if transcript in opendict:
				transcript_dict[transcript]["riboseq2"] += opendict[transcript]
			if transcript in covdict:
				transcript_dict[transcript]["riboseq2_cov"] += covdict[transcript]

	if rnaseq1_filepath:
		groupname += (rnaseq1_filepath.split("/")[-1]).replace(".sqlite","")+"_"
		if os.path.isfile(rnaseq1_filepath):
			sqlite_db = SqliteDict(rnaseq1_filepath, autocommit=False)
		else:
			return "error","File not found, please report this to tripsvizsite@gmail.com or via the contact page."
		if region == "fiveprime":
			opendict =sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict =sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict =sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict =sqlite_db["{}_all_totals".format(ambig_type)]
		if mapped_reads_norm == True:
			rnaseq1_tot_reads += float(sqlite_db["noncoding_counts"])
			rnaseq1_tot_reads += float(sqlite_db["coding_counts"])
		if min_cov > 0:
			if "{}_all_coverage".format(ambig_type) not in sqlite_db:
				calculate_coverages(sqlite_db,longest_tran_list,ambig_type, region,traninfo_dict)
			if region == "fiveprime":
				covdict =sqlite_db["{}_fiveprime_coverage".format(ambig_type)]
			elif region == "cds":
				covdict =sqlite_db["{}_cds_coverage".format(ambig_type)]
			elif region == "threeprime":
				covdict =sqlite_db["{}_threeprime_coverage".format(ambig_type)]
			elif region == "all":
				covdict =sqlite_db["{}_all_coverage".format(ambig_type)]
		sqlite_db.close()
		for transcript in longest_tran_list:
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001,
												"riboseq1_cov":0.0001,"riboseq2_cov":0.0001, "rnaseq1_cov":0.0001, "rnaseq2_cov":0.0001,
												"max_cov":0}
			if transcript in opendict:
				transcript_dict[transcript]["rnaseq1"] += opendict[transcript]
			if transcript in covdict:
				transcript_dict[transcript]["rnaseq1_cov"] += covdict[transcript]
	if rnaseq2_filepath:
		groupname += (rnaseq2_filepath.split("/")[-1]).replace(".sqlite","")+"_"
		if os.path.isfile(rnaseq2_filepath):
			sqlite_db = SqliteDict(rnaseq2_filepath, autocommit=False)
		else:
			return "error","File not found, please report this to tripsvizsite@gmail.com or via the contact page."
		if region == "fiveprime":
			opendict =sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict =sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict =sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict =sqlite_db["{}_all_totals".format(ambig_type)]
		if mapped_reads_norm == True:
			rnaseq2_tot_reads += float(sqlite_db["noncoding_counts"])
			rnaseq2_tot_reads += float(sqlite_db["coding_counts"])
		if min_cov > 0:
			if "{}_all_coverage".format(ambig_type) not in sqlite_db:
				calculate_coverages(sqlite_db,longest_tran_list,ambig_type, region,traninfo_dict)
			if region == "fiveprime":
				covdict =sqlite_db["{}_fiveprime_coverage".format(ambig_type)]
			elif region == "cds":
				covdict =sqlite_db["{}_cds_coverage".format(ambig_type)]
			elif region == "threeprime":
				covdict =sqlite_db["{}_threeprime_coverage".format(ambig_type)]
			elif region == "all":
				covdict =sqlite_db["{}_all_coverage".format(ambig_type)]	
		sqlite_db.close()
		for transcript in longest_tran_list:
			if transcript not in transcript_dict:
				transcript_dict[transcript] = {"riboseq1":0.0001,"riboseq2":0.0001, "rnaseq1":0.0001, "rnaseq2":0.0001,
												"riboseq1_cov":0.0001,"riboseq2_cov":0.0001, "rnaseq1_cov":0.0001, "rnaseq2_cov":0.0001,
												"max_cov":0}
			if transcript in opendict:
				transcript_dict[transcript]["rnaseq2"] += opendict[transcript]
			if transcript in covdict:
				transcript_dict[transcript]["rnaseq2_cov"] += covdict[transcript]


	for transcript in transcript_dict:
		if transcript in ["ENSMUST00000055032","ENSMUST00000091197","ENSMUST00000069309"]:
			print transcript,"ribo1:",transcript_dict[transcript]["riboseq1"],"ribo2:",transcript_dict[transcript]["riboseq2"]
	current_min_reads_list = []
	diff_expressed = []

	if mapped_reads_norm == True:
		ribo1_modifier = riboseq1_tot_reads/riboseq2_tot_reads
		if ribo1_modifier < 1:
			ribo1_modifier = 1

		ribo2_modifier = riboseq2_tot_reads/riboseq1_tot_reads
		if ribo2_modifier < 1:
			ribo2_modifier = 1

		rna1_modifier = rnaseq1_tot_reads/rnaseq2_tot_reads
		if rna1_modifier < 1:
			rna1_modifier = 1

		rna2_modifier = rnaseq2_tot_reads/rnaseq1_tot_reads
		if rna2_modifier < 1:
			rna2_modifier = 1
	else:
		ribo1_modifier = 1
		ribo2_modifier = 1
		rna1_modifier = 1
		rna2_modifier = 1
	
	del_list = []
	for transcript in transcript_dict:
		if label == "TE":
			current_min_reads = min(transcript_dict[transcript]["riboseq1"]/ribo1_modifier,transcript_dict[transcript]["riboseq2"]/ribo2_modifier,transcript_dict[transcript]["rnaseq1"]/rna1_modifier,transcript_dict[transcript]["rnaseq2"]/rna2_modifier)
		elif label == "Riboseq":
			current_min_reads = min(transcript_dict[transcript]["riboseq1"]/ribo1_modifier,transcript_dict[transcript]["riboseq2"]/ribo2_modifier)
		elif label == "Rnaseq":
			current_min_reads = min(transcript_dict[transcript]["rnaseq1"]/rna1_modifier,transcript_dict[transcript]["rnaseq2"]/rna2_modifier)
		if minreads != 0:
			#print "minreads not zero", current_min_reads
			if current_min_reads < minreads:
				#print "deleting transcript beacuase its reads {} were less than the min {}".format(current_min_reads, minreads)
				del_list.append(transcript)
	
	for transcript in transcript_dict:
		if min_cov > 0:
			max_cov = max(transcript_dict[transcript]["riboseq1_cov"],transcript_dict[transcript]["riboseq2_cov"],transcript_dict[transcript]["rnaseq1_cov"],transcript_dict[transcript]["rnaseq2_cov"])
			if max_cov > transcript_dict[transcript]["max_cov"]:
				transcript_dict[transcript]["max_cov"] = max_cov
		if transcript_dict[transcript]["riboseq1"] >= 1:
			transcript_dict[transcript]["riboseq1"] = log((transcript_dict[transcript]["riboseq1"]/ribo1_modifier),2)
		if transcript_dict[transcript]["riboseq2"] >= 1:
			transcript_dict[transcript]["riboseq2"] = log((transcript_dict[transcript]["riboseq2"]/ribo2_modifier),2)
		if transcript_dict[transcript]["rnaseq1"] >= 1:
			transcript_dict[transcript]["rnaseq1"] = log((transcript_dict[transcript]["rnaseq1"]/rna1_modifier),2)
		if transcript_dict[transcript]["rnaseq2"] >= 1:
			transcript_dict[transcript]["rnaseq2"] = log((transcript_dict[transcript]["rnaseq2"]/rna2_modifier),2)
	

	for transcript in master_dict:
		if transcript not in transcript_dict:
			del_list.append(transcript)
	for transcript in del_list:
		if transcript in transcript_dict:
			del transcript_dict[transcript]
		if transcript in master_dict:
			del master_dict[transcript]

	for transcript in transcript_dict:
		gene = "unknown"
		if transcript in traninfo_dict:
			if "gene" in traninfo_dict[transcript]:
				gene = traninfo_dict[transcript]["gene"]
		skip = False
		if label == "TE":
			product = abs(transcript_dict[transcript]["riboseq1"])*abs(transcript_dict[transcript]["riboseq2"])*abs(transcript_dict[transcript]["rnaseq1"])*abs(transcript_dict[transcript]["rnaseq2"])
			geometric_mean = product**(1/float(4))
			te1 = (float(transcript_dict[transcript]["riboseq1"]))-(float(transcript_dict[transcript]["rnaseq1"]))
			te2 = (float(transcript_dict[transcript]["riboseq2"]))-(float(transcript_dict[transcript]["rnaseq2"]))
			fold_change = te2-te1
			current_min_reads_list.append([transcript, geometric_mean, fold_change, gene,skip])
		
		elif label == "Riboseq":
			product = abs(transcript_dict[transcript]["riboseq1"])*abs(transcript_dict[transcript]["riboseq2"])
			geometric_mean = product**(1/float(2))
			ribo1 = float(transcript_dict[transcript]["riboseq1"])
			ribo2 = float(transcript_dict[transcript]["riboseq2"])
			fold_change = ribo2-ribo1
			current_min_reads_list.append([transcript, geometric_mean, fold_change, gene, skip])

		elif label == "Rnaseq":
			product = abs(transcript_dict[transcript]["rnaseq1"])*abs(transcript_dict[transcript]["rnaseq2"])
			geometric_mean = product**(1/float(2))
			rna1 = float(transcript_dict[transcript]["rnaseq1"])
			rna2 = float(transcript_dict[transcript]["rnaseq2"])
			fold_change = rna2-rna1
			current_min_reads_list.append([transcript, geometric_mean, fold_change, gene, skip])


	positive_fc_final_list = []
	negative_fc_final_list = []
	
	sorted_current_min_reads_list = sorted(current_min_reads_list, key=lambda x: x[1])
	bin_list = []
	for i in range(0, len(sorted_current_min_reads_list), 300):
		# if we are not near the end of the list calculate mean and std dev
		if i < (len(sorted_current_min_reads_list)-300):
			#for every transcript in this bin add the min exp to bin_values list
			bin_values = []
			for x in range(i,i+300):
				bin_values.append(sorted_current_min_reads_list[x][2])
			#work out the mean and standard deviation of bin_values
			mean = np.mean(bin_values)
			standard_dev = np.std(bin_values)
			# Append mean and std dev to sorted_current_min_reads_list so we can work out z-score for each gene later
			for x in range(i,i+300):
				sorted_current_min_reads_list[x].append(mean)
				sorted_current_min_reads_list[x].append(standard_dev)
			y = minzscore*(standard_dev)
			threshold = y+abs(mean)
			bin_list.append([mean, standard_dev,threshold])
		else:
			bin_list.append(bin_list[-1])
			# Append mean and std dev to sorted_current_min_reads_list so we can work out z-score for each gene later
			for x in range(i,len(sorted_current_min_reads_list)):
				sorted_current_min_reads_list[x].append(bin_list[-1][0])
				sorted_current_min_reads_list[x].append(bin_list[-1][1])
	for row in sorted_current_min_reads_list:
		tran = row[0]
		if tran not in master_dict:
			master_dict[tran] = {"skip":False}
		if row[4] == True:
			master_dict[tran]["skip"] = True
		else:
			master_dict[tran][groupname] = {"geometric_mean":row[1],"fold_change":row[2],"gene":row[3],"mean":row[5],"standard_dev":row[6]}
	transcript_dict["ribo1_modifier"] = ribo1_modifier
	transcript_dict["ribo2_modifier"] = ribo2_modifier
	transcript_dict["rna1_modifier"] = rna1_modifier
	transcript_dict["rna2_modifier"] = rna2_modifier
	return transcript_dict,groupname

