from flask import Blueprint, render_template, abort, request,jsonify, url_for
from flask import current_app as app
import sqlite3
from sqlitedict import SqliteDict
import ast
import os
import time
import numpy as np
import pandas as pd
from numpy.random import multinomial, random
from bisect import bisect_left
from random import shuffle
from math import log
import config
import subprocess
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code,fetch_user
import riboflask_diff
import collections
from flask_login import current_user
import json



# Differential expression page, used to find diffentially expressed genes via z-score
diff_plotpage_blueprint = Blueprint("diffpage", __name__, template_folder="templates")
@diff_plotpage_blueprint.route('/<organism>/<transcriptome>/differential/')
def diffpage(organism,transcriptome):
	#global user_short_passed
	user_short_passed = True
	global local
	try:
		print (local)
	except:
		local = False

	organism = str(organism)
	connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
	connection.text_factory = str
	cursor = connection.cursor()
	
	user,logged_in = fetch_user()	
	
	
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
				 "transcript_list":str(request.args.get('transcript_list')),
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


def prepare_return_str(input_string):
	return input_string
	

# Creates/serves the z-score plot for differential expression
diffquery_blueprint = Blueprint("diffquery", __name__, template_folder="templates")
@diffquery_blueprint.route('/diffquery', methods=['POST'])
def diffquery():
	#global user_short_passed
	user_short_passed = True
	data = json.loads(request.data)
	plottype = data["plottype"]

	genetype = data["genetype"]
	region = data["region"] #can be all, cds,fiveprime, or threeprime
	if genetype != "coding" and region != "all":
		return prepare_return_str( "If gene type is not set to 'coding' then region has to be set to 'all'")
	html_args = data["html_args"]
	organism = data["organism"]
	transcriptome = data["transcriptome"]
	store_de_results = data["store_de_results"]
	cond_desc = data["cond_desc"]
	gene_list = data["gene_list"]
	minreads = float(data["minreads"])
	transcript_list = ((data["transcript_list"].strip(" ")).replace(" ",",")).split(",")
	if data["min_cov"] != "undefined":
		min_cov = float(data["min_cov"])
	else:
		min_cov = 0
	if min_cov > 1:
		return prepare_return_str( "Minimum coverage should be a value between 0 and 1")
	filename = organism+"_differential_translation_"+str(time.time())+".csv"
	csv_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
	master_file_dict = data["master_file_dict"]
	master_transcript_dict = {}
	connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
	connection.text_factory = str
	cursor = connection.cursor()
	
	user,logged_in = fetch_user()	
	
	
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]

	if owner == 1:
		if os.path.isfile("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome)):
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			return prepare_return_str( "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome))
	else:
		transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
	trancursor = transhelve.cursor()
	if transcript_list == ['']:
		if genetype == "all":
			trancursor.execute("SELECT * from transcripts WHERE principal = 1")
		elif genetype == "coding":
			trancursor.execute("SELECT * from transcripts WHERE principal = 1 AND (tran_type = 'coding' OR tran_type = 1)")
		elif genetype == "noncoding":
			trancursor.execute("SELECT * from transcripts WHERE principal = 1 AND (tran_type = 'noncoding' OR tran_type = 0)")
	else:
		trancursor.execute("SELECT * from transcripts WHERE transcript IN ({});".format(str(transcript_list).strip("[]")))
	result =  trancursor.fetchall()
	traninfo_dict = {}
	for row in result:
		traninfo_dict[row[0]] = {"transcript":str(row[0]) , "gene":row[1], "length":row[2] , "cds_start":row[3] , "cds_stop":row[4] , "seq":row[5] ,
				"strand":row[6], "stop_list":row[7].split(","),"start_list":row[8].split(","), "exon_junctions":row[9].split(","),
				"tran_type":row[10], "principal":row[11]}
	transhelve.close()
	longest_tran_list = traninfo_dict.keys()

	#This will hold one or more z-scores for every transcript
	master_dict = {}

	#keys are genes, values are the x and y values 
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
	deseq_basename = "{}/static/tmp/{}".format(config.SCRIPT_LOC,filename)
	deseq_count_filename = "{}/static/tmp/{}_deseq_counts.csv".format(config.SCRIPT_LOC,filename)
	deseq_sample_filename = "{}/static/tmp/{}_deseq_sample_info.csv".format(config.SCRIPT_LOC,filename)
	deseq_count_file = open(deseq_count_filename,"w")
	deseq_sample_file = open(deseq_sample_filename,"w")
	anota2seq = False

	minzscore = float(data["minzscore"])
	if "mapped_reads_norm" in data:
		mapped_reads_norm = True
	else:
		mapped_reads_norm = False
	#DESeq2 requires all genes be included
	if plottype == "deseq2":
		minreads = 0
		mapped_reads_norm = False
	if plottype == "anota2seq":
		anota2seq = True
		minreads = 0
		mapped_reads_norm = False
	else:
		anota2seq = False
		
	if "ambiguous" in data:
		ambiguous = True
	else:
		ambiguous = False
	
	if len(master_file_dict["riboseq1"]["file_ids"]) == 0 and len(master_file_dict["riboseq2"]["file_ids"]) == 0 and len(master_file_dict["rnaseq1"]["file_ids"]) == 0 and len(master_file_dict["rnaseq2"]["file_ids"]) == 0:
		return prepare_return_str("Error: No files selected.")
	# User can decide to look at just riboseq fold-change, rnaseq fold-change or TE fold-change
	if len(master_file_dict["riboseq1"]["file_ids"]) != len(master_file_dict["riboseq2"]["file_ids"]):
		return prepare_return_str("Error: Both Riboseq Condition boxes need to have an equal number of files.")
	if len(master_file_dict["rnaseq1"]["file_ids"]) != len(master_file_dict["rnaseq2"]["file_ids"]):
		return prepare_return_str("Error: Both RNA-Seq Condition boxes need to have an equal number of files.")
	if len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and (len(master_file_dict["riboseq1"]["file_ids"]) != len(master_file_dict["rnaseq1"]["file_ids"])):
		return prepare_return_str( "Error: If RNA-Seq boxes are not empty they need to have an equal number of files as the riboseq boxes")
	if len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq2"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq2"]["file_ids"]) >= 1:
		label = "TE"
	elif len(master_file_dict["riboseq1"]["file_ids"]) >= 1 and len(master_file_dict["riboseq2"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq1"]["file_ids"]) == 0 and len(master_file_dict["rnaseq2"]["file_ids"]) == 0:
		label = "Riboseq"
	elif len(master_file_dict["riboseq1"]["file_ids"]) == 0 and len(master_file_dict["riboseq2"]["file_ids"]) == 0 and len(master_file_dict["rnaseq1"]["file_ids"]) >= 1 and len(master_file_dict["rnaseq2"]["file_ids"]) >= 1:
		label = "Rnaseq"
	else:
		label = "TE"
		return prepare_return_str("ERROR IMBALANCED OR NO FILES: Either all 4 boxes (RIBO-seq files 1, RIBO-seq files 2, mRNA-seq files 1, mRNA-seq files 2) must have a file associated with it OR both riboseq boxes OR both rnaseq boxes, you currently have {} files in RIBO-seq condition 1 files, {} in RIBO-seq condition 2 files, {} in mRNA-seq condition 1 files and {} in mRNA-seq condition 2 files".format(len(riboseq1_filepaths),len(riboseq2_filepaths),len(rnaseq1_filepaths),len(rnaseq2_filepaths)))

	no_groups = max(len(master_file_dict["riboseq1"]["file_ids"]),len(master_file_dict["rnaseq1"]["file_ids"]))
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
		#if plottype == "z_score":
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
			return prepare_return_str("Error: {}".format(groupname))
		master_transcript_dict[groupname] = transcript_dict

	deseq_sample_file.write(",Condition,SeqType\n")
	headers = ["Gene/Transcript"]
	sample_rows = []
	sample_data = []
	counts = []
	
	if plottype == "deseq2":
		if no_groups <= 1:
			return prepare_return_str("At least two replicates are required when using DESeq2")
		mapped_reads_norm = False
		for i in range(1,no_groups+1):
			if label in ["TE","Riboseq"]:
				deseq_sample_file.write("Ribo_cond1_count{},control,riboseq\n".format(i))
				headers.append("Ribo_cond1_count{}".format(i))
			sample_data.append(["control","riboseq"])
			sample_rows.append("Ribo_cond1_count{}".format(i))
			if label in ["TE","Riboseq"]:
				deseq_sample_file.write("Ribo_cond2_count{},treatment,riboseq\n".format(i))
				headers.append("Ribo_cond2_count{}".format(i))
			sample_data.append(["treatment","riboseq"])
			sample_rows.append("Ribo_cond2_count{}".format(i))
			if label in ["TE","Rnaseq"]:
				deseq_sample_file.write("mRNA_cond1_count{},control,rnaseq\n".format(i))
				headers.append("mRNA_cond1_count{}".format(i))
			sample_data.append(["control","rnaseq"])
			sample_rows.append("mRNA_cond1_count{}".format(i))
			if label in ["TE","Rnaseq"]:
				deseq_sample_file.write("mRNA_cond2_count{},treatment,rnaseq\n".format(i))
				headers.append("mRNA_cond2_count{}".format(i))
			sample_data.append(["treatment","rnaseq"])
			sample_rows.append("mRNA_cond2_count{}".format(i))
		deseq_sample_file.close()
		deseq_count_file.write("{}\n".format(str(headers).strip("[]").replace("'","")))
		deseq_dict = {}
		count_dict = {}
		sample_df = pd.DataFrame(sample_data, columns = ["Condition","SeqType"],index=sample_rows)
		for tran in master_dict:
			gene = (traninfo_dict[tran]["gene"]).replace(",","-")
			if tran not in count_dict:
				count_dict[tran] = {"ribo1":0,"ribo2":0,"rna1":0,"rna2":0}
			deseq_count_file.write("{}___{},".format(gene, tran))
			count_list = [tran]
			group_count = 1
			for group in master_transcript_dict:
				ribo1 = int((2**master_transcript_dict[group][tran]["riboseq1"])-0.0001)
				ribo2 = int((2**master_transcript_dict[group][tran]["riboseq2"])-0.0001)
				rna1 = int((2**master_transcript_dict[group][tran]["rnaseq1"])-0.0001)
				rna2 = int((2**master_transcript_dict[group][tran]["rnaseq2"])-0.0001)
				count_dict[tran]["ribo1"] += ribo1
				count_dict[tran]["rna1"] += rna1
				count_dict[tran]["ribo2"] += ribo2
				count_dict[tran]["rna2"] += rna2
				if label == "TE":
					deseq_count_file.write("{},{},{},{}".format(ribo1,ribo2,rna1,rna2))
				elif label == "Riboseq":
					deseq_count_file.write("{},{}".format(ribo1,ribo2))
				elif label == "Rnaseq":
					deseq_count_file.write("{},{}".format(rna1,rna2))
				if group_count < no_groups:
					deseq_count_file.write(",")
				group_count += 1
			counts.append(count_list)
			deseq_count_file.write("\n")
		deseq_count_file.close()
		subprocess.call("{}/R/deseq2_triplicates.R {} {} {} {} {}".format(config.SCRIPT_LOC, deseq_basename,deseq_count_filename, deseq_sample_filename, label,no_groups),shell=True)
		if label == "TE":
			subprocess.call("tar -C {1}/static/tmp/ -czvf {1}/static/tmp/{0}.tar.gz  {0}_deseq_counts.csv {0}_deseq_sample_info.csv {0}_DESeq2_RIBOSEQ.txt {0}_DESeq2_RNASEQ.txt {0}_DESeq2_TE.txt".format(filename,config.SCRIPT_LOC, deseq_count_filename, deseq_sample_filename),shell=True)
		elif label == "Riboseq":
			subprocess.call("tar -C {1}/static/tmp/ -czvf {1}/static/tmp/{0}.tar.gz  {0}_deseq_counts.csv {0}_deseq_sample_info.csv {0}_DESeq2_RIBOSEQ.txt".format(filename, config.SCRIPT_LOC, deseq_count_filename, deseq_sample_filename),shell=True)
		elif label == "Rnaseq":
			subprocess.call("tar -C {1}/static/tmp/ -czvf {1}/static/tmp/{0}.tar.gz  {0}_deseq_counts.csv {0}_deseq_sample_info.csv {0}_DESeq2_RNASEQ.txt".format(filename, config.SCRIPT_LOC, deseq_count_filename, deseq_sample_filename),shell=True)
		#Parse TE file to get the x,y values of the transcripts and get the genes below a signifigcance level.
		if os.path.isfile("{0}/static/tmp/{1}_DESeq2_TE.txt".format(config.SCRIPT_LOC,filename)):
			te_file = open("{0}/static/tmp/{1}_DESeq2_TE.txt".format(config.SCRIPT_LOC,filename)).readlines()
			for line in te_file[1:]:
				splitline = line.split(",")
				transcript = splitline[0].split("___")[1]
				gene = splitline[0].split("___")[0]
				basemean = splitline[1]
				log2fc = splitline[2]
				lfcSE = splitline[3]
				padj = splitline[6].replace("\n","")
				deseq_dict[gene] = {"te_padj":padj,"ribo_padj":"NA","rna_padj":"NA","ribo_fc":"NA","rna_fc":"NA","tran":transcript,
									"ribo1":count_dict[transcript]["ribo1"],"ribo2":count_dict[transcript]["ribo2"], "rna1":count_dict[transcript]["rna1"],
									"rna2":count_dict[transcript]["rna2"],"te_basemean":basemean,"te_lfcSE":lfcSE,"ribo_basemean":"NA","ribo_lfcSE":"NA","rna_basemean":"NA","rna_lfcSE":"NA"}
		if os.path.isfile("{0}/static/tmp/{1}_DESeq2_RIBOSEQ.txt".format(config.SCRIPT_LOC,filename)):
			ribo_file = open("{0}/static/tmp/{1}_DESeq2_RIBOSEQ.txt".format(config.SCRIPT_LOC,filename)).readlines()
			for line in ribo_file[1:]:
				splitline = line.split(",")
				transcript = splitline[0].split("___")[1]
				gene = splitline[0].split("___")[0]
				basemean = splitline[1]
				log2fc = splitline[2]
				lfcSE = splitline[3]
				padj = splitline[6].replace("\n","")
				if gene not in deseq_dict:
					deseq_dict[gene] = {"te_padj":"NA","ribo_padj":padj,"rna_padj":"NA","ribo_fc":log2fc,"rna_fc":"NA","tran":transcript,
						"ribo1":count_dict[transcript]["ribo1"],"ribo2":count_dict[transcript]["ribo2"], "rna1":count_dict[transcript]["rna1"],
						"rna2":count_dict[transcript]["rna2"],"te_basemean":"NA","te_lfcSE":"NA","ribo_basemean":basemean,"ribo_lfcSE":lfcSE,"rna_basemean":"NA","rna_lfcSE":"NA"}
				else:
					deseq_dict[gene]["ribo_padj"] = padj
					deseq_dict[gene]["ribo_fc"] = log2fc
		if os.path.isfile("{0}/static/tmp/{1}_DESeq2_RNASEQ.txt".format(config.SCRIPT_LOC,filename)):
			rna_file = open("{0}/static/tmp/{1}_DESeq2_RNASEQ.txt".format(config.SCRIPT_LOC,filename)).readlines()
			for line in rna_file[1:]:
				splitline = line.split(",")
				transcript = splitline[0].split("___")[1]
				gene = splitline[0].split("___")[0]
				basemean = splitline[1]
				log2fc = splitline[2]
				lfcSE = splitline[3]
				padj = splitline[6].replace("\n","")
				if gene not in deseq_dict:
					deseq_dict[gene] = {"te_padj":"NA","ribo_padj":"NA","rna_padj":padj,"ribo_fc":"NA","rna_fc":log2fc,"tran":transcript,
						"ribo1":count_dict[transcript]["ribo1"],"ribo2":count_dict[transcript]["ribo2"], "rna1":count_dict[transcript]["rna1"],
						"rna2":count_dict[transcript]["rna2"],"te_basemean":"NA","te_lfcSE":"NA","ribo_basemean":"NA","ribo_lfcSE":"NA","rna_basemean":basemean,"rna_lfcSE":lfcSE}
				else:
					deseq_dict[gene]["rna_padj"] = padj
					deseq_dict[gene]["rna_fc"] = log2fc
		
		
		#DELETE GENES FROM DESEQ_DICT THAT DON'T HAVE A padj for either rna and ribo (these have low counts)
		del_list = []
		for gene in deseq_dict:
			if deseq_dict[gene]["rna_padj"] == "NA" and  deseq_dict[gene]["ribo_padj"] == "NA":
				del_list.append(gene)
		
		for gene in del_list:
			del deseq_dict[gene]

		
		
		background_col = config.BACKGROUND_COL
		uga_col = config.UGA_COL
		uag_col = config.UAG_COL
		uaa_col = config.UAA_COL
		title_size = config.TITLE_SIZE
		subheading_size = config.SUBHEADING_SIZE
		axis_label_size = config.AXIS_LABEL_SIZE
		marker_size = config.MARKER_SIZE
		if current_user.is_authenticated:
			user_name = current_user.name
			cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user_name))
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
		
		if html_args["user_short"] == "None" or user_short_passed == True:
			short_code = generate_short_code(data,organism,html_args["transcriptome"],"differential")
		else:
			short_code = html_args["user_short"]
			user_short_passed = True
		
		tar_file = filename+".tar.gz"
		connection.close()
		
		return riboflask_diff.deseq2_plot(deseq_dict,
										organism,
										transcriptome,
										master_file_dict["riboseq1"]["file_ids"],
										master_file_dict["riboseq2"]["file_ids"],
										master_file_dict["rnaseq1"]["file_ids"],
										master_file_dict["rnaseq2"]["file_ids"],
										background_col,
										short_code,
										mapped_reads_norm,
										tar_file,
										no_groups,
										str(title_size)+"pt",
										str(axis_label_size)+"pt",
										str(subheading_size)+"pt",
										str(marker_size)+"pt",
										ambiguous,
										gene_list,
										label,
										minzscore)
	if anota2seq:
		if label != "TE" or no_groups <= 2:
			return prepare_return_str("At least 3 replicates are required for all groups when using Anota2seq")
		for i in range(1,no_groups+1):
			if label in ["TE","Riboseq"]:
				deseq_sample_file.write("Ribo_cond1_count{},control,riboseq\n".format(i))
				headers.append("Ribo_cond1_count{}".format(i))
			sample_data.append(["control","riboseq"])
			sample_rows.append("Ribo_cond1_count{}".format(i))
			if label in ["TE","Riboseq"]:
				deseq_sample_file.write("Ribo_cond2_count{},treatment,riboseq\n".format(i))
				headers.append("Ribo_cond2_count{}".format(i))
			sample_data.append(["treatment","riboseq"])
			sample_rows.append("Ribo_cond2_count{}".format(i))
			if label in ["TE","Rnaseq"]:
				deseq_sample_file.write("mRNA_cond1_count{},control,rnaseq\n".format(i))
				headers.append("mRNA_cond1_count{}".format(i))
			sample_data.append(["control","rnaseq"])
			sample_rows.append("mRNA_cond1_count{}".format(i))
			if label in ["TE","Rnaseq"]:
				deseq_sample_file.write("mRNA_cond2_count{},treatment,rnaseq\n".format(i))
				headers.append("mRNA_cond2_count{}".format(i))
			sample_data.append(["treatment","rnaseq"])
			sample_rows.append("mRNA_cond2_count{}".format(i))
		deseq_sample_file.close()
		deseq_count_file.write("{}\n".format(str(headers).strip("[]").replace("'","")))
		deseq_dict = {}
		count_dict = {}
		sample_df = pd.DataFrame(sample_data, columns = ["Condition","SeqType"],index=sample_rows)
		for tran in master_dict:
			gene = (traninfo_dict[tran]["gene"]).replace(",","-")
			if tran not in count_dict:
				count_dict[tran] = {"ribo1":0,"ribo2":0,"rna1":0,"rna2":0}
			deseq_count_file.write("{}___{},".format(gene, tran))
			count_list = [tran]
			group_count = 1
			for group in master_transcript_dict:
				ribo1 = int((2**master_transcript_dict[group][tran]["riboseq1"])-0.0001)
				ribo2 = int((2**master_transcript_dict[group][tran]["riboseq2"])-0.0001)
				rna1 = int((2**master_transcript_dict[group][tran]["rnaseq1"])-0.0001)
				rna2 = int((2**master_transcript_dict[group][tran]["rnaseq2"])-0.0001)
				count_dict[tran]["ribo1"] += ribo1
				count_dict[tran]["rna1"] += rna1
				count_dict[tran]["ribo2"] += ribo2
				count_dict[tran]["rna2"] += rna2
				if label == "TE":
					deseq_count_file.write("{},{},{},{}".format(ribo1,ribo2,rna1,rna2))
				elif label == "Riboseq":
					deseq_count_file.write("{},{}".format(ribo1,ribo2))
				elif label == "Rnaseq":
					deseq_count_file.write("{},{}".format(rna1,rna2))
				if group_count < no_groups:
					deseq_count_file.write(",")
				group_count += 1
			counts.append(count_list)
			deseq_count_file.write("\n")
		deseq_count_file.close()
		subprocess.call("{}/R/anota2seq.R {} {} {} {} {} {}".format(config.SCRIPT_LOC, deseq_basename,deseq_count_filename, deseq_sample_filename, label,no_groups, minzscore/100),shell=True)
		if label == "TE":
			subprocess.call("tar -C {1}/static/tmp/ -czvf {1}/static/tmp/{0}.tar.gz  {0}_deseq_counts.csv {0}_deseq_sample_info.csv {0}_anota2seq_RIBOSEQ.csv {0}_anota2seq_BUFFERED.csv {0}_anota2seq_mRNA.csv".format(filename,config.SCRIPT_LOC, deseq_count_filename, deseq_sample_filename),shell=True)
		elif label == "Riboseq":
			subprocess.call("tar -C {1}/static/tmp/ -czvf {1}/static/tmp/{0}.tar.gz  {0}_deseq_counts.csv {0}_deseq_sample_info.csv {0}_anota2seq_RIBOSEQ.csv".format(filename, config.SCRIPT_LOC, deseq_count_filename, deseq_sample_filename),shell=True)
		elif label == "Rnaseq":
			subprocess.call("tar -C {1}/static/tmp/ -czvf {1}/static/tmp/{0}.tar.gz  {0}_deseq_counts.csv {0}_deseq_sample_info.csv {0}_anota2seq_mRNA.csv".format(filename, config.SCRIPT_LOC, deseq_count_filename, deseq_sample_filename),shell=True)
		#Parse TE file to get the x,y values of the transcripts and get the genes below a signifigcance level.
		sig_translated = []
		sig_rna = []
		sig_buffered = []
		if os.path.isfile("{0}/static/tmp/{1}_anota2seq_sig_translated_genes.csv".format(config.SCRIPT_LOC,filename)):
			sig_translated_file = open("{0}/static/tmp/{1}_anota2seq_sig_translated_genes.csv".format(config.SCRIPT_LOC,filename)).readlines()
			for line in sig_translated_file[1:]:
				splitline = line.split(",")
				gene = splitline[1].split("___")[0]
				sig_translated.append(gene)
		if os.path.isfile("{0}/static/tmp/{1}_anota2seq_sig_rna_genes.csv".format(config.SCRIPT_LOC,filename)):
			sig_rna_file = open("{0}/static/tmp/{1}_anota2seq_sig_rna_genes.csv".format(config.SCRIPT_LOC,filename)).readlines()
			for line in sig_rna_file[1:]:
				splitline = line.split(",")
				gene = splitline[1].split("___")[0]
				sig_rna.append(gene)
		if os.path.isfile("{0}/static/tmp/{1}_anota2seq_sig_buffered_genes.csv".format(config.SCRIPT_LOC,filename)):
			sig_buffered_file = open("{0}/static/tmp/{1}_anota2seq_sig_buffered_genes.csv".format(config.SCRIPT_LOC,filename)).readlines()
			for line in sig_buffered_file[1:]:
				splitline = line.split(",")
				gene = splitline[1].split("___")[0]
				sig_buffered.append(gene)
		if os.path.isfile("{0}/static/tmp/{1}_anota2seq_RIBOSEQ.csv".format(config.SCRIPT_LOC,filename)):
			ribo_file = open("{0}/static/tmp/{1}_anota2seq_RIBOSEQ.csv".format(config.SCRIPT_LOC,filename)).readlines()
			for line in ribo_file[1:]:
				splitline = line.split(",")
				transcript = splitline[1].split("___")[1]
				gene = splitline[1].split("___")[0]
				basemean = splitline[1]
				log2fc = splitline[10]
				lfcSE = splitline[3]
				padj = splitline[9]
				if gene not in deseq_dict:
					deseq_dict[gene] = {"te_padj":"NA","ribo_padj":padj,"rna_padj":"NA","ribo_fc":log2fc,"rna_fc":"NA","tran":transcript,
						"ribo1":count_dict[transcript]["ribo1"],"ribo2":count_dict[transcript]["ribo2"], "rna1":count_dict[transcript]["rna1"],
						"rna2":count_dict[transcript]["rna2"],"te_basemean":"NA","te_lfcSE":"NA","ribo_basemean":basemean,"ribo_lfcSE":lfcSE,"rna_basemean":"NA","rna_lfcSE":"NA"}
				else:
					deseq_dict[gene]["ribo_padj"] = padj
					deseq_dict[gene]["ribo_fc"] = log2fc
		if os.path.isfile("{0}/static/tmp/{1}_anota2seq_mRNA.csv".format(config.SCRIPT_LOC,filename)):
			rna_file = open("{0}/static/tmp/{1}_anota2seq_mRNA.csv".format(config.SCRIPT_LOC,filename)).readlines()
			for line in rna_file[1:]:
				splitline = line.split(",")
				transcript = splitline[1].split("___")[1]
				gene = splitline[1].split("___")[0]
				basemean = splitline[1]
				log2fc = splitline[8]
				lfcSE = splitline[3]
				padj = splitline[7]
				if gene not in deseq_dict:
					deseq_dict[gene] = {"te_padj":"NA","ribo_padj":"NA","rna_padj":padj,"ribo_fc":"NA","rna_fc":log2fc,"tran":transcript,
						"ribo1":count_dict[transcript]["ribo1"],"ribo2":count_dict[transcript]["ribo2"], "rna1":count_dict[transcript]["rna1"],
						"rna2":count_dict[transcript]["rna2"],"te_basemean":"NA","te_lfcSE":"NA","ribo_basemean":"NA","ribo_lfcSE":"NA","rna_basemean":basemean,"rna_lfcSE":lfcSE}
				else:
					deseq_dict[gene]["rna_padj"] = padj
					deseq_dict[gene]["rna_fc"] = log2fc

		
		background_col = config.BACKGROUND_COL
		uga_col = config.UGA_COL
		uag_col = config.UAG_COL
		uaa_col = config.UAA_COL
		title_size = config.TITLE_SIZE
		subheading_size = config.SUBHEADING_SIZE
		axis_label_size = config.AXIS_LABEL_SIZE
		marker_size = config.MARKER_SIZE
		if current_user.is_authenticated:
			user_id = current_user.id
			#get a list of organism id's this user can access
			cursor.execute("SELECT background_col,title_size,subheading_size,axis_label_size,marker_size from user_settings WHERE user_id = '{}';".format(user_id))
			result = (cursor.fetchone())
			background_col = result[0]
			title_size = result[1]
			subheading_size = result[2]
			axis_label_size = result[3]
			marker_size = result[4]
			
		if html_args["user_short"] == "None" or user_short_passed == True:
			short_code = generate_short_code(data,organism,html_args["transcriptome"],"differential")
		else:
			short_code = html_args["user_short"]
			user_short_passed = True
		
		tar_file = filename+".tar.gz"
		connection.close()
			
		
		return riboflask_diff.anota2seq_plot(deseq_dict,
										organism,
										transcriptome,
										master_file_dict["riboseq1"]["file_ids"],
										master_file_dict["riboseq2"]["file_ids"],
										master_file_dict["rnaseq1"]["file_ids"],
										master_file_dict["rnaseq2"]["file_ids"],
										background_col,
										short_code,
										mapped_reads_norm,
										tar_file,
										no_groups,
										str(title_size)+"pt",
										str(axis_label_size)+"pt",
										str(subheading_size)+"pt",
										str(marker_size)+"pt",
										ambiguous,
										gene_list,
										label,
										minzscore,
										sig_translated,
										sig_rna,
										sig_buffered)

	
	
	
	
	
		'''
		
		read_delim = robjects.r('read.delim')
		as_data_frame = robjects.r('as.data.frame')
		as_factor = robjects.r('as.factor')
		rapply = robjects.r("apply")
		relevel = robjects.r("relevel")
		to_dataframe = robjects.r('function(x) data.frame(x)')
		deseq = importr('DESeq2')
		count_df = pd.read_csv(deseq_count_filename)
		design_df = pd.read_csv(deseq_sample_filename,index_col=0)
		count_df=count_df.dropna()
		design = Formula("~ Condition + SeqType + Condition:SeqType")
		gene_id = count_df["id"]
		count_matrix = pandas2ri.py2ri(count_df.drop("id",axis=1))
		design_matrix = pandas2ri.py2ri(design_df)
		dds = deseq.DESeqDataSetFromMatrix(countData=count_matrix,colData=design_matrix,design=design)
		dds = deseq.DESeq(dds)
		try:
			dollar = base.__dict__["$"]
			dds.SeqType = relevel(dds.SeqType,"rnaseq")
		except Exception as e:
			pass
		#normalized_count_matrix = deseq.counts(dds, normalized=True)
		comparison = deseq.resultsNames(dds)
		deseq_result = deseq.results(dds)
		deseq_result = to_dataframe(deseq_result)
		deseq_result = pandas2ri.ri2py(deseq_result) ## back to pandas dataframe
		deseq_result["id"] = gene_id.values
		
		return riboflask_diff.deseq2_plot(ribo_vs_rna_dict,
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
		'''

		





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
						fdr = round(valb/(vala),4)
					except:
						fdr = 0.00001
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
	background_col = config.BACKGROUND_COL
	uga_col = config.UGA_COL
	uag_col = config.UAG_COL
	uaa_col = config.UAA_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE
	if current_user.is_authenticated:
		user_id = current_user.id
		#get a list of organism id's this user can access
		cursor.execute("SELECT background_col,title_size,subheading_size,axis_label_size,marker_size from user_settings WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchone())
		background_col = result[0]
		title_size = result[1]
		subheading_size = result[2]
		axis_label_size = result[3]
		marker_size = result[4]
		
	if html_args["user_short"] == "None" or user_short_passed == True:
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"differential")
	else:
		short_code = html_args["user_short"]
		user_short_passed = True
	connection.close()
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
			return prepare_return_str("error","File not found: {}".format(groupname))
		if region == "fiveprime":
			opendict = sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict = sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict = sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict = sqlite_db["{}_all_totals".format(ambig_type)]
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
			return prepare_return_str("error","File not found, please report this to tripsvizsite@gmail.com or via the contact page.")
		if region == "fiveprime":
			opendict = sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict = sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict = sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict = sqlite_db["{}_all_totals".format(ambig_type)]
		if mapped_reads_norm == True:
			riboseq2_tot_reads += float(sqlite_db["noncoding_counts"])
			riboseq2_tot_reads += float(sqlite_db["coding_counts"])
		if min_cov > 0:
			if "{}_all_coverage".format(ambig_type) not in sqlite_db:
				calculate_coverages(sqlite_db,longest_tran_list,ambig_type, region,traninfo_dict)
			if region == "fiveprime":
				covdict = sqlite_db["{}_fiveprime_coverage".format(ambig_type)]
			elif region == "cds":
				covdict = sqlite_db["{}_cds_coverage".format(ambig_type)]
			elif region == "threeprime":
				covdict = sqlite_db["{}_threeprime_coverage".format(ambig_type)]
			elif region == "all":
				covdict = sqlite_db["{}_all_coverage".format(ambig_type)]
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
			return prepare_return_str("error","File not found, please report this to tripsvizsite@gmail.com or via the contact page.")
		if region == "fiveprime":
			opendict = sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict = sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict = sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict = sqlite_db["{}_all_totals".format(ambig_type)]
		if mapped_reads_norm == True:
			rnaseq1_tot_reads += float(sqlite_db["noncoding_counts"])
			rnaseq1_tot_reads += float(sqlite_db["coding_counts"])
		if min_cov > 0:
			if "{}_all_coverage".format(ambig_type) not in sqlite_db:
				calculate_coverages(sqlite_db,longest_tran_list,ambig_type, region,traninfo_dict)
			if region == "fiveprime":
				covdict = sqlite_db["{}_fiveprime_coverage".format(ambig_type)]
			elif region == "cds":
				covdict = sqlite_db["{}_cds_coverage".format(ambig_type)]
			elif region == "threeprime":
				covdict = sqlite_db["{}_threeprime_coverage".format(ambig_type)]
			elif region == "all":
				covdict = sqlite_db["{}_all_coverage".format(ambig_type)]
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
			return prepare_return_str("error","File not found, please report this to tripsvizsite@gmail.com or via the contact page.")
		if region == "fiveprime":
			opendict = sqlite_db["{}_fiveprime_totals".format(ambig_type)]
		elif region == "cds":
			opendict = sqlite_db["{}_cds_totals".format(ambig_type)]
		elif region == "threeprime":
			opendict = sqlite_db["{}_threeprime_totals".format(ambig_type)]
		elif region == "all":
			opendict = sqlite_db["{}_all_totals".format(ambig_type)]
		if mapped_reads_norm == True:
			rnaseq2_tot_reads += float(sqlite_db["noncoding_counts"])
			rnaseq2_tot_reads += float(sqlite_db["coding_counts"])
		if min_cov > 0:
			if "{}_all_coverage".format(ambig_type) not in sqlite_db:
				calculate_coverages(sqlite_db,longest_tran_list,ambig_type, region,traninfo_dict)
			if region == "fiveprime":
				covdict = sqlite_db["{}_fiveprime_coverage".format(ambig_type)]
			elif region == "cds":
				covdict = sqlite_db["{}_cds_coverage".format(ambig_type)]
			elif region == "threeprime":
				covdict = sqlite_db["{}_threeprime_coverage".format(ambig_type)]
			elif region == "all":
				covdict = sqlite_db["{}_all_coverage".format(ambig_type)]
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
			if current_min_reads < minreads:
				del_list.append(transcript)

	for transcript in transcript_dict:
		if min_cov > 0:
			max_cov = max(transcript_dict[transcript]["riboseq1_cov"],transcript_dict[transcript]["riboseq2_cov"],transcript_dict[transcript]["rnaseq1_cov"],transcript_dict[transcript]["rnaseq2_cov"])
			if max_cov > transcript_dict[transcript]["max_cov"]:
				transcript_dict[transcript]["max_cov"] = max_cov
		if transcript_dict[transcript]["riboseq1"] != 0:
			transcript_dict[transcript]["riboseq1"] = log((transcript_dict[transcript]["riboseq1"]/ribo1_modifier),2)
		if transcript_dict[transcript]["riboseq2"] != 0:
			transcript_dict[transcript]["riboseq2"] = log((transcript_dict[transcript]["riboseq2"]/ribo2_modifier),2)
		if transcript_dict[transcript]["rnaseq1"] != 0:
			transcript_dict[transcript]["rnaseq1"] = log((transcript_dict[transcript]["rnaseq1"]/rna1_modifier),2)
		if transcript_dict[transcript]["rnaseq2"] != 0:
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
			# Remove log2 tranformation before getting geometic mean, otherwise
			# counts just above 0 will have very large absolute values and
			# result in an incorrect geometric mean
			ribo1_count = transcript_dict[transcript]["riboseq1"]
			ribo2_count = transcript_dict[transcript]["riboseq2"]
			rna1_count = transcript_dict[transcript]["rnaseq1"]
			rna2_count = transcript_dict[transcript]["rnaseq2"]
			if ribo1_count < 1:
				ribo1_count = 1
			if ribo2_count < 1:
				ribo2_count = 1
			if rna1_count < 1:
				rna1_count = 1
			if rna2_count < 1:
				rna2_count = 1
			product = ribo1_count*ribo2_count*rna1_count*rna2_count
			geometric_mean = product**(1/float(4))
			te1 = (float(ribo1_count))-(float(ribo2_count))
			te2 = (float(rna1_count))-(float(rna2_count))
			fold_change = te2-te1
			current_min_reads_list.append([transcript, geometric_mean, fold_change, gene,skip])

		elif label == "Riboseq":
			# Remove log2 tranformation before getting geometic mean, otherwise
			# counts just above 0 will have very large absolute values and
			# result in an incorrect geometric mean
			ribo1_count = transcript_dict[transcript]["riboseq1"]
			ribo2_count = transcript_dict[transcript]["riboseq2"]
			if ribo1_count < 1:
				ribo1_count = 1
			if ribo2_count < 1:
				ribo2_count = 1
			product = ribo1_count*ribo2_count
			geometric_mean = product**(1/float(2))

			ribo1 = float(ribo1_count)
			ribo2 = float(ribo2_count)
			fold_change = ribo2-ribo1
			current_min_reads_list.append([transcript, geometric_mean, fold_change, gene, skip])

		elif label == "Rnaseq":
			# Remove log2 tranformation before getting geometic mean, otherwise
			# counts just above 0 will have very large absolute values and
			# result in an incorrect geometric mean
			rna1_count = transcript_dict[transcript]["rnaseq1"]
			rna2_count = transcript_dict[transcript]["rnaseq2"]
			if rna1_count < 1:
				rna1_count = 1
			if rna2_count < 1:
				rna2_count = 1
			product = rna1_count*rna2_count
			geometric_mean = product**(1/float(2))
			rna1 = float(rna1_count)
			rna2 = float(rna2_count)
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
