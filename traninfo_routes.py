from flask import Blueprint, render_template, abort, request
import sqlite3
from sqlitedict import SqliteDict
import ast
import os
import time
import re
import operator
from math import log
import config
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code, fetch_user
import traninfo_plots
import collections
from flask_login import current_user
import json



traninfo_plotpage_blueprint = Blueprint("traninfo_plotpage", __name__, template_folder="templates")
@traninfo_plotpage_blueprint.route('/<organism>/<transcriptome>/traninfo_plot/')
def traninfo_plotpage(organism, transcriptome):
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
	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict, accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)

	cursor.execute("SELECT gwips_clade,gwips_organism,gwips_database,default_transcript from organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism,transcriptome))
	result = (cursor.fetchone())
	studyinfo_dict = fetch_study_info(organism)
	gwips_clade = result[0]
	gwips_org = result[1]
	gwips_db = result[2]
	gwips_info = {"organism":gwips_org,
				  "clade": gwips_clade,
				  "database": gwips_db}
	default_tran = result[3]
	# holds all values the user could possibly pass in the url (keywords are after request.args.get), anything not passed by user will be a string: "None"
	html_args = {"user_short":str(request.args.get('short')),
				 "user_plot_type":str(request.args.get('plot')),
				 "gc_tranlist":str(request.args.get('gc_tranlist')),
				 "gc_tranlist2":str(request.args.get('gc_tranlist2')),
				 "gc_tranlist3":str(request.args.get('gc_tranlist3')),
				 "gc_tranlist4":str(request.args.get('gc_tranlist4')),
				 "plot_type":str(request.args.get('plot_type')),
				 "gc_location":str(request.args.get('gc_location')),
				 "nucleotide":str(request.args.get('nucleotide')),
				 "nuc_freq_plot_anchor":str(request.args.get('nuc_freq_plot_anchor')),
				 "nuc_freq_plot_window":str(request.args.get('nuc_freq_plot_window')),
				 "nuc_freq_plot_tranlist":str(request.args.get('nuc_freq_plot_tranlist')),
				 "orftype":str(request.args.get('orftype')),
				 "transcriptome":str(transcriptome),
				 "metagene_tranlist":str(request.args.get('metagene_tranlist'))}
	connection.close()
	return render_template('traninfo_index.html', gwips_clade=gwips_clade, gwips_org=gwips_org, gwips_db=gwips_db,transcriptome=transcriptome,organism=organism,default_tran=default_tran,current_username=user,local=local,
						   studies_dict=accepted_studies,accepted_files=accepted_files,html_args=html_args,
						   studyinfo_dict=studyinfo_dict,seq_types=seq_types)



# Used to create custom metagene plots on the traninformation plot page
def create_custom_metagene(custom_seq_list, exclude_first_val, exclude_last_val, include_first_val, include_last_val, custom_search_region, exclude_first, exclude_last, include_first, include_last,sqlite_db, organism,metagene_tranlist,metagene_frame):
	custom_metagene_id = "cmgc_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}_{}".format((custom_seq_list.upper()).replace(" ","").replace("T","U").replace(",","_"), custom_search_region, exclude_first, exclude_last, include_first, include_last,exclude_first_val,exclude_last_val,include_first_val,include_last_val,metagene_tranlist)
	#try:
	#	mgc = sqlite_db[custom_metagene_id]
	#except:
	#	pass
	transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
	cursor = transhelve.cursor()
	if metagene_tranlist == "":
		cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE principal = 1")
	else:
		metagene_tranlist = metagene_tranlist.split(",")
		cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE transcript IN ({})".format(str(metagene_tranlist).strip("[]")))

	result = cursor.fetchall()
	mgc = {"fiveprime":{},"threeprime":{}}
	iupac_dict =   {"R":"[AG]","Y":"[CU]","S":"[GC]","W":"[AU]","K":"[GU]",
					"M":"[AC]","B":"[CGU]","D":"[AGU]","D":"[AGU]",
					"H":"[ACU]","V":"[ACG]","N":"[AUGC]"}

	subseq_list = []
	for subseq in custom_seq_list.split(","):
		subseq = subseq.upper()
		subseq = subseq.replace("T","U").replace(" ","")
		for code in iupac_dict:
			subseq = subseq.replace(code,iupac_dict[code])
		subseq_list.append(subseq)
	for row in result:
		tran = row[0]
		cds_start = int(row[1])
		cds_stop = int(row[2])
		seq = row[3].replace("T","U")
		if custom_search_region == "whole_gene":
			min_pos = 0
			max_pos = len(seq)
		if custom_search_region == "five_leader":
			min_pos = 0
			max_pos = cds_start
		if custom_search_region == "cds":
			min_pos = cds_start
			max_pos = cds_stop
		if custom_search_region == "three_trailer":
			min_pos = cds_stop
			max_pos = len(seq)
		if include_first == True:
			max_pos = min_pos+include_first_val
		if include_last == True:
			min_pos = max_pos-include_last_val
		if exclude_first == True:
			min_pos = min_pos + exclude_first_val
		if exclude_last == True:
			max_pos = max_pos-exclude_last_val

		if cds_start != "NULL" and cds_start != None:
			seq_positions = []
			for subseq in subseq_list:
				pattern = re.compile(r"{}".format(subseq))
				m = ""
				search_pos = min_pos
				while m != None:
					m = pattern.search(seq,search_pos,max_pos+1)
					if m != None:
						position = m.span()[0]
						seq_positions.append(position)
						search_pos = position+1
		if seq_positions != []:
			profile = {"fiveprime":{},"threeprime":{}}
			try:
				tran_reads = sqlite_db[tran]["unambig"]
			except:
				tran_reads = {"unambig":{}}
			offsets = sqlite_db["offsets"]
			for readlen in tran_reads:
				if readlen in offsets["fiveprime"]["offsets"]:
					offset = offsets["fiveprime"]["offsets"][readlen]
				else:
					offset = 15
				profile["fiveprime"][readlen] = {}
				profile["threeprime"][readlen] = {}
				for pos in tran_reads[readlen]:
					relative_frame = ((pos+offset)-cds_start)%3
					if metagene_frame == "minus_frame":
						if relative_frame != 1:
							continue
					if metagene_frame == "in_frame":
						if relative_frame != 2:
							continue
					if metagene_frame == "plus_frame":
						if relative_frame != 0:
							continue
					three_pos = pos+readlen
					count = tran_reads[readlen][pos]
					try:
						profile["fiveprime"][readlen][pos] += count
					except:
						profile["fiveprime"][readlen][pos] = 0
						profile["fiveprime"][readlen][pos] += count
					try:
						profile["threeprime"][readlen][three_pos] += count
					except:
						profile["threeprime"][readlen][three_pos] = 0
						profile["threeprime"][readlen][three_pos] += count
			for readlen in profile["fiveprime"]:
				if readlen not in mgc["fiveprime"]:
					mgc["fiveprime"][readlen] = {}
					for i in range(-600,601):
						mgc["fiveprime"][readlen][i] = 0
				for pos in profile["fiveprime"][readlen]:
					count = profile["fiveprime"][readlen][pos]
					for seq_position in seq_positions:
						if seq_position >= pos-600 and seq_position <= pos+600:
							relative_seq_position = pos-seq_position
							mgc["fiveprime"][readlen][relative_seq_position] += count

			for readlen in profile["threeprime"]:
				if readlen not in mgc["threeprime"]:
					mgc["threeprime"][readlen] = {}
					for i in range(-600,601):
						mgc["threeprime"][readlen][i] = 0
				for pos in profile["threeprime"][readlen]:
					count = profile["threeprime"][readlen][pos]
					for seq_position in seq_positions:
						if seq_position >= pos-600 and seq_position <= pos+600:
							relative_seq_position = pos-seq_position
							mgc["threeprime"][readlen][relative_seq_position] += count
	sqlite_db[custom_metagene_id] = mgc
	sqlite_db.commit()
	return mgc







traninfoquery_blueprint = Blueprint("traninfoquery", __name__, template_folder="templates")
@traninfoquery_blueprint.route('/traninfoquery', methods=['POST'])
def traninfoquery():
	#global user_short_passed
	user_short_passed = True
	tran_dict = {}
	gene_dict = {}
	data = json.loads(request.data)
	plottype = data["plottype"]

	custom_seq_list = data["custom_seq_list"]
	exclude_first_val = int(data["exclude_first_val"])
	exclude_last_val = int(data["exclude_last_val"])
	include_first_val = int(data["include_first_val"])
	include_last_val = int(data["include_last_val"])
	metagene_tranlist = (data["metagene_tranlist"].strip(" ")).replace(" ",",")

	mismatch_minreadcount = int(data['mismatch_minreadcount'])
	mismatch_minper = int(data['mismatch_minper'])
	mismatch_maxper = int(data['mismatch_maxper'])
	mismatch_maxhit = int(data['mismatch_maxhit'])
	mismatch_region = data['mismatch_region']
	gc_tranlist = data['gc_tranlist'].upper()
	gc_tranlist2 = data['gc_tranlist2'].upper()
	gc_tranlist3 = data['gc_tranlist3'].upper()
	gc_tranlist4 = data['gc_tranlist4'].upper()
	nuc_freq_plot_tranlist = data['nuc_freq_plot_tranlist']
	nuc_freq_plot_window = int(data['nuc_freq_plot_window'])
	maxscaleval = data['maxscaleval']


	if maxscaleval != "None" and maxscaleval != "":
		try:
			maxscaleval = int(maxscaleval)
		except:
			maxscalval = "None"

	organism = data['organism']
	transcriptome = data['transcriptome']
	custom_search_region = data["custom_search_region"]
	contaminant_organism = data["contaminant_organism"]
	nuc_comp_type = data["nuc_comp_type"]
	gc_location = data["gc_location"]
	nucleotide = data["nucleotide"]
	plot_type = data["plot_type"]
	nuc_freq_plot_anchor = data["nuc_freq_plot_anchor"]
	orftype = data["orftype"]
	smooth_amount = int(data["smooth_amount"])
	html_args = data["html_args"]
	nuccomp_reads = data["nuccomp_reads"]
	corr_type = data["corr_type"]
	
	if "principal" in data:
		principal = True
	else:
		principal = False
	
	if "exons" in data:
		exons = True
	else:
		exons = False
	
	connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
	connection.text_factory = str
	cursor = connection.cursor()
	
	user,logged_in = fetch_user()
	
	background_col = config.BACKGROUND_COL
	readlength_col = config.READLENGTH_COL
	metagene_fiveprime_col = config.METAGENE_FIVEPRIME_COL
	metagene_threeprime_col = config.METAGENE_THREEPRIME_COL
	a_col = config.A_COL
	t_col = config.T_COL
	g_col = config.G_COL
	c_col = config.C_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE

	#get a list of organism id's this user can access
	if current_user.is_authenticated:
		#get user_id
		user_name = current_user.name
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user_name))
		result = (cursor.fetchone())
		user_id = result[0]
		cursor.execute("SELECT background_col,readlength_col,metagene_fiveprime_col,metagene_threeprime_col,nuc_comp_a_col,nuc_comp_t_col,nuc_comp_g_col,nuc_comp_c_col,title_size,subheading_size,axis_label_size,marker_size from user_settings WHERE user_id = '{}';".format(user_id))
		result = (cursor.fetchone())
		background_col = result[0]
		readlength_col = result[1]
		metagene_fiveprime_col = result[2]
		metagene_threeprime_col = result[3]
		a_col = result[4]
		t_col = result[5]
		g_col = result[6]
		c_col = result[7]
		title_size = result[8]
		subheading_size = result[9]
		axis_label_size = result[10]
		marker_size = result[11]

	if "readlen_ambig" in data:
		readlen_ambig = True
	else:
		readlen_ambig = False

	if "breakdown_per" in data:
		breakdown_per = True
	else:
		breakdown_per = False

	if "metagene_offsets" in data:
		metagene_offsets = True
	else:
		metagene_offsets = False

	if "metagene_aggregate" in data:
		metagene_aggregate = True
	else:
		metagene_aggregate = False

	if "metagene_normalise" in data:
		metagene_normalise = True
	else:
		metagene_normalise = False

	if "exclude_first" in data:
		exclude_first = True
	else:
		exclude_first = False

	if "exclude_last" in data:
		exclude_last = True
	else:
		exclude_last = False

	if "include_first" in data:
		include_first = True
	else:
		include_first = False

	if "include_last" in data:
		include_last = True
	else:
		include_last = False
	if "normalise" in data:
		normalise = True
	else:
		normalise = False



	if "log_scale" in data:
		log_scale = True
	else:
		log_scale = False

	if "mrna_dist_per" in data:
		mrna_dist_per = True
	else:
		mrna_dist_per = False

	if "md_start" in data:
		md_start = True
	else:
		md_start = False

	if "md_stop" in data:
		md_stop = True
	else:
		md_stop = False


	if "mrna_readlen_per" in data:
		mrna_readlen_per = True
	else:
		mrna_readlen_per = False
	if "count_agg" in data:
		count_agg = True
	else:
		count_agg = False

	if "mismatch_agg" in data:
		mismatch_agg = True
	else:
		mismatch_agg = False

	if "count_gc" in data:
		count_gc = True
	else:
		count_gc = False

	if "reverse_scale" in data:
		reverse_scale = True
	else:
		reverse_scale = False

	if html_args["user_short"] == "None" or user_short_passed == True:
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"traninfo_plot")
	else:
		short_code = html_args["user_short"]
		user_short_passed = True
	if plottype == "fetch_seq":
		filename = organism+"_sequences_"+str(time.time())+".fa"
		table_str = filename+"?~"
		splitlist = (gc_tranlist.replace(" ",",")).split(",")
		strlist = str(splitlist).strip("[]")
		tmp_fa_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
		connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
		connection.text_factory = str
		cursor = connection.cursor()
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		cursor = transhelve.cursor()
		
		
		cursor.execute("SELECT gene,transcript, sequence, cds_start, cds_stop FROM transcripts WHERE transcript IN ({})".format(strlist))

		#gc_location,exclude_first_val,exclude_last_val,include_first_val,include_last_val,gc_tranlist
		#keep track of transcript that have been written to file to avoid duplicates
		finished_trans = []
		result = cursor.fetchall()
		for row in result:
			gene = row[0]
			tran = row[1]
			if tran in finished_trans:
				continue
			seq = row[2]
			cds_start = row[3]
			cds_stop = row[4]
			if gc_location != "all" and cds_start == None:
				return "Error: {} is non-coding, region cannot be {}".format(tran, gc_location)
			if gc_location != "all":
				cds_start = int(cds_start)-1
				cds_stop = int(cds_stop)
			if gc_location == "five":
				seq = seq[:cds_start]
			elif gc_location == "cds":
				seq = seq[cds_start:cds_stop]
			elif gc_location == "three":
				seq = seq[cds_stop:]
			
			subseq = ""
			if sum([include_first_val, include_last_val,exclude_first_val,exclude_last_val]) == 0:
				subseq = seq
			else:
				if include_first_val != 0:
					subseq += seq[:include_first_val]
				if include_last_val != 0:
					subseq += seq[(include_last_val*-1):]
				if exclude_first_val != 0:
					seqlen = len(seq)
					subseq += seq[(seqlen-exclude_first_val)*-1:]
				if exclude_last_val != 0:
					seqlen = len(seq)
					subseq += seq[:(seqlen-exclude_last_val)]
			tmp_fa_file.write(">{}_{}\n{}\n".format(gene, tran, subseq))
			finished_trans.append(tran)
		tmp_fa_file.close()
		total_rows = len(finished_trans)
		table_str = "FA?~"+str(total_rows)+"?~"+table_str

		return table_str
	if plottype == "nuc_freq_plot":
		splitlist = (nuc_freq_plot_tranlist.replace(" ",",")).split(",")
		filename = "Sequences_{}.fa".format(time.time())
		outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
		connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
		connection.text_factory = str
		cursor = connection.cursor()
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		cursor = transhelve.cursor()
		if splitlist == []:
			cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts")
		else:
			cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE transcript IN ({})".format(str(splitlist).strip("[]")))
		result = cursor.fetchall()
		master_dict = {}
		for i in range(nuc_freq_plot_window*-1,nuc_freq_plot_window):
			master_dict[i] = {"A":0,"T":0,"G":0,"C":0}
		for row in result:
			tran = row[0]
			try:
				cds_start = int(row[1])
				cds_stop = int(row[2])
			except:
				cds_start = None
				cds_stop = None
			seq = row[3]
			seqlen = len(seq)
			if nuc_freq_plot_anchor == "tss":
				outfile.write(">{}\n{}\n".format(tran,seq[0:nuc_freq_plot_window]))
				for i in range(0,nuc_freq_plot_window):
					try:
						master_dict[i][seq[i]] += 1
					except:
						pass
			if nuc_freq_plot_anchor == "cds_start":
				if cds_start != None:
					outfile.write(">{}\n{}\n".format(tran,seq[cds_start+(nuc_freq_plot_window*-1):cds_start+nuc_freq_plot_window]))
					for i in range(nuc_freq_plot_window*-1,nuc_freq_plot_window):
						try:
							seqpos = cds_start+i
							if seqpos >= 0:
								master_dict[i][seq[seqpos]] += 1
						except:
							pass
			if nuc_freq_plot_anchor == "cds_stop":
				if cds_stop != None:
					outfile.write(">{}\n{}\n".format(tran,seq[cds_stop+(nuc_freq_plot_window*-1):cds_stop+nuc_freq_plot_window]))
					for i in range(nuc_freq_plot_window*-1,nuc_freq_plot_window):
						seqpos = cds_stop+i
						try:
							master_dict[i][seq[seqpos]] += 1
						except:
							pass
			if nuc_freq_plot_anchor == "tts":
				outfile.write(">{}\n{}\n".format(tran,seq[seqlen-nuc_freq_plot_window:seqlen]))
				for i in range(nuc_freq_plot_window*-1,0):
					try:
						master_dict[i][seq[seqlen+i]] += 1
					except:
						pass
		
		title = nuc_freq_plot_anchor
		return traninfo_plots.nuc_freq_plot(master_dict,title,short_code, background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size,filename)
		
	#Nucleotide composition (single transcript)
	if plottype == "nuc_comp_single":
		master_dict = {}
		metagene_tranlist = metagene_tranlist.split(",")
		if len(metagene_tranlist) == 1 and metagene_tranlist != ['']:
			tran = metagene_tranlist[0]
			connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
			connection.text_factory = str
			cursor = connection.cursor()
			cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
			owner = (cursor.fetchone())[0]
			if owner == 1:
				transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
			else:
				transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
			cursor = transhelve.cursor()
			cursor.execute("SELECT * from transcripts WHERE transcript = '{}'".format(tran.upper()))
			result = cursor.fetchone()
			traninfo = {"transcript":result[0] , "gene":result[1], "length":result[2], "cds_start":result[3] , "cds_stop":result[4] , "seq":result[5] ,
				"strand":result[6], "stop_list":result[7].split(","),"start_list":result[8].split(","), "exon_junctions":result[9].split(","),
				"tran_type":result[10], "principal":result[11]}

			title = tran
			connection.close()
			return traninfo_plots.nuc_comp_single(tran,master_dict,title,short_code, background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size,traninfo)
		else:
			connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
			connection.text_factory = str
			cursor = connection.cursor()
			cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
			owner = (cursor.fetchone())[0]
			if owner == 1:
				transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
			else:
				transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
			cursor = transhelve.cursor()
			if metagene_tranlist != ['']:
				cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE transcript IN ({})".format(str(metagene_tranlist).strip("[]").replace("'",'"')))
			else:
				cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE tran_type = 1 and principal = 1")
			result = cursor.fetchall()
			traninfo = []
			for row in result:
				traninfo.append([row[0],row[1],row[2],row[3]])
			title = "GC metagene of {} genes".format(len(traninfo))
			connection.close()
			return traninfo_plots.gc_metagene(title,short_code, background_col,readlength_col,title_size, axis_label_size, subheading_size,marker_size,traninfo)
		
	elif plottype == "orfstats":
		filename = organism+"_orfstats_"+str(time.time())+".csv"
		table_str = filename+"?~"
		tmp_te_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
		connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
		connection.text_factory = str
		cursor = connection.cursor()
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		cursor = transhelve.cursor()
		tmp_te_file.write("Gene,Tran,ORF Type, Start Codon, Length, Start position, Stop position, CDS coverage\n")
		if principal == True:
			cursor.execute("SELECT transcripts.gene, {0}.transcript, {0}.start_codon,{0}.length,{0}.start,{0}.stop, {0}.cds_coverage FROM {0} INNER JOIN transcripts ON transcripts.transcript = {0}.transcript WHERE transcripts.principal = 1".format(orftype))
		else:
			cursor.execute("SELECT transcripts.gene, {0}.transcript, {0}.start_codon,{0}.length,{0}.start,{0}.stop, {0}.cds_coverage FROM {0} INNER JOIN transcripts ON transcripts.transcript = {0}.transcript".format(orftype))
		input_list = []
		result = cursor.fetchall()
		for row in result:
			gene = row[0]
			tran = row[1]
			start_codon = row[2]
			length = row[3]
			start = row[4]
			stop = row[5]
			cds_cov = row[6]
			input_list.append([gene,tran,orftype,start_codon,length,start,stop,cds_cov])
			tmp_te_file.write("{},{},{},{},{},{},{},{}\n".format(gene,tran,orftype,start_codon,length,start,stop,cds_cov))
		tmp_te_file.close()
		total_rows = 0
		#(transcript VARCHAR(300), start_codon VARCHAR(10), length INT(6), start INT(6), stop INT(6), sequence VARCHAR(50000), cds_coverage FLOAT(20))
		#input_list = [["SRRBLAH", "ATP12","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"],["SRRBLAH", "ATP13","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"],["SRRBLAH", "ATP14","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"],["SRRBLAH", "ATP15","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"],["SRRBLAH", "ATP16","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"],["SRRBLAH", "ATP12","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"],["SRRBLAH", "ATP12","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"],["SRRBLAH", "ATP12","ENST00004483939",orftype, "0", "0", "0","https://trips.ucc.ie"]]	
		all_sorted_rows = input_list
		for row in all_sorted_rows:
			total_rows += 1
			if total_rows <= 1000:
				input_str = ""
				for item in row:
					input_str += "{}.;".format(item)
				input_str += "?~"
				table_str += input_str
		table_str = "TE?~"+str(total_rows)+"?~"+table_str
		return table_str
		
		
		
		return "ORF STATS"

	#Nucleotide composition (multiple transcripts)
	elif plottype == "nuc_comp_multi":
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		filename = "{}_{}_nucleotide_comp_{}.csv".format(organism,gc_location,str(time.time()))
		outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
		trancursor = transhelve.cursor()
		
		if plot_type == "scatter":
			master_dict =  {1:{"trans":[],"lengths":[]},
							2:{"trans":[],"lengths":[]},
							3:{"trans":[],"lengths":[]},
							4:{"trans":[],"lengths":[]}}
			gc_dict = {}
			if gc_tranlist != "":
				gc_tranlist = gc_tranlist.upper()
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				
				if gc_tranlist2 != "":
					splitlist2 = (gc_tranlist2.replace(" ",",")).split(",")
					for item in splitlist2:
						splitlist.append(item)
				if gc_tranlist3 != "":
					splitlist3 = (gc_tranlist3.replace(" ",",")).split(",")
					for item in splitlist3:
						splitlist.append(item)
				if gc_tranlist4 != "":
					splitlist4 = (gc_tranlist4.replace(" ",",")).split(",")
					for item in splitlist4:
						splitlist.append(item)
				strlist = str(splitlist).strip("[]")
				trancursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE transcript IN ({});".format(strlist.upper()))
			else:
				if gc_location == "all":
					trancursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE principal = 1;")
				else:
					trancursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE principal = 1 and tran_type = 1;")
			result = trancursor.fetchall()
			for row in result:
				tran = row[0]
				seq = row[1]
				if gc_location != "all":
					cds_start = int(row[2])
					cds_stop = int(row[3])
				if gc_location == "all":
					nuc_dict = calc_gc(seq)
					#outfile.write(">{}\n{}\n".format(tran,seq))
				elif gc_location == "five":
					nuc_dict = calc_gc(seq[:cds_start-1])
					#outfile.write("{>{}\n{}\n".format(tran,seq[:cds_start-1]))
				elif gc_location == "cds":
					nuc_dict = calc_gc(seq[cds_start-1:cds_stop+3])
					#outfile.write(">{}\n{}\n".format(tran,seq[cds_start-1:cds_stop+3]))
				elif gc_location == "three":
					nuc_dict = calc_gc(seq[cds_stop+3:])
					#outfile.write(">{}\n{}\n".format(tran,seq[cds_stop+3:]))
				tot = sum([nuc_dict["G"],nuc_dict["C"],nuc_dict["T"],nuc_dict["A"]])
				if tot > 0:
					if nucleotide == "GC":
						gc = ((nuc_dict["G"]+nuc_dict["C"])/tot)
						gc = gc*100
						gc_dict[tran]= gc
					elif nucleotide == "A":
						a = ((nuc_dict["A"])/tot)
						a = a*100
						gc_dict[tran]= a
					elif nucleotide == "T":
						t = ((nuc_dict["T"])/tot)
						t = t*100
						gc_dict[tran]= t
					elif nucleotide == "G":
						g = ((nuc_dict["G"])/tot)
						g = g*100
						gc_dict[tran]= g
					elif nucleotide == "C":
						c = ((nuc_dict["C"])/tot)
						c = c*100
						gc_dict[tran]= c
			if gc_tranlist == "":
				for tran in gc_dict:
					master_dict[1]["trans"].append(tran)
					master_dict[1]["lengths"].append(gc_dict[tran])
			else:
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				for tran in splitlist:
					if tran in gc_dict:
						master_dict[1]["trans"].append(tran)
						master_dict[1]["lengths"].append(gc_dict[tran])
				if gc_tranlist2 != "":
					splitlist = (gc_tranlist2.replace(" ",",")).split(",")
					for tran in splitlist2:
						if tran in gc_dict:
							master_dict[2]["trans"].append(tran)
							master_dict[2]["lengths"].append(gc_dict[tran])
				if gc_tranlist3 != "":
					splitlist = (gc_tranlist3.replace(" ",",")).split(",")
					for tran in splitlist3:
						if tran in gc_dict:
							master_dict[3]["trans"].append(tran)
							master_dict[3]["lengths"].append(gc_dict[tran])
				if gc_tranlist4 != "":
					splitlist = (gc_tranlist4.replace(" ",",")).split(",")
					for tran in splitlist4:
						if tran in gc_dict:
							master_dict[4]["trans"].append(tran)
							master_dict[4]["lengths"].append(gc_dict[tran])
			
			outfile.write("Group , GC\n")
			for group in master_dict:
				for gc in master_dict[group]["lengths"]:
					outfile.write("{}.{}\n".format(group, gc))
			
			outfile.close()
			return traninfo_plots.nuc_comp_scatter(master_dict,filename,str(title_size)+"pt",str(axis_label_size)+"pt",str(marker_size)+"pt",nucleotide,short_code)
		elif plot_type == "box":
			master_dict =  {1:{"trans":[],"gc":[],"a":[],"t":[],"g":[],"c":[]},
							2:{"trans":[],"gc":[],"a":[],"t":[],"g":[],"c":[]},
							3:{"trans":[],"gc":[],"a":[],"t":[],"g":[],"c":[]},
							4:{"trans":[],"gc":[],"a":[],"t":[],"g":[],"c":[]}}
			gc_dict = {}
			if gc_tranlist != "":
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				
				if gc_tranlist2 != "":
					splitlist2 = (gc_tranlist2.replace(" ",",")).split(",")
					for item in splitlist2:
						splitlist.append(item)
				if gc_tranlist3 != "":
					splitlist3 = (gc_tranlist3.replace(" ",",")).split(",")
					for item in splitlist3:
						splitlist.append(item)
				if gc_tranlist4 != "":
					splitlist4 = (gc_tranlist4.replace(" ",",")).split(",")
					for item in splitlist4:
						splitlist.append(item)
				strlist = str(splitlist).strip("[]")
				trancursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE transcript IN ({});".format(strlist.upper()))
			else:
				if gc_location == "all":
					trancursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE principal = 1;")
				else:
					trancursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE principal = 1 and tran_type = 1;")
			result = trancursor.fetchall()
			for row in result:
				tran = row[0]
				seq = row[1]
				if gc_location != "all":
					try:
						cds_start = int(row[2])
						cds_stop = int(row[3])
					except:
						return "Error: {} is non coding, remove this from the transcript list or change the region to 'All'".format(tran, row[2])
				if gc_location == "all":
					nuc_dict = calc_gc(seq)
					#outfile.write(">{}\n{}\n".format(tran,seq))
				elif gc_location == "five":
					nuc_dict = calc_gc(seq[:cds_start-1])
					outfile.write(">{}\n{}\n".format(tran,seq[:cds_start-1]))
				elif gc_location == "cds":
					nuc_dict = calc_gc(seq[cds_start-1:cds_stop+3])
					#outfile.write(">{}\n{}\n".format(tran,seq[cds_start-1:cds_stop+3]))
				elif gc_location == "three":
					nuc_dict = calc_gc(seq[cds_stop+3:])
					#outfile.write(">{}\n{}\n".format(tran,seq[cds_stop+3:]))
				tot = sum([nuc_dict["G"],nuc_dict["C"],nuc_dict["T"],nuc_dict["A"]])
				if tot > 0:
					gc = ((nuc_dict["G"]+nuc_dict["C"])/tot)
					a = ((nuc_dict["A"])/tot)
					t = ((nuc_dict["T"])/tot)
					g = ((nuc_dict["G"])/tot)
					c = ((nuc_dict["C"])/tot)
					
					gc = gc*100
					a = a*100
					t = t*100
					g = g*100
					c = c*100
					
					gc_dict[tran]= {"gc":gc,"a":a,"t":t,"g":g,"c":c}
			if gc_tranlist == "":
				for tran in gc_dict:
					master_dict[1]["trans"].append(tran)
					master_dict[1]["gc"].append(gc_dict[tran]["gc"])
					master_dict[1]["a"].append(gc_dict[tran]["a"])
					master_dict[1]["t"].append(gc_dict[tran]["t"])
					master_dict[1]["g"].append(gc_dict[tran]["g"])
					master_dict[1]["c"].append(gc_dict[tran]["c"])
			else:
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				for tran in splitlist:
					if tran in gc_dict:
						master_dict[1]["trans"].append(tran)
						master_dict[1]["gc"].append(gc_dict[tran]["gc"])
						master_dict[1]["a"].append(gc_dict[tran]["a"])
						master_dict[1]["t"].append(gc_dict[tran]["t"])
						master_dict[1]["g"].append(gc_dict[tran]["g"])
						master_dict[1]["c"].append(gc_dict[tran]["c"])
				if gc_tranlist2 != "":
					splitlist = (gc_tranlist2.replace(" ",",")).split(",")
					for tran in splitlist2:
						if tran in gc_dict:
							master_dict[2]["trans"].append(tran)
							master_dict[2]["gc"].append(gc_dict[tran]["gc"])
							master_dict[2]["a"].append(gc_dict[tran]["a"])
							master_dict[2]["t"].append(gc_dict[tran]["t"])
							master_dict[2]["g"].append(gc_dict[tran]["g"])
							master_dict[2]["c"].append(gc_dict[tran]["c"])
				if gc_tranlist3 != "":
					splitlist = (gc_tranlist3.replace(" ",",")).split(",")
					for tran in splitlist3:
						if tran in gc_dict:
							master_dict[3]["trans"].append(tran)
							master_dict[3]["gc"].append(gc_dict[tran]["gc"])
							master_dict[3]["a"].append(gc_dict[tran]["a"])
							master_dict[3]["t"].append(gc_dict[tran]["t"])
							master_dict[3]["g"].append(gc_dict[tran]["g"])
							master_dict[3]["c"].append(gc_dict[tran]["c"])
				if gc_tranlist4 != "":
					splitlist = (gc_tranlist4.replace(" ",",")).split(",")
					for tran in splitlist4:
						if tran in gc_dict:
							master_dict[4]["trans"].append(tran)
							master_dict[4]["gc"].append(gc_dict[tran]["gc"])
							master_dict[4]["a"].append(gc_dict[tran]["a"])
							master_dict[4]["t"].append(gc_dict[tran]["t"])
							master_dict[4]["g"].append(gc_dict[tran]["g"])
							master_dict[4]["c"].append(gc_dict[tran]["c"])
			outfile.write("Group ,GC\n")
			for group in master_dict:
				for gc in master_dict[group]["gc"]:
					outfile.write("{},{}\n".format(group, gc))				
			outfile.close()
			return traninfo_plots.nuc_comp_box(master_dict,filename,nucleotide,str(title_size)+"pt",config.A_COL,str(axis_label_size)+"pt",str(marker_size)+"pt",short_code)
			
		#return traninfo_plots.nuc_comp(master_dict, nuc_maxreadlen,title, nuc_comp_type,nuc_comp_direction,short_code,background_col,a_col,t_col,g_col,c_col,title_size, axis_label_size, subheading_size,marker_size)

	elif plottype == "codon_usage":
		aa_dict = { "TTT":"Phenylalanine", "TTC":"Phenylalanine", "TTA":"Leucine", "TTG":"Leucine",
					"TCT":"Serine", "TCC":"Serine", "TCA":"Serine", "TCG":"Serine",
					"TAT":"Tyrosine", "TAC":"Tyrosine", "TAA":"*", "TAG":"*",
					"TGT":"Cysteine", "TGC":"Cysteine", "TGA":"*", "TGG":"Tryptophan",
					"CTT":"Leucine", "CTC":"Leucine", "CTA":"Leucine", "CTG":"Leucine",

					"CCT":"Proline", "CCC":"Proline", "CCA":"Proline", "CCG":"Proline",
					"CAT":"Histidine", "CAC":"Histidine", "CAA":"Glutamine", "CAG":"Glutamine",
					"CGT":"Arginine", "CGC":"Arginine", "CGA":"Arginine", "CGG":"Arginine",
					"ATT":"Isoleucine", "ATC":"Isoleucine", "ATA":"Isoleucine", "ATG":"Methionine",

					"ACT":"Threonine", "ACC":"Threonine", "ACA":"Threonine", "ACG":"Threonine",
					"AAT":"Asparagine", "AAC":"Asparagine", "AAA":"Lysine", "AAG":"Lysine",
					"AGT":"Serine", "AGC":"Serine", "AGA":"Arginine", "AGG":"Arginine",
					"GTT":"Valine", "GTC":"Valine", "GTA":"Valine", "GTG":"Valine",

					"GCT":"Alanine", "GCC":"Alanine", "GCA":"Alanine", "GCG":"Alanine",
					"GAT":"Aspartic Acid", "GAC":"Aspartic Acid", "GAA":"Glutamic Acid", "GAG":"Glutamic Acid",
					"GGT":"Glycine", "GGC":"Glycine", "GGA":"Glycine", "GGG":"Glycine"}
		filename = organism+"_codon_usage_"+str(time.time())+".csv"
		cu_file = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))


		traninfo_cursor = transhelve.cursor()
		codon_dict = {}
		if gc_tranlist != "":
			splitlist = (gc_tranlist.replace(" ",",")).split(",")
			strlist = str(splitlist).strip("[]")
			traninfo_cursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE transcript IN ({});".format(strlist.upper()))
		else:
			traninfo_cursor.execute("SELECT transcript,sequence,cds_start,cds_stop from transcripts WHERE principal = 1 and tran_type = 1;")
		result = traninfo_cursor.fetchall()
		total_trans = 0
		total_length = 0
		for row in result:
			tran = row[0]
			seq = row[1]
			total_trans += 1
			if row[2] == None or row[3] == None:
				return "Enter coding transcripts only, {} is noncoding".format(row[0])
			start_pos = int(row[2])-1
			end_pos = int(row[3])+2
			total_length += (end_pos-start_pos)+1
			for i in range(start_pos,end_pos,3):
				try:
					codon = seq[i:i+3]
				except:
					pass
				if len(codon) != 3:
					continue
				if codon not in codon_dict:
					codon_dict[codon] = 0
				codon_dict[codon] += 1
		cu_file.write("Total_transcripts,{}\n".format(total_trans))
		cu_file.write("Total_length_(nts),{}\n".format(total_length))
		for codon in sorted(codon_dict.keys()):
			cu_file.write("{},{},{}\n".format(codon,aa_dict[codon], codon_dict[codon]))
				
		return traninfo_plots.codon_usage(codon_dict, short_code, str(title_size)+"pt", str(axis_label_size)+"pt", str(marker_size)+"pt",filename)
	#This is the lengths plot
	elif plottype == "lengths_plot":
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		
		trancursor = transhelve.cursor()
		
		if plot_type == "box":
			master_dict =  {1:{"trans":[],"lengths":[]},
							2:{"trans":[],"lengths":[]},
							3:{"trans":[],"lengths":[]},
							4:{"trans":[],"lengths":[]}}
			gc_dict = {}
			if gc_tranlist != "":
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				
				if gc_tranlist2 != "":
					splitlist2 = (gc_tranlist2.replace(" ",",")).split(",")
					for item in splitlist2:
						splitlist.append(item)
				if gc_tranlist3 != "":
					splitlist3 = (gc_tranlist3.replace(" ",",")).split(",")
					for item in splitlist3:
						splitlist.append(item)
				if gc_tranlist4 != "":
					splitlist4 = (gc_tranlist4.replace(" ",",")).split(",")
					for item in splitlist4:
						splitlist.append(item)
				strlist = str(splitlist).strip("[]")
				trancursor.execute("SELECT transcript,length,cds_start,cds_stop,exon_junctions from transcripts WHERE transcript IN ({});".format(strlist.upper()))
			else:
				if gc_location == "all":
					trancursor.execute("SELECT transcript,length,cds_start,cds_stop,exon_junctions from transcripts WHERE principal = 1;")
				else:
					trancursor.execute("SELECT transcript,length,cds_start,cds_stop,exon_junctions from transcripts WHERE principal = 1 and tran_type = 1;")
			result = trancursor.fetchall()
			if result == []:
				return "Could not find any info on given transcript list"
			for row in result:
				tran = row[0]
				tranlen = int(row[1])
				exon_junctions_raw = row[4].split(",")
				exon_junctions = []
				for item in exon_junctions_raw:
					if item != '':
						exon_junctions.append(int(item))
				if gc_location != "all":
					cds_start = int(row[2])
					cds_stop = int(row[3])
					
				if gc_location == "all":
					length = tranlen
					exon_no = length/(len(exon_junctions)+1)
				elif gc_location == "five":
					length = cds_start
					local_exons = sum(i < cds_start for i in exon_junctions)
					exon_no = length/(local_exons+1)
				elif gc_location == "cds":
					length = cds_stop-cds_start
					local_exons = sum((i > cds_start and i < cds_stop) for i in exon_junctions)
					exon_no = length/(local_exons+1)
				elif gc_location == "three":
					length = tranlen-cds_stop
					local_exons = sum(i > cds_stop for i in exon_junctions)
					exon_no = length/(local_exons+1)
				if exons == False:
					gc_dict[tran]= length
				else:
					gc_dict[tran] = exon_no
			filename = "Lengths_{}.csv".format(time.time())
			outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")	
			if gc_tranlist == "":
				for tran in gc_dict:
					master_dict[1]["trans"].append(tran)
					master_dict[1]["lengths"].append(gc_dict[tran])
			else:
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				for tran in splitlist:
					if tran in gc_dict:
						master_dict[1]["trans"].append(tran)
						master_dict[1]["lengths"].append(gc_dict[tran])
						outfile.write("Group1,{},{}\n".format(tran,gc_dict[tran]))
					else:
						return tran
				if gc_tranlist2 != "":
					splitlist = (gc_tranlist2.replace(" ",",")).split(",")
					for tran in splitlist2:
						if tran in gc_dict:
							master_dict[2]["trans"].append(tran)
							master_dict[2]["lengths"].append(gc_dict[tran])
							outfile.write("Group2,{},{}\n".format(tran,gc_dict[tran]))
				if gc_tranlist3 != "":
					splitlist = (gc_tranlist3.replace(" ",",")).split(",")
					for tran in splitlist3:
						if tran in gc_dict:
							master_dict[3]["trans"].append(tran)
							master_dict[3]["lengths"].append(gc_dict[tran])
							outfile.write("Group3,{},{}\n".format(tran,gc_dict[tran]))
				if gc_tranlist4 != "":
					splitlist = (gc_tranlist4.replace(" ",",")).split(",")
					for tran in splitlist4:
						if tran in gc_dict:
							master_dict[4]["trans"].append(tran)
							master_dict[4]["lengths"].append(gc_dict[tran])
							outfile.write("Group4,{},{}\n".format(tran,gc_dict[tran]))
			return traninfo_plots.lengths_box(master_dict,filename,config.A_COL,short_code,str(title_size)+"pt",str(marker_size)+"pt",str(axis_label_size)+"pt")
		elif plot_type == "scatter":
			master_dict =  {1:{"trans":[],"lengths":[]},
							2:{"trans":[],"lengths":[]},
							3:{"trans":[],"lengths":[]},
							4:{"trans":[],"lengths":[]}}
			gc_dict = {}
			if gc_tranlist != "":
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				
				if gc_tranlist2 != "":
					splitlist2 = (gc_tranlist2.replace(" ",",")).split(",")
					for item in splitlist2:
						splitlist.append(item)
				if gc_tranlist3 != "":
					splitlist3 = (gc_tranlist3.replace(" ",",")).split(",")
					for item in splitlist3:
						splitlist.append(item)
				if gc_tranlist4 != "":
					splitlist4 = (gc_tranlist4.replace(" ",",")).split(",")
					for item in splitlist4:
						splitlist.append(item)
				strlist = str(splitlist).strip("[]")
				trancursor.execute("SELECT transcript,length,cds_start,cds_stop from transcripts WHERE transcript IN ({});".format(strlist.upper()))
			else:
				if gc_location == "all":
					trancursor.execute("SELECT transcript,length,cds_start,cds_stop from transcripts WHERE principal = 1;")
				else:
					trancursor.execute("SELECT transcript,length,cds_start,cds_stop from transcripts WHERE principal = 1 and tran_type = 1;")
			result = trancursor.fetchall()
			for row in result:
				tran = row[0]
				tranlen = int(row[1])
				if gc_location != "all":
					cds_start = int(row[2])
					cds_stop = int(row[3])
				if gc_location == "all":
					length = tranlen
				elif gc_location == "five":
					length = cds_start
				elif gc_location == "cds":
					length = cds_stop-cds_start
				elif gc_location == "three":
					length = tranlen-cds_stop
				gc_dict[tran]= length
			filename = "Lengths_{}.csv".format(time.time())
			outfile = open("{}/static/tmp/{}".format(config.SCRIPT_LOC,filename),"w")	
						
			if gc_tranlist == "":
				for tran in gc_dict:
					master_dict[1]["trans"].append(tran)
					master_dict[1]["lengths"].append(gc_dict[tran])
			else:
				splitlist = (gc_tranlist.replace(" ",",")).split(",")
				for tran in splitlist:
					if tran in gc_dict:
						master_dict[1]["trans"].append(tran)
						master_dict[1]["lengths"].append(gc_dict[tran])
				if gc_tranlist2 != "":
					splitlist = (gc_tranlist2.replace(" ",",")).split(",")
					for tran in splitlist2:
						if tran in gc_dict:
							master_dict[2]["trans"].append(tran)
							master_dict[2]["lengths"].append(gc_dict[tran])
				if gc_tranlist3 != "":
					splitlist = (gc_tranlist3.replace(" ",",")).split(",")
					for tran in splitlist3:
						if tran in gc_dict:
							master_dict[3]["trans"].append(tran)
							master_dict[3]["lengths"].append(gc_dict[tran])
				if gc_tranlist4 != "":
					splitlist = (gc_tranlist4.replace(" ",",")).split(",")
					for tran in splitlist4:
						if tran in gc_dict:
							master_dict[4]["trans"].append(tran)
							master_dict[4]["lengths"].append(gc_dict[tran])
			outfile.close()
			return traninfo_plots.lengths_scatter(master_dict,filename,str(title_size)+"pt",str(axis_label_size)+"pt",str(marker_size)+"pt",short_code)

	elif plottype == "gene_count":
		cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
		owner = (cursor.fetchone())[0]
		if owner == 1:
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
		cursor = transhelve.cursor()
		cursor.execute("SELECT COUNT(*) FROM transcripts;")
		all_transcripts = cursor.fetchall()[0][0]
		cursor.execute("SELECT COUNT(*) FROM transcripts where tran_type = 1;")
		coding_transcripts = cursor.fetchall()[0][0]
		cursor.execute("SELECT COUNT(DISTINCT gene) FROM transcripts;")
		all_genes = cursor.fetchall()[0][0]
		cursor.execute("SELECT COUNT(DISTINCT gene) FROM transcripts WHERE gene_type = 1")
		coding_genes = cursor.fetchall()[0][0]
		coding = [1,coding_genes,coding_transcripts,1]
		noncoding = [1,all_genes-coding_genes,all_transcripts-coding_transcripts,1]
		return traninfo_plots.gene_count(short_code,background_col,title_size, axis_label_size, subheading_size,marker_size,coding,noncoding)


	else:
		if (plottype.strip(" ").replace("\n","")) not in ["replicate_comp"]:
			print ("Unknown plottype",plottype)
	return "Error, unknown plot type selected: {}".format(plottype)


def get_nuc_comp_reads(sqlite_db, nuccomp_reads, organism, transcriptome):
	cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (cursor.fetchone())[0]
	if owner == 1:
		transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
	else:
		transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
	cursor = transhelve.cursor()
	cursor.execute("SELECT transcript,cds_start,cds_stop,sequence from transcripts WHERE principal = 1")
	result = cursor.fetchall()
	master_dict = {}
	offsets = sqlite_db["offsets"]["fiveprime"]["offsets"]
	for row in result:
		tran = row[0]
		cds_start = int(row[1])
		cds_stop = int(row[2])
		seq = row[3].replace("T","U")
		if tran in sqlite_db:
			counts = sqlite_db[tran]["unambig"]
			for readlen in counts:
				if readlen not in master_dict:
					master_dict[readlen] = {}
					for i in range(0,int(readlen)):
						master_dict[readlen][i] = {"A":0,"T":0,"G":0,"C":0,"N":0}
				for pos in counts[readlen]:
					count = counts[readlen][pos]
					if readlen in offsets:
						offset_pos = pos+offsets[readlen]
					else:
						offset_pos = pos+15
					if offset_pos >= cds_start and offset_pos <= cds_stop:
						readframe = offset_pos%3
						cds_frame = (cds_start+2)%3
						if readframe == cds_frame:
							inframe = True
						else:
							inframe = False
						if nuccomp_reads == "inframe" and inframe == False:
							continue
						if nuccomp_reads == "offrame" and inframe == True:
							continue
						readseq = seq[pos:pos+readlen]
						for i in range(0,len(readseq)):
							char = readseq[i].replace("U","T")
							master_dict[readlen][i][char]+= count
	#save results so they won't have to be computed again later
	if "nuc_counts" in sqlite_db:
		new_nuc_counts = sqlite_db["nuc_counts"]
		new_nuc_counts[nuccomp_reads] = master_dict
		sqlite_db["nuc_counts"] = new_nuc_counts
	else:
		sqlite_db["nuc_counts"] = {nuccomp_reads:master_dict}


	sqlite_db.commit()
	return master_dict






def calc_gc(seq):
	count_dict = {"A":0.0,"T":0.0,"G":0.0,"C":0.0,"N":0.0}
	for char in seq:
		count_dict[char] += 1
	return count_dict
