from flask import Blueprint, render_template, abort, request, jsonify, url_for
from flask import current_app as app
import sqlite3
from sqlitedict import SqliteDict
import ast
import os
import config
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code,fetch_user
import riboflask_compare
import collections
from flask_login import current_user
import json



# Single transcript comparison page, user chooses a gene and groups of files to display
comparison_plotpage_blueprint = Blueprint("comparisonpage", __name__, template_folder="templates")
@comparison_plotpage_blueprint.route('/<organism>/<transcriptome>/comparison/')
def comparisonpage(organism, transcriptome):
	#global user_short_passed
	user_short_passed = False
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
	user_file_dict = {}
	if str(request.args.get('files')) != "None":
		colors = str(request.args.get('files')).split("_")
		for filelist in colors:
			all_items = filelist.split(",")
			files = []
			for item in all_items:
				if "#" in item:
					color = item
				else:
					files.append(item)
			user_file_dict[color] = files

	user_label_dict = {}
	if str(request.args.get('labels')) != "None":
		colors = str(request.args.get('labels')).split("_")
		for label in colors:
			all_items = label.split(",")
			for item in all_items:
				if "#" in item:
					color = item
				else:
					label = item
			user_label_dict[color] = label

	html_args = {"user_short":str(request.args.get('short')),
				 "user_file_dict":user_file_dict,
				 "user_label_dict":user_label_dict,
				 "transcript":str(request.args.get('transcript')),
				 "minread":str(request.args.get('minread')),
				 "maxread":str(request.args.get('maxread')),
				 "hili_start":str(request.args.get('hili_start')),
				 "hili_stop":str(request.args.get('hili_stop')),
				 "ambig":str(request.args.get('ambig')),
				 "cov":str(request.args.get('cov')),
				 "normalize":str(request.args.get('normalize')),
				 "transcriptome":str(transcriptome)}

	accepted_studies = fetch_studies(user, organism, transcriptome)
	file_id_to_name_dict,accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)
	print (local)
	connection.close()
	return render_template('index_compare.html', studies_dict=accepted_studies, accepted_files=accepted_files, gwips_info=gwips_info,
						gwips_clade=gwips_clade, gwips_org=gwips_org,gwips_db=gwips_db, organism=organism,transcriptome=transcriptome,
						default_tran=default_tran,local=local,html_args=html_args,file_id_to_name_dict=file_id_to_name_dict,studyinfo_dict=studyinfo_dict,
						seq_types=seq_types)




# Creates/serves the comparison plots
comparisonquery_blueprint = Blueprint("comparequery", __name__, template_folder="templates")
@comparisonquery_blueprint.route('/comparequery', methods=['POST'])
def comparequery():
	#global user_short_passed
	user_short_passed = False
	tran_dict = {}
	data = json.loads(request.data)
	tran = data['transcript'].upper().strip()
	organism = data['organism']
	transcriptome = data['transcriptome']
	
	user,logged_in = fetch_user()
	
	trips_connection = sqlite3.connect('{}/{}'.format(config.SCRIPT_LOC,config.DATABASE_NAME))
	trips_connection.text_factory = str
	trips_cursor = trips_connection.cursor()
	
	trips_cursor.execute("SELECT owner FROM organisms WHERE organism_name = '{}' and transcriptome_list = '{}';".format(organism, transcriptome))
	owner = (trips_cursor.fetchone())[0]
	if owner == 1:
		if os.path.isfile("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome)):
			transhelve = sqlite3.connect("{0}/{1}/{2}/{2}.{3}.sqlite".format(config.SCRIPT_LOC, config.ANNOTATION_DIR,organism,transcriptome))
		else:
			return "Cannot find annotation file {}.{}.sqlite".format(organism,transcriptome)
	else:
		transhelve = sqlite3.connect("{0}transcriptomes/{1}/{2}/{3}/{2}_{3}.sqlite".format(config.UPLOADS_DIR,owner,organism,transcriptome))
	cursor = transhelve.cursor()
	cursor.execute("SELECT * from transcripts WHERE transcript = '{}'".format(tran))
	result = cursor.fetchone()
	inputtran = True
	if result != None:
		newtran = result[0]
	else:
		inputtran = False
	if inputtran == False:
		cursor.execute("SELECT * from transcripts WHERE gene = '{}'".format(tran))
		result = cursor.fetchall()
		if result != []:
			if len(result) == 1:
				tran = str(result[0][0])
			else:
				return_str = "TRANSCRIPTS"
				for transcript in result:
					cursor.execute("SELECT length,cds_start,cds_stop,principal from transcripts WHERE transcript = '{}'".format(transcript[0]))
					tran_result = cursor.fetchone()
					tranlen = tran_result[0]
					cds_start = tran_result[1]
					cds_stop = tran_result[2]
					if tran_result[3] == 1:
						principal = "principal"
					else:
						principal = ""
					if cds_start == "NULL" or cds_start == None:
						cdslen = "NULL"
						threeutrlen = "NULL"
					else:
						cdslen = cds_stop-cds_start
						threeutrlen = tranlen - cds_stop
					return_str += (":{},{},{},{},{},{}".format(transcript[0], tranlen, cds_start, cdslen, threeutrlen,principal))
				tot_prog = 100
				if app.debug == True:
					return return_str, "NO_CELERY", {'Location': None}
				else:
					return jsonify({'current': 100, 'total': 100, 'status': 'tran_list','result': return_str}), 200, {'Location': ""}  
		else:
			return_str =  "ERROR! Could not find any transcript corresponding to {}".format(tran)
			if app.debug == True:
				return return_str, "NO_CELERY", {'Location': None}
			else:
				return jsonify({'current': 100, 'total': 100, 'status': 'return_str','result': return_str}), 200, {'Location': ""} 
	transhelve.close()
	minread = int(data['minread'])
	maxread = int(data['maxread'])
	hili_start = int(data['hili_start'])
	hili_stop = int(data['hili_stop'])
	master_filepath_dict = {}
	master_file_dict = data['master_file_dict']
	
	if master_file_dict == {}:
		return_str =  "Error: No files in the File list box. To add files to the file list box click on a study in the studies section above. This will populate the Ribo-seq and RNA-Seq sections with a list of files. Click on one of the files and then press the  Add button in the studies section. This will add the file to the File list box. Selecting another file and clicking Add again will add the new file to the same group in the File list. Alternatively to add a new group simply change the selected colour (by clicking on the coloured box in the studies section) and then click the Add file button."
		if app.debug == True:
			return return_str, "NO_CELERY", {'Location': None}
		else:
			return jsonify({'current': 100, 'total': 100, 'status': 'return_str','result': return_str}), 200, {'Location': ""}  


	# This section is purely to sort by label alphabetically
	for color in master_file_dict:
		master_filepath_dict[color] = {"filepaths":[],"file_ids":[],"file_names":[],"file_descs":[],"mapped_reads":0,"minread":minread,"maxread":maxread}
		# Overwrite the default minread and maxread with the minread/maxread values that are group specific, this allows users to easily visualise
		# how the profile of different readlenghts differs across a transcript
		if "minread" in master_file_dict[color]:
			master_filepath_dict[color]["minread"] = int(master_file_dict[color]["minread"])
			
		if "maxread" in master_file_dict[color]:
			master_filepath_dict[color]["maxread"] = int(master_file_dict[color]["maxread"])
			
		for file_id in master_file_dict[color]["file_ids"]:
			trips_cursor.execute("SELECT file_name,file_description,file_type from files WHERE file_id = {};".format(file_id))
			result = (trips_cursor.fetchone())
			file_name = master_file_dict[color]["label"] 
			file_paths = fetch_file_paths([file_id],organism)

			for filetype in file_paths:
				for file_id in file_paths[filetype]:
					filepath = file_paths[filetype][file_id]
					if os.path.isfile(filepath):
						sqlite_db = SqliteDict(filepath, autocommit=False)
					else:
						return_str =  "File not found, please report this to tripsvizsite@gmail.com or via the contact page."
						if app.debug == True:
							return return_str, "NO_CELERY", {'Location': None}
						else:
							return jsonify({'current': 100, 'total': 100, 'status': 'return_str','result': return_str}), 200, {'Location': ""} 
					#if maxread != 150:
					#	read_lengths = sqlite_db["read_lengths"]
					#	for i in range(master_filepath_dict[color]["minread"],master_filepath_dict[color]["maxread"]+1):
					#		master_filepath_dict[color]["mapped_reads"] += read_lengths[i]
					
					if "noncoding_counts" in sqlite_db and "coding_counts" in sqlite_db:
						master_filepath_dict[color]["mapped_reads"] += float(sqlite_db["noncoding_counts"])
						master_filepath_dict[color]["mapped_reads"] += float(sqlite_db["coding_counts"])
					else:
						if "normalize" in data:
							return_str = "One or more selected files is missing values for 'coding_counts' and 'non_coding_counts' so cannot normalize with these files, please report this to tripsvizsite@gmail.com or via the contact page."
							if app.debug == True:
								return return_str, "NO_CELERY", {'Location': None}
							else:
								return jsonify({'current': 100, 'total': 100, 'status': 'return_str','result': return_str}), 200, {'Location': ""}  
					master_filepath_dict[color]["filepaths"].append(filepath)
					master_filepath_dict[color]["file_ids"].append(file_id)
					master_filepath_dict[color]["file_names"].append(file_name)
					master_filepath_dict[color]["file_descs"].append(result[1])
					master_filepath_dict[color]["file_type"] = result[2]

	if 'ribocoverage' in data:
		ribocoverage = True
	else:
		ribocoverage = False
	if "ambiguous" in data:
		ambiguous = "ambig"
	else:
		ambiguous = "unambig"
		
	if "normalize" in data:
		normalize = True
	else:
		normalize = False
	html_args = data["html_args"]
	if html_args["user_short"] == "None" or user_short_passed == True:
		short_code = generate_short_code(data,organism,html_args["transcriptome"],"comparison")
	else:
		short_code = html_args["user_short"]
		user_short_passed = True


	#set colours to default values, if user logged in these will be overwritten
	background_col = config.BACKGROUND_COL
	comp_uga_col = config.UGA_COL
	comp_uag_col = config.UAG_COL
	comp_uaa_col = config.UAA_COL
	title_size = config.TITLE_SIZE
	subheading_size = config.SUBHEADING_SIZE
	axis_label_size = config.AXIS_LABEL_SIZE
	marker_size = config.MARKER_SIZE
	cds_marker_size = config.CDS_MARKER_SIZE
	cds_marker_colour = config.CDS_MARKER_COLOUR
	legend_size = config.LEGEND_SIZE
	if current_user.is_authenticated:
		user_name = current_user.name
		trips_cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user_name))
		result = (trips_cursor.fetchone())
		user_id = result[0]
		#get a list of organism id's this user can access
		trips_cursor.execute("SELECT background_col,comp_uga_col,comp_uag_col,comp_uaa_col,title_size,subheading_size,axis_label_size,marker_size,cds_marker_width,cds_marker_colour,legend_size from user_settings WHERE user_id = '{}';".format(user_id))
		result = (trips_cursor.fetchone())
		background_col = result[0]
		uga_col = result[1]
		uag_col = result[2]
		uaa_col = result[3]
		title_size = result[4]
		subheading_size = result[5]
		axis_label_size = result[6]
		marker_size = result[7]
		cds_marker_size = result[8]
		cds_marker_colour = result[9]
		legend_size = result[10]
	trips_connection.close()
		
	if tran != "":
		return riboflask_compare.generate_compare_plot(tran, ambiguous, minread, maxread, master_filepath_dict, "y", {}, 
											ribocoverage, organism,normalize,short_code,background_col,hili_start,
											hili_stop,comp_uag_col,comp_uga_col,comp_uaa_col,title_size, subheading_size,axis_label_size, marker_size,cds_marker_size,cds_marker_colour,
											legend_size,transcriptome)

			
	else:
		return_str =  "ERROR! Could not find any transcript corresponding to {}".format(tran)
		return return_str










