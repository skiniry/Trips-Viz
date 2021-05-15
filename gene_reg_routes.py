from flask import Blueprint, render_template, abort, request
import sqlite3
import ast
import os,math
import numpy as np
from scipy.stats import  mannwhitneyu,percentileofscore
from werkzeug import secure_filename
import config
from scipy.spatial.distance import cosine

from scipy.stats import spearmanr, pearsonr,kendalltau

gene_regulation_page = Blueprint("gene_regulation_page", __name__, template_folder="templates2")
@gene_regulation_page.route('/<organism>/<transcriptome>/gene_regulation/')
def gene_regulationpage(organism,transcriptome):
	#ip = request.environ['REMOTE_ADDR']
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

	connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips.sqlite")
	cursor = connection.cursor()
	advanced = False
	if user != None:
		cursor.execute("SELECT advanced from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		print "\n\n\n\nAdvanced {}".format(result)
		if result[0] == 1:
			advanced = True
		else:
			advanced = False
	gene_reg_tables = []
	for filename in os.listdir("/home/DATA/www/tripsviz/tripsviz/gene_reg_tables"):
		if filename.endswith(".csv"):
			gene_reg_tables.append(filename.strip(".csv"))

	return render_template('gene_regulation.html', studies_dict={}, accepted_files={},user=user,
						   organism=organism,default_tran="",local=local,transcriptome=transcriptome,advanced=advanced,
						   seq_types={},studyinfo_dict={},gene_reg_tables=gene_reg_tables)



def parse_csv_file(filepath,filename):
	openfile = open(filepath).readlines()
	table_html = "<table style  class='prediction_table hover'><tr><thead><th>File name {}</th></thead></tr>".format(filename)
	dict_user = {}
	count = 5
	for row in openfile:
		print "ROW", row
		try:
			splitrow = row.split("\t")
			dict_user[splitrow[0]] =  float(splitrow[1])
		except:
			count += -1
			if count> 0:
				table_html += "<tr><td>Skip row {}</td></tr>".format(row)
			if count == 0:
				table_html += "<tr><td>Skip ..</td></tr>"


	infileopen = open("pmid/ortho")
	dict_ortho = {}
	dict_ortho2 = {}
	for line in infileopen:
		linesplit = line[:-1].split()
		dict_ortho[linesplit[1]] = linesplit[2]
		dict_ortho2[linesplit[2]] = linesplit[1]
		#print linesplit

	keys_user = set(dict_user.keys())
	dict1 = {}
	files_human = os.listdir("/home/DATA/www/tripsviz/tripsviz/z_combine_ribo/human")
	files_mouse = os.listdir("/home/DATA/www/tripsviz/tripsviz/z_combine_ribo/mouse")

	keys_user2 = {}
	#print keys_user
	if len(set(dict_ortho.keys()).intersection(keys_user)) > len(set(dict_ortho.values()).intersection(keys_user)):
		mouse = 0
		table_html += "<tr><td>Genes identifed as Human</td></tr>"
		for filea in files_human:
			infileopen = open("/home/DATA/www/tripsviz/tripsviz/z_combine_ribo/human/%s"%filea)
			dict1[filea] = {}
			infileopen.readline()
			for line in infileopen:
				linesplit = line.split("\t")
				dict1[filea][linesplit[0]] = (float(linesplit[2]))
		for filea in files_mouse:
			infileopen = open("/home/DATA/www/tripsviz/tripsviz/z_combine_ribo/mouse/%s"%filea)
			dict1[filea] = {}
			infileopen.readline()
			for line in infileopen:
				linesplit = line.split("\t")
				if linesplit[1] != "":
					dict1[filea][linesplit[1]] = (float(linesplit[2]))   #for mouse data (i only use orthologues)
	else:
		mouse =1
		table_html += "<tr><td>Genes identifed as Mouse</td></tr>"
		for filea in files_human:
			infileopen = open("/home/DATA/www/tripsviz/tripsviz/z_combine_ribo/human/%s"%filea)
			dict1[filea] = {}
			infileopen.readline()
			for line in infileopen:
				linesplit = line.split("\t")
				#print filea, linesplit
				if linesplit[1] != "":
					dict1[filea][dict_ortho[linesplit[1]]] = (float(linesplit[2]))
		for filea in files_mouse:
			infileopen = open("/home/DATA/www/tripsviz/tripsviz/z_combine_ribo/mouse/%s"%filea)
			dict1[filea] = {}
			infileopen.readline()
			for line in infileopen:
				linesplit = line.split("\t")
				dict1[filea][linesplit[0]] = (float(linesplit[2]))

	out_data =[]
	genes_used = []
	genes_discarded = []
	for filea in dict1.keys():
		common_keys = set(dict1[filea].keys())
		common_keys = common_keys.intersection(keys_user)

		list1,list2 = [],[]
		genes_discarded.append(len( set(keys_user)) - len(common_keys))
		genes_used.append( len(common_keys))
		if len(common_keys) < 4000:
			print filea
			continue
#if len(common_keys1) > len(common_keys):
		for key  in common_keys:
			list1.append(dict1[filea][key])
			list2.append(dict_user[key])


		spt = cosine(list1,list2)
		list1c = []
		for n in list1:
			if n > 0:
				list1c.append(n**2)
			else:
				list1c.append(-1*(n**2))
				#list1c = []
		list2c = []
		for n in list2:
			if n > 0:
				list2c.append(n**2)
			else:
				list2c.append(-1*(n**2))
		spe = pearsonr(list1c,list2c)
		#if 1:
			#print dict1[filea].keys()
		#if abs(spt-1) > 0.20:
		out_data.append((abs(spt-1), spt, filea,abs(spe[0])))
			#if filea == "Jan14_Sec61":
				#print (abs(spt-1), spt, filea)

	if mouse:
		dict_user2 = {}
		for key,values in dict_ortho.items():
			if dict_user.has_key(values):
				dict_user2[key] = dict_user[values]
		dict_user = dict_user2
		keys_user = set(dict_user.keys())
	genes_discarded = np.mean(genes_used)
	genes_used = min(genes_used)
	if genes_used < 10000:
		table_html += "<tr><thead><th><font color='red'>{} <font color='red'>{}</th><th>Comparison may be unreliable</th></thead></tr>".format("Genes recognised",genes_used)
		table_html += "<tr><thead><th><font color='red'>{}  <font color='red'>{}</th></thead></tr>".format( "Gene discarded",genes_discarded)
	else:
		table_html += "<tr><thead><th>Gene recognised {}</th></thead></tr>".format(genes_used)

		table_html += "<tr><thead><th>Gene discarded {}</th></thead></tr>".format( genes_discarded)
	
				#keys_user2
	table_html += "<tr><thead><th>Gene Clusters</th></thead></tr>"

	infileopen = open("pmid/description")
	dict_des = {}
	for line in infileopen:
		linesplit = line[:-1].split("\t")
		print line
		dict_des[linesplit[0]] = linesplit[1]

	#table_html += "<br>"
	#table_html += "<br>"
	#table_html += "<tr><thead><td>Ribo_cluster</td></thead></tr>"
	out_data2 = []
	ddx = dict_user.values()
	ddx.sort()
	ddx2 = ddx
	ddx2.reverse()
	for file1 in ["out"]:
		infileopen = open("pmid/%s"%file1)
		for line in infileopen:
			#if len(line) == 1: continue
			line_list = line.split("\t")[1:]
			#if len(line_list) == 0: continue
			if len(line_list) < 3: continue

			common_keys = keys_user.intersection(line_list)
			if len(common_keys)> 2:
				listaa = [0,0]
				for gnee in common_keys:
					listaa.append(dict_user[gnee])
				pvale2a = mannwhitneyu(listaa, ddx)[1]

				#pvale2a = pvale[1]
				if pvale2a < 0.05:
					pbb = mannwhitneyu(ddx[:len(listaa)],ddx)[1]
					if pbb == 0: pbb = 1e-309
					pvale1a = (math.log10(pvale2a))/math.log10(pbb)
					#if pvale1a >= 0.1:
					if 1:
						line_list.sort()
						line_list = ",".join(line_list)
					#description_d = 
						out_data2.append((round(pvale1a,3), round(np.median(listaa),2), line_list,pvale2a,line.split("\t")[0]))
			#if set(gene_list) - set(line_list) == set([]):
			#	secondary_tablestr += "{}:{}.,/".format(file1,", ".join(line_list))
	out_data2.sort()
	out_data2.reverse()

	table_html += "<tr><thead><th></th><th>Mann-Whitney(P)</th><th>log10(P)/log10(max(P))</th><th>Median Z-score</th><th>Gene List</th></thead></tr>"
	skip = []
	for a, b,c,d,e in out_data2:
		#print a,  b
		if set(c)-set(skip) == set([]):continue
		skip.extend(c)
		if a>  0.50 and abs(b) > 1:
			table_html += "<tr><td><font color='blue'>{}</td><td><font color='blue'>{:.4e}</td><td><font color='blue'>{}</td><td><font color='blue'>{}</td><td><font color='blue'>{}</td></tr>".format(e,d,a,b,c)
		#elif a  >.50 :
			#table_html += "<tr><td></td><td>{}</td><td><font color='blue'>{}</td><td>{}</td><td>{}</td></tr>".format(d,a,b,c)
		#elif abs(b) > 1:
			#table_html += "<tr><td></td><td>{}</td><td>{}</td><td><font color='blue'>{}</td><td>{}</td></tr>".format(d,a,b,c)
		#table_html += "<tr><td>{}</td><td>{}</td></tr>".format(splitrow[0],splitrow[1])
		else:
			table_html += "<tr><td>{}</td><td>{:.4e}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(e, d,a,b,c)

	table_html += "<br>"
	table_html += "<br>"
	#table_html += "<tr><thead><td>Ribo_cluster</td></thead></tr>"
	table_html += "<tr><thead><th>Gene expresssion Responses</th>"
	table_html += "<tr><thead><th>Name</th><th>Mann-Whitney(P)</th><th>log10(P)/log10(max(P))</th><th>Median Z-score</th><th><align='left'> Gene List</th></thead></tr>"
	out_data2 = []
	for file1 in ["supervised"]:
		infileopen = open("pmid/%s"%file1)
		for line in infileopen:
			line_list = line.split("\t")
			#if len(line_list) == 0: continue
			gene_list = line_list[1].split(",")[:-1]
			common_keys = keys_user.intersection(gene_list)
			if len(common_keys)> 2:
				listaa = [0,0]
				for gnee in common_keys:
					listaa.append(dict_user[gnee])
				pvale2a = mannwhitneyu(listaa, dict_user.values())[1]
				if pvale2a < 0.05:
					pbb = mannwhitneyu(ddx[:len(listaa)],ddx)[1]
					if pbb == 0: pbb = 1e-309

					pvale1a = (math.log10(pvale2a))/math.log10(pbb)
					if pvale1a >= 0.1:
						gene_list = ",".join(gene_list)
						out_data2.append((round(pvale1a,3), round(np.median(listaa),2), gene_list,pvale2a,line_list[0]))
	out_data2.sort()
	out_data2.reverse()
	for d, b,c,a,e in out_data2:
		if d > 0.5 and abs(b) > 1:
			table_html += "<tr><td>	<font color='blue'>{}</td><td><font color='blue'>{:.2e}</td><td>	<font color='blue'>{}</td><td><font color='blue'>{}</td><td>	<font color='blue'>{}</td></tr>".format(e,a,d,b,c)
		else:
			table_html += "<tr><td>	{}</td><td>{:.2e}</td><td>{}</td><td>{}</td><td>{}</td></tr>".format(e,a,d,b,c)

	table_html += "<br>"
	table_html += "<br>"
	table_html += "<tr><thead><th>Individual Experiment Comparison</th>"
	table_html += "<tr><thead><th>{}</th>".format(genes_discarded)
	table_html += "<tr><thead><th>Experiment</th><th></th><th>abs(1-Cosine)</th><th>Direction</th><th>Study Description</th></thead></tr>"
	out_data.sort()
	out_data.reverse()
	for a,c, b, d in out_data:

		#table_html += "<tr><td>{}</td><td>{}</td></tr>".format(splitrow[0],splitrow[1])
		dict_des.setdefault(b, "")
		#if d < 50 and abs(b) > 1:
		if c < 1:
			#if a >= 0.3:
				#table_html += "<tr><td><font color='blue'>{}</td><td></td><td><font color='blue'>{}</td><td><font color='blue'>{}</td><td><font color='blue'>{}</td</tr>".format(b,round(a,4),round(c,4),dict_des[b])
			#else:
				#table_html += "<tr><td>{}</td><td></td><td>{}</td><td>{}</td><td>{}</td</tr>".format(b,round(a,4),round(c,4),dict_des[b])

			if a >= 0.3:
				table_html += "<tr><td><font color='blue'>{}</td><td>{}</td><td><font color='blue'>{}</td><td>+</td><td><font color='blue'>{}</td</tr>".format(b,round(a,4),round(d,4),dict_des[b])
			else:
				table_html += "<tr><td>{}</td><td>{}</td><td>{}</td><td>+</td><td>{}</td</tr>".format(b,round(a,4),round(d,4),dict_des[b])
		else:
			if a >= 0.3:
				table_html += "<tr><td><font color='blue'>{}</td><td>{}</td><td><font color='blue'>{}</td><td>-</td><td><font color='blue'>{}</td</tr>".format(b,round(a,4),round(d,4),dict_des[b])
			else:
				table_html += "<tr><td>{}</td><td>{}</td><td>{}</td><td>-</td><td>{}</td</tr>".format(b,round(a,4),round(d,4),dict_des[b])

	table_html += "</table>"
	#table_html +='<img src="/home/DATA/www/tripsviz/tripsviz/pmid/pic_trulli.jpg" alt="Trulli" width="500" height="333">'
	print "returning html"
	return table_html


gene_regulation_query = Blueprint("gene_regulation_query", __name__, template_folder="templates")
@gene_regulation_query.route('/gene_regulation_query', methods=['POST'])
def generegulationquery():
	print "gene regulation called"
	if request.method == 'POST':
		uploaded_files = request.files.getlist("file")
		for f in uploaded_files:
			connection = sqlite3.connect('{}/trips.sqlite'.format(config.SCRIPT_LOC))
			connection.text_factory = str
			cursor = connection.cursor()
			#user = current_user.name
			#cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
			#result = (cursor.fetchone())
			#user_id = result[0]
			filename = secure_filename(f.filename)
			upload_file_path = "{}/static/tmp/{}".format(config.SCRIPT_LOC,filename)
			f.save("{}".format(upload_file_path))
			tablehtml = parse_csv_file(upload_file_path,filename)
			print tablehtml
			return tablehtml
		selected_file = request.form.get("selected_file")
		#if selected_file != "None":
		try:
			#print selected_file
			print len(selected_file)
			#print "None" in selected_file
			tablehtml = parse_csv_file("/home/DATA/www/tripsviz/tripsviz/gene_reg_tables/{}.csv".format(selected_file),selected_file)
			return tablehtml
		except: pass
	connection = sqlite3.connect("/home/DATA/www/tripsviz/tripsviz/trips.sqlite")
	cursor = connection.cursor()
	script_location = "/home/DATA/www/tripsviz/tripsviz"
	acceptable = 0
	data = ast.literal_eval(request.data)

	try:
		user = current_user.name
		cursor.execute("SELECT user_id from users WHERE username = '{}';".format(user))
		result = (cursor.fetchone())
		user_id = result[0]
		#cursor.execute("SELECT user_id from users WHERE username = 'master';")
		#result = (cursor.fetchone())
	except:
		user = None

	organism = data["organism"]
	returnstr = ""
	gene_list = (data["tran_list"].replace(" ","")).split(",")
	#drop_file = data["drop_file"]
	print "gene list", gene_list
	#if gene_list == [""]:
	#	print "parsing premade csv file"
	#	tablehtml = parse_csv_file("/home/DATA/www/tripsviz/tripsviz/gene_reg_tables/{}.csv".format(drop_file),drop_file)
	#	print tablehtml
	#	return tablehtml
	transcriptome = data["transcriptome"]
	print "data", data
	gene_list = set([gxx.upper() for gxx in gene_list])
	if len(gene_list) == 0:
		#break
		pass
	else :
		secondary_tablestr = ""
		#for file1 in ["8_(high_confidence)","7_(high_confidence)","6_(high_confidence)","5_(low_confidence)","4_(low_confidence)","3_(low_confidence)", "2_(low_confidence)"]:
		if 0:
			for file1 in ["Ribo","TE","RNA"]:
				infileopen = open("pmid/%s"%file1)
				for line in infileopen:
					line_list = line.split("\t")[:-1]
					line_list.sort()

					if set(gene_list) - set(line_list) == set([]):
						secondary_tablestr += "{}:{}.,/".format(file1,", ".join(line_list))

		#secondary_tablestr += "{}:{}.,/".format(6,"ATF4,ATF5")
		#secondary_tablestr += "{}:{}.,/".format(7,"ATF4,ATF5,SLC35A4")
		#if len(gene_list) == 1:
		if 1:
			if 0:
				infileopen = open("pmid/dip_pair_top")
				start_dict = {}
				for line in infileopen:
					linesplit = line[:-1].split("\t")
					#print linesplit
					if linesplit[0] == linesplit[1]:
						start_dict[linesplit[0]] = float(linesplit[2])

			#infileopen = open("cluster_euclidean_average2b")
			dicta = {}
			dictb = {}
			dictc = {}
			dictd = {}
			used = []
			#infileopen = open("pmid/cluster_euclidean_average2b")
			##show = []
			##for line in infileopen:
				##if line[-1] != "\n": continue
				##linesplit = line[:-1].split(",")
				##list1 = linesplit[2].split("\t")
				##list1 = [nxx for nxx in list1 if nxx != ""]
				##list2 = linesplit[3].split("\t")
				##list2 = [nxx for nxx in list2 if nxx != ""]
				###list3.extend(linesplit[3].split("\t"))

				##list3 = linesplit[2].split("\t")
				##list3.extend(linesplit[3].split("\t"))
				##list3 = [nxx for nxx in list3 if nxx != ""]
				###stb = "".join([for nxx in list3])
				##stb = " ".join(list3)
				###listdd.append((float(linesplit[0]), list3))
				##if set(gene_list) - set(list3)  == set([]):
					##if float(linesplit[0]) > 3:
						###show = [linesplit[0], stb]

			#geen = gene_list[0]
			#dicta = {}
			#dictb = {}
			#dictc = {}
			#for line in infileopen:
				#linesplit = line[:-1].split(",")
				##print linesplit[2], 1, linesplit[3], 2
				##print
				#if line[-1] != "\n": continue
				#list1 = linesplit[2].split("\t")
				#list1 = [nxx for nxx in list1 if nxx != ""]
				#list2 = linesplit[3].split("\t")
				#list2 = [nxx for nxx in list2 if nxx != ""]
				##list3.extend(linesplit[3].split("\t"))

				#list3 = linesplit[2].split("\t")
				#list3.extend(linesplit[3].split("\t"))
				#list3 = [nxx for nxx in list3 if nxx != ""]
				##for key in list3:
					##dictb[key] = float(linesplit[0])
				##print list1, list2
				#if len(list1) == 1 and len(list2):
					#dicta[str(list3)] = "(%s,%s)"%(list1[0],list2[0])
					#dictb[str(list2)] = float(linesplit[0])
				#elif len(list1) ==1:
					#dictb[str(list1)] = float(linesplit[0])
					#dicta[str(list3)] = "(%s,%s)"%(list1[0],dicta[str(list2)])
				#elif len(list2) ==1:
					#dicta[str(list3)] = "(%s,%s)"%(dicta[str(list1)],list2[0])
				#else:
					#dicta[str(list3)] = "(%s,%s)"%(dicta[str(list1)],dicta[str(list2)])

				##dictb[str(list1)] = float(linesplit[0])
				##dictb[str(list2)] = float(linesplit[0])
				##dictb[str(list3)] =
				##print dictb
				#dictc[str(list3)] = float(linesplit[0])
				#if len(list1) != 1:
					#dictd[str(list1)] = dictc[str(list1)]- float(linesplit[0])
				#else:
					#dictd[str(list1)] = start_dict[list1[0]]- float(linesplit[0])

				#if len(list2) != 1:
					#dictd[str(list2)] = dictc[str(list2)]- float(linesplit[0])
				#else:
					#dictd[str(list2)] = start_dict[list2[0]]- float(linesplit[0])

				##print  dictb[str(list1)]- float(linesplit[0])
				##if len(list1) == 1 and len(list2):
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),list2[0],float(linesplit[0]))
				##elif len(list1) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])
				##elif len(list2) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),list2[0],float(linesplit[0]))
				##else:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])

				#if geen in list3:
					#if len(list3) > 500 :break

					##print dicta[str(list3)]
				#used.append(list3)
				#try:
					#used.remove(list1)
				#except:
					#continue

				#try:
					#used.remove(list2)
				#except:
					#continue


			#dicta1 = {}
			#infileopen.seek(0)
			#for line in infileopen:
				#linesplit = line[:-1].split(",")
				#if line[-1] != "\n": continue
				##print linesplit[2], 1, linesplit[3], 2
				##print
				#list1 = linesplit[2].split("\t")
				#list1 = [nxx for nxx in list1 if nxx != ""]
				#list2 = linesplit[3].split("\t")
				#list2 = [nxx for nxx in list2 if nxx != ""]
				##list3.extend(linesplit[3].split("\t"))

				#list3 = linesplit[2].split("\t")
				#list3.extend(linesplit[3].split("\t"))
				#list3 = [nxx for nxx in list3 if nxx != ""]
				##for key in list3:
					##dictb[key] = float(linesplit[0])
				##print list1, list2
				#if len(list1) == 1 and len(list2):
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],dictd[str(list1)],list2[0],dictd[str(list2)])
					##dictb[str(list2)] = float(linesplit[0])
				#elif len(list1) ==1:
					##dictb[str(list1)] = float(linesplit[0])
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],dictd[str(list1)],dicta1[str(list2)],dictd[str(list2)])
				#elif len(list2) ==1:
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(dicta1[str(list1)],dictd[str(list1)],list2[0],dictd[str(list2)])
				#else:
					#dicta1[str(list3)] = "(%s:%s,%s:%s)"%(dicta1[str(list1)],dictd[str(list1)],dicta1[str(list2)],dictd[str(list2)])

				##dictb[str(list1)] = float(linesplit[0])
				##dictb[str(list2)] = float(linesplit[0])
				##dictb[str(list3)] =
				##print dictb
				##dictc[str(list3)] = float(linesplit[0])
				##if len(list1) != 1:
					##dictd[str(list1)] = dictc[str(list1)]- float(linesplit[0])
				##else:
					##dictd[str(list1)] = 0

				##if len(list2) != 1:
					##dictd[str(list2)] = dictc[str(list2)]- float(linesplit[0])
				##else:
					##dictd[str(list2)] = 0

				##print  dictb[str(list1)]- float(linesplit[0])
				##if len(list1) == 1 and len(list2):
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),list2[0],float(linesplit[0]))
				##elif len(list1) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(list1[0],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])
				##elif len(list2) ==1:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),list2[0],float(linesplit[0]))
				##else:
					##dicta[str(list3)] = "(%s:%s,%s:%s)"%(dicta[str(list1)],float(linesplit[0]),dicta[str(list2),float(linesplit[0])])

				#if geen in list3:
					#if len(list3) > 500 :break

					##print dicta[str(list3)]
				#used.append(list3)
				#try:
					#used.remove(list1)
				#except:
					#continue

				#try:
					#used.remove(list2)
				#except:
					#continue

				##used.remove(list2)
				##else:
			#from ete3 import Tree
			##from ete3 import Tree, TreeStyle
			##circular_style = TreeStyle()
			##circular_style.mode = "c" # draw tree in circular mode
			##circular_style.scale = 20
			#for x in  used:
				#if geen in x:
					#xp = dicta1[str( x)]
					###print xp, "xp"
					#t = Tree("%s;"%xp, format=1)
					###t = Tree("((D:0.723274,F:0.567784)1.000000:0.067192,(B:0.279326,H:0.756049)1.000000:0.807788);",format = 0)
					###print xp
					###print t
					##returnstr += ",{},,,,,,.,/".format(xp)

			#if gene_list[0] in list3:
			#from ete3 import Tree
			#if show != []:
				#returnstr += "{},{},,,,,,.,/".format(show[0], show[1])

		infileopen = open("pmid/pmid")
		infileopen.readline()
		dict_pmid = {}
		for line in infileopen:
			linesplit = line[:-1].split("\t")
			dict_pmid[linesplit[0]] = linesplit[1]

		infileopen = open("pmid/description")
		#infileopen.readline()
		dict_des = {}
		for line in infileopen:
			linesplit = line[:-1].split("\t")
			dict_des[linesplit[0]] = "" #linesplit[0]
			if len(linesplit)> 1:
				dict_des[linesplit[0]] = linesplit[1]
		infileopen = open("pmid/new_name")
		#infileopen.readline()
		dict_name = {}
		for line in infileopen:
			linesplit = line[:-1].split("\t")
			dict_name[linesplit[0]] = linesplit[0]
			if len(linesplit)> 1:
				dict_name[linesplit[0]] = linesplit[1]

			#else:

		dict_pathway = {}
		if 0:
		#list_file2 =["mTOR_inhibition"]
			list_file2 = os.listdir("cluster/")
			list_out = []
			for filea in list_file2 :
				infileopen = open("cluster/%s"%filea)
				line = infileopen.readline()
				linesplit = line[:-1].split("\t")
				#dict_pathway[filea] =
				for key in linesplit[1:]:
					try:
						aba = float(key)
					except:
						dict_pathway.setdefault(key,[])
						if filea not in dict_pathway[key]:
								dict_pathway[key].append((0, filea))
								dict_pathway[key].sort()
								dict_pathway[key].reverse()
						#if dict_pathway[key] == "":
							#dict_pathway[key] = "{}".format(filea)
						#else:
							#dict_pathway[key] ="{}_{}".format(dict_pathway[key],filea)
				count = 0
				count1 = ""
			#for line in infileopen:
				#linesplit = line[:-1].split("\t")
				#if linesplit[0] in gene_list:
					##quer2.append(abs(float(linesplit2[1])))
					#count += 1
					##count1.append(linesplit[0])
					#if count1 == "":
						#count1 = linesplit[0]
					#else:
						#count1 = "{}_{}".format(count1, linesplit[0])

			##list_out.append((count, filea,count1))
		#list_out.sort()
		#list_out.reverse()
		#for z_score in list_out:
			#if z_score[0]> 0:
				#returnstr += "{}{},{},None,None.,/".format(z_score[0],z_score[1], z_score[2])
	#returnstr += "None,None,None,None.,/"
		#infileopen = open("patrick_regulated_in")
		#for line in infileopen:
			#linesplit = line.split("\t")
			#if gene_list[0] == linesplit[0]:
				#data = linesplit[1:]
				#data1 = []
				#data2 = []
				#for x_i, x in enumerate(data):
					#if x_i %2 ==1:
						#data1.append(x)
					#else:
						#data2.append(x)

				#break
		#for study, z_score in zip(data1,data2):
			##returnstr += "{},{},None,None.,/".format(z_score, study)
			#stra = ""
			#dict_pathway.setdefault(study,[])
			#for tuple1 in dict_pathway[study]:
				#if stra == "":
					#stra = "{}".format(tuple1[1])
				#else:
					#stra = "{}_{}".format(stra, tuple1[1])
			#annotation_popup = "Annotation description here_top"
			#returnstr += "{},{},None,{},{}.,/".format(study, z_score,stra,annotation_popup)
	#else:
	#if 1:
		list_file2 = os.listdir("z_all/")
		#list_file2 = ["Alvarez-Dominguez_24_v0hr_erythroid_differentiation", "Alvarez-Dominguez_33_v24hr_erythroid_differentiation", "Alvarez-Dominguez_48_v33hr_erythroid_differentiation", "Arango18_NAT10_ko", "Aviner17_hnRNP_c_kd", "Bennett_msi2_shRNA", "bercovich_kinori_influenza_2_4hr", "bercovich_kinori_influenza_2hr", "bercovich_kinori_influenza_4_8hr", "Blair_hesc_npc", "Blair_neu14_hesc", "Blair_neu14_npc", "Blanco16_Nsun1_ko", "castaneda_mael_ko", "Castelo17_zt0", "Castelo17_zt2", "Castelo17_zt4", "Castelo17_zt6", "Castelo17_zt8", "Castelo17_zt10", "Castelo17_zt12", "Castelo17_zt14", "Castelo17_zt16", "Castelo17_zt18", "Castelo17_zt20", "Castelo17_zt22", "dai_2hr_vaccinia_infection", "dai_2v4hr_vaccinia_infection", "dai_4v8hr_vaccinia_infection", "Diaz-Munoz_Hur_KO", "Diaz-Munoz_Hur_KO_LPS", "dima_aicar_f", "dimaelife_arsenite", "dimaGB_ogd20", "dimaGB_ogd20_40", "dimaGB_ogd40_60", "dima_mechano_h1", "dima_mechano_I1", "dima_myxo_anoxia", "dima_myxo_E1", "dima_myxo_g5_hct116sca", "dima_myxo_G1_hct116", "
#dima_myxo_H1", "dima_phendc3h1", "dima_phendc3k", "dima_pieri_f", "dima_pieri_g5_hct116sca", "dima_pieri_G1_hct116", "dima_pieri_H1", "dima_te_myxo_sco2_anoxia", "Eichhorn_mi155ko_2h", "Eichhorn_mi155ko_4h", "Eichhorn_miR-1", "fijalkowska_eif1_kd", "Fradejas_Secisbp2_ko", "Fradejas_Trsp_ko", "Gameiro_EtOH_NoBCAA", "Gameiro_EtOH_NoCys", "Gameiro_EtOH_NoG", "Gameiro_EtOH_NoQ", "Gameiro_TAM_NoBCAA", "Gameiro_TAM_NoCys", "Gameiro_TAM_NoG", "Gameiro_TAM_NoQ", "Gameiro_TAM_torin", "Gao17_liver_fasting", "Gao17_MEF_starvation", "gao_aminoastarvation", "gao_eif2phos_mimetic", "Ginossar12_HCMV_24v5hr", "Gonzalez14_Tumor_v_normal", "goodarzi16_argCCG_oe", "goodarzi16_gluUUC_oe", "Goodarzi_highly_v_lowly_metastatic", "Grabole_AZD-8055_tsc2-cellline", "Grabole_AZD-8055_tsc2+cellline", "Grow_Rec_OE", "Guo14_te_miR-1_v_control", "howard13_A37G_0.1_0_selenocysteine", "howard13_A37G_2_0.1_selenocysteine", "howard13_wt_0.1_0_selenocysteine", "howard13_wt_2_0.1_selenocysteine", "hsieh_pp242", "Irigoyen_coronavirus_1h", "
##Irigoyen_coronavirus_8h", "iwasaki_roca", "Jackson18_LPS_simulation", "Jakobson_KO_EEF1AKMT4", "jan_0", "jan_2", "jan_4", "jan_6", "jan_8", "jan_10", "jan_12", "jan_14", "jan_16", "jan_18", "jan_20", "jan_22", "Janich15_zt0", "Janich15_zt2", "Janich15_zt4", "Janich15_zt6", "Janich15_zt8", "Janich15_zt10", "Janich15_zt12", "Janich15_zt14", "Janich15_zt16", "Janich15_zt18", "Janich15_zt20", "Janich15_zt22", "Ji15_Src_induced_transformation_1h", "Ji15_Src_induced_transformation_4h", "Jiang17_hypoxia1h", "khajuria_rps5", "khajuria_rps19", "Laguesse_Elp3cKO", "Leshkowitz19_eIF1A_sil", "Li_cyscys_deprivation", "Li_glucose_deprivation", "Liu18_FMRP_ko", "Loayza-Puch13_Nutlin-3a_2hv0", "Loayza-Puch13_Nutlin-3a_4hv2h", "Loayza-Puch13_Nutlin-3a_6hv4h", "Loayza-Puch13_Nutlin-3a_6hv19h", "Loayza-Puch13_pre-senescence_senescence", "Loayza-Puch13_profileration_quiescence", "Loayza-Puch14_myc_v_control", "Loayza-Puch16_tumour_v_control", "Loayza-Puch_siSL3A2", "Loayza-Puch_TGFb", "Murat_dhx9_kd", "oh_arsenite_ddx3", "oh_
#arsenite_ddx3_r534h", "oh_arsenite_wt", "Paolini_tunicamycin", "park16_s_v_m", "park_torin24", "Razooky17_RP8_v_mock", "reid_0_6h_denv1", "reid_0_6h_denv2", "reid_6_12h_denv1", "reid_6_12h_denv2", "reid_12_24h_denv1", "reid_12_24h_denv2", "reid_24_48h_denv1", "reid_24_48h_denv2", "reid_48_72h_denv1", "reid_48_72h_denv2", "Reid13_cyto_Thapsigargin_0.5_v_0hr", "Reid13_cyto_Thapsigargin_1_v_0.5hr", "Reid13_cyto_Thapsigargin_2_v_1hr", "Reid13_cyto_Thapsigargin_4_v_2hr", "Reid13_er_Thapsigargin_0.5_v_0hr", "Reid13_er_Thapsigargin_1_v_0.5hr", "Reid13_er_Thapsigargin_2_v_1hr", "Reid13_er_Thapsigargin_4_v_2hr", "reid_cytosol_v_ER", "Ricci13_Stau1_overexpression_v_KnockDown", "rubio_silvestrol_2h", "Rutkowski_hsv1_2h", "Rutkowski_hsv1_4h", "Sendoel_SOX2_induction", "Shalgi_2hr_severe_heatshock", "Shalgi_8hr_mild_heatshock", "shi_siYTHDF3", "sid_tunamycin", "Sims14_tumour_normal", "Simsek17_siPKM", "Slobodin17_Campthotecin_v_control", "Slobodin17_Nutlin", "Stumpf_G1_v_S", "Stumpf_M_v_G1", "Stumpf_S_v_M", "Su15_IFNg",
#"#tanenbaum_g2_g1", "tanenbaum_g2_m", "tanenbaum_m_g1", "Thoreen_torin1", "Tichon_NORD_si", "Tirosh_HCMV_5h", "Tirosh_HCMV_12h", "Tirosh_HCMV_24h", "Tirosh_irradiated_HCMV", "Tirosh_type1_interferon", "wang_siYTHDF2", "Wein14_mut_wt", "werner_hES_shKBTBD8", "werner_ND1_shKBTBD8", "werner_ND3_shKBTBD8", "werner_ND6_shKBTBD8", "witta_bortezomib1.5_3h", "witta_bortezomib1.5h", "witta_bortezomib3_6h", "witta_bortezomib6_9h", "witta_bortezomib9_12h", "wolfe_silvestrol", "xu_correct", "xu_mutant", "you15_siKrr1", "zang_mi_oe", "Zhou15_heatshock", "Zhou15_heatshock_FTO", "Zhou15_heatshock_YTHDF2", "Zur_M_G1"]
		list_file = []
		for filea in list_file2 :
				if ".py" in filea: continue
				#try:
					#infileopen = open("z_combine_rna/%s"%filea)
				##dict2[filea] = 0
					#change2rna = [0,0, len(gene_list)]
					#dicta = {}
					#quer1, quer2 = [],[]
					#for line2 in infileopen:
						#linesplit2 = line2.split("\t")
						#quer1.append((float(linesplit2[2])))
						#if linesplit2[1] in gene_list:
							#quer2.append((float(linesplit2[2])))
							#if float(linesplit2[2]) >0:
								#change2rna[0] += 1
							#else:
								#change2rna[1] += 1
					#rna_score = mannwhitneyu(quer1, quer2)[1]
				#except :
					#rna_score = 1
					#change2rna = [0,0, len(gene_list)]

				#rna_score = np.median(quer2)

				infileopen = open("z_all/%s"%filea)
				infileopen.readline()
				#dict2[filea] = 0
				dicta = {}
				quer1, quer2 = [],[]
				quer1rna, quer2rna = [],[]
				quer1te, quer2te = [],[]
				change2ribo = [0,0, len(gene_list)]
				change2rna = [0,0, len(gene_list)]
				change2b = [0,0, len(gene_list)]
				for line2 in infileopen:
					linesplit2 = line2.split("\t")
					quer1.append((float(linesplit2[2])))
					try:
						quer1te.append((float(linesplit2[3])))
						quer1rna.append((float(linesplit2[4])))
					except:pass
					if linesplit2[1] in gene_list:
						quer2.append((float(linesplit2[2])))
						if float(linesplit2[2]) >0:
							change2ribo[0] += 1
						else:
							change2ribo[1] += 1

						try:
							quer2te.append((float(linesplit2[3])))
							quer2rna.append((float(linesplit2[4])))
							if float(linesplit2[3]) >0:
								change2b[0] += 1
							else:
								change2b[1] += 1

							if float(linesplit2[4]) >0:
								change2rna[0] += 1
							else:
								change2rna[1] += 1


						except:pass

				ribo_score = mannwhitneyu(quer1, quer2)[1]
				if len(quer2rna) == 0:
					te_score = 2
					rna_score = 2
				else:
					rna_score = mannwhitneyu(quer1rna, quer2rna)[1]
					te_score = mannwhitneyu(quer1te, quer2te)[1]

				#ribo_score = np.median(quer2)

				#try:
					#infileopen = open("z_combine/%s"%filea)
					##dict2[filea] = 0
					#dicta = {}
					#quer1, quer2 = [],[]
					#change2b = [0,0, len(gene_list)]
					#for line2 in infileopen:
						#linesplit2 = line2.split("\t")
						#quer1.append((float(linesplit2[2])))
						#if linesplit2[1] in gene_list:
							#quer2.append((float(linesplit2[2])))
							#if float(linesplit2[2]) >0:
								#change2b[0] += 1
							#else:
								#change2b[1] += 1
					#te_score = mannwhitneyu(quer1, quer2)[1]
				#except :
					#te_score = 1
					#change2b = [0,0, len(gene_list)]
						#if float(linesplit2[1]) > 0:
							#dict2[filea]  += 1
						#else:
							#dict2[filea]  += -1
				if len(quer2) > 0:
					#if np.mean(quer2) > 1:
					dop= [te_score,ribo_score,rna_score ]
					#dop= [np.median(quer2),ribo_score,rna_score ]
					#dop = mannwhitneyu(quer1, quer2)[1]
					quer1.sort()
					quer2a = quer1[:len(quer2)]
					quer2b = quer1[-len(quer2):]
					ideal_score = min([mannwhitneyu(quer1, quer2a)[1],mannwhitneyu(quer1, quer2b)[1]])
					#if dop < rna_score:
					if 1:

						list_file.append((max(dop), dop,  filea,round(percentileofscore(quer1, abs(max(quer2))),2),ideal_score,len(quer2),change2b, change2ribo, change2rna))
				#quer1.sort()
				#list_file_200[filea] = quer1[-1000]
		#list_file.sort()

		#returnstr += ",{},{},,,,.,/".format( "#Min Pvalue#",list_file[0][4])
		#returnstr += ",{},{},,.,/".format(10**(np.log10(list_file[0][-1])*4/5), "#STRONG SCORE#")

		for study in list_file:
			#if study[0] > 0.05: continue
			if study[3] < 70: continue
			#if np.log10(study[0])< np.log10(study[3])*1/2:
			if 1:


				dict_pathway.setdefault(study[2],[])
				#dict_pathway2.setdefault(study[1],"")
				stra = ""
				#for tuple1 in dict_pathway[study[1]]:
					#if stra == "":
						#stra = "{}".format(tuple1[1])
					#else:
						#stra = "{}_{}".format(stra, tuple1[1])
				dict_pmid.setdefault(study[2].split("_")[0], "na")
				dict_des.setdefault(study[2], "")
				dict_name.setdefault(study[2], study[2])
				annotation_popup = dict_des[study[2]]
				#returnstr += "{},{},{},{},{},.,/".format(study[1], study[0],dict_pmid[study[1].split("_")[0]],stra,annotation_popup)
				missing = study[6][2] - study[6][1]-study[6][0]
				changestr = "{}_{}_{}".format(study[6][0],study[6][1], missing)
				missing = study[7][2] - study[7][1]-study[7][0]
				changestr_ribo = "{}_{}_{}".format(study[7][0],study[7][1], missing)
				missing = study[8][2] - study[8][1]-study[8][0]
				changestr_rna = "{}_{}_{}".format(study[8][0],study[8][1], missing)

				#returnstr += "{},{},{},{},{},{},{}.,/".format(dict_name[study[1]], dict_pmid[study[1].split("_")[0]],study[0],study[2],changestr,stra,annotation_popup)
				returnstr += "{},{},{},{},{},{},{},{},{}.,/".format(dict_name[study[2]], dict_pmid[study[2].split("_")[0]],study[1][0],changestr,study[1][1],changestr_ribo,study[1][2],changestr_rna,annotation_popup)
				#print "return str noW", returnstr
	#print returnstr
	returnstr += "???"+secondary_tablestr
	return returnstr
