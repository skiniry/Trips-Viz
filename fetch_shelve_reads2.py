import os
import collections
from bokeh.palettes import all_palettes
from sqlitedict import SqliteDict

# Merge two dictionaries
def merge_dicts(dict1,dict2):
	for readlen in dict2:
		if readlen not in dict1:
			dict1[readlen] = dict2[readlen]
		else:
			for pos in dict2[readlen]:
				if pos in dict1[readlen]:
					dict1[readlen][pos] += dict2[readlen][pos]
				else:
					dict1[readlen][pos] = dict2[readlen][pos]
	return dict1


# Create dictionary of read counts at each position in a transcript
def get_reads(read_type, min_read, max_read, tran, user_files,tranlen,coverage, organism, subcodon, noisered, primetype, filetype, readscore, secondary_readscore=1,pcr=False,get_mismatches=False):
	mismatch_dict = collections.OrderedDict()
	mismatch_dict["A"] = {}
	mismatch_dict["T"] = {}
	mismatch_dict["G"] = {}
	mismatch_dict["C"] = {}
	if get_mismatches == True:
		for i in range(0,tranlen+1):
			mismatch_dict["A"][i] = 0
			mismatch_dict["T"][i] = 0
			mismatch_dict["G"][i] = 0
			mismatch_dict["C"][i] = 0


	master_dict = collections.OrderedDict()
	three_frame_dict = {0:{},1:{},2:{}}
	for i in range(0,tranlen+1):
		master_dict[i] = 0
	master_file_dict = {}
	list_master_dict = {}

	# first make a master dict consisting of all the read dicts from each filename
	offset_dict = {}
	secondary_offset_dict = {}

	if read_type == "unambig":
		if filetype in user_files:
			for file_id in user_files[filetype]:
				#print "fILE ID ", file_id
				filename = user_files[filetype][file_id]
				if os.path.isfile(filename):
					sqlite_db = SqliteDict(filename, autocommit=False)
				else:
					return "File not found",filename.split("/")[-1]
				#print "Filename", filename
				accepted_offsets = {}
				accepted_secondary_offsets = {}
				if "offsets" in sqlite_db:
					all_offsets = sqlite_db["offsets"][primetype]["offsets"]
					scores = sqlite_db["offsets"][primetype]["read_scores"]
				else:
					all_offsets = {}
					scores = {}
				if scores == {}:
					for i in range(min_read,max_read):
						scores[i] = 1
						all_offsets[i] = 15
				for rl in scores:
					if float(scores[rl]) >= float(readscore):
						if rl in all_offsets:
							accepted_offsets[rl] = all_offsets[rl]
						else:
							accepted_offsets[rl] = 15

				if accepted_offsets != {}:
					offset_dict[filename] = {}
					for length in accepted_offsets:
						#if length not in accepted_secondary_offsets:
						offset_dict[filename][length] = accepted_offsets[length]
						
				if get_mismatches == True:
					try:
						sqlite_db_seqvar = sqlite_db[tran]["seq"]
					
						for pos in sqlite_db_seqvar:
							#convert to one based
							fixed_pos = pos+1
							for char in sqlite_db_seqvar[pos]:
								#if char != "N":
								if fixed_pos not in mismatch_dict[char]:
									mismatch_dict[char][fixed_pos] = 0
								count = sqlite_db_seqvar[pos][char]
								mismatch_dict[char][fixed_pos] += count
					except Exception as e:
						pass

				try:
					alltrandict = sqlite_db[tran]
					#print "alltran dict", alltrandict
					unambig_tran_dict = alltrandict["unambig"]
					#print "unambig tran dict", unambig_tran_dict
					if pcr == True:
						if "unambig_pcr" in alltrandict:
							unambig_tran_dict = merge_dicts(unambig_tran_dict, alltrandict["unambig_pcr"])
					sqlite_db.close()
					master_file_dict[filename] = unambig_tran_dict
				except Exception as e:
					#print "error", e
					pass
	else:
		if filetype in user_files:
			for file_id in user_files[filetype]:
				filename = user_files[filetype][file_id]
				if os.path.isfile(filename):
					sqlite_db = SqliteDict(filename, autocommit=False)
				else:
					return "File not found",mismatch_dict
				accepted_offsets = {}
				accepted_secondary_offsets = {}
				if "offsets" in sqlite_db:
					all_offsets = sqlite_db["offsets"][primetype]["offsets"]
					scores = sqlite_db["offsets"][primetype]["read_scores"]
				else:
					all_offsets = {}
					scores = {}
				#oyster has no readscores so subcodon profiles not displaying, so i put this in to give every readlength a score of 1
				if scores == {}:
					for i in range(min_read,max_read):
						scores[i] = 1
						all_offsets[i] = 15

				for rl in scores:
					if float(scores[rl]) >= float(readscore):
						if rl in all_offsets:
							accepted_offsets[rl] = all_offsets[rl]
						else:
							accepted_offsets[rl] = 15

				if accepted_offsets != {}:
					offset_dict[filename] = {}
					for length in accepted_offsets:
						offset_dict[filename][length] = accepted_offsets[length]
				try:
					sqlite_db_seqvar = sqlite_db[tran]["seq"]
					for pos in sqlite_db_seqvar:
						#convert to one based
						fixed_pos = pos+1
						for char in sqlite_db_seqvar[pos]:
							if char != "N":
								if fixed_pos not in mismatch_dict[char]:
									mismatch_dict[char][fixed_pos] = 0
								count = sqlite_db_seqvar[pos][char]
								mismatch_dict[char][fixed_pos] += count
				except:
					pass

				try:
					alltrandict = sqlite_db[tran]
					sqlite_db.close()
					unambig_tran_dict = alltrandict["unambig"]
					if "ambig" in alltrandict:
						ambig_tran_dict = alltrandict["ambig"]
					else:
						ambig_tran_dict = {}
					# TODO: Change merge_dicts to take a list of dicts instead of two
					trandict = merge_dicts(unambig_tran_dict, ambig_tran_dict)
					if pcr == True:
						if "unambig_pcr" in alltrandict:
							trandict = merge_dicts(trandict, alltrandict["unambig_pcr"])
						if "ambig_pcr" in alltrandict:
							trandict = merge_dicts(trandict, alltrandict["ambig_pcr"])
					master_file_dict[filename] = trandict
				except Exception as e:
					#print "error: ", e
					pass	

	#Next check coverage, if that's true then calculate coverage for each rl and return dict
	if coverage == True and subcodon == False:
		for filename in master_file_dict:
			for readlen in master_file_dict[filename]:
				if readlen >= min_read and readlen <= max_read:
					for pos in master_file_dict[filename][readlen]:
						count = master_file_dict[filename][readlen][pos]
						if pos != 0 and pos-1 not in master_dict:
							master_dict[pos-1] = 0
						for i in range(pos,pos+(readlen+1)):
							if i in master_dict:
								master_dict[i] += count
							else:
								master_dict[i] = count
						#use this so line graph does not have 'ramps'
						if i+1 not in master_dict:
							master_dict[i+1] = 0
		return master_dict, mismatch_dict

	#Next check if subcodon is true, if not just give an offset of 15 to everything
	if subcodon == False:
		sorted_master_dict = collections.OrderedDict()
		for filename in master_file_dict:
			for readlen in master_file_dict[filename]:
				if readlen >= min_read and readlen <= max_read:
					for pos in master_file_dict[filename][readlen]:
						#use this so line graph does not have 'ramps'
						count = master_file_dict[filename][readlen][pos]
						offset_pos = pos+15
						master_dict[offset_pos] += count
						if offset_pos+1 not in master_dict:
							master_dict[offset_pos+1] = 0
		return master_dict, mismatch_dict
	

	#Fetching subcodon reads
	if subcodon == True:
		if coverage == False:
			for filename in master_file_dict:
				if filename not in offset_dict:
					continue
				for readlen in master_file_dict[filename]:
					if readlen >= min_read and readlen <= max_read:
						if readlen in offset_dict[filename]:
							offset = offset_dict[filename][readlen]
							for pos in master_file_dict[filename][readlen]:
								count = master_file_dict[filename][readlen][pos]
								if primetype == "threeprime":
									pos += readlen
								offset_pos = pos+offset
								try:
									master_dict[offset_pos] += count
								except Exception as e:
									print "Error tried adding to position {} but tranlen is only {}".format(e,tranlen)
									pass
		elif coverage == True:
			for filename in master_file_dict:
				if filename not in offset_dict:
					continue
				for readlen in master_file_dict[filename]:
					if readlen >= min_read and readlen <= max_read:
						if readlen in offset_dict[filename]:
							offset = offset_dict[filename][readlen]
							for pos in master_file_dict[filename][readlen]:
								count = master_file_dict[filename][readlen][pos]
								if primetype == "threeprime":
									pos += readlen							
								for i in range(0,readlen,3):
									new_offset_pos = (i+pos)+(offset%3)
									try:
										master_dict[new_offset_pos] += count
									except Exception as e:
										#print "Error tried adding to position {} but tranlen is only {}".format(e,tranlen)
										pass
	
		return master_dict, mismatch_dict

# Create dictionary of counts at each position, averged by readlength
def get_readlength_breakdown(read_type, min_read, max_read, tran, user_files,offset_dict,tranlen,coverage, organism, subcodon, noisered, primetype, preprocess, filetype, colorbar_minread,colorbar_maxread):
	master_dict = {}
	color_range = float(colorbar_maxread - colorbar_minread)
	color_list =  all_palettes["RdYlGn"][10]

	for i in range(0,tranlen+max_read):
		master_dict[i] = {}
		for x in range(min_read,max_read+1):
			master_dict[i][x] = 0
	# the keys of master readlen dict are readlengths the value is a dictionary of position:count, there is also a colour key
	colored_master_dict = {}
	master_file_dict = {}
	# first make a master dict consisting of all the read dicts from each filename
	if read_type == "unambig":
		if filetype in user_files:
			for file_id in user_files[filetype]:
				filename = user_files[filetype][file_id]
				if os.path.isfile(filename):
					openshelf = SqliteDict(filename)
				else:
					continue
				if tran in openshelf:
					alltrandict = dict(openshelf[tran])
					unambig_tran_dict = alltrandict["unambig"]
					openshelf.close()
					trandict = unambig_tran_dict
					master_file_dict[filename] = trandict
	else:
		if filetype in user_files:
			for file_id in user_files[filetype]:
				filename = user_files[filetype][file_id]
				if os.path.isfile(filename):
					openshelf = SqliteDict(filename)
				else:
					continue
				if tran in openshelf:
					alltrandict = dict(openshelf[tran])
					openshelf.close()
					unambig_tran_dict = alltrandict["unambig"]
					ambig_tran_dict = alltrandict["ambig"]
					# TODO: Change merge_dicts to take a list of dicts instead of two
					trandict = merge_dicts(unambig_tran_dict, ambig_tran_dict)
					master_file_dict[filename] = trandict

	for filename in master_file_dict:
		for readlen in master_file_dict[filename]:
			if readlen >= min_read and readlen <= max_read:
				if coverage == True:
					for pos in master_file_dict[filename][readlen]:
						count = master_file_dict[filename][readlen][pos]
						for i in range(pos,pos+(readlen+1)):
							master_dict[i][readlen] += count
				else:
					if os.path.isfile(filename):
						openshelf = SqliteDict(filename)
					else:
						continue
					offsets = openshelf["offsets"]["fiveprime"]["offsets"]
					if readlen in offsets:
						offset = offsets[readlen]
					else:
						offset = 15
					for pos in master_file_dict[filename][readlen]:
						count = master_file_dict[filename][readlen][pos]
						if (pos+offset) in master_dict:
							master_dict[pos+offset][readlen] += count

	sorted_master_dict = collections.OrderedDict()
	for key in sorted(master_dict.keys()):
		sorted_master_dict[key] = master_dict[key]

	for pos in sorted_master_dict:
		count = sum(sorted_master_dict[pos].values())
		tot_readlen = 0.0
		tot_count = 0.0001
		for readlen in sorted_master_dict[pos]:
			tot_count += sorted_master_dict[pos][readlen]
			tot_readlen += (sorted_master_dict[pos][readlen]*readlen)
		avg_readlen = int(tot_readlen/tot_count)
		if avg_readlen < colorbar_minread:
			avg_readlen = colorbar_minread
		if avg_readlen > colorbar_maxread:
			avg_readlen = colorbar_maxread
		# find where this avg readlen lies in the range of min readlen to max readlen and use that to assign a color

		y = avg_readlen - colorbar_minread
		per = y/color_range
		final_per = int(per*10)
		if final_per > 9:
			final_per = 9
		color = color_list[final_per]
		if color not in colored_master_dict:
			colored_master_dict[color] = collections.OrderedDict()
			for i in range(0,tranlen+max_read+1):
				colored_master_dict[color][i] = 0
		if pos not in colored_master_dict[color]:
			colored_master_dict[color][pos] = 0
		colored_master_dict[color][pos] += count
	return color_list, colored_master_dict



# Create a dictionary of mismatches at each position
def get_seq_var(user_files, tranlen, organism, tran):
	mismatch_dict = collections.OrderedDict()
	mismatch_dict["A"] = {}
	mismatch_dict["T"] = {}
	mismatch_dict["G"] = {}
	mismatch_dict["C"] = {}
	for i in range(0,tranlen+1):
		mismatch_dict["A"][i] = 0
		mismatch_dict["T"][i] = 0
		mismatch_dict["G"][i] = 0
		mismatch_dict["C"][i] = 0
	for filetype in ["riboseq","rnaseq"]:
		if filetype in user_files:
			for file_id in user_files[filetype]:
				filename = user_files[filetype][file_id]
				if os.path.isfile(filename):
					sqlite_db = SqliteDict(filename, autocommit=False)
				else:
					return "File not found"
				if tran in sqlite_db:
					if "seq" in sqlite_db[tran]:
						sqlite_db_seqvar = dict(sqlite_db[tran]["seq"])
						
						for pos in sqlite_db_seqvar:
							#convert to one based
							fixed_pos = pos+1
							for char in sqlite_db_seqvar[pos]:
								if char != "N":
									if fixed_pos not in mismatch_dict[char]:
										mismatch_dict[char][fixed_pos] = 0
									count = sqlite_db_seqvar[pos][char]
									mismatch_dict[char][fixed_pos] += count

					sqlite_db.close()


	return mismatch_dict




