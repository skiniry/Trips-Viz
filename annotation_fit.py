from sqlitedict import SqliteDict
import sqlite3
import bisect 
import collections
import time

def get_annotation_connection(annotation_path):
    '''
    Opne a sqlite connection with file at provided path. Returns connection
    '''
    conn = sqlite3.connect(annotation_path)
    return conn


def check_overlap(region1, region2):
    
    combined_list = []
    if region1[0] in range(region2[0], region2[1]) and region1[1] in range(region2[0], region2[1]):
        combined_list.append(region2)

    elif region2[0] in range(region1[0], region1[1]) and region2[1] in range(region1[0], region1[1]):
        combined_list.append(region1)

    elif region1[0] in range(region2[0], region2[1]) and region1[1] not in range(region2[0], region2[1]):
        combined_list.append((region2[0], region1[1]))

    elif region1[0] not in range(region2[0], region2[1]) and region1[1] in range(region2[0], region2[1]):
        combined_list.append((region1[0], region2[1]))

    else:
        combined_list.append(region2)
        combined_list.append(region1)


    return combined_list


def merge_list(coding_list):
    coding_list.sort()
    number_of_comparisons = len(coding_list) - 1
    comparison_count = 0 
    while comparison_count < number_of_comparisons - 1:
        combined_list = check_overlap(coding_list[comparison_count], coding_list[comparison_count + 1])
        second = coding_list[comparison_count + 1]
        if len(combined_list) == 1:
            coding_list.remove(coding_list[comparison_count])
            coding_list.remove(coding_list[comparison_count])
            coding_list.append(combined_list[0])
            coding_list.sort()

        else: 
            comparison_count += 1 
        if second == coding_list[-1]:
            comparison_count = number_of_comparisons 
    return coding_list


def get_merged_cds_coordinates(transcript, cds, connection):
    curs = connection.cursor()
    result = curs.execute("SELECT * FROM coding_regions WHERE transcript = ?;", (transcript,)).fetchall()
    coding_regions = []
    for entry in result:
        coding_regions.append((entry[1], entry[2]))
    coding_regions.append(cds)

    if len(coding_regions) == 1:
        return coding_regions
    elif len(coding_regions) == 2:
        new_coding_regions = check_overlap(coding_regions[0], coding_regions[1])
    else:
        new_coding_regions = merge_list(coding_regions)

    return new_coding_regions


def get_gene_info(gene_id, connection):
    '''
    Return a dictionary of gene info for a given gene_id using sqlite row factory.
    keys are the column names of the transcript sqlite table 
    keys = [
    transcript, gene, length, cds_start, cds_stop, sequence, strand, stop_list, start_list, 
    exon_junctions, tran_type, gene_type, principal,version,gc,five_gc, cds_gc, three_gc, chrom
    ]
    '''
    gene_info_dict = {}
    gene_id = gene_id.upper()
    connection.row_factory = sqlite3.Row 
    curs = connection.cursor()
    gene_info = curs.execute('SELECT * FROM transcripts WHERE gene = (:gene);', {'gene': gene_id}).fetchall()
    return gene_info


def get_gene_from_transcript(transcript, connection):
    '''
    return the gene name of a given transcript
    '''
    curs = connection.cursor()
    gene = curs.execute("SELECT gene FROM transcripts WHERE transcript = (:transcript);", {'transcript':transcript}).fetchall()[0][0]
    return gene


def transcript_architecture(gene_info):
    '''
    Input gene info dict. Return a nested dictionary with transcripts as keys of outer dict.
    Inner dict contains coordinates of 5' UTR CDS and 3' UTR
    '''
    architecture = {}
    for row in gene_info:
        if row['cds_start'] > 0 and row['cds_stop'] > 0:
            architecture[row['transcript']] = {
                'five_prime' : (0, row['cds_start']), 
                'coding':(row['cds_start'], row['cds_stop']),
                'three_prime': (row['cds_stop'] + 1 , row['length'])
                }
        else:
            architecture[row['transcript']] = {
                "noncoding":(0, row['length'])
            }
    
    return architecture


def get_reads_per_transcript_location(transcript_id, sqlite_path_reads):
    '''
    get the unambiguously mapped reads per transcript. Returns nest dictionary with lengths as keys and dictionary of positions:counts format 
    '''
    infile = SqliteDict(sqlite_path_reads)
    if transcript_id not in infile:
        print("No unambiguous reads support this gene", transcript_id)
        return None
    reads = infile[transcript_id]["unambig"]
    return reads


def get_offsets_for_read_file(sqlite_path_reads):
    '''
    get the five prime offset dict from the read file to fascilitate the calculation of a site predictions 

    '''
    infile = SqliteDict(sqlite_path_reads)
    return infile['offsets']['fiveprime']['offsets']


def get_index_for_range(positions, region):
    '''
    return the indexes to slice a list to contain the elements relevant to given region 
    takes sorted list and tuple of 2 integers 
    '''
    if region[0] in positions:
        region_start = positions.index(region[0])
    
    else: 
        bisect.insort(positions, region[0]) 
        region_start = positions.index(region[0])
        positions.pop(region_start)
    
    if region[1] in positions:
        region_stop = positions.index(region[1])
    
    else: 
        bisect.insort(positions, region[1]) 
        region_stop = positions.index(region[1])
        positions.pop(region_stop)

    return region_start, region_stop


def get_counts_per_region(trans_architecture, transcript_counts, offsets, coding_regions):
    '''
    calculate the asite counts per region of each transcript 

    '''
    count_dict = {'merged':0, 'five_prime':0, 'three_prime':0, 'coding':0}

    for length in transcript_counts:
        positions = sorted(transcript_counts[length].keys())
        
        for region in trans_architecture:
            updated_region = (trans_architecture[region][0] - offsets[length], trans_architecture[region][1] - offsets[length])
            region_start, region_stop = get_index_for_range(positions, updated_region)
            relevant_positions = positions[region_start:region_stop]
            read_dict = {k:transcript_counts[length][k] for k in relevant_positions if k in relevant_positions}
            count_dict[region] += sum(read_dict.values())
        
        for region in coding_regions:
            updated_region = (region[0] - offsets[length], region[1] - offsets[length])
            region_start, region_stop = get_index_for_range(positions, updated_region)
            relevant_positions = positions[region_start:region_stop]
            read_dict = {k:transcript_counts[length][k] for k in relevant_positions if k in relevant_positions}
            count_dict['merged'] += sum(read_dict.values())
    return count_dict

def get_range(d, begin, end):
    return {i:d[i] for i in d.keys() if begin <= i <= end}


def annotation_fit(reads, transcript, transhelve):
    '''
    main function that gets regions within transcript and counts the reads that map to those regions 
    Returns a dict with keys region_counts, region_lengths, merged_cds, total_reads which can be used to calcualte annotation fit
    reads: product of get_reads[0] 
    transcript: transcript ID 
    transhelve: connection to annotation dict 

    '''
    gene = get_gene_from_transcript(transcript, transhelve)
    gene_info = get_gene_info(gene, transhelve)
    transcript_architecture_dict = transcript_architecture(gene_info)

    annotation_fit = {}
    for region in transcript_architecture_dict[transcript]:
        reads_in_region = get_range(dict(reads[0]), transcript_architecture_dict[transcript][region][0], transcript_architecture_dict[transcript][region][1])
        length =  region + "_length"
        counts = region + "_counts"
        annotation_fit[counts] = sum(reads_in_region.values())
        annotation_fit[length] = transcript_architecture_dict[transcript][region][1] - transcript_architecture_dict[transcript][region][0]
        if region == 'coding':
            reads_in_frame = {i:reads_in_region[i] for i in reads_in_region.keys() if (i % 3) == (transcript_architecture_dict[transcript][region][0] % 3)}
            annotation_fit["in_frame_counts"] = sum(reads_in_frame.values())
            number_of_codons = float(transcript_architecture_dict[transcript][region][1] - transcript_architecture_dict[transcript][region][0])/3
            covered_codons = [position for position in reads_in_frame if reads_in_frame[position] >= 1]
            annotation_fit['in_frame_coverage'] = len(covered_codons)/number_of_codons

    coding_regions = get_merged_cds_coordinates(transcript, transcript_architecture_dict[transcript]['coding'], transhelve)

    merged_cds_counts = 0
    merged_length = 0
    for region in coding_regions:
        merged_length += region[1] - region[0]
        reads_in_region = get_range(dict(reads[0]), region[0], region[1])
        merged_cds_counts += sum(reads_in_region.values())


    annotation_fit["merged_counts"] = merged_cds_counts
    annotation_fit["merged_length"] = merged_length
    annotation_fit["total_reads"] = sum(reads[0].values())
    return annotation_fit

