import sqlite3
from sqlitedict import SqliteDict


def get_gene_info(gene_id, sqlite_path_organism):
    # get the transcript ID, exon junctions and transcript sequence from the specified sqlite file for a given gene ID
    gene_id = gene_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript, exon_junctions, sequence FROM transcripts WHERE gene = (:gene);', {'gene': gene_id})
    gene_info = c.fetchall()
    return gene_info


def get_transcript_length(gene, sqlite_path_organism):
    # get the transcript length from the specified sqlite file for a given gene ID
    gene = gene.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript, length FROM transcripts WHERE gene = (:gene);', {'gene': gene})
    transcript_length = c.fetchall()
    length = {transcript_len[0]:transcript_len[1] for transcript_len in transcript_length} 
    return length


def get_genes_principal(gene_id, sqlite_path_organism):
    # get the transcript ID for the given genes principal isoform as per APPRIS
    gene_id = gene_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript FROM transcripts '
              'WHERE gene = (:gene) AND principal = 1;', {'gene': gene_id})
    principal_transcript = c.fetchall()
    return principal_transcript[0][0]


def string_to_integer_list(lst):
    # Return the list of integers from a list of strings. If the list is empty then return 0
    if lst[0] == "":
        lst[0] = "0"
    new_lst = []
    for i in range(len(lst)):
        new_lst.append(int(lst[i]))
    return new_lst


def get_3prime_exon(junction_list, sequence):
    # Part of the process of producing the sequences of all exons in a transcript.
    # this function slices the 3' exon sequence from the transcript sequence returning
    # both the exon sequence and the remaining sequence of the transcript
    cut_off = junction_list[-1]
    exon = sequence[cut_off - 1:]
    seq_less_exon = sequence[:cut_off - 1]
    return exon, seq_less_exon


def get_transcript_info(transcript_id, sqlite_path_organism):
    # Returns the exon junctions and sequences for the given transcript
    transcript_id = transcript_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT exon_junctions, sequence FROM transcripts WHERE transcript = (:trans);', {'trans': transcript_id})
    transcript_info = c.fetchall()
    return transcript_info



def get_exon_coordinates_for_orf(transcript_id, sqlite_path_organism):
    # return the coordinates of the exons in a given transcript. 0 -> first junction, first junction -> second junction
    # last junction -> end
    transcript_info = get_transcript_info(transcript_id, sqlite_path_organism)
    exon_junctions = string_to_integer_list(str(transcript_info[0][0]).split(","))
    transcript_end = len(transcript_info[0][1])
    exon_coordinates = []
    for junction in exon_junctions:
        if exon_junctions.index(junction) == 0:
            exon_coordinates.append((0, junction))
        else:
            exon_coordinates.append((exon_junctions[exon_junctions.index(junction) - 1] + 1, junction))

    exon_coordinates.append((exon_junctions[-1] + 1, transcript_end))

    return exon_coordinates


def get_protein_coding_transcript_ids(gene, sqlite_path_organism):
    # get transcript IDs for protein coding genes for the given gene
    # returns a list of strings 
    gene_id = gene.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()

    c.execute('SELECT transcript, tran_type FROM transcripts WHERE gene = (:gene) AND tran_type = 1;', {'gene': gene_id})
    coding_info = c.fetchall()
    protein_coding_transcripts = [str(i[0]) for i in coding_info]
    return protein_coding_transcripts


def get_start_stop_codon_positions(transcript_id, sqlite_path_organism):
    # Get start and stop codon positions from annotation sqlite
    transcript_id = transcript_id.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT cds_start, cds_stop FROM transcripts WHERE transcript = (:transcript);', {'transcript': transcript_id})
    a = c.fetchone()
    return a


def get_orf_exon_structure(start_stop, exon_coordinates):
    # determine the stucture of each ORF. Returns coordinates in the form:
    # Initiation site to junction, junction to junction, junction to translation stop
    start, stop = start_stop

    started = False
    structure = []
    for exon in exon_coordinates:
        if start in range(exon[0], exon[1]):
            started = True
            structure.append((start, exon[1]))

        elif started and (stop in range(exon[0], exon[1])):
            structure.append((exon[0], stop))
        elif started:
            structure.append(exon)
    return structure


def get_coding_regions_for_genes_transcripts(gene, sqlite_path_organism):
    # Returns the coding region coordinates annotated on each transcript for a given gene
    gene_id = gene.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT * FROM coding_regions WHERE transcript IN (SELECT transcript FROM transcripts WHERE gene = (:gene));', {'gene': gene_id})
    coding_regions = c.fetchall()
    return coding_regions


def get_max_min_exon_genomic_positions(exon_info):
    # Of all exons starts return the lowest position and highest position of all stops
    starts = []
    for i, value in enumerate(exon_info):
        starts.append(exon_info[i][1])

    stops = []
    for i, value in enumerate(exon_info):
        stops.append(exon_info[i][2])
    return min(starts), max(stops)


def exons_of_transcript(transcript_id, sqlite_path_organism):
    # For a given transcript return the exon sequences in a list in 5' to 3' direction
    exon_lst = []
    trans_info = get_transcript_info(transcript_id, sqlite_path_organism)
    exon_junctions = trans_info[0][0].split(",")
    sequence = trans_info[0][1]

    exon_junct_int = string_to_integer_list(exon_junctions)

    while len(exon_junct_int) != 0:
        exon, sequence = get_3prime_exon(exon_junct_int, sequence)
        exon_lst.append(exon)

        exon_junct_int.pop(-1)
    exon_lst.append(sequence)
    return exon_lst[::-1]


def get_exon_coordinate_ranges(sequence, junctions):
    # Return list of transcript coordinate ranges (start position, stop position).
    junctions = [0] + junctions + [len(sequence) - 1]
    ranges = [(int(junctions[i]), int(junctions[(i + 1) % len(junctions)]))  for i in range(len(junctions)-1)]

    return ranges


def transcript_exon_coordinate_ranges(gene_info):
    # return dictionary of transcript exon coordinates for a given locus
    transcript_exon_coordinate_ranges = {} 
    for transcript in gene_info:
        if transcript[0] not in transcript_exon_coordinate_ranges:
            transcript_exon_coordinate_ranges[transcript[0]] = get_exon_coordinate_ranges(transcript[2], str(transcript[1]).split(","))
    return  transcript_exon_coordinate_ranges


def get_exon_info(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # Return the exon starts,stops and transcript for a given gene
    gene = gene.upper()
    conn = sqlite3.connect(sqlite_path_organism)
    c = conn.cursor()
    c.execute('SELECT transcript, exon_start, exon_stop FROM exons WHERE transcript IN '
              '(SELECT transcript FROM transcripts WHERE gene = (:gene));', {'gene': gene})
    exon_info = c.fetchall()
    if filter:
        supported_exon_info = []
        for exon in exon_info:
            if exon[0] in supported_transcripts:
                supported_exon_info.append(exon)
        return supported_exon_info
    return exon_info


def genomic_exon_coordinate_ranges(gene, sqlite_path_organism, supported_transcripts, filter=True):
    # create a dictionary with transcript_ids as keys and exon coordinates in tuple (start, stop) as values
    # subtract the minumum start codon position of any exon for the gene
    exon_info = get_exon_info(gene, sqlite_path_organism, supported_transcripts, filter=filter)
    minimum, _ = get_max_min_exon_genomic_positions(exon_info)
    genomic_coordinates_per_transcript = {}
    for exon in exon_info:
        if exon[0] not in genomic_coordinates_per_transcript:
            genomic_coordinates_per_transcript[exon[0]] = [(exon[1] - minimum, exon[2] - minimum)]
        else:
            genomic_coordinates_per_transcript[exon[0]].append((exon[1] - minimum, exon[2] - minimum))
    return genomic_coordinates_per_transcript



def genomic_orf_coordinate_ranges(gene, sqlite_path_organism, supported_transcripts, genomic_coordinates, filter=True):
    # Translate the orf structure to genomic coordinates using genomic coordinates from the exons table of sqlite file
    orf_structures = {}
    genomic_orf_structures = {}
    for transcript in supported_transcripts:
        start_stop = get_start_stop_codon_positions(transcript, sqlite_path_organism)
        transcript_coordinates = get_exon_coordinates_for_orf(transcript, sqlite_path_organism)
        orf_structures[transcript] = get_orf_exon_structure(start_stop, transcript_coordinates)

        for region in orf_structures[transcript]:
            for exon in transcript_coordinates:
                if region[1] == exon[1]:
                    difference = region[1] - region[0]
                    genomic_exon = genomic_coordinates[transcript][transcript_coordinates.index(exon)]
                    genomic_orf = (genomic_exon[1] - difference, genomic_exon[1])
                    if transcript in genomic_orf_structures:
                        genomic_orf_structures[transcript].append(genomic_orf)
                    else:
                        genomic_orf_structures[transcript] = [genomic_orf]

                elif region[0] == exon[0]:
                    difference = region[1] - region[0]

                    genomic_exon = genomic_coordinates[transcript][transcript_coordinates.index(exon)]
                    genomic_orf = (genomic_exon[0], genomic_exon[0] + difference)
                    if transcript in genomic_orf_structures:
                        genomic_orf_structures[transcript].append(genomic_orf)
                    else:
                        genomic_orf_structures[transcript] = [genomic_orf]

    return genomic_orf_structures


def genomic_junction_positions(gene, sqlite_path_organism, supported_transcripts, exons, filter=True):
    genomic_junctions = {}
    for transcript in exons:
        if transcript not in genomic_junctions:
            genomic_junctions[transcript] = []

        number_of_exons = len(exons[transcript])
        for index, exon in enumerate(exons[transcript]):
            if index < number_of_exons - 1:
                genomic_junctions[transcript].append((exons[transcript][index][1], exons[transcript][index + 1][0]))
    return genomic_junctions


def genomic_junction_scores(gene_name, sqlite_path_organism, sqlite_path_reads, supported_transcripts, exons, genomic_junctions, filter=True):

    scores = get_scores_per_exonjunction_for_gene(sqlite_path_organism, sqlite_path_reads, supported_transcripts)
    paired_dictionary = {}
    for transcript in genomic_junctions:
        sorted_scores = sorted(scores[transcript].keys())
        paired = zip(genomic_junctions[transcript], sorted_scores)

        paired_dictionary[transcript] = paired
    genomic_junctions_scores = {}
    for transcript in paired_dictionary:
        if transcript not in genomic_junctions_scores:
            genomic_junctions_scores[transcript] = {}

        for junction in paired_dictionary[transcript]:
            genomic_junctions_scores[transcript][junction[0]] = scores[transcript][junction[1]]

    return genomic_junctions_scores


def get_reads_per_transcript_location(transcript_id, sqlite_path_reads):
    infile = SqliteDict(sqlite_path_reads)
    if transcript_id not in infile:
        print("No unambiguous reads support this gene", transcript_id)
        return None
    reads = infile[transcript_id]["unambig"]
    return reads


def get_transcript_genomic_exon_dict(genomic, transcriptomic):
    genomic_transcriptomic_exon_dict = {}
    for transcript in transcriptomic: 
        t_list = sorted(transcriptomic[transcript])
        g_list = sorted(genomic[transcript])
        for pair in zip(t_list, g_list):
            if pair[0] not in genomic_transcriptomic_exon_dict:
                genomic_transcriptomic_exon_dict[pair[0]] = pair[1]

    return genomic_transcriptomic_exon_dict


def get_reads_per_genomic_location_asite(gene, sqlite_path_reads, sqlite_path_organism, supported_transcripts, genomic_exon_coordinates, filter=True):
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.
    gene_info = get_gene_info(gene, sqlite_path_organism)
    genomic_read_dictionary = {}
    for read_file in sqlite_path_reads:

        infile = SqliteDict(read_file)
        for transcript in gene_info:
            counted_reads = []

            if filter:
                if transcript[0] not in supported_transcripts:
                    print 'filtered'
                    continue
            if transcript[0] not in infile:
                print("No unambiguous reads support this transcript", transcript[0])
                continue
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = sorted(genomic_exon_coordinates[transcript[0]])

            exon_junctions = string_to_integer_list(transcript[1].split(","))

            sequence = transcript[2]
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(sequence, exon_junctions)

            for length in transcript_read_dictionary:
                for location in transcript_read_dictionary[length]:
                    position = location + infile["offsets"]["fiveprime"]["offsets"][length]

                    range_counter = 0

                    for exon in transcript_ranges:
                        if position in range(exon[0], exon[1]):
                            difference_between_read_position_and_exon_asite = position - exon[0]

                            genomic_asite = genomic_ranges[range_counter][
                                                0] + difference_between_read_position_and_exon_asite
                                    
                            if (length, genomic_asite) not in counted_reads:
                                if genomic_asite not in genomic_read_dictionary:
                                    genomic_read_dictionary[genomic_asite] = transcript_read_dictionary[length][location]
                                else:
                                    genomic_read_dictionary[genomic_asite] += transcript_read_dictionary[length][location]
                                counted_reads.append((length, genomic_asite))
                        range_counter += 1
    # print genomic_read_dictionary
    return genomic_read_dictionary


def get_reads_per_genomic_location_asite_faster(gene, sqlite_path_reads, sqlite_path_organism, supported_transcripts, genomic_exon_coordinates, transcript_exon_coordinate_ranges,filter=True):
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.
    gene_info = get_gene_info(gene, sqlite_path_organism)
    genomic_read_dictionary = {}
    for read_file in sqlite_path_reads:
        counted_exons = {}
        gene_info.sort(key=lambda a: len(a[2]))

        infile = SqliteDict(read_file)
        # for i in infile.keys(): 
        #     if not i.isupper():
        #         print i
        for transcript in gene_info:
            counted_reads = []

            if filter:
                if transcript[0] not in supported_transcripts:
                    print 'filtered'
                    continue
            if transcript[0] not in infile:
                print("No unambiguous reads support this transcript", transcript[0])
                continue
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = sorted(genomic_exon_coordinates[transcript[0]])
            genomic_transcriptomic_exon_dict = get_transcript_genomic_exon_dict(genomic_exon_coordinates, transcript_exon_coordinate_ranges)

            exon_junctions = string_to_integer_list(transcript[1].split(","))
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(transcript[2], exon_junctions)


            for transcript_exon in transcript_ranges:
                exon_count_total = 0
                if genomic_transcriptomic_exon_dict[transcript_exon] not in counted_exons:
                    for exon in counted_exons:
                        genomic_start = genomic_transcriptomic_exon_dict[transcript_exon][0]
                        genomic_stop = genomic_transcriptomic_exon_dict[transcript_exon][1]
                        if exon[0] in range(genomic_start, genomic_stop)  and exon[1] in range(genomic_start, genomic_stop):
                            print "Within ranges: ", exon, genomic_transcriptomic_exon_dict[transcript_exon]
                    for length in transcript_read_dictionary:
                        sorted_read_positions = sorted(transcript_read_dictionary[length])
                        relevant_reads = [i for i in sorted_read_positions if i in range(transcript_exon[0], transcript_exon[1])]
                        for location in relevant_reads:
                            position = location + infile["offsets"]["fiveprime"]["offsets"][length]
                            difference_between_read_position_and_exon_asite = position - transcript_exon[0]


                            genomic_asite = genomic_transcriptomic_exon_dict[transcript_exon][
                                                0] + difference_between_read_position_and_exon_asite

                            if (length, genomic_asite) not in counted_reads:
                                if genomic_asite not in genomic_read_dictionary:
                                    genomic_read_dictionary[genomic_asite] = transcript_read_dictionary[length][location]

                                else:
                                    genomic_read_dictionary[genomic_asite] += transcript_read_dictionary[length][location]

                                counted_reads.append((length, genomic_asite))
                                exon_count_total += transcript_read_dictionary[length][location]

                    counted_exons[genomic_transcriptomic_exon_dict[transcript_exon]] = exon_count_total
                else: 
                    a = 1
                    # print counted_exons
                    # print genomic_transcriptomic_exon_dict[transcript_exon], counted_exons[genomic_transcriptomic_exon_dict[transcript_exon]]

    return genomic_read_dictionary

def get_reads_per_genomic_location_fiveprime(gene, sqlite_path_reads, sqlite_path_organism, supported_transcripts, genomic_exon_coordinates,
                                             filter):
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.

    gene_info = get_gene_info(gene, sqlite_path_organism)
    genomic_read_dictionary = {}
    counted_reads = []
    for read_file in sqlite_path_reads:
        infile = SqliteDict(read_file)

        for transcript in gene_info:
            if filter:
                if transcript[0] not in supported_transcripts:
                    continue
            if transcript[0] not in infile:
                print("No unambiguous reads support this gene")
                continue
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = genomic_exon_coordinates[transcript[0]]

            exon_junctions = string_to_integer_list(transcript[1].split(","))
            sequence = transcript[2]
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(sequence, exon_junctions)

            for length in transcript_read_dictionary:
                for position in transcript_read_dictionary[length]:

                    range_counter = 0
                    for exon in transcript_ranges:
                        if position in range(exon[0], exon[1]):
                            difference_between_read_position_and_exon_start = position - exon[0]
                            genomic_start_pos = genomic_ranges[range_counter][
                                                    0] + difference_between_read_position_and_exon_start
                            if (length, genomic_start_pos) not in counted_reads:
                                if genomic_start_pos not in genomic_read_dictionary:
                                    genomic_read_dictionary[genomic_start_pos] = transcript_read_dictionary[length][
                                        position]
                                else:
                                    genomic_read_dictionary[genomic_start_pos] += transcript_read_dictionary[length][
                                        position]

                                counted_reads.append((length, genomic_start_pos))
                        range_counter += 1
    return genomic_read_dictionary


def get_read_ranges_genomic_location(gene, sqlite_path_reads, sqlite_path_organism, supported_transcripts, genomic_exon_coordinates, filter):
    # get the number of reads supporting each genomic position to be used in the display of support of the
    # supertranscript model. This function takes the reads mapped for each transcript of a gene and uses a combination
    # of genomic and transcriptomic ranges to translate each transcript position to a genomic one.

    gene_info = get_gene_info(gene, sqlite_path_organism)
    genomic_read_dictionary = {}

    for read_file in sqlite_path_reads:
        infile = SqliteDict(read_file)

        for transcript in gene_info:
            if filter:
                if transcript[0] not in supported_transcripts:
                    continue
            if transcript[0] not in infile:
                print("No unambiguous reads support this gene")
                return None
            transcript_read_dictionary = infile[transcript[0]]["unambig"]
            genomic_ranges = genomic_exon_coordinates[transcript[0]]

            exon_junctions = string_to_integer_list(transcript[1].split(","))
            sequence = transcript[2]
            exons = exons_of_transcript(transcript[0], sqlite_path_organism)
            transcript_ranges = get_exon_coordinate_ranges(sequence, exon_junctions)

            for length in transcript_read_dictionary:
                for position in transcript_read_dictionary[length]:
                    range_counter = 0
                    for exon in transcript_ranges:
                        if position in range(exon[0], exon[1]):
                            difference_between_read_position_and_exon_start = position - exon[0]
                            genomic_start_pos = genomic_ranges[range_counter][
                                                    0] + difference_between_read_position_and_exon_start
                            genomic_read_range = (genomic_start_pos, genomic_start_pos + length)

                            if genomic_read_range not in genomic_read_dictionary:
                                genomic_read_dictionary[genomic_read_range] = transcript_read_dictionary[length][
                                    position]
                            else:
                                genomic_read_dictionary[genomic_read_range] += transcript_read_dictionary[length][
                                    position]
                        range_counter += 1

    return genomic_read_dictionary


def get_exonjunction_pileup_for_transcript(transcript_id, sqlite_path_organism, sqlite_path_reads):
    # count the number of reads in the read file that span each exon-exon junction. for a given transcript
    # returns a dictionary with junctions as keys and counts as values d
    transcript_info = get_transcript_info(transcript_id, sqlite_path_organism)
    exon_junctions = string_to_integer_list(transcript_info[0][0].split(","))
    counts = {}

    for read_file in sqlite_path_reads:
        reads = get_reads_per_transcript_location(transcript_id, read_file)
        if reads == None:
            continue
        for junction in exon_junctions:
            for read_length in reads:
                for position in reads[read_length]:
                    if junction in range(position, position + read_length):

                        if junction in counts:
                            counts[junction] += reads[read_length][position]
                        else:
                            counts[junction] = reads[read_length][position]
                    else:
                        if junction in counts:
                            counts[junction] += 0
                        else:
                            counts[junction] = 0
    return counts


def get_scores_per_exonjunction_for_gene(sqlite_path_organism, sqlite_path_reads, supported):
    # count the reads in the reads file whos p sites lite within the exon sequence
    # returns a dictionary with all unique exons in the gene as keys and counts as values
    pileup = {}
    for trans in supported:
        pileup[trans] = get_exonjunction_pileup_for_transcript(trans, sqlite_path_organism, sqlite_path_reads)

    return pileup


