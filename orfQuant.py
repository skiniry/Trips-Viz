from tripsSplice import genomic_exon_coordinate_ranges
from tripsSplice import genomic_orf_coordinate_ranges
from tripsSplice import genomic_junction_positions
from tripsSplice import get_gene_info
from sqlitedict import SqliteDict


from tripsSplice import get_protein_coding_transcript_ids
from tripsSplice import get_reads_per_genomic_location_asite
from tripsSplice import get_reads_per_genomic_location_asite_faster

from tripsSplice import genomic_junction_scores
from tripsSplice import get_start_stop_codon_positions
from tripsSplice import transcript_exon_coordinate_ranges
import tripsSplice

from tripsCount import count_read_supporting_regions_per_transcript


def classify_regions_shared_unique(regions):
# Takes in dictionary of genomic coordinates that make up the features for each transcript
# for each region consider each transcript and count the number of occurences of this exact region in that transcript
# If this count is > 1 it is considered not unique. Otherwise unique 

    classified = {}
    for tran1 in regions:
        if tran1 not in classified:
            classified[tran1] = {"unique":[], "shared":[]}
        for region in regions[tran1]:
            counter = 0
            for tran2 in regions:
                counter += regions[tran2].count(region)

            if counter > 1:
                classified[tran1]["shared"].append(region)
            else:
                classified[tran1]["unique"].append(region)
    return classified


def get_orf_coordinate_junctions(orf_coordinates):
# Get the coordinates of junctions within ORFs. Returns a dictionary with transcripts as keys and a list of jucntions
# as values. Essentially reformats orf_coordinates to be junction oriented. 
    junctions = {}
    for transcript in orf_coordinates:
        if transcript not in junctions:
            junctions[transcript] = []
        for i in range(len(orf_coordinates[transcript]) - 1):
            junctions[transcript].append((orf_coordinates[transcript][i][1], orf_coordinates[transcript][i + 1][0] - 1))

    return junctions


def region_coverage(regions, counts):
# Returns a dictioanry of coverage of each coding region. Calculated as counts/length rounded to 4 places
    region_counts = {}
    for transcript in regions:
        region_counts_transcript = zip(regions[transcript], counts[transcript])
        if transcript not in region_counts:
            region_counts[transcript] = region_counts_transcript

    coverage_region_transcript = {}
    for transcript in region_counts:
        if transcript not in coverage_region_transcript:
            coverage_region_transcript[transcript] = {}

        for item in region_counts[transcript]:
            length = item[0][1] - item[0][0]
            coverage = round(float(item[1])/float(length), 4)
            coverage_region_transcript[transcript][item[0]] = coverage

    return coverage_region_transcript


def coverage_junction_transcript(junction_scores):
    #function to calculate the normalised coverage of junction features. returns a nested dict.
    # 60 = 30 + 30 which is the sum of two riboseq read lengths.
    coverage_junction = {}
    for transcript in junction_scores:
        if transcript not in coverage_junction:
            coverage_junction[transcript] = {}
        for junction in junction_scores[transcript]:
            coverage = float(junction_scores[transcript][junction]) / float(60)
            coverage_junction[transcript][junction] = coverage

    return coverage_junction


def average_unique_coverage(unique_shared_exons, unique_shared_junctions, coverage_junction, coverage_exons):
# Return the average coverage of all unique features. Sum of coverages over number of unique features

    average_features = {}

    for transcript in unique_shared_exons:
        sum = 0
        count = 0
        for item in unique_shared_exons[transcript]['unique']:
            sum += coverage_exons[transcript][item]
            count += 1

        for junction in unique_shared_junctions[transcript]['unique']:
            sum += coverage_junction[transcript][junction]
            count += 1

        if transcript not in average_features:
            if unique_shared_exons[transcript]['unique'] == [] and unique_shared_junctions[transcript]['unique'] == []:
                average_features[transcript] = 0
            else:
                average_features[transcript] = float(sum)/float(count)

    return average_features


def all_feature_average(coverage_junction, coverage_exons):
# Average of the coverage of all features. Sum of coverages over number of features
    average_features = {}

    for transcript in coverage_exons:
        sum = 0
        count = 0
        for exon in coverage_exons[transcript]:
            sum += coverage_exons[transcript][exon]
            count += 1

        for junction in coverage_junction[transcript]:
            sum += coverage_junction[transcript][junction]
            count += 1


        if transcript not in average_features:
            average_features[transcript] = float(sum) / float(count)
    return average_features


def cORF_ratio(average_unique, average_all):
# Calculate the ratio of average unique coverages to average all coverages
    cORF ={}
    for transcript in average_all:
        cORF[transcript] = average_unique[transcript] / average_all[transcript]

    return cORF


def adjusted_coverage_for_non_unique_orfs(transcript, coverage_exons, coverage_junctions, average_coverage_all):
# Calculate an adjusted value for cases where a transcript has no unique features
# This is coverage of exons - coverage of exons * sum of all coverages from that transcript
# Average is calculated of these adjusted values 

    adjusted_coverage_exons = {transcript: {}}
    for exon in coverage_exons[transcript]:
        signal = 0
        for trans in coverage_exons:
            if exon in coverage_exons[trans]:
                signal += coverage_exons[trans][exon]

        if exon not in adjusted_coverage_exons[transcript]:
            adjusted_coverage_exons[transcript][exon] = coverage_exons[transcript][exon] - (
                    coverage_exons[transcript][exon]  * signal)

    adjusted_coverage_junctions = {transcript: {}}
    for junction in coverage_junctions[transcript]:
        signal = 0
        for trans in coverage_junctions:
            if junction in coverage_junctions[trans]:
                signal += coverage_junctions[trans][junction]

        if junction not in adjusted_coverage_junctions[transcript]:
            adjusted_coverage_junctions[transcript][junction] = coverage_junctions[transcript][junction] - (
                        coverage_junctions[transcript][junction] * signal)
    sum = 0
    count = 0
    for region in adjusted_coverage_exons[transcript]:
        sum += adjusted_coverage_exons[transcript][region]
        count += 1

    for junction in adjusted_coverage_junctions[transcript]:
        sum += adjusted_coverage_junctions[transcript][junction]
        count += 1
    average_adjusted_coverage = float(sum)/float(count)

    cORF_transcript = average_adjusted_coverage / average_coverage_all[transcript]
    return cORF_transcript


def shared_coverage(coverage_exons, coverage_junctions):
# Calculate an adjusted average value where all features are shared between transcripts. 


    shared_exon_coverage = {}
    for transcript in coverage_exons:
        if transcript not in shared_exon_coverage:
            shared_exon_coverage[transcript] = {}
        for exon in coverage_exons[transcript]:
            number_ORF_over_F = 0
            for tran in coverage_exons:
                for exon2 in coverage_exons[tran]:
                    if transcript == tran:
                        continue

                    if exon == exon2:
                        number_ORF_over_F += 1
            number_ORF_over_F += 1

            shared_exon_coverage[transcript][exon] = coverage_exons[transcript][exon]/number_ORF_over_F

    shared_junction_coverage = {}
    for transcript in coverage_junctions:
        if transcript not in shared_junction_coverage:
            shared_junction_coverage[transcript] = {}
        for junction in coverage_junctions[transcript]:
            number_ORF_over_F = 0
            for tran in coverage_junctions:
                for junction2 in coverage_junctions[tran]:
                    if transcript == tran:
                        continue

                    if junction == junction2:
                        number_ORF_over_F += 1
            number_ORF_over_F += 1
            shared_junction_coverage[transcript][junction] = coverage_junctions[transcript][junction] / number_ORF_over_F

    shared_cORF = {}
    for transcript in shared_exon_coverage:
        sum = 0
        count = 0
        for exon in shared_exon_coverage[transcript]:
            sum += shared_exon_coverage[transcript][exon]
            count += 1

        for junction in shared_junction_coverage[transcript]:
            sum += shared_junction_coverage[transcript][junction]
            count += 1

        if transcript not in shared_cORF:
            shared_cORF[transcript] = float(sum)/float(count)

    return shared_cORF


def aORF(cORF, counts):
# Returns a dictionary of the number of a sites found in each transcript. Coverage by the sum of counts.
    a_sites_perORF = {}
    for transcript in cORF:
        if transcript not in a_sites_perORF:
            a_sites_perORF[transcript] = cORF[transcript] * sum(counts[transcript])

    return a_sites_perORF


def lORF(coding, sqlite_path_organism):
# Returns the lengths of each open reading frame in a dictionary. Transcript as keys, lengths as values
    lORFs = {}

    for transcript in coding:
        start, stop = get_start_stop_codon_positions(transcript, sqlite_path_organism)
        length = (stop + 1) - start
        if transcript not in lORFs:
            lORFs[transcript] = length
    return lORFs


def ORFs_per_million(aORF, lORF):
# Return the TPM like value OPM. Normalise coverage by length of each orf.
# One million divided by scaling factor is the sum of all normalised coverages 

    aORF_lORF = {}
    full_set = 0
    for transcript in aORF:
        if transcript not in aORF_lORF:
            aORF_lORF[transcript] = aORF[transcript]/lORF[transcript]
            full_set += aORF[transcript]/lORF[transcript]

    orfs_per_million = {}
    for transcript in aORF_lORF:
        if transcript not in orfs_per_million:
            orfs_per_million[transcript] = aORF_lORF[transcript] * (10**6 / full_set)
    return orfs_per_million


def pct_gene_signal_per_orf(aORF):
# Value related to OPM above but as a proportional of expression from the locus. between 0 & 1
    total_aORF = 0
    for transcript in aORF:
        total_aORF += aORF[transcript]

    pct_ORF = {}
    for transcript in aORF:
        if transcript not in pct_ORF:
            pct_ORF[transcript] = aORF[transcript] / total_aORF

    return pct_ORF




def orfQuant(gene, sqlite_path_organism, sqlite_path_reads, supported, counts, exons, filter=True):
# Main function for implementing method. Executes above functions and determines if all features are shared 
# Returms OPM values from ORFs_per_million function
    junctions = genomic_junction_positions(gene, sqlite_path_organism, supported, exons)
    orf_coordinates = genomic_orf_coordinate_ranges(gene, sqlite_path_organism, supported, exons)
    junction_scores = genomic_junction_scores(gene, sqlite_path_organism, sqlite_path_reads, supported, exons, junctions,
                               filter=True)

    unique_shared_exons = classify_regions_shared_unique(exons)
    coverage_exons = region_coverage(exons, counts)

    unique_shared_junctions = classify_regions_shared_unique(junctions)
    coverage_junctions = coverage_junction_transcript(junction_scores)


    average_unique = average_unique_coverage(unique_shared_exons, unique_shared_junctions, coverage_junctions, coverage_exons)
    average_all = all_feature_average(coverage_junctions, coverage_exons)

    cORF = cORF_ratio(average_unique, average_all)

    all_shared = True
    for transcript in supported:
        if unique_shared_junctions[transcript]["unique"] == [] and unique_shared_exons[transcript]["unique"] == []:
            cORF[transcript] = adjusted_coverage_for_non_unique_orfs(transcript, coverage_exons, coverage_junctions, average_all)
        else:
            all_shared = False

    if all_shared:
        cORF = shared_coverage(coverage_exons, coverage_junctions)
    else:
        cORF = cORF_ratio(average_unique, average_all)

    adjusted_a_sites = aORF(cORF, counts)
    orf_lengths = lORF(supported, sqlite_path_organism)

    orfs_per_million = ORFs_per_million(adjusted_a_sites, orf_lengths)


    return orfs_per_million



def incl_OPM_run_orfQuant(gene, sqlite_path_organism, sqlite_path_reads, force=False):
    # Run the ORFquant algorithm on the files that do not have a OPM value stored for the required trancripts
    # Store calculated OPMs for a file back in the read sqlitedict under "OPM"
    # Returns the average OPMs of the selected files for this locus 
    gene_info = get_gene_info(gene, sqlite_path_organism)
    transcripts = [transcript[0] for transcript in gene_info]
    coding = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
    exons = genomic_exon_coordinate_ranges(gene, sqlite_path_organism, coding)

    transcript_exon_coordinate_ranges = tripsSplice.transcript_exon_coordinate_ranges(gene_info)
    transcript_OPMs = {}
    read_files = []
    for file in sqlite_path_reads:
        if force:read_files.append(file)
        print file
        infile = SqliteDict(file)
        for transcript in coding:
            if transcript in infile:
				if "OPM" in infile[transcript]:
					if transcript in transcript_OPMs:
						transcript_OPMs[transcript].append(infile[transcript]["OPM"])
					else:
						transcript_OPMs[transcript] = [infile[transcript]["OPM"]]
                    

				else:
					if file not in read_files:
						read_files.append(file)
        infile.close()


    for file in read_files:
        genomic_read_positions_1 = get_reads_per_genomic_location_asite(gene, [file], sqlite_path_organism,
                                                                  coding, exons, filter=True)

        genomic_read_positions_2 = get_reads_per_genomic_location_asite_faster(gene, [file], sqlite_path_organism,
                                                                  coding, exons, transcript_exon_coordinate_ranges, filter=True)

        set1 = set(genomic_read_positions_1.items())
        set2 = set(genomic_read_positions_2.items())
        print set1 ^ set2
        print 
        print 

        print set1 - set2 

        genomic_read_positions = genomic_read_positions_1 
        counts = count_read_supporting_regions_per_transcript(exons, genomic_read_positions)

        if genomic_read_positions != {}: # don't run calculation if no reads map to gene
            orfQuant_res = orfQuant(gene, sqlite_path_organism, [file], coding, counts, exons, filter=True)

        infile = SqliteDict(file, autocommit=False)

        for transcript in transcripts:
            transcript_dict = infile[transcript].copy()

            if transcript in orfQuant_res:
                transcript_dict['OPM'] = orfQuant_res[transcript]
                if transcript not in transcript_OPMs:
                    transcript_OPMs[transcript] = [orfQuant_res[transcript]]
                else:
                    transcript_OPMs[transcript].append(orfQuant_res[transcript])
            
            infile[transcript] = transcript_dict


        infile.commit()
        infile.close()



    avg_opms_per_transcript = {}
    for transcript in transcript_OPMs:
        avg_opms_per_transcript[transcript] = round(sum(transcript_OPMs[transcript])/float(len(transcript_OPMs[transcript])), 2)
    return avg_opms_per_transcript


# if __name__ == "__main__":
#     # gene = "phpt1"
#     gene = "SARS-COV2"

#     # sqlite_path_organism = "/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_annotations/homo_sapiens/old_homo_sapiens.Gencode_v25.sqlite"
#     sqlite_path_organism = '/home/jack/projects/trips/Trips-Viz/trips_annotations_sample/sarscov2/sarscov2.sarscov2_homo_sapiens.sqlite'

#     # sqlite_path_reads = ["/home/jackt/Tools/Trips-Viz/Trips-Viz-master/trips_shelves/rnaseq/homo_sapiens/Park16/SRR3306574.sqlite"]
#     sqlite_path_reads =  ['/home/jack/projects/trips/Trips-Viz/trips_shelves/riboseq/sarscov2/Finkel20/SRR12216749.sqlite']
#     orfQuant_res = incl_OPM_run_orfQuant(gene, sqlite_path_organism, sqlite_path_reads, force=True)

#     print orfQuant_res
