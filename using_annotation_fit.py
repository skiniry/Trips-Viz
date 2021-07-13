from annotation_fit import * 
from fetch_shelve_reads2 import get_reads
from core_functions import fetch_studies, fetch_files,fetch_study_info,fetch_file_paths,generate_short_code,fetch_user
import config
import pandas as pd 
import matplotlib as mpl



def get_genes_for_organism(connection, coding=False):
    '''
    return list of gene names for organism sqlite 
    '''
    curs = connection.cursor()
    if coding:
        result = curs.execute("SELECT DISTINCT gene FROM transcripts WHERE gene_type == 1").fetchall()
    else:
        result = curs.execute("SELECT DISTINCT gene FROM transcripts").fetchall()

    genes = [i[0] for i in result]

    return genes

def get_file_paths(username, organism, transcriptome):
    accepted_studies = fetch_studies(username, organism, transcriptome)
    file_id_to_name_dict,accepted_studies,accepted_files,seq_types = fetch_files(accepted_studies)
    print accepted_files['riboseq'].keys()
    file_list = []
    for study in accepted_studies.keys():
        for file in accepted_files['riboseq'][study]:
            file_list.append(file)
    file_path_dict = fetch_file_paths(file_list,organism)
    return file_path_dict


def get_base_coverage(reads, cutoff_list, tranlen):
    '''
    given a list of cutoffs return the coverage scores from the read file that hit that cutoff
    '''
    coverage_scores = [] 
    for coverage_limit in cutoff_list:
        covered_positions = [position for position in reads[0] if reads[0][position] >= coverage_limit]
        proportion_covered = round(float(len(covered_positions))/float(tranlen), 2)
        coverage_scores.append((proportion_covered, coverage_limit))
    
    return coverage_scores


ambig = "unambig"
minread = 5
maxread = 150
ribocoverage = False
noisered = False
primetype = "fiveprime"
readscore = 0
secondary_readscore = 1
pcr = False

organism = "sarscov2"
transcriptome = "sarscov2_homo_sapiens"
# organism = "homo_sapiens"
# transcriptome = "Gencode_v25"


# ecoli_sqlite = '/home/jack/projects/trips/Trips-Viz/trips_annotations_sample/escherichia_coli/escherichia_coli.Ensembl_k_12_ASM584v2.sqlite'
sqlite = '/home/jack/projects/trips/Trips-Viz/trips_annotations_sample/sarscov2/sarscov2.sarscov2_homo_sapiens.sqlite'

def region_scores_calculation_to_csv(sqlite, output, ambig, minread, maxread, ribocoverage, organism, transcriptome, noisered, primetype, readscore, secondary_readscore, pcr):
    data = {'five_prime_counts': [], 'merged_length': [], 'coding_counts': [], 'three_prime_length': [], 'three_prime_counts': [], 'coding_length': [], 'five_prime_length': [], 'merged_counts': [], 'total_reads': []}
    df = pd.DataFrame(data)
    file_paths_dict = get_file_paths("test", organism, transcriptome)

    connection = get_annotation_connection(sqlite)
    gene_list = get_genes_for_organism(connection, coding=True)[-10:]
    total_genes = float(len(gene_list))

    for gene_number, gene in enumerate(gene_list):
        start = time.time()
        print gene , '\t',float(gene_number)/total_genes * 100, '%'

        gene_info = get_gene_info(gene, connection)
        transcript_architecture_dict = transcript_architecture(gene_info)
        locus_counts = {}
        if transcript_architecture_dict == {}:
            continue
        for transcript in gene_info:
            transcript_id = transcript[0]
            cds = transcript_architecture_dict[transcript_id]['coding']
            time1 = time.time()
            coding_regions = get_merged_cds_coordinates(transcript_id, cds, connection)
            coding_regions.sort()
            time2 = time.time()

            tranlen = len(transcript[5])

            reads = get_reads(ambig, minread, maxread, transcript[0], file_paths_dict,tranlen,ribocoverage, organism, True,noisered, primetype,"riboseq",readscore,secondary_readscore,pcr)
            annotation_fit_res = annotation_fit(reads, transcript[0], connection)
            annotation_fit_res["gene"] = gene 
            annotation_fit_res["transcript"] = transcript[0]
            try:
                annotation_fit_res["annotation_fit"] = round(float(annotation_fit_res['coding_counts'])/float(annotation_fit_res['total_reads']), 2) * 100
                annotation_fit_res["merged_fit"] = round(float(annotation_fit_res['merged_counts'])/float(annotation_fit_res['total_reads']), 2) * 100
            except ZeroDivisionError: 
                annotation_fit_res["annotation_fit"] = 0
                annotation_fit_res["merged_fit"] = 0   

            try: 
                annotation_fit_res['in_frame_score'] = float(annotation_fit_res["in_frame_counts"])/float(annotation_fit_res["coding_counts"]) 
            except ZeroDivisionError: 
                annotation_fit_res['in_frame_score'] = 0

            coverage_scores = get_base_coverage(reads, [1, 10, 100], tranlen)
            for proportion_limit_pair in coverage_scores:
                var_string = "covered_bases_limit_" + str(proportion_limit_pair[1])
                annotation_fit_res[var_string] = proportion_limit_pair[0]
            df = pd.DataFrame.append(df, annotation_fit_res, ignore_index=True)
                # print gene, " done: ", time.time() - start 
                # print 
            # except: 
            #     continue

    df.to_csv(output)

def csv_to_stats_df(csv_path):
    '''
    Take an input csv with counts per region and calculate statistics 
    Annotation fit = cds reads / all reads as % 
    merged_fit = merged reads / all reads as %
    '''
    df = pd.read_csv(csv_path) 
    genes = list(df["gene"].unique())
    annotation_fit_results = pd.DataFrame()
    genes = genes[-10:]
    ouptut_dict = {}
    # df.coding_counts.hist()
    # mpl.pyplot.show()
    for gene in genes:
        gene_df = df[df["gene"] == gene]
        gene_scores = {}
        for index, row in gene_df.iterrows():
            if row["transcript"] not in gene_scores:
                gene_scores[row['transcript']] = {"annotation_fit":0, "merged_fit":0}
                try:
                    gene_scores[row['transcript']]["annotation_fit"] = round(row['coding_counts']/row['total_reads'], 2) * 100
                except:
                    continue
                try:
                    gene_scores[row['transcript']]["merged_fit"] = round(row['merged_counts']/row['total_reads'], 2) * 100
                except:
                    continue
            if row['transcript'] not in ouptut_dict:
                ouptut_dict[row['transcript']] = dict(row)
                ouptut_dict[row['transcript']]['annotation_fit'] = gene_scores[row['transcript']]['annotation_fit']
                ouptut_dict[row['transcript']]['merged_fit'] = gene_scores[row['transcript']]['merged_fit'] 
    for i in ouptut_dict:
        print i, ouptut_dict[i]
            
            # print row["gene"], row["transcript"],row["coding_counts"]
        # print gene_df
    # for col in df:
    #     print df[col] 

# csv_to_stats_df('sarscov2_homo_sapiensout.csv')
region_scores_calculation_to_csv(sqlite, "test_out.csv", ambig, minread, maxread, ribocoverage, organism, transcriptome, noisered, primetype, readscore, secondary_readscore, pcr)