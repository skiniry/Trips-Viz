from tripsSplice import get_reads_per_transcript_location
from tripsSplice import get_protein_coding_transcript_ids
from tripsSplice import get_start_stop_codon_positions
from tripsSplice import get_gene_info
from tripsSplice import get_transcript_length

def get_counts_meanFLD(transcripts, read_file):
    length_freq = {}
    transcript_counts = {}
    for transcript in transcripts:
        if transcript not in transcript_counts:
            transcript_counts[transcript] = 0
        #print "transcript, read_file", transcript, read_file
        reads = get_reads_per_transcript_location(transcript, read_file)
        
        for length in reads:
            for position in reads[length]:
                if length not in length_freq:
                    length_freq[length] = reads[length][position]

                else:
                    length_freq[length] += reads[length][position]
                transcript_counts[transcript] += reads[length][position]
    sum = 0
    count = 0
    for length in length_freq:
        sum += length * length_freq[length]
        count += length_freq[length]
    mean_fld = round(float(sum) /float(count), 0) 
    return transcript_counts, mean_fld


def transcript_reads_per_kilobase(transcript_counts, cds_lengths, meanFLD):
    RPK = {}
    for transcript in cds_lengths:
        effective_length = cds_lengths[transcript] - meanFLD + 1
        effective_length_per_kilobase = effective_length/1000
        if transcript not in RPK:
            RPK[transcript] = transcript_counts[transcript]/effective_length_per_kilobase

    return RPK


def TPM(gene, sqlite_path_organism, sqlite_path_reads, type):
    # Calculate transcripts per million for a given read type 
    if type == "ribo":
        transcripts = get_protein_coding_transcript_ids(gene, sqlite_path_organism)
        start_stops = {transcript:get_start_stop_codon_positions(transcript, sqlite_path_organism) for transcript in transcripts}
        lengths = {transcript:start_stops[transcript][1] - start_stops[transcript][0] for transcript in start_stops}

    elif type == "rna":
        transcripts = [transcript[0] for transcript in get_gene_info(gene, sqlite_path_organism)]
        lengths = get_transcript_length(gene, sqlite_path_organism)


    all_TPMs = {transcript:[] for transcript in transcripts}

    for read_file in sqlite_path_reads:
        counts, meanFLD = get_counts_meanFLD(transcripts, read_file)


        RPK = transcript_reads_per_kilobase(counts, lengths, meanFLD)
        per_million_scaling_factor = sum(RPK.values())/1000000

        TPM = {transcript:round(RPK[transcript]/per_million_scaling_factor, 2) for transcript in transcripts}
        for transcript in TPM:
            all_TPMs[transcript].append(TPM[transcript])

    avg_TPM = {transcript:sum(all_TPMs[transcript])/len(all_TPMs[transcript]) for transcript in all_TPMs}

    return avg_TPM

