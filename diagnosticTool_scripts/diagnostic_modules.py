"""Kodoja pipeline."""
from __future__ import print_function

import subprocess
import pandas as pd
import random
import os
import pickle
import shutil
import json

from ete3 import NCBITaxa
from math import isnan

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# The user-facing scripts will all report this version number via --version:
version = "0.0.10"
ncbi = NCBITaxa()

def check_path(dirs):
    """Check if directory path has '/' at the end.

    Return value is either '/' or empty string ''.
    """
    if dirs[-1] != "/":
        return "/"
    else:
        return ""


def test_format(file1, user_format):
    """Check data format.

    Check if data is in the fasta or fastq format and
    assert the user has specified the correct format for
    the data provided.

    Return an assert stament and stop or continue.
    """
    with open(file1) as myfile:
        # Would have used xrange under Python 2, but want this to work
        # on both Python 2 and 3 and a list of 8 elements is tiny.
        small_file = [next(myfile) for x in range(8)]

    file_format = "not identified"

    if small_file[0][0] == "@" and small_file[4][0] == "@":
        file_format = "fastq"
    if small_file[0][0] == ">":
        file_format = "fasta"

    assert (file_format == "fasta") | (file_format == "fastq"), \
        "Cannot proceed with file as it is not in fasta or fastq format."
    assert user_format == file_format, \
        "File has been detected to be in " + file_format + \
        " format rather than " + user_format + " format."


def rename_seqIDs(input_file, out_dir, user_format, paired=False):
    """Rename sequence identifiers to just the read number.

    Write a new file where each sequence ID is replaced with
    the read number (counting from one).

    Does not attempt to include "/1" and "/2" name suffixes, nor
    include "1:" or "2:" in the description, for paired reads.

    Returns dictionary mapping the sequence number to the old
    identifier (first word only from the description line,
    and if paired without any "/1" or "/2" suffix).
    """
    if paired == 2:
        output_file = os.path.join(out_dir, "renamed_file_2." + user_format)
    elif paired == 1 or paired is False:
        output_file = os.path.join(out_dir, "renamed_file_1." + user_format)
    else:
        raise ValueError("Wanted 1, 2 or False - not %r" % paired)
    id_dict = {}
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        if user_format == 'fasta':
            for index, (title, seq) in enumerate(SimpleFastaParser(in_file)):
                name = title.split(None, 1)[0]
                if (paired == 1 and name.endswith("/1")) or (paired == 2 and name.endswith("/2")):
                    name = name[:-2]
                id_dict[index + 1] = name
                out_file.write(">%i\n%s\n" % (index + 1, seq))
        else:
            for index, (title, seq, qual) in enumerate(FastqGeneralIterator(in_file)):
                name = title.split(None, 1)[0]
                if (paired == 1 and name.endswith("/1")) or (paired == 2 and name.endswith("/2")):
                    name = name[:-2]
                id_dict[index + 1] = name
                out_file.write("@%i\n%s\n+\n%s\n" % (index + 1, seq, qual))
    return id_dict


def check_file(file1, out_dir, user_format, file2=False):
    """Rename sequnce ids and check PE files.

    Rename sequnce ids for SE or PE files to ensure
    consistency between kraken and kaiju (which modify
    id names). Create dictionaries containing real IDs and
    renamed version and pickle. If data is PE, assert
    paired files have the same number of entries and if
    the paired reads are matched by choosing random
    entries and confirming the IDs match (optionally
    with /1 and /2 suffixes).
    """
    if file2:
        ids1 = rename_seqIDs(file1, out_dir, user_format, paired=1)
        ids2 = rename_seqIDs(file2, out_dir, user_format, paired=2)
        with open(os.path.join(out_dir, 'ids2.pkl'), 'wb') as pkl_dict:
            pickle.dump(ids2, pkl_dict, protocol=pickle.HIGHEST_PROTOCOL)

        assert len(ids1) == len(ids2), \
            "Paired files have different number of reads"

        for values in range(1, 50):
            random_id = random.randint(1, len(ids1) - 1)
            id_1 = ids1[random_id]
            id_2 = ids2[random_id]
            assert id_1 == id_2, \
                ("Paired-end sequences don't match, e.g. %r vs %r"
                 % (id_1, id_2))
    else:
        ids1 = rename_seqIDs(file1, out_dir, user_format, paired=False)

    with open(os.path.join(out_dir, "log_file.txt"), "a") as log_file:
        log_file.write("Number of sequences = " + str(list(ids1)[-1]) + "\n")

    with open(os.path.join(out_dir, 'ids1.pkl'), 'wb') as pkl_dict:
        pickle.dump(ids1, pkl_dict, protocol=pickle.HIGHEST_PROTOCOL)


def fastqc_trim(out_dir, file1, trim_minlen, threads, adapter_file, file2=False):
    """Quality and adaptor trimming of fastq files.

    Takes fastq data (either single or paired), trims sequences using trimmomatic
    (in the case of paried end reads, it deletes extra files) and uses fastqc to
    show the user what the sequence quality looks like after trimming.

    Returns trimmed sequence files and fastq analysis files
    """
    trimAdapt_command = " LEADING:20 TRAILING:20 MINLEN:" + \
                        str(trim_minlen)
    if adapter_file:
        trimAdapt_command += " ILLUMINACLIP:" + adapter_file + ":2:30:10"

    if file2:
        PE_trim_command = "trimmomatic PE -threads " + str(threads) + " " + file1 + " " + file2 + \
            " " + os.path.join(out_dir, "trimmed_read1") + \
            " " + os.path.join(out_dir, "PE_trimmed_data_1U") + \
            " " + os.path.join(out_dir, "trimmed_read2") + \
            " " + os.path.join(out_dir, "PE_trimmed_data_2U") + trimAdapt_command

        subprocess.check_call(PE_trim_command, shell=True)
        os.remove(os.path.join(out_dir, "PE_trimmed_data_1U"))
        os.remove(os.path.join(out_dir, "PE_trimmed_data_2U"))
        subprocess.check_call("fastqc " + os.path.join(out_dir, "trimmed_read1") +
                              " -o " + out_dir, shell=True)
        subprocess.check_call("fastqc " + os.path.join(out_dir, "trimmed_read2") +
                              " -o " + out_dir, shell=True)
    else:
        subprocess.check_call("trimmomatic SE -threads " + str(threads) + " " + file1 +
                              " " + os.path.join(out_dir, "trimmed_read1") +
                              " " + trimAdapt_command, shell=True)
        subprocess.check_call("fastqc " + os.path.join(out_dir, "trimmed_read1") +
                              " -o " + out_dir, shell=True)


def kraken_classify(out_dir, kraken_file1, threads, kraken_db, memory, kraken_file2=False,
                    quick_minhits=False):
    """Kraken classification.

    Add appropiate switches for kraken command (format, minimum hits, memory
    if paired or single end) and call
    kraken command, followed by kraken-translate to get full taxonomy for each
    sequence based on thir sequence id (Seq_tax: d__superkingdom, k__kingdom,
    p__phylum, c__class, o__order, f__family, g__genus, s__species).

    Return kraken_table file with a row for each sequence and kraken classification
    (or unclassified) and kraken_labels file witha row for each sequence that was
    classified by kraken with full taxonomy.
    """

    kraken_command = "kraken2 "

    #kraken_command += "--threads " + str(threads) + " --report kraken.report --db " + kraken_db 
    kraken_command += "--threads " + str(threads) + " --report " + out_dir + "kraken.report --db " + kraken_db 
    
    if memory:
        kraken_command += " --memory-mapping "
        
    if quick_minhits:
        kraken_command += " --quick --min-hits " + str(quick_minhits)

    if kraken_file2:
        kraken_command += " --paired " + kraken_file1 + " " + \
                          kraken_file2 + " > " + os.path.join(out_dir, "kraken_table.txt")
    else:
        kraken_command += " " + kraken_file1 + " > " + os.path.join(out_dir, "kraken_table.txt")

    subprocess.check_call(kraken_command, shell=True)
    #subprocess.check_call("kraken-translate --mpa-format --db " + kraken_db +
    #                      " " + os.path.join(out_dir, "kraken_table.txt") + " > " +
    #                      os.path.join(out_dir, "kraken_labels.txt"), shell=True)

def filter_sequence_file(input_file, output_file, user_format, wanted,
                         ignore_suffix=None):
    """Create a subset of sequences based on sequence IDs.

    Writes a FASTA or FASTQ file in the output file specified.

    Argument wanted should be a Python set of identifers (with no white space,
    i.e. the first word only from the FASTA or FASTQ title lines).

    Optional argument ignore_suffix="/1" means remove any suffix "/1"
    from the input read names before matching to the wanted list.
    The suffix is retained in the output file.
    """
    print("Selecting %i unique identifiers" % len(wanted))
    if ignore_suffix:
        cut = len(ignore_suffix)
        records = (r for r in SeqIO.parse(input_file, user_format)
                   if (r.id[:-cut] if r.id.endswith(ignore_suffix) else r.id) in wanted)
    else:
        records = (r for r in SeqIO.parse(input_file, user_format)
                   if r.id in wanted)
    count = SeqIO.write(records, output_file, user_format)
    print("Saved %i records from %s to %s" % (count, input_file, output_file))
    if count < len(wanted):
        print("Warning %i IDs not found in %s" % (len(wanted) - count, input_file))

def kaiju_classify(kaiju_file1, threads, out_dir, kaiju_db, kaiju_minlen, kraken_db,
                   kaiju_file2=False, kaiju_mismatch=False, kaiju_score=False):
    """Run kaiju command for kaiju classification of sequences.

    It ensures that if mismatches are allowed that a score has also been provided.
    Once classification is complete, it uses kraken-translate (as in kraken_classify())
    to get full taxonomy names for each sequence that has been classified. It deletes the files
    used for this analysis.

    """
    #kaiju_nodes = kraken_db + "taxonomy/nodes.dmp"
    kaiju_nodes = kaiju_db + "../nodes.dmp"
    kaiju_fmi = kaiju_db + "kaiju_library.fmi"
    # kaiju_names = kaiju_db + "names.dmp"

    if kaiju_mismatch:
        assert(kaiju_score), "Set kaiju_score for greedy mode"
        mode = "greedy -e " + str(kaiju_mismatch) + " -s " + str(kaiju_score)
    else:
        mode = "mem"

    kaiju_command = "kaiju -z " + str(threads) + " -t " + kaiju_nodes + " -f " + kaiju_fmi + \
                    " -i " + kaiju_file1 + " -o " + os.path.join(out_dir, "kaiju_table.txt") + \
                    " -x -v -a " + mode + " -m " + str(kaiju_minlen)

    if kaiju_file2:
        kaiju_command += " -j " + kaiju_file2

    subprocess.check_call(kaiju_command, shell=True)
    #subprocess.check_call("kraken-translate --mpa-format --db " + kraken_db + " " +
    #                      os.path.join(out_dir, "kaiju_table.txt") + " > " +
    #                      os.path.join(out_dir, "kaiju_labels.txt"), shell=True)

    for dirs, sub_dirs, files in os.walk(out_dir):
        # Only delete file when it's in out_put
        for filenames in files:
            if kaiju_file1 == filenames:
                os.remove(kaiju_file1)
                if kaiju_file2:
                    os.remove(kaiju_file2)

def add_krona_representation(out_dir):
    """Make krona representation (html file) with kaiju_table.txt and kraken_table.txt
    """
    kraken_file = os.path.join(out_dir, "kraken_table.txt")
    kaiju_file = os.path.join(out_dir, "kaiju_table.txt")
    kraken_html_file = os.path.join(out_dir, "kraken.html")
    kaiju_html_file = os.path.join(out_dir, "kaiju.html")
    kraken_krona_command = "ktImportTaxonomy -o " + kraken_html_file + " -t 3 " + kraken_file
    kaiju_krona_command = "ktImportTaxonomy -o " + kaiju_html_file + " -t 3 " + kaiju_file

    subprocess.check_call(kraken_krona_command, shell=True)
    subprocess.check_call(kaiju_krona_command, shell=True)
    try:        
        shutil.rmtree(os.path.join(out_dir, "kraken.html.files"))
        shutil.rmtree(os.path.join(out_dir, "kaiju.html.files"))
    except FileNotFoundError:
        pass

def get_desired_ranks(taxid, desired_ranks):
    """
    give taxonomic name and ID for all rank in desired_ranks
    """
    # get lineage (current taxid)
    lineage = ncbi.get_lineage(taxid)   
    # get names
    names = ncbi.get_taxid_translator(lineage)
    # get current taxid : rank
    lineage2ranks = ncbi.get_rank(names)
    # {rank: {current taxid: name}}
    ranks2lineage = dict((rank, ncbi.get_taxid_translator([taxid])) for (taxid, rank) in lineage2ranks.items())
    # reformat like {taxid_initial: {rank: {current taxid: name}}}
    return {format(infos) : ranks2lineage.get(infos, 0) for infos in desired_ranks} 

def clean_dataset(kodoja_seq_data):
    """
    Get some number to calcul statistic later
    and remove row/read if both tools don't have result for that row/reads
    """ 
    # TODO store read name when no result for both tool to adress dark matter PBM
    #   
    # save the number of reads that have been processed (one pair is count once)
    total_reads_nb = len(kodoja_seq_data.index)
    # when kaiju don't have result, put 0 instead of Nan (NoneType)
    kodoja_seq_data["kaiju_Tax_ID"] = kodoja_seq_data["kaiju_Tax_ID"].fillna(0)
    # remove row if both kraken and kaiju don't have result
    kodoja_seq_data = kodoja_seq_data.drop(kodoja_seq_data[
        (kodoja_seq_data["kraken_Tax_ID"] == 0) & (kodoja_seq_data["kaiju_Tax_ID"] == 0)].index)
    # read nb that have been assigned by at least one tool
    total_known_read = len(kodoja_seq_data.index)
    # read nb that NOT have been assigned by at least one tool
    total_unknown_reads = total_reads_nb - total_known_read

    return total_reads_nb, total_known_read, total_unknown_reads

def get_combination(kodoja_seq_data):
    """
    put row/read to combined result is same result for both tool
    put in either if at least one of the tool classified the row/read (meanning that all reads in conbine are in either)
    """
    # check if reads have same assignation
    # TODO check all rank od interest not only the most precise one
    kodoja_seq_data["combined_result"] = kodoja_seq_data.kraken_Tax_ID[kodoja_seq_data["kraken_Tax_ID"] ==
        kodoja_seq_data["kaiju_Tax_ID"]]
    
    # when combined_result don"t have result, put 0 instead of Nan (NoneType)
    kodoja_seq_data["combined_result"] = kodoja_seq_data["combined_result"].fillna(0)
    # convert to int (default float)
    kodoja_seq_data[["combined_result"]] = kodoja_seq_data[["combined_result"]].astype(int)
    # Number of same sequences classified to same taxID by both tools
    combined_class = dict(kodoja_seq_data["combined_result"].value_counts())
    combined_class.pop(0)

    # get kraken nb of read by taxid
    kraken_class = dict(kodoja_seq_data["kraken_Tax_ID"].value_counts())
    # get kaiju nb of read by taxid
    kaiju_class = dict(kodoja_seq_data["kaiju_Tax_ID"].value_counts())
    # Number of sequences classified to taxID by at least one of the tools
    either_class = {k: kraken_class.get(k, 0) + kaiju_class.get(k, 0) for k in set(kraken_class) | set(kaiju_class)}
    either_class.pop(0, None)

    return either_class, combined_class, kraken_class, kaiju_class

def get_comparison_list(seq_orgData_list, desired_ranks2, tax_name_id_path_all, combined_class, 
    taxid, nbreads, kraken_class, kaiju_class, total_reads_nb, total_known_read):
    """
    writing in a list everything need to make an easy human readable result table (virus_table.tsv)
    """
    def add_NBreads(dictresult, taxid, results):
        """
        count nb of read in a taxid
        """
        switch=True
        for key, NBreads in dictresult.items():
            if key==taxid:
                results.append(str(NBreads))
                switch=False
                break
        if switch:
            results.append(str(0))
        return results
    
    seq_orgData_list.append(taxid)
    for rank in desired_ranks2:
        if tax_name_id_path_all[taxid][rank] != 0:
            seq_orgData_list.append(list(tax_name_id_path_all[taxid][rank].values())[0])
        else:
            seq_orgData_list.append("None")

    nb_common = int(add_NBreads(combined_class, taxid, list())[0])
    # print(nb_common)
    nb_either = nbreads - nb_common
    # NBreads_either
    seq_orgData_list.append(nb_either)
    # NBreads_common
    seq_orgData_list.append(nb_common)
    # percentage_common
    percentage_common = str((nb_common/nb_either)*100) + "%"
    seq_orgData_list.append(percentage_common)
    # NBreads_kaiju
    nbreads_kaiju = int(add_NBreads(kaiju_class, taxid, list())[0])
    seq_orgData_list.append(nbreads_kaiju)
    # NBreads_kraken
    nbreads_kraken = int(add_NBreads(kraken_class, taxid, list())[0])
    seq_orgData_list.append(nbreads_kraken)
    # either_absolute
    either_absolute = str((nb_either/total_reads_nb)*100) + "%"
    seq_orgData_list.append(either_absolute) #keep that one
    # either_logical
    either_logical = str((nb_either/total_known_read)*100) + "%"
    seq_orgData_list.append(either_logical)
    # common_absolute
    common_absolute = str((nb_common/total_reads_nb)*100) + "%"
    seq_orgData_list.append(common_absolute)
    # common_logical
    common_logical = str((nb_common/total_known_read)*100) + "%"
    seq_orgData_list.append(common_logical)   

    return seq_orgData_list

def format_result_table(kodoja_seq_data, host_subset):
    """
    Merge the classification data (either kraken or kaiju) with the 'label'
    data which has full taxonomy for the classified sequence.
    """
    label_colNames = ["Seq_ID", "Seq_ID_initial", "Tax_ID_kraken", "Seq_tax_kraken", 
        "Rank_kraken", "Tax_ID_kaiju", "Seq_tax_kaiju", "Rank_kaiju"]
    org_colNames = ["Original_query_taxid", "superkingdom", "family", "genus", 
        "species", "NBreads_either", "NBreads_common", "percentage_common", 
        "NBreads_kaiju", "NBreads_kraken", "pg_either_absolute", "pg_either_logical",
         "pg_common_absolute", "pg_common_logical"]
    desired_ranks = ["species","genus","family","superkingdom"]
    desired_ranks2 = ["superkingdom", "family", "genus", "species"] 
    seq_orgData = pd.DataFrame(columns=org_colNames)
    seq_labelDataList = []

    # clean main dataset and get raw number to make stats
    total_reads_nb, total_known_read, total_unknown_reads = clean_dataset(kodoja_seq_data)
    # Get reads assignation combination between kreaken and kaiju
    either_class, combined_class, kraken_class, kaiju_class = get_combination(kodoja_seq_data)

    #TODO check all rank of interest not only the most precise one
    #TODO Adding number of read for each rank to tax_name_id_path_all ? 
    tax_name_id_path_all={}   #12182: {'superkingdom': {10239: 'Viruses'}, ... 'species': {12182: 'Potato aucuba mosaic virus'}}
    tax_name_id_path={}   #12182: { ['species', 'Potato aucuba mosaic virus']}
    for taxid, nbreads in either_class.items():
        tax_name_id_path_all.update({taxid:get_desired_ranks(taxid, desired_ranks)})
        taxid2name = ncbi.get_taxid_translator([str(taxid)])
        taxid2rank = ncbi.get_rank([taxid])
        tax_name_id_path.update({taxid:[taxid2rank[taxid], taxid2name[taxid]]})

        seq_orgData_list_tmp = get_comparison_list(
            list(), desired_ranks2, tax_name_id_path_all, 
            combined_class, taxid, nbreads, kraken_class, kaiju_class,
            total_reads_nb, total_known_read)
             
        seq_orgData = seq_orgData.append(pd.Series(seq_orgData_list_tmp, index = seq_orgData.columns ), ignore_index=True)   
    #add unknown/dark matter count total_unknown_reads
    either_absolute = str((total_unknown_reads/total_reads_nb)*100) + "%"
    either_logical = str((total_unknown_reads/total_known_read)*100) + "%"
    seq_orgData_list_tmp = [-1, "Unassigned", "Unassigned", "Unassigned", "Unassigned", total_unknown_reads, 0, "0%", 0, 
    0, either_absolute, either_logical, "0%", "0%"]
    seq_orgData = seq_orgData.append(pd.Series(seq_orgData_list_tmp, index = seq_orgData.columns ), ignore_index=True)


    for element in kodoja_seq_data.itertuples():
        ### find taxonomic name and rank if possible
        # for kraken
        if element.kraken_Tax_ID !=0:
            # use the stored taxonomic name find the one belonging to each read 
            seq_tax_kraken = tax_name_id_path[element.kraken_Tax_ID][1]
            rank_kraken = tax_name_id_path[element.kraken_Tax_ID][0]
        else:
            seq_tax_kraken = 0
            rank_kraken = "None"
        # for kaiju
        if element.kaiju_Tax_ID !=0:
            # use the stored taxonomic name find the one belonging to each read 
            seq_tax_kaiju = tax_name_id_path[element.kaiju_Tax_ID][1]
            rank_kaiju = tax_name_id_path[element.kaiju_Tax_ID][0]
        else:
            seq_tax_kaiju = 0
            rank_kaiju = "None"     
        ### write taxonomic name and rank

        #####performance issue: 
        # https://stackoverflow.com/questions/36489576/why-does-concatenation-of-dataframes-get-exponentially-slower/36489724#36489724
        # seq_labelData = seq_labelData.append(pd.concat([pd.DataFrame([[ 
        #     element.Seq_ID, # Seq_ID
        #     element.Seq_ID_initial, # Seq_ID_initial
        #     element.kraken_Tax_ID, # Tax_ID_kraken
        #     seq_tax_kraken, # Seq_tax_kraken
        #     rank_kraken, # Rank_kraken
        #     element.kaiju_Tax_ID, # Tax_ID_kaiju
        #     seq_tax_kaiju, # Seq_tax_kaiju
        #     rank_kaiju # Rank_kaiju
        #     ]], columns=label_colNames)])) 

        seq_labelDataList.append((      
            element.Seq_ID, # Seq_ID
            element.Seq_ID_initial, # Seq_ID_initial
            element.kraken_Tax_ID, # Tax_ID_kraken
            seq_tax_kraken, # Seq_tax_kraken
            rank_kraken, # Rank_kraken
            element.kaiju_Tax_ID, # Tax_ID_kaiju
            seq_tax_kaiju, # Seq_tax_kaiju
            rank_kaiju # Rank_kaiju)
        ))
    seq_labelData = pd.DataFrame(seq_labelDataList, columns = label_colNames)            
    
    if host_subset:
        seq_labelData = seq_labelData[(seq_labelData["Tax_ID_kraken"] != float(host_subset)) &
                        (seq_labelData["Tax_ID_kaiju"] != float(host_subset))]

    if hasattr(seq_labelData, "sort_values"):
        # pandas 0.17 onwards
        seq_labelData.sort_values(["Seq_ID"], inplace=True)
    seq_labelData.reset_index(drop=True, inplace=True)

    return seq_labelData, seq_orgData, tax_name_id_path_all

def result_comparaison(out_dir, kraken_table,  kaiju_table, host_subset):
    """
    Open kraken/kaiju result file 
    launch comparative analaysis of the result
    clean result repository
    write output file 
    """

    kraken_colNames = ["kraken_classified", "Seq_ID", "kraken_Tax_ID", "kraken_length", "kraken_k-mer"]
    kaiju_colNames = ["kaiju_classified", "Seq_ID", "kaiju_Tax_ID", "kaiju_lenBest", "kaiju_tax_AN", 
        "kaiju_accession", "kaiju_fragment"]
    
    # open kraken result file
    kraken_seq_data = pd.read_csv(os.path.join(out_dir, kraken_table),
                            sep="\t", header=None, names=kraken_colNames,
                            index_col=False)
    # open kaiju result file
    kaiju_seq_data = pd.read_csv(os.path.join(out_dir, kaiju_table),
                            sep="\t", header=None, names=kaiju_colNames,
                            index_col=False)
    # check if same number of row a.k.a. number of reads beetween kraken kaiju
    # TODO if not, exit(1) with warning
    if len(kraken_seq_data.index) !=  len(kaiju_seq_data.index):
        print("WARNING: kraken and kaiju don't have the same number of reads processed")
    # all result on one table                        
    kodoja_seq_data = pd.merge(kraken_seq_data, kaiju_seq_data, on="Seq_ID", how="outer")
    
    # open pickle file
    with open(os.path.join(out_dir, "ids1.pkl"), "rb") as id_dict:
        ids1 = pickle.load(id_dict)
    # add the reads name in a new column
    kodoja_seq_data["Seq_ID_initial"] = kodoja_seq_data["Seq_ID"].map(ids1)

    kodoja_format_seq_data, virus_table, tax_name_id_path_all= format_result_table(kodoja_seq_data, host_subset)
   
    #TODO faire ca
    os.remove(os.path.join(out_dir, "ids1.pkl"))
    os.remove(os.path.join(out_dir, "ids2.pkl"))

    kodoja_format_seq_data.to_csv(os.path.join(out_dir, "kodoja_VRL.tsv"),
                   sep="\t", index=False)
    virus_table.to_csv(os.path.join(out_dir, "virus_table.tsv"),
                   sep="\t", index=False)

    # NOTE store this for further analysis
    #12182: {'superkingdom': {10239: 'Viruses'}, ... 'species': {12182: 'Potato aucuba mosaic virus'}}
    with open(os.path.join(out_dir, "Additional_taxonomic_prediction.json"), 'w') as file:
        file.write(json.dumps(tax_name_id_path_all)) # use `json.loads` to do the reverse     
    kodoja_seq_data.to_csv(os.path.join(out_dir, "Additional_all_data.tsv"),
                sep="\t", index=False)     

    #TODO write genus_taxid.pkl for kodoja_retrieve 
    # diagnostic_moduleOLD l443