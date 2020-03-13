"""Kodoja pipeline."""
from __future__ import print_function

import subprocess
import pandas as pd
import random
import os
import pickle
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
    kraken_command += "--threads " + str(threads) + " --report " + out_dir + "/kraken.report --db " + kraken_db 
    
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


def format_result_table(out_dir, data_table, table_colNames):
    """Merge classification and label data.

    Merge the classification data (either kraken or kaiju) with the 'label'
    data which has full taxonomy for the classified sequence.

    Return merged table
    """
    #label_colNames = ["Seq_ID", "Seq_tax"]
    label_colNames = ["Seq_ID", "Seq_tax", "Rank"]
    seq_data = pd.read_csv(os.path.join(out_dir, data_table),
                           sep="\t", header=None, names=table_colNames,
                           index_col=False)
    seq_data_clean = seq_data[[table_colNames[0], table_colNames[1], table_colNames[2]]].copy()
    
    seq_labelData = pd.concat([pd.DataFrame([[ el.Seq_ID, ncbi.get_taxid_translator([el.Tax_ID])[el.Tax_ID], 
                            ncbi.get_rank([el.Tax_ID])[el.Tax_ID]]], columns=label_colNames) 
                            for el in seq_data_clean.itertuples() if el.Tax_ID!=0])
    # give  a proper exit if no kraken or kaiju result ?
    # ValueError: No objects to concatenate
    seq_result = pd.merge(seq_data_clean, seq_labelData, on='Seq_ID', how='outer')
  
    return seq_result


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


def seq_reanalysis(kraken_table, kraken_labels, out_dir, user_format, forSubset_file1,
                   forSubset_file2=False):
    """Format table and subset sequences for kaiju analysis.

    Add label to kraken_table  using format_result_table() and write to disk
    (delete kraken_table).
    If subset = True, make a list of "Seq_ID" column value if sequence is unclassified
    in "Classified" column or classified as VRL (virus) in column "Div_ID". This list will be
    used to subset sequences using sequence_subset(). This should be used when the host plant
    genome is used to classify sequences.

    Return merged kraken tableresult tables and subsetted sequence files (i subset=True).
    """
    kraken_colNames = ["kraken_classified", "Seq_ID", "Tax_ID", "kraken_length",
                       "kraken_k-mer"]
    kraken_fullTable = format_result_table(out_dir, "kraken_table.txt", kraken_colNames)
    kraken_results = kraken_fullTable[["kraken_classified", "Seq_ID", "Tax_ID", "Seq_tax", "Rank"]]
    kraken_results.to_csv(os.path.join(out_dir, 'kraken_VRL.txt'),
                          sep='\t', index=False)

    with open(os.path.join(out_dir, 'ids1.pkl'), 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    kraken_fullTable["Seq_ID"] = kraken_fullTable["Seq_ID"].map(ids1)
    kraken_fullTable.to_csv(os.path.join(out_dir, "kraken_FormattedTable.txt"),
                            sep='\t', index=False)
    if os.path.isfile(os.path.join(out_dir, "kraken_FormattedTable.txt.gz")):
        os.remove(os.path.join(out_dir, "kraken_FormattedTable.txt.gz"))
    subprocess.check_call("gzip " + os.path.join(out_dir, "kraken_FormattedTable.txt"),
                          shell=True)
    
    #os.remove(os.path.join(out_dir, "kraken_table.txt"))



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
        os.rmdir(os.path.join(out_dir, "kraken.html.files"))
        os.rmdir(os.path.join(out_dir, "kaiju.html.files"))
    except FileNotFoundError:
        pass

def result_analysis(out_dir, kraken_VRL, host_subset):
    """Kodoja results table.

    Imports kraken results table, formats kaiju_table and merges
    kraken and kaiju results into one table (kodoja). It then makes a table with
    all identified species and count number of intances for each using virusSummary().
    """
    kraken_results = pd.read_csv(os.path.join(out_dir + kraken_VRL),
                                 header=0, sep='\t',
                                 dtype={"kraken_classified": str, "Seq_ID": int,
                                        "Tax_ID": int, "Seq_tax": str})

    kaiju_colNames = ["kaiju_classified", "Seq_ID", "Tax_ID", "kaiju_lenBest",
                      "kaiju_tax_AN", "kaiju_accession", "kaiju_fragment"]
    kaiju_fullTable = format_result_table(out_dir, "kaiju_table.txt", kaiju_colNames)
    # kaiju_fullTable['Seq_ID'] = kaiju_fullTable['Seq_ID'].astype(float)
    # kaiju_fullTable['Seq_ID'] = kaiju_fullTable['Seq_ID'].astype(int)
    kaiju_results = kaiju_fullTable[["kaiju_classified", "Seq_ID", "Tax_ID", "Seq_tax", "Rank"]]
    with open(os.path.join(out_dir, 'ids1.pkl'), 'rb') as id_dict:
        ids1 = pickle.load(id_dict)
    kaiju_fullTable["Seq_ID"] = kaiju_fullTable["Seq_ID"].map(ids1)
    kaiju_fullTable.to_csv(os.path.join(out_dir, 'kaiju_FormattedTable.txt'),
                           sep='\t', index=False)
    if os.path.isfile(os.path.join(out_dir, "kaiju_FormattedTable.txt.gz")):
        os.remove(os.path.join(out_dir, "kaiju_FormattedTable.txt.gz"))
    subprocess.check_call('gzip ' + os.path.join(out_dir, 'kaiju_FormattedTable.txt'),
                          shell=True)

    kodoja = pd.merge(kraken_results, kaiju_results, on='Seq_ID', how='outer')
    # REVIEW does it mean that there is always more kraken result than kaiju ?
    # to confirm with check default parameters for kaiju and kraken (kaiju seems more stringent than kraken) 
    assert len(kraken_results) == len(kodoja), \
        'ERROR: Kraken and Kaiju results not merged properly'
    if hasattr(kodoja, 'sort_values'):
        # pandas 0.17 onwards
        kodoja.sort_values(['Seq_ID'], inplace=True)
    else:
        kodoja.sort(['Seq_ID'], inplace=True)
    kodoja.reset_index(drop=True, inplace=True)
    kodoja.rename(columns={"Seq_tax_x": "kraken_seq_tax", "Seq_tax_y": "kaiju_seq_tax",
                           'Tax_ID_x': 'kraken_tax_ID', 'Tax_ID_y': 'kaiju_tax_ID', 'Rank_x': 'kraken_rank', 'Rank_y': 'kaiju_rank'}, inplace=True)

    kodoja["Seq_ID"] = kodoja["Seq_ID"].map(ids1)

    # TODO protein completness metrics with kaiju result
    # Not remove that file anymore, better to keep to be abbl to dig depper
    # example calculing protein completness metrics with kaiju result 
    #os.remove(os.path.join(out_dir, "kaiju_table.txt"))
    os.remove(os.path.join(out_dir, "kraken_VRL.txt"))
    # these pkl file are not useful anymore see  https://github.com/abaizan/kodoja/pull/28
    os.remove(os.path.join(out_dir, "ids1.pkl"))
    os.remove(os.path.join(out_dir, "ids2.pkl"))
    kodoja['combined_result'] = kodoja.kraken_tax_ID[kodoja['kraken_tax_ID'] == kodoja['kaiju_tax_ID']]
    if host_subset:
        kodoja = kodoja[(kodoja['kraken_tax_ID'] != float(host_subset)) &
                        (kodoja['kaiju_tax_ID'] != float(host_subset))]
    kodoja.to_csv(os.path.join(out_dir, 'kodoja_VRL.txt'),
                  sep='\t', index=False)

    def virusSummary(kodoja_data):
        """Merge tables to create summary table.

        Creates a summary table with virus species names, tax id, count of
        sequences by kraken, kaiju and sequences that were identified by both
        tools as belonging to that species.

        For each tax id, a sequence count for kraken, kaiju and the combined
        is made. '_levels' dict have all tax ids present in th table with the
        taxanomic 'labels' given by kraken-traslate.

        'associated_tax' dict, has tax ids which would be related to a species
        tax id, as they belong to taxa which are higher, and therefore if
        they could belong to a species but cannot be identified specifically
        (i.e. a sequence whih has been given the following label
        'd__Viruses|f__Closteroviridae|g__Ampelovirus' could be an unspecifically
        identified 'Grapevine_leafroll-associated_virus_4' the label for which is
        'd__Viruses|f__Closteroviridae|g__Ampelovirus|s__Grapevine_leafroll-associated_virus_4').
        """
        
        
        kraken_class = dict(kodoja_data['kraken_tax_ID'].value_counts())
        kraken_levels = pd.Series(kodoja_data.kraken_seq_tax.values,
                                  index=kodoja_data.kraken_tax_ID).to_dict()
        kaiju_class = dict(kodoja_data['kaiju_tax_ID'].value_counts())
        kaiju_levels = pd.Series(kodoja_data.kaiju_seq_tax.values,
                                 index=kodoja_data.kaiju_tax_ID).to_dict()

        # Number of sequences classified to taxID by both tools
        combined_class = dict(kodoja_data['combined_result'].value_counts())

        # Number of sequences classified to taxID by either tool
        either_class = kraken_class.copy()
        either_class.update(kaiju_class)
        either_class.pop(0, None)
        for key, value in either_class.items():
            if key in kraken_class:
                if key in kaiju_class:
                    either_class[key] = kraken_class[key] + kaiju_class[key]
                    if key in combined_class:
                        either_class[key] = either_class[key] - combined_class[key]
                else:
                    either_class[key] = kraken_class[key]
            else:
                either_class[key] = kaiju_class[key]

        levels_dict = kraken_levels.copy()
        levels_dict.update(kaiju_levels)
        levels_dict.pop(0, None)
        levels_dict = {k: levels_dict[k] for k in levels_dict if not isnan(k)}
        levels_tax = {key: list(map(str, value.split('|')))
                      for key, value in levels_dict.items()}

        LCA_tax = {}
        # NOTE NCBITaxa allow us to have the complete tree of tax_ID can we use that to improve this ?
        # example kaiju detect at genus level, kraken at species level and the two genus are concording

        # Iterate over a copy of the values as we may remove taxonomy entries
        for key, tax in list(levels_tax.items()):
            if tax[-1][0] != 's':
                LCA_tax[key] = tax[-1]
                levels_tax.pop(key)

        species_dict = {}
        for key in levels_tax:
            species_dict[key] = " ".join(levels_tax[key][-1][3:].split("_"))

        # Find the genus for each species
        genus_per_species = {}
        for key, value in levels_dict.items():
            if len(value.split('g__')) > 1:
                genus_per_species[key] = value.split('g__')[1].split('|')[0]
            else:
                genus_per_species[key] = ''

        # TaxID for genus
        genus_taxid = {}
        for key, value in LCA_tax.items():
            if value[0:3] == 'g__':
                genus = value[3:]
                if genus in genus_taxid:
                    genus_taxid[genus].append(key)
                else:
                    genus_taxid[genus] = [key]
        with open(os.path.join(out_dir, 'genus_taxid.pkl'), 'wb') as pkl_dict:
            pickle.dump(genus_taxid, pkl_dict, protocol=pickle.HIGHEST_PROTOCOL)

        # Number of sequences classified to genus level
        def genus_seq_count(dict_class):
            genus_dict = {}
            for key, value in genus_taxid.items():
                seq_sum = 0
                for taxid in value:
                    if taxid in dict_class:
                        seq_sum += dict_class[taxid]
                genus_dict[key] = seq_sum
            return genus_dict

        genus_either = genus_seq_count(either_class)
        genus_combined = genus_seq_count(combined_class)

        table_summary = pd.DataFrame(columns=['Species', 'Species TaxID',
                                              'Species sequences',
                                              'Species sequences (stringent)',
                                              'Genus',
                                              'Genus sequences',
                                              'Genus sequences (stringent)'])
        table_summary['Species TaxID'] = [int(key) for key in levels_tax]
        table_summary['Species sequences'] = table_summary['Species TaxID'].map(either_class)
        table_summary['Species sequences (stringent)'] = table_summary['Species TaxID'].map(combined_class)
        table_summary['Species'] = table_summary['Species TaxID'].map(species_dict)
        table_summary['Genus'] = table_summary['Species TaxID'].map(genus_per_species)
        # Using functions in map to set default value of 0,
        # can use a defaultdict or Counter if have pandas 0.20 onwards
        table_summary['Genus sequences'] = table_summary['Genus'].map(lambda g: genus_either.get(g, 0))
        table_summary['Genus sequences (stringent)'] = table_summary['Genus'].map(lambda g: genus_combined.get(g, 0))
        # TODO correct virus_table
        if hasattr(table_summary, 'sort_values'):
            # pandas 0.17 onwards
            table_summary.sort_values(['Species sequences (stringent)', 'Species sequences'],
                                      ascending=False, inplace=True)
        else:
            table_summary.sort(['Species sequences (stringent)', 'Species sequences'],
                               ascending=False, inplace=True)
        table_summary.to_csv(os.path.join(out_dir, 'virus_table.txt'),
                             sep='\t', index=False)

    virusSummary(kodoja)
