#!/usr/bin/env python
"""Script for running Kodoja pipeline modules."""

import argparse
import os
import sys
import time

from diagnostic_modules import version
from diagnostic_modules import check_path
from diagnostic_modules import test_format
from diagnostic_modules import check_file
from diagnostic_modules import fastqc_trim
from diagnostic_modules import kraken_classify
from diagnostic_modules import kaiju_classify
from diagnostic_modules import add_krona_representation
from diagnostic_modules import result_comparaison

help_text = """Kodoja Search is a tool intended to identify viral sequences
in a FASTQ/FASTA sequencing run by matching them against both Kraken and
Kaiju databases.
"""

help_epilog = """
The main output of ``kodoja_search.py`` is a file called ``virus_table.txt``
in the specified output directory. This is a plain text tab-separated table,
the columns are as follows:

1. Species name,
2. Species NCBI taxonomy identifier (TaxID),
3. Number of reads assigned by *either* Kraken or Kaiju to this species,
4. Number of Reads assigned by *both* Kraken and Kaiju to this species,
5. Genus name,
6. Number of reads assigned by *either* Kraken or Kaiju to this genus,
7. Number of reads assigned by *both* Kraken and Kaiju to this genus.

The output directory includes additional files, including ``kodoja_VRL.txt``
(a table listing the read identifiers used) which is intended mainly as
input to the ``kodoja_retrieve.py`` script.

See also https://github.com/abaizan/kodoja/wiki/Kodoja-Manual
"""

parser = argparse.ArgumentParser(description=help_text,
                                 epilog=help_epilog,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--version',
                    action='version',
                    version='Kodoja v' + version)
parser.add_argument('-o', '--output_dir', type=str, required=True,
                    help='Output directory path, required')
parser.add_argument('-d1', '--kraken_db', type=str, required=True,
                    help='Kraken database path, required')
parser.add_argument('-d2', '--kaiju_db', type=str, required=True,
                    help='Kaiju database path, required')
parser.add_argument('-r1', '--read1', type=str, required=True,
                    help='Read 1 file path, required')
parser.add_argument('-r2', '--read2', type=str, default=False,
                    help='Read 2 file path')
parser.add_argument('-f', '--data_format', type=str, default='fastq',
                    help='Sequence data format (default fastq)')
parser.add_argument('-t', '--threads', type=int, default=1,
                    help='Number of threads (default 1)')
parser.add_argument('-s', '--host_subset', type=int, default=False,
                    help='Subset sequences with this tax id from results')
parser.add_argument('-m', '--trim_minlen', type=int, default=50,
                    help='Trimmomatic minimum length')
parser.add_argument('-a', '--trim_adapt', type=str, default=False,
                    help='Illumina adapter sequence file')
parser.add_argument('-q', '--kraken_quick', type=int, default=False,
                    help='Number of minium hits by Kraken')
parser.add_argument('-c', '--kaiju_score', type=int, default=85,
                    help='Kaju alignment score')
parser.add_argument('-l', '--kaiju_minlen', type=int, default=15,
                    help='Kaju minimum length')
parser.add_argument('-i', '--kaiju_mismatch', type=int, default=1,
                    help='Kaju allowed mismatches')
parser.add_argument('-mm', '--memory', default=False, action='store_true',
                    help='kraken --memory-mapping for big database')                
args = parser.parse_args()

# Check that dirs have "/" at the end  
args.output_dir += check_path(args.output_dir)
args.kraken_db += check_path(args.kraken_db)
args.kaiju_db += check_path(args.kaiju_db)

# Check args.output_dir exits else make dir
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# os.chdir(args.output_dir)


# Write a log_file:
log_filename = os.path.join(args.output_dir, "log_file.txt")
with open(log_filename, "w") as log_file:
    log_file.write("General parameters:\n"
                   "file1 = %s\n"
                   "file2 = %s\n"
                   "output directory = %s\n"
                   "threads = %i\n"
                   "host_subset = %s\n"
                   "Trimmomatic parameters:\n"
                   "trim_minlen = %s\n"
                   "Kraken parameters:\n"
                   "kraken database = %s\n"
                   "quick_minhits = %s\n"
                   "Kaiju parameters:\n"
                   "args.kaiju_db = %s\n"
                   "kaiju_minlen = %i\n"
                   "kaiju_score = %i\n"
                   "kaiju_mismatch = %i\n"
                   "memory = %s\n"
                   % (args.read1, args.read2, args.output_dir,
                      args.threads, args.host_subset, args.trim_minlen,
                      args.kraken_db, args.kraken_quick,
                      args.kaiju_db, args.kaiju_minlen, args.kaiju_score,
                      args.kaiju_mismatch, args.memory))

# TODO: Review this and consider using Python standard library's logging module
def log(message):
    """Append the message to the log file."""
    with open(log_filename, "a") as log_file:
        log_file.write(message)


def main():
    """Run kodoja search using command line options."""
    t0 = time.time()

    # Test format, change seqIDs and check paired files are correct
    test_format(args.read1, args.data_format)
    check_file(args.read1, args.output_dir, args.data_format, args.read2)

    t1 = time.time()

    # Set all variables
    initial_file1 = os.path.join(args.output_dir, 'renamed_file_1.' + args.data_format)
    kraken_file1 = kaiju_file1 = os.path.join(args.output_dir, "trimmed_read1")

    if args.read2:
        # Set tool files
        kraken_file2 = kaiju_file2 = os.path.join(args.output_dir, "trimmed_read2")
        initial_file2 = os.path.join(args.output_dir, 'renamed_file_2.' + args.data_format)
    else:
        kraken_file2 = kaiju_file2 = False
        initial_file2 = False

    if args.data_format == "fastq":
        # fasta files cannot be QC'd - only for fastq files
        log("Starting FASTQ read trimming\n")
        fastqc_trim(args.output_dir, initial_file1, args.trim_minlen, args.threads,
                    args.trim_adapt, initial_file2)
    else:
        kraken_file1 = kaiju_file1 = initial_file1
        kraken_file2 = kaiju_file2 = initial_file2


    t2 = time.time()
    # Kraken classification
    log("Starting Kraken classification\n")
    kraken_classify(args.output_dir, kraken_file1, args.threads,
                    args.kraken_db, args.memory, kraken_file2,
                    args.kraken_quick)
    t3 = time.time()
    # Kaiju classification of all sequences or subset sequences
    log("Starting Kaiju classification\n")
    kaiju_classify(kaiju_file1, args.threads, args.output_dir,
                   args.kaiju_db, args.kaiju_minlen, args.kraken_db,
                   kaiju_file2, kaiju_mismatch=args.kaiju_mismatch,
                   kaiju_score=args.kaiju_score)
    t4 = time.time()
    # Make krona representation
    log("Make krona html display\n")
    add_krona_representation(args.output_dir)
    t5 = time.time()

    # Merge results
    log("Analyzing Kraken and Kaiju results\n")
    result_comparaison(args.output_dir, "kraken_table.txt",  "kaiju_table.txt", args.host_subset)
    t6 = time.time()

    # Create log file
    if args.host_subset:
        print_statment = "comparing result and subset sequences = %0.1f h\n" % ((t6 - t5) / 3600)
    else:
        print_statment = "comparing result data = %0.1f h" % ((t6 - t5) / 3600)

    log("Script timer:\n"
        "testing format/replace seqID = %0.1f s\n"
        "fastq and trim = %0.1f min\n"
        "kraken classification = %0.1f h\n"
        "kaiju classification = %0.1f h\n"
        "krona html representation = %0.1f m\n"
        "%s\n"
        "total = %0.1f h\n"
        % (t1 - t0,
           (t2 - t1) / 60,
           (t3 - t2) / 3600,
           (t4 - t3) / 3600,
           (t5 - t4) / 60,
           print_statment,
           (t6 - t0) / 3600))

    log("\nkodoja_search.py finished sucessfully.\n")


try:
    main()
except KeyboardInterrupt:
    msg = "Kodoja was interupted by the user.\n"
    log(msg)
    sys.exit(msg)
except Exception:
    import traceback
    msg = ("Kodoja failed unexpectedly with the following:\n"
           "\n"
           "%s\n"
           "\n"
           "If this happens consistently, please check you are using the\n"
           "latest version, then check the issue tracker to see if this is\n"
           "a known issue, and if not please report the problem:\n"
           "\n"
           "https://github.com/abaizan/kodoja/issues\n"
           "\n"
           "Kodoja aborted.\n" % traceback.format_exc())
    log(msg)
    sys.exit(msg)
