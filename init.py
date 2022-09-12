import os
from collections import defaultdict
from random import sample
import pandas as pd
import argparse

def add_sample(sample_dict, sample_id, header, path):
    if sample_id not in sample_dict:
        sample_dict[sample_id][header]=path

def get_samples(path):
    """
    create a table containing the fastq files
    """
    samples = defaultdict(dict)
    for root, dirs, fs in os.walk(os.path.abspath(path)):
        for fq in fs:

            # only check fq files
            if "fastq" in fq or "fq" in fq:
                fq_path = os.path.join(root, fq)
                sample_id = fq.split(".")[0]
                if "_R1" == sample_id[-3:]:
                    fq_path = fq_path + "," + fq_path.replace("_R1", "_R2")
                    sample_id = sample_id[:-3]
                if "_R2" == sample_id[-3:]:
                    fq_path = fq_path.replace("_R2", "_R1") + "," + fq_path
                    sample_id = sample_id[:-3]
                add_sample(samples, sample_id, "fqs", fq_path)
    samples_dt = pd.DataFrame(samples).T
    return samples_dt

def parse_arguments():
    """Read arguments from the console"""
    parser = argparse.ArgumentParser(
        description="Note: generate sample.tsv",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("-i", "--illumina", help='path to illumina fastq files')
    parser.add_argument("-n", "--nanopore", help='path to nanopore fastq files')
    parser.add_argument("-o", "--out", help='path to working directory', default=os.getcwd())

    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    # if args exists, then use it
    if args.illumina and not args.nanopore:
        sample_dt = get_samples(args.illunima)
        # add new column: type, illumina
        sample_dt["type"] = "illumina"
    if args.nanopore and not args.illumina:
        sample_dt = get_samples(args.nanopore)
        sample_dt["type"] = "nanopore"
    if args.illumina and args.nanopore:
        sample_dt_illu = get_samples(args.illumina)
        sample_dt_illu["type"] = "illumina"
        sample_dt_nano = get_samples(args.nanopore)
        sample_dt_nano["type"] = "nanopore"
        sample_dt = pd.concat([sample_dt_illu, sample_dt_nano])
    sample_dt.to_csv(args.out + "/samples.tsv", sep="\t")

if __name__ == "__main__":
    main()
