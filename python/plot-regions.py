#!/usr/bin/env python

import argparse
import gzip
import os
import fnmatch
import subprocess


def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


def detect_available_summary_statistics(root):
    indices = find("*.formatted*.tbi", root)
    filenames = [i[:-4] for i in indices]
    dataset_names = [i.split("/")[-1].split(".formatted")[0] for i in filenames]

    return {k: v for k, v in zip(dataset_names, filenames)}

def main():
    args = parse_args()

    os.mkdir(args.region)

    # The available summary stats are in /data/projects
    datasets = detect_available_summary_statistics(
        args.summary_statistics_root
    )

    for name, filename in datasets.items():
        # Prepare the output directory.
        output_dir = os.path.join(args.region, name)
        os.mkdir(output_dir)

        # Extract the region and prepare the filename.

        # Read the header or else it's thrown out by tabix.
        with gzip.open(filename, "rb") as f:
            header = f.readline()

        # Extract the statistics.
        stats_filename = os.path.join(
            output_dir, "statistics_{}.tsv".format(args.region)
        )

        with open(stats_filename, "wb") as f:
            f.write(header)

            with subprocess.Popen(
                ["tabix", filename, args.region],
                stdout=subprocess.PIPE
            ) as proc:
                f.write(proc.communicate()[0])

        # Launch region-plot on the extracted regions.
        command = [
            "launch-region-plot",
            "--assoc", stats_filename,
            "--genotypes", args.reference_genotypes,
            "--genotypes-format", args.reference_genotypes_format,
            "--plot-format", "png",
            "--significant", str(args.significance_threshold),
            "--genetic-map", "/shares/data/works/lemieuxl/data/genetic_map.txt.gz",
            "--genetic-chr-col", "Chromosome",
            "--genetic-pos-col", "Position(bp)",
            "--genetic-rate-col", "Rate(cM/Mb)",
            "--output-directory", output_dir
        ]

        print(command)

        subprocess.check_call(command)



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "region",
        help="A region of the form: CHR:START-END.",
        type=str
    )

    parser.add_argument(
        "--summary-statistics-root",
        help="Root directory where the look for summary statistics files.",
        default="/data/projects/summary_statistics"
    )

    parser.add_argument(
        "--reference-genotypes",
        help="Path to the reference panel used for LD computation.",
        default="/shares/data/works/lemieuxl/data/1000G_phase3/EUR.phase3_shapeit2_mvncall_integrated_v5.20130502.SNPs.founders",
    )

    parser.add_argument(
        "--reference-genotypes-format",
        help=("Format for the genotype file. This will be passed to "
              "launch-region-plot"),
        default="plink",
    )

    parser.add_argument(
        "--significance-threshold",
        default=0.001,
        type=float,
        help="This is the position of the significance threshold line."
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()
