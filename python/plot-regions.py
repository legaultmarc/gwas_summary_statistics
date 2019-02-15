#!/usr/bin/env python

import argparse
import gzip
import os
import fnmatch
import subprocess
import json
import urllib.request
from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


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


def get_genes_in_region(region):
    url = ("http://grch37.rest.ensembl.org/overlap/region/human/"
           "{}?feature=gene;content-type=application/json"
           "".format(region))

    try:
        with urllib.request.urlopen(url) as f:
            results = json.load(f)
    except:
        return df.DataFrame()

    df = pd.DataFrame.from_records(results)

    # ["symbol", "start", "end", "strand", "biotype"] are expected.
    df = df.rename(columns={
        "external_name": "symbol"
    })

    return df


def read_genetic_map(args):
    """Reads the genetic map."""
    data = pd.read_csv(
        "/data/projects/summary_statistics/utils/genetic_map.txt.gz",
        sep="\t", compression="gzip",
        dtype={"Chromosome": str}
    )

    # Sub-setting the data to get a region of X base pair on each side of the
    # hit
    chrom, start, end = parse_region(args.region)

    region = data["Chromosome"] == chrom
    region = region & (data["Position(bp)"] >= start)
    region = region & (data["Position(bp)"] <= end)

    data = data[region]

    return data


def plot_stats(stats_filename, output_dir, args):
    plot_filename = os.path.join(output_dir, "region.png")

    stats = pd.read_csv(stats_filename, sep="\t")

    fig = plt.figure(figsize=(11, 4))

    # The axes
    recomb_ax = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
    genes_ax = plt.subplot2grid((4, 1), (3, 0), sharex=recomb_ax)
    assoc_ax = recomb_ax.twinx()

    _plot_recomb(recomb_ax, args)
    _plot_genes(fig, genes_ax, args)
    _plot_assoc(assoc_ax, stats, args)

    plt.savefig(plot_filename, dpi=600, bbox_inches='tight')


def parse_region(region):
    chrom, tail = region.split(":")
    start, end = (int(i) for i in tail.split("-"))

    return chrom, start, end


def _plot_genes(fig, ax, args):
    chrom, start, end = parse_region(args.region)

    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("none")
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 9))
    ax.spines["right"].set_visible(False)
    ax.axes.get_yaxis().set_visible(False)

    ax.tick_params(axis="both", which="major", labelsize=8)

    genes = get_genes_in_region(args.region)
    genes = genes.sort_values(by=["start", "end"])

    renderer = fig.canvas.get_renderer()

    last_t_obj = {}
    last_end = defaultdict(int)

    for i in range(genes.shape[0]):
        gene_start = genes.iloc[i, :].start
        gene_end = genes.iloc[i, :].end
        gene_name = genes.iloc[i, :].symbol

        # Checking the starting position of the gene
        if gene_start < start:
            gene_start = start
        gene_start /= 1e6

        # Checking the ending position of the gene
        if gene_end > end:
            gene_end = end
        gene_end /= 1e6

        # Updating the gene label
        gene_label = None
        if genes.iloc[i, :].strand == 1:
            gene_label = gene_name + r"$\rightarrow$"
        else:
            gene_label = r"$\leftarrow$" + gene_name

        # We find the first j where we can put the line
        j = 0
        while True:
            if last_end[j] < gene_start:
                break
            j -= 1

        # Trying to put the label there
        t = ax.text((gene_start + gene_end) / 2, j - 0.15, gene_label,
                    fontsize=5, ha="center", va="top")

        # Is there a bbox in this location?
        if j in last_t_obj:
            # Getting the bbox
            bb = t.get_window_extent(renderer=renderer)
            last_bb = last_t_obj[j].get_window_extent(renderer=renderer)

            while last_bb.overlaps(bb):
                # BBoxes overlap
                j -= 1
                t.set_y(j - 0.15)

                if j not in last_t_obj:
                    break

                # Need to update both bboxes
                bb = t.get_window_extent(renderer=renderer)
                last_bb = last_t_obj[j].get_window_extent(renderer=renderer)

        marker = "-"
        other_param = {}
        if (gene_end - gene_start) < 3e-3:
            # Too small
            marker = "s"
            other_param["ms"] = 1.8
        ax.plot([gene_start, gene_end], [j, j], marker, lw=2,
                color="#000000", clip_on=False, **other_param)

        # Saving the last position (last end and bbox)
        last_end[j] = gene_end + 3e-3
        last_t_obj[j] = t

        min_y, max_y = ax.get_ylim()
        if min_y == -1:
            ax.set_ylim(-1.5, max_y)

        # Setting the ticks below the X axis for genes
        ax.get_xaxis().set_tick_params(direction='out')

        ax.set_xlabel("Position on chr{} (Mb)".format(chrom), fontsize=10,
                      weight="normal")


def _plot_recomb(ax, args):
    chrom, start, end = parse_region(args.region)

    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("right")
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_position(("outward", 9))
    ax.spines["bottom"].set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    ax.yaxis.set_label_position("right")

    ax.tick_params(axis="both", which="major", labelsize=8)

    ax.set_ylabel("Recombination Rate (cM/Mb)", fontsize=10,
                  weight="normal", rotation=270, va="bottom")

    genetic_map = read_genetic_map(args)

    ax.plot(genetic_map["Position(bp)"] / 1e6,
            genetic_map["Rate(cM/Mb)"],
            "-", lw=1, color="black", clip_on=False)

    ax.set_ylim(0, 100)
    ax.set_xlim(start/1e6, end/1e6)
    

def _plot_assoc(ax, stats, args):
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 9))
    ax.yaxis.set_label_position("left")

    ax.tick_params(axis='both', which='major', labelsize=8)

    ax.plot(
        stats.pos / 1e6,
        -np.log10(stats.p),
        ".", clip_on=False,
        ms=6
    )

    ax.axhline(-np.log10(args.significance_threshold), ls="--",
               color="#000000", lw=1)

    ax.set_ylabel("$-\log(p)$")


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

        plot_stats(stats_filename, output_dir, args)

        # Launch region-plot on the extracted regions.
        # command = [
        #     "launch-region-plot",
        #     "--assoc", stats_filename,
        #     "--genotypes", args.reference_genotypes,
        #     "--genotypes-format", args.reference_genotypes_format,
        #     "--plot-format", "png",
        #     "--significant", str(args.significance_threshold),
        #     "--snp-col", "study_marker_name",
        #     "--chr-col", "chrom",
        #     "--pos-col", "pos",
        #     "--p-col", "p",
        #     "--a1-col", "coded_allele",
        #     "--a2-col", "reference_allele",
        #     "--genetic-map", "/shares/data/works/lemieuxl/data/genetic_map.txt.gz",
        #     "--genetic-chr-col", "Chromosome",
        #     "--genetic-pos-col", "Position(bp)",
        #     "--genetic-rate-col", "Rate(cM/Mb)",
        #     "--output-directory", output_dir
        # ]

        # print(command)

        # subprocess.check_call(command)



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
