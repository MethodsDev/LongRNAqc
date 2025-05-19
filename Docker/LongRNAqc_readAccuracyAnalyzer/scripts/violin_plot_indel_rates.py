#!/usr/bin/env python3

import os
import sys
import pysam
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter
import math
from dataclasses import dataclass


@dataclass
class ReadAlignmentStats:
    cigar_length: int
    match_length: int
    mismatch_length: int
    insertion_count: int
    insertion_length: int
    deletion_count: int
    deletion_length: int
    soft_clip: int
    hard_clip: int


def parse_bam_file_paths(sample_names_str, bam_files_str):
    sample_names = sample_names_str.split(',')
    bam_files = bam_files_str.split(',')
    if len(sample_names) != len(bam_files):
        print("Error: The number of sample names must match the number of BAM files.")
        sys.exit(1)
    return {sample_name: bam_file for sample_name, bam_file in zip(sample_names, bam_files)}


def get_read_stats_from_bam(bam_file):
    read_to_stats = {}
    with pysam.AlignmentFile(bam_file, "rb") as inbam:
        for read in inbam.fetch():
            if read.is_unmapped or read.is_secondary:
                continue
            total_cigar_length = 0
            total_match_length = 0
            total_mismatch_length = 0
            total_insertion_count = 0
            total_insertion_length = 0
            total_deletion_count = 0
            total_deletion_length = 0
            total_soft_clip_length = 0
            total_hard_clip_length = 0
            
            for op, op_length in read.cigartuples:
                if op == 0 or op == 7:  # M or = in CIGAR
                    total_cigar_length += op_length
                    total_match_length += op_length
                elif op == 8:  # X in CIGAR
                    total_mismatch_length += op_length
                    total_cigar_length += op_length
                elif op == 1:  # I in CIGAR
                    total_insertion_length += op_length
                    total_insertion_count += 1
                    # total_cigar_length += op_length
                elif op == 2:  # D in CIGAR
                    total_deletion_count += 1
                    total_deletion_length += op_length
                elif op == 4:  # S in CIGAR
                    total_soft_clip_length += op_length
                    # total_cigar_length += op_length
                elif op == 5:  # H in CIGAR
                    total_hard_clip_length += op_length
            
            read_stats = ReadAlignmentStats(total_cigar_length, total_match_length, total_mismatch_length, total_insertion_count, total_insertion_length, total_deletion_count, total_deletion_length, total_soft_clip_length, total_hard_clip_length)
            yield read_stats


def weighted_boxplot_stats(sorted_values, sorted_counts):
    """
    From a sorted array of `sorted_values` and parallel `sorted_counts`,
    compute typical boxplot statistics:
      Q1, median, Q3, lower_whisker, upper_whisker, outliers

    Returns a dict with those keys.
    """
    total = sum(sorted_counts)
    if total == 0:
        # no data
        return {
            "q1": None, "median": None, "q3": None,
            "lower_whisker": None, "upper_whisker": None,
            "outliers": []
        }

    # find Q1, median, Q3 position along cumulative counts
    q1_cut = 0.25 * total
    median_cut = 0.50 * total
    q3_cut = 0.75 * total

    cumsum = 0
    q1 = None
    median = None
    q3 = None

    for val, cnt in zip(sorted_values, sorted_counts):
        cumsum_prev = cumsum
        cumsum += cnt

        # if we cross q1_cut in this segment, pick val as Q1
        if q1 is None and cumsum >= q1_cut:
            q1 = val
        # if we cross median_cut, pick val as median
        if median is None and cumsum >= median_cut:
            median = val
        # if we cross q3_cut, pick val as Q3
        if q3 is None and cumsum >= q3_cut:
            q3 = val
        if q1 is not None and median is not None and q3 is not None:
            break

    # compute interquartile range (iqr) and whiskers
    iqr = q3 - q1
    lower_whisker = q1 - 1.5 * iqr
    upper_whisker = q3 + 1.5 * iqr

    # we define outliers as values < lower_whisker or > upper_whisker.
    outliers = []
    for val, cnt in zip(sorted_values, sorted_counts):
        if val < lower_whisker or val > upper_whisker:
            # outliers.extend([val] * cnt)
            outliers.append(val)

    return {
        "q1": q1,
        "median": median,
        "q3": q3,
        "lower_whisker": lower_whisker,
        "upper_whisker": upper_whisker,
        "outliers": outliers
    }


def overlay_boxplot_on_violin(
    ax,
    sorted_values,
    sorted_counts,
    x_center=0.0,
    box_width=0.3,
    color="k",
    whisker_linewidth=1,
    median_linewidth=2,
    box_linewidth=1,
    outlier_marker_size=5,
    clip_y=None
):
    """
    Overlays a basic boxplot on the Axes at x_center, using
    the weighted distribution from (sorted_values, sorted_counts).
    
    - The box is drawn with no fill (facecolor="none"), 
      only an outline (edgecolor=color).
    - Median, whiskers, and outliers are drawn as opaque lines/points.

    clip_y : tuple or None
        If given, (clip_min, clip_max). Box, whiskers, and outliers
        will be clamped to this range.
    """
    stats = weighted_boxplot_stats(sorted_values, sorted_counts)
    if stats["q1"] is None:
        return  # no data to plot

    q1 = stats["q1"]
    med = stats["median"]
    q3 = stats["q3"]
    lw = stats["lower_whisker"]
    uw = stats["upper_whisker"]
    outliers = stats["outliers"]

    # if clip_y is set, clamp the box/whisker values and filter outliers
    if clip_y is not None:
        y_min, y_max = clip_y

        q1 = max(y_min, min(q1, y_max))
        med = max(y_min, min(med, y_max))
        q3 = max(y_min, min(q3, y_max))
        lw = max(y_min, min(lw, y_max))
        uw = max(y_min, min(uw, y_max))

        # filter outliers to those within [y_min, y_max]
        outliers = [val for val in outliers if y_min <= val <= y_max]


    # compute left/right x-coords of the box
    x_left = x_center - box_width/2
    x_right = x_center + box_width/2
    box_height = q3 - q1

    # draw the box (unfilled rectangle)
    rect = patches.Rectangle(
        (x_left, q1),
        box_width,        # width
        box_height,       # height
        linewidth=box_linewidth,
        edgecolor=color,
        facecolor="none"  # no fill (fully transparent)
    )
    ax.add_patch(rect)

    # whiskers
    #    vertical line from lower_whisker to Q1
    ax.plot(
        [x_center, x_center], 
        [lw, q1],
        color=color,
        linewidth=whisker_linewidth
    )
    #    vertical line from Q3 to upper_whisker
    ax.plot(
        [x_center, x_center],
        [q3, uw],
        color=color,
        linewidth=whisker_linewidth
    )
    #    small horizontal ticks at lw and uw
    ax.plot(
        [x_center - box_width/4, x_center + box_width/4],
        [lw, lw],
        color=color,
        linewidth=whisker_linewidth
    )
    ax.plot(
        [x_center - box_width/4, x_center + box_width/4],
        [uw, uw],
        color=color,
        linewidth=whisker_linewidth
    )
        
    # median line (horizontal inside the box)
    ax.plot(
        [x_left, x_right], 
        [med, med], 
        color="red", 
        linewidth=median_linewidth
    )

    # outliers as scatter points (opaque)
    if outliers:
        ax.scatter(
            [x_center]*len(outliers),
            outliers,
            color=color,
            s=outlier_marker_size,
            alpha=1.0  # fully opaque
        )


def plot_weighted_violin(ax, values, weights, x_center=0.0, width=0.6,
                         color="blue", alpha=0.4, clip=(0,1), bw_adjust=0.2, gridsize=1000):
    """
    Plot a 'violin' by using sns.kdeplot() with (values, weights)
    and mirroring the density around x_center.

    Parameters:
    -----------
    ax : Matplotlib Axes
        Where to draw the violin
    values : 1D array-like
        Data values in [clip[0], clip[1]]
    weights : 1D array-like
        Counts or weights for each value in `values`
    x_center : float
        The horizontal center of the violin
    width : float
        Maximum full width of the violin. The peak of the density
        will be scaled to width/2 in either direction.
    color : str
        Face color of the violin
    alpha : float
        Transparency
    clip : tuple (low, high)
        Clip range for the KDE
    bw_adjust : float
        Bandwidth adjustment for KDE (like in Seaborn)

    Returns:
    --------
    None, draws on `ax`.
    """

    # if there's only one unique bin
    if len(values) == 1:
        val = values[0]
        cnt = weights[0]

        # for a single-value distribution, draw a tiny rectangle.
        # half the normal "width" and a small vertical extent.
        half_width = width / 4.0
        half_height = 0.005

        # make a small rectangle around (x_center, val)
        x_left  = x_center - half_width
        x_right = x_center + half_width
        y_low   = val - half_height
        y_high  = val + half_height

        ax.fill_betweenx(
            [y_low, y_high],
            x_left,
            x_right,
            color=color,
            alpha=alpha
        )

        return

    # create a temporary figure/axes to run sns.kdeplot and extract the line data
    fig_tmp, ax_tmp = plt.subplots()
    
    # for vertical orientation `y=values`.
    p = sns.kdeplot(
        y=values,
        weights=weights,
        fill=False,          # fill manually
        clip=clip,
        bw_adjust=bw_adjust,
        gridsize=gridsize,   # higher => more resolution
        common_norm=False,   # each distribution handled separately
        cut=0,               # don't extend beyond data range
        ax=ax_tmp
    )

    # extract the line data for the KDE from the temporary Axes
    # usually p.lines[0] is the main density line
    if not p.lines:
        # if there's no density line (all zero?), just close and return
        plt.close(fig_tmp)
        return

    line = p.lines[0]
    x_data = line.get_xdata()  # "density" dimension
    y_data = line.get_ydata()  # "value" dimension

    plt.close(fig_tmp)  # don't need the temporary figure anymore

    # scale x_data so its maximum matches width/2
    max_density = x_data.max()
    if max_density > 0:
        scaling = (width / 2.0) / max_density
    else:
        scaling = 1.0
    x_data_scaled = x_data * scaling

    # on the real Axes, fill between x_center - x_data_scaled and x_center + x_data_scaled with y_data as the vertical axis
    ax.fill_betweenx(
        y_data,
        x_center - x_data_scaled,
        x_center + x_data_scaled,
        color=color,
        alpha=alpha,
    )


def plot_proportions_violin_plots(results, output_path_prefix):
    
    stats_order = [
        "match_proportion",
        "mismatch_proportion",
        "softclip_proportion",
        "insertion_rate",
        "deletion_rate",
        "indel_rate",
        "insertion_proportion",
        "deletion_proportion",
        "indel_proportion"
    ]
    
    titles = [
        "Match Proportion",
        "Mismatch Proportion",
        "Softclip Proportion",
        "Insertion Rate",
        "Deletion Rate",
        "Indel Rate",
        "Insertion Proportion",
        "Deletion Proportion",
        "Indel Proportion"
    ]
    
    # fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(18, 22))
    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(18, 33))
    axs = axs.flatten()  # make it easier to iterate through the 9 axes in a single list
    
    for ax, stat_name, title in zip(axs, stats_order, titles):
        # plot each run_name (BAM) as one violin along x = i
        for i, (run_name, stats_dict) in enumerate(results.items()):
            bin_dict = stats_dict[stat_name]
            # convert keys/values to sorted arrays
            keys_list = sorted(bin_dict.keys())
            values_list = [bin_dict[k] for k in keys_list]
            
            plot_weighted_violin(
                ax=ax,
                values=keys_list,
                weights=values_list,
                x_center=i,
                width=0.8,
                color="blue",
                alpha=0.4
            )

            overlay_boxplot_on_violin(
                ax=ax,
                sorted_values=keys_list,
                sorted_counts=values_list,
                x_center=i,
                box_width=0.3,  # smaller than the violin width
                color="black",
                clip_y=(0,1)
            )
        
        ax.set_title(title)
        # put run_name labels on the x-axis
        x_positions = range(len(results))
        run_names = list(results.keys())
        ax.set_xticks(x_positions)
        ax.set_xticklabels(run_names, rotation=90)
    
        # if stat_name in ("match_proportion", "mismatch_proportion", 
        #                  "softclip_proportion", "insertion_rate", 
        #                  "deletion_rate", "indel_rate"):
        if stat_name in ("match_proportion", "mismatch_proportion", "softclip_proportion"):
            ax.set_ylim(-0.05, 1.05)  
    
    plt.tight_layout()
    plt.savefig((output_path_prefix + '_violin_plots.png'), dpi=300)


def plot_phred_violin_plots(results, output_path_prefix):
    # fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    fig, axs = plt.subplots(2, 3, figsize=(18, 12))
    axs = axs.flatten()  # make it easier to iterate through the 6 axes in a single list
    
    phred_metrics = [
        "insertion_rate_phred",
        "deletion_rate_phred",
        "indel_rate_phred",
        "insertion_proportion_phred",
        "deletion_proportion_phred",
        "indel_proportion_phred"
    ]

    titles = [
        "Insertion Rate Phred Score",
        "Deletion Rate Phred Score",
        "Indel Rate Phred Score",
        "Insertion Proportion Phred Score",
        "Deletion Proportion Phred Score",
        "Indel Proportion Phred Score"
    ]

    for ax, metric, title in zip(axs, phred_metrics, titles):
        for i, (run_name, stats_dict) in enumerate(results.items()):
            phred_dict = stats_dict[metric]
            if not phred_dict:
                continue

            keys_list = sorted(phred_dict.keys())
            values_list = [phred_dict[k] for k in keys_list]

            plot_weighted_violin(
                ax=ax,
                values=keys_list,
                weights=values_list,
                x_center=i,
                width=0.8,
                color="blue",
                alpha=0.4,
                clip=(0, 40),
            )
            
            overlay_boxplot_on_violin(
                ax=ax,
                sorted_values=keys_list,
                sorted_counts=values_list,
                x_center=i,
                box_width=0.3,  # smaller than the violin width
                color="black",
                clip_y=(0,40)
            )

        ax.set_title(title)
        ax.set_ylim(0, 40)  # If you are capping at 40
        x_positions = range(len(results))
        ax.set_xticks(x_positions)
        ax.set_xticklabels(list(results.keys()), rotation=90)

    plt.tight_layout()
    plt.savefig((output_path_prefix + '_phred_violin_plots.png'), dpi=300)


def build_phred_dict(rate_dict, cap=40.0):
    """
    Given a dictionary of {rate -> count}, return a dictionary {phred -> count},
    capping at `cap` for the case rate=0 and for very small rates
    """
    phred_dict = Counter()
    for rate_val, cnt in rate_dict.items():
        if rate_val > 0:
            phred_score = -10.0 * math.log10(rate_val)
            if phred_score > cap:
                phred_score = cap
        else:
            phred_score = cap

        # round, though the input should already be rounded
        phred_score_rounded = round(phred_score, 2)
        phred_dict[phred_score_rounded] += cnt
    return phred_dict


def calculate_stats(bam_file_dict):

    results = {}

    for run_name, bam_file in bam_file_dict.items():

        match_proportion = Counter()
        mismatch_proportion = Counter()
        softclip_proportion = Counter()
        insertion_rate = Counter()
        insertion_proportion = Counter()
        deletion_rate = Counter()
        deletion_proportion = Counter()
        indel_rate = Counter()
        indel_proportion = Counter()

        for stats in get_read_stats_from_bam(bam_file):
            match_proportion[round((stats.match_length / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1
            mismatch_proportion[round((stats.mismatch_length / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1
            softclip_proportion[round((stats.soft_clip / (stats.cigar_length + stats.soft_clip) if stats.cigar_length > 0 else 0), 3)] += 1
            insertion_rate[round((stats.insertion_count / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1
            insertion_proportion[round((stats.insertion_length / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1
            deletion_rate[round((stats.deletion_count / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1
            deletion_proportion[round((stats.deletion_length / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1
            indel_rate[round(((stats.insertion_count + stats.deletion_count) / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1
            indel_proportion[round(((stats.insertion_length + stats.deletion_length) / stats.cigar_length if stats.cigar_length > 0 else 0), 3)] += 1

        insertion_rate_phred = build_phred_dict(insertion_rate, cap=40.0)
        insertion_proportion_phred = build_phred_dict(insertion_proportion, cap=40.0)
        deletion_rate_phred = build_phred_dict(deletion_rate, cap=40.0)
        deletion_proportion_phred = build_phred_dict(deletion_proportion, cap=40.0)
        indel_rate_phred = build_phred_dict(indel_rate, cap=40.0)
        indel_proportion_phred = build_phred_dict(indel_proportion, cap=40.0)

        results[run_name] = {
            "match_proportion": match_proportion,
            "mismatch_proportion": mismatch_proportion,
            "softclip_proportion": softclip_proportion,
            "insertion_rate": insertion_rate,
            "insertion_proportion": insertion_proportion,
            "deletion_rate": deletion_rate,
            "deletion_proportion": deletion_proportion,
            "indel_rate": indel_rate,
            "insertion_rate_phred": insertion_rate_phred,
            "insertion_proportion_phred": insertion_proportion_phred,
            "deletion_rate_phred": deletion_rate_phred,
            "deletion_proportion_phred": deletion_proportion_phred,
            "indel_rate_phred": indel_rate_phred,
            "indel_proportion_phred": indel_proportion_phred,
        }

    return results


def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py output_prefix <sample_name1,sample_name2,...> <bam_file1,bam_file2,...>")
        sys.exit(1)

    output_prefix_str, sample_names_str, bam_files_str = sys.argv[1], sys.argv[2], sys.argv[3]
    bam_file_dict = parse_bam_file_paths(sample_names_str, bam_files_str)

    results = calculate_stats(bam_file_dict)
    
    plot_proportions_violin_plots(results, output_prefix_str)
    plot_phred_violin_plots(results, output_prefix_str)


if __name__ == "__main__":
    main()

