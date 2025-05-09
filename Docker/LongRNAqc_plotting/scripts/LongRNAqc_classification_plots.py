#!/usr/bin/env python3

import argparse
import os
import re
import gzip

# import pgzip
from collections import defaultdict, Counter

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.stats import gaussian_kde
import seaborn as sns
import colorcet as cc

import math
import pandas as pd


def smart_open(file_path, mode="rt", thread=8):
    # Open file in binary mode and check the first two bytes.
    with open(file_path, "rb") as f:
        signature = f.read(2)
    # If the file starts with the gzip magic number, use (p)gzip.open().
    if signature == b"\x1f\x8b":
        # return pgzip.open(file_path, mode=mode, thread=thread)
        return gzip.open(file_path, mode=mode)
    else:
        return open(file_path, mode=mode)


def parse_sqanti_classification(sqanti_files, read_lengths_outfile):
    results = {}

    read_lens_ofh = open(read_lengths_outfile, "wt")

    for sample_name, sqanti_file in sqanti_files.items():
        with smart_open(sqanti_file, mode="rt") as f:
            FSM_counts = Counter()
            FSM_alt3_counts = Counter()
            FSM_alt5_counts = Counter()
            FSM_alt35_counts = Counter()
            FSM_ref_counts = Counter()
            FSM_monoexon_counts = Counter()
            ISM_counts = Counter()
            ISM_3frag_counts = Counter()
            ISM_5frag_counts = Counter()
            ISM_intfrag_counts = Counter()
            ISM_intron_retention_counts = Counter()
            ISM_monoexon_counts = Counter()
            NIC_counts = Counter()
            NIC_annotated_junctions_counts = Counter()
            NIC_annotated_splice_sites_counts = Counter()
            NIC_intron_retention_counts = Counter()
            NIC_monoexon_counts = Counter()
            NIC_monoexon_intron_retention_counts = Counter()
            NNC_counts = Counter()
            NNC_intron_retention_counts = Counter()
            NNC_partial_unk_splice_site_counts = Counter()
            genic_genomic_counts = Counter()
            antisense_counts = Counter()
            fusion_counts = Counter()
            fusion_intron_retention_counts = Counter()
            fusion_multiexon_counts = Counter()
            intergenic_counts = Counter()
            genic_intron_counts = Counter()

            f.readline()  # skip header

            for line in f:
                (
                    isoform,
                    chrom,
                    strand,
                    length,
                    exons,
                    structural_category,
                    associated_gene,
                    associated_transcript,
                    ref_length,
                    ref_exons,
                    diff_to_TSS,
                    diff_to_TTS,
                    diff_to_gene_TSS,
                    diff_to_gene_TTS,
                    subcategory,
                    RTS_stage,
                    all_canonical,
                    min_sample_cov,
                    min_cov,
                    min_cov_pos,
                    sd_cov,
                    FL,
                    n_indels,
                    n_indels_junc,
                    bite,
                    iso_exp,
                    gene_exp,
                    ratio_exp,
                    FSM_class,
                    coding,
                    ORF_length,
                    CDS_length,
                    CDS_start,
                    CDS_end,
                    CDS_genomic_start,
                    CDS_genomic_end,
                    predicted_NMD,
                    perc_A_downstream_TTS,
                    seq_A_downstream_TTS,
                    dist_to_CAGE_peak,
                    within_CAGE_peak,
                    dist_to_polyA_site,
                    within_polyA_site,
                    polyA_motif,
                    polyA_dist,
                    polyA_motif_found,
                    ORF_seq,
                    ratio_TSS,
                ) = line.split("\t")

                print("\t".join([sample_name, length]), file=read_lens_ofh)

                # since reused for all assignments, cast only once here
                length = int(length)

                if structural_category == "full-splice_match":
                    FSM_counts[length] += 1
                    if subcategory == "reference_match":
                        FSM_ref_counts[length] += 1
                    elif subcategory == "alternative_5end":
                        FSM_alt5_counts[length] += 1
                    elif subcategory == "alternative_3end5end":
                        FSM_alt35_counts[length] += 1
                    elif subcategory == "alternative_3end":
                        FSM_alt3_counts[length] += 1
                    elif subcategory == "mono-exon":
                        FSM_monoexon_counts[length] += 1
                elif structural_category == "incomplete-splice_match":
                    ISM_counts[length] += 1
                    if subcategory == "3prime_fragment":
                        ISM_3frag_counts[length] += 1
                    elif subcategory == "5prime_fragment":
                        ISM_5frag_counts[length] += 1
                    elif subcategory == "internal_fragment":
                        ISM_intfrag_counts[length] += 1
                    elif subcategory == "intron_retention":
                        ISM_intron_retention_counts[length] += 1
                    elif subcategory == "mono-exon":
                        ISM_monoexon_counts[length] += 1
                elif structural_category == "novel_in_catalog":
                    NIC_counts[length] += 1
                    if subcategory == "combination_of_known_junctions":
                        NIC_annotated_junctions_counts[length] += 1
                    elif subcategory == "combination_of_known_splicesites":
                        NIC_annotated_splice_sites_counts[length] += 1
                    elif subcategory == "intron_retention":
                        NIC_intron_retention_counts[length] += 1
                    elif subcategory == "mono-exon":
                        NIC_monoexon_counts[length] += 1
                    elif subcategory == "mono-exon_by_intron_retention":
                        NIC_monoexon_intron_retention_counts[length] += 1
                elif structural_category == "novel_not_in_catalog":
                    NNC_counts[length] += 1
                    if subcategory == "intron_retention":
                        NNC_intron_retention_counts[length] += 1
                    elif subcategory == "at_least_one_novel_splicesite":
                        NNC_partial_unk_splice_site_counts[length] += 1
                elif structural_category == "genic":
                    genic_genomic_counts[length] += 1
                elif structural_category == "antisense":
                    antisense_counts[length] += 1
                elif structural_category == "fusion":
                    fusion_counts[length] += 1
                    if subcategory == "intron_retention":
                        fusion_intron_retention_counts[length] += 1
                    elif subcategory == "multi-exon":
                        fusion_multiexon_counts[length] += 1
                elif structural_category == "intergenic":
                    intergenic_counts[length] += 1
                elif structural_category == "genic_intron":
                    genic_intron_counts[length] += 1

            total_counts = (
                FSM_counts
                + ISM_counts
                + NIC_counts
                + NNC_counts
                + genic_genomic_counts
                + antisense_counts
                + fusion_counts
                + intergenic_counts
                + genic_intron_counts
            )

            results[sample_name] = {
                "total_counts": total_counts,
                "total_sum": total_counts.total(),
                "FSM_counts": FSM_counts,
                "FSM_sum": FSM_counts.total(),
                "FSM_alt3_counts": FSM_alt3_counts,
                "FSM_alt5_counts": FSM_alt5_counts,
                "FSM_alt35_counts": FSM_alt35_counts,
                "FSM_ref_counts": FSM_ref_counts,
                "FSM_monoexon_counts": FSM_monoexon_counts,
                "ISM_counts": ISM_counts,
                "ISM_sum": ISM_counts.total(),
                "ISM_3frag_counts": ISM_3frag_counts,
                "ISM_5frag_counts": ISM_5frag_counts,
                "ISM_intfrag_counts": ISM_intfrag_counts,
                "ISM_intron_retention_counts": ISM_intron_retention_counts,
                "ISM_monoexon_counts": ISM_monoexon_counts,
                "NIC_counts": NIC_counts,
                "NIC_sum": NIC_counts.total(),
                "NIC_annotated_junctions_counts": NIC_annotated_junctions_counts,
                "NIC_annotated_splice_sites_counts": NIC_annotated_splice_sites_counts,
                "NIC_intron_retention_counts": NIC_intron_retention_counts,
                "NIC_monoexon_counts": NIC_monoexon_counts,
                "NIC_monoexon_intron_retention_counts": NIC_monoexon_intron_retention_counts,
                "NNC_counts": NNC_counts,
                "NNC_sum": NNC_counts.total(),
                "NNC_intron_retention_counts": NNC_intron_retention_counts,
                "NNC_partial_unk_splice_site_counts": NNC_partial_unk_splice_site_counts,
                "genic_genomic_counts": genic_genomic_counts,
                "genic_genomic_sum": genic_genomic_counts.total(),
                "antisense_counts": antisense_counts,
                "antisense_sum": antisense_counts.total(),
                "fusion_counts": fusion_counts,
                "fusion_sum": fusion_counts.total(),
                "fusion_intron_retention_counts": fusion_intron_retention_counts,
                "fusion_multiexon_counts": fusion_multiexon_counts,
                "intergenic_counts": intergenic_counts,
                "intergenic_sum": intergenic_counts.total(),
                "genic_intron_counts": genic_intron_counts,
                "genic_intron_sum": genic_intron_counts.total(),
            }

    read_lens_ofh.close()

    return results


def parse_lraa_sqantilike_classification(lraa_files, read_lengths_outfile):
    results = {}

    read_lens_ofh = open(read_lengths_outfile, "wt")

    for sample_name, sqanti_file in sqanti_files.items():
        with smart_open(sqanti_file, mode="rt") as f:
            FSM_counts = Counter()
            FSM_multi_exon_counts = Counter()
            FSM_single_exon_counts = Counter()
            ISM_counts = Counter()
            ISM_multi_exon_counts = Counter()
            ISM_single_exon_counts = Counter()
            NIC_counts = Counter()
            NNIC_counts = Counter()
            antisense_counts = Counter()
            antisense_multi_exon_counts = Counter()
            antisense_single_exon_counts = Counter()
            intergenic_counts = Counter()
            intergenic_multi_exon_counts = Counter()
            intergenic_single_exon_counts = Counter()
            intronic_counts = Counter()
            intronic_multi_exon_counts = Counter()
            intronic_single_exon_counts = Counter()
            genic_counts = Counter()
            genic_multi_exon_counts = Counter()
            genic_single_exon_counts = Counter()

            f.readline()  # skip header

            for line in f:
                (
                    feature_name,
                    sqanti_cat,
                    feature_length,
                    num_exon_segments,
                    structure,
                    matching_isoforms,
                ) = line.split("\t")

                print("\t".join([sample_name, feature_length]), file=read_lens_ofh)

                # since reused for all assignments, cast only once here
                length = int(feature_length)

                if sqanti_cat == "FSM":
                    FSM_multi_exon_counts[length] += 1
                elif sqanti_cat == "se_FSM":
                    FSM_single_exon_counts[length] += 1
                elif sqanti_cat == "ISM":
                    ISM_multi_exon_counts[length] += 1
                elif sqanti_cat == "se_ISM":
                    ISM_single_exon_counts[length] += 1
                elif sqanti_cat == "NIC":
                    NIC_counts[length] += 1
                elif sqanti_cat == "NNIC":
                    NNIC_counts[length] += 1
                elif sqanti_cat == "antisense":
                    antisense_multi_exon_counts[length] += 1
                elif sqanti_cat == "se_antisense":
                    antisense_single_exon_counts[length] += 1
                elif sqanti_cat == "antisense":
                    antisense_multi_exon_counts[length] += 1
                elif sqanti_cat == "se_antisense":
                    antisense_single_exon_counts[length] += 1
                elif sqanti_cat == "genic":
                    genic_multi_exon_counts[length] += 1
                elif sqanti_cat == "se_genic":
                    genic_single_exon_counts[length] += 1
                elif sqanti_cat == "intergenic":
                    intergenic_multi_exon_counts[length] += 1
                elif sqanti_cat == "se_intergenic":
                    intergenic_single_exon_counts[length] += 1
                elif sqanti_cat == "intronic":
                    intronic_multi_exon_counts[length] += 1
                elif sqanti_cat == "se_intronic":
                    intronic_single_exon_counts[length] += 1

            FSM_counts = FSM_multi_exon_counts + FSM_single_exon_counts
            ISM_counts = ISM_multi_exon_counts + ISM_single_exon_counts
            antisense_counts = antisense_multi_exon_counts + antisense_single_exon_counts
            intergenic_counts = intergenic_multi_exon_counts + intergenic_single_exon_counts
            intronic_counts = intronic_multi_exon_counts + intronic_single_exon_counts
            genic_counts = genic_multi_exon_counts + genic_single_exon_counts

            total_counts = (
                FSM_counts
                + ISM_counts
                + NIC_counts
                + NNIC_counts
                + antisense_counts
                + intergenic_counts
                + intronic_counts
                + genic_counts
            )

            results[sample_name] = {
                "total_counts": total_counts,
                "total_sum": total_counts.total(),
                "FSM_counts": FSM_counts,
                "FSM_sum": FSM_counts.total(),
                "FSM_multi_exon_counts": FSM_multi_exon_counts,
                "FSM_single_exon_counts": FSM_single_exon_counts,
                "ISM_counts": ISM_counts,
                "ISM_sum": ISM_counts.total(),
                "ISM_multi_exon_counts": ISM_multi_exon_counts,
                "ISM_single_exon_counts": ISM_single_exon_counts,
                "NIC_counts": NIC_counts,
                "NIC_sum": NIC_counts.total(),
                "NNIC_counts": NNIC_counts,
                "NNIC_sum": NNIC_counts.total(),
                "antisense_counts": antisense_counts,
                "antisense_sum": antisense_counts.total(),
                "antisense_multi_exon_counts": antisense_multi_exon_counts,
                "antisense_single_exon_counts": antisense_single_exon_counts,
                "intergenic_counts": intergenic_counts,
                "intergenic_sum": intergenic_counts.total(),
                "intergenic_multi_exon_counts": intergenic_multi_exon_counts,
                "intergenic_single_exon_counts": intergenic_single_exon_counts,
                "intergenic_counts": intergenic_counts,
                "intergenic_sum": intergenic_counts.total(),
                "intronic_counts": intronic_counts,
                "intronic_sum": intronic_counts.total(),
                "intronic_multi_exon_counts": intronic_multi_exon_counts,
                "intronic_single_exon_counts": intronic_single_exon_counts,
                "genic_counts": genic_counts,
                "genic_sum": genic_counts.total(),
                "genic_multi_exon_counts": genic_multi_exon_counts,
                "genic_single_exon_counts": genic_single_exon_counts,
            }

    read_lens_ofh.close()

    return results


def get_quantiles(results):
    quantile_one_percent = math.inf
    quantile_ninetynine_percent = -math.inf

    for sample_name, counts_dict in results.items():
        n_reads = counts_dict["total_counts"].total()
        one_percent_reads = math.floor(0.01 * n_reads)

        sorted_keys = sorted(counts_dict["total_counts"])
        current_count = 0
        for k in sorted_keys:
            if current_count + counts_dict["total_counts"][k] >= one_percent_reads:
                quantile_one_percent = (
                    k if k < quantile_one_percent else quantile_one_percent
                )
                break
            current_count += counts_dict["total_counts"][k]

        current_count = 0
        for k in sorted_keys[::-1]:
            if current_count + counts_dict["total_counts"][k] >= one_percent_reads:
                quantile_ninetynine_percent = (
                    k
                    if k > quantile_ninetynine_percent
                    else quantile_ninetynine_percent
                )
                break
            current_count += counts_dict["total_counts"][k]

    return (quantile_one_percent, quantile_ninetynine_percent)


def make_read_length_table(results, category):
    all_dfs = []
    for sample_name, counts_dict in results.items():
        df = pd.DataFrame.from_dict(
            counts_dict[category + "_counts"], orient="index"
        ).reset_index(names="read_length")
        df["sample_name"] = sample_name
        all_dfs.append(df)
    full_df = pd.concat(all_dfs)
    return full_df


# unused at the moment
def make_normalized_read_length_table(results, category):
    all_dfs = []
    for sample_name, counts_dict in results.items():
        df = pd.DataFrame.from_dict(
            counts_dict[category + "_counts"], orient="index"
        ).reset_index(names="read_length")
        df = pd.DataFrame.from_dict(
            {
                k: v / counts_dict[category + "_sum"]
                for (k, v) in counts_dict[category + "_counts"].items()
            },
            orient="index",
        ).reset_index(names="read_length")
        df["sample_name"] = sample_name
        all_dfs.append(df)
    full_df = pd.concat(all_dfs)
    return full_df


def make_total_normalized_read_length_table(results, category):
    all_dfs = []
    for sample_name, counts_dict in results.items():
        df = pd.DataFrame.from_dict(
            counts_dict[category + "_counts"], orient="index"
        ).reset_index(names="read_length")
        df = pd.DataFrame.from_dict(
            {
                k: v / counts_dict["total_sum"]
                for (k, v) in counts_dict[category + "_counts"].items()
            },
            orient="index",
        ).reset_index(names="read_length")
        df["sample_name"] = sample_name
        all_dfs.append(df)
    full_df = pd.concat(all_dfs)
    return full_df


def plot_categories_histogram(
    results,
    output_handle=None,
    categories_to_plot=[
        "FSM",
        "ISM",
        "NIC",
        "NNC",
        "genic_genomic",
        "antisense",
        "fusion",
        "intergenic",
        "genic_intron",
    ],
    categories_spacing=2,
    col_width=1,
):
    n_samples = len(results)
    n_categories = len(categories_to_plot)
    categories_group_spacing = (n_samples * col_width) + categories_spacing

    col_height = np.empty(n_categories * n_samples, dtype=np.float32)
    col_placement = np.empty(n_categories * n_samples, dtype=np.int32)
    col_labels = []
    palette = sns.color_palette(cc.glasbey, n_colors=n_samples)
    col_colors = palette * n_categories
    # col_colors = []

    for i, (sample_name, counts_dict) in enumerate(results.items()):
        col_labels.append(sample_name)
        # col_colors.append(custom_palette[sample_name])
        for j, category in enumerate(categories_to_plot):
            col_height[(n_samples * j) + i] = (
                counts_dict[category + "_sum"] / counts_dict["total_sum"] * 100
            )
            col_placement[(n_samples * j) + i] = (i * col_width) + (
                j * categories_group_spacing
            )

    fig, ax = plt.subplots(figsize=(24, 12))

    ax.bar(
        col_placement,
        col_height,
        width=col_width,
        facecolor=col_colors,
        edgecolor="white",
        linewidth=(0.7 * col_width),
    )
    ax.set_xlim(-1, (categories_group_spacing * n_categories) - categories_spacing)
    ax.set_xticks(
        np.arange(
            ((((n_samples - 1) * col_width)) / 2),
            (categories_group_spacing * n_categories),
            categories_group_spacing,
        )
    )
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 10))
    ax.set_yticks(np.arange(0, 101, 1), minor=True)
    ax.set_yticklabels(np.arange(0, 101, 10), fontsize=20)
    plt.tick_params(axis="y", which="both", color="black", left=True)
    ax.set_xticklabels(categories_to_plot, fontsize=20)
    ax.set_title(
        "Reads Distribution Across Structural Categories and Samples", fontsize=40
    )
    ax.legend(
        handles=[
            matplotlib.patches.Rectangle([0, 0], 2, 1, facecolor=i) for i in palette
        ],
        labels=col_labels,
        facecolor="white",
    )
    # ax.legend(handles=[matplotlib.patches.Rectangle([0,0], 2, 1, facecolor=i) for i in col_colors], labels=col_labels, facecolor="white")

    if output_handle:
        output_handle.savefig(fig)
        plt.close(fig)
    else:
        plt.show()


def plot_subcategories_histogram(
    results,
    category_to_plot,
    output_handle=None,
    subcategories_spacing=2,
    col_width=1,
):
    n_samples = len(results)
    subcategories_group_spacing = (n_samples * col_width) + subcategories_spacing
    subcategories_to_plot = []
    first_sample = next(iter(results))
    for subcategory in results[first_sample].keys():
        if subcategory == category_to_plot + "_counts":
            continue
        elif subcategory[-3:] == "sum":
            continue
        elif category_to_plot in subcategory:
            subcategories_to_plot.append(subcategory)
    n_categories = len(subcategories_to_plot)

    col_height = np.empty(n_categories * n_samples, dtype=np.float32)
    col_placement = np.empty(n_categories * n_samples, dtype=np.int32)
    col_labels = []
    palette = sns.color_palette(cc.glasbey, n_colors=n_samples)
    col_colors = palette * n_categories
    # col_colors = []

    for i, (sample_name, counts_dict) in enumerate(results.items()):
        col_labels.append(sample_name)
        # col_colors.append(custom_palette[sample_name])
        for j, category in enumerate(subcategories_to_plot):
            col_height[(n_samples * j) + i] = (
                counts_dict[category].total()
                / counts_dict[category_to_plot + "_sum"]
                * 100
            )
            col_placement[(n_samples * j) + i] = (i * col_width) + (
                j * subcategories_group_spacing
            )

    fig, ax = plt.subplots(figsize=(24, 12))

    ax.bar(
        col_placement,
        col_height,
        width=col_width,
        facecolor=col_colors,
        edgecolor="white",
        linewidth=(0.7 * col_width),
    )
    ax.set_xlim(
        -1, (subcategories_group_spacing * n_categories) - subcategories_spacing
    )
    ax.set_xticks(
        np.arange(
            ((((n_samples - 1) * col_width)) / 2),
            (subcategories_group_spacing * n_categories),
            subcategories_group_spacing,
        )
    )
    ax.set_ylim(0, 100)
    ax.set_yticks(np.arange(0, 101, 10))
    ax.set_yticks(np.arange(0, 101, 1), minor=True)
    ax.set_yticklabels(np.arange(0, 101, 10), fontsize=20)
    plt.tick_params(axis="y", which="both", color="black", left=True)
    ax.set_xticklabels(subcategories_to_plot, fontsize=16)
    ax.set_title(
        "Reads Distribution Across " + category_to_plot + " Subcategories and Samples",
        fontsize=40,
    )
    ax.legend(
        handles=[
            matplotlib.patches.Rectangle([0, 0], 2, 1, facecolor=i) for i in palette
        ],
        labels=col_labels,
        facecolor="white",
    )
    # ax.legend(handles=[matplotlib.patches.Rectangle([0,0], 2, 1, facecolor=i) for i in col_colors], labels=col_labels, facecolor="white")

    if output_handle:
        output_handle.savefig(fig)
        plt.close(fig)
    else:
        plt.show()


def plot_read_length_distribution_ridgeline(
    df,
    min_x,
    max_x,
    output_handle=None,
    share_y_axis=False,
    title="",
    global_y_max=None,
):

    matplotlib.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=cc.glasbey)

    n_samples = len(df["sample_name"].unique())

    base_height = 5
    target_aspect = 4

    scaled_height = (base_height + (0.3 * base_height * (n_samples - 1))) / n_samples
    scaled_aspect = (base_height * target_aspect) / scaled_height

    scaled_fontsize = 20

    g = sns.FacetGrid(
        df,
        row="sample_name",
        hue="sample_name",
        aspect=scaled_aspect,
        height=scaled_height,  # palette=custom_palette,
        legend_out=False,
        margin_titles=False,
        sharey=share_y_axis,
    )
    g.fig.patch.set_facecolor("none")
    for ax in g.axes.flat:
        ax.patch.set_facecolor("none")
        ax.patch.set_alpha(0)
    g.set_titles("")
    g.fig.suptitle(title, fontsize=scaled_fontsize)

    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    def weighted_hist(x, weights, **kwargs):
        sns.histplot(x=x, weights=weights, **kwargs)

    g.map(
        weighted_hist,
        "read_length",
        0,
        binwidth=20,
        binrange=(min_x, max_x),
        common_bins=True,
        log_scale=False,
        kde=True,
        kde_kws={"clip": (min_x, max_x), "bw_adjust": 0.2},
        line_kws={"alpha": 0.5},
    )

    for key in df.sample_name.unique():
        g.axes_dict[key].lines[0].set_color("black")

    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(
            1.02,
            0.01,
            label,
            fontsize=scaled_fontsize * 0.8,
            fontweight="bold",
            color=color,
            ha="left",
            va="center",
            transform=ax.transAxes,
        )

    g.map(label, 0)

    g.set(yticks=[], ylabel="")
    g.set_xlabels("Read length", fontsize=scaled_fontsize)
    g.tick_params(labelsize=scaled_fontsize * 0.8)
    for ax in g.axes.flat:
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_visible(False)

    if global_y_max is not None:
        for ax in g.axes.flat:
            ax.set_ylim(0, global_y_max)

    g.figure.subplots_adjust(hspace=-0.7)
    if output_handle:
        output_handle.savefig(bbox_inches="tight")
        plt.close(g.fig)
    else:
        plt.show()


def plot_read_length_distribution_ridgeline_manual_kde(
    df,
    min_x,
    max_x,
    output_handle=None,
    share_y_axis=False,
    title="",
    global_y_max=None,
):

    matplotlib.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=cc.glasbey)

    n_samples = len(df["sample_name"].unique())

    base_height = 5
    target_aspect = 4

    scaled_height = (base_height + (0.3 * base_height * (n_samples - 1))) / n_samples
    scaled_aspect = (base_height * target_aspect) / scaled_height

    scaled_fontsize = 20

    g = sns.FacetGrid(
        df,
        row="sample_name",
        hue="sample_name",
        aspect=scaled_aspect,
        height=scaled_height,  # palette=custom_palette,
        legend_out=False,
        margin_titles=False,
        sharey=share_y_axis,
    )
    g.fig.patch.set_facecolor("none")
    for ax in g.axes.flat:
        ax.patch.set_facecolor("none")
        ax.patch.set_alpha(0)
    g.set_titles("")
    g.fig.suptitle(title, fontsize=scaled_fontsize)

    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    def weighted_hist(x, weights, **kwargs):
        sns.histplot(x=x, weights=weights, **kwargs)

    g.map(
        weighted_hist,
        "read_length",
        0,
        binwidth=20,
        binrange=(min_x, max_x),
        alpha=0.5,
        common_bins=True,
        log_scale=False,
        kde=False,
    )

    # Manually computed and scaled KDE
    def manual_scaled_kde(
        data, color, min_x, max_x, binwidth, bw_adjust, grid_points, **kwargs
    ):
        # Extract x values and associated weights
        x_vals = data["read_length"].values
        weights = data[0].values  # Assumes the counts column is named 0
        total_count = weights.sum()
        # Create a grid over the x range
        x_grid = np.linspace(min_x, max_x, grid_points)
        # Compute the KDE using the data and weights
        kde = gaussian_kde(x_vals, weights=weights)
        kde.set_bandwidth(bw_adjust * kde.factor)
        density = kde(x_grid)
        # Scale the density so its area corresponds to total counts.
        # Multiply by (total_count * binwidth) so that the density is on the same scale as the histogram.
        scaled_density = density * total_count * binwidth
        plt.plot(x_grid, scaled_density, color=color, linewidth=2)

    # Overlay the manual KDE using map_dataframe
    gridpoints = int((max_x - min_x) / 2)  # /20 the bin width and *10
    g.map_dataframe(
        manual_scaled_kde,
        min_x=min_x,
        max_x=max_x,
        binwidth=20,
        bw_adjust=0.2,
        grid_points=gridpoints,
    )

    for key in df.sample_name.unique():
        g.axes_dict[key].lines[0].set_color("black")

    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    def label(x, color, label):
        ax = plt.gca()
        ax.text(
            1.02,
            0.01,
            label,
            fontsize=scaled_fontsize * 0.8,
            fontweight="bold",
            color=color,
            ha="left",
            va="center",
            transform=ax.transAxes,
        )

    g.map(label, 0)

    g.set(yticks=[], ylabel="")
    g.set_xlabels("Read length", fontsize=scaled_fontsize)
    g.tick_params(labelsize=scaled_fontsize * 0.8)
    for ax in g.axes.flat:
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_visible(False)

    if global_y_max is not None:
        for ax in g.axes.flat:
            ax.set_ylim(0, global_y_max)

    g.figure.subplots_adjust(hspace=-0.7)
    if output_handle:
        output_handle.savefig(bbox_inches="tight")
        plt.close(g.fig)
    else:
        plt.show()


def plot_trend_line_only(
    df,
    min_x,
    max_x,
    output_handle=None,
    share_y_axis=True,
    title="",
    global_y_max=None,
    ylabel="",
):

    matplotlib.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=cc.glasbey)

    n_samples = len(df["sample_name"].unique())

    g = sns.FacetGrid(
        df,
        row="sample_name",
        hue="sample_name",
        aspect=(2 * n_samples),
        height=(10 / n_samples),  # palette=custom_palette,
        legend_out=False,
        margin_titles=False,
        sharey=share_y_axis,
    )
    g.fig.patch.set_facecolor("none")
    for ax in g.axes.flat:
        ax.patch.set_facecolor("none")
        ax.patch.set_alpha(0)
    g.set_titles("")
    g.fig.suptitle(title, fontsize=20)

    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    def weighted_hist(x, weights, **kwargs):
        sns.histplot(x=x, weights=weights, **kwargs)

    g.map(
        weighted_hist,
        "read_length",
        0,
        binwidth=20,
        binrange=(min_x, max_x),
        common_bins=True,
        log_scale=False,
        kde=True,
        fill=False,
        linewidth=0,
        kde_kws={"clip": (min_x, max_x), "bw_adjust": 0.2},
        line_kws={"alpha": 0.5, "linewidth": 5},
    )

    g.refline(y=0, linewidth=2, linestyle="-", color="black", clip_on=False)

    g.set_ylabels(ylabel, fontsize=20)
    g.set_xlabels("Read length", fontsize=20)
    g.tick_params(labelsize=18)
    for ax in g.axes.flat:
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_visible(False)

    if global_y_max is not None:
        for ax in g.axes.flat:
            ax.set_ylim(0, global_y_max)

    # palette = sns.color_palette(cc.glasbey, n_colors=len(df))
    g.add_legend(
        handles=[
            matplotlib.patches.Rectangle([0, 0], 2, 1, facecolor=i) for i in g._colors
        ],
        labels=g.hue_names,
        fontsize=16,
        facecolor="white",
    )

    g.figure.subplots_adjust(hspace=-1)
    if output_handle:
        output_handle.savefig(bbox_inches="tight")
        plt.close(g.fig)
    else:
        plt.show()


def plot_trend_line_manual_kde(
    df,
    min_x,
    max_x,
    output_handle=None,
    share_y_axis=True,
    title="",
    global_y_max=None,
    ylabel="",
):
    matplotlib.rcParams["axes.prop_cycle"] = matplotlib.cycler(color=cc.glasbey)

    n_samples = len(df["sample_name"].unique())

    g = sns.FacetGrid(
        df,
        row="sample_name",
        hue="sample_name",
        aspect=(2 * n_samples),
        height=(10 / n_samples),  # palette=custom_palette,
        legend_out=False,
        margin_titles=False,
        sharey=share_y_axis,
    )
    g.fig.patch.set_facecolor("none")
    for ax in g.axes.flat:
        ax.patch.set_facecolor("none")
        ax.patch.set_alpha(0)
    g.set_titles("")
    g.fig.suptitle(title, fontsize=20)
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    # First, plot the histogram WITHOUT the built-in KDE.
    def weighted_hist(x, weights, **kwargs):
        sns.histplot(
            x=x,
            weights=weights,
            binwidth=20,
            binrange=(min_x, max_x),
            common_bins=True,
            log_scale=False,
            kde=False,
            fill=False,
            linewidth=0,
        )

    g.map(weighted_hist, "read_length", 0)

    # # Manually computed and scaled KDE
    def manual_scaled_kde(
        data, color, min_x, max_x, binwidth, bw_adjust, grid_points, **kwargs
    ):
        # Extract x values and associated weights
        x_vals = data["read_length"].values
        weights = data[0].values  # Assumes the counts column is named 0
        total_count = weights.sum()
        # Create a grid over the x range
        x_grid = np.linspace(min_x, max_x, grid_points)
        # Compute the KDE using the data and weights
        kde = gaussian_kde(x_vals, weights=weights)
        kde.set_bandwidth(bw_adjust * kde.factor)
        density = kde(x_grid)
        # Scale the density so its area corresponds to total counts.
        # Multiply by (total_count * binwidth) so that the density is on the same scale as the histogram.
        scaled_density = density * total_count * binwidth
        plt.plot(x_grid, scaled_density, color=color, linewidth=2)

    # Overlay the manual KDE using map_dataframe.
    gridpoints = int((max_x - min_x) / 2)  # /20 the bin width and *10
    g.map_dataframe(
        manual_scaled_kde,
        min_x=min_x,
        max_x=max_x,
        binwidth=20,
        bw_adjust=0.2,
        grid_points=gridpoints,
    )

    # Draw a reference line at y=0.
    g.refline(y=0, linewidth=2, linestyle="-", color=None, clip_on=False)

    g.set_ylabels(ylabel, fontsize=20)
    g.set_xlabels("Read length", fontsize=20)
    g.tick_params(labelsize=18)
    for ax in g.axes.flat:
        for spine in ["left", "bottom"]:
            ax.spines[spine].set_visible(False)

    if global_y_max is not None:
        for ax in g.axes.flat:
            ax.set_ylim(0, global_y_max)

    # Create a legend using your custom colors.
    # palette = list(custom_palette.values())
    g.add_legend(
        handles=[
            matplotlib.patches.Rectangle([0, 0], 2, 1, facecolor=i) for i in g._colors
        ],
        labels=g.hue_names,
        fontsize=16,
        facecolor="white",
    )

    g.figure.subplots_adjust(hspace=-1)

    if output_handle:
        # Use bbox_inches="tight" to avoid cutting off labels.
        output_handle.savefig(bbox_inches="tight")
        plt.close(g.fig)
    else:
        plt.show()


def report_sqanti_category_summary_table(
    results, distribution_categories, output_tsv_filename
):

    full_df = None
    for sample_name, sample_dict in results.items():

        total_sum = sample_dict["total_sum"]

        # restrict to subset of cateogries of interest.
        sample_dict = {
            k: sample_dict[k + "_sum"] / total_sum * 100
            for k in distribution_categories
            if (k + "_sum") in sample_dict
        }

        sample_df = pd.DataFrame(sample_dict, index=[0])
        sample_df["Sample"] = sample_name

        if full_df is None:
            full_df = sample_df
        else:
            full_df = pd.concat([full_df, sample_df])

    full_df.to_csv(output_tsv_filename, sep="\t", index=False)

    return


def parse_args():
    parser = argparse.ArgumentParser(
        description="Associate classification files with sample names."
    )
    parser.add_argument(
        "--classification_files",
        required=True,
        help="Comma-separated list of classification file paths.",
    )
    parser.add_argument(
        "--sample_names",
        required=True,
        help="Comma-separated list of sample names (in same order as classification_files).",
    )
    parser.add_argument(
        "--output", required=True, help="Output file or output string PREFIX."
    )
    parser.add_argument(
        "--type", required=True, help="Type of classification or output (as a string)."
    )
    return parser.parse_args()


def main():
    args = parse_args()

    classification_files = args.classification_files.split(",")
    sample_names = args.sample_names.split(",")

    if len(classification_files) != len(sample_names):
        raise ValueError(
            "The number of classification files must match the number of sample names."
        )

    sample_dict = {
        sample: file for sample, file in zip(sample_names, classification_files)
    }

    results = None
    read_length_categories = None
    distribution_categories = None
    with_subcategories = None
    if args.type == "sqanti":
        # include read length tsv file
        read_lengths_outfile = args.output + ".read_lengths.tsv"
        results = parse_sqanti_classification(sample_dict, read_lengths_outfile)

        read_length_categories = [
            "total",
            "FSM",
            "ISM",
            "NIC",
            "NNC",
            "antisense",
            "fusion",
            "genic_genomic",
            "genic_intron",
            "intergenic",
        ]
        distribution_categories = [
            "FSM",
            "ISM",
            "NIC",
            "NNC",
            "antisense",
            "fusion",
            "genic_genomic",
            "genic_intron",
            "intergenic",
        ]
        with_subcategories = ["FSM", "ISM", "NIC", "NNC", "fusion"]

        report_sqanti_category_summary_table(
            results, distribution_categories, args.output + ".categories.tsv"
        )
    elif args.type == "lraa_sqanti_like":
        read_lengths_outfile = args.output + ".read_lengths.tsv"
        results = parse_lraa_sqantilike_classification(sample_dict, read_lengths_outfile)

        read_length_categories = [
            "total",
            "FSM",
            "ISM",
            "NIC",
            "NNIC",
            "antisense",
            "genic",
            "intronic",
            "intergenic",
        ]
        distribution_categories = [
            "FSM",
            "ISM",
            "NIC",
            "NNIC",
            "antisense",
            "genic",
            "intronic",
            "intergenic",
        ]
        with_subcategories = ["FSM", "ISM", "antisense", "genic", "intronic", "intergenic"]

        report_sqanti_category_summary_table(
            results, distribution_categories, output_name + ".categories.tsv"
        )   

    (quantile_one_percent, quantile_ninetynine_percent) = get_quantiles(results)

    # histograms of categories split and of subcategories within each category
    pdf_handle = PdfPages(args.output + ".categories.pdf")
    plot_categories_histogram(
        results, output_handle=pdf_handle, categories_to_plot=distribution_categories
    )
    for category in with_subcategories:
        plot_subcategories_histogram(
            results, output_handle=pdf_handle, category_to_plot=category
        )
    pdf_handle.close()

    pdf_handle = PdfPages(args.output + ".read_length_distributions.pdf")
    fig, ax = plt.subplots()
    plt.close()
    for category in read_length_categories:
        full_df = make_read_length_table(results, category)
        # plot_read_length_distribution_ridgeline(
        #      full_df,
        #      min_x = quantile_one_percent,
        #      max_x = quantile_ninetynine_percent,
        #      output_handle = pdf_handle,
        #      share_y_axis = True,
        #      title = "Read length distribution of " + category + " reads in absolute counts",
        # )

        # plot_trend_line_only(
        #      full_df,
        #      min_x = quantile_one_percent,
        #      max_x = quantile_ninetynine_percent,
        #      output_handle = pdf_handle,
        #      share_y_axis = True,
        #      title = "Read length distribution of " + category + " reads in absolute counts",
        #      ylabel = "Absolute read counts"
        # )

        plot_read_length_distribution_ridgeline_manual_kde(
            full_df,
            min_x=quantile_one_percent,
            max_x=quantile_ninetynine_percent,
            output_handle=pdf_handle,
            share_y_axis=True,
            title="Read length distribution of "
            + category
            + " reads in absolute counts",
        )

        plot_trend_line_manual_kde(
            full_df,
            min_x=quantile_one_percent,
            max_x=quantile_ninetynine_percent,
            output_handle=pdf_handle,
            share_y_axis=True,
            title="Read length distribution of "
            + category
            + " reads in absolute counts",
            ylabel="Absolute read counts",
        )

        # total normalized distributions
        full_df = make_total_normalized_read_length_table(results, category)
        # plot_read_length_distribution_ridgeline(
        #      full_df,
        #      min_x = quantile_one_percent,
        #      max_x = quantile_ninetynine_percent,
        #      output_handle = pdf_handle,
        #      share_y_axis = True,
        #      title = "Normalized read length distribution of " + category + " reads",
        # )

        # plot_trend_line_only(
        #      full_df,
        #      min_x = quantile_one_percent,
        #      max_x = quantile_ninetynine_percent,
        #      output_handle = pdf_handle,
        #      share_y_axis = True,
        #      title = "Normalized read length distribution of " + category + " reads",
        #      ylabel = "Normalized read count"
        # )

        plot_read_length_distribution_ridgeline_manual_kde(
            full_df,
            min_x=quantile_one_percent,
            max_x=quantile_ninetynine_percent,
            output_handle=pdf_handle,
            share_y_axis=True,
            title="Normalized read length distribution of " + category + " reads",
        )

        plot_trend_line_manual_kde(
            full_df,
            min_x=quantile_one_percent,
            max_x=quantile_ninetynine_percent,
            output_handle=pdf_handle,
            share_y_axis=True,
            title="Normalized read length distribution of " + category + " reads",
            ylabel="Normalized read count",
        )
    pdf_handle.close()


if __name__ == "__main__":
    main()
