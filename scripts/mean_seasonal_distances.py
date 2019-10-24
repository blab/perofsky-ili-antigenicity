"""Calculate the distance between sequences between seasons.
"""
import argparse

from augur.distance import read_distance_map, get_distance_between_nodes
from augur.reconstruct_sequences import load_alignments
from augur.utils import annotate_parents_for_tree, read_node_data, write_json

import Bio
import Bio.Phylo
from collections import defaultdict
import copy
import json
import numpy as np
import pandas as pd
import pprint
import sys

from utils import get_seasons, get_nodes_by_season


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--alignment", nargs="+", help="sequence(s) to be used, supplied as FASTA files", required=True)
    parser.add_argument('--gene-names', nargs="+", type=str, help="names of the sequences in the alignment, same order assumed", required=True)
    parser.add_argument("--attribute-name", nargs="+", help="name to store distances associated with the given distance map; multiple attribute names are linked to corresponding positional comparison method and distance map arguments", required=True)
    parser.add_argument("--map", nargs="+", help="JSON providing the distance map between sites and, optionally, sequences present at those sites; the distance map JSON minimally requires a 'default' field defining a default numeric distance and a 'map' field defining a dictionary of genes and one-based coordinates", required=True)
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date", required=True)
    parser.add_argument("--start-date", help="date to start seasonal intervals (e.g., 2000-10-01)", required=True)
    parser.add_argument("--end-date", help="date to end seasonal intervals (e.g., 2010-10-01)", required=True)
    parser.add_argument("--interval", type=int, help="number of months per season", required=True)
    parser.add_argument("--seasons-away-to-compare", type=int, help="number of seasons away from the current season to compare (e.g., 2 for the current season and the two that follow it)", required=True)
    parser.add_argument("--output", help="tab-delimited file with mean distances by node name and attribute name", required=True)
    parser.add_argument("--strain-output", help="tab-delimited file with mean pairwise distance per strain by node names and attribute name")

    args = parser.parse_args()

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load sequences.
    alignments = load_alignments(args.alignment, args.gene_names)

    # Index sequences by node name and gene.
    sequences_by_node_and_gene = defaultdict(dict)
    for gene, alignment in alignments.items():
        for record in alignment:
            sequences_by_node_and_gene[record.name][gene] = str(record.seq)

    # Create season intervals.
    seasons = get_seasons(args.start_date, args.end_date, args.interval)

    # Load date annotations and annotate tree with them.
    date_annotations = read_node_data(args.date_annotations)
    for node in tree.find_clades():
        node.attr = date_annotations["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Assign tips to seasons.
    tips_by_season = get_nodes_by_season(tree, seasons, terminal_only=True)

    # Build a list of season keys for iteration downstream. These are pairs of
    # dates representing the start and end of each season.
    season_keys = sorted(tips_by_season.keys())

    # Store the mean pairwise distance between seasons for one or more distance
    # maps (e.g., epitope sites, non-epitope sites, etc.).
    mean_seasonal_distances = []
    pairwise_seasonal_distances = []
    for attribute, distance_map_file in zip(args.attribute_name, args.map):
        # Load the given distance map.
        distance_map = read_distance_map(distance_map_file)

        for current_season_index, (current_season_start_date, current_season_end_date) in enumerate(season_keys):
            other_season_start_index = current_season_index
            other_season_end_index = current_season_index + args.seasons_away_to_compare + 1

            for other_season_start_date, other_season_end_date in season_keys[other_season_start_index:other_season_end_index]:
                print("Compare %s to %s with %s map" % (current_season_start_date, other_season_start_date, distance_map["name"]), flush=True)
                total_distance = 0
                total_comparisons = 0
                total_other_season_tips = float(len(tips_by_season[(other_season_start_date, other_season_end_date)]))

                for current_season_tip in tips_by_season[(current_season_start_date, current_season_end_date)]:
                    current_tip_distance = 0

                    for other_season_tip in tips_by_season[(other_season_start_date, other_season_end_date)]:
                        current_tip_distance += get_distance_between_nodes(
                            sequences_by_node_and_gene[current_season_tip.name],
                            sequences_by_node_and_gene[other_season_tip.name],
                            distance_map
                        )
                        total_comparisons += 1

                    total_distance += current_tip_distance
                    if args.strain_output:
                        pairwise_seasonal_distances.append((
                            current_season_start_date,
                            current_season_end_date,
                            other_season_start_date,
                            other_season_end_date,
                            distance_map["name"],
                            current_season_tip.name,
                            current_tip_distance / total_other_season_tips
                        ))

                mean_seasonal_distances.append([
                    current_season_start_date,
                    current_season_end_date,
                    other_season_start_date,
                    other_season_end_date,
                    distance_map["name"],
                    total_distance,
                    total_comparisons,
                    total_distance / float(total_comparisons)
                ])

    # Export distances to a table.
    df = pd.DataFrame(
        mean_seasonal_distances,
        columns=["current_season_start", "current_season_end", "other_season_start", "other_season_end", "distance_map", "total_distance", "total_comparisons", "mean_distance"]
    )
    df.to_csv(args.output, sep="\t", float_format="%.2f", index=False)

    if args.strain_output:
        # Export all pairwise distances to a table.
        pairwise_df = pd.DataFrame(
            pairwise_seasonal_distances,
            columns=["current_season_start", "current_season_end", "other_season_start", "other_season_end", "distance_map", "current_strain", "mean_distance"]
        )
        pairwise_df.to_csv(args.strain_output, sep="\t", float_format="%.2f", index=False)
