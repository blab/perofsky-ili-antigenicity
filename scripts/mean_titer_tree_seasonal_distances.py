"""Calculate the distance between sequences between seasons.
"""
import argparse

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


def annotate_terminals_for_tree(tree):
    """Make a single pass through the tree in postorder to store a set of all
    terminals descending from each node. This uses more memory, but it allows faster
    identification of MRCAs between any pair of tips in the tree and speeds up
    pairwise distance calculations by orders of magnitude.

    """
    for node in tree.find_clades(order="postorder"):
        node.terminals = set()
        for child in node.clades:
            if child.is_terminal():
                node.terminals.add(child.name)
            else:
                node.terminals.update(child.terminals)

    return tree


def get_titer_distance_between_nodes(tree, past_node, current_node, titer_attr="dTiter"):
    # Find MRCA of tips from one tip up. Sum the titer attribute of interest
    # while walking up to the MRCA, to avoid an additional pass later. The loop
    # below stops when the past node is found in the list of the candidate
    # MRCA's terminals. This test should always evaluate to true when the MRCA
    # is the root node, so we should not have to worry about trying to find the
    # parent of the root.
    current_node_branch_sum = 0.0
    mrca = current_node
    while past_node.name not in mrca.terminals:
        current_node_branch_sum += mrca.attr[titer_attr]
        mrca = mrca.parent

    # Sum the node weights for the other tip from the bottom up until we reach
    # the MRCA. The value of the MRCA is intentionally excluded here, as it
    # would represent the branch leading to the MRCA and would be outside the
    # path between the two tips.
    past_node_branch_sum = 0.0
    current_node = past_node
    while current_node != mrca:
        past_node_branch_sum += current_node.attr[titer_attr]
        current_node = current_node.parent

    final_sum = past_node_branch_sum + current_node_branch_sum
    return final_sum


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree", required=True)
    parser.add_argument("--model", help="JSON of titer tree model with weights assigned to nodes of the given tree", required=True)
    parser.add_argument("--model-attribute-name", help="name of attribute to use from titer model file", default="dTiter")
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
    tree = annotate_terminals_for_tree(tree)

    # Create season intervals.
    seasons = get_seasons(args.start_date, args.end_date, args.interval)

    # Load node annotations and annotate tree with them.
    node_annotations = read_node_data([args.date_annotations, args.model])
    for node in tree.find_clades():
        node.attr = node_annotations["nodes"][node.name]
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
    for current_season_index, (current_season_start_date, current_season_end_date) in enumerate(season_keys):
        other_season_start_index = current_season_index
        other_season_end_index = current_season_index + args.seasons_away_to_compare + 1

        for other_season_start_date, other_season_end_date in season_keys[other_season_start_index:other_season_end_index]:
            print("Compare %s to %s" % (current_season_start_date, other_season_start_date), flush=True)
            total_distance = 0
            total_comparisons = 0
            total_other_season_tips = float(len(tips_by_season[(other_season_start_date, other_season_end_date)]))

            for current_season_tip in tips_by_season[(current_season_start_date, current_season_end_date)]:
                current_tip_distance = 0

                for other_season_tip in tips_by_season[(other_season_start_date, other_season_end_date)]:
                    current_tip_distance += get_titer_distance_between_nodes(
                        tree,
                        other_season_tip,
                        current_season_tip,
                        args.model_attribute_name,
                    )
                    total_comparisons += 1

                total_distance += current_tip_distance
                if args.strain_output:
                    pairwise_seasonal_distances.append((
                        current_season_start_date,
                        current_season_end_date,
                        other_season_start_date,
                        other_season_end_date,
                        "titer_tree_model",
                        current_season_tip.name,
                        current_tip_distance / total_other_season_tips
                    ))

            mean_seasonal_distances.append([
                current_season_start_date,
                current_season_end_date,
                other_season_start_date,
                other_season_end_date,
                "titer_tree_model",
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
