"""Calculate the distance between sequences between seasons.
"""
import argparse

from augur.distance import read_distance_map, get_distance_between_nodes
from augur.frequency_estimators import timestamp_to_float
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


def get_distances_to_last_ancestor_by_tips(tips, sequences_by_node_and_gene, distance_map, latest_date):
    """Calculate distances between each sample in the given sequences and its last
    ancestor in a previous season using the given distance map.

    Parameters
    ----------
    tips : list
        a list of Bio.Phylo nodes for tips to find the ancestral distance to

    sequences_by_node_and_gene : dict
        nucleotide or amino acid sequences by node name and gene

    distance_map : dict
        site-specific and, optionally, sequence-specific distances between two
        sequences

    latest_date : pandas.Timestamp
        latest date to consider a node as a potential ancestor of a given
        sample; used to define a previous season relative to the most recent
        samples in the given tree.

    Returns
    -------
    dict :
        distances calculated between each sample in the tree and its last
        ancestor sequence with distances indexed by node name

    """
    if latest_date is not None:
        latest_date = timestamp_to_float(latest_date)

    distances_by_node = {}

    # Calculate distance between each tip and its closest ancestor in the last
    # season as defined by the given latest date threshold.
    for node in tips:
        # If the given latest date is not None, skip nodes that were sampled
        # prior to this date.
        if latest_date is not None and node.attr["num_date"] < latest_date:
            continue

        # Find the closest ancestor of this node that was also sampled prior to
        # the given latest date. Stop searching once we reach the root. If the
        # latest date requested is None, the immediate parent of each node will
        # be used.
        parent = node.parent
        while parent != tree.root and latest_date is not None and parent.attr["num_date"] > latest_date:
            parent = parent.parent

        # Calculate distance between current node and its ancestor.
        distances_by_node[node.name] = get_distance_between_nodes(
            sequences_by_node_and_gene[parent.name],
            sequences_by_node_and_gene[node.name],
            distance_map
        )

    return distances_by_node


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
    parser.add_argument("--months-prior-to-season", type=int, help="number of months prior to each season to look for a tip's ancestor", required=True)
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)

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
    seasons = pd.date_range(args.start_date, args.end_date, freq="%iMS" % args.interval)

    # Load date annotations and annotate tree with them.
    date_annotations = read_node_data(args.date_annotations)
    for node in tree.find_clades():
        node.attr = date_annotations["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Assign tips to seasons.
    tips_assigned = set()
    tips_by_season = {}
    for season in seasons:
        season_date = season.strftime("%Y-%m-%d")
        season_float = timestamp_to_float(season)
        tips_by_season[season_date] = []

        for tip in tree.find_clades(terminal=True):
            if tip.name not in tips_assigned and tip.attr["num_date"] < season_float:
                tips_assigned.add(tip.name)
                tips_by_season[season_date].append(tip)

    final_distances_by_node = {}
    distance_map_names = []
    compare_to_list = ["ancestor"] * len(args.attribute_name)
    for compare_to, attribute, distance_map_file in zip(compare_to_list, args.attribute_name, args.map):
        # Load the given distance map.
        distance_map = read_distance_map(distance_map_file)
        distance_map_names.append(distance_map.get("name", distance_map_file))

        for season_date, tips in tips_by_season.items():
            # Use the distance map to calculate distances between all samples in the
            # given tree and their last ancestor in the previous season.
            latest_date = pd.to_datetime(season_date) - pd.DateOffset(months=args.months_prior_to_season)

            distances_by_node = get_distances_to_last_ancestor_by_tips(
                tips,
                sequences_by_node_and_gene,
                distance_map,
                latest_date
            )

            # Map distances to the requested attribute name.
            # Convert data like:
            # {
            #   "A/AbuDhabi/24/2017": 1
            # }
            # to data like:
            #
            # {
            #   "A/AbuDhabi/24/2017": {
            #     "ep": 1
            #   }
            # }
            #
            for node_name, values in distances_by_node.items():
                if node_name not in final_distances_by_node:
                    final_distances_by_node[node_name] = {}

                final_distances_by_node[node_name][attribute] = values

    # Prepare params for export.
    params = {
        "attribute": args.attribute_name,
        "map_name": distance_map_names,
        "start_date": args.start_date,
        "end_date": args.end_date,
        "interval": args.interval,
        "months_prior_to_season": args.months_prior_to_season
    }

    # Export distances to JSON.
    write_json({"params": params, "nodes": final_distances_by_node}, args.output)
