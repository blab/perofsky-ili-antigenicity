"""Calculate the distance between sequences between seasons.
"""
import argparse

from augur.distance import read_distance_map, get_distance_between_nodes
from augur.frequency_estimators import timestamp_to_float
from augur.reconstruct_sequences import load_alignments
from augur.utils import annotate_parents_for_tree, read_node_data, write_json

import Bio
import Bio.Phylo
from collections import defaultdict, OrderedDict
import copy
import datetime
import json
import numpy as np
import pandas as pd
import pprint
import sys
from treetime.utils import numeric_date


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


def get_titer_tree_distances_to_ancestor_by_tips(tips, tree, latest_date):
    """Calculate distances between each sample in the given list of tips and its
    last ancestor in a previous season using the titer tree model annotations in
    the given tree.

    Parameters
    ----------
    tips : list
        a list of Bio.Phylo nodes for tips to find the ancestral distance to

    tree : Bio.Phylo
        a tree annotated with titer tree model distances in the `dTiter` node attribute

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
        distances_by_node[node.name] = get_titer_distance_between_nodes(
            tree,
            node,
            parent
        )

    return distances_by_node


def get_vaccine_strain_for_date(vaccines, date):
    """Get the name of the vaccine strain whose start date immediately precedes the given date.

    Parameters
    ----------
    vaccines : pandas.DataFrame
        table of vaccine strains and their numeric vaccine start dates in float format (e.g., 2010.25).

    date : float
        numeric date (e.g., 2010.25).

    Returns
    -------
    str :
        name of the vaccine strain preceding the given date

    """
    # Find all vaccine strains with dates just before the given date.
    potential_strains = vaccines[vaccines["numeric_vaccine_start_date"] < date]["strain"].values

    # Select the last item in the sorted list.
    if len(potential_strains) > 0:
        return potential_strains[-1]
    else:
        return None


def get_distances_to_vaccine_by_tips(tips, sequences_by_node_and_gene, distance_map, vaccines):
    """Calculate distances between each sample in the given sequences and the
    vaccine that was available in the season those sequences were circulating using
    the given distance map.

    Parameters
    ----------
    tips : list
        a list of Bio.Phylo nodes for tips to find the ancestral distance to

    sequences_by_node_and_gene : dict
        nucleotide or amino acid sequences by node name and gene

    distance_map : dict
        site-specific and, optionally, sequence-specific distances between two
        sequences

    vaccines : pandas.DataFrame
        table of vaccine strains and their corresponding vaccine start dates

    Returns
    -------
    dict :
        distances calculated between each sample in the tree and its vaccine
        sequence with distances indexed by node name

    dict :
        vaccines strain used per sample

    """
    distances_by_node = {}
    vaccines_by_node = {}

    # Calculate distance between each tip and its vaccine from the same
    # season as defined by the given latest date threshold.
    for node in tips:
        # Find the appropriate vaccine for the given tip's date.
        vaccine_strain = get_vaccine_strain_for_date(vaccines, node.attr["num_date"])
        if vaccine_strain is None:
            continue

        print(f"Using {vaccine_strain} for {node.name} with date of {node.attr['num_date']}", flush=True)

        distances_by_node[node.name] = get_distance_between_nodes(
            sequences_by_node_and_gene[vaccine_strain],
            sequences_by_node_and_gene[node.name],
            distance_map
        )
        vaccines_by_node[node.name] = vaccine_strain

    return distances_by_node, vaccines_by_node


def get_titer_tree_distances_to_vaccine_by_tips(tips, tree, vaccines):
    """
    tips : list
        a list of Bio.Phylo nodes for tips to find the distance to

    tree : Bio.Phylo
        a tree annotated with titer tree model distances in the `dTiter` node attribute

    vaccines : pandas.DataFrame
        table of vaccine strains and their corresponding vaccine start dates
    """
    distances_by_node = {}
    vaccines_by_node = {}

    # Calculate titer tree distance between each tip and its vaccine from the
    # same season as defined by the given latest date threshold.
    for node in tips:
        # Find the appropriate vaccine for the given tip's date.
        vaccine_strain = get_vaccine_strain_for_date(vaccines, node.attr["num_date"])
        if vaccine_strain is None:
            continue

        print(f"Using {vaccine_strain} for {node.name} with date of {node.attr['num_date']}", flush=True)

        # Get the vaccine strain node.
        vaccine_node = [
            node
            for node in tree.find_clades(terminal=True)
            if node.name == vaccine_strain
        ][0]

        distances_by_node[node.name] = get_titer_distance_between_nodes(
            tree,
            vaccine_node,
            node
        )
        vaccines_by_node[node.name] = vaccine_strain

    return distances_by_node, vaccines_by_node


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
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--alignment", nargs="+", help="sequence(s) to be used, supplied as FASTA files")
    parser.add_argument('--gene-names', nargs="+", type=str, help="names of the sequences in the alignment, same order assumed")
    parser.add_argument("--attribute-name", nargs="+", help="name to store distances associated with the given distance map; multiple attribute names are linked to corresponding positional comparison method and distance map arguments", required=True)
    parser.add_argument("--map", nargs="+", help="JSON providing the distance map between sites and, optionally, sequences present at those sites; the distance map JSON minimally requires a 'default' field defining a default numeric distance and a 'map' field defining a dictionary of genes and one-based coordinates")
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date", required=True)
    parser.add_argument("--start-date", help="date to start seasonal intervals (e.g., 2000-10-01)", required=True)
    parser.add_argument("--end-date", help="date to end seasonal intervals (e.g., 2010-10-01)", required=True)
    parser.add_argument("--interval", type=int, help="number of months per season", required=True)
    parser.add_argument("--months-prior-to-season", type=int, help="number of months prior to each season to look for a tip's ancestor", required=True)
    parser.add_argument("--vaccines", help="table of vaccine strains annotated by selection and start date to use for seasonal distance calculations instead of the ancestor from the last season")
    parser.add_argument("--titer-tree-model", help="JSON representing titer tree model mapping weights to branches in the given tree")
    parser.add_argument("--output", help="JSON file with calculated distances stored by node name and attribute name", required=True)

    args = parser.parse_args()

    if args.alignment is None and args.titer_tree_model is None:
        print("Error: You must define either an alignment or titer tree model to use for distance calculations.", file=sys.stderr)
        sys.exit(1)

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Load sequences.
    if args.alignment and args.gene_names:
        alignments = load_alignments(args.alignment, args.gene_names)

        # Index sequences by node name and gene.
        sequences_by_node_and_gene = defaultdict(dict)
        for gene, alignment in alignments.items():
            for record in alignment:
                sequences_by_node_and_gene[record.name][gene] = str(record.seq)

    # Create season intervals.
    seasons = pd.date_range(args.start_date, args.end_date, freq="%iMS" % args.interval)
    print("Seasons:")
    print(seasons)

    # Load date annotations and annotate tree with them.
    date_annotations = read_node_data(args.date_annotations)
    for node in tree.find_clades():
        node.attr = date_annotations["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Load titer model annotations, if they have been provided.
    if args.titer_tree_model:
        titer_annotations = read_node_data(args.titer_tree_model)
        for node in tree.find_clades():
            node.attr["dTiter"] = titer_annotations["nodes"][node.name]["dTiter"]

        # Make a single pass through the tree in postorder to store a set of all
        # terminals descending from each node. This uses more memory, but it
        # allows faster identification of MRCAs between any pair of tips in the
        # tree and speeds up pairwise distance calculations by orders of
        # magnitude.
        for node in tree.find_clades(order="postorder"):
            node.terminals = set()
            for child in node.clades:
                if child.is_terminal():
                    node.terminals.add(child.name)
                else:
                    node.terminals.update(child.terminals)

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

    # If we're comparing to vaccine strains instead of last ancestors, assign
    # vaccines to seasons by their selection dates.
    if args.vaccines is not None:
        vaccines = pd.read_csv(
            args.vaccines,
            sep="\t",
            parse_dates=["selection_date", "vaccine_start_date"]
        )
        vaccines["numeric_vaccine_start_date"] = vaccines["vaccine_start_date"].apply(numeric_date)

        print(seasons)
        print(vaccines)
        compare_to_list = ["vaccine"] * len(args.attribute_name)
    else:
        compare_to_list = ["ancestor"] * len(args.attribute_name)

    if args.map is not None:
        maps = args.map
    else:
        maps = [None] * len(args.attribute_name)

    final_distances_by_node = {}
    distance_map_names = []
    for compare_to, attribute, distance_map_file in zip(compare_to_list, args.attribute_name, maps):
        if distance_map_file:
            # Load the given distance map.
            distance_map = read_distance_map(distance_map_file)
            distance_map_names.append(distance_map.get("name", distance_map_file))

        for season_date, tips in tips_by_season.items():
            vaccines_by_node = None

            if compare_to == "ancestor":
                # Use the distance map to calculate distances between all samples in the
                # given tree and their last ancestor in the previous season.
                latest_date = pd.to_datetime(season_date) - pd.DateOffset(months=args.months_prior_to_season)

                if args.titer_tree_model:
                    distances_by_node = get_titer_tree_distances_to_ancestor_by_tips(
                        tips,
                        tree,
                        latest_date
                    )
                else:
                    distances_by_node = get_distances_to_last_ancestor_by_tips(
                        tips,
                        sequences_by_node_and_gene,
                        distance_map,
                        latest_date
                    )
            elif compare_to == "vaccine":
                if args.titer_tree_model:
                    distances_by_node, vaccines_by_node = get_titer_tree_distances_to_vaccine_by_tips(
                        tips,
                        tree,
                        vaccines
                    )
                else:
                    distances_by_node, vaccines_by_node = get_distances_to_vaccine_by_tips(
                        tips,
                        sequences_by_node_and_gene,
                        distance_map,
                        vaccines
                    )
            else:
                print(f"Cannot compare tips with the requested method: {compare_to}", file=sys.stderr)
                sys.exit(1)

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

            if vaccines_by_node is not None:
                for node_name, values in vaccines_by_node.items():
                    if node_name not in final_distances_by_node:
                        final_distances_by_node[node_name] = {}

                    final_distances_by_node[node_name][f"vaccine_strain_for_{attribute}"] = values

    # Prepare params for export.
    params = {
        "attribute": args.attribute_name,
        "map_name": distance_map_names,
        "start_date": args.start_date,
        "end_date": args.end_date,
        "interval": args.interval,
        "months_prior_to_season": args.months_prior_to_season,
        "vaccines": args.vaccines is not None
    }

    # Export distances to JSON.
    write_json({"params": params, "nodes": final_distances_by_node}, args.output)
