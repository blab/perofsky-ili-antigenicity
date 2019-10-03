"""Calculate the distance between sequences between seasons.
"""
import argparse

from augur.lbi import calculate_LBI
from augur.utils import annotate_parents_for_tree, read_node_data, write_json

import Bio.Phylo
import pandas as pd

from utils import get_seasons, get_nodes_by_season


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree", help="Newick tree", required=True)
    parser.add_argument("--date-annotations", help="JSON of branch lengths and date annotations from augur refine for samples in the given tree; required for comparisons to earliest or latest date", required=True)
    parser.add_argument("--start-date", help="date to start seasonal intervals (e.g., 2000-10-01)", required=True)
    parser.add_argument("--end-date", help="date to end seasonal intervals (e.g., 2010-10-01)", required=True)
    parser.add_argument("--interval", type=int, help="number of months per season", required=True)
    parser.add_argument("--output", help="tab-delimited file with mean LBI by season", required=True)
    parser.add_argument("--strain-output", help="tab-delimited file with LBI per strain and season", required=True)

    args = parser.parse_args()

    # Load tree and annotate parents.
    tree = Bio.Phylo.read(args.tree, "newick")
    tree = annotate_parents_for_tree(tree)

    # Create season intervals.
    seasons = get_seasons(args.start_date, args.end_date, args.interval)

    # Load date annotations and annotate tree with them.
    date_annotations = read_node_data(args.date_annotations)
    for node in tree.find_clades():
        node.attr = date_annotations["nodes"][node.name]
        node.attr["num_date"] = node.attr["numdate"]

    # Assign nodes to seasons.
    nodes_by_season = get_nodes_by_season(tree, seasons, terminal_only=False)

    # Store the LBI per strain and season.
    strain_and_season_lbi = []

    for current_season in seasons:
        current_season_date = current_season.strftime("%Y-%m-%d")
        for node in tree.find_clades():
            # Nodes are alive if they are assigned to the current season.
            node.alive = node in nodes_by_season[current_season_date]

            # Clear any existing LBI annotations.
            if hasattr(node, "lbi"):
                del node.lbi

        # Calculate unnormalized LBI for the nodes that are alive.
        calculate_LBI(tree, normalize=False)

        # Collect LBI annotations for tips alive in this season.
        for node in nodes_by_season[current_season_date]:
            if node.is_terminal():
                strain_and_season_lbi.append({
                    "strain": node.name,
                    "season": current_season_date,
                    "lbi": node.lbi
                })

    # Export strain LBIs to a table.
    df = pd.DataFrame(strain_and_season_lbi)
    df.to_csv(args.strain_output, sep="\t", float_format="%.4f", index=False)

    # Summarize LBI by season.
    lbi_by_season = df.groupby("season")["lbi"].aggregate(["mean", "std", "count"]).reset_index()
    lbi_by_season.to_csv(args.output, sep="\t", float_format="%.4f", index=False)
