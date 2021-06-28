import argparse
from augur.utils import json_to_tree
import Bio
import Bio.Phylo
import json
import pandas as pd
import sys


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice tree JSON")
    parser.add_argument("output", help="tab-delimited file of attributes per node of the given tree")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data from internal nodes in output")
    parser.add_argument("--attributes", nargs="+", help="names of attributes to export from the given tree")
    parser.add_argument("--annotations", nargs="+", help="additional annotations to add to the output table in the format of 'key=value' pairs")

    args = parser.parse_args()

    # Load tree from JSON.
    with open(args.tree, "r") as fh:
        tree_json = json.load(fh)

    tree = json_to_tree(tree_json)

    if args.attributes:
        attributes = args.attributes
    else:
        attributes = sorted(list(tree.root.node_attr.keys()) + list(tree.root.branch_attrs.keys()))

    records = []
    for node in tree.find_clades():
        if node.is_terminal() or args.include_internal_nodes:
            record = {
                "name": node.name
            }

            for attribute in attributes:
                if attribute in node.node_attrs:
                    record[attribute] = node.node_attrs[attribute]["value"]
                elif attribute in node.branch_attrs:
                    record[attribute] = node.branch_attrs[attribute]["value"]
                else:
                    print(f"Attribute '{attribute}' missing from node '{node.name}'", file=sys.stderr)

            records.append(record)

    # Convert records to a data frame.
    df = pd.DataFrame(records)

    # Add any additional annotations requested by the user in the format of
    # "key=value" pairs where each key becomes a new column with the given
    # value.
    if args.annotations:
        annotation_columns = []
        for annotation in args.annotations:
            key, value = annotation.split("=")
            df[key] = value
            annotation_columns.append(key)

    # Save as a tab-delimited file.
    df.to_csv(
        args.output,
        sep="\t",
        header=True,
        index=False,
        columns=["name"] + list(args.attributes) + annotation_columns,
        float_format="%.2f"
    )
