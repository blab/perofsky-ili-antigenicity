import argparse
import Bio
import Bio.Phylo
import json
import pandas as pd


# Biopython's trees don't store links to node parents, so we need to build
# a map of each node to its parent.
# Code from the Bio.Phylo cookbook: http://biopython.org/wiki/Phylo_cookbook
def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents


def annotate_parents(tree):
    # Get all parent nodes by node.
    parents_by_node = all_parents(tree)

    # Next, annotate each node with its parent.
    for node in tree.find_clades():
        if node == tree.root:
            node.up = None
        else:
            node.up = parents_by_node[node]

    # Return the tree.
    return tree


def json_to_tree(json_dict, root=True):
    """Returns a Bio.Phylo tree corresponding to the given JSON dictionary exported
    by `tree_to_json`.

    Assigns links back to parent nodes for the root of the tree.

    >>> import json
    >>> json_fh = open("tests/data/json_tree_to_nexus/flu_h3n2_ha_3y_tree.json", "r")
    >>> json_dict = json.load(json_fh)
    >>> tree = json_to_tree(json_dict)
    >>> tree.name
    u'NODE_0002020'
    >>> len(tree.clades)
    2
    >>> tree.clades[0].name
    u'NODE_0001489'
    >>> hasattr(tree, "attr")
    True
    >>> "dTiter" in tree.attr
    True
    """
    node = Bio.Phylo.Newick.Clade()
    node.name = json_dict["strain"]

    if "children" in json_dict:
        # Recursively add children to the current node.
        node.clades = [json_to_tree(child, root=False) for child in json_dict["children"]]

    # Assign all non-children attributes.
    for attr, value in json_dict.items():
        if attr != "children":
            setattr(node, attr, value)

    node.numdate = node.attr.get("num_date")
    node.branch_length = node.attr.get("div")

    if "translations" in node.attr:
        node.translations = node.attr["translations"]

    if root:
        node = annotate_parents(node)

    return node


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice tree JSON")
    parser.add_argument("output", help="tab-delimited file of attributes per node of the given tree")
    parser.add_argument("--include-internal-nodes", action="store_true", help="include data from internal nodes in output")
    parser.add_argument("--attributes", nargs="+", help="names of attributes to export from the given tree")

    args = parser.parse_args()

    # Load tree from JSON.
    with open(args.tree, "r") as fh:
        tree_json = json.load(fh)

    tree = json_to_tree(tree_json)

    # Collect attributes per node from the tree to export.
    records = []
    for node in tree.find_clades():
        if node.is_terminal() or args.include_internal_nodes:
            record = {
                "name": node.name
            }

            for attribute in args.attributes:
                record[attribute] = node.attr.get(attribute, "N/A")

            records.append(record)

    # Convert records to a data frame and save as a tab-delimited file.
    df = pd.DataFrame(records)
    df.to_csv(args.output, sep="\t", header=True, index=False, columns=["name"] + list(args.attributes), float_format="%.2f")
