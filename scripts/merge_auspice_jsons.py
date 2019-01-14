"""
Merge metadata and tree JSONs from auspice into a single JSON for use in auspice static sites.
"""
import argparse
import json


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="auspice tree JSON")
    parser.add_argument("metadata", help="auspice metadata JSON")
    parser.add_argument("output", help="merged auspice JSON")

    args = parser.parse_args()

    # Load tree.
    with open(args.tree, "r") as fh:
        tree_json = json.load(fh)

    # Load metadata.
    with open(args.metadata, "r") as fh:
        metadata_json = json.load(fh)

    # Write a merged JSON.
    merged_json = {
        "meta": metadata_json,
        "tree": tree_json
    }

    with open(args.output, "w") as oh:
        json.dump(merged_json, oh, indent=1, sort_keys=True)
