import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--distances", required=True, help="table of distances to summarize across replicates")
    parser.add_argument("--output", required=True, help="table of mean, std, and count for distances across replicates")
    args = parser.parse_args()

    df = pd.read_csv(args.distances, sep="\t")

    group_columns = [
        "current_season_start",
        "current_season_end",
        "other_season_start",
        "other_season_end",
        "distance_map"
    ]
    grouped = df.groupby(
        group_columns
    ).aggregate({
        "mean_distance": ["mean", "std", "count"]
    }).reset_index()
    grouped.columns = group_columns + grouped.columns.get_level_values(1)[-3:].values.tolist()

    grouped.to_csv(
        args.output,
        sep="\t",
        index=False,
    )
