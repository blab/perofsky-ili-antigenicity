"""Concatenate two or more tables as data frames.
"""
import argparse
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables", nargs="+", help="tables to concatenate")
    parser.add_argument("--separator", default="\t", help="separator between columns in the given tables")
    parser.add_argument("--id-column", default="replicate", help="name of the column to store ids for the given tables")
    parser.add_argument("--ids", nargs="+", help="id to associated with records from each of the given tables. Defaults to the original table filename if no ids are given.")
    parser.add_argument("--sort-by", help="what to sort the dataframe by (optional)")
    parser.add_argument("--output", help="concatenated table")

    args = parser.parse_args()

    # Assign ids for each table to help distinguish sources in the merged output.
    ids = args.tables
    if args.ids:
        ids = args.ids

    id_by_table_number = dict(enumerate(ids))

    # Concatenate tables.
    df = pd.concat([
        pd.read_csv(table_file, sep=args.separator).assign(_table_number=table_number)
        for table_number, table_file in enumerate(args.tables)
    ], ignore_index=True, sort=True)

    # Map ids to the requested column by table number.
    df[args.id_column] = df["_table_number"].map(id_by_table_number)
    df = df.drop(columns=["_table_number"])

    if args.sort_by is not None:
        df = df.sort_values(by=[args.sort_by])
        cols_to_order = [args.sort_by]
        new_columns = cols_to_order + (df.columns.drop(cols_to_order).tolist())
        df = df[new_columns]

    df.to_csv(args.output, sep=args.separator, header=True, index=False)
