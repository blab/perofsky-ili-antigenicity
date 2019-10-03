"""
Utilities shared by scripts
"""
from augur.frequency_estimators import timestamp_to_float
import pandas as pd


def get_seasons(start_date, end_date, interval):
    """Return a list of datetime values corresponding to the seasons between the
    given start and end dates where the interval defines the season spacing.

    Arguments
    ---------
    start_date, end_date : str
        datetime formatted string like "YYYY-MM-DD"

    interval : int
        number of months between seasons

    Returns
    -------
    list :
        list of pandas datetime instances
    """
    # Create season intervals.
    seasons = pd.date_range(start_date, end_date, freq="%iMS" % interval)

    return seasons


def get_nodes_by_season(tree, seasons, terminal_only=True):
    """Assign the nodes in the given tree to the given seasons based on the
    collection date annotated in the "num_date" key of the `attr` attribute for
    each node.

    Arguments
    ---------
    tree : Bio.Phylo
        a phylogenetic tree instance from BioPython

    seasons : list
        a list of pandas datetime instances

    Returns
    -------
    dict :
        a map of node instances to season date strings (in YYYY-MM-DD format)
    """
    nodes_assigned = set()
    nodes_by_season = {}
    terminal = True if terminal_only else None

    for season in seasons:
        season_date = season.strftime("%Y-%m-%d")
        season_float = timestamp_to_float(season)
        nodes_by_season[season_date] = []

        for node in tree.find_clades(terminal=terminal):
            if node.name not in nodes_assigned and node.attr["num_date"] < season_float:
                nodes_assigned.add(node.name)
                nodes_by_season[season_date].append(node)

    return nodes_by_season
