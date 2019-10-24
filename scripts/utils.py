"""
Utilities shared by scripts
"""
from augur.frequency_estimators import timestamp_to_float
from collections import defaultdict
import numpy as np
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
    nodes_by_season = defaultdict(list)
    terminal = True if terminal_only else None

    # Create a list of season floating point dates to search.
    season_floats = [timestamp_to_float(season) for season in seasons]

    # Create a list of season keys to track nodes by including the start and end date.
    season_keys = [
        (season_start.strftime("%Y-%m-%d"), season_end.strftime("%Y-%m-%d"))
        for season_start, season_end in zip(seasons[:-1], seasons[1:])
    ]

    for node in tree.find_clades(terminal=terminal):
        # Find the season for the current node. This is the insertion index of
        # the node in the array of season float values minus one, so we get the
        # start date of the season.
        index = np.searchsorted(season_floats, node.attr["num_date"]) - 1
        if 0 <= index < len(season_keys):
            nodes_by_season[season_keys[index]].append(node)

    return nodes_by_season
