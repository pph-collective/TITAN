import random
from functools import wraps
from typing import TypeVar, Collection, Union, Iterable, Dict, Tuple, Set
from math import floor
import logging
import os
import csv
from datetime import datetime
import argparse

import oyaml as yaml  # type: ignore
import networkx as nx  # type: ignore

from . import distributions
from .parse_params import ObjMap


def memo(f):
    """
    Decorator to memoize a function
    (caches results given args, only use if deterministic)
    """
    cache = {}

    @wraps(f)
    def wrap(*arg):
        if arg not in cache:
            cache[arg] = f(*arg)
        return cache[arg]

    return wrap


def get_check_rand_int(seed: int) -> int:
    """
    Check the value passed of a seed, make sure it's an int, if 0, get a random seed

    args:
        seed: integer to check or replace with a seed

    returns:
        validated seed
    """
    if type(seed) is not int or seed < 0:
        raise ValueError("Random seed must be positive integer")
    elif seed == 0:
        return random.randint(1, 1000000)
    else:
        return seed


def safe_divide(numerator: int, denominator: int) -> float:
    """
    Divide two numbers, but default 0 if denominator is 0, otherwise divide as normal.

    args:
        numerator: number being divided
        denominator: number doing the dividing

    returns:
        resulting number
    """
    if denominator == 0:
        return 0.0
    else:
        return 1.0 * numerator / denominator


# Requirement for safe_random_choice function
T = TypeVar("T")


def safe_random_choice(seq, rand_gen, weights=None):
    """
    Return None or a random choice from a collection of items

    args:
        seq: collection to select a random item from
        rand_gen: random number generator
        weights: an optional collection of weights to use instead of a uniform distribution

    returns:
        an item, or `None` if the collection is empty
    """
    if not seq:
        return None

    if isinstance(seq, set):
        seq = tuple(seq)

    # don't call out to random choices if we don't need to (for performance)
    if len(seq) == 1:
        return seq[0]
    elif len(seq) == 2 and weights is None:
        return seq[0] if rand_gen.random() <= 0.5 else seq[1]

    choices = rand_gen.choices(seq, weights=weights)
    return choices[0]


def safe_random_int(start: int, stop: int, rand_gen) -> int:
    """
    Return an integer between [start, stop)

    args:
        start: start value
        stop: stop value
        rand_gen: random number generator

    returns:
        an item, or `None` if the collection is empty
    """
    return floor(rand_gen.random() * (stop - start) + start)


def safe_shuffle(seq: Collection[T], rand_gen) -> Iterable[T]:
    """
    Return None or a shuffled sequence

    args:
        seq: collection to shuffle
        rand_gen: random number generator

    returns:
        shuffled sequence, or `None` if empty
    """
    if seq:
        if isinstance(seq, set):
            seq = list(seq)
        rand_gen.shuffle(seq)
        return seq
    else:
        return []


@memo
def parse_var(dist_value, dist_type):
    type_caster = eval(dist_type)
    return type_caster(dist_value)


@memo
def get_dist(rand_gen, dist_type):
    if dist_type == "randint":
        return lambda *args: safe_random_int(*args, rand_gen)
    elif hasattr(rand_gen, dist_type):
        return getattr(rand_gen, dist_type)
    elif hasattr(distributions, dist_type):
        dist = getattr(distributions, dist_type)
        return lambda *args: dist(rand_gen, *args)
    else:
        raise AttributeError(f"Distribution type {dist_type} not found!")


def safe_dist(dist_info: ObjMap, rand_gen) -> Union[int, float]:
    """
    Draw a value from a distribution as defined in `dist_info`.

    args:
        dist_info: a definition of a distribution to use [params.classes.distributions]
        rand_gen: random number generator

    returns:
        a value drawn from the distribution
    """
    # gather arguments
    args = []
    for d in dist_info.vars.values():
        args.append(parse_var(d.value, d.value_type))

    dist = get_dist(rand_gen, dist_info.dist_type)
    value = dist(*args)

    if hasattr(value, "__iter__"):  # check if value is any type of sequence
        return value[0]
    else:
        return value


def binom_0(n: int, p: float):
    """
    Mirrors scipy binom.pmf as used in code
    """
    return (1 - p) ** n


def poisson(mu: float, np_rand):
    """
    Mirrors scipy poisson.rvs function as used in code
    """
    return np_rand.poisson(mu)


def get_param_from_path(params: ObjMap, param_path: str, delimiter: str):
    """
    Given a params object and a delimited path, get the leaf of the params tree
    and the last key to access it
    """
    path = param_path.split(delimiter)
    path_params = params
    for p in path[:-1]:
        try:
            path_params = path_params[p]
        except KeyError:
            path_params = path_params[int(p)]

    return path_params, path[-1]


def scale_param(params: ObjMap, param_path: str, scalar: float, delimiter="|"):
    """
    Given the params and a parameter path in the format prep|cap, scale the
    current value by the scalar
    """
    scaling_item, last_key = get_param_from_path(params, param_path, delimiter)

    old_val = scaling_item[last_key]
    logging.info(f"scaling - {param_path}: {old_val} => {old_val * scalar}")
    scaling_item[last_key] = old_val * scalar


def override_param(params: ObjMap, param_path: str, value, delimiter="|"):
    """
    Given the params and a parameter path in the format prep|cap, change the
    current value to new value
    """
    override_item, last_key = get_param_from_path(params, param_path, delimiter)
    try:
        old_val = override_item[last_key]
    except KeyError:
        last_key = int(last_key)
        old_val = override_item[last_key]

    logging.info(f"overriding - {param_path}: {old_val} => {value}")
    override_item[last_key] = value


def total_probability(p: float, num_acts: int) -> float:
    """
    Given a per act probability and a number of acts, return the total probability.

    args:
        p: the per act probability
        num_acts: the number of acts

    returns:
        the total probability
    """
    if num_acts == 1:
        return p
    elif num_acts >= 1:
        return 1.0 - binom_0(num_acts, p)
    else:
        return 0.0


def connected_components(graph):
    """
    Get connected components in graph

    args:
        graph: the model's underlying graph

    returns:
        list of connected components
    """
    return sorted(
        list(graph.subgraph(c) for c in nx.connected_components(graph)),
        key=len,
        reverse=True,
    )


def get_independent_bin(rand_gen, bin_def: ObjMap) -> int:
    """
    Get the bin key given independent bins.  A probability is selected at random, then each bin's `prob` is compared to it, the first bin that has a `prob` less than or equal to that probability is returned.

    args:
        rand_gen: A random number generator
        bin_def: The ObjMap containing the bins

    returns:
        The integer key of the matched bin (or last bin if no matches)
    """
    rand_val = rand_gen.random()
    for bin, fields in bin_def.items():
        if rand_val <= fields.prob:
            break

    return bin


def get_cumulative_bin(rand_gen, bin_def: ObjMap) -> int:
    """
    Get the bin key given cumulative bins.  A probability is selected at random, then each bin's `prob` is compared to it, the first bin that has a cumulative `prob` (e.g. for bin 2, the prob of bin 1 plus the prob of bin 2) less than or equal to that probability is returned.

    args:
        rand_gen: random number generator
        bin_def: ObjMap containing the bins

    returns:
        integer key of the matched bin (or last bin if no matches)
    """
    rand_val = rand_gen.random()
    p = 0.0
    for bin, fields in bin_def.items():
        p += fields.prob
        if rand_val <= p:
            break

    return bin


def set_up_logging(params):
    # set up logging
    if params.outputs.logging.destination == "file":
        log_msg_format = "{message:<64}  [{asctime}][{levelname}][{module}]"
        log_dt_format = "%Y%m%d%H%M%S"
        path = (
            os.getcwd()
            if params.outputs.logging.filepath == "__cwd__"
            else params.outputs.logging.filepath
        )
        logging.basicConfig(
            format=log_msg_format,
            style="{",
            datefmt=log_dt_format,
            filename=os.path.join(
                path, f"titan_log_{datetime.now().strftime(log_dt_format)}.txt"
            ),
            level=params.outputs.logging.level,
        )
    else:
        log_msg_format = "%(message)s"
        logging.basicConfig(
            format=log_msg_format,
            level=params.outputs.logging.level,
        )


# if edge exists, add it to edges - sorted by name to avoid duplicates as edges are undirected
def add_edge(edges, loc, edge):
    if edge != "" and edge != loc:
        edges.add(tuple(sorted([loc, edge])))


def grid_file_to_edges(file_path: str, diagonal_neighbors: bool = False) -> Dict:
    """
    Read a csv describing the layout of locations and return a dictionary describing the location edges, which can then be used in the params [location.edges].

    Sample csv:
    ```
    location_1,location_3,location_4
    location_1,location_3,location_4
    location_1,location_3,location_4
    location_1,location_2,
    location_2,location_2,
    ```

    Would generate the edges:
    * location_1, location_3
    * location_1, location_2
    * location_2, location_3
    * location_3, location_4

    If `diagonal_neighbors` were True, the edge [location_2, location_4] would also be returned.

    args:
        file_path: path to a csv file which contains a layout for the locations in the model.
        diagonal_neighbors: whether diagonally adjacent cells should be considered neighbors [default false]

    returns:
        A dictionary with generated edge names to locations
    """
    # read in the grid
    grid = []
    with open(file_path, newline="") as f:
        reader = csv.reader(f)
        for row in reader:
            grid.append(row)

    # generate edge pairs
    edges: Set[Tuple[str, str]] = set()
    for i in range(len(grid) - 1):
        for j in range(len(grid[0]) - 1):
            loc = grid[i][j]
            if loc == "":
                continue

            add_edge(edges, loc, grid[i + 1][j])
            add_edge(edges, loc, grid[i][j + 1])

            if diagonal_neighbors:
                add_edge(edges, loc, grid[i + 1][j + 1])
                if i >= 1:
                    add_edge(edges, loc, grid[i - 1][j + 1])
                if j >= 1:
                    add_edge(edges, loc, grid[i + 1][j - 1])

    res = {}
    for (i, edge) in enumerate(edges):
        res[f"edge_{i+1}"] = {"location_1": edge[0], "location_2": edge[1]}

    return res


def grid_file_to_edge_yml(
    file_path: str, outfile_path: str, diagonal_neighbors: bool = False
) -> None:
    """
    Read a csv describing the layout of locations and write the results to a yml file describing the location edges, which can then be used in the params [location.edges].

    Sample csv:
    ```
    location_1,location_3,location_4
    location_1,location_3,location_4
    location_1,location_3,location_4
    location_1,location_2,
    location_2,location_2,
    ```

    Would generate the edges:
    * location_1, location_3
    * location_1, location_2
    * location_2, location_3
    * location_3, location_4

    If `diagonal_neighbors` were True, the edge [location_2, location_4] would also be returned.

    args:
        file_path: path to a csv file which contains a layout for the locations in the model.
        outfile_path: path where the resulting yml file should be saved
        diagonal_neighbors: whether diagonally adjacent cells should be considered neighbors [default false]
    """
    edges = grid_file_to_edges(file_path, diagonal_neighbors=diagonal_neighbors)

    with open(outfile_path, "w") as f:
        yaml.dump({"edges": edges}, f)


# cli entry point
def grid_to_edges_cli():
    # set up args parsing
    parser = argparse.ArgumentParser(
        description="Convert CSV grid file to location edge yml"
    )
    parser.add_argument(
        "-i", "--infile", required=True, help="csv file with location grid"
    )

    parser.add_argument(
        "-o", "--outfile", required=True, help="yml file to save generated edges"
    )

    parser.add_argument(
        "-d",
        "--diagonals",
        action="store_true",
        help="Count diagonally adjacent cells as neighbors",
    )
    args = parser.parse_args()
    grid_file_to_edge_yml(
        args.infile.strip(), args.outfile.strip(), diagonal_neighbors=args.diagonals
    )
