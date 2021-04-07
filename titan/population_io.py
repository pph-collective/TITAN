import os
import csv
from typing import Dict, Any
from shutil import make_archive, unpack_archive
from tempfile import mkdtemp
import glob
import re
import logging

from .population import Population
from .agent import Agent, Relationship
from .parse_params import ObjMap
from .location import Location
from . import features
from . import exposures
from . import utils

agent_feature_attrs = [
    feature.name for feature in features.BaseFeature.__subclasses__()
]
agent_exposure_attrs = [
    exposure.name for exposure in exposures.BaseExposure.__subclasses__()
]

# these are functionally saved in the relationships or other files and complicate the agent file
agent_exclude_attrs = (
    {"partners", "relationships"}.union(agent_feature_attrs).union(agent_exposure_attrs)
)


def write(pop: Population, dir: str, compress: bool = True) -> str:
    """
    Write a non-empty Population to file.

    args:
        pop: a non-empty agent population
        dir: path to directory where files should be written
        compress: whether to compress and archive the csv

    returns:
        path, or archive name if compress is true
    """
    assert len(pop.relationships) > 0, "Can't write empty population"

    utils.set_up_logging(pop.params)

    # open agent file
    agent_file = os.path.join(dir, f"{pop.id}_agents.csv")

    a = next(iter(pop.all_agents))
    # get all attributes
    agent_attrs = [k for k in a.__dict__.keys() if k not in agent_exclude_attrs]

    write_class_file(agent_file, pop.all_agents, agent_attrs)

    extra_files = []

    # write agent extras (features, exposures) to their own files
    def write_extra_class(extra_attrs, extra_type):
        for extra in extra_attrs:
            extra_obj = getattr(a, extra)
            extra_attrs = list(extra_obj.__dict__.keys())
            extra_file = os.path.join(dir, f"{pop.id}_{extra_type}_{extra}.csv")
            extra_files.append(extra_file)
            write_extra_class_file(extra_file, pop.all_agents, extra, extra_attrs)

    write_extra_class(agent_feature_attrs, "feat")
    write_extra_class(agent_exposure_attrs, "exposure")

    # open relationship file
    rel_file = os.path.join(dir, f"{pop.id}_relationships.csv")

    r = next(iter(pop.relationships))
    rel_attrs = list(r.__dict__.keys())

    write_class_file(rel_file, pop.relationships, rel_attrs)

    if compress:
        archive_name = make_archive(
            os.path.join(dir, f"{pop.id}_pop"), "gztar", root_dir=dir, base_dir="."
        )
        os.remove(agent_file)
        os.remove(rel_file)
        for f in extra_files:
            os.remove(f)

        return archive_name
    else:
        return dir


def write_extra_class_file(file_name, collection, extra, attrs):
    logging.info(f"Creating {file_name}")
    with open(file_name, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=attrs)
        writer.writeheader()
        for item in collection:
            feat = getattr(item, extra)
            writer.writerow({attr: repr(getattr(feat, attr)) for attr in attrs})


def write_class_file(file_name, collection, attrs):
    logging.info(f"Creating {file_name}")
    with open(file_name, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=attrs)
        writer.writeheader()
        for item in collection:
            writer.writerow({attr: repr(getattr(item, attr)) for attr in attrs})


def read(params: ObjMap, path: str) -> Population:
    """
    Read a population from file and return a Population instance

    args:
        params: the parameters used for creating this popultation
        path: path where [id]_agents.csv and [id]_relationships.csv are or tar.gz file
            containing population
    returns:
        the re-constituted population
    """
    if os.path.isfile(path):
        dir = mkdtemp()
        unpack_archive(path, dir)
        path = dir

    agent_file = glob.glob(os.path.join(path, "*_agents.csv"))[0]
    rel_file = glob.glob(os.path.join(path, "*_relationships.csv"))[0]
    feat_files = glob.glob(os.path.join(path, "*_feat_*.csv"))
    exposure_files = glob.glob(os.path.join(path, "*_exposure_*.csv"))
    assert os.path.isfile(agent_file), f"can't find agents.csv in {dir}"
    assert os.path.isfile(rel_file), f"can't find relationships.csv in {dir}"

    _, agent_filename = os.path.split(agent_file)
    id = agent_filename[:8]

    # create feature dict
    agent_extras: Dict[str, Dict] = {}

    def update_agent_extras(files, extra_type):
        pattern = re.compile(f"^.*_{extra_type}_(.*)\\.csv$")
        for file in files:
            m = pattern.match(file)
            if m is not None:
                extra = m.group(1)
                agent_extras[extra] = {}
                with open(file, newline="") as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        agent_extras[extra][int(row["agent"])] = row

    update_agent_extras(feat_files, "feat")
    update_agent_extras(exposure_files, "exposure")

    # don't create any agents on init
    params.model.num_pop = 0
    pop = Population(params, id=id)

    # re-create all agents and add to population
    with open(agent_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            a = create_agent(
                row,
                params.classes.bond_types.keys(),
                pop.geography.locations,
                agent_extras,
            )
            pop.add_agent(a)

    # update num_pop to actual population
    params.model.num_pop = pop.all_agents.num_members()

    # re-create all relationships and add to population
    with open(rel_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = create_relationship(row, pop)
            pop.add_relationship(r)

    pop.update_agent_components()

    return pop


def create_agent(
    row: Dict[str, str],
    bond_types,
    locations: Dict[str, Location],
    agent_extras: Dict[str, Any],
) -> Agent:
    """
    Initialize an Agent from a row of the saved population
    """
    init_attrs = ["sex_type", "age", "race", "drug_type", "id", "location"]
    location = locations[eval(row["location"])]
    agent = Agent(
        eval(row["sex_type"]),
        eval(row["age"]),
        eval(row["race"]),
        eval(row["drug_type"]),
        location,
        eval(row["id"]),
    )

    for attr, val in row.items():
        if attr not in init_attrs:
            setattr(agent, attr, eval(val))

    for extra in agent_extras:
        extra_row = agent_extras[extra][agent.id]
        agent_extra = getattr(agent, extra)
        for attr, val in extra_row.items():
            if not attr == "agent":
                setattr(agent_extra, attr, eval(val))

        if agent_extra.active:
            agent_extra.add_agent(agent)

    agent.partners = {bond: set() for bond in bond_types}

    return agent


def create_relationship(row: Dict[str, str], pop: Population) -> Relationship:
    """
    Initialize a Relationship from a row of the saved population
    """
    init_attrs = ["agent1", "agent2", "duration", "bond_type", "id"]
    agent1 = find_agent(pop, row["agent1"])
    agent2 = find_agent(pop, row["agent2"])
    rel = Relationship(
        agent1,
        agent2,
        eval(row["duration"]),
        eval(row["bond_type"]),
        id=eval(row["id"]),
    )

    for attr, val in row.items():
        if attr not in init_attrs:
            setattr(rel, attr, eval(val))

    return rel


def find_agent(pop: Population, id_str: str) -> Agent:
    """
    Given a Population and an id (as a string), return the Agent with that id
    """
    id = eval(id_str)
    for a in pop.all_agents:
        if a.id == id:
            return a

    raise Exception(f"Agent {id_str} not found")
