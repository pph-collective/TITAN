import os
import csv
from typing import Dict
from shutil import make_archive, unpack_archive
from tempfile import mkdtemp
import glob

from .population import Population
from .agent import Agent, Relationship
from .parse_params import ObjMap
from .location import Location

# These attributes are the non-intervention attributes of an agent.  They are considered
# "core" as they are assigned in creating an agent and are stable over time (likely
# backwards compatible)
agent_core_attrs = [
    "id",
    "sex_type",
    "age",
    "age_bin",
    "race",
    "drug_type",
    "location",
    "msmw",
    "sex_role",
    "mean_num_partners",
    "target_partners",
    "hiv",
    "hiv_time",
    "hiv_dx",
    "aids",
]

# these are functionally saved in the relationships file and complicate the agent file
agent_exclude_attrs = ["partners", "relationships"]


def write(
    pop: Population, dir: str, intervention_attrs: bool = False, compress: bool = True
) -> str:
    """
    Write a non-empty Population to file.

    args:
        pop: a non-empty agent population
        dir: path to directory where files should be written
        intervention_attrs: whether to include intervention attributions in addition to
            core agent attributes (less likely to be backwards compatible if used with
            different versions of the model)
        compress: whether to compress and archive the csv

    returns:
        path, or archive name if compress is true
    """
    assert len(pop.relationships) > 0, "can't write empty population"

    # open agent file
    agent_file = os.path.join(dir, f"{pop.id}_agents.csv")

    if intervention_attrs:
        # get all attributes
        a = next(iter(pop.all_agents))
        agent_attrs = [k for k in a.__dict__.keys() if k not in agent_exclude_attrs]
    else:
        agent_attrs = agent_core_attrs

    with open(agent_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=agent_attrs)
        writer.writeheader()
        for agent in pop.all_agents:
            writer.writerow({attr: repr(getattr(agent, attr)) for attr in agent_attrs})

    # open relationship file
    rel_file = os.path.join(dir, f"{pop.id}_relationships.csv")

    r = next(iter(pop.relationships))
    rel_attrs = list(r.__dict__.keys())

    with open(rel_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rel_attrs)
        writer.writeheader()
        for rel in pop.relationships:
            writer.writerow({attr: repr(getattr(rel, attr)) for attr in rel_attrs})

    if compress:
        archive_name = make_archive(
            os.path.join(dir, f"{pop.id}_pop"), "gztar", root_dir=dir, base_dir="."
        )
        os.remove(agent_file)
        os.remove(rel_file)
        return archive_name
    else:
        return dir


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
    assert os.path.isfile(agent_file), f"can't find agents.csv in {dir}"
    assert os.path.isfile(rel_file), f"can't find relationships.csv in {dir}"

    _, agent_filename = os.path.split(agent_file)
    id = agent_filename[:8]

    # don't create any agents on init
    params.model.num_pop = 0
    pop = Population(params, id=id)

    # re-create all agents and add to population
    with open(agent_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            a = create_agent(
                row, params.classes.bond_types.keys(), pop.geography.locations
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

    return pop


def create_agent(
    row: Dict[str, str], bond_types, locations: Dict[str, Location]
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
