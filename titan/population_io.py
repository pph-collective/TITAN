import os
import csv
from typing import Dict, Optional
from shutil import make_archive, unpack_archive
from tempfile import mkdtemp
import glob

from .population import Population
from .agent import Agent, Relationship
from .parse_params import ObjMap

core_attrs = [
    "id",
    "so",
    "age",
    "age_bin",
    "race",
    "drug_use",
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
exclude_attrs = ["partners", "relationships"]


def write(
    pop: Population, dir: str, intervention_attrs: bool = False, compress: bool = True
):
    """
    Write a non-empty Population to file.

    :Input:
        pop: Population - a non-empty agent population
        dir: str - path to directory where files should be written
        intervention_attrs: boolean - whether to include intervention attributions in addition to core agent attributes (less likely to be backwards compatible if used with different versions of the model) [default False]
        compress: boolean - whether to compress and archive the csvs [default True]
    :Output:
        path: str - archive name if compress is true
    """
    # open agent file
    assert len(pop.relationships) > 0, "can't write empty population"

    agent_file = os.path.join(dir, f"{pop.id}_agents.csv")

    if intervention_attrs:
        # get all attributes
        a = next(iter(pop.all_agents))
        agent_attrs = [k for k in a.__dict__.keys() if k not in exclude_attrs]
    else:
        agent_attrs = core_attrs

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


def read(params: ObjMap, path: str) -> Population:
    """
    Read a population from file and return a Population instnace

    :Input:
        params: ObjMap - the parameters used for creating this popultation
        path: str - path where <id>_agents.csv and <id>_relationships.csv are or tar.gz file containing population
    :Output:
        pop : Population
    """
    if os.path.isfile(path):
        dir = mkdtemp()
        unpack_archive(path, dir)
        path = dir

    print(glob.glob(os.path.join(path, "*")))
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
            a = create_agent(row, params.classes.bond_types.keys())
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


def create_agent(row: Dict[str, str], bond_types) -> Agent:
    init_attrs = ["so", "age", "race", "drug_use", "id"]
    agent = Agent(
        eval(row["so"]),
        eval(row["age"]),
        eval(row["race"]),
        eval(row["drug_use"]),
        eval(row["id"]),
    )

    for attr, val in row.items():
        if attr not in init_attrs:
            setattr(agent, attr, eval(val))

    agent.partners = {bond: set() for bond in bond_types}

    return agent


def create_relationship(row: Dict[str, str], pop: Population) -> Relationship:
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
    id = eval(id_str)
    for a in pop.all_agents:
        if a.id == id:
            return a

    raise Exception("No agent found")
