import pytest
import os
import csv

from titan.population_io import write, read, core_attrs, exclude_attrs, find_agent


@pytest.mark.unit
def test_write_pop_core(tmpdir, make_population):
    pop = make_population(n=10)

    write(pop, tmpdir, compress=False)

    agent_file = os.path.join(tmpdir, f"{pop.id}_agents.csv")
    assert os.path.isfile(agent_file)

    rel_file = os.path.join(tmpdir, f"{pop.id}_relationships.csv")
    assert os.path.isfile(rel_file)

    # agents writing
    res = []
    with open(agent_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            res.append(row)

    for row, a in zip(res, pop.all_agents):
        for attr in core_attrs:
            assert repr(getattr(a, attr)) == row[attr]

    # relationship writing
    res = []
    with open(rel_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            res.append(row)

    assert len(res) > 0
    for row, r in zip(res, pop.relationships):
        for attr in r.__dict__.keys():
            assert repr(getattr(r, attr)) == row[attr]


@pytest.mark.unit
def test_write_pop_intervention(tmpdir, make_population):
    pop = make_population(n=10)

    write(pop, tmpdir, intervention_attrs=True, compress=False)

    agent_file = os.path.join(tmpdir, f"{pop.id}_agents.csv")
    assert os.path.isfile(agent_file)

    rel_file = os.path.join(tmpdir, f"{pop.id}_relationships.csv")
    assert os.path.isfile(rel_file)

    # agents writing
    res = []
    with open(agent_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            res.append(row)

    for row, a in zip(res, pop.all_agents):
        for attr in a.__dict__.keys():
            if attr in exclude_attrs:
                assert attr not in row
            else:
                assert repr(getattr(a, attr)) == row[attr]

    # relationship writing
    res = []
    with open(rel_file, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            res.append(row)

    for row, r in zip(res, pop.relationships):
        for attr in r.__dict__.keys():
            assert repr(getattr(r, attr)) == row[attr]


@pytest.mark.unit
def test_read_pop_core(tmpdir, make_population, params):
    pop = make_population(n=10)

    write(pop, tmpdir, compress=False)

    new_pop = read(params, tmpdir)

    assert pop.id == new_pop.id
    assert pop.all_agents.num_members() == new_pop.all_agents.num_members()
    assert len(pop.relationships) == len(new_pop.relationships)

    agent = next(iter(pop.all_agents))
    new_agent = find_agent(new_pop, str(agent.id))

    for attr in core_attrs:
        assert getattr(agent, attr) == getattr(new_agent, attr)

    assert agent.relationships == new_agent.relationships
    assert agent.partners == new_agent.partners


@pytest.mark.unit
def test_read_pop_intervention(tmpdir, make_population, params):
    pop = make_population(n=10)

    write(pop, tmpdir, intervention_attrs=True, compress=False)

    new_pop = read(params, tmpdir)

    assert pop.id == new_pop.id
    assert pop.all_agents.num_members() == new_pop.all_agents.num_members()
    assert len(pop.relationships) == len(new_pop.relationships)

    agent = next(iter(pop.all_agents))
    new_agent = find_agent(new_pop, str(agent.id))

    attrs = agent.__dict__.keys()

    for attr in attrs:
        assert getattr(agent, attr) == getattr(new_agent, attr)


@pytest.mark.unit
def test_write_read_pop_compressed(tmpdir, make_population, params):
    pop = make_population(n=10)

    archive = write(pop, tmpdir, intervention_attrs=True, compress=True)

    new_pop = read(params, archive)

    assert pop.id == new_pop.id
    assert pop.all_agents.num_members() == new_pop.all_agents.num_members()
    assert len(pop.relationships) == len(new_pop.relationships)

    agent = next(iter(pop.all_agents))
    new_agent = find_agent(new_pop, str(agent.id))

    attrs = agent.__dict__.keys()

    for attr in attrs:
        assert getattr(agent, attr) == getattr(new_agent, attr)
