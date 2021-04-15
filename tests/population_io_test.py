import pytest
import os
import csv
from copy import deepcopy

from titan.population_io import (
    write,
    read,
    agent_exclude_attrs,
    agent_feature_attrs,
    agent_exposure_attrs,
    find_agent,
)
from titan.features import Prep, BaseFeature
from titan.exposures import BaseExposure


@pytest.mark.unit
def test_write_pop(tmpdir, make_population):
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
        for attr in a.__dict__.keys():
            if any(
                attr in exclude_set
                for exclude_set in (
                    agent_exclude_attrs,
                    agent_exposure_attrs,
                    agent_feature_attrs,
                )
            ):
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
def test_read_pop(tmpdir, make_population, params):
    pop = make_population(n=10)
    prep_counts = deepcopy(Prep.counts)

    write(pop, tmpdir, compress=False)

    new_pop = read(params, tmpdir)

    assert pop.id == new_pop.id
    assert pop.all_agents.num_members() == new_pop.all_agents.num_members()
    assert len(pop.relationships) == len(new_pop.relationships)
    assert prep_counts == Prep.counts

    agent = next(iter(pop.all_agents))
    new_agent = find_agent(new_pop, str(agent.id))

    attrs = agent.__dict__.keys()

    for attr in attrs:
        if attr == "component":
            continue

        orig_attr = getattr(agent, attr)
        new_attr = getattr(new_agent, attr)
        if isinstance(orig_attr, BaseFeature):
            feat_attrs = orig_attr.__dict__.keys()
            for feat_attr in feat_attrs:
                assert getattr(orig_attr, feat_attr) == getattr(new_attr, feat_attr)
        elif isinstance(orig_attr, BaseExposure):
            expose_attrs = orig_attr.__dict__.keys()
            for expose_attr in expose_attrs:
                assert getattr(orig_attr, expose_attr) == getattr(new_attr, expose_attr)
        else:
            assert orig_attr == new_attr


@pytest.mark.unit
def test_write_read_pop_compressed(tmpdir, make_population, params):
    params.prep.cap = 0.5
    pop = make_population(n=10)
    prep_counts = deepcopy(Prep.counts)

    archive = write(pop, tmpdir, compress=True)

    new_pop = read(params, archive)

    assert pop.id == new_pop.id
    assert pop.all_agents.num_members() == new_pop.all_agents.num_members()
    assert len(pop.relationships) == len(new_pop.relationships)
    assert prep_counts == Prep.counts

    agent = next(iter(pop.all_agents))
    new_agent = find_agent(new_pop, str(agent.id))

    attrs = agent.__dict__.keys()

    for attr in attrs:
        orig_attr = getattr(agent, attr)
        new_attr = getattr(new_agent, attr)
        if isinstance(orig_attr, BaseFeature):
            feat_attrs = orig_attr.__dict__.keys()
            for feat_attr in feat_attrs:
                assert getattr(orig_attr, feat_attr) == getattr(new_attr, feat_attr)
        elif isinstance(orig_attr, BaseExposure):
            expose_attrs = orig_attr.__dict__.keys()
            for expose_attr in expose_attrs:
                assert getattr(orig_attr, expose_attr) == getattr(new_attr, expose_attr)
        else:
            assert orig_attr == new_attr
