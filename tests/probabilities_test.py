import pytest
import os

import titan.probabilities as probs


@pytest.mark.unit
def test_adherent_prob():
    # initiate result dict with 2 time steps
    assert probs.adherent_prob(True) == 0.0008
    assert probs.adherent_prob(False) == 0.00368


@pytest.mark.unit
def test_get_death_rate(make_population):
    pop = make_population(n=0)
    location = pop.geography.locations["world"]
    for hiv in [True, False]:
        for aids in [True, False]:
            for race in ["white", "black"]:
                for drug_type in ["None", "Inj"]:
                    for adh in [0, 5]:
                        assert (
                            probs.get_death_rate(
                                hiv,
                                aids,
                                drug_type,
                                adh,
                                race,
                                location,
                                location.params.model.time.steps_per_year,
                            )
                            > 0
                        )
