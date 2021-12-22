import pytest
import os

import titan.probabilities as probs


@pytest.mark.unit
def test_get_death_rate(make_population):
    pop = make_population(n=0)
    location = pop.geography.locations["world"]
    for hiv in [True, False]:
        for aids in [True, False]:
            for race in ["white", "black"]:
                for drug_type in ["None", "Inj"]:
                    for sex_type in ["MSM", "HM", "HF", "MTF", "WSW"]:
                        for adh in [0, 5]:
                            assert (
                                probs.get_death_rate(
                                    hiv,
                                    aids,
                                    drug_type,
                                    sex_type,
                                    adh,
                                    race,
                                    location,
                                    location.params.model.time.steps_per_year,
                                    "death",
                                )
                                > 0
                            )
