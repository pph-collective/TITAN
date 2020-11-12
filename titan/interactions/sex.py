from . import base_interaction
from .. import utils
from .. import model
from .. import agent


class Sex(base_interaction.BaseInteraction):

    name = "sex"

    @staticmethod
    def get_num_acts(model: "model.HIVModel", rel: "agent.Relationship") -> int:
        """
        Simulate random transmission of HIV between two agents through Sex. One of the agents must be HIV+.

        args:
            model: The model being run
            rel : Relationship
        """
        # unprotected sex probabilities for primary partnerships
        mean_sex_acts = (
            rel.agent1.get_number_of_sex_acts(model.np_random)
            * model.calibration.sex.act
        )
        total_sex_acts = utils.poisson(mean_sex_acts, model.np_random)

        # Get condom usage
        p_safe_sex = rel.agent1.location.params.demographics[rel.agent1.race][
            rel.agent1.sex_type
        ].safe_sex
        # increase condom usage if diagnosed
        if rel.agent1.hiv.dx or rel.agent1.hiv.dx:  # type: ignore[attr-defined]
            # Calculate probability of safe sex given risk reduction
            p_unsafe_sex = (1 - p_safe_sex) * (
                1 - model.params.hiv.dx.risk_reduction.sex
            )
            p_safe_sex = 1 - p_unsafe_sex

        # Reduction of risk acts between partners for condom usage
        unsafe_sex_acts = total_sex_acts
        for n in range(unsafe_sex_acts):
            if model.run_random.random() < p_safe_sex:
                unsafe_sex_acts -= 1

        rel.total_sex_acts += unsafe_sex_acts  # TO_REVIEW should this be total_sex_acts?, also this attribute seems to be unused...

        return unsafe_sex_acts
