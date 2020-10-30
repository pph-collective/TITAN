from . import base_interaction
from .. import utils
from .. import model
from .. import agent


class Sex(base_interaction.BaseInteraction):

    name = "sex"

    @staticmethod
    def interact(model: "model.HIVModel", rel: "agent.Relationship") -> bool:
        """
        Simulate random transmission of HIV between two agents through Sex. One of the agents must be HIV+.

        args:
            model: The model being run
            rel : Relationship

        returns:
            whether the agents interacted
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
            p_safe_sex *= 1 - p_unsafe_sex

        # Reduction of risk acts between partners for condom usage
        unsafe_sex_acts = total_sex_acts
        for n in range(unsafe_sex_acts):
            if model.run_random.random() < p_safe_sex:
                unsafe_sex_acts -= 1

        rel.total_sex_acts += unsafe_sex_acts

        if unsafe_sex_acts >= 1:
            for exposure in model.exposures:
                if model.time < model.params[exposure.name].start_time:
                    continue

                agent1_exposure = getattr(rel.agent1, exposure.name)
                agent2_exposure = getattr(rel.agent2, exposure.name)

                # Agent 1 is HIV+, Agent 2 is not, Agent 2 is succept
                if agent1_exposure.active and not agent2_exposure.active:  # type: ignore[attr-defined]
                    agent = rel.agent1
                    partner = rel.agent2
                # If Agent 2 is HIV and Agent 1 is not, Agent 1 is succept
                elif not agent1_exposure.active and agent2_exposure.active:  # type: ignore[attr-defined]
                    agent = rel.agent2
                    partner = rel.agent1
                else:  # neither agent is HIV or both are
                    continue

                # Everything from here is only run if one of them is exposure active
                agent_exposure = getattr(agent, exposure.name)
                p_per_act = agent_exposure.get_transmission_probability(
                    model, "sex", partner
                )

                p_total_transmission: float
                if unsafe_sex_acts == 1:
                    p_total_transmission = p_per_act
                else:
                    p_total_transmission = 1.0 - utils.binom_0(
                        unsafe_sex_acts, p_per_act
                    )

                if model.run_random.random() < p_total_transmission:
                    # if agent HIV+ partner becomes HIV+
                    partner_exposure = getattr(partner, exposure.name)
                    partner_exposure.convert(model)

        return True
