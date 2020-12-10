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
        if model.time < model.params.hiv.start_time:
            return False

        # Agent 1 is HIV+, Agent 2 is not, Agent 2 is succept
        if rel.agent1.hiv and not rel.agent2.hiv:
            agent = rel.agent1
            partner = rel.agent2
        # If Agent 2 is HIV and Agent 1 is not, Agent 1 is succept
        elif not rel.agent1.hiv and rel.agent2.hiv:
            agent = rel.agent2
            partner = rel.agent1
        else:  # neither agent is HIV or both are
            return False

        # Everything from here is only run if one of them is HIV+

        # unprotected sex probabilities for primary partnerships
        mean_sex_acts = (
            rel.get_number_of_sex_acts(model.np_random) * model.calibration.sex.act
        )
        total_sex_acts = utils.poisson(mean_sex_acts, model.np_random)

        # Get condom usage
        p_safe_sex = (
            agent.location.params.demographics[agent.race]
            .sex_type[agent.sex_type]
            .safe_sex
        )
        # increase condom usage if diagnosed
        if agent.hiv_dx or partner.hiv_dx:
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

        if unsafe_sex_acts >= 1:
            # agent is HIV+
            rel.total_sex_acts += unsafe_sex_acts
            p_per_act = model.get_transmission_probability("sex", agent, partner)

            p_total_transmission: float
            if unsafe_sex_acts == 1:
                p_total_transmission = p_per_act
            else:
                p_total_transmission = 1.0 - utils.binom_0(unsafe_sex_acts, p_per_act)

            if model.run_random.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                model.hiv_convert(partner)

        return True
