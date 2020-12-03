from . import base_interaction
from .. import utils
from .. import features
from .. import model
from .. import agent


class Injection(base_interaction.BaseInteraction):

    name = "injection"

    @staticmethod
    def interact(model: "model.HIVModel", rel: "agent.Relationship") -> bool:
        """
        Simulate random transmission of HIV between two PWID agents through injection.

        args:
            model: The currently running model
            rel: The relationship in which the interaction is happening

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

        # make sure both agents have Inj drug type, should only be possible for
        # the relationship to have the injection interaction type if both agents PWID
        assert agent.drug_type == "Inj"
        assert partner.drug_type == "Inj"

        agent_params = agent.location.params.demographics[agent.race][
            agent.sex_type
        ].injection

        mean_num_acts = agent_params.num_acts * model.calibration.injection.act
        share_acts = utils.poisson(mean_num_acts, model.np_random)

        if (
            agent.syringe_services.active or partner.syringe_services.active  # type: ignore[attr-defined]
        ):  # syringe services program risk
            p_unsafe_injection = features.SyringeServices.enrolled_risk
        else:
            # TO_REVIEW should this be outside the else?
            # If sharing, minimum of 1 share act
            if share_acts < 1:
                share_acts = 1

            p_unsafe_injection = agent_params.unsafe_prob

            if agent.hiv_dx or partner.hiv_dx:  # diagnosis risk reduction
                p_unsafe_injection *= 1 - model.params.hiv.dx.risk_reduction.injection

        for n in range(share_acts):
            if model.run_random.random() > p_unsafe_injection:
                share_acts -= 1

        if share_acts > 0:
            p = model.get_transmission_probability("injection", agent, partner)

            p_total_transmission: float
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1.0 - utils.binom_0(share_acts, p)

            if model.run_random.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                model.hiv_convert(partner)

        return True
