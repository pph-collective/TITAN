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

        # make sure both agents have Inj drug type, should only be possible for
        # the relationship to have the injection interaction type if both agents PWID
        assert rel.agent1.drug_type == "Inj"
        assert rel.agent2.drug_type == "Inj"

        agent_params = rel.agent1.location.params.demographics[rel.agent1.race][
            rel.agent1.sex_type
        ].injection

        mean_num_acts = agent_params.num_acts * model.calibration.injection.act
        share_acts = utils.poisson(mean_num_acts, model.np_random)

        # If sharing, minimum of 1 share act
        if share_acts < 1:
            share_acts = 1

        if (
            rel.agent1.syringe_services.active or rel.agent2.syringe_services.active  # type: ignore[attr-defined]
        ):  # syringe services program risk
            p_unsafe_injection = features.SyringeServices.enrolled_risk
        else:
            p_unsafe_injection = agent_params.unsafe_prob

            # diagnosis risk reduction
            if rel.agent1.hiv.dx or rel.agent1.hiv.dx:  # type: ignore[attr-defined]
                p_unsafe_injection *= 1 - model.params.hiv.dx.risk_reduction.injection

        for n in range(share_acts):
            if model.run_random.random() > p_unsafe_injection:
                share_acts -= 1

        if share_acts > 0:
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

                agent_exposure = getattr(agent, exposure.name)
                p_per_act = agent_exposure.get_transmission_probability(
                    model, "injection", partner
                )

                p_total_transmission: float
                if share_acts == 1:
                    p_total_transmission = p_per_act
                else:
                    p_total_transmission = 1.0 - utils.binom_0(share_acts, p_per_act)

                if model.run_random.random() < p_total_transmission:
                    # if agent HIV+ partner becomes HIV+
                    partner_exposure = getattr(partner, exposure.name)
                    partner_exposure.convert(model)

        return True
