from . import base_interaction
from .. import utils
from .. import features
from .. import model
from .. import agent


class Injection(base_interaction.BaseInteraction):

    name = "injection"

    @staticmethod
    def get_num_acts(model: "model.TITAN", rel: "agent.Relationship") -> int:
        """
        Simulate random transmission of HIV between two PWID agents through injection.

        args:
            model: The currently running model
            rel: The relationship in which the interaction is happening
        """

        # make sure both agents have Inj drug type, should only be possible for
        # the relationship to have the injection interaction type if both agents PWID
        assert rel.agent1.drug_type == "Inj"
        assert rel.agent2.drug_type == "Inj"

        agent_params = rel.agent1.location.params.demographics[rel.agent1.race][
            rel.agent1.sex_type
        ].injection

        # TO_REVIEW should this be looking to partnership.injection.frequency?
        mean_num_acts = agent_params.num_acts * model.calibration.injection.act
        share_acts = utils.poisson(mean_num_acts, model.np_random)

        # If sharing, minimum of 1 share act, TO_REVIEW should we allow 0?
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

        return share_acts
