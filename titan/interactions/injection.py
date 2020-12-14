from . import base_interaction
from .. import utils
from .. import features
from .. import model
from .. import agent


class Injection(base_interaction.BaseInteraction):

    name = "injection"

    @classmethod
    def get_num_acts(cls, model: "model.TITAN", rel: "agent.Relationship") -> int:
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

        agent_params = (
            rel.agent1.location.params.demographics[rel.agent1.race]
            .sex_type[rel.agent1.sex_type]
            .injection
        )
        partner_params = (
            rel.agent2.location.params.demographics[rel.agent2.race]
            .sex_type[rel.agent2.sex_type]
            .injection
        )

        mean_num_acts = (
            min(agent_params.num_acts, partner_params.num_acts)
            * model.calibration.injection.act
        )
        share_acts = utils.poisson(mean_num_acts, model.np_random)

        if share_acts < 1:
            return 0

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
