from . import base_feature
from .. import agent
from .. import model


class PartnerTracing(base_feature.BaseFeature):

    name = "partner_tracing"

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False
        self.time = None

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        If the agent is was diagnosed last time step, trace their partners. If the agent is traced but not diagnosed, stochastically diagnose.  If the agents tracking has expired, mark them as inactive.

        args:
            model: the instance of TITAN currently being run
        """
        params = self.agent.location.params.partner_tracing

        if model.time < params.start_time or model.time > params.stop_time:
            return

        agent_exposure = getattr(self.agent, params.exposure)

        if agent_exposure.active:
            # was this agent diagnosed with the target exposure last time step
            if agent_exposure.dx and agent_exposure.dx_time == model.time - 1:
                for ptnr in self.agent.get_partners(params.bond_type):
                    partner_exposure = getattr(ptnr, params.exposure)
                    if (
                        not partner_exposure.dx
                        and model.run_random.random() < params.prob
                    ):
                        ptnr.partner_tracing.active = True  # type: ignore[attr-defined]
                        ptnr.partner_tracing.time = model.time  # type: ignore[attr-defined]

            # second chance at diagnosis if traced
            if (
                self.active
                and self.time < model.time
                and not agent_exposure.dx
                and model.run_random.random() < params.dx_prob
            ):
                agent_exposure.diagnose(model)

        # stop tracing of this agent if time
        if self.active and model.time >= self.time + params.trace_duration:
            self.active = False
            self.time = None
