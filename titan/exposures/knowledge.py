from typing import List, Dict

import numpy as np  # type: ignore
import networkx as nx  # type: ignore

from . import base_exposure
from .. import agent as ag
from .. import population
from .. import model
from .. import utils


class Knowledge(base_exposure.BaseExposure):

    name: str = "knowledge"
    """Name of exposure in the params file.  Also used to name the attribute in Agent"""

    stats: List[str] = ["knowledge_aware"]
    """List of names of stats that come from this exposure (e.g. hiv.dx)"""

    def __init__(self, agent: "ag.Agent"):
        super().__init__(agent)

        self.active = False
        self.opinion = 0.0

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        knowledge_params = self.agent.location.params.knowledge
        if pop.pop_random.random() < knowledge_params.init:
            self.active = True

            attprob = pop.pop_random.random()
            pvalue = 0.0
            for bin, fields in knowledge_params.opinion.init.items():
                pvalue += fields.prob
                if attprob < pvalue:
                    self.opinion = bin
                break

    def update_agent(self, model: "model.HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        args:
            model: the instance of HIVModel currently being run
        """
        knowledge_params = self.agent.location.params.knowledge
        if (
            model.time >= knowledge_params.start_time
            and model.run_random.random() < knowledge_params.prob
        ):
            self.convert(model)

    def set_stats(self, stats: Dict[str, int], time: int):
        """
        Update the `stats` dictionary passed for this agent.  Called from `output.get_stats` for each enabled feature in the model.

        The stats to be updated must be declared in the class attribute `stats` to make sure the dictionary has the expected keys/counter value initialized.

        args:
            stats: the dictionary to update with this agent's feature statistics
            time: the time step of the model when the stats are set
        """
        if self.active:
            stats["knowledge_aware"] += 1

    @staticmethod
    def expose(
        model: "model.HIVModel",
        interaction: str,
        rel: "ag.Relationship",
        num_acts: int,
    ):
        assert (
            model.params.model.network.enable
        ), "Network must be enabled for knowledge exposure"

        # agent/partner ordering is irrelevant at this point for knowledge transmission
        p = rel.agent1.knowledge.get_transmission_probability(  # type: ignore[attr-defined]
            model, interaction, rel.agent2, num_acts
        )

        if model.run_random.random() < p:
            agent1_aware = rel.agent1.knowledge.active  # type: ignore[attr-defined]
            agent2_aware = rel.agent2.knowledge.active  # type: ignore[attr-defined]

            if agent1_aware and agent2_aware:
                influence(model, rel)
            elif agent1_aware:
                rel.agent2.knowledge.convert(model)  # type: ignore[attr-defined]
            elif agent2_aware:
                rel.agent1.knowledge.convert(model)  # type: ignore[attr-defined]

    def get_transmission_probability(
        self,
        model: "model.HIVModel",
        interaction: str,
        partner: "ag.Agent",
        num_acts: int,
    ) -> float:
        """
        Get the probability of knowledge/opinion transmission in this relationship

        args:
            model: The running model
            rel: The relationship where the knowledge interaction is happening
            num_acts: The number of interactions the agents had

        returns:
            the probability of knowledge/opinion transmission
        """
        if not interaction == "pca":
            return 0.0

        if self.active and partner.knowledge.active:  # type: ignore[attr-defined]
            p = model.params.knowledge.opinion.prob
        else:
            p = model.params.knowledge.prob

        return utils.total_probability(p, num_acts)

    def convert(self, model: "model.HIVModel"):
        """
        Make an agent aware, stochastically enroll in feature if their opinion meets the threshold

        args:
            model: The running model
        """
        params = self.agent.location.params.knowledge
        self.active = True  # type: ignore[attr-defined]
        if (
            self.opinion > params.opinion.threshold  # type: ignore[attr-defined]
            and model.run_random.random() < params.feature.prob
        ):
            agent_attr = getattr(self.agent, params.feature.name)
            agent_attr.initiate(model, force=True)


# ===================== HELPER FUNCTIONS ===================
def influence(model: "model.HIVModel", rel: "ag.Relationship"):
    """
    The higher influence agent chages the opinion of the other agent to be the mean of their opinions.  If the opinion excedes the threshold, initiate feature.

    args:
        model: The running model
        rel: a relationship where an agent is influencing their partner
    """

    def get_influence(agent):
        return nx.closeness_centrality(model.pop.graph, agent)

    if get_influence(rel.agent1) > get_influence(rel.agent2):
        agent = rel.agent1
        partner = rel.agent2
    else:
        agent = rel.agent2
        partner = rel.agent1

    partner_init_opinion = partner.knowledge.opinion  # type: ignore[attr-defined]
    partner.knowledge.opinion = np.mean([agent.knowledge.opinion, partner.knowledge.opinion])  # type: ignore[attr-defined]

    # partner crossed threshhold of opinion, stochastically enroll in feature
    params = partner.location.params.knowledge
    if (
        partner_init_opinion
        < params.opinion.threshold
        < partner.knowledge.opinion  # type: ignore[attr-defined]
    ):
        if model.run_random.random() < params.feature.prob:
            agent_attr = getattr(partner, params.feature.name)
            agent_attr.initiate(model, force=True)
