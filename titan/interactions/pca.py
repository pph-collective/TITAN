import numpy as np  # type: ignore
import networkx as nx  # type: ignore

from . import base_interaction
from .. import utils


class PCA(base_interaction.BaseInteraction):

    name = "pca"

    @staticmethod
    def interact(model, rel) -> bool:
        """
        Simulate peer change agent interactions. Knowledge if one agent is aware and one unaware,
            opinion if one agent swaying the other.

        args:
            model: The running model
            rel: The relationship PCA is happening in
        """
        if not model.params.features.pca or model.time < model.params.pca.start_time:
            return False

        assert (
            model.params.model.network.enable
        ), "Network must be enabled for pca interactions"

        acts_prob = model.run_random.random()
        acts_bin = 0
        current_p_value = 0.0
        params = model.params.partnership.pca.frequency[rel.bond_type]

        while acts_prob > current_p_value:
            acts_bin += 1
            current_p_value += params[acts_bin].prob

        min = params[acts_bin].min
        max = params[acts_bin].max

        if min == max:
            num_acts = min
        else:
            num_acts = model.run_random.randrange(min, max)

        if num_acts < 1:
            return False

        agent1_aware = rel.agent1.pca.awareness  # type: ignore[attr-defined]
        agent2_aware = rel.agent2.pca.awareness  # type: ignore[attr-defined]

        if model.run_random.random() < knowledge_transmission_probability(model, rel):
            if agent1_aware and agent2_aware:
                influence(model, rel.agent1, rel.agent2)
            elif agent1_aware:
                knowledge_dissemination(model, rel.agent2)
            elif agent2_aware:
                knowledge_dissemination(model, rel.agent1)


# ===================== HELPER FUNCTIONS ===================
def influence(model, rel):
    """
    The higher influence agent chages the opinion of the other agent to be the mean of their opinions.  If the opinion excedes the threshold, initiate prep.

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

    partner_init_opinion = partner.pca.opinion
    partner.pca.opinion = np.mean([agent.pca.opinion, partner.pca.opinion])

    # partner crossed threshhold of pca opinion, stochastically enroll in prep
    if (
        partner_init_opinion
        < partner.location.params.pca.opinion.threshold
        < partner.pca.opinion
    ):
        if model.run_random.random() < partner.location.params.pca.prep.prob:
            partner.prep.initiate(model, force=True)


def knowledge_dissemination(model, partner):
    """
    Make an agent pca aware, stochastically enroll in prep if their opinion meets the threshold

    args:
        model: The running model
        partner: The agent to whom knowledge is being disseminated
    """
    partner.pca.awareness = True
    if (
        partner.pca.opinion > partner.location.params.pca.opinion.threshold
        and model.run_random.random() < partner.location.params.pca.prep.prob
    ):
        partner.prep.initiate(model, force=True)


def knowledge_transmission_probability(model, rel):
    """
    Get the probability of knowledge/opinion transmission in this relationship

    args:
        model: The running model
        rel: The relationship where the pca interaction is happening
    """
    if rel.agent1.pca.awareness and rel.agent2.pca.awareness:
        p = model.params.pca.opinion.transmission
    else:
        p = model.params.pca.knowledge.transmission

    if num_acts == 1:
        p_total_transmission = p
    elif num_acts >= 1:
        p_total_transmission = 1.0 - utils.binom_0(num_acts, p)
    else:
        p_total_transmission = 0

    return p_total_transmission
