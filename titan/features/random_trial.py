from typing import Dict
import logging

from . import base_feature
from .. import utils
from .. import model

import networkx as nx  # type: ignore


class RandomTrial(base_feature.BaseFeature):

    name = "random_trial"
    stats = [
        "random_trial",
        "random_trial_treated",
        "random_trial_suitable",
        "random_trial_treated_hiv",
    ]
    """
        Random Trial collects the following stats:

        * random_trial - number of agents with active random_trial
        * random_trial_treated - number of active agents treated
        * random_trial_treated_hiv - number of HIV+ agents treated
        * random_trial_suitable - number of active agents suitable
    """

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False
        self.treated = False
        self.suitable = False

    @classmethod
    def update_pop(cls, model: "model.TITAN"):
        """
        Update the feature for the entire population (class method).

        Initialize a random trial in the population if time is the random trial start time.

        args:
            model: the instance of TITAN currently being run
        """
        rt_params = model.params.random_trial

        if not model.time == rt_params.start_time:
            return

        assert (
            model.params.model.network.enable
        ), "Network must be enabled for random trial"

        logging.info(f"Starting random trial ({rt_params.choice})")
        components = model.pop.connected_components()

        # set up helper methods based on params
        if rt_params.treatment == "prep":
            assert (
                model.params.features.prep
            ), "Prep feature must be enabled to use the prep random trial treatment"
            treat = treat_prep
            suitable = suitable_prep
        elif rt_params.treatment == "knowledge":
            assert (
                model.params.exposures.knowledge
            ), "Knowledge exposure must be enabled to use the knowledge random trial treatment"
            treat = treat_knowledge
            suitable = suitable_knowledge

        total_nodes = 0
        logging.info(
            f"Number of components {len([1 for comp in components if comp.number_of_nodes()])}",
        )
        for comp in components:
            total_nodes += comp.number_of_nodes()
            if model.run_random.random() < rt_params.prob:
                # Component selected as treatment pod!
                for agent in comp.nodes:
                    agent.random_trial.active = True

                # treat all agents
                if rt_params.choice == "all":
                    for agent in comp.nodes():
                        if suitable(agent, model):
                            treat(agent, model)
                            agent.random_trial.suitable = True
                            agent.random_trial.treated = True

                # chose an agent central to the component
                elif rt_params.choice == "eigenvector":
                    centrality = nx.algorithms.centrality.eigenvector_centrality(comp)
                    assert len(centrality) >= 1, "Empty centrality"
                    ordered_centrality = sorted(centrality, key=centrality.get)

                    # find the most central suitable agent, or if none, use most central
                    intervention_agent = ordered_centrality[0]
                    for agent in ordered_centrality:
                        if suitable(agent, model):
                            intervention_agent = agent
                            intervention_agent.random_trial.suitable = True
                            break

                    intervention_agent.random_trial.treated = True
                    treat(intervention_agent, model)

                # chose an agent that is a bridge in the network
                elif rt_params.choice == "bridge":
                    # list all edges that are bridges
                    all_bridges = list(nx.bridges(comp))
                    suitable_agents = [
                        agent
                        for agents in all_bridges
                        for agent in agents
                        if suitable(agent, model)
                    ]  # all suitable agents in bridges

                    chosen_agent = utils.safe_random_choice(
                        suitable_agents, model.run_random
                    )  # select change agent
                    if chosen_agent is not None:
                        chosen_agent.random_trial.suitable = True  # type: ignore[attr-defined]

                    else:  # if no suitable agents, mark a non-suitable agent
                        chosen_agent = utils.safe_random_choice(
                            list(comp.nodes), model.run_random
                        )

                    chosen_agent.random_trial.treated = True  # type: ignore[attr-defined]
                    treat(chosen_agent, model)

                # chose an agent from the component at random
                elif rt_params.choice == "random":
                    suitable_agents = [
                        agent for agent in comp.nodes if suitable(agent, model)
                    ]

                    # if there are agents who meet eligibility criteria,
                    # select one randomly
                    chosen_agent = utils.safe_random_choice(
                        suitable_agents, model.run_random
                    )

                    if chosen_agent is not None:
                        chosen_agent.random_trial.suitable = True
                    else:  # if no suitable agents, mark a non-suitable agent
                        chosen_agent = utils.safe_random_choice(
                            list(comp.nodes), model.run_random
                        )

                    chosen_agent.random_trial.treated = True  # type: ignore[attr-defined]
                    treat(chosen_agent, model)

        logging.info(f"Total agents in trial: {total_nodes}")

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.active:
            stats["random_trial"] += 1
            if self.treated:
                stats["random_trial_treated"] += 1
                if self.agent.hiv.active:  # type: ignore[attr-defined]
                    stats["random_trial_treated_hiv"] += 1
            if self.suitable:
                stats["random_trial_suitable"] += 1

    # ============= HELPER METHODS ====================


# ============= HELPER FUNCTIONS ==================


def treat_prep(agent, model):
    agent.prep.enroll(model.run_random, model.time)


def suitable_prep(agent, model) -> bool:
    if (
        not agent.hiv.active
        and not agent.prep.active
        and model.run_random.random() < agent.location.params.prep.cap
        and not agent.vaccine.active
    ):
        return True
    else:
        return False


def treat_knowledge(agent, model):
    agent.knowledge.convert(model)


def suitable_knowledge(agent, model) -> bool:
    return not agent.hiv.active
