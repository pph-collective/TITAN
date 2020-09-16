from .base_feature import BaseFeature
from .. import utils

import networkx as nx # type: ignore


class RandomTrial(BaseFeature):

    name = "random_trial"

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False
        self.treated = False

    @classmethod
    def update_pop(cls, model):
        if (
            model.params.features.prep
            and model.time == model.params.random_trial.start_time
        ):
            cls.initialize_random_trial(model)

    @classmethod
    def initialize_random_trial(cls, model):
        """
        Initialize a random trial in the population
        """
        assert (
            model.params.model.network.enable
        ), "Network must be enabled for random trial"

        print("Starting random trial")
        components = model.pop.connected_components()

        total_nodes = 0
        print(
            "Number of components",
            len([1 for comp in components if comp.number_of_nodes()]),
        )
        for comp in components:
            total_nodes += comp.number_of_nodes()
            if model.run_random.random() < model.params.random_trial.prob:
                # Component selected as treatment pod!
                for ag in comp.nodes:
                    ag.random_trial.active = True

                # prep case
                if not model.params.features.pca:
                    for ag in comp.nodes():
                        if not ag.hiv and not ag.prep.active:
                            ag.random_trial.treated = True
                            if (
                                model.run_random.random()
                                < ag.location.params.prep.target
                                and not ag.vaccine.active
                            ):
                                ag.prep.enroll(ag, model.run_random)

                # pca - eigenvector
                elif model.params.pca.choice == "eigenvector":
                    centrality = nx.algorithms.centrality.eigenvector_centrality(comp)
                    assert len(centrality) >= 1, "Empty centrality"
                    ordered_centrality = sorted(centrality, key=centrality.get)

                    # find the most central suitable agent, or if none, use most central
                    intervention_agent = ordered_centrality[0]
                    for ag in ordered_centrality:
                        if not ag.hiv:
                            intervention_agent = ag
                            intervention_agent.sutable = True
                            intervention_agent.random_trial.treated = True
                            break

                    intervention_agent.pca.awareness = True
                    intervention_agent.pca.active = True

                # pca - bridge
                elif model.params.pca.choice == "bridge":
                    # list all edges that are bridges
                    all_bridges = list(nx.bridges(comp))
                    comp_agents = [
                        agent
                        for agents in all_bridges
                        for agent in agents
                        if not agent.hiv
                    ]  # all suitable agents in bridges

                    if comp_agents:
                        chosen_agent = utils.safe_random_choice(
                            comp_agents, model.run_random
                        )  # select change agent
                        chosen_agent.pca.suitable = True
                        chosen_agent.random_trial.treated = True

                    else:  # if no suitable agents, mark a non-suitable agent
                        chosen_agent = utils.safe_random_choice(
                            list(comp.nodes), model.run_random
                        )

                    chosen_agent.pca.awareness = True  # make aware
                    chosen_agent.pca.active = True

                # pca - random
                elif model.params.pca.choice == "random":
                    suitable_agent_choices = [
                        agent for agent in comp.nodes if not agent.hiv
                    ]

                    if (
                        suitable_agent_choices
                    ):  # if there are agents who meet eligibility criteria,
                        # select one randomly
                        chosen_agent = utils.safe_random_choice(
                            suitable_agent_choices, model.run_random
                        )
                        chosen_agent.pca.suitable = True
                        chosen_agent.random_trial.treated = True
                    else:  # if no suitable agents, mark a non-suitable agent
                        chosen_agent = utils.safe_random_choice(
                            list(comp.nodes), model.run_random
                        )

                    chosen_agent.pca.awareness = True  # make aware
                    chosen_agent.pca.active = True

        print(("Total agents in trial: ", total_nodes))
