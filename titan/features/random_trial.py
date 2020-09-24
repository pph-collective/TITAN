from . import base_feature
from .. import utils
from .. import agent
from .. import population
from .. import model

import networkx as nx  # type: ignore


class RandomTrial(base_feature.BaseFeature):

    name = "random_trial"

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False
        self.treated = False

    @classmethod
    def update_pop(cls, model: "model.HIVModel"):
        """
        Update the feature for the entire population (class method).

        Initialize a random trial in the population if time is the random trial start time.

        args:
            model: the instance of HIVModel currently being run
        """
        if (
            model.params.features.prep
            and model.time == model.params.random_trial.start_time
        ):
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
                        centrality = nx.algorithms.centrality.eigenvector_centrality(
                            comp
                        )
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

                        chosen_agent = utils.safe_random_choice(
                            comp_agents, model.run_random
                        )  # select change agent
                        if chosen_agent is not None:
                            chosen_agent.pca.suitable = True  # type: ignore[attr-defined]
                            chosen_agent.random_trial.treated = True  # type: ignore[attr-defined]

                        else:  # if no suitable agents, mark a non-suitable agent
                            chosen_agent = utils.safe_random_choice(
                                list(comp.nodes), model.run_random
                            )

                        chosen_agent.pca.awareness = True  # type: ignore[attr-defined]
                        chosen_agent.pca.active = True  # type: ignore[attr-defined]

                    # pca - random
                    elif model.params.pca.choice == "random":
                        suitable_agent_choices = [
                            agent for agent in comp.nodes if not agent.hiv
                        ]

                        chosen_agent = utils.safe_random_choice(
                            suitable_agent_choices, model.run_random
                        )

                        if (
                            chosen_agent is not None
                        ):  # if there are agents who meet eligibility criteria,
                            # select one randomly

                            chosen_agent.pca.suitable = True
                            chosen_agent.random_trial.treated = True
                        else:  # if no suitable agents, mark a non-suitable agent
                            chosen_agent = utils.safe_random_choice(
                                list(comp.nodes), model.run_random
                            )

                        chosen_agent.pca.awareness = True  # make aware
                        chosen_agent.pca.active = True

            print(("Total agents in trial: ", total_nodes))

    # ============= HELPER METHODS ====================
