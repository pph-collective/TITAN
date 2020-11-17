from . import base_interaction
from .. import model
from .. import agent


class PCA(base_interaction.BaseInteraction):

    name = "pca"

    @staticmethod
    def get_num_acts(model: "model.TITAN", rel: "agent.Relationship") -> int:
        params = model.params.partnership.pca.frequency[rel.bond_type]
        acts_prob = model.run_random.random()
        acts_bin = 1
        current_p_value = params[acts_bin].prob

        while acts_prob > current_p_value:
            acts_bin += 1
            current_p_value += params[acts_bin].prob

        min = params[acts_bin].min
        max = params[acts_bin].max

        if min == max:
            num_acts = min
        else:
            num_acts = model.run_random.randrange(min, max)

        return num_acts
