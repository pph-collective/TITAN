from . import base_interaction
from .. import model
from .. import agent
from .. import utils


class PCA(base_interaction.BaseInteraction):

    name = "pca"

    @classmethod
    def get_num_acts(cls, model: "model.TITAN", rel: "agent.Relationship") -> int:
        params = model.params.partnership.pca.frequency[rel.bond_type]

        if params.type == "bins":
            params = params.bins
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
                num_acts = model.run_random.randint(min, max)
        elif params.type == "distribution":
            num_acts = round(utils.safe_dist(params.distribution, model.run_random))

        return num_acts
