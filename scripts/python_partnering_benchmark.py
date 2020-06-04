import os
from sys import path

path_to_this_dir = os.path.dirname(os.path.realpath(__file__))
path.append(os.path.join(path_to_this_dir, ".."))

from titan.partnering import select_partner
from titan.agent import Agent, AgentSet
from titan.parse_params import ObjMap
from titan import utils

import random
import time as time_mod

run_rand = random.Random(123)

num_agents = 1000000

agents = set()
for i in range(num_agents):
    new_agent = Agent("so", 40, "r", "du")
    new_agent.partners["test"] = set()
    new_agent.race = utils.safe_random_choice(["a", "b", "c"], run_rand)
    agents.add(new_agent)

params = ObjMap(
    {
        "classes": {"bond_types": {"test": {"acts_allowed": ["sex", "injection"]}}},
        "features": {"assort_mix": True},
        "assort_mix": {
            "assort_def": {
                "attribute": "race",
                "agent_value": "a",
                "partner_values": {"a": 0.3, "__other__": 0.7},
            },
        },
    }
)

agent = Agent("so", 40, "r", "du")
agent.partners["test"] = set()
agent.race = "a"

pwid_agents = AgentSet("pwid")
pwid_agents.members = agents

num_trials = 100
wct = []

for i in range(num_trials):
    tic = time_mod.time()
    p = select_partner(
        agent, agents, {"so": agents}, pwid_agents, params, run_rand, "test"
    )
    wct.append(time_mod.time() - tic)


def mean(seq):
    return sum(seq) / len(seq)


print(("\nSUMMARY:\nall tasks - mean: %8.6f seconds" % mean(wct)))
print(("all tasks - min:  %8.6f seconds" % min(wct)))
print(("all tasks - max:  %8.6f seconds" % max(wct)))
print(("all tasks - sum:  %8.6f seconds" % sum(wct)))
