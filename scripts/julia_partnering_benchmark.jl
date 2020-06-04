include(joinpath(@__DIR__, "..", "titan", "partnering.jl"))

using ProfileView
using Dates
using Statistics

import Profile

Profile.init(n = 10^7, delay = 0.005)

num_agents = 1000000

agents = Set(Agent[])
for i = 1:num_agents
    new_agent = Agent()
    new_agent.partners["test"] = Set(Agent[])
    new_agent.race = rand(["a", "b", "c"])
    push!(agents, new_agent)
end

run_rand = MersenneTwister(123)

params = DotMap(Dict(
    "classes" => Dict(
        "bond_types" => Dict(
            "test" => Dict("acts_allowed" => ["sex", "injection"]),
        ),
    ),
    "features" => Dict("assort_mix" => true),
    "assort_mix" => Dict(
        "assort_def" => Dict(
            "attribute" => "race",
            "agent_value" => "a",
            "partner_values" => Dict("a" => 0.3, "__other__" => 0.7),
        ),
    ),
))

agent = Agent()
agent.partners["test"] = Set(Agent[])
agent.race = "a"

function profile_test(num_trials, agents, run_rand, params, agent)
    wct = []

    for i = 1:num_trials
        tic = now()
        p = select_partner(
            agent,
            agents,
            Dict("so" => agents),
            agents,
            params,
            run_rand,
            "test",
        )
        push!(wct, Dates.toms(now() - tic))
    end

    println("\nSUMMARY:")
    println("all tasks - mean: $(mean(wct) / 1000) seconds")
    println("all tasks - min:  $(minimum(wct) / 1000) seconds")
    println("all tasks - max:  $(maximum(wct) / 1000) seconds")
    println("all tasks - sum:  $(sum(wct) / 1000) seconds")
end

profile_test(100, agents, run_rand, params, agent)

# @profview profile_test(1, agents, run_rand, params, agent)
# @profview profile_test(100, agents, run_rand, params, agent)
