using Random

include("utils.jl")
include("dotmap.jl")

mutable struct Agent
    id::Int
    so::String
    age::Int
    race::String
    drug_use::String
    partners::Dict{String,Set{Agent}}

    Agent() = new(rand(Int), "so", rand(20:80), "r", "du", Dict())
end

Agent(id, so, age, race, drug_use) = Agent(id, so, age, race, drug_use, Dict())


function assort!(eligible_partners::Set{Agent}, assort_params::DotMap, rand_gen::AbstractRNG)
    partner_values = collect(assort_params.partner_values)
    partner_types = first.(partner_values)
    partner_weights = last.(partner_values)
    partner_type = safe_random_choice(
        partner_types, rand_gen, partner_weights
    )
    attribute = Symbol(assort_params.attribute)

    if partner_type == "__other__"
        for p in partner_types
            if p != "__other__"
                filter!(partner -> string(getfield(partner, attribute)) != p, eligible_partners)
            end
        end
    else
        filter!(partner -> string(getfield(partner, attribute)) == partner_type, eligible_partners)
    end

    return eligible_partners
end

function select_partner(
    agent::Agent,
    partnerable_agents::Set{Agent},
    sex_partners::Dict{String,Set{Agent}},
    pwid_agents::Set{Agent},
    params::DotMap,
    rand_gen::AbstractRNG,
    bond_type::String,
)
    eligible = copy(partnerable_agents)
    setdiff!(eligible, Set([agent]))
    for bond in string.(keys(params.classes.bond_types))
        setdiff!(eligible, agent.partners[bond])
    end

    acts_allowed = params.classes.bond_types[bond_type].acts_allowed

    if "injection" in acts_allowed
        intersect!(eligible, pwid_agents)
    end

    if "sex" in acts_allowed
        intersect!(eligible, sex_partners[agent.so])
    end

    if !isempty(eligible) && params.features.assort_mix
        for assort_def in values(params.assort_mix)
            if string(getfield(agent, Symbol(assort_def.attribute))) == assort_def.agent_value
                assort!(eligible, assort_def, rand_gen)
            end
        end
    end

    return safe_random_choice(eligible, rand_gen)
end
