using Random
using StatsBase

function safe_random_choice(seq, rand_gen::MersenneTwister)
    if !isempty(seq)
        return rand(rand_gen, seq)
    end
end

function safe_random_choice(seq, rand_gen, weights)
    if !isempty(seq)
        return sample(rand_gen, seq, Weights(weights))
    end
end

safe_random_choice(seq::Set, rand_gen, weights) = safe_random_choice(collect(seq), rand_gen, weights)
