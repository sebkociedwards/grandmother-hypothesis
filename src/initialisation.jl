using StatsBase

# entry point to initlaise founding population
function init_population(life_table, params::SimParams)
    population = Agent[]

    w = stable_age_distribution(life_table, params.female_fraction)
    age_idxs = life_table.age_idx
    
    # initiate population
    for i in 1:params.n0
        age_idx = sample(age_idxs, Weights(w))

        help_alleles = (rand_allele(params.p_help_init), rand_allele(params.p_help_init))
        recog_alleles = (rand_allele(params.p_recog_init), rand_allele(params.p_recog_init))

        agent = Agent(i, age_idx, rand(Bool), nothing, nothing, help_alleles, recog_alleles)
        push!(population, agent)
    end

    next_id = params.n0 + 1
    return population, next_id
end