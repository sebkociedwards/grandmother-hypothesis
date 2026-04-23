using Distributions

function init_population_record()
    return PopulationRecord(
        Int[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Int[],
        Int[],
        Vector{Vector{Int}}()
    )
end


# helper functions

function is_dependent(a, params) 
    dep_age_idx = floor(Int, params.dep_age / params.age_bin_width)
    return a.age_idx < dep_age_idx
end

function age_distribution_counts(population, n_bins)
    counts = zeros(Int, n_bins)
    for agent in population
        idx = min(agent.age_idx, n_bins)
        counts[idx] += 1
    end
    return counts
end

function total_dependent(population, params::SimParams)
    return count(agent -> is_dependent(agent, params), population)
end

function total_helped(population, life_table, params::SimParams)
    id_map = build_id_map(population)
    return count(agent ->
        is_dependent(agent, params) &&
        has_active_helping_maternal_grandmother(agent, id_map, life_table),
        population
    )
end

function helping_grandmother_ids(population, life_table, params::SimParams)
    id_map = build_id_map(population)
    helpers = Set{Int}()

    for a in population
        if is_dependent(a, params)
            gm = get_maternal_grandmother(a, id_map)
            if gm !== nothing &&
               is_postreproductive(gm, life_table) &&
               is_helper(gm) &&
               has_recognition(gm)

                push!(helpers, gm.id)
            end
        end
    end

    return helpers
end   


# core functions

function record_population!(rec::PopulationRecord, population, life_table, params::SimParams)
    push!(rec.n, length(population))
    push!(rec.pf, calculate_pf(population))
    push!(rec.help_af, help_allele_frequency(population))
    push!(rec.helper_pf, helper_frequency(population))
    push!(rec.recog_af, recog_allele_frequency(population))
    push!(rec.recog_pf, recognition_frequency(population))

    dep_total = total_dependent(population, params)
    dep_helped = total_helped(population, life_table, params)

    push!(rec.dep_total, dep_total)
    push!(rec.dep_helped, dep_helped)
    push!(rec.age_dist, age_distribution_counts(population, nrow(life_table)))

    return rec
end

function fertility_step!(life_table, population, next_id, params::SimParams)
    newborns = Agent[]
    females = [agent for agent in population if agent.female && mx_at(agent, life_table) > 0.0]
    males = [agent for agent in population if is_reproductive_male(agent, params.age_bin_width)]

    (isempty(females) || isempty(males)) && return newborns, next_id

    for mother in females
        rate = mx_at(mother, life_table)
        rate <= 0.0 && continue

        n_births = rand(Poisson(rate))
        for _ in 1:n_births
            father = rand(males)

            child_help = (
                inherit_allele(mother.help_alleles),
                inherit_allele(father.help_alleles)
            )
            child_recog = (
                inherit_allele(mother.recog_alleles),
                inherit_allele(father.recog_alleles)
            )

            newborn = Agent(next_id, 1, rand(Bool), mother.id, father.id, child_help, child_recog)
            push!(newborns, newborn)
            next_id += 1
        end
    end

    return newborns, next_id
end

function mortality_step!(life_table, population, params::SimParams)
    survivors = Agent[]
    id_map = build_id_map(population)
    active_helper_ids = helping_grandmother_ids(population, life_table, params)

    for agent in population
        pr = qx_at(agent, life_table)

        if is_dependent(agent, params) &&
           has_active_helping_maternal_grandmother(agent, id_map, life_table)
            pr *= params.grandmother_multiplier
        end

        if agent.id in active_helper_ids
            pr *= params.helper_cost_multiplier
        end

        pr = clamp(pr, 0.0, 1.0)
        rand() > pr && push!(survivors, agent)
    end

    return survivors
end

function age_step!(population)
    for agent in population
        agent.age_idx += 1
    end
end

# to run every iteration of the simulation
# order or oporations:
# 1. babies born -> age and chance for death next step
# 2. death of adults
# 3. age adults
function step!(life_table, population, next_id, params::SimParams)
    newborns, next_id = fertility_step!(life_table, population, next_id, params)
    append!(population, newborns)

    survivors = mortality_step!(life_table, population, params)
    age_step!(survivors)

    empty!(population)
    append!(population, survivors)

    return next_id
end

function run_from_population(life_table, population, next_id, steps, params::SimParams)
    rec = init_population_record()
    record_population!(rec, population, life_table, params)

    for i in 1:steps
        next_id = step!(life_table, population, next_id, params)
        record_population!(rec, population, life_table, params)

        pct = round(100 * i / steps, digits=1)
        print("\rProgress: $pct%")
        flush(stdout)
    end

    println("")
    return rec, population, next_id
end

function run(life_table, params::SimParams)
    println("Initialising population...")
    population, next_id = init_population(life_table, params)

    R0 = calculate_R0(life_table.lx, life_table.mx, params.female_fraction)
    T = calculate_T(life_table, params.female_fraction)
    println("Initialised population: n0 = $(params.n0); R0 = $(round(R0, digits=2)); T = $(round(T, digits=2))")

    burnin_steps = ceil(Int, params.burnin_T * T / params.age_bin_width)
    println("Initiating burnin: $burnin_steps steps; $(burnin_steps * params.age_bin_width) years...")
    burnin_record, population, next_id = run_from_population(life_table, population, next_id, burnin_steps, params)

    sim_steps = ceil(Int, params.total_years / params.age_bin_width)
    println("Initiating simulation: $sim_steps steps; $(params.total_years) years...")
    sim_record, population, next_id = run_from_population(life_table, population, next_id, sim_steps, params)

    return burnin_record, sim_record
end