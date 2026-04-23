struct SimParams
    female_fraction::Float64
    age_bin_width::Int
    n0::Int
    burnin_T::Float64
    total_years::Int
    dep_age::Int
    grandmother_multiplier::Float64
    helper_cost_multiplier::Float64
    p_help_init::Float64
    p_recog_init::Float64
    use_custom_life_table::Bool
end

struct RunConfig
    country::String
    year::Int
    scale_mx::Bool
    fertility_file::String
    mortality_file::String
    output_dir::String
end

mutable struct Agent
    id::Int
    age_idx::Int
    female::Bool
    mother_id::Union{Int, Nothing}
    father_id::Union{Int, Nothing}
    help_alleles::NTuple{2, Bool}
    recog_alleles::NTuple{2, Bool}
end

mutable struct PopulationRecord
    n::Vector{Int}
    pf::Vector{Float64}
    help_af::Vector{Float64}
    helper_pf::Vector{Float64}
    recog_af::Vector{Float64}
    recog_pf::Vector{Float64}
    dep_total::Vector{Int}
    dep_helped::Vector{Int}
    age_dist::Vector{Vector{Int}}
end