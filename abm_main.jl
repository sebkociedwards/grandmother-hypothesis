using Random

Random.seed!(1)

include("src/types.jl")
include("src/genetics.jl")
include("src/demography.jl")
include("src/kinship.jl")
include("src/life_table.jl")
include("src/initialisation.jl")
include("src/simulation.jl")
include("src/plots.jl")

const FERTILITY_FILE = joinpath(@__DIR__, "data", "raw", "fertility.txt")
const MORTALITY_FILE = joinpath(@__DIR__, "data", "raw", "mortality_both.txt")
const OUTPUT_DIR = joinpath(@__DIR__, "data", "processed")

config = RunConfig(
    "SWE", 
    2024,
    true, # scale mx; false = raw max, true = scaled mx with R0 = 1
    FERTILITY_FILE, MORTALITY_FILE, OUTPUT_DIR
)

params = SimParams(
    0.5, # female fraction; use 0.5 if mx is total births and you want daughters
    1, # age bin width; 1 keeps original single-year resolution; e.g. 5 gives 0-4, 5-9 ...
    10_000, # n0
    2.0, # burnin T
    1, # total years
    15, # dependency age
    1.0, # grandmother multiplier; < 1.0 = benefit to child survival
    1.0, # helper cost multiplier; > 1.0 = cost to helping grandchild survival
    0.5, # help initial allele frequency
    0.5, # recog initial allele freqyency
    false, # use custom life table
)

life_table = compile_life_table(config, params)

# TODO needs to be debugged
if params.use_custom_life_table
    age_start = collect(0:params.age_bin_width:110)
    qx = make_type2_qx(age_start; qx_const=0.1)
    mx = life_table.mx

    life_table = build_custom_life_table(age_start, mx, qx, config, params)
end


# initiate the simulation
burnin_record, sim_record = run(life_table, params)


println("Generating plots...")

fig = Figure(size=(600*4, 400*4))
Label(fig[0, 2:3], "Simulation test", fontsize=28, font=:bold, tellwidth=false)

# figure 1: demography input
# fig_demography = Figure(size=(600, 400))
plot_life_table(fig[1, 1], life_table, config, params)
# display(fig_demography)

# fig_pop_dynamics = Figure(size=(600, 1200))
plot_n_history(fig[1, 2], burnin_record, sim_record, params.age_bin_width)
plot_pf_history(fig[2, 2], burnin_record, sim_record, params.age_bin_width)
plot_population_age_distribution(fig[3, 2], sim_record.age_dist[1], sim_record.age_dist[end], life_table)
# display(fig_pop_dynamics)

# fig_grandmother_mechanism = Figure(size=(600, 800))
plot_dependent(fig[1, 3], burnin_record, sim_record, params.age_bin_width)
plot_p_helped_dependent(fig[2, 3], burnin_record, sim_record, params.age_bin_width)
# display(fig_grandmother_mechanism)

# fig_genetics = Figure(size=(600, 1500))
plot_frequency_history(fig[1, 4], burnin_record.help_af, sim_record.help_af, params.age_bin_width; ylabel="p", title="Help allele frequency")
plot_frequency_history(fig[2, 4], burnin_record.helper_pf, sim_record.helper_pf, params.age_bin_width; ylabel="frequency", title="Helper phenotype frequency")
plot_frequency_history(fig[3, 4], burnin_record.recog_af, sim_record.recog_af, params.age_bin_width; ylabel="p", title="Recognition allele frequency")
plot_frequency_history(fig[4, 4], burnin_record.recog_pf, sim_record.recog_pf, params.age_bin_width; ylabel="frequency", title="Recognition phenotype frequency")
# display(fig_genetics)

display(fig)


test_fig = Figure()
plot_life_table(test_fig[1, 1], life_table, config, params)
display(test_fig)

return