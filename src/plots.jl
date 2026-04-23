using CairoMakie


# figure 1

function plot_life_table(fig, df, config::RunConfig, params::SimParams; series=[:lx, :qx, :mx])
    ax = Axis(
        fig, 
        xlabel="Age",
        ylabel="Value",
        title="$(config.country) $(config.year) bin$(params.age_bin_width) $(config.scale_mx ? "scaled" : "raw") life table"
    )

    for param in series
        lines!(ax, df.age_start, df[!, param], label=String(param))
    end

    axislegend(ax)
    return ax 
end


function history_x(record, age_bin_width; burnin=false)
    if burnin
        return (-(length(record.n)-1):0) .* age_bin_width
    else
        return (0:length(record.n)-1) .* age_bin_width
    end
end


# figure 2

function plot_n_history(fig, burnin_record, sim_record, age_bin_width)
    ax = Axis(fig, xlabel="Year", ylabel="N", title="N by year", ytickformat="{:,.0f}")
    lines!(ax, history_x(burnin_record, age_bin_width; burnin=true), burnin_record.n, color=:black, label="burnin")
    lines!(ax, history_x(sim_record, age_bin_width), sim_record.n, color=:seagreen4, label="sim")
    return ax
end

function plot_pf_history(fig, burnin_record, sim_record, age_bin_width)
    ax = Axis(fig, xlabel="Year", ylabel="pf", title="pf by year")
    ylims!(ax, 0.45, 0.55)
    lines!(ax, history_x(burnin_record, age_bin_width; burnin=true), burnin_record.pf, color=:black, label="burnin")
    lines!(ax, history_x(sim_record, age_bin_width), sim_record.pf, color=:seagreen4, label="sim")
    return ax
end

function plot_population_age_distribution(fig, start_counts, end_counts, life_table)
    x = life_table.age_start

    start_props = start_counts ./ sum(start_counts)
    end_props = end_counts ./ sum(end_counts)

    width = first(life_table.age_end .- life_table.age_start) + 1

    ax = Axis(fig, xlabel="Age", ylabel="Proportion", title="Population age distribution (standardised)")
    barplot!(ax, x, start_props; width=width, gap=0, color=(:red, 0.5), label="start")
    barplot!(ax, x, end_props; width=width, gap=0, color=(:blue, 0.5), label="end")
    axislegend(ax)
    return ax
end


# figure 3

function plot_dependent(fig, burnin_record, sim_record, age_bin_width)
    ax = Axis(fig, xlabel="Year", ylabel="N", title="Dependent individuals per year")
    lines!(ax, history_x(burnin_record, age_bin_width; burnin=true), burnin_record.dep_total, color=:black)
    lines!(ax, history_x(burnin_record, age_bin_width; burnin=true), burnin_record.dep_helped, color=:black)
    lines!(ax, history_x(sim_record, age_bin_width), sim_record.dep_total, color=:seagreen4, label="total")
    lines!(ax, history_x(sim_record, age_bin_width), sim_record.dep_helped, color=:sienna2, label="helped")
    axislegend(ax)
    return ax
end

function plot_p_helped_dependent(fig, burnin_record, sim_record, age_bin_width)
    burnin_p = [t == 0 ? 0.0 : h / t for (t, h) in zip(burnin_record.dep_total, burnin_record.dep_helped)]
    sim_p = [t == 0 ? 0.0 : h / t for (t, h) in zip(sim_record.dep_total, sim_record.dep_helped)]

    ax = Axis(fig, xlabel="Year", ylabel="p", title="Proportion of dependent individuals helped")
    ylims!(ax, 0, 1)

    lines!(ax, history_x(burnin_record, age_bin_width; burnin=true), burnin_p, color=:black, label="burnin")
    lines!(ax, history_x(sim_record, age_bin_width), sim_p, color=:seagreen4, label="sim")

    return ax
end


# figure 4

function plot_frequency_history(fig, burnin_history, sim_history, age_bin_width; ylabel, title)
    ax = Axis(fig, xlabel="Year", ylabel=ylabel, title=title)
    ylims!(ax, 0, 1)

    burnin_x = (-(length(burnin_history)-1):0) .* age_bin_width
    sim_x = (0:length(sim_history)-1) .* age_bin_width

    lines!(ax, burnin_x, burnin_history, label="burnin", color=:black)
    lines!(ax, sim_x, sim_history, label="sim", color=:seagreen4)

    return ax
end
