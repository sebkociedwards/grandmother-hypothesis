using LinearAlgebra

safe_age_idx(a, life_table) = min(a.age_idx, nrow(life_table))

mx_at(a, life_table) = life_table.mx[safe_age_idx(a, life_table)]
qx_at(a, life_table) = life_table.qx[safe_age_idx(a, life_table)]
age_years(a, life_table) = life_table.age_start[safe_age_idx(a, life_table)]

function years_to_age_idx(age_years, age_bin_width)
    return fld(age_years, age_bin_width) + 1
end

function calculate_R0(lx, mx, female_fraction)
    return sum(lx .* (female_fraction .* mx))
end

function calculate_T(life_table, female_fraction)
    lx = life_table.lx
    mx = life_table.mx
    ages_mid = life_table.age_start .+ ((life_table.age_end .- life_table.age_start) ./ 2)

    R0 = calculate_R0(lx, mx, female_fraction)
    return sum(ages_mid .* lx .* (female_fraction .* mx)) / R0
end

function stable_age_distribution(life_table, female_fraction)
    n = nrow(life_table)
    sx = clamp.(life_table.sx, 0.0, 1.0)

    L = zeros(Float64, n, n)
    L[1, :] .= female_fraction .* life_table.mx

    for i in 1:(n - 1)
        L[i + 1, i] = sx[i]
    end

    vals, vecs = eigen(L)
    idx = argmax(real(vals))
    w = abs.(real(vecs[:, idx]))
    w ./= sum(w)

    return w
end

function calculate_pf(population)
    isempty(population) && return NaN
    return count(agent -> agent.female, population) / length(population)
end

# return true if female mx <= 0
function is_postreproductive(agent, life_table)
    return agent.female && mx_at(agent, life_table) <= 0.0
end

function is_reproductive_male(a, age_bin_width; min_age=15, max_age=80)
    min_idx = years_to_age_idx(min_age, age_bin_width)
    max_idx = years_to_age_idx(max_age, age_bin_width)
    return !a.female && (min_idx <= a.age_idx <= max_idx)
end