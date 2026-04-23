using DelimitedFiles
using DataFrames
using CSV

# read data from fertility/mortality database, assuming the first 3 rows are metadata
function read_raw_table(file)
    lines = readlines(file)
    return [split(line) for line in lines[4:end] if !isempty(strip(line))]
end

# filter raw table rows matching "country" and "year", assuming that column 1 is the country code, and column 2 is the year
function filter_raw_rows(rows, country, year) 
    year_str = string(year)
    return [row for row in rows if row[1] == country && row[2] == year_str]
end

# convert the age lables like "12", "12-", or "110+" into an integer age key
function parse_age(age_label)
    cleaned = strip(age_label)
    numeric = replace(cleaned, "+" => "", "-" => "")
    return parse(Int, numeric)
end

function build_age_value_dict(rows, specs::NamedTuple)
    values = Dict{Int, NamedTuple}()

    for row in rows
        age = parse_age(row[3])
        values[age] = (; (name => parse(Float64, row[col]) for (name, col) in pairs(specs))...)
    end
    return values
end

function merge_fertility_mortality_rows(fertility_rows, mortality_rows)
    fert = build_age_value_dict(fertility_rows, (mx = 4, ))
    mort = build_age_value_dict(mortality_rows, (qx = 5, lx = 7))

    ages = sort(collect(keys(mort)))

    df = DataFrame(
        age_start = ages, # age start is the lower bound of each age interval
        mx = [get(fert, a, (mx = 0.0, )).mx for a in ages],
        qx = [mort[a].qx for a in ages],
        lx = [mort[a].lx for a in ages]
    )

    # normalise lx so lx0 = 1
    df.lx ./= first(df.lx)
    
    return df
end

function compute_lx_from_qx(qx::Vector{Float64})
    lx = zeros(Float64, length(qx))
    lx[1] = 1.0
    for i in 2:length(qx)
        lx[i] = lx[i-1] * (1.0 - qx[i-1])
    end
    return lx
end

function rebin_life_table(df, age_bin_width)
    df = copy(df)
    df.age_idx = fld.(df.age_start, age_bin_width)

    rebinned = combine(groupby(df, :age_idx),
        :age_start => minimum => :age_start,
        :mx => sum => :mx,
        :qx => (qx -> 1 .- prod(1 .- qx)) => :qx,
        :lx => first => :lx
    )

    sort!(rebinned, :age_start)
    return rebinned
end

function finalise_life_table!(df, age_bin_width)
    sort!(df, :age_start)
    df.age_idx = collect(1:nrow(df))
    df.age_end = df.age_start .+ age_bin_width .- 1
    df.sx = 1 .- df.qx

    select!(df, [:age_idx, :age_start, :age_end, :lx, :qx, :sx, :mx])

    return df
end

function scale_mx!(df, female_fraction)
    R0 = sum(df.lx .* (female_fraction .* df.mx))
    R0 == 0.0 && error("R0 is zero; cannot scale mx")
    df.mx ./= R0
    return df
end

function make_type2_qx(age_start; qx_const=0.05)
    return fill(qx_const, length(age_start))
end

function build_custom_life_table(age_start::Vector{Int}, mx::Vector{Float64}, qx::Vector{Float64}, config::RunConfig, params::SimParams)
    length(age_start) == length(mx) == length(qx) || error("age_start, mx, qx must have same length")

    df = DataFrame(
        age_start = age_start,
        mx = mx,
        qx = qx,
        lx = compute_lx_from_qx(qx)
    )

    finalise_life_table!(df, params.age_bin_width)

    if config.scale_mx
        scale_mx!(df, params.female_fraction)
    end

    return df
end

function compile_life_table(config::RunConfig, params::SimParams)
    fertility_rows = filter_raw_rows(read_raw_table(config.fertility_file), config.country, config.year)
    mortality_rows = filter_raw_rows(read_raw_table(config.mortality_file), config.country, config.year)

    isempty(fertility_rows) && error("No fertility rows found for $(config.country) $(config.year)")
    isempty(mortality_rows) && error("No mortality rows found for $(config.country) $(config.year)")

    df = merge_fertility_mortality_rows(fertility_rows, mortality_rows)

    if params.age_bin_width > 1
        df = rebin_life_table(df, params.age_bin_width)
    end

    finalise_life_table!(df, params.age_bin_width)

    if config.scale_mx
        scale_mx!(df, params.female_fraction)
    end

    mkpath(config.output_dir)
    outname = "$(config.country)_$(config.year)_bin$(params.age_bin_width)_$(config.scale_mx ? "scaled" : "raw").csv"
    CSV.write(joinpath(config.output_dir, outname), df)

    return df
end