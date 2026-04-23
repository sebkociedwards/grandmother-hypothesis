function rand_allele(p)
    return rand() < p
end

function inherit_allele(alleles)
    return rand(alleles)
end


# genotype to phenotype logic

function is_helper(agent)
    return agent.help_alleles[1] || agent.help_alleles[2]
end

function has_recognition(agent)
    return agent.recog_alleles[1] || agent.recog_alleles[2]
end


# measure genetic data

function help_allele_frequency(population)
    isempty(population) && return NaN
    total = 2 * length(population)
    n_help = sum(Int(a.help_alleles[1]) + Int(a.help_alleles[2]) for a in population)
    return n_help / total
end

function helper_frequency(population)
    isempty(population) && return NaN
    return count(is_helper, population) / length(population)
end

function recog_allele_frequency(population)
    isempty(population) && return NaN
    total = 2 * length(population)
    n_recog = sum(Int(a.recog_alleles[1]) + Int(a.recog_alleles[2]) for a in population)
    return n_recog / total
end

function recognition_frequency(population)
    isempty(population) && return NaN
    return count(has_recognition, population) / length(population)
end