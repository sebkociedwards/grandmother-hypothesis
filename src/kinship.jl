# assign population to a dictionary for kin lookup; id => agent
function build_id_map(population)
    return Dict(agent.id => agent for agent in population)
end

function get_mother(agent, id_map)
    isnothing(agent) && return nothing
    isnothing(agent.mother_id) && return nothing
    return get(id_map, agent.mother_id, nothing)
end


function get_maternal_grandmother(agent, id_map)
    mother = get_mother(agent, id_map)
    return get_mother(mother, id_map)
end

function has_active_helping_maternal_grandmother(agent, id_map, life_table)
    gm = get_maternal_grandmother(agent, id_map)
    gm === nothing && return false
    return is_postreproductive(gm, life_table) &&
        is_helper(gm) &&
        has_recognition(gm)
end









