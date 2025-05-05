

function get_ecosystem(
    local_niche::AbstractVector{Int},
    dist::AbstractVector{Int},
    cell_labels::Dict,
    base_cells::BaseCellType)::AbstractMatrix
    ecosystem =  zeros(length(base_cells.types), base_cells.wave_diffusion)
    for cell in eachindex(local_niche)
        depth_index = dist[cell]
        depth_index == 0 && continue
        cell_index = cell_labels[local_niche[cell]]
        ecosystem[cell_index, depth_index] += 1
    end
    return ecosystem
end

function apply_shift(
    shifts::AbstractVector{Int},
    noise::Float64,
    damp::AbstractVector{Float64})
    for i in eachindex(shifts)
        local_shift = 0
        for d in damp
            local_shift += apply_noise(shifts[i] * d, noise, 1)
        end
            new_shifts[i] = sum(local_shift)
        end
        
    return new_shifts
end

function apply_type_shift(
    type_shifts::Dict{K, Vector{T}},
    noise::Float64) where {K, T}
    new_shifts = Dict(k => Vector{Int}() for k in keys(type_shifts))
    for k in keys(type_shifts)
        new_shifts[k] = apply_shift(type_shifts, noise, 1.0)
    end
    return new_shifts
end


function apply_contact_shift(
    ecosystem::AbstractMatrix,
    contact_shifts::Dict{K, Vector{T}},
    noise::Float64,
    density_damp::Float64)
    contact_cells = ecosystem[:,1]
    new_shifts = Dict(k => Vector{Int}() for k in keys(type_shifts))
    for k in keys(type_shifts)
        local_damp = damp_decay(contact_cells[k], density_damp)
        new_shifts[k] = apply_shift(contact_shifts, noise, local_damp)
    end
    return new_shifts
end

function apply_wave_shift(
    ecosystem::AbstractMatrix,
    wave_shifts::Dict,
    density_damp::Float64,
    wave_damp::Float64)
    new_shifts = Dict(k => Vector{Int}() for k in keys(wave_shifts))
    for k in keys(wave_shifts)
        local_damp = damp_decay(contact_cells[k], density_damp)
    end
end


function ecosystem_shifts(
    ecosystem::AbstractMatrix,
    base_cells::BaseCellType)
    type_shifts = base_cells.type_shifts
    contact_shifts = base_cells.contact_shifts
    density_damp = base_cells.density_damp
    wave_damp = base_cells.wave_damp
    for cell in 1:size(ecosystem,1)
        for wave in 1:size(ecosystem,2)
            
        end
    end
    
end

function add_ecosystems(
    tissue::Tissue;)::Tissue
    base_cells = tissue.base_cell_types
    ecosystems = Dict(k => Dict() for k in keys(tissue.cell_graph))
    for eco in keys(ecosystems)
        local_niche, dist = get_neighborhood_graph(
            tissue.cell_graph,
            eco,
            depth = base_cells.wave_diffusion)
        local_eco = get_ecosystem(
            local_niche,
            dist,
            tissue.cell_labels,
            base_cells)
        local_shift = ecosystem_shifts(
            local_eco,
            base_cells)
    end
end