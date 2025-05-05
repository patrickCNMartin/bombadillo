

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

# function apply_shift(
#     shifts::AbstractVector{Int},
#     noise::Float64,
#     damp::AbstractVector{Float64})
#     for i in eachindex(shifts)
#         local_shift = 0
#         for d in damp
#             local_shift += apply_noise(shifts[i] * d, noise)
#         end
#             new_shifts[i] = sum(local_shift)
#         end
        
#     return new_shifts
# end

function compute_type_shift(
    type_shifts::Dict,
    noise::Float64)
    new_shifts = Dict(k => Vector{Int}() for k in keys(type_shifts))
    for (k, v) in type_shifts
        local_shift = similar(v)
        for s in eachindex(local_shift)
            local_shift[s] = compute_shifts(v[s], noise)
        end
        new_shifts[k] = local_shift
    end
    return new_shifts
end


function compute_contact_shift(
    ecosystem::AbstractMatrix,
    contact_shifts::Dict,
    noise::Float64,
    density_damp::Float64)
    contact_cells = ecosystem[:,1]
    new_shifts = Dict(k => Vector{Int}() for k in keys(contact_shifts))
    for (k, v) in contact_shifts
        local_shift = similar(v)
        for s in eachindex(local_shift)
            damp = density_decay(contact_cells[k], density_damp)
            for d in damp
                if round(v[s] * d) == 0
                    break
                end
                local_shift[s] += compute_shifts(v[s] * d, noise)
            end
        end
        new_shifts[k] = local_shift
    end
    return new_shifts
end

function compute_wave_shift(
    ecosystem::AbstractMatrix,
    wave_shifts::Dict,
    noise::Float64,
    density_damp::Float64,
    wave_damp::Float64)
    new_shifts = Dict(k => Vector{Int}() for k in keys(wave_shifts))
    for (k, v) in wave_shifts
        local_shift = similar(v)
        for depth in axes(ecosystem, 2)
            cell_dense = ecosystem[k,depth]
            damp_w = wave_decay(depth, wave_damp)
            for s in eachindex(local_shift)
                damp_d = density_decay(cell_dense, density_damp)
                for d in damp_d
                    if round(v[s] * (d * damp_w)) == 0
                        break
                    end
                    local_shift[s] += compute_shifts(v[s] * (d * damp_w) , noise)
                end
            end
        end
        new_shifts[k] = local_shift
    end
    return new_shifts
end

function apply_shift(
    gene_vec::AbstractVector{Int},
    gene_set::Dict{K, Vector{Int}},
    gene_shift::Dict{K, Vector{Int}}) where K
    
    for cell in keys(gene_set)
        genes = gene_set[cell]
        shifts = gene_shift[cell]
        for (gene_array, shift_array) in zip(genes, shifts)
            local_index = gene_array
            local_shift = shift_array
            if local_index + local_shift < 1
                shift_array = -local_index + 1
            end
            if local_index + local_shift > maximum(gene_vec)
                local_shift =  maximum(gene_vec) - local_index
            end
            gene_vec[local_index] = local_index + local_shift
        end
    end
    
    return gene_vec
end

function aggregate_shift(
    type_shift::Dict,
    type_set::Dict,
    contact_shift::Dict,
    contact_set::Dict,
    wave_shift::Dict,
    wave_set::Dict,
    n_genes::Int64)
    gene_vector = collect(1:n_genes)
    gene_vector = apply_shift(gene_vector, type_set, type_shift)
    gene_vector = apply_shift(gene_vector, contact_set, contact_shift)
    gene_vector = apply_shift(gene_vector, wave_set, wave_shift)
    return gene_vector
end

function ecosystem_shifts(
    ecosystem::AbstractMatrix,
    base_cells::BaseCellType)
    noise = base_cells.noise
    density_damp = base_cells.density_damp
    wave_damp = base_cells.wave_damp
    n_genes = base_cells.n_genes

    type_shift = base_cells.type_shifts
    type_set = base_cells.type_sets
    type_shift = compute_type_shift(type_shift, noise)


    contact_shift = base_cells.contact_shifts
    contact_set = base_cells.contact_sets
    contact_shift = compute_contact_shift(ecosystem,
        contact_shift,
        noise,
        density_damp)

    wave_shift = base_cells.wave_shifts
    wave_set = base_cells.wave_sets
    wave_shift = compute_wave_shift(ecosystem,
        wave_shift,
        noise,
        density_damp,
        wave_damp)

    cell_shift = aggregate_shift(type_shift,
        type_set,
        contact_shift,
        contact_set,
        wave_shift,
        wave_set,
        n_genes)
    return cell_shift
    
end

function collapse_ecosystem(
    ecosystem::Matrix{Float64},
    collapse_to::Union{Int64, String})
    if typeof(collapse_to) == Int
        return round.(ecosystem ./ collapse_to) .* collapse_to
    else
        return map(x -> x > 0 ? one(x) : zero(x), ecosystem)
    end
    
end

function add_ecosystems(
    tissue::Tissue;
    collapse_to::Union{Int64, String} = 1)::Tissue
    base_cells = tissue.base_cell_types
    ecosystems = Dict(k => (zeros(0, 0), Float64[]) for k in keys(tissue.cell_graph))
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
        if !isnothing(collapse_to)
            local_eco = collapse_ecosystem(local_eco, collapse_to)
        end
        local_shift = ecosystem_shifts(
            local_eco,
            base_cells)
        ecosystems[eco] = (local_eco,local_shift)
    end
    new_labels = ecosystem_cell_types(ecosystems,tissue.cell_labels)
    tissue.ecosystems = ecosystems
    tissue.cell_labels = new_labels
    return tissue
end