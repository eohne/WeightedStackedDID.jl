function compute_weights!(dataset::DataFrame; treatedVar::Symbol, eventTimeVar::Symbol, subexpVar::Symbol)
    # Step 1: Compute stack - time counts for treated and control
    transform!(dataset, eventTimeVar => ByRow( x -> begin
        total = count(==(x), dataset[!, eventTimeVar])
        treated_count = sum(dataset[!, treatedVar] .== 1)
        control_count = sum(dataset[!, treatedVar] .== 0)
        return (total, treated_count, control_count)
    end) => [:stack_n, :stack_treat_n, :stack_control_n])

    # Step 2: Compute sub_exp-level counts
    transform!(groupby(dataset, [subexpVar, eventTimeVar]), treatedVar => length => :sub_n )
    transform!(groupby(dataset, [subexpVar, eventTimeVar]), treatedVar => (x-> sum(isequal.(x,1))) => :sub_treat_n )
    transform!(groupby(dataset, [subexpVar, eventTimeVar]), treatedVar => (x-> sum(isequal.(x,0))) => :sub_control_n )


    # Step 3: Compute sub-experiment share of totals
    transform!(dataset, [:sub_n, :stack_n] => ByRow((x,y)->x/y) => :sub_share)
    transform!(dataset, [:sub_treat_n, :stack_treat_n] => ByRow((x,y)->x/y) => :sub_treat_share)
    transform!(dataset, [:sub_control_n, :stack_control_n] => ByRow((x,y)->x/y) => :sub_control_share)

    # Step 4: Compute weights for treated and control groups
    dataset[ !,:stack_weight] .= 1.
    dataset[ isequal.(dataset[:, treatedVar], 0), :stack_weight] .= dataset[isequal.(dataset[:, treatedVar], 0), :sub_treat_share] ./ dataset[isequal.(dataset[:, treatedVar], 0), :sub_control_share]

    return dataset
end