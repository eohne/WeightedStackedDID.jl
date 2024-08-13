using CSV, DataFrames, FixedEffectModels
data = CSV.File(raw"data\acs1860_unins_2008_2021.csv") |> DataFrame;

larger_na_false(x::Number,y::Number) = x>y
larger_na_false(x::Missing,y::Number) = false

function create_sub_exp(dataset::DataFrame; timeID::Symbol, groupID::Symbol, adoptionTime::Symbol, focalAdoptionTime::Int, kappa_pre::Int, kappa_post::Int)

    # Copy dataset
    dt_temp = deepcopy(dataset)

    # Determine earliest and latest time in the data.
    # Used for feasibility check later
    minTime = minimum(dt_temp[!, timeID])
    maxTime = maximum(dt_temp[!, timeID])

    # Include only the treated groups and the clean controls that adopt at least kappa_post periods after the focal atime.
    subset!(dt_temp, adoptionTime => ByRow(x-> (isequal.(x,focalAdoptionTime)) | larger_na_false(x, focalAdoptionTime + kappa_post) | isequal(x,true) | ismissing(x) ) )
    
    # Limit to time periods inside the event window defined by the kappas
    subset!(dt_temp, timeID => ByRow(x-> (focalAdoptionTime-kappa_pre) <= x <= (focalAdoptionTime + kappa_post))) 

    # Make treatment group dummy
    dt_temp[:, :treat] .= 0;
    dt_temp[ isequal.(dt_temp[:, adoptionTime], focalAdoptionTime), :treat] .= 1;

    # Make a post variable
    dt_temp[:, :post] .= ifelse.(dt_temp[:, timeID] .>= focalAdoptionTime, 1, 0)

    # Make event time variable
    dt_temp[:, :event_time] .= dt_temp[:, timeID] .- focalAdoptionTime

    # Create a feasible variable
    dt_temp[:, :feasible] .= ifelse((focalAdoptionTime - kappa_pre) >= minTime && 
                                    (focalAdoptionTime + kappa_post) <= maxTime, 1, 0)

    # Make a sub experiment ID
    dt_temp[:, :sub_exp] .= focalAdoptionTime

    return dt_temp
end;

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


"""
    stacked_did_reg(data::DataFrame;yvar::Symbol, timeID::Symbol, groupID::Symbol, adoptionTime::Symbol, kappa_pre::Int, kappa_post::Int,x_vars::Vector{Symbol}=Symbol[],fes::Vector{Symbol}=Symbol[],cluster::Union{Nothing,Symbol}=nothing,contrasts::Dict{Symbol, DummyCoding}=Dict{Symbol, DummyCoding}())

Implements the Weighted Stacked Difference-in-Differences of Wing, Freedman, and Hollingsworth (2024).

"""
function stacked_did_reg(data::DataFrame;yvar::Symbol, timeID::Symbol, groupID::Symbol, adoptionTime::Symbol, kappa_pre::Int, kappa_post::Int,x_vars::Vector{Symbol}=Symbol[],fes::Vector{Symbol}=Symbol[],cluster::Union{Nothing,Symbol}=nothing,contrasts::Dict{Symbol, DummyCoding}=Dict{Symbol, DummyCoding}())
    dataset = select(data, unique(vcat(yvar,x_vars,fes,timeID,groupID,adoptionTime)))    
    dropmissing!(dataset)
    
    _formula=term(yvar)~sum(vcat(term(:event_time_treat),term.(x_vars),fe.(term.(vcat(fes,:treat,:event_time)))))
    if isnothing(cluster)
        cluster = groupID
        @info "\nCluster variable was not provided. Clustering will happen on the $(string(groupID)) level!\n"
    end
    if size(dataset,1) != size(data,1)
        @info "\n$(size(data,1) - size(dataset,1)) observations were dropped due to missingness!\n"
    end
    events = unique(skipmissing(dataset[:, adoptionTime]))
    # Initialize a list (Vector) to store the sub experiments
    allowmissing!(dataset, adoptionTime)
    dataset[:,adoptionTime][isequal.(data[:,adoptionTime],0)] .=missing
    stacked_dtc = DataFrame()
    # Loop over the events and create a data set for each one
    for j in events
    append!(stacked_dtc,
        create_sub_exp(
                dataset,
                timeID=timeID,
                groupID=groupID, 
                adoptionTime=adoptionTime, 
                focalAdoptionTime=j,
                kappa_pre=kappa_pre,
                kappa_post=kappa_post
            ))
    end
    # Remove the sub-experiments that are not feasible
    subset!(stacked_dtc, :feasible =>ByRow( x-> x==1));

    compute_weights!(
        stacked_dtc,
        treatedVar = :treat,
        eventTimeVar = :event_time,
        subexpVar = :sub_exp)
        
    transform!(stacked_dtc, [:event_time,:treat] =>  ByRow(*) =>:event_time_treat )

    weight_stack = reg(stacked_dtc, _formula, 
                            Vcov.cluster(cluster),
                            weights =:stack_weight, contrasts = merge(Dict(:event_time_treat => DummyCoding(base=-1)),contrasts))
    coefnames(weight_stack) .= replace.(coefnames(weight_stack), "event_time_treat: "=>"ttt::")

    return weight_stack
end

using CSV, DataFrames, FixedEffectModels
data = CSV.File(raw"data\acs1860_unins_2008_2021.csv") |> DataFrame;
data.adopt_year[ismissing.(data.adopt_year)] .=0; # Never Treated Adoption time must be set to zero
stacked_did_reg(data,
    yvar=:unins,
    x_vars=Symbol[], # not needed is automatically set to empty
    fes = Symbol[], # not needed is automatically set to empty
    timeID=:year,
    groupID=:statefip,
    adoptionTime=:adopt_year,
    kappa_pre=3,
    kappa_post=2,
    cluster=:statefip)

