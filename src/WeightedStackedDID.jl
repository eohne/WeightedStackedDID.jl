
module WeightedStackedDID
    
using DataFrames, FixedEffectModels

export stacked_did_reg, compute_weights!, create_sub_exp

include("compute_weights.jl");
include("create_sub_experiment.jl");

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
end