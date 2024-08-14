module WeightedStackedDID
    
using DataFrames, FixedEffectModels, Distributions
using PrecompileTools: @setup_workload, @compile_workload

export stacked_did_reg, compute_weights!, create_sub_exp, agg_to_ATT, WSDID_AGG

include("compute_weights.jl");
include("create_sub_experiment.jl");
include("ATT.jl")

"""
    stacked_did_reg(data::DataFrame;yvar::Symbol, timeID::Symbol, unitID::Symbol, cohort::Symbol, kappa_pre::Int, kappa_post::Int,x_vars::Vector{Symbol}=Symbol[],fes::Vector{Symbol}=Symbol[],cluster::Union{Nothing,Symbol}=nothing,contrasts::Dict{Symbol, DummyCoding}=Dict{Symbol, DummyCoding}(),multi_thread::Bool=true)

Implements the Weighted Stacked Difference-in-Differences of Wing, Freedman, and Hollingsworth (2024).

# Arguments  

 * data: The DataFrame
 * yvar: Y variable
 * timeID: The calander time variable (assumes an integer not a date variable)
 * unitID: The unit/individual ID
 * cohort: Time period of first treatment (assumes an integer not a date variable)
 * kappa_pre: Number of pre periods
 * kappa_post: Number of post periods
 * x_vars: Vector of control variables
 * fes: Vector of additional fixed effects
 * cluster: The cluster variable (currently only supports one). If not supplied it is set at the level of the unitID
 * contrasts: contrasts for DummyVariables. Required if you included a Dummy Variable in x_vars
 * multi_thread: Whether to use multithreading or not in building the stacks

 # Output

 FixedEffectModel 

"""
function stacked_did_reg(data::DataFrame;yvar::Symbol, timeID::Symbol, unitID::Symbol, cohort::Symbol, kappa_pre::Int, kappa_post::Int,x_vars::Vector{Symbol}=Symbol[],fes::Vector{Symbol}=Symbol[],cluster::Union{Nothing,Symbol}=nothing,contrasts::Dict{Symbol, DummyCoding}=Dict{Symbol, DummyCoding}(), multi_thread::Bool=true)
    dataset = select(data, unique(vcat(yvar,x_vars,fes,timeID,unitID,cohort)))    
    dropmissing!(dataset)
    
    _formula=term(yvar)~sum(vcat(term(:event_time_treat),term.(x_vars),fe.(term.(vcat(fes,:treat,:event_time)))))
    if isnothing(cluster)
        cluster = unitID
        @info "\nCluster variable was not provided. Clustering will happen on the $(string(unitID)) level!\n"
    end
    if size(dataset,1) != size(data,1)
        @info "\n$(size(data,1) - size(dataset,1)) observations were dropped due to missingness!\n"
    end
    # Initialize a list (Vector) to store the sub experiments
    allowmissing!(dataset, cohort)
    dataset[!,cohort][isequal.(data[:,cohort],0)] .=missing
    events = unique(skipmissing(dataset[:, cohort]))
    if multi_thread
        stacked_dtc = [DataFrame() for i in 1:Threads.nthreads()]
        # Loop over the events and create a data set for each one
        Threads.@threads for j in events
        append!(stacked_dtc[Threads.threadid()],
            create_sub_exp(
                    dataset,
                    timeID=timeID,
                    unitID=unitID, 
                    cohort=cohort, 
                    focalcohort=j,
                    kappa_pre=kappa_pre,
                    kappa_post=kappa_post
                ))
        end
        # Remove the sub-experiments that are not feasible
        stacked_dtc = vcat(stacked_dtc...)
    else
        stacked_dtc = DataFrame()
        # Loop over the events and create a data set for each one
        for j in events
        append!(stacked_dtc,
            create_sub_exp(
                    dataset,
                    timeID=timeID,
                    unitID=unitID, 
                    cohort=cohort, 
                    focalcohort=j,
                    kappa_pre=kappa_pre,
                    kappa_post=kappa_post
                ))
        end
    end
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


@setup_workload begin
    using CSV
    data = CSV.File(raw"data\acs1860_unins_2008_2021.csv") |> DataFrame
    data[!,:x1] = rand(size(data,1))
    data.adopt_year[ismissing.(data.adopt_year)].=0
    @compile_workload begin
        res =stacked_did_reg(data;yvar=:unins, timeID=:year, unitID=:statefip, cohort=:adopt_year, kappa_pre=3, kappa_post=2,cluster=:statefip);
        res =stacked_did_reg(data;yvar=:unins, timeID=:year, unitID=:statefip, cohort=:adopt_year, kappa_pre=3, kappa_post=2,cluster=:statefip,multi_thread=false);
        res =stacked_did_reg(data;yvar=:unins, timeID=:year, unitID=:statefip, cohort=:adopt_year, kappa_pre=3, kappa_post=2,cluster=:statefip,fes=[:statefip],x_vars=[:x1]);
        agg_to_ATT(res);
    end
end

end