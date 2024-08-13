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