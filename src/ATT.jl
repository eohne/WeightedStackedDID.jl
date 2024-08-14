function _lh_test_stats_F(hyp_m::Tuple,β,Σ,dof)
    lhs,rhs = first(hyp_m),last(hyp_m)
    valHyp = lhs' * β - rhs
    vcovHyp = lhs' * Σ * lhs
    SSH = valHyp * inv(vcovHyp) * valHyp
    SSH = SSH/1
    std =  sqrt(vcovHyp)
    p = ccdf(FDist(1,dof),SSH)
    return (;st_err = std, FVal=SSH, PVal = p)
    end

mutable struct WSDID_AGG 
    post_val::Float64
    post_se::Float64
    post_f::Float64
    post_p::Float64

    pre_val::Float64
    pre_se::Float64
    pre_f::Float64
    pre_p::Float64
    dof::Int
end
function Base.show(io::IO, obj::WSDID_AGG)
res="\n             Trimmed Aggregate ATT
           (average post-treatment effect)
=====================================================
Dof   : $(round(obj.dof,digits=0))
T-ATT : $(round(obj.post_val,digits=4)) \t\t Std. Error: $(round(obj.post_se,digits=4))
FVal  : $(round(obj.post_f,digits=4)) \t\t P-Value   : $(round(obj.post_p,digits=4))
-----------------------------------------------------
Average pre-treatment effect:
Dof   : $(round(obj.dof,digits=0))
Effect: $(round(obj.pre_val,digits=4)) \t\t\t Std. Error: $(round(obj.pre_se,digits=4))
FVal  : $(round(obj.pre_f,digits=4)) \t\t\t P-Value   : $(round(obj.pre_p,digits=4))
=====================================================\n"
println(io,res)
end

"""
    agg_to_ATT(model::FixedEffectModel)

Calculates the Trimmed Aggregate ATT (average post-treatment effect) of a weighted stacked DiD regression.
Also calculates the average pre-treatment effect.
"""
function agg_to_ATT(model::FixedEffectModel)
    cnames = coefnames(model)
    pre_hm = occursin.(r"\:\:\-[0-9]",cnames)
    pre_hm = pre_hm./sum(pre_hm)
    post_hm = occursin.(r"\:\:[0-9]",cnames)
    post_hm = post_hm./sum(post_hm)
    β = coef(model)
    Σ = vcov(model)
    deg_free = dof_residual(model)
    pre_val = pre_hm' * β
    post_val = post_hm' * β
    pre_test = _lh_test_stats_F((pre_hm,0.),β,Σ,deg_free)
    post_test = _lh_test_stats_F((post_hm,0.),β,Σ,deg_free)
    return WSDID_AGG(post_val, post_test.st_err, post_test.FVal, post_test.PVal, pre_val,pre_test.st_err,pre_test.FVal,pre_test.PVal,deg_free)
end


