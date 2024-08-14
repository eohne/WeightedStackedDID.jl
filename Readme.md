Weighted Stacked Regressions (Wing, Freedman, and Hollingsworth, 2024) - [NBER Paper](https://www.nber.org/system/files/working_papers/w32054/w32054.pdf)  

Following this code example: [Link](https://github.com/hollina/stacked-did-weights/tree/main?tab=readme-ov-file)

```julia
julia> using CSV, DataFrames, FixedEffectModels

julia> data = CSV.File(raw"data\acs1860_unins_2008_2021.csv") |> DataFrame;

julia> data.adopt_year[ismissing.(data.adopt_year)] .=0; # Never Treated Adoption time must be set to zero

julia> wsd_res = stacked_did_reg(data,                                                                                                                  
           yvar=:unins,                                                                                                                       
           x_vars=Symbol[], # not needed is automatically set to empty                                                                        
           fes = Symbol[], # not needed is automatically set to empty                                                                         
           timeID=:year,                                                                                                                      
           groupID=:statefip,                                                                                                                 
           adoptionTime=:adopt_year,                                                                                                          
           kappa_pre=3,                                                                                                                       
           kappa_post=2,                                                                                                                      
           cluster=:statefip)
                                FixedEffectModel
================================================================================
Number of obs:                      336  Converged:                         true
dof (model):                          5  dof (residuals):                     39
R²:                               0.391  R² adjusted:                      0.369
F-statistic:                    5.57619  P-value:                          0.001
R² within:                        0.018  Iterations:                           2
================================================================================
             Estimate  Std. Error     t-stat  Pr(>|t|)    Lower 95%    Upper 95%
────────────────────────────────────────────────────────────────────────────────
ttt::-3   0.000638982  0.00451112   0.141646    0.8881  -0.00848562   0.00976358
ttt::-2  -0.00374609   0.00317994  -1.17804     0.2459  -0.0101781    0.00268595
ttt::0   -0.0181662    0.00417973  -4.34627     <1e-04  -0.0266205   -0.00971193
ttt::1   -0.0260809    0.00648439  -4.0221      0.0003  -0.0391968   -0.012965
ttt::2   -0.0303825    0.00808953  -3.75578     0.0006  -0.0467451   -0.0140199
================================================================================
```

Calculate the Trimmed Aggregate ATT (average post-treatment effect). Also shows the average pre-treatment effect.  

```julia
julia> agg_to_ATT(wsd_res)

             Trimmed Aggregate ATT
           (average post-treatment effect)
=====================================================
Dof   : 51.0
T-ATT : -0.0219                  Std. Error: 0.0056
FVal  : 15.0847                  P-Value   : 0.0003
-----------------------------------------------------
Average pre-treatment effect:
Dof   : 51.0
Effect: -0.002                   Std. Error: 0.0032
FVal  : 0.4137                   P-Value   : 0.523
=====================================================
```