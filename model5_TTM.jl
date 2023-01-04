# Model 4 + constrain outbound volume by plants
#   
# Assume CPU has the same freight penalty as PPD, and then zero out in post processing
#
#using Pkg
#Pkg.add("Cbc")

using DataFrames # Includes PrettyTables
using Pipe
using Printf
using XLSX
using Statistics
using VegaLite
using VegaDatasets
#using JSON3
using JSON


include("read_geodata.jl")
using .LoadGeoData
varinfo(LoadGeoData)

include("get_forecast.jl")
using .LoadForecast
varinfo(LoadForecast)

@pipe combine(fcst_df, :OLMUnits => sum => :OLMUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#Copy this to the clipboard and paste into a file for the capacity to use in the model
#clipboard(sprint(show, "text/tab-separated-values", fcst_df))

maxship = Dict("KENTON" => 4600000.0,
     "VISALIA" => 6400000.0,
     "SHELBYVILLE" => 11100000.0,
     "CLARKSVILLE" => 1600000.0,
     "PITTSTON" => 1400000.0)

shipsites = keys(maxship)
plants = shipsites

const DATA_DIR = joinpath(@__DIR__, "Data")

#Get the historical shipments to determine the ship-to location splits by part type
#Save all this for when we run the model from the forecast...
        #= ship_df = DataFrames.DataFrame(
            XLSX.readtable(joinpath(DATA_DIR, "ship_history_last6.xlsx"), "Sheet 1"; infer_eltypes=true)...
        )
        describe(ship_df) # Note missing values for MfgSite
        dropmissing!(ship_df, :MfgSiteName)
        describe(ship_df) # Note missing values for MfgSite

        # Convert MUnits from Any to Float64
        ship_df[!, "MUnits"] = convert.(Float64, ship_df[!, "MUnits"])

        # Filter to just desired Print Group
        #ship_df = filter(:PrintGroup => ==("COKECORP"), ship_df)

        # Combine the City and State to a Loc
        ship_df.Loc = string.(ship_df.MSTCity, ", ", ship_df.MSTState)

        # Create a "Lid" or "Cup" Part type in the TypeFix column
        # Not needed any more, but we'll add it anyway
        ship_df.PTypeFix = string.(first.(ship_df.PartType))
        ship_df.PTypeFix = ifelse.(string.(first.(ship_df.PartType)) .== "L", "LID", "CUP") 
        ship_df.PTypeFix = Symbol.(ship_df.PTypeFix)
        # check that it works
        ship_df[ship_df[:, :PTypeFix] .== :LID, :]

        # Get rid of the PlasticInenuity
        ship_df = ship_df[ship_df[:, :ShipSiteName] .!= "PLASTICS INGENUITY", :]

        @pipe combine(ship_df, :MUnits => sum => :MUnits) |> 
            PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
        #dropmissing!(demand_df, :MSTZip)

        # Just include where we have distances
        ship_df = innerjoin(ship_df, DataFrames.select(distw_df, :Dest), on = [:Loc => :Dest])

        @pipe combine(ship_df, :MUnits => sum => :MUnits) |> 
            PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

        # How much is PPD freight
        @pipe combine(filter(:FrtTermsFix => x -> x == "PPD", ship_df), :MUnits => sum => :MUnits) |> 
            PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

        # For the first pass, make everything PPD.  Deal with CPU by lane later.
        ship_df[!, :FrtTermsFix] .= "PPD"

        #Check
        @pipe combine(filter(:FrtTermsFix => x -> x == "PPD", ship_df), :MUnits => sum => :MUnits) |> 
            PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

        # consolidate
        ship_df = combine(groupby(ship_df, [:MSTName, :Loc, :FrtTermsFix, :PrintGroup, :PartType]), :MUnits => sum => :MUnits)

        rename!(ship_df, :FrtTermsFix => :FTerms)

        # Make sure we haven't lost anything
        @pipe combine(ship_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

        #TODO: Need a check somewhere that we don't have demand for part types not in 
        #       the history

        # Do special stuff here like changing part types, etc.
        # 
        ship_df.PartType = replace.(ship_df.PartType, ("LCRS-032" => "LCRSSL32"))
        ship_df.PartType = replace.(ship_df.PartType, ("DMSL0120" => "DMS-0120"))


        #recombine to eliminate duplicates
        ship_df = combine(groupby(ship_df, Between(:MSTName, :PartType)), :MUnits => sum => :MUnits)

        #Check the smallest
        ship_df[partialsortperm(ship_df.MUnits, 1:10, rev=false), :]

        # Calculate the percent of each PartType at each DC
        @pipe transform!(groupby(ship_df, :PartType), :MUnits => sum) |>
            transform!(_, [:MUnits, :MUnits_sum] => ByRow(/) => :Pct) |>
            select!(_, [:MSTName, :Loc, :FTerms, :PrintGroup, :PartType, :Pct])

        # Copy to the clipboard and review
        clipboard(sprint(show, "text/tab-separated-values", ship_df))

        # Aggregate the forecast by parttype
        dmd_df = combine(groupby(fcst_df, [:PrintGroup, :PartType]), :OLMUnits => sum => :MUnits)

        # More special stuff here

        # Replace the LCRS-032 with the LCRSL032 coming from Pittston
        #dmd_df.PartType = replace.(dmd_df.PartType, ("LCRS-032" => "LCRSSL32"))

        # Aggregate again
        dmd_df = combine(groupby(dmd_df, [:PrintGroup, :PartType]), :MUnits => sum => :MUnits)

        dmd_df = leftjoin(ship_df, dmd_df, on = [:PrintGroup, :PartType])

        # Remove any that no longer are in the forecast
        dropmissing!(dmd_df, :MUnits)    

        # Calculate the volume to each location
        @pipe transform!(dmd_df, [:Pct, :MUnits] => ByRow(*) => :MUnits) |>
            select!(_, Not([:Pct]))

        # Did we lose anything?
        @pipe combine(dmd_df, :MUnits => sum => :MUnits) |> 
            PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

        # Copy to the clipboard
        clipboard(sprint(show, "text/tab-separated-values", dmd_df))
 =#

ship_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_history.xlsx"), "Sheet 1"; infer_eltypes=true)
)
describe(ship_df) # Note missing values for MfgSite
dropmissing!(ship_df, :MfgSiteName)
describe(ship_df) # Note no missing values for MfgSite

# Convert MUnits from Any to Float64
ship_df[!, "MUnits"] = convert.(Float64, ship_df[!, "MUnits"])

# Combine the City and State to a Loc
ship_df.Loc = string.(ship_df.MSTCity, ", ", ship_df.MSTState)

# Create a "Lid" or "Cup" Part type in the TypeFix column
# Not needed any more, but we'll add it anyway
ship_df.PTypeFix = string.(first.(ship_df.PartType))
ship_df.PTypeFix = ifelse.(string.(first.(ship_df.PartType)) .== "L", "LID", "CUP") 
ship_df.PTypeFix = Symbol.(ship_df.PTypeFix)
# check that it works
ship_df[ship_df[:, :PTypeFix] .== :LID, :]

# Get rid of the PlasticIngenuity
ship_df = ship_df[ship_df[:, :ShipSiteName] .!= "PLASTICS INGENUITY", :]

@pipe combine(ship_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
#dropmissing!(demand_df, :MSTZip)

# Get the root part type Tables
rootpt_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "CupPlantRootPartTypes.xlsx"), 
        "RootPartType"; infer_eltypes=true)
)

rename!(rootpt_df, "Part Type" => :PartType, "Root Part Type" => :RootPartType)

rootpt_df = select(rootpt_df, :PartType, :RootPartType)

rootpt_df = unique(rootpt_df, :PartType)

ship_df = innerjoin(ship_df, rootpt_df, on = :PartType)

#How many missing CupPlantRootPartTypes
describe(ship_df)

#TODO:  Add logic to create a new root part type field.  If the RootPart type is a LID,
#    use the part type instead (i.e., 'Hot PS')

ship_df.RPType = ifelse.(ship_df.PTypeFix .== :LID, ship_df.PartType, ship_df.RootPartType)

describe(ship_df)

#clipboard(sprint(show, "text/tab-separated-values", ship_df))

# Just include where we have distances
ship_df = innerjoin(ship_df, DataFrames.select(distw_df, :Dest), on = [:Loc => :Dest])

@pipe combine(ship_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# How much is PPD freight
@pipe combine(filter(:FrtTermsFix => x -> x == "PPD", ship_df), :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Keep both PPD and CPU in the demand file to get total demand, not just PPD
demand_df = combine(groupby(ship_df, 
    [:MSTName, :Loc, :FrtTermsFix, :RPType]), :MUnits => sum => :MUnits)
rename!(demand_df, :FrtTermsFix => :FTerms)

# Make sure we haven't lost anything
@pipe combine(demand_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Build the various data structures to hold the supply and demand info

#supply_df = DataFrames.DataFrame(
#    XLSX.readtable(joinpath(DATA_DIR, "COKE", "capacity.xlsx"), "Plan"; infer_eltypes=true)...
#)

supply_df = combine(groupby(ship_df, [:MfgSiteName, :RPType]), :MUnits => sum => :MUnits)
@pipe supply_df |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
clipboard(sprint(show, "text/tab-separated-values", supply_df))

dropmissing!(supply_df)
rename!(supply_df, :MUnits => :Capacity, :MfgSiteName => :MfgSite)
supply_df[!, :Capacity] = convert.(Float64, supply_df[!, :Capacity])
#select!(supplypep_df, Not([:Shipped]))
subset!(supply_df, :Capacity => ByRow(!=(0)))

#Aggregate
supply_df = combine(groupby(supply_df, [:MfgSite, :RPType]), :Capacity => sum => :Capacity)

@pipe combine(supply_df, :Capacity => sum => :Capacity) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
    
rptypes = unique(supply_df.RPType)

struct Supply
    plant::String
    rptype::String
end

#[
#:LID; :CUP
#]

supply_dict = Dict(Supply.(supply_df.MfgSite, supply_df.RPType) .=> supply_df.Capacity)
supplies  = keys(supply_dict)
for s in supplies
    println(s.plant, " ", s.rptype, " ", supply_dict[s])
end

# Don't do the following - use the sparse matrix logic from the following link
# https://shuvomoy.github.io/blogs/posts/Solving-transportation-problem-in-Julia-Jump/ 
#= # Add missing elements with a 0 supply
for p in plants
    for t in ptypes
        #println(Demand(c, p))
        # if !(Supply(p, t) in keys(supply_dict)) 
        #     supply_dict[Supply(p, t)] = 0
        # end
        supply_dict[Supply(p, t)] = get(supply_dict, Supply(p, t), 0)
    end
end
 =#
supply_dict

# Customers as a Named tuple
cust_nt = Tables.rowtable(DataFrames.select(demand_df, :MSTName, :Loc, :FTerms))
unique(cust_nt)

# Better to use a Struct?
struct Cust
    mstname::String 
    loc::String 
    fterms::String 
end

# populate with the data from the dataframe
customers = Cust.(demand_df.MSTName, demand_df.Loc, demand_df.FTerms)
customers = unique(customers)

eltype(customers)

# build the demand dict using the Cust structure
#demand_dict = Dict(structfromnt.(Cust, Tables.rowtable(DataFrames.select(dmd_df, :MSTName, :Loc, :FTerms))) .=> dmd_df.MUnits)
# Another way to do it
#demand_dict = Dict(Pair.(structfromnt.(Cust, Tables.rowtable(DataFrames.select(dmd_df, :MSTName, :Loc, :FTerms))), dmd_df.MUnits))

# Use this as the key for a demand dict
struct Demand
    customer::Cust
    rptype::String
end

# Solution from Bogumit Kaminski from Julialang@discoursemail.com

demand_dict = Dict(Demand.(Cust.(demand_df.MSTName, demand_df.Loc, demand_df.FTerms), 
    demand_df.RPType) .=> demand_df.MUnits)
@printf("%d", sum(values(demand_dict)))
demands = keys(demand_dict)

#for d in demands
#    println(d.customer, " ", d.ptype, " ",  demand_dict[d])
#end

#demand_key = Demand.(Cust.(dmd_df.MSTName, dmd_df.Loc, dmd_df.FTerms), dmd_df.PartType)

#keys(demand_dict)
#demand_dict[Demand(Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD"), :LID)]

# Don't do this either.  See above for the supplies
#= # Add missing elements with a 0 demand
# Can I skip this and use the get() function with a zero default value
for c in customers
    for t in ptypes
        #println(Demand(c, p))
        #if !( Demand(c, t) in keys(demand_dict)) 
        #    demand_dict[Demand(c, t)] = 0
        #end
        demand_dict[Demand(c, t)] = get(demand_dict, Demand(c, t), 0)
    end
end
 =#

 # sum the demand values 
@printf("%d", sum(values(demand_dict)))
@pipe [sum(values(demand_dict))] |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
# make sure it's the same
@pipe combine(demand_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))


# create a dict storing for each customer the products they take
custprod = Dict{Cust, Array}()
for c in customers
    dummy_array = String[]
    for d in demands
        if d.customer == c
            #for p in ptypes
                #println(i, j, products[i] == routes[j].p)
            #if p == d.ptype
            push!(dummy_array, d.rptype)
            #println(orig[products[i]])
        else
        #println("Oops, something is not right")
    end #if
  end #for
  custprod[c]=unique(dummy_array)
end #for

# Set of Origins
# plants = String[]
# for (i, j) in supply_dict
#     push!(plants, i.plant)
# end
# plants = unique(plants)

# Build the lanes (PLANT, DEST pairs)
struct Lane
    shipsite::String
    cust::Cust
end

lanes = Lane[]

for s in shipsites
    for c in customers
        push!(lanes, Lane(s, c))
    end
end

length(lanes)

# Then make a dict Lane => dist
# use 'only' function in the next step to return miles as a scalar, not a 1-d vector
dist_dict = Dict(l => only(dist_df[(dist_df[:, :Orig].== l.shipsite) .& 
    (dist_df[:, :Dest].== l.cust.loc), :Miles]) for l in lanes)

#values(dist_dict)

# Build the Routes: Lane and a product - only where the customer has demand for a product
struct Route 
    lane::Lane
    rptype::String
end

routes = Route[]
for l in lanes
    for p in custprod[l.cust]
        push!(routes, Route(l, p))
    end
end
#clipboard(sprint(show, "text/tab-separated-values", enumerate(routes)))
#enumerate(routes)

# Build the BF dictionary with costs per M for each O/D Pair
bf_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "bf_lanes.xlsx"), "Sheet1"; infer_eltypes=true)
)

struct BFLane
    plant::String
    shipsite::String
end

bf_lanes = BFLane.(bf_df.Plant, bf_df.ShipSite)
bf_dict = Dict(bf_lanes .=> bf_df.CostPerMUnit)

bf_dict[BFLane("CLARKSVILLE", "SHELBYVILLE")]

#Overwite with large freight penalty to force shipping from the mfg site
# for p in plants
#     for s in shipsites
#         if p != s
#             bf_dict[BFLane(p,s)] = 100.0
#         end
#     end
# end


# Scalar f  freight in dollars per case/M per thousand miles 
# Jan - Nov Munits * Miles / 1000 = 7,599,747 [1M Miles * MUnits Shipped] 
# Jan - Nov PPD customer Freight = $42,067,887
# 
f = 6.15 # freight in dollars per M per thousand miles 
 
# Parameter c(l)  transport cost in thousands of dollars per case on each lane ;
#            c(l) = f * dist(l) / 1000 ;
# TODO:  There might be a way to reduce the number of lanes by removing the
#  FTerms component.  Freight will e the same for PPD and CPU
# We first declare an empty dictionary and then we fill it with the values
cost = Dict{Lane,Float64}() # transport cost in thousands of dollars per M ;
#[ c[p,c] = f * dist_dict[p,c] / 1000 for p in PLANTS, c in CUST]
[ cost[l] = f * dist_dict[l] / 1000 for l in lanes ]
cost

# Define the Mz for limiting the z sites shipping to each customer
#Mz(s,c,t) = min(supply(s),demand(c,t))

# Mz = []
#     for c in customers
#         for t in ptypes
#             push!(Mz, min(demand_dict[Demand(c, t)]))
#         end
#     end

#demand_dict[Demand(Cust("PEPSI BOTTLING GROUP", "WINNIPEG, MB", "PPD"), "DMSL0120")]


using JuMP
using Gurobi
#using GLPK
#using Cbc
#using HiGHS
#model = Model(HiGHS.Optimizer)
#model = Model(Cbc.Optimizer)
#model = Model(GLPK.Optimizer) 
model = Model(Gurobi.Optimizer)
# These for Cbc
#set_optimizer_attribute(model, "logLevel", 1)
#set_optimizer_attribute(model, "allowableFractionGap", 0.1)
#set_optimizer_attribute(model, "maxSolutions", 2)
#set_optimizer_attribute(model, "maxNodes", 100)

# These for HiGHS
#set_optimizer_attribute(model, "threads", 4)
#set_optimizer_attribute(model, "mip_max_nodes", 2147483647)
#set_optimizer_attribute(model, "ipm_optimality_tolerance", 1e-08)
#set_optimizer_attribute(model, "mip_max_improving_sols", 2)
#set_optimizer_attribute(model, "ipm_iteration_limit", 1)

#unset_silent(model)

#MUnits for Mfg Plant to Ship site. 
@variable(model, y[supplies, shipsites] >= 0)

# Define variable for Ship Site Plant to customer
# route is a structure consisting of Lane (O/D Pair) and a part type
@variable(model, x[routes])

# Binary variable if Site S ships to Cust C = 1, 0 otherwise
@variable(model, z[lanes], Bin)

# Binary variable if site S is open
#@variable(model, r[shipsites], Bin)

# Constraint that each plant can ship no more than its capacity:
@constraint(model,
    supply_con[s in supplies], sum(y[Supply(s.plant, s.rptype), d] for d in shipsites) <= 
        supply_dict[Supply(s.plant, s.rptype)])

# Customer Demand constraint
@constraint(model,
    demand_con[d in demands], sum(x[Route(Lane(site, d.customer), d.rptype)] for site in shipsites) ==
    demand_dict[Demand(d.customer, d.rptype)]) 

# shipsite inbound = outbound for net of 0
@constraint(model, inout_con[s in shipsites, t in rptypes], 
    sum(y[Supply(p.plant, p.rptype), s] for p in supplies if p.rptype == t) == 
    sum(x[Route(Lane(s, c), t)] for c in customers if t in custprod[c]))

#limit the ship sites to 1 per customer
@constraint(model, maxsite_con[c in customers], sum(z[Lane(s, c)] for s in shipsites) == 1)

#limit the amount shipping through each ship site
@constraint(model, maxship_con[s in shipsites], 
    sum(x[Route(Lane(s, d.customer), d.rptype)]
    for d in demands) <= 
        maxship[s])

# Constraint to Limit the number of ship sites - use the last Big M for the relaxion help
#@constraint(model, maxopen_con, sum(r[:]) <= 5)

# Big Mz constraint for the number of sites servicing each DC
# See p 449 of AMPL
@constraint(model, avail_con[s in shipsites, d in demands], 
    x[Route(Lane(s, d.customer), d.rptype)] <= 
        min(demand_dict[Demand(d.customer, d.rptype)]) * z[Lane(s, d.customer)])

#@constraint(model, avail_con1[s in supplies],
#    y[supplies, shipsites] <= sum((supply_dict[]) * z[Lane(s, d.customer)])

#Big M for the total number of sites
#@constraint(model, openship_con[s in shipsites, d in demands], x[Lane(s, d.customer), d.ptype] <= 
#    min(demand_dict[Demand(d.customer, d.ptype)]) * r[s])

# Objective
@objective(model, Min, sum(bf_dict[BFLane(s.plant, ss)] * 
                        y[Supply(s.plant, s.rptype), ss] for s in supplies, ss in shipsites) + 
    sum(cost[r.lane] * x[r] for r in routes)) #, t in custprod[l.cust]))

#print(model)
#write_to_file(model, "model5_TTM.lp")
#latex_formulation(model)

#set_optimizer(model, HiGHS.Optimizer)

optimize!(model)

solution_summary(model) 

println("RESULTS:")
for p in supplies, s in plants
        println("Supply: $p to Site: $s => ", value(y[p, s]))
end

for r in routes #l in lanes, t in custprod[l.cust] 
    #if value(x[Lane(s, cu), t]) > 0
        println("Route: $r = " , value(x[r]), "  Freight = ",  value(x[r]) * cost[r.lane])
    #end
end

for c in customers, s in shipsites
    println("Binary Variable z[$s, $c] = ", value(z[Lane(s, c)]))
end

for s in shipsites
    println("Binary Variable r[$s] = ", value(r[s]))
end

@pipe [objective_value(model)] |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# # This is a test for the output
# for l in lanes
#         println(" $(l) = $(value(x[l]))")
#         println(l.Cust.FTerms)
# end

# Create the output DataFrame - OutBound Freight from the ship site
# recall that x[Lane(s, c), t] and lanes are Lane(shipsite, cust) demands are Demand(c, t)
# customers =  Cust(mstname, loc, fterms)
outob_df = DataFrame(ShipSite = String[], Cust = String[], Dest = String[], 
    FTerms = String[], PType = String[], MUnits = Float64[], Freight = Float64[])

for s in shipsites
    for c in customers
        for p in custprod[c]
            #println("Site: ", s, " Cust: ", c, " Prod: ", p)
            if !isapprox(value(x[Route(Lane(s, c), p)]), 0.0, atol = 1e-10, rtol = 0) #&& l.Cust.FTerms == "PPD"
                frt = if c.fterms == "PPD" value(x[Route(Lane(s, c), p)]) * cost[Lane(s, c)] else 0 end
                tupl = (s, c.mstname, c.loc, c.fterms, p, value(x[Route(Lane(s, c), p)]), frt)
                #push!(outob_df.ShipSite, s)
                #push!(outob_df.Cust, c.mstname)
                #push!(outob_df.Dest, c.loc)
                #push!(outob_df.FTerms, c.fterms)
                #push!(outob_df.PType, p)
                #push!(outob_df.MUnits, value(x[Lane(s, c), p]))
                #if c.fterms == "PPD" 
                #    push!(outob_df.Freight, value(x[Lane(s, c), p]) * cost[Lane(s, c)])
                #else
                #    push!(outob_df.Freight, 0.0)
                #end
                push!(outob_df, tupl)
            end
        end
    end
end


@pipe combine(outob_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

outob_df = innerjoin(outob_df, dist_df, on = [:ShipSite => :Orig, :Dest])
#transform!(out_df, [:MUnits, :Miles] => ByRow((a, b) -> (a * b * f / 1000)) => :Freight1)

@pipe combine(outob_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#aggregate munits by plant
@pipe combine(groupby(outob_df, :ShipSite), :MUnits => sum => :SumMUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

@pipe combine(groupby(outob_df, :ShipSite), :MUnits => sum => :SumMUnits, :Freight => sum => :Freight) |> 
    clipboard(sprint(show, "text/tab-separated-values", _))


# out_df[!, :MUnitMiles] = out_df[!, :MUnits] .* out_df[!, :Miles]  # This also works

const OUT_DIR = joinpath(@__DIR__, "Data", "TTM")

XLSX.writetable(joinpath(OUT_DIR, "model5TTMCb30Iterations.xlsx"), overwrite = true, 
    collect(DataFrames.eachcol(outob_df)), DataFrames.names(outob_df))

# agregate by lane
outplotob_df = combine(groupby(outob_df, [:ShipSite, :Dest]), :MUnits => sum => :MUnits)
#geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)
outplotob_df = innerjoin(outplotob_df, geocodes_fix, on = :Dest)

# Create the output DataFrame - BF Freight 
# y[supplies, plant], Supply = plant, ptype
outbf_df = DataFrame(MfgSite = String[], ShipSite = String[], 
    PType = String[], MUnits = Float64[], Freight = Float64[])
for s in supplies
    for p in shipsites
		if value(y[s, p]) > 0 && (s.plant != p)
			push!(outbf_df.MfgSite, s.plant)
            push!(outbf_df.ShipSite, p)
            push!(outbf_df.PType, s.rptype)
			push!(outbf_df.MUnits, value(y[s, p]))
            push!(outbf_df.Freight, value(y[s, p]) * bf_dict[BFLane(s.plant, p)]) 
		end
    end
end

@pipe combine(outbf_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
@pipe combine(outbf_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#aggregate munits by plant
@pipe combine(groupby(outbf_df, :MfgSite), :MUnits => sum => :SumMUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

@pipe combine(groupby(outbf_df, :MfgSite), :MUnits => sum => :SumMUnits, :Freight => sum => :Freight) |> 
    clipboard(sprint(show, "text/tab-separated-values", _))

# out_df[!, :MUnitMiles] = out_df[!, :MUnits] .* out_df[!, :Miles]  # This also works

XLSX.writetable(joinpath(OUT_DIR, "model5TTMCbcbf30Iterations.xlsx"), overwrite = true, 
    collect(DataFrames.eachcol(outbf_df)), DataFrames.names(outbf_df))

# agregate by lane
outplotbf_df = combine(groupby(outbf_df, [:MfgSite, :ShipSite]), :MUnits => sum => :MUnits)
#geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)
#outplot5Gbf_df = innerjoin(outplot5Gbf_df, geocodes_fix, on = :Dest)


na = JSON.parsefile(joinpath(DATA_DIR, "na.topojson"))
d = VegaDatasets.VegaJSONDataset(na, joinpath(DATA_DIR, "na.topojson"))
@vlplot(width=800, height=500) +
@vlplot(
    mark={
        :geoshape,
        fill=:lightgray,
        stroke=:white
    },
    data={
        values=d,
        format={
            type=:topojson,
            feature=:na
        }
        },
    projection={
        type=:albers,
        scale = 800
    }
) + @vlplot(
    mark={:circle,
        fill=:cyan,
        stroke=:black
    },
    data = combine(groupby(outplotob_df, [:Dest, :lat, :lon]), :MUnits => sum => :MUnits),
    projection={
        type=:albers,
        scale = 800
    },
    longitude="lon:q",
    latitude="lat:q",
    size={
        field=:MUnits,
        legend={title="MUnits"},
        scale={type=:quantize,
        nice=:true},
        bin={
            step=100000
        },
    }
) + @vlplot(
    mark={:rule,
        opacity = .6,
        size = 1},
        color=:ShipSite,
        data = outplotob_df,
        transform=[
        {
            lookup=:ShipSite,
            from={
                data=plants_df,
                key=:Name,
                fields=["Latitude", "Longitude"]
            },
            as=["origin_latitude", "origin_longitude"]
        },
        {
            lookup=:Dest,
            from={
                data=geocodes_fix,
                key=:Dest,
                fields=["lat", "lon"]
                },
            as=["dest_latitude", "dest_longitude"]
        }
    ],
    projection={
        type=:albers,
        scale = 800
    },
    longitude="origin_longitude:q",
    latitude="origin_latitude:q",
    longitude2="dest_longitude:q",
    latitude2="dest_latitude:q",
    
)

# BF Freight Plot
@vlplot(width=800, height=500) +
@vlplot(
    mark={
        :geoshape,
        fill=:lightgray,
        stroke=:white
    },
    data={
        values=d,
        format={
            type=:topojson,
            feature=:na
        }
        },
    projection={
        type=:albers,
        scale = 800
    }
)  + @vlplot(
    mark={:rule,
        opacity = .6},
        color = :MfgSite,
        size = :MUnits,
        data=outplot5Gbf_df,
        transform=[
        {
            joinaggregate = [{
                op= :sum,
                field= :MUnits,
                as= :MUnits
            }],
            groupby= ["Dest"]
        },
        {
            lookup=:MfgSite,
            from={
                data=plants_df,
                key=:Name,
                fields=["Latitude", "Longitude"]
            },
            as=["origin_latitude", "origin_longitude"]
        },
        {
            lookup=:ShipSite,
            from={
                data=plants_df,
                key=:Name,
                fields=["Latitude", "Longitude"]
                },
            as=["dest_latitude", "dest_longitude"]
        }
    ],
    projection={
        type=:albers,
        scale = 800
    },
    longitude="origin_longitude:q",
    latitude="origin_latitude:q",
    longitude2="dest_longitude:q",
    latitude2="dest_latitude:q",
    size={
        field=:MUnits,
        legend={title="BF MUnits"},
        scale={type=:quantize,
        nice=:true},
        bin={
            step=10000
        }
    }
)

