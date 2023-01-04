# Model 4 + actual part types by ship site
#   allow BF freight and control the number of ship sites to each
#   customer
# Assume CPU has the same freight penalty as PPD, and then zero out in post processing
#
#using Pkg
#Pkg.add("DataFrames")
using DataFrames # Includes PrettyTables
using Pipe
using Printf
using XLSX
using Statistics
using VegaLite
using VegaDatasets
using JSON3
using JSON

data = JSON3.read("""
{
    "shipsites": ["KENTON",
                 "VISALIA",
                 "SHELBYVILLE",
                 "CLARKSVILLE",
                 "PITTSTON"
                 ] 
 
}
""")

shipsites = collect(values(data.shipsites))
plants = shipsites

include("read_geodata.jl")
using .LoadGeoData
varinfo(LoadGeoData)

const DATA_DIR = joinpath(@__DIR__, "Data")

ship_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_history.xlsx"), "Sheet 1"; infer_eltypes=true)...
)
describe(ship_df) # Note missing values for MfgSite
dropmissing!(ship_df, :MfgSiteName)
describe(ship_df) # Note missing values for MfgSite

# Convert MUnits from Any to Fload64
ship_df[!, "MUnits"] = convert.(Float64, ship_df[!, "MUnits"])

# Filter to just FIVEGUYS Print Group
ship5_df = filter(:PrintGroup => ==("FIVEGUYS"), ship_df)

# Combine the City and State to a Loc
ship5_df.Loc = string.(ship5_df.MSTCity, ", ", ship5_df.MSTState)

# Create a "Lid" or "Cup" Part type in the TypeFix column
# Not needed any more, but we'll add it anyway
ship5_df.PTypeFix = string.(first.(ship5_df.PartType))
ship5_df.PTypeFix = ifelse.(string.(first.(ship5_df.PartType)) .== "L", "LID", "CUP") 
ship5_df.PTypeFix = Symbol.(ship5_df.PTypeFix)
# check that it works
ship5_df[ship5_df[:, :PTypeFix] .== :LID, :]

# Get rid of the PlasticInenuity
ship5_df = ship5_df[ship5_df[:, :ShipSiteName] .!= "PLASTICS INGENUITY", :]

@pipe combine(ship5_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
#dropmissing!(demand_df, :MSTZip)

# Just include where we have distances
ship5_df = innerjoin(ship5_df, DataFrames.select(distw_df, :Dest), on = [:Loc => :Dest])

# Test with two parts and locations
#ship5_df = ship5_df[((ship5_df.PartType .== "SDRA0160") .| (ship5_df.PartType .== "DMR-0320")) .& 
#    ((ship5_df.MSTCity .== "FONTANA") .| (ship5_df.MSTCity .== "COAL TOWNSHIP")), :]

@pipe combine(ship5_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# How much is PPD freight
@pipe combine(filter(:FrtTermsFix => x -> x == "PPD", ship5_df), :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Keep both PPD and CPU in the demand file to get total demand, not just PPD
demand5_df = combine(groupby(ship5_df, [:MSTName, :Loc, :FrtTermsFix, :PartType]), :MUnits => sum => :MUnits)
rename!(demand5_df, :FrtTermsFix => :FTerms)

# Make sure we haven't lost anything
@pipe combine(demand5_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))


# Build the various data structures to hold the supply and demand info

supply5_df = combine(groupby(ship5_df, [:MfgSiteName, :PartType]), :MUnits => sum => :MUnits)
@pipe supply5_df |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

supply5_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "fiveguys_capacity.xlsx"), "Sheet1"; infer_eltypes=true)...
)
dropmissing!(supply5_df)
supply5_df[!, :Capacity] = convert.(Float64, supply5_df[!, :Capacity])
select!(supply5_df, Not([:Shipped]))
subset!(supply5_df, :Capacity => ByRow(!=(0)))
#supply5_df = supply5_df[((supply5_df.PartType .== "SDRA0160") .| (supply5_df.PartType .== "DMR-0320")), :]

# Copy to the clipboard
clipboard(sprint(show, "text/tab-separated-values", supply5_df))
clipboard(sprint(show, "text/tab-separated-values", demand5_df))
clipboard(sprint(show, "text/tab-separated-values", dist_dict))

struct Supply
    plant::String
    ptype::String
end
ptypes = unique(supply5_df.PartType)

#[
#:LID; :CUP
#]

supply_dict = Dict(Supply.(supply5_df.MfgSite, supply5_df.PartType) .=> supply5_df.Capacity)
# Add missing elements with a 0 supply
for p in plants
    for t in ptypes
        #println(Demand(c, p))
        if !(Supply(p, t) in keys(supply_dict)) 
            supply_dict[Supply(p, t)] = 0
        end
    end
end

supply_dict

# Customers as a Named tuple
cust_nt = Tables.rowtable(DataFrames.select(demand5_df, :MSTName, :Loc, :FTerms))
unique(cust_nt)

# Better to use a Struct?
struct Cust
    mstname::String 
    loc::String 
    fterms::String 
end

# populate with the data from the dataframe
customers = Cust.(demand5_df.MSTName, demand5_df.Loc, demand5_df.FTerms)
customers = unique(customers)

eltype(customers)

# build the demand dict using the Cust structure
#demand_dict = Dict(structfromnt.(Cust, Tables.rowtable(DataFrames.select(dmd_df, :MSTName, :Loc, :FTerms))) .=> dmd_df.MUnits)
# Another way to do it
#demand_dict = Dict(Pair.(structfromnt.(Cust, Tables.rowtable(DataFrames.select(dmd_df, :MSTName, :Loc, :FTerms))), dmd_df.MUnits))

# Use this as the key for a demand dict
struct Demand
    customer::Cust
    ptype::String
end

# Solution from Bogumit Kaminski from Julialang@discoursemail.com
demand_key = Demand.(Cust.(demand5_df.MSTName, demand5_df.Loc, demand5_df.FTerms), demand5_df.PartType)

demand_dict = Dict(demand_key .=> demand5_df.MUnits)
keys(demand_dict)
#demand_dict[Demand(Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD"), :LID)]
filter(k -> k.first == Demand(Cust("BEN E. KEITH", "MISSOURI CITY, TX", "PPD"), "DMMS0240"), demand_dict)

# Add missing elements with a 0 demand
for c in customers
    for t in ptypes
        #println(Demand(c, p))
        if !( Demand(c, t) in keys(demand_dict)) 
            demand_dict[Demand(c, t)] = 0
        end
    end
end

# sum the demand values 
@printf("%d", sum(values(demand_dict)))
@pipe [sum(values(demand_dict))] |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
# make sure it's the same
@pipe combine(demand5_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

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

for s in data.shipsites
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

# Build the BF dictionary with costs per M for each O/D Pair
bf_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "bf_lanes.xlsx"), "Sheet1"; infer_eltypes=true)...
)

struct BFLane
    plant::String
    shipsite::String
end

bf_lanes = BFLane.(bf_df.Plant, bf_df.ShipSite)
bf_dict = Dict(bf_lanes .=> bf_df.CostPerM)

#Overwite with large freight penalty to force shipping from the mfg site
for p in plants
    for s in shipsites
        if (p != s) & (s == "KENTON")
            bf_dict[BFLane(p,s)] = 100.0
        end
    end
end


# Scalar f  freight in dollars per case/M per thousand miles 
# Jan - Nov Munits * Miles / 1000 = 7,599,747 [1M Miles * MUnits Shipped] 
# Jan - Nov PPD customer Freight = $42,067,887
# 
f = 5.77 # freight in dollars per M per thousand miles 
 
# Parameter c(l)  transport cost in thousands of dollars per case on each lane ;
#            c(l) = f * dist(l) / 1000 ;
# TODO:  There might be a way to reduce the number of lanes by removing the
#  FTerms component.  Freight will e the same for PPD and CPU
# We first declare an empty dictionary and then we fill it with the values
c = Dict{Lane,Float64}() # transport cost in thousands of dollars per M ;
#[ c[p,c] = f * dist_dict[p,c] / 1000 for p in PLANTS, c in CUST]
[ c[l] = f * dist_dict[l] / 1000 for l in lanes ]
c


using JuMP
using GLPK
model = Model(GLPK.Optimizer)

#MUnits for Plant to Ship site. 
@variable(model, y[bf_lanes, ptypes] >= 0)
#@variable(model5G, t[plants, shipsites, ptypes] >= 0)

# Define variable for Plant to customer
@variable(model, x[lanes, ptypes] >= 0)

# Binary variable if Site S ships to Cust C = 1, 0 otherwise
@variable(model, z[lanes], Bin)

# Binary variable if site S is open
@variable(model, r[shipsites], Bin)

# We need a constraint that each plant can ship no more than its capacity:
@constraint(model,
    [p in plants, t in ptypes], sum(y[BFLane(p, s), t] for s in shipsites) <= 
        supply_dict[Supply(p, t)])

# Customer Demand constraint
@constraint(model,
    [c in customers, t in ptypes], sum(x[Lane(s, c), t] for s in shipsites) == 
        demand_dict[Demand(c, t)])

# shipsite inbound = outbound for net of 0
@constraint(model, [s in shipsites, t in ptypes], sum(y[BFLane(p, s), t] for p in plants) == 
    sum(x[Lane(s, c), t] for c in customers))

#limit the ship sites to 1 per customer
@constraint(model, [c in customers], sum(z[Lane(s, c)] for s in shipsites) == 1)

# Constraint to Limit the number of ship sites
@constraint(model, sum(r[:]) == 2)

# Big Mz constraint for the number of sites servicing each DC
#@constraint(model, [s in shipsites, c in customers, t in ptypes], x[Lane(s, c), t] <= 8000 * z[Lane(s, c)])
@constraint(model, [s in shipsites, c in customers, t in ptypes], x[Lane(s, c), t] <= 
    min(demand_dict[Demand(c, t)]) * z[Lane(s, c)])

@constraint(model, [s in shipsites, c in customers, t in ptypes], x[Lane(s, c), t] <= min(demand_dict[Demand(c, t)]) * r[s])

# Objective
#@objective(model5G, Min, sum(bf_dict[b] * y[b, t] + c[l] * x[l, t] for b in bf_lanes, l in lanes, t in ptypes))
@objective(model, Min, sum(bf_dict[b] * y[b, t] for b in bf_lanes, t in ptypes) + sum(c[l] * x[l, t] for l in lanes, t in ptypes))

optimize!(model)

solution_summary(model)

println("RESULTS:")
for p in plants, s in shipsites, t in ptypes
        println("Plant: $p to Site: $s Product: $t => ", value(y[BFLane(p, s), t]))
end

for s in shipsites, cu in customers, t in ptypes
    if value(x[Lane(s, cu), t]) > 0
        println("Site: $s to $cu Product: $t MUnits =  ", value(x[Lane(s, cu), t]), "  Freight = ",  value(x[Lane(s, cu), t]) * c[Lane(s, cu)])
    end
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
out5Gob_df = DataFrame(ShipSite = String[], Cust = String[], Dest = String[], 
    FTerms = String[], PType = String[], MUnits = Float64[], Freight = Float64[])
for l in lanes
    for p in ptypes
		if !isapprox(value(x[l,p]), 0.0, atol = 1e-10, rtol = 0) #&& l.Cust.FTerms == "PPD"
			push!(out5Gob_df.ShipSite, l.shipsite)
            push!(out5Gob_df.Cust, l.cust.mstname)
            push!(out5Gob_df.Dest, l.cust.loc)
            push!(out5Gob_df.FTerms, l.cust.fterms)
            push!(out5Gob_df.PType, p)
			push!(out5Gob_df.MUnits, value(x[l,p]))
            if l.cust.fterms == "PPD" 
                push!(out5Gob_df.Freight, value(x[l,p]) * c[l])
            else
                push!(out5Gob_df.Freight, 0.0)
            end
		end
    end
end

@pipe combine(out5Gob_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

out5Gob_df = innerjoin(out5Gob_df, dist_df, on = [:ShipSite => :Orig, :Dest])
#transform!(out_df, [:MUnits, :Miles] => ByRow((a, b) -> (a * b * f / 1000)) => :Freight1)

@pipe combine(out5Gob_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#aggregate munits by plant
@pipe combine(groupby(out5Gob_df, :ShipSite), :MUnits => sum => :SumMUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# out_df[!, :MUnitMiles] = out_df[!, :MUnits] .* out_df[!, :Miles]  # This also works

const OUT_DIR = joinpath(@__DIR__, "FiveGuys")

XLSX.writetable(joinpath(OUT_DIR, "model5Gob_output.xlsx"), overwrite = true, 
    collect(DataFrames.eachcol(out5Gob_df)), DataFrames.names(out5Gob_df))

# agregate by lane
outplot5Gob_df = combine(groupby(out5Gob_df, [:ShipSite, :Dest]), :MUnits => sum => :MUnits)
#geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)
outplot5Gob_df = innerjoin(outplot5Gob_df, geocodes_fix, on = :Dest)

# Create the output DataFrame - BF Freight 
out5Gbf_df = DataFrame(MfgSite = String[], ShipSite = String[], 
    PType = String[], MUnits = Float64[], Freight = Float64[])
for b in bf_lanes
    for t in ptypes
		if value(y[b, t]) > 0 && (b.shipsite != b.plant)
			push!(out5Gbf_df.MfgSite, b.plant)
            push!(out5Gbf_df.ShipSite, b.shipsite)
            push!(out5Gbf_df.PType, t)
			push!(out5Gbf_df.MUnits, value(y[b, t]))
            push!(out5Gbf_df.Freight, value(y[b, t]) * bf_dict[b]) 
		end
    end
end

@pipe combine(out5Gbf_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
@pipe combine(out5Gbf_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#aggregate munits by plant
@pipe combine(groupby(out5Gbf_df, [:MfgSite, :ShipSite]), :MUnits => sum => :SumMUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# out_df[!, :MUnitMiles] = out_df[!, :MUnits] .* out_df[!, :Miles]  # This also works

XLSX.writetable(joinpath(OUT_DIR, "model5Gbf_output.xlsx"), overwrite = true, 
    collect(DataFrames.eachcol(out5Gbf_df)), DataFrames.names(out5Gbf_df))

# agregate by lane
outplot5Gbf_df = combine(groupby(out5Gbf_df, [:MfgSite, :ShipSite]), :MUnits => sum => :MUnits)
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
    data = combine(groupby(outplot5Gob_df, [:Dest, :lat, :lon]), :MUnits => sum => :MUnits),
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
            step=5000
        },
    }
) + @vlplot(
    mark={:rule,
        opacity = .6,
        size = 2},
        color=:ShipSite,
        data = outplot5Gob_df,
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

