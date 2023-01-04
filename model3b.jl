# Model 3 + actual part types by ship site
# Assume CPU has the same freight penalty as PPD, and then zero out in post processing
#
#using Pkg
#Pkg.add("DataFrames")
using DataFrames # Includes PrettyTables
using Pipe
using Printf
using XLSX
#using StructArrays
using Statistics
#using NamedTupleTools
using VegaLite
using VegaDatasets
using JSON

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

# These do the same thing
ship_df.Loc = string.(ship_df.MSTCity, ", ", ship_df.MSTState)
ship_df.Loc = ship_df.MSTCity .* ", " .* ship_df.MSTState

# Create a "Lid" or "Cup" Part type in the TypeFix column
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

# Just the PPD freight
@pipe combine(filter(:FrtTermsFix => x -> x == "PPD", ship_df), :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Keep both PPD and CPU in the demand file to get total demand, not just PPD
demand_df = combine(groupby(ship_df, [:MSTName, :Loc, :FrtTermsFix, :PartType]), :MUnits => sum => :MUnits)
rename!(demand_df, :FrtTermsFix => :FTerms)

# Make sure we haven't lost anything
@pipe combine(demand_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Build the supply table
#ship_df.PartType = Symbol.(ship_df.PartType)
supply_df = combine(groupby(ship_df, [:MfgSiteName, :PartType]), :MUnits => sum => :MUnits)
@pipe supply_df |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

struct Supply
    plant::String
    ptype::String
end
ptypes = unique(supply_df.PartType)
plants = unique(supply_df.MfgSiteName)
#[
#:LID; :CUP
#]

supply_dict = Dict(Supply.(supply_df.MfgSiteName, supply_df.PartType) .=> supply_df.MUnits)
# Add missing elements with a 0 supply
for s in plants
    for t in ptypes
        #println(Demand(c, p))
        if !(Supply(s, t) in keys(supply_dict)) 
            supply_dict[Supply(s, t)] = 0
        end
    end
end

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
    ptype::String
end

# Solution from Bogumit Kaminski from Julialang@discoursemail.com
demand_key = Demand.(Cust.(demand_df.MSTName, demand_df.Loc, demand_df.FTerms), demand_df.PartType)

demand_dict = Dict(demand_key .=> demand_df.MUnits)
keys(demand_dict)
#demand_dict[Demand(Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD"), :LID)]

# Add missing elements with a 0 demand
for c in customers
    for p in ptypes
        #println(Demand(c, p))
        if !( Demand(c, p) in keys(demand_dict)) 
            demand_dict[Demand(c, p)] = 0
        end
    end
end

# sum the demand values 
@printf("%d", sum(values(demand_dict)))
@pipe [sum(values(demand_dict))] |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
# make sure it's the same
@pipe combine(demand_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Set of Origins
# plants = String[]
# for (i, j) in supply_dict
#     push!(plants, i.plant)
# end
# plants = unique(plants)

# Build the lanes (PLANT, DEST pairs)
struct Lane
    Plant::String
    Cust::Cust
end

lanes = Lane[]

for p in plants
    for c in customers
        push!(lanes, Lane(p, c))
    end
end

length(lanes)

# Then make a dict Lane => dist
# use 'only' function in the next step to return miles as a scalar, not a 1-d vector
dist_dict = Dict(l => only(dist_df[(dist_df[:, :Orig].== l.Plant) .& 
    (dist_df[:, :Dest].== l.Cust.loc), :Miles]) for l in lanes)

#values(dist_dict)

lanes[1]
dist_dict[Lane("KENTON", Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD"))]
dist_dict[lanes[1]]

# Scalar f  freight in dollars per case/M per thousand miles 
# Jan - Nov Munits * Miles / 1000 = 7,599,747 [1M Miles * MUnits Shipped] 
# Jan - Nov PPD customer Freight = $42,067,887
# 
f = 5.565 # freight in dollars per M per thousand miles 
 
# Parameter c(l)  transport cost in thousands of dollars per case on each lane ;
#            c(l) = f * dist(l) / 1000 ;
# TODO:  There might be a way to reduce the number of lanes by removing the
#  FTerms component.  Freight will e the same for PPD and CPU
# We first declare an empty dictionary and then we fill it with the values
c = Dict{Lane,Float64}() # transport cost in thousands of dollars per M ;
#[ c[p,c] = f * dist_dict[p,c] / 1000 for p in PLANTS, c in CUST]
[ c[l] = f * dist_dict[l] / 1000 for l in lanes ]
c
# Check one
dist_dict[Lane("SHELBYVILLE", Cust("MILLENNIUM MIDWEST", "LA GRANGE PARK, IL", "PPD"))]

#cap_dict = Dict(PLANTS .=> [3000000.0, 1500000.0, 8000000.0, 3500000.0, 1800000.0])

#----------------------------------------------
using JuMP
import Test
import GLPK

function solve_trans_1(lanes::Vector{Lane}, 
                        ptypes::Vector, 
                        customers::Vector{Cust}, 
                        plants::Vector)
    model3b = Model(GLPK.Optimizer)
    set_silent(model3b)
    
    # Define variables
    @variable(model3b, x[lanes, ptypes] >= 0)
    
    # Demand constraint
    @constraints model3b begin 
        dmd[c in customers, t in ptypes], sum(x[Lane(p, c), t] for p in plants) == 
           demand_dict[Demand(c, t)]
    end
    
    #Supply Constraint
    @constraints model3b begin
        sply[p in plants, t in ptypes], sum(x[Lane(p, c), t] for c in customers) <= 
            supply_dict[Supply(p, t)]
    end
    
    # Objective
    @objective(model3b, Min, sum(c[l] * x[l, t] for l in lanes, t in ptypes))
    
    optimize!(model3b)

    return value.(x)
end

res = solve_trans_1(lanes, ptypes, customers, plants)

res
value.(res)

#has_lower_bound(x[Lane("KENTON", Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD")), :LID])
#lower_bound.(x)
#model3b[:x]

#x[Lane("KENTON", Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD")), :LID]

model3b[:dmd]
model3b[:sply]


#print(model2_output)

status = termination_status(model3b)

if (status == MOI.OPTIMAL) || (status == MOI.LOCALLY_SOLVED)
    println("Objective value: ",objective_value(model3b))
    println(value.(x))
    #println("Shadow prices of supply:")
    #[println("$p = $(dual(supply[p]))") for p in plants]
    #println("Shadow prices of demand:")
    #[println("$c = $(dual(dmd[c]))") for c in customers]
else
    println("Model didn't solved")
    println(status)
end

Test.@test termination_status(model3b) == MOI.OPTIMAL
Test.@test primal_status(model3b) == MOI.FEASIBLE_POINT
println("the objective is: ")
println(objective_value(model3b))
@printf "%d" objective_value(model3b)
@pipe [objective_value(model3b)] |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# # This is a test for the output
# for l in lanes
#         println(" $(l) = $(value(x[l]))")
#         println(l.Cust.FTerms)
# end

# Create the output DataFrame
out3b_df = DataFrame(Plant = String[], Cust = String[], Dest = String[], 
    FTerms = String[], PType = String[], MUnits = Float64[], Freight = Float64[])
for l in lanes
    for p in ptypes
		if value(res[l,p]) > 0 #&& l.Cust.FTerms == "PPD"
			push!(out3b_df.Plant, l.Plant)
            push!(out3b_df.Cust, l.Cust.mstname)
            push!(out3b_df.Dest, l.Cust.loc)
            push!(out3b_df.FTerms, l.Cust.fterms)
            push!(out3b_df.PType, p)
			push!(out3b_df.MUnits, value(res[l,p]))
            if l.Cust.fterms == "PPD" 
                push!(out3b_df.Freight, value(res[l,p]) * c[l])
            else
                push!(out3b_df.Freight, 0.0)
            end
		end
    end
end

@pipe combine(out3b_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

out3b_df = innerjoin(out3b_df, dist_df, on = [:Plant => :Orig, :Dest])
#transform!(out_df, [:MUnits, :Miles] => ByRow((a, b) -> (a * b * f / 1000)) => :Freight1)

@pipe combine(out3b_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#aggregate munits by plant
@pipe combine(groupby(out3b_df, :Plant), :MUnits => sum => :SumMUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# out_df[!, :MUnitMiles] = out_df[!, :MUnits] .* out_df[!, :Miles]  # This also works

XLSX.writetable("model3b_output.xlsx", overwrite = true, 
    collect(DataFrames.eachcol(out3b_df)), DataFrames.names(out3b_df))

# agregate by lane
outplot3b_df = combine(groupby(out3b_df, [:Plant, :Dest]), :MUnits => sum => :MUnits)
#geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)
outplot3b_df = innerjoin(outplot3b_df, geocodes_fix, on = :Dest)

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
    data=outplot3b_df,
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
            lookup = "Dest",
            from = {
                data = geocodes_df,
                key = "Dest",
                fields=["lat", "lon"]
            },
            as= ["lat", "lon"]
        }
    ],
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
            step=70000
        },
    }
) + @vlplot(
    mark={:rule,
        opacity = .5},
        color=:Plant,
        data=outplot3b_df,
        transform=[
    #{filter={field=:ShipSiteName ,equal=:SHELBYVILLE}},
        {
            lookup=:Plant,
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
