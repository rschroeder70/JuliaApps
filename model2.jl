# Modele- Add Total Capacity Constraints by plant.  Equivalent to a Single Part.
# Need to deal with CPU and PPD.  The same Ship-To can be both
# Assume CPU has the same freight penalty as PPD, and then zero out in post processing
# Configure the destinations by customer, not just city state as in Model2
#
#using Pkg
#Pkg.add("DataFrames")
using DataFrames # Includes PrettyTables
# using PrettyTables
using Pipe
#using InvertedIndices #Need this for the 'Not' operator
using Printf
using XLSX
using StructArrays
using NamedTupleTools
using VegaLite
using VegaDatasets
using JSON

include("read_geodata.jl")
using .LoadGeoData

const DATA_DIR = joinpath(@__DIR__, "Data")

ship_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_history.xlsx"), "Sheet 1"; infer_eltypes=true)...
)

ship_df.Loc = string.(ship_df.MSTCity, ", ", ship_df.MSTState)

@pipe combine(ship_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
#dropmissing!(demand_df, :MSTZip)

# Keep both PPD and CPU in the demand file
demand_df = combine(groupby(ship_df, [:MSTName, :Loc, :FrtTermsFix]), :MUnits => sum => :MUnits)
rename!(demand_df, :FrtTermsFix => :FTerms)

# Make sure we haven't lost anything
@pipe combine(demand_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Just include where we have distances
dmd_df = innerjoin(DataFrames.select(demand_df, :MSTName, :Loc, :FTerms, :MUnits), DataFrames.select(distw_df, :Dest), on = [:Loc => :Dest])

# Named tuple
cust_nt = Tables.rowtable(DataFrames.select(dmd_df, :MSTName, :Loc, :FTerms))

# Better to use a Struct?
struct Cust
    MSTName::String 
    Loc::String 
    FTerms::String 
end

# convert the named tuple to a struct
cust_st = structfromnt.(Cust, cust_nt)
cust_st[1].MSTName

# build the demand dict using the Cust structure
demand_dict = Dict(structfromnt.(Cust, Tables.rowtable(DataFrames.select(dmd_df, :MSTName, :Loc, :FTerms))) .=> dmd_df.MUnits)
# Another way to do it
#demand_dict = Dict(Pair.(structfromnt.(Cust, Tables.rowtable(DataFrames.select(dmd_df, :MSTName, :Loc, :FTerms))), dmd_df.MUnits))

# sum the demand values - not sure how to format with commas
@printf("%d", sum(values(demand_dict)))
@pipe [sum(values(demand_dict))] |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
# make sure it's the same
@pipe combine(dmd_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# build the Set of Destinations
# DEST = dmd_df.Dest
CUST = cust_st
eltype(CUST)

# Set of Origins
PLANTS = DataFrames.names(distw_df)
PLANTS = PLANTS[2:6]
#ORIG = unique(dist_df.Orig)

# Build the lanes (PLANT, DEST pairs)
struct Lane
    Plant::String
    Cust::Cust
end

lanes = Lane[]

for p in PLANTS
    for c in CUST
        push!(lanes, Lane(p, c))
    end
end

length(lanes)

# Then make a dict Lane => dist
# use 'only' function in the next step to return miles as a scalar, not a 1-d vector
dist_dict = Dict(l => only(dist_df[(dist_df[:, :Orig].== l.Plant) .& (dist_df[:, :Dest].== l.Cust.Loc), :Miles]) for l in lanes)
values(dist_dict)


lanes[1]
dist_dict[Lane("KENTON", Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD"))]
dist_dict[lanes[1]]


# Scalar f  freight in dollars per case/M per thousand miles  /90/ ;
f = 5.64 # freight in dollars per M per thousand miles 
 
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

cap_dict = Dict(PLANTS .=> [3000000.0, 1500000.0, 8000000.0, 3500000.0, 1800000.0])

using JuMP
import Test
import GLPK

model2 = Model(GLPK.Optimizer)

# Define variables
@variable(model2, x[lanes] >= 0)
has_lower_bound(x[Lane("KENTON", Cust("ACME PAPER & SUPPLY, INC.", "JESSUP, MD", "PPD"))])
lower_bound.(x)
model2[:x]


# Demand constraint
@constraints model2 begin 
    dmd[c in CUST], sum(x[Lane(p, c)] for p in PLANTS) == demand_dict[c]
end

# Supply Constraint
@constraints model2 begin
    sply[p in PLANTS], sum(x[Lane(p, c)] for c in CUST) <= cap_dict[p]
end

model2[:dmd]
model2[:sply]

# Objective
@objective(model2, Min, sum(c[l] * x[l] for l in lanes))

optimize!(model2)

#print(model2_output)

status = termination_status(model2)
if (status == MOI.OPTIMAL) || (status == MOI.LOCALLY_SOLVED)
    println("Objective value: ",objective_value(model2))
    println(value.(x))
    #println("Shadow prices of supply:")
    #[println("$p = $(dual(supply[p]))") for p in plants]
    println("Shadow prices of demand:")
    [println("$c = $(dual(dmd[c]))") for c in CUST]
else
    println("Model didn't solved")
    println(status)
end

Test.@test termination_status(model2) == MOI.OPTIMAL
Test.@test primal_status(model2) == MOI.FEASIBLE_POINT
println("the objective is: ")
println(objective_value(model2))
@printf "%d" objective_value(model2)
@pipe [objective_value(model2)] |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# # This is a test for the output
# for l in lanes
#         println(" $(l) = $(value(x[l]))")
#         println(l.Cust.FTerms)
# end

lanes[1].Cust.FTerms

println("The optimal solution is: ")
println(value.(x))

# Create the output DataFrame
out_df = DataFrame(Plant = String[], Cust = String[], Dest = String[], FTerms = String[], MUnits = Float64[], Freight = Float64[])
for l in lanes
		if value(x[l]) > 0 #&& l.Cust.FTerms == "PPD"
			push!(out_df.Plant, l.Plant)
            push!(out_df.Cust, l.Cust.MSTName)
            push!(out_df.Dest, l.Cust.Loc)
            push!(out_df.FTerms, l.Cust.FTerms)
			push!(out_df.MUnits, value(x[l]))
            if l.Cust.FTerms == "PPD" 
                push!(out_df.Freight, value(x[l]) * c[l])
            else
                push!(out_df.Freight, 0.0)
            end
		end
end

@pipe combine(out_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

out_df = innerjoin(out_df, dist_df, on = [:Plant => :Orig, :Dest])
#transform!(out_df, [:MUnits, :Miles] => ByRow((a, b) -> (a * b * f / 1000)) => :Freight1)

@pipe combine(out_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#aggregate munits by plant
@pipe combine(groupby(out_df, :Plant), :MUnits => sum => :SumMUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# out_df[!, :MUnitMiles] = out_df[!, :MUnits] .* out_df[!, :Miles]  # This also works

XLSX.writetable("model2_output.xlsx", overwrite = true, collect(DataFrames.eachcol(out_df)), DataFrames.names(out_df))

outplot_df = innerjoin(out_df, geocodes_df, on = :Dest)

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
    data=out_df,
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
            step=50000
        },
    }
) + @vlplot(
    mark={:rule},
        color=:Plant,
    data=out_df,
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
                data=geocodes_df,
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
