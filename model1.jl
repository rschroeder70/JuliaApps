# Model2 - No capacity constraints.  Single Part.
#using Pkg
#Pkg.add("DataFrames")
using DataFrames
using Pipe
using InvertedIndices #Need this for the 'Not' operator
using Printf
using JSON
using VegaLite
using VegaDatasets

#Pkg.add("XLSX")
using XLSX

include("read_geodata.jl")
using .LoadGeoData

# will try to do this with a dataframe
#dist_m = Matrix(dist_df[:, 2:6])
#parse.(Float64, string.(dist))

#XLSX.writetable("temp.xlsx", collect(DataFrames.eachcol(demandf_df)), DataFrames.names(demandf_df))

#Create a dataframe of the demand for PPD shipments to each Destinations
const DATA_DIR = joinpath(@__DIR__, "Data")

ship_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_history.xlsx"), "Sheet 1")...
)

ship_df.Dest = string.(ship_df.MSTCity, ", ", ship_df.MSTState)

@pipe combine(ship_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
#dropmissing!(demand_df, :MSTZip)

names(ship_df)

# Keep both PPD and CPU in the demand file
demand_df = combine(groupby(ship_df, [:Dest, :FrtTermsFix]), :MUnits => sum => :MUnits)

# Make sure we haven't lost anything
@pipe combine(demand_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Total PPD Demand
ppddmd_df = demand_df[demand_df.FrtTermsFix .== "PPD", :]
@pipe combine(ppddmd_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Just include where we have distances
ppddmd_df = innerjoin(DataFrames.select(ppddmd_df, :Dest, :MUnits), DataFrames.select(distw_df, :Dest), on = :Dest)
@pipe combine(ppddmd_df, :MUnits => sum => :MUnits) |> PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# build the Set of Destinations from the distance table
DEST = ppddmd_df.Dest
#eltype.(Dest)

# Set of Origins
PLANTS = DataFrames.names(distw_df)
PLANTS = PLANTS[2:6]
#ORIG = unique(dist_df.Orig)

# Convert the demand data frame to a Dict
demand = Dict(Pair.(ppddmd_df.Dest, ppddmd_df.MUnits))

# only need distances for cities in the dmdppd_df
distwppd_df = innerjoin(distw_df, ppddmd_df, on = :Dest)

# transpose the distance dataframe. The first row is a duplicate of the column headings - don't want as a row
#distw_df = DataFrame([[names(distw_df)]; collect.(eachrow(distw_df))], [:column; Symbol.(axes(distw_df, 1))])
distwppdt_df = DataFrame([[names(distwppd_df)]; collect.(eachrow(distwppd_df))], [:column; Symbol.(distwppd_df[:, 1])])[Not([1]), :]
distwppdt_df[!, 2:end] = convert.(Float64, distwppdt_df[:, 2:end])
distwppdt_df[!, 1] = convert.(String, distwppdt_df[:, 1])
distwppdt_df = rename(distwppdt_df, :column => :PLANTS)

dist = Dict((r[:PLANTS], d) => r[Symbol(d)] for r in eachrow(distwppdt_df), d in DEST)
# Here we are converting the table in a "(plant, dest) => distance" dictionary
# r[:plants]:   the first key, row field using a fixed header
# m:            the second key
# r[Symbol(m)]: the value, the row field with a dynamic header

# check
dist["VISALIA", "MILFORD, MA"]

# Scalar f  freight in dollars per case/M per thousand miles  /90/ ;
f = 5.64 # freight in dollars per M per thousand miles 
 
# Parameter c(i,j)  transport cost in thousands of dollars per case ;
#            c(i,j) = f * d(i,j) / 1000 ;
# We first declare an empty dictionary and then we fill it with the values
c = Dict() # transport cost in thousands of dollars per M ;
[ c[p,d] = f * dist[p,d] / 1000 for p in PLANTS, d in DEST]

using JuMP
import Test
import GLPK

model1 = Model(GLPK.Optimizer)

# Define variables ##
#  Variables
#       x(i,j)  shipment quantities in MUnits
#       z       total transportation costs in thousands of dollars ;
#  Positive Variable x ;
@variable(model1, x[p in PLANTS, d in DEST] >= 0)

## Define contrains ##
# supply(i)   observe supply limit at plant i
# supply(i) .. sum (j, x(i,j)) =l= a(i)
# demand(j)   satisfy demand at dest j ;  
# demand(j) .. sum(i, x(i,j)) =g= b(j);
@constraints model1 begin
    #supply[p in plants],   # observe supply limit at plant p
    #    sum(x[p,m] for m in markets)  <=  a[p]
    dmd[d in DEST],  # satisfy demand at dest d
        sum(x[p,d] for p in PLANTS)   >=  demand[d]
end

# Objective
@objective model1 Min begin
    sum(c[p,d]*x[p,d] for p in PLANTS, d in DEST)
end

print(model1)

optimize!(model1)

status = termination_status(model1)
if (status == MOI.OPTIMAL) || (status == MOI.LOCALLY_SOLVED)
    println("Objective value: ",objective_value(model1))
    println(value.(x))
    #println("Shadow prices of supply:")
    #[println("$p = $(dual(supply[p]))") for p in plants]
    println("Shadow prices of demand:")
    [println("$d = $(dual(dmd[d]))") for d in DEST]
else
    println("Model didn't solved")
    println(status)
end

Test.@test termination_status(model1) == MOI.OPTIMAL
Test.@test primal_status(model1) == MOI.FEASIBLE_POINT
println("the objective is: ")
println(objective_value(model1))
@printf "%d" objective_value(model1)

for p in PLANTS
    for d in DEST
        println(" $(p) $(d) = $(value(x[p, d]))")
    end
end

println("The optimal solution is: ")
println(value.(x))

# Create the output DataFrame
out_df = DataFrame(Orig = String[], Dest = String[], MUnits = Float64[])
for p in PLANTS
	for d in DEST
		if value(x[p,d]) > 0 	
			push!(out_df.Orig, p)
			push!(out_df.Dest, d)
			push!(out_df.MUnits, value(x[p,d]))
		end
	end
end

@pipe combine(out_df, :MUnits => sum => :MUnits) |> pretty_table(_; formatters = ft_printf("%'d"))

out_df = innerjoin(out_df, dist_df, on = [:Orig, :Dest])
transform!(out_df, [:MUnits, :Miles] => ByRow((a, b) -> (a * b * f / 1000)) => :Freight)

@pipe combine(out_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> pretty_table(_; formatters = ft_printf("%'d"))

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
