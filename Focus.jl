using DataFrames
using Pipe
using ODBC
using Printf
using XLSX
using Statistics
using VegaLite
using VegaDatasets
using JSON
using JSON3

const DATA_DIR = abspath(joinpath(joinpath(@__DIR__, "..", "GPI Foodservice", "Projects", "Focus Brands")))
const GEO_DIR = joinpath(@__DIR__, "Data")

include("read_geodata.jl")
using .LoadGeoData
varinfo(LoadGeoData)

#Get the historical shipments to determine the ship-to location splits by part type
supply_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "SupplyDemandData.xlsx"), "Source"; infer_eltypes=true)...
)
describe(supply_df) # Note missing values for MfgSite
dropmissing!(supply_df, :PrintGroup)

# Two ways of creating a new column.  Add 10 to make supply > demand
#transform(supply_df, :MMUnits =>  ByRow(x -> x * 1000) => :MUnits)
supply_df.MUnits = supply_df.MMUnits .* 1000 .+ 10
clipboard(sprint(show, "text/tab-separated-values", supply_df))

@pipe combine(supply_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))


# Read the demand info
dmd_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "SupplyDemandData.xlsx"), "Raw Demand"; infer_eltypes=true)...
)

describe(dmd_df)
dmd_df = stack(dmd_df, 7:84, variable_name=:MST, value_name=:MUnits)
select!(dmd_df, [:ITEM, :Print, :MST, :MUnits])

clipboard(sprint(show, "text/tab-separated-values", dmd_df))

dmd_df.MST = SubString.(dmd_df.MST, 1, 8)
dropmissing!(dmd_df, :ITEM)
dmd_df = dmd_df[dmd_df.ITEM .!= "ITEM", :]
# clipboard(sprint(show, "text/tab-separated-values", dmd_df))

dmd_df.MUnits = parse.(Float64, string.(dmd_df.MUnits))
dmd_df = dmd_df[dmd_df.MUnits .> 0.0, :]

#Agregate the Lineage DC in Wilmington, Includes
dmd_df = combine(groupby(dmd_df, [:Print, :ITEM, :MST]), :MUnits => sum => :MUnits)

@pipe combine(dmd_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

#Get the MST Locations
conn = ODBC.Connection("CUPS", "crystal", "Bobject\$")
query = """ 
    SELECT DISTINCT [MST]
    ,[MSTName]
    --,[MSTSetup]
    ,[MSTZip]
    ,[MSTLine1]
    ,[MSTLine2]
    ,[MSTCity]
    ,[MSTState]
    ,[MSTCountry]
    FROM [LIDS].[dbo].[MasterShipTo]
    WHERE MST in ('98034627','98015247','98033229','98013040','98026320','98036225','98035410','98033630',
    '98026924','98033795','98018011','98030490','98033836','98030119','98032733','98034700','98028496',
    '95000320','59257002','98014527','98032192','98012210','98035724','98000882','98000971','79470679', 
    '98020147')
"""
mst_df = DBInterface.execute(conn, query) |> DataFrames.DataFrame

DBInterface.close!(conn)

clipboard(sprint(show, "text/tab-separated-values", mst_df))

dmd_df = leftjoin(dmd_df, select(mst_df, [:MST, :MSTName, :MSTCity, :MSTState, :MSTZip]), on = [:MST])

clipboard(sprint(show, "text/tab-separated-values", dmd_df))

# Combine the City and State to a Loc
dmd_df.Loc = string.(dmd_df.MSTCity, ", ", dmd_df.MSTState)
dmd_df.FTerms .= "PPD"

#Aggregate to ensure unique keys
dmd_df = combine(groupby(dmd_df, [:MSTName, :Loc, :ITEM, :FTerms]), :MUnits => sum => :MUnits)

# Build the customer structure
# Better to use a Struct?
struct Cust
    mstname::String 
    loc::String 
    fterms::String 
end

# populate with the data from the dataframe
customers = Cust.(dmd_df.MSTName, dmd_df.Loc, dmd_df.FTerms)
customers = unique(customers)

struct Demand
    customer::Cust
    ptype::String
end

# Solution from Bogumit Kaminski from Julialang@discoursemail.com
dmd_dict = Dict(Demand.(Cust.(dmd_df.MSTName, dmd_df.Loc, dmd_df.FTerms), dmd_df.ITEM) .=> dmd_df.MUnits)
@printf("%d", sum(values(dmd_dict)))

ptypes = unique(supply_df.PartType)

# Add missing elements with a 0 demand
# Can I skip this and use the get() function with a zero default value
for c in customers
    for t in ptypes
        #println(Demand(c, p))
        #if !( Demand(c, t) in keys(demand_dict)) 
        #    demand_dict[Demand(c, t)] = 0
        #end
        dmd_dict[Demand(c, t)] = get(dmd_dict, Demand(c, t), 0)
    end
end

# sum the demand values 
@printf("%d", sum(values(dmd_dict)))


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

#Aggregate
supply_df = combine(groupby(supply_df, [:MfgSiteName, :PartType]), :MUnits => sum => :Capacity)
rename!(supply_df, :MfgSiteName => :MfgSite)


@pipe combine(supply_df, :Capacity => sum => :Capacity) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

struct Supply
    plant::String
    ptype::String
end

#[
#:LID; :CUP
#]

supply_dict = Dict(Supply.(supply_df.MfgSite, supply_df.PartType) .=> supply_df.Capacity)
@printf("%d", sum(values(supply_dict)))

# Add missing elements with a 0 supply
for p in plants
    for t in ptypes
        #println(Demand(c, p))
        # if !(Supply(p, t) in keys(supply_dict)) 
        #     supply_dict[Supply(p, t)] = 0
        # end
        supply_dict[Supply(p, t)] = get(supply_dict, Supply(p, t), 0)
    end
end

supply_dict

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
    XLSX.readtable(joinpath(GEO_DIR, "bf_lanes.xlsx"), "Sheet1"; infer_eltypes=true)...
)

struct BFLane
    plant::String
    shipsite::String
end

bf_lanes = BFLane.(bf_df.Plant, bf_df.ShipSite)
bf_dict = Dict(bf_lanes .=> bf_df.CostPerMUnit)

bf_dict[BFLane("CLARKSVILLE", "SHELBYVILLE")]

f = 6.15 # freight in dollars per M per thousand miles 
 
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
using HiGHS
import Test
model = Model(HiGHS.Optimizer)

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
        dmd_dict[Demand(c, t)])

# shipsite inbound = outbound for net of 0
@constraint(model, [s in shipsites, t in ptypes], sum(y[BFLane(p, s), t] for p in plants) == 
    sum(x[Lane(s, c), t] for c in customers))

#limit the ship sites to 1 per customer
@constraint(model, [c in customers], sum(z[Lane(s, c)] for s in shipsites) == 1)

# Constraint to Limit the number of ship sites
@constraint(model, sum(r[:]) <= 5)
@constraint(model, r["VISALIA"] == 0)
@constraint(model, r["PITTSTON"] == 0)
@constraint(model, r["CLARKSVILLE"] == 0)

# Big Mz constraint for the number of sites servicing each DC
#@constraint(model, [s in shipsites, c in customers, t in ptypes], x[Lane(s, c), t] <= 8000 * z[Lane(s, c)])
@constraint(model, [s in shipsites, c in customers, t in ptypes], x[Lane(s, c), t] <= 
    min(dmd_dict[Demand(c, t)]) * z[Lane(s, c)])

#Big M for the total number of sites
@constraint(model, [s in shipsites, c in customers, t in ptypes], x[Lane(s, c), t] <= min(dmd_dict[Demand(c, t)]) * r[s])

# Objective
#@objective(model5G, Min, sum(bf_dict[b] * y[b, t] + c[l] * x[l, t] for b in bf_lanes, l in lanes, t in ptypes))
@objective(model, Min, sum(bf_dict[b] * y[b, t] for b in bf_lanes, t in ptypes) + sum(c[l] * x[l, t] for l in lanes, t in ptypes))

optimize!(model)

Test.@test termination_status(model) == OPTIMAL

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
outob_df = DataFrame(ShipSite = String[], Cust = String[], Dest = String[], 
    FTerms = String[], PType = String[], MUnits = Float64[], Freight = Float64[])
for l in lanes
    for p in ptypes
		#if value(x[l,p]) > 0 #&& l.Cust.FTerms == "PPD"
        if !isapprox(value(x[l,p]), 0.0, atol = 1e-10, rtol = 0) #&& l.Cust.FTerms == "PPD"
			push!(outob_df.ShipSite, l.shipsite)
            push!(outob_df.Cust, l.cust.mstname)
            push!(outob_df.Dest, l.cust.loc)
            push!(outob_df.FTerms, l.cust.fterms)
            push!(outob_df.PType, p)
			push!(outob_df.MUnits, value(x[l,p]))
            if l.cust.fterms == "PPD" 
                push!(outob_df.Freight, value(x[l,p]) * c[l])
            else
                push!(outob_df.Freight, 0.0)
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

#const OUT_DIR = joinpath(@__DIR__, "Data", "BK")
const OUT_DIR = abspath(joinpath(joinpath(@__DIR__, "..", "GPI Foodservice", "Projects", "Focus Brands")))

XLSX.writetable(joinpath(OUT_DIR, "model4ob_output.xlsx"), overwrite = true, 
    collect(DataFrames.eachcol(outob_df)), DataFrames.names(outob_df))

# agregate by lane
outplotob_df = combine(groupby(outob_df, [:ShipSite, :Dest]), :MUnits => sum => :MUnits)
#geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)
outplotob_df = innerjoin(outplotob_df, geocodes_fix, on = :Dest)

# Create the output DataFrame - BF Freight 
outbf_df = DataFrame(MfgSite = String[], ShipSite = String[], 
    PType = String[], MUnits = Float64[], Freight = Float64[])
for b in bf_lanes
    for t in ptypes
		if value(y[b, t]) > 0 && (b.shipsite != b.plant)
			push!(outbf_df.MfgSite, b.plant)
            push!(outbf_df.ShipSite, b.shipsite)
            push!(outbf_df.PType, t)
			push!(outbf_df.MUnits, value(y[b, t]))
            push!(outbf_df.Freight, value(y[b, t]) * bf_dict[b]) 
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

XLSX.writetable(joinpath(OUT_DIR, "model4bf_output.xlsx"), overwrite = true, 
    collect(DataFrames.eachcol(outbf_df)), DataFrames.names(outbf_df))

# agregate by lane
outplotbf_df = combine(groupby(outbf_df, [:MfgSite, :ShipSite]), :MUnits => sum => :MUnits)
#geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)
#outplot5Gbf_df = innerjoin(outplot5Gbf_df, geocodes_fix, on = :Dest)


na = JSON.parsefile(joinpath(GEO_DIR, "na.topojson"))
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
            step=50000
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

