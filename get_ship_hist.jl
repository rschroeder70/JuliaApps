using ODBC
using TableView
using DataFrames
using XLSX
using JLD  # to salve the workspace
using VegaLite
using VegaDatasets
using JSON

include("read_distances.jl")

# Need to first load the geocode info from the read_distances.jl file

const DATA_DIR = joinpath(@__DIR__, "Data")

conn = ODBC.Connection("CUPS", "crystal", "Bobject\$")
#df = DBInterface.execute(conn, "SELECT TOP (100) * FROM [LIDS].[dbo].[FinanceForecastV1]") |> DataFrames.DataFrame
#query = """ SELECT TOP (100) * FROM [LIDS].[dbo].[FinanceForecastV1];"""
#df = DBInterface.execute(conn, query) |> DataFrames.DataFrame

frt_sql = raw"C:\Users\SchroedR\OneDrive - Graphic Packaging International, Inc\Documents\MyDocuments\GPI Foodservice\SQL\Freight\ship_history_from_idh_r.sql"
#sql_file = joinpath(homedir(), "OneDrive - Graphic Packaging International, Inc", "Documents", "MyDocuments", "GPI Foodservice", "SQL", "Freight", "ship_history_from_idh_r.sql")

f = open("ship_history_from_idh_r.sql", "r")
                                 
f = open(frt_sql, "r")

sql = read(f, String)

ship_df = DBInterface.execute(conn, sql) |> DataFrames.DataFrame

#VSCodeServer.vscodedisplay(ship_df)
DataFrames.describe(ship_df)

demand_df = DataFrames.combine(DataFrames.groupby(ship_df, [:MSTCity, :MSTState, :FrtTermsFix]), :MUnits => sum, renamecols = false)

demand_df.MUnits = Float64.(demand_df.MUnits)

# remove EX rows and PR
filter!(:MSTState => x -> !(x in(["EX", "PR", "HI"])), demand_df)
filter!(:MSTCity => x -> x != "GNANGARA", demand_df)

demand_df.Dest = string.(demand_df.MSTCity, ", ", demand_df.MSTState)

demand_df = leftjoin(demand_df, geocodes_df, on = [:MSTCity, :MSTState])

dropmissing!(demand_df, :MSTZip)

#DataFrames.select!(demand_df, :Dest, :MUnits)

#VSCodeServer.vscodedisplay(demand_df)

# Filter out the CPU for the baseline dataframe
baseline_df = DataFrames.filter(:FrtTermsFix => !=("CPU"), ship_df)
baseline_df = DataFrames.combine(DataFrames.groupby(baseline_df, [:ShipSiteName, :MSTCity, :MSTState, :FrtTermsFix]), :MUnits => sum, renamecols = false)
DataFrames.describe(baseline_df)
baseline_df.MUnits = Float64.(baseline_df.MUnits)
baseline_df[ismissing.(baseline_df.ShipSiteName) .| ismissing.(baseline_df.MSTCity) .| ismissing.(baseline_df.MSTState), :]
dropmissing!(baseline_df, disallowmissing = true)

baseline_df.Dest = string.(baseline_df.MSTCity, ", ", baseline_df.MSTState)

#DataFrames.select!(baseline_df, :ShipSiteName, :Dest, :MUnits)

baseline_df = leftjoin(baseline_df, dist_df, on = [:ShipSiteName => :Orig, :Dest])
# TODO: show what is missing
dropmissing!(baseline_df, disallowmissing = true)

transform!(baseline_df, [:MUnits, :Miles] => ByRow((a, b) -> (a * b)) => :MUnitMiles)

# out_df[!, :MUnitMiles] = out_df[!, :MUnits] .* out_df[!, :Miles]  # This also works

XLSX.writetable("baseline.xlsx", overwrite = true, collect(DataFrames.eachcol(baseline_df)), DataFrames.names(baseline_df))

#@save "model2.jld"

# VSCodeServer.vscodedisplay(baseline_df)


# Read na topojson data map created in R (US, Canada & Mexico) plot PPD demand on it
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
    data=demand_df,
    transform=[
        {filter="datum.FrtTermsFix != 'CPU'"}
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
        color=:ShipSiteName,
    data=baseline_df,
    transform=[
    #{filter={field=:ShipSiteName ,equal=:SHELBYVILLE}},
        {
            lookup=:ShipSiteName,
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
                data=demand_df,
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

