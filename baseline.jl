#using ODBC
using TableView
using DataFrames
using XLSX
#using JLD  # to save the workspace
using VegaLite
using VegaDatasets
using JSON
using Pipe
using PrettyTables
using Printf

# Need to first load the geocode info from the read_geodata.jl file
include("read_geodata.jl")
using .LoadGeoData

#= If we want to read data directly via SQL Server
    conn = ODBC.Connection("CUPS", "crystal", "Bobject\$")
    #df = DBInterface.execute(conn, "SELECT TOP (100) * FROM [LIDS].[dbo].[FinanceForecastV1]") |> DataFrames.DataFrame
    #query = """ SELECT TOP (100) * FROM [LIDS].[dbo].[FinanceForecastV1];"""
    #df = DBInterface.execute(conn, query) |> DataFrames.DataFrame

    frt_sql = raw"C:\Users\SchroedR\OneDrive - Graphic Packaging International, Inc\Documents\MyDocuments\GPI Foodservice\SQL\Freight\ship_history_from_idh_r.sql"
    f = open("ship_history_from_idh_r.sql", "r")          
    f = open(frt_sql, "r")
    sql = read(f, String)
    ship_df = DBInterface.execute(conn, sql) |> DataFrames.DataFrame
=#

const DATA_DIR = joinpath(@__DIR__, "Data")

baseship_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_history.xlsx"), "Sheet 1")...
)

#VSCodeServer.vscodedisplay(ship_df)
DataFrames.describe(baseship_df)

# These do the same thing
baseship_df.Loc = string.(baseship_df.MSTCity, ", ", baseship_df.MSTState)
baseship_df.Loc = baseship_df.MSTCity .* ", " .* baseship_df.MSTState

#baseship_df.PTypeFix = Symbol.(baseship_df.PTypeFix)

# Get rid of the PlasticInenuity
baseship_df = baseship_df[baseship_df[:, :ShipSiteName] .!= "PLASTICS INGENUITY", :]

@pipe combine(baseship_df, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))
#dropmissing!(demand_df, :MSTZip)

#baseline_df = filter(:PrintGroup => ==("FIVEGUYS"), baseship_df)
dropmissing!(baseship_df, :MfgSiteName)

#get the BF volume.  Keep CPU here.
basebf_df = combine(groupby(baseship_df, 
    [:MfgSiteName, :ShipSiteName]), :MUnits => sum, renamecols = false)

bf_df = DataFrames.DataFrame(
        XLSX.readtable(joinpath(DATA_DIR, "bf_lanes.xlsx"), "Sheet1"; infer_eltypes=true)...
    )
basebf_df = leftjoin(basebf_df, bf_df, on = [:ShipSiteName => :ShipSite, :MfgSiteName => :Plant])

transform!(basebf_df, [:MUnits, :CostPerMUnit] => ByRow((a, b) -> (a * b)) => :Freight)
@pipe combine(basebf_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

# Copy to the clipboard and paste in the excel file
clipboard(sprint(show, "text/tab-separated-values", basebf_df))

baseprod_df = combine(groupby(baseship_df, 
    [:MfgSiteName]), :MUnits => sum, renamecols = false)
clipboard(sprint(show, "text/tab-separated-values", baseprod_df))


# Keep the location for the CPUs, but make the freight 0
#baseline_df = DataFrames.filter(:FrtTermsFix => !=("CPU"), baseline_df)
baseship_df.Loc = string.(baseship_df.MSTCity, ", ", baseship_df.MSTState)
baseship_df = DataFrames.combine(DataFrames.groupby(baseship_df, 
    [:ShipSiteName, :Loc, :FrtTermsFix]), :MUnits => sum, renamecols = false)
DataFrames.describe(baseship_df)

baseship_df[ismissing.(baseship_df.ShipSiteName) .| 
    ismissing.(baseship_df.Loc), :]
dropmissing!(baseship_df, disallowmissing = true)

baseship_df = leftjoin(baseship_df, dist_df, on = [:ShipSiteName => :Orig, :Loc => :Dest])
# TODO: show what is missing
dropmissing!(baseship_df, disallowmissing = true)

# Scalar f  freight in dollars per case/M per thousand miles  /90/ ;
f = 6.15 # freight in dollars per M per thousand miles 

g(a, b, c) = (c == "CPU" ? 0 : a * b * f/1000)
baseship_df.Freight = g.(baseship_df.MUnits, baseship_df.Miles, baseship_df.FrtTermsFix)

#transform!(baseline_df, [:MUnits, :Miles] => ByRow((a, b) -> (a * b * f / 1000)) => :Freight)

@pipe combine(baseship_df, :Freight => sum => :Freight, :MUnits => sum => :MUnits) |> 
    PrettyTables.pretty_table(_; formatters = PrettyTables.ft_printf("%'d"))

const OUT_DIR = joinpath(@__DIR__, "Data", "Baseline")
XLSX.writetable(joinpath(OUT_DIR, "baseship.xlsx"), overwrite = true, 
        collect(DataFrames.eachcol(baseship_df)), DataFrames.names(baseship_df))

# VSCodeServer.vscodedisplay(baseline_df)
#geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)

baseship_dmd = leftjoin(combine(groupby(baseship_df, :Loc), 
    :MUnits => sum => :MUnits), geocodes_fix, on = :Loc => :Dest)

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
    data=baseship_dmd,
    #data = combine(groupby(baseline_dmd, [:Dest, :lat, :lon]), :MUnits => sum => :MUnits),
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
        color=:ShipSiteName,
    data=baseship_df,
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
            lookup=:Loc,
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

