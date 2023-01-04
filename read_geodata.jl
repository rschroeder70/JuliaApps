module LoadGeoData

export dist_df, distw_df, geocodes_df, plants_df, geocodes_fix

using XLSX
using DataFrames
using Statistics

#=
    Read all the geocode stuff
=#

const DATA_DIR = joinpath(@__DIR__, "Data")

dist_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_to_dist.xlsx"), "Sheet 1"; infer_eltypes=true)
)

DataFrames.describe(dist_df)
dist_df.Miles = passmissing(convert).(Float64, dist_df.Miles)
# An alternative Here
#dist_df[!, :Miles] = convert.(Float64, dist_df[!, :Miles])

dist_df[ismissing.(dist_df."Distance.Distance") .| ismissing.(dist_df.Miles), :]

# remove EX rows and PR
filter!(row -> row."de.State" != "EX", dist_df)
filter!(row -> row."de.State" != "PR", dist_df)
filter!(row -> row."de.State" != "HI", dist_df)
filter!(row -> row."de.State" != "XX", dist_df)
filter!(row -> !((row."de.City" == "GNANGARA") & (row."de.State" == "WA")), dist_df)

dist_df[ismissing.(dist_df."Distance.Distance") .| ismissing.(dist_df.Miles), :]

# remove or.State
dist_df = dist_df[:, InvertedIndices.Not("or.State")] 
dist_df = dist_df[:, InvertedIndices.Not("Distance.Distance")]

dist_df.Dest = string.(dist_df."de.City", ", ", dist_df."de.State")
DataFrames.rename!(dist_df, "or.City" => "Orig")

DataFrames.select!(dist_df, :Orig, :Dest, :Miles)

# Change to wide format
distw_df = DataFrames.unstack(dist_df, :Orig, :Miles)

# Fix datatypes in the dist_df
dist_df[!, :Miles] = convert.(Float64, dist_df[:, :Miles])

# Read the geocodes
geocodes_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "ship_to_geocodes.xlsx"), "Sheet 1"; infer_eltypes=true)
)

geocodes_df.Dest = string.(geocodes_df.MSTCity, ", ", geocodes_df.MSTState)

DataFrames.describe(geocodes_df)

geocodes_df[ismissing.(geocodes_df.lat) .| ismissing.(geocodes_df.lon), :]
geocodes_df.lat = convert.(Float64, geocodes_df.lat)
geocodes_df.lon = convert.(Float64, geocodes_df.lon)
geocodes_fix = combine(groupby(geocodes_df, [:Dest]), :lat => mean => :lat, :lon => mean => :lon)


# Read the plant locations
plants_df = DataFrames.DataFrame(
    XLSX.readtable(joinpath(DATA_DIR, "Facilities.xlsx"), "Facilities"; infer_eltypes=true)
)

DataFrames.select!(plants_df, :Name, :Latitude, :Longitude)

end 