using VegaLite
using VegaDatasets
using JSON3
using JSON
using URIs

const DATA_DIR = joinpath(@__DIR__, "Data")
# Read us states topojson data
test = JSON.parsefile(joinpath(DATA_DIR, "test.json"))
JSON3.pretty(test)
print(json(test, 4))
d=VegaDatasets.VegaJSONDataset(test, joinpath(DATA_DIR, "test.json"))
@vlplot(width=800, height=500) +
@vlplot(
    mark={
        :geoshape,
        fill=:lightgray,
        stroke=:green
    },
    data={
        values=d,
        format={
            type=:topojson, 
            feature=:"two-squares"
        }
    }
) + @vlplot(
    mark={
        :geoshape,
        stroke=:red
    },
    data={
        values=d,
        format={
            type=:topojson, 
            feature=:"one-line"
        }
    }
) + @vlplot(
    mark={
        :geoshape,
        stroke=:blue
    },
    data={
        values=d,
        format={
            type=:topojson, 
            feature=:"two-places"
        }
    }
)
