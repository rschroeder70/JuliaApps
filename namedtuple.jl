using DataFrames
using NamedTupleTools
df = DataFrame(City=["Boston", "Dallas", "Chicago"], State=["MA", "TX", "IL"], Part = [:X, :Y, :Z])

#= 
Desired Output
(Loc = (City = "Boston", State = "MA"), Part = :X)
(Loc = (City = "Dallas", State = "TX"), Part = :Y)
(Loc = (City = "Chicago", State = "IL"), Part = :Z)
=#

struct Loc
    City::String
    State::String
end

# Loc as a named tuple
loc_nt = Tables.rowtable(DataFrames.select(df, :City, :State))
# convert to a structure
loc = structfromnt.(Loc, loc_nt)

struct Dmd
    Dest::Loc
    Part::Symbol
end

Demand.(Loc.(df.City, df.State), df.Part)

dmd_st = Dmd[]
for i in 1 : nrow(df)
    #println(cust_st[i])
    push!(dmd_st, Dmd(loc[i], df.Part[i]))
end

# Solution from Bogumit Kaminski from Julialang@discoursemail.com
Dmd.(Loc.(df.City, df.State), df.Part)

using DataFramesMeta
@with(df, Dmd.(Loc.(:City, :State), :Part))