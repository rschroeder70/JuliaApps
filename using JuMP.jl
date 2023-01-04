using JuMP
import GLPK
import JSON

data = JSON.parse("""
{
    "plants": {
        "Seattle": {"capacity": 350},
        "San-Diego": {"capacity": 600}
    },
    "markets": {
        "New-York": {"demand": 300},
        "Chicago": {"demand": 300},
        "Topeka": {"demand": 300}
    },
    "distances": {
        "Seattle => New-York": 2.5,
        "Seattle => Chicago": 1.7,
        "Seattle => Topeka": 1.8,
        "San-Diego => New-York": 2.5,
        "San-Diego => Chicago": 1.8,
        "San-Diego => Topeka": 1.4
    }
}
""")

P = collect(keys(data["plants"]))


M = collect(keys(data["markets"]))

distance(p::String, m::String) = data["distances"]["$(p) => $(m)"]

model = Model(GLPK.Optimizer)

@variable(model, x[P, M] >= 0)

@constraint(model, [p in P], sum(x[p, :]) <= data["plants"][p]["capacity"])