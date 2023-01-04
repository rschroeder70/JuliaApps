# Model 4:  Allow trans shipments
# Assume CPU has the same freight penalty as PPD, and then zero out in post processing
#
# This from StackOverflow showing the big M and the binary variable:
# min sum((c,s,k), T(c,s,k)*x(c,s,k))
#     sum((c,k), x(c,s,k)) ≤ supply(s)    ∀s
#     sum(s, x(c,s,k)) = demand(c,k)      ∀c,k
#     x(c,s,k) ≤ M(c,s,k)*y(s)            ∀c,s,k     "if closed don't ship from DC" 
#     sum(s,y(s)) ≤ limit                            "limit number of open DCs" 
#     y(s) ∈ {0,1}                                   "DC is open or closed"
#     x(c,s,k) ≥ 0 
#     M(c,s,k) = min(supply(s),demand(c,k))          "constants" 

using JSON3

function pretty_print(d::Dict, pre=1)
    for (k,v) in d
        if typeof(v) <: Dict
            s = "$(repr(k)) => "
            println(join(fill(" ", pre)) * s)
            pretty_print(v, pre+1+length(s))
        else
            println(join(fill(" ", pre)) * "$(repr(k)) => $(repr(v))")
        end
    end
    nothing
end

data = JSON3.read("""
{
    "plants": {
        "P1": {"capacity": 200},
        "P2": {"capacity": 250},
        "P3": {"capacity": 300},
        "P4": {"capacity": 450}
    },
    "shipsites": {
        "S1": {"net": 0},
        "S2": {"net": 0},
        "S3": {"net": 0},
        "S4": {"net": 0}
    },
    "customers": {
        "C1": {"demand": 600},
        "C2": {"demand": 600}
    },
    "distances": {
        "P1 => S1": 0,
        "P1 => S2": 6,
        "P1 => S3": 24,
        "P1 => S4": 7,
        "S1 => C1": 24,
        "S1 => C2": 10,
        "P2 => S1": 6,
        "P2 => S2": 0,
        "P2 => S3": 10,
        "P2 => S4": 12,
        "S2 => C1": 5,
        "S2 => C2": 20,
        "P3 => S1": 24,
        "P3 => S2": 20,
        "P3 => S3": 0,
        "P3 => S4": 8,
        "S3 => C1": 45,
        "S3 => C2": 7,
        "P4 => S1": 7,
        "P4 => S2": 12,
        "P4 => S3": 8,
        "P4 => S4": 0,
        "S4 => C1": 30,
        "S4 => C2": 6
    }
}
""")

demand = JSON3.read("""
{
    "customers": {
        "C1": {"products" : {"T1": {"demand": 600},
                             "T2": {"demand": 500}
                            }
               },
        "C2": {"products" : {"T1": {"demand": 800},
                             "T2": {"demand": 200}
                             }
               }
        }
}
""")
print(json(demand, 4))
pretty_print(demand)

capacity = JSON3.read("""
{
    "plants": {
        "P1": {"products" : {"T1": {"capacity": 600},
                             "T2": {"capacity": 350}
                             }
              },
        "P2": {"products" : {"T1": {"capacity": 100},
                             "T2": {"capacity": 150}
                            }
                },
        "P3": {"products" : {"T1": {"capacity": 500},
                             "T2": {"capacity": 150}
                             }
                },
        "P4": {"products" : {"T1": {"capacity": 300},
                             "T2": {"capacity": 150}
                            }
                }
    }
}
""")
JSON3.pretty(JSON3.write(capacity))

C = collect(keys(demand["customers"]))
T = collect(keys(capacity["plants"]["P1"]["products"]))
P = collect(keys(capacity["plants"]))
S = collect(keys(data["shipsites"]))

distance(p::Symbol, m::Symbol) = data["distances"]["$(p) => $(m)"]

# Check for the Capacity by product
for t in T
    cap = 0
    for p in P
        cap = cap + values(capacity["plants"][p]["products"][t]["capacity"])
    end
    println("Capacity: $t = ", cap)
end

# Check for the Demand by product
for t in T
    dmd = 0
    for c in C
        dmd = dmd + values(demand["customers"][c]["products"][t]["demand"])
    end
    println("Demand: $t = ", dmd)
end

using JuMP
using GLPK
model = Model(GLPK.Optimizer)
#Plant to Ship site
@variable(model, y[P, S, T] >= 0)

#Ship Site to customer
@variable(model, x[S, C, T] >= 0)

# Binary variable if Site S ships to Cust C = 1, 0 otherwise
@variable(model, z[S, C], Bin)

# We need a constraint that each plant can ship no more than its capacity:
@constraint(model, [p in P, t in T], sum(y[p, :, t]) <= 
    capacity["plants"][p]["products"][t]["capacity"])

# customer demand
@constraint(model, [c in C, t in T], sum(x[:, c, t]) == 
    demand["customers"][c]["products"][t]["demand"])

# shipsite net = 0
@constraint(model, [s in S, t in T], sum(y[:, s, t]) == sum(x[s, :, t]))

#limit the ship sites to 1
@constraint(model, [c in C], sum(z[:, c]) == 1)

# Indicator constraint - not working
#@constraint(model, [s in S, c in C, t in T], !z[s, c] => {x[s, c, t] == 0})
@constraint(model, [s in S, c in C, t in T], x[s, c, t] <= 1000 * z[s, c])

@objective(model, Min, sum(distance(p, s) * y[p, s, t] for p in P, s in S, t in T) + 
    sum(distance(s, c) * x[s, c, t] for s in S, c in C, t in T))

optimize!(model)

solution_summary(model)

println("RESULTS:")
for p in P, s in S, t in T
    println("Plant: $p to Site: $s Product: $t => ", value(y[p, s, t]))
end

for s in S, c in C, t in T
    println("Site: $s to Customer: $c Product: $t => ", value(x[s, c, t]))
end

for s in S, c in C
    println("Binary Variable z[$s, $c] = ", value(z[s, c]))
end

value(x)