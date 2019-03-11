using Weave
path = dirname(@__FILE__)
examples=[
    "benders_quatile_regression.jl"
]
cd(path)
for ex in examples
    weave(ex)
    weave(ex, doctype = "github")
    # weave("benders_quatile_regression.jl", doctype = "md2pdf")
end