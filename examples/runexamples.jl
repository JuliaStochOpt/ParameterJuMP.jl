using Weave
path = dirname(@__FILE__)
examples=[
    "benders_quatile_regression.jl"
]
cd(path)
for ex in examples
    weave(ex)
    # weave(ex, doctype = "github")
    # weave(ex, doctype = "notebook")
    # Weave.convert_doc(ex, ex[1:end-3]*".ipynb")
    # weave(ex, doctype = "md2pdf")
end