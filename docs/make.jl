using ForestMensuration
using Documenter

DocMeta.setdocmeta!(ForestMensuration, :DocTestSetup, :(using ForestMensuration); recursive=true)

makedocs(;
    modules=[ForestMensuration],
    authors="Marcos Daniel da Silva <marcosdasilva@5a.tec.br> and contributors",
    sitename="ForestMensuration.jl",
    format=Documenter.HTML(;
        canonical="https://marcosdanieldasilva.github.io/ForestMensuration.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/marcosdanieldasilva/ForestMensuration.jl",
    devbranch="main",
)
