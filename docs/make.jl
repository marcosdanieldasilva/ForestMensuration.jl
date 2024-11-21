using ForestMensuration
using Documenter

DocMeta.setdocmeta!(ForestMensuration, :DocTestSetup, :(using ForestMensuration); recursive=true)

makedocs(;
    modules=[ForestMensuration],
    doctest=true,
    # linkcheck = true,
    authors="Marcos Daniel da Silva <marcosdasilva@5a.tec.br> and contributors",
    sitename="ForestMensuration.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://marcosdanieldasilva.github.io/ForestMensuration.jl",
        edit_link="main",
        assets=[
            joinpath("assets", "favicon.ico"),
            joinpath("assets", "style.css")
        ]
    ),
    checkdocs=:exports,
    pages=[
        "ForestMensuration Package" => "forestmensuration.md",
        "Getting Started" => "tutorial.md",
        "API Reference" => "reference.md",
        "Bibliography" => "bibliography.md",
        "Index" => "index.md"
    ]
)

deploydocs(;
    repo="github.com/marcosdanieldasilva/ForestMensuration.jl",
    devbranch="main",
    push_preview=true
)
