using ForestMensuration
using ForestModeling
using Documenter

DocMeta.setdocmeta!(ForestMensuration, :DocTestSetup, :(using ForestMensuration); recursive=true)

makedocs(;
    modules=[ForestMensuration, ForestModeling],
    doctest=true,
    # linkcheck = true,
    authors="Marcos Daniel da Silva <marcosdasilva@5a.tec.br> and contributors",
    sitename="ForestMensuration.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaForests.github.io/ForestMensuration.jl",
        edit_link="main",
        assets=[
            joinpath("assets", "favicon.ico"),
            joinpath("assets", "style.css")
        ],
        # reference.md now documents both ForestMensuration and the re-exported
        # ForestModeling API on one page, comfortably past Documenter's default
        # 200 KiB size warning/error threshold for a single generated page.
        size_threshold=400 * 1024,
        size_threshold_warn=250 * 1024,
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
    repo="github.com/JuliaForests/ForestMensuration.jl",
    devbranch="main",
    push_preview=true
)
