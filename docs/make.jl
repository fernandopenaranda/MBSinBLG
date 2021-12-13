using MBSinBLG
using Documenter

DocMeta.setdocmeta!(MBSinBLG, :DocTestSetup, :(using MBSinBLG); recursive=true)

makedocs(;
    modules=[MBSinBLG],
    authors="Fernando Peñaranda",
    repo="https://github.com/fernandopenaranda/MBSinBLG.jl/blob/{commit}{path}#{line}",
    sitename="MBSinBLG.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://fernandopenaranda.github.io/MBSinBLG.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/fernandopenaranda/MBSinBLG.jl",
    devbranch="main",
)
