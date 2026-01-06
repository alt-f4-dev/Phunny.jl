using Phunny
using Documenter

DocMeta.setdocmeta!(Phunny, :DocTestSetup, :(using Phunny); recursive=true)

makedocs(;
    modules=[Phunny],
    authors="Isaac C. Ownby, Immanuel Schmidt",
    sitename="Phunny.jl",
    format=Documenter.HTML(;
    	prettyurls = get(ENV, "CI", "false") == "true",
        canonical="https://alt-f4-dev.github.io/Phunny.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => "tutorials.md",
        "Conventions"=> "conventions.md",
        "References" => "refs.md"
    ],
)

deploydocs(;
    repo="github.com/alt-f4-dev/Phunny.jl.git",
    devbranch="main",
)
