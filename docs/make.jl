using logic_cooling
using Documenter

DocMeta.setdocmeta!(logic_cooling, :DocTestSetup, :(using logic_cooling); recursive=true)

makedocs(;
    modules=[logic_cooling],
    authors="Rose-max111 <luyimingboy@163.com> and contributors",
    sitename="logic_cooling.jl",
    format=Documenter.HTML(;
        canonical="https://Rose_max111.github.io/logic_cooling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Rose_max111/logic_cooling.jl",
    devbranch="main",
)
