using QuantumAdiabaticAnnealing
using Documenter

DocMeta.setdocmeta!(QuantumAdiabaticAnnealing, :DocTestSetup, :(using QuantumAdiabaticAnnealing); recursive=true)

makedocs(;
    modules=[QuantumAdiabaticAnnealing],
    authors="Rose-max111 <luyimingboy@163.com> and contributors",
    sitename="QuantumAdiabaticAnnealing.jl",
    format=Documenter.HTML(;
        canonical="https://Rose_max111.github.io/QuantumAdiabaticAnnealing.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Rose_max111/QuantumAdiabaticAnnealing.jl",
    devbranch="main",
)
