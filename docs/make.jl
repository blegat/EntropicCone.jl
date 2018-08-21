using Documenter, EntropicCone

makedocs(
    format = :html,
    sitename = "EntropicCone",
    pages = [
        "Index" => "index.md",
        "Introduction" => "intro.md",
        "Entropic Vector" => "vector.md",
    ],
    # The following ensures that we only include the docstrings from
    # this module for functions define in Base that we overwrite.
    modules = [EntropicCone]
)

deploydocs(
    repo   = "github.com/blegat/EntropicCone.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing
)
