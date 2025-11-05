using Pkg
Pkg.activate("docs")
Pkg.instantiate()

using Documenter, OpticalRayTracing

module MakeDocs

function setup(file)
    init = read(file, String)
    for page_file in pages
        page_file = "docs/src/" * page_file
        page = read(page_file, String)
        page = replace(page, r"```@(example sample|setup \b\w+\b).*?```"s =>
                             SubstitutionString("```@\\1\n$init```"))
        write(page_file, page)
    end
end

const branch = get(ENV, "GITHUB_REF_NAME", "")

const pages = [
    "index.md",
    "Ray Tracing.md",
    "Vignetting Analysis.md",
    "Seidel Aberrations.md",
]

end # module MakeDocs

MakeDocs.setup("docs/setup.jl")

DocMeta.setdocmeta!(
    OpticalRayTracing,
    :DocTestSetup,
    :(using OpticalRayTracing);
    recursive = true
)

makedocs(
    sitename = "OpticalRayTracing.jl",
    pages = [
        MakeDocs.pages;
        "Plotting Examples.md";
        "API" => "API.md"
    ],
    modules = [OpticalRayTracing],
    format = Documenter.HTML(warn_outdated = false)
)

if !isempty(MakeDocs.branch)
    deploydocs(
        repo = "github.com/Sagnac/OpticalRayTracing.jl.git",
        devbranch = MakeDocs.branch,
        devurl = MakeDocs.branch,
        versions = ["stable" => "master", "v^", "dev" => "dev"]
    )
end
