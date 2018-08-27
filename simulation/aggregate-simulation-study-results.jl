using ArgParse
using CSV
using DataFrames

function parse_commandline()
    s = ArgParseSettings("Aggregate simulation results from MCMC.")
    @add_arg_table s begin
        "directory"
            help = "Directory containing all simulation result tsv files."
    end
    return parse_args(s)
end

args = parse_commandline()

directory = abspath(args["directory"])
@assert isdir(directory)

global_test_dfs = DataFrame[]
local_estimates_dfs = DataFrame[] 
for (root, dirs, files) in walkdir(directory)
    for file in files
        if file in ["global-test.tsv", "local-estimates.tsv"]
            df = CSV.read(joinpath(root, file), delim='\t')
            colnames = names(df)
            path, rep = splitdir(root)
            path, setting = splitdir(path)
            path, model = splitdir(path)
            df[:model] = model
            df[:setting] = setting
            df[:rep] = rep
            df = df[vcat([:model, :setting, :rep], colnames)]
            if file == "global-test.tsv"
                push!(global_test_dfs, df)
            end
            if file == "local-estimates.tsv"
                push!(local_estimates_dfs, df)
            end
        end
    end
end

global_test_df = vcat(global_test_dfs...)
local_estimates_df = vcat(local_estimates_dfs...)

CSV.write(joinpath(directory, "global-test-results.tsv"), global_test_df, delim='\t')
CSV.write(joinpath(directory, "local-estimates-results.tsv"), local_estimates_df, delim='\t')