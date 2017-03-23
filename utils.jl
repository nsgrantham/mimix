### Utility functions for simulation study

import Iterators: product

function simulate{T <: Function}(generate::Function,
                  factors::Associative{Symbol, Any},
                  metrics::Associative{Symbol, String},
                  methods::Associative{Symbol, T};
                  reps::Int=1, dir=pwd(), parallel::Bool=true)
    mkdir(dir)
    empty_metrics = typeof(metrics)(zip(keys(metrics), fill("", length(metrics))))
    header = [:setting; :rep; :method; collect(keys(factors)); collect(keys(empty_metrics))]
    settings = dictproduct(factors)
    for (i, setting) in enumerate(settings)
        println("Simulating setting $i...")
        args = [((i - 1) * reps + j, setting, generate, methods) for j in 1:reps]
        results = parallel ? pmap(simulate_worker, args) : map(simulate_worker, args)
        open(joinpath(dir, "results-$i.csv"), "w") do f
            writecsvrow(f, header)
            for (j, result) in enumerate(results)
                for (method, output) in result
                    row = [i; j; method; collect(values(setting)); collect(values(leftmerge(empty_metrics, output)))]
                    writecsvrow(f, row)
                end
            end
        end
        println("Setting $i complete.")
    end
    println("Simulation study complete!")
end

function simulate_worker(args)
    seed, setting, generate, methods = args
    Y, X, Z, β = generate(; setting...)
    results = Dict{Symbol, Any}()
    shuffled_methods = shuffle(collect(keys(methods)))
    for method in shuffled_methods
        f = methods[method]
        srand(seed)
        results[method] = f(Y, X, Z, β)
    end
    results
end

function writecsvrow(f::IO, x)
    write(f, join(x, ","), "\n")
end

function dictproduct(dict::Associative)
    [typeof(dict)(zip(keys(dict), x)) for x in product(values(dict)...)]
end

function leftmerge(dict1::Associative, dict2::Associative)
    d2 = deepcopy(dict2)
    for key in keys(d2)
        if !(key in keys(dict1))
            delete!(d2, key)
        end
    end
    merge(dict1, d2)
end
