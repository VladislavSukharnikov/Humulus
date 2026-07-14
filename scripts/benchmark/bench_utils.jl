using Humulus
using BenchmarkTools
using JLD2

# =============================================================================
# Utilities.
# =============================================================================

function print_header(title::AbstractString)
    println("\n")
    println("="^length(title))
    println(title)
    println("="^length(title))
    println()
end;

function report_evaluation(bench)
    allocated_bytes = median(bench).memory
    @info "Memory allocations:" allocated_bytes
end;

function benchmark_construction(object, bench)
    display(bench)

    allocated_bytes = median(bench).memory
    object_size_bytes = Base.summarysize(object)

    println()
    @info "Construction:" allocated_bytes object_size_bytes ratio = allocated_bytes / object_size_bytes
end

function benchmark_evaluation(bench)
    display(bench)

    allocated_bytes = median(bench).memory

    println()
    @info "Evaluation:" allocated_bytes
end