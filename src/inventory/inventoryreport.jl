"""
    SamplingReport

Wraps the set of tables produced by a multi-table forest-inventory sampling design (for
example a per-stratum auxiliary table alongside the final result table) so they print and
travel together instead of as a bare `Vector`/`Tuple`. Designs that only ever produce one
table (e.g. [`simplecasualsampling`](@ref)) return a plain `DataFrame` instead —
wrapping a single table would add nothing.

# Fields
- `tables::NamedTuple`: the report's tables, keyed by name, e.g.
  `(anova=..., auxiliary_table=..., result_table=...)` for a stratified design, or
  `(occasion1=..., occasion2=..., change=...)` for a successive-occasions design.

# Examples

Individual tables are reachable both through `.tables` and directly as properties:

```julia-repl
julia> report = stratifiedsampling(:stratum, :volume, plot_area, total_area, data);
julia> report.result_table   # shorthand for report.tables.result_table
julia> report.auxiliary_table
```
"""
struct SamplingReport
  tables::NamedTuple
end

Base.propertynames(r::SamplingReport) = (:tables, propertynames(getfield(r, :tables))...)

function Base.getproperty(r::SamplingReport, name::Symbol)
  name === :tables && return getfield(r, :tables)
  return getfield(r, :tables)[name]
end

function Base.show(io::IO, r::SamplingReport)
  tables = getfield(r, :tables)
  names = propertynames(tables)
  for (i, name) in enumerate(names)
    println(io, titlecase(replace(string(name), "_" => " ")))
    show(io, tables[name])
    i < length(names) && (println(io); println(io))
  end
end
