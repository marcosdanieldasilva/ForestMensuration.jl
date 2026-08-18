"""
    SamplingReport

Wraps the set of tables produced by a multi-table forest-inventory sampling design (for
example a per-stratum auxiliary table alongside the final result table) so they print and
travel together instead of as a bare `Vector`/`Tuple`. Designs that only ever produce one
table (e.g. [`simplecasualsampling`](@ref)) return a plain `DataFrame` instead —
wrapping a single table would add nothing.

# Fields
- `tables::NamedTuple`: the report's tables, keyed by name, e.g.
  `(anova=..., auxiliaryTable=..., resultTable=...)` for a stratified design, or
  `(occasion1=..., occasion2=..., change=...)` for a successive-occasions design.

# Examples

Individual tables are reachable both through `.tables` and directly as properties:

```julia-repl
julia> report = stratifiedsampling(:stratum, :volume, plot_area, total_area, data);
julia> report.resultTable   # shorthand for report.tables.resultTable
julia> report.auxiliaryTable
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
    println(io, titlecase(replace(string(name), r"(?<=[a-z])(?=[A-Z])" => " ")))
    show(io, tables[name])
    i < length(names) && (println(io); println(io))
  end
end

"""
    resultTable(report::SamplingReport)

The final one-row, one-column-per-statistic result table — same as `report.resultTable`.
Produced by every multi-table design except the four successive-occasions designs
(see [`occasion1`](@ref)/[`occasion2`](@ref)/[`change`](@ref) instead). The function-call
equivalent of dot access, matching how `AllometricModel`'s own results are reached
through functions (`coef`, `confint`, ...) rather than raw field access elsewhere in
this ecosystem; calling it on a report that has no such table errors the same way
`report.resultTable` already would.

# Examples
```julia-repl
julia> report = clustersampling(:cluster, :volume, plot_area, total_area, data);

julia> resultTable(report)   # same as report.resultTable
```
"""
resultTable(r::SamplingReport) = r.resultTable

"""
    clusterTable(report::SamplingReport)

Per-cluster descriptive statistics (`n`, mean, variance) from [`clustersampling`](@ref) —
same as `report.clusterTable`.
"""
clusterTable(r::SamplingReport) = r.clusterTable

"""
    auxiliaryTable(report::SamplingReport)

Per-stratum `n`, mean, variance, and allocation weights from
[`stratifiedsampling`](@ref) — same as `report.auxiliaryTable`.
"""
auxiliaryTable(r::SamplingReport) = r.auxiliaryTable

"""
    pointTable(report::SamplingReport)

Per-point tree count, basal area/ha, and volume/ha from
[`horizontalpointsampling`](@ref) (only points with at least one counted tree) — same
as `report.pointTable`.
"""
pointTable(r::SamplingReport) = r.pointTable

"""
    startTable(report::SamplingReport)

Per-start descriptive statistics from [`multistartsystematicsampling`](@ref) — same as
`report.startTable`.
"""
startTable(r::SamplingReport) = r.startTable

"""
    primaryTable(report::SamplingReport)

Per-primary-unit descriptive statistics from [`twostagesampling`](@ref) — same as
`report.primaryTable`.
"""
primaryTable(r::SamplingReport) = r.primaryTable

"""
    anova(report::SamplingReport)

The between/within-strata analysis-of-variance table from [`stratifiedsampling`](@ref) —
same as `report.anova`.
"""
anova(r::SamplingReport) = r.anova

"""
    occasion1(report::SamplingReport)

The first-occasion simple-random-sampling result table, from any of the four
sampling-on-successive-occasions designs — same as `report.occasion1`.
"""
occasion1(r::SamplingReport) = r.occasion1

"""
    occasion2(report::SamplingReport)

The second-occasion result table, from any of the four sampling-on-successive-occasions
designs — same as `report.occasion2`.
"""
occasion2(r::SamplingReport) = r.occasion2

"""
    change(report::SamplingReport)

The estimated growth (mean and total change) between occasions, from any of the four
sampling-on-successive-occasions designs — same as `report.change`.
"""
change(r::SamplingReport) = r.change
