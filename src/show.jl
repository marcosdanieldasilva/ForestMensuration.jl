import Base: show

# Define the show method for custom display of ModelEquation
function show(io::IO, model::ModelEquation)
  print(io, model.output)
end

# Custom show method for SiteAnalysis to display  the site_table and site_plot
function show(io::IO, analysis::SiteAnalysis)
  show(io, analysis.site_table)
  display(analysis.site_plot)
end