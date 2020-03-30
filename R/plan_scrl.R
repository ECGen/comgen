plan_scrl <- drake_plan(
### Load data
### Wrangling
### Modeling
### Analyses
### Analytical Checks
### Figures
    ## xg_size.pdf = plot_xg_size(xgs.data, file = "results/xg_size.pdf"),
### Tables
### Tables and Figures for Manuscript
    ## tables_figures = list(
    ##     ),
### Supplementary
### Generate the manuscript
    ## update.manuscript = update_manuscript(
    ##     files = tables_figures, 
    ##     dir = "docs/lcn_manuscript", 
    ##    file.tex = "main.tex")
)
