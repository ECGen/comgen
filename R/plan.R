plan <- drake_plan(
    ## Load data
    raw_data = lcn_lichen("data"),
    ## Data wrangling
    data = raw_data %>%
        mutate(Species = forcats::fct_inorder(Species)),
    ## Analyses
    fit = lm(Sepal.Width ~ Petal.Width + Species, data),
    ## Plots
    hist = create_plot(data),
    ## Generate a report
    report = rmarkdown::render(
         knitr_in("src/lcn_notebook.Rmd"),
         output_file = file_out("results/lcn_notebook.html"),
         quiet = TRUE
 )
)
