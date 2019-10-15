plan <- drake_plan(

    ## Load data
    
    ## Data wrangling

    ## Analyses
                                        # 1. network~geno+ permanova
    
                                        # 2. species centrality anova
    
                                        # 3. REMLs: cen, etc.
                                        # 4. Cen ~ SR, SE, CT
                                        # 5. Net ~ traits
    ## Plots
                                        # fig: cn_onc
                                        # fig: cn:nms
                                        # fig: cen by geno
    ## Tables
                                        # table:onc_d_Ftable
                                        # H2 table all
    ## Generate a report
    report = rmarkdown::render(
         knitr_in("R/report_lcn.Rmd"),
         output_file = file_out("results/report_lcn.html"),
         quiet = TRUE
 )
)




