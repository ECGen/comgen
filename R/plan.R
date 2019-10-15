plan <- drake_plan(

    ## Load data
    ## Trees were removed from the analysis genotype RL6 and N1.31
    ## No physciods
    ## Lecanoras merged

    xgal.size <- read.csv("../data/lcn/ONC_Xgal_SizeData_May2011.csv"),
    garden.data <- read.csv("../data/lcn/LCO_data_ONC_PIT.csv"),
    rough <- read.csv("../data/lcn/ONC_raw_roughness.csv"),

    ## Data wrangling
                                        # rm genotype RL6 and N1.31
    garden.data <- garden.data[garden.data$Geno != "RL6", ],
    garden.data <- garden.data[garden.data$Tree != "N1.31", ],

    


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




