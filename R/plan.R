plan <- drake_plan(
    ## Load data
    ## Trees were removed from the analysis genotype RL6 and N1.31
    ## No physciods
    ## Lecanoras merged
    ## xgal.size.in = read.csv("./data/lcn/ONC_Xgal_SizeData_May2011.csv"),
    garden.data.in = read.csv("./data/lcn/LCO_data_ONC_PIT.csv"),
    rough.in = read.csv("./data/lcn/ONC_raw_roughness.csv"),
    onc.nc.in = read.csv("./data/lcn/ONC_phytochem_NC.csv"), 
    onc.tan.in = read.csv("./data/lcn/ONC_phytochem_tannin.csv"),
    onc.ph.in = read.csv("./data/lcn/ONC_Bark_lichen_pH_data.csv"),
    ## wild.dat.in = read.csv("./data/lcn/lco_Apr2012.csv"),
    ## env.in = read.csv("./data/lcn/Uinta2012_all_data_from_Lamit.csv"), 
    ## age.in = read.csv(
    ##     "./data/lcn/UintaMaster_LichenHeritNL_FallSpring_2012_ForLau.csv"),
    ## Data wrangling
    garden.data = proc_garden_data(garden.data.in),
    ## pit = proc_pit(garden.data.in),
    onc = proc_onc(garden.data),
    onc.q = proc_onc_q(onc),
    onc.dat = proc_onc_dat(garden.data, rough.in, 
                           onc, onc.q, 
                           onc.nc.in, onc.tan.in, onc.ph.in),
    cn.onc = proc_cn_onc(onc),
    onc.ns = proc_onc_ns(cn.onc),
    cn.d.onc = proc_cn_d_onc(cn.onc, onc.dat, rm.na = TRUE),
    onc.com = proc_onc_com(garden.data, onc, onc.q, onc.dat, rm.na = TRUE),
    onc.com.rel = proc_onc_com_rel(onc.com, onc.dat),
    ## Analyses
    reml.results = run_reml(onc.dat),
    perm.results = run_perm(onc.dat, onc.com.rel, cn.d.onc),
    ## Plots
                                        # fig: cn_onc
                                        # fig: cn:nms
                                        # fig: cen by geno
    ## Tables
                                        # table:onc_d_Ftable
                                        # H2 table all
    tables = make_tables(onc.dat, reml.results, perm.results, digits = 3),
    ## Update lichen manuscript tables and figures
    print(xtable(tables[["h2_reml"]], type = "latex"), 
          file = "results/h2_reml.tex", 
          include.rownames = FALSE,
          include.colnames = TRUE),
    print(xtable(tables[["cn"]], type = "latex"), 
          file = "results/cn_perm.tex", 
          include.rownames = TRUE,
          include.colnames = TRUE),
    print(xtable(tables[["com"]], type = "latex"), 
          file = "results/com_perm.tex", 
          include.rownames = TRUE,
          include.colnames = TRUE),
    ## Generate a report
    report.md = rmarkdown::render(
      knitr_in("R/report_lcn.Rmd"),
      output_dir = file_out("results/"),
      output_format = "md_document",
      quiet = TRUE
      ),
    report.html = rmarkdown::render(
      knitr_in("R/report_lcn.Rmd"),
      output_dir = file_out("results/"),
      output_format = "html_document",
      quiet = TRUE
      ),
    report.pdf = rmarkdown::render(
      knitr_in("R/report_lcn.Rmd"),
      output_dir = file_out("results/"),
      output_format = "pdf_document",
      quiet = TRUE
      )
)
