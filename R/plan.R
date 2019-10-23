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
    onc.ph = proc_onc_ph(garden.data, 
                         rough.in, 
                         onc, onc.q, 
                         onc.nc.in, onc.tan.in, onc.ph.in),
    onc.dat = proc_onc_dat(garden.data, rough.in, 
                           onc, onc.q, 
                           onc.nc.in, onc.tan.in, 
                           onc.ph),
    cn.onc = proc_cn_onc(onc),
    onc.ns = proc_onc_ns(cn.onc),
    cn.d.onc.na = proc_cn_d_onc(cn.onc, onc.dat, rm.na = TRUE),
    onc.com = proc_onc_com(garden.data, onc, onc.q),
    onc.com.rel = proc_onc_com_rel(onc.com),
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
    tabs = make_tables(onc.dat, cn.d.onc.na, onc.ns, onc.com.rel),
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
