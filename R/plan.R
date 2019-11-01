plan <- drake_plan(
### Load data
    ## Lichen Quadrat Data from Ogden Nature Center
    garden.data.in = read.csv("./data/lcn/LCO_data_ONC_PIT.csv"),
    ## Roughness Data
    rough.in = read.csv("./data/lcn/ONC_raw_roughness.csv"),
    ## Chemistry Data
    onc.nc.in = read.csv("./data/lcn/ONC_phytochem_NC.csv"), 
    onc.tan.in = read.csv("./data/lcn/ONC_phytochem_tannin.csv"),
    onc.ph.in = read.csv("./data/lcn/ONC_Bark_lichen_pH_data.csv"),
    ## X. galericulata size estimates
    xgal.size.in = read.csv("./data/lcn/ONC_Xgal_SizeData_May2011.csv"),
### Data wrangling
    ## Genotypes and Trees to Remove
    rm.geno = c("RL6", "T6", "1007"),
    rm.tree = c("N1.31"),
    garden.data = proc_garden_data(garden.data.in, 
                                   rm.geno,
                                   rm.tree),
    ## Lichen Quadrats with Cell-wise Observations
    onc = proc_onc(garden.data),
    onc.q = proc_onc_q(onc),
    ## Tree Information
    onc.dat = proc_onc_dat(garden.data, rough.in, 
        onc, onc.q, 
        onc.nc.in, onc.tan.in, onc.ph.in),
    ## Size data 
    xgs.data = proc_size(xgal.size.in), 
### Modeling
    ## Lichen Network Models
    cn.onc = proc_cn_onc(onc),
    ## Lichen Network Metrics
    onc.ns = proc_onc_ns(cn.onc),
    ## Lichen Network Model Similarity
    cn.d.onc = proc_cn_d_onc(cn.onc, onc.dat, rm.na = TRUE),
    ## Lichen Community Matrix
    onc.com = proc_onc_com(garden.data, onc, onc.q, onc.dat, rm.na = TRUE),
    ## Relativized Lichen Community Matrix
    onc.com.rel = proc_onc_com_rel(onc.com, onc.dat),
### Analyses
    ## REML for Genotype Effects
    ## These use Fourth Root Transformations
    reml.reml = run_reml(onc.dat, raw = TRUE), 
    reml.results = run_reml(onc.dat),
    ## Analytical Checks
    check.shapiro = check_shapiro(reml.reml),
    check.fligner = check_fligner(onc.dat),
    ## PERMANOVAs for Network and Community Similarity
    perm.results = run_perm(onc.dat, onc.com.rel, cn.d.onc),
    ## Network Ordination
    cn.ord = run_nms(cn.d.onc, onc.dat[, c("CT", "BR", "PC", "SR")]),
    ## Size analysis
    ## Size is square-rooted
    xg.reml = run_xgsize(xgs.data), 
### Plots
    ## fig: cn_onc
    plot_nets(cn.onc, onc.dat, file = "./results/cn_onc.pdf"), 
    ## fig:cn_chplot
    plot_netsim(cn.ord, onc.dat, file = "./results/ch_plot.pdf"),
    ## fig: heritable traits
    plot_mdc(onc.dat, file = "./results/cn_metrics.pdf"),
    ## SUPPLEMENTARY
    plot_xg_size(xgs.data, file = "./results/xg_size.pdf"),
### Tables
    ## tab:h2_table = Heritability table 
    ## tab:cn_perm_table = network similarity PERMANOVA
    ## tab:com_perm_table = community similarity PERMANOVA
    ## H2 table all
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
