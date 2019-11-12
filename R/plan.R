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
    cn.ord = run_nms(cn.d.onc, onc.dat[, c(8, 15)]),
    ## Size analysis
    ## Size is square-rooted
    xg.reml = run_xgsize(xgs.data), 
    ## Test correlation between traits with genetic basis
    reg.trait = summary(lm(I(BR^(1/4))~CT, data = onc.dat)),
    ## Test correlation between networks and community comp
    cor.cn.com = vegan::mantel(cn.d.onc, vegdist(onc.com.rel)),
    ## Test correlation between traits and network metrics
    reg.trait.nm = run_trait_nm(onc.dat), 
    ## Get correlation values
    cor.trait.nm = cor(onc.dat[, -c(1, 6, 7)]),
    ## Species centrality analysis
    spp.cen = run_spp_centrality(cn.onc, onc.dat),
### Plots
    ## fig:cn_onc
    cn_onc.pdf = plot_nets(cn.onc, onc.dat, file = "results/cn_onc.pdf"),
    ## fig:h2_plot
    h2_plot.pdf = plot_h2(cn.ord, 
        onc.dat, sig.alpha = 0.15, 
        plot.vectors = TRUE,
        file = "results/h2_plot.pdf"),
    ## fig:spp_cen
    spp_cen.pdf = plot_sppcen(spp.cen, file = "results/spp_cen.pdf"),
    ## SUPPLEMENTARY
    xg_size.pdf = plot_xg_size(xgs.data, file = "results/xg_size.pdf"),
### Tables
    ## tab:h2_table = Heritability table 
    ## tab:cn_perm_table = network similarity PERMANOVA
    ## tab:com_perm_table = community similarity PERMANOVA
    ## H2 table all
    xtab = make_tables(onc.dat, reml.results, perm.results, digits = 4),
    ## Update lichen manuscript tables and figures
    h2_reml.tex = print(
        xtab[["h2_reml"]], 
        file = "results/h2_reml.tex", 
        include.rownames = FALSE,
        include.colnames = TRUE),
    cn_perm.tex = print(
        xtab[["cn"]], 
        file = "results/cn_perm.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
    com_perm.tex = print(
        xtab[["com"]], 
        file = "results/com_perm.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
    ## Tables and Figures for Manuscript
    tables_figures = list(
        h2_reml.tex = h2_reml.tex,
        cn_perm.tex = cn_perm.tex,
        com_perm.tex = com_perm.tex,
        cn_onc.pdf = cn_onc.pdf,
        h2_plot.pdf = h2_plot.pdf,
        spp_cen.pdf = spp_cen.pdf, 
        xg_size.pdf = xg_size.pdf
        ),
    ## Generate the manuscript
    update.manuscript = update_manuscript(
        files = tables_figures, 
        dir = "docs/lcn_manuscript", 
        file.tex = "main.tex")
)
