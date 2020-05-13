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
    rm.geno = c("RL6", "T6", "999"),
    rm.tree = c("N1.31"),
    garden.data = proc_garden_data(garden.data.in, 
                                   rm.geno,
                                   rm.tree),
    ## Lichen Quadrats with Cell-wise Observations
    onc.q = proc_onc_q(garden.data, rm.zero = TRUE, rm.nless = 2),
    ## Tree Information
    onc.dat = proc_onc_dat(garden.data, rough.in, 
        onc.q, onc.nc.in, onc.tan.in, onc.ph.in),
    ## Size data 
    xgs.data = proc_size(xgal.size.in), 
### Modeling
    ## Lichen Network Models
    cn.onc = proc_cn_onc(onc.q),
    ## Lichen Network Model Similarity
    cn.d.onc = proc_cn_d_onc(cn.onc, onc.dat, 
        method = "euclidean", rm.na = TRUE),
    ## Lichen Community Matrix
    onc.com = proc_onc_com(garden.data, onc.q, onc.dat, rm.na = TRUE),
    ## Relativized Lichen Community Matrix
    onc.com.rel = proc_onc_com_rel(onc.com, onc.dat),
### Analyses
    ## REML for Genotype Effects
    ## These use tranformations, see R/functions/run_reml for details
    trait.results = run_trait_path(onc.dat),
    reml.results = run_reml(onc.dat, trait.results),
### Analytical Checks
    check.shapiro = check_shapiro(run_reml(onc.dat, trait.results, raw = TRUE)),
    check.fligner = check_fligner(onc.dat),
    ## PERMANOVAs for Network and Community Similarity
    ## NOTE: run_perm may transform the response distances, see R/functions
    perm.results = run_perm(onc.dat, onc.com.rel, cn.d.onc),
### Network Ordination
    cn.ord = run_nms(cn.d.onc, onc.dat[, c("BR", "L", "Cen", "AMI")]),
### Size analysis
    ## Size is square-rooted
    ## xg.reml = run_xgsize(xgs.data), 
    ## Test correlation between traits with genetic basis
    ## reg.trait = summary(lm(I(BR^(1/4))~CT, data = onc.dat)),
    ## Test correlation between networks and community comp
    ## cor.cn.com = vegan::mantel(cn.d.onc, vegdist(onc.com.rel)),
    ## Test correlation between traits and network metrics
    ## reg.trait.nm = run_trait_nm(onc.dat), 
    ## Get correlation values
    ## cor.trait.nm = cor(onc.dat[, -c(1, 6, 7)]),
### Species centrality analysis
    spp.cen = run_spp_centrality(cn.onc, onc.dat, cmode = "freeman", digits = 4),
    spp.cen.pos.in = run_spp_centrality(cn.onc, onc.dat, cmode = "in", type = "pos", digits = 4),
    spp.cen.pos.out = run_spp_centrality(cn.onc, onc.dat, cmode = "out", type = "pos", digits = 4),
    spp.cen.neg.in = run_spp_centrality(cn.onc, onc.dat, cmode = "in", type = "neg", digits = 4),
    spp.cen.neg.out = run_spp_centrality(cn.onc, onc.dat, cmode = "out", type = "neg", digits = 4),
### Species centrality table
    sppcen_xtab = make_table_sppcen(spp.cen.pos.in, spp.cen.pos.out, 
                                    spp.cen.neg.in, spp.cen.neg.out, digits = 4),
### correlation matrix for lichen and network
    cormat_xtab = cormat_tab(onc.dat),
### species area curves by genotype
    spac.g = spac_geno(onc.q, onc.dat),
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
    spp_cen_in.pdf = plot_sppcen(spp.cen.pos.in, file = "results/spp_cen_in.pdf"),
    spp_cen_out.pdf = plot_sppcen(spp.cen.pos.out, file = "results/spp_cen_out.pdf"),
    ## SUPPLEMENTARY
    xg_size.pdf = plot_xg_size(xgs.data, file = "results/xg_size.pdf"),
    spac_geno.pdf = plot_spag(spac.g, file = "results/spac_geno.pdf"),
### Tables
    ## tab:h2_table = Heritability table 
    ## tab:cn_perm_table = network similarity PERMANOVA
    ## tab:com_perm_table = community similarity PERMANOVA
    ## H2 table all
    xtab = make_tables(onc.dat, reml.results, perm.results, digits = 4),
    trait_path_xtab = make_table_path(trait.results, onc.dat),
    vec_tab = make_table_vectors(cn.ord[["vec"]]),
    vec_xtab = xtable(vec_tab, digits = 3),
    ## Update lichen manuscript tables and figures
    h2_reml.tex = print(
        xtab[["h2_reml"]], 
        file = "results/h2_reml.tex", 
        include.rownames = FALSE,
        include.colnames = TRUE),
    h2_reml_net.tex = print(
        xtab[["h2_net"]], 
        file = "results/h2_reml_net.tex", 
        include.rownames = FALSE,
        include.colnames = TRUE),
    h2_reml_trait.tex = print(
        xtab[["h2_trait"]], 
        file = "results/h2_reml_trait.tex", 
        include.rownames = FALSE,
        include.colnames = TRUE),
    trait_path.tex = print(
        trait_path_xtab, 
        file = "results/trait_path.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
    cn_perm.tex = print(
        xtab[["cn"]], 
        file = "results/cn_perm.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
    cn_trait_perm.tex = print(
        xtab[["cn_trait"]], 
        file = "results/cn_trait_perm.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
    com_perm.tex = print(
        xtab[["com"]], 
        file = "results/com_perm.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
    vec.tex = print(
        vec_xtab, 
        file = "results/vec.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
    sppcen.tex = print(
        sppcen_xtab, 
        file = "results/sppcen.tex", 
        include.rownames = FALSE,
        include.colnames = TRUE),
    cormat.tex = print(
        cormat_xtab, 
        file = "results/cormat.tex", 
        include.rownames = TRUE,
        include.colnames = TRUE),
### Tables and Figures for Manuscript
    tables_figures = list(
        h2_reml.tex = h2_reml.tex,
        h2_reml_net.tex = h2_reml_net.tex,
        h2_reml_trait.tex = h2_reml_trait.tex,
        trait_path.tex = trait_path.tex,
        cn_perm.tex = cn_perm.tex,
        cn_trait_perm.tex = cn_trait_perm.tex,
        com_perm.tex = com_perm.tex,
        vec.tex = vec.tex,
        cormat.tex = cormat.tex,
        cn_onc.pdf = cn_onc.pdf,
        h2_plot.pdf = h2_plot.pdf,
        spp_cen.pdf = spp_cen.pdf, 
        spp_cen_in.pdf = spp_cen_in.pdf,
        spp_cen_out.pdf = spp_cen_out.pdf,
        sppcen.tex = sppcen.tex,
        spac_geno.pdf = spac_geno.pdf,
        xg_size.pdf = xg_size.pdf
        ),
### Generate the manuscript
    update.manuscript = update_manuscript(
        files = tables_figures, 
        dir = "docs/lcn_manuscript", 
        file.tex = "main.tex",
        render = FALSE)
)
