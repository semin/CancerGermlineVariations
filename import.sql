DROP TABLE IF EXISTS snp_annotations;
CREATE TABLE snp_annotations AS 
  SELECT
    m.chrom, m.pos, m.id, m.qual, m.ref, m.alt, m.alt1, m.alt2, m.alt3, m.rsid1, m.rsid2, m.rsid3, m.all_ac1, m.all_ac2, m.all_ac3, m.all_af1, m.all_af2, m.all_af3, m.eur_ac1, m.eur_ac2, m.eur_ac3, m.eur_af1, m.eur_af2, m.eur_af3, m.tg_all_af1, m.tg_all_af2, m.tg_all_af3, m.tg_asn_af1, m.tg_asn_af2, m.tg_asn_af3, m.tg_amr_af1, m.tg_amr_af2, m.tg_amr_af3, m.tg_afr_af1, m.tg_afr_af2, m.tg_afr_af3, m.tg_eur_af1, m.tg_eur_af2, m.tg_eur_af3, m.chisq_pval_all_ac1,  m.chisq_pval_all_ac2,  m.chisq_pval_all_ac3, m.fisher_pval_all_ac1, m.fisher_pval_all_ac2, m.fisher_pval_all_ac3, m.chisq_pval_eur_ac1,  m.chisq_pval_eur_ac2,  m.chisq_pval_eur_ac3, m.fisher_pval_eur_ac1, m.fisher_pval_eur_ac2, m.fisher_pval_eur_ac3, m.cadd_raw1, m.cadd_raw2, m.cadd_raw3, m.cadd_phred1, m.cadd_phred2, m.cadd_phred3,
    v.region AS region, v.feature AS feature,
    e.exonic_vartype AS exonic_vartype, e.exonic_change AS exonic_change,
    s1.rsid AS dbsnp137, s2.rsid AS dbsnp138,
    t1.af AS av_tg_all_af, t2.af AS av_tg_eur_af, t3.af AS av_tg_afr_af, t4.af AS av_tg_amr_af, t5.af AS av_tg_asn_af,
    e1.af AS esp6500_all_af, e2.af AS esp6500_ea_af, e3.af AS esp6500_aa_af,
    cs.id AS cosmic,
    sf.score AS sift,
    pd.score AS pp2hdiv,
    pv.score AS pp2hvar,
    ma.score AS ma,
    mt.score AS mt,
    g.phenotype AS gwas_catalog,
    c1.score AS phast_cons_elements46way,
    c2.score AS phast_cons_elements46way_placental,
    c3.score AS phast_cons_elements46way_primates,
    dc.num_of_celllines AS wg_encode_reg_dnase_clustered_v2,
    t.tfs AS tfbs_cons_sites,
    tc.tfs AS wg_encode_reg_tfbs_clustered_v3,
    h1.state AS wg_encode_broad_hmm_gm12878,
    h2.state AS wg_encode_broad_hmm_h1hesc,
    h3.state AS wg_encode_broad_hmm_hepg2,
    h4.state AS wg_encode_broad_hmm_hmec,
    h5.state AS wg_encode_broad_hmm_hsmm,
    h6.state AS wg_encode_broad_hmm_huvec,
    h7.state AS wg_encode_broad_hmm_k562,
    h8.state AS wg_encode_broad_hmm_nhek,
    h9.state AS wg_encode_broad_hmm_nhlf,
    r.type AS repeat_masker,
    s.type AS simple_repeat,
    n.type AS nested_repeat,
    wr.rna AS wg_rna,
    ts.target AS mirna_target,
    d.var AS dgv,
    ms.type AS microsatelite,
    d1.num_of_celllines AS combined_dhs_peak,
    d2.genes AS distal_dhs_to_promoter_dhs,
    d3.genes AS dhs_to_gene_expression,
    p1.num_of_celllines AS promoter_dhs_master_known,
    p2.num_of_celllines AS promoter_dhs_master_novel,
    quote(m.TCGA_A6_2670_10A) AS TCGA_A6_2670_10A, quote(m.TCGA_A6_2671_10A) AS TCGA_A6_2671_10A, quote(m.TCGA_A6_2674_10A) AS TCGA_A6_2674_10A, quote(m.TCGA_A6_2676_10A) AS TCGA_A6_2676_10A, quote(m.TCGA_A6_2678_10A) AS TCGA_A6_2678_10A, quote(m.TCGA_A6_2679_10A) AS TCGA_A6_2679_10A, quote(m.TCGA_A6_2680_10A) AS TCGA_A6_2680_10A, quote(m.TCGA_A6_2681_10A) AS TCGA_A6_2681_10A, quote(m.TCGA_A6_2682_10A) AS TCGA_A6_2682_10A, quote(m.TCGA_A6_2683_10A) AS TCGA_A6_2683_10A, quote(m.TCGA_A6_2684_10A) AS TCGA_A6_2684_10A, quote(m.TCGA_A6_3807_11A) AS TCGA_A6_3807_11A, quote(m.TCGA_A6_3808_11A) AS TCGA_A6_3808_11A, quote(m.TCGA_A6_3810_11A) AS TCGA_A6_3810_11A, quote(m.TCGA_AA_3494_11A) AS TCGA_AA_3494_11A, quote(m.TCGA_AA_3495_11A) AS TCGA_AA_3495_11A, quote(m.TCGA_AA_3502_11A) AS TCGA_AA_3502_11A, quote(m.TCGA_AA_3509_11A) AS TCGA_AA_3509_11A, quote(m.TCGA_AA_3510_11A) AS TCGA_AA_3510_11A, quote(m.TCGA_AA_3514_11A) AS TCGA_AA_3514_11A, quote(m.TCGA_AA_3516_10A) AS TCGA_AA_3516_10A, quote(m.TCGA_AA_3529_11A) AS TCGA_AA_3529_11A, quote(m.TCGA_AA_3548_10A) AS TCGA_AA_3548_10A, quote(m.TCGA_AA_3549_10A) AS TCGA_AA_3549_10A, quote(m.TCGA_AA_3553_10A) AS TCGA_AA_3553_10A, quote(m.TCGA_AA_3555_10A) AS TCGA_AA_3555_10A, quote(m.TCGA_AA_3558_10A) AS TCGA_AA_3558_10A, quote(m.TCGA_AA_3664_10A) AS TCGA_AA_3664_10A, quote(m.TCGA_AA_3666_10A) AS TCGA_AA_3666_10A, quote(m.TCGA_AA_3672_10A) AS TCGA_AA_3672_10A, quote(m.TCGA_AA_3675_10A) AS TCGA_AA_3675_10A, quote(m.TCGA_AA_3681_10A) AS TCGA_AA_3681_10A, quote(m.TCGA_AA_3685_10A) AS TCGA_AA_3685_10A, quote(m.TCGA_AA_3688_10A) AS TCGA_AA_3688_10A, quote(m.TCGA_AA_3692_10A) AS TCGA_AA_3692_10A, quote(m.TCGA_AA_3693_10A) AS TCGA_AA_3693_10A, quote(m.TCGA_AA_3715_10A) AS TCGA_AA_3715_10A, quote(m.TCGA_AA_3812_10A) AS TCGA_AA_3812_10A, quote(m.TCGA_AA_3861_10A) AS TCGA_AA_3861_10A, quote(m.TCGA_AA_3956_10A) AS TCGA_AA_3956_10A, quote(m.TCGA_AA_3966_10A) AS TCGA_AA_3966_10A, quote(m.TCGA_AA_3968_10A) AS TCGA_AA_3968_10A, quote(m.TCGA_AA_3970_10A) AS TCGA_AA_3970_10A, quote(m.TCGA_AA_3994_10A) AS TCGA_AA_3994_10A, quote(m.TCGA_AA_A00Z_10A) AS TCGA_AA_A00Z_10A, quote(m.TCGA_AA_A01C_10A) AS TCGA_AA_A01C_10A, quote(m.TCGA_AA_A01D_10A) AS TCGA_AA_A01D_10A, quote(m.TCGA_AA_A01G_10A) AS TCGA_AA_A01G_10A, quote(m.TCGA_AA_A01I_10A) AS TCGA_AA_A01I_10A, quote(m.TCGA_AA_A01K_10A) AS TCGA_AA_A01K_10A, quote(m.TCGA_AA_A01P_11A) AS TCGA_AA_A01P_11A, quote(m.TCGA_AA_A01Q_10A) AS TCGA_AA_A01Q_10A, quote(m.TCGA_AA_A01R_11A) AS TCGA_AA_A01R_11A, quote(m.TCGA_AA_A01S_11A) AS TCGA_AA_A01S_11A, quote(m.TCGA_AA_A01T_11A) AS TCGA_AA_A01T_11A, quote(m.TCGA_AA_A01V_11A) AS TCGA_AA_A01V_11A, quote(m.TCGA_AA_A01X_11A) AS TCGA_AA_A01X_11A, quote(m.TCGA_AA_A02K_10A) AS TCGA_AA_A02K_10A, quote(m.TCGA_AA_A02O_11A) AS TCGA_AA_A02O_11A, quote(m.TCGA_AA_A02R_10A) AS TCGA_AA_A02R_10A, quote(m.TCGA_AA_A02W_10A) AS TCGA_AA_A02W_10A, quote(m.TCGA_AA_A02Y_10A) AS TCGA_AA_A02Y_10A, quote(m.TCGA_AA_A03F_11A) AS TCGA_AA_A03F_11A, quote(m.TCGA_AA_A03J_11A) AS TCGA_AA_A03J_11A, quote(m.TCGA_AF_2689_10A) AS TCGA_AF_2689_10A, quote(m.TCGA_AF_2691_10A) AS TCGA_AF_2691_10A, quote(m.TCGA_AF_2692_10A) AS TCGA_AF_2692_10A, quote(m.TCGA_AF_3913_11A) AS TCGA_AF_3913_11A, quote(m.TCGA_AG_3574_10A) AS TCGA_AG_3574_10A, quote(m.TCGA_AG_3582_10A) AS TCGA_AG_3582_10A, quote(m.TCGA_AG_3661_10A) AS TCGA_AG_3661_10A, quote(m.TCGA_AG_3664_10A) AS TCGA_AG_3664_10A, quote(m.TCGA_AG_3727_10A) AS TCGA_AG_3727_10A, quote(m.TCGA_AG_3728_10A) AS TCGA_AG_3728_10A, quote(m.TCGA_AG_3878_10A) AS TCGA_AG_3878_10A, quote(m.TCGA_AG_3881_10A) AS TCGA_AG_3881_10A, quote(m.TCGA_AG_3882_10A) AS TCGA_AG_3882_10A, quote(m.TCGA_AG_3885_10A) AS TCGA_AG_3885_10A, quote(m.TCGA_AG_3887_10A) AS TCGA_AG_3887_10A, quote(m.TCGA_AG_3890_10A) AS TCGA_AG_3890_10A, quote(m.TCGA_AG_3892_10A) AS TCGA_AG_3892_10A, quote(m.TCGA_AG_3893_10A) AS TCGA_AG_3893_10A, quote(m.TCGA_AG_3894_10A) AS TCGA_AG_3894_10A, quote(m.TCGA_AG_3896_10A) AS TCGA_AG_3896_10A, quote(m.TCGA_AG_3898_10A) AS TCGA_AG_3898_10A, quote(m.TCGA_AG_3901_10A) AS TCGA_AG_3901_10A, quote(m.TCGA_AG_3902_10A) AS TCGA_AG_3902_10A, quote(m.TCGA_AG_3909_10A) AS TCGA_AG_3909_10A, quote(m.TCGA_AG_3999_10A) AS TCGA_AG_3999_10A, quote(m.TCGA_AG_4001_10A) AS TCGA_AG_4001_10A, quote(m.TCGA_AG_4005_10A) AS TCGA_AG_4005_10A, quote(m.TCGA_AG_4007_10A) AS TCGA_AG_4007_10A, quote(m.TCGA_AG_4008_10A) AS TCGA_AG_4008_10A, quote(m.TCGA_AG_4015_10A) AS TCGA_AG_4015_10A, quote(m.TCGA_AG_A00C_10A) AS TCGA_AG_A00C_10A, quote(m.TCGA_AG_A00H_10A) AS TCGA_AG_A00H_10A, quote(m.TCGA_AG_A00Y_10A) AS TCGA_AG_A00Y_10A, quote(m.TCGA_AG_A011_10A) AS TCGA_AG_A011_10A, quote(m.TCGA_AG_A014_10A) AS TCGA_AG_A014_10A, quote(m.TCGA_AG_A015_10A) AS TCGA_AG_A015_10A, quote(m.TCGA_AG_A01J_10A) AS TCGA_AG_A01J_10A, quote(m.TCGA_AG_A01L_10A) AS TCGA_AG_A01L_10A, quote(m.TCGA_AG_A01N_10A) AS TCGA_AG_A01N_10A, quote(m.TCGA_AG_A01W_11A) AS TCGA_AG_A01W_11A, quote(m.TCGA_AG_A01Y_11A) AS TCGA_AG_A01Y_11A, quote(m.TCGA_AG_A020_11A) AS TCGA_AG_A020_11A, quote(m.TCGA_AG_A032_10A) AS TCGA_AG_A032_10A, quote(m.TCGA_AY_4070_10A) AS TCGA_AY_4070_10A, quote(m.TCGA_AY_4071_10A) AS TCGA_AY_4071_10A, quote(m.TCGA_AZ_4308_10A) AS TCGA_AZ_4308_10A, quote(m.TCGA_AZ_4313_10A) AS TCGA_AZ_4313_10A, quote(m.TCGA_AZ_4315_10A) AS TCGA_AZ_4315_10A, quote(m.TCGA_AZ_4614_10A) AS TCGA_AZ_4614_10A, quote(m.TCGA_AZ_4615_10A) AS TCGA_AZ_4615_10A, quote(m.TCGA_AZ_4681_10A) AS TCGA_AZ_4681_10A, quote(m.TCGA_AZ_4682_10A) AS TCGA_AZ_4682_10A, quote(m.TCGA_AZ_4684_10A) AS TCGA_AZ_4684_10A, quote(m.TCGA_CA_5256_10A) AS TCGA_CA_5256_10A, quote(m.TCGA_CK_4951_10A) AS TCGA_CK_4951_10A, quote(m.TCGA_CM_4744_10A) AS TCGA_CM_4744_10A, quote(m.TCGA_CM_4746_10A) AS TCGA_CM_4746_10A, quote(m.TCGA_CM_4747_10A) AS TCGA_CM_4747_10A, quote(m.TCGA_CM_4748_10A) AS TCGA_CM_4748_10A, quote(m.TCGA_CM_4750_10A) AS TCGA_CM_4750_10A, quote(m.TCGA_CM_4752_10A) AS TCGA_CM_4752_10A
    FROM snps AS m
    LEFT JOIN snp_variant_functions AS v ON m.chrom = v.vchrom AND m.pos = v.pos AND m.ref = v.ref AND m.alt1 = v.alt
    LEFT JOIN snp_exonic_variant_functions AS e ON m.chrom = e.vchrom AND m.pos = e.pos AND m.ref = e.ref AND m.alt1 = e.alt
    LEFT JOIN snp_dbsnp137 AS s1 ON m.chrom = s1.vchrom AND m.pos = s1.pos AND m.ref = s1.ref AND m.alt1 = s1.alt
    LEFT JOIN snp_dbsnp138 AS s2 ON m.chrom = s2.vchrom AND m.pos = s2.pos AND m.ref = s2.ref AND m.alt1 = s2.alt
    LEFT JOIN snp_tg_all_afs AS t1 ON m.chrom = t1.vchrom AND m.pos = t1.pos AND m.ref = t1.ref AND m.alt1 = t1.alt
    LEFT JOIN snp_tg_eur_afs AS t2 ON m.chrom = t2.vchrom AND m.pos = t2.pos AND m.ref = t2.ref AND m.alt1 = t2.alt
    LEFT JOIN snp_tg_afr_afs AS t3 ON m.chrom = t3.vchrom AND m.pos = t3.pos AND m.ref = t3.ref AND m.alt1 = t3.alt
    LEFT JOIN snp_tg_amr_afs AS t4 ON m.chrom = t4.vchrom AND m.pos = t4.pos AND m.ref = t4.ref AND m.alt1 = t4.alt
    LEFT JOIN snp_tg_asn_afs AS t5 ON m.chrom = t5.vchrom AND m.pos = t5.pos AND m.ref = t5.ref AND m.alt1 = t5.alt
    LEFT JOIN snp_esp_all_afs AS e1 ON m.chrom = e1.vchrom AND m.pos = e1.pos AND m.ref = e1.ref AND m.alt1 = e1.alt
    LEFT JOIN snp_esp_ea_afs AS e2 ON m.chrom = e2.vchrom AND m.pos = e2.pos AND m.ref = e2.ref AND m.alt1 = e2.alt
    LEFT JOIN snp_esp_aa_afs AS e3 ON m.chrom = e3.vchrom AND m.pos = e3.pos AND m.ref = e3.ref AND m.alt1 = e3.alt
    LEFT JOIN snp_cosmic68 AS cs ON m.chrom = cs.vchrom AND m.pos = cs.pos AND m.ref = cs.ref AND m.alt1 = cs.alt
    LEFT JOIN snp_ljb23_sift AS sf ON m.chrom = sf.vchrom AND m.pos = sf.pos AND m.ref = sf.ref AND m.alt1 = sf.alt
    LEFT JOIN snp_ljb23_pp2hdiv AS pd ON m.chrom = pd.vchrom AND m.pos = pd.pos AND m.ref = pd.ref AND m.alt1 = pd.alt
    LEFT JOIN snp_ljb23_pp2hvar AS pv ON m.chrom = pv.vchrom AND m.pos = pv.pos AND m.ref = pv.ref AND m.alt1 = pv.alt
    LEFT JOIN snp_ljb23_ma AS ma ON m.chrom = ma.vchrom AND m.pos = ma.pos AND m.ref = ma.ref AND m.alt1 = ma.alt
    LEFT JOIN snp_ljb23_mt AS mt ON m.chrom = mt.vchrom AND m.pos = mt.pos AND m.ref = mt.ref AND m.alt1 = mt.alt
    LEFT JOIN snp_gwas_catalogs AS g ON m.chrom = g.vchrom AND m.pos = g.pos AND m.ref = g.ref AND m.alt1 = g.alt
    LEFT JOIN snp_phast_cons_elements46way AS c1 ON m.chrom = c1.vchrom AND m.pos = c1.pos AND m.ref = c1.ref AND m.alt1 = c1.alt
    LEFT JOIN snp_phast_cons_elements46way_placental AS c2 ON m.chrom = c2.vchrom AND m.pos = c2.pos AND m.ref = c2.ref AND m.alt1 = c2.alt
    LEFT JOIN snp_phast_cons_elements46way_primates AS c3 ON m.chrom = c3.vchrom AND m.pos = c3.pos AND m.ref = c3.ref AND m.alt1 = c3.alt
    LEFT JOIN snp_tfbs_cons_sites AS t ON m.chrom = t.vchrom AND m.pos = t.pos AND m.ref = t.ref AND m.alt1 = t.alt
    LEFT JOIN snp_wg_rna AS wr ON m.chrom = wr.vchrom AND m.pos = wr.pos AND m.ref = wr.ref AND m.alt1 = wr.alt
    LEFT JOIN snp_target_scans AS ts ON m.chrom = ts.vchrom AND m.pos = ts.pos AND m.ref = ts.ref AND m.alt1 = ts.alt
    LEFT JOIN snp_dgvs AS d ON m.chrom = d.vchrom AND m.pos = d.pos AND m.ref = d.ref AND m.alt1 = d.alt
    LEFT JOIN snp_wg_encode_reg_dnase_clustered_v2 AS dc ON m.chrom = dc.vchrom AND m.pos = dc.pos AND m.ref = dc.ref AND m.alt1 = dc.alt
    LEFT JOIN snp_wg_encode_reg_tfbs_clustered_v3 AS tc ON m.chrom = tc.vchrom AND m.pos = tc.pos AND m.ref = tc.ref AND m.alt1 = tc.alt
    LEFT JOIN snp_wg_encode_broad_hmm_gm12878 AS h1 ON m.chrom = h1.vchrom AND m.pos = h1.pos AND m.ref = h1.ref AND m.alt1 = h1.alt
    LEFT JOIN snp_wg_encode_broad_hmm_h1hesc AS h2 ON m.chrom = h2.vchrom AND m.pos = h2.pos AND m.ref = h2.ref AND m.alt1 = h2.alt
    LEFT JOIN snp_wg_encode_broad_hmm_hepg2 AS h3 ON m.chrom = h3.vchrom AND m.pos = h3.pos AND m.ref = h3.ref AND m.alt1 = h3.alt
    LEFT JOIN snp_wg_encode_broad_hmm_hmec AS h4 ON m.chrom = h4.vchrom AND m.pos = h4.pos AND m.ref = h4.ref AND m.alt1 = h4.alt
    LEFT JOIN snp_wg_encode_broad_hmm_hsmm AS h5 ON m.chrom = h5.vchrom AND m.pos = h5.pos AND m.ref = h5.ref AND m.alt1 = h5.alt
    LEFT JOIN snp_wg_encode_broad_hmm_huvec AS h6 ON m.chrom = h6.vchrom AND m.pos = h6.pos AND m.ref = h6.ref AND m.alt1 = h6.alt
    LEFT JOIN snp_wg_encode_broad_hmm_k562 AS h7 ON m.chrom = h7.vchrom AND m.pos = h7.pos AND m.ref = h7.ref AND m.alt1 = h7.alt
    LEFT JOIN snp_wg_encode_broad_hmm_nhek AS h8 ON m.chrom = h8.vchrom AND m.pos = h8.pos AND m.ref = h8.ref AND m.alt1 = h8.alt
    LEFT JOIN snp_wg_encode_broad_hmm_nhlf AS h9 ON m.chrom = h9.vchrom AND m.pos = h9.pos AND m.ref = h9.ref AND m.alt1 = h9.alt
    LEFT JOIN snp_rmsk AS r ON m.chrom = r.vchrom AND m.pos = r.pos AND m.ref = r.ref AND m.alt1 = r.alt
    LEFT JOIN snp_simple_repeats AS s ON m.chrom = s.vchrom AND m.pos = s.pos AND m.ref = s.ref AND m.alt1 = s.alt
    LEFT JOIN snp_nested_repeats AS n ON m.chrom = n.vchrom AND m.pos = n.pos AND m.ref = n.ref AND m.alt1 = n.alt
    LEFT JOIN snp_microsatelites AS ms ON m.chrom = ms.vchrom AND m.pos = ms.pos AND m.ref = ms.ref AND m.alt1 = ms.alt
    LEFT JOIN snp_combined_dhs_peaks AS d1 ON m.chrom = d1.vchrom AND m.pos = d1.pos AND m.ref = d1.ref AND m.alt1 = d1.alt
    LEFT JOIN snp_distal_dhs_to_promoter_dhs AS d2 ON m.chrom = d2.vchrom AND m.pos = d2.pos AND m.ref = d2.ref AND m.alt1 = d2.alt
    LEFT JOIN snp_dhs_to_gene_expression AS d3 ON m.chrom = d3.vchrom AND m.pos = d3.pos AND m.ref = d3.ref AND m.alt1 = d3.alt
    LEFT JOIN snp_promoter_dhs_master_known AS p1 ON m.chrom = p1.vchrom AND m.pos = p1.pos AND m.ref = p1.ref AND m.alt1 = p1.alt
    LEFT JOIN snp_promoter_dhs_master_novel AS p2 ON m.chrom = p2.vchrom AND m.pos = p2.pos AND m.ref = p2.ref AND m.alt1 = p2.alt;

