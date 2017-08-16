\Problem name: benchmark_LP

Maximize
 obj: biomass_SC5_notrace
Subject To
val-L_c: 1.0*LMPD_72_val-L_c - 1.0*LMPD_72_val-L_c_reverse_768ff + 1.0*VALt2r - 1.0*VALt2r_reverse_e3d24 - 0.2646*biomass_SC5_notrace + 0.2646*biomass_SC5_notrace_reverse_e32ff = 0.0
13BDglcn_c: 1.0*LMPD_1_13BDglcn_c - 1.0*LMPD_1_13BDglcn_c_reverse_e2bdf - 1.1348*biomass_SC5_notrace + 1.1348*biomass_SC5_notrace_reverse_e32ff = 0.0
pro-L_c: 1.0*LMPD_53_pro-L_c - 1.0*LMPD_53_pro-L_c_reverse_830d7 + 1.0*LMPD_54_pro-L_c - 1.0*LMPD_54_pro-L_c_reverse_a6168 + 1.0*PROt2r - 1.0*PROt2r_reverse_f59df - 0.1647*biomass_SC5_notrace + 0.1647*biomass_SC5_notrace_reverse_e32ff = 0.0
ps_SC_c: 1.0*LMPD_48_pc_SC_c - 1.0*LMPD_48_pc_SC_c_reverse_be709 + 1.0*LMPD_49_pc_SC_c - 1.0*LMPD_49_pc_SC_c_reverse_55f95 + 1.0*LMPD_51_pe_SC_c - 1.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*LMPD_56_ps_SC_c - 1.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 1.0*PStm_SC + 1.0*PStm_SC_reverse_a026e - 1.7e-5*biomass_SC5_notrace + 1.7e-5*biomass_SC5_notrace_reverse_e32ff = 0.0
triglyc_SC_c: 1.0*LMPD_65_triglyc_SC_c - 1.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 1.0*LMPD_66_triglyc_SC_c - 1.0*LMPD_66_triglyc_SC_c_reverse_42d98 - 6.6e-5*biomass_SC5_notrace + 6.6e-5*biomass_SC5_notrace_reverse_e32ff = 0.0
damp_c: 1.0*LMPD_13_damp_c - 1.0*LMPD_13_damp_c_reverse_e961c + 1.0*LMPD_14_damp_c - 1.0*LMPD_14_damp_c_reverse_d5737 + 1.0*LMPD_15_damp_c - 1.0*LMPD_15_damp_c_reverse_33cc4 + 1.0*LMPD_16_damp_c - 1.0*LMPD_16_damp_c_reverse_be99f - 0.0036*biomass_SC5_notrace + 0.0036*biomass_SC5_notrace_reverse_e32ff = 0.0
trp-L_c: 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*TRPt2r - 1.0*TRPt2r_reverse_159c4 - 0.0284*biomass_SC5_notrace + 0.0284*biomass_SC5_notrace_reverse_e32ff = 0.0
gly_c: 1.0*GLYt2r - 1.0*GLYt2r_reverse_d768d + 1.0*LMPD_13_damp_c - 1.0*LMPD_13_damp_c_reverse_e961c + 1.0*LMPD_14_damp_c - 1.0*LMPD_14_damp_c_reverse_d5737 + 1.0*LMPD_15_damp_c - 1.0*LMPD_15_damp_c_reverse_33cc4 + 1.0*LMPD_16_damp_c - 1.0*LMPD_16_damp_c_reverse_be99f + 1.0*LMPD_19_dgmp_c - 1.0*LMPD_19_dgmp_c_reverse_7f2b1 + 1.0*LMPD_20_dgmp_c - 1.0*LMPD_20_dgmp_c_reverse_e4f3a + 1.0*LMPD_21_dgmp_c - 1.0*LMPD_21_dgmp_c_reverse_65045 + 1.0*LMPD_22_dgmp_c - 1.0*LMPD_22_dgmp_c_reverse_ebc77 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 1.0*LMPD_28_gly_c - 1.0*LMPD_28_gly_c_reverse_77428 + 1.0*LMPD_2_amp_c - 1.0*LMPD_2_amp_c_reverse_bde95 + 1.0*LMPD_30_gmp_c - 1.0*LMPD_30_gmp_c_reverse_1f79a + 1.0*LMPD_31_gmp_c - 1.0*LMPD_31_gmp_c_reverse_d0c02 + 1.0*LMPD_32_his-L_c - 1.0*LMPD_32_his-L_c_reverse_c5967 + 1.0*LMPD_33_his-L_c - 1.0*LMPD_33_his-L_c_reverse_0eb8b + 1.0*LMPD_34_his-L_c - 1.0*LMPD_34_his-L_c_reverse_e16a9 + 1.0*LMPD_35_his-L_c - 1.0*LMPD_35_his-L_c_reverse_5e623 + 1.0*LMPD_3_amp_c - 1.0*LMPD_3_amp_c_reverse_96e30 + 1.0*LMPD_44_met-L_c - 1.0*LMPD_44_met-L_c_reverse_1646f + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 + 300.0*LMPD_48_pc_SC_c - 300.0*LMPD_48_pc_SC_c_reverse_be709 + 300.0*LMPD_49_pc_SC_c - 300.0*LMPD_49_pc_SC_c_reverse_55f95 + 1.0*LMPD_4_amp_c - 1.0*LMPD_4_amp_c_reverse_7925f + 1.0*LMPD_5_amp_c - 1.0*LMPD_5_amp_c_reverse_8e974 - 0.2904*biomass_SC5_notrace + 0.2904*biomass_SC5_notrace_reverse_e32ff = 0.0
tyr-L_c: 1.0*LMPD_68_tyr-L_c - 1.0*LMPD_68_tyr-L_c_reverse_6983e + 1.0*LMPD_69_tyr-L_c - 1.0*LMPD_69_tyr-L_c_reverse_03f73 + 1.0*TYRt2r - 1.0*TYRt2r_reverse_6f73b - 0.102*biomass_SC5_notrace + 0.102*biomass_SC5_notrace_reverse_e32ff = 0.0
dcmp_c: 1.0*LMPD_17_dcmp_c - 1.0*LMPD_17_dcmp_c_reverse_5c6b0 + 1.0*LMPD_18_dcmp_c - 1.0*LMPD_18_dcmp_c_reverse_15021 - 0.0024*biomass_SC5_notrace + 0.0024*biomass_SC5_notrace_reverse_e32ff = 0.0
ump_c: 1.0*LMPD_70_ump_c - 1.0*LMPD_70_ump_c_reverse_3dcab + 1.0*LMPD_71_ump_c - 1.0*LMPD_71_ump_c_reverse_dc9d0 - 0.0599*biomass_SC5_notrace + 0.0599*biomass_SC5_notrace_reverse_e32ff = 0.0
dgmp_c: 1.0*LMPD_19_dgmp_c - 1.0*LMPD_19_dgmp_c_reverse_7f2b1 + 1.0*LMPD_20_dgmp_c - 1.0*LMPD_20_dgmp_c_reverse_e4f3a + 1.0*LMPD_21_dgmp_c - 1.0*LMPD_21_dgmp_c_reverse_65045 + 1.0*LMPD_22_dgmp_c - 1.0*LMPD_22_dgmp_c_reverse_ebc77 - 0.0024*biomass_SC5_notrace + 0.0024*biomass_SC5_notrace_reverse_e32ff = 0.0
mannan_c: -1.0*MANNANter + 1.0*MANNANter_reverse_6b4bd - 0.8079*biomass_SC5_notrace + 0.8079*biomass_SC5_notrace_reverse_e32ff = 0.0
ala-L_c: -1.0*ALATA_L + 1.0*ALATA_L_reverse_e54ff + 1.0*ALAt2r - 1.0*ALAt2r_reverse_552b7 + 1.0*ALAtmi - 1.0*ALAtmi_reverse_c89ef - 1.0*LMPD_28_gly_c + 1.0*LMPD_28_gly_c_reverse_77428 - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 - 0.4588*biomass_SC5_notrace + 0.4588*biomass_SC5_notrace_reverse_e32ff = 0.0
ptd1ino_SC_c: 1.0*LMPD_57_ptd1ino_SC_c - 1.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 1.0*LMPD_58_ptd1ino_SC_c - 1.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 5.3e-5*biomass_SC5_notrace + 5.3e-5*biomass_SC5_notrace_reverse_e32ff = 0.0
amp_c: 1.0*ACS - 1.0*ACS_reverse_37635 - 1.0*ADK1 + 1.0*ADK1_reverse_a6f90 + 1.0*ATPH1 - 1.0*ATPH1_reverse_23c0c + 1.0*LMPD_10_cmp_c - 1.0*LMPD_10_cmp_c_reverse_3acfa + 1.0*LMPD_11_cmp_c - 1.0*LMPD_11_cmp_c_reverse_98fa8 + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a + 2.0*LMPD_13_damp_c - 2.0*LMPD_13_damp_c_reverse_e961c + 3.0*LMPD_14_damp_c - 3.0*LMPD_14_damp_c_reverse_d5737 + 2.0*LMPD_15_damp_c - 2.0*LMPD_15_damp_c_reverse_33cc4 + 3.0*LMPD_16_damp_c - 3.0*LMPD_16_damp_c_reverse_be99f + 1.0*LMPD_17_dcmp_c - 1.0*LMPD_17_dcmp_c_reverse_5c6b0 + 1.0*LMPD_18_dcmp_c - 1.0*LMPD_18_dcmp_c_reverse_15021 + 2.0*LMPD_19_dgmp_c - 2.0*LMPD_19_dgmp_c_reverse_7f2b1 + 2.0*LMPD_20_dgmp_c - 2.0*LMPD_20_dgmp_c_reverse_e4f3a + 2.0*LMPD_21_dgmp_c - 2.0*LMPD_21_dgmp_c_reverse_65045 + 2.0*LMPD_22_dgmp_c - 2.0*LMPD_22_dgmp_c_reverse_ebc77 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 1.0*LMPD_25_ergst_c - 1.0*LMPD_25_ergst_c_reverse_96537 + 1.0*LMPD_26_ergst_c - 1.0*LMPD_26_ergst_c_reverse_a59e0 + 2.0*LMPD_2_amp_c - 2.0*LMPD_2_amp_c_reverse_bde95 + 2.0*LMPD_30_gmp_c - 2.0*LMPD_30_gmp_c_reverse_1f79a + 2.0*LMPD_31_gmp_c - 2.0*LMPD_31_gmp_c_reverse_d0c02 + 3.0*LMPD_32_his-L_c - 3.0*LMPD_32_his-L_c_reverse_c5967 + 2.0*LMPD_33_his-L_c - 2.0*LMPD_33_his-L_c_reverse_0eb8b + 2.0*LMPD_34_his-L_c - 2.0*LMPD_34_his-L_c_reverse_e16a9 + 3.0*LMPD_35_his-L_c - 3.0*LMPD_35_his-L_c_reverse_5e623 + 2.0*LMPD_3_amp_c - 2.0*LMPD_3_amp_c_reverse_96e30 + 1.0*LMPD_40_lys-L_c - 1.0*LMPD_40_lys-L_c_reverse_4d3a6 + 1.0*LMPD_41_lys-L_c - 1.0*LMPD_41_lys-L_c_reverse_c021a + 1.0*LMPD_43_mannan_c - 1.0*LMPD_43_mannan_c_reverse_7dd82 + 1.0*LMPD_44_met-L_c - 1.0*LMPD_44_met-L_c_reverse_1646f + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 + 300.0*LMPD_48_pc_SC_c - 300.0*LMPD_48_pc_SC_c_reverse_be709 + 300.0*LMPD_49_pc_SC_c - 300.0*LMPD_49_pc_SC_c_reverse_55f95 + 3.0*LMPD_4_amp_c - 3.0*LMPD_4_amp_c_reverse_7925f + 3.0*LMPD_59_ribflv_c - 3.0*LMPD_59_ribflv_c_reverse_01d9b + 3.0*LMPD_5_amp_c - 3.0*LMPD_5_amp_c_reverse_8e974 + 2.0*LMPD_60_ribflv_c - 2.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 + 1.0*LMPD_70_ump_c - 1.0*LMPD_70_ump_c_reverse_3dcab + 1.0*LMPD_71_ump_c - 1.0*LMPD_71_ump_c_reverse_dc9d0 + 1.0*LMPD_7_asn-L_c - 1.0*LMPD_7_asn-L_c_reverse_b41c7 + 1.0*LMPD_8_cmp_c - 1.0*LMPD_8_cmp_c_reverse_23cbb + 1.0*LMPD_9_cmp_c - 1.0*LMPD_9_cmp_c_reverse_85279 + 1.0*NDP1 - 1.0*NDP1_reverse_58712 - 0.046*biomass_SC5_notrace + 0.046*biomass_SC5_notrace_reverse_e32ff = 0.0
pe_SC_c: -1.0*LMPD_48_pc_SC_c + 1.0*LMPD_48_pc_SC_c_reverse_be709 - 1.0*LMPD_49_pc_SC_c + 1.0*LMPD_49_pc_SC_c_reverse_55f95 - 1.0*PEtm_SC + 1.0*PEtm_SC_reverse_a93d3 - 4.5e-5*biomass_SC5_notrace + 4.5e-5*biomass_SC5_notrace_reverse_e32ff = 0.0
arg-L_c: 1.0*ARGt2r - 1.0*ARGt2r_reverse_57e3f + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 - 0.1607*biomass_SC5_notrace + 0.1607*biomass_SC5_notrace_reverse_e32ff = 0.0
met-L_c: 1.0*LMPD_44_met-L_c - 1.0*LMPD_44_met-L_c_reverse_1646f + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 + 1.0*METt2r - 1.0*METt2r_reverse_22908 - 0.0507*biomass_SC5_notrace + 0.0507*biomass_SC5_notrace_reverse_e32ff = 0.0
pa_SC_c: 1.0*LMPD_47_pa_SC_c - 1.0*LMPD_47_pa_SC_c_reverse_4049b - 6.0e-6*biomass_SC5_notrace + 6.0e-6*biomass_SC5_notrace_reverse_e32ff = 0.0
ribflv_c: 1.0*LMPD_59_ribflv_c - 1.0*LMPD_59_ribflv_c_reverse_01d9b + 1.0*LMPD_60_ribflv_c - 1.0*LMPD_60_ribflv_c_reverse_b7a89 - 0.00099*biomass_SC5_notrace + 0.00099*biomass_SC5_notrace_reverse_e32ff = 0.0
dtmp_c: 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 - 0.0036*biomass_SC5_notrace + 0.0036*biomass_SC5_notrace_reverse_e32ff = 0.0
ile-L_c: 1.0*ILEt2r - 1.0*ILEt2r_reverse_9891b + 1.0*ILEtmi - 1.0*ILEtmi_reverse_5c3e7 - 0.1927*biomass_SC5_notrace + 0.1927*biomass_SC5_notrace_reverse_e32ff = 0.0
ergst_c: 1.0*ERGSTt - 1.0*ERGSTt_reverse_c595b + 1.0*ERGSTter - 1.0*ERGSTter_reverse_23bab - 0.0007*biomass_SC5_notrace + 0.0007*biomass_SC5_notrace_reverse_e32ff = 0.0
gmp_c: 1.0*LMPD_30_gmp_c - 1.0*LMPD_30_gmp_c_reverse_1f79a + 1.0*LMPD_31_gmp_c - 1.0*LMPD_31_gmp_c_reverse_d0c02 - 0.046*biomass_SC5_notrace + 0.046*biomass_SC5_notrace_reverse_e32ff = 0.0
gln-L_c: 1.0*GLNt2r - 1.0*GLNt2r_reverse_35d46 + 1.0*LMPD_27_gln-L_c - 1.0*LMPD_27_gln-L_c_reverse_77507 - 0.1054*biomass_SC5_notrace + 0.1054*biomass_SC5_notrace_reverse_e32ff = 0.0
leu-L_c: 1.0*LEUt2r - 1.0*LEUt2r_reverse_4a477 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_39_leu-L_c - 1.0*LMPD_39_leu-L_c_reverse_73152 - 0.2964*biomass_SC5_notrace + 0.2964*biomass_SC5_notrace_reverse_e32ff = 0.0
glu-L_c: 1.0*ALATA_L - 1.0*ALATA_L_reverse_e54ff - 1.0*ASPGLU2m + 1.0*ASPGLU2m_reverse_efaea + 1.0*ASPTA - 1.0*ASPTA_reverse_36525 - 1.0*GLUDC + 1.0*GLUDC_reverse_c4aae - 1.0*GLUDxi + 1.0*GLUDxi_reverse_6fbe2 + 1.0*GLUDyi - 1.0*GLUDyi_reverse_a2351 + 1.0*GLUt2r - 1.0*GLUt2r_reverse_3e88a - 1.0*GLUt7m + 1.0*GLUt7m_reverse_26492 - 1.0*LMPD_12_cys-L_c + 1.0*LMPD_12_cys-L_c_reverse_6dc0a - 2.0*LMPD_13_damp_c + 2.0*LMPD_13_damp_c_reverse_e961c - 2.0*LMPD_14_damp_c + 2.0*LMPD_14_damp_c_reverse_d5737 - 2.0*LMPD_15_damp_c + 2.0*LMPD_15_damp_c_reverse_33cc4 - 2.0*LMPD_16_damp_c + 2.0*LMPD_16_damp_c_reverse_be99f - 2.0*LMPD_19_dgmp_c + 2.0*LMPD_19_dgmp_c_reverse_7f2b1 - 2.0*LMPD_20_dgmp_c + 2.0*LMPD_20_dgmp_c_reverse_e4f3a - 2.0*LMPD_21_dgmp_c + 2.0*LMPD_21_dgmp_c_reverse_65045 - 2.0*LMPD_22_dgmp_c + 2.0*LMPD_22_dgmp_c_reverse_ebc77 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 - 1.0*LMPD_27_gln-L_c + 1.0*LMPD_27_gln-L_c_reverse_77507 - 2.0*LMPD_2_amp_c + 2.0*LMPD_2_amp_c_reverse_bde95 - 2.0*LMPD_30_gmp_c + 2.0*LMPD_30_gmp_c_reverse_1f79a - 2.0*LMPD_31_gmp_c + 2.0*LMPD_31_gmp_c_reverse_d0c02 - 2.0*LMPD_32_his-L_c + 2.0*LMPD_32_his-L_c_reverse_c5967 - 2.0*LMPD_33_his-L_c + 2.0*LMPD_33_his-L_c_reverse_0eb8b - 2.0*LMPD_34_his-L_c + 2.0*LMPD_34_his-L_c_reverse_e16a9 - 2.0*LMPD_35_his-L_c + 2.0*LMPD_35_his-L_c_reverse_5e623 - 1.0*LMPD_38_leu-L_c + 1.0*LMPD_38_leu-L_c_reverse_0e66c - 1.0*LMPD_39_leu-L_c + 1.0*LMPD_39_leu-L_c_reverse_73152 - 2.0*LMPD_3_amp_c + 2.0*LMPD_3_amp_c_reverse_96e30 - 2.0*LMPD_40_lys-L_c + 2.0*LMPD_40_lys-L_c_reverse_4d3a6 - 2.0*LMPD_41_lys-L_c + 2.0*LMPD_41_lys-L_c_reverse_c021a - 1.0*LMPD_44_met-L_c + 1.0*LMPD_44_met-L_c_reverse_1646f - 1.0*LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c_reverse_912f1 - 400.0*LMPD_48_pc_SC_c + 400.0*LMPD_48_pc_SC_c_reverse_be709 - 400.0*LMPD_49_pc_SC_c + 400.0*LMPD_49_pc_SC_c_reverse_55f95 - 2.0*LMPD_4_amp_c + 2.0*LMPD_4_amp_c_reverse_7925f - 100.0*LMPD_51_pe_SC_c + 100.0*LMPD_51_pe_SC_c_reverse_09777 - 1.0*LMPD_52_phe-L_c + 1.0*LMPD_52_phe-L_c_reverse_909a1 - 1.0*LMPD_53_pro-L_c + 1.0*LMPD_53_pro-L_c_reverse_830d7 - 1.0*LMPD_54_pro-L_c + 1.0*LMPD_54_pro-L_c_reverse_a6168 - 100.0*LMPD_56_ps_SC_c + 100.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 2.0*LMPD_5_amp_c + 2.0*LMPD_5_amp_c_reverse_8e974 - 1.0*LMPD_61_ser-L_c + 1.0*LMPD_61_ser-L_c_reverse_40666 - 1.0*LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_68_tyr-L_c + 1.0*LMPD_68_tyr-L_c_reverse_6983e - 1.0*LMPD_69_tyr-L_c + 1.0*LMPD_69_tyr-L_c_reverse_03f73 - 1.0*LMPD_72_val-L_c + 1.0*LMPD_72_val-L_c_reverse_768ff - 0.3018*biomass_SC5_notrace + 0.3018*biomass_SC5_notrace_reverse_e32ff = 0.0
lys-L_c: 1.0*LMPD_40_lys-L_c - 1.0*LMPD_40_lys-L_c_reverse_4d3a6 + 1.0*LMPD_41_lys-L_c - 1.0*LMPD_41_lys-L_c_reverse_c021a + 1.0*LYSt2r - 1.0*LYSt2r_reverse_7e0ef - 0.2862*biomass_SC5_notrace + 0.2862*biomass_SC5_notrace_reverse_e32ff = 0.0
zymst_c: 1.0*LMPD_73_zymst_c - 1.0*LMPD_73_zymst_c_reverse_e6142 + 1.0*LMPD_74_zymst_c - 1.0*LMPD_74_zymst_c_reverse_84617 + 1.0*ZYMSTt - 1.0*ZYMSTt_reverse_ab42f - 0.0015*biomass_SC5_notrace + 0.0015*biomass_SC5_notrace_reverse_e32ff = 0.0
ser-L_c: 1.0*LMPD_61_ser-L_c - 1.0*LMPD_61_ser-L_c_reverse_40666 + 1.0*SERt2r - 1.0*SERt2r_reverse_1b436 - 0.1854*biomass_SC5_notrace + 0.1854*biomass_SC5_notrace_reverse_e32ff = 0.0
asn-L_c: 1.0*ASNt2r - 1.0*ASNt2r_reverse_5bcd2 + 1.0*LMPD_7_asn-L_c - 1.0*LMPD_7_asn-L_c_reverse_b41c7 - 0.1017*biomass_SC5_notrace + 0.1017*biomass_SC5_notrace_reverse_e32ff = 0.0
asp-L_c: 1.0*ASPGLU2m - 1.0*ASPGLU2m_reverse_efaea - 1.0*ASPTA + 1.0*ASPTA_reverse_36525 - 1.0*ASPt2m + 1.0*ASPt2m_reverse_c9585 + 1.0*ASPt2r - 1.0*ASPt2r_reverse_3a681 - 1.0*LMPD_10_cmp_c + 1.0*LMPD_10_cmp_c_reverse_3acfa - 1.0*LMPD_11_cmp_c + 1.0*LMPD_11_cmp_c_reverse_98fa8 - 2.0*LMPD_13_damp_c + 2.0*LMPD_13_damp_c_reverse_e961c - 2.0*LMPD_14_damp_c + 2.0*LMPD_14_damp_c_reverse_d5737 - 2.0*LMPD_15_damp_c + 2.0*LMPD_15_damp_c_reverse_33cc4 - 2.0*LMPD_16_damp_c + 2.0*LMPD_16_damp_c_reverse_be99f - 1.0*LMPD_17_dcmp_c + 1.0*LMPD_17_dcmp_c_reverse_5c6b0 - 1.0*LMPD_18_dcmp_c + 1.0*LMPD_18_dcmp_c_reverse_15021 - 1.0*LMPD_19_dgmp_c + 1.0*LMPD_19_dgmp_c_reverse_7f2b1 - 1.0*LMPD_20_dgmp_c + 1.0*LMPD_20_dgmp_c_reverse_e4f3a - 1.0*LMPD_21_dgmp_c + 1.0*LMPD_21_dgmp_c_reverse_65045 - 1.0*LMPD_22_dgmp_c + 1.0*LMPD_22_dgmp_c_reverse_ebc77 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 - 2.0*LMPD_2_amp_c + 2.0*LMPD_2_amp_c_reverse_bde95 - 1.0*LMPD_30_gmp_c + 1.0*LMPD_30_gmp_c_reverse_1f79a - 1.0*LMPD_31_gmp_c + 1.0*LMPD_31_gmp_c_reverse_d0c02 - 1.0*LMPD_32_his-L_c + 1.0*LMPD_32_his-L_c_reverse_c5967 - 1.0*LMPD_33_his-L_c + 1.0*LMPD_33_his-L_c_reverse_0eb8b - 1.0*LMPD_34_his-L_c + 1.0*LMPD_34_his-L_c_reverse_e16a9 - 1.0*LMPD_35_his-L_c + 1.0*LMPD_35_his-L_c_reverse_5e623 - 1.0*LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 - 2.0*LMPD_3_amp_c + 2.0*LMPD_3_amp_c_reverse_96e30 - 1.0*LMPD_44_met-L_c + 1.0*LMPD_44_met-L_c_reverse_1646f - 1.0*LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c_reverse_912f1 - 2.0*LMPD_4_amp_c + 2.0*LMPD_4_amp_c_reverse_7925f - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 2.0*LMPD_5_amp_c + 2.0*LMPD_5_amp_c_reverse_8e974 - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 - 1.0*LMPD_62_thr-L_c + 1.0*LMPD_62_thr-L_c_reverse_e4d09 - 1.0*LMPD_63_thr-L_c + 1.0*LMPD_63_thr-L_c_reverse_52245 - 1.0*LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c_reverse_33f31 - 1.0*LMPD_70_ump_c + 1.0*LMPD_70_ump_c_reverse_3dcab - 1.0*LMPD_71_ump_c + 1.0*LMPD_71_ump_c_reverse_dc9d0 - 1.0*LMPD_7_asn-L_c + 1.0*LMPD_7_asn-L_c_reverse_b41c7 - 1.0*LMPD_8_cmp_c + 1.0*LMPD_8_cmp_c_reverse_23cbb - 1.0*LMPD_9_cmp_c + 1.0*LMPD_9_cmp_c_reverse_85279 - 0.2975*biomass_SC5_notrace + 0.2975*biomass_SC5_notrace_reverse_e32ff = 0.0
atp_c: -1.0*ACS + 1.0*ACS_reverse_37635 - 1.0*ADK1 + 1.0*ADK1_reverse_a6f90 - 1.0*ATPH1 + 1.0*ATPH1_reverse_23c0c + 1.0*ATPM_reverse_5b752 - 1.0*ATPS + 1.0*ATPS_reverse_64286 + 1.0*ATPtm-H - 1.0*ATPtm-H_reverse_ab831 - 1.0*GLUK + 1.0*GLUK_reverse_ecb52 - 1.0*HEX1 + 1.0*HEX1_reverse_25efa - 5.0*LMPD_10_cmp_c + 5.0*LMPD_10_cmp_c_reverse_3acfa - 5.0*LMPD_11_cmp_c + 5.0*LMPD_11_cmp_c_reverse_98fa8 - 2.0*LMPD_12_cys-L_c + 2.0*LMPD_12_cys-L_c_reverse_6dc0a - 7.0*LMPD_13_damp_c + 7.0*LMPD_13_damp_c_reverse_e961c - 6.0*LMPD_14_damp_c + 6.0*LMPD_14_damp_c_reverse_d5737 - 7.0*LMPD_15_damp_c + 7.0*LMPD_15_damp_c_reverse_33cc4 - 6.0*LMPD_16_damp_c + 6.0*LMPD_16_damp_c_reverse_be99f - 4.0*LMPD_17_dcmp_c + 4.0*LMPD_17_dcmp_c_reverse_5c6b0 - 4.0*LMPD_18_dcmp_c + 4.0*LMPD_18_dcmp_c_reverse_15021 - 10.0*LMPD_19_dgmp_c + 10.0*LMPD_19_dgmp_c_reverse_7f2b1 - 1.0*LMPD_1_13BDglcn_c + 1.0*LMPD_1_13BDglcn_c_reverse_e2bdf - 9.0*LMPD_20_dgmp_c + 9.0*LMPD_20_dgmp_c_reverse_e4f3a - 9.0*LMPD_21_dgmp_c + 9.0*LMPD_21_dgmp_c_reverse_65045 - 10.0*LMPD_22_dgmp_c + 10.0*LMPD_22_dgmp_c_reverse_ebc77 - 4.0*LMPD_23_dtmp_c + 4.0*LMPD_23_dtmp_c_reverse_ed12d - 4.0*LMPD_24_dtmp_c + 4.0*LMPD_24_dtmp_c_reverse_5d241 - 21.0*LMPD_25_ergst_c + 21.0*LMPD_25_ergst_c_reverse_96537 - 21.0*LMPD_26_ergst_c + 21.0*LMPD_26_ergst_c_reverse_a59e0 - 1.0*LMPD_27_gln-L_c + 1.0*LMPD_27_gln-L_c_reverse_77507 - 1.0*LMPD_29_glycogen_c + 1.0*LMPD_29_glycogen_c_reverse_68cb3 - 8.0*LMPD_2_amp_c + 8.0*LMPD_2_amp_c_reverse_bde95 - 9.0*LMPD_30_gmp_c + 9.0*LMPD_30_gmp_c_reverse_1f79a - 9.0*LMPD_31_gmp_c + 9.0*LMPD_31_gmp_c_reverse_d0c02 - 3.0*LMPD_32_his-L_c + 3.0*LMPD_32_his-L_c_reverse_c5967 - 4.0*LMPD_33_his-L_c + 4.0*LMPD_33_his-L_c_reverse_0eb8b - 4.0*LMPD_34_his-L_c + 4.0*LMPD_34_his-L_c_reverse_e16a9 - 3.0*LMPD_35_his-L_c + 3.0*LMPD_35_his-L_c_reverse_5e623 - 2.0*LMPD_36_ile-L_c + 2.0*LMPD_36_ile-L_c_reverse_d2bf7 - 2.0*LMPD_37_ile-L_c + 2.0*LMPD_37_ile-L_c_reverse_1ea52 - 8.0*LMPD_3_amp_c + 8.0*LMPD_3_amp_c_reverse_96e30 - 1.0*LMPD_40_lys-L_c + 1.0*LMPD_40_lys-L_c_reverse_4d3a6 - 1.0*LMPD_41_lys-L_c + 1.0*LMPD_41_lys-L_c_reverse_c021a - 1.0*LMPD_42_mannan_c + 1.0*LMPD_42_mannan_c_reverse_9f626 - 3.0*LMPD_44_met-L_c + 3.0*LMPD_44_met-L_c_reverse_1646f - 3.0*LMPD_45_met-L_c + 3.0*LMPD_45_met-L_c_reverse_912f1 - 1420.0*LMPD_47_pa_SC_c + 1420.0*LMPD_47_pa_SC_c_reverse_4049b - 2220.0*LMPD_48_pc_SC_c + 2220.0*LMPD_48_pc_SC_c_reverse_be709 - 2220.0*LMPD_49_pc_SC_c + 2220.0*LMPD_49_pc_SC_c_reverse_55f95 - 7.0*LMPD_4_amp_c + 7.0*LMPD_4_amp_c_reverse_7925f - 1620.0*LMPD_51_pe_SC_c + 1620.0*LMPD_51_pe_SC_c_reverse_09777 - 1.0*LMPD_52_phe-L_c + 1.0*LMPD_52_phe-L_c_reverse_909a1 - 1.0*LMPD_53_pro-L_c + 1.0*LMPD_53_pro-L_c_reverse_830d7 - 1.0*LMPD_54_pro-L_c + 1.0*LMPD_54_pro-L_c_reverse_a6168 - 1620.0*LMPD_56_ps_SC_c + 1620.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 1620.0*LMPD_57_ptd1ino_SC_c + 1620.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 1620.0*LMPD_58_ptd1ino_SC_c + 1620.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 12.0*LMPD_59_ribflv_c + 12.0*LMPD_59_ribflv_c_reverse_01d9b - 7.0*LMPD_5_amp_c + 7.0*LMPD_5_amp_c_reverse_8e974 - 13.0*LMPD_60_ribflv_c + 13.0*LMPD_60_ribflv_c_reverse_b7a89 - 2.0*LMPD_62_thr-L_c + 2.0*LMPD_62_thr-L_c_reverse_e4d09 - 2.0*LMPD_63_thr-L_c + 2.0*LMPD_63_thr-L_c_reverse_52245 - 1.0*LMPD_64_tre_c + 1.0*LMPD_64_tre_c_reverse_64f54 - 2130.0*LMPD_65_triglyc_SC_c + 2130.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 2130.0*LMPD_66_triglyc_SC_c + 2130.0*LMPD_66_triglyc_SC_c_reverse_42d98 - 3.0*LMPD_67_trp-L_c + 3.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_68_tyr-L_c + 1.0*LMPD_68_tyr-L_c_reverse_6983e - 1.0*LMPD_69_tyr-L_c + 1.0*LMPD_69_tyr-L_c_reverse_03f73 - 4.0*LMPD_6_arg-L_c + 4.0*LMPD_6_arg-L_c_reverse_33f31 - 4.0*LMPD_70_ump_c + 4.0*LMPD_70_ump_c_reverse_3dcab - 4.0*LMPD_71_ump_c + 4.0*LMPD_71_ump_c_reverse_dc9d0 - 18.0*LMPD_73_zymst_c + 18.0*LMPD_73_zymst_c_reverse_e6142 - 18.0*LMPD_74_zymst_c + 18.0*LMPD_74_zymst_c_reverse_84617 - 2.0*LMPD_7_asn-L_c + 2.0*LMPD_7_asn-L_c_reverse_b41c7 - 6.0*LMPD_8_cmp_c + 6.0*LMPD_8_cmp_c_reverse_23cbb - 6.0*LMPD_9_cmp_c + 6.0*LMPD_9_cmp_c_reverse_85279 - 1.0*NADK + 1.0*NADK_reverse_bba52 - 1.0*PC + 1.0*PC_reverse_88dba - 1.0*PFK - 1.0*PFK_3 + 1.0*PFK_3_reverse_037e0 + 1.0*PFK_reverse_d24a6 - 1.0*PGK + 1.0*PGK_reverse_02696 - 1.0*PPCK + 1.0*PPCK_reverse_2557d + 1.0*PYK - 1.0*PYK_reverse_bc8ff - 59.276*biomass_SC5_notrace + 59.276*biomass_SC5_notrace_reverse_e32ff - 1.0*ATPM = 0.0
pc_SC_c: 1.0*LMPD_48_pc_SC_c - 1.0*LMPD_48_pc_SC_c_reverse_be709 + 1.0*LMPD_49_pc_SC_c - 1.0*LMPD_49_pc_SC_c_reverse_55f95 - 6.0e-5*biomass_SC5_notrace + 6.0e-5*biomass_SC5_notrace_reverse_e32ff = 0.0
glycogen_c: 1.0*LMPD_29_glycogen_c - 1.0*LMPD_29_glycogen_c_reverse_68cb3 - 0.5185*biomass_SC5_notrace + 0.5185*biomass_SC5_notrace_reverse_e32ff = 0.0
phe-L_c: 1.0*LMPD_52_phe-L_c - 1.0*LMPD_52_phe-L_c_reverse_909a1 + 1.0*PHEt2r - 1.0*PHEt2r_reverse_19e89 - 0.1339*biomass_SC5_notrace + 0.1339*biomass_SC5_notrace_reverse_e32ff = 0.0
h_c: 1.0*ABUTt2r - 1.0*ABUTt2r_reverse_d925e - 1.0*ACS + 1.0*ACS_reverse_37635 + 1.0*ACt2r - 1.0*ACt2r_reverse_64e05 + 1.0*AKGt2r - 1.0*AKGt2r_reverse_5d500 + 1.0*ALAt2r - 1.0*ALAt2r_reverse_552b7 + 1.0*ALCD2if - 1.0*ALCD2if_reverse_85957 - 1.0*ALCD2ir + 1.0*ALCD2ir_reverse_ba067 + 1.0*ALCD2x - 1.0*ALCD2x_reverse_5d107 + 2.0*ALDD2y - 2.0*ALDD2y_reverse_03afb + 1.0*ARGt2r - 1.0*ARGt2r_reverse_57e3f + 1.0*ASNt2r - 1.0*ASNt2r_reverse_5bcd2 - 1.0*ASPt2m + 1.0*ASPt2m_reverse_c9585 + 1.0*ASPt2r - 1.0*ASPt2r_reverse_3a681 + 2.0*ATPH1 - 2.0*ATPH1_reverse_23c0c - 1.0*ATPM_reverse_5b752 - 3.0*ATPS3m + 3.0*ATPS3m_reverse_5d8d3 - 1.0*ATPtm-H + 1.0*ATPtm-H_reverse_ab831 + 1.0*CITt2r - 1.0*CITt2r_reverse_9f606 + 2.0*CYOOm - 2.0*CYOOm_reverse_d730e + 3.5*CYOR_u6m - 3.5*CYOR_u6m_reverse_ba2b8 + 1.0*CYSt2r - 1.0*CYSt2r_reverse_6167b - 1.0*D-LACtm + 1.0*D-LACtm_reverse_08640 + 1.0*DB4PS - 1.0*DB4PS_reverse_43dd1 + 1.0*FBP - 1.0*FBP_reverse_bf2c9 - 1.0*FDNG + 1.0*FDNG_reverse_f3615 + 1.0*FUMt2r - 1.0*FUMt2r_reverse_d93c0 + 1.0*G6PDH2 - 1.0*G6PDH2_reverse_47bd4 + 1.0*GLNt2r - 1.0*GLNt2r_reverse_35d46 - 1.0*GLUDC + 1.0*GLUDC_reverse_c4aae + 1.0*GLUDxi - 1.0*GLUDxi_reverse_6fbe2 - 1.0*GLUDyi + 1.0*GLUDyi_reverse_a2351 + 1.0*GLUt2r - 1.0*GLUt2r_reverse_3e88a + 1.0*GLYt2r - 1.0*GLYt2r_reverse_d768d + 1.0*HCO3E - 1.0*HCO3E_reverse_97ea5 + 1.0*HISt2r - 1.0*HISt2r_reverse_b00a3 + 1.0*ILEt2r - 1.0*ILEt2r_reverse_9891b + 1.0*LEUt2r - 1.0*LEUt2r_reverse_4a477 + 4.0*LMPD_10_cmp_c - 4.0*LMPD_10_cmp_c_reverse_3acfa + 4.0*LMPD_11_cmp_c - 4.0*LMPD_11_cmp_c_reverse_98fa8 - 2.0*LMPD_12_cys-L_c + 2.0*LMPD_12_cys-L_c_reverse_6dc0a + 16.0*LMPD_13_damp_c - 16.0*LMPD_13_damp_c_reverse_e961c + 16.0*LMPD_14_damp_c - 16.0*LMPD_14_damp_c_reverse_d5737 + 16.0*LMPD_15_damp_c - 16.0*LMPD_15_damp_c_reverse_33cc4 + 16.0*LMPD_16_damp_c - 16.0*LMPD_16_damp_c_reverse_be99f + 2.0*LMPD_17_dcmp_c - 2.0*LMPD_17_dcmp_c_reverse_5c6b0 + 2.0*LMPD_18_dcmp_c - 2.0*LMPD_18_dcmp_c_reverse_15021 + 18.0*LMPD_19_dgmp_c - 18.0*LMPD_19_dgmp_c_reverse_7f2b1 + 1.0*LMPD_1_13BDglcn_c - 1.0*LMPD_1_13BDglcn_c_reverse_e2bdf + 17.0*LMPD_20_dgmp_c - 17.0*LMPD_20_dgmp_c_reverse_e4f3a + 17.0*LMPD_21_dgmp_c - 17.0*LMPD_21_dgmp_c_reverse_65045 + 18.0*LMPD_22_dgmp_c - 18.0*LMPD_22_dgmp_c_reverse_ebc77 + 2.0*LMPD_23_dtmp_c - 2.0*LMPD_23_dtmp_c_reverse_ed12d + 2.0*LMPD_24_dtmp_c - 2.0*LMPD_24_dtmp_c_reverse_5d241 - 17.0*LMPD_25_ergst_c + 17.0*LMPD_25_ergst_c_reverse_96537 - 17.0*LMPD_26_ergst_c + 17.0*LMPD_26_ergst_c_reverse_a59e0 + 1.0*LMPD_27_gln-L_c - 1.0*LMPD_27_gln-L_c_reverse_77507 + 1.0*LMPD_29_glycogen_c - 1.0*LMPD_29_glycogen_c_reverse_68cb3 + 17.0*LMPD_2_amp_c - 17.0*LMPD_2_amp_c_reverse_bde95 + 18.0*LMPD_30_gmp_c - 18.0*LMPD_30_gmp_c_reverse_1f79a + 18.0*LMPD_31_gmp_c - 18.0*LMPD_31_gmp_c_reverse_d0c02 + 11.0*LMPD_32_his-L_c - 11.0*LMPD_32_his-L_c_reverse_c5967 + 11.0*LMPD_33_his-L_c - 11.0*LMPD_33_his-L_c_reverse_0eb8b + 11.0*LMPD_34_his-L_c - 11.0*LMPD_34_his-L_c_reverse_e16a9 + 11.0*LMPD_35_his-L_c - 11.0*LMPD_35_his-L_c_reverse_5e623 - 1.0*LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 17.0*LMPD_3_amp_c - 17.0*LMPD_3_amp_c_reverse_96e30 - 2.0*LMPD_40_lys-L_c + 2.0*LMPD_40_lys-L_c_reverse_4d3a6 - 2.0*LMPD_41_lys-L_c + 2.0*LMPD_41_lys-L_c_reverse_c021a + 1.0*LMPD_42_mannan_c - 1.0*LMPD_42_mannan_c_reverse_9f626 + 1.0*LMPD_43_mannan_c - 1.0*LMPD_43_mannan_c_reverse_7dd82 - 5.0*LMPD_44_met-L_c + 5.0*LMPD_44_met-L_c_reverse_1646f - 5.0*LMPD_45_met-L_c + 5.0*LMPD_45_met-L_c_reverse_912f1 - 3058.0*LMPD_47_pa_SC_c + 3058.0*LMPD_47_pa_SC_c_reverse_4049b - 2258.0*LMPD_48_pc_SC_c + 2258.0*LMPD_48_pc_SC_c_reverse_be709 - 2258.0*LMPD_49_pc_SC_c + 2258.0*LMPD_49_pc_SC_c_reverse_55f95 + 17.0*LMPD_4_amp_c - 17.0*LMPD_4_amp_c_reverse_7925f - 2858.0*LMPD_51_pe_SC_c + 2858.0*LMPD_51_pe_SC_c_reverse_09777 + 2.0*LMPD_52_phe-L_c - 2.0*LMPD_52_phe-L_c_reverse_909a1 - 1.0*LMPD_53_pro-L_c + 1.0*LMPD_53_pro-L_c_reverse_830d7 - 1.0*LMPD_54_pro-L_c + 1.0*LMPD_54_pro-L_c_reverse_a6168 - 2858.0*LMPD_56_ps_SC_c + 2858.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 2958.0*LMPD_57_ptd1ino_SC_c + 2958.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 2958.0*LMPD_58_ptd1ino_SC_c + 2958.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 + 14.0*LMPD_59_ribflv_c - 14.0*LMPD_59_ribflv_c_reverse_01d9b + 17.0*LMPD_5_amp_c - 17.0*LMPD_5_amp_c_reverse_8e974 + 14.0*LMPD_60_ribflv_c - 14.0*LMPD_60_ribflv_c_reverse_b7a89 + 2.0*LMPD_61_ser-L_c - 2.0*LMPD_61_ser-L_c_reverse_40666 - 1.0*LMPD_62_thr-L_c + 1.0*LMPD_62_thr-L_c_reverse_e4d09 - 1.0*LMPD_63_thr-L_c + 1.0*LMPD_63_thr-L_c_reverse_52245 + 1.0*LMPD_64_tre_c - 1.0*LMPD_64_tre_c_reverse_64f54 - 4537.0*LMPD_65_triglyc_SC_c + 4537.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 4537.0*LMPD_66_triglyc_SC_c + 4537.0*LMPD_66_triglyc_SC_c_reverse_42d98 + 6.0*LMPD_67_trp-L_c - 6.0*LMPD_67_trp-L_c_reverse_17a73 + 3.0*LMPD_68_tyr-L_c - 3.0*LMPD_68_tyr-L_c_reverse_6983e + 3.0*LMPD_69_tyr-L_c - 3.0*LMPD_69_tyr-L_c_reverse_03f73 + 3.0*LMPD_6_arg-L_c - 3.0*LMPD_6_arg-L_c_reverse_33f31 + 2.0*LMPD_70_ump_c - 2.0*LMPD_70_ump_c_reverse_3dcab + 2.0*LMPD_71_ump_c - 2.0*LMPD_71_ump_c_reverse_dc9d0 - 16.0*LMPD_73_zymst_c + 16.0*LMPD_73_zymst_c_reverse_e6142 - 16.0*LMPD_74_zymst_c + 16.0*LMPD_74_zymst_c_reverse_84617 + 1.0*LMPD_7_asn-L_c - 1.0*LMPD_7_asn-L_c_reverse_b41c7 + 5.0*LMPD_8_cmp_c - 5.0*LMPD_8_cmp_c_reverse_23cbb + 5.0*LMPD_9_cmp_c - 5.0*LMPD_9_cmp_c_reverse_85279 + 1.0*LYSt2r - 1.0*LYSt2r_reverse_7e0ef + 1.0*MALt2r - 1.0*MALt2r_reverse_bfe2f + 1.0*MDH - 1.0*MDH_reverse_ee52c + 1.0*METt2r - 1.0*METt2r_reverse_22908 - 1.0*NADH2-u6cm + 1.0*NADH2-u6cm_reverse_07733 + 1.0*NADPPPS - 1.0*NADPPPS_reverse_c9c40 + 1.0*NDP1 - 1.0*NDP1_reverse_58712 - 1.0*OAAt2m + 1.0*OAAt2m_reverse_42ec8 + 1.0*PC - 1.0*PC_reverse_88dba - 1.0*PGK + 1.0*PGK_reverse_02696 + 1.0*PGL - 1.0*PGL_reverse_2bb6b + 1.0*PHEt2r - 1.0*PHEt2r_reverse_19e89 - 1.0*PIt2m + 1.0*PIt2m_reverse_cc401 + 1.0*PIt2r - 1.0*PIt2r_reverse_1cd61 + 2.0*PPA - 2.0*PPA_reverse_c5293 - 1.0*PPCK + 1.0*PPCK_reverse_2557d + 1.0*PROt2r - 1.0*PROt2r_reverse_f59df - 1.0*PYRDC + 1.0*PYRDC_reverse_51259 - 1.0*PYRt2m + 1.0*PYRt2m_reverse_fb762 + 1.0*SERt2r - 1.0*SERt2r_reverse_1b436 + 1.0*SUCCt2r - 1.0*SUCCt2r_reverse_e3abe - 1.0*THRt2m + 1.0*THRt2m_reverse_fd5a8 + 1.0*THRt2r - 1.0*THRt2r_reverse_721ff + 1.0*TRPt2r - 1.0*TRPt2r_reverse_159c4 + 1.0*TYRt2r - 1.0*TYRt2r_reverse_6f73b + 1.0*VALt2r - 1.0*VALt2r_reverse_e3d24 + 58.70001*biomass_SC5_notrace - 58.70001*biomass_SC5_notrace_reverse_e32ff + 1.0*ATPM = 0.0
so4_c: -1.0*LMPD_12_cys-L_c + 1.0*LMPD_12_cys-L_c_reverse_6dc0a - 1.0*LMPD_44_met-L_c + 1.0*LMPD_44_met-L_c_reverse_1646f - 1.0*LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c_reverse_912f1 + 1.0*SO4ti - 1.0*SO4ti_reverse_4a93d - 0.02*biomass_SC5_notrace + 0.02*biomass_SC5_notrace_reverse_e32ff = 0.0
adp_c: 2.0*ADK1 - 2.0*ADK1_reverse_a6f90 - 1.0*ATPM_reverse_5b752 + 1.0*ATPS - 1.0*ATPS_reverse_64286 - 1.0*ATPtm-H + 1.0*ATPtm-H_reverse_ab831 + 1.0*GLUK - 1.0*GLUK_reverse_ecb52 + 1.0*HEX1 - 1.0*HEX1_reverse_25efa + 4.0*LMPD_10_cmp_c - 4.0*LMPD_10_cmp_c_reverse_3acfa + 4.0*LMPD_11_cmp_c - 4.0*LMPD_11_cmp_c_reverse_98fa8 + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a + 5.0*LMPD_13_damp_c - 5.0*LMPD_13_damp_c_reverse_e961c + 3.0*LMPD_14_damp_c - 3.0*LMPD_14_damp_c_reverse_d5737 + 5.0*LMPD_15_damp_c - 5.0*LMPD_15_damp_c_reverse_33cc4 + 3.0*LMPD_16_damp_c - 3.0*LMPD_16_damp_c_reverse_be99f + 3.0*LMPD_17_dcmp_c - 3.0*LMPD_17_dcmp_c_reverse_5c6b0 + 3.0*LMPD_18_dcmp_c - 3.0*LMPD_18_dcmp_c_reverse_15021 + 8.0*LMPD_19_dgmp_c - 8.0*LMPD_19_dgmp_c_reverse_7f2b1 + 1.0*LMPD_1_13BDglcn_c - 1.0*LMPD_1_13BDglcn_c_reverse_e2bdf + 7.0*LMPD_20_dgmp_c - 7.0*LMPD_20_dgmp_c_reverse_e4f3a + 7.0*LMPD_21_dgmp_c - 7.0*LMPD_21_dgmp_c_reverse_65045 + 8.0*LMPD_22_dgmp_c - 8.0*LMPD_22_dgmp_c_reverse_ebc77 + 3.0*LMPD_23_dtmp_c - 3.0*LMPD_23_dtmp_c_reverse_ed12d + 3.0*LMPD_24_dtmp_c - 3.0*LMPD_24_dtmp_c_reverse_5d241 + 20.0*LMPD_25_ergst_c - 20.0*LMPD_25_ergst_c_reverse_96537 + 20.0*LMPD_26_ergst_c - 20.0*LMPD_26_ergst_c_reverse_a59e0 + 1.0*LMPD_27_gln-L_c - 1.0*LMPD_27_gln-L_c_reverse_77507 + 1.0*LMPD_29_glycogen_c - 1.0*LMPD_29_glycogen_c_reverse_68cb3 + 7.0*LMPD_2_amp_c - 7.0*LMPD_2_amp_c_reverse_bde95 + 7.0*LMPD_30_gmp_c - 7.0*LMPD_30_gmp_c_reverse_1f79a + 7.0*LMPD_31_gmp_c - 7.0*LMPD_31_gmp_c_reverse_d0c02 + 2.0*LMPD_33_his-L_c - 2.0*LMPD_33_his-L_c_reverse_0eb8b + 2.0*LMPD_34_his-L_c - 2.0*LMPD_34_his-L_c_reverse_e16a9 + 2.0*LMPD_36_ile-L_c - 2.0*LMPD_36_ile-L_c_reverse_d2bf7 + 2.0*LMPD_37_ile-L_c - 2.0*LMPD_37_ile-L_c_reverse_1ea52 + 7.0*LMPD_3_amp_c - 7.0*LMPD_3_amp_c_reverse_96e30 + 1.0*LMPD_42_mannan_c - 1.0*LMPD_42_mannan_c_reverse_9f626 - 1.0*LMPD_43_mannan_c + 1.0*LMPD_43_mannan_c_reverse_7dd82 + 2.0*LMPD_44_met-L_c - 2.0*LMPD_44_met-L_c_reverse_1646f + 2.0*LMPD_45_met-L_c - 2.0*LMPD_45_met-L_c_reverse_912f1 + 1420.0*LMPD_47_pa_SC_c - 1420.0*LMPD_47_pa_SC_c_reverse_4049b + 1920.0*LMPD_48_pc_SC_c - 1920.0*LMPD_48_pc_SC_c_reverse_be709 + 1920.0*LMPD_49_pc_SC_c - 1920.0*LMPD_49_pc_SC_c_reverse_55f95 + 5.0*LMPD_4_amp_c - 5.0*LMPD_4_amp_c_reverse_7925f + 1620.0*LMPD_51_pe_SC_c - 1620.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*LMPD_52_phe-L_c - 1.0*LMPD_52_phe-L_c_reverse_909a1 + 1.0*LMPD_53_pro-L_c - 1.0*LMPD_53_pro-L_c_reverse_830d7 + 1.0*LMPD_54_pro-L_c - 1.0*LMPD_54_pro-L_c_reverse_a6168 + 1620.0*LMPD_56_ps_SC_c - 1620.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 1620.0*LMPD_57_ptd1ino_SC_c - 1620.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 1620.0*LMPD_58_ptd1ino_SC_c - 1620.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 + 9.0*LMPD_59_ribflv_c - 9.0*LMPD_59_ribflv_c_reverse_01d9b + 5.0*LMPD_5_amp_c - 5.0*LMPD_5_amp_c_reverse_8e974 + 11.0*LMPD_60_ribflv_c - 11.0*LMPD_60_ribflv_c_reverse_b7a89 + 2.0*LMPD_62_thr-L_c - 2.0*LMPD_62_thr-L_c_reverse_e4d09 + 2.0*LMPD_63_thr-L_c - 2.0*LMPD_63_thr-L_c_reverse_52245 + 1.0*LMPD_64_tre_c - 1.0*LMPD_64_tre_c_reverse_64f54 + 2130.0*LMPD_65_triglyc_SC_c - 2130.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 2130.0*LMPD_66_triglyc_SC_c - 2130.0*LMPD_66_triglyc_SC_c_reverse_42d98 + 2.0*LMPD_67_trp-L_c - 2.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*LMPD_68_tyr-L_c - 1.0*LMPD_68_tyr-L_c_reverse_6983e + 1.0*LMPD_69_tyr-L_c - 1.0*LMPD_69_tyr-L_c_reverse_03f73 + 3.0*LMPD_6_arg-L_c - 3.0*LMPD_6_arg-L_c_reverse_33f31 + 3.0*LMPD_70_ump_c - 3.0*LMPD_70_ump_c_reverse_3dcab + 3.0*LMPD_71_ump_c - 3.0*LMPD_71_ump_c_reverse_dc9d0 + 18.0*LMPD_73_zymst_c - 18.0*LMPD_73_zymst_c_reverse_e6142 + 18.0*LMPD_74_zymst_c - 18.0*LMPD_74_zymst_c_reverse_84617 + 1.0*LMPD_7_asn-L_c - 1.0*LMPD_7_asn-L_c_reverse_b41c7 + 5.0*LMPD_8_cmp_c - 5.0*LMPD_8_cmp_c_reverse_23cbb + 5.0*LMPD_9_cmp_c - 5.0*LMPD_9_cmp_c_reverse_85279 + 1.0*NADK - 1.0*NADK_reverse_bba52 - 1.0*NDP1 + 1.0*NDP1_reverse_58712 + 1.0*PC - 1.0*PC_reverse_88dba + 1.0*PFK + 1.0*PFK_3 - 1.0*PFK_3_reverse_037e0 - 1.0*PFK_reverse_d24a6 + 1.0*PGK - 1.0*PGK_reverse_02696 + 1.0*PPCK - 1.0*PPCK_reverse_2557d - 1.0*PYK + 1.0*PYK_reverse_bc8ff + 59.276*biomass_SC5_notrace - 59.276*biomass_SC5_notrace_reverse_e32ff + 1.0*ATPM = 0.0
thr-L_c: 1.0*LMPD_36_ile-L_c - 1.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*LMPD_62_thr-L_c - 1.0*LMPD_62_thr-L_c_reverse_e4d09 + 1.0*LMPD_63_thr-L_c - 1.0*LMPD_63_thr-L_c_reverse_52245 - 1.0*THRt2m + 1.0*THRt2m_reverse_fd5a8 + 1.0*THRt2r - 1.0*THRt2r_reverse_721ff - 0.1914*biomass_SC5_notrace + 0.1914*biomass_SC5_notrace_reverse_e32ff = 0.0
cmp_c: 1.0*LMPD_10_cmp_c - 1.0*LMPD_10_cmp_c_reverse_3acfa + 1.0*LMPD_11_cmp_c - 1.0*LMPD_11_cmp_c_reverse_98fa8 + 1.0*LMPD_8_cmp_c - 1.0*LMPD_8_cmp_c_reverse_23cbb + 1.0*LMPD_9_cmp_c - 1.0*LMPD_9_cmp_c_reverse_85279 - 0.0447*biomass_SC5_notrace + 0.0447*biomass_SC5_notrace_reverse_e32ff = 0.0
tre_c: 1.0*LMPD_64_tre_c - 1.0*LMPD_64_tre_c_reverse_64f54 - 0.0234*biomass_SC5_notrace + 0.0234*biomass_SC5_notrace_reverse_e32ff = 0.0
his-L_c: 1.0*HISt2r - 1.0*HISt2r_reverse_b00a3 + 1.0*LMPD_32_his-L_c - 1.0*LMPD_32_his-L_c_reverse_c5967 + 1.0*LMPD_33_his-L_c - 1.0*LMPD_33_his-L_c_reverse_0eb8b + 1.0*LMPD_34_his-L_c - 1.0*LMPD_34_his-L_c_reverse_e16a9 + 1.0*LMPD_35_his-L_c - 1.0*LMPD_35_his-L_c_reverse_5e623 - 0.0663*biomass_SC5_notrace + 0.0663*biomass_SC5_notrace_reverse_e32ff = 0.0
h2o_c: -1.0*ALDD2y + 1.0*ALDD2y_reverse_03afb - 2.0*ATPH1 + 2.0*ATPH1_reverse_23c0c + 1.0*ATPM_reverse_5b752 - 1.0*ATPS + 1.0*ATPS_reverse_64286 + 1.0*ENO - 1.0*ENO_reverse_40eea - 1.0*FBP + 1.0*FBP_reverse_bf2c9 - 1.0*FUM + 1.0*FUM_reverse_d3642 - 1.0*GLUDxi + 1.0*GLUDxi_reverse_6fbe2 + 1.0*GLUDyi - 1.0*GLUDyi_reverse_a2351 + 1.0*H2Ot - 1.0*H2Ot_reverse_aa560 - 1.0*H2Oter + 1.0*H2Oter_reverse_c874d - 1.0*H2Otm + 1.0*H2Otm_reverse_d5f1a - 1.0*HCO3E + 1.0*HCO3E_reverse_97ea5 + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a - 2.0*LMPD_13_damp_c + 2.0*LMPD_13_damp_c_reverse_e961c - 2.0*LMPD_14_damp_c + 2.0*LMPD_14_damp_c_reverse_d5737 - 2.0*LMPD_15_damp_c + 2.0*LMPD_15_damp_c_reverse_33cc4 - 2.0*LMPD_16_damp_c + 2.0*LMPD_16_damp_c_reverse_be99f + 2.0*LMPD_17_dcmp_c - 2.0*LMPD_17_dcmp_c_reverse_5c6b0 + 2.0*LMPD_18_dcmp_c - 2.0*LMPD_18_dcmp_c_reverse_15021 - 5.0*LMPD_19_dgmp_c + 5.0*LMPD_19_dgmp_c_reverse_7f2b1 - 4.0*LMPD_20_dgmp_c + 4.0*LMPD_20_dgmp_c_reverse_e4f3a - 4.0*LMPD_21_dgmp_c + 4.0*LMPD_21_dgmp_c_reverse_65045 - 5.0*LMPD_22_dgmp_c + 5.0*LMPD_22_dgmp_c_reverse_ebc77 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 7.0*LMPD_25_ergst_c - 7.0*LMPD_25_ergst_c_reverse_96537 + 7.0*LMPD_26_ergst_c - 7.0*LMPD_26_ergst_c_reverse_a59e0 - 3.0*LMPD_2_amp_c + 3.0*LMPD_2_amp_c_reverse_bde95 - 5.0*LMPD_30_gmp_c + 5.0*LMPD_30_gmp_c_reverse_1f79a - 5.0*LMPD_31_gmp_c + 5.0*LMPD_31_gmp_c_reverse_d0c02 - 3.0*LMPD_32_his-L_c + 3.0*LMPD_32_his-L_c_reverse_c5967 - 3.0*LMPD_33_his-L_c + 3.0*LMPD_33_his-L_c_reverse_0eb8b - 3.0*LMPD_34_his-L_c + 3.0*LMPD_34_his-L_c_reverse_e16a9 - 3.0*LMPD_35_his-L_c + 3.0*LMPD_35_his-L_c_reverse_5e623 - 1.0*LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 - 1.0*LMPD_38_leu-L_c + 1.0*LMPD_38_leu-L_c_reverse_0e66c - 3.0*LMPD_3_amp_c + 3.0*LMPD_3_amp_c_reverse_96e30 + 2.0*LMPD_44_met-L_c - 2.0*LMPD_44_met-L_c_reverse_1646f + 2.0*LMPD_45_met-L_c - 2.0*LMPD_45_met-L_c_reverse_912f1 + 1656.0*LMPD_47_pa_SC_c - 1656.0*LMPD_47_pa_SC_c_reverse_4049b + 956.0*LMPD_48_pc_SC_c - 956.0*LMPD_48_pc_SC_c_reverse_be709 + 956.0*LMPD_49_pc_SC_c - 956.0*LMPD_49_pc_SC_c_reverse_55f95 - 3.0*LMPD_4_amp_c + 3.0*LMPD_4_amp_c_reverse_7925f + 1556.0*LMPD_51_pe_SC_c - 1556.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*LMPD_52_phe-L_c - 1.0*LMPD_52_phe-L_c_reverse_909a1 + 1.0*LMPD_53_pro-L_c - 1.0*LMPD_53_pro-L_c_reverse_830d7 + 1.0*LMPD_54_pro-L_c - 1.0*LMPD_54_pro-L_c_reverse_a6168 + 1556.0*LMPD_56_ps_SC_c - 1556.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 1556.0*LMPD_57_ptd1ino_SC_c - 1556.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 1556.0*LMPD_58_ptd1ino_SC_c - 1556.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 4.0*LMPD_59_ribflv_c + 4.0*LMPD_59_ribflv_c_reverse_01d9b - 3.0*LMPD_5_amp_c + 3.0*LMPD_5_amp_c_reverse_8e974 - 4.0*LMPD_60_ribflv_c + 4.0*LMPD_60_ribflv_c_reverse_b7a89 - 1.0*LMPD_61_ser-L_c + 1.0*LMPD_61_ser-L_c_reverse_40666 - 1.0*LMPD_62_thr-L_c + 1.0*LMPD_62_thr-L_c_reverse_e4d09 - 1.0*LMPD_63_thr-L_c + 1.0*LMPD_63_thr-L_c_reverse_52245 - 1.0*LMPD_64_tre_c + 1.0*LMPD_64_tre_c_reverse_64f54 + 2384.0*LMPD_65_triglyc_SC_c - 2384.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 2384.0*LMPD_66_triglyc_SC_c - 2384.0*LMPD_66_triglyc_SC_c_reverse_42d98 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c_reverse_33f31 + 6.0*LMPD_73_zymst_c - 6.0*LMPD_73_zymst_c_reverse_e6142 + 6.0*LMPD_74_zymst_c - 6.0*LMPD_74_zymst_c_reverse_84617 - 1.0*LMPD_7_asn-L_c + 1.0*LMPD_7_asn-L_c_reverse_b41c7 - 1.0*LMPD_8_cmp_c + 1.0*LMPD_8_cmp_c_reverse_23cbb - 1.0*LMPD_9_cmp_c + 1.0*LMPD_9_cmp_c_reverse_85279 - 1.0*NADPPPS + 1.0*NADPPPS_reverse_c9c40 - 1.0*NDP1 + 1.0*NDP1_reverse_58712 - 1.0*PGL + 1.0*PGL_reverse_2bb6b - 1.0*PPA + 1.0*PPA_reverse_c5293 - 59.276*biomass_SC5_notrace + 59.276*biomass_SC5_notrace_reverse_e32ff - 1.0*ATPM = 0.0
cys-L_c: 1.0*CYSt2r - 1.0*CYSt2r_reverse_6167b + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a - 0.0066*biomass_SC5_notrace + 0.0066*biomass_SC5_notrace_reverse_e32ff = 0.0
pi_c: 2.0*ATPH1 - 2.0*ATPH1_reverse_23c0c - 1.0*ATPM_reverse_5b752 + 1.0*ATPS - 1.0*ATPS_reverse_64286 + 1.0*FBP - 1.0*FBP_reverse_bf2c9 - 1.0*GAPD + 1.0*GAPD_reverse_459c1 + 4.0*LMPD_10_cmp_c - 4.0*LMPD_10_cmp_c_reverse_3acfa + 4.0*LMPD_11_cmp_c - 4.0*LMPD_11_cmp_c_reverse_98fa8 + 2.0*LMPD_12_cys-L_c - 2.0*LMPD_12_cys-L_c_reverse_6dc0a + 9.0*LMPD_13_damp_c - 9.0*LMPD_13_damp_c_reverse_e961c + 9.0*LMPD_14_damp_c - 9.0*LMPD_14_damp_c_reverse_d5737 + 9.0*LMPD_15_damp_c - 9.0*LMPD_15_damp_c_reverse_33cc4 + 9.0*LMPD_16_damp_c - 9.0*LMPD_16_damp_c_reverse_be99f + 3.0*LMPD_17_dcmp_c - 3.0*LMPD_17_dcmp_c_reverse_5c6b0 + 3.0*LMPD_18_dcmp_c - 3.0*LMPD_18_dcmp_c_reverse_15021 + 10.0*LMPD_19_dgmp_c - 10.0*LMPD_19_dgmp_c_reverse_7f2b1 + 9.0*LMPD_20_dgmp_c - 9.0*LMPD_20_dgmp_c_reverse_e4f3a + 9.0*LMPD_21_dgmp_c - 9.0*LMPD_21_dgmp_c_reverse_65045 + 10.0*LMPD_22_dgmp_c - 10.0*LMPD_22_dgmp_c_reverse_ebc77 + 4.0*LMPD_23_dtmp_c - 4.0*LMPD_23_dtmp_c_reverse_ed12d + 4.0*LMPD_24_dtmp_c - 4.0*LMPD_24_dtmp_c_reverse_5d241 + 8.0*LMPD_25_ergst_c - 8.0*LMPD_25_ergst_c_reverse_96537 + 8.0*LMPD_26_ergst_c - 8.0*LMPD_26_ergst_c_reverse_a59e0 + 1.0*LMPD_27_gln-L_c - 1.0*LMPD_27_gln-L_c_reverse_77507 + 9.0*LMPD_2_amp_c - 9.0*LMPD_2_amp_c_reverse_bde95 + 9.0*LMPD_30_gmp_c - 9.0*LMPD_30_gmp_c_reverse_1f79a + 9.0*LMPD_31_gmp_c - 9.0*LMPD_31_gmp_c_reverse_d0c02 + 4.0*LMPD_32_his-L_c - 4.0*LMPD_32_his-L_c_reverse_c5967 + 4.0*LMPD_33_his-L_c - 4.0*LMPD_33_his-L_c_reverse_0eb8b + 4.0*LMPD_34_his-L_c - 4.0*LMPD_34_his-L_c_reverse_e16a9 + 4.0*LMPD_35_his-L_c - 4.0*LMPD_35_his-L_c_reverse_5e623 + 2.0*LMPD_36_ile-L_c - 2.0*LMPD_36_ile-L_c_reverse_d2bf7 + 2.0*LMPD_37_ile-L_c - 2.0*LMPD_37_ile-L_c_reverse_1ea52 + 9.0*LMPD_3_amp_c - 9.0*LMPD_3_amp_c_reverse_96e30 + 3.0*LMPD_44_met-L_c - 3.0*LMPD_44_met-L_c_reverse_1646f + 3.0*LMPD_45_met-L_c - 3.0*LMPD_45_met-L_c_reverse_912f1 + 1420.0*LMPD_47_pa_SC_c - 1420.0*LMPD_47_pa_SC_c_reverse_4049b + 2120.0*LMPD_48_pc_SC_c - 2120.0*LMPD_48_pc_SC_c_reverse_be709 + 2120.0*LMPD_49_pc_SC_c - 2120.0*LMPD_49_pc_SC_c_reverse_55f95 + 9.0*LMPD_4_amp_c - 9.0*LMPD_4_amp_c_reverse_7925f + 1520.0*LMPD_51_pe_SC_c - 1520.0*LMPD_51_pe_SC_c_reverse_09777 + 4.0*LMPD_52_phe-L_c - 4.0*LMPD_52_phe-L_c_reverse_909a1 + 1.0*LMPD_53_pro-L_c - 1.0*LMPD_53_pro-L_c_reverse_830d7 + 1.0*LMPD_54_pro-L_c - 1.0*LMPD_54_pro-L_c_reverse_a6168 + 1520.0*LMPD_56_ps_SC_c - 1520.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 1520.0*LMPD_57_ptd1ino_SC_c - 1520.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 1520.0*LMPD_58_ptd1ino_SC_c - 1520.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 + 12.0*LMPD_59_ribflv_c - 12.0*LMPD_59_ribflv_c_reverse_01d9b + 9.0*LMPD_5_amp_c - 9.0*LMPD_5_amp_c_reverse_8e974 + 12.0*LMPD_60_ribflv_c - 12.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_61_ser-L_c - 1.0*LMPD_61_ser-L_c_reverse_40666 + 2.0*LMPD_62_thr-L_c - 2.0*LMPD_62_thr-L_c_reverse_e4d09 + 2.0*LMPD_63_thr-L_c - 2.0*LMPD_63_thr-L_c_reverse_52245 + 1.0*LMPD_64_tre_c - 1.0*LMPD_64_tre_c_reverse_64f54 + 2230.0*LMPD_65_triglyc_SC_c - 2230.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 2230.0*LMPD_66_triglyc_SC_c - 2230.0*LMPD_66_triglyc_SC_c_reverse_42d98 + 6.0*LMPD_67_trp-L_c - 6.0*LMPD_67_trp-L_c_reverse_17a73 + 4.0*LMPD_68_tyr-L_c - 4.0*LMPD_68_tyr-L_c_reverse_6983e + 4.0*LMPD_69_tyr-L_c - 4.0*LMPD_69_tyr-L_c_reverse_03f73 + 3.0*LMPD_6_arg-L_c - 3.0*LMPD_6_arg-L_c_reverse_33f31 + 3.0*LMPD_70_ump_c - 3.0*LMPD_70_ump_c_reverse_3dcab + 3.0*LMPD_71_ump_c - 3.0*LMPD_71_ump_c_reverse_dc9d0 + 6.0*LMPD_73_zymst_c - 6.0*LMPD_73_zymst_c_reverse_e6142 + 6.0*LMPD_74_zymst_c - 6.0*LMPD_74_zymst_c_reverse_84617 + 1.0*LMPD_7_asn-L_c - 1.0*LMPD_7_asn-L_c_reverse_b41c7 + 5.0*LMPD_8_cmp_c - 5.0*LMPD_8_cmp_c_reverse_23cbb + 5.0*LMPD_9_cmp_c - 5.0*LMPD_9_cmp_c_reverse_85279 + 1.0*MALtm - 1.0*MALtm_reverse_f798a + 1.0*NADPPPS - 1.0*NADPPPS_reverse_c9c40 + 1.0*NDP1 - 1.0*NDP1_reverse_58712 + 1.0*PC - 1.0*PC_reverse_88dba - 1.0*PIt2m + 1.0*PIt2m_reverse_cc401 + 1.0*PIt2r - 1.0*PIt2r_reverse_1cd61 + 2.0*PPA - 2.0*PPA_reverse_c5293 + 1.0*SUCCtm - 1.0*SUCCtm_reverse_5f36b + 59.305*biomass_SC5_notrace - 59.305*biomass_SC5_notrace_reverse_e32ff + 1.0*ATPM = 0.0
4abut_e: -1.0*ABUTt2r + 1.0*ABUTt2r_reverse_d925e - 1.0*EX_4abut_e + 1.0*EX_4abut_e_reverse_82295 = 0.0
ac_e: -1.0*ACt2r + 1.0*ACt2r_reverse_64e05 - 1.0*ACtr + 1.0*ACtr_reverse_b9231 - 1.0*EX_ac_e + 1.0*EX_ac_e_reverse_0be96 = 0.0
acald_e: -1.0*ACALDt + 1.0*ACALDt_reverse_858fa - 1.0*EX_acald_e + 1.0*EX_acald_e_reverse_c096e = 0.0
akg_e: 1.0*AKGMAL - 1.0*AKGMAL_reverse_76d29 - 1.0*AKGt2r + 1.0*AKGt2r_reverse_5d500 - 1.0*EX_akg_e + 1.0*EX_akg_e_reverse_70d85 = 0.0
ala-L_e: -1.0*ALAt2r + 1.0*ALAt2r_reverse_552b7 - 1.0*EX_ala-L_e + 1.0*EX_ala-L_e_reverse_45174 = 0.0
arg-L_e: -1.0*ARGt2r + 1.0*ARGt2r_reverse_57e3f - 1.0*EX_arg-L_e + 1.0*EX_arg-L_e_reverse_2768b = 0.0
asn-L_e: -1.0*ASNt2r + 1.0*ASNt2r_reverse_5bcd2 - 1.0*EX_asn-L_e + 1.0*EX_asn-L_e_reverse_73ed2 = 0.0
asp-L_e: -1.0*ASPt2r + 1.0*ASPt2r_reverse_3a681 - 1.0*EX_asp-L_e + 1.0*EX_asp-L_e_reverse_ceee6 = 0.0
cit_e: -1.0*CITt2r + 1.0*CITt2r_reverse_9f606 - 1.0*EX_cit_e + 1.0*EX_cit_e_reverse_0835e = 0.0
co2_e: -1.0*CO2t + 1.0*CO2t_reverse_7c42f - 1.0*EX_co2_e + 1.0*EX_co2_e_reverse_d0466 = 0.0
cys-L_e: -1.0*CYSt2r + 1.0*CYSt2r_reverse_6167b - 1.0*EX_cys-L_e + 1.0*EX_cys-L_e_reverse_c0cbc = 0.0
ergst_e: -1.0*ERGSTt + 1.0*ERGSTt_reverse_c595b - 1.0*EX_ergst_e + 1.0*EX_ergst_e_reverse_92efd = 0.0
etoh_e: -1.0*ETOHt + 1.0*ETOHt_reverse_a82a5 - 1.0*EX_etoh_e + 1.0*EX_etoh_e_reverse_cc64f = 0.0
for_e: -1.0*EX_for_e + 1.0*EX_for_e_reverse_23269 - 1.0*FORt + 1.0*FORt_reverse_40f9f = 0.0
fum_e: -1.0*EX_fum_e + 1.0*EX_fum_e_reverse_e3432 - 1.0*FUMt2r + 1.0*FUMt2r_reverse_d93c0 = 0.0
glc-D_e: -1.0*EX_glc_e + 1.0*EX_glc_e_reverse_6226d - 1.0*GLCt1 + 1.0*GLCt1_reverse_3bcc5 = 0.0
gln-L_e: -1.0*EX_gln-L_e + 1.0*EX_gln-L_e_reverse_1be13 - 1.0*GLNt2r + 1.0*GLNt2r_reverse_35d46 = 0.0
glu-L_e: -1.0*EX_glu-L_e + 1.0*EX_glu-L_e_reverse_727e6 - 1.0*GLUt2r + 1.0*GLUt2r_reverse_3e88a = 0.0
gly_e: -1.0*EX_gly_e + 1.0*EX_gly_e_reverse_6956b - 1.0*GLYt2r + 1.0*GLYt2r_reverse_d768d = 0.0
h_e: -1.0*ABUTt2r + 1.0*ABUTt2r_reverse_d925e - 1.0*ACt2r + 1.0*ACt2r_reverse_64e05 - 1.0*AKGt2r + 1.0*AKGt2r_reverse_5d500 - 1.0*ALAt2r + 1.0*ALAt2r_reverse_552b7 - 1.0*ARGt2r + 1.0*ARGt2r_reverse_57e3f - 1.0*ASNt2r + 1.0*ASNt2r_reverse_5bcd2 - 1.0*ASPt2r + 1.0*ASPt2r_reverse_3a681 + 1.0*ATPS - 1.0*ATPS_reverse_64286 - 1.0*CITt2r + 1.0*CITt2r_reverse_9f606 - 1.0*CYSt2r + 1.0*CYSt2r_reverse_6167b - 1.0*EX_h_e + 1.0*EX_h_e_reverse_3e0c5 - 1.0*FUMt2r + 1.0*FUMt2r_reverse_d93c0 - 1.0*GLNt2r + 1.0*GLNt2r_reverse_35d46 - 1.0*GLUt2r + 1.0*GLUt2r_reverse_3e88a - 1.0*GLYt2r + 1.0*GLYt2r_reverse_d768d - 1.0*HISt2r + 1.0*HISt2r_reverse_b00a3 - 1.0*ILEt2r + 1.0*ILEt2r_reverse_9891b - 1.0*LEUt2r + 1.0*LEUt2r_reverse_4a477 - 1.0*LYSt2r + 1.0*LYSt2r_reverse_7e0ef - 1.0*MALt2r + 1.0*MALt2r_reverse_bfe2f - 1.0*METt2r + 1.0*METt2r_reverse_22908 - 1.0*PHEt2r + 1.0*PHEt2r_reverse_19e89 - 1.0*PIt2r + 1.0*PIt2r_reverse_1cd61 - 1.0*PROt2r + 1.0*PROt2r_reverse_f59df - 1.0*SERt2r + 1.0*SERt2r_reverse_1b436 - 1.0*SUCCt2r + 1.0*SUCCt2r_reverse_e3abe - 1.0*THRt2r + 1.0*THRt2r_reverse_721ff - 1.0*TRPt2r + 1.0*TRPt2r_reverse_159c4 - 1.0*TYRt2r + 1.0*TYRt2r_reverse_6f73b - 1.0*VALt2r + 1.0*VALt2r_reverse_e3d24 = 0.0
h2o_e: -1.0*EX_h2o_e + 1.0*EX_h2o_e_reverse_3ced4 - 1.0*H2Ot + 1.0*H2Ot_reverse_aa560 = 0.0
his-L_e: -1.0*EX_his-L_e + 1.0*EX_his-L_e_reverse_e1c17 - 1.0*HISt2r + 1.0*HISt2r_reverse_b00a3 = 0.0
ile-L_e: -1.0*EX_ile-L_e + 1.0*EX_ile-L_e_reverse_e9c53 - 1.0*ILEt2r + 1.0*ILEt2r_reverse_9891b = 0.0
leu-L_e: -1.0*EX_leu-L_e + 1.0*EX_leu-L_e_reverse_9947f - 1.0*LEUt2r + 1.0*LEUt2r_reverse_4a477 = 0.0
lys-L_e: -1.0*EX_lys-L_e + 1.0*EX_lys-L_e_reverse_5779d - 1.0*LYSt2r + 1.0*LYSt2r_reverse_7e0ef = 0.0
mal-L_e: -1.0*AKGMAL + 1.0*AKGMAL_reverse_76d29 - 1.0*EX_mal-L_e + 1.0*EX_mal-L_e_reverse_00553 - 1.0*MALt2r + 1.0*MALt2r_reverse_bfe2f = 0.0
met-L_e: -1.0*EX_met-L_e + 1.0*EX_met-L_e_reverse_f8e12 - 1.0*METt2r + 1.0*METt2r_reverse_22908 = 0.0
nh4_e: -1.0*EX_nh4_e + 1.0*EX_nh4_e_reverse_f9cc6 - 1.0*NH4t + 1.0*NH4t_reverse_551ee = 0.0
o2_e: -1.0*EX_o2_e + 1.0*EX_o2_e_reverse_efa94 - 1.0*O2t + 1.0*O2t_reverse_4d957 = 0.0
oaa_e: -1.0*EX_oaa_e + 1.0*EX_oaa_e_reverse_9b144 + 1.0*OAAt - 1.0*OAAt_reverse_d0cc5 = 0.0
phe-L_e: -1.0*EX_phe-L_e + 1.0*EX_phe-L_e_reverse_155f7 - 1.0*PHEt2r + 1.0*PHEt2r_reverse_19e89 = 0.0
pi_e: -1.0*EX_pi_e + 1.0*EX_pi_e_reverse_1fb09 - 1.0*PIt2r + 1.0*PIt2r_reverse_1cd61 = 0.0
pro-L_e: -1.0*EX_pro-L_e + 1.0*EX_pro-L_e_reverse_26581 - 1.0*PROt2r + 1.0*PROt2r_reverse_f59df = 0.0
pyr_e: -1.0*EX_pyr_e + 1.0*EX_pyr_e_reverse_1f6de + 1.0*PYRt - 1.0*PYRt_reverse_6fee7 = 0.0
ser-L_e: -1.0*EX_ser-L_e + 1.0*EX_ser-L_e_reverse_6b4e0 - 1.0*SERt2r + 1.0*SERt2r_reverse_1b436 = 0.0
so4_e: -1.0*EX_so4_e + 1.0*EX_so4_e_reverse_5c8ed - 1.0*SO4ti + 1.0*SO4ti_reverse_4a93d = 0.0
succ_e: -1.0*EX_succ_e + 1.0*EX_succ_e_reverse_a9039 - 1.0*SUCCt2r + 1.0*SUCCt2r_reverse_e3abe = 0.0
thr-L_e: -1.0*EX_thr-L_e + 1.0*EX_thr-L_e_reverse_9a854 - 1.0*THRt2r + 1.0*THRt2r_reverse_721ff = 0.0
trp-L_e: -1.0*EX_trp-L_e + 1.0*EX_trp-L_e_reverse_273bc - 1.0*TRPt2r + 1.0*TRPt2r_reverse_159c4 = 0.0
tyr-L_e: -1.0*EX_tyr-L_e + 1.0*EX_tyr-L_e_reverse_08f6e - 1.0*TYRt2r + 1.0*TYRt2r_reverse_6f73b = 0.0
val-L_e: -1.0*EX_val-L_e + 1.0*EX_val-L_e_reverse_4c155 - 1.0*VALt2r + 1.0*VALt2r_reverse_e3d24 = 0.0
zymst_e: -1.0*EX_zymst_e + 1.0*EX_zymst_e_reverse_78276 - 1.0*ZYMSTt + 1.0*ZYMSTt_reverse_ab42f = 0.0
6pgl_r: 1.0*6PGLter - 1.0*6PGLter_reverse_0e388 + 1.0*G6PDH2er - 1.0*G6PDH2er_reverse_9f5fd = 0.0
6pgl_c: -1.0*6PGLter + 1.0*6PGLter_reverse_0e388 + 1.0*G6PDH2 - 1.0*G6PDH2_reverse_47bd4 - 1.0*PGL + 1.0*PGL_reverse_2bb6b = 0.0
4abut_c: 1.0*ABUTt2r - 1.0*ABUTt2r_reverse_d925e + 1.0*GLUDC - 1.0*GLUDC_reverse_c4aae = 0.0
acald_c: 1.0*ACALDt - 1.0*ACALDt_reverse_858fa + 1.0*ACALDtm - 1.0*ACALDtm_reverse_c2973 + 1.0*ALCD2if - 1.0*ALCD2if_reverse_85957 - 1.0*ALCD2ir + 1.0*ALCD2ir_reverse_ba067 + 1.0*ALCD2x - 1.0*ALCD2x_reverse_5d107 - 1.0*ALDD2y + 1.0*ALDD2y_reverse_03afb + 1.0*PYRDC - 1.0*PYRDC_reverse_51259 = 0.0
acald_m: -1.0*ACALDtm + 1.0*ACALDtm_reverse_c2973 - 1.0*ALCD2irm + 1.0*ALCD2irm_reverse_2557e - 1.0*ALDD2xm + 1.0*ALDD2xm_reverse_c93a0 - 1.0*ALDD2ym + 1.0*ALDD2ym_reverse_cfc40 = 0.0
h_m: 1.0*ACOAHim - 1.0*ACOAHim_reverse_069ca - 1.0*ACSm + 1.0*ACSm_reverse_5e300 - 1.0*AKGDam + 1.0*AKGDam_reverse_a5f55 - 1.0*ALCD2irm + 1.0*ALCD2irm_reverse_2557e + 2.0*ALDD2xm - 2.0*ALDD2xm_reverse_c93a0 + 2.0*ALDD2ym - 2.0*ALDD2ym_reverse_cfc40 + 1.0*ASPt2m - 1.0*ASPt2m_reverse_c9585 + 2.0*ATPS3m - 2.0*ATPS3m_reverse_5d8d3 + 1.0*ATPtm-H - 1.0*ATPtm-H_reverse_ab831 + 1.0*CSm - 1.0*CSm_reverse_bb3a1 - 6.0*CYOOm + 6.0*CYOOm_reverse_d730e - 1.5*CYOR_u6m + 1.5*CYOR_u6m_reverse_ba2b8 + 1.0*D-LACtm - 1.0*D-LACtm_reverse_08640 + 1.0*GCC2cm - 1.0*GCC2cm_reverse_e896f + 1.0*HCITSm - 1.0*HCITSm_reverse_145c0 - 2.0*LMPD_36_ile-L_c + 2.0*LMPD_36_ile-L_c_reverse_d2bf7 - 2.0*LMPD_37_ile-L_c + 2.0*LMPD_37_ile-L_c_reverse_1ea52 - 2.0*LMPD_38_leu-L_c + 2.0*LMPD_38_leu-L_c_reverse_0e66c - 1.0*LMPD_39_leu-L_c + 1.0*LMPD_39_leu-L_c_reverse_73152 - 100.0*LMPD_48_pc_SC_c + 100.0*LMPD_48_pc_SC_c_reverse_be709 - 100.0*LMPD_49_pc_SC_c + 100.0*LMPD_49_pc_SC_c_reverse_55f95 - 100.0*LMPD_51_pe_SC_c + 100.0*LMPD_51_pe_SC_c_reverse_09777 - 2.0*LMPD_72_val-L_c + 2.0*LMPD_72_val-L_c_reverse_768ff + 1.0*MDHm - 1.0*MDHm_reverse_80605 - 1.0*NADH2-u6m + 1.0*NADH2-u6m_reverse_18109 + 1.0*NADPPPSm - 1.0*NADPPPSm_reverse_64d27 + 1.0*OAAt2m - 1.0*OAAt2m_reverse_42ec8 + 1.0*PDHcm - 1.0*PDHcm_reverse_40ed3 + 1.0*PIt2m - 1.0*PIt2m_reverse_cc401 + 2.0*PPAm - 2.0*PPAm_reverse_6c806 + 1.0*PYRt2m - 1.0*PYRt2m_reverse_fb762 + 1.0*THRt2m - 1.0*THRt2m_reverse_fd5a8 = 0.0
coa_m: 1.0*ACOAHim - 1.0*ACOAHim_reverse_069ca - 1.0*ACSm + 1.0*ACSm_reverse_5e300 - 1.0*AKGDbm + 1.0*AKGDbm_reverse_7cf54 + 1.0*CSm - 1.0*CSm_reverse_bb3a1 + 1.0*HCITSm - 1.0*HCITSm_reverse_145c0 + 1.0*LMPD_39_leu-L_c - 1.0*LMPD_39_leu-L_c_reverse_73152 - 1.0*PDHm + 1.0*PDHm_reverse_31d17 - 1.0*SUCOASm + 1.0*SUCOASm_reverse_87069 = 0.0
accoa_m: -1.0*ACOAHim + 1.0*ACOAHim_reverse_069ca + 1.0*ACSm - 1.0*ACSm_reverse_5e300 - 1.0*CSm + 1.0*CSm_reverse_bb3a1 - 1.0*HCITSm + 1.0*HCITSm_reverse_145c0 - 1.0*LMPD_39_leu-L_c + 1.0*LMPD_39_leu-L_c_reverse_73152 + 1.0*PDHm - 1.0*PDHm_reverse_31d17 = 0.0
ac_m: 1.0*ACOAHim - 1.0*ACOAHim_reverse_069ca - 1.0*ACSm + 1.0*ACSm_reverse_5e300 + 1.0*ACtm - 1.0*ACtm_reverse_5029e + 1.0*ALDD2xm - 1.0*ALDD2xm_reverse_c93a0 + 1.0*ALDD2ym - 1.0*ALDD2ym_reverse_cfc40 = 0.0
h2o_m: -1.0*ACOAHim + 1.0*ACOAHim_reverse_069ca - 1.0*ALDD2xm + 1.0*ALDD2xm_reverse_c93a0 - 1.0*ALDD2ym + 1.0*ALDD2ym_reverse_cfc40 + 1.0*ATPS3m - 1.0*ATPS3m_reverse_5d8d3 - 1.0*CSm + 1.0*CSm_reverse_bb3a1 + 2.0*CYOOm - 2.0*CYOOm_reverse_d730e - 1.0*FUMm + 1.0*FUMm_reverse_5e620 + 1.0*H2Otm - 1.0*H2Otm_reverse_d5f1a - 1.0*HCITSm + 1.0*HCITSm_reverse_145c0 + 1.0*LMPD_36_ile-L_c - 1.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_72_val-L_c - 1.0*LMPD_72_val-L_c_reverse_768ff - 1.0*NADPPPSm + 1.0*NADPPPSm_reverse_64d27 - 1.0*PPAm + 1.0*PPAm_reverse_6c806 = 0.0
cit_c: -1.0*ACONT + 1.0*ACONT_reverse_7c2d5 + 1.0*CITt2r - 1.0*CITt2r_reverse_9f606 - 1.0*CITtam + 1.0*CITtam_reverse_5230c - 1.0*CITtcm + 1.0*CITtcm_reverse_df725 = 0.0
icit_c: 1.0*ACONT - 1.0*ACONT_reverse_7c2d5 + 1.0*CITtcm - 1.0*CITtcm_reverse_df725 - 1.0*ICDHy + 1.0*ICDHy_reverse_c12f8 - 1.0*LMPD_28_gly_c + 1.0*LMPD_28_gly_c_reverse_77428 - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 = 0.0
cit_m: -1.0*ACONTm + 1.0*ACONTm_reverse_8079b + 1.0*CITtam - 1.0*CITtam_reverse_5230c + 1.0*CITtcm - 1.0*CITtcm_reverse_df725 + 1.0*CSm - 1.0*CSm_reverse_bb3a1 = 0.0
icit_m: 1.0*ACONTm - 1.0*ACONTm_reverse_8079b - 1.0*CITtcm + 1.0*CITtcm_reverse_df725 - 1.0*ICDHxm + 1.0*ICDHxm_reverse_14a38 - 1.0*ICDHym + 1.0*ICDHym_reverse_38fb7 = 0.0
ac_c: -1.0*ACS + 1.0*ACS_reverse_37635 + 1.0*ACt2r - 1.0*ACt2r_reverse_64e05 - 1.0*ACtm + 1.0*ACtm_reverse_5029e + 1.0*ACtr - 1.0*ACtr_reverse_b9231 + 1.0*ALDD2y - 1.0*ALDD2y_reverse_03afb + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a + 1.0*LMPD_44_met-L_c - 1.0*LMPD_44_met-L_c_reverse_1646f + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 = 0.0
ppi_c: 1.0*ACS - 1.0*ACS_reverse_37635 + 1.0*LMPD_10_cmp_c - 1.0*LMPD_10_cmp_c_reverse_3acfa + 1.0*LMPD_11_cmp_c - 1.0*LMPD_11_cmp_c_reverse_98fa8 + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a + 1.0*LMPD_13_damp_c - 1.0*LMPD_13_damp_c_reverse_e961c + 1.0*LMPD_14_damp_c - 1.0*LMPD_14_damp_c_reverse_d5737 + 1.0*LMPD_15_damp_c - 1.0*LMPD_15_damp_c_reverse_33cc4 + 1.0*LMPD_16_damp_c - 1.0*LMPD_16_damp_c_reverse_be99f + 1.0*LMPD_17_dcmp_c - 1.0*LMPD_17_dcmp_c_reverse_5c6b0 + 1.0*LMPD_18_dcmp_c - 1.0*LMPD_18_dcmp_c_reverse_15021 + 2.0*LMPD_19_dgmp_c - 2.0*LMPD_19_dgmp_c_reverse_7f2b1 + 1.0*LMPD_1_13BDglcn_c - 1.0*LMPD_1_13BDglcn_c_reverse_e2bdf + 2.0*LMPD_20_dgmp_c - 2.0*LMPD_20_dgmp_c_reverse_e4f3a + 2.0*LMPD_21_dgmp_c - 2.0*LMPD_21_dgmp_c_reverse_65045 + 2.0*LMPD_22_dgmp_c - 2.0*LMPD_22_dgmp_c_reverse_ebc77 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 7.0*LMPD_25_ergst_c - 7.0*LMPD_25_ergst_c_reverse_96537 + 7.0*LMPD_26_ergst_c - 7.0*LMPD_26_ergst_c_reverse_a59e0 + 1.0*LMPD_29_glycogen_c - 1.0*LMPD_29_glycogen_c_reverse_68cb3 + 1.0*LMPD_2_amp_c - 1.0*LMPD_2_amp_c_reverse_bde95 + 2.0*LMPD_30_gmp_c - 2.0*LMPD_30_gmp_c_reverse_1f79a + 2.0*LMPD_31_gmp_c - 2.0*LMPD_31_gmp_c_reverse_d0c02 + 2.0*LMPD_32_his-L_c - 2.0*LMPD_32_his-L_c_reverse_c5967 + 2.0*LMPD_33_his-L_c - 2.0*LMPD_33_his-L_c_reverse_0eb8b + 2.0*LMPD_34_his-L_c - 2.0*LMPD_34_his-L_c_reverse_e16a9 + 2.0*LMPD_35_his-L_c - 2.0*LMPD_35_his-L_c_reverse_5e623 + 1.0*LMPD_3_amp_c - 1.0*LMPD_3_amp_c_reverse_96e30 + 1.0*LMPD_40_lys-L_c - 1.0*LMPD_40_lys-L_c_reverse_4d3a6 + 1.0*LMPD_41_lys-L_c - 1.0*LMPD_41_lys-L_c_reverse_c021a + 1.0*LMPD_42_mannan_c - 1.0*LMPD_42_mannan_c_reverse_9f626 + 1.0*LMPD_43_mannan_c - 1.0*LMPD_43_mannan_c_reverse_7dd82 + 1.0*LMPD_44_met-L_c - 1.0*LMPD_44_met-L_c_reverse_1646f + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 + 400.0*LMPD_48_pc_SC_c - 400.0*LMPD_48_pc_SC_c_reverse_be709 + 400.0*LMPD_49_pc_SC_c - 400.0*LMPD_49_pc_SC_c_reverse_55f95 + 1.0*LMPD_4_amp_c - 1.0*LMPD_4_amp_c_reverse_7925f + 100.0*LMPD_51_pe_SC_c - 100.0*LMPD_51_pe_SC_c_reverse_09777 + 100.0*LMPD_56_ps_SC_c - 100.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 100.0*LMPD_57_ptd1ino_SC_c - 100.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 100.0*LMPD_58_ptd1ino_SC_c - 100.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 + 3.0*LMPD_59_ribflv_c - 3.0*LMPD_59_ribflv_c_reverse_01d9b + 1.0*LMPD_5_amp_c - 1.0*LMPD_5_amp_c_reverse_8e974 + 3.0*LMPD_60_ribflv_c - 3.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_64_tre_c - 1.0*LMPD_64_tre_c_reverse_64f54 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 + 1.0*LMPD_70_ump_c - 1.0*LMPD_70_ump_c_reverse_3dcab + 1.0*LMPD_71_ump_c - 1.0*LMPD_71_ump_c_reverse_dc9d0 + 6.0*LMPD_73_zymst_c - 6.0*LMPD_73_zymst_c_reverse_e6142 + 6.0*LMPD_74_zymst_c - 6.0*LMPD_74_zymst_c_reverse_84617 + 1.0*LMPD_7_asn-L_c - 1.0*LMPD_7_asn-L_c_reverse_b41c7 + 1.0*LMPD_8_cmp_c - 1.0*LMPD_8_cmp_c_reverse_23cbb + 1.0*LMPD_9_cmp_c - 1.0*LMPD_9_cmp_c_reverse_85279 - 1.0*PPA + 1.0*PPA_reverse_c5293 = 0.0
accoa_c: 1.0*ACS - 1.0*ACS_reverse_37635 - 1.0*LMPD_12_cys-L_c + 1.0*LMPD_12_cys-L_c_reverse_6dc0a - 18.0*LMPD_25_ergst_c + 18.0*LMPD_25_ergst_c_reverse_96537 - 18.0*LMPD_26_ergst_c + 18.0*LMPD_26_ergst_c_reverse_a59e0 - 1.0*LMPD_38_leu-L_c + 1.0*LMPD_38_leu-L_c_reverse_0e66c - 1.0*LMPD_44_met-L_c + 1.0*LMPD_44_met-L_c_reverse_1646f - 1.0*LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c_reverse_912f1 - 1620.0*LMPD_47_pa_SC_c + 1620.0*LMPD_47_pa_SC_c_reverse_4049b - 1620.0*LMPD_48_pc_SC_c + 1620.0*LMPD_48_pc_SC_c_reverse_be709 - 1620.0*LMPD_49_pc_SC_c + 1620.0*LMPD_49_pc_SC_c_reverse_55f95 - 1620.0*LMPD_51_pe_SC_c + 1620.0*LMPD_51_pe_SC_c_reverse_09777 - 1620.0*LMPD_56_ps_SC_c + 1620.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 1620.0*LMPD_57_ptd1ino_SC_c + 1620.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 1620.0*LMPD_58_ptd1ino_SC_c + 1620.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 2430.0*LMPD_65_triglyc_SC_c + 2430.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 2430.0*LMPD_66_triglyc_SC_c + 2430.0*LMPD_66_triglyc_SC_c_reverse_42d98 - 18.0*LMPD_73_zymst_c + 18.0*LMPD_73_zymst_c_reverse_e6142 - 18.0*LMPD_74_zymst_c + 18.0*LMPD_74_zymst_c_reverse_84617 = 0.0
coa_c: -1.0*ACS + 1.0*ACS_reverse_37635 + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a + 18.0*LMPD_25_ergst_c - 18.0*LMPD_25_ergst_c_reverse_96537 + 18.0*LMPD_26_ergst_c - 18.0*LMPD_26_ergst_c_reverse_a59e0 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_44_met-L_c - 1.0*LMPD_44_met-L_c_reverse_1646f + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 + 1620.0*LMPD_47_pa_SC_c - 1620.0*LMPD_47_pa_SC_c_reverse_4049b + 1620.0*LMPD_48_pc_SC_c - 1620.0*LMPD_48_pc_SC_c_reverse_be709 + 1620.0*LMPD_49_pc_SC_c - 1620.0*LMPD_49_pc_SC_c_reverse_55f95 + 1620.0*LMPD_51_pe_SC_c - 1620.0*LMPD_51_pe_SC_c_reverse_09777 + 1620.0*LMPD_56_ps_SC_c - 1620.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 1620.0*LMPD_57_ptd1ino_SC_c - 1620.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 1620.0*LMPD_58_ptd1ino_SC_c - 1620.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 + 2430.0*LMPD_65_triglyc_SC_c - 2430.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 2430.0*LMPD_66_triglyc_SC_c - 2430.0*LMPD_66_triglyc_SC_c_reverse_42d98 + 18.0*LMPD_73_zymst_c - 18.0*LMPD_73_zymst_c_reverse_e6142 + 18.0*LMPD_74_zymst_c - 18.0*LMPD_74_zymst_c_reverse_84617 = 0.0
atp_m: -1.0*ACSm + 1.0*ACSm_reverse_5e300 - 1.0*ADK1m + 1.0*ADK1m_reverse_59ef9 + 1.0*ATPS3m - 1.0*ATPS3m_reverse_5d8d3 - 1.0*ATPtm-H + 1.0*ATPtm-H_reverse_ab831 - 1.0*LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c_reverse_33f31 - 1.0*NADKm + 1.0*NADKm_reverse_f0454 - 1.0*SUCOASm + 1.0*SUCOASm_reverse_87069 = 0.0
ppi_m: 1.0*ACSm - 1.0*ACSm_reverse_5e300 - 1.0*PPAm + 1.0*PPAm_reverse_6c806 = 0.0
amp_m: 1.0*ACSm - 1.0*ACSm_reverse_5e300 - 1.0*ADK1m + 1.0*ADK1m_reverse_59ef9 = 0.0
adp_m: 2.0*ADK1m - 2.0*ADK1m_reverse_59ef9 - 1.0*ATPS3m + 1.0*ATPS3m_reverse_5d8d3 + 1.0*ATPtm-H - 1.0*ATPtm-H_reverse_ab831 + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 + 1.0*NADKm - 1.0*NADKm_reverse_f0454 + 1.0*SUCOASm - 1.0*SUCOASm_reverse_87069 = 0.0
co2_m: 1.0*AKGDam - 1.0*AKGDam_reverse_a5f55 + 1.0*CO2tm - 1.0*CO2tm_reverse_48c28 + 1.0*ICDHxm - 1.0*ICDHxm_reverse_14a38 + 1.0*ICDHym - 1.0*ICDHym_reverse_38fb7 + 1.0*LMPD_36_ile-L_c - 1.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_39_leu-L_c - 1.0*LMPD_39_leu-L_c_reverse_73152 + 1.0*LMPD_40_lys-L_c - 1.0*LMPD_40_lys-L_c_reverse_4d3a6 + 1.0*LMPD_41_lys-L_c - 1.0*LMPD_41_lys-L_c_reverse_c021a + 100.0*LMPD_48_pc_SC_c - 100.0*LMPD_48_pc_SC_c_reverse_be709 + 100.0*LMPD_49_pc_SC_c - 100.0*LMPD_49_pc_SC_c_reverse_55f95 + 100.0*LMPD_51_pe_SC_c - 100.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*LMPD_72_val-L_c - 1.0*LMPD_72_val-L_c_reverse_768ff + 1.0*ME1m - 1.0*ME1m_reverse_6868b + 1.0*ME2m - 1.0*ME2m_reverse_f04a4 + 1.0*PDHm - 1.0*PDHm_reverse_31d17 = 0.0
lpam_m: -1.0*AKGDam + 1.0*AKGDam_reverse_a5f55 + 1.0*GCC2cm - 1.0*GCC2cm_reverse_e896f + 1.0*PDHcm - 1.0*PDHcm_reverse_40ed3 = 0.0
sdhlam_m: 1.0*AKGDam - 1.0*AKGDam_reverse_a5f55 - 1.0*AKGDbm + 1.0*AKGDbm_reverse_7cf54 = 0.0
akg_m: -1.0*AKGDam + 1.0*AKGDam_reverse_a5f55 - 1.0*ALATA_Lm + 1.0*ALATA_Lm_reverse_c845b - 1.0*ASPTAm + 1.0*ASPTAm_reverse_dc3b5 - 1.0*HCITSm + 1.0*HCITSm_reverse_145c0 + 1.0*ICDHxm - 1.0*ICDHxm_reverse_14a38 + 1.0*ICDHym - 1.0*ICDHym_reverse_38fb7 + 1.0*LMPD_36_ile-L_c - 1.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 = 0.0
succoa_m: 1.0*AKGDbm - 1.0*AKGDbm_reverse_7cf54 + 1.0*SUCOASm - 1.0*SUCOASm_reverse_87069 = 0.0
dhlam_m: 1.0*AKGDbm - 1.0*AKGDbm_reverse_7cf54 - 1.0*GCC2cm + 1.0*GCC2cm_reverse_e896f - 1.0*PDHcm + 1.0*PDHcm_reverse_40ed3 = 0.0
akg_c: -1.0*AKGMAL + 1.0*AKGMAL_reverse_76d29 + 1.0*AKGt2r - 1.0*AKGt2r_reverse_5d500 - 1.0*ALATA_L + 1.0*ALATA_L_reverse_e54ff - 1.0*ASPTA + 1.0*ASPTA_reverse_36525 + 1.0*GLUDxi - 1.0*GLUDxi_reverse_6fbe2 - 1.0*GLUDyi + 1.0*GLUDyi_reverse_a2351 + 1.0*ICDHy - 1.0*ICDHy_reverse_c12f8 + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a + 2.0*LMPD_13_damp_c - 2.0*LMPD_13_damp_c_reverse_e961c + 2.0*LMPD_14_damp_c - 2.0*LMPD_14_damp_c_reverse_d5737 + 2.0*LMPD_15_damp_c - 2.0*LMPD_15_damp_c_reverse_33cc4 + 2.0*LMPD_16_damp_c - 2.0*LMPD_16_damp_c_reverse_be99f + 2.0*LMPD_19_dgmp_c - 2.0*LMPD_19_dgmp_c_reverse_7f2b1 + 2.0*LMPD_20_dgmp_c - 2.0*LMPD_20_dgmp_c_reverse_e4f3a + 2.0*LMPD_21_dgmp_c - 2.0*LMPD_21_dgmp_c_reverse_65045 + 2.0*LMPD_22_dgmp_c - 2.0*LMPD_22_dgmp_c_reverse_ebc77 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 2.0*LMPD_2_amp_c - 2.0*LMPD_2_amp_c_reverse_bde95 + 2.0*LMPD_30_gmp_c - 2.0*LMPD_30_gmp_c_reverse_1f79a + 2.0*LMPD_31_gmp_c - 2.0*LMPD_31_gmp_c_reverse_d0c02 + 2.0*LMPD_32_his-L_c - 2.0*LMPD_32_his-L_c_reverse_c5967 + 2.0*LMPD_33_his-L_c - 2.0*LMPD_33_his-L_c_reverse_0eb8b + 2.0*LMPD_34_his-L_c - 2.0*LMPD_34_his-L_c_reverse_e16a9 + 2.0*LMPD_35_his-L_c - 2.0*LMPD_35_his-L_c_reverse_5e623 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_39_leu-L_c - 1.0*LMPD_39_leu-L_c_reverse_73152 + 2.0*LMPD_3_amp_c - 2.0*LMPD_3_amp_c_reverse_96e30 + 2.0*LMPD_40_lys-L_c - 2.0*LMPD_40_lys-L_c_reverse_4d3a6 + 2.0*LMPD_41_lys-L_c - 2.0*LMPD_41_lys-L_c_reverse_c021a + 1.0*LMPD_44_met-L_c - 1.0*LMPD_44_met-L_c_reverse_1646f + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 + 400.0*LMPD_48_pc_SC_c - 400.0*LMPD_48_pc_SC_c_reverse_be709 + 400.0*LMPD_49_pc_SC_c - 400.0*LMPD_49_pc_SC_c_reverse_55f95 + 2.0*LMPD_4_amp_c - 2.0*LMPD_4_amp_c_reverse_7925f + 100.0*LMPD_51_pe_SC_c - 100.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*LMPD_52_phe-L_c - 1.0*LMPD_52_phe-L_c_reverse_909a1 + 100.0*LMPD_56_ps_SC_c - 100.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 2.0*LMPD_5_amp_c - 2.0*LMPD_5_amp_c_reverse_8e974 + 1.0*LMPD_61_ser-L_c - 1.0*LMPD_61_ser-L_c_reverse_40666 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*LMPD_68_tyr-L_c - 1.0*LMPD_68_tyr-L_c_reverse_6983e + 1.0*LMPD_69_tyr-L_c - 1.0*LMPD_69_tyr-L_c_reverse_03f73 + 1.0*LMPD_72_val-L_c - 1.0*LMPD_72_val-L_c_reverse_768ff = 0.0
mal-L_c: 1.0*AKGMAL - 1.0*AKGMAL_reverse_76d29 + 1.0*CITtam - 1.0*CITtam_reverse_5230c + 1.0*FUM - 1.0*FUM_reverse_d3642 + 1.0*MALt2r - 1.0*MALt2r_reverse_bfe2f - 1.0*MALtm + 1.0*MALtm_reverse_f798a - 1.0*MDH + 1.0*MDH_reverse_ee52c = 0.0
pyr_c: 1.0*ALATA_L - 1.0*ALATA_L_reverse_e54ff + 1.0*D-LACt2m - 1.0*D-LACt2m_reverse_bb9a8 + 1.0*LMPD_28_gly_c - 1.0*LMPD_28_gly_c_reverse_77428 + 1.0*LMPD_59_ribflv_c - 1.0*LMPD_59_ribflv_c_reverse_01d9b + 1.0*LMPD_60_ribflv_c - 1.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*PC + 1.0*PC_reverse_88dba + 1.0*PYK - 1.0*PYK_reverse_bc8ff - 1.0*PYRDC + 1.0*PYRDC_reverse_51259 - 1.0*PYRt - 1.0*PYRt2m + 1.0*PYRt2m_reverse_fb762 + 1.0*PYRt_reverse_6fee7 = 0.0
ala-L_m: -1.0*ALATA_Lm + 1.0*ALATA_Lm_reverse_c845b - 1.0*ALAtmi + 1.0*ALAtmi_reverse_c89ef = 0.0
glu-L_m: 1.0*ALATA_Lm - 1.0*ALATA_Lm_reverse_c845b + 1.0*ASPGLU2m - 1.0*ASPGLU2m_reverse_efaea + 1.0*ASPTAm - 1.0*ASPTAm_reverse_dc3b5 + 1.0*GLUt7m - 1.0*GLUt7m_reverse_26492 - 1.0*LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 - 2.0*LMPD_6_arg-L_c + 2.0*LMPD_6_arg-L_c_reverse_33f31 = 0.0
pyr_m: 1.0*ALATA_Lm - 1.0*ALATA_Lm_reverse_c845b - 1.0*D-LACt2m + 1.0*D-LACt2m_reverse_bb9a8 - 1.0*LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 - 2.0*LMPD_38_leu-L_c + 2.0*LMPD_38_leu-L_c_reverse_0e66c - 2.0*LMPD_39_leu-L_c + 2.0*LMPD_39_leu-L_c_reverse_73152 - 2.0*LMPD_72_val-L_c + 2.0*LMPD_72_val-L_c_reverse_768ff + 1.0*ME1m - 1.0*ME1m_reverse_6868b + 1.0*ME2m - 1.0*ME2m_reverse_f04a4 - 1.0*PDHm + 1.0*PDHm_reverse_31d17 + 1.0*PYRt2m - 1.0*PYRt2m_reverse_fb762 = 0.0
nadh_c: 1.0*ALCD2if - 1.0*ALCD2if_reverse_85957 - 1.0*ALCD2ir + 1.0*ALCD2ir_reverse_ba067 + 1.0*ALCD2x - 1.0*ALCD2x_reverse_5d107 + 1.0*FDH - 1.0*FDH_reverse_06346 + 1.0*GAPD - 1.0*GAPD_reverse_459c1 + 1.0*GLUDxi - 1.0*GLUDxi_reverse_6fbe2 + 1.0*LMPD_12_cys-L_c - 1.0*LMPD_12_cys-L_c_reverse_6dc0a + 4.0*LMPD_13_damp_c - 4.0*LMPD_13_damp_c_reverse_e961c + 4.0*LMPD_14_damp_c - 4.0*LMPD_14_damp_c_reverse_d5737 + 2.0*LMPD_15_damp_c - 2.0*LMPD_15_damp_c_reverse_33cc4 + 2.0*LMPD_16_damp_c - 2.0*LMPD_16_damp_c_reverse_be99f + 5.0*LMPD_19_dgmp_c - 5.0*LMPD_19_dgmp_c_reverse_7f2b1 + 5.0*LMPD_20_dgmp_c - 5.0*LMPD_20_dgmp_c_reverse_e4f3a + 3.0*LMPD_21_dgmp_c - 3.0*LMPD_21_dgmp_c_reverse_65045 + 3.0*LMPD_22_dgmp_c - 3.0*LMPD_22_dgmp_c_reverse_ebc77 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 2.0*LMPD_25_ergst_c - 2.0*LMPD_25_ergst_c_reverse_96537 - 1.0*LMPD_26_ergst_c + 1.0*LMPD_26_ergst_c_reverse_a59e0 + 2.0*LMPD_2_amp_c - 2.0*LMPD_2_amp_c_reverse_bde95 + 5.0*LMPD_30_gmp_c - 5.0*LMPD_30_gmp_c_reverse_1f79a + 3.0*LMPD_31_gmp_c - 3.0*LMPD_31_gmp_c_reverse_d0c02 + 3.0*LMPD_32_his-L_c - 3.0*LMPD_32_his-L_c_reverse_c5967 + 4.0*LMPD_33_his-L_c - 4.0*LMPD_33_his-L_c_reverse_0eb8b + 3.0*LMPD_34_his-L_c - 3.0*LMPD_34_his-L_c_reverse_e16a9 + 4.0*LMPD_35_his-L_c - 4.0*LMPD_35_his-L_c_reverse_5e623 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_39_leu-L_c - 1.0*LMPD_39_leu-L_c_reverse_73152 + 4.0*LMPD_3_amp_c - 4.0*LMPD_3_amp_c_reverse_96e30 + 1.0*LMPD_40_lys-L_c - 1.0*LMPD_40_lys-L_c_reverse_4d3a6 + 1.0*LMPD_45_met-L_c - 1.0*LMPD_45_met-L_c_reverse_912f1 - 100.0*LMPD_47_pa_SC_c + 100.0*LMPD_47_pa_SC_c_reverse_4049b + 400.0*LMPD_48_pc_SC_c - 400.0*LMPD_48_pc_SC_c_reverse_be709 + 300.0*LMPD_49_pc_SC_c - 300.0*LMPD_49_pc_SC_c_reverse_55f95 + 2.0*LMPD_4_amp_c - 2.0*LMPD_4_amp_c_reverse_7925f - 1.0*LMPD_53_pro-L_c + 1.0*LMPD_53_pro-L_c_reverse_830d7 - 100.0*LMPD_58_ptd1ino_SC_c + 100.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 + 1.0*LMPD_59_ribflv_c - 1.0*LMPD_59_ribflv_c_reverse_01d9b + 4.0*LMPD_5_amp_c - 4.0*LMPD_5_amp_c_reverse_8e974 + 1.0*LMPD_60_ribflv_c - 1.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_61_ser-L_c - 1.0*LMPD_61_ser-L_c_reverse_40666 - 1.0*LMPD_62_thr-L_c + 1.0*LMPD_62_thr-L_c_reverse_e4d09 - 100.0*LMPD_65_triglyc_SC_c + 100.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*LMPD_68_tyr-L_c - 1.0*LMPD_68_tyr-L_c_reverse_6983e - 1.0*LMPD_73_zymst_c + 1.0*LMPD_73_zymst_c_reverse_e6142 + 2.0*LMPD_74_zymst_c - 2.0*LMPD_74_zymst_c_reverse_84617 + 1.0*MDH - 1.0*MDH_reverse_ee52c - 1.0*NADH2-u6cm + 1.0*NADH2-u6cm_reverse_07733 = 0.0
nad_c: -1.0*ALCD2if + 1.0*ALCD2if_reverse_85957 + 1.0*ALCD2ir - 1.0*ALCD2ir_reverse_ba067 - 1.0*ALCD2x + 1.0*ALCD2x_reverse_5d107 - 1.0*FDH + 1.0*FDH_reverse_06346 - 1.0*GAPD + 1.0*GAPD_reverse_459c1 - 1.0*GLUDxi + 1.0*GLUDxi_reverse_6fbe2 - 1.0*LMPD_12_cys-L_c + 1.0*LMPD_12_cys-L_c_reverse_6dc0a - 4.0*LMPD_13_damp_c + 4.0*LMPD_13_damp_c_reverse_e961c - 4.0*LMPD_14_damp_c + 4.0*LMPD_14_damp_c_reverse_d5737 - 2.0*LMPD_15_damp_c + 2.0*LMPD_15_damp_c_reverse_33cc4 - 2.0*LMPD_16_damp_c + 2.0*LMPD_16_damp_c_reverse_be99f - 5.0*LMPD_19_dgmp_c + 5.0*LMPD_19_dgmp_c_reverse_7f2b1 - 5.0*LMPD_20_dgmp_c + 5.0*LMPD_20_dgmp_c_reverse_e4f3a - 3.0*LMPD_21_dgmp_c + 3.0*LMPD_21_dgmp_c_reverse_65045 - 3.0*LMPD_22_dgmp_c + 3.0*LMPD_22_dgmp_c_reverse_ebc77 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 - 2.0*LMPD_25_ergst_c + 2.0*LMPD_25_ergst_c_reverse_96537 + 1.0*LMPD_26_ergst_c - 1.0*LMPD_26_ergst_c_reverse_a59e0 - 2.0*LMPD_2_amp_c + 2.0*LMPD_2_amp_c_reverse_bde95 - 5.0*LMPD_30_gmp_c + 5.0*LMPD_30_gmp_c_reverse_1f79a - 3.0*LMPD_31_gmp_c + 3.0*LMPD_31_gmp_c_reverse_d0c02 - 3.0*LMPD_32_his-L_c + 3.0*LMPD_32_his-L_c_reverse_c5967 - 4.0*LMPD_33_his-L_c + 4.0*LMPD_33_his-L_c_reverse_0eb8b - 3.0*LMPD_34_his-L_c + 3.0*LMPD_34_his-L_c_reverse_e16a9 - 4.0*LMPD_35_his-L_c + 4.0*LMPD_35_his-L_c_reverse_5e623 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 - 1.0*LMPD_38_leu-L_c + 1.0*LMPD_38_leu-L_c_reverse_0e66c - 1.0*LMPD_39_leu-L_c + 1.0*LMPD_39_leu-L_c_reverse_73152 - 4.0*LMPD_3_amp_c + 4.0*LMPD_3_amp_c_reverse_96e30 - 1.0*LMPD_40_lys-L_c + 1.0*LMPD_40_lys-L_c_reverse_4d3a6 - 1.0*LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c_reverse_912f1 + 100.0*LMPD_47_pa_SC_c - 100.0*LMPD_47_pa_SC_c_reverse_4049b - 400.0*LMPD_48_pc_SC_c + 400.0*LMPD_48_pc_SC_c_reverse_be709 - 300.0*LMPD_49_pc_SC_c + 300.0*LMPD_49_pc_SC_c_reverse_55f95 - 2.0*LMPD_4_amp_c + 2.0*LMPD_4_amp_c_reverse_7925f + 1.0*LMPD_53_pro-L_c - 1.0*LMPD_53_pro-L_c_reverse_830d7 + 100.0*LMPD_58_ptd1ino_SC_c - 100.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 4.0*LMPD_5_amp_c + 4.0*LMPD_5_amp_c_reverse_8e974 - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 - 1.0*LMPD_61_ser-L_c + 1.0*LMPD_61_ser-L_c_reverse_40666 + 1.0*LMPD_62_thr-L_c - 1.0*LMPD_62_thr-L_c_reverse_e4d09 + 100.0*LMPD_65_triglyc_SC_c - 100.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 1.0*LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_68_tyr-L_c + 1.0*LMPD_68_tyr-L_c_reverse_6983e + 1.0*LMPD_73_zymst_c - 1.0*LMPD_73_zymst_c_reverse_e6142 - 2.0*LMPD_74_zymst_c + 2.0*LMPD_74_zymst_c_reverse_84617 - 1.0*MDH + 1.0*MDH_reverse_ee52c + 1.0*NADH2-u6cm - 1.0*NADH2-u6cm_reverse_07733 - 1.0*NADK + 1.0*NADK_reverse_bba52 + 1.0*NADPPPS - 1.0*NADPPPS_reverse_c9c40 = 0.0
etoh_c: -1.0*ALCD2if + 1.0*ALCD2if_reverse_85957 + 1.0*ALCD2ir - 1.0*ALCD2ir_reverse_ba067 - 1.0*ALCD2x + 1.0*ALCD2x_reverse_5d107 + 1.0*ETOHt - 1.0*ETOHt_reverse_a82a5 - 1.0*ETOHtm + 1.0*ETOHtm_reverse_9d4eb = 0.0
nadh_m: -1.0*ALCD2irm + 1.0*ALCD2irm_reverse_2557e + 1.0*ALDD2xm - 1.0*ALDD2xm_reverse_c93a0 + 1.0*GCC2cm - 1.0*GCC2cm_reverse_e896f + 1.0*ICDHxm - 1.0*ICDHxm_reverse_14a38 + 1.0*LMPD_40_lys-L_c - 1.0*LMPD_40_lys-L_c_reverse_4d3a6 + 1.0*LMPD_41_lys-L_c - 1.0*LMPD_41_lys-L_c_reverse_c021a + 1.0*MDHm - 1.0*MDHm_reverse_80605 + 1.0*ME1m - 1.0*ME1m_reverse_6868b - 1.0*NADH2-u6m + 1.0*NADH2-u6m_reverse_18109 + 1.0*PDHcm - 1.0*PDHcm_reverse_40ed3 + 1.0*PDHm - 1.0*PDHm_reverse_31d17 = 0.0
nad_m: 1.0*ALCD2irm - 1.0*ALCD2irm_reverse_2557e - 1.0*ALDD2xm + 1.0*ALDD2xm_reverse_c93a0 - 1.0*GCC2cm + 1.0*GCC2cm_reverse_e896f - 1.0*ICDHxm + 1.0*ICDHxm_reverse_14a38 - 1.0*LMPD_40_lys-L_c + 1.0*LMPD_40_lys-L_c_reverse_4d3a6 - 1.0*LMPD_41_lys-L_c + 1.0*LMPD_41_lys-L_c_reverse_c021a - 1.0*MDHm + 1.0*MDHm_reverse_80605 - 1.0*ME1m + 1.0*ME1m_reverse_6868b + 1.0*NADH2-u6m - 1.0*NADH2-u6m_reverse_18109 - 1.0*NADKm + 1.0*NADKm_reverse_f0454 + 1.0*NADPPPSm - 1.0*NADPPPSm_reverse_64d27 - 1.0*PDHcm + 1.0*PDHcm_reverse_40ed3 - 1.0*PDHm + 1.0*PDHm_reverse_31d17 = 0.0
etoh_m: 1.0*ALCD2irm - 1.0*ALCD2irm_reverse_2557e + 1.0*ETOHtm - 1.0*ETOHtm_reverse_9d4eb = 0.0
nadph_c: 1.0*ALDD2y - 1.0*ALDD2y_reverse_03afb + 1.0*G6PDH2 - 1.0*G6PDH2_reverse_47bd4 - 1.0*GLUDyi + 1.0*GLUDyi_reverse_a2351 + 1.0*GND - 1.0*GND_reverse_eec5c + 1.0*ICDHy - 1.0*ICDHy_reverse_c12f8 - 4.0*LMPD_12_cys-L_c + 4.0*LMPD_12_cys-L_c_reverse_6dc0a - 1.0*LMPD_13_damp_c + 1.0*LMPD_13_damp_c_reverse_e961c - 1.0*LMPD_14_damp_c + 1.0*LMPD_14_damp_c_reverse_d5737 + 1.0*LMPD_15_damp_c - 1.0*LMPD_15_damp_c_reverse_33cc4 + 1.0*LMPD_16_damp_c - 1.0*LMPD_16_damp_c_reverse_be99f - 1.0*LMPD_17_dcmp_c + 1.0*LMPD_17_dcmp_c_reverse_5c6b0 - 1.0*LMPD_18_dcmp_c + 1.0*LMPD_18_dcmp_c_reverse_15021 - 1.0*LMPD_19_dgmp_c + 1.0*LMPD_19_dgmp_c_reverse_7f2b1 - 1.0*LMPD_20_dgmp_c + 1.0*LMPD_20_dgmp_c_reverse_e4f3a + 1.0*LMPD_21_dgmp_c - 1.0*LMPD_21_dgmp_c_reverse_65045 + 1.0*LMPD_22_dgmp_c - 1.0*LMPD_22_dgmp_c_reverse_ebc77 - 2.0*LMPD_23_dtmp_c + 2.0*LMPD_23_dtmp_c_reverse_ed12d - 2.0*LMPD_24_dtmp_c + 2.0*LMPD_24_dtmp_c_reverse_5d241 - 27.0*LMPD_25_ergst_c + 27.0*LMPD_25_ergst_c_reverse_96537 - 24.0*LMPD_26_ergst_c + 24.0*LMPD_26_ergst_c_reverse_a59e0 + 2.0*LMPD_2_amp_c - 2.0*LMPD_2_amp_c_reverse_bde95 + 2.0*LMPD_31_gmp_c - 2.0*LMPD_31_gmp_c_reverse_d0c02 + 1.0*LMPD_32_his-L_c - 1.0*LMPD_32_his-L_c_reverse_c5967 + 1.0*LMPD_34_his-L_c - 1.0*LMPD_34_his-L_c_reverse_e16a9 - 2.0*LMPD_36_ile-L_c + 2.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 - 2.0*LMPD_40_lys-L_c + 2.0*LMPD_40_lys-L_c_reverse_4d3a6 - 1.0*LMPD_41_lys-L_c + 1.0*LMPD_41_lys-L_c_reverse_c021a - 6.0*LMPD_44_met-L_c + 6.0*LMPD_44_met-L_c_reverse_1646f - 7.0*LMPD_45_met-L_c + 7.0*LMPD_45_met-L_c_reverse_912f1 - 2958.0*LMPD_47_pa_SC_c + 2958.0*LMPD_47_pa_SC_c_reverse_4049b - 3358.0*LMPD_48_pc_SC_c + 3358.0*LMPD_48_pc_SC_c_reverse_be709 - 3258.0*LMPD_49_pc_SC_c + 3258.0*LMPD_49_pc_SC_c_reverse_55f95 + 2.0*LMPD_4_amp_c - 2.0*LMPD_4_amp_c_reverse_7925f - 2958.0*LMPD_51_pe_SC_c + 2958.0*LMPD_51_pe_SC_c_reverse_09777 - 1.0*LMPD_52_phe-L_c + 1.0*LMPD_52_phe-L_c_reverse_909a1 - 1.0*LMPD_53_pro-L_c + 1.0*LMPD_53_pro-L_c_reverse_830d7 - 2.0*LMPD_54_pro-L_c + 2.0*LMPD_54_pro-L_c_reverse_a6168 - 2958.0*LMPD_56_ps_SC_c + 2958.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 3058.0*LMPD_57_ptd1ino_SC_c + 3058.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 2958.0*LMPD_58_ptd1ino_SC_c + 2958.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 - 1.0*LMPD_62_thr-L_c + 1.0*LMPD_62_thr-L_c_reverse_e4d09 - 2.0*LMPD_63_thr-L_c + 2.0*LMPD_63_thr-L_c_reverse_52245 - 4437.0*LMPD_65_triglyc_SC_c + 4437.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 4537.0*LMPD_66_triglyc_SC_c + 4537.0*LMPD_66_triglyc_SC_c_reverse_42d98 - 1.0*LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_68_tyr-L_c + 1.0*LMPD_68_tyr-L_c_reverse_6983e - 22.0*LMPD_73_zymst_c + 22.0*LMPD_73_zymst_c_reverse_e6142 - 25.0*LMPD_74_zymst_c + 25.0*LMPD_74_zymst_c_reverse_84617 = 0.0
nadp_c: -1.0*ALDD2y + 1.0*ALDD2y_reverse_03afb - 1.0*G6PDH2 + 1.0*G6PDH2_reverse_47bd4 + 1.0*GLUDyi - 1.0*GLUDyi_reverse_a2351 - 1.0*GND + 1.0*GND_reverse_eec5c - 1.0*ICDHy + 1.0*ICDHy_reverse_c12f8 + 4.0*LMPD_12_cys-L_c - 4.0*LMPD_12_cys-L_c_reverse_6dc0a + 1.0*LMPD_13_damp_c - 1.0*LMPD_13_damp_c_reverse_e961c + 1.0*LMPD_14_damp_c - 1.0*LMPD_14_damp_c_reverse_d5737 - 1.0*LMPD_15_damp_c + 1.0*LMPD_15_damp_c_reverse_33cc4 - 1.0*LMPD_16_damp_c + 1.0*LMPD_16_damp_c_reverse_be99f + 1.0*LMPD_17_dcmp_c - 1.0*LMPD_17_dcmp_c_reverse_5c6b0 + 1.0*LMPD_18_dcmp_c - 1.0*LMPD_18_dcmp_c_reverse_15021 + 1.0*LMPD_19_dgmp_c - 1.0*LMPD_19_dgmp_c_reverse_7f2b1 + 1.0*LMPD_20_dgmp_c - 1.0*LMPD_20_dgmp_c_reverse_e4f3a - 1.0*LMPD_21_dgmp_c + 1.0*LMPD_21_dgmp_c_reverse_65045 - 1.0*LMPD_22_dgmp_c + 1.0*LMPD_22_dgmp_c_reverse_ebc77 + 2.0*LMPD_23_dtmp_c - 2.0*LMPD_23_dtmp_c_reverse_ed12d + 2.0*LMPD_24_dtmp_c - 2.0*LMPD_24_dtmp_c_reverse_5d241 + 27.0*LMPD_25_ergst_c - 27.0*LMPD_25_ergst_c_reverse_96537 + 24.0*LMPD_26_ergst_c - 24.0*LMPD_26_ergst_c_reverse_a59e0 - 2.0*LMPD_2_amp_c + 2.0*LMPD_2_amp_c_reverse_bde95 - 2.0*LMPD_31_gmp_c + 2.0*LMPD_31_gmp_c_reverse_d0c02 - 1.0*LMPD_32_his-L_c + 1.0*LMPD_32_his-L_c_reverse_c5967 - 1.0*LMPD_34_his-L_c + 1.0*LMPD_34_his-L_c_reverse_e16a9 + 2.0*LMPD_36_ile-L_c - 2.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 2.0*LMPD_40_lys-L_c - 2.0*LMPD_40_lys-L_c_reverse_4d3a6 + 1.0*LMPD_41_lys-L_c - 1.0*LMPD_41_lys-L_c_reverse_c021a + 6.0*LMPD_44_met-L_c - 6.0*LMPD_44_met-L_c_reverse_1646f + 7.0*LMPD_45_met-L_c - 7.0*LMPD_45_met-L_c_reverse_912f1 + 2958.0*LMPD_47_pa_SC_c - 2958.0*LMPD_47_pa_SC_c_reverse_4049b + 3358.0*LMPD_48_pc_SC_c - 3358.0*LMPD_48_pc_SC_c_reverse_be709 + 3258.0*LMPD_49_pc_SC_c - 3258.0*LMPD_49_pc_SC_c_reverse_55f95 - 2.0*LMPD_4_amp_c + 2.0*LMPD_4_amp_c_reverse_7925f + 2958.0*LMPD_51_pe_SC_c - 2958.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*LMPD_52_phe-L_c - 1.0*LMPD_52_phe-L_c_reverse_909a1 + 1.0*LMPD_53_pro-L_c - 1.0*LMPD_53_pro-L_c_reverse_830d7 + 2.0*LMPD_54_pro-L_c - 2.0*LMPD_54_pro-L_c_reverse_a6168 + 2958.0*LMPD_56_ps_SC_c - 2958.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 3058.0*LMPD_57_ptd1ino_SC_c - 3058.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 2958.0*LMPD_58_ptd1ino_SC_c - 2958.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 + 1.0*LMPD_59_ribflv_c - 1.0*LMPD_59_ribflv_c_reverse_01d9b + 1.0*LMPD_60_ribflv_c - 1.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_62_thr-L_c - 1.0*LMPD_62_thr-L_c_reverse_e4d09 + 2.0*LMPD_63_thr-L_c - 2.0*LMPD_63_thr-L_c_reverse_52245 + 4437.0*LMPD_65_triglyc_SC_c - 4437.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 4537.0*LMPD_66_triglyc_SC_c - 4537.0*LMPD_66_triglyc_SC_c_reverse_42d98 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*LMPD_68_tyr-L_c - 1.0*LMPD_68_tyr-L_c_reverse_6983e + 22.0*LMPD_73_zymst_c - 22.0*LMPD_73_zymst_c_reverse_e6142 + 25.0*LMPD_74_zymst_c - 25.0*LMPD_74_zymst_c_reverse_84617 + 1.0*NADK - 1.0*NADK_reverse_bba52 - 1.0*NADPPPS + 1.0*NADPPPS_reverse_c9c40 = 0.0
nadp_m: -1.0*ALDD2ym + 1.0*ALDD2ym_reverse_cfc40 - 1.0*ICDHym + 1.0*ICDHym_reverse_38fb7 + 1.0*LMPD_36_ile-L_c - 1.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_39_leu-L_c - 1.0*LMPD_39_leu-L_c_reverse_73152 + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 + 1.0*LMPD_72_val-L_c - 1.0*LMPD_72_val-L_c_reverse_768ff - 1.0*ME2m + 1.0*ME2m_reverse_f04a4 + 1.0*NADKm - 1.0*NADKm_reverse_f0454 - 1.0*NADPPPSm + 1.0*NADPPPSm_reverse_64d27 = 0.0
nadph_m: 1.0*ALDD2ym - 1.0*ALDD2ym_reverse_cfc40 + 1.0*ICDHym - 1.0*ICDHym_reverse_38fb7 - 1.0*LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 - 1.0*LMPD_38_leu-L_c + 1.0*LMPD_38_leu-L_c_reverse_0e66c - 1.0*LMPD_39_leu-L_c + 1.0*LMPD_39_leu-L_c_reverse_73152 - 1.0*LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c_reverse_33f31 - 1.0*LMPD_72_val-L_c + 1.0*LMPD_72_val-L_c_reverse_768ff + 1.0*ME2m - 1.0*ME2m_reverse_f04a4 = 0.0
asp-L_m: -1.0*ASPGLU2m + 1.0*ASPGLU2m_reverse_efaea - 1.0*ASPTAm + 1.0*ASPTAm_reverse_dc3b5 + 1.0*ASPt2m - 1.0*ASPt2m_reverse_c9585 = 0.0
oaa_c: 1.0*ASPTA - 1.0*ASPTA_reverse_36525 + 1.0*MDH - 1.0*MDH_reverse_ee52c - 1.0*OAAt - 1.0*OAAt2m + 1.0*OAAt2m_reverse_42ec8 + 1.0*OAAt_reverse_d0cc5 + 1.0*PC - 1.0*PC_reverse_88dba - 1.0*PPCK + 1.0*PPCK_reverse_2557d = 0.0
oaa_m: 1.0*ASPTAm - 1.0*ASPTAm_reverse_dc3b5 - 1.0*CSm + 1.0*CSm_reverse_bb3a1 + 1.0*MDHm - 1.0*MDHm_reverse_80605 + 1.0*OAAt2m - 1.0*OAAt2m_reverse_42ec8 = 0.0
pi_m: -1.0*ATPS3m + 1.0*ATPS3m_reverse_5d8d3 + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 - 1.0*MALtm + 1.0*MALtm_reverse_f798a + 1.0*NADPPPSm - 1.0*NADPPPSm_reverse_64d27 + 1.0*PIt2m - 1.0*PIt2m_reverse_cc401 + 2.0*PPAm - 2.0*PPAm_reverse_6c806 - 1.0*SUCCtm + 1.0*SUCCtm_reverse_5f36b + 1.0*SUCOASm - 1.0*SUCOASm_reverse_87069 = 0.0
mal-L_m: -1.0*CITtam + 1.0*CITtam_reverse_5230c + 1.0*FUMm - 1.0*FUMm_reverse_5e620 + 1.0*MALtm - 1.0*MALtm_reverse_f798a - 1.0*MDHm + 1.0*MDHm_reverse_80605 - 1.0*ME1m + 1.0*ME1m_reverse_6868b - 1.0*ME2m + 1.0*ME2m_reverse_f04a4 = 0.0
co2_c: 1.0*CO2t - 1.0*CO2t_reverse_7c42f - 1.0*CO2tm + 1.0*CO2tm_reverse_48c28 + 1.0*FDH - 1.0*FDH_reverse_06346 + 1.0*FDNG - 1.0*FDNG_reverse_f3615 + 1.0*GLUDC - 1.0*GLUDC_reverse_c4aae + 1.0*GND - 1.0*GND_reverse_eec5c - 1.0*HCO3E + 1.0*HCO3E_reverse_97ea5 + 1.0*ICDHy - 1.0*ICDHy_reverse_c12f8 + 1.0*LMPD_10_cmp_c - 1.0*LMPD_10_cmp_c_reverse_3acfa + 1.0*LMPD_11_cmp_c - 1.0*LMPD_11_cmp_c_reverse_98fa8 - 1.0*LMPD_13_damp_c + 1.0*LMPD_13_damp_c_reverse_e961c - 1.0*LMPD_14_damp_c + 1.0*LMPD_14_damp_c_reverse_d5737 - 1.0*LMPD_15_damp_c + 1.0*LMPD_15_damp_c_reverse_33cc4 - 1.0*LMPD_16_damp_c + 1.0*LMPD_16_damp_c_reverse_be99f + 1.0*LMPD_17_dcmp_c - 1.0*LMPD_17_dcmp_c_reverse_5c6b0 + 1.0*LMPD_18_dcmp_c - 1.0*LMPD_18_dcmp_c_reverse_15021 - 1.0*LMPD_19_dgmp_c + 1.0*LMPD_19_dgmp_c_reverse_7f2b1 - 1.0*LMPD_20_dgmp_c + 1.0*LMPD_20_dgmp_c_reverse_e4f3a - 1.0*LMPD_21_dgmp_c + 1.0*LMPD_21_dgmp_c_reverse_65045 - 1.0*LMPD_22_dgmp_c + 1.0*LMPD_22_dgmp_c_reverse_ebc77 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 8.0*LMPD_25_ergst_c - 8.0*LMPD_25_ergst_c_reverse_96537 + 8.0*LMPD_26_ergst_c - 8.0*LMPD_26_ergst_c_reverse_a59e0 - 1.0*LMPD_2_amp_c + 1.0*LMPD_2_amp_c_reverse_bde95 - 1.0*LMPD_30_gmp_c + 1.0*LMPD_30_gmp_c_reverse_1f79a - 1.0*LMPD_31_gmp_c + 1.0*LMPD_31_gmp_c_reverse_d0c02 + 1.0*LMPD_38_leu-L_c - 1.0*LMPD_38_leu-L_c_reverse_0e66c + 1.0*LMPD_39_leu-L_c - 1.0*LMPD_39_leu-L_c_reverse_73152 - 1.0*LMPD_3_amp_c + 1.0*LMPD_3_amp_c_reverse_96e30 + 1420.0*LMPD_47_pa_SC_c - 1420.0*LMPD_47_pa_SC_c_reverse_4049b + 1420.0*LMPD_48_pc_SC_c - 1420.0*LMPD_48_pc_SC_c_reverse_be709 + 1420.0*LMPD_49_pc_SC_c - 1420.0*LMPD_49_pc_SC_c_reverse_55f95 - 1.0*LMPD_4_amp_c + 1.0*LMPD_4_amp_c_reverse_7925f + 1420.0*LMPD_51_pe_SC_c - 1420.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*LMPD_52_phe-L_c - 1.0*LMPD_52_phe-L_c_reverse_909a1 + 1420.0*LMPD_56_ps_SC_c - 1420.0*LMPD_56_ps_SC_c_reverse_7e9d7 + 1420.0*LMPD_57_ptd1ino_SC_c - 1420.0*LMPD_57_ptd1ino_SC_c_reverse_65089 + 1420.0*LMPD_58_ptd1ino_SC_c - 1420.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 1.0*LMPD_5_amp_c + 1.0*LMPD_5_amp_c_reverse_8e974 - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 + 2130.0*LMPD_65_triglyc_SC_c - 2130.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 + 2130.0*LMPD_66_triglyc_SC_c - 2130.0*LMPD_66_triglyc_SC_c_reverse_42d98 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 + 1.0*LMPD_68_tyr-L_c - 1.0*LMPD_68_tyr-L_c_reverse_6983e + 1.0*LMPD_69_tyr-L_c - 1.0*LMPD_69_tyr-L_c_reverse_03f73 + 1.0*LMPD_70_ump_c - 1.0*LMPD_70_ump_c_reverse_3dcab + 1.0*LMPD_71_ump_c - 1.0*LMPD_71_ump_c_reverse_dc9d0 + 8.0*LMPD_73_zymst_c - 8.0*LMPD_73_zymst_c_reverse_e6142 + 8.0*LMPD_74_zymst_c - 8.0*LMPD_74_zymst_c_reverse_84617 + 1.0*LMPD_8_cmp_c - 1.0*LMPD_8_cmp_c_reverse_23cbb + 1.0*LMPD_9_cmp_c - 1.0*LMPD_9_cmp_c_reverse_85279 + 1.0*PPCK - 1.0*PPCK_reverse_2557d + 1.0*PYRDC - 1.0*PYRDC_reverse_51259 = 0.0
ficytc_m: 4.0*CYOOm - 4.0*CYOOm_reverse_d730e - 2.0*CYOR_u6m + 2.0*CYOR_u6m_reverse_ba2b8 = 0.0
focytc_m: -4.0*CYOOm + 4.0*CYOOm_reverse_d730e + 2.0*CYOR_u6m - 2.0*CYOR_u6m_reverse_ba2b8 = 0.0
o2_m: -1.0*CYOOm + 1.0*CYOOm_reverse_d730e + 1.0*O2tm - 1.0*O2tm_reverse_fe4d6 = 0.0
q6h2_m: -1.0*CYOR_u6m + 1.0*CYOR_u6m_reverse_ba2b8 + 1.0*FDNG - 1.0*FDNG_reverse_f3615 + 1.0*LMPD_10_cmp_c - 1.0*LMPD_10_cmp_c_reverse_3acfa + 1.0*LMPD_18_dcmp_c - 1.0*LMPD_18_dcmp_c_reverse_15021 + 1.0*LMPD_23_dtmp_c - 1.0*LMPD_23_dtmp_c_reverse_ed12d + 1.0*LMPD_71_ump_c - 1.0*LMPD_71_ump_c_reverse_dc9d0 + 1.0*LMPD_8_cmp_c - 1.0*LMPD_8_cmp_c_reverse_23cbb + 1.0*NADH2-u6cm - 1.0*NADH2-u6cm_reverse_07733 + 1.0*NADH2-u6m - 1.0*NADH2-u6m_reverse_18109 + 1.0*SUCD2_u6m - 1.0*SUCD2_u6m_reverse_b6f79 + 1.0*SUCD3-u6m - 1.0*SUCD3-u6m_reverse_0e152 = 0.0
q6_m: 1.0*CYOR_u6m - 1.0*CYOR_u6m_reverse_ba2b8 - 1.0*FDNG + 1.0*FDNG_reverse_f3615 - 1.0*LMPD_10_cmp_c + 1.0*LMPD_10_cmp_c_reverse_3acfa - 1.0*LMPD_18_dcmp_c + 1.0*LMPD_18_dcmp_c_reverse_15021 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_71_ump_c + 1.0*LMPD_71_ump_c_reverse_dc9d0 - 1.0*LMPD_8_cmp_c + 1.0*LMPD_8_cmp_c_reverse_23cbb - 1.0*NADH2-u6cm + 1.0*NADH2-u6cm_reverse_07733 - 1.0*NADH2-u6m + 1.0*NADH2-u6m_reverse_18109 - 1.0*SUCD2_u6m + 1.0*SUCD2_u6m_reverse_b6f79 - 1.0*SUCD3-u6m + 1.0*SUCD3-u6m_reverse_0e152 = 0.0
lac-D_c: -1.0*D-LACt2m + 1.0*D-LACt2m_reverse_bb9a8 - 1.0*D-LACtm + 1.0*D-LACtm_reverse_08640 = 0.0
lac-D_m: 1.0*D-LACt2m - 1.0*D-LACt2m_reverse_bb9a8 + 1.0*D-LACtm - 1.0*D-LACtm_reverse_08640 = 0.0
db4p_c: 1.0*DB4PS - 1.0*DB4PS_reverse_43dd1 - 2.0*LMPD_59_ribflv_c + 2.0*LMPD_59_ribflv_c_reverse_01d9b - 2.0*LMPD_60_ribflv_c + 2.0*LMPD_60_ribflv_c_reverse_b7a89 = 0.0
for_c: 1.0*DB4PS - 1.0*DB4PS_reverse_43dd1 - 1.0*FDH + 1.0*FDH_reverse_06346 - 1.0*FDNG + 1.0*FDNG_reverse_f3615 + 1.0*FORt - 1.0*FORt_reverse_40f9f - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_73_zymst_c - 1.0*LMPD_73_zymst_c_reverse_e6142 + 1.0*LMPD_74_zymst_c - 1.0*LMPD_74_zymst_c_reverse_84617 = 0.0
ru5p-D_c: -1.0*DB4PS + 1.0*DB4PS_reverse_43dd1 + 1.0*GND - 1.0*GND_reverse_eec5c - 1.0*RPE + 1.0*RPE_reverse_a1b04 + 1.0*RPI - 1.0*RPI_reverse_853a1 = 0.0
pep_c: 1.0*ENO - 1.0*ENO_reverse_40eea - 2.0*LMPD_52_phe-L_c + 2.0*LMPD_52_phe-L_c_reverse_909a1 - 2.0*LMPD_67_trp-L_c + 2.0*LMPD_67_trp-L_c_reverse_17a73 - 2.0*LMPD_68_tyr-L_c + 2.0*LMPD_68_tyr-L_c_reverse_6983e - 2.0*LMPD_69_tyr-L_c + 2.0*LMPD_69_tyr-L_c_reverse_03f73 + 1.0*PPCK - 1.0*PPCK_reverse_2557d - 1.0*PYK + 1.0*PYK_reverse_bc8ff = 0.0
2pg_c: -1.0*ENO + 1.0*ENO_reverse_40eea - 1.0*PGM + 1.0*PGM_reverse_fc9af = 0.0
ergst_r: -1.0*ERGSTter + 1.0*ERGSTter_reverse_23bab + 1.0*LMPD_25_ergst_c - 1.0*LMPD_25_ergst_c_reverse_96537 + 1.0*LMPD_26_ergst_c - 1.0*LMPD_26_ergst_c_reverse_a59e0 = 0.0
dhap_c: 1.0*FBA + 1.0*FBA3 - 1.0*FBA3_reverse_0d49f - 1.0*FBA_reverse_84806 - 100.0*LMPD_47_pa_SC_c + 100.0*LMPD_47_pa_SC_c_reverse_4049b - 100.0*LMPD_48_pc_SC_c + 100.0*LMPD_48_pc_SC_c_reverse_be709 - 100.0*LMPD_49_pc_SC_c + 100.0*LMPD_49_pc_SC_c_reverse_55f95 - 100.0*LMPD_51_pe_SC_c + 100.0*LMPD_51_pe_SC_c_reverse_09777 - 100.0*LMPD_56_ps_SC_c + 100.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 100.0*LMPD_57_ptd1ino_SC_c + 100.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 100.0*LMPD_58_ptd1ino_SC_c + 100.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 100.0*LMPD_65_triglyc_SC_c + 100.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 100.0*LMPD_66_triglyc_SC_c + 100.0*LMPD_66_triglyc_SC_c_reverse_42d98 - 1.0*TPI + 1.0*TPI_reverse_c2c3b = 0.0
fdp_c: -1.0*FBA + 1.0*FBA_reverse_84806 - 1.0*FBP + 1.0*FBP_reverse_bf2c9 + 1.0*PFK - 1.0*PFK_reverse_d24a6 = 0.0
g3p_c: 1.0*FBA - 1.0*FBA_reverse_84806 - 1.0*GAPD + 1.0*GAPD_reverse_459c1 + 1.0*LMPD_67_trp-L_c - 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*TALA + 1.0*TALA_reverse_adfda + 1.0*TKT1 - 1.0*TKT1_reverse_a1021 + 1.0*TKT2 - 1.0*TKT2_reverse_7ebc7 + 1.0*TPI - 1.0*TPI_reverse_c2c3b = 0.0
s17bp_c: -1.0*FBA3 + 1.0*FBA3_reverse_0d49f + 1.0*PFK_3 - 1.0*PFK_3_reverse_037e0 = 0.0
e4p_c: 1.0*FBA3 - 1.0*FBA3_reverse_0d49f - 1.0*LMPD_52_phe-L_c + 1.0*LMPD_52_phe-L_c_reverse_909a1 - 1.0*LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_68_tyr-L_c + 1.0*LMPD_68_tyr-L_c_reverse_6983e - 1.0*LMPD_69_tyr-L_c + 1.0*LMPD_69_tyr-L_c_reverse_03f73 + 1.0*TALA - 1.0*TALA_reverse_adfda - 1.0*TKT2 + 1.0*TKT2_reverse_7ebc7 = 0.0
f6p_c: 1.0*FBP - 1.0*FBP_reverse_bf2c9 + 1.0*G6PI3 - 1.0*G6PI3_reverse_734c4 + 1.0*MAN6PI - 1.0*MAN6PI_reverse_d96f0 - 1.0*PFK + 1.0*PFK_reverse_d24a6 + 1.0*PGI - 1.0*PGI_reverse_27efc + 1.0*TALA - 1.0*TALA_reverse_adfda + 1.0*TKT2 - 1.0*TKT2_reverse_7ebc7 = 0.0
fum_c: -1.0*FRDcm + 1.0*FRDcm_reverse_6cbef - 1.0*FUM + 1.0*FUM_reverse_d3642 + 1.0*FUMt2r - 1.0*FUMt2r_reverse_d93c0 - 1.0*LMPD_11_cmp_c + 1.0*LMPD_11_cmp_c_reverse_98fa8 + 2.0*LMPD_13_damp_c - 2.0*LMPD_13_damp_c_reverse_e961c + 2.0*LMPD_14_damp_c - 2.0*LMPD_14_damp_c_reverse_d5737 + 2.0*LMPD_15_damp_c - 2.0*LMPD_15_damp_c_reverse_33cc4 + 2.0*LMPD_16_damp_c - 2.0*LMPD_16_damp_c_reverse_be99f - 1.0*LMPD_17_dcmp_c + 1.0*LMPD_17_dcmp_c_reverse_5c6b0 + 1.0*LMPD_19_dgmp_c - 1.0*LMPD_19_dgmp_c_reverse_7f2b1 + 1.0*LMPD_20_dgmp_c - 1.0*LMPD_20_dgmp_c_reverse_e4f3a + 1.0*LMPD_21_dgmp_c - 1.0*LMPD_21_dgmp_c_reverse_65045 + 1.0*LMPD_22_dgmp_c - 1.0*LMPD_22_dgmp_c_reverse_ebc77 - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 + 2.0*LMPD_2_amp_c - 2.0*LMPD_2_amp_c_reverse_bde95 + 1.0*LMPD_30_gmp_c - 1.0*LMPD_30_gmp_c_reverse_1f79a + 1.0*LMPD_31_gmp_c - 1.0*LMPD_31_gmp_c_reverse_d0c02 + 1.0*LMPD_32_his-L_c - 1.0*LMPD_32_his-L_c_reverse_c5967 + 1.0*LMPD_33_his-L_c - 1.0*LMPD_33_his-L_c_reverse_0eb8b + 1.0*LMPD_34_his-L_c - 1.0*LMPD_34_his-L_c_reverse_e16a9 + 1.0*LMPD_35_his-L_c - 1.0*LMPD_35_his-L_c_reverse_5e623 + 2.0*LMPD_3_amp_c - 2.0*LMPD_3_amp_c_reverse_96e30 + 2.0*LMPD_4_amp_c - 2.0*LMPD_4_amp_c_reverse_7925f + 1.0*LMPD_59_ribflv_c - 1.0*LMPD_59_ribflv_c_reverse_01d9b + 2.0*LMPD_5_amp_c - 2.0*LMPD_5_amp_c_reverse_8e974 + 1.0*LMPD_60_ribflv_c - 1.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_6_arg-L_c - 1.0*LMPD_6_arg-L_c_reverse_33f31 - 1.0*LMPD_70_ump_c + 1.0*LMPD_70_ump_c_reverse_3dcab - 1.0*LMPD_9_cmp_c + 1.0*LMPD_9_cmp_c_reverse_85279 + 1.0*SUCFUMtm - 1.0*SUCFUMtm_reverse_99c3a = 0.0
fad_m: 1.0*FRDcm - 1.0*FRDcm_reverse_6cbef + 1.0*FRDm - 1.0*FRDm_reverse_ef6ae - 1.0*SUCD1m + 1.0*SUCD1m_reverse_60ba1 + 1.0*SUCD3-u6m - 1.0*SUCD3-u6m_reverse_0e152 = 0.0
succ_c: 1.0*FRDcm - 1.0*FRDcm_reverse_6cbef + 1.0*LMPD_11_cmp_c - 1.0*LMPD_11_cmp_c_reverse_98fa8 + 1.0*LMPD_17_dcmp_c - 1.0*LMPD_17_dcmp_c_reverse_5c6b0 + 1.0*LMPD_24_dtmp_c - 1.0*LMPD_24_dtmp_c_reverse_5d241 + 1.0*LMPD_28_gly_c - 1.0*LMPD_28_gly_c_reverse_77428 + 1.0*LMPD_59_ribflv_c - 1.0*LMPD_59_ribflv_c_reverse_01d9b + 1.0*LMPD_60_ribflv_c - 1.0*LMPD_60_ribflv_c_reverse_b7a89 + 1.0*LMPD_70_ump_c - 1.0*LMPD_70_ump_c_reverse_3dcab + 1.0*LMPD_9_cmp_c - 1.0*LMPD_9_cmp_c_reverse_85279 + 1.0*SUCCt2r - 1.0*SUCCt2r_reverse_e3abe - 1.0*SUCCtm + 1.0*SUCCtm_reverse_5f36b - 1.0*SUCFUMtm + 1.0*SUCFUMtm_reverse_99c3a = 0.0
fadh2_m: -1.0*FRDcm + 1.0*FRDcm_reverse_6cbef - 1.0*FRDm + 1.0*FRDm_reverse_ef6ae + 1.0*SUCD1m - 1.0*SUCD1m_reverse_60ba1 - 1.0*SUCD3-u6m + 1.0*SUCD3-u6m_reverse_0e152 = 0.0
succ_m: 1.0*FRDm - 1.0*FRDm_reverse_ef6ae + 1.0*SUCCtm - 1.0*SUCCtm_reverse_5f36b - 1.0*SUCD1m + 1.0*SUCD1m_reverse_60ba1 - 1.0*SUCD2_u6m + 1.0*SUCD2_u6m_reverse_b6f79 + 1.0*SUCFUMtm - 1.0*SUCFUMtm_reverse_99c3a - 1.0*SUCOASm + 1.0*SUCOASm_reverse_87069 = 0.0
fum_m: -1.0*FRDm + 1.0*FRDm_reverse_ef6ae - 1.0*FUMm + 1.0*FUMm_reverse_5e620 + 1.0*SUCD1m - 1.0*SUCD1m_reverse_60ba1 + 1.0*SUCD2_u6m - 1.0*SUCD2_u6m_reverse_b6f79 - 1.0*SUCFUMtm + 1.0*SUCFUMtm_reverse_99c3a = 0.0
g6p_c: -1.0*G6PDH2 + 1.0*G6PDH2_reverse_47bd4 - 1.0*G6PI + 1.0*G6PI_reverse_834e6 - 1.0*G6Pter + 1.0*G6Pter_reverse_b4f95 + 1.0*HEX1 - 1.0*HEX1_reverse_25efa - 100.0*LMPD_57_ptd1ino_SC_c + 100.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 100.0*LMPD_58_ptd1ino_SC_c + 100.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 1.0*LMPD_64_tre_c + 1.0*LMPD_64_tre_c_reverse_64f54 - 1.0*PGI + 1.0*PGI_reverse_27efc + 1.0*PGMT - 1.0*PGMT_reverse_5bcdd = 0.0
g6p_r: -1.0*G6PDH2er + 1.0*G6PDH2er_reverse_9f5fd + 1.0*G6Pter - 1.0*G6Pter_reverse_b4f95 = 0.0
h_r: 1.0*G6PDH2er - 1.0*G6PDH2er_reverse_9f5fd - 2.0*LMPD_25_ergst_c + 2.0*LMPD_25_ergst_c_reverse_96537 - 2.0*LMPD_26_ergst_c + 2.0*LMPD_26_ergst_c_reverse_a59e0 - 1.0*LMPD_73_zymst_c + 1.0*LMPD_73_zymst_c_reverse_e6142 - 1.0*LMPD_74_zymst_c + 1.0*LMPD_74_zymst_c_reverse_84617 = 0.0
nadp_r: -1.0*G6PDH2er + 1.0*G6PDH2er_reverse_9f5fd + 2.0*LMPD_25_ergst_c - 2.0*LMPD_25_ergst_c_reverse_96537 + 2.0*LMPD_26_ergst_c - 2.0*LMPD_26_ergst_c_reverse_a59e0 + 1.0*LMPD_73_zymst_c - 1.0*LMPD_73_zymst_c_reverse_e6142 + 1.0*LMPD_74_zymst_c - 1.0*LMPD_74_zymst_c_reverse_84617 = 0.0
nadph_r: 1.0*G6PDH2er - 1.0*G6PDH2er_reverse_9f5fd - 2.0*LMPD_25_ergst_c + 2.0*LMPD_25_ergst_c_reverse_96537 - 2.0*LMPD_26_ergst_c + 2.0*LMPD_26_ergst_c_reverse_a59e0 - 1.0*LMPD_73_zymst_c + 1.0*LMPD_73_zymst_c_reverse_e6142 - 1.0*LMPD_74_zymst_c + 1.0*LMPD_74_zymst_c_reverse_84617 = 0.0
g6p-B_c: 1.0*G6PI - 1.0*G6PI3 + 1.0*G6PI3_reverse_734c4 - 1.0*G6PI_reverse_834e6 + 1.0*GLUK - 1.0*GLUK_reverse_ecb52 = 0.0
13dpg_c: 1.0*GAPD - 1.0*GAPD_reverse_459c1 + 1.0*PGK - 1.0*PGK_reverse_02696 = 0.0
glc-D_c: 1.0*GLCt1 - 1.0*GLCt1_reverse_3bcc5 - 1.0*GLUK + 1.0*GLUK_reverse_ecb52 - 1.0*HEX1 + 1.0*HEX1_reverse_25efa = 0.0
nh4_c: 1.0*GLUDxi - 1.0*GLUDxi_reverse_6fbe2 - 1.0*GLUDyi + 1.0*GLUDyi_reverse_a2351 - 2.0*LMPD_10_cmp_c + 2.0*LMPD_10_cmp_c_reverse_3acfa - 2.0*LMPD_11_cmp_c + 2.0*LMPD_11_cmp_c_reverse_98fa8 - 2.0*LMPD_13_damp_c + 2.0*LMPD_13_damp_c_reverse_e961c - 2.0*LMPD_14_damp_c + 2.0*LMPD_14_damp_c_reverse_d5737 - 2.0*LMPD_15_damp_c + 2.0*LMPD_15_damp_c_reverse_33cc4 - 2.0*LMPD_16_damp_c + 2.0*LMPD_16_damp_c_reverse_be99f - 2.0*LMPD_17_dcmp_c + 2.0*LMPD_17_dcmp_c_reverse_5c6b0 - 2.0*LMPD_18_dcmp_c + 2.0*LMPD_18_dcmp_c_reverse_15021 - 3.0*LMPD_19_dgmp_c + 3.0*LMPD_19_dgmp_c_reverse_7f2b1 - 3.0*LMPD_20_dgmp_c + 3.0*LMPD_20_dgmp_c_reverse_e4f3a - 3.0*LMPD_21_dgmp_c + 3.0*LMPD_21_dgmp_c_reverse_65045 - 3.0*LMPD_22_dgmp_c + 3.0*LMPD_22_dgmp_c_reverse_ebc77 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 - 1.0*LMPD_27_gln-L_c + 1.0*LMPD_27_gln-L_c_reverse_77507 - 2.0*LMPD_2_amp_c + 2.0*LMPD_2_amp_c_reverse_bde95 - 3.0*LMPD_30_gmp_c + 3.0*LMPD_30_gmp_c_reverse_1f79a - 3.0*LMPD_31_gmp_c + 3.0*LMPD_31_gmp_c_reverse_d0c02 - 1.0*LMPD_32_his-L_c + 1.0*LMPD_32_his-L_c_reverse_c5967 - 1.0*LMPD_33_his-L_c + 1.0*LMPD_33_his-L_c_reverse_0eb8b - 1.0*LMPD_34_his-L_c + 1.0*LMPD_34_his-L_c_reverse_e16a9 - 1.0*LMPD_35_his-L_c + 1.0*LMPD_35_his-L_c_reverse_5e623 - 2.0*LMPD_3_amp_c + 2.0*LMPD_3_amp_c_reverse_96e30 - 2.0*LMPD_4_amp_c + 2.0*LMPD_4_amp_c_reverse_7925f - 2.0*LMPD_59_ribflv_c + 2.0*LMPD_59_ribflv_c_reverse_01d9b - 2.0*LMPD_5_amp_c + 2.0*LMPD_5_amp_c_reverse_8e974 - 2.0*LMPD_60_ribflv_c + 2.0*LMPD_60_ribflv_c_reverse_b7a89 - 1.0*LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c_reverse_33f31 - 1.0*LMPD_70_ump_c + 1.0*LMPD_70_ump_c_reverse_3dcab - 1.0*LMPD_71_ump_c + 1.0*LMPD_71_ump_c_reverse_dc9d0 - 1.0*LMPD_7_asn-L_c + 1.0*LMPD_7_asn-L_c_reverse_b41c7 - 2.0*LMPD_8_cmp_c + 2.0*LMPD_8_cmp_c_reverse_23cbb - 2.0*LMPD_9_cmp_c + 2.0*LMPD_9_cmp_c_reverse_85279 + 1.0*NH4t - 1.0*NH4t_reverse_551ee - 1.0*NH4tm + 1.0*NH4tm_reverse_e28df = 0.0
6pgc_c: -1.0*GND + 1.0*GND_reverse_eec5c + 1.0*PGL - 1.0*PGL_reverse_2bb6b = 0.0
h2o_r: 1.0*H2Oter - 1.0*H2Oter_reverse_c874d + 1.0*LMPD_25_ergst_c - 1.0*LMPD_25_ergst_c_reverse_96537 + 1.0*LMPD_26_ergst_c - 1.0*LMPD_26_ergst_c_reverse_a59e0 + 1.0*LMPD_73_zymst_c - 1.0*LMPD_73_zymst_c_reverse_e6142 + 1.0*LMPD_74_zymst_c - 1.0*LMPD_74_zymst_c_reverse_84617 = 0.0
hcit_m: 1.0*HCITSm - 1.0*HCITSm_reverse_145c0 - 1.0*LMPD_40_lys-L_c + 1.0*LMPD_40_lys-L_c_reverse_4d3a6 - 1.0*LMPD_41_lys-L_c + 1.0*LMPD_41_lys-L_c_reverse_c021a = 0.0
hco3_c: 1.0*HCO3E - 1.0*HCO3E_reverse_97ea5 - 1.0*LMPD_10_cmp_c + 1.0*LMPD_10_cmp_c_reverse_3acfa - 1.0*LMPD_11_cmp_c + 1.0*LMPD_11_cmp_c_reverse_98fa8 - 1.0*LMPD_17_dcmp_c + 1.0*LMPD_17_dcmp_c_reverse_5c6b0 - 1.0*LMPD_18_dcmp_c + 1.0*LMPD_18_dcmp_c_reverse_15021 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 - 1420.0*LMPD_47_pa_SC_c + 1420.0*LMPD_47_pa_SC_c_reverse_4049b - 1420.0*LMPD_48_pc_SC_c + 1420.0*LMPD_48_pc_SC_c_reverse_be709 - 1420.0*LMPD_49_pc_SC_c + 1420.0*LMPD_49_pc_SC_c_reverse_55f95 - 1420.0*LMPD_51_pe_SC_c + 1420.0*LMPD_51_pe_SC_c_reverse_09777 - 1420.0*LMPD_56_ps_SC_c + 1420.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 1420.0*LMPD_57_ptd1ino_SC_c + 1420.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 1420.0*LMPD_58_ptd1ino_SC_c + 1420.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 2130.0*LMPD_65_triglyc_SC_c + 2130.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 2130.0*LMPD_66_triglyc_SC_c + 2130.0*LMPD_66_triglyc_SC_c_reverse_42d98 - 1.0*LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c_reverse_33f31 - 1.0*LMPD_70_ump_c + 1.0*LMPD_70_ump_c_reverse_3dcab - 1.0*LMPD_71_ump_c + 1.0*LMPD_71_ump_c_reverse_dc9d0 - 1.0*LMPD_8_cmp_c + 1.0*LMPD_8_cmp_c_reverse_23cbb - 1.0*LMPD_9_cmp_c + 1.0*LMPD_9_cmp_c_reverse_85279 - 1.0*PC + 1.0*PC_reverse_88dba = 0.0
ile-L_m: -1.0*ILEtmi + 1.0*ILEtmi_reverse_5c3e7 + 1.0*LMPD_36_ile-L_c - 1.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 = 0.0
man6p_c: -1.0*LMPD_42_mannan_c + 1.0*LMPD_42_mannan_c_reverse_9f626 - 1.0*LMPD_43_mannan_c + 1.0*LMPD_43_mannan_c_reverse_7dd82 - 1.0*MAN6PI + 1.0*MAN6PI_reverse_d96f0 = 0.0
mannan_r: 1.0*LMPD_42_mannan_c - 1.0*LMPD_42_mannan_c_reverse_9f626 + 1.0*LMPD_43_mannan_c - 1.0*LMPD_43_mannan_c_reverse_7dd82 + 1.0*MANNANter - 1.0*MANNANter_reverse_6b4bd = 0.0
nh4_m: 1.0*LMPD_36_ile-L_c - 1.0*LMPD_36_ile-L_c_reverse_d2bf7 + 1.0*LMPD_37_ile-L_c - 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*NH4tm - 1.0*NH4tm_reverse_e28df = 0.0
o2_c: -10.0*LMPD_25_ergst_c + 10.0*LMPD_25_ergst_c_reverse_96537 - 10.0*LMPD_26_ergst_c + 10.0*LMPD_26_ergst_c_reverse_a59e0 - 118.0*LMPD_47_pa_SC_c + 118.0*LMPD_47_pa_SC_c_reverse_4049b - 118.0*LMPD_48_pc_SC_c + 118.0*LMPD_48_pc_SC_c_reverse_be709 - 118.0*LMPD_49_pc_SC_c + 118.0*LMPD_49_pc_SC_c_reverse_55f95 - 118.0*LMPD_51_pe_SC_c + 118.0*LMPD_51_pe_SC_c_reverse_09777 - 118.0*LMPD_56_ps_SC_c + 118.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 118.0*LMPD_57_ptd1ino_SC_c + 118.0*LMPD_57_ptd1ino_SC_c_reverse_65089 - 118.0*LMPD_58_ptd1ino_SC_c + 118.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 - 177.0*LMPD_65_triglyc_SC_c + 177.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 - 177.0*LMPD_66_triglyc_SC_c + 177.0*LMPD_66_triglyc_SC_c_reverse_42d98 - 9.0*LMPD_73_zymst_c + 9.0*LMPD_73_zymst_c_reverse_e6142 - 9.0*LMPD_74_zymst_c + 9.0*LMPD_74_zymst_c_reverse_84617 + 1.0*O2t - 1.0*O2t_reverse_4d957 - 1.0*O2ter + 1.0*O2ter_reverse_ea1ad - 1.0*O2tm + 1.0*O2tm_reverse_fe4d6 = 0.0
o2_r: -1.0*LMPD_25_ergst_c + 1.0*LMPD_25_ergst_c_reverse_96537 - 1.0*LMPD_26_ergst_c + 1.0*LMPD_26_ergst_c_reverse_a59e0 - 1.0*LMPD_73_zymst_c + 1.0*LMPD_73_zymst_c_reverse_e6142 - 1.0*LMPD_74_zymst_c + 1.0*LMPD_74_zymst_c_reverse_84617 + 1.0*O2ter - 1.0*O2ter_reverse_ea1ad = 0.0
pe_SC_m: 1.0*LMPD_48_pc_SC_c - 1.0*LMPD_48_pc_SC_c_reverse_be709 + 1.0*LMPD_49_pc_SC_c - 1.0*LMPD_49_pc_SC_c_reverse_55f95 + 1.0*LMPD_51_pe_SC_c - 1.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*PEtm_SC - 1.0*PEtm_SC_reverse_a93d3 = 0.0
s7p_c: -1.0*PFK_3 + 1.0*PFK_3_reverse_037e0 - 1.0*TALA + 1.0*TALA_reverse_adfda + 1.0*TKT1 - 1.0*TKT1_reverse_a1021 = 0.0
3pg_c: -1.0*LMPD_12_cys-L_c + 1.0*LMPD_12_cys-L_c_reverse_6dc0a - 2.0*LMPD_13_damp_c + 2.0*LMPD_13_damp_c_reverse_e961c - 2.0*LMPD_14_damp_c + 2.0*LMPD_14_damp_c_reverse_d5737 - 2.0*LMPD_15_damp_c + 2.0*LMPD_15_damp_c_reverse_33cc4 - 2.0*LMPD_16_damp_c + 2.0*LMPD_16_damp_c_reverse_be99f - 2.0*LMPD_19_dgmp_c + 2.0*LMPD_19_dgmp_c_reverse_7f2b1 - 2.0*LMPD_20_dgmp_c + 2.0*LMPD_20_dgmp_c_reverse_e4f3a - 2.0*LMPD_21_dgmp_c + 2.0*LMPD_21_dgmp_c_reverse_65045 - 2.0*LMPD_22_dgmp_c + 2.0*LMPD_22_dgmp_c_reverse_ebc77 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 - 2.0*LMPD_2_amp_c + 2.0*LMPD_2_amp_c_reverse_bde95 - 2.0*LMPD_30_gmp_c + 2.0*LMPD_30_gmp_c_reverse_1f79a - 2.0*LMPD_31_gmp_c + 2.0*LMPD_31_gmp_c_reverse_d0c02 - 1.0*LMPD_32_his-L_c + 1.0*LMPD_32_his-L_c_reverse_c5967 - 1.0*LMPD_33_his-L_c + 1.0*LMPD_33_his-L_c_reverse_0eb8b - 1.0*LMPD_34_his-L_c + 1.0*LMPD_34_his-L_c_reverse_e16a9 - 1.0*LMPD_35_his-L_c + 1.0*LMPD_35_his-L_c_reverse_5e623 - 2.0*LMPD_3_amp_c + 2.0*LMPD_3_amp_c_reverse_96e30 - 1.0*LMPD_44_met-L_c + 1.0*LMPD_44_met-L_c_reverse_1646f - 1.0*LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c_reverse_912f1 - 400.0*LMPD_48_pc_SC_c + 400.0*LMPD_48_pc_SC_c_reverse_be709 - 400.0*LMPD_49_pc_SC_c + 400.0*LMPD_49_pc_SC_c_reverse_55f95 - 2.0*LMPD_4_amp_c + 2.0*LMPD_4_amp_c_reverse_7925f - 100.0*LMPD_51_pe_SC_c + 100.0*LMPD_51_pe_SC_c_reverse_09777 - 100.0*LMPD_56_ps_SC_c + 100.0*LMPD_56_ps_SC_c_reverse_7e9d7 - 2.0*LMPD_5_amp_c + 2.0*LMPD_5_amp_c_reverse_8e974 - 1.0*LMPD_61_ser-L_c + 1.0*LMPD_61_ser-L_c_reverse_40666 - 1.0*LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*PGK + 1.0*PGK_reverse_02696 + 1.0*PGM - 1.0*PGM_reverse_fc9af = 0.0
g1p_c: -1.0*LMPD_1_13BDglcn_c + 1.0*LMPD_1_13BDglcn_c_reverse_e2bdf - 1.0*LMPD_29_glycogen_c + 1.0*LMPD_29_glycogen_c_reverse_68cb3 - 1.0*LMPD_64_tre_c + 1.0*LMPD_64_tre_c_reverse_64f54 - 1.0*PGMT + 1.0*PGMT_reverse_5bcdd = 0.0
ps_SC_m: -1.0*LMPD_48_pc_SC_c + 1.0*LMPD_48_pc_SC_c_reverse_be709 - 1.0*LMPD_49_pc_SC_c + 1.0*LMPD_49_pc_SC_c_reverse_55f95 - 1.0*LMPD_51_pe_SC_c + 1.0*LMPD_51_pe_SC_c_reverse_09777 + 1.0*PStm_SC - 1.0*PStm_SC_reverse_a026e = 0.0
xu5p-D_c: 1.0*RPE - 1.0*RPE_reverse_a1b04 - 1.0*TKT1 + 1.0*TKT1_reverse_a1021 - 1.0*TKT2 + 1.0*TKT2_reverse_7ebc7 = 0.0
r5p_c: -1.0*LMPD_10_cmp_c + 1.0*LMPD_10_cmp_c_reverse_3acfa - 1.0*LMPD_11_cmp_c + 1.0*LMPD_11_cmp_c_reverse_98fa8 - 1.0*LMPD_13_damp_c + 1.0*LMPD_13_damp_c_reverse_e961c - 1.0*LMPD_14_damp_c + 1.0*LMPD_14_damp_c_reverse_d5737 - 1.0*LMPD_15_damp_c + 1.0*LMPD_15_damp_c_reverse_33cc4 - 1.0*LMPD_16_damp_c + 1.0*LMPD_16_damp_c_reverse_be99f - 1.0*LMPD_17_dcmp_c + 1.0*LMPD_17_dcmp_c_reverse_5c6b0 - 1.0*LMPD_18_dcmp_c + 1.0*LMPD_18_dcmp_c_reverse_15021 - 1.0*LMPD_19_dgmp_c + 1.0*LMPD_19_dgmp_c_reverse_7f2b1 - 1.0*LMPD_20_dgmp_c + 1.0*LMPD_20_dgmp_c_reverse_e4f3a - 1.0*LMPD_21_dgmp_c + 1.0*LMPD_21_dgmp_c_reverse_65045 - 1.0*LMPD_22_dgmp_c + 1.0*LMPD_22_dgmp_c_reverse_ebc77 - 1.0*LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d - 1.0*LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 - 1.0*LMPD_2_amp_c + 1.0*LMPD_2_amp_c_reverse_bde95 - 1.0*LMPD_30_gmp_c + 1.0*LMPD_30_gmp_c_reverse_1f79a - 1.0*LMPD_31_gmp_c + 1.0*LMPD_31_gmp_c_reverse_d0c02 - 1.0*LMPD_32_his-L_c + 1.0*LMPD_32_his-L_c_reverse_c5967 - 1.0*LMPD_33_his-L_c + 1.0*LMPD_33_his-L_c_reverse_0eb8b - 1.0*LMPD_34_his-L_c + 1.0*LMPD_34_his-L_c_reverse_e16a9 - 1.0*LMPD_35_his-L_c + 1.0*LMPD_35_his-L_c_reverse_5e623 - 1.0*LMPD_3_amp_c + 1.0*LMPD_3_amp_c_reverse_96e30 - 1.0*LMPD_4_amp_c + 1.0*LMPD_4_amp_c_reverse_7925f - 1.0*LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b - 1.0*LMPD_5_amp_c + 1.0*LMPD_5_amp_c_reverse_8e974 - 1.0*LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 - 1.0*LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 - 1.0*LMPD_70_ump_c + 1.0*LMPD_70_ump_c_reverse_3dcab - 1.0*LMPD_71_ump_c + 1.0*LMPD_71_ump_c_reverse_dc9d0 - 1.0*LMPD_8_cmp_c + 1.0*LMPD_8_cmp_c_reverse_23cbb - 1.0*LMPD_9_cmp_c + 1.0*LMPD_9_cmp_c_reverse_85279 - 1.0*RPI + 1.0*RPI_reverse_853a1 - 1.0*TKT1 + 1.0*TKT1_reverse_a1021 = 0.0
thr-L_m: -1.0*LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 - 1.0*LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 + 1.0*THRt2m - 1.0*THRt2m_reverse_fd5a8 = 0.0
SU_biomass_SC5_notrace: 0 < 1.0*BU_biomass_SC5_notrace + 1.0*FU_biomass_SC5_notrace < 1
UF_biomass_SC5_notrace: -1000.0*FU_biomass_SC5_notrace + 1.0*biomass_SC5_notrace < 0
UR_biomass_SC5_notrace: -1000.0*BU_biomass_SC5_notrace + 1.0*biomass_SC5_notrace_reverse_e32ff < 0
SU_EX_4abut_e: 0 < 1.0*BU_EX_4abut_e + 1.0*FU_EX_4abut_e < 1
UF_EX_4abut_e: -1000.0*FU_EX_4abut_e + 1.0*EX_4abut_e < 0
UR_EX_4abut_e: -1000.0*BU_EX_4abut_e + 1.0*EX_4abut_e_reverse_82295 < 0
SU_EX_ac_e: 0 < 1.0*BU_EX_ac_e + 1.0*FU_EX_ac_e < 1
UF_EX_ac_e: -1000.0*FU_EX_ac_e + 1.0*EX_ac_e < 0
UR_EX_ac_e: -1000.0*BU_EX_ac_e + 1.0*EX_ac_e_reverse_0be96 < 0
SU_EX_acald_e: 0 < 1.0*BU_EX_acald_e + 1.0*FU_EX_acald_e < 1
UF_EX_acald_e: -1000.0*FU_EX_acald_e + 1.0*EX_acald_e < 0
UR_EX_acald_e: -1000.0*BU_EX_acald_e + 1.0*EX_acald_e_reverse_c096e < 0
SU_EX_akg_e: 0 < 1.0*BU_EX_akg_e + 1.0*FU_EX_akg_e < 1
UF_EX_akg_e: -1000.0*FU_EX_akg_e + 1.0*EX_akg_e < 0
UR_EX_akg_e: -1000.0*BU_EX_akg_e + 1.0*EX_akg_e_reverse_70d85 < 0
SU_EX_ala-L_e: 0 < 1.0*BU_EX_ala-L_e + 1.0*FU_EX_ala-L_e < 1
UF_EX_ala-L_e: -1000.0*FU_EX_ala-L_e + 1.0*EX_ala-L_e < 0
UR_EX_ala-L_e: -1000.0*BU_EX_ala-L_e + 1.0*EX_ala-L_e_reverse_45174 < 0
SU_EX_arg-L_e: 0 < 1.0*BU_EX_arg-L_e + 1.0*FU_EX_arg-L_e < 1
UF_EX_arg-L_e: -1000.0*FU_EX_arg-L_e + 1.0*EX_arg-L_e < 0
UR_EX_arg-L_e: -1000.0*BU_EX_arg-L_e + 1.0*EX_arg-L_e_reverse_2768b < 0
SU_EX_asn-L_e: 0 < 1.0*BU_EX_asn-L_e + 1.0*FU_EX_asn-L_e < 1
UF_EX_asn-L_e: -1000.0*FU_EX_asn-L_e + 1.0*EX_asn-L_e < 0
UR_EX_asn-L_e: -1000.0*BU_EX_asn-L_e + 1.0*EX_asn-L_e_reverse_73ed2 < 0
SU_EX_asp-L_e: 0 < 1.0*BU_EX_asp-L_e + 1.0*FU_EX_asp-L_e < 1
UF_EX_asp-L_e: -1000.0*FU_EX_asp-L_e + 1.0*EX_asp-L_e < 0
UR_EX_asp-L_e: -1000.0*BU_EX_asp-L_e + 1.0*EX_asp-L_e_reverse_ceee6 < 0
SU_EX_cit_e: 0 < 1.0*BU_EX_cit_e + 1.0*FU_EX_cit_e < 1
UF_EX_cit_e: -1000.0*FU_EX_cit_e + 1.0*EX_cit_e < 0
UR_EX_cit_e: -1000.0*BU_EX_cit_e + 1.0*EX_cit_e_reverse_0835e < 0
SU_EX_co2_e: 0 < 1.0*BU_EX_co2_e + 1.0*FU_EX_co2_e < 1
UF_EX_co2_e: -1000.0*FU_EX_co2_e + 1.0*EX_co2_e < 0
UR_EX_co2_e: -1000.0*BU_EX_co2_e + 1.0*EX_co2_e_reverse_d0466 < 0
SU_EX_cys-L_e: 0 < 1.0*BU_EX_cys-L_e + 1.0*FU_EX_cys-L_e < 1
UF_EX_cys-L_e: -1000.0*FU_EX_cys-L_e + 1.0*EX_cys-L_e < 0
UR_EX_cys-L_e: -1000.0*BU_EX_cys-L_e + 1.0*EX_cys-L_e_reverse_c0cbc < 0
SU_EX_ergst_e: 0 < 1.0*BU_EX_ergst_e + 1.0*FU_EX_ergst_e < 1
UF_EX_ergst_e: -1000.0*FU_EX_ergst_e + 1.0*EX_ergst_e < 0
UR_EX_ergst_e: -1000.0*BU_EX_ergst_e + 1.0*EX_ergst_e_reverse_92efd < 0
SU_EX_etoh_e: 0 < 1.0*BU_EX_etoh_e + 1.0*FU_EX_etoh_e < 1
UF_EX_etoh_e: -1000.0*FU_EX_etoh_e + 1.0*EX_etoh_e < 0
UR_EX_etoh_e: -1000.0*BU_EX_etoh_e + 1.0*EX_etoh_e_reverse_cc64f < 0
SU_EX_for_e: 0 < 1.0*BU_EX_for_e + 1.0*FU_EX_for_e < 1
UF_EX_for_e: -1000.0*FU_EX_for_e + 1.0*EX_for_e < 0
UR_EX_for_e: -1000.0*BU_EX_for_e + 1.0*EX_for_e_reverse_23269 < 0
SU_EX_fum_e: 0 < 1.0*BU_EX_fum_e + 1.0*FU_EX_fum_e < 1
UF_EX_fum_e: -1000.0*FU_EX_fum_e + 1.0*EX_fum_e < 0
UR_EX_fum_e: -1000.0*BU_EX_fum_e + 1.0*EX_fum_e_reverse_e3432 < 0
SU_EX_glc_e: 0 < 1.0*BU_EX_glc_e + 1.0*FU_EX_glc_e < 1
UF_EX_glc_e: -1000.0*FU_EX_glc_e + 1.0*EX_glc_e < 0
UR_EX_glc_e: -1000.0*BU_EX_glc_e + 1.0*EX_glc_e_reverse_6226d < 0
SU_EX_gln-L_e: 0 < 1.0*BU_EX_gln-L_e + 1.0*FU_EX_gln-L_e < 1
UF_EX_gln-L_e: -1000.0*FU_EX_gln-L_e + 1.0*EX_gln-L_e < 0
UR_EX_gln-L_e: -1000.0*BU_EX_gln-L_e + 1.0*EX_gln-L_e_reverse_1be13 < 0
SU_EX_glu-L_e: 0 < 1.0*BU_EX_glu-L_e + 1.0*FU_EX_glu-L_e < 1
UF_EX_glu-L_e: -1000.0*FU_EX_glu-L_e + 1.0*EX_glu-L_e < 0
UR_EX_glu-L_e: -1000.0*BU_EX_glu-L_e + 1.0*EX_glu-L_e_reverse_727e6 < 0
SU_EX_gly_e: 0 < 1.0*BU_EX_gly_e + 1.0*FU_EX_gly_e < 1
UF_EX_gly_e: -1000.0*FU_EX_gly_e + 1.0*EX_gly_e < 0
UR_EX_gly_e: -1000.0*BU_EX_gly_e + 1.0*EX_gly_e_reverse_6956b < 0
SU_EX_h_e: 0 < 1.0*BU_EX_h_e + 1.0*FU_EX_h_e < 1
UF_EX_h_e: -1000.0*FU_EX_h_e + 1.0*EX_h_e < 0
UR_EX_h_e: -1000.0*BU_EX_h_e + 1.0*EX_h_e_reverse_3e0c5 < 0
SU_EX_h2o_e: 0 < 1.0*BU_EX_h2o_e + 1.0*FU_EX_h2o_e < 1
UF_EX_h2o_e: -1000.0*FU_EX_h2o_e + 1.0*EX_h2o_e < 0
UR_EX_h2o_e: -1000.0*BU_EX_h2o_e + 1.0*EX_h2o_e_reverse_3ced4 < 0
SU_EX_his-L_e: 0 < 1.0*BU_EX_his-L_e + 1.0*FU_EX_his-L_e < 1
UF_EX_his-L_e: -1000.0*FU_EX_his-L_e + 1.0*EX_his-L_e < 0
UR_EX_his-L_e: -1000.0*BU_EX_his-L_e + 1.0*EX_his-L_e_reverse_e1c17 < 0
SU_EX_ile-L_e: 0 < 1.0*BU_EX_ile-L_e + 1.0*FU_EX_ile-L_e < 1
UF_EX_ile-L_e: -1000.0*FU_EX_ile-L_e + 1.0*EX_ile-L_e < 0
UR_EX_ile-L_e: -1000.0*BU_EX_ile-L_e + 1.0*EX_ile-L_e_reverse_e9c53 < 0
SU_EX_leu-L_e: 0 < 1.0*BU_EX_leu-L_e + 1.0*FU_EX_leu-L_e < 1
UF_EX_leu-L_e: -1000.0*FU_EX_leu-L_e + 1.0*EX_leu-L_e < 0
UR_EX_leu-L_e: -1000.0*BU_EX_leu-L_e + 1.0*EX_leu-L_e_reverse_9947f < 0
SU_EX_lys-L_e: 0 < 1.0*BU_EX_lys-L_e + 1.0*FU_EX_lys-L_e < 1
UF_EX_lys-L_e: -1000.0*FU_EX_lys-L_e + 1.0*EX_lys-L_e < 0
UR_EX_lys-L_e: -1000.0*BU_EX_lys-L_e + 1.0*EX_lys-L_e_reverse_5779d < 0
SU_EX_mal-L_e: 0 < 1.0*BU_EX_mal-L_e + 1.0*FU_EX_mal-L_e < 1
UF_EX_mal-L_e: -1000.0*FU_EX_mal-L_e + 1.0*EX_mal-L_e < 0
UR_EX_mal-L_e: -1000.0*BU_EX_mal-L_e + 1.0*EX_mal-L_e_reverse_00553 < 0
SU_EX_met-L_e: 0 < 1.0*BU_EX_met-L_e + 1.0*FU_EX_met-L_e < 1
UF_EX_met-L_e: -1000.0*FU_EX_met-L_e + 1.0*EX_met-L_e < 0
UR_EX_met-L_e: -1000.0*BU_EX_met-L_e + 1.0*EX_met-L_e_reverse_f8e12 < 0
SU_EX_nh4_e: 0 < 1.0*BU_EX_nh4_e + 1.0*FU_EX_nh4_e < 1
UF_EX_nh4_e: -1000.0*FU_EX_nh4_e + 1.0*EX_nh4_e < 0
UR_EX_nh4_e: -1000.0*BU_EX_nh4_e + 1.0*EX_nh4_e_reverse_f9cc6 < 0
SU_EX_o2_e: 0 < 1.0*BU_EX_o2_e + 1.0*FU_EX_o2_e < 1
UF_EX_o2_e: -1000.0*FU_EX_o2_e + 1.0*EX_o2_e < 0
UR_EX_o2_e: -1000.0*BU_EX_o2_e + 1.0*EX_o2_e_reverse_efa94 < 0
SU_EX_oaa_e: 0 < 1.0*BU_EX_oaa_e + 1.0*FU_EX_oaa_e < 1
UF_EX_oaa_e: -1000.0*FU_EX_oaa_e + 1.0*EX_oaa_e < 0
UR_EX_oaa_e: -1000.0*BU_EX_oaa_e + 1.0*EX_oaa_e_reverse_9b144 < 0
SU_EX_phe-L_e: 0 < 1.0*BU_EX_phe-L_e + 1.0*FU_EX_phe-L_e < 1
UF_EX_phe-L_e: -1000.0*FU_EX_phe-L_e + 1.0*EX_phe-L_e < 0
UR_EX_phe-L_e: -1000.0*BU_EX_phe-L_e + 1.0*EX_phe-L_e_reverse_155f7 < 0
SU_EX_pi_e: 0 < 1.0*BU_EX_pi_e + 1.0*FU_EX_pi_e < 1
UF_EX_pi_e: -1000.0*FU_EX_pi_e + 1.0*EX_pi_e < 0
UR_EX_pi_e: -1000.0*BU_EX_pi_e + 1.0*EX_pi_e_reverse_1fb09 < 0
SU_EX_pro-L_e: 0 < 1.0*BU_EX_pro-L_e + 1.0*FU_EX_pro-L_e < 1
UF_EX_pro-L_e: -1000.0*FU_EX_pro-L_e + 1.0*EX_pro-L_e < 0
UR_EX_pro-L_e: -1000.0*BU_EX_pro-L_e + 1.0*EX_pro-L_e_reverse_26581 < 0
SU_EX_pyr_e: 0 < 1.0*BU_EX_pyr_e + 1.0*FU_EX_pyr_e < 1
UF_EX_pyr_e: -1000.0*FU_EX_pyr_e + 1.0*EX_pyr_e < 0
UR_EX_pyr_e: -1000.0*BU_EX_pyr_e + 1.0*EX_pyr_e_reverse_1f6de < 0
SU_EX_ser-L_e: 0 < 1.0*BU_EX_ser-L_e + 1.0*FU_EX_ser-L_e < 1
UF_EX_ser-L_e: -1000.0*FU_EX_ser-L_e + 1.0*EX_ser-L_e < 0
UR_EX_ser-L_e: -1000.0*BU_EX_ser-L_e + 1.0*EX_ser-L_e_reverse_6b4e0 < 0
SU_EX_so4_e: 0 < 1.0*BU_EX_so4_e + 1.0*FU_EX_so4_e < 1
UF_EX_so4_e: -1000.0*FU_EX_so4_e + 1.0*EX_so4_e < 0
UR_EX_so4_e: -1000.0*BU_EX_so4_e + 1.0*EX_so4_e_reverse_5c8ed < 0
SU_EX_succ_e: 0 < 1.0*BU_EX_succ_e + 1.0*FU_EX_succ_e < 1
UF_EX_succ_e: -1000.0*FU_EX_succ_e + 1.0*EX_succ_e < 0
UR_EX_succ_e: -1000.0*BU_EX_succ_e + 1.0*EX_succ_e_reverse_a9039 < 0
SU_EX_thr-L_e: 0 < 1.0*BU_EX_thr-L_e + 1.0*FU_EX_thr-L_e < 1
UF_EX_thr-L_e: -1000.0*FU_EX_thr-L_e + 1.0*EX_thr-L_e < 0
UR_EX_thr-L_e: -1000.0*BU_EX_thr-L_e + 1.0*EX_thr-L_e_reverse_9a854 < 0
SU_EX_trp-L_e: 0 < 1.0*BU_EX_trp-L_e + 1.0*FU_EX_trp-L_e < 1
UF_EX_trp-L_e: -1000.0*FU_EX_trp-L_e + 1.0*EX_trp-L_e < 0
UR_EX_trp-L_e: -1000.0*BU_EX_trp-L_e + 1.0*EX_trp-L_e_reverse_273bc < 0
SU_EX_tyr-L_e: 0 < 1.0*BU_EX_tyr-L_e + 1.0*FU_EX_tyr-L_e < 1
UF_EX_tyr-L_e: -1000.0*FU_EX_tyr-L_e + 1.0*EX_tyr-L_e < 0
UR_EX_tyr-L_e: -1000.0*BU_EX_tyr-L_e + 1.0*EX_tyr-L_e_reverse_08f6e < 0
SU_EX_val-L_e: 0 < 1.0*BU_EX_val-L_e + 1.0*FU_EX_val-L_e < 1
UF_EX_val-L_e: -1000.0*FU_EX_val-L_e + 1.0*EX_val-L_e < 0
UR_EX_val-L_e: -1000.0*BU_EX_val-L_e + 1.0*EX_val-L_e_reverse_4c155 < 0
SU_EX_zymst_e: 0 < 1.0*BU_EX_zymst_e + 1.0*FU_EX_zymst_e < 1
UF_EX_zymst_e: -1000.0*FU_EX_zymst_e + 1.0*EX_zymst_e < 0
UR_EX_zymst_e: -1000.0*BU_EX_zymst_e + 1.0*EX_zymst_e_reverse_78276 < 0
G_6PGLter: -1.0*DG_6PGLter - 0.592089376625*LC_6pgl_c + 0.592089376625*LC_6pgl_r + 1.0*DGo_6PGLter = 0
FU_6PGLter: 1.0*DG_6PGLter + 1000.0*FU_6PGLter < 999.999999
BU_6PGLter: -1.0*DG_6PGLter + 1000.0*BU_6PGLter < 999.999999
SU_6PGLter: 0 < 1.0*BU_6PGLter + 1.0*FU_6PGLter < 1
UF_6PGLter: -1000.0*FU_6PGLter + 1.0*6PGLter < 0
UR_6PGLter: -1000.0*BU_6PGLter + 1.0*6PGLter_reverse_0e388 < 0
G_ABUTt2r: -1.0*DG_ABUTt2r + 0.592089376625*LC_4abut_c - 0.592089376625*LC_4abut_e + 1.0*DGo_ABUTt2r = 0
FU_ABUTt2r: 1.0*DG_ABUTt2r + 1000.0*FU_ABUTt2r < 999.999999
BU_ABUTt2r: -1.0*DG_ABUTt2r + 1000.0*BU_ABUTt2r < 999.999999
SU_ABUTt2r: 0 < 1.0*BU_ABUTt2r + 1.0*FU_ABUTt2r < 1
UF_ABUTt2r: -1000.0*FU_ABUTt2r + 1.0*ABUTt2r < 0
UR_ABUTt2r: -1000.0*BU_ABUTt2r + 1.0*ABUTt2r_reverse_d925e < 0
G_ACALDt: 1.0*DGo_ACALDt - 1.0*DG_ACALDt + 0.592089376625*LC_acald_c - 0.592089376625*LC_acald_e = 0
FU_ACALDt: 1.0*DG_ACALDt + 1000.0*FU_ACALDt < 999.999999
BU_ACALDt: -1.0*DG_ACALDt + 1000.0*BU_ACALDt < 999.999999
SU_ACALDt: 0 < 1.0*BU_ACALDt + 1.0*FU_ACALDt < 1
UF_ACALDt: -1000.0*FU_ACALDt + 1.0*ACALDt < 0
UR_ACALDt: -1000.0*BU_ACALDt + 1.0*ACALDt_reverse_858fa < 0
G_ACALDtm: -1.0*DG_ACALDtm + 0.592089376625*LC_acald_c - 0.592089376625*LC_acald_m + 1.0*DGo_ACALDtm = 0
FU_ACALDtm: 1.0*DG_ACALDtm + 1000.0*FU_ACALDtm < 999.999999
BU_ACALDtm: -1.0*DG_ACALDtm + 1000.0*BU_ACALDtm < 999.999999
SU_ACALDtm: 0 < 1.0*BU_ACALDtm + 1.0*FU_ACALDtm < 1
UF_ACALDtm: -1000.0*FU_ACALDtm + 1.0*ACALDtm < 0
UR_ACALDtm: -1000.0*BU_ACALDtm + 1.0*ACALDtm_reverse_c2973 < 0
G_ACOAHim: -1.0*DG_ACOAHim + 0.592089376625*LC_ac_m - 0.592089376625*LC_accoa_m + 0.592089376625*LC_coa_m + 1.0*DGo_ACOAHim = 0
FU_ACOAHim: 1.0*DG_ACOAHim + 1000.0*FU_ACOAHim < 999.999999
BU_ACOAHim: -1.0*DG_ACOAHim + 1000.0*BU_ACOAHim < 999.999999
SU_ACOAHim: 0 < 1.0*BU_ACOAHim + 1.0*FU_ACOAHim < 1
UF_ACOAHim: -1000.0*FU_ACOAHim + 1.0*ACOAHim < 0
UR_ACOAHim: -1000.0*BU_ACOAHim + 1.0*ACOAHim_reverse_069ca < 0
G_ACONT: -1.0*DG_ACONT - 0.592089376625*LC_cit_c + 0.592089376625*LC_icit_c + 1.0*DGo_ACONT = 0
FU_ACONT: 1.0*DG_ACONT + 1000.0*FU_ACONT < 999.999999
BU_ACONT: -1.0*DG_ACONT + 1000.0*BU_ACONT < 999.999999
SU_ACONT: 0 < 1.0*BU_ACONT + 1.0*FU_ACONT < 1
UF_ACONT: -1000.0*FU_ACONT + 1.0*ACONT < 0
UR_ACONT: -1000.0*BU_ACONT + 1.0*ACONT_reverse_7c2d5 < 0
G_ACONTm: -1.0*DG_ACONTm - 0.592089376625*LC_cit_m + 0.592089376625*LC_icit_m + 1.0*DGo_ACONTm = 0
FU_ACONTm: 1.0*DG_ACONTm + 1000.0*FU_ACONTm < 999.999999
BU_ACONTm: -1.0*DG_ACONTm + 1000.0*BU_ACONTm < 999.999999
SU_ACONTm: 0 < 1.0*BU_ACONTm + 1.0*FU_ACONTm < 1
UF_ACONTm: -1000.0*FU_ACONTm + 1.0*ACONTm < 0
UR_ACONTm: -1000.0*BU_ACONTm + 1.0*ACONTm_reverse_8079b < 0
G_ACS: -1.0*DG_ACS - 0.592089376625*LC_ac_c + 0.592089376625*LC_accoa_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_atp_c - 0.592089376625*LC_coa_c + 0.592089376625*LC_ppi_c + 1.0*DGo_ACS = 0
FU_ACS: 1.0*DG_ACS + 1000.0*FU_ACS < 999.999999
BU_ACS: -1.0*DG_ACS + 1000.0*BU_ACS < 999.999999
SU_ACS: 0 < 1.0*BU_ACS + 1.0*FU_ACS < 1
UF_ACS: -1000.0*FU_ACS + 1.0*ACS < 0
UR_ACS: -1000.0*BU_ACS + 1.0*ACS_reverse_37635 < 0
G_ACSm: -1.0*DG_ACSm - 0.592089376625*LC_ac_m + 0.592089376625*LC_accoa_m + 0.592089376625*LC_amp_m - 0.592089376625*LC_atp_m - 0.592089376625*LC_coa_m + 0.592089376625*LC_ppi_m + 1.0*DGo_ACSm = 0
FU_ACSm: 1.0*DG_ACSm + 1000.0*FU_ACSm < 999.999999
BU_ACSm: -1.0*DG_ACSm + 1000.0*BU_ACSm < 999.999999
SU_ACSm: 0 < 1.0*BU_ACSm + 1.0*FU_ACSm < 1
UF_ACSm: -1000.0*FU_ACSm + 1.0*ACSm < 0
UR_ACSm: -1000.0*BU_ACSm + 1.0*ACSm_reverse_5e300 < 0
G_ACt2r: -1.0*DG_ACt2r + 0.592089376625*LC_ac_c - 0.592089376625*LC_ac_e + 1.0*DGo_ACt2r = 0
FU_ACt2r: 1.0*DG_ACt2r + 1000.0*FU_ACt2r < 999.999999
BU_ACt2r: -1.0*DG_ACt2r + 1000.0*BU_ACt2r < 999.999999
SU_ACt2r: 0 < 1.0*BU_ACt2r + 1.0*FU_ACt2r < 1
UF_ACt2r: -1000.0*FU_ACt2r + 1.0*ACt2r < 0
UR_ACt2r: -1000.0*BU_ACt2r + 1.0*ACt2r_reverse_64e05 < 0
G_ACtm: -1.0*DG_ACtm - 0.592089376625*LC_ac_c + 0.592089376625*LC_ac_m + 1.0*DGo_ACtm = 0
FU_ACtm: 1.0*DG_ACtm + 1000.0*FU_ACtm < 999.999999
BU_ACtm: -1.0*DG_ACtm + 1000.0*BU_ACtm < 999.999999
SU_ACtm: 0 < 1.0*BU_ACtm + 1.0*FU_ACtm < 1
UF_ACtm: -1000.0*FU_ACtm + 1.0*ACtm < 0
UR_ACtm: -1000.0*BU_ACtm + 1.0*ACtm_reverse_5029e < 0
G_ACtr: -1.0*DG_ACtr + 0.592089376625*LC_ac_c - 0.592089376625*LC_ac_e + 1.0*DGo_ACtr = 0
FU_ACtr: 1.0*DG_ACtr + 1000.0*FU_ACtr < 999.999999
BU_ACtr: -1.0*DG_ACtr + 1000.0*BU_ACtr < 999.999999
SU_ACtr: 0 < 1.0*BU_ACtr + 1.0*FU_ACtr < 1
UF_ACtr: -1000.0*FU_ACtr + 1.0*ACtr < 0
UR_ACtr: -1000.0*BU_ACtr + 1.0*ACtr_reverse_b9231 < 0
G_ADK1: 1.0*DGo_ADK1 - 1.0*DG_ADK1 + 1.18417875325*LC_adp_c - 0.592089376625*LC_amp_c - 0.592089376625*LC_atp_c = 0
FU_ADK1: 1.0*DG_ADK1 + 1000.0*FU_ADK1 < 999.999999
BU_ADK1: -1.0*DG_ADK1 + 1000.0*BU_ADK1 < 999.999999
SU_ADK1: 0 < 1.0*BU_ADK1 + 1.0*FU_ADK1 < 1
UF_ADK1: -1000.0*FU_ADK1 + 1.0*ADK1 < 0
UR_ADK1: -1000.0*BU_ADK1 + 1.0*ADK1_reverse_a6f90 < 0
G_ADK1m: 1.0*DGo_ADK1m - 1.0*DG_ADK1m + 1.18417875325*LC_adp_m - 0.592089376625*LC_amp_m - 0.592089376625*LC_atp_m = 0
FU_ADK1m: 1.0*DG_ADK1m + 1000.0*FU_ADK1m < 999.999999
BU_ADK1m: -1.0*DG_ADK1m + 1000.0*BU_ADK1m < 999.999999
SU_ADK1m: 0 < 1.0*BU_ADK1m + 1.0*FU_ADK1m < 1
UF_ADK1m: -1000.0*FU_ADK1m + 1.0*ADK1m < 0
UR_ADK1m: -1000.0*BU_ADK1m + 1.0*ADK1m_reverse_59ef9 < 0
G_AKGDam: 1.0*DGo_AKGDam - 1.0*DG_AKGDam - 0.592089376625*LC_akg_m + 0.592089376625*LC_co2_m - 0.592089376625*LC_lpam_m + 0.592089376625*LC_sdhlam_m = 0
FU_AKGDam: 1.0*DG_AKGDam + 1000.0*FU_AKGDam < 999.999999
BU_AKGDam: -1.0*DG_AKGDam + 1000.0*BU_AKGDam < 999.999999
SU_AKGDam: 0 < 1.0*BU_AKGDam + 1.0*FU_AKGDam < 1
UF_AKGDam: -1000.0*FU_AKGDam + 1.0*AKGDam < 0
UR_AKGDam: -1000.0*BU_AKGDam + 1.0*AKGDam_reverse_a5f55 < 0
G_AKGDbm: -1.0*DG_AKGDbm - 0.592089376625*LC_coa_m + 0.592089376625*LC_dhlam_m - 0.592089376625*LC_sdhlam_m + 0.592089376625*LC_succoa_m + 1.0*DGo_AKGDbm = 0
FU_AKGDbm: 1.0*DG_AKGDbm + 1000.0*FU_AKGDbm < 999.999999
BU_AKGDbm: -1.0*DG_AKGDbm + 1000.0*BU_AKGDbm < 999.999999
SU_AKGDbm: 0 < 1.0*BU_AKGDbm + 1.0*FU_AKGDbm < 1
UF_AKGDbm: -1000.0*FU_AKGDbm + 1.0*AKGDbm < 0
UR_AKGDbm: -1000.0*BU_AKGDbm + 1.0*AKGDbm_reverse_7cf54 < 0
G_AKGMAL: 1.0*DGo_AKGMAL - 1.0*DG_AKGMAL - 0.592089376625*LC_akg_c + 0.592089376625*LC_akg_e + 0.592089376625*LC_mal-L_c - 0.592089376625*LC_mal-L_e = 0
FU_AKGMAL: 1.0*DG_AKGMAL + 1000.0*FU_AKGMAL < 999.999999
BU_AKGMAL: -1.0*DG_AKGMAL + 1000.0*BU_AKGMAL < 999.999999
SU_AKGMAL: 0 < 1.0*BU_AKGMAL + 1.0*FU_AKGMAL < 1
UF_AKGMAL: -1000.0*FU_AKGMAL + 1.0*AKGMAL < 0
UR_AKGMAL: -1000.0*BU_AKGMAL + 1.0*AKGMAL_reverse_76d29 < 0
G_AKGt2r: -1.0*DG_AKGt2r + 0.592089376625*LC_akg_c - 0.592089376625*LC_akg_e + 1.0*DGo_AKGt2r = 0
FU_AKGt2r: 1.0*DG_AKGt2r + 1000.0*FU_AKGt2r < 999.999999
BU_AKGt2r: -1.0*DG_AKGt2r + 1000.0*BU_AKGt2r < 999.999999
SU_AKGt2r: 0 < 1.0*BU_AKGt2r + 1.0*FU_AKGt2r < 1
UF_AKGt2r: -1000.0*FU_AKGt2r + 1.0*AKGt2r < 0
UR_AKGt2r: -1000.0*BU_AKGt2r + 1.0*AKGt2r_reverse_5d500 < 0
G_ALATA_L: -1.0*DG_ALATA_L - 0.592089376625*LC_akg_c - 0.592089376625*LC_ala-L_c + 0.592089376625*LC_glu-L_c + 0.592089376625*LC_pyr_c + 1.0*DGo_ALATA_L = 0
FU_ALATA_L: 1.0*DG_ALATA_L + 1000.0*FU_ALATA_L < 999.999999
BU_ALATA_L: -1.0*DG_ALATA_L + 1000.0*BU_ALATA_L < 999.999999
SU_ALATA_L: 0 < 1.0*BU_ALATA_L + 1.0*FU_ALATA_L < 1
UF_ALATA_L: -1000.0*FU_ALATA_L + 1.0*ALATA_L < 0
UR_ALATA_L: -1000.0*BU_ALATA_L + 1.0*ALATA_L_reverse_e54ff < 0
G_ALATA_Lm: -1.0*DG_ALATA_Lm - 0.592089376625*LC_akg_m - 0.592089376625*LC_ala-L_m + 0.592089376625*LC_glu-L_m + 0.592089376625*LC_pyr_m + 1.0*DGo_ALATA_Lm = 0
FU_ALATA_Lm: 1.0*DG_ALATA_Lm + 1000.0*FU_ALATA_Lm < 999.999999
BU_ALATA_Lm: -1.0*DG_ALATA_Lm + 1000.0*BU_ALATA_Lm < 999.999999
SU_ALATA_Lm: 0 < 1.0*BU_ALATA_Lm + 1.0*FU_ALATA_Lm < 1
UF_ALATA_Lm: -1000.0*FU_ALATA_Lm + 1.0*ALATA_Lm < 0
UR_ALATA_Lm: -1000.0*BU_ALATA_Lm + 1.0*ALATA_Lm_reverse_c845b < 0
G_ALAt2r: -1.0*DG_ALAt2r + 0.592089376625*LC_ala-L_c - 0.592089376625*LC_ala-L_e + 1.0*DGo_ALAt2r = 0
FU_ALAt2r: 1.0*DG_ALAt2r + 1000.0*FU_ALAt2r < 999.999999
BU_ALAt2r: -1.0*DG_ALAt2r + 1000.0*BU_ALAt2r < 999.999999
SU_ALAt2r: 0 < 1.0*BU_ALAt2r + 1.0*FU_ALAt2r < 1
UF_ALAt2r: -1000.0*FU_ALAt2r + 1.0*ALAt2r < 0
UR_ALAt2r: -1000.0*BU_ALAt2r + 1.0*ALAt2r_reverse_552b7 < 0
G_ALAtmi: -1.0*DG_ALAtmi + 0.592089376625*LC_ala-L_c - 0.592089376625*LC_ala-L_m + 1.0*DGo_ALAtmi = 0
FU_ALAtmi: 1.0*DG_ALAtmi + 1000.0*FU_ALAtmi < 999.999999
BU_ALAtmi: -1.0*DG_ALAtmi + 1000.0*BU_ALAtmi < 999.999999
SU_ALAtmi: 0 < 1.0*BU_ALAtmi + 1.0*FU_ALAtmi < 1
UF_ALAtmi: -1000.0*FU_ALAtmi + 1.0*ALAtmi < 0
UR_ALAtmi: -1000.0*BU_ALAtmi + 1.0*ALAtmi_reverse_c89ef < 0
G_ALCD2if: -1.0*DG_ALCD2if + 0.592089376625*LC_acald_c - 0.592089376625*LC_etoh_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 1.0*DGo_ALCD2if = 0
FU_ALCD2if: 1.0*DG_ALCD2if + 1000.0*FU_ALCD2if < 999.999999
BU_ALCD2if: -1.0*DG_ALCD2if + 1000.0*BU_ALCD2if < 999.999999
SU_ALCD2if: 0 < 1.0*BU_ALCD2if + 1.0*FU_ALCD2if < 1
UF_ALCD2if: -1000.0*FU_ALCD2if + 1.0*ALCD2if < 0
UR_ALCD2if: -1000.0*BU_ALCD2if + 1.0*ALCD2if_reverse_85957 < 0
G_ALCD2ir: -1.0*DG_ALCD2ir - 0.592089376625*LC_acald_c + 0.592089376625*LC_etoh_c + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadh_c + 1.0*DGo_ALCD2ir = 0
FU_ALCD2ir: 1.0*DG_ALCD2ir + 1000.0*FU_ALCD2ir < 999.999999
BU_ALCD2ir: -1.0*DG_ALCD2ir + 1000.0*BU_ALCD2ir < 999.999999
SU_ALCD2ir: 0 < 1.0*BU_ALCD2ir + 1.0*FU_ALCD2ir < 1
UF_ALCD2ir: -1000.0*FU_ALCD2ir + 1.0*ALCD2ir < 0
UR_ALCD2ir: -1000.0*BU_ALCD2ir + 1.0*ALCD2ir_reverse_ba067 < 0
G_ALCD2irm: -1.0*DG_ALCD2irm - 0.592089376625*LC_acald_m + 0.592089376625*LC_etoh_m + 0.592089376625*LC_nad_m - 0.592089376625*LC_nadh_m + 1.0*DGo_ALCD2irm = 0
FU_ALCD2irm: 1.0*DG_ALCD2irm + 1000.0*FU_ALCD2irm < 999.999999
BU_ALCD2irm: -1.0*DG_ALCD2irm + 1000.0*BU_ALCD2irm < 999.999999
SU_ALCD2irm: 0 < 1.0*BU_ALCD2irm + 1.0*FU_ALCD2irm < 1
UF_ALCD2irm: -1000.0*FU_ALCD2irm + 1.0*ALCD2irm < 0
UR_ALCD2irm: -1000.0*BU_ALCD2irm + 1.0*ALCD2irm_reverse_2557e < 0
G_ALCD2x: -1.0*DG_ALCD2x + 0.592089376625*LC_acald_c - 0.592089376625*LC_etoh_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 1.0*DGo_ALCD2x = 0
FU_ALCD2x: 1.0*DG_ALCD2x + 1000.0*FU_ALCD2x < 999.999999
BU_ALCD2x: -1.0*DG_ALCD2x + 1000.0*BU_ALCD2x < 999.999999
SU_ALCD2x: 0 < 1.0*BU_ALCD2x + 1.0*FU_ALCD2x < 1
UF_ALCD2x: -1000.0*FU_ALCD2x + 1.0*ALCD2x < 0
UR_ALCD2x: -1000.0*BU_ALCD2x + 1.0*ALCD2x_reverse_5d107 < 0
G_ALDD2xm: -1.0*DG_ALDD2xm + 1.0*DGo_ALDD2xm + 0.592089376625*LC_ac_m - 0.592089376625*LC_acald_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m = 0
FU_ALDD2xm: 1.0*DG_ALDD2xm + 1000.0*FU_ALDD2xm < 999.999999
BU_ALDD2xm: -1.0*DG_ALDD2xm + 1000.0*BU_ALDD2xm < 999.999999
SU_ALDD2xm: 0 < 1.0*BU_ALDD2xm + 1.0*FU_ALDD2xm < 1
UF_ALDD2xm: -1000.0*FU_ALDD2xm + 1.0*ALDD2xm < 0
UR_ALDD2xm: -1000.0*BU_ALDD2xm + 1.0*ALDD2xm_reverse_c93a0 < 0
G_ALDD2y: -1.0*DG_ALDD2y + 1.0*DGo_ALDD2y + 0.592089376625*LC_ac_c - 0.592089376625*LC_acald_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c = 0
FU_ALDD2y: 1.0*DG_ALDD2y + 1000.0*FU_ALDD2y < 999.999999
BU_ALDD2y: -1.0*DG_ALDD2y + 1000.0*BU_ALDD2y < 999.999999
SU_ALDD2y: 0 < 1.0*BU_ALDD2y + 1.0*FU_ALDD2y < 1
UF_ALDD2y: -1000.0*FU_ALDD2y + 1.0*ALDD2y < 0
UR_ALDD2y: -1000.0*BU_ALDD2y + 1.0*ALDD2y_reverse_03afb < 0
G_ALDD2ym: -1.0*DG_ALDD2ym + 0.592089376625*LC_ac_m - 0.592089376625*LC_acald_m - 0.592089376625*LC_nadp_m + 0.592089376625*LC_nadph_m + 1.0*DGo_ALDD2ym = 0
FU_ALDD2ym: 1.0*DG_ALDD2ym + 1000.0*FU_ALDD2ym < 999.999999
BU_ALDD2ym: -1.0*DG_ALDD2ym + 1000.0*BU_ALDD2ym < 999.999999
SU_ALDD2ym: 0 < 1.0*BU_ALDD2ym + 1.0*FU_ALDD2ym < 1
UF_ALDD2ym: -1000.0*FU_ALDD2ym + 1.0*ALDD2ym < 0
UR_ALDD2ym: -1000.0*BU_ALDD2ym + 1.0*ALDD2ym_reverse_cfc40 < 0
G_ARGt2r: -1.0*DG_ARGt2r + 0.592089376625*LC_arg-L_c - 0.592089376625*LC_arg-L_e + 1.0*DGo_ARGt2r = 0
FU_ARGt2r: 1.0*DG_ARGt2r + 1000.0*FU_ARGt2r < 999.999999
BU_ARGt2r: -1.0*DG_ARGt2r + 1000.0*BU_ARGt2r < 999.999999
SU_ARGt2r: 0 < 1.0*BU_ARGt2r + 1.0*FU_ARGt2r < 1
UF_ARGt2r: -1000.0*FU_ARGt2r + 1.0*ARGt2r < 0
UR_ARGt2r: -1000.0*BU_ARGt2r + 1.0*ARGt2r_reverse_57e3f < 0
G_ASNt2r: -1.0*DG_ASNt2r + 0.592089376625*LC_asn-L_c - 0.592089376625*LC_asn-L_e + 1.0*DGo_ASNt2r = 0
FU_ASNt2r: 1.0*DG_ASNt2r + 1000.0*FU_ASNt2r < 999.999999
BU_ASNt2r: -1.0*DG_ASNt2r + 1000.0*BU_ASNt2r < 999.999999
SU_ASNt2r: 0 < 1.0*BU_ASNt2r + 1.0*FU_ASNt2r < 1
UF_ASNt2r: -1000.0*FU_ASNt2r + 1.0*ASNt2r < 0
UR_ASNt2r: -1000.0*BU_ASNt2r + 1.0*ASNt2r_reverse_5bcd2 < 0
G_ASPGLU2m: -1.0*DG_ASPGLU2m + 0.592089376625*LC_asp-L_c - 0.592089376625*LC_asp-L_m - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_glu-L_m + 1.0*DGo_ASPGLU2m = 0
FU_ASPGLU2m: 1.0*DG_ASPGLU2m + 1000.0*FU_ASPGLU2m < 999.999999
BU_ASPGLU2m: -1.0*DG_ASPGLU2m + 1000.0*BU_ASPGLU2m < 999.999999
SU_ASPGLU2m: 0 < 1.0*BU_ASPGLU2m + 1.0*FU_ASPGLU2m < 1
UF_ASPGLU2m: -1000.0*FU_ASPGLU2m + 1.0*ASPGLU2m < 0
UR_ASPGLU2m: -1000.0*BU_ASPGLU2m + 1.0*ASPGLU2m_reverse_efaea < 0
G_ASPTA: -1.0*DG_ASPTA - 0.592089376625*LC_akg_c - 0.592089376625*LC_asp-L_c + 0.592089376625*LC_glu-L_c + 0.592089376625*LC_oaa_c + 1.0*DGo_ASPTA = 0
FU_ASPTA: 1.0*DG_ASPTA + 1000.0*FU_ASPTA < 999.999999
BU_ASPTA: -1.0*DG_ASPTA + 1000.0*BU_ASPTA < 999.999999
SU_ASPTA: 0 < 1.0*BU_ASPTA + 1.0*FU_ASPTA < 1
UF_ASPTA: -1000.0*FU_ASPTA + 1.0*ASPTA < 0
UR_ASPTA: -1000.0*BU_ASPTA + 1.0*ASPTA_reverse_36525 < 0
G_ASPTAm: -1.0*DG_ASPTAm - 0.592089376625*LC_akg_m - 0.592089376625*LC_asp-L_m + 0.592089376625*LC_glu-L_m + 0.592089376625*LC_oaa_m + 1.0*DGo_ASPTAm = 0
FU_ASPTAm: 1.0*DG_ASPTAm + 1000.0*FU_ASPTAm < 999.999999
BU_ASPTAm: -1.0*DG_ASPTAm + 1000.0*BU_ASPTAm < 999.999999
SU_ASPTAm: 0 < 1.0*BU_ASPTAm + 1.0*FU_ASPTAm < 1
UF_ASPTAm: -1000.0*FU_ASPTAm + 1.0*ASPTAm < 0
UR_ASPTAm: -1000.0*BU_ASPTAm + 1.0*ASPTAm_reverse_dc3b5 < 0
G_ASPt2m: -1.0*DG_ASPt2m - 0.592089376625*LC_asp-L_c + 0.592089376625*LC_asp-L_m + 1.0*DGo_ASPt2m = 0
FU_ASPt2m: 1.0*DG_ASPt2m + 1000.0*FU_ASPt2m < 999.999999
BU_ASPt2m: -1.0*DG_ASPt2m + 1000.0*BU_ASPt2m < 999.999999
SU_ASPt2m: 0 < 1.0*BU_ASPt2m + 1.0*FU_ASPt2m < 1
UF_ASPt2m: -1000.0*FU_ASPt2m + 1.0*ASPt2m < 0
UR_ASPt2m: -1000.0*BU_ASPt2m + 1.0*ASPt2m_reverse_c9585 < 0
G_ASPt2r: -1.0*DG_ASPt2r + 0.592089376625*LC_asp-L_c - 0.592089376625*LC_asp-L_e + 1.0*DGo_ASPt2r = 0
FU_ASPt2r: 1.0*DG_ASPt2r + 1000.0*FU_ASPt2r < 999.999999
BU_ASPt2r: -1.0*DG_ASPt2r + 1000.0*BU_ASPt2r < 999.999999
SU_ASPt2r: 0 < 1.0*BU_ASPt2r + 1.0*FU_ASPt2r < 1
UF_ASPt2r: -1000.0*FU_ASPt2r + 1.0*ASPt2r < 0
UR_ASPt2r: -1000.0*BU_ASPt2r + 1.0*ASPt2r_reverse_3a681 < 0
G_ATPH1: -1.0*DG_ATPH1 + 0.592089376625*LC_amp_c - 0.592089376625*LC_atp_c + 1.18417875325*LC_pi_c + 1.0*DGo_ATPH1 = 0
FU_ATPH1: 1.0*DG_ATPH1 + 1000.0*FU_ATPH1 < 999.999999
BU_ATPH1: -1.0*DG_ATPH1 + 1000.0*BU_ATPH1 < 999.999999
SU_ATPH1: 0 < 1.0*BU_ATPH1 + 1.0*FU_ATPH1 < 1
UF_ATPH1: -1000.0*FU_ATPH1 + 1.0*ATPH1 < 0
UR_ATPH1: -1000.0*BU_ATPH1 + 1.0*ATPH1_reverse_23c0c < 0
G_ATPM: -1.0*DG_ATPM + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_pi_c + 1.0*DGo_ATPM = 0
FU_ATPM: 1.0*DG_ATPM + 1000.0*FU_ATPM < 999.999999
BU_ATPM: -1.0*DG_ATPM + 1000.0*BU_ATPM < 999.999999
SU_ATPM: 0 < 1.0*BU_ATPM + 1.0*FU_ATPM < 1
UF_ATPM: -1000.0*FU_ATPM + 1.0*ATPM < 0
UR_ATPM: -1000.0*BU_ATPM + 1.0*ATPM_reverse_5b752 < 0
G_ATPS: -1.0*DG_ATPS + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_pi_c + 1.0*DGo_ATPS = 0
FU_ATPS: 1.0*DG_ATPS + 1000.0*FU_ATPS < 999.999999
BU_ATPS: -1.0*DG_ATPS + 1000.0*BU_ATPS < 999.999999
SU_ATPS: 0 < 1.0*BU_ATPS + 1.0*FU_ATPS < 1
UF_ATPS: -1000.0*FU_ATPS + 1.0*ATPS < 0
UR_ATPS: -1000.0*BU_ATPS + 1.0*ATPS_reverse_64286 < 0
G_ATPS3m: -1.0*DG_ATPS3m - 0.592089376625*LC_adp_m + 0.592089376625*LC_atp_m - 0.592089376625*LC_pi_m + 1.0*DGo_ATPS3m = 0
FU_ATPS3m: 1.0*DG_ATPS3m + 1000.0*FU_ATPS3m < 999.999999
BU_ATPS3m: -1.0*DG_ATPS3m + 1000.0*BU_ATPS3m < 999.999999
SU_ATPS3m: 0 < 1.0*BU_ATPS3m + 1.0*FU_ATPS3m < 1
UF_ATPS3m: -1000.0*FU_ATPS3m + 1.0*ATPS3m < 0
UR_ATPS3m: -1000.0*BU_ATPS3m + 1.0*ATPS3m_reverse_5d8d3 < 0
G_ATPtm-H: -1.0*DG_ATPtm-H + 1.0*DGo_ATPtm-H - 0.592089376625*LC_adp_c + 0.592089376625*LC_adp_m + 0.592089376625*LC_atp_c - 0.592089376625*LC_atp_m = 0
FU_ATPtm-H: 1.0*DG_ATPtm-H + 1000.0*FU_ATPtm-H < 999.999999
BU_ATPtm-H: -1.0*DG_ATPtm-H + 1000.0*BU_ATPtm-H < 999.999999
SU_ATPtm-H: 0 < 1.0*BU_ATPtm-H + 1.0*FU_ATPtm-H < 1
UF_ATPtm-H: -1000.0*FU_ATPtm-H + 1.0*ATPtm-H < 0
UR_ATPtm-H: -1000.0*BU_ATPtm-H + 1.0*ATPtm-H_reverse_ab831 < 0
G_CITt2r: 1.0*DGo_CITt2r - 1.0*DG_CITt2r + 0.592089376625*LC_cit_c - 0.592089376625*LC_cit_e = 0
FU_CITt2r: 1.0*DG_CITt2r + 1000.0*FU_CITt2r < 999.999999
BU_CITt2r: -1.0*DG_CITt2r + 1000.0*BU_CITt2r < 999.999999
SU_CITt2r: 0 < 1.0*BU_CITt2r + 1.0*FU_CITt2r < 1
UF_CITt2r: -1000.0*FU_CITt2r + 1.0*CITt2r < 0
UR_CITt2r: -1000.0*BU_CITt2r + 1.0*CITt2r_reverse_9f606 < 0
G_CITtam: -1.0*DG_CITtam - 0.592089376625*LC_cit_c + 0.592089376625*LC_cit_m + 0.592089376625*LC_mal-L_c - 0.592089376625*LC_mal-L_m + 1.0*DGo_CITtam = 0
FU_CITtam: 1.0*DG_CITtam + 1000.0*FU_CITtam < 999.999999
BU_CITtam: -1.0*DG_CITtam + 1000.0*BU_CITtam < 999.999999
SU_CITtam: 0 < 1.0*BU_CITtam + 1.0*FU_CITtam < 1
UF_CITtam: -1000.0*FU_CITtam + 1.0*CITtam < 0
UR_CITtam: -1000.0*BU_CITtam + 1.0*CITtam_reverse_5230c < 0
G_CITtcm: -1.0*DG_CITtcm - 0.592089376625*LC_cit_c + 0.592089376625*LC_cit_m + 0.592089376625*LC_icit_c - 0.592089376625*LC_icit_m + 1.0*DGo_CITtcm = 0
FU_CITtcm: 1.0*DG_CITtcm + 1000.0*FU_CITtcm < 999.999999
BU_CITtcm: -1.0*DG_CITtcm + 1000.0*BU_CITtcm < 999.999999
SU_CITtcm: 0 < 1.0*BU_CITtcm + 1.0*FU_CITtcm < 1
UF_CITtcm: -1000.0*FU_CITtcm + 1.0*CITtcm < 0
UR_CITtcm: -1000.0*BU_CITtcm + 1.0*CITtcm_reverse_df725 < 0
G_CO2t: -1.0*DG_CO2t + 0.592089376625*LC_co2_c - 0.592089376625*LC_co2_e + 1.0*DGo_CO2t = 0
FU_CO2t: 1.0*DG_CO2t + 1000.0*FU_CO2t < 999.999999
BU_CO2t: -1.0*DG_CO2t + 1000.0*BU_CO2t < 999.999999
SU_CO2t: 0 < 1.0*BU_CO2t + 1.0*FU_CO2t < 1
UF_CO2t: -1000.0*FU_CO2t + 1.0*CO2t < 0
UR_CO2t: -1000.0*BU_CO2t + 1.0*CO2t_reverse_7c42f < 0
G_CO2tm: -1.0*DG_CO2tm - 0.592089376625*LC_co2_c + 0.592089376625*LC_co2_m + 1.0*DGo_CO2tm = 0
FU_CO2tm: 1.0*DG_CO2tm + 1000.0*FU_CO2tm < 999.999999
BU_CO2tm: -1.0*DG_CO2tm + 1000.0*BU_CO2tm < 999.999999
SU_CO2tm: 0 < 1.0*BU_CO2tm + 1.0*FU_CO2tm < 1
UF_CO2tm: -1000.0*FU_CO2tm + 1.0*CO2tm < 0
UR_CO2tm: -1000.0*BU_CO2tm + 1.0*CO2tm_reverse_48c28 < 0
G_CSm: 1.0*DGo_CSm - 1.0*DG_CSm - 0.592089376625*LC_accoa_m + 0.592089376625*LC_cit_m + 0.592089376625*LC_coa_m - 0.592089376625*LC_oaa_m = 0
FU_CSm: 1.0*DG_CSm + 1000.0*FU_CSm < 999.999999
BU_CSm: -1.0*DG_CSm + 1000.0*BU_CSm < 999.999999
SU_CSm: 0 < 1.0*BU_CSm + 1.0*FU_CSm < 1
UF_CSm: -1000.0*FU_CSm + 1.0*CSm < 0
UR_CSm: -1000.0*BU_CSm + 1.0*CSm_reverse_bb3a1 < 0
G_CYOOm: -1.0*DG_CYOOm + 2.3683575065*LC_ficytc_m - 2.3683575065*LC_focytc_m - 0.592089376625*LC_o2_m + 1.0*DGo_CYOOm = 0
FU_CYOOm: 1.0*DG_CYOOm + 1000.0*FU_CYOOm < 999.999999
BU_CYOOm: -1.0*DG_CYOOm + 1000.0*BU_CYOOm < 999.999999
SU_CYOOm: 0 < 1.0*BU_CYOOm + 1.0*FU_CYOOm < 1
UF_CYOOm: -1000.0*FU_CYOOm + 1.0*CYOOm < 0
UR_CYOOm: -1000.0*BU_CYOOm + 1.0*CYOOm_reverse_d730e < 0
G_CYOR_u6m: -1.0*DG_CYOR_u6m - 1.18417875325*LC_ficytc_m + 1.18417875325*LC_focytc_m + 0.592089376625*LC_q6_m - 0.592089376625*LC_q6h2_m + 1.0*DGo_CYOR_u6m = 0
FU_CYOR_u6m: 1.0*DG_CYOR_u6m + 1000.0*FU_CYOR_u6m < 999.999999
BU_CYOR_u6m: -1.0*DG_CYOR_u6m + 1000.0*BU_CYOR_u6m < 999.999999
SU_CYOR_u6m: 0 < 1.0*BU_CYOR_u6m + 1.0*FU_CYOR_u6m < 1
UF_CYOR_u6m: -1000.0*FU_CYOR_u6m + 1.0*CYOR_u6m < 0
UR_CYOR_u6m: -1000.0*BU_CYOR_u6m + 1.0*CYOR_u6m_reverse_ba2b8 < 0
G_CYSt2r: -1.0*DG_CYSt2r + 0.592089376625*LC_cys-L_c - 0.592089376625*LC_cys-L_e + 1.0*DGo_CYSt2r = 0
FU_CYSt2r: 1.0*DG_CYSt2r + 1000.0*FU_CYSt2r < 999.999999
BU_CYSt2r: -1.0*DG_CYSt2r + 1000.0*BU_CYSt2r < 999.999999
SU_CYSt2r: 0 < 1.0*BU_CYSt2r + 1.0*FU_CYSt2r < 1
UF_CYSt2r: -1000.0*FU_CYSt2r + 1.0*CYSt2r < 0
UR_CYSt2r: -1000.0*BU_CYSt2r + 1.0*CYSt2r_reverse_6167b < 0
G_D-LACt2m: 1.0*DGo_D-LACt2m - 1.0*DG_D-LACt2m - 0.592089376625*LC_lac-D_c + 0.592089376625*LC_lac-D_m + 0.592089376625*LC_pyr_c - 0.592089376625*LC_pyr_m = 0
FU_D-LACt2m: 1.0*DG_D-LACt2m + 1000.0*FU_D-LACt2m < 999.999999
BU_D-LACt2m: -1.0*DG_D-LACt2m + 1000.0*BU_D-LACt2m < 999.999999
SU_D-LACt2m: 0 < 1.0*BU_D-LACt2m + 1.0*FU_D-LACt2m < 1
UF_D-LACt2m: -1000.0*FU_D-LACt2m + 1.0*D-LACt2m < 0
UR_D-LACt2m: -1000.0*BU_D-LACt2m + 1.0*D-LACt2m_reverse_bb9a8 < 0
G_D-LACtm: -1.0*DG_D-LACtm - 0.592089376625*LC_lac-D_c + 0.592089376625*LC_lac-D_m + 1.0*DGo_D-LACtm = 0
FU_D-LACtm: 1.0*DG_D-LACtm + 1000.0*FU_D-LACtm < 999.999999
BU_D-LACtm: -1.0*DG_D-LACtm + 1000.0*BU_D-LACtm < 999.999999
SU_D-LACtm: 0 < 1.0*BU_D-LACtm + 1.0*FU_D-LACtm < 1
UF_D-LACtm: -1000.0*FU_D-LACtm + 1.0*D-LACtm < 0
UR_D-LACtm: -1000.0*BU_D-LACtm + 1.0*D-LACtm_reverse_08640 < 0
G_DB4PS: -1.0*DG_DB4PS + 0.592089376625*LC_db4p_c + 0.592089376625*LC_for_c - 0.592089376625*LC_ru5p-D_c + 1.0*DGo_DB4PS = 0
FU_DB4PS: 1.0*DG_DB4PS + 1000.0*FU_DB4PS < 999.999999
BU_DB4PS: -1.0*DG_DB4PS + 1000.0*BU_DB4PS < 999.999999
SU_DB4PS: 0 < 1.0*BU_DB4PS + 1.0*FU_DB4PS < 1
UF_DB4PS: -1000.0*FU_DB4PS + 1.0*DB4PS < 0
UR_DB4PS: -1000.0*BU_DB4PS + 1.0*DB4PS_reverse_43dd1 < 0
G_ENO: -1.0*DG_ENO - 0.592089376625*LC_2pg_c + 0.592089376625*LC_pep_c + 1.0*DGo_ENO = 0
FU_ENO: 1.0*DG_ENO + 1000.0*FU_ENO < 999.999999
BU_ENO: -1.0*DG_ENO + 1000.0*BU_ENO < 999.999999
SU_ENO: 0 < 1.0*BU_ENO + 1.0*FU_ENO < 1
UF_ENO: -1000.0*FU_ENO + 1.0*ENO < 0
UR_ENO: -1000.0*BU_ENO + 1.0*ENO_reverse_40eea < 0
G_ERGSTt: -1.0*DG_ERGSTt + 0.592089376625*LC_ergst_c - 0.592089376625*LC_ergst_e + 1.0*DGo_ERGSTt = 0
FU_ERGSTt: 1.0*DG_ERGSTt + 1000.0*FU_ERGSTt < 999.999999
BU_ERGSTt: -1.0*DG_ERGSTt + 1000.0*BU_ERGSTt < 999.999999
SU_ERGSTt: 0 < 1.0*BU_ERGSTt + 1.0*FU_ERGSTt < 1
UF_ERGSTt: -1000.0*FU_ERGSTt + 1.0*ERGSTt < 0
UR_ERGSTt: -1000.0*BU_ERGSTt + 1.0*ERGSTt_reverse_c595b < 0
G_ERGSTter: -1.0*DG_ERGSTter + 0.592089376625*LC_ergst_c - 0.592089376625*LC_ergst_r + 1.0*DGo_ERGSTter = 0
FU_ERGSTter: 1.0*DG_ERGSTter + 1000.0*FU_ERGSTter < 999.999999
BU_ERGSTter: -1.0*DG_ERGSTter + 1000.0*BU_ERGSTter < 999.999999
SU_ERGSTter: 0 < 1.0*BU_ERGSTter + 1.0*FU_ERGSTter < 1
UF_ERGSTter: -1000.0*FU_ERGSTter + 1.0*ERGSTter < 0
UR_ERGSTter: -1000.0*BU_ERGSTter + 1.0*ERGSTter_reverse_23bab < 0
G_ETOHt: 1.0*DGo_ETOHt - 1.0*DG_ETOHt + 0.592089376625*LC_etoh_c - 0.592089376625*LC_etoh_e = 0
FU_ETOHt: 1.0*DG_ETOHt + 1000.0*FU_ETOHt < 999.999999
BU_ETOHt: -1.0*DG_ETOHt + 1000.0*BU_ETOHt < 999.999999
SU_ETOHt: 0 < 1.0*BU_ETOHt + 1.0*FU_ETOHt < 1
UF_ETOHt: -1000.0*FU_ETOHt + 1.0*ETOHt < 0
UR_ETOHt: -1000.0*BU_ETOHt + 1.0*ETOHt_reverse_a82a5 < 0
G_ETOHtm: -1.0*DG_ETOHtm - 0.592089376625*LC_etoh_c + 0.592089376625*LC_etoh_m + 1.0*DGo_ETOHtm = 0
FU_ETOHtm: 1.0*DG_ETOHtm + 1000.0*FU_ETOHtm < 999.999999
BU_ETOHtm: -1.0*DG_ETOHtm + 1000.0*BU_ETOHtm < 999.999999
SU_ETOHtm: 0 < 1.0*BU_ETOHtm + 1.0*FU_ETOHtm < 1
UF_ETOHtm: -1000.0*FU_ETOHtm + 1.0*ETOHtm < 0
UR_ETOHtm: -1000.0*BU_ETOHtm + 1.0*ETOHtm_reverse_9d4eb < 0
G_FBA: -1.0*DG_FBA + 0.592089376625*LC_dhap_c - 0.592089376625*LC_fdp_c + 0.592089376625*LC_g3p_c + 1.0*DGo_FBA = 0
FU_FBA: 1.0*DG_FBA + 1000.0*FU_FBA < 999.999999
BU_FBA: -1.0*DG_FBA + 1000.0*BU_FBA < 999.999999
SU_FBA: 0 < 1.0*BU_FBA + 1.0*FU_FBA < 1
UF_FBA: -1000.0*FU_FBA + 1.0*FBA < 0
UR_FBA: -1000.0*BU_FBA + 1.0*FBA_reverse_84806 < 0
G_FBA3: -1.0*DG_FBA3 + 0.592089376625*LC_dhap_c + 0.592089376625*LC_e4p_c - 0.592089376625*LC_s17bp_c + 1.0*DGo_FBA3 = 0
FU_FBA3: 1.0*DG_FBA3 + 1000.0*FU_FBA3 < 999.999999
BU_FBA3: -1.0*DG_FBA3 + 1000.0*BU_FBA3 < 999.999999
SU_FBA3: 0 < 1.0*BU_FBA3 + 1.0*FU_FBA3 < 1
UF_FBA3: -1000.0*FU_FBA3 + 1.0*FBA3 < 0
UR_FBA3: -1000.0*BU_FBA3 + 1.0*FBA3_reverse_0d49f < 0
G_FBP: 1.0*DGo_FBP - 1.0*DG_FBP + 0.592089376625*LC_f6p_c - 0.592089376625*LC_fdp_c + 0.592089376625*LC_pi_c = 0
FU_FBP: 1.0*DG_FBP + 1000.0*FU_FBP < 999.999999
BU_FBP: -1.0*DG_FBP + 1000.0*BU_FBP < 999.999999
SU_FBP: 0 < 1.0*BU_FBP + 1.0*FU_FBP < 1
UF_FBP: -1000.0*FU_FBP + 1.0*FBP < 0
UR_FBP: -1000.0*BU_FBP + 1.0*FBP_reverse_bf2c9 < 0
G_FDH: -1.0*DG_FDH + 0.592089376625*LC_co2_c - 0.592089376625*LC_for_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 1.0*DGo_FDH = 0
FU_FDH: 1.0*DG_FDH + 1000.0*FU_FDH < 999.999999
BU_FDH: -1.0*DG_FDH + 1000.0*BU_FDH < 999.999999
SU_FDH: 0 < 1.0*BU_FDH + 1.0*FU_FDH < 1
UF_FDH: -1000.0*FU_FDH + 1.0*FDH < 0
UR_FDH: -1000.0*BU_FDH + 1.0*FDH_reverse_06346 < 0
G_FDNG: -1.0*DG_FDNG + 0.592089376625*LC_co2_c - 0.592089376625*LC_for_c - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m + 1.0*DGo_FDNG = 0
FU_FDNG: 1.0*DG_FDNG + 1000.0*FU_FDNG < 999.999999
BU_FDNG: -1.0*DG_FDNG + 1000.0*BU_FDNG < 999.999999
SU_FDNG: 0 < 1.0*BU_FDNG + 1.0*FU_FDNG < 1
UF_FDNG: -1000.0*FU_FDNG + 1.0*FDNG < 0
UR_FDNG: -1000.0*BU_FDNG + 1.0*FDNG_reverse_f3615 < 0
G_FORt: 1.0*DGo_FORt - 1.0*DG_FORt + 0.592089376625*LC_for_c - 0.592089376625*LC_for_e = 0
FU_FORt: 1.0*DG_FORt + 1000.0*FU_FORt < 999.999999
BU_FORt: -1.0*DG_FORt + 1000.0*BU_FORt < 999.999999
SU_FORt: 0 < 1.0*BU_FORt + 1.0*FU_FORt < 1
UF_FORt: -1000.0*FU_FORt + 1.0*FORt < 0
UR_FORt: -1000.0*BU_FORt + 1.0*FORt_reverse_40f9f < 0
G_FRDcm: -1.0*DG_FRDcm + 0.592089376625*LC_fad_m - 0.592089376625*LC_fadh2_m - 0.592089376625*LC_fum_c + 0.592089376625*LC_succ_c + 1.0*DGo_FRDcm = 0
FU_FRDcm: 1.0*DG_FRDcm + 1000.0*FU_FRDcm < 999.999999
BU_FRDcm: -1.0*DG_FRDcm + 1000.0*BU_FRDcm < 999.999999
SU_FRDcm: 0 < 1.0*BU_FRDcm + 1.0*FU_FRDcm < 1
UF_FRDcm: -1000.0*FU_FRDcm + 1.0*FRDcm < 0
UR_FRDcm: -1000.0*BU_FRDcm + 1.0*FRDcm_reverse_6cbef < 0
G_FRDm: -1.0*DG_FRDm + 1.0*DGo_FRDm + 0.592089376625*LC_fad_m - 0.592089376625*LC_fadh2_m - 0.592089376625*LC_fum_m + 0.592089376625*LC_succ_m = 0
FU_FRDm: 1.0*DG_FRDm + 1000.0*FU_FRDm < 999.999999
BU_FRDm: -1.0*DG_FRDm + 1000.0*BU_FRDm < 999.999999
SU_FRDm: 0 < 1.0*BU_FRDm + 1.0*FU_FRDm < 1
UF_FRDm: -1000.0*FU_FRDm + 1.0*FRDm < 0
UR_FRDm: -1000.0*BU_FRDm + 1.0*FRDm_reverse_ef6ae < 0
G_FUM: 1.0*DGo_FUM - 1.0*DG_FUM - 0.592089376625*LC_fum_c + 0.592089376625*LC_mal-L_c = 0
FU_FUM: 1.0*DG_FUM + 1000.0*FU_FUM < 999.999999
BU_FUM: -1.0*DG_FUM + 1000.0*BU_FUM < 999.999999
SU_FUM: 0 < 1.0*BU_FUM + 1.0*FU_FUM < 1
UF_FUM: -1000.0*FU_FUM + 1.0*FUM < 0
UR_FUM: -1000.0*BU_FUM + 1.0*FUM_reverse_d3642 < 0
G_FUMm: 1.0*DGo_FUMm - 1.0*DG_FUMm - 0.592089376625*LC_fum_m + 0.592089376625*LC_mal-L_m = 0
FU_FUMm: 1.0*DG_FUMm + 1000.0*FU_FUMm < 999.999999
BU_FUMm: -1.0*DG_FUMm + 1000.0*BU_FUMm < 999.999999
SU_FUMm: 0 < 1.0*BU_FUMm + 1.0*FU_FUMm < 1
UF_FUMm: -1000.0*FU_FUMm + 1.0*FUMm < 0
UR_FUMm: -1000.0*BU_FUMm + 1.0*FUMm_reverse_5e620 < 0
G_FUMt2r: -1.0*DG_FUMt2r + 0.592089376625*LC_fum_c - 0.592089376625*LC_fum_e + 1.0*DGo_FUMt2r = 0
FU_FUMt2r: 1.0*DG_FUMt2r + 1000.0*FU_FUMt2r < 999.999999
BU_FUMt2r: -1.0*DG_FUMt2r + 1000.0*BU_FUMt2r < 999.999999
SU_FUMt2r: 0 < 1.0*BU_FUMt2r + 1.0*FU_FUMt2r < 1
UF_FUMt2r: -1000.0*FU_FUMt2r + 1.0*FUMt2r < 0
UR_FUMt2r: -1000.0*BU_FUMt2r + 1.0*FUMt2r_reverse_d93c0 < 0
G_G6PDH2: -1.0*DG_G6PDH2 + 0.592089376625*LC_6pgl_c - 0.592089376625*LC_g6p_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c + 1.0*DGo_G6PDH2 = 0
FU_G6PDH2: 1.0*DG_G6PDH2 + 1000.0*FU_G6PDH2 < 999.999999
BU_G6PDH2: -1.0*DG_G6PDH2 + 1000.0*BU_G6PDH2 < 999.999999
SU_G6PDH2: 0 < 1.0*BU_G6PDH2 + 1.0*FU_G6PDH2 < 1
UF_G6PDH2: -1000.0*FU_G6PDH2 + 1.0*G6PDH2 < 0
UR_G6PDH2: -1000.0*BU_G6PDH2 + 1.0*G6PDH2_reverse_47bd4 < 0
G_G6PDH2er: -1.0*DG_G6PDH2er + 0.592089376625*LC_6pgl_r - 0.592089376625*LC_g6p_r - 0.592089376625*LC_nadp_r + 0.592089376625*LC_nadph_r + 1.0*DGo_G6PDH2er = 0
FU_G6PDH2er: 1.0*DG_G6PDH2er + 1000.0*FU_G6PDH2er < 999.999999
BU_G6PDH2er: -1.0*DG_G6PDH2er + 1000.0*BU_G6PDH2er < 999.999999
SU_G6PDH2er: 0 < 1.0*BU_G6PDH2er + 1.0*FU_G6PDH2er < 1
UF_G6PDH2er: -1000.0*FU_G6PDH2er + 1.0*G6PDH2er < 0
UR_G6PDH2er: -1000.0*BU_G6PDH2er + 1.0*G6PDH2er_reverse_9f5fd < 0
G_G6PI: 1.0*DGo_G6PI - 1.0*DG_G6PI + 0.592089376625*LC_g6p-B_c - 0.592089376625*LC_g6p_c = 0
FU_G6PI: 1.0*DG_G6PI + 1000.0*FU_G6PI < 999.999999
BU_G6PI: -1.0*DG_G6PI + 1000.0*BU_G6PI < 999.999999
SU_G6PI: 0 < 1.0*BU_G6PI + 1.0*FU_G6PI < 1
UF_G6PI: -1000.0*FU_G6PI + 1.0*G6PI < 0
UR_G6PI: -1000.0*BU_G6PI + 1.0*G6PI_reverse_834e6 < 0
G_G6PI3: -1.0*DG_G6PI3 + 0.592089376625*LC_f6p_c - 0.592089376625*LC_g6p-B_c + 1.0*DGo_G6PI3 = 0
FU_G6PI3: 1.0*DG_G6PI3 + 1000.0*FU_G6PI3 < 999.999999
BU_G6PI3: -1.0*DG_G6PI3 + 1000.0*BU_G6PI3 < 999.999999
SU_G6PI3: 0 < 1.0*BU_G6PI3 + 1.0*FU_G6PI3 < 1
UF_G6PI3: -1000.0*FU_G6PI3 + 1.0*G6PI3 < 0
UR_G6PI3: -1000.0*BU_G6PI3 + 1.0*G6PI3_reverse_734c4 < 0
G_G6Pter: -1.0*DG_G6Pter - 0.592089376625*LC_g6p_c + 0.592089376625*LC_g6p_r + 1.0*DGo_G6Pter = 0
FU_G6Pter: 1.0*DG_G6Pter + 1000.0*FU_G6Pter < 999.999999
BU_G6Pter: -1.0*DG_G6Pter + 1000.0*BU_G6Pter < 999.999999
SU_G6Pter: 0 < 1.0*BU_G6Pter + 1.0*FU_G6Pter < 1
UF_G6Pter: -1000.0*FU_G6Pter + 1.0*G6Pter < 0
UR_G6Pter: -1000.0*BU_G6Pter + 1.0*G6Pter_reverse_b4f95 < 0
G_GAPD: -1.0*DG_GAPD + 0.592089376625*LC_13dpg_c - 0.592089376625*LC_g3p_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c - 0.592089376625*LC_pi_c + 1.0*DGo_GAPD = 0
FU_GAPD: 1.0*DG_GAPD + 1000.0*FU_GAPD < 999.999999
BU_GAPD: -1.0*DG_GAPD + 1000.0*BU_GAPD < 999.999999
SU_GAPD: 0 < 1.0*BU_GAPD + 1.0*FU_GAPD < 1
UF_GAPD: -1000.0*FU_GAPD + 1.0*GAPD < 0
UR_GAPD: -1000.0*BU_GAPD + 1.0*GAPD_reverse_459c1 < 0
G_GCC2cm: -1.0*DG_GCC2cm - 0.592089376625*LC_dhlam_m + 0.592089376625*LC_lpam_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m + 1.0*DGo_GCC2cm = 0
FU_GCC2cm: 1.0*DG_GCC2cm + 1000.0*FU_GCC2cm < 999.999999
BU_GCC2cm: -1.0*DG_GCC2cm + 1000.0*BU_GCC2cm < 999.999999
SU_GCC2cm: 0 < 1.0*BU_GCC2cm + 1.0*FU_GCC2cm < 1
UF_GCC2cm: -1000.0*FU_GCC2cm + 1.0*GCC2cm < 0
UR_GCC2cm: -1000.0*BU_GCC2cm + 1.0*GCC2cm_reverse_e896f < 0
G_GLCt1: -1.0*DG_GLCt1 + 0.592089376625*LC_glc-D_c - 0.592089376625*LC_glc-D_e + 1.0*DGo_GLCt1 = 0
FU_GLCt1: 1.0*DG_GLCt1 + 1000.0*FU_GLCt1 < 999.999999
BU_GLCt1: -1.0*DG_GLCt1 + 1000.0*BU_GLCt1 < 999.999999
SU_GLCt1: 0 < 1.0*BU_GLCt1 + 1.0*FU_GLCt1 < 1
UF_GLCt1: -1000.0*FU_GLCt1 + 1.0*GLCt1 < 0
UR_GLCt1: -1000.0*BU_GLCt1 + 1.0*GLCt1_reverse_3bcc5 < 0
G_GLNt2r: -1.0*DG_GLNt2r + 0.592089376625*LC_gln-L_c - 0.592089376625*LC_gln-L_e + 1.0*DGo_GLNt2r = 0
FU_GLNt2r: 1.0*DG_GLNt2r + 1000.0*FU_GLNt2r < 999.999999
BU_GLNt2r: -1.0*DG_GLNt2r + 1000.0*BU_GLNt2r < 999.999999
SU_GLNt2r: 0 < 1.0*BU_GLNt2r + 1.0*FU_GLNt2r < 1
UF_GLNt2r: -1000.0*FU_GLNt2r + 1.0*GLNt2r < 0
UR_GLNt2r: -1000.0*BU_GLNt2r + 1.0*GLNt2r_reverse_35d46 < 0
G_GLUDC: -1.0*DG_GLUDC + 0.592089376625*LC_4abut_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_glu-L_c + 1.0*DGo_GLUDC = 0
FU_GLUDC: 1.0*DG_GLUDC + 1000.0*FU_GLUDC < 999.999999
BU_GLUDC: -1.0*DG_GLUDC + 1000.0*BU_GLUDC < 999.999999
SU_GLUDC: 0 < 1.0*BU_GLUDC + 1.0*FU_GLUDC < 1
UF_GLUDC: -1000.0*FU_GLUDC + 1.0*GLUDC < 0
UR_GLUDC: -1000.0*BU_GLUDC + 1.0*GLUDC_reverse_c4aae < 0
G_GLUDxi: -1.0*DG_GLUDxi + 0.592089376625*LC_akg_c - 0.592089376625*LC_glu-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nh4_c + 1.0*DGo_GLUDxi = 0
FU_GLUDxi: 1.0*DG_GLUDxi + 1000.0*FU_GLUDxi < 999.999999
BU_GLUDxi: -1.0*DG_GLUDxi + 1000.0*BU_GLUDxi < 999.999999
SU_GLUDxi: 0 < 1.0*BU_GLUDxi + 1.0*FU_GLUDxi < 1
UF_GLUDxi: -1000.0*FU_GLUDxi + 1.0*GLUDxi < 0
UR_GLUDxi: -1000.0*BU_GLUDxi + 1.0*GLUDxi_reverse_6fbe2 < 0
G_GLUDyi: -1.0*DG_GLUDyi + 1.0*DGo_GLUDyi - 0.592089376625*LC_akg_c + 0.592089376625*LC_glu-L_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 0.592089376625*LC_nh4_c = 0
FU_GLUDyi: 1.0*DG_GLUDyi + 1000.0*FU_GLUDyi < 999.999999
BU_GLUDyi: -1.0*DG_GLUDyi + 1000.0*BU_GLUDyi < 999.999999
SU_GLUDyi: 0 < 1.0*BU_GLUDyi + 1.0*FU_GLUDyi < 1
UF_GLUDyi: -1000.0*FU_GLUDyi + 1.0*GLUDyi < 0
UR_GLUDyi: -1000.0*BU_GLUDyi + 1.0*GLUDyi_reverse_a2351 < 0
G_GLUK: -1.0*DG_GLUK + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_g6p-B_c - 0.592089376625*LC_glc-D_c + 1.0*DGo_GLUK = 0
FU_GLUK: 1.0*DG_GLUK + 1000.0*FU_GLUK < 999.999999
BU_GLUK: -1.0*DG_GLUK + 1000.0*BU_GLUK < 999.999999
SU_GLUK: 0 < 1.0*BU_GLUK + 1.0*FU_GLUK < 1
UF_GLUK: -1000.0*FU_GLUK + 1.0*GLUK < 0
UR_GLUK: -1000.0*BU_GLUK + 1.0*GLUK_reverse_ecb52 < 0
G_GLUt2r: -1.0*DG_GLUt2r + 0.592089376625*LC_glu-L_c - 0.592089376625*LC_glu-L_e + 1.0*DGo_GLUt2r = 0
FU_GLUt2r: 1.0*DG_GLUt2r + 1000.0*FU_GLUt2r < 999.999999
BU_GLUt2r: -1.0*DG_GLUt2r + 1000.0*BU_GLUt2r < 999.999999
SU_GLUt2r: 0 < 1.0*BU_GLUt2r + 1.0*FU_GLUt2r < 1
UF_GLUt2r: -1000.0*FU_GLUt2r + 1.0*GLUt2r < 0
UR_GLUt2r: -1000.0*BU_GLUt2r + 1.0*GLUt2r_reverse_3e88a < 0
G_GLUt7m: -1.0*DG_GLUt7m - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_glu-L_m + 1.0*DGo_GLUt7m = 0
FU_GLUt7m: 1.0*DG_GLUt7m + 1000.0*FU_GLUt7m < 999.999999
BU_GLUt7m: -1.0*DG_GLUt7m + 1000.0*BU_GLUt7m < 999.999999
SU_GLUt7m: 0 < 1.0*BU_GLUt7m + 1.0*FU_GLUt7m < 1
UF_GLUt7m: -1000.0*FU_GLUt7m + 1.0*GLUt7m < 0
UR_GLUt7m: -1000.0*BU_GLUt7m + 1.0*GLUt7m_reverse_26492 < 0
G_GLYt2r: -1.0*DG_GLYt2r + 0.592089376625*LC_gly_c - 0.592089376625*LC_gly_e + 1.0*DGo_GLYt2r = 0
FU_GLYt2r: 1.0*DG_GLYt2r + 1000.0*FU_GLYt2r < 999.999999
BU_GLYt2r: -1.0*DG_GLYt2r + 1000.0*BU_GLYt2r < 999.999999
SU_GLYt2r: 0 < 1.0*BU_GLYt2r + 1.0*FU_GLYt2r < 1
UF_GLYt2r: -1000.0*FU_GLYt2r + 1.0*GLYt2r < 0
UR_GLYt2r: -1000.0*BU_GLYt2r + 1.0*GLYt2r_reverse_d768d < 0
G_GND: 1.0*DGo_GND - 1.0*DG_GND - 0.592089376625*LC_6pgc_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c + 0.592089376625*LC_ru5p-D_c = 0
FU_GND: 1.0*DG_GND + 1000.0*FU_GND < 999.999999
BU_GND: -1.0*DG_GND + 1000.0*BU_GND < 999.999999
SU_GND: 0 < 1.0*BU_GND + 1.0*FU_GND < 1
UF_GND: -1000.0*FU_GND + 1.0*GND < 0
UR_GND: -1000.0*BU_GND + 1.0*GND_reverse_eec5c < 0
SU_H2Ot: 0 < 1.0*BU_H2Ot + 1.0*FU_H2Ot < 1
UF_H2Ot: -1000.0*FU_H2Ot + 1.0*H2Ot < 0
UR_H2Ot: -1000.0*BU_H2Ot + 1.0*H2Ot_reverse_aa560 < 0
SU_H2Oter: 0 < 1.0*BU_H2Oter + 1.0*FU_H2Oter < 1
UF_H2Oter: -1000.0*FU_H2Oter + 1.0*H2Oter < 0
UR_H2Oter: -1000.0*BU_H2Oter + 1.0*H2Oter_reverse_c874d < 0
SU_H2Otm: 0 < 1.0*BU_H2Otm + 1.0*FU_H2Otm < 1
UF_H2Otm: -1000.0*FU_H2Otm + 1.0*H2Otm < 0
UR_H2Otm: -1000.0*BU_H2Otm + 1.0*H2Otm_reverse_d5f1a < 0
G_HCITSm: 1.0*DGo_HCITSm - 1.0*DG_HCITSm - 0.592089376625*LC_accoa_m - 0.592089376625*LC_akg_m + 0.592089376625*LC_coa_m + 0.592089376625*LC_hcit_m = 0
FU_HCITSm: 1.0*DG_HCITSm + 1000.0*FU_HCITSm < 999.999999
BU_HCITSm: -1.0*DG_HCITSm + 1000.0*BU_HCITSm < 999.999999
SU_HCITSm: 0 < 1.0*BU_HCITSm + 1.0*FU_HCITSm < 1
UF_HCITSm: -1000.0*FU_HCITSm + 1.0*HCITSm < 0
UR_HCITSm: -1000.0*BU_HCITSm + 1.0*HCITSm_reverse_145c0 < 0
G_HCO3E: -1.0*DG_HCO3E - 0.592089376625*LC_co2_c + 0.592089376625*LC_hco3_c + 1.0*DGo_HCO3E = 0
FU_HCO3E: 1.0*DG_HCO3E + 1000.0*FU_HCO3E < 999.999999
BU_HCO3E: -1.0*DG_HCO3E + 1000.0*BU_HCO3E < 999.999999
SU_HCO3E: 0 < 1.0*BU_HCO3E + 1.0*FU_HCO3E < 1
UF_HCO3E: -1000.0*FU_HCO3E + 1.0*HCO3E < 0
UR_HCO3E: -1000.0*BU_HCO3E + 1.0*HCO3E_reverse_97ea5 < 0
G_HEX1: -1.0*DG_HEX1 + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_g6p_c - 0.592089376625*LC_glc-D_c + 1.0*DGo_HEX1 = 0
FU_HEX1: 1.0*DG_HEX1 + 1000.0*FU_HEX1 < 999.999999
BU_HEX1: -1.0*DG_HEX1 + 1000.0*BU_HEX1 < 999.999999
SU_HEX1: 0 < 1.0*BU_HEX1 + 1.0*FU_HEX1 < 1
UF_HEX1: -1000.0*FU_HEX1 + 1.0*HEX1 < 0
UR_HEX1: -1000.0*BU_HEX1 + 1.0*HEX1_reverse_25efa < 0
G_HISt2r: -1.0*DG_HISt2r + 0.592089376625*LC_his-L_c - 0.592089376625*LC_his-L_e + 1.0*DGo_HISt2r = 0
FU_HISt2r: 1.0*DG_HISt2r + 1000.0*FU_HISt2r < 999.999999
BU_HISt2r: -1.0*DG_HISt2r + 1000.0*BU_HISt2r < 999.999999
SU_HISt2r: 0 < 1.0*BU_HISt2r + 1.0*FU_HISt2r < 1
UF_HISt2r: -1000.0*FU_HISt2r + 1.0*HISt2r < 0
UR_HISt2r: -1000.0*BU_HISt2r + 1.0*HISt2r_reverse_b00a3 < 0
G_ICDHxm: -1.0*DG_ICDHxm + 0.592089376625*LC_akg_m + 0.592089376625*LC_co2_m - 0.592089376625*LC_icit_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m + 1.0*DGo_ICDHxm = 0
FU_ICDHxm: 1.0*DG_ICDHxm + 1000.0*FU_ICDHxm < 999.999999
BU_ICDHxm: -1.0*DG_ICDHxm + 1000.0*BU_ICDHxm < 999.999999
SU_ICDHxm: 0 < 1.0*BU_ICDHxm + 1.0*FU_ICDHxm < 1
UF_ICDHxm: -1000.0*FU_ICDHxm + 1.0*ICDHxm < 0
UR_ICDHxm: -1000.0*BU_ICDHxm + 1.0*ICDHxm_reverse_14a38 < 0
G_ICDHy: -1.0*DG_ICDHy + 0.592089376625*LC_akg_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_icit_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c + 1.0*DGo_ICDHy = 0
FU_ICDHy: 1.0*DG_ICDHy + 1000.0*FU_ICDHy < 999.999999
BU_ICDHy: -1.0*DG_ICDHy + 1000.0*BU_ICDHy < 999.999999
SU_ICDHy: 0 < 1.0*BU_ICDHy + 1.0*FU_ICDHy < 1
UF_ICDHy: -1000.0*FU_ICDHy + 1.0*ICDHy < 0
UR_ICDHy: -1000.0*BU_ICDHy + 1.0*ICDHy_reverse_c12f8 < 0
G_ICDHym: -1.0*DG_ICDHym + 0.592089376625*LC_akg_m + 0.592089376625*LC_co2_m - 0.592089376625*LC_icit_m - 0.592089376625*LC_nadp_m + 0.592089376625*LC_nadph_m + 1.0*DGo_ICDHym = 0
FU_ICDHym: 1.0*DG_ICDHym + 1000.0*FU_ICDHym < 999.999999
BU_ICDHym: -1.0*DG_ICDHym + 1000.0*BU_ICDHym < 999.999999
SU_ICDHym: 0 < 1.0*BU_ICDHym + 1.0*FU_ICDHym < 1
UF_ICDHym: -1000.0*FU_ICDHym + 1.0*ICDHym < 0
UR_ICDHym: -1000.0*BU_ICDHym + 1.0*ICDHym_reverse_38fb7 < 0
G_ILEt2r: -1.0*DG_ILEt2r + 0.592089376625*LC_ile-L_c - 0.592089376625*LC_ile-L_e + 1.0*DGo_ILEt2r = 0
FU_ILEt2r: 1.0*DG_ILEt2r + 1000.0*FU_ILEt2r < 999.999999
BU_ILEt2r: -1.0*DG_ILEt2r + 1000.0*BU_ILEt2r < 999.999999
SU_ILEt2r: 0 < 1.0*BU_ILEt2r + 1.0*FU_ILEt2r < 1
UF_ILEt2r: -1000.0*FU_ILEt2r + 1.0*ILEt2r < 0
UR_ILEt2r: -1000.0*BU_ILEt2r + 1.0*ILEt2r_reverse_9891b < 0
G_ILEtmi: -1.0*DG_ILEtmi + 0.592089376625*LC_ile-L_c - 0.592089376625*LC_ile-L_m + 1.0*DGo_ILEtmi = 0
FU_ILEtmi: 1.0*DG_ILEtmi + 1000.0*FU_ILEtmi < 999.999999
BU_ILEtmi: -1.0*DG_ILEtmi + 1000.0*BU_ILEtmi < 999.999999
SU_ILEtmi: 0 < 1.0*BU_ILEtmi + 1.0*FU_ILEtmi < 1
UF_ILEtmi: -1000.0*FU_ILEtmi + 1.0*ILEtmi < 0
UR_ILEtmi: -1000.0*BU_ILEtmi + 1.0*ILEtmi_reverse_5c3e7 < 0
G_LEUt2r: -1.0*DG_LEUt2r + 0.592089376625*LC_leu-L_c - 0.592089376625*LC_leu-L_e + 1.0*DGo_LEUt2r = 0
FU_LEUt2r: 1.0*DG_LEUt2r + 1000.0*FU_LEUt2r < 999.999999
BU_LEUt2r: -1.0*DG_LEUt2r + 1000.0*BU_LEUt2r < 999.999999
SU_LEUt2r: 0 < 1.0*BU_LEUt2r + 1.0*FU_LEUt2r < 1
UF_LEUt2r: -1000.0*FU_LEUt2r + 1.0*LEUt2r < 0
UR_LEUt2r: -1000.0*BU_LEUt2r + 1.0*LEUt2r_reverse_4a477 < 0
G_LYSt2r: -1.0*DG_LYSt2r + 0.592089376625*LC_lys-L_c - 0.592089376625*LC_lys-L_e + 1.0*DGo_LYSt2r = 0
FU_LYSt2r: 1.0*DG_LYSt2r + 1000.0*FU_LYSt2r < 999.999999
BU_LYSt2r: -1.0*DG_LYSt2r + 1000.0*BU_LYSt2r < 999.999999
SU_LYSt2r: 0 < 1.0*BU_LYSt2r + 1.0*FU_LYSt2r < 1
UF_LYSt2r: -1000.0*FU_LYSt2r + 1.0*LYSt2r < 0
UR_LYSt2r: -1000.0*BU_LYSt2r + 1.0*LYSt2r_reverse_7e0ef < 0
G_MALt2r: -1.0*DG_MALt2r + 0.592089376625*LC_mal-L_c - 0.592089376625*LC_mal-L_e + 1.0*DGo_MALt2r = 0
FU_MALt2r: 1.0*DG_MALt2r + 1000.0*FU_MALt2r < 999.999999
BU_MALt2r: -1.0*DG_MALt2r + 1000.0*BU_MALt2r < 999.999999
SU_MALt2r: 0 < 1.0*BU_MALt2r + 1.0*FU_MALt2r < 1
UF_MALt2r: -1000.0*FU_MALt2r + 1.0*MALt2r < 0
UR_MALt2r: -1000.0*BU_MALt2r + 1.0*MALt2r_reverse_bfe2f < 0
G_MALtm: -1.0*DG_MALtm - 0.592089376625*LC_mal-L_c + 0.592089376625*LC_mal-L_m + 0.592089376625*LC_pi_c - 0.592089376625*LC_pi_m + 1.0*DGo_MALtm = 0
FU_MALtm: 1.0*DG_MALtm + 1000.0*FU_MALtm < 999.999999
BU_MALtm: -1.0*DG_MALtm + 1000.0*BU_MALtm < 999.999999
SU_MALtm: 0 < 1.0*BU_MALtm + 1.0*FU_MALtm < 1
UF_MALtm: -1000.0*FU_MALtm + 1.0*MALtm < 0
UR_MALtm: -1000.0*BU_MALtm + 1.0*MALtm_reverse_f798a < 0
G_MAN6PI: -1.0*DG_MAN6PI + 0.592089376625*LC_f6p_c - 0.592089376625*LC_man6p_c + 1.0*DGo_MAN6PI = 0
FU_MAN6PI: 1.0*DG_MAN6PI + 1000.0*FU_MAN6PI < 999.999999
BU_MAN6PI: -1.0*DG_MAN6PI + 1000.0*BU_MAN6PI < 999.999999
SU_MAN6PI: 0 < 1.0*BU_MAN6PI + 1.0*FU_MAN6PI < 1
UF_MAN6PI: -1000.0*FU_MAN6PI + 1.0*MAN6PI < 0
UR_MAN6PI: -1000.0*BU_MAN6PI + 1.0*MAN6PI_reverse_d96f0 < 0
G_MANNANter: -1.0*DG_MANNANter - 0.592089376625*LC_mannan_c + 0.592089376625*LC_mannan_r + 1.0*DGo_MANNANter = 0
FU_MANNANter: 1.0*DG_MANNANter + 1000.0*FU_MANNANter < 999.999999
BU_MANNANter: -1.0*DG_MANNANter + 1000.0*BU_MANNANter < 999.999999
SU_MANNANter: 0 < 1.0*BU_MANNANter + 1.0*FU_MANNANter < 1
UF_MANNANter: -1000.0*FU_MANNANter + 1.0*MANNANter < 0
UR_MANNANter: -1000.0*BU_MANNANter + 1.0*MANNANter_reverse_6b4bd < 0
G_MDH: -1.0*DG_MDH - 0.592089376625*LC_mal-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_oaa_c + 1.0*DGo_MDH = 0
FU_MDH: 1.0*DG_MDH + 1000.0*FU_MDH < 999.999999
BU_MDH: -1.0*DG_MDH + 1000.0*BU_MDH < 999.999999
SU_MDH: 0 < 1.0*BU_MDH + 1.0*FU_MDH < 1
UF_MDH: -1000.0*FU_MDH + 1.0*MDH < 0
UR_MDH: -1000.0*BU_MDH + 1.0*MDH_reverse_ee52c < 0
G_MDHm: -1.0*DG_MDHm - 0.592089376625*LC_mal-L_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m + 0.592089376625*LC_oaa_m + 1.0*DGo_MDHm = 0
FU_MDHm: 1.0*DG_MDHm + 1000.0*FU_MDHm < 999.999999
BU_MDHm: -1.0*DG_MDHm + 1000.0*BU_MDHm < 999.999999
SU_MDHm: 0 < 1.0*BU_MDHm + 1.0*FU_MDHm < 1
UF_MDHm: -1000.0*FU_MDHm + 1.0*MDHm < 0
UR_MDHm: -1000.0*BU_MDHm + 1.0*MDHm_reverse_80605 < 0
G_ME1m: -1.0*DG_ME1m + 0.592089376625*LC_co2_m - 0.592089376625*LC_mal-L_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m + 0.592089376625*LC_pyr_m + 1.0*DGo_ME1m = 0
FU_ME1m: 1.0*DG_ME1m + 1000.0*FU_ME1m < 999.999999
BU_ME1m: -1.0*DG_ME1m + 1000.0*BU_ME1m < 999.999999
SU_ME1m: 0 < 1.0*BU_ME1m + 1.0*FU_ME1m < 1
UF_ME1m: -1000.0*FU_ME1m + 1.0*ME1m < 0
UR_ME1m: -1000.0*BU_ME1m + 1.0*ME1m_reverse_6868b < 0
G_ME2m: 1.0*DGo_ME2m - 1.0*DG_ME2m + 0.592089376625*LC_co2_m - 0.592089376625*LC_mal-L_m - 0.592089376625*LC_nadp_m + 0.592089376625*LC_nadph_m + 0.592089376625*LC_pyr_m = 0
FU_ME2m: 1.0*DG_ME2m + 1000.0*FU_ME2m < 999.999999
BU_ME2m: -1.0*DG_ME2m + 1000.0*BU_ME2m < 999.999999
SU_ME2m: 0 < 1.0*BU_ME2m + 1.0*FU_ME2m < 1
UF_ME2m: -1000.0*FU_ME2m + 1.0*ME2m < 0
UR_ME2m: -1000.0*BU_ME2m + 1.0*ME2m_reverse_f04a4 < 0
G_METt2r: -1.0*DG_METt2r + 0.592089376625*LC_met-L_c - 0.592089376625*LC_met-L_e + 1.0*DGo_METt2r = 0
FU_METt2r: 1.0*DG_METt2r + 1000.0*FU_METt2r < 999.999999
BU_METt2r: -1.0*DG_METt2r + 1000.0*BU_METt2r < 999.999999
SU_METt2r: 0 < 1.0*BU_METt2r + 1.0*FU_METt2r < 1
UF_METt2r: -1000.0*FU_METt2r + 1.0*METt2r < 0
UR_METt2r: -1000.0*BU_METt2r + 1.0*METt2r_reverse_22908 < 0
G_NADH2-u6cm: -1.0*DG_NADH2-u6cm + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadh_c - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m + 1.0*DGo_NADH2-u6cm = 0
FU_NADH2-u6cm: 1.0*DG_NADH2-u6cm + 1000.0*FU_NADH2-u6cm < 999.999999
BU_NADH2-u6cm: -1.0*DG_NADH2-u6cm + 1000.0*BU_NADH2-u6cm < 999.999999
SU_NADH2-u6cm: 0 < 1.0*BU_NADH2-u6cm + 1.0*FU_NADH2-u6cm < 1
UF_NADH2-u6cm: -1000.0*FU_NADH2-u6cm + 1.0*NADH2-u6cm < 0
UR_NADH2-u6cm: -1000.0*BU_NADH2-u6cm + 1.0*NADH2-u6cm_reverse_07733 < 0
G_NADH2-u6m: -1.0*DG_NADH2-u6m + 0.592089376625*LC_nad_m - 0.592089376625*LC_nadh_m - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m + 1.0*DGo_NADH2-u6m = 0
FU_NADH2-u6m: 1.0*DG_NADH2-u6m + 1000.0*FU_NADH2-u6m < 999.999999
BU_NADH2-u6m: -1.0*DG_NADH2-u6m + 1000.0*BU_NADH2-u6m < 999.999999
SU_NADH2-u6m: 0 < 1.0*BU_NADH2-u6m + 1.0*FU_NADH2-u6m < 1
UF_NADH2-u6m: -1000.0*FU_NADH2-u6m + 1.0*NADH2-u6m < 0
UR_NADH2-u6m: -1000.0*BU_NADH2-u6m + 1.0*NADH2-u6m_reverse_18109 < 0
G_NADK: -1.0*DG_NADK + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadp_c + 1.0*DGo_NADK = 0
FU_NADK: 1.0*DG_NADK + 1000.0*FU_NADK < 999.999999
BU_NADK: -1.0*DG_NADK + 1000.0*BU_NADK < 999.999999
SU_NADK: 0 < 1.0*BU_NADK + 1.0*FU_NADK < 1
UF_NADK: -1000.0*FU_NADK + 1.0*NADK < 0
UR_NADK: -1000.0*BU_NADK + 1.0*NADK_reverse_bba52 < 0
G_NADKm: -1.0*DG_NADKm + 0.592089376625*LC_adp_m - 0.592089376625*LC_atp_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadp_m + 1.0*DGo_NADKm = 0
FU_NADKm: 1.0*DG_NADKm + 1000.0*FU_NADKm < 999.999999
BU_NADKm: -1.0*DG_NADKm + 1000.0*BU_NADKm < 999.999999
SU_NADKm: 0 < 1.0*BU_NADKm + 1.0*FU_NADKm < 1
UF_NADKm: -1000.0*FU_NADKm + 1.0*NADKm < 0
UR_NADKm: -1000.0*BU_NADKm + 1.0*NADKm_reverse_f0454 < 0
G_NADPPPS: -1.0*DG_NADPPPS + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_pi_c + 1.0*DGo_NADPPPS = 0
FU_NADPPPS: 1.0*DG_NADPPPS + 1000.0*FU_NADPPPS < 999.999999
BU_NADPPPS: -1.0*DG_NADPPPS + 1000.0*BU_NADPPPS < 999.999999
SU_NADPPPS: 0 < 1.0*BU_NADPPPS + 1.0*FU_NADPPPS < 1
UF_NADPPPS: -1000.0*FU_NADPPPS + 1.0*NADPPPS < 0
UR_NADPPPS: -1000.0*BU_NADPPPS + 1.0*NADPPPS_reverse_c9c40 < 0
G_NADPPPSm: -1.0*DG_NADPPPSm + 0.592089376625*LC_nad_m - 0.592089376625*LC_nadp_m + 0.592089376625*LC_pi_m + 1.0*DGo_NADPPPSm = 0
FU_NADPPPSm: 1.0*DG_NADPPPSm + 1000.0*FU_NADPPPSm < 999.999999
BU_NADPPPSm: -1.0*DG_NADPPPSm + 1000.0*BU_NADPPPSm < 999.999999
SU_NADPPPSm: 0 < 1.0*BU_NADPPPSm + 1.0*FU_NADPPPSm < 1
UF_NADPPPSm: -1000.0*FU_NADPPPSm + 1.0*NADPPPSm < 0
UR_NADPPPSm: -1000.0*BU_NADPPPSm + 1.0*NADPPPSm_reverse_64d27 < 0
G_NDP1: -1.0*DG_NDP1 - 0.592089376625*LC_adp_c + 0.592089376625*LC_amp_c + 0.592089376625*LC_pi_c + 1.0*DGo_NDP1 = 0
FU_NDP1: 1.0*DG_NDP1 + 1000.0*FU_NDP1 < 999.999999
BU_NDP1: -1.0*DG_NDP1 + 1000.0*BU_NDP1 < 999.999999
SU_NDP1: 0 < 1.0*BU_NDP1 + 1.0*FU_NDP1 < 1
UF_NDP1: -1000.0*FU_NDP1 + 1.0*NDP1 < 0
UR_NDP1: -1000.0*BU_NDP1 + 1.0*NDP1_reverse_58712 < 0
G_NH4t: -1.0*DG_NH4t + 0.592089376625*LC_nh4_c - 0.592089376625*LC_nh4_e + 1.0*DGo_NH4t = 0
FU_NH4t: 1.0*DG_NH4t + 1000.0*FU_NH4t < 999.999999
BU_NH4t: -1.0*DG_NH4t + 1000.0*BU_NH4t < 999.999999
SU_NH4t: 0 < 1.0*BU_NH4t + 1.0*FU_NH4t < 1
UF_NH4t: -1000.0*FU_NH4t + 1.0*NH4t < 0
UR_NH4t: -1000.0*BU_NH4t + 1.0*NH4t_reverse_551ee < 0
G_NH4tm: -1.0*DG_NH4tm - 0.592089376625*LC_nh4_c + 0.592089376625*LC_nh4_m + 1.0*DGo_NH4tm = 0
FU_NH4tm: 1.0*DG_NH4tm + 1000.0*FU_NH4tm < 999.999999
BU_NH4tm: -1.0*DG_NH4tm + 1000.0*BU_NH4tm < 999.999999
SU_NH4tm: 0 < 1.0*BU_NH4tm + 1.0*FU_NH4tm < 1
UF_NH4tm: -1000.0*FU_NH4tm + 1.0*NH4tm < 0
UR_NH4tm: -1000.0*BU_NH4tm + 1.0*NH4tm_reverse_e28df < 0
G_O2t: -1.0*DG_O2t + 0.592089376625*LC_o2_c - 0.592089376625*LC_o2_e + 1.0*DGo_O2t = 0
FU_O2t: 1.0*DG_O2t + 1000.0*FU_O2t < 999.999999
BU_O2t: -1.0*DG_O2t + 1000.0*BU_O2t < 999.999999
SU_O2t: 0 < 1.0*BU_O2t + 1.0*FU_O2t < 1
UF_O2t: -1000.0*FU_O2t + 1.0*O2t < 0
UR_O2t: -1000.0*BU_O2t + 1.0*O2t_reverse_4d957 < 0
G_O2ter: -1.0*DG_O2ter - 0.592089376625*LC_o2_c + 0.592089376625*LC_o2_r + 1.0*DGo_O2ter = 0
FU_O2ter: 1.0*DG_O2ter + 1000.0*FU_O2ter < 999.999999
BU_O2ter: -1.0*DG_O2ter + 1000.0*BU_O2ter < 999.999999
SU_O2ter: 0 < 1.0*BU_O2ter + 1.0*FU_O2ter < 1
UF_O2ter: -1000.0*FU_O2ter + 1.0*O2ter < 0
UR_O2ter: -1000.0*BU_O2ter + 1.0*O2ter_reverse_ea1ad < 0
G_O2tm: -1.0*DG_O2tm - 0.592089376625*LC_o2_c + 0.592089376625*LC_o2_m + 1.0*DGo_O2tm = 0
FU_O2tm: 1.0*DG_O2tm + 1000.0*FU_O2tm < 999.999999
BU_O2tm: -1.0*DG_O2tm + 1000.0*BU_O2tm < 999.999999
SU_O2tm: 0 < 1.0*BU_O2tm + 1.0*FU_O2tm < 1
UF_O2tm: -1000.0*FU_O2tm + 1.0*O2tm < 0
UR_O2tm: -1000.0*BU_O2tm + 1.0*O2tm_reverse_fe4d6 < 0
G_OAAt: -1.0*DG_OAAt - 0.592089376625*LC_oaa_c + 0.592089376625*LC_oaa_e + 1.0*DGo_OAAt = 0
FU_OAAt: 1.0*DG_OAAt + 1000.0*FU_OAAt < 999.999999
BU_OAAt: -1.0*DG_OAAt + 1000.0*BU_OAAt < 999.999999
SU_OAAt: 0 < 1.0*BU_OAAt + 1.0*FU_OAAt < 1
UF_OAAt: -1000.0*FU_OAAt + 1.0*OAAt < 0
UR_OAAt: -1000.0*BU_OAAt + 1.0*OAAt_reverse_d0cc5 < 0
G_OAAt2m: -1.0*DG_OAAt2m - 0.592089376625*LC_oaa_c + 0.592089376625*LC_oaa_m + 1.0*DGo_OAAt2m = 0
FU_OAAt2m: 1.0*DG_OAAt2m + 1000.0*FU_OAAt2m < 999.999999
BU_OAAt2m: -1.0*DG_OAAt2m + 1000.0*BU_OAAt2m < 999.999999
SU_OAAt2m: 0 < 1.0*BU_OAAt2m + 1.0*FU_OAAt2m < 1
UF_OAAt2m: -1000.0*FU_OAAt2m + 1.0*OAAt2m < 0
UR_OAAt2m: -1000.0*BU_OAAt2m + 1.0*OAAt2m_reverse_42ec8 < 0
G_PC: -1.0*DG_PC + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c - 0.592089376625*LC_hco3_c + 0.592089376625*LC_oaa_c + 0.592089376625*LC_pi_c - 0.592089376625*LC_pyr_c + 1.0*DGo_PC = 0
FU_PC: 1.0*DG_PC + 1000.0*FU_PC < 999.999999
BU_PC: -1.0*DG_PC + 1000.0*BU_PC < 999.999999
SU_PC: 0 < 1.0*BU_PC + 1.0*FU_PC < 1
UF_PC: -1000.0*FU_PC + 1.0*PC < 0
UR_PC: -1000.0*BU_PC + 1.0*PC_reverse_88dba < 0
G_PDHcm: -1.0*DG_PDHcm - 0.592089376625*LC_dhlam_m + 0.592089376625*LC_lpam_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m + 1.0*DGo_PDHcm = 0
FU_PDHcm: 1.0*DG_PDHcm + 1000.0*FU_PDHcm < 999.999999
BU_PDHcm: -1.0*DG_PDHcm + 1000.0*BU_PDHcm < 999.999999
SU_PDHcm: 0 < 1.0*BU_PDHcm + 1.0*FU_PDHcm < 1
UF_PDHcm: -1000.0*FU_PDHcm + 1.0*PDHcm < 0
UR_PDHcm: -1000.0*BU_PDHcm + 1.0*PDHcm_reverse_40ed3 < 0
G_PDHm: -1.0*DG_PDHm + 1.0*DGo_PDHm + 0.592089376625*LC_accoa_m + 0.592089376625*LC_co2_m - 0.592089376625*LC_coa_m - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m - 0.592089376625*LC_pyr_m = 0
FU_PDHm: 1.0*DG_PDHm + 1000.0*FU_PDHm < 999.999999
BU_PDHm: -1.0*DG_PDHm + 1000.0*BU_PDHm < 999.999999
SU_PDHm: 0 < 1.0*BU_PDHm + 1.0*FU_PDHm < 1
UF_PDHm: -1000.0*FU_PDHm + 1.0*PDHm < 0
UR_PDHm: -1000.0*BU_PDHm + 1.0*PDHm_reverse_31d17 < 0
SU_PEtm_SC: 0 < 1.0*BU_PEtm_SC + 1.0*FU_PEtm_SC < 1
UF_PEtm_SC: -1000.0*FU_PEtm_SC + 1.0*PEtm_SC < 0
UR_PEtm_SC: -1000.0*BU_PEtm_SC + 1.0*PEtm_SC_reverse_a93d3 < 0
G_PFK: -1.0*DG_PFK + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c - 0.592089376625*LC_f6p_c + 0.592089376625*LC_fdp_c + 1.0*DGo_PFK = 0
FU_PFK: 1.0*DG_PFK + 1000.0*FU_PFK < 999.999999
BU_PFK: -1.0*DG_PFK + 1000.0*BU_PFK < 999.999999
SU_PFK: 0 < 1.0*BU_PFK + 1.0*FU_PFK < 1
UF_PFK: -1000.0*FU_PFK + 1.0*PFK < 0
UR_PFK: -1000.0*BU_PFK + 1.0*PFK_reverse_d24a6 < 0
G_PFK_3: -1.0*DG_PFK_3 + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_s17bp_c - 0.592089376625*LC_s7p_c + 1.0*DGo_PFK_3 = 0
FU_PFK_3: 1.0*DG_PFK_3 + 1000.0*FU_PFK_3 < 999.999999
BU_PFK_3: -1.0*DG_PFK_3 + 1000.0*BU_PFK_3 < 999.999999
SU_PFK_3: 0 < 1.0*BU_PFK_3 + 1.0*FU_PFK_3 < 1
UF_PFK_3: -1000.0*FU_PFK_3 + 1.0*PFK_3 < 0
UR_PFK_3: -1000.0*BU_PFK_3 + 1.0*PFK_3_reverse_037e0 < 0
G_PGI: 1.0*DGo_PGI - 1.0*DG_PGI + 0.592089376625*LC_f6p_c - 0.592089376625*LC_g6p_c = 0
FU_PGI: 1.0*DG_PGI + 1000.0*FU_PGI < 999.999999
BU_PGI: -1.0*DG_PGI + 1000.0*BU_PGI < 999.999999
SU_PGI: 0 < 1.0*BU_PGI + 1.0*FU_PGI < 1
UF_PGI: -1000.0*FU_PGI + 1.0*PGI < 0
UR_PGI: -1000.0*BU_PGI + 1.0*PGI_reverse_27efc < 0
G_PGK: -1.0*DG_PGK + 0.592089376625*LC_13dpg_c - 0.592089376625*LC_3pg_c + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 1.0*DGo_PGK = 0
FU_PGK: 1.0*DG_PGK + 1000.0*FU_PGK < 999.999999
BU_PGK: -1.0*DG_PGK + 1000.0*BU_PGK < 999.999999
SU_PGK: 0 < 1.0*BU_PGK + 1.0*FU_PGK < 1
UF_PGK: -1000.0*FU_PGK + 1.0*PGK < 0
UR_PGK: -1000.0*BU_PGK + 1.0*PGK_reverse_02696 < 0
G_PGL: -1.0*DG_PGL + 0.592089376625*LC_6pgc_c - 0.592089376625*LC_6pgl_c + 1.0*DGo_PGL = 0
FU_PGL: 1.0*DG_PGL + 1000.0*FU_PGL < 999.999999
BU_PGL: -1.0*DG_PGL + 1000.0*BU_PGL < 999.999999
SU_PGL: 0 < 1.0*BU_PGL + 1.0*FU_PGL < 1
UF_PGL: -1000.0*FU_PGL + 1.0*PGL < 0
UR_PGL: -1000.0*BU_PGL + 1.0*PGL_reverse_2bb6b < 0
G_PGM: -1.0*DG_PGM - 0.592089376625*LC_2pg_c + 0.592089376625*LC_3pg_c + 1.0*DGo_PGM = 0
FU_PGM: 1.0*DG_PGM + 1000.0*FU_PGM < 999.999999
BU_PGM: -1.0*DG_PGM + 1000.0*BU_PGM < 999.999999
SU_PGM: 0 < 1.0*BU_PGM + 1.0*FU_PGM < 1
UF_PGM: -1000.0*FU_PGM + 1.0*PGM < 0
UR_PGM: -1000.0*BU_PGM + 1.0*PGM_reverse_fc9af < 0
G_PGMT: -1.0*DG_PGMT - 0.592089376625*LC_g1p_c + 0.592089376625*LC_g6p_c + 1.0*DGo_PGMT = 0
FU_PGMT: 1.0*DG_PGMT + 1000.0*FU_PGMT < 999.999999
BU_PGMT: -1.0*DG_PGMT + 1000.0*BU_PGMT < 999.999999
SU_PGMT: 0 < 1.0*BU_PGMT + 1.0*FU_PGMT < 1
UF_PGMT: -1000.0*FU_PGMT + 1.0*PGMT < 0
UR_PGMT: -1000.0*BU_PGMT + 1.0*PGMT_reverse_5bcdd < 0
G_PHEt2r: -1.0*DG_PHEt2r + 0.592089376625*LC_phe-L_c - 0.592089376625*LC_phe-L_e + 1.0*DGo_PHEt2r = 0
FU_PHEt2r: 1.0*DG_PHEt2r + 1000.0*FU_PHEt2r < 999.999999
BU_PHEt2r: -1.0*DG_PHEt2r + 1000.0*BU_PHEt2r < 999.999999
SU_PHEt2r: 0 < 1.0*BU_PHEt2r + 1.0*FU_PHEt2r < 1
UF_PHEt2r: -1000.0*FU_PHEt2r + 1.0*PHEt2r < 0
UR_PHEt2r: -1000.0*BU_PHEt2r + 1.0*PHEt2r_reverse_19e89 < 0
G_PIt2m: -1.0*DG_PIt2m - 0.592089376625*LC_pi_c + 0.592089376625*LC_pi_m + 1.0*DGo_PIt2m = 0
FU_PIt2m: 1.0*DG_PIt2m + 1000.0*FU_PIt2m < 999.999999
BU_PIt2m: -1.0*DG_PIt2m + 1000.0*BU_PIt2m < 999.999999
SU_PIt2m: 0 < 1.0*BU_PIt2m + 1.0*FU_PIt2m < 1
UF_PIt2m: -1000.0*FU_PIt2m + 1.0*PIt2m < 0
UR_PIt2m: -1000.0*BU_PIt2m + 1.0*PIt2m_reverse_cc401 < 0
G_PIt2r: 1.0*DGo_PIt2r - 1.0*DG_PIt2r + 0.592089376625*LC_pi_c - 0.592089376625*LC_pi_e = 0
FU_PIt2r: 1.0*DG_PIt2r + 1000.0*FU_PIt2r < 999.999999
BU_PIt2r: -1.0*DG_PIt2r + 1000.0*BU_PIt2r < 999.999999
SU_PIt2r: 0 < 1.0*BU_PIt2r + 1.0*FU_PIt2r < 1
UF_PIt2r: -1000.0*FU_PIt2r + 1.0*PIt2r < 0
UR_PIt2r: -1000.0*BU_PIt2r + 1.0*PIt2r_reverse_1cd61 < 0
G_PPA: -1.0*DG_PPA + 1.18417875325*LC_pi_c - 0.592089376625*LC_ppi_c + 1.0*DGo_PPA = 0
FU_PPA: 1.0*DG_PPA + 1000.0*FU_PPA < 999.999999
BU_PPA: -1.0*DG_PPA + 1000.0*BU_PPA < 999.999999
SU_PPA: 0 < 1.0*BU_PPA + 1.0*FU_PPA < 1
UF_PPA: -1000.0*FU_PPA + 1.0*PPA < 0
UR_PPA: -1000.0*BU_PPA + 1.0*PPA_reverse_c5293 < 0
G_PPAm: -1.0*DG_PPAm + 1.18417875325*LC_pi_m - 0.592089376625*LC_ppi_m + 1.0*DGo_PPAm = 0
FU_PPAm: 1.0*DG_PPAm + 1000.0*FU_PPAm < 999.999999
BU_PPAm: -1.0*DG_PPAm + 1000.0*BU_PPAm < 999.999999
SU_PPAm: 0 < 1.0*BU_PPAm + 1.0*FU_PPAm < 1
UF_PPAm: -1000.0*FU_PPAm + 1.0*PPAm < 0
UR_PPAm: -1000.0*BU_PPAm + 1.0*PPAm_reverse_6c806 < 0
G_PPCK: 1.0*DGo_PPCK - 1.0*DG_PPCK + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_oaa_c + 0.592089376625*LC_pep_c = 0
FU_PPCK: 1.0*DG_PPCK + 1000.0*FU_PPCK < 999.999999
BU_PPCK: -1.0*DG_PPCK + 1000.0*BU_PPCK < 999.999999
SU_PPCK: 0 < 1.0*BU_PPCK + 1.0*FU_PPCK < 1
UF_PPCK: -1000.0*FU_PPCK + 1.0*PPCK < 0
UR_PPCK: -1000.0*BU_PPCK + 1.0*PPCK_reverse_2557d < 0
G_PROt2r: -1.0*DG_PROt2r + 0.592089376625*LC_pro-L_c - 0.592089376625*LC_pro-L_e + 1.0*DGo_PROt2r = 0
FU_PROt2r: 1.0*DG_PROt2r + 1000.0*FU_PROt2r < 999.999999
BU_PROt2r: -1.0*DG_PROt2r + 1000.0*BU_PROt2r < 999.999999
SU_PROt2r: 0 < 1.0*BU_PROt2r + 1.0*FU_PROt2r < 1
UF_PROt2r: -1000.0*FU_PROt2r + 1.0*PROt2r < 0
UR_PROt2r: -1000.0*BU_PROt2r + 1.0*PROt2r_reverse_f59df < 0
SU_PStm_SC: 0 < 1.0*BU_PStm_SC + 1.0*FU_PStm_SC < 1
UF_PStm_SC: -1000.0*FU_PStm_SC + 1.0*PStm_SC < 0
UR_PStm_SC: -1000.0*BU_PStm_SC + 1.0*PStm_SC_reverse_a026e < 0
G_PYK: -1.0*DG_PYK - 0.592089376625*LC_adp_c + 0.592089376625*LC_atp_c - 0.592089376625*LC_pep_c + 0.592089376625*LC_pyr_c + 1.0*DGo_PYK = 0
FU_PYK: 1.0*DG_PYK + 1000.0*FU_PYK < 999.999999
BU_PYK: -1.0*DG_PYK + 1000.0*BU_PYK < 999.999999
SU_PYK: 0 < 1.0*BU_PYK + 1.0*FU_PYK < 1
UF_PYK: -1000.0*FU_PYK + 1.0*PYK < 0
UR_PYK: -1000.0*BU_PYK + 1.0*PYK_reverse_bc8ff < 0
G_PYRDC: -1.0*DG_PYRDC + 0.592089376625*LC_acald_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_pyr_c + 1.0*DGo_PYRDC = 0
FU_PYRDC: 1.0*DG_PYRDC + 1000.0*FU_PYRDC < 999.999999
BU_PYRDC: -1.0*DG_PYRDC + 1000.0*BU_PYRDC < 999.999999
SU_PYRDC: 0 < 1.0*BU_PYRDC + 1.0*FU_PYRDC < 1
UF_PYRDC: -1000.0*FU_PYRDC + 1.0*PYRDC < 0
UR_PYRDC: -1000.0*BU_PYRDC + 1.0*PYRDC_reverse_51259 < 0
G_PYRt: -1.0*DG_PYRt - 0.592089376625*LC_pyr_c + 0.592089376625*LC_pyr_e + 1.0*DGo_PYRt = 0
FU_PYRt: 1.0*DG_PYRt + 1000.0*FU_PYRt < 999.999999
BU_PYRt: -1.0*DG_PYRt + 1000.0*BU_PYRt < 999.999999
SU_PYRt: 0 < 1.0*BU_PYRt + 1.0*FU_PYRt < 1
UF_PYRt: -1000.0*FU_PYRt + 1.0*PYRt < 0
UR_PYRt: -1000.0*BU_PYRt + 1.0*PYRt_reverse_6fee7 < 0
G_PYRt2m: -1.0*DG_PYRt2m - 0.592089376625*LC_pyr_c + 0.592089376625*LC_pyr_m + 1.0*DGo_PYRt2m = 0
FU_PYRt2m: 1.0*DG_PYRt2m + 1000.0*FU_PYRt2m < 999.999999
BU_PYRt2m: -1.0*DG_PYRt2m + 1000.0*BU_PYRt2m < 999.999999
SU_PYRt2m: 0 < 1.0*BU_PYRt2m + 1.0*FU_PYRt2m < 1
UF_PYRt2m: -1000.0*FU_PYRt2m + 1.0*PYRt2m < 0
UR_PYRt2m: -1000.0*BU_PYRt2m + 1.0*PYRt2m_reverse_fb762 < 0
G_RPE: -1.0*DG_RPE - 0.592089376625*LC_ru5p-D_c + 0.592089376625*LC_xu5p-D_c + 1.0*DGo_RPE = 0
FU_RPE: 1.0*DG_RPE + 1000.0*FU_RPE < 999.999999
BU_RPE: -1.0*DG_RPE + 1000.0*BU_RPE < 999.999999
SU_RPE: 0 < 1.0*BU_RPE + 1.0*FU_RPE < 1
UF_RPE: -1000.0*FU_RPE + 1.0*RPE < 0
UR_RPE: -1000.0*BU_RPE + 1.0*RPE_reverse_a1b04 < 0
G_RPI: 1.0*DGo_RPI - 1.0*DG_RPI - 0.592089376625*LC_r5p_c + 0.592089376625*LC_ru5p-D_c = 0
FU_RPI: 1.0*DG_RPI + 1000.0*FU_RPI < 999.999999
BU_RPI: -1.0*DG_RPI + 1000.0*BU_RPI < 999.999999
SU_RPI: 0 < 1.0*BU_RPI + 1.0*FU_RPI < 1
UF_RPI: -1000.0*FU_RPI + 1.0*RPI < 0
UR_RPI: -1000.0*BU_RPI + 1.0*RPI_reverse_853a1 < 0
G_SERt2r: -1.0*DG_SERt2r + 0.592089376625*LC_ser-L_c - 0.592089376625*LC_ser-L_e + 1.0*DGo_SERt2r = 0
FU_SERt2r: 1.0*DG_SERt2r + 1000.0*FU_SERt2r < 999.999999
BU_SERt2r: -1.0*DG_SERt2r + 1000.0*BU_SERt2r < 999.999999
SU_SERt2r: 0 < 1.0*BU_SERt2r + 1.0*FU_SERt2r < 1
UF_SERt2r: -1000.0*FU_SERt2r + 1.0*SERt2r < 0
UR_SERt2r: -1000.0*BU_SERt2r + 1.0*SERt2r_reverse_1b436 < 0
G_SO4ti: 1.0*DGo_SO4ti - 1.0*DG_SO4ti + 0.592089376625*LC_so4_c - 0.592089376625*LC_so4_e = 0
FU_SO4ti: 1.0*DG_SO4ti + 1000.0*FU_SO4ti < 999.999999
BU_SO4ti: -1.0*DG_SO4ti + 1000.0*BU_SO4ti < 999.999999
SU_SO4ti: 0 < 1.0*BU_SO4ti + 1.0*FU_SO4ti < 1
UF_SO4ti: -1000.0*FU_SO4ti + 1.0*SO4ti < 0
UR_SO4ti: -1000.0*BU_SO4ti + 1.0*SO4ti_reverse_4a93d < 0
G_SUCCt2r: -1.0*DG_SUCCt2r + 0.592089376625*LC_succ_c - 0.592089376625*LC_succ_e + 1.0*DGo_SUCCt2r = 0
FU_SUCCt2r: 1.0*DG_SUCCt2r + 1000.0*FU_SUCCt2r < 999.999999
BU_SUCCt2r: -1.0*DG_SUCCt2r + 1000.0*BU_SUCCt2r < 999.999999
SU_SUCCt2r: 0 < 1.0*BU_SUCCt2r + 1.0*FU_SUCCt2r < 1
UF_SUCCt2r: -1000.0*FU_SUCCt2r + 1.0*SUCCt2r < 0
UR_SUCCt2r: -1000.0*BU_SUCCt2r + 1.0*SUCCt2r_reverse_e3abe < 0
G_SUCCtm: -1.0*DG_SUCCtm + 0.592089376625*LC_pi_c - 0.592089376625*LC_pi_m - 0.592089376625*LC_succ_c + 0.592089376625*LC_succ_m + 1.0*DGo_SUCCtm = 0
FU_SUCCtm: 1.0*DG_SUCCtm + 1000.0*FU_SUCCtm < 999.999999
BU_SUCCtm: -1.0*DG_SUCCtm + 1000.0*BU_SUCCtm < 999.999999
SU_SUCCtm: 0 < 1.0*BU_SUCCtm + 1.0*FU_SUCCtm < 1
UF_SUCCtm: -1000.0*FU_SUCCtm + 1.0*SUCCtm < 0
UR_SUCCtm: -1000.0*BU_SUCCtm + 1.0*SUCCtm_reverse_5f36b < 0
G_SUCD1m: -1.0*DG_SUCD1m - 0.592089376625*LC_fad_m + 0.592089376625*LC_fadh2_m + 0.592089376625*LC_fum_m - 0.592089376625*LC_succ_m + 1.0*DGo_SUCD1m = 0
FU_SUCD1m: 1.0*DG_SUCD1m + 1000.0*FU_SUCD1m < 999.999999
BU_SUCD1m: -1.0*DG_SUCD1m + 1000.0*BU_SUCD1m < 999.999999
SU_SUCD1m: 0 < 1.0*BU_SUCD1m + 1.0*FU_SUCD1m < 1
UF_SUCD1m: -1000.0*FU_SUCD1m + 1.0*SUCD1m < 0
UR_SUCD1m: -1000.0*BU_SUCD1m + 1.0*SUCD1m_reverse_60ba1 < 0
G_SUCD2_u6m: -1.0*DG_SUCD2_u6m + 0.592089376625*LC_fum_m - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m - 0.592089376625*LC_succ_m + 1.0*DGo_SUCD2_u6m = 0
FU_SUCD2_u6m: 1.0*DG_SUCD2_u6m + 1000.0*FU_SUCD2_u6m < 999.999999
BU_SUCD2_u6m: -1.0*DG_SUCD2_u6m + 1000.0*BU_SUCD2_u6m < 999.999999
SU_SUCD2_u6m: 0 < 1.0*BU_SUCD2_u6m + 1.0*FU_SUCD2_u6m < 1
UF_SUCD2_u6m: -1000.0*FU_SUCD2_u6m + 1.0*SUCD2_u6m < 0
UR_SUCD2_u6m: -1000.0*BU_SUCD2_u6m + 1.0*SUCD2_u6m_reverse_b6f79 < 0
G_SUCD3-u6m: -1.0*DG_SUCD3-u6m + 0.592089376625*LC_fad_m - 0.592089376625*LC_fadh2_m - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m + 1.0*DGo_SUCD3-u6m = 0
FU_SUCD3-u6m: 1.0*DG_SUCD3-u6m + 1000.0*FU_SUCD3-u6m < 999.999999
BU_SUCD3-u6m: -1.0*DG_SUCD3-u6m + 1000.0*BU_SUCD3-u6m < 999.999999
SU_SUCD3-u6m: 0 < 1.0*BU_SUCD3-u6m + 1.0*FU_SUCD3-u6m < 1
UF_SUCD3-u6m: -1000.0*FU_SUCD3-u6m + 1.0*SUCD3-u6m < 0
UR_SUCD3-u6m: -1000.0*BU_SUCD3-u6m + 1.0*SUCD3-u6m_reverse_0e152 < 0
G_SUCFUMtm: 1.0*DGo_SUCFUMtm - 1.0*DG_SUCFUMtm + 0.592089376625*LC_fum_c - 0.592089376625*LC_fum_m - 0.592089376625*LC_succ_c + 0.592089376625*LC_succ_m = 0
FU_SUCFUMtm: 1.0*DG_SUCFUMtm + 1000.0*FU_SUCFUMtm < 999.999999
BU_SUCFUMtm: -1.0*DG_SUCFUMtm + 1000.0*BU_SUCFUMtm < 999.999999
SU_SUCFUMtm: 0 < 1.0*BU_SUCFUMtm + 1.0*FU_SUCFUMtm < 1
UF_SUCFUMtm: -1000.0*FU_SUCFUMtm + 1.0*SUCFUMtm < 0
UR_SUCFUMtm: -1000.0*BU_SUCFUMtm + 1.0*SUCFUMtm_reverse_99c3a < 0
G_SUCOASm: -1.0*DG_SUCOASm + 0.592089376625*LC_adp_m - 0.592089376625*LC_atp_m - 0.592089376625*LC_coa_m + 0.592089376625*LC_pi_m - 0.592089376625*LC_succ_m + 0.592089376625*LC_succoa_m + 1.0*DGo_SUCOASm = 0
FU_SUCOASm: 1.0*DG_SUCOASm + 1000.0*FU_SUCOASm < 999.999999
BU_SUCOASm: -1.0*DG_SUCOASm + 1000.0*BU_SUCOASm < 999.999999
SU_SUCOASm: 0 < 1.0*BU_SUCOASm + 1.0*FU_SUCOASm < 1
UF_SUCOASm: -1000.0*FU_SUCOASm + 1.0*SUCOASm < 0
UR_SUCOASm: -1000.0*BU_SUCOASm + 1.0*SUCOASm_reverse_87069 < 0
G_TALA: -1.0*DG_TALA + 0.592089376625*LC_e4p_c + 0.592089376625*LC_f6p_c - 0.592089376625*LC_g3p_c - 0.592089376625*LC_s7p_c + 1.0*DGo_TALA = 0
FU_TALA: 1.0*DG_TALA + 1000.0*FU_TALA < 999.999999
BU_TALA: -1.0*DG_TALA + 1000.0*BU_TALA < 999.999999
SU_TALA: 0 < 1.0*BU_TALA + 1.0*FU_TALA < 1
UF_TALA: -1000.0*FU_TALA + 1.0*TALA < 0
UR_TALA: -1000.0*BU_TALA + 1.0*TALA_reverse_adfda < 0
G_THRt2m: -1.0*DG_THRt2m - 0.592089376625*LC_thr-L_c + 0.592089376625*LC_thr-L_m + 1.0*DGo_THRt2m = 0
FU_THRt2m: 1.0*DG_THRt2m + 1000.0*FU_THRt2m < 999.999999
BU_THRt2m: -1.0*DG_THRt2m + 1000.0*BU_THRt2m < 999.999999
SU_THRt2m: 0 < 1.0*BU_THRt2m + 1.0*FU_THRt2m < 1
UF_THRt2m: -1000.0*FU_THRt2m + 1.0*THRt2m < 0
UR_THRt2m: -1000.0*BU_THRt2m + 1.0*THRt2m_reverse_fd5a8 < 0
G_THRt2r: -1.0*DG_THRt2r + 0.592089376625*LC_thr-L_c - 0.592089376625*LC_thr-L_e + 1.0*DGo_THRt2r = 0
FU_THRt2r: 1.0*DG_THRt2r + 1000.0*FU_THRt2r < 999.999999
BU_THRt2r: -1.0*DG_THRt2r + 1000.0*BU_THRt2r < 999.999999
SU_THRt2r: 0 < 1.0*BU_THRt2r + 1.0*FU_THRt2r < 1
UF_THRt2r: -1000.0*FU_THRt2r + 1.0*THRt2r < 0
UR_THRt2r: -1000.0*BU_THRt2r + 1.0*THRt2r_reverse_721ff < 0
G_TKT1: 1.0*DGo_TKT1 - 1.0*DG_TKT1 + 0.592089376625*LC_g3p_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_s7p_c - 0.592089376625*LC_xu5p-D_c = 0
FU_TKT1: 1.0*DG_TKT1 + 1000.0*FU_TKT1 < 999.999999
BU_TKT1: -1.0*DG_TKT1 + 1000.0*BU_TKT1 < 999.999999
SU_TKT1: 0 < 1.0*BU_TKT1 + 1.0*FU_TKT1 < 1
UF_TKT1: -1000.0*FU_TKT1 + 1.0*TKT1 < 0
UR_TKT1: -1000.0*BU_TKT1 + 1.0*TKT1_reverse_a1021 < 0
G_TKT2: -1.0*DG_TKT2 - 0.592089376625*LC_e4p_c + 0.592089376625*LC_f6p_c + 0.592089376625*LC_g3p_c - 0.592089376625*LC_xu5p-D_c + 1.0*DGo_TKT2 = 0
FU_TKT2: 1.0*DG_TKT2 + 1000.0*FU_TKT2 < 999.999999
BU_TKT2: -1.0*DG_TKT2 + 1000.0*BU_TKT2 < 999.999999
SU_TKT2: 0 < 1.0*BU_TKT2 + 1.0*FU_TKT2 < 1
UF_TKT2: -1000.0*FU_TKT2 + 1.0*TKT2 < 0
UR_TKT2: -1000.0*BU_TKT2 + 1.0*TKT2_reverse_7ebc7 < 0
G_TPI: -1.0*DG_TPI - 0.592089376625*LC_dhap_c + 0.592089376625*LC_g3p_c + 1.0*DGo_TPI = 0
FU_TPI: 1.0*DG_TPI + 1000.0*FU_TPI < 999.999999
BU_TPI: -1.0*DG_TPI + 1000.0*BU_TPI < 999.999999
SU_TPI: 0 < 1.0*BU_TPI + 1.0*FU_TPI < 1
UF_TPI: -1000.0*FU_TPI + 1.0*TPI < 0
UR_TPI: -1000.0*BU_TPI + 1.0*TPI_reverse_c2c3b < 0
G_TRPt2r: -1.0*DG_TRPt2r + 0.592089376625*LC_trp-L_c - 0.592089376625*LC_trp-L_e + 1.0*DGo_TRPt2r = 0
FU_TRPt2r: 1.0*DG_TRPt2r + 1000.0*FU_TRPt2r < 999.999999
BU_TRPt2r: -1.0*DG_TRPt2r + 1000.0*BU_TRPt2r < 999.999999
SU_TRPt2r: 0 < 1.0*BU_TRPt2r + 1.0*FU_TRPt2r < 1
UF_TRPt2r: -1000.0*FU_TRPt2r + 1.0*TRPt2r < 0
UR_TRPt2r: -1000.0*BU_TRPt2r + 1.0*TRPt2r_reverse_159c4 < 0
G_TYRt2r: -1.0*DG_TYRt2r + 0.592089376625*LC_tyr-L_c - 0.592089376625*LC_tyr-L_e + 1.0*DGo_TYRt2r = 0
FU_TYRt2r: 1.0*DG_TYRt2r + 1000.0*FU_TYRt2r < 999.999999
BU_TYRt2r: -1.0*DG_TYRt2r + 1000.0*BU_TYRt2r < 999.999999
SU_TYRt2r: 0 < 1.0*BU_TYRt2r + 1.0*FU_TYRt2r < 1
UF_TYRt2r: -1000.0*FU_TYRt2r + 1.0*TYRt2r < 0
UR_TYRt2r: -1000.0*BU_TYRt2r + 1.0*TYRt2r_reverse_6f73b < 0
G_VALt2r: -1.0*DG_VALt2r + 0.592089376625*LC_val-L_c - 0.592089376625*LC_val-L_e + 1.0*DGo_VALt2r = 0
FU_VALt2r: 1.0*DG_VALt2r + 1000.0*FU_VALt2r < 999.999999
BU_VALt2r: -1.0*DG_VALt2r + 1000.0*BU_VALt2r < 999.999999
SU_VALt2r: 0 < 1.0*BU_VALt2r + 1.0*FU_VALt2r < 1
UF_VALt2r: -1000.0*FU_VALt2r + 1.0*VALt2r < 0
UR_VALt2r: -1000.0*BU_VALt2r + 1.0*VALt2r_reverse_e3d24 < 0
G_ZYMSTt: -1.0*DG_ZYMSTt + 0.592089376625*LC_zymst_c - 0.592089376625*LC_zymst_e + 1.0*DGo_ZYMSTt = 0
FU_ZYMSTt: 1.0*DG_ZYMSTt + 1000.0*FU_ZYMSTt < 999.999999
BU_ZYMSTt: -1.0*DG_ZYMSTt + 1000.0*BU_ZYMSTt < 999.999999
SU_ZYMSTt: 0 < 1.0*BU_ZYMSTt + 1.0*FU_ZYMSTt < 1
UF_ZYMSTt: -1000.0*FU_ZYMSTt + 1.0*ZYMSTt < 0
UR_ZYMSTt: -1000.0*BU_ZYMSTt + 1.0*ZYMSTt_reverse_ab42f < 0
SU_LMPD_1_13BDglcn_c: 0 < 1.0*BU_LMPD_1_13BDglcn_c + 1.0*FU_LMPD_1_13BDglcn_c < 1
UF_LMPD_1_13BDglcn_c: -1000.0*FU_LMPD_1_13BDglcn_c + 1.0*LMPD_1_13BDglcn_c < 0
UR_LMPD_1_13BDglcn_c: -1000.0*BU_LMPD_1_13BDglcn_c + 1.0*LMPD_1_13BDglcn_c_reverse_e2bdf < 0
G_LMPD_2_amp_c: -1.0*DG_LMPD_2_amp_c - 1.18417875325*LC_3pg_c + 4.144625636375*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 1.18417875325*LC_asp-L_c - 4.736715013*LC_atp_c - 0.592089376625*LC_co2_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 1.18417875325*LC_nad_c + 1.18417875325*LC_nadh_c - 1.18417875325*LC_nadp_c + 1.18417875325*LC_nadph_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_2_amp_c = 0
FU_LMPD_2_amp_c: 1.0*DG_LMPD_2_amp_c + 1000.0*FU_LMPD_2_amp_c < 999.999999
BU_LMPD_2_amp_c: -1.0*DG_LMPD_2_amp_c + 1000.0*BU_LMPD_2_amp_c < 999.999999
SU_LMPD_2_amp_c: 0 < 1.0*BU_LMPD_2_amp_c + 1.0*FU_LMPD_2_amp_c < 1
UF_LMPD_2_amp_c: -1000.0*FU_LMPD_2_amp_c + 1.0*LMPD_2_amp_c < 0
UR_LMPD_2_amp_c: -1000.0*BU_LMPD_2_amp_c + 1.0*LMPD_2_amp_c_reverse_bde95 < 0
G_LMPD_3_amp_c: -1.0*DG_LMPD_3_amp_c - 1.18417875325*LC_3pg_c + 4.144625636375*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 1.18417875325*LC_asp-L_c - 4.736715013*LC_atp_c - 0.592089376625*LC_co2_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 2.3683575065*LC_nad_c + 2.3683575065*LC_nadh_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_3_amp_c = 0
FU_LMPD_3_amp_c: 1.0*DG_LMPD_3_amp_c + 1000.0*FU_LMPD_3_amp_c < 999.999999
BU_LMPD_3_amp_c: -1.0*DG_LMPD_3_amp_c + 1000.0*BU_LMPD_3_amp_c < 999.999999
SU_LMPD_3_amp_c: 0 < 1.0*BU_LMPD_3_amp_c + 1.0*FU_LMPD_3_amp_c < 1
UF_LMPD_3_amp_c: -1000.0*FU_LMPD_3_amp_c + 1.0*LMPD_3_amp_c < 0
UR_LMPD_3_amp_c: -1000.0*BU_LMPD_3_amp_c + 1.0*LMPD_3_amp_c_reverse_96e30 < 0
G_LMPD_4_amp_c: -1.0*DG_LMPD_4_amp_c - 1.18417875325*LC_3pg_c + 2.960446883125*LC_adp_c + 1.18417875325*LC_akg_c + 1.776268129875*LC_amp_c - 1.18417875325*LC_asp-L_c - 4.144625636375*LC_atp_c - 0.592089376625*LC_co2_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 1.18417875325*LC_nad_c + 1.18417875325*LC_nadh_c - 1.18417875325*LC_nadp_c + 1.18417875325*LC_nadph_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_4_amp_c = 0
FU_LMPD_4_amp_c: 1.0*DG_LMPD_4_amp_c + 1000.0*FU_LMPD_4_amp_c < 999.999999
BU_LMPD_4_amp_c: -1.0*DG_LMPD_4_amp_c + 1000.0*BU_LMPD_4_amp_c < 999.999999
SU_LMPD_4_amp_c: 0 < 1.0*BU_LMPD_4_amp_c + 1.0*FU_LMPD_4_amp_c < 1
UF_LMPD_4_amp_c: -1000.0*FU_LMPD_4_amp_c + 1.0*LMPD_4_amp_c < 0
UR_LMPD_4_amp_c: -1000.0*BU_LMPD_4_amp_c + 1.0*LMPD_4_amp_c_reverse_7925f < 0
G_LMPD_5_amp_c: -1.0*DG_LMPD_5_amp_c - 1.18417875325*LC_3pg_c + 2.960446883125*LC_adp_c + 1.18417875325*LC_akg_c + 1.776268129875*LC_amp_c - 1.18417875325*LC_asp-L_c - 4.144625636375*LC_atp_c - 0.592089376625*LC_co2_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 2.3683575065*LC_nad_c + 2.3683575065*LC_nadh_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_5_amp_c = 0
FU_LMPD_5_amp_c: 1.0*DG_LMPD_5_amp_c + 1000.0*FU_LMPD_5_amp_c < 999.999999
BU_LMPD_5_amp_c: -1.0*DG_LMPD_5_amp_c + 1000.0*BU_LMPD_5_amp_c < 999.999999
SU_LMPD_5_amp_c: 0 < 1.0*BU_LMPD_5_amp_c + 1.0*FU_LMPD_5_amp_c < 1
UF_LMPD_5_amp_c: -1000.0*FU_LMPD_5_amp_c + 1.0*LMPD_5_amp_c < 0
UR_LMPD_5_amp_c: -1000.0*BU_LMPD_5_amp_c + 1.0*LMPD_5_amp_c_reverse_8e974 < 0
G_LMPD_6_arg-L_c: 1.0*DGo_LMPD_6_arg-L_c - 1.0*DG_LMPD_6_arg-L_c + 1.776268129875*LC_adp_c + 0.592089376625*LC_adp_m + 0.592089376625*LC_akg_m + 0.592089376625*LC_amp_c + 0.592089376625*LC_arg-L_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c - 0.592089376625*LC_atp_m + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_m - 0.592089376625*LC_hco3_c + 0.592089376625*LC_nadp_m - 0.592089376625*LC_nadph_m - 0.592089376625*LC_nh4_c + 1.776268129875*LC_pi_c + 0.592089376625*LC_pi_m + 0.592089376625*LC_ppi_c = 0
FU_LMPD_6_arg-L_c: 1.0*DG_LMPD_6_arg-L_c + 1000.0*FU_LMPD_6_arg-L_c < 999.999999
BU_LMPD_6_arg-L_c: -1.0*DG_LMPD_6_arg-L_c + 1000.0*BU_LMPD_6_arg-L_c < 999.999999
SU_LMPD_6_arg-L_c: 0 < 1.0*BU_LMPD_6_arg-L_c + 1.0*FU_LMPD_6_arg-L_c < 1
UF_LMPD_6_arg-L_c: -1000.0*FU_LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c < 0
UR_LMPD_6_arg-L_c: -1000.0*BU_LMPD_6_arg-L_c + 1.0*LMPD_6_arg-L_c_reverse_33f31 < 0
G_LMPD_7_asn-L_c: -1.0*DG_LMPD_7_asn-L_c + 0.592089376625*LC_adp_c + 0.592089376625*LC_amp_c + 0.592089376625*LC_asn-L_c - 0.592089376625*LC_asp-L_c - 1.18417875325*LC_atp_c - 0.592089376625*LC_nh4_c + 0.592089376625*LC_pi_c + 0.592089376625*LC_ppi_c + 1.0*DGo_LMPD_7_asn-L_c = 0
FU_LMPD_7_asn-L_c: 1.0*DG_LMPD_7_asn-L_c + 1000.0*FU_LMPD_7_asn-L_c < 999.999999
BU_LMPD_7_asn-L_c: -1.0*DG_LMPD_7_asn-L_c + 1000.0*BU_LMPD_7_asn-L_c < 999.999999
SU_LMPD_7_asn-L_c: 0 < 1.0*BU_LMPD_7_asn-L_c + 1.0*FU_LMPD_7_asn-L_c < 1
UF_LMPD_7_asn-L_c: -1000.0*FU_LMPD_7_asn-L_c + 1.0*LMPD_7_asn-L_c < 0
UR_LMPD_7_asn-L_c: -1000.0*BU_LMPD_7_asn-L_c + 1.0*LMPD_7_asn-L_c_reverse_b41c7 < 0
G_LMPD_8_cmp_c: -1.0*DG_LMPD_8_cmp_c + 2.960446883125*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 3.55253625975*LC_atp_c + 0.592089376625*LC_cmp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_hco3_c - 1.18417875325*LC_nh4_c + 2.960446883125*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_8_cmp_c = 0
FU_LMPD_8_cmp_c: 1.0*DG_LMPD_8_cmp_c + 1000.0*FU_LMPD_8_cmp_c < 999.999999
BU_LMPD_8_cmp_c: -1.0*DG_LMPD_8_cmp_c + 1000.0*BU_LMPD_8_cmp_c < 999.999999
SU_LMPD_8_cmp_c: 0 < 1.0*BU_LMPD_8_cmp_c + 1.0*FU_LMPD_8_cmp_c < 1
UF_LMPD_8_cmp_c: -1000.0*FU_LMPD_8_cmp_c + 1.0*LMPD_8_cmp_c < 0
UR_LMPD_8_cmp_c: -1000.0*BU_LMPD_8_cmp_c + 1.0*LMPD_8_cmp_c_reverse_23cbb < 0
G_LMPD_9_cmp_c: -1.0*DG_LMPD_9_cmp_c + 2.960446883125*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 3.55253625975*LC_atp_c + 0.592089376625*LC_cmp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_fum_c - 0.592089376625*LC_hco3_c - 1.18417875325*LC_nh4_c + 2.960446883125*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_succ_c + 1.0*DGo_LMPD_9_cmp_c = 0
FU_LMPD_9_cmp_c: 1.0*DG_LMPD_9_cmp_c + 1000.0*FU_LMPD_9_cmp_c < 999.999999
BU_LMPD_9_cmp_c: -1.0*DG_LMPD_9_cmp_c + 1000.0*BU_LMPD_9_cmp_c < 999.999999
SU_LMPD_9_cmp_c: 0 < 1.0*BU_LMPD_9_cmp_c + 1.0*FU_LMPD_9_cmp_c < 1
UF_LMPD_9_cmp_c: -1000.0*FU_LMPD_9_cmp_c + 1.0*LMPD_9_cmp_c < 0
UR_LMPD_9_cmp_c: -1000.0*BU_LMPD_9_cmp_c + 1.0*LMPD_9_cmp_c_reverse_85279 < 0
G_LMPD_10_cmp_c: -1.0*DG_LMPD_10_cmp_c + 2.3683575065*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.960446883125*LC_atp_c + 0.592089376625*LC_cmp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_hco3_c - 1.18417875325*LC_nh4_c + 2.3683575065*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_10_cmp_c = 0
FU_LMPD_10_cmp_c: 1.0*DG_LMPD_10_cmp_c + 1000.0*FU_LMPD_10_cmp_c < 999.999999
BU_LMPD_10_cmp_c: -1.0*DG_LMPD_10_cmp_c + 1000.0*BU_LMPD_10_cmp_c < 999.999999
SU_LMPD_10_cmp_c: 0 < 1.0*BU_LMPD_10_cmp_c + 1.0*FU_LMPD_10_cmp_c < 1
UF_LMPD_10_cmp_c: -1000.0*FU_LMPD_10_cmp_c + 1.0*LMPD_10_cmp_c < 0
UR_LMPD_10_cmp_c: -1000.0*BU_LMPD_10_cmp_c + 1.0*LMPD_10_cmp_c_reverse_3acfa < 0
G_LMPD_11_cmp_c: -1.0*DG_LMPD_11_cmp_c + 2.3683575065*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.960446883125*LC_atp_c + 0.592089376625*LC_cmp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_fum_c - 0.592089376625*LC_hco3_c - 1.18417875325*LC_nh4_c + 2.3683575065*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_succ_c + 1.0*DGo_LMPD_11_cmp_c = 0
FU_LMPD_11_cmp_c: 1.0*DG_LMPD_11_cmp_c + 1000.0*FU_LMPD_11_cmp_c < 999.999999
BU_LMPD_11_cmp_c: -1.0*DG_LMPD_11_cmp_c + 1000.0*BU_LMPD_11_cmp_c < 999.999999
SU_LMPD_11_cmp_c: 0 < 1.0*BU_LMPD_11_cmp_c + 1.0*FU_LMPD_11_cmp_c < 1
UF_LMPD_11_cmp_c: -1000.0*FU_LMPD_11_cmp_c + 1.0*LMPD_11_cmp_c < 0
UR_LMPD_11_cmp_c: -1000.0*BU_LMPD_11_cmp_c + 1.0*LMPD_11_cmp_c_reverse_98fa8 < 0
G_LMPD_12_cys-L_c: -1.0*DG_LMPD_12_cys-L_c - 0.592089376625*LC_3pg_c + 0.592089376625*LC_ac_c - 0.592089376625*LC_accoa_c + 0.592089376625*LC_adp_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_amp_c - 1.18417875325*LC_atp_c + 0.592089376625*LC_coa_c + 0.592089376625*LC_cys-L_c - 0.592089376625*LC_glu-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 2.3683575065*LC_nadp_c - 2.3683575065*LC_nadph_c + 1.18417875325*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_so4_c + 1.0*DGo_LMPD_12_cys-L_c = 0
FU_LMPD_12_cys-L_c: 1.0*DG_LMPD_12_cys-L_c + 1000.0*FU_LMPD_12_cys-L_c < 999.999999
BU_LMPD_12_cys-L_c: -1.0*DG_LMPD_12_cys-L_c + 1000.0*BU_LMPD_12_cys-L_c < 999.999999
SU_LMPD_12_cys-L_c: 0 < 1.0*BU_LMPD_12_cys-L_c + 1.0*FU_LMPD_12_cys-L_c < 1
UF_LMPD_12_cys-L_c: -1000.0*FU_LMPD_12_cys-L_c + 1.0*LMPD_12_cys-L_c < 0
UR_LMPD_12_cys-L_c: -1000.0*BU_LMPD_12_cys-L_c + 1.0*LMPD_12_cys-L_c_reverse_6dc0a < 0
G_LMPD_13_damp_c: -1.0*DG_LMPD_13_damp_c - 1.18417875325*LC_3pg_c + 2.960446883125*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 1.18417875325*LC_asp-L_c - 4.144625636375*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_damp_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 2.3683575065*LC_nad_c + 2.3683575065*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_13_damp_c = 0
FU_LMPD_13_damp_c: 1.0*DG_LMPD_13_damp_c + 1000.0*FU_LMPD_13_damp_c < 999.999999
BU_LMPD_13_damp_c: -1.0*DG_LMPD_13_damp_c + 1000.0*BU_LMPD_13_damp_c < 999.999999
SU_LMPD_13_damp_c: 0 < 1.0*BU_LMPD_13_damp_c + 1.0*FU_LMPD_13_damp_c < 1
UF_LMPD_13_damp_c: -1000.0*FU_LMPD_13_damp_c + 1.0*LMPD_13_damp_c < 0
UR_LMPD_13_damp_c: -1000.0*BU_LMPD_13_damp_c + 1.0*LMPD_13_damp_c_reverse_e961c < 0
G_LMPD_14_damp_c: -1.0*DG_LMPD_14_damp_c - 1.18417875325*LC_3pg_c + 1.776268129875*LC_adp_c + 1.18417875325*LC_akg_c + 1.776268129875*LC_amp_c - 1.18417875325*LC_asp-L_c - 3.55253625975*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_damp_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 2.3683575065*LC_nad_c + 2.3683575065*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_14_damp_c = 0
FU_LMPD_14_damp_c: 1.0*DG_LMPD_14_damp_c + 1000.0*FU_LMPD_14_damp_c < 999.999999
BU_LMPD_14_damp_c: -1.0*DG_LMPD_14_damp_c + 1000.0*BU_LMPD_14_damp_c < 999.999999
SU_LMPD_14_damp_c: 0 < 1.0*BU_LMPD_14_damp_c + 1.0*FU_LMPD_14_damp_c < 1
UF_LMPD_14_damp_c: -1000.0*FU_LMPD_14_damp_c + 1.0*LMPD_14_damp_c < 0
UR_LMPD_14_damp_c: -1000.0*BU_LMPD_14_damp_c + 1.0*LMPD_14_damp_c_reverse_d5737 < 0
G_LMPD_15_damp_c: -1.0*DG_LMPD_15_damp_c - 1.18417875325*LC_3pg_c + 2.960446883125*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 1.18417875325*LC_asp-L_c - 4.144625636375*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_damp_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 1.18417875325*LC_nad_c + 1.18417875325*LC_nadh_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_15_damp_c = 0
FU_LMPD_15_damp_c: 1.0*DG_LMPD_15_damp_c + 1000.0*FU_LMPD_15_damp_c < 999.999999
BU_LMPD_15_damp_c: -1.0*DG_LMPD_15_damp_c + 1000.0*BU_LMPD_15_damp_c < 999.999999
SU_LMPD_15_damp_c: 0 < 1.0*BU_LMPD_15_damp_c + 1.0*FU_LMPD_15_damp_c < 1
UF_LMPD_15_damp_c: -1000.0*FU_LMPD_15_damp_c + 1.0*LMPD_15_damp_c < 0
UR_LMPD_15_damp_c: -1000.0*BU_LMPD_15_damp_c + 1.0*LMPD_15_damp_c_reverse_33cc4 < 0
G_LMPD_16_damp_c: -1.0*DG_LMPD_16_damp_c - 1.18417875325*LC_3pg_c + 1.776268129875*LC_adp_c + 1.18417875325*LC_akg_c + 1.776268129875*LC_amp_c - 1.18417875325*LC_asp-L_c - 3.55253625975*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_damp_c + 1.18417875325*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 1.18417875325*LC_nad_c + 1.18417875325*LC_nadh_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 5.328804389625*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_16_damp_c = 0
FU_LMPD_16_damp_c: 1.0*DG_LMPD_16_damp_c + 1000.0*FU_LMPD_16_damp_c < 999.999999
BU_LMPD_16_damp_c: -1.0*DG_LMPD_16_damp_c + 1000.0*BU_LMPD_16_damp_c < 999.999999
SU_LMPD_16_damp_c: 0 < 1.0*BU_LMPD_16_damp_c + 1.0*FU_LMPD_16_damp_c < 1
UF_LMPD_16_damp_c: -1000.0*FU_LMPD_16_damp_c + 1.0*LMPD_16_damp_c < 0
UR_LMPD_16_damp_c: -1000.0*BU_LMPD_16_damp_c + 1.0*LMPD_16_damp_c_reverse_be99f < 0
G_LMPD_17_dcmp_c: -1.0*DG_LMPD_17_dcmp_c + 1.776268129875*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_co2_c + 0.592089376625*LC_dcmp_c - 0.592089376625*LC_fum_c - 0.592089376625*LC_hco3_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 1.776268129875*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_succ_c + 1.0*DGo_LMPD_17_dcmp_c = 0
FU_LMPD_17_dcmp_c: 1.0*DG_LMPD_17_dcmp_c + 1000.0*FU_LMPD_17_dcmp_c < 999.999999
BU_LMPD_17_dcmp_c: -1.0*DG_LMPD_17_dcmp_c + 1000.0*BU_LMPD_17_dcmp_c < 999.999999
SU_LMPD_17_dcmp_c: 0 < 1.0*BU_LMPD_17_dcmp_c + 1.0*FU_LMPD_17_dcmp_c < 1
UF_LMPD_17_dcmp_c: -1000.0*FU_LMPD_17_dcmp_c + 1.0*LMPD_17_dcmp_c < 0
UR_LMPD_17_dcmp_c: -1000.0*BU_LMPD_17_dcmp_c + 1.0*LMPD_17_dcmp_c_reverse_5c6b0 < 0
G_LMPD_18_dcmp_c: -1.0*DG_LMPD_18_dcmp_c + 1.776268129875*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_co2_c + 0.592089376625*LC_dcmp_c - 0.592089376625*LC_hco3_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 1.776268129875*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_18_dcmp_c = 0
FU_LMPD_18_dcmp_c: 1.0*DG_LMPD_18_dcmp_c + 1000.0*FU_LMPD_18_dcmp_c < 999.999999
BU_LMPD_18_dcmp_c: -1.0*DG_LMPD_18_dcmp_c + 1000.0*BU_LMPD_18_dcmp_c < 999.999999
SU_LMPD_18_dcmp_c: 0 < 1.0*BU_LMPD_18_dcmp_c + 1.0*FU_LMPD_18_dcmp_c < 1
UF_LMPD_18_dcmp_c: -1000.0*FU_LMPD_18_dcmp_c + 1.0*LMPD_18_dcmp_c < 0
UR_LMPD_18_dcmp_c: -1000.0*BU_LMPD_18_dcmp_c + 1.0*LMPD_18_dcmp_c_reverse_15021 < 0
G_LMPD_19_dgmp_c: -1.0*DG_LMPD_19_dgmp_c - 1.18417875325*LC_3pg_c + 4.736715013*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 5.92089376625*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_dgmp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 2.960446883125*LC_nad_c + 2.960446883125*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.776268129875*LC_nh4_c + 5.92089376625*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_19_dgmp_c = 0
FU_LMPD_19_dgmp_c: 1.0*DG_LMPD_19_dgmp_c + 1000.0*FU_LMPD_19_dgmp_c < 999.999999
BU_LMPD_19_dgmp_c: -1.0*DG_LMPD_19_dgmp_c + 1000.0*BU_LMPD_19_dgmp_c < 999.999999
SU_LMPD_19_dgmp_c: 0 < 1.0*BU_LMPD_19_dgmp_c + 1.0*FU_LMPD_19_dgmp_c < 1
UF_LMPD_19_dgmp_c: -1000.0*FU_LMPD_19_dgmp_c + 1.0*LMPD_19_dgmp_c < 0
UR_LMPD_19_dgmp_c: -1000.0*BU_LMPD_19_dgmp_c + 1.0*LMPD_19_dgmp_c_reverse_7f2b1 < 0
G_LMPD_20_dgmp_c: -1.0*DG_LMPD_20_dgmp_c - 1.18417875325*LC_3pg_c + 4.144625636375*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 5.328804389625*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_dgmp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 2.960446883125*LC_nad_c + 2.960446883125*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.776268129875*LC_nh4_c + 5.328804389625*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_20_dgmp_c = 0
FU_LMPD_20_dgmp_c: 1.0*DG_LMPD_20_dgmp_c + 1000.0*FU_LMPD_20_dgmp_c < 999.999999
BU_LMPD_20_dgmp_c: -1.0*DG_LMPD_20_dgmp_c + 1000.0*BU_LMPD_20_dgmp_c < 999.999999
SU_LMPD_20_dgmp_c: 0 < 1.0*BU_LMPD_20_dgmp_c + 1.0*FU_LMPD_20_dgmp_c < 1
UF_LMPD_20_dgmp_c: -1000.0*FU_LMPD_20_dgmp_c + 1.0*LMPD_20_dgmp_c < 0
UR_LMPD_20_dgmp_c: -1000.0*BU_LMPD_20_dgmp_c + 1.0*LMPD_20_dgmp_c_reverse_e4f3a < 0
G_LMPD_21_dgmp_c: -1.0*DG_LMPD_21_dgmp_c - 1.18417875325*LC_3pg_c + 4.144625636375*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 5.328804389625*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_dgmp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 1.776268129875*LC_nad_c + 1.776268129875*LC_nadh_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c - 1.776268129875*LC_nh4_c + 5.328804389625*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_21_dgmp_c = 0
FU_LMPD_21_dgmp_c: 1.0*DG_LMPD_21_dgmp_c + 1000.0*FU_LMPD_21_dgmp_c < 999.999999
BU_LMPD_21_dgmp_c: -1.0*DG_LMPD_21_dgmp_c + 1000.0*BU_LMPD_21_dgmp_c < 999.999999
SU_LMPD_21_dgmp_c: 0 < 1.0*BU_LMPD_21_dgmp_c + 1.0*FU_LMPD_21_dgmp_c < 1
UF_LMPD_21_dgmp_c: -1000.0*FU_LMPD_21_dgmp_c + 1.0*LMPD_21_dgmp_c < 0
UR_LMPD_21_dgmp_c: -1000.0*BU_LMPD_21_dgmp_c + 1.0*LMPD_21_dgmp_c_reverse_65045 < 0
G_LMPD_22_dgmp_c: -1.0*DG_LMPD_22_dgmp_c - 1.18417875325*LC_3pg_c + 4.736715013*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 5.92089376625*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_dgmp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c - 1.776268129875*LC_nad_c + 1.776268129875*LC_nadh_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c - 1.776268129875*LC_nh4_c + 5.92089376625*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_22_dgmp_c = 0
FU_LMPD_22_dgmp_c: 1.0*DG_LMPD_22_dgmp_c + 1000.0*FU_LMPD_22_dgmp_c < 999.999999
BU_LMPD_22_dgmp_c: -1.0*DG_LMPD_22_dgmp_c + 1000.0*BU_LMPD_22_dgmp_c < 999.999999
SU_LMPD_22_dgmp_c: 0 < 1.0*BU_LMPD_22_dgmp_c + 1.0*FU_LMPD_22_dgmp_c < 1
UF_LMPD_22_dgmp_c: -1000.0*FU_LMPD_22_dgmp_c + 1.0*LMPD_22_dgmp_c < 0
UR_LMPD_22_dgmp_c: -1000.0*BU_LMPD_22_dgmp_c + 1.0*LMPD_22_dgmp_c_reverse_ebc77 < 0
G_LMPD_23_dtmp_c: -1.0*DG_LMPD_23_dtmp_c - 0.592089376625*LC_3pg_c + 1.776268129875*LC_adp_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_co2_c + 0.592089376625*LC_dtmp_c - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_gly_c - 0.592089376625*LC_hco3_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 1.18417875325*LC_nadp_c - 1.18417875325*LC_nadph_c - 0.592089376625*LC_nh4_c + 2.3683575065*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_23_dtmp_c = 0
FU_LMPD_23_dtmp_c: 1.0*DG_LMPD_23_dtmp_c + 1000.0*FU_LMPD_23_dtmp_c < 999.999999
BU_LMPD_23_dtmp_c: -1.0*DG_LMPD_23_dtmp_c + 1000.0*BU_LMPD_23_dtmp_c < 999.999999
SU_LMPD_23_dtmp_c: 0 < 1.0*BU_LMPD_23_dtmp_c + 1.0*FU_LMPD_23_dtmp_c < 1
UF_LMPD_23_dtmp_c: -1000.0*FU_LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c < 0
UR_LMPD_23_dtmp_c: -1000.0*BU_LMPD_23_dtmp_c + 1.0*LMPD_23_dtmp_c_reverse_ed12d < 0
G_LMPD_24_dtmp_c: -1.0*DG_LMPD_24_dtmp_c - 0.592089376625*LC_3pg_c + 1.776268129875*LC_adp_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_co2_c + 0.592089376625*LC_dtmp_c - 0.592089376625*LC_fum_c - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_gly_c - 0.592089376625*LC_hco3_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 1.18417875325*LC_nadp_c - 1.18417875325*LC_nadph_c - 0.592089376625*LC_nh4_c + 2.3683575065*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_succ_c + 1.0*DGo_LMPD_24_dtmp_c = 0
FU_LMPD_24_dtmp_c: 1.0*DG_LMPD_24_dtmp_c + 1000.0*FU_LMPD_24_dtmp_c < 999.999999
BU_LMPD_24_dtmp_c: -1.0*DG_LMPD_24_dtmp_c + 1000.0*BU_LMPD_24_dtmp_c < 999.999999
SU_LMPD_24_dtmp_c: 0 < 1.0*BU_LMPD_24_dtmp_c + 1.0*FU_LMPD_24_dtmp_c < 1
UF_LMPD_24_dtmp_c: -1000.0*FU_LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c < 0
UR_LMPD_24_dtmp_c: -1000.0*BU_LMPD_24_dtmp_c + 1.0*LMPD_24_dtmp_c_reverse_5d241 < 0
G_LMPD_25_ergst_c: -1.0*DG_LMPD_25_ergst_c + 1.0*DGo_LMPD_25_ergst_c - 10.65760877925*LC_accoa_c + 11.8417875325*LC_adp_c + 0.592089376625*LC_amp_c - 12.433876909125*LC_atp_c + 4.736715013*LC_co2_c + 10.65760877925*LC_coa_c + 0.592089376625*LC_ergst_r - 1.18417875325*LC_nad_c + 1.18417875325*LC_nadh_c + 15.986413168875*LC_nadp_c + 1.18417875325*LC_nadp_r - 15.986413168875*LC_nadph_c - 1.18417875325*LC_nadph_r - 5.92089376625*LC_o2_c - 0.592089376625*LC_o2_r + 4.736715013*LC_pi_c + 4.144625636375*LC_ppi_c = 0
FU_LMPD_25_ergst_c: 1.0*DG_LMPD_25_ergst_c + 1000.0*FU_LMPD_25_ergst_c < 999.999999
BU_LMPD_25_ergst_c: -1.0*DG_LMPD_25_ergst_c + 1000.0*BU_LMPD_25_ergst_c < 999.999999
SU_LMPD_25_ergst_c: 0 < 1.0*BU_LMPD_25_ergst_c + 1.0*FU_LMPD_25_ergst_c < 1
UF_LMPD_25_ergst_c: -1000.0*FU_LMPD_25_ergst_c + 1.0*LMPD_25_ergst_c < 0
UR_LMPD_25_ergst_c: -1000.0*BU_LMPD_25_ergst_c + 1.0*LMPD_25_ergst_c_reverse_96537 < 0
G_LMPD_26_ergst_c: -1.0*DG_LMPD_26_ergst_c + 1.0*DGo_LMPD_26_ergst_c - 10.65760877925*LC_accoa_c + 11.8417875325*LC_adp_c + 0.592089376625*LC_amp_c - 12.433876909125*LC_atp_c + 4.736715013*LC_co2_c + 10.65760877925*LC_coa_c + 0.592089376625*LC_ergst_r + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadh_c + 14.210145039*LC_nadp_c + 1.18417875325*LC_nadp_r - 14.210145039*LC_nadph_c - 1.18417875325*LC_nadph_r - 5.92089376625*LC_o2_c - 0.592089376625*LC_o2_r + 4.736715013*LC_pi_c + 4.144625636375*LC_ppi_c = 0
FU_LMPD_26_ergst_c: 1.0*DG_LMPD_26_ergst_c + 1000.0*FU_LMPD_26_ergst_c < 999.999999
BU_LMPD_26_ergst_c: -1.0*DG_LMPD_26_ergst_c + 1000.0*BU_LMPD_26_ergst_c < 999.999999
SU_LMPD_26_ergst_c: 0 < 1.0*BU_LMPD_26_ergst_c + 1.0*FU_LMPD_26_ergst_c < 1
UF_LMPD_26_ergst_c: -1000.0*FU_LMPD_26_ergst_c + 1.0*LMPD_26_ergst_c < 0
UR_LMPD_26_ergst_c: -1000.0*BU_LMPD_26_ergst_c + 1.0*LMPD_26_ergst_c_reverse_a59e0 < 0
G_LMPD_27_gln-L_c: -1.0*DG_LMPD_27_gln-L_c + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_gln-L_c - 0.592089376625*LC_glu-L_c - 0.592089376625*LC_nh4_c + 0.592089376625*LC_pi_c + 1.0*DGo_LMPD_27_gln-L_c = 0
FU_LMPD_27_gln-L_c: 1.0*DG_LMPD_27_gln-L_c + 1000.0*FU_LMPD_27_gln-L_c < 999.999999
BU_LMPD_27_gln-L_c: -1.0*DG_LMPD_27_gln-L_c + 1000.0*BU_LMPD_27_gln-L_c < 999.999999
SU_LMPD_27_gln-L_c: 0 < 1.0*BU_LMPD_27_gln-L_c + 1.0*FU_LMPD_27_gln-L_c < 1
UF_LMPD_27_gln-L_c: -1000.0*FU_LMPD_27_gln-L_c + 1.0*LMPD_27_gln-L_c < 0
UR_LMPD_27_gln-L_c: -1000.0*BU_LMPD_27_gln-L_c + 1.0*LMPD_27_gln-L_c_reverse_77507 < 0
G_LMPD_28_gly_c: -1.0*DG_LMPD_28_gly_c - 0.592089376625*LC_ala-L_c + 0.592089376625*LC_gly_c - 0.592089376625*LC_icit_c + 0.592089376625*LC_pyr_c + 0.592089376625*LC_succ_c + 1.0*DGo_LMPD_28_gly_c = 0
FU_LMPD_28_gly_c: 1.0*DG_LMPD_28_gly_c + 1000.0*FU_LMPD_28_gly_c < 999.999999
BU_LMPD_28_gly_c: -1.0*DG_LMPD_28_gly_c + 1000.0*BU_LMPD_28_gly_c < 999.999999
SU_LMPD_28_gly_c: 0 < 1.0*BU_LMPD_28_gly_c + 1.0*FU_LMPD_28_gly_c < 1
UF_LMPD_28_gly_c: -1000.0*FU_LMPD_28_gly_c + 1.0*LMPD_28_gly_c < 0
UR_LMPD_28_gly_c: -1000.0*BU_LMPD_28_gly_c + 1.0*LMPD_28_gly_c_reverse_77428 < 0
SU_LMPD_29_glycogen_c: 0 < 1.0*BU_LMPD_29_glycogen_c + 1.0*FU_LMPD_29_glycogen_c < 1
UF_LMPD_29_glycogen_c: -1000.0*FU_LMPD_29_glycogen_c + 1.0*LMPD_29_glycogen_c < 0
UR_LMPD_29_glycogen_c: -1000.0*BU_LMPD_29_glycogen_c + 1.0*LMPD_29_glycogen_c_reverse_68cb3 < 0
G_LMPD_30_gmp_c: -1.0*DG_LMPD_30_gmp_c - 1.18417875325*LC_3pg_c + 4.144625636375*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 5.328804389625*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_gmp_c - 2.960446883125*LC_nad_c + 2.960446883125*LC_nadh_c - 1.776268129875*LC_nh4_c + 5.328804389625*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_30_gmp_c = 0
FU_LMPD_30_gmp_c: 1.0*DG_LMPD_30_gmp_c + 1000.0*FU_LMPD_30_gmp_c < 999.999999
BU_LMPD_30_gmp_c: -1.0*DG_LMPD_30_gmp_c + 1000.0*BU_LMPD_30_gmp_c < 999.999999
SU_LMPD_30_gmp_c: 0 < 1.0*BU_LMPD_30_gmp_c + 1.0*FU_LMPD_30_gmp_c < 1
UF_LMPD_30_gmp_c: -1000.0*FU_LMPD_30_gmp_c + 1.0*LMPD_30_gmp_c < 0
UR_LMPD_30_gmp_c: -1000.0*BU_LMPD_30_gmp_c + 1.0*LMPD_30_gmp_c_reverse_1f79a < 0
G_LMPD_31_gmp_c: -1.0*DG_LMPD_31_gmp_c - 1.18417875325*LC_3pg_c + 4.144625636375*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 5.328804389625*LC_atp_c - 0.592089376625*LC_co2_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_gmp_c - 1.776268129875*LC_nad_c + 1.776268129875*LC_nadh_c - 1.18417875325*LC_nadp_c + 1.18417875325*LC_nadph_c - 1.776268129875*LC_nh4_c + 5.328804389625*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_31_gmp_c = 0
FU_LMPD_31_gmp_c: 1.0*DG_LMPD_31_gmp_c + 1000.0*FU_LMPD_31_gmp_c < 999.999999
BU_LMPD_31_gmp_c: -1.0*DG_LMPD_31_gmp_c + 1000.0*BU_LMPD_31_gmp_c < 999.999999
SU_LMPD_31_gmp_c: 0 < 1.0*BU_LMPD_31_gmp_c + 1.0*FU_LMPD_31_gmp_c < 1
UF_LMPD_31_gmp_c: -1000.0*FU_LMPD_31_gmp_c + 1.0*LMPD_31_gmp_c < 0
UR_LMPD_31_gmp_c: -1000.0*BU_LMPD_31_gmp_c + 1.0*LMPD_31_gmp_c_reverse_d0c02 < 0
G_LMPD_32_his-L_c: -1.0*DG_LMPD_32_his-L_c - 0.592089376625*LC_3pg_c + 1.18417875325*LC_akg_c + 1.776268129875*LC_amp_c - 0.592089376625*LC_asp-L_c - 1.776268129875*LC_atp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_his-L_c - 1.776268129875*LC_nad_c + 1.776268129875*LC_nadh_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c - 0.592089376625*LC_nh4_c + 2.3683575065*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_32_his-L_c = 0
FU_LMPD_32_his-L_c: 1.0*DG_LMPD_32_his-L_c + 1000.0*FU_LMPD_32_his-L_c < 999.999999
BU_LMPD_32_his-L_c: -1.0*DG_LMPD_32_his-L_c + 1000.0*BU_LMPD_32_his-L_c < 999.999999
SU_LMPD_32_his-L_c: 0 < 1.0*BU_LMPD_32_his-L_c + 1.0*FU_LMPD_32_his-L_c < 1
UF_LMPD_32_his-L_c: -1000.0*FU_LMPD_32_his-L_c + 1.0*LMPD_32_his-L_c < 0
UR_LMPD_32_his-L_c: -1000.0*BU_LMPD_32_his-L_c + 1.0*LMPD_32_his-L_c_reverse_c5967 < 0
G_LMPD_33_his-L_c: -1.0*DG_LMPD_33_his-L_c - 0.592089376625*LC_3pg_c + 1.18417875325*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_his-L_c - 2.3683575065*LC_nad_c + 2.3683575065*LC_nadh_c - 0.592089376625*LC_nh4_c + 2.3683575065*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_33_his-L_c = 0
FU_LMPD_33_his-L_c: 1.0*DG_LMPD_33_his-L_c + 1000.0*FU_LMPD_33_his-L_c < 999.999999
BU_LMPD_33_his-L_c: -1.0*DG_LMPD_33_his-L_c + 1000.0*BU_LMPD_33_his-L_c < 999.999999
SU_LMPD_33_his-L_c: 0 < 1.0*BU_LMPD_33_his-L_c + 1.0*FU_LMPD_33_his-L_c < 1
UF_LMPD_33_his-L_c: -1000.0*FU_LMPD_33_his-L_c + 1.0*LMPD_33_his-L_c < 0
UR_LMPD_33_his-L_c: -1000.0*BU_LMPD_33_his-L_c + 1.0*LMPD_33_his-L_c_reverse_0eb8b < 0
G_LMPD_34_his-L_c: -1.0*DG_LMPD_34_his-L_c - 0.592089376625*LC_3pg_c + 1.18417875325*LC_adp_c + 1.18417875325*LC_akg_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_his-L_c - 1.776268129875*LC_nad_c + 1.776268129875*LC_nadh_c - 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadph_c - 0.592089376625*LC_nh4_c + 2.3683575065*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_34_his-L_c = 0
FU_LMPD_34_his-L_c: 1.0*DG_LMPD_34_his-L_c + 1000.0*FU_LMPD_34_his-L_c < 999.999999
BU_LMPD_34_his-L_c: -1.0*DG_LMPD_34_his-L_c + 1000.0*BU_LMPD_34_his-L_c < 999.999999
SU_LMPD_34_his-L_c: 0 < 1.0*BU_LMPD_34_his-L_c + 1.0*FU_LMPD_34_his-L_c < 1
UF_LMPD_34_his-L_c: -1000.0*FU_LMPD_34_his-L_c + 1.0*LMPD_34_his-L_c < 0
UR_LMPD_34_his-L_c: -1000.0*BU_LMPD_34_his-L_c + 1.0*LMPD_34_his-L_c_reverse_e16a9 < 0
G_LMPD_35_his-L_c: -1.0*DG_LMPD_35_his-L_c - 0.592089376625*LC_3pg_c + 1.18417875325*LC_akg_c + 1.776268129875*LC_amp_c - 0.592089376625*LC_asp-L_c - 1.776268129875*LC_atp_c + 0.592089376625*LC_fum_c - 1.18417875325*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_his-L_c - 2.3683575065*LC_nad_c + 2.3683575065*LC_nadh_c - 0.592089376625*LC_nh4_c + 2.3683575065*LC_pi_c + 1.18417875325*LC_ppi_c - 0.592089376625*LC_r5p_c + 1.0*DGo_LMPD_35_his-L_c = 0
FU_LMPD_35_his-L_c: 1.0*DG_LMPD_35_his-L_c + 1000.0*FU_LMPD_35_his-L_c < 999.999999
BU_LMPD_35_his-L_c: -1.0*DG_LMPD_35_his-L_c + 1000.0*BU_LMPD_35_his-L_c < 999.999999
SU_LMPD_35_his-L_c: 0 < 1.0*BU_LMPD_35_his-L_c + 1.0*FU_LMPD_35_his-L_c < 1
UF_LMPD_35_his-L_c: -1000.0*FU_LMPD_35_his-L_c + 1.0*LMPD_35_his-L_c < 0
UR_LMPD_35_his-L_c: -1000.0*BU_LMPD_35_his-L_c + 1.0*LMPD_35_his-L_c_reverse_5e623 < 0
G_LMPD_36_ile-L_c: -1.0*DG_LMPD_36_ile-L_c + 1.18417875325*LC_adp_c + 0.592089376625*LC_akg_m - 0.592089376625*LC_asp-L_c - 1.18417875325*LC_atp_c + 0.592089376625*LC_co2_m - 0.592089376625*LC_glu-L_m + 0.592089376625*LC_ile-L_m + 1.18417875325*LC_nadp_c + 0.592089376625*LC_nadp_m - 1.18417875325*LC_nadph_c - 0.592089376625*LC_nadph_m + 0.592089376625*LC_nh4_m + 1.18417875325*LC_pi_c - 0.592089376625*LC_pyr_m + 0.592089376625*LC_thr-L_c - 0.592089376625*LC_thr-L_m + 1.0*DGo_LMPD_36_ile-L_c - 0.592089376625*LC_h2o_c + 0.592089376625*LC_h2o_m = 0
FU_LMPD_36_ile-L_c: 1.0*DG_LMPD_36_ile-L_c + 1000.0*FU_LMPD_36_ile-L_c < 999.999999
BU_LMPD_36_ile-L_c: -1.0*DG_LMPD_36_ile-L_c + 1000.0*BU_LMPD_36_ile-L_c < 999.999999
SU_LMPD_36_ile-L_c: 0 < 1.0*BU_LMPD_36_ile-L_c + 1.0*FU_LMPD_36_ile-L_c < 1
UF_LMPD_36_ile-L_c: -1000.0*FU_LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c < 0
UR_LMPD_36_ile-L_c: -1000.0*BU_LMPD_36_ile-L_c + 1.0*LMPD_36_ile-L_c_reverse_d2bf7 < 0
G_LMPD_37_ile-L_c: -1.0*DG_LMPD_37_ile-L_c + 1.18417875325*LC_adp_c + 0.592089376625*LC_akg_m - 0.592089376625*LC_asp-L_c - 1.18417875325*LC_atp_c + 0.592089376625*LC_co2_m - 0.592089376625*LC_glu-L_m + 0.592089376625*LC_ile-L_m + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_c + 0.592089376625*LC_nadp_m - 0.592089376625*LC_nadph_c - 0.592089376625*LC_nadph_m + 0.592089376625*LC_nh4_m + 1.18417875325*LC_pi_c - 0.592089376625*LC_pyr_m + 0.592089376625*LC_thr-L_c - 0.592089376625*LC_thr-L_m + 1.0*DGo_LMPD_37_ile-L_c - 0.592089376625*LC_h2o_c + 0.592089376625*LC_h2o_m = 0
FU_LMPD_37_ile-L_c: 1.0*DG_LMPD_37_ile-L_c + 1000.0*FU_LMPD_37_ile-L_c < 999.999999
BU_LMPD_37_ile-L_c: -1.0*DG_LMPD_37_ile-L_c + 1000.0*BU_LMPD_37_ile-L_c < 999.999999
SU_LMPD_37_ile-L_c: 0 < 1.0*BU_LMPD_37_ile-L_c + 1.0*FU_LMPD_37_ile-L_c < 1
UF_LMPD_37_ile-L_c: -1000.0*FU_LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c < 0
UR_LMPD_37_ile-L_c: -1000.0*BU_LMPD_37_ile-L_c + 1.0*LMPD_37_ile-L_c_reverse_1ea52 < 0
G_LMPD_38_leu-L_c: -1.0*DG_LMPD_38_leu-L_c + 1.0*DGo_LMPD_38_leu-L_c - 0.592089376625*LC_accoa_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_co2_c + 0.592089376625*LC_co2_m + 0.592089376625*LC_coa_c - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_leu-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_m - 0.592089376625*LC_nadph_m - 1.18417875325*LC_pyr_m - 0.592089376625*LC_h2o_c + 0.592089376625*LC_h2o_m = 0
FU_LMPD_38_leu-L_c: 1.0*DG_LMPD_38_leu-L_c + 1000.0*FU_LMPD_38_leu-L_c < 999.999999
BU_LMPD_38_leu-L_c: -1.0*DG_LMPD_38_leu-L_c + 1000.0*BU_LMPD_38_leu-L_c < 999.999999
SU_LMPD_38_leu-L_c: 0 < 1.0*BU_LMPD_38_leu-L_c + 1.0*FU_LMPD_38_leu-L_c < 1
UF_LMPD_38_leu-L_c: -1000.0*FU_LMPD_38_leu-L_c + 1.0*LMPD_38_leu-L_c < 0
UR_LMPD_38_leu-L_c: -1000.0*BU_LMPD_38_leu-L_c + 1.0*LMPD_38_leu-L_c_reverse_0e66c < 0
G_LMPD_39_leu-L_c: -1.0*DG_LMPD_39_leu-L_c - 0.592089376625*LC_accoa_m + 0.592089376625*LC_akg_c + 0.592089376625*LC_co2_c + 0.592089376625*LC_co2_m + 0.592089376625*LC_coa_m - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_leu-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_m - 0.592089376625*LC_nadph_m - 1.18417875325*LC_pyr_m + 1.0*DGo_LMPD_39_leu-L_c = 0
FU_LMPD_39_leu-L_c: 1.0*DG_LMPD_39_leu-L_c + 1000.0*FU_LMPD_39_leu-L_c < 999.999999
BU_LMPD_39_leu-L_c: -1.0*DG_LMPD_39_leu-L_c + 1000.0*BU_LMPD_39_leu-L_c < 999.999999
SU_LMPD_39_leu-L_c: 0 < 1.0*BU_LMPD_39_leu-L_c + 1.0*FU_LMPD_39_leu-L_c < 1
UF_LMPD_39_leu-L_c: -1000.0*FU_LMPD_39_leu-L_c + 1.0*LMPD_39_leu-L_c < 0
UR_LMPD_39_leu-L_c: -1000.0*BU_LMPD_39_leu-L_c + 1.0*LMPD_39_leu-L_c_reverse_73152 < 0
G_LMPD_40_lys-L_c: -1.0*DG_LMPD_40_lys-L_c + 1.18417875325*LC_akg_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_co2_m - 1.18417875325*LC_glu-L_c - 0.592089376625*LC_hcit_m + 0.592089376625*LC_lys-L_c - 0.592089376625*LC_nad_c - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadh_m + 1.18417875325*LC_nadp_c - 1.18417875325*LC_nadph_c + 0.592089376625*LC_ppi_c + 1.0*DGo_LMPD_40_lys-L_c = 0
FU_LMPD_40_lys-L_c: 1.0*DG_LMPD_40_lys-L_c + 1000.0*FU_LMPD_40_lys-L_c < 999.999999
BU_LMPD_40_lys-L_c: -1.0*DG_LMPD_40_lys-L_c + 1000.0*BU_LMPD_40_lys-L_c < 999.999999
SU_LMPD_40_lys-L_c: 0 < 1.0*BU_LMPD_40_lys-L_c + 1.0*FU_LMPD_40_lys-L_c < 1
UF_LMPD_40_lys-L_c: -1000.0*FU_LMPD_40_lys-L_c + 1.0*LMPD_40_lys-L_c < 0
UR_LMPD_40_lys-L_c: -1000.0*BU_LMPD_40_lys-L_c + 1.0*LMPD_40_lys-L_c_reverse_4d3a6 < 0
G_LMPD_41_lys-L_c: -1.0*DG_LMPD_41_lys-L_c + 1.18417875325*LC_akg_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_co2_m - 1.18417875325*LC_glu-L_c - 0.592089376625*LC_hcit_m + 0.592089376625*LC_lys-L_c - 0.592089376625*LC_nad_m + 0.592089376625*LC_nadh_m + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c + 0.592089376625*LC_ppi_c + 1.0*DGo_LMPD_41_lys-L_c = 0
FU_LMPD_41_lys-L_c: 1.0*DG_LMPD_41_lys-L_c + 1000.0*FU_LMPD_41_lys-L_c < 999.999999
BU_LMPD_41_lys-L_c: -1.0*DG_LMPD_41_lys-L_c + 1000.0*BU_LMPD_41_lys-L_c < 999.999999
SU_LMPD_41_lys-L_c: 0 < 1.0*BU_LMPD_41_lys-L_c + 1.0*FU_LMPD_41_lys-L_c < 1
UF_LMPD_41_lys-L_c: -1000.0*FU_LMPD_41_lys-L_c + 1.0*LMPD_41_lys-L_c < 0
UR_LMPD_41_lys-L_c: -1000.0*BU_LMPD_41_lys-L_c + 1.0*LMPD_41_lys-L_c_reverse_c021a < 0
SU_LMPD_42_mannan_c: 0 < 1.0*BU_LMPD_42_mannan_c + 1.0*FU_LMPD_42_mannan_c < 1
UF_LMPD_42_mannan_c: -1000.0*FU_LMPD_42_mannan_c + 1.0*LMPD_42_mannan_c < 0
UR_LMPD_42_mannan_c: -1000.0*BU_LMPD_42_mannan_c + 1.0*LMPD_42_mannan_c_reverse_9f626 < 0
SU_LMPD_43_mannan_c: 0 < 1.0*BU_LMPD_43_mannan_c + 1.0*FU_LMPD_43_mannan_c < 1
UF_LMPD_43_mannan_c: -1000.0*FU_LMPD_43_mannan_c + 1.0*LMPD_43_mannan_c < 0
UR_LMPD_43_mannan_c: -1000.0*BU_LMPD_43_mannan_c + 1.0*LMPD_43_mannan_c_reverse_7dd82 < 0
G_LMPD_44_met-L_c: -1.0*DG_LMPD_44_met-L_c - 0.592089376625*LC_3pg_c + 0.592089376625*LC_ac_c - 0.592089376625*LC_accoa_c + 1.18417875325*LC_adp_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 1.776268129875*LC_atp_c + 0.592089376625*LC_coa_c - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_met-L_c + 3.55253625975*LC_nadp_c - 3.55253625975*LC_nadph_c + 1.776268129875*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_so4_c + 1.0*DGo_LMPD_44_met-L_c = 0
FU_LMPD_44_met-L_c: 1.0*DG_LMPD_44_met-L_c + 1000.0*FU_LMPD_44_met-L_c < 999.999999
BU_LMPD_44_met-L_c: -1.0*DG_LMPD_44_met-L_c + 1000.0*BU_LMPD_44_met-L_c < 999.999999
SU_LMPD_44_met-L_c: 0 < 1.0*BU_LMPD_44_met-L_c + 1.0*FU_LMPD_44_met-L_c < 1
UF_LMPD_44_met-L_c: -1000.0*FU_LMPD_44_met-L_c + 1.0*LMPD_44_met-L_c < 0
UR_LMPD_44_met-L_c: -1000.0*BU_LMPD_44_met-L_c + 1.0*LMPD_44_met-L_c_reverse_1646f < 0
G_LMPD_45_met-L_c: -1.0*DG_LMPD_45_met-L_c - 0.592089376625*LC_3pg_c + 0.592089376625*LC_ac_c - 0.592089376625*LC_accoa_c + 1.18417875325*LC_adp_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 1.776268129875*LC_atp_c + 0.592089376625*LC_coa_c - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_gly_c + 0.592089376625*LC_met-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 4.144625636375*LC_nadp_c - 4.144625636375*LC_nadph_c + 1.776268129875*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_so4_c + 1.0*DGo_LMPD_45_met-L_c = 0
FU_LMPD_45_met-L_c: 1.0*DG_LMPD_45_met-L_c + 1000.0*FU_LMPD_45_met-L_c < 999.999999
BU_LMPD_45_met-L_c: -1.0*DG_LMPD_45_met-L_c + 1000.0*BU_LMPD_45_met-L_c < 999.999999
SU_LMPD_45_met-L_c: 0 < 1.0*BU_LMPD_45_met-L_c + 1.0*FU_LMPD_45_met-L_c < 1
UF_LMPD_45_met-L_c: -1000.0*FU_LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c < 0
UR_LMPD_45_met-L_c: -1000.0*BU_LMPD_45_met-L_c + 1.0*LMPD_45_met-L_c_reverse_912f1 < 0
SU_LMPD_47_pa_SC_c: 0 < 1.0*BU_LMPD_47_pa_SC_c + 1.0*FU_LMPD_47_pa_SC_c < 1
UF_LMPD_47_pa_SC_c: -1000.0*FU_LMPD_47_pa_SC_c + 1.0*LMPD_47_pa_SC_c < 0
UR_LMPD_47_pa_SC_c: -1000.0*BU_LMPD_47_pa_SC_c + 1.0*LMPD_47_pa_SC_c_reverse_4049b < 0
SU_LMPD_48_pc_SC_c: 0 < 1.0*BU_LMPD_48_pc_SC_c + 1.0*FU_LMPD_48_pc_SC_c < 1
UF_LMPD_48_pc_SC_c: -1000.0*FU_LMPD_48_pc_SC_c + 1.0*LMPD_48_pc_SC_c < 0
UR_LMPD_48_pc_SC_c: -1000.0*BU_LMPD_48_pc_SC_c + 1.0*LMPD_48_pc_SC_c_reverse_be709 < 0
SU_LMPD_49_pc_SC_c: 0 < 1.0*BU_LMPD_49_pc_SC_c + 1.0*FU_LMPD_49_pc_SC_c < 1
UF_LMPD_49_pc_SC_c: -1000.0*FU_LMPD_49_pc_SC_c + 1.0*LMPD_49_pc_SC_c < 0
UR_LMPD_49_pc_SC_c: -1000.0*BU_LMPD_49_pc_SC_c + 1.0*LMPD_49_pc_SC_c_reverse_55f95 < 0
SU_LMPD_51_pe_SC_c: 0 < 1.0*BU_LMPD_51_pe_SC_c + 1.0*FU_LMPD_51_pe_SC_c < 1
UF_LMPD_51_pe_SC_c: -1000.0*FU_LMPD_51_pe_SC_c + 1.0*LMPD_51_pe_SC_c < 0
UR_LMPD_51_pe_SC_c: -1000.0*BU_LMPD_51_pe_SC_c + 1.0*LMPD_51_pe_SC_c_reverse_09777 < 0
G_LMPD_52_phe-L_c: -1.0*DG_LMPD_52_phe-L_c + 0.592089376625*LC_adp_c + 0.592089376625*LC_akg_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_e4p_c - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_pep_c + 0.592089376625*LC_phe-L_c + 2.3683575065*LC_pi_c + 1.0*DGo_LMPD_52_phe-L_c = 0
FU_LMPD_52_phe-L_c: 1.0*DG_LMPD_52_phe-L_c + 1000.0*FU_LMPD_52_phe-L_c < 999.999999
BU_LMPD_52_phe-L_c: -1.0*DG_LMPD_52_phe-L_c + 1000.0*BU_LMPD_52_phe-L_c < 999.999999
SU_LMPD_52_phe-L_c: 0 < 1.0*BU_LMPD_52_phe-L_c + 1.0*FU_LMPD_52_phe-L_c < 1
UF_LMPD_52_phe-L_c: -1000.0*FU_LMPD_52_phe-L_c + 1.0*LMPD_52_phe-L_c < 0
UR_LMPD_52_phe-L_c: -1000.0*BU_LMPD_52_phe-L_c + 1.0*LMPD_52_phe-L_c_reverse_909a1 < 0
G_LMPD_53_pro-L_c: -1.0*DG_LMPD_53_pro-L_c + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c + 0.592089376625*LC_pi_c + 0.592089376625*LC_pro-L_c + 1.0*DGo_LMPD_53_pro-L_c = 0
FU_LMPD_53_pro-L_c: 1.0*DG_LMPD_53_pro-L_c + 1000.0*FU_LMPD_53_pro-L_c < 999.999999
BU_LMPD_53_pro-L_c: -1.0*DG_LMPD_53_pro-L_c + 1000.0*BU_LMPD_53_pro-L_c < 999.999999
SU_LMPD_53_pro-L_c: 0 < 1.0*BU_LMPD_53_pro-L_c + 1.0*FU_LMPD_53_pro-L_c < 1
UF_LMPD_53_pro-L_c: -1000.0*FU_LMPD_53_pro-L_c + 1.0*LMPD_53_pro-L_c < 0
UR_LMPD_53_pro-L_c: -1000.0*BU_LMPD_53_pro-L_c + 1.0*LMPD_53_pro-L_c_reverse_830d7 < 0
G_LMPD_54_pro-L_c: -1.0*DG_LMPD_54_pro-L_c + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c - 0.592089376625*LC_glu-L_c + 1.18417875325*LC_nadp_c - 1.18417875325*LC_nadph_c + 0.592089376625*LC_pi_c + 0.592089376625*LC_pro-L_c + 1.0*DGo_LMPD_54_pro-L_c = 0
FU_LMPD_54_pro-L_c: 1.0*DG_LMPD_54_pro-L_c + 1000.0*FU_LMPD_54_pro-L_c < 999.999999
BU_LMPD_54_pro-L_c: -1.0*DG_LMPD_54_pro-L_c + 1000.0*BU_LMPD_54_pro-L_c < 999.999999
SU_LMPD_54_pro-L_c: 0 < 1.0*BU_LMPD_54_pro-L_c + 1.0*FU_LMPD_54_pro-L_c < 1
UF_LMPD_54_pro-L_c: -1000.0*FU_LMPD_54_pro-L_c + 1.0*LMPD_54_pro-L_c < 0
UR_LMPD_54_pro-L_c: -1000.0*BU_LMPD_54_pro-L_c + 1.0*LMPD_54_pro-L_c_reverse_a6168 < 0
SU_LMPD_56_ps_SC_c: 0 < 1.0*BU_LMPD_56_ps_SC_c + 1.0*FU_LMPD_56_ps_SC_c < 1
UF_LMPD_56_ps_SC_c: -1000.0*FU_LMPD_56_ps_SC_c + 1.0*LMPD_56_ps_SC_c < 0
UR_LMPD_56_ps_SC_c: -1000.0*BU_LMPD_56_ps_SC_c + 1.0*LMPD_56_ps_SC_c_reverse_7e9d7 < 0
SU_LMPD_57_ptd1ino_SC_c: 0 < 1.0*BU_LMPD_57_ptd1ino_SC_c + 1.0*FU_LMPD_57_ptd1ino_SC_c < 1
UF_LMPD_57_ptd1ino_SC_c: -1000.0*FU_LMPD_57_ptd1ino_SC_c + 1.0*LMPD_57_ptd1ino_SC_c < 0
UR_LMPD_57_ptd1ino_SC_c: -1000.0*BU_LMPD_57_ptd1ino_SC_c + 1.0*LMPD_57_ptd1ino_SC_c_reverse_65089 < 0
SU_LMPD_58_ptd1ino_SC_c: 0 < 1.0*BU_LMPD_58_ptd1ino_SC_c + 1.0*FU_LMPD_58_ptd1ino_SC_c < 1
UF_LMPD_58_ptd1ino_SC_c: -1000.0*FU_LMPD_58_ptd1ino_SC_c + 1.0*LMPD_58_ptd1ino_SC_c < 0
UR_LMPD_58_ptd1ino_SC_c: -1000.0*BU_LMPD_58_ptd1ino_SC_c + 1.0*LMPD_58_ptd1ino_SC_c_reverse_d40c6 < 0
G_LMPD_59_ribflv_c: -1.0*DG_LMPD_59_ribflv_c + 1.0*DGo_LMPD_59_ribflv_c + 5.328804389625*LC_adp_c - 0.592089376625*LC_ala-L_c + 1.776268129875*LC_amp_c - 0.592089376625*LC_asp-L_c - 7.1050725195*LC_atp_c - 0.592089376625*LC_co2_c - 1.18417875325*LC_db4p_c - 0.592089376625*LC_for_c + 0.592089376625*LC_fum_c - 0.592089376625*LC_icit_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 7.1050725195*LC_pi_c + 1.776268129875*LC_ppi_c + 0.592089376625*LC_pyr_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_ribflv_c + 0.592089376625*LC_succ_c = 0
FU_LMPD_59_ribflv_c: 1.0*DG_LMPD_59_ribflv_c + 1000.0*FU_LMPD_59_ribflv_c < 999.999999
BU_LMPD_59_ribflv_c: -1.0*DG_LMPD_59_ribflv_c + 1000.0*BU_LMPD_59_ribflv_c < 999.999999
SU_LMPD_59_ribflv_c: 0 < 1.0*BU_LMPD_59_ribflv_c + 1.0*FU_LMPD_59_ribflv_c < 1
UF_LMPD_59_ribflv_c: -1000.0*FU_LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c < 0
UR_LMPD_59_ribflv_c: -1000.0*BU_LMPD_59_ribflv_c + 1.0*LMPD_59_ribflv_c_reverse_01d9b < 0
G_LMPD_60_ribflv_c: -1.0*DG_LMPD_60_ribflv_c + 1.0*DGo_LMPD_60_ribflv_c + 6.512983142875*LC_adp_c - 0.592089376625*LC_ala-L_c + 1.18417875325*LC_amp_c - 0.592089376625*LC_asp-L_c - 7.697161896125*LC_atp_c - 0.592089376625*LC_co2_c - 1.18417875325*LC_db4p_c - 0.592089376625*LC_for_c + 0.592089376625*LC_fum_c - 0.592089376625*LC_icit_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_nh4_c + 7.1050725195*LC_pi_c + 1.776268129875*LC_ppi_c + 0.592089376625*LC_pyr_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_ribflv_c + 0.592089376625*LC_succ_c = 0
FU_LMPD_60_ribflv_c: 1.0*DG_LMPD_60_ribflv_c + 1000.0*FU_LMPD_60_ribflv_c < 999.999999
BU_LMPD_60_ribflv_c: -1.0*DG_LMPD_60_ribflv_c + 1000.0*BU_LMPD_60_ribflv_c < 999.999999
SU_LMPD_60_ribflv_c: 0 < 1.0*BU_LMPD_60_ribflv_c + 1.0*FU_LMPD_60_ribflv_c < 1
UF_LMPD_60_ribflv_c: -1000.0*FU_LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c < 0
UR_LMPD_60_ribflv_c: -1000.0*BU_LMPD_60_ribflv_c + 1.0*LMPD_60_ribflv_c_reverse_b7a89 < 0
G_LMPD_61_ser-L_c: -1.0*DG_LMPD_61_ser-L_c - 0.592089376625*LC_3pg_c + 0.592089376625*LC_akg_c - 0.592089376625*LC_glu-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_pi_c + 0.592089376625*LC_ser-L_c + 1.0*DGo_LMPD_61_ser-L_c = 0
FU_LMPD_61_ser-L_c: 1.0*DG_LMPD_61_ser-L_c + 1000.0*FU_LMPD_61_ser-L_c < 999.999999
BU_LMPD_61_ser-L_c: -1.0*DG_LMPD_61_ser-L_c + 1000.0*BU_LMPD_61_ser-L_c < 999.999999
SU_LMPD_61_ser-L_c: 0 < 1.0*BU_LMPD_61_ser-L_c + 1.0*FU_LMPD_61_ser-L_c < 1
UF_LMPD_61_ser-L_c: -1000.0*FU_LMPD_61_ser-L_c + 1.0*LMPD_61_ser-L_c < 0
UR_LMPD_61_ser-L_c: -1000.0*BU_LMPD_61_ser-L_c + 1.0*LMPD_61_ser-L_c_reverse_40666 < 0
G_LMPD_62_thr-L_c: -1.0*DG_LMPD_62_thr-L_c + 1.18417875325*LC_adp_c - 0.592089376625*LC_asp-L_c - 1.18417875325*LC_atp_c + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c + 1.18417875325*LC_pi_c + 0.592089376625*LC_thr-L_c + 1.0*DGo_LMPD_62_thr-L_c = 0
FU_LMPD_62_thr-L_c: 1.0*DG_LMPD_62_thr-L_c + 1000.0*FU_LMPD_62_thr-L_c < 999.999999
BU_LMPD_62_thr-L_c: -1.0*DG_LMPD_62_thr-L_c + 1000.0*BU_LMPD_62_thr-L_c < 999.999999
SU_LMPD_62_thr-L_c: 0 < 1.0*BU_LMPD_62_thr-L_c + 1.0*FU_LMPD_62_thr-L_c < 1
UF_LMPD_62_thr-L_c: -1000.0*FU_LMPD_62_thr-L_c + 1.0*LMPD_62_thr-L_c < 0
UR_LMPD_62_thr-L_c: -1000.0*BU_LMPD_62_thr-L_c + 1.0*LMPD_62_thr-L_c_reverse_e4d09 < 0
G_LMPD_63_thr-L_c: -1.0*DG_LMPD_63_thr-L_c + 1.18417875325*LC_adp_c - 0.592089376625*LC_asp-L_c - 1.18417875325*LC_atp_c + 1.18417875325*LC_nadp_c - 1.18417875325*LC_nadph_c + 1.18417875325*LC_pi_c + 0.592089376625*LC_thr-L_c + 1.0*DGo_LMPD_63_thr-L_c = 0
FU_LMPD_63_thr-L_c: 1.0*DG_LMPD_63_thr-L_c + 1000.0*FU_LMPD_63_thr-L_c < 999.999999
BU_LMPD_63_thr-L_c: -1.0*DG_LMPD_63_thr-L_c + 1000.0*BU_LMPD_63_thr-L_c < 999.999999
SU_LMPD_63_thr-L_c: 0 < 1.0*BU_LMPD_63_thr-L_c + 1.0*FU_LMPD_63_thr-L_c < 1
UF_LMPD_63_thr-L_c: -1000.0*FU_LMPD_63_thr-L_c + 1.0*LMPD_63_thr-L_c < 0
UR_LMPD_63_thr-L_c: -1000.0*BU_LMPD_63_thr-L_c + 1.0*LMPD_63_thr-L_c_reverse_52245 < 0
G_LMPD_64_tre_c: -1.0*DG_LMPD_64_tre_c + 0.592089376625*LC_adp_c - 0.592089376625*LC_atp_c - 0.592089376625*LC_g1p_c - 0.592089376625*LC_g6p_c + 0.592089376625*LC_pi_c + 0.592089376625*LC_ppi_c + 0.592089376625*LC_tre_c + 1.0*DGo_LMPD_64_tre_c = 0
FU_LMPD_64_tre_c: 1.0*DG_LMPD_64_tre_c + 1000.0*FU_LMPD_64_tre_c < 999.999999
BU_LMPD_64_tre_c: -1.0*DG_LMPD_64_tre_c + 1000.0*BU_LMPD_64_tre_c < 999.999999
SU_LMPD_64_tre_c: 0 < 1.0*BU_LMPD_64_tre_c + 1.0*FU_LMPD_64_tre_c < 1
UF_LMPD_64_tre_c: -1000.0*FU_LMPD_64_tre_c + 1.0*LMPD_64_tre_c < 0
UR_LMPD_64_tre_c: -1000.0*BU_LMPD_64_tre_c + 1.0*LMPD_64_tre_c_reverse_64f54 < 0
SU_LMPD_65_triglyc_SC_c: 0 < 1.0*BU_LMPD_65_triglyc_SC_c + 1.0*FU_LMPD_65_triglyc_SC_c < 1
UF_LMPD_65_triglyc_SC_c: -1000.0*FU_LMPD_65_triglyc_SC_c + 1.0*LMPD_65_triglyc_SC_c < 0
UR_LMPD_65_triglyc_SC_c: -1000.0*BU_LMPD_65_triglyc_SC_c + 1.0*LMPD_65_triglyc_SC_c_reverse_7f8d4 < 0
SU_LMPD_66_triglyc_SC_c: 0 < 1.0*BU_LMPD_66_triglyc_SC_c + 1.0*FU_LMPD_66_triglyc_SC_c < 1
UF_LMPD_66_triglyc_SC_c: -1000.0*FU_LMPD_66_triglyc_SC_c + 1.0*LMPD_66_triglyc_SC_c < 0
UR_LMPD_66_triglyc_SC_c: -1000.0*BU_LMPD_66_triglyc_SC_c + 1.0*LMPD_66_triglyc_SC_c_reverse_42d98 < 0
G_LMPD_67_trp-L_c: -1.0*DG_LMPD_67_trp-L_c + 1.0*DGo_LMPD_67_trp-L_c - 0.592089376625*LC_3pg_c + 1.18417875325*LC_adp_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_amp_c - 1.776268129875*LC_atp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_e4p_c + 0.592089376625*LC_g3p_c - 0.592089376625*LC_glu-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 0.592089376625*LC_nh4_c - 1.18417875325*LC_pep_c + 3.55253625975*LC_pi_c + 0.592089376625*LC_ppi_c + 0.592089376625*LC_pyr_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_trp-L_c = 0
FU_LMPD_67_trp-L_c: 1.0*DG_LMPD_67_trp-L_c + 1000.0*FU_LMPD_67_trp-L_c < 999.999999
BU_LMPD_67_trp-L_c: -1.0*DG_LMPD_67_trp-L_c + 1000.0*BU_LMPD_67_trp-L_c < 999.999999
SU_LMPD_67_trp-L_c: 0 < 1.0*BU_LMPD_67_trp-L_c + 1.0*FU_LMPD_67_trp-L_c < 1
UF_LMPD_67_trp-L_c: -1000.0*FU_LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c < 0
UR_LMPD_67_trp-L_c: -1000.0*BU_LMPD_67_trp-L_c + 1.0*LMPD_67_trp-L_c_reverse_17a73 < 0
G_LMPD_68_tyr-L_c: -1.0*DG_LMPD_68_tyr-L_c + 0.592089376625*LC_adp_c + 0.592089376625*LC_akg_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_e4p_c - 0.592089376625*LC_glu-L_c - 0.592089376625*LC_nad_c + 0.592089376625*LC_nadh_c + 0.592089376625*LC_nadp_c - 0.592089376625*LC_nadph_c - 1.18417875325*LC_pep_c + 2.3683575065*LC_pi_c + 0.592089376625*LC_tyr-L_c + 1.0*DGo_LMPD_68_tyr-L_c = 0
FU_LMPD_68_tyr-L_c: 1.0*DG_LMPD_68_tyr-L_c + 1000.0*FU_LMPD_68_tyr-L_c < 999.999999
BU_LMPD_68_tyr-L_c: -1.0*DG_LMPD_68_tyr-L_c + 1000.0*BU_LMPD_68_tyr-L_c < 999.999999
SU_LMPD_68_tyr-L_c: 0 < 1.0*BU_LMPD_68_tyr-L_c + 1.0*FU_LMPD_68_tyr-L_c < 1
UF_LMPD_68_tyr-L_c: -1000.0*FU_LMPD_68_tyr-L_c + 1.0*LMPD_68_tyr-L_c < 0
UR_LMPD_68_tyr-L_c: -1000.0*BU_LMPD_68_tyr-L_c + 1.0*LMPD_68_tyr-L_c_reverse_6983e < 0
G_LMPD_69_tyr-L_c: -1.0*DG_LMPD_69_tyr-L_c + 0.592089376625*LC_adp_c + 0.592089376625*LC_akg_c - 0.592089376625*LC_atp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_e4p_c - 0.592089376625*LC_glu-L_c - 1.18417875325*LC_pep_c + 2.3683575065*LC_pi_c + 0.592089376625*LC_tyr-L_c + 1.0*DGo_LMPD_69_tyr-L_c = 0
FU_LMPD_69_tyr-L_c: 1.0*DG_LMPD_69_tyr-L_c + 1000.0*FU_LMPD_69_tyr-L_c < 999.999999
BU_LMPD_69_tyr-L_c: -1.0*DG_LMPD_69_tyr-L_c + 1000.0*BU_LMPD_69_tyr-L_c < 999.999999
SU_LMPD_69_tyr-L_c: 0 < 1.0*BU_LMPD_69_tyr-L_c + 1.0*FU_LMPD_69_tyr-L_c < 1
UF_LMPD_69_tyr-L_c: -1000.0*FU_LMPD_69_tyr-L_c + 1.0*LMPD_69_tyr-L_c < 0
UR_LMPD_69_tyr-L_c: -1000.0*BU_LMPD_69_tyr-L_c + 1.0*LMPD_69_tyr-L_c_reverse_03f73 < 0
G_LMPD_70_ump_c: -1.0*DG_LMPD_70_ump_c + 1.776268129875*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_fum_c - 0.592089376625*LC_hco3_c - 0.592089376625*LC_nh4_c + 1.776268129875*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_r5p_c + 0.592089376625*LC_succ_c + 0.592089376625*LC_ump_c + 1.0*DGo_LMPD_70_ump_c = 0
FU_LMPD_70_ump_c: 1.0*DG_LMPD_70_ump_c + 1000.0*FU_LMPD_70_ump_c < 999.999999
BU_LMPD_70_ump_c: -1.0*DG_LMPD_70_ump_c + 1000.0*BU_LMPD_70_ump_c < 999.999999
SU_LMPD_70_ump_c: 0 < 1.0*BU_LMPD_70_ump_c + 1.0*FU_LMPD_70_ump_c < 1
UF_LMPD_70_ump_c: -1000.0*FU_LMPD_70_ump_c + 1.0*LMPD_70_ump_c < 0
UR_LMPD_70_ump_c: -1000.0*BU_LMPD_70_ump_c + 1.0*LMPD_70_ump_c_reverse_3dcab < 0
G_LMPD_71_ump_c: -1.0*DG_LMPD_71_ump_c + 1.776268129875*LC_adp_c + 0.592089376625*LC_amp_c - 0.592089376625*LC_asp-L_c - 2.3683575065*LC_atp_c + 0.592089376625*LC_co2_c - 0.592089376625*LC_hco3_c - 0.592089376625*LC_nh4_c + 1.776268129875*LC_pi_c + 0.592089376625*LC_ppi_c - 0.592089376625*LC_q6_m + 0.592089376625*LC_q6h2_m - 0.592089376625*LC_r5p_c + 0.592089376625*LC_ump_c + 1.0*DGo_LMPD_71_ump_c = 0
FU_LMPD_71_ump_c: 1.0*DG_LMPD_71_ump_c + 1000.0*FU_LMPD_71_ump_c < 999.999999
BU_LMPD_71_ump_c: -1.0*DG_LMPD_71_ump_c + 1000.0*BU_LMPD_71_ump_c < 999.999999
SU_LMPD_71_ump_c: 0 < 1.0*BU_LMPD_71_ump_c + 1.0*FU_LMPD_71_ump_c < 1
UF_LMPD_71_ump_c: -1000.0*FU_LMPD_71_ump_c + 1.0*LMPD_71_ump_c < 0
UR_LMPD_71_ump_c: -1000.0*BU_LMPD_71_ump_c + 1.0*LMPD_71_ump_c_reverse_dc9d0 < 0
G_LMPD_72_val-L_c: -1.0*DG_LMPD_72_val-L_c + 0.592089376625*LC_akg_c + 0.592089376625*LC_co2_m - 0.592089376625*LC_glu-L_c + 0.592089376625*LC_nadp_m - 0.592089376625*LC_nadph_m - 1.18417875325*LC_pyr_m + 0.592089376625*LC_val-L_c + 1.0*DGo_LMPD_72_val-L_c = 0
FU_LMPD_72_val-L_c: 1.0*DG_LMPD_72_val-L_c + 1000.0*FU_LMPD_72_val-L_c < 999.999999
BU_LMPD_72_val-L_c: -1.0*DG_LMPD_72_val-L_c + 1000.0*BU_LMPD_72_val-L_c < 999.999999
SU_LMPD_72_val-L_c: 0 < 1.0*BU_LMPD_72_val-L_c + 1.0*FU_LMPD_72_val-L_c < 1
UF_LMPD_72_val-L_c: -1000.0*FU_LMPD_72_val-L_c + 1.0*LMPD_72_val-L_c < 0
UR_LMPD_72_val-L_c: -1000.0*BU_LMPD_72_val-L_c + 1.0*LMPD_72_val-L_c_reverse_768ff < 0
G_LMPD_73_zymst_c: -1.0*DG_LMPD_73_zymst_c + 1.0*DGo_LMPD_73_zymst_c - 10.65760877925*LC_accoa_c + 10.65760877925*LC_adp_c - 10.65760877925*LC_atp_c + 4.736715013*LC_co2_c + 10.65760877925*LC_coa_c + 0.592089376625*LC_for_c + 0.592089376625*LC_nad_c - 0.592089376625*LC_nadh_c + 13.02596628575*LC_nadp_c + 0.592089376625*LC_nadp_r - 13.02596628575*LC_nadph_c - 0.592089376625*LC_nadph_r - 5.328804389625*LC_o2_c - 0.592089376625*LC_o2_r + 3.55253625975*LC_pi_c + 3.55253625975*LC_ppi_c + 0.592089376625*LC_zymst_c = 0
FU_LMPD_73_zymst_c: 1.0*DG_LMPD_73_zymst_c + 1000.0*FU_LMPD_73_zymst_c < 999.999999
BU_LMPD_73_zymst_c: -1.0*DG_LMPD_73_zymst_c + 1000.0*BU_LMPD_73_zymst_c < 999.999999
SU_LMPD_73_zymst_c: 0 < 1.0*BU_LMPD_73_zymst_c + 1.0*FU_LMPD_73_zymst_c < 1
UF_LMPD_73_zymst_c: -1000.0*FU_LMPD_73_zymst_c + 1.0*LMPD_73_zymst_c < 0
UR_LMPD_73_zymst_c: -1000.0*BU_LMPD_73_zymst_c + 1.0*LMPD_73_zymst_c_reverse_e6142 < 0
G_LMPD_74_zymst_c: -1.0*DG_LMPD_74_zymst_c + 1.0*DGo_LMPD_74_zymst_c - 10.65760877925*LC_accoa_c + 10.65760877925*LC_adp_c - 10.65760877925*LC_atp_c + 4.736715013*LC_co2_c + 10.65760877925*LC_coa_c + 0.592089376625*LC_for_c - 1.18417875325*LC_nad_c + 1.18417875325*LC_nadh_c + 14.802234415625*LC_nadp_c + 0.592089376625*LC_nadp_r - 14.802234415625*LC_nadph_c - 0.592089376625*LC_nadph_r - 5.328804389625*LC_o2_c - 0.592089376625*LC_o2_r + 3.55253625975*LC_pi_c + 3.55253625975*LC_ppi_c + 0.592089376625*LC_zymst_c = 0
FU_LMPD_74_zymst_c: 1.0*DG_LMPD_74_zymst_c + 1000.0*FU_LMPD_74_zymst_c < 999.999999
BU_LMPD_74_zymst_c: -1.0*DG_LMPD_74_zymst_c + 1000.0*BU_LMPD_74_zymst_c < 999.999999
SU_LMPD_74_zymst_c: 0 < 1.0*BU_LMPD_74_zymst_c + 1.0*FU_LMPD_74_zymst_c < 1
UF_LMPD_74_zymst_c: -1000.0*FU_LMPD_74_zymst_c + 1.0*LMPD_74_zymst_c < 0
UR_LMPD_74_zymst_c: -1000.0*BU_LMPD_74_zymst_c + 1.0*LMPD_74_zymst_c_reverse_84617 < 0
Bounds
0.0 <= biomass_SC5_notrace <= 100.0
0.0 <= biomass_SC5_notrace_reverse_e32ff <= 0.0
0.0 <= EX_4abut_e <= 100.0
0.0 <= EX_4abut_e_reverse_82295 <= 0.0
0.0 <= EX_ac_e <= 100.0
0.0 <= EX_ac_e_reverse_0be96 <= 0.0
0.0 <= EX_acald_e <= 100.0
0.0 <= EX_acald_e_reverse_c096e <= 0.0
0.0 <= EX_akg_e <= 100.0
0.0 <= EX_akg_e_reverse_70d85 <= 0.0
0.0 <= EX_ala-L_e <= 100.0
0.0 <= EX_ala-L_e_reverse_45174 <= 0.0
0.0 <= EX_arg-L_e <= 100.0
0.0 <= EX_arg-L_e_reverse_2768b <= 0.0
0.0 <= EX_asn-L_e <= 100.0
0.0 <= EX_asn-L_e_reverse_73ed2 <= 0.0
0.0 <= EX_asp-L_e <= 100.0
0.0 <= EX_asp-L_e_reverse_ceee6 <= 0.0
0.0 <= EX_cit_e <= 100.0
0.0 <= EX_cit_e_reverse_0835e <= 0.0
0.0 <= EX_co2_e <= 100.0
0.0 <= EX_co2_e_reverse_d0466 <= 0.0
0.0 <= EX_cys-L_e <= 100.0
0.0 <= EX_cys-L_e_reverse_c0cbc <= 0.0
0.0 <= EX_ergst_e <= 100.0
0.0 <= EX_ergst_e_reverse_92efd <= 0.0
0.0 <= EX_etoh_e <= 100.0
0.0 <= EX_etoh_e_reverse_cc64f <= 0.0
0.0 <= EX_for_e <= 100.0
0.0 <= EX_for_e_reverse_23269 <= 0.0
0.0 <= EX_fum_e <= 100.0
0.0 <= EX_fum_e_reverse_e3432 <= 0.0
0.0 <= EX_glc_e <= 100.0
0.0 <= EX_glc_e_reverse_6226d <= 10.0
0.0 <= EX_gln-L_e <= 100.0
0.0 <= EX_gln-L_e_reverse_1be13 <= 0.0
0.0 <= EX_glu-L_e <= 100.0
0.0 <= EX_glu-L_e_reverse_727e6 <= 0.0
0.0 <= EX_gly_e <= 100.0
0.0 <= EX_gly_e_reverse_6956b <= 0.0
0.0 <= EX_h_e <= 100.0
0.0 <= EX_h_e_reverse_3e0c5 <= 100.0
0.0 <= EX_h2o_e <= 100.0
0.0 <= EX_h2o_e_reverse_3ced4 <= 100.0
0.0 <= EX_his-L_e <= 100.0
0.0 <= EX_his-L_e_reverse_e1c17 <= 0.0
0.0 <= EX_ile-L_e <= 100.0
0.0 <= EX_ile-L_e_reverse_e9c53 <= 0.0
0.0 <= EX_leu-L_e <= 100.0
0.0 <= EX_leu-L_e_reverse_9947f <= 0.0
0.0 <= EX_lys-L_e <= 100.0
0.0 <= EX_lys-L_e_reverse_5779d <= 0.0
0.0 <= EX_mal-L_e <= 100.0
0.0 <= EX_mal-L_e_reverse_00553 <= 0.0
0.0 <= EX_met-L_e <= 100.0
0.0 <= EX_met-L_e_reverse_f8e12 <= 0.0
0.0 <= EX_nh4_e <= 100.0
0.0 <= EX_nh4_e_reverse_f9cc6 <= 100.0
0.0 <= EX_o2_e <= 100.0
0.0 <= EX_o2_e_reverse_efa94 <= 2.0
0.0 <= EX_oaa_e <= 100.0
0.0 <= EX_oaa_e_reverse_9b144 <= 0.0
0.0 <= EX_phe-L_e <= 100.0
0.0 <= EX_phe-L_e_reverse_155f7 <= 0.0
0.0 <= EX_pi_e <= 100.0
0.0 <= EX_pi_e_reverse_1fb09 <= 100.0
0.0 <= EX_pro-L_e <= 100.0
0.0 <= EX_pro-L_e_reverse_26581 <= 0.0
0.0 <= EX_pyr_e <= 100.0
0.0 <= EX_pyr_e_reverse_1f6de <= 0.0
0.0 <= EX_ser-L_e <= 100.0
0.0 <= EX_ser-L_e_reverse_6b4e0 <= 0.0
0.0 <= EX_so4_e <= 100.0
0.0 <= EX_so4_e_reverse_5c8ed <= 100.0
0.0 <= EX_succ_e <= 100.0
0.0 <= EX_succ_e_reverse_a9039 <= 0.0
0.0 <= EX_thr-L_e <= 100.0
0.0 <= EX_thr-L_e_reverse_9a854 <= 0.0
0.0 <= EX_trp-L_e <= 100.0
0.0 <= EX_trp-L_e_reverse_273bc <= 0.0
0.0 <= EX_tyr-L_e <= 100.0
0.0 <= EX_tyr-L_e_reverse_08f6e <= 0.0
0.0 <= EX_val-L_e <= 100.0
0.0 <= EX_val-L_e_reverse_4c155 <= 0.0
0.0 <= EX_zymst_e <= 100.0
0.0 <= EX_zymst_e_reverse_78276 <= 0.0
0.0 <= 6PGLter <= 100.0
0.0 <= 6PGLter_reverse_0e388 <= 100.0
0.0 <= ABUTt2r <= 100.0
0.0 <= ABUTt2r_reverse_d925e <= 100.0
0.0 <= ACALDt <= 100.0
0.0 <= ACALDt_reverse_858fa <= 100.0
0.0 <= ACALDtm <= 100.0
0.0 <= ACALDtm_reverse_c2973 <= 100.0
0.0 <= ACOAHim <= 100.0
0.0 <= ACOAHim_reverse_069ca <= 0.0
0.0 <= ACONT <= 100.0
0.0 <= ACONT_reverse_7c2d5 <= 100.0
0.0 <= ACONTm <= 100.0
0.0 <= ACONTm_reverse_8079b <= 100.0
0.0 <= ACS <= 100.0
0.0 <= ACS_reverse_37635 <= 0.0
0.0 <= ACSm <= 100.0
0.0 <= ACSm_reverse_5e300 <= 0.0
0.0 <= ACt2r <= 100.0
0.0 <= ACt2r_reverse_64e05 <= 100.0
0.0 <= ACtm <= 100.0
0.0 <= ACtm_reverse_5029e <= 100.0
0.0 <= ACtr <= 100.0
0.0 <= ACtr_reverse_b9231 <= 100.0
0.0 <= ADK1 <= 100.0
0.0 <= ADK1_reverse_a6f90 <= 100.0
0.0 <= ADK1m <= 100.0
0.0 <= ADK1m_reverse_59ef9 <= 100.0
0.0 <= AKGDam <= 100.0
0.0 <= AKGDam_reverse_a5f55 <= 100.0
0.0 <= AKGDbm <= 100.0
0.0 <= AKGDbm_reverse_7cf54 <= 0.0
0.0 <= AKGMAL <= 100.0
0.0 <= AKGMAL_reverse_76d29 <= 100.0
0.0 <= AKGt2r <= 100.0
0.0 <= AKGt2r_reverse_5d500 <= 100.0
0.0 <= ALATA_L <= 100.0
0.0 <= ALATA_L_reverse_e54ff <= 100.0
0.0 <= ALATA_Lm <= 100.0
0.0 <= ALATA_Lm_reverse_c845b <= 100.0
0.0 <= ALAt2r <= 100.0
0.0 <= ALAt2r_reverse_552b7 <= 100.0
0.0 <= ALAtmi <= 100.0
0.0 <= ALAtmi_reverse_c89ef <= 0.0
0.0 <= ALCD2if <= 100.0
0.0 <= ALCD2if_reverse_85957 <= 0.0
0.0 <= ALCD2ir <= 100.0
0.0 <= ALCD2ir_reverse_ba067 <= 0.0
0.0 <= ALCD2irm <= 100.0
0.0 <= ALCD2irm_reverse_2557e <= 0.0
0.0 <= ALCD2x <= 100.0
0.0 <= ALCD2x_reverse_5d107 <= 100.0
0.0 <= ALDD2xm <= 100.0
0.0 <= ALDD2xm_reverse_c93a0 <= 0.0
0.0 <= ALDD2y <= 100.0
0.0 <= ALDD2y_reverse_03afb <= 0.0
0.0 <= ALDD2ym <= 100.0
0.0 <= ALDD2ym_reverse_cfc40 <= 0.0
0.0 <= ARGt2r <= 100.0
0.0 <= ARGt2r_reverse_57e3f <= 100.0
0.0 <= ASNt2r <= 100.0
0.0 <= ASNt2r_reverse_5bcd2 <= 100.0
0.0 <= ASPGLU2m <= 100.0
0.0 <= ASPGLU2m_reverse_efaea <= 0.0
0.0 <= ASPTA <= 100.0
0.0 <= ASPTA_reverse_36525 <= 100.0
0.0 <= ASPTAm <= 100.0
0.0 <= ASPTAm_reverse_dc3b5 <= 100.0
0.0 <= ASPt2m <= 100.0
0.0 <= ASPt2m_reverse_c9585 <= 100.0
0.0 <= ASPt2r <= 100.0
0.0 <= ASPt2r_reverse_3a681 <= 100.0
0.0 <= ATPH1 <= 100.0
0.0 <= ATPH1_reverse_23c0c <= 0.0
1.0 <= ATPM <= 1.0
0.0 <= ATPM_reverse_5b752 <= 0.0
0.0 <= ATPS <= 100.0
0.0 <= ATPS_reverse_64286 <= 0.0
0.0 <= ATPS3m <= 100.0
0.0 <= ATPS3m_reverse_5d8d3 <= 0.0
0.0 <= ATPtm-H <= 100.0
0.0 <= ATPtm-H_reverse_ab831 <= 0.0
0.0 <= CITt2r <= 100.0
0.0 <= CITt2r_reverse_9f606 <= 100.0
0.0 <= CITtam <= 100.0
0.0 <= CITtam_reverse_5230c <= 100.0
0.0 <= CITtcm <= 100.0
0.0 <= CITtcm_reverse_df725 <= 100.0
0.0 <= CO2t <= 100.0
0.0 <= CO2t_reverse_7c42f <= 100.0
0.0 <= CO2tm <= 100.0
0.0 <= CO2tm_reverse_48c28 <= 100.0
0.0 <= CSm <= 100.0
0.0 <= CSm_reverse_bb3a1 <= 0.0
0.0 <= CYOOm <= 100.0
0.0 <= CYOOm_reverse_d730e <= 0.0
0.0 <= CYOR_u6m <= 100.0
0.0 <= CYOR_u6m_reverse_ba2b8 <= 0.0
0.0 <= CYSt2r <= 100.0
0.0 <= CYSt2r_reverse_6167b <= 100.0
0.0 <= D-LACt2m <= 100.0
0.0 <= D-LACt2m_reverse_bb9a8 <= 100.0
0.0 <= D-LACtm <= 100.0
0.0 <= D-LACtm_reverse_08640 <= 100.0
0.0 <= DB4PS <= 100.0
0.0 <= DB4PS_reverse_43dd1 <= 0.0
0.0 <= ENO <= 100.0
0.0 <= ENO_reverse_40eea <= 100.0
0.0 <= ERGSTt <= 100.0
0.0 <= ERGSTt_reverse_c595b <= 100.0
0.0 <= ERGSTter <= 100.0
0.0 <= ERGSTter_reverse_23bab <= 100.0
0.0 <= ETOHt <= 100.0
0.0 <= ETOHt_reverse_a82a5 <= 100.0
0.0 <= ETOHtm <= 100.0
0.0 <= ETOHtm_reverse_9d4eb <= 100.0
0.0 <= FBA <= 100.0
0.0 <= FBA_reverse_84806 <= 100.0
0.0 <= FBA3 <= 100.0
0.0 <= FBA3_reverse_0d49f <= 100.0
0.0 <= FBP <= 100.0
0.0 <= FBP_reverse_bf2c9 <= 0.0
0.0 <= FDH <= 100.0
0.0 <= FDH_reverse_06346 <= 0.0
0.0 <= FDNG <= 100.0
0.0 <= FDNG_reverse_f3615 <= 0.0
0.0 <= FORt <= 100.0
0.0 <= FORt_reverse_40f9f <= 100.0
0.0 <= FRDcm <= 100.0
0.0 <= FRDcm_reverse_6cbef <= 0.0
0.0 <= FRDm <= 100.0
0.0 <= FRDm_reverse_ef6ae <= 0.0
0.0 <= FUM <= 100.0
0.0 <= FUM_reverse_d3642 <= 100.0
0.0 <= FUMm <= 100.0
0.0 <= FUMm_reverse_5e620 <= 100.0
0.0 <= FUMt2r <= 100.0
0.0 <= FUMt2r_reverse_d93c0 <= 100.0
0.0 <= G6PDH2 <= 100.0
0.0 <= G6PDH2_reverse_47bd4 <= 0.0
0.0 <= G6PDH2er <= 100.0
0.0 <= G6PDH2er_reverse_9f5fd <= 0.0
0.0 <= G6PI <= 100.0
0.0 <= G6PI_reverse_834e6 <= 100.0
0.0 <= G6PI3 <= 100.0
0.0 <= G6PI3_reverse_734c4 <= 100.0
0.0 <= G6Pter <= 100.0
0.0 <= G6Pter_reverse_b4f95 <= 100.0
0.0 <= GAPD <= 100.0
0.0 <= GAPD_reverse_459c1 <= 100.0
0.0 <= GCC2cm <= 100.0
0.0 <= GCC2cm_reverse_e896f <= 100.0
0.0 <= GLCt1 <= 100.0
0.0 <= GLCt1_reverse_3bcc5 <= 0.0
0.0 <= GLNt2r <= 100.0
0.0 <= GLNt2r_reverse_35d46 <= 100.0
0.0 <= GLUDC <= 100.0
0.0 <= GLUDC_reverse_c4aae <= 0.0
0.0 <= GLUDxi <= 100.0
0.0 <= GLUDxi_reverse_6fbe2 <= 0.0
0.0 <= GLUDyi <= 100.0
0.0 <= GLUDyi_reverse_a2351 <= 0.0
0.0 <= GLUK <= 100.0
0.0 <= GLUK_reverse_ecb52 <= 0.0
0.0 <= GLUt2r <= 100.0
0.0 <= GLUt2r_reverse_3e88a <= 100.0
0.0 <= GLUt7m <= 100.0
0.0 <= GLUt7m_reverse_26492 <= 0.0
0.0 <= GLYt2r <= 100.0
0.0 <= GLYt2r_reverse_d768d <= 100.0
0.0 <= GND <= 100.0
0.0 <= GND_reverse_eec5c <= 0.0
0.0 <= H2Ot <= 100.0
0.0 <= H2Ot_reverse_aa560 <= 100.0
0.0 <= H2Oter <= 100.0
0.0 <= H2Oter_reverse_c874d <= 100.0
0.0 <= H2Otm <= 100.0
0.0 <= H2Otm_reverse_d5f1a <= 100.0
0.0 <= HCITSm <= 100.0
0.0 <= HCITSm_reverse_145c0 <= 0.0
0.0 <= HCO3E <= 100.0
0.0 <= HCO3E_reverse_97ea5 <= 100.0
0.0 <= HEX1 <= 100.0
0.0 <= HEX1_reverse_25efa <= 0.0
0.0 <= HISt2r <= 100.0
0.0 <= HISt2r_reverse_b00a3 <= 100.0
0.0 <= ICDHxm <= 100.0
0.0 <= ICDHxm_reverse_14a38 <= 0.0
0.0 <= ICDHy <= 100.0
0.0 <= ICDHy_reverse_c12f8 <= 0.0
0.0 <= ICDHym <= 100.0
0.0 <= ICDHym_reverse_38fb7 <= 0.0
0.0 <= ILEt2r <= 100.0
0.0 <= ILEt2r_reverse_9891b <= 100.0
0.0 <= ILEtmi <= 100.0
0.0 <= ILEtmi_reverse_5c3e7 <= 0.0
0.0 <= LEUt2r <= 100.0
0.0 <= LEUt2r_reverse_4a477 <= 100.0
0.0 <= LYSt2r <= 100.0
0.0 <= LYSt2r_reverse_7e0ef <= 100.0
0.0 <= MALt2r <= 100.0
0.0 <= MALt2r_reverse_bfe2f <= 100.0
0.0 <= MALtm <= 100.0
0.0 <= MALtm_reverse_f798a <= 100.0
0.0 <= MAN6PI <= 100.0
0.0 <= MAN6PI_reverse_d96f0 <= 100.0
0.0 <= MANNANter <= 100.0
0.0 <= MANNANter_reverse_6b4bd <= 100.0
0.0 <= MDH <= 100.0
0.0 <= MDH_reverse_ee52c <= 100.0
0.0 <= MDHm <= 100.0
0.0 <= MDHm_reverse_80605 <= 100.0
0.0 <= ME1m <= 100.0
0.0 <= ME1m_reverse_6868b <= 0.0
0.0 <= ME2m <= 100.0
0.0 <= ME2m_reverse_f04a4 <= 0.0
0.0 <= METt2r <= 100.0
0.0 <= METt2r_reverse_22908 <= 100.0
0.0 <= NADH2-u6cm <= 100.0
0.0 <= NADH2-u6cm_reverse_07733 <= 0.0
0.0 <= NADH2-u6m <= 100.0
0.0 <= NADH2-u6m_reverse_18109 <= 0.0
0.0 <= NADK <= 100.0
0.0 <= NADK_reverse_bba52 <= 0.0
0.0 <= NADKm <= 100.0
0.0 <= NADKm_reverse_f0454 <= 0.0
0.0 <= NADPPPS <= 100.0
0.0 <= NADPPPS_reverse_c9c40 <= 0.0
0.0 <= NADPPPSm <= 100.0
0.0 <= NADPPPSm_reverse_64d27 <= 0.0
0.0 <= NDP1 <= 100.0
0.0 <= NDP1_reverse_58712 <= 0.0
0.0 <= NH4t <= 100.0
0.0 <= NH4t_reverse_551ee <= 100.0
0.0 <= NH4tm <= 100.0
0.0 <= NH4tm_reverse_e28df <= 100.0
0.0 <= O2t <= 100.0
0.0 <= O2t_reverse_4d957 <= 100.0
0.0 <= O2ter <= 100.0
0.0 <= O2ter_reverse_ea1ad <= 100.0
0.0 <= O2tm <= 100.0
0.0 <= O2tm_reverse_fe4d6 <= 100.0
0.0 <= OAAt <= 100.0
0.0 <= OAAt_reverse_d0cc5 <= 100.0
0.0 <= OAAt2m <= 100.0
0.0 <= OAAt2m_reverse_42ec8 <= 100.0
0.0 <= PC <= 100.0
0.0 <= PC_reverse_88dba <= 0.0
0.0 <= PDHcm <= 100.0
0.0 <= PDHcm_reverse_40ed3 <= 0.0
0.0 <= PDHm <= 100.0
0.0 <= PDHm_reverse_31d17 <= 0.0
0.0 <= PEtm_SC <= 100.0
0.0 <= PEtm_SC_reverse_a93d3 <= 100.0
0.0 <= PFK <= 100.0
0.0 <= PFK_reverse_d24a6 <= 0.0
0.0 <= PFK_3 <= 100.0
0.0 <= PFK_3_reverse_037e0 <= 0.0
0.0 <= PGI <= 100.0
0.0 <= PGI_reverse_27efc <= 100.0
0.0 <= PGK <= 100.0
0.0 <= PGK_reverse_02696 <= 100.0
0.0 <= PGL <= 100.0
0.0 <= PGL_reverse_2bb6b <= 0.0
0.0 <= PGM <= 100.0
0.0 <= PGM_reverse_fc9af <= 100.0
0.0 <= PGMT <= 100.0
0.0 <= PGMT_reverse_5bcdd <= 100.0
0.0 <= PHEt2r <= 100.0
0.0 <= PHEt2r_reverse_19e89 <= 100.0
0.0 <= PIt2m <= 100.0
0.0 <= PIt2m_reverse_cc401 <= 100.0
0.0 <= PIt2r <= 100.0
0.0 <= PIt2r_reverse_1cd61 <= 100.0
0.0 <= PPA <= 100.0
0.0 <= PPA_reverse_c5293 <= 0.0
0.0 <= PPAm <= 100.0
0.0 <= PPAm_reverse_6c806 <= 0.0
0.0 <= PPCK <= 100.0
0.0 <= PPCK_reverse_2557d <= 0.0
0.0 <= PROt2r <= 100.0
0.0 <= PROt2r_reverse_f59df <= 100.0
0.0 <= PStm_SC <= 100.0
0.0 <= PStm_SC_reverse_a026e <= 100.0
0.0 <= PYK <= 100.0
0.0 <= PYK_reverse_bc8ff <= 0.0
0.0 <= PYRDC <= 100.0
0.0 <= PYRDC_reverse_51259 <= 0.0
0.0 <= PYRt <= 100.0
0.0 <= PYRt_reverse_6fee7 <= 0.0
0.0 <= PYRt2m <= 100.0
0.0 <= PYRt2m_reverse_fb762 <= 100.0
0.0 <= RPE <= 100.0
0.0 <= RPE_reverse_a1b04 <= 100.0
0.0 <= RPI <= 100.0
0.0 <= RPI_reverse_853a1 <= 100.0
0.0 <= SERt2r <= 100.0
0.0 <= SERt2r_reverse_1b436 <= 100.0
0.0 <= SO4ti <= 100.0
0.0 <= SO4ti_reverse_4a93d <= 0.0
0.0 <= SUCCt2r <= 100.0
0.0 <= SUCCt2r_reverse_e3abe <= 100.0
0.0 <= SUCCtm <= 100.0
0.0 <= SUCCtm_reverse_5f36b <= 0.0
0.0 <= SUCD1m <= 100.0
0.0 <= SUCD1m_reverse_60ba1 <= 100.0
0.0 <= SUCD2_u6m <= 100.0
0.0 <= SUCD2_u6m_reverse_b6f79 <= 100.0
0.0 <= SUCD3-u6m <= 100.0
0.0 <= SUCD3-u6m_reverse_0e152 <= 100.0
0.0 <= SUCFUMtm <= 100.0
0.0 <= SUCFUMtm_reverse_99c3a <= 0.0
0.0 <= SUCOASm <= 100.0
0.0 <= SUCOASm_reverse_87069 <= 100.0
0.0 <= TALA <= 100.0
0.0 <= TALA_reverse_adfda <= 100.0
0.0 <= THRt2m <= 100.0
0.0 <= THRt2m_reverse_fd5a8 <= 100.0
0.0 <= THRt2r <= 100.0
0.0 <= THRt2r_reverse_721ff <= 100.0
0.0 <= TKT1 <= 100.0
0.0 <= TKT1_reverse_a1021 <= 100.0
0.0 <= TKT2 <= 100.0
0.0 <= TKT2_reverse_7ebc7 <= 100.0
0.0 <= TPI <= 100.0
0.0 <= TPI_reverse_c2c3b <= 100.0
0.0 <= TRPt2r <= 100.0
0.0 <= TRPt2r_reverse_159c4 <= 100.0
0.0 <= TYRt2r <= 100.0
0.0 <= TYRt2r_reverse_6f73b <= 100.0
0.0 <= VALt2r <= 100.0
0.0 <= VALt2r_reverse_e3d24 <= 100.0
0.0 <= ZYMSTt <= 100.0
0.0 <= ZYMSTt_reverse_ab42f <= 100.0
0.0 <= LMPD_1_13BDglcn_c <= 100.0
0.0 <= LMPD_1_13BDglcn_c_reverse_e2bdf <= 0.0
0.0 <= LMPD_2_amp_c <= 100.0
0.0 <= LMPD_2_amp_c_reverse_bde95 <= 0.0
0.0 <= LMPD_3_amp_c <= 100.0
0.0 <= LMPD_3_amp_c_reverse_96e30 <= 0.0
0.0 <= LMPD_4_amp_c <= 100.0
0.0 <= LMPD_4_amp_c_reverse_7925f <= 0.0
0.0 <= LMPD_5_amp_c <= 100.0
0.0 <= LMPD_5_amp_c_reverse_8e974 <= 0.0
0.0 <= LMPD_6_arg-L_c <= 100.0
0.0 <= LMPD_6_arg-L_c_reverse_33f31 <= 0.0
0.0 <= LMPD_7_asn-L_c <= 100.0
0.0 <= LMPD_7_asn-L_c_reverse_b41c7 <= 0.0
0.0 <= LMPD_8_cmp_c <= 100.0
0.0 <= LMPD_8_cmp_c_reverse_23cbb <= 0.0
0.0 <= LMPD_9_cmp_c <= 100.0
0.0 <= LMPD_9_cmp_c_reverse_85279 <= 0.0
0.0 <= LMPD_10_cmp_c <= 100.0
0.0 <= LMPD_10_cmp_c_reverse_3acfa <= 0.0
0.0 <= LMPD_11_cmp_c <= 100.0
0.0 <= LMPD_11_cmp_c_reverse_98fa8 <= 0.0
0.0 <= LMPD_12_cys-L_c <= 100.0
0.0 <= LMPD_12_cys-L_c_reverse_6dc0a <= 0.0
0.0 <= LMPD_13_damp_c <= 100.0
0.0 <= LMPD_13_damp_c_reverse_e961c <= 0.0
0.0 <= LMPD_14_damp_c <= 100.0
0.0 <= LMPD_14_damp_c_reverse_d5737 <= 0.0
0.0 <= LMPD_15_damp_c <= 100.0
0.0 <= LMPD_15_damp_c_reverse_33cc4 <= 0.0
0.0 <= LMPD_16_damp_c <= 100.0
0.0 <= LMPD_16_damp_c_reverse_be99f <= 0.0
0.0 <= LMPD_17_dcmp_c <= 100.0
0.0 <= LMPD_17_dcmp_c_reverse_5c6b0 <= 0.0
0.0 <= LMPD_18_dcmp_c <= 100.0
0.0 <= LMPD_18_dcmp_c_reverse_15021 <= 0.0
0.0 <= LMPD_19_dgmp_c <= 100.0
0.0 <= LMPD_19_dgmp_c_reverse_7f2b1 <= 0.0
0.0 <= LMPD_20_dgmp_c <= 100.0
0.0 <= LMPD_20_dgmp_c_reverse_e4f3a <= 0.0
0.0 <= LMPD_21_dgmp_c <= 100.0
0.0 <= LMPD_21_dgmp_c_reverse_65045 <= 0.0
0.0 <= LMPD_22_dgmp_c <= 100.0
0.0 <= LMPD_22_dgmp_c_reverse_ebc77 <= 0.0
0.0 <= LMPD_23_dtmp_c <= 100.0
0.0 <= LMPD_23_dtmp_c_reverse_ed12d <= 0.0
0.0 <= LMPD_24_dtmp_c <= 100.0
0.0 <= LMPD_24_dtmp_c_reverse_5d241 <= 0.0
0.0 <= LMPD_25_ergst_c <= 100.0
0.0 <= LMPD_25_ergst_c_reverse_96537 <= 0.0
0.0 <= LMPD_26_ergst_c <= 100.0
0.0 <= LMPD_26_ergst_c_reverse_a59e0 <= 0.0
0.0 <= LMPD_27_gln-L_c <= 100.0
0.0 <= LMPD_27_gln-L_c_reverse_77507 <= 0.0
0.0 <= LMPD_28_gly_c <= 100.0
0.0 <= LMPD_28_gly_c_reverse_77428 <= 0.0
0.0 <= LMPD_29_glycogen_c <= 100.0
0.0 <= LMPD_29_glycogen_c_reverse_68cb3 <= 0.0
0.0 <= LMPD_30_gmp_c <= 100.0
0.0 <= LMPD_30_gmp_c_reverse_1f79a <= 0.0
0.0 <= LMPD_31_gmp_c <= 100.0
0.0 <= LMPD_31_gmp_c_reverse_d0c02 <= 0.0
0.0 <= LMPD_32_his-L_c <= 100.0
0.0 <= LMPD_32_his-L_c_reverse_c5967 <= 0.0
0.0 <= LMPD_33_his-L_c <= 100.0
0.0 <= LMPD_33_his-L_c_reverse_0eb8b <= 0.0
0.0 <= LMPD_34_his-L_c <= 100.0
0.0 <= LMPD_34_his-L_c_reverse_e16a9 <= 0.0
0.0 <= LMPD_35_his-L_c <= 100.0
0.0 <= LMPD_35_his-L_c_reverse_5e623 <= 0.0
0.0 <= LMPD_36_ile-L_c <= 100.0
0.0 <= LMPD_36_ile-L_c_reverse_d2bf7 <= 0.0
0.0 <= LMPD_37_ile-L_c <= 100.0
0.0 <= LMPD_37_ile-L_c_reverse_1ea52 <= 0.0
0.0 <= LMPD_38_leu-L_c <= 100.0
0.0 <= LMPD_38_leu-L_c_reverse_0e66c <= 0.0
0.0 <= LMPD_39_leu-L_c <= 100.0
0.0 <= LMPD_39_leu-L_c_reverse_73152 <= 0.0
0.0 <= LMPD_40_lys-L_c <= 100.0
0.0 <= LMPD_40_lys-L_c_reverse_4d3a6 <= 0.0
0.0 <= LMPD_41_lys-L_c <= 100.0
0.0 <= LMPD_41_lys-L_c_reverse_c021a <= 0.0
0.0 <= LMPD_42_mannan_c <= 100.0
0.0 <= LMPD_42_mannan_c_reverse_9f626 <= 0.0
0.0 <= LMPD_43_mannan_c <= 100.0
0.0 <= LMPD_43_mannan_c_reverse_7dd82 <= 0.0
0.0 <= LMPD_44_met-L_c <= 100.0
0.0 <= LMPD_44_met-L_c_reverse_1646f <= 0.0
0.0 <= LMPD_45_met-L_c <= 100.0
0.0 <= LMPD_45_met-L_c_reverse_912f1 <= 0.0
0.0 <= LMPD_47_pa_SC_c <= 100.0
0.0 <= LMPD_47_pa_SC_c_reverse_4049b <= 0.0
0.0 <= LMPD_48_pc_SC_c <= 100.0
0.0 <= LMPD_48_pc_SC_c_reverse_be709 <= 0.0
0.0 <= LMPD_49_pc_SC_c <= 100.0
0.0 <= LMPD_49_pc_SC_c_reverse_55f95 <= 0.0
0.0 <= LMPD_51_pe_SC_c <= 100.0
0.0 <= LMPD_51_pe_SC_c_reverse_09777 <= 0.0
0.0 <= LMPD_52_phe-L_c <= 100.0
0.0 <= LMPD_52_phe-L_c_reverse_909a1 <= 0.0
0.0 <= LMPD_53_pro-L_c <= 100.0
0.0 <= LMPD_53_pro-L_c_reverse_830d7 <= 0.0
0.0 <= LMPD_54_pro-L_c <= 100.0
0.0 <= LMPD_54_pro-L_c_reverse_a6168 <= 0.0
0.0 <= LMPD_56_ps_SC_c <= 100.0
0.0 <= LMPD_56_ps_SC_c_reverse_7e9d7 <= 0.0
0.0 <= LMPD_57_ptd1ino_SC_c <= 100.0
0.0 <= LMPD_57_ptd1ino_SC_c_reverse_65089 <= 0.0
0.0 <= LMPD_58_ptd1ino_SC_c <= 100.0
0.0 <= LMPD_58_ptd1ino_SC_c_reverse_d40c6 <= 0.0
0.0 <= LMPD_59_ribflv_c <= 100.0
0.0 <= LMPD_59_ribflv_c_reverse_01d9b <= 0.0
0.0 <= LMPD_60_ribflv_c <= 100.0
0.0 <= LMPD_60_ribflv_c_reverse_b7a89 <= 0.0
0.0 <= LMPD_61_ser-L_c <= 100.0
0.0 <= LMPD_61_ser-L_c_reverse_40666 <= 0.0
0.0 <= LMPD_62_thr-L_c <= 100.0
0.0 <= LMPD_62_thr-L_c_reverse_e4d09 <= 0.0
0.0 <= LMPD_63_thr-L_c <= 100.0
0.0 <= LMPD_63_thr-L_c_reverse_52245 <= 0.0
0.0 <= LMPD_64_tre_c <= 100.0
0.0 <= LMPD_64_tre_c_reverse_64f54 <= 0.0
0.0 <= LMPD_65_triglyc_SC_c <= 100.0
0.0 <= LMPD_65_triglyc_SC_c_reverse_7f8d4 <= 0.0
0.0 <= LMPD_66_triglyc_SC_c <= 100.0
0.0 <= LMPD_66_triglyc_SC_c_reverse_42d98 <= 0.0
0.0 <= LMPD_67_trp-L_c <= 100.0
0.0 <= LMPD_67_trp-L_c_reverse_17a73 <= 0.0
0.0 <= LMPD_68_tyr-L_c <= 100.0
0.0 <= LMPD_68_tyr-L_c_reverse_6983e <= 0.0
0.0 <= LMPD_69_tyr-L_c <= 100.0
0.0 <= LMPD_69_tyr-L_c_reverse_03f73 <= 0.0
0.0 <= LMPD_70_ump_c <= 100.0
0.0 <= LMPD_70_ump_c_reverse_3dcab <= 0.0
0.0 <= LMPD_71_ump_c <= 100.0
0.0 <= LMPD_71_ump_c_reverse_dc9d0 <= 0.0
0.0 <= LMPD_72_val-L_c <= 100.0
0.0 <= LMPD_72_val-L_c_reverse_768ff <= 0.0
0.0 <= LMPD_73_zymst_c <= 100.0
0.0 <= LMPD_73_zymst_c_reverse_e6142 <= 0.0
0.0 <= LMPD_74_zymst_c <= 100.0
0.0 <= LMPD_74_zymst_c_reverse_84617 <= 0.0
-16.11809565095832 <= LC_val-L_c <= -4.605170185988091
-16.11809565095832 <= LC_13BDglcn_c <= -4.605170185988091
-16.11809565095832 <= LC_pro-L_c <= -4.605170185988091
-16.11809565095832 <= LC_damp_c <= -4.605170185988091
-16.11809565095832 <= LC_trp-L_c <= -4.605170185988091
-16.11809565095832 <= LC_gly_c <= -4.605170185988091
-16.11809565095832 <= LC_tyr-L_c <= -4.605170185988091
-16.11809565095832 <= LC_dcmp_c <= -4.605170185988091
-16.11809565095832 <= LC_ump_c <= -4.605170185988091
-16.11809565095832 <= LC_dgmp_c <= -4.605170185988091
-16.11809565095832 <= LC_mannan_c <= -4.605170185988091
-16.11809565095832 <= LC_ala-L_c <= -4.605170185988091
-16.11809565095832 <= LC_amp_c <= -4.605170185988091
-16.11809565095832 <= LC_arg-L_c <= -4.605170185988091
-16.11809565095832 <= LC_met-L_c <= -4.605170185988091
-16.11809565095832 <= LC_ribflv_c <= -4.605170185988091
-16.11809565095832 <= LC_dtmp_c <= -4.605170185988091
-16.11809565095832 <= LC_ile-L_c <= -4.605170185988091
-16.11809565095832 <= LC_ergst_c <= -4.605170185988091
-16.11809565095832 <= LC_gmp_c <= -4.605170185988091
-16.11809565095832 <= LC_gln-L_c <= -4.605170185988091
-16.11809565095832 <= LC_leu-L_c <= -4.605170185988091
-16.11809565095832 <= LC_glu-L_c <= -4.605170185988091
-16.11809565095832 <= LC_lys-L_c <= -4.605170185988091
-16.11809565095832 <= LC_zymst_c <= -4.605170185988091
-16.11809565095832 <= LC_ser-L_c <= -4.605170185988091
-16.11809565095832 <= LC_asn-L_c <= -4.605170185988091
-16.11809565095832 <= LC_asp-L_c <= -4.605170185988091
-16.11809565095832 <= LC_atp_c <= -4.605170185988091
-16.11809565095832 <= LC_glycogen_c <= -4.605170185988091
-16.11809565095832 <= LC_phe-L_c <= -4.605170185988091
-16.11809565095832 <= LC_h_c <= -16.11809565095832
-16.11809565095832 <= LC_so4_c <= -4.605170185988091
-16.11809565095832 <= LC_adp_c <= -4.605170185988091
-16.11809565095832 <= LC_thr-L_c <= -4.605170185988091
-16.11809565095832 <= LC_cmp_c <= -4.605170185988091
-16.11809565095832 <= LC_tre_c <= -4.605170185988091
-16.11809565095832 <= LC_his-L_c <= -4.605170185988091
0 <= LC_h2o_c <= 0
-16.11809565095832 <= LC_cys-L_c <= -4.605170185988091
-16.11809565095832 <= LC_pi_c <= -4.605170185988091
-16.11809565095832 <= LC_4abut_e <= -4.605170185988091
-16.11809565095832 <= LC_ac_e <= -4.605170185988091
-16.11809565095832 <= LC_acald_e <= -4.605170185988091
-16.11809565095832 <= LC_akg_e <= -4.605170185988091
-16.11809565095832 <= LC_ala-L_e <= -4.605170185988091
-16.11809565095832 <= LC_arg-L_e <= -4.605170185988091
-16.11809565095832 <= LC_asn-L_e <= -4.605170185988091
-16.11809565095832 <= LC_asp-L_e <= -4.605170185988091
-16.11809565095832 <= LC_cit_e <= -4.605170185988091
-16.11809565095832 <= LC_co2_e <= -4.605170185988091
-16.11809565095832 <= LC_cys-L_e <= -4.605170185988091
-16.11809565095832 <= LC_ergst_e <= -4.605170185988091
-16.11809565095832 <= LC_etoh_e <= -4.605170185988091
-16.11809565095832 <= LC_for_e <= -4.605170185988091
-16.11809565095832 <= LC_fum_e <= -4.605170185988091
-16.11809565095832 <= LC_glc-D_e <= -4.605170185988091
-16.11809565095832 <= LC_gln-L_e <= -4.605170185988091
-16.11809565095832 <= LC_glu-L_e <= -4.605170185988091
-16.11809565095832 <= LC_gly_e <= -4.605170185988091
-11.512925464970229 <= LC_h_e <= -11.512925464970229
0 <= LC_h2o_e <= 0
-16.11809565095832 <= LC_his-L_e <= -4.605170185988091
-16.11809565095832 <= LC_ile-L_e <= -4.605170185988091
-16.11809565095832 <= LC_leu-L_e <= -4.605170185988091
-16.11809565095832 <= LC_lys-L_e <= -4.605170185988091
-16.11809565095832 <= LC_mal-L_e <= -4.605170185988091
-16.11809565095832 <= LC_met-L_e <= -4.605170185988091
-16.11809565095832 <= LC_nh4_e <= -4.605170185988091
-16.11809565095832 <= LC_o2_e <= -4.605170185988091
-16.11809565095832 <= LC_oaa_e <= -4.605170185988091
-16.11809565095832 <= LC_phe-L_e <= -4.605170185988091
-16.11809565095832 <= LC_pi_e <= -4.605170185988091
-16.11809565095832 <= LC_pro-L_e <= -4.605170185988091
-16.11809565095832 <= LC_pyr_e <= -4.605170185988091
-16.11809565095832 <= LC_ser-L_e <= -4.605170185988091
-16.11809565095832 <= LC_so4_e <= -4.605170185988091
-16.11809565095832 <= LC_succ_e <= -4.605170185988091
-16.11809565095832 <= LC_thr-L_e <= -4.605170185988091
-16.11809565095832 <= LC_trp-L_e <= -4.605170185988091
-16.11809565095832 <= LC_tyr-L_e <= -4.605170185988091
-16.11809565095832 <= LC_val-L_e <= -4.605170185988091
-16.11809565095832 <= LC_zymst_e <= -4.605170185988091
-16.11809565095832 <= LC_6pgl_r <= -4.605170185988091
-16.11809565095832 <= LC_6pgl_c <= -4.605170185988091
-16.11809565095832 <= LC_4abut_c <= -4.605170185988091
-16.11809565095832 <= LC_acald_c <= -4.605170185988091
-16.11809565095832 <= LC_acald_m <= -4.605170185988091
-17.269388197455342 <= LC_h_m <= -17.269388197455342
-16.11809565095832 <= LC_coa_m <= -4.605170185988091
-16.11809565095832 <= LC_accoa_m <= -4.605170185988091
-16.11809565095832 <= LC_ac_m <= -4.605170185988091
0 <= LC_h2o_m <= 0
-16.11809565095832 <= LC_cit_c <= -4.605170185988091
-16.11809565095832 <= LC_icit_c <= -4.605170185988091
-16.11809565095832 <= LC_cit_m <= -4.605170185988091
-16.11809565095832 <= LC_icit_m <= -4.605170185988091
-16.11809565095832 <= LC_ac_c <= -4.605170185988091
-16.11809565095832 <= LC_ppi_c <= -4.605170185988091
-16.11809565095832 <= LC_accoa_c <= -4.605170185988091
-16.11809565095832 <= LC_coa_c <= -4.605170185988091
-16.11809565095832 <= LC_atp_m <= -4.605170185988091
-16.11809565095832 <= LC_ppi_m <= -4.605170185988091
-16.11809565095832 <= LC_amp_m <= -4.605170185988091
-16.11809565095832 <= LC_adp_m <= -4.605170185988091
-16.11809565095832 <= LC_co2_m <= -4.605170185988091
-16.11809565095832 <= LC_lpam_m <= -4.605170185988091
-16.11809565095832 <= LC_sdhlam_m <= -4.605170185988091
-16.11809565095832 <= LC_akg_m <= -4.605170185988091
-16.11809565095832 <= LC_succoa_m <= -4.605170185988091
-16.11809565095832 <= LC_dhlam_m <= -4.605170185988091
-16.11809565095832 <= LC_akg_c <= -4.605170185988091
-16.11809565095832 <= LC_mal-L_c <= -4.605170185988091
-16.11809565095832 <= LC_pyr_c <= -4.605170185988091
-16.11809565095832 <= LC_ala-L_m <= -4.605170185988091
-16.11809565095832 <= LC_glu-L_m <= -4.605170185988091
-16.11809565095832 <= LC_pyr_m <= -4.605170185988091
-16.11809565095832 <= LC_nadh_c <= -4.605170185988091
-16.11809565095832 <= LC_nad_c <= -4.605170185988091
-16.11809565095832 <= LC_etoh_c <= -4.605170185988091
-16.11809565095832 <= LC_nadh_m <= -4.605170185988091
-16.11809565095832 <= LC_nad_m <= -4.605170185988091
-16.11809565095832 <= LC_etoh_m <= -4.605170185988091
-16.11809565095832 <= LC_nadph_c <= -4.605170185988091
-16.11809565095832 <= LC_nadp_c <= -4.605170185988091
-16.11809565095832 <= LC_nadp_m <= -4.605170185988091
-16.11809565095832 <= LC_nadph_m <= -4.605170185988091
-16.11809565095832 <= LC_asp-L_m <= -4.605170185988091
-16.11809565095832 <= LC_oaa_c <= -4.605170185988091
-16.11809565095832 <= LC_oaa_m <= -4.605170185988091
-16.11809565095832 <= LC_pi_m <= -4.605170185988091
-16.11809565095832 <= LC_mal-L_m <= -4.605170185988091
-16.11809565095832 <= LC_co2_c <= -4.605170185988091
-16.11809565095832 <= LC_ficytc_m <= -4.605170185988091
-16.11809565095832 <= LC_focytc_m <= -4.605170185988091
-16.11809565095832 <= LC_o2_m <= -4.605170185988091
-16.11809565095832 <= LC_q6h2_m <= -4.605170185988091
-16.11809565095832 <= LC_q6_m <= -4.605170185988091
-16.11809565095832 <= LC_lac-D_c <= -4.605170185988091
-16.11809565095832 <= LC_lac-D_m <= -4.605170185988091
-16.11809565095832 <= LC_db4p_c <= -4.605170185988091
-16.11809565095832 <= LC_for_c <= -4.605170185988091
-16.11809565095832 <= LC_ru5p-D_c <= -4.605170185988091
-16.11809565095832 <= LC_pep_c <= -4.605170185988091
-16.11809565095832 <= LC_2pg_c <= -4.605170185988091
-16.11809565095832 <= LC_ergst_r <= -4.605170185988091
-16.11809565095832 <= LC_dhap_c <= -4.605170185988091
-16.11809565095832 <= LC_fdp_c <= -4.605170185988091
-16.11809565095832 <= LC_g3p_c <= -4.605170185988091
-16.11809565095832 <= LC_s17bp_c <= -4.605170185988091
-16.11809565095832 <= LC_e4p_c <= -4.605170185988091
-16.11809565095832 <= LC_f6p_c <= -4.605170185988091
-16.11809565095832 <= LC_fum_c <= -4.605170185988091
-16.11809565095832 <= LC_fad_m <= -4.605170185988091
-16.11809565095832 <= LC_succ_c <= -4.605170185988091
-16.11809565095832 <= LC_fadh2_m <= -4.605170185988091
-16.11809565095832 <= LC_succ_m <= -4.605170185988091
-16.11809565095832 <= LC_fum_m <= -4.605170185988091
-16.11809565095832 <= LC_g6p_c <= -4.605170185988091
-16.11809565095832 <= LC_g6p_r <= -4.605170185988091
-16.11809565095832 <= LC_h_r <= -16.11809565095832
-16.11809565095832 <= LC_nadp_r <= -4.605170185988091
-16.11809565095832 <= LC_nadph_r <= -4.605170185988091
-16.11809565095832 <= LC_g6p-B_c <= -4.605170185988091
-16.11809565095832 <= LC_13dpg_c <= -4.605170185988091
-16.11809565095832 <= LC_glc-D_c <= -4.605170185988091
-16.11809565095832 <= LC_nh4_c <= -4.605170185988091
-16.11809565095832 <= LC_6pgc_c <= -4.605170185988091
0 <= LC_h2o_r <= 0
-16.11809565095832 <= LC_hcit_m <= -4.605170185988091
-16.11809565095832 <= LC_hco3_c <= -4.605170185988091
-16.11809565095832 <= LC_ile-L_m <= -4.605170185988091
-16.11809565095832 <= LC_man6p_c <= -4.605170185988091
-16.11809565095832 <= LC_mannan_r <= -4.605170185988091
-16.11809565095832 <= LC_nh4_m <= -4.605170185988091
-16.11809565095832 <= LC_o2_c <= -4.605170185988091
-16.11809565095832 <= LC_o2_r <= -4.605170185988091
-16.11809565095832 <= LC_s7p_c <= -4.605170185988091
-16.11809565095832 <= LC_3pg_c <= -4.605170185988091
-16.11809565095832 <= LC_g1p_c <= -4.605170185988091
-16.11809565095832 <= LC_xu5p-D_c <= -4.605170185988091
-16.11809565095832 <= LC_r5p_c <= -4.605170185988091
-16.11809565095832 <= LC_thr-L_m <= -4.605170185988091
0 <= FU_biomass_SC5_notrace <= 1
0 <= BU_biomass_SC5_notrace <= 1
0 <= FU_EX_4abut_e <= 1
0 <= BU_EX_4abut_e <= 1
0 <= FU_EX_ac_e <= 1
0 <= BU_EX_ac_e <= 1
0 <= FU_EX_acald_e <= 1
0 <= BU_EX_acald_e <= 1
0 <= FU_EX_akg_e <= 1
0 <= BU_EX_akg_e <= 1
0 <= FU_EX_ala-L_e <= 1
0 <= BU_EX_ala-L_e <= 1
0 <= FU_EX_arg-L_e <= 1
0 <= BU_EX_arg-L_e <= 1
0 <= FU_EX_asn-L_e <= 1
0 <= BU_EX_asn-L_e <= 1
0 <= FU_EX_asp-L_e <= 1
0 <= BU_EX_asp-L_e <= 1
0 <= FU_EX_cit_e <= 1
0 <= BU_EX_cit_e <= 1
0 <= FU_EX_co2_e <= 1
0 <= BU_EX_co2_e <= 1
0 <= FU_EX_cys-L_e <= 1
0 <= BU_EX_cys-L_e <= 1
0 <= FU_EX_ergst_e <= 1
0 <= BU_EX_ergst_e <= 1
0 <= FU_EX_etoh_e <= 1
0 <= BU_EX_etoh_e <= 1
0 <= FU_EX_for_e <= 1
0 <= BU_EX_for_e <= 1
0 <= FU_EX_fum_e <= 1
0 <= BU_EX_fum_e <= 1
0 <= FU_EX_glc_e <= 1
0 <= BU_EX_glc_e <= 1
0 <= FU_EX_gln-L_e <= 1
0 <= BU_EX_gln-L_e <= 1
0 <= FU_EX_glu-L_e <= 1
0 <= BU_EX_glu-L_e <= 1
0 <= FU_EX_gly_e <= 1
0 <= BU_EX_gly_e <= 1
0 <= FU_EX_h_e <= 1
0 <= BU_EX_h_e <= 1
0 <= FU_EX_h2o_e <= 1
0 <= BU_EX_h2o_e <= 1
0 <= FU_EX_his-L_e <= 1
0 <= BU_EX_his-L_e <= 1
0 <= FU_EX_ile-L_e <= 1
0 <= BU_EX_ile-L_e <= 1
0 <= FU_EX_leu-L_e <= 1
0 <= BU_EX_leu-L_e <= 1
0 <= FU_EX_lys-L_e <= 1
0 <= BU_EX_lys-L_e <= 1
0 <= FU_EX_mal-L_e <= 1
0 <= BU_EX_mal-L_e <= 1
0 <= FU_EX_met-L_e <= 1
0 <= BU_EX_met-L_e <= 1
0 <= FU_EX_nh4_e <= 1
0 <= BU_EX_nh4_e <= 1
0 <= FU_EX_o2_e <= 1
0 <= BU_EX_o2_e <= 1
0 <= FU_EX_oaa_e <= 1
0 <= BU_EX_oaa_e <= 1
0 <= FU_EX_phe-L_e <= 1
0 <= BU_EX_phe-L_e <= 1
0 <= FU_EX_pi_e <= 1
0 <= BU_EX_pi_e <= 1
0 <= FU_EX_pro-L_e <= 1
0 <= BU_EX_pro-L_e <= 1
0 <= FU_EX_pyr_e <= 1
0 <= BU_EX_pyr_e <= 1
0 <= FU_EX_ser-L_e <= 1
0 <= BU_EX_ser-L_e <= 1
0 <= FU_EX_so4_e <= 1
0 <= BU_EX_so4_e <= 1
0 <= FU_EX_succ_e <= 1
0 <= BU_EX_succ_e <= 1
0 <= FU_EX_thr-L_e <= 1
0 <= BU_EX_thr-L_e <= 1
0 <= FU_EX_trp-L_e <= 1
0 <= BU_EX_trp-L_e <= 1
0 <= FU_EX_tyr-L_e <= 1
0 <= BU_EX_tyr-L_e <= 1
0 <= FU_EX_val-L_e <= 1
0 <= BU_EX_val-L_e <= 1
0 <= FU_EX_zymst_e <= 1
0 <= BU_EX_zymst_e <= 1
-1000 <= DG_6PGLter <= 1000
-3.15412860233 <= DGo_6PGLter <= 0.845871397668
0 <= FU_6PGLter <= 1
0 <= BU_6PGLter <= 1
-1000 <= DG_ABUTt2r <= 1000
-7.31065797989 <= DGo_ABUTt2r <= -3.31065797989
0 <= FU_ABUTt2r <= 1
0 <= BU_ABUTt2r <= 1
-1000 <= DG_ACALDt <= 1000
-1.22593477799 <= DGo_ACALDt <= 2.77406522201
0 <= FU_ACALDt <= 1
0 <= BU_ACALDt <= 1
-1000 <= DG_ACALDtm <= 1000
-2.0 <= DGo_ACALDtm <= 2.0
0 <= FU_ACALDtm <= 1
0 <= BU_ACALDtm <= 1
-1000 <= DG_ACOAHim <= 1000
-7.67768989177 <= DGo_ACOAHim <= -5.06460129778
0 <= FU_ACOAHim <= 1
0 <= BU_ACOAHim <= 1
-1000 <= DG_ACONT <= 1000
0.983477834035 <= DGo_ACONT <= 1.85559901828
0 <= FU_ACONT <= 1
0 <= BU_ACONT <= 1
-1000 <= DG_ACONTm <= 1000
0.983790142837 <= DGo_ACONTm <= 1.85591132708
0 <= FU_ACONTm <= 1
0 <= BU_ACONTm <= 1
-1000 <= DG_ACS <= 1000
-4.7093997787 <= DGo_ACS <= -2.05113417873
0 <= FU_ACS <= 1
0 <= BU_ACS <= 1
-1000 <= DG_ACSm <= 1000
-5.05148318662 <= DGo_ACSm <= -2.39321758664
0 <= FU_ACSm <= 1
0 <= BU_ACSm <= 1
-1000 <= DG_ACt2r <= 1000
-4.0771880555 <= DGo_ACt2r <= -0.0771880555013
0 <= FU_ACt2r <= 1
0 <= BU_ACt2r <= 1
-1000 <= DG_ACtm <= 1000
2.15216240891 <= DGo_ACtm <= 6.15216240891
0 <= FU_ACtm <= 1
0 <= BU_ACtm <= 1
-1000 <= DG_ACtr <= 1000
0.0331442891724 <= DGo_ACtr <= 4.03314428917
0 <= FU_ACtr <= 1
0 <= BU_ACtr <= 1
-1000 <= DG_ADK1 <= 1000
-1.56896108821 <= DGo_ADK1 <= 2.43103891179
0 <= FU_ADK1 <= 1
0 <= BU_ADK1 <= 1
-1000 <= DG_ADK1m <= 1000
-1.57045564321 <= DGo_ADK1m <= 2.42954435679
0 <= FU_ADK1m <= 1
0 <= BU_ADK1m <= 1
-1000 <= DG_AKGDam <= 1000
-10.1636954099 <= DGo_AKGDam <= -6.26637973881
0 <= FU_AKGDam <= 1
0 <= BU_AKGDam <= 1
-1000 <= DG_AKGDbm <= 1000
-5.55300094771 <= DGo_AKGDbm <= -1.55300094771
0 <= FU_AKGDbm <= 1
0 <= BU_AKGDbm <= 1
-1000 <= DG_AKGMAL <= 1000
-1.42532982439 <= DGo_AKGMAL <= 2.57467017561
0 <= FU_AKGMAL <= 1
0 <= BU_AKGMAL <= 1
-1000 <= DG_AKGt2r <= 1000
-3.30078969523 <= DGo_AKGt2r <= 0.699210304768
0 <= FU_AKGt2r <= 1
0 <= BU_AKGt2r <= 1
-1000 <= DG_ALATA_L <= 1000
-6.10911998456 <= DGo_ALATA_L <= -2.10911998456
0 <= FU_ALATA_L <= 1
0 <= BU_ALATA_L <= 1
-1000 <= DG_ALATA_Lm <= 1000
-6.79066526314 <= DGo_ALATA_Lm <= -2.79066526314
0 <= FU_ALATA_Lm <= 1
0 <= BU_ALATA_Lm <= 1
-1000 <= DG_ALAt2r <= 1000
-4.75571820616 <= DGo_ALAt2r <= -0.755718206157
0 <= FU_ALAt2r <= 1
0 <= BU_ALAt2r <= 1
-1000 <= DG_ALAtmi <= 1000
-2.0 <= DGo_ALAtmi <= 2.0
0 <= FU_ALAtmi <= 1
0 <= BU_ALAtmi <= 1
-1000 <= DG_ALCD2if <= 1000
3.73002385018 <= DGo_ALCD2if <= 6.25513929309
0 <= FU_ALCD2if <= 1
0 <= BU_ALCD2if <= 1
-1000 <= DG_ALCD2ir <= 1000
-6.25513929309 <= DGo_ALCD2ir <= -3.73002385018
0 <= FU_ALCD2ir <= 1
0 <= BU_ALCD2ir <= 1
-1000 <= DG_ALCD2irm <= 1000
-5.57347120692 <= DGo_ALCD2irm <= -3.04835576401
0 <= FU_ALCD2irm <= 1
0 <= BU_ALCD2irm <= 1
-1000 <= DG_ALCD2x <= 1000
3.73002385018 <= DGo_ALCD2x <= 6.25513929309
0 <= FU_ALCD2x <= 1
0 <= BU_ALCD2x <= 1
-1000 <= DG_ALDD2xm <= 1000
-15.9026468148 <= DGo_ALDD2xm <= -13.1867292903
0 <= FU_ALDD2xm <= 1
0 <= BU_ALDD2xm <= 1
-1000 <= DG_ALDD2y <= 1000
-14.9275256624 <= DGo_ALDD2y <= -12.2116081379
0 <= FU_ALDD2y <= 1
0 <= BU_ALDD2y <= 1
-1000 <= DG_ALDD2ym <= 1000
-16.2896794258 <= DGo_ALDD2ym <= -13.5737619013
0 <= FU_ALDD2ym <= 1
0 <= BU_ALDD2ym <= 1
-1000 <= DG_ARGt2r <= 1000
-4.78476406764 <= DGo_ARGt2r <= -0.78476406764
0 <= FU_ARGt2r <= 1
0 <= BU_ARGt2r <= 1
-1000 <= DG_ASNt2r <= 1000
-4.56220190065 <= DGo_ASNt2r <= -0.562201900654
0 <= FU_ASNt2r <= 1
0 <= BU_ASNt2r <= 1
-1000 <= DG_ASPGLU2m <= 1000
-2.00120494009 <= DGo_ASPGLU2m <= 1.99879505991
0 <= FU_ASPGLU2m <= 1
0 <= BU_ASPGLU2m <= 1
-1000 <= DG_ASPTA <= 1000
-3.1443785382 <= DGo_ASPTA <= 0.855621461795
0 <= FU_ASPTA <= 1
0 <= BU_ASPTA <= 1
-1000 <= DG_ASPTAm <= 1000
-3.14494130332 <= DGo_ASPTAm <= 0.855058696684
0 <= FU_ASPTAm <= 1
0 <= BU_ASPTAm <= 1
-1000 <= DG_ASPt2m <= 1000
-3.36192716967 <= DGo_ASPt2m <= 0.638072830334
0 <= FU_ASPt2m <= 1
0 <= BU_ASPt2m <= 1
-1000 <= DG_ASPt2r <= 1000
-6.77176618096 <= DGo_ASPt2r <= -2.77176618096
0 <= FU_ASPt2r <= 1
0 <= BU_ASPt2r <= 1
-1000 <= DG_ATPH1 <= 1000
-16.2630685005 <= DGo_ATPH1 <= -13.3928517839
0 <= FU_ATPH1 <= 1
0 <= BU_ATPH1 <= 1
-1000 <= DG_ATPM <= 1000
-7.91601479433 <= DGo_ATPM <= -6.48090643606
0 <= FU_ATPM <= 1
0 <= BU_ATPM <= 1
-1000 <= DG_ATPS <= 1000
-7.91601479433 <= DGo_ATPS <= -6.48090643606
0 <= FU_ATPS <= 1
0 <= BU_ATPS <= 1
-1000 <= DG_ATPS3m <= 1000
-7.4528727438 <= DGo_ATPS3m <= -6.01776438552
0 <= FU_ATPS3m <= 1
0 <= BU_ATPS3m <= 1
-1000 <= DG_ATPtm-H <= 1000
-11.0071959583 <= DGo_ATPtm-H <= -7.00719595833
0 <= FU_ATPtm-H <= 1
0 <= BU_ATPtm-H <= 1
-1000 <= DG_CITt2r <= 1000
-0.571942053328 <= DGo_CITt2r <= 3.42805794667
0 <= FU_CITt2r <= 1
0 <= BU_CITt2r <= 1
-1000 <= DG_CITtam <= 1000
2.16271075289 <= DGo_CITtam <= 6.16271075289
0 <= FU_CITtam <= 1
0 <= BU_CITtam <= 1
-1000 <= DG_CITtcm <= 1000
-2.0003123088 <= DGo_CITtcm <= 1.9996876912
0 <= FU_CITtcm <= 1
0 <= BU_CITtcm <= 1
-1000 <= DG_CO2t <= 1000
-2.0 <= DGo_CO2t <= 2.0
0 <= FU_CO2t <= 1
0 <= BU_CO2t <= 1
-1000 <= DG_CO2tm <= 1000
-2.0 <= DGo_CO2tm <= 2.0
0 <= FU_CO2tm <= 1
0 <= BU_CO2tm <= 1
-1000 <= DG_CSm <= 1000
-10.0893809582 <= DGo_CSm <= -7.47216347393
0 <= FU_CSm <= 1
0 <= BU_CSm <= 1
-1000 <= DG_CYOOm <= 1000
-31.2761597068 <= DGo_CYOOm <= -29.0400917293
0 <= FU_CYOOm <= 1
0 <= BU_CYOOm <= 1
-1000 <= DG_CYOR_u6m <= 1000
6.90818794877 <= DGo_CYOR_u6m <= 12.1806628894
0 <= FU_CYOR_u6m <= 1
0 <= BU_CYOR_u6m <= 1
-1000 <= DG_CYSt2r <= 1000
-7.95448045664 <= DGo_CYSt2r <= -3.95448045664
0 <= FU_CYSt2r <= 1
0 <= BU_CYSt2r <= 1
-1000 <= DG_D-LACt2m <= 1000
-1.99987308484 <= DGo_D-LACt2m <= 2.00012691516
0 <= FU_D-LACt2m <= 1
0 <= BU_D-LACt2m <= 1
-1000 <= DG_D-LACtm <= 1000
-2.68154117101 <= DGo_D-LACtm <= 1.31845882899
0 <= FU_D-LACtm <= 1
0 <= BU_D-LACtm <= 1
-1000 <= DG_DB4PS <= 1000
-21.7780378364 <= DGo_DB4PS <= -20.5838806608
0 <= FU_DB4PS <= 1
0 <= BU_DB4PS <= 1
-1000 <= DG_ENO <= 1000
-2.3129561054 <= DGo_ENO <= -0.754528898242
0 <= FU_ENO <= 1
0 <= BU_ENO <= 1
-1000 <= DG_ERGSTt <= 1000
6.51471744211 <= DGo_ERGSTt <= 10.5147174421
0 <= FU_ERGSTt <= 1
0 <= BU_ERGSTt <= 1
-1000 <= DG_ERGSTter <= 1000
6.51471744211 <= DGo_ERGSTter <= 10.5147174421
0 <= FU_ERGSTter <= 1
0 <= BU_ERGSTter <= 1
-1000 <= DG_ETOHt <= 1000
-0.838902166985 <= DGo_ETOHt <= 3.16109783301
0 <= FU_ETOHt <= 1
0 <= BU_ETOHt <= 1
-1000 <= DG_ETOHtm <= 1000
-2.0 <= DGo_ETOHtm <= 2.0
0 <= FU_ETOHtm <= 1
0 <= BU_ETOHtm <= 1
-1000 <= DG_FBA <= 1000
2.85786845753 <= DGo_FBA <= 5.32054816874
0 <= FU_FBA <= 1
0 <= BU_FBA <= 1
-1000 <= DG_FBA3 <= 1000
2.27728470129 <= DGo_FBA3 <= 3.03653917773
0 <= FU_FBA3 <= 1
0 <= BU_FBA3 <= 1
-1000 <= DG_FBP <= 1000
-0.970246492784 <= DGo_FBP <= 0.501020987996
0 <= FU_FBP <= 1
0 <= BU_FBP <= 1
-1000 <= DG_FDH <= 1000
-4.73901690126 <= DGo_FDH <= -2.21390145836
0 <= FU_FDH <= 1
0 <= BU_FDH <= 1
-1000 <= DG_FDNG <= 1000
-21.4951272287 <= DGo_FDNG <= -16.0362819295
0 <= FU_FDNG <= 1
0 <= BU_FDNG <= 1
-1000 <= DG_FORt <= 1000
-0.515806907414 <= DGo_FORt <= 3.48419309259
0 <= FU_FORt <= 1
0 <= BU_FORt <= 1
-1000 <= DG_FRDcm <= 1000
-16.6409444938 <= DGo_FRDcm <= -12.8021158391
0 <= FU_FRDcm <= 1
0 <= BU_FRDcm <= 1
-1000 <= DG_FRDm <= 1000
-15.2732908525 <= DGo_FRDm <= -11.4344621978
0 <= FU_FRDm <= 1
0 <= BU_FRDm <= 1
-1000 <= DG_FUM <= 1000
-1.48799804487 <= DGo_FUM <= 0.257453002853
0 <= FU_FUM <= 1
0 <= BU_FUM <= 1
-1000 <= DG_FUMm <= 1000
-1.48644558836 <= DGo_FUMm <= 0.259005459369
0 <= FU_FUMm <= 1
0 <= BU_FUMm <= 1
-1000 <= DG_FUMt2r <= 1000
-3.53019177534 <= DGo_FUMt2r <= 0.469808224659
0 <= FU_FUMt2r <= 1
0 <= BU_FUMt2r <= 1
-1000 <= DG_G6PDH2 <= 1000
-4.35335229512 <= DGo_G6PDH2 <= -1.13493218581
0 <= FU_G6PDH2 <= 1
0 <= BU_G6PDH2 <= 1
-1000 <= DG_G6PDH2er <= 1000
-3.28127740622 <= DGo_G6PDH2er <= -0.0628572969084
0 <= FU_G6PDH2er <= 1
0 <= BU_G6PDH2er <= 1
-1000 <= DG_G6PI <= 1000
-1.5696912011 <= DGo_G6PI <= 2.4303087989
0 <= FU_G6PI <= 1
0 <= BU_G6PI <= 1
-1000 <= DG_G6PI3 <= 1000
-2.12207759428 <= DGo_G6PI3 <= -0.904320902649
0 <= FU_G6PI3 <= 1
0 <= BU_G6PI3 <= 1
-1000 <= DG_G6Pter <= 1000
-3.45213826923 <= DGo_G6Pter <= 0.547861730772
0 <= FU_G6Pter <= 1
0 <= BU_G6Pter <= 1
-1000 <= DG_GAPD <= 1000
-2.70951368986 <= DGo_GAPD <= -0.319476870056
0 <= FU_GAPD <= 1
0 <= BU_GAPD <= 1
-1000 <= DG_GCC2cm <= 1000
-2.22211940435 <= DGo_GCC2cm <= 2.37994637528
0 <= FU_GCC2cm <= 1
0 <= BU_GCC2cm <= 1
-1000 <= DG_GLCt1 <= 1000
0.322195666029 <= DGo_GLCt1 <= 4.32219566603
0 <= FU_GLCt1 <= 1
0 <= BU_GLCt1 <= 1
-1000 <= DG_GLNt2r <= 1000
-4.17516928965 <= DGo_GLNt2r <= -0.175169289649
0 <= FU_GLNt2r <= 1
0 <= BU_GLNt2r <= 1
-1000 <= DG_GLUDC <= 1000
-2.49055871979 <= DGo_GLUDC <= -1.4071588121
0 <= FU_GLUDC <= 1
0 <= BU_GLUDC <= 1
-1000 <= DG_GLUDxi <= 1000
10.3095081319 <= DGo_GLUDxi <= 12.9338297585
0 <= FU_GLUDxi <= 1
0 <= BU_GLUDxi <= 1
-1000 <= DG_GLUDyi <= 1000
-12.5467971475 <= DGo_GLUDyi <= -9.92247552087
0 <= FU_GLUDyi <= 1
0 <= BU_GLUDyi <= 1
-1000 <= DG_GLUK <= 1000
-5.28625518608 <= DGo_GLUK <= -4.81279557851
0 <= FU_GLUK <= 1
0 <= BU_GLUK <= 1
-1000 <= DG_GLUt2r <= 1000
-6.77198770934 <= DGo_GLUt2r <= -2.77198770934
0 <= FU_GLUt2r <= 1
0 <= BU_GLUt2r <= 1
-1000 <= DG_GLUt7m <= 1000
1.46951597642 <= DGo_GLUt7m <= 5.46951597642
0 <= FU_GLUt7m <= 1
0 <= BU_GLUt7m <= 1
-1000 <= DG_GLYt2r <= 1000
-5.14275081716 <= DGo_GLYt2r <= -1.14275081716
0 <= FU_GLYt2r <= 1
0 <= BU_GLYt2r <= 1
-1000 <= DG_GND <= 1000
-0.709932431852 <= DGo_GND <= 1.74444638162
0 <= FU_GND <= 1
0 <= BU_GND <= 1
0 <= FU_H2Ot <= 1
0 <= BU_H2Ot <= 1
0 <= FU_H2Oter <= 1
0 <= BU_H2Oter <= 1
0 <= FU_H2Otm <= 1
0 <= BU_H2Otm <= 1
-1000 <= DG_HCITSm <= 1000
-10.0837033284 <= DGo_HCITSm <= -7.46648584412
0 <= FU_HCITSm <= 1
0 <= BU_HCITSm <= 1
-1000 <= DG_HCO3E <= 1000
-2.17756006993 <= DGo_HCO3E <= -0.445509262362
0 <= FU_HCO3E <= 1
0 <= BU_HCO3E <= 1
-1000 <= DG_HEX1 <= 1000
-5.71656398499 <= DGo_HEX1 <= -5.24310437742
0 <= FU_HEX1 <= 1
0 <= BU_HEX1 <= 1
-1000 <= DG_HISt2r <= 1000
-4.36868559515 <= DGo_HISt2r <= -0.368685595152
0 <= FU_HISt2r <= 1
0 <= BU_HISt2r <= 1
-1000 <= DG_ICDHxm <= 1000
2.28282921784 <= DGo_ICDHxm <= 4.76434717264
0 <= FU_ICDHxm <= 1
0 <= BU_ICDHxm <= 1
-1000 <= DG_ICDHy <= 1000
1.90974702024 <= DGo_ICDHy <= 4.39126497504
0 <= FU_ICDHy <= 1
0 <= BU_ICDHy <= 1
-1000 <= DG_ICDHym <= 1000
1.89579660684 <= DGo_ICDHym <= 4.37731456164
0 <= FU_ICDHym <= 1
0 <= BU_ICDHym <= 1
-1000 <= DG_ILEt2r <= 1000
-3.59462037314 <= DGo_ILEt2r <= 0.405379626858
0 <= FU_ILEt2r <= 1
0 <= BU_ILEt2r <= 1
-1000 <= DG_ILEtmi <= 1000
-2.0 <= DGo_ILEtmi <= 2.0
0 <= FU_ILEtmi <= 1
0 <= BU_ILEtmi <= 1
-1000 <= DG_LEUt2r <= 1000
-3.59462037314 <= DGo_LEUt2r <= 0.405379626858
0 <= FU_LEUt2r <= 1
0 <= BU_LEUt2r <= 1
-1000 <= DG_LYSt2r <= 1000
-4.78476406764 <= DGo_LYSt2r <= -0.78476406764
0 <= FU_LYSt2r <= 1
0 <= BU_LYSt2r <= 1
-1000 <= DG_MALt2r <= 1000
-2.72611951963 <= DGo_MALt2r <= 1.27388048037
0 <= FU_MALt2r <= 1
0 <= BU_MALt2r <= 1
-1000 <= DG_MALtm <= 1000
-2.13908227298 <= DGo_MALtm <= 1.86091772702
0 <= FU_MALtm <= 1
0 <= BU_MALtm <= 1
-1000 <= DG_MAN6PI <= 1000
-2.12207759428 <= DGo_MAN6PI <= -0.904320902649
0 <= FU_MAN6PI <= 1
0 <= BU_MAN6PI <= 1
-1000 <= DG_MANNANter <= 1000
-8.19252177608 <= DGo_MANNANter <= -4.19252177608
0 <= FU_MANNANter <= 1
0 <= BU_MANNANter <= 1
-1000 <= DG_MDH <= 1000
4.33934197271 <= DGo_MDH <= 6.52952552666
0 <= FU_MDH <= 1
0 <= BU_MDH <= 1
-1000 <= DG_MDHm <= 1000
3.65640870981 <= DGo_MDHm <= 5.84659226376
0 <= FU_MDHm <= 1
0 <= BU_MDHm <= 1
-1000 <= DG_ME1m <= 1000
0.103902214234 <= DGo_ME1m <= 2.53794152369
0 <= FU_ME1m <= 1
0 <= BU_ME1m <= 1
-1000 <= DG_ME2m <= 1000
-0.283130396771 <= DGo_ME2m <= 2.15090891269
0 <= FU_ME2m <= 1
0 <= BU_ME2m <= 1
-1000 <= DG_METt2r <= 1000
-3.98165298415 <= DGo_METt2r <= 0.0183470158531
0 <= FU_METt2r <= 1
0 <= BU_METt2r <= 1
-1000 <= DG_NADH2-u6cm <= 1000
-18.1254047742 <= DGo_NADH2-u6cm <= -12.4530860244
0 <= FU_NADH2-u6cm <= 1
0 <= BU_NADH2-u6cm <= 1
-1000 <= DG_NADH2-u6m <= 1000
-18.8070728604 <= DGo_NADH2-u6m <= -13.1347541105
0 <= FU_NADH2-u6m <= 1
0 <= BU_NADH2-u6m <= 1
-1000 <= DG_NADK <= 1000
-4.49806110571 <= DGo_NADK <= -4.02460149814
0 <= FU_NADK <= 1
0 <= BU_NADK <= 1
-1000 <= DG_NADKm <= 1000
-4.52162897788 <= DGo_NADKm <= -4.04816937031
0 <= FU_NADKm <= 1
0 <= BU_NADKm <= 1
-1000 <= DG_NADPPPS <= 1000
-3.67276305366 <= DGo_NADPPPS <= -2.20149557288
0 <= FU_NADPPPS <= 1
0 <= BU_NADPPPS <= 1
-1000 <= DG_NADPPPSm <= 1000
-4.21336026014 <= DGo_NADPPPSm <= -2.74209277936
0 <= FU_NADPPPSm <= 1
0 <= BU_NADPPPSm <= 1
-1000 <= DG_NDP1 <= 1000
-8.34705370612 <= DGo_NDP1 <= -6.91194534784
0 <= FU_NDP1 <= 1
0 <= BU_NDP1 <= 1
-1000 <= DG_NH4t <= 1000
-2.80311108349 <= DGo_NH4t <= 1.19688891651
0 <= FU_NH4t <= 1
0 <= BU_NH4t <= 1
-1000 <= DG_NH4tm <= 1000
-6.15098 <= DGo_NH4tm <= -2.15098
0 <= FU_NH4tm <= 1
0 <= BU_NH4tm <= 1
-1000 <= DG_O2t <= 1000
-2.0 <= DGo_O2t <= 2.0
0 <= FU_O2t <= 1
0 <= BU_O2t <= 1
-1000 <= DG_O2ter <= 1000
-2.0 <= DGo_O2ter <= 2.0
0 <= FU_O2ter <= 1
0 <= BU_O2ter <= 1
-1000 <= DG_O2tm <= 1000
-2.0 <= DGo_O2tm <= 2.0
0 <= FU_O2tm <= 1
0 <= BU_O2tm <= 1
-1000 <= DG_OAAt <= 1000
-4.67955769725 <= DGo_OAAt <= -0.679557697254
0 <= FU_OAAt <= 1
0 <= BU_OAAt <= 1
-1000 <= DG_OAAt2m <= 1000
1.47003534381 <= DGo_OAAt2m <= 5.47003534381
0 <= FU_OAAt2m <= 1
0 <= BU_OAAt2m <= 1
-1000 <= DG_PC <= 1000
-2.51485844788 <= DGo_PC <= -1.03594690219
0 <= FU_PC <= 1
0 <= BU_PC <= 1
-1000 <= DG_PDHcm <= 1000
-2.22211940435 <= DGo_PDHcm <= 2.37994637528
0 <= FU_PDHcm <= 1
0 <= BU_PDHcm <= 1
-1000 <= DG_PDHm <= 1000
-12.1291865858 <= DGo_PDHm <= -8.93379052269
0 <= FU_PDHm <= 1
0 <= BU_PDHm <= 1
0 <= FU_PEtm_SC <= 1
0 <= BU_PEtm_SC <= 1
-1000 <= DG_PFK <= 1000
-7.20057766658 <= DGo_PFK <= -6.72711805902
0 <= FU_PFK <= 1
0 <= BU_PFK <= 1
-1000 <= DG_PFK_3 <= 1000
-6.7623238335 <= DGo_PFK_3 <= -6.28886422593
0 <= FU_PFK_3 <= 1
0 <= BU_PFK_3 <= 1
-1000 <= DG_PGI <= 1000
-1.69176879537 <= DGo_PGI <= -0.474012103745
0 <= FU_PGI <= 1
0 <= BU_PGI <= 1
-1000 <= DG_PGK <= 1000
0.489750617092 <= DGo_PGK <= 1.07052940458
0 <= FU_PGK <= 1
0 <= BU_PGK <= 1
-1000 <= DG_PGL <= 1000
-7.92548456248 <= DGo_PGL <= -5.08393146602
0 <= FU_PGL <= 1
0 <= BU_PGL <= 1
-1000 <= DG_PGM <= 1000
-2.09723204 <= DGo_PGM <= 1.90276796
0 <= FU_PGM <= 1
0 <= BU_PGM <= 1
-1000 <= DG_PGMT <= 1000
-2.10115769483 <= DGo_PGMT <= 1.89884230517
0 <= FU_PGMT <= 1
0 <= BU_PGMT <= 1
-1000 <= DG_PHEt2r <= 1000
-3.98165298415 <= DGo_PHEt2r <= 0.0183470158532
0 <= FU_PHEt2r <= 1
0 <= BU_PHEt2r <= 1
-1000 <= DG_PIt2m <= 1000
1.61038279352 <= DGo_PIt2m <= 5.61038279352
0 <= FU_PIt2m <= 1
0 <= BU_PIt2m <= 1
-1000 <= DG_PIt2r <= 1000
-1.11314146895 <= DGo_PIt2r <= 2.88685853105
0 <= FU_PIt2r <= 1
0 <= BU_PIt2r <= 1
-1000 <= DG_PPA <= 1000
-6.65796480856 <= DGo_PPA <= -4.20847506577
0 <= FU_PPA <= 1
0 <= BU_PPA <= 1
-1000 <= DG_PPAm <= 1000
-7.08604463446 <= DGo_PPAm <= -4.63655489168
0 <= FU_PPAm <= 1
0 <= BU_PPAm <= 1
-1000 <= DG_PPCK <= 1000
-1.33041791598 <= DGo_PPCK <= 0.326402943373
0 <= FU_PPCK <= 1
0 <= BU_PPCK <= 1
-1000 <= DG_PROt2r <= 1000
-3.37205820616 <= DGo_PROt2r <= 0.627941793843
0 <= FU_PROt2r <= 1
0 <= BU_PROt2r <= 1
0 <= FU_PStm_SC <= 1
0 <= BU_PStm_SC <= 1
-1000 <= DG_PYK <= 1000
-4.26382076247 <= DGo_PYK <= -2.95521081294
0 <= FU_PYK <= 1
0 <= BU_PYK <= 1
-1000 <= DG_PYRDC <= 1000
-3.80778799109 <= DGo_PYRDC <= -2.27144037418
0 <= FU_PYRDC <= 1
0 <= BU_PYRDC <= 1
-1000 <= DG_PYRt <= 1000
-3.770692611 <= DGo_PYRt <= 0.229307388995
0 <= FU_PYRt <= 1
0 <= BU_PYRt <= 1
-1000 <= DG_PYRt2m <= 1000
-2.68166808617 <= DGo_PYRt2m <= 1.31833191383
0 <= FU_PYRt2m <= 1
0 <= BU_PYRt2m <= 1
-1000 <= DG_RPE <= 1000
-2.0 <= DGo_RPE <= 2.0
0 <= FU_RPE <= 1
0 <= BU_RPE <= 1
-1000 <= DG_RPI <= 1000
-0.586270264425 <= DGo_RPI <= 2.01668562456
0 <= FU_RPI <= 1
0 <= BU_RPI <= 1
-1000 <= DG_SERt2r <= 1000
-4.75571820616 <= DGo_SERt2r <= -0.755718206157
0 <= FU_SERt2r <= 1
0 <= BU_SERt2r <= 1
-1000 <= DG_SO4ti <= 1000
-0.00674522200977 <= DGo_SO4ti <= 3.99325477799
0 <= FU_SO4ti <= 1
0 <= BU_SO4ti <= 1
-1000 <= DG_SUCCt2r <= 1000
-2.2819233204 <= DGo_SUCCt2r <= 1.7180766796
0 <= FU_SUCCt2r <= 1
0 <= BU_SUCCt2r <= 1
-1000 <= DG_SUCCtm <= 1000
-2.13631726053 <= DGo_SUCCtm <= 1.86368273947
0 <= FU_SUCCtm <= 1
0 <= BU_SUCCtm <= 1
-1000 <= DG_SUCD1m <= 1000
11.4344621978 <= DGo_SUCD1m <= 15.2732908525
0 <= FU_SUCD1m <= 1
0 <= BU_SUCD1m <= 1
-1000 <= DG_SUCD2_u6m <= 1000
-4.99213387661 <= DGo_SUCD2_u6m <= 0.116149341852
0 <= FU_SUCD2_u6m <= 1
0 <= BU_SUCD2_u6m <= 1
-1000 <= DG_SUCD3-u6m <= 1000
-19.0119321459 <= DGo_SUCD3-u6m <= -12.5718054391
0 <= FU_SUCD3-u6m <= 1
0 <= BU_SUCD3-u6m <= 1
-1000 <= DG_SUCFUMtm <= 1000
-1.99568253104 <= DGo_SUCFUMtm <= 2.00431746896
0 <= FU_SUCFUMtm <= 1
0 <= BU_SUCFUMtm <= 1
-1000 <= DG_SUCOASm <= 1000
-5.30090974893 <= DGo_SUCOASm <= -2.87407477308
0 <= FU_SUCOASm <= 1
0 <= BU_SUCOASm <= 1
-1000 <= DG_TALA <= 1000
-2.27972873866 <= DGo_TALA <= 0.291643657556
0 <= FU_TALA <= 1
0 <= BU_TALA <= 1
-1000 <= DG_THRt2m <= 1000
-6.83264808617 <= DGo_THRt2m <= -2.83264808617
0 <= FU_THRt2m <= 1
0 <= BU_THRt2m <= 1
-1000 <= DG_THRt2r <= 1000
-4.36868559515 <= DGo_THRt2r <= -0.368685595152
0 <= FU_THRt2r <= 1
0 <= BU_THRt2r <= 1
-1000 <= DG_TKT1 <= 1000
-0.544723803861 <= DGo_TKT1 <= 2.13682574817
0 <= FU_TKT1 <= 1
0 <= BU_TKT1 <= 1
-1000 <= DG_TKT2 <= 1000
-5.03804215448 <= DGo_TKT2 <= -2.46666975827
0 <= FU_TKT2 <= 1
0 <= BU_TKT2 <= 1
-1000 <= DG_TPI <= 1000
1.15931105774 <= DGo_TPI <= 1.76068894226
0 <= FU_TPI <= 1
0 <= BU_TPI <= 1
-1000 <= DG_TRPt2r <= 1000
-3.78813667864 <= DGo_TRPt2r <= 0.211863321356
0 <= FU_TRPt2r <= 1
0 <= BU_TRPt2r <= 1
-1000 <= DG_TYRt2r <= 1000
-3.98165298415 <= DGo_TYRt2r <= 0.0183470158532
0 <= FU_TYRt2r <= 1
0 <= BU_TYRt2r <= 1
-1000 <= DG_VALt2r <= 1000
-3.98165298415 <= DGo_VALt2r <= 0.0183470158532
0 <= FU_VALt2r <= 1
0 <= BU_VALt2r <= 1
-1000 <= DG_ZYMSTt <= 1000
6.51471744211 <= DGo_ZYMSTt <= 10.5147174421
0 <= FU_ZYMSTt <= 1
0 <= BU_ZYMSTt <= 1
0 <= FU_LMPD_1_13BDglcn_c <= 1
0 <= BU_LMPD_1_13BDglcn_c <= 1
-1000 <= DG_LMPD_2_amp_c <= 1000
-22.0724501305 <= DGo_LMPD_2_amp_c <= -7.41304515483
0 <= FU_LMPD_2_amp_c <= 1
0 <= BU_LMPD_2_amp_c <= 1
-1000 <= DG_LMPD_3_amp_c <= 1000
-21.2983849085 <= DGo_LMPD_3_amp_c <= -6.63897993282
0 <= FU_LMPD_3_amp_c <= 1
0 <= BU_LMPD_3_amp_c <= 1
-1000 <= DG_LMPD_4_amp_c <= 1000
-22.5034890423 <= DGo_LMPD_4_amp_c <= -7.84408406662
0 <= FU_LMPD_4_amp_c <= 1
0 <= BU_LMPD_4_amp_c <= 1
-1000 <= DG_LMPD_5_amp_c <= 1000
-21.7294238203 <= DGo_LMPD_5_amp_c <= -7.07001884461
0 <= FU_LMPD_5_amp_c <= 1
0 <= BU_LMPD_5_amp_c <= 1
-1000 <= DG_LMPD_6_arg-L_c <= 1000
-10.8104595959 <= DGo_LMPD_6_arg-L_c <= -4.57540899471
0 <= FU_LMPD_6_arg-L_c <= 1
0 <= BU_LMPD_6_arg-L_c <= 1
-1000 <= DG_LMPD_7_asn-L_c <= 1000
-9.31749296261 <= DGo_LMPD_7_asn-L_c <= -6.94317874583
0 <= FU_LMPD_7_asn-L_c <= 1
0 <= BU_LMPD_7_asn-L_c <= 1
-1000 <= DG_LMPD_8_cmp_c <= 1000
-48.515606034 <= DGo_LMPD_8_cmp_c <= -40.3612999716
0 <= FU_LMPD_8_cmp_c <= 1
0 <= BU_LMPD_8_cmp_c <= 1
-1000 <= DG_LMPD_9_cmp_c <= 1000
-46.772669684 <= DGo_LMPD_9_cmp_c <= -39.9635590695
0 <= FU_LMPD_9_cmp_c <= 1
0 <= BU_LMPD_9_cmp_c <= 1
-1000 <= DG_LMPD_10_cmp_c <= 1000
-40.972241688 <= DGo_LMPD_10_cmp_c <= -33.5077430873
0 <= FU_LMPD_10_cmp_c <= 1
0 <= BU_LMPD_10_cmp_c <= 1
-1000 <= DG_LMPD_11_cmp_c <= 1000
-39.1525223205 <= DGo_LMPD_11_cmp_c <= -33.1867852026
0 <= FU_LMPD_11_cmp_c <= 1
0 <= BU_LMPD_11_cmp_c <= 1
-1000 <= DG_LMPD_12_cys-L_c <= 1000
-48.6722025076 <= DGo_LMPD_12_cys-L_c <= -41.0651991554
0 <= FU_LMPD_12_cys-L_c <= 1
0 <= BU_LMPD_12_cys-L_c <= 1
-1000 <= DG_LMPD_13_damp_c <= 1000
-32.5060287755 <= DGo_LMPD_13_damp_c <= -19.0237845172
0 <= FU_LMPD_13_damp_c <= 1
0 <= BU_LMPD_13_damp_c <= 1
-1000 <= DG_LMPD_14_damp_c <= 1000
-32.9370676873 <= DGo_LMPD_14_damp_c <= -19.454823429
0 <= FU_LMPD_14_damp_c <= 1
0 <= BU_LMPD_14_damp_c <= 1
-1000 <= DG_LMPD_15_damp_c <= 1000
-33.2800939975 <= DGo_LMPD_15_damp_c <= -19.7978497392
0 <= FU_LMPD_15_damp_c <= 1
0 <= BU_LMPD_15_damp_c <= 1
-1000 <= DG_LMPD_16_damp_c <= 1000
-33.7111329093 <= DGo_LMPD_16_damp_c <= -20.228888651
0 <= FU_LMPD_16_damp_c <= 1
0 <= BU_LMPD_16_damp_c <= 1
-1000 <= DG_LMPD_17_dcmp_c <= 1000
-43.3321528382 <= DGo_LMPD_17_dcmp_c <= -37.4353313757
0 <= FU_LMPD_17_dcmp_c <= 1
0 <= BU_LMPD_17_dcmp_c <= 1
-1000 <= DG_LMPD_18_dcmp_c <= 1000
-45.1502682066 <= DGo_LMPD_18_dcmp_c <= -37.7578932594
0 <= FU_LMPD_18_dcmp_c <= 1
0 <= BU_LMPD_18_dcmp_c <= 1
-1000 <= DG_LMPD_19_dgmp_c <= 1000
-54.4194762275 <= DGo_LMPD_19_dgmp_c <= -38.3717560954
0 <= FU_LMPD_19_dgmp_c <= 1
0 <= BU_LMPD_19_dgmp_c <= 1
-1000 <= DG_LMPD_20_dgmp_c <= 1000
-46.7496384731 <= DGo_LMPD_20_dgmp_c <= -31.6446726194
0 <= FU_LMPD_20_dgmp_c <= 1
0 <= BU_LMPD_20_dgmp_c <= 1
-1000 <= DG_LMPD_21_dgmp_c <= 1000
-47.5237036951 <= DGo_LMPD_21_dgmp_c <= -32.4187378414
0 <= FU_LMPD_21_dgmp_c <= 1
0 <= BU_LMPD_21_dgmp_c <= 1
-1000 <= DG_LMPD_22_dgmp_c <= 1000
-55.1935414495 <= DGo_LMPD_22_dgmp_c <= -39.1458213174
0 <= FU_LMPD_22_dgmp_c <= 1
0 <= BU_LMPD_22_dgmp_c <= 1
-1000 <= DG_LMPD_23_dtmp_c <= 1000
-58.1570297347 <= DGo_LMPD_23_dtmp_c <= -50.947093857
0 <= FU_LMPD_23_dtmp_c <= 1
0 <= BU_LMPD_23_dtmp_c <= 1
-1000 <= DG_LMPD_24_dtmp_c <= 1000
-56.4248639819 <= DGo_LMPD_24_dtmp_c <= -50.5385823576
0 <= FU_LMPD_24_dtmp_c <= 1
0 <= BU_LMPD_24_dtmp_c <= 1
-1000 <= DG_LMPD_25_ergst_c <= 1000
-1413.32563784 <= DGo_LMPD_25_ergst_c <= -1340.74716751
0 <= FU_LMPD_25_ergst_c <= 1
0 <= BU_LMPD_25_ergst_c <= 1
-1000 <= DG_LMPD_26_ergst_c <= 1000
-1414.48673568 <= DGo_LMPD_26_ergst_c <= -1341.90826534
0 <= FU_LMPD_26_ergst_c <= 1
0 <= BU_LMPD_26_ergst_c <= 1
-1000 <= DG_LMPD_27_gln-L_c <= 1000
1.51887211296 <= DGo_LMPD_27_gln-L_c <= 3.2968147442
0 <= FU_LMPD_27_gln-L_c <= 1
0 <= BU_LMPD_27_gln-L_c <= 1
-1000 <= DG_LMPD_28_gly_c <= 1000
3.7791485331 <= DGo_LMPD_28_gly_c <= 4.91630526811
0 <= FU_LMPD_28_gly_c <= 1
0 <= BU_LMPD_28_gly_c <= 1
0 <= FU_LMPD_29_glycogen_c <= 1
0 <= BU_LMPD_29_glycogen_c <= 1
-1000 <= DG_LMPD_30_gmp_c <= 1000
-36.0785462623 <= DGo_LMPD_30_gmp_c <= -19.490666909
0 <= FU_LMPD_30_gmp_c <= 1
0 <= BU_LMPD_30_gmp_c <= 1
-1000 <= DG_LMPD_31_gmp_c <= 1000
-36.8526114843 <= DGo_LMPD_31_gmp_c <= -20.264732131
0 <= FU_LMPD_31_gmp_c <= 1
0 <= BU_LMPD_31_gmp_c <= 1
-1000 <= DG_LMPD_32_his-L_c <= 1000
-52.8877489465 <= DGo_LMPD_32_his-L_c <= -42.4264198072
0 <= FU_LMPD_32_his-L_c <= 1
0 <= BU_LMPD_32_his-L_c <= 1
-1000 <= DG_LMPD_33_his-L_c <= 1000
-52.0696774237 <= DGo_LMPD_33_his-L_c <= -41.6083482845
0 <= FU_LMPD_33_his-L_c <= 1
0 <= BU_LMPD_33_his-L_c <= 1
-1000 <= DG_LMPD_34_his-L_c <= 1000
-52.4567100347 <= DGo_LMPD_34_his-L_c <= -41.9953808955
0 <= FU_LMPD_34_his-L_c <= 1
0 <= BU_LMPD_34_his-L_c <= 1
-1000 <= DG_LMPD_35_his-L_c <= 1000
-52.5007163355 <= DGo_LMPD_35_his-L_c <= -42.0393871962
0 <= FU_LMPD_35_his-L_c <= 1
0 <= BU_LMPD_35_his-L_c <= 1
-1000 <= DG_LMPD_36_ile-L_c <= 1000
-23.5373070498 <= DGo_LMPD_36_ile-L_c <= -16.75333034
0 <= FU_LMPD_36_ile-L_c <= 1
0 <= BU_LMPD_36_ile-L_c <= 1
-1000 <= DG_LMPD_37_ile-L_c <= 1000
-23.9243396608 <= DGo_LMPD_37_ile-L_c <= -17.140362951
0 <= FU_LMPD_37_ile-L_c <= 1
0 <= BU_LMPD_37_ile-L_c <= 1
-1000 <= DG_LMPD_38_leu-L_c <= 1000
-15.4823085306 <= DGo_LMPD_38_leu-L_c <= -12.3960902363
0 <= FU_LMPD_38_leu-L_c <= 1
0 <= BU_LMPD_38_leu-L_c <= 1
-1000 <= DG_LMPD_39_leu-L_c <= 1000
-26.1871275665 <= DGo_LMPD_39_leu-L_c <= -23.1009092722
0 <= FU_LMPD_39_leu-L_c <= 1
0 <= BU_LMPD_39_leu-L_c <= 1
-1000 <= DG_LMPD_40_lys-L_c <= 1000
11.4317837095 <= DGo_LMPD_40_lys-L_c <= 13.4399409145
0 <= FU_LMPD_40_lys-L_c <= 1
0 <= BU_LMPD_40_lys-L_c <= 1
-1000 <= DG_LMPD_41_lys-L_c <= 1000
11.0447510985 <= DGo_LMPD_41_lys-L_c <= 13.0529083035
0 <= FU_LMPD_41_lys-L_c <= 1
0 <= BU_LMPD_41_lys-L_c <= 1
0 <= FU_LMPD_42_mannan_c <= 1
0 <= BU_LMPD_42_mannan_c <= 1
0 <= FU_LMPD_43_mannan_c <= 1
0 <= BU_LMPD_43_mannan_c <= 1
-1000 <= DG_LMPD_44_met-L_c <= 1000
-50.5911689628 <= DGo_LMPD_44_met-L_c <= -37.2603768123
0 <= FU_LMPD_44_met-L_c <= 1
0 <= BU_LMPD_44_met-L_c <= 1
-1000 <= DG_LMPD_45_met-L_c <= 1000
-50.2041363518 <= DGo_LMPD_45_met-L_c <= -36.8733442013
0 <= FU_LMPD_45_met-L_c <= 1
0 <= BU_LMPD_45_met-L_c <= 1
0 <= FU_LMPD_47_pa_SC_c <= 1
0 <= BU_LMPD_47_pa_SC_c <= 1
0 <= FU_LMPD_48_pc_SC_c <= 1
0 <= BU_LMPD_48_pc_SC_c <= 1
0 <= FU_LMPD_49_pc_SC_c <= 1
0 <= BU_LMPD_49_pc_SC_c <= 1
0 <= FU_LMPD_51_pe_SC_c <= 1
0 <= BU_LMPD_51_pe_SC_c <= 1
-1000 <= DG_LMPD_52_phe-L_c <= 1000
-97.3224585258 <= DGo_LMPD_52_phe-L_c <= -91.6863250994
0 <= FU_LMPD_52_phe-L_c <= 1
0 <= BU_LMPD_52_phe-L_c <= 1
-1000 <= DG_LMPD_53_pro-L_c <= 1000
-6.27207771765 <= DGo_LMPD_53_pro-L_c <= -1.73656171172
0 <= FU_LMPD_53_pro-L_c <= 1
0 <= BU_LMPD_53_pro-L_c <= 1
-1000 <= DG_LMPD_54_pro-L_c <= 1000
-5.88504510664 <= DGo_LMPD_54_pro-L_c <= -1.34952910072
0 <= FU_LMPD_54_pro-L_c <= 1
0 <= BU_LMPD_54_pro-L_c <= 1
0 <= FU_LMPD_56_ps_SC_c <= 1
0 <= BU_LMPD_56_ps_SC_c <= 1
0 <= FU_LMPD_57_ptd1ino_SC_c <= 1
0 <= BU_LMPD_57_ptd1ino_SC_c <= 1
0 <= FU_LMPD_58_ptd1ino_SC_c <= 1
0 <= BU_LMPD_58_ptd1ino_SC_c <= 1
-1000 <= DG_LMPD_59_ribflv_c <= 1000
-132.196559443 <= DGo_LMPD_59_ribflv_c <= -117.012484385
0 <= FU_LMPD_59_ribflv_c <= 1
0 <= BU_LMPD_59_ribflv_c <= 1
-1000 <= DG_LMPD_60_ribflv_c <= 1000
-131.765520531 <= DGo_LMPD_60_ribflv_c <= -116.581445473
0 <= FU_LMPD_60_ribflv_c <= 1
0 <= BU_LMPD_60_ribflv_c <= 1
-1000 <= DG_LMPD_61_ser-L_c <= 1000
5.96706856345 <= DGo_LMPD_61_ser-L_c <= 8.58136081949
0 <= FU_LMPD_61_ser-L_c <= 1
0 <= BU_LMPD_61_ser-L_c <= 1
-1000 <= DG_LMPD_62_thr-L_c <= 1000
-8.30819597227 <= DGo_LMPD_62_thr-L_c <= -3.50947120834
0 <= FU_LMPD_62_thr-L_c <= 1
0 <= BU_LMPD_62_thr-L_c <= 1
-1000 <= DG_LMPD_63_thr-L_c <= 1000
-7.92116336127 <= DGo_LMPD_63_thr-L_c <= -3.12243859733
0 <= FU_LMPD_63_thr-L_c <= 1
0 <= BU_LMPD_63_thr-L_c <= 1
-1000 <= DG_LMPD_64_tre_c <= 1000
-3.2067575345 <= DGo_LMPD_64_tre_c <= -1.1745449472
0 <= FU_LMPD_64_tre_c <= 1
0 <= BU_LMPD_64_tre_c <= 1
0 <= FU_LMPD_65_triglyc_SC_c <= 1
0 <= BU_LMPD_65_triglyc_SC_c <= 1
0 <= FU_LMPD_66_triglyc_SC_c <= 1
0 <= BU_LMPD_66_triglyc_SC_c <= 1
-1000 <= DG_LMPD_67_trp-L_c <= 1000
-136.730104125 <= DGo_LMPD_67_trp-L_c <= -128.571496195
0 <= FU_LMPD_67_trp-L_c <= 1
0 <= BU_LMPD_67_trp-L_c <= 1
-1000 <= DG_LMPD_68_tyr-L_c <= 1000
-84.9259439066 <= DGo_LMPD_68_tyr-L_c <= -79.7636765753
0 <= FU_LMPD_68_tyr-L_c <= 1
0 <= BU_LMPD_68_tyr-L_c <= 1
-1000 <= DG_LMPD_69_tyr-L_c <= 1000
-85.3129765176 <= DGo_LMPD_69_tyr-L_c <= -80.1507091863
0 <= FU_LMPD_69_tyr-L_c <= 1
0 <= BU_LMPD_69_tyr-L_c <= 1
-1000 <= DG_LMPD_70_ump_c <= 1000
-36.5259648522 <= DGo_LMPD_70_ump_c <= -31.6757150277
0 <= FU_LMPD_70_ump_c <= 1
0 <= BU_LMPD_70_ump_c <= 1
-1000 <= DG_LMPD_71_ump_c <= 1000
-38.5461314309 <= DGo_LMPD_71_ump_c <= -31.7962257012
0 <= FU_LMPD_71_ump_c <= 1
0 <= BU_LMPD_71_ump_c <= 1
-1000 <= DG_LMPD_72_val-L_c <= 1000
-20.0635629167 <= DGo_LMPD_72_val-L_c <= -17.4592042624
0 <= FU_LMPD_72_val-L_c <= 1
0 <= BU_LMPD_72_val-L_c <= 1
-1000 <= DG_LMPD_73_zymst_c <= 1000
-1264.11708091 <= DGo_LMPD_73_zymst_c <= -1197.27025665
0 <= FU_LMPD_73_zymst_c <= 1
0 <= BU_LMPD_73_zymst_c <= 1
-1000 <= DG_LMPD_74_zymst_c <= 1000
-1262.95598307 <= DGo_LMPD_74_zymst_c <= -1196.10915882
0 <= FU_LMPD_74_zymst_c <= 1
0 <= BU_LMPD_74_zymst_c <= 1
Binaries
FU_biomass_SC5_notrace	BU_biomass_SC5_notrace	FU_EX_4abut_e	BU_EX_4abut_e	FU_EX_ac_e	BU_EX_ac_e	
FU_EX_acald_e	BU_EX_acald_e	FU_EX_akg_e	BU_EX_akg_e	FU_EX_ala-L_e	BU_EX_ala-L_e	
FU_EX_arg-L_e	BU_EX_arg-L_e	FU_EX_asn-L_e	BU_EX_asn-L_e	FU_EX_asp-L_e	BU_EX_asp-L_e	
FU_EX_cit_e	BU_EX_cit_e	FU_EX_co2_e	BU_EX_co2_e	FU_EX_cys-L_e	BU_EX_cys-L_e	
FU_EX_ergst_e	BU_EX_ergst_e	FU_EX_etoh_e	BU_EX_etoh_e	FU_EX_for_e	BU_EX_for_e	
FU_EX_fum_e	BU_EX_fum_e	FU_EX_glc_e	BU_EX_glc_e	FU_EX_gln-L_e	BU_EX_gln-L_e	
FU_EX_glu-L_e	BU_EX_glu-L_e	FU_EX_gly_e	BU_EX_gly_e	FU_EX_h_e	BU_EX_h_e	
FU_EX_h2o_e	BU_EX_h2o_e	FU_EX_his-L_e	BU_EX_his-L_e	FU_EX_ile-L_e	BU_EX_ile-L_e	
FU_EX_leu-L_e	BU_EX_leu-L_e	FU_EX_lys-L_e	BU_EX_lys-L_e	FU_EX_mal-L_e	BU_EX_mal-L_e	
FU_EX_met-L_e	BU_EX_met-L_e	FU_EX_nh4_e	BU_EX_nh4_e	FU_EX_o2_e	BU_EX_o2_e	
FU_EX_oaa_e	BU_EX_oaa_e	FU_EX_phe-L_e	BU_EX_phe-L_e	FU_EX_pi_e	BU_EX_pi_e	
FU_EX_pro-L_e	BU_EX_pro-L_e	FU_EX_pyr_e	BU_EX_pyr_e	FU_EX_ser-L_e	BU_EX_ser-L_e	
FU_EX_so4_e	BU_EX_so4_e	FU_EX_succ_e	BU_EX_succ_e	FU_EX_thr-L_e	BU_EX_thr-L_e	
FU_EX_trp-L_e	BU_EX_trp-L_e	FU_EX_tyr-L_e	BU_EX_tyr-L_e	FU_EX_val-L_e	BU_EX_val-L_e	
FU_EX_zymst_e	BU_EX_zymst_e	FU_6PGLter	BU_6PGLter	FU_ABUTt2r	BU_ABUTt2r	
FU_ACALDt	BU_ACALDt	FU_ACALDtm	BU_ACALDtm	FU_ACOAHim	BU_ACOAHim	
FU_ACONT	BU_ACONT	FU_ACONTm	BU_ACONTm	FU_ACS	BU_ACS	
FU_ACSm	BU_ACSm	FU_ACt2r	BU_ACt2r	FU_ACtm	BU_ACtm	
FU_ACtr	BU_ACtr	FU_ADK1	BU_ADK1	FU_ADK1m	BU_ADK1m	
FU_AKGDam	BU_AKGDam	FU_AKGDbm	BU_AKGDbm	FU_AKGMAL	BU_AKGMAL	
FU_AKGt2r	BU_AKGt2r	FU_ALATA_L	BU_ALATA_L	FU_ALATA_Lm	BU_ALATA_Lm	
FU_ALAt2r	BU_ALAt2r	FU_ALAtmi	BU_ALAtmi	FU_ALCD2if	BU_ALCD2if	
FU_ALCD2ir	BU_ALCD2ir	FU_ALCD2irm	BU_ALCD2irm	FU_ALCD2x	BU_ALCD2x	
FU_ALDD2xm	BU_ALDD2xm	FU_ALDD2y	BU_ALDD2y	FU_ALDD2ym	BU_ALDD2ym	
FU_ARGt2r	BU_ARGt2r	FU_ASNt2r	BU_ASNt2r	FU_ASPGLU2m	BU_ASPGLU2m	
FU_ASPTA	BU_ASPTA	FU_ASPTAm	BU_ASPTAm	FU_ASPt2m	BU_ASPt2m	
FU_ASPt2r	BU_ASPt2r	FU_ATPH1	BU_ATPH1	FU_ATPM	BU_ATPM	
FU_ATPS	BU_ATPS	FU_ATPS3m	BU_ATPS3m	FU_ATPtm-H	BU_ATPtm-H	
FU_CITt2r	BU_CITt2r	FU_CITtam	BU_CITtam	FU_CITtcm	BU_CITtcm	
FU_CO2t	BU_CO2t	FU_CO2tm	BU_CO2tm	FU_CSm	BU_CSm	
FU_CYOOm	BU_CYOOm	FU_CYOR_u6m	BU_CYOR_u6m	FU_CYSt2r	BU_CYSt2r	
FU_D-LACt2m	BU_D-LACt2m	FU_D-LACtm	BU_D-LACtm	FU_DB4PS	BU_DB4PS	
FU_ENO	BU_ENO	FU_ERGSTt	BU_ERGSTt	FU_ERGSTter	BU_ERGSTter	
FU_ETOHt	BU_ETOHt	FU_ETOHtm	BU_ETOHtm	FU_FBA	BU_FBA	
FU_FBA3	BU_FBA3	FU_FBP	BU_FBP	FU_FDH	BU_FDH	
FU_FDNG	BU_FDNG	FU_FORt	BU_FORt	FU_FRDcm	BU_FRDcm	
FU_FRDm	BU_FRDm	FU_FUM	BU_FUM	FU_FUMm	BU_FUMm	
FU_FUMt2r	BU_FUMt2r	FU_G6PDH2	BU_G6PDH2	FU_G6PDH2er	BU_G6PDH2er	
FU_G6PI	BU_G6PI	FU_G6PI3	BU_G6PI3	FU_G6Pter	BU_G6Pter	
FU_GAPD	BU_GAPD	FU_GCC2cm	BU_GCC2cm	FU_GLCt1	BU_GLCt1	
FU_GLNt2r	BU_GLNt2r	FU_GLUDC	BU_GLUDC	FU_GLUDxi	BU_GLUDxi	
FU_GLUDyi	BU_GLUDyi	FU_GLUK	BU_GLUK	FU_GLUt2r	BU_GLUt2r	
FU_GLUt7m	BU_GLUt7m	FU_GLYt2r	BU_GLYt2r	FU_GND	BU_GND	
FU_H2Ot	BU_H2Ot	FU_H2Oter	BU_H2Oter	FU_H2Otm	BU_H2Otm	
FU_HCITSm	BU_HCITSm	FU_HCO3E	BU_HCO3E	FU_HEX1	BU_HEX1	
FU_HISt2r	BU_HISt2r	FU_ICDHxm	BU_ICDHxm	FU_ICDHy	BU_ICDHy	
FU_ICDHym	BU_ICDHym	FU_ILEt2r	BU_ILEt2r	FU_ILEtmi	BU_ILEtmi	
FU_LEUt2r	BU_LEUt2r	FU_LYSt2r	BU_LYSt2r	FU_MALt2r	BU_MALt2r	
FU_MALtm	BU_MALtm	FU_MAN6PI	BU_MAN6PI	FU_MANNANter	BU_MANNANter	
FU_MDH	BU_MDH	FU_MDHm	BU_MDHm	FU_ME1m	BU_ME1m	
FU_ME2m	BU_ME2m	FU_METt2r	BU_METt2r	FU_NADH2-u6cm	BU_NADH2-u6cm	
FU_NADH2-u6m	BU_NADH2-u6m	FU_NADK	BU_NADK	FU_NADKm	BU_NADKm	
FU_NADPPPS	BU_NADPPPS	FU_NADPPPSm	BU_NADPPPSm	FU_NDP1	BU_NDP1	
FU_NH4t	BU_NH4t	FU_NH4tm	BU_NH4tm	FU_O2t	BU_O2t	
FU_O2ter	BU_O2ter	FU_O2tm	BU_O2tm	FU_OAAt	BU_OAAt	
FU_OAAt2m	BU_OAAt2m	FU_PC	BU_PC	FU_PDHcm	BU_PDHcm	
FU_PDHm	BU_PDHm	FU_PEtm_SC	BU_PEtm_SC	FU_PFK	BU_PFK	
FU_PFK_3	BU_PFK_3	FU_PGI	BU_PGI	FU_PGK	BU_PGK	
FU_PGL	BU_PGL	FU_PGM	BU_PGM	FU_PGMT	BU_PGMT	
FU_PHEt2r	BU_PHEt2r	FU_PIt2m	BU_PIt2m	FU_PIt2r	BU_PIt2r	
FU_PPA	BU_PPA	FU_PPAm	BU_PPAm	FU_PPCK	BU_PPCK	
FU_PROt2r	BU_PROt2r	FU_PStm_SC	BU_PStm_SC	FU_PYK	BU_PYK	
FU_PYRDC	BU_PYRDC	FU_PYRt	BU_PYRt	FU_PYRt2m	BU_PYRt2m	
FU_RPE	BU_RPE	FU_RPI	BU_RPI	FU_SERt2r	BU_SERt2r	
FU_SO4ti	BU_SO4ti	FU_SUCCt2r	BU_SUCCt2r	FU_SUCCtm	BU_SUCCtm	
FU_SUCD1m	BU_SUCD1m	FU_SUCD2_u6m	BU_SUCD2_u6m	FU_SUCD3-u6m	BU_SUCD3-u6m	
FU_SUCFUMtm	BU_SUCFUMtm	FU_SUCOASm	BU_SUCOASm	FU_TALA	BU_TALA	
FU_THRt2m	BU_THRt2m	FU_THRt2r	BU_THRt2r	FU_TKT1	BU_TKT1	
FU_TKT2	BU_TKT2	FU_TPI	BU_TPI	FU_TRPt2r	BU_TRPt2r	
FU_TYRt2r	BU_TYRt2r	FU_VALt2r	BU_VALt2r	FU_ZYMSTt	BU_ZYMSTt	
FU_LMPD_1_13BDglcn_c	BU_LMPD_1_13BDglcn_c	FU_LMPD_2_amp_c	BU_LMPD_2_amp_c	FU_LMPD_3_amp_c	BU_LMPD_3_amp_c	
FU_LMPD_4_amp_c	BU_LMPD_4_amp_c	FU_LMPD_5_amp_c	BU_LMPD_5_amp_c	FU_LMPD_6_arg-L_c	BU_LMPD_6_arg-L_c	
FU_LMPD_7_asn-L_c	BU_LMPD_7_asn-L_c	FU_LMPD_8_cmp_c	BU_LMPD_8_cmp_c	FU_LMPD_9_cmp_c	BU_LMPD_9_cmp_c	
FU_LMPD_10_cmp_c	BU_LMPD_10_cmp_c	FU_LMPD_11_cmp_c	BU_LMPD_11_cmp_c	FU_LMPD_12_cys-L_c	BU_LMPD_12_cys-L_c	
FU_LMPD_13_damp_c	BU_LMPD_13_damp_c	FU_LMPD_14_damp_c	BU_LMPD_14_damp_c	FU_LMPD_15_damp_c	BU_LMPD_15_damp_c	
FU_LMPD_16_damp_c	BU_LMPD_16_damp_c	FU_LMPD_17_dcmp_c	BU_LMPD_17_dcmp_c	FU_LMPD_18_dcmp_c	BU_LMPD_18_dcmp_c	
FU_LMPD_19_dgmp_c	BU_LMPD_19_dgmp_c	FU_LMPD_20_dgmp_c	BU_LMPD_20_dgmp_c	FU_LMPD_21_dgmp_c	BU_LMPD_21_dgmp_c	
FU_LMPD_22_dgmp_c	BU_LMPD_22_dgmp_c	FU_LMPD_23_dtmp_c	BU_LMPD_23_dtmp_c	FU_LMPD_24_dtmp_c	BU_LMPD_24_dtmp_c	
FU_LMPD_25_ergst_c	BU_LMPD_25_ergst_c	FU_LMPD_26_ergst_c	BU_LMPD_26_ergst_c	FU_LMPD_27_gln-L_c	BU_LMPD_27_gln-L_c	
FU_LMPD_28_gly_c	BU_LMPD_28_gly_c	FU_LMPD_29_glycogen_c	BU_LMPD_29_glycogen_c	FU_LMPD_30_gmp_c	BU_LMPD_30_gmp_c	
FU_LMPD_31_gmp_c	BU_LMPD_31_gmp_c	FU_LMPD_32_his-L_c	BU_LMPD_32_his-L_c	FU_LMPD_33_his-L_c	BU_LMPD_33_his-L_c	
FU_LMPD_34_his-L_c	BU_LMPD_34_his-L_c	FU_LMPD_35_his-L_c	BU_LMPD_35_his-L_c	FU_LMPD_36_ile-L_c	BU_LMPD_36_ile-L_c	
FU_LMPD_37_ile-L_c	BU_LMPD_37_ile-L_c	FU_LMPD_38_leu-L_c	BU_LMPD_38_leu-L_c	FU_LMPD_39_leu-L_c	BU_LMPD_39_leu-L_c	
FU_LMPD_40_lys-L_c	BU_LMPD_40_lys-L_c	FU_LMPD_41_lys-L_c	BU_LMPD_41_lys-L_c	FU_LMPD_42_mannan_c	BU_LMPD_42_mannan_c	
FU_LMPD_43_mannan_c	BU_LMPD_43_mannan_c	FU_LMPD_44_met-L_c	BU_LMPD_44_met-L_c	FU_LMPD_45_met-L_c	BU_LMPD_45_met-L_c	
FU_LMPD_47_pa_SC_c	BU_LMPD_47_pa_SC_c	FU_LMPD_48_pc_SC_c	BU_LMPD_48_pc_SC_c	FU_LMPD_49_pc_SC_c	BU_LMPD_49_pc_SC_c	
FU_LMPD_51_pe_SC_c	BU_LMPD_51_pe_SC_c	FU_LMPD_52_phe-L_c	BU_LMPD_52_phe-L_c	FU_LMPD_53_pro-L_c	BU_LMPD_53_pro-L_c	
FU_LMPD_54_pro-L_c	BU_LMPD_54_pro-L_c	FU_LMPD_56_ps_SC_c	BU_LMPD_56_ps_SC_c	FU_LMPD_57_ptd1ino_SC_c	BU_LMPD_57_ptd1ino_SC_c	
FU_LMPD_58_ptd1ino_SC_c	BU_LMPD_58_ptd1ino_SC_c	FU_LMPD_59_ribflv_c	BU_LMPD_59_ribflv_c	FU_LMPD_60_ribflv_c	BU_LMPD_60_ribflv_c	
FU_LMPD_61_ser-L_c	BU_LMPD_61_ser-L_c	FU_LMPD_62_thr-L_c	BU_LMPD_62_thr-L_c	FU_LMPD_63_thr-L_c	BU_LMPD_63_thr-L_c	
FU_LMPD_64_tre_c	BU_LMPD_64_tre_c	FU_LMPD_65_triglyc_SC_c	BU_LMPD_65_triglyc_SC_c	FU_LMPD_66_triglyc_SC_c	BU_LMPD_66_triglyc_SC_c	
FU_LMPD_67_trp-L_c	BU_LMPD_67_trp-L_c	FU_LMPD_68_tyr-L_c	BU_LMPD_68_tyr-L_c	FU_LMPD_69_tyr-L_c	BU_LMPD_69_tyr-L_c	
FU_LMPD_70_ump_c	BU_LMPD_70_ump_c	FU_LMPD_71_ump_c	BU_LMPD_71_ump_c	FU_LMPD_72_val-L_c	BU_LMPD_72_val-L_c	
FU_LMPD_73_zymst_c	BU_LMPD_73_zymst_c	FU_LMPD_74_zymst_c	BU_LMPD_74_zymst_c	
End