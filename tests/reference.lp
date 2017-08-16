\Problem name: matlab_new_formulation_LP

Maximize
 obj: F_biomass_SC5_notrace
Subject To
M_13BDglcn_c: -1.1348 F_biomass_SC5_notrace +1.1348 R_biomass_SC5_notrace + F_LMPD_1_13BDglcn_c - R_LMPD_1_13BDglcn_c = 0
M_13dpg_c: F_GAPD - R_GAPD + F_PGK - R_PGK = 0
M_2pg_c: - F_ENO + R_ENO - F_PGM + R_PGM = 0
M_3pg_c: - F_PGK + R_PGK + F_PGM - R_PGM -2 F_LMPD_2_amp_c +2 R_LMPD_2_amp_c -2 F_LMPD_3_amp_c +2 R_LMPD_3_amp_c -2 F_LMPD_4_amp_c +2 R_LMPD_4_amp_c -2 F_LMPD_5_amp_c +2 R_LMPD_5_amp_c - F_LMPD_12_cys-L_c + R_LMPD_12_cys-L_c -2 F_LMPD_13_damp_c +2 R_LMPD_13_damp_c -2 F_LMPD_14_damp_c +2 R_LMPD_14_damp_c -2 F_LMPD_15_damp_c +2 R_LMPD_15_damp_c -2 F_LMPD_16_damp_c +2 R_LMPD_16_damp_c -2 F_LMPD_19_dgmp_c +2 R_LMPD_19_dgmp_c -2 F_LMPD_20_dgmp_c +2 R_LMPD_20_dgmp_c -2 F_LMPD_21_dgmp_c +2 R_LMPD_21_dgmp_c -2 F_LMPD_22_dgmp_c +2 R_LMPD_22_dgmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c -2 F_LMPD_30_gmp_c +2 R_LMPD_30_gmp_c -2 F_LMPD_31_gmp_c +2 R_LMPD_31_gmp_c - F_LMPD_32_his-L_c + R_LMPD_32_his-L_c - F_LMPD_33_his-L_c + R_LMPD_33_his-L_c - F_LMPD_34_his-L_c + R_LMPD_34_his-L_c - F_LMPD_35_his-L_c + R_LMPD_35_his-L_c - F_LMPD_44_met-L_c + R_LMPD_44_met-L_c - F_LMPD_45_met-L_c + R_LMPD_45_met-L_c -400 F_LMPD_48_pc_SC_c +400 R_LMPD_48_pc_SC_c -400 F_LMPD_49_pc_SC_c +400 R_LMPD_49_pc_SC_c -100 F_LMPD_51_pe_SC_c +100 R_LMPD_51_pe_SC_c -100 F_LMPD_56_ps_SC_c +100 R_LMPD_56_ps_SC_c - F_LMPD_61_ser-L_c + R_LMPD_61_ser-L_c - F_LMPD_67_trp-L_c + R_LMPD_67_trp-L_c = 0
M_4abut_c: F_ABUTt2r - R_ABUTt2r + F_GLUDC - R_GLUDC = 0
M_6pgc_c: - F_GND + R_GND + F_PGL - R_PGL = 0
M_6pgl_c: - F_6PGLter + R_6PGLter + F_G6PDH2 - R_G6PDH2 - F_PGL + R_PGL = 0
M_ac_c: - F_ACS + R_ACS + F_ACt2r - R_ACt2r - F_ACtm + R_ACtm + F_ACtr - R_ACtr + F_ALDD2y - R_ALDD2y + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c + F_LMPD_44_met-L_c - R_LMPD_44_met-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c = 0
M_acald_c: F_ACALDt - R_ACALDt + F_ACALDtm - R_ACALDtm + F_ALCD2if - R_ALCD2if - F_ALCD2ir + R_ALCD2ir + F_ALCD2x - R_ALCD2x - F_ALDD2y + R_ALDD2y + F_PYRDC - R_PYRDC = 0
M_accoa_c: F_ACS - R_ACS - F_LMPD_12_cys-L_c + R_LMPD_12_cys-L_c -18 F_LMPD_25_ergst_c +18 R_LMPD_25_ergst_c -18 F_LMPD_26_ergst_c +18 R_LMPD_26_ergst_c - F_LMPD_38_leu-L_c + R_LMPD_38_leu-L_c - F_LMPD_44_met-L_c + R_LMPD_44_met-L_c - F_LMPD_45_met-L_c + R_LMPD_45_met-L_c -1620 F_LMPD_47_pa_SC_c +1620 R_LMPD_47_pa_SC_c -1620 F_LMPD_48_pc_SC_c +1620 R_LMPD_48_pc_SC_c -1620 F_LMPD_49_pc_SC_c +1620 R_LMPD_49_pc_SC_c -1620 F_LMPD_51_pe_SC_c +1620 R_LMPD_51_pe_SC_c -1620 F_LMPD_56_ps_SC_c +1620 R_LMPD_56_ps_SC_c -1620 F_LMPD_57_ptd1ino_SC_c +1620 R_LMPD_57_ptd1ino_SC_c -1620 F_LMPD_58_ptd1ino_SC_c +1620 R_LMPD_58_ptd1ino_SC_c -2430 F_LMPD_65_triglyc_SC_c +2430 R_LMPD_65_triglyc_SC_c -2430 F_LMPD_66_triglyc_SC_c +2430 R_LMPD_66_triglyc_SC_c -18 F_LMPD_73_zymst_c +18 R_LMPD_73_zymst_c -18 F_LMPD_74_zymst_c +18 R_LMPD_74_zymst_c = 0
M_adp_c: +59.276 F_biomass_SC5_notrace -59.276 R_biomass_SC5_notrace +2 F_ADK1 -2 R_ADK1 + F_ATPM - R_ATPM + F_ATPS - R_ATPS - F_ATPtm-H + R_ATPtm-H + F_GLUK - R_GLUK + F_HEX1 - R_HEX1 + F_NADK - R_NADK - F_NDP1 + R_NDP1 + F_PC - R_PC + F_PFK - R_PFK + F_PFK_3 - R_PFK_3 + F_PGK - R_PGK + F_PPCK - R_PPCK - F_PYK + R_PYK + F_LMPD_1_13BDglcn_c - R_LMPD_1_13BDglcn_c +7 F_LMPD_2_amp_c -7 R_LMPD_2_amp_c +7 F_LMPD_3_amp_c -7 R_LMPD_3_amp_c +5 F_LMPD_4_amp_c -5 R_LMPD_4_amp_c +5 F_LMPD_5_amp_c -5 R_LMPD_5_amp_c +3 F_LMPD_6_arg-L_c -3 R_LMPD_6_arg-L_c + F_LMPD_7_asn-L_c - R_LMPD_7_asn-L_c +5 F_LMPD_8_cmp_c -5 R_LMPD_8_cmp_c +5 F_LMPD_9_cmp_c -5 R_LMPD_9_cmp_c +4 F_LMPD_10_cmp_c -4 R_LMPD_10_cmp_c +4 F_LMPD_11_cmp_c -4 R_LMPD_11_cmp_c + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c +5 F_LMPD_13_damp_c -5 R_LMPD_13_damp_c +3 F_LMPD_14_damp_c -3 R_LMPD_14_damp_c +5 F_LMPD_15_damp_c -5 R_LMPD_15_damp_c +3 F_LMPD_16_damp_c -3 R_LMPD_16_damp_c +3 F_LMPD_17_dcmp_c -3 R_LMPD_17_dcmp_c +3 F_LMPD_18_dcmp_c -3 R_LMPD_18_dcmp_c +8 F_LMPD_19_dgmp_c -8 R_LMPD_19_dgmp_c +7 F_LMPD_20_dgmp_c -7 R_LMPD_20_dgmp_c +7 F_LMPD_21_dgmp_c -7 R_LMPD_21_dgmp_c +8 F_LMPD_22_dgmp_c -8 R_LMPD_22_dgmp_c +3 F_LMPD_23_dtmp_c -3 R_LMPD_23_dtmp_c +3 F_LMPD_24_dtmp_c -3 R_LMPD_24_dtmp_c +20 F_LMPD_25_ergst_c -20 R_LMPD_25_ergst_c +20 F_LMPD_26_ergst_c -20 R_LMPD_26_ergst_c + F_LMPD_27_gln-L_c - R_LMPD_27_gln-L_c + F_LMPD_29_glycogen_c - R_LMPD_29_glycogen_c +7 F_LMPD_30_gmp_c -7 R_LMPD_30_gmp_c +7 F_LMPD_31_gmp_c -7 R_LMPD_31_gmp_c +2 F_LMPD_33_his-L_c -2 R_LMPD_33_his-L_c +2 F_LMPD_34_his-L_c -2 R_LMPD_34_his-L_c +2 F_LMPD_36_ile-L_c -2 R_LMPD_36_ile-L_c +2 F_LMPD_37_ile-L_c -2 R_LMPD_37_ile-L_c + F_LMPD_42_mannan_c - R_LMPD_42_mannan_c - F_LMPD_43_mannan_c + R_LMPD_43_mannan_c +2 F_LMPD_44_met-L_c -2 R_LMPD_44_met-L_c +2 F_LMPD_45_met-L_c -2 R_LMPD_45_met-L_c +1420 F_LMPD_47_pa_SC_c -1420 R_LMPD_47_pa_SC_c +1920 F_LMPD_48_pc_SC_c -1920 R_LMPD_48_pc_SC_c +1920 F_LMPD_49_pc_SC_c -1920 R_LMPD_49_pc_SC_c +1620 F_LMPD_51_pe_SC_c -1620 R_LMPD_51_pe_SC_c + F_LMPD_52_phe-L_c - R_LMPD_52_phe-L_c + F_LMPD_53_pro-L_c - R_LMPD_53_pro-L_c + F_LMPD_54_pro-L_c - R_LMPD_54_pro-L_c +1620 F_LMPD_56_ps_SC_c -1620 R_LMPD_56_ps_SC_c +1620 F_LMPD_57_ptd1ino_SC_c -1620 R_LMPD_57_ptd1ino_SC_c +1620 F_LMPD_58_ptd1ino_SC_c -1620 R_LMPD_58_ptd1ino_SC_c +9 F_LMPD_59_ribflv_c -9 R_LMPD_59_ribflv_c +11 F_LMPD_60_ribflv_c -11 R_LMPD_60_ribflv_c +2 F_LMPD_62_thr-L_c -2 R_LMPD_62_thr-L_c +2 F_LMPD_63_thr-L_c -2 R_LMPD_63_thr-L_c + F_LMPD_64_tre_c - R_LMPD_64_tre_c +2130 F_LMPD_65_triglyc_SC_c -2130 R_LMPD_65_triglyc_SC_c +2130 F_LMPD_66_triglyc_SC_c -2130 R_LMPD_66_triglyc_SC_c +2 F_LMPD_67_trp-L_c -2 R_LMPD_67_trp-L_c + F_LMPD_68_tyr-L_c - R_LMPD_68_tyr-L_c + F_LMPD_69_tyr-L_c - R_LMPD_69_tyr-L_c +3 F_LMPD_70_ump_c -3 R_LMPD_70_ump_c +3 F_LMPD_71_ump_c -3 R_LMPD_71_ump_c +18 F_LMPD_73_zymst_c -18 R_LMPD_73_zymst_c +18 F_LMPD_74_zymst_c -18 R_LMPD_74_zymst_c = 0
M_akg_c: - F_AKGMAL + R_AKGMAL + F_AKGt2r - R_AKGt2r - F_ALATA_L + R_ALATA_L - F_ASPTA + R_ASPTA + F_GLUDxi - R_GLUDxi - F_GLUDyi + R_GLUDyi + F_ICDHy - R_ICDHy +2 F_LMPD_2_amp_c -2 R_LMPD_2_amp_c +2 F_LMPD_3_amp_c -2 R_LMPD_3_amp_c +2 F_LMPD_4_amp_c -2 R_LMPD_4_amp_c +2 F_LMPD_5_amp_c -2 R_LMPD_5_amp_c + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c +2 F_LMPD_13_damp_c -2 R_LMPD_13_damp_c +2 F_LMPD_14_damp_c -2 R_LMPD_14_damp_c +2 F_LMPD_15_damp_c -2 R_LMPD_15_damp_c +2 F_LMPD_16_damp_c -2 R_LMPD_16_damp_c +2 F_LMPD_19_dgmp_c -2 R_LMPD_19_dgmp_c +2 F_LMPD_20_dgmp_c -2 R_LMPD_20_dgmp_c +2 F_LMPD_21_dgmp_c -2 R_LMPD_21_dgmp_c +2 F_LMPD_22_dgmp_c -2 R_LMPD_22_dgmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c +2 F_LMPD_30_gmp_c -2 R_LMPD_30_gmp_c +2 F_LMPD_31_gmp_c -2 R_LMPD_31_gmp_c +2 F_LMPD_32_his-L_c -2 R_LMPD_32_his-L_c +2 F_LMPD_33_his-L_c -2 R_LMPD_33_his-L_c +2 F_LMPD_34_his-L_c -2 R_LMPD_34_his-L_c +2 F_LMPD_35_his-L_c -2 R_LMPD_35_his-L_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_39_leu-L_c - R_LMPD_39_leu-L_c +2 F_LMPD_40_lys-L_c -2 R_LMPD_40_lys-L_c +2 F_LMPD_41_lys-L_c -2 R_LMPD_41_lys-L_c + F_LMPD_44_met-L_c - R_LMPD_44_met-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c +400 F_LMPD_48_pc_SC_c -400 R_LMPD_48_pc_SC_c +400 F_LMPD_49_pc_SC_c -400 R_LMPD_49_pc_SC_c +100 F_LMPD_51_pe_SC_c -100 R_LMPD_51_pe_SC_c + F_LMPD_52_phe-L_c - R_LMPD_52_phe-L_c +100 F_LMPD_56_ps_SC_c -100 R_LMPD_56_ps_SC_c + F_LMPD_61_ser-L_c - R_LMPD_61_ser-L_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c + F_LMPD_68_tyr-L_c - R_LMPD_68_tyr-L_c + F_LMPD_69_tyr-L_c - R_LMPD_69_tyr-L_c + F_LMPD_72_val-L_c - R_LMPD_72_val-L_c = 0
M_ala-L_c: -0.4588 F_biomass_SC5_notrace +0.4588 R_biomass_SC5_notrace - F_ALATA_L + R_ALATA_L + F_ALAt2r - R_ALAt2r + F_ALAtmi - R_ALAtmi - F_LMPD_28_gly_c + R_LMPD_28_gly_c - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c = 0
M_amp_c: -0.046 F_biomass_SC5_notrace +0.046 R_biomass_SC5_notrace + F_ACS - R_ACS - F_ADK1 + R_ADK1 + F_ATPH1 - R_ATPH1 + F_NDP1 - R_NDP1 +2 F_LMPD_2_amp_c -2 R_LMPD_2_amp_c +2 F_LMPD_3_amp_c -2 R_LMPD_3_amp_c +3 F_LMPD_4_amp_c -3 R_LMPD_4_amp_c +3 F_LMPD_5_amp_c -3 R_LMPD_5_amp_c + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c + F_LMPD_7_asn-L_c - R_LMPD_7_asn-L_c + F_LMPD_8_cmp_c - R_LMPD_8_cmp_c + F_LMPD_9_cmp_c - R_LMPD_9_cmp_c + F_LMPD_10_cmp_c - R_LMPD_10_cmp_c + F_LMPD_11_cmp_c - R_LMPD_11_cmp_c + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c +2 F_LMPD_13_damp_c -2 R_LMPD_13_damp_c +3 F_LMPD_14_damp_c -3 R_LMPD_14_damp_c +2 F_LMPD_15_damp_c -2 R_LMPD_15_damp_c +3 F_LMPD_16_damp_c -3 R_LMPD_16_damp_c + F_LMPD_17_dcmp_c - R_LMPD_17_dcmp_c + F_LMPD_18_dcmp_c - R_LMPD_18_dcmp_c +2 F_LMPD_19_dgmp_c -2 R_LMPD_19_dgmp_c +2 F_LMPD_20_dgmp_c -2 R_LMPD_20_dgmp_c +2 F_LMPD_21_dgmp_c -2 R_LMPD_21_dgmp_c +2 F_LMPD_22_dgmp_c -2 R_LMPD_22_dgmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c + F_LMPD_25_ergst_c - R_LMPD_25_ergst_c + F_LMPD_26_ergst_c - R_LMPD_26_ergst_c +2 F_LMPD_30_gmp_c -2 R_LMPD_30_gmp_c +2 F_LMPD_31_gmp_c -2 R_LMPD_31_gmp_c +3 F_LMPD_32_his-L_c -3 R_LMPD_32_his-L_c +2 F_LMPD_33_his-L_c -2 R_LMPD_33_his-L_c +2 F_LMPD_34_his-L_c -2 R_LMPD_34_his-L_c +3 F_LMPD_35_his-L_c -3 R_LMPD_35_his-L_c + F_LMPD_40_lys-L_c - R_LMPD_40_lys-L_c + F_LMPD_41_lys-L_c - R_LMPD_41_lys-L_c + F_LMPD_43_mannan_c - R_LMPD_43_mannan_c + F_LMPD_44_met-L_c - R_LMPD_44_met-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c +300 F_LMPD_48_pc_SC_c -300 R_LMPD_48_pc_SC_c +300 F_LMPD_49_pc_SC_c -300 R_LMPD_49_pc_SC_c +3 F_LMPD_59_ribflv_c -3 R_LMPD_59_ribflv_c +2 F_LMPD_60_ribflv_c -2 R_LMPD_60_ribflv_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c + F_LMPD_70_ump_c - R_LMPD_70_ump_c + F_LMPD_71_ump_c - R_LMPD_71_ump_c = 0
M_arg-L_c: -0.1607 F_biomass_SC5_notrace +0.1607 R_biomass_SC5_notrace + F_ARGt2r - R_ARGt2r + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c = 0
M_asn-L_c: -0.1017 F_biomass_SC5_notrace +0.1017 R_biomass_SC5_notrace + F_ASNt2r - R_ASNt2r + F_LMPD_7_asn-L_c - R_LMPD_7_asn-L_c = 0
M_asp-L_c: -0.2975 F_biomass_SC5_notrace +0.2975 R_biomass_SC5_notrace + F_ASPGLU2m - R_ASPGLU2m - F_ASPTA + R_ASPTA - F_ASPt2m + R_ASPt2m + F_ASPt2r - R_ASPt2r -2 F_LMPD_2_amp_c +2 R_LMPD_2_amp_c -2 F_LMPD_3_amp_c +2 R_LMPD_3_amp_c -2 F_LMPD_4_amp_c +2 R_LMPD_4_amp_c -2 F_LMPD_5_amp_c +2 R_LMPD_5_amp_c - F_LMPD_6_arg-L_c + R_LMPD_6_arg-L_c - F_LMPD_7_asn-L_c + R_LMPD_7_asn-L_c - F_LMPD_8_cmp_c + R_LMPD_8_cmp_c - F_LMPD_9_cmp_c + R_LMPD_9_cmp_c - F_LMPD_10_cmp_c + R_LMPD_10_cmp_c - F_LMPD_11_cmp_c + R_LMPD_11_cmp_c -2 F_LMPD_13_damp_c +2 R_LMPD_13_damp_c -2 F_LMPD_14_damp_c +2 R_LMPD_14_damp_c -2 F_LMPD_15_damp_c +2 R_LMPD_15_damp_c -2 F_LMPD_16_damp_c +2 R_LMPD_16_damp_c - F_LMPD_17_dcmp_c + R_LMPD_17_dcmp_c - F_LMPD_18_dcmp_c + R_LMPD_18_dcmp_c - F_LMPD_19_dgmp_c + R_LMPD_19_dgmp_c - F_LMPD_20_dgmp_c + R_LMPD_20_dgmp_c - F_LMPD_21_dgmp_c + R_LMPD_21_dgmp_c - F_LMPD_22_dgmp_c + R_LMPD_22_dgmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c - F_LMPD_30_gmp_c + R_LMPD_30_gmp_c - F_LMPD_31_gmp_c + R_LMPD_31_gmp_c - F_LMPD_32_his-L_c + R_LMPD_32_his-L_c - F_LMPD_33_his-L_c + R_LMPD_33_his-L_c - F_LMPD_34_his-L_c + R_LMPD_34_his-L_c - F_LMPD_35_his-L_c + R_LMPD_35_his-L_c - F_LMPD_36_ile-L_c + R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c - F_LMPD_44_met-L_c + R_LMPD_44_met-L_c - F_LMPD_45_met-L_c + R_LMPD_45_met-L_c - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c - F_LMPD_62_thr-L_c + R_LMPD_62_thr-L_c - F_LMPD_63_thr-L_c + R_LMPD_63_thr-L_c - F_LMPD_70_ump_c + R_LMPD_70_ump_c - F_LMPD_71_ump_c + R_LMPD_71_ump_c = 0
M_atp_c: -59.276 F_biomass_SC5_notrace +59.276 R_biomass_SC5_notrace - F_ACS + R_ACS - F_ADK1 + R_ADK1 - F_ATPH1 + R_ATPH1 - F_ATPM + R_ATPM - F_ATPS + R_ATPS + F_ATPtm-H - R_ATPtm-H - F_GLUK + R_GLUK - F_HEX1 + R_HEX1 - F_NADK + R_NADK - F_PC + R_PC - F_PFK + R_PFK - F_PFK_3 + R_PFK_3 - F_PGK + R_PGK - F_PPCK + R_PPCK + F_PYK - R_PYK - F_LMPD_1_13BDglcn_c + R_LMPD_1_13BDglcn_c -8 F_LMPD_2_amp_c +8 R_LMPD_2_amp_c -8 F_LMPD_3_amp_c +8 R_LMPD_3_amp_c -7 F_LMPD_4_amp_c +7 R_LMPD_4_amp_c -7 F_LMPD_5_amp_c +7 R_LMPD_5_amp_c -4 F_LMPD_6_arg-L_c +4 R_LMPD_6_arg-L_c -2 F_LMPD_7_asn-L_c +2 R_LMPD_7_asn-L_c -6 F_LMPD_8_cmp_c +6 R_LMPD_8_cmp_c -6 F_LMPD_9_cmp_c +6 R_LMPD_9_cmp_c -5 F_LMPD_10_cmp_c +5 R_LMPD_10_cmp_c -5 F_LMPD_11_cmp_c +5 R_LMPD_11_cmp_c -2 F_LMPD_12_cys-L_c +2 R_LMPD_12_cys-L_c -7 F_LMPD_13_damp_c +7 R_LMPD_13_damp_c -6 F_LMPD_14_damp_c +6 R_LMPD_14_damp_c -7 F_LMPD_15_damp_c +7 R_LMPD_15_damp_c -6 F_LMPD_16_damp_c +6 R_LMPD_16_damp_c -4 F_LMPD_17_dcmp_c +4 R_LMPD_17_dcmp_c -4 F_LMPD_18_dcmp_c +4 R_LMPD_18_dcmp_c -10 F_LMPD_19_dgmp_c +10 R_LMPD_19_dgmp_c -9 F_LMPD_20_dgmp_c +9 R_LMPD_20_dgmp_c -9 F_LMPD_21_dgmp_c +9 R_LMPD_21_dgmp_c -10 F_LMPD_22_dgmp_c +10 R_LMPD_22_dgmp_c -4 F_LMPD_23_dtmp_c +4 R_LMPD_23_dtmp_c -4 F_LMPD_24_dtmp_c +4 R_LMPD_24_dtmp_c -21 F_LMPD_25_ergst_c +21 R_LMPD_25_ergst_c -21 F_LMPD_26_ergst_c +21 R_LMPD_26_ergst_c - F_LMPD_27_gln-L_c + R_LMPD_27_gln-L_c - F_LMPD_29_glycogen_c + R_LMPD_29_glycogen_c -9 F_LMPD_30_gmp_c +9 R_LMPD_30_gmp_c -9 F_LMPD_31_gmp_c +9 R_LMPD_31_gmp_c -3 F_LMPD_32_his-L_c +3 R_LMPD_32_his-L_c -4 F_LMPD_33_his-L_c +4 R_LMPD_33_his-L_c -4 F_LMPD_34_his-L_c +4 R_LMPD_34_his-L_c -3 F_LMPD_35_his-L_c +3 R_LMPD_35_his-L_c -2 F_LMPD_36_ile-L_c +2 R_LMPD_36_ile-L_c -2 F_LMPD_37_ile-L_c +2 R_LMPD_37_ile-L_c - F_LMPD_40_lys-L_c + R_LMPD_40_lys-L_c - F_LMPD_41_lys-L_c + R_LMPD_41_lys-L_c - F_LMPD_42_mannan_c + R_LMPD_42_mannan_c -3 F_LMPD_44_met-L_c +3 R_LMPD_44_met-L_c -3 F_LMPD_45_met-L_c +3 R_LMPD_45_met-L_c -1420 F_LMPD_47_pa_SC_c +1420 R_LMPD_47_pa_SC_c -2220 F_LMPD_48_pc_SC_c +2220 R_LMPD_48_pc_SC_c -2220 F_LMPD_49_pc_SC_c +2220 R_LMPD_49_pc_SC_c -1620 F_LMPD_51_pe_SC_c +1620 R_LMPD_51_pe_SC_c - F_LMPD_52_phe-L_c + R_LMPD_52_phe-L_c - F_LMPD_53_pro-L_c + R_LMPD_53_pro-L_c - F_LMPD_54_pro-L_c + R_LMPD_54_pro-L_c -1620 F_LMPD_56_ps_SC_c +1620 R_LMPD_56_ps_SC_c -1620 F_LMPD_57_ptd1ino_SC_c +1620 R_LMPD_57_ptd1ino_SC_c -1620 F_LMPD_58_ptd1ino_SC_c +1620 R_LMPD_58_ptd1ino_SC_c -12 F_LMPD_59_ribflv_c +12 R_LMPD_59_ribflv_c -13 F_LMPD_60_ribflv_c +13 R_LMPD_60_ribflv_c -2 F_LMPD_62_thr-L_c +2 R_LMPD_62_thr-L_c -2 F_LMPD_63_thr-L_c +2 R_LMPD_63_thr-L_c - F_LMPD_64_tre_c + R_LMPD_64_tre_c -2130 F_LMPD_65_triglyc_SC_c +2130 R_LMPD_65_triglyc_SC_c -2130 F_LMPD_66_triglyc_SC_c +2130 R_LMPD_66_triglyc_SC_c -3 F_LMPD_67_trp-L_c +3 R_LMPD_67_trp-L_c - F_LMPD_68_tyr-L_c + R_LMPD_68_tyr-L_c - F_LMPD_69_tyr-L_c + R_LMPD_69_tyr-L_c -4 F_LMPD_70_ump_c +4 R_LMPD_70_ump_c -4 F_LMPD_71_ump_c +4 R_LMPD_71_ump_c -18 F_LMPD_73_zymst_c +18 R_LMPD_73_zymst_c -18 F_LMPD_74_zymst_c +18 R_LMPD_74_zymst_c = 0
M_cit_c: - F_ACONT + R_ACONT + F_CITt2r - R_CITt2r - F_CITtam + R_CITtam - F_CITtcm + R_CITtcm = 0
M_cmp_c: -0.0447 F_biomass_SC5_notrace +0.0447 R_biomass_SC5_notrace + F_LMPD_8_cmp_c - R_LMPD_8_cmp_c + F_LMPD_9_cmp_c - R_LMPD_9_cmp_c + F_LMPD_10_cmp_c - R_LMPD_10_cmp_c + F_LMPD_11_cmp_c - R_LMPD_11_cmp_c = 0
M_co2_c: F_CO2t - R_CO2t - F_CO2tm + R_CO2tm + F_FDH - R_FDH + F_FDNG - R_FDNG + F_GLUDC - R_GLUDC + F_GND - R_GND - F_HCO3E + R_HCO3E + F_ICDHy - R_ICDHy + F_PPCK - R_PPCK + F_PYRDC - R_PYRDC - F_LMPD_2_amp_c + R_LMPD_2_amp_c - F_LMPD_3_amp_c + R_LMPD_3_amp_c - F_LMPD_4_amp_c + R_LMPD_4_amp_c - F_LMPD_5_amp_c + R_LMPD_5_amp_c + F_LMPD_8_cmp_c - R_LMPD_8_cmp_c + F_LMPD_9_cmp_c - R_LMPD_9_cmp_c + F_LMPD_10_cmp_c - R_LMPD_10_cmp_c + F_LMPD_11_cmp_c - R_LMPD_11_cmp_c - F_LMPD_13_damp_c + R_LMPD_13_damp_c - F_LMPD_14_damp_c + R_LMPD_14_damp_c - F_LMPD_15_damp_c + R_LMPD_15_damp_c - F_LMPD_16_damp_c + R_LMPD_16_damp_c + F_LMPD_17_dcmp_c - R_LMPD_17_dcmp_c + F_LMPD_18_dcmp_c - R_LMPD_18_dcmp_c - F_LMPD_19_dgmp_c + R_LMPD_19_dgmp_c - F_LMPD_20_dgmp_c + R_LMPD_20_dgmp_c - F_LMPD_21_dgmp_c + R_LMPD_21_dgmp_c - F_LMPD_22_dgmp_c + R_LMPD_22_dgmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c +8 F_LMPD_25_ergst_c -8 R_LMPD_25_ergst_c +8 F_LMPD_26_ergst_c -8 R_LMPD_26_ergst_c - F_LMPD_30_gmp_c + R_LMPD_30_gmp_c - F_LMPD_31_gmp_c + R_LMPD_31_gmp_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_39_leu-L_c - R_LMPD_39_leu-L_c +1420 F_LMPD_47_pa_SC_c -1420 R_LMPD_47_pa_SC_c +1420 F_LMPD_48_pc_SC_c -1420 R_LMPD_48_pc_SC_c +1420 F_LMPD_49_pc_SC_c -1420 R_LMPD_49_pc_SC_c +1420 F_LMPD_51_pe_SC_c -1420 R_LMPD_51_pe_SC_c + F_LMPD_52_phe-L_c - R_LMPD_52_phe-L_c +1420 F_LMPD_56_ps_SC_c -1420 R_LMPD_56_ps_SC_c +1420 F_LMPD_57_ptd1ino_SC_c -1420 R_LMPD_57_ptd1ino_SC_c +1420 F_LMPD_58_ptd1ino_SC_c -1420 R_LMPD_58_ptd1ino_SC_c - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c +2130 F_LMPD_65_triglyc_SC_c -2130 R_LMPD_65_triglyc_SC_c +2130 F_LMPD_66_triglyc_SC_c -2130 R_LMPD_66_triglyc_SC_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c + F_LMPD_68_tyr-L_c - R_LMPD_68_tyr-L_c + F_LMPD_69_tyr-L_c - R_LMPD_69_tyr-L_c + F_LMPD_70_ump_c - R_LMPD_70_ump_c + F_LMPD_71_ump_c - R_LMPD_71_ump_c +8 F_LMPD_73_zymst_c -8 R_LMPD_73_zymst_c +8 F_LMPD_74_zymst_c -8 R_LMPD_74_zymst_c = 0
M_coa_c: - F_ACS + R_ACS + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c +18 F_LMPD_25_ergst_c -18 R_LMPD_25_ergst_c +18 F_LMPD_26_ergst_c -18 R_LMPD_26_ergst_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_44_met-L_c - R_LMPD_44_met-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c +1620 F_LMPD_47_pa_SC_c -1620 R_LMPD_47_pa_SC_c +1620 F_LMPD_48_pc_SC_c -1620 R_LMPD_48_pc_SC_c +1620 F_LMPD_49_pc_SC_c -1620 R_LMPD_49_pc_SC_c +1620 F_LMPD_51_pe_SC_c -1620 R_LMPD_51_pe_SC_c +1620 F_LMPD_56_ps_SC_c -1620 R_LMPD_56_ps_SC_c +1620 F_LMPD_57_ptd1ino_SC_c -1620 R_LMPD_57_ptd1ino_SC_c +1620 F_LMPD_58_ptd1ino_SC_c -1620 R_LMPD_58_ptd1ino_SC_c +2430 F_LMPD_65_triglyc_SC_c -2430 R_LMPD_65_triglyc_SC_c +2430 F_LMPD_66_triglyc_SC_c -2430 R_LMPD_66_triglyc_SC_c +18 F_LMPD_73_zymst_c -18 R_LMPD_73_zymst_c +18 F_LMPD_74_zymst_c -18 R_LMPD_74_zymst_c = 0
M_cys-L_c: -0.0066 F_biomass_SC5_notrace +0.0066 R_biomass_SC5_notrace + F_CYSt2r - R_CYSt2r + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c = 0
M_damp_c: -0.0036 F_biomass_SC5_notrace +0.0036 R_biomass_SC5_notrace + F_LMPD_13_damp_c - R_LMPD_13_damp_c + F_LMPD_14_damp_c - R_LMPD_14_damp_c + F_LMPD_15_damp_c - R_LMPD_15_damp_c + F_LMPD_16_damp_c - R_LMPD_16_damp_c = 0
M_db4p_c: F_DB4PS - R_DB4PS -2 F_LMPD_59_ribflv_c +2 R_LMPD_59_ribflv_c -2 F_LMPD_60_ribflv_c +2 R_LMPD_60_ribflv_c = 0
M_dcmp_c: -0.0024 F_biomass_SC5_notrace +0.0024 R_biomass_SC5_notrace + F_LMPD_17_dcmp_c - R_LMPD_17_dcmp_c + F_LMPD_18_dcmp_c - R_LMPD_18_dcmp_c = 0
M_dgmp_c: -0.0024 F_biomass_SC5_notrace +0.0024 R_biomass_SC5_notrace + F_LMPD_19_dgmp_c - R_LMPD_19_dgmp_c + F_LMPD_20_dgmp_c - R_LMPD_20_dgmp_c + F_LMPD_21_dgmp_c - R_LMPD_21_dgmp_c + F_LMPD_22_dgmp_c - R_LMPD_22_dgmp_c = 0
M_dhap_c: F_FBA - R_FBA + F_FBA3 - R_FBA3 - F_TPI + R_TPI -100 F_LMPD_47_pa_SC_c +100 R_LMPD_47_pa_SC_c -100 F_LMPD_48_pc_SC_c +100 R_LMPD_48_pc_SC_c -100 F_LMPD_49_pc_SC_c +100 R_LMPD_49_pc_SC_c -100 F_LMPD_51_pe_SC_c +100 R_LMPD_51_pe_SC_c -100 F_LMPD_56_ps_SC_c +100 R_LMPD_56_ps_SC_c -100 F_LMPD_57_ptd1ino_SC_c +100 R_LMPD_57_ptd1ino_SC_c -100 F_LMPD_58_ptd1ino_SC_c +100 R_LMPD_58_ptd1ino_SC_c -100 F_LMPD_65_triglyc_SC_c +100 R_LMPD_65_triglyc_SC_c -100 F_LMPD_66_triglyc_SC_c +100 R_LMPD_66_triglyc_SC_c = 0
M_dtmp_c: -0.0036 F_biomass_SC5_notrace +0.0036 R_biomass_SC5_notrace + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c = 0
M_e4p_c: F_FBA3 - R_FBA3 + F_TALA - R_TALA - F_TKT2 + R_TKT2 - F_LMPD_52_phe-L_c + R_LMPD_52_phe-L_c - F_LMPD_67_trp-L_c + R_LMPD_67_trp-L_c - F_LMPD_68_tyr-L_c + R_LMPD_68_tyr-L_c - F_LMPD_69_tyr-L_c + R_LMPD_69_tyr-L_c = 0
M_ergst_c: -0.0007 F_biomass_SC5_notrace +0.0007 R_biomass_SC5_notrace + F_ERGSTt - R_ERGSTt + F_ERGSTter - R_ERGSTter = 0
M_etoh_c: - F_ALCD2if + R_ALCD2if + F_ALCD2ir - R_ALCD2ir - F_ALCD2x + R_ALCD2x + F_ETOHt - R_ETOHt - F_ETOHtm + R_ETOHtm = 0
M_f6p_c: F_FBP - R_FBP + F_G6PI3 - R_G6PI3 + F_MAN6PI - R_MAN6PI - F_PFK + R_PFK + F_PGI - R_PGI + F_TALA - R_TALA + F_TKT2 - R_TKT2 = 0
M_fdp_c: - F_FBA + R_FBA - F_FBP + R_FBP + F_PFK - R_PFK = 0
M_for_c: F_DB4PS - R_DB4PS - F_FDH + R_FDH - F_FDNG + R_FDNG + F_FORt - R_FORt - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c + F_LMPD_73_zymst_c - R_LMPD_73_zymst_c + F_LMPD_74_zymst_c - R_LMPD_74_zymst_c = 0
M_fum_c: - F_FRDcm + R_FRDcm - F_FUM + R_FUM + F_FUMt2r - R_FUMt2r + F_SUCFUMtm - R_SUCFUMtm +2 F_LMPD_2_amp_c -2 R_LMPD_2_amp_c +2 F_LMPD_3_amp_c -2 R_LMPD_3_amp_c +2 F_LMPD_4_amp_c -2 R_LMPD_4_amp_c +2 F_LMPD_5_amp_c -2 R_LMPD_5_amp_c + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c - F_LMPD_9_cmp_c + R_LMPD_9_cmp_c - F_LMPD_11_cmp_c + R_LMPD_11_cmp_c +2 F_LMPD_13_damp_c -2 R_LMPD_13_damp_c +2 F_LMPD_14_damp_c -2 R_LMPD_14_damp_c +2 F_LMPD_15_damp_c -2 R_LMPD_15_damp_c +2 F_LMPD_16_damp_c -2 R_LMPD_16_damp_c - F_LMPD_17_dcmp_c + R_LMPD_17_dcmp_c + F_LMPD_19_dgmp_c - R_LMPD_19_dgmp_c + F_LMPD_20_dgmp_c - R_LMPD_20_dgmp_c + F_LMPD_21_dgmp_c - R_LMPD_21_dgmp_c + F_LMPD_22_dgmp_c - R_LMPD_22_dgmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c + F_LMPD_30_gmp_c - R_LMPD_30_gmp_c + F_LMPD_31_gmp_c - R_LMPD_31_gmp_c + F_LMPD_32_his-L_c - R_LMPD_32_his-L_c + F_LMPD_33_his-L_c - R_LMPD_33_his-L_c + F_LMPD_34_his-L_c - R_LMPD_34_his-L_c + F_LMPD_35_his-L_c - R_LMPD_35_his-L_c + F_LMPD_59_ribflv_c - R_LMPD_59_ribflv_c + F_LMPD_60_ribflv_c - R_LMPD_60_ribflv_c - F_LMPD_70_ump_c + R_LMPD_70_ump_c = 0
M_g1p_c: - F_PGMT + R_PGMT - F_LMPD_1_13BDglcn_c + R_LMPD_1_13BDglcn_c - F_LMPD_29_glycogen_c + R_LMPD_29_glycogen_c - F_LMPD_64_tre_c + R_LMPD_64_tre_c = 0
M_g3p_c: F_FBA - R_FBA - F_GAPD + R_GAPD - F_TALA + R_TALA + F_TKT1 - R_TKT1 + F_TKT2 - R_TKT2 + F_TPI - R_TPI + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c = 0
M_g6p_c: - F_G6PDH2 + R_G6PDH2 - F_G6PI + R_G6PI - F_G6Pter + R_G6Pter + F_HEX1 - R_HEX1 - F_PGI + R_PGI + F_PGMT - R_PGMT -100 F_LMPD_57_ptd1ino_SC_c +100 R_LMPD_57_ptd1ino_SC_c -100 F_LMPD_58_ptd1ino_SC_c +100 R_LMPD_58_ptd1ino_SC_c - F_LMPD_64_tre_c + R_LMPD_64_tre_c = 0
M_g6p-B_c: F_G6PI - R_G6PI - F_G6PI3 + R_G6PI3 + F_GLUK - R_GLUK = 0
M_glc-D_c: F_GLCt1 - R_GLCt1 - F_GLUK + R_GLUK - F_HEX1 + R_HEX1 = 0
M_gln-L_c: -0.1054 F_biomass_SC5_notrace +0.1054 R_biomass_SC5_notrace + F_GLNt2r - R_GLNt2r + F_LMPD_27_gln-L_c - R_LMPD_27_gln-L_c = 0
M_glu-L_c: -0.3018 F_biomass_SC5_notrace +0.3018 R_biomass_SC5_notrace + F_ALATA_L - R_ALATA_L - F_ASPGLU2m + R_ASPGLU2m + F_ASPTA - R_ASPTA - F_GLUDC + R_GLUDC - F_GLUDxi + R_GLUDxi + F_GLUDyi - R_GLUDyi + F_GLUt2r - R_GLUt2r - F_GLUt7m + R_GLUt7m -2 F_LMPD_2_amp_c +2 R_LMPD_2_amp_c -2 F_LMPD_3_amp_c +2 R_LMPD_3_amp_c -2 F_LMPD_4_amp_c +2 R_LMPD_4_amp_c -2 F_LMPD_5_amp_c +2 R_LMPD_5_amp_c - F_LMPD_12_cys-L_c + R_LMPD_12_cys-L_c -2 F_LMPD_13_damp_c +2 R_LMPD_13_damp_c -2 F_LMPD_14_damp_c +2 R_LMPD_14_damp_c -2 F_LMPD_15_damp_c +2 R_LMPD_15_damp_c -2 F_LMPD_16_damp_c +2 R_LMPD_16_damp_c -2 F_LMPD_19_dgmp_c +2 R_LMPD_19_dgmp_c -2 F_LMPD_20_dgmp_c +2 R_LMPD_20_dgmp_c -2 F_LMPD_21_dgmp_c +2 R_LMPD_21_dgmp_c -2 F_LMPD_22_dgmp_c +2 R_LMPD_22_dgmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c - F_LMPD_27_gln-L_c + R_LMPD_27_gln-L_c -2 F_LMPD_30_gmp_c +2 R_LMPD_30_gmp_c -2 F_LMPD_31_gmp_c +2 R_LMPD_31_gmp_c -2 F_LMPD_32_his-L_c +2 R_LMPD_32_his-L_c -2 F_LMPD_33_his-L_c +2 R_LMPD_33_his-L_c -2 F_LMPD_34_his-L_c +2 R_LMPD_34_his-L_c -2 F_LMPD_35_his-L_c +2 R_LMPD_35_his-L_c - F_LMPD_38_leu-L_c + R_LMPD_38_leu-L_c - F_LMPD_39_leu-L_c + R_LMPD_39_leu-L_c -2 F_LMPD_40_lys-L_c +2 R_LMPD_40_lys-L_c -2 F_LMPD_41_lys-L_c +2 R_LMPD_41_lys-L_c - F_LMPD_44_met-L_c + R_LMPD_44_met-L_c - F_LMPD_45_met-L_c + R_LMPD_45_met-L_c -400 F_LMPD_48_pc_SC_c +400 R_LMPD_48_pc_SC_c -400 F_LMPD_49_pc_SC_c +400 R_LMPD_49_pc_SC_c -100 F_LMPD_51_pe_SC_c +100 R_LMPD_51_pe_SC_c - F_LMPD_52_phe-L_c + R_LMPD_52_phe-L_c - F_LMPD_53_pro-L_c + R_LMPD_53_pro-L_c - F_LMPD_54_pro-L_c + R_LMPD_54_pro-L_c -100 F_LMPD_56_ps_SC_c +100 R_LMPD_56_ps_SC_c - F_LMPD_61_ser-L_c + R_LMPD_61_ser-L_c - F_LMPD_67_trp-L_c + R_LMPD_67_trp-L_c - F_LMPD_68_tyr-L_c + R_LMPD_68_tyr-L_c - F_LMPD_69_tyr-L_c + R_LMPD_69_tyr-L_c - F_LMPD_72_val-L_c + R_LMPD_72_val-L_c = 0
M_gly_c: -0.2904 F_biomass_SC5_notrace +0.2904 R_biomass_SC5_notrace + F_GLYt2r - R_GLYt2r + F_LMPD_2_amp_c - R_LMPD_2_amp_c + F_LMPD_3_amp_c - R_LMPD_3_amp_c + F_LMPD_4_amp_c - R_LMPD_4_amp_c + F_LMPD_5_amp_c - R_LMPD_5_amp_c + F_LMPD_13_damp_c - R_LMPD_13_damp_c + F_LMPD_14_damp_c - R_LMPD_14_damp_c + F_LMPD_15_damp_c - R_LMPD_15_damp_c + F_LMPD_16_damp_c - R_LMPD_16_damp_c + F_LMPD_19_dgmp_c - R_LMPD_19_dgmp_c + F_LMPD_20_dgmp_c - R_LMPD_20_dgmp_c + F_LMPD_21_dgmp_c - R_LMPD_21_dgmp_c + F_LMPD_22_dgmp_c - R_LMPD_22_dgmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c + F_LMPD_28_gly_c - R_LMPD_28_gly_c + F_LMPD_30_gmp_c - R_LMPD_30_gmp_c + F_LMPD_31_gmp_c - R_LMPD_31_gmp_c + F_LMPD_32_his-L_c - R_LMPD_32_his-L_c + F_LMPD_33_his-L_c - R_LMPD_33_his-L_c + F_LMPD_34_his-L_c - R_LMPD_34_his-L_c + F_LMPD_35_his-L_c - R_LMPD_35_his-L_c + F_LMPD_44_met-L_c - R_LMPD_44_met-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c +300 F_LMPD_48_pc_SC_c -300 R_LMPD_48_pc_SC_c +300 F_LMPD_49_pc_SC_c -300 R_LMPD_49_pc_SC_c = 0
M_glycogen_c: -0.5185 F_biomass_SC5_notrace +0.5185 R_biomass_SC5_notrace + F_LMPD_29_glycogen_c - R_LMPD_29_glycogen_c = 0
M_gmp_c: -0.046 F_biomass_SC5_notrace +0.046 R_biomass_SC5_notrace + F_LMPD_30_gmp_c - R_LMPD_30_gmp_c + F_LMPD_31_gmp_c - R_LMPD_31_gmp_c = 0
M_h_c: +58.7 F_biomass_SC5_notrace -58.7 R_biomass_SC5_notrace + F_ABUTt2r - R_ABUTt2r - F_ACS + R_ACS + F_ACt2r - R_ACt2r + F_AKGt2r - R_AKGt2r + F_ALAt2r - R_ALAt2r + F_ALCD2if - R_ALCD2if - F_ALCD2ir + R_ALCD2ir + F_ALCD2x - R_ALCD2x +2 F_ALDD2y -2 R_ALDD2y + F_ARGt2r - R_ARGt2r + F_ASNt2r - R_ASNt2r - F_ASPt2m + R_ASPt2m + F_ASPt2r - R_ASPt2r +2 F_ATPH1 -2 R_ATPH1 + F_ATPM - R_ATPM -3 F_ATPS3m +3 R_ATPS3m - F_ATPtm-H + R_ATPtm-H + F_CITt2r - R_CITt2r +2 F_CYOOm -2 R_CYOOm +3.5 F_CYOR_u6m -3.5 R_CYOR_u6m + F_CYSt2r - R_CYSt2r - F_D-LACtm + R_D-LACtm + F_DB4PS - R_DB4PS + F_FBP - R_FBP - F_FDNG + R_FDNG + F_FUMt2r - R_FUMt2r + F_G6PDH2 - R_G6PDH2 + F_GLNt2r - R_GLNt2r - F_GLUDC + R_GLUDC + F_GLUDxi - R_GLUDxi - F_GLUDyi + R_GLUDyi + F_GLUt2r - R_GLUt2r + F_GLYt2r - R_GLYt2r + F_HCO3E - R_HCO3E + F_HISt2r - R_HISt2r + F_ILEt2r - R_ILEt2r + F_LEUt2r - R_LEUt2r + F_LYSt2r - R_LYSt2r + F_MALt2r - R_MALt2r + F_MDH - R_MDH + F_METt2r - R_METt2r - F_NADH2-u6cm + R_NADH2-u6cm + F_NADPPPS - R_NADPPPS + F_NDP1 - R_NDP1 - F_OAAt2m + R_OAAt2m + F_PC - R_PC - F_PGK + R_PGK + F_PGL - R_PGL + F_PHEt2r - R_PHEt2r - F_PIt2m + R_PIt2m + F_PIt2r - R_PIt2r +2 F_PPA -2 R_PPA - F_PPCK + R_PPCK + F_PROt2r - R_PROt2r - F_PYRDC + R_PYRDC - F_PYRt2m + R_PYRt2m + F_SERt2r - R_SERt2r + F_SUCCt2r - R_SUCCt2r - F_THRt2m + R_THRt2m + F_THRt2r - R_THRt2r + F_TRPt2r - R_TRPt2r + F_TYRt2r - R_TYRt2r + F_VALt2r - R_VALt2r + F_LMPD_1_13BDglcn_c - R_LMPD_1_13BDglcn_c +17 F_LMPD_2_amp_c -17 R_LMPD_2_amp_c +17 F_LMPD_3_amp_c -17 R_LMPD_3_amp_c +17 F_LMPD_4_amp_c -17 R_LMPD_4_amp_c +17 F_LMPD_5_amp_c -17 R_LMPD_5_amp_c +3 F_LMPD_6_arg-L_c -3 R_LMPD_6_arg-L_c + F_LMPD_7_asn-L_c - R_LMPD_7_asn-L_c +5 F_LMPD_8_cmp_c -5 R_LMPD_8_cmp_c +5 F_LMPD_9_cmp_c -5 R_LMPD_9_cmp_c +4 F_LMPD_10_cmp_c -4 R_LMPD_10_cmp_c +4 F_LMPD_11_cmp_c -4 R_LMPD_11_cmp_c -2 F_LMPD_12_cys-L_c +2 R_LMPD_12_cys-L_c +16 F_LMPD_13_damp_c -16 R_LMPD_13_damp_c +16 F_LMPD_14_damp_c -16 R_LMPD_14_damp_c +16 F_LMPD_15_damp_c -16 R_LMPD_15_damp_c +16 F_LMPD_16_damp_c -16 R_LMPD_16_damp_c +2 F_LMPD_17_dcmp_c -2 R_LMPD_17_dcmp_c +2 F_LMPD_18_dcmp_c -2 R_LMPD_18_dcmp_c +18 F_LMPD_19_dgmp_c -18 R_LMPD_19_dgmp_c +17 F_LMPD_20_dgmp_c -17 R_LMPD_20_dgmp_c +17 F_LMPD_21_dgmp_c -17 R_LMPD_21_dgmp_c +18 F_LMPD_22_dgmp_c -18 R_LMPD_22_dgmp_c +2 F_LMPD_23_dtmp_c -2 R_LMPD_23_dtmp_c +2 F_LMPD_24_dtmp_c -2 R_LMPD_24_dtmp_c -17 F_LMPD_25_ergst_c +17 R_LMPD_25_ergst_c -17 F_LMPD_26_ergst_c +17 R_LMPD_26_ergst_c + F_LMPD_27_gln-L_c - R_LMPD_27_gln-L_c + F_LMPD_29_glycogen_c - R_LMPD_29_glycogen_c +18 F_LMPD_30_gmp_c -18 R_LMPD_30_gmp_c +18 F_LMPD_31_gmp_c -18 R_LMPD_31_gmp_c +11 F_LMPD_32_his-L_c -11 R_LMPD_32_his-L_c +11 F_LMPD_33_his-L_c -11 R_LMPD_33_his-L_c +11 F_LMPD_34_his-L_c -11 R_LMPD_34_his-L_c +11 F_LMPD_35_his-L_c -11 R_LMPD_35_his-L_c - F_LMPD_36_ile-L_c + R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c -2 F_LMPD_40_lys-L_c +2 R_LMPD_40_lys-L_c -2 F_LMPD_41_lys-L_c +2 R_LMPD_41_lys-L_c + F_LMPD_42_mannan_c - R_LMPD_42_mannan_c + F_LMPD_43_mannan_c - R_LMPD_43_mannan_c -5 F_LMPD_44_met-L_c +5 R_LMPD_44_met-L_c -5 F_LMPD_45_met-L_c +5 R_LMPD_45_met-L_c -3058 F_LMPD_47_pa_SC_c +3058 R_LMPD_47_pa_SC_c -2258 F_LMPD_48_pc_SC_c +2258 R_LMPD_48_pc_SC_c -2258 F_LMPD_49_pc_SC_c +2258 R_LMPD_49_pc_SC_c -2858 F_LMPD_51_pe_SC_c +2858 R_LMPD_51_pe_SC_c +2 F_LMPD_52_phe-L_c -2 R_LMPD_52_phe-L_c - F_LMPD_53_pro-L_c + R_LMPD_53_pro-L_c - F_LMPD_54_pro-L_c + R_LMPD_54_pro-L_c -2858 F_LMPD_56_ps_SC_c +2858 R_LMPD_56_ps_SC_c -2958 F_LMPD_57_ptd1ino_SC_c +2958 R_LMPD_57_ptd1ino_SC_c -2958 F_LMPD_58_ptd1ino_SC_c +2958 R_LMPD_58_ptd1ino_SC_c +14 F_LMPD_59_ribflv_c -14 R_LMPD_59_ribflv_c +14 F_LMPD_60_ribflv_c -14 R_LMPD_60_ribflv_c +2 F_LMPD_61_ser-L_c -2 R_LMPD_61_ser-L_c - F_LMPD_62_thr-L_c + R_LMPD_62_thr-L_c - F_LMPD_63_thr-L_c + R_LMPD_63_thr-L_c + F_LMPD_64_tre_c - R_LMPD_64_tre_c -4537 F_LMPD_65_triglyc_SC_c +4537 R_LMPD_65_triglyc_SC_c -4537 F_LMPD_66_triglyc_SC_c +4537 R_LMPD_66_triglyc_SC_c +6 F_LMPD_67_trp-L_c -6 R_LMPD_67_trp-L_c +3 F_LMPD_68_tyr-L_c -3 R_LMPD_68_tyr-L_c +3 F_LMPD_69_tyr-L_c -3 R_LMPD_69_tyr-L_c +2 F_LMPD_70_ump_c -2 R_LMPD_70_ump_c +2 F_LMPD_71_ump_c -2 R_LMPD_71_ump_c -16 F_LMPD_73_zymst_c +16 R_LMPD_73_zymst_c -16 F_LMPD_74_zymst_c +16 R_LMPD_74_zymst_c = 0
M_h2o_c: -59.276 F_biomass_SC5_notrace +59.276 R_biomass_SC5_notrace - F_ALDD2y + R_ALDD2y -2 F_ATPH1 +2 R_ATPH1 - F_ATPM + R_ATPM - F_ATPS + R_ATPS + F_ENO - R_ENO - F_FBP + R_FBP - F_FUM + R_FUM - F_GLUDxi + R_GLUDxi + F_GLUDyi - R_GLUDyi + F_H2Ot - R_H2Ot - F_H2Oter + R_H2Oter - F_H2Otm + R_H2Otm - F_HCO3E + R_HCO3E - F_NADPPPS + R_NADPPPS - F_NDP1 + R_NDP1 - F_PGL + R_PGL - F_PPA + R_PPA -3 F_LMPD_2_amp_c +3 R_LMPD_2_amp_c -3 F_LMPD_3_amp_c +3 R_LMPD_3_amp_c -3 F_LMPD_4_amp_c +3 R_LMPD_4_amp_c -3 F_LMPD_5_amp_c +3 R_LMPD_5_amp_c - F_LMPD_6_arg-L_c + R_LMPD_6_arg-L_c - F_LMPD_7_asn-L_c + R_LMPD_7_asn-L_c - F_LMPD_8_cmp_c + R_LMPD_8_cmp_c - F_LMPD_9_cmp_c + R_LMPD_9_cmp_c + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c -2 F_LMPD_13_damp_c +2 R_LMPD_13_damp_c -2 F_LMPD_14_damp_c +2 R_LMPD_14_damp_c -2 F_LMPD_15_damp_c +2 R_LMPD_15_damp_c -2 F_LMPD_16_damp_c +2 R_LMPD_16_damp_c +2 F_LMPD_17_dcmp_c -2 R_LMPD_17_dcmp_c +2 F_LMPD_18_dcmp_c -2 R_LMPD_18_dcmp_c -5 F_LMPD_19_dgmp_c +5 R_LMPD_19_dgmp_c -4 F_LMPD_20_dgmp_c +4 R_LMPD_20_dgmp_c -4 F_LMPD_21_dgmp_c +4 R_LMPD_21_dgmp_c -5 F_LMPD_22_dgmp_c +5 R_LMPD_22_dgmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c +7 F_LMPD_25_ergst_c -7 R_LMPD_25_ergst_c +7 F_LMPD_26_ergst_c -7 R_LMPD_26_ergst_c -5 F_LMPD_30_gmp_c +5 R_LMPD_30_gmp_c -5 F_LMPD_31_gmp_c +5 R_LMPD_31_gmp_c -3 F_LMPD_32_his-L_c +3 R_LMPD_32_his-L_c -3 F_LMPD_33_his-L_c +3 R_LMPD_33_his-L_c -3 F_LMPD_34_his-L_c +3 R_LMPD_34_his-L_c -3 F_LMPD_35_his-L_c +3 R_LMPD_35_his-L_c - F_LMPD_36_ile-L_c + R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c - F_LMPD_38_leu-L_c + R_LMPD_38_leu-L_c +2 F_LMPD_44_met-L_c -2 R_LMPD_44_met-L_c +2 F_LMPD_45_met-L_c -2 R_LMPD_45_met-L_c +1656 F_LMPD_47_pa_SC_c -1656 R_LMPD_47_pa_SC_c +956 F_LMPD_48_pc_SC_c -956 R_LMPD_48_pc_SC_c +956 F_LMPD_49_pc_SC_c -956 R_LMPD_49_pc_SC_c +1556 F_LMPD_51_pe_SC_c -1556 R_LMPD_51_pe_SC_c + F_LMPD_52_phe-L_c - R_LMPD_52_phe-L_c + F_LMPD_53_pro-L_c - R_LMPD_53_pro-L_c + F_LMPD_54_pro-L_c - R_LMPD_54_pro-L_c +1556 F_LMPD_56_ps_SC_c -1556 R_LMPD_56_ps_SC_c +1556 F_LMPD_57_ptd1ino_SC_c -1556 R_LMPD_57_ptd1ino_SC_c +1556 F_LMPD_58_ptd1ino_SC_c -1556 R_LMPD_58_ptd1ino_SC_c -4 F_LMPD_59_ribflv_c +4 R_LMPD_59_ribflv_c -4 F_LMPD_60_ribflv_c +4 R_LMPD_60_ribflv_c - F_LMPD_61_ser-L_c + R_LMPD_61_ser-L_c - F_LMPD_62_thr-L_c + R_LMPD_62_thr-L_c - F_LMPD_63_thr-L_c + R_LMPD_63_thr-L_c - F_LMPD_64_tre_c + R_LMPD_64_tre_c +2384 F_LMPD_65_triglyc_SC_c -2384 R_LMPD_65_triglyc_SC_c +2384 F_LMPD_66_triglyc_SC_c -2384 R_LMPD_66_triglyc_SC_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c +6 F_LMPD_73_zymst_c -6 R_LMPD_73_zymst_c +6 F_LMPD_74_zymst_c -6 R_LMPD_74_zymst_c = 0
M_hco3_c: F_HCO3E - R_HCO3E - F_PC + R_PC - F_LMPD_6_arg-L_c + R_LMPD_6_arg-L_c - F_LMPD_8_cmp_c + R_LMPD_8_cmp_c - F_LMPD_9_cmp_c + R_LMPD_9_cmp_c - F_LMPD_10_cmp_c + R_LMPD_10_cmp_c - F_LMPD_11_cmp_c + R_LMPD_11_cmp_c - F_LMPD_17_dcmp_c + R_LMPD_17_dcmp_c - F_LMPD_18_dcmp_c + R_LMPD_18_dcmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c -1420 F_LMPD_47_pa_SC_c +1420 R_LMPD_47_pa_SC_c -1420 F_LMPD_48_pc_SC_c +1420 R_LMPD_48_pc_SC_c -1420 F_LMPD_49_pc_SC_c +1420 R_LMPD_49_pc_SC_c -1420 F_LMPD_51_pe_SC_c +1420 R_LMPD_51_pe_SC_c -1420 F_LMPD_56_ps_SC_c +1420 R_LMPD_56_ps_SC_c -1420 F_LMPD_57_ptd1ino_SC_c +1420 R_LMPD_57_ptd1ino_SC_c -1420 F_LMPD_58_ptd1ino_SC_c +1420 R_LMPD_58_ptd1ino_SC_c -2130 F_LMPD_65_triglyc_SC_c +2130 R_LMPD_65_triglyc_SC_c -2130 F_LMPD_66_triglyc_SC_c +2130 R_LMPD_66_triglyc_SC_c - F_LMPD_70_ump_c + R_LMPD_70_ump_c - F_LMPD_71_ump_c + R_LMPD_71_ump_c = 0
M_his-L_c: -0.0663 F_biomass_SC5_notrace +0.0663 R_biomass_SC5_notrace + F_HISt2r - R_HISt2r + F_LMPD_32_his-L_c - R_LMPD_32_his-L_c + F_LMPD_33_his-L_c - R_LMPD_33_his-L_c + F_LMPD_34_his-L_c - R_LMPD_34_his-L_c + F_LMPD_35_his-L_c - R_LMPD_35_his-L_c = 0
M_icit_c: F_ACONT - R_ACONT + F_CITtcm - R_CITtcm - F_ICDHy + R_ICDHy - F_LMPD_28_gly_c + R_LMPD_28_gly_c - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c = 0
M_ile-L_c: -0.1927 F_biomass_SC5_notrace +0.1927 R_biomass_SC5_notrace + F_ILEt2r - R_ILEt2r + F_ILEtmi - R_ILEtmi = 0
M_lac-D_c: - F_D-LACt2m + R_D-LACt2m - F_D-LACtm + R_D-LACtm = 0
M_leu-L_c: -0.2964 F_biomass_SC5_notrace +0.2964 R_biomass_SC5_notrace + F_LEUt2r - R_LEUt2r + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_39_leu-L_c - R_LMPD_39_leu-L_c = 0
M_lys-L_c: -0.2862 F_biomass_SC5_notrace +0.2862 R_biomass_SC5_notrace + F_LYSt2r - R_LYSt2r + F_LMPD_40_lys-L_c - R_LMPD_40_lys-L_c + F_LMPD_41_lys-L_c - R_LMPD_41_lys-L_c = 0
M_mal-L_c: F_AKGMAL - R_AKGMAL + F_CITtam - R_CITtam + F_FUM - R_FUM + F_MALt2r - R_MALt2r - F_MALtm + R_MALtm - F_MDH + R_MDH = 0
M_man6p_c: - F_MAN6PI + R_MAN6PI - F_LMPD_42_mannan_c + R_LMPD_42_mannan_c - F_LMPD_43_mannan_c + R_LMPD_43_mannan_c = 0
M_mannan_c: -0.8079 F_biomass_SC5_notrace +0.8079 R_biomass_SC5_notrace - F_MANNANter + R_MANNANter = 0
M_met-L_c: -0.0507 F_biomass_SC5_notrace +0.0507 R_biomass_SC5_notrace + F_METt2r - R_METt2r + F_LMPD_44_met-L_c - R_LMPD_44_met-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c = 0
M_nad_c: - F_ALCD2if + R_ALCD2if + F_ALCD2ir - R_ALCD2ir - F_ALCD2x + R_ALCD2x - F_FDH + R_FDH - F_GAPD + R_GAPD - F_GLUDxi + R_GLUDxi - F_MDH + R_MDH + F_NADH2-u6cm - R_NADH2-u6cm - F_NADK + R_NADK + F_NADPPPS - R_NADPPPS -2 F_LMPD_2_amp_c +2 R_LMPD_2_amp_c -4 F_LMPD_3_amp_c +4 R_LMPD_3_amp_c -2 F_LMPD_4_amp_c +2 R_LMPD_4_amp_c -4 F_LMPD_5_amp_c +4 R_LMPD_5_amp_c - F_LMPD_12_cys-L_c + R_LMPD_12_cys-L_c -4 F_LMPD_13_damp_c +4 R_LMPD_13_damp_c -4 F_LMPD_14_damp_c +4 R_LMPD_14_damp_c -2 F_LMPD_15_damp_c +2 R_LMPD_15_damp_c -2 F_LMPD_16_damp_c +2 R_LMPD_16_damp_c -5 F_LMPD_19_dgmp_c +5 R_LMPD_19_dgmp_c -5 F_LMPD_20_dgmp_c +5 R_LMPD_20_dgmp_c -3 F_LMPD_21_dgmp_c +3 R_LMPD_21_dgmp_c -3 F_LMPD_22_dgmp_c +3 R_LMPD_22_dgmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c -2 F_LMPD_25_ergst_c +2 R_LMPD_25_ergst_c + F_LMPD_26_ergst_c - R_LMPD_26_ergst_c -5 F_LMPD_30_gmp_c +5 R_LMPD_30_gmp_c -3 F_LMPD_31_gmp_c +3 R_LMPD_31_gmp_c -3 F_LMPD_32_his-L_c +3 R_LMPD_32_his-L_c -4 F_LMPD_33_his-L_c +4 R_LMPD_33_his-L_c -3 F_LMPD_34_his-L_c +3 R_LMPD_34_his-L_c -4 F_LMPD_35_his-L_c +4 R_LMPD_35_his-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c - F_LMPD_38_leu-L_c + R_LMPD_38_leu-L_c - F_LMPD_39_leu-L_c + R_LMPD_39_leu-L_c - F_LMPD_40_lys-L_c + R_LMPD_40_lys-L_c - F_LMPD_45_met-L_c + R_LMPD_45_met-L_c +100 F_LMPD_47_pa_SC_c -100 R_LMPD_47_pa_SC_c -400 F_LMPD_48_pc_SC_c +400 R_LMPD_48_pc_SC_c -300 F_LMPD_49_pc_SC_c +300 R_LMPD_49_pc_SC_c + F_LMPD_53_pro-L_c - R_LMPD_53_pro-L_c +100 F_LMPD_58_ptd1ino_SC_c -100 R_LMPD_58_ptd1ino_SC_c - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c - F_LMPD_61_ser-L_c + R_LMPD_61_ser-L_c + F_LMPD_62_thr-L_c - R_LMPD_62_thr-L_c +100 F_LMPD_65_triglyc_SC_c -100 R_LMPD_65_triglyc_SC_c - F_LMPD_67_trp-L_c + R_LMPD_67_trp-L_c - F_LMPD_68_tyr-L_c + R_LMPD_68_tyr-L_c + F_LMPD_73_zymst_c - R_LMPD_73_zymst_c -2 F_LMPD_74_zymst_c +2 R_LMPD_74_zymst_c = 0
M_nadh_c: F_ALCD2if - R_ALCD2if - F_ALCD2ir + R_ALCD2ir + F_ALCD2x - R_ALCD2x + F_FDH - R_FDH + F_GAPD - R_GAPD + F_GLUDxi - R_GLUDxi + F_MDH - R_MDH - F_NADH2-u6cm + R_NADH2-u6cm +2 F_LMPD_2_amp_c -2 R_LMPD_2_amp_c +4 F_LMPD_3_amp_c -4 R_LMPD_3_amp_c +2 F_LMPD_4_amp_c -2 R_LMPD_4_amp_c +4 F_LMPD_5_amp_c -4 R_LMPD_5_amp_c + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c +4 F_LMPD_13_damp_c -4 R_LMPD_13_damp_c +4 F_LMPD_14_damp_c -4 R_LMPD_14_damp_c +2 F_LMPD_15_damp_c -2 R_LMPD_15_damp_c +2 F_LMPD_16_damp_c -2 R_LMPD_16_damp_c +5 F_LMPD_19_dgmp_c -5 R_LMPD_19_dgmp_c +5 F_LMPD_20_dgmp_c -5 R_LMPD_20_dgmp_c +3 F_LMPD_21_dgmp_c -3 R_LMPD_21_dgmp_c +3 F_LMPD_22_dgmp_c -3 R_LMPD_22_dgmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c +2 F_LMPD_25_ergst_c -2 R_LMPD_25_ergst_c - F_LMPD_26_ergst_c + R_LMPD_26_ergst_c +5 F_LMPD_30_gmp_c -5 R_LMPD_30_gmp_c +3 F_LMPD_31_gmp_c -3 R_LMPD_31_gmp_c +3 F_LMPD_32_his-L_c -3 R_LMPD_32_his-L_c +4 F_LMPD_33_his-L_c -4 R_LMPD_33_his-L_c +3 F_LMPD_34_his-L_c -3 R_LMPD_34_his-L_c +4 F_LMPD_35_his-L_c -4 R_LMPD_35_his-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_39_leu-L_c - R_LMPD_39_leu-L_c + F_LMPD_40_lys-L_c - R_LMPD_40_lys-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c -100 F_LMPD_47_pa_SC_c +100 R_LMPD_47_pa_SC_c +400 F_LMPD_48_pc_SC_c -400 R_LMPD_48_pc_SC_c +300 F_LMPD_49_pc_SC_c -300 R_LMPD_49_pc_SC_c - F_LMPD_53_pro-L_c + R_LMPD_53_pro-L_c -100 F_LMPD_58_ptd1ino_SC_c +100 R_LMPD_58_ptd1ino_SC_c + F_LMPD_59_ribflv_c - R_LMPD_59_ribflv_c + F_LMPD_60_ribflv_c - R_LMPD_60_ribflv_c + F_LMPD_61_ser-L_c - R_LMPD_61_ser-L_c - F_LMPD_62_thr-L_c + R_LMPD_62_thr-L_c -100 F_LMPD_65_triglyc_SC_c +100 R_LMPD_65_triglyc_SC_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c + F_LMPD_68_tyr-L_c - R_LMPD_68_tyr-L_c - F_LMPD_73_zymst_c + R_LMPD_73_zymst_c +2 F_LMPD_74_zymst_c -2 R_LMPD_74_zymst_c = 0
M_nadp_c: - F_ALDD2y + R_ALDD2y - F_G6PDH2 + R_G6PDH2 + F_GLUDyi - R_GLUDyi - F_GND + R_GND - F_ICDHy + R_ICDHy + F_NADK - R_NADK - F_NADPPPS + R_NADPPPS -2 F_LMPD_2_amp_c +2 R_LMPD_2_amp_c -2 F_LMPD_4_amp_c +2 R_LMPD_4_amp_c +4 F_LMPD_12_cys-L_c -4 R_LMPD_12_cys-L_c + F_LMPD_13_damp_c - R_LMPD_13_damp_c + F_LMPD_14_damp_c - R_LMPD_14_damp_c - F_LMPD_15_damp_c + R_LMPD_15_damp_c - F_LMPD_16_damp_c + R_LMPD_16_damp_c + F_LMPD_17_dcmp_c - R_LMPD_17_dcmp_c + F_LMPD_18_dcmp_c - R_LMPD_18_dcmp_c + F_LMPD_19_dgmp_c - R_LMPD_19_dgmp_c + F_LMPD_20_dgmp_c - R_LMPD_20_dgmp_c - F_LMPD_21_dgmp_c + R_LMPD_21_dgmp_c - F_LMPD_22_dgmp_c + R_LMPD_22_dgmp_c +2 F_LMPD_23_dtmp_c -2 R_LMPD_23_dtmp_c +2 F_LMPD_24_dtmp_c -2 R_LMPD_24_dtmp_c +27 F_LMPD_25_ergst_c -27 R_LMPD_25_ergst_c +24 F_LMPD_26_ergst_c -24 R_LMPD_26_ergst_c -2 F_LMPD_31_gmp_c +2 R_LMPD_31_gmp_c - F_LMPD_32_his-L_c + R_LMPD_32_his-L_c - F_LMPD_34_his-L_c + R_LMPD_34_his-L_c +2 F_LMPD_36_ile-L_c -2 R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c +2 F_LMPD_40_lys-L_c -2 R_LMPD_40_lys-L_c + F_LMPD_41_lys-L_c - R_LMPD_41_lys-L_c +6 F_LMPD_44_met-L_c -6 R_LMPD_44_met-L_c +7 F_LMPD_45_met-L_c -7 R_LMPD_45_met-L_c +2958 F_LMPD_47_pa_SC_c -2958 R_LMPD_47_pa_SC_c +3358 F_LMPD_48_pc_SC_c -3358 R_LMPD_48_pc_SC_c +3258 F_LMPD_49_pc_SC_c -3258 R_LMPD_49_pc_SC_c +2958 F_LMPD_51_pe_SC_c -2958 R_LMPD_51_pe_SC_c + F_LMPD_52_phe-L_c - R_LMPD_52_phe-L_c + F_LMPD_53_pro-L_c - R_LMPD_53_pro-L_c +2 F_LMPD_54_pro-L_c -2 R_LMPD_54_pro-L_c +2958 F_LMPD_56_ps_SC_c -2958 R_LMPD_56_ps_SC_c +3058 F_LMPD_57_ptd1ino_SC_c -3058 R_LMPD_57_ptd1ino_SC_c +2958 F_LMPD_58_ptd1ino_SC_c -2958 R_LMPD_58_ptd1ino_SC_c + F_LMPD_59_ribflv_c - R_LMPD_59_ribflv_c + F_LMPD_60_ribflv_c - R_LMPD_60_ribflv_c + F_LMPD_62_thr-L_c - R_LMPD_62_thr-L_c +2 F_LMPD_63_thr-L_c -2 R_LMPD_63_thr-L_c +4437 F_LMPD_65_triglyc_SC_c -4437 R_LMPD_65_triglyc_SC_c +4537 F_LMPD_66_triglyc_SC_c -4537 R_LMPD_66_triglyc_SC_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c + F_LMPD_68_tyr-L_c - R_LMPD_68_tyr-L_c +22 F_LMPD_73_zymst_c -22 R_LMPD_73_zymst_c +25 F_LMPD_74_zymst_c -25 R_LMPD_74_zymst_c = 0
M_nadph_c: F_ALDD2y - R_ALDD2y + F_G6PDH2 - R_G6PDH2 - F_GLUDyi + R_GLUDyi + F_GND - R_GND + F_ICDHy - R_ICDHy +2 F_LMPD_2_amp_c -2 R_LMPD_2_amp_c +2 F_LMPD_4_amp_c -2 R_LMPD_4_amp_c -4 F_LMPD_12_cys-L_c +4 R_LMPD_12_cys-L_c - F_LMPD_13_damp_c + R_LMPD_13_damp_c - F_LMPD_14_damp_c + R_LMPD_14_damp_c + F_LMPD_15_damp_c - R_LMPD_15_damp_c + F_LMPD_16_damp_c - R_LMPD_16_damp_c - F_LMPD_17_dcmp_c + R_LMPD_17_dcmp_c - F_LMPD_18_dcmp_c + R_LMPD_18_dcmp_c - F_LMPD_19_dgmp_c + R_LMPD_19_dgmp_c - F_LMPD_20_dgmp_c + R_LMPD_20_dgmp_c + F_LMPD_21_dgmp_c - R_LMPD_21_dgmp_c + F_LMPD_22_dgmp_c - R_LMPD_22_dgmp_c -2 F_LMPD_23_dtmp_c +2 R_LMPD_23_dtmp_c -2 F_LMPD_24_dtmp_c +2 R_LMPD_24_dtmp_c -27 F_LMPD_25_ergst_c +27 R_LMPD_25_ergst_c -24 F_LMPD_26_ergst_c +24 R_LMPD_26_ergst_c +2 F_LMPD_31_gmp_c -2 R_LMPD_31_gmp_c + F_LMPD_32_his-L_c - R_LMPD_32_his-L_c + F_LMPD_34_his-L_c - R_LMPD_34_his-L_c -2 F_LMPD_36_ile-L_c +2 R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c -2 F_LMPD_40_lys-L_c +2 R_LMPD_40_lys-L_c - F_LMPD_41_lys-L_c + R_LMPD_41_lys-L_c -6 F_LMPD_44_met-L_c +6 R_LMPD_44_met-L_c -7 F_LMPD_45_met-L_c +7 R_LMPD_45_met-L_c -2958 F_LMPD_47_pa_SC_c +2958 R_LMPD_47_pa_SC_c -3358 F_LMPD_48_pc_SC_c +3358 R_LMPD_48_pc_SC_c -3258 F_LMPD_49_pc_SC_c +3258 R_LMPD_49_pc_SC_c -2958 F_LMPD_51_pe_SC_c +2958 R_LMPD_51_pe_SC_c - F_LMPD_52_phe-L_c + R_LMPD_52_phe-L_c - F_LMPD_53_pro-L_c + R_LMPD_53_pro-L_c -2 F_LMPD_54_pro-L_c +2 R_LMPD_54_pro-L_c -2958 F_LMPD_56_ps_SC_c +2958 R_LMPD_56_ps_SC_c -3058 F_LMPD_57_ptd1ino_SC_c +3058 R_LMPD_57_ptd1ino_SC_c -2958 F_LMPD_58_ptd1ino_SC_c +2958 R_LMPD_58_ptd1ino_SC_c - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c - F_LMPD_62_thr-L_c + R_LMPD_62_thr-L_c -2 F_LMPD_63_thr-L_c +2 R_LMPD_63_thr-L_c -4437 F_LMPD_65_triglyc_SC_c +4437 R_LMPD_65_triglyc_SC_c -4537 F_LMPD_66_triglyc_SC_c +4537 R_LMPD_66_triglyc_SC_c - F_LMPD_67_trp-L_c + R_LMPD_67_trp-L_c - F_LMPD_68_tyr-L_c + R_LMPD_68_tyr-L_c -22 F_LMPD_73_zymst_c +22 R_LMPD_73_zymst_c -25 F_LMPD_74_zymst_c +25 R_LMPD_74_zymst_c = 0
M_nh4_c: F_GLUDxi - R_GLUDxi - F_GLUDyi + R_GLUDyi + F_NH4t - R_NH4t - F_NH4tm + R_NH4tm -2 F_LMPD_2_amp_c +2 R_LMPD_2_amp_c -2 F_LMPD_3_amp_c +2 R_LMPD_3_amp_c -2 F_LMPD_4_amp_c +2 R_LMPD_4_amp_c -2 F_LMPD_5_amp_c +2 R_LMPD_5_amp_c - F_LMPD_6_arg-L_c + R_LMPD_6_arg-L_c - F_LMPD_7_asn-L_c + R_LMPD_7_asn-L_c -2 F_LMPD_8_cmp_c +2 R_LMPD_8_cmp_c -2 F_LMPD_9_cmp_c +2 R_LMPD_9_cmp_c -2 F_LMPD_10_cmp_c +2 R_LMPD_10_cmp_c -2 F_LMPD_11_cmp_c +2 R_LMPD_11_cmp_c -2 F_LMPD_13_damp_c +2 R_LMPD_13_damp_c -2 F_LMPD_14_damp_c +2 R_LMPD_14_damp_c -2 F_LMPD_15_damp_c +2 R_LMPD_15_damp_c -2 F_LMPD_16_damp_c +2 R_LMPD_16_damp_c -2 F_LMPD_17_dcmp_c +2 R_LMPD_17_dcmp_c -2 F_LMPD_18_dcmp_c +2 R_LMPD_18_dcmp_c -3 F_LMPD_19_dgmp_c +3 R_LMPD_19_dgmp_c -3 F_LMPD_20_dgmp_c +3 R_LMPD_20_dgmp_c -3 F_LMPD_21_dgmp_c +3 R_LMPD_21_dgmp_c -3 F_LMPD_22_dgmp_c +3 R_LMPD_22_dgmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c - F_LMPD_27_gln-L_c + R_LMPD_27_gln-L_c -3 F_LMPD_30_gmp_c +3 R_LMPD_30_gmp_c -3 F_LMPD_31_gmp_c +3 R_LMPD_31_gmp_c - F_LMPD_32_his-L_c + R_LMPD_32_his-L_c - F_LMPD_33_his-L_c + R_LMPD_33_his-L_c - F_LMPD_34_his-L_c + R_LMPD_34_his-L_c - F_LMPD_35_his-L_c + R_LMPD_35_his-L_c -2 F_LMPD_59_ribflv_c +2 R_LMPD_59_ribflv_c -2 F_LMPD_60_ribflv_c +2 R_LMPD_60_ribflv_c - F_LMPD_67_trp-L_c + R_LMPD_67_trp-L_c - F_LMPD_70_ump_c + R_LMPD_70_ump_c - F_LMPD_71_ump_c + R_LMPD_71_ump_c = 0
M_o2_c: F_O2t - R_O2t - F_O2ter + R_O2ter - F_O2tm + R_O2tm -10 F_LMPD_25_ergst_c +10 R_LMPD_25_ergst_c -10 F_LMPD_26_ergst_c +10 R_LMPD_26_ergst_c -118 F_LMPD_47_pa_SC_c +118 R_LMPD_47_pa_SC_c -118 F_LMPD_48_pc_SC_c +118 R_LMPD_48_pc_SC_c -118 F_LMPD_49_pc_SC_c +118 R_LMPD_49_pc_SC_c -118 F_LMPD_51_pe_SC_c +118 R_LMPD_51_pe_SC_c -118 F_LMPD_56_ps_SC_c +118 R_LMPD_56_ps_SC_c -118 F_LMPD_57_ptd1ino_SC_c +118 R_LMPD_57_ptd1ino_SC_c -118 F_LMPD_58_ptd1ino_SC_c +118 R_LMPD_58_ptd1ino_SC_c -177 F_LMPD_65_triglyc_SC_c +177 R_LMPD_65_triglyc_SC_c -177 F_LMPD_66_triglyc_SC_c +177 R_LMPD_66_triglyc_SC_c -9 F_LMPD_73_zymst_c +9 R_LMPD_73_zymst_c -9 F_LMPD_74_zymst_c +9 R_LMPD_74_zymst_c = 0
M_oaa_c: F_ASPTA - R_ASPTA + F_MDH - R_MDH - F_OAAt + R_OAAt - F_OAAt2m + R_OAAt2m + F_PC - R_PC - F_PPCK + R_PPCK = 0
M_pa_SC_c: -6e-06 F_biomass_SC5_notrace +6e-06 R_biomass_SC5_notrace + F_LMPD_47_pa_SC_c - R_LMPD_47_pa_SC_c = 0
M_pc_SC_c: -6e-05 F_biomass_SC5_notrace +6e-05 R_biomass_SC5_notrace + F_LMPD_48_pc_SC_c - R_LMPD_48_pc_SC_c + F_LMPD_49_pc_SC_c - R_LMPD_49_pc_SC_c = 0
M_pe_SC_c: -4.5e-05 F_biomass_SC5_notrace +4.5e-05 R_biomass_SC5_notrace - F_PEtm_SC + R_PEtm_SC - F_LMPD_48_pc_SC_c + R_LMPD_48_pc_SC_c - F_LMPD_49_pc_SC_c + R_LMPD_49_pc_SC_c = 0
M_pep_c: F_ENO - R_ENO + F_PPCK - R_PPCK - F_PYK + R_PYK -2 F_LMPD_52_phe-L_c +2 R_LMPD_52_phe-L_c -2 F_LMPD_67_trp-L_c +2 R_LMPD_67_trp-L_c -2 F_LMPD_68_tyr-L_c +2 R_LMPD_68_tyr-L_c -2 F_LMPD_69_tyr-L_c +2 R_LMPD_69_tyr-L_c = 0
M_phe-L_c: -0.1339 F_biomass_SC5_notrace +0.1339 R_biomass_SC5_notrace + F_PHEt2r - R_PHEt2r + F_LMPD_52_phe-L_c - R_LMPD_52_phe-L_c = 0
M_pi_c: +59.305 F_biomass_SC5_notrace -59.305 R_biomass_SC5_notrace +2 F_ATPH1 -2 R_ATPH1 + F_ATPM - R_ATPM + F_ATPS - R_ATPS + F_FBP - R_FBP - F_GAPD + R_GAPD + F_MALtm - R_MALtm + F_NADPPPS - R_NADPPPS + F_NDP1 - R_NDP1 + F_PC - R_PC - F_PIt2m + R_PIt2m + F_PIt2r - R_PIt2r +2 F_PPA -2 R_PPA + F_SUCCtm - R_SUCCtm +9 F_LMPD_2_amp_c -9 R_LMPD_2_amp_c +9 F_LMPD_3_amp_c -9 R_LMPD_3_amp_c +9 F_LMPD_4_amp_c -9 R_LMPD_4_amp_c +9 F_LMPD_5_amp_c -9 R_LMPD_5_amp_c +3 F_LMPD_6_arg-L_c -3 R_LMPD_6_arg-L_c + F_LMPD_7_asn-L_c - R_LMPD_7_asn-L_c +5 F_LMPD_8_cmp_c -5 R_LMPD_8_cmp_c +5 F_LMPD_9_cmp_c -5 R_LMPD_9_cmp_c +4 F_LMPD_10_cmp_c -4 R_LMPD_10_cmp_c +4 F_LMPD_11_cmp_c -4 R_LMPD_11_cmp_c +2 F_LMPD_12_cys-L_c -2 R_LMPD_12_cys-L_c +9 F_LMPD_13_damp_c -9 R_LMPD_13_damp_c +9 F_LMPD_14_damp_c -9 R_LMPD_14_damp_c +9 F_LMPD_15_damp_c -9 R_LMPD_15_damp_c +9 F_LMPD_16_damp_c -9 R_LMPD_16_damp_c +3 F_LMPD_17_dcmp_c -3 R_LMPD_17_dcmp_c +3 F_LMPD_18_dcmp_c -3 R_LMPD_18_dcmp_c +10 F_LMPD_19_dgmp_c -10 R_LMPD_19_dgmp_c +9 F_LMPD_20_dgmp_c -9 R_LMPD_20_dgmp_c +9 F_LMPD_21_dgmp_c -9 R_LMPD_21_dgmp_c +10 F_LMPD_22_dgmp_c -10 R_LMPD_22_dgmp_c +4 F_LMPD_23_dtmp_c -4 R_LMPD_23_dtmp_c +4 F_LMPD_24_dtmp_c -4 R_LMPD_24_dtmp_c +8 F_LMPD_25_ergst_c -8 R_LMPD_25_ergst_c +8 F_LMPD_26_ergst_c -8 R_LMPD_26_ergst_c + F_LMPD_27_gln-L_c - R_LMPD_27_gln-L_c +9 F_LMPD_30_gmp_c -9 R_LMPD_30_gmp_c +9 F_LMPD_31_gmp_c -9 R_LMPD_31_gmp_c +4 F_LMPD_32_his-L_c -4 R_LMPD_32_his-L_c +4 F_LMPD_33_his-L_c -4 R_LMPD_33_his-L_c +4 F_LMPD_34_his-L_c -4 R_LMPD_34_his-L_c +4 F_LMPD_35_his-L_c -4 R_LMPD_35_his-L_c +2 F_LMPD_36_ile-L_c -2 R_LMPD_36_ile-L_c +2 F_LMPD_37_ile-L_c -2 R_LMPD_37_ile-L_c +3 F_LMPD_44_met-L_c -3 R_LMPD_44_met-L_c +3 F_LMPD_45_met-L_c -3 R_LMPD_45_met-L_c +1420 F_LMPD_47_pa_SC_c -1420 R_LMPD_47_pa_SC_c +2120 F_LMPD_48_pc_SC_c -2120 R_LMPD_48_pc_SC_c +2120 F_LMPD_49_pc_SC_c -2120 R_LMPD_49_pc_SC_c +1520 F_LMPD_51_pe_SC_c -1520 R_LMPD_51_pe_SC_c +4 F_LMPD_52_phe-L_c -4 R_LMPD_52_phe-L_c + F_LMPD_53_pro-L_c - R_LMPD_53_pro-L_c + F_LMPD_54_pro-L_c - R_LMPD_54_pro-L_c +1520 F_LMPD_56_ps_SC_c -1520 R_LMPD_56_ps_SC_c +1520 F_LMPD_57_ptd1ino_SC_c -1520 R_LMPD_57_ptd1ino_SC_c +1520 F_LMPD_58_ptd1ino_SC_c -1520 R_LMPD_58_ptd1ino_SC_c +12 F_LMPD_59_ribflv_c -12 R_LMPD_59_ribflv_c +12 F_LMPD_60_ribflv_c -12 R_LMPD_60_ribflv_c + F_LMPD_61_ser-L_c - R_LMPD_61_ser-L_c +2 F_LMPD_62_thr-L_c -2 R_LMPD_62_thr-L_c +2 F_LMPD_63_thr-L_c -2 R_LMPD_63_thr-L_c + F_LMPD_64_tre_c - R_LMPD_64_tre_c +2230 F_LMPD_65_triglyc_SC_c -2230 R_LMPD_65_triglyc_SC_c +2230 F_LMPD_66_triglyc_SC_c -2230 R_LMPD_66_triglyc_SC_c +6 F_LMPD_67_trp-L_c -6 R_LMPD_67_trp-L_c +4 F_LMPD_68_tyr-L_c -4 R_LMPD_68_tyr-L_c +4 F_LMPD_69_tyr-L_c -4 R_LMPD_69_tyr-L_c +3 F_LMPD_70_ump_c -3 R_LMPD_70_ump_c +3 F_LMPD_71_ump_c -3 R_LMPD_71_ump_c +6 F_LMPD_73_zymst_c -6 R_LMPD_73_zymst_c +6 F_LMPD_74_zymst_c -6 R_LMPD_74_zymst_c = 0
M_ppi_c: F_ACS - R_ACS - F_PPA + R_PPA + F_LMPD_1_13BDglcn_c - R_LMPD_1_13BDglcn_c + F_LMPD_2_amp_c - R_LMPD_2_amp_c + F_LMPD_3_amp_c - R_LMPD_3_amp_c + F_LMPD_4_amp_c - R_LMPD_4_amp_c + F_LMPD_5_amp_c - R_LMPD_5_amp_c + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c + F_LMPD_7_asn-L_c - R_LMPD_7_asn-L_c + F_LMPD_8_cmp_c - R_LMPD_8_cmp_c + F_LMPD_9_cmp_c - R_LMPD_9_cmp_c + F_LMPD_10_cmp_c - R_LMPD_10_cmp_c + F_LMPD_11_cmp_c - R_LMPD_11_cmp_c + F_LMPD_12_cys-L_c - R_LMPD_12_cys-L_c + F_LMPD_13_damp_c - R_LMPD_13_damp_c + F_LMPD_14_damp_c - R_LMPD_14_damp_c + F_LMPD_15_damp_c - R_LMPD_15_damp_c + F_LMPD_16_damp_c - R_LMPD_16_damp_c + F_LMPD_17_dcmp_c - R_LMPD_17_dcmp_c + F_LMPD_18_dcmp_c - R_LMPD_18_dcmp_c +2 F_LMPD_19_dgmp_c -2 R_LMPD_19_dgmp_c +2 F_LMPD_20_dgmp_c -2 R_LMPD_20_dgmp_c +2 F_LMPD_21_dgmp_c -2 R_LMPD_21_dgmp_c +2 F_LMPD_22_dgmp_c -2 R_LMPD_22_dgmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c +7 F_LMPD_25_ergst_c -7 R_LMPD_25_ergst_c +7 F_LMPD_26_ergst_c -7 R_LMPD_26_ergst_c + F_LMPD_29_glycogen_c - R_LMPD_29_glycogen_c +2 F_LMPD_30_gmp_c -2 R_LMPD_30_gmp_c +2 F_LMPD_31_gmp_c -2 R_LMPD_31_gmp_c +2 F_LMPD_32_his-L_c -2 R_LMPD_32_his-L_c +2 F_LMPD_33_his-L_c -2 R_LMPD_33_his-L_c +2 F_LMPD_34_his-L_c -2 R_LMPD_34_his-L_c +2 F_LMPD_35_his-L_c -2 R_LMPD_35_his-L_c + F_LMPD_40_lys-L_c - R_LMPD_40_lys-L_c + F_LMPD_41_lys-L_c - R_LMPD_41_lys-L_c + F_LMPD_42_mannan_c - R_LMPD_42_mannan_c + F_LMPD_43_mannan_c - R_LMPD_43_mannan_c + F_LMPD_44_met-L_c - R_LMPD_44_met-L_c + F_LMPD_45_met-L_c - R_LMPD_45_met-L_c +400 F_LMPD_48_pc_SC_c -400 R_LMPD_48_pc_SC_c +400 F_LMPD_49_pc_SC_c -400 R_LMPD_49_pc_SC_c +100 F_LMPD_51_pe_SC_c -100 R_LMPD_51_pe_SC_c +100 F_LMPD_56_ps_SC_c -100 R_LMPD_56_ps_SC_c +100 F_LMPD_57_ptd1ino_SC_c -100 R_LMPD_57_ptd1ino_SC_c +100 F_LMPD_58_ptd1ino_SC_c -100 R_LMPD_58_ptd1ino_SC_c +3 F_LMPD_59_ribflv_c -3 R_LMPD_59_ribflv_c +3 F_LMPD_60_ribflv_c -3 R_LMPD_60_ribflv_c + F_LMPD_64_tre_c - R_LMPD_64_tre_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c + F_LMPD_70_ump_c - R_LMPD_70_ump_c + F_LMPD_71_ump_c - R_LMPD_71_ump_c +6 F_LMPD_73_zymst_c -6 R_LMPD_73_zymst_c +6 F_LMPD_74_zymst_c -6 R_LMPD_74_zymst_c = 0
M_pro-L_c: -0.1647 F_biomass_SC5_notrace +0.1647 R_biomass_SC5_notrace + F_PROt2r - R_PROt2r + F_LMPD_53_pro-L_c - R_LMPD_53_pro-L_c + F_LMPD_54_pro-L_c - R_LMPD_54_pro-L_c = 0
M_ps_SC_c: -1.7e-05 F_biomass_SC5_notrace +1.7e-05 R_biomass_SC5_notrace - F_PStm_SC + R_PStm_SC + F_LMPD_48_pc_SC_c - R_LMPD_48_pc_SC_c + F_LMPD_49_pc_SC_c - R_LMPD_49_pc_SC_c + F_LMPD_51_pe_SC_c - R_LMPD_51_pe_SC_c + F_LMPD_56_ps_SC_c - R_LMPD_56_ps_SC_c = 0
M_ptd1ino_SC_c: -5.3e-05 F_biomass_SC5_notrace +5.3e-05 R_biomass_SC5_notrace + F_LMPD_57_ptd1ino_SC_c - R_LMPD_57_ptd1ino_SC_c + F_LMPD_58_ptd1ino_SC_c - R_LMPD_58_ptd1ino_SC_c = 0
M_pyr_c: F_ALATA_L - R_ALATA_L + F_D-LACt2m - R_D-LACt2m - F_PC + R_PC + F_PYK - R_PYK - F_PYRDC + R_PYRDC - F_PYRt + R_PYRt - F_PYRt2m + R_PYRt2m + F_LMPD_28_gly_c - R_LMPD_28_gly_c + F_LMPD_59_ribflv_c - R_LMPD_59_ribflv_c + F_LMPD_60_ribflv_c - R_LMPD_60_ribflv_c + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c = 0
M_r5p_c: - F_RPI + R_RPI - F_TKT1 + R_TKT1 - F_LMPD_2_amp_c + R_LMPD_2_amp_c - F_LMPD_3_amp_c + R_LMPD_3_amp_c - F_LMPD_4_amp_c + R_LMPD_4_amp_c - F_LMPD_5_amp_c + R_LMPD_5_amp_c - F_LMPD_8_cmp_c + R_LMPD_8_cmp_c - F_LMPD_9_cmp_c + R_LMPD_9_cmp_c - F_LMPD_10_cmp_c + R_LMPD_10_cmp_c - F_LMPD_11_cmp_c + R_LMPD_11_cmp_c - F_LMPD_13_damp_c + R_LMPD_13_damp_c - F_LMPD_14_damp_c + R_LMPD_14_damp_c - F_LMPD_15_damp_c + R_LMPD_15_damp_c - F_LMPD_16_damp_c + R_LMPD_16_damp_c - F_LMPD_17_dcmp_c + R_LMPD_17_dcmp_c - F_LMPD_18_dcmp_c + R_LMPD_18_dcmp_c - F_LMPD_19_dgmp_c + R_LMPD_19_dgmp_c - F_LMPD_20_dgmp_c + R_LMPD_20_dgmp_c - F_LMPD_21_dgmp_c + R_LMPD_21_dgmp_c - F_LMPD_22_dgmp_c + R_LMPD_22_dgmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_24_dtmp_c + R_LMPD_24_dtmp_c - F_LMPD_30_gmp_c + R_LMPD_30_gmp_c - F_LMPD_31_gmp_c + R_LMPD_31_gmp_c - F_LMPD_32_his-L_c + R_LMPD_32_his-L_c - F_LMPD_33_his-L_c + R_LMPD_33_his-L_c - F_LMPD_34_his-L_c + R_LMPD_34_his-L_c - F_LMPD_35_his-L_c + R_LMPD_35_his-L_c - F_LMPD_59_ribflv_c + R_LMPD_59_ribflv_c - F_LMPD_60_ribflv_c + R_LMPD_60_ribflv_c - F_LMPD_67_trp-L_c + R_LMPD_67_trp-L_c - F_LMPD_70_ump_c + R_LMPD_70_ump_c - F_LMPD_71_ump_c + R_LMPD_71_ump_c = 0
M_ribflv_c: -0.00099 F_biomass_SC5_notrace +0.00099 R_biomass_SC5_notrace + F_LMPD_59_ribflv_c - R_LMPD_59_ribflv_c + F_LMPD_60_ribflv_c - R_LMPD_60_ribflv_c = 0
M_ru5p-D_c: - F_DB4PS + R_DB4PS + F_GND - R_GND - F_RPE + R_RPE + F_RPI - R_RPI = 0
M_s17bp_c: - F_FBA3 + R_FBA3 + F_PFK_3 - R_PFK_3 = 0
M_s7p_c: - F_PFK_3 + R_PFK_3 - F_TALA + R_TALA + F_TKT1 - R_TKT1 = 0
M_ser-L_c: -0.1854 F_biomass_SC5_notrace +0.1854 R_biomass_SC5_notrace + F_SERt2r - R_SERt2r + F_LMPD_61_ser-L_c - R_LMPD_61_ser-L_c = 0
M_so4_c: -0.02 F_biomass_SC5_notrace +0.02 R_biomass_SC5_notrace + F_SO4ti - R_SO4ti - F_LMPD_12_cys-L_c + R_LMPD_12_cys-L_c - F_LMPD_44_met-L_c + R_LMPD_44_met-L_c - F_LMPD_45_met-L_c + R_LMPD_45_met-L_c = 0
M_succ_c: F_FRDcm - R_FRDcm + F_SUCCt2r - R_SUCCt2r - F_SUCCtm + R_SUCCtm - F_SUCFUMtm + R_SUCFUMtm + F_LMPD_9_cmp_c - R_LMPD_9_cmp_c + F_LMPD_11_cmp_c - R_LMPD_11_cmp_c + F_LMPD_17_dcmp_c - R_LMPD_17_dcmp_c + F_LMPD_24_dtmp_c - R_LMPD_24_dtmp_c + F_LMPD_28_gly_c - R_LMPD_28_gly_c + F_LMPD_59_ribflv_c - R_LMPD_59_ribflv_c + F_LMPD_60_ribflv_c - R_LMPD_60_ribflv_c + F_LMPD_70_ump_c - R_LMPD_70_ump_c = 0
M_thr-L_c: -0.1914 F_biomass_SC5_notrace +0.1914 R_biomass_SC5_notrace - F_THRt2m + R_THRt2m + F_THRt2r - R_THRt2r + F_LMPD_36_ile-L_c - R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c + F_LMPD_62_thr-L_c - R_LMPD_62_thr-L_c + F_LMPD_63_thr-L_c - R_LMPD_63_thr-L_c = 0
M_tre_c: -0.0234 F_biomass_SC5_notrace +0.0234 R_biomass_SC5_notrace + F_LMPD_64_tre_c - R_LMPD_64_tre_c = 0
M_triglyc_SC_c: -6.6e-05 F_biomass_SC5_notrace +6.6e-05 R_biomass_SC5_notrace + F_LMPD_65_triglyc_SC_c - R_LMPD_65_triglyc_SC_c + F_LMPD_66_triglyc_SC_c - R_LMPD_66_triglyc_SC_c = 0
M_trp-L_c: -0.0284 F_biomass_SC5_notrace +0.0284 R_biomass_SC5_notrace + F_TRPt2r - R_TRPt2r + F_LMPD_67_trp-L_c - R_LMPD_67_trp-L_c = 0
M_tyr-L_c: -0.102 F_biomass_SC5_notrace +0.102 R_biomass_SC5_notrace + F_TYRt2r - R_TYRt2r + F_LMPD_68_tyr-L_c - R_LMPD_68_tyr-L_c + F_LMPD_69_tyr-L_c - R_LMPD_69_tyr-L_c = 0
M_ump_c: -0.0599 F_biomass_SC5_notrace +0.0599 R_biomass_SC5_notrace + F_LMPD_70_ump_c - R_LMPD_70_ump_c + F_LMPD_71_ump_c - R_LMPD_71_ump_c = 0
M_val-L_c: -0.2646 F_biomass_SC5_notrace +0.2646 R_biomass_SC5_notrace + F_VALt2r - R_VALt2r + F_LMPD_72_val-L_c - R_LMPD_72_val-L_c = 0
M_xu5p-D_c: F_RPE - R_RPE - F_TKT1 + R_TKT1 - F_TKT2 + R_TKT2 = 0
M_zymst_c: -0.0015 F_biomass_SC5_notrace +0.0015 R_biomass_SC5_notrace + F_ZYMSTt - R_ZYMSTt + F_LMPD_73_zymst_c - R_LMPD_73_zymst_c + F_LMPD_74_zymst_c - R_LMPD_74_zymst_c = 0
M_4abut_e: - F_EX_4abut_e + R_EX_4abut_e - F_ABUTt2r + R_ABUTt2r = 0
M_ac_e: - F_EX_ac_e + R_EX_ac_e - F_ACt2r + R_ACt2r - F_ACtr + R_ACtr = 0
M_acald_e: - F_EX_acald_e + R_EX_acald_e - F_ACALDt + R_ACALDt = 0
M_akg_e: - F_EX_akg_e + R_EX_akg_e + F_AKGMAL - R_AKGMAL - F_AKGt2r + R_AKGt2r = 0
M_ala-L_e: - F_EX_ala-L_e + R_EX_ala-L_e - F_ALAt2r + R_ALAt2r = 0
M_arg-L_e: - F_EX_arg-L_e + R_EX_arg-L_e - F_ARGt2r + R_ARGt2r = 0
M_asn-L_e: - F_EX_asn-L_e + R_EX_asn-L_e - F_ASNt2r + R_ASNt2r = 0
M_asp-L_e: - F_EX_asp-L_e + R_EX_asp-L_e - F_ASPt2r + R_ASPt2r = 0
M_cit_e: - F_EX_cit_e + R_EX_cit_e - F_CITt2r + R_CITt2r = 0
M_co2_e: - F_EX_co2_e + R_EX_co2_e - F_CO2t + R_CO2t = 0
M_cys-L_e: - F_EX_cys-L_e + R_EX_cys-L_e - F_CYSt2r + R_CYSt2r = 0
M_ergst_e: - F_EX_ergst_e + R_EX_ergst_e - F_ERGSTt + R_ERGSTt = 0
M_etoh_e: - F_EX_etoh_e + R_EX_etoh_e - F_ETOHt + R_ETOHt = 0
M_for_e: - F_EX_for_e + R_EX_for_e - F_FORt + R_FORt = 0
M_fum_e: - F_EX_fum_e + R_EX_fum_e - F_FUMt2r + R_FUMt2r = 0
M_glc-D_e: - F_EX_glc_e + R_EX_glc_e - F_GLCt1 + R_GLCt1 = 0
M_gln-L_e: - F_EX_gln-L_e + R_EX_gln-L_e - F_GLNt2r + R_GLNt2r = 0
M_glu-L_e: - F_EX_glu-L_e + R_EX_glu-L_e - F_GLUt2r + R_GLUt2r = 0
M_gly_e: - F_EX_gly_e + R_EX_gly_e - F_GLYt2r + R_GLYt2r = 0
M_h_e: - F_EX_h_e + R_EX_h_e - F_ABUTt2r + R_ABUTt2r - F_ACt2r + R_ACt2r - F_AKGt2r + R_AKGt2r - F_ALAt2r + R_ALAt2r - F_ARGt2r + R_ARGt2r - F_ASNt2r + R_ASNt2r - F_ASPt2r + R_ASPt2r + F_ATPS - R_ATPS - F_CITt2r + R_CITt2r - F_CYSt2r + R_CYSt2r - F_FUMt2r + R_FUMt2r - F_GLNt2r + R_GLNt2r - F_GLUt2r + R_GLUt2r - F_GLYt2r + R_GLYt2r - F_HISt2r + R_HISt2r - F_ILEt2r + R_ILEt2r - F_LEUt2r + R_LEUt2r - F_LYSt2r + R_LYSt2r - F_MALt2r + R_MALt2r - F_METt2r + R_METt2r - F_PHEt2r + R_PHEt2r - F_PIt2r + R_PIt2r - F_PROt2r + R_PROt2r - F_SERt2r + R_SERt2r - F_SUCCt2r + R_SUCCt2r - F_THRt2r + R_THRt2r - F_TRPt2r + R_TRPt2r - F_TYRt2r + R_TYRt2r - F_VALt2r + R_VALt2r = 0
M_h2o_e: - F_EX_h2o_e + R_EX_h2o_e - F_H2Ot + R_H2Ot = 0
M_his-L_e: - F_EX_his-L_e + R_EX_his-L_e - F_HISt2r + R_HISt2r = 0
M_ile-L_e: - F_EX_ile-L_e + R_EX_ile-L_e - F_ILEt2r + R_ILEt2r = 0
M_leu-L_e: - F_EX_leu-L_e + R_EX_leu-L_e - F_LEUt2r + R_LEUt2r = 0
M_lys-L_e: - F_EX_lys-L_e + R_EX_lys-L_e - F_LYSt2r + R_LYSt2r = 0
M_mal-L_e: - F_EX_mal-L_e + R_EX_mal-L_e - F_AKGMAL + R_AKGMAL - F_MALt2r + R_MALt2r = 0
M_met-L_e: - F_EX_met-L_e + R_EX_met-L_e - F_METt2r + R_METt2r = 0
M_nh4_e: - F_EX_nh4_e + R_EX_nh4_e - F_NH4t + R_NH4t = 0
M_o2_e: - F_EX_o2_e + R_EX_o2_e - F_O2t + R_O2t = 0
M_oaa_e: - F_EX_oaa_e + R_EX_oaa_e + F_OAAt - R_OAAt = 0
M_phe-L_e: - F_EX_phe-L_e + R_EX_phe-L_e - F_PHEt2r + R_PHEt2r = 0
M_pi_e: - F_EX_pi_e + R_EX_pi_e - F_PIt2r + R_PIt2r = 0
M_pro-L_e: - F_EX_pro-L_e + R_EX_pro-L_e - F_PROt2r + R_PROt2r = 0
M_pyr_e: - F_EX_pyr_e + R_EX_pyr_e + F_PYRt - R_PYRt = 0
M_ser-L_e: - F_EX_ser-L_e + R_EX_ser-L_e - F_SERt2r + R_SERt2r = 0
M_so4_e: - F_EX_so4_e + R_EX_so4_e - F_SO4ti + R_SO4ti = 0
M_succ_e: - F_EX_succ_e + R_EX_succ_e - F_SUCCt2r + R_SUCCt2r = 0
M_thr-L_e: - F_EX_thr-L_e + R_EX_thr-L_e - F_THRt2r + R_THRt2r = 0
M_trp-L_e: - F_EX_trp-L_e + R_EX_trp-L_e - F_TRPt2r + R_TRPt2r = 0
M_tyr-L_e: - F_EX_tyr-L_e + R_EX_tyr-L_e - F_TYRt2r + R_TYRt2r = 0
M_val-L_e: - F_EX_val-L_e + R_EX_val-L_e - F_VALt2r + R_VALt2r = 0
M_zymst_e: - F_EX_zymst_e + R_EX_zymst_e - F_ZYMSTt + R_ZYMSTt = 0
M_ac_m: F_ACOAHim - R_ACOAHim - F_ACSm + R_ACSm + F_ACtm - R_ACtm + F_ALDD2xm - R_ALDD2xm + F_ALDD2ym - R_ALDD2ym = 0
M_acald_m: - F_ACALDtm + R_ACALDtm - F_ALCD2irm + R_ALCD2irm - F_ALDD2xm + R_ALDD2xm - F_ALDD2ym + R_ALDD2ym = 0
M_accoa_m: - F_ACOAHim + R_ACOAHim + F_ACSm - R_ACSm - F_CSm + R_CSm - F_HCITSm + R_HCITSm + F_PDHm - R_PDHm - F_LMPD_39_leu-L_c + R_LMPD_39_leu-L_c = 0
M_adp_m: +2 F_ADK1m -2 R_ADK1m - F_ATPS3m + R_ATPS3m + F_ATPtm-H - R_ATPtm-H + F_NADKm - R_NADKm + F_SUCOASm - R_SUCOASm + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c = 0
M_akg_m: - F_AKGDam + R_AKGDam - F_ALATA_Lm + R_ALATA_Lm - F_ASPTAm + R_ASPTAm - F_HCITSm + R_HCITSm + F_ICDHxm - R_ICDHxm + F_ICDHym - R_ICDHym + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c + F_LMPD_36_ile-L_c - R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c = 0
M_ala-L_m: - F_ALATA_Lm + R_ALATA_Lm - F_ALAtmi + R_ALAtmi = 0
M_amp_m: F_ACSm - R_ACSm - F_ADK1m + R_ADK1m = 0
M_asp-L_m: - F_ASPGLU2m + R_ASPGLU2m - F_ASPTAm + R_ASPTAm + F_ASPt2m - R_ASPt2m = 0
M_atp_m: - F_ACSm + R_ACSm - F_ADK1m + R_ADK1m + F_ATPS3m - R_ATPS3m - F_ATPtm-H + R_ATPtm-H - F_NADKm + R_NADKm - F_SUCOASm + R_SUCOASm - F_LMPD_6_arg-L_c + R_LMPD_6_arg-L_c = 0
M_cit_m: - F_ACONTm + R_ACONTm + F_CITtam - R_CITtam + F_CITtcm - R_CITtcm + F_CSm - R_CSm = 0
M_co2_m: F_AKGDam - R_AKGDam + F_CO2tm - R_CO2tm + F_ICDHxm - R_ICDHxm + F_ICDHym - R_ICDHym + F_ME1m - R_ME1m + F_ME2m - R_ME2m + F_PDHm - R_PDHm + F_LMPD_36_ile-L_c - R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_39_leu-L_c - R_LMPD_39_leu-L_c + F_LMPD_40_lys-L_c - R_LMPD_40_lys-L_c + F_LMPD_41_lys-L_c - R_LMPD_41_lys-L_c +100 F_LMPD_48_pc_SC_c -100 R_LMPD_48_pc_SC_c +100 F_LMPD_49_pc_SC_c -100 R_LMPD_49_pc_SC_c +100 F_LMPD_51_pe_SC_c -100 R_LMPD_51_pe_SC_c + F_LMPD_72_val-L_c - R_LMPD_72_val-L_c = 0
M_coa_m: F_ACOAHim - R_ACOAHim - F_ACSm + R_ACSm - F_AKGDbm + R_AKGDbm + F_CSm - R_CSm + F_HCITSm - R_HCITSm - F_PDHm + R_PDHm - F_SUCOASm + R_SUCOASm + F_LMPD_39_leu-L_c - R_LMPD_39_leu-L_c = 0
M_dhlam_m: F_AKGDbm - R_AKGDbm - F_GCC2cm + R_GCC2cm - F_PDHcm + R_PDHcm = 0
M_etoh_m: F_ALCD2irm - R_ALCD2irm + F_ETOHtm - R_ETOHtm = 0
M_fad_m: F_FRDcm - R_FRDcm + F_FRDm - R_FRDm - F_SUCD1m + R_SUCD1m + F_SUCD3-u6m - R_SUCD3-u6m = 0
M_fadh2_m: - F_FRDcm + R_FRDcm - F_FRDm + R_FRDm + F_SUCD1m - R_SUCD1m - F_SUCD3-u6m + R_SUCD3-u6m = 0
M_ficytc_m: +4 F_CYOOm -4 R_CYOOm -2 F_CYOR_u6m +2 R_CYOR_u6m = 0
M_focytc_m: -4 F_CYOOm +4 R_CYOOm +2 F_CYOR_u6m -2 R_CYOR_u6m = 0
M_fum_m: - F_FRDm + R_FRDm - F_FUMm + R_FUMm + F_SUCD1m - R_SUCD1m + F_SUCD2_u6m - R_SUCD2_u6m - F_SUCFUMtm + R_SUCFUMtm = 0
M_glu-L_m: F_ALATA_Lm - R_ALATA_Lm + F_ASPGLU2m - R_ASPGLU2m + F_ASPTAm - R_ASPTAm + F_GLUt7m - R_GLUt7m -2 F_LMPD_6_arg-L_c +2 R_LMPD_6_arg-L_c - F_LMPD_36_ile-L_c + R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c = 0
M_h_m: F_ACOAHim - R_ACOAHim - F_ACSm + R_ACSm - F_AKGDam + R_AKGDam - F_ALCD2irm + R_ALCD2irm +2 F_ALDD2xm -2 R_ALDD2xm +2 F_ALDD2ym -2 R_ALDD2ym + F_ASPt2m - R_ASPt2m +2 F_ATPS3m -2 R_ATPS3m + F_ATPtm-H - R_ATPtm-H + F_CSm - R_CSm -6 F_CYOOm +6 R_CYOOm -1.5 F_CYOR_u6m +1.5 R_CYOR_u6m + F_D-LACtm - R_D-LACtm + F_GCC2cm - R_GCC2cm + F_HCITSm - R_HCITSm + F_MDHm - R_MDHm - F_NADH2-u6m + R_NADH2-u6m + F_NADPPPSm - R_NADPPPSm + F_OAAt2m - R_OAAt2m + F_PDHcm - R_PDHcm + F_PIt2m - R_PIt2m +2 F_PPAm -2 R_PPAm + F_PYRt2m - R_PYRt2m + F_THRt2m - R_THRt2m -2 F_LMPD_36_ile-L_c +2 R_LMPD_36_ile-L_c -2 F_LMPD_37_ile-L_c +2 R_LMPD_37_ile-L_c -2 F_LMPD_38_leu-L_c +2 R_LMPD_38_leu-L_c - F_LMPD_39_leu-L_c + R_LMPD_39_leu-L_c -100 F_LMPD_48_pc_SC_c +100 R_LMPD_48_pc_SC_c -100 F_LMPD_49_pc_SC_c +100 R_LMPD_49_pc_SC_c -100 F_LMPD_51_pe_SC_c +100 R_LMPD_51_pe_SC_c -2 F_LMPD_72_val-L_c +2 R_LMPD_72_val-L_c = 0
M_h2o_m: - F_ACOAHim + R_ACOAHim - F_ALDD2xm + R_ALDD2xm - F_ALDD2ym + R_ALDD2ym + F_ATPS3m - R_ATPS3m - F_CSm + R_CSm +2 F_CYOOm -2 R_CYOOm - F_FUMm + R_FUMm + F_H2Otm - R_H2Otm - F_HCITSm + R_HCITSm - F_NADPPPSm + R_NADPPPSm - F_PPAm + R_PPAm + F_LMPD_36_ile-L_c - R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_72_val-L_c - R_LMPD_72_val-L_c = 0
M_hcit_m: F_HCITSm - R_HCITSm - F_LMPD_40_lys-L_c + R_LMPD_40_lys-L_c - F_LMPD_41_lys-L_c + R_LMPD_41_lys-L_c = 0
M_icit_m: F_ACONTm - R_ACONTm - F_CITtcm + R_CITtcm - F_ICDHxm + R_ICDHxm - F_ICDHym + R_ICDHym = 0
M_ile-L_m: - F_ILEtmi + R_ILEtmi + F_LMPD_36_ile-L_c - R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c = 0
M_lac-D_m: F_D-LACt2m - R_D-LACt2m + F_D-LACtm - R_D-LACtm = 0
M_lpam_m: - F_AKGDam + R_AKGDam + F_GCC2cm - R_GCC2cm + F_PDHcm - R_PDHcm = 0
M_mal-L_m: - F_CITtam + R_CITtam + F_FUMm - R_FUMm + F_MALtm - R_MALtm - F_MDHm + R_MDHm - F_ME1m + R_ME1m - F_ME2m + R_ME2m = 0
M_nad_m: F_ALCD2irm - R_ALCD2irm - F_ALDD2xm + R_ALDD2xm - F_GCC2cm + R_GCC2cm - F_ICDHxm + R_ICDHxm - F_MDHm + R_MDHm - F_ME1m + R_ME1m + F_NADH2-u6m - R_NADH2-u6m - F_NADKm + R_NADKm + F_NADPPPSm - R_NADPPPSm - F_PDHcm + R_PDHcm - F_PDHm + R_PDHm - F_LMPD_40_lys-L_c + R_LMPD_40_lys-L_c - F_LMPD_41_lys-L_c + R_LMPD_41_lys-L_c = 0
M_nadh_m: - F_ALCD2irm + R_ALCD2irm + F_ALDD2xm - R_ALDD2xm + F_GCC2cm - R_GCC2cm + F_ICDHxm - R_ICDHxm + F_MDHm - R_MDHm + F_ME1m - R_ME1m - F_NADH2-u6m + R_NADH2-u6m + F_PDHcm - R_PDHcm + F_PDHm - R_PDHm + F_LMPD_40_lys-L_c - R_LMPD_40_lys-L_c + F_LMPD_41_lys-L_c - R_LMPD_41_lys-L_c = 0
M_nadp_m: - F_ALDD2ym + R_ALDD2ym - F_ICDHym + R_ICDHym - F_ME2m + R_ME2m + F_NADKm - R_NADKm - F_NADPPPSm + R_NADPPPSm + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c + F_LMPD_36_ile-L_c - R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c + F_LMPD_38_leu-L_c - R_LMPD_38_leu-L_c + F_LMPD_39_leu-L_c - R_LMPD_39_leu-L_c + F_LMPD_72_val-L_c - R_LMPD_72_val-L_c = 0
M_nadph_m: F_ALDD2ym - R_ALDD2ym + F_ICDHym - R_ICDHym + F_ME2m - R_ME2m - F_LMPD_6_arg-L_c + R_LMPD_6_arg-L_c - F_LMPD_36_ile-L_c + R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c - F_LMPD_38_leu-L_c + R_LMPD_38_leu-L_c - F_LMPD_39_leu-L_c + R_LMPD_39_leu-L_c - F_LMPD_72_val-L_c + R_LMPD_72_val-L_c = 0
M_nh4_m: F_NH4tm - R_NH4tm + F_LMPD_36_ile-L_c - R_LMPD_36_ile-L_c + F_LMPD_37_ile-L_c - R_LMPD_37_ile-L_c = 0
M_o2_m: - F_CYOOm + R_CYOOm + F_O2tm - R_O2tm = 0
M_oaa_m: F_ASPTAm - R_ASPTAm - F_CSm + R_CSm + F_MDHm - R_MDHm + F_OAAt2m - R_OAAt2m = 0
M_pe_SC_m: F_PEtm_SC - R_PEtm_SC + F_LMPD_48_pc_SC_c - R_LMPD_48_pc_SC_c + F_LMPD_49_pc_SC_c - R_LMPD_49_pc_SC_c + F_LMPD_51_pe_SC_c - R_LMPD_51_pe_SC_c = 0
M_pi_m: - F_ATPS3m + R_ATPS3m - F_MALtm + R_MALtm + F_NADPPPSm - R_NADPPPSm + F_PIt2m - R_PIt2m +2 F_PPAm -2 R_PPAm - F_SUCCtm + R_SUCCtm + F_SUCOASm - R_SUCOASm + F_LMPD_6_arg-L_c - R_LMPD_6_arg-L_c = 0
M_ppi_m: F_ACSm - R_ACSm - F_PPAm + R_PPAm = 0
M_ps_SC_m: F_PStm_SC - R_PStm_SC - F_LMPD_48_pc_SC_c + R_LMPD_48_pc_SC_c - F_LMPD_49_pc_SC_c + R_LMPD_49_pc_SC_c - F_LMPD_51_pe_SC_c + R_LMPD_51_pe_SC_c = 0
M_pyr_m: F_ALATA_Lm - R_ALATA_Lm - F_D-LACt2m + R_D-LACt2m + F_ME1m - R_ME1m + F_ME2m - R_ME2m - F_PDHm + R_PDHm + F_PYRt2m - R_PYRt2m - F_LMPD_36_ile-L_c + R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c -2 F_LMPD_38_leu-L_c +2 R_LMPD_38_leu-L_c -2 F_LMPD_39_leu-L_c +2 R_LMPD_39_leu-L_c -2 F_LMPD_72_val-L_c +2 R_LMPD_72_val-L_c = 0
M_q6_m: F_CYOR_u6m - R_CYOR_u6m - F_FDNG + R_FDNG - F_NADH2-u6cm + R_NADH2-u6cm - F_NADH2-u6m + R_NADH2-u6m - F_SUCD2_u6m + R_SUCD2_u6m - F_SUCD3-u6m + R_SUCD3-u6m - F_LMPD_8_cmp_c + R_LMPD_8_cmp_c - F_LMPD_10_cmp_c + R_LMPD_10_cmp_c - F_LMPD_18_dcmp_c + R_LMPD_18_dcmp_c - F_LMPD_23_dtmp_c + R_LMPD_23_dtmp_c - F_LMPD_71_ump_c + R_LMPD_71_ump_c = 0
M_q6h2_m: - F_CYOR_u6m + R_CYOR_u6m + F_FDNG - R_FDNG + F_NADH2-u6cm - R_NADH2-u6cm + F_NADH2-u6m - R_NADH2-u6m + F_SUCD2_u6m - R_SUCD2_u6m + F_SUCD3-u6m - R_SUCD3-u6m + F_LMPD_8_cmp_c - R_LMPD_8_cmp_c + F_LMPD_10_cmp_c - R_LMPD_10_cmp_c + F_LMPD_18_dcmp_c - R_LMPD_18_dcmp_c + F_LMPD_23_dtmp_c - R_LMPD_23_dtmp_c + F_LMPD_71_ump_c - R_LMPD_71_ump_c = 0
M_sdhlam_m: F_AKGDam - R_AKGDam - F_AKGDbm + R_AKGDbm = 0
M_succ_m: F_FRDm - R_FRDm + F_SUCCtm - R_SUCCtm - F_SUCD1m + R_SUCD1m - F_SUCD2_u6m + R_SUCD2_u6m + F_SUCFUMtm - R_SUCFUMtm - F_SUCOASm + R_SUCOASm = 0
M_succoa_m: F_AKGDbm - R_AKGDbm + F_SUCOASm - R_SUCOASm = 0
M_thr-L_m: F_THRt2m - R_THRt2m - F_LMPD_36_ile-L_c + R_LMPD_36_ile-L_c - F_LMPD_37_ile-L_c + R_LMPD_37_ile-L_c = 0
M_6pgl_r: F_6PGLter - R_6PGLter + F_G6PDH2er - R_G6PDH2er = 0
M_ergst_r: - F_ERGSTter + R_ERGSTter + F_LMPD_25_ergst_c - R_LMPD_25_ergst_c + F_LMPD_26_ergst_c - R_LMPD_26_ergst_c = 0
M_g6p_r: - F_G6PDH2er + R_G6PDH2er + F_G6Pter - R_G6Pter = 0
M_h_r: F_G6PDH2er - R_G6PDH2er -2 F_LMPD_25_ergst_c +2 R_LMPD_25_ergst_c -2 F_LMPD_26_ergst_c +2 R_LMPD_26_ergst_c - F_LMPD_73_zymst_c + R_LMPD_73_zymst_c - F_LMPD_74_zymst_c + R_LMPD_74_zymst_c = 0
M_h2o_r: F_H2Oter - R_H2Oter + F_LMPD_25_ergst_c - R_LMPD_25_ergst_c + F_LMPD_26_ergst_c - R_LMPD_26_ergst_c + F_LMPD_73_zymst_c - R_LMPD_73_zymst_c + F_LMPD_74_zymst_c - R_LMPD_74_zymst_c = 0
M_mannan_r: F_MANNANter - R_MANNANter + F_LMPD_42_mannan_c - R_LMPD_42_mannan_c + F_LMPD_43_mannan_c - R_LMPD_43_mannan_c = 0
M_nadp_r: - F_G6PDH2er + R_G6PDH2er +2 F_LMPD_25_ergst_c -2 R_LMPD_25_ergst_c +2 F_LMPD_26_ergst_c -2 R_LMPD_26_ergst_c + F_LMPD_73_zymst_c - R_LMPD_73_zymst_c + F_LMPD_74_zymst_c - R_LMPD_74_zymst_c = 0
M_nadph_r: F_G6PDH2er - R_G6PDH2er -2 F_LMPD_25_ergst_c +2 R_LMPD_25_ergst_c -2 F_LMPD_26_ergst_c +2 R_LMPD_26_ergst_c - F_LMPD_73_zymst_c + R_LMPD_73_zymst_c - F_LMPD_74_zymst_c + R_LMPD_74_zymst_c = 0
M_o2_r: F_O2ter - R_O2ter - F_LMPD_25_ergst_c + R_LMPD_25_ergst_c - F_LMPD_26_ergst_c + R_LMPD_26_ergst_c - F_LMPD_73_zymst_c + R_LMPD_73_zymst_c - F_LMPD_74_zymst_c + R_LMPD_74_zymst_c = 0
SU_biomass_SC5_notrace: FU_biomass_SC5_notrace + BU_biomass_SC5_notrace < 1
UF_biomass_SC5_notrace: F_biomass_SC5_notrace -1000 FU_biomass_SC5_notrace < 0
UR_biomass_SC5_notrace: R_biomass_SC5_notrace -1000 BU_biomass_SC5_notrace < 0
SU_EX_4abut_e: FU_EX_4abut_e + BU_EX_4abut_e < 1
UF_EX_4abut_e: F_EX_4abut_e -1000 FU_EX_4abut_e < 0
UR_EX_4abut_e: R_EX_4abut_e -1000 BU_EX_4abut_e < 0
SU_EX_ac_e: FU_EX_ac_e + BU_EX_ac_e < 1
UF_EX_ac_e: F_EX_ac_e -1000 FU_EX_ac_e < 0
UR_EX_ac_e: R_EX_ac_e -1000 BU_EX_ac_e < 0
SU_EX_acald_e: FU_EX_acald_e + BU_EX_acald_e < 1
UF_EX_acald_e: F_EX_acald_e -1000 FU_EX_acald_e < 0
UR_EX_acald_e: R_EX_acald_e -1000 BU_EX_acald_e < 0
SU_EX_akg_e: FU_EX_akg_e + BU_EX_akg_e < 1
UF_EX_akg_e: F_EX_akg_e -1000 FU_EX_akg_e < 0
UR_EX_akg_e: R_EX_akg_e -1000 BU_EX_akg_e < 0
SU_EX_ala-L_e: FU_EX_ala-L_e + BU_EX_ala-L_e < 1
UF_EX_ala-L_e: F_EX_ala-L_e -1000 FU_EX_ala-L_e < 0
UR_EX_ala-L_e: R_EX_ala-L_e -1000 BU_EX_ala-L_e < 0
SU_EX_arg-L_e: FU_EX_arg-L_e + BU_EX_arg-L_e < 1
UF_EX_arg-L_e: F_EX_arg-L_e -1000 FU_EX_arg-L_e < 0
UR_EX_arg-L_e: R_EX_arg-L_e -1000 BU_EX_arg-L_e < 0
SU_EX_asn-L_e: FU_EX_asn-L_e + BU_EX_asn-L_e < 1
UF_EX_asn-L_e: F_EX_asn-L_e -1000 FU_EX_asn-L_e < 0
UR_EX_asn-L_e: R_EX_asn-L_e -1000 BU_EX_asn-L_e < 0
SU_EX_asp-L_e: FU_EX_asp-L_e + BU_EX_asp-L_e < 1
UF_EX_asp-L_e: F_EX_asp-L_e -1000 FU_EX_asp-L_e < 0
UR_EX_asp-L_e: R_EX_asp-L_e -1000 BU_EX_asp-L_e < 0
SU_EX_cit_e: FU_EX_cit_e + BU_EX_cit_e < 1
UF_EX_cit_e: F_EX_cit_e -1000 FU_EX_cit_e < 0
UR_EX_cit_e: R_EX_cit_e -1000 BU_EX_cit_e < 0
SU_EX_co2_e: FU_EX_co2_e + BU_EX_co2_e < 1
UF_EX_co2_e: F_EX_co2_e -1000 FU_EX_co2_e < 0
UR_EX_co2_e: R_EX_co2_e -1000 BU_EX_co2_e < 0
SU_EX_cys-L_e: FU_EX_cys-L_e + BU_EX_cys-L_e < 1
UF_EX_cys-L_e: F_EX_cys-L_e -1000 FU_EX_cys-L_e < 0
UR_EX_cys-L_e: R_EX_cys-L_e -1000 BU_EX_cys-L_e < 0
SU_EX_ergst_e: FU_EX_ergst_e + BU_EX_ergst_e < 1
UF_EX_ergst_e: F_EX_ergst_e -1000 FU_EX_ergst_e < 0
UR_EX_ergst_e: R_EX_ergst_e -1000 BU_EX_ergst_e < 0
SU_EX_etoh_e: FU_EX_etoh_e + BU_EX_etoh_e < 1
UF_EX_etoh_e: F_EX_etoh_e -1000 FU_EX_etoh_e < 0
UR_EX_etoh_e: R_EX_etoh_e -1000 BU_EX_etoh_e < 0
SU_EX_for_e: FU_EX_for_e + BU_EX_for_e < 1
UF_EX_for_e: F_EX_for_e -1000 FU_EX_for_e < 0
UR_EX_for_e: R_EX_for_e -1000 BU_EX_for_e < 0
SU_EX_fum_e: FU_EX_fum_e + BU_EX_fum_e < 1
UF_EX_fum_e: F_EX_fum_e -1000 FU_EX_fum_e < 0
UR_EX_fum_e: R_EX_fum_e -1000 BU_EX_fum_e < 0
SU_EX_glc_e: FU_EX_glc_e + BU_EX_glc_e < 1
UF_EX_glc_e: F_EX_glc_e -1000 FU_EX_glc_e < 0
UR_EX_glc_e: R_EX_glc_e -1000 BU_EX_glc_e < 0
SU_EX_gln-L_e: FU_EX_gln-L_e + BU_EX_gln-L_e < 1
UF_EX_gln-L_e: F_EX_gln-L_e -1000 FU_EX_gln-L_e < 0
UR_EX_gln-L_e: R_EX_gln-L_e -1000 BU_EX_gln-L_e < 0
SU_EX_glu-L_e: FU_EX_glu-L_e + BU_EX_glu-L_e < 1
UF_EX_glu-L_e: F_EX_glu-L_e -1000 FU_EX_glu-L_e < 0
UR_EX_glu-L_e: R_EX_glu-L_e -1000 BU_EX_glu-L_e < 0
SU_EX_gly_e: FU_EX_gly_e + BU_EX_gly_e < 1
UF_EX_gly_e: F_EX_gly_e -1000 FU_EX_gly_e < 0
UR_EX_gly_e: R_EX_gly_e -1000 BU_EX_gly_e < 0
SU_EX_h_e: FU_EX_h_e + BU_EX_h_e < 1
UF_EX_h_e: F_EX_h_e -1000 FU_EX_h_e < 0
UR_EX_h_e: R_EX_h_e -1000 BU_EX_h_e < 0
SU_EX_h2o_e: FU_EX_h2o_e + BU_EX_h2o_e < 1
UF_EX_h2o_e: F_EX_h2o_e -1000 FU_EX_h2o_e < 0
UR_EX_h2o_e: R_EX_h2o_e -1000 BU_EX_h2o_e < 0
SU_EX_his-L_e: FU_EX_his-L_e + BU_EX_his-L_e < 1
UF_EX_his-L_e: F_EX_his-L_e -1000 FU_EX_his-L_e < 0
UR_EX_his-L_e: R_EX_his-L_e -1000 BU_EX_his-L_e < 0
SU_EX_ile-L_e: FU_EX_ile-L_e + BU_EX_ile-L_e < 1
UF_EX_ile-L_e: F_EX_ile-L_e -1000 FU_EX_ile-L_e < 0
UR_EX_ile-L_e: R_EX_ile-L_e -1000 BU_EX_ile-L_e < 0
SU_EX_leu-L_e: FU_EX_leu-L_e + BU_EX_leu-L_e < 1
UF_EX_leu-L_e: F_EX_leu-L_e -1000 FU_EX_leu-L_e < 0
UR_EX_leu-L_e: R_EX_leu-L_e -1000 BU_EX_leu-L_e < 0
SU_EX_lys-L_e: FU_EX_lys-L_e + BU_EX_lys-L_e < 1
UF_EX_lys-L_e: F_EX_lys-L_e -1000 FU_EX_lys-L_e < 0
UR_EX_lys-L_e: R_EX_lys-L_e -1000 BU_EX_lys-L_e < 0
SU_EX_mal-L_e: FU_EX_mal-L_e + BU_EX_mal-L_e < 1
UF_EX_mal-L_e: F_EX_mal-L_e -1000 FU_EX_mal-L_e < 0
UR_EX_mal-L_e: R_EX_mal-L_e -1000 BU_EX_mal-L_e < 0
SU_EX_met-L_e: FU_EX_met-L_e + BU_EX_met-L_e < 1
UF_EX_met-L_e: F_EX_met-L_e -1000 FU_EX_met-L_e < 0
UR_EX_met-L_e: R_EX_met-L_e -1000 BU_EX_met-L_e < 0
SU_EX_nh4_e: FU_EX_nh4_e + BU_EX_nh4_e < 1
UF_EX_nh4_e: F_EX_nh4_e -1000 FU_EX_nh4_e < 0
UR_EX_nh4_e: R_EX_nh4_e -1000 BU_EX_nh4_e < 0
SU_EX_o2_e: FU_EX_o2_e + BU_EX_o2_e < 1
UF_EX_o2_e: F_EX_o2_e -1000 FU_EX_o2_e < 0
UR_EX_o2_e: R_EX_o2_e -1000 BU_EX_o2_e < 0
SU_EX_oaa_e: FU_EX_oaa_e + BU_EX_oaa_e < 1
UF_EX_oaa_e: F_EX_oaa_e -1000 FU_EX_oaa_e < 0
UR_EX_oaa_e: R_EX_oaa_e -1000 BU_EX_oaa_e < 0
SU_EX_phe-L_e: FU_EX_phe-L_e + BU_EX_phe-L_e < 1
UF_EX_phe-L_e: F_EX_phe-L_e -1000 FU_EX_phe-L_e < 0
UR_EX_phe-L_e: R_EX_phe-L_e -1000 BU_EX_phe-L_e < 0
SU_EX_pi_e: FU_EX_pi_e + BU_EX_pi_e < 1
UF_EX_pi_e: F_EX_pi_e -1000 FU_EX_pi_e < 0
UR_EX_pi_e: R_EX_pi_e -1000 BU_EX_pi_e < 0
SU_EX_pro-L_e: FU_EX_pro-L_e + BU_EX_pro-L_e < 1
UF_EX_pro-L_e: F_EX_pro-L_e -1000 FU_EX_pro-L_e < 0
UR_EX_pro-L_e: R_EX_pro-L_e -1000 BU_EX_pro-L_e < 0
SU_EX_pyr_e: FU_EX_pyr_e + BU_EX_pyr_e < 1
UF_EX_pyr_e: F_EX_pyr_e -1000 FU_EX_pyr_e < 0
UR_EX_pyr_e: R_EX_pyr_e -1000 BU_EX_pyr_e < 0
SU_EX_ser-L_e: FU_EX_ser-L_e + BU_EX_ser-L_e < 1
UF_EX_ser-L_e: F_EX_ser-L_e -1000 FU_EX_ser-L_e < 0
UR_EX_ser-L_e: R_EX_ser-L_e -1000 BU_EX_ser-L_e < 0
SU_EX_so4_e: FU_EX_so4_e + BU_EX_so4_e < 1
UF_EX_so4_e: F_EX_so4_e -1000 FU_EX_so4_e < 0
UR_EX_so4_e: R_EX_so4_e -1000 BU_EX_so4_e < 0
SU_EX_succ_e: FU_EX_succ_e + BU_EX_succ_e < 1
UF_EX_succ_e: F_EX_succ_e -1000 FU_EX_succ_e < 0
UR_EX_succ_e: R_EX_succ_e -1000 BU_EX_succ_e < 0
SU_EX_thr-L_e: FU_EX_thr-L_e + BU_EX_thr-L_e < 1
UF_EX_thr-L_e: F_EX_thr-L_e -1000 FU_EX_thr-L_e < 0
UR_EX_thr-L_e: R_EX_thr-L_e -1000 BU_EX_thr-L_e < 0
SU_EX_trp-L_e: FU_EX_trp-L_e + BU_EX_trp-L_e < 1
UF_EX_trp-L_e: F_EX_trp-L_e -1000 FU_EX_trp-L_e < 0
UR_EX_trp-L_e: R_EX_trp-L_e -1000 BU_EX_trp-L_e < 0
SU_EX_tyr-L_e: FU_EX_tyr-L_e + BU_EX_tyr-L_e < 1
UF_EX_tyr-L_e: F_EX_tyr-L_e -1000 FU_EX_tyr-L_e < 0
UR_EX_tyr-L_e: R_EX_tyr-L_e -1000 BU_EX_tyr-L_e < 0
SU_EX_val-L_e: FU_EX_val-L_e + BU_EX_val-L_e < 1
UF_EX_val-L_e: F_EX_val-L_e -1000 FU_EX_val-L_e < 0
UR_EX_val-L_e: R_EX_val-L_e -1000 BU_EX_val-L_e < 0
SU_EX_zymst_e: FU_EX_zymst_e + BU_EX_zymst_e < 1
UF_EX_zymst_e: F_EX_zymst_e -1000 FU_EX_zymst_e < 0
UR_EX_zymst_e: R_EX_zymst_e -1000 BU_EX_zymst_e < 0
G_6PGLter: -0.59209 LC_6pgl_c +0.59209 LC_6pgl_r - DG_6PGLter + DGo_6PGLter = 0
FU_6PGLter: DG_6PGLter +1000 FU_6PGLter < 1000
BU_6PGLter: - DG_6PGLter +1000 BU_6PGLter < 1000
SU_6PGLter: FU_6PGLter + BU_6PGLter < 1
UF_6PGLter: F_6PGLter -1000 FU_6PGLter < 0
UR_6PGLter: R_6PGLter -1000 BU_6PGLter < 0
G_ABUTt2r: +0.59209 LC_4abut_c -0.59209 LC_4abut_e - DG_ABUTt2r + DGo_ABUTt2r = 0
FU_ABUTt2r: DG_ABUTt2r +1000 FU_ABUTt2r < 1000
BU_ABUTt2r: - DG_ABUTt2r +1000 BU_ABUTt2r < 1000
SU_ABUTt2r: FU_ABUTt2r + BU_ABUTt2r < 1
UF_ABUTt2r: F_ABUTt2r -1000 FU_ABUTt2r < 0
UR_ABUTt2r: R_ABUTt2r -1000 BU_ABUTt2r < 0
G_ACALDt: +0.59209 LC_acald_c -0.59209 LC_acald_e - DG_ACALDt + DGo_ACALDt = 0
FU_ACALDt: DG_ACALDt +1000 FU_ACALDt < 1000
BU_ACALDt: - DG_ACALDt +1000 BU_ACALDt < 1000
SU_ACALDt: FU_ACALDt + BU_ACALDt < 1
UF_ACALDt: F_ACALDt -1000 FU_ACALDt < 0
UR_ACALDt: R_ACALDt -1000 BU_ACALDt < 0
G_ACALDtm: +0.59209 LC_acald_c -0.59209 LC_acald_m - DG_ACALDtm + DGo_ACALDtm = 0
FU_ACALDtm: DG_ACALDtm +1000 FU_ACALDtm < 1000
BU_ACALDtm: - DG_ACALDtm +1000 BU_ACALDtm < 1000
SU_ACALDtm: FU_ACALDtm + BU_ACALDtm < 1
UF_ACALDtm: F_ACALDtm -1000 FU_ACALDtm < 0
UR_ACALDtm: R_ACALDtm -1000 BU_ACALDtm < 0
G_ACOAHim: +0.59209 LC_ac_m -0.59209 LC_accoa_m +0.59209 LC_coa_m - DG_ACOAHim + DGo_ACOAHim = 0
FU_ACOAHim: DG_ACOAHim +1000 FU_ACOAHim < 1000
BU_ACOAHim: - DG_ACOAHim +1000 BU_ACOAHim < 1000
SU_ACOAHim: FU_ACOAHim + BU_ACOAHim < 1
UF_ACOAHim: F_ACOAHim -1000 FU_ACOAHim < 0
UR_ACOAHim: R_ACOAHim -1000 BU_ACOAHim < 0
G_ACONT: -0.59209 LC_cit_c +0.59209 LC_icit_c - DG_ACONT + DGo_ACONT = 0
FU_ACONT: DG_ACONT +1000 FU_ACONT < 1000
BU_ACONT: - DG_ACONT +1000 BU_ACONT < 1000
SU_ACONT: FU_ACONT + BU_ACONT < 1
UF_ACONT: F_ACONT -1000 FU_ACONT < 0
UR_ACONT: R_ACONT -1000 BU_ACONT < 0
G_ACONTm: -0.59209 LC_cit_m +0.59209 LC_icit_m - DG_ACONTm + DGo_ACONTm = 0
FU_ACONTm: DG_ACONTm +1000 FU_ACONTm < 1000
BU_ACONTm: - DG_ACONTm +1000 BU_ACONTm < 1000
SU_ACONTm: FU_ACONTm + BU_ACONTm < 1
UF_ACONTm: F_ACONTm -1000 FU_ACONTm < 0
UR_ACONTm: R_ACONTm -1000 BU_ACONTm < 0
G_ACS: -0.59209 LC_ac_c +0.59209 LC_accoa_c +0.59209 LC_amp_c -0.59209 LC_atp_c -0.59209 LC_coa_c +0.59209 LC_ppi_c - DG_ACS + DGo_ACS = 0
FU_ACS: DG_ACS +1000 FU_ACS < 1000
BU_ACS: - DG_ACS +1000 BU_ACS < 1000
SU_ACS: FU_ACS + BU_ACS < 1
UF_ACS: F_ACS -1000 FU_ACS < 0
UR_ACS: R_ACS -1000 BU_ACS < 0
G_ACSm: -0.59209 LC_ac_m +0.59209 LC_accoa_m +0.59209 LC_amp_m -0.59209 LC_atp_m -0.59209 LC_coa_m +0.59209 LC_ppi_m - DG_ACSm + DGo_ACSm = 0
FU_ACSm: DG_ACSm +1000 FU_ACSm < 1000
BU_ACSm: - DG_ACSm +1000 BU_ACSm < 1000
SU_ACSm: FU_ACSm + BU_ACSm < 1
UF_ACSm: F_ACSm -1000 FU_ACSm < 0
UR_ACSm: R_ACSm -1000 BU_ACSm < 0
G_ACt2r: +0.59209 LC_ac_c -0.59209 LC_ac_e - DG_ACt2r + DGo_ACt2r = 0
FU_ACt2r: DG_ACt2r +1000 FU_ACt2r < 1000
BU_ACt2r: - DG_ACt2r +1000 BU_ACt2r < 1000
SU_ACt2r: FU_ACt2r + BU_ACt2r < 1
UF_ACt2r: F_ACt2r -1000 FU_ACt2r < 0
UR_ACt2r: R_ACt2r -1000 BU_ACt2r < 0
G_ACtm: -0.59209 LC_ac_c +0.59209 LC_ac_m - DG_ACtm + DGo_ACtm = 0
FU_ACtm: DG_ACtm +1000 FU_ACtm < 1000
BU_ACtm: - DG_ACtm +1000 BU_ACtm < 1000
SU_ACtm: FU_ACtm + BU_ACtm < 1
UF_ACtm: F_ACtm -1000 FU_ACtm < 0
UR_ACtm: R_ACtm -1000 BU_ACtm < 0
G_ACtr: +0.59209 LC_ac_c -0.59209 LC_ac_e - DG_ACtr + DGo_ACtr = 0
FU_ACtr: DG_ACtr +1000 FU_ACtr < 1000
BU_ACtr: - DG_ACtr +1000 BU_ACtr < 1000
SU_ACtr: FU_ACtr + BU_ACtr < 1
UF_ACtr: F_ACtr -1000 FU_ACtr < 0
UR_ACtr: R_ACtr -1000 BU_ACtr < 0
G_ADK1: +1.1842 LC_adp_c -0.59209 LC_amp_c -0.59209 LC_atp_c - DG_ADK1 + DGo_ADK1 = 0
FU_ADK1: DG_ADK1 +1000 FU_ADK1 < 1000
BU_ADK1: - DG_ADK1 +1000 BU_ADK1 < 1000
SU_ADK1: FU_ADK1 + BU_ADK1 < 1
UF_ADK1: F_ADK1 -1000 FU_ADK1 < 0
UR_ADK1: R_ADK1 -1000 BU_ADK1 < 0
G_ADK1m: +1.1842 LC_adp_m -0.59209 LC_amp_m -0.59209 LC_atp_m - DG_ADK1m + DGo_ADK1m = 0
FU_ADK1m: DG_ADK1m +1000 FU_ADK1m < 1000
BU_ADK1m: - DG_ADK1m +1000 BU_ADK1m < 1000
SU_ADK1m: FU_ADK1m + BU_ADK1m < 1
UF_ADK1m: F_ADK1m -1000 FU_ADK1m < 0
UR_ADK1m: R_ADK1m -1000 BU_ADK1m < 0
G_AKGDam: -0.59209 LC_akg_m +0.59209 LC_co2_m -0.59209 LC_lpam_m +0.59209 LC_sdhlam_m - DG_AKGDam + DGo_AKGDam = 0
FU_AKGDam: DG_AKGDam +1000 FU_AKGDam < 1000
BU_AKGDam: - DG_AKGDam +1000 BU_AKGDam < 1000
SU_AKGDam: FU_AKGDam + BU_AKGDam < 1
UF_AKGDam: F_AKGDam -1000 FU_AKGDam < 0
UR_AKGDam: R_AKGDam -1000 BU_AKGDam < 0
G_AKGDbm: -0.59209 LC_coa_m +0.59209 LC_dhlam_m -0.59209 LC_sdhlam_m +0.59209 LC_succoa_m - DG_AKGDbm + DGo_AKGDbm = 0
FU_AKGDbm: DG_AKGDbm +1000 FU_AKGDbm < 1000
BU_AKGDbm: - DG_AKGDbm +1000 BU_AKGDbm < 1000
SU_AKGDbm: FU_AKGDbm + BU_AKGDbm < 1
UF_AKGDbm: F_AKGDbm -1000 FU_AKGDbm < 0
UR_AKGDbm: R_AKGDbm -1000 BU_AKGDbm < 0
G_AKGMAL: -0.59209 LC_akg_c +0.59209 LC_mal-L_c +0.59209 LC_akg_e -0.59209 LC_mal-L_e - DG_AKGMAL + DGo_AKGMAL = 0
FU_AKGMAL: DG_AKGMAL +1000 FU_AKGMAL < 1000
BU_AKGMAL: - DG_AKGMAL +1000 BU_AKGMAL < 1000
SU_AKGMAL: FU_AKGMAL + BU_AKGMAL < 1
UF_AKGMAL: F_AKGMAL -1000 FU_AKGMAL < 0
UR_AKGMAL: R_AKGMAL -1000 BU_AKGMAL < 0
G_AKGt2r: +0.59209 LC_akg_c -0.59209 LC_akg_e - DG_AKGt2r + DGo_AKGt2r = 0
FU_AKGt2r: DG_AKGt2r +1000 FU_AKGt2r < 1000
BU_AKGt2r: - DG_AKGt2r +1000 BU_AKGt2r < 1000
SU_AKGt2r: FU_AKGt2r + BU_AKGt2r < 1
UF_AKGt2r: F_AKGt2r -1000 FU_AKGt2r < 0
UR_AKGt2r: R_AKGt2r -1000 BU_AKGt2r < 0
G_ALATA_L: -0.59209 LC_akg_c -0.59209 LC_ala-L_c +0.59209 LC_glu-L_c +0.59209 LC_pyr_c - DG_ALATA_L + DGo_ALATA_L = 0
FU_ALATA_L: DG_ALATA_L +1000 FU_ALATA_L < 1000
BU_ALATA_L: - DG_ALATA_L +1000 BU_ALATA_L < 1000
SU_ALATA_L: FU_ALATA_L + BU_ALATA_L < 1
UF_ALATA_L: F_ALATA_L -1000 FU_ALATA_L < 0
UR_ALATA_L: R_ALATA_L -1000 BU_ALATA_L < 0
G_ALATA_Lm: -0.59209 LC_akg_m -0.59209 LC_ala-L_m +0.59209 LC_glu-L_m +0.59209 LC_pyr_m - DG_ALATA_Lm + DGo_ALATA_Lm = 0
FU_ALATA_Lm: DG_ALATA_Lm +1000 FU_ALATA_Lm < 1000
BU_ALATA_Lm: - DG_ALATA_Lm +1000 BU_ALATA_Lm < 1000
SU_ALATA_Lm: FU_ALATA_Lm + BU_ALATA_Lm < 1
UF_ALATA_Lm: F_ALATA_Lm -1000 FU_ALATA_Lm < 0
UR_ALATA_Lm: R_ALATA_Lm -1000 BU_ALATA_Lm < 0
G_ALAt2r: +0.59209 LC_ala-L_c -0.59209 LC_ala-L_e - DG_ALAt2r + DGo_ALAt2r = 0
FU_ALAt2r: DG_ALAt2r +1000 FU_ALAt2r < 1000
BU_ALAt2r: - DG_ALAt2r +1000 BU_ALAt2r < 1000
SU_ALAt2r: FU_ALAt2r + BU_ALAt2r < 1
UF_ALAt2r: F_ALAt2r -1000 FU_ALAt2r < 0
UR_ALAt2r: R_ALAt2r -1000 BU_ALAt2r < 0
G_ALAtmi: +0.59209 LC_ala-L_c -0.59209 LC_ala-L_m - DG_ALAtmi + DGo_ALAtmi = 0
FU_ALAtmi: DG_ALAtmi +1000 FU_ALAtmi < 1000
BU_ALAtmi: - DG_ALAtmi +1000 BU_ALAtmi < 1000
SU_ALAtmi: FU_ALAtmi + BU_ALAtmi < 1
UF_ALAtmi: F_ALAtmi -1000 FU_ALAtmi < 0
UR_ALAtmi: R_ALAtmi -1000 BU_ALAtmi < 0
G_ALCD2if: +0.59209 LC_acald_c -0.59209 LC_etoh_c -0.59209 LC_nad_c +0.59209 LC_nadh_c - DG_ALCD2if + DGo_ALCD2if = 0
FU_ALCD2if: DG_ALCD2if +1000 FU_ALCD2if < 1000
BU_ALCD2if: - DG_ALCD2if +1000 BU_ALCD2if < 1000
SU_ALCD2if: FU_ALCD2if + BU_ALCD2if < 1
UF_ALCD2if: F_ALCD2if -1000 FU_ALCD2if < 0
UR_ALCD2if: R_ALCD2if -1000 BU_ALCD2if < 0
G_ALCD2ir: -0.59209 LC_acald_c +0.59209 LC_etoh_c +0.59209 LC_nad_c -0.59209 LC_nadh_c - DG_ALCD2ir + DGo_ALCD2ir = 0
FU_ALCD2ir: DG_ALCD2ir +1000 FU_ALCD2ir < 1000
BU_ALCD2ir: - DG_ALCD2ir +1000 BU_ALCD2ir < 1000
SU_ALCD2ir: FU_ALCD2ir + BU_ALCD2ir < 1
UF_ALCD2ir: F_ALCD2ir -1000 FU_ALCD2ir < 0
UR_ALCD2ir: R_ALCD2ir -1000 BU_ALCD2ir < 0
G_ALCD2irm: -0.59209 LC_acald_m +0.59209 LC_etoh_m +0.59209 LC_nad_m -0.59209 LC_nadh_m - DG_ALCD2irm + DGo_ALCD2irm = 0
FU_ALCD2irm: DG_ALCD2irm +1000 FU_ALCD2irm < 1000
BU_ALCD2irm: - DG_ALCD2irm +1000 BU_ALCD2irm < 1000
SU_ALCD2irm: FU_ALCD2irm + BU_ALCD2irm < 1
UF_ALCD2irm: F_ALCD2irm -1000 FU_ALCD2irm < 0
UR_ALCD2irm: R_ALCD2irm -1000 BU_ALCD2irm < 0
G_ALCD2x: +0.59209 LC_acald_c -0.59209 LC_etoh_c -0.59209 LC_nad_c +0.59209 LC_nadh_c - DG_ALCD2x + DGo_ALCD2x = 0
FU_ALCD2x: DG_ALCD2x +1000 FU_ALCD2x < 1000
BU_ALCD2x: - DG_ALCD2x +1000 BU_ALCD2x < 1000
SU_ALCD2x: FU_ALCD2x + BU_ALCD2x < 1
UF_ALCD2x: F_ALCD2x -1000 FU_ALCD2x < 0
UR_ALCD2x: R_ALCD2x -1000 BU_ALCD2x < 0
G_ALDD2xm: +0.59209 LC_ac_m -0.59209 LC_acald_m -0.59209 LC_nad_m +0.59209 LC_nadh_m - DG_ALDD2xm + DGo_ALDD2xm = 0
FU_ALDD2xm: DG_ALDD2xm +1000 FU_ALDD2xm < 1000
BU_ALDD2xm: - DG_ALDD2xm +1000 BU_ALDD2xm < 1000
SU_ALDD2xm: FU_ALDD2xm + BU_ALDD2xm < 1
UF_ALDD2xm: F_ALDD2xm -1000 FU_ALDD2xm < 0
UR_ALDD2xm: R_ALDD2xm -1000 BU_ALDD2xm < 0
G_ALDD2y: +0.59209 LC_ac_c -0.59209 LC_acald_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_ALDD2y + DGo_ALDD2y = 0
FU_ALDD2y: DG_ALDD2y +1000 FU_ALDD2y < 1000
BU_ALDD2y: - DG_ALDD2y +1000 BU_ALDD2y < 1000
SU_ALDD2y: FU_ALDD2y + BU_ALDD2y < 1
UF_ALDD2y: F_ALDD2y -1000 FU_ALDD2y < 0
UR_ALDD2y: R_ALDD2y -1000 BU_ALDD2y < 0
G_ALDD2ym: +0.59209 LC_ac_m -0.59209 LC_acald_m -0.59209 LC_nadp_m +0.59209 LC_nadph_m - DG_ALDD2ym + DGo_ALDD2ym = 0
FU_ALDD2ym: DG_ALDD2ym +1000 FU_ALDD2ym < 1000
BU_ALDD2ym: - DG_ALDD2ym +1000 BU_ALDD2ym < 1000
SU_ALDD2ym: FU_ALDD2ym + BU_ALDD2ym < 1
UF_ALDD2ym: F_ALDD2ym -1000 FU_ALDD2ym < 0
UR_ALDD2ym: R_ALDD2ym -1000 BU_ALDD2ym < 0
G_ARGt2r: +0.59209 LC_arg-L_c -0.59209 LC_arg-L_e - DG_ARGt2r + DGo_ARGt2r = 0
FU_ARGt2r: DG_ARGt2r +1000 FU_ARGt2r < 1000
BU_ARGt2r: - DG_ARGt2r +1000 BU_ARGt2r < 1000
SU_ARGt2r: FU_ARGt2r + BU_ARGt2r < 1
UF_ARGt2r: F_ARGt2r -1000 FU_ARGt2r < 0
UR_ARGt2r: R_ARGt2r -1000 BU_ARGt2r < 0
G_ASNt2r: +0.59209 LC_asn-L_c -0.59209 LC_asn-L_e - DG_ASNt2r + DGo_ASNt2r = 0
FU_ASNt2r: DG_ASNt2r +1000 FU_ASNt2r < 1000
BU_ASNt2r: - DG_ASNt2r +1000 BU_ASNt2r < 1000
SU_ASNt2r: FU_ASNt2r + BU_ASNt2r < 1
UF_ASNt2r: F_ASNt2r -1000 FU_ASNt2r < 0
UR_ASNt2r: R_ASNt2r -1000 BU_ASNt2r < 0
G_ASPGLU2m: +0.59209 LC_asp-L_c -0.59209 LC_glu-L_c -0.59209 LC_asp-L_m +0.59209 LC_glu-L_m - DG_ASPGLU2m + DGo_ASPGLU2m = 0
FU_ASPGLU2m: DG_ASPGLU2m +1000 FU_ASPGLU2m < 1000
BU_ASPGLU2m: - DG_ASPGLU2m +1000 BU_ASPGLU2m < 1000
SU_ASPGLU2m: FU_ASPGLU2m + BU_ASPGLU2m < 1
UF_ASPGLU2m: F_ASPGLU2m -1000 FU_ASPGLU2m < 0
UR_ASPGLU2m: R_ASPGLU2m -1000 BU_ASPGLU2m < 0
G_ASPTA: -0.59209 LC_akg_c -0.59209 LC_asp-L_c +0.59209 LC_glu-L_c +0.59209 LC_oaa_c - DG_ASPTA + DGo_ASPTA = 0
FU_ASPTA: DG_ASPTA +1000 FU_ASPTA < 1000
BU_ASPTA: - DG_ASPTA +1000 BU_ASPTA < 1000
SU_ASPTA: FU_ASPTA + BU_ASPTA < 1
UF_ASPTA: F_ASPTA -1000 FU_ASPTA < 0
UR_ASPTA: R_ASPTA -1000 BU_ASPTA < 0
G_ASPTAm: -0.59209 LC_akg_m -0.59209 LC_asp-L_m +0.59209 LC_glu-L_m +0.59209 LC_oaa_m - DG_ASPTAm + DGo_ASPTAm = 0
FU_ASPTAm: DG_ASPTAm +1000 FU_ASPTAm < 1000
BU_ASPTAm: - DG_ASPTAm +1000 BU_ASPTAm < 1000
SU_ASPTAm: FU_ASPTAm + BU_ASPTAm < 1
UF_ASPTAm: F_ASPTAm -1000 FU_ASPTAm < 0
UR_ASPTAm: R_ASPTAm -1000 BU_ASPTAm < 0
G_ASPt2m: -0.59209 LC_asp-L_c +0.59209 LC_asp-L_m - DG_ASPt2m + DGo_ASPt2m = 0
FU_ASPt2m: DG_ASPt2m +1000 FU_ASPt2m < 1000
BU_ASPt2m: - DG_ASPt2m +1000 BU_ASPt2m < 1000
SU_ASPt2m: FU_ASPt2m + BU_ASPt2m < 1
UF_ASPt2m: F_ASPt2m -1000 FU_ASPt2m < 0
UR_ASPt2m: R_ASPt2m -1000 BU_ASPt2m < 0
G_ASPt2r: +0.59209 LC_asp-L_c -0.59209 LC_asp-L_e - DG_ASPt2r + DGo_ASPt2r = 0
FU_ASPt2r: DG_ASPt2r +1000 FU_ASPt2r < 1000
BU_ASPt2r: - DG_ASPt2r +1000 BU_ASPt2r < 1000
SU_ASPt2r: FU_ASPt2r + BU_ASPt2r < 1
UF_ASPt2r: F_ASPt2r -1000 FU_ASPt2r < 0
UR_ASPt2r: R_ASPt2r -1000 BU_ASPt2r < 0
G_ATPH1: +0.59209 LC_amp_c -0.59209 LC_atp_c +1.1842 LC_pi_c - DG_ATPH1 + DGo_ATPH1 = 0
FU_ATPH1: DG_ATPH1 +1000 FU_ATPH1 < 1000
BU_ATPH1: - DG_ATPH1 +1000 BU_ATPH1 < 1000
SU_ATPH1: FU_ATPH1 + BU_ATPH1 < 1
UF_ATPH1: F_ATPH1 -1000 FU_ATPH1 < 0
UR_ATPH1: R_ATPH1 -1000 BU_ATPH1 < 0
G_ATPM: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_pi_c - DG_ATPM + DGo_ATPM = 0
FU_ATPM: DG_ATPM +1000 FU_ATPM < 1000
BU_ATPM: - DG_ATPM +1000 BU_ATPM < 1000
SU_ATPM: FU_ATPM + BU_ATPM < 1
UF_ATPM: F_ATPM -1000 FU_ATPM < 0
UR_ATPM: R_ATPM -1000 BU_ATPM < 0
G_ATPS: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_pi_c - DG_ATPS + DGo_ATPS = 0
FU_ATPS: DG_ATPS +1000 FU_ATPS < 1000
BU_ATPS: - DG_ATPS +1000 BU_ATPS < 1000
SU_ATPS: FU_ATPS + BU_ATPS < 1
UF_ATPS: F_ATPS -1000 FU_ATPS < 0
UR_ATPS: R_ATPS -1000 BU_ATPS < 0
G_ATPS3m: -0.59209 LC_adp_m +0.59209 LC_atp_m -0.59209 LC_pi_m - DG_ATPS3m + DGo_ATPS3m = 0
FU_ATPS3m: DG_ATPS3m +1000 FU_ATPS3m < 1000
BU_ATPS3m: - DG_ATPS3m +1000 BU_ATPS3m < 1000
SU_ATPS3m: FU_ATPS3m + BU_ATPS3m < 1
UF_ATPS3m: F_ATPS3m -1000 FU_ATPS3m < 0
UR_ATPS3m: R_ATPS3m -1000 BU_ATPS3m < 0
G_ATPtm-H: -0.59209 LC_adp_c +0.59209 LC_atp_c +0.59209 LC_adp_m -0.59209 LC_atp_m - DG_ATPtm-H + DGo_ATPtm-H = 0
FU_ATPtm-H: DG_ATPtm-H +1000 FU_ATPtm-H < 1000
BU_ATPtm-H: - DG_ATPtm-H +1000 BU_ATPtm-H < 1000
SU_ATPtm-H: FU_ATPtm-H + BU_ATPtm-H < 1
UF_ATPtm-H: F_ATPtm-H -1000 FU_ATPtm-H < 0
UR_ATPtm-H: R_ATPtm-H -1000 BU_ATPtm-H < 0
G_CITt2r: +0.59209 LC_cit_c -0.59209 LC_cit_e - DG_CITt2r + DGo_CITt2r = 0
FU_CITt2r: DG_CITt2r +1000 FU_CITt2r < 1000
BU_CITt2r: - DG_CITt2r +1000 BU_CITt2r < 1000
SU_CITt2r: FU_CITt2r + BU_CITt2r < 1
UF_CITt2r: F_CITt2r -1000 FU_CITt2r < 0
UR_CITt2r: R_CITt2r -1000 BU_CITt2r < 0
G_CITtam: -0.59209 LC_cit_c +0.59209 LC_mal-L_c +0.59209 LC_cit_m -0.59209 LC_mal-L_m - DG_CITtam + DGo_CITtam = 0
FU_CITtam: DG_CITtam +1000 FU_CITtam < 1000
BU_CITtam: - DG_CITtam +1000 BU_CITtam < 1000
SU_CITtam: FU_CITtam + BU_CITtam < 1
UF_CITtam: F_CITtam -1000 FU_CITtam < 0
UR_CITtam: R_CITtam -1000 BU_CITtam < 0
G_CITtcm: -0.59209 LC_cit_c +0.59209 LC_icit_c +0.59209 LC_cit_m -0.59209 LC_icit_m - DG_CITtcm + DGo_CITtcm = 0
FU_CITtcm: DG_CITtcm +1000 FU_CITtcm < 1000
BU_CITtcm: - DG_CITtcm +1000 BU_CITtcm < 1000
SU_CITtcm: FU_CITtcm + BU_CITtcm < 1
UF_CITtcm: F_CITtcm -1000 FU_CITtcm < 0
UR_CITtcm: R_CITtcm -1000 BU_CITtcm < 0
G_CO2t: +0.59209 LC_co2_c -0.59209 LC_co2_e - DG_CO2t + DGo_CO2t = 0
FU_CO2t: DG_CO2t +1000 FU_CO2t < 1000
BU_CO2t: - DG_CO2t +1000 BU_CO2t < 1000
SU_CO2t: FU_CO2t + BU_CO2t < 1
UF_CO2t: F_CO2t -1000 FU_CO2t < 0
UR_CO2t: R_CO2t -1000 BU_CO2t < 0
G_CO2tm: -0.59209 LC_co2_c +0.59209 LC_co2_m - DG_CO2tm + DGo_CO2tm = 0
FU_CO2tm: DG_CO2tm +1000 FU_CO2tm < 1000
BU_CO2tm: - DG_CO2tm +1000 BU_CO2tm < 1000
SU_CO2tm: FU_CO2tm + BU_CO2tm < 1
UF_CO2tm: F_CO2tm -1000 FU_CO2tm < 0
UR_CO2tm: R_CO2tm -1000 BU_CO2tm < 0
G_CSm: -0.59209 LC_accoa_m +0.59209 LC_cit_m +0.59209 LC_coa_m -0.59209 LC_oaa_m - DG_CSm + DGo_CSm = 0
FU_CSm: DG_CSm +1000 FU_CSm < 1000
BU_CSm: - DG_CSm +1000 BU_CSm < 1000
SU_CSm: FU_CSm + BU_CSm < 1
UF_CSm: F_CSm -1000 FU_CSm < 0
UR_CSm: R_CSm -1000 BU_CSm < 0
G_CYOOm: +2.3684 LC_ficytc_m -2.3684 LC_focytc_m -0.59209 LC_o2_m - DG_CYOOm + DGo_CYOOm = 0
FU_CYOOm: DG_CYOOm +1000 FU_CYOOm < 1000
BU_CYOOm: - DG_CYOOm +1000 BU_CYOOm < 1000
SU_CYOOm: FU_CYOOm + BU_CYOOm < 1
UF_CYOOm: F_CYOOm -1000 FU_CYOOm < 0
UR_CYOOm: R_CYOOm -1000 BU_CYOOm < 0
G_CYOR_u6m: -1.1842 LC_ficytc_m +1.1842 LC_focytc_m +0.59209 LC_q6_m -0.59209 LC_q6h2_m - DG_CYOR_u6m + DGo_CYOR_u6m = 0
FU_CYOR_u6m: DG_CYOR_u6m +1000 FU_CYOR_u6m < 1000
BU_CYOR_u6m: - DG_CYOR_u6m +1000 BU_CYOR_u6m < 1000
SU_CYOR_u6m: FU_CYOR_u6m + BU_CYOR_u6m < 1
UF_CYOR_u6m: F_CYOR_u6m -1000 FU_CYOR_u6m < 0
UR_CYOR_u6m: R_CYOR_u6m -1000 BU_CYOR_u6m < 0
G_CYSt2r: +0.59209 LC_cys-L_c -0.59209 LC_cys-L_e - DG_CYSt2r + DGo_CYSt2r = 0
FU_CYSt2r: DG_CYSt2r +1000 FU_CYSt2r < 1000
BU_CYSt2r: - DG_CYSt2r +1000 BU_CYSt2r < 1000
SU_CYSt2r: FU_CYSt2r + BU_CYSt2r < 1
UF_CYSt2r: F_CYSt2r -1000 FU_CYSt2r < 0
UR_CYSt2r: R_CYSt2r -1000 BU_CYSt2r < 0
G_D-LACt2m: -0.59209 LC_lac-D_c +0.59209 LC_pyr_c +0.59209 LC_lac-D_m -0.59209 LC_pyr_m - DG_D-LACt2m + DGo_D-LACt2m = 0
FU_D-LACt2m: DG_D-LACt2m +1000 FU_D-LACt2m < 1000
BU_D-LACt2m: - DG_D-LACt2m +1000 BU_D-LACt2m < 1000
SU_D-LACt2m: FU_D-LACt2m + BU_D-LACt2m < 1
UF_D-LACt2m: F_D-LACt2m -1000 FU_D-LACt2m < 0
UR_D-LACt2m: R_D-LACt2m -1000 BU_D-LACt2m < 0
G_D-LACtm: -0.59209 LC_lac-D_c +0.59209 LC_lac-D_m - DG_D-LACtm + DGo_D-LACtm = 0
FU_D-LACtm: DG_D-LACtm +1000 FU_D-LACtm < 1000
BU_D-LACtm: - DG_D-LACtm +1000 BU_D-LACtm < 1000
SU_D-LACtm: FU_D-LACtm + BU_D-LACtm < 1
UF_D-LACtm: F_D-LACtm -1000 FU_D-LACtm < 0
UR_D-LACtm: R_D-LACtm -1000 BU_D-LACtm < 0
G_DB4PS: +0.59209 LC_db4p_c +0.59209 LC_for_c -0.59209 LC_ru5p-D_c - DG_DB4PS + DGo_DB4PS = 0
FU_DB4PS: DG_DB4PS +1000 FU_DB4PS < 1000
BU_DB4PS: - DG_DB4PS +1000 BU_DB4PS < 1000
SU_DB4PS: FU_DB4PS + BU_DB4PS < 1
UF_DB4PS: F_DB4PS -1000 FU_DB4PS < 0
UR_DB4PS: R_DB4PS -1000 BU_DB4PS < 0
G_ENO: -0.59209 LC_2pg_c +0.59209 LC_pep_c - DG_ENO + DGo_ENO = 0
FU_ENO: DG_ENO +1000 FU_ENO < 1000
BU_ENO: - DG_ENO +1000 BU_ENO < 1000
SU_ENO: FU_ENO + BU_ENO < 1
UF_ENO: F_ENO -1000 FU_ENO < 0
UR_ENO: R_ENO -1000 BU_ENO < 0
G_ERGSTt: +0.59209 LC_ergst_c -0.59209 LC_ergst_e - DG_ERGSTt + DGo_ERGSTt = 0
FU_ERGSTt: DG_ERGSTt +1000 FU_ERGSTt < 1000
BU_ERGSTt: - DG_ERGSTt +1000 BU_ERGSTt < 1000
SU_ERGSTt: FU_ERGSTt + BU_ERGSTt < 1
UF_ERGSTt: F_ERGSTt -1000 FU_ERGSTt < 0
UR_ERGSTt: R_ERGSTt -1000 BU_ERGSTt < 0
G_ERGSTter: +0.59209 LC_ergst_c -0.59209 LC_ergst_r - DG_ERGSTter + DGo_ERGSTter = 0
FU_ERGSTter: DG_ERGSTter +1000 FU_ERGSTter < 1000
BU_ERGSTter: - DG_ERGSTter +1000 BU_ERGSTter < 1000
SU_ERGSTter: FU_ERGSTter + BU_ERGSTter < 1
UF_ERGSTter: F_ERGSTter -1000 FU_ERGSTter < 0
UR_ERGSTter: R_ERGSTter -1000 BU_ERGSTter < 0
G_ETOHt: +0.59209 LC_etoh_c -0.59209 LC_etoh_e - DG_ETOHt + DGo_ETOHt = 0
FU_ETOHt: DG_ETOHt +1000 FU_ETOHt < 1000
BU_ETOHt: - DG_ETOHt +1000 BU_ETOHt < 1000
SU_ETOHt: FU_ETOHt + BU_ETOHt < 1
UF_ETOHt: F_ETOHt -1000 FU_ETOHt < 0
UR_ETOHt: R_ETOHt -1000 BU_ETOHt < 0
G_ETOHtm: -0.59209 LC_etoh_c +0.59209 LC_etoh_m - DG_ETOHtm + DGo_ETOHtm = 0
FU_ETOHtm: DG_ETOHtm +1000 FU_ETOHtm < 1000
BU_ETOHtm: - DG_ETOHtm +1000 BU_ETOHtm < 1000
SU_ETOHtm: FU_ETOHtm + BU_ETOHtm < 1
UF_ETOHtm: F_ETOHtm -1000 FU_ETOHtm < 0
UR_ETOHtm: R_ETOHtm -1000 BU_ETOHtm < 0
G_FBA: +0.59209 LC_dhap_c -0.59209 LC_fdp_c +0.59209 LC_g3p_c - DG_FBA + DGo_FBA = 0
FU_FBA: DG_FBA +1000 FU_FBA < 1000
BU_FBA: - DG_FBA +1000 BU_FBA < 1000
SU_FBA: FU_FBA + BU_FBA < 1
UF_FBA: F_FBA -1000 FU_FBA < 0
UR_FBA: R_FBA -1000 BU_FBA < 0
G_FBA3: +0.59209 LC_dhap_c +0.59209 LC_e4p_c -0.59209 LC_s17bp_c - DG_FBA3 + DGo_FBA3 = 0
FU_FBA3: DG_FBA3 +1000 FU_FBA3 < 1000
BU_FBA3: - DG_FBA3 +1000 BU_FBA3 < 1000
SU_FBA3: FU_FBA3 + BU_FBA3 < 1
UF_FBA3: F_FBA3 -1000 FU_FBA3 < 0
UR_FBA3: R_FBA3 -1000 BU_FBA3 < 0
G_FBP: +0.59209 LC_f6p_c -0.59209 LC_fdp_c +0.59209 LC_pi_c - DG_FBP + DGo_FBP = 0
FU_FBP: DG_FBP +1000 FU_FBP < 1000
BU_FBP: - DG_FBP +1000 BU_FBP < 1000
SU_FBP: FU_FBP + BU_FBP < 1
UF_FBP: F_FBP -1000 FU_FBP < 0
UR_FBP: R_FBP -1000 BU_FBP < 0
G_FDH: +0.59209 LC_co2_c -0.59209 LC_for_c -0.59209 LC_nad_c +0.59209 LC_nadh_c - DG_FDH + DGo_FDH = 0
FU_FDH: DG_FDH +1000 FU_FDH < 1000
BU_FDH: - DG_FDH +1000 BU_FDH < 1000
SU_FDH: FU_FDH + BU_FDH < 1
UF_FDH: F_FDH -1000 FU_FDH < 0
UR_FDH: R_FDH -1000 BU_FDH < 0
G_FDNG: +0.59209 LC_co2_c -0.59209 LC_for_c -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_FDNG + DGo_FDNG = 0
FU_FDNG: DG_FDNG +1000 FU_FDNG < 1000
BU_FDNG: - DG_FDNG +1000 BU_FDNG < 1000
SU_FDNG: FU_FDNG + BU_FDNG < 1
UF_FDNG: F_FDNG -1000 FU_FDNG < 0
UR_FDNG: R_FDNG -1000 BU_FDNG < 0
G_FORt: +0.59209 LC_for_c -0.59209 LC_for_e - DG_FORt + DGo_FORt = 0
FU_FORt: DG_FORt +1000 FU_FORt < 1000
BU_FORt: - DG_FORt +1000 BU_FORt < 1000
SU_FORt: FU_FORt + BU_FORt < 1
UF_FORt: F_FORt -1000 FU_FORt < 0
UR_FORt: R_FORt -1000 BU_FORt < 0
G_FRDcm: -0.59209 LC_fum_c +0.59209 LC_succ_c +0.59209 LC_fad_m -0.59209 LC_fadh2_m - DG_FRDcm + DGo_FRDcm = 0
FU_FRDcm: DG_FRDcm +1000 FU_FRDcm < 1000
BU_FRDcm: - DG_FRDcm +1000 BU_FRDcm < 1000
SU_FRDcm: FU_FRDcm + BU_FRDcm < 1
UF_FRDcm: F_FRDcm -1000 FU_FRDcm < 0
UR_FRDcm: R_FRDcm -1000 BU_FRDcm < 0
G_FRDm: +0.59209 LC_fad_m -0.59209 LC_fadh2_m -0.59209 LC_fum_m +0.59209 LC_succ_m - DG_FRDm + DGo_FRDm = 0
FU_FRDm: DG_FRDm +1000 FU_FRDm < 1000
BU_FRDm: - DG_FRDm +1000 BU_FRDm < 1000
SU_FRDm: FU_FRDm + BU_FRDm < 1
UF_FRDm: F_FRDm -1000 FU_FRDm < 0
UR_FRDm: R_FRDm -1000 BU_FRDm < 0
G_FUM: -0.59209 LC_fum_c +0.59209 LC_mal-L_c - DG_FUM + DGo_FUM = 0
FU_FUM: DG_FUM +1000 FU_FUM < 1000
BU_FUM: - DG_FUM +1000 BU_FUM < 1000
SU_FUM: FU_FUM + BU_FUM < 1
UF_FUM: F_FUM -1000 FU_FUM < 0
UR_FUM: R_FUM -1000 BU_FUM < 0
G_FUMm: -0.59209 LC_fum_m +0.59209 LC_mal-L_m - DG_FUMm + DGo_FUMm = 0
FU_FUMm: DG_FUMm +1000 FU_FUMm < 1000
BU_FUMm: - DG_FUMm +1000 BU_FUMm < 1000
SU_FUMm: FU_FUMm + BU_FUMm < 1
UF_FUMm: F_FUMm -1000 FU_FUMm < 0
UR_FUMm: R_FUMm -1000 BU_FUMm < 0
G_FUMt2r: +0.59209 LC_fum_c -0.59209 LC_fum_e - DG_FUMt2r + DGo_FUMt2r = 0
FU_FUMt2r: DG_FUMt2r +1000 FU_FUMt2r < 1000
BU_FUMt2r: - DG_FUMt2r +1000 BU_FUMt2r < 1000
SU_FUMt2r: FU_FUMt2r + BU_FUMt2r < 1
UF_FUMt2r: F_FUMt2r -1000 FU_FUMt2r < 0
UR_FUMt2r: R_FUMt2r -1000 BU_FUMt2r < 0
G_G6PDH2: +0.59209 LC_6pgl_c -0.59209 LC_g6p_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_G6PDH2 + DGo_G6PDH2 = 0
FU_G6PDH2: DG_G6PDH2 +1000 FU_G6PDH2 < 1000
BU_G6PDH2: - DG_G6PDH2 +1000 BU_G6PDH2 < 1000
SU_G6PDH2: FU_G6PDH2 + BU_G6PDH2 < 1
UF_G6PDH2: F_G6PDH2 -1000 FU_G6PDH2 < 0
UR_G6PDH2: R_G6PDH2 -1000 BU_G6PDH2 < 0
G_G6PDH2er: +0.59209 LC_6pgl_r -0.59209 LC_g6p_r -0.59209 LC_nadp_r +0.59209 LC_nadph_r - DG_G6PDH2er + DGo_G6PDH2er = 0
FU_G6PDH2er: DG_G6PDH2er +1000 FU_G6PDH2er < 1000
BU_G6PDH2er: - DG_G6PDH2er +1000 BU_G6PDH2er < 1000
SU_G6PDH2er: FU_G6PDH2er + BU_G6PDH2er < 1
UF_G6PDH2er: F_G6PDH2er -1000 FU_G6PDH2er < 0
UR_G6PDH2er: R_G6PDH2er -1000 BU_G6PDH2er < 0
G_G6PI: -0.59209 LC_g6p_c +0.59209 LC_g6p-B_c - DG_G6PI + DGo_G6PI = 0
FU_G6PI: DG_G6PI +1000 FU_G6PI < 1000
BU_G6PI: - DG_G6PI +1000 BU_G6PI < 1000
SU_G6PI: FU_G6PI + BU_G6PI < 1
UF_G6PI: F_G6PI -1000 FU_G6PI < 0
UR_G6PI: R_G6PI -1000 BU_G6PI < 0
G_G6PI3: +0.59209 LC_f6p_c -0.59209 LC_g6p-B_c - DG_G6PI3 + DGo_G6PI3 = 0
FU_G6PI3: DG_G6PI3 +1000 FU_G6PI3 < 1000
BU_G6PI3: - DG_G6PI3 +1000 BU_G6PI3 < 1000
SU_G6PI3: FU_G6PI3 + BU_G6PI3 < 1
UF_G6PI3: F_G6PI3 -1000 FU_G6PI3 < 0
UR_G6PI3: R_G6PI3 -1000 BU_G6PI3 < 0
G_G6Pter: -0.59209 LC_g6p_c +0.59209 LC_g6p_r - DG_G6Pter + DGo_G6Pter = 0
FU_G6Pter: DG_G6Pter +1000 FU_G6Pter < 1000
BU_G6Pter: - DG_G6Pter +1000 BU_G6Pter < 1000
SU_G6Pter: FU_G6Pter + BU_G6Pter < 1
UF_G6Pter: F_G6Pter -1000 FU_G6Pter < 0
UR_G6Pter: R_G6Pter -1000 BU_G6Pter < 0
G_GAPD: +0.59209 LC_13dpg_c -0.59209 LC_g3p_c -0.59209 LC_nad_c +0.59209 LC_nadh_c -0.59209 LC_pi_c - DG_GAPD + DGo_GAPD = 0
FU_GAPD: DG_GAPD +1000 FU_GAPD < 1000
BU_GAPD: - DG_GAPD +1000 BU_GAPD < 1000
SU_GAPD: FU_GAPD + BU_GAPD < 1
UF_GAPD: F_GAPD -1000 FU_GAPD < 0
UR_GAPD: R_GAPD -1000 BU_GAPD < 0
G_GCC2cm: -0.59209 LC_dhlam_m +0.59209 LC_lpam_m -0.59209 LC_nad_m +0.59209 LC_nadh_m - DG_GCC2cm + DGo_GCC2cm = 0
FU_GCC2cm: DG_GCC2cm +1000 FU_GCC2cm < 1000
BU_GCC2cm: - DG_GCC2cm +1000 BU_GCC2cm < 1000
SU_GCC2cm: FU_GCC2cm + BU_GCC2cm < 1
UF_GCC2cm: F_GCC2cm -1000 FU_GCC2cm < 0
UR_GCC2cm: R_GCC2cm -1000 BU_GCC2cm < 0
G_GLCt1: +0.59209 LC_glc-D_c -0.59209 LC_glc-D_e - DG_GLCt1 + DGo_GLCt1 = 0
FU_GLCt1: DG_GLCt1 +1000 FU_GLCt1 < 1000
BU_GLCt1: - DG_GLCt1 +1000 BU_GLCt1 < 1000
SU_GLCt1: FU_GLCt1 + BU_GLCt1 < 1
UF_GLCt1: F_GLCt1 -1000 FU_GLCt1 < 0
UR_GLCt1: R_GLCt1 -1000 BU_GLCt1 < 0
G_GLNt2r: +0.59209 LC_gln-L_c -0.59209 LC_gln-L_e - DG_GLNt2r + DGo_GLNt2r = 0
FU_GLNt2r: DG_GLNt2r +1000 FU_GLNt2r < 1000
BU_GLNt2r: - DG_GLNt2r +1000 BU_GLNt2r < 1000
SU_GLNt2r: FU_GLNt2r + BU_GLNt2r < 1
UF_GLNt2r: F_GLNt2r -1000 FU_GLNt2r < 0
UR_GLNt2r: R_GLNt2r -1000 BU_GLNt2r < 0
G_GLUDC: +0.59209 LC_4abut_c +0.59209 LC_co2_c -0.59209 LC_glu-L_c - DG_GLUDC + DGo_GLUDC = 0
FU_GLUDC: DG_GLUDC +1000 FU_GLUDC < 1000
BU_GLUDC: - DG_GLUDC +1000 BU_GLUDC < 1000
SU_GLUDC: FU_GLUDC + BU_GLUDC < 1
UF_GLUDC: F_GLUDC -1000 FU_GLUDC < 0
UR_GLUDC: R_GLUDC -1000 BU_GLUDC < 0
G_GLUDxi: +0.59209 LC_akg_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nh4_c - DG_GLUDxi + DGo_GLUDxi = 0
FU_GLUDxi: DG_GLUDxi +1000 FU_GLUDxi < 1000
BU_GLUDxi: - DG_GLUDxi +1000 BU_GLUDxi < 1000
SU_GLUDxi: FU_GLUDxi + BU_GLUDxi < 1
UF_GLUDxi: F_GLUDxi -1000 FU_GLUDxi < 0
UR_GLUDxi: R_GLUDxi -1000 BU_GLUDxi < 0
G_GLUDyi: -0.59209 LC_akg_c +0.59209 LC_glu-L_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -0.59209 LC_nh4_c - DG_GLUDyi + DGo_GLUDyi = 0
FU_GLUDyi: DG_GLUDyi +1000 FU_GLUDyi < 1000
BU_GLUDyi: - DG_GLUDyi +1000 BU_GLUDyi < 1000
SU_GLUDyi: FU_GLUDyi + BU_GLUDyi < 1
UF_GLUDyi: F_GLUDyi -1000 FU_GLUDyi < 0
UR_GLUDyi: R_GLUDyi -1000 BU_GLUDyi < 0
G_GLUK: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_g6p-B_c -0.59209 LC_glc-D_c - DG_GLUK + DGo_GLUK = 0
FU_GLUK: DG_GLUK +1000 FU_GLUK < 1000
BU_GLUK: - DG_GLUK +1000 BU_GLUK < 1000
SU_GLUK: FU_GLUK + BU_GLUK < 1
UF_GLUK: F_GLUK -1000 FU_GLUK < 0
UR_GLUK: R_GLUK -1000 BU_GLUK < 0
G_GLUt2r: +0.59209 LC_glu-L_c -0.59209 LC_glu-L_e - DG_GLUt2r + DGo_GLUt2r = 0
FU_GLUt2r: DG_GLUt2r +1000 FU_GLUt2r < 1000
BU_GLUt2r: - DG_GLUt2r +1000 BU_GLUt2r < 1000
SU_GLUt2r: FU_GLUt2r + BU_GLUt2r < 1
UF_GLUt2r: F_GLUt2r -1000 FU_GLUt2r < 0
UR_GLUt2r: R_GLUt2r -1000 BU_GLUt2r < 0
G_GLUt7m: -0.59209 LC_glu-L_c +0.59209 LC_glu-L_m - DG_GLUt7m + DGo_GLUt7m = 0
FU_GLUt7m: DG_GLUt7m +1000 FU_GLUt7m < 1000
BU_GLUt7m: - DG_GLUt7m +1000 BU_GLUt7m < 1000
SU_GLUt7m: FU_GLUt7m + BU_GLUt7m < 1
UF_GLUt7m: F_GLUt7m -1000 FU_GLUt7m < 0
UR_GLUt7m: R_GLUt7m -1000 BU_GLUt7m < 0
G_GLYt2r: +0.59209 LC_gly_c -0.59209 LC_gly_e - DG_GLYt2r + DGo_GLYt2r = 0
FU_GLYt2r: DG_GLYt2r +1000 FU_GLYt2r < 1000
BU_GLYt2r: - DG_GLYt2r +1000 BU_GLYt2r < 1000
SU_GLYt2r: FU_GLYt2r + BU_GLYt2r < 1
UF_GLYt2r: F_GLYt2r -1000 FU_GLYt2r < 0
UR_GLYt2r: R_GLYt2r -1000 BU_GLYt2r < 0
G_GND: -0.59209 LC_6pgc_c +0.59209 LC_co2_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c +0.59209 LC_ru5p-D_c - DG_GND + DGo_GND = 0
FU_GND: DG_GND +1000 FU_GND < 1000
BU_GND: - DG_GND +1000 BU_GND < 1000
SU_GND: FU_GND + BU_GND < 1
UF_GND: F_GND -1000 FU_GND < 0
UR_GND: R_GND -1000 BU_GND < 0
SU_H2Ot: FU_H2Ot + BU_H2Ot < 1
UF_H2Ot: F_H2Ot -1000 FU_H2Ot < 0
UR_H2Ot: R_H2Ot -1000 BU_H2Ot < 0
SU_H2Oter: FU_H2Oter + BU_H2Oter < 1
UF_H2Oter: F_H2Oter -1000 FU_H2Oter < 0
UR_H2Oter: R_H2Oter -1000 BU_H2Oter < 0
SU_H2Otm: FU_H2Otm + BU_H2Otm < 1
UF_H2Otm: F_H2Otm -1000 FU_H2Otm < 0
UR_H2Otm: R_H2Otm -1000 BU_H2Otm < 0
G_HCITSm: -0.59209 LC_accoa_m -0.59209 LC_akg_m +0.59209 LC_coa_m +0.59209 LC_hcit_m - DG_HCITSm + DGo_HCITSm = 0
FU_HCITSm: DG_HCITSm +1000 FU_HCITSm < 1000
BU_HCITSm: - DG_HCITSm +1000 BU_HCITSm < 1000
SU_HCITSm: FU_HCITSm + BU_HCITSm < 1
UF_HCITSm: F_HCITSm -1000 FU_HCITSm < 0
UR_HCITSm: R_HCITSm -1000 BU_HCITSm < 0
G_HCO3E: -0.59209 LC_co2_c +0.59209 LC_hco3_c - DG_HCO3E + DGo_HCO3E = 0
FU_HCO3E: DG_HCO3E +1000 FU_HCO3E < 1000
BU_HCO3E: - DG_HCO3E +1000 BU_HCO3E < 1000
SU_HCO3E: FU_HCO3E + BU_HCO3E < 1
UF_HCO3E: F_HCO3E -1000 FU_HCO3E < 0
UR_HCO3E: R_HCO3E -1000 BU_HCO3E < 0
G_HEX1: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_g6p_c -0.59209 LC_glc-D_c - DG_HEX1 + DGo_HEX1 = 0
FU_HEX1: DG_HEX1 +1000 FU_HEX1 < 1000
BU_HEX1: - DG_HEX1 +1000 BU_HEX1 < 1000
SU_HEX1: FU_HEX1 + BU_HEX1 < 1
UF_HEX1: F_HEX1 -1000 FU_HEX1 < 0
UR_HEX1: R_HEX1 -1000 BU_HEX1 < 0
G_HISt2r: +0.59209 LC_his-L_c -0.59209 LC_his-L_e - DG_HISt2r + DGo_HISt2r = 0
FU_HISt2r: DG_HISt2r +1000 FU_HISt2r < 1000
BU_HISt2r: - DG_HISt2r +1000 BU_HISt2r < 1000
SU_HISt2r: FU_HISt2r + BU_HISt2r < 1
UF_HISt2r: F_HISt2r -1000 FU_HISt2r < 0
UR_HISt2r: R_HISt2r -1000 BU_HISt2r < 0
G_ICDHxm: +0.59209 LC_akg_m +0.59209 LC_co2_m -0.59209 LC_icit_m -0.59209 LC_nad_m +0.59209 LC_nadh_m - DG_ICDHxm + DGo_ICDHxm = 0
FU_ICDHxm: DG_ICDHxm +1000 FU_ICDHxm < 1000
BU_ICDHxm: - DG_ICDHxm +1000 BU_ICDHxm < 1000
SU_ICDHxm: FU_ICDHxm + BU_ICDHxm < 1
UF_ICDHxm: F_ICDHxm -1000 FU_ICDHxm < 0
UR_ICDHxm: R_ICDHxm -1000 BU_ICDHxm < 0
G_ICDHy: +0.59209 LC_akg_c +0.59209 LC_co2_c -0.59209 LC_icit_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_ICDHy + DGo_ICDHy = 0
FU_ICDHy: DG_ICDHy +1000 FU_ICDHy < 1000
BU_ICDHy: - DG_ICDHy +1000 BU_ICDHy < 1000
SU_ICDHy: FU_ICDHy + BU_ICDHy < 1
UF_ICDHy: F_ICDHy -1000 FU_ICDHy < 0
UR_ICDHy: R_ICDHy -1000 BU_ICDHy < 0
G_ICDHym: +0.59209 LC_akg_m +0.59209 LC_co2_m -0.59209 LC_icit_m -0.59209 LC_nadp_m +0.59209 LC_nadph_m - DG_ICDHym + DGo_ICDHym = 0
FU_ICDHym: DG_ICDHym +1000 FU_ICDHym < 1000
BU_ICDHym: - DG_ICDHym +1000 BU_ICDHym < 1000
SU_ICDHym: FU_ICDHym + BU_ICDHym < 1
UF_ICDHym: F_ICDHym -1000 FU_ICDHym < 0
UR_ICDHym: R_ICDHym -1000 BU_ICDHym < 0
G_ILEt2r: +0.59209 LC_ile-L_c -0.59209 LC_ile-L_e - DG_ILEt2r + DGo_ILEt2r = 0
FU_ILEt2r: DG_ILEt2r +1000 FU_ILEt2r < 1000
BU_ILEt2r: - DG_ILEt2r +1000 BU_ILEt2r < 1000
SU_ILEt2r: FU_ILEt2r + BU_ILEt2r < 1
UF_ILEt2r: F_ILEt2r -1000 FU_ILEt2r < 0
UR_ILEt2r: R_ILEt2r -1000 BU_ILEt2r < 0
G_ILEtmi: +0.59209 LC_ile-L_c -0.59209 LC_ile-L_m - DG_ILEtmi + DGo_ILEtmi = 0
FU_ILEtmi: DG_ILEtmi +1000 FU_ILEtmi < 1000
BU_ILEtmi: - DG_ILEtmi +1000 BU_ILEtmi < 1000
SU_ILEtmi: FU_ILEtmi + BU_ILEtmi < 1
UF_ILEtmi: F_ILEtmi -1000 FU_ILEtmi < 0
UR_ILEtmi: R_ILEtmi -1000 BU_ILEtmi < 0
G_LEUt2r: +0.59209 LC_leu-L_c -0.59209 LC_leu-L_e - DG_LEUt2r + DGo_LEUt2r = 0
FU_LEUt2r: DG_LEUt2r +1000 FU_LEUt2r < 1000
BU_LEUt2r: - DG_LEUt2r +1000 BU_LEUt2r < 1000
SU_LEUt2r: FU_LEUt2r + BU_LEUt2r < 1
UF_LEUt2r: F_LEUt2r -1000 FU_LEUt2r < 0
UR_LEUt2r: R_LEUt2r -1000 BU_LEUt2r < 0
G_LYSt2r: +0.59209 LC_lys-L_c -0.59209 LC_lys-L_e - DG_LYSt2r + DGo_LYSt2r = 0
FU_LYSt2r: DG_LYSt2r +1000 FU_LYSt2r < 1000
BU_LYSt2r: - DG_LYSt2r +1000 BU_LYSt2r < 1000
SU_LYSt2r: FU_LYSt2r + BU_LYSt2r < 1
UF_LYSt2r: F_LYSt2r -1000 FU_LYSt2r < 0
UR_LYSt2r: R_LYSt2r -1000 BU_LYSt2r < 0
G_MALt2r: +0.59209 LC_mal-L_c -0.59209 LC_mal-L_e - DG_MALt2r + DGo_MALt2r = 0
FU_MALt2r: DG_MALt2r +1000 FU_MALt2r < 1000
BU_MALt2r: - DG_MALt2r +1000 BU_MALt2r < 1000
SU_MALt2r: FU_MALt2r + BU_MALt2r < 1
UF_MALt2r: F_MALt2r -1000 FU_MALt2r < 0
UR_MALt2r: R_MALt2r -1000 BU_MALt2r < 0
G_MALtm: -0.59209 LC_mal-L_c +0.59209 LC_pi_c +0.59209 LC_mal-L_m -0.59209 LC_pi_m - DG_MALtm + DGo_MALtm = 0
FU_MALtm: DG_MALtm +1000 FU_MALtm < 1000
BU_MALtm: - DG_MALtm +1000 BU_MALtm < 1000
SU_MALtm: FU_MALtm + BU_MALtm < 1
UF_MALtm: F_MALtm -1000 FU_MALtm < 0
UR_MALtm: R_MALtm -1000 BU_MALtm < 0
G_MAN6PI: +0.59209 LC_f6p_c -0.59209 LC_man6p_c - DG_MAN6PI + DGo_MAN6PI = 0
FU_MAN6PI: DG_MAN6PI +1000 FU_MAN6PI < 1000
BU_MAN6PI: - DG_MAN6PI +1000 BU_MAN6PI < 1000
SU_MAN6PI: FU_MAN6PI + BU_MAN6PI < 1
UF_MAN6PI: F_MAN6PI -1000 FU_MAN6PI < 0
UR_MAN6PI: R_MAN6PI -1000 BU_MAN6PI < 0
G_MANNANter: -0.59209 LC_mannan_c +0.59209 LC_mannan_r - DG_MANNANter + DGo_MANNANter = 0
FU_MANNANter: DG_MANNANter +1000 FU_MANNANter < 1000
BU_MANNANter: - DG_MANNANter +1000 BU_MANNANter < 1000
SU_MANNANter: FU_MANNANter + BU_MANNANter < 1
UF_MANNANter: F_MANNANter -1000 FU_MANNANter < 0
UR_MANNANter: R_MANNANter -1000 BU_MANNANter < 0
G_MDH: -0.59209 LC_mal-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_oaa_c - DG_MDH + DGo_MDH = 0
FU_MDH: DG_MDH +1000 FU_MDH < 1000
BU_MDH: - DG_MDH +1000 BU_MDH < 1000
SU_MDH: FU_MDH + BU_MDH < 1
UF_MDH: F_MDH -1000 FU_MDH < 0
UR_MDH: R_MDH -1000 BU_MDH < 0
G_MDHm: -0.59209 LC_mal-L_m -0.59209 LC_nad_m +0.59209 LC_nadh_m +0.59209 LC_oaa_m - DG_MDHm + DGo_MDHm = 0
FU_MDHm: DG_MDHm +1000 FU_MDHm < 1000
BU_MDHm: - DG_MDHm +1000 BU_MDHm < 1000
SU_MDHm: FU_MDHm + BU_MDHm < 1
UF_MDHm: F_MDHm -1000 FU_MDHm < 0
UR_MDHm: R_MDHm -1000 BU_MDHm < 0
G_ME1m: +0.59209 LC_co2_m -0.59209 LC_mal-L_m -0.59209 LC_nad_m +0.59209 LC_nadh_m +0.59209 LC_pyr_m - DG_ME1m + DGo_ME1m = 0
FU_ME1m: DG_ME1m +1000 FU_ME1m < 1000
BU_ME1m: - DG_ME1m +1000 BU_ME1m < 1000
SU_ME1m: FU_ME1m + BU_ME1m < 1
UF_ME1m: F_ME1m -1000 FU_ME1m < 0
UR_ME1m: R_ME1m -1000 BU_ME1m < 0
G_ME2m: +0.59209 LC_co2_m -0.59209 LC_mal-L_m -0.59209 LC_nadp_m +0.59209 LC_nadph_m +0.59209 LC_pyr_m - DG_ME2m + DGo_ME2m = 0
FU_ME2m: DG_ME2m +1000 FU_ME2m < 1000
BU_ME2m: - DG_ME2m +1000 BU_ME2m < 1000
SU_ME2m: FU_ME2m + BU_ME2m < 1
UF_ME2m: F_ME2m -1000 FU_ME2m < 0
UR_ME2m: R_ME2m -1000 BU_ME2m < 0
G_METt2r: +0.59209 LC_met-L_c -0.59209 LC_met-L_e - DG_METt2r + DGo_METt2r = 0
FU_METt2r: DG_METt2r +1000 FU_METt2r < 1000
BU_METt2r: - DG_METt2r +1000 BU_METt2r < 1000
SU_METt2r: FU_METt2r + BU_METt2r < 1
UF_METt2r: F_METt2r -1000 FU_METt2r < 0
UR_METt2r: R_METt2r -1000 BU_METt2r < 0
G_NADH2-u6cm: +0.59209 LC_nad_c -0.59209 LC_nadh_c -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_NADH2-u6cm + DGo_NADH2-u6cm = 0
FU_NADH2-u6cm: DG_NADH2-u6cm +1000 FU_NADH2-u6cm < 1000
BU_NADH2-u6cm: - DG_NADH2-u6cm +1000 BU_NADH2-u6cm < 1000
SU_NADH2-u6cm: FU_NADH2-u6cm + BU_NADH2-u6cm < 1
UF_NADH2-u6cm: F_NADH2-u6cm -1000 FU_NADH2-u6cm < 0
UR_NADH2-u6cm: R_NADH2-u6cm -1000 BU_NADH2-u6cm < 0
G_NADH2-u6m: +0.59209 LC_nad_m -0.59209 LC_nadh_m -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_NADH2-u6m + DGo_NADH2-u6m = 0
FU_NADH2-u6m: DG_NADH2-u6m +1000 FU_NADH2-u6m < 1000
BU_NADH2-u6m: - DG_NADH2-u6m +1000 BU_NADH2-u6m < 1000
SU_NADH2-u6m: FU_NADH2-u6m + BU_NADH2-u6m < 1
UF_NADH2-u6m: F_NADH2-u6m -1000 FU_NADH2-u6m < 0
UR_NADH2-u6m: R_NADH2-u6m -1000 BU_NADH2-u6m < 0
G_NADK: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_nad_c +0.59209 LC_nadp_c - DG_NADK + DGo_NADK = 0
FU_NADK: DG_NADK +1000 FU_NADK < 1000
BU_NADK: - DG_NADK +1000 BU_NADK < 1000
SU_NADK: FU_NADK + BU_NADK < 1
UF_NADK: F_NADK -1000 FU_NADK < 0
UR_NADK: R_NADK -1000 BU_NADK < 0
G_NADKm: +0.59209 LC_adp_m -0.59209 LC_atp_m -0.59209 LC_nad_m +0.59209 LC_nadp_m - DG_NADKm + DGo_NADKm = 0
FU_NADKm: DG_NADKm +1000 FU_NADKm < 1000
BU_NADKm: - DG_NADKm +1000 BU_NADKm < 1000
SU_NADKm: FU_NADKm + BU_NADKm < 1
UF_NADKm: F_NADKm -1000 FU_NADKm < 0
UR_NADKm: R_NADKm -1000 BU_NADKm < 0
G_NADPPPS: +0.59209 LC_nad_c -0.59209 LC_nadp_c +0.59209 LC_pi_c - DG_NADPPPS + DGo_NADPPPS = 0
FU_NADPPPS: DG_NADPPPS +1000 FU_NADPPPS < 1000
BU_NADPPPS: - DG_NADPPPS +1000 BU_NADPPPS < 1000
SU_NADPPPS: FU_NADPPPS + BU_NADPPPS < 1
UF_NADPPPS: F_NADPPPS -1000 FU_NADPPPS < 0
UR_NADPPPS: R_NADPPPS -1000 BU_NADPPPS < 0
G_NADPPPSm: +0.59209 LC_nad_m -0.59209 LC_nadp_m +0.59209 LC_pi_m - DG_NADPPPSm + DGo_NADPPPSm = 0
FU_NADPPPSm: DG_NADPPPSm +1000 FU_NADPPPSm < 1000
BU_NADPPPSm: - DG_NADPPPSm +1000 BU_NADPPPSm < 1000
SU_NADPPPSm: FU_NADPPPSm + BU_NADPPPSm < 1
UF_NADPPPSm: F_NADPPPSm -1000 FU_NADPPPSm < 0
UR_NADPPPSm: R_NADPPPSm -1000 BU_NADPPPSm < 0
G_NDP1: -0.59209 LC_adp_c +0.59209 LC_amp_c +0.59209 LC_pi_c - DG_NDP1 + DGo_NDP1 = 0
FU_NDP1: DG_NDP1 +1000 FU_NDP1 < 1000
BU_NDP1: - DG_NDP1 +1000 BU_NDP1 < 1000
SU_NDP1: FU_NDP1 + BU_NDP1 < 1
UF_NDP1: F_NDP1 -1000 FU_NDP1 < 0
UR_NDP1: R_NDP1 -1000 BU_NDP1 < 0
G_NH4t: +0.59209 LC_nh4_c -0.59209 LC_nh4_e - DG_NH4t + DGo_NH4t = 0
FU_NH4t: DG_NH4t +1000 FU_NH4t < 1000
BU_NH4t: - DG_NH4t +1000 BU_NH4t < 1000
SU_NH4t: FU_NH4t + BU_NH4t < 1
UF_NH4t: F_NH4t -1000 FU_NH4t < 0
UR_NH4t: R_NH4t -1000 BU_NH4t < 0
G_NH4tm: -0.59209 LC_nh4_c +0.59209 LC_nh4_m - DG_NH4tm + DGo_NH4tm = 0
FU_NH4tm: DG_NH4tm +1000 FU_NH4tm < 1000
BU_NH4tm: - DG_NH4tm +1000 BU_NH4tm < 1000
SU_NH4tm: FU_NH4tm + BU_NH4tm < 1
UF_NH4tm: F_NH4tm -1000 FU_NH4tm < 0
UR_NH4tm: R_NH4tm -1000 BU_NH4tm < 0
G_O2t: +0.59209 LC_o2_c -0.59209 LC_o2_e - DG_O2t + DGo_O2t = 0
FU_O2t: DG_O2t +1000 FU_O2t < 1000
BU_O2t: - DG_O2t +1000 BU_O2t < 1000
SU_O2t: FU_O2t + BU_O2t < 1
UF_O2t: F_O2t -1000 FU_O2t < 0
UR_O2t: R_O2t -1000 BU_O2t < 0
G_O2ter: -0.59209 LC_o2_c +0.59209 LC_o2_r - DG_O2ter + DGo_O2ter = 0
FU_O2ter: DG_O2ter +1000 FU_O2ter < 1000
BU_O2ter: - DG_O2ter +1000 BU_O2ter < 1000
SU_O2ter: FU_O2ter + BU_O2ter < 1
UF_O2ter: F_O2ter -1000 FU_O2ter < 0
UR_O2ter: R_O2ter -1000 BU_O2ter < 0
G_O2tm: -0.59209 LC_o2_c +0.59209 LC_o2_m - DG_O2tm + DGo_O2tm = 0
FU_O2tm: DG_O2tm +1000 FU_O2tm < 1000
BU_O2tm: - DG_O2tm +1000 BU_O2tm < 1000
SU_O2tm: FU_O2tm + BU_O2tm < 1
UF_O2tm: F_O2tm -1000 FU_O2tm < 0
UR_O2tm: R_O2tm -1000 BU_O2tm < 0
G_OAAt: -0.59209 LC_oaa_c +0.59209 LC_oaa_e - DG_OAAt + DGo_OAAt = 0
FU_OAAt: DG_OAAt +1000 FU_OAAt < 1000
BU_OAAt: - DG_OAAt +1000 BU_OAAt < 1000
SU_OAAt: FU_OAAt + BU_OAAt < 1
UF_OAAt: F_OAAt -1000 FU_OAAt < 0
UR_OAAt: R_OAAt -1000 BU_OAAt < 0
G_OAAt2m: -0.59209 LC_oaa_c +0.59209 LC_oaa_m - DG_OAAt2m + DGo_OAAt2m = 0
FU_OAAt2m: DG_OAAt2m +1000 FU_OAAt2m < 1000
BU_OAAt2m: - DG_OAAt2m +1000 BU_OAAt2m < 1000
SU_OAAt2m: FU_OAAt2m + BU_OAAt2m < 1
UF_OAAt2m: F_OAAt2m -1000 FU_OAAt2m < 0
UR_OAAt2m: R_OAAt2m -1000 BU_OAAt2m < 0
G_PC: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_hco3_c +0.59209 LC_oaa_c +0.59209 LC_pi_c -0.59209 LC_pyr_c - DG_PC + DGo_PC = 0
FU_PC: DG_PC +1000 FU_PC < 1000
BU_PC: - DG_PC +1000 BU_PC < 1000
SU_PC: FU_PC + BU_PC < 1
UF_PC: F_PC -1000 FU_PC < 0
UR_PC: R_PC -1000 BU_PC < 0
G_PDHcm: -0.59209 LC_dhlam_m +0.59209 LC_lpam_m -0.59209 LC_nad_m +0.59209 LC_nadh_m - DG_PDHcm + DGo_PDHcm = 0
FU_PDHcm: DG_PDHcm +1000 FU_PDHcm < 1000
BU_PDHcm: - DG_PDHcm +1000 BU_PDHcm < 1000
SU_PDHcm: FU_PDHcm + BU_PDHcm < 1
UF_PDHcm: F_PDHcm -1000 FU_PDHcm < 0
UR_PDHcm: R_PDHcm -1000 BU_PDHcm < 0
G_PDHm: +0.59209 LC_accoa_m +0.59209 LC_co2_m -0.59209 LC_coa_m -0.59209 LC_nad_m +0.59209 LC_nadh_m -0.59209 LC_pyr_m - DG_PDHm + DGo_PDHm = 0
FU_PDHm: DG_PDHm +1000 FU_PDHm < 1000
BU_PDHm: - DG_PDHm +1000 BU_PDHm < 1000
SU_PDHm: FU_PDHm + BU_PDHm < 1
UF_PDHm: F_PDHm -1000 FU_PDHm < 0
UR_PDHm: R_PDHm -1000 BU_PDHm < 0
SU_PEtm_SC: FU_PEtm_SC + BU_PEtm_SC < 1
UF_PEtm_SC: F_PEtm_SC -1000 FU_PEtm_SC < 0
UR_PEtm_SC: R_PEtm_SC -1000 BU_PEtm_SC < 0
G_PFK: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_f6p_c +0.59209 LC_fdp_c - DG_PFK + DGo_PFK = 0
FU_PFK: DG_PFK +1000 FU_PFK < 1000
BU_PFK: - DG_PFK +1000 BU_PFK < 1000
SU_PFK: FU_PFK + BU_PFK < 1
UF_PFK: F_PFK -1000 FU_PFK < 0
UR_PFK: R_PFK -1000 BU_PFK < 0
G_PFK_3: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_s17bp_c -0.59209 LC_s7p_c - DG_PFK_3 + DGo_PFK_3 = 0
FU_PFK_3: DG_PFK_3 +1000 FU_PFK_3 < 1000
BU_PFK_3: - DG_PFK_3 +1000 BU_PFK_3 < 1000
SU_PFK_3: FU_PFK_3 + BU_PFK_3 < 1
UF_PFK_3: F_PFK_3 -1000 FU_PFK_3 < 0
UR_PFK_3: R_PFK_3 -1000 BU_PFK_3 < 0
G_PGI: +0.59209 LC_f6p_c -0.59209 LC_g6p_c - DG_PGI + DGo_PGI = 0
FU_PGI: DG_PGI +1000 FU_PGI < 1000
BU_PGI: - DG_PGI +1000 BU_PGI < 1000
SU_PGI: FU_PGI + BU_PGI < 1
UF_PGI: F_PGI -1000 FU_PGI < 0
UR_PGI: R_PGI -1000 BU_PGI < 0
G_PGK: +0.59209 LC_13dpg_c -0.59209 LC_3pg_c +0.59209 LC_adp_c -0.59209 LC_atp_c - DG_PGK + DGo_PGK = 0
FU_PGK: DG_PGK +1000 FU_PGK < 1000
BU_PGK: - DG_PGK +1000 BU_PGK < 1000
SU_PGK: FU_PGK + BU_PGK < 1
UF_PGK: F_PGK -1000 FU_PGK < 0
UR_PGK: R_PGK -1000 BU_PGK < 0
G_PGL: +0.59209 LC_6pgc_c -0.59209 LC_6pgl_c - DG_PGL + DGo_PGL = 0
FU_PGL: DG_PGL +1000 FU_PGL < 1000
BU_PGL: - DG_PGL +1000 BU_PGL < 1000
SU_PGL: FU_PGL + BU_PGL < 1
UF_PGL: F_PGL -1000 FU_PGL < 0
UR_PGL: R_PGL -1000 BU_PGL < 0
G_PGM: -0.59209 LC_2pg_c +0.59209 LC_3pg_c - DG_PGM + DGo_PGM = 0
FU_PGM: DG_PGM +1000 FU_PGM < 1000
BU_PGM: - DG_PGM +1000 BU_PGM < 1000
SU_PGM: FU_PGM + BU_PGM < 1
UF_PGM: F_PGM -1000 FU_PGM < 0
UR_PGM: R_PGM -1000 BU_PGM < 0
G_PGMT: -0.59209 LC_g1p_c +0.59209 LC_g6p_c - DG_PGMT + DGo_PGMT = 0
FU_PGMT: DG_PGMT +1000 FU_PGMT < 1000
BU_PGMT: - DG_PGMT +1000 BU_PGMT < 1000
SU_PGMT: FU_PGMT + BU_PGMT < 1
UF_PGMT: F_PGMT -1000 FU_PGMT < 0
UR_PGMT: R_PGMT -1000 BU_PGMT < 0
G_PHEt2r: +0.59209 LC_phe-L_c -0.59209 LC_phe-L_e - DG_PHEt2r + DGo_PHEt2r = 0
FU_PHEt2r: DG_PHEt2r +1000 FU_PHEt2r < 1000
BU_PHEt2r: - DG_PHEt2r +1000 BU_PHEt2r < 1000
SU_PHEt2r: FU_PHEt2r + BU_PHEt2r < 1
UF_PHEt2r: F_PHEt2r -1000 FU_PHEt2r < 0
UR_PHEt2r: R_PHEt2r -1000 BU_PHEt2r < 0
G_PIt2m: -0.59209 LC_pi_c +0.59209 LC_pi_m - DG_PIt2m + DGo_PIt2m = 0
FU_PIt2m: DG_PIt2m +1000 FU_PIt2m < 1000
BU_PIt2m: - DG_PIt2m +1000 BU_PIt2m < 1000
SU_PIt2m: FU_PIt2m + BU_PIt2m < 1
UF_PIt2m: F_PIt2m -1000 FU_PIt2m < 0
UR_PIt2m: R_PIt2m -1000 BU_PIt2m < 0
G_PIt2r: +0.59209 LC_pi_c -0.59209 LC_pi_e - DG_PIt2r + DGo_PIt2r = 0
FU_PIt2r: DG_PIt2r +1000 FU_PIt2r < 1000
BU_PIt2r: - DG_PIt2r +1000 BU_PIt2r < 1000
SU_PIt2r: FU_PIt2r + BU_PIt2r < 1
UF_PIt2r: F_PIt2r -1000 FU_PIt2r < 0
UR_PIt2r: R_PIt2r -1000 BU_PIt2r < 0
G_PPA: +1.1842 LC_pi_c -0.59209 LC_ppi_c - DG_PPA + DGo_PPA = 0
FU_PPA: DG_PPA +1000 FU_PPA < 1000
BU_PPA: - DG_PPA +1000 BU_PPA < 1000
SU_PPA: FU_PPA + BU_PPA < 1
UF_PPA: F_PPA -1000 FU_PPA < 0
UR_PPA: R_PPA -1000 BU_PPA < 0
G_PPAm: +1.1842 LC_pi_m -0.59209 LC_ppi_m - DG_PPAm + DGo_PPAm = 0
FU_PPAm: DG_PPAm +1000 FU_PPAm < 1000
BU_PPAm: - DG_PPAm +1000 BU_PPAm < 1000
SU_PPAm: FU_PPAm + BU_PPAm < 1
UF_PPAm: F_PPAm -1000 FU_PPAm < 0
UR_PPAm: R_PPAm -1000 BU_PPAm < 0
G_PPCK: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_oaa_c +0.59209 LC_pep_c - DG_PPCK + DGo_PPCK = 0
FU_PPCK: DG_PPCK +1000 FU_PPCK < 1000
BU_PPCK: - DG_PPCK +1000 BU_PPCK < 1000
SU_PPCK: FU_PPCK + BU_PPCK < 1
UF_PPCK: F_PPCK -1000 FU_PPCK < 0
UR_PPCK: R_PPCK -1000 BU_PPCK < 0
G_PROt2r: +0.59209 LC_pro-L_c -0.59209 LC_pro-L_e - DG_PROt2r + DGo_PROt2r = 0
FU_PROt2r: DG_PROt2r +1000 FU_PROt2r < 1000
BU_PROt2r: - DG_PROt2r +1000 BU_PROt2r < 1000
SU_PROt2r: FU_PROt2r + BU_PROt2r < 1
UF_PROt2r: F_PROt2r -1000 FU_PROt2r < 0
UR_PROt2r: R_PROt2r -1000 BU_PROt2r < 0
SU_PStm_SC: FU_PStm_SC + BU_PStm_SC < 1
UF_PStm_SC: F_PStm_SC -1000 FU_PStm_SC < 0
UR_PStm_SC: R_PStm_SC -1000 BU_PStm_SC < 0
G_PYK: -0.59209 LC_adp_c +0.59209 LC_atp_c -0.59209 LC_pep_c +0.59209 LC_pyr_c - DG_PYK + DGo_PYK = 0
FU_PYK: DG_PYK +1000 FU_PYK < 1000
BU_PYK: - DG_PYK +1000 BU_PYK < 1000
SU_PYK: FU_PYK + BU_PYK < 1
UF_PYK: F_PYK -1000 FU_PYK < 0
UR_PYK: R_PYK -1000 BU_PYK < 0
G_PYRDC: +0.59209 LC_acald_c +0.59209 LC_co2_c -0.59209 LC_pyr_c - DG_PYRDC + DGo_PYRDC = 0
FU_PYRDC: DG_PYRDC +1000 FU_PYRDC < 1000
BU_PYRDC: - DG_PYRDC +1000 BU_PYRDC < 1000
SU_PYRDC: FU_PYRDC + BU_PYRDC < 1
UF_PYRDC: F_PYRDC -1000 FU_PYRDC < 0
UR_PYRDC: R_PYRDC -1000 BU_PYRDC < 0
G_PYRt: -0.59209 LC_pyr_c +0.59209 LC_pyr_e - DG_PYRt + DGo_PYRt = 0
FU_PYRt: DG_PYRt +1000 FU_PYRt < 1000
BU_PYRt: - DG_PYRt +1000 BU_PYRt < 1000
SU_PYRt: FU_PYRt + BU_PYRt < 1
UF_PYRt: F_PYRt -1000 FU_PYRt < 0
UR_PYRt: R_PYRt -1000 BU_PYRt < 0
G_PYRt2m: -0.59209 LC_pyr_c +0.59209 LC_pyr_m - DG_PYRt2m + DGo_PYRt2m = 0
FU_PYRt2m: DG_PYRt2m +1000 FU_PYRt2m < 1000
BU_PYRt2m: - DG_PYRt2m +1000 BU_PYRt2m < 1000
SU_PYRt2m: FU_PYRt2m + BU_PYRt2m < 1
UF_PYRt2m: F_PYRt2m -1000 FU_PYRt2m < 0
UR_PYRt2m: R_PYRt2m -1000 BU_PYRt2m < 0
G_RPE: -0.59209 LC_ru5p-D_c +0.59209 LC_xu5p-D_c - DG_RPE + DGo_RPE = 0
FU_RPE: DG_RPE +1000 FU_RPE < 1000
BU_RPE: - DG_RPE +1000 BU_RPE < 1000
SU_RPE: FU_RPE + BU_RPE < 1
UF_RPE: F_RPE -1000 FU_RPE < 0
UR_RPE: R_RPE -1000 BU_RPE < 0
G_RPI: -0.59209 LC_r5p_c +0.59209 LC_ru5p-D_c - DG_RPI + DGo_RPI = 0
FU_RPI: DG_RPI +1000 FU_RPI < 1000
BU_RPI: - DG_RPI +1000 BU_RPI < 1000
SU_RPI: FU_RPI + BU_RPI < 1
UF_RPI: F_RPI -1000 FU_RPI < 0
UR_RPI: R_RPI -1000 BU_RPI < 0
G_SERt2r: +0.59209 LC_ser-L_c -0.59209 LC_ser-L_e - DG_SERt2r + DGo_SERt2r = 0
FU_SERt2r: DG_SERt2r +1000 FU_SERt2r < 1000
BU_SERt2r: - DG_SERt2r +1000 BU_SERt2r < 1000
SU_SERt2r: FU_SERt2r + BU_SERt2r < 1
UF_SERt2r: F_SERt2r -1000 FU_SERt2r < 0
UR_SERt2r: R_SERt2r -1000 BU_SERt2r < 0
G_SO4ti: +0.59209 LC_so4_c -0.59209 LC_so4_e - DG_SO4ti + DGo_SO4ti = 0
FU_SO4ti: DG_SO4ti +1000 FU_SO4ti < 1000
BU_SO4ti: - DG_SO4ti +1000 BU_SO4ti < 1000
SU_SO4ti: FU_SO4ti + BU_SO4ti < 1
UF_SO4ti: F_SO4ti -1000 FU_SO4ti < 0
UR_SO4ti: R_SO4ti -1000 BU_SO4ti < 0
G_SUCCt2r: +0.59209 LC_succ_c -0.59209 LC_succ_e - DG_SUCCt2r + DGo_SUCCt2r = 0
FU_SUCCt2r: DG_SUCCt2r +1000 FU_SUCCt2r < 1000
BU_SUCCt2r: - DG_SUCCt2r +1000 BU_SUCCt2r < 1000
SU_SUCCt2r: FU_SUCCt2r + BU_SUCCt2r < 1
UF_SUCCt2r: F_SUCCt2r -1000 FU_SUCCt2r < 0
UR_SUCCt2r: R_SUCCt2r -1000 BU_SUCCt2r < 0
G_SUCCtm: +0.59209 LC_pi_c -0.59209 LC_succ_c -0.59209 LC_pi_m +0.59209 LC_succ_m - DG_SUCCtm + DGo_SUCCtm = 0
FU_SUCCtm: DG_SUCCtm +1000 FU_SUCCtm < 1000
BU_SUCCtm: - DG_SUCCtm +1000 BU_SUCCtm < 1000
SU_SUCCtm: FU_SUCCtm + BU_SUCCtm < 1
UF_SUCCtm: F_SUCCtm -1000 FU_SUCCtm < 0
UR_SUCCtm: R_SUCCtm -1000 BU_SUCCtm < 0
G_SUCD1m: -0.59209 LC_fad_m +0.59209 LC_fadh2_m +0.59209 LC_fum_m -0.59209 LC_succ_m - DG_SUCD1m + DGo_SUCD1m = 0
FU_SUCD1m: DG_SUCD1m +1000 FU_SUCD1m < 1000
BU_SUCD1m: - DG_SUCD1m +1000 BU_SUCD1m < 1000
SU_SUCD1m: FU_SUCD1m + BU_SUCD1m < 1
UF_SUCD1m: F_SUCD1m -1000 FU_SUCD1m < 0
UR_SUCD1m: R_SUCD1m -1000 BU_SUCD1m < 0
G_SUCD2_u6m: +0.59209 LC_fum_m -0.59209 LC_q6_m +0.59209 LC_q6h2_m -0.59209 LC_succ_m - DG_SUCD2_u6m + DGo_SUCD2_u6m = 0
FU_SUCD2_u6m: DG_SUCD2_u6m +1000 FU_SUCD2_u6m < 1000
BU_SUCD2_u6m: - DG_SUCD2_u6m +1000 BU_SUCD2_u6m < 1000
SU_SUCD2_u6m: FU_SUCD2_u6m + BU_SUCD2_u6m < 1
UF_SUCD2_u6m: F_SUCD2_u6m -1000 FU_SUCD2_u6m < 0
UR_SUCD2_u6m: R_SUCD2_u6m -1000 BU_SUCD2_u6m < 0
G_SUCD3-u6m: +0.59209 LC_fad_m -0.59209 LC_fadh2_m -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_SUCD3-u6m + DGo_SUCD3-u6m = 0
FU_SUCD3-u6m: DG_SUCD3-u6m +1000 FU_SUCD3-u6m < 1000
BU_SUCD3-u6m: - DG_SUCD3-u6m +1000 BU_SUCD3-u6m < 1000
SU_SUCD3-u6m: FU_SUCD3-u6m + BU_SUCD3-u6m < 1
UF_SUCD3-u6m: F_SUCD3-u6m -1000 FU_SUCD3-u6m < 0
UR_SUCD3-u6m: R_SUCD3-u6m -1000 BU_SUCD3-u6m < 0
G_SUCFUMtm: +0.59209 LC_fum_c -0.59209 LC_succ_c -0.59209 LC_fum_m +0.59209 LC_succ_m - DG_SUCFUMtm + DGo_SUCFUMtm = 0
FU_SUCFUMtm: DG_SUCFUMtm +1000 FU_SUCFUMtm < 1000
BU_SUCFUMtm: - DG_SUCFUMtm +1000 BU_SUCFUMtm < 1000
SU_SUCFUMtm: FU_SUCFUMtm + BU_SUCFUMtm < 1
UF_SUCFUMtm: F_SUCFUMtm -1000 FU_SUCFUMtm < 0
UR_SUCFUMtm: R_SUCFUMtm -1000 BU_SUCFUMtm < 0
G_SUCOASm: +0.59209 LC_adp_m -0.59209 LC_atp_m -0.59209 LC_coa_m +0.59209 LC_pi_m -0.59209 LC_succ_m +0.59209 LC_succoa_m - DG_SUCOASm + DGo_SUCOASm = 0
FU_SUCOASm: DG_SUCOASm +1000 FU_SUCOASm < 1000
BU_SUCOASm: - DG_SUCOASm +1000 BU_SUCOASm < 1000
SU_SUCOASm: FU_SUCOASm + BU_SUCOASm < 1
UF_SUCOASm: F_SUCOASm -1000 FU_SUCOASm < 0
UR_SUCOASm: R_SUCOASm -1000 BU_SUCOASm < 0
G_TALA: +0.59209 LC_e4p_c +0.59209 LC_f6p_c -0.59209 LC_g3p_c -0.59209 LC_s7p_c - DG_TALA + DGo_TALA = 0
FU_TALA: DG_TALA +1000 FU_TALA < 1000
BU_TALA: - DG_TALA +1000 BU_TALA < 1000
SU_TALA: FU_TALA + BU_TALA < 1
UF_TALA: F_TALA -1000 FU_TALA < 0
UR_TALA: R_TALA -1000 BU_TALA < 0
G_THRt2m: -0.59209 LC_thr-L_c +0.59209 LC_thr-L_m - DG_THRt2m + DGo_THRt2m = 0
FU_THRt2m: DG_THRt2m +1000 FU_THRt2m < 1000
BU_THRt2m: - DG_THRt2m +1000 BU_THRt2m < 1000
SU_THRt2m: FU_THRt2m + BU_THRt2m < 1
UF_THRt2m: F_THRt2m -1000 FU_THRt2m < 0
UR_THRt2m: R_THRt2m -1000 BU_THRt2m < 0
G_THRt2r: +0.59209 LC_thr-L_c -0.59209 LC_thr-L_e - DG_THRt2r + DGo_THRt2r = 0
FU_THRt2r: DG_THRt2r +1000 FU_THRt2r < 1000
BU_THRt2r: - DG_THRt2r +1000 BU_THRt2r < 1000
SU_THRt2r: FU_THRt2r + BU_THRt2r < 1
UF_THRt2r: F_THRt2r -1000 FU_THRt2r < 0
UR_THRt2r: R_THRt2r -1000 BU_THRt2r < 0
G_TKT1: +0.59209 LC_g3p_c -0.59209 LC_r5p_c +0.59209 LC_s7p_c -0.59209 LC_xu5p-D_c - DG_TKT1 + DGo_TKT1 = 0
FU_TKT1: DG_TKT1 +1000 FU_TKT1 < 1000
BU_TKT1: - DG_TKT1 +1000 BU_TKT1 < 1000
SU_TKT1: FU_TKT1 + BU_TKT1 < 1
UF_TKT1: F_TKT1 -1000 FU_TKT1 < 0
UR_TKT1: R_TKT1 -1000 BU_TKT1 < 0
G_TKT2: -0.59209 LC_e4p_c +0.59209 LC_f6p_c +0.59209 LC_g3p_c -0.59209 LC_xu5p-D_c - DG_TKT2 + DGo_TKT2 = 0
FU_TKT2: DG_TKT2 +1000 FU_TKT2 < 1000
BU_TKT2: - DG_TKT2 +1000 BU_TKT2 < 1000
SU_TKT2: FU_TKT2 + BU_TKT2 < 1
UF_TKT2: F_TKT2 -1000 FU_TKT2 < 0
UR_TKT2: R_TKT2 -1000 BU_TKT2 < 0
G_TPI: -0.59209 LC_dhap_c +0.59209 LC_g3p_c - DG_TPI + DGo_TPI = 0
FU_TPI: DG_TPI +1000 FU_TPI < 1000
BU_TPI: - DG_TPI +1000 BU_TPI < 1000
SU_TPI: FU_TPI + BU_TPI < 1
UF_TPI: F_TPI -1000 FU_TPI < 0
UR_TPI: R_TPI -1000 BU_TPI < 0
G_TRPt2r: +0.59209 LC_trp-L_c -0.59209 LC_trp-L_e - DG_TRPt2r + DGo_TRPt2r = 0
FU_TRPt2r: DG_TRPt2r +1000 FU_TRPt2r < 1000
BU_TRPt2r: - DG_TRPt2r +1000 BU_TRPt2r < 1000
SU_TRPt2r: FU_TRPt2r + BU_TRPt2r < 1
UF_TRPt2r: F_TRPt2r -1000 FU_TRPt2r < 0
UR_TRPt2r: R_TRPt2r -1000 BU_TRPt2r < 0
G_TYRt2r: +0.59209 LC_tyr-L_c -0.59209 LC_tyr-L_e - DG_TYRt2r + DGo_TYRt2r = 0
FU_TYRt2r: DG_TYRt2r +1000 FU_TYRt2r < 1000
BU_TYRt2r: - DG_TYRt2r +1000 BU_TYRt2r < 1000
SU_TYRt2r: FU_TYRt2r + BU_TYRt2r < 1
UF_TYRt2r: F_TYRt2r -1000 FU_TYRt2r < 0
UR_TYRt2r: R_TYRt2r -1000 BU_TYRt2r < 0
G_VALt2r: +0.59209 LC_val-L_c -0.59209 LC_val-L_e - DG_VALt2r + DGo_VALt2r = 0
FU_VALt2r: DG_VALt2r +1000 FU_VALt2r < 1000
BU_VALt2r: - DG_VALt2r +1000 BU_VALt2r < 1000
SU_VALt2r: FU_VALt2r + BU_VALt2r < 1
UF_VALt2r: F_VALt2r -1000 FU_VALt2r < 0
UR_VALt2r: R_VALt2r -1000 BU_VALt2r < 0
G_ZYMSTt: +0.59209 LC_zymst_c -0.59209 LC_zymst_e - DG_ZYMSTt + DGo_ZYMSTt = 0
FU_ZYMSTt: DG_ZYMSTt +1000 FU_ZYMSTt < 1000
BU_ZYMSTt: - DG_ZYMSTt +1000 BU_ZYMSTt < 1000
SU_ZYMSTt: FU_ZYMSTt + BU_ZYMSTt < 1
UF_ZYMSTt: F_ZYMSTt -1000 FU_ZYMSTt < 0
UR_ZYMSTt: R_ZYMSTt -1000 BU_ZYMSTt < 0
SU_LMPD_1_13BDglcn_c: FU_LMPD_1_13BDglcn_c + BU_LMPD_1_13BDglcn_c < 1
UF_LMPD_1_13BDglcn_c: F_LMPD_1_13BDglcn_c -1000 FU_LMPD_1_13BDglcn_c < 0
UR_LMPD_1_13BDglcn_c: R_LMPD_1_13BDglcn_c -1000 BU_LMPD_1_13BDglcn_c < 0
G_LMPD_2_amp_c: -1.1842 LC_3pg_c +4.1446 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -1.1842 LC_asp-L_c -4.7367 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_2_amp_c + DGo_LMPD_2_amp_c = 0
FU_LMPD_2_amp_c: DG_LMPD_2_amp_c +1000 FU_LMPD_2_amp_c < 1000
BU_LMPD_2_amp_c: - DG_LMPD_2_amp_c +1000 BU_LMPD_2_amp_c < 1000
SU_LMPD_2_amp_c: FU_LMPD_2_amp_c + BU_LMPD_2_amp_c < 1
UF_LMPD_2_amp_c: F_LMPD_2_amp_c -1000 FU_LMPD_2_amp_c < 0
UR_LMPD_2_amp_c: R_LMPD_2_amp_c -1000 BU_LMPD_2_amp_c < 0
G_LMPD_3_amp_c: -1.1842 LC_3pg_c +4.1446 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -1.1842 LC_asp-L_c -4.7367 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_3_amp_c + DGo_LMPD_3_amp_c = 0
FU_LMPD_3_amp_c: DG_LMPD_3_amp_c +1000 FU_LMPD_3_amp_c < 1000
BU_LMPD_3_amp_c: - DG_LMPD_3_amp_c +1000 BU_LMPD_3_amp_c < 1000
SU_LMPD_3_amp_c: FU_LMPD_3_amp_c + BU_LMPD_3_amp_c < 1
UF_LMPD_3_amp_c: F_LMPD_3_amp_c -1000 FU_LMPD_3_amp_c < 0
UR_LMPD_3_amp_c: R_LMPD_3_amp_c -1000 BU_LMPD_3_amp_c < 0
G_LMPD_4_amp_c: -1.1842 LC_3pg_c +2.9604 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -1.1842 LC_asp-L_c -4.1446 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_4_amp_c + DGo_LMPD_4_amp_c = 0
FU_LMPD_4_amp_c: DG_LMPD_4_amp_c +1000 FU_LMPD_4_amp_c < 1000
BU_LMPD_4_amp_c: - DG_LMPD_4_amp_c +1000 BU_LMPD_4_amp_c < 1000
SU_LMPD_4_amp_c: FU_LMPD_4_amp_c + BU_LMPD_4_amp_c < 1
UF_LMPD_4_amp_c: F_LMPD_4_amp_c -1000 FU_LMPD_4_amp_c < 0
UR_LMPD_4_amp_c: R_LMPD_4_amp_c -1000 BU_LMPD_4_amp_c < 0
G_LMPD_5_amp_c: -1.1842 LC_3pg_c +2.9604 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -1.1842 LC_asp-L_c -4.1446 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_5_amp_c + DGo_LMPD_5_amp_c = 0
FU_LMPD_5_amp_c: DG_LMPD_5_amp_c +1000 FU_LMPD_5_amp_c < 1000
BU_LMPD_5_amp_c: - DG_LMPD_5_amp_c +1000 BU_LMPD_5_amp_c < 1000
SU_LMPD_5_amp_c: FU_LMPD_5_amp_c + BU_LMPD_5_amp_c < 1
UF_LMPD_5_amp_c: F_LMPD_5_amp_c -1000 FU_LMPD_5_amp_c < 0
UR_LMPD_5_amp_c: R_LMPD_5_amp_c -1000 BU_LMPD_5_amp_c < 0
G_LMPD_6_arg-L_c: +1.7763 LC_adp_c +0.59209 LC_amp_c +0.59209 LC_arg-L_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_fum_c -0.59209 LC_hco3_c -0.59209 LC_nh4_c +1.7763 LC_pi_c +0.59209 LC_ppi_c +0.59209 LC_adp_m +0.59209 LC_akg_m -0.59209 LC_atp_m -1.1842 LC_glu-L_m +0.59209 LC_nadp_m -0.59209 LC_nadph_m +0.59209 LC_pi_m - DG_LMPD_6_arg-L_c + DGo_LMPD_6_arg-L_c = 0
FU_LMPD_6_arg-L_c: DG_LMPD_6_arg-L_c +1000 FU_LMPD_6_arg-L_c < 1000
BU_LMPD_6_arg-L_c: - DG_LMPD_6_arg-L_c +1000 BU_LMPD_6_arg-L_c < 1000
SU_LMPD_6_arg-L_c: FU_LMPD_6_arg-L_c + BU_LMPD_6_arg-L_c < 1
UF_LMPD_6_arg-L_c: F_LMPD_6_arg-L_c -1000 FU_LMPD_6_arg-L_c < 0
UR_LMPD_6_arg-L_c: R_LMPD_6_arg-L_c -1000 BU_LMPD_6_arg-L_c < 0
G_LMPD_7_asn-L_c: +0.59209 LC_adp_c +0.59209 LC_amp_c +0.59209 LC_asn-L_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c -0.59209 LC_nh4_c +0.59209 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_7_asn-L_c + DGo_LMPD_7_asn-L_c = 0
FU_LMPD_7_asn-L_c: DG_LMPD_7_asn-L_c +1000 FU_LMPD_7_asn-L_c < 1000
BU_LMPD_7_asn-L_c: - DG_LMPD_7_asn-L_c +1000 BU_LMPD_7_asn-L_c < 1000
SU_LMPD_7_asn-L_c: FU_LMPD_7_asn-L_c + BU_LMPD_7_asn-L_c < 1
UF_LMPD_7_asn-L_c: F_LMPD_7_asn-L_c -1000 FU_LMPD_7_asn-L_c < 0
UR_LMPD_7_asn-L_c: R_LMPD_7_asn-L_c -1000 BU_LMPD_7_asn-L_c < 0
G_LMPD_8_cmp_c: +2.9604 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -3.5525 LC_atp_c +0.59209 LC_cmp_c +0.59209 LC_co2_c -0.59209 LC_hco3_c -1.1842 LC_nh4_c +2.9604 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_LMPD_8_cmp_c + DGo_LMPD_8_cmp_c = 0
FU_LMPD_8_cmp_c: DG_LMPD_8_cmp_c +1000 FU_LMPD_8_cmp_c < 1000
BU_LMPD_8_cmp_c: - DG_LMPD_8_cmp_c +1000 BU_LMPD_8_cmp_c < 1000
SU_LMPD_8_cmp_c: FU_LMPD_8_cmp_c + BU_LMPD_8_cmp_c < 1
UF_LMPD_8_cmp_c: F_LMPD_8_cmp_c -1000 FU_LMPD_8_cmp_c < 0
UR_LMPD_8_cmp_c: R_LMPD_8_cmp_c -1000 BU_LMPD_8_cmp_c < 0
G_LMPD_9_cmp_c: +2.9604 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -3.5525 LC_atp_c +0.59209 LC_cmp_c +0.59209 LC_co2_c -0.59209 LC_fum_c -0.59209 LC_hco3_c -1.1842 LC_nh4_c +2.9604 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_succ_c - DG_LMPD_9_cmp_c + DGo_LMPD_9_cmp_c = 0
FU_LMPD_9_cmp_c: DG_LMPD_9_cmp_c +1000 FU_LMPD_9_cmp_c < 1000
BU_LMPD_9_cmp_c: - DG_LMPD_9_cmp_c +1000 BU_LMPD_9_cmp_c < 1000
SU_LMPD_9_cmp_c: FU_LMPD_9_cmp_c + BU_LMPD_9_cmp_c < 1
UF_LMPD_9_cmp_c: F_LMPD_9_cmp_c -1000 FU_LMPD_9_cmp_c < 0
UR_LMPD_9_cmp_c: R_LMPD_9_cmp_c -1000 BU_LMPD_9_cmp_c < 0
G_LMPD_10_cmp_c: +2.3684 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.9604 LC_atp_c +0.59209 LC_cmp_c +0.59209 LC_co2_c -0.59209 LC_hco3_c -1.1842 LC_nh4_c +2.3684 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_LMPD_10_cmp_c + DGo_LMPD_10_cmp_c = 0
FU_LMPD_10_cmp_c: DG_LMPD_10_cmp_c +1000 FU_LMPD_10_cmp_c < 1000
BU_LMPD_10_cmp_c: - DG_LMPD_10_cmp_c +1000 BU_LMPD_10_cmp_c < 1000
SU_LMPD_10_cmp_c: FU_LMPD_10_cmp_c + BU_LMPD_10_cmp_c < 1
UF_LMPD_10_cmp_c: F_LMPD_10_cmp_c -1000 FU_LMPD_10_cmp_c < 0
UR_LMPD_10_cmp_c: R_LMPD_10_cmp_c -1000 BU_LMPD_10_cmp_c < 0
G_LMPD_11_cmp_c: +2.3684 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.9604 LC_atp_c +0.59209 LC_cmp_c +0.59209 LC_co2_c -0.59209 LC_fum_c -0.59209 LC_hco3_c -1.1842 LC_nh4_c +2.3684 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_succ_c - DG_LMPD_11_cmp_c + DGo_LMPD_11_cmp_c = 0
FU_LMPD_11_cmp_c: DG_LMPD_11_cmp_c +1000 FU_LMPD_11_cmp_c < 1000
BU_LMPD_11_cmp_c: - DG_LMPD_11_cmp_c +1000 BU_LMPD_11_cmp_c < 1000
SU_LMPD_11_cmp_c: FU_LMPD_11_cmp_c + BU_LMPD_11_cmp_c < 1
UF_LMPD_11_cmp_c: F_LMPD_11_cmp_c -1000 FU_LMPD_11_cmp_c < 0
UR_LMPD_11_cmp_c: R_LMPD_11_cmp_c -1000 BU_LMPD_11_cmp_c < 0
G_LMPD_12_cys-L_c: -0.59209 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +0.59209 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c -1.1842 LC_atp_c +0.59209 LC_coa_c +0.59209 LC_cys-L_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +2.3684 LC_nadp_c -2.3684 LC_nadph_c +1.1842 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_so4_c - DG_LMPD_12_cys-L_c + DGo_LMPD_12_cys-L_c = 0
FU_LMPD_12_cys-L_c: DG_LMPD_12_cys-L_c +1000 FU_LMPD_12_cys-L_c < 1000
BU_LMPD_12_cys-L_c: - DG_LMPD_12_cys-L_c +1000 BU_LMPD_12_cys-L_c < 1000
SU_LMPD_12_cys-L_c: FU_LMPD_12_cys-L_c + BU_LMPD_12_cys-L_c < 1
UF_LMPD_12_cys-L_c: F_LMPD_12_cys-L_c -1000 FU_LMPD_12_cys-L_c < 0
UR_LMPD_12_cys-L_c: R_LMPD_12_cys-L_c -1000 BU_LMPD_12_cys-L_c < 0
G_LMPD_13_damp_c: -1.1842 LC_3pg_c +2.9604 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -1.1842 LC_asp-L_c -4.1446 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_damp_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -2.3684 LC_nad_c +2.3684 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_13_damp_c + DGo_LMPD_13_damp_c = 0
FU_LMPD_13_damp_c: DG_LMPD_13_damp_c +1000 FU_LMPD_13_damp_c < 1000
BU_LMPD_13_damp_c: - DG_LMPD_13_damp_c +1000 BU_LMPD_13_damp_c < 1000
SU_LMPD_13_damp_c: FU_LMPD_13_damp_c + BU_LMPD_13_damp_c < 1
UF_LMPD_13_damp_c: F_LMPD_13_damp_c -1000 FU_LMPD_13_damp_c < 0
UR_LMPD_13_damp_c: R_LMPD_13_damp_c -1000 BU_LMPD_13_damp_c < 0
G_LMPD_14_damp_c: -1.1842 LC_3pg_c +1.7763 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -1.1842 LC_asp-L_c -3.5525 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_damp_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -2.3684 LC_nad_c +2.3684 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_14_damp_c + DGo_LMPD_14_damp_c = 0
FU_LMPD_14_damp_c: DG_LMPD_14_damp_c +1000 FU_LMPD_14_damp_c < 1000
BU_LMPD_14_damp_c: - DG_LMPD_14_damp_c +1000 BU_LMPD_14_damp_c < 1000
SU_LMPD_14_damp_c: FU_LMPD_14_damp_c + BU_LMPD_14_damp_c < 1
UF_LMPD_14_damp_c: F_LMPD_14_damp_c -1000 FU_LMPD_14_damp_c < 0
UR_LMPD_14_damp_c: R_LMPD_14_damp_c -1000 BU_LMPD_14_damp_c < 0
G_LMPD_15_damp_c: -1.1842 LC_3pg_c +2.9604 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -1.1842 LC_asp-L_c -4.1446 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_damp_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_15_damp_c + DGo_LMPD_15_damp_c = 0
FU_LMPD_15_damp_c: DG_LMPD_15_damp_c +1000 FU_LMPD_15_damp_c < 1000
BU_LMPD_15_damp_c: - DG_LMPD_15_damp_c +1000 BU_LMPD_15_damp_c < 1000
SU_LMPD_15_damp_c: FU_LMPD_15_damp_c + BU_LMPD_15_damp_c < 1
UF_LMPD_15_damp_c: F_LMPD_15_damp_c -1000 FU_LMPD_15_damp_c < 0
UR_LMPD_15_damp_c: R_LMPD_15_damp_c -1000 BU_LMPD_15_damp_c < 0
G_LMPD_16_damp_c: -1.1842 LC_3pg_c +1.7763 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -1.1842 LC_asp-L_c -3.5525 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_damp_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.1842 LC_nh4_c +5.3288 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_16_damp_c + DGo_LMPD_16_damp_c = 0
FU_LMPD_16_damp_c: DG_LMPD_16_damp_c +1000 FU_LMPD_16_damp_c < 1000
BU_LMPD_16_damp_c: - DG_LMPD_16_damp_c +1000 BU_LMPD_16_damp_c < 1000
SU_LMPD_16_damp_c: FU_LMPD_16_damp_c + BU_LMPD_16_damp_c < 1
UF_LMPD_16_damp_c: F_LMPD_16_damp_c -1000 FU_LMPD_16_damp_c < 0
UR_LMPD_16_damp_c: R_LMPD_16_damp_c -1000 BU_LMPD_16_damp_c < 0
G_LMPD_17_dcmp_c: +1.7763 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_co2_c +0.59209 LC_dcmp_c -0.59209 LC_fum_c -0.59209 LC_hco3_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_nh4_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_succ_c - DG_LMPD_17_dcmp_c + DGo_LMPD_17_dcmp_c = 0
FU_LMPD_17_dcmp_c: DG_LMPD_17_dcmp_c +1000 FU_LMPD_17_dcmp_c < 1000
BU_LMPD_17_dcmp_c: - DG_LMPD_17_dcmp_c +1000 BU_LMPD_17_dcmp_c < 1000
SU_LMPD_17_dcmp_c: FU_LMPD_17_dcmp_c + BU_LMPD_17_dcmp_c < 1
UF_LMPD_17_dcmp_c: F_LMPD_17_dcmp_c -1000 FU_LMPD_17_dcmp_c < 0
UR_LMPD_17_dcmp_c: R_LMPD_17_dcmp_c -1000 BU_LMPD_17_dcmp_c < 0
G_LMPD_18_dcmp_c: +1.7763 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_co2_c +0.59209 LC_dcmp_c -0.59209 LC_hco3_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_nh4_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_LMPD_18_dcmp_c + DGo_LMPD_18_dcmp_c = 0
FU_LMPD_18_dcmp_c: DG_LMPD_18_dcmp_c +1000 FU_LMPD_18_dcmp_c < 1000
BU_LMPD_18_dcmp_c: - DG_LMPD_18_dcmp_c +1000 BU_LMPD_18_dcmp_c < 1000
SU_LMPD_18_dcmp_c: FU_LMPD_18_dcmp_c + BU_LMPD_18_dcmp_c < 1
UF_LMPD_18_dcmp_c: F_LMPD_18_dcmp_c -1000 FU_LMPD_18_dcmp_c < 0
UR_LMPD_18_dcmp_c: R_LMPD_18_dcmp_c -1000 BU_LMPD_18_dcmp_c < 0
G_LMPD_19_dgmp_c: -1.1842 LC_3pg_c +4.7367 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -5.9209 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_dgmp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -2.9604 LC_nad_c +2.9604 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_19_dgmp_c + DGo_LMPD_19_dgmp_c = 0
FU_LMPD_19_dgmp_c: DG_LMPD_19_dgmp_c +1000 FU_LMPD_19_dgmp_c < 1000
BU_LMPD_19_dgmp_c: - DG_LMPD_19_dgmp_c +1000 BU_LMPD_19_dgmp_c < 1000
SU_LMPD_19_dgmp_c: FU_LMPD_19_dgmp_c + BU_LMPD_19_dgmp_c < 1
UF_LMPD_19_dgmp_c: F_LMPD_19_dgmp_c -1000 FU_LMPD_19_dgmp_c < 0
UR_LMPD_19_dgmp_c: R_LMPD_19_dgmp_c -1000 BU_LMPD_19_dgmp_c < 0
G_LMPD_20_dgmp_c: -1.1842 LC_3pg_c +4.1446 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -5.3288 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_dgmp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -2.9604 LC_nad_c +2.9604 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.7763 LC_nh4_c +5.3288 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_20_dgmp_c + DGo_LMPD_20_dgmp_c = 0
FU_LMPD_20_dgmp_c: DG_LMPD_20_dgmp_c +1000 FU_LMPD_20_dgmp_c < 1000
BU_LMPD_20_dgmp_c: - DG_LMPD_20_dgmp_c +1000 BU_LMPD_20_dgmp_c < 1000
SU_LMPD_20_dgmp_c: FU_LMPD_20_dgmp_c + BU_LMPD_20_dgmp_c < 1
UF_LMPD_20_dgmp_c: F_LMPD_20_dgmp_c -1000 FU_LMPD_20_dgmp_c < 0
UR_LMPD_20_dgmp_c: R_LMPD_20_dgmp_c -1000 BU_LMPD_20_dgmp_c < 0
G_LMPD_21_dgmp_c: -1.1842 LC_3pg_c +4.1446 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -5.3288 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_dgmp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.7763 LC_nh4_c +5.3288 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_21_dgmp_c + DGo_LMPD_21_dgmp_c = 0
FU_LMPD_21_dgmp_c: DG_LMPD_21_dgmp_c +1000 FU_LMPD_21_dgmp_c < 1000
BU_LMPD_21_dgmp_c: - DG_LMPD_21_dgmp_c +1000 BU_LMPD_21_dgmp_c < 1000
SU_LMPD_21_dgmp_c: FU_LMPD_21_dgmp_c + BU_LMPD_21_dgmp_c < 1
UF_LMPD_21_dgmp_c: F_LMPD_21_dgmp_c -1000 FU_LMPD_21_dgmp_c < 0
UR_LMPD_21_dgmp_c: R_LMPD_21_dgmp_c -1000 BU_LMPD_21_dgmp_c < 0
G_LMPD_22_dgmp_c: -1.1842 LC_3pg_c +4.7367 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -5.9209 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_dgmp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_22_dgmp_c + DGo_LMPD_22_dgmp_c = 0
FU_LMPD_22_dgmp_c: DG_LMPD_22_dgmp_c +1000 FU_LMPD_22_dgmp_c < 1000
BU_LMPD_22_dgmp_c: - DG_LMPD_22_dgmp_c +1000 BU_LMPD_22_dgmp_c < 1000
SU_LMPD_22_dgmp_c: FU_LMPD_22_dgmp_c + BU_LMPD_22_dgmp_c < 1
UF_LMPD_22_dgmp_c: F_LMPD_22_dgmp_c -1000 FU_LMPD_22_dgmp_c < 0
UR_LMPD_22_dgmp_c: R_LMPD_22_dgmp_c -1000 BU_LMPD_22_dgmp_c < 0
G_LMPD_23_dtmp_c: -0.59209 LC_3pg_c +1.7763 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_co2_c +0.59209 LC_dtmp_c -0.59209 LC_glu-L_c +0.59209 LC_gly_c -0.59209 LC_hco3_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_LMPD_23_dtmp_c + DGo_LMPD_23_dtmp_c = 0
FU_LMPD_23_dtmp_c: DG_LMPD_23_dtmp_c +1000 FU_LMPD_23_dtmp_c < 1000
BU_LMPD_23_dtmp_c: - DG_LMPD_23_dtmp_c +1000 BU_LMPD_23_dtmp_c < 1000
SU_LMPD_23_dtmp_c: FU_LMPD_23_dtmp_c + BU_LMPD_23_dtmp_c < 1
UF_LMPD_23_dtmp_c: F_LMPD_23_dtmp_c -1000 FU_LMPD_23_dtmp_c < 0
UR_LMPD_23_dtmp_c: R_LMPD_23_dtmp_c -1000 BU_LMPD_23_dtmp_c < 0
G_LMPD_24_dtmp_c: -0.59209 LC_3pg_c +1.7763 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_co2_c +0.59209 LC_dtmp_c -0.59209 LC_fum_c -0.59209 LC_glu-L_c +0.59209 LC_gly_c -0.59209 LC_hco3_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_succ_c - DG_LMPD_24_dtmp_c + DGo_LMPD_24_dtmp_c = 0
FU_LMPD_24_dtmp_c: DG_LMPD_24_dtmp_c +1000 FU_LMPD_24_dtmp_c < 1000
BU_LMPD_24_dtmp_c: - DG_LMPD_24_dtmp_c +1000 BU_LMPD_24_dtmp_c < 1000
SU_LMPD_24_dtmp_c: FU_LMPD_24_dtmp_c + BU_LMPD_24_dtmp_c < 1
UF_LMPD_24_dtmp_c: F_LMPD_24_dtmp_c -1000 FU_LMPD_24_dtmp_c < 0
UR_LMPD_24_dtmp_c: R_LMPD_24_dtmp_c -1000 BU_LMPD_24_dtmp_c < 0
G_LMPD_25_ergst_c: -10.6576 LC_accoa_c +11.8418 LC_adp_c +0.59209 LC_amp_c -12.4339 LC_atp_c +4.7367 LC_co2_c +10.6576 LC_coa_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +15.9864 LC_nadp_c -15.9864 LC_nadph_c -5.9209 LC_o2_c +4.7367 LC_pi_c +4.1446 LC_ppi_c +0.59209 LC_ergst_r +1.1842 LC_nadp_r -1.1842 LC_nadph_r -0.59209 LC_o2_r - DG_LMPD_25_ergst_c + DGo_LMPD_25_ergst_c = 0
FU_LMPD_25_ergst_c: DG_LMPD_25_ergst_c +1000 FU_LMPD_25_ergst_c < 1000
BU_LMPD_25_ergst_c: - DG_LMPD_25_ergst_c +1000 BU_LMPD_25_ergst_c < 1000
SU_LMPD_25_ergst_c: FU_LMPD_25_ergst_c + BU_LMPD_25_ergst_c < 1
UF_LMPD_25_ergst_c: F_LMPD_25_ergst_c -1000 FU_LMPD_25_ergst_c < 0
UR_LMPD_25_ergst_c: R_LMPD_25_ergst_c -1000 BU_LMPD_25_ergst_c < 0
G_LMPD_26_ergst_c: -10.6576 LC_accoa_c +11.8418 LC_adp_c +0.59209 LC_amp_c -12.4339 LC_atp_c +4.7367 LC_co2_c +10.6576 LC_coa_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +14.2101 LC_nadp_c -14.2101 LC_nadph_c -5.9209 LC_o2_c +4.7367 LC_pi_c +4.1446 LC_ppi_c +0.59209 LC_ergst_r +1.1842 LC_nadp_r -1.1842 LC_nadph_r -0.59209 LC_o2_r - DG_LMPD_26_ergst_c + DGo_LMPD_26_ergst_c = 0
FU_LMPD_26_ergst_c: DG_LMPD_26_ergst_c +1000 FU_LMPD_26_ergst_c < 1000
BU_LMPD_26_ergst_c: - DG_LMPD_26_ergst_c +1000 BU_LMPD_26_ergst_c < 1000
SU_LMPD_26_ergst_c: FU_LMPD_26_ergst_c + BU_LMPD_26_ergst_c < 1
UF_LMPD_26_ergst_c: F_LMPD_26_ergst_c -1000 FU_LMPD_26_ergst_c < 0
UR_LMPD_26_ergst_c: R_LMPD_26_ergst_c -1000 BU_LMPD_26_ergst_c < 0
G_LMPD_27_gln-L_c: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_gln-L_c -0.59209 LC_glu-L_c -0.59209 LC_nh4_c +0.59209 LC_pi_c - DG_LMPD_27_gln-L_c + DGo_LMPD_27_gln-L_c = 0
FU_LMPD_27_gln-L_c: DG_LMPD_27_gln-L_c +1000 FU_LMPD_27_gln-L_c < 1000
BU_LMPD_27_gln-L_c: - DG_LMPD_27_gln-L_c +1000 BU_LMPD_27_gln-L_c < 1000
SU_LMPD_27_gln-L_c: FU_LMPD_27_gln-L_c + BU_LMPD_27_gln-L_c < 1
UF_LMPD_27_gln-L_c: F_LMPD_27_gln-L_c -1000 FU_LMPD_27_gln-L_c < 0
UR_LMPD_27_gln-L_c: R_LMPD_27_gln-L_c -1000 BU_LMPD_27_gln-L_c < 0
G_LMPD_28_gly_c: -0.59209 LC_ala-L_c +0.59209 LC_gly_c -0.59209 LC_icit_c +0.59209 LC_pyr_c +0.59209 LC_succ_c - DG_LMPD_28_gly_c + DGo_LMPD_28_gly_c = 0
FU_LMPD_28_gly_c: DG_LMPD_28_gly_c +1000 FU_LMPD_28_gly_c < 1000
BU_LMPD_28_gly_c: - DG_LMPD_28_gly_c +1000 BU_LMPD_28_gly_c < 1000
SU_LMPD_28_gly_c: FU_LMPD_28_gly_c + BU_LMPD_28_gly_c < 1
UF_LMPD_28_gly_c: F_LMPD_28_gly_c -1000 FU_LMPD_28_gly_c < 0
UR_LMPD_28_gly_c: R_LMPD_28_gly_c -1000 BU_LMPD_28_gly_c < 0
SU_LMPD_29_glycogen_c: FU_LMPD_29_glycogen_c + BU_LMPD_29_glycogen_c < 1
UF_LMPD_29_glycogen_c: F_LMPD_29_glycogen_c -1000 FU_LMPD_29_glycogen_c < 0
UR_LMPD_29_glycogen_c: R_LMPD_29_glycogen_c -1000 BU_LMPD_29_glycogen_c < 0
G_LMPD_30_gmp_c: -1.1842 LC_3pg_c +4.1446 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -5.3288 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_gmp_c -2.9604 LC_nad_c +2.9604 LC_nadh_c -1.7763 LC_nh4_c +5.3288 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_30_gmp_c + DGo_LMPD_30_gmp_c = 0
FU_LMPD_30_gmp_c: DG_LMPD_30_gmp_c +1000 FU_LMPD_30_gmp_c < 1000
BU_LMPD_30_gmp_c: - DG_LMPD_30_gmp_c +1000 BU_LMPD_30_gmp_c < 1000
SU_LMPD_30_gmp_c: FU_LMPD_30_gmp_c + BU_LMPD_30_gmp_c < 1
UF_LMPD_30_gmp_c: F_LMPD_30_gmp_c -1000 FU_LMPD_30_gmp_c < 0
UR_LMPD_30_gmp_c: R_LMPD_30_gmp_c -1000 BU_LMPD_30_gmp_c < 0
G_LMPD_31_gmp_c: -1.1842 LC_3pg_c +4.1446 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -5.3288 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_gmp_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.7763 LC_nh4_c +5.3288 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_31_gmp_c + DGo_LMPD_31_gmp_c = 0
FU_LMPD_31_gmp_c: DG_LMPD_31_gmp_c +1000 FU_LMPD_31_gmp_c < 1000
BU_LMPD_31_gmp_c: - DG_LMPD_31_gmp_c +1000 BU_LMPD_31_gmp_c < 1000
SU_LMPD_31_gmp_c: FU_LMPD_31_gmp_c + BU_LMPD_31_gmp_c < 1
UF_LMPD_31_gmp_c: F_LMPD_31_gmp_c -1000 FU_LMPD_31_gmp_c < 0
UR_LMPD_31_gmp_c: R_LMPD_31_gmp_c -1000 BU_LMPD_31_gmp_c < 0
G_LMPD_32_his-L_c: -0.59209 LC_3pg_c +1.1842 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -1.7763 LC_atp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_his-L_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_32_his-L_c + DGo_LMPD_32_his-L_c = 0
FU_LMPD_32_his-L_c: DG_LMPD_32_his-L_c +1000 FU_LMPD_32_his-L_c < 1000
BU_LMPD_32_his-L_c: - DG_LMPD_32_his-L_c +1000 BU_LMPD_32_his-L_c < 1000
SU_LMPD_32_his-L_c: FU_LMPD_32_his-L_c + BU_LMPD_32_his-L_c < 1
UF_LMPD_32_his-L_c: F_LMPD_32_his-L_c -1000 FU_LMPD_32_his-L_c < 0
UR_LMPD_32_his-L_c: R_LMPD_32_his-L_c -1000 BU_LMPD_32_his-L_c < 0
G_LMPD_33_his-L_c: -0.59209 LC_3pg_c +1.1842 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_his-L_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_33_his-L_c + DGo_LMPD_33_his-L_c = 0
FU_LMPD_33_his-L_c: DG_LMPD_33_his-L_c +1000 FU_LMPD_33_his-L_c < 1000
BU_LMPD_33_his-L_c: - DG_LMPD_33_his-L_c +1000 BU_LMPD_33_his-L_c < 1000
SU_LMPD_33_his-L_c: FU_LMPD_33_his-L_c + BU_LMPD_33_his-L_c < 1
UF_LMPD_33_his-L_c: F_LMPD_33_his-L_c -1000 FU_LMPD_33_his-L_c < 0
UR_LMPD_33_his-L_c: R_LMPD_33_his-L_c -1000 BU_LMPD_33_his-L_c < 0
G_LMPD_34_his-L_c: -0.59209 LC_3pg_c +1.1842 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_his-L_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_34_his-L_c + DGo_LMPD_34_his-L_c = 0
FU_LMPD_34_his-L_c: DG_LMPD_34_his-L_c +1000 FU_LMPD_34_his-L_c < 1000
BU_LMPD_34_his-L_c: - DG_LMPD_34_his-L_c +1000 BU_LMPD_34_his-L_c < 1000
SU_LMPD_34_his-L_c: FU_LMPD_34_his-L_c + BU_LMPD_34_his-L_c < 1
UF_LMPD_34_his-L_c: F_LMPD_34_his-L_c -1000 FU_LMPD_34_his-L_c < 0
UR_LMPD_34_his-L_c: R_LMPD_34_his-L_c -1000 BU_LMPD_34_his-L_c < 0
G_LMPD_35_his-L_c: -0.59209 LC_3pg_c +1.1842 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -1.7763 LC_atp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_his-L_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_35_his-L_c + DGo_LMPD_35_his-L_c = 0
FU_LMPD_35_his-L_c: DG_LMPD_35_his-L_c +1000 FU_LMPD_35_his-L_c < 1000
BU_LMPD_35_his-L_c: - DG_LMPD_35_his-L_c +1000 BU_LMPD_35_his-L_c < 1000
SU_LMPD_35_his-L_c: FU_LMPD_35_his-L_c + BU_LMPD_35_his-L_c < 1
UF_LMPD_35_his-L_c: F_LMPD_35_his-L_c -1000 FU_LMPD_35_his-L_c < 0
UR_LMPD_35_his-L_c: R_LMPD_35_his-L_c -1000 BU_LMPD_35_his-L_c < 0
G_LMPD_36_ile-L_c: +1.1842 LC_adp_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c -0.59209 LC_h2o_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +1.1842 LC_pi_c +0.59209 LC_akg_m +0.59209 LC_co2_m -0.59209 LC_glu-L_m +0.59209 LC_h2o_m +0.59209 LC_ile-L_m +0.59209 LC_nadp_m -0.59209 LC_nadph_m +0.59209 LC_nh4_m -0.59209 LC_pyr_m - DG_LMPD_36_ile-L_c + DGo_LMPD_36_ile-L_c = 0
FU_LMPD_36_ile-L_c: DG_LMPD_36_ile-L_c +1000 FU_LMPD_36_ile-L_c < 1000
BU_LMPD_36_ile-L_c: - DG_LMPD_36_ile-L_c +1000 BU_LMPD_36_ile-L_c < 1000
SU_LMPD_36_ile-L_c: FU_LMPD_36_ile-L_c + BU_LMPD_36_ile-L_c < 1
UF_LMPD_36_ile-L_c: F_LMPD_36_ile-L_c -1000 FU_LMPD_36_ile-L_c < 0
UR_LMPD_36_ile-L_c: R_LMPD_36_ile-L_c -1000 BU_LMPD_36_ile-L_c < 0
G_LMPD_37_ile-L_c: +1.1842 LC_adp_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c -0.59209 LC_h2o_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c +1.1842 LC_pi_c +0.59209 LC_akg_m +0.59209 LC_co2_m -0.59209 LC_glu-L_m +0.59209 LC_h2o_m +0.59209 LC_ile-L_m +0.59209 LC_nadp_m -0.59209 LC_nadph_m +0.59209 LC_nh4_m -0.59209 LC_pyr_m - DG_LMPD_37_ile-L_c + DGo_LMPD_37_ile-L_c = 0
FU_LMPD_37_ile-L_c: DG_LMPD_37_ile-L_c +1000 FU_LMPD_37_ile-L_c < 1000
BU_LMPD_37_ile-L_c: - DG_LMPD_37_ile-L_c +1000 BU_LMPD_37_ile-L_c < 1000
SU_LMPD_37_ile-L_c: FU_LMPD_37_ile-L_c + BU_LMPD_37_ile-L_c < 1
UF_LMPD_37_ile-L_c: F_LMPD_37_ile-L_c -1000 FU_LMPD_37_ile-L_c < 0
UR_LMPD_37_ile-L_c: R_LMPD_37_ile-L_c -1000 BU_LMPD_37_ile-L_c < 0
G_LMPD_38_leu-L_c: -0.59209 LC_accoa_c +0.59209 LC_akg_c +0.59209 LC_co2_c +0.59209 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_h2o_c +0.59209 LC_leu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_co2_m +0.59209 LC_h2o_m +0.59209 LC_nadp_m -0.59209 LC_nadph_m -1.1842 LC_pyr_m - DG_LMPD_38_leu-L_c + DGo_LMPD_38_leu-L_c = 0
FU_LMPD_38_leu-L_c: DG_LMPD_38_leu-L_c +1000 FU_LMPD_38_leu-L_c < 1000
BU_LMPD_38_leu-L_c: - DG_LMPD_38_leu-L_c +1000 BU_LMPD_38_leu-L_c < 1000
SU_LMPD_38_leu-L_c: FU_LMPD_38_leu-L_c + BU_LMPD_38_leu-L_c < 1
UF_LMPD_38_leu-L_c: F_LMPD_38_leu-L_c -1000 FU_LMPD_38_leu-L_c < 0
UR_LMPD_38_leu-L_c: R_LMPD_38_leu-L_c -1000 BU_LMPD_38_leu-L_c < 0
G_LMPD_39_leu-L_c: +0.59209 LC_akg_c +0.59209 LC_co2_c -0.59209 LC_glu-L_c +0.59209 LC_leu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c -0.59209 LC_accoa_m +0.59209 LC_co2_m +0.59209 LC_coa_m +0.59209 LC_nadp_m -0.59209 LC_nadph_m -1.1842 LC_pyr_m - DG_LMPD_39_leu-L_c + DGo_LMPD_39_leu-L_c = 0
FU_LMPD_39_leu-L_c: DG_LMPD_39_leu-L_c +1000 FU_LMPD_39_leu-L_c < 1000
BU_LMPD_39_leu-L_c: - DG_LMPD_39_leu-L_c +1000 BU_LMPD_39_leu-L_c < 1000
SU_LMPD_39_leu-L_c: FU_LMPD_39_leu-L_c + BU_LMPD_39_leu-L_c < 1
UF_LMPD_39_leu-L_c: F_LMPD_39_leu-L_c -1000 FU_LMPD_39_leu-L_c < 0
UR_LMPD_39_leu-L_c: R_LMPD_39_leu-L_c -1000 BU_LMPD_39_leu-L_c < 0
G_LMPD_40_lys-L_c: +1.1842 LC_akg_c +0.59209 LC_amp_c -0.59209 LC_atp_c -1.1842 LC_glu-L_c +0.59209 LC_lys-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +0.59209 LC_ppi_c +0.59209 LC_co2_m -0.59209 LC_hcit_m -0.59209 LC_nad_m +0.59209 LC_nadh_m - DG_LMPD_40_lys-L_c + DGo_LMPD_40_lys-L_c = 0
FU_LMPD_40_lys-L_c: DG_LMPD_40_lys-L_c +1000 FU_LMPD_40_lys-L_c < 1000
BU_LMPD_40_lys-L_c: - DG_LMPD_40_lys-L_c +1000 BU_LMPD_40_lys-L_c < 1000
SU_LMPD_40_lys-L_c: FU_LMPD_40_lys-L_c + BU_LMPD_40_lys-L_c < 1
UF_LMPD_40_lys-L_c: F_LMPD_40_lys-L_c -1000 FU_LMPD_40_lys-L_c < 0
UR_LMPD_40_lys-L_c: R_LMPD_40_lys-L_c -1000 BU_LMPD_40_lys-L_c < 0
G_LMPD_41_lys-L_c: +1.1842 LC_akg_c +0.59209 LC_amp_c -0.59209 LC_atp_c -1.1842 LC_glu-L_c +0.59209 LC_lys-L_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c +0.59209 LC_ppi_c +0.59209 LC_co2_m -0.59209 LC_hcit_m -0.59209 LC_nad_m +0.59209 LC_nadh_m - DG_LMPD_41_lys-L_c + DGo_LMPD_41_lys-L_c = 0
FU_LMPD_41_lys-L_c: DG_LMPD_41_lys-L_c +1000 FU_LMPD_41_lys-L_c < 1000
BU_LMPD_41_lys-L_c: - DG_LMPD_41_lys-L_c +1000 BU_LMPD_41_lys-L_c < 1000
SU_LMPD_41_lys-L_c: FU_LMPD_41_lys-L_c + BU_LMPD_41_lys-L_c < 1
UF_LMPD_41_lys-L_c: F_LMPD_41_lys-L_c -1000 FU_LMPD_41_lys-L_c < 0
UR_LMPD_41_lys-L_c: R_LMPD_41_lys-L_c -1000 BU_LMPD_41_lys-L_c < 0
SU_LMPD_42_mannan_c: FU_LMPD_42_mannan_c + BU_LMPD_42_mannan_c < 1
UF_LMPD_42_mannan_c: F_LMPD_42_mannan_c -1000 FU_LMPD_42_mannan_c < 0
UR_LMPD_42_mannan_c: R_LMPD_42_mannan_c -1000 BU_LMPD_42_mannan_c < 0
SU_LMPD_43_mannan_c: FU_LMPD_43_mannan_c + BU_LMPD_43_mannan_c < 1
UF_LMPD_43_mannan_c: F_LMPD_43_mannan_c -1000 FU_LMPD_43_mannan_c < 0
UR_LMPD_43_mannan_c: R_LMPD_43_mannan_c -1000 BU_LMPD_43_mannan_c < 0
G_LMPD_44_met-L_c: -0.59209 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -1.7763 LC_atp_c +0.59209 LC_coa_c -0.59209 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_met-L_c +3.5525 LC_nadp_c -3.5525 LC_nadph_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_so4_c - DG_LMPD_44_met-L_c + DGo_LMPD_44_met-L_c = 0
FU_LMPD_44_met-L_c: DG_LMPD_44_met-L_c +1000 FU_LMPD_44_met-L_c < 1000
BU_LMPD_44_met-L_c: - DG_LMPD_44_met-L_c +1000 BU_LMPD_44_met-L_c < 1000
SU_LMPD_44_met-L_c: FU_LMPD_44_met-L_c + BU_LMPD_44_met-L_c < 1
UF_LMPD_44_met-L_c: F_LMPD_44_met-L_c -1000 FU_LMPD_44_met-L_c < 0
UR_LMPD_44_met-L_c: R_LMPD_44_met-L_c -1000 BU_LMPD_44_met-L_c < 0
G_LMPD_45_met-L_c: -0.59209 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -1.7763 LC_atp_c +0.59209 LC_coa_c -0.59209 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_met-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +4.1446 LC_nadp_c -4.1446 LC_nadph_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_so4_c - DG_LMPD_45_met-L_c + DGo_LMPD_45_met-L_c = 0
FU_LMPD_45_met-L_c: DG_LMPD_45_met-L_c +1000 FU_LMPD_45_met-L_c < 1000
BU_LMPD_45_met-L_c: - DG_LMPD_45_met-L_c +1000 BU_LMPD_45_met-L_c < 1000
SU_LMPD_45_met-L_c: FU_LMPD_45_met-L_c + BU_LMPD_45_met-L_c < 1
UF_LMPD_45_met-L_c: F_LMPD_45_met-L_c -1000 FU_LMPD_45_met-L_c < 0
UR_LMPD_45_met-L_c: R_LMPD_45_met-L_c -1000 BU_LMPD_45_met-L_c < 0
SU_LMPD_47_pa_SC_c: FU_LMPD_47_pa_SC_c + BU_LMPD_47_pa_SC_c < 1
UF_LMPD_47_pa_SC_c: F_LMPD_47_pa_SC_c -1000 FU_LMPD_47_pa_SC_c < 0
UR_LMPD_47_pa_SC_c: R_LMPD_47_pa_SC_c -1000 BU_LMPD_47_pa_SC_c < 0
SU_LMPD_48_pc_SC_c: FU_LMPD_48_pc_SC_c + BU_LMPD_48_pc_SC_c < 1
UF_LMPD_48_pc_SC_c: F_LMPD_48_pc_SC_c -1000 FU_LMPD_48_pc_SC_c < 0
UR_LMPD_48_pc_SC_c: R_LMPD_48_pc_SC_c -1000 BU_LMPD_48_pc_SC_c < 0
SU_LMPD_49_pc_SC_c: FU_LMPD_49_pc_SC_c + BU_LMPD_49_pc_SC_c < 1
UF_LMPD_49_pc_SC_c: F_LMPD_49_pc_SC_c -1000 FU_LMPD_49_pc_SC_c < 0
UR_LMPD_49_pc_SC_c: R_LMPD_49_pc_SC_c -1000 BU_LMPD_49_pc_SC_c < 0
SU_LMPD_51_pe_SC_c: FU_LMPD_51_pe_SC_c + BU_LMPD_51_pe_SC_c < 1
UF_LMPD_51_pe_SC_c: F_LMPD_51_pe_SC_c -1000 FU_LMPD_51_pe_SC_c < 0
UR_LMPD_51_pe_SC_c: R_LMPD_51_pe_SC_c -1000 BU_LMPD_51_pe_SC_c < 0
G_LMPD_52_phe-L_c: +0.59209 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_glu-L_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_pep_c +0.59209 LC_phe-L_c +2.3684 LC_pi_c - DG_LMPD_52_phe-L_c + DGo_LMPD_52_phe-L_c = 0
FU_LMPD_52_phe-L_c: DG_LMPD_52_phe-L_c +1000 FU_LMPD_52_phe-L_c < 1000
BU_LMPD_52_phe-L_c: - DG_LMPD_52_phe-L_c +1000 BU_LMPD_52_phe-L_c < 1000
SU_LMPD_52_phe-L_c: FU_LMPD_52_phe-L_c + BU_LMPD_52_phe-L_c < 1
UF_LMPD_52_phe-L_c: F_LMPD_52_phe-L_c -1000 FU_LMPD_52_phe-L_c < 0
UR_LMPD_52_phe-L_c: R_LMPD_52_phe-L_c -1000 BU_LMPD_52_phe-L_c < 0
G_LMPD_53_pro-L_c: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_glu-L_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c +0.59209 LC_pi_c +0.59209 LC_pro-L_c - DG_LMPD_53_pro-L_c + DGo_LMPD_53_pro-L_c = 0
FU_LMPD_53_pro-L_c: DG_LMPD_53_pro-L_c +1000 FU_LMPD_53_pro-L_c < 1000
BU_LMPD_53_pro-L_c: - DG_LMPD_53_pro-L_c +1000 BU_LMPD_53_pro-L_c < 1000
SU_LMPD_53_pro-L_c: FU_LMPD_53_pro-L_c + BU_LMPD_53_pro-L_c < 1
UF_LMPD_53_pro-L_c: F_LMPD_53_pro-L_c -1000 FU_LMPD_53_pro-L_c < 0
UR_LMPD_53_pro-L_c: R_LMPD_53_pro-L_c -1000 BU_LMPD_53_pro-L_c < 0
G_LMPD_54_pro-L_c: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_glu-L_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +0.59209 LC_pi_c +0.59209 LC_pro-L_c - DG_LMPD_54_pro-L_c + DGo_LMPD_54_pro-L_c = 0
FU_LMPD_54_pro-L_c: DG_LMPD_54_pro-L_c +1000 FU_LMPD_54_pro-L_c < 1000
BU_LMPD_54_pro-L_c: - DG_LMPD_54_pro-L_c +1000 BU_LMPD_54_pro-L_c < 1000
SU_LMPD_54_pro-L_c: FU_LMPD_54_pro-L_c + BU_LMPD_54_pro-L_c < 1
UF_LMPD_54_pro-L_c: F_LMPD_54_pro-L_c -1000 FU_LMPD_54_pro-L_c < 0
UR_LMPD_54_pro-L_c: R_LMPD_54_pro-L_c -1000 BU_LMPD_54_pro-L_c < 0
SU_LMPD_56_ps_SC_c: FU_LMPD_56_ps_SC_c + BU_LMPD_56_ps_SC_c < 1
UF_LMPD_56_ps_SC_c: F_LMPD_56_ps_SC_c -1000 FU_LMPD_56_ps_SC_c < 0
UR_LMPD_56_ps_SC_c: R_LMPD_56_ps_SC_c -1000 BU_LMPD_56_ps_SC_c < 0
SU_LMPD_57_ptd1ino_SC_c: FU_LMPD_57_ptd1ino_SC_c + BU_LMPD_57_ptd1ino_SC_c < 1
UF_LMPD_57_ptd1ino_SC_c: F_LMPD_57_ptd1ino_SC_c -1000 FU_LMPD_57_ptd1ino_SC_c < 0
UR_LMPD_57_ptd1ino_SC_c: R_LMPD_57_ptd1ino_SC_c -1000 BU_LMPD_57_ptd1ino_SC_c < 0
SU_LMPD_58_ptd1ino_SC_c: FU_LMPD_58_ptd1ino_SC_c + BU_LMPD_58_ptd1ino_SC_c < 1
UF_LMPD_58_ptd1ino_SC_c: F_LMPD_58_ptd1ino_SC_c -1000 FU_LMPD_58_ptd1ino_SC_c < 0
UR_LMPD_58_ptd1ino_SC_c: R_LMPD_58_ptd1ino_SC_c -1000 BU_LMPD_58_ptd1ino_SC_c < 0
G_LMPD_59_ribflv_c: +5.3288 LC_adp_c -0.59209 LC_ala-L_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_co2_c -1.1842 LC_db4p_c -0.59209 LC_for_c +0.59209 LC_fum_c -0.59209 LC_icit_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_nh4_c +7.1051 LC_pi_c +1.7763 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c +0.59209 LC_ribflv_c +0.59209 LC_succ_c - DG_LMPD_59_ribflv_c + DGo_LMPD_59_ribflv_c = 0
FU_LMPD_59_ribflv_c: DG_LMPD_59_ribflv_c +1000 FU_LMPD_59_ribflv_c < 1000
BU_LMPD_59_ribflv_c: - DG_LMPD_59_ribflv_c +1000 BU_LMPD_59_ribflv_c < 1000
SU_LMPD_59_ribflv_c: FU_LMPD_59_ribflv_c + BU_LMPD_59_ribflv_c < 1
UF_LMPD_59_ribflv_c: F_LMPD_59_ribflv_c -1000 FU_LMPD_59_ribflv_c < 0
UR_LMPD_59_ribflv_c: R_LMPD_59_ribflv_c -1000 BU_LMPD_59_ribflv_c < 0
G_LMPD_60_ribflv_c: +6.513 LC_adp_c -0.59209 LC_ala-L_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -7.6972 LC_atp_c -0.59209 LC_co2_c -1.1842 LC_db4p_c -0.59209 LC_for_c +0.59209 LC_fum_c -0.59209 LC_icit_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_nh4_c +7.1051 LC_pi_c +1.7763 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c +0.59209 LC_ribflv_c +0.59209 LC_succ_c - DG_LMPD_60_ribflv_c + DGo_LMPD_60_ribflv_c = 0
FU_LMPD_60_ribflv_c: DG_LMPD_60_ribflv_c +1000 FU_LMPD_60_ribflv_c < 1000
BU_LMPD_60_ribflv_c: - DG_LMPD_60_ribflv_c +1000 BU_LMPD_60_ribflv_c < 1000
SU_LMPD_60_ribflv_c: FU_LMPD_60_ribflv_c + BU_LMPD_60_ribflv_c < 1
UF_LMPD_60_ribflv_c: F_LMPD_60_ribflv_c -1000 FU_LMPD_60_ribflv_c < 0
UR_LMPD_60_ribflv_c: R_LMPD_60_ribflv_c -1000 BU_LMPD_60_ribflv_c < 0
G_LMPD_61_ser-L_c: -0.59209 LC_3pg_c +0.59209 LC_akg_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_pi_c +0.59209 LC_ser-L_c - DG_LMPD_61_ser-L_c + DGo_LMPD_61_ser-L_c = 0
FU_LMPD_61_ser-L_c: DG_LMPD_61_ser-L_c +1000 FU_LMPD_61_ser-L_c < 1000
BU_LMPD_61_ser-L_c: - DG_LMPD_61_ser-L_c +1000 BU_LMPD_61_ser-L_c < 1000
SU_LMPD_61_ser-L_c: FU_LMPD_61_ser-L_c + BU_LMPD_61_ser-L_c < 1
UF_LMPD_61_ser-L_c: F_LMPD_61_ser-L_c -1000 FU_LMPD_61_ser-L_c < 0
UR_LMPD_61_ser-L_c: R_LMPD_61_ser-L_c -1000 BU_LMPD_61_ser-L_c < 0
G_LMPD_62_thr-L_c: +1.1842 LC_adp_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c +1.1842 LC_pi_c +0.59209 LC_thr-L_c - DG_LMPD_62_thr-L_c + DGo_LMPD_62_thr-L_c = 0
FU_LMPD_62_thr-L_c: DG_LMPD_62_thr-L_c +1000 FU_LMPD_62_thr-L_c < 1000
BU_LMPD_62_thr-L_c: - DG_LMPD_62_thr-L_c +1000 BU_LMPD_62_thr-L_c < 1000
SU_LMPD_62_thr-L_c: FU_LMPD_62_thr-L_c + BU_LMPD_62_thr-L_c < 1
UF_LMPD_62_thr-L_c: F_LMPD_62_thr-L_c -1000 FU_LMPD_62_thr-L_c < 0
UR_LMPD_62_thr-L_c: R_LMPD_62_thr-L_c -1000 BU_LMPD_62_thr-L_c < 0
G_LMPD_63_thr-L_c: +1.1842 LC_adp_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +1.1842 LC_pi_c +0.59209 LC_thr-L_c - DG_LMPD_63_thr-L_c + DGo_LMPD_63_thr-L_c = 0
FU_LMPD_63_thr-L_c: DG_LMPD_63_thr-L_c +1000 FU_LMPD_63_thr-L_c < 1000
BU_LMPD_63_thr-L_c: - DG_LMPD_63_thr-L_c +1000 BU_LMPD_63_thr-L_c < 1000
SU_LMPD_63_thr-L_c: FU_LMPD_63_thr-L_c + BU_LMPD_63_thr-L_c < 1
UF_LMPD_63_thr-L_c: F_LMPD_63_thr-L_c -1000 FU_LMPD_63_thr-L_c < 0
UR_LMPD_63_thr-L_c: R_LMPD_63_thr-L_c -1000 BU_LMPD_63_thr-L_c < 0
G_LMPD_64_tre_c: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_g1p_c -0.59209 LC_g6p_c +0.59209 LC_pi_c +0.59209 LC_ppi_c +0.59209 LC_tre_c - DG_LMPD_64_tre_c + DGo_LMPD_64_tre_c = 0
FU_LMPD_64_tre_c: DG_LMPD_64_tre_c +1000 FU_LMPD_64_tre_c < 1000
BU_LMPD_64_tre_c: - DG_LMPD_64_tre_c +1000 BU_LMPD_64_tre_c < 1000
SU_LMPD_64_tre_c: FU_LMPD_64_tre_c + BU_LMPD_64_tre_c < 1
UF_LMPD_64_tre_c: F_LMPD_64_tre_c -1000 FU_LMPD_64_tre_c < 0
UR_LMPD_64_tre_c: R_LMPD_64_tre_c -1000 BU_LMPD_64_tre_c < 0
SU_LMPD_65_triglyc_SC_c: FU_LMPD_65_triglyc_SC_c + BU_LMPD_65_triglyc_SC_c < 1
UF_LMPD_65_triglyc_SC_c: F_LMPD_65_triglyc_SC_c -1000 FU_LMPD_65_triglyc_SC_c < 0
UR_LMPD_65_triglyc_SC_c: R_LMPD_65_triglyc_SC_c -1000 BU_LMPD_65_triglyc_SC_c < 0
SU_LMPD_66_triglyc_SC_c: FU_LMPD_66_triglyc_SC_c + BU_LMPD_66_triglyc_SC_c < 1
UF_LMPD_66_triglyc_SC_c: F_LMPD_66_triglyc_SC_c -1000 FU_LMPD_66_triglyc_SC_c < 0
UR_LMPD_66_triglyc_SC_c: R_LMPD_66_triglyc_SC_c -1000 BU_LMPD_66_triglyc_SC_c < 0
G_LMPD_67_trp-L_c: -0.59209 LC_3pg_c +1.1842 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c -1.7763 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_g3p_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -0.59209 LC_nh4_c -1.1842 LC_pep_c +3.5525 LC_pi_c +0.59209 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c +0.59209 LC_trp-L_c - DG_LMPD_67_trp-L_c + DGo_LMPD_67_trp-L_c = 0
FU_LMPD_67_trp-L_c: DG_LMPD_67_trp-L_c +1000 FU_LMPD_67_trp-L_c < 1000
BU_LMPD_67_trp-L_c: - DG_LMPD_67_trp-L_c +1000 BU_LMPD_67_trp-L_c < 1000
SU_LMPD_67_trp-L_c: FU_LMPD_67_trp-L_c + BU_LMPD_67_trp-L_c < 1
UF_LMPD_67_trp-L_c: F_LMPD_67_trp-L_c -1000 FU_LMPD_67_trp-L_c < 0
UR_LMPD_67_trp-L_c: R_LMPD_67_trp-L_c -1000 BU_LMPD_67_trp-L_c < 0
G_LMPD_68_tyr-L_c: +0.59209 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_pep_c +2.3684 LC_pi_c +0.59209 LC_tyr-L_c - DG_LMPD_68_tyr-L_c + DGo_LMPD_68_tyr-L_c = 0
FU_LMPD_68_tyr-L_c: DG_LMPD_68_tyr-L_c +1000 FU_LMPD_68_tyr-L_c < 1000
BU_LMPD_68_tyr-L_c: - DG_LMPD_68_tyr-L_c +1000 BU_LMPD_68_tyr-L_c < 1000
SU_LMPD_68_tyr-L_c: FU_LMPD_68_tyr-L_c + BU_LMPD_68_tyr-L_c < 1
UF_LMPD_68_tyr-L_c: F_LMPD_68_tyr-L_c -1000 FU_LMPD_68_tyr-L_c < 0
UR_LMPD_68_tyr-L_c: R_LMPD_68_tyr-L_c -1000 BU_LMPD_68_tyr-L_c < 0
G_LMPD_69_tyr-L_c: +0.59209 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_glu-L_c -1.1842 LC_pep_c +2.3684 LC_pi_c +0.59209 LC_tyr-L_c - DG_LMPD_69_tyr-L_c + DGo_LMPD_69_tyr-L_c = 0
FU_LMPD_69_tyr-L_c: DG_LMPD_69_tyr-L_c +1000 FU_LMPD_69_tyr-L_c < 1000
BU_LMPD_69_tyr-L_c: - DG_LMPD_69_tyr-L_c +1000 BU_LMPD_69_tyr-L_c < 1000
SU_LMPD_69_tyr-L_c: FU_LMPD_69_tyr-L_c + BU_LMPD_69_tyr-L_c < 1
UF_LMPD_69_tyr-L_c: F_LMPD_69_tyr-L_c -1000 FU_LMPD_69_tyr-L_c < 0
UR_LMPD_69_tyr-L_c: R_LMPD_69_tyr-L_c -1000 BU_LMPD_69_tyr-L_c < 0
G_LMPD_70_ump_c: +1.7763 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_fum_c -0.59209 LC_hco3_c -0.59209 LC_nh4_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_succ_c +0.59209 LC_ump_c - DG_LMPD_70_ump_c + DGo_LMPD_70_ump_c = 0
FU_LMPD_70_ump_c: DG_LMPD_70_ump_c +1000 FU_LMPD_70_ump_c < 1000
BU_LMPD_70_ump_c: - DG_LMPD_70_ump_c +1000 BU_LMPD_70_ump_c < 1000
SU_LMPD_70_ump_c: FU_LMPD_70_ump_c + BU_LMPD_70_ump_c < 1
UF_LMPD_70_ump_c: F_LMPD_70_ump_c -1000 FU_LMPD_70_ump_c < 0
UR_LMPD_70_ump_c: R_LMPD_70_ump_c -1000 BU_LMPD_70_ump_c < 0
G_LMPD_71_ump_c: +1.7763 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_hco3_c -0.59209 LC_nh4_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_ump_c -0.59209 LC_q6_m +0.59209 LC_q6h2_m - DG_LMPD_71_ump_c + DGo_LMPD_71_ump_c = 0
FU_LMPD_71_ump_c: DG_LMPD_71_ump_c +1000 FU_LMPD_71_ump_c < 1000
BU_LMPD_71_ump_c: - DG_LMPD_71_ump_c +1000 BU_LMPD_71_ump_c < 1000
SU_LMPD_71_ump_c: FU_LMPD_71_ump_c + BU_LMPD_71_ump_c < 1
UF_LMPD_71_ump_c: F_LMPD_71_ump_c -1000 FU_LMPD_71_ump_c < 0
UR_LMPD_71_ump_c: R_LMPD_71_ump_c -1000 BU_LMPD_71_ump_c < 0
G_LMPD_72_val-L_c: +0.59209 LC_akg_c -0.59209 LC_glu-L_c +0.59209 LC_val-L_c +0.59209 LC_co2_m +0.59209 LC_nadp_m -0.59209 LC_nadph_m -1.1842 LC_pyr_m - DG_LMPD_72_val-L_c + DGo_LMPD_72_val-L_c = 0
FU_LMPD_72_val-L_c: DG_LMPD_72_val-L_c +1000 FU_LMPD_72_val-L_c < 1000
BU_LMPD_72_val-L_c: - DG_LMPD_72_val-L_c +1000 BU_LMPD_72_val-L_c < 1000
SU_LMPD_72_val-L_c: FU_LMPD_72_val-L_c + BU_LMPD_72_val-L_c < 1
UF_LMPD_72_val-L_c: F_LMPD_72_val-L_c -1000 FU_LMPD_72_val-L_c < 0
UR_LMPD_72_val-L_c: R_LMPD_72_val-L_c -1000 BU_LMPD_72_val-L_c < 0
G_LMPD_73_zymst_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +4.7367 LC_co2_c +10.6576 LC_coa_c +0.59209 LC_for_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +13.026 LC_nadp_c -13.026 LC_nadph_c -5.3288 LC_o2_c +3.5525 LC_pi_c +3.5525 LC_ppi_c +0.59209 LC_zymst_c +0.59209 LC_nadp_r -0.59209 LC_nadph_r -0.59209 LC_o2_r - DG_LMPD_73_zymst_c + DGo_LMPD_73_zymst_c = 0
FU_LMPD_73_zymst_c: DG_LMPD_73_zymst_c +1000 FU_LMPD_73_zymst_c < 1000
BU_LMPD_73_zymst_c: - DG_LMPD_73_zymst_c +1000 BU_LMPD_73_zymst_c < 1000
SU_LMPD_73_zymst_c: FU_LMPD_73_zymst_c + BU_LMPD_73_zymst_c < 1
UF_LMPD_73_zymst_c: F_LMPD_73_zymst_c -1000 FU_LMPD_73_zymst_c < 0
UR_LMPD_73_zymst_c: R_LMPD_73_zymst_c -1000 BU_LMPD_73_zymst_c < 0
G_LMPD_74_zymst_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +4.7367 LC_co2_c +10.6576 LC_coa_c +0.59209 LC_for_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +14.8022 LC_nadp_c -14.8022 LC_nadph_c -5.3288 LC_o2_c +3.5525 LC_pi_c +3.5525 LC_ppi_c +0.59209 LC_zymst_c +0.59209 LC_nadp_r -0.59209 LC_nadph_r -0.59209 LC_o2_r - DG_LMPD_74_zymst_c + DGo_LMPD_74_zymst_c = 0
FU_LMPD_74_zymst_c: DG_LMPD_74_zymst_c +1000 FU_LMPD_74_zymst_c < 1000
BU_LMPD_74_zymst_c: - DG_LMPD_74_zymst_c +1000 BU_LMPD_74_zymst_c < 1000
SU_LMPD_74_zymst_c: FU_LMPD_74_zymst_c + BU_LMPD_74_zymst_c < 1
UF_LMPD_74_zymst_c: F_LMPD_74_zymst_c -1000 FU_LMPD_74_zymst_c < 0
UR_LMPD_74_zymst_c: R_LMPD_74_zymst_c -1000 BU_LMPD_74_zymst_c < 0
Bounds
 0  <= F_biomass_SC5_notrace  <= 100
 0  <= R_biomass_SC5_notrace  <= 0
 0  <= F_EX_4abut_e  <= 100
 0  <= R_EX_4abut_e  <= 0
 0  <= F_EX_ac_e  <= 100
 0  <= R_EX_ac_e  <= 0
 0  <= F_EX_acald_e  <= 100
 0  <= R_EX_acald_e  <= 0
 0  <= F_EX_akg_e  <= 100
 0  <= R_EX_akg_e  <= 0
 0  <= F_EX_ala-L_e  <= 100
 0  <= R_EX_ala-L_e  <= 0
 0  <= F_EX_arg-L_e  <= 100
 0  <= R_EX_arg-L_e  <= 0
 0  <= F_EX_asn-L_e  <= 100
 0  <= R_EX_asn-L_e  <= 0
 0  <= F_EX_asp-L_e  <= 100
 0  <= R_EX_asp-L_e  <= 0
 0  <= F_EX_cit_e  <= 100
 0  <= R_EX_cit_e  <= 0
 0  <= F_EX_co2_e  <= 100
 0  <= R_EX_co2_e  <= 0
 0  <= F_EX_cys-L_e  <= 100
 0  <= R_EX_cys-L_e  <= 0
 0  <= F_EX_ergst_e  <= 100
 0  <= R_EX_ergst_e  <= 0
 0  <= F_EX_etoh_e  <= 100
 0  <= R_EX_etoh_e  <= 0
 0  <= F_EX_for_e  <= 100
 0  <= R_EX_for_e  <= 0
 0  <= F_EX_fum_e  <= 100
 0  <= R_EX_fum_e  <= 0
 0  <= F_EX_glc_e  <= 100
 0  <= R_EX_glc_e  <= 10
 0  <= F_EX_gln-L_e  <= 100
 0  <= R_EX_gln-L_e  <= 0
 0  <= F_EX_glu-L_e  <= 100
 0  <= R_EX_glu-L_e  <= 0
 0  <= F_EX_gly_e  <= 100
 0  <= R_EX_gly_e  <= 0
 0  <= F_EX_h_e  <= 100
 0  <= R_EX_h_e  <= 100
 0  <= F_EX_h2o_e  <= 100
 0  <= R_EX_h2o_e  <= 100
 0  <= F_EX_his-L_e  <= 100
 0  <= R_EX_his-L_e  <= 0
 0  <= F_EX_ile-L_e  <= 100
 0  <= R_EX_ile-L_e  <= 0
 0  <= F_EX_leu-L_e  <= 100
 0  <= R_EX_leu-L_e  <= 0
 0  <= F_EX_lys-L_e  <= 100
 0  <= R_EX_lys-L_e  <= 0
 0  <= F_EX_mal-L_e  <= 100
 0  <= R_EX_mal-L_e  <= 0
 0  <= F_EX_met-L_e  <= 100
 0  <= R_EX_met-L_e  <= 0
 0  <= F_EX_nh4_e  <= 100
 0  <= R_EX_nh4_e  <= 100
 0  <= F_EX_o2_e  <= 100
 0  <= R_EX_o2_e  <= 2
 0  <= F_EX_oaa_e  <= 100
 0  <= R_EX_oaa_e  <= 0
 0  <= F_EX_phe-L_e  <= 100
 0  <= R_EX_phe-L_e  <= 0
 0  <= F_EX_pi_e  <= 100
 0  <= R_EX_pi_e  <= 100
 0  <= F_EX_pro-L_e  <= 100
 0  <= R_EX_pro-L_e  <= 0
 0  <= F_EX_pyr_e  <= 100
 0  <= R_EX_pyr_e  <= 0
 0  <= F_EX_ser-L_e  <= 100
 0  <= R_EX_ser-L_e  <= 0
 0  <= F_EX_so4_e  <= 100
 0  <= R_EX_so4_e  <= 100
 0  <= F_EX_succ_e  <= 100
 0  <= R_EX_succ_e  <= 0
 0  <= F_EX_thr-L_e  <= 100
 0  <= R_EX_thr-L_e  <= 0
 0  <= F_EX_trp-L_e  <= 100
 0  <= R_EX_trp-L_e  <= 0
 0  <= F_EX_tyr-L_e  <= 100
 0  <= R_EX_tyr-L_e  <= 0
 0  <= F_EX_val-L_e  <= 100
 0  <= R_EX_val-L_e  <= 0
 0  <= F_EX_zymst_e  <= 100
 0  <= R_EX_zymst_e  <= 0
 0  <= F_6PGLter  <= 100
 0  <= R_6PGLter  <= 100
 0  <= F_ABUTt2r  <= 100
 0  <= R_ABUTt2r  <= 100
 0  <= F_ACALDt  <= 100
 0  <= R_ACALDt  <= 100
 0  <= F_ACALDtm  <= 100
 0  <= R_ACALDtm  <= 100
 0  <= F_ACOAHim  <= 100
 0  <= R_ACOAHim  <= 0
 0  <= F_ACONT  <= 100
 0  <= R_ACONT  <= 100
 0  <= F_ACONTm  <= 100
 0  <= R_ACONTm  <= 100
 0  <= F_ACS  <= 100
 0  <= R_ACS  <= 0
 0  <= F_ACSm  <= 100
 0  <= R_ACSm  <= 0
 0  <= F_ACt2r  <= 100
 0  <= R_ACt2r  <= 100
 0  <= F_ACtm  <= 100
 0  <= R_ACtm  <= 100
 0  <= F_ACtr  <= 100
 0  <= R_ACtr  <= 100
 0  <= F_ADK1  <= 100
 0  <= R_ADK1  <= 100
 0  <= F_ADK1m  <= 100
 0  <= R_ADK1m  <= 100
 0  <= F_AKGDam  <= 100
 0  <= R_AKGDam  <= 100
 0  <= F_AKGDbm  <= 100
 0  <= R_AKGDbm  <= 0
 0  <= F_AKGMAL  <= 100
 0  <= R_AKGMAL  <= 100
 0  <= F_AKGt2r  <= 100
 0  <= R_AKGt2r  <= 100
 0  <= F_ALATA_L  <= 100
 0  <= R_ALATA_L  <= 100
 0  <= F_ALATA_Lm  <= 100
 0  <= R_ALATA_Lm  <= 100
 0  <= F_ALAt2r  <= 100
 0  <= R_ALAt2r  <= 100
 0  <= F_ALAtmi  <= 100
 0  <= R_ALAtmi  <= 0
 0  <= F_ALCD2if  <= 100
 0  <= R_ALCD2if  <= 0
 0  <= F_ALCD2ir  <= 100
 0  <= R_ALCD2ir  <= 0
 0  <= F_ALCD2irm  <= 100
 0  <= R_ALCD2irm  <= 0
 0  <= F_ALCD2x  <= 100
 0  <= R_ALCD2x  <= 100
 0  <= F_ALDD2xm  <= 100
 0  <= R_ALDD2xm  <= 0
 0  <= F_ALDD2y  <= 100
 0  <= R_ALDD2y  <= 0
 0  <= F_ALDD2ym  <= 100
 0  <= R_ALDD2ym  <= 0
 0  <= F_ARGt2r  <= 100
 0  <= R_ARGt2r  <= 100
 0  <= F_ASNt2r  <= 100
 0  <= R_ASNt2r  <= 100
 0  <= F_ASPGLU2m  <= 100
 0  <= R_ASPGLU2m  <= 0
 0  <= F_ASPTA  <= 100
 0  <= R_ASPTA  <= 100
 0  <= F_ASPTAm  <= 100
 0  <= R_ASPTAm  <= 100
 0  <= F_ASPt2m  <= 100
 0  <= R_ASPt2m  <= 100
 0  <= F_ASPt2r  <= 100
 0  <= R_ASPt2r  <= 100
 0  <= F_ATPH1  <= 100
 0  <= R_ATPH1  <= 0
 1  <= F_ATPM  <= 1
 0  <= R_ATPM  <= 0
 0  <= F_ATPS  <= 100
 0  <= R_ATPS  <= 0
 0  <= F_ATPS3m  <= 100
 0  <= R_ATPS3m  <= 0
 0  <= F_ATPtm-H  <= 100
 0  <= R_ATPtm-H  <= 0
 0  <= F_CITt2r  <= 100
 0  <= R_CITt2r  <= 100
 0  <= F_CITtam  <= 100
 0  <= R_CITtam  <= 100
 0  <= F_CITtcm  <= 100
 0  <= R_CITtcm  <= 100
 0  <= F_CO2t  <= 100
 0  <= R_CO2t  <= 100
 0  <= F_CO2tm  <= 100
 0  <= R_CO2tm  <= 100
 0  <= F_CSm  <= 100
 0  <= R_CSm  <= 0
 0  <= F_CYOOm  <= 100
 0  <= R_CYOOm  <= 0
 0  <= F_CYOR_u6m  <= 100
 0  <= R_CYOR_u6m  <= 0
 0  <= F_CYSt2r  <= 100
 0  <= R_CYSt2r  <= 100
 0  <= F_D-LACt2m  <= 100
 0  <= R_D-LACt2m  <= 100
 0  <= F_D-LACtm  <= 100
 0  <= R_D-LACtm  <= 100
 0  <= F_DB4PS  <= 100
 0  <= R_DB4PS  <= 0
 0  <= F_ENO  <= 100
 0  <= R_ENO  <= 100
 0  <= F_ERGSTt  <= 100
 0  <= R_ERGSTt  <= 100
 0  <= F_ERGSTter  <= 100
 0  <= R_ERGSTter  <= 100
 0  <= F_ETOHt  <= 100
 0  <= R_ETOHt  <= 100
 0  <= F_ETOHtm  <= 100
 0  <= R_ETOHtm  <= 100
 0  <= F_FBA  <= 100
 0  <= R_FBA  <= 100
 0  <= F_FBA3  <= 100
 0  <= R_FBA3  <= 100
 0  <= F_FBP  <= 100
 0  <= R_FBP  <= 0
 0  <= F_FDH  <= 100
 0  <= R_FDH  <= 0
 0  <= F_FDNG  <= 100
 0  <= R_FDNG  <= 0
 0  <= F_FORt  <= 100
 0  <= R_FORt  <= 100
 0  <= F_FRDcm  <= 100
 0  <= R_FRDcm  <= 0
 0  <= F_FRDm  <= 100
 0  <= R_FRDm  <= 0
 0  <= F_FUM  <= 100
 0  <= R_FUM  <= 100
 0  <= F_FUMm  <= 100
 0  <= R_FUMm  <= 100
 0  <= F_FUMt2r  <= 100
 0  <= R_FUMt2r  <= 100
 0  <= F_G6PDH2  <= 100
 0  <= R_G6PDH2  <= 0
 0  <= F_G6PDH2er  <= 100
 0  <= R_G6PDH2er  <= 0
 0  <= F_G6PI  <= 100
 0  <= R_G6PI  <= 100
 0  <= F_G6PI3  <= 100
 0  <= R_G6PI3  <= 100
 0  <= F_G6Pter  <= 100
 0  <= R_G6Pter  <= 100
 0  <= F_GAPD  <= 100
 0  <= R_GAPD  <= 100
 0  <= F_GCC2cm  <= 100
 0  <= R_GCC2cm  <= 100
 0  <= F_GLCt1  <= 100
 0  <= R_GLCt1  <= 0
 0  <= F_GLNt2r  <= 100
 0  <= R_GLNt2r  <= 100
 0  <= F_GLUDC  <= 100
 0  <= R_GLUDC  <= 0
 0  <= F_GLUDxi  <= 100
 0  <= R_GLUDxi  <= 0
 0  <= F_GLUDyi  <= 100
 0  <= R_GLUDyi  <= 0
 0  <= F_GLUK  <= 100
 0  <= R_GLUK  <= 0
 0  <= F_GLUt2r  <= 100
 0  <= R_GLUt2r  <= 100
 0  <= F_GLUt7m  <= 100
 0  <= R_GLUt7m  <= 0
 0  <= F_GLYt2r  <= 100
 0  <= R_GLYt2r  <= 100
 0  <= F_GND  <= 100
 0  <= R_GND  <= 0
 0  <= F_H2Ot  <= 100
 0  <= R_H2Ot  <= 100
 0  <= F_H2Oter  <= 100
 0  <= R_H2Oter  <= 100
 0  <= F_H2Otm  <= 100
 0  <= R_H2Otm  <= 100
 0  <= F_HCITSm  <= 100
 0  <= R_HCITSm  <= 0
 0  <= F_HCO3E  <= 100
 0  <= R_HCO3E  <= 100
 0  <= F_HEX1  <= 100
 0  <= R_HEX1  <= 0
 0  <= F_HISt2r  <= 100
 0  <= R_HISt2r  <= 100
 0  <= F_ICDHxm  <= 100
 0  <= R_ICDHxm  <= 0
 0  <= F_ICDHy  <= 100
 0  <= R_ICDHy  <= 0
 0  <= F_ICDHym  <= 100
 0  <= R_ICDHym  <= 0
 0  <= F_ILEt2r  <= 100
 0  <= R_ILEt2r  <= 100
 0  <= F_ILEtmi  <= 100
 0  <= R_ILEtmi  <= 0
 0  <= F_LEUt2r  <= 100
 0  <= R_LEUt2r  <= 100
 0  <= F_LYSt2r  <= 100
 0  <= R_LYSt2r  <= 100
 0  <= F_MALt2r  <= 100
 0  <= R_MALt2r  <= 100
 0  <= F_MALtm  <= 100
 0  <= R_MALtm  <= 100
 0  <= F_MAN6PI  <= 100
 0  <= R_MAN6PI  <= 100
 0  <= F_MANNANter  <= 100
 0  <= R_MANNANter  <= 100
 0  <= F_MDH  <= 100
 0  <= R_MDH  <= 100
 0  <= F_MDHm  <= 100
 0  <= R_MDHm  <= 100
 0  <= F_ME1m  <= 100
 0  <= R_ME1m  <= 0
 0  <= F_ME2m  <= 100
 0  <= R_ME2m  <= 0
 0  <= F_METt2r  <= 100
 0  <= R_METt2r  <= 100
 0  <= F_NADH2-u6cm  <= 100
 0  <= R_NADH2-u6cm  <= 0
 0  <= F_NADH2-u6m  <= 100
 0  <= R_NADH2-u6m  <= 0
 0  <= F_NADK  <= 100
 0  <= R_NADK  <= 0
 0  <= F_NADKm  <= 100
 0  <= R_NADKm  <= 0
 0  <= F_NADPPPS  <= 100
 0  <= R_NADPPPS  <= 0
 0  <= F_NADPPPSm  <= 100
 0  <= R_NADPPPSm  <= 0
 0  <= F_NDP1  <= 100
 0  <= R_NDP1  <= 0
 0  <= F_NH4t  <= 100
 0  <= R_NH4t  <= 100
 0  <= F_NH4tm  <= 100
 0  <= R_NH4tm  <= 100
 0  <= F_O2t  <= 100
 0  <= R_O2t  <= 100
 0  <= F_O2ter  <= 100
 0  <= R_O2ter  <= 100
 0  <= F_O2tm  <= 100
 0  <= R_O2tm  <= 100
 0  <= F_OAAt  <= 100
 0  <= R_OAAt  <= 100
 0  <= F_OAAt2m  <= 100
 0  <= R_OAAt2m  <= 100
 0  <= F_PC  <= 100
 0  <= R_PC  <= 0
 0  <= F_PDHcm  <= 100
 0  <= R_PDHcm  <= 0
 0  <= F_PDHm  <= 100
 0  <= R_PDHm  <= 0
 0  <= F_PEtm_SC  <= 100
 0  <= R_PEtm_SC  <= 100
 0  <= F_PFK  <= 100
 0  <= R_PFK  <= 0
 0  <= F_PFK_3  <= 100
 0  <= R_PFK_3  <= 0
 0  <= F_PGI  <= 100
 0  <= R_PGI  <= 100
 0  <= F_PGK  <= 100
 0  <= R_PGK  <= 100
 0  <= F_PGL  <= 100
 0  <= R_PGL  <= 0
 0  <= F_PGM  <= 100
 0  <= R_PGM  <= 100
 0  <= F_PGMT  <= 100
 0  <= R_PGMT  <= 100
 0  <= F_PHEt2r  <= 100
 0  <= R_PHEt2r  <= 100
 0  <= F_PIt2m  <= 100
 0  <= R_PIt2m  <= 100
 0  <= F_PIt2r  <= 100
 0  <= R_PIt2r  <= 100
 0  <= F_PPA  <= 100
 0  <= R_PPA  <= 0
 0  <= F_PPAm  <= 100
 0  <= R_PPAm  <= 0
 0  <= F_PPCK  <= 100
 0  <= R_PPCK  <= 0
 0  <= F_PROt2r  <= 100
 0  <= R_PROt2r  <= 100
 0  <= F_PStm_SC  <= 100
 0  <= R_PStm_SC  <= 100
 0  <= F_PYK  <= 100
 0  <= R_PYK  <= 0
 0  <= F_PYRDC  <= 100
 0  <= R_PYRDC  <= 0
 0  <= F_PYRt  <= 100
 0  <= R_PYRt  <= 0
 0  <= F_PYRt2m  <= 100
 0  <= R_PYRt2m  <= 100
 0  <= F_RPE  <= 100
 0  <= R_RPE  <= 100
 0  <= F_RPI  <= 100
 0  <= R_RPI  <= 100
 0  <= F_SERt2r  <= 100
 0  <= R_SERt2r  <= 100
 0  <= F_SO4ti  <= 100
 0  <= R_SO4ti  <= 0
 0  <= F_SUCCt2r  <= 100
 0  <= R_SUCCt2r  <= 100
 0  <= F_SUCCtm  <= 100
 0  <= R_SUCCtm  <= 0
 0  <= F_SUCD1m  <= 100
 0  <= R_SUCD1m  <= 100
 0  <= F_SUCD2_u6m  <= 100
 0  <= R_SUCD2_u6m  <= 100
 0  <= F_SUCD3-u6m  <= 100
 0  <= R_SUCD3-u6m  <= 100
 0  <= F_SUCFUMtm  <= 100
 0  <= R_SUCFUMtm  <= 0
 0  <= F_SUCOASm  <= 100
 0  <= R_SUCOASm  <= 100
 0  <= F_TALA  <= 100
 0  <= R_TALA  <= 100
 0  <= F_THRt2m  <= 100
 0  <= R_THRt2m  <= 100
 0  <= F_THRt2r  <= 100
 0  <= R_THRt2r  <= 100
 0  <= F_TKT1  <= 100
 0  <= R_TKT1  <= 100
 0  <= F_TKT2  <= 100
 0  <= R_TKT2  <= 100
 0  <= F_TPI  <= 100
 0  <= R_TPI  <= 100
 0  <= F_TRPt2r  <= 100
 0  <= R_TRPt2r  <= 100
 0  <= F_TYRt2r  <= 100
 0  <= R_TYRt2r  <= 100
 0  <= F_VALt2r  <= 100
 0  <= R_VALt2r  <= 100
 0  <= F_ZYMSTt  <= 100
 0  <= R_ZYMSTt  <= 100
 0  <= F_LMPD_1_13BDglcn_c  <= 100
 0  <= R_LMPD_1_13BDglcn_c  <= 0
 0  <= F_LMPD_2_amp_c  <= 100
 0  <= R_LMPD_2_amp_c  <= 0
 0  <= F_LMPD_3_amp_c  <= 100
 0  <= R_LMPD_3_amp_c  <= 0
 0  <= F_LMPD_4_amp_c  <= 100
 0  <= R_LMPD_4_amp_c  <= 0
 0  <= F_LMPD_5_amp_c  <= 100
 0  <= R_LMPD_5_amp_c  <= 0
 0  <= F_LMPD_6_arg-L_c  <= 100
 0  <= R_LMPD_6_arg-L_c  <= 0
 0  <= F_LMPD_7_asn-L_c  <= 100
 0  <= R_LMPD_7_asn-L_c  <= 0
 0  <= F_LMPD_8_cmp_c  <= 100
 0  <= R_LMPD_8_cmp_c  <= 0
 0  <= F_LMPD_9_cmp_c  <= 100
 0  <= R_LMPD_9_cmp_c  <= 0
 0  <= F_LMPD_10_cmp_c  <= 100
 0  <= R_LMPD_10_cmp_c  <= 0
 0  <= F_LMPD_11_cmp_c  <= 100
 0  <= R_LMPD_11_cmp_c  <= 0
 0  <= F_LMPD_12_cys-L_c  <= 100
 0  <= R_LMPD_12_cys-L_c  <= 0
 0  <= F_LMPD_13_damp_c  <= 100
 0  <= R_LMPD_13_damp_c  <= 0
 0  <= F_LMPD_14_damp_c  <= 100
 0  <= R_LMPD_14_damp_c  <= 0
 0  <= F_LMPD_15_damp_c  <= 100
 0  <= R_LMPD_15_damp_c  <= 0
 0  <= F_LMPD_16_damp_c  <= 100
 0  <= R_LMPD_16_damp_c  <= 0
 0  <= F_LMPD_17_dcmp_c  <= 100
 0  <= R_LMPD_17_dcmp_c  <= 0
 0  <= F_LMPD_18_dcmp_c  <= 100
 0  <= R_LMPD_18_dcmp_c  <= 0
 0  <= F_LMPD_19_dgmp_c  <= 100
 0  <= R_LMPD_19_dgmp_c  <= 0
 0  <= F_LMPD_20_dgmp_c  <= 100
 0  <= R_LMPD_20_dgmp_c  <= 0
 0  <= F_LMPD_21_dgmp_c  <= 100
 0  <= R_LMPD_21_dgmp_c  <= 0
 0  <= F_LMPD_22_dgmp_c  <= 100
 0  <= R_LMPD_22_dgmp_c  <= 0
 0  <= F_LMPD_23_dtmp_c  <= 100
 0  <= R_LMPD_23_dtmp_c  <= 0
 0  <= F_LMPD_24_dtmp_c  <= 100
 0  <= R_LMPD_24_dtmp_c  <= 0
 0  <= F_LMPD_25_ergst_c  <= 100
 0  <= R_LMPD_25_ergst_c  <= 0
 0  <= F_LMPD_26_ergst_c  <= 100
 0  <= R_LMPD_26_ergst_c  <= 0
 0  <= F_LMPD_27_gln-L_c  <= 100
 0  <= R_LMPD_27_gln-L_c  <= 0
 0  <= F_LMPD_28_gly_c  <= 100
 0  <= R_LMPD_28_gly_c  <= 0
 0  <= F_LMPD_29_glycogen_c  <= 100
 0  <= R_LMPD_29_glycogen_c  <= 0
 0  <= F_LMPD_30_gmp_c  <= 100
 0  <= R_LMPD_30_gmp_c  <= 0
 0  <= F_LMPD_31_gmp_c  <= 100
 0  <= R_LMPD_31_gmp_c  <= 0
 0  <= F_LMPD_32_his-L_c  <= 100
 0  <= R_LMPD_32_his-L_c  <= 0
 0  <= F_LMPD_33_his-L_c  <= 100
 0  <= R_LMPD_33_his-L_c  <= 0
 0  <= F_LMPD_34_his-L_c  <= 100
 0  <= R_LMPD_34_his-L_c  <= 0
 0  <= F_LMPD_35_his-L_c  <= 100
 0  <= R_LMPD_35_his-L_c  <= 0
 0  <= F_LMPD_36_ile-L_c  <= 100
 0  <= R_LMPD_36_ile-L_c  <= 0
 0  <= F_LMPD_37_ile-L_c  <= 100
 0  <= R_LMPD_37_ile-L_c  <= 0
 0  <= F_LMPD_38_leu-L_c  <= 100
 0  <= R_LMPD_38_leu-L_c  <= 0
 0  <= F_LMPD_39_leu-L_c  <= 100
 0  <= R_LMPD_39_leu-L_c  <= 0
 0  <= F_LMPD_40_lys-L_c  <= 100
 0  <= R_LMPD_40_lys-L_c  <= 0
 0  <= F_LMPD_41_lys-L_c  <= 100
 0  <= R_LMPD_41_lys-L_c  <= 0
 0  <= F_LMPD_42_mannan_c  <= 100
 0  <= R_LMPD_42_mannan_c  <= 0
 0  <= F_LMPD_43_mannan_c  <= 100
 0  <= R_LMPD_43_mannan_c  <= 0
 0  <= F_LMPD_44_met-L_c  <= 100
 0  <= R_LMPD_44_met-L_c  <= 0
 0  <= F_LMPD_45_met-L_c  <= 100
 0  <= R_LMPD_45_met-L_c  <= 0
 0  <= F_LMPD_47_pa_SC_c  <= 100
 0  <= R_LMPD_47_pa_SC_c  <= 0
 0  <= F_LMPD_48_pc_SC_c  <= 100
 0  <= R_LMPD_48_pc_SC_c  <= 0
 0  <= F_LMPD_49_pc_SC_c  <= 100
 0  <= R_LMPD_49_pc_SC_c  <= 0
 0  <= F_LMPD_51_pe_SC_c  <= 100
 0  <= R_LMPD_51_pe_SC_c  <= 0
 0  <= F_LMPD_52_phe-L_c  <= 100
 0  <= R_LMPD_52_phe-L_c  <= 0
 0  <= F_LMPD_53_pro-L_c  <= 100
 0  <= R_LMPD_53_pro-L_c  <= 0
 0  <= F_LMPD_54_pro-L_c  <= 100
 0  <= R_LMPD_54_pro-L_c  <= 0
 0  <= F_LMPD_56_ps_SC_c  <= 100
 0  <= R_LMPD_56_ps_SC_c  <= 0
 0  <= F_LMPD_57_ptd1ino_SC_c  <= 100
 0  <= R_LMPD_57_ptd1ino_SC_c  <= 0
 0  <= F_LMPD_58_ptd1ino_SC_c  <= 100
 0  <= R_LMPD_58_ptd1ino_SC_c  <= 0
 0  <= F_LMPD_59_ribflv_c  <= 100
 0  <= R_LMPD_59_ribflv_c  <= 0
 0  <= F_LMPD_60_ribflv_c  <= 100
 0  <= R_LMPD_60_ribflv_c  <= 0
 0  <= F_LMPD_61_ser-L_c  <= 100
 0  <= R_LMPD_61_ser-L_c  <= 0
 0  <= F_LMPD_62_thr-L_c  <= 100
 0  <= R_LMPD_62_thr-L_c  <= 0
 0  <= F_LMPD_63_thr-L_c  <= 100
 0  <= R_LMPD_63_thr-L_c  <= 0
 0  <= F_LMPD_64_tre_c  <= 100
 0  <= R_LMPD_64_tre_c  <= 0
 0  <= F_LMPD_65_triglyc_SC_c  <= 100
 0  <= R_LMPD_65_triglyc_SC_c  <= 0
 0  <= F_LMPD_66_triglyc_SC_c  <= 100
 0  <= R_LMPD_66_triglyc_SC_c  <= 0
 0  <= F_LMPD_67_trp-L_c  <= 100
 0  <= R_LMPD_67_trp-L_c  <= 0
 0  <= F_LMPD_68_tyr-L_c  <= 100
 0  <= R_LMPD_68_tyr-L_c  <= 0
 0  <= F_LMPD_69_tyr-L_c  <= 100
 0  <= R_LMPD_69_tyr-L_c  <= 0
 0  <= F_LMPD_70_ump_c  <= 100
 0  <= R_LMPD_70_ump_c  <= 0
 0  <= F_LMPD_71_ump_c  <= 100
 0  <= R_LMPD_71_ump_c  <= 0
 0  <= F_LMPD_72_val-L_c  <= 100
 0  <= R_LMPD_72_val-L_c  <= 0
 0  <= F_LMPD_73_zymst_c  <= 100
 0  <= R_LMPD_73_zymst_c  <= 0
 0  <= F_LMPD_74_zymst_c  <= 100
 0  <= R_LMPD_74_zymst_c  <= 0
 -1.611810e+01  <= LC_13BDglcn_c  <= -4.605170e+00
 -1.611810e+01  <= LC_13dpg_c  <= -4.605170e+00
 -1.611810e+01  <= LC_2pg_c  <= -4.605170e+00
 -1.611810e+01  <= LC_3pg_c  <= -4.605170e+00
 -1.611810e+01  <= LC_4abut_c  <= -4.605170e+00
 -1.611810e+01  <= LC_6pgc_c  <= -4.605170e+00
 -1.611810e+01  <= LC_6pgl_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ac_c  <= -4.605170e+00
 -1.611810e+01  <= LC_acald_c  <= -4.605170e+00
 -1.611810e+01  <= LC_accoa_c  <= -4.605170e+00
 -1.611810e+01  <= LC_adp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_akg_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ala-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_amp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_arg-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_asn-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_asp-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_atp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_cit_c  <= -4.605170e+00
 -1.611810e+01  <= LC_cmp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_co2_c  <= -4.605170e+00
 -1.611810e+01  <= LC_coa_c  <= -4.605170e+00
 -1.611810e+01  <= LC_cys-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_damp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_db4p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_dcmp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_dgmp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_dhap_c  <= -4.605170e+00
 -1.611810e+01  <= LC_dtmp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_e4p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ergst_c  <= -4.605170e+00
 -1.611810e+01  <= LC_etoh_c  <= -4.605170e+00
 -1.611810e+01  <= LC_f6p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_fdp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_for_c  <= -4.605170e+00
 -1.611810e+01  <= LC_fum_c  <= -4.605170e+00
 -1.611810e+01  <= LC_g1p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_g3p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_g6p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_g6p-B_c  <= -4.605170e+00
 -1.611810e+01  <= LC_glc-D_c  <= -4.605170e+00
 -1.611810e+01  <= LC_gln-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_glu-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_gly_c  <= -4.605170e+00
 -1.611810e+01  <= LC_glycogen_c  <= -4.605170e+00
 -1.611810e+01  <= LC_gmp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_h_c  <= -1.611810e+01
 0  <= LC_h2o_c  <= 0
 -1.611810e+01  <= LC_hco3_c  <= -4.605170e+00
 -1.611810e+01  <= LC_his-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_icit_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ile-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_lac-D_c  <= -4.605170e+00
 -1.611810e+01  <= LC_leu-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_lys-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_mal-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_man6p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_mannan_c  <= -4.605170e+00
 -1.611810e+01  <= LC_met-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_nad_c  <= -4.605170e+00
 -1.611810e+01  <= LC_nadh_c  <= -4.605170e+00
 -1.611810e+01  <= LC_nadp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_nadph_c  <= -4.605170e+00
 -1.611810e+01  <= LC_nh4_c  <= -4.605170e+00
 -1.611810e+01  <= LC_o2_c  <= -4.605170e+00
 -1.611810e+01  <= LC_oaa_c  <= -4.605170e+00
 -1.611810e+01  <= LC_pep_c  <= -4.605170e+00
 -1.611810e+01  <= LC_phe-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_pi_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ppi_c  <= -4.605170e+00
 -1.611810e+01  <= LC_pro-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_pyr_c  <= -4.605170e+00
 -1.611810e+01  <= LC_r5p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ribflv_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ru5p-D_c  <= -4.605170e+00
 -1.611810e+01  <= LC_s17bp_c  <= -4.605170e+00
 -1.611810e+01  <= LC_s7p_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ser-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_so4_c  <= -4.605170e+00
 -1.611810e+01  <= LC_succ_c  <= -4.605170e+00
 -1.611810e+01  <= LC_thr-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_tre_c  <= -4.605170e+00
 -1.611810e+01  <= LC_trp-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_tyr-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_ump_c  <= -4.605170e+00
 -1.611810e+01  <= LC_val-L_c  <= -4.605170e+00
 -1.611810e+01  <= LC_xu5p-D_c  <= -4.605170e+00
 -1.611810e+01  <= LC_zymst_c  <= -4.605170e+00
 -1.611810e+01  <= LC_4abut_e  <= -4.605170e+00
 -1.611810e+01  <= LC_ac_e  <= -4.605170e+00
 -1.611810e+01  <= LC_acald_e  <= -4.605170e+00
 -1.611810e+01  <= LC_akg_e  <= -4.605170e+00
 -1.611810e+01  <= LC_ala-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_arg-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_asn-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_asp-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_cit_e  <= -4.605170e+00
 -1.611810e+01  <= LC_co2_e  <= -4.605170e+00
 -1.611810e+01  <= LC_cys-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_ergst_e  <= -4.605170e+00
 -1.611810e+01  <= LC_etoh_e  <= -4.605170e+00
 -1.611810e+01  <= LC_for_e  <= -4.605170e+00
 -1.611810e+01  <= LC_fum_e  <= -4.605170e+00
 -1.611810e+01  <= LC_glc-D_e  <= -4.605170e+00
 -1.611810e+01  <= LC_gln-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_glu-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_gly_e  <= -4.605170e+00
 -1.151293e+01  <= LC_h_e  <= -1.151293e+01
 0  <= LC_h2o_e  <= 0
 -1.611810e+01  <= LC_his-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_ile-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_leu-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_lys-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_mal-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_met-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_nh4_e  <= -4.605170e+00
 -1.611810e+01  <= LC_o2_e  <= -4.605170e+00
 -1.611810e+01  <= LC_oaa_e  <= -4.605170e+00
 -1.611810e+01  <= LC_phe-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_pi_e  <= -4.605170e+00
 -1.611810e+01  <= LC_pro-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_pyr_e  <= -4.605170e+00
 -1.611810e+01  <= LC_ser-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_so4_e  <= -4.605170e+00
 -1.611810e+01  <= LC_succ_e  <= -4.605170e+00
 -1.611810e+01  <= LC_thr-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_trp-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_tyr-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_val-L_e  <= -4.605170e+00
 -1.611810e+01  <= LC_zymst_e  <= -4.605170e+00
 -1.611810e+01  <= LC_ac_m  <= -4.605170e+00
 -1.611810e+01  <= LC_acald_m  <= -4.605170e+00
 -1.611810e+01  <= LC_accoa_m  <= -4.605170e+00
 -1.611810e+01  <= LC_adp_m  <= -4.605170e+00
 -1.611810e+01  <= LC_akg_m  <= -4.605170e+00
 -1.611810e+01  <= LC_ala-L_m  <= -4.605170e+00
 -1.611810e+01  <= LC_amp_m  <= -4.605170e+00
 -1.611810e+01  <= LC_asp-L_m  <= -4.605170e+00
 -1.611810e+01  <= LC_atp_m  <= -4.605170e+00
 -1.611810e+01  <= LC_cit_m  <= -4.605170e+00
 -1.611810e+01  <= LC_co2_m  <= -4.605170e+00
 -1.611810e+01  <= LC_coa_m  <= -4.605170e+00
 -1.611810e+01  <= LC_dhlam_m  <= -4.605170e+00
 -1.611810e+01  <= LC_etoh_m  <= -4.605170e+00
 -1.611810e+01  <= LC_fad_m  <= -4.605170e+00
 -1.611810e+01  <= LC_fadh2_m  <= -4.605170e+00
 -1.611810e+01  <= LC_ficytc_m  <= -4.605170e+00
 -1.611810e+01  <= LC_focytc_m  <= -4.605170e+00
 -1.611810e+01  <= LC_fum_m  <= -4.605170e+00
 -1.611810e+01  <= LC_glu-L_m  <= -4.605170e+00
 -1.726939e+01  <= LC_h_m  <= -1.726939e+01
 0  <= LC_h2o_m  <= 0
 -1.611810e+01  <= LC_hcit_m  <= -4.605170e+00
 -1.611810e+01  <= LC_icit_m  <= -4.605170e+00
 -1.611810e+01  <= LC_ile-L_m  <= -4.605170e+00
 -1.611810e+01  <= LC_lac-D_m  <= -4.605170e+00
 -1.611810e+01  <= LC_lpam_m  <= -4.605170e+00
 -1.611810e+01  <= LC_mal-L_m  <= -4.605170e+00
 -1.611810e+01  <= LC_nad_m  <= -4.605170e+00
 -1.611810e+01  <= LC_nadh_m  <= -4.605170e+00
 -1.611810e+01  <= LC_nadp_m  <= -4.605170e+00
 -1.611810e+01  <= LC_nadph_m  <= -4.605170e+00
 -1.611810e+01  <= LC_nh4_m  <= -4.605170e+00
 -1.611810e+01  <= LC_o2_m  <= -4.605170e+00
 -1.611810e+01  <= LC_oaa_m  <= -4.605170e+00
 -1.611810e+01  <= LC_pi_m  <= -4.605170e+00
 -1.611810e+01  <= LC_ppi_m  <= -4.605170e+00
 -1.611810e+01  <= LC_pyr_m  <= -4.605170e+00
 -1.611810e+01  <= LC_q6_m  <= -4.605170e+00
 -1.611810e+01  <= LC_q6h2_m  <= -4.605170e+00
 -1.611810e+01  <= LC_sdhlam_m  <= -4.605170e+00
 -1.611810e+01  <= LC_succ_m  <= -4.605170e+00
 -1.611810e+01  <= LC_succoa_m  <= -4.605170e+00
 -1.611810e+01  <= LC_thr-L_m  <= -4.605170e+00
 -1.611810e+01  <= LC_6pgl_r  <= -4.605170e+00
 -1.611810e+01  <= LC_ergst_r  <= -4.605170e+00
 -1.611810e+01  <= LC_g6p_r  <= -4.605170e+00
 -1.611810e+01  <= LC_h_r  <= -1.611810e+01
 0  <= LC_h2o_r  <= 0
 -1.611810e+01  <= LC_mannan_r  <= -4.605170e+00
 -1.611810e+01  <= LC_nadp_r  <= -4.605170e+00
 -1.611810e+01  <= LC_nadph_r  <= -4.605170e+00
 -1.611810e+01  <= LC_o2_r  <= -4.605170e+00
 0  <= FU_biomass_SC5_notrace  <= 1
 0  <= BU_biomass_SC5_notrace  <= 1
 0  <= FU_EX_4abut_e  <= 1
 0  <= BU_EX_4abut_e  <= 1
 0  <= FU_EX_ac_e  <= 1
 0  <= BU_EX_ac_e  <= 1
 0  <= FU_EX_acald_e  <= 1
 0  <= BU_EX_acald_e  <= 1
 0  <= FU_EX_akg_e  <= 1
 0  <= BU_EX_akg_e  <= 1
 0  <= FU_EX_ala-L_e  <= 1
 0  <= BU_EX_ala-L_e  <= 1
 0  <= FU_EX_arg-L_e  <= 1
 0  <= BU_EX_arg-L_e  <= 1
 0  <= FU_EX_asn-L_e  <= 1
 0  <= BU_EX_asn-L_e  <= 1
 0  <= FU_EX_asp-L_e  <= 1
 0  <= BU_EX_asp-L_e  <= 1
 0  <= FU_EX_cit_e  <= 1
 0  <= BU_EX_cit_e  <= 1
 0  <= FU_EX_co2_e  <= 1
 0  <= BU_EX_co2_e  <= 1
 0  <= FU_EX_cys-L_e  <= 1
 0  <= BU_EX_cys-L_e  <= 1
 0  <= FU_EX_ergst_e  <= 1
 0  <= BU_EX_ergst_e  <= 1
 0  <= FU_EX_etoh_e  <= 1
 0  <= BU_EX_etoh_e  <= 1
 0  <= FU_EX_for_e  <= 1
 0  <= BU_EX_for_e  <= 1
 0  <= FU_EX_fum_e  <= 1
 0  <= BU_EX_fum_e  <= 1
 0  <= FU_EX_glc_e  <= 1
 0  <= BU_EX_glc_e  <= 1
 0  <= FU_EX_gln-L_e  <= 1
 0  <= BU_EX_gln-L_e  <= 1
 0  <= FU_EX_glu-L_e  <= 1
 0  <= BU_EX_glu-L_e  <= 1
 0  <= FU_EX_gly_e  <= 1
 0  <= BU_EX_gly_e  <= 1
 0  <= FU_EX_h_e  <= 1
 0  <= BU_EX_h_e  <= 1
 0  <= FU_EX_h2o_e  <= 1
 0  <= BU_EX_h2o_e  <= 1
 0  <= FU_EX_his-L_e  <= 1
 0  <= BU_EX_his-L_e  <= 1
 0  <= FU_EX_ile-L_e  <= 1
 0  <= BU_EX_ile-L_e  <= 1
 0  <= FU_EX_leu-L_e  <= 1
 0  <= BU_EX_leu-L_e  <= 1
 0  <= FU_EX_lys-L_e  <= 1
 0  <= BU_EX_lys-L_e  <= 1
 0  <= FU_EX_mal-L_e  <= 1
 0  <= BU_EX_mal-L_e  <= 1
 0  <= FU_EX_met-L_e  <= 1
 0  <= BU_EX_met-L_e  <= 1
 0  <= FU_EX_nh4_e  <= 1
 0  <= BU_EX_nh4_e  <= 1
 0  <= FU_EX_o2_e  <= 1
 0  <= BU_EX_o2_e  <= 1
 0  <= FU_EX_oaa_e  <= 1
 0  <= BU_EX_oaa_e  <= 1
 0  <= FU_EX_phe-L_e  <= 1
 0  <= BU_EX_phe-L_e  <= 1
 0  <= FU_EX_pi_e  <= 1
 0  <= BU_EX_pi_e  <= 1
 0  <= FU_EX_pro-L_e  <= 1
 0  <= BU_EX_pro-L_e  <= 1
 0  <= FU_EX_pyr_e  <= 1
 0  <= BU_EX_pyr_e  <= 1
 0  <= FU_EX_ser-L_e  <= 1
 0  <= BU_EX_ser-L_e  <= 1
 0  <= FU_EX_so4_e  <= 1
 0  <= BU_EX_so4_e  <= 1
 0  <= FU_EX_succ_e  <= 1
 0  <= BU_EX_succ_e  <= 1
 0  <= FU_EX_thr-L_e  <= 1
 0  <= BU_EX_thr-L_e  <= 1
 0  <= FU_EX_trp-L_e  <= 1
 0  <= BU_EX_trp-L_e  <= 1
 0  <= FU_EX_tyr-L_e  <= 1
 0  <= BU_EX_tyr-L_e  <= 1
 0  <= FU_EX_val-L_e  <= 1
 0  <= BU_EX_val-L_e  <= 1
 0  <= FU_EX_zymst_e  <= 1
 0  <= BU_EX_zymst_e  <= 1
 -1000  <= DG_6PGLter  <= 1000
 -3.154129e+00  <= DGo_6PGLter  <= 8.458714e-01
 0  <= FU_6PGLter  <= 1
 0  <= BU_6PGLter  <= 1
 -1000  <= DG_ABUTt2r  <= 1000
 -7.310658e+00  <= DGo_ABUTt2r  <= -3.310658e+00
 0  <= FU_ABUTt2r  <= 1
 0  <= BU_ABUTt2r  <= 1
 -1000  <= DG_ACALDt  <= 1000
 -1.225929e+00  <= DGo_ACALDt  <= 2.774071e+00
 0  <= FU_ACALDt  <= 1
 0  <= BU_ACALDt  <= 1
 -1000  <= DG_ACALDtm  <= 1000
 -2.000000e+00  <= DGo_ACALDtm  <= 2.000000e+00
 0  <= FU_ACALDtm  <= 1
 0  <= BU_ACALDtm  <= 1
 -1000  <= DG_ACOAHim  <= 1000
 -7.677690e+00  <= DGo_ACOAHim  <= -5.064601e+00
 0  <= FU_ACOAHim  <= 1
 0  <= BU_ACOAHim  <= 1
 -1000  <= DG_ACONT  <= 1000
 9.834778e-01  <= DGo_ACONT  <= 1.855599e+00
 0  <= FU_ACONT  <= 1
 0  <= BU_ACONT  <= 1
 -1000  <= DG_ACONTm  <= 1000
 9.837901e-01  <= DGo_ACONTm  <= 1.855911e+00
 0  <= FU_ACONTm  <= 1
 0  <= BU_ACONTm  <= 1
 -1000  <= DG_ACS  <= 1000
 -4.709400e+00  <= DGo_ACS  <= -2.051134e+00
 0  <= FU_ACS  <= 1
 0  <= BU_ACS  <= 1
 -1000  <= DG_ACSm  <= 1000
 -5.051483e+00  <= DGo_ACSm  <= -2.393218e+00
 0  <= FU_ACSm  <= 1
 0  <= BU_ACSm  <= 1
 -1000  <= DG_ACt2r  <= 1000
 -4.077188e+00  <= DGo_ACt2r  <= -7.718806e-02
 0  <= FU_ACt2r  <= 1
 0  <= BU_ACt2r  <= 1
 -1000  <= DG_ACtm  <= 1000
 2.152162e+00  <= DGo_ACtm  <= 6.152162e+00
 0  <= FU_ACtm  <= 1
 0  <= BU_ACtm  <= 1
 -1000  <= DG_ACtr  <= 1000
 3.314429e-02  <= DGo_ACtr  <= 4.033144e+00
 0  <= FU_ACtr  <= 1
 0  <= BU_ACtr  <= 1
 -1000  <= DG_ADK1  <= 1000
 -1.568961e+00  <= DGo_ADK1  <= 2.431039e+00
 0  <= FU_ADK1  <= 1
 0  <= BU_ADK1  <= 1
 -1000  <= DG_ADK1m  <= 1000
 -1.570456e+00  <= DGo_ADK1m  <= 2.429544e+00
 0  <= FU_ADK1m  <= 1
 0  <= BU_ADK1m  <= 1
 -1000  <= DG_AKGDam  <= 1000
 -1.016370e+01  <= DGo_AKGDam  <= -6.266380e+00
 0  <= FU_AKGDam  <= 1
 0  <= BU_AKGDam  <= 1
 -1000  <= DG_AKGDbm  <= 1000
 -5.553001e+00  <= DGo_AKGDbm  <= -1.553001e+00
 0  <= FU_AKGDbm  <= 1
 0  <= BU_AKGDbm  <= 1
 -1000  <= DG_AKGMAL  <= 1000
 -1.425330e+00  <= DGo_AKGMAL  <= 2.574670e+00
 0  <= FU_AKGMAL  <= 1
 0  <= BU_AKGMAL  <= 1
 -1000  <= DG_AKGt2r  <= 1000
 -3.300790e+00  <= DGo_AKGt2r  <= 6.992103e-01
 0  <= FU_AKGt2r  <= 1
 0  <= BU_AKGt2r  <= 1
 -1000  <= DG_ALATA_L  <= 1000
 -6.109120e+00  <= DGo_ALATA_L  <= -2.109120e+00
 0  <= FU_ALATA_L  <= 1
 0  <= BU_ALATA_L  <= 1
 -1000  <= DG_ALATA_Lm  <= 1000
 -6.790665e+00  <= DGo_ALATA_Lm  <= -2.790665e+00
 0  <= FU_ALATA_Lm  <= 1
 0  <= BU_ALATA_Lm  <= 1
 -1000  <= DG_ALAt2r  <= 1000
 -4.755712e+00  <= DGo_ALAt2r  <= -7.557123e-01
 0  <= FU_ALAt2r  <= 1
 0  <= BU_ALAt2r  <= 1
 -1000  <= DG_ALAtmi  <= 1000
 -2.000000e+00  <= DGo_ALAtmi  <= 2.000000e+00
 0  <= FU_ALAtmi  <= 1
 0  <= BU_ALAtmi  <= 1
 -1000  <= DG_ALCD2if  <= 1000
 3.730024e+00  <= DGo_ALCD2if  <= 6.255139e+00
 0  <= FU_ALCD2if  <= 1
 0  <= BU_ALCD2if  <= 1
 -1000  <= DG_ALCD2ir  <= 1000
 -6.255139e+00  <= DGo_ALCD2ir  <= -3.730024e+00
 0  <= FU_ALCD2ir  <= 1
 0  <= BU_ALCD2ir  <= 1
 -1000  <= DG_ALCD2irm  <= 1000
 -5.573471e+00  <= DGo_ALCD2irm  <= -3.048356e+00
 0  <= FU_ALCD2irm  <= 1
 0  <= BU_ALCD2irm  <= 1
 -1000  <= DG_ALCD2x  <= 1000
 3.730024e+00  <= DGo_ALCD2x  <= 6.255139e+00
 0  <= FU_ALCD2x  <= 1
 0  <= BU_ALCD2x  <= 1
 -1000  <= DG_ALDD2xm  <= 1000
 -1.590265e+01  <= DGo_ALDD2xm  <= -1.318673e+01
 0  <= FU_ALDD2xm  <= 1
 0  <= BU_ALDD2xm  <= 1
 -1000  <= DG_ALDD2y  <= 1000
 -1.492753e+01  <= DGo_ALDD2y  <= -1.221161e+01
 0  <= FU_ALDD2y  <= 1
 0  <= BU_ALDD2y  <= 1
 -1000  <= DG_ALDD2ym  <= 1000
 -1.628968e+01  <= DGo_ALDD2ym  <= -1.357376e+01
 0  <= FU_ALDD2ym  <= 1
 0  <= BU_ALDD2ym  <= 1
 -1000  <= DG_ARGt2r  <= 1000
 -4.784758e+00  <= DGo_ARGt2r  <= -7.847582e-01
 0  <= FU_ARGt2r  <= 1
 0  <= BU_ARGt2r  <= 1
 -1000  <= DG_ASNt2r  <= 1000
 -4.562196e+00  <= DGo_ASNt2r  <= -5.621960e-01
 0  <= FU_ASNt2r  <= 1
 0  <= BU_ASNt2r  <= 1
 -1000  <= DG_ASPGLU2m  <= 1000
 -2.001205e+00  <= DGo_ASPGLU2m  <= 1.998795e+00
 0  <= FU_ASPGLU2m  <= 1
 0  <= BU_ASPGLU2m  <= 1
 -1000  <= DG_ASPTA  <= 1000
 -3.144379e+00  <= DGo_ASPTA  <= 8.556215e-01
 0  <= FU_ASPTA  <= 1
 0  <= BU_ASPTA  <= 1
 -1000  <= DG_ASPTAm  <= 1000
 -3.144941e+00  <= DGo_ASPTAm  <= 8.550587e-01
 0  <= FU_ASPTAm  <= 1
 0  <= BU_ASPTAm  <= 1
 -1000  <= DG_ASPt2m  <= 1000
 -3.361927e+00  <= DGo_ASPt2m  <= 6.380728e-01
 0  <= FU_ASPt2m  <= 1
 0  <= BU_ASPt2m  <= 1
 -1000  <= DG_ASPt2r  <= 1000
 -6.771766e+00  <= DGo_ASPt2r  <= -2.771766e+00
 0  <= FU_ASPt2r  <= 1
 0  <= BU_ASPt2r  <= 1
 -1000  <= DG_ATPH1  <= 1000
 -1.626307e+01  <= DGo_ATPH1  <= -1.339285e+01
 0  <= FU_ATPH1  <= 1
 0  <= BU_ATPH1  <= 1
 -1000  <= DG_ATPM  <= 1000
 -7.916015e+00  <= DGo_ATPM  <= -6.480906e+00
 0  <= FU_ATPM  <= 1
 0  <= BU_ATPM  <= 1
 -1000  <= DG_ATPS  <= 1000
 -7.916015e+00  <= DGo_ATPS  <= -6.480906e+00
 0  <= FU_ATPS  <= 1
 0  <= BU_ATPS  <= 1
 -1000  <= DG_ATPS3m  <= 1000
 -7.452873e+00  <= DGo_ATPS3m  <= -6.017764e+00
 0  <= FU_ATPS3m  <= 1
 0  <= BU_ATPS3m  <= 1
 -1000  <= DG_ATPtm-H  <= 1000
 -1.100720e+01  <= DGo_ATPtm-H  <= -7.007196e+00
 0  <= FU_ATPtm-H  <= 1
 0  <= BU_ATPtm-H  <= 1
 -1000  <= DG_CITt2r  <= 1000
 -5.719421e-01  <= DGo_CITt2r  <= 3.428058e+00
 0  <= FU_CITt2r  <= 1
 0  <= BU_CITt2r  <= 1
 -1000  <= DG_CITtam  <= 1000
 2.162711e+00  <= DGo_CITtam  <= 6.162711e+00
 0  <= FU_CITtam  <= 1
 0  <= BU_CITtam  <= 1
 -1000  <= DG_CITtcm  <= 1000
 -2.000312e+00  <= DGo_CITtcm  <= 1.999688e+00
 0  <= FU_CITtcm  <= 1
 0  <= BU_CITtcm  <= 1
 -1000  <= DG_CO2t  <= 1000
 -1.999994e+00  <= DGo_CO2t  <= 2.000006e+00
 0  <= FU_CO2t  <= 1
 0  <= BU_CO2t  <= 1
 -1000  <= DG_CO2tm  <= 1000
 -2.000000e+00  <= DGo_CO2tm  <= 2.000000e+00
 0  <= FU_CO2tm  <= 1
 0  <= BU_CO2tm  <= 1
 -1000  <= DG_CSm  <= 1000
 -1.008938e+01  <= DGo_CSm  <= -7.472163e+00
 0  <= FU_CSm  <= 1
 0  <= BU_CSm  <= 1
 -1000  <= DG_CYOOm  <= 1000
 -3.127616e+01  <= DGo_CYOOm  <= -2.904009e+01
 0  <= FU_CYOOm  <= 1
 0  <= BU_CYOOm  <= 1
 -1000  <= DG_CYOR_u6m  <= 1000
 6.908188e+00  <= DGo_CYOR_u6m  <= 1.218066e+01
 0  <= FU_CYOR_u6m  <= 1
 0  <= BU_CYOR_u6m  <= 1
 -1000  <= DG_CYSt2r  <= 1000
 -7.954480e+00  <= DGo_CYSt2r  <= -3.954480e+00
 0  <= FU_CYSt2r  <= 1
 0  <= BU_CYSt2r  <= 1
 -1000  <= DG_D-LACt2m  <= 1000
 -1.999873e+00  <= DGo_D-LACt2m  <= 2.000127e+00
 0  <= FU_D-LACt2m  <= 1
 0  <= BU_D-LACt2m  <= 1
 -1000  <= DG_D-LACtm  <= 1000
 -2.681541e+00  <= DGo_D-LACtm  <= 1.318459e+00
 0  <= FU_D-LACtm  <= 1
 0  <= BU_D-LACtm  <= 1
 -1000  <= DG_DB4PS  <= 1000
 -2.177804e+01  <= DGo_DB4PS  <= -2.058388e+01
 0  <= FU_DB4PS  <= 1
 0  <= BU_DB4PS  <= 1
 -1000  <= DG_ENO  <= 1000
 -2.312956e+00  <= DGo_ENO  <= -7.545290e-01
 0  <= FU_ENO  <= 1
 0  <= BU_ENO  <= 1
 -1000  <= DG_ERGSTt  <= 1000
 6.514723e+00  <= DGo_ERGSTt  <= 1.051472e+01
 0  <= FU_ERGSTt  <= 1
 0  <= BU_ERGSTt  <= 1
 -1000  <= DG_ERGSTter  <= 1000
 6.514717e+00  <= DGo_ERGSTter  <= 1.051472e+01
 0  <= FU_ERGSTter  <= 1
 0  <= BU_ERGSTter  <= 1
 -1000  <= DG_ETOHt  <= 1000
 -8.388963e-01  <= DGo_ETOHt  <= 3.161104e+00
 0  <= FU_ETOHt  <= 1
 0  <= BU_ETOHt  <= 1
 -1000  <= DG_ETOHtm  <= 1000
 -2.000000e+00  <= DGo_ETOHtm  <= 2.000000e+00
 0  <= FU_ETOHtm  <= 1
 0  <= BU_ETOHtm  <= 1
 -1000  <= DG_FBA  <= 1000
 2.857868e+00  <= DGo_FBA  <= 5.320548e+00
 0  <= FU_FBA  <= 1
 0  <= BU_FBA  <= 1
 -1000  <= DG_FBA3  <= 1000
 2.277285e+00  <= DGo_FBA3  <= 3.036539e+00
 0  <= FU_FBA3  <= 1
 0  <= BU_FBA3  <= 1
 -1000  <= DG_FBP  <= 1000
 -9.702464e-01  <= DGo_FBP  <= 5.010210e-01
 0  <= FU_FBP  <= 1
 0  <= BU_FBP  <= 1
 -1000  <= DG_FDH  <= 1000
 -4.739017e+00  <= DGo_FDH  <= -2.213902e+00
 0  <= FU_FDH  <= 1
 0  <= BU_FDH  <= 1
 -1000  <= DG_FDNG  <= 1000
 -2.149513e+01  <= DGo_FDNG  <= -1.603628e+01
 0  <= FU_FDNG  <= 1
 0  <= BU_FDNG  <= 1
 -1000  <= DG_FORt  <= 1000
 -5.158069e-01  <= DGo_FORt  <= 3.484193e+00
 0  <= FU_FORt  <= 1
 0  <= BU_FORt  <= 1
 -1000  <= DG_FRDcm  <= 1000
 -1.664094e+01  <= DGo_FRDcm  <= -1.280212e+01
 0  <= FU_FRDcm  <= 1
 0  <= BU_FRDcm  <= 1
 -1000  <= DG_FRDm  <= 1000
 -1.527329e+01  <= DGo_FRDm  <= -1.143446e+01
 0  <= FU_FRDm  <= 1
 0  <= BU_FRDm  <= 1
 -1000  <= DG_FUM  <= 1000
 -1.487998e+00  <= DGo_FUM  <= 2.574531e-01
 0  <= FU_FUM  <= 1
 0  <= BU_FUM  <= 1
 -1000  <= DG_FUMm  <= 1000
 -1.486446e+00  <= DGo_FUMm  <= 2.590055e-01
 0  <= FU_FUMm  <= 1
 0  <= BU_FUMm  <= 1
 -1000  <= DG_FUMt2r  <= 1000
 -3.530192e+00  <= DGo_FUMt2r  <= 4.698082e-01
 0  <= FU_FUMt2r  <= 1
 0  <= BU_FUMt2r  <= 1
 -1000  <= DG_G6PDH2  <= 1000
 -4.353352e+00  <= DGo_G6PDH2  <= -1.134932e+00
 0  <= FU_G6PDH2  <= 1
 0  <= BU_G6PDH2  <= 1
 -1000  <= DG_G6PDH2er  <= 1000
 -3.281277e+00  <= DGo_G6PDH2er  <= -6.285730e-02
 0  <= FU_G6PDH2er  <= 1
 0  <= BU_G6PDH2er  <= 1
 -1000  <= DG_G6PI  <= 1000
 -1.569691e+00  <= DGo_G6PI  <= 2.430309e+00
 0  <= FU_G6PI  <= 1
 0  <= BU_G6PI  <= 1
 -1000  <= DG_G6PI3  <= 1000
 -2.122078e+00  <= DGo_G6PI3  <= -9.043209e-01
 0  <= FU_G6PI3  <= 1
 0  <= BU_G6PI3  <= 1
 -1000  <= DG_G6Pter  <= 1000
 -3.452138e+00  <= DGo_G6Pter  <= 5.478617e-01
 0  <= FU_G6Pter  <= 1
 0  <= BU_G6Pter  <= 1
 -1000  <= DG_GAPD  <= 1000
 -2.709514e+00  <= DGo_GAPD  <= -3.194769e-01
 0  <= FU_GAPD  <= 1
 0  <= BU_GAPD  <= 1
 -1000  <= DG_GCC2cm  <= 1000
 -2.222119e+00  <= DGo_GCC2cm  <= 2.379946e+00
 0  <= FU_GCC2cm  <= 1
 0  <= BU_GCC2cm  <= 1
 -1000  <= DG_GLCt1  <= 1000
 3.222015e-01  <= DGo_GLCt1  <= 4.322202e+00
 0  <= FU_GLCt1  <= 1
 0  <= BU_GLCt1  <= 1
 -1000  <= DG_GLNt2r  <= 1000
 -4.175163e+00  <= DGo_GLNt2r  <= -1.751634e-01
 0  <= FU_GLNt2r  <= 1
 0  <= BU_GLNt2r  <= 1
 -1000  <= DG_GLUDC  <= 1000
 -2.490559e+00  <= DGo_GLUDC  <= -1.407159e+00
 0  <= FU_GLUDC  <= 1
 0  <= BU_GLUDC  <= 1
 -1000  <= DG_GLUDxi  <= 1000
 1.030951e+01  <= DGo_GLUDxi  <= 1.293383e+01
 0  <= FU_GLUDxi  <= 1
 0  <= BU_GLUDxi  <= 1
 -1000  <= DG_GLUDyi  <= 1000
 -1.254680e+01  <= DGo_GLUDyi  <= -9.922476e+00
 0  <= FU_GLUDyi  <= 1
 0  <= BU_GLUDyi  <= 1
 -1000  <= DG_GLUK  <= 1000
 -5.286255e+00  <= DGo_GLUK  <= -4.812796e+00
 0  <= FU_GLUK  <= 1
 0  <= BU_GLUK  <= 1
 -1000  <= DG_GLUt2r  <= 1000
 -6.771988e+00  <= DGo_GLUt2r  <= -2.771988e+00
 0  <= FU_GLUt2r  <= 1
 0  <= BU_GLUt2r  <= 1
 -1000  <= DG_GLUt7m  <= 1000
 1.469516e+00  <= DGo_GLUt7m  <= 5.469516e+00
 0  <= FU_GLUt7m  <= 1
 0  <= BU_GLUt7m  <= 1
 -1000  <= DG_GLYt2r  <= 1000
 -5.142745e+00  <= DGo_GLYt2r  <= -1.142745e+00
 0  <= FU_GLYt2r  <= 1
 0  <= BU_GLYt2r  <= 1
 -1000  <= DG_GND  <= 1000
 -7.099325e-01  <= DGo_GND  <= 1.744446e+00
 0  <= FU_GND  <= 1
 0  <= BU_GND  <= 1
 0  <= FU_H2Ot  <= 1
 0  <= BU_H2Ot  <= 1
 0  <= FU_H2Oter  <= 1
 0  <= BU_H2Oter  <= 1
 0  <= FU_H2Otm  <= 1
 0  <= BU_H2Otm  <= 1
 -1000  <= DG_HCITSm  <= 1000
 -1.008370e+01  <= DGo_HCITSm  <= -7.466486e+00
 0  <= FU_HCITSm  <= 1
 0  <= BU_HCITSm  <= 1
 -1000  <= DG_HCO3E  <= 1000
 -2.177560e+00  <= DGo_HCO3E  <= -4.455091e-01
 0  <= FU_HCO3E  <= 1
 0  <= BU_HCO3E  <= 1
 -1000  <= DG_HEX1  <= 1000
 -5.716564e+00  <= DGo_HEX1  <= -5.243104e+00
 0  <= FU_HEX1  <= 1
 0  <= BU_HEX1  <= 1
 -1000  <= DG_HISt2r  <= 1000
 -4.368680e+00  <= DGo_HISt2r  <= -3.686797e-01
 0  <= FU_HISt2r  <= 1
 0  <= BU_HISt2r  <= 1
 -1000  <= DG_ICDHxm  <= 1000
 2.282829e+00  <= DGo_ICDHxm  <= 4.764347e+00
 0  <= FU_ICDHxm  <= 1
 0  <= BU_ICDHxm  <= 1
 -1000  <= DG_ICDHy  <= 1000
 1.909747e+00  <= DGo_ICDHy  <= 4.391265e+00
 0  <= FU_ICDHy  <= 1
 0  <= BU_ICDHy  <= 1
 -1000  <= DG_ICDHym  <= 1000
 1.895797e+00  <= DGo_ICDHym  <= 4.377315e+00
 0  <= FU_ICDHym  <= 1
 0  <= BU_ICDHym  <= 1
 -1000  <= DG_ILEt2r  <= 1000
 -3.594615e+00  <= DGo_ILEt2r  <= 4.053855e-01
 0  <= FU_ILEt2r  <= 1
 0  <= BU_ILEt2r  <= 1
 -1000  <= DG_ILEtmi  <= 1000
 -2.000000e+00  <= DGo_ILEtmi  <= 2.000000e+00
 0  <= FU_ILEtmi  <= 1
 0  <= BU_ILEtmi  <= 1
 -1000  <= DG_LEUt2r  <= 1000
 -3.594615e+00  <= DGo_LEUt2r  <= 4.053855e-01
 0  <= FU_LEUt2r  <= 1
 0  <= BU_LEUt2r  <= 1
 -1000  <= DG_LYSt2r  <= 1000
 -4.784758e+00  <= DGo_LYSt2r  <= -7.847582e-01
 0  <= FU_LYSt2r  <= 1
 0  <= BU_LYSt2r  <= 1
 -1000  <= DG_MALt2r  <= 1000
 -2.726120e+00  <= DGo_MALt2r  <= 1.273880e+00
 0  <= FU_MALt2r  <= 1
 0  <= BU_MALt2r  <= 1
 -1000  <= DG_MALtm  <= 1000
 -2.139082e+00  <= DGo_MALtm  <= 1.860918e+00
 0  <= FU_MALtm  <= 1
 0  <= BU_MALtm  <= 1
 -1000  <= DG_MAN6PI  <= 1000
 -2.122078e+00  <= DGo_MAN6PI  <= -9.043209e-01
 0  <= FU_MAN6PI  <= 1
 0  <= BU_MAN6PI  <= 1
 -1000  <= DG_MANNANter  <= 1000
 -8.192522e+00  <= DGo_MANNANter  <= -4.192522e+00
 0  <= FU_MANNANter  <= 1
 0  <= BU_MANNANter  <= 1
 -1000  <= DG_MDH  <= 1000
 4.339342e+00  <= DGo_MDH  <= 6.529526e+00
 0  <= FU_MDH  <= 1
 0  <= BU_MDH  <= 1
 -1000  <= DG_MDHm  <= 1000
 3.656409e+00  <= DGo_MDHm  <= 5.846592e+00
 0  <= FU_MDHm  <= 1
 0  <= BU_MDHm  <= 1
 -1000  <= DG_ME1m  <= 1000
 1.039022e-01  <= DGo_ME1m  <= 2.537941e+00
 0  <= FU_ME1m  <= 1
 0  <= BU_ME1m  <= 1
 -1000  <= DG_ME2m  <= 1000
 -2.831304e-01  <= DGo_ME2m  <= 2.150909e+00
 0  <= FU_ME2m  <= 1
 0  <= BU_ME2m  <= 1
 -1000  <= DG_METt2r  <= 1000
 -3.981647e+00  <= DGo_METt2r  <= 1.835288e-02
 0  <= FU_METt2r  <= 1
 0  <= BU_METt2r  <= 1
 -1000  <= DG_NADH2-u6cm  <= 1000
 -1.812540e+01  <= DGo_NADH2-u6cm  <= -1.245309e+01
 0  <= FU_NADH2-u6cm  <= 1
 0  <= BU_NADH2-u6cm  <= 1
 -1000  <= DG_NADH2-u6m  <= 1000
 -1.880707e+01  <= DGo_NADH2-u6m  <= -1.313475e+01
 0  <= FU_NADH2-u6m  <= 1
 0  <= BU_NADH2-u6m  <= 1
 -1000  <= DG_NADK  <= 1000
 -4.498061e+00  <= DGo_NADK  <= -4.024601e+00
 0  <= FU_NADK  <= 1
 0  <= BU_NADK  <= 1
 -1000  <= DG_NADKm  <= 1000
 -4.521629e+00  <= DGo_NADKm  <= -4.048169e+00
 0  <= FU_NADKm  <= 1
 0  <= BU_NADKm  <= 1
 -1000  <= DG_NADPPPS  <= 1000
 -3.672763e+00  <= DGo_NADPPPS  <= -2.201496e+00
 0  <= FU_NADPPPS  <= 1
 0  <= BU_NADPPPS  <= 1
 -1000  <= DG_NADPPPSm  <= 1000
 -4.213360e+00  <= DGo_NADPPPSm  <= -2.742093e+00
 0  <= FU_NADPPPSm  <= 1
 0  <= BU_NADPPPSm  <= 1
 -1000  <= DG_NDP1  <= 1000
 -8.347054e+00  <= DGo_NDP1  <= -6.911945e+00
 0  <= FU_NDP1  <= 1
 0  <= BU_NDP1  <= 1
 -1000  <= DG_NH4t  <= 1000
 -2.803105e+00  <= DGo_NH4t  <= 1.196895e+00
 0  <= FU_NH4t  <= 1
 0  <= BU_NH4t  <= 1
 -1000  <= DG_NH4tm  <= 1000
 -6.150980e+00  <= DGo_NH4tm  <= -2.150980e+00
 0  <= FU_NH4tm  <= 1
 0  <= BU_NH4tm  <= 1
 -1000  <= DG_O2t  <= 1000
 -1.999994e+00  <= DGo_O2t  <= 2.000006e+00
 0  <= FU_O2t  <= 1
 0  <= BU_O2t  <= 1
 -1000  <= DG_O2ter  <= 1000
 -2  <= DGo_O2ter  <= 2
 0  <= FU_O2ter  <= 1
 0  <= BU_O2ter  <= 1
 -1000  <= DG_O2tm  <= 1000
 -2.000000e+00  <= DGo_O2tm  <= 2.000000e+00
 0  <= FU_O2tm  <= 1
 0  <= BU_O2tm  <= 1
 -1000  <= DG_OAAt  <= 1000
 -4.679558e+00  <= DGo_OAAt  <= -6.795577e-01
 0  <= FU_OAAt  <= 1
 0  <= BU_OAAt  <= 1
 -1000  <= DG_OAAt2m  <= 1000
 1.470035e+00  <= DGo_OAAt2m  <= 5.470035e+00
 0  <= FU_OAAt2m  <= 1
 0  <= BU_OAAt2m  <= 1
 -1000  <= DG_PC  <= 1000
 -2.514858e+00  <= DGo_PC  <= -1.035947e+00
 0  <= FU_PC  <= 1
 0  <= BU_PC  <= 1
 -1000  <= DG_PDHcm  <= 1000
 -2.222119e+00  <= DGo_PDHcm  <= 2.379946e+00
 0  <= FU_PDHcm  <= 1
 0  <= BU_PDHcm  <= 1
 -1000  <= DG_PDHm  <= 1000
 -1.212919e+01  <= DGo_PDHm  <= -8.933791e+00
 0  <= FU_PDHm  <= 1
 0  <= BU_PDHm  <= 1
 0  <= FU_PEtm_SC  <= 1
 0  <= BU_PEtm_SC  <= 1
 -1000  <= DG_PFK  <= 1000
 -7.200578e+00  <= DGo_PFK  <= -6.727118e+00
 0  <= FU_PFK  <= 1
 0  <= BU_PFK  <= 1
 -1000  <= DG_PFK_3  <= 1000
 -6.762324e+00  <= DGo_PFK_3  <= -6.288864e+00
 0  <= FU_PFK_3  <= 1
 0  <= BU_PFK_3  <= 1
 -1000  <= DG_PGI  <= 1000
 -1.691769e+00  <= DGo_PGI  <= -4.740121e-01
 0  <= FU_PGI  <= 1
 0  <= BU_PGI  <= 1
 -1000  <= DG_PGK  <= 1000
 4.897506e-01  <= DGo_PGK  <= 1.070529e+00
 0  <= FU_PGK  <= 1
 0  <= BU_PGK  <= 1
 -1000  <= DG_PGL  <= 1000
 -7.925485e+00  <= DGo_PGL  <= -5.083931e+00
 0  <= FU_PGL  <= 1
 0  <= BU_PGL  <= 1
 -1000  <= DG_PGM  <= 1000
 -2.097232e+00  <= DGo_PGM  <= 1.902768e+00
 0  <= FU_PGM  <= 1
 0  <= BU_PGM  <= 1
 -1000  <= DG_PGMT  <= 1000
 -2.101158e+00  <= DGo_PGMT  <= 1.898842e+00
 0  <= FU_PGMT  <= 1
 0  <= BU_PGMT  <= 1
 -1000  <= DG_PHEt2r  <= 1000
 -3.981647e+00  <= DGo_PHEt2r  <= 1.835288e-02
 0  <= FU_PHEt2r  <= 1
 0  <= BU_PHEt2r  <= 1
 -1000  <= DG_PIt2m  <= 1000
 1.610383e+00  <= DGo_PIt2m  <= 5.610383e+00
 0  <= FU_PIt2m  <= 1
 0  <= BU_PIt2m  <= 1
 -1000  <= DG_PIt2r  <= 1000
 -1.113141e+00  <= DGo_PIt2r  <= 2.886859e+00
 0  <= FU_PIt2r  <= 1
 0  <= BU_PIt2r  <= 1
 -1000  <= DG_PPA  <= 1000
 -6.657965e+00  <= DGo_PPA  <= -4.208475e+00
 0  <= FU_PPA  <= 1
 0  <= BU_PPA  <= 1
 -1000  <= DG_PPAm  <= 1000
 -7.086045e+00  <= DGo_PPAm  <= -4.636555e+00
 0  <= FU_PPAm  <= 1
 0  <= BU_PPAm  <= 1
 -1000  <= DG_PPCK  <= 1000
 -1.330418e+00  <= DGo_PPCK  <= 3.264029e-01
 0  <= FU_PPCK  <= 1
 0  <= BU_PPCK  <= 1
 -1000  <= DG_PROt2r  <= 1000
 -3.372052e+00  <= DGo_PROt2r  <= 6.279477e-01
 0  <= FU_PROt2r  <= 1
 0  <= BU_PROt2r  <= 1
 0  <= FU_PStm_SC  <= 1
 0  <= BU_PStm_SC  <= 1
 -1000  <= DG_PYK  <= 1000
 -4.263821e+00  <= DGo_PYK  <= -2.955211e+00
 0  <= FU_PYK  <= 1
 0  <= BU_PYK  <= 1
 -1000  <= DG_PYRDC  <= 1000
 -3.807788e+00  <= DGo_PYRDC  <= -2.271440e+00
 0  <= FU_PYRDC  <= 1
 0  <= BU_PYRDC  <= 1
 -1000  <= DG_PYRt  <= 1000
 -3.770698e+00  <= DGo_PYRt  <= 2.293015e-01
 0  <= FU_PYRt  <= 1
 0  <= BU_PYRt  <= 1
 -1000  <= DG_PYRt2m  <= 1000
 -2.681668e+00  <= DGo_PYRt2m  <= 1.318332e+00
 0  <= FU_PYRt2m  <= 1
 0  <= BU_PYRt2m  <= 1
 -1000  <= DG_RPE  <= 1000
 -2  <= DGo_RPE  <= 2
 0  <= FU_RPE  <= 1
 0  <= BU_RPE  <= 1
 -1000  <= DG_RPI  <= 1000
 -5.862703e-01  <= DGo_RPI  <= 2.016686e+00
 0  <= FU_RPI  <= 1
 0  <= BU_RPI  <= 1
 -1000  <= DG_SERt2r  <= 1000
 -4.755712e+00  <= DGo_SERt2r  <= -7.557123e-01
 0  <= FU_SERt2r  <= 1
 0  <= BU_SERt2r  <= 1
 -1000  <= DG_SO4ti  <= 1000
 -6.739360e-03  <= DGo_SO4ti  <= 3.993261e+00
 0  <= FU_SO4ti  <= 1
 0  <= BU_SO4ti  <= 1
 -1000  <= DG_SUCCt2r  <= 1000
 -2.281923e+00  <= DGo_SUCCt2r  <= 1.718077e+00
 0  <= FU_SUCCt2r  <= 1
 0  <= BU_SUCCt2r  <= 1
 -1000  <= DG_SUCCtm  <= 1000
 -2.136317e+00  <= DGo_SUCCtm  <= 1.863683e+00
 0  <= FU_SUCCtm  <= 1
 0  <= BU_SUCCtm  <= 1
 -1000  <= DG_SUCD1m  <= 1000
 1.143446e+01  <= DGo_SUCD1m  <= 1.527329e+01
 0  <= FU_SUCD1m  <= 1
 0  <= BU_SUCD1m  <= 1
 -1000  <= DG_SUCD2_u6m  <= 1000
 -4.992134e+00  <= DGo_SUCD2_u6m  <= 1.161493e-01
 0  <= FU_SUCD2_u6m  <= 1
 0  <= BU_SUCD2_u6m  <= 1
 -1000  <= DG_SUCD3-u6m  <= 1000
 -1.901193e+01  <= DGo_SUCD3-u6m  <= -1.257181e+01
 0  <= FU_SUCD3-u6m  <= 1
 0  <= BU_SUCD3-u6m  <= 1
 -1000  <= DG_SUCFUMtm  <= 1000
 -1.995683e+00  <= DGo_SUCFUMtm  <= 2.004317e+00
 0  <= FU_SUCFUMtm  <= 1
 0  <= BU_SUCFUMtm  <= 1
 -1000  <= DG_SUCOASm  <= 1000
 -5.300910e+00  <= DGo_SUCOASm  <= -2.874075e+00
 0  <= FU_SUCOASm  <= 1
 0  <= BU_SUCOASm  <= 1
 -1000  <= DG_TALA  <= 1000
 -2.279729e+00  <= DGo_TALA  <= 2.916437e-01
 0  <= FU_TALA  <= 1
 0  <= BU_TALA  <= 1
 -1000  <= DG_THRt2m  <= 1000
 -6.832648e+00  <= DGo_THRt2m  <= -2.832648e+00
 0  <= FU_THRt2m  <= 1
 0  <= BU_THRt2m  <= 1
 -1000  <= DG_THRt2r  <= 1000
 -4.368680e+00  <= DGo_THRt2r  <= -3.686797e-01
 0  <= FU_THRt2r  <= 1
 0  <= BU_THRt2r  <= 1
 -1000  <= DG_TKT1  <= 1000
 -5.447238e-01  <= DGo_TKT1  <= 2.136826e+00
 0  <= FU_TKT1  <= 1
 0  <= BU_TKT1  <= 1
 -1000  <= DG_TKT2  <= 1000
 -5.038042e+00  <= DGo_TKT2  <= -2.466670e+00
 0  <= FU_TKT2  <= 1
 0  <= BU_TKT2  <= 1
 -1000  <= DG_TPI  <= 1000
 1.159311e+00  <= DGo_TPI  <= 1.760689e+00
 0  <= FU_TPI  <= 1
 0  <= BU_TPI  <= 1
 -1000  <= DG_TRPt2r  <= 1000
 -3.788131e+00  <= DGo_TRPt2r  <= 2.118692e-01
 0  <= FU_TRPt2r  <= 1
 0  <= BU_TRPt2r  <= 1
 -1000  <= DG_TYRt2r  <= 1000
 -3.981647e+00  <= DGo_TYRt2r  <= 1.835288e-02
 0  <= FU_TYRt2r  <= 1
 0  <= BU_TYRt2r  <= 1
 -1000  <= DG_VALt2r  <= 1000
 -3.981647e+00  <= DGo_VALt2r  <= 1.835288e-02
 0  <= FU_VALt2r  <= 1
 0  <= BU_VALt2r  <= 1
 -1000  <= DG_ZYMSTt  <= 1000
 6.514723e+00  <= DGo_ZYMSTt  <= 1.051472e+01
 0  <= FU_ZYMSTt  <= 1
 0  <= BU_ZYMSTt  <= 1
 0  <= FU_LMPD_1_13BDglcn_c  <= 1
 0  <= BU_LMPD_1_13BDglcn_c  <= 1
 -1000  <= DG_LMPD_2_amp_c  <= 1000
 -2.207245e+01  <= DGo_LMPD_2_amp_c  <= -7.413045e+00
 0  <= FU_LMPD_2_amp_c  <= 1
 0  <= BU_LMPD_2_amp_c  <= 1
 -1000  <= DG_LMPD_3_amp_c  <= 1000
 -2.129838e+01  <= DGo_LMPD_3_amp_c  <= -6.638980e+00
 0  <= FU_LMPD_3_amp_c  <= 1
 0  <= BU_LMPD_3_amp_c  <= 1
 -1000  <= DG_LMPD_4_amp_c  <= 1000
 -2.250349e+01  <= DGo_LMPD_4_amp_c  <= -7.844084e+00
 0  <= FU_LMPD_4_amp_c  <= 1
 0  <= BU_LMPD_4_amp_c  <= 1
 -1000  <= DG_LMPD_5_amp_c  <= 1000
 -2.172942e+01  <= DGo_LMPD_5_amp_c  <= -7.070019e+00
 0  <= FU_LMPD_5_amp_c  <= 1
 0  <= BU_LMPD_5_amp_c  <= 1
 -1000  <= DG_LMPD_6_arg-L_c  <= 1000
 -1.081046e+01  <= DGo_LMPD_6_arg-L_c  <= -4.575409e+00
 0  <= FU_LMPD_6_arg-L_c  <= 1
 0  <= BU_LMPD_6_arg-L_c  <= 1
 -1000  <= DG_LMPD_7_asn-L_c  <= 1000
 -9.317493e+00  <= DGo_LMPD_7_asn-L_c  <= -6.943179e+00
 0  <= FU_LMPD_7_asn-L_c  <= 1
 0  <= BU_LMPD_7_asn-L_c  <= 1
 -1000  <= DG_LMPD_8_cmp_c  <= 1000
 -4.851561e+01  <= DGo_LMPD_8_cmp_c  <= -4.036130e+01
 0  <= FU_LMPD_8_cmp_c  <= 1
 0  <= BU_LMPD_8_cmp_c  <= 1
 -1000  <= DG_LMPD_9_cmp_c  <= 1000
 -4.677267e+01  <= DGo_LMPD_9_cmp_c  <= -3.996356e+01
 0  <= FU_LMPD_9_cmp_c  <= 1
 0  <= BU_LMPD_9_cmp_c  <= 1
 -1000  <= DG_LMPD_10_cmp_c  <= 1000
 -4.097224e+01  <= DGo_LMPD_10_cmp_c  <= -3.350774e+01
 0  <= FU_LMPD_10_cmp_c  <= 1
 0  <= BU_LMPD_10_cmp_c  <= 1
 -1000  <= DG_LMPD_11_cmp_c  <= 1000
 -3.915252e+01  <= DGo_LMPD_11_cmp_c  <= -3.318679e+01
 0  <= FU_LMPD_11_cmp_c  <= 1
 0  <= BU_LMPD_11_cmp_c  <= 1
 -1000  <= DG_LMPD_12_cys-L_c  <= 1000
 -4.867220e+01  <= DGo_LMPD_12_cys-L_c  <= -4.106520e+01
 0  <= FU_LMPD_12_cys-L_c  <= 1
 0  <= BU_LMPD_12_cys-L_c  <= 1
 -1000  <= DG_LMPD_13_damp_c  <= 1000
 -3.250603e+01  <= DGo_LMPD_13_damp_c  <= -1.902378e+01
 0  <= FU_LMPD_13_damp_c  <= 1
 0  <= BU_LMPD_13_damp_c  <= 1
 -1000  <= DG_LMPD_14_damp_c  <= 1000
 -3.293707e+01  <= DGo_LMPD_14_damp_c  <= -1.945482e+01
 0  <= FU_LMPD_14_damp_c  <= 1
 0  <= BU_LMPD_14_damp_c  <= 1
 -1000  <= DG_LMPD_15_damp_c  <= 1000
 -3.328009e+01  <= DGo_LMPD_15_damp_c  <= -1.979785e+01
 0  <= FU_LMPD_15_damp_c  <= 1
 0  <= BU_LMPD_15_damp_c  <= 1
 -1000  <= DG_LMPD_16_damp_c  <= 1000
 -3.371113e+01  <= DGo_LMPD_16_damp_c  <= -2.022889e+01
 0  <= FU_LMPD_16_damp_c  <= 1
 0  <= BU_LMPD_16_damp_c  <= 1
 -1000  <= DG_LMPD_17_dcmp_c  <= 1000
 -4.333215e+01  <= DGo_LMPD_17_dcmp_c  <= -3.743533e+01
 0  <= FU_LMPD_17_dcmp_c  <= 1
 0  <= BU_LMPD_17_dcmp_c  <= 1
 -1000  <= DG_LMPD_18_dcmp_c  <= 1000
 -4.515027e+01  <= DGo_LMPD_18_dcmp_c  <= -3.775789e+01
 0  <= FU_LMPD_18_dcmp_c  <= 1
 0  <= BU_LMPD_18_dcmp_c  <= 1
 -1000  <= DG_LMPD_19_dgmp_c  <= 1000
 -5.441948e+01  <= DGo_LMPD_19_dgmp_c  <= -3.837176e+01
 0  <= FU_LMPD_19_dgmp_c  <= 1
 0  <= BU_LMPD_19_dgmp_c  <= 1
 -1000  <= DG_LMPD_20_dgmp_c  <= 1000
 -4.674964e+01  <= DGo_LMPD_20_dgmp_c  <= -3.164467e+01
 0  <= FU_LMPD_20_dgmp_c  <= 1
 0  <= BU_LMPD_20_dgmp_c  <= 1
 -1000  <= DG_LMPD_21_dgmp_c  <= 1000
 -4.752370e+01  <= DGo_LMPD_21_dgmp_c  <= -3.241874e+01
 0  <= FU_LMPD_21_dgmp_c  <= 1
 0  <= BU_LMPD_21_dgmp_c  <= 1
 -1000  <= DG_LMPD_22_dgmp_c  <= 1000
 -5.519354e+01  <= DGo_LMPD_22_dgmp_c  <= -3.914582e+01
 0  <= FU_LMPD_22_dgmp_c  <= 1
 0  <= BU_LMPD_22_dgmp_c  <= 1
 -1000  <= DG_LMPD_23_dtmp_c  <= 1000
 -5.815703e+01  <= DGo_LMPD_23_dtmp_c  <= -5.094709e+01
 0  <= FU_LMPD_23_dtmp_c  <= 1
 0  <= BU_LMPD_23_dtmp_c  <= 1
 -1000  <= DG_LMPD_24_dtmp_c  <= 1000
 -5.642486e+01  <= DGo_LMPD_24_dtmp_c  <= -5.053858e+01
 0  <= FU_LMPD_24_dtmp_c  <= 1
 0  <= BU_LMPD_24_dtmp_c  <= 1
 -1000  <= DG_LMPD_25_ergst_c  <= 1000
 -1.413326e+03  <= DGo_LMPD_25_ergst_c  <= -1.340747e+03
 0  <= FU_LMPD_25_ergst_c  <= 1
 0  <= BU_LMPD_25_ergst_c  <= 1
 -1000  <= DG_LMPD_26_ergst_c  <= 1000
 -1.414487e+03  <= DGo_LMPD_26_ergst_c  <= -1.341908e+03
 0  <= FU_LMPD_26_ergst_c  <= 1
 0  <= BU_LMPD_26_ergst_c  <= 1
 -1000  <= DG_LMPD_27_gln-L_c  <= 1000
 1.518872e+00  <= DGo_LMPD_27_gln-L_c  <= 3.296815e+00
 0  <= FU_LMPD_27_gln-L_c  <= 1
 0  <= BU_LMPD_27_gln-L_c  <= 1
 -1000  <= DG_LMPD_28_gly_c  <= 1000
 3.779148e+00  <= DGo_LMPD_28_gly_c  <= 4.916305e+00
 0  <= FU_LMPD_28_gly_c  <= 1
 0  <= BU_LMPD_28_gly_c  <= 1
 0  <= FU_LMPD_29_glycogen_c  <= 1
 0  <= BU_LMPD_29_glycogen_c  <= 1
 -1000  <= DG_LMPD_30_gmp_c  <= 1000
 -3.607855e+01  <= DGo_LMPD_30_gmp_c  <= -1.949067e+01
 0  <= FU_LMPD_30_gmp_c  <= 1
 0  <= BU_LMPD_30_gmp_c  <= 1
 -1000  <= DG_LMPD_31_gmp_c  <= 1000
 -3.685261e+01  <= DGo_LMPD_31_gmp_c  <= -2.026473e+01
 0  <= FU_LMPD_31_gmp_c  <= 1
 0  <= BU_LMPD_31_gmp_c  <= 1
 -1000  <= DG_LMPD_32_his-L_c  <= 1000
 -5.288775e+01  <= DGo_LMPD_32_his-L_c  <= -4.242642e+01
 0  <= FU_LMPD_32_his-L_c  <= 1
 0  <= BU_LMPD_32_his-L_c  <= 1
 -1000  <= DG_LMPD_33_his-L_c  <= 1000
 -5.206968e+01  <= DGo_LMPD_33_his-L_c  <= -4.160835e+01
 0  <= FU_LMPD_33_his-L_c  <= 1
 0  <= BU_LMPD_33_his-L_c  <= 1
 -1000  <= DG_LMPD_34_his-L_c  <= 1000
 -5.245671e+01  <= DGo_LMPD_34_his-L_c  <= -4.199538e+01
 0  <= FU_LMPD_34_his-L_c  <= 1
 0  <= BU_LMPD_34_his-L_c  <= 1
 -1000  <= DG_LMPD_35_his-L_c  <= 1000
 -5.250072e+01  <= DGo_LMPD_35_his-L_c  <= -4.203939e+01
 0  <= FU_LMPD_35_his-L_c  <= 1
 0  <= BU_LMPD_35_his-L_c  <= 1
 -1000  <= DG_LMPD_36_ile-L_c  <= 1000
 -2.353731e+01  <= DGo_LMPD_36_ile-L_c  <= -1.675333e+01
 0  <= FU_LMPD_36_ile-L_c  <= 1
 0  <= BU_LMPD_36_ile-L_c  <= 1
 -1000  <= DG_LMPD_37_ile-L_c  <= 1000
 -2.392434e+01  <= DGo_LMPD_37_ile-L_c  <= -1.714036e+01
 0  <= FU_LMPD_37_ile-L_c  <= 1
 0  <= BU_LMPD_37_ile-L_c  <= 1
 -1000  <= DG_LMPD_38_leu-L_c  <= 1000
 -1.548231e+01  <= DGo_LMPD_38_leu-L_c  <= -1.239609e+01
 0  <= FU_LMPD_38_leu-L_c  <= 1
 0  <= BU_LMPD_38_leu-L_c  <= 1
 -1000  <= DG_LMPD_39_leu-L_c  <= 1000
 -2.618713e+01  <= DGo_LMPD_39_leu-L_c  <= -2.310091e+01
 0  <= FU_LMPD_39_leu-L_c  <= 1
 0  <= BU_LMPD_39_leu-L_c  <= 1
 -1000  <= DG_LMPD_40_lys-L_c  <= 1000
 1.143178e+01  <= DGo_LMPD_40_lys-L_c  <= 1.343994e+01
 0  <= FU_LMPD_40_lys-L_c  <= 1
 0  <= BU_LMPD_40_lys-L_c  <= 1
 -1000  <= DG_LMPD_41_lys-L_c  <= 1000
 1.104475e+01  <= DGo_LMPD_41_lys-L_c  <= 1.305291e+01
 0  <= FU_LMPD_41_lys-L_c  <= 1
 0  <= BU_LMPD_41_lys-L_c  <= 1
 0  <= FU_LMPD_42_mannan_c  <= 1
 0  <= BU_LMPD_42_mannan_c  <= 1
 0  <= FU_LMPD_43_mannan_c  <= 1
 0  <= BU_LMPD_43_mannan_c  <= 1
 -1000  <= DG_LMPD_44_met-L_c  <= 1000
 -5.059117e+01  <= DGo_LMPD_44_met-L_c  <= -3.726038e+01
 0  <= FU_LMPD_44_met-L_c  <= 1
 0  <= BU_LMPD_44_met-L_c  <= 1
 -1000  <= DG_LMPD_45_met-L_c  <= 1000
 -5.020414e+01  <= DGo_LMPD_45_met-L_c  <= -3.687334e+01
 0  <= FU_LMPD_45_met-L_c  <= 1
 0  <= BU_LMPD_45_met-L_c  <= 1
 0  <= FU_LMPD_47_pa_SC_c  <= 1
 0  <= BU_LMPD_47_pa_SC_c  <= 1
 0  <= FU_LMPD_48_pc_SC_c  <= 1
 0  <= BU_LMPD_48_pc_SC_c  <= 1
 0  <= FU_LMPD_49_pc_SC_c  <= 1
 0  <= BU_LMPD_49_pc_SC_c  <= 1
 0  <= FU_LMPD_51_pe_SC_c  <= 1
 0  <= BU_LMPD_51_pe_SC_c  <= 1
 -1000  <= DG_LMPD_52_phe-L_c  <= 1000
 -9.732246e+01  <= DGo_LMPD_52_phe-L_c  <= -9.168633e+01
 0  <= FU_LMPD_52_phe-L_c  <= 1
 0  <= BU_LMPD_52_phe-L_c  <= 1
 -1000  <= DG_LMPD_53_pro-L_c  <= 1000
 -6.272078e+00  <= DGo_LMPD_53_pro-L_c  <= -1.736562e+00
 0  <= FU_LMPD_53_pro-L_c  <= 1
 0  <= BU_LMPD_53_pro-L_c  <= 1
 -1000  <= DG_LMPD_54_pro-L_c  <= 1000
 -5.885045e+00  <= DGo_LMPD_54_pro-L_c  <= -1.349529e+00
 0  <= FU_LMPD_54_pro-L_c  <= 1
 0  <= BU_LMPD_54_pro-L_c  <= 1
 0  <= FU_LMPD_56_ps_SC_c  <= 1
 0  <= BU_LMPD_56_ps_SC_c  <= 1
 0  <= FU_LMPD_57_ptd1ino_SC_c  <= 1
 0  <= BU_LMPD_57_ptd1ino_SC_c  <= 1
 0  <= FU_LMPD_58_ptd1ino_SC_c  <= 1
 0  <= BU_LMPD_58_ptd1ino_SC_c  <= 1
 -1000  <= DG_LMPD_59_ribflv_c  <= 1000
 -1.321966e+02  <= DGo_LMPD_59_ribflv_c  <= -1.170125e+02
 0  <= FU_LMPD_59_ribflv_c  <= 1
 0  <= BU_LMPD_59_ribflv_c  <= 1
 -1000  <= DG_LMPD_60_ribflv_c  <= 1000
 -1.317655e+02  <= DGo_LMPD_60_ribflv_c  <= -1.165814e+02
 0  <= FU_LMPD_60_ribflv_c  <= 1
 0  <= BU_LMPD_60_ribflv_c  <= 1
 -1000  <= DG_LMPD_61_ser-L_c  <= 1000
 5.967069e+00  <= DGo_LMPD_61_ser-L_c  <= 8.581361e+00
 0  <= FU_LMPD_61_ser-L_c  <= 1
 0  <= BU_LMPD_61_ser-L_c  <= 1
 -1000  <= DG_LMPD_62_thr-L_c  <= 1000
 -8.308196e+00  <= DGo_LMPD_62_thr-L_c  <= -3.509471e+00
 0  <= FU_LMPD_62_thr-L_c  <= 1
 0  <= BU_LMPD_62_thr-L_c  <= 1
 -1000  <= DG_LMPD_63_thr-L_c  <= 1000
 -7.921163e+00  <= DGo_LMPD_63_thr-L_c  <= -3.122439e+00
 0  <= FU_LMPD_63_thr-L_c  <= 1
 0  <= BU_LMPD_63_thr-L_c  <= 1
 -1000  <= DG_LMPD_64_tre_c  <= 1000
 -3.206758e+00  <= DGo_LMPD_64_tre_c  <= -1.174545e+00
 0  <= FU_LMPD_64_tre_c  <= 1
 0  <= BU_LMPD_64_tre_c  <= 1
 0  <= FU_LMPD_65_triglyc_SC_c  <= 1
 0  <= BU_LMPD_65_triglyc_SC_c  <= 1
 0  <= FU_LMPD_66_triglyc_SC_c  <= 1
 0  <= BU_LMPD_66_triglyc_SC_c  <= 1
 -1000  <= DG_LMPD_67_trp-L_c  <= 1000
 -1.367301e+02  <= DGo_LMPD_67_trp-L_c  <= -1.285715e+02
 0  <= FU_LMPD_67_trp-L_c  <= 1
 0  <= BU_LMPD_67_trp-L_c  <= 1
 -1000  <= DG_LMPD_68_tyr-L_c  <= 1000
 -8.492594e+01  <= DGo_LMPD_68_tyr-L_c  <= -7.976368e+01
 0  <= FU_LMPD_68_tyr-L_c  <= 1
 0  <= BU_LMPD_68_tyr-L_c  <= 1
 -1000  <= DG_LMPD_69_tyr-L_c  <= 1000
 -8.531298e+01  <= DGo_LMPD_69_tyr-L_c  <= -8.015071e+01
 0  <= FU_LMPD_69_tyr-L_c  <= 1
 0  <= BU_LMPD_69_tyr-L_c  <= 1
 -1000  <= DG_LMPD_70_ump_c  <= 1000
 -3.652596e+01  <= DGo_LMPD_70_ump_c  <= -3.167572e+01
 0  <= FU_LMPD_70_ump_c  <= 1
 0  <= BU_LMPD_70_ump_c  <= 1
 -1000  <= DG_LMPD_71_ump_c  <= 1000
 -3.854613e+01  <= DGo_LMPD_71_ump_c  <= -3.179623e+01
 0  <= FU_LMPD_71_ump_c  <= 1
 0  <= BU_LMPD_71_ump_c  <= 1
 -1000  <= DG_LMPD_72_val-L_c  <= 1000
 -2.006356e+01  <= DGo_LMPD_72_val-L_c  <= -1.745920e+01
 0  <= FU_LMPD_72_val-L_c  <= 1
 0  <= BU_LMPD_72_val-L_c  <= 1
 -1000  <= DG_LMPD_73_zymst_c  <= 1000
 -1.264117e+03  <= DGo_LMPD_73_zymst_c  <= -1.197270e+03
 0  <= FU_LMPD_73_zymst_c  <= 1
 0  <= BU_LMPD_73_zymst_c  <= 1
 -1000  <= DG_LMPD_74_zymst_c  <= 1000
 -1.262956e+03  <= DGo_LMPD_74_zymst_c  <= -1.196109e+03
 0  <= FU_LMPD_74_zymst_c  <= 1
 0  <= BU_LMPD_74_zymst_c  <= 1
Binaries
 FU_biomass_SC5_notrace	 BU_biomass_SC5_notrace	 FU_EX_4abut_e	 BU_EX_4abut_e	 FU_EX_ac_e	 BU_EX_ac_e	
 FU_EX_acald_e	 BU_EX_acald_e	 FU_EX_akg_e	 BU_EX_akg_e	 FU_EX_ala-L_e	 BU_EX_ala-L_e	
 FU_EX_arg-L_e	 BU_EX_arg-L_e	 FU_EX_asn-L_e	 BU_EX_asn-L_e	 FU_EX_asp-L_e	 BU_EX_asp-L_e	
 FU_EX_cit_e	 BU_EX_cit_e	 FU_EX_co2_e	 BU_EX_co2_e	 FU_EX_cys-L_e	 BU_EX_cys-L_e	
 FU_EX_ergst_e	 BU_EX_ergst_e	 FU_EX_etoh_e	 BU_EX_etoh_e	 FU_EX_for_e	 BU_EX_for_e	
 FU_EX_fum_e	 BU_EX_fum_e	 FU_EX_glc_e	 BU_EX_glc_e	 FU_EX_gln-L_e	 BU_EX_gln-L_e	
 FU_EX_glu-L_e	 BU_EX_glu-L_e	 FU_EX_gly_e	 BU_EX_gly_e	 FU_EX_h_e	 BU_EX_h_e	
 FU_EX_h2o_e	 BU_EX_h2o_e	 FU_EX_his-L_e	 BU_EX_his-L_e	 FU_EX_ile-L_e	 BU_EX_ile-L_e	
 FU_EX_leu-L_e	 BU_EX_leu-L_e	 FU_EX_lys-L_e	 BU_EX_lys-L_e	 FU_EX_mal-L_e	 BU_EX_mal-L_e	
 FU_EX_met-L_e	 BU_EX_met-L_e	 FU_EX_nh4_e	 BU_EX_nh4_e	 FU_EX_o2_e	 BU_EX_o2_e	
 FU_EX_oaa_e	 BU_EX_oaa_e	 FU_EX_phe-L_e	 BU_EX_phe-L_e	 FU_EX_pi_e	 BU_EX_pi_e	
 FU_EX_pro-L_e	 BU_EX_pro-L_e	 FU_EX_pyr_e	 BU_EX_pyr_e	 FU_EX_ser-L_e	 BU_EX_ser-L_e	
 FU_EX_so4_e	 BU_EX_so4_e	 FU_EX_succ_e	 BU_EX_succ_e	 FU_EX_thr-L_e	 BU_EX_thr-L_e	
 FU_EX_trp-L_e	 BU_EX_trp-L_e	 FU_EX_tyr-L_e	 BU_EX_tyr-L_e	 FU_EX_val-L_e	 BU_EX_val-L_e	
 FU_EX_zymst_e	 BU_EX_zymst_e	 FU_6PGLter	 BU_6PGLter	 FU_ABUTt2r	 BU_ABUTt2r	
 FU_ACALDt	 BU_ACALDt	 FU_ACALDtm	 BU_ACALDtm	 FU_ACOAHim	 BU_ACOAHim	
 FU_ACONT	 BU_ACONT	 FU_ACONTm	 BU_ACONTm	 FU_ACS	 BU_ACS	
 FU_ACSm	 BU_ACSm	 FU_ACt2r	 BU_ACt2r	 FU_ACtm	 BU_ACtm	
 FU_ACtr	 BU_ACtr	 FU_ADK1	 BU_ADK1	 FU_ADK1m	 BU_ADK1m	
 FU_AKGDam	 BU_AKGDam	 FU_AKGDbm	 BU_AKGDbm	 FU_AKGMAL	 BU_AKGMAL	
 FU_AKGt2r	 BU_AKGt2r	 FU_ALATA_L	 BU_ALATA_L	 FU_ALATA_Lm	 BU_ALATA_Lm	
 FU_ALAt2r	 BU_ALAt2r	 FU_ALAtmi	 BU_ALAtmi	 FU_ALCD2if	 BU_ALCD2if	
 FU_ALCD2ir	 BU_ALCD2ir	 FU_ALCD2irm	 BU_ALCD2irm	 FU_ALCD2x	 BU_ALCD2x	
 FU_ALDD2xm	 BU_ALDD2xm	 FU_ALDD2y	 BU_ALDD2y	 FU_ALDD2ym	 BU_ALDD2ym	
 FU_ARGt2r	 BU_ARGt2r	 FU_ASNt2r	 BU_ASNt2r	 FU_ASPGLU2m	 BU_ASPGLU2m	
 FU_ASPTA	 BU_ASPTA	 FU_ASPTAm	 BU_ASPTAm	 FU_ASPt2m	 BU_ASPt2m	
 FU_ASPt2r	 BU_ASPt2r	 FU_ATPH1	 BU_ATPH1	 FU_ATPM	 BU_ATPM	
 FU_ATPS	 BU_ATPS	 FU_ATPS3m	 BU_ATPS3m	 FU_ATPtm-H	 BU_ATPtm-H	
 FU_CITt2r	 BU_CITt2r	 FU_CITtam	 BU_CITtam	 FU_CITtcm	 BU_CITtcm	
 FU_CO2t	 BU_CO2t	 FU_CO2tm	 BU_CO2tm	 FU_CSm	 BU_CSm	
 FU_CYOOm	 BU_CYOOm	 FU_CYOR_u6m	 BU_CYOR_u6m	 FU_CYSt2r	 BU_CYSt2r	
 FU_D-LACt2m	 BU_D-LACt2m	 FU_D-LACtm	 BU_D-LACtm	 FU_DB4PS	 BU_DB4PS	
 FU_ENO	 BU_ENO	 FU_ERGSTt	 BU_ERGSTt	 FU_ERGSTter	 BU_ERGSTter	
 FU_ETOHt	 BU_ETOHt	 FU_ETOHtm	 BU_ETOHtm	 FU_FBA	 BU_FBA	
 FU_FBA3	 BU_FBA3	 FU_FBP	 BU_FBP	 FU_FDH	 BU_FDH	
 FU_FDNG	 BU_FDNG	 FU_FORt	 BU_FORt	 FU_FRDcm	 BU_FRDcm	
 FU_FRDm	 BU_FRDm	 FU_FUM	 BU_FUM	 FU_FUMm	 BU_FUMm	
 FU_FUMt2r	 BU_FUMt2r	 FU_G6PDH2	 BU_G6PDH2	 FU_G6PDH2er	 BU_G6PDH2er	
 FU_G6PI	 BU_G6PI	 FU_G6PI3	 BU_G6PI3	 FU_G6Pter	 BU_G6Pter	
 FU_GAPD	 BU_GAPD	 FU_GCC2cm	 BU_GCC2cm	 FU_GLCt1	 BU_GLCt1	
 FU_GLNt2r	 BU_GLNt2r	 FU_GLUDC	 BU_GLUDC	 FU_GLUDxi	 BU_GLUDxi	
 FU_GLUDyi	 BU_GLUDyi	 FU_GLUK	 BU_GLUK	 FU_GLUt2r	 BU_GLUt2r	
 FU_GLUt7m	 BU_GLUt7m	 FU_GLYt2r	 BU_GLYt2r	 FU_GND	 BU_GND	
 FU_H2Ot	 BU_H2Ot	 FU_H2Oter	 BU_H2Oter	 FU_H2Otm	 BU_H2Otm	
 FU_HCITSm	 BU_HCITSm	 FU_HCO3E	 BU_HCO3E	 FU_HEX1	 BU_HEX1	
 FU_HISt2r	 BU_HISt2r	 FU_ICDHxm	 BU_ICDHxm	 FU_ICDHy	 BU_ICDHy	
 FU_ICDHym	 BU_ICDHym	 FU_ILEt2r	 BU_ILEt2r	 FU_ILEtmi	 BU_ILEtmi	
 FU_LEUt2r	 BU_LEUt2r	 FU_LYSt2r	 BU_LYSt2r	 FU_MALt2r	 BU_MALt2r	
 FU_MALtm	 BU_MALtm	 FU_MAN6PI	 BU_MAN6PI	 FU_MANNANter	 BU_MANNANter	
 FU_MDH	 BU_MDH	 FU_MDHm	 BU_MDHm	 FU_ME1m	 BU_ME1m	
 FU_ME2m	 BU_ME2m	 FU_METt2r	 BU_METt2r	 FU_NADH2-u6cm	 BU_NADH2-u6cm	
 FU_NADH2-u6m	 BU_NADH2-u6m	 FU_NADK	 BU_NADK	 FU_NADKm	 BU_NADKm	
 FU_NADPPPS	 BU_NADPPPS	 FU_NADPPPSm	 BU_NADPPPSm	 FU_NDP1	 BU_NDP1	
 FU_NH4t	 BU_NH4t	 FU_NH4tm	 BU_NH4tm	 FU_O2t	 BU_O2t	
 FU_O2ter	 BU_O2ter	 FU_O2tm	 BU_O2tm	 FU_OAAt	 BU_OAAt	
 FU_OAAt2m	 BU_OAAt2m	 FU_PC	 BU_PC	 FU_PDHcm	 BU_PDHcm	
 FU_PDHm	 BU_PDHm	 FU_PEtm_SC	 BU_PEtm_SC	 FU_PFK	 BU_PFK	
 FU_PFK_3	 BU_PFK_3	 FU_PGI	 BU_PGI	 FU_PGK	 BU_PGK	
 FU_PGL	 BU_PGL	 FU_PGM	 BU_PGM	 FU_PGMT	 BU_PGMT	
 FU_PHEt2r	 BU_PHEt2r	 FU_PIt2m	 BU_PIt2m	 FU_PIt2r	 BU_PIt2r	
 FU_PPA	 BU_PPA	 FU_PPAm	 BU_PPAm	 FU_PPCK	 BU_PPCK	
 FU_PROt2r	 BU_PROt2r	 FU_PStm_SC	 BU_PStm_SC	 FU_PYK	 BU_PYK	
 FU_PYRDC	 BU_PYRDC	 FU_PYRt	 BU_PYRt	 FU_PYRt2m	 BU_PYRt2m	
 FU_RPE	 BU_RPE	 FU_RPI	 BU_RPI	 FU_SERt2r	 BU_SERt2r	
 FU_SO4ti	 BU_SO4ti	 FU_SUCCt2r	 BU_SUCCt2r	 FU_SUCCtm	 BU_SUCCtm	
 FU_SUCD1m	 BU_SUCD1m	 FU_SUCD2_u6m	 BU_SUCD2_u6m	 FU_SUCD3-u6m	 BU_SUCD3-u6m	
 FU_SUCFUMtm	 BU_SUCFUMtm	 FU_SUCOASm	 BU_SUCOASm	 FU_TALA	 BU_TALA	
 FU_THRt2m	 BU_THRt2m	 FU_THRt2r	 BU_THRt2r	 FU_TKT1	 BU_TKT1	
 FU_TKT2	 BU_TKT2	 FU_TPI	 BU_TPI	 FU_TRPt2r	 BU_TRPt2r	
 FU_TYRt2r	 BU_TYRt2r	 FU_VALt2r	 BU_VALt2r	 FU_ZYMSTt	 BU_ZYMSTt	
 FU_LMPD_1_13BDglcn_c	 BU_LMPD_1_13BDglcn_c	 FU_LMPD_2_amp_c	 BU_LMPD_2_amp_c	 FU_LMPD_3_amp_c	 BU_LMPD_3_amp_c	
 FU_LMPD_4_amp_c	 BU_LMPD_4_amp_c	 FU_LMPD_5_amp_c	 BU_LMPD_5_amp_c	 FU_LMPD_6_arg-L_c	 BU_LMPD_6_arg-L_c	
 FU_LMPD_7_asn-L_c	 BU_LMPD_7_asn-L_c	 FU_LMPD_8_cmp_c	 BU_LMPD_8_cmp_c	 FU_LMPD_9_cmp_c	 BU_LMPD_9_cmp_c	
 FU_LMPD_10_cmp_c	 BU_LMPD_10_cmp_c	 FU_LMPD_11_cmp_c	 BU_LMPD_11_cmp_c	 FU_LMPD_12_cys-L_c	 BU_LMPD_12_cys-L_c	
 FU_LMPD_13_damp_c	 BU_LMPD_13_damp_c	 FU_LMPD_14_damp_c	 BU_LMPD_14_damp_c	 FU_LMPD_15_damp_c	 BU_LMPD_15_damp_c	
 FU_LMPD_16_damp_c	 BU_LMPD_16_damp_c	 FU_LMPD_17_dcmp_c	 BU_LMPD_17_dcmp_c	 FU_LMPD_18_dcmp_c	 BU_LMPD_18_dcmp_c	
 FU_LMPD_19_dgmp_c	 BU_LMPD_19_dgmp_c	 FU_LMPD_20_dgmp_c	 BU_LMPD_20_dgmp_c	 FU_LMPD_21_dgmp_c	 BU_LMPD_21_dgmp_c	
 FU_LMPD_22_dgmp_c	 BU_LMPD_22_dgmp_c	 FU_LMPD_23_dtmp_c	 BU_LMPD_23_dtmp_c	 FU_LMPD_24_dtmp_c	 BU_LMPD_24_dtmp_c	
 FU_LMPD_25_ergst_c	 BU_LMPD_25_ergst_c	 FU_LMPD_26_ergst_c	 BU_LMPD_26_ergst_c	 FU_LMPD_27_gln-L_c	 BU_LMPD_27_gln-L_c	
 FU_LMPD_28_gly_c	 BU_LMPD_28_gly_c	 FU_LMPD_29_glycogen_c	 BU_LMPD_29_glycogen_c	 FU_LMPD_30_gmp_c	 BU_LMPD_30_gmp_c	
 FU_LMPD_31_gmp_c	 BU_LMPD_31_gmp_c	 FU_LMPD_32_his-L_c	 BU_LMPD_32_his-L_c	 FU_LMPD_33_his-L_c	 BU_LMPD_33_his-L_c	
 FU_LMPD_34_his-L_c	 BU_LMPD_34_his-L_c	 FU_LMPD_35_his-L_c	 BU_LMPD_35_his-L_c	 FU_LMPD_36_ile-L_c	 BU_LMPD_36_ile-L_c	
 FU_LMPD_37_ile-L_c	 BU_LMPD_37_ile-L_c	 FU_LMPD_38_leu-L_c	 BU_LMPD_38_leu-L_c	 FU_LMPD_39_leu-L_c	 BU_LMPD_39_leu-L_c	
 FU_LMPD_40_lys-L_c	 BU_LMPD_40_lys-L_c	 FU_LMPD_41_lys-L_c	 BU_LMPD_41_lys-L_c	 FU_LMPD_42_mannan_c	 BU_LMPD_42_mannan_c	
 FU_LMPD_43_mannan_c	 BU_LMPD_43_mannan_c	 FU_LMPD_44_met-L_c	 BU_LMPD_44_met-L_c	 FU_LMPD_45_met-L_c	 BU_LMPD_45_met-L_c	
 FU_LMPD_47_pa_SC_c	 BU_LMPD_47_pa_SC_c	 FU_LMPD_48_pc_SC_c	 BU_LMPD_48_pc_SC_c	 FU_LMPD_49_pc_SC_c	 BU_LMPD_49_pc_SC_c	
 FU_LMPD_51_pe_SC_c	 BU_LMPD_51_pe_SC_c	 FU_LMPD_52_phe-L_c	 BU_LMPD_52_phe-L_c	 FU_LMPD_53_pro-L_c	 BU_LMPD_53_pro-L_c	
 FU_LMPD_54_pro-L_c	 BU_LMPD_54_pro-L_c	 FU_LMPD_56_ps_SC_c	 BU_LMPD_56_ps_SC_c	 FU_LMPD_57_ptd1ino_SC_c	 BU_LMPD_57_ptd1ino_SC_c	
 FU_LMPD_58_ptd1ino_SC_c	 BU_LMPD_58_ptd1ino_SC_c	 FU_LMPD_59_ribflv_c	 BU_LMPD_59_ribflv_c	 FU_LMPD_60_ribflv_c	 BU_LMPD_60_ribflv_c	
 FU_LMPD_61_ser-L_c	 BU_LMPD_61_ser-L_c	 FU_LMPD_62_thr-L_c	 BU_LMPD_62_thr-L_c	 FU_LMPD_63_thr-L_c	 BU_LMPD_63_thr-L_c	
 FU_LMPD_64_tre_c	 BU_LMPD_64_tre_c	 FU_LMPD_65_triglyc_SC_c	 BU_LMPD_65_triglyc_SC_c	 FU_LMPD_66_triglyc_SC_c	 BU_LMPD_66_triglyc_SC_c	
 FU_LMPD_67_trp-L_c	 BU_LMPD_67_trp-L_c	 FU_LMPD_68_tyr-L_c	 BU_LMPD_68_tyr-L_c	 FU_LMPD_69_tyr-L_c	 BU_LMPD_69_tyr-L_c	
 FU_LMPD_70_ump_c	 BU_LMPD_70_ump_c	 FU_LMPD_71_ump_c	 BU_LMPD_71_ump_c	 FU_LMPD_72_val-L_c	 BU_LMPD_72_val-L_c	
 FU_LMPD_73_zymst_c	 BU_LMPD_73_zymst_c	 FU_LMPD_74_zymst_c	 BU_LMPD_74_zymst_c	
End