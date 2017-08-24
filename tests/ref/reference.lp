\Problem name: matlab_new_formulation_LP

Maximize
 obj: F_Ec_biomass_iJO1366_WT_53p95M
Subject To
M_10fthf_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_1_10fthf_c - R_LMPD_1_10fthf_c + F_LMPD_2_10fthf_c - R_LMPD_2_10fthf_c + F_LMPD_3_10fthf_c - R_LMPD_3_10fthf_c + F_LMPD_4_10fthf_c - R_LMPD_4_10fthf_c = 0
M_13dpg_c: F_GAPD - R_GAPD + F_PGK - R_PGK = 0
M_2ddg6p_c: - F_EDA + R_EDA + F_EDD - R_EDD = 0
M_2dmmq8_c: F_FRD3 - R_FRD3 - F_G3PD7 + R_G3PD7 - F_GLYCTO4 + R_GLYCTO4 - F_NADH18pp + R_NADH18pp - F_NADH9 + R_NADH9 - F_NADPHQR4 + R_NADPHQR4 = 0
M_2dmmql8_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M - F_FRD3 + R_FRD3 + F_G3PD7 - R_G3PD7 + F_GLYCTO4 - R_GLYCTO4 + F_NADH18pp - R_NADH18pp + F_NADH9 - R_NADH9 + F_NADPHQR4 - R_NADPHQR4 + F_LMPD_5_2dmmql8_c - R_LMPD_5_2dmmql8_c = 0
M_2fe2s_c: -2.5e-05 F_Ec_biomass_iJO1366_WT_53p95M +2.5e-05 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_6_2fe2s_c - R_LMPD_6_2fe2s_c + F_LMPD_7_2fe2s_c - R_LMPD_7_2fe2s_c + F_LMPD_8_2fe2s_c - R_LMPD_8_2fe2s_c + F_LMPD_9_2fe2s_c - R_LMPD_9_2fe2s_c + F_LMPD_10_2fe2s_c - R_LMPD_10_2fe2s_c + F_LMPD_11_2fe2s_c - R_LMPD_11_2fe2s_c + F_LMPD_12_2fe2s_c - R_LMPD_12_2fe2s_c + F_LMPD_13_2fe2s_c - R_LMPD_13_2fe2s_c + F_LMPD_130_lipopb_c - R_LMPD_130_lipopb_c + F_LMPD_131_lipopb_c - R_LMPD_131_lipopb_c + F_LMPD_132_lipopb_c - R_LMPD_132_lipopb_c + F_LMPD_133_lipopb_c - R_LMPD_133_lipopb_c = 0
M_2h3oppan_c: F_GLXCL - R_GLXCL - F_TRSARr + R_TRSARr = 0
M_2pg_c: - F_ENO + R_ENO + F_GLYCK2 - R_GLYCK2 - F_PGM + R_PGM = 0
M_3pg_c: F_GLYCK - R_GLYCK - F_PGK + R_PGK + F_PGM - R_PGM -2 F_LMPD_1_10fthf_c +2 R_LMPD_1_10fthf_c -2 F_LMPD_2_10fthf_c +2 R_LMPD_2_10fthf_c -3 F_LMPD_3_10fthf_c +3 R_LMPD_3_10fthf_c -3 F_LMPD_4_10fthf_c +3 R_LMPD_4_10fthf_c -2 F_LMPD_6_2fe2s_c +2 R_LMPD_6_2fe2s_c -2 F_LMPD_7_2fe2s_c +2 R_LMPD_7_2fe2s_c -2 F_LMPD_8_2fe2s_c +2 R_LMPD_8_2fe2s_c -2 F_LMPD_9_2fe2s_c +2 R_LMPD_9_2fe2s_c -2 F_LMPD_10_2fe2s_c +2 R_LMPD_10_2fe2s_c -2 F_LMPD_11_2fe2s_c +2 R_LMPD_11_2fe2s_c -2 F_LMPD_12_2fe2s_c +2 R_LMPD_12_2fe2s_c -2 F_LMPD_13_2fe2s_c +2 R_LMPD_13_2fe2s_c -4 F_LMPD_14_4fe4s_c +4 R_LMPD_14_4fe4s_c -4 F_LMPD_15_4fe4s_c +4 R_LMPD_15_4fe4s_c -4 F_LMPD_16_4fe4s_c +4 R_LMPD_16_4fe4s_c -4 F_LMPD_17_4fe4s_c +4 R_LMPD_17_4fe4s_c -4 F_LMPD_18_4fe4s_c +4 R_LMPD_18_4fe4s_c -4 F_LMPD_19_4fe4s_c +4 R_LMPD_19_4fe4s_c -4 F_LMPD_20_4fe4s_c +4 R_LMPD_20_4fe4s_c -4 F_LMPD_21_4fe4s_c +4 R_LMPD_21_4fe4s_c - F_LMPD_22_5mthf_c + R_LMPD_22_5mthf_c - F_LMPD_23_5mthf_c + R_LMPD_23_5mthf_c -2 F_LMPD_24_accoa_c +2 R_LMPD_24_accoa_c -2 F_LMPD_25_accoa_c +2 R_LMPD_25_accoa_c -2 F_LMPD_26_adocbl_c +2 R_LMPD_26_adocbl_c -2 F_LMPD_27_adocbl_c +2 R_LMPD_27_adocbl_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c -2 F_LMPD_31_amet_c +2 R_LMPD_31_amet_c -2 F_LMPD_32_amet_c +2 R_LMPD_32_amet_c -12 F_LMPD_35_bmocogdp_c +12 R_LMPD_35_bmocogdp_c -12 F_LMPD_36_bmocogdp_c +12 R_LMPD_36_bmocogdp_c -8 F_LMPD_37_bmocogdp_c +8 R_LMPD_37_bmocogdp_c -8 F_LMPD_38_bmocogdp_c +8 R_LMPD_38_bmocogdp_c -4 F_LMPD_39_btn_c +4 R_LMPD_39_btn_c -4 F_LMPD_40_btn_c +4 R_LMPD_40_btn_c -4 F_LMPD_41_btn_c +4 R_LMPD_41_btn_c -4 F_LMPD_42_btn_c +4 R_LMPD_42_btn_c -6 F_LMPD_43_btn_c +6 R_LMPD_43_btn_c -6 F_LMPD_44_btn_c +6 R_LMPD_44_btn_c -6 F_LMPD_45_btn_c +6 R_LMPD_45_btn_c -6 F_LMPD_47_btn_c +6 R_LMPD_47_btn_c - F_LMPD_48_ca2_c + R_LMPD_48_ca2_c - F_LMPD_49_ca2_c + R_LMPD_49_ca2_c - F_LMPD_92_cys-L_c + R_LMPD_92_cys-L_c - F_LMPD_93_cys-L_c + R_LMPD_93_cys-L_c -2 F_LMPD_94_datp_c +2 R_LMPD_94_datp_c -2 F_LMPD_95_datp_c +2 R_LMPD_95_datp_c - F_LMPD_96_datp_c + R_LMPD_96_datp_c - F_LMPD_97_datp_c + R_LMPD_97_datp_c -2 F_LMPD_99_dgtp_c +2 R_LMPD_99_dgtp_c -2 F_LMPD_100_dgtp_c +2 R_LMPD_100_dgtp_c -2 F_LMPD_101_dgtp_c +2 R_LMPD_101_dgtp_c - F_LMPD_102_dgtp_c + R_LMPD_102_dgtp_c - F_LMPD_103_dgtp_c + R_LMPD_103_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c -3 F_LMPD_107_enter_c +3 R_LMPD_107_enter_c -2 F_LMPD_108_fad_c +2 R_LMPD_108_fad_c -2 F_LMPD_109_fad_c +2 R_LMPD_109_fad_c -4 F_LMPD_110_fad_c +4 R_LMPD_110_fad_c -4 F_LMPD_111_fad_c +4 R_LMPD_111_fad_c -2 F_LMPD_112_GAM_c +2 R_LMPD_112_GAM_c -2 F_LMPD_113_GAM_c +2 R_LMPD_113_GAM_c - F_LMPD_114_GAM_c + R_LMPD_114_GAM_c - F_LMPD_115_GAM_c + R_LMPD_115_GAM_c -2 F_LMPD_118_gthrd_c +2 R_LMPD_118_gthrd_c -2 F_LMPD_119_gthrd_c +2 R_LMPD_119_gthrd_c -2 F_LMPD_120_gtp_c +2 R_LMPD_120_gtp_c -2 F_LMPD_121_gtp_c +2 R_LMPD_121_gtp_c - F_LMPD_122_gtp_c + R_LMPD_122_gtp_c - F_LMPD_123_gtp_c + R_LMPD_123_gtp_c -4 F_LMPD_130_lipopb_c +4 R_LMPD_130_lipopb_c -4 F_LMPD_131_lipopb_c +4 R_LMPD_131_lipopb_c -4 F_LMPD_132_lipopb_c +4 R_LMPD_132_lipopb_c -4 F_LMPD_133_lipopb_c +4 R_LMPD_133_lipopb_c -2 F_LMPD_135_malcoa_c +2 R_LMPD_135_malcoa_c -2 F_LMPD_136_malcoa_c +2 R_LMPD_136_malcoa_c -2 F_LMPD_137_met-L_c +2 R_LMPD_137_met-L_c -2 F_LMPD_138_met-L_c +2 R_LMPD_138_met-L_c - F_LMPD_139_mlthf_c + R_LMPD_139_mlthf_c - F_LMPD_140_mlthf_c + R_LMPD_140_mlthf_c -3 F_LMPD_141_mococdp_c +3 R_LMPD_141_mococdp_c -3 F_LMPD_142_mococdp_c +3 R_LMPD_142_mococdp_c -3 F_LMPD_143_mococdp_c +3 R_LMPD_143_mococdp_c -3 F_LMPD_144_mococdp_c +3 R_LMPD_144_mococdp_c -4 F_LMPD_145_mococdp_c +4 R_LMPD_145_mococdp_c -4 F_LMPD_146_mococdp_c +4 R_LMPD_146_mococdp_c -4 F_LMPD_147_mococdp_c +4 R_LMPD_147_mococdp_c -4 F_LMPD_148_mococdp_c +4 R_LMPD_148_mococdp_c -4 F_LMPD_149_mocogdp_c +4 R_LMPD_149_mocogdp_c -4 F_LMPD_150_mocogdp_c +4 R_LMPD_150_mocogdp_c -6 F_LMPD_151_mocogdp_c +6 R_LMPD_151_mocogdp_c -6 F_LMPD_152_mocogdp_c +6 R_LMPD_152_mocogdp_c -2 F_LMPD_159_nad_c +2 R_LMPD_159_nad_c -2 F_LMPD_160_nad_c +2 R_LMPD_160_nad_c - F_LMPD_161_nad_c + R_LMPD_161_nad_c - F_LMPD_162_nad_c + R_LMPD_162_nad_c - F_LMPD_163_pe160_c + R_LMPD_163_pe160_c - F_LMPD_164_pe160_c + R_LMPD_164_pe160_c - F_LMPD_165_pe160_c + R_LMPD_165_pe160_c - F_LMPD_166_pe160_c + R_LMPD_166_pe160_c - F_LMPD_167_pe160_c + R_LMPD_167_pe160_c - F_LMPD_168_pe160_c + R_LMPD_168_pe160_c - F_LMPD_169_pe160_c + R_LMPD_169_pe160_c - F_LMPD_170_pe160_c + R_LMPD_170_pe160_c - F_LMPD_171_pe160_c + R_LMPD_171_pe160_c - F_LMPD_172_pe160_c + R_LMPD_172_pe160_c - F_LMPD_173_pe161_c + R_LMPD_173_pe161_c - F_LMPD_174_pe161_c + R_LMPD_174_pe161_c - F_LMPD_175_pe161_c + R_LMPD_175_pe161_c - F_LMPD_176_pe161_c + R_LMPD_176_pe161_c - F_LMPD_177_pe161_c + R_LMPD_177_pe161_c - F_LMPD_178_pe161_c + R_LMPD_178_pe161_c - F_LMPD_179_pe161_c + R_LMPD_179_pe161_c - F_LMPD_180_pe181_c + R_LMPD_180_pe181_c - F_LMPD_181_pe181_c + R_LMPD_181_pe181_c - F_LMPD_182_pe181_c + R_LMPD_182_pe181_c - F_LMPD_183_pe181_c + R_LMPD_183_pe181_c - F_LMPD_184_pe181_c + R_LMPD_184_pe181_c - F_LMPD_185_pe181_c + R_LMPD_185_pe181_c - F_LMPD_186_pe181_c + R_LMPD_186_pe181_c - F_LMPD_187_pe181_c + R_LMPD_187_pe181_c - F_LMPD_188_pe181_c + R_LMPD_188_pe181_c - F_LMPD_189_pe181_c + R_LMPD_189_pe181_c -3 F_LMPD_224_q8h2_c +3 R_LMPD_224_q8h2_c -3 F_LMPD_225_q8h2_c +3 R_LMPD_225_q8h2_c -3 F_LMPD_226_q8h2_c +3 R_LMPD_226_q8h2_c -6 F_LMPD_232_q8h2_c +6 R_LMPD_232_q8h2_c -2 F_LMPD_235_ribflv_c +2 R_LMPD_235_ribflv_c -2 F_LMPD_236_ribflv_c +2 R_LMPD_236_ribflv_c - F_LMPD_237_ribflv_c + R_LMPD_237_ribflv_c - F_LMPD_238_ribflv_c + R_LMPD_238_ribflv_c - F_LMPD_239_ser-L_c + R_LMPD_239_ser-L_c -2 F_LMPD_241_spmd_c +2 R_LMPD_241_spmd_c -2 F_LMPD_242_spmd_c +2 R_LMPD_242_spmd_c -2 F_LMPD_243_thf_c +2 R_LMPD_243_thf_c -2 F_LMPD_244_thf_c +2 R_LMPD_244_thf_c - F_LMPD_245_thf_c + R_LMPD_245_thf_c - F_LMPD_246_thf_c + R_LMPD_246_thf_c -2 F_LMPD_247_thmpp_c +2 R_LMPD_247_thmpp_c -2 F_LMPD_248_thmpp_c +2 R_LMPD_248_thmpp_c = 0
M_4crsol_c: - F_DM_4CRSOL + R_DM_4CRSOL + F_LMPD_247_thmpp_c - R_LMPD_247_thmpp_c + F_LMPD_248_thmpp_c - R_LMPD_248_thmpp_c = 0
M_4fe4s_c: -0.000248 F_Ec_biomass_iJO1366_WT_53p95M +0.000248 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_14_4fe4s_c - R_LMPD_14_4fe4s_c + F_LMPD_15_4fe4s_c - R_LMPD_15_4fe4s_c + F_LMPD_16_4fe4s_c - R_LMPD_16_4fe4s_c + F_LMPD_17_4fe4s_c - R_LMPD_17_4fe4s_c + F_LMPD_18_4fe4s_c - R_LMPD_18_4fe4s_c + F_LMPD_19_4fe4s_c - R_LMPD_19_4fe4s_c + F_LMPD_20_4fe4s_c - R_LMPD_20_4fe4s_c + F_LMPD_21_4fe4s_c - R_LMPD_21_4fe4s_c = 0
M_5drib_c: - F_DM_5DRIB + R_DM_5DRIB + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c + F_LMPD_43_btn_c - R_LMPD_43_btn_c + F_LMPD_44_btn_c - R_LMPD_44_btn_c + F_LMPD_45_btn_c - R_LMPD_45_btn_c + F_LMPD_47_btn_c - R_LMPD_47_btn_c +2 F_LMPD_125_hemeO_c -2 R_LMPD_125_hemeO_c +2 F_LMPD_130_lipopb_c -2 R_LMPD_130_lipopb_c +2 F_LMPD_131_lipopb_c -2 R_LMPD_131_lipopb_c +2 F_LMPD_132_lipopb_c -2 R_LMPD_132_lipopb_c +2 F_LMPD_133_lipopb_c -2 R_LMPD_133_lipopb_c +2 F_LMPD_219_pheme_c -2 R_LMPD_219_pheme_c + F_LMPD_247_thmpp_c - R_LMPD_247_thmpp_c + F_LMPD_248_thmpp_c - R_LMPD_248_thmpp_c = 0
M_5mthf_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_22_5mthf_c - R_LMPD_22_5mthf_c + F_LMPD_23_5mthf_c - R_LMPD_23_5mthf_c = 0
M_6pgc_c: - F_EDD + R_EDD - F_GND + R_GND + F_PGL - R_PGL = 0
M_6pgl_c: F_G6PDH2r - R_G6PDH2r - F_PGL + R_PGL = 0
M_ac_c: - F_ACKr + R_ACKr - F_ACS + R_ACS + F_ACt2rpp - R_ACt2rpp + F_ALDD2x - R_ALDD2x + F_ALDD2y - R_ALDD2y + F_CITL - R_CITL + F_POX - R_POX +2 F_LMPD_6_2fe2s_c -2 R_LMPD_6_2fe2s_c +2 F_LMPD_7_2fe2s_c -2 R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c +2 F_LMPD_9_2fe2s_c -2 R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c +2 F_LMPD_11_2fe2s_c -2 R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +4 F_LMPD_14_4fe4s_c -4 R_LMPD_14_4fe4s_c +4 F_LMPD_15_4fe4s_c -4 R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c +4 F_LMPD_17_4fe4s_c -4 R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c +4 F_LMPD_19_4fe4s_c -4 R_LMPD_19_4fe4s_c +4 F_LMPD_20_4fe4s_c -4 R_LMPD_20_4fe4s_c +4 F_LMPD_21_4fe4s_c -4 R_LMPD_21_4fe4s_c + F_LMPD_24_accoa_c - R_LMPD_24_accoa_c + F_LMPD_25_accoa_c - R_LMPD_25_accoa_c + F_LMPD_31_amet_c - R_LMPD_31_amet_c + F_LMPD_32_amet_c - R_LMPD_32_amet_c + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c +4 F_LMPD_35_bmocogdp_c -4 R_LMPD_35_bmocogdp_c +4 F_LMPD_36_bmocogdp_c -4 R_LMPD_36_bmocogdp_c +4 F_LMPD_37_bmocogdp_c -4 R_LMPD_37_bmocogdp_c +4 F_LMPD_38_bmocogdp_c -4 R_LMPD_38_bmocogdp_c +2 F_LMPD_39_btn_c -2 R_LMPD_39_btn_c +2 F_LMPD_40_btn_c -2 R_LMPD_40_btn_c +2 F_LMPD_41_btn_c -2 R_LMPD_41_btn_c +2 F_LMPD_42_btn_c -2 R_LMPD_42_btn_c +3 F_LMPD_43_btn_c -3 R_LMPD_43_btn_c +3 F_LMPD_44_btn_c -3 R_LMPD_44_btn_c +3 F_LMPD_45_btn_c -3 R_LMPD_45_btn_c +3 F_LMPD_47_btn_c -3 R_LMPD_47_btn_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c +2 F_LMPD_79_colipa_e -2 R_LMPD_79_colipa_e +2 F_LMPD_80_colipa_e -2 R_LMPD_80_colipa_e +2 F_LMPD_81_colipa_e -2 R_LMPD_81_colipa_e +2 F_LMPD_82_colipa_e -2 R_LMPD_82_colipa_e +2 F_LMPD_83_colipa_e -2 R_LMPD_83_colipa_e +2 F_LMPD_84_colipa_e -2 R_LMPD_84_colipa_e +2 F_LMPD_85_colipa_e -2 R_LMPD_85_colipa_e +2 F_LMPD_86_colipa_e -2 R_LMPD_86_colipa_e +2 F_LMPD_87_colipa_e -2 R_LMPD_87_colipa_e +2 F_LMPD_88_colipa_e -2 R_LMPD_88_colipa_e +2 F_LMPD_89_colipa_e -2 R_LMPD_89_colipa_e +2 F_LMPD_90_colipa_e -2 R_LMPD_90_colipa_e + F_LMPD_92_cys-L_c - R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c + F_LMPD_118_gthrd_c - R_LMPD_118_gthrd_c + F_LMPD_119_gthrd_c - R_LMPD_119_gthrd_c +4 F_LMPD_130_lipopb_c -4 R_LMPD_130_lipopb_c +4 F_LMPD_131_lipopb_c -4 R_LMPD_131_lipopb_c +4 F_LMPD_132_lipopb_c -4 R_LMPD_132_lipopb_c +4 F_LMPD_133_lipopb_c -4 R_LMPD_133_lipopb_c + F_LMPD_135_malcoa_c - R_LMPD_135_malcoa_c + F_LMPD_136_malcoa_c - R_LMPD_136_malcoa_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c +2 F_LMPD_141_mococdp_c -2 R_LMPD_141_mococdp_c +2 F_LMPD_142_mococdp_c -2 R_LMPD_142_mococdp_c +2 F_LMPD_143_mococdp_c -2 R_LMPD_143_mococdp_c +2 F_LMPD_144_mococdp_c -2 R_LMPD_144_mococdp_c +2 F_LMPD_145_mococdp_c -2 R_LMPD_145_mococdp_c +2 F_LMPD_146_mococdp_c -2 R_LMPD_146_mococdp_c +2 F_LMPD_147_mococdp_c -2 R_LMPD_147_mococdp_c +2 F_LMPD_148_mococdp_c -2 R_LMPD_148_mococdp_c +2 F_LMPD_149_mocogdp_c -2 R_LMPD_149_mocogdp_c +2 F_LMPD_150_mocogdp_c -2 R_LMPD_150_mocogdp_c +2 F_LMPD_151_mocogdp_c -2 R_LMPD_151_mocogdp_c +2 F_LMPD_152_mocogdp_c -2 R_LMPD_152_mocogdp_c + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c +2 F_LMPD_241_spmd_c -2 R_LMPD_241_spmd_c +2 F_LMPD_242_spmd_c -2 R_LMPD_242_spmd_c + F_LMPD_247_thmpp_c - R_LMPD_247_thmpp_c + F_LMPD_248_thmpp_c - R_LMPD_248_thmpp_c = 0
M_acald_c: - F_ACALD + R_ACALD + F_ACALDtpp - R_ACALDtpp + F_ALCD2x - R_ALCD2x - F_ALDD2x + R_ALDD2x - F_ALDD2y + R_ALDD2y + F_LMPD_117_gly_c - R_LMPD_117_gly_c = 0
M_accoa_c: -0.000279 F_Ec_biomass_iJO1366_WT_53p95M +0.000279 R_Ec_biomass_iJO1366_WT_53p95M + F_ACALD - R_ACALD + F_ACS - R_ACS - F_CS + R_CS - F_MALS + R_MALS + F_PDH - R_PDH + F_PFL - R_PFL + F_POR5 - R_POR5 - F_PTAr + R_PTAr -2 F_LMPD_6_2fe2s_c +2 R_LMPD_6_2fe2s_c -2 F_LMPD_7_2fe2s_c +2 R_LMPD_7_2fe2s_c -2 F_LMPD_8_2fe2s_c +2 R_LMPD_8_2fe2s_c -2 F_LMPD_9_2fe2s_c +2 R_LMPD_9_2fe2s_c -2 F_LMPD_10_2fe2s_c +2 R_LMPD_10_2fe2s_c -2 F_LMPD_11_2fe2s_c +2 R_LMPD_11_2fe2s_c -2 F_LMPD_12_2fe2s_c +2 R_LMPD_12_2fe2s_c -2 F_LMPD_13_2fe2s_c +2 R_LMPD_13_2fe2s_c -4 F_LMPD_14_4fe4s_c +4 R_LMPD_14_4fe4s_c -4 F_LMPD_15_4fe4s_c +4 R_LMPD_15_4fe4s_c -4 F_LMPD_16_4fe4s_c +4 R_LMPD_16_4fe4s_c -4 F_LMPD_17_4fe4s_c +4 R_LMPD_17_4fe4s_c -4 F_LMPD_18_4fe4s_c +4 R_LMPD_18_4fe4s_c -4 F_LMPD_19_4fe4s_c +4 R_LMPD_19_4fe4s_c -4 F_LMPD_20_4fe4s_c +4 R_LMPD_20_4fe4s_c -4 F_LMPD_21_4fe4s_c +4 R_LMPD_21_4fe4s_c - F_LMPD_24_accoa_c + R_LMPD_24_accoa_c - F_LMPD_25_accoa_c + R_LMPD_25_accoa_c - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c - F_LMPD_33_arg-L_c + R_LMPD_33_arg-L_c -4 F_LMPD_35_bmocogdp_c +4 R_LMPD_35_bmocogdp_c -4 F_LMPD_36_bmocogdp_c +4 R_LMPD_36_bmocogdp_c -4 F_LMPD_37_bmocogdp_c +4 R_LMPD_37_bmocogdp_c -4 F_LMPD_38_bmocogdp_c +4 R_LMPD_38_bmocogdp_c -5 F_LMPD_39_btn_c +5 R_LMPD_39_btn_c -5 F_LMPD_40_btn_c +5 R_LMPD_40_btn_c -5 F_LMPD_41_btn_c +5 R_LMPD_41_btn_c -5 F_LMPD_42_btn_c +5 R_LMPD_42_btn_c -6 F_LMPD_43_btn_c +6 R_LMPD_43_btn_c -6 F_LMPD_44_btn_c +6 R_LMPD_44_btn_c -6 F_LMPD_45_btn_c +6 R_LMPD_45_btn_c -6 F_LMPD_47_btn_c +6 R_LMPD_47_btn_c - F_LMPD_48_ca2_c + R_LMPD_48_ca2_c - F_LMPD_49_ca2_c + R_LMPD_49_ca2_c -32 F_LMPD_52_clpn160_p +32 R_LMPD_52_clpn160_p -32 F_LMPD_53_clpn160_p +32 R_LMPD_53_clpn160_p -32 F_LMPD_54_clpn160_p +32 R_LMPD_54_clpn160_p -32 F_LMPD_55_clpn160_p +32 R_LMPD_55_clpn160_p -32 F_LMPD_56_clpn160_p +32 R_LMPD_56_clpn160_p -32 F_LMPD_57_clpn160_p +32 R_LMPD_57_clpn160_p -32 F_LMPD_58_clpn160_p +32 R_LMPD_58_clpn160_p -32 F_LMPD_59_clpn160_p +32 R_LMPD_59_clpn160_p -32 F_LMPD_60_clpn160_p +32 R_LMPD_60_clpn160_p -32 F_LMPD_61_clpn160_p +32 R_LMPD_61_clpn160_p -32 F_LMPD_62_clpn161_p +32 R_LMPD_62_clpn161_p -32 F_LMPD_63_clpn161_p +32 R_LMPD_63_clpn161_p -32 F_LMPD_64_clpn161_p +32 R_LMPD_64_clpn161_p -32 F_LMPD_65_clpn161_p +32 R_LMPD_65_clpn161_p -32 F_LMPD_66_clpn161_p +32 R_LMPD_66_clpn161_p -32 F_LMPD_67_clpn161_p +32 R_LMPD_67_clpn161_p -32 F_LMPD_68_clpn161_p +32 R_LMPD_68_clpn161_p -36 F_LMPD_69_clpn181_p +36 R_LMPD_69_clpn181_p -36 F_LMPD_70_clpn181_p +36 R_LMPD_70_clpn181_p -36 F_LMPD_71_clpn181_p +36 R_LMPD_71_clpn181_p -36 F_LMPD_72_clpn181_p +36 R_LMPD_72_clpn181_p -36 F_LMPD_73_clpn181_p +36 R_LMPD_73_clpn181_p -36 F_LMPD_74_clpn181_p +36 R_LMPD_74_clpn181_p -36 F_LMPD_75_clpn181_p +36 R_LMPD_75_clpn181_p -36 F_LMPD_76_clpn181_p +36 R_LMPD_76_clpn181_p -36 F_LMPD_77_clpn181_p +36 R_LMPD_77_clpn181_p -36 F_LMPD_78_clpn181_p +36 R_LMPD_78_clpn181_p -43 F_LMPD_79_colipa_e +43 R_LMPD_79_colipa_e -43 F_LMPD_80_colipa_e +43 R_LMPD_80_colipa_e -43 F_LMPD_81_colipa_e +43 R_LMPD_81_colipa_e -43 F_LMPD_82_colipa_e +43 R_LMPD_82_colipa_e -43 F_LMPD_83_colipa_e +43 R_LMPD_83_colipa_e -43 F_LMPD_84_colipa_e +43 R_LMPD_84_colipa_e -43 F_LMPD_85_colipa_e +43 R_LMPD_85_colipa_e -43 F_LMPD_86_colipa_e +43 R_LMPD_86_colipa_e -43 F_LMPD_87_colipa_e +43 R_LMPD_87_colipa_e -43 F_LMPD_88_colipa_e +43 R_LMPD_88_colipa_e -43 F_LMPD_89_colipa_e +43 R_LMPD_89_colipa_e -43 F_LMPD_90_colipa_e +43 R_LMPD_90_colipa_e - F_LMPD_92_cys-L_c + R_LMPD_92_cys-L_c - F_LMPD_93_cys-L_c + R_LMPD_93_cys-L_c - F_LMPD_118_gthrd_c + R_LMPD_118_gthrd_c - F_LMPD_119_gthrd_c + R_LMPD_119_gthrd_c - F_LMPD_129_leu-L_c + R_LMPD_129_leu-L_c -8 F_LMPD_130_lipopb_c +8 R_LMPD_130_lipopb_c -8 F_LMPD_131_lipopb_c +8 R_LMPD_131_lipopb_c -8 F_LMPD_132_lipopb_c +8 R_LMPD_132_lipopb_c -8 F_LMPD_133_lipopb_c +8 R_LMPD_133_lipopb_c -2 F_LMPD_135_malcoa_c +2 R_LMPD_135_malcoa_c -2 F_LMPD_136_malcoa_c +2 R_LMPD_136_malcoa_c - F_LMPD_137_met-L_c + R_LMPD_137_met-L_c - F_LMPD_138_met-L_c + R_LMPD_138_met-L_c -2 F_LMPD_141_mococdp_c +2 R_LMPD_141_mococdp_c -2 F_LMPD_142_mococdp_c +2 R_LMPD_142_mococdp_c -2 F_LMPD_143_mococdp_c +2 R_LMPD_143_mococdp_c -2 F_LMPD_144_mococdp_c +2 R_LMPD_144_mococdp_c -2 F_LMPD_145_mococdp_c +2 R_LMPD_145_mococdp_c -2 F_LMPD_146_mococdp_c +2 R_LMPD_146_mococdp_c -2 F_LMPD_147_mococdp_c +2 R_LMPD_147_mococdp_c -2 F_LMPD_148_mococdp_c +2 R_LMPD_148_mococdp_c -2 F_LMPD_149_mocogdp_c +2 R_LMPD_149_mocogdp_c -2 F_LMPD_150_mocogdp_c +2 R_LMPD_150_mocogdp_c -2 F_LMPD_151_mocogdp_c +2 R_LMPD_151_mocogdp_c -2 F_LMPD_152_mocogdp_c +2 R_LMPD_152_mocogdp_c -4 F_LMPD_154_murein3p3p_p +4 R_LMPD_154_murein3p3p_p -4 F_LMPD_155_murein3px4p_p +4 R_LMPD_155_murein3px4p_p -4 F_LMPD_156_murein4p4p_p +4 R_LMPD_156_murein4p4p_p -4 F_LMPD_157_murein4px4p_p +4 R_LMPD_157_murein4px4p_p -6 F_LMPD_158_murein4px4px4p_p +6 R_LMPD_158_murein4px4px4p_p -16 F_LMPD_163_pe160_c +16 R_LMPD_163_pe160_c -16 F_LMPD_164_pe160_c +16 R_LMPD_164_pe160_c -16 F_LMPD_165_pe160_c +16 R_LMPD_165_pe160_c -16 F_LMPD_166_pe160_c +16 R_LMPD_166_pe160_c -16 F_LMPD_167_pe160_c +16 R_LMPD_167_pe160_c -16 F_LMPD_168_pe160_c +16 R_LMPD_168_pe160_c -16 F_LMPD_169_pe160_c +16 R_LMPD_169_pe160_c -16 F_LMPD_170_pe160_c +16 R_LMPD_170_pe160_c -16 F_LMPD_171_pe160_c +16 R_LMPD_171_pe160_c -16 F_LMPD_172_pe160_c +16 R_LMPD_172_pe160_c -16 F_LMPD_173_pe161_c +16 R_LMPD_173_pe161_c -16 F_LMPD_174_pe161_c +16 R_LMPD_174_pe161_c -16 F_LMPD_175_pe161_c +16 R_LMPD_175_pe161_c -16 F_LMPD_176_pe161_c +16 R_LMPD_176_pe161_c -16 F_LMPD_177_pe161_c +16 R_LMPD_177_pe161_c -16 F_LMPD_178_pe161_c +16 R_LMPD_178_pe161_c -16 F_LMPD_179_pe161_c +16 R_LMPD_179_pe161_c -18 F_LMPD_180_pe181_c +18 R_LMPD_180_pe181_c -18 F_LMPD_181_pe181_c +18 R_LMPD_181_pe181_c -18 F_LMPD_182_pe181_c +18 R_LMPD_182_pe181_c -18 F_LMPD_183_pe181_c +18 R_LMPD_183_pe181_c -18 F_LMPD_184_pe181_c +18 R_LMPD_184_pe181_c -18 F_LMPD_185_pe181_c +18 R_LMPD_185_pe181_c -18 F_LMPD_186_pe181_c +18 R_LMPD_186_pe181_c -18 F_LMPD_187_pe181_c +18 R_LMPD_187_pe181_c -18 F_LMPD_188_pe181_c +18 R_LMPD_188_pe181_c -18 F_LMPD_189_pe181_c +18 R_LMPD_189_pe181_c -16 F_LMPD_190_pg160_c +16 R_LMPD_190_pg160_c -16 F_LMPD_191_pg160_c +16 R_LMPD_191_pg160_c -16 F_LMPD_192_pg160_c +16 R_LMPD_192_pg160_c -16 F_LMPD_193_pg160_c +16 R_LMPD_193_pg160_c -16 F_LMPD_194_pg160_c +16 R_LMPD_194_pg160_c -16 F_LMPD_195_pg160_c +16 R_LMPD_195_pg160_c -16 F_LMPD_196_pg160_c +16 R_LMPD_196_pg160_c -16 F_LMPD_197_pg160_c +16 R_LMPD_197_pg160_c -16 F_LMPD_198_pg160_c +16 R_LMPD_198_pg160_c -16 F_LMPD_199_pg160_c +16 R_LMPD_199_pg160_c -16 F_LMPD_200_pg161_c +16 R_LMPD_200_pg161_c -16 F_LMPD_201_pg161_c +16 R_LMPD_201_pg161_c -16 F_LMPD_202_pg161_c +16 R_LMPD_202_pg161_c -16 F_LMPD_203_pg161_c +16 R_LMPD_203_pg161_c -16 F_LMPD_204_pg161_c +16 R_LMPD_204_pg161_c -16 F_LMPD_205_pg161_c +16 R_LMPD_205_pg161_c -16 F_LMPD_206_pg161_c +16 R_LMPD_206_pg161_c -18 F_LMPD_207_pg181_c +18 R_LMPD_207_pg181_c -18 F_LMPD_208_pg181_c +18 R_LMPD_208_pg181_c -18 F_LMPD_209_pg181_c +18 R_LMPD_209_pg181_c -18 F_LMPD_210_pg181_c +18 R_LMPD_210_pg181_c -18 F_LMPD_211_pg181_c +18 R_LMPD_211_pg181_c -18 F_LMPD_212_pg181_c +18 R_LMPD_212_pg181_c -18 F_LMPD_213_pg181_c +18 R_LMPD_213_pg181_c -18 F_LMPD_214_pg181_c +18 R_LMPD_214_pg181_c -18 F_LMPD_215_pg181_c +18 R_LMPD_215_pg181_c -18 F_LMPD_216_pg181_c +18 R_LMPD_216_pg181_c - F_LMPD_221_ptrc_c + R_LMPD_221_ptrc_c -2 F_LMPD_241_spmd_c +2 R_LMPD_241_spmd_c -2 F_LMPD_242_spmd_c +2 R_LMPD_242_spmd_c - F_LMPD_247_thmpp_c + R_LMPD_247_thmpp_c - F_LMPD_248_thmpp_c + R_LMPD_248_thmpp_c = 0
M_acon-C_c: F_ACONTa - R_ACONTa - F_ACONTb + R_ACONTb = 0
M_actp_c: F_ACKr - R_ACKr + F_PTAr - R_PTAr = 0
M_adocbl_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_26_adocbl_c - R_LMPD_26_adocbl_c + F_LMPD_27_adocbl_c - R_LMPD_27_adocbl_c + F_LMPD_28_adocbl_c - R_LMPD_28_adocbl_c + F_LMPD_29_adocbl_c - R_LMPD_29_adocbl_c = 0
M_adp_c: +53.95 F_Ec_biomass_iJO1366_WT_53p95M -53.95 R_Ec_biomass_iJO1366_WT_53p95M + F_ACKr - R_ACKr +2 F_ADK1 -2 R_ADK1 + F_ATPM - R_ATPM - F_ATPS4rpp + R_ATPS4rpp + F_COLIPAabcpp - R_COLIPAabcpp + F_COLIPAabctex - R_COLIPAabctex + F_CYSabc2pp - R_CYSabc2pp + F_FE2abcpp - R_FE2abcpp + F_FE3abcpp - R_FE3abcpp + F_GLCS1 - R_GLCS1 + F_GLCabcpp - R_GLCabcpp + F_GLUabcpp - R_GLUabcpp + F_GLYCK - R_GLYCK + F_GLYCK2 - R_GLYCK2 + F_GLYK - R_GLYK + F_GTHRDabc2pp - R_GTHRDabc2pp + F_HEX1 - R_HEX1 + F_Kabcpp - R_Kabcpp + F_MG2uabcpp - R_MG2uabcpp + F_MOBDabcpp - R_MOBDabcpp + F_NADK - R_NADK + F_NI2uabcpp - R_NI2uabcpp + F_NTP1 - R_NTP1 + F_PE160abcpp - R_PE160abcpp + F_PE161abcpp - R_PE161abcpp + F_PE181abcpp - R_PE181abcpp + F_PFK - R_PFK + F_PFK_3 - R_PFK_3 + F_PG160abcpp - R_PG160abcpp + F_PG161abcpp - R_PG161abcpp + F_PG181abcpp - R_PG181abcpp + F_PGK - R_PGK + F_PHEMEabcpp - R_PHEMEabcpp + F_PIuabcpp - R_PIuabcpp + F_PPCK - R_PPCK + F_PPK2r - R_PPK2r + F_PPKr - R_PPKr - F_PYK + R_PYK + F_SUCOAS - R_SUCOAS + F_SULabcpp - R_SULabcpp + F_ZNabcpp - R_ZNabcpp +12 F_LMPD_1_10fthf_c -12 R_LMPD_1_10fthf_c +14 F_LMPD_2_10fthf_c -14 R_LMPD_2_10fthf_c +11 F_LMPD_3_10fthf_c -11 R_LMPD_3_10fthf_c +13 F_LMPD_4_10fthf_c -13 R_LMPD_4_10fthf_c +25 F_LMPD_5_2dmmql8_c -25 R_LMPD_5_2dmmql8_c + F_LMPD_8_2fe2s_c - R_LMPD_8_2fe2s_c + F_LMPD_9_2fe2s_c - R_LMPD_9_2fe2s_c +4 F_LMPD_10_2fe2s_c -4 R_LMPD_10_2fe2s_c +4 F_LMPD_11_2fe2s_c -4 R_LMPD_11_2fe2s_c +5 F_LMPD_12_2fe2s_c -5 R_LMPD_12_2fe2s_c +5 F_LMPD_13_2fe2s_c -5 R_LMPD_13_2fe2s_c +2 F_LMPD_14_4fe4s_c -2 R_LMPD_14_4fe4s_c +2 F_LMPD_15_4fe4s_c -2 R_LMPD_15_4fe4s_c +8 F_LMPD_16_4fe4s_c -8 R_LMPD_16_4fe4s_c +8 F_LMPD_17_4fe4s_c -8 R_LMPD_17_4fe4s_c +10 F_LMPD_18_4fe4s_c -10 R_LMPD_18_4fe4s_c +10 F_LMPD_19_4fe4s_c -10 R_LMPD_19_4fe4s_c +13 F_LMPD_22_5mthf_c -13 R_LMPD_22_5mthf_c +15 F_LMPD_23_5mthf_c -15 R_LMPD_23_5mthf_c +12 F_LMPD_24_accoa_c -12 R_LMPD_24_accoa_c +16 F_LMPD_25_accoa_c -16 R_LMPD_25_accoa_c +7 F_LMPD_26_adocbl_c -7 R_LMPD_26_adocbl_c +9 F_LMPD_27_adocbl_c -9 R_LMPD_27_adocbl_c +8 F_LMPD_28_adocbl_c -8 R_LMPD_28_adocbl_c +10 F_LMPD_29_adocbl_c -10 R_LMPD_29_adocbl_c +9 F_LMPD_31_amet_c -9 R_LMPD_31_amet_c +13 F_LMPD_32_amet_c -13 R_LMPD_32_amet_c +2 F_LMPD_33_arg-L_c -2 R_LMPD_33_arg-L_c +32 F_LMPD_35_bmocogdp_c -32 R_LMPD_35_bmocogdp_c +48 F_LMPD_36_bmocogdp_c -48 R_LMPD_36_bmocogdp_c +36 F_LMPD_37_bmocogdp_c -36 R_LMPD_37_bmocogdp_c +52 F_LMPD_38_bmocogdp_c -52 R_LMPD_38_bmocogdp_c +13 F_LMPD_39_btn_c -13 R_LMPD_39_btn_c +14 F_LMPD_40_btn_c -14 R_LMPD_40_btn_c +19 F_LMPD_41_btn_c -19 R_LMPD_41_btn_c +20 F_LMPD_42_btn_c -20 R_LMPD_42_btn_c +30 F_LMPD_43_btn_c -30 R_LMPD_43_btn_c +31 F_LMPD_44_btn_c -31 R_LMPD_44_btn_c +42 F_LMPD_45_btn_c -42 R_LMPD_45_btn_c +43 F_LMPD_47_btn_c -43 R_LMPD_47_btn_c +3 F_LMPD_48_ca2_c -3 R_LMPD_48_ca2_c +3 F_LMPD_49_ca2_c -3 R_LMPD_49_ca2_c +3 F_LMPD_50_ca2_c -3 R_LMPD_50_ca2_c + F_LMPD_51_chor_c - R_LMPD_51_chor_c +32 F_LMPD_52_clpn160_p -32 R_LMPD_52_clpn160_p +32 F_LMPD_53_clpn160_p -32 R_LMPD_53_clpn160_p +32 F_LMPD_54_clpn160_p -32 R_LMPD_54_clpn160_p +32 F_LMPD_55_clpn160_p -32 R_LMPD_55_clpn160_p +32 F_LMPD_56_clpn160_p -32 R_LMPD_56_clpn160_p +32 F_LMPD_57_clpn160_p -32 R_LMPD_57_clpn160_p +32 F_LMPD_58_clpn160_p -32 R_LMPD_58_clpn160_p +32 F_LMPD_59_clpn160_p -32 R_LMPD_59_clpn160_p +4 F_LMPD_60_clpn160_p -4 R_LMPD_60_clpn160_p +4 F_LMPD_61_clpn160_p -4 R_LMPD_61_clpn160_p +32 F_LMPD_62_clpn161_p -32 R_LMPD_62_clpn161_p +32 F_LMPD_63_clpn161_p -32 R_LMPD_63_clpn161_p +32 F_LMPD_64_clpn161_p -32 R_LMPD_64_clpn161_p +32 F_LMPD_65_clpn161_p -32 R_LMPD_65_clpn161_p +32 F_LMPD_66_clpn161_p -32 R_LMPD_66_clpn161_p +32 F_LMPD_67_clpn161_p -32 R_LMPD_67_clpn161_p +32 F_LMPD_68_clpn161_p -32 R_LMPD_68_clpn161_p +36 F_LMPD_69_clpn181_p -36 R_LMPD_69_clpn181_p +36 F_LMPD_70_clpn181_p -36 R_LMPD_70_clpn181_p +36 F_LMPD_71_clpn181_p -36 R_LMPD_71_clpn181_p +36 F_LMPD_72_clpn181_p -36 R_LMPD_72_clpn181_p +36 F_LMPD_73_clpn181_p -36 R_LMPD_73_clpn181_p +36 F_LMPD_74_clpn181_p -36 R_LMPD_74_clpn181_p +36 F_LMPD_75_clpn181_p -36 R_LMPD_75_clpn181_p +36 F_LMPD_76_clpn181_p -36 R_LMPD_76_clpn181_p +4 F_LMPD_77_clpn181_p -4 R_LMPD_77_clpn181_p +4 F_LMPD_78_clpn181_p -4 R_LMPD_78_clpn181_p +62 F_LMPD_79_colipa_e -62 R_LMPD_79_colipa_e +62 F_LMPD_80_colipa_e -62 R_LMPD_80_colipa_e +62 F_LMPD_81_colipa_e -62 R_LMPD_81_colipa_e +62 F_LMPD_82_colipa_e -62 R_LMPD_82_colipa_e +62 F_LMPD_83_colipa_e -62 R_LMPD_83_colipa_e +62 F_LMPD_84_colipa_e -62 R_LMPD_84_colipa_e +62 F_LMPD_85_colipa_e -62 R_LMPD_85_colipa_e +62 F_LMPD_86_colipa_e -62 R_LMPD_86_colipa_e +62 F_LMPD_87_colipa_e -62 R_LMPD_87_colipa_e +62 F_LMPD_88_colipa_e -62 R_LMPD_88_colipa_e +62 F_LMPD_89_colipa_e -62 R_LMPD_89_colipa_e +62 F_LMPD_90_colipa_e -62 R_LMPD_90_colipa_e +5 F_LMPD_91_ctp_c -5 R_LMPD_91_ctp_c +2 F_LMPD_93_cys-L_c -2 R_LMPD_93_cys-L_c +8 F_LMPD_94_datp_c -8 R_LMPD_94_datp_c +6 F_LMPD_95_datp_c -6 R_LMPD_95_datp_c +7 F_LMPD_96_datp_c -7 R_LMPD_96_datp_c +9 F_LMPD_97_datp_c -9 R_LMPD_97_datp_c +5 F_LMPD_98_dctp_c -5 R_LMPD_98_dctp_c +8 F_LMPD_99_dgtp_c -8 R_LMPD_99_dgtp_c +10 F_LMPD_100_dgtp_c -10 R_LMPD_100_dgtp_c +10 F_LMPD_101_dgtp_c -10 R_LMPD_101_dgtp_c +9 F_LMPD_102_dgtp_c -9 R_LMPD_102_dgtp_c +11 F_LMPD_103_dgtp_c -11 R_LMPD_103_dgtp_c +11 F_LMPD_104_dgtp_c -11 R_LMPD_104_dgtp_c +6 F_LMPD_105_dttp_c -6 R_LMPD_105_dttp_c +6 F_LMPD_106_dttp_c -6 R_LMPD_106_dttp_c +3 F_LMPD_107_enter_c -3 R_LMPD_107_enter_c +17 F_LMPD_108_fad_c -17 R_LMPD_108_fad_c +21 F_LMPD_109_fad_c -21 R_LMPD_109_fad_c +15 F_LMPD_110_fad_c -15 R_LMPD_110_fad_c +19 F_LMPD_111_fad_c -19 R_LMPD_111_fad_c +6 F_LMPD_112_GAM_c -6 R_LMPD_112_GAM_c +8 F_LMPD_113_GAM_c -8 R_LMPD_113_GAM_c +9 F_LMPD_114_GAM_c -9 R_LMPD_114_GAM_c +7 F_LMPD_115_GAM_c -7 R_LMPD_115_GAM_c + F_LMPD_116_gln-L_c - R_LMPD_116_gln-L_c +2 F_LMPD_117_gly_c -2 R_LMPD_117_gly_c +2 F_LMPD_118_gthrd_c -2 R_LMPD_118_gthrd_c +4 F_LMPD_119_gthrd_c -4 R_LMPD_119_gthrd_c +8 F_LMPD_120_gtp_c -8 R_LMPD_120_gtp_c +10 F_LMPD_121_gtp_c -10 R_LMPD_121_gtp_c +9 F_LMPD_122_gtp_c -9 R_LMPD_122_gtp_c +11 F_LMPD_123_gtp_c -11 R_LMPD_123_gtp_c +9 F_LMPD_124_hemeO_c -9 R_LMPD_124_hemeO_c +9 F_LMPD_125_hemeO_c -9 R_LMPD_125_hemeO_c + F_LMPD_126_his-L_c - R_LMPD_126_his-L_c +3 F_LMPD_127_his-L_c -3 R_LMPD_127_his-L_c +2 F_LMPD_128_ile-L_c -2 R_LMPD_128_ile-L_c +8 F_LMPD_131_lipopb_c -8 R_LMPD_131_lipopb_c +10 F_LMPD_132_lipopb_c -10 R_LMPD_132_lipopb_c +8 F_LMPD_133_lipopb_c -8 R_LMPD_133_lipopb_c + F_LMPD_134_lys-L_c - R_LMPD_134_lys-L_c +13 F_LMPD_135_malcoa_c -13 R_LMPD_135_malcoa_c +17 F_LMPD_136_malcoa_c -17 R_LMPD_136_malcoa_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c +3 F_LMPD_138_met-L_c -3 R_LMPD_138_met-L_c +13 F_LMPD_139_mlthf_c -13 R_LMPD_139_mlthf_c +15 F_LMPD_140_mlthf_c -15 R_LMPD_140_mlthf_c +14 F_LMPD_141_mococdp_c -14 R_LMPD_141_mococdp_c +16 F_LMPD_142_mococdp_c -16 R_LMPD_142_mococdp_c +20 F_LMPD_143_mococdp_c -20 R_LMPD_143_mococdp_c +22 F_LMPD_144_mococdp_c -22 R_LMPD_144_mococdp_c +13 F_LMPD_145_mococdp_c -13 R_LMPD_145_mococdp_c +15 F_LMPD_146_mococdp_c -15 R_LMPD_146_mococdp_c +19 F_LMPD_147_mococdp_c -19 R_LMPD_147_mococdp_c +21 F_LMPD_148_mococdp_c -21 R_LMPD_148_mococdp_c +18 F_LMPD_149_mocogdp_c -18 R_LMPD_149_mocogdp_c +26 F_LMPD_150_mocogdp_c -26 R_LMPD_150_mocogdp_c +16 F_LMPD_151_mocogdp_c -16 R_LMPD_151_mocogdp_c +24 F_LMPD_152_mocogdp_c -24 R_LMPD_152_mocogdp_c +26 F_LMPD_153_mql8_c -26 R_LMPD_153_mql8_c +22 F_LMPD_154_murein3p3p_p -22 R_LMPD_154_murein3p3p_p +22 F_LMPD_155_murein3px4p_p -22 R_LMPD_155_murein3px4p_p +22 F_LMPD_156_murein4p4p_p -22 R_LMPD_156_murein4p4p_p +22 F_LMPD_157_murein4px4p_p -22 R_LMPD_157_murein4px4p_p +33 F_LMPD_158_murein4px4px4p_p -33 R_LMPD_158_murein4px4px4p_p +6 F_LMPD_159_nad_c -6 R_LMPD_159_nad_c +8 F_LMPD_160_nad_c -8 R_LMPD_160_nad_c +7 F_LMPD_161_nad_c -7 R_LMPD_161_nad_c +9 F_LMPD_162_nad_c -9 R_LMPD_162_nad_c +16 F_LMPD_163_pe160_c -16 R_LMPD_163_pe160_c +16 F_LMPD_164_pe160_c -16 R_LMPD_164_pe160_c +16 F_LMPD_165_pe160_c -16 R_LMPD_165_pe160_c +16 F_LMPD_166_pe160_c -16 R_LMPD_166_pe160_c +16 F_LMPD_167_pe160_c -16 R_LMPD_167_pe160_c +16 F_LMPD_168_pe160_c -16 R_LMPD_168_pe160_c +16 F_LMPD_169_pe160_c -16 R_LMPD_169_pe160_c +16 F_LMPD_170_pe160_c -16 R_LMPD_170_pe160_c +2 F_LMPD_171_pe160_c -2 R_LMPD_171_pe160_c +2 F_LMPD_172_pe160_c -2 R_LMPD_172_pe160_c +16 F_LMPD_173_pe161_c -16 R_LMPD_173_pe161_c +16 F_LMPD_174_pe161_c -16 R_LMPD_174_pe161_c +16 F_LMPD_175_pe161_c -16 R_LMPD_175_pe161_c +16 F_LMPD_176_pe161_c -16 R_LMPD_176_pe161_c +16 F_LMPD_177_pe161_c -16 R_LMPD_177_pe161_c +16 F_LMPD_178_pe161_c -16 R_LMPD_178_pe161_c +16 F_LMPD_179_pe161_c -16 R_LMPD_179_pe161_c +18 F_LMPD_180_pe181_c -18 R_LMPD_180_pe181_c +18 F_LMPD_181_pe181_c -18 R_LMPD_181_pe181_c +18 F_LMPD_182_pe181_c -18 R_LMPD_182_pe181_c +18 F_LMPD_183_pe181_c -18 R_LMPD_183_pe181_c +18 F_LMPD_184_pe181_c -18 R_LMPD_184_pe181_c +18 F_LMPD_185_pe181_c -18 R_LMPD_185_pe181_c +18 F_LMPD_186_pe181_c -18 R_LMPD_186_pe181_c +18 F_LMPD_187_pe181_c -18 R_LMPD_187_pe181_c +2 F_LMPD_188_pe181_c -2 R_LMPD_188_pe181_c +2 F_LMPD_189_pe181_c -2 R_LMPD_189_pe181_c +16 F_LMPD_190_pg160_c -16 R_LMPD_190_pg160_c +16 F_LMPD_191_pg160_c -16 R_LMPD_191_pg160_c +16 F_LMPD_192_pg160_c -16 R_LMPD_192_pg160_c +16 F_LMPD_193_pg160_c -16 R_LMPD_193_pg160_c +16 F_LMPD_194_pg160_c -16 R_LMPD_194_pg160_c +16 F_LMPD_195_pg160_c -16 R_LMPD_195_pg160_c +16 F_LMPD_196_pg160_c -16 R_LMPD_196_pg160_c +16 F_LMPD_197_pg160_c -16 R_LMPD_197_pg160_c +2 F_LMPD_198_pg160_c -2 R_LMPD_198_pg160_c +2 F_LMPD_199_pg160_c -2 R_LMPD_199_pg160_c +16 F_LMPD_200_pg161_c -16 R_LMPD_200_pg161_c +16 F_LMPD_201_pg161_c -16 R_LMPD_201_pg161_c +16 F_LMPD_202_pg161_c -16 R_LMPD_202_pg161_c +16 F_LMPD_203_pg161_c -16 R_LMPD_203_pg161_c +16 F_LMPD_204_pg161_c -16 R_LMPD_204_pg161_c +16 F_LMPD_205_pg161_c -16 R_LMPD_205_pg161_c +16 F_LMPD_206_pg161_c -16 R_LMPD_206_pg161_c +18 F_LMPD_207_pg181_c -18 R_LMPD_207_pg181_c +18 F_LMPD_208_pg181_c -18 R_LMPD_208_pg181_c +18 F_LMPD_209_pg181_c -18 R_LMPD_209_pg181_c +18 F_LMPD_210_pg181_c -18 R_LMPD_210_pg181_c +18 F_LMPD_211_pg181_c -18 R_LMPD_211_pg181_c +18 F_LMPD_212_pg181_c -18 R_LMPD_212_pg181_c +18 F_LMPD_213_pg181_c -18 R_LMPD_213_pg181_c +18 F_LMPD_214_pg181_c -18 R_LMPD_214_pg181_c +2 F_LMPD_215_pg181_c -2 R_LMPD_215_pg181_c +2 F_LMPD_216_pg181_c -2 R_LMPD_216_pg181_c + F_LMPD_217_phe-L_c - R_LMPD_217_phe-L_c + F_LMPD_220_pro-L_c - R_LMPD_220_pro-L_c + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c +25 F_LMPD_224_q8h2_c -25 R_LMPD_224_q8h2_c +27 F_LMPD_225_q8h2_c -27 R_LMPD_225_q8h2_c +31 F_LMPD_226_q8h2_c -31 R_LMPD_226_q8h2_c +28 F_LMPD_227_q8h2_c -28 R_LMPD_227_q8h2_c +33 F_LMPD_228_q8h2_c -33 R_LMPD_228_q8h2_c +30 F_LMPD_229_q8h2_c -30 R_LMPD_229_q8h2_c +37 F_LMPD_230_q8h2_c -37 R_LMPD_230_q8h2_c +34 F_LMPD_231_q8h2_c -34 R_LMPD_231_q8h2_c +58 F_LMPD_232_q8h2_c -58 R_LMPD_232_q8h2_c +70 F_LMPD_233_q8h2_c -70 R_LMPD_233_q8h2_c +64 F_LMPD_234_q8h2_c -64 R_LMPD_234_q8h2_c +8 F_LMPD_235_ribflv_c -8 R_LMPD_235_ribflv_c +10 F_LMPD_236_ribflv_c -10 R_LMPD_236_ribflv_c +9 F_LMPD_237_ribflv_c -9 R_LMPD_237_ribflv_c +11 F_LMPD_238_ribflv_c -11 R_LMPD_238_ribflv_c +2 F_LMPD_240_sheme_c -2 R_LMPD_240_sheme_c +2 F_LMPD_241_spmd_c -2 R_LMPD_241_spmd_c +4 F_LMPD_242_spmd_c -4 R_LMPD_242_spmd_c +11 F_LMPD_243_thf_c -11 R_LMPD_243_thf_c +13 F_LMPD_244_thf_c -13 R_LMPD_244_thf_c +12 F_LMPD_245_thf_c -12 R_LMPD_245_thf_c +14 F_LMPD_246_thf_c -14 R_LMPD_246_thf_c +8 F_LMPD_247_thmpp_c -8 R_LMPD_247_thmpp_c +10 F_LMPD_248_thmpp_c -10 R_LMPD_248_thmpp_c +2 F_LMPD_249_thr-L_c -2 R_LMPD_249_thr-L_c +2 F_LMPD_250_trp-L_c -2 R_LMPD_250_trp-L_c + F_LMPD_251_tyr-L_c - R_LMPD_251_tyr-L_c +33 F_LMPD_252_udcpdp_c -33 R_LMPD_252_udcpdp_c +3 F_LMPD_253_utp_c -3 R_LMPD_253_utp_c = 0
M_adpglc_c: - F_GLCS1 + R_GLCS1 + F_GLGC - R_GLGC = 0
M_akg_c: - F_AKGDH + R_AKGDH + F_AKGt2rpp - R_AKGt2rpp - F_ASPTA + R_ASPTA + F_GLUDy - R_GLUDy + F_ICDHyr - R_ICDHyr +2 F_LMPD_1_10fthf_c -2 R_LMPD_1_10fthf_c +2 F_LMPD_2_10fthf_c -2 R_LMPD_2_10fthf_c +3 F_LMPD_3_10fthf_c -3 R_LMPD_3_10fthf_c +3 F_LMPD_4_10fthf_c -3 R_LMPD_4_10fthf_c - F_LMPD_5_2dmmql8_c + R_LMPD_5_2dmmql8_c +2 F_LMPD_6_2fe2s_c -2 R_LMPD_6_2fe2s_c +2 F_LMPD_7_2fe2s_c -2 R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c +2 F_LMPD_9_2fe2s_c -2 R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c +2 F_LMPD_11_2fe2s_c -2 R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +4 F_LMPD_14_4fe4s_c -4 R_LMPD_14_4fe4s_c +4 F_LMPD_15_4fe4s_c -4 R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c +4 F_LMPD_17_4fe4s_c -4 R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c +4 F_LMPD_19_4fe4s_c -4 R_LMPD_19_4fe4s_c +4 F_LMPD_20_4fe4s_c -4 R_LMPD_20_4fe4s_c +4 F_LMPD_21_4fe4s_c -4 R_LMPD_21_4fe4s_c + F_LMPD_22_5mthf_c - R_LMPD_22_5mthf_c + F_LMPD_23_5mthf_c - R_LMPD_23_5mthf_c +2 F_LMPD_24_accoa_c -2 R_LMPD_24_accoa_c +2 F_LMPD_25_accoa_c -2 R_LMPD_25_accoa_c +2 F_LMPD_26_adocbl_c -2 R_LMPD_26_adocbl_c +2 F_LMPD_27_adocbl_c -2 R_LMPD_27_adocbl_c + F_LMPD_28_adocbl_c - R_LMPD_28_adocbl_c + F_LMPD_29_adocbl_c - R_LMPD_29_adocbl_c + F_LMPD_30_ala-L_c - R_LMPD_30_ala-L_c +2 F_LMPD_31_amet_c -2 R_LMPD_31_amet_c +2 F_LMPD_32_amet_c -2 R_LMPD_32_amet_c + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c +12 F_LMPD_35_bmocogdp_c -12 R_LMPD_35_bmocogdp_c +12 F_LMPD_36_bmocogdp_c -12 R_LMPD_36_bmocogdp_c +8 F_LMPD_37_bmocogdp_c -8 R_LMPD_37_bmocogdp_c +8 F_LMPD_38_bmocogdp_c -8 R_LMPD_38_bmocogdp_c +4 F_LMPD_39_btn_c -4 R_LMPD_39_btn_c +4 F_LMPD_40_btn_c -4 R_LMPD_40_btn_c +4 F_LMPD_41_btn_c -4 R_LMPD_41_btn_c +4 F_LMPD_42_btn_c -4 R_LMPD_42_btn_c +6 F_LMPD_43_btn_c -6 R_LMPD_43_btn_c +6 F_LMPD_44_btn_c -6 R_LMPD_44_btn_c +6 F_LMPD_45_btn_c -6 R_LMPD_45_btn_c +6 F_LMPD_47_btn_c -6 R_LMPD_47_btn_c + F_LMPD_48_ca2_c - R_LMPD_48_ca2_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c + F_LMPD_92_cys-L_c - R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c +2 F_LMPD_94_datp_c -2 R_LMPD_94_datp_c +2 F_LMPD_95_datp_c -2 R_LMPD_95_datp_c + F_LMPD_96_datp_c - R_LMPD_96_datp_c + F_LMPD_97_datp_c - R_LMPD_97_datp_c +2 F_LMPD_99_dgtp_c -2 R_LMPD_99_dgtp_c +2 F_LMPD_100_dgtp_c -2 R_LMPD_100_dgtp_c +2 F_LMPD_101_dgtp_c -2 R_LMPD_101_dgtp_c + F_LMPD_102_dgtp_c - R_LMPD_102_dgtp_c + F_LMPD_103_dgtp_c - R_LMPD_103_dgtp_c + F_LMPD_104_dgtp_c - R_LMPD_104_dgtp_c +3 F_LMPD_107_enter_c -3 R_LMPD_107_enter_c +2 F_LMPD_108_fad_c -2 R_LMPD_108_fad_c +2 F_LMPD_109_fad_c -2 R_LMPD_109_fad_c +4 F_LMPD_110_fad_c -4 R_LMPD_110_fad_c +4 F_LMPD_111_fad_c -4 R_LMPD_111_fad_c +2 F_LMPD_112_GAM_c -2 R_LMPD_112_GAM_c +2 F_LMPD_113_GAM_c -2 R_LMPD_113_GAM_c + F_LMPD_114_GAM_c - R_LMPD_114_GAM_c + F_LMPD_115_GAM_c - R_LMPD_115_GAM_c +2 F_LMPD_118_gthrd_c -2 R_LMPD_118_gthrd_c +2 F_LMPD_119_gthrd_c -2 R_LMPD_119_gthrd_c +2 F_LMPD_120_gtp_c -2 R_LMPD_120_gtp_c +2 F_LMPD_121_gtp_c -2 R_LMPD_121_gtp_c + F_LMPD_122_gtp_c - R_LMPD_122_gtp_c + F_LMPD_123_gtp_c - R_LMPD_123_gtp_c + F_LMPD_126_his-L_c - R_LMPD_126_his-L_c + F_LMPD_127_his-L_c - R_LMPD_127_his-L_c + F_LMPD_128_ile-L_c - R_LMPD_128_ile-L_c + F_LMPD_129_leu-L_c - R_LMPD_129_leu-L_c +4 F_LMPD_130_lipopb_c -4 R_LMPD_130_lipopb_c +4 F_LMPD_131_lipopb_c -4 R_LMPD_131_lipopb_c +4 F_LMPD_132_lipopb_c -4 R_LMPD_132_lipopb_c +4 F_LMPD_133_lipopb_c -4 R_LMPD_133_lipopb_c + F_LMPD_134_lys-L_c - R_LMPD_134_lys-L_c +2 F_LMPD_135_malcoa_c -2 R_LMPD_135_malcoa_c +2 F_LMPD_136_malcoa_c -2 R_LMPD_136_malcoa_c +2 F_LMPD_137_met-L_c -2 R_LMPD_137_met-L_c +2 F_LMPD_138_met-L_c -2 R_LMPD_138_met-L_c + F_LMPD_139_mlthf_c - R_LMPD_139_mlthf_c + F_LMPD_140_mlthf_c - R_LMPD_140_mlthf_c +3 F_LMPD_141_mococdp_c -3 R_LMPD_141_mococdp_c +3 F_LMPD_142_mococdp_c -3 R_LMPD_142_mococdp_c +3 F_LMPD_143_mococdp_c -3 R_LMPD_143_mococdp_c +3 F_LMPD_144_mococdp_c -3 R_LMPD_144_mococdp_c +4 F_LMPD_145_mococdp_c -4 R_LMPD_145_mococdp_c +4 F_LMPD_146_mococdp_c -4 R_LMPD_146_mococdp_c +4 F_LMPD_147_mococdp_c -4 R_LMPD_147_mococdp_c +4 F_LMPD_148_mococdp_c -4 R_LMPD_148_mococdp_c +4 F_LMPD_149_mocogdp_c -4 R_LMPD_149_mocogdp_c +4 F_LMPD_150_mocogdp_c -4 R_LMPD_150_mocogdp_c +6 F_LMPD_151_mocogdp_c -6 R_LMPD_151_mocogdp_c +6 F_LMPD_152_mocogdp_c -6 R_LMPD_152_mocogdp_c - F_LMPD_153_mql8_c + R_LMPD_153_mql8_c +4 F_LMPD_154_murein3p3p_p -4 R_LMPD_154_murein3p3p_p +8 F_LMPD_155_murein3px4p_p -8 R_LMPD_155_murein3px4p_p +8 F_LMPD_156_murein4p4p_p -8 R_LMPD_156_murein4p4p_p +8 F_LMPD_157_murein4px4p_p -8 R_LMPD_157_murein4px4p_p +12 F_LMPD_158_murein4px4px4p_p -12 R_LMPD_158_murein4px4px4p_p +2 F_LMPD_159_nad_c -2 R_LMPD_159_nad_c +2 F_LMPD_160_nad_c -2 R_LMPD_160_nad_c + F_LMPD_161_nad_c - R_LMPD_161_nad_c + F_LMPD_162_nad_c - R_LMPD_162_nad_c + F_LMPD_163_pe160_c - R_LMPD_163_pe160_c + F_LMPD_164_pe160_c - R_LMPD_164_pe160_c + F_LMPD_165_pe160_c - R_LMPD_165_pe160_c + F_LMPD_166_pe160_c - R_LMPD_166_pe160_c + F_LMPD_167_pe160_c - R_LMPD_167_pe160_c + F_LMPD_168_pe160_c - R_LMPD_168_pe160_c + F_LMPD_169_pe160_c - R_LMPD_169_pe160_c + F_LMPD_170_pe160_c - R_LMPD_170_pe160_c + F_LMPD_171_pe160_c - R_LMPD_171_pe160_c + F_LMPD_172_pe160_c - R_LMPD_172_pe160_c + F_LMPD_173_pe161_c - R_LMPD_173_pe161_c + F_LMPD_174_pe161_c - R_LMPD_174_pe161_c + F_LMPD_175_pe161_c - R_LMPD_175_pe161_c + F_LMPD_176_pe161_c - R_LMPD_176_pe161_c + F_LMPD_177_pe161_c - R_LMPD_177_pe161_c + F_LMPD_178_pe161_c - R_LMPD_178_pe161_c + F_LMPD_179_pe161_c - R_LMPD_179_pe161_c + F_LMPD_180_pe181_c - R_LMPD_180_pe181_c + F_LMPD_181_pe181_c - R_LMPD_181_pe181_c + F_LMPD_182_pe181_c - R_LMPD_182_pe181_c + F_LMPD_183_pe181_c - R_LMPD_183_pe181_c + F_LMPD_184_pe181_c - R_LMPD_184_pe181_c + F_LMPD_185_pe181_c - R_LMPD_185_pe181_c + F_LMPD_186_pe181_c - R_LMPD_186_pe181_c + F_LMPD_187_pe181_c - R_LMPD_187_pe181_c + F_LMPD_188_pe181_c - R_LMPD_188_pe181_c + F_LMPD_189_pe181_c - R_LMPD_189_pe181_c + F_LMPD_217_phe-L_c - R_LMPD_217_phe-L_c + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c + F_LMPD_222_pydx5p_c - R_LMPD_222_pydx5p_c + F_LMPD_223_pydx5p_c - R_LMPD_223_pydx5p_c +3 F_LMPD_224_q8h2_c -3 R_LMPD_224_q8h2_c +3 F_LMPD_225_q8h2_c -3 R_LMPD_225_q8h2_c +3 F_LMPD_226_q8h2_c -3 R_LMPD_226_q8h2_c +6 F_LMPD_232_q8h2_c -6 R_LMPD_232_q8h2_c +2 F_LMPD_235_ribflv_c -2 R_LMPD_235_ribflv_c +2 F_LMPD_236_ribflv_c -2 R_LMPD_236_ribflv_c + F_LMPD_237_ribflv_c - R_LMPD_237_ribflv_c + F_LMPD_238_ribflv_c - R_LMPD_238_ribflv_c + F_LMPD_239_ser-L_c - R_LMPD_239_ser-L_c +3 F_LMPD_241_spmd_c -3 R_LMPD_241_spmd_c +3 F_LMPD_242_spmd_c -3 R_LMPD_242_spmd_c +2 F_LMPD_243_thf_c -2 R_LMPD_243_thf_c +2 F_LMPD_244_thf_c -2 R_LMPD_244_thf_c + F_LMPD_245_thf_c - R_LMPD_245_thf_c + F_LMPD_246_thf_c - R_LMPD_246_thf_c +3 F_LMPD_247_thmpp_c -3 R_LMPD_247_thmpp_c +3 F_LMPD_248_thmpp_c -3 R_LMPD_248_thmpp_c + F_LMPD_251_tyr-L_c - R_LMPD_251_tyr-L_c + F_LMPD_254_val-L_c - R_LMPD_254_val-L_c = 0
M_ala-L_c: -0.49915 F_Ec_biomass_iJO1366_WT_53p95M +0.49915 R_Ec_biomass_iJO1366_WT_53p95M + F_ALAt2rpp - R_ALAt2rpp +2 F_LMPD_6_2fe2s_c -2 R_LMPD_6_2fe2s_c +2 F_LMPD_7_2fe2s_c -2 R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c +2 F_LMPD_9_2fe2s_c -2 R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c +2 F_LMPD_11_2fe2s_c -2 R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +4 F_LMPD_14_4fe4s_c -4 R_LMPD_14_4fe4s_c +4 F_LMPD_15_4fe4s_c -4 R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c +4 F_LMPD_17_4fe4s_c -4 R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c +4 F_LMPD_19_4fe4s_c -4 R_LMPD_19_4fe4s_c +4 F_LMPD_20_4fe4s_c -4 R_LMPD_20_4fe4s_c +4 F_LMPD_21_4fe4s_c -4 R_LMPD_21_4fe4s_c + F_LMPD_30_ala-L_c - R_LMPD_30_ala-L_c +4 F_LMPD_35_bmocogdp_c -4 R_LMPD_35_bmocogdp_c +4 F_LMPD_36_bmocogdp_c -4 R_LMPD_36_bmocogdp_c +4 F_LMPD_37_bmocogdp_c -4 R_LMPD_37_bmocogdp_c +4 F_LMPD_38_bmocogdp_c -4 R_LMPD_38_bmocogdp_c +4 F_LMPD_130_lipopb_c -4 R_LMPD_130_lipopb_c +4 F_LMPD_131_lipopb_c -4 R_LMPD_131_lipopb_c +4 F_LMPD_132_lipopb_c -4 R_LMPD_132_lipopb_c +4 F_LMPD_133_lipopb_c -4 R_LMPD_133_lipopb_c +2 F_LMPD_141_mococdp_c -2 R_LMPD_141_mococdp_c +2 F_LMPD_142_mococdp_c -2 R_LMPD_142_mococdp_c +2 F_LMPD_143_mococdp_c -2 R_LMPD_143_mococdp_c +2 F_LMPD_144_mococdp_c -2 R_LMPD_144_mococdp_c +2 F_LMPD_145_mococdp_c -2 R_LMPD_145_mococdp_c +2 F_LMPD_146_mococdp_c -2 R_LMPD_146_mococdp_c +2 F_LMPD_147_mococdp_c -2 R_LMPD_147_mococdp_c +2 F_LMPD_148_mococdp_c -2 R_LMPD_148_mococdp_c +2 F_LMPD_149_mocogdp_c -2 R_LMPD_149_mocogdp_c +2 F_LMPD_150_mocogdp_c -2 R_LMPD_150_mocogdp_c +2 F_LMPD_151_mocogdp_c -2 R_LMPD_151_mocogdp_c +2 F_LMPD_152_mocogdp_c -2 R_LMPD_152_mocogdp_c + F_LMPD_247_thmpp_c - R_LMPD_247_thmpp_c + F_LMPD_248_thmpp_c - R_LMPD_248_thmpp_c = 0
M_amet_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_31_amet_c - R_LMPD_31_amet_c + F_LMPD_32_amet_c - R_LMPD_32_amet_c = 0
M_amob_c: - F_DM_AMOB + R_DM_AMOB + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c + F_LMPD_43_btn_c - R_LMPD_43_btn_c + F_LMPD_44_btn_c - R_LMPD_44_btn_c + F_LMPD_45_btn_c - R_LMPD_45_btn_c + F_LMPD_47_btn_c - R_LMPD_47_btn_c = 0
M_amp_c: F_ACS - R_ACS - F_ADK1 + R_ADK1 + F_NTPP6 - R_NTPP6 + F_PPS - R_PPS +4 F_LMPD_1_10fthf_c -4 R_LMPD_1_10fthf_c +3 F_LMPD_2_10fthf_c -3 R_LMPD_2_10fthf_c +4 F_LMPD_3_10fthf_c -4 R_LMPD_3_10fthf_c +3 F_LMPD_4_10fthf_c -3 R_LMPD_4_10fthf_c + F_LMPD_5_2dmmql8_c - R_LMPD_5_2dmmql8_c +4 F_LMPD_6_2fe2s_c -4 R_LMPD_6_2fe2s_c +4 F_LMPD_7_2fe2s_c -4 R_LMPD_7_2fe2s_c +4 F_LMPD_8_2fe2s_c -4 R_LMPD_8_2fe2s_c +4 F_LMPD_9_2fe2s_c -4 R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c +2 F_LMPD_11_2fe2s_c -2 R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +8 F_LMPD_14_4fe4s_c -8 R_LMPD_14_4fe4s_c +8 F_LMPD_15_4fe4s_c -8 R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c +4 F_LMPD_17_4fe4s_c -4 R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c +4 F_LMPD_19_4fe4s_c -4 R_LMPD_19_4fe4s_c +8 F_LMPD_20_4fe4s_c -8 R_LMPD_20_4fe4s_c +8 F_LMPD_21_4fe4s_c -8 R_LMPD_21_4fe4s_c +4 F_LMPD_22_5mthf_c -4 R_LMPD_22_5mthf_c +3 F_LMPD_23_5mthf_c -3 R_LMPD_23_5mthf_c +6 F_LMPD_24_accoa_c -6 R_LMPD_24_accoa_c +4 F_LMPD_25_accoa_c -4 R_LMPD_25_accoa_c +3 F_LMPD_26_adocbl_c -3 R_LMPD_26_adocbl_c +2 F_LMPD_27_adocbl_c -2 R_LMPD_27_adocbl_c +3 F_LMPD_28_adocbl_c -3 R_LMPD_28_adocbl_c +2 F_LMPD_29_adocbl_c -2 R_LMPD_29_adocbl_c +5 F_LMPD_31_amet_c -5 R_LMPD_31_amet_c +3 F_LMPD_32_amet_c -3 R_LMPD_32_amet_c + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c + F_LMPD_34_asn-L_c - R_LMPD_34_asn-L_c +26 F_LMPD_35_bmocogdp_c -26 R_LMPD_35_bmocogdp_c +18 F_LMPD_36_bmocogdp_c -18 R_LMPD_36_bmocogdp_c +26 F_LMPD_37_bmocogdp_c -26 R_LMPD_37_bmocogdp_c +18 F_LMPD_38_bmocogdp_c -18 R_LMPD_38_bmocogdp_c +11 F_LMPD_39_btn_c -11 R_LMPD_39_btn_c +11 F_LMPD_40_btn_c -11 R_LMPD_40_btn_c +8 F_LMPD_41_btn_c -8 R_LMPD_41_btn_c +8 F_LMPD_42_btn_c -8 R_LMPD_42_btn_c +15 F_LMPD_43_btn_c -15 R_LMPD_43_btn_c +15 F_LMPD_44_btn_c -15 R_LMPD_44_btn_c +9 F_LMPD_45_btn_c -9 R_LMPD_45_btn_c +9 F_LMPD_47_btn_c -9 R_LMPD_47_btn_c +4 F_LMPD_60_clpn160_p -4 R_LMPD_60_clpn160_p +4 F_LMPD_61_clpn160_p -4 R_LMPD_61_clpn160_p +4 F_LMPD_77_clpn181_p -4 R_LMPD_77_clpn181_p +4 F_LMPD_78_clpn181_p -4 R_LMPD_78_clpn181_p + F_LMPD_91_ctp_c - R_LMPD_91_ctp_c +2 F_LMPD_92_cys-L_c -2 R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c +2 F_LMPD_94_datp_c -2 R_LMPD_94_datp_c +3 F_LMPD_95_datp_c -3 R_LMPD_95_datp_c +3 F_LMPD_96_datp_c -3 R_LMPD_96_datp_c +2 F_LMPD_97_datp_c -2 R_LMPD_97_datp_c + F_LMPD_98_dctp_c - R_LMPD_98_dctp_c +3 F_LMPD_99_dgtp_c -3 R_LMPD_99_dgtp_c +2 F_LMPD_100_dgtp_c -2 R_LMPD_100_dgtp_c +2 F_LMPD_101_dgtp_c -2 R_LMPD_101_dgtp_c +3 F_LMPD_102_dgtp_c -3 R_LMPD_102_dgtp_c +2 F_LMPD_103_dgtp_c -2 R_LMPD_103_dgtp_c +2 F_LMPD_104_dgtp_c -2 R_LMPD_104_dgtp_c + F_LMPD_105_dttp_c - R_LMPD_105_dttp_c + F_LMPD_106_dttp_c - R_LMPD_106_dttp_c +6 F_LMPD_107_enter_c -6 R_LMPD_107_enter_c +6 F_LMPD_108_fad_c -6 R_LMPD_108_fad_c +4 F_LMPD_109_fad_c -4 R_LMPD_109_fad_c +6 F_LMPD_110_fad_c -6 R_LMPD_110_fad_c +4 F_LMPD_111_fad_c -4 R_LMPD_111_fad_c +3 F_LMPD_112_GAM_c -3 R_LMPD_112_GAM_c +2 F_LMPD_113_GAM_c -2 R_LMPD_113_GAM_c +2 F_LMPD_114_GAM_c -2 R_LMPD_114_GAM_c +3 F_LMPD_115_GAM_c -3 R_LMPD_115_GAM_c +2 F_LMPD_118_gthrd_c -2 R_LMPD_118_gthrd_c + F_LMPD_119_gthrd_c - R_LMPD_119_gthrd_c +3 F_LMPD_120_gtp_c -3 R_LMPD_120_gtp_c +2 F_LMPD_121_gtp_c -2 R_LMPD_121_gtp_c +3 F_LMPD_122_gtp_c -3 R_LMPD_122_gtp_c +2 F_LMPD_123_gtp_c -2 R_LMPD_123_gtp_c +8 F_LMPD_124_hemeO_c -8 R_LMPD_124_hemeO_c +12 F_LMPD_125_hemeO_c -12 R_LMPD_125_hemeO_c +3 F_LMPD_126_his-L_c -3 R_LMPD_126_his-L_c +2 F_LMPD_127_his-L_c -2 R_LMPD_127_his-L_c +13 F_LMPD_130_lipopb_c -13 R_LMPD_130_lipopb_c +9 F_LMPD_131_lipopb_c -9 R_LMPD_131_lipopb_c +9 F_LMPD_132_lipopb_c -9 R_LMPD_132_lipopb_c +9 F_LMPD_133_lipopb_c -9 R_LMPD_133_lipopb_c +6 F_LMPD_135_malcoa_c -6 R_LMPD_135_malcoa_c +4 F_LMPD_136_malcoa_c -4 R_LMPD_136_malcoa_c +2 F_LMPD_137_met-L_c -2 R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c +4 F_LMPD_139_mlthf_c -4 R_LMPD_139_mlthf_c +3 F_LMPD_140_mlthf_c -3 R_LMPD_140_mlthf_c +11 F_LMPD_141_mococdp_c -11 R_LMPD_141_mococdp_c +11 F_LMPD_142_mococdp_c -11 R_LMPD_142_mococdp_c +8 F_LMPD_143_mococdp_c -8 R_LMPD_143_mococdp_c +8 F_LMPD_144_mococdp_c -8 R_LMPD_144_mococdp_c +11 F_LMPD_145_mococdp_c -11 R_LMPD_145_mococdp_c +11 F_LMPD_146_mococdp_c -11 R_LMPD_146_mococdp_c +8 F_LMPD_147_mococdp_c -8 R_LMPD_147_mococdp_c +8 F_LMPD_148_mococdp_c -8 R_LMPD_148_mococdp_c +13 F_LMPD_149_mocogdp_c -13 R_LMPD_149_mocogdp_c +9 F_LMPD_150_mocogdp_c -9 R_LMPD_150_mocogdp_c +13 F_LMPD_151_mocogdp_c -13 R_LMPD_151_mocogdp_c +9 F_LMPD_152_mocogdp_c -9 R_LMPD_152_mocogdp_c +3 F_LMPD_153_mql8_c -3 R_LMPD_153_mql8_c +5 F_LMPD_159_nad_c -5 R_LMPD_159_nad_c +4 F_LMPD_160_nad_c -4 R_LMPD_160_nad_c +5 F_LMPD_161_nad_c -5 R_LMPD_161_nad_c +4 F_LMPD_162_nad_c -4 R_LMPD_162_nad_c +2 F_LMPD_171_pe160_c -2 R_LMPD_171_pe160_c +2 F_LMPD_172_pe160_c -2 R_LMPD_172_pe160_c +2 F_LMPD_188_pe181_c -2 R_LMPD_188_pe181_c +2 F_LMPD_189_pe181_c -2 R_LMPD_189_pe181_c +2 F_LMPD_198_pg160_c -2 R_LMPD_198_pg160_c +2 F_LMPD_199_pg160_c -2 R_LMPD_199_pg160_c +2 F_LMPD_215_pg181_c -2 R_LMPD_215_pg181_c +2 F_LMPD_216_pg181_c -2 R_LMPD_216_pg181_c +8 F_LMPD_218_pheme_c -8 R_LMPD_218_pheme_c +12 F_LMPD_219_pheme_c -12 R_LMPD_219_pheme_c +6 F_LMPD_224_q8h2_c -6 R_LMPD_224_q8h2_c +6 F_LMPD_225_q8h2_c -6 R_LMPD_225_q8h2_c +6 F_LMPD_226_q8h2_c -6 R_LMPD_226_q8h2_c +6 F_LMPD_227_q8h2_c -6 R_LMPD_227_q8h2_c +3 F_LMPD_228_q8h2_c -3 R_LMPD_228_q8h2_c +6 F_LMPD_229_q8h2_c -6 R_LMPD_229_q8h2_c +3 F_LMPD_230_q8h2_c -3 R_LMPD_230_q8h2_c +6 F_LMPD_231_q8h2_c -6 R_LMPD_231_q8h2_c +12 F_LMPD_232_q8h2_c -12 R_LMPD_232_q8h2_c +6 F_LMPD_233_q8h2_c -6 R_LMPD_233_q8h2_c +12 F_LMPD_234_q8h2_c -12 R_LMPD_234_q8h2_c +3 F_LMPD_235_ribflv_c -3 R_LMPD_235_ribflv_c +2 F_LMPD_236_ribflv_c -2 R_LMPD_236_ribflv_c +3 F_LMPD_237_ribflv_c -3 R_LMPD_237_ribflv_c +2 F_LMPD_238_ribflv_c -2 R_LMPD_238_ribflv_c +12 F_LMPD_240_sheme_c -12 R_LMPD_240_sheme_c +4 F_LMPD_241_spmd_c -4 R_LMPD_241_spmd_c +3 F_LMPD_242_spmd_c -3 R_LMPD_242_spmd_c +4 F_LMPD_243_thf_c -4 R_LMPD_243_thf_c +3 F_LMPD_244_thf_c -3 R_LMPD_244_thf_c +4 F_LMPD_245_thf_c -4 R_LMPD_245_thf_c +3 F_LMPD_246_thf_c -3 R_LMPD_246_thf_c +6 F_LMPD_247_thmpp_c -6 R_LMPD_247_thmpp_c +5 F_LMPD_248_thmpp_c -5 R_LMPD_248_thmpp_c + F_LMPD_250_trp-L_c - R_LMPD_250_trp-L_c + F_LMPD_253_utp_c - R_LMPD_253_utp_c = 0
M_arg-L_c: -0.28742 F_Ec_biomass_iJO1366_WT_53p95M +0.28742 R_Ec_biomass_iJO1366_WT_53p95M - F_ARGt3pp + R_ARGt3pp + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c = 0
M_asn-L_c: -0.23423 F_Ec_biomass_iJO1366_WT_53p95M +0.23423 R_Ec_biomass_iJO1366_WT_53p95M + F_ASNt2rpp - R_ASNt2rpp + F_LMPD_34_asn-L_c - R_LMPD_34_asn-L_c = 0
M_asp-L_c: -0.23423 F_Ec_biomass_iJO1366_WT_53p95M +0.23423 R_Ec_biomass_iJO1366_WT_53p95M - F_ASPO3 + R_ASPO3 - F_ASPO4 + R_ASPO4 - F_ASPO5 + R_ASPO5 - F_ASPO6 + R_ASPO6 - F_ASPT + R_ASPT - F_ASPTA + R_ASPTA + F_ASPt2rpp - R_ASPt2rpp - F_LMPD_1_10fthf_c + R_LMPD_1_10fthf_c - F_LMPD_2_10fthf_c + R_LMPD_2_10fthf_c - F_LMPD_3_10fthf_c + R_LMPD_3_10fthf_c - F_LMPD_4_10fthf_c + R_LMPD_4_10fthf_c - F_LMPD_22_5mthf_c + R_LMPD_22_5mthf_c - F_LMPD_23_5mthf_c + R_LMPD_23_5mthf_c -3 F_LMPD_24_accoa_c +3 R_LMPD_24_accoa_c -3 F_LMPD_25_accoa_c +3 R_LMPD_25_accoa_c -2 F_LMPD_26_adocbl_c +2 R_LMPD_26_adocbl_c -2 F_LMPD_27_adocbl_c +2 R_LMPD_27_adocbl_c -2 F_LMPD_28_adocbl_c +2 R_LMPD_28_adocbl_c -2 F_LMPD_29_adocbl_c +2 R_LMPD_29_adocbl_c -3 F_LMPD_31_amet_c +3 R_LMPD_31_amet_c -3 F_LMPD_32_amet_c +3 R_LMPD_32_amet_c - F_LMPD_33_arg-L_c + R_LMPD_33_arg-L_c - F_LMPD_34_asn-L_c + R_LMPD_34_asn-L_c -4 F_LMPD_35_bmocogdp_c +4 R_LMPD_35_bmocogdp_c -4 F_LMPD_36_bmocogdp_c +4 R_LMPD_36_bmocogdp_c -4 F_LMPD_37_bmocogdp_c +4 R_LMPD_37_bmocogdp_c -4 F_LMPD_38_bmocogdp_c +4 R_LMPD_38_bmocogdp_c -3 F_LMPD_39_btn_c +3 R_LMPD_39_btn_c -3 F_LMPD_40_btn_c +3 R_LMPD_40_btn_c -3 F_LMPD_41_btn_c +3 R_LMPD_41_btn_c -3 F_LMPD_42_btn_c +3 R_LMPD_42_btn_c -8 F_LMPD_43_btn_c +8 R_LMPD_43_btn_c -8 F_LMPD_44_btn_c +8 R_LMPD_44_btn_c -8 F_LMPD_45_btn_c +8 R_LMPD_45_btn_c -8 F_LMPD_47_btn_c +8 R_LMPD_47_btn_c - F_LMPD_91_ctp_c + R_LMPD_91_ctp_c -2 F_LMPD_94_datp_c +2 R_LMPD_94_datp_c -2 F_LMPD_95_datp_c +2 R_LMPD_95_datp_c -2 F_LMPD_96_datp_c +2 R_LMPD_96_datp_c -2 F_LMPD_97_datp_c +2 R_LMPD_97_datp_c - F_LMPD_98_dctp_c + R_LMPD_98_dctp_c - F_LMPD_99_dgtp_c + R_LMPD_99_dgtp_c - F_LMPD_100_dgtp_c + R_LMPD_100_dgtp_c - F_LMPD_101_dgtp_c + R_LMPD_101_dgtp_c - F_LMPD_102_dgtp_c + R_LMPD_102_dgtp_c - F_LMPD_103_dgtp_c + R_LMPD_103_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c - F_LMPD_105_dttp_c + R_LMPD_105_dttp_c - F_LMPD_106_dttp_c + R_LMPD_106_dttp_c -3 F_LMPD_108_fad_c +3 R_LMPD_108_fad_c -3 F_LMPD_109_fad_c +3 R_LMPD_109_fad_c -3 F_LMPD_110_fad_c +3 R_LMPD_110_fad_c -3 F_LMPD_111_fad_c +3 R_LMPD_111_fad_c -2 F_LMPD_112_GAM_c +2 R_LMPD_112_GAM_c -2 F_LMPD_113_GAM_c +2 R_LMPD_113_GAM_c -2 F_LMPD_114_GAM_c +2 R_LMPD_114_GAM_c -2 F_LMPD_115_GAM_c +2 R_LMPD_115_GAM_c - F_LMPD_117_gly_c + R_LMPD_117_gly_c - F_LMPD_120_gtp_c + R_LMPD_120_gtp_c - F_LMPD_121_gtp_c + R_LMPD_121_gtp_c - F_LMPD_122_gtp_c + R_LMPD_122_gtp_c - F_LMPD_123_gtp_c + R_LMPD_123_gtp_c - F_LMPD_126_his-L_c + R_LMPD_126_his-L_c - F_LMPD_127_his-L_c + R_LMPD_127_his-L_c - F_LMPD_128_ile-L_c + R_LMPD_128_ile-L_c - F_LMPD_134_lys-L_c + R_LMPD_134_lys-L_c -3 F_LMPD_135_malcoa_c +3 R_LMPD_135_malcoa_c -3 F_LMPD_136_malcoa_c +3 R_LMPD_136_malcoa_c - F_LMPD_137_met-L_c + R_LMPD_137_met-L_c - F_LMPD_138_met-L_c + R_LMPD_138_met-L_c - F_LMPD_139_mlthf_c + R_LMPD_139_mlthf_c - F_LMPD_140_mlthf_c + R_LMPD_140_mlthf_c -2 F_LMPD_141_mococdp_c +2 R_LMPD_141_mococdp_c -2 F_LMPD_142_mococdp_c +2 R_LMPD_142_mococdp_c -2 F_LMPD_143_mococdp_c +2 R_LMPD_143_mococdp_c -2 F_LMPD_144_mococdp_c +2 R_LMPD_144_mococdp_c -2 F_LMPD_145_mococdp_c +2 R_LMPD_145_mococdp_c -2 F_LMPD_146_mococdp_c +2 R_LMPD_146_mococdp_c -2 F_LMPD_147_mococdp_c +2 R_LMPD_147_mococdp_c -2 F_LMPD_148_mococdp_c +2 R_LMPD_148_mococdp_c -2 F_LMPD_149_mocogdp_c +2 R_LMPD_149_mocogdp_c -2 F_LMPD_150_mocogdp_c +2 R_LMPD_150_mocogdp_c -2 F_LMPD_151_mocogdp_c +2 R_LMPD_151_mocogdp_c -2 F_LMPD_152_mocogdp_c +2 R_LMPD_152_mocogdp_c -2 F_LMPD_154_murein3p3p_p +2 R_LMPD_154_murein3p3p_p -2 F_LMPD_155_murein3px4p_p +2 R_LMPD_155_murein3px4p_p -2 F_LMPD_156_murein4p4p_p +2 R_LMPD_156_murein4p4p_p -2 F_LMPD_157_murein4px4p_p +2 R_LMPD_157_murein4px4p_p -3 F_LMPD_158_murein4px4px4p_p +3 R_LMPD_158_murein4px4px4p_p -2 F_LMPD_159_nad_c +2 R_LMPD_159_nad_c -2 F_LMPD_160_nad_c +2 R_LMPD_160_nad_c -2 F_LMPD_161_nad_c +2 R_LMPD_161_nad_c -2 F_LMPD_162_nad_c +2 R_LMPD_162_nad_c - F_LMPD_235_ribflv_c + R_LMPD_235_ribflv_c - F_LMPD_236_ribflv_c + R_LMPD_236_ribflv_c - F_LMPD_237_ribflv_c + R_LMPD_237_ribflv_c - F_LMPD_238_ribflv_c + R_LMPD_238_ribflv_c - F_LMPD_241_spmd_c + R_LMPD_241_spmd_c - F_LMPD_242_spmd_c + R_LMPD_242_spmd_c - F_LMPD_243_thf_c + R_LMPD_243_thf_c - F_LMPD_244_thf_c + R_LMPD_244_thf_c - F_LMPD_245_thf_c + R_LMPD_245_thf_c - F_LMPD_246_thf_c + R_LMPD_246_thf_c - F_LMPD_249_thr-L_c + R_LMPD_249_thr-L_c - F_LMPD_253_utp_c + R_LMPD_253_utp_c = 0
M_atp_c: -54.12 F_Ec_biomass_iJO1366_WT_53p95M +54.12 R_Ec_biomass_iJO1366_WT_53p95M - F_ACKr + R_ACKr - F_ACS + R_ACS - F_ADK1 + R_ADK1 - F_ATPM + R_ATPM + F_ATPS4rpp - R_ATPS4rpp - F_COLIPAabcpp + R_COLIPAabcpp - F_COLIPAabctex + R_COLIPAabctex - F_CYSabc2pp + R_CYSabc2pp - F_FE2abcpp + R_FE2abcpp - F_FE3abcpp + R_FE3abcpp - F_GLCabcpp + R_GLCabcpp - F_GLGC + R_GLGC - F_GLUabcpp + R_GLUabcpp - F_GLYCK + R_GLYCK - F_GLYCK2 + R_GLYCK2 - F_GLYK + R_GLYK - F_GTHRDabc2pp + R_GTHRDabc2pp - F_HEX1 + R_HEX1 - F_Kabcpp + R_Kabcpp - F_MG2uabcpp + R_MG2uabcpp - F_MOBDabcpp + R_MOBDabcpp - F_NADK + R_NADK - F_NI2uabcpp + R_NI2uabcpp - F_NTP1 + R_NTP1 - F_NTPP6 + R_NTPP6 - F_PE160abcpp + R_PE160abcpp - F_PE161abcpp + R_PE161abcpp - F_PE181abcpp + R_PE181abcpp - F_PFK + R_PFK - F_PFK_3 + R_PFK_3 - F_PG160abcpp + R_PG160abcpp - F_PG161abcpp + R_PG161abcpp - F_PG181abcpp + R_PG181abcpp - F_PGK + R_PGK - F_PHEMEabcpp + R_PHEMEabcpp - F_PIuabcpp + R_PIuabcpp - F_PPCK + R_PPCK - F_PPK2r + R_PPK2r - F_PPKr + R_PPKr - F_PPS + R_PPS + F_PYK - R_PYK - F_SUCOAS + R_SUCOAS - F_SULabcpp + R_SULabcpp - F_ZNabcpp + R_ZNabcpp -16 F_LMPD_1_10fthf_c +16 R_LMPD_1_10fthf_c -17 F_LMPD_2_10fthf_c +17 R_LMPD_2_10fthf_c -15 F_LMPD_3_10fthf_c +15 R_LMPD_3_10fthf_c -16 F_LMPD_4_10fthf_c +16 R_LMPD_4_10fthf_c -26 F_LMPD_5_2dmmql8_c +26 R_LMPD_5_2dmmql8_c -4 F_LMPD_6_2fe2s_c +4 R_LMPD_6_2fe2s_c -4 F_LMPD_7_2fe2s_c +4 R_LMPD_7_2fe2s_c -5 F_LMPD_8_2fe2s_c +5 R_LMPD_8_2fe2s_c -5 F_LMPD_9_2fe2s_c +5 R_LMPD_9_2fe2s_c -6 F_LMPD_10_2fe2s_c +6 R_LMPD_10_2fe2s_c -6 F_LMPD_11_2fe2s_c +6 R_LMPD_11_2fe2s_c -7 F_LMPD_12_2fe2s_c +7 R_LMPD_12_2fe2s_c -7 F_LMPD_13_2fe2s_c +7 R_LMPD_13_2fe2s_c -10 F_LMPD_14_4fe4s_c +10 R_LMPD_14_4fe4s_c -10 F_LMPD_15_4fe4s_c +10 R_LMPD_15_4fe4s_c -12 F_LMPD_16_4fe4s_c +12 R_LMPD_16_4fe4s_c -12 F_LMPD_17_4fe4s_c +12 R_LMPD_17_4fe4s_c -14 F_LMPD_18_4fe4s_c +14 R_LMPD_18_4fe4s_c -14 F_LMPD_19_4fe4s_c +14 R_LMPD_19_4fe4s_c -8 F_LMPD_20_4fe4s_c +8 R_LMPD_20_4fe4s_c -8 F_LMPD_21_4fe4s_c +8 R_LMPD_21_4fe4s_c -17 F_LMPD_22_5mthf_c +17 R_LMPD_22_5mthf_c -18 F_LMPD_23_5mthf_c +18 R_LMPD_23_5mthf_c -18 F_LMPD_24_accoa_c +18 R_LMPD_24_accoa_c -20 F_LMPD_25_accoa_c +20 R_LMPD_25_accoa_c -10 F_LMPD_26_adocbl_c +10 R_LMPD_26_adocbl_c -11 F_LMPD_27_adocbl_c +11 R_LMPD_27_adocbl_c -11 F_LMPD_28_adocbl_c +11 R_LMPD_28_adocbl_c -12 F_LMPD_29_adocbl_c +12 R_LMPD_29_adocbl_c -14 F_LMPD_31_amet_c +14 R_LMPD_31_amet_c -16 F_LMPD_32_amet_c +16 R_LMPD_32_amet_c -3 F_LMPD_33_arg-L_c +3 R_LMPD_33_arg-L_c - F_LMPD_34_asn-L_c + R_LMPD_34_asn-L_c -58 F_LMPD_35_bmocogdp_c +58 R_LMPD_35_bmocogdp_c -66 F_LMPD_36_bmocogdp_c +66 R_LMPD_36_bmocogdp_c -62 F_LMPD_37_bmocogdp_c +62 R_LMPD_37_bmocogdp_c -70 F_LMPD_38_bmocogdp_c +70 R_LMPD_38_bmocogdp_c -24 F_LMPD_39_btn_c +24 R_LMPD_39_btn_c -25 F_LMPD_40_btn_c +25 R_LMPD_40_btn_c -27 F_LMPD_41_btn_c +27 R_LMPD_41_btn_c -28 F_LMPD_42_btn_c +28 R_LMPD_42_btn_c -45 F_LMPD_43_btn_c +45 R_LMPD_43_btn_c -46 F_LMPD_44_btn_c +46 R_LMPD_44_btn_c -51 F_LMPD_45_btn_c +51 R_LMPD_45_btn_c -52 F_LMPD_47_btn_c +52 R_LMPD_47_btn_c -3 F_LMPD_48_ca2_c +3 R_LMPD_48_ca2_c -3 F_LMPD_49_ca2_c +3 R_LMPD_49_ca2_c -3 F_LMPD_50_ca2_c +3 R_LMPD_50_ca2_c - F_LMPD_51_chor_c + R_LMPD_51_chor_c -32 F_LMPD_52_clpn160_p +32 R_LMPD_52_clpn160_p -32 F_LMPD_53_clpn160_p +32 R_LMPD_53_clpn160_p -32 F_LMPD_54_clpn160_p +32 R_LMPD_54_clpn160_p -32 F_LMPD_55_clpn160_p +32 R_LMPD_55_clpn160_p -32 F_LMPD_56_clpn160_p +32 R_LMPD_56_clpn160_p -32 F_LMPD_57_clpn160_p +32 R_LMPD_57_clpn160_p -32 F_LMPD_58_clpn160_p +32 R_LMPD_58_clpn160_p -32 F_LMPD_59_clpn160_p +32 R_LMPD_59_clpn160_p -8 F_LMPD_60_clpn160_p +8 R_LMPD_60_clpn160_p -8 F_LMPD_61_clpn160_p +8 R_LMPD_61_clpn160_p -32 F_LMPD_62_clpn161_p +32 R_LMPD_62_clpn161_p -32 F_LMPD_63_clpn161_p +32 R_LMPD_63_clpn161_p -32 F_LMPD_64_clpn161_p +32 R_LMPD_64_clpn161_p -32 F_LMPD_65_clpn161_p +32 R_LMPD_65_clpn161_p -32 F_LMPD_66_clpn161_p +32 R_LMPD_66_clpn161_p -32 F_LMPD_67_clpn161_p +32 R_LMPD_67_clpn161_p -32 F_LMPD_68_clpn161_p +32 R_LMPD_68_clpn161_p -36 F_LMPD_69_clpn181_p +36 R_LMPD_69_clpn181_p -36 F_LMPD_70_clpn181_p +36 R_LMPD_70_clpn181_p -36 F_LMPD_71_clpn181_p +36 R_LMPD_71_clpn181_p -36 F_LMPD_72_clpn181_p +36 R_LMPD_72_clpn181_p -36 F_LMPD_73_clpn181_p +36 R_LMPD_73_clpn181_p -36 F_LMPD_74_clpn181_p +36 R_LMPD_74_clpn181_p -36 F_LMPD_75_clpn181_p +36 R_LMPD_75_clpn181_p -36 F_LMPD_76_clpn181_p +36 R_LMPD_76_clpn181_p -8 F_LMPD_77_clpn181_p +8 R_LMPD_77_clpn181_p -8 F_LMPD_78_clpn181_p +8 R_LMPD_78_clpn181_p -62 F_LMPD_79_colipa_e +62 R_LMPD_79_colipa_e -62 F_LMPD_80_colipa_e +62 R_LMPD_80_colipa_e -62 F_LMPD_81_colipa_e +62 R_LMPD_81_colipa_e -62 F_LMPD_82_colipa_e +62 R_LMPD_82_colipa_e -62 F_LMPD_83_colipa_e +62 R_LMPD_83_colipa_e -62 F_LMPD_84_colipa_e +62 R_LMPD_84_colipa_e -62 F_LMPD_85_colipa_e +62 R_LMPD_85_colipa_e -62 F_LMPD_86_colipa_e +62 R_LMPD_86_colipa_e -62 F_LMPD_87_colipa_e +62 R_LMPD_87_colipa_e -62 F_LMPD_88_colipa_e +62 R_LMPD_88_colipa_e -62 F_LMPD_89_colipa_e +62 R_LMPD_89_colipa_e -62 F_LMPD_90_colipa_e +62 R_LMPD_90_colipa_e -6 F_LMPD_91_ctp_c +6 R_LMPD_91_ctp_c -2 F_LMPD_92_cys-L_c +2 R_LMPD_92_cys-L_c -3 F_LMPD_93_cys-L_c +3 R_LMPD_93_cys-L_c -10 F_LMPD_94_datp_c +10 R_LMPD_94_datp_c -9 F_LMPD_95_datp_c +9 R_LMPD_95_datp_c -10 F_LMPD_96_datp_c +10 R_LMPD_96_datp_c -11 F_LMPD_97_datp_c +11 R_LMPD_97_datp_c -6 F_LMPD_98_dctp_c +6 R_LMPD_98_dctp_c -11 F_LMPD_99_dgtp_c +11 R_LMPD_99_dgtp_c -12 F_LMPD_100_dgtp_c +12 R_LMPD_100_dgtp_c -12 F_LMPD_101_dgtp_c +12 R_LMPD_101_dgtp_c -12 F_LMPD_102_dgtp_c +12 R_LMPD_102_dgtp_c -13 F_LMPD_103_dgtp_c +13 R_LMPD_103_dgtp_c -13 F_LMPD_104_dgtp_c +13 R_LMPD_104_dgtp_c -7 F_LMPD_105_dttp_c +7 R_LMPD_105_dttp_c -7 F_LMPD_106_dttp_c +7 R_LMPD_106_dttp_c -9 F_LMPD_107_enter_c +9 R_LMPD_107_enter_c -23 F_LMPD_108_fad_c +23 R_LMPD_108_fad_c -25 F_LMPD_109_fad_c +25 R_LMPD_109_fad_c -21 F_LMPD_110_fad_c +21 R_LMPD_110_fad_c -23 F_LMPD_111_fad_c +23 R_LMPD_111_fad_c -8 F_LMPD_112_GAM_c +8 R_LMPD_112_GAM_c -9 F_LMPD_113_GAM_c +9 R_LMPD_113_GAM_c -10 F_LMPD_114_GAM_c +10 R_LMPD_114_GAM_c -9 F_LMPD_115_GAM_c +9 R_LMPD_115_GAM_c - F_LMPD_116_gln-L_c + R_LMPD_116_gln-L_c -2 F_LMPD_117_gly_c +2 R_LMPD_117_gly_c -4 F_LMPD_118_gthrd_c +4 R_LMPD_118_gthrd_c -5 F_LMPD_119_gthrd_c +5 R_LMPD_119_gthrd_c -11 F_LMPD_120_gtp_c +11 R_LMPD_120_gtp_c -12 F_LMPD_121_gtp_c +12 R_LMPD_121_gtp_c -12 F_LMPD_122_gtp_c +12 R_LMPD_122_gtp_c -13 F_LMPD_123_gtp_c +13 R_LMPD_123_gtp_c -17 F_LMPD_124_hemeO_c +17 R_LMPD_124_hemeO_c -21 F_LMPD_125_hemeO_c +21 R_LMPD_125_hemeO_c -4 F_LMPD_126_his-L_c +4 R_LMPD_126_his-L_c -5 F_LMPD_127_his-L_c +5 R_LMPD_127_his-L_c -2 F_LMPD_128_ile-L_c +2 R_LMPD_128_ile-L_c -13 F_LMPD_130_lipopb_c +13 R_LMPD_130_lipopb_c -17 F_LMPD_131_lipopb_c +17 R_LMPD_131_lipopb_c -19 F_LMPD_132_lipopb_c +19 R_LMPD_132_lipopb_c -17 F_LMPD_133_lipopb_c +17 R_LMPD_133_lipopb_c - F_LMPD_134_lys-L_c + R_LMPD_134_lys-L_c -19 F_LMPD_135_malcoa_c +19 R_LMPD_135_malcoa_c -21 F_LMPD_136_malcoa_c +21 R_LMPD_136_malcoa_c -3 F_LMPD_137_met-L_c +3 R_LMPD_137_met-L_c -4 F_LMPD_138_met-L_c +4 R_LMPD_138_met-L_c -17 F_LMPD_139_mlthf_c +17 R_LMPD_139_mlthf_c -18 F_LMPD_140_mlthf_c +18 R_LMPD_140_mlthf_c -25 F_LMPD_141_mococdp_c +25 R_LMPD_141_mococdp_c -27 F_LMPD_142_mococdp_c +27 R_LMPD_142_mococdp_c -28 F_LMPD_143_mococdp_c +28 R_LMPD_143_mococdp_c -30 F_LMPD_144_mococdp_c +30 R_LMPD_144_mococdp_c -24 F_LMPD_145_mococdp_c +24 R_LMPD_145_mococdp_c -26 F_LMPD_146_mococdp_c +26 R_LMPD_146_mococdp_c -27 F_LMPD_147_mococdp_c +27 R_LMPD_147_mococdp_c -29 F_LMPD_148_mococdp_c +29 R_LMPD_148_mococdp_c -31 F_LMPD_149_mocogdp_c +31 R_LMPD_149_mocogdp_c -35 F_LMPD_150_mocogdp_c +35 R_LMPD_150_mocogdp_c -29 F_LMPD_151_mocogdp_c +29 R_LMPD_151_mocogdp_c -33 F_LMPD_152_mocogdp_c +33 R_LMPD_152_mocogdp_c -29 F_LMPD_153_mql8_c +29 R_LMPD_153_mql8_c -22 F_LMPD_154_murein3p3p_p +22 R_LMPD_154_murein3p3p_p -22 F_LMPD_155_murein3px4p_p +22 R_LMPD_155_murein3px4p_p -22 F_LMPD_156_murein4p4p_p +22 R_LMPD_156_murein4p4p_p -22 F_LMPD_157_murein4px4p_p +22 R_LMPD_157_murein4px4p_p -33 F_LMPD_158_murein4px4px4p_p +33 R_LMPD_158_murein4px4px4p_p -11 F_LMPD_159_nad_c +11 R_LMPD_159_nad_c -12 F_LMPD_160_nad_c +12 R_LMPD_160_nad_c -12 F_LMPD_161_nad_c +12 R_LMPD_161_nad_c -13 F_LMPD_162_nad_c +13 R_LMPD_162_nad_c -16 F_LMPD_163_pe160_c +16 R_LMPD_163_pe160_c -16 F_LMPD_164_pe160_c +16 R_LMPD_164_pe160_c -16 F_LMPD_165_pe160_c +16 R_LMPD_165_pe160_c -16 F_LMPD_166_pe160_c +16 R_LMPD_166_pe160_c -16 F_LMPD_167_pe160_c +16 R_LMPD_167_pe160_c -16 F_LMPD_168_pe160_c +16 R_LMPD_168_pe160_c -16 F_LMPD_169_pe160_c +16 R_LMPD_169_pe160_c -16 F_LMPD_170_pe160_c +16 R_LMPD_170_pe160_c -4 F_LMPD_171_pe160_c +4 R_LMPD_171_pe160_c -4 F_LMPD_172_pe160_c +4 R_LMPD_172_pe160_c -16 F_LMPD_173_pe161_c +16 R_LMPD_173_pe161_c -16 F_LMPD_174_pe161_c +16 R_LMPD_174_pe161_c -16 F_LMPD_175_pe161_c +16 R_LMPD_175_pe161_c -16 F_LMPD_176_pe161_c +16 R_LMPD_176_pe161_c -16 F_LMPD_177_pe161_c +16 R_LMPD_177_pe161_c -16 F_LMPD_178_pe161_c +16 R_LMPD_178_pe161_c -16 F_LMPD_179_pe161_c +16 R_LMPD_179_pe161_c -18 F_LMPD_180_pe181_c +18 R_LMPD_180_pe181_c -18 F_LMPD_181_pe181_c +18 R_LMPD_181_pe181_c -18 F_LMPD_182_pe181_c +18 R_LMPD_182_pe181_c -18 F_LMPD_183_pe181_c +18 R_LMPD_183_pe181_c -18 F_LMPD_184_pe181_c +18 R_LMPD_184_pe181_c -18 F_LMPD_185_pe181_c +18 R_LMPD_185_pe181_c -18 F_LMPD_186_pe181_c +18 R_LMPD_186_pe181_c -18 F_LMPD_187_pe181_c +18 R_LMPD_187_pe181_c -4 F_LMPD_188_pe181_c +4 R_LMPD_188_pe181_c -4 F_LMPD_189_pe181_c +4 R_LMPD_189_pe181_c -16 F_LMPD_190_pg160_c +16 R_LMPD_190_pg160_c -16 F_LMPD_191_pg160_c +16 R_LMPD_191_pg160_c -16 F_LMPD_192_pg160_c +16 R_LMPD_192_pg160_c -16 F_LMPD_193_pg160_c +16 R_LMPD_193_pg160_c -16 F_LMPD_194_pg160_c +16 R_LMPD_194_pg160_c -16 F_LMPD_195_pg160_c +16 R_LMPD_195_pg160_c -16 F_LMPD_196_pg160_c +16 R_LMPD_196_pg160_c -16 F_LMPD_197_pg160_c +16 R_LMPD_197_pg160_c -4 F_LMPD_198_pg160_c +4 R_LMPD_198_pg160_c -4 F_LMPD_199_pg160_c +4 R_LMPD_199_pg160_c -16 F_LMPD_200_pg161_c +16 R_LMPD_200_pg161_c -16 F_LMPD_201_pg161_c +16 R_LMPD_201_pg161_c -16 F_LMPD_202_pg161_c +16 R_LMPD_202_pg161_c -16 F_LMPD_203_pg161_c +16 R_LMPD_203_pg161_c -16 F_LMPD_204_pg161_c +16 R_LMPD_204_pg161_c -16 F_LMPD_205_pg161_c +16 R_LMPD_205_pg161_c -16 F_LMPD_206_pg161_c +16 R_LMPD_206_pg161_c -18 F_LMPD_207_pg181_c +18 R_LMPD_207_pg181_c -18 F_LMPD_208_pg181_c +18 R_LMPD_208_pg181_c -18 F_LMPD_209_pg181_c +18 R_LMPD_209_pg181_c -18 F_LMPD_210_pg181_c +18 R_LMPD_210_pg181_c -18 F_LMPD_211_pg181_c +18 R_LMPD_211_pg181_c -18 F_LMPD_212_pg181_c +18 R_LMPD_212_pg181_c -18 F_LMPD_213_pg181_c +18 R_LMPD_213_pg181_c -18 F_LMPD_214_pg181_c +18 R_LMPD_214_pg181_c -4 F_LMPD_215_pg181_c +4 R_LMPD_215_pg181_c -4 F_LMPD_216_pg181_c +4 R_LMPD_216_pg181_c - F_LMPD_217_phe-L_c + R_LMPD_217_phe-L_c -8 F_LMPD_218_pheme_c +8 R_LMPD_218_pheme_c -12 F_LMPD_219_pheme_c +12 R_LMPD_219_pheme_c - F_LMPD_220_pro-L_c + R_LMPD_220_pro-L_c - F_LMPD_221_ptrc_c + R_LMPD_221_ptrc_c -31 F_LMPD_224_q8h2_c +31 R_LMPD_224_q8h2_c -33 F_LMPD_225_q8h2_c +33 R_LMPD_225_q8h2_c -37 F_LMPD_226_q8h2_c +37 R_LMPD_226_q8h2_c -34 F_LMPD_227_q8h2_c +34 R_LMPD_227_q8h2_c -36 F_LMPD_228_q8h2_c +36 R_LMPD_228_q8h2_c -36 F_LMPD_229_q8h2_c +36 R_LMPD_229_q8h2_c -40 F_LMPD_230_q8h2_c +40 R_LMPD_230_q8h2_c -40 F_LMPD_231_q8h2_c +40 R_LMPD_231_q8h2_c -70 F_LMPD_232_q8h2_c +70 R_LMPD_232_q8h2_c -76 F_LMPD_233_q8h2_c +76 R_LMPD_233_q8h2_c -76 F_LMPD_234_q8h2_c +76 R_LMPD_234_q8h2_c -11 F_LMPD_235_ribflv_c +11 R_LMPD_235_ribflv_c -12 F_LMPD_236_ribflv_c +12 R_LMPD_236_ribflv_c -12 F_LMPD_237_ribflv_c +12 R_LMPD_237_ribflv_c -13 F_LMPD_238_ribflv_c +13 R_LMPD_238_ribflv_c -14 F_LMPD_240_sheme_c +14 R_LMPD_240_sheme_c -6 F_LMPD_241_spmd_c +6 R_LMPD_241_spmd_c -7 F_LMPD_242_spmd_c +7 R_LMPD_242_spmd_c -15 F_LMPD_243_thf_c +15 R_LMPD_243_thf_c -16 F_LMPD_244_thf_c +16 R_LMPD_244_thf_c -16 F_LMPD_245_thf_c +16 R_LMPD_245_thf_c -17 F_LMPD_246_thf_c +17 R_LMPD_246_thf_c -14 F_LMPD_247_thmpp_c +14 R_LMPD_247_thmpp_c -15 F_LMPD_248_thmpp_c +15 R_LMPD_248_thmpp_c -2 F_LMPD_249_thr-L_c +2 R_LMPD_249_thr-L_c -3 F_LMPD_250_trp-L_c +3 R_LMPD_250_trp-L_c - F_LMPD_251_tyr-L_c + R_LMPD_251_tyr-L_c -33 F_LMPD_252_udcpdp_c +33 R_LMPD_252_udcpdp_c -4 F_LMPD_253_utp_c +4 R_LMPD_253_utp_c = 0
M_bglycogen_c: F_GLBRAN2 - R_GLBRAN2 - F_GLDBRAN2 + R_GLDBRAN2 = 0
M_bmocogdp_c: -0.000116 F_Ec_biomass_iJO1366_WT_53p95M +0.000116 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_35_bmocogdp_c - R_LMPD_35_bmocogdp_c + F_LMPD_36_bmocogdp_c - R_LMPD_36_bmocogdp_c + F_LMPD_37_bmocogdp_c - R_LMPD_37_bmocogdp_c + F_LMPD_38_bmocogdp_c - R_LMPD_38_bmocogdp_c = 0
M_btn_c: -2e-06 F_Ec_biomass_iJO1366_WT_53p95M +2e-06 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c + F_LMPD_43_btn_c - R_LMPD_43_btn_c + F_LMPD_44_btn_c - R_LMPD_44_btn_c + F_LMPD_45_btn_c - R_LMPD_45_btn_c + F_LMPD_47_btn_c - R_LMPD_47_btn_c = 0
M_ca2_c: -0.004952 F_Ec_biomass_iJO1366_WT_53p95M +0.004952 R_Ec_biomass_iJO1366_WT_53p95M - F_CAt6pp + R_CAt6pp = 0
M_chor_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_51_chor_c - R_LMPD_51_chor_c = 0
M_cit_c: - F_ACONTa + R_ACONTa - F_CITL + R_CITL - F_CITt3pp + R_CITt3pp + F_CITt7pp - R_CITt7pp + F_CS - R_CS = 0
M_cl_c: -0.004952 F_Ec_biomass_iJO1366_WT_53p95M +0.004952 R_Ec_biomass_iJO1366_WT_53p95M +2 F_CLt3_2pp -2 R_CLt3_2pp = 0
M_co2_c: F_AKGDH - R_AKGDH + F_CO2tpp - R_CO2tpp + F_GLXCL - R_GLXCL + F_GND - R_GND + F_ICDHyr - R_ICDHyr + F_ME1 - R_ME1 + F_ME2 - R_ME2 + F_OAADC - R_OAADC + F_PDH - R_PDH + F_POR5 - R_POR5 + F_POX - R_POX - F_PPC + R_PPC + F_PPCK - R_PPCK - F_LMPD_1_10fthf_c + R_LMPD_1_10fthf_c - F_LMPD_2_10fthf_c + R_LMPD_2_10fthf_c - F_LMPD_3_10fthf_c + R_LMPD_3_10fthf_c - F_LMPD_4_10fthf_c + R_LMPD_4_10fthf_c +10 F_LMPD_5_2dmmql8_c -10 R_LMPD_5_2dmmql8_c - F_LMPD_22_5mthf_c + R_LMPD_22_5mthf_c - F_LMPD_23_5mthf_c + R_LMPD_23_5mthf_c +2 F_LMPD_24_accoa_c -2 R_LMPD_24_accoa_c +2 F_LMPD_25_accoa_c -2 R_LMPD_25_accoa_c - F_LMPD_26_adocbl_c + R_LMPD_26_adocbl_c - F_LMPD_27_adocbl_c + R_LMPD_27_adocbl_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c - F_LMPD_33_arg-L_c + R_LMPD_33_arg-L_c -4 F_LMPD_35_bmocogdp_c +4 R_LMPD_35_bmocogdp_c -4 F_LMPD_36_bmocogdp_c +4 R_LMPD_36_bmocogdp_c -4 F_LMPD_37_bmocogdp_c +4 R_LMPD_37_bmocogdp_c -4 F_LMPD_38_bmocogdp_c +4 R_LMPD_38_bmocogdp_c -2 F_LMPD_39_btn_c +2 R_LMPD_39_btn_c -2 F_LMPD_40_btn_c +2 R_LMPD_40_btn_c -2 F_LMPD_41_btn_c +2 R_LMPD_41_btn_c -2 F_LMPD_42_btn_c +2 R_LMPD_42_btn_c -4 F_LMPD_43_btn_c +4 R_LMPD_43_btn_c -4 F_LMPD_44_btn_c +4 R_LMPD_44_btn_c -4 F_LMPD_45_btn_c +4 R_LMPD_45_btn_c -4 F_LMPD_47_btn_c +4 R_LMPD_47_btn_c - F_LMPD_94_datp_c + R_LMPD_94_datp_c - F_LMPD_95_datp_c + R_LMPD_95_datp_c - F_LMPD_96_datp_c + R_LMPD_96_datp_c - F_LMPD_97_datp_c + R_LMPD_97_datp_c - F_LMPD_99_dgtp_c + R_LMPD_99_dgtp_c - F_LMPD_100_dgtp_c + R_LMPD_100_dgtp_c - F_LMPD_101_dgtp_c + R_LMPD_101_dgtp_c - F_LMPD_102_dgtp_c + R_LMPD_102_dgtp_c - F_LMPD_103_dgtp_c + R_LMPD_103_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c -2 F_LMPD_108_fad_c +2 R_LMPD_108_fad_c -2 F_LMPD_109_fad_c +2 R_LMPD_109_fad_c -2 F_LMPD_110_fad_c +2 R_LMPD_110_fad_c -2 F_LMPD_111_fad_c +2 R_LMPD_111_fad_c - F_LMPD_112_GAM_c + R_LMPD_112_GAM_c - F_LMPD_113_GAM_c + R_LMPD_113_GAM_c - F_LMPD_114_GAM_c + R_LMPD_114_GAM_c - F_LMPD_115_GAM_c + R_LMPD_115_GAM_c - F_LMPD_120_gtp_c + R_LMPD_120_gtp_c - F_LMPD_121_gtp_c + R_LMPD_121_gtp_c - F_LMPD_122_gtp_c + R_LMPD_122_gtp_c - F_LMPD_123_gtp_c + R_LMPD_123_gtp_c +9 F_LMPD_124_hemeO_c -9 R_LMPD_124_hemeO_c +9 F_LMPD_125_hemeO_c -9 R_LMPD_125_hemeO_c + F_LMPD_128_ile-L_c - R_LMPD_128_ile-L_c +2 F_LMPD_129_leu-L_c -2 R_LMPD_129_leu-L_c + F_LMPD_134_lys-L_c - R_LMPD_134_lys-L_c + F_LMPD_135_malcoa_c - R_LMPD_135_malcoa_c + F_LMPD_136_malcoa_c - R_LMPD_136_malcoa_c - F_LMPD_139_mlthf_c + R_LMPD_139_mlthf_c - F_LMPD_140_mlthf_c + R_LMPD_140_mlthf_c - F_LMPD_141_mococdp_c + R_LMPD_141_mococdp_c - F_LMPD_142_mococdp_c + R_LMPD_142_mococdp_c - F_LMPD_143_mococdp_c + R_LMPD_143_mococdp_c - F_LMPD_144_mococdp_c + R_LMPD_144_mococdp_c - F_LMPD_145_mococdp_c + R_LMPD_145_mococdp_c - F_LMPD_146_mococdp_c + R_LMPD_146_mococdp_c - F_LMPD_147_mococdp_c + R_LMPD_147_mococdp_c - F_LMPD_148_mococdp_c + R_LMPD_148_mococdp_c -2 F_LMPD_149_mocogdp_c +2 R_LMPD_149_mocogdp_c -2 F_LMPD_150_mocogdp_c +2 R_LMPD_150_mocogdp_c -2 F_LMPD_151_mocogdp_c +2 R_LMPD_151_mocogdp_c -2 F_LMPD_152_mocogdp_c +2 R_LMPD_152_mocogdp_c +10 F_LMPD_153_mql8_c -10 R_LMPD_153_mql8_c + F_LMPD_163_pe160_c - R_LMPD_163_pe160_c + F_LMPD_164_pe160_c - R_LMPD_164_pe160_c + F_LMPD_165_pe160_c - R_LMPD_165_pe160_c + F_LMPD_166_pe160_c - R_LMPD_166_pe160_c + F_LMPD_167_pe160_c - R_LMPD_167_pe160_c + F_LMPD_168_pe160_c - R_LMPD_168_pe160_c + F_LMPD_169_pe160_c - R_LMPD_169_pe160_c + F_LMPD_170_pe160_c - R_LMPD_170_pe160_c + F_LMPD_171_pe160_c - R_LMPD_171_pe160_c + F_LMPD_172_pe160_c - R_LMPD_172_pe160_c + F_LMPD_173_pe161_c - R_LMPD_173_pe161_c + F_LMPD_174_pe161_c - R_LMPD_174_pe161_c + F_LMPD_175_pe161_c - R_LMPD_175_pe161_c + F_LMPD_176_pe161_c - R_LMPD_176_pe161_c + F_LMPD_177_pe161_c - R_LMPD_177_pe161_c + F_LMPD_178_pe161_c - R_LMPD_178_pe161_c + F_LMPD_179_pe161_c - R_LMPD_179_pe161_c + F_LMPD_180_pe181_c - R_LMPD_180_pe181_c + F_LMPD_181_pe181_c - R_LMPD_181_pe181_c + F_LMPD_182_pe181_c - R_LMPD_182_pe181_c + F_LMPD_183_pe181_c - R_LMPD_183_pe181_c + F_LMPD_184_pe181_c - R_LMPD_184_pe181_c + F_LMPD_185_pe181_c - R_LMPD_185_pe181_c + F_LMPD_186_pe181_c - R_LMPD_186_pe181_c + F_LMPD_187_pe181_c - R_LMPD_187_pe181_c + F_LMPD_188_pe181_c - R_LMPD_188_pe181_c + F_LMPD_189_pe181_c - R_LMPD_189_pe181_c + F_LMPD_217_phe-L_c - R_LMPD_217_phe-L_c +6 F_LMPD_218_pheme_c -6 R_LMPD_218_pheme_c +6 F_LMPD_219_pheme_c -6 R_LMPD_219_pheme_c + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c +2 F_LMPD_222_pydx5p_c -2 R_LMPD_222_pydx5p_c +2 F_LMPD_223_pydx5p_c -2 R_LMPD_223_pydx5p_c +9 F_LMPD_224_q8h2_c -9 R_LMPD_224_q8h2_c +9 F_LMPD_225_q8h2_c -9 R_LMPD_225_q8h2_c +9 F_LMPD_226_q8h2_c -9 R_LMPD_226_q8h2_c +9 F_LMPD_227_q8h2_c -9 R_LMPD_227_q8h2_c +9 F_LMPD_228_q8h2_c -9 R_LMPD_228_q8h2_c +9 F_LMPD_229_q8h2_c -9 R_LMPD_229_q8h2_c +9 F_LMPD_230_q8h2_c -9 R_LMPD_230_q8h2_c +9 F_LMPD_231_q8h2_c -9 R_LMPD_231_q8h2_c +18 F_LMPD_232_q8h2_c -18 R_LMPD_232_q8h2_c +18 F_LMPD_233_q8h2_c -18 R_LMPD_233_q8h2_c +18 F_LMPD_234_q8h2_c -18 R_LMPD_234_q8h2_c - F_LMPD_235_ribflv_c + R_LMPD_235_ribflv_c - F_LMPD_236_ribflv_c + R_LMPD_236_ribflv_c - F_LMPD_237_ribflv_c + R_LMPD_237_ribflv_c - F_LMPD_238_ribflv_c + R_LMPD_238_ribflv_c +2 F_LMPD_241_spmd_c -2 R_LMPD_241_spmd_c +2 F_LMPD_242_spmd_c -2 R_LMPD_242_spmd_c - F_LMPD_243_thf_c + R_LMPD_243_thf_c - F_LMPD_244_thf_c + R_LMPD_244_thf_c - F_LMPD_245_thf_c + R_LMPD_245_thf_c - F_LMPD_246_thf_c + R_LMPD_246_thf_c +3 F_LMPD_247_thmpp_c -3 R_LMPD_247_thmpp_c +3 F_LMPD_248_thmpp_c -3 R_LMPD_248_thmpp_c + F_LMPD_250_trp-L_c - R_LMPD_250_trp-L_c + F_LMPD_251_tyr-L_c - R_LMPD_251_tyr-L_c +11 F_LMPD_252_udcpdp_c -11 R_LMPD_252_udcpdp_c + F_LMPD_254_val-L_c - R_LMPD_254_val-L_c = 0
M_coa_c: -0.000168 F_Ec_biomass_iJO1366_WT_53p95M +0.000168 R_Ec_biomass_iJO1366_WT_53p95M - F_ACALD + R_ACALD - F_ACS + R_ACS - F_AKGDH + R_AKGDH + F_CS - R_CS + F_MALS - R_MALS - F_PDH + R_PDH - F_PFL + R_PFL - F_POR5 + R_POR5 + F_PTAr - R_PTAr - F_SUCOAS + R_SUCOAS +2 F_LMPD_6_2fe2s_c -2 R_LMPD_6_2fe2s_c +2 F_LMPD_7_2fe2s_c -2 R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c +2 F_LMPD_9_2fe2s_c -2 R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c +2 F_LMPD_11_2fe2s_c -2 R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +4 F_LMPD_14_4fe4s_c -4 R_LMPD_14_4fe4s_c +4 F_LMPD_15_4fe4s_c -4 R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c +4 F_LMPD_17_4fe4s_c -4 R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c +4 F_LMPD_19_4fe4s_c -4 R_LMPD_19_4fe4s_c +4 F_LMPD_20_4fe4s_c -4 R_LMPD_20_4fe4s_c +4 F_LMPD_21_4fe4s_c -4 R_LMPD_21_4fe4s_c +2 F_LMPD_24_accoa_c -2 R_LMPD_24_accoa_c +2 F_LMPD_25_accoa_c -2 R_LMPD_25_accoa_c +2 F_LMPD_31_amet_c -2 R_LMPD_31_amet_c +2 F_LMPD_32_amet_c -2 R_LMPD_32_amet_c + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c +4 F_LMPD_35_bmocogdp_c -4 R_LMPD_35_bmocogdp_c +4 F_LMPD_36_bmocogdp_c -4 R_LMPD_36_bmocogdp_c +4 F_LMPD_37_bmocogdp_c -4 R_LMPD_37_bmocogdp_c +4 F_LMPD_38_bmocogdp_c -4 R_LMPD_38_bmocogdp_c +6 F_LMPD_39_btn_c -6 R_LMPD_39_btn_c +6 F_LMPD_40_btn_c -6 R_LMPD_40_btn_c +6 F_LMPD_41_btn_c -6 R_LMPD_41_btn_c +6 F_LMPD_42_btn_c -6 R_LMPD_42_btn_c +8 F_LMPD_43_btn_c -8 R_LMPD_43_btn_c +8 F_LMPD_44_btn_c -8 R_LMPD_44_btn_c +8 F_LMPD_45_btn_c -8 R_LMPD_45_btn_c +8 F_LMPD_47_btn_c -8 R_LMPD_47_btn_c + F_LMPD_48_ca2_c - R_LMPD_48_ca2_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c +32 F_LMPD_52_clpn160_p -32 R_LMPD_52_clpn160_p +32 F_LMPD_53_clpn160_p -32 R_LMPD_53_clpn160_p +32 F_LMPD_54_clpn160_p -32 R_LMPD_54_clpn160_p +32 F_LMPD_55_clpn160_p -32 R_LMPD_55_clpn160_p +32 F_LMPD_56_clpn160_p -32 R_LMPD_56_clpn160_p +32 F_LMPD_57_clpn160_p -32 R_LMPD_57_clpn160_p +32 F_LMPD_58_clpn160_p -32 R_LMPD_58_clpn160_p +32 F_LMPD_59_clpn160_p -32 R_LMPD_59_clpn160_p +32 F_LMPD_60_clpn160_p -32 R_LMPD_60_clpn160_p +32 F_LMPD_61_clpn160_p -32 R_LMPD_61_clpn160_p +32 F_LMPD_62_clpn161_p -32 R_LMPD_62_clpn161_p +32 F_LMPD_63_clpn161_p -32 R_LMPD_63_clpn161_p +32 F_LMPD_64_clpn161_p -32 R_LMPD_64_clpn161_p +32 F_LMPD_65_clpn161_p -32 R_LMPD_65_clpn161_p +32 F_LMPD_66_clpn161_p -32 R_LMPD_66_clpn161_p +32 F_LMPD_67_clpn161_p -32 R_LMPD_67_clpn161_p +32 F_LMPD_68_clpn161_p -32 R_LMPD_68_clpn161_p +36 F_LMPD_69_clpn181_p -36 R_LMPD_69_clpn181_p +36 F_LMPD_70_clpn181_p -36 R_LMPD_70_clpn181_p +36 F_LMPD_71_clpn181_p -36 R_LMPD_71_clpn181_p +36 F_LMPD_72_clpn181_p -36 R_LMPD_72_clpn181_p +36 F_LMPD_73_clpn181_p -36 R_LMPD_73_clpn181_p +36 F_LMPD_74_clpn181_p -36 R_LMPD_74_clpn181_p +36 F_LMPD_75_clpn181_p -36 R_LMPD_75_clpn181_p +36 F_LMPD_76_clpn181_p -36 R_LMPD_76_clpn181_p +36 F_LMPD_77_clpn181_p -36 R_LMPD_77_clpn181_p +36 F_LMPD_78_clpn181_p -36 R_LMPD_78_clpn181_p +43 F_LMPD_79_colipa_e -43 R_LMPD_79_colipa_e +43 F_LMPD_80_colipa_e -43 R_LMPD_80_colipa_e +43 F_LMPD_81_colipa_e -43 R_LMPD_81_colipa_e +43 F_LMPD_82_colipa_e -43 R_LMPD_82_colipa_e +43 F_LMPD_83_colipa_e -43 R_LMPD_83_colipa_e +43 F_LMPD_84_colipa_e -43 R_LMPD_84_colipa_e +43 F_LMPD_85_colipa_e -43 R_LMPD_85_colipa_e +43 F_LMPD_86_colipa_e -43 R_LMPD_86_colipa_e +43 F_LMPD_87_colipa_e -43 R_LMPD_87_colipa_e +43 F_LMPD_88_colipa_e -43 R_LMPD_88_colipa_e +43 F_LMPD_89_colipa_e -43 R_LMPD_89_colipa_e +43 F_LMPD_90_colipa_e -43 R_LMPD_90_colipa_e + F_LMPD_92_cys-L_c - R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c + F_LMPD_118_gthrd_c - R_LMPD_118_gthrd_c + F_LMPD_119_gthrd_c - R_LMPD_119_gthrd_c + F_LMPD_129_leu-L_c - R_LMPD_129_leu-L_c +8 F_LMPD_130_lipopb_c -8 R_LMPD_130_lipopb_c +8 F_LMPD_131_lipopb_c -8 R_LMPD_131_lipopb_c +8 F_LMPD_132_lipopb_c -8 R_LMPD_132_lipopb_c +8 F_LMPD_133_lipopb_c -8 R_LMPD_133_lipopb_c + F_LMPD_134_lys-L_c - R_LMPD_134_lys-L_c +2 F_LMPD_135_malcoa_c -2 R_LMPD_135_malcoa_c +2 F_LMPD_136_malcoa_c -2 R_LMPD_136_malcoa_c +2 F_LMPD_137_met-L_c -2 R_LMPD_137_met-L_c +2 F_LMPD_138_met-L_c -2 R_LMPD_138_met-L_c +2 F_LMPD_141_mococdp_c -2 R_LMPD_141_mococdp_c +2 F_LMPD_142_mococdp_c -2 R_LMPD_142_mococdp_c +2 F_LMPD_143_mococdp_c -2 R_LMPD_143_mococdp_c +2 F_LMPD_144_mococdp_c -2 R_LMPD_144_mococdp_c +2 F_LMPD_145_mococdp_c -2 R_LMPD_145_mococdp_c +2 F_LMPD_146_mococdp_c -2 R_LMPD_146_mococdp_c +2 F_LMPD_147_mococdp_c -2 R_LMPD_147_mococdp_c +2 F_LMPD_148_mococdp_c -2 R_LMPD_148_mococdp_c +2 F_LMPD_149_mocogdp_c -2 R_LMPD_149_mocogdp_c +2 F_LMPD_150_mocogdp_c -2 R_LMPD_150_mocogdp_c +2 F_LMPD_151_mocogdp_c -2 R_LMPD_151_mocogdp_c +2 F_LMPD_152_mocogdp_c -2 R_LMPD_152_mocogdp_c +6 F_LMPD_154_murein3p3p_p -6 R_LMPD_154_murein3p3p_p +6 F_LMPD_155_murein3px4p_p -6 R_LMPD_155_murein3px4p_p +6 F_LMPD_156_murein4p4p_p -6 R_LMPD_156_murein4p4p_p +6 F_LMPD_157_murein4px4p_p -6 R_LMPD_157_murein4px4p_p +9 F_LMPD_158_murein4px4px4p_p -9 R_LMPD_158_murein4px4px4p_p +16 F_LMPD_163_pe160_c -16 R_LMPD_163_pe160_c +16 F_LMPD_164_pe160_c -16 R_LMPD_164_pe160_c +16 F_LMPD_165_pe160_c -16 R_LMPD_165_pe160_c +16 F_LMPD_166_pe160_c -16 R_LMPD_166_pe160_c +16 F_LMPD_167_pe160_c -16 R_LMPD_167_pe160_c +16 F_LMPD_168_pe160_c -16 R_LMPD_168_pe160_c +16 F_LMPD_169_pe160_c -16 R_LMPD_169_pe160_c +16 F_LMPD_170_pe160_c -16 R_LMPD_170_pe160_c +16 F_LMPD_171_pe160_c -16 R_LMPD_171_pe160_c +16 F_LMPD_172_pe160_c -16 R_LMPD_172_pe160_c +16 F_LMPD_173_pe161_c -16 R_LMPD_173_pe161_c +16 F_LMPD_174_pe161_c -16 R_LMPD_174_pe161_c +16 F_LMPD_175_pe161_c -16 R_LMPD_175_pe161_c +16 F_LMPD_176_pe161_c -16 R_LMPD_176_pe161_c +16 F_LMPD_177_pe161_c -16 R_LMPD_177_pe161_c +16 F_LMPD_178_pe161_c -16 R_LMPD_178_pe161_c +16 F_LMPD_179_pe161_c -16 R_LMPD_179_pe161_c +18 F_LMPD_180_pe181_c -18 R_LMPD_180_pe181_c +18 F_LMPD_181_pe181_c -18 R_LMPD_181_pe181_c +18 F_LMPD_182_pe181_c -18 R_LMPD_182_pe181_c +18 F_LMPD_183_pe181_c -18 R_LMPD_183_pe181_c +18 F_LMPD_184_pe181_c -18 R_LMPD_184_pe181_c +18 F_LMPD_185_pe181_c -18 R_LMPD_185_pe181_c +18 F_LMPD_186_pe181_c -18 R_LMPD_186_pe181_c +18 F_LMPD_187_pe181_c -18 R_LMPD_187_pe181_c +18 F_LMPD_188_pe181_c -18 R_LMPD_188_pe181_c +18 F_LMPD_189_pe181_c -18 R_LMPD_189_pe181_c +16 F_LMPD_190_pg160_c -16 R_LMPD_190_pg160_c +16 F_LMPD_191_pg160_c -16 R_LMPD_191_pg160_c +16 F_LMPD_192_pg160_c -16 R_LMPD_192_pg160_c +16 F_LMPD_193_pg160_c -16 R_LMPD_193_pg160_c +16 F_LMPD_194_pg160_c -16 R_LMPD_194_pg160_c +16 F_LMPD_195_pg160_c -16 R_LMPD_195_pg160_c +16 F_LMPD_196_pg160_c -16 R_LMPD_196_pg160_c +16 F_LMPD_197_pg160_c -16 R_LMPD_197_pg160_c +16 F_LMPD_198_pg160_c -16 R_LMPD_198_pg160_c +16 F_LMPD_199_pg160_c -16 R_LMPD_199_pg160_c +16 F_LMPD_200_pg161_c -16 R_LMPD_200_pg161_c +16 F_LMPD_201_pg161_c -16 R_LMPD_201_pg161_c +16 F_LMPD_202_pg161_c -16 R_LMPD_202_pg161_c +16 F_LMPD_203_pg161_c -16 R_LMPD_203_pg161_c +16 F_LMPD_204_pg161_c -16 R_LMPD_204_pg161_c +16 F_LMPD_205_pg161_c -16 R_LMPD_205_pg161_c +16 F_LMPD_206_pg161_c -16 R_LMPD_206_pg161_c +18 F_LMPD_207_pg181_c -18 R_LMPD_207_pg181_c +18 F_LMPD_208_pg181_c -18 R_LMPD_208_pg181_c +18 F_LMPD_209_pg181_c -18 R_LMPD_209_pg181_c +18 F_LMPD_210_pg181_c -18 R_LMPD_210_pg181_c +18 F_LMPD_211_pg181_c -18 R_LMPD_211_pg181_c +18 F_LMPD_212_pg181_c -18 R_LMPD_212_pg181_c +18 F_LMPD_213_pg181_c -18 R_LMPD_213_pg181_c +18 F_LMPD_214_pg181_c -18 R_LMPD_214_pg181_c +18 F_LMPD_215_pg181_c -18 R_LMPD_215_pg181_c +18 F_LMPD_216_pg181_c -18 R_LMPD_216_pg181_c + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c +3 F_LMPD_241_spmd_c -3 R_LMPD_241_spmd_c +3 F_LMPD_242_spmd_c -3 R_LMPD_242_spmd_c + F_LMPD_247_thmpp_c - R_LMPD_247_thmpp_c + F_LMPD_248_thmpp_c - R_LMPD_248_thmpp_c = 0
M_cobalt2_c: -2.4e-05 F_Ec_biomass_iJO1366_WT_53p95M +2.4e-05 R_Ec_biomass_iJO1366_WT_53p95M + F_COBALT2tpp - R_COBALT2tpp = 0
M_colipa_c: - F_COLIPAabcpp + R_COLIPAabcpp + F_LMPD_79_colipa_e - R_LMPD_79_colipa_e + F_LMPD_80_colipa_e - R_LMPD_80_colipa_e + F_LMPD_81_colipa_e - R_LMPD_81_colipa_e + F_LMPD_82_colipa_e - R_LMPD_82_colipa_e + F_LMPD_83_colipa_e - R_LMPD_83_colipa_e + F_LMPD_84_colipa_e - R_LMPD_84_colipa_e + F_LMPD_85_colipa_e - R_LMPD_85_colipa_e + F_LMPD_86_colipa_e - R_LMPD_86_colipa_e + F_LMPD_87_colipa_e - R_LMPD_87_colipa_e + F_LMPD_88_colipa_e - R_LMPD_88_colipa_e + F_LMPD_89_colipa_e - R_LMPD_89_colipa_e + F_LMPD_90_colipa_e - R_LMPD_90_colipa_e = 0
M_ctp_c: -0.1298 F_Ec_biomass_iJO1366_WT_53p95M +0.1298 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_91_ctp_c - R_LMPD_91_ctp_c = 0
M_cu2_c: -0.000674 F_Ec_biomass_iJO1366_WT_53p95M +0.000674 R_Ec_biomass_iJO1366_WT_53p95M + F_CU2tpp - R_CU2tpp = 0
M_cys-L_c: -0.088988 F_Ec_biomass_iJO1366_WT_53p95M +0.088988 R_Ec_biomass_iJO1366_WT_53p95M - F_CYSabc2pp + R_CYSabc2pp - F_CYStpp + R_CYStpp + F_LMPD_92_cys-L_c - R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c = 0
M_datp_c: -0.024805 F_Ec_biomass_iJO1366_WT_53p95M +0.024805 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_94_datp_c - R_LMPD_94_datp_c + F_LMPD_95_datp_c - R_LMPD_95_datp_c + F_LMPD_96_datp_c - R_LMPD_96_datp_c + F_LMPD_97_datp_c - R_LMPD_97_datp_c = 0
M_db4p_c: F_DB4PS - R_DB4PS -2 F_LMPD_108_fad_c +2 R_LMPD_108_fad_c -2 F_LMPD_109_fad_c +2 R_LMPD_109_fad_c -2 F_LMPD_110_fad_c +2 R_LMPD_110_fad_c -2 F_LMPD_111_fad_c +2 R_LMPD_111_fad_c -2 F_LMPD_235_ribflv_c +2 R_LMPD_235_ribflv_c -2 F_LMPD_236_ribflv_c +2 R_LMPD_236_ribflv_c -2 F_LMPD_237_ribflv_c +2 R_LMPD_237_ribflv_c -2 F_LMPD_238_ribflv_c +2 R_LMPD_238_ribflv_c = 0
M_dctp_c: -0.025612 F_Ec_biomass_iJO1366_WT_53p95M +0.025612 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_98_dctp_c - R_LMPD_98_dctp_c = 0
M_dgtp_c: -0.025612 F_Ec_biomass_iJO1366_WT_53p95M +0.025612 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_99_dgtp_c - R_LMPD_99_dgtp_c + F_LMPD_100_dgtp_c - R_LMPD_100_dgtp_c + F_LMPD_101_dgtp_c - R_LMPD_101_dgtp_c + F_LMPD_102_dgtp_c - R_LMPD_102_dgtp_c + F_LMPD_103_dgtp_c - R_LMPD_103_dgtp_c + F_LMPD_104_dgtp_c - R_LMPD_104_dgtp_c = 0
M_dha_c: - F_DHAPT + R_DHAPT + F_DHAtpp - R_DHAtpp + F_F6PA - R_F6PA + F_GLYCDx - R_GLYCDx = 0
M_dhap_c: F_DHAPT - R_DHAPT + F_FBA - R_FBA + F_FBA3 - R_FBA3 + F_G3PD2 - R_G3PD2 + F_G3PD5 - R_G3PD5 + F_G3PD6 - R_G3PD6 + F_G3PD7 - R_G3PD7 - F_TPI + R_TPI - F_LMPD_159_nad_c + R_LMPD_159_nad_c - F_LMPD_160_nad_c + R_LMPD_160_nad_c - F_LMPD_161_nad_c + R_LMPD_161_nad_c - F_LMPD_162_nad_c + R_LMPD_162_nad_c = 0
M_dhor-S_c: - F_DHORD2 + R_DHORD2 - F_DHORD5 + R_DHORD5 - F_DHORDfum + R_DHORDfum + F_LMPD_91_ctp_c - R_LMPD_91_ctp_c + F_LMPD_98_dctp_c - R_LMPD_98_dctp_c + F_LMPD_105_dttp_c - R_LMPD_105_dttp_c + F_LMPD_106_dttp_c - R_LMPD_106_dttp_c + F_LMPD_141_mococdp_c - R_LMPD_141_mococdp_c + F_LMPD_142_mococdp_c - R_LMPD_142_mococdp_c + F_LMPD_143_mococdp_c - R_LMPD_143_mococdp_c + F_LMPD_144_mococdp_c - R_LMPD_144_mococdp_c + F_LMPD_145_mococdp_c - R_LMPD_145_mococdp_c + F_LMPD_146_mococdp_c - R_LMPD_146_mococdp_c + F_LMPD_147_mococdp_c - R_LMPD_147_mococdp_c + F_LMPD_148_mococdp_c - R_LMPD_148_mococdp_c + F_LMPD_253_utp_c - R_LMPD_253_utp_c = 0
M_dttp_c: -0.024805 F_Ec_biomass_iJO1366_WT_53p95M +0.024805 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_105_dttp_c - R_LMPD_105_dttp_c + F_LMPD_106_dttp_c - R_LMPD_106_dttp_c = 0
M_e4p_c: F_FBA3 - R_FBA3 + F_TALA - R_TALA - F_TKT2 + R_TKT2 - F_LMPD_1_10fthf_c + R_LMPD_1_10fthf_c - F_LMPD_2_10fthf_c + R_LMPD_2_10fthf_c - F_LMPD_3_10fthf_c + R_LMPD_3_10fthf_c - F_LMPD_4_10fthf_c + R_LMPD_4_10fthf_c - F_LMPD_5_2dmmql8_c + R_LMPD_5_2dmmql8_c - F_LMPD_22_5mthf_c + R_LMPD_22_5mthf_c - F_LMPD_23_5mthf_c + R_LMPD_23_5mthf_c - F_LMPD_51_chor_c + R_LMPD_51_chor_c -3 F_LMPD_107_enter_c +3 R_LMPD_107_enter_c - F_LMPD_139_mlthf_c + R_LMPD_139_mlthf_c - F_LMPD_140_mlthf_c + R_LMPD_140_mlthf_c - F_LMPD_153_mql8_c + R_LMPD_153_mql8_c - F_LMPD_217_phe-L_c + R_LMPD_217_phe-L_c - F_LMPD_222_pydx5p_c + R_LMPD_222_pydx5p_c - F_LMPD_223_pydx5p_c + R_LMPD_223_pydx5p_c - F_LMPD_224_q8h2_c + R_LMPD_224_q8h2_c - F_LMPD_225_q8h2_c + R_LMPD_225_q8h2_c - F_LMPD_226_q8h2_c + R_LMPD_226_q8h2_c - F_LMPD_227_q8h2_c + R_LMPD_227_q8h2_c - F_LMPD_228_q8h2_c + R_LMPD_228_q8h2_c - F_LMPD_229_q8h2_c + R_LMPD_229_q8h2_c - F_LMPD_230_q8h2_c + R_LMPD_230_q8h2_c - F_LMPD_231_q8h2_c + R_LMPD_231_q8h2_c -2 F_LMPD_232_q8h2_c +2 R_LMPD_232_q8h2_c -2 F_LMPD_233_q8h2_c +2 R_LMPD_233_q8h2_c -2 F_LMPD_234_q8h2_c +2 R_LMPD_234_q8h2_c - F_LMPD_243_thf_c + R_LMPD_243_thf_c - F_LMPD_244_thf_c + R_LMPD_244_thf_c - F_LMPD_245_thf_c + R_LMPD_245_thf_c - F_LMPD_246_thf_c + R_LMPD_246_thf_c - F_LMPD_247_thmpp_c + R_LMPD_247_thmpp_c - F_LMPD_248_thmpp_c + R_LMPD_248_thmpp_c - F_LMPD_250_trp-L_c + R_LMPD_250_trp-L_c - F_LMPD_251_tyr-L_c + R_LMPD_251_tyr-L_c = 0
M_enter_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M - F_FEENTERtpp + R_FEENTERtpp + F_LMPD_107_enter_c - R_LMPD_107_enter_c = 0
M_etoh_c: - F_ALCD2x + R_ALCD2x + F_ETOHtrpp - R_ETOHtrpp = 0
M_f6p_c: - F_F6PA + R_F6PA + F_FBP - R_FBP - F_PFK + R_PFK + F_PGI - R_PGI + F_TALA - R_TALA + F_TKT2 - R_TKT2 -2 F_LMPD_79_colipa_e +2 R_LMPD_79_colipa_e -2 F_LMPD_80_colipa_e +2 R_LMPD_80_colipa_e -2 F_LMPD_81_colipa_e +2 R_LMPD_81_colipa_e -2 F_LMPD_82_colipa_e +2 R_LMPD_82_colipa_e -2 F_LMPD_83_colipa_e +2 R_LMPD_83_colipa_e -2 F_LMPD_84_colipa_e +2 R_LMPD_84_colipa_e -2 F_LMPD_85_colipa_e +2 R_LMPD_85_colipa_e -2 F_LMPD_86_colipa_e +2 R_LMPD_86_colipa_e -2 F_LMPD_87_colipa_e +2 R_LMPD_87_colipa_e -2 F_LMPD_88_colipa_e +2 R_LMPD_88_colipa_e -2 F_LMPD_89_colipa_e +2 R_LMPD_89_colipa_e -2 F_LMPD_90_colipa_e +2 R_LMPD_90_colipa_e -4 F_LMPD_154_murein3p3p_p +4 R_LMPD_154_murein3p3p_p -4 F_LMPD_155_murein3px4p_p +4 R_LMPD_155_murein3px4p_p -4 F_LMPD_156_murein4p4p_p +4 R_LMPD_156_murein4p4p_p -4 F_LMPD_157_murein4px4p_p +4 R_LMPD_157_murein4px4p_p -6 F_LMPD_158_murein4px4px4p_p +6 R_LMPD_158_murein4px4px4p_p = 0
M_fad_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_108_fad_c - R_LMPD_108_fad_c + F_LMPD_109_fad_c - R_LMPD_109_fad_c + F_LMPD_110_fad_c - R_LMPD_110_fad_c + F_LMPD_111_fad_c - R_LMPD_111_fad_c = 0
M_fdp_c: - F_FBA + R_FBA - F_FBP + R_FBP + F_PFK - R_PFK = 0
M_fe2_c: -0.006388 F_Ec_biomass_iJO1366_WT_53p95M +0.006388 R_Ec_biomass_iJO1366_WT_53p95M + F_FE2abcpp - R_FE2abcpp + F_FE2t2pp - R_FE2t2pp + F_FE2tpp - R_FE2tpp -2 F_LMPD_6_2fe2s_c +2 R_LMPD_6_2fe2s_c -2 F_LMPD_7_2fe2s_c +2 R_LMPD_7_2fe2s_c -2 F_LMPD_8_2fe2s_c +2 R_LMPD_8_2fe2s_c -2 F_LMPD_9_2fe2s_c +2 R_LMPD_9_2fe2s_c -2 F_LMPD_10_2fe2s_c +2 R_LMPD_10_2fe2s_c -2 F_LMPD_11_2fe2s_c +2 R_LMPD_11_2fe2s_c -2 F_LMPD_12_2fe2s_c +2 R_LMPD_12_2fe2s_c -2 F_LMPD_13_2fe2s_c +2 R_LMPD_13_2fe2s_c -4 F_LMPD_14_4fe4s_c +4 R_LMPD_14_4fe4s_c -4 F_LMPD_15_4fe4s_c +4 R_LMPD_15_4fe4s_c -4 F_LMPD_16_4fe4s_c +4 R_LMPD_16_4fe4s_c -4 F_LMPD_17_4fe4s_c +4 R_LMPD_17_4fe4s_c -4 F_LMPD_18_4fe4s_c +4 R_LMPD_18_4fe4s_c -4 F_LMPD_19_4fe4s_c +4 R_LMPD_19_4fe4s_c -4 F_LMPD_20_4fe4s_c +4 R_LMPD_20_4fe4s_c -4 F_LMPD_21_4fe4s_c +4 R_LMPD_21_4fe4s_c - F_LMPD_124_hemeO_c + R_LMPD_124_hemeO_c - F_LMPD_125_hemeO_c + R_LMPD_125_hemeO_c -2 F_LMPD_130_lipopb_c +2 R_LMPD_130_lipopb_c -2 F_LMPD_131_lipopb_c +2 R_LMPD_131_lipopb_c -2 F_LMPD_132_lipopb_c +2 R_LMPD_132_lipopb_c -2 F_LMPD_133_lipopb_c +2 R_LMPD_133_lipopb_c - F_LMPD_218_pheme_c + R_LMPD_218_pheme_c - F_LMPD_219_pheme_c + R_LMPD_219_pheme_c - F_LMPD_240_sheme_c + R_LMPD_240_sheme_c = 0
M_fe3_c: -0.007428 F_Ec_biomass_iJO1366_WT_53p95M +0.007428 R_Ec_biomass_iJO1366_WT_53p95M + F_FE3abcpp - R_FE3abcpp = 0
M_flxr_c: +2 F_FLDR2 -2 R_FLDR2 +2 F_POR5 -2 R_POR5 -16 F_LMPD_5_2dmmql8_c +16 R_LMPD_5_2dmmql8_c -2 F_LMPD_94_datp_c +2 R_LMPD_94_datp_c -2 F_LMPD_95_datp_c +2 R_LMPD_95_datp_c -2 F_LMPD_96_datp_c +2 R_LMPD_96_datp_c -2 F_LMPD_97_datp_c +2 R_LMPD_97_datp_c -2 F_LMPD_98_dctp_c +2 R_LMPD_98_dctp_c -2 F_LMPD_99_dgtp_c +2 R_LMPD_99_dgtp_c -2 F_LMPD_100_dgtp_c +2 R_LMPD_100_dgtp_c -2 F_LMPD_102_dgtp_c +2 R_LMPD_102_dgtp_c -2 F_LMPD_103_dgtp_c +2 R_LMPD_103_dgtp_c -2 F_LMPD_105_dttp_c +2 R_LMPD_105_dttp_c -6 F_LMPD_124_hemeO_c +6 R_LMPD_124_hemeO_c -6 F_LMPD_125_hemeO_c +6 R_LMPD_125_hemeO_c -16 F_LMPD_153_mql8_c +16 R_LMPD_153_mql8_c -16 F_LMPD_224_q8h2_c +16 R_LMPD_224_q8h2_c -16 F_LMPD_225_q8h2_c +16 R_LMPD_225_q8h2_c -16 F_LMPD_226_q8h2_c +16 R_LMPD_226_q8h2_c -16 F_LMPD_227_q8h2_c +16 R_LMPD_227_q8h2_c -16 F_LMPD_228_q8h2_c +16 R_LMPD_228_q8h2_c -16 F_LMPD_229_q8h2_c +16 R_LMPD_229_q8h2_c -16 F_LMPD_230_q8h2_c +16 R_LMPD_230_q8h2_c -16 F_LMPD_231_q8h2_c +16 R_LMPD_231_q8h2_c -32 F_LMPD_232_q8h2_c +32 R_LMPD_232_q8h2_c -32 F_LMPD_233_q8h2_c +32 R_LMPD_233_q8h2_c -32 F_LMPD_234_q8h2_c +32 R_LMPD_234_q8h2_c -22 F_LMPD_252_udcpdp_c +22 R_LMPD_252_udcpdp_c = 0
M_flxso_c: -2 F_FLDR2 +2 R_FLDR2 -2 F_POR5 +2 R_POR5 +16 F_LMPD_5_2dmmql8_c -16 R_LMPD_5_2dmmql8_c +2 F_LMPD_94_datp_c -2 R_LMPD_94_datp_c +2 F_LMPD_95_datp_c -2 R_LMPD_95_datp_c +2 F_LMPD_96_datp_c -2 R_LMPD_96_datp_c +2 F_LMPD_97_datp_c -2 R_LMPD_97_datp_c +2 F_LMPD_98_dctp_c -2 R_LMPD_98_dctp_c +2 F_LMPD_99_dgtp_c -2 R_LMPD_99_dgtp_c +2 F_LMPD_100_dgtp_c -2 R_LMPD_100_dgtp_c +2 F_LMPD_102_dgtp_c -2 R_LMPD_102_dgtp_c +2 F_LMPD_103_dgtp_c -2 R_LMPD_103_dgtp_c +2 F_LMPD_105_dttp_c -2 R_LMPD_105_dttp_c +6 F_LMPD_124_hemeO_c -6 R_LMPD_124_hemeO_c +6 F_LMPD_125_hemeO_c -6 R_LMPD_125_hemeO_c +16 F_LMPD_153_mql8_c -16 R_LMPD_153_mql8_c +16 F_LMPD_224_q8h2_c -16 R_LMPD_224_q8h2_c +16 F_LMPD_225_q8h2_c -16 R_LMPD_225_q8h2_c +16 F_LMPD_226_q8h2_c -16 R_LMPD_226_q8h2_c +16 F_LMPD_227_q8h2_c -16 R_LMPD_227_q8h2_c +16 F_LMPD_228_q8h2_c -16 R_LMPD_228_q8h2_c +16 F_LMPD_229_q8h2_c -16 R_LMPD_229_q8h2_c +16 F_LMPD_230_q8h2_c -16 R_LMPD_230_q8h2_c +16 F_LMPD_231_q8h2_c -16 R_LMPD_231_q8h2_c +32 F_LMPD_232_q8h2_c -32 R_LMPD_232_q8h2_c +32 F_LMPD_233_q8h2_c -32 R_LMPD_233_q8h2_c +32 F_LMPD_234_q8h2_c -32 R_LMPD_234_q8h2_c +22 F_LMPD_252_udcpdp_c -22 R_LMPD_252_udcpdp_c = 0
M_for_c: F_DB4PS - R_DB4PS - F_FORtppi + R_FORtppi + F_PFL - R_PFL + F_LMPD_3_10fthf_c - R_LMPD_3_10fthf_c + F_LMPD_4_10fthf_c - R_LMPD_4_10fthf_c - F_LMPD_22_5mthf_c + R_LMPD_22_5mthf_c - F_LMPD_23_5mthf_c + R_LMPD_23_5mthf_c -2 F_LMPD_24_accoa_c +2 R_LMPD_24_accoa_c -2 F_LMPD_25_accoa_c +2 R_LMPD_25_accoa_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c -2 F_LMPD_31_amet_c +2 R_LMPD_31_amet_c -2 F_LMPD_32_amet_c +2 R_LMPD_32_amet_c -4 F_LMPD_37_bmocogdp_c +4 R_LMPD_37_bmocogdp_c -4 F_LMPD_38_bmocogdp_c +4 R_LMPD_38_bmocogdp_c -2 F_LMPD_39_btn_c +2 R_LMPD_39_btn_c -2 F_LMPD_40_btn_c +2 R_LMPD_40_btn_c -2 F_LMPD_41_btn_c +2 R_LMPD_41_btn_c -2 F_LMPD_42_btn_c +2 R_LMPD_42_btn_c -6 F_LMPD_43_btn_c +6 R_LMPD_43_btn_c -6 F_LMPD_44_btn_c +6 R_LMPD_44_btn_c -6 F_LMPD_45_btn_c +6 R_LMPD_45_btn_c -6 F_LMPD_47_btn_c +6 R_LMPD_47_btn_c - F_LMPD_96_datp_c + R_LMPD_96_datp_c - F_LMPD_97_datp_c + R_LMPD_97_datp_c - F_LMPD_102_dgtp_c + R_LMPD_102_dgtp_c - F_LMPD_103_dgtp_c + R_LMPD_103_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c - F_LMPD_105_dttp_c + R_LMPD_105_dttp_c - F_LMPD_106_dttp_c + R_LMPD_106_dttp_c - F_LMPD_108_fad_c + R_LMPD_108_fad_c - F_LMPD_109_fad_c + R_LMPD_109_fad_c + F_LMPD_110_fad_c - R_LMPD_110_fad_c + F_LMPD_111_fad_c - R_LMPD_111_fad_c - F_LMPD_114_GAM_c + R_LMPD_114_GAM_c - F_LMPD_115_GAM_c + R_LMPD_115_GAM_c + F_LMPD_118_gthrd_c - R_LMPD_118_gthrd_c + F_LMPD_119_gthrd_c - R_LMPD_119_gthrd_c - F_LMPD_122_gtp_c + R_LMPD_122_gtp_c - F_LMPD_123_gtp_c + R_LMPD_123_gtp_c - F_LMPD_126_his-L_c + R_LMPD_126_his-L_c - F_LMPD_127_his-L_c + R_LMPD_127_his-L_c -2 F_LMPD_135_malcoa_c +2 R_LMPD_135_malcoa_c -2 F_LMPD_136_malcoa_c +2 R_LMPD_136_malcoa_c - F_LMPD_139_mlthf_c + R_LMPD_139_mlthf_c - F_LMPD_140_mlthf_c + R_LMPD_140_mlthf_c - F_LMPD_141_mococdp_c + R_LMPD_141_mococdp_c - F_LMPD_142_mococdp_c + R_LMPD_142_mococdp_c - F_LMPD_143_mococdp_c + R_LMPD_143_mococdp_c - F_LMPD_144_mococdp_c + R_LMPD_144_mococdp_c -2 F_LMPD_149_mocogdp_c +2 R_LMPD_149_mocogdp_c -2 F_LMPD_150_mocogdp_c +2 R_LMPD_150_mocogdp_c - F_LMPD_153_mql8_c + R_LMPD_153_mql8_c - F_LMPD_161_nad_c + R_LMPD_161_nad_c - F_LMPD_162_nad_c + R_LMPD_162_nad_c -3 F_LMPD_227_q8h2_c +3 R_LMPD_227_q8h2_c -3 F_LMPD_228_q8h2_c +3 R_LMPD_228_q8h2_c -3 F_LMPD_229_q8h2_c +3 R_LMPD_229_q8h2_c -3 F_LMPD_230_q8h2_c +3 R_LMPD_230_q8h2_c -3 F_LMPD_231_q8h2_c +3 R_LMPD_231_q8h2_c -6 F_LMPD_233_q8h2_c +6 R_LMPD_233_q8h2_c -6 F_LMPD_234_q8h2_c +6 R_LMPD_234_q8h2_c + F_LMPD_235_ribflv_c - R_LMPD_235_ribflv_c + F_LMPD_236_ribflv_c - R_LMPD_236_ribflv_c -2 F_LMPD_240_sheme_c +2 R_LMPD_240_sheme_c + F_LMPD_243_thf_c - R_LMPD_243_thf_c + F_LMPD_244_thf_c - R_LMPD_244_thf_c +2 F_LMPD_247_thmpp_c -2 R_LMPD_247_thmpp_c +2 F_LMPD_248_thmpp_c -2 R_LMPD_248_thmpp_c = 0
M_fum_c: - F_ASPO5 + R_ASPO5 + F_ASPT - R_ASPT - F_DHORDfum + R_DHORDfum - F_FRD3 + R_FRD3 - F_FUM + R_FUM -3 F_PPPGO3 +3 R_PPPGO3 + F_SUCDi - R_SUCDi + F_LMPD_1_10fthf_c - R_LMPD_1_10fthf_c + F_LMPD_2_10fthf_c - R_LMPD_2_10fthf_c + F_LMPD_3_10fthf_c - R_LMPD_3_10fthf_c + F_LMPD_4_10fthf_c - R_LMPD_4_10fthf_c + F_LMPD_22_5mthf_c - R_LMPD_22_5mthf_c + F_LMPD_23_5mthf_c - R_LMPD_23_5mthf_c +2 F_LMPD_24_accoa_c -2 R_LMPD_24_accoa_c +2 F_LMPD_25_accoa_c -2 R_LMPD_25_accoa_c +2 F_LMPD_26_adocbl_c -2 R_LMPD_26_adocbl_c +2 F_LMPD_27_adocbl_c -2 R_LMPD_27_adocbl_c +2 F_LMPD_28_adocbl_c -2 R_LMPD_28_adocbl_c +2 F_LMPD_29_adocbl_c -2 R_LMPD_29_adocbl_c +2 F_LMPD_31_amet_c -2 R_LMPD_31_amet_c +2 F_LMPD_32_amet_c -2 R_LMPD_32_amet_c + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c +4 F_LMPD_35_bmocogdp_c -4 R_LMPD_35_bmocogdp_c +4 F_LMPD_36_bmocogdp_c -4 R_LMPD_36_bmocogdp_c +4 F_LMPD_37_bmocogdp_c -4 R_LMPD_37_bmocogdp_c +4 F_LMPD_38_bmocogdp_c -4 R_LMPD_38_bmocogdp_c +2 F_LMPD_39_btn_c -2 R_LMPD_39_btn_c +2 F_LMPD_40_btn_c -2 R_LMPD_40_btn_c +2 F_LMPD_41_btn_c -2 R_LMPD_41_btn_c +2 F_LMPD_42_btn_c -2 R_LMPD_42_btn_c +6 F_LMPD_43_btn_c -6 R_LMPD_43_btn_c +6 F_LMPD_44_btn_c -6 R_LMPD_44_btn_c +6 F_LMPD_45_btn_c -6 R_LMPD_45_btn_c +6 F_LMPD_47_btn_c -6 R_LMPD_47_btn_c +2 F_LMPD_94_datp_c -2 R_LMPD_94_datp_c +2 F_LMPD_95_datp_c -2 R_LMPD_95_datp_c +2 F_LMPD_96_datp_c -2 R_LMPD_96_datp_c +2 F_LMPD_97_datp_c -2 R_LMPD_97_datp_c + F_LMPD_99_dgtp_c - R_LMPD_99_dgtp_c + F_LMPD_100_dgtp_c - R_LMPD_100_dgtp_c + F_LMPD_101_dgtp_c - R_LMPD_101_dgtp_c + F_LMPD_102_dgtp_c - R_LMPD_102_dgtp_c + F_LMPD_103_dgtp_c - R_LMPD_103_dgtp_c + F_LMPD_104_dgtp_c - R_LMPD_104_dgtp_c +3 F_LMPD_108_fad_c -3 R_LMPD_108_fad_c +3 F_LMPD_109_fad_c -3 R_LMPD_109_fad_c +3 F_LMPD_110_fad_c -3 R_LMPD_110_fad_c +3 F_LMPD_111_fad_c -3 R_LMPD_111_fad_c +2 F_LMPD_112_GAM_c -2 R_LMPD_112_GAM_c +2 F_LMPD_113_GAM_c -2 R_LMPD_113_GAM_c +2 F_LMPD_114_GAM_c -2 R_LMPD_114_GAM_c +2 F_LMPD_115_GAM_c -2 R_LMPD_115_GAM_c + F_LMPD_120_gtp_c - R_LMPD_120_gtp_c + F_LMPD_121_gtp_c - R_LMPD_121_gtp_c + F_LMPD_122_gtp_c - R_LMPD_122_gtp_c + F_LMPD_123_gtp_c - R_LMPD_123_gtp_c + F_LMPD_126_his-L_c - R_LMPD_126_his-L_c + F_LMPD_127_his-L_c - R_LMPD_127_his-L_c +2 F_LMPD_135_malcoa_c -2 R_LMPD_135_malcoa_c +2 F_LMPD_136_malcoa_c -2 R_LMPD_136_malcoa_c + F_LMPD_139_mlthf_c - R_LMPD_139_mlthf_c + F_LMPD_140_mlthf_c - R_LMPD_140_mlthf_c + F_LMPD_141_mococdp_c - R_LMPD_141_mococdp_c + F_LMPD_142_mococdp_c - R_LMPD_142_mococdp_c + F_LMPD_143_mococdp_c - R_LMPD_143_mococdp_c + F_LMPD_144_mococdp_c - R_LMPD_144_mococdp_c + F_LMPD_145_mococdp_c - R_LMPD_145_mococdp_c + F_LMPD_146_mococdp_c - R_LMPD_146_mococdp_c + F_LMPD_147_mococdp_c - R_LMPD_147_mococdp_c + F_LMPD_148_mococdp_c - R_LMPD_148_mococdp_c +2 F_LMPD_149_mocogdp_c -2 R_LMPD_149_mocogdp_c +2 F_LMPD_150_mocogdp_c -2 R_LMPD_150_mocogdp_c +2 F_LMPD_151_mocogdp_c -2 R_LMPD_151_mocogdp_c +2 F_LMPD_152_mocogdp_c -2 R_LMPD_152_mocogdp_c +2 F_LMPD_159_nad_c -2 R_LMPD_159_nad_c +2 F_LMPD_160_nad_c -2 R_LMPD_160_nad_c +2 F_LMPD_161_nad_c -2 R_LMPD_161_nad_c +2 F_LMPD_162_nad_c -2 R_LMPD_162_nad_c + F_LMPD_235_ribflv_c - R_LMPD_235_ribflv_c + F_LMPD_236_ribflv_c - R_LMPD_236_ribflv_c + F_LMPD_237_ribflv_c - R_LMPD_237_ribflv_c + F_LMPD_238_ribflv_c - R_LMPD_238_ribflv_c + F_LMPD_243_thf_c - R_LMPD_243_thf_c + F_LMPD_244_thf_c - R_LMPD_244_thf_c + F_LMPD_245_thf_c - R_LMPD_245_thf_c + F_LMPD_246_thf_c - R_LMPD_246_thf_c = 0
M_g1p_c: F_GLCP - R_GLCP - F_GLGC + R_GLGC - F_PGMT + R_PGMT -5 F_LMPD_79_colipa_e +5 R_LMPD_79_colipa_e -5 F_LMPD_80_colipa_e +5 R_LMPD_80_colipa_e -5 F_LMPD_81_colipa_e +5 R_LMPD_81_colipa_e -5 F_LMPD_82_colipa_e +5 R_LMPD_82_colipa_e -5 F_LMPD_83_colipa_e +5 R_LMPD_83_colipa_e -5 F_LMPD_84_colipa_e +5 R_LMPD_84_colipa_e -5 F_LMPD_85_colipa_e +5 R_LMPD_85_colipa_e -5 F_LMPD_86_colipa_e +5 R_LMPD_86_colipa_e -5 F_LMPD_87_colipa_e +5 R_LMPD_87_colipa_e -5 F_LMPD_88_colipa_e +5 R_LMPD_88_colipa_e -5 F_LMPD_89_colipa_e +5 R_LMPD_89_colipa_e -5 F_LMPD_90_colipa_e +5 R_LMPD_90_colipa_e = 0
M_g3p_c: F_EDA - R_EDA + F_F6PA - R_F6PA + F_FBA - R_FBA - F_GAPD + R_GAPD - F_TALA + R_TALA + F_TKT1 - R_TKT1 + F_TKT2 - R_TKT2 + F_TPI - R_TPI -8 F_LMPD_5_2dmmql8_c +8 R_LMPD_5_2dmmql8_c -3 F_LMPD_124_hemeO_c +3 R_LMPD_124_hemeO_c -3 F_LMPD_125_hemeO_c +3 R_LMPD_125_hemeO_c -8 F_LMPD_153_mql8_c +8 R_LMPD_153_mql8_c - F_LMPD_222_pydx5p_c + R_LMPD_222_pydx5p_c - F_LMPD_223_pydx5p_c + R_LMPD_223_pydx5p_c -8 F_LMPD_224_q8h2_c +8 R_LMPD_224_q8h2_c -8 F_LMPD_225_q8h2_c +8 R_LMPD_225_q8h2_c -8 F_LMPD_226_q8h2_c +8 R_LMPD_226_q8h2_c -8 F_LMPD_227_q8h2_c +8 R_LMPD_227_q8h2_c -8 F_LMPD_228_q8h2_c +8 R_LMPD_228_q8h2_c -8 F_LMPD_229_q8h2_c +8 R_LMPD_229_q8h2_c -8 F_LMPD_230_q8h2_c +8 R_LMPD_230_q8h2_c -8 F_LMPD_231_q8h2_c +8 R_LMPD_231_q8h2_c -16 F_LMPD_232_q8h2_c +16 R_LMPD_232_q8h2_c -16 F_LMPD_233_q8h2_c +16 R_LMPD_233_q8h2_c -16 F_LMPD_234_q8h2_c +16 R_LMPD_234_q8h2_c - F_LMPD_247_thmpp_c + R_LMPD_247_thmpp_c - F_LMPD_248_thmpp_c + R_LMPD_248_thmpp_c + F_LMPD_250_trp-L_c - R_LMPD_250_trp-L_c -11 F_LMPD_252_udcpdp_c +11 R_LMPD_252_udcpdp_c = 0
M_g6p_c: - F_G6PDH2r + R_G6PDH2r - F_G6PP + R_G6PP + F_GLCptspp - R_GLCptspp + F_HEX1 - R_HEX1 - F_PGI + R_PGI + F_PGMT - R_PGMT = 0
M_glc-D_c: F_G6PP - R_G6PP + F_GLCabcpp - R_GLCabcpp + F_GLCt2pp - R_GLCt2pp - F_HEX1 + R_HEX1 = 0
M_gln-L_c: -0.25571 F_Ec_biomass_iJO1366_WT_53p95M +0.25571 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_116_gln-L_c - R_LMPD_116_gln-L_c = 0
M_glu-L_c: -0.25571 F_Ec_biomass_iJO1366_WT_53p95M +0.25571 R_Ec_biomass_iJO1366_WT_53p95M + F_ASPTA - R_ASPTA - F_GLUDy + R_GLUDy + F_GLUabcpp - R_GLUabcpp + F_GLUt2rpp - R_GLUt2rpp + F_GLUt4pp - R_GLUt4pp -3 F_LMPD_1_10fthf_c +3 R_LMPD_1_10fthf_c -3 F_LMPD_2_10fthf_c +3 R_LMPD_2_10fthf_c -4 F_LMPD_3_10fthf_c +4 R_LMPD_3_10fthf_c -4 F_LMPD_4_10fthf_c +4 R_LMPD_4_10fthf_c -2 F_LMPD_6_2fe2s_c +2 R_LMPD_6_2fe2s_c -2 F_LMPD_7_2fe2s_c +2 R_LMPD_7_2fe2s_c -2 F_LMPD_8_2fe2s_c +2 R_LMPD_8_2fe2s_c -2 F_LMPD_9_2fe2s_c +2 R_LMPD_9_2fe2s_c -2 F_LMPD_10_2fe2s_c +2 R_LMPD_10_2fe2s_c -2 F_LMPD_11_2fe2s_c +2 R_LMPD_11_2fe2s_c -2 F_LMPD_12_2fe2s_c +2 R_LMPD_12_2fe2s_c -2 F_LMPD_13_2fe2s_c +2 R_LMPD_13_2fe2s_c -4 F_LMPD_14_4fe4s_c +4 R_LMPD_14_4fe4s_c -4 F_LMPD_15_4fe4s_c +4 R_LMPD_15_4fe4s_c -4 F_LMPD_16_4fe4s_c +4 R_LMPD_16_4fe4s_c -4 F_LMPD_17_4fe4s_c +4 R_LMPD_17_4fe4s_c -4 F_LMPD_18_4fe4s_c +4 R_LMPD_18_4fe4s_c -4 F_LMPD_19_4fe4s_c +4 R_LMPD_19_4fe4s_c -4 F_LMPD_20_4fe4s_c +4 R_LMPD_20_4fe4s_c -4 F_LMPD_21_4fe4s_c +4 R_LMPD_21_4fe4s_c -2 F_LMPD_22_5mthf_c +2 R_LMPD_22_5mthf_c -2 F_LMPD_23_5mthf_c +2 R_LMPD_23_5mthf_c -2 F_LMPD_24_accoa_c +2 R_LMPD_24_accoa_c -2 F_LMPD_25_accoa_c +2 R_LMPD_25_accoa_c -2 F_LMPD_26_adocbl_c +2 R_LMPD_26_adocbl_c -2 F_LMPD_27_adocbl_c +2 R_LMPD_27_adocbl_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c - F_LMPD_30_ala-L_c + R_LMPD_30_ala-L_c -2 F_LMPD_31_amet_c +2 R_LMPD_31_amet_c -2 F_LMPD_32_amet_c +2 R_LMPD_32_amet_c -2 F_LMPD_33_arg-L_c +2 R_LMPD_33_arg-L_c -12 F_LMPD_35_bmocogdp_c +12 R_LMPD_35_bmocogdp_c -12 F_LMPD_36_bmocogdp_c +12 R_LMPD_36_bmocogdp_c -8 F_LMPD_37_bmocogdp_c +8 R_LMPD_37_bmocogdp_c -8 F_LMPD_38_bmocogdp_c +8 R_LMPD_38_bmocogdp_c -4 F_LMPD_39_btn_c +4 R_LMPD_39_btn_c -4 F_LMPD_40_btn_c +4 R_LMPD_40_btn_c -4 F_LMPD_41_btn_c +4 R_LMPD_41_btn_c -4 F_LMPD_42_btn_c +4 R_LMPD_42_btn_c -6 F_LMPD_43_btn_c +6 R_LMPD_43_btn_c -6 F_LMPD_44_btn_c +6 R_LMPD_44_btn_c -6 F_LMPD_45_btn_c +6 R_LMPD_45_btn_c -6 F_LMPD_47_btn_c +6 R_LMPD_47_btn_c -2 F_LMPD_48_ca2_c +2 R_LMPD_48_ca2_c -2 F_LMPD_49_ca2_c +2 R_LMPD_49_ca2_c - F_LMPD_50_ca2_c + R_LMPD_50_ca2_c - F_LMPD_92_cys-L_c + R_LMPD_92_cys-L_c - F_LMPD_93_cys-L_c + R_LMPD_93_cys-L_c -2 F_LMPD_94_datp_c +2 R_LMPD_94_datp_c -2 F_LMPD_95_datp_c +2 R_LMPD_95_datp_c - F_LMPD_96_datp_c + R_LMPD_96_datp_c - F_LMPD_97_datp_c + R_LMPD_97_datp_c -2 F_LMPD_99_dgtp_c +2 R_LMPD_99_dgtp_c -2 F_LMPD_100_dgtp_c +2 R_LMPD_100_dgtp_c -2 F_LMPD_101_dgtp_c +2 R_LMPD_101_dgtp_c - F_LMPD_102_dgtp_c + R_LMPD_102_dgtp_c - F_LMPD_103_dgtp_c + R_LMPD_103_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c -3 F_LMPD_107_enter_c +3 R_LMPD_107_enter_c -2 F_LMPD_108_fad_c +2 R_LMPD_108_fad_c -2 F_LMPD_109_fad_c +2 R_LMPD_109_fad_c -4 F_LMPD_110_fad_c +4 R_LMPD_110_fad_c -4 F_LMPD_111_fad_c +4 R_LMPD_111_fad_c -2 F_LMPD_112_GAM_c +2 R_LMPD_112_GAM_c -2 F_LMPD_113_GAM_c +2 R_LMPD_113_GAM_c - F_LMPD_114_GAM_c + R_LMPD_114_GAM_c - F_LMPD_115_GAM_c + R_LMPD_115_GAM_c - F_LMPD_116_gln-L_c + R_LMPD_116_gln-L_c -3 F_LMPD_118_gthrd_c +3 R_LMPD_118_gthrd_c -3 F_LMPD_119_gthrd_c +3 R_LMPD_119_gthrd_c -2 F_LMPD_120_gtp_c +2 R_LMPD_120_gtp_c -2 F_LMPD_121_gtp_c +2 R_LMPD_121_gtp_c - F_LMPD_122_gtp_c + R_LMPD_122_gtp_c - F_LMPD_123_gtp_c + R_LMPD_123_gtp_c -8 F_LMPD_124_hemeO_c +8 R_LMPD_124_hemeO_c -8 F_LMPD_125_hemeO_c +8 R_LMPD_125_hemeO_c - F_LMPD_126_his-L_c + R_LMPD_126_his-L_c - F_LMPD_127_his-L_c + R_LMPD_127_his-L_c - F_LMPD_128_ile-L_c + R_LMPD_128_ile-L_c - F_LMPD_129_leu-L_c + R_LMPD_129_leu-L_c -4 F_LMPD_130_lipopb_c +4 R_LMPD_130_lipopb_c -4 F_LMPD_131_lipopb_c +4 R_LMPD_131_lipopb_c -4 F_LMPD_132_lipopb_c +4 R_LMPD_132_lipopb_c -4 F_LMPD_133_lipopb_c +4 R_LMPD_133_lipopb_c - F_LMPD_134_lys-L_c + R_LMPD_134_lys-L_c -2 F_LMPD_135_malcoa_c +2 R_LMPD_135_malcoa_c -2 F_LMPD_136_malcoa_c +2 R_LMPD_136_malcoa_c -2 F_LMPD_137_met-L_c +2 R_LMPD_137_met-L_c -2 F_LMPD_138_met-L_c +2 R_LMPD_138_met-L_c -2 F_LMPD_139_mlthf_c +2 R_LMPD_139_mlthf_c -2 F_LMPD_140_mlthf_c +2 R_LMPD_140_mlthf_c -3 F_LMPD_141_mococdp_c +3 R_LMPD_141_mococdp_c -3 F_LMPD_142_mococdp_c +3 R_LMPD_142_mococdp_c -3 F_LMPD_143_mococdp_c +3 R_LMPD_143_mococdp_c -3 F_LMPD_144_mococdp_c +3 R_LMPD_144_mococdp_c -4 F_LMPD_145_mococdp_c +4 R_LMPD_145_mococdp_c -4 F_LMPD_146_mococdp_c +4 R_LMPD_146_mococdp_c -4 F_LMPD_147_mococdp_c +4 R_LMPD_147_mococdp_c -4 F_LMPD_148_mococdp_c +4 R_LMPD_148_mococdp_c -4 F_LMPD_149_mocogdp_c +4 R_LMPD_149_mocogdp_c -4 F_LMPD_150_mocogdp_c +4 R_LMPD_150_mocogdp_c -6 F_LMPD_151_mocogdp_c +6 R_LMPD_151_mocogdp_c -6 F_LMPD_152_mocogdp_c +6 R_LMPD_152_mocogdp_c -6 F_LMPD_154_murein3p3p_p +6 R_LMPD_154_murein3p3p_p -10 F_LMPD_155_murein3px4p_p +10 R_LMPD_155_murein3px4p_p -10 F_LMPD_156_murein4p4p_p +10 R_LMPD_156_murein4p4p_p -10 F_LMPD_157_murein4px4p_p +10 R_LMPD_157_murein4px4p_p -15 F_LMPD_158_murein4px4px4p_p +15 R_LMPD_158_murein4px4px4p_p -2 F_LMPD_159_nad_c +2 R_LMPD_159_nad_c -2 F_LMPD_160_nad_c +2 R_LMPD_160_nad_c - F_LMPD_161_nad_c + R_LMPD_161_nad_c - F_LMPD_162_nad_c + R_LMPD_162_nad_c - F_LMPD_163_pe160_c + R_LMPD_163_pe160_c - F_LMPD_164_pe160_c + R_LMPD_164_pe160_c - F_LMPD_165_pe160_c + R_LMPD_165_pe160_c - F_LMPD_166_pe160_c + R_LMPD_166_pe160_c - F_LMPD_167_pe160_c + R_LMPD_167_pe160_c - F_LMPD_168_pe160_c + R_LMPD_168_pe160_c - F_LMPD_169_pe160_c + R_LMPD_169_pe160_c - F_LMPD_170_pe160_c + R_LMPD_170_pe160_c - F_LMPD_171_pe160_c + R_LMPD_171_pe160_c - F_LMPD_172_pe160_c + R_LMPD_172_pe160_c - F_LMPD_173_pe161_c + R_LMPD_173_pe161_c - F_LMPD_174_pe161_c + R_LMPD_174_pe161_c - F_LMPD_175_pe161_c + R_LMPD_175_pe161_c - F_LMPD_176_pe161_c + R_LMPD_176_pe161_c - F_LMPD_177_pe161_c + R_LMPD_177_pe161_c - F_LMPD_178_pe161_c + R_LMPD_178_pe161_c - F_LMPD_179_pe161_c + R_LMPD_179_pe161_c - F_LMPD_180_pe181_c + R_LMPD_180_pe181_c - F_LMPD_181_pe181_c + R_LMPD_181_pe181_c - F_LMPD_182_pe181_c + R_LMPD_182_pe181_c - F_LMPD_183_pe181_c + R_LMPD_183_pe181_c - F_LMPD_184_pe181_c + R_LMPD_184_pe181_c - F_LMPD_185_pe181_c + R_LMPD_185_pe181_c - F_LMPD_186_pe181_c + R_LMPD_186_pe181_c - F_LMPD_187_pe181_c + R_LMPD_187_pe181_c - F_LMPD_188_pe181_c + R_LMPD_188_pe181_c - F_LMPD_189_pe181_c + R_LMPD_189_pe181_c - F_LMPD_217_phe-L_c + R_LMPD_217_phe-L_c -8 F_LMPD_218_pheme_c +8 R_LMPD_218_pheme_c -8 F_LMPD_219_pheme_c +8 R_LMPD_219_pheme_c - F_LMPD_220_pro-L_c + R_LMPD_220_pro-L_c -2 F_LMPD_221_ptrc_c +2 R_LMPD_221_ptrc_c - F_LMPD_222_pydx5p_c + R_LMPD_222_pydx5p_c - F_LMPD_223_pydx5p_c + R_LMPD_223_pydx5p_c -3 F_LMPD_224_q8h2_c +3 R_LMPD_224_q8h2_c -3 F_LMPD_225_q8h2_c +3 R_LMPD_225_q8h2_c -3 F_LMPD_226_q8h2_c +3 R_LMPD_226_q8h2_c -6 F_LMPD_232_q8h2_c +6 R_LMPD_232_q8h2_c -2 F_LMPD_235_ribflv_c +2 R_LMPD_235_ribflv_c -2 F_LMPD_236_ribflv_c +2 R_LMPD_236_ribflv_c - F_LMPD_237_ribflv_c + R_LMPD_237_ribflv_c - F_LMPD_238_ribflv_c + R_LMPD_238_ribflv_c - F_LMPD_239_ser-L_c + R_LMPD_239_ser-L_c -8 F_LMPD_240_sheme_c +8 R_LMPD_240_sheme_c -4 F_LMPD_241_spmd_c +4 R_LMPD_241_spmd_c -4 F_LMPD_242_spmd_c +4 R_LMPD_242_spmd_c -3 F_LMPD_243_thf_c +3 R_LMPD_243_thf_c -3 F_LMPD_244_thf_c +3 R_LMPD_244_thf_c -2 F_LMPD_245_thf_c +2 R_LMPD_245_thf_c -2 F_LMPD_246_thf_c +2 R_LMPD_246_thf_c -3 F_LMPD_247_thmpp_c +3 R_LMPD_247_thmpp_c -3 F_LMPD_248_thmpp_c +3 R_LMPD_248_thmpp_c - F_LMPD_251_tyr-L_c + R_LMPD_251_tyr-L_c - F_LMPD_254_val-L_c + R_LMPD_254_val-L_c = 0
M_glx_c: -2 F_GLXCL +2 R_GLXCL - F_GLYCLTDx + R_GLYCLTDx - F_GLYCLTDy + R_GLYCLTDy + F_GLYCTO2 - R_GLYCTO2 + F_GLYCTO3 - R_GLYCTO3 + F_GLYCTO4 - R_GLYCTO4 + F_ICL - R_ICL - F_MALS + R_MALS = 0
M_gly_c: -0.5953 F_Ec_biomass_iJO1366_WT_53p95M +0.5953 R_Ec_biomass_iJO1366_WT_53p95M + F_GLYt2rpp - R_GLYt2rpp + F_LMPD_1_10fthf_c - R_LMPD_1_10fthf_c + F_LMPD_2_10fthf_c - R_LMPD_2_10fthf_c +2 F_LMPD_3_10fthf_c -2 R_LMPD_3_10fthf_c +2 F_LMPD_4_10fthf_c -2 R_LMPD_4_10fthf_c + F_LMPD_26_adocbl_c - R_LMPD_26_adocbl_c + F_LMPD_27_adocbl_c - R_LMPD_27_adocbl_c +4 F_LMPD_35_bmocogdp_c -4 R_LMPD_35_bmocogdp_c +4 F_LMPD_36_bmocogdp_c -4 R_LMPD_36_bmocogdp_c + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c + F_LMPD_94_datp_c - R_LMPD_94_datp_c + F_LMPD_95_datp_c - R_LMPD_95_datp_c + F_LMPD_99_dgtp_c - R_LMPD_99_dgtp_c + F_LMPD_100_dgtp_c - R_LMPD_100_dgtp_c + F_LMPD_101_dgtp_c - R_LMPD_101_dgtp_c +2 F_LMPD_110_fad_c -2 R_LMPD_110_fad_c +2 F_LMPD_111_fad_c -2 R_LMPD_111_fad_c + F_LMPD_112_GAM_c - R_LMPD_112_GAM_c + F_LMPD_113_GAM_c - R_LMPD_113_GAM_c + F_LMPD_117_gly_c - R_LMPD_117_gly_c + F_LMPD_120_gtp_c - R_LMPD_120_gtp_c + F_LMPD_121_gtp_c - R_LMPD_121_gtp_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c + F_LMPD_145_mococdp_c - R_LMPD_145_mococdp_c + F_LMPD_146_mococdp_c - R_LMPD_146_mococdp_c + F_LMPD_147_mococdp_c - R_LMPD_147_mococdp_c + F_LMPD_148_mococdp_c - R_LMPD_148_mococdp_c +2 F_LMPD_151_mocogdp_c -2 R_LMPD_151_mocogdp_c +2 F_LMPD_152_mocogdp_c -2 R_LMPD_152_mocogdp_c + F_LMPD_159_nad_c - R_LMPD_159_nad_c + F_LMPD_160_nad_c - R_LMPD_160_nad_c +3 F_LMPD_224_q8h2_c -3 R_LMPD_224_q8h2_c +3 F_LMPD_225_q8h2_c -3 R_LMPD_225_q8h2_c +3 F_LMPD_226_q8h2_c -3 R_LMPD_226_q8h2_c +6 F_LMPD_232_q8h2_c -6 R_LMPD_232_q8h2_c + F_LMPD_235_ribflv_c - R_LMPD_235_ribflv_c + F_LMPD_236_ribflv_c - R_LMPD_236_ribflv_c + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c + F_LMPD_243_thf_c - R_LMPD_243_thf_c + F_LMPD_244_thf_c - R_LMPD_244_thf_c = 0
M_glyc_c: F_G3PT - R_G3PT - F_GLYCDx + R_GLYCDx - F_GLYCtpp + R_GLYCtpp - F_GLYK + R_GLYK = 0
M_glyc-R_c: F_GLYCAt2rpp - R_GLYCAt2rpp - F_GLYCK + R_GLYCK - F_GLYCK2 + R_GLYCK2 + F_TRSARr - R_TRSARr = 0
M_glyc3p_c: - F_G3PD2 + R_G3PD2 - F_G3PD5 + R_G3PD5 - F_G3PD6 + R_G3PD6 - F_G3PD7 + R_G3PD7 - F_G3PT + R_G3PT + F_GLYK - R_GLYK -4 F_LMPD_52_clpn160_p +4 R_LMPD_52_clpn160_p -4 F_LMPD_53_clpn160_p +4 R_LMPD_53_clpn160_p -4 F_LMPD_54_clpn160_p +4 R_LMPD_54_clpn160_p -4 F_LMPD_55_clpn160_p +4 R_LMPD_55_clpn160_p -4 F_LMPD_56_clpn160_p +4 R_LMPD_56_clpn160_p -4 F_LMPD_57_clpn160_p +4 R_LMPD_57_clpn160_p -4 F_LMPD_58_clpn160_p +4 R_LMPD_58_clpn160_p -4 F_LMPD_59_clpn160_p +4 R_LMPD_59_clpn160_p -4 F_LMPD_60_clpn160_p +4 R_LMPD_60_clpn160_p -4 F_LMPD_61_clpn160_p +4 R_LMPD_61_clpn160_p -4 F_LMPD_62_clpn161_p +4 R_LMPD_62_clpn161_p -4 F_LMPD_63_clpn161_p +4 R_LMPD_63_clpn161_p -4 F_LMPD_64_clpn161_p +4 R_LMPD_64_clpn161_p -4 F_LMPD_65_clpn161_p +4 R_LMPD_65_clpn161_p -4 F_LMPD_66_clpn161_p +4 R_LMPD_66_clpn161_p -4 F_LMPD_67_clpn161_p +4 R_LMPD_67_clpn161_p -4 F_LMPD_68_clpn161_p +4 R_LMPD_68_clpn161_p -4 F_LMPD_69_clpn181_p +4 R_LMPD_69_clpn181_p -4 F_LMPD_70_clpn181_p +4 R_LMPD_70_clpn181_p -4 F_LMPD_71_clpn181_p +4 R_LMPD_71_clpn181_p -4 F_LMPD_72_clpn181_p +4 R_LMPD_72_clpn181_p -4 F_LMPD_73_clpn181_p +4 R_LMPD_73_clpn181_p -4 F_LMPD_74_clpn181_p +4 R_LMPD_74_clpn181_p -4 F_LMPD_75_clpn181_p +4 R_LMPD_75_clpn181_p -4 F_LMPD_76_clpn181_p +4 R_LMPD_76_clpn181_p -4 F_LMPD_77_clpn181_p +4 R_LMPD_77_clpn181_p -4 F_LMPD_78_clpn181_p +4 R_LMPD_78_clpn181_p - F_LMPD_163_pe160_c + R_LMPD_163_pe160_c - F_LMPD_164_pe160_c + R_LMPD_164_pe160_c - F_LMPD_165_pe160_c + R_LMPD_165_pe160_c - F_LMPD_166_pe160_c + R_LMPD_166_pe160_c - F_LMPD_167_pe160_c + R_LMPD_167_pe160_c - F_LMPD_168_pe160_c + R_LMPD_168_pe160_c - F_LMPD_169_pe160_c + R_LMPD_169_pe160_c - F_LMPD_170_pe160_c + R_LMPD_170_pe160_c - F_LMPD_171_pe160_c + R_LMPD_171_pe160_c - F_LMPD_172_pe160_c + R_LMPD_172_pe160_c - F_LMPD_173_pe161_c + R_LMPD_173_pe161_c - F_LMPD_174_pe161_c + R_LMPD_174_pe161_c - F_LMPD_175_pe161_c + R_LMPD_175_pe161_c - F_LMPD_176_pe161_c + R_LMPD_176_pe161_c - F_LMPD_177_pe161_c + R_LMPD_177_pe161_c - F_LMPD_178_pe161_c + R_LMPD_178_pe161_c - F_LMPD_179_pe161_c + R_LMPD_179_pe161_c - F_LMPD_180_pe181_c + R_LMPD_180_pe181_c - F_LMPD_181_pe181_c + R_LMPD_181_pe181_c - F_LMPD_182_pe181_c + R_LMPD_182_pe181_c - F_LMPD_183_pe181_c + R_LMPD_183_pe181_c - F_LMPD_184_pe181_c + R_LMPD_184_pe181_c - F_LMPD_185_pe181_c + R_LMPD_185_pe181_c - F_LMPD_186_pe181_c + R_LMPD_186_pe181_c - F_LMPD_187_pe181_c + R_LMPD_187_pe181_c - F_LMPD_188_pe181_c + R_LMPD_188_pe181_c - F_LMPD_189_pe181_c + R_LMPD_189_pe181_c -2 F_LMPD_190_pg160_c +2 R_LMPD_190_pg160_c -2 F_LMPD_191_pg160_c +2 R_LMPD_191_pg160_c -2 F_LMPD_192_pg160_c +2 R_LMPD_192_pg160_c -2 F_LMPD_193_pg160_c +2 R_LMPD_193_pg160_c -2 F_LMPD_194_pg160_c +2 R_LMPD_194_pg160_c -2 F_LMPD_195_pg160_c +2 R_LMPD_195_pg160_c -2 F_LMPD_196_pg160_c +2 R_LMPD_196_pg160_c -2 F_LMPD_197_pg160_c +2 R_LMPD_197_pg160_c -2 F_LMPD_198_pg160_c +2 R_LMPD_198_pg160_c -2 F_LMPD_199_pg160_c +2 R_LMPD_199_pg160_c -2 F_LMPD_200_pg161_c +2 R_LMPD_200_pg161_c -2 F_LMPD_201_pg161_c +2 R_LMPD_201_pg161_c -2 F_LMPD_202_pg161_c +2 R_LMPD_202_pg161_c -2 F_LMPD_203_pg161_c +2 R_LMPD_203_pg161_c -2 F_LMPD_204_pg161_c +2 R_LMPD_204_pg161_c -2 F_LMPD_205_pg161_c +2 R_LMPD_205_pg161_c -2 F_LMPD_206_pg161_c +2 R_LMPD_206_pg161_c -2 F_LMPD_207_pg181_c +2 R_LMPD_207_pg181_c -2 F_LMPD_208_pg181_c +2 R_LMPD_208_pg181_c -2 F_LMPD_209_pg181_c +2 R_LMPD_209_pg181_c -2 F_LMPD_210_pg181_c +2 R_LMPD_210_pg181_c -2 F_LMPD_211_pg181_c +2 R_LMPD_211_pg181_c -2 F_LMPD_212_pg181_c +2 R_LMPD_212_pg181_c -2 F_LMPD_213_pg181_c +2 R_LMPD_213_pg181_c -2 F_LMPD_214_pg181_c +2 R_LMPD_214_pg181_c -2 F_LMPD_215_pg181_c +2 R_LMPD_215_pg181_c -2 F_LMPD_216_pg181_c +2 R_LMPD_216_pg181_c = 0
M_glyclt_c: F_GLYCLTDx - R_GLYCLTDx + F_GLYCLTDy - R_GLYCLTDy + F_GLYCLTt2rpp - R_GLYCLTt2rpp - F_GLYCTO2 + R_GLYCTO2 - F_GLYCTO3 + R_GLYCTO3 - F_GLYCTO4 + R_GLYCTO4 + F_LMPD_1_10fthf_c - R_LMPD_1_10fthf_c + F_LMPD_2_10fthf_c - R_LMPD_2_10fthf_c + F_LMPD_3_10fthf_c - R_LMPD_3_10fthf_c + F_LMPD_4_10fthf_c - R_LMPD_4_10fthf_c + F_LMPD_22_5mthf_c - R_LMPD_22_5mthf_c + F_LMPD_23_5mthf_c - R_LMPD_23_5mthf_c + F_LMPD_139_mlthf_c - R_LMPD_139_mlthf_c + F_LMPD_140_mlthf_c - R_LMPD_140_mlthf_c + F_LMPD_243_thf_c - R_LMPD_243_thf_c + F_LMPD_244_thf_c - R_LMPD_244_thf_c + F_LMPD_245_thf_c - R_LMPD_245_thf_c + F_LMPD_246_thf_c - R_LMPD_246_thf_c = 0
M_glycogen_c: -0.15419 F_Ec_biomass_iJO1366_WT_53p95M +0.15419 R_Ec_biomass_iJO1366_WT_53p95M - F_GLBRAN2 + R_GLBRAN2 - F_GLCP + R_GLCP + F_GLCS1 - R_GLCS1 + F_GLDBRAN2 - R_GLDBRAN2 = 0
M_glycogenn1_c: - F_DM_glycogenn1_c + R_DM_glycogenn1_c + F_GLCP - R_GLCP - F_GLCS1 + R_GLCS1 = 0
M_gthrd_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M - F_GTHRDabc2pp + R_GTHRDabc2pp + F_LMPD_48_ca2_c - R_LMPD_48_ca2_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c + F_LMPD_50_ca2_c - R_LMPD_50_ca2_c + F_LMPD_118_gthrd_c - R_LMPD_118_gthrd_c + F_LMPD_119_gthrd_c - R_LMPD_119_gthrd_c = 0
M_gtp_c: -0.20912 F_Ec_biomass_iJO1366_WT_53p95M +0.20912 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_120_gtp_c - R_LMPD_120_gtp_c + F_LMPD_121_gtp_c - R_LMPD_121_gtp_c + F_LMPD_122_gtp_c - R_LMPD_122_gtp_c + F_LMPD_123_gtp_c - R_LMPD_123_gtp_c = 0
M_h_c: +53.95 F_Ec_biomass_iJO1366_WT_53p95M -53.95 R_Ec_biomass_iJO1366_WT_53p95M + F_ACALD - R_ACALD - F_ACKr + R_ACKr - F_ACS + R_ACS + F_ACt2rpp - R_ACt2rpp + F_AKGt2rpp - R_AKGt2rpp + F_ALAt2rpp - R_ALAt2rpp + F_ALCD2x - R_ALCD2x +2 F_ALDD2x -2 R_ALDD2x +2 F_ALDD2y -2 R_ALDD2y + F_ARGt3pp - R_ARGt3pp + F_ASNt2rpp - R_ASNt2rpp + F_ASPO3 - R_ASPO3 + F_ASPO4 - R_ASPO4 + F_ASPO5 - R_ASPO5 + F_ASPO6 - R_ASPO6 + F_ASPt2rpp - R_ASPt2rpp + F_ATPM - R_ATPM +3 F_ATPS4rpp -3 R_ATPS4rpp + F_CBL1tonex - R_CBL1tonex + F_CITt3pp - R_CITt3pp - F_CLt3_2pp + R_CLt3_2pp + F_COLIPAabcpp - R_COLIPAabcpp + F_COLIPAabctex - R_COLIPAabctex + F_CS - R_CS + F_CYSabc2pp - R_CYSabc2pp -2 F_CYTBD2pp +2 R_CYTBD2pp -2 F_CYTBDpp +2 R_CYTBDpp -4 F_CYTBO3_4pp +4 R_CYTBO3_4pp + F_D-LACt2pp - R_D-LACt2pp + F_DB4PS - R_DB4PS + F_FBP - R_FBP -2 F_FDH4pp +2 R_FDH4pp -2 F_FDH5pp +2 R_FDH5pp + F_FE2abcpp - R_FE2abcpp + F_FE2t2pp - R_FE2t2pp + F_FE3abcpp - R_FE3abcpp + F_FEENTERtex - R_FEENTERtex + F_FEENTERtpp - R_FEENTERtpp - F_FLDR2 + R_FLDR2 + F_G3PD2 - R_G3PD2 + F_G3PT - R_G3PT + F_G6PDH2r - R_G6PDH2r + F_G6PP - R_G6PP - F_GLCP + R_GLCP + F_GLCabcpp - R_GLCabcpp + F_GLCt2pp - R_GLCt2pp + F_GLUDy - R_GLUDy + F_GLUabcpp - R_GLUabcpp + F_GLUt2rpp - R_GLUt2rpp - F_GLXCL + R_GLXCL + F_GLYCAt2rpp - R_GLYCAt2rpp + F_GLYCDx - R_GLYCDx - F_GLYCLTDx + R_GLYCLTDx - F_GLYCLTDy + R_GLYCLTDy + F_GLYCLTt2rpp - R_GLYCLTt2rpp + F_GLYt2rpp - R_GLYt2rpp + F_GTHRDabc2pp - R_GTHRDabc2pp + F_HISt2rpp - R_HISt2rpp + F_ILEt2rpp - R_ILEt2rpp + F_Kabcpp - R_Kabcpp + F_Kt2pp - R_Kt2pp + F_LDH_D - R_LDH_D + F_LEUt2rpp - R_LEUt2rpp + F_LYSt3pp - R_LYSt3pp + F_MALS - R_MALS + F_MALt3pp - R_MALt3pp + F_MDH - R_MDH -2 F_MG2t3_2pp +2 R_MG2t3_2pp + F_MG2uabcpp - R_MG2uabcpp + F_MNt2pp - R_MNt2pp + F_MOBDabcpp - R_MOBDabcpp - F_NADH10 + R_NADH10 -4 F_NADH16pp +4 R_NADH16pp -4 F_NADH17pp +4 R_NADH17pp -4 F_NADH18pp +4 R_NADH18pp - F_NADH5 + R_NADH5 - F_NADH9 + R_NADH9 - F_NADPHQR2 + R_NADPHQR2 - F_NADPHQR3 + R_NADPHQR3 - F_NADPHQR4 + R_NADPHQR4 + F_NADPPPS - R_NADPPPS +3 F_NAt3_1p5pp -3 R_NAt3_1p5pp +2 F_NAt3_2pp -2 R_NAt3_2pp + F_NAt3pp - R_NAt3pp + F_NI2uabcpp - R_NI2uabcpp + F_NTP1 - R_NTP1 - F_OAADC + R_OAADC + F_PE160abcpp - R_PE160abcpp + F_PE161abcpp - R_PE161abcpp + F_PE181abcpp - R_PE181abcpp + F_PG160abcpp - R_PG160abcpp + F_PG161abcpp - R_PG161abcpp + F_PG181abcpp - R_PG181abcpp - F_PGK + R_PGK + F_PGL - R_PGL + F_PHEMEabcpp - R_PHEMEabcpp + F_PHEt2rpp - R_PHEt2rpp + F_PIt2rpp - R_PIt2rpp + F_PIuabcpp - R_PIuabcpp - F_POR5 + R_POR5 +2 F_PPA -2 R_PPA +2 F_PPC -2 R_PPC - F_PPCK + R_PPCK + F_PPK2r - R_PPK2r - F_PPKr + R_PPKr + F_PPS - R_PPS + F_PROt2rpp - R_PROt2rpp - F_PTAr + R_PTAr + F_PYRt2rpp - R_PYRt2rpp + F_SERt2rpp - R_SERt2rpp + F_SO4t2pp - R_SO4t2pp + F_SPMDt3pp - R_SPMDt3pp +2 F_SUCCt2_2pp -2 R_SUCCt2_2pp +3 F_SUCCt2_3pp -3 R_SUCCt2_3pp + F_SUCCt3pp - R_SUCCt3pp + F_SULabcpp - R_SULabcpp +2 F_THD2pp -2 R_THD2pp + F_THRt2pp - R_THRt2pp + F_THRt2rpp - R_THRt2rpp - F_TRDR + R_TRDR + F_TRPt2rpp - R_TRPt2rpp - F_TRSARr + R_TRSARr + F_TYRt2rpp - R_TYRt2rpp + F_VALt2rpp - R_VALt2rpp + F_ZNabcpp - R_ZNabcpp +28 F_LMPD_1_10fthf_c -28 R_LMPD_1_10fthf_c +28 F_LMPD_2_10fthf_c -28 R_LMPD_2_10fthf_c +31 F_LMPD_3_10fthf_c -31 R_LMPD_3_10fthf_c +31 F_LMPD_4_10fthf_c -31 R_LMPD_4_10fthf_c -24 F_LMPD_5_2dmmql8_c +24 R_LMPD_5_2dmmql8_c + F_LMPD_6_2fe2s_c - R_LMPD_6_2fe2s_c + F_LMPD_7_2fe2s_c - R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c +2 F_LMPD_9_2fe2s_c -2 R_LMPD_9_2fe2s_c + F_LMPD_10_2fe2s_c - R_LMPD_10_2fe2s_c + F_LMPD_11_2fe2s_c - R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +5 F_LMPD_14_4fe4s_c -5 R_LMPD_14_4fe4s_c +5 F_LMPD_15_4fe4s_c -5 R_LMPD_15_4fe4s_c +3 F_LMPD_16_4fe4s_c -3 R_LMPD_16_4fe4s_c +3 F_LMPD_17_4fe4s_c -3 R_LMPD_17_4fe4s_c +5 F_LMPD_18_4fe4s_c -5 R_LMPD_18_4fe4s_c +5 F_LMPD_19_4fe4s_c -5 R_LMPD_19_4fe4s_c +3 F_LMPD_20_4fe4s_c -3 R_LMPD_20_4fe4s_c +3 F_LMPD_21_4fe4s_c -3 R_LMPD_21_4fe4s_c +23 F_LMPD_22_5mthf_c -23 R_LMPD_22_5mthf_c +23 F_LMPD_23_5mthf_c -23 R_LMPD_23_5mthf_c +9 F_LMPD_24_accoa_c -9 R_LMPD_24_accoa_c +9 F_LMPD_25_accoa_c -9 R_LMPD_25_accoa_c +18 F_LMPD_26_adocbl_c -18 R_LMPD_26_adocbl_c +18 F_LMPD_27_adocbl_c -18 R_LMPD_27_adocbl_c +15 F_LMPD_28_adocbl_c -15 R_LMPD_28_adocbl_c +15 F_LMPD_29_adocbl_c -15 R_LMPD_29_adocbl_c +12 F_LMPD_31_amet_c -12 R_LMPD_31_amet_c +12 F_LMPD_32_amet_c -12 R_LMPD_32_amet_c +3 F_LMPD_33_arg-L_c -3 R_LMPD_33_arg-L_c +76 F_LMPD_35_bmocogdp_c -76 R_LMPD_35_bmocogdp_c +76 F_LMPD_36_bmocogdp_c -76 R_LMPD_36_bmocogdp_c +64 F_LMPD_37_bmocogdp_c -64 R_LMPD_37_bmocogdp_c +64 F_LMPD_38_bmocogdp_c -64 R_LMPD_38_bmocogdp_c +17 F_LMPD_39_btn_c -17 R_LMPD_39_btn_c +18 F_LMPD_40_btn_c -18 R_LMPD_40_btn_c +17 F_LMPD_41_btn_c -17 R_LMPD_41_btn_c +18 F_LMPD_42_btn_c -18 R_LMPD_42_btn_c +39 F_LMPD_43_btn_c -39 R_LMPD_43_btn_c +40 F_LMPD_44_btn_c -40 R_LMPD_44_btn_c +39 F_LMPD_45_btn_c -39 R_LMPD_45_btn_c +40 F_LMPD_47_btn_c -40 R_LMPD_47_btn_c +5 F_LMPD_48_ca2_c -5 R_LMPD_48_ca2_c +6 F_LMPD_49_ca2_c -6 R_LMPD_49_ca2_c +3 F_LMPD_50_ca2_c -3 R_LMPD_50_ca2_c +3 F_LMPD_51_chor_c -3 R_LMPD_51_chor_c -26 F_LMPD_52_clpn160_p +26 R_LMPD_52_clpn160_p -26 F_LMPD_53_clpn160_p +26 R_LMPD_53_clpn160_p -26 F_LMPD_54_clpn160_p +26 R_LMPD_54_clpn160_p -26 F_LMPD_55_clpn160_p +26 R_LMPD_55_clpn160_p -26 F_LMPD_56_clpn160_p +26 R_LMPD_56_clpn160_p -26 F_LMPD_57_clpn160_p +26 R_LMPD_57_clpn160_p -26 F_LMPD_58_clpn160_p +26 R_LMPD_58_clpn160_p -26 F_LMPD_59_clpn160_p +26 R_LMPD_59_clpn160_p -54 F_LMPD_60_clpn160_p +54 R_LMPD_60_clpn160_p -54 F_LMPD_61_clpn160_p +54 R_LMPD_61_clpn160_p -22 F_LMPD_62_clpn161_p +22 R_LMPD_62_clpn161_p -22 F_LMPD_63_clpn161_p +22 R_LMPD_63_clpn161_p -22 F_LMPD_64_clpn161_p +22 R_LMPD_64_clpn161_p -22 F_LMPD_65_clpn161_p +22 R_LMPD_65_clpn161_p -22 F_LMPD_66_clpn161_p +22 R_LMPD_66_clpn161_p -22 F_LMPD_67_clpn161_p +22 R_LMPD_67_clpn161_p -22 F_LMPD_68_clpn161_p +22 R_LMPD_68_clpn161_p -26 F_LMPD_69_clpn181_p +26 R_LMPD_69_clpn181_p -26 F_LMPD_70_clpn181_p +26 R_LMPD_70_clpn181_p -26 F_LMPD_71_clpn181_p +26 R_LMPD_71_clpn181_p -26 F_LMPD_72_clpn181_p +26 R_LMPD_72_clpn181_p -26 F_LMPD_73_clpn181_p +26 R_LMPD_73_clpn181_p -26 F_LMPD_74_clpn181_p +26 R_LMPD_74_clpn181_p -26 F_LMPD_75_clpn181_p +26 R_LMPD_75_clpn181_p -26 F_LMPD_76_clpn181_p +26 R_LMPD_76_clpn181_p -58 F_LMPD_77_clpn181_p +58 R_LMPD_77_clpn181_p -58 F_LMPD_78_clpn181_p +58 R_LMPD_78_clpn181_p -16 F_LMPD_79_colipa_e +16 R_LMPD_79_colipa_e -16 F_LMPD_80_colipa_e +16 R_LMPD_80_colipa_e -16 F_LMPD_81_colipa_e +16 R_LMPD_81_colipa_e -16 F_LMPD_82_colipa_e +16 R_LMPD_82_colipa_e -16 F_LMPD_83_colipa_e +16 R_LMPD_83_colipa_e -16 F_LMPD_84_colipa_e +16 R_LMPD_84_colipa_e -16 F_LMPD_85_colipa_e +16 R_LMPD_85_colipa_e -16 F_LMPD_86_colipa_e +16 R_LMPD_86_colipa_e -16 F_LMPD_87_colipa_e +16 R_LMPD_87_colipa_e -16 F_LMPD_88_colipa_e +16 R_LMPD_88_colipa_e -16 F_LMPD_89_colipa_e +16 R_LMPD_89_colipa_e -16 F_LMPD_90_colipa_e +16 R_LMPD_90_colipa_e +4 F_LMPD_91_ctp_c -4 R_LMPD_91_ctp_c +18 F_LMPD_94_datp_c -18 R_LMPD_94_datp_c +18 F_LMPD_95_datp_c -18 R_LMPD_95_datp_c +15 F_LMPD_96_datp_c -15 R_LMPD_96_datp_c +15 F_LMPD_97_datp_c -15 R_LMPD_97_datp_c +4 F_LMPD_98_dctp_c -4 R_LMPD_98_dctp_c +19 F_LMPD_99_dgtp_c -19 R_LMPD_99_dgtp_c +19 F_LMPD_100_dgtp_c -19 R_LMPD_100_dgtp_c +19 F_LMPD_101_dgtp_c -19 R_LMPD_101_dgtp_c +16 F_LMPD_102_dgtp_c -16 R_LMPD_102_dgtp_c +16 F_LMPD_103_dgtp_c -16 R_LMPD_103_dgtp_c +16 F_LMPD_104_dgtp_c -16 R_LMPD_104_dgtp_c - F_LMPD_105_dttp_c + R_LMPD_105_dttp_c - F_LMPD_106_dttp_c + R_LMPD_106_dttp_c +15 F_LMPD_107_enter_c -15 R_LMPD_107_enter_c +33 F_LMPD_108_fad_c -33 R_LMPD_108_fad_c +33 F_LMPD_109_fad_c -33 R_LMPD_109_fad_c +39 F_LMPD_110_fad_c -39 R_LMPD_110_fad_c +39 F_LMPD_111_fad_c -39 R_LMPD_111_fad_c +18 F_LMPD_112_GAM_c -18 R_LMPD_112_GAM_c +18 F_LMPD_113_GAM_c -18 R_LMPD_113_GAM_c +15 F_LMPD_114_GAM_c -15 R_LMPD_114_GAM_c +15 F_LMPD_115_GAM_c -15 R_LMPD_115_GAM_c + F_LMPD_116_gln-L_c - R_LMPD_116_gln-L_c - F_LMPD_117_gly_c + R_LMPD_117_gly_c +6 F_LMPD_118_gthrd_c -6 R_LMPD_118_gthrd_c +6 F_LMPD_119_gthrd_c -6 R_LMPD_119_gthrd_c +19 F_LMPD_120_gtp_c -19 R_LMPD_120_gtp_c +19 F_LMPD_121_gtp_c -19 R_LMPD_121_gtp_c +16 F_LMPD_122_gtp_c -16 R_LMPD_122_gtp_c +16 F_LMPD_123_gtp_c -16 R_LMPD_123_gtp_c -25 F_LMPD_124_hemeO_c +25 R_LMPD_124_hemeO_c -23 F_LMPD_125_hemeO_c +23 R_LMPD_125_hemeO_c +8 F_LMPD_126_his-L_c -8 R_LMPD_126_his-L_c +8 F_LMPD_127_his-L_c -8 R_LMPD_127_his-L_c -3 F_LMPD_128_ile-L_c +3 R_LMPD_128_ile-L_c - F_LMPD_129_leu-L_c + R_LMPD_129_leu-L_c -4 F_LMPD_130_lipopb_c +4 R_LMPD_130_lipopb_c -4 F_LMPD_131_lipopb_c +4 R_LMPD_131_lipopb_c -2 F_LMPD_132_lipopb_c +2 R_LMPD_132_lipopb_c -4 F_LMPD_133_lipopb_c +4 R_LMPD_133_lipopb_c -2 F_LMPD_134_lys-L_c +2 R_LMPD_134_lys-L_c +11 F_LMPD_135_malcoa_c -11 R_LMPD_135_malcoa_c +11 F_LMPD_136_malcoa_c -11 R_LMPD_136_malcoa_c +24 F_LMPD_139_mlthf_c -24 R_LMPD_139_mlthf_c +24 F_LMPD_140_mlthf_c -24 R_LMPD_140_mlthf_c +19 F_LMPD_141_mococdp_c -19 R_LMPD_141_mococdp_c +21 F_LMPD_142_mococdp_c -21 R_LMPD_142_mococdp_c +19 F_LMPD_143_mococdp_c -19 R_LMPD_143_mococdp_c +21 F_LMPD_144_mococdp_c -21 R_LMPD_144_mococdp_c +22 F_LMPD_145_mococdp_c -22 R_LMPD_145_mococdp_c +24 F_LMPD_146_mococdp_c -24 R_LMPD_146_mococdp_c +22 F_LMPD_147_mococdp_c -22 R_LMPD_147_mococdp_c +24 F_LMPD_148_mococdp_c -24 R_LMPD_148_mococdp_c +31 F_LMPD_149_mocogdp_c -31 R_LMPD_149_mocogdp_c +31 F_LMPD_150_mocogdp_c -31 R_LMPD_150_mocogdp_c +37 F_LMPD_151_mocogdp_c -37 R_LMPD_151_mocogdp_c +37 F_LMPD_152_mocogdp_c -37 R_LMPD_152_mocogdp_c -25 F_LMPD_153_mql8_c +25 R_LMPD_153_mql8_c +18 F_LMPD_154_murein3p3p_p -18 R_LMPD_154_murein3p3p_p +18 F_LMPD_155_murein3px4p_p -18 R_LMPD_155_murein3px4p_p +18 F_LMPD_156_murein4p4p_p -18 R_LMPD_156_murein4p4p_p +18 F_LMPD_157_murein4px4p_p -18 R_LMPD_157_murein4px4p_p +27 F_LMPD_158_murein4px4px4p_p -27 R_LMPD_158_murein4px4px4p_p +17 F_LMPD_159_nad_c -17 R_LMPD_159_nad_c +17 F_LMPD_160_nad_c -17 R_LMPD_160_nad_c +14 F_LMPD_161_nad_c -14 R_LMPD_161_nad_c +14 F_LMPD_162_nad_c -14 R_LMPD_162_nad_c -13 F_LMPD_163_pe160_c +13 R_LMPD_163_pe160_c -13 F_LMPD_164_pe160_c +13 R_LMPD_164_pe160_c -13 F_LMPD_165_pe160_c +13 R_LMPD_165_pe160_c -13 F_LMPD_166_pe160_c +13 R_LMPD_166_pe160_c -13 F_LMPD_167_pe160_c +13 R_LMPD_167_pe160_c -13 F_LMPD_168_pe160_c +13 R_LMPD_168_pe160_c -13 F_LMPD_169_pe160_c +13 R_LMPD_169_pe160_c -13 F_LMPD_170_pe160_c +13 R_LMPD_170_pe160_c -27 F_LMPD_171_pe160_c +27 R_LMPD_171_pe160_c -27 F_LMPD_172_pe160_c +27 R_LMPD_172_pe160_c -11 F_LMPD_173_pe161_c +11 R_LMPD_173_pe161_c -11 F_LMPD_174_pe161_c +11 R_LMPD_174_pe161_c -11 F_LMPD_175_pe161_c +11 R_LMPD_175_pe161_c -11 F_LMPD_176_pe161_c +11 R_LMPD_176_pe161_c -11 F_LMPD_177_pe161_c +11 R_LMPD_177_pe161_c -11 F_LMPD_178_pe161_c +11 R_LMPD_178_pe161_c -11 F_LMPD_179_pe161_c +11 R_LMPD_179_pe161_c -13 F_LMPD_180_pe181_c +13 R_LMPD_180_pe181_c -13 F_LMPD_181_pe181_c +13 R_LMPD_181_pe181_c -13 F_LMPD_182_pe181_c +13 R_LMPD_182_pe181_c -13 F_LMPD_183_pe181_c +13 R_LMPD_183_pe181_c -13 F_LMPD_184_pe181_c +13 R_LMPD_184_pe181_c -13 F_LMPD_185_pe181_c +13 R_LMPD_185_pe181_c -13 F_LMPD_186_pe181_c +13 R_LMPD_186_pe181_c -13 F_LMPD_187_pe181_c +13 R_LMPD_187_pe181_c -29 F_LMPD_188_pe181_c +29 R_LMPD_188_pe181_c -29 F_LMPD_189_pe181_c +29 R_LMPD_189_pe181_c -13 F_LMPD_190_pg160_c +13 R_LMPD_190_pg160_c -13 F_LMPD_191_pg160_c +13 R_LMPD_191_pg160_c -13 F_LMPD_192_pg160_c +13 R_LMPD_192_pg160_c -13 F_LMPD_193_pg160_c +13 R_LMPD_193_pg160_c -13 F_LMPD_194_pg160_c +13 R_LMPD_194_pg160_c -13 F_LMPD_195_pg160_c +13 R_LMPD_195_pg160_c -13 F_LMPD_196_pg160_c +13 R_LMPD_196_pg160_c -13 F_LMPD_197_pg160_c +13 R_LMPD_197_pg160_c -27 F_LMPD_198_pg160_c +27 R_LMPD_198_pg160_c -27 F_LMPD_199_pg160_c +27 R_LMPD_199_pg160_c -11 F_LMPD_200_pg161_c +11 R_LMPD_200_pg161_c -11 F_LMPD_201_pg161_c +11 R_LMPD_201_pg161_c -11 F_LMPD_202_pg161_c +11 R_LMPD_202_pg161_c -11 F_LMPD_203_pg161_c +11 R_LMPD_203_pg161_c -11 F_LMPD_204_pg161_c +11 R_LMPD_204_pg161_c -11 F_LMPD_205_pg161_c +11 R_LMPD_205_pg161_c -11 F_LMPD_206_pg161_c +11 R_LMPD_206_pg161_c -13 F_LMPD_207_pg181_c +13 R_LMPD_207_pg181_c -13 F_LMPD_208_pg181_c +13 R_LMPD_208_pg181_c -13 F_LMPD_209_pg181_c +13 R_LMPD_209_pg181_c -13 F_LMPD_210_pg181_c +13 R_LMPD_210_pg181_c -13 F_LMPD_211_pg181_c +13 R_LMPD_211_pg181_c -13 F_LMPD_212_pg181_c +13 R_LMPD_212_pg181_c -13 F_LMPD_213_pg181_c +13 R_LMPD_213_pg181_c -13 F_LMPD_214_pg181_c +13 R_LMPD_214_pg181_c -29 F_LMPD_215_pg181_c +29 R_LMPD_215_pg181_c -29 F_LMPD_216_pg181_c +29 R_LMPD_216_pg181_c +2 F_LMPD_217_phe-L_c -2 R_LMPD_217_phe-L_c -16 F_LMPD_218_pheme_c +16 R_LMPD_218_pheme_c -14 F_LMPD_219_pheme_c +14 R_LMPD_219_pheme_c - F_LMPD_220_pro-L_c + R_LMPD_220_pro-L_c - F_LMPD_221_ptrc_c + R_LMPD_221_ptrc_c +4 F_LMPD_222_pydx5p_c -4 R_LMPD_222_pydx5p_c +5 F_LMPD_223_pydx5p_c -5 R_LMPD_223_pydx5p_c -16 F_LMPD_224_q8h2_c +16 R_LMPD_224_q8h2_c -13 F_LMPD_225_q8h2_c +13 R_LMPD_225_q8h2_c -7 F_LMPD_226_q8h2_c +7 R_LMPD_226_q8h2_c -25 F_LMPD_227_q8h2_c +25 R_LMPD_227_q8h2_c -19 F_LMPD_228_q8h2_c +19 R_LMPD_228_q8h2_c -22 F_LMPD_229_q8h2_c +22 R_LMPD_229_q8h2_c -13 F_LMPD_230_q8h2_c +13 R_LMPD_230_q8h2_c -16 F_LMPD_231_q8h2_c +16 R_LMPD_231_q8h2_c -20 F_LMPD_232_q8h2_c +20 R_LMPD_232_q8h2_c -32 F_LMPD_233_q8h2_c +32 R_LMPD_233_q8h2_c -38 F_LMPD_234_q8h2_c +38 R_LMPD_234_q8h2_c +21 F_LMPD_235_ribflv_c -21 R_LMPD_235_ribflv_c +21 F_LMPD_236_ribflv_c -21 R_LMPD_236_ribflv_c +18 F_LMPD_237_ribflv_c -18 R_LMPD_237_ribflv_c +18 F_LMPD_238_ribflv_c -18 R_LMPD_238_ribflv_c +2 F_LMPD_239_ser-L_c -2 R_LMPD_239_ser-L_c -11 F_LMPD_240_sheme_c +11 R_LMPD_240_sheme_c -2 F_LMPD_241_spmd_c +2 R_LMPD_241_spmd_c -2 F_LMPD_242_spmd_c +2 R_LMPD_242_spmd_c +28 F_LMPD_243_thf_c -28 R_LMPD_243_thf_c +28 F_LMPD_244_thf_c -28 R_LMPD_244_thf_c +25 F_LMPD_245_thf_c -25 R_LMPD_245_thf_c +25 F_LMPD_246_thf_c -25 R_LMPD_246_thf_c +13 F_LMPD_247_thmpp_c -13 R_LMPD_247_thmpp_c +13 F_LMPD_248_thmpp_c -13 R_LMPD_248_thmpp_c - F_LMPD_249_thr-L_c + R_LMPD_249_thr-L_c +4 F_LMPD_250_trp-L_c -4 R_LMPD_250_trp-L_c +3 F_LMPD_251_tyr-L_c -3 R_LMPD_251_tyr-L_c -33 F_LMPD_252_udcpdp_c +33 R_LMPD_252_udcpdp_c + F_LMPD_253_utp_c - R_LMPD_253_utp_c -2 F_LMPD_254_val-L_c +2 R_LMPD_254_val-L_c = 0
M_h2o_c: -48.7529 F_Ec_biomass_iJO1366_WT_53p95M +48.7529 R_Ec_biomass_iJO1366_WT_53p95M + F_ACONTa - R_ACONTa - F_ACONTb + R_ACONTb - F_ALDD2x + R_ALDD2x - F_ALDD2y + R_ALDD2y - F_ATPM + R_ATPM + F_ATPS4rpp - R_ATPS4rpp - F_COLIPAabcpp + R_COLIPAabcpp - F_COLIPAabctex + R_COLIPAabctex - F_CS + R_CS - F_CYSabc2pp + R_CYSabc2pp + F_CYTBD2pp - R_CYTBD2pp + F_CYTBDpp - R_CYTBDpp + F_CYTBO3_4pp - R_CYTBO3_4pp + F_EDD - R_EDD + F_ENO - R_ENO - F_FBP + R_FBP - F_FE2abcpp + R_FE2abcpp - F_FE3abcpp + R_FE3abcpp - F_FUM + R_FUM - F_G3PT + R_G3PT - F_G6PP + R_G6PP - F_GLCabcpp + R_GLCabcpp - F_GLUDy + R_GLUDy - F_GLUabcpp + R_GLUabcpp - F_GTHRDabc2pp + R_GTHRDabc2pp + F_H2Otpp - R_H2Otpp - F_Kabcpp + R_Kabcpp - F_MALS + R_MALS - F_MG2uabcpp + R_MG2uabcpp - F_MOBDabcpp + R_MOBDabcpp - F_NADPPPS + R_NADPPPS - F_NI2uabcpp + R_NI2uabcpp - F_NTP1 + R_NTP1 - F_NTPP6 + R_NTPP6 - F_PE160abcpp + R_PE160abcpp - F_PE161abcpp + R_PE161abcpp - F_PE181abcpp + R_PE181abcpp - F_PG160abcpp + R_PG160abcpp - F_PG161abcpp + R_PG161abcpp - F_PG181abcpp + R_PG181abcpp - F_PGL + R_PGL - F_PHEMEabcpp + R_PHEMEabcpp - F_PIuabcpp + R_PIuabcpp - F_POX + R_POX - F_PPA + R_PPA - F_PPA2 + R_PPA2 - F_PPC + R_PPC +3 F_PPPGO -3 R_PPPGO - F_PPS + R_PPS - F_SULabcpp + R_SULabcpp +2 F_THIORDXi -2 R_THIORDXi - F_ZNabcpp + R_ZNabcpp -10 F_LMPD_1_10fthf_c +10 R_LMPD_1_10fthf_c -10 F_LMPD_2_10fthf_c +10 R_LMPD_2_10fthf_c -11 F_LMPD_3_10fthf_c +11 R_LMPD_3_10fthf_c -11 F_LMPD_4_10fthf_c +11 R_LMPD_4_10fthf_c +17 F_LMPD_5_2dmmql8_c -17 R_LMPD_5_2dmmql8_c - F_LMPD_8_2fe2s_c + R_LMPD_8_2fe2s_c - F_LMPD_9_2fe2s_c + R_LMPD_9_2fe2s_c - F_LMPD_12_2fe2s_c + R_LMPD_12_2fe2s_c - F_LMPD_13_2fe2s_c + R_LMPD_13_2fe2s_c -2 F_LMPD_14_4fe4s_c +2 R_LMPD_14_4fe4s_c -2 F_LMPD_15_4fe4s_c +2 R_LMPD_15_4fe4s_c -2 F_LMPD_18_4fe4s_c +2 R_LMPD_18_4fe4s_c -2 F_LMPD_19_4fe4s_c +2 R_LMPD_19_4fe4s_c -8 F_LMPD_22_5mthf_c +8 R_LMPD_22_5mthf_c -8 F_LMPD_23_5mthf_c +8 R_LMPD_23_5mthf_c -2 F_LMPD_24_accoa_c +2 R_LMPD_24_accoa_c -2 F_LMPD_25_accoa_c +2 R_LMPD_25_accoa_c -5 F_LMPD_26_adocbl_c +5 R_LMPD_26_adocbl_c -5 F_LMPD_27_adocbl_c +5 R_LMPD_27_adocbl_c -4 F_LMPD_28_adocbl_c +4 R_LMPD_28_adocbl_c -4 F_LMPD_29_adocbl_c +4 R_LMPD_29_adocbl_c -4 F_LMPD_31_amet_c +4 R_LMPD_31_amet_c -4 F_LMPD_32_amet_c +4 R_LMPD_32_amet_c - F_LMPD_33_arg-L_c + R_LMPD_33_arg-L_c -25 F_LMPD_35_bmocogdp_c +25 R_LMPD_35_bmocogdp_c -25 F_LMPD_36_bmocogdp_c +25 R_LMPD_36_bmocogdp_c -21 F_LMPD_37_bmocogdp_c +21 R_LMPD_37_bmocogdp_c -21 F_LMPD_38_bmocogdp_c +21 R_LMPD_38_bmocogdp_c -11 F_LMPD_39_btn_c +11 R_LMPD_39_btn_c -12 F_LMPD_40_btn_c +12 R_LMPD_40_btn_c -11 F_LMPD_41_btn_c +11 R_LMPD_41_btn_c -12 F_LMPD_42_btn_c +12 R_LMPD_42_btn_c -16 F_LMPD_43_btn_c +16 R_LMPD_43_btn_c -17 F_LMPD_44_btn_c +17 R_LMPD_44_btn_c -16 F_LMPD_45_btn_c +16 R_LMPD_45_btn_c -17 F_LMPD_47_btn_c +17 R_LMPD_47_btn_c -3 F_LMPD_48_ca2_c +3 R_LMPD_48_ca2_c -4 F_LMPD_49_ca2_c +4 R_LMPD_49_ca2_c -2 F_LMPD_50_ca2_c +2 R_LMPD_50_ca2_c -2 F_LMPD_52_clpn160_p +2 R_LMPD_52_clpn160_p -2 F_LMPD_53_clpn160_p +2 R_LMPD_53_clpn160_p -2 F_LMPD_54_clpn160_p +2 R_LMPD_54_clpn160_p -2 F_LMPD_55_clpn160_p +2 R_LMPD_55_clpn160_p -2 F_LMPD_56_clpn160_p +2 R_LMPD_56_clpn160_p -2 F_LMPD_57_clpn160_p +2 R_LMPD_57_clpn160_p -2 F_LMPD_58_clpn160_p +2 R_LMPD_58_clpn160_p -2 F_LMPD_59_clpn160_p +2 R_LMPD_59_clpn160_p +22 F_LMPD_60_clpn160_p -22 R_LMPD_60_clpn160_p +22 F_LMPD_61_clpn160_p -22 R_LMPD_61_clpn160_p -2 F_LMPD_62_clpn161_p +2 R_LMPD_62_clpn161_p -2 F_LMPD_63_clpn161_p +2 R_LMPD_63_clpn161_p -2 F_LMPD_64_clpn161_p +2 R_LMPD_64_clpn161_p -2 F_LMPD_65_clpn161_p +2 R_LMPD_65_clpn161_p -2 F_LMPD_66_clpn161_p +2 R_LMPD_66_clpn161_p -2 F_LMPD_67_clpn161_p +2 R_LMPD_67_clpn161_p -2 F_LMPD_68_clpn161_p +2 R_LMPD_68_clpn161_p -2 F_LMPD_69_clpn181_p +2 R_LMPD_69_clpn181_p -2 F_LMPD_70_clpn181_p +2 R_LMPD_70_clpn181_p -2 F_LMPD_71_clpn181_p +2 R_LMPD_71_clpn181_p -2 F_LMPD_72_clpn181_p +2 R_LMPD_72_clpn181_p -2 F_LMPD_73_clpn181_p +2 R_LMPD_73_clpn181_p -2 F_LMPD_74_clpn181_p +2 R_LMPD_74_clpn181_p -2 F_LMPD_75_clpn181_p +2 R_LMPD_75_clpn181_p -2 F_LMPD_76_clpn181_p +2 R_LMPD_76_clpn181_p +26 F_LMPD_77_clpn181_p -26 R_LMPD_77_clpn181_p +26 F_LMPD_78_clpn181_p -26 R_LMPD_78_clpn181_p -16 F_LMPD_79_colipa_e +16 R_LMPD_79_colipa_e -16 F_LMPD_80_colipa_e +16 R_LMPD_80_colipa_e -16 F_LMPD_81_colipa_e +16 R_LMPD_81_colipa_e -16 F_LMPD_82_colipa_e +16 R_LMPD_82_colipa_e -16 F_LMPD_83_colipa_e +16 R_LMPD_83_colipa_e -16 F_LMPD_84_colipa_e +16 R_LMPD_84_colipa_e -16 F_LMPD_85_colipa_e +16 R_LMPD_85_colipa_e -16 F_LMPD_86_colipa_e +16 R_LMPD_86_colipa_e -16 F_LMPD_87_colipa_e +16 R_LMPD_87_colipa_e -16 F_LMPD_88_colipa_e +16 R_LMPD_88_colipa_e -16 F_LMPD_89_colipa_e +16 R_LMPD_89_colipa_e -16 F_LMPD_90_colipa_e +16 R_LMPD_90_colipa_e -3 F_LMPD_94_datp_c +3 R_LMPD_94_datp_c -3 F_LMPD_95_datp_c +3 R_LMPD_95_datp_c -2 F_LMPD_96_datp_c +2 R_LMPD_96_datp_c -2 F_LMPD_97_datp_c +2 R_LMPD_97_datp_c + F_LMPD_98_dctp_c - R_LMPD_98_dctp_c -5 F_LMPD_99_dgtp_c +5 R_LMPD_99_dgtp_c -5 F_LMPD_100_dgtp_c +5 R_LMPD_100_dgtp_c -5 F_LMPD_101_dgtp_c +5 R_LMPD_101_dgtp_c -4 F_LMPD_102_dgtp_c +4 R_LMPD_102_dgtp_c -4 F_LMPD_103_dgtp_c +4 R_LMPD_103_dgtp_c -4 F_LMPD_104_dgtp_c +4 R_LMPD_104_dgtp_c +2 F_LMPD_105_dttp_c -2 R_LMPD_105_dttp_c +2 F_LMPD_106_dttp_c -2 R_LMPD_106_dttp_c -6 F_LMPD_107_enter_c +6 R_LMPD_107_enter_c -9 F_LMPD_108_fad_c +9 R_LMPD_108_fad_c -9 F_LMPD_109_fad_c +9 R_LMPD_109_fad_c -11 F_LMPD_110_fad_c +11 R_LMPD_110_fad_c -11 F_LMPD_111_fad_c +11 R_LMPD_111_fad_c -4 F_LMPD_112_GAM_c +4 R_LMPD_112_GAM_c -4 F_LMPD_113_GAM_c +4 R_LMPD_113_GAM_c -3 F_LMPD_114_GAM_c +3 R_LMPD_114_GAM_c -3 F_LMPD_115_GAM_c +3 R_LMPD_115_GAM_c - F_LMPD_117_gly_c + R_LMPD_117_gly_c -2 F_LMPD_118_gthrd_c +2 R_LMPD_118_gthrd_c -2 F_LMPD_119_gthrd_c +2 R_LMPD_119_gthrd_c -6 F_LMPD_120_gtp_c +6 R_LMPD_120_gtp_c -6 F_LMPD_121_gtp_c +6 R_LMPD_121_gtp_c -5 F_LMPD_122_gtp_c +5 R_LMPD_122_gtp_c -5 F_LMPD_123_gtp_c +5 R_LMPD_123_gtp_c +15 F_LMPD_124_hemeO_c -15 R_LMPD_124_hemeO_c +9 F_LMPD_125_hemeO_c -9 R_LMPD_125_hemeO_c -2 F_LMPD_126_his-L_c +2 R_LMPD_126_his-L_c -2 F_LMPD_127_his-L_c +2 R_LMPD_127_his-L_c -2 F_LMPD_130_lipopb_c +2 R_LMPD_130_lipopb_c -2 F_LMPD_131_lipopb_c +2 R_LMPD_131_lipopb_c -4 F_LMPD_132_lipopb_c +4 R_LMPD_132_lipopb_c -2 F_LMPD_133_lipopb_c +2 R_LMPD_133_lipopb_c -3 F_LMPD_135_malcoa_c +3 R_LMPD_135_malcoa_c -3 F_LMPD_136_malcoa_c +3 R_LMPD_136_malcoa_c - F_LMPD_137_met-L_c + R_LMPD_137_met-L_c - F_LMPD_138_met-L_c + R_LMPD_138_met-L_c -8 F_LMPD_139_mlthf_c +8 R_LMPD_139_mlthf_c -8 F_LMPD_140_mlthf_c +8 R_LMPD_140_mlthf_c -5 F_LMPD_141_mococdp_c +5 R_LMPD_141_mococdp_c -7 F_LMPD_142_mococdp_c +7 R_LMPD_142_mococdp_c -5 F_LMPD_143_mococdp_c +5 R_LMPD_143_mococdp_c -7 F_LMPD_144_mococdp_c +7 R_LMPD_144_mococdp_c -6 F_LMPD_145_mococdp_c +6 R_LMPD_145_mococdp_c -8 F_LMPD_146_mococdp_c +8 R_LMPD_146_mococdp_c -6 F_LMPD_147_mococdp_c +6 R_LMPD_147_mococdp_c -8 F_LMPD_148_mococdp_c +8 R_LMPD_148_mococdp_c -10 F_LMPD_149_mocogdp_c +10 R_LMPD_149_mocogdp_c -10 F_LMPD_150_mocogdp_c +10 R_LMPD_150_mocogdp_c -12 F_LMPD_151_mocogdp_c +12 R_LMPD_151_mocogdp_c -12 F_LMPD_152_mocogdp_c +12 R_LMPD_152_mocogdp_c +15 F_LMPD_153_mql8_c -15 R_LMPD_153_mql8_c -4 F_LMPD_154_murein3p3p_p +4 R_LMPD_154_murein3p3p_p -2 F_LMPD_155_murein3px4p_p +2 R_LMPD_155_murein3px4p_p -2 F_LMPD_156_murein4p4p_p +2 R_LMPD_156_murein4p4p_p -2 F_LMPD_157_murein4px4p_p +2 R_LMPD_157_murein4px4p_p -3 F_LMPD_158_murein4px4px4p_p +3 R_LMPD_158_murein4px4px4p_p -2 F_LMPD_159_nad_c +2 R_LMPD_159_nad_c -2 F_LMPD_160_nad_c +2 R_LMPD_160_nad_c - F_LMPD_161_nad_c + R_LMPD_161_nad_c - F_LMPD_162_nad_c + R_LMPD_162_nad_c - F_LMPD_163_pe160_c + R_LMPD_163_pe160_c - F_LMPD_164_pe160_c + R_LMPD_164_pe160_c - F_LMPD_165_pe160_c + R_LMPD_165_pe160_c - F_LMPD_166_pe160_c + R_LMPD_166_pe160_c - F_LMPD_167_pe160_c + R_LMPD_167_pe160_c - F_LMPD_168_pe160_c + R_LMPD_168_pe160_c - F_LMPD_169_pe160_c + R_LMPD_169_pe160_c - F_LMPD_170_pe160_c + R_LMPD_170_pe160_c +11 F_LMPD_171_pe160_c -11 R_LMPD_171_pe160_c +11 F_LMPD_172_pe160_c -11 R_LMPD_172_pe160_c - F_LMPD_173_pe161_c + R_LMPD_173_pe161_c - F_LMPD_174_pe161_c + R_LMPD_174_pe161_c - F_LMPD_175_pe161_c + R_LMPD_175_pe161_c - F_LMPD_176_pe161_c + R_LMPD_176_pe161_c - F_LMPD_177_pe161_c + R_LMPD_177_pe161_c - F_LMPD_178_pe161_c + R_LMPD_178_pe161_c - F_LMPD_179_pe161_c + R_LMPD_179_pe161_c - F_LMPD_180_pe181_c + R_LMPD_180_pe181_c - F_LMPD_181_pe181_c + R_LMPD_181_pe181_c - F_LMPD_182_pe181_c + R_LMPD_182_pe181_c - F_LMPD_183_pe181_c + R_LMPD_183_pe181_c - F_LMPD_184_pe181_c + R_LMPD_184_pe181_c - F_LMPD_185_pe181_c + R_LMPD_185_pe181_c - F_LMPD_186_pe181_c + R_LMPD_186_pe181_c - F_LMPD_187_pe181_c + R_LMPD_187_pe181_c +13 F_LMPD_188_pe181_c -13 R_LMPD_188_pe181_c +13 F_LMPD_189_pe181_c -13 R_LMPD_189_pe181_c - F_LMPD_190_pg160_c + R_LMPD_190_pg160_c - F_LMPD_191_pg160_c + R_LMPD_191_pg160_c - F_LMPD_192_pg160_c + R_LMPD_192_pg160_c - F_LMPD_193_pg160_c + R_LMPD_193_pg160_c - F_LMPD_194_pg160_c + R_LMPD_194_pg160_c - F_LMPD_195_pg160_c + R_LMPD_195_pg160_c - F_LMPD_196_pg160_c + R_LMPD_196_pg160_c - F_LMPD_197_pg160_c + R_LMPD_197_pg160_c +11 F_LMPD_198_pg160_c -11 R_LMPD_198_pg160_c +11 F_LMPD_199_pg160_c -11 R_LMPD_199_pg160_c - F_LMPD_200_pg161_c + R_LMPD_200_pg161_c - F_LMPD_201_pg161_c + R_LMPD_201_pg161_c - F_LMPD_202_pg161_c + R_LMPD_202_pg161_c - F_LMPD_203_pg161_c + R_LMPD_203_pg161_c - F_LMPD_204_pg161_c + R_LMPD_204_pg161_c - F_LMPD_205_pg161_c + R_LMPD_205_pg161_c - F_LMPD_206_pg161_c + R_LMPD_206_pg161_c - F_LMPD_207_pg181_c + R_LMPD_207_pg181_c - F_LMPD_208_pg181_c + R_LMPD_208_pg181_c - F_LMPD_209_pg181_c + R_LMPD_209_pg181_c - F_LMPD_210_pg181_c + R_LMPD_210_pg181_c - F_LMPD_211_pg181_c + R_LMPD_211_pg181_c - F_LMPD_212_pg181_c + R_LMPD_212_pg181_c - F_LMPD_213_pg181_c + R_LMPD_213_pg181_c - F_LMPD_214_pg181_c + R_LMPD_214_pg181_c +13 F_LMPD_215_pg181_c -13 R_LMPD_215_pg181_c +13 F_LMPD_216_pg181_c -13 R_LMPD_216_pg181_c + F_LMPD_217_phe-L_c - R_LMPD_217_phe-L_c +10 F_LMPD_218_pheme_c -10 R_LMPD_218_pheme_c +4 F_LMPD_219_pheme_c -4 R_LMPD_219_pheme_c + F_LMPD_220_pro-L_c - R_LMPD_220_pro-L_c - F_LMPD_221_ptrc_c + R_LMPD_221_ptrc_c + F_LMPD_222_pydx5p_c - R_LMPD_222_pydx5p_c + F_LMPD_223_pydx5p_c - R_LMPD_223_pydx5p_c +7 F_LMPD_224_q8h2_c -7 R_LMPD_224_q8h2_c +4 F_LMPD_225_q8h2_c -4 R_LMPD_225_q8h2_c -2 F_LMPD_226_q8h2_c +2 R_LMPD_226_q8h2_c +10 F_LMPD_227_q8h2_c -10 R_LMPD_227_q8h2_c +7 F_LMPD_228_q8h2_c -7 R_LMPD_228_q8h2_c +7 F_LMPD_229_q8h2_c -7 R_LMPD_229_q8h2_c + F_LMPD_230_q8h2_c - R_LMPD_230_q8h2_c + F_LMPD_231_q8h2_c - R_LMPD_231_q8h2_c +2 F_LMPD_232_q8h2_c -2 R_LMPD_232_q8h2_c +8 F_LMPD_233_q8h2_c -8 R_LMPD_233_q8h2_c +8 F_LMPD_234_q8h2_c -8 R_LMPD_234_q8h2_c -7 F_LMPD_235_ribflv_c +7 R_LMPD_235_ribflv_c -7 F_LMPD_236_ribflv_c +7 R_LMPD_236_ribflv_c -6 F_LMPD_237_ribflv_c +6 R_LMPD_237_ribflv_c -6 F_LMPD_238_ribflv_c +6 R_LMPD_238_ribflv_c - F_LMPD_239_ser-L_c + R_LMPD_239_ser-L_c +4 F_LMPD_240_sheme_c -4 R_LMPD_240_sheme_c -4 F_LMPD_241_spmd_c +4 R_LMPD_241_spmd_c -4 F_LMPD_242_spmd_c +4 R_LMPD_242_spmd_c -10 F_LMPD_243_thf_c +10 R_LMPD_243_thf_c -10 F_LMPD_244_thf_c +10 R_LMPD_244_thf_c -9 F_LMPD_245_thf_c +9 R_LMPD_245_thf_c -9 F_LMPD_246_thf_c +9 R_LMPD_246_thf_c -4 F_LMPD_247_thmpp_c +4 R_LMPD_247_thmpp_c -4 F_LMPD_248_thmpp_c +4 R_LMPD_248_thmpp_c - F_LMPD_249_thr-L_c + R_LMPD_249_thr-L_c +2 F_LMPD_250_trp-L_c -2 R_LMPD_250_trp-L_c +22 F_LMPD_252_udcpdp_c -22 R_LMPD_252_udcpdp_c + F_LMPD_253_utp_c - R_LMPD_253_utp_c + F_LMPD_254_val-L_c - R_LMPD_254_val-L_c = 0
M_h2o2_c: F_ASPO6 - R_ASPO6 + F_MOX - R_MOX - F_THIORDXi + R_THIORDXi + F_LMPD_222_pydx5p_c - R_LMPD_222_pydx5p_c = 0
M_hemeO_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_124_hemeO_c - R_LMPD_124_hemeO_c + F_LMPD_125_hemeO_c - R_LMPD_125_hemeO_c = 0
M_his-L_c: -0.092056 F_Ec_biomass_iJO1366_WT_53p95M +0.092056 R_Ec_biomass_iJO1366_WT_53p95M + F_HISt2rpp - R_HISt2rpp + F_LMPD_126_his-L_c - R_LMPD_126_his-L_c + F_LMPD_127_his-L_c - R_LMPD_127_his-L_c = 0
M_iasp_c: F_ASPO3 - R_ASPO3 + F_ASPO4 - R_ASPO4 + F_ASPO5 - R_ASPO5 + F_ASPO6 - R_ASPO6 - F_LMPD_159_nad_c + R_LMPD_159_nad_c - F_LMPD_160_nad_c + R_LMPD_160_nad_c - F_LMPD_161_nad_c + R_LMPD_161_nad_c - F_LMPD_162_nad_c + R_LMPD_162_nad_c = 0
M_icit_c: F_ACONTb - R_ACONTb - F_ICDHyr + R_ICDHyr - F_ICL + R_ICL = 0
M_ile-L_c: -0.28231 F_Ec_biomass_iJO1366_WT_53p95M +0.28231 R_Ec_biomass_iJO1366_WT_53p95M + F_ILEt2rpp - R_ILEt2rpp + F_LMPD_128_ile-L_c - R_LMPD_128_ile-L_c = 0
M_k_c: -0.18569 F_Ec_biomass_iJO1366_WT_53p95M +0.18569 R_Ec_biomass_iJO1366_WT_53p95M + F_Kabcpp - R_Kabcpp + F_Kt2pp - R_Kt2pp = 0
M_lac-D_c: F_D-LACt2pp - R_D-LACt2pp - F_LDH_D + R_LDH_D - F_LDH_D2 + R_LDH_D2 = 0
M_leu-L_c: -0.43778 F_Ec_biomass_iJO1366_WT_53p95M +0.43778 R_Ec_biomass_iJO1366_WT_53p95M + F_LEUt2rpp - R_LEUt2rpp + F_LMPD_129_leu-L_c - R_LMPD_129_leu-L_c = 0
M_lipopb_c: -3e-06 F_Ec_biomass_iJO1366_WT_53p95M +3e-06 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_130_lipopb_c - R_LMPD_130_lipopb_c + F_LMPD_131_lipopb_c - R_LMPD_131_lipopb_c + F_LMPD_132_lipopb_c - R_LMPD_132_lipopb_c + F_LMPD_133_lipopb_c - R_LMPD_133_lipopb_c = 0
M_lys-L_c: -0.33345 F_Ec_biomass_iJO1366_WT_53p95M +0.33345 R_Ec_biomass_iJO1366_WT_53p95M - F_LYSt3pp + R_LYSt3pp + F_LMPD_134_lys-L_c - R_LMPD_134_lys-L_c = 0
M_mal-L_c: F_FUM - R_FUM + F_MALS - R_MALS - F_MALt3pp + R_MALt3pp - F_MDH + R_MDH - F_MDH2 + R_MDH2 - F_MDH3 + R_MDH3 - F_ME1 + R_ME1 - F_ME2 + R_ME2 - F_MOX + R_MOX = 0
M_malcoa_c: -3.1e-05 F_Ec_biomass_iJO1366_WT_53p95M +3.1e-05 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_135_malcoa_c - R_LMPD_135_malcoa_c + F_LMPD_136_malcoa_c - R_LMPD_136_malcoa_c = 0
M_met-L_c: -0.14934 F_Ec_biomass_iJO1366_WT_53p95M +0.14934 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_43_btn_c - R_LMPD_43_btn_c + F_LMPD_44_btn_c - R_LMPD_44_btn_c + F_LMPD_45_btn_c - R_LMPD_45_btn_c + F_LMPD_47_btn_c - R_LMPD_47_btn_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c = 0
M_mg2_c: -0.008253 F_Ec_biomass_iJO1366_WT_53p95M +0.008253 R_Ec_biomass_iJO1366_WT_53p95M + F_MG2t3_2pp - R_MG2t3_2pp + F_MG2tpp - R_MG2tpp + F_MG2uabcpp - R_MG2uabcpp = 0
M_mlthf_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_139_mlthf_c - R_LMPD_139_mlthf_c + F_LMPD_140_mlthf_c - R_LMPD_140_mlthf_c = 0
M_mn2_c: -0.000658 F_Ec_biomass_iJO1366_WT_53p95M +0.000658 R_Ec_biomass_iJO1366_WT_53p95M + F_MN2tpp - R_MN2tpp + F_MNt2pp - R_MNt2pp = 0
M_mobd_c: -7e-06 F_Ec_biomass_iJO1366_WT_53p95M +7e-06 R_Ec_biomass_iJO1366_WT_53p95M + F_MOBDabcpp - R_MOBDabcpp - F_LMPD_35_bmocogdp_c + R_LMPD_35_bmocogdp_c - F_LMPD_36_bmocogdp_c + R_LMPD_36_bmocogdp_c - F_LMPD_37_bmocogdp_c + R_LMPD_37_bmocogdp_c - F_LMPD_38_bmocogdp_c + R_LMPD_38_bmocogdp_c - F_LMPD_141_mococdp_c + R_LMPD_141_mococdp_c - F_LMPD_142_mococdp_c + R_LMPD_142_mococdp_c - F_LMPD_143_mococdp_c + R_LMPD_143_mococdp_c - F_LMPD_144_mococdp_c + R_LMPD_144_mococdp_c - F_LMPD_145_mococdp_c + R_LMPD_145_mococdp_c - F_LMPD_146_mococdp_c + R_LMPD_146_mococdp_c - F_LMPD_147_mococdp_c + R_LMPD_147_mococdp_c - F_LMPD_148_mococdp_c + R_LMPD_148_mococdp_c - F_LMPD_149_mocogdp_c + R_LMPD_149_mocogdp_c - F_LMPD_150_mocogdp_c + R_LMPD_150_mocogdp_c - F_LMPD_151_mocogdp_c + R_LMPD_151_mocogdp_c - F_LMPD_152_mocogdp_c + R_LMPD_152_mocogdp_c = 0
M_mococdp_c: -7e-06 F_Ec_biomass_iJO1366_WT_53p95M +7e-06 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_141_mococdp_c - R_LMPD_141_mococdp_c + F_LMPD_142_mococdp_c - R_LMPD_142_mococdp_c + F_LMPD_143_mococdp_c - R_LMPD_143_mococdp_c + F_LMPD_144_mococdp_c - R_LMPD_144_mococdp_c + F_LMPD_145_mococdp_c - R_LMPD_145_mococdp_c + F_LMPD_146_mococdp_c - R_LMPD_146_mococdp_c + F_LMPD_147_mococdp_c - R_LMPD_147_mococdp_c + F_LMPD_148_mococdp_c - R_LMPD_148_mococdp_c = 0
M_mocogdp_c: -7e-06 F_Ec_biomass_iJO1366_WT_53p95M +7e-06 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_149_mocogdp_c - R_LMPD_149_mocogdp_c + F_LMPD_150_mocogdp_c - R_LMPD_150_mocogdp_c + F_LMPD_151_mocogdp_c - R_LMPD_151_mocogdp_c + F_LMPD_152_mocogdp_c - R_LMPD_152_mocogdp_c = 0
M_mql8_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_ASPO4 - R_ASPO4 - F_CYTBD2pp + R_CYTBD2pp + F_DHORD5 - R_DHORD5 + F_FDH5pp - R_FDH5pp + F_G3PD6 - R_G3PD6 + F_GLYCTO3 - R_GLYCTO3 + F_MDH3 - R_MDH3 + F_NADH10 - R_NADH10 + F_NADH17pp - R_NADH17pp + F_NADPHQR3 - R_NADPHQR3 + F_LMPD_153_mql8_c - R_LMPD_153_mql8_c = 0
M_mqn8_c: - F_ASPO4 + R_ASPO4 + F_CYTBD2pp - R_CYTBD2pp - F_DHORD5 + R_DHORD5 - F_FDH5pp + R_FDH5pp - F_G3PD6 + R_G3PD6 - F_GLYCTO3 + R_GLYCTO3 - F_MDH3 + R_MDH3 - F_NADH10 + R_NADH10 - F_NADH17pp + R_NADH17pp - F_NADPHQR3 + R_NADPHQR3 = 0
M_mththf_c: - F_DM_MTHTHF + R_DM_MTHTHF + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c + F_LMPD_43_btn_c - R_LMPD_43_btn_c + F_LMPD_44_btn_c - R_LMPD_44_btn_c + F_LMPD_45_btn_c - R_LMPD_45_btn_c + F_LMPD_47_btn_c - R_LMPD_47_btn_c + F_LMPD_153_mql8_c - R_LMPD_153_mql8_c +3 F_LMPD_224_q8h2_c -3 R_LMPD_224_q8h2_c +3 F_LMPD_225_q8h2_c -3 R_LMPD_225_q8h2_c +3 F_LMPD_226_q8h2_c -3 R_LMPD_226_q8h2_c +3 F_LMPD_227_q8h2_c -3 R_LMPD_227_q8h2_c +3 F_LMPD_228_q8h2_c -3 R_LMPD_228_q8h2_c +3 F_LMPD_229_q8h2_c -3 R_LMPD_229_q8h2_c +3 F_LMPD_230_q8h2_c -3 R_LMPD_230_q8h2_c +3 F_LMPD_231_q8h2_c -3 R_LMPD_231_q8h2_c +6 F_LMPD_232_q8h2_c -6 R_LMPD_232_q8h2_c +6 F_LMPD_233_q8h2_c -6 R_LMPD_233_q8h2_c +6 F_LMPD_234_q8h2_c -6 R_LMPD_234_q8h2_c +2 F_LMPD_240_sheme_c -2 R_LMPD_240_sheme_c = 0
M_na1_c: F_CAt6pp - R_CAt6pp + F_GLUt4pp - R_GLUt4pp -2 F_NAt3_1p5pp +2 R_NAt3_1p5pp - F_NAt3_2pp + R_NAt3_2pp - F_NAt3pp + R_NAt3pp = 0
M_nad_c: -0.001787 F_Ec_biomass_iJO1366_WT_53p95M +0.001787 R_Ec_biomass_iJO1366_WT_53p95M - F_ACALD + R_ACALD - F_AKGDH + R_AKGDH - F_ALCD2x + R_ALCD2x - F_ALDD2x + R_ALDD2x - F_GAPD + R_GAPD - F_GLYCDx + R_GLYCDx + F_GLYCLTDx - R_GLYCLTDx - F_LDH_D + R_LDH_D - F_MDH + R_MDH - F_ME1 + R_ME1 + F_NADH10 - R_NADH10 + F_NADH16pp - R_NADH16pp + F_NADH17pp - R_NADH17pp + F_NADH18pp - R_NADH18pp + F_NADH5 - R_NADH5 + F_NADH9 - R_NADH9 - F_NADK + R_NADK + F_NADPPPS - R_NADPPPS - F_NADTRHD + R_NADTRHD - F_PDH + R_PDH + F_THD2pp - R_THD2pp + F_TRSARr - R_TRSARr -4 F_LMPD_1_10fthf_c +4 R_LMPD_1_10fthf_c -4 F_LMPD_2_10fthf_c +4 R_LMPD_2_10fthf_c -5 F_LMPD_3_10fthf_c +5 R_LMPD_3_10fthf_c -5 F_LMPD_4_10fthf_c +5 R_LMPD_4_10fthf_c +8 F_LMPD_5_2dmmql8_c -8 R_LMPD_5_2dmmql8_c -2 F_LMPD_6_2fe2s_c +2 R_LMPD_6_2fe2s_c - F_LMPD_7_2fe2s_c + R_LMPD_7_2fe2s_c -2 F_LMPD_8_2fe2s_c +2 R_LMPD_8_2fe2s_c - F_LMPD_9_2fe2s_c + R_LMPD_9_2fe2s_c -2 F_LMPD_10_2fe2s_c +2 R_LMPD_10_2fe2s_c - F_LMPD_11_2fe2s_c + R_LMPD_11_2fe2s_c -2 F_LMPD_12_2fe2s_c +2 R_LMPD_12_2fe2s_c - F_LMPD_13_2fe2s_c + R_LMPD_13_2fe2s_c -4 F_LMPD_14_4fe4s_c +4 R_LMPD_14_4fe4s_c - F_LMPD_15_4fe4s_c + R_LMPD_15_4fe4s_c -4 F_LMPD_16_4fe4s_c +4 R_LMPD_16_4fe4s_c - F_LMPD_17_4fe4s_c + R_LMPD_17_4fe4s_c -4 F_LMPD_18_4fe4s_c +4 R_LMPD_18_4fe4s_c - F_LMPD_19_4fe4s_c + R_LMPD_19_4fe4s_c -4 F_LMPD_20_4fe4s_c +4 R_LMPD_20_4fe4s_c - F_LMPD_21_4fe4s_c + R_LMPD_21_4fe4s_c -2 F_LMPD_22_5mthf_c +2 R_LMPD_22_5mthf_c -2 F_LMPD_23_5mthf_c +2 R_LMPD_23_5mthf_c -2 F_LMPD_24_accoa_c +2 R_LMPD_24_accoa_c -2 F_LMPD_25_accoa_c +2 R_LMPD_25_accoa_c -2 F_LMPD_26_adocbl_c +2 R_LMPD_26_adocbl_c -2 F_LMPD_27_adocbl_c +2 R_LMPD_27_adocbl_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c -12 F_LMPD_35_bmocogdp_c +12 R_LMPD_35_bmocogdp_c -12 F_LMPD_36_bmocogdp_c +12 R_LMPD_36_bmocogdp_c -8 F_LMPD_37_bmocogdp_c +8 R_LMPD_37_bmocogdp_c -8 F_LMPD_38_bmocogdp_c +8 R_LMPD_38_bmocogdp_c -2 F_LMPD_39_btn_c +2 R_LMPD_39_btn_c -2 F_LMPD_40_btn_c +2 R_LMPD_40_btn_c -2 F_LMPD_41_btn_c +2 R_LMPD_41_btn_c -2 F_LMPD_42_btn_c +2 R_LMPD_42_btn_c -3 F_LMPD_43_btn_c +3 R_LMPD_43_btn_c -3 F_LMPD_44_btn_c +3 R_LMPD_44_btn_c -3 F_LMPD_45_btn_c +3 R_LMPD_45_btn_c -3 F_LMPD_47_btn_c +3 R_LMPD_47_btn_c - F_LMPD_48_ca2_c + R_LMPD_48_ca2_c - F_LMPD_49_ca2_c + R_LMPD_49_ca2_c +12 F_LMPD_52_clpn160_p -12 R_LMPD_52_clpn160_p +16 F_LMPD_53_clpn160_p -16 R_LMPD_53_clpn160_p +20 F_LMPD_54_clpn160_p -20 R_LMPD_54_clpn160_p +24 F_LMPD_55_clpn160_p -24 R_LMPD_55_clpn160_p +28 F_LMPD_56_clpn160_p -28 R_LMPD_56_clpn160_p +4 F_LMPD_57_clpn160_p -4 R_LMPD_57_clpn160_p +8 F_LMPD_59_clpn160_p -8 R_LMPD_59_clpn160_p +28 F_LMPD_60_clpn160_p -28 R_LMPD_60_clpn160_p +56 F_LMPD_61_clpn160_p -56 R_LMPD_61_clpn160_p +12 F_LMPD_62_clpn161_p -12 R_LMPD_62_clpn161_p +16 F_LMPD_63_clpn161_p -16 R_LMPD_63_clpn161_p +20 F_LMPD_64_clpn161_p -20 R_LMPD_64_clpn161_p +24 F_LMPD_65_clpn161_p -24 R_LMPD_65_clpn161_p +4 F_LMPD_66_clpn161_p -4 R_LMPD_66_clpn161_p +8 F_LMPD_68_clpn161_p -8 R_LMPD_68_clpn161_p +12 F_LMPD_69_clpn181_p -12 R_LMPD_69_clpn181_p +16 F_LMPD_70_clpn181_p -16 R_LMPD_70_clpn181_p +20 F_LMPD_71_clpn181_p -20 R_LMPD_71_clpn181_p +24 F_LMPD_72_clpn181_p -24 R_LMPD_72_clpn181_p +28 F_LMPD_73_clpn181_p -28 R_LMPD_73_clpn181_p +4 F_LMPD_74_clpn181_p -4 R_LMPD_74_clpn181_p +8 F_LMPD_76_clpn181_p -8 R_LMPD_76_clpn181_p +32 F_LMPD_77_clpn181_p -32 R_LMPD_77_clpn181_p +60 F_LMPD_78_clpn181_p -60 R_LMPD_78_clpn181_p +12 F_LMPD_79_colipa_e -12 R_LMPD_79_colipa_e +13 F_LMPD_80_colipa_e -13 R_LMPD_80_colipa_e +18 F_LMPD_81_colipa_e -18 R_LMPD_81_colipa_e +19 F_LMPD_82_colipa_e -19 R_LMPD_82_colipa_e +24 F_LMPD_83_colipa_e -24 R_LMPD_83_colipa_e +25 F_LMPD_84_colipa_e -25 R_LMPD_84_colipa_e +30 F_LMPD_85_colipa_e -30 R_LMPD_85_colipa_e +31 F_LMPD_86_colipa_e -31 R_LMPD_86_colipa_e +6 F_LMPD_87_colipa_e -6 R_LMPD_87_colipa_e +7 F_LMPD_89_colipa_e -7 R_LMPD_89_colipa_e + F_LMPD_90_colipa_e - R_LMPD_90_colipa_e - F_LMPD_92_cys-L_c + R_LMPD_92_cys-L_c - F_LMPD_93_cys-L_c + R_LMPD_93_cys-L_c -2 F_LMPD_94_datp_c +2 R_LMPD_94_datp_c -2 F_LMPD_95_datp_c +2 R_LMPD_95_datp_c - F_LMPD_96_datp_c + R_LMPD_96_datp_c - F_LMPD_97_datp_c + R_LMPD_97_datp_c -3 F_LMPD_99_dgtp_c +3 R_LMPD_99_dgtp_c -3 F_LMPD_100_dgtp_c +3 R_LMPD_100_dgtp_c -3 F_LMPD_101_dgtp_c +3 R_LMPD_101_dgtp_c -2 F_LMPD_102_dgtp_c +2 R_LMPD_102_dgtp_c -2 F_LMPD_103_dgtp_c +2 R_LMPD_103_dgtp_c -2 F_LMPD_104_dgtp_c +2 R_LMPD_104_dgtp_c -6 F_LMPD_107_enter_c +6 R_LMPD_107_enter_c -3 F_LMPD_108_fad_c +3 R_LMPD_108_fad_c -3 F_LMPD_109_fad_c +3 R_LMPD_109_fad_c -5 F_LMPD_110_fad_c +5 R_LMPD_110_fad_c -5 F_LMPD_111_fad_c +5 R_LMPD_111_fad_c -2 F_LMPD_112_GAM_c +2 R_LMPD_112_GAM_c -2 F_LMPD_113_GAM_c +2 R_LMPD_113_GAM_c - F_LMPD_114_GAM_c + R_LMPD_114_GAM_c - F_LMPD_115_GAM_c + R_LMPD_115_GAM_c -2 F_LMPD_118_gthrd_c +2 R_LMPD_118_gthrd_c -2 F_LMPD_119_gthrd_c +2 R_LMPD_119_gthrd_c -3 F_LMPD_120_gtp_c +3 R_LMPD_120_gtp_c -3 F_LMPD_121_gtp_c +3 R_LMPD_121_gtp_c -2 F_LMPD_122_gtp_c +2 R_LMPD_122_gtp_c -2 F_LMPD_123_gtp_c +2 R_LMPD_123_gtp_c +3 F_LMPD_124_hemeO_c -3 R_LMPD_124_hemeO_c +3 F_LMPD_125_hemeO_c -3 R_LMPD_125_hemeO_c -2 F_LMPD_126_his-L_c +2 R_LMPD_126_his-L_c -2 F_LMPD_127_his-L_c +2 R_LMPD_127_his-L_c - F_LMPD_129_leu-L_c + R_LMPD_129_leu-L_c +4 F_LMPD_130_lipopb_c -4 R_LMPD_130_lipopb_c +4 F_LMPD_131_lipopb_c -4 R_LMPD_131_lipopb_c +4 F_LMPD_132_lipopb_c -4 R_LMPD_132_lipopb_c -2 F_LMPD_133_lipopb_c +2 R_LMPD_133_lipopb_c -2 F_LMPD_135_malcoa_c +2 R_LMPD_135_malcoa_c -2 F_LMPD_136_malcoa_c +2 R_LMPD_136_malcoa_c - F_LMPD_137_met-L_c + R_LMPD_137_met-L_c - F_LMPD_138_met-L_c + R_LMPD_138_met-L_c -3 F_LMPD_139_mlthf_c +3 R_LMPD_139_mlthf_c -3 F_LMPD_140_mlthf_c +3 R_LMPD_140_mlthf_c -2 F_LMPD_141_mococdp_c +2 R_LMPD_141_mococdp_c -2 F_LMPD_142_mococdp_c +2 R_LMPD_142_mococdp_c -2 F_LMPD_143_mococdp_c +2 R_LMPD_143_mococdp_c -2 F_LMPD_144_mococdp_c +2 R_LMPD_144_mococdp_c -3 F_LMPD_145_mococdp_c +3 R_LMPD_145_mococdp_c -3 F_LMPD_146_mococdp_c +3 R_LMPD_146_mococdp_c -3 F_LMPD_147_mococdp_c +3 R_LMPD_147_mococdp_c -3 F_LMPD_148_mococdp_c +3 R_LMPD_148_mococdp_c -4 F_LMPD_149_mocogdp_c +4 R_LMPD_149_mocogdp_c -4 F_LMPD_150_mocogdp_c +4 R_LMPD_150_mocogdp_c -6 F_LMPD_151_mocogdp_c +6 R_LMPD_151_mocogdp_c -6 F_LMPD_152_mocogdp_c +6 R_LMPD_152_mocogdp_c +9 F_LMPD_153_mql8_c -9 R_LMPD_153_mql8_c - F_LMPD_159_nad_c + R_LMPD_159_nad_c - F_LMPD_160_nad_c + R_LMPD_160_nad_c +11 F_LMPD_163_pe160_c -11 R_LMPD_163_pe160_c +13 F_LMPD_164_pe160_c -13 R_LMPD_164_pe160_c +3 F_LMPD_165_pe160_c -3 R_LMPD_165_pe160_c +5 F_LMPD_166_pe160_c -5 R_LMPD_166_pe160_c +7 F_LMPD_167_pe160_c -7 R_LMPD_167_pe160_c +9 F_LMPD_168_pe160_c -9 R_LMPD_168_pe160_c - F_LMPD_169_pe160_c + R_LMPD_169_pe160_c + F_LMPD_170_pe160_c - R_LMPD_170_pe160_c +13 F_LMPD_171_pe160_c -13 R_LMPD_171_pe160_c +27 F_LMPD_172_pe160_c -27 R_LMPD_172_pe160_c +11 F_LMPD_173_pe161_c -11 R_LMPD_173_pe161_c +3 F_LMPD_174_pe161_c -3 R_LMPD_174_pe161_c +5 F_LMPD_175_pe161_c -5 R_LMPD_175_pe161_c +7 F_LMPD_176_pe161_c -7 R_LMPD_176_pe161_c +9 F_LMPD_177_pe161_c -9 R_LMPD_177_pe161_c - F_LMPD_178_pe161_c + R_LMPD_178_pe161_c + F_LMPD_179_pe161_c - R_LMPD_179_pe161_c +11 F_LMPD_180_pe181_c -11 R_LMPD_180_pe181_c +13 F_LMPD_181_pe181_c -13 R_LMPD_181_pe181_c +3 F_LMPD_182_pe181_c -3 R_LMPD_182_pe181_c +5 F_LMPD_183_pe181_c -5 R_LMPD_183_pe181_c +7 F_LMPD_184_pe181_c -7 R_LMPD_184_pe181_c +9 F_LMPD_185_pe181_c -9 R_LMPD_185_pe181_c - F_LMPD_186_pe181_c + R_LMPD_186_pe181_c + F_LMPD_187_pe181_c - R_LMPD_187_pe181_c +15 F_LMPD_188_pe181_c -15 R_LMPD_188_pe181_c +29 F_LMPD_189_pe181_c -29 R_LMPD_189_pe181_c +10 F_LMPD_190_pg160_c -10 R_LMPD_190_pg160_c +12 F_LMPD_191_pg160_c -12 R_LMPD_191_pg160_c +14 F_LMPD_192_pg160_c -14 R_LMPD_192_pg160_c +2 F_LMPD_193_pg160_c -2 R_LMPD_193_pg160_c +4 F_LMPD_195_pg160_c -4 R_LMPD_195_pg160_c +6 F_LMPD_196_pg160_c -6 R_LMPD_196_pg160_c +8 F_LMPD_197_pg160_c -8 R_LMPD_197_pg160_c +14 F_LMPD_198_pg160_c -14 R_LMPD_198_pg160_c +28 F_LMPD_199_pg160_c -28 R_LMPD_199_pg160_c +10 F_LMPD_200_pg161_c -10 R_LMPD_200_pg161_c +12 F_LMPD_201_pg161_c -12 R_LMPD_201_pg161_c +2 F_LMPD_202_pg161_c -2 R_LMPD_202_pg161_c +4 F_LMPD_204_pg161_c -4 R_LMPD_204_pg161_c +6 F_LMPD_205_pg161_c -6 R_LMPD_205_pg161_c +8 F_LMPD_206_pg161_c -8 R_LMPD_206_pg161_c +10 F_LMPD_207_pg181_c -10 R_LMPD_207_pg181_c +12 F_LMPD_208_pg181_c -12 R_LMPD_208_pg181_c +14 F_LMPD_209_pg181_c -14 R_LMPD_209_pg181_c +2 F_LMPD_210_pg181_c -2 R_LMPD_210_pg181_c +4 F_LMPD_212_pg181_c -4 R_LMPD_212_pg181_c +6 F_LMPD_213_pg181_c -6 R_LMPD_213_pg181_c +8 F_LMPD_214_pg181_c -8 R_LMPD_214_pg181_c +16 F_LMPD_215_pg181_c -16 R_LMPD_215_pg181_c +30 F_LMPD_216_pg181_c -30 R_LMPD_216_pg181_c -3 F_LMPD_222_pydx5p_c +3 R_LMPD_222_pydx5p_c -4 F_LMPD_223_pydx5p_c +4 R_LMPD_223_pydx5p_c +8 F_LMPD_224_q8h2_c -8 R_LMPD_224_q8h2_c +7 F_LMPD_225_q8h2_c -7 R_LMPD_225_q8h2_c +5 F_LMPD_226_q8h2_c -5 R_LMPD_226_q8h2_c +11 F_LMPD_227_q8h2_c -11 R_LMPD_227_q8h2_c +10 F_LMPD_228_q8h2_c -10 R_LMPD_228_q8h2_c +10 F_LMPD_229_q8h2_c -10 R_LMPD_229_q8h2_c +8 F_LMPD_230_q8h2_c -8 R_LMPD_230_q8h2_c +8 F_LMPD_231_q8h2_c -8 R_LMPD_231_q8h2_c +12 F_LMPD_232_q8h2_c -12 R_LMPD_232_q8h2_c +18 F_LMPD_233_q8h2_c -18 R_LMPD_233_q8h2_c +18 F_LMPD_234_q8h2_c -18 R_LMPD_234_q8h2_c -3 F_LMPD_235_ribflv_c +3 R_LMPD_235_ribflv_c -3 F_LMPD_236_ribflv_c +3 R_LMPD_236_ribflv_c -2 F_LMPD_237_ribflv_c +2 R_LMPD_237_ribflv_c -2 F_LMPD_238_ribflv_c +2 R_LMPD_238_ribflv_c - F_LMPD_239_ser-L_c + R_LMPD_239_ser-L_c + F_LMPD_240_sheme_c - R_LMPD_240_sheme_c - F_LMPD_241_spmd_c + R_LMPD_241_spmd_c - F_LMPD_242_spmd_c + R_LMPD_242_spmd_c -4 F_LMPD_243_thf_c +4 R_LMPD_243_thf_c -4 F_LMPD_244_thf_c +4 R_LMPD_244_thf_c -3 F_LMPD_245_thf_c +3 R_LMPD_245_thf_c -3 F_LMPD_246_thf_c +3 R_LMPD_246_thf_c -4 F_LMPD_247_thmpp_c +4 R_LMPD_247_thmpp_c -4 F_LMPD_248_thmpp_c +4 R_LMPD_248_thmpp_c - F_LMPD_251_tyr-L_c + R_LMPD_251_tyr-L_c +11 F_LMPD_252_udcpdp_c -11 R_LMPD_252_udcpdp_c = 0
M_nadh_c: -4.5e-05 F_Ec_biomass_iJO1366_WT_53p95M +4.5e-05 R_Ec_biomass_iJO1366_WT_53p95M + F_ACALD - R_ACALD + F_AKGDH - R_AKGDH + F_ALCD2x - R_ALCD2x + F_ALDD2x - R_ALDD2x + F_GAPD - R_GAPD + F_GLYCDx - R_GLYCDx - F_GLYCLTDx + R_GLYCLTDx + F_LDH_D - R_LDH_D + F_MDH - R_MDH + F_ME1 - R_ME1 - F_NADH10 + R_NADH10 - F_NADH16pp + R_NADH16pp - F_NADH17pp + R_NADH17pp - F_NADH18pp + R_NADH18pp - F_NADH5 + R_NADH5 - F_NADH9 + R_NADH9 + F_NADTRHD - R_NADTRHD + F_PDH - R_PDH - F_THD2pp + R_THD2pp - F_TRSARr + R_TRSARr +4 F_LMPD_1_10fthf_c -4 R_LMPD_1_10fthf_c +4 F_LMPD_2_10fthf_c -4 R_LMPD_2_10fthf_c +5 F_LMPD_3_10fthf_c -5 R_LMPD_3_10fthf_c +5 F_LMPD_4_10fthf_c -5 R_LMPD_4_10fthf_c -8 F_LMPD_5_2dmmql8_c +8 R_LMPD_5_2dmmql8_c +2 F_LMPD_6_2fe2s_c -2 R_LMPD_6_2fe2s_c + F_LMPD_7_2fe2s_c - R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c + F_LMPD_9_2fe2s_c - R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c + F_LMPD_11_2fe2s_c - R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c + F_LMPD_13_2fe2s_c - R_LMPD_13_2fe2s_c +4 F_LMPD_14_4fe4s_c -4 R_LMPD_14_4fe4s_c + F_LMPD_15_4fe4s_c - R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c + F_LMPD_17_4fe4s_c - R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c + F_LMPD_19_4fe4s_c - R_LMPD_19_4fe4s_c +4 F_LMPD_20_4fe4s_c -4 R_LMPD_20_4fe4s_c + F_LMPD_21_4fe4s_c - R_LMPD_21_4fe4s_c +2 F_LMPD_22_5mthf_c -2 R_LMPD_22_5mthf_c +2 F_LMPD_23_5mthf_c -2 R_LMPD_23_5mthf_c +2 F_LMPD_24_accoa_c -2 R_LMPD_24_accoa_c +2 F_LMPD_25_accoa_c -2 R_LMPD_25_accoa_c +2 F_LMPD_26_adocbl_c -2 R_LMPD_26_adocbl_c +2 F_LMPD_27_adocbl_c -2 R_LMPD_27_adocbl_c + F_LMPD_28_adocbl_c - R_LMPD_28_adocbl_c + F_LMPD_29_adocbl_c - R_LMPD_29_adocbl_c + F_LMPD_31_amet_c - R_LMPD_31_amet_c + F_LMPD_32_amet_c - R_LMPD_32_amet_c +12 F_LMPD_35_bmocogdp_c -12 R_LMPD_35_bmocogdp_c +12 F_LMPD_36_bmocogdp_c -12 R_LMPD_36_bmocogdp_c +8 F_LMPD_37_bmocogdp_c -8 R_LMPD_37_bmocogdp_c +8 F_LMPD_38_bmocogdp_c -8 R_LMPD_38_bmocogdp_c +2 F_LMPD_39_btn_c -2 R_LMPD_39_btn_c +2 F_LMPD_40_btn_c -2 R_LMPD_40_btn_c +2 F_LMPD_41_btn_c -2 R_LMPD_41_btn_c +2 F_LMPD_42_btn_c -2 R_LMPD_42_btn_c +3 F_LMPD_43_btn_c -3 R_LMPD_43_btn_c +3 F_LMPD_44_btn_c -3 R_LMPD_44_btn_c +3 F_LMPD_45_btn_c -3 R_LMPD_45_btn_c +3 F_LMPD_47_btn_c -3 R_LMPD_47_btn_c + F_LMPD_48_ca2_c - R_LMPD_48_ca2_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c -12 F_LMPD_52_clpn160_p +12 R_LMPD_52_clpn160_p -16 F_LMPD_53_clpn160_p +16 R_LMPD_53_clpn160_p -20 F_LMPD_54_clpn160_p +20 R_LMPD_54_clpn160_p -24 F_LMPD_55_clpn160_p +24 R_LMPD_55_clpn160_p -28 F_LMPD_56_clpn160_p +28 R_LMPD_56_clpn160_p -4 F_LMPD_57_clpn160_p +4 R_LMPD_57_clpn160_p -8 F_LMPD_59_clpn160_p +8 R_LMPD_59_clpn160_p -28 F_LMPD_60_clpn160_p +28 R_LMPD_60_clpn160_p -56 F_LMPD_61_clpn160_p +56 R_LMPD_61_clpn160_p -12 F_LMPD_62_clpn161_p +12 R_LMPD_62_clpn161_p -16 F_LMPD_63_clpn161_p +16 R_LMPD_63_clpn161_p -20 F_LMPD_64_clpn161_p +20 R_LMPD_64_clpn161_p -24 F_LMPD_65_clpn161_p +24 R_LMPD_65_clpn161_p -4 F_LMPD_66_clpn161_p +4 R_LMPD_66_clpn161_p -8 F_LMPD_68_clpn161_p +8 R_LMPD_68_clpn161_p -12 F_LMPD_69_clpn181_p +12 R_LMPD_69_clpn181_p -16 F_LMPD_70_clpn181_p +16 R_LMPD_70_clpn181_p -20 F_LMPD_71_clpn181_p +20 R_LMPD_71_clpn181_p -24 F_LMPD_72_clpn181_p +24 R_LMPD_72_clpn181_p -28 F_LMPD_73_clpn181_p +28 R_LMPD_73_clpn181_p -4 F_LMPD_74_clpn181_p +4 R_LMPD_74_clpn181_p -8 F_LMPD_76_clpn181_p +8 R_LMPD_76_clpn181_p -32 F_LMPD_77_clpn181_p +32 R_LMPD_77_clpn181_p -60 F_LMPD_78_clpn181_p +60 R_LMPD_78_clpn181_p -12 F_LMPD_79_colipa_e +12 R_LMPD_79_colipa_e -13 F_LMPD_80_colipa_e +13 R_LMPD_80_colipa_e -18 F_LMPD_81_colipa_e +18 R_LMPD_81_colipa_e -19 F_LMPD_82_colipa_e +19 R_LMPD_82_colipa_e -24 F_LMPD_83_colipa_e +24 R_LMPD_83_colipa_e -25 F_LMPD_84_colipa_e +25 R_LMPD_84_colipa_e -30 F_LMPD_85_colipa_e +30 R_LMPD_85_colipa_e -31 F_LMPD_86_colipa_e +31 R_LMPD_86_colipa_e -6 F_LMPD_87_colipa_e +6 R_LMPD_87_colipa_e -7 F_LMPD_89_colipa_e +7 R_LMPD_89_colipa_e - F_LMPD_90_colipa_e + R_LMPD_90_colipa_e + F_LMPD_92_cys-L_c - R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c +2 F_LMPD_94_datp_c -2 R_LMPD_94_datp_c +2 F_LMPD_95_datp_c -2 R_LMPD_95_datp_c + F_LMPD_96_datp_c - R_LMPD_96_datp_c + F_LMPD_97_datp_c - R_LMPD_97_datp_c +3 F_LMPD_99_dgtp_c -3 R_LMPD_99_dgtp_c +3 F_LMPD_100_dgtp_c -3 R_LMPD_100_dgtp_c +3 F_LMPD_101_dgtp_c -3 R_LMPD_101_dgtp_c +2 F_LMPD_102_dgtp_c -2 R_LMPD_102_dgtp_c +2 F_LMPD_103_dgtp_c -2 R_LMPD_103_dgtp_c +2 F_LMPD_104_dgtp_c -2 R_LMPD_104_dgtp_c +6 F_LMPD_107_enter_c -6 R_LMPD_107_enter_c +3 F_LMPD_108_fad_c -3 R_LMPD_108_fad_c +3 F_LMPD_109_fad_c -3 R_LMPD_109_fad_c +5 F_LMPD_110_fad_c -5 R_LMPD_110_fad_c +5 F_LMPD_111_fad_c -5 R_LMPD_111_fad_c +2 F_LMPD_112_GAM_c -2 R_LMPD_112_GAM_c +2 F_LMPD_113_GAM_c -2 R_LMPD_113_GAM_c + F_LMPD_114_GAM_c - R_LMPD_114_GAM_c + F_LMPD_115_GAM_c - R_LMPD_115_GAM_c +2 F_LMPD_118_gthrd_c -2 R_LMPD_118_gthrd_c +2 F_LMPD_119_gthrd_c -2 R_LMPD_119_gthrd_c +3 F_LMPD_120_gtp_c -3 R_LMPD_120_gtp_c +3 F_LMPD_121_gtp_c -3 R_LMPD_121_gtp_c +2 F_LMPD_122_gtp_c -2 R_LMPD_122_gtp_c +2 F_LMPD_123_gtp_c -2 R_LMPD_123_gtp_c -3 F_LMPD_124_hemeO_c +3 R_LMPD_124_hemeO_c -3 F_LMPD_125_hemeO_c +3 R_LMPD_125_hemeO_c +2 F_LMPD_126_his-L_c -2 R_LMPD_126_his-L_c +2 F_LMPD_127_his-L_c -2 R_LMPD_127_his-L_c + F_LMPD_129_leu-L_c - R_LMPD_129_leu-L_c -4 F_LMPD_130_lipopb_c +4 R_LMPD_130_lipopb_c -4 F_LMPD_131_lipopb_c +4 R_LMPD_131_lipopb_c -4 F_LMPD_132_lipopb_c +4 R_LMPD_132_lipopb_c +2 F_LMPD_133_lipopb_c -2 R_LMPD_133_lipopb_c +2 F_LMPD_135_malcoa_c -2 R_LMPD_135_malcoa_c +2 F_LMPD_136_malcoa_c -2 R_LMPD_136_malcoa_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c +3 F_LMPD_139_mlthf_c -3 R_LMPD_139_mlthf_c +3 F_LMPD_140_mlthf_c -3 R_LMPD_140_mlthf_c +2 F_LMPD_141_mococdp_c -2 R_LMPD_141_mococdp_c +2 F_LMPD_142_mococdp_c -2 R_LMPD_142_mococdp_c +2 F_LMPD_143_mococdp_c -2 R_LMPD_143_mococdp_c +2 F_LMPD_144_mococdp_c -2 R_LMPD_144_mococdp_c +3 F_LMPD_145_mococdp_c -3 R_LMPD_145_mococdp_c +3 F_LMPD_146_mococdp_c -3 R_LMPD_146_mococdp_c +3 F_LMPD_147_mococdp_c -3 R_LMPD_147_mococdp_c +3 F_LMPD_148_mococdp_c -3 R_LMPD_148_mococdp_c +4 F_LMPD_149_mocogdp_c -4 R_LMPD_149_mocogdp_c +4 F_LMPD_150_mocogdp_c -4 R_LMPD_150_mocogdp_c +6 F_LMPD_151_mocogdp_c -6 R_LMPD_151_mocogdp_c +6 F_LMPD_152_mocogdp_c -6 R_LMPD_152_mocogdp_c -9 F_LMPD_153_mql8_c +9 R_LMPD_153_mql8_c +2 F_LMPD_159_nad_c -2 R_LMPD_159_nad_c +2 F_LMPD_160_nad_c -2 R_LMPD_160_nad_c + F_LMPD_161_nad_c - R_LMPD_161_nad_c + F_LMPD_162_nad_c - R_LMPD_162_nad_c -11 F_LMPD_163_pe160_c +11 R_LMPD_163_pe160_c -13 F_LMPD_164_pe160_c +13 R_LMPD_164_pe160_c -3 F_LMPD_165_pe160_c +3 R_LMPD_165_pe160_c -5 F_LMPD_166_pe160_c +5 R_LMPD_166_pe160_c -7 F_LMPD_167_pe160_c +7 R_LMPD_167_pe160_c -9 F_LMPD_168_pe160_c +9 R_LMPD_168_pe160_c + F_LMPD_169_pe160_c - R_LMPD_169_pe160_c - F_LMPD_170_pe160_c + R_LMPD_170_pe160_c -13 F_LMPD_171_pe160_c +13 R_LMPD_171_pe160_c -27 F_LMPD_172_pe160_c +27 R_LMPD_172_pe160_c -11 F_LMPD_173_pe161_c +11 R_LMPD_173_pe161_c -3 F_LMPD_174_pe161_c +3 R_LMPD_174_pe161_c -5 F_LMPD_175_pe161_c +5 R_LMPD_175_pe161_c -7 F_LMPD_176_pe161_c +7 R_LMPD_176_pe161_c -9 F_LMPD_177_pe161_c +9 R_LMPD_177_pe161_c + F_LMPD_178_pe161_c - R_LMPD_178_pe161_c - F_LMPD_179_pe161_c + R_LMPD_179_pe161_c -11 F_LMPD_180_pe181_c +11 R_LMPD_180_pe181_c -13 F_LMPD_181_pe181_c +13 R_LMPD_181_pe181_c -3 F_LMPD_182_pe181_c +3 R_LMPD_182_pe181_c -5 F_LMPD_183_pe181_c +5 R_LMPD_183_pe181_c -7 F_LMPD_184_pe181_c +7 R_LMPD_184_pe181_c -9 F_LMPD_185_pe181_c +9 R_LMPD_185_pe181_c + F_LMPD_186_pe181_c - R_LMPD_186_pe181_c - F_LMPD_187_pe181_c + R_LMPD_187_pe181_c -15 F_LMPD_188_pe181_c +15 R_LMPD_188_pe181_c -29 F_LMPD_189_pe181_c +29 R_LMPD_189_pe181_c -10 F_LMPD_190_pg160_c +10 R_LMPD_190_pg160_c -12 F_LMPD_191_pg160_c +12 R_LMPD_191_pg160_c -14 F_LMPD_192_pg160_c +14 R_LMPD_192_pg160_c -2 F_LMPD_193_pg160_c +2 R_LMPD_193_pg160_c -4 F_LMPD_195_pg160_c +4 R_LMPD_195_pg160_c -6 F_LMPD_196_pg160_c +6 R_LMPD_196_pg160_c -8 F_LMPD_197_pg160_c +8 R_LMPD_197_pg160_c -14 F_LMPD_198_pg160_c +14 R_LMPD_198_pg160_c -28 F_LMPD_199_pg160_c +28 R_LMPD_199_pg160_c -10 F_LMPD_200_pg161_c +10 R_LMPD_200_pg161_c -12 F_LMPD_201_pg161_c +12 R_LMPD_201_pg161_c -2 F_LMPD_202_pg161_c +2 R_LMPD_202_pg161_c -4 F_LMPD_204_pg161_c +4 R_LMPD_204_pg161_c -6 F_LMPD_205_pg161_c +6 R_LMPD_205_pg161_c -8 F_LMPD_206_pg161_c +8 R_LMPD_206_pg161_c -10 F_LMPD_207_pg181_c +10 R_LMPD_207_pg181_c -12 F_LMPD_208_pg181_c +12 R_LMPD_208_pg181_c -14 F_LMPD_209_pg181_c +14 R_LMPD_209_pg181_c -2 F_LMPD_210_pg181_c +2 R_LMPD_210_pg181_c -4 F_LMPD_212_pg181_c +4 R_LMPD_212_pg181_c -6 F_LMPD_213_pg181_c +6 R_LMPD_213_pg181_c -8 F_LMPD_214_pg181_c +8 R_LMPD_214_pg181_c -16 F_LMPD_215_pg181_c +16 R_LMPD_215_pg181_c -30 F_LMPD_216_pg181_c +30 R_LMPD_216_pg181_c +3 F_LMPD_222_pydx5p_c -3 R_LMPD_222_pydx5p_c +4 F_LMPD_223_pydx5p_c -4 R_LMPD_223_pydx5p_c -8 F_LMPD_224_q8h2_c +8 R_LMPD_224_q8h2_c -7 F_LMPD_225_q8h2_c +7 R_LMPD_225_q8h2_c -5 F_LMPD_226_q8h2_c +5 R_LMPD_226_q8h2_c -11 F_LMPD_227_q8h2_c +11 R_LMPD_227_q8h2_c -10 F_LMPD_228_q8h2_c +10 R_LMPD_228_q8h2_c -10 F_LMPD_229_q8h2_c +10 R_LMPD_229_q8h2_c -8 F_LMPD_230_q8h2_c +8 R_LMPD_230_q8h2_c -8 F_LMPD_231_q8h2_c +8 R_LMPD_231_q8h2_c -12 F_LMPD_232_q8h2_c +12 R_LMPD_232_q8h2_c -18 F_LMPD_233_q8h2_c +18 R_LMPD_233_q8h2_c -18 F_LMPD_234_q8h2_c +18 R_LMPD_234_q8h2_c +3 F_LMPD_235_ribflv_c -3 R_LMPD_235_ribflv_c +3 F_LMPD_236_ribflv_c -3 R_LMPD_236_ribflv_c +2 F_LMPD_237_ribflv_c -2 R_LMPD_237_ribflv_c +2 F_LMPD_238_ribflv_c -2 R_LMPD_238_ribflv_c + F_LMPD_239_ser-L_c - R_LMPD_239_ser-L_c - F_LMPD_240_sheme_c + R_LMPD_240_sheme_c + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c +4 F_LMPD_243_thf_c -4 R_LMPD_243_thf_c +4 F_LMPD_244_thf_c -4 R_LMPD_244_thf_c +3 F_LMPD_245_thf_c -3 R_LMPD_245_thf_c +3 F_LMPD_246_thf_c -3 R_LMPD_246_thf_c +4 F_LMPD_247_thmpp_c -4 R_LMPD_247_thmpp_c +4 F_LMPD_248_thmpp_c -4 R_LMPD_248_thmpp_c + F_LMPD_251_tyr-L_c - R_LMPD_251_tyr-L_c -11 F_LMPD_252_udcpdp_c +11 R_LMPD_252_udcpdp_c = 0
M_nadp_c: -0.000112 F_Ec_biomass_iJO1366_WT_53p95M +0.000112 R_Ec_biomass_iJO1366_WT_53p95M - F_ALDD2y + R_ALDD2y + F_FLDR2 - R_FLDR2 - F_G3PD2 + R_G3PD2 - F_G6PDH2r + R_G6PDH2r - F_GLUDy + R_GLUDy + F_GLYCLTDy - R_GLYCLTDy - F_GND + R_GND - F_ICDHyr + R_ICDHyr - F_ME2 + R_ME2 + F_NADK - R_NADK + F_NADPHQR2 - R_NADPHQR2 + F_NADPHQR3 - R_NADPHQR3 + F_NADPHQR4 - R_NADPHQR4 - F_NADPPPS + R_NADPPPS + F_NADTRHD - R_NADTRHD - F_THD2pp + R_THD2pp + F_TRDR - R_TRDR - F_LMPD_3_10fthf_c + R_LMPD_3_10fthf_c - F_LMPD_4_10fthf_c + R_LMPD_4_10fthf_c +9 F_LMPD_5_2dmmql8_c -9 R_LMPD_5_2dmmql8_c +7 F_LMPD_6_2fe2s_c -7 R_LMPD_6_2fe2s_c +6 F_LMPD_7_2fe2s_c -6 R_LMPD_7_2fe2s_c +7 F_LMPD_8_2fe2s_c -7 R_LMPD_8_2fe2s_c +6 F_LMPD_9_2fe2s_c -6 R_LMPD_9_2fe2s_c +7 F_LMPD_10_2fe2s_c -7 R_LMPD_10_2fe2s_c +6 F_LMPD_11_2fe2s_c -6 R_LMPD_11_2fe2s_c +7 F_LMPD_12_2fe2s_c -7 R_LMPD_12_2fe2s_c +6 F_LMPD_13_2fe2s_c -6 R_LMPD_13_2fe2s_c +15 F_LMPD_14_4fe4s_c -15 R_LMPD_14_4fe4s_c +12 F_LMPD_15_4fe4s_c -12 R_LMPD_15_4fe4s_c +15 F_LMPD_16_4fe4s_c -15 R_LMPD_16_4fe4s_c +12 F_LMPD_17_4fe4s_c -12 R_LMPD_17_4fe4s_c +15 F_LMPD_18_4fe4s_c -15 R_LMPD_18_4fe4s_c +12 F_LMPD_19_4fe4s_c -12 R_LMPD_19_4fe4s_c +15 F_LMPD_20_4fe4s_c -15 R_LMPD_20_4fe4s_c +12 F_LMPD_21_4fe4s_c -12 R_LMPD_21_4fe4s_c +2 F_LMPD_22_5mthf_c -2 R_LMPD_22_5mthf_c +2 F_LMPD_23_5mthf_c -2 R_LMPD_23_5mthf_c +5 F_LMPD_24_accoa_c -5 R_LMPD_24_accoa_c +5 F_LMPD_25_accoa_c -5 R_LMPD_25_accoa_c -2 F_LMPD_26_adocbl_c +2 R_LMPD_26_adocbl_c -2 F_LMPD_27_adocbl_c +2 R_LMPD_27_adocbl_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c +5 F_LMPD_31_amet_c -5 R_LMPD_31_amet_c +5 F_LMPD_32_amet_c -5 R_LMPD_32_amet_c + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c +4 F_LMPD_35_bmocogdp_c -4 R_LMPD_35_bmocogdp_c +4 F_LMPD_36_bmocogdp_c -4 R_LMPD_36_bmocogdp_c +8 F_LMPD_37_bmocogdp_c -8 R_LMPD_37_bmocogdp_c +8 F_LMPD_38_bmocogdp_c -8 R_LMPD_38_bmocogdp_c +12 F_LMPD_39_btn_c -12 R_LMPD_39_btn_c +12 F_LMPD_40_btn_c -12 R_LMPD_40_btn_c +12 F_LMPD_41_btn_c -12 R_LMPD_41_btn_c +12 F_LMPD_42_btn_c -12 R_LMPD_42_btn_c +17 F_LMPD_43_btn_c -17 R_LMPD_43_btn_c +17 F_LMPD_44_btn_c -17 R_LMPD_44_btn_c +17 F_LMPD_45_btn_c -17 R_LMPD_45_btn_c +17 F_LMPD_47_btn_c -17 R_LMPD_47_btn_c + F_LMPD_51_chor_c - R_LMPD_51_chor_c +44 F_LMPD_52_clpn160_p -44 R_LMPD_52_clpn160_p +40 F_LMPD_53_clpn160_p -40 R_LMPD_53_clpn160_p +36 F_LMPD_54_clpn160_p -36 R_LMPD_54_clpn160_p +32 F_LMPD_55_clpn160_p -32 R_LMPD_55_clpn160_p +28 F_LMPD_56_clpn160_p -28 R_LMPD_56_clpn160_p +52 F_LMPD_57_clpn160_p -52 R_LMPD_57_clpn160_p +56 F_LMPD_58_clpn160_p -56 R_LMPD_58_clpn160_p +48 F_LMPD_59_clpn160_p -48 R_LMPD_59_clpn160_p +28 F_LMPD_60_clpn160_p -28 R_LMPD_60_clpn160_p +40 F_LMPD_62_clpn161_p -40 R_LMPD_62_clpn161_p +36 F_LMPD_63_clpn161_p -36 R_LMPD_63_clpn161_p +32 F_LMPD_64_clpn161_p -32 R_LMPD_64_clpn161_p +28 F_LMPD_65_clpn161_p -28 R_LMPD_65_clpn161_p +48 F_LMPD_66_clpn161_p -48 R_LMPD_66_clpn161_p +52 F_LMPD_67_clpn161_p -52 R_LMPD_67_clpn161_p +44 F_LMPD_68_clpn161_p -44 R_LMPD_68_clpn161_p +48 F_LMPD_69_clpn181_p -48 R_LMPD_69_clpn181_p +44 F_LMPD_70_clpn181_p -44 R_LMPD_70_clpn181_p +40 F_LMPD_71_clpn181_p -40 R_LMPD_71_clpn181_p +36 F_LMPD_72_clpn181_p -36 R_LMPD_72_clpn181_p +32 F_LMPD_73_clpn181_p -32 R_LMPD_73_clpn181_p +56 F_LMPD_74_clpn181_p -56 R_LMPD_74_clpn181_p +60 F_LMPD_75_clpn181_p -60 R_LMPD_75_clpn181_p +52 F_LMPD_76_clpn181_p -52 R_LMPD_76_clpn181_p +28 F_LMPD_77_clpn181_p -28 R_LMPD_77_clpn181_p +55 F_LMPD_79_colipa_e -55 R_LMPD_79_colipa_e +54 F_LMPD_80_colipa_e -54 R_LMPD_80_colipa_e +49 F_LMPD_81_colipa_e -49 R_LMPD_81_colipa_e +48 F_LMPD_82_colipa_e -48 R_LMPD_82_colipa_e +43 F_LMPD_83_colipa_e -43 R_LMPD_83_colipa_e +42 F_LMPD_84_colipa_e -42 R_LMPD_84_colipa_e +37 F_LMPD_85_colipa_e -37 R_LMPD_85_colipa_e +36 F_LMPD_86_colipa_e -36 R_LMPD_86_colipa_e +61 F_LMPD_87_colipa_e -61 R_LMPD_87_colipa_e +67 F_LMPD_88_colipa_e -67 R_LMPD_88_colipa_e +60 F_LMPD_89_colipa_e -60 R_LMPD_89_colipa_e +66 F_LMPD_90_colipa_e -66 R_LMPD_90_colipa_e +3 F_LMPD_92_cys-L_c -3 R_LMPD_92_cys-L_c +3 F_LMPD_93_cys-L_c -3 R_LMPD_93_cys-L_c -2 F_LMPD_94_datp_c +2 R_LMPD_94_datp_c -2 F_LMPD_95_datp_c +2 R_LMPD_95_datp_c - F_LMPD_96_datp_c + R_LMPD_96_datp_c - F_LMPD_97_datp_c + R_LMPD_97_datp_c -2 F_LMPD_99_dgtp_c +2 R_LMPD_99_dgtp_c -2 F_LMPD_100_dgtp_c +2 R_LMPD_100_dgtp_c -2 F_LMPD_101_dgtp_c +2 R_LMPD_101_dgtp_c - F_LMPD_102_dgtp_c + R_LMPD_102_dgtp_c - F_LMPD_103_dgtp_c + R_LMPD_103_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c +2 F_LMPD_105_dttp_c -2 R_LMPD_105_dttp_c +2 F_LMPD_106_dttp_c -2 R_LMPD_106_dttp_c +3 F_LMPD_107_enter_c -3 R_LMPD_107_enter_c - F_LMPD_108_fad_c + R_LMPD_108_fad_c - F_LMPD_109_fad_c + R_LMPD_109_fad_c -3 F_LMPD_110_fad_c +3 R_LMPD_110_fad_c -3 F_LMPD_111_fad_c +3 R_LMPD_111_fad_c -2 F_LMPD_112_GAM_c +2 R_LMPD_112_GAM_c -2 F_LMPD_113_GAM_c +2 R_LMPD_113_GAM_c - F_LMPD_114_GAM_c + R_LMPD_114_GAM_c - F_LMPD_115_GAM_c + R_LMPD_115_GAM_c +2 F_LMPD_117_gly_c -2 R_LMPD_117_gly_c +2 F_LMPD_118_gthrd_c -2 R_LMPD_118_gthrd_c +2 F_LMPD_119_gthrd_c -2 R_LMPD_119_gthrd_c -2 F_LMPD_120_gtp_c +2 R_LMPD_120_gtp_c -2 F_LMPD_121_gtp_c +2 R_LMPD_121_gtp_c - F_LMPD_122_gtp_c + R_LMPD_122_gtp_c - F_LMPD_123_gtp_c + R_LMPD_123_gtp_c +11 F_LMPD_124_hemeO_c -11 R_LMPD_124_hemeO_c +11 F_LMPD_125_hemeO_c -11 R_LMPD_125_hemeO_c +3 F_LMPD_128_ile-L_c -3 R_LMPD_128_ile-L_c + F_LMPD_129_leu-L_c - R_LMPD_129_leu-L_c +12 F_LMPD_130_lipopb_c -12 R_LMPD_130_lipopb_c +12 F_LMPD_131_lipopb_c -12 R_LMPD_131_lipopb_c +12 F_LMPD_132_lipopb_c -12 R_LMPD_132_lipopb_c +18 F_LMPD_133_lipopb_c -18 R_LMPD_133_lipopb_c +2 F_LMPD_134_lys-L_c -2 R_LMPD_134_lys-L_c +5 F_LMPD_135_malcoa_c -5 R_LMPD_135_malcoa_c +5 F_LMPD_136_malcoa_c -5 R_LMPD_136_malcoa_c +5 F_LMPD_137_met-L_c -5 R_LMPD_137_met-L_c +5 F_LMPD_138_met-L_c -5 R_LMPD_138_met-L_c +2 F_LMPD_139_mlthf_c -2 R_LMPD_139_mlthf_c +2 F_LMPD_140_mlthf_c -2 R_LMPD_140_mlthf_c +5 F_LMPD_141_mococdp_c -5 R_LMPD_141_mococdp_c +5 F_LMPD_142_mococdp_c -5 R_LMPD_142_mococdp_c +5 F_LMPD_143_mococdp_c -5 R_LMPD_143_mococdp_c +5 F_LMPD_144_mococdp_c -5 R_LMPD_144_mococdp_c +4 F_LMPD_145_mococdp_c -4 R_LMPD_145_mococdp_c +4 F_LMPD_146_mococdp_c -4 R_LMPD_146_mococdp_c +4 F_LMPD_147_mococdp_c -4 R_LMPD_147_mococdp_c +4 F_LMPD_148_mococdp_c -4 R_LMPD_148_mococdp_c +4 F_LMPD_149_mocogdp_c -4 R_LMPD_149_mocogdp_c +4 F_LMPD_150_mocogdp_c -4 R_LMPD_150_mocogdp_c +2 F_LMPD_151_mocogdp_c -2 R_LMPD_151_mocogdp_c +2 F_LMPD_152_mocogdp_c -2 R_LMPD_152_mocogdp_c +10 F_LMPD_153_mql8_c -10 R_LMPD_153_mql8_c +6 F_LMPD_154_murein3p3p_p -6 R_LMPD_154_murein3p3p_p +6 F_LMPD_155_murein3px4p_p -6 R_LMPD_155_murein3px4p_p +6 F_LMPD_156_murein4p4p_p -6 R_LMPD_156_murein4p4p_p +6 F_LMPD_157_murein4px4p_p -6 R_LMPD_157_murein4px4p_p +9 F_LMPD_158_murein4px4px4p_p -9 R_LMPD_158_murein4px4px4p_p -2 F_LMPD_159_nad_c +2 R_LMPD_159_nad_c -2 F_LMPD_160_nad_c +2 R_LMPD_160_nad_c - F_LMPD_161_nad_c + R_LMPD_161_nad_c - F_LMPD_162_nad_c + R_LMPD_162_nad_c +16 F_LMPD_163_pe160_c -16 R_LMPD_163_pe160_c +14 F_LMPD_164_pe160_c -14 R_LMPD_164_pe160_c +24 F_LMPD_165_pe160_c -24 R_LMPD_165_pe160_c +22 F_LMPD_166_pe160_c -22 R_LMPD_166_pe160_c +20 F_LMPD_167_pe160_c -20 R_LMPD_167_pe160_c +18 F_LMPD_168_pe160_c -18 R_LMPD_168_pe160_c +28 F_LMPD_169_pe160_c -28 R_LMPD_169_pe160_c +26 F_LMPD_170_pe160_c -26 R_LMPD_170_pe160_c +14 F_LMPD_171_pe160_c -14 R_LMPD_171_pe160_c +14 F_LMPD_173_pe161_c -14 R_LMPD_173_pe161_c +22 F_LMPD_174_pe161_c -22 R_LMPD_174_pe161_c +20 F_LMPD_175_pe161_c -20 R_LMPD_175_pe161_c +18 F_LMPD_176_pe161_c -18 R_LMPD_176_pe161_c +16 F_LMPD_177_pe161_c -16 R_LMPD_177_pe161_c +26 F_LMPD_178_pe161_c -26 R_LMPD_178_pe161_c +24 F_LMPD_179_pe161_c -24 R_LMPD_179_pe161_c +18 F_LMPD_180_pe181_c -18 R_LMPD_180_pe181_c +16 F_LMPD_181_pe181_c -16 R_LMPD_181_pe181_c +26 F_LMPD_182_pe181_c -26 R_LMPD_182_pe181_c +24 F_LMPD_183_pe181_c -24 R_LMPD_183_pe181_c +22 F_LMPD_184_pe181_c -22 R_LMPD_184_pe181_c +20 F_LMPD_185_pe181_c -20 R_LMPD_185_pe181_c +30 F_LMPD_186_pe181_c -30 R_LMPD_186_pe181_c +28 F_LMPD_187_pe181_c -28 R_LMPD_187_pe181_c +14 F_LMPD_188_pe181_c -14 R_LMPD_188_pe181_c +18 F_LMPD_190_pg160_c -18 R_LMPD_190_pg160_c +16 F_LMPD_191_pg160_c -16 R_LMPD_191_pg160_c +14 F_LMPD_192_pg160_c -14 R_LMPD_192_pg160_c +26 F_LMPD_193_pg160_c -26 R_LMPD_193_pg160_c +28 F_LMPD_194_pg160_c -28 R_LMPD_194_pg160_c +24 F_LMPD_195_pg160_c -24 R_LMPD_195_pg160_c +22 F_LMPD_196_pg160_c -22 R_LMPD_196_pg160_c +20 F_LMPD_197_pg160_c -20 R_LMPD_197_pg160_c +14 F_LMPD_198_pg160_c -14 R_LMPD_198_pg160_c +16 F_LMPD_200_pg161_c -16 R_LMPD_200_pg161_c +14 F_LMPD_201_pg161_c -14 R_LMPD_201_pg161_c +24 F_LMPD_202_pg161_c -24 R_LMPD_202_pg161_c +26 F_LMPD_203_pg161_c -26 R_LMPD_203_pg161_c +22 F_LMPD_204_pg161_c -22 R_LMPD_204_pg161_c +20 F_LMPD_205_pg161_c -20 R_LMPD_205_pg161_c +18 F_LMPD_206_pg161_c -18 R_LMPD_206_pg161_c +20 F_LMPD_207_pg181_c -20 R_LMPD_207_pg181_c +18 F_LMPD_208_pg181_c -18 R_LMPD_208_pg181_c +16 F_LMPD_209_pg181_c -16 R_LMPD_209_pg181_c +28 F_LMPD_210_pg181_c -28 R_LMPD_210_pg181_c +30 F_LMPD_211_pg181_c -30 R_LMPD_211_pg181_c +26 F_LMPD_212_pg181_c -26 R_LMPD_212_pg181_c +24 F_LMPD_213_pg181_c -24 R_LMPD_213_pg181_c +22 F_LMPD_214_pg181_c -22 R_LMPD_214_pg181_c +14 F_LMPD_215_pg181_c -14 R_LMPD_215_pg181_c + F_LMPD_217_phe-L_c - R_LMPD_217_phe-L_c +8 F_LMPD_218_pheme_c -8 R_LMPD_218_pheme_c +8 F_LMPD_219_pheme_c -8 R_LMPD_219_pheme_c +2 F_LMPD_220_pro-L_c -2 R_LMPD_220_pro-L_c + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c +9 F_LMPD_224_q8h2_c -9 R_LMPD_224_q8h2_c +9 F_LMPD_225_q8h2_c -9 R_LMPD_225_q8h2_c +9 F_LMPD_226_q8h2_c -9 R_LMPD_226_q8h2_c +12 F_LMPD_227_q8h2_c -12 R_LMPD_227_q8h2_c +12 F_LMPD_228_q8h2_c -12 R_LMPD_228_q8h2_c +12 F_LMPD_229_q8h2_c -12 R_LMPD_229_q8h2_c +12 F_LMPD_230_q8h2_c -12 R_LMPD_230_q8h2_c +12 F_LMPD_231_q8h2_c -12 R_LMPD_231_q8h2_c +18 F_LMPD_232_q8h2_c -18 R_LMPD_232_q8h2_c +24 F_LMPD_233_q8h2_c -24 R_LMPD_233_q8h2_c +24 F_LMPD_234_q8h2_c -24 R_LMPD_234_q8h2_c - F_LMPD_235_ribflv_c + R_LMPD_235_ribflv_c - F_LMPD_236_ribflv_c + R_LMPD_236_ribflv_c +10 F_LMPD_240_sheme_c -10 R_LMPD_240_sheme_c +6 F_LMPD_241_spmd_c -6 R_LMPD_241_spmd_c +6 F_LMPD_242_spmd_c -6 R_LMPD_242_spmd_c + F_LMPD_245_thf_c - R_LMPD_245_thf_c + F_LMPD_246_thf_c - R_LMPD_246_thf_c +5 F_LMPD_247_thmpp_c -5 R_LMPD_247_thmpp_c +5 F_LMPD_248_thmpp_c -5 R_LMPD_248_thmpp_c +2 F_LMPD_249_thr-L_c -2 R_LMPD_249_thr-L_c + F_LMPD_250_trp-L_c - R_LMPD_250_trp-L_c + F_LMPD_251_tyr-L_c - R_LMPD_251_tyr-L_c +11 F_LMPD_252_udcpdp_c -11 R_LMPD_252_udcpdp_c + F_LMPD_254_val-L_c - R_LMPD_254_val-L_c = 0
M_nadph_c: -0.000335 F_Ec_biomass_iJO1366_WT_53p95M +0.000335 R_Ec_biomass_iJO1366_WT_53p95M + F_ALDD2y - R_ALDD2y - F_FLDR2 + R_FLDR2 + F_G3PD2 - R_G3PD2 + F_G6PDH2r - R_G6PDH2r + F_GLUDy - R_GLUDy - F_GLYCLTDy + R_GLYCLTDy + F_GND - R_GND + F_ICDHyr - R_ICDHyr + F_ME2 - R_ME2 - F_NADPHQR2 + R_NADPHQR2 - F_NADPHQR3 + R_NADPHQR3 - F_NADPHQR4 + R_NADPHQR4 - F_NADTRHD + R_NADTRHD + F_THD2pp - R_THD2pp - F_TRDR + R_TRDR + F_LMPD_3_10fthf_c - R_LMPD_3_10fthf_c + F_LMPD_4_10fthf_c - R_LMPD_4_10fthf_c -9 F_LMPD_5_2dmmql8_c +9 R_LMPD_5_2dmmql8_c -7 F_LMPD_6_2fe2s_c +7 R_LMPD_6_2fe2s_c -6 F_LMPD_7_2fe2s_c +6 R_LMPD_7_2fe2s_c -7 F_LMPD_8_2fe2s_c +7 R_LMPD_8_2fe2s_c -6 F_LMPD_9_2fe2s_c +6 R_LMPD_9_2fe2s_c -7 F_LMPD_10_2fe2s_c +7 R_LMPD_10_2fe2s_c -6 F_LMPD_11_2fe2s_c +6 R_LMPD_11_2fe2s_c -7 F_LMPD_12_2fe2s_c +7 R_LMPD_12_2fe2s_c -6 F_LMPD_13_2fe2s_c +6 R_LMPD_13_2fe2s_c -15 F_LMPD_14_4fe4s_c +15 R_LMPD_14_4fe4s_c -12 F_LMPD_15_4fe4s_c +12 R_LMPD_15_4fe4s_c -15 F_LMPD_16_4fe4s_c +15 R_LMPD_16_4fe4s_c -12 F_LMPD_17_4fe4s_c +12 R_LMPD_17_4fe4s_c -15 F_LMPD_18_4fe4s_c +15 R_LMPD_18_4fe4s_c -12 F_LMPD_19_4fe4s_c +12 R_LMPD_19_4fe4s_c -15 F_LMPD_20_4fe4s_c +15 R_LMPD_20_4fe4s_c -12 F_LMPD_21_4fe4s_c +12 R_LMPD_21_4fe4s_c -2 F_LMPD_22_5mthf_c +2 R_LMPD_22_5mthf_c -2 F_LMPD_23_5mthf_c +2 R_LMPD_23_5mthf_c -5 F_LMPD_24_accoa_c +5 R_LMPD_24_accoa_c -5 F_LMPD_25_accoa_c +5 R_LMPD_25_accoa_c +2 F_LMPD_26_adocbl_c -2 R_LMPD_26_adocbl_c +2 F_LMPD_27_adocbl_c -2 R_LMPD_27_adocbl_c + F_LMPD_28_adocbl_c - R_LMPD_28_adocbl_c + F_LMPD_29_adocbl_c - R_LMPD_29_adocbl_c -5 F_LMPD_31_amet_c +5 R_LMPD_31_amet_c -5 F_LMPD_32_amet_c +5 R_LMPD_32_amet_c - F_LMPD_33_arg-L_c + R_LMPD_33_arg-L_c -4 F_LMPD_35_bmocogdp_c +4 R_LMPD_35_bmocogdp_c -4 F_LMPD_36_bmocogdp_c +4 R_LMPD_36_bmocogdp_c -8 F_LMPD_37_bmocogdp_c +8 R_LMPD_37_bmocogdp_c -8 F_LMPD_38_bmocogdp_c +8 R_LMPD_38_bmocogdp_c -12 F_LMPD_39_btn_c +12 R_LMPD_39_btn_c -12 F_LMPD_40_btn_c +12 R_LMPD_40_btn_c -12 F_LMPD_41_btn_c +12 R_LMPD_41_btn_c -12 F_LMPD_42_btn_c +12 R_LMPD_42_btn_c -17 F_LMPD_43_btn_c +17 R_LMPD_43_btn_c -17 F_LMPD_44_btn_c +17 R_LMPD_44_btn_c -17 F_LMPD_45_btn_c +17 R_LMPD_45_btn_c -17 F_LMPD_47_btn_c +17 R_LMPD_47_btn_c - F_LMPD_51_chor_c + R_LMPD_51_chor_c -44 F_LMPD_52_clpn160_p +44 R_LMPD_52_clpn160_p -40 F_LMPD_53_clpn160_p +40 R_LMPD_53_clpn160_p -36 F_LMPD_54_clpn160_p +36 R_LMPD_54_clpn160_p -32 F_LMPD_55_clpn160_p +32 R_LMPD_55_clpn160_p -28 F_LMPD_56_clpn160_p +28 R_LMPD_56_clpn160_p -52 F_LMPD_57_clpn160_p +52 R_LMPD_57_clpn160_p -56 F_LMPD_58_clpn160_p +56 R_LMPD_58_clpn160_p -48 F_LMPD_59_clpn160_p +48 R_LMPD_59_clpn160_p -28 F_LMPD_60_clpn160_p +28 R_LMPD_60_clpn160_p -40 F_LMPD_62_clpn161_p +40 R_LMPD_62_clpn161_p -36 F_LMPD_63_clpn161_p +36 R_LMPD_63_clpn161_p -32 F_LMPD_64_clpn161_p +32 R_LMPD_64_clpn161_p -28 F_LMPD_65_clpn161_p +28 R_LMPD_65_clpn161_p -48 F_LMPD_66_clpn161_p +48 R_LMPD_66_clpn161_p -52 F_LMPD_67_clpn161_p +52 R_LMPD_67_clpn161_p -44 F_LMPD_68_clpn161_p +44 R_LMPD_68_clpn161_p -48 F_LMPD_69_clpn181_p +48 R_LMPD_69_clpn181_p -44 F_LMPD_70_clpn181_p +44 R_LMPD_70_clpn181_p -40 F_LMPD_71_clpn181_p +40 R_LMPD_71_clpn181_p -36 F_LMPD_72_clpn181_p +36 R_LMPD_72_clpn181_p -32 F_LMPD_73_clpn181_p +32 R_LMPD_73_clpn181_p -56 F_LMPD_74_clpn181_p +56 R_LMPD_74_clpn181_p -60 F_LMPD_75_clpn181_p +60 R_LMPD_75_clpn181_p -52 F_LMPD_76_clpn181_p +52 R_LMPD_76_clpn181_p -28 F_LMPD_77_clpn181_p +28 R_LMPD_77_clpn181_p -55 F_LMPD_79_colipa_e +55 R_LMPD_79_colipa_e -54 F_LMPD_80_colipa_e +54 R_LMPD_80_colipa_e -49 F_LMPD_81_colipa_e +49 R_LMPD_81_colipa_e -48 F_LMPD_82_colipa_e +48 R_LMPD_82_colipa_e -43 F_LMPD_83_colipa_e +43 R_LMPD_83_colipa_e -42 F_LMPD_84_colipa_e +42 R_LMPD_84_colipa_e -37 F_LMPD_85_colipa_e +37 R_LMPD_85_colipa_e -36 F_LMPD_86_colipa_e +36 R_LMPD_86_colipa_e -61 F_LMPD_87_colipa_e +61 R_LMPD_87_colipa_e -67 F_LMPD_88_colipa_e +67 R_LMPD_88_colipa_e -60 F_LMPD_89_colipa_e +60 R_LMPD_89_colipa_e -66 F_LMPD_90_colipa_e +66 R_LMPD_90_colipa_e -3 F_LMPD_92_cys-L_c +3 R_LMPD_92_cys-L_c -3 F_LMPD_93_cys-L_c +3 R_LMPD_93_cys-L_c +2 F_LMPD_94_datp_c -2 R_LMPD_94_datp_c +2 F_LMPD_95_datp_c -2 R_LMPD_95_datp_c + F_LMPD_96_datp_c - R_LMPD_96_datp_c + F_LMPD_97_datp_c - R_LMPD_97_datp_c +2 F_LMPD_99_dgtp_c -2 R_LMPD_99_dgtp_c +2 F_LMPD_100_dgtp_c -2 R_LMPD_100_dgtp_c +2 F_LMPD_101_dgtp_c -2 R_LMPD_101_dgtp_c + F_LMPD_102_dgtp_c - R_LMPD_102_dgtp_c + F_LMPD_103_dgtp_c - R_LMPD_103_dgtp_c + F_LMPD_104_dgtp_c - R_LMPD_104_dgtp_c -2 F_LMPD_105_dttp_c +2 R_LMPD_105_dttp_c -2 F_LMPD_106_dttp_c +2 R_LMPD_106_dttp_c -3 F_LMPD_107_enter_c +3 R_LMPD_107_enter_c + F_LMPD_108_fad_c - R_LMPD_108_fad_c + F_LMPD_109_fad_c - R_LMPD_109_fad_c +3 F_LMPD_110_fad_c -3 R_LMPD_110_fad_c +3 F_LMPD_111_fad_c -3 R_LMPD_111_fad_c +2 F_LMPD_112_GAM_c -2 R_LMPD_112_GAM_c +2 F_LMPD_113_GAM_c -2 R_LMPD_113_GAM_c + F_LMPD_114_GAM_c - R_LMPD_114_GAM_c + F_LMPD_115_GAM_c - R_LMPD_115_GAM_c -2 F_LMPD_117_gly_c +2 R_LMPD_117_gly_c -2 F_LMPD_118_gthrd_c +2 R_LMPD_118_gthrd_c -2 F_LMPD_119_gthrd_c +2 R_LMPD_119_gthrd_c +2 F_LMPD_120_gtp_c -2 R_LMPD_120_gtp_c +2 F_LMPD_121_gtp_c -2 R_LMPD_121_gtp_c + F_LMPD_122_gtp_c - R_LMPD_122_gtp_c + F_LMPD_123_gtp_c - R_LMPD_123_gtp_c -11 F_LMPD_124_hemeO_c +11 R_LMPD_124_hemeO_c -11 F_LMPD_125_hemeO_c +11 R_LMPD_125_hemeO_c -3 F_LMPD_128_ile-L_c +3 R_LMPD_128_ile-L_c - F_LMPD_129_leu-L_c + R_LMPD_129_leu-L_c -12 F_LMPD_130_lipopb_c +12 R_LMPD_130_lipopb_c -12 F_LMPD_131_lipopb_c +12 R_LMPD_131_lipopb_c -12 F_LMPD_132_lipopb_c +12 R_LMPD_132_lipopb_c -18 F_LMPD_133_lipopb_c +18 R_LMPD_133_lipopb_c -2 F_LMPD_134_lys-L_c +2 R_LMPD_134_lys-L_c -5 F_LMPD_135_malcoa_c +5 R_LMPD_135_malcoa_c -5 F_LMPD_136_malcoa_c +5 R_LMPD_136_malcoa_c -5 F_LMPD_137_met-L_c +5 R_LMPD_137_met-L_c -5 F_LMPD_138_met-L_c +5 R_LMPD_138_met-L_c -2 F_LMPD_139_mlthf_c +2 R_LMPD_139_mlthf_c -2 F_LMPD_140_mlthf_c +2 R_LMPD_140_mlthf_c -5 F_LMPD_141_mococdp_c +5 R_LMPD_141_mococdp_c -5 F_LMPD_142_mococdp_c +5 R_LMPD_142_mococdp_c -5 F_LMPD_143_mococdp_c +5 R_LMPD_143_mococdp_c -5 F_LMPD_144_mococdp_c +5 R_LMPD_144_mococdp_c -4 F_LMPD_145_mococdp_c +4 R_LMPD_145_mococdp_c -4 F_LMPD_146_mococdp_c +4 R_LMPD_146_mococdp_c -4 F_LMPD_147_mococdp_c +4 R_LMPD_147_mococdp_c -4 F_LMPD_148_mococdp_c +4 R_LMPD_148_mococdp_c -4 F_LMPD_149_mocogdp_c +4 R_LMPD_149_mocogdp_c -4 F_LMPD_150_mocogdp_c +4 R_LMPD_150_mocogdp_c -2 F_LMPD_151_mocogdp_c +2 R_LMPD_151_mocogdp_c -2 F_LMPD_152_mocogdp_c +2 R_LMPD_152_mocogdp_c -10 F_LMPD_153_mql8_c +10 R_LMPD_153_mql8_c -6 F_LMPD_154_murein3p3p_p +6 R_LMPD_154_murein3p3p_p -6 F_LMPD_155_murein3px4p_p +6 R_LMPD_155_murein3px4p_p -6 F_LMPD_156_murein4p4p_p +6 R_LMPD_156_murein4p4p_p -6 F_LMPD_157_murein4px4p_p +6 R_LMPD_157_murein4px4p_p -9 F_LMPD_158_murein4px4px4p_p +9 R_LMPD_158_murein4px4px4p_p +2 F_LMPD_159_nad_c -2 R_LMPD_159_nad_c +2 F_LMPD_160_nad_c -2 R_LMPD_160_nad_c + F_LMPD_161_nad_c - R_LMPD_161_nad_c + F_LMPD_162_nad_c - R_LMPD_162_nad_c -16 F_LMPD_163_pe160_c +16 R_LMPD_163_pe160_c -14 F_LMPD_164_pe160_c +14 R_LMPD_164_pe160_c -24 F_LMPD_165_pe160_c +24 R_LMPD_165_pe160_c -22 F_LMPD_166_pe160_c +22 R_LMPD_166_pe160_c -20 F_LMPD_167_pe160_c +20 R_LMPD_167_pe160_c -18 F_LMPD_168_pe160_c +18 R_LMPD_168_pe160_c -28 F_LMPD_169_pe160_c +28 R_LMPD_169_pe160_c -26 F_LMPD_170_pe160_c +26 R_LMPD_170_pe160_c -14 F_LMPD_171_pe160_c +14 R_LMPD_171_pe160_c -14 F_LMPD_173_pe161_c +14 R_LMPD_173_pe161_c -22 F_LMPD_174_pe161_c +22 R_LMPD_174_pe161_c -20 F_LMPD_175_pe161_c +20 R_LMPD_175_pe161_c -18 F_LMPD_176_pe161_c +18 R_LMPD_176_pe161_c -16 F_LMPD_177_pe161_c +16 R_LMPD_177_pe161_c -26 F_LMPD_178_pe161_c +26 R_LMPD_178_pe161_c -24 F_LMPD_179_pe161_c +24 R_LMPD_179_pe161_c -18 F_LMPD_180_pe181_c +18 R_LMPD_180_pe181_c -16 F_LMPD_181_pe181_c +16 R_LMPD_181_pe181_c -26 F_LMPD_182_pe181_c +26 R_LMPD_182_pe181_c -24 F_LMPD_183_pe181_c +24 R_LMPD_183_pe181_c -22 F_LMPD_184_pe181_c +22 R_LMPD_184_pe181_c -20 F_LMPD_185_pe181_c +20 R_LMPD_185_pe181_c -30 F_LMPD_186_pe181_c +30 R_LMPD_186_pe181_c -28 F_LMPD_187_pe181_c +28 R_LMPD_187_pe181_c -14 F_LMPD_188_pe181_c +14 R_LMPD_188_pe181_c -18 F_LMPD_190_pg160_c +18 R_LMPD_190_pg160_c -16 F_LMPD_191_pg160_c +16 R_LMPD_191_pg160_c -14 F_LMPD_192_pg160_c +14 R_LMPD_192_pg160_c -26 F_LMPD_193_pg160_c +26 R_LMPD_193_pg160_c -28 F_LMPD_194_pg160_c +28 R_LMPD_194_pg160_c -24 F_LMPD_195_pg160_c +24 R_LMPD_195_pg160_c -22 F_LMPD_196_pg160_c +22 R_LMPD_196_pg160_c -20 F_LMPD_197_pg160_c +20 R_LMPD_197_pg160_c -14 F_LMPD_198_pg160_c +14 R_LMPD_198_pg160_c -16 F_LMPD_200_pg161_c +16 R_LMPD_200_pg161_c -14 F_LMPD_201_pg161_c +14 R_LMPD_201_pg161_c -24 F_LMPD_202_pg161_c +24 R_LMPD_202_pg161_c -26 F_LMPD_203_pg161_c +26 R_LMPD_203_pg161_c -22 F_LMPD_204_pg161_c +22 R_LMPD_204_pg161_c -20 F_LMPD_205_pg161_c +20 R_LMPD_205_pg161_c -18 F_LMPD_206_pg161_c +18 R_LMPD_206_pg161_c -20 F_LMPD_207_pg181_c +20 R_LMPD_207_pg181_c -18 F_LMPD_208_pg181_c +18 R_LMPD_208_pg181_c -16 F_LMPD_209_pg181_c +16 R_LMPD_209_pg181_c -28 F_LMPD_210_pg181_c +28 R_LMPD_210_pg181_c -30 F_LMPD_211_pg181_c +30 R_LMPD_211_pg181_c -26 F_LMPD_212_pg181_c +26 R_LMPD_212_pg181_c -24 F_LMPD_213_pg181_c +24 R_LMPD_213_pg181_c -22 F_LMPD_214_pg181_c +22 R_LMPD_214_pg181_c -14 F_LMPD_215_pg181_c +14 R_LMPD_215_pg181_c - F_LMPD_217_phe-L_c + R_LMPD_217_phe-L_c -8 F_LMPD_218_pheme_c +8 R_LMPD_218_pheme_c -8 F_LMPD_219_pheme_c +8 R_LMPD_219_pheme_c -2 F_LMPD_220_pro-L_c +2 R_LMPD_220_pro-L_c - F_LMPD_221_ptrc_c + R_LMPD_221_ptrc_c -9 F_LMPD_224_q8h2_c +9 R_LMPD_224_q8h2_c -9 F_LMPD_225_q8h2_c +9 R_LMPD_225_q8h2_c -9 F_LMPD_226_q8h2_c +9 R_LMPD_226_q8h2_c -12 F_LMPD_227_q8h2_c +12 R_LMPD_227_q8h2_c -12 F_LMPD_228_q8h2_c +12 R_LMPD_228_q8h2_c -12 F_LMPD_229_q8h2_c +12 R_LMPD_229_q8h2_c -12 F_LMPD_230_q8h2_c +12 R_LMPD_230_q8h2_c -12 F_LMPD_231_q8h2_c +12 R_LMPD_231_q8h2_c -18 F_LMPD_232_q8h2_c +18 R_LMPD_232_q8h2_c -24 F_LMPD_233_q8h2_c +24 R_LMPD_233_q8h2_c -24 F_LMPD_234_q8h2_c +24 R_LMPD_234_q8h2_c + F_LMPD_235_ribflv_c - R_LMPD_235_ribflv_c + F_LMPD_236_ribflv_c - R_LMPD_236_ribflv_c -10 F_LMPD_240_sheme_c +10 R_LMPD_240_sheme_c -6 F_LMPD_241_spmd_c +6 R_LMPD_241_spmd_c -6 F_LMPD_242_spmd_c +6 R_LMPD_242_spmd_c - F_LMPD_245_thf_c + R_LMPD_245_thf_c - F_LMPD_246_thf_c + R_LMPD_246_thf_c -5 F_LMPD_247_thmpp_c +5 R_LMPD_247_thmpp_c -5 F_LMPD_248_thmpp_c +5 R_LMPD_248_thmpp_c -2 F_LMPD_249_thr-L_c +2 R_LMPD_249_thr-L_c - F_LMPD_250_trp-L_c + R_LMPD_250_trp-L_c - F_LMPD_251_tyr-L_c + R_LMPD_251_tyr-L_c -11 F_LMPD_252_udcpdp_c +11 R_LMPD_252_udcpdp_c - F_LMPD_254_val-L_c + R_LMPD_254_val-L_c = 0
M_nh4_c: -0.012379 F_Ec_biomass_iJO1366_WT_53p95M +0.012379 R_Ec_biomass_iJO1366_WT_53p95M + F_ASPT - R_ASPT + F_GLUDy - R_GLUDy + F_NH4tpp - R_NH4tpp -4 F_LMPD_1_10fthf_c +4 R_LMPD_1_10fthf_c -4 F_LMPD_2_10fthf_c +4 R_LMPD_2_10fthf_c -4 F_LMPD_3_10fthf_c +4 R_LMPD_3_10fthf_c -4 F_LMPD_4_10fthf_c +4 R_LMPD_4_10fthf_c -4 F_LMPD_22_5mthf_c +4 R_LMPD_22_5mthf_c -4 F_LMPD_23_5mthf_c +4 R_LMPD_23_5mthf_c -2 F_LMPD_24_accoa_c +2 R_LMPD_24_accoa_c -2 F_LMPD_25_accoa_c +2 R_LMPD_25_accoa_c -2 F_LMPD_26_adocbl_c +2 R_LMPD_26_adocbl_c -2 F_LMPD_27_adocbl_c +2 R_LMPD_27_adocbl_c -2 F_LMPD_28_adocbl_c +2 R_LMPD_28_adocbl_c -2 F_LMPD_29_adocbl_c +2 R_LMPD_29_adocbl_c - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c - F_LMPD_33_arg-L_c + R_LMPD_33_arg-L_c - F_LMPD_34_asn-L_c + R_LMPD_34_asn-L_c -12 F_LMPD_35_bmocogdp_c +12 R_LMPD_35_bmocogdp_c -12 F_LMPD_36_bmocogdp_c +12 R_LMPD_36_bmocogdp_c -12 F_LMPD_37_bmocogdp_c +12 R_LMPD_37_bmocogdp_c -12 F_LMPD_38_bmocogdp_c +12 R_LMPD_38_bmocogdp_c - F_LMPD_39_btn_c + R_LMPD_39_btn_c - F_LMPD_40_btn_c + R_LMPD_40_btn_c - F_LMPD_41_btn_c + R_LMPD_41_btn_c - F_LMPD_42_btn_c + R_LMPD_42_btn_c -4 F_LMPD_43_btn_c +4 R_LMPD_43_btn_c -4 F_LMPD_44_btn_c +4 R_LMPD_44_btn_c -4 F_LMPD_45_btn_c +4 R_LMPD_45_btn_c -4 F_LMPD_47_btn_c +4 R_LMPD_47_btn_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c -2 F_LMPD_79_colipa_e +2 R_LMPD_79_colipa_e -2 F_LMPD_80_colipa_e +2 R_LMPD_80_colipa_e -2 F_LMPD_81_colipa_e +2 R_LMPD_81_colipa_e -2 F_LMPD_82_colipa_e +2 R_LMPD_82_colipa_e -2 F_LMPD_83_colipa_e +2 R_LMPD_83_colipa_e -2 F_LMPD_84_colipa_e +2 R_LMPD_84_colipa_e -2 F_LMPD_85_colipa_e +2 R_LMPD_85_colipa_e -2 F_LMPD_86_colipa_e +2 R_LMPD_86_colipa_e -2 F_LMPD_87_colipa_e +2 R_LMPD_87_colipa_e -2 F_LMPD_88_colipa_e +2 R_LMPD_88_colipa_e -2 F_LMPD_89_colipa_e +2 R_LMPD_89_colipa_e -2 F_LMPD_90_colipa_e +2 R_LMPD_90_colipa_e -2 F_LMPD_91_ctp_c +2 R_LMPD_91_ctp_c -2 F_LMPD_94_datp_c +2 R_LMPD_94_datp_c -2 F_LMPD_95_datp_c +2 R_LMPD_95_datp_c -2 F_LMPD_96_datp_c +2 R_LMPD_96_datp_c -2 F_LMPD_97_datp_c +2 R_LMPD_97_datp_c -2 F_LMPD_98_dctp_c +2 R_LMPD_98_dctp_c -3 F_LMPD_99_dgtp_c +3 R_LMPD_99_dgtp_c -3 F_LMPD_100_dgtp_c +3 R_LMPD_100_dgtp_c -3 F_LMPD_101_dgtp_c +3 R_LMPD_101_dgtp_c -3 F_LMPD_102_dgtp_c +3 R_LMPD_102_dgtp_c -3 F_LMPD_103_dgtp_c +3 R_LMPD_103_dgtp_c -3 F_LMPD_104_dgtp_c +3 R_LMPD_104_dgtp_c - F_LMPD_105_dttp_c + R_LMPD_105_dttp_c - F_LMPD_106_dttp_c + R_LMPD_106_dttp_c -4 F_LMPD_108_fad_c +4 R_LMPD_108_fad_c -4 F_LMPD_109_fad_c +4 R_LMPD_109_fad_c -4 F_LMPD_110_fad_c +4 R_LMPD_110_fad_c -4 F_LMPD_111_fad_c +4 R_LMPD_111_fad_c -2 F_LMPD_112_GAM_c +2 R_LMPD_112_GAM_c -2 F_LMPD_113_GAM_c +2 R_LMPD_113_GAM_c -2 F_LMPD_114_GAM_c +2 R_LMPD_114_GAM_c -2 F_LMPD_115_GAM_c +2 R_LMPD_115_GAM_c - F_LMPD_116_gln-L_c + R_LMPD_116_gln-L_c -3 F_LMPD_120_gtp_c +3 R_LMPD_120_gtp_c -3 F_LMPD_121_gtp_c +3 R_LMPD_121_gtp_c -3 F_LMPD_122_gtp_c +3 R_LMPD_122_gtp_c -3 F_LMPD_123_gtp_c +3 R_LMPD_123_gtp_c +4 F_LMPD_124_hemeO_c -4 R_LMPD_124_hemeO_c +4 F_LMPD_125_hemeO_c -4 R_LMPD_125_hemeO_c - F_LMPD_126_his-L_c + R_LMPD_126_his-L_c - F_LMPD_127_his-L_c + R_LMPD_127_his-L_c + F_LMPD_128_ile-L_c - R_LMPD_128_ile-L_c -2 F_LMPD_135_malcoa_c +2 R_LMPD_135_malcoa_c -2 F_LMPD_136_malcoa_c +2 R_LMPD_136_malcoa_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c -4 F_LMPD_139_mlthf_c +4 R_LMPD_139_mlthf_c -4 F_LMPD_140_mlthf_c +4 R_LMPD_140_mlthf_c -5 F_LMPD_141_mococdp_c +5 R_LMPD_141_mococdp_c -5 F_LMPD_142_mococdp_c +5 R_LMPD_142_mococdp_c -5 F_LMPD_143_mococdp_c +5 R_LMPD_143_mococdp_c -5 F_LMPD_144_mococdp_c +5 R_LMPD_144_mococdp_c -5 F_LMPD_145_mococdp_c +5 R_LMPD_145_mococdp_c -5 F_LMPD_146_mococdp_c +5 R_LMPD_146_mococdp_c -5 F_LMPD_147_mococdp_c +5 R_LMPD_147_mococdp_c -5 F_LMPD_148_mococdp_c +5 R_LMPD_148_mococdp_c -6 F_LMPD_149_mocogdp_c +6 R_LMPD_149_mocogdp_c -6 F_LMPD_150_mocogdp_c +6 R_LMPD_150_mocogdp_c -6 F_LMPD_151_mocogdp_c +6 R_LMPD_151_mocogdp_c -6 F_LMPD_152_mocogdp_c +6 R_LMPD_152_mocogdp_c -4 F_LMPD_154_murein3p3p_p +4 R_LMPD_154_murein3p3p_p -4 F_LMPD_155_murein3px4p_p +4 R_LMPD_155_murein3px4p_p -4 F_LMPD_156_murein4p4p_p +4 R_LMPD_156_murein4p4p_p -4 F_LMPD_157_murein4px4p_p +4 R_LMPD_157_murein4px4p_p -6 F_LMPD_158_murein4px4px4p_p +6 R_LMPD_158_murein4px4px4p_p -3 F_LMPD_159_nad_c +3 R_LMPD_159_nad_c -3 F_LMPD_160_nad_c +3 R_LMPD_160_nad_c -3 F_LMPD_161_nad_c +3 R_LMPD_161_nad_c -3 F_LMPD_162_nad_c +3 R_LMPD_162_nad_c +4 F_LMPD_218_pheme_c -4 R_LMPD_218_pheme_c +4 F_LMPD_219_pheme_c -4 R_LMPD_219_pheme_c -2 F_LMPD_235_ribflv_c +2 R_LMPD_235_ribflv_c -2 F_LMPD_236_ribflv_c +2 R_LMPD_236_ribflv_c -2 F_LMPD_237_ribflv_c +2 R_LMPD_237_ribflv_c -2 F_LMPD_238_ribflv_c +2 R_LMPD_238_ribflv_c +4 F_LMPD_240_sheme_c -4 R_LMPD_240_sheme_c + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c -4 F_LMPD_243_thf_c +4 R_LMPD_243_thf_c -4 F_LMPD_244_thf_c +4 R_LMPD_244_thf_c -4 F_LMPD_245_thf_c +4 R_LMPD_245_thf_c -4 F_LMPD_246_thf_c +4 R_LMPD_246_thf_c -2 F_LMPD_247_thmpp_c +2 R_LMPD_247_thmpp_c -2 F_LMPD_248_thmpp_c +2 R_LMPD_248_thmpp_c -2 F_LMPD_250_trp-L_c +2 R_LMPD_250_trp-L_c - F_LMPD_253_utp_c + R_LMPD_253_utp_c = 0
M_ni2_c: -0.000307 F_Ec_biomass_iJO1366_WT_53p95M +0.000307 R_Ec_biomass_iJO1366_WT_53p95M + F_NI2tpp - R_NI2tpp + F_NI2uabcpp - R_NI2uabcpp = 0
M_o2_c: - F_ASPO6 + R_ASPO6 -0.5 F_CYTBD2pp +0.5 R_CYTBD2pp -0.5 F_CYTBDpp +0.5 R_CYTBDpp -0.5 F_CYTBO3_4pp +0.5 R_CYTBO3_4pp - F_MOX + R_MOX + F_O2tpp - R_O2tpp -1.5 F_PPPGO +1.5 R_PPPGO - F_LMPD_124_hemeO_c + R_LMPD_124_hemeO_c - F_LMPD_218_pheme_c + R_LMPD_218_pheme_c - F_LMPD_222_pydx5p_c + R_LMPD_222_pydx5p_c -1.5 F_LMPD_224_q8h2_c +1.5 R_LMPD_224_q8h2_c - F_LMPD_225_q8h2_c + R_LMPD_225_q8h2_c -1.5 F_LMPD_227_q8h2_c +1.5 R_LMPD_227_q8h2_c - F_LMPD_228_q8h2_c + R_LMPD_228_q8h2_c - F_LMPD_229_q8h2_c + R_LMPD_229_q8h2_c - F_LMPD_232_q8h2_c + R_LMPD_232_q8h2_c - F_LMPD_233_q8h2_c + R_LMPD_233_q8h2_c - F_LMPD_234_q8h2_c + R_LMPD_234_q8h2_c = 0
M_oaa_c: F_ASPTA - R_ASPTA + F_CITL - R_CITL - F_CS + R_CS + F_MDH - R_MDH + F_MDH2 - R_MDH2 + F_MDH3 - R_MDH3 + F_MOX - R_MOX - F_OAADC + R_OAADC + F_PPC - R_PPC - F_PPCK + R_PPCK = 0
M_orot_c: F_DHORD2 - R_DHORD2 + F_DHORD5 - R_DHORD5 + F_DHORDfum - R_DHORDfum - F_LMPD_91_ctp_c + R_LMPD_91_ctp_c - F_LMPD_98_dctp_c + R_LMPD_98_dctp_c - F_LMPD_105_dttp_c + R_LMPD_105_dttp_c - F_LMPD_106_dttp_c + R_LMPD_106_dttp_c - F_LMPD_141_mococdp_c + R_LMPD_141_mococdp_c - F_LMPD_142_mococdp_c + R_LMPD_142_mococdp_c - F_LMPD_143_mococdp_c + R_LMPD_143_mococdp_c - F_LMPD_144_mococdp_c + R_LMPD_144_mococdp_c - F_LMPD_145_mococdp_c + R_LMPD_145_mococdp_c - F_LMPD_146_mococdp_c + R_LMPD_146_mococdp_c - F_LMPD_147_mococdp_c + R_LMPD_147_mococdp_c - F_LMPD_148_mococdp_c + R_LMPD_148_mococdp_c - F_LMPD_253_utp_c + R_LMPD_253_utp_c = 0
M_pe160_c: -0.012366 F_Ec_biomass_iJO1366_WT_53p95M +0.012366 R_Ec_biomass_iJO1366_WT_53p95M - F_PE160abcpp + R_PE160abcpp + F_LMPD_163_pe160_c - R_LMPD_163_pe160_c + F_LMPD_164_pe160_c - R_LMPD_164_pe160_c + F_LMPD_165_pe160_c - R_LMPD_165_pe160_c + F_LMPD_166_pe160_c - R_LMPD_166_pe160_c + F_LMPD_167_pe160_c - R_LMPD_167_pe160_c + F_LMPD_168_pe160_c - R_LMPD_168_pe160_c + F_LMPD_169_pe160_c - R_LMPD_169_pe160_c + F_LMPD_170_pe160_c - R_LMPD_170_pe160_c + F_LMPD_171_pe160_c - R_LMPD_171_pe160_c + F_LMPD_172_pe160_c - R_LMPD_172_pe160_c = 0
M_pe161_c: -0.009618 F_Ec_biomass_iJO1366_WT_53p95M +0.009618 R_Ec_biomass_iJO1366_WT_53p95M - F_PE161abcpp + R_PE161abcpp + F_LMPD_173_pe161_c - R_LMPD_173_pe161_c + F_LMPD_174_pe161_c - R_LMPD_174_pe161_c + F_LMPD_175_pe161_c - R_LMPD_175_pe161_c + F_LMPD_176_pe161_c - R_LMPD_176_pe161_c + F_LMPD_177_pe161_c - R_LMPD_177_pe161_c + F_LMPD_178_pe161_c - R_LMPD_178_pe161_c + F_LMPD_179_pe161_c - R_LMPD_179_pe161_c = 0
M_pe181_c: -0.004957 F_Ec_biomass_iJO1366_WT_53p95M +0.004957 R_Ec_biomass_iJO1366_WT_53p95M - F_PE181abcpp + R_PE181abcpp + F_LMPD_180_pe181_c - R_LMPD_180_pe181_c + F_LMPD_181_pe181_c - R_LMPD_181_pe181_c + F_LMPD_182_pe181_c - R_LMPD_182_pe181_c + F_LMPD_183_pe181_c - R_LMPD_183_pe181_c + F_LMPD_184_pe181_c - R_LMPD_184_pe181_c + F_LMPD_185_pe181_c - R_LMPD_185_pe181_c + F_LMPD_186_pe181_c - R_LMPD_186_pe181_c + F_LMPD_187_pe181_c - R_LMPD_187_pe181_c + F_LMPD_188_pe181_c - R_LMPD_188_pe181_c + F_LMPD_189_pe181_c - R_LMPD_189_pe181_c = 0
M_pep_c: - F_DHAPT + R_DHAPT + F_ENO - R_ENO - F_GLCptspp + R_GLCptspp - F_PPC + R_PPC + F_PPCK - R_PPCK + F_PPS - R_PPS - F_PYK + R_PYK -2 F_LMPD_1_10fthf_c +2 R_LMPD_1_10fthf_c -2 F_LMPD_2_10fthf_c +2 R_LMPD_2_10fthf_c -2 F_LMPD_3_10fthf_c +2 R_LMPD_3_10fthf_c -2 F_LMPD_4_10fthf_c +2 R_LMPD_4_10fthf_c -2 F_LMPD_5_2dmmql8_c +2 R_LMPD_5_2dmmql8_c -2 F_LMPD_22_5mthf_c +2 R_LMPD_22_5mthf_c -2 F_LMPD_23_5mthf_c +2 R_LMPD_23_5mthf_c -2 F_LMPD_51_chor_c +2 R_LMPD_51_chor_c -3 F_LMPD_79_colipa_e +3 R_LMPD_79_colipa_e -3 F_LMPD_80_colipa_e +3 R_LMPD_80_colipa_e -3 F_LMPD_81_colipa_e +3 R_LMPD_81_colipa_e -3 F_LMPD_82_colipa_e +3 R_LMPD_82_colipa_e -3 F_LMPD_83_colipa_e +3 R_LMPD_83_colipa_e -3 F_LMPD_84_colipa_e +3 R_LMPD_84_colipa_e -3 F_LMPD_85_colipa_e +3 R_LMPD_85_colipa_e -3 F_LMPD_86_colipa_e +3 R_LMPD_86_colipa_e -3 F_LMPD_87_colipa_e +3 R_LMPD_87_colipa_e -3 F_LMPD_88_colipa_e +3 R_LMPD_88_colipa_e -3 F_LMPD_89_colipa_e +3 R_LMPD_89_colipa_e -3 F_LMPD_90_colipa_e +3 R_LMPD_90_colipa_e -6 F_LMPD_107_enter_c +6 R_LMPD_107_enter_c -2 F_LMPD_139_mlthf_c +2 R_LMPD_139_mlthf_c -2 F_LMPD_140_mlthf_c +2 R_LMPD_140_mlthf_c -2 F_LMPD_153_mql8_c +2 R_LMPD_153_mql8_c -2 F_LMPD_154_murein3p3p_p +2 R_LMPD_154_murein3p3p_p -2 F_LMPD_155_murein3px4p_p +2 R_LMPD_155_murein3px4p_p -2 F_LMPD_156_murein4p4p_p +2 R_LMPD_156_murein4p4p_p -2 F_LMPD_157_murein4px4p_p +2 R_LMPD_157_murein4px4p_p -3 F_LMPD_158_murein4px4px4p_p +3 R_LMPD_158_murein4px4px4p_p -2 F_LMPD_217_phe-L_c +2 R_LMPD_217_phe-L_c -2 F_LMPD_224_q8h2_c +2 R_LMPD_224_q8h2_c -2 F_LMPD_225_q8h2_c +2 R_LMPD_225_q8h2_c -2 F_LMPD_226_q8h2_c +2 R_LMPD_226_q8h2_c -2 F_LMPD_227_q8h2_c +2 R_LMPD_227_q8h2_c -2 F_LMPD_228_q8h2_c +2 R_LMPD_228_q8h2_c -2 F_LMPD_229_q8h2_c +2 R_LMPD_229_q8h2_c -2 F_LMPD_230_q8h2_c +2 R_LMPD_230_q8h2_c -2 F_LMPD_231_q8h2_c +2 R_LMPD_231_q8h2_c -4 F_LMPD_232_q8h2_c +4 R_LMPD_232_q8h2_c -4 F_LMPD_233_q8h2_c +4 R_LMPD_233_q8h2_c -4 F_LMPD_234_q8h2_c +4 R_LMPD_234_q8h2_c -2 F_LMPD_243_thf_c +2 R_LMPD_243_thf_c -2 F_LMPD_244_thf_c +2 R_LMPD_244_thf_c -2 F_LMPD_245_thf_c +2 R_LMPD_245_thf_c -2 F_LMPD_246_thf_c +2 R_LMPD_246_thf_c -2 F_LMPD_247_thmpp_c +2 R_LMPD_247_thmpp_c -2 F_LMPD_248_thmpp_c +2 R_LMPD_248_thmpp_c -2 F_LMPD_250_trp-L_c +2 R_LMPD_250_trp-L_c -2 F_LMPD_251_tyr-L_c +2 R_LMPD_251_tyr-L_c = 0
M_pg160_c: -0.005707 F_Ec_biomass_iJO1366_WT_53p95M +0.005707 R_Ec_biomass_iJO1366_WT_53p95M - F_PG160abcpp + R_PG160abcpp +2 F_LMPD_52_clpn160_p -2 R_LMPD_52_clpn160_p +2 F_LMPD_53_clpn160_p -2 R_LMPD_53_clpn160_p +2 F_LMPD_54_clpn160_p -2 R_LMPD_54_clpn160_p +2 F_LMPD_55_clpn160_p -2 R_LMPD_55_clpn160_p +2 F_LMPD_56_clpn160_p -2 R_LMPD_56_clpn160_p +2 F_LMPD_57_clpn160_p -2 R_LMPD_57_clpn160_p +2 F_LMPD_58_clpn160_p -2 R_LMPD_58_clpn160_p +2 F_LMPD_59_clpn160_p -2 R_LMPD_59_clpn160_p +2 F_LMPD_60_clpn160_p -2 R_LMPD_60_clpn160_p +2 F_LMPD_61_clpn160_p -2 R_LMPD_61_clpn160_p + F_LMPD_190_pg160_c - R_LMPD_190_pg160_c + F_LMPD_191_pg160_c - R_LMPD_191_pg160_c + F_LMPD_192_pg160_c - R_LMPD_192_pg160_c + F_LMPD_193_pg160_c - R_LMPD_193_pg160_c + F_LMPD_194_pg160_c - R_LMPD_194_pg160_c + F_LMPD_195_pg160_c - R_LMPD_195_pg160_c + F_LMPD_196_pg160_c - R_LMPD_196_pg160_c + F_LMPD_197_pg160_c - R_LMPD_197_pg160_c + F_LMPD_198_pg160_c - R_LMPD_198_pg160_c + F_LMPD_199_pg160_c - R_LMPD_199_pg160_c = 0
M_pg161_c: -0.004439 F_Ec_biomass_iJO1366_WT_53p95M +0.004439 R_Ec_biomass_iJO1366_WT_53p95M - F_PG161abcpp + R_PG161abcpp +2 F_LMPD_62_clpn161_p -2 R_LMPD_62_clpn161_p +2 F_LMPD_63_clpn161_p -2 R_LMPD_63_clpn161_p +2 F_LMPD_64_clpn161_p -2 R_LMPD_64_clpn161_p +2 F_LMPD_65_clpn161_p -2 R_LMPD_65_clpn161_p +2 F_LMPD_66_clpn161_p -2 R_LMPD_66_clpn161_p +2 F_LMPD_67_clpn161_p -2 R_LMPD_67_clpn161_p +2 F_LMPD_68_clpn161_p -2 R_LMPD_68_clpn161_p + F_LMPD_200_pg161_c - R_LMPD_200_pg161_c + F_LMPD_201_pg161_c - R_LMPD_201_pg161_c + F_LMPD_202_pg161_c - R_LMPD_202_pg161_c + F_LMPD_203_pg161_c - R_LMPD_203_pg161_c + F_LMPD_204_pg161_c - R_LMPD_204_pg161_c + F_LMPD_205_pg161_c - R_LMPD_205_pg161_c + F_LMPD_206_pg161_c - R_LMPD_206_pg161_c = 0
M_pg181_c: -0.002288 F_Ec_biomass_iJO1366_WT_53p95M +0.002288 R_Ec_biomass_iJO1366_WT_53p95M - F_PG181abcpp + R_PG181abcpp +2 F_LMPD_69_clpn181_p -2 R_LMPD_69_clpn181_p +2 F_LMPD_70_clpn181_p -2 R_LMPD_70_clpn181_p +2 F_LMPD_71_clpn181_p -2 R_LMPD_71_clpn181_p +2 F_LMPD_72_clpn181_p -2 R_LMPD_72_clpn181_p +2 F_LMPD_73_clpn181_p -2 R_LMPD_73_clpn181_p +2 F_LMPD_74_clpn181_p -2 R_LMPD_74_clpn181_p +2 F_LMPD_75_clpn181_p -2 R_LMPD_75_clpn181_p +2 F_LMPD_76_clpn181_p -2 R_LMPD_76_clpn181_p +2 F_LMPD_77_clpn181_p -2 R_LMPD_77_clpn181_p +2 F_LMPD_78_clpn181_p -2 R_LMPD_78_clpn181_p + F_LMPD_207_pg181_c - R_LMPD_207_pg181_c + F_LMPD_208_pg181_c - R_LMPD_208_pg181_c + F_LMPD_209_pg181_c - R_LMPD_209_pg181_c + F_LMPD_210_pg181_c - R_LMPD_210_pg181_c + F_LMPD_211_pg181_c - R_LMPD_211_pg181_c + F_LMPD_212_pg181_c - R_LMPD_212_pg181_c + F_LMPD_213_pg181_c - R_LMPD_213_pg181_c + F_LMPD_214_pg181_c - R_LMPD_214_pg181_c + F_LMPD_215_pg181_c - R_LMPD_215_pg181_c + F_LMPD_216_pg181_c - R_LMPD_216_pg181_c = 0
M_phe-L_c: -0.18002 F_Ec_biomass_iJO1366_WT_53p95M +0.18002 R_Ec_biomass_iJO1366_WT_53p95M + F_PHEt2rpp - R_PHEt2rpp + F_LMPD_217_phe-L_c - R_LMPD_217_phe-L_c = 0
M_pheme_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M - F_PHEMEabcpp + R_PHEMEabcpp + F_LMPD_218_pheme_c - R_LMPD_218_pheme_c + F_LMPD_219_pheme_c - R_LMPD_219_pheme_c = 0
M_pi_c: +53.9459 F_Ec_biomass_iJO1366_WT_53p95M -53.9459 R_Ec_biomass_iJO1366_WT_53p95M + F_ATPM - R_ATPM - F_ATPS4rpp + R_ATPS4rpp + F_COLIPAabcpp - R_COLIPAabcpp + F_COLIPAabctex - R_COLIPAabctex + F_CYSabc2pp - R_CYSabc2pp + F_FBP - R_FBP + F_FE2abcpp - R_FE2abcpp + F_FE3abcpp - R_FE3abcpp + F_G3PT - R_G3PT + F_G6PP - R_G6PP - F_GAPD + R_GAPD - F_GLCP + R_GLCP + F_GLCabcpp - R_GLCabcpp + F_GLUabcpp - R_GLUabcpp + F_GTHRDabc2pp - R_GTHRDabc2pp + F_Kabcpp - R_Kabcpp + F_MG2uabcpp - R_MG2uabcpp + F_MOBDabcpp - R_MOBDabcpp + F_NADPPPS - R_NADPPPS + F_NI2uabcpp - R_NI2uabcpp + F_NTP1 - R_NTP1 + F_PE160abcpp - R_PE160abcpp + F_PE161abcpp - R_PE161abcpp + F_PE181abcpp - R_PE181abcpp + F_PG160abcpp - R_PG160abcpp + F_PG161abcpp - R_PG161abcpp + F_PG181abcpp - R_PG181abcpp + F_PHEMEabcpp - R_PHEMEabcpp + F_PIt2rpp - R_PIt2rpp +2 F_PIuabcpp -2 R_PIuabcpp +2 F_PPA -2 R_PPA + F_PPA2 - R_PPA2 + F_PPC - R_PPC - F_PPKr + R_PPKr + F_PPS - R_PPS - F_PTAr + R_PTAr + F_SUCOAS - R_SUCOAS + F_SULabcpp - R_SULabcpp + F_ZNabcpp - R_ZNabcpp +18 F_LMPD_1_10fthf_c -18 R_LMPD_1_10fthf_c +18 F_LMPD_2_10fthf_c -18 R_LMPD_2_10fthf_c +18 F_LMPD_3_10fthf_c -18 R_LMPD_3_10fthf_c +18 F_LMPD_4_10fthf_c -18 R_LMPD_4_10fthf_c +4 F_LMPD_5_2dmmql8_c -4 R_LMPD_5_2dmmql8_c +6 F_LMPD_6_2fe2s_c -6 R_LMPD_6_2fe2s_c +6 F_LMPD_7_2fe2s_c -6 R_LMPD_7_2fe2s_c +7 F_LMPD_8_2fe2s_c -7 R_LMPD_8_2fe2s_c +7 F_LMPD_9_2fe2s_c -7 R_LMPD_9_2fe2s_c +6 F_LMPD_10_2fe2s_c -6 R_LMPD_10_2fe2s_c +6 F_LMPD_11_2fe2s_c -6 R_LMPD_11_2fe2s_c +7 F_LMPD_12_2fe2s_c -7 R_LMPD_12_2fe2s_c +7 F_LMPD_13_2fe2s_c -7 R_LMPD_13_2fe2s_c +14 F_LMPD_14_4fe4s_c -14 R_LMPD_14_4fe4s_c +14 F_LMPD_15_4fe4s_c -14 R_LMPD_15_4fe4s_c +12 F_LMPD_16_4fe4s_c -12 R_LMPD_16_4fe4s_c +12 F_LMPD_17_4fe4s_c -12 R_LMPD_17_4fe4s_c +14 F_LMPD_18_4fe4s_c -14 R_LMPD_18_4fe4s_c +14 F_LMPD_19_4fe4s_c -14 R_LMPD_19_4fe4s_c +12 F_LMPD_20_4fe4s_c -12 R_LMPD_20_4fe4s_c +12 F_LMPD_21_4fe4s_c -12 R_LMPD_21_4fe4s_c +18 F_LMPD_22_5mthf_c -18 R_LMPD_22_5mthf_c +18 F_LMPD_23_5mthf_c -18 R_LMPD_23_5mthf_c +14 F_LMPD_24_accoa_c -14 R_LMPD_24_accoa_c +14 F_LMPD_25_accoa_c -14 R_LMPD_25_accoa_c +11 F_LMPD_26_adocbl_c -11 R_LMPD_26_adocbl_c +11 F_LMPD_27_adocbl_c -11 R_LMPD_27_adocbl_c +11 F_LMPD_28_adocbl_c -11 R_LMPD_28_adocbl_c +11 F_LMPD_29_adocbl_c -11 R_LMPD_29_adocbl_c +16 F_LMPD_31_amet_c -16 R_LMPD_31_amet_c +16 F_LMPD_32_amet_c -16 R_LMPD_32_amet_c +2 F_LMPD_33_arg-L_c -2 R_LMPD_33_arg-L_c +52 F_LMPD_35_bmocogdp_c -52 R_LMPD_35_bmocogdp_c +52 F_LMPD_36_bmocogdp_c -52 R_LMPD_36_bmocogdp_c +52 F_LMPD_37_bmocogdp_c -52 R_LMPD_37_bmocogdp_c +52 F_LMPD_38_bmocogdp_c -52 R_LMPD_38_bmocogdp_c +26 F_LMPD_39_btn_c -26 R_LMPD_39_btn_c +27 F_LMPD_40_btn_c -27 R_LMPD_40_btn_c +26 F_LMPD_41_btn_c -26 R_LMPD_41_btn_c +27 F_LMPD_42_btn_c -27 R_LMPD_42_btn_c +51 F_LMPD_43_btn_c -51 R_LMPD_43_btn_c +52 F_LMPD_44_btn_c -52 R_LMPD_44_btn_c +51 F_LMPD_45_btn_c -51 R_LMPD_45_btn_c +52 F_LMPD_47_btn_c -52 R_LMPD_47_btn_c +4 F_LMPD_48_ca2_c -4 R_LMPD_48_ca2_c +4 F_LMPD_49_ca2_c -4 R_LMPD_49_ca2_c +3 F_LMPD_50_ca2_c -3 R_LMPD_50_ca2_c +4 F_LMPD_51_chor_c -4 R_LMPD_51_chor_c +30 F_LMPD_52_clpn160_p -30 R_LMPD_52_clpn160_p +30 F_LMPD_53_clpn160_p -30 R_LMPD_53_clpn160_p +30 F_LMPD_54_clpn160_p -30 R_LMPD_54_clpn160_p +30 F_LMPD_55_clpn160_p -30 R_LMPD_55_clpn160_p +30 F_LMPD_56_clpn160_p -30 R_LMPD_56_clpn160_p +30 F_LMPD_57_clpn160_p -30 R_LMPD_57_clpn160_p +30 F_LMPD_58_clpn160_p -30 R_LMPD_58_clpn160_p +30 F_LMPD_59_clpn160_p -30 R_LMPD_59_clpn160_p +2 F_LMPD_60_clpn160_p -2 R_LMPD_60_clpn160_p +2 F_LMPD_61_clpn160_p -2 R_LMPD_61_clpn160_p +30 F_LMPD_62_clpn161_p -30 R_LMPD_62_clpn161_p +30 F_LMPD_63_clpn161_p -30 R_LMPD_63_clpn161_p +30 F_LMPD_64_clpn161_p -30 R_LMPD_64_clpn161_p +30 F_LMPD_65_clpn161_p -30 R_LMPD_65_clpn161_p +30 F_LMPD_66_clpn161_p -30 R_LMPD_66_clpn161_p +30 F_LMPD_67_clpn161_p -30 R_LMPD_67_clpn161_p +30 F_LMPD_68_clpn161_p -30 R_LMPD_68_clpn161_p +34 F_LMPD_69_clpn181_p -34 R_LMPD_69_clpn181_p +34 F_LMPD_70_clpn181_p -34 R_LMPD_70_clpn181_p +34 F_LMPD_71_clpn181_p -34 R_LMPD_71_clpn181_p +34 F_LMPD_72_clpn181_p -34 R_LMPD_72_clpn181_p +34 F_LMPD_73_clpn181_p -34 R_LMPD_73_clpn181_p +34 F_LMPD_74_clpn181_p -34 R_LMPD_74_clpn181_p +34 F_LMPD_75_clpn181_p -34 R_LMPD_75_clpn181_p +34 F_LMPD_76_clpn181_p -34 R_LMPD_76_clpn181_p +2 F_LMPD_77_clpn181_p -2 R_LMPD_77_clpn181_p +2 F_LMPD_78_clpn181_p -2 R_LMPD_78_clpn181_p +47 F_LMPD_79_colipa_e -47 R_LMPD_79_colipa_e +47 F_LMPD_80_colipa_e -47 R_LMPD_80_colipa_e +47 F_LMPD_81_colipa_e -47 R_LMPD_81_colipa_e +47 F_LMPD_82_colipa_e -47 R_LMPD_82_colipa_e +47 F_LMPD_83_colipa_e -47 R_LMPD_83_colipa_e +47 F_LMPD_84_colipa_e -47 R_LMPD_84_colipa_e +47 F_LMPD_85_colipa_e -47 R_LMPD_85_colipa_e +47 F_LMPD_86_colipa_e -47 R_LMPD_86_colipa_e +47 F_LMPD_87_colipa_e -47 R_LMPD_87_colipa_e +47 F_LMPD_88_colipa_e -47 R_LMPD_88_colipa_e +47 F_LMPD_89_colipa_e -47 R_LMPD_89_colipa_e +47 F_LMPD_90_colipa_e -47 R_LMPD_90_colipa_e +3 F_LMPD_91_ctp_c -3 R_LMPD_91_ctp_c +3 F_LMPD_92_cys-L_c -3 R_LMPD_92_cys-L_c +3 F_LMPD_93_cys-L_c -3 R_LMPD_93_cys-L_c +10 F_LMPD_94_datp_c -10 R_LMPD_94_datp_c +10 F_LMPD_95_datp_c -10 R_LMPD_95_datp_c +10 F_LMPD_96_datp_c -10 R_LMPD_96_datp_c +10 F_LMPD_97_datp_c -10 R_LMPD_97_datp_c +3 F_LMPD_98_dctp_c -3 R_LMPD_98_dctp_c +10 F_LMPD_99_dgtp_c -10 R_LMPD_99_dgtp_c +10 F_LMPD_100_dgtp_c -10 R_LMPD_100_dgtp_c +10 F_LMPD_101_dgtp_c -10 R_LMPD_101_dgtp_c +10 F_LMPD_102_dgtp_c -10 R_LMPD_102_dgtp_c +10 F_LMPD_103_dgtp_c -10 R_LMPD_103_dgtp_c +10 F_LMPD_104_dgtp_c -10 R_LMPD_104_dgtp_c +2 F_LMPD_105_dttp_c -2 R_LMPD_105_dttp_c +2 F_LMPD_106_dttp_c -2 R_LMPD_106_dttp_c +15 F_LMPD_107_enter_c -15 R_LMPD_107_enter_c +23 F_LMPD_108_fad_c -23 R_LMPD_108_fad_c +23 F_LMPD_109_fad_c -23 R_LMPD_109_fad_c +23 F_LMPD_110_fad_c -23 R_LMPD_110_fad_c +23 F_LMPD_111_fad_c -23 R_LMPD_111_fad_c +10 F_LMPD_112_GAM_c -10 R_LMPD_112_GAM_c +10 F_LMPD_113_GAM_c -10 R_LMPD_113_GAM_c +10 F_LMPD_114_GAM_c -10 R_LMPD_114_GAM_c +10 F_LMPD_115_GAM_c -10 R_LMPD_115_GAM_c + F_LMPD_116_gln-L_c - R_LMPD_116_gln-L_c +2 F_LMPD_117_gly_c -2 R_LMPD_117_gly_c +6 F_LMPD_118_gthrd_c -6 R_LMPD_118_gthrd_c +6 F_LMPD_119_gthrd_c -6 R_LMPD_119_gthrd_c +10 F_LMPD_120_gtp_c -10 R_LMPD_120_gtp_c +10 F_LMPD_121_gtp_c -10 R_LMPD_121_gtp_c +10 F_LMPD_122_gtp_c -10 R_LMPD_122_gtp_c +10 F_LMPD_123_gtp_c -10 R_LMPD_123_gtp_c +2 F_LMPD_125_hemeO_c -2 R_LMPD_125_hemeO_c +4 F_LMPD_126_his-L_c -4 R_LMPD_126_his-L_c +4 F_LMPD_127_his-L_c -4 R_LMPD_127_his-L_c +2 F_LMPD_128_ile-L_c -2 R_LMPD_128_ile-L_c +14 F_LMPD_130_lipopb_c -14 R_LMPD_130_lipopb_c +14 F_LMPD_131_lipopb_c -14 R_LMPD_131_lipopb_c +16 F_LMPD_132_lipopb_c -16 R_LMPD_132_lipopb_c +14 F_LMPD_133_lipopb_c -14 R_LMPD_133_lipopb_c + F_LMPD_134_lys-L_c - R_LMPD_134_lys-L_c +15 F_LMPD_135_malcoa_c -15 R_LMPD_135_malcoa_c +15 F_LMPD_136_malcoa_c -15 R_LMPD_136_malcoa_c +5 F_LMPD_137_met-L_c -5 R_LMPD_137_met-L_c +5 F_LMPD_138_met-L_c -5 R_LMPD_138_met-L_c +18 F_LMPD_139_mlthf_c -18 R_LMPD_139_mlthf_c +18 F_LMPD_140_mlthf_c -18 R_LMPD_140_mlthf_c +19 F_LMPD_141_mococdp_c -19 R_LMPD_141_mococdp_c +21 F_LMPD_142_mococdp_c -21 R_LMPD_142_mococdp_c +19 F_LMPD_143_mococdp_c -19 R_LMPD_143_mococdp_c +21 F_LMPD_144_mococdp_c -21 R_LMPD_144_mococdp_c +19 F_LMPD_145_mococdp_c -19 R_LMPD_145_mococdp_c +21 F_LMPD_146_mococdp_c -21 R_LMPD_146_mococdp_c +19 F_LMPD_147_mococdp_c -19 R_LMPD_147_mococdp_c +21 F_LMPD_148_mococdp_c -21 R_LMPD_148_mococdp_c +26 F_LMPD_149_mocogdp_c -26 R_LMPD_149_mocogdp_c +26 F_LMPD_150_mocogdp_c -26 R_LMPD_150_mocogdp_c +26 F_LMPD_151_mocogdp_c -26 R_LMPD_151_mocogdp_c +26 F_LMPD_152_mocogdp_c -26 R_LMPD_152_mocogdp_c +6 F_LMPD_153_mql8_c -6 R_LMPD_153_mql8_c +20 F_LMPD_154_murein3p3p_p -20 R_LMPD_154_murein3p3p_p +20 F_LMPD_155_murein3px4p_p -20 R_LMPD_155_murein3px4p_p +20 F_LMPD_156_murein4p4p_p -20 R_LMPD_156_murein4p4p_p +20 F_LMPD_157_murein4px4p_p -20 R_LMPD_157_murein4px4p_p +30 F_LMPD_158_murein4px4px4p_p -30 R_LMPD_158_murein4px4px4p_p +11 F_LMPD_159_nad_c -11 R_LMPD_159_nad_c +11 F_LMPD_160_nad_c -11 R_LMPD_160_nad_c +11 F_LMPD_161_nad_c -11 R_LMPD_161_nad_c +11 F_LMPD_162_nad_c -11 R_LMPD_162_nad_c +15 F_LMPD_163_pe160_c -15 R_LMPD_163_pe160_c +15 F_LMPD_164_pe160_c -15 R_LMPD_164_pe160_c +15 F_LMPD_165_pe160_c -15 R_LMPD_165_pe160_c +15 F_LMPD_166_pe160_c -15 R_LMPD_166_pe160_c +15 F_LMPD_167_pe160_c -15 R_LMPD_167_pe160_c +15 F_LMPD_168_pe160_c -15 R_LMPD_168_pe160_c +15 F_LMPD_169_pe160_c -15 R_LMPD_169_pe160_c +15 F_LMPD_170_pe160_c -15 R_LMPD_170_pe160_c + F_LMPD_171_pe160_c - R_LMPD_171_pe160_c + F_LMPD_172_pe160_c - R_LMPD_172_pe160_c +15 F_LMPD_173_pe161_c -15 R_LMPD_173_pe161_c +15 F_LMPD_174_pe161_c -15 R_LMPD_174_pe161_c +15 F_LMPD_175_pe161_c -15 R_LMPD_175_pe161_c +15 F_LMPD_176_pe161_c -15 R_LMPD_176_pe161_c +15 F_LMPD_177_pe161_c -15 R_LMPD_177_pe161_c +15 F_LMPD_178_pe161_c -15 R_LMPD_178_pe161_c +15 F_LMPD_179_pe161_c -15 R_LMPD_179_pe161_c +17 F_LMPD_180_pe181_c -17 R_LMPD_180_pe181_c +17 F_LMPD_181_pe181_c -17 R_LMPD_181_pe181_c +17 F_LMPD_182_pe181_c -17 R_LMPD_182_pe181_c +17 F_LMPD_183_pe181_c -17 R_LMPD_183_pe181_c +17 F_LMPD_184_pe181_c -17 R_LMPD_184_pe181_c +17 F_LMPD_185_pe181_c -17 R_LMPD_185_pe181_c +17 F_LMPD_186_pe181_c -17 R_LMPD_186_pe181_c +17 F_LMPD_187_pe181_c -17 R_LMPD_187_pe181_c + F_LMPD_188_pe181_c - R_LMPD_188_pe181_c + F_LMPD_189_pe181_c - R_LMPD_189_pe181_c +15 F_LMPD_190_pg160_c -15 R_LMPD_190_pg160_c +15 F_LMPD_191_pg160_c -15 R_LMPD_191_pg160_c +15 F_LMPD_192_pg160_c -15 R_LMPD_192_pg160_c +15 F_LMPD_193_pg160_c -15 R_LMPD_193_pg160_c +15 F_LMPD_194_pg160_c -15 R_LMPD_194_pg160_c +15 F_LMPD_195_pg160_c -15 R_LMPD_195_pg160_c +15 F_LMPD_196_pg160_c -15 R_LMPD_196_pg160_c +15 F_LMPD_197_pg160_c -15 R_LMPD_197_pg160_c + F_LMPD_198_pg160_c - R_LMPD_198_pg160_c + F_LMPD_199_pg160_c - R_LMPD_199_pg160_c +15 F_LMPD_200_pg161_c -15 R_LMPD_200_pg161_c +15 F_LMPD_201_pg161_c -15 R_LMPD_201_pg161_c +15 F_LMPD_202_pg161_c -15 R_LMPD_202_pg161_c +15 F_LMPD_203_pg161_c -15 R_LMPD_203_pg161_c +15 F_LMPD_204_pg161_c -15 R_LMPD_204_pg161_c +15 F_LMPD_205_pg161_c -15 R_LMPD_205_pg161_c +15 F_LMPD_206_pg161_c -15 R_LMPD_206_pg161_c +17 F_LMPD_207_pg181_c -17 R_LMPD_207_pg181_c +17 F_LMPD_208_pg181_c -17 R_LMPD_208_pg181_c +17 F_LMPD_209_pg181_c -17 R_LMPD_209_pg181_c +17 F_LMPD_210_pg181_c -17 R_LMPD_210_pg181_c +17 F_LMPD_211_pg181_c -17 R_LMPD_211_pg181_c +17 F_LMPD_212_pg181_c -17 R_LMPD_212_pg181_c +17 F_LMPD_213_pg181_c -17 R_LMPD_213_pg181_c +17 F_LMPD_214_pg181_c -17 R_LMPD_214_pg181_c + F_LMPD_215_pg181_c - R_LMPD_215_pg181_c + F_LMPD_216_pg181_c - R_LMPD_216_pg181_c +4 F_LMPD_217_phe-L_c -4 R_LMPD_217_phe-L_c +2 F_LMPD_219_pheme_c -2 R_LMPD_219_pheme_c + F_LMPD_220_pro-L_c - R_LMPD_220_pro-L_c + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c + F_LMPD_222_pydx5p_c - R_LMPD_222_pydx5p_c + F_LMPD_223_pydx5p_c - R_LMPD_223_pydx5p_c +10 F_LMPD_224_q8h2_c -10 R_LMPD_224_q8h2_c +12 F_LMPD_225_q8h2_c -12 R_LMPD_225_q8h2_c +16 F_LMPD_226_q8h2_c -16 R_LMPD_226_q8h2_c +10 F_LMPD_227_q8h2_c -10 R_LMPD_227_q8h2_c +15 F_LMPD_228_q8h2_c -15 R_LMPD_228_q8h2_c +12 F_LMPD_229_q8h2_c -12 R_LMPD_229_q8h2_c +19 F_LMPD_230_q8h2_c -19 R_LMPD_230_q8h2_c +16 F_LMPD_231_q8h2_c -16 R_LMPD_231_q8h2_c +28 F_LMPD_232_q8h2_c -28 R_LMPD_232_q8h2_c +34 F_LMPD_233_q8h2_c -34 R_LMPD_233_q8h2_c +28 F_LMPD_234_q8h2_c -28 R_LMPD_234_q8h2_c +13 F_LMPD_235_ribflv_c -13 R_LMPD_235_ribflv_c +13 F_LMPD_236_ribflv_c -13 R_LMPD_236_ribflv_c +13 F_LMPD_237_ribflv_c -13 R_LMPD_237_ribflv_c +13 F_LMPD_238_ribflv_c -13 R_LMPD_238_ribflv_c + F_LMPD_239_ser-L_c - R_LMPD_239_ser-L_c +4 F_LMPD_240_sheme_c -4 R_LMPD_240_sheme_c +7 F_LMPD_241_spmd_c -7 R_LMPD_241_spmd_c +7 F_LMPD_242_spmd_c -7 R_LMPD_242_spmd_c +17 F_LMPD_243_thf_c -17 R_LMPD_243_thf_c +17 F_LMPD_244_thf_c -17 R_LMPD_244_thf_c +17 F_LMPD_245_thf_c -17 R_LMPD_245_thf_c +17 F_LMPD_246_thf_c -17 R_LMPD_246_thf_c +14 F_LMPD_247_thmpp_c -14 R_LMPD_247_thmpp_c +14 F_LMPD_248_thmpp_c -14 R_LMPD_248_thmpp_c +2 F_LMPD_249_thr-L_c -2 R_LMPD_249_thr-L_c +5 F_LMPD_250_trp-L_c -5 R_LMPD_250_trp-L_c +4 F_LMPD_251_tyr-L_c -4 R_LMPD_251_tyr-L_c + F_LMPD_253_utp_c - R_LMPD_253_utp_c = 0
M_ppi_c: +0.74983 F_Ec_biomass_iJO1366_WT_53p95M -0.74983 R_Ec_biomass_iJO1366_WT_53p95M + F_ACS - R_ACS + F_GLGC - R_GLGC + F_NTPP6 - R_NTPP6 - F_PPA + R_PPA + F_PPA2 - R_PPA2 - F_PPK2r + R_PPK2r + F_PPKr - R_PPKr +4 F_LMPD_1_10fthf_c -4 R_LMPD_1_10fthf_c +4 F_LMPD_2_10fthf_c -4 R_LMPD_2_10fthf_c +4 F_LMPD_3_10fthf_c -4 R_LMPD_3_10fthf_c +4 F_LMPD_4_10fthf_c -4 R_LMPD_4_10fthf_c +17 F_LMPD_5_2dmmql8_c -17 R_LMPD_5_2dmmql8_c +2 F_LMPD_6_2fe2s_c -2 R_LMPD_6_2fe2s_c +2 F_LMPD_7_2fe2s_c -2 R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c +2 F_LMPD_9_2fe2s_c -2 R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c +2 F_LMPD_11_2fe2s_c -2 R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +4 F_LMPD_14_4fe4s_c -4 R_LMPD_14_4fe4s_c +4 F_LMPD_15_4fe4s_c -4 R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c +4 F_LMPD_17_4fe4s_c -4 R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c +4 F_LMPD_19_4fe4s_c -4 R_LMPD_19_4fe4s_c +4 F_LMPD_20_4fe4s_c -4 R_LMPD_20_4fe4s_c +4 F_LMPD_21_4fe4s_c -4 R_LMPD_21_4fe4s_c +4 F_LMPD_22_5mthf_c -4 R_LMPD_22_5mthf_c +4 F_LMPD_23_5mthf_c -4 R_LMPD_23_5mthf_c +5 F_LMPD_24_accoa_c -5 R_LMPD_24_accoa_c +5 F_LMPD_25_accoa_c -5 R_LMPD_25_accoa_c + F_LMPD_26_adocbl_c - R_LMPD_26_adocbl_c + F_LMPD_27_adocbl_c - R_LMPD_27_adocbl_c + F_LMPD_28_adocbl_c - R_LMPD_28_adocbl_c + F_LMPD_29_adocbl_c - R_LMPD_29_adocbl_c +3 F_LMPD_31_amet_c -3 R_LMPD_31_amet_c +3 F_LMPD_32_amet_c -3 R_LMPD_32_amet_c + F_LMPD_33_arg-L_c - R_LMPD_33_arg-L_c + F_LMPD_34_asn-L_c - R_LMPD_34_asn-L_c +22 F_LMPD_35_bmocogdp_c -22 R_LMPD_35_bmocogdp_c +22 F_LMPD_36_bmocogdp_c -22 R_LMPD_36_bmocogdp_c +22 F_LMPD_37_bmocogdp_c -22 R_LMPD_37_bmocogdp_c +22 F_LMPD_38_bmocogdp_c -22 R_LMPD_38_bmocogdp_c +8 F_LMPD_39_btn_c -8 R_LMPD_39_btn_c +8 F_LMPD_40_btn_c -8 R_LMPD_40_btn_c +8 F_LMPD_41_btn_c -8 R_LMPD_41_btn_c +8 F_LMPD_42_btn_c -8 R_LMPD_42_btn_c +9 F_LMPD_43_btn_c -9 R_LMPD_43_btn_c +9 F_LMPD_44_btn_c -9 R_LMPD_44_btn_c +9 F_LMPD_45_btn_c -9 R_LMPD_45_btn_c +9 F_LMPD_47_btn_c -9 R_LMPD_47_btn_c +2 F_LMPD_52_clpn160_p -2 R_LMPD_52_clpn160_p +2 F_LMPD_53_clpn160_p -2 R_LMPD_53_clpn160_p +2 F_LMPD_54_clpn160_p -2 R_LMPD_54_clpn160_p +2 F_LMPD_55_clpn160_p -2 R_LMPD_55_clpn160_p +2 F_LMPD_56_clpn160_p -2 R_LMPD_56_clpn160_p +2 F_LMPD_57_clpn160_p -2 R_LMPD_57_clpn160_p +2 F_LMPD_58_clpn160_p -2 R_LMPD_58_clpn160_p +2 F_LMPD_59_clpn160_p -2 R_LMPD_59_clpn160_p +6 F_LMPD_60_clpn160_p -6 R_LMPD_60_clpn160_p +6 F_LMPD_61_clpn160_p -6 R_LMPD_61_clpn160_p +2 F_LMPD_62_clpn161_p -2 R_LMPD_62_clpn161_p +2 F_LMPD_63_clpn161_p -2 R_LMPD_63_clpn161_p +2 F_LMPD_64_clpn161_p -2 R_LMPD_64_clpn161_p +2 F_LMPD_65_clpn161_p -2 R_LMPD_65_clpn161_p +2 F_LMPD_66_clpn161_p -2 R_LMPD_66_clpn161_p +2 F_LMPD_67_clpn161_p -2 R_LMPD_67_clpn161_p +2 F_LMPD_68_clpn161_p -2 R_LMPD_68_clpn161_p +2 F_LMPD_69_clpn181_p -2 R_LMPD_69_clpn181_p +2 F_LMPD_70_clpn181_p -2 R_LMPD_70_clpn181_p +2 F_LMPD_71_clpn181_p -2 R_LMPD_71_clpn181_p +2 F_LMPD_72_clpn181_p -2 R_LMPD_72_clpn181_p +2 F_LMPD_73_clpn181_p -2 R_LMPD_73_clpn181_p +2 F_LMPD_74_clpn181_p -2 R_LMPD_74_clpn181_p +2 F_LMPD_75_clpn181_p -2 R_LMPD_75_clpn181_p +2 F_LMPD_76_clpn181_p -2 R_LMPD_76_clpn181_p +6 F_LMPD_77_clpn181_p -6 R_LMPD_77_clpn181_p +6 F_LMPD_78_clpn181_p -6 R_LMPD_78_clpn181_p +14 F_LMPD_79_colipa_e -14 R_LMPD_79_colipa_e +14 F_LMPD_80_colipa_e -14 R_LMPD_80_colipa_e +14 F_LMPD_81_colipa_e -14 R_LMPD_81_colipa_e +14 F_LMPD_82_colipa_e -14 R_LMPD_82_colipa_e +14 F_LMPD_83_colipa_e -14 R_LMPD_83_colipa_e +14 F_LMPD_84_colipa_e -14 R_LMPD_84_colipa_e +14 F_LMPD_85_colipa_e -14 R_LMPD_85_colipa_e +14 F_LMPD_86_colipa_e -14 R_LMPD_86_colipa_e +14 F_LMPD_87_colipa_e -14 R_LMPD_87_colipa_e +14 F_LMPD_88_colipa_e -14 R_LMPD_88_colipa_e +14 F_LMPD_89_colipa_e -14 R_LMPD_89_colipa_e +14 F_LMPD_90_colipa_e -14 R_LMPD_90_colipa_e + F_LMPD_91_ctp_c - R_LMPD_91_ctp_c + F_LMPD_92_cys-L_c - R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c + F_LMPD_94_datp_c - R_LMPD_94_datp_c + F_LMPD_95_datp_c - R_LMPD_95_datp_c + F_LMPD_96_datp_c - R_LMPD_96_datp_c + F_LMPD_97_datp_c - R_LMPD_97_datp_c + F_LMPD_98_dctp_c - R_LMPD_98_dctp_c +2 F_LMPD_99_dgtp_c -2 R_LMPD_99_dgtp_c +2 F_LMPD_100_dgtp_c -2 R_LMPD_100_dgtp_c +2 F_LMPD_101_dgtp_c -2 R_LMPD_101_dgtp_c +2 F_LMPD_102_dgtp_c -2 R_LMPD_102_dgtp_c +2 F_LMPD_103_dgtp_c -2 R_LMPD_103_dgtp_c +2 F_LMPD_104_dgtp_c -2 R_LMPD_104_dgtp_c +2 F_LMPD_105_dttp_c -2 R_LMPD_105_dttp_c +2 F_LMPD_106_dttp_c -2 R_LMPD_106_dttp_c +6 F_LMPD_107_enter_c -6 R_LMPD_107_enter_c +5 F_LMPD_108_fad_c -5 R_LMPD_108_fad_c +5 F_LMPD_109_fad_c -5 R_LMPD_109_fad_c +5 F_LMPD_110_fad_c -5 R_LMPD_110_fad_c +5 F_LMPD_111_fad_c -5 R_LMPD_111_fad_c + F_LMPD_112_GAM_c - R_LMPD_112_GAM_c + F_LMPD_113_GAM_c - R_LMPD_113_GAM_c + F_LMPD_114_GAM_c - R_LMPD_114_GAM_c + F_LMPD_115_GAM_c - R_LMPD_115_GAM_c + F_LMPD_118_gthrd_c - R_LMPD_118_gthrd_c + F_LMPD_119_gthrd_c - R_LMPD_119_gthrd_c +2 F_LMPD_120_gtp_c -2 R_LMPD_120_gtp_c +2 F_LMPD_121_gtp_c -2 R_LMPD_121_gtp_c +2 F_LMPD_122_gtp_c -2 R_LMPD_122_gtp_c +2 F_LMPD_123_gtp_c -2 R_LMPD_123_gtp_c +14 F_LMPD_124_hemeO_c -14 R_LMPD_124_hemeO_c +18 F_LMPD_125_hemeO_c -18 R_LMPD_125_hemeO_c +2 F_LMPD_126_his-L_c -2 R_LMPD_126_his-L_c +2 F_LMPD_127_his-L_c -2 R_LMPD_127_his-L_c +9 F_LMPD_130_lipopb_c -9 R_LMPD_130_lipopb_c +9 F_LMPD_131_lipopb_c -9 R_LMPD_131_lipopb_c +9 F_LMPD_132_lipopb_c -9 R_LMPD_132_lipopb_c +9 F_LMPD_133_lipopb_c -9 R_LMPD_133_lipopb_c +5 F_LMPD_135_malcoa_c -5 R_LMPD_135_malcoa_c +5 F_LMPD_136_malcoa_c -5 R_LMPD_136_malcoa_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c +4 F_LMPD_139_mlthf_c -4 R_LMPD_139_mlthf_c +4 F_LMPD_140_mlthf_c -4 R_LMPD_140_mlthf_c +10 F_LMPD_141_mococdp_c -10 R_LMPD_141_mococdp_c +10 F_LMPD_142_mococdp_c -10 R_LMPD_142_mococdp_c +10 F_LMPD_143_mococdp_c -10 R_LMPD_143_mococdp_c +10 F_LMPD_144_mococdp_c -10 R_LMPD_144_mococdp_c +10 F_LMPD_145_mococdp_c -10 R_LMPD_145_mococdp_c +10 F_LMPD_146_mococdp_c -10 R_LMPD_146_mococdp_c +10 F_LMPD_147_mococdp_c -10 R_LMPD_147_mococdp_c +10 F_LMPD_148_mococdp_c -10 R_LMPD_148_mococdp_c +11 F_LMPD_149_mocogdp_c -11 R_LMPD_149_mocogdp_c +11 F_LMPD_150_mocogdp_c -11 R_LMPD_150_mocogdp_c +11 F_LMPD_151_mocogdp_c -11 R_LMPD_151_mocogdp_c +11 F_LMPD_152_mocogdp_c -11 R_LMPD_152_mocogdp_c +19 F_LMPD_153_mql8_c -19 R_LMPD_153_mql8_c +4 F_LMPD_154_murein3p3p_p -4 R_LMPD_154_murein3p3p_p +4 F_LMPD_155_murein3px4p_p -4 R_LMPD_155_murein3px4p_p +4 F_LMPD_156_murein4p4p_p -4 R_LMPD_156_murein4p4p_p +4 F_LMPD_157_murein4px4p_p -4 R_LMPD_157_murein4px4p_p +6 F_LMPD_158_murein4px4px4p_p -6 R_LMPD_158_murein4px4px4p_p +4 F_LMPD_159_nad_c -4 R_LMPD_159_nad_c +4 F_LMPD_160_nad_c -4 R_LMPD_160_nad_c +4 F_LMPD_161_nad_c -4 R_LMPD_161_nad_c +4 F_LMPD_162_nad_c -4 R_LMPD_162_nad_c + F_LMPD_163_pe160_c - R_LMPD_163_pe160_c + F_LMPD_164_pe160_c - R_LMPD_164_pe160_c + F_LMPD_165_pe160_c - R_LMPD_165_pe160_c + F_LMPD_166_pe160_c - R_LMPD_166_pe160_c + F_LMPD_167_pe160_c - R_LMPD_167_pe160_c + F_LMPD_168_pe160_c - R_LMPD_168_pe160_c + F_LMPD_169_pe160_c - R_LMPD_169_pe160_c + F_LMPD_170_pe160_c - R_LMPD_170_pe160_c +3 F_LMPD_171_pe160_c -3 R_LMPD_171_pe160_c +3 F_LMPD_172_pe160_c -3 R_LMPD_172_pe160_c + F_LMPD_173_pe161_c - R_LMPD_173_pe161_c + F_LMPD_174_pe161_c - R_LMPD_174_pe161_c + F_LMPD_175_pe161_c - R_LMPD_175_pe161_c + F_LMPD_176_pe161_c - R_LMPD_176_pe161_c + F_LMPD_177_pe161_c - R_LMPD_177_pe161_c + F_LMPD_178_pe161_c - R_LMPD_178_pe161_c + F_LMPD_179_pe161_c - R_LMPD_179_pe161_c + F_LMPD_180_pe181_c - R_LMPD_180_pe181_c + F_LMPD_181_pe181_c - R_LMPD_181_pe181_c + F_LMPD_182_pe181_c - R_LMPD_182_pe181_c + F_LMPD_183_pe181_c - R_LMPD_183_pe181_c + F_LMPD_184_pe181_c - R_LMPD_184_pe181_c + F_LMPD_185_pe181_c - R_LMPD_185_pe181_c + F_LMPD_186_pe181_c - R_LMPD_186_pe181_c + F_LMPD_187_pe181_c - R_LMPD_187_pe181_c +3 F_LMPD_188_pe181_c -3 R_LMPD_188_pe181_c +3 F_LMPD_189_pe181_c -3 R_LMPD_189_pe181_c + F_LMPD_190_pg160_c - R_LMPD_190_pg160_c + F_LMPD_191_pg160_c - R_LMPD_191_pg160_c + F_LMPD_192_pg160_c - R_LMPD_192_pg160_c + F_LMPD_193_pg160_c - R_LMPD_193_pg160_c + F_LMPD_194_pg160_c - R_LMPD_194_pg160_c + F_LMPD_195_pg160_c - R_LMPD_195_pg160_c + F_LMPD_196_pg160_c - R_LMPD_196_pg160_c + F_LMPD_197_pg160_c - R_LMPD_197_pg160_c +3 F_LMPD_198_pg160_c -3 R_LMPD_198_pg160_c +3 F_LMPD_199_pg160_c -3 R_LMPD_199_pg160_c + F_LMPD_200_pg161_c - R_LMPD_200_pg161_c + F_LMPD_201_pg161_c - R_LMPD_201_pg161_c + F_LMPD_202_pg161_c - R_LMPD_202_pg161_c + F_LMPD_203_pg161_c - R_LMPD_203_pg161_c + F_LMPD_204_pg161_c - R_LMPD_204_pg161_c + F_LMPD_205_pg161_c - R_LMPD_205_pg161_c + F_LMPD_206_pg161_c - R_LMPD_206_pg161_c + F_LMPD_207_pg181_c - R_LMPD_207_pg181_c + F_LMPD_208_pg181_c - R_LMPD_208_pg181_c + F_LMPD_209_pg181_c - R_LMPD_209_pg181_c + F_LMPD_210_pg181_c - R_LMPD_210_pg181_c + F_LMPD_211_pg181_c - R_LMPD_211_pg181_c + F_LMPD_212_pg181_c - R_LMPD_212_pg181_c + F_LMPD_213_pg181_c - R_LMPD_213_pg181_c + F_LMPD_214_pg181_c - R_LMPD_214_pg181_c +3 F_LMPD_215_pg181_c -3 R_LMPD_215_pg181_c +3 F_LMPD_216_pg181_c -3 R_LMPD_216_pg181_c +8 F_LMPD_218_pheme_c -8 R_LMPD_218_pheme_c +12 F_LMPD_219_pheme_c -12 R_LMPD_219_pheme_c +22 F_LMPD_224_q8h2_c -22 R_LMPD_224_q8h2_c +22 F_LMPD_225_q8h2_c -22 R_LMPD_225_q8h2_c +22 F_LMPD_226_q8h2_c -22 R_LMPD_226_q8h2_c +22 F_LMPD_227_q8h2_c -22 R_LMPD_227_q8h2_c +19 F_LMPD_228_q8h2_c -19 R_LMPD_228_q8h2_c +22 F_LMPD_229_q8h2_c -22 R_LMPD_229_q8h2_c +19 F_LMPD_230_q8h2_c -19 R_LMPD_230_q8h2_c +22 F_LMPD_231_q8h2_c -22 R_LMPD_231_q8h2_c +44 F_LMPD_232_q8h2_c -44 R_LMPD_232_q8h2_c +38 F_LMPD_233_q8h2_c -38 R_LMPD_233_q8h2_c +44 F_LMPD_234_q8h2_c -44 R_LMPD_234_q8h2_c +3 F_LMPD_235_ribflv_c -3 R_LMPD_235_ribflv_c +3 F_LMPD_236_ribflv_c -3 R_LMPD_236_ribflv_c +3 F_LMPD_237_ribflv_c -3 R_LMPD_237_ribflv_c +3 F_LMPD_238_ribflv_c -3 R_LMPD_238_ribflv_c +12 F_LMPD_240_sheme_c -12 R_LMPD_240_sheme_c +3 F_LMPD_241_spmd_c -3 R_LMPD_241_spmd_c +3 F_LMPD_242_spmd_c -3 R_LMPD_242_spmd_c +4 F_LMPD_243_thf_c -4 R_LMPD_243_thf_c +4 F_LMPD_244_thf_c -4 R_LMPD_244_thf_c +4 F_LMPD_245_thf_c -4 R_LMPD_245_thf_c +4 F_LMPD_246_thf_c -4 R_LMPD_246_thf_c +6 F_LMPD_247_thmpp_c -6 R_LMPD_247_thmpp_c +6 F_LMPD_248_thmpp_c -6 R_LMPD_248_thmpp_c + F_LMPD_250_trp-L_c - R_LMPD_250_trp-L_c +21 F_LMPD_252_udcpdp_c -21 R_LMPD_252_udcpdp_c + F_LMPD_253_utp_c - R_LMPD_253_utp_c = 0
M_ppp9_c: F_PPPGO - R_PPPGO + F_PPPGO3 - R_PPPGO3 - F_LMPD_124_hemeO_c + R_LMPD_124_hemeO_c - F_LMPD_125_hemeO_c + R_LMPD_125_hemeO_c - F_LMPD_218_pheme_c + R_LMPD_218_pheme_c - F_LMPD_219_pheme_c + R_LMPD_219_pheme_c = 0
M_pppg9_c: - F_PPPGO + R_PPPGO - F_PPPGO3 + R_PPPGO3 + F_LMPD_124_hemeO_c - R_LMPD_124_hemeO_c + F_LMPD_125_hemeO_c - R_LMPD_125_hemeO_c + F_LMPD_218_pheme_c - R_LMPD_218_pheme_c + F_LMPD_219_pheme_c - R_LMPD_219_pheme_c = 0
M_pppi_c: - F_PPA2 + R_PPA2 + F_PPK2r - R_PPK2r + F_LMPD_26_adocbl_c - R_LMPD_26_adocbl_c + F_LMPD_27_adocbl_c - R_LMPD_27_adocbl_c + F_LMPD_28_adocbl_c - R_LMPD_28_adocbl_c + F_LMPD_29_adocbl_c - R_LMPD_29_adocbl_c = 0
M_pro-L_c: -0.2148 F_Ec_biomass_iJO1366_WT_53p95M +0.2148 R_Ec_biomass_iJO1366_WT_53p95M + F_PROt2rpp - R_PROt2rpp + F_LMPD_220_pro-L_c - R_LMPD_220_pro-L_c = 0
M_ptrc_c: -0.03327 F_Ec_biomass_iJO1366_WT_53p95M +0.03327 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_221_ptrc_c - R_LMPD_221_ptrc_c = 0
M_pydx5p_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_222_pydx5p_c - R_LMPD_222_pydx5p_c + F_LMPD_223_pydx5p_c - R_LMPD_223_pydx5p_c = 0
M_pyr_c: F_DHAPT - R_DHAPT + F_EDA - R_EDA + F_GLCptspp - R_GLCptspp + F_LDH_D - R_LDH_D + F_LDH_D2 - R_LDH_D2 + F_ME1 - R_ME1 + F_ME2 - R_ME2 + F_OAADC - R_OAADC - F_PDH + R_PDH - F_PFL + R_PFL - F_POR5 + R_POR5 - F_POX + R_POX - F_PPS + R_PPS + F_PYK - R_PYK + F_PYRt2rpp - R_PYRt2rpp + F_LMPD_1_10fthf_c - R_LMPD_1_10fthf_c + F_LMPD_2_10fthf_c - R_LMPD_2_10fthf_c + F_LMPD_3_10fthf_c - R_LMPD_3_10fthf_c + F_LMPD_4_10fthf_c - R_LMPD_4_10fthf_c -7 F_LMPD_5_2dmmql8_c +7 R_LMPD_5_2dmmql8_c + F_LMPD_22_5mthf_c - R_LMPD_22_5mthf_c + F_LMPD_23_5mthf_c - R_LMPD_23_5mthf_c -2 F_LMPD_24_accoa_c +2 R_LMPD_24_accoa_c -2 F_LMPD_25_accoa_c +2 R_LMPD_25_accoa_c - F_LMPD_30_ala-L_c + R_LMPD_30_ala-L_c + F_LMPD_31_amet_c - R_LMPD_31_amet_c + F_LMPD_32_amet_c - R_LMPD_32_amet_c + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c +2 F_LMPD_43_btn_c -2 R_LMPD_43_btn_c +2 F_LMPD_44_btn_c -2 R_LMPD_44_btn_c +2 F_LMPD_45_btn_c -2 R_LMPD_45_btn_c +2 F_LMPD_47_btn_c -2 R_LMPD_47_btn_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c +3 F_LMPD_107_enter_c -3 R_LMPD_107_enter_c -3 F_LMPD_124_hemeO_c +3 R_LMPD_124_hemeO_c -3 F_LMPD_125_hemeO_c +3 R_LMPD_125_hemeO_c - F_LMPD_128_ile-L_c + R_LMPD_128_ile-L_c -2 F_LMPD_129_leu-L_c +2 R_LMPD_129_leu-L_c - F_LMPD_134_lys-L_c + R_LMPD_134_lys-L_c -2 F_LMPD_135_malcoa_c +2 R_LMPD_135_malcoa_c -2 F_LMPD_136_malcoa_c +2 R_LMPD_136_malcoa_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c + F_LMPD_139_mlthf_c - R_LMPD_139_mlthf_c + F_LMPD_140_mlthf_c - R_LMPD_140_mlthf_c -7 F_LMPD_153_mql8_c +7 R_LMPD_153_mql8_c -4 F_LMPD_154_murein3p3p_p +4 R_LMPD_154_murein3p3p_p -8 F_LMPD_155_murein3px4p_p +8 R_LMPD_155_murein3px4p_p -8 F_LMPD_156_murein4p4p_p +8 R_LMPD_156_murein4p4p_p -8 F_LMPD_157_murein4px4p_p +8 R_LMPD_157_murein4px4p_p -12 F_LMPD_158_murein4px4px4p_p +12 R_LMPD_158_murein4px4px4p_p - F_LMPD_222_pydx5p_c + R_LMPD_222_pydx5p_c - F_LMPD_223_pydx5p_c + R_LMPD_223_pydx5p_c -7 F_LMPD_224_q8h2_c +7 R_LMPD_224_q8h2_c -7 F_LMPD_225_q8h2_c +7 R_LMPD_225_q8h2_c -7 F_LMPD_226_q8h2_c +7 R_LMPD_226_q8h2_c -7 F_LMPD_227_q8h2_c +7 R_LMPD_227_q8h2_c -7 F_LMPD_228_q8h2_c +7 R_LMPD_228_q8h2_c -7 F_LMPD_229_q8h2_c +7 R_LMPD_229_q8h2_c -7 F_LMPD_230_q8h2_c +7 R_LMPD_230_q8h2_c -7 F_LMPD_231_q8h2_c +7 R_LMPD_231_q8h2_c -14 F_LMPD_232_q8h2_c +14 R_LMPD_232_q8h2_c -14 F_LMPD_233_q8h2_c +14 R_LMPD_233_q8h2_c -14 F_LMPD_234_q8h2_c +14 R_LMPD_234_q8h2_c + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c + F_LMPD_243_thf_c - R_LMPD_243_thf_c + F_LMPD_244_thf_c - R_LMPD_244_thf_c + F_LMPD_245_thf_c - R_LMPD_245_thf_c + F_LMPD_246_thf_c - R_LMPD_246_thf_c - F_LMPD_247_thmpp_c + R_LMPD_247_thmpp_c - F_LMPD_248_thmpp_c + R_LMPD_248_thmpp_c -11 F_LMPD_252_udcpdp_c +11 R_LMPD_252_udcpdp_c -2 F_LMPD_254_val-L_c +2 R_LMPD_254_val-L_c = 0
M_q8_c: - F_ASPO3 + R_ASPO3 + F_CYTBDpp - R_CYTBDpp + F_CYTBO3_4pp - R_CYTBO3_4pp - F_DHORD2 + R_DHORD2 - F_FDH4pp + R_FDH4pp - F_G3PD5 + R_G3PD5 - F_GLCDpp + R_GLCDpp - F_GLYCTO2 + R_GLYCTO2 - F_LDH_D2 + R_LDH_D2 - F_MDH2 + R_MDH2 - F_NADH16pp + R_NADH16pp - F_NADH5 + R_NADH5 - F_NADPHQR2 + R_NADPHQR2 - F_POX + R_POX - F_SUCDi + R_SUCDi = 0
M_q8h2_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_ASPO3 - R_ASPO3 - F_CYTBDpp + R_CYTBDpp - F_CYTBO3_4pp + R_CYTBO3_4pp + F_DHORD2 - R_DHORD2 + F_FDH4pp - R_FDH4pp + F_G3PD5 - R_G3PD5 + F_GLCDpp - R_GLCDpp + F_GLYCTO2 - R_GLYCTO2 + F_LDH_D2 - R_LDH_D2 + F_MDH2 - R_MDH2 + F_NADH16pp - R_NADH16pp + F_NADH5 - R_NADH5 + F_NADPHQR2 - R_NADPHQR2 + F_POX - R_POX + F_SUCDi - R_SUCDi + F_LMPD_224_q8h2_c - R_LMPD_224_q8h2_c + F_LMPD_225_q8h2_c - R_LMPD_225_q8h2_c + F_LMPD_226_q8h2_c - R_LMPD_226_q8h2_c + F_LMPD_227_q8h2_c - R_LMPD_227_q8h2_c + F_LMPD_228_q8h2_c - R_LMPD_228_q8h2_c + F_LMPD_229_q8h2_c - R_LMPD_229_q8h2_c + F_LMPD_230_q8h2_c - R_LMPD_230_q8h2_c + F_LMPD_231_q8h2_c - R_LMPD_231_q8h2_c +2 F_LMPD_232_q8h2_c -2 R_LMPD_232_q8h2_c +2 F_LMPD_233_q8h2_c -2 R_LMPD_233_q8h2_c +2 F_LMPD_234_q8h2_c -2 R_LMPD_234_q8h2_c = 0
M_r5p_c: - F_RPI + R_RPI - F_TKT1 + R_TKT1 - F_LMPD_1_10fthf_c + R_LMPD_1_10fthf_c - F_LMPD_2_10fthf_c + R_LMPD_2_10fthf_c - F_LMPD_3_10fthf_c + R_LMPD_3_10fthf_c - F_LMPD_4_10fthf_c + R_LMPD_4_10fthf_c - F_LMPD_22_5mthf_c + R_LMPD_22_5mthf_c - F_LMPD_23_5mthf_c + R_LMPD_23_5mthf_c - F_LMPD_24_accoa_c + R_LMPD_24_accoa_c - F_LMPD_25_accoa_c + R_LMPD_25_accoa_c - F_LMPD_26_adocbl_c + R_LMPD_26_adocbl_c - F_LMPD_27_adocbl_c + R_LMPD_27_adocbl_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c -4 F_LMPD_35_bmocogdp_c +4 R_LMPD_35_bmocogdp_c -4 F_LMPD_36_bmocogdp_c +4 R_LMPD_36_bmocogdp_c -4 F_LMPD_37_bmocogdp_c +4 R_LMPD_37_bmocogdp_c -4 F_LMPD_38_bmocogdp_c +4 R_LMPD_38_bmocogdp_c -3 F_LMPD_39_btn_c +3 R_LMPD_39_btn_c -3 F_LMPD_40_btn_c +3 R_LMPD_40_btn_c -3 F_LMPD_41_btn_c +3 R_LMPD_41_btn_c -3 F_LMPD_42_btn_c +3 R_LMPD_42_btn_c -3 F_LMPD_43_btn_c +3 R_LMPD_43_btn_c -3 F_LMPD_44_btn_c +3 R_LMPD_44_btn_c -3 F_LMPD_45_btn_c +3 R_LMPD_45_btn_c -3 F_LMPD_47_btn_c +3 R_LMPD_47_btn_c - F_LMPD_91_ctp_c + R_LMPD_91_ctp_c - F_LMPD_94_datp_c + R_LMPD_94_datp_c - F_LMPD_95_datp_c + R_LMPD_95_datp_c - F_LMPD_96_datp_c + R_LMPD_96_datp_c - F_LMPD_97_datp_c + R_LMPD_97_datp_c - F_LMPD_98_dctp_c + R_LMPD_98_dctp_c - F_LMPD_99_dgtp_c + R_LMPD_99_dgtp_c - F_LMPD_100_dgtp_c + R_LMPD_100_dgtp_c - F_LMPD_101_dgtp_c + R_LMPD_101_dgtp_c - F_LMPD_102_dgtp_c + R_LMPD_102_dgtp_c - F_LMPD_103_dgtp_c + R_LMPD_103_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c - F_LMPD_105_dttp_c + R_LMPD_105_dttp_c - F_LMPD_106_dttp_c + R_LMPD_106_dttp_c -2 F_LMPD_108_fad_c +2 R_LMPD_108_fad_c -2 F_LMPD_109_fad_c +2 R_LMPD_109_fad_c -2 F_LMPD_110_fad_c +2 R_LMPD_110_fad_c -2 F_LMPD_111_fad_c +2 R_LMPD_111_fad_c - F_LMPD_112_GAM_c + R_LMPD_112_GAM_c - F_LMPD_113_GAM_c + R_LMPD_113_GAM_c - F_LMPD_114_GAM_c + R_LMPD_114_GAM_c - F_LMPD_115_GAM_c + R_LMPD_115_GAM_c - F_LMPD_120_gtp_c + R_LMPD_120_gtp_c - F_LMPD_121_gtp_c + R_LMPD_121_gtp_c - F_LMPD_122_gtp_c + R_LMPD_122_gtp_c - F_LMPD_123_gtp_c + R_LMPD_123_gtp_c -2 F_LMPD_125_hemeO_c +2 R_LMPD_125_hemeO_c - F_LMPD_126_his-L_c + R_LMPD_126_his-L_c - F_LMPD_127_his-L_c + R_LMPD_127_his-L_c -2 F_LMPD_130_lipopb_c +2 R_LMPD_130_lipopb_c -2 F_LMPD_131_lipopb_c +2 R_LMPD_131_lipopb_c -2 F_LMPD_132_lipopb_c +2 R_LMPD_132_lipopb_c -2 F_LMPD_133_lipopb_c +2 R_LMPD_133_lipopb_c - F_LMPD_135_malcoa_c + R_LMPD_135_malcoa_c - F_LMPD_136_malcoa_c + R_LMPD_136_malcoa_c - F_LMPD_139_mlthf_c + R_LMPD_139_mlthf_c - F_LMPD_140_mlthf_c + R_LMPD_140_mlthf_c -2 F_LMPD_141_mococdp_c +2 R_LMPD_141_mococdp_c -2 F_LMPD_142_mococdp_c +2 R_LMPD_142_mococdp_c -2 F_LMPD_143_mococdp_c +2 R_LMPD_143_mococdp_c -2 F_LMPD_144_mococdp_c +2 R_LMPD_144_mococdp_c -2 F_LMPD_145_mococdp_c +2 R_LMPD_145_mococdp_c -2 F_LMPD_146_mococdp_c +2 R_LMPD_146_mococdp_c -2 F_LMPD_147_mococdp_c +2 R_LMPD_147_mococdp_c -2 F_LMPD_148_mococdp_c +2 R_LMPD_148_mococdp_c -2 F_LMPD_149_mocogdp_c +2 R_LMPD_149_mocogdp_c -2 F_LMPD_150_mocogdp_c +2 R_LMPD_150_mocogdp_c -2 F_LMPD_151_mocogdp_c +2 R_LMPD_151_mocogdp_c -2 F_LMPD_152_mocogdp_c +2 R_LMPD_152_mocogdp_c - F_LMPD_153_mql8_c + R_LMPD_153_mql8_c -2 F_LMPD_159_nad_c +2 R_LMPD_159_nad_c -2 F_LMPD_160_nad_c +2 R_LMPD_160_nad_c -2 F_LMPD_161_nad_c +2 R_LMPD_161_nad_c -2 F_LMPD_162_nad_c +2 R_LMPD_162_nad_c -2 F_LMPD_219_pheme_c +2 R_LMPD_219_pheme_c -3 F_LMPD_224_q8h2_c +3 R_LMPD_224_q8h2_c -3 F_LMPD_225_q8h2_c +3 R_LMPD_225_q8h2_c -3 F_LMPD_226_q8h2_c +3 R_LMPD_226_q8h2_c -3 F_LMPD_227_q8h2_c +3 R_LMPD_227_q8h2_c -3 F_LMPD_228_q8h2_c +3 R_LMPD_228_q8h2_c -3 F_LMPD_229_q8h2_c +3 R_LMPD_229_q8h2_c -3 F_LMPD_230_q8h2_c +3 R_LMPD_230_q8h2_c -3 F_LMPD_231_q8h2_c +3 R_LMPD_231_q8h2_c -6 F_LMPD_232_q8h2_c +6 R_LMPD_232_q8h2_c -6 F_LMPD_233_q8h2_c +6 R_LMPD_233_q8h2_c -6 F_LMPD_234_q8h2_c +6 R_LMPD_234_q8h2_c - F_LMPD_235_ribflv_c + R_LMPD_235_ribflv_c - F_LMPD_236_ribflv_c + R_LMPD_236_ribflv_c - F_LMPD_237_ribflv_c + R_LMPD_237_ribflv_c - F_LMPD_238_ribflv_c + R_LMPD_238_ribflv_c -2 F_LMPD_240_sheme_c +2 R_LMPD_240_sheme_c - F_LMPD_241_spmd_c + R_LMPD_241_spmd_c - F_LMPD_242_spmd_c + R_LMPD_242_spmd_c - F_LMPD_243_thf_c + R_LMPD_243_thf_c - F_LMPD_244_thf_c + R_LMPD_244_thf_c - F_LMPD_245_thf_c + R_LMPD_245_thf_c - F_LMPD_246_thf_c + R_LMPD_246_thf_c -2 F_LMPD_247_thmpp_c +2 R_LMPD_247_thmpp_c -2 F_LMPD_248_thmpp_c +2 R_LMPD_248_thmpp_c - F_LMPD_250_trp-L_c + R_LMPD_250_trp-L_c - F_LMPD_253_utp_c + R_LMPD_253_utp_c = 0
M_ribflv_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_235_ribflv_c - R_LMPD_235_ribflv_c + F_LMPD_236_ribflv_c - R_LMPD_236_ribflv_c + F_LMPD_237_ribflv_c - R_LMPD_237_ribflv_c + F_LMPD_238_ribflv_c - R_LMPD_238_ribflv_c = 0
M_ru5p-D_c: - F_DB4PS + R_DB4PS + F_GND - R_GND - F_RPE + R_RPE + F_RPI - R_RPI -3 F_LMPD_79_colipa_e +3 R_LMPD_79_colipa_e -3 F_LMPD_80_colipa_e +3 R_LMPD_80_colipa_e -3 F_LMPD_81_colipa_e +3 R_LMPD_81_colipa_e -3 F_LMPD_82_colipa_e +3 R_LMPD_82_colipa_e -3 F_LMPD_83_colipa_e +3 R_LMPD_83_colipa_e -3 F_LMPD_84_colipa_e +3 R_LMPD_84_colipa_e -3 F_LMPD_85_colipa_e +3 R_LMPD_85_colipa_e -3 F_LMPD_86_colipa_e +3 R_LMPD_86_colipa_e -3 F_LMPD_87_colipa_e +3 R_LMPD_87_colipa_e -3 F_LMPD_88_colipa_e +3 R_LMPD_88_colipa_e -3 F_LMPD_89_colipa_e +3 R_LMPD_89_colipa_e -3 F_LMPD_90_colipa_e +3 R_LMPD_90_colipa_e = 0
M_s17bp_c: - F_FBA3 + R_FBA3 + F_PFK_3 - R_PFK_3 = 0
M_s7p_c: - F_PFK_3 + R_PFK_3 - F_TALA + R_TALA + F_TKT1 - R_TKT1 -4 F_LMPD_79_colipa_e +4 R_LMPD_79_colipa_e -4 F_LMPD_80_colipa_e +4 R_LMPD_80_colipa_e -4 F_LMPD_81_colipa_e +4 R_LMPD_81_colipa_e -4 F_LMPD_82_colipa_e +4 R_LMPD_82_colipa_e -4 F_LMPD_83_colipa_e +4 R_LMPD_83_colipa_e -4 F_LMPD_84_colipa_e +4 R_LMPD_84_colipa_e -4 F_LMPD_85_colipa_e +4 R_LMPD_85_colipa_e -4 F_LMPD_86_colipa_e +4 R_LMPD_86_colipa_e -4 F_LMPD_87_colipa_e +4 R_LMPD_87_colipa_e -4 F_LMPD_88_colipa_e +4 R_LMPD_88_colipa_e -4 F_LMPD_89_colipa_e +4 R_LMPD_89_colipa_e -4 F_LMPD_90_colipa_e +4 R_LMPD_90_colipa_e = 0
M_ser-L_c: -0.20968 F_Ec_biomass_iJO1366_WT_53p95M +0.20968 R_Ec_biomass_iJO1366_WT_53p95M + F_SERt2rpp - R_SERt2rpp + F_LMPD_239_ser-L_c - R_LMPD_239_ser-L_c = 0
M_sheme_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_240_sheme_c - R_LMPD_240_sheme_c = 0
M_so4_c: -0.004126 F_Ec_biomass_iJO1366_WT_53p95M +0.004126 R_Ec_biomass_iJO1366_WT_53p95M + F_SO4t2pp - R_SO4t2pp + F_SULabcpp - R_SULabcpp -2 F_LMPD_6_2fe2s_c +2 R_LMPD_6_2fe2s_c -2 F_LMPD_7_2fe2s_c +2 R_LMPD_7_2fe2s_c -2 F_LMPD_8_2fe2s_c +2 R_LMPD_8_2fe2s_c -2 F_LMPD_9_2fe2s_c +2 R_LMPD_9_2fe2s_c -2 F_LMPD_10_2fe2s_c +2 R_LMPD_10_2fe2s_c -2 F_LMPD_11_2fe2s_c +2 R_LMPD_11_2fe2s_c -2 F_LMPD_12_2fe2s_c +2 R_LMPD_12_2fe2s_c -2 F_LMPD_13_2fe2s_c +2 R_LMPD_13_2fe2s_c -4 F_LMPD_14_4fe4s_c +4 R_LMPD_14_4fe4s_c -4 F_LMPD_15_4fe4s_c +4 R_LMPD_15_4fe4s_c -4 F_LMPD_16_4fe4s_c +4 R_LMPD_16_4fe4s_c -4 F_LMPD_17_4fe4s_c +4 R_LMPD_17_4fe4s_c -4 F_LMPD_18_4fe4s_c +4 R_LMPD_18_4fe4s_c -4 F_LMPD_19_4fe4s_c +4 R_LMPD_19_4fe4s_c -4 F_LMPD_20_4fe4s_c +4 R_LMPD_20_4fe4s_c -4 F_LMPD_21_4fe4s_c +4 R_LMPD_21_4fe4s_c - F_LMPD_24_accoa_c + R_LMPD_24_accoa_c - F_LMPD_25_accoa_c + R_LMPD_25_accoa_c - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c -4 F_LMPD_35_bmocogdp_c +4 R_LMPD_35_bmocogdp_c -4 F_LMPD_36_bmocogdp_c +4 R_LMPD_36_bmocogdp_c -4 F_LMPD_37_bmocogdp_c +4 R_LMPD_37_bmocogdp_c -4 F_LMPD_38_bmocogdp_c +4 R_LMPD_38_bmocogdp_c -2 F_LMPD_39_btn_c +2 R_LMPD_39_btn_c -2 F_LMPD_40_btn_c +2 R_LMPD_40_btn_c -2 F_LMPD_41_btn_c +2 R_LMPD_41_btn_c -2 F_LMPD_42_btn_c +2 R_LMPD_42_btn_c -3 F_LMPD_43_btn_c +3 R_LMPD_43_btn_c -3 F_LMPD_44_btn_c +3 R_LMPD_44_btn_c -3 F_LMPD_45_btn_c +3 R_LMPD_45_btn_c -3 F_LMPD_47_btn_c +3 R_LMPD_47_btn_c - F_LMPD_92_cys-L_c + R_LMPD_92_cys-L_c - F_LMPD_93_cys-L_c + R_LMPD_93_cys-L_c - F_LMPD_118_gthrd_c + R_LMPD_118_gthrd_c - F_LMPD_119_gthrd_c + R_LMPD_119_gthrd_c -4 F_LMPD_130_lipopb_c +4 R_LMPD_130_lipopb_c -4 F_LMPD_131_lipopb_c +4 R_LMPD_131_lipopb_c -4 F_LMPD_132_lipopb_c +4 R_LMPD_132_lipopb_c -4 F_LMPD_133_lipopb_c +4 R_LMPD_133_lipopb_c - F_LMPD_135_malcoa_c + R_LMPD_135_malcoa_c - F_LMPD_136_malcoa_c + R_LMPD_136_malcoa_c - F_LMPD_137_met-L_c + R_LMPD_137_met-L_c - F_LMPD_138_met-L_c + R_LMPD_138_met-L_c -2 F_LMPD_141_mococdp_c +2 R_LMPD_141_mococdp_c -2 F_LMPD_142_mococdp_c +2 R_LMPD_142_mococdp_c -2 F_LMPD_143_mococdp_c +2 R_LMPD_143_mococdp_c -2 F_LMPD_144_mococdp_c +2 R_LMPD_144_mococdp_c -2 F_LMPD_145_mococdp_c +2 R_LMPD_145_mococdp_c -2 F_LMPD_146_mococdp_c +2 R_LMPD_146_mococdp_c -2 F_LMPD_147_mococdp_c +2 R_LMPD_147_mococdp_c -2 F_LMPD_148_mococdp_c +2 R_LMPD_148_mococdp_c -2 F_LMPD_149_mocogdp_c +2 R_LMPD_149_mocogdp_c -2 F_LMPD_150_mocogdp_c +2 R_LMPD_150_mocogdp_c -2 F_LMPD_151_mocogdp_c +2 R_LMPD_151_mocogdp_c -2 F_LMPD_152_mocogdp_c +2 R_LMPD_152_mocogdp_c - F_LMPD_241_spmd_c + R_LMPD_241_spmd_c - F_LMPD_242_spmd_c + R_LMPD_242_spmd_c - F_LMPD_247_thmpp_c + R_LMPD_247_thmpp_c - F_LMPD_248_thmpp_c + R_LMPD_248_thmpp_c = 0
M_spmd_c: -0.006744 F_Ec_biomass_iJO1366_WT_53p95M +0.006744 R_Ec_biomass_iJO1366_WT_53p95M - F_SPMDt3pp + R_SPMDt3pp + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c = 0
M_succ_c: F_ASPO5 - R_ASPO5 - F_CITt7pp + R_CITt7pp + F_DHORDfum - R_DHORDfum + F_FRD3 - R_FRD3 + F_ICL - R_ICL +3 F_PPPGO3 -3 R_PPPGO3 + F_SUCCt2_2pp - R_SUCCt2_2pp + F_SUCCt2_3pp - R_SUCCt2_3pp - F_SUCCt3pp + R_SUCCt3pp - F_SUCDi + R_SUCDi - F_SUCOAS + R_SUCOAS + F_LMPD_31_amet_c - R_LMPD_31_amet_c + F_LMPD_32_amet_c - R_LMPD_32_amet_c + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c +2 F_LMPD_43_btn_c -2 R_LMPD_43_btn_c +2 F_LMPD_44_btn_c -2 R_LMPD_44_btn_c +2 F_LMPD_45_btn_c -2 R_LMPD_45_btn_c +2 F_LMPD_47_btn_c -2 R_LMPD_47_btn_c + F_LMPD_134_lys-L_c - R_LMPD_134_lys-L_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c +2 F_LMPD_154_murein3p3p_p -2 R_LMPD_154_murein3p3p_p +2 F_LMPD_155_murein3px4p_p -2 R_LMPD_155_murein3px4p_p +2 F_LMPD_156_murein4p4p_p -2 R_LMPD_156_murein4p4p_p +2 F_LMPD_157_murein4px4p_p -2 R_LMPD_157_murein4px4p_p +3 F_LMPD_158_murein4px4px4p_p -3 R_LMPD_158_murein4px4px4p_p + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c = 0
M_succoa_c: -9.8e-05 F_Ec_biomass_iJO1366_WT_53p95M +9.8e-05 R_Ec_biomass_iJO1366_WT_53p95M + F_AKGDH - R_AKGDH + F_SUCOAS - R_SUCOAS - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c - F_LMPD_39_btn_c + R_LMPD_39_btn_c - F_LMPD_40_btn_c + R_LMPD_40_btn_c - F_LMPD_41_btn_c + R_LMPD_41_btn_c - F_LMPD_42_btn_c + R_LMPD_42_btn_c -2 F_LMPD_43_btn_c +2 R_LMPD_43_btn_c -2 F_LMPD_44_btn_c +2 R_LMPD_44_btn_c -2 F_LMPD_45_btn_c +2 R_LMPD_45_btn_c -2 F_LMPD_47_btn_c +2 R_LMPD_47_btn_c - F_LMPD_134_lys-L_c + R_LMPD_134_lys-L_c - F_LMPD_137_met-L_c + R_LMPD_137_met-L_c - F_LMPD_138_met-L_c + R_LMPD_138_met-L_c -2 F_LMPD_154_murein3p3p_p +2 R_LMPD_154_murein3p3p_p -2 F_LMPD_155_murein3px4p_p +2 R_LMPD_155_murein3px4p_p -2 F_LMPD_156_murein4p4p_p +2 R_LMPD_156_murein4p4p_p -2 F_LMPD_157_murein4px4p_p +2 R_LMPD_157_murein4px4p_p -3 F_LMPD_158_murein4px4px4p_p +3 R_LMPD_158_murein4px4px4p_p - F_LMPD_241_spmd_c + R_LMPD_241_spmd_c - F_LMPD_242_spmd_c + R_LMPD_242_spmd_c = 0
M_thf_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_243_thf_c - R_LMPD_243_thf_c + F_LMPD_244_thf_c - R_LMPD_244_thf_c + F_LMPD_245_thf_c - R_LMPD_245_thf_c + F_LMPD_246_thf_c - R_LMPD_246_thf_c = 0
M_thmpp_c: -0.000223 F_Ec_biomass_iJO1366_WT_53p95M +0.000223 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_247_thmpp_c - R_LMPD_247_thmpp_c + F_LMPD_248_thmpp_c - R_LMPD_248_thmpp_c = 0
M_thr-L_c: -0.24651 F_Ec_biomass_iJO1366_WT_53p95M +0.24651 R_Ec_biomass_iJO1366_WT_53p95M - F_THRt2pp + R_THRt2pp + F_THRt2rpp - R_THRt2rpp + F_LMPD_249_thr-L_c - R_LMPD_249_thr-L_c = 0
M_trdox_c: F_THIORDXi - R_THIORDXi - F_TRDR + R_TRDR +2 F_LMPD_6_2fe2s_c -2 R_LMPD_6_2fe2s_c +2 F_LMPD_7_2fe2s_c -2 R_LMPD_7_2fe2s_c +2 F_LMPD_8_2fe2s_c -2 R_LMPD_8_2fe2s_c +2 F_LMPD_9_2fe2s_c -2 R_LMPD_9_2fe2s_c +2 F_LMPD_10_2fe2s_c -2 R_LMPD_10_2fe2s_c +2 F_LMPD_11_2fe2s_c -2 R_LMPD_11_2fe2s_c +2 F_LMPD_12_2fe2s_c -2 R_LMPD_12_2fe2s_c +2 F_LMPD_13_2fe2s_c -2 R_LMPD_13_2fe2s_c +4 F_LMPD_14_4fe4s_c -4 R_LMPD_14_4fe4s_c +4 F_LMPD_15_4fe4s_c -4 R_LMPD_15_4fe4s_c +4 F_LMPD_16_4fe4s_c -4 R_LMPD_16_4fe4s_c +4 F_LMPD_17_4fe4s_c -4 R_LMPD_17_4fe4s_c +4 F_LMPD_18_4fe4s_c -4 R_LMPD_18_4fe4s_c +4 F_LMPD_19_4fe4s_c -4 R_LMPD_19_4fe4s_c +4 F_LMPD_20_4fe4s_c -4 R_LMPD_20_4fe4s_c +4 F_LMPD_21_4fe4s_c -4 R_LMPD_21_4fe4s_c + F_LMPD_24_accoa_c - R_LMPD_24_accoa_c + F_LMPD_25_accoa_c - R_LMPD_25_accoa_c + F_LMPD_31_amet_c - R_LMPD_31_amet_c + F_LMPD_32_amet_c - R_LMPD_32_amet_c +4 F_LMPD_35_bmocogdp_c -4 R_LMPD_35_bmocogdp_c +4 F_LMPD_36_bmocogdp_c -4 R_LMPD_36_bmocogdp_c +4 F_LMPD_37_bmocogdp_c -4 R_LMPD_37_bmocogdp_c +4 F_LMPD_38_bmocogdp_c -4 R_LMPD_38_bmocogdp_c +2 F_LMPD_39_btn_c -2 R_LMPD_39_btn_c +2 F_LMPD_40_btn_c -2 R_LMPD_40_btn_c +2 F_LMPD_41_btn_c -2 R_LMPD_41_btn_c +2 F_LMPD_42_btn_c -2 R_LMPD_42_btn_c +3 F_LMPD_43_btn_c -3 R_LMPD_43_btn_c +3 F_LMPD_44_btn_c -3 R_LMPD_44_btn_c +3 F_LMPD_45_btn_c -3 R_LMPD_45_btn_c +3 F_LMPD_47_btn_c -3 R_LMPD_47_btn_c + F_LMPD_92_cys-L_c - R_LMPD_92_cys-L_c + F_LMPD_93_cys-L_c - R_LMPD_93_cys-L_c + F_LMPD_101_dgtp_c - R_LMPD_101_dgtp_c + F_LMPD_104_dgtp_c - R_LMPD_104_dgtp_c + F_LMPD_106_dttp_c - R_LMPD_106_dttp_c + F_LMPD_118_gthrd_c - R_LMPD_118_gthrd_c + F_LMPD_119_gthrd_c - R_LMPD_119_gthrd_c +4 F_LMPD_130_lipopb_c -4 R_LMPD_130_lipopb_c +4 F_LMPD_131_lipopb_c -4 R_LMPD_131_lipopb_c +4 F_LMPD_132_lipopb_c -4 R_LMPD_132_lipopb_c +4 F_LMPD_133_lipopb_c -4 R_LMPD_133_lipopb_c + F_LMPD_135_malcoa_c - R_LMPD_135_malcoa_c + F_LMPD_136_malcoa_c - R_LMPD_136_malcoa_c + F_LMPD_137_met-L_c - R_LMPD_137_met-L_c + F_LMPD_138_met-L_c - R_LMPD_138_met-L_c +2 F_LMPD_141_mococdp_c -2 R_LMPD_141_mococdp_c +2 F_LMPD_142_mococdp_c -2 R_LMPD_142_mococdp_c +2 F_LMPD_143_mococdp_c -2 R_LMPD_143_mococdp_c +2 F_LMPD_144_mococdp_c -2 R_LMPD_144_mococdp_c +2 F_LMPD_145_mococdp_c -2 R_LMPD_145_mococdp_c +2 F_LMPD_146_mococdp_c -2 R_LMPD_146_mococdp_c +2 F_LMPD_147_mococdp_c -2 R_LMPD_147_mococdp_c +2 F_LMPD_148_mococdp_c -2 R_LMPD_148_mococdp_c +2 F_LMPD_149_mocogdp_c -2 R_LMPD_149_mocogdp_c +2 F_LMPD_150_mocogdp_c -2 R_LMPD_150_mocogdp_c +2 F_LMPD_151_mocogdp_c -2 R_LMPD_151_mocogdp_c +2 F_LMPD_152_mocogdp_c -2 R_LMPD_152_mocogdp_c + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c + F_LMPD_247_thmpp_c - R_LMPD_247_thmpp_c + F_LMPD_248_thmpp_c - R_LMPD_248_thmpp_c = 0
M_trdrd_c: - F_THIORDXi + R_THIORDXi + F_TRDR - R_TRDR -2 F_LMPD_6_2fe2s_c +2 R_LMPD_6_2fe2s_c -2 F_LMPD_7_2fe2s_c +2 R_LMPD_7_2fe2s_c -2 F_LMPD_8_2fe2s_c +2 R_LMPD_8_2fe2s_c -2 F_LMPD_9_2fe2s_c +2 R_LMPD_9_2fe2s_c -2 F_LMPD_10_2fe2s_c +2 R_LMPD_10_2fe2s_c -2 F_LMPD_11_2fe2s_c +2 R_LMPD_11_2fe2s_c -2 F_LMPD_12_2fe2s_c +2 R_LMPD_12_2fe2s_c -2 F_LMPD_13_2fe2s_c +2 R_LMPD_13_2fe2s_c -4 F_LMPD_14_4fe4s_c +4 R_LMPD_14_4fe4s_c -4 F_LMPD_15_4fe4s_c +4 R_LMPD_15_4fe4s_c -4 F_LMPD_16_4fe4s_c +4 R_LMPD_16_4fe4s_c -4 F_LMPD_17_4fe4s_c +4 R_LMPD_17_4fe4s_c -4 F_LMPD_18_4fe4s_c +4 R_LMPD_18_4fe4s_c -4 F_LMPD_19_4fe4s_c +4 R_LMPD_19_4fe4s_c -4 F_LMPD_20_4fe4s_c +4 R_LMPD_20_4fe4s_c -4 F_LMPD_21_4fe4s_c +4 R_LMPD_21_4fe4s_c - F_LMPD_24_accoa_c + R_LMPD_24_accoa_c - F_LMPD_25_accoa_c + R_LMPD_25_accoa_c - F_LMPD_31_amet_c + R_LMPD_31_amet_c - F_LMPD_32_amet_c + R_LMPD_32_amet_c -4 F_LMPD_35_bmocogdp_c +4 R_LMPD_35_bmocogdp_c -4 F_LMPD_36_bmocogdp_c +4 R_LMPD_36_bmocogdp_c -4 F_LMPD_37_bmocogdp_c +4 R_LMPD_37_bmocogdp_c -4 F_LMPD_38_bmocogdp_c +4 R_LMPD_38_bmocogdp_c -2 F_LMPD_39_btn_c +2 R_LMPD_39_btn_c -2 F_LMPD_40_btn_c +2 R_LMPD_40_btn_c -2 F_LMPD_41_btn_c +2 R_LMPD_41_btn_c -2 F_LMPD_42_btn_c +2 R_LMPD_42_btn_c -3 F_LMPD_43_btn_c +3 R_LMPD_43_btn_c -3 F_LMPD_44_btn_c +3 R_LMPD_44_btn_c -3 F_LMPD_45_btn_c +3 R_LMPD_45_btn_c -3 F_LMPD_47_btn_c +3 R_LMPD_47_btn_c - F_LMPD_92_cys-L_c + R_LMPD_92_cys-L_c - F_LMPD_93_cys-L_c + R_LMPD_93_cys-L_c - F_LMPD_101_dgtp_c + R_LMPD_101_dgtp_c - F_LMPD_104_dgtp_c + R_LMPD_104_dgtp_c - F_LMPD_106_dttp_c + R_LMPD_106_dttp_c - F_LMPD_118_gthrd_c + R_LMPD_118_gthrd_c - F_LMPD_119_gthrd_c + R_LMPD_119_gthrd_c -4 F_LMPD_130_lipopb_c +4 R_LMPD_130_lipopb_c -4 F_LMPD_131_lipopb_c +4 R_LMPD_131_lipopb_c -4 F_LMPD_132_lipopb_c +4 R_LMPD_132_lipopb_c -4 F_LMPD_133_lipopb_c +4 R_LMPD_133_lipopb_c - F_LMPD_135_malcoa_c + R_LMPD_135_malcoa_c - F_LMPD_136_malcoa_c + R_LMPD_136_malcoa_c - F_LMPD_137_met-L_c + R_LMPD_137_met-L_c - F_LMPD_138_met-L_c + R_LMPD_138_met-L_c -2 F_LMPD_141_mococdp_c +2 R_LMPD_141_mococdp_c -2 F_LMPD_142_mococdp_c +2 R_LMPD_142_mococdp_c -2 F_LMPD_143_mococdp_c +2 R_LMPD_143_mococdp_c -2 F_LMPD_144_mococdp_c +2 R_LMPD_144_mococdp_c -2 F_LMPD_145_mococdp_c +2 R_LMPD_145_mococdp_c -2 F_LMPD_146_mococdp_c +2 R_LMPD_146_mococdp_c -2 F_LMPD_147_mococdp_c +2 R_LMPD_147_mococdp_c -2 F_LMPD_148_mococdp_c +2 R_LMPD_148_mococdp_c -2 F_LMPD_149_mocogdp_c +2 R_LMPD_149_mocogdp_c -2 F_LMPD_150_mocogdp_c +2 R_LMPD_150_mocogdp_c -2 F_LMPD_151_mocogdp_c +2 R_LMPD_151_mocogdp_c -2 F_LMPD_152_mocogdp_c +2 R_LMPD_152_mocogdp_c - F_LMPD_241_spmd_c + R_LMPD_241_spmd_c - F_LMPD_242_spmd_c + R_LMPD_242_spmd_c - F_LMPD_247_thmpp_c + R_LMPD_247_thmpp_c - F_LMPD_248_thmpp_c + R_LMPD_248_thmpp_c = 0
M_trp-L_c: -0.055234 F_Ec_biomass_iJO1366_WT_53p95M +0.055234 R_Ec_biomass_iJO1366_WT_53p95M + F_TRPt2rpp - R_TRPt2rpp + F_LMPD_250_trp-L_c - R_LMPD_250_trp-L_c = 0
M_tyr-L_c: -0.13399 F_Ec_biomass_iJO1366_WT_53p95M +0.13399 R_Ec_biomass_iJO1366_WT_53p95M + F_TYRt2rpp - R_TYRt2rpp + F_LMPD_251_tyr-L_c - R_LMPD_251_tyr-L_c = 0
M_udcpdp_c: -5.5e-05 F_Ec_biomass_iJO1366_WT_53p95M +5.5e-05 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_252_udcpdp_c - R_LMPD_252_udcpdp_c = 0
M_utp_c: -0.1401 F_Ec_biomass_iJO1366_WT_53p95M +0.1401 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_253_utp_c - R_LMPD_253_utp_c = 0
M_val-L_c: -0.41118 F_Ec_biomass_iJO1366_WT_53p95M +0.41118 R_Ec_biomass_iJO1366_WT_53p95M + F_VALt2rpp - R_VALt2rpp + F_LMPD_254_val-L_c - R_LMPD_254_val-L_c = 0
M_xu5p-D_c: F_RPE - R_RPE - F_TKT1 + R_TKT1 - F_TKT2 + R_TKT2 = 0
M_zn2_c: -0.000324 F_Ec_biomass_iJO1366_WT_53p95M +0.000324 R_Ec_biomass_iJO1366_WT_53p95M + F_ZN2tpp - R_ZN2tpp + F_ZNabcpp - R_ZNabcpp = 0
M_5mtr_e: - F_DM_5mtr_e + R_DM_5mtr_e - F_5MTRtex + R_5MTRtex = 0
M_ac_e: - F_DM_ac_e + R_DM_ac_e - F_ACtex + R_ACtex = 0
M_acald_e: - F_DM_acald_e + R_DM_acald_e - F_ACALDtex + R_ACALDtex = 0
M_acser_e: - F_DM_acser_e + R_DM_acser_e - F_ACSERtex + R_ACSERtex = 0
M_akg_e: - F_DM_akg_e + R_DM_akg_e - F_AKGtex + R_AKGtex = 0
M_ala-D_e: - F_DM_ala-D_e + R_DM_ala-D_e - F_DALAtex + R_DALAtex = 0
M_ala-L_e: - F_DM_ala-L_e + R_DM_ala-L_e - F_ALAtex + R_ALAtex = 0
M_alaala_e: - F_DM_alaala_e + R_DM_alaala_e - F_ALAALAtex + R_ALAALAtex = 0
M_arg-L_e: - F_DM_arg-L_e + R_DM_arg-L_e - F_ARGtex + R_ARGtex = 0
M_asn-L_e: - F_DM_asn-L_e + R_DM_asn-L_e - F_ASNtex + R_ASNtex = 0
M_asp-L_e: - F_DM_asp-L_e + R_DM_asp-L_e - F_ASPtex + R_ASPtex = 0
M_ca2_e: - F_DM_ca2_e + R_DM_ca2_e - F_CA2tex + R_CA2tex = 0
M_cbl1_e: - F_DM_cbl1_e + R_DM_cbl1_e - F_CBL1tonex + R_CBL1tonex = 0
M_cit_e: - F_DM_cit_e + R_DM_cit_e - F_CITtex + R_CITtex = 0
M_cl_e: - F_DM_cl_e + R_DM_cl_e - F_CLtex + R_CLtex = 0
M_co2_e: - F_DM_co2_e + R_DM_co2_e - F_CO2tex + R_CO2tex = 0
M_cobalt2_e: - F_DM_cobalt2_e + R_DM_cobalt2_e - F_COBALT2tex + R_COBALT2tex = 0
M_colipa_e: -0.008151 F_Ec_biomass_iJO1366_WT_53p95M +0.008151 R_Ec_biomass_iJO1366_WT_53p95M - F_DM_colipa_e + R_DM_colipa_e + F_COLIPAabctex - R_COLIPAabctex = 0
M_cu2_e: - F_DM_cu2_e + R_DM_cu2_e - F_CU2tex + R_CU2tex = 0
M_cys-L_e: - F_DM_cys-L_e + R_DM_cys-L_e - F_CYStex + R_CYStex = 0
M_dha_e: - F_DM_dha_e + R_DM_dha_e - F_DHAtex + R_DHAtex = 0
M_enter_e: - F_DM_enter_e + R_DM_enter_e + F_FEENTERtex - R_FEENTERtex = 0
M_etoh_e: - F_DM_etoh_e + R_DM_etoh_e - F_ETOHtex + R_ETOHtex = 0
M_fe2_e: - F_DM_fe2_e + R_DM_fe2_e - F_FE2tex + R_FE2tex = 0
M_fe3_e: - F_DM_fe3_e + R_DM_fe3_e - F_FE3tex + R_FE3tex = 0
M_for_e: - F_DM_for_e + R_DM_for_e - F_FORtex + R_FORtex = 0
M_glc-D_e: - F_DM_glc_e + R_DM_glc_e - F_GLCtex + R_GLCtex - F_GLCtexi + R_GLCtexi = 0
M_glcn_e: - F_DM_glcn_e + R_DM_glcn_e - F_GLCNtex + R_GLCNtex = 0
M_glu-L_e: - F_DM_glu-L_e + R_DM_glu-L_e - F_GLUtex + R_GLUtex = 0
M_gly_e: - F_DM_gly_e + R_DM_gly_e - F_GLYtex + R_GLYtex = 0
M_glyc_e: - F_DM_glyc_e + R_DM_glyc_e - F_GLYCtex + R_GLYCtex = 0
M_glyc-R_e: - F_DM_glyc-R_e + R_DM_glyc-R_e - F_GLYCAtex + R_GLYCAtex = 0
M_glyclt_e: - F_DM_glyclt_e + R_DM_glyclt_e - F_GLYCLTtex + R_GLYCLTtex = 0
M_gthrd_e: - F_DM_gthrd_e + R_DM_gthrd_e - F_GTHRDtex + R_GTHRDtex = 0
M_h_e: - F_DM_h_e + R_DM_h_e - F_Htex + R_Htex = 0
M_h2o_e: - F_DM_h2o_e + R_DM_h2o_e - F_H2Otex + R_H2Otex = 0
M_his-L_e: - F_DM_his-L_e + R_DM_his-L_e - F_HIStex + R_HIStex = 0
M_ile-L_e: - F_DM_ile-L_e + R_DM_ile-L_e - F_ILEtex + R_ILEtex = 0
M_k_e: - F_DM_k_e + R_DM_k_e - F_Ktex + R_Ktex = 0
M_lac-D_e: - F_DM_lac-D_e + R_DM_lac-D_e - F_D-LACtex + R_D-LACtex = 0
M_leu-L_e: - F_DM_leu-L_e + R_DM_leu-L_e - F_LEUtex + R_LEUtex = 0
M_lys-L_e: - F_DM_lys-L_e + R_DM_lys-L_e - F_LYStex + R_LYStex = 0
M_mal-L_e: - F_DM_mal-L_e + R_DM_mal-L_e - F_MALtex + R_MALtex = 0
M_meoh_e: - F_DM_meoh_e + R_DM_meoh_e - F_MEOHtex + R_MEOHtex = 0
M_mg2_e: - F_DM_mg2_e + R_DM_mg2_e - F_MG2tex + R_MG2tex = 0
M_mn2_e: - F_DM_mn2_e + R_DM_mn2_e - F_MNtex + R_MNtex = 0
M_mobd_e: - F_DM_mobd_e + R_DM_mobd_e - F_MOBDtex + R_MOBDtex = 0
M_nh4_e: - F_DM_nh4_e + R_DM_nh4_e - F_NH4tex + R_NH4tex = 0
M_ni2_e: - F_DM_ni2_e + R_DM_ni2_e - F_NI2tex + R_NI2tex = 0
M_o2_e: - F_DM_o2_e + R_DM_o2_e - F_O2tex + R_O2tex = 0
M_phe-L_e: - F_DM_phe-L_e + R_DM_phe-L_e - F_PHEtex + R_PHEtex = 0
M_pheme_e: - F_DM_pheme_e + R_DM_pheme_e + F_PHEMEtiex - R_PHEMEtiex = 0
M_pi_e: - F_DM_pi_e + R_DM_pi_e - F_PItex + R_PItex = 0
M_pro-L_e: - F_DM_pro-L_e + R_DM_pro-L_e - F_PROtex + R_PROtex = 0
M_pyr_e: - F_DM_pyr_e + R_DM_pyr_e - F_PYRtex + R_PYRtex = 0
M_ser-L_e: - F_DM_ser-L_e + R_DM_ser-L_e - F_SERtex + R_SERtex = 0
M_so4_e: - F_DM_so4_e + R_DM_so4_e - F_SO4tex + R_SO4tex = 0
M_spmd_e: - F_DM_spmd_e + R_DM_spmd_e - F_SPMDtex + R_SPMDtex = 0
M_succ_e: - F_DM_succ_e + R_DM_succ_e - F_SUCCtex + R_SUCCtex = 0
M_thr-L_e: - F_DM_thr-L_e + R_DM_thr-L_e - F_THRtex + R_THRtex = 0
M_trp-L_e: - F_DM_trp-L_e + R_DM_trp-L_e - F_TRPtex + R_TRPtex = 0
M_tyr-L_e: - F_DM_tyr-L_e + R_DM_tyr-L_e - F_TYRtex + R_TYRtex = 0
M_val-L_e: - F_DM_val-L_e + R_DM_val-L_e - F_VALtex + R_VALtex = 0
M_zn2_e: - F_DM_zn2_e + R_DM_zn2_e - F_Zn2tex + R_Zn2tex = 0
M_5mtr_p: F_5MTRtex - R_5MTRtex + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c = 0
M_ac_p: - F_ACt2rpp + R_ACt2rpp + F_ACtex - R_ACtex = 0
M_acald_p: F_ACALDtex - R_ACALDtex - F_ACALDtpp + R_ACALDtpp = 0
M_acser_p: F_ACSERtex - R_ACSERtex + F_LMPD_48_ca2_c - R_LMPD_48_ca2_c = 0
M_akg_p: - F_AKGt2rpp + R_AKGt2rpp + F_AKGtex - R_AKGtex = 0
M_ala-D_p: F_DALAtex - R_DALAtex + F_LMPD_155_murein3px4p_p - R_LMPD_155_murein3px4p_p +2 F_LMPD_156_murein4p4p_p -2 R_LMPD_156_murein4p4p_p +2 F_LMPD_157_murein4px4p_p -2 R_LMPD_157_murein4px4p_p +3 F_LMPD_158_murein4px4px4p_p -3 R_LMPD_158_murein4px4px4p_p = 0
M_ala-L_p: - F_ALAt2rpp + R_ALAt2rpp + F_ALAtex - R_ALAtex = 0
M_alaala_p: F_ALAALAtex - R_ALAALAtex + F_LMPD_155_murein3px4p_p - R_LMPD_155_murein3px4p_p = 0
M_arg-L_p: F_ARGt3pp - R_ARGt3pp + F_ARGtex - R_ARGtex = 0
M_asn-L_p: - F_ASNt2rpp + R_ASNt2rpp + F_ASNtex - R_ASNtex = 0
M_asp-L_p: - F_ASPt2rpp + R_ASPt2rpp + F_ASPtex - R_ASPtex = 0
M_ca2_p: F_CA2tex - R_CA2tex + F_CAt6pp - R_CAt6pp = 0
M_cbl1_p: F_CBL1tonex - R_CBL1tonex - F_LMPD_26_adocbl_c + R_LMPD_26_adocbl_c - F_LMPD_27_adocbl_c + R_LMPD_27_adocbl_c - F_LMPD_28_adocbl_c + R_LMPD_28_adocbl_c - F_LMPD_29_adocbl_c + R_LMPD_29_adocbl_c = 0
M_cit_p: F_CITt3pp - R_CITt3pp - F_CITt7pp + R_CITt7pp + F_CITtex - R_CITtex = 0
M_cl_p: -2 F_CLt3_2pp +2 R_CLt3_2pp + F_CLtex - R_CLtex = 0
M_clpn160_p: -0.002944 F_Ec_biomass_iJO1366_WT_53p95M +0.002944 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_52_clpn160_p - R_LMPD_52_clpn160_p + F_LMPD_53_clpn160_p - R_LMPD_53_clpn160_p + F_LMPD_54_clpn160_p - R_LMPD_54_clpn160_p + F_LMPD_55_clpn160_p - R_LMPD_55_clpn160_p + F_LMPD_56_clpn160_p - R_LMPD_56_clpn160_p + F_LMPD_57_clpn160_p - R_LMPD_57_clpn160_p + F_LMPD_58_clpn160_p - R_LMPD_58_clpn160_p + F_LMPD_59_clpn160_p - R_LMPD_59_clpn160_p + F_LMPD_60_clpn160_p - R_LMPD_60_clpn160_p + F_LMPD_61_clpn160_p - R_LMPD_61_clpn160_p = 0
M_clpn161_p: -0.00229 F_Ec_biomass_iJO1366_WT_53p95M +0.00229 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_62_clpn161_p - R_LMPD_62_clpn161_p + F_LMPD_63_clpn161_p - R_LMPD_63_clpn161_p + F_LMPD_64_clpn161_p - R_LMPD_64_clpn161_p + F_LMPD_65_clpn161_p - R_LMPD_65_clpn161_p + F_LMPD_66_clpn161_p - R_LMPD_66_clpn161_p + F_LMPD_67_clpn161_p - R_LMPD_67_clpn161_p + F_LMPD_68_clpn161_p - R_LMPD_68_clpn161_p = 0
M_clpn181_p: -0.00118 F_Ec_biomass_iJO1366_WT_53p95M +0.00118 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_69_clpn181_p - R_LMPD_69_clpn181_p + F_LMPD_70_clpn181_p - R_LMPD_70_clpn181_p + F_LMPD_71_clpn181_p - R_LMPD_71_clpn181_p + F_LMPD_72_clpn181_p - R_LMPD_72_clpn181_p + F_LMPD_73_clpn181_p - R_LMPD_73_clpn181_p + F_LMPD_74_clpn181_p - R_LMPD_74_clpn181_p + F_LMPD_75_clpn181_p - R_LMPD_75_clpn181_p + F_LMPD_76_clpn181_p - R_LMPD_76_clpn181_p + F_LMPD_77_clpn181_p - R_LMPD_77_clpn181_p + F_LMPD_78_clpn181_p - R_LMPD_78_clpn181_p = 0
M_co2_p: F_CO2tex - R_CO2tex - F_CO2tpp + R_CO2tpp + F_FDH4pp - R_FDH4pp + F_FDH5pp - R_FDH5pp = 0
M_cobalt2_p: F_COBALT2tex - R_COBALT2tex - F_COBALT2tpp + R_COBALT2tpp = 0
M_colipa_p: F_COLIPAabcpp - R_COLIPAabcpp - F_COLIPAabctex + R_COLIPAabctex = 0
M_cu2_p: F_CU2tex - R_CU2tex - F_CU2tpp + R_CU2tpp = 0
M_cys-L_p: F_CYSabc2pp - R_CYSabc2pp + F_CYStex - R_CYStex + F_CYStpp - R_CYStpp = 0
M_dha_p: F_DHAtex - R_DHAtex - F_DHAtpp + R_DHAtpp = 0
M_enter_p: - F_FEENTERtex + R_FEENTERtex + F_FEENTERtpp - R_FEENTERtpp = 0
M_etoh_p: F_ETOHtex - R_ETOHtex - F_ETOHtrpp + R_ETOHtrpp = 0
M_fe2_p: - F_FE2abcpp + R_FE2abcpp - F_FE2t2pp + R_FE2t2pp + F_FE2tex - R_FE2tex - F_FE2tpp + R_FE2tpp = 0
M_fe3_p: - F_FE3abcpp + R_FE3abcpp + F_FE3tex - R_FE3tex = 0
M_for_p: - F_FDH4pp + R_FDH4pp - F_FDH5pp + R_FDH5pp + F_FORtex - R_FORtex + F_FORtppi - R_FORtppi = 0
M_glc-D_p: - F_GLCDpp + R_GLCDpp - F_GLCabcpp + R_GLCabcpp - F_GLCptspp + R_GLCptspp - F_GLCt2pp + R_GLCt2pp + F_GLCtex - R_GLCtex + F_GLCtexi - R_GLCtexi = 0
M_glcn_p: F_GLCDpp - R_GLCDpp + F_GLCNtex - R_GLCNtex = 0
M_glu-L_p: - F_GLUabcpp + R_GLUabcpp - F_GLUt2rpp + R_GLUt2rpp - F_GLUt4pp + R_GLUt4pp + F_GLUtex - R_GLUtex + F_LMPD_48_ca2_c - R_LMPD_48_ca2_c + F_LMPD_49_ca2_c - R_LMPD_49_ca2_c + F_LMPD_50_ca2_c - R_LMPD_50_ca2_c = 0
M_gly_p: - F_GLYt2rpp + R_GLYt2rpp + F_GLYtex - R_GLYtex = 0
M_glyc_p: F_GLYCtex - R_GLYCtex + F_GLYCtpp - R_GLYCtpp + F_LMPD_52_clpn160_p - R_LMPD_52_clpn160_p + F_LMPD_53_clpn160_p - R_LMPD_53_clpn160_p + F_LMPD_54_clpn160_p - R_LMPD_54_clpn160_p + F_LMPD_55_clpn160_p - R_LMPD_55_clpn160_p + F_LMPD_56_clpn160_p - R_LMPD_56_clpn160_p + F_LMPD_57_clpn160_p - R_LMPD_57_clpn160_p + F_LMPD_58_clpn160_p - R_LMPD_58_clpn160_p + F_LMPD_59_clpn160_p - R_LMPD_59_clpn160_p + F_LMPD_60_clpn160_p - R_LMPD_60_clpn160_p + F_LMPD_61_clpn160_p - R_LMPD_61_clpn160_p + F_LMPD_62_clpn161_p - R_LMPD_62_clpn161_p + F_LMPD_63_clpn161_p - R_LMPD_63_clpn161_p + F_LMPD_64_clpn161_p - R_LMPD_64_clpn161_p + F_LMPD_65_clpn161_p - R_LMPD_65_clpn161_p + F_LMPD_66_clpn161_p - R_LMPD_66_clpn161_p + F_LMPD_67_clpn161_p - R_LMPD_67_clpn161_p + F_LMPD_68_clpn161_p - R_LMPD_68_clpn161_p + F_LMPD_69_clpn181_p - R_LMPD_69_clpn181_p + F_LMPD_70_clpn181_p - R_LMPD_70_clpn181_p + F_LMPD_71_clpn181_p - R_LMPD_71_clpn181_p + F_LMPD_72_clpn181_p - R_LMPD_72_clpn181_p + F_LMPD_73_clpn181_p - R_LMPD_73_clpn181_p + F_LMPD_74_clpn181_p - R_LMPD_74_clpn181_p + F_LMPD_75_clpn181_p - R_LMPD_75_clpn181_p + F_LMPD_76_clpn181_p - R_LMPD_76_clpn181_p + F_LMPD_77_clpn181_p - R_LMPD_77_clpn181_p + F_LMPD_78_clpn181_p - R_LMPD_78_clpn181_p = 0
M_glyc-R_p: - F_GLYCAt2rpp + R_GLYCAt2rpp + F_GLYCAtex - R_GLYCAtex = 0
M_glyclt_p: - F_GLYCLTt2rpp + R_GLYCLTt2rpp + F_GLYCLTtex - R_GLYCLTtex = 0
M_gthrd_p: F_GTHRDabc2pp - R_GTHRDabc2pp + F_GTHRDtex - R_GTHRDtex - F_LMPD_48_ca2_c + R_LMPD_48_ca2_c - F_LMPD_49_ca2_c + R_LMPD_49_ca2_c - F_LMPD_50_ca2_c + R_LMPD_50_ca2_c = 0
M_h_p: - F_ACt2rpp + R_ACt2rpp - F_AKGt2rpp + R_AKGt2rpp - F_ALAt2rpp + R_ALAt2rpp - F_ARGt3pp + R_ARGt3pp - F_ASNt2rpp + R_ASNt2rpp - F_ASPt2rpp + R_ASPt2rpp -4 F_ATPS4rpp +4 R_ATPS4rpp - F_CBL1tonex + R_CBL1tonex - F_CITt3pp + R_CITt3pp + F_CLt3_2pp - R_CLt3_2pp +2 F_CYTBD2pp -2 R_CYTBD2pp +2 F_CYTBDpp -2 R_CYTBDpp +4 F_CYTBO3_4pp -4 R_CYTBO3_4pp - F_D-LACt2pp + R_D-LACt2pp + F_FDH4pp - R_FDH4pp + F_FDH5pp - R_FDH5pp - F_FE2t2pp + R_FE2t2pp - F_FEENTERtex + R_FEENTERtex - F_FEENTERtpp + R_FEENTERtpp + F_GLCDpp - R_GLCDpp - F_GLCt2pp + R_GLCt2pp - F_GLUt2rpp + R_GLUt2rpp - F_GLYCAt2rpp + R_GLYCAt2rpp - F_GLYCLTt2rpp + R_GLYCLTt2rpp - F_GLYt2rpp + R_GLYt2rpp - F_HISt2rpp + R_HISt2rpp + F_Htex - R_Htex - F_ILEt2rpp + R_ILEt2rpp - F_Kt2pp + R_Kt2pp - F_LEUt2rpp + R_LEUt2rpp - F_LYSt3pp + R_LYSt3pp - F_MALt3pp + R_MALt3pp +2 F_MG2t3_2pp -2 R_MG2t3_2pp - F_MNt2pp + R_MNt2pp +3 F_NADH16pp -3 R_NADH16pp +3 F_NADH17pp -3 R_NADH17pp +3 F_NADH18pp -3 R_NADH18pp -3 F_NAt3_1p5pp +3 R_NAt3_1p5pp -2 F_NAt3_2pp +2 R_NAt3_2pp - F_NAt3pp + R_NAt3pp - F_PHEt2rpp + R_PHEt2rpp - F_PIt2rpp + R_PIt2rpp - F_PROt2rpp + R_PROt2rpp - F_PYRt2rpp + R_PYRt2rpp - F_SERt2rpp + R_SERt2rpp - F_SO4t2pp + R_SO4t2pp - F_SPMDt3pp + R_SPMDt3pp -2 F_SUCCt2_2pp +2 R_SUCCt2_2pp -3 F_SUCCt2_3pp +3 R_SUCCt2_3pp - F_SUCCt3pp + R_SUCCt3pp -2 F_THD2pp +2 R_THD2pp - F_THRt2pp + R_THRt2pp - F_THRt2rpp + R_THRt2rpp - F_TRPt2rpp + R_TRPt2rpp - F_TYRt2rpp + R_TYRt2rpp - F_VALt2rpp + R_VALt2rpp +2 F_LMPD_43_btn_c -2 R_LMPD_43_btn_c +2 F_LMPD_44_btn_c -2 R_LMPD_44_btn_c +2 F_LMPD_45_btn_c -2 R_LMPD_45_btn_c +2 F_LMPD_47_btn_c -2 R_LMPD_47_btn_c + F_LMPD_241_spmd_c - R_LMPD_241_spmd_c + F_LMPD_242_spmd_c - R_LMPD_242_spmd_c = 0
M_h2o_p: - F_GLCDpp + R_GLCDpp + F_H2Otex - R_H2Otex - F_H2Otpp + R_H2Otpp - F_LMPD_48_ca2_c + R_LMPD_48_ca2_c - F_LMPD_49_ca2_c + R_LMPD_49_ca2_c - F_LMPD_50_ca2_c + R_LMPD_50_ca2_c -2 F_LMPD_154_murein3p3p_p +2 R_LMPD_154_murein3p3p_p - F_LMPD_155_murein3px4p_p + R_LMPD_155_murein3px4p_p -2 F_LMPD_156_murein4p4p_p +2 R_LMPD_156_murein4p4p_p - F_LMPD_157_murein4px4p_p + R_LMPD_157_murein4px4p_p - F_LMPD_158_murein4px4px4p_p + R_LMPD_158_murein4px4px4p_p = 0
M_his-L_p: - F_HISt2rpp + R_HISt2rpp + F_HIStex - R_HIStex = 0
M_ile-L_p: - F_ILEt2rpp + R_ILEt2rpp + F_ILEtex - R_ILEtex = 0
M_k_p: - F_Kabcpp + R_Kabcpp - F_Kt2pp + R_Kt2pp + F_Ktex - R_Ktex = 0
M_lac-D_p: - F_D-LACt2pp + R_D-LACt2pp + F_D-LACtex - R_D-LACtex = 0
M_leu-L_p: - F_LEUt2rpp + R_LEUt2rpp + F_LEUtex - R_LEUtex = 0
M_lys-L_p: F_LYSt3pp - R_LYSt3pp + F_LYStex - R_LYStex = 0
M_mal-L_p: F_MALt3pp - R_MALt3pp + F_MALtex - R_MALtex = 0
M_meoh_p: F_MEOHtex - R_MEOHtex + F_LMPD_39_btn_c - R_LMPD_39_btn_c + F_LMPD_40_btn_c - R_LMPD_40_btn_c + F_LMPD_41_btn_c - R_LMPD_41_btn_c + F_LMPD_42_btn_c - R_LMPD_42_btn_c + F_LMPD_43_btn_c - R_LMPD_43_btn_c + F_LMPD_44_btn_c - R_LMPD_44_btn_c + F_LMPD_45_btn_c - R_LMPD_45_btn_c + F_LMPD_47_btn_c - R_LMPD_47_btn_c = 0
M_mg2_p: - F_MG2t3_2pp + R_MG2t3_2pp + F_MG2tex - R_MG2tex - F_MG2tpp + R_MG2tpp - F_MG2uabcpp + R_MG2uabcpp = 0
M_mn2_p: - F_MN2tpp + R_MN2tpp - F_MNt2pp + R_MNt2pp + F_MNtex - R_MNtex = 0
M_mobd_p: - F_MOBDabcpp + R_MOBDabcpp + F_MOBDtex - R_MOBDtex = 0
M_murein3p3p_p: -0.001345 F_Ec_biomass_iJO1366_WT_53p95M +0.001345 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_154_murein3p3p_p - R_LMPD_154_murein3p3p_p = 0
M_murein3px4p_p: -0.000605 F_Ec_biomass_iJO1366_WT_53p95M +0.000605 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_155_murein3px4p_p - R_LMPD_155_murein3px4p_p = 0
M_murein4p4p_p: -0.005381 F_Ec_biomass_iJO1366_WT_53p95M +0.005381 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_156_murein4p4p_p - R_LMPD_156_murein4p4p_p = 0
M_murein4px4p_p: -0.005448 F_Ec_biomass_iJO1366_WT_53p95M +0.005448 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_157_murein4px4p_p - R_LMPD_157_murein4px4p_p = 0
M_murein4px4px4p_p: -0.000673 F_Ec_biomass_iJO1366_WT_53p95M +0.000673 R_Ec_biomass_iJO1366_WT_53p95M + F_LMPD_158_murein4px4px4p_p - R_LMPD_158_murein4px4px4p_p = 0
M_na1_p: - F_CAt6pp + R_CAt6pp - F_GLUt4pp + R_GLUt4pp +2 F_NAt3_1p5pp -2 R_NAt3_1p5pp + F_NAt3_2pp - R_NAt3_2pp + F_NAt3pp - R_NAt3pp = 0
M_nh4_p: F_NH4tex - R_NH4tex - F_NH4tpp + R_NH4tpp = 0
M_ni2_p: F_NI2tex - R_NI2tex - F_NI2tpp + R_NI2tpp - F_NI2uabcpp + R_NI2uabcpp = 0
M_o2_p: F_O2tex - R_O2tex - F_O2tpp + R_O2tpp = 0
M_pe160_p: -0.031798 F_Ec_biomass_iJO1366_WT_53p95M +0.031798 R_Ec_biomass_iJO1366_WT_53p95M + F_PE160abcpp - R_PE160abcpp = 0
M_pe161_p: -0.024732 F_Ec_biomass_iJO1366_WT_53p95M +0.024732 R_Ec_biomass_iJO1366_WT_53p95M + F_PE161abcpp - R_PE161abcpp = 0
M_pe181_p: -0.012747 F_Ec_biomass_iJO1366_WT_53p95M +0.012747 R_Ec_biomass_iJO1366_WT_53p95M + F_PE181abcpp - R_PE181abcpp = 0
M_pg160_p: -0.004892 F_Ec_biomass_iJO1366_WT_53p95M +0.004892 R_Ec_biomass_iJO1366_WT_53p95M + F_PG160abcpp - R_PG160abcpp -2 F_LMPD_52_clpn160_p +2 R_LMPD_52_clpn160_p -2 F_LMPD_53_clpn160_p +2 R_LMPD_53_clpn160_p -2 F_LMPD_54_clpn160_p +2 R_LMPD_54_clpn160_p -2 F_LMPD_55_clpn160_p +2 R_LMPD_55_clpn160_p -2 F_LMPD_56_clpn160_p +2 R_LMPD_56_clpn160_p -2 F_LMPD_57_clpn160_p +2 R_LMPD_57_clpn160_p -2 F_LMPD_58_clpn160_p +2 R_LMPD_58_clpn160_p -2 F_LMPD_59_clpn160_p +2 R_LMPD_59_clpn160_p -2 F_LMPD_60_clpn160_p +2 R_LMPD_60_clpn160_p -2 F_LMPD_61_clpn160_p +2 R_LMPD_61_clpn160_p = 0
M_pg161_p: -0.003805 F_Ec_biomass_iJO1366_WT_53p95M +0.003805 R_Ec_biomass_iJO1366_WT_53p95M + F_PG161abcpp - R_PG161abcpp -2 F_LMPD_62_clpn161_p +2 R_LMPD_62_clpn161_p -2 F_LMPD_63_clpn161_p +2 R_LMPD_63_clpn161_p -2 F_LMPD_64_clpn161_p +2 R_LMPD_64_clpn161_p -2 F_LMPD_65_clpn161_p +2 R_LMPD_65_clpn161_p -2 F_LMPD_66_clpn161_p +2 R_LMPD_66_clpn161_p -2 F_LMPD_67_clpn161_p +2 R_LMPD_67_clpn161_p -2 F_LMPD_68_clpn161_p +2 R_LMPD_68_clpn161_p = 0
M_pg181_p: -0.001961 F_Ec_biomass_iJO1366_WT_53p95M +0.001961 R_Ec_biomass_iJO1366_WT_53p95M + F_PG181abcpp - R_PG181abcpp -2 F_LMPD_69_clpn181_p +2 R_LMPD_69_clpn181_p -2 F_LMPD_70_clpn181_p +2 R_LMPD_70_clpn181_p -2 F_LMPD_71_clpn181_p +2 R_LMPD_71_clpn181_p -2 F_LMPD_72_clpn181_p +2 R_LMPD_72_clpn181_p -2 F_LMPD_73_clpn181_p +2 R_LMPD_73_clpn181_p -2 F_LMPD_74_clpn181_p +2 R_LMPD_74_clpn181_p -2 F_LMPD_75_clpn181_p +2 R_LMPD_75_clpn181_p -2 F_LMPD_76_clpn181_p +2 R_LMPD_76_clpn181_p -2 F_LMPD_77_clpn181_p +2 R_LMPD_77_clpn181_p -2 F_LMPD_78_clpn181_p +2 R_LMPD_78_clpn181_p = 0
M_phe-L_p: - F_PHEt2rpp + R_PHEt2rpp + F_PHEtex - R_PHEtex = 0
M_pheme_p: F_PHEMEabcpp - R_PHEMEabcpp - F_PHEMEtiex + R_PHEMEtiex = 0
M_pi_p: - F_PIt2rpp + R_PIt2rpp + F_PItex - R_PItex - F_PIuabcpp + R_PIuabcpp = 0
M_pro-L_p: - F_PROt2rpp + R_PROt2rpp + F_PROtex - R_PROtex = 0
M_pyr_p: - F_PYRt2rpp + R_PYRt2rpp + F_PYRtex - R_PYRtex = 0
M_ser-L_p: - F_SERt2rpp + R_SERt2rpp + F_SERtex - R_SERtex = 0
M_so4_p: - F_SO4t2pp + R_SO4t2pp + F_SO4tex - R_SO4tex - F_SULabcpp + R_SULabcpp = 0
M_spmd_p: F_SPMDt3pp - R_SPMDt3pp + F_SPMDtex - R_SPMDtex = 0
M_succ_p: F_CITt7pp - R_CITt7pp - F_SUCCt2_2pp + R_SUCCt2_2pp - F_SUCCt2_3pp + R_SUCCt2_3pp + F_SUCCt3pp - R_SUCCt3pp + F_SUCCtex - R_SUCCtex = 0
M_thr-L_p: F_THRt2pp - R_THRt2pp - F_THRt2rpp + R_THRt2rpp + F_THRtex - R_THRtex = 0
M_trp-L_p: - F_TRPt2rpp + R_TRPt2rpp + F_TRPtex - R_TRPtex = 0
M_tyr-L_p: - F_TYRt2rpp + R_TYRt2rpp + F_TYRtex - R_TYRtex = 0
M_val-L_p: - F_VALt2rpp + R_VALt2rpp + F_VALtex - R_VALtex = 0
M_zn2_p: - F_ZN2tpp + R_ZN2tpp - F_ZNabcpp + R_ZNabcpp + F_Zn2tex - R_Zn2tex = 0
M_ade_p: +2 F_LMPD_43_btn_c -2 R_LMPD_43_btn_c +2 F_LMPD_44_btn_c -2 R_LMPD_44_btn_c +2 F_LMPD_45_btn_c -2 R_LMPD_45_btn_c +2 F_LMPD_47_btn_c -2 R_LMPD_47_btn_c = 0
SU_DM_4CRSOL: FU_DM_4CRSOL + BU_DM_4CRSOL < 1
UF_DM_4CRSOL: F_DM_4CRSOL -1000 FU_DM_4CRSOL < 0
UR_DM_4CRSOL: R_DM_4CRSOL -1000 BU_DM_4CRSOL < 0
SU_DM_5DRIB: FU_DM_5DRIB + BU_DM_5DRIB < 1
UF_DM_5DRIB: F_DM_5DRIB -1000 FU_DM_5DRIB < 0
UR_DM_5DRIB: R_DM_5DRIB -1000 BU_DM_5DRIB < 0
SU_DM_AMOB: FU_DM_AMOB + BU_DM_AMOB < 1
UF_DM_AMOB: F_DM_AMOB -1000 FU_DM_AMOB < 0
UR_DM_AMOB: R_DM_AMOB -1000 BU_DM_AMOB < 0
SU_DM_MTHTHF: FU_DM_MTHTHF + BU_DM_MTHTHF < 1
UF_DM_MTHTHF: F_DM_MTHTHF -1000 FU_DM_MTHTHF < 0
UR_DM_MTHTHF: R_DM_MTHTHF -1000 BU_DM_MTHTHF < 0
SU_Ec_biomass_iJO1366_WT_53p95M: FU_Ec_biomass_iJO1366_WT_53p95M + BU_Ec_biomass_iJO1366_WT_53p95M < 1
UF_Ec_biomass_iJO1366_WT_53p95M: F_Ec_biomass_iJO1366_WT_53p95M -1000 FU_Ec_biomass_iJO1366_WT_53p95M < 0
UR_Ec_biomass_iJO1366_WT_53p95M: R_Ec_biomass_iJO1366_WT_53p95M -1000 BU_Ec_biomass_iJO1366_WT_53p95M < 0
SU_DM_5mtr_e: FU_DM_5mtr_e + BU_DM_5mtr_e < 1
UF_DM_5mtr_e: F_DM_5mtr_e -1000 FU_DM_5mtr_e < 0
UR_DM_5mtr_e: R_DM_5mtr_e -1000 BU_DM_5mtr_e < 0
SU_DM_ac_e: FU_DM_ac_e + BU_DM_ac_e < 1
UF_DM_ac_e: F_DM_ac_e -1000 FU_DM_ac_e < 0
UR_DM_ac_e: R_DM_ac_e -1000 BU_DM_ac_e < 0
SU_DM_acald_e: FU_DM_acald_e + BU_DM_acald_e < 1
UF_DM_acald_e: F_DM_acald_e -1000 FU_DM_acald_e < 0
UR_DM_acald_e: R_DM_acald_e -1000 BU_DM_acald_e < 0
SU_DM_acser_e: FU_DM_acser_e + BU_DM_acser_e < 1
UF_DM_acser_e: F_DM_acser_e -1000 FU_DM_acser_e < 0
UR_DM_acser_e: R_DM_acser_e -1000 BU_DM_acser_e < 0
SU_DM_akg_e: FU_DM_akg_e + BU_DM_akg_e < 1
UF_DM_akg_e: F_DM_akg_e -1000 FU_DM_akg_e < 0
UR_DM_akg_e: R_DM_akg_e -1000 BU_DM_akg_e < 0
SU_DM_ala-D_e: FU_DM_ala-D_e + BU_DM_ala-D_e < 1
UF_DM_ala-D_e: F_DM_ala-D_e -1000 FU_DM_ala-D_e < 0
UR_DM_ala-D_e: R_DM_ala-D_e -1000 BU_DM_ala-D_e < 0
SU_DM_ala-L_e: FU_DM_ala-L_e + BU_DM_ala-L_e < 1
UF_DM_ala-L_e: F_DM_ala-L_e -1000 FU_DM_ala-L_e < 0
UR_DM_ala-L_e: R_DM_ala-L_e -1000 BU_DM_ala-L_e < 0
SU_DM_alaala_e: FU_DM_alaala_e + BU_DM_alaala_e < 1
UF_DM_alaala_e: F_DM_alaala_e -1000 FU_DM_alaala_e < 0
UR_DM_alaala_e: R_DM_alaala_e -1000 BU_DM_alaala_e < 0
SU_DM_arg-L_e: FU_DM_arg-L_e + BU_DM_arg-L_e < 1
UF_DM_arg-L_e: F_DM_arg-L_e -1000 FU_DM_arg-L_e < 0
UR_DM_arg-L_e: R_DM_arg-L_e -1000 BU_DM_arg-L_e < 0
SU_DM_asn-L_e: FU_DM_asn-L_e + BU_DM_asn-L_e < 1
UF_DM_asn-L_e: F_DM_asn-L_e -1000 FU_DM_asn-L_e < 0
UR_DM_asn-L_e: R_DM_asn-L_e -1000 BU_DM_asn-L_e < 0
SU_DM_asp-L_e: FU_DM_asp-L_e + BU_DM_asp-L_e < 1
UF_DM_asp-L_e: F_DM_asp-L_e -1000 FU_DM_asp-L_e < 0
UR_DM_asp-L_e: R_DM_asp-L_e -1000 BU_DM_asp-L_e < 0
SU_DM_ca2_e: FU_DM_ca2_e + BU_DM_ca2_e < 1
UF_DM_ca2_e: F_DM_ca2_e -1000 FU_DM_ca2_e < 0
UR_DM_ca2_e: R_DM_ca2_e -1000 BU_DM_ca2_e < 0
SU_DM_cbl1_e: FU_DM_cbl1_e + BU_DM_cbl1_e < 1
UF_DM_cbl1_e: F_DM_cbl1_e -1000 FU_DM_cbl1_e < 0
UR_DM_cbl1_e: R_DM_cbl1_e -1000 BU_DM_cbl1_e < 0
SU_DM_cit_e: FU_DM_cit_e + BU_DM_cit_e < 1
UF_DM_cit_e: F_DM_cit_e -1000 FU_DM_cit_e < 0
UR_DM_cit_e: R_DM_cit_e -1000 BU_DM_cit_e < 0
SU_DM_cl_e: FU_DM_cl_e + BU_DM_cl_e < 1
UF_DM_cl_e: F_DM_cl_e -1000 FU_DM_cl_e < 0
UR_DM_cl_e: R_DM_cl_e -1000 BU_DM_cl_e < 0
SU_DM_co2_e: FU_DM_co2_e + BU_DM_co2_e < 1
UF_DM_co2_e: F_DM_co2_e -1000 FU_DM_co2_e < 0
UR_DM_co2_e: R_DM_co2_e -1000 BU_DM_co2_e < 0
SU_DM_cobalt2_e: FU_DM_cobalt2_e + BU_DM_cobalt2_e < 1
UF_DM_cobalt2_e: F_DM_cobalt2_e -1000 FU_DM_cobalt2_e < 0
UR_DM_cobalt2_e: R_DM_cobalt2_e -1000 BU_DM_cobalt2_e < 0
SU_DM_colipa_e: FU_DM_colipa_e + BU_DM_colipa_e < 1
UF_DM_colipa_e: F_DM_colipa_e -1000 FU_DM_colipa_e < 0
UR_DM_colipa_e: R_DM_colipa_e -1000 BU_DM_colipa_e < 0
SU_DM_cu2_e: FU_DM_cu2_e + BU_DM_cu2_e < 1
UF_DM_cu2_e: F_DM_cu2_e -1000 FU_DM_cu2_e < 0
UR_DM_cu2_e: R_DM_cu2_e -1000 BU_DM_cu2_e < 0
SU_DM_cys-L_e: FU_DM_cys-L_e + BU_DM_cys-L_e < 1
UF_DM_cys-L_e: F_DM_cys-L_e -1000 FU_DM_cys-L_e < 0
UR_DM_cys-L_e: R_DM_cys-L_e -1000 BU_DM_cys-L_e < 0
SU_DM_dha_e: FU_DM_dha_e + BU_DM_dha_e < 1
UF_DM_dha_e: F_DM_dha_e -1000 FU_DM_dha_e < 0
UR_DM_dha_e: R_DM_dha_e -1000 BU_DM_dha_e < 0
SU_DM_enter_e: FU_DM_enter_e + BU_DM_enter_e < 1
UF_DM_enter_e: F_DM_enter_e -1000 FU_DM_enter_e < 0
UR_DM_enter_e: R_DM_enter_e -1000 BU_DM_enter_e < 0
SU_DM_etoh_e: FU_DM_etoh_e + BU_DM_etoh_e < 1
UF_DM_etoh_e: F_DM_etoh_e -1000 FU_DM_etoh_e < 0
UR_DM_etoh_e: R_DM_etoh_e -1000 BU_DM_etoh_e < 0
SU_DM_fe2_e: FU_DM_fe2_e + BU_DM_fe2_e < 1
UF_DM_fe2_e: F_DM_fe2_e -1000 FU_DM_fe2_e < 0
UR_DM_fe2_e: R_DM_fe2_e -1000 BU_DM_fe2_e < 0
SU_DM_fe3_e: FU_DM_fe3_e + BU_DM_fe3_e < 1
UF_DM_fe3_e: F_DM_fe3_e -1000 FU_DM_fe3_e < 0
UR_DM_fe3_e: R_DM_fe3_e -1000 BU_DM_fe3_e < 0
SU_DM_for_e: FU_DM_for_e + BU_DM_for_e < 1
UF_DM_for_e: F_DM_for_e -1000 FU_DM_for_e < 0
UR_DM_for_e: R_DM_for_e -1000 BU_DM_for_e < 0
SU_DM_glc_e: FU_DM_glc_e + BU_DM_glc_e < 1
UF_DM_glc_e: F_DM_glc_e -1000 FU_DM_glc_e < 0
UR_DM_glc_e: R_DM_glc_e -1000 BU_DM_glc_e < 0
SU_DM_glcn_e: FU_DM_glcn_e + BU_DM_glcn_e < 1
UF_DM_glcn_e: F_DM_glcn_e -1000 FU_DM_glcn_e < 0
UR_DM_glcn_e: R_DM_glcn_e -1000 BU_DM_glcn_e < 0
SU_DM_glu-L_e: FU_DM_glu-L_e + BU_DM_glu-L_e < 1
UF_DM_glu-L_e: F_DM_glu-L_e -1000 FU_DM_glu-L_e < 0
UR_DM_glu-L_e: R_DM_glu-L_e -1000 BU_DM_glu-L_e < 0
SU_DM_gly_e: FU_DM_gly_e + BU_DM_gly_e < 1
UF_DM_gly_e: F_DM_gly_e -1000 FU_DM_gly_e < 0
UR_DM_gly_e: R_DM_gly_e -1000 BU_DM_gly_e < 0
SU_DM_glyc_e: FU_DM_glyc_e + BU_DM_glyc_e < 1
UF_DM_glyc_e: F_DM_glyc_e -1000 FU_DM_glyc_e < 0
UR_DM_glyc_e: R_DM_glyc_e -1000 BU_DM_glyc_e < 0
SU_DM_glyc-R_e: FU_DM_glyc-R_e + BU_DM_glyc-R_e < 1
UF_DM_glyc-R_e: F_DM_glyc-R_e -1000 FU_DM_glyc-R_e < 0
UR_DM_glyc-R_e: R_DM_glyc-R_e -1000 BU_DM_glyc-R_e < 0
SU_DM_glycogenn1_c: FU_DM_glycogenn1_c + BU_DM_glycogenn1_c < 1
UF_DM_glycogenn1_c: F_DM_glycogenn1_c -1000 FU_DM_glycogenn1_c < 0
UR_DM_glycogenn1_c: R_DM_glycogenn1_c -1000 BU_DM_glycogenn1_c < 0
SU_DM_glyclt_e: FU_DM_glyclt_e + BU_DM_glyclt_e < 1
UF_DM_glyclt_e: F_DM_glyclt_e -1000 FU_DM_glyclt_e < 0
UR_DM_glyclt_e: R_DM_glyclt_e -1000 BU_DM_glyclt_e < 0
SU_DM_gthrd_e: FU_DM_gthrd_e + BU_DM_gthrd_e < 1
UF_DM_gthrd_e: F_DM_gthrd_e -1000 FU_DM_gthrd_e < 0
UR_DM_gthrd_e: R_DM_gthrd_e -1000 BU_DM_gthrd_e < 0
SU_DM_h_e: FU_DM_h_e + BU_DM_h_e < 1
UF_DM_h_e: F_DM_h_e -1000 FU_DM_h_e < 0
UR_DM_h_e: R_DM_h_e -1000 BU_DM_h_e < 0
SU_DM_h2o_e: FU_DM_h2o_e + BU_DM_h2o_e < 1
UF_DM_h2o_e: F_DM_h2o_e -1000 FU_DM_h2o_e < 0
UR_DM_h2o_e: R_DM_h2o_e -1000 BU_DM_h2o_e < 0
SU_DM_his-L_e: FU_DM_his-L_e + BU_DM_his-L_e < 1
UF_DM_his-L_e: F_DM_his-L_e -1000 FU_DM_his-L_e < 0
UR_DM_his-L_e: R_DM_his-L_e -1000 BU_DM_his-L_e < 0
SU_DM_ile-L_e: FU_DM_ile-L_e + BU_DM_ile-L_e < 1
UF_DM_ile-L_e: F_DM_ile-L_e -1000 FU_DM_ile-L_e < 0
UR_DM_ile-L_e: R_DM_ile-L_e -1000 BU_DM_ile-L_e < 0
SU_DM_k_e: FU_DM_k_e + BU_DM_k_e < 1
UF_DM_k_e: F_DM_k_e -1000 FU_DM_k_e < 0
UR_DM_k_e: R_DM_k_e -1000 BU_DM_k_e < 0
SU_DM_lac-D_e: FU_DM_lac-D_e + BU_DM_lac-D_e < 1
UF_DM_lac-D_e: F_DM_lac-D_e -1000 FU_DM_lac-D_e < 0
UR_DM_lac-D_e: R_DM_lac-D_e -1000 BU_DM_lac-D_e < 0
SU_DM_leu-L_e: FU_DM_leu-L_e + BU_DM_leu-L_e < 1
UF_DM_leu-L_e: F_DM_leu-L_e -1000 FU_DM_leu-L_e < 0
UR_DM_leu-L_e: R_DM_leu-L_e -1000 BU_DM_leu-L_e < 0
SU_DM_lys-L_e: FU_DM_lys-L_e + BU_DM_lys-L_e < 1
UF_DM_lys-L_e: F_DM_lys-L_e -1000 FU_DM_lys-L_e < 0
UR_DM_lys-L_e: R_DM_lys-L_e -1000 BU_DM_lys-L_e < 0
SU_DM_mal-L_e: FU_DM_mal-L_e + BU_DM_mal-L_e < 1
UF_DM_mal-L_e: F_DM_mal-L_e -1000 FU_DM_mal-L_e < 0
UR_DM_mal-L_e: R_DM_mal-L_e -1000 BU_DM_mal-L_e < 0
SU_DM_meoh_e: FU_DM_meoh_e + BU_DM_meoh_e < 1
UF_DM_meoh_e: F_DM_meoh_e -1000 FU_DM_meoh_e < 0
UR_DM_meoh_e: R_DM_meoh_e -1000 BU_DM_meoh_e < 0
SU_DM_mg2_e: FU_DM_mg2_e + BU_DM_mg2_e < 1
UF_DM_mg2_e: F_DM_mg2_e -1000 FU_DM_mg2_e < 0
UR_DM_mg2_e: R_DM_mg2_e -1000 BU_DM_mg2_e < 0
SU_DM_mn2_e: FU_DM_mn2_e + BU_DM_mn2_e < 1
UF_DM_mn2_e: F_DM_mn2_e -1000 FU_DM_mn2_e < 0
UR_DM_mn2_e: R_DM_mn2_e -1000 BU_DM_mn2_e < 0
SU_DM_mobd_e: FU_DM_mobd_e + BU_DM_mobd_e < 1
UF_DM_mobd_e: F_DM_mobd_e -1000 FU_DM_mobd_e < 0
UR_DM_mobd_e: R_DM_mobd_e -1000 BU_DM_mobd_e < 0
SU_DM_nh4_e: FU_DM_nh4_e + BU_DM_nh4_e < 1
UF_DM_nh4_e: F_DM_nh4_e -1000 FU_DM_nh4_e < 0
UR_DM_nh4_e: R_DM_nh4_e -1000 BU_DM_nh4_e < 0
SU_DM_ni2_e: FU_DM_ni2_e + BU_DM_ni2_e < 1
UF_DM_ni2_e: F_DM_ni2_e -1000 FU_DM_ni2_e < 0
UR_DM_ni2_e: R_DM_ni2_e -1000 BU_DM_ni2_e < 0
SU_DM_o2_e: FU_DM_o2_e + BU_DM_o2_e < 1
UF_DM_o2_e: F_DM_o2_e -1000 FU_DM_o2_e < 0
UR_DM_o2_e: R_DM_o2_e -1000 BU_DM_o2_e < 0
SU_DM_phe-L_e: FU_DM_phe-L_e + BU_DM_phe-L_e < 1
UF_DM_phe-L_e: F_DM_phe-L_e -1000 FU_DM_phe-L_e < 0
UR_DM_phe-L_e: R_DM_phe-L_e -1000 BU_DM_phe-L_e < 0
SU_DM_pheme_e: FU_DM_pheme_e + BU_DM_pheme_e < 1
UF_DM_pheme_e: F_DM_pheme_e -1000 FU_DM_pheme_e < 0
UR_DM_pheme_e: R_DM_pheme_e -1000 BU_DM_pheme_e < 0
SU_DM_pi_e: FU_DM_pi_e + BU_DM_pi_e < 1
UF_DM_pi_e: F_DM_pi_e -1000 FU_DM_pi_e < 0
UR_DM_pi_e: R_DM_pi_e -1000 BU_DM_pi_e < 0
SU_DM_pro-L_e: FU_DM_pro-L_e + BU_DM_pro-L_e < 1
UF_DM_pro-L_e: F_DM_pro-L_e -1000 FU_DM_pro-L_e < 0
UR_DM_pro-L_e: R_DM_pro-L_e -1000 BU_DM_pro-L_e < 0
SU_DM_pyr_e: FU_DM_pyr_e + BU_DM_pyr_e < 1
UF_DM_pyr_e: F_DM_pyr_e -1000 FU_DM_pyr_e < 0
UR_DM_pyr_e: R_DM_pyr_e -1000 BU_DM_pyr_e < 0
SU_DM_ser-L_e: FU_DM_ser-L_e + BU_DM_ser-L_e < 1
UF_DM_ser-L_e: F_DM_ser-L_e -1000 FU_DM_ser-L_e < 0
UR_DM_ser-L_e: R_DM_ser-L_e -1000 BU_DM_ser-L_e < 0
SU_DM_so4_e: FU_DM_so4_e + BU_DM_so4_e < 1
UF_DM_so4_e: F_DM_so4_e -1000 FU_DM_so4_e < 0
UR_DM_so4_e: R_DM_so4_e -1000 BU_DM_so4_e < 0
SU_DM_spmd_e: FU_DM_spmd_e + BU_DM_spmd_e < 1
UF_DM_spmd_e: F_DM_spmd_e -1000 FU_DM_spmd_e < 0
UR_DM_spmd_e: R_DM_spmd_e -1000 BU_DM_spmd_e < 0
SU_DM_succ_e: FU_DM_succ_e + BU_DM_succ_e < 1
UF_DM_succ_e: F_DM_succ_e -1000 FU_DM_succ_e < 0
UR_DM_succ_e: R_DM_succ_e -1000 BU_DM_succ_e < 0
SU_DM_thr-L_e: FU_DM_thr-L_e + BU_DM_thr-L_e < 1
UF_DM_thr-L_e: F_DM_thr-L_e -1000 FU_DM_thr-L_e < 0
UR_DM_thr-L_e: R_DM_thr-L_e -1000 BU_DM_thr-L_e < 0
SU_DM_trp-L_e: FU_DM_trp-L_e + BU_DM_trp-L_e < 1
UF_DM_trp-L_e: F_DM_trp-L_e -1000 FU_DM_trp-L_e < 0
UR_DM_trp-L_e: R_DM_trp-L_e -1000 BU_DM_trp-L_e < 0
SU_DM_tyr-L_e: FU_DM_tyr-L_e + BU_DM_tyr-L_e < 1
UF_DM_tyr-L_e: F_DM_tyr-L_e -1000 FU_DM_tyr-L_e < 0
UR_DM_tyr-L_e: R_DM_tyr-L_e -1000 BU_DM_tyr-L_e < 0
SU_DM_val-L_e: FU_DM_val-L_e + BU_DM_val-L_e < 1
UF_DM_val-L_e: F_DM_val-L_e -1000 FU_DM_val-L_e < 0
UR_DM_val-L_e: R_DM_val-L_e -1000 BU_DM_val-L_e < 0
SU_DM_zn2_e: FU_DM_zn2_e + BU_DM_zn2_e < 1
UF_DM_zn2_e: F_DM_zn2_e -1000 FU_DM_zn2_e < 0
UR_DM_zn2_e: R_DM_zn2_e -1000 BU_DM_zn2_e < 0
G_5MTRtex: -0.59209 LC_5mtr_e +0.59209 LC_5mtr_p - DG_5MTRtex + DGo_5MTRtex = 0
FU_5MTRtex: DG_5MTRtex +1000 FU_5MTRtex < 1000
BU_5MTRtex: - DG_5MTRtex +1000 BU_5MTRtex < 1000
SU_5MTRtex: FU_5MTRtex + BU_5MTRtex < 1
UF_5MTRtex: F_5MTRtex -1000 FU_5MTRtex < 0
UR_5MTRtex: R_5MTRtex -1000 BU_5MTRtex < 0
G_ACALD: -0.59209 LC_acald_c +0.59209 LC_accoa_c -0.59209 LC_coa_c -0.59209 LC_nad_c +0.59209 LC_nadh_c - DG_ACALD + DGo_ACALD = 0
FU_ACALD: DG_ACALD +1000 FU_ACALD < 1000
BU_ACALD: - DG_ACALD +1000 BU_ACALD < 1000
SU_ACALD: FU_ACALD + BU_ACALD < 1
UF_ACALD: F_ACALD -1000 FU_ACALD < 0
UR_ACALD: R_ACALD -1000 BU_ACALD < 0
G_ACALDtex: -0.59209 LC_acald_e +0.59209 LC_acald_p - DG_ACALDtex + DGo_ACALDtex = 0
FU_ACALDtex: DG_ACALDtex +1000 FU_ACALDtex < 1000
BU_ACALDtex: - DG_ACALDtex +1000 BU_ACALDtex < 1000
SU_ACALDtex: FU_ACALDtex + BU_ACALDtex < 1
UF_ACALDtex: F_ACALDtex -1000 FU_ACALDtex < 0
UR_ACALDtex: R_ACALDtex -1000 BU_ACALDtex < 0
G_ACALDtpp: +0.59209 LC_acald_c -0.59209 LC_acald_p - DG_ACALDtpp + DGo_ACALDtpp = 0
FU_ACALDtpp: DG_ACALDtpp +1000 FU_ACALDtpp < 1000
BU_ACALDtpp: - DG_ACALDtpp +1000 BU_ACALDtpp < 1000
SU_ACALDtpp: FU_ACALDtpp + BU_ACALDtpp < 1
UF_ACALDtpp: F_ACALDtpp -1000 FU_ACALDtpp < 0
UR_ACALDtpp: R_ACALDtpp -1000 BU_ACALDtpp < 0
G_ACKr: -0.59209 LC_ac_c +0.59209 LC_actp_c +0.59209 LC_adp_c -0.59209 LC_atp_c - DG_ACKr + DGo_ACKr = 0
FU_ACKr: DG_ACKr +1000 FU_ACKr < 1000
BU_ACKr: - DG_ACKr +1000 BU_ACKr < 1000
SU_ACKr: FU_ACKr + BU_ACKr < 1
UF_ACKr: F_ACKr -1000 FU_ACKr < 0
UR_ACKr: R_ACKr -1000 BU_ACKr < 0
G_ACONTa: +0.59209 LC_acon-C_c -0.59209 LC_cit_c - DG_ACONTa + DGo_ACONTa = 0
FU_ACONTa: DG_ACONTa +1000 FU_ACONTa < 1000
BU_ACONTa: - DG_ACONTa +1000 BU_ACONTa < 1000
SU_ACONTa: FU_ACONTa + BU_ACONTa < 1
UF_ACONTa: F_ACONTa -1000 FU_ACONTa < 0
UR_ACONTa: R_ACONTa -1000 BU_ACONTa < 0
G_ACONTb: -0.59209 LC_acon-C_c +0.59209 LC_icit_c - DG_ACONTb + DGo_ACONTb = 0
FU_ACONTb: DG_ACONTb +1000 FU_ACONTb < 1000
BU_ACONTb: - DG_ACONTb +1000 BU_ACONTb < 1000
SU_ACONTb: FU_ACONTb + BU_ACONTb < 1
UF_ACONTb: F_ACONTb -1000 FU_ACONTb < 0
UR_ACONTb: R_ACONTb -1000 BU_ACONTb < 0
G_ACS: -0.59209 LC_ac_c +0.59209 LC_accoa_c +0.59209 LC_amp_c -0.59209 LC_atp_c -0.59209 LC_coa_c +0.59209 LC_ppi_c - DG_ACS + DGo_ACS = 0
FU_ACS: DG_ACS +1000 FU_ACS < 1000
BU_ACS: - DG_ACS +1000 BU_ACS < 1000
SU_ACS: FU_ACS + BU_ACS < 1
UF_ACS: F_ACS -1000 FU_ACS < 0
UR_ACS: R_ACS -1000 BU_ACS < 0
G_ACSERtex: -0.59209 LC_acser_e +0.59209 LC_acser_p - DG_ACSERtex + DGo_ACSERtex = 0
FU_ACSERtex: DG_ACSERtex +1000 FU_ACSERtex < 1000
BU_ACSERtex: - DG_ACSERtex +1000 BU_ACSERtex < 1000
SU_ACSERtex: FU_ACSERtex + BU_ACSERtex < 1
UF_ACSERtex: F_ACSERtex -1000 FU_ACSERtex < 0
UR_ACSERtex: R_ACSERtex -1000 BU_ACSERtex < 0
G_ACt2rpp: +0.59209 LC_ac_c -0.59209 LC_ac_p - DG_ACt2rpp + DGo_ACt2rpp = 0
FU_ACt2rpp: DG_ACt2rpp +1000 FU_ACt2rpp < 1000
BU_ACt2rpp: - DG_ACt2rpp +1000 BU_ACt2rpp < 1000
SU_ACt2rpp: FU_ACt2rpp + BU_ACt2rpp < 1
UF_ACt2rpp: F_ACt2rpp -1000 FU_ACt2rpp < 0
UR_ACt2rpp: R_ACt2rpp -1000 BU_ACt2rpp < 0
G_ACtex: -0.59209 LC_ac_e +0.59209 LC_ac_p - DG_ACtex + DGo_ACtex = 0
FU_ACtex: DG_ACtex +1000 FU_ACtex < 1000
BU_ACtex: - DG_ACtex +1000 BU_ACtex < 1000
SU_ACtex: FU_ACtex + BU_ACtex < 1
UF_ACtex: F_ACtex -1000 FU_ACtex < 0
UR_ACtex: R_ACtex -1000 BU_ACtex < 0
G_ADK1: +1.1842 LC_adp_c -0.59209 LC_amp_c -0.59209 LC_atp_c - DG_ADK1 + DGo_ADK1 = 0
FU_ADK1: DG_ADK1 +1000 FU_ADK1 < 1000
BU_ADK1: - DG_ADK1 +1000 BU_ADK1 < 1000
SU_ADK1: FU_ADK1 + BU_ADK1 < 1
UF_ADK1: F_ADK1 -1000 FU_ADK1 < 0
UR_ADK1: R_ADK1 -1000 BU_ADK1 < 0
G_AKGDH: -0.59209 LC_akg_c +0.59209 LC_co2_c -0.59209 LC_coa_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_succoa_c - DG_AKGDH + DGo_AKGDH = 0
FU_AKGDH: DG_AKGDH +1000 FU_AKGDH < 1000
BU_AKGDH: - DG_AKGDH +1000 BU_AKGDH < 1000
SU_AKGDH: FU_AKGDH + BU_AKGDH < 1
UF_AKGDH: F_AKGDH -1000 FU_AKGDH < 0
UR_AKGDH: R_AKGDH -1000 BU_AKGDH < 0
G_AKGt2rpp: +0.59209 LC_akg_c -0.59209 LC_akg_p - DG_AKGt2rpp + DGo_AKGt2rpp = 0
FU_AKGt2rpp: DG_AKGt2rpp +1000 FU_AKGt2rpp < 1000
BU_AKGt2rpp: - DG_AKGt2rpp +1000 BU_AKGt2rpp < 1000
SU_AKGt2rpp: FU_AKGt2rpp + BU_AKGt2rpp < 1
UF_AKGt2rpp: F_AKGt2rpp -1000 FU_AKGt2rpp < 0
UR_AKGt2rpp: R_AKGt2rpp -1000 BU_AKGt2rpp < 0
G_AKGtex: -0.59209 LC_akg_e +0.59209 LC_akg_p - DG_AKGtex + DGo_AKGtex = 0
FU_AKGtex: DG_AKGtex +1000 FU_AKGtex < 1000
BU_AKGtex: - DG_AKGtex +1000 BU_AKGtex < 1000
SU_AKGtex: FU_AKGtex + BU_AKGtex < 1
UF_AKGtex: F_AKGtex -1000 FU_AKGtex < 0
UR_AKGtex: R_AKGtex -1000 BU_AKGtex < 0
G_ALAALAtex: -0.59209 LC_alaala_e +0.59209 LC_alaala_p - DG_ALAALAtex + DGo_ALAALAtex = 0
FU_ALAALAtex: DG_ALAALAtex +1000 FU_ALAALAtex < 1000
BU_ALAALAtex: - DG_ALAALAtex +1000 BU_ALAALAtex < 1000
SU_ALAALAtex: FU_ALAALAtex + BU_ALAALAtex < 1
UF_ALAALAtex: F_ALAALAtex -1000 FU_ALAALAtex < 0
UR_ALAALAtex: R_ALAALAtex -1000 BU_ALAALAtex < 0
G_ALAt2rpp: +0.59209 LC_ala-L_c -0.59209 LC_ala-L_p - DG_ALAt2rpp + DGo_ALAt2rpp = 0
FU_ALAt2rpp: DG_ALAt2rpp +1000 FU_ALAt2rpp < 1000
BU_ALAt2rpp: - DG_ALAt2rpp +1000 BU_ALAt2rpp < 1000
SU_ALAt2rpp: FU_ALAt2rpp + BU_ALAt2rpp < 1
UF_ALAt2rpp: F_ALAt2rpp -1000 FU_ALAt2rpp < 0
UR_ALAt2rpp: R_ALAt2rpp -1000 BU_ALAt2rpp < 0
G_ALAtex: -0.59209 LC_ala-L_e +0.59209 LC_ala-L_p - DG_ALAtex + DGo_ALAtex = 0
FU_ALAtex: DG_ALAtex +1000 FU_ALAtex < 1000
BU_ALAtex: - DG_ALAtex +1000 BU_ALAtex < 1000
SU_ALAtex: FU_ALAtex + BU_ALAtex < 1
UF_ALAtex: F_ALAtex -1000 FU_ALAtex < 0
UR_ALAtex: R_ALAtex -1000 BU_ALAtex < 0
G_ALCD2x: +0.59209 LC_acald_c -0.59209 LC_etoh_c -0.59209 LC_nad_c +0.59209 LC_nadh_c - DG_ALCD2x + DGo_ALCD2x = 0
FU_ALCD2x: DG_ALCD2x +1000 FU_ALCD2x < 1000
BU_ALCD2x: - DG_ALCD2x +1000 BU_ALCD2x < 1000
SU_ALCD2x: FU_ALCD2x + BU_ALCD2x < 1
UF_ALCD2x: F_ALCD2x -1000 FU_ALCD2x < 0
UR_ALCD2x: R_ALCD2x -1000 BU_ALCD2x < 0
G_ALDD2x: +0.59209 LC_ac_c -0.59209 LC_acald_c -0.59209 LC_nad_c +0.59209 LC_nadh_c - DG_ALDD2x + DGo_ALDD2x = 0
FU_ALDD2x: DG_ALDD2x +1000 FU_ALDD2x < 1000
BU_ALDD2x: - DG_ALDD2x +1000 BU_ALDD2x < 1000
SU_ALDD2x: FU_ALDD2x + BU_ALDD2x < 1
UF_ALDD2x: F_ALDD2x -1000 FU_ALDD2x < 0
UR_ALDD2x: R_ALDD2x -1000 BU_ALDD2x < 0
G_ALDD2y: +0.59209 LC_ac_c -0.59209 LC_acald_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_ALDD2y + DGo_ALDD2y = 0
FU_ALDD2y: DG_ALDD2y +1000 FU_ALDD2y < 1000
BU_ALDD2y: - DG_ALDD2y +1000 BU_ALDD2y < 1000
SU_ALDD2y: FU_ALDD2y + BU_ALDD2y < 1
UF_ALDD2y: F_ALDD2y -1000 FU_ALDD2y < 0
UR_ALDD2y: R_ALDD2y -1000 BU_ALDD2y < 0
G_ARGt3pp: -0.59209 LC_arg-L_c +0.59209 LC_arg-L_p - DG_ARGt3pp + DGo_ARGt3pp = 0
FU_ARGt3pp: DG_ARGt3pp +1000 FU_ARGt3pp < 1000
BU_ARGt3pp: - DG_ARGt3pp +1000 BU_ARGt3pp < 1000
SU_ARGt3pp: FU_ARGt3pp + BU_ARGt3pp < 1
UF_ARGt3pp: F_ARGt3pp -1000 FU_ARGt3pp < 0
UR_ARGt3pp: R_ARGt3pp -1000 BU_ARGt3pp < 0
G_ARGtex: -0.59209 LC_arg-L_e +0.59209 LC_arg-L_p - DG_ARGtex + DGo_ARGtex = 0
FU_ARGtex: DG_ARGtex +1000 FU_ARGtex < 1000
BU_ARGtex: - DG_ARGtex +1000 BU_ARGtex < 1000
SU_ARGtex: FU_ARGtex + BU_ARGtex < 1
UF_ARGtex: F_ARGtex -1000 FU_ARGtex < 0
UR_ARGtex: R_ARGtex -1000 BU_ARGtex < 0
G_ASNt2rpp: +0.59209 LC_asn-L_c -0.59209 LC_asn-L_p - DG_ASNt2rpp + DGo_ASNt2rpp = 0
FU_ASNt2rpp: DG_ASNt2rpp +1000 FU_ASNt2rpp < 1000
BU_ASNt2rpp: - DG_ASNt2rpp +1000 BU_ASNt2rpp < 1000
SU_ASNt2rpp: FU_ASNt2rpp + BU_ASNt2rpp < 1
UF_ASNt2rpp: F_ASNt2rpp -1000 FU_ASNt2rpp < 0
UR_ASNt2rpp: R_ASNt2rpp -1000 BU_ASNt2rpp < 0
G_ASNtex: -0.59209 LC_asn-L_e +0.59209 LC_asn-L_p - DG_ASNtex + DGo_ASNtex = 0
FU_ASNtex: DG_ASNtex +1000 FU_ASNtex < 1000
BU_ASNtex: - DG_ASNtex +1000 BU_ASNtex < 1000
SU_ASNtex: FU_ASNtex + BU_ASNtex < 1
UF_ASNtex: F_ASNtex -1000 FU_ASNtex < 0
UR_ASNtex: R_ASNtex -1000 BU_ASNtex < 0
G_ASPO3: -0.59209 LC_asp-L_c +0.59209 LC_iasp_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_ASPO3 + DGo_ASPO3 = 0
FU_ASPO3: DG_ASPO3 +1000 FU_ASPO3 < 1000
BU_ASPO3: - DG_ASPO3 +1000 BU_ASPO3 < 1000
SU_ASPO3: FU_ASPO3 + BU_ASPO3 < 1
UF_ASPO3: F_ASPO3 -1000 FU_ASPO3 < 0
UR_ASPO3: R_ASPO3 -1000 BU_ASPO3 < 0
G_ASPO4: -0.59209 LC_asp-L_c +0.59209 LC_iasp_c +0.59209 LC_mql8_c -0.59209 LC_mqn8_c - DG_ASPO4 + DGo_ASPO4 = 0
FU_ASPO4: DG_ASPO4 +1000 FU_ASPO4 < 1000
BU_ASPO4: - DG_ASPO4 +1000 BU_ASPO4 < 1000
SU_ASPO4: FU_ASPO4 + BU_ASPO4 < 1
UF_ASPO4: F_ASPO4 -1000 FU_ASPO4 < 0
UR_ASPO4: R_ASPO4 -1000 BU_ASPO4 < 0
G_ASPO5: -0.59209 LC_asp-L_c -0.59209 LC_fum_c +0.59209 LC_iasp_c +0.59209 LC_succ_c - DG_ASPO5 + DGo_ASPO5 = 0
FU_ASPO5: DG_ASPO5 +1000 FU_ASPO5 < 1000
BU_ASPO5: - DG_ASPO5 +1000 BU_ASPO5 < 1000
SU_ASPO5: FU_ASPO5 + BU_ASPO5 < 1
UF_ASPO5: F_ASPO5 -1000 FU_ASPO5 < 0
UR_ASPO5: R_ASPO5 -1000 BU_ASPO5 < 0
G_ASPO6: -0.59209 LC_asp-L_c +0.59209 LC_h2o2_c +0.59209 LC_iasp_c -0.59209 LC_o2_c - DG_ASPO6 + DGo_ASPO6 = 0
FU_ASPO6: DG_ASPO6 +1000 FU_ASPO6 < 1000
BU_ASPO6: - DG_ASPO6 +1000 BU_ASPO6 < 1000
SU_ASPO6: FU_ASPO6 + BU_ASPO6 < 1
UF_ASPO6: F_ASPO6 -1000 FU_ASPO6 < 0
UR_ASPO6: R_ASPO6 -1000 BU_ASPO6 < 0
G_ASPT: -0.59209 LC_asp-L_c +0.59209 LC_fum_c +0.59209 LC_nh4_c - DG_ASPT + DGo_ASPT = 0
FU_ASPT: DG_ASPT +1000 FU_ASPT < 1000
BU_ASPT: - DG_ASPT +1000 BU_ASPT < 1000
SU_ASPT: FU_ASPT + BU_ASPT < 1
UF_ASPT: F_ASPT -1000 FU_ASPT < 0
UR_ASPT: R_ASPT -1000 BU_ASPT < 0
G_ASPTA: -0.59209 LC_akg_c -0.59209 LC_asp-L_c +0.59209 LC_glu-L_c +0.59209 LC_oaa_c - DG_ASPTA + DGo_ASPTA = 0
FU_ASPTA: DG_ASPTA +1000 FU_ASPTA < 1000
BU_ASPTA: - DG_ASPTA +1000 BU_ASPTA < 1000
SU_ASPTA: FU_ASPTA + BU_ASPTA < 1
UF_ASPTA: F_ASPTA -1000 FU_ASPTA < 0
UR_ASPTA: R_ASPTA -1000 BU_ASPTA < 0
G_ASPt2rpp: +0.59209 LC_asp-L_c -0.59209 LC_asp-L_p - DG_ASPt2rpp + DGo_ASPt2rpp = 0
FU_ASPt2rpp: DG_ASPt2rpp +1000 FU_ASPt2rpp < 1000
BU_ASPt2rpp: - DG_ASPt2rpp +1000 BU_ASPt2rpp < 1000
SU_ASPt2rpp: FU_ASPt2rpp + BU_ASPt2rpp < 1
UF_ASPt2rpp: F_ASPt2rpp -1000 FU_ASPt2rpp < 0
UR_ASPt2rpp: R_ASPt2rpp -1000 BU_ASPt2rpp < 0
G_ASPtex: -0.59209 LC_asp-L_e +0.59209 LC_asp-L_p - DG_ASPtex + DGo_ASPtex = 0
FU_ASPtex: DG_ASPtex +1000 FU_ASPtex < 1000
BU_ASPtex: - DG_ASPtex +1000 BU_ASPtex < 1000
SU_ASPtex: FU_ASPtex + BU_ASPtex < 1
UF_ASPtex: F_ASPtex -1000 FU_ASPtex < 0
UR_ASPtex: R_ASPtex -1000 BU_ASPtex < 0
G_ATPM: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_pi_c - DG_ATPM + DGo_ATPM = 0
FU_ATPM: DG_ATPM +1000 FU_ATPM < 1000
BU_ATPM: - DG_ATPM +1000 BU_ATPM < 1000
SU_ATPM: FU_ATPM + BU_ATPM < 1
UF_ATPM: F_ATPM -1000 FU_ATPM < 0
UR_ATPM: R_ATPM -1000 BU_ATPM < 0
G_ATPS4rpp: -0.59209 LC_adp_c +0.59209 LC_atp_c -0.59209 LC_pi_c - DG_ATPS4rpp + DGo_ATPS4rpp = 0
FU_ATPS4rpp: DG_ATPS4rpp +1000 FU_ATPS4rpp < 1000
BU_ATPS4rpp: - DG_ATPS4rpp +1000 BU_ATPS4rpp < 1000
SU_ATPS4rpp: FU_ATPS4rpp + BU_ATPS4rpp < 1
UF_ATPS4rpp: F_ATPS4rpp -1000 FU_ATPS4rpp < 0
UR_ATPS4rpp: R_ATPS4rpp -1000 BU_ATPS4rpp < 0
SU_CA2tex: FU_CA2tex + BU_CA2tex < 1
UF_CA2tex: F_CA2tex -1000 FU_CA2tex < 0
UR_CA2tex: R_CA2tex -1000 BU_CA2tex < 0
SU_CAt6pp: FU_CAt6pp + BU_CAt6pp < 1
UF_CAt6pp: F_CAt6pp -1000 FU_CAt6pp < 0
UR_CAt6pp: R_CAt6pp -1000 BU_CAt6pp < 0
SU_CBL1tonex: FU_CBL1tonex + BU_CBL1tonex < 1
UF_CBL1tonex: F_CBL1tonex -1000 FU_CBL1tonex < 0
UR_CBL1tonex: R_CBL1tonex -1000 BU_CBL1tonex < 0
G_CITL: +0.59209 LC_ac_c -0.59209 LC_cit_c +0.59209 LC_oaa_c - DG_CITL + DGo_CITL = 0
FU_CITL: DG_CITL +1000 FU_CITL < 1000
BU_CITL: - DG_CITL +1000 BU_CITL < 1000
SU_CITL: FU_CITL + BU_CITL < 1
UF_CITL: F_CITL -1000 FU_CITL < 0
UR_CITL: R_CITL -1000 BU_CITL < 0
G_CITt3pp: -0.59209 LC_cit_c +0.59209 LC_cit_p - DG_CITt3pp + DGo_CITt3pp = 0
FU_CITt3pp: DG_CITt3pp +1000 FU_CITt3pp < 1000
BU_CITt3pp: - DG_CITt3pp +1000 BU_CITt3pp < 1000
SU_CITt3pp: FU_CITt3pp + BU_CITt3pp < 1
UF_CITt3pp: F_CITt3pp -1000 FU_CITt3pp < 0
UR_CITt3pp: R_CITt3pp -1000 BU_CITt3pp < 0
G_CITt7pp: +0.59209 LC_cit_c -0.59209 LC_succ_c -0.59209 LC_cit_p +0.59209 LC_succ_p - DG_CITt7pp + DGo_CITt7pp = 0
FU_CITt7pp: DG_CITt7pp +1000 FU_CITt7pp < 1000
BU_CITt7pp: - DG_CITt7pp +1000 BU_CITt7pp < 1000
SU_CITt7pp: FU_CITt7pp + BU_CITt7pp < 1
UF_CITt7pp: F_CITt7pp -1000 FU_CITt7pp < 0
UR_CITt7pp: R_CITt7pp -1000 BU_CITt7pp < 0
G_CITtex: -0.59209 LC_cit_e +0.59209 LC_cit_p - DG_CITtex + DGo_CITtex = 0
FU_CITtex: DG_CITtex +1000 FU_CITtex < 1000
BU_CITtex: - DG_CITtex +1000 BU_CITtex < 1000
SU_CITtex: FU_CITtex + BU_CITtex < 1
UF_CITtex: F_CITtex -1000 FU_CITtex < 0
UR_CITtex: R_CITtex -1000 BU_CITtex < 0
SU_CLt3_2pp: FU_CLt3_2pp + BU_CLt3_2pp < 1
UF_CLt3_2pp: F_CLt3_2pp -1000 FU_CLt3_2pp < 0
UR_CLt3_2pp: R_CLt3_2pp -1000 BU_CLt3_2pp < 0
SU_CLtex: FU_CLtex + BU_CLtex < 1
UF_CLtex: F_CLtex -1000 FU_CLtex < 0
UR_CLtex: R_CLtex -1000 BU_CLtex < 0
G_CO2tex: -0.59209 LC_co2_e +0.59209 LC_co2_p - DG_CO2tex + DGo_CO2tex = 0
FU_CO2tex: DG_CO2tex +1000 FU_CO2tex < 1000
BU_CO2tex: - DG_CO2tex +1000 BU_CO2tex < 1000
SU_CO2tex: FU_CO2tex + BU_CO2tex < 1
UF_CO2tex: F_CO2tex -1000 FU_CO2tex < 0
UR_CO2tex: R_CO2tex -1000 BU_CO2tex < 0
G_CO2tpp: +0.59209 LC_co2_c -0.59209 LC_co2_p - DG_CO2tpp + DGo_CO2tpp = 0
FU_CO2tpp: DG_CO2tpp +1000 FU_CO2tpp < 1000
BU_CO2tpp: - DG_CO2tpp +1000 BU_CO2tpp < 1000
SU_CO2tpp: FU_CO2tpp + BU_CO2tpp < 1
UF_CO2tpp: F_CO2tpp -1000 FU_CO2tpp < 0
UR_CO2tpp: R_CO2tpp -1000 BU_CO2tpp < 0
SU_COBALT2tex: FU_COBALT2tex + BU_COBALT2tex < 1
UF_COBALT2tex: F_COBALT2tex -1000 FU_COBALT2tex < 0
UR_COBALT2tex: R_COBALT2tex -1000 BU_COBALT2tex < 0
SU_COBALT2tpp: FU_COBALT2tpp + BU_COBALT2tpp < 1
UF_COBALT2tpp: F_COBALT2tpp -1000 FU_COBALT2tpp < 0
UR_COBALT2tpp: R_COBALT2tpp -1000 BU_COBALT2tpp < 0
G_COLIPAabcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_colipa_c +0.59209 LC_pi_c +0.59209 LC_colipa_p - DG_COLIPAabcpp + DGo_COLIPAabcpp = 0
FU_COLIPAabcpp: DG_COLIPAabcpp +1000 FU_COLIPAabcpp < 1000
BU_COLIPAabcpp: - DG_COLIPAabcpp +1000 BU_COLIPAabcpp < 1000
SU_COLIPAabcpp: FU_COLIPAabcpp + BU_COLIPAabcpp < 1
UF_COLIPAabcpp: F_COLIPAabcpp -1000 FU_COLIPAabcpp < 0
UR_COLIPAabcpp: R_COLIPAabcpp -1000 BU_COLIPAabcpp < 0
G_COLIPAabctex: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_pi_c +0.59209 LC_colipa_e -0.59209 LC_colipa_p - DG_COLIPAabctex + DGo_COLIPAabctex = 0
FU_COLIPAabctex: DG_COLIPAabctex +1000 FU_COLIPAabctex < 1000
BU_COLIPAabctex: - DG_COLIPAabctex +1000 BU_COLIPAabctex < 1000
SU_COLIPAabctex: FU_COLIPAabctex + BU_COLIPAabctex < 1
UF_COLIPAabctex: F_COLIPAabctex -1000 FU_COLIPAabctex < 0
UR_COLIPAabctex: R_COLIPAabctex -1000 BU_COLIPAabctex < 0
G_CS: -0.59209 LC_accoa_c +0.59209 LC_cit_c +0.59209 LC_coa_c -0.59209 LC_oaa_c - DG_CS + DGo_CS = 0
FU_CS: DG_CS +1000 FU_CS < 1000
BU_CS: - DG_CS +1000 BU_CS < 1000
SU_CS: FU_CS + BU_CS < 1
UF_CS: F_CS -1000 FU_CS < 0
UR_CS: R_CS -1000 BU_CS < 0
SU_CU2tex: FU_CU2tex + BU_CU2tex < 1
UF_CU2tex: F_CU2tex -1000 FU_CU2tex < 0
UR_CU2tex: R_CU2tex -1000 BU_CU2tex < 0
SU_CU2tpp: FU_CU2tpp + BU_CU2tpp < 1
UF_CU2tpp: F_CU2tpp -1000 FU_CU2tpp < 0
UR_CU2tpp: R_CU2tpp -1000 BU_CU2tpp < 0
G_CYSabc2pp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_cys-L_c +0.59209 LC_pi_c +0.59209 LC_cys-L_p - DG_CYSabc2pp + DGo_CYSabc2pp = 0
FU_CYSabc2pp: DG_CYSabc2pp +1000 FU_CYSabc2pp < 1000
BU_CYSabc2pp: - DG_CYSabc2pp +1000 BU_CYSabc2pp < 1000
SU_CYSabc2pp: FU_CYSabc2pp + BU_CYSabc2pp < 1
UF_CYSabc2pp: F_CYSabc2pp -1000 FU_CYSabc2pp < 0
UR_CYSabc2pp: R_CYSabc2pp -1000 BU_CYSabc2pp < 0
G_CYStex: -0.59209 LC_cys-L_e +0.59209 LC_cys-L_p - DG_CYStex + DGo_CYStex = 0
FU_CYStex: DG_CYStex +1000 FU_CYStex < 1000
BU_CYStex: - DG_CYStex +1000 BU_CYStex < 1000
SU_CYStex: FU_CYStex + BU_CYStex < 1
UF_CYStex: F_CYStex -1000 FU_CYStex < 0
UR_CYStex: R_CYStex -1000 BU_CYStex < 0
G_CYStpp: -0.59209 LC_cys-L_c +0.59209 LC_cys-L_p - DG_CYStpp + DGo_CYStpp = 0
FU_CYStpp: DG_CYStpp +1000 FU_CYStpp < 1000
BU_CYStpp: - DG_CYStpp +1000 BU_CYStpp < 1000
SU_CYStpp: FU_CYStpp + BU_CYStpp < 1
UF_CYStpp: F_CYStpp -1000 FU_CYStpp < 0
UR_CYStpp: R_CYStpp -1000 BU_CYStpp < 0
G_CYTBD2pp: -0.59209 LC_mql8_c +0.59209 LC_mqn8_c -0.29604 LC_o2_c - DG_CYTBD2pp + DGo_CYTBD2pp = 0
FU_CYTBD2pp: DG_CYTBD2pp +1000 FU_CYTBD2pp < 1000
BU_CYTBD2pp: - DG_CYTBD2pp +1000 BU_CYTBD2pp < 1000
SU_CYTBD2pp: FU_CYTBD2pp + BU_CYTBD2pp < 1
UF_CYTBD2pp: F_CYTBD2pp -1000 FU_CYTBD2pp < 0
UR_CYTBD2pp: R_CYTBD2pp -1000 BU_CYTBD2pp < 0
G_CYTBDpp: -0.29604 LC_o2_c +0.59209 LC_q8_c -0.59209 LC_q8h2_c - DG_CYTBDpp + DGo_CYTBDpp = 0
FU_CYTBDpp: DG_CYTBDpp +1000 FU_CYTBDpp < 1000
BU_CYTBDpp: - DG_CYTBDpp +1000 BU_CYTBDpp < 1000
SU_CYTBDpp: FU_CYTBDpp + BU_CYTBDpp < 1
UF_CYTBDpp: F_CYTBDpp -1000 FU_CYTBDpp < 0
UR_CYTBDpp: R_CYTBDpp -1000 BU_CYTBDpp < 0
G_CYTBO3_4pp: -0.29604 LC_o2_c +0.59209 LC_q8_c -0.59209 LC_q8h2_c - DG_CYTBO3_4pp + DGo_CYTBO3_4pp = 0
FU_CYTBO3_4pp: DG_CYTBO3_4pp +1000 FU_CYTBO3_4pp < 1000
BU_CYTBO3_4pp: - DG_CYTBO3_4pp +1000 BU_CYTBO3_4pp < 1000
SU_CYTBO3_4pp: FU_CYTBO3_4pp + BU_CYTBO3_4pp < 1
UF_CYTBO3_4pp: F_CYTBO3_4pp -1000 FU_CYTBO3_4pp < 0
UR_CYTBO3_4pp: R_CYTBO3_4pp -1000 BU_CYTBO3_4pp < 0
G_D-LACt2pp: +0.59209 LC_lac-D_c -0.59209 LC_lac-D_p - DG_D-LACt2pp + DGo_D-LACt2pp = 0
FU_D-LACt2pp: DG_D-LACt2pp +1000 FU_D-LACt2pp < 1000
BU_D-LACt2pp: - DG_D-LACt2pp +1000 BU_D-LACt2pp < 1000
SU_D-LACt2pp: FU_D-LACt2pp + BU_D-LACt2pp < 1
UF_D-LACt2pp: F_D-LACt2pp -1000 FU_D-LACt2pp < 0
UR_D-LACt2pp: R_D-LACt2pp -1000 BU_D-LACt2pp < 0
G_D-LACtex: -0.59209 LC_lac-D_e +0.59209 LC_lac-D_p - DG_D-LACtex + DGo_D-LACtex = 0
FU_D-LACtex: DG_D-LACtex +1000 FU_D-LACtex < 1000
BU_D-LACtex: - DG_D-LACtex +1000 BU_D-LACtex < 1000
SU_D-LACtex: FU_D-LACtex + BU_D-LACtex < 1
UF_D-LACtex: F_D-LACtex -1000 FU_D-LACtex < 0
UR_D-LACtex: R_D-LACtex -1000 BU_D-LACtex < 0
G_DALAtex: -0.59209 LC_ala-D_e +0.59209 LC_ala-D_p - DG_DALAtex + DGo_DALAtex = 0
FU_DALAtex: DG_DALAtex +1000 FU_DALAtex < 1000
BU_DALAtex: - DG_DALAtex +1000 BU_DALAtex < 1000
SU_DALAtex: FU_DALAtex + BU_DALAtex < 1
UF_DALAtex: F_DALAtex -1000 FU_DALAtex < 0
UR_DALAtex: R_DALAtex -1000 BU_DALAtex < 0
G_DB4PS: +0.59209 LC_db4p_c +0.59209 LC_for_c -0.59209 LC_ru5p-D_c - DG_DB4PS + DGo_DB4PS = 0
FU_DB4PS: DG_DB4PS +1000 FU_DB4PS < 1000
BU_DB4PS: - DG_DB4PS +1000 BU_DB4PS < 1000
SU_DB4PS: FU_DB4PS + BU_DB4PS < 1
UF_DB4PS: F_DB4PS -1000 FU_DB4PS < 0
UR_DB4PS: R_DB4PS -1000 BU_DB4PS < 0
G_DHAPT: -0.59209 LC_dha_c +0.59209 LC_dhap_c -0.59209 LC_pep_c +0.59209 LC_pyr_c - DG_DHAPT + DGo_DHAPT = 0
FU_DHAPT: DG_DHAPT +1000 FU_DHAPT < 1000
BU_DHAPT: - DG_DHAPT +1000 BU_DHAPT < 1000
SU_DHAPT: FU_DHAPT + BU_DHAPT < 1
UF_DHAPT: F_DHAPT -1000 FU_DHAPT < 0
UR_DHAPT: R_DHAPT -1000 BU_DHAPT < 0
G_DHAtex: -0.59209 LC_dha_e +0.59209 LC_dha_p - DG_DHAtex + DGo_DHAtex = 0
FU_DHAtex: DG_DHAtex +1000 FU_DHAtex < 1000
BU_DHAtex: - DG_DHAtex +1000 BU_DHAtex < 1000
SU_DHAtex: FU_DHAtex + BU_DHAtex < 1
UF_DHAtex: F_DHAtex -1000 FU_DHAtex < 0
UR_DHAtex: R_DHAtex -1000 BU_DHAtex < 0
G_DHAtpp: +0.59209 LC_dha_c -0.59209 LC_dha_p - DG_DHAtpp + DGo_DHAtpp = 0
FU_DHAtpp: DG_DHAtpp +1000 FU_DHAtpp < 1000
BU_DHAtpp: - DG_DHAtpp +1000 BU_DHAtpp < 1000
SU_DHAtpp: FU_DHAtpp + BU_DHAtpp < 1
UF_DHAtpp: F_DHAtpp -1000 FU_DHAtpp < 0
UR_DHAtpp: R_DHAtpp -1000 BU_DHAtpp < 0
G_DHORD2: -0.59209 LC_dhor-S_c +0.59209 LC_orot_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_DHORD2 + DGo_DHORD2 = 0
FU_DHORD2: DG_DHORD2 +1000 FU_DHORD2 < 1000
BU_DHORD2: - DG_DHORD2 +1000 BU_DHORD2 < 1000
SU_DHORD2: FU_DHORD2 + BU_DHORD2 < 1
UF_DHORD2: F_DHORD2 -1000 FU_DHORD2 < 0
UR_DHORD2: R_DHORD2 -1000 BU_DHORD2 < 0
G_DHORD5: -0.59209 LC_dhor-S_c +0.59209 LC_mql8_c -0.59209 LC_mqn8_c +0.59209 LC_orot_c - DG_DHORD5 + DGo_DHORD5 = 0
FU_DHORD5: DG_DHORD5 +1000 FU_DHORD5 < 1000
BU_DHORD5: - DG_DHORD5 +1000 BU_DHORD5 < 1000
SU_DHORD5: FU_DHORD5 + BU_DHORD5 < 1
UF_DHORD5: F_DHORD5 -1000 FU_DHORD5 < 0
UR_DHORD5: R_DHORD5 -1000 BU_DHORD5 < 0
G_DHORDfum: -0.59209 LC_dhor-S_c -0.59209 LC_fum_c +0.59209 LC_orot_c +0.59209 LC_succ_c - DG_DHORDfum + DGo_DHORDfum = 0
FU_DHORDfum: DG_DHORDfum +1000 FU_DHORDfum < 1000
BU_DHORDfum: - DG_DHORDfum +1000 BU_DHORDfum < 1000
SU_DHORDfum: FU_DHORDfum + BU_DHORDfum < 1
UF_DHORDfum: F_DHORDfum -1000 FU_DHORDfum < 0
UR_DHORDfum: R_DHORDfum -1000 BU_DHORDfum < 0
G_EDA: -0.59209 LC_2ddg6p_c +0.59209 LC_g3p_c +0.59209 LC_pyr_c - DG_EDA + DGo_EDA = 0
FU_EDA: DG_EDA +1000 FU_EDA < 1000
BU_EDA: - DG_EDA +1000 BU_EDA < 1000
SU_EDA: FU_EDA + BU_EDA < 1
UF_EDA: F_EDA -1000 FU_EDA < 0
UR_EDA: R_EDA -1000 BU_EDA < 0
G_EDD: +0.59209 LC_2ddg6p_c -0.59209 LC_6pgc_c - DG_EDD + DGo_EDD = 0
FU_EDD: DG_EDD +1000 FU_EDD < 1000
BU_EDD: - DG_EDD +1000 BU_EDD < 1000
SU_EDD: FU_EDD + BU_EDD < 1
UF_EDD: F_EDD -1000 FU_EDD < 0
UR_EDD: R_EDD -1000 BU_EDD < 0
G_ENO: -0.59209 LC_2pg_c +0.59209 LC_pep_c - DG_ENO + DGo_ENO = 0
FU_ENO: DG_ENO +1000 FU_ENO < 1000
BU_ENO: - DG_ENO +1000 BU_ENO < 1000
SU_ENO: FU_ENO + BU_ENO < 1
UF_ENO: F_ENO -1000 FU_ENO < 0
UR_ENO: R_ENO -1000 BU_ENO < 0
G_ETOHtex: -0.59209 LC_etoh_e +0.59209 LC_etoh_p - DG_ETOHtex + DGo_ETOHtex = 0
FU_ETOHtex: DG_ETOHtex +1000 FU_ETOHtex < 1000
BU_ETOHtex: - DG_ETOHtex +1000 BU_ETOHtex < 1000
SU_ETOHtex: FU_ETOHtex + BU_ETOHtex < 1
UF_ETOHtex: F_ETOHtex -1000 FU_ETOHtex < 0
UR_ETOHtex: R_ETOHtex -1000 BU_ETOHtex < 0
G_ETOHtrpp: +0.59209 LC_etoh_c -0.59209 LC_etoh_p - DG_ETOHtrpp + DGo_ETOHtrpp = 0
FU_ETOHtrpp: DG_ETOHtrpp +1000 FU_ETOHtrpp < 1000
BU_ETOHtrpp: - DG_ETOHtrpp +1000 BU_ETOHtrpp < 1000
SU_ETOHtrpp: FU_ETOHtrpp + BU_ETOHtrpp < 1
UF_ETOHtrpp: F_ETOHtrpp -1000 FU_ETOHtrpp < 0
UR_ETOHtrpp: R_ETOHtrpp -1000 BU_ETOHtrpp < 0
G_F6PA: +0.59209 LC_dha_c -0.59209 LC_f6p_c +0.59209 LC_g3p_c - DG_F6PA + DGo_F6PA = 0
FU_F6PA: DG_F6PA +1000 FU_F6PA < 1000
BU_F6PA: - DG_F6PA +1000 BU_F6PA < 1000
SU_F6PA: FU_F6PA + BU_F6PA < 1
UF_F6PA: F_F6PA -1000 FU_F6PA < 0
UR_F6PA: R_F6PA -1000 BU_F6PA < 0
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
G_FDH4pp: -0.59209 LC_q8_c +0.59209 LC_q8h2_c +0.59209 LC_co2_p -0.59209 LC_for_p - DG_FDH4pp + DGo_FDH4pp = 0
FU_FDH4pp: DG_FDH4pp +1000 FU_FDH4pp < 1000
BU_FDH4pp: - DG_FDH4pp +1000 BU_FDH4pp < 1000
SU_FDH4pp: FU_FDH4pp + BU_FDH4pp < 1
UF_FDH4pp: F_FDH4pp -1000 FU_FDH4pp < 0
UR_FDH4pp: R_FDH4pp -1000 BU_FDH4pp < 0
G_FDH5pp: +0.59209 LC_mql8_c -0.59209 LC_mqn8_c +0.59209 LC_co2_p -0.59209 LC_for_p - DG_FDH5pp + DGo_FDH5pp = 0
FU_FDH5pp: DG_FDH5pp +1000 FU_FDH5pp < 1000
BU_FDH5pp: - DG_FDH5pp +1000 BU_FDH5pp < 1000
SU_FDH5pp: FU_FDH5pp + BU_FDH5pp < 1
UF_FDH5pp: F_FDH5pp -1000 FU_FDH5pp < 0
UR_FDH5pp: R_FDH5pp -1000 BU_FDH5pp < 0
G_FE2abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_fe2_c +0.59209 LC_pi_c -0.59209 LC_fe2_p - DG_FE2abcpp + DGo_FE2abcpp = 0
FU_FE2abcpp: DG_FE2abcpp +1000 FU_FE2abcpp < 1000
BU_FE2abcpp: - DG_FE2abcpp +1000 BU_FE2abcpp < 1000
SU_FE2abcpp: FU_FE2abcpp + BU_FE2abcpp < 1
UF_FE2abcpp: F_FE2abcpp -1000 FU_FE2abcpp < 0
UR_FE2abcpp: R_FE2abcpp -1000 BU_FE2abcpp < 0
G_FE2t2pp: +0.59209 LC_fe2_c -0.59209 LC_fe2_p - DG_FE2t2pp + DGo_FE2t2pp = 0
FU_FE2t2pp: DG_FE2t2pp +1000 FU_FE2t2pp < 1000
BU_FE2t2pp: - DG_FE2t2pp +1000 BU_FE2t2pp < 1000
SU_FE2t2pp: FU_FE2t2pp + BU_FE2t2pp < 1
UF_FE2t2pp: F_FE2t2pp -1000 FU_FE2t2pp < 0
UR_FE2t2pp: R_FE2t2pp -1000 BU_FE2t2pp < 0
G_FE2tex: -0.59209 LC_fe2_e +0.59209 LC_fe2_p - DG_FE2tex + DGo_FE2tex = 0
FU_FE2tex: DG_FE2tex +1000 FU_FE2tex < 1000
BU_FE2tex: - DG_FE2tex +1000 BU_FE2tex < 1000
SU_FE2tex: FU_FE2tex + BU_FE2tex < 1
UF_FE2tex: F_FE2tex -1000 FU_FE2tex < 0
UR_FE2tex: R_FE2tex -1000 BU_FE2tex < 0
G_FE2tpp: +0.59209 LC_fe2_c -0.59209 LC_fe2_p - DG_FE2tpp + DGo_FE2tpp = 0
FU_FE2tpp: DG_FE2tpp +1000 FU_FE2tpp < 1000
BU_FE2tpp: - DG_FE2tpp +1000 BU_FE2tpp < 1000
SU_FE2tpp: FU_FE2tpp + BU_FE2tpp < 1
UF_FE2tpp: F_FE2tpp -1000 FU_FE2tpp < 0
UR_FE2tpp: R_FE2tpp -1000 BU_FE2tpp < 0
G_FE3abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_fe3_c +0.59209 LC_pi_c -0.59209 LC_fe3_p - DG_FE3abcpp + DGo_FE3abcpp = 0
FU_FE3abcpp: DG_FE3abcpp +1000 FU_FE3abcpp < 1000
BU_FE3abcpp: - DG_FE3abcpp +1000 BU_FE3abcpp < 1000
SU_FE3abcpp: FU_FE3abcpp + BU_FE3abcpp < 1
UF_FE3abcpp: F_FE3abcpp -1000 FU_FE3abcpp < 0
UR_FE3abcpp: R_FE3abcpp -1000 BU_FE3abcpp < 0
G_FE3tex: -0.59209 LC_fe3_e +0.59209 LC_fe3_p - DG_FE3tex + DGo_FE3tex = 0
FU_FE3tex: DG_FE3tex +1000 FU_FE3tex < 1000
BU_FE3tex: - DG_FE3tex +1000 BU_FE3tex < 1000
SU_FE3tex: FU_FE3tex + BU_FE3tex < 1
UF_FE3tex: F_FE3tex -1000 FU_FE3tex < 0
UR_FE3tex: R_FE3tex -1000 BU_FE3tex < 0
SU_FEENTERtex: FU_FEENTERtex + BU_FEENTERtex < 1
UF_FEENTERtex: F_FEENTERtex -1000 FU_FEENTERtex < 0
UR_FEENTERtex: R_FEENTERtex -1000 BU_FEENTERtex < 0
SU_FEENTERtpp: FU_FEENTERtpp + BU_FEENTERtpp < 1
UF_FEENTERtpp: F_FEENTERtpp -1000 FU_FEENTERtpp < 0
UR_FEENTERtpp: R_FEENTERtpp -1000 BU_FEENTERtpp < 0
SU_FLDR2: FU_FLDR2 + BU_FLDR2 < 1
UF_FLDR2: F_FLDR2 -1000 FU_FLDR2 < 0
UR_FLDR2: R_FLDR2 -1000 BU_FLDR2 < 0
G_FORtex: -0.59209 LC_for_e +0.59209 LC_for_p - DG_FORtex + DGo_FORtex = 0
FU_FORtex: DG_FORtex +1000 FU_FORtex < 1000
BU_FORtex: - DG_FORtex +1000 BU_FORtex < 1000
SU_FORtex: FU_FORtex + BU_FORtex < 1
UF_FORtex: F_FORtex -1000 FU_FORtex < 0
UR_FORtex: R_FORtex -1000 BU_FORtex < 0
G_FORtppi: -0.59209 LC_for_c +0.59209 LC_for_p - DG_FORtppi + DGo_FORtppi = 0
FU_FORtppi: DG_FORtppi +1000 FU_FORtppi < 1000
BU_FORtppi: - DG_FORtppi +1000 BU_FORtppi < 1000
SU_FORtppi: FU_FORtppi + BU_FORtppi < 1
UF_FORtppi: F_FORtppi -1000 FU_FORtppi < 0
UR_FORtppi: R_FORtppi -1000 BU_FORtppi < 0
G_FRD3: +0.59209 LC_2dmmq8_c -0.59209 LC_2dmmql8_c -0.59209 LC_fum_c +0.59209 LC_succ_c - DG_FRD3 + DGo_FRD3 = 0
FU_FRD3: DG_FRD3 +1000 FU_FRD3 < 1000
BU_FRD3: - DG_FRD3 +1000 BU_FRD3 < 1000
SU_FRD3: FU_FRD3 + BU_FRD3 < 1
UF_FRD3: F_FRD3 -1000 FU_FRD3 < 0
UR_FRD3: R_FRD3 -1000 BU_FRD3 < 0
G_FUM: -0.59209 LC_fum_c +0.59209 LC_mal-L_c - DG_FUM + DGo_FUM = 0
FU_FUM: DG_FUM +1000 FU_FUM < 1000
BU_FUM: - DG_FUM +1000 BU_FUM < 1000
SU_FUM: FU_FUM + BU_FUM < 1
UF_FUM: F_FUM -1000 FU_FUM < 0
UR_FUM: R_FUM -1000 BU_FUM < 0
G_G3PD2: +0.59209 LC_dhap_c -0.59209 LC_glyc3p_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_G3PD2 + DGo_G3PD2 = 0
FU_G3PD2: DG_G3PD2 +1000 FU_G3PD2 < 1000
BU_G3PD2: - DG_G3PD2 +1000 BU_G3PD2 < 1000
SU_G3PD2: FU_G3PD2 + BU_G3PD2 < 1
UF_G3PD2: F_G3PD2 -1000 FU_G3PD2 < 0
UR_G3PD2: R_G3PD2 -1000 BU_G3PD2 < 0
G_G3PD5: +0.59209 LC_dhap_c -0.59209 LC_glyc3p_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_G3PD5 + DGo_G3PD5 = 0
FU_G3PD5: DG_G3PD5 +1000 FU_G3PD5 < 1000
BU_G3PD5: - DG_G3PD5 +1000 BU_G3PD5 < 1000
SU_G3PD5: FU_G3PD5 + BU_G3PD5 < 1
UF_G3PD5: F_G3PD5 -1000 FU_G3PD5 < 0
UR_G3PD5: R_G3PD5 -1000 BU_G3PD5 < 0
G_G3PD6: +0.59209 LC_dhap_c -0.59209 LC_glyc3p_c +0.59209 LC_mql8_c -0.59209 LC_mqn8_c - DG_G3PD6 + DGo_G3PD6 = 0
FU_G3PD6: DG_G3PD6 +1000 FU_G3PD6 < 1000
BU_G3PD6: - DG_G3PD6 +1000 BU_G3PD6 < 1000
SU_G3PD6: FU_G3PD6 + BU_G3PD6 < 1
UF_G3PD6: F_G3PD6 -1000 FU_G3PD6 < 0
UR_G3PD6: R_G3PD6 -1000 BU_G3PD6 < 0
G_G3PD7: -0.59209 LC_2dmmq8_c +0.59209 LC_2dmmql8_c +0.59209 LC_dhap_c -0.59209 LC_glyc3p_c - DG_G3PD7 + DGo_G3PD7 = 0
FU_G3PD7: DG_G3PD7 +1000 FU_G3PD7 < 1000
BU_G3PD7: - DG_G3PD7 +1000 BU_G3PD7 < 1000
SU_G3PD7: FU_G3PD7 + BU_G3PD7 < 1
UF_G3PD7: F_G3PD7 -1000 FU_G3PD7 < 0
UR_G3PD7: R_G3PD7 -1000 BU_G3PD7 < 0
G_G3PT: +0.59209 LC_glyc_c -0.59209 LC_glyc3p_c +0.59209 LC_pi_c - DG_G3PT + DGo_G3PT = 0
FU_G3PT: DG_G3PT +1000 FU_G3PT < 1000
BU_G3PT: - DG_G3PT +1000 BU_G3PT < 1000
SU_G3PT: FU_G3PT + BU_G3PT < 1
UF_G3PT: F_G3PT -1000 FU_G3PT < 0
UR_G3PT: R_G3PT -1000 BU_G3PT < 0
G_G6PDH2r: +0.59209 LC_6pgl_c -0.59209 LC_g6p_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_G6PDH2r + DGo_G6PDH2r = 0
FU_G6PDH2r: DG_G6PDH2r +1000 FU_G6PDH2r < 1000
BU_G6PDH2r: - DG_G6PDH2r +1000 BU_G6PDH2r < 1000
SU_G6PDH2r: FU_G6PDH2r + BU_G6PDH2r < 1
UF_G6PDH2r: F_G6PDH2r -1000 FU_G6PDH2r < 0
UR_G6PDH2r: R_G6PDH2r -1000 BU_G6PDH2r < 0
G_G6PP: -0.59209 LC_g6p_c +0.59209 LC_glc-D_c +0.59209 LC_pi_c - DG_G6PP + DGo_G6PP = 0
FU_G6PP: DG_G6PP +1000 FU_G6PP < 1000
BU_G6PP: - DG_G6PP +1000 BU_G6PP < 1000
SU_G6PP: FU_G6PP + BU_G6PP < 1
UF_G6PP: F_G6PP -1000 FU_G6PP < 0
UR_G6PP: R_G6PP -1000 BU_G6PP < 0
G_GAPD: +0.59209 LC_13dpg_c -0.59209 LC_g3p_c -0.59209 LC_nad_c +0.59209 LC_nadh_c -0.59209 LC_pi_c - DG_GAPD + DGo_GAPD = 0
FU_GAPD: DG_GAPD +1000 FU_GAPD < 1000
BU_GAPD: - DG_GAPD +1000 BU_GAPD < 1000
SU_GAPD: FU_GAPD + BU_GAPD < 1
UF_GAPD: F_GAPD -1000 FU_GAPD < 0
UR_GAPD: R_GAPD -1000 BU_GAPD < 0
G_GLBRAN2: +0.59209 LC_bglycogen_c -0.59209 LC_glycogen_c - DG_GLBRAN2 + DGo_GLBRAN2 = 0
FU_GLBRAN2: DG_GLBRAN2 +1000 FU_GLBRAN2 < 1000
BU_GLBRAN2: - DG_GLBRAN2 +1000 BU_GLBRAN2 < 1000
SU_GLBRAN2: FU_GLBRAN2 + BU_GLBRAN2 < 1
UF_GLBRAN2: F_GLBRAN2 -1000 FU_GLBRAN2 < 0
UR_GLBRAN2: R_GLBRAN2 -1000 BU_GLBRAN2 < 0
G_GLCDpp: -0.59209 LC_q8_c +0.59209 LC_q8h2_c -0.59209 LC_glc-D_p +0.59209 LC_glcn_p - DG_GLCDpp + DGo_GLCDpp = 0
FU_GLCDpp: DG_GLCDpp +1000 FU_GLCDpp < 1000
BU_GLCDpp: - DG_GLCDpp +1000 BU_GLCDpp < 1000
SU_GLCDpp: FU_GLCDpp + BU_GLCDpp < 1
UF_GLCDpp: F_GLCDpp -1000 FU_GLCDpp < 0
UR_GLCDpp: R_GLCDpp -1000 BU_GLCDpp < 0
G_GLCNtex: -0.59209 LC_glcn_e +0.59209 LC_glcn_p - DG_GLCNtex + DGo_GLCNtex = 0
FU_GLCNtex: DG_GLCNtex +1000 FU_GLCNtex < 1000
BU_GLCNtex: - DG_GLCNtex +1000 BU_GLCNtex < 1000
SU_GLCNtex: FU_GLCNtex + BU_GLCNtex < 1
UF_GLCNtex: F_GLCNtex -1000 FU_GLCNtex < 0
UR_GLCNtex: R_GLCNtex -1000 BU_GLCNtex < 0
G_GLCP: +0.59209 LC_g1p_c -0.59209 LC_glycogen_c +0.59209 LC_glycogenn1_c -0.59209 LC_pi_c - DG_GLCP + DGo_GLCP = 0
FU_GLCP: DG_GLCP +1000 FU_GLCP < 1000
BU_GLCP: - DG_GLCP +1000 BU_GLCP < 1000
SU_GLCP: FU_GLCP + BU_GLCP < 1
UF_GLCP: F_GLCP -1000 FU_GLCP < 0
UR_GLCP: R_GLCP -1000 BU_GLCP < 0
G_GLCS1: +0.59209 LC_adp_c -0.59209 LC_adpglc_c +0.59209 LC_glycogen_c -0.59209 LC_glycogenn1_c - DG_GLCS1 + DGo_GLCS1 = 0
FU_GLCS1: DG_GLCS1 +1000 FU_GLCS1 < 1000
BU_GLCS1: - DG_GLCS1 +1000 BU_GLCS1 < 1000
SU_GLCS1: FU_GLCS1 + BU_GLCS1 < 1
UF_GLCS1: F_GLCS1 -1000 FU_GLCS1 < 0
UR_GLCS1: R_GLCS1 -1000 BU_GLCS1 < 0
G_GLCabcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_glc-D_c +0.59209 LC_pi_c -0.59209 LC_glc-D_p - DG_GLCabcpp + DGo_GLCabcpp = 0
FU_GLCabcpp: DG_GLCabcpp +1000 FU_GLCabcpp < 1000
BU_GLCabcpp: - DG_GLCabcpp +1000 BU_GLCabcpp < 1000
SU_GLCabcpp: FU_GLCabcpp + BU_GLCabcpp < 1
UF_GLCabcpp: F_GLCabcpp -1000 FU_GLCabcpp < 0
UR_GLCabcpp: R_GLCabcpp -1000 BU_GLCabcpp < 0
G_GLCptspp: +0.59209 LC_g6p_c -0.59209 LC_pep_c +0.59209 LC_pyr_c -0.59209 LC_glc-D_p - DG_GLCptspp + DGo_GLCptspp = 0
FU_GLCptspp: DG_GLCptspp +1000 FU_GLCptspp < 1000
BU_GLCptspp: - DG_GLCptspp +1000 BU_GLCptspp < 1000
SU_GLCptspp: FU_GLCptspp + BU_GLCptspp < 1
UF_GLCptspp: F_GLCptspp -1000 FU_GLCptspp < 0
UR_GLCptspp: R_GLCptspp -1000 BU_GLCptspp < 0
G_GLCt2pp: +0.59209 LC_glc-D_c -0.59209 LC_glc-D_p - DG_GLCt2pp + DGo_GLCt2pp = 0
FU_GLCt2pp: DG_GLCt2pp +1000 FU_GLCt2pp < 1000
BU_GLCt2pp: - DG_GLCt2pp +1000 BU_GLCt2pp < 1000
SU_GLCt2pp: FU_GLCt2pp + BU_GLCt2pp < 1
UF_GLCt2pp: F_GLCt2pp -1000 FU_GLCt2pp < 0
UR_GLCt2pp: R_GLCt2pp -1000 BU_GLCt2pp < 0
G_GLCtex: -0.59209 LC_glc-D_e +0.59209 LC_glc-D_p - DG_GLCtex + DGo_GLCtex = 0
FU_GLCtex: DG_GLCtex +1000 FU_GLCtex < 1000
BU_GLCtex: - DG_GLCtex +1000 BU_GLCtex < 1000
SU_GLCtex: FU_GLCtex + BU_GLCtex < 1
UF_GLCtex: F_GLCtex -1000 FU_GLCtex < 0
UR_GLCtex: R_GLCtex -1000 BU_GLCtex < 0
G_GLCtexi: -0.59209 LC_glc-D_e +0.59209 LC_glc-D_p - DG_GLCtexi + DGo_GLCtexi = 0
FU_GLCtexi: DG_GLCtexi +1000 FU_GLCtexi < 1000
BU_GLCtexi: - DG_GLCtexi +1000 BU_GLCtexi < 1000
SU_GLCtexi: FU_GLCtexi + BU_GLCtexi < 1
UF_GLCtexi: F_GLCtexi -1000 FU_GLCtexi < 0
UR_GLCtexi: R_GLCtexi -1000 BU_GLCtexi < 0
G_GLDBRAN2: -0.59209 LC_bglycogen_c +0.59209 LC_glycogen_c - DG_GLDBRAN2 + DGo_GLDBRAN2 = 0
FU_GLDBRAN2: DG_GLDBRAN2 +1000 FU_GLDBRAN2 < 1000
BU_GLDBRAN2: - DG_GLDBRAN2 +1000 BU_GLDBRAN2 < 1000
SU_GLDBRAN2: FU_GLDBRAN2 + BU_GLDBRAN2 < 1
UF_GLDBRAN2: F_GLDBRAN2 -1000 FU_GLDBRAN2 < 0
UR_GLDBRAN2: R_GLDBRAN2 -1000 BU_GLDBRAN2 < 0
G_GLGC: +0.59209 LC_adpglc_c -0.59209 LC_atp_c -0.59209 LC_g1p_c +0.59209 LC_ppi_c - DG_GLGC + DGo_GLGC = 0
FU_GLGC: DG_GLGC +1000 FU_GLGC < 1000
BU_GLGC: - DG_GLGC +1000 BU_GLGC < 1000
SU_GLGC: FU_GLGC + BU_GLGC < 1
UF_GLGC: F_GLGC -1000 FU_GLGC < 0
UR_GLGC: R_GLGC -1000 BU_GLGC < 0
G_GLUDy: +0.59209 LC_akg_c -0.59209 LC_glu-L_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c +0.59209 LC_nh4_c - DG_GLUDy + DGo_GLUDy = 0
FU_GLUDy: DG_GLUDy +1000 FU_GLUDy < 1000
BU_GLUDy: - DG_GLUDy +1000 BU_GLUDy < 1000
SU_GLUDy: FU_GLUDy + BU_GLUDy < 1
UF_GLUDy: F_GLUDy -1000 FU_GLUDy < 0
UR_GLUDy: R_GLUDy -1000 BU_GLUDy < 0
G_GLUabcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_glu-L_c +0.59209 LC_pi_c -0.59209 LC_glu-L_p - DG_GLUabcpp + DGo_GLUabcpp = 0
FU_GLUabcpp: DG_GLUabcpp +1000 FU_GLUabcpp < 1000
BU_GLUabcpp: - DG_GLUabcpp +1000 BU_GLUabcpp < 1000
SU_GLUabcpp: FU_GLUabcpp + BU_GLUabcpp < 1
UF_GLUabcpp: F_GLUabcpp -1000 FU_GLUabcpp < 0
UR_GLUabcpp: R_GLUabcpp -1000 BU_GLUabcpp < 0
G_GLUt2rpp: +0.59209 LC_glu-L_c -0.59209 LC_glu-L_p - DG_GLUt2rpp + DGo_GLUt2rpp = 0
FU_GLUt2rpp: DG_GLUt2rpp +1000 FU_GLUt2rpp < 1000
BU_GLUt2rpp: - DG_GLUt2rpp +1000 BU_GLUt2rpp < 1000
SU_GLUt2rpp: FU_GLUt2rpp + BU_GLUt2rpp < 1
UF_GLUt2rpp: F_GLUt2rpp -1000 FU_GLUt2rpp < 0
UR_GLUt2rpp: R_GLUt2rpp -1000 BU_GLUt2rpp < 0
SU_GLUt4pp: FU_GLUt4pp + BU_GLUt4pp < 1
UF_GLUt4pp: F_GLUt4pp -1000 FU_GLUt4pp < 0
UR_GLUt4pp: R_GLUt4pp -1000 BU_GLUt4pp < 0
G_GLUtex: -0.59209 LC_glu-L_e +0.59209 LC_glu-L_p - DG_GLUtex + DGo_GLUtex = 0
FU_GLUtex: DG_GLUtex +1000 FU_GLUtex < 1000
BU_GLUtex: - DG_GLUtex +1000 BU_GLUtex < 1000
SU_GLUtex: FU_GLUtex + BU_GLUtex < 1
UF_GLUtex: F_GLUtex -1000 FU_GLUtex < 0
UR_GLUtex: R_GLUtex -1000 BU_GLUtex < 0
G_GLXCL: +0.59209 LC_2h3oppan_c +0.59209 LC_co2_c -1.1842 LC_glx_c - DG_GLXCL + DGo_GLXCL = 0
FU_GLXCL: DG_GLXCL +1000 FU_GLXCL < 1000
BU_GLXCL: - DG_GLXCL +1000 BU_GLXCL < 1000
SU_GLXCL: FU_GLXCL + BU_GLXCL < 1
UF_GLXCL: F_GLXCL -1000 FU_GLXCL < 0
UR_GLXCL: R_GLXCL -1000 BU_GLXCL < 0
G_GLYCAt2rpp: +0.59209 LC_glyc-R_c -0.59209 LC_glyc-R_p - DG_GLYCAt2rpp + DGo_GLYCAt2rpp = 0
FU_GLYCAt2rpp: DG_GLYCAt2rpp +1000 FU_GLYCAt2rpp < 1000
BU_GLYCAt2rpp: - DG_GLYCAt2rpp +1000 BU_GLYCAt2rpp < 1000
SU_GLYCAt2rpp: FU_GLYCAt2rpp + BU_GLYCAt2rpp < 1
UF_GLYCAt2rpp: F_GLYCAt2rpp -1000 FU_GLYCAt2rpp < 0
UR_GLYCAt2rpp: R_GLYCAt2rpp -1000 BU_GLYCAt2rpp < 0
G_GLYCAtex: -0.59209 LC_glyc-R_e +0.59209 LC_glyc-R_p - DG_GLYCAtex + DGo_GLYCAtex = 0
FU_GLYCAtex: DG_GLYCAtex +1000 FU_GLYCAtex < 1000
BU_GLYCAtex: - DG_GLYCAtex +1000 BU_GLYCAtex < 1000
SU_GLYCAtex: FU_GLYCAtex + BU_GLYCAtex < 1
UF_GLYCAtex: F_GLYCAtex -1000 FU_GLYCAtex < 0
UR_GLYCAtex: R_GLYCAtex -1000 BU_GLYCAtex < 0
G_GLYCDx: +0.59209 LC_dha_c -0.59209 LC_glyc_c -0.59209 LC_nad_c +0.59209 LC_nadh_c - DG_GLYCDx + DGo_GLYCDx = 0
FU_GLYCDx: DG_GLYCDx +1000 FU_GLYCDx < 1000
BU_GLYCDx: - DG_GLYCDx +1000 BU_GLYCDx < 1000
SU_GLYCDx: FU_GLYCDx + BU_GLYCDx < 1
UF_GLYCDx: F_GLYCDx -1000 FU_GLYCDx < 0
UR_GLYCDx: R_GLYCDx -1000 BU_GLYCDx < 0
G_GLYCK: +0.59209 LC_3pg_c +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_glyc-R_c - DG_GLYCK + DGo_GLYCK = 0
FU_GLYCK: DG_GLYCK +1000 FU_GLYCK < 1000
BU_GLYCK: - DG_GLYCK +1000 BU_GLYCK < 1000
SU_GLYCK: FU_GLYCK + BU_GLYCK < 1
UF_GLYCK: F_GLYCK -1000 FU_GLYCK < 0
UR_GLYCK: R_GLYCK -1000 BU_GLYCK < 0
G_GLYCK2: +0.59209 LC_2pg_c +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_glyc-R_c - DG_GLYCK2 + DGo_GLYCK2 = 0
FU_GLYCK2: DG_GLYCK2 +1000 FU_GLYCK2 < 1000
BU_GLYCK2: - DG_GLYCK2 +1000 BU_GLYCK2 < 1000
SU_GLYCK2: FU_GLYCK2 + BU_GLYCK2 < 1
UF_GLYCK2: F_GLYCK2 -1000 FU_GLYCK2 < 0
UR_GLYCK2: R_GLYCK2 -1000 BU_GLYCK2 < 0
G_GLYCLTDx: -0.59209 LC_glx_c +0.59209 LC_glyclt_c +0.59209 LC_nad_c -0.59209 LC_nadh_c - DG_GLYCLTDx + DGo_GLYCLTDx = 0
FU_GLYCLTDx: DG_GLYCLTDx +1000 FU_GLYCLTDx < 1000
BU_GLYCLTDx: - DG_GLYCLTDx +1000 BU_GLYCLTDx < 1000
SU_GLYCLTDx: FU_GLYCLTDx + BU_GLYCLTDx < 1
UF_GLYCLTDx: F_GLYCLTDx -1000 FU_GLYCLTDx < 0
UR_GLYCLTDx: R_GLYCLTDx -1000 BU_GLYCLTDx < 0
G_GLYCLTDy: -0.59209 LC_glx_c +0.59209 LC_glyclt_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c - DG_GLYCLTDy + DGo_GLYCLTDy = 0
FU_GLYCLTDy: DG_GLYCLTDy +1000 FU_GLYCLTDy < 1000
BU_GLYCLTDy: - DG_GLYCLTDy +1000 BU_GLYCLTDy < 1000
SU_GLYCLTDy: FU_GLYCLTDy + BU_GLYCLTDy < 1
UF_GLYCLTDy: F_GLYCLTDy -1000 FU_GLYCLTDy < 0
UR_GLYCLTDy: R_GLYCLTDy -1000 BU_GLYCLTDy < 0
G_GLYCLTt2rpp: +0.59209 LC_glyclt_c -0.59209 LC_glyclt_p - DG_GLYCLTt2rpp + DGo_GLYCLTt2rpp = 0
FU_GLYCLTt2rpp: DG_GLYCLTt2rpp +1000 FU_GLYCLTt2rpp < 1000
BU_GLYCLTt2rpp: - DG_GLYCLTt2rpp +1000 BU_GLYCLTt2rpp < 1000
SU_GLYCLTt2rpp: FU_GLYCLTt2rpp + BU_GLYCLTt2rpp < 1
UF_GLYCLTt2rpp: F_GLYCLTt2rpp -1000 FU_GLYCLTt2rpp < 0
UR_GLYCLTt2rpp: R_GLYCLTt2rpp -1000 BU_GLYCLTt2rpp < 0
G_GLYCLTtex: -0.59209 LC_glyclt_e +0.59209 LC_glyclt_p - DG_GLYCLTtex + DGo_GLYCLTtex = 0
FU_GLYCLTtex: DG_GLYCLTtex +1000 FU_GLYCLTtex < 1000
BU_GLYCLTtex: - DG_GLYCLTtex +1000 BU_GLYCLTtex < 1000
SU_GLYCLTtex: FU_GLYCLTtex + BU_GLYCLTtex < 1
UF_GLYCLTtex: F_GLYCLTtex -1000 FU_GLYCLTtex < 0
UR_GLYCLTtex: R_GLYCLTtex -1000 BU_GLYCLTtex < 0
G_GLYCTO2: +0.59209 LC_glx_c -0.59209 LC_glyclt_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_GLYCTO2 + DGo_GLYCTO2 = 0
FU_GLYCTO2: DG_GLYCTO2 +1000 FU_GLYCTO2 < 1000
BU_GLYCTO2: - DG_GLYCTO2 +1000 BU_GLYCTO2 < 1000
SU_GLYCTO2: FU_GLYCTO2 + BU_GLYCTO2 < 1
UF_GLYCTO2: F_GLYCTO2 -1000 FU_GLYCTO2 < 0
UR_GLYCTO2: R_GLYCTO2 -1000 BU_GLYCTO2 < 0
G_GLYCTO3: +0.59209 LC_glx_c -0.59209 LC_glyclt_c +0.59209 LC_mql8_c -0.59209 LC_mqn8_c - DG_GLYCTO3 + DGo_GLYCTO3 = 0
FU_GLYCTO3: DG_GLYCTO3 +1000 FU_GLYCTO3 < 1000
BU_GLYCTO3: - DG_GLYCTO3 +1000 BU_GLYCTO3 < 1000
SU_GLYCTO3: FU_GLYCTO3 + BU_GLYCTO3 < 1
UF_GLYCTO3: F_GLYCTO3 -1000 FU_GLYCTO3 < 0
UR_GLYCTO3: R_GLYCTO3 -1000 BU_GLYCTO3 < 0
G_GLYCTO4: -0.59209 LC_2dmmq8_c +0.59209 LC_2dmmql8_c +0.59209 LC_glx_c -0.59209 LC_glyclt_c - DG_GLYCTO4 + DGo_GLYCTO4 = 0
FU_GLYCTO4: DG_GLYCTO4 +1000 FU_GLYCTO4 < 1000
BU_GLYCTO4: - DG_GLYCTO4 +1000 BU_GLYCTO4 < 1000
SU_GLYCTO4: FU_GLYCTO4 + BU_GLYCTO4 < 1
UF_GLYCTO4: F_GLYCTO4 -1000 FU_GLYCTO4 < 0
UR_GLYCTO4: R_GLYCTO4 -1000 BU_GLYCTO4 < 0
G_GLYCtex: -0.59209 LC_glyc_e +0.59209 LC_glyc_p - DG_GLYCtex + DGo_GLYCtex = 0
FU_GLYCtex: DG_GLYCtex +1000 FU_GLYCtex < 1000
BU_GLYCtex: - DG_GLYCtex +1000 BU_GLYCtex < 1000
SU_GLYCtex: FU_GLYCtex + BU_GLYCtex < 1
UF_GLYCtex: F_GLYCtex -1000 FU_GLYCtex < 0
UR_GLYCtex: R_GLYCtex -1000 BU_GLYCtex < 0
G_GLYCtpp: -0.59209 LC_glyc_c +0.59209 LC_glyc_p - DG_GLYCtpp + DGo_GLYCtpp = 0
FU_GLYCtpp: DG_GLYCtpp +1000 FU_GLYCtpp < 1000
BU_GLYCtpp: - DG_GLYCtpp +1000 BU_GLYCtpp < 1000
SU_GLYCtpp: FU_GLYCtpp + BU_GLYCtpp < 1
UF_GLYCtpp: F_GLYCtpp -1000 FU_GLYCtpp < 0
UR_GLYCtpp: R_GLYCtpp -1000 BU_GLYCtpp < 0
G_GLYK: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_glyc_c +0.59209 LC_glyc3p_c - DG_GLYK + DGo_GLYK = 0
FU_GLYK: DG_GLYK +1000 FU_GLYK < 1000
BU_GLYK: - DG_GLYK +1000 BU_GLYK < 1000
SU_GLYK: FU_GLYK + BU_GLYK < 1
UF_GLYK: F_GLYK -1000 FU_GLYK < 0
UR_GLYK: R_GLYK -1000 BU_GLYK < 0
G_GLYt2rpp: +0.59209 LC_gly_c -0.59209 LC_gly_p - DG_GLYt2rpp + DGo_GLYt2rpp = 0
FU_GLYt2rpp: DG_GLYt2rpp +1000 FU_GLYt2rpp < 1000
BU_GLYt2rpp: - DG_GLYt2rpp +1000 BU_GLYt2rpp < 1000
SU_GLYt2rpp: FU_GLYt2rpp + BU_GLYt2rpp < 1
UF_GLYt2rpp: F_GLYt2rpp -1000 FU_GLYt2rpp < 0
UR_GLYt2rpp: R_GLYt2rpp -1000 BU_GLYt2rpp < 0
G_GLYtex: -0.59209 LC_gly_e +0.59209 LC_gly_p - DG_GLYtex + DGo_GLYtex = 0
FU_GLYtex: DG_GLYtex +1000 FU_GLYtex < 1000
BU_GLYtex: - DG_GLYtex +1000 BU_GLYtex < 1000
SU_GLYtex: FU_GLYtex + BU_GLYtex < 1
UF_GLYtex: F_GLYtex -1000 FU_GLYtex < 0
UR_GLYtex: R_GLYtex -1000 BU_GLYtex < 0
G_GND: -0.59209 LC_6pgc_c +0.59209 LC_co2_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c +0.59209 LC_ru5p-D_c - DG_GND + DGo_GND = 0
FU_GND: DG_GND +1000 FU_GND < 1000
BU_GND: - DG_GND +1000 BU_GND < 1000
SU_GND: FU_GND + BU_GND < 1
UF_GND: F_GND -1000 FU_GND < 0
UR_GND: R_GND -1000 BU_GND < 0
G_GTHRDabc2pp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_gthrd_c +0.59209 LC_pi_c +0.59209 LC_gthrd_p - DG_GTHRDabc2pp + DGo_GTHRDabc2pp = 0
FU_GTHRDabc2pp: DG_GTHRDabc2pp +1000 FU_GTHRDabc2pp < 1000
BU_GTHRDabc2pp: - DG_GTHRDabc2pp +1000 BU_GTHRDabc2pp < 1000
SU_GTHRDabc2pp: FU_GTHRDabc2pp + BU_GTHRDabc2pp < 1
UF_GTHRDabc2pp: F_GTHRDabc2pp -1000 FU_GTHRDabc2pp < 0
UR_GTHRDabc2pp: R_GTHRDabc2pp -1000 BU_GTHRDabc2pp < 0
G_GTHRDtex: -0.59209 LC_gthrd_e +0.59209 LC_gthrd_p - DG_GTHRDtex + DGo_GTHRDtex = 0
FU_GTHRDtex: DG_GTHRDtex +1000 FU_GTHRDtex < 1000
BU_GTHRDtex: - DG_GTHRDtex +1000 BU_GTHRDtex < 1000
SU_GTHRDtex: FU_GTHRDtex + BU_GTHRDtex < 1
UF_GTHRDtex: F_GTHRDtex -1000 FU_GTHRDtex < 0
UR_GTHRDtex: R_GTHRDtex -1000 BU_GTHRDtex < 0
SU_H2Otex: FU_H2Otex + BU_H2Otex < 1
UF_H2Otex: F_H2Otex -1000 FU_H2Otex < 0
UR_H2Otex: R_H2Otex -1000 BU_H2Otex < 0
SU_H2Otpp: FU_H2Otpp + BU_H2Otpp < 1
UF_H2Otpp: F_H2Otpp -1000 FU_H2Otpp < 0
UR_H2Otpp: R_H2Otpp -1000 BU_H2Otpp < 0
G_HEX1: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_g6p_c -0.59209 LC_glc-D_c - DG_HEX1 + DGo_HEX1 = 0
FU_HEX1: DG_HEX1 +1000 FU_HEX1 < 1000
BU_HEX1: - DG_HEX1 +1000 BU_HEX1 < 1000
SU_HEX1: FU_HEX1 + BU_HEX1 < 1
UF_HEX1: F_HEX1 -1000 FU_HEX1 < 0
UR_HEX1: R_HEX1 -1000 BU_HEX1 < 0
G_HISt2rpp: +0.59209 LC_his-L_c -0.59209 LC_his-L_p - DG_HISt2rpp + DGo_HISt2rpp = 0
FU_HISt2rpp: DG_HISt2rpp +1000 FU_HISt2rpp < 1000
BU_HISt2rpp: - DG_HISt2rpp +1000 BU_HISt2rpp < 1000
SU_HISt2rpp: FU_HISt2rpp + BU_HISt2rpp < 1
UF_HISt2rpp: F_HISt2rpp -1000 FU_HISt2rpp < 0
UR_HISt2rpp: R_HISt2rpp -1000 BU_HISt2rpp < 0
G_HIStex: -0.59209 LC_his-L_e +0.59209 LC_his-L_p - DG_HIStex + DGo_HIStex = 0
FU_HIStex: DG_HIStex +1000 FU_HIStex < 1000
BU_HIStex: - DG_HIStex +1000 BU_HIStex < 1000
SU_HIStex: FU_HIStex + BU_HIStex < 1
UF_HIStex: F_HIStex -1000 FU_HIStex < 0
UR_HIStex: R_HIStex -1000 BU_HIStex < 0
G_Htex: - DG_Htex + DGo_Htex = 0
FU_Htex: DG_Htex +1000 FU_Htex < 1000
BU_Htex: - DG_Htex +1000 BU_Htex < 1000
SU_Htex: FU_Htex + BU_Htex < 1
UF_Htex: F_Htex -1000 FU_Htex < 0
UR_Htex: R_Htex -1000 BU_Htex < 0
G_ICDHyr: +0.59209 LC_akg_c +0.59209 LC_co2_c -0.59209 LC_icit_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_ICDHyr + DGo_ICDHyr = 0
FU_ICDHyr: DG_ICDHyr +1000 FU_ICDHyr < 1000
BU_ICDHyr: - DG_ICDHyr +1000 BU_ICDHyr < 1000
SU_ICDHyr: FU_ICDHyr + BU_ICDHyr < 1
UF_ICDHyr: F_ICDHyr -1000 FU_ICDHyr < 0
UR_ICDHyr: R_ICDHyr -1000 BU_ICDHyr < 0
G_ICL: +0.59209 LC_glx_c -0.59209 LC_icit_c +0.59209 LC_succ_c - DG_ICL + DGo_ICL = 0
FU_ICL: DG_ICL +1000 FU_ICL < 1000
BU_ICL: - DG_ICL +1000 BU_ICL < 1000
SU_ICL: FU_ICL + BU_ICL < 1
UF_ICL: F_ICL -1000 FU_ICL < 0
UR_ICL: R_ICL -1000 BU_ICL < 0
G_ILEt2rpp: +0.59209 LC_ile-L_c -0.59209 LC_ile-L_p - DG_ILEt2rpp + DGo_ILEt2rpp = 0
FU_ILEt2rpp: DG_ILEt2rpp +1000 FU_ILEt2rpp < 1000
BU_ILEt2rpp: - DG_ILEt2rpp +1000 BU_ILEt2rpp < 1000
SU_ILEt2rpp: FU_ILEt2rpp + BU_ILEt2rpp < 1
UF_ILEt2rpp: F_ILEt2rpp -1000 FU_ILEt2rpp < 0
UR_ILEt2rpp: R_ILEt2rpp -1000 BU_ILEt2rpp < 0
G_ILEtex: -0.59209 LC_ile-L_e +0.59209 LC_ile-L_p - DG_ILEtex + DGo_ILEtex = 0
FU_ILEtex: DG_ILEtex +1000 FU_ILEtex < 1000
BU_ILEtex: - DG_ILEtex +1000 BU_ILEtex < 1000
SU_ILEtex: FU_ILEtex + BU_ILEtex < 1
UF_ILEtex: F_ILEtex -1000 FU_ILEtex < 0
UR_ILEtex: R_ILEtex -1000 BU_ILEtex < 0
SU_Kabcpp: FU_Kabcpp + BU_Kabcpp < 1
UF_Kabcpp: F_Kabcpp -1000 FU_Kabcpp < 0
UR_Kabcpp: R_Kabcpp -1000 BU_Kabcpp < 0
SU_Kt2pp: FU_Kt2pp + BU_Kt2pp < 1
UF_Kt2pp: F_Kt2pp -1000 FU_Kt2pp < 0
UR_Kt2pp: R_Kt2pp -1000 BU_Kt2pp < 0
SU_Ktex: FU_Ktex + BU_Ktex < 1
UF_Ktex: F_Ktex -1000 FU_Ktex < 0
UR_Ktex: R_Ktex -1000 BU_Ktex < 0
G_LDH_D: -0.59209 LC_lac-D_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_pyr_c - DG_LDH_D + DGo_LDH_D = 0
FU_LDH_D: DG_LDH_D +1000 FU_LDH_D < 1000
BU_LDH_D: - DG_LDH_D +1000 BU_LDH_D < 1000
SU_LDH_D: FU_LDH_D + BU_LDH_D < 1
UF_LDH_D: F_LDH_D -1000 FU_LDH_D < 0
UR_LDH_D: R_LDH_D -1000 BU_LDH_D < 0
G_LDH_D2: -0.59209 LC_lac-D_c +0.59209 LC_pyr_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_LDH_D2 + DGo_LDH_D2 = 0
FU_LDH_D2: DG_LDH_D2 +1000 FU_LDH_D2 < 1000
BU_LDH_D2: - DG_LDH_D2 +1000 BU_LDH_D2 < 1000
SU_LDH_D2: FU_LDH_D2 + BU_LDH_D2 < 1
UF_LDH_D2: F_LDH_D2 -1000 FU_LDH_D2 < 0
UR_LDH_D2: R_LDH_D2 -1000 BU_LDH_D2 < 0
G_LEUt2rpp: +0.59209 LC_leu-L_c -0.59209 LC_leu-L_p - DG_LEUt2rpp + DGo_LEUt2rpp = 0
FU_LEUt2rpp: DG_LEUt2rpp +1000 FU_LEUt2rpp < 1000
BU_LEUt2rpp: - DG_LEUt2rpp +1000 BU_LEUt2rpp < 1000
SU_LEUt2rpp: FU_LEUt2rpp + BU_LEUt2rpp < 1
UF_LEUt2rpp: F_LEUt2rpp -1000 FU_LEUt2rpp < 0
UR_LEUt2rpp: R_LEUt2rpp -1000 BU_LEUt2rpp < 0
G_LEUtex: -0.59209 LC_leu-L_e +0.59209 LC_leu-L_p - DG_LEUtex + DGo_LEUtex = 0
FU_LEUtex: DG_LEUtex +1000 FU_LEUtex < 1000
BU_LEUtex: - DG_LEUtex +1000 BU_LEUtex < 1000
SU_LEUtex: FU_LEUtex + BU_LEUtex < 1
UF_LEUtex: F_LEUtex -1000 FU_LEUtex < 0
UR_LEUtex: R_LEUtex -1000 BU_LEUtex < 0
G_LYSt3pp: -0.59209 LC_lys-L_c +0.59209 LC_lys-L_p - DG_LYSt3pp + DGo_LYSt3pp = 0
FU_LYSt3pp: DG_LYSt3pp +1000 FU_LYSt3pp < 1000
BU_LYSt3pp: - DG_LYSt3pp +1000 BU_LYSt3pp < 1000
SU_LYSt3pp: FU_LYSt3pp + BU_LYSt3pp < 1
UF_LYSt3pp: F_LYSt3pp -1000 FU_LYSt3pp < 0
UR_LYSt3pp: R_LYSt3pp -1000 BU_LYSt3pp < 0
G_LYStex: -0.59209 LC_lys-L_e +0.59209 LC_lys-L_p - DG_LYStex + DGo_LYStex = 0
FU_LYStex: DG_LYStex +1000 FU_LYStex < 1000
BU_LYStex: - DG_LYStex +1000 BU_LYStex < 1000
SU_LYStex: FU_LYStex + BU_LYStex < 1
UF_LYStex: F_LYStex -1000 FU_LYStex < 0
UR_LYStex: R_LYStex -1000 BU_LYStex < 0
G_MALS: -0.59209 LC_accoa_c +0.59209 LC_coa_c -0.59209 LC_glx_c +0.59209 LC_mal-L_c - DG_MALS + DGo_MALS = 0
FU_MALS: DG_MALS +1000 FU_MALS < 1000
BU_MALS: - DG_MALS +1000 BU_MALS < 1000
SU_MALS: FU_MALS + BU_MALS < 1
UF_MALS: F_MALS -1000 FU_MALS < 0
UR_MALS: R_MALS -1000 BU_MALS < 0
G_MALt3pp: -0.59209 LC_mal-L_c +0.59209 LC_mal-L_p - DG_MALt3pp + DGo_MALt3pp = 0
FU_MALt3pp: DG_MALt3pp +1000 FU_MALt3pp < 1000
BU_MALt3pp: - DG_MALt3pp +1000 BU_MALt3pp < 1000
SU_MALt3pp: FU_MALt3pp + BU_MALt3pp < 1
UF_MALt3pp: F_MALt3pp -1000 FU_MALt3pp < 0
UR_MALt3pp: R_MALt3pp -1000 BU_MALt3pp < 0
G_MALtex: -0.59209 LC_mal-L_e +0.59209 LC_mal-L_p - DG_MALtex + DGo_MALtex = 0
FU_MALtex: DG_MALtex +1000 FU_MALtex < 1000
BU_MALtex: - DG_MALtex +1000 BU_MALtex < 1000
SU_MALtex: FU_MALtex + BU_MALtex < 1
UF_MALtex: F_MALtex -1000 FU_MALtex < 0
UR_MALtex: R_MALtex -1000 BU_MALtex < 0
G_MDH: -0.59209 LC_mal-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_oaa_c - DG_MDH + DGo_MDH = 0
FU_MDH: DG_MDH +1000 FU_MDH < 1000
BU_MDH: - DG_MDH +1000 BU_MDH < 1000
SU_MDH: FU_MDH + BU_MDH < 1
UF_MDH: F_MDH -1000 FU_MDH < 0
UR_MDH: R_MDH -1000 BU_MDH < 0
G_MDH2: -0.59209 LC_mal-L_c +0.59209 LC_oaa_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_MDH2 + DGo_MDH2 = 0
FU_MDH2: DG_MDH2 +1000 FU_MDH2 < 1000
BU_MDH2: - DG_MDH2 +1000 BU_MDH2 < 1000
SU_MDH2: FU_MDH2 + BU_MDH2 < 1
UF_MDH2: F_MDH2 -1000 FU_MDH2 < 0
UR_MDH2: R_MDH2 -1000 BU_MDH2 < 0
G_MDH3: -0.59209 LC_mal-L_c +0.59209 LC_mql8_c -0.59209 LC_mqn8_c +0.59209 LC_oaa_c - DG_MDH3 + DGo_MDH3 = 0
FU_MDH3: DG_MDH3 +1000 FU_MDH3 < 1000
BU_MDH3: - DG_MDH3 +1000 BU_MDH3 < 1000
SU_MDH3: FU_MDH3 + BU_MDH3 < 1
UF_MDH3: F_MDH3 -1000 FU_MDH3 < 0
UR_MDH3: R_MDH3 -1000 BU_MDH3 < 0
G_ME1: +0.59209 LC_co2_c -0.59209 LC_mal-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_pyr_c - DG_ME1 + DGo_ME1 = 0
FU_ME1: DG_ME1 +1000 FU_ME1 < 1000
BU_ME1: - DG_ME1 +1000 BU_ME1 < 1000
SU_ME1: FU_ME1 + BU_ME1 < 1
UF_ME1: F_ME1 -1000 FU_ME1 < 0
UR_ME1: R_ME1 -1000 BU_ME1 < 0
G_ME2: +0.59209 LC_co2_c -0.59209 LC_mal-L_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c +0.59209 LC_pyr_c - DG_ME2 + DGo_ME2 = 0
FU_ME2: DG_ME2 +1000 FU_ME2 < 1000
BU_ME2: - DG_ME2 +1000 BU_ME2 < 1000
SU_ME2: FU_ME2 + BU_ME2 < 1
UF_ME2: F_ME2 -1000 FU_ME2 < 0
UR_ME2: R_ME2 -1000 BU_ME2 < 0
G_MEOHtex: -0.59209 LC_meoh_e +0.59209 LC_meoh_p - DG_MEOHtex + DGo_MEOHtex = 0
FU_MEOHtex: DG_MEOHtex +1000 FU_MEOHtex < 1000
BU_MEOHtex: - DG_MEOHtex +1000 BU_MEOHtex < 1000
SU_MEOHtex: FU_MEOHtex + BU_MEOHtex < 1
UF_MEOHtex: F_MEOHtex -1000 FU_MEOHtex < 0
UR_MEOHtex: R_MEOHtex -1000 BU_MEOHtex < 0
SU_MG2t3_2pp: FU_MG2t3_2pp + BU_MG2t3_2pp < 1
UF_MG2t3_2pp: F_MG2t3_2pp -1000 FU_MG2t3_2pp < 0
UR_MG2t3_2pp: R_MG2t3_2pp -1000 BU_MG2t3_2pp < 0
SU_MG2tex: FU_MG2tex + BU_MG2tex < 1
UF_MG2tex: F_MG2tex -1000 FU_MG2tex < 0
UR_MG2tex: R_MG2tex -1000 BU_MG2tex < 0
SU_MG2tpp: FU_MG2tpp + BU_MG2tpp < 1
UF_MG2tpp: F_MG2tpp -1000 FU_MG2tpp < 0
UR_MG2tpp: R_MG2tpp -1000 BU_MG2tpp < 0
SU_MG2uabcpp: FU_MG2uabcpp + BU_MG2uabcpp < 1
UF_MG2uabcpp: F_MG2uabcpp -1000 FU_MG2uabcpp < 0
UR_MG2uabcpp: R_MG2uabcpp -1000 BU_MG2uabcpp < 0
SU_MN2tpp: FU_MN2tpp + BU_MN2tpp < 1
UF_MN2tpp: F_MN2tpp -1000 FU_MN2tpp < 0
UR_MN2tpp: R_MN2tpp -1000 BU_MN2tpp < 0
SU_MNt2pp: FU_MNt2pp + BU_MNt2pp < 1
UF_MNt2pp: F_MNt2pp -1000 FU_MNt2pp < 0
UR_MNt2pp: R_MNt2pp -1000 BU_MNt2pp < 0
SU_MNtex: FU_MNtex + BU_MNtex < 1
UF_MNtex: F_MNtex -1000 FU_MNtex < 0
UR_MNtex: R_MNtex -1000 BU_MNtex < 0
SU_MOBDabcpp: FU_MOBDabcpp + BU_MOBDabcpp < 1
UF_MOBDabcpp: F_MOBDabcpp -1000 FU_MOBDabcpp < 0
UR_MOBDabcpp: R_MOBDabcpp -1000 BU_MOBDabcpp < 0
SU_MOBDtex: FU_MOBDtex + BU_MOBDtex < 1
UF_MOBDtex: F_MOBDtex -1000 FU_MOBDtex < 0
UR_MOBDtex: R_MOBDtex -1000 BU_MOBDtex < 0
G_MOX: +0.59209 LC_h2o2_c -0.59209 LC_mal-L_c -0.59209 LC_o2_c +0.59209 LC_oaa_c - DG_MOX + DGo_MOX = 0
FU_MOX: DG_MOX +1000 FU_MOX < 1000
BU_MOX: - DG_MOX +1000 BU_MOX < 1000
SU_MOX: FU_MOX + BU_MOX < 1
UF_MOX: F_MOX -1000 FU_MOX < 0
UR_MOX: R_MOX -1000 BU_MOX < 0
G_NADH10: +0.59209 LC_mql8_c -0.59209 LC_mqn8_c +0.59209 LC_nad_c -0.59209 LC_nadh_c - DG_NADH10 + DGo_NADH10 = 0
FU_NADH10: DG_NADH10 +1000 FU_NADH10 < 1000
BU_NADH10: - DG_NADH10 +1000 BU_NADH10 < 1000
SU_NADH10: FU_NADH10 + BU_NADH10 < 1
UF_NADH10: F_NADH10 -1000 FU_NADH10 < 0
UR_NADH10: R_NADH10 -1000 BU_NADH10 < 0
G_NADH16pp: +0.59209 LC_nad_c -0.59209 LC_nadh_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_NADH16pp + DGo_NADH16pp = 0
FU_NADH16pp: DG_NADH16pp +1000 FU_NADH16pp < 1000
BU_NADH16pp: - DG_NADH16pp +1000 BU_NADH16pp < 1000
SU_NADH16pp: FU_NADH16pp + BU_NADH16pp < 1
UF_NADH16pp: F_NADH16pp -1000 FU_NADH16pp < 0
UR_NADH16pp: R_NADH16pp -1000 BU_NADH16pp < 0
G_NADH17pp: +0.59209 LC_mql8_c -0.59209 LC_mqn8_c +0.59209 LC_nad_c -0.59209 LC_nadh_c - DG_NADH17pp + DGo_NADH17pp = 0
FU_NADH17pp: DG_NADH17pp +1000 FU_NADH17pp < 1000
BU_NADH17pp: - DG_NADH17pp +1000 BU_NADH17pp < 1000
SU_NADH17pp: FU_NADH17pp + BU_NADH17pp < 1
UF_NADH17pp: F_NADH17pp -1000 FU_NADH17pp < 0
UR_NADH17pp: R_NADH17pp -1000 BU_NADH17pp < 0
G_NADH18pp: -0.59209 LC_2dmmq8_c +0.59209 LC_2dmmql8_c +0.59209 LC_nad_c -0.59209 LC_nadh_c - DG_NADH18pp + DGo_NADH18pp = 0
FU_NADH18pp: DG_NADH18pp +1000 FU_NADH18pp < 1000
BU_NADH18pp: - DG_NADH18pp +1000 BU_NADH18pp < 1000
SU_NADH18pp: FU_NADH18pp + BU_NADH18pp < 1
UF_NADH18pp: F_NADH18pp -1000 FU_NADH18pp < 0
UR_NADH18pp: R_NADH18pp -1000 BU_NADH18pp < 0
G_NADH5: +0.59209 LC_nad_c -0.59209 LC_nadh_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_NADH5 + DGo_NADH5 = 0
FU_NADH5: DG_NADH5 +1000 FU_NADH5 < 1000
BU_NADH5: - DG_NADH5 +1000 BU_NADH5 < 1000
SU_NADH5: FU_NADH5 + BU_NADH5 < 1
UF_NADH5: F_NADH5 -1000 FU_NADH5 < 0
UR_NADH5: R_NADH5 -1000 BU_NADH5 < 0
G_NADH9: -0.59209 LC_2dmmq8_c +0.59209 LC_2dmmql8_c +0.59209 LC_nad_c -0.59209 LC_nadh_c - DG_NADH9 + DGo_NADH9 = 0
FU_NADH9: DG_NADH9 +1000 FU_NADH9 < 1000
BU_NADH9: - DG_NADH9 +1000 BU_NADH9 < 1000
SU_NADH9: FU_NADH9 + BU_NADH9 < 1
UF_NADH9: F_NADH9 -1000 FU_NADH9 < 0
UR_NADH9: R_NADH9 -1000 BU_NADH9 < 0
G_NADK: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_nad_c +0.59209 LC_nadp_c - DG_NADK + DGo_NADK = 0
FU_NADK: DG_NADK +1000 FU_NADK < 1000
BU_NADK: - DG_NADK +1000 BU_NADK < 1000
SU_NADK: FU_NADK + BU_NADK < 1
UF_NADK: F_NADK -1000 FU_NADK < 0
UR_NADK: R_NADK -1000 BU_NADK < 0
G_NADPHQR2: +0.59209 LC_nadp_c -0.59209 LC_nadph_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_NADPHQR2 + DGo_NADPHQR2 = 0
FU_NADPHQR2: DG_NADPHQR2 +1000 FU_NADPHQR2 < 1000
BU_NADPHQR2: - DG_NADPHQR2 +1000 BU_NADPHQR2 < 1000
SU_NADPHQR2: FU_NADPHQR2 + BU_NADPHQR2 < 1
UF_NADPHQR2: F_NADPHQR2 -1000 FU_NADPHQR2 < 0
UR_NADPHQR2: R_NADPHQR2 -1000 BU_NADPHQR2 < 0
G_NADPHQR3: +0.59209 LC_mql8_c -0.59209 LC_mqn8_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c - DG_NADPHQR3 + DGo_NADPHQR3 = 0
FU_NADPHQR3: DG_NADPHQR3 +1000 FU_NADPHQR3 < 1000
BU_NADPHQR3: - DG_NADPHQR3 +1000 BU_NADPHQR3 < 1000
SU_NADPHQR3: FU_NADPHQR3 + BU_NADPHQR3 < 1
UF_NADPHQR3: F_NADPHQR3 -1000 FU_NADPHQR3 < 0
UR_NADPHQR3: R_NADPHQR3 -1000 BU_NADPHQR3 < 0
G_NADPHQR4: -0.59209 LC_2dmmq8_c +0.59209 LC_2dmmql8_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c - DG_NADPHQR4 + DGo_NADPHQR4 = 0
FU_NADPHQR4: DG_NADPHQR4 +1000 FU_NADPHQR4 < 1000
BU_NADPHQR4: - DG_NADPHQR4 +1000 BU_NADPHQR4 < 1000
SU_NADPHQR4: FU_NADPHQR4 + BU_NADPHQR4 < 1
UF_NADPHQR4: F_NADPHQR4 -1000 FU_NADPHQR4 < 0
UR_NADPHQR4: R_NADPHQR4 -1000 BU_NADPHQR4 < 0
G_NADPPPS: +0.59209 LC_nad_c -0.59209 LC_nadp_c +0.59209 LC_pi_c - DG_NADPPPS + DGo_NADPPPS = 0
FU_NADPPPS: DG_NADPPPS +1000 FU_NADPPPS < 1000
BU_NADPPPS: - DG_NADPPPS +1000 BU_NADPPPS < 1000
SU_NADPPPS: FU_NADPPPS + BU_NADPPPS < 1
UF_NADPPPS: F_NADPPPS -1000 FU_NADPPPS < 0
UR_NADPPPS: R_NADPPPS -1000 BU_NADPPPS < 0
G_NADTRHD: -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c - DG_NADTRHD + DGo_NADTRHD = 0
FU_NADTRHD: DG_NADTRHD +1000 FU_NADTRHD < 1000
BU_NADTRHD: - DG_NADTRHD +1000 BU_NADTRHD < 1000
SU_NADTRHD: FU_NADTRHD + BU_NADTRHD < 1
UF_NADTRHD: F_NADTRHD -1000 FU_NADTRHD < 0
UR_NADTRHD: R_NADTRHD -1000 BU_NADTRHD < 0
SU_NAt3_1p5pp: FU_NAt3_1p5pp + BU_NAt3_1p5pp < 1
UF_NAt3_1p5pp: F_NAt3_1p5pp -1000 FU_NAt3_1p5pp < 0
UR_NAt3_1p5pp: R_NAt3_1p5pp -1000 BU_NAt3_1p5pp < 0
SU_NAt3_2pp: FU_NAt3_2pp + BU_NAt3_2pp < 1
UF_NAt3_2pp: F_NAt3_2pp -1000 FU_NAt3_2pp < 0
UR_NAt3_2pp: R_NAt3_2pp -1000 BU_NAt3_2pp < 0
SU_NAt3pp: FU_NAt3pp + BU_NAt3pp < 1
UF_NAt3pp: F_NAt3pp -1000 FU_NAt3pp < 0
UR_NAt3pp: R_NAt3pp -1000 BU_NAt3pp < 0
G_NH4tex: -0.59209 LC_nh4_e +0.59209 LC_nh4_p - DG_NH4tex + DGo_NH4tex = 0
FU_NH4tex: DG_NH4tex +1000 FU_NH4tex < 1000
BU_NH4tex: - DG_NH4tex +1000 BU_NH4tex < 1000
SU_NH4tex: FU_NH4tex + BU_NH4tex < 1
UF_NH4tex: F_NH4tex -1000 FU_NH4tex < 0
UR_NH4tex: R_NH4tex -1000 BU_NH4tex < 0
G_NH4tpp: +0.59209 LC_nh4_c -0.59209 LC_nh4_p - DG_NH4tpp + DGo_NH4tpp = 0
FU_NH4tpp: DG_NH4tpp +1000 FU_NH4tpp < 1000
BU_NH4tpp: - DG_NH4tpp +1000 BU_NH4tpp < 1000
SU_NH4tpp: FU_NH4tpp + BU_NH4tpp < 1
UF_NH4tpp: F_NH4tpp -1000 FU_NH4tpp < 0
UR_NH4tpp: R_NH4tpp -1000 BU_NH4tpp < 0
SU_NI2tex: FU_NI2tex + BU_NI2tex < 1
UF_NI2tex: F_NI2tex -1000 FU_NI2tex < 0
UR_NI2tex: R_NI2tex -1000 BU_NI2tex < 0
SU_NI2tpp: FU_NI2tpp + BU_NI2tpp < 1
UF_NI2tpp: F_NI2tpp -1000 FU_NI2tpp < 0
UR_NI2tpp: R_NI2tpp -1000 BU_NI2tpp < 0
SU_NI2uabcpp: FU_NI2uabcpp + BU_NI2uabcpp < 1
UF_NI2uabcpp: F_NI2uabcpp -1000 FU_NI2uabcpp < 0
UR_NI2uabcpp: R_NI2uabcpp -1000 BU_NI2uabcpp < 0
G_NTP1: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_pi_c - DG_NTP1 + DGo_NTP1 = 0
FU_NTP1: DG_NTP1 +1000 FU_NTP1 < 1000
BU_NTP1: - DG_NTP1 +1000 BU_NTP1 < 1000
SU_NTP1: FU_NTP1 + BU_NTP1 < 1
UF_NTP1: F_NTP1 -1000 FU_NTP1 < 0
UR_NTP1: R_NTP1 -1000 BU_NTP1 < 0
G_NTPP6: +0.59209 LC_amp_c -0.59209 LC_atp_c +0.59209 LC_ppi_c - DG_NTPP6 + DGo_NTPP6 = 0
FU_NTPP6: DG_NTPP6 +1000 FU_NTPP6 < 1000
BU_NTPP6: - DG_NTPP6 +1000 BU_NTPP6 < 1000
SU_NTPP6: FU_NTPP6 + BU_NTPP6 < 1
UF_NTPP6: F_NTPP6 -1000 FU_NTPP6 < 0
UR_NTPP6: R_NTPP6 -1000 BU_NTPP6 < 0
G_O2tex: -0.59209 LC_o2_e +0.59209 LC_o2_p - DG_O2tex + DGo_O2tex = 0
FU_O2tex: DG_O2tex +1000 FU_O2tex < 1000
BU_O2tex: - DG_O2tex +1000 BU_O2tex < 1000
SU_O2tex: FU_O2tex + BU_O2tex < 1
UF_O2tex: F_O2tex -1000 FU_O2tex < 0
UR_O2tex: R_O2tex -1000 BU_O2tex < 0
G_O2tpp: +0.59209 LC_o2_c -0.59209 LC_o2_p - DG_O2tpp + DGo_O2tpp = 0
FU_O2tpp: DG_O2tpp +1000 FU_O2tpp < 1000
BU_O2tpp: - DG_O2tpp +1000 BU_O2tpp < 1000
SU_O2tpp: FU_O2tpp + BU_O2tpp < 1
UF_O2tpp: F_O2tpp -1000 FU_O2tpp < 0
UR_O2tpp: R_O2tpp -1000 BU_O2tpp < 0
G_OAADC: +0.59209 LC_co2_c -0.59209 LC_oaa_c +0.59209 LC_pyr_c - DG_OAADC + DGo_OAADC = 0
FU_OAADC: DG_OAADC +1000 FU_OAADC < 1000
BU_OAADC: - DG_OAADC +1000 BU_OAADC < 1000
SU_OAADC: FU_OAADC + BU_OAADC < 1
UF_OAADC: F_OAADC -1000 FU_OAADC < 0
UR_OAADC: R_OAADC -1000 BU_OAADC < 0
G_PDH: +0.59209 LC_accoa_c +0.59209 LC_co2_c -0.59209 LC_coa_c -0.59209 LC_nad_c +0.59209 LC_nadh_c -0.59209 LC_pyr_c - DG_PDH + DGo_PDH = 0
FU_PDH: DG_PDH +1000 FU_PDH < 1000
BU_PDH: - DG_PDH +1000 BU_PDH < 1000
SU_PDH: FU_PDH + BU_PDH < 1
UF_PDH: F_PDH -1000 FU_PDH < 0
UR_PDH: R_PDH -1000 BU_PDH < 0
G_PE160abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_pe160_c +0.59209 LC_pi_c +0.59209 LC_pe160_p - DG_PE160abcpp + DGo_PE160abcpp = 0
FU_PE160abcpp: DG_PE160abcpp +1000 FU_PE160abcpp < 1000
BU_PE160abcpp: - DG_PE160abcpp +1000 BU_PE160abcpp < 1000
SU_PE160abcpp: FU_PE160abcpp + BU_PE160abcpp < 1
UF_PE160abcpp: F_PE160abcpp -1000 FU_PE160abcpp < 0
UR_PE160abcpp: R_PE160abcpp -1000 BU_PE160abcpp < 0
G_PE161abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_pe161_c +0.59209 LC_pi_c +0.59209 LC_pe161_p - DG_PE161abcpp + DGo_PE161abcpp = 0
FU_PE161abcpp: DG_PE161abcpp +1000 FU_PE161abcpp < 1000
BU_PE161abcpp: - DG_PE161abcpp +1000 BU_PE161abcpp < 1000
SU_PE161abcpp: FU_PE161abcpp + BU_PE161abcpp < 1
UF_PE161abcpp: F_PE161abcpp -1000 FU_PE161abcpp < 0
UR_PE161abcpp: R_PE161abcpp -1000 BU_PE161abcpp < 0
G_PE181abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_pe181_c +0.59209 LC_pi_c +0.59209 LC_pe181_p - DG_PE181abcpp + DGo_PE181abcpp = 0
FU_PE181abcpp: DG_PE181abcpp +1000 FU_PE181abcpp < 1000
BU_PE181abcpp: - DG_PE181abcpp +1000 BU_PE181abcpp < 1000
SU_PE181abcpp: FU_PE181abcpp + BU_PE181abcpp < 1
UF_PE181abcpp: F_PE181abcpp -1000 FU_PE181abcpp < 0
UR_PE181abcpp: R_PE181abcpp -1000 BU_PE181abcpp < 0
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
G_PFL: +0.59209 LC_accoa_c -0.59209 LC_coa_c +0.59209 LC_for_c -0.59209 LC_pyr_c - DG_PFL + DGo_PFL = 0
FU_PFL: DG_PFL +1000 FU_PFL < 1000
BU_PFL: - DG_PFL +1000 BU_PFL < 1000
SU_PFL: FU_PFL + BU_PFL < 1
UF_PFL: F_PFL -1000 FU_PFL < 0
UR_PFL: R_PFL -1000 BU_PFL < 0
G_PG160abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_pg160_c +0.59209 LC_pi_c +0.59209 LC_pg160_p - DG_PG160abcpp + DGo_PG160abcpp = 0
FU_PG160abcpp: DG_PG160abcpp +1000 FU_PG160abcpp < 1000
BU_PG160abcpp: - DG_PG160abcpp +1000 BU_PG160abcpp < 1000
SU_PG160abcpp: FU_PG160abcpp + BU_PG160abcpp < 1
UF_PG160abcpp: F_PG160abcpp -1000 FU_PG160abcpp < 0
UR_PG160abcpp: R_PG160abcpp -1000 BU_PG160abcpp < 0
G_PG161abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_pg161_c +0.59209 LC_pi_c +0.59209 LC_pg161_p - DG_PG161abcpp + DGo_PG161abcpp = 0
FU_PG161abcpp: DG_PG161abcpp +1000 FU_PG161abcpp < 1000
BU_PG161abcpp: - DG_PG161abcpp +1000 BU_PG161abcpp < 1000
SU_PG161abcpp: FU_PG161abcpp + BU_PG161abcpp < 1
UF_PG161abcpp: F_PG161abcpp -1000 FU_PG161abcpp < 0
UR_PG161abcpp: R_PG161abcpp -1000 BU_PG161abcpp < 0
G_PG181abcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_pg181_c +0.59209 LC_pi_c +0.59209 LC_pg181_p - DG_PG181abcpp + DGo_PG181abcpp = 0
FU_PG181abcpp: DG_PG181abcpp +1000 FU_PG181abcpp < 1000
BU_PG181abcpp: - DG_PG181abcpp +1000 BU_PG181abcpp < 1000
SU_PG181abcpp: FU_PG181abcpp + BU_PG181abcpp < 1
UF_PG181abcpp: F_PG181abcpp -1000 FU_PG181abcpp < 0
UR_PG181abcpp: R_PG181abcpp -1000 BU_PG181abcpp < 0
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
SU_PHEMEabcpp: FU_PHEMEabcpp + BU_PHEMEabcpp < 1
UF_PHEMEabcpp: F_PHEMEabcpp -1000 FU_PHEMEabcpp < 0
UR_PHEMEabcpp: R_PHEMEabcpp -1000 BU_PHEMEabcpp < 0
SU_PHEMEtiex: FU_PHEMEtiex + BU_PHEMEtiex < 1
UF_PHEMEtiex: F_PHEMEtiex -1000 FU_PHEMEtiex < 0
UR_PHEMEtiex: R_PHEMEtiex -1000 BU_PHEMEtiex < 0
G_PHEt2rpp: +0.59209 LC_phe-L_c -0.59209 LC_phe-L_p - DG_PHEt2rpp + DGo_PHEt2rpp = 0
FU_PHEt2rpp: DG_PHEt2rpp +1000 FU_PHEt2rpp < 1000
BU_PHEt2rpp: - DG_PHEt2rpp +1000 BU_PHEt2rpp < 1000
SU_PHEt2rpp: FU_PHEt2rpp + BU_PHEt2rpp < 1
UF_PHEt2rpp: F_PHEt2rpp -1000 FU_PHEt2rpp < 0
UR_PHEt2rpp: R_PHEt2rpp -1000 BU_PHEt2rpp < 0
G_PHEtex: -0.59209 LC_phe-L_e +0.59209 LC_phe-L_p - DG_PHEtex + DGo_PHEtex = 0
FU_PHEtex: DG_PHEtex +1000 FU_PHEtex < 1000
BU_PHEtex: - DG_PHEtex +1000 BU_PHEtex < 1000
SU_PHEtex: FU_PHEtex + BU_PHEtex < 1
UF_PHEtex: F_PHEtex -1000 FU_PHEtex < 0
UR_PHEtex: R_PHEtex -1000 BU_PHEtex < 0
G_PIt2rpp: +0.59209 LC_pi_c -0.59209 LC_pi_p - DG_PIt2rpp + DGo_PIt2rpp = 0
FU_PIt2rpp: DG_PIt2rpp +1000 FU_PIt2rpp < 1000
BU_PIt2rpp: - DG_PIt2rpp +1000 BU_PIt2rpp < 1000
SU_PIt2rpp: FU_PIt2rpp + BU_PIt2rpp < 1
UF_PIt2rpp: F_PIt2rpp -1000 FU_PIt2rpp < 0
UR_PIt2rpp: R_PIt2rpp -1000 BU_PIt2rpp < 0
G_PItex: -0.59209 LC_pi_e +0.59209 LC_pi_p - DG_PItex + DGo_PItex = 0
FU_PItex: DG_PItex +1000 FU_PItex < 1000
BU_PItex: - DG_PItex +1000 BU_PItex < 1000
SU_PItex: FU_PItex + BU_PItex < 1
UF_PItex: F_PItex -1000 FU_PItex < 0
UR_PItex: R_PItex -1000 BU_PItex < 0
G_PIuabcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c +1.1842 LC_pi_c - DG_PIuabcpp + DGo_PIuabcpp = 0
FU_PIuabcpp: DG_PIuabcpp +1000 FU_PIuabcpp < 1000
BU_PIuabcpp: - DG_PIuabcpp +1000 BU_PIuabcpp < 1000
SU_PIuabcpp: FU_PIuabcpp + BU_PIuabcpp < 1
UF_PIuabcpp: F_PIuabcpp -1000 FU_PIuabcpp < 0
UR_PIuabcpp: R_PIuabcpp -1000 BU_PIuabcpp < 0
SU_POR5: FU_POR5 + BU_POR5 < 1
UF_POR5: F_POR5 -1000 FU_POR5 < 0
UR_POR5: R_POR5 -1000 BU_POR5 < 0
G_POX: +0.59209 LC_ac_c +0.59209 LC_co2_c -0.59209 LC_pyr_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c - DG_POX + DGo_POX = 0
FU_POX: DG_POX +1000 FU_POX < 1000
BU_POX: - DG_POX +1000 BU_POX < 1000
SU_POX: FU_POX + BU_POX < 1
UF_POX: F_POX -1000 FU_POX < 0
UR_POX: R_POX -1000 BU_POX < 0
G_PPA: +1.1842 LC_pi_c -0.59209 LC_ppi_c - DG_PPA + DGo_PPA = 0
FU_PPA: DG_PPA +1000 FU_PPA < 1000
BU_PPA: - DG_PPA +1000 BU_PPA < 1000
SU_PPA: FU_PPA + BU_PPA < 1
UF_PPA: F_PPA -1000 FU_PPA < 0
UR_PPA: R_PPA -1000 BU_PPA < 0
SU_PPA2: FU_PPA2 + BU_PPA2 < 1
UF_PPA2: F_PPA2 -1000 FU_PPA2 < 0
UR_PPA2: R_PPA2 -1000 BU_PPA2 < 0
G_PPC: -0.59209 LC_co2_c +0.59209 LC_oaa_c -0.59209 LC_pep_c +0.59209 LC_pi_c - DG_PPC + DGo_PPC = 0
FU_PPC: DG_PPC +1000 FU_PPC < 1000
BU_PPC: - DG_PPC +1000 BU_PPC < 1000
SU_PPC: FU_PPC + BU_PPC < 1
UF_PPC: F_PPC -1000 FU_PPC < 0
UR_PPC: R_PPC -1000 BU_PPC < 0
G_PPCK: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_oaa_c +0.59209 LC_pep_c - DG_PPCK + DGo_PPCK = 0
FU_PPCK: DG_PPCK +1000 FU_PPCK < 1000
BU_PPCK: - DG_PPCK +1000 BU_PPCK < 1000
SU_PPCK: FU_PPCK + BU_PPCK < 1
UF_PPCK: F_PPCK -1000 FU_PPCK < 0
UR_PPCK: R_PPCK -1000 BU_PPCK < 0
SU_PPK2r: FU_PPK2r + BU_PPK2r < 1
UF_PPK2r: F_PPK2r -1000 FU_PPK2r < 0
UR_PPK2r: R_PPK2r -1000 BU_PPK2r < 0
G_PPKr: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_pi_c +0.59209 LC_ppi_c - DG_PPKr + DGo_PPKr = 0
FU_PPKr: DG_PPKr +1000 FU_PPKr < 1000
BU_PPKr: - DG_PPKr +1000 BU_PPKr < 1000
SU_PPKr: FU_PPKr + BU_PPKr < 1
UF_PPKr: F_PPKr -1000 FU_PPKr < 0
UR_PPKr: R_PPKr -1000 BU_PPKr < 0
SU_PPPGO: FU_PPPGO + BU_PPPGO < 1
UF_PPPGO: F_PPPGO -1000 FU_PPPGO < 0
UR_PPPGO: R_PPPGO -1000 BU_PPPGO < 0
SU_PPPGO3: FU_PPPGO3 + BU_PPPGO3 < 1
UF_PPPGO3: F_PPPGO3 -1000 FU_PPPGO3 < 0
UR_PPPGO3: R_PPPGO3 -1000 BU_PPPGO3 < 0
G_PPS: +0.59209 LC_amp_c -0.59209 LC_atp_c +0.59209 LC_pep_c +0.59209 LC_pi_c -0.59209 LC_pyr_c - DG_PPS + DGo_PPS = 0
FU_PPS: DG_PPS +1000 FU_PPS < 1000
BU_PPS: - DG_PPS +1000 BU_PPS < 1000
SU_PPS: FU_PPS + BU_PPS < 1
UF_PPS: F_PPS -1000 FU_PPS < 0
UR_PPS: R_PPS -1000 BU_PPS < 0
G_PROt2rpp: +0.59209 LC_pro-L_c -0.59209 LC_pro-L_p - DG_PROt2rpp + DGo_PROt2rpp = 0
FU_PROt2rpp: DG_PROt2rpp +1000 FU_PROt2rpp < 1000
BU_PROt2rpp: - DG_PROt2rpp +1000 BU_PROt2rpp < 1000
SU_PROt2rpp: FU_PROt2rpp + BU_PROt2rpp < 1
UF_PROt2rpp: F_PROt2rpp -1000 FU_PROt2rpp < 0
UR_PROt2rpp: R_PROt2rpp -1000 BU_PROt2rpp < 0
G_PROtex: -0.59209 LC_pro-L_e +0.59209 LC_pro-L_p - DG_PROtex + DGo_PROtex = 0
FU_PROtex: DG_PROtex +1000 FU_PROtex < 1000
BU_PROtex: - DG_PROtex +1000 BU_PROtex < 1000
SU_PROtex: FU_PROtex + BU_PROtex < 1
UF_PROtex: F_PROtex -1000 FU_PROtex < 0
UR_PROtex: R_PROtex -1000 BU_PROtex < 0
G_PTAr: -0.59209 LC_accoa_c +0.59209 LC_actp_c +0.59209 LC_coa_c -0.59209 LC_pi_c - DG_PTAr + DGo_PTAr = 0
FU_PTAr: DG_PTAr +1000 FU_PTAr < 1000
BU_PTAr: - DG_PTAr +1000 BU_PTAr < 1000
SU_PTAr: FU_PTAr + BU_PTAr < 1
UF_PTAr: F_PTAr -1000 FU_PTAr < 0
UR_PTAr: R_PTAr -1000 BU_PTAr < 0
G_PYK: -0.59209 LC_adp_c +0.59209 LC_atp_c -0.59209 LC_pep_c +0.59209 LC_pyr_c - DG_PYK + DGo_PYK = 0
FU_PYK: DG_PYK +1000 FU_PYK < 1000
BU_PYK: - DG_PYK +1000 BU_PYK < 1000
SU_PYK: FU_PYK + BU_PYK < 1
UF_PYK: F_PYK -1000 FU_PYK < 0
UR_PYK: R_PYK -1000 BU_PYK < 0
G_PYRt2rpp: +0.59209 LC_pyr_c -0.59209 LC_pyr_p - DG_PYRt2rpp + DGo_PYRt2rpp = 0
FU_PYRt2rpp: DG_PYRt2rpp +1000 FU_PYRt2rpp < 1000
BU_PYRt2rpp: - DG_PYRt2rpp +1000 BU_PYRt2rpp < 1000
SU_PYRt2rpp: FU_PYRt2rpp + BU_PYRt2rpp < 1
UF_PYRt2rpp: F_PYRt2rpp -1000 FU_PYRt2rpp < 0
UR_PYRt2rpp: R_PYRt2rpp -1000 BU_PYRt2rpp < 0
G_PYRtex: -0.59209 LC_pyr_e +0.59209 LC_pyr_p - DG_PYRtex + DGo_PYRtex = 0
FU_PYRtex: DG_PYRtex +1000 FU_PYRtex < 1000
BU_PYRtex: - DG_PYRtex +1000 BU_PYRtex < 1000
SU_PYRtex: FU_PYRtex + BU_PYRtex < 1
UF_PYRtex: F_PYRtex -1000 FU_PYRtex < 0
UR_PYRtex: R_PYRtex -1000 BU_PYRtex < 0
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
G_SERt2rpp: +0.59209 LC_ser-L_c -0.59209 LC_ser-L_p - DG_SERt2rpp + DGo_SERt2rpp = 0
FU_SERt2rpp: DG_SERt2rpp +1000 FU_SERt2rpp < 1000
BU_SERt2rpp: - DG_SERt2rpp +1000 BU_SERt2rpp < 1000
SU_SERt2rpp: FU_SERt2rpp + BU_SERt2rpp < 1
UF_SERt2rpp: F_SERt2rpp -1000 FU_SERt2rpp < 0
UR_SERt2rpp: R_SERt2rpp -1000 BU_SERt2rpp < 0
G_SERtex: -0.59209 LC_ser-L_e +0.59209 LC_ser-L_p - DG_SERtex + DGo_SERtex = 0
FU_SERtex: DG_SERtex +1000 FU_SERtex < 1000
BU_SERtex: - DG_SERtex +1000 BU_SERtex < 1000
SU_SERtex: FU_SERtex + BU_SERtex < 1
UF_SERtex: F_SERtex -1000 FU_SERtex < 0
UR_SERtex: R_SERtex -1000 BU_SERtex < 0
G_SO4t2pp: +0.59209 LC_so4_c -0.59209 LC_so4_p - DG_SO4t2pp + DGo_SO4t2pp = 0
FU_SO4t2pp: DG_SO4t2pp +1000 FU_SO4t2pp < 1000
BU_SO4t2pp: - DG_SO4t2pp +1000 BU_SO4t2pp < 1000
SU_SO4t2pp: FU_SO4t2pp + BU_SO4t2pp < 1
UF_SO4t2pp: F_SO4t2pp -1000 FU_SO4t2pp < 0
UR_SO4t2pp: R_SO4t2pp -1000 BU_SO4t2pp < 0
G_SO4tex: -0.59209 LC_so4_e +0.59209 LC_so4_p - DG_SO4tex + DGo_SO4tex = 0
FU_SO4tex: DG_SO4tex +1000 FU_SO4tex < 1000
BU_SO4tex: - DG_SO4tex +1000 BU_SO4tex < 1000
SU_SO4tex: FU_SO4tex + BU_SO4tex < 1
UF_SO4tex: F_SO4tex -1000 FU_SO4tex < 0
UR_SO4tex: R_SO4tex -1000 BU_SO4tex < 0
G_SPMDt3pp: -0.59209 LC_spmd_c +0.59209 LC_spmd_p - DG_SPMDt3pp + DGo_SPMDt3pp = 0
FU_SPMDt3pp: DG_SPMDt3pp +1000 FU_SPMDt3pp < 1000
BU_SPMDt3pp: - DG_SPMDt3pp +1000 BU_SPMDt3pp < 1000
SU_SPMDt3pp: FU_SPMDt3pp + BU_SPMDt3pp < 1
UF_SPMDt3pp: F_SPMDt3pp -1000 FU_SPMDt3pp < 0
UR_SPMDt3pp: R_SPMDt3pp -1000 BU_SPMDt3pp < 0
G_SPMDtex: -0.59209 LC_spmd_e +0.59209 LC_spmd_p - DG_SPMDtex + DGo_SPMDtex = 0
FU_SPMDtex: DG_SPMDtex +1000 FU_SPMDtex < 1000
BU_SPMDtex: - DG_SPMDtex +1000 BU_SPMDtex < 1000
SU_SPMDtex: FU_SPMDtex + BU_SPMDtex < 1
UF_SPMDtex: F_SPMDtex -1000 FU_SPMDtex < 0
UR_SPMDtex: R_SPMDtex -1000 BU_SPMDtex < 0
G_SUCCt2_2pp: +0.59209 LC_succ_c -0.59209 LC_succ_p - DG_SUCCt2_2pp + DGo_SUCCt2_2pp = 0
FU_SUCCt2_2pp: DG_SUCCt2_2pp +1000 FU_SUCCt2_2pp < 1000
BU_SUCCt2_2pp: - DG_SUCCt2_2pp +1000 BU_SUCCt2_2pp < 1000
SU_SUCCt2_2pp: FU_SUCCt2_2pp + BU_SUCCt2_2pp < 1
UF_SUCCt2_2pp: F_SUCCt2_2pp -1000 FU_SUCCt2_2pp < 0
UR_SUCCt2_2pp: R_SUCCt2_2pp -1000 BU_SUCCt2_2pp < 0
G_SUCCt2_3pp: +0.59209 LC_succ_c -0.59209 LC_succ_p - DG_SUCCt2_3pp + DGo_SUCCt2_3pp = 0
FU_SUCCt2_3pp: DG_SUCCt2_3pp +1000 FU_SUCCt2_3pp < 1000
BU_SUCCt2_3pp: - DG_SUCCt2_3pp +1000 BU_SUCCt2_3pp < 1000
SU_SUCCt2_3pp: FU_SUCCt2_3pp + BU_SUCCt2_3pp < 1
UF_SUCCt2_3pp: F_SUCCt2_3pp -1000 FU_SUCCt2_3pp < 0
UR_SUCCt2_3pp: R_SUCCt2_3pp -1000 BU_SUCCt2_3pp < 0
G_SUCCt3pp: -0.59209 LC_succ_c +0.59209 LC_succ_p - DG_SUCCt3pp + DGo_SUCCt3pp = 0
FU_SUCCt3pp: DG_SUCCt3pp +1000 FU_SUCCt3pp < 1000
BU_SUCCt3pp: - DG_SUCCt3pp +1000 BU_SUCCt3pp < 1000
SU_SUCCt3pp: FU_SUCCt3pp + BU_SUCCt3pp < 1
UF_SUCCt3pp: F_SUCCt3pp -1000 FU_SUCCt3pp < 0
UR_SUCCt3pp: R_SUCCt3pp -1000 BU_SUCCt3pp < 0
G_SUCCtex: -0.59209 LC_succ_e +0.59209 LC_succ_p - DG_SUCCtex + DGo_SUCCtex = 0
FU_SUCCtex: DG_SUCCtex +1000 FU_SUCCtex < 1000
BU_SUCCtex: - DG_SUCCtex +1000 BU_SUCCtex < 1000
SU_SUCCtex: FU_SUCCtex + BU_SUCCtex < 1
UF_SUCCtex: F_SUCCtex -1000 FU_SUCCtex < 0
UR_SUCCtex: R_SUCCtex -1000 BU_SUCCtex < 0
G_SUCDi: +0.59209 LC_fum_c -0.59209 LC_q8_c +0.59209 LC_q8h2_c -0.59209 LC_succ_c - DG_SUCDi + DGo_SUCDi = 0
FU_SUCDi: DG_SUCDi +1000 FU_SUCDi < 1000
BU_SUCDi: - DG_SUCDi +1000 BU_SUCDi < 1000
SU_SUCDi: FU_SUCDi + BU_SUCDi < 1
UF_SUCDi: F_SUCDi -1000 FU_SUCDi < 0
UR_SUCDi: R_SUCDi -1000 BU_SUCDi < 0
G_SUCOAS: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_coa_c +0.59209 LC_pi_c -0.59209 LC_succ_c +0.59209 LC_succoa_c - DG_SUCOAS + DGo_SUCOAS = 0
FU_SUCOAS: DG_SUCOAS +1000 FU_SUCOAS < 1000
BU_SUCOAS: - DG_SUCOAS +1000 BU_SUCOAS < 1000
SU_SUCOAS: FU_SUCOAS + BU_SUCOAS < 1
UF_SUCOAS: F_SUCOAS -1000 FU_SUCOAS < 0
UR_SUCOAS: R_SUCOAS -1000 BU_SUCOAS < 0
G_SULabcpp: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_pi_c +0.59209 LC_so4_c -0.59209 LC_so4_p - DG_SULabcpp + DGo_SULabcpp = 0
FU_SULabcpp: DG_SULabcpp +1000 FU_SULabcpp < 1000
BU_SULabcpp: - DG_SULabcpp +1000 BU_SULabcpp < 1000
SU_SULabcpp: FU_SULabcpp + BU_SULabcpp < 1
UF_SULabcpp: F_SULabcpp -1000 FU_SULabcpp < 0
UR_SULabcpp: R_SULabcpp -1000 BU_SULabcpp < 0
G_TALA: +0.59209 LC_e4p_c +0.59209 LC_f6p_c -0.59209 LC_g3p_c -0.59209 LC_s7p_c - DG_TALA + DGo_TALA = 0
FU_TALA: DG_TALA +1000 FU_TALA < 1000
BU_TALA: - DG_TALA +1000 BU_TALA < 1000
SU_TALA: FU_TALA + BU_TALA < 1
UF_TALA: F_TALA -1000 FU_TALA < 0
UR_TALA: R_TALA -1000 BU_TALA < 0
G_THD2pp: +0.59209 LC_nad_c -0.59209 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c - DG_THD2pp + DGo_THD2pp = 0
FU_THD2pp: DG_THD2pp +1000 FU_THD2pp < 1000
BU_THD2pp: - DG_THD2pp +1000 BU_THD2pp < 1000
SU_THD2pp: FU_THD2pp + BU_THD2pp < 1
UF_THD2pp: F_THD2pp -1000 FU_THD2pp < 0
UR_THD2pp: R_THD2pp -1000 BU_THD2pp < 0
G_THIORDXi: -0.59209 LC_h2o2_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_THIORDXi + DGo_THIORDXi = 0
FU_THIORDXi: DG_THIORDXi +1000 FU_THIORDXi < 1000
BU_THIORDXi: - DG_THIORDXi +1000 BU_THIORDXi < 1000
SU_THIORDXi: FU_THIORDXi + BU_THIORDXi < 1
UF_THIORDXi: F_THIORDXi -1000 FU_THIORDXi < 0
UR_THIORDXi: R_THIORDXi -1000 BU_THIORDXi < 0
G_THRt2pp: -0.59209 LC_thr-L_c +0.59209 LC_thr-L_p - DG_THRt2pp + DGo_THRt2pp = 0
FU_THRt2pp: DG_THRt2pp +1000 FU_THRt2pp < 1000
BU_THRt2pp: - DG_THRt2pp +1000 BU_THRt2pp < 1000
SU_THRt2pp: FU_THRt2pp + BU_THRt2pp < 1
UF_THRt2pp: F_THRt2pp -1000 FU_THRt2pp < 0
UR_THRt2pp: R_THRt2pp -1000 BU_THRt2pp < 0
G_THRt2rpp: +0.59209 LC_thr-L_c -0.59209 LC_thr-L_p - DG_THRt2rpp + DGo_THRt2rpp = 0
FU_THRt2rpp: DG_THRt2rpp +1000 FU_THRt2rpp < 1000
BU_THRt2rpp: - DG_THRt2rpp +1000 BU_THRt2rpp < 1000
SU_THRt2rpp: FU_THRt2rpp + BU_THRt2rpp < 1
UF_THRt2rpp: F_THRt2rpp -1000 FU_THRt2rpp < 0
UR_THRt2rpp: R_THRt2rpp -1000 BU_THRt2rpp < 0
G_THRtex: -0.59209 LC_thr-L_e +0.59209 LC_thr-L_p - DG_THRtex + DGo_THRtex = 0
FU_THRtex: DG_THRtex +1000 FU_THRtex < 1000
BU_THRtex: - DG_THRtex +1000 BU_THRtex < 1000
SU_THRtex: FU_THRtex + BU_THRtex < 1
UF_THRtex: F_THRtex -1000 FU_THRtex < 0
UR_THRtex: R_THRtex -1000 BU_THRtex < 0
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
G_TRDR: +0.59209 LC_nadp_c -0.59209 LC_nadph_c -0.59209 LC_trdox_c +0.59209 LC_trdrd_c - DG_TRDR + DGo_TRDR = 0
FU_TRDR: DG_TRDR +1000 FU_TRDR < 1000
BU_TRDR: - DG_TRDR +1000 BU_TRDR < 1000
SU_TRDR: FU_TRDR + BU_TRDR < 1
UF_TRDR: F_TRDR -1000 FU_TRDR < 0
UR_TRDR: R_TRDR -1000 BU_TRDR < 0
G_TRPt2rpp: +0.59209 LC_trp-L_c -0.59209 LC_trp-L_p - DG_TRPt2rpp + DGo_TRPt2rpp = 0
FU_TRPt2rpp: DG_TRPt2rpp +1000 FU_TRPt2rpp < 1000
BU_TRPt2rpp: - DG_TRPt2rpp +1000 BU_TRPt2rpp < 1000
SU_TRPt2rpp: FU_TRPt2rpp + BU_TRPt2rpp < 1
UF_TRPt2rpp: F_TRPt2rpp -1000 FU_TRPt2rpp < 0
UR_TRPt2rpp: R_TRPt2rpp -1000 BU_TRPt2rpp < 0
G_TRPtex: -0.59209 LC_trp-L_e +0.59209 LC_trp-L_p - DG_TRPtex + DGo_TRPtex = 0
FU_TRPtex: DG_TRPtex +1000 FU_TRPtex < 1000
BU_TRPtex: - DG_TRPtex +1000 BU_TRPtex < 1000
SU_TRPtex: FU_TRPtex + BU_TRPtex < 1
UF_TRPtex: F_TRPtex -1000 FU_TRPtex < 0
UR_TRPtex: R_TRPtex -1000 BU_TRPtex < 0
G_TRSARr: -0.59209 LC_2h3oppan_c +0.59209 LC_glyc-R_c +0.59209 LC_nad_c -0.59209 LC_nadh_c - DG_TRSARr + DGo_TRSARr = 0
FU_TRSARr: DG_TRSARr +1000 FU_TRSARr < 1000
BU_TRSARr: - DG_TRSARr +1000 BU_TRSARr < 1000
SU_TRSARr: FU_TRSARr + BU_TRSARr < 1
UF_TRSARr: F_TRSARr -1000 FU_TRSARr < 0
UR_TRSARr: R_TRSARr -1000 BU_TRSARr < 0
G_TYRt2rpp: +0.59209 LC_tyr-L_c -0.59209 LC_tyr-L_p - DG_TYRt2rpp + DGo_TYRt2rpp = 0
FU_TYRt2rpp: DG_TYRt2rpp +1000 FU_TYRt2rpp < 1000
BU_TYRt2rpp: - DG_TYRt2rpp +1000 BU_TYRt2rpp < 1000
SU_TYRt2rpp: FU_TYRt2rpp + BU_TYRt2rpp < 1
UF_TYRt2rpp: F_TYRt2rpp -1000 FU_TYRt2rpp < 0
UR_TYRt2rpp: R_TYRt2rpp -1000 BU_TYRt2rpp < 0
G_TYRtex: -0.59209 LC_tyr-L_e +0.59209 LC_tyr-L_p - DG_TYRtex + DGo_TYRtex = 0
FU_TYRtex: DG_TYRtex +1000 FU_TYRtex < 1000
BU_TYRtex: - DG_TYRtex +1000 BU_TYRtex < 1000
SU_TYRtex: FU_TYRtex + BU_TYRtex < 1
UF_TYRtex: F_TYRtex -1000 FU_TYRtex < 0
UR_TYRtex: R_TYRtex -1000 BU_TYRtex < 0
G_VALt2rpp: +0.59209 LC_val-L_c -0.59209 LC_val-L_p - DG_VALt2rpp + DGo_VALt2rpp = 0
FU_VALt2rpp: DG_VALt2rpp +1000 FU_VALt2rpp < 1000
BU_VALt2rpp: - DG_VALt2rpp +1000 BU_VALt2rpp < 1000
SU_VALt2rpp: FU_VALt2rpp + BU_VALt2rpp < 1
UF_VALt2rpp: F_VALt2rpp -1000 FU_VALt2rpp < 0
UR_VALt2rpp: R_VALt2rpp -1000 BU_VALt2rpp < 0
G_VALtex: -0.59209 LC_val-L_e +0.59209 LC_val-L_p - DG_VALtex + DGo_VALtex = 0
FU_VALtex: DG_VALtex +1000 FU_VALtex < 1000
BU_VALtex: - DG_VALtex +1000 BU_VALtex < 1000
SU_VALtex: FU_VALtex + BU_VALtex < 1
UF_VALtex: F_VALtex -1000 FU_VALtex < 0
UR_VALtex: R_VALtex -1000 BU_VALtex < 0
SU_ZN2tpp: FU_ZN2tpp + BU_ZN2tpp < 1
UF_ZN2tpp: F_ZN2tpp -1000 FU_ZN2tpp < 0
UR_ZN2tpp: R_ZN2tpp -1000 BU_ZN2tpp < 0
SU_ZNabcpp: FU_ZNabcpp + BU_ZNabcpp < 1
UF_ZNabcpp: F_ZNabcpp -1000 FU_ZNabcpp < 0
UR_ZNabcpp: R_ZNabcpp -1000 BU_ZNabcpp < 0
SU_Zn2tex: FU_Zn2tex + BU_Zn2tex < 1
UF_Zn2tex: F_Zn2tex -1000 FU_Zn2tex < 0
UR_Zn2tex: R_Zn2tex -1000 BU_Zn2tex < 0
G_LMPD_1_10fthf_c: +0.59209 LC_10fthf_c -1.1842 LC_3pg_c +7.1051 LC_adp_c +1.1842 LC_akg_c +2.3684 LC_amp_c -0.59209 LC_asp-L_c -9.4734 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_fum_c -1.7763 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_glyclt_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_1_10fthf_c + DGo_LMPD_1_10fthf_c = 0
FU_LMPD_1_10fthf_c: DG_LMPD_1_10fthf_c +1000 FU_LMPD_1_10fthf_c < 1000
BU_LMPD_1_10fthf_c: - DG_LMPD_1_10fthf_c +1000 BU_LMPD_1_10fthf_c < 1000
SU_LMPD_1_10fthf_c: FU_LMPD_1_10fthf_c + BU_LMPD_1_10fthf_c < 1
UF_LMPD_1_10fthf_c: F_LMPD_1_10fthf_c -1000 FU_LMPD_1_10fthf_c < 0
UR_LMPD_1_10fthf_c: R_LMPD_1_10fthf_c -1000 BU_LMPD_1_10fthf_c < 0
G_LMPD_2_10fthf_c: +0.59209 LC_10fthf_c -1.1842 LC_3pg_c +8.2893 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -10.0655 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_fum_c -1.7763 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_glyclt_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_2_10fthf_c + DGo_LMPD_2_10fthf_c = 0
FU_LMPD_2_10fthf_c: DG_LMPD_2_10fthf_c +1000 FU_LMPD_2_10fthf_c < 1000
BU_LMPD_2_10fthf_c: - DG_LMPD_2_10fthf_c +1000 BU_LMPD_2_10fthf_c < 1000
SU_LMPD_2_10fthf_c: FU_LMPD_2_10fthf_c + BU_LMPD_2_10fthf_c < 1
UF_LMPD_2_10fthf_c: F_LMPD_2_10fthf_c -1000 FU_LMPD_2_10fthf_c < 0
UR_LMPD_2_10fthf_c: R_LMPD_2_10fthf_c -1000 BU_LMPD_2_10fthf_c < 0
G_LMPD_3_10fthf_c: +0.59209 LC_10fthf_c -1.7763 LC_3pg_c +6.513 LC_adp_c +1.7763 LC_akg_c +2.3684 LC_amp_c -0.59209 LC_asp-L_c -8.8813 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_for_c +0.59209 LC_fum_c -2.3684 LC_glu-L_c +1.1842 LC_gly_c +0.59209 LC_glyclt_c -2.9604 LC_nad_c +2.9604 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_3_10fthf_c + DGo_LMPD_3_10fthf_c = 0
FU_LMPD_3_10fthf_c: DG_LMPD_3_10fthf_c +1000 FU_LMPD_3_10fthf_c < 1000
BU_LMPD_3_10fthf_c: - DG_LMPD_3_10fthf_c +1000 BU_LMPD_3_10fthf_c < 1000
SU_LMPD_3_10fthf_c: FU_LMPD_3_10fthf_c + BU_LMPD_3_10fthf_c < 1
UF_LMPD_3_10fthf_c: F_LMPD_3_10fthf_c -1000 FU_LMPD_3_10fthf_c < 0
UR_LMPD_3_10fthf_c: R_LMPD_3_10fthf_c -1000 BU_LMPD_3_10fthf_c < 0
G_LMPD_4_10fthf_c: +0.59209 LC_10fthf_c -1.7763 LC_3pg_c +7.6972 LC_adp_c +1.7763 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -9.4734 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_for_c +0.59209 LC_fum_c -2.3684 LC_glu-L_c +1.1842 LC_gly_c +0.59209 LC_glyclt_c -2.9604 LC_nad_c +2.9604 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_4_10fthf_c + DGo_LMPD_4_10fthf_c = 0
FU_LMPD_4_10fthf_c: DG_LMPD_4_10fthf_c +1000 FU_LMPD_4_10fthf_c < 1000
BU_LMPD_4_10fthf_c: - DG_LMPD_4_10fthf_c +1000 BU_LMPD_4_10fthf_c < 1000
SU_LMPD_4_10fthf_c: FU_LMPD_4_10fthf_c + BU_LMPD_4_10fthf_c < 1
UF_LMPD_4_10fthf_c: F_LMPD_4_10fthf_c -1000 FU_LMPD_4_10fthf_c < 0
UR_LMPD_4_10fthf_c: R_LMPD_4_10fthf_c -1000 BU_LMPD_4_10fthf_c < 0
SU_LMPD_5_2dmmql8_c: FU_LMPD_5_2dmmql8_c + BU_LMPD_5_2dmmql8_c < 1
UF_LMPD_5_2dmmql8_c: F_LMPD_5_2dmmql8_c -1000 FU_LMPD_5_2dmmql8_c < 0
UR_LMPD_5_2dmmql8_c: R_LMPD_5_2dmmql8_c -1000 BU_LMPD_5_2dmmql8_c < 0
SU_LMPD_6_2fe2s_c: FU_LMPD_6_2fe2s_c + BU_LMPD_6_2fe2s_c < 1
UF_LMPD_6_2fe2s_c: F_LMPD_6_2fe2s_c -1000 FU_LMPD_6_2fe2s_c < 0
UR_LMPD_6_2fe2s_c: R_LMPD_6_2fe2s_c -1000 BU_LMPD_6_2fe2s_c < 0
SU_LMPD_7_2fe2s_c: FU_LMPD_7_2fe2s_c + BU_LMPD_7_2fe2s_c < 1
UF_LMPD_7_2fe2s_c: F_LMPD_7_2fe2s_c -1000 FU_LMPD_7_2fe2s_c < 0
UR_LMPD_7_2fe2s_c: R_LMPD_7_2fe2s_c -1000 BU_LMPD_7_2fe2s_c < 0
SU_LMPD_8_2fe2s_c: FU_LMPD_8_2fe2s_c + BU_LMPD_8_2fe2s_c < 1
UF_LMPD_8_2fe2s_c: F_LMPD_8_2fe2s_c -1000 FU_LMPD_8_2fe2s_c < 0
UR_LMPD_8_2fe2s_c: R_LMPD_8_2fe2s_c -1000 BU_LMPD_8_2fe2s_c < 0
SU_LMPD_9_2fe2s_c: FU_LMPD_9_2fe2s_c + BU_LMPD_9_2fe2s_c < 1
UF_LMPD_9_2fe2s_c: F_LMPD_9_2fe2s_c -1000 FU_LMPD_9_2fe2s_c < 0
UR_LMPD_9_2fe2s_c: R_LMPD_9_2fe2s_c -1000 BU_LMPD_9_2fe2s_c < 0
SU_LMPD_10_2fe2s_c: FU_LMPD_10_2fe2s_c + BU_LMPD_10_2fe2s_c < 1
UF_LMPD_10_2fe2s_c: F_LMPD_10_2fe2s_c -1000 FU_LMPD_10_2fe2s_c < 0
UR_LMPD_10_2fe2s_c: R_LMPD_10_2fe2s_c -1000 BU_LMPD_10_2fe2s_c < 0
SU_LMPD_11_2fe2s_c: FU_LMPD_11_2fe2s_c + BU_LMPD_11_2fe2s_c < 1
UF_LMPD_11_2fe2s_c: F_LMPD_11_2fe2s_c -1000 FU_LMPD_11_2fe2s_c < 0
UR_LMPD_11_2fe2s_c: R_LMPD_11_2fe2s_c -1000 BU_LMPD_11_2fe2s_c < 0
SU_LMPD_12_2fe2s_c: FU_LMPD_12_2fe2s_c + BU_LMPD_12_2fe2s_c < 1
UF_LMPD_12_2fe2s_c: F_LMPD_12_2fe2s_c -1000 FU_LMPD_12_2fe2s_c < 0
UR_LMPD_12_2fe2s_c: R_LMPD_12_2fe2s_c -1000 BU_LMPD_12_2fe2s_c < 0
SU_LMPD_13_2fe2s_c: FU_LMPD_13_2fe2s_c + BU_LMPD_13_2fe2s_c < 1
UF_LMPD_13_2fe2s_c: F_LMPD_13_2fe2s_c -1000 FU_LMPD_13_2fe2s_c < 0
UR_LMPD_13_2fe2s_c: R_LMPD_13_2fe2s_c -1000 BU_LMPD_13_2fe2s_c < 0
SU_LMPD_14_4fe4s_c: FU_LMPD_14_4fe4s_c + BU_LMPD_14_4fe4s_c < 1
UF_LMPD_14_4fe4s_c: F_LMPD_14_4fe4s_c -1000 FU_LMPD_14_4fe4s_c < 0
UR_LMPD_14_4fe4s_c: R_LMPD_14_4fe4s_c -1000 BU_LMPD_14_4fe4s_c < 0
SU_LMPD_15_4fe4s_c: FU_LMPD_15_4fe4s_c + BU_LMPD_15_4fe4s_c < 1
UF_LMPD_15_4fe4s_c: F_LMPD_15_4fe4s_c -1000 FU_LMPD_15_4fe4s_c < 0
UR_LMPD_15_4fe4s_c: R_LMPD_15_4fe4s_c -1000 BU_LMPD_15_4fe4s_c < 0
SU_LMPD_16_4fe4s_c: FU_LMPD_16_4fe4s_c + BU_LMPD_16_4fe4s_c < 1
UF_LMPD_16_4fe4s_c: F_LMPD_16_4fe4s_c -1000 FU_LMPD_16_4fe4s_c < 0
UR_LMPD_16_4fe4s_c: R_LMPD_16_4fe4s_c -1000 BU_LMPD_16_4fe4s_c < 0
SU_LMPD_17_4fe4s_c: FU_LMPD_17_4fe4s_c + BU_LMPD_17_4fe4s_c < 1
UF_LMPD_17_4fe4s_c: F_LMPD_17_4fe4s_c -1000 FU_LMPD_17_4fe4s_c < 0
UR_LMPD_17_4fe4s_c: R_LMPD_17_4fe4s_c -1000 BU_LMPD_17_4fe4s_c < 0
SU_LMPD_18_4fe4s_c: FU_LMPD_18_4fe4s_c + BU_LMPD_18_4fe4s_c < 1
UF_LMPD_18_4fe4s_c: F_LMPD_18_4fe4s_c -1000 FU_LMPD_18_4fe4s_c < 0
UR_LMPD_18_4fe4s_c: R_LMPD_18_4fe4s_c -1000 BU_LMPD_18_4fe4s_c < 0
SU_LMPD_19_4fe4s_c: FU_LMPD_19_4fe4s_c + BU_LMPD_19_4fe4s_c < 1
UF_LMPD_19_4fe4s_c: F_LMPD_19_4fe4s_c -1000 FU_LMPD_19_4fe4s_c < 0
UR_LMPD_19_4fe4s_c: R_LMPD_19_4fe4s_c -1000 BU_LMPD_19_4fe4s_c < 0
SU_LMPD_20_4fe4s_c: FU_LMPD_20_4fe4s_c + BU_LMPD_20_4fe4s_c < 1
UF_LMPD_20_4fe4s_c: F_LMPD_20_4fe4s_c -1000 FU_LMPD_20_4fe4s_c < 0
UR_LMPD_20_4fe4s_c: R_LMPD_20_4fe4s_c -1000 BU_LMPD_20_4fe4s_c < 0
SU_LMPD_21_4fe4s_c: FU_LMPD_21_4fe4s_c + BU_LMPD_21_4fe4s_c < 1
UF_LMPD_21_4fe4s_c: F_LMPD_21_4fe4s_c -1000 FU_LMPD_21_4fe4s_c < 0
UR_LMPD_21_4fe4s_c: R_LMPD_21_4fe4s_c -1000 BU_LMPD_21_4fe4s_c < 0
G_LMPD_22_5mthf_c: -0.59209 LC_3pg_c +0.59209 LC_5mthf_c +7.6972 LC_adp_c +0.59209 LC_akg_c +2.3684 LC_amp_c -0.59209 LC_asp-L_c -10.0655 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_for_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_glyclt_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_22_5mthf_c + DGo_LMPD_22_5mthf_c = 0
FU_LMPD_22_5mthf_c: DG_LMPD_22_5mthf_c +1000 FU_LMPD_22_5mthf_c < 1000
BU_LMPD_22_5mthf_c: - DG_LMPD_22_5mthf_c +1000 BU_LMPD_22_5mthf_c < 1000
SU_LMPD_22_5mthf_c: FU_LMPD_22_5mthf_c + BU_LMPD_22_5mthf_c < 1
UF_LMPD_22_5mthf_c: F_LMPD_22_5mthf_c -1000 FU_LMPD_22_5mthf_c < 0
UR_LMPD_22_5mthf_c: R_LMPD_22_5mthf_c -1000 BU_LMPD_22_5mthf_c < 0
G_LMPD_23_5mthf_c: -0.59209 LC_3pg_c +0.59209 LC_5mthf_c +8.8813 LC_adp_c +0.59209 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -10.6576 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_for_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_glyclt_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_23_5mthf_c + DGo_LMPD_23_5mthf_c = 0
FU_LMPD_23_5mthf_c: DG_LMPD_23_5mthf_c +1000 FU_LMPD_23_5mthf_c < 1000
BU_LMPD_23_5mthf_c: - DG_LMPD_23_5mthf_c +1000 BU_LMPD_23_5mthf_c < 1000
SU_LMPD_23_5mthf_c: FU_LMPD_23_5mthf_c + BU_LMPD_23_5mthf_c < 1
UF_LMPD_23_5mthf_c: F_LMPD_23_5mthf_c -1000 FU_LMPD_23_5mthf_c < 0
UR_LMPD_23_5mthf_c: R_LMPD_23_5mthf_c -1000 BU_LMPD_23_5mthf_c < 0
G_LMPD_24_accoa_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +7.1051 LC_adp_c +1.1842 LC_akg_c +3.5525 LC_amp_c -1.7763 LC_asp-L_c -10.6576 LC_atp_c +1.1842 LC_co2_c +1.1842 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -1.1842 LC_nh4_c +8.2893 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_24_accoa_c + DGo_LMPD_24_accoa_c = 0
FU_LMPD_24_accoa_c: DG_LMPD_24_accoa_c +1000 FU_LMPD_24_accoa_c < 1000
BU_LMPD_24_accoa_c: - DG_LMPD_24_accoa_c +1000 BU_LMPD_24_accoa_c < 1000
SU_LMPD_24_accoa_c: FU_LMPD_24_accoa_c + BU_LMPD_24_accoa_c < 1
UF_LMPD_24_accoa_c: F_LMPD_24_accoa_c -1000 FU_LMPD_24_accoa_c < 0
UR_LMPD_24_accoa_c: R_LMPD_24_accoa_c -1000 BU_LMPD_24_accoa_c < 0
G_LMPD_25_accoa_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +9.4734 LC_adp_c +1.1842 LC_akg_c +2.3684 LC_amp_c -1.7763 LC_asp-L_c -11.8418 LC_atp_c +1.1842 LC_co2_c +1.1842 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -1.1842 LC_nh4_c +8.2893 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_25_accoa_c + DGo_LMPD_25_accoa_c = 0
FU_LMPD_25_accoa_c: DG_LMPD_25_accoa_c +1000 FU_LMPD_25_accoa_c < 1000
BU_LMPD_25_accoa_c: - DG_LMPD_25_accoa_c +1000 BU_LMPD_25_accoa_c < 1000
SU_LMPD_25_accoa_c: FU_LMPD_25_accoa_c + BU_LMPD_25_accoa_c < 1
UF_LMPD_25_accoa_c: F_LMPD_25_accoa_c -1000 FU_LMPD_25_accoa_c < 0
UR_LMPD_25_accoa_c: R_LMPD_25_accoa_c -1000 BU_LMPD_25_accoa_c < 0
SU_LMPD_26_adocbl_c: FU_LMPD_26_adocbl_c + BU_LMPD_26_adocbl_c < 1
UF_LMPD_26_adocbl_c: F_LMPD_26_adocbl_c -1000 FU_LMPD_26_adocbl_c < 0
UR_LMPD_26_adocbl_c: R_LMPD_26_adocbl_c -1000 BU_LMPD_26_adocbl_c < 0
SU_LMPD_27_adocbl_c: FU_LMPD_27_adocbl_c + BU_LMPD_27_adocbl_c < 1
UF_LMPD_27_adocbl_c: F_LMPD_27_adocbl_c -1000 FU_LMPD_27_adocbl_c < 0
UR_LMPD_27_adocbl_c: R_LMPD_27_adocbl_c -1000 BU_LMPD_27_adocbl_c < 0
SU_LMPD_28_adocbl_c: FU_LMPD_28_adocbl_c + BU_LMPD_28_adocbl_c < 1
UF_LMPD_28_adocbl_c: F_LMPD_28_adocbl_c -1000 FU_LMPD_28_adocbl_c < 0
UR_LMPD_28_adocbl_c: R_LMPD_28_adocbl_c -1000 BU_LMPD_28_adocbl_c < 0
SU_LMPD_29_adocbl_c: FU_LMPD_29_adocbl_c + BU_LMPD_29_adocbl_c < 1
UF_LMPD_29_adocbl_c: F_LMPD_29_adocbl_c -1000 FU_LMPD_29_adocbl_c < 0
UR_LMPD_29_adocbl_c: R_LMPD_29_adocbl_c -1000 BU_LMPD_29_adocbl_c < 0
G_LMPD_30_ala-L_c: +0.59209 LC_akg_c +0.59209 LC_ala-L_c -0.59209 LC_glu-L_c -0.59209 LC_pyr_c - DG_LMPD_30_ala-L_c + DGo_LMPD_30_ala-L_c = 0
FU_LMPD_30_ala-L_c: DG_LMPD_30_ala-L_c +1000 FU_LMPD_30_ala-L_c < 1000
BU_LMPD_30_ala-L_c: - DG_LMPD_30_ala-L_c +1000 BU_LMPD_30_ala-L_c < 1000
SU_LMPD_30_ala-L_c: FU_LMPD_30_ala-L_c + BU_LMPD_30_ala-L_c < 1
UF_LMPD_30_ala-L_c: F_LMPD_30_ala-L_c -1000 FU_LMPD_30_ala-L_c < 0
UR_LMPD_30_ala-L_c: R_LMPD_30_ala-L_c -1000 BU_LMPD_30_ala-L_c < 0
G_LMPD_31_amet_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +5.3288 LC_adp_c +1.1842 LC_akg_c +0.59209 LC_amet_c +2.9604 LC_amp_c -1.7763 LC_asp-L_c -8.2893 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -0.59209 LC_nh4_c +9.4734 LC_pi_c +1.7763 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_31_amet_c + DGo_LMPD_31_amet_c = 0
FU_LMPD_31_amet_c: DG_LMPD_31_amet_c +1000 FU_LMPD_31_amet_c < 1000
BU_LMPD_31_amet_c: - DG_LMPD_31_amet_c +1000 BU_LMPD_31_amet_c < 1000
SU_LMPD_31_amet_c: FU_LMPD_31_amet_c + BU_LMPD_31_amet_c < 1
UF_LMPD_31_amet_c: F_LMPD_31_amet_c -1000 FU_LMPD_31_amet_c < 0
UR_LMPD_31_amet_c: R_LMPD_31_amet_c -1000 BU_LMPD_31_amet_c < 0
G_LMPD_32_amet_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +7.6972 LC_adp_c +1.1842 LC_akg_c +0.59209 LC_amet_c +1.7763 LC_amp_c -1.7763 LC_asp-L_c -9.4734 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -0.59209 LC_nh4_c +9.4734 LC_pi_c +1.7763 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_32_amet_c + DGo_LMPD_32_amet_c = 0
FU_LMPD_32_amet_c: DG_LMPD_32_amet_c +1000 FU_LMPD_32_amet_c < 1000
BU_LMPD_32_amet_c: - DG_LMPD_32_amet_c +1000 BU_LMPD_32_amet_c < 1000
SU_LMPD_32_amet_c: FU_LMPD_32_amet_c + BU_LMPD_32_amet_c < 1
UF_LMPD_32_amet_c: F_LMPD_32_amet_c -1000 FU_LMPD_32_amet_c < 0
UR_LMPD_32_amet_c: R_LMPD_32_amet_c -1000 BU_LMPD_32_amet_c < 0
G_LMPD_33_arg-L_c: +0.59209 LC_ac_c -0.59209 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c +0.59209 LC_arg-L_c -0.59209 LC_asp-L_c -1.7763 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_coa_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -0.59209 LC_nh4_c +1.1842 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_33_arg-L_c + DGo_LMPD_33_arg-L_c = 0
FU_LMPD_33_arg-L_c: DG_LMPD_33_arg-L_c +1000 FU_LMPD_33_arg-L_c < 1000
BU_LMPD_33_arg-L_c: - DG_LMPD_33_arg-L_c +1000 BU_LMPD_33_arg-L_c < 1000
SU_LMPD_33_arg-L_c: FU_LMPD_33_arg-L_c + BU_LMPD_33_arg-L_c < 1
UF_LMPD_33_arg-L_c: F_LMPD_33_arg-L_c -1000 FU_LMPD_33_arg-L_c < 0
UR_LMPD_33_arg-L_c: R_LMPD_33_arg-L_c -1000 BU_LMPD_33_arg-L_c < 0
G_LMPD_34_asn-L_c: +0.59209 LC_amp_c +0.59209 LC_asn-L_c -0.59209 LC_asp-L_c -0.59209 LC_atp_c -0.59209 LC_nh4_c +0.59209 LC_ppi_c - DG_LMPD_34_asn-L_c + DGo_LMPD_34_asn-L_c = 0
FU_LMPD_34_asn-L_c: DG_LMPD_34_asn-L_c +1000 FU_LMPD_34_asn-L_c < 1000
BU_LMPD_34_asn-L_c: - DG_LMPD_34_asn-L_c +1000 BU_LMPD_34_asn-L_c < 1000
SU_LMPD_34_asn-L_c: FU_LMPD_34_asn-L_c + BU_LMPD_34_asn-L_c < 1
UF_LMPD_34_asn-L_c: F_LMPD_34_asn-L_c -1000 FU_LMPD_34_asn-L_c < 0
UR_LMPD_34_asn-L_c: R_LMPD_34_asn-L_c -1000 BU_LMPD_34_asn-L_c < 0
SU_LMPD_35_bmocogdp_c: FU_LMPD_35_bmocogdp_c + BU_LMPD_35_bmocogdp_c < 1
UF_LMPD_35_bmocogdp_c: F_LMPD_35_bmocogdp_c -1000 FU_LMPD_35_bmocogdp_c < 0
UR_LMPD_35_bmocogdp_c: R_LMPD_35_bmocogdp_c -1000 BU_LMPD_35_bmocogdp_c < 0
SU_LMPD_36_bmocogdp_c: FU_LMPD_36_bmocogdp_c + BU_LMPD_36_bmocogdp_c < 1
UF_LMPD_36_bmocogdp_c: F_LMPD_36_bmocogdp_c -1000 FU_LMPD_36_bmocogdp_c < 0
UR_LMPD_36_bmocogdp_c: R_LMPD_36_bmocogdp_c -1000 BU_LMPD_36_bmocogdp_c < 0
SU_LMPD_37_bmocogdp_c: FU_LMPD_37_bmocogdp_c + BU_LMPD_37_bmocogdp_c < 1
UF_LMPD_37_bmocogdp_c: F_LMPD_37_bmocogdp_c -1000 FU_LMPD_37_bmocogdp_c < 0
UR_LMPD_37_bmocogdp_c: R_LMPD_37_bmocogdp_c -1000 BU_LMPD_37_bmocogdp_c < 0
SU_LMPD_38_bmocogdp_c: FU_LMPD_38_bmocogdp_c + BU_LMPD_38_bmocogdp_c < 1
UF_LMPD_38_bmocogdp_c: F_LMPD_38_bmocogdp_c -1000 FU_LMPD_38_bmocogdp_c < 0
UR_LMPD_38_bmocogdp_c: R_LMPD_38_bmocogdp_c -1000 BU_LMPD_38_bmocogdp_c < 0
G_LMPD_39_btn_c: -2.3684 LC_3pg_c +0.59209 LC_5drib_c +1.1842 LC_ac_c -2.9604 LC_accoa_c +7.6972 LC_adp_c +2.3684 LC_akg_c +0.59209 LC_amob_c +6.513 LC_amp_c -1.7763 LC_asp-L_c -14.2101 LC_atp_c +0.59209 LC_btn_c -1.1842 LC_co2_c +3.5525 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -2.3684 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_mththf_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +7.1051 LC_nadp_c -7.1051 LC_nadph_c -0.59209 LC_nh4_c +15.3943 LC_pi_c +4.7367 LC_ppi_c +0.59209 LC_pyr_c -1.7763 LC_r5p_c -1.1842 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +1.1842 LC_trdox_c -1.1842 LC_trdrd_c +0.59209 LC_meoh_p - DG_LMPD_39_btn_c + DGo_LMPD_39_btn_c = 0
FU_LMPD_39_btn_c: DG_LMPD_39_btn_c +1000 FU_LMPD_39_btn_c < 1000
BU_LMPD_39_btn_c: - DG_LMPD_39_btn_c +1000 BU_LMPD_39_btn_c < 1000
SU_LMPD_39_btn_c: FU_LMPD_39_btn_c + BU_LMPD_39_btn_c < 1
UF_LMPD_39_btn_c: F_LMPD_39_btn_c -1000 FU_LMPD_39_btn_c < 0
UR_LMPD_39_btn_c: R_LMPD_39_btn_c -1000 BU_LMPD_39_btn_c < 0
G_LMPD_40_btn_c: -2.3684 LC_3pg_c +0.59209 LC_5drib_c +1.1842 LC_ac_c -2.9604 LC_accoa_c +8.2893 LC_adp_c +2.3684 LC_akg_c +0.59209 LC_amob_c +6.513 LC_amp_c -1.7763 LC_asp-L_c -14.8022 LC_atp_c +0.59209 LC_btn_c -1.1842 LC_co2_c +3.5525 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -2.3684 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_mththf_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +7.1051 LC_nadp_c -7.1051 LC_nadph_c -0.59209 LC_nh4_c +15.9864 LC_pi_c +4.7367 LC_ppi_c +0.59209 LC_pyr_c -1.7763 LC_r5p_c -1.1842 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +1.1842 LC_trdox_c -1.1842 LC_trdrd_c +0.59209 LC_meoh_p - DG_LMPD_40_btn_c + DGo_LMPD_40_btn_c = 0
FU_LMPD_40_btn_c: DG_LMPD_40_btn_c +1000 FU_LMPD_40_btn_c < 1000
BU_LMPD_40_btn_c: - DG_LMPD_40_btn_c +1000 BU_LMPD_40_btn_c < 1000
SU_LMPD_40_btn_c: FU_LMPD_40_btn_c + BU_LMPD_40_btn_c < 1
UF_LMPD_40_btn_c: F_LMPD_40_btn_c -1000 FU_LMPD_40_btn_c < 0
UR_LMPD_40_btn_c: R_LMPD_40_btn_c -1000 BU_LMPD_40_btn_c < 0
G_LMPD_41_btn_c: -2.3684 LC_3pg_c +0.59209 LC_5drib_c +1.1842 LC_ac_c -2.9604 LC_accoa_c +11.2497 LC_adp_c +2.3684 LC_akg_c +0.59209 LC_amob_c +4.7367 LC_amp_c -1.7763 LC_asp-L_c -15.9864 LC_atp_c +0.59209 LC_btn_c -1.1842 LC_co2_c +3.5525 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -2.3684 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_mththf_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +7.1051 LC_nadp_c -7.1051 LC_nadph_c -0.59209 LC_nh4_c +15.3943 LC_pi_c +4.7367 LC_ppi_c +0.59209 LC_pyr_c -1.7763 LC_r5p_c -1.1842 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +1.1842 LC_trdox_c -1.1842 LC_trdrd_c +0.59209 LC_meoh_p - DG_LMPD_41_btn_c + DGo_LMPD_41_btn_c = 0
FU_LMPD_41_btn_c: DG_LMPD_41_btn_c +1000 FU_LMPD_41_btn_c < 1000
BU_LMPD_41_btn_c: - DG_LMPD_41_btn_c +1000 BU_LMPD_41_btn_c < 1000
SU_LMPD_41_btn_c: FU_LMPD_41_btn_c + BU_LMPD_41_btn_c < 1
UF_LMPD_41_btn_c: F_LMPD_41_btn_c -1000 FU_LMPD_41_btn_c < 0
UR_LMPD_41_btn_c: R_LMPD_41_btn_c -1000 BU_LMPD_41_btn_c < 0
G_LMPD_42_btn_c: -2.3684 LC_3pg_c +0.59209 LC_5drib_c +1.1842 LC_ac_c -2.9604 LC_accoa_c +11.8418 LC_adp_c +2.3684 LC_akg_c +0.59209 LC_amob_c +4.7367 LC_amp_c -1.7763 LC_asp-L_c -16.5785 LC_atp_c +0.59209 LC_btn_c -1.1842 LC_co2_c +3.5525 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -2.3684 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_mththf_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +7.1051 LC_nadp_c -7.1051 LC_nadph_c -0.59209 LC_nh4_c +15.9864 LC_pi_c +4.7367 LC_ppi_c +0.59209 LC_pyr_c -1.7763 LC_r5p_c -1.1842 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +1.1842 LC_trdox_c -1.1842 LC_trdrd_c +0.59209 LC_meoh_p - DG_LMPD_42_btn_c + DGo_LMPD_42_btn_c = 0
FU_LMPD_42_btn_c: DG_LMPD_42_btn_c +1000 FU_LMPD_42_btn_c < 1000
BU_LMPD_42_btn_c: - DG_LMPD_42_btn_c +1000 BU_LMPD_42_btn_c < 1000
SU_LMPD_42_btn_c: FU_LMPD_42_btn_c + BU_LMPD_42_btn_c < 1
UF_LMPD_42_btn_c: F_LMPD_42_btn_c -1000 FU_LMPD_42_btn_c < 0
UR_LMPD_42_btn_c: R_LMPD_42_btn_c -1000 BU_LMPD_42_btn_c < 0
G_LMPD_43_btn_c: -3.5525 LC_3pg_c +0.59209 LC_5drib_c +1.7763 LC_ac_c -3.5525 LC_accoa_c +17.7627 LC_adp_c +3.5525 LC_akg_c +0.59209 LC_amob_c +8.8813 LC_amp_c -4.7367 LC_asp-L_c -26.644 LC_atp_c +0.59209 LC_btn_c -2.3684 LC_co2_c +4.7367 LC_coa_c -3.5525 LC_for_c +3.5525 LC_fum_c -3.5525 LC_glu-L_c +0.59209 LC_met-L_c +0.59209 LC_mththf_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +10.0655 LC_nadp_c -10.0655 LC_nadph_c -2.3684 LC_nh4_c +30.1966 LC_pi_c +5.3288 LC_ppi_c +1.1842 LC_pyr_c -1.7763 LC_r5p_c -1.7763 LC_so4_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +1.7763 LC_trdox_c -1.7763 LC_trdrd_c +0.59209 LC_meoh_p +1.1842 LC_ade_p - DG_LMPD_43_btn_c + DGo_LMPD_43_btn_c = 0
FU_LMPD_43_btn_c: DG_LMPD_43_btn_c +1000 FU_LMPD_43_btn_c < 1000
BU_LMPD_43_btn_c: - DG_LMPD_43_btn_c +1000 BU_LMPD_43_btn_c < 1000
SU_LMPD_43_btn_c: FU_LMPD_43_btn_c + BU_LMPD_43_btn_c < 1
UF_LMPD_43_btn_c: F_LMPD_43_btn_c -1000 FU_LMPD_43_btn_c < 0
UR_LMPD_43_btn_c: R_LMPD_43_btn_c -1000 BU_LMPD_43_btn_c < 0
G_LMPD_44_btn_c: -3.5525 LC_3pg_c +0.59209 LC_5drib_c +1.7763 LC_ac_c -3.5525 LC_accoa_c +18.3548 LC_adp_c +3.5525 LC_akg_c +0.59209 LC_amob_c +8.8813 LC_amp_c -4.7367 LC_asp-L_c -27.2361 LC_atp_c +0.59209 LC_btn_c -2.3684 LC_co2_c +4.7367 LC_coa_c -3.5525 LC_for_c +3.5525 LC_fum_c -3.5525 LC_glu-L_c +0.59209 LC_met-L_c +0.59209 LC_mththf_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +10.0655 LC_nadp_c -10.0655 LC_nadph_c -2.3684 LC_nh4_c +30.7886 LC_pi_c +5.3288 LC_ppi_c +1.1842 LC_pyr_c -1.7763 LC_r5p_c -1.7763 LC_so4_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +1.7763 LC_trdox_c -1.7763 LC_trdrd_c +0.59209 LC_meoh_p +1.1842 LC_ade_p - DG_LMPD_44_btn_c + DGo_LMPD_44_btn_c = 0
FU_LMPD_44_btn_c: DG_LMPD_44_btn_c +1000 FU_LMPD_44_btn_c < 1000
BU_LMPD_44_btn_c: - DG_LMPD_44_btn_c +1000 BU_LMPD_44_btn_c < 1000
SU_LMPD_44_btn_c: FU_LMPD_44_btn_c + BU_LMPD_44_btn_c < 1
UF_LMPD_44_btn_c: F_LMPD_44_btn_c -1000 FU_LMPD_44_btn_c < 0
UR_LMPD_44_btn_c: R_LMPD_44_btn_c -1000 BU_LMPD_44_btn_c < 0
G_LMPD_45_btn_c: -3.5525 LC_3pg_c +0.59209 LC_5drib_c +1.7763 LC_ac_c -3.5525 LC_accoa_c +24.8678 LC_adp_c +3.5525 LC_akg_c +0.59209 LC_amob_c +5.3288 LC_amp_c -4.7367 LC_asp-L_c -30.1966 LC_atp_c +0.59209 LC_btn_c -2.3684 LC_co2_c +4.7367 LC_coa_c -3.5525 LC_for_c +3.5525 LC_fum_c -3.5525 LC_glu-L_c +0.59209 LC_met-L_c +0.59209 LC_mththf_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +10.0655 LC_nadp_c -10.0655 LC_nadph_c -2.3684 LC_nh4_c +30.1966 LC_pi_c +5.3288 LC_ppi_c +1.1842 LC_pyr_c -1.7763 LC_r5p_c -1.7763 LC_so4_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +1.7763 LC_trdox_c -1.7763 LC_trdrd_c +0.59209 LC_meoh_p +1.1842 LC_ade_p - DG_LMPD_45_btn_c + DGo_LMPD_45_btn_c = 0
FU_LMPD_45_btn_c: DG_LMPD_45_btn_c +1000 FU_LMPD_45_btn_c < 1000
BU_LMPD_45_btn_c: - DG_LMPD_45_btn_c +1000 BU_LMPD_45_btn_c < 1000
SU_LMPD_45_btn_c: FU_LMPD_45_btn_c + BU_LMPD_45_btn_c < 1
UF_LMPD_45_btn_c: F_LMPD_45_btn_c -1000 FU_LMPD_45_btn_c < 0
UR_LMPD_45_btn_c: R_LMPD_45_btn_c -1000 BU_LMPD_45_btn_c < 0
G_LMPD_47_btn_c: -3.5525 LC_3pg_c +0.59209 LC_5drib_c +1.7763 LC_ac_c -3.5525 LC_accoa_c +25.4598 LC_adp_c +3.5525 LC_akg_c +0.59209 LC_amob_c +5.3288 LC_amp_c -4.7367 LC_asp-L_c -30.7886 LC_atp_c +0.59209 LC_btn_c -2.3684 LC_co2_c +4.7367 LC_coa_c -3.5525 LC_for_c +3.5525 LC_fum_c -3.5525 LC_glu-L_c +0.59209 LC_met-L_c +0.59209 LC_mththf_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +10.0655 LC_nadp_c -10.0655 LC_nadph_c -2.3684 LC_nh4_c +30.7886 LC_pi_c +5.3288 LC_ppi_c +1.1842 LC_pyr_c -1.7763 LC_r5p_c -1.7763 LC_so4_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +1.7763 LC_trdox_c -1.7763 LC_trdrd_c +0.59209 LC_meoh_p +1.1842 LC_ade_p - DG_LMPD_47_btn_c + DGo_LMPD_47_btn_c = 0
FU_LMPD_47_btn_c: DG_LMPD_47_btn_c +1000 FU_LMPD_47_btn_c < 1000
BU_LMPD_47_btn_c: - DG_LMPD_47_btn_c +1000 BU_LMPD_47_btn_c < 1000
SU_LMPD_47_btn_c: FU_LMPD_47_btn_c + BU_LMPD_47_btn_c < 1
UF_LMPD_47_btn_c: F_LMPD_47_btn_c -1000 FU_LMPD_47_btn_c < 0
UR_LMPD_47_btn_c: R_LMPD_47_btn_c -1000 BU_LMPD_47_btn_c < 0
G_LMPD_48_ca2_c: -0.59209 LC_3pg_c -0.59209 LC_accoa_c +1.7763 LC_adp_c +0.59209 LC_akg_c -1.7763 LC_atp_c +0.59209 LC_coa_c -1.1842 LC_glu-L_c +0.59209 LC_gthrd_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +2.3684 LC_pi_c +0.59209 LC_acser_p -0.59209 LC_gthrd_p - DG_LMPD_48_ca2_c + DGo_LMPD_48_ca2_c = 0
FU_LMPD_48_ca2_c: DG_LMPD_48_ca2_c +1000 FU_LMPD_48_ca2_c < 1000
BU_LMPD_48_ca2_c: - DG_LMPD_48_ca2_c +1000 BU_LMPD_48_ca2_c < 1000
SU_LMPD_48_ca2_c: FU_LMPD_48_ca2_c + BU_LMPD_48_ca2_c < 1
UF_LMPD_48_ca2_c: F_LMPD_48_ca2_c -1000 FU_LMPD_48_ca2_c < 0
UR_LMPD_48_ca2_c: R_LMPD_48_ca2_c -1000 BU_LMPD_48_ca2_c < 0
G_LMPD_49_ca2_c: -0.59209 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +1.7763 LC_adp_c +0.59209 LC_akg_c -1.7763 LC_atp_c +0.59209 LC_coa_c -1.1842 LC_glu-L_c +0.59209 LC_gthrd_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nh4_c +2.3684 LC_pi_c +0.59209 LC_pyr_c -0.59209 LC_gthrd_p - DG_LMPD_49_ca2_c + DGo_LMPD_49_ca2_c = 0
FU_LMPD_49_ca2_c: DG_LMPD_49_ca2_c +1000 FU_LMPD_49_ca2_c < 1000
BU_LMPD_49_ca2_c: - DG_LMPD_49_ca2_c +1000 BU_LMPD_49_ca2_c < 1000
SU_LMPD_49_ca2_c: FU_LMPD_49_ca2_c + BU_LMPD_49_ca2_c < 1
UF_LMPD_49_ca2_c: F_LMPD_49_ca2_c -1000 FU_LMPD_49_ca2_c < 0
UR_LMPD_49_ca2_c: R_LMPD_49_ca2_c -1000 BU_LMPD_49_ca2_c < 0
G_LMPD_50_ca2_c: +1.7763 LC_adp_c -1.7763 LC_atp_c -0.59209 LC_glu-L_c +0.59209 LC_gthrd_c +1.7763 LC_pi_c +0.59209 LC_glu-L_p -0.59209 LC_gthrd_p - DG_LMPD_50_ca2_c + DGo_LMPD_50_ca2_c = 0
FU_LMPD_50_ca2_c: DG_LMPD_50_ca2_c +1000 FU_LMPD_50_ca2_c < 1000
BU_LMPD_50_ca2_c: - DG_LMPD_50_ca2_c +1000 BU_LMPD_50_ca2_c < 1000
SU_LMPD_50_ca2_c: FU_LMPD_50_ca2_c + BU_LMPD_50_ca2_c < 1
UF_LMPD_50_ca2_c: F_LMPD_50_ca2_c -1000 FU_LMPD_50_ca2_c < 0
UR_LMPD_50_ca2_c: R_LMPD_50_ca2_c -1000 BU_LMPD_50_ca2_c < 0
G_LMPD_51_chor_c: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_chor_c -0.59209 LC_e4p_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_pep_c +2.3684 LC_pi_c - DG_LMPD_51_chor_c + DGo_LMPD_51_chor_c = 0
FU_LMPD_51_chor_c: DG_LMPD_51_chor_c +1000 FU_LMPD_51_chor_c < 1000
BU_LMPD_51_chor_c: - DG_LMPD_51_chor_c +1000 BU_LMPD_51_chor_c < 1000
SU_LMPD_51_chor_c: FU_LMPD_51_chor_c + BU_LMPD_51_chor_c < 1
UF_LMPD_51_chor_c: F_LMPD_51_chor_c -1000 FU_LMPD_51_chor_c < 0
UR_LMPD_51_chor_c: R_LMPD_51_chor_c -1000 BU_LMPD_51_chor_c < 0
G_LMPD_52_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +7.1051 LC_nad_c -7.1051 LC_nadh_c +26.0519 LC_nadp_c -26.0519 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_52_clpn160_p + DGo_LMPD_52_clpn160_p = 0
FU_LMPD_52_clpn160_p: DG_LMPD_52_clpn160_p +1000 FU_LMPD_52_clpn160_p < 1000
BU_LMPD_52_clpn160_p: - DG_LMPD_52_clpn160_p +1000 BU_LMPD_52_clpn160_p < 1000
SU_LMPD_52_clpn160_p: FU_LMPD_52_clpn160_p + BU_LMPD_52_clpn160_p < 1
UF_LMPD_52_clpn160_p: F_LMPD_52_clpn160_p -1000 FU_LMPD_52_clpn160_p < 0
UR_LMPD_52_clpn160_p: R_LMPD_52_clpn160_p -1000 BU_LMPD_52_clpn160_p < 0
G_LMPD_53_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +9.4734 LC_nad_c -9.4734 LC_nadh_c +23.6836 LC_nadp_c -23.6836 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_53_clpn160_p + DGo_LMPD_53_clpn160_p = 0
FU_LMPD_53_clpn160_p: DG_LMPD_53_clpn160_p +1000 FU_LMPD_53_clpn160_p < 1000
BU_LMPD_53_clpn160_p: - DG_LMPD_53_clpn160_p +1000 BU_LMPD_53_clpn160_p < 1000
SU_LMPD_53_clpn160_p: FU_LMPD_53_clpn160_p + BU_LMPD_53_clpn160_p < 1
UF_LMPD_53_clpn160_p: F_LMPD_53_clpn160_p -1000 FU_LMPD_53_clpn160_p < 0
UR_LMPD_53_clpn160_p: R_LMPD_53_clpn160_p -1000 BU_LMPD_53_clpn160_p < 0
G_LMPD_54_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +11.8418 LC_nad_c -11.8418 LC_nadh_c +21.3152 LC_nadp_c -21.3152 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_54_clpn160_p + DGo_LMPD_54_clpn160_p = 0
FU_LMPD_54_clpn160_p: DG_LMPD_54_clpn160_p +1000 FU_LMPD_54_clpn160_p < 1000
BU_LMPD_54_clpn160_p: - DG_LMPD_54_clpn160_p +1000 BU_LMPD_54_clpn160_p < 1000
SU_LMPD_54_clpn160_p: FU_LMPD_54_clpn160_p + BU_LMPD_54_clpn160_p < 1
UF_LMPD_54_clpn160_p: F_LMPD_54_clpn160_p -1000 FU_LMPD_54_clpn160_p < 0
UR_LMPD_54_clpn160_p: R_LMPD_54_clpn160_p -1000 BU_LMPD_54_clpn160_p < 0
G_LMPD_55_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +14.2101 LC_nad_c -14.2101 LC_nadh_c +18.9469 LC_nadp_c -18.9469 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_55_clpn160_p + DGo_LMPD_55_clpn160_p = 0
FU_LMPD_55_clpn160_p: DG_LMPD_55_clpn160_p +1000 FU_LMPD_55_clpn160_p < 1000
BU_LMPD_55_clpn160_p: - DG_LMPD_55_clpn160_p +1000 BU_LMPD_55_clpn160_p < 1000
SU_LMPD_55_clpn160_p: FU_LMPD_55_clpn160_p + BU_LMPD_55_clpn160_p < 1
UF_LMPD_55_clpn160_p: F_LMPD_55_clpn160_p -1000 FU_LMPD_55_clpn160_p < 0
UR_LMPD_55_clpn160_p: R_LMPD_55_clpn160_p -1000 BU_LMPD_55_clpn160_p < 0
G_LMPD_56_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +16.5785 LC_nad_c -16.5785 LC_nadh_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_56_clpn160_p + DGo_LMPD_56_clpn160_p = 0
FU_LMPD_56_clpn160_p: DG_LMPD_56_clpn160_p +1000 FU_LMPD_56_clpn160_p < 1000
BU_LMPD_56_clpn160_p: - DG_LMPD_56_clpn160_p +1000 BU_LMPD_56_clpn160_p < 1000
SU_LMPD_56_clpn160_p: FU_LMPD_56_clpn160_p + BU_LMPD_56_clpn160_p < 1
UF_LMPD_56_clpn160_p: F_LMPD_56_clpn160_p -1000 FU_LMPD_56_clpn160_p < 0
UR_LMPD_56_clpn160_p: R_LMPD_56_clpn160_p -1000 BU_LMPD_56_clpn160_p < 0
G_LMPD_57_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +2.3684 LC_nad_c -2.3684 LC_nadh_c +30.7886 LC_nadp_c -30.7886 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_57_clpn160_p + DGo_LMPD_57_clpn160_p = 0
FU_LMPD_57_clpn160_p: DG_LMPD_57_clpn160_p +1000 FU_LMPD_57_clpn160_p < 1000
BU_LMPD_57_clpn160_p: - DG_LMPD_57_clpn160_p +1000 BU_LMPD_57_clpn160_p < 1000
SU_LMPD_57_clpn160_p: FU_LMPD_57_clpn160_p + BU_LMPD_57_clpn160_p < 1
UF_LMPD_57_clpn160_p: F_LMPD_57_clpn160_p -1000 FU_LMPD_57_clpn160_p < 0
UR_LMPD_57_clpn160_p: R_LMPD_57_clpn160_p -1000 BU_LMPD_57_clpn160_p < 0
G_LMPD_58_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +33.157 LC_nadp_c -33.157 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_58_clpn160_p + DGo_LMPD_58_clpn160_p = 0
FU_LMPD_58_clpn160_p: DG_LMPD_58_clpn160_p +1000 FU_LMPD_58_clpn160_p < 1000
BU_LMPD_58_clpn160_p: - DG_LMPD_58_clpn160_p +1000 BU_LMPD_58_clpn160_p < 1000
SU_LMPD_58_clpn160_p: FU_LMPD_58_clpn160_p + BU_LMPD_58_clpn160_p < 1
UF_LMPD_58_clpn160_p: F_LMPD_58_clpn160_p -1000 FU_LMPD_58_clpn160_p < 0
UR_LMPD_58_clpn160_p: R_LMPD_58_clpn160_p -1000 BU_LMPD_58_clpn160_p < 0
G_LMPD_59_clpn160_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +4.7367 LC_nad_c -4.7367 LC_nadh_c +28.4203 LC_nadp_c -28.4203 LC_nadph_c +1.1842 LC_pg160_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_59_clpn160_p + DGo_LMPD_59_clpn160_p = 0
FU_LMPD_59_clpn160_p: DG_LMPD_59_clpn160_p +1000 FU_LMPD_59_clpn160_p < 1000
BU_LMPD_59_clpn160_p: - DG_LMPD_59_clpn160_p +1000 BU_LMPD_59_clpn160_p < 1000
SU_LMPD_59_clpn160_p: FU_LMPD_59_clpn160_p + BU_LMPD_59_clpn160_p < 1
UF_LMPD_59_clpn160_p: F_LMPD_59_clpn160_p -1000 FU_LMPD_59_clpn160_p < 0
UR_LMPD_59_clpn160_p: R_LMPD_59_clpn160_p -1000 BU_LMPD_59_clpn160_p < 0
G_LMPD_60_clpn160_p: -18.9469 LC_accoa_c +2.3684 LC_adp_c +2.3684 LC_amp_c -4.7367 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +16.5785 LC_nad_c -16.5785 LC_nadh_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +1.1842 LC_pg160_c +1.1842 LC_pi_c +3.5525 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_60_clpn160_p + DGo_LMPD_60_clpn160_p = 0
FU_LMPD_60_clpn160_p: DG_LMPD_60_clpn160_p +1000 FU_LMPD_60_clpn160_p < 1000
BU_LMPD_60_clpn160_p: - DG_LMPD_60_clpn160_p +1000 BU_LMPD_60_clpn160_p < 1000
SU_LMPD_60_clpn160_p: FU_LMPD_60_clpn160_p + BU_LMPD_60_clpn160_p < 1
UF_LMPD_60_clpn160_p: F_LMPD_60_clpn160_p -1000 FU_LMPD_60_clpn160_p < 0
UR_LMPD_60_clpn160_p: R_LMPD_60_clpn160_p -1000 BU_LMPD_60_clpn160_p < 0
G_LMPD_61_clpn160_p: -18.9469 LC_accoa_c +2.3684 LC_adp_c +2.3684 LC_amp_c -4.7367 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +33.157 LC_nad_c -33.157 LC_nadh_c +1.1842 LC_pg160_c +1.1842 LC_pi_c +3.5525 LC_ppi_c +0.59209 LC_clpn160_p +0.59209 LC_glyc_p -1.1842 LC_pg160_p - DG_LMPD_61_clpn160_p + DGo_LMPD_61_clpn160_p = 0
FU_LMPD_61_clpn160_p: DG_LMPD_61_clpn160_p +1000 FU_LMPD_61_clpn160_p < 1000
BU_LMPD_61_clpn160_p: - DG_LMPD_61_clpn160_p +1000 BU_LMPD_61_clpn160_p < 1000
SU_LMPD_61_clpn160_p: FU_LMPD_61_clpn160_p + BU_LMPD_61_clpn160_p < 1
UF_LMPD_61_clpn160_p: F_LMPD_61_clpn160_p -1000 FU_LMPD_61_clpn160_p < 0
UR_LMPD_61_clpn160_p: R_LMPD_61_clpn160_p -1000 BU_LMPD_61_clpn160_p < 0
G_LMPD_62_clpn161_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +7.1051 LC_nad_c -7.1051 LC_nadh_c +23.6836 LC_nadp_c -23.6836 LC_nadph_c +1.1842 LC_pg161_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn161_p +0.59209 LC_glyc_p -1.1842 LC_pg161_p - DG_LMPD_62_clpn161_p + DGo_LMPD_62_clpn161_p = 0
FU_LMPD_62_clpn161_p: DG_LMPD_62_clpn161_p +1000 FU_LMPD_62_clpn161_p < 1000
BU_LMPD_62_clpn161_p: - DG_LMPD_62_clpn161_p +1000 BU_LMPD_62_clpn161_p < 1000
SU_LMPD_62_clpn161_p: FU_LMPD_62_clpn161_p + BU_LMPD_62_clpn161_p < 1
UF_LMPD_62_clpn161_p: F_LMPD_62_clpn161_p -1000 FU_LMPD_62_clpn161_p < 0
UR_LMPD_62_clpn161_p: R_LMPD_62_clpn161_p -1000 BU_LMPD_62_clpn161_p < 0
G_LMPD_63_clpn161_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +9.4734 LC_nad_c -9.4734 LC_nadh_c +21.3152 LC_nadp_c -21.3152 LC_nadph_c +1.1842 LC_pg161_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn161_p +0.59209 LC_glyc_p -1.1842 LC_pg161_p - DG_LMPD_63_clpn161_p + DGo_LMPD_63_clpn161_p = 0
FU_LMPD_63_clpn161_p: DG_LMPD_63_clpn161_p +1000 FU_LMPD_63_clpn161_p < 1000
BU_LMPD_63_clpn161_p: - DG_LMPD_63_clpn161_p +1000 BU_LMPD_63_clpn161_p < 1000
SU_LMPD_63_clpn161_p: FU_LMPD_63_clpn161_p + BU_LMPD_63_clpn161_p < 1
UF_LMPD_63_clpn161_p: F_LMPD_63_clpn161_p -1000 FU_LMPD_63_clpn161_p < 0
UR_LMPD_63_clpn161_p: R_LMPD_63_clpn161_p -1000 BU_LMPD_63_clpn161_p < 0
G_LMPD_64_clpn161_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +11.8418 LC_nad_c -11.8418 LC_nadh_c +18.9469 LC_nadp_c -18.9469 LC_nadph_c +1.1842 LC_pg161_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn161_p +0.59209 LC_glyc_p -1.1842 LC_pg161_p - DG_LMPD_64_clpn161_p + DGo_LMPD_64_clpn161_p = 0
FU_LMPD_64_clpn161_p: DG_LMPD_64_clpn161_p +1000 FU_LMPD_64_clpn161_p < 1000
BU_LMPD_64_clpn161_p: - DG_LMPD_64_clpn161_p +1000 BU_LMPD_64_clpn161_p < 1000
SU_LMPD_64_clpn161_p: FU_LMPD_64_clpn161_p + BU_LMPD_64_clpn161_p < 1
UF_LMPD_64_clpn161_p: F_LMPD_64_clpn161_p -1000 FU_LMPD_64_clpn161_p < 0
UR_LMPD_64_clpn161_p: R_LMPD_64_clpn161_p -1000 BU_LMPD_64_clpn161_p < 0
G_LMPD_65_clpn161_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +14.2101 LC_nad_c -14.2101 LC_nadh_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +1.1842 LC_pg161_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn161_p +0.59209 LC_glyc_p -1.1842 LC_pg161_p - DG_LMPD_65_clpn161_p + DGo_LMPD_65_clpn161_p = 0
FU_LMPD_65_clpn161_p: DG_LMPD_65_clpn161_p +1000 FU_LMPD_65_clpn161_p < 1000
BU_LMPD_65_clpn161_p: - DG_LMPD_65_clpn161_p +1000 BU_LMPD_65_clpn161_p < 1000
SU_LMPD_65_clpn161_p: FU_LMPD_65_clpn161_p + BU_LMPD_65_clpn161_p < 1
UF_LMPD_65_clpn161_p: F_LMPD_65_clpn161_p -1000 FU_LMPD_65_clpn161_p < 0
UR_LMPD_65_clpn161_p: R_LMPD_65_clpn161_p -1000 BU_LMPD_65_clpn161_p < 0
G_LMPD_66_clpn161_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +2.3684 LC_nad_c -2.3684 LC_nadh_c +28.4203 LC_nadp_c -28.4203 LC_nadph_c +1.1842 LC_pg161_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn161_p +0.59209 LC_glyc_p -1.1842 LC_pg161_p - DG_LMPD_66_clpn161_p + DGo_LMPD_66_clpn161_p = 0
FU_LMPD_66_clpn161_p: DG_LMPD_66_clpn161_p +1000 FU_LMPD_66_clpn161_p < 1000
BU_LMPD_66_clpn161_p: - DG_LMPD_66_clpn161_p +1000 BU_LMPD_66_clpn161_p < 1000
SU_LMPD_66_clpn161_p: FU_LMPD_66_clpn161_p + BU_LMPD_66_clpn161_p < 1
UF_LMPD_66_clpn161_p: F_LMPD_66_clpn161_p -1000 FU_LMPD_66_clpn161_p < 0
UR_LMPD_66_clpn161_p: R_LMPD_66_clpn161_p -1000 BU_LMPD_66_clpn161_p < 0
G_LMPD_67_clpn161_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +30.7886 LC_nadp_c -30.7886 LC_nadph_c +1.1842 LC_pg161_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn161_p +0.59209 LC_glyc_p -1.1842 LC_pg161_p - DG_LMPD_67_clpn161_p + DGo_LMPD_67_clpn161_p = 0
FU_LMPD_67_clpn161_p: DG_LMPD_67_clpn161_p +1000 FU_LMPD_67_clpn161_p < 1000
BU_LMPD_67_clpn161_p: - DG_LMPD_67_clpn161_p +1000 BU_LMPD_67_clpn161_p < 1000
SU_LMPD_67_clpn161_p: FU_LMPD_67_clpn161_p + BU_LMPD_67_clpn161_p < 1
UF_LMPD_67_clpn161_p: F_LMPD_67_clpn161_p -1000 FU_LMPD_67_clpn161_p < 0
UR_LMPD_67_clpn161_p: R_LMPD_67_clpn161_p -1000 BU_LMPD_67_clpn161_p < 0
G_LMPD_68_clpn161_p: -18.9469 LC_accoa_c +18.9469 LC_adp_c -18.9469 LC_atp_c +18.9469 LC_coa_c -2.3684 LC_glyc3p_c +4.7367 LC_nad_c -4.7367 LC_nadh_c +26.0519 LC_nadp_c -26.0519 LC_nadph_c +1.1842 LC_pg161_c +17.7627 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn161_p +0.59209 LC_glyc_p -1.1842 LC_pg161_p - DG_LMPD_68_clpn161_p + DGo_LMPD_68_clpn161_p = 0
FU_LMPD_68_clpn161_p: DG_LMPD_68_clpn161_p +1000 FU_LMPD_68_clpn161_p < 1000
BU_LMPD_68_clpn161_p: - DG_LMPD_68_clpn161_p +1000 BU_LMPD_68_clpn161_p < 1000
SU_LMPD_68_clpn161_p: FU_LMPD_68_clpn161_p + BU_LMPD_68_clpn161_p < 1
UF_LMPD_68_clpn161_p: F_LMPD_68_clpn161_p -1000 FU_LMPD_68_clpn161_p < 0
UR_LMPD_68_clpn161_p: R_LMPD_68_clpn161_p -1000 BU_LMPD_68_clpn161_p < 0
G_LMPD_69_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +7.1051 LC_nad_c -7.1051 LC_nadh_c +28.4203 LC_nadp_c -28.4203 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_69_clpn181_p + DGo_LMPD_69_clpn181_p = 0
FU_LMPD_69_clpn181_p: DG_LMPD_69_clpn181_p +1000 FU_LMPD_69_clpn181_p < 1000
BU_LMPD_69_clpn181_p: - DG_LMPD_69_clpn181_p +1000 BU_LMPD_69_clpn181_p < 1000
SU_LMPD_69_clpn181_p: FU_LMPD_69_clpn181_p + BU_LMPD_69_clpn181_p < 1
UF_LMPD_69_clpn181_p: F_LMPD_69_clpn181_p -1000 FU_LMPD_69_clpn181_p < 0
UR_LMPD_69_clpn181_p: R_LMPD_69_clpn181_p -1000 BU_LMPD_69_clpn181_p < 0
G_LMPD_70_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +9.4734 LC_nad_c -9.4734 LC_nadh_c +26.0519 LC_nadp_c -26.0519 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_70_clpn181_p + DGo_LMPD_70_clpn181_p = 0
FU_LMPD_70_clpn181_p: DG_LMPD_70_clpn181_p +1000 FU_LMPD_70_clpn181_p < 1000
BU_LMPD_70_clpn181_p: - DG_LMPD_70_clpn181_p +1000 BU_LMPD_70_clpn181_p < 1000
SU_LMPD_70_clpn181_p: FU_LMPD_70_clpn181_p + BU_LMPD_70_clpn181_p < 1
UF_LMPD_70_clpn181_p: F_LMPD_70_clpn181_p -1000 FU_LMPD_70_clpn181_p < 0
UR_LMPD_70_clpn181_p: R_LMPD_70_clpn181_p -1000 BU_LMPD_70_clpn181_p < 0
G_LMPD_71_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +11.8418 LC_nad_c -11.8418 LC_nadh_c +23.6836 LC_nadp_c -23.6836 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_71_clpn181_p + DGo_LMPD_71_clpn181_p = 0
FU_LMPD_71_clpn181_p: DG_LMPD_71_clpn181_p +1000 FU_LMPD_71_clpn181_p < 1000
BU_LMPD_71_clpn181_p: - DG_LMPD_71_clpn181_p +1000 BU_LMPD_71_clpn181_p < 1000
SU_LMPD_71_clpn181_p: FU_LMPD_71_clpn181_p + BU_LMPD_71_clpn181_p < 1
UF_LMPD_71_clpn181_p: F_LMPD_71_clpn181_p -1000 FU_LMPD_71_clpn181_p < 0
UR_LMPD_71_clpn181_p: R_LMPD_71_clpn181_p -1000 BU_LMPD_71_clpn181_p < 0
G_LMPD_72_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +14.2101 LC_nad_c -14.2101 LC_nadh_c +21.3152 LC_nadp_c -21.3152 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_72_clpn181_p + DGo_LMPD_72_clpn181_p = 0
FU_LMPD_72_clpn181_p: DG_LMPD_72_clpn181_p +1000 FU_LMPD_72_clpn181_p < 1000
BU_LMPD_72_clpn181_p: - DG_LMPD_72_clpn181_p +1000 BU_LMPD_72_clpn181_p < 1000
SU_LMPD_72_clpn181_p: FU_LMPD_72_clpn181_p + BU_LMPD_72_clpn181_p < 1
UF_LMPD_72_clpn181_p: F_LMPD_72_clpn181_p -1000 FU_LMPD_72_clpn181_p < 0
UR_LMPD_72_clpn181_p: R_LMPD_72_clpn181_p -1000 BU_LMPD_72_clpn181_p < 0
G_LMPD_73_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +16.5785 LC_nad_c -16.5785 LC_nadh_c +18.9469 LC_nadp_c -18.9469 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_73_clpn181_p + DGo_LMPD_73_clpn181_p = 0
FU_LMPD_73_clpn181_p: DG_LMPD_73_clpn181_p +1000 FU_LMPD_73_clpn181_p < 1000
BU_LMPD_73_clpn181_p: - DG_LMPD_73_clpn181_p +1000 BU_LMPD_73_clpn181_p < 1000
SU_LMPD_73_clpn181_p: FU_LMPD_73_clpn181_p + BU_LMPD_73_clpn181_p < 1
UF_LMPD_73_clpn181_p: F_LMPD_73_clpn181_p -1000 FU_LMPD_73_clpn181_p < 0
UR_LMPD_73_clpn181_p: R_LMPD_73_clpn181_p -1000 BU_LMPD_73_clpn181_p < 0
G_LMPD_74_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +2.3684 LC_nad_c -2.3684 LC_nadh_c +33.157 LC_nadp_c -33.157 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_74_clpn181_p + DGo_LMPD_74_clpn181_p = 0
FU_LMPD_74_clpn181_p: DG_LMPD_74_clpn181_p +1000 FU_LMPD_74_clpn181_p < 1000
BU_LMPD_74_clpn181_p: - DG_LMPD_74_clpn181_p +1000 BU_LMPD_74_clpn181_p < 1000
SU_LMPD_74_clpn181_p: FU_LMPD_74_clpn181_p + BU_LMPD_74_clpn181_p < 1
UF_LMPD_74_clpn181_p: F_LMPD_74_clpn181_p -1000 FU_LMPD_74_clpn181_p < 0
UR_LMPD_74_clpn181_p: R_LMPD_74_clpn181_p -1000 BU_LMPD_74_clpn181_p < 0
G_LMPD_75_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +35.5254 LC_nadp_c -35.5254 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_75_clpn181_p + DGo_LMPD_75_clpn181_p = 0
FU_LMPD_75_clpn181_p: DG_LMPD_75_clpn181_p +1000 FU_LMPD_75_clpn181_p < 1000
BU_LMPD_75_clpn181_p: - DG_LMPD_75_clpn181_p +1000 BU_LMPD_75_clpn181_p < 1000
SU_LMPD_75_clpn181_p: FU_LMPD_75_clpn181_p + BU_LMPD_75_clpn181_p < 1
UF_LMPD_75_clpn181_p: F_LMPD_75_clpn181_p -1000 FU_LMPD_75_clpn181_p < 0
UR_LMPD_75_clpn181_p: R_LMPD_75_clpn181_p -1000 BU_LMPD_75_clpn181_p < 0
G_LMPD_76_clpn181_p: -21.3152 LC_accoa_c +21.3152 LC_adp_c -21.3152 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +4.7367 LC_nad_c -4.7367 LC_nadh_c +30.7886 LC_nadp_c -30.7886 LC_nadph_c +1.1842 LC_pg181_c +20.131 LC_pi_c +1.1842 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_76_clpn181_p + DGo_LMPD_76_clpn181_p = 0
FU_LMPD_76_clpn181_p: DG_LMPD_76_clpn181_p +1000 FU_LMPD_76_clpn181_p < 1000
BU_LMPD_76_clpn181_p: - DG_LMPD_76_clpn181_p +1000 BU_LMPD_76_clpn181_p < 1000
SU_LMPD_76_clpn181_p: FU_LMPD_76_clpn181_p + BU_LMPD_76_clpn181_p < 1
UF_LMPD_76_clpn181_p: F_LMPD_76_clpn181_p -1000 FU_LMPD_76_clpn181_p < 0
UR_LMPD_76_clpn181_p: R_LMPD_76_clpn181_p -1000 BU_LMPD_76_clpn181_p < 0
G_LMPD_77_clpn181_p: -21.3152 LC_accoa_c +2.3684 LC_adp_c +2.3684 LC_amp_c -4.7367 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +18.9469 LC_nad_c -18.9469 LC_nadh_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +1.1842 LC_pg181_c +1.1842 LC_pi_c +3.5525 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_77_clpn181_p + DGo_LMPD_77_clpn181_p = 0
FU_LMPD_77_clpn181_p: DG_LMPD_77_clpn181_p +1000 FU_LMPD_77_clpn181_p < 1000
BU_LMPD_77_clpn181_p: - DG_LMPD_77_clpn181_p +1000 BU_LMPD_77_clpn181_p < 1000
SU_LMPD_77_clpn181_p: FU_LMPD_77_clpn181_p + BU_LMPD_77_clpn181_p < 1
UF_LMPD_77_clpn181_p: F_LMPD_77_clpn181_p -1000 FU_LMPD_77_clpn181_p < 0
UR_LMPD_77_clpn181_p: R_LMPD_77_clpn181_p -1000 BU_LMPD_77_clpn181_p < 0
G_LMPD_78_clpn181_p: -21.3152 LC_accoa_c +2.3684 LC_adp_c +2.3684 LC_amp_c -4.7367 LC_atp_c +21.3152 LC_coa_c -2.3684 LC_glyc3p_c +35.5254 LC_nad_c -35.5254 LC_nadh_c +1.1842 LC_pg181_c +1.1842 LC_pi_c +3.5525 LC_ppi_c +0.59209 LC_clpn181_p +0.59209 LC_glyc_p -1.1842 LC_pg181_p - DG_LMPD_78_clpn181_p + DGo_LMPD_78_clpn181_p = 0
FU_LMPD_78_clpn181_p: DG_LMPD_78_clpn181_p +1000 FU_LMPD_78_clpn181_p < 1000
BU_LMPD_78_clpn181_p: - DG_LMPD_78_clpn181_p +1000 BU_LMPD_78_clpn181_p < 1000
SU_LMPD_78_clpn181_p: FU_LMPD_78_clpn181_p + BU_LMPD_78_clpn181_p < 1
UF_LMPD_78_clpn181_p: F_LMPD_78_clpn181_p -1000 FU_LMPD_78_clpn181_p < 0
UR_LMPD_78_clpn181_p: R_LMPD_78_clpn181_p -1000 BU_LMPD_78_clpn181_p < 0
G_LMPD_79_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +7.1051 LC_nad_c -7.1051 LC_nadh_c +32.5649 LC_nadp_c -32.5649 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_79_colipa_e + DGo_LMPD_79_colipa_e = 0
FU_LMPD_79_colipa_e: DG_LMPD_79_colipa_e +1000 FU_LMPD_79_colipa_e < 1000
BU_LMPD_79_colipa_e: - DG_LMPD_79_colipa_e +1000 BU_LMPD_79_colipa_e < 1000
SU_LMPD_79_colipa_e: FU_LMPD_79_colipa_e + BU_LMPD_79_colipa_e < 1
UF_LMPD_79_colipa_e: F_LMPD_79_colipa_e -1000 FU_LMPD_79_colipa_e < 0
UR_LMPD_79_colipa_e: R_LMPD_79_colipa_e -1000 BU_LMPD_79_colipa_e < 0
G_LMPD_80_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +7.6972 LC_nad_c -7.6972 LC_nadh_c +31.9728 LC_nadp_c -31.9728 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_80_colipa_e + DGo_LMPD_80_colipa_e = 0
FU_LMPD_80_colipa_e: DG_LMPD_80_colipa_e +1000 FU_LMPD_80_colipa_e < 1000
BU_LMPD_80_colipa_e: - DG_LMPD_80_colipa_e +1000 BU_LMPD_80_colipa_e < 1000
SU_LMPD_80_colipa_e: FU_LMPD_80_colipa_e + BU_LMPD_80_colipa_e < 1
UF_LMPD_80_colipa_e: F_LMPD_80_colipa_e -1000 FU_LMPD_80_colipa_e < 0
UR_LMPD_80_colipa_e: R_LMPD_80_colipa_e -1000 BU_LMPD_80_colipa_e < 0
G_LMPD_81_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +10.6576 LC_nad_c -10.6576 LC_nadh_c +29.0124 LC_nadp_c -29.0124 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_81_colipa_e + DGo_LMPD_81_colipa_e = 0
FU_LMPD_81_colipa_e: DG_LMPD_81_colipa_e +1000 FU_LMPD_81_colipa_e < 1000
BU_LMPD_81_colipa_e: - DG_LMPD_81_colipa_e +1000 BU_LMPD_81_colipa_e < 1000
SU_LMPD_81_colipa_e: FU_LMPD_81_colipa_e + BU_LMPD_81_colipa_e < 1
UF_LMPD_81_colipa_e: F_LMPD_81_colipa_e -1000 FU_LMPD_81_colipa_e < 0
UR_LMPD_81_colipa_e: R_LMPD_81_colipa_e -1000 BU_LMPD_81_colipa_e < 0
G_LMPD_82_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +11.2497 LC_nad_c -11.2497 LC_nadh_c +28.4203 LC_nadp_c -28.4203 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_82_colipa_e + DGo_LMPD_82_colipa_e = 0
FU_LMPD_82_colipa_e: DG_LMPD_82_colipa_e +1000 FU_LMPD_82_colipa_e < 1000
BU_LMPD_82_colipa_e: - DG_LMPD_82_colipa_e +1000 BU_LMPD_82_colipa_e < 1000
SU_LMPD_82_colipa_e: FU_LMPD_82_colipa_e + BU_LMPD_82_colipa_e < 1
UF_LMPD_82_colipa_e: F_LMPD_82_colipa_e -1000 FU_LMPD_82_colipa_e < 0
UR_LMPD_82_colipa_e: R_LMPD_82_colipa_e -1000 BU_LMPD_82_colipa_e < 0
G_LMPD_83_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +14.2101 LC_nad_c -14.2101 LC_nadh_c +25.4598 LC_nadp_c -25.4598 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_83_colipa_e + DGo_LMPD_83_colipa_e = 0
FU_LMPD_83_colipa_e: DG_LMPD_83_colipa_e +1000 FU_LMPD_83_colipa_e < 1000
BU_LMPD_83_colipa_e: - DG_LMPD_83_colipa_e +1000 BU_LMPD_83_colipa_e < 1000
SU_LMPD_83_colipa_e: FU_LMPD_83_colipa_e + BU_LMPD_83_colipa_e < 1
UF_LMPD_83_colipa_e: F_LMPD_83_colipa_e -1000 FU_LMPD_83_colipa_e < 0
UR_LMPD_83_colipa_e: R_LMPD_83_colipa_e -1000 BU_LMPD_83_colipa_e < 0
G_LMPD_84_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +14.8022 LC_nad_c -14.8022 LC_nadh_c +24.8678 LC_nadp_c -24.8678 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_84_colipa_e + DGo_LMPD_84_colipa_e = 0
FU_LMPD_84_colipa_e: DG_LMPD_84_colipa_e +1000 FU_LMPD_84_colipa_e < 1000
BU_LMPD_84_colipa_e: - DG_LMPD_84_colipa_e +1000 BU_LMPD_84_colipa_e < 1000
SU_LMPD_84_colipa_e: FU_LMPD_84_colipa_e + BU_LMPD_84_colipa_e < 1
UF_LMPD_84_colipa_e: F_LMPD_84_colipa_e -1000 FU_LMPD_84_colipa_e < 0
UR_LMPD_84_colipa_e: R_LMPD_84_colipa_e -1000 BU_LMPD_84_colipa_e < 0
G_LMPD_85_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +17.7627 LC_nad_c -17.7627 LC_nadh_c +21.9073 LC_nadp_c -21.9073 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_85_colipa_e + DGo_LMPD_85_colipa_e = 0
FU_LMPD_85_colipa_e: DG_LMPD_85_colipa_e +1000 FU_LMPD_85_colipa_e < 1000
BU_LMPD_85_colipa_e: - DG_LMPD_85_colipa_e +1000 BU_LMPD_85_colipa_e < 1000
SU_LMPD_85_colipa_e: FU_LMPD_85_colipa_e + BU_LMPD_85_colipa_e < 1
UF_LMPD_85_colipa_e: F_LMPD_85_colipa_e -1000 FU_LMPD_85_colipa_e < 0
UR_LMPD_85_colipa_e: R_LMPD_85_colipa_e -1000 BU_LMPD_85_colipa_e < 0
G_LMPD_86_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +18.3548 LC_nad_c -18.3548 LC_nadh_c +21.3152 LC_nadp_c -21.3152 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_86_colipa_e + DGo_LMPD_86_colipa_e = 0
FU_LMPD_86_colipa_e: DG_LMPD_86_colipa_e +1000 FU_LMPD_86_colipa_e < 1000
BU_LMPD_86_colipa_e: - DG_LMPD_86_colipa_e +1000 BU_LMPD_86_colipa_e < 1000
SU_LMPD_86_colipa_e: FU_LMPD_86_colipa_e + BU_LMPD_86_colipa_e < 1
UF_LMPD_86_colipa_e: F_LMPD_86_colipa_e -1000 FU_LMPD_86_colipa_e < 0
UR_LMPD_86_colipa_e: R_LMPD_86_colipa_e -1000 BU_LMPD_86_colipa_e < 0
G_LMPD_87_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +3.5525 LC_nad_c -3.5525 LC_nadh_c +36.1175 LC_nadp_c -36.1175 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_87_colipa_e + DGo_LMPD_87_colipa_e = 0
FU_LMPD_87_colipa_e: DG_LMPD_87_colipa_e +1000 FU_LMPD_87_colipa_e < 1000
BU_LMPD_87_colipa_e: - DG_LMPD_87_colipa_e +1000 BU_LMPD_87_colipa_e < 1000
SU_LMPD_87_colipa_e: FU_LMPD_87_colipa_e + BU_LMPD_87_colipa_e < 1
UF_LMPD_87_colipa_e: F_LMPD_87_colipa_e -1000 FU_LMPD_87_colipa_e < 0
UR_LMPD_87_colipa_e: R_LMPD_87_colipa_e -1000 BU_LMPD_87_colipa_e < 0
G_LMPD_88_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +39.67 LC_nadp_c -39.67 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_88_colipa_e + DGo_LMPD_88_colipa_e = 0
FU_LMPD_88_colipa_e: DG_LMPD_88_colipa_e +1000 FU_LMPD_88_colipa_e < 1000
BU_LMPD_88_colipa_e: - DG_LMPD_88_colipa_e +1000 BU_LMPD_88_colipa_e < 1000
SU_LMPD_88_colipa_e: FU_LMPD_88_colipa_e + BU_LMPD_88_colipa_e < 1
UF_LMPD_88_colipa_e: F_LMPD_88_colipa_e -1000 FU_LMPD_88_colipa_e < 0
UR_LMPD_88_colipa_e: R_LMPD_88_colipa_e -1000 BU_LMPD_88_colipa_e < 0
G_LMPD_89_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +4.1446 LC_nad_c -4.1446 LC_nadh_c +35.5254 LC_nadp_c -35.5254 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_89_colipa_e + DGo_LMPD_89_colipa_e = 0
FU_LMPD_89_colipa_e: DG_LMPD_89_colipa_e +1000 FU_LMPD_89_colipa_e < 1000
BU_LMPD_89_colipa_e: - DG_LMPD_89_colipa_e +1000 BU_LMPD_89_colipa_e < 1000
SU_LMPD_89_colipa_e: FU_LMPD_89_colipa_e + BU_LMPD_89_colipa_e < 1
UF_LMPD_89_colipa_e: F_LMPD_89_colipa_e -1000 FU_LMPD_89_colipa_e < 0
UR_LMPD_89_colipa_e: R_LMPD_89_colipa_e -1000 BU_LMPD_89_colipa_e < 0
G_LMPD_90_colipa_e: +1.1842 LC_ac_c -25.4598 LC_accoa_c +36.7095 LC_adp_c -36.7095 LC_atp_c +25.4598 LC_coa_c +0.59209 LC_colipa_c -1.1842 LC_f6p_c -2.9604 LC_g1p_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +39.0779 LC_nadp_c -39.0779 LC_nadph_c -1.1842 LC_nh4_c -1.7763 LC_pep_c +27.8282 LC_pi_c +8.2893 LC_ppi_c -1.7763 LC_ru5p-D_c -2.3684 LC_s7p_c - DG_LMPD_90_colipa_e + DGo_LMPD_90_colipa_e = 0
FU_LMPD_90_colipa_e: DG_LMPD_90_colipa_e +1000 FU_LMPD_90_colipa_e < 1000
BU_LMPD_90_colipa_e: - DG_LMPD_90_colipa_e +1000 BU_LMPD_90_colipa_e < 1000
SU_LMPD_90_colipa_e: FU_LMPD_90_colipa_e + BU_LMPD_90_colipa_e < 1
UF_LMPD_90_colipa_e: F_LMPD_90_colipa_e -1000 FU_LMPD_90_colipa_e < 0
UR_LMPD_90_colipa_e: R_LMPD_90_colipa_e -1000 BU_LMPD_90_colipa_e < 0
G_LMPD_91_ctp_c: +2.9604 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -3.5525 LC_atp_c +0.59209 LC_ctp_c +0.59209 LC_dhor-S_c -1.1842 LC_nh4_c -0.59209 LC_orot_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_91_ctp_c + DGo_LMPD_91_ctp_c = 0
FU_LMPD_91_ctp_c: DG_LMPD_91_ctp_c +1000 FU_LMPD_91_ctp_c < 1000
BU_LMPD_91_ctp_c: - DG_LMPD_91_ctp_c +1000 BU_LMPD_91_ctp_c < 1000
SU_LMPD_91_ctp_c: FU_LMPD_91_ctp_c + BU_LMPD_91_ctp_c < 1
UF_LMPD_91_ctp_c: F_LMPD_91_ctp_c -1000 FU_LMPD_91_ctp_c < 0
UR_LMPD_91_ctp_c: R_LMPD_91_ctp_c -1000 BU_LMPD_91_ctp_c < 0
G_LMPD_92_cys-L_c: -0.59209 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +0.59209 LC_akg_c +1.1842 LC_amp_c -1.1842 LC_atp_c +0.59209 LC_coa_c +0.59209 LC_cys-L_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +1.7763 LC_nadp_c -1.7763 LC_nadph_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_92_cys-L_c + DGo_LMPD_92_cys-L_c = 0
FU_LMPD_92_cys-L_c: DG_LMPD_92_cys-L_c +1000 FU_LMPD_92_cys-L_c < 1000
BU_LMPD_92_cys-L_c: - DG_LMPD_92_cys-L_c +1000 BU_LMPD_92_cys-L_c < 1000
SU_LMPD_92_cys-L_c: FU_LMPD_92_cys-L_c + BU_LMPD_92_cys-L_c < 1
UF_LMPD_92_cys-L_c: F_LMPD_92_cys-L_c -1000 FU_LMPD_92_cys-L_c < 0
UR_LMPD_92_cys-L_c: R_LMPD_92_cys-L_c -1000 BU_LMPD_92_cys-L_c < 0
G_LMPD_93_cys-L_c: -0.59209 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +0.59209 LC_amp_c -1.7763 LC_atp_c +0.59209 LC_coa_c +0.59209 LC_cys-L_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +1.7763 LC_nadp_c -1.7763 LC_nadph_c +1.7763 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_93_cys-L_c + DGo_LMPD_93_cys-L_c = 0
FU_LMPD_93_cys-L_c: DG_LMPD_93_cys-L_c +1000 FU_LMPD_93_cys-L_c < 1000
BU_LMPD_93_cys-L_c: - DG_LMPD_93_cys-L_c +1000 BU_LMPD_93_cys-L_c < 1000
SU_LMPD_93_cys-L_c: FU_LMPD_93_cys-L_c + BU_LMPD_93_cys-L_c < 1
UF_LMPD_93_cys-L_c: F_LMPD_93_cys-L_c -1000 FU_LMPD_93_cys-L_c < 0
UR_LMPD_93_cys-L_c: R_LMPD_93_cys-L_c -1000 BU_LMPD_93_cys-L_c < 0
SU_LMPD_94_datp_c: FU_LMPD_94_datp_c + BU_LMPD_94_datp_c < 1
UF_LMPD_94_datp_c: F_LMPD_94_datp_c -1000 FU_LMPD_94_datp_c < 0
UR_LMPD_94_datp_c: R_LMPD_94_datp_c -1000 BU_LMPD_94_datp_c < 0
SU_LMPD_95_datp_c: FU_LMPD_95_datp_c + BU_LMPD_95_datp_c < 1
UF_LMPD_95_datp_c: F_LMPD_95_datp_c -1000 FU_LMPD_95_datp_c < 0
UR_LMPD_95_datp_c: R_LMPD_95_datp_c -1000 BU_LMPD_95_datp_c < 0
SU_LMPD_96_datp_c: FU_LMPD_96_datp_c + BU_LMPD_96_datp_c < 1
UF_LMPD_96_datp_c: F_LMPD_96_datp_c -1000 FU_LMPD_96_datp_c < 0
UR_LMPD_96_datp_c: R_LMPD_96_datp_c -1000 BU_LMPD_96_datp_c < 0
SU_LMPD_97_datp_c: FU_LMPD_97_datp_c + BU_LMPD_97_datp_c < 1
UF_LMPD_97_datp_c: F_LMPD_97_datp_c -1000 FU_LMPD_97_datp_c < 0
UR_LMPD_97_datp_c: R_LMPD_97_datp_c -1000 BU_LMPD_97_datp_c < 0
SU_LMPD_98_dctp_c: FU_LMPD_98_dctp_c + BU_LMPD_98_dctp_c < 1
UF_LMPD_98_dctp_c: F_LMPD_98_dctp_c -1000 FU_LMPD_98_dctp_c < 0
UR_LMPD_98_dctp_c: R_LMPD_98_dctp_c -1000 BU_LMPD_98_dctp_c < 0
SU_LMPD_99_dgtp_c: FU_LMPD_99_dgtp_c + BU_LMPD_99_dgtp_c < 1
UF_LMPD_99_dgtp_c: F_LMPD_99_dgtp_c -1000 FU_LMPD_99_dgtp_c < 0
UR_LMPD_99_dgtp_c: R_LMPD_99_dgtp_c -1000 BU_LMPD_99_dgtp_c < 0
SU_LMPD_100_dgtp_c: FU_LMPD_100_dgtp_c + BU_LMPD_100_dgtp_c < 1
UF_LMPD_100_dgtp_c: F_LMPD_100_dgtp_c -1000 FU_LMPD_100_dgtp_c < 0
UR_LMPD_100_dgtp_c: R_LMPD_100_dgtp_c -1000 BU_LMPD_100_dgtp_c < 0
G_LMPD_101_dgtp_c: -1.1842 LC_3pg_c +5.9209 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_dgtp_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_101_dgtp_c + DGo_LMPD_101_dgtp_c = 0
FU_LMPD_101_dgtp_c: DG_LMPD_101_dgtp_c +1000 FU_LMPD_101_dgtp_c < 1000
BU_LMPD_101_dgtp_c: - DG_LMPD_101_dgtp_c +1000 BU_LMPD_101_dgtp_c < 1000
SU_LMPD_101_dgtp_c: FU_LMPD_101_dgtp_c + BU_LMPD_101_dgtp_c < 1
UF_LMPD_101_dgtp_c: F_LMPD_101_dgtp_c -1000 FU_LMPD_101_dgtp_c < 0
UR_LMPD_101_dgtp_c: R_LMPD_101_dgtp_c -1000 BU_LMPD_101_dgtp_c < 0
SU_LMPD_102_dgtp_c: FU_LMPD_102_dgtp_c + BU_LMPD_102_dgtp_c < 1
UF_LMPD_102_dgtp_c: F_LMPD_102_dgtp_c -1000 FU_LMPD_102_dgtp_c < 0
UR_LMPD_102_dgtp_c: R_LMPD_102_dgtp_c -1000 BU_LMPD_102_dgtp_c < 0
SU_LMPD_103_dgtp_c: FU_LMPD_103_dgtp_c + BU_LMPD_103_dgtp_c < 1
UF_LMPD_103_dgtp_c: F_LMPD_103_dgtp_c -1000 FU_LMPD_103_dgtp_c < 0
UR_LMPD_103_dgtp_c: R_LMPD_103_dgtp_c -1000 BU_LMPD_103_dgtp_c < 0
G_LMPD_104_dgtp_c: -0.59209 LC_3pg_c +6.513 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -7.6972 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_dgtp_c -0.59209 LC_for_c +0.59209 LC_fum_c -0.59209 LC_glu-L_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_104_dgtp_c + DGo_LMPD_104_dgtp_c = 0
FU_LMPD_104_dgtp_c: DG_LMPD_104_dgtp_c +1000 FU_LMPD_104_dgtp_c < 1000
BU_LMPD_104_dgtp_c: - DG_LMPD_104_dgtp_c +1000 BU_LMPD_104_dgtp_c < 1000
SU_LMPD_104_dgtp_c: FU_LMPD_104_dgtp_c + BU_LMPD_104_dgtp_c < 1
UF_LMPD_104_dgtp_c: F_LMPD_104_dgtp_c -1000 FU_LMPD_104_dgtp_c < 0
UR_LMPD_104_dgtp_c: R_LMPD_104_dgtp_c -1000 BU_LMPD_104_dgtp_c < 0
SU_LMPD_105_dttp_c: FU_LMPD_105_dttp_c + BU_LMPD_105_dttp_c < 1
UF_LMPD_105_dttp_c: F_LMPD_105_dttp_c -1000 FU_LMPD_105_dttp_c < 0
UR_LMPD_105_dttp_c: R_LMPD_105_dttp_c -1000 BU_LMPD_105_dttp_c < 0
G_LMPD_106_dttp_c: +3.5525 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -4.1446 LC_atp_c +0.59209 LC_dhor-S_c +0.59209 LC_dttp_c -0.59209 LC_for_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c -0.59209 LC_nh4_c -0.59209 LC_orot_c +1.1842 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_106_dttp_c + DGo_LMPD_106_dttp_c = 0
FU_LMPD_106_dttp_c: DG_LMPD_106_dttp_c +1000 FU_LMPD_106_dttp_c < 1000
BU_LMPD_106_dttp_c: - DG_LMPD_106_dttp_c +1000 BU_LMPD_106_dttp_c < 1000
SU_LMPD_106_dttp_c: FU_LMPD_106_dttp_c + BU_LMPD_106_dttp_c < 1
UF_LMPD_106_dttp_c: F_LMPD_106_dttp_c -1000 FU_LMPD_106_dttp_c < 0
UR_LMPD_106_dttp_c: R_LMPD_106_dttp_c -1000 BU_LMPD_106_dttp_c < 0
SU_LMPD_107_enter_c: FU_LMPD_107_enter_c + BU_LMPD_107_enter_c < 1
UF_LMPD_107_enter_c: F_LMPD_107_enter_c -1000 FU_LMPD_107_enter_c < 0
UR_LMPD_107_enter_c: R_LMPD_107_enter_c -1000 BU_LMPD_107_enter_c < 0
G_LMPD_108_fad_c: -1.1842 LC_3pg_c +10.0655 LC_adp_c +1.1842 LC_akg_c +3.5525 LC_amp_c -1.7763 LC_asp-L_c -13.6181 LC_atp_c -1.1842 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_fad_c -0.59209 LC_for_c +1.7763 LC_fum_c -1.1842 LC_glu-L_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -2.3684 LC_nh4_c +13.6181 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_108_fad_c + DGo_LMPD_108_fad_c = 0
FU_LMPD_108_fad_c: DG_LMPD_108_fad_c +1000 FU_LMPD_108_fad_c < 1000
BU_LMPD_108_fad_c: - DG_LMPD_108_fad_c +1000 BU_LMPD_108_fad_c < 1000
SU_LMPD_108_fad_c: FU_LMPD_108_fad_c + BU_LMPD_108_fad_c < 1
UF_LMPD_108_fad_c: F_LMPD_108_fad_c -1000 FU_LMPD_108_fad_c < 0
UR_LMPD_108_fad_c: R_LMPD_108_fad_c -1000 BU_LMPD_108_fad_c < 0
G_LMPD_109_fad_c: -1.1842 LC_3pg_c +12.4339 LC_adp_c +1.1842 LC_akg_c +2.3684 LC_amp_c -1.7763 LC_asp-L_c -14.8022 LC_atp_c -1.1842 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_fad_c -0.59209 LC_for_c +1.7763 LC_fum_c -1.1842 LC_glu-L_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -2.3684 LC_nh4_c +13.6181 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_109_fad_c + DGo_LMPD_109_fad_c = 0
FU_LMPD_109_fad_c: DG_LMPD_109_fad_c +1000 FU_LMPD_109_fad_c < 1000
BU_LMPD_109_fad_c: - DG_LMPD_109_fad_c +1000 BU_LMPD_109_fad_c < 1000
SU_LMPD_109_fad_c: FU_LMPD_109_fad_c + BU_LMPD_109_fad_c < 1
UF_LMPD_109_fad_c: F_LMPD_109_fad_c -1000 FU_LMPD_109_fad_c < 0
UR_LMPD_109_fad_c: R_LMPD_109_fad_c -1000 BU_LMPD_109_fad_c < 0
G_LMPD_110_fad_c: -2.3684 LC_3pg_c +8.8813 LC_adp_c +2.3684 LC_akg_c +3.5525 LC_amp_c -1.7763 LC_asp-L_c -12.4339 LC_atp_c -1.1842 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_fad_c +0.59209 LC_for_c +1.7763 LC_fum_c -2.3684 LC_glu-L_c +1.1842 LC_gly_c -2.9604 LC_nad_c +2.9604 LC_nadh_c -1.7763 LC_nadp_c +1.7763 LC_nadph_c -2.3684 LC_nh4_c +13.6181 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_110_fad_c + DGo_LMPD_110_fad_c = 0
FU_LMPD_110_fad_c: DG_LMPD_110_fad_c +1000 FU_LMPD_110_fad_c < 1000
BU_LMPD_110_fad_c: - DG_LMPD_110_fad_c +1000 BU_LMPD_110_fad_c < 1000
SU_LMPD_110_fad_c: FU_LMPD_110_fad_c + BU_LMPD_110_fad_c < 1
UF_LMPD_110_fad_c: F_LMPD_110_fad_c -1000 FU_LMPD_110_fad_c < 0
UR_LMPD_110_fad_c: R_LMPD_110_fad_c -1000 BU_LMPD_110_fad_c < 0
G_LMPD_111_fad_c: -2.3684 LC_3pg_c +11.2497 LC_adp_c +2.3684 LC_akg_c +2.3684 LC_amp_c -1.7763 LC_asp-L_c -13.6181 LC_atp_c -1.1842 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_fad_c +0.59209 LC_for_c +1.7763 LC_fum_c -2.3684 LC_glu-L_c +1.1842 LC_gly_c -2.9604 LC_nad_c +2.9604 LC_nadh_c -1.7763 LC_nadp_c +1.7763 LC_nadph_c -2.3684 LC_nh4_c +13.6181 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_111_fad_c + DGo_LMPD_111_fad_c = 0
FU_LMPD_111_fad_c: DG_LMPD_111_fad_c +1000 FU_LMPD_111_fad_c < 1000
BU_LMPD_111_fad_c: - DG_LMPD_111_fad_c +1000 BU_LMPD_111_fad_c < 1000
SU_LMPD_111_fad_c: FU_LMPD_111_fad_c + BU_LMPD_111_fad_c < 1
UF_LMPD_111_fad_c: F_LMPD_111_fad_c -1000 FU_LMPD_111_fad_c < 0
UR_LMPD_111_fad_c: R_LMPD_111_fad_c -1000 BU_LMPD_111_fad_c < 0
G_LMPD_112_GAM_c: -1.1842 LC_3pg_c +3.5525 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -1.1842 LC_asp-L_c -4.7367 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.1842 LC_nh4_c +5.9209 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_112_GAM_c + DGo_LMPD_112_GAM_c = 0
FU_LMPD_112_GAM_c: DG_LMPD_112_GAM_c +1000 FU_LMPD_112_GAM_c < 1000
BU_LMPD_112_GAM_c: - DG_LMPD_112_GAM_c +1000 BU_LMPD_112_GAM_c < 1000
SU_LMPD_112_GAM_c: FU_LMPD_112_GAM_c + BU_LMPD_112_GAM_c < 1
UF_LMPD_112_GAM_c: F_LMPD_112_GAM_c -1000 FU_LMPD_112_GAM_c < 0
UR_LMPD_112_GAM_c: R_LMPD_112_GAM_c -1000 BU_LMPD_112_GAM_c < 0
G_LMPD_113_GAM_c: -1.1842 LC_3pg_c +4.7367 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -1.1842 LC_asp-L_c -5.3288 LC_atp_c -0.59209 LC_co2_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.1842 LC_nh4_c +5.9209 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_113_GAM_c + DGo_LMPD_113_GAM_c = 0
FU_LMPD_113_GAM_c: DG_LMPD_113_GAM_c +1000 FU_LMPD_113_GAM_c < 1000
BU_LMPD_113_GAM_c: - DG_LMPD_113_GAM_c +1000 BU_LMPD_113_GAM_c < 1000
SU_LMPD_113_GAM_c: FU_LMPD_113_GAM_c + BU_LMPD_113_GAM_c < 1
UF_LMPD_113_GAM_c: F_LMPD_113_GAM_c -1000 FU_LMPD_113_GAM_c < 0
UR_LMPD_113_GAM_c: R_LMPD_113_GAM_c -1000 BU_LMPD_113_GAM_c < 0
G_LMPD_114_GAM_c: -0.59209 LC_3pg_c +5.3288 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -1.1842 LC_asp-L_c -5.9209 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_for_c +1.1842 LC_fum_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.1842 LC_nh4_c +5.9209 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_114_GAM_c + DGo_LMPD_114_GAM_c = 0
FU_LMPD_114_GAM_c: DG_LMPD_114_GAM_c +1000 FU_LMPD_114_GAM_c < 1000
BU_LMPD_114_GAM_c: - DG_LMPD_114_GAM_c +1000 BU_LMPD_114_GAM_c < 1000
SU_LMPD_114_GAM_c: FU_LMPD_114_GAM_c + BU_LMPD_114_GAM_c < 1
UF_LMPD_114_GAM_c: F_LMPD_114_GAM_c -1000 FU_LMPD_114_GAM_c < 0
UR_LMPD_114_GAM_c: R_LMPD_114_GAM_c -1000 BU_LMPD_114_GAM_c < 0
G_LMPD_115_GAM_c: -0.59209 LC_3pg_c +4.1446 LC_adp_c +0.59209 LC_akg_c +1.7763 LC_amp_c -1.1842 LC_asp-L_c -5.3288 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_for_c +1.1842 LC_fum_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.1842 LC_nh4_c +5.9209 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_115_GAM_c + DGo_LMPD_115_GAM_c = 0
FU_LMPD_115_GAM_c: DG_LMPD_115_GAM_c +1000 FU_LMPD_115_GAM_c < 1000
BU_LMPD_115_GAM_c: - DG_LMPD_115_GAM_c +1000 BU_LMPD_115_GAM_c < 1000
SU_LMPD_115_GAM_c: FU_LMPD_115_GAM_c + BU_LMPD_115_GAM_c < 1
UF_LMPD_115_GAM_c: F_LMPD_115_GAM_c -1000 FU_LMPD_115_GAM_c < 0
UR_LMPD_115_GAM_c: R_LMPD_115_GAM_c -1000 BU_LMPD_115_GAM_c < 0
G_LMPD_116_gln-L_c: +0.59209 LC_adp_c -0.59209 LC_atp_c +0.59209 LC_gln-L_c -0.59209 LC_glu-L_c -0.59209 LC_nh4_c +0.59209 LC_pi_c - DG_LMPD_116_gln-L_c + DGo_LMPD_116_gln-L_c = 0
FU_LMPD_116_gln-L_c: DG_LMPD_116_gln-L_c +1000 FU_LMPD_116_gln-L_c < 1000
BU_LMPD_116_gln-L_c: - DG_LMPD_116_gln-L_c +1000 BU_LMPD_116_gln-L_c < 1000
SU_LMPD_116_gln-L_c: FU_LMPD_116_gln-L_c + BU_LMPD_116_gln-L_c < 1
UF_LMPD_116_gln-L_c: F_LMPD_116_gln-L_c -1000 FU_LMPD_116_gln-L_c < 0
UR_LMPD_116_gln-L_c: R_LMPD_116_gln-L_c -1000 BU_LMPD_116_gln-L_c < 0
G_LMPD_117_gly_c: +0.59209 LC_acald_c +1.1842 LC_adp_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c +0.59209 LC_gly_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +1.1842 LC_pi_c - DG_LMPD_117_gly_c + DGo_LMPD_117_gly_c = 0
FU_LMPD_117_gly_c: DG_LMPD_117_gly_c +1000 FU_LMPD_117_gly_c < 1000
BU_LMPD_117_gly_c: - DG_LMPD_117_gly_c +1000 BU_LMPD_117_gly_c < 1000
SU_LMPD_117_gly_c: FU_LMPD_117_gly_c + BU_LMPD_117_gly_c < 1
UF_LMPD_117_gly_c: F_LMPD_117_gly_c -1000 FU_LMPD_117_gly_c < 0
UR_LMPD_117_gly_c: R_LMPD_117_gly_c -1000 BU_LMPD_117_gly_c < 0
G_LMPD_118_gthrd_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +1.1842 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -2.3684 LC_atp_c +0.59209 LC_coa_c +0.59209 LC_for_c -1.7763 LC_glu-L_c +0.59209 LC_gthrd_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +3.5525 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_118_gthrd_c + DGo_LMPD_118_gthrd_c = 0
FU_LMPD_118_gthrd_c: DG_LMPD_118_gthrd_c +1000 FU_LMPD_118_gthrd_c < 1000
BU_LMPD_118_gthrd_c: - DG_LMPD_118_gthrd_c +1000 BU_LMPD_118_gthrd_c < 1000
SU_LMPD_118_gthrd_c: FU_LMPD_118_gthrd_c + BU_LMPD_118_gthrd_c < 1
UF_LMPD_118_gthrd_c: F_LMPD_118_gthrd_c -1000 FU_LMPD_118_gthrd_c < 0
UR_LMPD_118_gthrd_c: R_LMPD_118_gthrd_c -1000 BU_LMPD_118_gthrd_c < 0
G_LMPD_119_gthrd_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +2.3684 LC_adp_c +1.1842 LC_akg_c +0.59209 LC_amp_c -2.9604 LC_atp_c +0.59209 LC_coa_c +0.59209 LC_for_c -1.7763 LC_glu-L_c +0.59209 LC_gthrd_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +3.5525 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_119_gthrd_c + DGo_LMPD_119_gthrd_c = 0
FU_LMPD_119_gthrd_c: DG_LMPD_119_gthrd_c +1000 FU_LMPD_119_gthrd_c < 1000
BU_LMPD_119_gthrd_c: - DG_LMPD_119_gthrd_c +1000 BU_LMPD_119_gthrd_c < 1000
SU_LMPD_119_gthrd_c: FU_LMPD_119_gthrd_c + BU_LMPD_119_gthrd_c < 1
UF_LMPD_119_gthrd_c: F_LMPD_119_gthrd_c -1000 FU_LMPD_119_gthrd_c < 0
UR_LMPD_119_gthrd_c: R_LMPD_119_gthrd_c -1000 BU_LMPD_119_gthrd_c < 0
G_LMPD_120_gtp_c: -1.1842 LC_3pg_c +4.7367 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -6.513 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_gtp_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_120_gtp_c + DGo_LMPD_120_gtp_c = 0
FU_LMPD_120_gtp_c: DG_LMPD_120_gtp_c +1000 FU_LMPD_120_gtp_c < 1000
BU_LMPD_120_gtp_c: - DG_LMPD_120_gtp_c +1000 BU_LMPD_120_gtp_c < 1000
SU_LMPD_120_gtp_c: FU_LMPD_120_gtp_c + BU_LMPD_120_gtp_c < 1
UF_LMPD_120_gtp_c: F_LMPD_120_gtp_c -1000 FU_LMPD_120_gtp_c < 0
UR_LMPD_120_gtp_c: R_LMPD_120_gtp_c -1000 BU_LMPD_120_gtp_c < 0
G_LMPD_121_gtp_c: -1.1842 LC_3pg_c +5.9209 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_co2_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_gtp_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_121_gtp_c + DGo_LMPD_121_gtp_c = 0
FU_LMPD_121_gtp_c: DG_LMPD_121_gtp_c +1000 FU_LMPD_121_gtp_c < 1000
BU_LMPD_121_gtp_c: - DG_LMPD_121_gtp_c +1000 BU_LMPD_121_gtp_c < 1000
SU_LMPD_121_gtp_c: FU_LMPD_121_gtp_c + BU_LMPD_121_gtp_c < 1
UF_LMPD_121_gtp_c: F_LMPD_121_gtp_c -1000 FU_LMPD_121_gtp_c < 0
UR_LMPD_121_gtp_c: R_LMPD_121_gtp_c -1000 BU_LMPD_121_gtp_c < 0
G_LMPD_122_gtp_c: -0.59209 LC_3pg_c +5.3288 LC_adp_c +0.59209 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_for_c +0.59209 LC_fum_c -0.59209 LC_glu-L_c +0.59209 LC_gtp_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_122_gtp_c + DGo_LMPD_122_gtp_c = 0
FU_LMPD_122_gtp_c: DG_LMPD_122_gtp_c +1000 FU_LMPD_122_gtp_c < 1000
BU_LMPD_122_gtp_c: - DG_LMPD_122_gtp_c +1000 BU_LMPD_122_gtp_c < 1000
SU_LMPD_122_gtp_c: FU_LMPD_122_gtp_c + BU_LMPD_122_gtp_c < 1
UF_LMPD_122_gtp_c: F_LMPD_122_gtp_c -1000 FU_LMPD_122_gtp_c < 0
UR_LMPD_122_gtp_c: R_LMPD_122_gtp_c -1000 BU_LMPD_122_gtp_c < 0
G_LMPD_123_gtp_c: -0.59209 LC_3pg_c +6.513 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -7.6972 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_for_c +0.59209 LC_fum_c -0.59209 LC_glu-L_c +0.59209 LC_gtp_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.7763 LC_nh4_c +5.9209 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_123_gtp_c + DGo_LMPD_123_gtp_c = 0
FU_LMPD_123_gtp_c: DG_LMPD_123_gtp_c +1000 FU_LMPD_123_gtp_c < 1000
BU_LMPD_123_gtp_c: - DG_LMPD_123_gtp_c +1000 BU_LMPD_123_gtp_c < 1000
SU_LMPD_123_gtp_c: FU_LMPD_123_gtp_c + BU_LMPD_123_gtp_c < 1
UF_LMPD_123_gtp_c: F_LMPD_123_gtp_c -1000 FU_LMPD_123_gtp_c < 0
UR_LMPD_123_gtp_c: R_LMPD_123_gtp_c -1000 BU_LMPD_123_gtp_c < 0
SU_LMPD_124_hemeO_c: FU_LMPD_124_hemeO_c + BU_LMPD_124_hemeO_c < 1
UF_LMPD_124_hemeO_c: F_LMPD_124_hemeO_c -1000 FU_LMPD_124_hemeO_c < 0
UR_LMPD_124_hemeO_c: R_LMPD_124_hemeO_c -1000 BU_LMPD_124_hemeO_c < 0
SU_LMPD_125_hemeO_c: FU_LMPD_125_hemeO_c + BU_LMPD_125_hemeO_c < 1
UF_LMPD_125_hemeO_c: F_LMPD_125_hemeO_c -1000 FU_LMPD_125_hemeO_c < 0
UR_LMPD_125_hemeO_c: R_LMPD_125_hemeO_c -1000 BU_LMPD_125_hemeO_c < 0
G_LMPD_126_his-L_c: +0.59209 LC_adp_c +0.59209 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c -0.59209 LC_for_c +0.59209 LC_fum_c -0.59209 LC_glu-L_c +0.59209 LC_his-L_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_126_his-L_c + DGo_LMPD_126_his-L_c = 0
FU_LMPD_126_his-L_c: DG_LMPD_126_his-L_c +1000 FU_LMPD_126_his-L_c < 1000
BU_LMPD_126_his-L_c: - DG_LMPD_126_his-L_c +1000 BU_LMPD_126_his-L_c < 1000
SU_LMPD_126_his-L_c: FU_LMPD_126_his-L_c + BU_LMPD_126_his-L_c < 1
UF_LMPD_126_his-L_c: F_LMPD_126_his-L_c -1000 FU_LMPD_126_his-L_c < 0
UR_LMPD_126_his-L_c: R_LMPD_126_his-L_c -1000 BU_LMPD_126_his-L_c < 0
G_LMPD_127_his-L_c: +1.7763 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -2.9604 LC_atp_c -0.59209 LC_for_c +0.59209 LC_fum_c -0.59209 LC_glu-L_c +0.59209 LC_his-L_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -0.59209 LC_nh4_c +2.3684 LC_pi_c +1.1842 LC_ppi_c -0.59209 LC_r5p_c - DG_LMPD_127_his-L_c + DGo_LMPD_127_his-L_c = 0
FU_LMPD_127_his-L_c: DG_LMPD_127_his-L_c +1000 FU_LMPD_127_his-L_c < 1000
BU_LMPD_127_his-L_c: - DG_LMPD_127_his-L_c +1000 BU_LMPD_127_his-L_c < 1000
SU_LMPD_127_his-L_c: FU_LMPD_127_his-L_c + BU_LMPD_127_his-L_c < 1
UF_LMPD_127_his-L_c: F_LMPD_127_his-L_c -1000 FU_LMPD_127_his-L_c < 0
UR_LMPD_127_his-L_c: R_LMPD_127_his-L_c -1000 BU_LMPD_127_his-L_c < 0
G_LMPD_128_ile-L_c: +1.1842 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_glu-L_c +0.59209 LC_ile-L_c +1.7763 LC_nadp_c -1.7763 LC_nadph_c +0.59209 LC_nh4_c +1.1842 LC_pi_c -0.59209 LC_pyr_c - DG_LMPD_128_ile-L_c + DGo_LMPD_128_ile-L_c = 0
FU_LMPD_128_ile-L_c: DG_LMPD_128_ile-L_c +1000 FU_LMPD_128_ile-L_c < 1000
BU_LMPD_128_ile-L_c: - DG_LMPD_128_ile-L_c +1000 BU_LMPD_128_ile-L_c < 1000
SU_LMPD_128_ile-L_c: FU_LMPD_128_ile-L_c + BU_LMPD_128_ile-L_c < 1
UF_LMPD_128_ile-L_c: F_LMPD_128_ile-L_c -1000 FU_LMPD_128_ile-L_c < 0
UR_LMPD_128_ile-L_c: R_LMPD_128_ile-L_c -1000 BU_LMPD_128_ile-L_c < 0
G_LMPD_129_leu-L_c: -0.59209 LC_accoa_c +0.59209 LC_akg_c +1.1842 LC_co2_c +0.59209 LC_coa_c -0.59209 LC_glu-L_c +0.59209 LC_leu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_pyr_c - DG_LMPD_129_leu-L_c + DGo_LMPD_129_leu-L_c = 0
FU_LMPD_129_leu-L_c: DG_LMPD_129_leu-L_c +1000 FU_LMPD_129_leu-L_c < 1000
BU_LMPD_129_leu-L_c: - DG_LMPD_129_leu-L_c +1000 BU_LMPD_129_leu-L_c < 1000
SU_LMPD_129_leu-L_c: FU_LMPD_129_leu-L_c + BU_LMPD_129_leu-L_c < 1
UF_LMPD_129_leu-L_c: F_LMPD_129_leu-L_c -1000 FU_LMPD_129_leu-L_c < 0
UR_LMPD_129_leu-L_c: R_LMPD_129_leu-L_c -1000 BU_LMPD_129_leu-L_c < 0
SU_LMPD_130_lipopb_c: FU_LMPD_130_lipopb_c + BU_LMPD_130_lipopb_c < 1
UF_LMPD_130_lipopb_c: F_LMPD_130_lipopb_c -1000 FU_LMPD_130_lipopb_c < 0
UR_LMPD_130_lipopb_c: R_LMPD_130_lipopb_c -1000 BU_LMPD_130_lipopb_c < 0
SU_LMPD_131_lipopb_c: FU_LMPD_131_lipopb_c + BU_LMPD_131_lipopb_c < 1
UF_LMPD_131_lipopb_c: F_LMPD_131_lipopb_c -1000 FU_LMPD_131_lipopb_c < 0
UR_LMPD_131_lipopb_c: R_LMPD_131_lipopb_c -1000 BU_LMPD_131_lipopb_c < 0
SU_LMPD_132_lipopb_c: FU_LMPD_132_lipopb_c + BU_LMPD_132_lipopb_c < 1
UF_LMPD_132_lipopb_c: F_LMPD_132_lipopb_c -1000 FU_LMPD_132_lipopb_c < 0
UR_LMPD_132_lipopb_c: R_LMPD_132_lipopb_c -1000 BU_LMPD_132_lipopb_c < 0
SU_LMPD_133_lipopb_c: FU_LMPD_133_lipopb_c + BU_LMPD_133_lipopb_c < 1
UF_LMPD_133_lipopb_c: F_LMPD_133_lipopb_c -1000 FU_LMPD_133_lipopb_c < 0
UR_LMPD_133_lipopb_c: R_LMPD_133_lipopb_c -1000 BU_LMPD_133_lipopb_c < 0
G_LMPD_134_lys-L_c: +0.59209 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_asp-L_c -0.59209 LC_atp_c +0.59209 LC_co2_c +0.59209 LC_coa_c -0.59209 LC_glu-L_c +0.59209 LC_lys-L_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +0.59209 LC_pi_c -0.59209 LC_pyr_c +0.59209 LC_succ_c -0.59209 LC_succoa_c - DG_LMPD_134_lys-L_c + DGo_LMPD_134_lys-L_c = 0
FU_LMPD_134_lys-L_c: DG_LMPD_134_lys-L_c +1000 FU_LMPD_134_lys-L_c < 1000
BU_LMPD_134_lys-L_c: - DG_LMPD_134_lys-L_c +1000 BU_LMPD_134_lys-L_c < 1000
SU_LMPD_134_lys-L_c: FU_LMPD_134_lys-L_c + BU_LMPD_134_lys-L_c < 1
UF_LMPD_134_lys-L_c: F_LMPD_134_lys-L_c -1000 FU_LMPD_134_lys-L_c < 0
UR_LMPD_134_lys-L_c: R_LMPD_134_lys-L_c -1000 BU_LMPD_134_lys-L_c < 0
G_LMPD_135_malcoa_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -1.1842 LC_accoa_c +7.6972 LC_adp_c +1.1842 LC_akg_c +3.5525 LC_amp_c -1.7763 LC_asp-L_c -11.2497 LC_atp_c +0.59209 LC_co2_c +1.1842 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_malcoa_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -1.1842 LC_nh4_c +8.8813 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_135_malcoa_c + DGo_LMPD_135_malcoa_c = 0
FU_LMPD_135_malcoa_c: DG_LMPD_135_malcoa_c +1000 FU_LMPD_135_malcoa_c < 1000
BU_LMPD_135_malcoa_c: - DG_LMPD_135_malcoa_c +1000 BU_LMPD_135_malcoa_c < 1000
SU_LMPD_135_malcoa_c: FU_LMPD_135_malcoa_c + BU_LMPD_135_malcoa_c < 1
UF_LMPD_135_malcoa_c: F_LMPD_135_malcoa_c -1000 FU_LMPD_135_malcoa_c < 0
UR_LMPD_135_malcoa_c: R_LMPD_135_malcoa_c -1000 BU_LMPD_135_malcoa_c < 0
G_LMPD_136_malcoa_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -1.1842 LC_accoa_c +10.0655 LC_adp_c +1.1842 LC_akg_c +2.3684 LC_amp_c -1.7763 LC_asp-L_c -12.4339 LC_atp_c +0.59209 LC_co2_c +1.1842 LC_coa_c -1.1842 LC_for_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_malcoa_c -1.1842 LC_nad_c +1.1842 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -1.1842 LC_nh4_c +8.8813 LC_pi_c +2.9604 LC_ppi_c -1.1842 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_136_malcoa_c + DGo_LMPD_136_malcoa_c = 0
FU_LMPD_136_malcoa_c: DG_LMPD_136_malcoa_c +1000 FU_LMPD_136_malcoa_c < 1000
BU_LMPD_136_malcoa_c: - DG_LMPD_136_malcoa_c +1000 BU_LMPD_136_malcoa_c < 1000
SU_LMPD_136_malcoa_c: FU_LMPD_136_malcoa_c + BU_LMPD_136_malcoa_c < 1
UF_LMPD_136_malcoa_c: F_LMPD_136_malcoa_c -1000 FU_LMPD_136_malcoa_c < 0
UR_LMPD_136_malcoa_c: R_LMPD_136_malcoa_c -1000 BU_LMPD_136_malcoa_c < 0
G_LMPD_137_met-L_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +0.59209 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -1.7763 LC_atp_c +1.1842 LC_coa_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_met-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c +0.59209 LC_nh4_c +2.9604 LC_pi_c +0.59209 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_137_met-L_c + DGo_LMPD_137_met-L_c = 0
FU_LMPD_137_met-L_c: DG_LMPD_137_met-L_c +1000 FU_LMPD_137_met-L_c < 1000
BU_LMPD_137_met-L_c: - DG_LMPD_137_met-L_c +1000 BU_LMPD_137_met-L_c < 1000
SU_LMPD_137_met-L_c: FU_LMPD_137_met-L_c + BU_LMPD_137_met-L_c < 1
UF_LMPD_137_met-L_c: F_LMPD_137_met-L_c -1000 FU_LMPD_137_met-L_c < 0
UR_LMPD_137_met-L_c: R_LMPD_137_met-L_c -1000 BU_LMPD_137_met-L_c < 0
G_LMPD_138_met-L_c: -1.1842 LC_3pg_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +1.7763 LC_adp_c +1.1842 LC_akg_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +1.1842 LC_coa_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_met-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c +0.59209 LC_nh4_c +2.9604 LC_pi_c +0.59209 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_so4_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_138_met-L_c + DGo_LMPD_138_met-L_c = 0
FU_LMPD_138_met-L_c: DG_LMPD_138_met-L_c +1000 FU_LMPD_138_met-L_c < 1000
BU_LMPD_138_met-L_c: - DG_LMPD_138_met-L_c +1000 BU_LMPD_138_met-L_c < 1000
SU_LMPD_138_met-L_c: FU_LMPD_138_met-L_c + BU_LMPD_138_met-L_c < 1
UF_LMPD_138_met-L_c: F_LMPD_138_met-L_c -1000 FU_LMPD_138_met-L_c < 0
UR_LMPD_138_met-L_c: R_LMPD_138_met-L_c -1000 BU_LMPD_138_met-L_c < 0
G_LMPD_139_mlthf_c: -0.59209 LC_3pg_c +7.6972 LC_adp_c +0.59209 LC_akg_c +2.3684 LC_amp_c -0.59209 LC_asp-L_c -10.0655 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_for_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_glyclt_c +0.59209 LC_mlthf_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_139_mlthf_c + DGo_LMPD_139_mlthf_c = 0
FU_LMPD_139_mlthf_c: DG_LMPD_139_mlthf_c +1000 FU_LMPD_139_mlthf_c < 1000
BU_LMPD_139_mlthf_c: - DG_LMPD_139_mlthf_c +1000 BU_LMPD_139_mlthf_c < 1000
SU_LMPD_139_mlthf_c: FU_LMPD_139_mlthf_c + BU_LMPD_139_mlthf_c < 1
UF_LMPD_139_mlthf_c: F_LMPD_139_mlthf_c -1000 FU_LMPD_139_mlthf_c < 0
UR_LMPD_139_mlthf_c: R_LMPD_139_mlthf_c -1000 BU_LMPD_139_mlthf_c < 0
G_LMPD_140_mlthf_c: -0.59209 LC_3pg_c +8.8813 LC_adp_c +0.59209 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -10.6576 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_for_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_glyclt_c +0.59209 LC_mlthf_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.6576 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c - DG_LMPD_140_mlthf_c + DGo_LMPD_140_mlthf_c = 0
FU_LMPD_140_mlthf_c: DG_LMPD_140_mlthf_c +1000 FU_LMPD_140_mlthf_c < 1000
BU_LMPD_140_mlthf_c: - DG_LMPD_140_mlthf_c +1000 BU_LMPD_140_mlthf_c < 1000
SU_LMPD_140_mlthf_c: FU_LMPD_140_mlthf_c + BU_LMPD_140_mlthf_c < 1
UF_LMPD_140_mlthf_c: F_LMPD_140_mlthf_c -1000 FU_LMPD_140_mlthf_c < 0
UR_LMPD_140_mlthf_c: R_LMPD_140_mlthf_c -1000 BU_LMPD_140_mlthf_c < 0
SU_LMPD_141_mococdp_c: FU_LMPD_141_mococdp_c + BU_LMPD_141_mococdp_c < 1
UF_LMPD_141_mococdp_c: F_LMPD_141_mococdp_c -1000 FU_LMPD_141_mococdp_c < 0
UR_LMPD_141_mococdp_c: R_LMPD_141_mococdp_c -1000 BU_LMPD_141_mococdp_c < 0
SU_LMPD_142_mococdp_c: FU_LMPD_142_mococdp_c + BU_LMPD_142_mococdp_c < 1
UF_LMPD_142_mococdp_c: F_LMPD_142_mococdp_c -1000 FU_LMPD_142_mococdp_c < 0
UR_LMPD_142_mococdp_c: R_LMPD_142_mococdp_c -1000 BU_LMPD_142_mococdp_c < 0
SU_LMPD_143_mococdp_c: FU_LMPD_143_mococdp_c + BU_LMPD_143_mococdp_c < 1
UF_LMPD_143_mococdp_c: F_LMPD_143_mococdp_c -1000 FU_LMPD_143_mococdp_c < 0
UR_LMPD_143_mococdp_c: R_LMPD_143_mococdp_c -1000 BU_LMPD_143_mococdp_c < 0
SU_LMPD_144_mococdp_c: FU_LMPD_144_mococdp_c + BU_LMPD_144_mococdp_c < 1
UF_LMPD_144_mococdp_c: F_LMPD_144_mococdp_c -1000 FU_LMPD_144_mococdp_c < 0
UR_LMPD_144_mococdp_c: R_LMPD_144_mococdp_c -1000 BU_LMPD_144_mococdp_c < 0
SU_LMPD_145_mococdp_c: FU_LMPD_145_mococdp_c + BU_LMPD_145_mococdp_c < 1
UF_LMPD_145_mococdp_c: F_LMPD_145_mococdp_c -1000 FU_LMPD_145_mococdp_c < 0
UR_LMPD_145_mococdp_c: R_LMPD_145_mococdp_c -1000 BU_LMPD_145_mococdp_c < 0
SU_LMPD_146_mococdp_c: FU_LMPD_146_mococdp_c + BU_LMPD_146_mococdp_c < 1
UF_LMPD_146_mococdp_c: F_LMPD_146_mococdp_c -1000 FU_LMPD_146_mococdp_c < 0
UR_LMPD_146_mococdp_c: R_LMPD_146_mococdp_c -1000 BU_LMPD_146_mococdp_c < 0
SU_LMPD_147_mococdp_c: FU_LMPD_147_mococdp_c + BU_LMPD_147_mococdp_c < 1
UF_LMPD_147_mococdp_c: F_LMPD_147_mococdp_c -1000 FU_LMPD_147_mococdp_c < 0
UR_LMPD_147_mococdp_c: R_LMPD_147_mococdp_c -1000 BU_LMPD_147_mococdp_c < 0
SU_LMPD_148_mococdp_c: FU_LMPD_148_mococdp_c + BU_LMPD_148_mococdp_c < 1
UF_LMPD_148_mococdp_c: F_LMPD_148_mococdp_c -1000 FU_LMPD_148_mococdp_c < 0
UR_LMPD_148_mococdp_c: R_LMPD_148_mococdp_c -1000 BU_LMPD_148_mococdp_c < 0
SU_LMPD_149_mocogdp_c: FU_LMPD_149_mocogdp_c + BU_LMPD_149_mocogdp_c < 1
UF_LMPD_149_mocogdp_c: F_LMPD_149_mocogdp_c -1000 FU_LMPD_149_mocogdp_c < 0
UR_LMPD_149_mocogdp_c: R_LMPD_149_mocogdp_c -1000 BU_LMPD_149_mocogdp_c < 0
SU_LMPD_150_mocogdp_c: FU_LMPD_150_mocogdp_c + BU_LMPD_150_mocogdp_c < 1
UF_LMPD_150_mocogdp_c: F_LMPD_150_mocogdp_c -1000 FU_LMPD_150_mocogdp_c < 0
UR_LMPD_150_mocogdp_c: R_LMPD_150_mocogdp_c -1000 BU_LMPD_150_mocogdp_c < 0
SU_LMPD_151_mocogdp_c: FU_LMPD_151_mocogdp_c + BU_LMPD_151_mocogdp_c < 1
UF_LMPD_151_mocogdp_c: F_LMPD_151_mocogdp_c -1000 FU_LMPD_151_mocogdp_c < 0
UR_LMPD_151_mocogdp_c: R_LMPD_151_mocogdp_c -1000 BU_LMPD_151_mocogdp_c < 0
SU_LMPD_152_mocogdp_c: FU_LMPD_152_mocogdp_c + BU_LMPD_152_mocogdp_c < 1
UF_LMPD_152_mocogdp_c: F_LMPD_152_mocogdp_c -1000 FU_LMPD_152_mocogdp_c < 0
UR_LMPD_152_mocogdp_c: R_LMPD_152_mocogdp_c -1000 BU_LMPD_152_mocogdp_c < 0
SU_LMPD_153_mql8_c: FU_LMPD_153_mql8_c + BU_LMPD_153_mql8_c < 1
UF_LMPD_153_mql8_c: F_LMPD_153_mql8_c -1000 FU_LMPD_153_mql8_c < 0
UR_LMPD_153_mql8_c: R_LMPD_153_mql8_c -1000 BU_LMPD_153_mql8_c < 0
G_LMPD_154_murein3p3p_p: -2.3684 LC_accoa_c +13.026 LC_adp_c +2.3684 LC_akg_c -1.1842 LC_asp-L_c -13.026 LC_atp_c +3.5525 LC_coa_c -2.3684 LC_f6p_c -3.5525 LC_glu-L_c +3.5525 LC_nadp_c -3.5525 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +11.8418 LC_pi_c +2.3684 LC_ppi_c -2.3684 LC_pyr_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +0.59209 LC_murein3p3p_p - DG_LMPD_154_murein3p3p_p + DGo_LMPD_154_murein3p3p_p = 0
FU_LMPD_154_murein3p3p_p: DG_LMPD_154_murein3p3p_p +1000 FU_LMPD_154_murein3p3p_p < 1000
BU_LMPD_154_murein3p3p_p: - DG_LMPD_154_murein3p3p_p +1000 BU_LMPD_154_murein3p3p_p < 1000
SU_LMPD_154_murein3p3p_p: FU_LMPD_154_murein3p3p_p + BU_LMPD_154_murein3p3p_p < 1
UF_LMPD_154_murein3p3p_p: F_LMPD_154_murein3p3p_p -1000 FU_LMPD_154_murein3p3p_p < 0
UR_LMPD_154_murein3p3p_p: R_LMPD_154_murein3p3p_p -1000 BU_LMPD_154_murein3p3p_p < 0
G_LMPD_155_murein3px4p_p: -2.3684 LC_accoa_c +13.026 LC_adp_c +4.7367 LC_akg_c -1.1842 LC_asp-L_c -13.026 LC_atp_c +3.5525 LC_coa_c -2.3684 LC_f6p_c -5.9209 LC_glu-L_c +3.5525 LC_nadp_c -3.5525 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +11.8418 LC_pi_c +2.3684 LC_ppi_c -4.7367 LC_pyr_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +0.59209 LC_ala-D_p +0.59209 LC_alaala_p +0.59209 LC_murein3px4p_p - DG_LMPD_155_murein3px4p_p + DGo_LMPD_155_murein3px4p_p = 0
FU_LMPD_155_murein3px4p_p: DG_LMPD_155_murein3px4p_p +1000 FU_LMPD_155_murein3px4p_p < 1000
BU_LMPD_155_murein3px4p_p: - DG_LMPD_155_murein3px4p_p +1000 BU_LMPD_155_murein3px4p_p < 1000
SU_LMPD_155_murein3px4p_p: FU_LMPD_155_murein3px4p_p + BU_LMPD_155_murein3px4p_p < 1
UF_LMPD_155_murein3px4p_p: F_LMPD_155_murein3px4p_p -1000 FU_LMPD_155_murein3px4p_p < 0
UR_LMPD_155_murein3px4p_p: R_LMPD_155_murein3px4p_p -1000 BU_LMPD_155_murein3px4p_p < 0
G_LMPD_156_murein4p4p_p: -2.3684 LC_accoa_c +13.026 LC_adp_c +4.7367 LC_akg_c -1.1842 LC_asp-L_c -13.026 LC_atp_c +3.5525 LC_coa_c -2.3684 LC_f6p_c -5.9209 LC_glu-L_c +3.5525 LC_nadp_c -3.5525 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +11.8418 LC_pi_c +2.3684 LC_ppi_c -4.7367 LC_pyr_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +1.1842 LC_ala-D_p +0.59209 LC_murein4p4p_p - DG_LMPD_156_murein4p4p_p + DGo_LMPD_156_murein4p4p_p = 0
FU_LMPD_156_murein4p4p_p: DG_LMPD_156_murein4p4p_p +1000 FU_LMPD_156_murein4p4p_p < 1000
BU_LMPD_156_murein4p4p_p: - DG_LMPD_156_murein4p4p_p +1000 BU_LMPD_156_murein4p4p_p < 1000
SU_LMPD_156_murein4p4p_p: FU_LMPD_156_murein4p4p_p + BU_LMPD_156_murein4p4p_p < 1
UF_LMPD_156_murein4p4p_p: F_LMPD_156_murein4p4p_p -1000 FU_LMPD_156_murein4p4p_p < 0
UR_LMPD_156_murein4p4p_p: R_LMPD_156_murein4p4p_p -1000 BU_LMPD_156_murein4p4p_p < 0
G_LMPD_157_murein4px4p_p: -2.3684 LC_accoa_c +13.026 LC_adp_c +4.7367 LC_akg_c -1.1842 LC_asp-L_c -13.026 LC_atp_c +3.5525 LC_coa_c -2.3684 LC_f6p_c -5.9209 LC_glu-L_c +3.5525 LC_nadp_c -3.5525 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +11.8418 LC_pi_c +2.3684 LC_ppi_c -4.7367 LC_pyr_c +1.1842 LC_succ_c -1.1842 LC_succoa_c +1.1842 LC_ala-D_p +0.59209 LC_murein4px4p_p - DG_LMPD_157_murein4px4p_p + DGo_LMPD_157_murein4px4p_p = 0
FU_LMPD_157_murein4px4p_p: DG_LMPD_157_murein4px4p_p +1000 FU_LMPD_157_murein4px4p_p < 1000
BU_LMPD_157_murein4px4p_p: - DG_LMPD_157_murein4px4p_p +1000 BU_LMPD_157_murein4px4p_p < 1000
SU_LMPD_157_murein4px4p_p: FU_LMPD_157_murein4px4p_p + BU_LMPD_157_murein4px4p_p < 1
UF_LMPD_157_murein4px4p_p: F_LMPD_157_murein4px4p_p -1000 FU_LMPD_157_murein4px4p_p < 0
UR_LMPD_157_murein4px4p_p: R_LMPD_157_murein4px4p_p -1000 BU_LMPD_157_murein4px4p_p < 0
G_LMPD_158_murein4px4px4p_p: -3.5525 LC_accoa_c +19.5389 LC_adp_c +7.1051 LC_akg_c -1.7763 LC_asp-L_c -19.5389 LC_atp_c +5.3288 LC_coa_c -3.5525 LC_f6p_c -8.8813 LC_glu-L_c +5.3288 LC_nadp_c -5.3288 LC_nadph_c -3.5525 LC_nh4_c -1.7763 LC_pep_c +17.7627 LC_pi_c +3.5525 LC_ppi_c -7.1051 LC_pyr_c +1.7763 LC_succ_c -1.7763 LC_succoa_c +1.7763 LC_ala-D_p +0.59209 LC_murein4px4px4p_p - DG_LMPD_158_murein4px4px4p_p + DGo_LMPD_158_murein4px4px4p_p = 0
FU_LMPD_158_murein4px4px4p_p: DG_LMPD_158_murein4px4px4p_p +1000 FU_LMPD_158_murein4px4px4p_p < 1000
BU_LMPD_158_murein4px4px4p_p: - DG_LMPD_158_murein4px4px4p_p +1000 BU_LMPD_158_murein4px4px4p_p < 1000
SU_LMPD_158_murein4px4px4p_p: FU_LMPD_158_murein4px4px4p_p + BU_LMPD_158_murein4px4px4p_p < 1
UF_LMPD_158_murein4px4px4p_p: F_LMPD_158_murein4px4px4p_p -1000 FU_LMPD_158_murein4px4px4p_p < 0
UR_LMPD_158_murein4px4px4p_p: R_LMPD_158_murein4px4px4p_p -1000 BU_LMPD_158_murein4px4px4p_p < 0
G_LMPD_159_nad_c: -1.1842 LC_3pg_c +3.5525 LC_adp_c +1.1842 LC_akg_c +2.9604 LC_amp_c -1.1842 LC_asp-L_c -6.513 LC_atp_c -0.59209 LC_dhap_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -0.59209 LC_iasp_c -0.59209 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.7763 LC_nh4_c +6.513 LC_pi_c +2.3684 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_159_nad_c + DGo_LMPD_159_nad_c = 0
FU_LMPD_159_nad_c: DG_LMPD_159_nad_c +1000 FU_LMPD_159_nad_c < 1000
BU_LMPD_159_nad_c: - DG_LMPD_159_nad_c +1000 BU_LMPD_159_nad_c < 1000
SU_LMPD_159_nad_c: FU_LMPD_159_nad_c + BU_LMPD_159_nad_c < 1
UF_LMPD_159_nad_c: F_LMPD_159_nad_c -1000 FU_LMPD_159_nad_c < 0
UR_LMPD_159_nad_c: R_LMPD_159_nad_c -1000 BU_LMPD_159_nad_c < 0
G_LMPD_160_nad_c: -1.1842 LC_3pg_c +4.7367 LC_adp_c +1.1842 LC_akg_c +2.3684 LC_amp_c -1.1842 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_dhap_c +1.1842 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -0.59209 LC_iasp_c -0.59209 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nadp_c +1.1842 LC_nadph_c -1.7763 LC_nh4_c +6.513 LC_pi_c +2.3684 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_160_nad_c + DGo_LMPD_160_nad_c = 0
FU_LMPD_160_nad_c: DG_LMPD_160_nad_c +1000 FU_LMPD_160_nad_c < 1000
BU_LMPD_160_nad_c: - DG_LMPD_160_nad_c +1000 BU_LMPD_160_nad_c < 1000
SU_LMPD_160_nad_c: FU_LMPD_160_nad_c + BU_LMPD_160_nad_c < 1
UF_LMPD_160_nad_c: F_LMPD_160_nad_c -1000 FU_LMPD_160_nad_c < 0
UR_LMPD_160_nad_c: R_LMPD_160_nad_c -1000 BU_LMPD_160_nad_c < 0
G_LMPD_161_nad_c: -0.59209 LC_3pg_c +4.1446 LC_adp_c +0.59209 LC_akg_c +2.9604 LC_amp_c -1.1842 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_dhap_c -0.59209 LC_for_c +1.1842 LC_fum_c -0.59209 LC_glu-L_c -0.59209 LC_iasp_c +0.59209 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.7763 LC_nh4_c +6.513 LC_pi_c +2.3684 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_161_nad_c + DGo_LMPD_161_nad_c = 0
FU_LMPD_161_nad_c: DG_LMPD_161_nad_c +1000 FU_LMPD_161_nad_c < 1000
BU_LMPD_161_nad_c: - DG_LMPD_161_nad_c +1000 BU_LMPD_161_nad_c < 1000
SU_LMPD_161_nad_c: FU_LMPD_161_nad_c + BU_LMPD_161_nad_c < 1
UF_LMPD_161_nad_c: F_LMPD_161_nad_c -1000 FU_LMPD_161_nad_c < 0
UR_LMPD_161_nad_c: R_LMPD_161_nad_c -1000 BU_LMPD_161_nad_c < 0
G_LMPD_162_nad_c: -0.59209 LC_3pg_c +5.3288 LC_adp_c +0.59209 LC_akg_c +2.3684 LC_amp_c -1.1842 LC_asp-L_c -7.6972 LC_atp_c -0.59209 LC_dhap_c -0.59209 LC_for_c +1.1842 LC_fum_c -0.59209 LC_glu-L_c -0.59209 LC_iasp_c +0.59209 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.7763 LC_nh4_c +6.513 LC_pi_c +2.3684 LC_ppi_c -1.1842 LC_r5p_c - DG_LMPD_162_nad_c + DGo_LMPD_162_nad_c = 0
FU_LMPD_162_nad_c: DG_LMPD_162_nad_c +1000 FU_LMPD_162_nad_c < 1000
BU_LMPD_162_nad_c: - DG_LMPD_162_nad_c +1000 BU_LMPD_162_nad_c < 1000
SU_LMPD_162_nad_c: FU_LMPD_162_nad_c + BU_LMPD_162_nad_c < 1
UF_LMPD_162_nad_c: F_LMPD_162_nad_c -1000 FU_LMPD_162_nad_c < 0
UR_LMPD_162_nad_c: R_LMPD_162_nad_c -1000 BU_LMPD_162_nad_c < 0
G_LMPD_163_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +6.513 LC_nad_c -6.513 LC_nadh_c +9.4734 LC_nadp_c -9.4734 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_163_pe160_c + DGo_LMPD_163_pe160_c = 0
FU_LMPD_163_pe160_c: DG_LMPD_163_pe160_c +1000 FU_LMPD_163_pe160_c < 1000
BU_LMPD_163_pe160_c: - DG_LMPD_163_pe160_c +1000 BU_LMPD_163_pe160_c < 1000
SU_LMPD_163_pe160_c: FU_LMPD_163_pe160_c + BU_LMPD_163_pe160_c < 1
UF_LMPD_163_pe160_c: F_LMPD_163_pe160_c -1000 FU_LMPD_163_pe160_c < 0
UR_LMPD_163_pe160_c: R_LMPD_163_pe160_c -1000 BU_LMPD_163_pe160_c < 0
G_LMPD_164_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +7.6972 LC_nad_c -7.6972 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_164_pe160_c + DGo_LMPD_164_pe160_c = 0
FU_LMPD_164_pe160_c: DG_LMPD_164_pe160_c +1000 FU_LMPD_164_pe160_c < 1000
BU_LMPD_164_pe160_c: - DG_LMPD_164_pe160_c +1000 BU_LMPD_164_pe160_c < 1000
SU_LMPD_164_pe160_c: FU_LMPD_164_pe160_c + BU_LMPD_164_pe160_c < 1
UF_LMPD_164_pe160_c: F_LMPD_164_pe160_c -1000 FU_LMPD_164_pe160_c < 0
UR_LMPD_164_pe160_c: R_LMPD_164_pe160_c -1000 BU_LMPD_164_pe160_c < 0
G_LMPD_165_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +1.7763 LC_nad_c -1.7763 LC_nadh_c +14.2101 LC_nadp_c -14.2101 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_165_pe160_c + DGo_LMPD_165_pe160_c = 0
FU_LMPD_165_pe160_c: DG_LMPD_165_pe160_c +1000 FU_LMPD_165_pe160_c < 1000
BU_LMPD_165_pe160_c: - DG_LMPD_165_pe160_c +1000 BU_LMPD_165_pe160_c < 1000
SU_LMPD_165_pe160_c: FU_LMPD_165_pe160_c + BU_LMPD_165_pe160_c < 1
UF_LMPD_165_pe160_c: F_LMPD_165_pe160_c -1000 FU_LMPD_165_pe160_c < 0
UR_LMPD_165_pe160_c: R_LMPD_165_pe160_c -1000 BU_LMPD_165_pe160_c < 0
G_LMPD_166_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +2.9604 LC_nad_c -2.9604 LC_nadh_c +13.026 LC_nadp_c -13.026 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_166_pe160_c + DGo_LMPD_166_pe160_c = 0
FU_LMPD_166_pe160_c: DG_LMPD_166_pe160_c +1000 FU_LMPD_166_pe160_c < 1000
BU_LMPD_166_pe160_c: - DG_LMPD_166_pe160_c +1000 BU_LMPD_166_pe160_c < 1000
SU_LMPD_166_pe160_c: FU_LMPD_166_pe160_c + BU_LMPD_166_pe160_c < 1
UF_LMPD_166_pe160_c: F_LMPD_166_pe160_c -1000 FU_LMPD_166_pe160_c < 0
UR_LMPD_166_pe160_c: R_LMPD_166_pe160_c -1000 BU_LMPD_166_pe160_c < 0
G_LMPD_167_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +4.1446 LC_nad_c -4.1446 LC_nadh_c +11.8418 LC_nadp_c -11.8418 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_167_pe160_c + DGo_LMPD_167_pe160_c = 0
FU_LMPD_167_pe160_c: DG_LMPD_167_pe160_c +1000 FU_LMPD_167_pe160_c < 1000
BU_LMPD_167_pe160_c: - DG_LMPD_167_pe160_c +1000 BU_LMPD_167_pe160_c < 1000
SU_LMPD_167_pe160_c: FU_LMPD_167_pe160_c + BU_LMPD_167_pe160_c < 1
UF_LMPD_167_pe160_c: F_LMPD_167_pe160_c -1000 FU_LMPD_167_pe160_c < 0
UR_LMPD_167_pe160_c: R_LMPD_167_pe160_c -1000 BU_LMPD_167_pe160_c < 0
G_LMPD_168_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +5.3288 LC_nad_c -5.3288 LC_nadh_c +10.6576 LC_nadp_c -10.6576 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_168_pe160_c + DGo_LMPD_168_pe160_c = 0
FU_LMPD_168_pe160_c: DG_LMPD_168_pe160_c +1000 FU_LMPD_168_pe160_c < 1000
BU_LMPD_168_pe160_c: - DG_LMPD_168_pe160_c +1000 BU_LMPD_168_pe160_c < 1000
SU_LMPD_168_pe160_c: FU_LMPD_168_pe160_c + BU_LMPD_168_pe160_c < 1
UF_LMPD_168_pe160_c: F_LMPD_168_pe160_c -1000 FU_LMPD_168_pe160_c < 0
UR_LMPD_168_pe160_c: R_LMPD_168_pe160_c -1000 BU_LMPD_168_pe160_c < 0
G_LMPD_169_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_169_pe160_c + DGo_LMPD_169_pe160_c = 0
FU_LMPD_169_pe160_c: DG_LMPD_169_pe160_c +1000 FU_LMPD_169_pe160_c < 1000
BU_LMPD_169_pe160_c: - DG_LMPD_169_pe160_c +1000 BU_LMPD_169_pe160_c < 1000
SU_LMPD_169_pe160_c: FU_LMPD_169_pe160_c + BU_LMPD_169_pe160_c < 1
UF_LMPD_169_pe160_c: F_LMPD_169_pe160_c -1000 FU_LMPD_169_pe160_c < 0
UR_LMPD_169_pe160_c: R_LMPD_169_pe160_c -1000 BU_LMPD_169_pe160_c < 0
G_LMPD_170_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +15.3943 LC_nadp_c -15.3943 LC_nadph_c +0.59209 LC_pe160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_170_pe160_c + DGo_LMPD_170_pe160_c = 0
FU_LMPD_170_pe160_c: DG_LMPD_170_pe160_c +1000 FU_LMPD_170_pe160_c < 1000
BU_LMPD_170_pe160_c: - DG_LMPD_170_pe160_c +1000 BU_LMPD_170_pe160_c < 1000
SU_LMPD_170_pe160_c: FU_LMPD_170_pe160_c + BU_LMPD_170_pe160_c < 1
UF_LMPD_170_pe160_c: F_LMPD_170_pe160_c -1000 FU_LMPD_170_pe160_c < 0
UR_LMPD_170_pe160_c: R_LMPD_170_pe160_c -1000 BU_LMPD_170_pe160_c < 0
G_LMPD_171_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -2.3684 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +7.6972 LC_nad_c -7.6972 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pe160_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_171_pe160_c + DGo_LMPD_171_pe160_c = 0
FU_LMPD_171_pe160_c: DG_LMPD_171_pe160_c +1000 FU_LMPD_171_pe160_c < 1000
BU_LMPD_171_pe160_c: - DG_LMPD_171_pe160_c +1000 BU_LMPD_171_pe160_c < 1000
SU_LMPD_171_pe160_c: FU_LMPD_171_pe160_c + BU_LMPD_171_pe160_c < 1
UF_LMPD_171_pe160_c: F_LMPD_171_pe160_c -1000 FU_LMPD_171_pe160_c < 0
UR_LMPD_171_pe160_c: R_LMPD_171_pe160_c -1000 BU_LMPD_171_pe160_c < 0
G_LMPD_172_pe160_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -2.3684 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +15.9864 LC_nad_c -15.9864 LC_nadh_c +0.59209 LC_pe160_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_172_pe160_c + DGo_LMPD_172_pe160_c = 0
FU_LMPD_172_pe160_c: DG_LMPD_172_pe160_c +1000 FU_LMPD_172_pe160_c < 1000
BU_LMPD_172_pe160_c: - DG_LMPD_172_pe160_c +1000 BU_LMPD_172_pe160_c < 1000
SU_LMPD_172_pe160_c: FU_LMPD_172_pe160_c + BU_LMPD_172_pe160_c < 1
UF_LMPD_172_pe160_c: F_LMPD_172_pe160_c -1000 FU_LMPD_172_pe160_c < 0
UR_LMPD_172_pe160_c: R_LMPD_172_pe160_c -1000 BU_LMPD_172_pe160_c < 0
G_LMPD_173_pe161_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +6.513 LC_nad_c -6.513 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pe161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_173_pe161_c + DGo_LMPD_173_pe161_c = 0
FU_LMPD_173_pe161_c: DG_LMPD_173_pe161_c +1000 FU_LMPD_173_pe161_c < 1000
BU_LMPD_173_pe161_c: - DG_LMPD_173_pe161_c +1000 BU_LMPD_173_pe161_c < 1000
SU_LMPD_173_pe161_c: FU_LMPD_173_pe161_c + BU_LMPD_173_pe161_c < 1
UF_LMPD_173_pe161_c: F_LMPD_173_pe161_c -1000 FU_LMPD_173_pe161_c < 0
UR_LMPD_173_pe161_c: R_LMPD_173_pe161_c -1000 BU_LMPD_173_pe161_c < 0
G_LMPD_174_pe161_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +1.7763 LC_nad_c -1.7763 LC_nadh_c +13.026 LC_nadp_c -13.026 LC_nadph_c +0.59209 LC_pe161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_174_pe161_c + DGo_LMPD_174_pe161_c = 0
FU_LMPD_174_pe161_c: DG_LMPD_174_pe161_c +1000 FU_LMPD_174_pe161_c < 1000
BU_LMPD_174_pe161_c: - DG_LMPD_174_pe161_c +1000 BU_LMPD_174_pe161_c < 1000
SU_LMPD_174_pe161_c: FU_LMPD_174_pe161_c + BU_LMPD_174_pe161_c < 1
UF_LMPD_174_pe161_c: F_LMPD_174_pe161_c -1000 FU_LMPD_174_pe161_c < 0
UR_LMPD_174_pe161_c: R_LMPD_174_pe161_c -1000 BU_LMPD_174_pe161_c < 0
G_LMPD_175_pe161_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +2.9604 LC_nad_c -2.9604 LC_nadh_c +11.8418 LC_nadp_c -11.8418 LC_nadph_c +0.59209 LC_pe161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_175_pe161_c + DGo_LMPD_175_pe161_c = 0
FU_LMPD_175_pe161_c: DG_LMPD_175_pe161_c +1000 FU_LMPD_175_pe161_c < 1000
BU_LMPD_175_pe161_c: - DG_LMPD_175_pe161_c +1000 BU_LMPD_175_pe161_c < 1000
SU_LMPD_175_pe161_c: FU_LMPD_175_pe161_c + BU_LMPD_175_pe161_c < 1
UF_LMPD_175_pe161_c: F_LMPD_175_pe161_c -1000 FU_LMPD_175_pe161_c < 0
UR_LMPD_175_pe161_c: R_LMPD_175_pe161_c -1000 BU_LMPD_175_pe161_c < 0
G_LMPD_176_pe161_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +4.1446 LC_nad_c -4.1446 LC_nadh_c +10.6576 LC_nadp_c -10.6576 LC_nadph_c +0.59209 LC_pe161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_176_pe161_c + DGo_LMPD_176_pe161_c = 0
FU_LMPD_176_pe161_c: DG_LMPD_176_pe161_c +1000 FU_LMPD_176_pe161_c < 1000
BU_LMPD_176_pe161_c: - DG_LMPD_176_pe161_c +1000 BU_LMPD_176_pe161_c < 1000
SU_LMPD_176_pe161_c: FU_LMPD_176_pe161_c + BU_LMPD_176_pe161_c < 1
UF_LMPD_176_pe161_c: F_LMPD_176_pe161_c -1000 FU_LMPD_176_pe161_c < 0
UR_LMPD_176_pe161_c: R_LMPD_176_pe161_c -1000 BU_LMPD_176_pe161_c < 0
G_LMPD_177_pe161_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +5.3288 LC_nad_c -5.3288 LC_nadh_c +9.4734 LC_nadp_c -9.4734 LC_nadph_c +0.59209 LC_pe161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_177_pe161_c + DGo_LMPD_177_pe161_c = 0
FU_LMPD_177_pe161_c: DG_LMPD_177_pe161_c +1000 FU_LMPD_177_pe161_c < 1000
BU_LMPD_177_pe161_c: - DG_LMPD_177_pe161_c +1000 BU_LMPD_177_pe161_c < 1000
SU_LMPD_177_pe161_c: FU_LMPD_177_pe161_c + BU_LMPD_177_pe161_c < 1
UF_LMPD_177_pe161_c: F_LMPD_177_pe161_c -1000 FU_LMPD_177_pe161_c < 0
UR_LMPD_177_pe161_c: R_LMPD_177_pe161_c -1000 BU_LMPD_177_pe161_c < 0
G_LMPD_178_pe161_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +15.3943 LC_nadp_c -15.3943 LC_nadph_c +0.59209 LC_pe161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_178_pe161_c + DGo_LMPD_178_pe161_c = 0
FU_LMPD_178_pe161_c: DG_LMPD_178_pe161_c +1000 FU_LMPD_178_pe161_c < 1000
BU_LMPD_178_pe161_c: - DG_LMPD_178_pe161_c +1000 BU_LMPD_178_pe161_c < 1000
SU_LMPD_178_pe161_c: FU_LMPD_178_pe161_c + BU_LMPD_178_pe161_c < 1
UF_LMPD_178_pe161_c: F_LMPD_178_pe161_c -1000 FU_LMPD_178_pe161_c < 0
UR_LMPD_178_pe161_c: R_LMPD_178_pe161_c -1000 BU_LMPD_178_pe161_c < 0
G_LMPD_179_pe161_c: -0.59209 LC_3pg_c -9.4734 LC_accoa_c +9.4734 LC_adp_c +0.59209 LC_akg_c -9.4734 LC_atp_c +0.59209 LC_co2_c +9.4734 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +14.2101 LC_nadp_c -14.2101 LC_nadph_c +0.59209 LC_pe161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_179_pe161_c + DGo_LMPD_179_pe161_c = 0
FU_LMPD_179_pe161_c: DG_LMPD_179_pe161_c +1000 FU_LMPD_179_pe161_c < 1000
BU_LMPD_179_pe161_c: - DG_LMPD_179_pe161_c +1000 BU_LMPD_179_pe161_c < 1000
SU_LMPD_179_pe161_c: FU_LMPD_179_pe161_c + BU_LMPD_179_pe161_c < 1
UF_LMPD_179_pe161_c: F_LMPD_179_pe161_c -1000 FU_LMPD_179_pe161_c < 0
UR_LMPD_179_pe161_c: R_LMPD_179_pe161_c -1000 BU_LMPD_179_pe161_c < 0
G_LMPD_180_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +6.513 LC_nad_c -6.513 LC_nadh_c +10.6576 LC_nadp_c -10.6576 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_180_pe181_c + DGo_LMPD_180_pe181_c = 0
FU_LMPD_180_pe181_c: DG_LMPD_180_pe181_c +1000 FU_LMPD_180_pe181_c < 1000
BU_LMPD_180_pe181_c: - DG_LMPD_180_pe181_c +1000 BU_LMPD_180_pe181_c < 1000
SU_LMPD_180_pe181_c: FU_LMPD_180_pe181_c + BU_LMPD_180_pe181_c < 1
UF_LMPD_180_pe181_c: F_LMPD_180_pe181_c -1000 FU_LMPD_180_pe181_c < 0
UR_LMPD_180_pe181_c: R_LMPD_180_pe181_c -1000 BU_LMPD_180_pe181_c < 0
G_LMPD_181_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +7.6972 LC_nad_c -7.6972 LC_nadh_c +9.4734 LC_nadp_c -9.4734 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_181_pe181_c + DGo_LMPD_181_pe181_c = 0
FU_LMPD_181_pe181_c: DG_LMPD_181_pe181_c +1000 FU_LMPD_181_pe181_c < 1000
BU_LMPD_181_pe181_c: - DG_LMPD_181_pe181_c +1000 BU_LMPD_181_pe181_c < 1000
SU_LMPD_181_pe181_c: FU_LMPD_181_pe181_c + BU_LMPD_181_pe181_c < 1
UF_LMPD_181_pe181_c: F_LMPD_181_pe181_c -1000 FU_LMPD_181_pe181_c < 0
UR_LMPD_181_pe181_c: R_LMPD_181_pe181_c -1000 BU_LMPD_181_pe181_c < 0
G_LMPD_182_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +1.7763 LC_nad_c -1.7763 LC_nadh_c +15.3943 LC_nadp_c -15.3943 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_182_pe181_c + DGo_LMPD_182_pe181_c = 0
FU_LMPD_182_pe181_c: DG_LMPD_182_pe181_c +1000 FU_LMPD_182_pe181_c < 1000
BU_LMPD_182_pe181_c: - DG_LMPD_182_pe181_c +1000 BU_LMPD_182_pe181_c < 1000
SU_LMPD_182_pe181_c: FU_LMPD_182_pe181_c + BU_LMPD_182_pe181_c < 1
UF_LMPD_182_pe181_c: F_LMPD_182_pe181_c -1000 FU_LMPD_182_pe181_c < 0
UR_LMPD_182_pe181_c: R_LMPD_182_pe181_c -1000 BU_LMPD_182_pe181_c < 0
G_LMPD_183_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +2.9604 LC_nad_c -2.9604 LC_nadh_c +14.2101 LC_nadp_c -14.2101 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_183_pe181_c + DGo_LMPD_183_pe181_c = 0
FU_LMPD_183_pe181_c: DG_LMPD_183_pe181_c +1000 FU_LMPD_183_pe181_c < 1000
BU_LMPD_183_pe181_c: - DG_LMPD_183_pe181_c +1000 BU_LMPD_183_pe181_c < 1000
SU_LMPD_183_pe181_c: FU_LMPD_183_pe181_c + BU_LMPD_183_pe181_c < 1
UF_LMPD_183_pe181_c: F_LMPD_183_pe181_c -1000 FU_LMPD_183_pe181_c < 0
UR_LMPD_183_pe181_c: R_LMPD_183_pe181_c -1000 BU_LMPD_183_pe181_c < 0
G_LMPD_184_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +4.1446 LC_nad_c -4.1446 LC_nadh_c +13.026 LC_nadp_c -13.026 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_184_pe181_c + DGo_LMPD_184_pe181_c = 0
FU_LMPD_184_pe181_c: DG_LMPD_184_pe181_c +1000 FU_LMPD_184_pe181_c < 1000
BU_LMPD_184_pe181_c: - DG_LMPD_184_pe181_c +1000 BU_LMPD_184_pe181_c < 1000
SU_LMPD_184_pe181_c: FU_LMPD_184_pe181_c + BU_LMPD_184_pe181_c < 1
UF_LMPD_184_pe181_c: F_LMPD_184_pe181_c -1000 FU_LMPD_184_pe181_c < 0
UR_LMPD_184_pe181_c: R_LMPD_184_pe181_c -1000 BU_LMPD_184_pe181_c < 0
G_LMPD_185_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +5.3288 LC_nad_c -5.3288 LC_nadh_c +11.8418 LC_nadp_c -11.8418 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_185_pe181_c + DGo_LMPD_185_pe181_c = 0
FU_LMPD_185_pe181_c: DG_LMPD_185_pe181_c +1000 FU_LMPD_185_pe181_c < 1000
BU_LMPD_185_pe181_c: - DG_LMPD_185_pe181_c +1000 BU_LMPD_185_pe181_c < 1000
SU_LMPD_185_pe181_c: FU_LMPD_185_pe181_c + BU_LMPD_185_pe181_c < 1
UF_LMPD_185_pe181_c: F_LMPD_185_pe181_c -1000 FU_LMPD_185_pe181_c < 0
UR_LMPD_185_pe181_c: R_LMPD_185_pe181_c -1000 BU_LMPD_185_pe181_c < 0
G_LMPD_186_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +17.7627 LC_nadp_c -17.7627 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_186_pe181_c + DGo_LMPD_186_pe181_c = 0
FU_LMPD_186_pe181_c: DG_LMPD_186_pe181_c +1000 FU_LMPD_186_pe181_c < 1000
BU_LMPD_186_pe181_c: - DG_LMPD_186_pe181_c +1000 BU_LMPD_186_pe181_c < 1000
SU_LMPD_186_pe181_c: FU_LMPD_186_pe181_c + BU_LMPD_186_pe181_c < 1
UF_LMPD_186_pe181_c: F_LMPD_186_pe181_c -1000 FU_LMPD_186_pe181_c < 0
UR_LMPD_186_pe181_c: R_LMPD_186_pe181_c -1000 BU_LMPD_186_pe181_c < 0
G_LMPD_187_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +10.6576 LC_adp_c +0.59209 LC_akg_c -10.6576 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +0.59209 LC_nad_c -0.59209 LC_nadh_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +0.59209 LC_pe181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_187_pe181_c + DGo_LMPD_187_pe181_c = 0
FU_LMPD_187_pe181_c: DG_LMPD_187_pe181_c +1000 FU_LMPD_187_pe181_c < 1000
BU_LMPD_187_pe181_c: - DG_LMPD_187_pe181_c +1000 BU_LMPD_187_pe181_c < 1000
SU_LMPD_187_pe181_c: FU_LMPD_187_pe181_c + BU_LMPD_187_pe181_c < 1
UF_LMPD_187_pe181_c: F_LMPD_187_pe181_c -1000 FU_LMPD_187_pe181_c < 0
UR_LMPD_187_pe181_c: R_LMPD_187_pe181_c -1000 BU_LMPD_187_pe181_c < 0
G_LMPD_188_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -2.3684 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +8.8813 LC_nad_c -8.8813 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pe181_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_188_pe181_c + DGo_LMPD_188_pe181_c = 0
FU_LMPD_188_pe181_c: DG_LMPD_188_pe181_c +1000 FU_LMPD_188_pe181_c < 1000
BU_LMPD_188_pe181_c: - DG_LMPD_188_pe181_c +1000 BU_LMPD_188_pe181_c < 1000
SU_LMPD_188_pe181_c: FU_LMPD_188_pe181_c + BU_LMPD_188_pe181_c < 1
UF_LMPD_188_pe181_c: F_LMPD_188_pe181_c -1000 FU_LMPD_188_pe181_c < 0
UR_LMPD_188_pe181_c: R_LMPD_188_pe181_c -1000 BU_LMPD_188_pe181_c < 0
G_LMPD_189_pe181_c: -0.59209 LC_3pg_c -10.6576 LC_accoa_c +1.1842 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -2.3684 LC_atp_c +0.59209 LC_co2_c +10.6576 LC_coa_c -0.59209 LC_glu-L_c -0.59209 LC_glyc3p_c +17.1706 LC_nad_c -17.1706 LC_nadh_c +0.59209 LC_pe181_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_189_pe181_c + DGo_LMPD_189_pe181_c = 0
FU_LMPD_189_pe181_c: DG_LMPD_189_pe181_c +1000 FU_LMPD_189_pe181_c < 1000
BU_LMPD_189_pe181_c: - DG_LMPD_189_pe181_c +1000 BU_LMPD_189_pe181_c < 1000
SU_LMPD_189_pe181_c: FU_LMPD_189_pe181_c + BU_LMPD_189_pe181_c < 1
UF_LMPD_189_pe181_c: F_LMPD_189_pe181_c -1000 FU_LMPD_189_pe181_c < 0
UR_LMPD_189_pe181_c: R_LMPD_189_pe181_c -1000 BU_LMPD_189_pe181_c < 0
G_LMPD_190_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +5.9209 LC_nad_c -5.9209 LC_nadh_c +10.6576 LC_nadp_c -10.6576 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_190_pg160_c + DGo_LMPD_190_pg160_c = 0
FU_LMPD_190_pg160_c: DG_LMPD_190_pg160_c +1000 FU_LMPD_190_pg160_c < 1000
BU_LMPD_190_pg160_c: - DG_LMPD_190_pg160_c +1000 BU_LMPD_190_pg160_c < 1000
SU_LMPD_190_pg160_c: FU_LMPD_190_pg160_c + BU_LMPD_190_pg160_c < 1
UF_LMPD_190_pg160_c: F_LMPD_190_pg160_c -1000 FU_LMPD_190_pg160_c < 0
UR_LMPD_190_pg160_c: R_LMPD_190_pg160_c -1000 BU_LMPD_190_pg160_c < 0
G_LMPD_191_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +7.1051 LC_nad_c -7.1051 LC_nadh_c +9.4734 LC_nadp_c -9.4734 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_191_pg160_c + DGo_LMPD_191_pg160_c = 0
FU_LMPD_191_pg160_c: DG_LMPD_191_pg160_c +1000 FU_LMPD_191_pg160_c < 1000
BU_LMPD_191_pg160_c: - DG_LMPD_191_pg160_c +1000 BU_LMPD_191_pg160_c < 1000
SU_LMPD_191_pg160_c: FU_LMPD_191_pg160_c + BU_LMPD_191_pg160_c < 1
UF_LMPD_191_pg160_c: F_LMPD_191_pg160_c -1000 FU_LMPD_191_pg160_c < 0
UR_LMPD_191_pg160_c: R_LMPD_191_pg160_c -1000 BU_LMPD_191_pg160_c < 0
G_LMPD_192_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +8.2893 LC_nad_c -8.2893 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_192_pg160_c + DGo_LMPD_192_pg160_c = 0
FU_LMPD_192_pg160_c: DG_LMPD_192_pg160_c +1000 FU_LMPD_192_pg160_c < 1000
BU_LMPD_192_pg160_c: - DG_LMPD_192_pg160_c +1000 BU_LMPD_192_pg160_c < 1000
SU_LMPD_192_pg160_c: FU_LMPD_192_pg160_c + BU_LMPD_192_pg160_c < 1
UF_LMPD_192_pg160_c: F_LMPD_192_pg160_c -1000 FU_LMPD_192_pg160_c < 0
UR_LMPD_192_pg160_c: R_LMPD_192_pg160_c -1000 BU_LMPD_192_pg160_c < 0
G_LMPD_193_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +1.1842 LC_nad_c -1.1842 LC_nadh_c +15.3943 LC_nadp_c -15.3943 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_193_pg160_c + DGo_LMPD_193_pg160_c = 0
FU_LMPD_193_pg160_c: DG_LMPD_193_pg160_c +1000 FU_LMPD_193_pg160_c < 1000
BU_LMPD_193_pg160_c: - DG_LMPD_193_pg160_c +1000 BU_LMPD_193_pg160_c < 1000
SU_LMPD_193_pg160_c: FU_LMPD_193_pg160_c + BU_LMPD_193_pg160_c < 1
UF_LMPD_193_pg160_c: F_LMPD_193_pg160_c -1000 FU_LMPD_193_pg160_c < 0
UR_LMPD_193_pg160_c: R_LMPD_193_pg160_c -1000 BU_LMPD_193_pg160_c < 0
G_LMPD_194_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_194_pg160_c + DGo_LMPD_194_pg160_c = 0
FU_LMPD_194_pg160_c: DG_LMPD_194_pg160_c +1000 FU_LMPD_194_pg160_c < 1000
BU_LMPD_194_pg160_c: - DG_LMPD_194_pg160_c +1000 BU_LMPD_194_pg160_c < 1000
SU_LMPD_194_pg160_c: FU_LMPD_194_pg160_c + BU_LMPD_194_pg160_c < 1
UF_LMPD_194_pg160_c: F_LMPD_194_pg160_c -1000 FU_LMPD_194_pg160_c < 0
UR_LMPD_194_pg160_c: R_LMPD_194_pg160_c -1000 BU_LMPD_194_pg160_c < 0
G_LMPD_195_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +2.3684 LC_nad_c -2.3684 LC_nadh_c +14.2101 LC_nadp_c -14.2101 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_195_pg160_c + DGo_LMPD_195_pg160_c = 0
FU_LMPD_195_pg160_c: DG_LMPD_195_pg160_c +1000 FU_LMPD_195_pg160_c < 1000
BU_LMPD_195_pg160_c: - DG_LMPD_195_pg160_c +1000 BU_LMPD_195_pg160_c < 1000
SU_LMPD_195_pg160_c: FU_LMPD_195_pg160_c + BU_LMPD_195_pg160_c < 1
UF_LMPD_195_pg160_c: F_LMPD_195_pg160_c -1000 FU_LMPD_195_pg160_c < 0
UR_LMPD_195_pg160_c: R_LMPD_195_pg160_c -1000 BU_LMPD_195_pg160_c < 0
G_LMPD_196_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +3.5525 LC_nad_c -3.5525 LC_nadh_c +13.026 LC_nadp_c -13.026 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_196_pg160_c + DGo_LMPD_196_pg160_c = 0
FU_LMPD_196_pg160_c: DG_LMPD_196_pg160_c +1000 FU_LMPD_196_pg160_c < 1000
BU_LMPD_196_pg160_c: - DG_LMPD_196_pg160_c +1000 BU_LMPD_196_pg160_c < 1000
SU_LMPD_196_pg160_c: FU_LMPD_196_pg160_c + BU_LMPD_196_pg160_c < 1
UF_LMPD_196_pg160_c: F_LMPD_196_pg160_c -1000 FU_LMPD_196_pg160_c < 0
UR_LMPD_196_pg160_c: R_LMPD_196_pg160_c -1000 BU_LMPD_196_pg160_c < 0
G_LMPD_197_pg160_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +4.7367 LC_nad_c -4.7367 LC_nadh_c +11.8418 LC_nadp_c -11.8418 LC_nadph_c +0.59209 LC_pg160_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_197_pg160_c + DGo_LMPD_197_pg160_c = 0
FU_LMPD_197_pg160_c: DG_LMPD_197_pg160_c +1000 FU_LMPD_197_pg160_c < 1000
BU_LMPD_197_pg160_c: - DG_LMPD_197_pg160_c +1000 BU_LMPD_197_pg160_c < 1000
SU_LMPD_197_pg160_c: FU_LMPD_197_pg160_c + BU_LMPD_197_pg160_c < 1
UF_LMPD_197_pg160_c: F_LMPD_197_pg160_c -1000 FU_LMPD_197_pg160_c < 0
UR_LMPD_197_pg160_c: R_LMPD_197_pg160_c -1000 BU_LMPD_197_pg160_c < 0
G_LMPD_198_pg160_c: -9.4734 LC_accoa_c +1.1842 LC_adp_c +1.1842 LC_amp_c -2.3684 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +8.2893 LC_nad_c -8.2893 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pg160_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_198_pg160_c + DGo_LMPD_198_pg160_c = 0
FU_LMPD_198_pg160_c: DG_LMPD_198_pg160_c +1000 FU_LMPD_198_pg160_c < 1000
BU_LMPD_198_pg160_c: - DG_LMPD_198_pg160_c +1000 BU_LMPD_198_pg160_c < 1000
SU_LMPD_198_pg160_c: FU_LMPD_198_pg160_c + BU_LMPD_198_pg160_c < 1
UF_LMPD_198_pg160_c: F_LMPD_198_pg160_c -1000 FU_LMPD_198_pg160_c < 0
UR_LMPD_198_pg160_c: R_LMPD_198_pg160_c -1000 BU_LMPD_198_pg160_c < 0
G_LMPD_199_pg160_c: -9.4734 LC_accoa_c +1.1842 LC_adp_c +1.1842 LC_amp_c -2.3684 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +16.5785 LC_nad_c -16.5785 LC_nadh_c +0.59209 LC_pg160_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_199_pg160_c + DGo_LMPD_199_pg160_c = 0
FU_LMPD_199_pg160_c: DG_LMPD_199_pg160_c +1000 FU_LMPD_199_pg160_c < 1000
BU_LMPD_199_pg160_c: - DG_LMPD_199_pg160_c +1000 BU_LMPD_199_pg160_c < 1000
SU_LMPD_199_pg160_c: FU_LMPD_199_pg160_c + BU_LMPD_199_pg160_c < 1
UF_LMPD_199_pg160_c: F_LMPD_199_pg160_c -1000 FU_LMPD_199_pg160_c < 0
UR_LMPD_199_pg160_c: R_LMPD_199_pg160_c -1000 BU_LMPD_199_pg160_c < 0
G_LMPD_200_pg161_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +5.9209 LC_nad_c -5.9209 LC_nadh_c +9.4734 LC_nadp_c -9.4734 LC_nadph_c +0.59209 LC_pg161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_200_pg161_c + DGo_LMPD_200_pg161_c = 0
FU_LMPD_200_pg161_c: DG_LMPD_200_pg161_c +1000 FU_LMPD_200_pg161_c < 1000
BU_LMPD_200_pg161_c: - DG_LMPD_200_pg161_c +1000 BU_LMPD_200_pg161_c < 1000
SU_LMPD_200_pg161_c: FU_LMPD_200_pg161_c + BU_LMPD_200_pg161_c < 1
UF_LMPD_200_pg161_c: F_LMPD_200_pg161_c -1000 FU_LMPD_200_pg161_c < 0
UR_LMPD_200_pg161_c: R_LMPD_200_pg161_c -1000 BU_LMPD_200_pg161_c < 0
G_LMPD_201_pg161_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +7.1051 LC_nad_c -7.1051 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pg161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_201_pg161_c + DGo_LMPD_201_pg161_c = 0
FU_LMPD_201_pg161_c: DG_LMPD_201_pg161_c +1000 FU_LMPD_201_pg161_c < 1000
BU_LMPD_201_pg161_c: - DG_LMPD_201_pg161_c +1000 BU_LMPD_201_pg161_c < 1000
SU_LMPD_201_pg161_c: FU_LMPD_201_pg161_c + BU_LMPD_201_pg161_c < 1
UF_LMPD_201_pg161_c: F_LMPD_201_pg161_c -1000 FU_LMPD_201_pg161_c < 0
UR_LMPD_201_pg161_c: R_LMPD_201_pg161_c -1000 BU_LMPD_201_pg161_c < 0
G_LMPD_202_pg161_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +1.1842 LC_nad_c -1.1842 LC_nadh_c +14.2101 LC_nadp_c -14.2101 LC_nadph_c +0.59209 LC_pg161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_202_pg161_c + DGo_LMPD_202_pg161_c = 0
FU_LMPD_202_pg161_c: DG_LMPD_202_pg161_c +1000 FU_LMPD_202_pg161_c < 1000
BU_LMPD_202_pg161_c: - DG_LMPD_202_pg161_c +1000 BU_LMPD_202_pg161_c < 1000
SU_LMPD_202_pg161_c: FU_LMPD_202_pg161_c + BU_LMPD_202_pg161_c < 1
UF_LMPD_202_pg161_c: F_LMPD_202_pg161_c -1000 FU_LMPD_202_pg161_c < 0
UR_LMPD_202_pg161_c: R_LMPD_202_pg161_c -1000 BU_LMPD_202_pg161_c < 0
G_LMPD_203_pg161_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +15.3943 LC_nadp_c -15.3943 LC_nadph_c +0.59209 LC_pg161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_203_pg161_c + DGo_LMPD_203_pg161_c = 0
FU_LMPD_203_pg161_c: DG_LMPD_203_pg161_c +1000 FU_LMPD_203_pg161_c < 1000
BU_LMPD_203_pg161_c: - DG_LMPD_203_pg161_c +1000 BU_LMPD_203_pg161_c < 1000
SU_LMPD_203_pg161_c: FU_LMPD_203_pg161_c + BU_LMPD_203_pg161_c < 1
UF_LMPD_203_pg161_c: F_LMPD_203_pg161_c -1000 FU_LMPD_203_pg161_c < 0
UR_LMPD_203_pg161_c: R_LMPD_203_pg161_c -1000 BU_LMPD_203_pg161_c < 0
G_LMPD_204_pg161_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +2.3684 LC_nad_c -2.3684 LC_nadh_c +13.026 LC_nadp_c -13.026 LC_nadph_c +0.59209 LC_pg161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_204_pg161_c + DGo_LMPD_204_pg161_c = 0
FU_LMPD_204_pg161_c: DG_LMPD_204_pg161_c +1000 FU_LMPD_204_pg161_c < 1000
BU_LMPD_204_pg161_c: - DG_LMPD_204_pg161_c +1000 BU_LMPD_204_pg161_c < 1000
SU_LMPD_204_pg161_c: FU_LMPD_204_pg161_c + BU_LMPD_204_pg161_c < 1
UF_LMPD_204_pg161_c: F_LMPD_204_pg161_c -1000 FU_LMPD_204_pg161_c < 0
UR_LMPD_204_pg161_c: R_LMPD_204_pg161_c -1000 BU_LMPD_204_pg161_c < 0
G_LMPD_205_pg161_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +3.5525 LC_nad_c -3.5525 LC_nadh_c +11.8418 LC_nadp_c -11.8418 LC_nadph_c +0.59209 LC_pg161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_205_pg161_c + DGo_LMPD_205_pg161_c = 0
FU_LMPD_205_pg161_c: DG_LMPD_205_pg161_c +1000 FU_LMPD_205_pg161_c < 1000
BU_LMPD_205_pg161_c: - DG_LMPD_205_pg161_c +1000 BU_LMPD_205_pg161_c < 1000
SU_LMPD_205_pg161_c: FU_LMPD_205_pg161_c + BU_LMPD_205_pg161_c < 1
UF_LMPD_205_pg161_c: F_LMPD_205_pg161_c -1000 FU_LMPD_205_pg161_c < 0
UR_LMPD_205_pg161_c: R_LMPD_205_pg161_c -1000 BU_LMPD_205_pg161_c < 0
G_LMPD_206_pg161_c: -9.4734 LC_accoa_c +9.4734 LC_adp_c -9.4734 LC_atp_c +9.4734 LC_coa_c -1.1842 LC_glyc3p_c +4.7367 LC_nad_c -4.7367 LC_nadh_c +10.6576 LC_nadp_c -10.6576 LC_nadph_c +0.59209 LC_pg161_c +8.8813 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_206_pg161_c + DGo_LMPD_206_pg161_c = 0
FU_LMPD_206_pg161_c: DG_LMPD_206_pg161_c +1000 FU_LMPD_206_pg161_c < 1000
BU_LMPD_206_pg161_c: - DG_LMPD_206_pg161_c +1000 BU_LMPD_206_pg161_c < 1000
SU_LMPD_206_pg161_c: FU_LMPD_206_pg161_c + BU_LMPD_206_pg161_c < 1
UF_LMPD_206_pg161_c: F_LMPD_206_pg161_c -1000 FU_LMPD_206_pg161_c < 0
UR_LMPD_206_pg161_c: R_LMPD_206_pg161_c -1000 BU_LMPD_206_pg161_c < 0
G_LMPD_207_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +5.9209 LC_nad_c -5.9209 LC_nadh_c +11.8418 LC_nadp_c -11.8418 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_207_pg181_c + DGo_LMPD_207_pg181_c = 0
FU_LMPD_207_pg181_c: DG_LMPD_207_pg181_c +1000 FU_LMPD_207_pg181_c < 1000
BU_LMPD_207_pg181_c: - DG_LMPD_207_pg181_c +1000 BU_LMPD_207_pg181_c < 1000
SU_LMPD_207_pg181_c: FU_LMPD_207_pg181_c + BU_LMPD_207_pg181_c < 1
UF_LMPD_207_pg181_c: F_LMPD_207_pg181_c -1000 FU_LMPD_207_pg181_c < 0
UR_LMPD_207_pg181_c: R_LMPD_207_pg181_c -1000 BU_LMPD_207_pg181_c < 0
G_LMPD_208_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +7.1051 LC_nad_c -7.1051 LC_nadh_c +10.6576 LC_nadp_c -10.6576 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_208_pg181_c + DGo_LMPD_208_pg181_c = 0
FU_LMPD_208_pg181_c: DG_LMPD_208_pg181_c +1000 FU_LMPD_208_pg181_c < 1000
BU_LMPD_208_pg181_c: - DG_LMPD_208_pg181_c +1000 BU_LMPD_208_pg181_c < 1000
SU_LMPD_208_pg181_c: FU_LMPD_208_pg181_c + BU_LMPD_208_pg181_c < 1
UF_LMPD_208_pg181_c: F_LMPD_208_pg181_c -1000 FU_LMPD_208_pg181_c < 0
UR_LMPD_208_pg181_c: R_LMPD_208_pg181_c -1000 BU_LMPD_208_pg181_c < 0
G_LMPD_209_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +8.2893 LC_nad_c -8.2893 LC_nadh_c +9.4734 LC_nadp_c -9.4734 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_209_pg181_c + DGo_LMPD_209_pg181_c = 0
FU_LMPD_209_pg181_c: DG_LMPD_209_pg181_c +1000 FU_LMPD_209_pg181_c < 1000
BU_LMPD_209_pg181_c: - DG_LMPD_209_pg181_c +1000 BU_LMPD_209_pg181_c < 1000
SU_LMPD_209_pg181_c: FU_LMPD_209_pg181_c + BU_LMPD_209_pg181_c < 1
UF_LMPD_209_pg181_c: F_LMPD_209_pg181_c -1000 FU_LMPD_209_pg181_c < 0
UR_LMPD_209_pg181_c: R_LMPD_209_pg181_c -1000 BU_LMPD_209_pg181_c < 0
G_LMPD_210_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +1.1842 LC_nad_c -1.1842 LC_nadh_c +16.5785 LC_nadp_c -16.5785 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_210_pg181_c + DGo_LMPD_210_pg181_c = 0
FU_LMPD_210_pg181_c: DG_LMPD_210_pg181_c +1000 FU_LMPD_210_pg181_c < 1000
BU_LMPD_210_pg181_c: - DG_LMPD_210_pg181_c +1000 BU_LMPD_210_pg181_c < 1000
SU_LMPD_210_pg181_c: FU_LMPD_210_pg181_c + BU_LMPD_210_pg181_c < 1
UF_LMPD_210_pg181_c: F_LMPD_210_pg181_c -1000 FU_LMPD_210_pg181_c < 0
UR_LMPD_210_pg181_c: R_LMPD_210_pg181_c -1000 BU_LMPD_210_pg181_c < 0
G_LMPD_211_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +17.7627 LC_nadp_c -17.7627 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_211_pg181_c + DGo_LMPD_211_pg181_c = 0
FU_LMPD_211_pg181_c: DG_LMPD_211_pg181_c +1000 FU_LMPD_211_pg181_c < 1000
BU_LMPD_211_pg181_c: - DG_LMPD_211_pg181_c +1000 BU_LMPD_211_pg181_c < 1000
SU_LMPD_211_pg181_c: FU_LMPD_211_pg181_c + BU_LMPD_211_pg181_c < 1
UF_LMPD_211_pg181_c: F_LMPD_211_pg181_c -1000 FU_LMPD_211_pg181_c < 0
UR_LMPD_211_pg181_c: R_LMPD_211_pg181_c -1000 BU_LMPD_211_pg181_c < 0
G_LMPD_212_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +2.3684 LC_nad_c -2.3684 LC_nadh_c +15.3943 LC_nadp_c -15.3943 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_212_pg181_c + DGo_LMPD_212_pg181_c = 0
FU_LMPD_212_pg181_c: DG_LMPD_212_pg181_c +1000 FU_LMPD_212_pg181_c < 1000
BU_LMPD_212_pg181_c: - DG_LMPD_212_pg181_c +1000 BU_LMPD_212_pg181_c < 1000
SU_LMPD_212_pg181_c: FU_LMPD_212_pg181_c + BU_LMPD_212_pg181_c < 1
UF_LMPD_212_pg181_c: F_LMPD_212_pg181_c -1000 FU_LMPD_212_pg181_c < 0
UR_LMPD_212_pg181_c: R_LMPD_212_pg181_c -1000 BU_LMPD_212_pg181_c < 0
G_LMPD_213_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +3.5525 LC_nad_c -3.5525 LC_nadh_c +14.2101 LC_nadp_c -14.2101 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_213_pg181_c + DGo_LMPD_213_pg181_c = 0
FU_LMPD_213_pg181_c: DG_LMPD_213_pg181_c +1000 FU_LMPD_213_pg181_c < 1000
BU_LMPD_213_pg181_c: - DG_LMPD_213_pg181_c +1000 BU_LMPD_213_pg181_c < 1000
SU_LMPD_213_pg181_c: FU_LMPD_213_pg181_c + BU_LMPD_213_pg181_c < 1
UF_LMPD_213_pg181_c: F_LMPD_213_pg181_c -1000 FU_LMPD_213_pg181_c < 0
UR_LMPD_213_pg181_c: R_LMPD_213_pg181_c -1000 BU_LMPD_213_pg181_c < 0
G_LMPD_214_pg181_c: -10.6576 LC_accoa_c +10.6576 LC_adp_c -10.6576 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +4.7367 LC_nad_c -4.7367 LC_nadh_c +13.026 LC_nadp_c -13.026 LC_nadph_c +0.59209 LC_pg181_c +10.0655 LC_pi_c +0.59209 LC_ppi_c - DG_LMPD_214_pg181_c + DGo_LMPD_214_pg181_c = 0
FU_LMPD_214_pg181_c: DG_LMPD_214_pg181_c +1000 FU_LMPD_214_pg181_c < 1000
BU_LMPD_214_pg181_c: - DG_LMPD_214_pg181_c +1000 BU_LMPD_214_pg181_c < 1000
SU_LMPD_214_pg181_c: FU_LMPD_214_pg181_c + BU_LMPD_214_pg181_c < 1
UF_LMPD_214_pg181_c: F_LMPD_214_pg181_c -1000 FU_LMPD_214_pg181_c < 0
UR_LMPD_214_pg181_c: R_LMPD_214_pg181_c -1000 BU_LMPD_214_pg181_c < 0
G_LMPD_215_pg181_c: -10.6576 LC_accoa_c +1.1842 LC_adp_c +1.1842 LC_amp_c -2.3684 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +9.4734 LC_nad_c -9.4734 LC_nadh_c +8.2893 LC_nadp_c -8.2893 LC_nadph_c +0.59209 LC_pg181_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_215_pg181_c + DGo_LMPD_215_pg181_c = 0
FU_LMPD_215_pg181_c: DG_LMPD_215_pg181_c +1000 FU_LMPD_215_pg181_c < 1000
BU_LMPD_215_pg181_c: - DG_LMPD_215_pg181_c +1000 BU_LMPD_215_pg181_c < 1000
SU_LMPD_215_pg181_c: FU_LMPD_215_pg181_c + BU_LMPD_215_pg181_c < 1
UF_LMPD_215_pg181_c: F_LMPD_215_pg181_c -1000 FU_LMPD_215_pg181_c < 0
UR_LMPD_215_pg181_c: R_LMPD_215_pg181_c -1000 BU_LMPD_215_pg181_c < 0
G_LMPD_216_pg181_c: -10.6576 LC_accoa_c +1.1842 LC_adp_c +1.1842 LC_amp_c -2.3684 LC_atp_c +10.6576 LC_coa_c -1.1842 LC_glyc3p_c +17.7627 LC_nad_c -17.7627 LC_nadh_c +0.59209 LC_pg181_c +0.59209 LC_pi_c +1.7763 LC_ppi_c - DG_LMPD_216_pg181_c + DGo_LMPD_216_pg181_c = 0
FU_LMPD_216_pg181_c: DG_LMPD_216_pg181_c +1000 FU_LMPD_216_pg181_c < 1000
BU_LMPD_216_pg181_c: - DG_LMPD_216_pg181_c +1000 BU_LMPD_216_pg181_c < 1000
SU_LMPD_216_pg181_c: FU_LMPD_216_pg181_c + BU_LMPD_216_pg181_c < 1
UF_LMPD_216_pg181_c: F_LMPD_216_pg181_c -1000 FU_LMPD_216_pg181_c < 0
UR_LMPD_216_pg181_c: R_LMPD_216_pg181_c -1000 BU_LMPD_216_pg181_c < 0
G_LMPD_217_phe-L_c: +0.59209 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_glu-L_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_pep_c +0.59209 LC_phe-L_c +2.3684 LC_pi_c - DG_LMPD_217_phe-L_c + DGo_LMPD_217_phe-L_c = 0
FU_LMPD_217_phe-L_c: DG_LMPD_217_phe-L_c +1000 FU_LMPD_217_phe-L_c < 1000
BU_LMPD_217_phe-L_c: - DG_LMPD_217_phe-L_c +1000 BU_LMPD_217_phe-L_c < 1000
SU_LMPD_217_phe-L_c: FU_LMPD_217_phe-L_c + BU_LMPD_217_phe-L_c < 1
UF_LMPD_217_phe-L_c: F_LMPD_217_phe-L_c -1000 FU_LMPD_217_phe-L_c < 0
UR_LMPD_217_phe-L_c: R_LMPD_217_phe-L_c -1000 BU_LMPD_217_phe-L_c < 0
SU_LMPD_218_pheme_c: FU_LMPD_218_pheme_c + BU_LMPD_218_pheme_c < 1
UF_LMPD_218_pheme_c: F_LMPD_218_pheme_c -1000 FU_LMPD_218_pheme_c < 0
UR_LMPD_218_pheme_c: R_LMPD_218_pheme_c -1000 BU_LMPD_218_pheme_c < 0
SU_LMPD_219_pheme_c: FU_LMPD_219_pheme_c + BU_LMPD_219_pheme_c < 1
UF_LMPD_219_pheme_c: F_LMPD_219_pheme_c -1000 FU_LMPD_219_pheme_c < 0
UR_LMPD_219_pheme_c: R_LMPD_219_pheme_c -1000 BU_LMPD_219_pheme_c < 0
G_LMPD_220_pro-L_c: +0.59209 LC_adp_c -0.59209 LC_atp_c -0.59209 LC_glu-L_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +0.59209 LC_pi_c +0.59209 LC_pro-L_c - DG_LMPD_220_pro-L_c + DGo_LMPD_220_pro-L_c = 0
FU_LMPD_220_pro-L_c: DG_LMPD_220_pro-L_c +1000 FU_LMPD_220_pro-L_c < 1000
BU_LMPD_220_pro-L_c: - DG_LMPD_220_pro-L_c +1000 BU_LMPD_220_pro-L_c < 1000
SU_LMPD_220_pro-L_c: FU_LMPD_220_pro-L_c + BU_LMPD_220_pro-L_c < 1
UF_LMPD_220_pro-L_c: F_LMPD_220_pro-L_c -1000 FU_LMPD_220_pro-L_c < 0
UR_LMPD_220_pro-L_c: R_LMPD_220_pro-L_c -1000 BU_LMPD_220_pro-L_c < 0
G_LMPD_221_ptrc_c: +0.59209 LC_ac_c -0.59209 LC_accoa_c +0.59209 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_atp_c +0.59209 LC_co2_c +0.59209 LC_coa_c -1.1842 LC_glu-L_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c +0.59209 LC_pi_c +0.59209 LC_ptrc_c - DG_LMPD_221_ptrc_c + DGo_LMPD_221_ptrc_c = 0
FU_LMPD_221_ptrc_c: DG_LMPD_221_ptrc_c +1000 FU_LMPD_221_ptrc_c < 1000
BU_LMPD_221_ptrc_c: - DG_LMPD_221_ptrc_c +1000 BU_LMPD_221_ptrc_c < 1000
SU_LMPD_221_ptrc_c: FU_LMPD_221_ptrc_c + BU_LMPD_221_ptrc_c < 1
UF_LMPD_221_ptrc_c: F_LMPD_221_ptrc_c -1000 FU_LMPD_221_ptrc_c < 0
UR_LMPD_221_ptrc_c: R_LMPD_221_ptrc_c -1000 BU_LMPD_221_ptrc_c < 0
G_LMPD_222_pydx5p_c: +0.59209 LC_akg_c +1.1842 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_g3p_c -0.59209 LC_glu-L_c +0.59209 LC_h2o2_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_o2_c +0.59209 LC_pi_c +0.59209 LC_pydx5p_c -0.59209 LC_pyr_c - DG_LMPD_222_pydx5p_c + DGo_LMPD_222_pydx5p_c = 0
FU_LMPD_222_pydx5p_c: DG_LMPD_222_pydx5p_c +1000 FU_LMPD_222_pydx5p_c < 1000
BU_LMPD_222_pydx5p_c: - DG_LMPD_222_pydx5p_c +1000 BU_LMPD_222_pydx5p_c < 1000
SU_LMPD_222_pydx5p_c: FU_LMPD_222_pydx5p_c + BU_LMPD_222_pydx5p_c < 1
UF_LMPD_222_pydx5p_c: F_LMPD_222_pydx5p_c -1000 FU_LMPD_222_pydx5p_c < 0
UR_LMPD_222_pydx5p_c: R_LMPD_222_pydx5p_c -1000 BU_LMPD_222_pydx5p_c < 0
G_LMPD_223_pydx5p_c: +0.59209 LC_akg_c +1.1842 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_g3p_c -0.59209 LC_glu-L_c -2.3684 LC_nad_c +2.3684 LC_nadh_c +0.59209 LC_pi_c +0.59209 LC_pydx5p_c -0.59209 LC_pyr_c - DG_LMPD_223_pydx5p_c + DGo_LMPD_223_pydx5p_c = 0
FU_LMPD_223_pydx5p_c: DG_LMPD_223_pydx5p_c +1000 FU_LMPD_223_pydx5p_c < 1000
BU_LMPD_223_pydx5p_c: - DG_LMPD_223_pydx5p_c +1000 BU_LMPD_223_pydx5p_c < 1000
SU_LMPD_223_pydx5p_c: FU_LMPD_223_pydx5p_c + BU_LMPD_223_pydx5p_c < 1
UF_LMPD_223_pydx5p_c: F_LMPD_223_pydx5p_c -1000 FU_LMPD_223_pydx5p_c < 0
UR_LMPD_223_pydx5p_c: R_LMPD_223_pydx5p_c -1000 BU_LMPD_223_pydx5p_c < 0
SU_LMPD_224_q8h2_c: FU_LMPD_224_q8h2_c + BU_LMPD_224_q8h2_c < 1
UF_LMPD_224_q8h2_c: F_LMPD_224_q8h2_c -1000 FU_LMPD_224_q8h2_c < 0
UR_LMPD_224_q8h2_c: R_LMPD_224_q8h2_c -1000 BU_LMPD_224_q8h2_c < 0
SU_LMPD_225_q8h2_c: FU_LMPD_225_q8h2_c + BU_LMPD_225_q8h2_c < 1
UF_LMPD_225_q8h2_c: F_LMPD_225_q8h2_c -1000 FU_LMPD_225_q8h2_c < 0
UR_LMPD_225_q8h2_c: R_LMPD_225_q8h2_c -1000 BU_LMPD_225_q8h2_c < 0
SU_LMPD_226_q8h2_c: FU_LMPD_226_q8h2_c + BU_LMPD_226_q8h2_c < 1
UF_LMPD_226_q8h2_c: F_LMPD_226_q8h2_c -1000 FU_LMPD_226_q8h2_c < 0
UR_LMPD_226_q8h2_c: R_LMPD_226_q8h2_c -1000 BU_LMPD_226_q8h2_c < 0
SU_LMPD_227_q8h2_c: FU_LMPD_227_q8h2_c + BU_LMPD_227_q8h2_c < 1
UF_LMPD_227_q8h2_c: F_LMPD_227_q8h2_c -1000 FU_LMPD_227_q8h2_c < 0
UR_LMPD_227_q8h2_c: R_LMPD_227_q8h2_c -1000 BU_LMPD_227_q8h2_c < 0
SU_LMPD_228_q8h2_c: FU_LMPD_228_q8h2_c + BU_LMPD_228_q8h2_c < 1
UF_LMPD_228_q8h2_c: F_LMPD_228_q8h2_c -1000 FU_LMPD_228_q8h2_c < 0
UR_LMPD_228_q8h2_c: R_LMPD_228_q8h2_c -1000 BU_LMPD_228_q8h2_c < 0
SU_LMPD_229_q8h2_c: FU_LMPD_229_q8h2_c + BU_LMPD_229_q8h2_c < 1
UF_LMPD_229_q8h2_c: F_LMPD_229_q8h2_c -1000 FU_LMPD_229_q8h2_c < 0
UR_LMPD_229_q8h2_c: R_LMPD_229_q8h2_c -1000 BU_LMPD_229_q8h2_c < 0
SU_LMPD_230_q8h2_c: FU_LMPD_230_q8h2_c + BU_LMPD_230_q8h2_c < 1
UF_LMPD_230_q8h2_c: F_LMPD_230_q8h2_c -1000 FU_LMPD_230_q8h2_c < 0
UR_LMPD_230_q8h2_c: R_LMPD_230_q8h2_c -1000 BU_LMPD_230_q8h2_c < 0
SU_LMPD_231_q8h2_c: FU_LMPD_231_q8h2_c + BU_LMPD_231_q8h2_c < 1
UF_LMPD_231_q8h2_c: F_LMPD_231_q8h2_c -1000 FU_LMPD_231_q8h2_c < 0
UR_LMPD_231_q8h2_c: R_LMPD_231_q8h2_c -1000 BU_LMPD_231_q8h2_c < 0
SU_LMPD_232_q8h2_c: FU_LMPD_232_q8h2_c + BU_LMPD_232_q8h2_c < 1
UF_LMPD_232_q8h2_c: F_LMPD_232_q8h2_c -1000 FU_LMPD_232_q8h2_c < 0
UR_LMPD_232_q8h2_c: R_LMPD_232_q8h2_c -1000 BU_LMPD_232_q8h2_c < 0
SU_LMPD_233_q8h2_c: FU_LMPD_233_q8h2_c + BU_LMPD_233_q8h2_c < 1
UF_LMPD_233_q8h2_c: F_LMPD_233_q8h2_c -1000 FU_LMPD_233_q8h2_c < 0
UR_LMPD_233_q8h2_c: R_LMPD_233_q8h2_c -1000 BU_LMPD_233_q8h2_c < 0
SU_LMPD_234_q8h2_c: FU_LMPD_234_q8h2_c + BU_LMPD_234_q8h2_c < 1
UF_LMPD_234_q8h2_c: F_LMPD_234_q8h2_c -1000 FU_LMPD_234_q8h2_c < 0
UR_LMPD_234_q8h2_c: R_LMPD_234_q8h2_c -1000 BU_LMPD_234_q8h2_c < 0
G_LMPD_235_ribflv_c: -1.1842 LC_3pg_c +4.7367 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -6.513 LC_atp_c -0.59209 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_for_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.1842 LC_nh4_c +7.6972 LC_pi_c +1.7763 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_ribflv_c - DG_LMPD_235_ribflv_c + DGo_LMPD_235_ribflv_c = 0
FU_LMPD_235_ribflv_c: DG_LMPD_235_ribflv_c +1000 FU_LMPD_235_ribflv_c < 1000
BU_LMPD_235_ribflv_c: - DG_LMPD_235_ribflv_c +1000 BU_LMPD_235_ribflv_c < 1000
SU_LMPD_235_ribflv_c: FU_LMPD_235_ribflv_c + BU_LMPD_235_ribflv_c < 1
UF_LMPD_235_ribflv_c: F_LMPD_235_ribflv_c -1000 FU_LMPD_235_ribflv_c < 0
UR_LMPD_235_ribflv_c: R_LMPD_235_ribflv_c -1000 BU_LMPD_235_ribflv_c < 0
G_LMPD_236_ribflv_c: -1.1842 LC_3pg_c +5.9209 LC_adp_c +1.1842 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_for_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_gly_c -1.7763 LC_nad_c +1.7763 LC_nadh_c -0.59209 LC_nadp_c +0.59209 LC_nadph_c -1.1842 LC_nh4_c +7.6972 LC_pi_c +1.7763 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_ribflv_c - DG_LMPD_236_ribflv_c + DGo_LMPD_236_ribflv_c = 0
FU_LMPD_236_ribflv_c: DG_LMPD_236_ribflv_c +1000 FU_LMPD_236_ribflv_c < 1000
BU_LMPD_236_ribflv_c: - DG_LMPD_236_ribflv_c +1000 BU_LMPD_236_ribflv_c < 1000
SU_LMPD_236_ribflv_c: FU_LMPD_236_ribflv_c + BU_LMPD_236_ribflv_c < 1
UF_LMPD_236_ribflv_c: F_LMPD_236_ribflv_c -1000 FU_LMPD_236_ribflv_c < 0
UR_LMPD_236_ribflv_c: R_LMPD_236_ribflv_c -1000 BU_LMPD_236_ribflv_c < 0
G_LMPD_237_ribflv_c: -0.59209 LC_3pg_c +5.3288 LC_adp_c +0.59209 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -7.1051 LC_atp_c -0.59209 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_fum_c -0.59209 LC_glu-L_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nh4_c +7.6972 LC_pi_c +1.7763 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_ribflv_c - DG_LMPD_237_ribflv_c + DGo_LMPD_237_ribflv_c = 0
FU_LMPD_237_ribflv_c: DG_LMPD_237_ribflv_c +1000 FU_LMPD_237_ribflv_c < 1000
BU_LMPD_237_ribflv_c: - DG_LMPD_237_ribflv_c +1000 BU_LMPD_237_ribflv_c < 1000
SU_LMPD_237_ribflv_c: FU_LMPD_237_ribflv_c + BU_LMPD_237_ribflv_c < 1
UF_LMPD_237_ribflv_c: F_LMPD_237_ribflv_c -1000 FU_LMPD_237_ribflv_c < 0
UR_LMPD_237_ribflv_c: R_LMPD_237_ribflv_c -1000 BU_LMPD_237_ribflv_c < 0
G_LMPD_238_ribflv_c: -0.59209 LC_3pg_c +6.513 LC_adp_c +0.59209 LC_akg_c +1.1842 LC_amp_c -0.59209 LC_asp-L_c -7.6972 LC_atp_c -0.59209 LC_co2_c -1.1842 LC_db4p_c +0.59209 LC_fum_c -0.59209 LC_glu-L_c -1.1842 LC_nad_c +1.1842 LC_nadh_c -1.1842 LC_nh4_c +7.6972 LC_pi_c +1.7763 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_ribflv_c - DG_LMPD_238_ribflv_c + DGo_LMPD_238_ribflv_c = 0
FU_LMPD_238_ribflv_c: DG_LMPD_238_ribflv_c +1000 FU_LMPD_238_ribflv_c < 1000
BU_LMPD_238_ribflv_c: - DG_LMPD_238_ribflv_c +1000 BU_LMPD_238_ribflv_c < 1000
SU_LMPD_238_ribflv_c: FU_LMPD_238_ribflv_c + BU_LMPD_238_ribflv_c < 1
UF_LMPD_238_ribflv_c: F_LMPD_238_ribflv_c -1000 FU_LMPD_238_ribflv_c < 0
UR_LMPD_238_ribflv_c: R_LMPD_238_ribflv_c -1000 BU_LMPD_238_ribflv_c < 0
G_LMPD_239_ser-L_c: -0.59209 LC_3pg_c +0.59209 LC_akg_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_pi_c +0.59209 LC_ser-L_c - DG_LMPD_239_ser-L_c + DGo_LMPD_239_ser-L_c = 0
FU_LMPD_239_ser-L_c: DG_LMPD_239_ser-L_c +1000 FU_LMPD_239_ser-L_c < 1000
BU_LMPD_239_ser-L_c: - DG_LMPD_239_ser-L_c +1000 BU_LMPD_239_ser-L_c < 1000
SU_LMPD_239_ser-L_c: FU_LMPD_239_ser-L_c + BU_LMPD_239_ser-L_c < 1
UF_LMPD_239_ser-L_c: F_LMPD_239_ser-L_c -1000 FU_LMPD_239_ser-L_c < 0
UR_LMPD_239_ser-L_c: R_LMPD_239_ser-L_c -1000 BU_LMPD_239_ser-L_c < 0
SU_LMPD_240_sheme_c: FU_LMPD_240_sheme_c + BU_LMPD_240_sheme_c < 1
UF_LMPD_240_sheme_c: F_LMPD_240_sheme_c -1000 FU_LMPD_240_sheme_c < 0
UR_LMPD_240_sheme_c: R_LMPD_240_sheme_c -1000 BU_LMPD_240_sheme_c < 0
G_LMPD_241_spmd_c: -1.1842 LC_3pg_c +1.1842 LC_ac_c -1.1842 LC_accoa_c +1.1842 LC_adp_c +1.7763 LC_akg_c +2.3684 LC_amp_c -0.59209 LC_asp-L_c -3.5525 LC_atp_c +1.1842 LC_co2_c +1.7763 LC_coa_c -2.3684 LC_glu-L_c +0.59209 LC_gly_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +3.5525 LC_nadp_c -3.5525 LC_nadph_c +0.59209 LC_nh4_c +4.1446 LC_pi_c +1.7763 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_spmd_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c +0.59209 LC_5mtr_p - DG_LMPD_241_spmd_c + DGo_LMPD_241_spmd_c = 0
FU_LMPD_241_spmd_c: DG_LMPD_241_spmd_c +1000 FU_LMPD_241_spmd_c < 1000
BU_LMPD_241_spmd_c: - DG_LMPD_241_spmd_c +1000 BU_LMPD_241_spmd_c < 1000
SU_LMPD_241_spmd_c: FU_LMPD_241_spmd_c + BU_LMPD_241_spmd_c < 1
UF_LMPD_241_spmd_c: F_LMPD_241_spmd_c -1000 FU_LMPD_241_spmd_c < 0
UR_LMPD_241_spmd_c: R_LMPD_241_spmd_c -1000 BU_LMPD_241_spmd_c < 0
G_LMPD_242_spmd_c: -1.1842 LC_3pg_c +1.1842 LC_ac_c -1.1842 LC_accoa_c +2.3684 LC_adp_c +1.7763 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -4.1446 LC_atp_c +1.1842 LC_co2_c +1.7763 LC_coa_c -2.3684 LC_glu-L_c +0.59209 LC_gly_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +3.5525 LC_nadp_c -3.5525 LC_nadph_c +0.59209 LC_nh4_c +4.1446 LC_pi_c +1.7763 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_spmd_c +0.59209 LC_succ_c -0.59209 LC_succoa_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c +0.59209 LC_5mtr_p - DG_LMPD_242_spmd_c + DGo_LMPD_242_spmd_c = 0
FU_LMPD_242_spmd_c: DG_LMPD_242_spmd_c +1000 FU_LMPD_242_spmd_c < 1000
BU_LMPD_242_spmd_c: - DG_LMPD_242_spmd_c +1000 BU_LMPD_242_spmd_c < 1000
SU_LMPD_242_spmd_c: FU_LMPD_242_spmd_c + BU_LMPD_242_spmd_c < 1
UF_LMPD_242_spmd_c: F_LMPD_242_spmd_c -1000 FU_LMPD_242_spmd_c < 0
UR_LMPD_242_spmd_c: R_LMPD_242_spmd_c -1000 BU_LMPD_242_spmd_c < 0
G_LMPD_243_thf_c: -1.1842 LC_3pg_c +6.513 LC_adp_c +1.1842 LC_akg_c +2.3684 LC_amp_c -0.59209 LC_asp-L_c -8.8813 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_for_c +0.59209 LC_fum_c -1.7763 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_glyclt_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.0655 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c +0.59209 LC_thf_c - DG_LMPD_243_thf_c + DGo_LMPD_243_thf_c = 0
FU_LMPD_243_thf_c: DG_LMPD_243_thf_c +1000 FU_LMPD_243_thf_c < 1000
BU_LMPD_243_thf_c: - DG_LMPD_243_thf_c +1000 BU_LMPD_243_thf_c < 1000
SU_LMPD_243_thf_c: FU_LMPD_243_thf_c + BU_LMPD_243_thf_c < 1
UF_LMPD_243_thf_c: F_LMPD_243_thf_c -1000 FU_LMPD_243_thf_c < 0
UR_LMPD_243_thf_c: R_LMPD_243_thf_c -1000 BU_LMPD_243_thf_c < 0
G_LMPD_244_thf_c: -1.1842 LC_3pg_c +7.6972 LC_adp_c +1.1842 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -9.4734 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_for_c +0.59209 LC_fum_c -1.7763 LC_glu-L_c +0.59209 LC_gly_c +0.59209 LC_glyclt_c -2.3684 LC_nad_c +2.3684 LC_nadh_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.0655 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c +0.59209 LC_thf_c - DG_LMPD_244_thf_c + DGo_LMPD_244_thf_c = 0
FU_LMPD_244_thf_c: DG_LMPD_244_thf_c +1000 FU_LMPD_244_thf_c < 1000
BU_LMPD_244_thf_c: - DG_LMPD_244_thf_c +1000 BU_LMPD_244_thf_c < 1000
SU_LMPD_244_thf_c: FU_LMPD_244_thf_c + BU_LMPD_244_thf_c < 1
UF_LMPD_244_thf_c: F_LMPD_244_thf_c -1000 FU_LMPD_244_thf_c < 0
UR_LMPD_244_thf_c: R_LMPD_244_thf_c -1000 BU_LMPD_244_thf_c < 0
G_LMPD_245_thf_c: -0.59209 LC_3pg_c +7.1051 LC_adp_c +0.59209 LC_akg_c +2.3684 LC_amp_c -0.59209 LC_asp-L_c -9.4734 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_glyclt_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.0655 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c +0.59209 LC_thf_c - DG_LMPD_245_thf_c + DGo_LMPD_245_thf_c = 0
FU_LMPD_245_thf_c: DG_LMPD_245_thf_c +1000 FU_LMPD_245_thf_c < 1000
BU_LMPD_245_thf_c: - DG_LMPD_245_thf_c +1000 BU_LMPD_245_thf_c < 1000
SU_LMPD_245_thf_c: FU_LMPD_245_thf_c + BU_LMPD_245_thf_c < 1
UF_LMPD_245_thf_c: F_LMPD_245_thf_c -1000 FU_LMPD_245_thf_c < 0
UR_LMPD_245_thf_c: R_LMPD_245_thf_c -1000 BU_LMPD_245_thf_c < 0
G_LMPD_246_thf_c: -0.59209 LC_3pg_c +8.2893 LC_adp_c +0.59209 LC_akg_c +1.7763 LC_amp_c -0.59209 LC_asp-L_c -10.0655 LC_atp_c -0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_fum_c -1.1842 LC_glu-L_c +0.59209 LC_glyclt_c -1.7763 LC_nad_c +1.7763 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -2.3684 LC_nh4_c -1.1842 LC_pep_c +10.0655 LC_pi_c +2.3684 LC_ppi_c +0.59209 LC_pyr_c -0.59209 LC_r5p_c +0.59209 LC_thf_c - DG_LMPD_246_thf_c + DGo_LMPD_246_thf_c = 0
FU_LMPD_246_thf_c: DG_LMPD_246_thf_c +1000 FU_LMPD_246_thf_c < 1000
BU_LMPD_246_thf_c: - DG_LMPD_246_thf_c +1000 BU_LMPD_246_thf_c < 1000
SU_LMPD_246_thf_c: FU_LMPD_246_thf_c + BU_LMPD_246_thf_c < 1
UF_LMPD_246_thf_c: F_LMPD_246_thf_c -1000 FU_LMPD_246_thf_c < 0
UR_LMPD_246_thf_c: R_LMPD_246_thf_c -1000 BU_LMPD_246_thf_c < 0
G_LMPD_247_thmpp_c: -1.1842 LC_3pg_c +0.59209 LC_4crsol_c +0.59209 LC_5drib_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +4.7367 LC_adp_c +1.7763 LC_akg_c +0.59209 LC_ala-L_c +3.5525 LC_amp_c -8.2893 LC_atp_c +1.7763 LC_co2_c +0.59209 LC_coa_c -0.59209 LC_e4p_c +1.1842 LC_for_c -0.59209 LC_g3p_c -1.7763 LC_glu-L_c -2.3684 LC_nad_c +2.3684 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -1.1842 LC_nh4_c -1.1842 LC_pep_c +8.2893 LC_pi_c +3.5525 LC_ppi_c -0.59209 LC_pyr_c -1.1842 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_thmpp_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_247_thmpp_c + DGo_LMPD_247_thmpp_c = 0
FU_LMPD_247_thmpp_c: DG_LMPD_247_thmpp_c +1000 FU_LMPD_247_thmpp_c < 1000
BU_LMPD_247_thmpp_c: - DG_LMPD_247_thmpp_c +1000 BU_LMPD_247_thmpp_c < 1000
SU_LMPD_247_thmpp_c: FU_LMPD_247_thmpp_c + BU_LMPD_247_thmpp_c < 1
UF_LMPD_247_thmpp_c: F_LMPD_247_thmpp_c -1000 FU_LMPD_247_thmpp_c < 0
UR_LMPD_247_thmpp_c: R_LMPD_247_thmpp_c -1000 BU_LMPD_247_thmpp_c < 0
G_LMPD_248_thmpp_c: -1.1842 LC_3pg_c +0.59209 LC_4crsol_c +0.59209 LC_5drib_c +0.59209 LC_ac_c -0.59209 LC_accoa_c +5.9209 LC_adp_c +1.7763 LC_akg_c +0.59209 LC_ala-L_c +2.9604 LC_amp_c -8.8813 LC_atp_c +1.7763 LC_co2_c +0.59209 LC_coa_c -0.59209 LC_e4p_c +1.1842 LC_for_c -0.59209 LC_g3p_c -1.7763 LC_glu-L_c -2.3684 LC_nad_c +2.3684 LC_nadh_c +2.9604 LC_nadp_c -2.9604 LC_nadph_c -1.1842 LC_nh4_c -1.1842 LC_pep_c +8.2893 LC_pi_c +3.5525 LC_ppi_c -0.59209 LC_pyr_c -1.1842 LC_r5p_c -0.59209 LC_so4_c +0.59209 LC_thmpp_c +0.59209 LC_trdox_c -0.59209 LC_trdrd_c - DG_LMPD_248_thmpp_c + DGo_LMPD_248_thmpp_c = 0
FU_LMPD_248_thmpp_c: DG_LMPD_248_thmpp_c +1000 FU_LMPD_248_thmpp_c < 1000
BU_LMPD_248_thmpp_c: - DG_LMPD_248_thmpp_c +1000 BU_LMPD_248_thmpp_c < 1000
SU_LMPD_248_thmpp_c: FU_LMPD_248_thmpp_c + BU_LMPD_248_thmpp_c < 1
UF_LMPD_248_thmpp_c: F_LMPD_248_thmpp_c -1000 FU_LMPD_248_thmpp_c < 0
UR_LMPD_248_thmpp_c: R_LMPD_248_thmpp_c -1000 BU_LMPD_248_thmpp_c < 0
G_LMPD_249_thr-L_c: +1.1842 LC_adp_c -0.59209 LC_asp-L_c -1.1842 LC_atp_c +1.1842 LC_nadp_c -1.1842 LC_nadph_c +1.1842 LC_pi_c +0.59209 LC_thr-L_c - DG_LMPD_249_thr-L_c + DGo_LMPD_249_thr-L_c = 0
FU_LMPD_249_thr-L_c: DG_LMPD_249_thr-L_c +1000 FU_LMPD_249_thr-L_c < 1000
BU_LMPD_249_thr-L_c: - DG_LMPD_249_thr-L_c +1000 BU_LMPD_249_thr-L_c < 1000
SU_LMPD_249_thr-L_c: FU_LMPD_249_thr-L_c + BU_LMPD_249_thr-L_c < 1
UF_LMPD_249_thr-L_c: F_LMPD_249_thr-L_c -1000 FU_LMPD_249_thr-L_c < 0
UR_LMPD_249_thr-L_c: R_LMPD_249_thr-L_c -1000 BU_LMPD_249_thr-L_c < 0
G_LMPD_250_trp-L_c: +1.1842 LC_adp_c +0.59209 LC_amp_c -1.7763 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_e4p_c +0.59209 LC_g3p_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_nh4_c -1.1842 LC_pep_c +2.9604 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_trp-L_c - DG_LMPD_250_trp-L_c + DGo_LMPD_250_trp-L_c = 0
FU_LMPD_250_trp-L_c: DG_LMPD_250_trp-L_c +1000 FU_LMPD_250_trp-L_c < 1000
BU_LMPD_250_trp-L_c: - DG_LMPD_250_trp-L_c +1000 BU_LMPD_250_trp-L_c < 1000
SU_LMPD_250_trp-L_c: FU_LMPD_250_trp-L_c + BU_LMPD_250_trp-L_c < 1
UF_LMPD_250_trp-L_c: F_LMPD_250_trp-L_c -1000 FU_LMPD_250_trp-L_c < 0
UR_LMPD_250_trp-L_c: R_LMPD_250_trp-L_c -1000 BU_LMPD_250_trp-L_c < 0
G_LMPD_251_tyr-L_c: +0.59209 LC_adp_c +0.59209 LC_akg_c -0.59209 LC_atp_c +0.59209 LC_co2_c -0.59209 LC_e4p_c -0.59209 LC_glu-L_c -0.59209 LC_nad_c +0.59209 LC_nadh_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_pep_c +2.3684 LC_pi_c +0.59209 LC_tyr-L_c - DG_LMPD_251_tyr-L_c + DGo_LMPD_251_tyr-L_c = 0
FU_LMPD_251_tyr-L_c: DG_LMPD_251_tyr-L_c +1000 FU_LMPD_251_tyr-L_c < 1000
BU_LMPD_251_tyr-L_c: - DG_LMPD_251_tyr-L_c +1000 BU_LMPD_251_tyr-L_c < 1000
SU_LMPD_251_tyr-L_c: FU_LMPD_251_tyr-L_c + BU_LMPD_251_tyr-L_c < 1
UF_LMPD_251_tyr-L_c: F_LMPD_251_tyr-L_c -1000 FU_LMPD_251_tyr-L_c < 0
UR_LMPD_251_tyr-L_c: R_LMPD_251_tyr-L_c -1000 BU_LMPD_251_tyr-L_c < 0
SU_LMPD_252_udcpdp_c: FU_LMPD_252_udcpdp_c + BU_LMPD_252_udcpdp_c < 1
UF_LMPD_252_udcpdp_c: F_LMPD_252_udcpdp_c -1000 FU_LMPD_252_udcpdp_c < 0
UR_LMPD_252_udcpdp_c: R_LMPD_252_udcpdp_c -1000 BU_LMPD_252_udcpdp_c < 0
G_LMPD_253_utp_c: +1.7763 LC_adp_c +0.59209 LC_amp_c -0.59209 LC_asp-L_c -2.3684 LC_atp_c +0.59209 LC_dhor-S_c -0.59209 LC_nh4_c -0.59209 LC_orot_c +0.59209 LC_pi_c +0.59209 LC_ppi_c -0.59209 LC_r5p_c +0.59209 LC_utp_c - DG_LMPD_253_utp_c + DGo_LMPD_253_utp_c = 0
FU_LMPD_253_utp_c: DG_LMPD_253_utp_c +1000 FU_LMPD_253_utp_c < 1000
BU_LMPD_253_utp_c: - DG_LMPD_253_utp_c +1000 BU_LMPD_253_utp_c < 1000
SU_LMPD_253_utp_c: FU_LMPD_253_utp_c + BU_LMPD_253_utp_c < 1
UF_LMPD_253_utp_c: F_LMPD_253_utp_c -1000 FU_LMPD_253_utp_c < 0
UR_LMPD_253_utp_c: R_LMPD_253_utp_c -1000 BU_LMPD_253_utp_c < 0
G_LMPD_254_val-L_c: +0.59209 LC_akg_c +0.59209 LC_co2_c -0.59209 LC_glu-L_c +0.59209 LC_nadp_c -0.59209 LC_nadph_c -1.1842 LC_pyr_c +0.59209 LC_val-L_c - DG_LMPD_254_val-L_c + DGo_LMPD_254_val-L_c = 0
FU_LMPD_254_val-L_c: DG_LMPD_254_val-L_c +1000 FU_LMPD_254_val-L_c < 1000
BU_LMPD_254_val-L_c: - DG_LMPD_254_val-L_c +1000 BU_LMPD_254_val-L_c < 1000
SU_LMPD_254_val-L_c: FU_LMPD_254_val-L_c + BU_LMPD_254_val-L_c < 1
UF_LMPD_254_val-L_c: F_LMPD_254_val-L_c -1000 FU_LMPD_254_val-L_c < 0
UR_LMPD_254_val-L_c: R_LMPD_254_val-L_c -1000 BU_LMPD_254_val-L_c < 0
Bounds
 0  <= F_DM_4CRSOL  <= 1000
 0  <= R_DM_4CRSOL  <= 0
 0  <= F_DM_5DRIB  <= 1000
 0  <= R_DM_5DRIB  <= 0
 0  <= F_DM_AMOB  <= 1000
 0  <= R_DM_AMOB  <= 0
 0  <= F_DM_MTHTHF  <= 1000
 0  <= R_DM_MTHTHF  <= 0
 0  <= F_Ec_biomass_iJO1366_WT_53p95M  <= 1000
 0  <= R_Ec_biomass_iJO1366_WT_53p95M  <= 0
 0  <= F_DM_5mtr_e  <= 1000
 0  <= R_DM_5mtr_e  <= 0
 0  <= F_DM_ac_e  <= 1000
 0  <= R_DM_ac_e  <= 0
 0  <= F_DM_acald_e  <= 1000
 0  <= R_DM_acald_e  <= 0
 0  <= F_DM_acser_e  <= 1000
 0  <= R_DM_acser_e  <= 0
 0  <= F_DM_akg_e  <= 1000
 0  <= R_DM_akg_e  <= 0
 0  <= F_DM_ala-D_e  <= 1000
 0  <= R_DM_ala-D_e  <= 0
 0  <= F_DM_ala-L_e  <= 1000
 0  <= R_DM_ala-L_e  <= 0
 0  <= F_DM_alaala_e  <= 1000
 0  <= R_DM_alaala_e  <= 0
 0  <= F_DM_arg-L_e  <= 1000
 0  <= R_DM_arg-L_e  <= 0
 0  <= F_DM_asn-L_e  <= 1000
 0  <= R_DM_asn-L_e  <= 0
 0  <= F_DM_asp-L_e  <= 1000
 0  <= R_DM_asp-L_e  <= 0
 0  <= F_DM_ca2_e  <= 1000
 0  <= R_DM_ca2_e  <= 1000
 0  <= F_DM_cbl1_e  <= 1000
 0  <= R_DM_cbl1_e  <= 1.000000e-02
 0  <= F_DM_cit_e  <= 1000
 0  <= R_DM_cit_e  <= 0
 0  <= F_DM_cl_e  <= 1000
 0  <= R_DM_cl_e  <= 1000
 0  <= F_DM_co2_e  <= 1000
 0  <= R_DM_co2_e  <= 1000
 0  <= F_DM_cobalt2_e  <= 1000
 0  <= R_DM_cobalt2_e  <= 1000
 0  <= F_DM_colipa_e  <= 1000
 0  <= R_DM_colipa_e  <= 0
 0  <= F_DM_cu2_e  <= 1000
 0  <= R_DM_cu2_e  <= 1000
 0  <= F_DM_cys-L_e  <= 1000
 0  <= R_DM_cys-L_e  <= 0
 0  <= F_DM_dha_e  <= 1000
 0  <= R_DM_dha_e  <= 0
 0  <= F_DM_enter_e  <= 1000
 0  <= R_DM_enter_e  <= 0
 0  <= F_DM_etoh_e  <= 1000
 0  <= R_DM_etoh_e  <= 0
 0  <= F_DM_fe2_e  <= 1000
 0  <= R_DM_fe2_e  <= 1000
 0  <= F_DM_fe3_e  <= 1000
 0  <= R_DM_fe3_e  <= 1000
 0  <= F_DM_for_e  <= 1000
 0  <= R_DM_for_e  <= 0
 0  <= F_DM_glc_e  <= 0
 0  <= R_DM_glc_e  <= 8.160000e+00
 0  <= F_DM_glcn_e  <= 1000
 0  <= R_DM_glcn_e  <= 0
 0  <= F_DM_glu-L_e  <= 1000
 0  <= R_DM_glu-L_e  <= 0
 0  <= F_DM_gly_e  <= 1000
 0  <= R_DM_gly_e  <= 0
 0  <= F_DM_glyc_e  <= 1000
 0  <= R_DM_glyc_e  <= 0
 0  <= F_DM_glyc-R_e  <= 1000
 0  <= R_DM_glyc-R_e  <= 0
 0  <= F_DM_glycogenn1_c  <= 0
 0  <= R_DM_glycogenn1_c  <= 10
 0  <= F_DM_glyclt_e  <= 1000
 0  <= R_DM_glyclt_e  <= 0
 0  <= F_DM_gthrd_e  <= 1000
 0  <= R_DM_gthrd_e  <= 0
 0  <= F_DM_h_e  <= 1000
 0  <= R_DM_h_e  <= 1000
 0  <= F_DM_h2o_e  <= 1000
 0  <= R_DM_h2o_e  <= 1000
 0  <= F_DM_his-L_e  <= 1000
 0  <= R_DM_his-L_e  <= 0
 0  <= F_DM_ile-L_e  <= 1000
 0  <= R_DM_ile-L_e  <= 0
 0  <= F_DM_k_e  <= 1000
 0  <= R_DM_k_e  <= 1000
 0  <= F_DM_lac-D_e  <= 1000
 0  <= R_DM_lac-D_e  <= 0
 0  <= F_DM_leu-L_e  <= 1000
 0  <= R_DM_leu-L_e  <= 0
 0  <= F_DM_lys-L_e  <= 1000
 0  <= R_DM_lys-L_e  <= 0
 0  <= F_DM_mal-L_e  <= 1000
 0  <= R_DM_mal-L_e  <= 0
 0  <= F_DM_meoh_e  <= 1000
 0  <= R_DM_meoh_e  <= 0
 0  <= F_DM_mg2_e  <= 1000
 0  <= R_DM_mg2_e  <= 1000
 0  <= F_DM_mn2_e  <= 1000
 0  <= R_DM_mn2_e  <= 1000
 0  <= F_DM_mobd_e  <= 1000
 0  <= R_DM_mobd_e  <= 1000
 0  <= F_DM_nh4_e  <= 1000
 0  <= R_DM_nh4_e  <= 1000
 0  <= F_DM_ni2_e  <= 1000
 0  <= R_DM_ni2_e  <= 1000
 0  <= F_DM_o2_e  <= 1000
 0  <= R_DM_o2_e  <= 1000
 0  <= F_DM_phe-L_e  <= 1000
 0  <= R_DM_phe-L_e  <= 0
 0  <= F_DM_pheme_e  <= 1000
 0  <= R_DM_pheme_e  <= 0
 0  <= F_DM_pi_e  <= 1000
 0  <= R_DM_pi_e  <= 1000
 0  <= F_DM_pro-L_e  <= 1000
 0  <= R_DM_pro-L_e  <= 0
 0  <= F_DM_pyr_e  <= 1000
 0  <= R_DM_pyr_e  <= 0
 0  <= F_DM_ser-L_e  <= 1000
 0  <= R_DM_ser-L_e  <= 0
 0  <= F_DM_so4_e  <= 1000
 0  <= R_DM_so4_e  <= 1000
 0  <= F_DM_spmd_e  <= 1000
 0  <= R_DM_spmd_e  <= 0
 0  <= F_DM_succ_e  <= 1000
 0  <= R_DM_succ_e  <= 0
 0  <= F_DM_thr-L_e  <= 1000
 0  <= R_DM_thr-L_e  <= 0
 0  <= F_DM_trp-L_e  <= 1000
 0  <= R_DM_trp-L_e  <= 0
 0  <= F_DM_tyr-L_e  <= 1000
 0  <= R_DM_tyr-L_e  <= 0
 0  <= F_DM_val-L_e  <= 1000
 0  <= R_DM_val-L_e  <= 0
 0  <= F_DM_zn2_e  <= 1000
 0  <= R_DM_zn2_e  <= 1000
 0  <= F_5MTRtex  <= 1000
 0  <= R_5MTRtex  <= 1000
 0  <= F_ACALD  <= 1000
 0  <= R_ACALD  <= 1000
 0  <= F_ACALDtex  <= 1000
 0  <= R_ACALDtex  <= 1000
 0  <= F_ACALDtpp  <= 1000
 0  <= R_ACALDtpp  <= 1000
 0  <= F_ACKr  <= 1000
 0  <= R_ACKr  <= 1000
 0  <= F_ACONTa  <= 1000
 0  <= R_ACONTa  <= 1000
 0  <= F_ACONTb  <= 1000
 0  <= R_ACONTb  <= 1000
 0  <= F_ACS  <= 1000
 0  <= R_ACS  <= 0
 0  <= F_ACSERtex  <= 1000
 0  <= R_ACSERtex  <= 1000
 0  <= F_ACt2rpp  <= 1000
 0  <= R_ACt2rpp  <= 1000
 0  <= F_ACtex  <= 1000
 0  <= R_ACtex  <= 1000
 0  <= F_ADK1  <= 1000
 0  <= R_ADK1  <= 1000
 0  <= F_AKGDH  <= 1000
 0  <= R_AKGDH  <= 0
 0  <= F_AKGt2rpp  <= 1000
 0  <= R_AKGt2rpp  <= 1000
 0  <= F_AKGtex  <= 1000
 0  <= R_AKGtex  <= 1000
 0  <= F_ALAALAtex  <= 1000
 0  <= R_ALAALAtex  <= 1000
 0  <= F_ALAt2rpp  <= 1000
 0  <= R_ALAt2rpp  <= 1000
 0  <= F_ALAtex  <= 1000
 0  <= R_ALAtex  <= 1000
 0  <= F_ALCD2x  <= 1000
 0  <= R_ALCD2x  <= 1000
 0  <= F_ALDD2x  <= 1000
 0  <= R_ALDD2x  <= 0
 0  <= F_ALDD2y  <= 1000
 0  <= R_ALDD2y  <= 0
 0  <= F_ARGt3pp  <= 1000
 0  <= R_ARGt3pp  <= 0
 0  <= F_ARGtex  <= 1000
 0  <= R_ARGtex  <= 1000
 0  <= F_ASNt2rpp  <= 1000
 0  <= R_ASNt2rpp  <= 1000
 0  <= F_ASNtex  <= 1000
 0  <= R_ASNtex  <= 1000
 0  <= F_ASPO3  <= 1000
 0  <= R_ASPO3  <= 0
 0  <= F_ASPO4  <= 1000
 0  <= R_ASPO4  <= 0
 0  <= F_ASPO5  <= 1000
 0  <= R_ASPO5  <= 0
 0  <= F_ASPO6  <= 1000
 0  <= R_ASPO6  <= 0
 0  <= F_ASPT  <= 1000
 0  <= R_ASPT  <= 0
 0  <= F_ASPTA  <= 1000
 0  <= R_ASPTA  <= 1000
 0  <= F_ASPt2rpp  <= 1000
 0  <= R_ASPt2rpp  <= 1000
 0  <= F_ASPtex  <= 1000
 0  <= R_ASPtex  <= 1000
 0  <= F_ATPM  <= 1000
 0  <= R_ATPM  <= 0
 0  <= F_ATPS4rpp  <= 1000
 0  <= R_ATPS4rpp  <= 1000
 0  <= F_CA2tex  <= 1000
 0  <= R_CA2tex  <= 1000
 0  <= F_CAt6pp  <= 1000
 0  <= R_CAt6pp  <= 1000
 0  <= F_CBL1tonex  <= 1000
 0  <= R_CBL1tonex  <= 0
 0  <= F_CITL  <= 1000
 0  <= R_CITL  <= 0
 0  <= F_CITt3pp  <= 1000
 0  <= R_CITt3pp  <= 0
 0  <= F_CITt7pp  <= 1000
 0  <= R_CITt7pp  <= 0
 0  <= F_CITtex  <= 1000
 0  <= R_CITtex  <= 1000
 0  <= F_CLt3_2pp  <= 1000
 0  <= R_CLt3_2pp  <= 0
 0  <= F_CLtex  <= 1000
 0  <= R_CLtex  <= 1000
 0  <= F_CO2tex  <= 1000
 0  <= R_CO2tex  <= 1000
 0  <= F_CO2tpp  <= 1000
 0  <= R_CO2tpp  <= 1000
 0  <= F_COBALT2tex  <= 1000
 0  <= R_COBALT2tex  <= 1000
 0  <= F_COBALT2tpp  <= 1000
 0  <= R_COBALT2tpp  <= 0
 0  <= F_COLIPAabcpp  <= 1000
 0  <= R_COLIPAabcpp  <= 0
 0  <= F_COLIPAabctex  <= 1000
 0  <= R_COLIPAabctex  <= 0
 0  <= F_CS  <= 1000
 0  <= R_CS  <= 0
 0  <= F_CU2tex  <= 1000
 0  <= R_CU2tex  <= 1000
 0  <= F_CU2tpp  <= 1000
 0  <= R_CU2tpp  <= 0
 0  <= F_CYSabc2pp  <= 1000
 0  <= R_CYSabc2pp  <= 0
 0  <= F_CYStex  <= 1000
 0  <= R_CYStex  <= 1000
 0  <= F_CYStpp  <= 1000
 0  <= R_CYStpp  <= 0
 0  <= F_CYTBD2pp  <= 1000
 0  <= R_CYTBD2pp  <= 0
 0  <= F_CYTBDpp  <= 1000
 0  <= R_CYTBDpp  <= 0
 0  <= F_CYTBO3_4pp  <= 1000
 0  <= R_CYTBO3_4pp  <= 0
 0  <= F_D-LACt2pp  <= 1000
 0  <= R_D-LACt2pp  <= 1000
 0  <= F_D-LACtex  <= 1000
 0  <= R_D-LACtex  <= 1000
 0  <= F_DALAtex  <= 1000
 0  <= R_DALAtex  <= 1000
 0  <= F_DB4PS  <= 1000
 0  <= R_DB4PS  <= 0
 0  <= F_DHAPT  <= 1000
 0  <= R_DHAPT  <= 0
 0  <= F_DHAtex  <= 1000
 0  <= R_DHAtex  <= 1000
 0  <= F_DHAtpp  <= 1000
 0  <= R_DHAtpp  <= 1000
 0  <= F_DHORD2  <= 1000
 0  <= R_DHORD2  <= 0
 0  <= F_DHORD5  <= 1000
 0  <= R_DHORD5  <= 0
 0  <= F_DHORDfum  <= 1000
 0  <= R_DHORDfum  <= 0
 0  <= F_EDA  <= 1000
 0  <= R_EDA  <= 0
 0  <= F_EDD  <= 1000
 0  <= R_EDD  <= 0
 0  <= F_ENO  <= 1000
 0  <= R_ENO  <= 1000
 0  <= F_ETOHtex  <= 1000
 0  <= R_ETOHtex  <= 1000
 0  <= F_ETOHtrpp  <= 1000
 0  <= R_ETOHtrpp  <= 1000
 0  <= F_F6PA  <= 1000
 0  <= R_F6PA  <= 1000
 0  <= F_FBA  <= 1000
 0  <= R_FBA  <= 1000
 0  <= F_FBA3  <= 1000
 0  <= R_FBA3  <= 1000
 0  <= F_FBP  <= 1000
 0  <= R_FBP  <= 0
 0  <= F_FDH4pp  <= 1000
 0  <= R_FDH4pp  <= 0
 0  <= F_FDH5pp  <= 1000
 0  <= R_FDH5pp  <= 0
 0  <= F_FE2abcpp  <= 1000
 0  <= R_FE2abcpp  <= 0
 0  <= F_FE2t2pp  <= 1000
 0  <= R_FE2t2pp  <= 0
 0  <= F_FE2tex  <= 1000
 0  <= R_FE2tex  <= 1000
 0  <= F_FE2tpp  <= 1000
 0  <= R_FE2tpp  <= 0
 0  <= F_FE3abcpp  <= 1000
 0  <= R_FE3abcpp  <= 0
 0  <= F_FE3tex  <= 1000
 0  <= R_FE3tex  <= 1000
 0  <= F_FEENTERtex  <= 1000
 0  <= R_FEENTERtex  <= 0
 0  <= F_FEENTERtpp  <= 1000
 0  <= R_FEENTERtpp  <= 0
 0  <= F_FLDR2  <= 1000
 0  <= R_FLDR2  <= 0
 0  <= F_FORtex  <= 1000
 0  <= R_FORtex  <= 1000
 0  <= F_FORtppi  <= 1000
 0  <= R_FORtppi  <= 0
 0  <= F_FRD3  <= 1000
 0  <= R_FRD3  <= 0
 0  <= F_FUM  <= 1000
 0  <= R_FUM  <= 1000
 0  <= F_G3PD2  <= 1000
 0  <= R_G3PD2  <= 1000
 0  <= F_G3PD5  <= 1000
 0  <= R_G3PD5  <= 0
 0  <= F_G3PD6  <= 1000
 0  <= R_G3PD6  <= 0
 0  <= F_G3PD7  <= 1000
 0  <= R_G3PD7  <= 0
 0  <= F_G3PT  <= 1000
 0  <= R_G3PT  <= 0
 0  <= F_G6PDH2r  <= 1000
 0  <= R_G6PDH2r  <= 1000
 0  <= F_G6PP  <= 1000
 0  <= R_G6PP  <= 0
 0  <= F_GAPD  <= 1000
 0  <= R_GAPD  <= 1000
 0  <= F_GLBRAN2  <= 1000
 0  <= R_GLBRAN2  <= 0
 0  <= F_GLCDpp  <= 1000
 0  <= R_GLCDpp  <= 0
 0  <= F_GLCNtex  <= 1000
 0  <= R_GLCNtex  <= 1000
 0  <= F_GLCP  <= 1000
 0  <= R_GLCP  <= 0
 0  <= F_GLCS1  <= 1000
 0  <= R_GLCS1  <= 0
 0  <= F_GLCabcpp  <= 1000
 0  <= R_GLCabcpp  <= 0
 0  <= F_GLCptspp  <= 1000
 0  <= R_GLCptspp  <= 0
 0  <= F_GLCt2pp  <= 1000
 0  <= R_GLCt2pp  <= 0
 0  <= F_GLCtex  <= 1000
 0  <= R_GLCtex  <= 1000
 0  <= F_GLCtexi  <= 1000
 0  <= R_GLCtexi  <= 0
 0  <= F_GLDBRAN2  <= 1000
 0  <= R_GLDBRAN2  <= 0
 0  <= F_GLGC  <= 1000
 0  <= R_GLGC  <= 0
 0  <= F_GLUDy  <= 1000
 0  <= R_GLUDy  <= 1000
 0  <= F_GLUabcpp  <= 1000
 0  <= R_GLUabcpp  <= 0
 0  <= F_GLUt2rpp  <= 1000
 0  <= R_GLUt2rpp  <= 1000
 0  <= F_GLUt4pp  <= 1000
 0  <= R_GLUt4pp  <= 0
 0  <= F_GLUtex  <= 1000
 0  <= R_GLUtex  <= 1000
 0  <= F_GLXCL  <= 1000
 0  <= R_GLXCL  <= 0
 0  <= F_GLYCAt2rpp  <= 1000
 0  <= R_GLYCAt2rpp  <= 1000
 0  <= F_GLYCAtex  <= 1000
 0  <= R_GLYCAtex  <= 1000
 0  <= F_GLYCDx  <= 1000
 0  <= R_GLYCDx  <= 0
 0  <= F_GLYCK  <= 1000
 0  <= R_GLYCK  <= 0
 0  <= F_GLYCK2  <= 1000
 0  <= R_GLYCK2  <= 0
 0  <= F_GLYCLTDx  <= 1000
 0  <= R_GLYCLTDx  <= 0
 0  <= F_GLYCLTDy  <= 1000
 0  <= R_GLYCLTDy  <= 0
 0  <= F_GLYCLTt2rpp  <= 1000
 0  <= R_GLYCLTt2rpp  <= 1000
 0  <= F_GLYCLTtex  <= 1000
 0  <= R_GLYCLTtex  <= 1000
 0  <= F_GLYCTO2  <= 1000
 0  <= R_GLYCTO2  <= 0
 0  <= F_GLYCTO3  <= 1000
 0  <= R_GLYCTO3  <= 0
 0  <= F_GLYCTO4  <= 1000
 0  <= R_GLYCTO4  <= 0
 0  <= F_GLYCtex  <= 1000
 0  <= R_GLYCtex  <= 1000
 0  <= F_GLYCtpp  <= 1000
 0  <= R_GLYCtpp  <= 1000
 0  <= F_GLYK  <= 1000
 0  <= R_GLYK  <= 0
 0  <= F_GLYt2rpp  <= 1000
 0  <= R_GLYt2rpp  <= 1000
 0  <= F_GLYtex  <= 1000
 0  <= R_GLYtex  <= 1000
 0  <= F_GND  <= 1000
 0  <= R_GND  <= 0
 0  <= F_GTHRDabc2pp  <= 1000
 0  <= R_GTHRDabc2pp  <= 0
 0  <= F_GTHRDtex  <= 1000
 0  <= R_GTHRDtex  <= 1000
 0  <= F_H2Otex  <= 1000
 0  <= R_H2Otex  <= 1000
 0  <= F_H2Otpp  <= 1000
 0  <= R_H2Otpp  <= 1000
 0  <= F_HEX1  <= 1000
 0  <= R_HEX1  <= 0
 0  <= F_HISt2rpp  <= 1000
 0  <= R_HISt2rpp  <= 1000
 0  <= F_HIStex  <= 1000
 0  <= R_HIStex  <= 1000
 0  <= F_Htex  <= 1000
 0  <= R_Htex  <= 1000
 0  <= F_ICDHyr  <= 1000
 0  <= R_ICDHyr  <= 1000
 0  <= F_ICL  <= 1000
 0  <= R_ICL  <= 0
 0  <= F_ILEt2rpp  <= 1000
 0  <= R_ILEt2rpp  <= 1000
 0  <= F_ILEtex  <= 1000
 0  <= R_ILEtex  <= 1000
 0  <= F_Kabcpp  <= 1000
 0  <= R_Kabcpp  <= 0
 0  <= F_Kt2pp  <= 1000
 0  <= R_Kt2pp  <= 0
 0  <= F_Ktex  <= 1000
 0  <= R_Ktex  <= 1000
 0  <= F_LDH_D  <= 1000
 0  <= R_LDH_D  <= 1000
 0  <= F_LDH_D2  <= 1000
 0  <= R_LDH_D2  <= 0
 0  <= F_LEUt2rpp  <= 1000
 0  <= R_LEUt2rpp  <= 1000
 0  <= F_LEUtex  <= 1000
 0  <= R_LEUtex  <= 1000
 0  <= F_LYSt3pp  <= 1000
 0  <= R_LYSt3pp  <= 0
 0  <= F_LYStex  <= 1000
 0  <= R_LYStex  <= 1000
 0  <= F_MALS  <= 1000
 0  <= R_MALS  <= 0
 0  <= F_MALt3pp  <= 1000
 0  <= R_MALt3pp  <= 0
 0  <= F_MALtex  <= 1000
 0  <= R_MALtex  <= 1000
 0  <= F_MDH  <= 1000
 0  <= R_MDH  <= 1000
 0  <= F_MDH2  <= 1000
 0  <= R_MDH2  <= 0
 0  <= F_MDH3  <= 1000
 0  <= R_MDH3  <= 0
 0  <= F_ME1  <= 1000
 0  <= R_ME1  <= 0
 0  <= F_ME2  <= 1000
 0  <= R_ME2  <= 0
 0  <= F_MEOHtex  <= 1000
 0  <= R_MEOHtex  <= 1000
 0  <= F_MG2t3_2pp  <= 1000
 0  <= R_MG2t3_2pp  <= 1000
 0  <= F_MG2tex  <= 1000
 0  <= R_MG2tex  <= 1000
 0  <= F_MG2tpp  <= 1000
 0  <= R_MG2tpp  <= 0
 0  <= F_MG2uabcpp  <= 1000
 0  <= R_MG2uabcpp  <= 0
 0  <= F_MN2tpp  <= 1000
 0  <= R_MN2tpp  <= 0
 0  <= F_MNt2pp  <= 1000
 0  <= R_MNt2pp  <= 0
 0  <= F_MNtex  <= 1000
 0  <= R_MNtex  <= 1000
 0  <= F_MOBDabcpp  <= 1000
 0  <= R_MOBDabcpp  <= 0
 0  <= F_MOBDtex  <= 1000
 0  <= R_MOBDtex  <= 1000
 0  <= F_MOX  <= 1000
 0  <= R_MOX  <= 1000
 0  <= F_NADH10  <= 1000
 0  <= R_NADH10  <= 0
 0  <= F_NADH16pp  <= 1000
 0  <= R_NADH16pp  <= 0
 0  <= F_NADH17pp  <= 1000
 0  <= R_NADH17pp  <= 0
 0  <= F_NADH18pp  <= 1000
 0  <= R_NADH18pp  <= 0
 0  <= F_NADH5  <= 1000
 0  <= R_NADH5  <= 0
 0  <= F_NADH9  <= 1000
 0  <= R_NADH9  <= 0
 0  <= F_NADK  <= 1000
 0  <= R_NADK  <= 0
 0  <= F_NADPHQR2  <= 1000
 0  <= R_NADPHQR2  <= 0
 0  <= F_NADPHQR3  <= 1000
 0  <= R_NADPHQR3  <= 0
 0  <= F_NADPHQR4  <= 1000
 0  <= R_NADPHQR4  <= 0
 0  <= F_NADPPPS  <= 1000
 0  <= R_NADPPPS  <= 0
 0  <= F_NADTRHD  <= 1000
 0  <= R_NADTRHD  <= 0
 0  <= F_NAt3_1p5pp  <= 1000
 0  <= R_NAt3_1p5pp  <= 0
 0  <= F_NAt3_2pp  <= 1000
 0  <= R_NAt3_2pp  <= 0
 0  <= F_NAt3pp  <= 1000
 0  <= R_NAt3pp  <= 0
 0  <= F_NH4tex  <= 1000
 0  <= R_NH4tex  <= 1000
 0  <= F_NH4tpp  <= 1000
 0  <= R_NH4tpp  <= 1000
 0  <= F_NI2tex  <= 1000
 0  <= R_NI2tex  <= 1000
 0  <= F_NI2tpp  <= 1000
 0  <= R_NI2tpp  <= 0
 0  <= F_NI2uabcpp  <= 1000
 0  <= R_NI2uabcpp  <= 0
 0  <= F_NTP1  <= 1000
 0  <= R_NTP1  <= 0
 0  <= F_NTPP6  <= 1000
 0  <= R_NTPP6  <= 0
 0  <= F_O2tex  <= 1000
 0  <= R_O2tex  <= 1000
 0  <= F_O2tpp  <= 1000
 0  <= R_O2tpp  <= 1000
 0  <= F_OAADC  <= 1000
 0  <= R_OAADC  <= 0
 0  <= F_PDH  <= 1000
 0  <= R_PDH  <= 0
 0  <= F_PE160abcpp  <= 1000
 0  <= R_PE160abcpp  <= 0
 0  <= F_PE161abcpp  <= 1000
 0  <= R_PE161abcpp  <= 0
 0  <= F_PE181abcpp  <= 1000
 0  <= R_PE181abcpp  <= 0
 0  <= F_PFK  <= 1000
 0  <= R_PFK  <= 0
 0  <= F_PFK_3  <= 1000
 0  <= R_PFK_3  <= 0
 0  <= F_PFL  <= 1000
 0  <= R_PFL  <= 0
 0  <= F_PG160abcpp  <= 1000
 0  <= R_PG160abcpp  <= 0
 0  <= F_PG161abcpp  <= 1000
 0  <= R_PG161abcpp  <= 0
 0  <= F_PG181abcpp  <= 1000
 0  <= R_PG181abcpp  <= 0
 0  <= F_PGI  <= 1000
 0  <= R_PGI  <= 1000
 0  <= F_PGK  <= 1000
 0  <= R_PGK  <= 1000
 0  <= F_PGL  <= 1000
 0  <= R_PGL  <= 0
 0  <= F_PGM  <= 1000
 0  <= R_PGM  <= 1000
 0  <= F_PGMT  <= 1000
 0  <= R_PGMT  <= 1000
 0  <= F_PHEMEabcpp  <= 1000
 0  <= R_PHEMEabcpp  <= 0
 0  <= F_PHEMEtiex  <= 1000
 0  <= R_PHEMEtiex  <= 0
 0  <= F_PHEt2rpp  <= 1000
 0  <= R_PHEt2rpp  <= 1000
 0  <= F_PHEtex  <= 1000
 0  <= R_PHEtex  <= 1000
 0  <= F_PIt2rpp  <= 1000
 0  <= R_PIt2rpp  <= 1000
 0  <= F_PItex  <= 1000
 0  <= R_PItex  <= 1000
 0  <= F_PIuabcpp  <= 1000
 0  <= R_PIuabcpp  <= 0
 0  <= F_POR5  <= 1000
 0  <= R_POR5  <= 1000
 0  <= F_POX  <= 1000
 0  <= R_POX  <= 0
 0  <= F_PPA  <= 1000
 0  <= R_PPA  <= 0
 0  <= F_PPA2  <= 1000
 0  <= R_PPA2  <= 0
 0  <= F_PPC  <= 1000
 0  <= R_PPC  <= 0
 0  <= F_PPCK  <= 1000
 0  <= R_PPCK  <= 0
 0  <= F_PPK2r  <= 1000
 0  <= R_PPK2r  <= 1000
 0  <= F_PPKr  <= 1000
 0  <= R_PPKr  <= 1000
 0  <= F_PPPGO  <= 1000
 0  <= R_PPPGO  <= 0
 0  <= F_PPPGO3  <= 1000
 0  <= R_PPPGO3  <= 0
 0  <= F_PPS  <= 1000
 0  <= R_PPS  <= 0
 0  <= F_PROt2rpp  <= 1000
 0  <= R_PROt2rpp  <= 1000
 0  <= F_PROtex  <= 1000
 0  <= R_PROtex  <= 1000
 0  <= F_PTAr  <= 1000
 0  <= R_PTAr  <= 1000
 0  <= F_PYK  <= 1000
 0  <= R_PYK  <= 0
 0  <= F_PYRt2rpp  <= 1000
 0  <= R_PYRt2rpp  <= 1000
 0  <= F_PYRtex  <= 1000
 0  <= R_PYRtex  <= 1000
 0  <= F_RPE  <= 1000
 0  <= R_RPE  <= 1000
 0  <= F_RPI  <= 1000
 0  <= R_RPI  <= 1000
 0  <= F_SERt2rpp  <= 1000
 0  <= R_SERt2rpp  <= 1000
 0  <= F_SERtex  <= 1000
 0  <= R_SERtex  <= 1000
 0  <= F_SO4t2pp  <= 1000
 0  <= R_SO4t2pp  <= 0
 0  <= F_SO4tex  <= 1000
 0  <= R_SO4tex  <= 1000
 0  <= F_SPMDt3pp  <= 1000
 0  <= R_SPMDt3pp  <= 0
 0  <= F_SPMDtex  <= 1000
 0  <= R_SPMDtex  <= 1000
 0  <= F_SUCCt2_2pp  <= 1000
 0  <= R_SUCCt2_2pp  <= 0
 0  <= F_SUCCt2_3pp  <= 1000
 0  <= R_SUCCt2_3pp  <= 0
 0  <= F_SUCCt3pp  <= 1000
 0  <= R_SUCCt3pp  <= 0
 0  <= F_SUCCtex  <= 1000
 0  <= R_SUCCtex  <= 1000
 0  <= F_SUCDi  <= 1000
 0  <= R_SUCDi  <= 0
 0  <= F_SUCOAS  <= 1000
 0  <= R_SUCOAS  <= 1000
 0  <= F_SULabcpp  <= 1000
 0  <= R_SULabcpp  <= 0
 0  <= F_TALA  <= 1000
 0  <= R_TALA  <= 1000
 0  <= F_THD2pp  <= 1000
 0  <= R_THD2pp  <= 0
 0  <= F_THIORDXi  <= 1000
 0  <= R_THIORDXi  <= 0
 0  <= F_THRt2pp  <= 1000
 0  <= R_THRt2pp  <= 0
 0  <= F_THRt2rpp  <= 1000
 0  <= R_THRt2rpp  <= 1000
 0  <= F_THRtex  <= 1000
 0  <= R_THRtex  <= 1000
 0  <= F_TKT1  <= 1000
 0  <= R_TKT1  <= 1000
 0  <= F_TKT2  <= 1000
 0  <= R_TKT2  <= 1000
 0  <= F_TPI  <= 1000
 0  <= R_TPI  <= 1000
 0  <= F_TRDR  <= 1000
 0  <= R_TRDR  <= 0
 0  <= F_TRPt2rpp  <= 1000
 0  <= R_TRPt2rpp  <= 1000
 0  <= F_TRPtex  <= 1000
 0  <= R_TRPtex  <= 1000
 0  <= F_TRSARr  <= 1000
 0  <= R_TRSARr  <= 1000
 0  <= F_TYRt2rpp  <= 1000
 0  <= R_TYRt2rpp  <= 1000
 0  <= F_TYRtex  <= 1000
 0  <= R_TYRtex  <= 1000
 0  <= F_VALt2rpp  <= 1000
 0  <= R_VALt2rpp  <= 1000
 0  <= F_VALtex  <= 1000
 0  <= R_VALtex  <= 1000
 0  <= F_ZN2tpp  <= 1000
 0  <= R_ZN2tpp  <= 0
 0  <= F_ZNabcpp  <= 1000
 0  <= R_ZNabcpp  <= 0
 0  <= F_Zn2tex  <= 1000
 0  <= R_Zn2tex  <= 1000
 0  <= F_LMPD_1_10fthf_c  <= 1000
 0  <= R_LMPD_1_10fthf_c  <= 0
 0  <= F_LMPD_2_10fthf_c  <= 1000
 0  <= R_LMPD_2_10fthf_c  <= 0
 0  <= F_LMPD_3_10fthf_c  <= 1000
 0  <= R_LMPD_3_10fthf_c  <= 0
 0  <= F_LMPD_4_10fthf_c  <= 1000
 0  <= R_LMPD_4_10fthf_c  <= 0
 0  <= F_LMPD_5_2dmmql8_c  <= 1000
 0  <= R_LMPD_5_2dmmql8_c  <= 0
 0  <= F_LMPD_6_2fe2s_c  <= 1000
 0  <= R_LMPD_6_2fe2s_c  <= 0
 0  <= F_LMPD_7_2fe2s_c  <= 1000
 0  <= R_LMPD_7_2fe2s_c  <= 0
 0  <= F_LMPD_8_2fe2s_c  <= 1000
 0  <= R_LMPD_8_2fe2s_c  <= 0
 0  <= F_LMPD_9_2fe2s_c  <= 1000
 0  <= R_LMPD_9_2fe2s_c  <= 0
 0  <= F_LMPD_10_2fe2s_c  <= 1000
 0  <= R_LMPD_10_2fe2s_c  <= 0
 0  <= F_LMPD_11_2fe2s_c  <= 1000
 0  <= R_LMPD_11_2fe2s_c  <= 0
 0  <= F_LMPD_12_2fe2s_c  <= 1000
 0  <= R_LMPD_12_2fe2s_c  <= 0
 0  <= F_LMPD_13_2fe2s_c  <= 1000
 0  <= R_LMPD_13_2fe2s_c  <= 0
 0  <= F_LMPD_14_4fe4s_c  <= 1000
 0  <= R_LMPD_14_4fe4s_c  <= 0
 0  <= F_LMPD_15_4fe4s_c  <= 1000
 0  <= R_LMPD_15_4fe4s_c  <= 0
 0  <= F_LMPD_16_4fe4s_c  <= 1000
 0  <= R_LMPD_16_4fe4s_c  <= 0
 0  <= F_LMPD_17_4fe4s_c  <= 1000
 0  <= R_LMPD_17_4fe4s_c  <= 0
 0  <= F_LMPD_18_4fe4s_c  <= 1000
 0  <= R_LMPD_18_4fe4s_c  <= 0
 0  <= F_LMPD_19_4fe4s_c  <= 1000
 0  <= R_LMPD_19_4fe4s_c  <= 0
 0  <= F_LMPD_20_4fe4s_c  <= 1000
 0  <= R_LMPD_20_4fe4s_c  <= 0
 0  <= F_LMPD_21_4fe4s_c  <= 1000
 0  <= R_LMPD_21_4fe4s_c  <= 0
 0  <= F_LMPD_22_5mthf_c  <= 1000
 0  <= R_LMPD_22_5mthf_c  <= 0
 0  <= F_LMPD_23_5mthf_c  <= 1000
 0  <= R_LMPD_23_5mthf_c  <= 0
 0  <= F_LMPD_24_accoa_c  <= 1000
 0  <= R_LMPD_24_accoa_c  <= 0
 0  <= F_LMPD_25_accoa_c  <= 1000
 0  <= R_LMPD_25_accoa_c  <= 0
 0  <= F_LMPD_26_adocbl_c  <= 1000
 0  <= R_LMPD_26_adocbl_c  <= 0
 0  <= F_LMPD_27_adocbl_c  <= 1000
 0  <= R_LMPD_27_adocbl_c  <= 0
 0  <= F_LMPD_28_adocbl_c  <= 1000
 0  <= R_LMPD_28_adocbl_c  <= 0
 0  <= F_LMPD_29_adocbl_c  <= 1000
 0  <= R_LMPD_29_adocbl_c  <= 0
 0  <= F_LMPD_30_ala-L_c  <= 1000
 0  <= R_LMPD_30_ala-L_c  <= 0
 0  <= F_LMPD_31_amet_c  <= 1000
 0  <= R_LMPD_31_amet_c  <= 0
 0  <= F_LMPD_32_amet_c  <= 1000
 0  <= R_LMPD_32_amet_c  <= 0
 0  <= F_LMPD_33_arg-L_c  <= 1000
 0  <= R_LMPD_33_arg-L_c  <= 0
 0  <= F_LMPD_34_asn-L_c  <= 1000
 0  <= R_LMPD_34_asn-L_c  <= 0
 0  <= F_LMPD_35_bmocogdp_c  <= 1000
 0  <= R_LMPD_35_bmocogdp_c  <= 0
 0  <= F_LMPD_36_bmocogdp_c  <= 1000
 0  <= R_LMPD_36_bmocogdp_c  <= 0
 0  <= F_LMPD_37_bmocogdp_c  <= 1000
 0  <= R_LMPD_37_bmocogdp_c  <= 0
 0  <= F_LMPD_38_bmocogdp_c  <= 1000
 0  <= R_LMPD_38_bmocogdp_c  <= 0
 0  <= F_LMPD_39_btn_c  <= 1000
 0  <= R_LMPD_39_btn_c  <= 0
 0  <= F_LMPD_40_btn_c  <= 1000
 0  <= R_LMPD_40_btn_c  <= 0
 0  <= F_LMPD_41_btn_c  <= 1000
 0  <= R_LMPD_41_btn_c  <= 0
 0  <= F_LMPD_42_btn_c  <= 1000
 0  <= R_LMPD_42_btn_c  <= 0
 0  <= F_LMPD_43_btn_c  <= 1000
 0  <= R_LMPD_43_btn_c  <= 0
 0  <= F_LMPD_44_btn_c  <= 1000
 0  <= R_LMPD_44_btn_c  <= 0
 0  <= F_LMPD_45_btn_c  <= 1000
 0  <= R_LMPD_45_btn_c  <= 0
 0  <= F_LMPD_47_btn_c  <= 1000
 0  <= R_LMPD_47_btn_c  <= 0
 0  <= F_LMPD_48_ca2_c  <= 1000
 0  <= R_LMPD_48_ca2_c  <= 0
 0  <= F_LMPD_49_ca2_c  <= 1000
 0  <= R_LMPD_49_ca2_c  <= 0
 0  <= F_LMPD_50_ca2_c  <= 1000
 0  <= R_LMPD_50_ca2_c  <= 0
 0  <= F_LMPD_51_chor_c  <= 1000
 0  <= R_LMPD_51_chor_c  <= 0
 0  <= F_LMPD_52_clpn160_p  <= 1000
 0  <= R_LMPD_52_clpn160_p  <= 0
 0  <= F_LMPD_53_clpn160_p  <= 1000
 0  <= R_LMPD_53_clpn160_p  <= 0
 0  <= F_LMPD_54_clpn160_p  <= 1000
 0  <= R_LMPD_54_clpn160_p  <= 0
 0  <= F_LMPD_55_clpn160_p  <= 1000
 0  <= R_LMPD_55_clpn160_p  <= 0
 0  <= F_LMPD_56_clpn160_p  <= 1000
 0  <= R_LMPD_56_clpn160_p  <= 0
 0  <= F_LMPD_57_clpn160_p  <= 1000
 0  <= R_LMPD_57_clpn160_p  <= 0
 0  <= F_LMPD_58_clpn160_p  <= 1000
 0  <= R_LMPD_58_clpn160_p  <= 0
 0  <= F_LMPD_59_clpn160_p  <= 1000
 0  <= R_LMPD_59_clpn160_p  <= 0
 0  <= F_LMPD_60_clpn160_p  <= 1000
 0  <= R_LMPD_60_clpn160_p  <= 0
 0  <= F_LMPD_61_clpn160_p  <= 1000
 0  <= R_LMPD_61_clpn160_p  <= 0
 0  <= F_LMPD_62_clpn161_p  <= 1000
 0  <= R_LMPD_62_clpn161_p  <= 0
 0  <= F_LMPD_63_clpn161_p  <= 1000
 0  <= R_LMPD_63_clpn161_p  <= 0
 0  <= F_LMPD_64_clpn161_p  <= 1000
 0  <= R_LMPD_64_clpn161_p  <= 0
 0  <= F_LMPD_65_clpn161_p  <= 1000
 0  <= R_LMPD_65_clpn161_p  <= 0
 0  <= F_LMPD_66_clpn161_p  <= 1000
 0  <= R_LMPD_66_clpn161_p  <= 0
 0  <= F_LMPD_67_clpn161_p  <= 1000
 0  <= R_LMPD_67_clpn161_p  <= 0
 0  <= F_LMPD_68_clpn161_p  <= 1000
 0  <= R_LMPD_68_clpn161_p  <= 0
 0  <= F_LMPD_69_clpn181_p  <= 1000
 0  <= R_LMPD_69_clpn181_p  <= 0
 0  <= F_LMPD_70_clpn181_p  <= 1000
 0  <= R_LMPD_70_clpn181_p  <= 0
 0  <= F_LMPD_71_clpn181_p  <= 1000
 0  <= R_LMPD_71_clpn181_p  <= 0
 0  <= F_LMPD_72_clpn181_p  <= 1000
 0  <= R_LMPD_72_clpn181_p  <= 0
 0  <= F_LMPD_73_clpn181_p  <= 1000
 0  <= R_LMPD_73_clpn181_p  <= 0
 0  <= F_LMPD_74_clpn181_p  <= 1000
 0  <= R_LMPD_74_clpn181_p  <= 0
 0  <= F_LMPD_75_clpn181_p  <= 1000
 0  <= R_LMPD_75_clpn181_p  <= 0
 0  <= F_LMPD_76_clpn181_p  <= 1000
 0  <= R_LMPD_76_clpn181_p  <= 0
 0  <= F_LMPD_77_clpn181_p  <= 1000
 0  <= R_LMPD_77_clpn181_p  <= 0
 0  <= F_LMPD_78_clpn181_p  <= 1000
 0  <= R_LMPD_78_clpn181_p  <= 0
 0  <= F_LMPD_79_colipa_e  <= 1000
 0  <= R_LMPD_79_colipa_e  <= 0
 0  <= F_LMPD_80_colipa_e  <= 1000
 0  <= R_LMPD_80_colipa_e  <= 0
 0  <= F_LMPD_81_colipa_e  <= 1000
 0  <= R_LMPD_81_colipa_e  <= 0
 0  <= F_LMPD_82_colipa_e  <= 1000
 0  <= R_LMPD_82_colipa_e  <= 0
 0  <= F_LMPD_83_colipa_e  <= 1000
 0  <= R_LMPD_83_colipa_e  <= 0
 0  <= F_LMPD_84_colipa_e  <= 1000
 0  <= R_LMPD_84_colipa_e  <= 0
 0  <= F_LMPD_85_colipa_e  <= 1000
 0  <= R_LMPD_85_colipa_e  <= 0
 0  <= F_LMPD_86_colipa_e  <= 1000
 0  <= R_LMPD_86_colipa_e  <= 0
 0  <= F_LMPD_87_colipa_e  <= 1000
 0  <= R_LMPD_87_colipa_e  <= 0
 0  <= F_LMPD_88_colipa_e  <= 1000
 0  <= R_LMPD_88_colipa_e  <= 0
 0  <= F_LMPD_89_colipa_e  <= 1000
 0  <= R_LMPD_89_colipa_e  <= 0
 0  <= F_LMPD_90_colipa_e  <= 1000
 0  <= R_LMPD_90_colipa_e  <= 0
 0  <= F_LMPD_91_ctp_c  <= 1000
 0  <= R_LMPD_91_ctp_c  <= 0
 0  <= F_LMPD_92_cys-L_c  <= 1000
 0  <= R_LMPD_92_cys-L_c  <= 0
 0  <= F_LMPD_93_cys-L_c  <= 1000
 0  <= R_LMPD_93_cys-L_c  <= 0
 0  <= F_LMPD_94_datp_c  <= 1000
 0  <= R_LMPD_94_datp_c  <= 0
 0  <= F_LMPD_95_datp_c  <= 1000
 0  <= R_LMPD_95_datp_c  <= 0
 0  <= F_LMPD_96_datp_c  <= 1000
 0  <= R_LMPD_96_datp_c  <= 0
 0  <= F_LMPD_97_datp_c  <= 1000
 0  <= R_LMPD_97_datp_c  <= 0
 0  <= F_LMPD_98_dctp_c  <= 1000
 0  <= R_LMPD_98_dctp_c  <= 0
 0  <= F_LMPD_99_dgtp_c  <= 1000
 0  <= R_LMPD_99_dgtp_c  <= 0
 0  <= F_LMPD_100_dgtp_c  <= 1000
 0  <= R_LMPD_100_dgtp_c  <= 0
 0  <= F_LMPD_101_dgtp_c  <= 1000
 0  <= R_LMPD_101_dgtp_c  <= 0
 0  <= F_LMPD_102_dgtp_c  <= 1000
 0  <= R_LMPD_102_dgtp_c  <= 0
 0  <= F_LMPD_103_dgtp_c  <= 1000
 0  <= R_LMPD_103_dgtp_c  <= 0
 0  <= F_LMPD_104_dgtp_c  <= 1000
 0  <= R_LMPD_104_dgtp_c  <= 0
 0  <= F_LMPD_105_dttp_c  <= 1000
 0  <= R_LMPD_105_dttp_c  <= 0
 0  <= F_LMPD_106_dttp_c  <= 1000
 0  <= R_LMPD_106_dttp_c  <= 0
 0  <= F_LMPD_107_enter_c  <= 1000
 0  <= R_LMPD_107_enter_c  <= 0
 0  <= F_LMPD_108_fad_c  <= 1000
 0  <= R_LMPD_108_fad_c  <= 0
 0  <= F_LMPD_109_fad_c  <= 1000
 0  <= R_LMPD_109_fad_c  <= 0
 0  <= F_LMPD_110_fad_c  <= 1000
 0  <= R_LMPD_110_fad_c  <= 0
 0  <= F_LMPD_111_fad_c  <= 1000
 0  <= R_LMPD_111_fad_c  <= 0
 0  <= F_LMPD_112_GAM_c  <= 1000
 0  <= R_LMPD_112_GAM_c  <= 0
 0  <= F_LMPD_113_GAM_c  <= 1000
 0  <= R_LMPD_113_GAM_c  <= 0
 0  <= F_LMPD_114_GAM_c  <= 1000
 0  <= R_LMPD_114_GAM_c  <= 0
 0  <= F_LMPD_115_GAM_c  <= 1000
 0  <= R_LMPD_115_GAM_c  <= 0
 0  <= F_LMPD_116_gln-L_c  <= 1000
 0  <= R_LMPD_116_gln-L_c  <= 0
 0  <= F_LMPD_117_gly_c  <= 1000
 0  <= R_LMPD_117_gly_c  <= 0
 0  <= F_LMPD_118_gthrd_c  <= 1000
 0  <= R_LMPD_118_gthrd_c  <= 0
 0  <= F_LMPD_119_gthrd_c  <= 1000
 0  <= R_LMPD_119_gthrd_c  <= 0
 0  <= F_LMPD_120_gtp_c  <= 1000
 0  <= R_LMPD_120_gtp_c  <= 0
 0  <= F_LMPD_121_gtp_c  <= 1000
 0  <= R_LMPD_121_gtp_c  <= 0
 0  <= F_LMPD_122_gtp_c  <= 1000
 0  <= R_LMPD_122_gtp_c  <= 0
 0  <= F_LMPD_123_gtp_c  <= 1000
 0  <= R_LMPD_123_gtp_c  <= 0
 0  <= F_LMPD_124_hemeO_c  <= 1000
 0  <= R_LMPD_124_hemeO_c  <= 0
 0  <= F_LMPD_125_hemeO_c  <= 1000
 0  <= R_LMPD_125_hemeO_c  <= 0
 0  <= F_LMPD_126_his-L_c  <= 1000
 0  <= R_LMPD_126_his-L_c  <= 0
 0  <= F_LMPD_127_his-L_c  <= 1000
 0  <= R_LMPD_127_his-L_c  <= 0
 0  <= F_LMPD_128_ile-L_c  <= 1000
 0  <= R_LMPD_128_ile-L_c  <= 0
 0  <= F_LMPD_129_leu-L_c  <= 1000
 0  <= R_LMPD_129_leu-L_c  <= 0
 0  <= F_LMPD_130_lipopb_c  <= 1000
 0  <= R_LMPD_130_lipopb_c  <= 0
 0  <= F_LMPD_131_lipopb_c  <= 1000
 0  <= R_LMPD_131_lipopb_c  <= 0
 0  <= F_LMPD_132_lipopb_c  <= 1000
 0  <= R_LMPD_132_lipopb_c  <= 0
 0  <= F_LMPD_133_lipopb_c  <= 1000
 0  <= R_LMPD_133_lipopb_c  <= 0
 0  <= F_LMPD_134_lys-L_c  <= 1000
 0  <= R_LMPD_134_lys-L_c  <= 0
 0  <= F_LMPD_135_malcoa_c  <= 1000
 0  <= R_LMPD_135_malcoa_c  <= 0
 0  <= F_LMPD_136_malcoa_c  <= 1000
 0  <= R_LMPD_136_malcoa_c  <= 0
 0  <= F_LMPD_137_met-L_c  <= 1000
 0  <= R_LMPD_137_met-L_c  <= 0
 0  <= F_LMPD_138_met-L_c  <= 1000
 0  <= R_LMPD_138_met-L_c  <= 0
 0  <= F_LMPD_139_mlthf_c  <= 1000
 0  <= R_LMPD_139_mlthf_c  <= 0
 0  <= F_LMPD_140_mlthf_c  <= 1000
 0  <= R_LMPD_140_mlthf_c  <= 0
 0  <= F_LMPD_141_mococdp_c  <= 1000
 0  <= R_LMPD_141_mococdp_c  <= 0
 0  <= F_LMPD_142_mococdp_c  <= 1000
 0  <= R_LMPD_142_mococdp_c  <= 0
 0  <= F_LMPD_143_mococdp_c  <= 1000
 0  <= R_LMPD_143_mococdp_c  <= 0
 0  <= F_LMPD_144_mococdp_c  <= 1000
 0  <= R_LMPD_144_mococdp_c  <= 0
 0  <= F_LMPD_145_mococdp_c  <= 1000
 0  <= R_LMPD_145_mococdp_c  <= 0
 0  <= F_LMPD_146_mococdp_c  <= 1000
 0  <= R_LMPD_146_mococdp_c  <= 0
 0  <= F_LMPD_147_mococdp_c  <= 1000
 0  <= R_LMPD_147_mococdp_c  <= 0
 0  <= F_LMPD_148_mococdp_c  <= 1000
 0  <= R_LMPD_148_mococdp_c  <= 0
 0  <= F_LMPD_149_mocogdp_c  <= 1000
 0  <= R_LMPD_149_mocogdp_c  <= 0
 0  <= F_LMPD_150_mocogdp_c  <= 1000
 0  <= R_LMPD_150_mocogdp_c  <= 0
 0  <= F_LMPD_151_mocogdp_c  <= 1000
 0  <= R_LMPD_151_mocogdp_c  <= 0
 0  <= F_LMPD_152_mocogdp_c  <= 1000
 0  <= R_LMPD_152_mocogdp_c  <= 0
 0  <= F_LMPD_153_mql8_c  <= 1000
 0  <= R_LMPD_153_mql8_c  <= 0
 0  <= F_LMPD_154_murein3p3p_p  <= 1000
 0  <= R_LMPD_154_murein3p3p_p  <= 0
 0  <= F_LMPD_155_murein3px4p_p  <= 1000
 0  <= R_LMPD_155_murein3px4p_p  <= 0
 0  <= F_LMPD_156_murein4p4p_p  <= 1000
 0  <= R_LMPD_156_murein4p4p_p  <= 0
 0  <= F_LMPD_157_murein4px4p_p  <= 1000
 0  <= R_LMPD_157_murein4px4p_p  <= 0
 0  <= F_LMPD_158_murein4px4px4p_p  <= 1000
 0  <= R_LMPD_158_murein4px4px4p_p  <= 0
 0  <= F_LMPD_159_nad_c  <= 1000
 0  <= R_LMPD_159_nad_c  <= 0
 0  <= F_LMPD_160_nad_c  <= 1000
 0  <= R_LMPD_160_nad_c  <= 0
 0  <= F_LMPD_161_nad_c  <= 1000
 0  <= R_LMPD_161_nad_c  <= 0
 0  <= F_LMPD_162_nad_c  <= 1000
 0  <= R_LMPD_162_nad_c  <= 0
 0  <= F_LMPD_163_pe160_c  <= 1000
 0  <= R_LMPD_163_pe160_c  <= 0
 0  <= F_LMPD_164_pe160_c  <= 1000
 0  <= R_LMPD_164_pe160_c  <= 0
 0  <= F_LMPD_165_pe160_c  <= 1000
 0  <= R_LMPD_165_pe160_c  <= 0
 0  <= F_LMPD_166_pe160_c  <= 1000
 0  <= R_LMPD_166_pe160_c  <= 0
 0  <= F_LMPD_167_pe160_c  <= 1000
 0  <= R_LMPD_167_pe160_c  <= 0
 0  <= F_LMPD_168_pe160_c  <= 1000
 0  <= R_LMPD_168_pe160_c  <= 0
 0  <= F_LMPD_169_pe160_c  <= 1000
 0  <= R_LMPD_169_pe160_c  <= 0
 0  <= F_LMPD_170_pe160_c  <= 1000
 0  <= R_LMPD_170_pe160_c  <= 0
 0  <= F_LMPD_171_pe160_c  <= 1000
 0  <= R_LMPD_171_pe160_c  <= 0
 0  <= F_LMPD_172_pe160_c  <= 1000
 0  <= R_LMPD_172_pe160_c  <= 0
 0  <= F_LMPD_173_pe161_c  <= 1000
 0  <= R_LMPD_173_pe161_c  <= 0
 0  <= F_LMPD_174_pe161_c  <= 1000
 0  <= R_LMPD_174_pe161_c  <= 0
 0  <= F_LMPD_175_pe161_c  <= 1000
 0  <= R_LMPD_175_pe161_c  <= 0
 0  <= F_LMPD_176_pe161_c  <= 1000
 0  <= R_LMPD_176_pe161_c  <= 0
 0  <= F_LMPD_177_pe161_c  <= 1000
 0  <= R_LMPD_177_pe161_c  <= 0
 0  <= F_LMPD_178_pe161_c  <= 1000
 0  <= R_LMPD_178_pe161_c  <= 0
 0  <= F_LMPD_179_pe161_c  <= 1000
 0  <= R_LMPD_179_pe161_c  <= 0
 0  <= F_LMPD_180_pe181_c  <= 1000
 0  <= R_LMPD_180_pe181_c  <= 0
 0  <= F_LMPD_181_pe181_c  <= 1000
 0  <= R_LMPD_181_pe181_c  <= 0
 0  <= F_LMPD_182_pe181_c  <= 1000
 0  <= R_LMPD_182_pe181_c  <= 0
 0  <= F_LMPD_183_pe181_c  <= 1000
 0  <= R_LMPD_183_pe181_c  <= 0
 0  <= F_LMPD_184_pe181_c  <= 1000
 0  <= R_LMPD_184_pe181_c  <= 0
 0  <= F_LMPD_185_pe181_c  <= 1000
 0  <= R_LMPD_185_pe181_c  <= 0
 0  <= F_LMPD_186_pe181_c  <= 1000
 0  <= R_LMPD_186_pe181_c  <= 0
 0  <= F_LMPD_187_pe181_c  <= 1000
 0  <= R_LMPD_187_pe181_c  <= 0
 0  <= F_LMPD_188_pe181_c  <= 1000
 0  <= R_LMPD_188_pe181_c  <= 0
 0  <= F_LMPD_189_pe181_c  <= 1000
 0  <= R_LMPD_189_pe181_c  <= 0
 0  <= F_LMPD_190_pg160_c  <= 1000
 0  <= R_LMPD_190_pg160_c  <= 0
 0  <= F_LMPD_191_pg160_c  <= 1000
 0  <= R_LMPD_191_pg160_c  <= 0
 0  <= F_LMPD_192_pg160_c  <= 1000
 0  <= R_LMPD_192_pg160_c  <= 0
 0  <= F_LMPD_193_pg160_c  <= 1000
 0  <= R_LMPD_193_pg160_c  <= 0
 0  <= F_LMPD_194_pg160_c  <= 1000
 0  <= R_LMPD_194_pg160_c  <= 0
 0  <= F_LMPD_195_pg160_c  <= 1000
 0  <= R_LMPD_195_pg160_c  <= 0
 0  <= F_LMPD_196_pg160_c  <= 1000
 0  <= R_LMPD_196_pg160_c  <= 0
 0  <= F_LMPD_197_pg160_c  <= 1000
 0  <= R_LMPD_197_pg160_c  <= 0
 0  <= F_LMPD_198_pg160_c  <= 1000
 0  <= R_LMPD_198_pg160_c  <= 0
 0  <= F_LMPD_199_pg160_c  <= 1000
 0  <= R_LMPD_199_pg160_c  <= 0
 0  <= F_LMPD_200_pg161_c  <= 1000
 0  <= R_LMPD_200_pg161_c  <= 0
 0  <= F_LMPD_201_pg161_c  <= 1000
 0  <= R_LMPD_201_pg161_c  <= 0
 0  <= F_LMPD_202_pg161_c  <= 1000
 0  <= R_LMPD_202_pg161_c  <= 0
 0  <= F_LMPD_203_pg161_c  <= 1000
 0  <= R_LMPD_203_pg161_c  <= 0
 0  <= F_LMPD_204_pg161_c  <= 1000
 0  <= R_LMPD_204_pg161_c  <= 0
 0  <= F_LMPD_205_pg161_c  <= 1000
 0  <= R_LMPD_205_pg161_c  <= 0
 0  <= F_LMPD_206_pg161_c  <= 1000
 0  <= R_LMPD_206_pg161_c  <= 0
 0  <= F_LMPD_207_pg181_c  <= 1000
 0  <= R_LMPD_207_pg181_c  <= 0
 0  <= F_LMPD_208_pg181_c  <= 1000
 0  <= R_LMPD_208_pg181_c  <= 0
 0  <= F_LMPD_209_pg181_c  <= 1000
 0  <= R_LMPD_209_pg181_c  <= 0
 0  <= F_LMPD_210_pg181_c  <= 1000
 0  <= R_LMPD_210_pg181_c  <= 0
 0  <= F_LMPD_211_pg181_c  <= 1000
 0  <= R_LMPD_211_pg181_c  <= 0
 0  <= F_LMPD_212_pg181_c  <= 1000
 0  <= R_LMPD_212_pg181_c  <= 0
 0  <= F_LMPD_213_pg181_c  <= 1000
 0  <= R_LMPD_213_pg181_c  <= 0
 0  <= F_LMPD_214_pg181_c  <= 1000
 0  <= R_LMPD_214_pg181_c  <= 0
 0  <= F_LMPD_215_pg181_c  <= 1000
 0  <= R_LMPD_215_pg181_c  <= 0
 0  <= F_LMPD_216_pg181_c  <= 1000
 0  <= R_LMPD_216_pg181_c  <= 0
 0  <= F_LMPD_217_phe-L_c  <= 1000
 0  <= R_LMPD_217_phe-L_c  <= 0
 0  <= F_LMPD_218_pheme_c  <= 1000
 0  <= R_LMPD_218_pheme_c  <= 0
 0  <= F_LMPD_219_pheme_c  <= 1000
 0  <= R_LMPD_219_pheme_c  <= 0
 0  <= F_LMPD_220_pro-L_c  <= 1000
 0  <= R_LMPD_220_pro-L_c  <= 0
 0  <= F_LMPD_221_ptrc_c  <= 1000
 0  <= R_LMPD_221_ptrc_c  <= 0
 0  <= F_LMPD_222_pydx5p_c  <= 1000
 0  <= R_LMPD_222_pydx5p_c  <= 0
 0  <= F_LMPD_223_pydx5p_c  <= 1000
 0  <= R_LMPD_223_pydx5p_c  <= 0
 0  <= F_LMPD_224_q8h2_c  <= 1000
 0  <= R_LMPD_224_q8h2_c  <= 0
 0  <= F_LMPD_225_q8h2_c  <= 1000
 0  <= R_LMPD_225_q8h2_c  <= 0
 0  <= F_LMPD_226_q8h2_c  <= 1000
 0  <= R_LMPD_226_q8h2_c  <= 0
 0  <= F_LMPD_227_q8h2_c  <= 1000
 0  <= R_LMPD_227_q8h2_c  <= 0
 0  <= F_LMPD_228_q8h2_c  <= 1000
 0  <= R_LMPD_228_q8h2_c  <= 0
 0  <= F_LMPD_229_q8h2_c  <= 1000
 0  <= R_LMPD_229_q8h2_c  <= 0
 0  <= F_LMPD_230_q8h2_c  <= 1000
 0  <= R_LMPD_230_q8h2_c  <= 0
 0  <= F_LMPD_231_q8h2_c  <= 1000
 0  <= R_LMPD_231_q8h2_c  <= 0
 0  <= F_LMPD_232_q8h2_c  <= 1000
 0  <= R_LMPD_232_q8h2_c  <= 0
 0  <= F_LMPD_233_q8h2_c  <= 1000
 0  <= R_LMPD_233_q8h2_c  <= 0
 0  <= F_LMPD_234_q8h2_c  <= 1000
 0  <= R_LMPD_234_q8h2_c  <= 0
 0  <= F_LMPD_235_ribflv_c  <= 1000
 0  <= R_LMPD_235_ribflv_c  <= 0
 0  <= F_LMPD_236_ribflv_c  <= 1000
 0  <= R_LMPD_236_ribflv_c  <= 0
 0  <= F_LMPD_237_ribflv_c  <= 1000
 0  <= R_LMPD_237_ribflv_c  <= 0
 0  <= F_LMPD_238_ribflv_c  <= 1000
 0  <= R_LMPD_238_ribflv_c  <= 0
 0  <= F_LMPD_239_ser-L_c  <= 1000
 0  <= R_LMPD_239_ser-L_c  <= 0
 0  <= F_LMPD_240_sheme_c  <= 1000
 0  <= R_LMPD_240_sheme_c  <= 0
 0  <= F_LMPD_241_spmd_c  <= 1000
 0  <= R_LMPD_241_spmd_c  <= 0
 0  <= F_LMPD_242_spmd_c  <= 1000
 0  <= R_LMPD_242_spmd_c  <= 0
 0  <= F_LMPD_243_thf_c  <= 1000
 0  <= R_LMPD_243_thf_c  <= 0
 0  <= F_LMPD_244_thf_c  <= 1000
 0  <= R_LMPD_244_thf_c  <= 0
 0  <= F_LMPD_245_thf_c  <= 1000
 0  <= R_LMPD_245_thf_c  <= 0
 0  <= F_LMPD_246_thf_c  <= 1000
 0  <= R_LMPD_246_thf_c  <= 0
 0  <= F_LMPD_247_thmpp_c  <= 1000
 0  <= R_LMPD_247_thmpp_c  <= 0
 0  <= F_LMPD_248_thmpp_c  <= 1000
 0  <= R_LMPD_248_thmpp_c  <= 0
 0  <= F_LMPD_249_thr-L_c  <= 1000
 0  <= R_LMPD_249_thr-L_c  <= 0
 0  <= F_LMPD_250_trp-L_c  <= 1000
 0  <= R_LMPD_250_trp-L_c  <= 0
 0  <= F_LMPD_251_tyr-L_c  <= 1000
 0  <= R_LMPD_251_tyr-L_c  <= 0
 0  <= F_LMPD_252_udcpdp_c  <= 1000
 0  <= R_LMPD_252_udcpdp_c  <= 0
 0  <= F_LMPD_253_utp_c  <= 1000
 0  <= R_LMPD_253_utp_c  <= 0
 0  <= F_LMPD_254_val-L_c  <= 1000
 0  <= R_LMPD_254_val-L_c  <= 0
 -1.381551e+01  <= LC_10fthf_c  <= -2.995732e+00
 -1.381551e+01  <= LC_13dpg_c  <= -2.995732e+00
 -1.381551e+01  <= LC_2ddg6p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_2dmmq8_c  <= -2.995732e+00
 -1.381551e+01  <= LC_2dmmql8_c  <= -2.995732e+00
 -1.381551e+01  <= LC_2h3oppan_c  <= -2.995732e+00
 -1.381551e+01  <= LC_2pg_c  <= -2.995732e+00
 -1.381551e+01  <= LC_3pg_c  <= -2.995732e+00
 -1.381551e+01  <= LC_4crsol_c  <= -2.995732e+00
 -1.381551e+01  <= LC_5drib_c  <= -2.995732e+00
 -1.381551e+01  <= LC_5mthf_c  <= -2.995732e+00
 -1.381551e+01  <= LC_6pgc_c  <= -2.995732e+00
 -1.381551e+01  <= LC_6pgl_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ac_c  <= -2.995732e+00
 -1.381551e+01  <= LC_acald_c  <= -2.995732e+00
 -1.381551e+01  <= LC_accoa_c  <= -2.995732e+00
 -1.381551e+01  <= LC_acon-C_c  <= -2.995732e+00
 -1.381551e+01  <= LC_actp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_adp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_adpglc_c  <= -2.995732e+00
 -1.381551e+01  <= LC_akg_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ala-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_amet_c  <= -2.995732e+00
 -1.381551e+01  <= LC_amob_c  <= -2.995732e+00
 -1.381551e+01  <= LC_amp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_arg-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_asn-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_asp-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_atp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_bglycogen_c  <= -2.995732e+00
 -1.381551e+01  <= LC_btn_c  <= -2.995732e+00
 -1.381551e+01  <= LC_chor_c  <= -2.995732e+00
 -1.381551e+01  <= LC_cit_c  <= -2.995732e+00
 -1.381551e+01  <= LC_co2_c  <= -2.995732e+00
 -1.381551e+01  <= LC_coa_c  <= -2.995732e+00
 -1.381551e+01  <= LC_colipa_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ctp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_cys-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_datp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_db4p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_dctp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_dgtp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_dha_c  <= -2.995732e+00
 -1.381551e+01  <= LC_dhap_c  <= -2.995732e+00
 -1.381551e+01  <= LC_dhor-S_c  <= -2.995732e+00
 -1.381551e+01  <= LC_dttp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_e4p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_etoh_c  <= -2.995732e+00
 -1.381551e+01  <= LC_f6p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_fad_c  <= -2.995732e+00
 -1.381551e+01  <= LC_fdp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_fe2_c  <= -2.995732e+00
 -1.381551e+01  <= LC_fe3_c  <= -2.995732e+00
 -1.381551e+01  <= LC_for_c  <= -2.995732e+00
 -1.381551e+01  <= LC_fum_c  <= -2.995732e+00
 -1.381551e+01  <= LC_g1p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_g3p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_g6p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glc-D_c  <= -2.995732e+00
 -1.381551e+01  <= LC_gln-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glu-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glx_c  <= -2.995732e+00
 -1.381551e+01  <= LC_gly_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glyc_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glyc-R_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glyc3p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glyclt_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glycogen_c  <= -2.995732e+00
 -1.381551e+01  <= LC_glycogenn1_c  <= -2.995732e+00
 -1.381551e+01  <= LC_gthrd_c  <= -2.995732e+00
 -1.381551e+01  <= LC_gtp_c  <= -2.995732e+00
 -1.726939e+01  <= LC_h_c  <= -1.726939e+01
 0  <= LC_h2o_c  <= 0
 -1.381551e+01  <= LC_h2o2_c  <= -2.995732e+00
 -1.381551e+01  <= LC_his-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_iasp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_icit_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ile-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_lac-D_c  <= -2.995732e+00
 -1.381551e+01  <= LC_leu-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_lys-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_mal-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_malcoa_c  <= -2.995732e+00
 -1.381551e+01  <= LC_met-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_mlthf_c  <= -2.995732e+00
 -1.381551e+01  <= LC_mql8_c  <= -2.995732e+00
 -1.381551e+01  <= LC_mqn8_c  <= -2.995732e+00
 -1.381551e+01  <= LC_mththf_c  <= -2.995732e+00
 -1.381551e+01  <= LC_nad_c  <= -2.995732e+00
 -1.381551e+01  <= LC_nadh_c  <= -2.995732e+00
 -1.381551e+01  <= LC_nadp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_nadph_c  <= -2.995732e+00
 -1.381551e+01  <= LC_nh4_c  <= -2.995732e+00
 -1.381551e+01  <= LC_o2_c  <= -2.995732e+00
 -1.381551e+01  <= LC_oaa_c  <= -2.995732e+00
 -1.381551e+01  <= LC_orot_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pe160_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pe161_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pe181_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pep_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pg160_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pg161_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pg181_c  <= -2.995732e+00
 -1.381551e+01  <= LC_phe-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pi_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ppi_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pro-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ptrc_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pydx5p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_pyr_c  <= -2.995732e+00
 -1.381551e+01  <= LC_q8_c  <= -2.995732e+00
 -1.381551e+01  <= LC_q8h2_c  <= -2.995732e+00
 -1.381551e+01  <= LC_r5p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ribflv_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ru5p-D_c  <= -2.995732e+00
 -1.381551e+01  <= LC_s17bp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_s7p_c  <= -2.995732e+00
 -1.381551e+01  <= LC_ser-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_so4_c  <= -2.995732e+00
 -1.381551e+01  <= LC_spmd_c  <= -2.995732e+00
 -1.381551e+01  <= LC_succ_c  <= -2.995732e+00
 -1.381551e+01  <= LC_succoa_c  <= -2.995732e+00
 -1.381551e+01  <= LC_thf_c  <= -2.995732e+00
 -1.381551e+01  <= LC_thmpp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_thr-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_trdox_c  <= -2.995732e+00
 -1.381551e+01  <= LC_trdrd_c  <= -2.995732e+00
 -1.381551e+01  <= LC_trp-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_tyr-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_udcpdp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_utp_c  <= -2.995732e+00
 -1.381551e+01  <= LC_val-L_c  <= -2.995732e+00
 -1.381551e+01  <= LC_xu5p-D_c  <= -2.995732e+00
 -1.842068e+01  <= LC_5mtr_e  <= -2.302585e+00
 -1.842068e+01  <= LC_ac_e  <= -2.302585e+00
 -1.842068e+01  <= LC_acald_e  <= -2.302585e+00
 -1.842068e+01  <= LC_acser_e  <= -2.302585e+00
 -1.842068e+01  <= LC_akg_e  <= -2.302585e+00
 -1.842068e+01  <= LC_ala-D_e  <= -2.302585e+00
 -1.842068e+01  <= LC_ala-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_alaala_e  <= -2.302585e+00
 -1.842068e+01  <= LC_arg-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_asn-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_asp-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_cit_e  <= -2.302585e+00
 -1.842068e+01  <= LC_co2_e  <= -2.302585e+00
 -1.842068e+01  <= LC_colipa_e  <= -2.302585e+00
 -1.842068e+01  <= LC_cys-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_dha_e  <= -2.302585e+00
 -1.842068e+01  <= LC_etoh_e  <= -2.302585e+00
 -1.842068e+01  <= LC_fe2_e  <= -2.302585e+00
 -1.842068e+01  <= LC_fe3_e  <= -2.302585e+00
 -1.842068e+01  <= LC_for_e  <= -2.302585e+00
 -1.842068e+01  <= LC_glc-D_e  <= -2.302585e+00
 -1.842068e+01  <= LC_glcn_e  <= -2.302585e+00
 -1.842068e+01  <= LC_glu-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_gly_e  <= -2.302585e+00
 -1.842068e+01  <= LC_glyc_e  <= -2.302585e+00
 -1.842068e+01  <= LC_glyc-R_e  <= -2.302585e+00
 -1.842068e+01  <= LC_glyclt_e  <= -2.302585e+00
 -1.842068e+01  <= LC_gthrd_e  <= -2.302585e+00
 -1.611810e+01  <= LC_h_e  <= -1.611810e+01
 0  <= LC_h2o_e  <= 0
 -1.842068e+01  <= LC_his-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_ile-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_lac-D_e  <= -2.302585e+00
 -1.842068e+01  <= LC_leu-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_lys-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_mal-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_meoh_e  <= -2.302585e+00
 -1.842068e+01  <= LC_nh4_e  <= -2.302585e+00
 -1.842068e+01  <= LC_o2_e  <= -2.302585e+00
 -1.842068e+01  <= LC_phe-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_pi_e  <= -2.302585e+00
 -1.842068e+01  <= LC_pro-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_pyr_e  <= -2.302585e+00
 -1.842068e+01  <= LC_ser-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_so4_e  <= -2.302585e+00
 -1.842068e+01  <= LC_spmd_e  <= -2.302585e+00
 -1.842068e+01  <= LC_succ_e  <= -2.302585e+00
 -1.842068e+01  <= LC_thr-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_trp-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_tyr-L_e  <= -2.302585e+00
 -1.842068e+01  <= LC_val-L_e  <= -2.302585e+00
 -1.381551e+01  <= LC_5mtr_p  <= -2.995732e+00
 -1.381551e+01  <= LC_ac_p  <= -2.995732e+00
 -1.381551e+01  <= LC_acald_p  <= -2.995732e+00
 -1.381551e+01  <= LC_acser_p  <= -2.995732e+00
 -1.381551e+01  <= LC_akg_p  <= -2.995732e+00
 -1.381551e+01  <= LC_ala-D_p  <= -2.995732e+00
 -1.381551e+01  <= LC_ala-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_alaala_p  <= -2.995732e+00
 -1.381551e+01  <= LC_arg-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_asn-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_asp-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_cit_p  <= -2.995732e+00
 -1.381551e+01  <= LC_clpn160_p  <= -2.995732e+00
 -1.381551e+01  <= LC_clpn161_p  <= -2.995732e+00
 -1.381551e+01  <= LC_clpn181_p  <= -2.995732e+00
 -1.381551e+01  <= LC_co2_p  <= -2.995732e+00
 -1.381551e+01  <= LC_colipa_p  <= -2.995732e+00
 -1.381551e+01  <= LC_cys-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_dha_p  <= -2.995732e+00
 -1.381551e+01  <= LC_etoh_p  <= -2.995732e+00
 -1.381551e+01  <= LC_fe2_p  <= -2.995732e+00
 -1.381551e+01  <= LC_fe3_p  <= -2.995732e+00
 -1.381551e+01  <= LC_for_p  <= -2.995732e+00
 -1.381551e+01  <= LC_glc-D_p  <= -2.995732e+00
 -1.381551e+01  <= LC_glcn_p  <= -2.995732e+00
 -1.381551e+01  <= LC_glu-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_gly_p  <= -2.995732e+00
 -1.381551e+01  <= LC_glyc_p  <= -2.995732e+00
 -1.381551e+01  <= LC_glyc-R_p  <= -2.995732e+00
 -1.381551e+01  <= LC_glyclt_p  <= -2.995732e+00
 -1.381551e+01  <= LC_gthrd_p  <= -2.995732e+00
 -1.611810e+01  <= LC_h_p  <= -1.611810e+01
 0  <= LC_h2o_p  <= 0
 -1.381551e+01  <= LC_his-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_ile-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_lac-D_p  <= -2.995732e+00
 -1.381551e+01  <= LC_leu-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_lys-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_mal-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_meoh_p  <= -2.995732e+00
 -1.381551e+01  <= LC_murein3p3p_p  <= -2.995732e+00
 -1.381551e+01  <= LC_murein3px4p_p  <= -2.995732e+00
 -1.381551e+01  <= LC_murein4p4p_p  <= -2.995732e+00
 -1.381551e+01  <= LC_murein4px4p_p  <= -2.995732e+00
 -1.381551e+01  <= LC_murein4px4px4p_p  <= -2.995732e+00
 -1.381551e+01  <= LC_nh4_p  <= -2.995732e+00
 -1.381551e+01  <= LC_o2_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pe160_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pe161_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pe181_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pg160_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pg161_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pg181_p  <= -2.995732e+00
 -1.381551e+01  <= LC_phe-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pi_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pro-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_pyr_p  <= -2.995732e+00
 -1.381551e+01  <= LC_ser-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_so4_p  <= -2.995732e+00
 -1.381551e+01  <= LC_spmd_p  <= -2.995732e+00
 -1.381551e+01  <= LC_succ_p  <= -2.995732e+00
 -1.381551e+01  <= LC_thr-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_trp-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_tyr-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_val-L_p  <= -2.995732e+00
 -1.381551e+01  <= LC_ade_p  <= -2.995732e+00
 0  <= FU_DM_4CRSOL  <= 1
 0  <= BU_DM_4CRSOL  <= 1
 0  <= FU_DM_5DRIB  <= 1
 0  <= BU_DM_5DRIB  <= 1
 0  <= FU_DM_AMOB  <= 1
 0  <= BU_DM_AMOB  <= 1
 0  <= FU_DM_MTHTHF  <= 1
 0  <= BU_DM_MTHTHF  <= 1
 0  <= FU_Ec_biomass_iJO1366_WT_53p95M  <= 1
 0  <= BU_Ec_biomass_iJO1366_WT_53p95M  <= 1
 0  <= FU_DM_5mtr_e  <= 1
 0  <= BU_DM_5mtr_e  <= 1
 0  <= FU_DM_ac_e  <= 1
 0  <= BU_DM_ac_e  <= 1
 0  <= FU_DM_acald_e  <= 1
 0  <= BU_DM_acald_e  <= 1
 0  <= FU_DM_acser_e  <= 1
 0  <= BU_DM_acser_e  <= 1
 0  <= FU_DM_akg_e  <= 1
 0  <= BU_DM_akg_e  <= 1
 0  <= FU_DM_ala-D_e  <= 1
 0  <= BU_DM_ala-D_e  <= 1
 0  <= FU_DM_ala-L_e  <= 1
 0  <= BU_DM_ala-L_e  <= 1
 0  <= FU_DM_alaala_e  <= 1
 0  <= BU_DM_alaala_e  <= 1
 0  <= FU_DM_arg-L_e  <= 1
 0  <= BU_DM_arg-L_e  <= 1
 0  <= FU_DM_asn-L_e  <= 1
 0  <= BU_DM_asn-L_e  <= 1
 0  <= FU_DM_asp-L_e  <= 1
 0  <= BU_DM_asp-L_e  <= 1
 0  <= FU_DM_ca2_e  <= 1
 0  <= BU_DM_ca2_e  <= 1
 0  <= FU_DM_cbl1_e  <= 1
 0  <= BU_DM_cbl1_e  <= 1
 0  <= FU_DM_cit_e  <= 1
 0  <= BU_DM_cit_e  <= 1
 0  <= FU_DM_cl_e  <= 1
 0  <= BU_DM_cl_e  <= 1
 0  <= FU_DM_co2_e  <= 1
 0  <= BU_DM_co2_e  <= 1
 0  <= FU_DM_cobalt2_e  <= 1
 0  <= BU_DM_cobalt2_e  <= 1
 0  <= FU_DM_colipa_e  <= 1
 0  <= BU_DM_colipa_e  <= 1
 0  <= FU_DM_cu2_e  <= 1
 0  <= BU_DM_cu2_e  <= 1
 0  <= FU_DM_cys-L_e  <= 1
 0  <= BU_DM_cys-L_e  <= 1
 0  <= FU_DM_dha_e  <= 1
 0  <= BU_DM_dha_e  <= 1
 0  <= FU_DM_enter_e  <= 1
 0  <= BU_DM_enter_e  <= 1
 0  <= FU_DM_etoh_e  <= 1
 0  <= BU_DM_etoh_e  <= 1
 0  <= FU_DM_fe2_e  <= 1
 0  <= BU_DM_fe2_e  <= 1
 0  <= FU_DM_fe3_e  <= 1
 0  <= BU_DM_fe3_e  <= 1
 0  <= FU_DM_for_e  <= 1
 0  <= BU_DM_for_e  <= 1
 0  <= FU_DM_glc_e  <= 1
 0  <= BU_DM_glc_e  <= 1
 0  <= FU_DM_glcn_e  <= 1
 0  <= BU_DM_glcn_e  <= 1
 0  <= FU_DM_glu-L_e  <= 1
 0  <= BU_DM_glu-L_e  <= 1
 0  <= FU_DM_gly_e  <= 1
 0  <= BU_DM_gly_e  <= 1
 0  <= FU_DM_glyc_e  <= 1
 0  <= BU_DM_glyc_e  <= 1
 0  <= FU_DM_glyc-R_e  <= 1
 0  <= BU_DM_glyc-R_e  <= 1
 0  <= FU_DM_glycogenn1_c  <= 1
 0  <= BU_DM_glycogenn1_c  <= 1
 0  <= FU_DM_glyclt_e  <= 1
 0  <= BU_DM_glyclt_e  <= 1
 0  <= FU_DM_gthrd_e  <= 1
 0  <= BU_DM_gthrd_e  <= 1
 0  <= FU_DM_h_e  <= 1
 0  <= BU_DM_h_e  <= 1
 0  <= FU_DM_h2o_e  <= 1
 0  <= BU_DM_h2o_e  <= 1
 0  <= FU_DM_his-L_e  <= 1
 0  <= BU_DM_his-L_e  <= 1
 0  <= FU_DM_ile-L_e  <= 1
 0  <= BU_DM_ile-L_e  <= 1
 0  <= FU_DM_k_e  <= 1
 0  <= BU_DM_k_e  <= 1
 0  <= FU_DM_lac-D_e  <= 1
 0  <= BU_DM_lac-D_e  <= 1
 0  <= FU_DM_leu-L_e  <= 1
 0  <= BU_DM_leu-L_e  <= 1
 0  <= FU_DM_lys-L_e  <= 1
 0  <= BU_DM_lys-L_e  <= 1
 0  <= FU_DM_mal-L_e  <= 1
 0  <= BU_DM_mal-L_e  <= 1
 0  <= FU_DM_meoh_e  <= 1
 0  <= BU_DM_meoh_e  <= 1
 0  <= FU_DM_mg2_e  <= 1
 0  <= BU_DM_mg2_e  <= 1
 0  <= FU_DM_mn2_e  <= 1
 0  <= BU_DM_mn2_e  <= 1
 0  <= FU_DM_mobd_e  <= 1
 0  <= BU_DM_mobd_e  <= 1
 0  <= FU_DM_nh4_e  <= 1
 0  <= BU_DM_nh4_e  <= 1
 0  <= FU_DM_ni2_e  <= 1
 0  <= BU_DM_ni2_e  <= 1
 0  <= FU_DM_o2_e  <= 1
 0  <= BU_DM_o2_e  <= 1
 0  <= FU_DM_phe-L_e  <= 1
 0  <= BU_DM_phe-L_e  <= 1
 0  <= FU_DM_pheme_e  <= 1
 0  <= BU_DM_pheme_e  <= 1
 0  <= FU_DM_pi_e  <= 1
 0  <= BU_DM_pi_e  <= 1
 0  <= FU_DM_pro-L_e  <= 1
 0  <= BU_DM_pro-L_e  <= 1
 0  <= FU_DM_pyr_e  <= 1
 0  <= BU_DM_pyr_e  <= 1
 0  <= FU_DM_ser-L_e  <= 1
 0  <= BU_DM_ser-L_e  <= 1
 0  <= FU_DM_so4_e  <= 1
 0  <= BU_DM_so4_e  <= 1
 0  <= FU_DM_spmd_e  <= 1
 0  <= BU_DM_spmd_e  <= 1
 0  <= FU_DM_succ_e  <= 1
 0  <= BU_DM_succ_e  <= 1
 0  <= FU_DM_thr-L_e  <= 1
 0  <= BU_DM_thr-L_e  <= 1
 0  <= FU_DM_trp-L_e  <= 1
 0  <= BU_DM_trp-L_e  <= 1
 0  <= FU_DM_tyr-L_e  <= 1
 0  <= BU_DM_tyr-L_e  <= 1
 0  <= FU_DM_val-L_e  <= 1
 0  <= BU_DM_val-L_e  <= 1
 0  <= FU_DM_zn2_e  <= 1
 0  <= BU_DM_zn2_e  <= 1
 -1000  <= DG_5MTRtex  <= 1000
 -2  <= DGo_5MTRtex  <= 2
 0  <= FU_5MTRtex  <= 1
 0  <= BU_5MTRtex  <= 1
 -1000  <= DG_ACALD  <= 1000
 -9.770762e+00  <= DGo_ACALD  <= -6.576323e+00
 0  <= FU_ACALD  <= 1
 0  <= BU_ACALD  <= 1
 -1000  <= DG_ACALDtex  <= 1000
 -2  <= DGo_ACALDtex  <= 2
 0  <= FU_ACALDtex  <= 1
 0  <= BU_ACALDtex  <= 1
 -1000  <= DG_ACALDtpp  <= 1000
 -1.225935e+00  <= DGo_ACALDtpp  <= 2.774065e+00
 0  <= FU_ACALDtpp  <= 1
 0  <= BU_ACALDtpp  <= 1
 -1000  <= DG_ACKr  <= 1000
 5.148750e+00  <= DGo_ACKr  <= 6.304408e+00
 0  <= FU_ACKr  <= 1
 0  <= BU_ACKr  <= 1
 -1000  <= DG_ACONTa  <= 1000
 5.804269e-01  <= DGo_ACONTa  <= 2.376633e+00
 0  <= FU_ACONTa  <= 1
 0  <= BU_ACONTa  <= 1
 -1000  <= DG_ACONTb  <= 1000
 -9.551789e-01  <= DGo_ACONTb  <= 8.378208e-01
 0  <= FU_ACONTb  <= 1
 0  <= BU_ACONTb  <= 1
 -1000  <= DG_ACS  <= 1000
 -5.051483e+00  <= DGo_ACS  <= -2.393218e+00
 0  <= FU_ACS  <= 1
 0  <= BU_ACS  <= 1
 -1000  <= DG_ACSERtex  <= 1000
 -2  <= DGo_ACSERtex  <= 2
 0  <= FU_ACSERtex  <= 1
 0  <= BU_ACSERtex  <= 1
 -1000  <= DG_ACt2rpp  <= 1000
 -2.291863e+00  <= DGo_ACt2rpp  <= 1.708137e+00
 0  <= FU_ACt2rpp  <= 1
 0  <= BU_ACt2rpp  <= 1
 -1000  <= DG_ACtex  <= 1000
 -2  <= DGo_ACtex  <= 2
 0  <= FU_ACtex  <= 1
 0  <= BU_ACtex  <= 1
 -1000  <= DG_ADK1  <= 1000
 -1.570456e+00  <= DGo_ADK1  <= 2.429544e+00
 0  <= FU_ADK1  <= 1
 0  <= BU_ADK1  <= 1
 -1000  <= DG_AKGDH  <= 1000
 -1.328682e+01  <= DGo_AKGDH  <= -1.009143e+01
 0  <= FU_AKGDH  <= 1
 0  <= BU_AKGDH  <= 1
 -1000  <= DG_AKGt2rpp  <= 1000
 7.778831e-01  <= DGo_AKGt2rpp  <= 4.777883e+00
 0  <= FU_AKGt2rpp  <= 1
 0  <= BU_AKGt2rpp  <= 1
 -1000  <= DG_AKGtex  <= 1000
 -2  <= DGo_AKGtex  <= 2
 0  <= FU_AKGtex  <= 1
 0  <= BU_AKGtex  <= 1
 -1000  <= DG_ALAALAtex  <= 1000
 -2  <= DGo_ALAALAtex  <= 2
 0  <= FU_ALAALAtex  <= 1
 0  <= BU_ALAALAtex  <= 1
 -1000  <= DG_ALAt2rpp  <= 1000
 -4.786204e+00  <= DGo_ALAt2rpp  <= -7.862039e-01
 0  <= FU_ALAt2rpp  <= 1
 0  <= BU_ALAt2rpp  <= 1
 -1000  <= DG_ALAtex  <= 1000
 -2  <= DGo_ALAtex  <= 2
 0  <= FU_ALAtex  <= 1
 0  <= BU_ALAtex  <= 1
 -1000  <= DG_ALCD2x  <= 1000
 3.048356e+00  <= DGo_ALCD2x  <= 5.573471e+00
 0  <= FU_ALCD2x  <= 1
 0  <= BU_ALCD2x  <= 1
 -1000  <= DG_ALDD2x  <= 1000
 -1.590265e+01  <= DGo_ALDD2x  <= -1.318673e+01
 0  <= FU_ALDD2x  <= 1
 0  <= BU_ALDD2x  <= 1
 -1000  <= DG_ALDD2y  <= 1000
 -1.628968e+01  <= DGo_ALDD2y  <= -1.357376e+01
 0  <= FU_ALDD2y  <= 1
 0  <= BU_ALDD2y  <= 1
 -1000  <= DG_ARGt3pp  <= 1000
 -5.390896e+00  <= DGo_ARGt3pp  <= -1.390896e+00
 0  <= FU_ARGt3pp  <= 1
 0  <= BU_ARGt3pp  <= 1
 -1000  <= DG_ARGtex  <= 1000
 -2  <= DGo_ARGtex  <= 2
 0  <= FU_ARGtex  <= 1
 0  <= BU_ARGtex  <= 1
 -1000  <= DG_ASNt2rpp  <= 1000
 -4.592688e+00  <= DGo_ASNt2rpp  <= -5.926876e-01
 0  <= FU_ASNt2rpp  <= 1
 0  <= BU_ASNt2rpp  <= 1
 -1000  <= DG_ASNtex  <= 1000
 -2  <= DGo_ASNtex  <= 2
 0  <= FU_ASNtex  <= 1
 0  <= BU_ASNtex  <= 1
 -1000  <= DG_ASPO3  <= 1000
 -3.964964e+01  <= DGo_ASPO3  <= -3.338357e+01
 0  <= FU_ASPO3  <= 1
 0  <= BU_ASPO3  <= 1
 -1000  <= DG_ASPO4  <= 1000
 -5.358119e+01  <= DGo_ASPO4  <= -4.688936e+01
 0  <= FU_ASPO4  <= 1
 0  <= BU_ASPO4  <= 1
 -1000  <= DG_ASPO5  <= 1000
 -3.591210e+01  <= DGo_ASPO5  <= -3.224512e+01
 0  <= FU_ASPO5  <= 1
 0  <= BU_ASPO5  <= 1
 -1000  <= DG_ASPO6  <= 1000
 -5.262102e+01  <= DGo_ASPO6  <= -4.895159e+01
 0  <= FU_ASPO6  <= 1
 0  <= BU_ASPO6  <= 1
 -1000  <= DG_ASPT  <= 1000
 5.453279e+00  <= DGo_ASPT  <= 7.224369e+00
 0  <= FU_ASPT  <= 1
 0  <= BU_ASPT  <= 1
 -1000  <= DG_ASPTA  <= 1000
 -3.144941e+00  <= DGo_ASPTA  <= 8.550587e-01
 0  <= FU_ASPTA  <= 1
 0  <= BU_ASPTA  <= 1
 -1000  <= DG_ASPt2rpp  <= 1000
 -3.162894e+00  <= DGo_ASPt2rpp  <= 8.371062e-01
 0  <= FU_ASPt2rpp  <= 1
 0  <= BU_ASPt2rpp  <= 1
 -1000  <= DG_ASPtex  <= 1000
 -2  <= DGo_ASPtex  <= 2
 0  <= FU_ASPtex  <= 1
 0  <= BU_ASPtex  <= 1
 -1000  <= DG_ATPM  <= 1000
 -8.480180e+00  <= DGo_ATPM  <= -7.045071e+00
 0  <= FU_ATPM  <= 1
 0  <= BU_ATPM  <= 1
 -1000  <= DG_ATPS4rpp  <= 1000
 -9.518201e+00  <= DGo_ATPS4rpp  <= -8.083092e+00
 0  <= FU_ATPS4rpp  <= 1
 0  <= BU_ATPS4rpp  <= 1
 0  <= FU_CA2tex  <= 1
 0  <= BU_CA2tex  <= 1
 0  <= FU_CAt6pp  <= 1
 0  <= BU_CAt6pp  <= 1
 0  <= FU_CBL1tonex  <= 1
 0  <= BU_CBL1tonex  <= 1
 -1000  <= DG_CITL  <= 1000
 1.745850e+00  <= DGo_CITL  <= 3.073403e+00
 0  <= FU_CITL  <= 1
 0  <= BU_CITL  <= 1
 -1000  <= DG_CITt3pp  <= 1000
 -1.586840e+01  <= DGo_CITt3pp  <= -1.186840e+01
 0  <= FU_CITt3pp  <= 1
 0  <= BU_CITt3pp  <= 1
 -1000  <= DG_CITt7pp  <= 1000
 7.861557e-01  <= DGo_CITt7pp  <= 4.786156e+00
 0  <= FU_CITt7pp  <= 1
 0  <= BU_CITt7pp  <= 1
 -1000  <= DG_CITtex  <= 1000
 -2  <= DGo_CITtex  <= 2
 0  <= FU_CITtex  <= 1
 0  <= BU_CITtex  <= 1
 0  <= FU_CLt3_2pp  <= 1
 0  <= BU_CLt3_2pp  <= 1
 0  <= FU_CLtex  <= 1
 0  <= BU_CLtex  <= 1
 -1000  <= DG_CO2tex  <= 1000
 -2  <= DGo_CO2tex  <= 2
 0  <= FU_CO2tex  <= 1
 0  <= BU_CO2tex  <= 1
 -1000  <= DG_CO2tpp  <= 1000
 -2.000000e+00  <= DGo_CO2tpp  <= 2.000000e+00
 0  <= FU_CO2tpp  <= 1
 0  <= BU_CO2tpp  <= 1
 0  <= FU_COBALT2tex  <= 1
 0  <= BU_COBALT2tex  <= 1
 0  <= FU_COBALT2tpp  <= 1
 0  <= BU_COBALT2tpp  <= 1
 -1000  <= DG_COLIPAabcpp  <= 1000
 -8.262144e+01  <= DGo_COLIPAabcpp  <= -8.118633e+01
 0  <= FU_COLIPAabcpp  <= 1
 0  <= BU_COLIPAabcpp  <= 1
 -1000  <= DG_COLIPAabctex  <= 1000
 -8.480180e+00  <= DGo_COLIPAabctex  <= -7.045071e+00
 0  <= FU_COLIPAabctex  <= 1
 0  <= BU_COLIPAabctex  <= 1
 -1000  <= DG_CS  <= 1000
 -1.008938e+01  <= DGo_CS  <= -7.472163e+00
 0  <= FU_CS  <= 1
 0  <= BU_CS  <= 1
 0  <= FU_CU2tex  <= 1
 0  <= BU_CU2tex  <= 1
 0  <= FU_CU2tpp  <= 1
 0  <= BU_CU2tpp  <= 1
 -1000  <= DG_CYSabc2pp  <= 1000
 -8.554899e+00  <= DGo_CYSabc2pp  <= -7.119790e+00
 0  <= FU_CYSabc2pp  <= 1
 0  <= BU_CYSabc2pp  <= 1
 -1000  <= DG_CYStex  <= 1000
 -2  <= DGo_CYStex  <= 2
 0  <= FU_CYStex  <= 1
 0  <= BU_CYStex  <= 1
 -1000  <= DG_CYStpp  <= 1000
 -2.074719e+00  <= DGo_CYStpp  <= 1.925281e+00
 0  <= FU_CYStpp  <= 1
 0  <= BU_CYStpp  <= 1
 -1000  <= DG_CYTBD2pp  <= 1000
 -1.788618e+01  <= DGo_CYTBD2pp  <= -1.200691e+01
 0  <= FU_CYTBD2pp  <= 1
 0  <= BU_CYTBD2pp  <= 1
 -1000  <= DG_CYTBDpp  <= 1000
 -3.136007e+01  <= DGo_CYTBDpp  <= -2.597036e+01
 0  <= FU_CYTBDpp  <= 1
 0  <= BU_CYTBDpp  <= 1
 -1000  <= DG_CYTBO3_4pp  <= 1000
 -2.307843e+01  <= DGo_CYTBO3_4pp  <= -1.768872e+01
 0  <= FU_CYTBO3_4pp  <= 1
 0  <= BU_CYTBO3_4pp  <= 1
 -1000  <= DG_D-LACt2pp  <= 1000
 -1.907305e+00  <= DGo_D-LACt2pp  <= 2.092695e+00
 0  <= FU_D-LACt2pp  <= 1
 0  <= BU_D-LACt2pp  <= 1
 -1000  <= DG_D-LACtex  <= 1000
 -2  <= DGo_D-LACtex  <= 2
 0  <= FU_D-LACtex  <= 1
 0  <= BU_D-LACtex  <= 1
 -1000  <= DG_DALAtex  <= 1000
 -2  <= DGo_DALAtex  <= 2
 0  <= FU_DALAtex  <= 1
 0  <= BU_DALAtex  <= 1
 -1000  <= DG_DB4PS  <= 1000
 -2.245931e+01  <= DGo_DB4PS  <= -2.126516e+01
 0  <= FU_DB4PS  <= 1
 0  <= BU_DB4PS  <= 1
 -1000  <= DG_DHAPT  <= 1000
 -1.072047e+01  <= DGo_DHAPT  <= -9.449523e+00
 0  <= FU_DHAPT  <= 1
 0  <= BU_DHAPT  <= 1
 -1000  <= DG_DHAtex  <= 1000
 -2  <= DGo_DHAtex  <= 2
 0  <= FU_DHAtex  <= 1
 0  <= BU_DHAtex  <= 1
 -1000  <= DG_DHAtpp  <= 1000
 -8.389021e-01  <= DGo_DHAtpp  <= 3.161098e+00
 0  <= FU_DHAtpp  <= 1
 0  <= BU_DHAtpp  <= 1
 -1000  <= DG_DHORD2  <= 1000
 -1.521840e+01  <= DGo_DHORD2  <= -1.051318e+01
 0  <= FU_DHORD2  <= 1
 0  <= BU_DHORD2  <= 1
 -1000  <= DG_DHORD5  <= 1000
 -2.931175e+01  <= DGo_DHORD5  <= -2.385717e+01
 0  <= FU_DHORD5  <= 1
 0  <= BU_DHORD5  <= 1
 -1000  <= DG_DHORDfum  <= 1000
 -1.119568e+01  <= DGo_DHORDfum  <= -9.659911e+00
 0  <= FU_DHORDfum  <= 1
 0  <= BU_DHORDfum  <= 1
 -1000  <= DG_EDA  <= 1000
 -3.471994e-01  <= DGo_EDA  <= 2.385338e-01
 0  <= FU_EDA  <= 1
 0  <= BU_EDA  <= 1
 -1000  <= DG_EDD  <= 1000
 -8.847150e+00  <= DGo_EDD  <= -7.460445e+00
 0  <= FU_EDD  <= 1
 0  <= BU_EDD  <= 1
 -1000  <= DG_ENO  <= 1000
 -2.439328e+00  <= DGo_ENO  <= -8.809013e-01
 0  <= FU_ENO  <= 1
 0  <= BU_ENO  <= 1
 -1000  <= DG_ETOHtex  <= 1000
 -2  <= DGo_ETOHtex  <= 2
 0  <= FU_ETOHtex  <= 1
 0  <= BU_ETOHtex  <= 1
 -1000  <= DG_ETOHtrpp  <= 1000
 -8.389021e-01  <= DGo_ETOHtrpp  <= 3.161098e+00
 0  <= FU_ETOHtrpp  <= 1
 0  <= BU_ETOHtrpp  <= 1
 -1000  <= DG_F6PA  <= 1000
 2.308389e+00  <= DGo_F6PA  <= 4.771069e+00
 0  <= FU_F6PA  <= 1
 0  <= BU_F6PA  <= 1
 -1000  <= DG_FBA  <= 1000
 2.855513e+00  <= DGo_FBA  <= 5.318193e+00
 0  <= FU_FBA  <= 1
 0  <= BU_FBA  <= 1
 -1000  <= DG_FBA3  <= 1000
 2.308001e+00  <= DGo_FBA3  <= 3.067255e+00
 0  <= FU_FBA3  <= 1
 0  <= BU_FBA3  <= 1
 -1000  <= DG_FBP  <= 1000
 -8.226627e-01  <= DGo_FBP  <= 6.486048e-01
 0  <= FU_FBP  <= 1
 0  <= BU_FBP  <= 1
 -1000  <= DG_FDH4pp  <= 1000
 -1.321296e+01  <= DGo_FDH4pp  <= -7.754118e+00
 0  <= FU_FDH4pp  <= 1
 0  <= BU_FDH4pp  <= 1
 -1000  <= DG_FDH5pp  <= 1000
 -2.717357e+01  <= DGo_FDH5pp  <= -2.123086e+01
 0  <= FU_FDH5pp  <= 1
 0  <= BU_FDH5pp  <= 1
 -1000  <= DG_FE2abcpp  <= 1000
 -1.617255e+01  <= DGo_FE2abcpp  <= -1.473744e+01
 0  <= FU_FE2abcpp  <= 1
 0  <= BU_FE2abcpp  <= 1
 -1000  <= DG_FE2t2pp  <= 1000
 -1.383318e+01  <= DGo_FE2t2pp  <= -9.833183e+00
 0  <= FU_FE2t2pp  <= 1
 0  <= BU_FE2t2pp  <= 1
 -1000  <= DG_FE2tex  <= 1000
 -2  <= DGo_FE2tex  <= 2
 0  <= FU_FE2tex  <= 1
 0  <= BU_FE2tex  <= 1
 -1000  <= DG_FE2tpp  <= 1000
 -9.692365e+00  <= DGo_FE2tpp  <= -5.692365e+00
 0  <= FU_FE2tpp  <= 1
 0  <= BU_FE2tpp  <= 1
 -1000  <= DG_FE3abcpp  <= 1000
 -2.059928e+01  <= DGo_FE3abcpp  <= -1.916417e+01
 0  <= FU_FE3abcpp  <= 1
 0  <= BU_FE3abcpp  <= 1
 -1000  <= DG_FE3tex  <= 1000
 -2  <= DGo_FE3tex  <= 2
 0  <= FU_FE3tex  <= 1
 0  <= BU_FE3tex  <= 1
 0  <= FU_FEENTERtex  <= 1
 0  <= BU_FEENTERtex  <= 1
 0  <= FU_FEENTERtpp  <= 1
 0  <= BU_FEENTERtpp  <= 1
 0  <= FU_FLDR2  <= 1
 0  <= BU_FLDR2  <= 1
 -1000  <= DG_FORtex  <= 1000
 -2  <= DGo_FORtex  <= 2
 0  <= FU_FORtex  <= 1
 0  <= BU_FORtex  <= 1
 -1000  <= DG_FORtppi  <= 1000
 -5.460070e+00  <= DGo_FORtppi  <= -1.460070e+00
 0  <= FU_FORtppi  <= 1
 0  <= BU_FORtppi  <= 1
 -1000  <= DG_FRD3  <= 1000
 2.674514e+00  <= DGo_FRD3  <= 7.939470e+00
 0  <= FU_FRD3  <= 1
 0  <= BU_FRD3  <= 1
 -1000  <= DG_FUM  <= 1000
 -1.486446e+00  <= DGo_FUM  <= 2.590055e-01
 0  <= FU_FUM  <= 1
 0  <= BU_FUM  <= 1
 -1000  <= DG_G3PD2  <= 1000
 -6.079920e-01  <= DGo_G3PD2  <= 1.551394e+00
 0  <= FU_G3PD2  <= 1
 0  <= BU_G3PD2  <= 1
 -1000  <= DG_G3PD5  <= 1000
 -1.774997e+01  <= DGo_G3PD5  <= -1.247439e+01
 0  <= FU_G3PD5  <= 1
 0  <= BU_G3PD5  <= 1
 -1000  <= DG_G3PD6  <= 1000
 -3.171826e+01  <= DGo_G3PD6  <= -2.594344e+01
 0  <= FU_G3PD6  <= 1
 0  <= BU_G3PD6  <= 1
 -1000  <= DG_G3PD7  <= 1000
 -2.069489e+01  <= DGo_G3PD7  <= -1.526747e+01
 0  <= FU_G3PD7  <= 1
 0  <= BU_G3PD7  <= 1
 -1000  <= DG_G3PT  <= 1000
 -2.801967e+00  <= DGo_G3PT  <= -1.330699e+00
 0  <= FU_G3PT  <= 1
 0  <= BU_G3PT  <= 1
 -1000  <= DG_G6PDH2r  <= 1000
 -5.004100e+00  <= DGo_G6PDH2r  <= -1.785680e+00
 0  <= FU_G6PDH2r  <= 1
 0  <= BU_G6PDH2r  <= 1
 -1000  <= DG_G6PP  <= 1000
 -2.340737e+00  <= DGo_G6PP  <= -8.694700e-01
 0  <= FU_G6PP  <= 1
 0  <= BU_G6PP  <= 1
 -1000  <= DG_GAPD  <= 1000
 -3.346819e+00  <= DGo_GAPD  <= -9.567827e-01
 0  <= FU_GAPD  <= 1
 0  <= BU_GAPD  <= 1
 -1000  <= DG_GLBRAN2  <= 1000
 -2  <= DGo_GLBRAN2  <= 2
 0  <= FU_GLBRAN2  <= 1
 0  <= BU_GLBRAN2  <= 1
 -1000  <= DG_GLCDpp  <= 1000
 -2.627912e+01  <= DGo_GLCDpp  <= -2.009313e+01
 0  <= FU_GLCDpp  <= 1
 0  <= BU_GLCDpp  <= 1
 -1000  <= DG_GLCNtex  <= 1000
 -2  <= DGo_GLCNtex  <= 2
 0  <= FU_GLCNtex  <= 1
 0  <= BU_GLCNtex  <= 1
 -1000  <= DG_GLCP  <= 1000
 -2.077938e+00  <= DGo_GLCP  <= -7.249845e-01
 0  <= FU_GLCP  <= 1
 0  <= BU_GLCP  <= 1
 -1000  <= DG_GLCS1  <= 1000
 -3.968046e+00  <= DGo_GLCS1  <= -2.703982e+00
 0  <= FU_GLCS1  <= 1
 0  <= BU_GLCS1  <= 1
 -1000  <= DG_GLCabcpp  <= 1000
 -6.157984e+00  <= DGo_GLCabcpp  <= -4.722876e+00
 0  <= FU_GLCabcpp  <= 1
 0  <= BU_GLCabcpp  <= 1
 -1000  <= DG_GLCptspp  <= 1000
 7.526943e-01  <= DGo_GLCptspp  <= 2.023640e+00
 0  <= FU_GLCptspp  <= 1
 0  <= BU_GLCptspp  <= 1
 -1000  <= DG_GLCt2pp  <= 1000
 -3.818622e+00  <= DGo_GLCt2pp  <= 1.813776e-01
 0  <= FU_GLCt2pp  <= 1
 0  <= BU_GLCt2pp  <= 1
 -1000  <= DG_GLCtex  <= 1000
 -2  <= DGo_GLCtex  <= 2
 0  <= FU_GLCtex  <= 1
 0  <= BU_GLCtex  <= 1
 -1000  <= DG_GLCtexi  <= 1000
 -2  <= DGo_GLCtexi  <= 2
 0  <= FU_GLCtexi  <= 1
 0  <= BU_GLCtexi  <= 1
 -1000  <= DG_GLDBRAN2  <= 1000
 -2  <= DGo_GLDBRAN2  <= 2
 0  <= FU_GLDBRAN2  <= 1
 0  <= BU_GLDBRAN2  <= 1
 -1000  <= DG_GLGC  <= 1000
 2.089255e+00  <= DGo_GLGC  <= 3.583044e+00
 0  <= FU_GLGC  <= 1
 0  <= BU_GLGC  <= 1
 -1000  <= DG_GLUDy  <= 1000
 9.922353e+00  <= DGo_GLUDy  <= 1.254667e+01
 0  <= FU_GLUDy  <= 1
 0  <= BU_GLUDy  <= 1
 -1000  <= DG_GLUabcpp  <= 1000
 -5.121142e+00  <= DGo_GLUabcpp  <= -3.686034e+00
 0  <= FU_GLUabcpp  <= 1
 0  <= BU_GLUabcpp  <= 1
 -1000  <= DG_GLUt2rpp  <= 1000
 -2.781780e+00  <= DGo_GLUt2rpp  <= 1.218220e+00
 0  <= FU_GLUt2rpp  <= 1
 0  <= BU_GLUt2rpp  <= 1
 0  <= FU_GLUt4pp  <= 1
 0  <= BU_GLUt4pp  <= 1
 -1000  <= DG_GLUtex  <= 1000
 -2  <= DGo_GLUtex  <= 2
 0  <= FU_GLUtex  <= 1
 0  <= BU_GLUtex  <= 1
 -1000  <= DG_GLXCL  <= 1000
 -1.016800e+01  <= DGo_GLXCL  <= -8.807914e+00
 0  <= FU_GLXCL  <= 1
 0  <= BU_GLXCL  <= 1
 -1000  <= DG_GLYCAt2rpp  <= 1000
 -1.907473e+00  <= DGo_GLYCAt2rpp  <= 2.092527e+00
 0  <= FU_GLYCAt2rpp  <= 1
 0  <= BU_GLYCAt2rpp  <= 1
 -1000  <= DG_GLYCAtex  <= 1000
 -2  <= DGo_GLYCAtex  <= 2
 0  <= FU_GLYCAtex  <= 1
 0  <= BU_GLYCAtex  <= 1
 -1000  <= DG_GLYCDx  <= 1000
 1.211220e+00  <= DGo_GLYCDx  <= 3.370606e+00
 0  <= FU_GLYCDx  <= 1
 0  <= BU_GLYCDx  <= 1
 -1000  <= DG_GLYCK  <= 1000
 -5.725058e+00  <= DGo_GLYCK  <= -5.251598e+00
 0  <= FU_GLYCK  <= 1
 0  <= BU_GLYCK  <= 1
 -1000  <= DG_GLYCK2  <= 1000
 -5.597066e+00  <= DGo_GLYCK2  <= -5.123606e+00
 0  <= FU_GLYCK2  <= 1
 0  <= BU_GLYCK2  <= 1
 -1000  <= DG_GLYCLTDx  <= 1000
 -7.296314e+00  <= DGo_GLYCLTDx  <= -5.125579e+00
 0  <= FU_GLYCLTDx  <= 1
 0  <= BU_GLYCLTDx  <= 1
 -1000  <= DG_GLYCLTDy  <= 1000
 -6.909281e+00  <= DGo_GLYCLTDy  <= -4.738547e+00
 0  <= FU_GLYCLTDy  <= 1
 0  <= BU_GLYCLTDy  <= 1
 -1000  <= DG_GLYCLTt2rpp  <= 1000
 -2.294468e+00  <= DGo_GLYCLTt2rpp  <= 1.705532e+00
 0  <= FU_GLYCLTt2rpp  <= 1
 0  <= BU_GLYCLTt2rpp  <= 1
 -1000  <= DG_GLYCLTtex  <= 1000
 -2  <= DGo_GLYCLTtex  <= 2
 0  <= FU_GLYCLTtex  <= 1
 0  <= BU_GLYCLTtex  <= 1
 -1000  <= DG_GLYCTO2  <= 1000
 -1.240008e+01  <= DGo_GLYCTO2  <= -7.119851e+00
 0  <= FU_GLYCTO2  <= 1
 0  <= BU_GLYCTO2  <= 1
 -1000  <= DG_GLYCTO3  <= 1000
 -2.636817e+01  <= DGo_GLYCTO3  <= -2.058910e+01
 0  <= FU_GLYCTO3  <= 1
 0  <= BU_GLYCTO3  <= 1
 -1000  <= DG_GLYCTO4  <= 1000
 -1.534494e+01  <= DGo_GLYCTO4  <= -9.912993e+00
 0  <= FU_GLYCTO4  <= 1
 0  <= BU_GLYCTO4  <= 1
 -1000  <= DG_GLYCtex  <= 1000
 -2  <= DGo_GLYCtex  <= 2
 0  <= FU_GLYCtex  <= 1
 0  <= BU_GLYCtex  <= 1
 -1000  <= DG_GLYCtpp  <= 1000
 -3.548130e+00  <= DGo_GLYCtpp  <= 4.518695e-01
 0  <= FU_GLYCtpp  <= 1
 0  <= BU_GLYCtpp  <= 1
 -1000  <= DG_GLYK  <= 1000
 -5.933022e+00  <= DGo_GLYK  <= -5.459563e+00
 0  <= FU_GLYK  <= 1
 0  <= BU_GLYK  <= 1
 -1000  <= DG_GLYt2rpp  <= 1000
 -5.173237e+00  <= DGo_GLYt2rpp  <= -1.173237e+00
 0  <= FU_GLYt2rpp  <= 1
 0  <= BU_GLYt2rpp  <= 1
 -1000  <= DG_GLYtex  <= 1000
 -2  <= DGo_GLYtex  <= 2
 0  <= FU_GLYtex  <= 1
 0  <= BU_GLYtex  <= 1
 -1000  <= DG_GND  <= 1000
 -6.902048e-01  <= DGo_GND  <= 1.764174e+00
 0  <= FU_GND  <= 1
 0  <= BU_GND  <= 1
 -1000  <= DG_GTHRDabc2pp  <= 1000
 -1.299411e+01  <= DGo_GTHRDabc2pp  <= -1.155900e+01
 0  <= FU_GTHRDabc2pp  <= 1
 0  <= BU_GTHRDabc2pp  <= 1
 -1000  <= DG_GTHRDtex  <= 1000
 -2  <= DGo_GTHRDtex  <= 2
 0  <= FU_GTHRDtex  <= 1
 0  <= BU_GTHRDtex  <= 1
 0  <= FU_H2Otex  <= 1
 0  <= BU_H2Otex  <= 1
 0  <= FU_H2Otpp  <= 1
 0  <= BU_H2Otpp  <= 1
 -1000  <= DG_HEX1  <= 1000
 -6.394252e+00  <= DGo_HEX1  <= -5.920792e+00
 0  <= FU_HEX1  <= 1
 0  <= BU_HEX1  <= 1
 -1000  <= DG_HISt2rpp  <= 1000
 -4.399171e+00  <= DGo_HISt2rpp  <= -3.991713e-01
 0  <= FU_HISt2rpp  <= 1
 0  <= BU_HISt2rpp  <= 1
 -1000  <= DG_HIStex  <= 1000
 -2  <= DGo_HIStex  <= 2
 0  <= FU_HIStex  <= 1
 0  <= BU_HIStex  <= 1
 -1000  <= DG_Htex  <= 1000
 -2  <= DGo_Htex  <= 2
 0  <= FU_Htex  <= 1
 0  <= BU_Htex  <= 1
 -1000  <= DG_ICDHyr  <= 1000
 1.895797e+00  <= DGo_ICDHyr  <= 4.377315e+00
 0  <= FU_ICDHyr  <= 1
 0  <= BU_ICDHyr  <= 1
 -1000  <= DG_ICL  <= 1000
 5.377016e+00  <= DGo_ICL  <= 6.219882e+00
 0  <= FU_ICL  <= 1
 0  <= BU_ICL  <= 1
 -1000  <= DG_ILEt2rpp  <= 1000
 -3.625106e+00  <= DGo_ILEt2rpp  <= 3.748939e-01
 0  <= FU_ILEt2rpp  <= 1
 0  <= BU_ILEt2rpp  <= 1
 -1000  <= DG_ILEtex  <= 1000
 -2  <= DGo_ILEtex  <= 2
 0  <= FU_ILEtex  <= 1
 0  <= BU_ILEtex  <= 1
 0  <= FU_Kabcpp  <= 1
 0  <= BU_Kabcpp  <= 1
 0  <= FU_Kt2pp  <= 1
 0  <= BU_Kt2pp  <= 1
 0  <= FU_Ktex  <= 1
 0  <= BU_Ktex  <= 1
 -1000  <= DG_LDH_D  <= 1000
 3.655880e+00  <= DGo_LDH_D  <= 5.846064e+00
 0  <= FU_LDH_D  <= 1
 0  <= BU_LDH_D  <= 1
 -1000  <= DG_LDH_D2  <= 1000
 -1.386407e+01  <= DGo_LDH_D2  <= -8.575812e+00
 0  <= FU_LDH_D2  <= 1
 0  <= BU_LDH_D2  <= 1
 -1000  <= DG_LEUt2rpp  <= 1000
 -3.625106e+00  <= DGo_LEUt2rpp  <= 3.748939e-01
 0  <= FU_LEUt2rpp  <= 1
 0  <= BU_LEUt2rpp  <= 1
 -1000  <= DG_LEUtex  <= 1000
 -2  <= DGo_LEUtex  <= 2
 0  <= FU_LEUtex  <= 1
 0  <= BU_LEUtex  <= 1
 -1000  <= DG_LYSt3pp  <= 1000
 -5.390896e+00  <= DGo_LYSt3pp  <= -1.390896e+00
 0  <= FU_LYSt3pp  <= 1
 0  <= BU_LYSt3pp  <= 1
 -1000  <= DG_LYStex  <= 1000
 -2  <= DGo_LYStex  <= 2
 0  <= FU_LYStex  <= 1
 0  <= BU_LYStex  <= 1
 -1000  <= DG_MALS  <= 1000
 -9.684243e+00  <= DGo_MALS  <= -7.172862e+00
 0  <= FU_MALS  <= 1
 0  <= BU_MALS  <= 1
 -1000  <= DG_MALt3pp  <= 1000
 -1.306888e+01  <= DGo_MALt3pp  <= -9.068875e+00
 0  <= FU_MALt3pp  <= 1
 0  <= BU_MALt3pp  <= 1
 -1000  <= DG_MALtex  <= 1000
 -2  <= DGo_MALtex  <= 2
 0  <= FU_MALtex  <= 1
 0  <= BU_MALtex  <= 1
 -1000  <= DG_MDH  <= 1000
 3.656409e+00  <= DGo_MDH  <= 5.846592e+00
 0  <= FU_MDH  <= 1
 0  <= BU_MDH  <= 1
 -1000  <= DG_MDH2  <= 1000
 -1.386354e+01  <= DGo_MDH2  <= -8.575284e+00
 0  <= FU_MDH2  <= 1
 0  <= BU_MDH2  <= 1
 -1000  <= DG_MDH3  <= 1000
 -2.783129e+01  <= DGo_MDH3  <= -2.204488e+01
 0  <= FU_MDH3  <= 1
 0  <= BU_MDH3  <= 1
 -1000  <= DG_ME1  <= 1000
 1.039022e-01  <= DGo_ME1  <= 2.537941e+00
 0  <= FU_ME1  <= 1
 0  <= BU_ME1  <= 1
 -1000  <= DG_ME2  <= 1000
 -2.831304e-01  <= DGo_ME2  <= 2.150909e+00
 0  <= FU_ME2  <= 1
 0  <= BU_ME2  <= 1
 -1000  <= DG_MEOHtex  <= 1000
 -2  <= DGo_MEOHtex  <= 2
 0  <= FU_MEOHtex  <= 1
 0  <= BU_MEOHtex  <= 1
 0  <= FU_MG2t3_2pp  <= 1
 0  <= BU_MG2t3_2pp  <= 1
 0  <= FU_MG2tex  <= 1
 0  <= BU_MG2tex  <= 1
 0  <= FU_MG2tpp  <= 1
 0  <= BU_MG2tpp  <= 1
 0  <= FU_MG2uabcpp  <= 1
 0  <= BU_MG2uabcpp  <= 1
 0  <= FU_MN2tpp  <= 1
 0  <= BU_MN2tpp  <= 1
 0  <= FU_MNt2pp  <= 1
 0  <= BU_MNt2pp  <= 1
 0  <= FU_MNtex  <= 1
 0  <= BU_MNtex  <= 1
 0  <= FU_MOBDabcpp  <= 1
 0  <= BU_MOBDabcpp  <= 1
 0  <= FU_MOBDtex  <= 1
 0  <= BU_MOBDtex  <= 1
 -1000  <= DG_MOX  <= 1000
 -2.626704e+01  <= DGo_MOX  <= -2.471118e+01
 0  <= FU_MOX  <= 1
 0  <= BU_MOX  <= 1
 -1000  <= DG_NADH10  <= 1000
 -3.275927e+01  <= DGo_NADH10  <= -2.661989e+01
 0  <= FU_NADH10  <= 1
 0  <= BU_NADH10  <= 1
 -1000  <= DG_NADH16pp  <= 1000
 -2.243801e+00  <= DGo_NADH16pp  <= 3.428518e+00
 0  <= FU_NADH16pp  <= 1
 0  <= BU_NADH16pp  <= 1
 -1000  <= DG_NADH17pp  <= 1000
 -1.619600e+01  <= DGo_NADH17pp  <= -1.005662e+01
 0  <= FU_NADH17pp  <= 1
 0  <= BU_NADH17pp  <= 1
 -1000  <= DG_NADH18pp  <= 1000
 -5.153863e+00  <= DGo_NADH18pp  <= 6.005804e-01
 0  <= FU_NADH18pp  <= 1
 0  <= BU_NADH18pp  <= 1
 -1000  <= DG_NADH5  <= 1000
 -1.880707e+01  <= DGo_NADH5  <= -1.313475e+01
 0  <= FU_NADH5  <= 1
 0  <= BU_NADH5  <= 1
 -1000  <= DG_NADH9  <= 1000
 -2.171714e+01  <= DGo_NADH9  <= -1.596269e+01
 0  <= FU_NADH9  <= 1
 0  <= BU_NADH9  <= 1
 -1000  <= DG_NADK  <= 1000
 -4.521629e+00  <= DGo_NADK  <= -4.048169e+00
 0  <= FU_NADK  <= 1
 0  <= BU_NADK  <= 1
 -1000  <= DG_NADPHQR2  <= 1000
 -1.842004e+01  <= DGo_NADPHQR2  <= -1.274772e+01
 0  <= FU_NADPHQR2  <= 1
 0  <= BU_NADPHQR2  <= 1
 -1000  <= DG_NADPHQR3  <= 1000
 -3.237224e+01  <= DGo_NADPHQR3  <= -2.623286e+01
 0  <= FU_NADPHQR3  <= 1
 0  <= BU_NADPHQR3  <= 1
 -1000  <= DG_NADPHQR4  <= 1000
 -2.133010e+01  <= DGo_NADPHQR4  <= -1.557566e+01
 0  <= FU_NADPHQR4  <= 1
 0  <= BU_NADPHQR4  <= 1
 -1000  <= DG_NADPPPS  <= 1000
 -4.213360e+00  <= DGo_NADPPPS  <= -2.742093e+00
 0  <= FU_NADPPPS  <= 1
 0  <= BU_NADPPPS  <= 1
 -1000  <= DG_NADTRHD  <= 1000
 -1.612967e+00  <= DGo_NADTRHD  <= 2.387033e+00
 0  <= FU_NADTRHD  <= 1
 0  <= BU_NADTRHD  <= 1
 0  <= FU_NAt3_1p5pp  <= 1
 0  <= BU_NAt3_1p5pp  <= 1
 0  <= FU_NAt3_2pp  <= 1
 0  <= BU_NAt3_2pp  <= 1
 0  <= FU_NAt3pp  <= 1
 0  <= BU_NAt3pp  <= 1
 -1000  <= DG_NH4tex  <= 1000
 -2  <= DGo_NH4tex  <= 2
 0  <= FU_NH4tex  <= 1
 0  <= BU_NH4tex  <= 1
 -1000  <= DG_NH4tpp  <= 1000
 -4.878601e+00  <= DGo_NH4tpp  <= -8.786010e-01
 0  <= FU_NH4tpp  <= 1
 0  <= BU_NH4tpp  <= 1
 0  <= FU_NI2tex  <= 1
 0  <= BU_NI2tex  <= 1
 0  <= FU_NI2tpp  <= 1
 0  <= BU_NI2tpp  <= 1
 0  <= FU_NI2uabcpp  <= 1
 0  <= BU_NI2uabcpp  <= 1
 -1000  <= DG_NTP1  <= 1000
 -8.480180e+00  <= DGo_NTP1  <= -7.045071e+00
 0  <= FU_NTP1  <= 1
 0  <= BU_NTP1  <= 1
 -1000  <= DG_NTPP6  <= 1000
 -1.084152e+01  <= DGo_NTPP6  <= -9.345475e+00
 0  <= FU_NTPP6  <= 1
 0  <= BU_NTPP6  <= 1
 -1000  <= DG_O2tex  <= 1000
 -2  <= DGo_O2tex  <= 2
 0  <= FU_O2tex  <= 1
 0  <= BU_O2tex  <= 1
 -1000  <= DG_O2tpp  <= 1000
 -2.000000e+00  <= DGo_O2tpp  <= 2.000000e+00
 0  <= FU_O2tpp  <= 1
 0  <= BU_O2tpp  <= 1
 -1000  <= DG_OAADC  <= 1000
 -3.961531e+00  <= DGo_OAADC  <= -2.899626e+00
 0  <= FU_OAADC  <= 1
 0  <= BU_OAADC  <= 1
 -1000  <= DG_PDH  <= 1000
 -1.212919e+01  <= DGo_PDH  <= -8.933791e+00
 0  <= FU_PDH  <= 1
 0  <= BU_PDH  <= 1
 -1000  <= DG_PE160abcpp  <= 1000
 -2.280039e+01  <= DGo_PE160abcpp  <= -2.136528e+01
 0  <= FU_PE160abcpp  <= 1
 0  <= BU_PE160abcpp  <= 1
 -1000  <= DG_PE161abcpp  <= 1000
 -2.202632e+01  <= DGo_PE161abcpp  <= -2.059121e+01
 0  <= FU_PE161abcpp  <= 1
 0  <= BU_PE161abcpp  <= 1
 -1000  <= DG_PE181abcpp  <= 1000
 -2.357445e+01  <= DGo_PE181abcpp  <= -2.213934e+01
 0  <= FU_PE181abcpp  <= 1
 0  <= BU_PE181abcpp  <= 1
 -1000  <= DG_PFK  <= 1000
 -7.912327e+00  <= DGo_PFK  <= -7.438867e+00
 0  <= FU_PFK  <= 1
 0  <= BU_PFK  <= 1
 -1000  <= DG_PFK_3  <= 1000
 -7.492830e+00  <= DGo_PFK_3  <= -7.019370e+00
 0  <= FU_PFK_3  <= 1
 0  <= BU_PFK_3  <= 1
 -1000  <= DG_PFL  <= 1000
 -8.262357e+00  <= DGo_PFL  <= -5.846918e+00
 0  <= FU_PFL  <= 1
 0  <= BU_PFL  <= 1
 -1000  <= DG_PG160abcpp  <= 1000
 -2.606602e+01  <= DGo_PG160abcpp  <= -2.463091e+01
 0  <= FU_PG160abcpp  <= 1
 0  <= BU_PG160abcpp  <= 1
 -1000  <= DG_PG161abcpp  <= 1000
 -2.529195e+01  <= DGo_PG161abcpp  <= -2.385685e+01
 0  <= FU_PG161abcpp  <= 1
 0  <= BU_PG161abcpp  <= 1
 -1000  <= DG_PG181abcpp  <= 1000
 -2.380039e+01  <= DGo_PG181abcpp  <= -2.236528e+01
 0  <= FU_PG181abcpp  <= 1
 0  <= BU_PG181abcpp  <= 1
 -1000  <= DG_PGI  <= 1000
 -1.699558e+00  <= DGo_PGI  <= -4.818010e-01
 0  <= FU_PGI  <= 1
 0  <= BU_PGI  <= 1
 -1000  <= DG_PGK  <= 1000
 5.092056e-01  <= DGo_PGK  <= 1.089984e+00
 0  <= FU_PGK  <= 1
 0  <= BU_PGK  <= 1
 -1000  <= DG_PGL  <= 1000
 -8.626880e+00  <= DGo_PGL  <= -5.785327e+00
 0  <= FU_PGL  <= 1
 0  <= BU_PGL  <= 1
 -1000  <= DG_PGM  <= 1000
 -2.127992e+00  <= DGo_PGM  <= 1.872008e+00
 0  <= FU_PGM  <= 1
 0  <= BU_PGM  <= 1
 -1000  <= DG_PGMT  <= 1000
 -2.106435e+00  <= DGo_PGMT  <= 1.893565e+00
 0  <= FU_PGMT  <= 1
 0  <= BU_PGMT  <= 1
 0  <= FU_PHEMEabcpp  <= 1
 0  <= BU_PHEMEabcpp  <= 1
 0  <= FU_PHEMEtiex  <= 1
 0  <= BU_PHEMEtiex  <= 1
 -1000  <= DG_PHEt2rpp  <= 1000
 -4.012139e+00  <= DGo_PHEt2rpp  <= -1.213869e-02
 0  <= FU_PHEt2rpp  <= 1
 0  <= BU_PHEt2rpp  <= 1
 -1000  <= DG_PHEtex  <= 1000
 -2  <= DGo_PHEtex  <= 2
 0  <= FU_PHEtex  <= 1
 0  <= BU_PHEtex  <= 1
 -1000  <= DG_PIt2rpp  <= 1000
 6.974249e-01  <= DGo_PIt2rpp  <= 4.697425e+00
 0  <= FU_PIt2rpp  <= 1
 0  <= BU_PIt2rpp  <= 1
 -1000  <= DG_PItex  <= 1000
 -2  <= DGo_PItex  <= 2
 0  <= FU_PItex  <= 1
 0  <= BU_PItex  <= 1
 -1000  <= DG_PIuabcpp  <= 1000
 4.594695e+00  <= DGo_PIuabcpp  <= 6.029803e+00
 0  <= FU_PIuabcpp  <= 1
 0  <= BU_PIuabcpp  <= 1
 0  <= FU_POR5  <= 1
 0  <= BU_POR5  <= 1
 -1000  <= DG_POX  <= 1000
 -3.566458e+01  <= DGo_POX  <= -3.008252e+01
 0  <= FU_POX  <= 1
 0  <= BU_POX  <= 1
 -1000  <= DG_PPA  <= 1000
 -7.086045e+00  <= DGo_PPA  <= -4.636555e+00
 0  <= FU_PPA  <= 1
 0  <= BU_PPA  <= 1
 0  <= FU_PPA2  <= 1
 0  <= BU_PPA2  <= 1
 -1000  <= DG_PPC  <= 1000
 -8.370873e+00  <= DGo_PPC  <= -6.206267e+00
 0  <= FU_PPC  <= 1
 0  <= BU_PPC  <= 1
 -1000  <= DG_PPCK  <= 1000
 -1.302466e+00  <= DGo_PPCK  <= 3.543552e-01
 0  <= FU_PPCK  <= 1
 0  <= BU_PPCK  <= 1
 0  <= FU_PPK2r  <= 1
 0  <= BU_PPK2r  <= 1
 -1000  <= DG_PPKr  <= 1000
 -2.618880e+00  <= DGo_PPKr  <= -1.183772e+00
 0  <= FU_PPKr  <= 1
 0  <= BU_PPKr  <= 1
 0  <= FU_PPPGO  <= 1
 0  <= BU_PPPGO  <= 1
 0  <= FU_PPPGO3  <= 1
 0  <= BU_PPPGO3  <= 1
 -1000  <= DG_PPS  <= 1000
 -6.221936e+00  <= DGo_PPS  <= -4.249357e+00
 0  <= FU_PPS  <= 1
 0  <= BU_PPS  <= 1
 -1000  <= DG_PROt2rpp  <= 1000
 -1.327054e+00  <= DGo_PROt2rpp  <= 2.672946e+00
 0  <= FU_PROt2rpp  <= 1
 0  <= BU_PROt2rpp  <= 1
 -1000  <= DG_PROtex  <= 1000
 -2  <= DGo_PROtex  <= 2
 0  <= FU_PROtex  <= 1
 0  <= BU_PROtex  <= 1
 -1000  <= DG_PTAr  <= 1000
 5.892377e+00  <= DGo_PTAr  <= 8.343742e+00
 0  <= FU_PTAr  <= 1
 0  <= BU_PTAr  <= 1
 -1000  <= DG_PYK  <= 1000
 -3.610828e+00  <= DGo_PYK  <= -2.302218e+00
 0  <= FU_PYK  <= 1
 0  <= BU_PYK  <= 1
 -1000  <= DG_PYRt2rpp  <= 1000
 -2.294635e+00  <= DGo_PYRt2rpp  <= 1.705365e+00
 0  <= FU_PYRt2rpp  <= 1
 0  <= BU_PYRt2rpp  <= 1
 -1000  <= DG_PYRtex  <= 1000
 -2  <= DGo_PYRtex  <= 2
 0  <= FU_PYRtex  <= 1
 0  <= BU_PYRtex  <= 1
 -1000  <= DG_RPE  <= 1000
 -2  <= DGo_RPE  <= 2
 0  <= FU_RPE  <= 1
 0  <= BU_RPE  <= 1
 -1000  <= DG_RPI  <= 1000
 -5.772118e-01  <= DGo_RPI  <= 2.025744e+00
 0  <= FU_RPI  <= 1
 0  <= BU_RPI  <= 1
 -1000  <= DG_SERt2rpp  <= 1000
 -4.786204e+00  <= DGo_SERt2rpp  <= -7.862039e-01
 0  <= FU_SERt2rpp  <= 1
 0  <= BU_SERt2rpp  <= 1
 -1000  <= DG_SERtex  <= 1000
 -2  <= DGo_SERtex  <= 2
 0  <= FU_SERtex  <= 1
 0  <= BU_SERtex  <= 1
 -1000  <= DG_SO4t2pp  <= 1000
 3.416732e-03  <= DGo_SO4t2pp  <= 4.003417e+00
 0  <= FU_SO4t2pp  <= 1
 0  <= BU_SO4t2pp  <= 1
 -1000  <= DG_SO4tex  <= 1000
 -2  <= DGo_SO4tex  <= 2
 0  <= FU_SO4tex  <= 1
 0  <= BU_SO4tex  <= 1
 -1000  <= DG_SPMDt3pp  <= 1000
 1.720920e+00  <= DGo_SPMDt3pp  <= 5.720920e+00
 0  <= FU_SPMDt3pp  <= 1
 0  <= BU_SPMDt3pp  <= 1
 -1000  <= DG_SPMDtex  <= 1000
 -2  <= DGo_SPMDtex  <= 2
 0  <= FU_SPMDtex  <= 1
 0  <= BU_SPMDtex  <= 1
 -1000  <= DG_SUCCt2_2pp  <= 1000
 -3.340208e+00  <= DGo_SUCCt2_2pp  <= 6.597919e-01
 0  <= FU_SUCCt2_2pp  <= 1
 0  <= BU_SUCCt2_2pp  <= 1
 -1000  <= DG_SUCCt2_3pp  <= 1000
 -7.481026e+00  <= DGo_SUCCt2_3pp  <= -3.481026e+00
 0  <= FU_SUCCt2_3pp  <= 1
 0  <= BU_SUCCt2_3pp  <= 1
 -1000  <= DG_SUCCt3pp  <= 1000
 -1.308225e+01  <= DGo_SUCCt3pp  <= -9.082246e+00
 0  <= FU_SUCCt3pp  <= 1
 0  <= BU_SUCCt3pp  <= 1
 -1000  <= DG_SUCCtex  <= 1000
 -2  <= DGo_SUCCtex  <= 2
 0  <= FU_SUCCtex  <= 1
 0  <= BU_SUCCtex  <= 1
 -1000  <= DG_SUCDi  <= 1000
 -4.992134e+00  <= DGo_SUCDi  <= 1.161493e-01
 0  <= FU_SUCDi  <= 1
 0  <= BU_SUCDi  <= 1
 -1000  <= DG_SUCOAS  <= 1000
 -5.300910e+00  <= DGo_SUCOAS  <= -2.874075e+00
 0  <= FU_SUCOAS  <= 1
 0  <= BU_SUCOAS  <= 1
 -1000  <= DG_SULabcpp  <= 1000
 -2.335945e+00  <= DGo_SULabcpp  <= -9.008367e-01
 0  <= FU_SULabcpp  <= 1
 0  <= BU_SULabcpp  <= 1
 -1000  <= DG_TALA  <= 1000
 -2.265415e+00  <= DGo_TALA  <= 3.059573e-01
 0  <= FU_TALA  <= 1
 0  <= BU_TALA  <= 1
 -1000  <= DG_THD2pp  <= 1000
 -1.066867e+01  <= DGo_THD2pp  <= -6.668669e+00
 0  <= FU_THD2pp  <= 1
 0  <= BU_THD2pp  <= 1
 -1000  <= DG_THIORDXi  <= 1000
 -7.965497e+01  <= DGo_THIORDXi  <= -7.472232e+01
 0  <= FU_THIORDXi  <= 1
 0  <= BU_THIORDXi  <= 1
 -1000  <= DG_THRt2pp  <= 1000
 -7.882465e+00  <= DGo_THRt2pp  <= -3.882465e+00
 0  <= FU_THRt2pp  <= 1
 0  <= BU_THRt2pp  <= 1
 -1000  <= DG_THRt2rpp  <= 1000
 -4.399171e+00  <= DGo_THRt2rpp  <= -3.991713e-01
 0  <= FU_THRt2rpp  <= 1
 0  <= BU_THRt2rpp  <= 1
 -1000  <= DG_THRtex  <= 1000
 -2  <= DGo_THRtex  <= 2
 0  <= FU_THRtex  <= 1
 0  <= BU_THRtex  <= 1
 -1000  <= DG_TKT1  <= 1000
 -5.886886e-01  <= DGo_TKT1  <= 2.092861e+00
 0  <= FU_TKT1  <= 1
 0  <= BU_TKT1  <= 1
 -1000  <= DG_TKT2  <= 1000
 -5.129775e+00  <= DGo_TKT2  <= -2.558403e+00
 0  <= FU_TKT2  <= 1
 0  <= BU_TKT2  <= 1
 -1000  <= DG_TPI  <= 1000
 1.159311e+00  <= DGo_TPI  <= 1.760689e+00
 0  <= FU_TPI  <= 1
 0  <= BU_TPI  <= 1
 -1000  <= DG_TRDR  <= 1000
 -5.033569e-01  <= DGo_TRDR  <= 4.464886e+00
 0  <= FU_TRDR  <= 1
 0  <= BU_TRDR  <= 1
 -1000  <= DG_TRPt2rpp  <= 1000
 -3.818622e+00  <= DGo_TRPt2rpp  <= 1.813776e-01
 0  <= FU_TRPt2rpp  <= 1
 0  <= BU_TRPt2rpp  <= 1
 -1000  <= DG_TRPtex  <= 1000
 -2  <= DGo_TRPtex  <= 2
 0  <= FU_TRPtex  <= 1
 0  <= BU_TRPtex  <= 1
 -1000  <= DG_TRSARr  <= 1000
 -4.820757e+00  <= DGo_TRSARr  <= -2.681100e+00
 0  <= FU_TRSARr  <= 1
 0  <= BU_TRSARr  <= 1
 -1000  <= DG_TYRt2rpp  <= 1000
 -4.012139e+00  <= DGo_TYRt2rpp  <= -1.213869e-02
 0  <= FU_TYRt2rpp  <= 1
 0  <= BU_TYRt2rpp  <= 1
 -1000  <= DG_TYRtex  <= 1000
 -2  <= DGo_TYRtex  <= 2
 0  <= FU_TYRtex  <= 1
 0  <= BU_TYRtex  <= 1
 -1000  <= DG_VALt2rpp  <= 1000
 -4.012139e+00  <= DGo_VALt2rpp  <= -1.213869e-02
 0  <= FU_VALt2rpp  <= 1
 0  <= BU_VALt2rpp  <= 1
 -1000  <= DG_VALtex  <= 1000
 -2  <= DGo_VALtex  <= 2
 0  <= FU_VALtex  <= 1
 0  <= BU_VALtex  <= 1
 0  <= FU_ZN2tpp  <= 1
 0  <= BU_ZN2tpp  <= 1
 0  <= FU_ZNabcpp  <= 1
 0  <= BU_ZNabcpp  <= 1
 0  <= FU_Zn2tex  <= 1
 0  <= BU_Zn2tex  <= 1
 -1000  <= DG_LMPD_1_10fthf_c  <= 1000
 -2.177097e+02  <= DGo_LMPD_1_10fthf_c  <= -1.928452e+02
 0  <= FU_LMPD_1_10fthf_c  <= 1
 0  <= BU_LMPD_1_10fthf_c  <= 1
 -1000  <= DG_LMPD_2_10fthf_c  <= 1000
 -2.172801e+02  <= DGo_LMPD_2_10fthf_c  <= -1.924157e+02
 0  <= FU_LMPD_2_10fthf_c  <= 1
 0  <= BU_LMPD_2_10fthf_c  <= 1
 -1000  <= DG_LMPD_3_10fthf_c  <= 1000
 -2.113119e+02  <= DGo_LMPD_3_10fthf_c  <= -1.842258e+02
 0  <= FU_LMPD_3_10fthf_c  <= 1
 0  <= BU_LMPD_3_10fthf_c  <= 1
 -1000  <= DG_LMPD_4_10fthf_c  <= 1000
 -2.108824e+02  <= DGo_LMPD_4_10fthf_c  <= -1.837962e+02
 0  <= FU_LMPD_4_10fthf_c  <= 1
 0  <= BU_LMPD_4_10fthf_c  <= 1
 0  <= FU_LMPD_5_2dmmql8_c  <= 1
 0  <= BU_LMPD_5_2dmmql8_c  <= 1
 0  <= FU_LMPD_6_2fe2s_c  <= 1
 0  <= BU_LMPD_6_2fe2s_c  <= 1
 0  <= FU_LMPD_7_2fe2s_c  <= 1
 0  <= BU_LMPD_7_2fe2s_c  <= 1
 0  <= FU_LMPD_8_2fe2s_c  <= 1
 0  <= BU_LMPD_8_2fe2s_c  <= 1
 0  <= FU_LMPD_9_2fe2s_c  <= 1
 0  <= BU_LMPD_9_2fe2s_c  <= 1
 0  <= FU_LMPD_10_2fe2s_c  <= 1
 0  <= BU_LMPD_10_2fe2s_c  <= 1
 0  <= FU_LMPD_11_2fe2s_c  <= 1
 0  <= BU_LMPD_11_2fe2s_c  <= 1
 0  <= FU_LMPD_12_2fe2s_c  <= 1
 0  <= BU_LMPD_12_2fe2s_c  <= 1
 0  <= FU_LMPD_13_2fe2s_c  <= 1
 0  <= BU_LMPD_13_2fe2s_c  <= 1
 0  <= FU_LMPD_14_4fe4s_c  <= 1
 0  <= BU_LMPD_14_4fe4s_c  <= 1
 0  <= FU_LMPD_15_4fe4s_c  <= 1
 0  <= BU_LMPD_15_4fe4s_c  <= 1
 0  <= FU_LMPD_16_4fe4s_c  <= 1
 0  <= BU_LMPD_16_4fe4s_c  <= 1
 0  <= FU_LMPD_17_4fe4s_c  <= 1
 0  <= BU_LMPD_17_4fe4s_c  <= 1
 0  <= FU_LMPD_18_4fe4s_c  <= 1
 0  <= BU_LMPD_18_4fe4s_c  <= 1
 0  <= FU_LMPD_19_4fe4s_c  <= 1
 0  <= BU_LMPD_19_4fe4s_c  <= 1
 0  <= FU_LMPD_20_4fe4s_c  <= 1
 0  <= BU_LMPD_20_4fe4s_c  <= 1
 0  <= FU_LMPD_21_4fe4s_c  <= 1
 0  <= BU_LMPD_21_4fe4s_c  <= 1
 -1000  <= DG_LMPD_22_5mthf_c  <= 1000
 -2.257353e+02  <= DGo_LMPD_22_5mthf_c  <= -2.033803e+02
 0  <= FU_LMPD_22_5mthf_c  <= 1
 0  <= BU_LMPD_22_5mthf_c  <= 1
 -1000  <= DG_LMPD_23_5mthf_c  <= 1000
 -2.253058e+02  <= DGo_LMPD_23_5mthf_c  <= -2.029507e+02
 0  <= FU_LMPD_23_5mthf_c  <= 1
 0  <= BU_LMPD_23_5mthf_c  <= 1
 -1000  <= DG_LMPD_24_accoa_c  <= 1000
 -1.447283e+02  <= DGo_LMPD_24_accoa_c  <= -1.253338e+02
 0  <= FU_LMPD_24_accoa_c  <= 1
 0  <= BU_LMPD_24_accoa_c  <= 1
 -1000  <= DG_LMPD_25_accoa_c  <= 1000
 -1.438692e+02  <= DGo_LMPD_25_accoa_c  <= -1.244747e+02
 0  <= FU_LMPD_25_accoa_c  <= 1
 0  <= BU_LMPD_25_accoa_c  <= 1
 0  <= FU_LMPD_26_adocbl_c  <= 1
 0  <= BU_LMPD_26_adocbl_c  <= 1
 0  <= FU_LMPD_27_adocbl_c  <= 1
 0  <= BU_LMPD_27_adocbl_c  <= 1
 0  <= FU_LMPD_28_adocbl_c  <= 1
 0  <= BU_LMPD_28_adocbl_c  <= 1
 0  <= FU_LMPD_29_adocbl_c  <= 1
 0  <= BU_LMPD_29_adocbl_c  <= 1
 -1000  <= DG_LMPD_30_ala-L_c  <= 1000
 2.790665e+00  <= DGo_LMPD_30_ala-L_c  <= 6.790665e+00
 0  <= FU_LMPD_30_ala-L_c  <= 1
 0  <= BU_LMPD_30_ala-L_c  <= 1
 -1000  <= DG_LMPD_31_amet_c  <= 1000
 -1.131925e+02  <= DGo_LMPD_31_amet_c  <= -9.155766e+01
 0  <= FU_LMPD_31_amet_c  <= 1
 0  <= BU_LMPD_31_amet_c  <= 1
 -1000  <= DG_LMPD_32_amet_c  <= 1000
 -1.123334e+02  <= DGo_LMPD_32_amet_c  <= -9.069857e+01
 0  <= FU_LMPD_32_amet_c  <= 1
 0  <= BU_LMPD_32_amet_c  <= 1
 -1000  <= DG_LMPD_33_arg-L_c  <= 1000
 3.142159e+00  <= DGo_LMPD_33_arg-L_c  <= 8.732701e+00
 0  <= FU_LMPD_33_arg-L_c  <= 1
 0  <= BU_LMPD_33_arg-L_c  <= 1
 -1000  <= DG_LMPD_34_asn-L_c  <= 1000
 -1.864113e+00  <= DGo_LMPD_34_asn-L_c  <= -3.663127e-02
 0  <= FU_LMPD_34_asn-L_c  <= 1
 0  <= BU_LMPD_34_asn-L_c  <= 1
 0  <= FU_LMPD_35_bmocogdp_c  <= 1
 0  <= BU_LMPD_35_bmocogdp_c  <= 1
 0  <= FU_LMPD_36_bmocogdp_c  <= 1
 0  <= BU_LMPD_36_bmocogdp_c  <= 1
 0  <= FU_LMPD_37_bmocogdp_c  <= 1
 0  <= BU_LMPD_37_bmocogdp_c  <= 1
 0  <= FU_LMPD_38_bmocogdp_c  <= 1
 0  <= BU_LMPD_38_bmocogdp_c  <= 1
 -1000  <= DG_LMPD_39_btn_c  <= 1000
 -3.136089e+02  <= DGo_LMPD_39_btn_c  <= -2.727382e+02
 0  <= FU_LMPD_39_btn_c  <= 1
 0  <= BU_LMPD_39_btn_c  <= 1
 -1000  <= DG_LMPD_40_btn_c  <= 1000
 -3.218565e+02  <= DGo_LMPD_40_btn_c  <= -2.800159e+02
 0  <= FU_LMPD_40_btn_c  <= 1
 0  <= BU_LMPD_40_btn_c  <= 1
 -1000  <= DG_LMPD_41_btn_c  <= 1000
 -3.123203e+02  <= DGo_LMPD_41_btn_c  <= -2.714496e+02
 0  <= FU_LMPD_41_btn_c  <= 1
 0  <= BU_LMPD_41_btn_c  <= 1
 -1000  <= DG_LMPD_42_btn_c  <= 1000
 -3.205679e+02  <= DGo_LMPD_42_btn_c  <= -2.787272e+02
 0  <= FU_LMPD_42_btn_c  <= 1
 0  <= BU_LMPD_42_btn_c  <= 1
 -1000  <= DG_LMPD_43_btn_c  <= 1000
 -4.622946e+02  <= DGo_LMPD_43_btn_c  <= -3.918699e+02
 0  <= FU_LMPD_43_btn_c  <= 1
 0  <= BU_LMPD_43_btn_c  <= 1
 -1000  <= DG_LMPD_44_btn_c  <= 1000
 -4.705619e+02  <= DGo_LMPD_44_btn_c  <= -3.991278e+02
 0  <= FU_LMPD_44_btn_c  <= 1
 0  <= BU_LMPD_44_btn_c  <= 1
 -1000  <= DG_LMPD_45_btn_c  <= 1000
 -4.597173e+02  <= DGo_LMPD_45_btn_c  <= -3.892926e+02
 0  <= FU_LMPD_45_btn_c  <= 1
 0  <= BU_LMPD_45_btn_c  <= 1
 -1000  <= DG_LMPD_47_btn_c  <= 1000
 -4.679847e+02  <= DGo_LMPD_47_btn_c  <= -3.965505e+02
 0  <= FU_LMPD_47_btn_c  <= 1
 0  <= BU_LMPD_47_btn_c  <= 1
 -1000  <= DG_LMPD_48_ca2_c  <= 1000
 -2.386934e+01  <= DGo_LMPD_48_ca2_c  <= -1.733729e+01
 0  <= FU_LMPD_48_ca2_c  <= 1
 0  <= BU_LMPD_48_ca2_c  <= 1
 -1000  <= DG_LMPD_49_ca2_c  <= 1000
 -2.906190e+01  <= DGo_LMPD_49_ca2_c  <= -2.171668e+01
 0  <= FU_LMPD_49_ca2_c  <= 1
 0  <= BU_LMPD_49_ca2_c  <= 1
 -1000  <= DG_LMPD_50_ca2_c  <= 1000
 -2.253528e+01  <= DGo_LMPD_50_ca2_c  <= -1.822995e+01
 0  <= FU_LMPD_50_ca2_c  <= 1
 0  <= BU_LMPD_50_ca2_c  <= 1
 -1000  <= DG_LMPD_51_chor_c  <= 1000
 -5.604764e+01  <= DGo_LMPD_51_chor_c  <= -5.050486e+01
 0  <= FU_LMPD_51_chor_c  <= 1
 0  <= BU_LMPD_51_chor_c  <= 1
 -1000  <= DG_LMPD_52_clpn160_p  <= 1000
 -5.532432e+02  <= DGo_LMPD_52_clpn160_p  <= -4.127467e+02
 0  <= FU_LMPD_52_clpn160_p  <= 1
 0  <= BU_LMPD_52_clpn160_p  <= 1
 -1000  <= DG_LMPD_53_clpn160_p  <= 1000
 -5.547913e+02  <= DGo_LMPD_53_clpn160_p  <= -4.142948e+02
 0  <= FU_LMPD_53_clpn160_p  <= 1
 0  <= BU_LMPD_53_clpn160_p  <= 1
 -1000  <= DG_LMPD_54_clpn160_p  <= 1000
 -5.563394e+02  <= DGo_LMPD_54_clpn160_p  <= -4.158430e+02
 0  <= FU_LMPD_54_clpn160_p  <= 1
 0  <= BU_LMPD_54_clpn160_p  <= 1
 -1000  <= DG_LMPD_55_clpn160_p  <= 1000
 -5.578876e+02  <= DGo_LMPD_55_clpn160_p  <= -4.173911e+02
 0  <= FU_LMPD_55_clpn160_p  <= 1
 0  <= BU_LMPD_55_clpn160_p  <= 1
 -1000  <= DG_LMPD_56_clpn160_p  <= 1000
 -5.594357e+02  <= DGo_LMPD_56_clpn160_p  <= -4.189392e+02
 0  <= FU_LMPD_56_clpn160_p  <= 1
 0  <= BU_LMPD_56_clpn160_p  <= 1
 -1000  <= DG_LMPD_57_clpn160_p  <= 1000
 -5.501469e+02  <= DGo_LMPD_57_clpn160_p  <= -4.096504e+02
 0  <= FU_LMPD_57_clpn160_p  <= 1
 0  <= BU_LMPD_57_clpn160_p  <= 1
 -1000  <= DG_LMPD_58_clpn160_p  <= 1000
 -5.485988e+02  <= DGo_LMPD_58_clpn160_p  <= -4.081023e+02
 0  <= FU_LMPD_58_clpn160_p  <= 1
 0  <= BU_LMPD_58_clpn160_p  <= 1
 -1000  <= DG_LMPD_59_clpn160_p  <= 1000
 -5.516950e+02  <= DGo_LMPD_59_clpn160_p  <= -4.111986e+02
 0  <= FU_LMPD_59_clpn160_p  <= 1
 0  <= BU_LMPD_59_clpn160_p  <= 1
 -1000  <= DG_LMPD_60_clpn160_p  <= 1000
 -3.816753e+02  <= DGo_LMPD_60_clpn160_p  <= -2.427406e+02
 0  <= FU_LMPD_60_clpn160_p  <= 1
 0  <= BU_LMPD_60_clpn160_p  <= 1
 -1000  <= DG_LMPD_61_clpn160_p  <= 1000
 -3.925122e+02  <= DGo_LMPD_61_clpn160_p  <= -2.535775e+02
 0  <= FU_LMPD_61_clpn160_p  <= 1
 0  <= BU_LMPD_61_clpn160_p  <= 1
 -1000  <= DG_LMPD_62_clpn161_p  <= 1000
 -4.793457e+02  <= DGo_LMPD_62_clpn161_p  <= -3.457865e+02
 0  <= FU_LMPD_62_clpn161_p  <= 1
 0  <= BU_LMPD_62_clpn161_p  <= 1
 -1000  <= DG_LMPD_63_clpn161_p  <= 1000
 -4.808938e+02  <= DGo_LMPD_63_clpn161_p  <= -3.473346e+02
 0  <= FU_LMPD_63_clpn161_p  <= 1
 0  <= BU_LMPD_63_clpn161_p  <= 1
 -1000  <= DG_LMPD_64_clpn161_p  <= 1000
 -4.824420e+02  <= DGo_LMPD_64_clpn161_p  <= -3.488827e+02
 0  <= FU_LMPD_64_clpn161_p  <= 1
 0  <= BU_LMPD_64_clpn161_p  <= 1
 -1000  <= DG_LMPD_65_clpn161_p  <= 1000
 -4.839901e+02  <= DGo_LMPD_65_clpn161_p  <= -3.504308e+02
 0  <= FU_LMPD_65_clpn161_p  <= 1
 0  <= BU_LMPD_65_clpn161_p  <= 1
 -1000  <= DG_LMPD_66_clpn161_p  <= 1000
 -4.762494e+02  <= DGo_LMPD_66_clpn161_p  <= -3.426902e+02
 0  <= FU_LMPD_66_clpn161_p  <= 1
 0  <= BU_LMPD_66_clpn161_p  <= 1
 -1000  <= DG_LMPD_67_clpn161_p  <= 1000
 -4.747013e+02  <= DGo_LMPD_67_clpn161_p  <= -3.411421e+02
 0  <= FU_LMPD_67_clpn161_p  <= 1
 0  <= BU_LMPD_67_clpn161_p  <= 1
 -1000  <= DG_LMPD_68_clpn161_p  <= 1000
 -4.777976e+02  <= DGo_LMPD_68_clpn161_p  <= -3.442383e+02
 0  <= FU_LMPD_68_clpn161_p  <= 1
 0  <= BU_LMPD_68_clpn161_p  <= 1
 -1000  <= DG_LMPD_69_clpn181_p  <= 1000
 -5.585181e+02  <= DGo_LMPD_69_clpn181_p  <= -4.056497e+02
 0  <= FU_LMPD_69_clpn181_p  <= 1
 0  <= BU_LMPD_69_clpn181_p  <= 1
 -1000  <= DG_LMPD_70_clpn181_p  <= 1000
 -5.600662e+02  <= DGo_LMPD_70_clpn181_p  <= -4.071978e+02
 0  <= FU_LMPD_70_clpn181_p  <= 1
 0  <= BU_LMPD_70_clpn181_p  <= 1
 -1000  <= DG_LMPD_71_clpn181_p  <= 1000
 -5.616144e+02  <= DGo_LMPD_71_clpn181_p  <= -4.087459e+02
 0  <= FU_LMPD_71_clpn181_p  <= 1
 0  <= BU_LMPD_71_clpn181_p  <= 1
 -1000  <= DG_LMPD_72_clpn181_p  <= 1000
 -5.631625e+02  <= DGo_LMPD_72_clpn181_p  <= -4.102940e+02
 0  <= FU_LMPD_72_clpn181_p  <= 1
 0  <= BU_LMPD_72_clpn181_p  <= 1
 -1000  <= DG_LMPD_73_clpn181_p  <= 1000
 -5.647106e+02  <= DGo_LMPD_73_clpn181_p  <= -4.118422e+02
 0  <= FU_LMPD_73_clpn181_p  <= 1
 0  <= BU_LMPD_73_clpn181_p  <= 1
 -1000  <= DG_LMPD_74_clpn181_p  <= 1000
 -5.554218e+02  <= DGo_LMPD_74_clpn181_p  <= -4.025534e+02
 0  <= FU_LMPD_74_clpn181_p  <= 1
 0  <= BU_LMPD_74_clpn181_p  <= 1
 -1000  <= DG_LMPD_75_clpn181_p  <= 1000
 -5.538737e+02  <= DGo_LMPD_75_clpn181_p  <= -4.010053e+02
 0  <= FU_LMPD_75_clpn181_p  <= 1
 0  <= BU_LMPD_75_clpn181_p  <= 1
 -1000  <= DG_LMPD_76_clpn181_p  <= 1000
 -5.569700e+02  <= DGo_LMPD_76_clpn181_p  <= -4.041015e+02
 0  <= FU_LMPD_76_clpn181_p  <= 1
 0  <= BU_LMPD_76_clpn181_p  <= 1
 -1000  <= DG_LMPD_77_clpn181_p  <= 1000
 -3.573792e+02  <= DGo_LMPD_77_clpn181_p  <= -2.062098e+02
 0  <= FU_LMPD_77_clpn181_p  <= 1
 0  <= BU_LMPD_77_clpn181_p  <= 1
 -1000  <= DG_LMPD_78_clpn181_p  <= 1000
 -3.682161e+02  <= DGo_LMPD_78_clpn181_p  <= -2.170467e+02
 0  <= FU_LMPD_78_clpn181_p  <= 1
 0  <= BU_LMPD_78_clpn181_p  <= 1
 -1000  <= DG_LMPD_79_colipa_e  <= 1000
 -6.952912e+02  <= DGo_LMPD_79_colipa_e  <= -5.162580e+02
 0  <= FU_LMPD_79_colipa_e  <= 1
 0  <= BU_LMPD_79_colipa_e  <= 1
 -1000  <= DG_LMPD_80_colipa_e  <= 1000
 -6.956783e+02  <= DGo_LMPD_80_colipa_e  <= -5.166450e+02
 0  <= FU_LMPD_80_colipa_e  <= 1
 0  <= BU_LMPD_80_colipa_e  <= 1
 -1000  <= DG_LMPD_81_colipa_e  <= 1000
 -6.976134e+02  <= DGo_LMPD_81_colipa_e  <= -5.185802e+02
 0  <= FU_LMPD_81_colipa_e  <= 1
 0  <= BU_LMPD_81_colipa_e  <= 1
 -1000  <= DG_LMPD_82_colipa_e  <= 1000
 -6.980005e+02  <= DGo_LMPD_82_colipa_e  <= -5.189672e+02
 0  <= FU_LMPD_82_colipa_e  <= 1
 0  <= BU_LMPD_82_colipa_e  <= 1
 -1000  <= DG_LMPD_83_colipa_e  <= 1000
 -6.999356e+02  <= DGo_LMPD_83_colipa_e  <= -5.209024e+02
 0  <= FU_LMPD_83_colipa_e  <= 1
 0  <= BU_LMPD_83_colipa_e  <= 1
 -1000  <= DG_LMPD_84_colipa_e  <= 1000
 -7.003226e+02  <= DGo_LMPD_84_colipa_e  <= -5.212894e+02
 0  <= FU_LMPD_84_colipa_e  <= 1
 0  <= BU_LMPD_84_colipa_e  <= 1
 -1000  <= DG_LMPD_85_colipa_e  <= 1000
 -7.022578e+02  <= DGo_LMPD_85_colipa_e  <= -5.232246e+02
 0  <= FU_LMPD_85_colipa_e  <= 1
 0  <= BU_LMPD_85_colipa_e  <= 1
 -1000  <= DG_LMPD_86_colipa_e  <= 1000
 -7.026448e+02  <= DGo_LMPD_86_colipa_e  <= -5.236116e+02
 0  <= FU_LMPD_86_colipa_e  <= 1
 0  <= BU_LMPD_86_colipa_e  <= 1
 -1000  <= DG_LMPD_87_colipa_e  <= 1000
 -6.929690e+02  <= DGo_LMPD_87_colipa_e  <= -5.139358e+02
 0  <= FU_LMPD_87_colipa_e  <= 1
 0  <= BU_LMPD_87_colipa_e  <= 1
 -1000  <= DG_LMPD_88_colipa_e  <= 1000
 -6.906468e+02  <= DGo_LMPD_88_colipa_e  <= -5.116136e+02
 0  <= FU_LMPD_88_colipa_e  <= 1
 0  <= BU_LMPD_88_colipa_e  <= 1
 -1000  <= DG_LMPD_89_colipa_e  <= 1000
 -6.933561e+02  <= DGo_LMPD_89_colipa_e  <= -5.143228e+02
 0  <= FU_LMPD_89_colipa_e  <= 1
 0  <= BU_LMPD_89_colipa_e  <= 1
 -1000  <= DG_LMPD_90_colipa_e  <= 1000
 -6.910339e+02  <= DGo_LMPD_90_colipa_e  <= -5.120006e+02
 0  <= FU_LMPD_90_colipa_e  <= 1
 0  <= BU_LMPD_90_colipa_e  <= 1
 -1000  <= DG_LMPD_91_ctp_c  <= 1000
 -2.415188e+01  <= DGo_LMPD_91_ctp_c  <= -1.926573e+01
 0  <= FU_LMPD_91_ctp_c  <= 1
 0  <= BU_LMPD_91_ctp_c  <= 1
 -1000  <= DG_LMPD_92_cys-L_c  <= 1000
 -5.677185e+01  <= DGo_LMPD_92_cys-L_c  <= -5.004789e+01
 0  <= FU_LMPD_92_cys-L_c  <= 1
 0  <= BU_LMPD_92_cys-L_c  <= 1
 -1000  <= DG_LMPD_93_cys-L_c  <= 1000
 -5.634231e+01  <= DGo_LMPD_93_cys-L_c  <= -4.961834e+01
 0  <= FU_LMPD_93_cys-L_c  <= 1
 0  <= BU_LMPD_93_cys-L_c  <= 1
 0  <= FU_LMPD_94_datp_c  <= 1
 0  <= BU_LMPD_94_datp_c  <= 1
 0  <= FU_LMPD_95_datp_c  <= 1
 0  <= BU_LMPD_95_datp_c  <= 1
 0  <= FU_LMPD_96_datp_c  <= 1
 0  <= BU_LMPD_96_datp_c  <= 1
 0  <= FU_LMPD_97_datp_c  <= 1
 0  <= BU_LMPD_97_datp_c  <= 1
 0  <= FU_LMPD_98_dctp_c  <= 1
 0  <= BU_LMPD_98_dctp_c  <= 1
 0  <= FU_LMPD_99_dgtp_c  <= 1
 0  <= BU_LMPD_99_dgtp_c  <= 1
 0  <= FU_LMPD_100_dgtp_c  <= 1
 0  <= BU_LMPD_100_dgtp_c  <= 1
 -1000  <= DG_LMPD_101_dgtp_c  <= 1000
 -6.700300e+01  <= DGo_LMPD_101_dgtp_c  <= -4.889688e+01
 0  <= FU_LMPD_101_dgtp_c  <= 1
 0  <= BU_LMPD_101_dgtp_c  <= 1
 0  <= FU_LMPD_102_dgtp_c  <= 1
 0  <= BU_LMPD_102_dgtp_c  <= 1
 0  <= FU_LMPD_103_dgtp_c  <= 1
 0  <= BU_LMPD_103_dgtp_c  <= 1
 -1000  <= DG_LMPD_104_dgtp_c  <= 1000
 -7.323616e+01  <= DGo_LMPD_104_dgtp_c  <= -5.768088e+01
 0  <= FU_LMPD_104_dgtp_c  <= 1
 0  <= BU_LMPD_104_dgtp_c  <= 1
 0  <= FU_LMPD_105_dttp_c  <= 1
 0  <= BU_LMPD_105_dttp_c  <= 1
 -1000  <= DG_LMPD_106_dttp_c  <= 1000
 -5.159644e+01  <= DGo_LMPD_106_dttp_c  <= -4.404646e+01
 0  <= FU_LMPD_106_dttp_c  <= 1
 0  <= BU_LMPD_106_dttp_c  <= 1
 0  <= FU_LMPD_107_enter_c  <= 1
 0  <= BU_LMPD_107_enter_c  <= 1
 -1000  <= DG_LMPD_108_fad_c  <= 1000
 -2.088134e+02  <= DGo_LMPD_108_fad_c  <= -1.783480e+02
 0  <= FU_LMPD_108_fad_c  <= 1
 0  <= BU_LMPD_108_fad_c  <= 1
 -1000  <= DG_LMPD_109_fad_c  <= 1000
 -2.079543e+02  <= DGo_LMPD_109_fad_c  <= -1.774889e+02
 0  <= FU_LMPD_109_fad_c  <= 1
 0  <= BU_LMPD_109_fad_c  <= 1
 -1000  <= DG_LMPD_110_fad_c  <= 1000
 -1.958534e+02  <= DGo_LMPD_110_fad_c  <= -1.612737e+02
 0  <= FU_LMPD_110_fad_c  <= 1
 0  <= BU_LMPD_110_fad_c  <= 1
 -1000  <= DG_LMPD_111_fad_c  <= 1000
 -1.949944e+02  <= DGo_LMPD_111_fad_c  <= -1.604146e+02
 0  <= FU_LMPD_111_fad_c  <= 1
 0  <= BU_LMPD_111_fad_c  <= 1
 -1000  <= DG_LMPD_112_GAM_c  <= 1000
 -3.809952e+01  <= DGo_LMPD_112_GAM_c  <= -2.254223e+01
 0  <= FU_LMPD_112_GAM_c  <= 1
 0  <= BU_LMPD_112_GAM_c  <= 1
 -1000  <= DG_LMPD_113_GAM_c  <= 1000
 -3.766997e+01  <= DGo_LMPD_113_GAM_c  <= -2.211268e+01
 0  <= FU_LMPD_113_GAM_c  <= 1
 0  <= BU_LMPD_113_GAM_c  <= 1
 -1000  <= DG_LMPD_114_GAM_c  <= 1000
 -4.414184e+01  <= DGo_LMPD_114_GAM_c  <= -3.065798e+01
 0  <= FU_LMPD_114_GAM_c  <= 1
 0  <= BU_LMPD_114_GAM_c  <= 1
 -1000  <= DG_LMPD_115_GAM_c  <= 1000
 -4.457138e+01  <= DGo_LMPD_115_GAM_c  <= -3.108752e+01
 0  <= FU_LMPD_115_GAM_c  <= 1
 0  <= BU_LMPD_115_GAM_c  <= 1
 -1000  <= DG_LMPD_116_gln-L_c  <= 1000
 1.636171e+00  <= DGo_LMPD_116_gln-L_c  <= 3.414114e+00
 0  <= FU_LMPD_116_gln-L_c  <= 1
 0  <= BU_LMPD_116_gln-L_c  <= 1
 -1000  <= DG_LMPD_117_gly_c  <= 1000
 -3.158528e+00  <= DGo_LMPD_117_gly_c  <= 1.728792e+00
 0  <= FU_LMPD_117_gly_c  <= 1
 0  <= BU_LMPD_117_gly_c  <= 1
 -1000  <= DG_LMPD_118_gthrd_c  <= 1000
 -6.656355e+01  <= DGo_LMPD_118_gthrd_c  <= -5.829421e+01
 0  <= FU_LMPD_118_gthrd_c  <= 1
 0  <= BU_LMPD_118_gthrd_c  <= 1
 -1000  <= DG_LMPD_119_gthrd_c  <= 1000
 -6.613401e+01  <= DGo_LMPD_119_gthrd_c  <= -5.786467e+01
 0  <= FU_LMPD_119_gthrd_c  <= 1
 0  <= BU_LMPD_119_gthrd_c  <= 1
 -1000  <= DG_LMPD_120_gtp_c  <= 1000
 -5.442179e+01  <= DGo_LMPD_120_gtp_c  <= -3.691388e+01
 0  <= FU_LMPD_120_gtp_c  <= 1
 0  <= BU_LMPD_120_gtp_c  <= 1
 -1000  <= DG_LMPD_121_gtp_c  <= 1000
 -5.399225e+01  <= DGo_LMPD_121_gtp_c  <= -3.648434e+01
 0  <= FU_LMPD_121_gtp_c  <= 1
 0  <= BU_LMPD_121_gtp_c  <= 1
 -1000  <= DG_LMPD_122_gtp_c  <= 1000
 -6.062973e+01  <= DGo_LMPD_122_gtp_c  <= -4.572310e+01
 0  <= FU_LMPD_122_gtp_c  <= 1
 0  <= BU_LMPD_122_gtp_c  <= 1
 -1000  <= DG_LMPD_123_gtp_c  <= 1000
 -6.020019e+01  <= DGo_LMPD_123_gtp_c  <= -4.529355e+01
 0  <= FU_LMPD_123_gtp_c  <= 1
 0  <= BU_LMPD_123_gtp_c  <= 1
 0  <= FU_LMPD_124_hemeO_c  <= 1
 0  <= BU_LMPD_124_hemeO_c  <= 1
 0  <= FU_LMPD_125_hemeO_c  <= 1
 0  <= BU_LMPD_125_hemeO_c  <= 1
 -1000  <= DG_LMPD_126_his-L_c  <= 1000
 -6.407459e+01  <= DGo_LMPD_126_his-L_c  <= -5.659146e+01
 0  <= FU_LMPD_126_his-L_c  <= 1
 0  <= BU_LMPD_126_his-L_c  <= 1
 -1000  <= DG_LMPD_127_his-L_c  <= 1000
 -6.364504e+01  <= DGo_LMPD_127_his-L_c  <= -5.616192e+01
 0  <= FU_LMPD_127_his-L_c  <= 1
 0  <= BU_LMPD_127_his-L_c  <= 1
 -1000  <= DG_LMPD_128_ile-L_c  <= 1000
 -2.671205e+01  <= DGo_LMPD_128_ile-L_c  <= -1.992807e+01
 0  <= FU_LMPD_128_ile-L_c  <= 1
 0  <= BU_LMPD_128_ile-L_c  <= 1
 -1000  <= DG_LMPD_129_leu-L_c  <= 1000
 -1.868890e+01  <= DGo_LMPD_129_leu-L_c  <= -1.560268e+01
 0  <= FU_LMPD_129_leu-L_c  <= 1
 0  <= BU_LMPD_129_leu-L_c  <= 1
 0  <= FU_LMPD_130_lipopb_c  <= 1
 0  <= BU_LMPD_130_lipopb_c  <= 1
 0  <= FU_LMPD_131_lipopb_c  <= 1
 0  <= BU_LMPD_131_lipopb_c  <= 1
 0  <= FU_LMPD_132_lipopb_c  <= 1
 0  <= BU_LMPD_132_lipopb_c  <= 1
 0  <= FU_LMPD_133_lipopb_c  <= 1
 0  <= BU_LMPD_133_lipopb_c  <= 1
 -1000  <= DG_LMPD_134_lys-L_c  <= 1000
 -8.039323e+00  <= DGo_LMPD_134_lys-L_c  <= -3.058055e+00
 0  <= FU_LMPD_134_lys-L_c  <= 1
 0  <= BU_LMPD_134_lys-L_c  <= 1
 -1000  <= DG_LMPD_135_malcoa_c  <= 1000
 -1.506379e+02  <= DGo_LMPD_135_malcoa_c  <= -1.304030e+02
 0  <= FU_LMPD_135_malcoa_c  <= 1
 0  <= BU_LMPD_135_malcoa_c  <= 1
 -1000  <= DG_LMPD_136_malcoa_c  <= 1000
 -1.497788e+02  <= DGo_LMPD_136_malcoa_c  <= -1.295439e+02
 0  <= FU_LMPD_136_malcoa_c  <= 1
 0  <= BU_LMPD_136_malcoa_c  <= 1
 -1000  <= DG_LMPD_137_met-L_c  <= 1000
 -5.859629e+01  <= DGo_LMPD_137_met-L_c  <= -4.759830e+01
 0  <= FU_LMPD_137_met-L_c  <= 1
 0  <= BU_LMPD_137_met-L_c  <= 1
 -1000  <= DG_LMPD_138_met-L_c  <= 1000
 -5.816674e+01  <= DGo_LMPD_138_met-L_c  <= -4.716876e+01
 0  <= FU_LMPD_138_met-L_c  <= 1
 0  <= BU_LMPD_138_met-L_c  <= 1
 -1000  <= DG_LMPD_139_mlthf_c  <= 1000
 -2.235946e+02  <= DGo_LMPD_139_mlthf_c  <= -2.008792e+02
 0  <= FU_LMPD_139_mlthf_c  <= 1
 0  <= BU_LMPD_139_mlthf_c  <= 1
 -1000  <= DG_LMPD_140_mlthf_c  <= 1000
 -2.231651e+02  <= DGo_LMPD_140_mlthf_c  <= -2.004496e+02
 0  <= FU_LMPD_140_mlthf_c  <= 1
 0  <= BU_LMPD_140_mlthf_c  <= 1
 0  <= FU_LMPD_141_mococdp_c  <= 1
 0  <= BU_LMPD_141_mococdp_c  <= 1
 0  <= FU_LMPD_142_mococdp_c  <= 1
 0  <= BU_LMPD_142_mococdp_c  <= 1
 0  <= FU_LMPD_143_mococdp_c  <= 1
 0  <= BU_LMPD_143_mococdp_c  <= 1
 0  <= FU_LMPD_144_mococdp_c  <= 1
 0  <= BU_LMPD_144_mococdp_c  <= 1
 0  <= FU_LMPD_145_mococdp_c  <= 1
 0  <= BU_LMPD_145_mococdp_c  <= 1
 0  <= FU_LMPD_146_mococdp_c  <= 1
 0  <= BU_LMPD_146_mococdp_c  <= 1
 0  <= FU_LMPD_147_mococdp_c  <= 1
 0  <= BU_LMPD_147_mococdp_c  <= 1
 0  <= FU_LMPD_148_mococdp_c  <= 1
 0  <= BU_LMPD_148_mococdp_c  <= 1
 0  <= FU_LMPD_149_mocogdp_c  <= 1
 0  <= BU_LMPD_149_mocogdp_c  <= 1
 0  <= FU_LMPD_150_mocogdp_c  <= 1
 0  <= BU_LMPD_150_mocogdp_c  <= 1
 0  <= FU_LMPD_151_mocogdp_c  <= 1
 0  <= BU_LMPD_151_mocogdp_c  <= 1
 0  <= FU_LMPD_152_mocogdp_c  <= 1
 0  <= BU_LMPD_152_mocogdp_c  <= 1
 0  <= FU_LMPD_153_mql8_c  <= 1
 0  <= BU_LMPD_153_mql8_c  <= 1
 -1000  <= DG_LMPD_154_murein3p3p_p  <= 1000
 -2.689941e+02  <= DGo_LMPD_154_murein3p3p_p  <= -2.367532e+02
 0  <= FU_LMPD_154_murein3p3p_p  <= 1
 0  <= BU_LMPD_154_murein3p3p_p  <= 1
 -1000  <= DG_LMPD_155_murein3px4p_p  <= 1000
 -2.795191e+02  <= DGo_LMPD_155_murein3px4p_p  <= -2.458663e+02
 0  <= FU_LMPD_155_murein3px4p_p  <= 1
 0  <= BU_LMPD_155_murein3px4p_p  <= 1
 -1000  <= DG_LMPD_156_murein4p4p_p  <= 1000
 -2.751931e+02  <= DGo_LMPD_156_murein4p4p_p  <= -2.420911e+02
 0  <= FU_LMPD_156_murein4p4p_p  <= 1
 0  <= BU_LMPD_156_murein4p4p_p  <= 1
 -1000  <= DG_LMPD_157_murein4px4p_p  <= 1000
 -2.749375e+02  <= DGo_LMPD_157_murein4px4p_p  <= -2.412847e+02
 0  <= FU_LMPD_157_murein4px4p_p  <= 1
 0  <= BU_LMPD_157_murein4px4p_p  <= 1
 -1000  <= DG_LMPD_158_murein4px4px4p_p  <= 1000
 -3.959061e+02  <= DGo_LMPD_158_murein4px4px4p_p  <= -3.450730e+02
 0  <= FU_LMPD_158_murein4px4px4p_p  <= 1
 0  <= BU_LMPD_158_murein4px4px4p_p  <= 1
 -1000  <= DG_LMPD_159_nad_c  <= 1000
 -5.320017e+01  <= DGo_LMPD_159_nad_c  <= -3.631099e+01
 0  <= FU_LMPD_159_nad_c  <= 1
 0  <= BU_LMPD_159_nad_c  <= 1
 -1000  <= DG_LMPD_160_nad_c  <= 1000
 -5.277063e+01  <= DGo_LMPD_160_nad_c  <= -3.588145e+01
 0  <= FU_LMPD_160_nad_c  <= 1
 0  <= BU_LMPD_160_nad_c  <= 1
 -1000  <= DG_LMPD_161_nad_c  <= 1000
 -6.006356e+01  <= DGo_LMPD_161_nad_c  <= -4.446476e+01
 0  <= FU_LMPD_161_nad_c  <= 1
 0  <= BU_LMPD_161_nad_c  <= 1
 -1000  <= DG_LMPD_162_nad_c  <= 1000
 -5.963402e+01  <= DGo_LMPD_162_nad_c  <= -4.403521e+01
 0  <= FU_LMPD_162_nad_c  <= 1
 0  <= BU_LMPD_162_nad_c  <= 1
 -1000  <= DG_LMPD_163_pe160_c  <= 1000
 -2.237651e+02  <= DGo_LMPD_163_pe160_c  <= -1.552543e+02
 0  <= FU_LMPD_163_pe160_c  <= 1
 0  <= BU_LMPD_163_pe160_c  <= 1
 -1000  <= DG_LMPD_164_pe160_c  <= 1000
 -2.245392e+02  <= DGo_LMPD_164_pe160_c  <= -1.560283e+02
 0  <= FU_LMPD_164_pe160_c  <= 1
 0  <= BU_LMPD_164_pe160_c  <= 1
 -1000  <= DG_LMPD_165_pe160_c  <= 1000
 -2.206689e+02  <= DGo_LMPD_165_pe160_c  <= -1.521580e+02
 0  <= FU_LMPD_165_pe160_c  <= 1
 0  <= BU_LMPD_165_pe160_c  <= 1
 -1000  <= DG_LMPD_166_pe160_c  <= 1000
 -2.214429e+02  <= DGo_LMPD_166_pe160_c  <= -1.529321e+02
 0  <= FU_LMPD_166_pe160_c  <= 1
 0  <= BU_LMPD_166_pe160_c  <= 1
 -1000  <= DG_LMPD_167_pe160_c  <= 1000
 -2.222170e+02  <= DGo_LMPD_167_pe160_c  <= -1.537061e+02
 0  <= FU_LMPD_167_pe160_c  <= 1
 0  <= BU_LMPD_167_pe160_c  <= 1
 -1000  <= DG_LMPD_168_pe160_c  <= 1000
 -2.229911e+02  <= DGo_LMPD_168_pe160_c  <= -1.544802e+02
 0  <= FU_LMPD_168_pe160_c  <= 1
 0  <= BU_LMPD_168_pe160_c  <= 1
 -1000  <= DG_LMPD_169_pe160_c  <= 1000
 -2.191207e+02  <= DGo_LMPD_169_pe160_c  <= -1.506099e+02
 0  <= FU_LMPD_169_pe160_c  <= 1
 0  <= BU_LMPD_169_pe160_c  <= 1
 -1000  <= DG_LMPD_170_pe160_c  <= 1000
 -2.198948e+02  <= DGo_LMPD_170_pe160_c  <= -1.513839e+02
 0  <= FU_LMPD_170_pe160_c  <= 1
 0  <= BU_LMPD_170_pe160_c  <= 1
 -1000  <= DG_LMPD_171_pe160_c  <= 1000
 -1.356490e+02  <= DGo_LMPD_171_pe160_c  <= -6.793900e+01
 0  <= FU_LMPD_171_pe160_c  <= 1
 0  <= BU_LMPD_171_pe160_c  <= 1
 -1000  <= DG_LMPD_172_pe160_c  <= 1000
 -1.410674e+02  <= DGo_LMPD_172_pe160_c  <= -7.335745e+01
 0  <= FU_LMPD_172_pe160_c  <= 1
 0  <= BU_LMPD_172_pe160_c  <= 1
 -1000  <= DG_LMPD_173_pe161_c  <= 1000
 -1.895631e+02  <= DGo_LMPD_173_pe161_c  <= -1.244808e+02
 0  <= FU_LMPD_173_pe161_c  <= 1
 0  <= BU_LMPD_173_pe161_c  <= 1
 -1000  <= DG_LMPD_174_pe161_c  <= 1000
 -1.864668e+02  <= DGo_LMPD_174_pe161_c  <= -1.213845e+02
 0  <= FU_LMPD_174_pe161_c  <= 1
 0  <= BU_LMPD_174_pe161_c  <= 1
 -1000  <= DG_LMPD_175_pe161_c  <= 1000
 -1.872409e+02  <= DGo_LMPD_175_pe161_c  <= -1.221586e+02
 0  <= FU_LMPD_175_pe161_c  <= 1
 0  <= BU_LMPD_175_pe161_c  <= 1
 -1000  <= DG_LMPD_176_pe161_c  <= 1000
 -1.880149e+02  <= DGo_LMPD_176_pe161_c  <= -1.229327e+02
 0  <= FU_LMPD_176_pe161_c  <= 1
 0  <= BU_LMPD_176_pe161_c  <= 1
 -1000  <= DG_LMPD_177_pe161_c  <= 1000
 -1.887890e+02  <= DGo_LMPD_177_pe161_c  <= -1.237067e+02
 0  <= FU_LMPD_177_pe161_c  <= 1
 0  <= BU_LMPD_177_pe161_c  <= 1
 -1000  <= DG_LMPD_178_pe161_c  <= 1000
 -1.849187e+02  <= DGo_LMPD_178_pe161_c  <= -1.198364e+02
 0  <= FU_LMPD_178_pe161_c  <= 1
 0  <= BU_LMPD_178_pe161_c  <= 1
 -1000  <= DG_LMPD_179_pe161_c  <= 1000
 -1.856927e+02  <= DGo_LMPD_179_pe161_c  <= -1.206105e+02
 0  <= FU_LMPD_179_pe161_c  <= 1
 0  <= BU_LMPD_179_pe161_c  <= 1
 -1000  <= DG_LMPD_180_pe181_c  <= 1000
 -2.206478e+02  <= DGo_LMPD_180_pe181_c  <= -1.459277e+02
 0  <= FU_LMPD_180_pe181_c  <= 1
 0  <= BU_LMPD_180_pe181_c  <= 1
 -1000  <= DG_LMPD_181_pe181_c  <= 1000
 -2.214219e+02  <= DGo_LMPD_181_pe181_c  <= -1.467018e+02
 0  <= FU_LMPD_181_pe181_c  <= 1
 0  <= BU_LMPD_181_pe181_c  <= 1
 -1000  <= DG_LMPD_182_pe181_c  <= 1000
 -2.175516e+02  <= DGo_LMPD_182_pe181_c  <= -1.428315e+02
 0  <= FU_LMPD_182_pe181_c  <= 1
 0  <= BU_LMPD_182_pe181_c  <= 1
 -1000  <= DG_LMPD_183_pe181_c  <= 1000
 -2.183256e+02  <= DGo_LMPD_183_pe181_c  <= -1.436055e+02
 0  <= FU_LMPD_183_pe181_c  <= 1
 0  <= BU_LMPD_183_pe181_c  <= 1
 -1000  <= DG_LMPD_184_pe181_c  <= 1000
 -2.190997e+02  <= DGo_LMPD_184_pe181_c  <= -1.443796e+02
 0  <= FU_LMPD_184_pe181_c  <= 1
 0  <= BU_LMPD_184_pe181_c  <= 1
 -1000  <= DG_LMPD_185_pe181_c  <= 1000
 -2.198738e+02  <= DGo_LMPD_185_pe181_c  <= -1.451537e+02
 0  <= FU_LMPD_185_pe181_c  <= 1
 0  <= BU_LMPD_185_pe181_c  <= 1
 -1000  <= DG_LMPD_186_pe181_c  <= 1000
 -2.160034e+02  <= DGo_LMPD_186_pe181_c  <= -1.412834e+02
 0  <= FU_LMPD_186_pe181_c  <= 1
 0  <= BU_LMPD_186_pe181_c  <= 1
 -1000  <= DG_LMPD_187_pe181_c  <= 1000
 -2.167775e+02  <= DGo_LMPD_187_pe181_c  <= -1.420574e+02
 0  <= FU_LMPD_187_pe181_c  <= 1
 0  <= BU_LMPD_187_pe181_c  <= 1
 -1000  <= DG_LMPD_188_pe181_c  <= 1000
 -1.177463e+02  <= DGo_LMPD_188_pe181_c  <= -4.389548e+01
 0  <= FU_LMPD_188_pe181_c  <= 1
 0  <= BU_LMPD_188_pe181_c  <= 1
 -1000  <= DG_LMPD_189_pe181_c  <= 1000
 -1.231648e+02  <= DGo_LMPD_189_pe181_c  <= -4.931394e+01
 0  <= FU_LMPD_189_pe181_c  <= 1
 0  <= BU_LMPD_189_pe181_c  <= 1
 -1000  <= DG_LMPD_190_pg160_c  <= 1000
 -2.311854e+02  <= DGo_LMPD_190_pg160_c  <= -1.609372e+02
 0  <= FU_LMPD_190_pg160_c  <= 1
 0  <= BU_LMPD_190_pg160_c  <= 1
 -1000  <= DG_LMPD_191_pg160_c  <= 1000
 -2.319595e+02  <= DGo_LMPD_191_pg160_c  <= -1.617113e+02
 0  <= FU_LMPD_191_pg160_c  <= 1
 0  <= BU_LMPD_191_pg160_c  <= 1
 -1000  <= DG_LMPD_192_pg160_c  <= 1000
 -2.327336e+02  <= DGo_LMPD_192_pg160_c  <= -1.624853e+02
 0  <= FU_LMPD_192_pg160_c  <= 1
 0  <= BU_LMPD_192_pg160_c  <= 1
 -1000  <= DG_LMPD_193_pg160_c  <= 1000
 -2.280892e+02  <= DGo_LMPD_193_pg160_c  <= -1.578409e+02
 0  <= FU_LMPD_193_pg160_c  <= 1
 0  <= BU_LMPD_193_pg160_c  <= 1
 -1000  <= DG_LMPD_194_pg160_c  <= 1000
 -2.273151e+02  <= DGo_LMPD_194_pg160_c  <= -1.570669e+02
 0  <= FU_LMPD_194_pg160_c  <= 1
 0  <= BU_LMPD_194_pg160_c  <= 1
 -1000  <= DG_LMPD_195_pg160_c  <= 1000
 -2.288632e+02  <= DGo_LMPD_195_pg160_c  <= -1.586150e+02
 0  <= FU_LMPD_195_pg160_c  <= 1
 0  <= BU_LMPD_195_pg160_c  <= 1
 -1000  <= DG_LMPD_196_pg160_c  <= 1000
 -2.296373e+02  <= DGo_LMPD_196_pg160_c  <= -1.593891e+02
 0  <= FU_LMPD_196_pg160_c  <= 1
 0  <= BU_LMPD_196_pg160_c  <= 1
 -1000  <= DG_LMPD_197_pg160_c  <= 1000
 -2.304114e+02  <= DGo_LMPD_197_pg160_c  <= -1.601631e+02
 0  <= FU_LMPD_197_pg160_c  <= 1
 0  <= BU_LMPD_197_pg160_c  <= 1
 -1000  <= DG_LMPD_198_pg160_c  <= 1000
 -1.438534e+02  <= DGo_LMPD_198_pg160_c  <= -7.438599e+01
 0  <= FU_LMPD_198_pg160_c  <= 1
 0  <= BU_LMPD_198_pg160_c  <= 1
 -1000  <= DG_LMPD_199_pg160_c  <= 1000
 -1.492718e+02  <= DGo_LMPD_199_pg160_c  <= -7.980445e+01
 0  <= FU_LMPD_199_pg160_c  <= 1
 0  <= BU_LMPD_199_pg160_c  <= 1
 -1000  <= DG_LMPD_200_pg161_c  <= 1000
 -1.969634e+02  <= DGo_LMPD_200_pg161_c  <= -1.301837e+02
 0  <= FU_LMPD_200_pg161_c  <= 1
 0  <= BU_LMPD_200_pg161_c  <= 1
 -1000  <= DG_LMPD_201_pg161_c  <= 1000
 -1.977374e+02  <= DGo_LMPD_201_pg161_c  <= -1.309578e+02
 0  <= FU_LMPD_201_pg161_c  <= 1
 0  <= BU_LMPD_201_pg161_c  <= 1
 -1000  <= DG_LMPD_202_pg161_c  <= 1000
 -1.938671e+02  <= DGo_LMPD_202_pg161_c  <= -1.270875e+02
 0  <= FU_LMPD_202_pg161_c  <= 1
 0  <= BU_LMPD_202_pg161_c  <= 1
 -1000  <= DG_LMPD_203_pg161_c  <= 1000
 -1.930930e+02  <= DGo_LMPD_203_pg161_c  <= -1.263134e+02
 0  <= FU_LMPD_203_pg161_c  <= 1
 0  <= BU_LMPD_203_pg161_c  <= 1
 -1000  <= DG_LMPD_204_pg161_c  <= 1000
 -1.946412e+02  <= DGo_LMPD_204_pg161_c  <= -1.278615e+02
 0  <= FU_LMPD_204_pg161_c  <= 1
 0  <= BU_LMPD_204_pg161_c  <= 1
 -1000  <= DG_LMPD_205_pg161_c  <= 1000
 -1.954152e+02  <= DGo_LMPD_205_pg161_c  <= -1.286356e+02
 0  <= FU_LMPD_205_pg161_c  <= 1
 0  <= BU_LMPD_205_pg161_c  <= 1
 -1000  <= DG_LMPD_206_pg161_c  <= 1000
 -1.961893e+02  <= DGo_LMPD_206_pg161_c  <= -1.294097e+02
 0  <= FU_LMPD_206_pg161_c  <= 1
 0  <= BU_LMPD_206_pg161_c  <= 1
 -1000  <= DG_LMPD_207_pg181_c  <= 1000
 -2.383097e+02  <= DGo_LMPD_207_pg181_c  <= -1.618755e+02
 0  <= FU_LMPD_207_pg181_c  <= 1
 0  <= BU_LMPD_207_pg181_c  <= 1
 -1000  <= DG_LMPD_208_pg181_c  <= 1000
 -2.390838e+02  <= DGo_LMPD_208_pg181_c  <= -1.626496e+02
 0  <= FU_LMPD_208_pg181_c  <= 1
 0  <= BU_LMPD_208_pg181_c  <= 1
 -1000  <= DG_LMPD_209_pg181_c  <= 1000
 -2.398579e+02  <= DGo_LMPD_209_pg181_c  <= -1.634236e+02
 0  <= FU_LMPD_209_pg181_c  <= 1
 0  <= BU_LMPD_209_pg181_c  <= 1
 -1000  <= DG_LMPD_210_pg181_c  <= 1000
 -2.352135e+02  <= DGo_LMPD_210_pg181_c  <= -1.587793e+02
 0  <= FU_LMPD_210_pg181_c  <= 1
 0  <= BU_LMPD_210_pg181_c  <= 1
 -1000  <= DG_LMPD_211_pg181_c  <= 1000
 -2.344394e+02  <= DGo_LMPD_211_pg181_c  <= -1.580052e+02
 0  <= FU_LMPD_211_pg181_c  <= 1
 0  <= BU_LMPD_211_pg181_c  <= 1
 -1000  <= DG_LMPD_212_pg181_c  <= 1000
 -2.359875e+02  <= DGo_LMPD_212_pg181_c  <= -1.595533e+02
 0  <= FU_LMPD_212_pg181_c  <= 1
 0  <= BU_LMPD_212_pg181_c  <= 1
 -1000  <= DG_LMPD_213_pg181_c  <= 1000
 -2.367616e+02  <= DGo_LMPD_213_pg181_c  <= -1.603274e+02
 0  <= FU_LMPD_213_pg181_c  <= 1
 0  <= BU_LMPD_213_pg181_c  <= 1
 -1000  <= DG_LMPD_214_pg181_c  <= 1000
 -2.375357e+02  <= DGo_LMPD_214_pg181_c  <= -1.611015e+02
 0  <= FU_LMPD_214_pg181_c  <= 1
 0  <= BU_LMPD_214_pg181_c  <= 1
 -1000  <= DG_LMPD_215_pg181_c  <= 1000
 -1.361921e+02  <= DGo_LMPD_215_pg181_c  <= -6.060747e+01
 0  <= FU_LMPD_215_pg181_c  <= 1
 0  <= BU_LMPD_215_pg181_c  <= 1
 -1000  <= DG_LMPD_216_pg181_c  <= 1000
 -1.416106e+02  <= DGo_LMPD_216_pg181_c  <= -6.602592e+01
 0  <= FU_LMPD_216_pg181_c  <= 1
 0  <= BU_LMPD_216_pg181_c  <= 1
 -1000  <= DG_LMPD_217_phe-L_c  <= 1000
 -9.558149e+01  <= DGo_LMPD_217_phe-L_c  <= -8.994535e+01
 0  <= FU_LMPD_217_phe-L_c  <= 1
 0  <= BU_LMPD_217_phe-L_c  <= 1
 0  <= FU_LMPD_218_pheme_c  <= 1
 0  <= BU_LMPD_218_pheme_c  <= 1
 0  <= FU_LMPD_219_pheme_c  <= 1
 0  <= BU_LMPD_219_pheme_c  <= 1
 -1000  <= DG_LMPD_220_pro-L_c  <= 1000
 -4.404410e+00  <= DGo_LMPD_220_pro-L_c  <= 1.311060e-01
 0  <= FU_LMPD_220_pro-L_c  <= 1
 0  <= BU_LMPD_220_pro-L_c  <= 1
 -1000  <= DG_LMPD_221_ptrc_c  <= 1000
 4.538380e+00  <= DGo_LMPD_221_ptrc_c  <= 8.275137e+00
 0  <= FU_LMPD_221_ptrc_c  <= 1
 0  <= BU_LMPD_221_ptrc_c  <= 1
 -1000  <= DG_LMPD_222_pydx5p_c  <= 1000
 -7.331219e+01  <= DGo_LMPD_222_pydx5p_c  <= -6.586296e+01
 0  <= FU_LMPD_222_pydx5p_c  <= 1
 0  <= BU_LMPD_222_pydx5p_c  <= 1
 -1000  <= DG_LMPD_223_pydx5p_c  <= 1000
 -4.388894e+01  <= DGo_LMPD_223_pydx5p_c  <= -3.480498e+01
 0  <= FU_LMPD_223_pydx5p_c  <= 1
 0  <= BU_LMPD_223_pydx5p_c  <= 1
 0  <= FU_LMPD_224_q8h2_c  <= 1
 0  <= BU_LMPD_224_q8h2_c  <= 1
 0  <= FU_LMPD_225_q8h2_c  <= 1
 0  <= BU_LMPD_225_q8h2_c  <= 1
 0  <= FU_LMPD_226_q8h2_c  <= 1
 0  <= BU_LMPD_226_q8h2_c  <= 1
 0  <= FU_LMPD_227_q8h2_c  <= 1
 0  <= BU_LMPD_227_q8h2_c  <= 1
 0  <= FU_LMPD_228_q8h2_c  <= 1
 0  <= BU_LMPD_228_q8h2_c  <= 1
 0  <= FU_LMPD_229_q8h2_c  <= 1
 0  <= BU_LMPD_229_q8h2_c  <= 1
 0  <= FU_LMPD_230_q8h2_c  <= 1
 0  <= BU_LMPD_230_q8h2_c  <= 1
 0  <= FU_LMPD_231_q8h2_c  <= 1
 0  <= BU_LMPD_231_q8h2_c  <= 1
 0  <= FU_LMPD_232_q8h2_c  <= 1
 0  <= BU_LMPD_232_q8h2_c  <= 1
 0  <= FU_LMPD_233_q8h2_c  <= 1
 0  <= BU_LMPD_233_q8h2_c  <= 1
 0  <= FU_LMPD_234_q8h2_c  <= 1
 0  <= BU_LMPD_234_q8h2_c  <= 1
 -1000  <= DG_LMPD_235_ribflv_c  <= 1000
 -1.396057e+02  <= DGo_LMPD_235_ribflv_c  <= -1.203148e+02
 0  <= FU_LMPD_235_ribflv_c  <= 1
 0  <= BU_LMPD_235_ribflv_c  <= 1
 -1000  <= DG_LMPD_236_ribflv_c  <= 1000
 -1.391762e+02  <= DGo_LMPD_236_ribflv_c  <= -1.198852e+02
 0  <= FU_LMPD_236_ribflv_c  <= 1
 0  <= BU_LMPD_236_ribflv_c  <= 1
 -1000  <= DG_LMPD_237_ribflv_c  <= 1000
 -1.461094e+02  <= DGo_LMPD_237_ribflv_c  <= -1.288283e+02
 0  <= FU_LMPD_237_ribflv_c  <= 1
 0  <= BU_LMPD_237_ribflv_c  <= 1
 -1000  <= DG_LMPD_238_ribflv_c  <= 1000
 -1.456799e+02  <= DGo_LMPD_238_ribflv_c  <= -1.283987e+02
 0  <= FU_LMPD_238_ribflv_c  <= 1
 0  <= BU_LMPD_238_ribflv_c  <= 1
 -1000  <= DG_LMPD_239_ser-L_c  <= 1000
 5.960161e+00  <= DGo_LMPD_239_ser-L_c  <= 8.574453e+00
 0  <= FU_LMPD_239_ser-L_c  <= 1
 0  <= BU_LMPD_239_ser-L_c  <= 1
 0  <= FU_LMPD_240_sheme_c  <= 1
 0  <= BU_LMPD_240_sheme_c  <= 1
 -1000  <= DG_LMPD_241_spmd_c  <= 1000
 -8.496048e+01  <= DGo_LMPD_241_spmd_c  <= -6.947194e+01
 0  <= FU_LMPD_241_spmd_c  <= 1
 0  <= BU_LMPD_241_spmd_c  <= 1
 -1000  <= DG_LMPD_242_spmd_c  <= 1000
 -8.453094e+01  <= DGo_LMPD_242_spmd_c  <= -6.904240e+01
 0  <= FU_LMPD_242_spmd_c  <= 1
 0  <= BU_LMPD_242_spmd_c  <= 1
 -1000  <= DG_LMPD_243_thf_c  <= 1000
 -2.165566e+02  <= DGo_LMPD_243_thf_c  <= -1.925235e+02
 0  <= FU_LMPD_243_thf_c  <= 1
 0  <= BU_LMPD_243_thf_c  <= 1
 -1000  <= DG_LMPD_244_thf_c  <= 1000
 -2.161270e+02  <= DGo_LMPD_244_thf_c  <= -1.920939e+02
 0  <= FU_LMPD_244_thf_c  <= 1
 0  <= BU_LMPD_244_thf_c  <= 1
 -1000  <= DG_LMPD_245_thf_c  <= 1000
 -2.232377e+02  <= DGo_LMPD_245_thf_c  <= -2.008595e+02
 0  <= FU_LMPD_245_thf_c  <= 1
 0  <= BU_LMPD_245_thf_c  <= 1
 -1000  <= DG_LMPD_246_thf_c  <= 1000
 -2.228082e+02  <= DGo_LMPD_246_thf_c  <= -2.004300e+02
 0  <= FU_LMPD_246_thf_c  <= 1
 0  <= BU_LMPD_246_thf_c  <= 1
 -1000  <= DG_LMPD_247_thmpp_c  <= 1000
 -3.423789e+02  <= DGo_LMPD_247_thmpp_c  <= -3.228940e+02
 0  <= FU_LMPD_247_thmpp_c  <= 1
 0  <= BU_LMPD_247_thmpp_c  <= 1
 -1000  <= DG_LMPD_248_thmpp_c  <= 1000
 -3.419493e+02  <= DGo_LMPD_248_thmpp_c  <= -3.224645e+02
 0  <= FU_LMPD_248_thmpp_c  <= 1
 0  <= BU_LMPD_248_thmpp_c  <= 1
 -1000  <= DG_LMPD_249_thr-L_c  <= 1000
 -6.324230e+00  <= DGo_LMPD_249_thr-L_c  <= -1.525505e+00
 0  <= FU_LMPD_249_thr-L_c  <= 1
 0  <= BU_LMPD_249_thr-L_c  <= 1
 -1000  <= DG_LMPD_250_trp-L_c  <= 1000
 -1.369855e+02  <= DGo_LMPD_250_trp-L_c  <= -1.289287e+02
 0  <= FU_LMPD_250_trp-L_c  <= 1
 0  <= BU_LMPD_250_trp-L_c  <= 1
 -1000  <= DG_LMPD_251_tyr-L_c  <= 1000
 -8.386664e+01  <= DGo_LMPD_251_tyr-L_c  <= -7.870437e+01
 0  <= FU_LMPD_251_tyr-L_c  <= 1
 0  <= BU_LMPD_251_tyr-L_c  <= 1
 0  <= FU_LMPD_252_udcpdp_c  <= 1
 0  <= BU_LMPD_252_udcpdp_c  <= 1
 -1000  <= DG_LMPD_253_utp_c  <= 1000
 -1.241570e+01  <= DGo_LMPD_253_utp_c  <= -8.871699e+00
 0  <= FU_LMPD_253_utp_c  <= 1
 0  <= BU_LMPD_253_utp_c  <= 1
 -1000  <= DG_LMPD_254_val-L_c  <= 1000
 -1.461034e+01  <= DGo_LMPD_254_val-L_c  <= -1.200598e+01
 0  <= FU_LMPD_254_val-L_c  <= 1
 0  <= BU_LMPD_254_val-L_c  <= 1
Binaries
 FU_DM_4CRSOL	 BU_DM_4CRSOL	 FU_DM_5DRIB	 BU_DM_5DRIB	 FU_DM_AMOB	 BU_DM_AMOB	
 FU_DM_MTHTHF	 BU_DM_MTHTHF	 FU_Ec_biomass_iJO1366_WT_53p95M	 BU_Ec_biomass_iJO1366_WT_53p95M	 FU_DM_5mtr_e	 BU_DM_5mtr_e	
 FU_DM_ac_e	 BU_DM_ac_e	 FU_DM_acald_e	 BU_DM_acald_e	 FU_DM_acser_e	 BU_DM_acser_e	
 FU_DM_akg_e	 BU_DM_akg_e	 FU_DM_ala-D_e	 BU_DM_ala-D_e	 FU_DM_ala-L_e	 BU_DM_ala-L_e	
 FU_DM_alaala_e	 BU_DM_alaala_e	 FU_DM_arg-L_e	 BU_DM_arg-L_e	 FU_DM_asn-L_e	 BU_DM_asn-L_e	
 FU_DM_asp-L_e	 BU_DM_asp-L_e	 FU_DM_ca2_e	 BU_DM_ca2_e	 FU_DM_cbl1_e	 BU_DM_cbl1_e	
 FU_DM_cit_e	 BU_DM_cit_e	 FU_DM_cl_e	 BU_DM_cl_e	 FU_DM_co2_e	 BU_DM_co2_e	
 FU_DM_cobalt2_e	 BU_DM_cobalt2_e	 FU_DM_colipa_e	 BU_DM_colipa_e	 FU_DM_cu2_e	 BU_DM_cu2_e	
 FU_DM_cys-L_e	 BU_DM_cys-L_e	 FU_DM_dha_e	 BU_DM_dha_e	 FU_DM_enter_e	 BU_DM_enter_e	
 FU_DM_etoh_e	 BU_DM_etoh_e	 FU_DM_fe2_e	 BU_DM_fe2_e	 FU_DM_fe3_e	 BU_DM_fe3_e	
 FU_DM_for_e	 BU_DM_for_e	 FU_DM_glc_e	 BU_DM_glc_e	 FU_DM_glcn_e	 BU_DM_glcn_e	
 FU_DM_glu-L_e	 BU_DM_glu-L_e	 FU_DM_gly_e	 BU_DM_gly_e	 FU_DM_glyc_e	 BU_DM_glyc_e	
 FU_DM_glyc-R_e	 BU_DM_glyc-R_e	 FU_DM_glycogenn1_c	 BU_DM_glycogenn1_c	 FU_DM_glyclt_e	 BU_DM_glyclt_e	
 FU_DM_gthrd_e	 BU_DM_gthrd_e	 FU_DM_h_e	 BU_DM_h_e	 FU_DM_h2o_e	 BU_DM_h2o_e	
 FU_DM_his-L_e	 BU_DM_his-L_e	 FU_DM_ile-L_e	 BU_DM_ile-L_e	 FU_DM_k_e	 BU_DM_k_e	
 FU_DM_lac-D_e	 BU_DM_lac-D_e	 FU_DM_leu-L_e	 BU_DM_leu-L_e	 FU_DM_lys-L_e	 BU_DM_lys-L_e	
 FU_DM_mal-L_e	 BU_DM_mal-L_e	 FU_DM_meoh_e	 BU_DM_meoh_e	 FU_DM_mg2_e	 BU_DM_mg2_e	
 FU_DM_mn2_e	 BU_DM_mn2_e	 FU_DM_mobd_e	 BU_DM_mobd_e	 FU_DM_nh4_e	 BU_DM_nh4_e	
 FU_DM_ni2_e	 BU_DM_ni2_e	 FU_DM_o2_e	 BU_DM_o2_e	 FU_DM_phe-L_e	 BU_DM_phe-L_e	
 FU_DM_pheme_e	 BU_DM_pheme_e	 FU_DM_pi_e	 BU_DM_pi_e	 FU_DM_pro-L_e	 BU_DM_pro-L_e	
 FU_DM_pyr_e	 BU_DM_pyr_e	 FU_DM_ser-L_e	 BU_DM_ser-L_e	 FU_DM_so4_e	 BU_DM_so4_e	
 FU_DM_spmd_e	 BU_DM_spmd_e	 FU_DM_succ_e	 BU_DM_succ_e	 FU_DM_thr-L_e	 BU_DM_thr-L_e	
 FU_DM_trp-L_e	 BU_DM_trp-L_e	 FU_DM_tyr-L_e	 BU_DM_tyr-L_e	 FU_DM_val-L_e	 BU_DM_val-L_e	
 FU_DM_zn2_e	 BU_DM_zn2_e	 FU_5MTRtex	 BU_5MTRtex	 FU_ACALD	 BU_ACALD	
 FU_ACALDtex	 BU_ACALDtex	 FU_ACALDtpp	 BU_ACALDtpp	 FU_ACKr	 BU_ACKr	
 FU_ACONTa	 BU_ACONTa	 FU_ACONTb	 BU_ACONTb	 FU_ACS	 BU_ACS	
 FU_ACSERtex	 BU_ACSERtex	 FU_ACt2rpp	 BU_ACt2rpp	 FU_ACtex	 BU_ACtex	
 FU_ADK1	 BU_ADK1	 FU_AKGDH	 BU_AKGDH	 FU_AKGt2rpp	 BU_AKGt2rpp	
 FU_AKGtex	 BU_AKGtex	 FU_ALAALAtex	 BU_ALAALAtex	 FU_ALAt2rpp	 BU_ALAt2rpp	
 FU_ALAtex	 BU_ALAtex	 FU_ALCD2x	 BU_ALCD2x	 FU_ALDD2x	 BU_ALDD2x	
 FU_ALDD2y	 BU_ALDD2y	 FU_ARGt3pp	 BU_ARGt3pp	 FU_ARGtex	 BU_ARGtex	
 FU_ASNt2rpp	 BU_ASNt2rpp	 FU_ASNtex	 BU_ASNtex	 FU_ASPO3	 BU_ASPO3	
 FU_ASPO4	 BU_ASPO4	 FU_ASPO5	 BU_ASPO5	 FU_ASPO6	 BU_ASPO6	
 FU_ASPT	 BU_ASPT	 FU_ASPTA	 BU_ASPTA	 FU_ASPt2rpp	 BU_ASPt2rpp	
 FU_ASPtex	 BU_ASPtex	 FU_ATPM	 BU_ATPM	 FU_ATPS4rpp	 BU_ATPS4rpp	
 FU_CA2tex	 BU_CA2tex	 FU_CAt6pp	 BU_CAt6pp	 FU_CBL1tonex	 BU_CBL1tonex	
 FU_CITL	 BU_CITL	 FU_CITt3pp	 BU_CITt3pp	 FU_CITt7pp	 BU_CITt7pp	
 FU_CITtex	 BU_CITtex	 FU_CLt3_2pp	 BU_CLt3_2pp	 FU_CLtex	 BU_CLtex	
 FU_CO2tex	 BU_CO2tex	 FU_CO2tpp	 BU_CO2tpp	 FU_COBALT2tex	 BU_COBALT2tex	
 FU_COBALT2tpp	 BU_COBALT2tpp	 FU_COLIPAabcpp	 BU_COLIPAabcpp	 FU_COLIPAabctex	 BU_COLIPAabctex	
 FU_CS	 BU_CS	 FU_CU2tex	 BU_CU2tex	 FU_CU2tpp	 BU_CU2tpp	
 FU_CYSabc2pp	 BU_CYSabc2pp	 FU_CYStex	 BU_CYStex	 FU_CYStpp	 BU_CYStpp	
 FU_CYTBD2pp	 BU_CYTBD2pp	 FU_CYTBDpp	 BU_CYTBDpp	 FU_CYTBO3_4pp	 BU_CYTBO3_4pp	
 FU_D-LACt2pp	 BU_D-LACt2pp	 FU_D-LACtex	 BU_D-LACtex	 FU_DALAtex	 BU_DALAtex	
 FU_DB4PS	 BU_DB4PS	 FU_DHAPT	 BU_DHAPT	 FU_DHAtex	 BU_DHAtex	
 FU_DHAtpp	 BU_DHAtpp	 FU_DHORD2	 BU_DHORD2	 FU_DHORD5	 BU_DHORD5	
 FU_DHORDfum	 BU_DHORDfum	 FU_EDA	 BU_EDA	 FU_EDD	 BU_EDD	
 FU_ENO	 BU_ENO	 FU_ETOHtex	 BU_ETOHtex	 FU_ETOHtrpp	 BU_ETOHtrpp	
 FU_F6PA	 BU_F6PA	 FU_FBA	 BU_FBA	 FU_FBA3	 BU_FBA3	
 FU_FBP	 BU_FBP	 FU_FDH4pp	 BU_FDH4pp	 FU_FDH5pp	 BU_FDH5pp	
 FU_FE2abcpp	 BU_FE2abcpp	 FU_FE2t2pp	 BU_FE2t2pp	 FU_FE2tex	 BU_FE2tex	
 FU_FE2tpp	 BU_FE2tpp	 FU_FE3abcpp	 BU_FE3abcpp	 FU_FE3tex	 BU_FE3tex	
 FU_FEENTERtex	 BU_FEENTERtex	 FU_FEENTERtpp	 BU_FEENTERtpp	 FU_FLDR2	 BU_FLDR2	
 FU_FORtex	 BU_FORtex	 FU_FORtppi	 BU_FORtppi	 FU_FRD3	 BU_FRD3	
 FU_FUM	 BU_FUM	 FU_G3PD2	 BU_G3PD2	 FU_G3PD5	 BU_G3PD5	
 FU_G3PD6	 BU_G3PD6	 FU_G3PD7	 BU_G3PD7	 FU_G3PT	 BU_G3PT	
 FU_G6PDH2r	 BU_G6PDH2r	 FU_G6PP	 BU_G6PP	 FU_GAPD	 BU_GAPD	
 FU_GLBRAN2	 BU_GLBRAN2	 FU_GLCDpp	 BU_GLCDpp	 FU_GLCNtex	 BU_GLCNtex	
 FU_GLCP	 BU_GLCP	 FU_GLCS1	 BU_GLCS1	 FU_GLCabcpp	 BU_GLCabcpp	
 FU_GLCptspp	 BU_GLCptspp	 FU_GLCt2pp	 BU_GLCt2pp	 FU_GLCtex	 BU_GLCtex	
 FU_GLCtexi	 BU_GLCtexi	 FU_GLDBRAN2	 BU_GLDBRAN2	 FU_GLGC	 BU_GLGC	
 FU_GLUDy	 BU_GLUDy	 FU_GLUabcpp	 BU_GLUabcpp	 FU_GLUt2rpp	 BU_GLUt2rpp	
 FU_GLUt4pp	 BU_GLUt4pp	 FU_GLUtex	 BU_GLUtex	 FU_GLXCL	 BU_GLXCL	
 FU_GLYCAt2rpp	 BU_GLYCAt2rpp	 FU_GLYCAtex	 BU_GLYCAtex	 FU_GLYCDx	 BU_GLYCDx	
 FU_GLYCK	 BU_GLYCK	 FU_GLYCK2	 BU_GLYCK2	 FU_GLYCLTDx	 BU_GLYCLTDx	
 FU_GLYCLTDy	 BU_GLYCLTDy	 FU_GLYCLTt2rpp	 BU_GLYCLTt2rpp	 FU_GLYCLTtex	 BU_GLYCLTtex	
 FU_GLYCTO2	 BU_GLYCTO2	 FU_GLYCTO3	 BU_GLYCTO3	 FU_GLYCTO4	 BU_GLYCTO4	
 FU_GLYCtex	 BU_GLYCtex	 FU_GLYCtpp	 BU_GLYCtpp	 FU_GLYK	 BU_GLYK	
 FU_GLYt2rpp	 BU_GLYt2rpp	 FU_GLYtex	 BU_GLYtex	 FU_GND	 BU_GND	
 FU_GTHRDabc2pp	 BU_GTHRDabc2pp	 FU_GTHRDtex	 BU_GTHRDtex	 FU_H2Otex	 BU_H2Otex	
 FU_H2Otpp	 BU_H2Otpp	 FU_HEX1	 BU_HEX1	 FU_HISt2rpp	 BU_HISt2rpp	
 FU_HIStex	 BU_HIStex	 FU_Htex	 BU_Htex	 FU_ICDHyr	 BU_ICDHyr	
 FU_ICL	 BU_ICL	 FU_ILEt2rpp	 BU_ILEt2rpp	 FU_ILEtex	 BU_ILEtex	
 FU_Kabcpp	 BU_Kabcpp	 FU_Kt2pp	 BU_Kt2pp	 FU_Ktex	 BU_Ktex	
 FU_LDH_D	 BU_LDH_D	 FU_LDH_D2	 BU_LDH_D2	 FU_LEUt2rpp	 BU_LEUt2rpp	
 FU_LEUtex	 BU_LEUtex	 FU_LYSt3pp	 BU_LYSt3pp	 FU_LYStex	 BU_LYStex	
 FU_MALS	 BU_MALS	 FU_MALt3pp	 BU_MALt3pp	 FU_MALtex	 BU_MALtex	
 FU_MDH	 BU_MDH	 FU_MDH2	 BU_MDH2	 FU_MDH3	 BU_MDH3	
 FU_ME1	 BU_ME1	 FU_ME2	 BU_ME2	 FU_MEOHtex	 BU_MEOHtex	
 FU_MG2t3_2pp	 BU_MG2t3_2pp	 FU_MG2tex	 BU_MG2tex	 FU_MG2tpp	 BU_MG2tpp	
 FU_MG2uabcpp	 BU_MG2uabcpp	 FU_MN2tpp	 BU_MN2tpp	 FU_MNt2pp	 BU_MNt2pp	
 FU_MNtex	 BU_MNtex	 FU_MOBDabcpp	 BU_MOBDabcpp	 FU_MOBDtex	 BU_MOBDtex	
 FU_MOX	 BU_MOX	 FU_NADH10	 BU_NADH10	 FU_NADH16pp	 BU_NADH16pp	
 FU_NADH17pp	 BU_NADH17pp	 FU_NADH18pp	 BU_NADH18pp	 FU_NADH5	 BU_NADH5	
 FU_NADH9	 BU_NADH9	 FU_NADK	 BU_NADK	 FU_NADPHQR2	 BU_NADPHQR2	
 FU_NADPHQR3	 BU_NADPHQR3	 FU_NADPHQR4	 BU_NADPHQR4	 FU_NADPPPS	 BU_NADPPPS	
 FU_NADTRHD	 BU_NADTRHD	 FU_NAt3_1p5pp	 BU_NAt3_1p5pp	 FU_NAt3_2pp	 BU_NAt3_2pp	
 FU_NAt3pp	 BU_NAt3pp	 FU_NH4tex	 BU_NH4tex	 FU_NH4tpp	 BU_NH4tpp	
 FU_NI2tex	 BU_NI2tex	 FU_NI2tpp	 BU_NI2tpp	 FU_NI2uabcpp	 BU_NI2uabcpp	
 FU_NTP1	 BU_NTP1	 FU_NTPP6	 BU_NTPP6	 FU_O2tex	 BU_O2tex	
 FU_O2tpp	 BU_O2tpp	 FU_OAADC	 BU_OAADC	 FU_PDH	 BU_PDH	
 FU_PE160abcpp	 BU_PE160abcpp	 FU_PE161abcpp	 BU_PE161abcpp	 FU_PE181abcpp	 BU_PE181abcpp	
 FU_PFK	 BU_PFK	 FU_PFK_3	 BU_PFK_3	 FU_PFL	 BU_PFL	
 FU_PG160abcpp	 BU_PG160abcpp	 FU_PG161abcpp	 BU_PG161abcpp	 FU_PG181abcpp	 BU_PG181abcpp	
 FU_PGI	 BU_PGI	 FU_PGK	 BU_PGK	 FU_PGL	 BU_PGL	
 FU_PGM	 BU_PGM	 FU_PGMT	 BU_PGMT	 FU_PHEMEabcpp	 BU_PHEMEabcpp	
 FU_PHEMEtiex	 BU_PHEMEtiex	 FU_PHEt2rpp	 BU_PHEt2rpp	 FU_PHEtex	 BU_PHEtex	
 FU_PIt2rpp	 BU_PIt2rpp	 FU_PItex	 BU_PItex	 FU_PIuabcpp	 BU_PIuabcpp	
 FU_POR5	 BU_POR5	 FU_POX	 BU_POX	 FU_PPA	 BU_PPA	
 FU_PPA2	 BU_PPA2	 FU_PPC	 BU_PPC	 FU_PPCK	 BU_PPCK	
 FU_PPK2r	 BU_PPK2r	 FU_PPKr	 BU_PPKr	 FU_PPPGO	 BU_PPPGO	
 FU_PPPGO3	 BU_PPPGO3	 FU_PPS	 BU_PPS	 FU_PROt2rpp	 BU_PROt2rpp	
 FU_PROtex	 BU_PROtex	 FU_PTAr	 BU_PTAr	 FU_PYK	 BU_PYK	
 FU_PYRt2rpp	 BU_PYRt2rpp	 FU_PYRtex	 BU_PYRtex	 FU_RPE	 BU_RPE	
 FU_RPI	 BU_RPI	 FU_SERt2rpp	 BU_SERt2rpp	 FU_SERtex	 BU_SERtex	
 FU_SO4t2pp	 BU_SO4t2pp	 FU_SO4tex	 BU_SO4tex	 FU_SPMDt3pp	 BU_SPMDt3pp	
 FU_SPMDtex	 BU_SPMDtex	 FU_SUCCt2_2pp	 BU_SUCCt2_2pp	 FU_SUCCt2_3pp	 BU_SUCCt2_3pp	
 FU_SUCCt3pp	 BU_SUCCt3pp	 FU_SUCCtex	 BU_SUCCtex	 FU_SUCDi	 BU_SUCDi	
 FU_SUCOAS	 BU_SUCOAS	 FU_SULabcpp	 BU_SULabcpp	 FU_TALA	 BU_TALA	
 FU_THD2pp	 BU_THD2pp	 FU_THIORDXi	 BU_THIORDXi	 FU_THRt2pp	 BU_THRt2pp	
 FU_THRt2rpp	 BU_THRt2rpp	 FU_THRtex	 BU_THRtex	 FU_TKT1	 BU_TKT1	
 FU_TKT2	 BU_TKT2	 FU_TPI	 BU_TPI	 FU_TRDR	 BU_TRDR	
 FU_TRPt2rpp	 BU_TRPt2rpp	 FU_TRPtex	 BU_TRPtex	 FU_TRSARr	 BU_TRSARr	
 FU_TYRt2rpp	 BU_TYRt2rpp	 FU_TYRtex	 BU_TYRtex	 FU_VALt2rpp	 BU_VALt2rpp	
 FU_VALtex	 BU_VALtex	 FU_ZN2tpp	 BU_ZN2tpp	 FU_ZNabcpp	 BU_ZNabcpp	
 FU_Zn2tex	 BU_Zn2tex	 FU_LMPD_1_10fthf_c	 BU_LMPD_1_10fthf_c	 FU_LMPD_2_10fthf_c	 BU_LMPD_2_10fthf_c	
 FU_LMPD_3_10fthf_c	 BU_LMPD_3_10fthf_c	 FU_LMPD_4_10fthf_c	 BU_LMPD_4_10fthf_c	 FU_LMPD_5_2dmmql8_c	 BU_LMPD_5_2dmmql8_c	
 FU_LMPD_6_2fe2s_c	 BU_LMPD_6_2fe2s_c	 FU_LMPD_7_2fe2s_c	 BU_LMPD_7_2fe2s_c	 FU_LMPD_8_2fe2s_c	 BU_LMPD_8_2fe2s_c	
 FU_LMPD_9_2fe2s_c	 BU_LMPD_9_2fe2s_c	 FU_LMPD_10_2fe2s_c	 BU_LMPD_10_2fe2s_c	 FU_LMPD_11_2fe2s_c	 BU_LMPD_11_2fe2s_c	
 FU_LMPD_12_2fe2s_c	 BU_LMPD_12_2fe2s_c	 FU_LMPD_13_2fe2s_c	 BU_LMPD_13_2fe2s_c	 FU_LMPD_14_4fe4s_c	 BU_LMPD_14_4fe4s_c	
 FU_LMPD_15_4fe4s_c	 BU_LMPD_15_4fe4s_c	 FU_LMPD_16_4fe4s_c	 BU_LMPD_16_4fe4s_c	 FU_LMPD_17_4fe4s_c	 BU_LMPD_17_4fe4s_c	
 FU_LMPD_18_4fe4s_c	 BU_LMPD_18_4fe4s_c	 FU_LMPD_19_4fe4s_c	 BU_LMPD_19_4fe4s_c	 FU_LMPD_20_4fe4s_c	 BU_LMPD_20_4fe4s_c	
 FU_LMPD_21_4fe4s_c	 BU_LMPD_21_4fe4s_c	 FU_LMPD_22_5mthf_c	 BU_LMPD_22_5mthf_c	 FU_LMPD_23_5mthf_c	 BU_LMPD_23_5mthf_c	
 FU_LMPD_24_accoa_c	 BU_LMPD_24_accoa_c	 FU_LMPD_25_accoa_c	 BU_LMPD_25_accoa_c	 FU_LMPD_26_adocbl_c	 BU_LMPD_26_adocbl_c	
 FU_LMPD_27_adocbl_c	 BU_LMPD_27_adocbl_c	 FU_LMPD_28_adocbl_c	 BU_LMPD_28_adocbl_c	 FU_LMPD_29_adocbl_c	 BU_LMPD_29_adocbl_c	
 FU_LMPD_30_ala-L_c	 BU_LMPD_30_ala-L_c	 FU_LMPD_31_amet_c	 BU_LMPD_31_amet_c	 FU_LMPD_32_amet_c	 BU_LMPD_32_amet_c	
 FU_LMPD_33_arg-L_c	 BU_LMPD_33_arg-L_c	 FU_LMPD_34_asn-L_c	 BU_LMPD_34_asn-L_c	 FU_LMPD_35_bmocogdp_c	 BU_LMPD_35_bmocogdp_c	
 FU_LMPD_36_bmocogdp_c	 BU_LMPD_36_bmocogdp_c	 FU_LMPD_37_bmocogdp_c	 BU_LMPD_37_bmocogdp_c	 FU_LMPD_38_bmocogdp_c	 BU_LMPD_38_bmocogdp_c	
 FU_LMPD_39_btn_c	 BU_LMPD_39_btn_c	 FU_LMPD_40_btn_c	 BU_LMPD_40_btn_c	 FU_LMPD_41_btn_c	 BU_LMPD_41_btn_c	
 FU_LMPD_42_btn_c	 BU_LMPD_42_btn_c	 FU_LMPD_43_btn_c	 BU_LMPD_43_btn_c	 FU_LMPD_44_btn_c	 BU_LMPD_44_btn_c	
 FU_LMPD_45_btn_c	 BU_LMPD_45_btn_c	 FU_LMPD_47_btn_c	 BU_LMPD_47_btn_c	 FU_LMPD_48_ca2_c	 BU_LMPD_48_ca2_c	
 FU_LMPD_49_ca2_c	 BU_LMPD_49_ca2_c	 FU_LMPD_50_ca2_c	 BU_LMPD_50_ca2_c	 FU_LMPD_51_chor_c	 BU_LMPD_51_chor_c	
 FU_LMPD_52_clpn160_p	 BU_LMPD_52_clpn160_p	 FU_LMPD_53_clpn160_p	 BU_LMPD_53_clpn160_p	 FU_LMPD_54_clpn160_p	 BU_LMPD_54_clpn160_p	
 FU_LMPD_55_clpn160_p	 BU_LMPD_55_clpn160_p	 FU_LMPD_56_clpn160_p	 BU_LMPD_56_clpn160_p	 FU_LMPD_57_clpn160_p	 BU_LMPD_57_clpn160_p	
 FU_LMPD_58_clpn160_p	 BU_LMPD_58_clpn160_p	 FU_LMPD_59_clpn160_p	 BU_LMPD_59_clpn160_p	 FU_LMPD_60_clpn160_p	 BU_LMPD_60_clpn160_p	
 FU_LMPD_61_clpn160_p	 BU_LMPD_61_clpn160_p	 FU_LMPD_62_clpn161_p	 BU_LMPD_62_clpn161_p	 FU_LMPD_63_clpn161_p	 BU_LMPD_63_clpn161_p	
 FU_LMPD_64_clpn161_p	 BU_LMPD_64_clpn161_p	 FU_LMPD_65_clpn161_p	 BU_LMPD_65_clpn161_p	 FU_LMPD_66_clpn161_p	 BU_LMPD_66_clpn161_p	
 FU_LMPD_67_clpn161_p	 BU_LMPD_67_clpn161_p	 FU_LMPD_68_clpn161_p	 BU_LMPD_68_clpn161_p	 FU_LMPD_69_clpn181_p	 BU_LMPD_69_clpn181_p	
 FU_LMPD_70_clpn181_p	 BU_LMPD_70_clpn181_p	 FU_LMPD_71_clpn181_p	 BU_LMPD_71_clpn181_p	 FU_LMPD_72_clpn181_p	 BU_LMPD_72_clpn181_p	
 FU_LMPD_73_clpn181_p	 BU_LMPD_73_clpn181_p	 FU_LMPD_74_clpn181_p	 BU_LMPD_74_clpn181_p	 FU_LMPD_75_clpn181_p	 BU_LMPD_75_clpn181_p	
 FU_LMPD_76_clpn181_p	 BU_LMPD_76_clpn181_p	 FU_LMPD_77_clpn181_p	 BU_LMPD_77_clpn181_p	 FU_LMPD_78_clpn181_p	 BU_LMPD_78_clpn181_p	
 FU_LMPD_79_colipa_e	 BU_LMPD_79_colipa_e	 FU_LMPD_80_colipa_e	 BU_LMPD_80_colipa_e	 FU_LMPD_81_colipa_e	 BU_LMPD_81_colipa_e	
 FU_LMPD_82_colipa_e	 BU_LMPD_82_colipa_e	 FU_LMPD_83_colipa_e	 BU_LMPD_83_colipa_e	 FU_LMPD_84_colipa_e	 BU_LMPD_84_colipa_e	
 FU_LMPD_85_colipa_e	 BU_LMPD_85_colipa_e	 FU_LMPD_86_colipa_e	 BU_LMPD_86_colipa_e	 FU_LMPD_87_colipa_e	 BU_LMPD_87_colipa_e	
 FU_LMPD_88_colipa_e	 BU_LMPD_88_colipa_e	 FU_LMPD_89_colipa_e	 BU_LMPD_89_colipa_e	 FU_LMPD_90_colipa_e	 BU_LMPD_90_colipa_e	
 FU_LMPD_91_ctp_c	 BU_LMPD_91_ctp_c	 FU_LMPD_92_cys-L_c	 BU_LMPD_92_cys-L_c	 FU_LMPD_93_cys-L_c	 BU_LMPD_93_cys-L_c	
 FU_LMPD_94_datp_c	 BU_LMPD_94_datp_c	 FU_LMPD_95_datp_c	 BU_LMPD_95_datp_c	 FU_LMPD_96_datp_c	 BU_LMPD_96_datp_c	
 FU_LMPD_97_datp_c	 BU_LMPD_97_datp_c	 FU_LMPD_98_dctp_c	 BU_LMPD_98_dctp_c	 FU_LMPD_99_dgtp_c	 BU_LMPD_99_dgtp_c	
 FU_LMPD_100_dgtp_c	 BU_LMPD_100_dgtp_c	 FU_LMPD_101_dgtp_c	 BU_LMPD_101_dgtp_c	 FU_LMPD_102_dgtp_c	 BU_LMPD_102_dgtp_c	
 FU_LMPD_103_dgtp_c	 BU_LMPD_103_dgtp_c	 FU_LMPD_104_dgtp_c	 BU_LMPD_104_dgtp_c	 FU_LMPD_105_dttp_c	 BU_LMPD_105_dttp_c	
 FU_LMPD_106_dttp_c	 BU_LMPD_106_dttp_c	 FU_LMPD_107_enter_c	 BU_LMPD_107_enter_c	 FU_LMPD_108_fad_c	 BU_LMPD_108_fad_c	
 FU_LMPD_109_fad_c	 BU_LMPD_109_fad_c	 FU_LMPD_110_fad_c	 BU_LMPD_110_fad_c	 FU_LMPD_111_fad_c	 BU_LMPD_111_fad_c	
 FU_LMPD_112_GAM_c	 BU_LMPD_112_GAM_c	 FU_LMPD_113_GAM_c	 BU_LMPD_113_GAM_c	 FU_LMPD_114_GAM_c	 BU_LMPD_114_GAM_c	
 FU_LMPD_115_GAM_c	 BU_LMPD_115_GAM_c	 FU_LMPD_116_gln-L_c	 BU_LMPD_116_gln-L_c	 FU_LMPD_117_gly_c	 BU_LMPD_117_gly_c	
 FU_LMPD_118_gthrd_c	 BU_LMPD_118_gthrd_c	 FU_LMPD_119_gthrd_c	 BU_LMPD_119_gthrd_c	 FU_LMPD_120_gtp_c	 BU_LMPD_120_gtp_c	
 FU_LMPD_121_gtp_c	 BU_LMPD_121_gtp_c	 FU_LMPD_122_gtp_c	 BU_LMPD_122_gtp_c	 FU_LMPD_123_gtp_c	 BU_LMPD_123_gtp_c	
 FU_LMPD_124_hemeO_c	 BU_LMPD_124_hemeO_c	 FU_LMPD_125_hemeO_c	 BU_LMPD_125_hemeO_c	 FU_LMPD_126_his-L_c	 BU_LMPD_126_his-L_c	
 FU_LMPD_127_his-L_c	 BU_LMPD_127_his-L_c	 FU_LMPD_128_ile-L_c	 BU_LMPD_128_ile-L_c	 FU_LMPD_129_leu-L_c	 BU_LMPD_129_leu-L_c	
 FU_LMPD_130_lipopb_c	 BU_LMPD_130_lipopb_c	 FU_LMPD_131_lipopb_c	 BU_LMPD_131_lipopb_c	 FU_LMPD_132_lipopb_c	 BU_LMPD_132_lipopb_c	
 FU_LMPD_133_lipopb_c	 BU_LMPD_133_lipopb_c	 FU_LMPD_134_lys-L_c	 BU_LMPD_134_lys-L_c	 FU_LMPD_135_malcoa_c	 BU_LMPD_135_malcoa_c	
 FU_LMPD_136_malcoa_c	 BU_LMPD_136_malcoa_c	 FU_LMPD_137_met-L_c	 BU_LMPD_137_met-L_c	 FU_LMPD_138_met-L_c	 BU_LMPD_138_met-L_c	
 FU_LMPD_139_mlthf_c	 BU_LMPD_139_mlthf_c	 FU_LMPD_140_mlthf_c	 BU_LMPD_140_mlthf_c	 FU_LMPD_141_mococdp_c	 BU_LMPD_141_mococdp_c	
 FU_LMPD_142_mococdp_c	 BU_LMPD_142_mococdp_c	 FU_LMPD_143_mococdp_c	 BU_LMPD_143_mococdp_c	 FU_LMPD_144_mococdp_c	 BU_LMPD_144_mococdp_c	
 FU_LMPD_145_mococdp_c	 BU_LMPD_145_mococdp_c	 FU_LMPD_146_mococdp_c	 BU_LMPD_146_mococdp_c	 FU_LMPD_147_mococdp_c	 BU_LMPD_147_mococdp_c	
 FU_LMPD_148_mococdp_c	 BU_LMPD_148_mococdp_c	 FU_LMPD_149_mocogdp_c	 BU_LMPD_149_mocogdp_c	 FU_LMPD_150_mocogdp_c	 BU_LMPD_150_mocogdp_c	
 FU_LMPD_151_mocogdp_c	 BU_LMPD_151_mocogdp_c	 FU_LMPD_152_mocogdp_c	 BU_LMPD_152_mocogdp_c	 FU_LMPD_153_mql8_c	 BU_LMPD_153_mql8_c	
 FU_LMPD_154_murein3p3p_p	 BU_LMPD_154_murein3p3p_p	 FU_LMPD_155_murein3px4p_p	 BU_LMPD_155_murein3px4p_p	 FU_LMPD_156_murein4p4p_p	 BU_LMPD_156_murein4p4p_p	
 FU_LMPD_157_murein4px4p_p	 BU_LMPD_157_murein4px4p_p	 FU_LMPD_158_murein4px4px4p_p	 BU_LMPD_158_murein4px4px4p_p	 FU_LMPD_159_nad_c	 BU_LMPD_159_nad_c	
 FU_LMPD_160_nad_c	 BU_LMPD_160_nad_c	 FU_LMPD_161_nad_c	 BU_LMPD_161_nad_c	 FU_LMPD_162_nad_c	 BU_LMPD_162_nad_c	
 FU_LMPD_163_pe160_c	 BU_LMPD_163_pe160_c	 FU_LMPD_164_pe160_c	 BU_LMPD_164_pe160_c	 FU_LMPD_165_pe160_c	 BU_LMPD_165_pe160_c	
 FU_LMPD_166_pe160_c	 BU_LMPD_166_pe160_c	 FU_LMPD_167_pe160_c	 BU_LMPD_167_pe160_c	 FU_LMPD_168_pe160_c	 BU_LMPD_168_pe160_c	
 FU_LMPD_169_pe160_c	 BU_LMPD_169_pe160_c	 FU_LMPD_170_pe160_c	 BU_LMPD_170_pe160_c	 FU_LMPD_171_pe160_c	 BU_LMPD_171_pe160_c	
 FU_LMPD_172_pe160_c	 BU_LMPD_172_pe160_c	 FU_LMPD_173_pe161_c	 BU_LMPD_173_pe161_c	 FU_LMPD_174_pe161_c	 BU_LMPD_174_pe161_c	
 FU_LMPD_175_pe161_c	 BU_LMPD_175_pe161_c	 FU_LMPD_176_pe161_c	 BU_LMPD_176_pe161_c	 FU_LMPD_177_pe161_c	 BU_LMPD_177_pe161_c	
 FU_LMPD_178_pe161_c	 BU_LMPD_178_pe161_c	 FU_LMPD_179_pe161_c	 BU_LMPD_179_pe161_c	 FU_LMPD_180_pe181_c	 BU_LMPD_180_pe181_c	
 FU_LMPD_181_pe181_c	 BU_LMPD_181_pe181_c	 FU_LMPD_182_pe181_c	 BU_LMPD_182_pe181_c	 FU_LMPD_183_pe181_c	 BU_LMPD_183_pe181_c	
 FU_LMPD_184_pe181_c	 BU_LMPD_184_pe181_c	 FU_LMPD_185_pe181_c	 BU_LMPD_185_pe181_c	 FU_LMPD_186_pe181_c	 BU_LMPD_186_pe181_c	
 FU_LMPD_187_pe181_c	 BU_LMPD_187_pe181_c	 FU_LMPD_188_pe181_c	 BU_LMPD_188_pe181_c	 FU_LMPD_189_pe181_c	 BU_LMPD_189_pe181_c	
 FU_LMPD_190_pg160_c	 BU_LMPD_190_pg160_c	 FU_LMPD_191_pg160_c	 BU_LMPD_191_pg160_c	 FU_LMPD_192_pg160_c	 BU_LMPD_192_pg160_c	
 FU_LMPD_193_pg160_c	 BU_LMPD_193_pg160_c	 FU_LMPD_194_pg160_c	 BU_LMPD_194_pg160_c	 FU_LMPD_195_pg160_c	 BU_LMPD_195_pg160_c	
 FU_LMPD_196_pg160_c	 BU_LMPD_196_pg160_c	 FU_LMPD_197_pg160_c	 BU_LMPD_197_pg160_c	 FU_LMPD_198_pg160_c	 BU_LMPD_198_pg160_c	
 FU_LMPD_199_pg160_c	 BU_LMPD_199_pg160_c	 FU_LMPD_200_pg161_c	 BU_LMPD_200_pg161_c	 FU_LMPD_201_pg161_c	 BU_LMPD_201_pg161_c	
 FU_LMPD_202_pg161_c	 BU_LMPD_202_pg161_c	 FU_LMPD_203_pg161_c	 BU_LMPD_203_pg161_c	 FU_LMPD_204_pg161_c	 BU_LMPD_204_pg161_c	
 FU_LMPD_205_pg161_c	 BU_LMPD_205_pg161_c	 FU_LMPD_206_pg161_c	 BU_LMPD_206_pg161_c	 FU_LMPD_207_pg181_c	 BU_LMPD_207_pg181_c	
 FU_LMPD_208_pg181_c	 BU_LMPD_208_pg181_c	 FU_LMPD_209_pg181_c	 BU_LMPD_209_pg181_c	 FU_LMPD_210_pg181_c	 BU_LMPD_210_pg181_c	
 FU_LMPD_211_pg181_c	 BU_LMPD_211_pg181_c	 FU_LMPD_212_pg181_c	 BU_LMPD_212_pg181_c	 FU_LMPD_213_pg181_c	 BU_LMPD_213_pg181_c	
 FU_LMPD_214_pg181_c	 BU_LMPD_214_pg181_c	 FU_LMPD_215_pg181_c	 BU_LMPD_215_pg181_c	 FU_LMPD_216_pg181_c	 BU_LMPD_216_pg181_c	
 FU_LMPD_217_phe-L_c	 BU_LMPD_217_phe-L_c	 FU_LMPD_218_pheme_c	 BU_LMPD_218_pheme_c	 FU_LMPD_219_pheme_c	 BU_LMPD_219_pheme_c	
 FU_LMPD_220_pro-L_c	 BU_LMPD_220_pro-L_c	 FU_LMPD_221_ptrc_c	 BU_LMPD_221_ptrc_c	 FU_LMPD_222_pydx5p_c	 BU_LMPD_222_pydx5p_c	
 FU_LMPD_223_pydx5p_c	 BU_LMPD_223_pydx5p_c	 FU_LMPD_224_q8h2_c	 BU_LMPD_224_q8h2_c	 FU_LMPD_225_q8h2_c	 BU_LMPD_225_q8h2_c	
 FU_LMPD_226_q8h2_c	 BU_LMPD_226_q8h2_c	 FU_LMPD_227_q8h2_c	 BU_LMPD_227_q8h2_c	 FU_LMPD_228_q8h2_c	 BU_LMPD_228_q8h2_c	
 FU_LMPD_229_q8h2_c	 BU_LMPD_229_q8h2_c	 FU_LMPD_230_q8h2_c	 BU_LMPD_230_q8h2_c	 FU_LMPD_231_q8h2_c	 BU_LMPD_231_q8h2_c	
 FU_LMPD_232_q8h2_c	 BU_LMPD_232_q8h2_c	 FU_LMPD_233_q8h2_c	 BU_LMPD_233_q8h2_c	 FU_LMPD_234_q8h2_c	 BU_LMPD_234_q8h2_c	
 FU_LMPD_235_ribflv_c	 BU_LMPD_235_ribflv_c	 FU_LMPD_236_ribflv_c	 BU_LMPD_236_ribflv_c	 FU_LMPD_237_ribflv_c	 BU_LMPD_237_ribflv_c	
 FU_LMPD_238_ribflv_c	 BU_LMPD_238_ribflv_c	 FU_LMPD_239_ser-L_c	 BU_LMPD_239_ser-L_c	 FU_LMPD_240_sheme_c	 BU_LMPD_240_sheme_c	
 FU_LMPD_241_spmd_c	 BU_LMPD_241_spmd_c	 FU_LMPD_242_spmd_c	 BU_LMPD_242_spmd_c	 FU_LMPD_243_thf_c	 BU_LMPD_243_thf_c	
 FU_LMPD_244_thf_c	 BU_LMPD_244_thf_c	 FU_LMPD_245_thf_c	 BU_LMPD_245_thf_c	 FU_LMPD_246_thf_c	 BU_LMPD_246_thf_c	
 FU_LMPD_247_thmpp_c	 BU_LMPD_247_thmpp_c	 FU_LMPD_248_thmpp_c	 BU_LMPD_248_thmpp_c	 FU_LMPD_249_thr-L_c	 BU_LMPD_249_thr-L_c	
 FU_LMPD_250_trp-L_c	 BU_LMPD_250_trp-L_c	 FU_LMPD_251_tyr-L_c	 BU_LMPD_251_tyr-L_c	 FU_LMPD_252_udcpdp_c	 BU_LMPD_252_udcpdp_c	
 FU_LMPD_253_utp_c	 BU_LMPD_253_utp_c	 FU_LMPD_254_val-L_c	 BU_LMPD_254_val-L_c	
End