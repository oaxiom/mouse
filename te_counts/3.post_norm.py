"""

post EDASeq clean-up and annotation.


"""

import sys, os, glob
from glbase3 import *

user_path = os.path.expanduser("~")
ensg = glload(os.path.join("../mm10/mm10_ensembl_v95_ensg_tes.glb"))

raw_expn = expression(filename="rawtags_gc_normed.tsv", format={"force_tsv": True, "skiplines": 0, "ensg": 0}, expn="column[1:]")
raw_expn.sort_conditions()

arr = raw_expn.mean_replicates(
    ['adipocytes_brown_rp1', 'adipocytes_brown_rp2'],
    ['adipose_rp1', 'adipose_rp2'],
    ['B_resting_rp1', 'B_resting_rp2', 'B_resting_rp3', 'B_resting_rp4'],
    ['basophilic_erythroblast_rp1', 'basophilic_erythroblast_rp2', 'basophilic_erythroblast_rp3'],
    ['BM_monocyte_rp1', 'BM_monocyte_rp2', 'BM_monocyte_rp3'],
    ['BM_macrophages_rp1', 'BM_macrophages_rp2', 'BM_macrophages_rp3', 'BM_macrophages_rp4', ],
    ['CLP_rp1', 'CLP_rp2', 'CLP_rp3', 'CLP_rp4' ],
    ['CMP_rp1', 'CMP_rp2', 'CMP_rp3', 'CMP_rp4' ],
    ['CD4T_cells_rp1', 'CD4T_cells_rp2', 'CD4T_cells_rp3', 'CD4T_cells_rp4'],
    ['CD8T_nkg2dn_rp1', 'CD8T_nkg2dn_rp2', 'CD8T_nkg2dn_rp3', ],
    ['CD8T_nkg2dp_rp1', 'CD8T_nkg2dp_rp2', 'CD8T_nkg2dp_rp3', ],
    ['cortex_neuron_rp1', 'cortex_neuron_rp2', 'cortex_neuron_rp4', 'cortex_neuron_rp5', 'cortex_neuron_rp6', 'cortex_neuron_rp7', 'cortex_neuron_rp8'],
    ['dentate_gyrus_rp1', 'dentate_gyrus_rp2', 'dentate_gyrus_rp3', 'dentate_gyrus_rp4', 'dentate_gyrus_rp5'],
    ['DN1_rp1', 'DN1_rp2'],
    ['DN2a_rp1', 'DN2a_rp2'],
    ['DN2b_rp1', 'DN2b_rp2', 'DN2b_rp3'],
    ['DN3_rp1', 'DN3_rp2'],
    ['DP_rp1', 'DP_rp2', 'DP_rp3', 'DP_rp4'],
    ['eos_im_lm_rp1', 'eos_im_lm_rp2'],
    ['eos_im_lp_rp1', 'eos_im_lp_rp2'],
    ['epidermis_basal_E14_rp1', 'epidermis_basal_E14_rp2'],
    ['epidermis_suprabasal_E14_rp1', 'epidermis_suprabasal_E14_rp2'],
    ['EpiSC_rp1', 'EpiSC_rp2', 'EpiSC_rp3', 'EpiSC_rp4', 'EpiSC_rp6', 'EpiSC_DFN_rp1', 'EpiSC_DFN_rp2', 'EpiSC_KFN_rp1', 'EpiSC_KFN_rp2', 'EpiSC_KFN_rp3', 'EpiSC_KFN_rp4'],
    ['erythrocytes_typeA_rp1', 'erythrocytes_typeA_rp2'],
    ['ESC_2iL_HVn_rp1', 'ESC_2iL_HVn_rp2'],
    ['ESC_2iL_HVp_rp1', 'ESC_2iL_HVp_rp2'],
    ['ESC_V65_rp1', 'ESC_V65_rp2'],
    ['ESC_rp1', 'ESC_rp2', 'ESC_E14_rp1', 'ESC_E14_rp2', 'ESC_E14_rp3', 'ESC_E14_rp4', 'ESC_E14_rp5'],
    ['ESC_groundstate_line1_rp1', 'ESC_groundstate_line2_rp1', 'ESC_groundstate_line3_rp1', 'ESC_groundstate_line4_rp1'],
    ['ESC_naive_2iL_rp1','ESC_naive_2iL_rp2','ESC_naive_2iL_rp3',],
    ['ESC_naive_2iLGO_rp1','ESC_naive_2iLGO_rp2','ESC_naive_2iLGO_rp3',],
    ['ESC_PDL_rp1','ESC_PDL_rp2','ESC_PDL_rp3',],
    ['frontal_cortex_rp1', 'frontal_cortex_rp2', 'frontal_cortex_rp3', 'frontal_cortex_rp4', 'frontal_cortex_rp5', 'frontal_cortex_rp6', 'frontal_cortex_rp7', 'frontal_cortex_rp8', 'frontal_cortex_rp9', 'frontal_cortex_rp10', 'frontal_cortex_rp11', 'frontal_cortex_rp12'],
    ['GMP_rp1', 'GMP_rp2', 'GMP_rp3', 'GMP_rp4'],
    ['heart_rp1', 'heart_rp2', 'heart_rp3', 'heart_rp4'],
    ['hair_follicle_stem_cells_rp1','hair_follicle_stem_cells_rp2',],
    ['hair_follicle_transit_amplifying_cells_rp1', 'hair_follicle_transit_amplifying_cells_rp2'],
    ['hippocampus_neuron_rp1', 'hippocampus_neuron_rp2', 'hippocampus_neuron_rp3', 'hippocampus_neuron_rp4', 'hippocampus_neuron_rp5', 'hippocampus_neuron_rp6', 'hippocampus_neuron_rp7', 'hippocampus_neuron_rp8', 'hippocampus_neuron_rp9', 'hippocampus_neuron_rp10', 'hippocampus_neuron_rp11', 'hippocampus_neuron_rp15'],
    ['hippocampus_tissue_rp1', 'hippocampus_tissue_rp2', 'hippocampus_tissue_rp3'],
    ['HSC_rp1', 'HSC_rp2', 'HSC_rp3', 'HSC_rp4', 'HSC_rp5', 'HSC_rp6', 'HSC_rp7', 'HSC_rp8', 'HSC_rp9', 'HSC_rp10', 'HSC_rp11', 'HSC_rp12', 'HSC_rp13', 'HSC_rp14', 'HSC_rp15', 'HSC_rp16'], # deleted: , 'HSC_rp17', 'HSC_rp18', 'HSC_rp19', 'HSC_rp20'],
    ['HSC_MPP1_rp1', 'HSC_MPP1_rp2', 'HSC_MPP1_rp3',],
    ['HSC_MPP2_rp1', 'HSC_MPP2_rp2', 'HSC_MPP2_rp3',],
    ['HSC_MPP3_rp1', 'HSC_MPP3_rp2', 'HSC_MPP3_rp3',],
    ['HSC_MPP4_rp1', 'HSC_MPP4_rp2', 'HSC_MPP4_rp3',],
    ['ICM_E35_rp1', 'ICM_E35_rp2'],
    ['ILC2_rp1', 'ILC2_rp2'],
    ['ileum_rp1', 'ileum_rp2'],
    ['intestine_rp1', 'intestine_rp2'],
    ['keratinocytes_rp1', 'keratinocytes_rp2'],
    ['kidney_rp1', 'kidney_rp2', 'kidney_rp3', 'kidney_rp4', 'kidney_rp5', 'kidney_rp6'],
    ['large_intestine_rp1', 'large_intestine_rp2', 'large_intestine_rp3', 'large_intestine_rp4', 'large_intestine_rp5', 'large_intestine_rp6', 'large_intestine_rp7', 'large_intestine_rp8', 'large_intestine_rp9', 'large_intestine_rp10', 'large_intestine_rp11', 'large_intestine_rp12', 'large_intestine_rp13', 'large_intestine_rp14', 'large_intestine_rp15', 'large_intestine_rp16', 'large_intestine_rp17', 'large_intestine_rp18', ],
    ['liver_rp1', 'liver_rp2', 'liver_rp3', 'liver_rp4', 'liver_rp5'],
    ['lung_rp1', 'lung_rp2', 'lung_rp3', 'lung_rp4', 'lung_rp5'],
    ['lungFibroblasts_rp1', 'lungFibroblasts_rp2'],
    ['mac_im_lm_rp1', 'mac_im_lm_rp2'],
    ['mac_im_lp_rp1', 'mac_im_lp_rp2'],
    ['mac_ip_lm_rp1', 'mac_ip_lm_rp2'],
    ['mac_ip_lp_rp1', 'mac_ip_lp_rp2'],
    ['mast_cells_cd25n_rp1', 'mast_cells_cd25n_rp2'],
    ['mast_cells_cd25p_rp1', 'mast_cells_cd25p_rp2'],
    ['MEF_rp1', 'MEF_rp2', 'MEF_rp3', 'MEF_rp4'],
    ['morula_rp1', 'morula_rp2', 'morula_rp3'],
    ['mst_im_lm_rp1', 'mst_im_lm_rp2'],
    ['mst_im_lp_rp1', 'mst_im_lp_rp2'],
    ['mst_ip_lm_rp1', 'mst_ip_lm_rp2'],
    ['mst_ip_lp_rp1', 'mst_ip_lp_rp2'],
    ['monocytes_rp1', 'monocytes_rp2', 'monocytes_rp3', 'monocytes_rp4',],
    ['myeloid_rp1', 'myeloid_rp2', 'myeloid_rp3'],
    ['NK_cells_rp1', 'NK_cells_rp2'],
    ['neocortex_cortical_plate_rp1', 'neocortex_cortical_plate_rp2', 'neocortex_cortical_plate_rp3', 'neocortex_cortical_plate_rp4', 'neocortex_cortical_plate_rp5'],
    ['neocortex_subventricular_zone_rp1', 'neocortex_subventricular_zone_rp2', 'neocortex_subventricular_zone_rp3', 'neocortex_subventricular_zone_rp4', 'neocortex_subventricular_zone_rp5'],
    ['neocortex_ventricular_zone_rp1', 'neocortex_ventricular_zone_rp2', 'neocortex_ventricular_zone_rp3', 'neocortex_ventricular_zone_rp4', 'neocortex_ventricular_zone_rp5'],
    ['neu_im_lm_rp1', 'neu_im_lm_rp2'],
    ['neu_im_lp_rp1', 'neu_im_lp_rp2'],
    ['neu_ip_lm_rp1', 'neu_ip_lm_rp2'],
    ['neu_ip_lp_rp1', 'neu_ip_lp_rp2'],
    ['NPC_rp1', 'NPC_rp2'],
    ['nucleus_accumbens_rp1', 'nucleus_accumbens_rp2', 'nucleus_accumbens_rp3'],
    ['olfactory_bulb_rp1', 'olfactory_bulb_rp2', 'olfactory_bulb_rp3', 'olfactorybulb_rp4', 'olfactorybulb_rp5'],
    ['oocytes_rp1', 'oocytes_rp2', 'oocytes_rp3', 'oocyte_rp4', 'oocyte_rp5', 'oocyte_rp6', 'oocyte_rp7'],
    ['orthochromatic_erythroblast_rp1', 'orthochromatic_erythroblast_rp2', 'orthochromatic_erythroblast_rp3'],
    ['polychromatic_erythroblast_rp1', 'polychromatic_erythroblast_rp2', 'polychromatic_erythroblast_rp3'],
    ['preadipocytes_brown_rp1', 'preadipocytes_brown_rp2'],
    #['prim_typeA_spermatogonia_rp1', 'prim_typeA_spermatogonia_rp2'],
    ['proerythroblast_rp1', 'proerythroblast_rp2', 'proerythroblast_rp3'],
    ['pronuclei_rp1', 'pronuclei_rp2', 'pronuclei_rp3'],
    ['prostate_basal_cells_rp1', 'prostate_basal_cells_rp2', 'prostate_basal_cells_rp3', 'prostate_basal_cells_rp4', 'prostate_basal_cells_rp5', 'prostate_basal_cells_rp6'],
    ['rsEpiSC_rp1', 'rsEpiSC_rp2'],
    ['sdc_im_lm_rp1', 'sdc_im_lm_rp2'],
    ['sdc_im_lp_rp1', 'sdc_im_lp_rp2'],
    ['sdc_ip_lm_rp1', 'sdc_ip_lm_rp2'],
    ['sdc_ip_lp_rp1', 'sdc_ip_lp_rp2'],
    ['skeletal_muscle_rp1', 'skeletal_muscle_rp2', 'skeletal_muscle_rp3'],
    ['striatum_rp1', 'striatum_rp2'],
    ['substantia_nigra_rp1', 'substantia_nigra_rp2', 'substantia_nigra_rp3'],
    ['subventricularzone_rp1', 'subventricularzone_rp2', 'subventricularzone_rp3', 'subventricularzone_rp4'],
    #['typeA_spermatogonia_rp1', 'typeA_spermatogonia_rp2'],
    #['typeB_spermatogonia_rp1', 'typeB_spermatogonia_rp2'],
    ['ventral_tagmental_rp1', 'ventral_tagmental_rp2', 'ventral_tagmental_rp3'],
    ['white_matter_glia_rp1', 'white_matter_glia_rp2'],
    ['proB_fracB_rp1', 'proB_fracB_rp2'],
    ['proB_fracCC_rp1', 'proB_fracCC_rp2'],
    ['proB_rp1', 'proB_rp2'],
    ['spermatids_rp1', 'spermatids_rp2', 'spermatids_rp3', 'spermatids_rp4', 'spermatids_rp5'],
    ['spermatocytes_rp1', 'spermatocytes_rp2', 'spermatocytes_rp3', 'spermatocytes_rp4', 'spermatocytes_rp5'],
    ['spinalcord_rp1', 'spinalcord_rp2', 'spinal_cord_rp3', 'spinal_cord_rp4'],
    ['B_germinal_center_rp1', 'B_germinal_center_rp2'],
    ['B_active_rp1', 'B_active_rp2', 'B_active_rp3', 'B_active_rp4'],
    ['BMDM_rp1', 'BMDM_rp2', 'BMDM_rp3'],
    ['B_naive_rp1', 'B_naive_rp2', 'B_naive_rp3', 'B_naive_rp4'],
    ['caecum_rp1', 'caecum_rp2', 'caecum_rp3'],
    ['CD4T_iTreg_rp1', 'CD4T_iTreg_rp2'],
    ['CD4T_Th1_rp1', 'CD4T_Th1_rp2', 'CD4T_Th1_rp3', 'CD4T_Th1_rp4', 'CD4T_Th1_rp5', 'CD4T_Th1_rp6'],
    ['CD4T_Th17_rp1' ,'CD4T_Th17_rp2', 'CD4T_Th17_rp3'],
    ['CD4T_Th2_rp1', 'CD4T_Th2_rp2', 'CD4T_Th2_rp3', 'CD4T_Th2_rp4', 'CD4T_Th2_rp5', 'CD4T_Th2_rp6'],
    ['CD8T_rp1', 'CD8T_cells_rp1', 'CD8T_cells_rp2', 'CD8T_cells_rp3', 'CD8T_cells_rp4', 'CD8T_untre_rp1'],
    ['cerebellar_gran_neurons_rp1', 'cerebellar_gran_neurons_rp2', 'cerebellar_gran_neurons_rp3', 'cerebellar_gran_neurons_rp4'],
    ['cerebellum_rp1', 'cerebellum_rp2', 'cerebellum_rp3', 'cerebellum_rp4'],
    ['dorsal_rootganglia_rp1', 'dorsal_rootganglia_rp2', 'dorsal_rootganglia_rp3', 'dorsal_rootganglia_rp4', 'dorsal_rootganglia_rp5', 'dorsal_rootganglia_rp6'],
    ['E165_cortical_neurons_rp1', 'E165_cortical_neurons_rp2', 'E165_cortical_neurons_rp3'],
    ['E165_skin_rp1', 'E165_skin_rp2'],
    ['erythroblast_rp1', 'erythroblast_rp2', 'erythroblast_rp3', 'erythroblast_rp4'],
    ['eye_rp1', 'eye_rp2', 'eye_rp3', 'eye_rp4', 'eye_rp5', 'eye_rp6'],
    ['megakryocyte_erythroid_progenitor_rp1', 'megakryocyte_erythroid_progenitor_rp2', 'megakryocyte_erythroid_progenitor_rp3', 'megakryocyte_erythroid_progenitor_rp4', 'megakryocyte_erythroid_progenitor_rp5', 'megakryocyte_erythroid_progenitor_rp6'],
    ['motor_neurons_rp1', 'motor_neurons_rp2', 'motor_neurons_rp3', 'motor_neurons_rp4'],
    ['myoblast_rp1', 'myoblast_rp2'],
    ['neonatal_tail_fibroblast_rp1', 'neonatal_tail_fibroblast_rp2'],
    ['pineal_gland_rp1', 'pineal_gland_rp2'],
    ['placenta_rp1', 'placenta_rp2', 'placenta_rp3', 'placenta_rp4', 'placenta_rp5'],
    ['retina_rp1', 'retina_rp2'],
    ['spleen_rp1', 'spleen_rp2', 'spleen_rp3', 'spleen_rp4', 'spleen_rp5'],
    ['splenicB_cells_rp1', 'splenicB_cells_rp2', 'splenicB_cells_rp3', 'splenicB_cells_rp4'],
    ['striated_muscle_rp1', 'striated_muscle_rp2'],
    ['telencephalon_rp1', 'telencephalon_rp2', 'telencephalon_rp3'],
    ['testicular_rp1', 'testicular_rp2', 'testicular_rp3', 'testicular_rp4'],
    ['X2clc_rp1', 'X2clc_rp2'],
    ['X2clc_D1_rp3', 'X2clc_D1_rp4', 'X2clc_D2_rp5', 'X2clc_D2_rp6'], # very different...
    ['X2C_embryo_rp1', 'X2C_embryo_rp2', 'X2C_embryo_rp3', 'embryo_2C_rp1', 'embryo_2C_rp2', 'embryo_2C_rp3'],
    ['X4C_embryo_rp1', 'X4C_embryo_rp2', 'X4C_embryo_rp3'],
    ['X8C_embryo_rp1', 'X8C_embryo_rp2', 'X8C_embryo_rp3'],
    ['MSC_bone_rp1', 'MSC_bone_rp2', 'MSC_bone_rp3', 'MSC_bone_rp4', 'MSC_bone_rp5', 'MSC_bone_rp6', 'MSC_bone_rp7' ,'MSC_bone_rp8', 'MSC_bone_rp9'],
    ['MSC_skin_rp1', 'MSC_skin_rp2', 'MSC_skin_rp3', 'MSC_skin_rp4', 'MSC_skin_rp5', 'MSC_skin_rp6', 'MSC_skin_rp7', 'MSC_skin_rp8', 'MSC_skin_rp9', 'MSC_skin_rp10', 'MSC_skin_rp11', 'MSC_skin_rp12', 'MSC_skin_rp13'],
    ['MSC_thymus_rp1', 'MSC_thymus_rp2', 'MSC_thymus_rp3', 'MSC_thymus_rp4', 'MSC_thymus_rp5', 'MSC_thymus_rp6', 'MSC_thymus_rp7', 'MSC_thymus_rp8'],
    ['patski_rp1', 'patski_rp2', 'patski_rp3', 'patski_rp4', 'patski_rp5', 'patski_rp6', 'patski_rp7'],
    #['genital_fatpad_rp1', 'genital_fatpad_rp2', 'genital_fatpad_rp3', 'genital_fatpad_rp4', 'genital_fatpad_rp5', 'genital_fatpad_rp6'],
    ['NKT_stage1_rp1', 'NKT_stage1_rp2' ,'NKT_stage1_rp3', 'NKT_stage1_rp4', 'NKT_stage1_rp5'],
    ['NKT_stage2_rp1', 'NKT_stage2_rp2', 'NKT_stage2_rp3', 'NKT_stage2_rp4', 'NKT_stage2_rp5', 'NKT_stage2_rp6', 'NKT_stage2_rp7'],
    ['NKT_stage3_rp1', 'NKT_stage3_rp2', 'NKT_stage3_rp3', 'NKT_stage3_rp4', 'NKT_stage3_rp5', 'NKT_stage3_rp6'],
    ['B_cells_CD19p_rp1', 'B_cells_CD19p_rp2', 'B_cells_CD19p_rp3', 'B_cells_CD19p_rp4', 'B_cells_CD19p_rp5', 'B_cells_CD19p_rp6'],
    ['B_cells_CD43m_rp1', 'B_cells_CD43m_rp2', 'B_cells_CD43m_rp3', 'B_cells_CD43m_rp4', 'B_cells_CD43m_rp5', 'B_cells_CD43m_rp6', 'B_cells_CD43m_rp7', 'B_cells_CD43m_rp8', 'B_cells_CD43m_rp9', 'B_cells_CD43m_rp10', 'B_cells_CD43m_rp11', 'B_cells_CD43m_rp12'],
    ['XEN_rp1', 'XEN_rp2', 'XEN_rp3', 'XEN_rp4', 'XEN_rp5', 'XEN_rp6'],
    ['astrocyte_rp1', 'astrocyte_rp2'],
    #['cerebral_cortex_endothelium_rp1', 'cerebral_cortex_endothelium_rp2'],
    ['chondrocyte_rib_rp1', 'chondrocyte_rib_rp2', 'chondrocyte_rib_rp3', 'chondrocyte_rib_rp4'],
    ['CD4T_DPCD3hi_rp1', 'CD4T_DPCD3hi_rp2'],
    ['CD4T_DPCD3lo_rp1', 'CD4T_DPCD3lo_rp2'],
    ['oligodendrocyte_myelinating_rp1', 'oligodendrocyte_myelinating_rp2'],
    ['oligodendrocyte_newly_formed_rp1', 'oligodendrocyte_newly_formed_rp2'],
    ['oligodendrocyte_precursor_rp1', 'oligodendrocyte_precursor_rp2'],
    ['adrenal_medulla_rp1', 'adrenal_medulla_rp2', 'adrenal_medulla_rp3'],
    ['carotid_body_rp1', 'carotid_body_rp2', 'carotid_body_rp3'],
    #['cranial_neural_crest_mandible_rp1', 'cranial_neural_crest_mandible_rp2', 'cranial_neural_crest_mandible_rp3'],
    ['E125_atrioventricular_canal_rp1', 'E125_atrioventricular_canal_rp2', 'E125_atrioventricular_canal_rp3', 'E125_atrioventricular_canal_rp4'],
    ['pancreatic_alpha_rp1', 'pancreatic_alpha_rp2'],
    ['pancreatic_beta_rp1', 'pancreatic_beta_rp2'],
    ['cortical_thymic_epithelial_rp1', 'cortical_thymic_epithelial_cell_rp1', 'cortical_thymic_epithelial_cell_rp2'],
    ['prostate_luminal_cells_rp1', 'prostate_luminal_cells_rp2', 'prostate_luminal_cells_rp3', 'prostate_luminal_cells_rp4', 'prostate_luminal_cells_rp5'],
    ['brain_rp1', 'brain_rp2', 'brain_rp3'],
    ['BM_rp1', 'BM_rp2'],
    ['pgc_male_rp1', 'pgc_male_rp2'],
    ['ectoplacental_cone_rp1', 'ectoplacental_cone_rp2'],
    ['TSC_rp1', 'TSC_rp2'],
    ['medullary_thymic_epithelial_rp1', 'medullary_thymic_epithelial_cell_rp1', 'medullary_thymic_epithelial_cell_rp2'],
    ['medullary_thymic_epithelial_cell_immature_rp1', 'medullary_thymic_epithelial_cell_immature_rp2'],
    ['medullary_thymic_epithelial_cell_mature_rp1', 'medullary_thymic_epithelial_cell_mature_rp2'],
    ['intestinal_stem_cell_rp1', 'intestinal_stem_cell_rp2'],
    ['small_intestine_cells_rp1', 'small_intestine_cells_rp2'],
    ['colon_epithelial_cells_rp1', 'colon_epithelial_cells_rp2', 'colon_epithelial_cells_rp3'],
    #['leptotene_spermatocytes_rp1', 'leptotene_spermatocytes_rp2'],
    #['pachytene_spermatocytes_rp1', 'pachytene_spermatocytes_rp2'],
    ['neuron_rp1', 'neuron_rp2'],
    ['DN1eP_rp1', 'DN1eP_rp2', 'DN1eP_rp3', 'DN1eP_rp4'],
    ['DN1ePNKT_rp1', 'DN1ePNKT_rp2', 'DN1ePNKT_rp3', 'DN1ePNKT_rp4'],
    ['round_spermatids_rp1', 'round_spermatids_rp2'],
    ['elongative_spermatids_rp1', 'elongative_spermatids_rp2'],
    ['fatpad_rp1', 'fatpad_rp2'],
    ['granulocytes_rp1', 'granulocyte_rp2', 'granulocyte_rp3', 'granulocyte_rp4', 'granulocyte_rp5', 'granulocyte_rp6'],
    ['colon_rp1', 'colon_rp2', 'colon_rp3'],
    ['MPP_rp1', 'MPP_rp2'],
    ['cortex_rp1', 'cortex_rp2', 'cortex_rp3', 'cortex_rp4', 'cortex_rp5'],
    ['microglia_rp1', 'microglia_rp2', 'microglia_rp3'],
    # Single cell embryo data:
    ["SS_Zygote_E1_C1", "SS_Zygote_E2_C1", "SS_Zygote_E3_C1" , "SS_Zygote_E4_C1"],
    ["SS_Embryo2C_early_E1_C2", "SS_Embryo2C_early_E2_C1", "SS_Embryo2C_early_E3_C2"],
    ["SS_Embryo2C_late_E6_C2", "SS_Embryo2C_late_E7_C1", "SS_Embryo2C_late_E7_C1", "SS_Embryo2C_late_E7_C2", "SS_Embryo2C_late_E8_C2", "SS_Embryo2C_late_E9_C1", "SS_Embryo2C_late_E9_C2"],
    ["SS_Embryo2C_mid_E5_C2", "SS_Embryo2C_mid_E6_C1", "SS_Embryo2C_mid_E6_C2", "SS_Embryo2C_mid_E7_C1", "SS_Embryo2C_mid_E7_C2"],
    ["SS_Embryo4C_E1_C2", "SS_Embryo4C_E1_C4", "SS_Embryo4C_E2_C1", "SS_Embryo4C_E2_C2", "SS_Embryo4C_E2_C3", "SS_Embryo4C_E2_C4", "SS_Embryo4C_E3_C1", "SS_Embryo4C_E3_C3", "SS_Embryo4C_E3_C4", "SS_Embryo4C_E4_C1", "SS_Embryo4C_E4_C2", "SS_Embryo4C_E4_C3", "SS_Embryo4C_E4_C4"],
    ["SS_Embryo8C_E1_C1", "SS_Embryo8C_E1_C2","SS_Embryo8C_E1_C6","SS_Embryo8C_E1_C7", "SS_Embryo8C_E1_C8", "SS_Embryo8C_E2_C1", "SS_Embryo8C_E2_C2",
    "SS_Embryo8C_E2_C3", "SS_Embryo8C_E2_C4", "SS_Embryo8C_E2_C8","SS_Embryo8C_E5_C1","SS_Embryo8C_E5_C2","SS_Embryo8C_E5_C3","SS_Embryo8C_E5_C6",
    "SS_Embryo8C_E5_C8","SS_Embryo8C_E8_C1","SS_Embryo8C_E8_C2","SS_Embryo8C_E8_C6","SS_Embryo8C_E8_C7","SS_Embryo8C_E8_C8",
    'X8C_embryo_SS_rp1', 'X8C_embryo_SS_rp2', 'X8C_embryo_SS_rp3', 'X8C_embryo_SS_rp4', 'X8C_embryo_SS_rp5',
    'X8C_embryo_SS_rp6', 'X8C_embryo_SS_rp7', 'X8C_embryo_SS_rp8', 'X8C_embryo_SS_rp9', 'X8C_embryo_SS_rp10',
    'X8C_embryo_SS_rp11', 'X8C_embryo_SS_rp12'],
    ["SS_Morula_E1_C10", "SS_Morula_E1_C11", "SS_Morula_E1_C14", "SS_Morula_E1_C15", "SS_Morula_E1_C2", "SS_Morula_E1_C3", "SS_Morula_E1_C4", "SS_Morula_E1_C5", "SS_Morula_E1_C6", "SS_Morula_E1_C7",
    "SS_Morula_E1_C8", "SS_Morula_E1_C9", "SS_Morula_E4_C1", "SS_Morula_E4_C3", "SS_Morula_E4_C4", "SS_Morula_E4_C5", "SS_Morula_E4_C6", "SS_Morula_E4_C7", "SS_Morula_E5_C1", "SS_Morula_E5_C10",
    "SS_Morula_E5_C11", "SS_Morula_E5_C12", "SS_Morula_E5_C13", "SS_Morula_E5_C2", "SS_Morula_E5_C5", "SS_Morula_E5_C6", "SS_Morula_E5_C7", "SS_Morula_E5_C8", "SS_Morula_E5_C9", "SS_Morula_E6_C1",
    "SS_Morula_E6_C12", "SS_Morula_E6_C2", "SS_Morula_E6_C3", "SS_Morula_E6_C4", "SS_Morula_E6_C6", "SS_Morula_E6_C7", "SS_Morula_E6_C8"],
    # blastocyst early TE
    ["SS_Early_blastocyst_E4_C8", "SS_Early_blastocyst_E2_C5", "SS_Early_blastocyst_E2_C4", "SS_Early_blastocyst_E2_C3", "SS_Early_blastocyst_E2_C2", "SS_Early_blastocyst_E3_C13",
    "SS_Early_blastocyst_E2_C7", "SS_Early_blastocyst_E4_C9", "SS_Early_blastocyst_E4_C18", "SS_Early_blastocyst_E4_C13", "SS_Early_blastocyst_E2_C1"],
    # blastocyst mid TE
    ["SS_Mid_blastocyst_E2_C17", "SS_Mid_blastocyst_E3_C7", "SS_Mid_blastocyst_E3_C8", "SS_Mid_blastocyst_E2_C12", "SS_Mid_blastocyst_E3_C6", "SS_Mid_blastocyst_E1_C14", "SS_Mid_blastocyst_E3_C5",
    "SS_Mid_blastocyst_E3_C4", "SS_Mid_blastocyst_E3_C2", "SS_Mid_blastocyst_E1_C15", "SS_Mid_blastocyst_E1_C10", "SS_Mid_blastocyst_E3_C23", "SS_Mid_blastocyst_E3_C18"],
    # blastocyst late TE
    ["SS_Late_blastocyst_E2_C8", "SS_Late_blastocyst_E2_C3", "SS_Late_blastocyst_E1_C26", "SS_Late_blastocyst_E1_C9", "SS_Late_blastocyst_E1_C19", "SS_Late_blastocyst_E1_C16"],
    # blastocyst early ICM
    ["SS_Early_blastocyst_E2_C16", "SS_Early_blastocyst_E3_C6", "SS_Early_blastocyst_E3_C4", "SS_Early_blastocyst_E3_C1", "SS_Early_blastocyst_E3_C3", "SS_Early_blastocyst_E3_C2",
    "SS_Early_blastocyst_E2_C22", "SS_Early_blastocyst_E2_C19", "SS_Early_blastocyst_E4_C16", "SS_Early_blastocyst_E4_C12", "SS_Early_blastocyst_E3_C9", "SS_Early_blastocyst_E4_C17",
    "SS_Early_blastocyst_E4_C14", "SS_Early_blastocyst_E2_C17"],
    # blastocyst mid ICM
    ["SS_Mid_blastocyst_E2_C9", "SS_Mid_blastocyst_E3_C13", "SS_Mid_blastocyst_E1_C20", "SS_Mid_blastocyst_E1_C13", "SS_Mid_blastocyst_E1_C23", "SS_Mid_blastocyst_E2_C6", "SS_Mid_blastocyst_E2_C4", "SS_Mid_blastocyst_E2_C3",
    "SS_Mid_blastocyst_E2_C1", "SS_Mid_blastocyst_E3_C11", "SS_Mid_blastocyst_E2_C2", "SS_Mid_blastocyst_E1_C8", "SS_Mid_blastocyst_E1_C24", "SS_Mid_blastocyst_E1_C11", "SS_Mid_blastocyst_E1_C9",
    "SS_Mid_blastocyst_E1_C12", "SS_Mid_blastocyst_E2_C15", "SS_Mid_blastocyst_E2_C18", "SS_Mid_blastocyst_E2_C16", "SS_Mid_blastocyst_E2_C5", "SS_Mid_blastocyst_E2_C14", "SS_Mid_blastocyst_E2_C24",
    "SS_Mid_blastocyst_E2_C13", "SS_Mid_blastocyst_E2_C23", "SS_Mid_blastocyst_E1_C17", "SS_Mid_blastocyst_E2_C10", "SS_Mid_blastocyst_E1_C19", "SS_Mid_blastocyst_E1_C5"],
    # blastocyst PrE
    ["SS_Late_blastocyst_E2_C2", "SS_Late_blastocyst_E2_C9", "SS_Late_blastocyst_E1_C8", "SS_Late_blastocyst_E1_C7", "SS_Late_blastocyst_E2_C14", "SS_Late_blastocyst_E1_C27",
    "SS_Late_blastocyst_E1_C10", "SS_Late_blastocyst_E1_C24", "SS_Late_blastocyst_E1_C11"],
    output_pears="Pearson_correlations.tsv",
    pearson_hist="Pearson_histogram.png",
    threshold=0.6,
    _ignore_missing_samples=True)

# pretify the condition names here:
pretty = []
for i in arr.getConditionNames():
    newi = i.replace("_rp1", "").replace("_", " ")
    if newi.endswith(" 1"):
        newi = newi.replace(" 1", "")
    pretty.append(newi)
arr.setConditionNames(pretty)

# It is possible to make it here with a mean_replicates genes that actually have zero in all columns.
# I thought that sort of case would be really rare. But unfortunately it's surprisingly common.
arr = arr.filter_low_expressed(10**1.6, 1)

mapped = ensg.map(key="ensg", genelist=arr)
mapped.strip_errs() # I probably want a error key containing version too...
print(list(mapped.keys()))
mapped.saveTSV("genes_ntc_expression.tsv", key_order=['ensg', 'name'])

# Clean out hanger on data
# ['loc', 'name', 'ensg', 'tss_loc', 'conditions', 'strand']
f = {"force_tsv": True, 'ensg': 0, 'name': 1}
mapped = expression(filename='genes_ntc_expression.tsv', format=f, expn='column[2:]')
mapped.save("genes_ntc_expression.glb") # These are the final tables.

mapped = ensg.map(key="ensg", genelist=raw_expn)
mapped.save("genes_ntc_expression_nonmerged.glb")
mapped.saveTSV("genes_ntc_expression_nonmerged.tsv")

arr.tree(filename="tree.png", color_threshold=0.0, label_size=5, size=(5,14))
print()
print("Num conditions =", len(arr.getConditionNames()))
print()
