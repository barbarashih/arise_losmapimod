# plot change in specific GO terms
library(biomaRt)
#library(org.Hs.eg.db)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

GO_interest <- list (
					leukocyte_migration = "GO:0050900",
					T_cell_mediated_cytotoxicity = "GO:0001913",
					type_I_interferon_signaling_pathway = "GO:0060337",
					type_I_interferon_biosynthetic_process ="GO:0045351",
					response_to_interferon_gamma = "GO:0034341",
					T_cell_activation ="GO:0042110",
					regulation_of_wound_healing ="GO:0061041",
					response_to_cytokine_stimulus ="GO:0034097",
					regulation_of_leucocyte_activation = "GO:0002694",
					MAPK_cascade = "GO:0000165",
					leucocyte_activation= "GO:0045321",
					innate_immune_response = "GO:0002226",
					inflammatory_response = "GO:0006954",
					interleukin17_production = "GO:0032620",
					cellular_defense_response = "GO:0006968",
					cytokine_production = "GO:0001816",
					cellular_response_to_stress = "GO:0033554",
					cell_activation = "GO:0001775",
					positive_regulation_of_MAPK_cascade = "GO:0043410",
					stress_activated_MAPK_cascade = "GO:0051403",
					positive_regulation_of_stress_activated_MAPK_cascade = "GO:0032874"
					)
GO_interest_S6 <- list (
					lymphocyte_chemotaxis='GO:0048247',
					monocyte_chemotaxis='GO:0002548',
					chemokine_mediated_signaling_pathway='GO:0070098',
					neutrophil_chemotaxis='GO:0030593',
					positive_regulation_of_cytokine_biosynthetic_process='GO:0042108',
					positive_regulation_of_T_cell_proliferation='GO:0042102',
					cellular_response_to_interleukin_1='GO:0071347',
					cellular_response_to_interferon_gamma='GO:0071346',
					cellular_response_to_tumor_necrosis_factor='GO:0071356',
					positive_regulation_of_ERK1_and_ERK2_cascade='GO:0070374',
					cellular_response_to_lipopolysaccharide='GO:0071222',
					inflammatory_response='GO:0006954',
					negative_regulation_of_cell_population_proliferation='GO:0008285',
					leukocyte_activation='GO:0045321')
		# set of genes that is differentially expressed in S6 in both improver and nonimprovers
GO_interest_S6 <- list (					
	neutrophil_chemotaxis='GO:0030593',
	positive_regulation_of_leukocyte_migration='GO:0002687',
	myeloid_leukocyte_differentiation='GO:0002573',
	cellular_response_to_interferon_gamma='GO:0071346',
	cellular_response_to_lipopolysaccharide='GO:0071222',
	cellular_response_to_tumor_necrosis_factor='GO:0071356',
	positive_regulation_of_leukocyte_cell_cell_adhesion='GO:1903039',
	extracellular_matrix_organization='GO:0030198',
	response_to_virus='GO:0009615',
	positive_regulation_of_response_to_external_stimulus='GO:0032103',
	neutrophil_degranulation='GO:0043312',
	positive_regulation_of_leukocyte_activation='GO:0002696',
	positive_regulation_of_cell_population_proliferation='GO:0008284')
GO_interest_S6 <- list (
					neutrophil_degranulation='GO:0043312',
					neutrophil_chemotaxis='GO:0030593',
					cellular_response_to_lipopolysaccharide='GO:0071222',
					chemokine_mediated_signaling_pathway='GO:0070098',
					extracellular_matrix_organization='GO:0030198',
					regulation_of_signaling_receptor_activity='GO:0010469',
					positive_regulation_of_cytokine_secretion='GO:0050715',
					acute_inflammatory_response='GO:0002526',
					positive_regulation_of_acute_inflammatory_response_='GO:0002675',
					positive_regulation_of_apoptotic_process='GO:0043065',
					response_to_toxic_substance='GO:0009636',
					cellular_response_to_tumor_necrosis_factor='GO:0071356',
					response_to_virus='GO:0009615',
					lymphocyte_chemotaxis='GO:0048247',
					regulation_of_leukocyte_apoptotic_process='GO:2000106',
					positive_regulation_of_ERK1_and_ERK2_cascade='GO:0070374',
					angiogenesis='GO:0001525',
					response_to_inorganic_substance='GO:0010035',
					cellular_defense_response='GO:0006968',
					monocyte_chemotaxis='GO:0002548',
					negative_regulation_of_interleukin_12_production_='GO:0032695',
					antimicrobial_humoral_immune_response_mediated_by_antimicrobial_peptide='GO:0061844',
					positive_regulation_of_cytosolic_calcium_ion_concentration='GO:0007204')		
GO_genes_S6 <- lapply(GO_interest_S6, FUN=function(x)unique(getBM(attributes = c('entrezgene','hgnc_symbol'), 
												  filters = 'go', 
												  values = x, 
												  mart = mart)$hgnc_symbol))

GO_genes_S6 <- lapply(GO_genes_S6, FUN=function(x)unique(x))

GO_interest_D2 <- list (
					positive_regulation_of_cell_killing='GO:0031343',
					positive_regulation_of_chemotaxis='GO:0050921',
					positive_regulation_of_immune_effector_process='GO:0002699',
					cellular_response_to_interferon_gamma='GO:0071346',
					cellular_response_to_tumor_necrosis_factor='GO:0071356',
					positive_regulation_of_innate_immune_response='GO:0045089',
					inflammatory_response='GO:0006954',
					regulation_of_T_cell_activation='GO:0050863',
					cytokine_mediated_signaling_pathway='GO:0019221',
					adaptive_immune_response='GO:0002250',
					leukocyte_activation_involved_in_immune_response='GO:0002366',
					leukocyte_mediated_immunity='GO:0002443',
					response_to_other_organism='GO:0051707')
GO_interest_D2 <- list (
					positive_regulation_of_alpha_beta_T_cell_proliferation='GO:0046641',	
					negative_regulation_of_interleukin_6_production='GO:0032715',	
					regulation_of_mononuclear_cell_migration='GO:0071675',	
					T_cell_costimulation='GO:0031295',	
					negative_regulation_of_viral_genome_replication='GO:0045071',	
					cellular_defense_response='GO:0006968',	
					positive_regulation_of_interferon_gamma_production='GO:0032729',	
					type_I_interferon_signaling_pathway='GO:0060337',	
					interferon_gamma_mediated_signaling_pathway='GO:0060333',	
					positive_regulation_of_tumor_necrosis_factor_production='GO:0032760',	
					chemokine_mediated_signaling_pathway='GO:0070098',	
					neutrophil_chemotaxis='GO:0030593',	
					myeloid_leukocyte_differentiation='GO:0002573',	
					T_cell_receptor_signaling_pathway='GO:0050852',	
					positive_regulation_of_cytokine_secretion='GO:0050715',	
					positive_regulation_of_NF_kappaB_transcription_factor_activity='GO:0051092',	
					regulation_of_cytokine_mediated_signaling_pathway='GO:0001959',	
					positive_regulation_of_I_kappaB_kinase_NF_kappaB_signaling='GO:0043123',	
					defense_response_to_virus='GO:0051607',	
					positive_regulation_of_ERK1_and_ERK2_cascade='GO:0070374',	
					positive_regulation_of_cytosolic_calcium_ion_concentration='GO:0007204',	
					inflammatory_response='GO:0006954',	
					neutrophil_degranulation='GO:0043312')
GO_genes_D2 <- lapply(GO_interest_D2, FUN=function(x)unique(getBM(attributes = c('entrezgene','hgnc_symbol'), 
												  filters = 'go', 
												  values = x, 
												  mart = mart)$hgnc_symbol))

GO_genes_D2 <- lapply(GO_genes_D2, FUN=function(x)unique(x))
