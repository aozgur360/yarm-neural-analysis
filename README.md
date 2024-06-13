"# yarm-neural-analysis" 

svm_allmice_averages_main:
to get main decoding using a support vector machine fed spikerate sums, ex. fig23_svm_decoding, can change window in get_spikerate_sums_per_cell_per_trial function, ex. fig24_svm_decoding_window_timebins _wstats
![fig23_svm_decoding](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/1db8d79a-ad15-4a14-87c1-d54cd493f05b)

![fig24_svm_decoding_window_timebins _wstats](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/422a4781-7e0a-46af-9ba3-4c6494c3ebd6)

clustering_main:
for clustering thousands of neurons across mice to find similarities in the process of learning, datasets used: mice_datasets\mice_params/mice_INIT_TF_params_406065_231128.csv (INIT, 0.2 threshold), mice_datasets\mice_params/mice_FTP_TF_params_406065_231130.csv' (FTP, 0.2 threshold), mice_datasets\crossreg\all_cell_matching_240130\init_n_e1 (supplemental INIT), mice_datasets\crossreg\all_cell_matching_240130\ftp_n_e1_nobias (supplemental FTP)
![fig27_init_clustering](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/dfb6ebf6-0471-4128-903d-d95ba3c0007c)

![fig28_init_clustering_stats](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/3d466ea3-9d4f-419c-ad14-c4d4bb72cab8)

![fig31_init_clustering_supp_sizes](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/e634668c-24b1-4f85-aea4-a779e306b3ec)


cell_matching_analysis_graphs_and_decoding_main:
for active cell turnover rates and cell matched decoding across sessions, datasets used: mice_datasets\crossreg\all_cell_matching_240130\[INSERST DESIRED: ex. init_e1_e2]\{mouse_id}\mappings.pkl'
![fig34_cellmatching_stability](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/176b2928-c084-4dca-b704-78ce6f3b63d5)

![fig35_cellmatching_decoding](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/6c8cebff-1285-40af-a670-33b91008d15c)

![fig36_initeeres_indiv](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/4e61aa3e-0b6c-4af6-b3b6-4746638d976e)

![fig37_initee_indiv_mice_scatters](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/5077dd53-2479-4170-9af8-31a278e637b0)

cross_reg_clustering_and_heatmaps_main: 
for cell matched clustering and heatmaps, datasets used: mice_datasets/mice_params/1to1_INIT_TF_ee.csv" (INIT ee), mice_datasets/mice_params/1to1_INIT_TF_ne.csv" (INIT ne)
![fig38_cell_matched_clustering](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/dad42f28-1c27-421e-a9f2-d6287cd52b00)

![fig39_heatmaps_indiv_mice_INIT_ee](https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/b9bcdb6b-f22b-4cb6-893d-95f056669a11)

https://github.com/aozgur360/yarm-neural-analysis/assets/77759136/b76869f6-dddb-493f-afaa-1aafea79778c

mice used when not specified by csv dataset: = ['m2r','m9r','m11n','m12lr','m13l','m14r','m15n','m16l','m17r','m19n','m21r','m22l','m23r','m25lr','m26r','m29n', 'm30lr', 'm31r','m32l','m33r','m34n','m36rr','m37l','m38n','m39rr','m40l', 'm41r','min4r','sn','snr']
