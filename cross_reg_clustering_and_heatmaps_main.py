#%%
#NOTE: MAKE SURE TO SPECIFY TURN DIRECTION VS. VIS STIM FOR LR (comment & uncomment in functions section accordingly.. this code only looks at 1 of these at a time, not both)

#%%
def average_absolute_difference(trace1, trace2):
    # Check if the traces have the same length
    if len(trace1) != len(trace2):
        raise ValueError("Traces must have the same length")

    # Calculate absolute differences for each index
    differences = [abs(trace1[i] - trace2[i]) for i in range(len(trace1))]

    # Calculate the average of the absolute differences
    average_difference = sum(differences) / len(differences)

    return average_difference

# Example usage:
#trace1 = [1, 2, 3, 4, 5]
#trace2 = [5, 4, 3, 2, 1]
#result = average_absolute_difference(trace1, trace2)
#print("Average Absolute Difference:", result)

#%%
import pandas as pd

# Assuming TF_params is your DataFrame

# Iterate through pairs of indices
for i in range(0, len(TF_params), 2):
    # Select rows at current and next index
    pair = TF_params.iloc[i:i+2]  # This will select the row at index i and i+1
    # Do something with the pair
    print(pair)
    # If you want to perform operations on this pair, you can do so here


#%%
#CURRENTLY USING TF_PARAMS_PATH FOR FTP
stage = 's2'

#0 for no, 1 for yes
if stage == 's2':
  init_only = 1
else:
  init_only = 0
#0 for turn decoding or visual stimuli decoding, 1 for correctness decoding
# IGNORE TRIAL TYPE
trial_type = 0


import pandas as pd
import os
import glob
import pandas as pd
import math
import numpy as np
import re
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.pyplot import cm
from scipy.stats import linregress

#TF_params_path = r"E:/main_data/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/1to1_INIT_TF_ee.csv" #INIT ee
TF_params_path = r"E:/main_data/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/1to1_INIT_TF_ne.csv" #INIT ne
#TF_params_path = r"E:/main_data/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/mice_INIT_FTF_params_231012.csv" #INIT

#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/mice_FTP_TF_custom7m_1n1e_231009.csv' #FTP
#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/mice_INIT_FTF_params_231012.csv' #INIT
#TF_params_path = r"E:/main_data/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/edrive_test_mice_INIT_FTF_params_231024.csv" #INIT lenient
#TF_params_path = r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params\for_clustering_sessions_for_cell_matching_init_n_e1_240130.csv"

#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/init_1n1e_permouse.csv' #5 init mice, 1n/e (no bias, >20 clean trials)
#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/mice_INIT_TF_params.csv' #all init (no bias, >20 clean trials)
#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/mice_FTP_TF_params.csv' #all ftp (no bias, >20 clean trials)
TF_params = pd.read_csv(TF_params_path)

#################################################
# # #### if want to look at only single mouse

# pair_number = 5  # For the 3rd pair, pair = mouse number
# start_index = 2 * (pair_number - 1)

# # Overwrite TF_params with just the 3rd pair
# TF_params = TF_params.iloc[start_index:start_index + 2]

# # Assuming you have already overwritten TF_params with just the 3rd pair
# # Now reset the indices to start at 0
# TF_params = TF_params.reset_index(drop=True)

# # TF_params now only contains the data for the 3rd pair, with indices starting at 0
# print(TF_params)
#################################################


# Initialize a list to store sessionID values
selected_sessions = []
selected_sessions_expert = []
selected_sessions_naive = []

#crossreg_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/m32l/crossreg_230731_230814/mappings.pkl'
#mappings = pd.read_pickle(crossreg_path)
#pairings = mappings.iloc[:, [1, 2]].values.tolist()
#pairings = [pair for pair in pairings if not (np.isnan(pair[0]) or np.isnan(pair[1]))]
#pairings = [[round(pair[0]), round(pair[1])] for pair in pairings]

reject_dict = {}
reject_dict_naive_sessions = {}
reject_dict_expert_sessions = {}
# Iterate through the rows of the new DataFrame 'new_df'
for index, row in TF_params.iterrows():
    bias = row['bias']
    enough_trials = row['enough_trials']
    expert = row['expert']
    sessionID = row['sessionID']

    # Check if the conditions are met
    if not bias and enough_trials and not pd.isna(expert):
        selected_sessions.append(sessionID)

    # get above conditions but also only expert
    if not bias and enough_trials and not pd.isna(expert) and expert:
        selected_sessions_expert.append(sessionID)

    # get above conditions but also only naive
    if not bias and enough_trials and not pd.isna(expert) and not expert:
        selected_sessions_naive.append(sessionID)

# Print the list of selected sessionIDs
#print(selected_sessions)
#########################################
temp = selected_sessions

### IF USING DESKTOP THEN UNCOMMENT BELOW:
#temp = [path.replace('/content/drive/MyDrive', 'E:/main_data') for path in temp]
#temp = [os.path.join(path.replace('/content/drive/MyDrive', 'E:/main_data')) for path in temp]



behavior = []
folders_analyzed = []
temp.sort()
param_paths = []
cal_paths_raw = []
cal_paths_spike = []
rejected_cells_paths = []


for file in temp:
    param_paths.append(glob.glob(file + '/trial_parameters' + '*.csv'))
    #cal_paths_spike.append(glob.glob(file + '/spikerate' + '*.csv'))
    cal_paths_raw.append(glob.glob(file + '/calcium' + '*.csv'))
    rejected_cells_paths.append(glob.glob(file + '/rejected' + '*.csv'))

data_transposed_C_arrays_expert = []
data_transposed_W_arrays_expert = []
data_transposed_L_arrays_expert = []
data_transposed_R_arrays_expert = []
data_transposed_C_arrays_naive = []
data_transposed_W_arrays_naive = []
data_transposed_L_arrays_naive = []
data_transposed_R_arrays_naive = []


transposed_data_arrays = []
session_info_list = []
for i in range(len(param_paths)):
#for i in range(0,10):
#for i in range(8,11):
    #try:
      current_folder = temp[i]
      cal = pd.read_csv(cal_paths_raw[i][0])
      #cal_spike = pd.read_csv(cal_paths_spike[i][0])
      params = pd.read_csv(param_paths[i][0])
      reject = pd.read_csv(rejected_cells_paths[i][0])
      reject = [item for sublist in reject.values.tolist() for item in sublist]
      reject_dict[current_folder] = reject
      if current_folder in selected_sessions_naive:
        reject_dict_naive_sessions[current_folder] = reject
      if current_folder in selected_sessions_expert:
        reject_dict_expert_sessions[current_folder] = reject


      del cal[cal.columns[0]]

      ##blood vessels remove
      #cal.drop(cal.columns[reject], axis=1, inplace=True)
      ##reset column indices / cell ids
      #cal = cal.transpose()
      #cal = cal.reset_index(drop=True)
      #cal = cal.transpose()

      #del cal_spike[cal_spike.columns[0]]
      ##blood vessels remove
      #cal_spike.drop(cal_spike.columns[reject], axis=1, inplace=True)
      ###reset column indices / cell ids
      #cal_spike = cal_spike.transpose()
      #cal_spike = cal_spike.reset_index(drop=True)
      #cal_spike = cal_spike.transpose()

      #get df/f new approach
      f0_values = cal.quantile(0.1, axis=0)
      f0_values = np.where(f0_values == 0, 0.000000000001, f0_values)
      f0_values = np.nan_to_num(f0_values, nan=0.000000000001)
      cal = cal.div(f0_values, axis=1)


      # z score cal
      #from scipy.stats import zscore
      #cal = cal.apply(zscore)

      # Apply a 30-row moving average to each column, excluding NaN values
      #window_size = 30
      #cal = cal.apply(lambda x: x.rolling(window=2 * window_size + 1, center=True, min_periods=1).mean())

      parameters = get_parameters(cal,params)
      trace_average_per_cell = get_trace_average_per_cell(cal,parameters)

      print(current_folder)
      ################################
      #t_trace_average_per_cell = np.transpose(trace_average_per_cell)
      data = trace_average_per_cell

      from scipy.stats import zscore #MAIN KEEP
      data = zscore(data, axis=0)

      window_size = 6 #6 def
      # Create an array to store the result
      smoothed_data = np.zeros_like(data)
      # Apply the moving average to each column
      for i in range(data.shape[1]):
          column = data[:, i]
          for t in range(data.shape[0]):
              start_idx = max(0, t - window_size // 2)
              end_idx = min(data.shape[0], t + window_size // 2 + 1)
              smoothed_data[t, i] = np.mean(column[start_idx:end_idx])
      data = smoothed_data
      #print('gen')

      ################
      parameters_C = get_parameters_C(cal,params)
      if current_folder in selected_sessions_expert:
        print('ex')
        C_cal_expert = get_trace_average_per_cell(cal,parameters_C)
        C_cal_expert = zscore(C_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_C_cal_expert = np.zeros_like(C_cal_expert)
        # Apply the moving average to each column
        for i in range(C_cal_expert.shape[1]):
            column = C_cal_expert[:, i]
            for t in range(C_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(C_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_C_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        C_cal_expert = smoothed_C_cal_expert
      elif current_folder in selected_sessions_naive:
        print('na')
        C_cal_naive = get_trace_average_per_cell(cal,parameters_C)
        C_cal_naive  = zscore(C_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_C_cal_naive = np.zeros_like(C_cal_naive)
        # Apply the moving average to each column
        for i in range(C_cal_naive.shape[1]):
            column = C_cal_naive[:, i]
            for t in range(C_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(C_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_C_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        C_cal_naive = smoothed_C_cal_naive

      parameters_W = get_parameters_W(cal,params)
      if current_folder in selected_sessions_expert:
        W_cal_expert = get_trace_average_per_cell(cal,parameters_W)
        W_cal_expert = zscore(W_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_W_cal_expert = np.zeros_like(W_cal_expert)
        # Apply the moving average to each column
        for i in range(W_cal_expert.shape[1]):
            column = W_cal_expert[:, i]
            for t in range(W_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(W_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_W_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        W_cal_expert = smoothed_W_cal_expert
      elif current_folder in selected_sessions_naive:
        W_cal_naive = get_trace_average_per_cell(cal,parameters_W)
        W_cal_naive  = zscore(W_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_W_cal_naive = np.zeros_like(W_cal_naive)
        # Apply the moving average to each column
        for i in range(W_cal_naive.shape[1]):
            column = W_cal_naive[:, i]
            for t in range(W_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(W_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_W_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        W_cal_naive = smoothed_W_cal_naive


      parameters_L = get_parameters_L(cal,params)
      if current_folder in selected_sessions_expert:
        L_cal_expert = get_trace_average_per_cell(cal,parameters_L)
        L_cal_expert = zscore(L_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_L_cal_expert = np.zeros_like(L_cal_expert)
        # Apply the moving average to each column
        for i in range(L_cal_expert.shape[1]):
            column = L_cal_expert[:, i]
            for t in range(L_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(L_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_L_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        L_cal_expert = smoothed_L_cal_expert
      elif current_folder in selected_sessions_naive:
        L_cal_naive = get_trace_average_per_cell(cal,parameters_L)
        L_cal_naive  = zscore(L_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_L_cal_naive = np.zeros_like(L_cal_naive)
        # Apply the moving average to each column
        for i in range(L_cal_naive.shape[1]):
            column = L_cal_naive[:, i]
            for t in range(L_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(L_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_L_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        L_cal_naive = smoothed_L_cal_naive

      parameters_R = get_parameters_R(cal,params)
      if current_folder in selected_sessions_expert:
        R_cal_expert = get_trace_average_per_cell(cal,parameters_R)
        R_cal_expert = zscore(R_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_R_cal_expert = np.zeros_like(R_cal_expert)
        # Apply the moving average to each column
        for i in range(R_cal_expert.shape[1]):
            column = R_cal_expert[:, i]
            for t in range(R_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(R_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_R_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        R_cal_expert = smoothed_R_cal_expert
      elif current_folder in selected_sessions_naive:
        R_cal_naive = get_trace_average_per_cell(cal,parameters_R)
        R_cal_naive  = zscore(R_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_R_cal_naive = np.zeros_like(R_cal_naive)
        # Apply the moving average to each column
        for i in range(R_cal_naive.shape[1]):
            column = R_cal_naive[:, i]
            for t in range(R_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(R_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_R_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        R_cal_naive = smoothed_R_cal_naive

      #######################################################################

      # Transpose the data array to have cell IDs as rows and frames as columns
      data_transposed = data.T
      transposed_data_arrays.append(data_transposed)


      if current_folder in selected_sessions_expert:
          data_transposed_C_expert = C_cal_expert.T
          data_transposed_W_expert = W_cal_expert.T
          data_transposed_L_expert = L_cal_expert.T
          data_transposed_R_expert = R_cal_expert.T
          data_transposed_C_arrays_expert.append(data_transposed_C_expert)
          data_transposed_W_arrays_expert.append(data_transposed_W_expert)
          data_transposed_L_arrays_expert.append(data_transposed_L_expert)
          data_transposed_R_arrays_expert.append(data_transposed_R_expert)

          copy_shape = data_transposed_C_expert.shape
          filler = "filler"
          filler_array = np.full(copy_shape, filler)
          data_transposed_C_arrays_naive.append(filler_array)
          data_transposed_W_arrays_naive.append(filler_array)
          data_transposed_L_arrays_naive.append(filler_array)
          data_transposed_R_arrays_naive.append(filler_array)

      elif current_folder in selected_sessions_naive:
          data_transposed_C_naive = C_cal_naive.T
          data_transposed_W_naive = W_cal_naive.T
          data_transposed_L_naive = L_cal_naive.T
          data_transposed_R_naive = R_cal_naive.T
          data_transposed_C_arrays_naive.append(data_transposed_C_naive)
          data_transposed_W_arrays_naive.append(data_transposed_W_naive)
          data_transposed_L_arrays_naive.append(data_transposed_L_naive)
          data_transposed_R_arrays_naive.append(data_transposed_R_naive)

          copy_shape = data_transposed_C_naive.shape
          filler = "filler"
          filler_array = np.full(copy_shape, filler)
          data_transposed_C_arrays_expert.append(filler_array)
          data_transposed_W_arrays_expert.append(filler_array)
          data_transposed_L_arrays_expert.append(filler_array)
          data_transposed_R_arrays_expert.append(filler_array)

    #except:
    #  pass


concatenated_data = np.vstack(transposed_data_arrays)

concatenated_data_C_expert = np.vstack(data_transposed_C_arrays_expert)
concatenated_data_W_expert = np.vstack(data_transposed_W_arrays_expert)
concatenated_data_L_expert = np.vstack(data_transposed_L_arrays_expert)
concatenated_data_R_expert = np.vstack(data_transposed_R_arrays_expert)

concatenated_data_C_naive = np.vstack(data_transposed_C_arrays_naive)
concatenated_data_W_naive = np.vstack(data_transposed_W_arrays_naive)
concatenated_data_L_naive = np.vstack(data_transposed_L_arrays_naive)
concatenated_data_R_naive = np.vstack(data_transposed_R_arrays_naive)

# Initialize a list to store the information for each row in concatenated_data
data_info = []

# Initialize variables to keep track of original row indices and file names
current_file_index = 0
original_row_index_offset = 0

# Iterate through the rows of concatenated_data and keep track of original row indices and file names
for i, concatenated_row in enumerate(concatenated_data):
    file_name = temp[current_file_index]
    original_row_index = i - original_row_index_offset
    data_info.append([original_row_index, i, file_name])

    # Check if we have reached the end of the current file's data_transposed
    if original_row_index == len(data_transposed_C_arrays_expert[current_file_index]) - 1:
        current_file_index += 1
        original_row_index_offset = i + 1

# Now data_info contains pairs of original row indices from each file and corresponding row indices in concatenated_data

# Determine which rows have NaNs
nan_rows = []
nan_indices = np.isnan(concatenated_data).any(axis=1)

# Determine which rows have NaNs
for index, has_nan in enumerate(nan_indices):
    if has_nan:
        nan_rows.append(index)

# Remove rows with NaNs
concatenated_data = np.delete(concatenated_data, nan_rows, axis=0)
concatenated_data_C_expert = np.delete(concatenated_data_C_expert, nan_rows, axis=0)
concatenated_data_W_expert = np.delete(concatenated_data_W_expert, nan_rows, axis=0)
concatenated_data_L_expert = np.delete(concatenated_data_L_expert, nan_rows, axis=0)
concatenated_data_R_expert = np.delete(concatenated_data_R_expert, nan_rows, axis=0)
concatenated_data_C_naive = np.delete(concatenated_data_C_naive, nan_rows, axis=0)
concatenated_data_W_naive = np.delete(concatenated_data_W_naive, nan_rows, axis=0)
concatenated_data_L_naive = np.delete(concatenated_data_L_naive, nan_rows, axis=0)
concatenated_data_R_naive = np.delete(concatenated_data_R_naive, nan_rows, axis=0)


expert_cellids = []
naive_cellids = []

for index, row in enumerate(concatenated_data_C_expert):
    if "filler" not in row:
        expert_cellids.append(index)
for index, row in enumerate(concatenated_data_C_naive):
    if "filler" not in row:
        naive_cellids.append(index)

#%%
import pandas as pd
import numpy as np
import os

#crossreg_paths = ['/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/FTP_initd1e_1n1e_231009/m11n/crossreg_220506_220520/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/FTP_initd1e_1n1e_231009/m12lr/crossreg_220505_220512/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/FTP_initd1e_1n1e_231009/m13l/crossreg_220505_220512/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/FTP_initd1e_1n1e_231009/m14r/crossreg_220603_220613/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/FTP_initd1e_1n1e_231009/m25lr/crossreg_221107_221110/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/FTP_initd1e_1n1e_231009/m26r/crossreg_221107_221213/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/FTP_initd1e_1n1e_231009/m9r/crossreg_220505_220518/mappings.pkl']

#crossreg_paths = ['/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/INIT_1n1e_231012/m21r/crossreg_221013_221110/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/INIT_1n1e_231012/m32l/crossreg_230731_230815/mappings.pkl',
#                  '/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/INIT_1n1e_231012/m33r/crossreg_230808_230901/mappings.pkl']

# crossreg_paths = ['E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m15n/crossreg_220615_220617/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m16l/crossreg_220921_220927/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m19n/crossreg_220919_221007/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m22l/crossreg_221028_221219/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m25lr/crossreg_221205_221221/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m26r/crossreg_221214_221219/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m30lr/crossreg_230801_230816/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m31r/crossreg_230731_230801/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m32l/crossreg_230731_230816/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m33r/crossreg_230808_230901/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m40l/crossreg_230929_231004/mappings.pkl',
#                   'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\INIT_1n1e_231024/m41r/crossreg_230928_231003/mappings.pkl']

# crossreg_paths = [r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\all_cell_matching_240130\init_n_e1/m21r/mappings.pkl",
#                   r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\all_cell_matching_240130\init_n_e1/m22l/mappings.pkl",
#                   r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\all_cell_matching_240130\init_n_e1/m32l/mappings.pkl",
#                   r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\all_cell_matching_240130\init_n_e1/m33r/mappings.pkl",
#                   r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\all_cell_matching_240130\init_n_e1/m41r/mappings.pkl"]

#for 1to1_INIT_TF_ne
crossreg_paths = [r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_n_e1/m21r/mappings.pkl",
                  r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_n_e1/m22l/mappings.pkl",
                  r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_n_e1/m32l/mappings.pkl",
                  r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_n_e1/m33r/mappings.pkl",
                  r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_n_e1/m41r/mappings.pkl"]

# # #for 1to1_INIT_TF_ee
# crossreg_paths = [r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m11n/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m12lr/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m13l/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m14r/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m15n/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m16l/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m25lr/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m26r/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m30lr/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m31r/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m32l/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m33r/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m38n/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m39rr/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m40l/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m41r/mappings.pkl",
#                   r"E:/main_data\Calcium_imaging_yarm_LurLab\mice_datasets/crossreg/all_cell_matching_240130/init_e1_e2/m9r/mappings.pkl"]







########
import re

reject_dict_naive = {}
reject_dict_expert = {}

# ##for old structure
# for key, value in reject_dict_naive_sessions.items():
#     folder_name = key.split('/')[-3]  # Extract the desired folder name
#     reject_dict_naive[folder_name] = value
# for key, value in reject_dict_expert_sessions.items():
#     folder_name = key.split('/')[-3]  # Extract the desired folder name
#     reject_dict_expert[folder_name] = value

# for new e drive structure
for key, value in reject_dict_naive_sessions.items():
    #print(key,value)
    folder_name = key.split('\\')[-3]  # Extract the desired folder name
    #print(folder_name)
    reject_dict_naive[folder_name] = value
for key, value in reject_dict_expert_sessions.items():
    folder_name = key.split('\\')[-3]  # Extract the desired folder name
    reject_dict_expert[folder_name] = value

# for value in crossreg_paths:
#     folder_name = os.path.basename(os.path.dirname(value))
#     reject_dict_naive[value] = folder_name
#     reject_dict_expert[value] = folder_name

# Initialize dictionaries to store the lists
expert_matched_cellids_dict = {}
naive_matched_cellids_dict = {}
expert_unmatched_cellids_dict = {}
naive_unmatched_cellids_dict = {}

for file in crossreg_paths:
    # Extract the 'm11n', 'm12lr', etc., part from the filename
    match = re.search(r'm\d+[a-z]+', file)
    #print(match)
    if match:
        key = match.group()  # Use the matched part as the key

        mappings = pd.read_pickle(file)
        n_expert_cells = mappings.iloc[:, 2].count()
        expert_unmatched_cells = list(range(n_expert_cells))
        n_naive_cells = mappings.iloc[:, 1].count()
        naive_unmatched_cells = list(range(n_naive_cells))


        pairings = mappings.iloc[:, [1, 2]].values.tolist()
        pairings = [pair for pair in pairings if not (np.isnan(pair[0]) or np.isnan(pair[1]))]
        pairings = [[round(pair[0]), round(pair[1])] for pair in pairings]
        #print(pairings)
        expert_matched_cellids = []
        naive_matched_cellids = []

        for pairing in pairings:
            expert_matched_cellids.append(pairing[1])
            naive_matched_cellids.append(pairing[0])

        expert_unmatched_cellids = [cellid for cellid in expert_unmatched_cells if cellid not in expert_matched_cellids]
        naive_unmatched_cellids = [cellid for cellid in naive_unmatched_cells if cellid not in naive_matched_cellids]

        # Save the lists in the dictionaries with the appropriate key
        expert_matched_cellids_dict[key] = expert_matched_cellids
        naive_matched_cellids_dict[key] = naive_matched_cellids
        expert_unmatched_cellids_dict[key] = expert_unmatched_cellids
        naive_unmatched_cellids_dict[key] = naive_unmatched_cellids
    else:
        print(f"Skipping file '{file}' as it doesn't have the expected format.")
######
#this will give the final cellids in the massive array for matched or unmatched cells in expert sessions

data_info_expert_sessions = []
data_info_naive_sessions = []
for item in data_info:
    if item[2] in selected_sessions_expert:
        data_info_expert_sessions.append(item)
    elif item[2] in selected_sessions_naive:
        data_info_naive_sessions.append(item)


#mice_ids = ['m21r','m22l','m32l','m33r','m41r'] # for INIT_TF_ne
#mice_ids = ['m11n',''] # for INIT_TF_ee
# Assuming TF_params is your DataFrame
mice_ids = TF_params['mouseID'].unique().tolist()

# Now mice_ids contains all unique mouse IDs
print(mice_ids)


expert_matched_cells_final = []
expert_unmatched_cells_final = []
naive_matched_cells_final = []
naive_unmatched_cells_final = []

# ##for old structure
# for m in mice_ids:
#   for d in data_info_expert_sessions:
#     if m in d[2]:
#         if any(value == d[0] for value in expert_matched_cellids_dict[m]):
#           if d[0] not in reject_dict_expert[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#             expert_matched_cells_final.append(d[1]) #d1 is the cell ID within the massive collective array
#             #print(d[0]) #d0 is the original cell ID for that session
#         if any(value == d[0] for value in expert_unmatched_cellids_dict[m]):
#           if d[0] not in reject_dict_expert[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#             expert_unmatched_cells_final.append(d[1])
#   for d in data_info_naive_sessions:
#     if m in d[2]:
#       if any(value == d[0] for value in naive_matched_cellids_dict[m]):
#         if d[0] not in reject_dict_naive[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#           naive_matched_cells_final.append(d[1])
#       if any(value == d[0] for value in naive_unmatched_cellids_dict[m]):
#         if d[0] not in reject_dict_naive[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#           naive_unmatched_cells_final.append(d[1])

#for new edrive structure
for m in mice_ids:
  for d in data_info_expert_sessions:
    if m in d[2]:
        if any(value == d[0] for value in expert_matched_cellids_dict[m]):
          if d[0] not in reject_dict_expert[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
            expert_matched_cells_final.append(d[1]) #d1 is the cell ID within the massive collective array
            #print(d[0]) #d0 is the original cell ID for that session
        if any(value == d[0] for value in expert_unmatched_cellids_dict[m]):
          if d[0] not in reject_dict_expert[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
            expert_unmatched_cells_final.append(d[1])
  for d in data_info_naive_sessions:
    if m in d[2]:
      if any(value == d[0] for value in naive_matched_cellids_dict[m]):
        if d[0] not in reject_dict_naive[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
          naive_matched_cells_final.append(d[1])
      if any(value == d[0] for value in naive_unmatched_cellids_dict[m]):
        if d[0] not in reject_dict_naive[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
          naive_unmatched_cells_final.append(d[1])

# #for new edrive structure with IGNORE rejects
# for m in mice_ids:
#   for d in data_info_expert_sessions:
#     if m in d[2]:
#         if any(value == d[0] for value in expert_matched_cellids_dict[m]):
#           #if d[0] not in reject_dict_expert[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#             expert_matched_cells_final.append(d[1]) #d1 is the cell ID within the massive collective array
#             #print(d[0]) #d0 is the original cell ID for that session
#         if any(value == d[0] for value in expert_unmatched_cellids_dict[m]):
#           #if d[0] not in reject_dict_expert[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#             expert_unmatched_cells_final.append(d[1])
#   for d in data_info_naive_sessions:
#     if m in d[2]:
#       if any(value == d[0] for value in naive_matched_cellids_dict[m]):
#         #if d[0] not in reject_dict_naive[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#           naive_matched_cells_final.append(d[1])
#       if any(value == d[0] for value in naive_unmatched_cellids_dict[m]):
#         #if d[0] not in reject_dict_naive[m]: #this is saying to only continue with the cells in the session that were not rejected due to size.
#           naive_unmatched_cells_final.append(d[1])


#%%
key

#%%
# remove rejected cells from dict, in order to match cell ids after clustering, per mouse for decoding analysis
for key in expert_matched_cellids_dict:
    expert_matched_cellids_dict[key] = [x for x in expert_matched_cellids_dict[key] if x not in reject_dict_expert.get(key, [])]
for key in expert_unmatched_cellids_dict:
    expert_unmatched_cellids_dict[key] = [x for x in expert_unmatched_cellids_dict[key] if x not in reject_dict_expert.get(key, [])]
for key in naive_matched_cellids_dict:
    naive_matched_cellids_dict[key] = [x for x in naive_matched_cellids_dict[key] if x not in reject_dict_naive.get(key, [])]
for key in naive_unmatched_cellids_dict:
    naive_unmatched_cellids_dict[key] = [x for x in naive_unmatched_cellids_dict[key] if x not in reject_dict_naive.get(key, [])]

#%%
# printout = expert_matched_cellids_dict['m22l']
# print('[' + ', '.join(map(str, printout)) + ']')

#%%
printout = expert_matched_cellids_dict['m11n']
print('[' + ', '.join(map(str, printout)) + ']')
print(len(printout))

#%%
printout = naive_matched_cellids_dict['m11n']
print('[' + ', '.join(map(str, printout)) + ']')
print(len(printout))

#%%
len(expert_matched_cells_final)

#%%
len(naive_matched_cells_final)

#%%
len(printout)

#%%
expert_matched_cells_final_array

#%%


#%%
# to check step , confirmed
#expert_matched_cellids_dict
#naive_matched_cellids_dict
#data_info
#expert_matched_cells_final
#naive_matched_cells_final

#%%
#total_count = sum(len(values) for values in expert_matched_cellids_dict.values())
#print("Total number of values across all keys:", total_count)

#%%

expert_matched_cells_final_array = concatenated_data[expert_matched_cells_final]

# cluster using all data to get cellIDs
###### cluster
# Calculate the linkage matrix using 'ward' method
#Z = linkage(data_transposed, method='ward')
Z = linkage(expert_matched_cells_final_array, method='ward')

# Set the percentage threshold for determining the number of clusters
percentage_threshold = 0.25  # Adjust this value as needed, 0.2

# Calculate the height threshold based on the percentage
max_d = percentage_threshold * np.max(Z[:, 2])

# Retrieve the clusters based on the height threshold
clusters = fcluster(Z, max_d, criterion='distance')
num_clusters = np.max(clusters)

# Generate a random color for each cluster
cmap = cm.rainbow(np.linspace(0, 1, num_clusters))
colors = [mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap]

# Create a color dictionary for each cluster
color_dict = {cluster_id: color for cluster_id, color in zip(range(1, num_clusters + 1), colors)}

# Assign colors to the dendrogram
set_link_color_palette(colors)

# Plot the dendrogram with colored clusters
#plt.figure(figsize=(10, 6))
#dendrogram(Z, color_threshold=max_d, above_threshold_color='gray')

# Assign unique integer labels to leaf nodes
leaf_labels = dendrogram(Z, no_plot=True)['ivl']
label_map = {label: idx for idx, label in enumerate(leaf_labels)}
leaf_labels_int = [label_map[label] for label in leaf_labels]

## Apply colors to the leaf nodes (cell IDs) based on the cluster they belong to
#leaf_colors = [color_dict[clusters[label]] for label in leaf_labels_int]
#ax = plt.gca()
#ax.set_xticklabels(ax.get_xticks(), rotation=90)  # Rotate x-axis labels for better readability
#plt.bar(range(len(leaf_labels_int)), np.zeros(len(leaf_labels_int)), color=leaf_colors, align='center', width=0.5)
## Plot a horizontal line at the height threshold
#plt.axhline(y=max_d, c='red', linestyle='dashed')
#plt.title('Hierarchical Clustering Dendrogram')
#plt.xlabel('Cell IDs')
#plt.ylabel('Distance')

# Dictionary to store cluster indices
#cluster_indices_dict = {}
#ncells_per_cluster = []

# Print the cell IDs for each cluster in list format and store indices
#for cluster_id in range(1, num_clusters + 1):
#    cluster_indices = np.where(clusters == cluster_id)[0]
#    cell_ids = [int(leaf_labels_int[idx]) for idx in cluster_indices]
#    cluster_name = f"cluster{cluster_id}"
#    cluster_indices_dict[cluster_name] = cell_ids
#    print(f"{cluster_name} = {cell_ids}")
#    ncells_per_cluster.append(len(cell_ids))
#cluster_expert_dict = cluster_indices_dict

## Plot the representative cell traces for each cluster on the same graph
#plt.figure(figsize=(10, 6))
#for cluster_id in range(1, num_clusters + 1):
#    cluster_indices = np.where(clusters == cluster_id)[0]
#    cluster_data = concatenated_data[cluster_indices]
#    representative_trace = np.mean(cluster_data, axis=0)
#    plt.plot(representative_trace, color=color_dict[cluster_id], label=f'Cluster {cluster_id}')
#plt.title('Representative Cell Traces for Each Cluster')
#plt.xlabel('Frames')
#plt.ylabel('Amplitude')
#plt.legend()
#plt.show()


### Create separate graphs for individual cell traces in each cluster
#for cluster_id in range(1, num_clusters + 1):
#    cluster_indices = np.where(clusters == cluster_id)[0]
#    cluster_data = concatenated_data[cluster_indices]
#    plt.figure(figsize=(10, 6))
#    for cell_trace in cluster_data:
#        plt.plot(cell_trace, color=color_dict[cluster_id], alpha=0.3)
#    plt.title(f'Cluster {cluster_id} - Individual Cell Traces')
#    plt.xlabel('Frames')
#    plt.ylabel('Amplitude')
#    plt.show()
# Create dictionaries to store the cluster-specific expert and naive lists
cluster_indices_dict = {}
ncells_per_cluster = []
for cluster_id in range(1, num_clusters + 1):
    cluster_indices = np.where(clusters == cluster_id)[0]
    cell_ids = [expert_matched_cells_final[idx] for idx in cluster_indices]  # Replace with values from the original list
    cluster_name = f"cluster{cluster_id}"
    cluster_indices_dict[cluster_name] = cell_ids
    #print(f"{cluster_name} = {cell_ids}")
    ncells_per_cluster.append(len(cell_ids))
cluster_expert_dict = cluster_indices_dict

# #for each value in cluster_expert_dict, first find that value in the list expert_matched_cells_final, and determine its index.
# #then use that index to find what value corresponds to the index in naive_matched_cells_final. save these values in a new dictionary called cluster_naive_dict
################################################################################################ OLD
# cluster_naive_dict = {}
# # Iterate through the clusters in cluster_expert_dict
# for cluster_name, cluster in cluster_expert_dict.items():
#     cluster_naive_dict[cluster_name] = []
#     for value in cluster:
#         # Find the index of the value in expert_matched_cells_final
#         expert_index = expert_matched_cells_final.index(value)
#         # Use the expert_index to get the corresponding value from naive_matched_cells_final
#         naive_value = naive_matched_cells_final[expert_index]
#         # Append the naive_value to the cluster_naive_dict for the current cluster
#         cluster_naive_dict[cluster_name].append(naive_value)


# # Print the results
# print("Expert Matched Cellids:")
# for key, value in cluster_expert_dict.items():
#     print(f'{key}: {value}')

# print("\nNaive Cellids:")
# for key, value in cluster_naive_dict.items():
#     print(f'{key}: {value}')
################################################################################################ OLD
cluster_naive_dict = {}
expert_cluster_lengths = {}
naive_cluster_lengths = {}

# Iterate through the clusters in cluster_expert_dict
for cluster_name, cluster in cluster_expert_dict.items():
    cluster_naive_dict[cluster_name] = []
    expert_cluster_lengths[cluster_name] = len(cluster)
    for value in cluster:
        try:
            # Find the index of the value in expert_matched_cells_final
            expert_index = expert_matched_cells_final.index(value)
            # Use the expert_index to get the corresponding value from naive_matched_cells_final
            naive_value = naive_matched_cells_final[expert_index]
            # Append the naive_value to the cluster_naive_dict for the current cluster
            cluster_naive_dict[cluster_name].append(naive_value)
        except IndexError:
            # Handle the case where the index is out of range in naive_matched_cells_final
            # Skip this index and continue to the next iteration
            pass
    naive_cluster_lengths[cluster_name] = len(cluster_naive_dict[cluster_name])

# # Print the lengths of expert clusters
# print("Length of Expert Clusters:")
# for key, value in expert_cluster_lengths.items():
#     print(f'{key}: {value}')

# # Print the lengths of naive clusters
# print("\nLength of Naive Clusters:")
# for key, value in naive_cluster_lengths.items():
#     print(f'{key}: {value}')

# # Print the results
# print("\nExpert Matched Cellids:")
# for key, value in cluster_expert_dict.items():
#     print(f'{key}: {value}')

# print("\nNaive Cellids:")
# for key, value in cluster_naive_dict.items():
#     print(f'{key}: {value}')

###################################
# Remove values from expert cluster that do not have a match in naive cluster
for cluster_name, cluster in cluster_expert_dict.items():
    updated_cluster = []
    for value in cluster:
        try:
            # Find the index of the value in expert_matched_cells_final
            expert_index = expert_matched_cells_final.index(value)
            # Use the expert_index to get the corresponding value from naive_matched_cells_final
            naive_value = naive_matched_cells_final[expert_index]
            # Append the value to the updated_cluster if it has a match in naive cluster
            updated_cluster.append(value)
        except IndexError:
            # Handle the case where the index is out of range in naive_matched_cells_final
            # Skip this index and continue to the next iteration
            pass
    # Update the expert cluster with the filtered values
    cluster_expert_dict[cluster_name] = updated_cluster

# Calculate lengths of expert clusters
expert_cluster_lengths = {cluster_name: len(cluster) for cluster_name, cluster in cluster_expert_dict.items()}

# Calculate lengths of naive clusters
naive_cluster_lengths = {cluster_name: len(cluster) for cluster_name, cluster in cluster_naive_dict.items()}

# Print the lengths of expert clusters
print("Length of Expert Clusters:")
for key, value in expert_cluster_lengths.items():
    print(f'{key}: {value}')

# Print the lengths of naive clusters
print("\nLength of Naive Clusters:")
for key, value in naive_cluster_lengths.items():
    print(f'{key}: {value}')

# Print the results
print("\nExpert Matched Cellids:")
for key, value in cluster_expert_dict.items():
    print(f'{key}: {value}')

print("\nNaive Cellids:")
for key, value in cluster_naive_dict.items():
    print(f'{key}: {value}')


#%%

# expert_matched_cells_final_array = concatenated_data[expert_matched_cells_final]

# # cluster using all data to get cellIDs
# ###### cluster
# # Calculate the linkage matrix using 'ward' method
# #Z = linkage(data_transposed, method='ward')
# Z = linkage(expert_matched_cells_final_array, method='ward')

# # Set the percentage threshold for determining the number of clusters
# percentage_threshold = 0.2  # Adjust this value as needed, 0.2

# # Calculate the height threshold based on the percentage
# max_d = percentage_threshold * np.max(Z[:, 2])

# # Retrieve the clusters based on the height threshold
# clusters = fcluster(Z, max_d, criterion='distance')
# num_clusters = np.max(clusters)

# # Generate a random color for each cluster
# cmap = cm.rainbow(np.linspace(0, 1, num_clusters))
# colors = [mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap]

# # Create a color dictionary for each cluster
# color_dict = {cluster_id: color for cluster_id, color in zip(range(1, num_clusters + 1), colors)}

# # Assign colors to the dendrogram
# set_link_color_palette(colors)

# # Plot the dendrogram with colored clusters
# #plt.figure(figsize=(10, 6))
# #dendrogram(Z, color_threshold=max_d, above_threshold_color='gray')

# # Assign unique integer labels to leaf nodes
# leaf_labels = dendrogram(Z, no_plot=True)['ivl']
# label_map = {label: idx for idx, label in enumerate(leaf_labels)}
# leaf_labels_int = [label_map[label] for label in leaf_labels]

# ## Apply colors to the leaf nodes (cell IDs) based on the cluster they belong to
# #leaf_colors = [color_dict[clusters[label]] for label in leaf_labels_int]
# #ax = plt.gca()
# #ax.set_xticklabels(ax.get_xticks(), rotation=90)  # Rotate x-axis labels for better readability
# #plt.bar(range(len(leaf_labels_int)), np.zeros(len(leaf_labels_int)), color=leaf_colors, align='center', width=0.5)
# ## Plot a horizontal line at the height threshold
# #plt.axhline(y=max_d, c='red', linestyle='dashed')
# #plt.title('Hierarchical Clustering Dendrogram')
# #plt.xlabel('Cell IDs')
# #plt.ylabel('Distance')

# # Dictionary to store cluster indices
# #cluster_indices_dict = {}
# #ncells_per_cluster = []

# # Print the cell IDs for each cluster in list format and store indices
# #for cluster_id in range(1, num_clusters + 1):
# #    cluster_indices = np.where(clusters == cluster_id)[0]
# #    cell_ids = [int(leaf_labels_int[idx]) for idx in cluster_indices]
# #    cluster_name = f"cluster{cluster_id}"
# #    cluster_indices_dict[cluster_name] = cell_ids
# #    print(f"{cluster_name} = {cell_ids}")
# #    ncells_per_cluster.append(len(cell_ids))
# #cluster_expert_dict = cluster_indices_dict

# ## Plot the representative cell traces for each cluster on the same graph
# #plt.figure(figsize=(10, 6))
# #for cluster_id in range(1, num_clusters + 1):
# #    cluster_indices = np.where(clusters == cluster_id)[0]
# #    cluster_data = concatenated_data[cluster_indices]
# #    representative_trace = np.mean(cluster_data, axis=0)
# #    plt.plot(representative_trace, color=color_dict[cluster_id], label=f'Cluster {cluster_id}')
# #plt.title('Representative Cell Traces for Each Cluster')
# #plt.xlabel('Frames')
# #plt.ylabel('Amplitude')
# #plt.legend()
# #plt.show()


# ### Create separate graphs for individual cell traces in each cluster
# #for cluster_id in range(1, num_clusters + 1):
# #    cluster_indices = np.where(clusters == cluster_id)[0]
# #    cluster_data = concatenated_data[cluster_indices]
# #    plt.figure(figsize=(10, 6))
# #    for cell_trace in cluster_data:
# #        plt.plot(cell_trace, color=color_dict[cluster_id], alpha=0.3)
# #    plt.title(f'Cluster {cluster_id} - Individual Cell Traces')
# #    plt.xlabel('Frames')
# #    plt.ylabel('Amplitude')
# #    plt.show()
# # Create dictionaries to store the cluster-specific expert and naive lists
# cluster_indices_dict = {}
# ncells_per_cluster = []
# for cluster_id in range(1, num_clusters + 1):
#     cluster_indices = np.where(clusters == cluster_id)[0]
#     cell_ids = [expert_matched_cells_final[idx] for idx in cluster_indices]  # Replace with values from the original list
#     cluster_name = f"cluster{cluster_id}"
#     cluster_indices_dict[cluster_name] = cell_ids
#     #print(f"{cluster_name} = {cell_ids}")
#     ncells_per_cluster.append(len(cell_ids))
# cluster_expert_dict = cluster_indices_dict

# # #for each value in cluster_expert_dict, first find that value in the list expert_matched_cells_final, and determine its index.
# # #then use that index to find what value corresponds to the index in naive_matched_cells_final. save these values in a new dictionary called cluster_naive_dict
# ################################################################################################ OLD
# # cluster_naive_dict = {}
# # # Iterate through the clusters in cluster_expert_dict
# # for cluster_name, cluster in cluster_expert_dict.items():
# #     cluster_naive_dict[cluster_name] = []
# #     for value in cluster:
# #         # Find the index of the value in expert_matched_cells_final
# #         expert_index = expert_matched_cells_final.index(value)
# #         # Use the expert_index to get the corresponding value from naive_matched_cells_final
# #         naive_value = naive_matched_cells_final[expert_index]
# #         # Append the naive_value to the cluster_naive_dict for the current cluster
# #         cluster_naive_dict[cluster_name].append(naive_value)


# # # Print the results
# # print("Expert Matched Cellids:")
# # for key, value in cluster_expert_dict.items():
# #     print(f'{key}: {value}')

# # print("\nNaive Cellids:")
# # for key, value in cluster_naive_dict.items():
# #     print(f'{key}: {value}')
# ################################################################################################ OLD
# cluster_naive_dict = {}
# expert_cluster_lengths = {}
# naive_cluster_lengths = {}

# # Iterate through the clusters in cluster_expert_dict
# for cluster_name, cluster in cluster_expert_dict.items():
#     cluster_naive_dict[cluster_name] = []
#     expert_cluster_lengths[cluster_name] = len(cluster)
#     for value in cluster:
#         try:
#             # Find the index of the value in expert_matched_cells_final
#             expert_index = expert_matched_cells_final.index(value)
#             # Use the expert_index to get the corresponding value from naive_matched_cells_final
#             naive_value = naive_matched_cells_final[expert_index]
#             # Append the naive_value to the cluster_naive_dict for the current cluster
#             cluster_naive_dict[cluster_name].append(naive_value)
#         except IndexError:
#             # Handle the case where the index is out of range in naive_matched_cells_final
#             # Skip this index and continue to the next iteration
#             pass
#     naive_cluster_lengths[cluster_name] = len(cluster_naive_dict[cluster_name])

# # Print the lengths of expert clusters
# print("Length of Expert Clusters:")
# for key, value in expert_cluster_lengths.items():
#     print(f'{key}: {value}')

# # Print the lengths of naive clusters
# print("\nLength of Naive Clusters:")
# for key, value in naive_cluster_lengths.items():
#     print(f'{key}: {value}')

# # Print the results
# print("\nExpert Matched Cellids:")
# for key, value in cluster_expert_dict.items():
#     print(f'{key}: {value}')

# print("\nNaive Cellids:")
# for key, value in cluster_naive_dict.items():
#     print(f'{key}: {value}')





#%%
id_matches = {}  # Use a dictionary to store pairs by key

# Initialize an index to track the position in expert_matched_cells_final
index = 0

# Iterate through the dictionary and create pairs with available values
for key, values in expert_matched_cellids_dict.items():
    pairs = []
    for value in values:
        if index < len(expert_matched_cells_final):
            pairs.append((value, expert_matched_cells_final[index]))
            index += 1
    id_matches[key] = pairs

cellids_percluster_permouse = {}

for key, pairs in id_matches.items():
    cluster_matches = {}
    for cluster, values in cluster_expert_dict.items():
        first_values = [pair[0] for pair in pairs if pair[1] in values]
        if first_values:
            cluster_matches[cluster] = first_values
    cellids_percluster_permouse[key] = cluster_matches

print(cellids_percluster_permouse)


#%%
#mice_ids = ['m15n', 'm16l', 'm19n', 'm22l', 'm25lr', 'm26r', 'm30lr', 'm31r', 'm32l', 'm33r', 'm40l', 'm41r']

for m in mice_ids:
    # Check if the mouse ID exists in cellids_percluster_permouse
    if m in cellids_percluster_permouse:
      try:
        # Get the cluster1 data for the current mouse
        selected_cellids = cellids_percluster_permouse[m]['cluster8']

        # Convert the list to a string with elements separated by a space
        selected_cellids_str = '[' + ', '.join(map(str, selected_cellids)) + ']'

        # Print the formatted output for the current mouse
        #print(f'if m == \'{m}\':')
        #print(f'    selected_cellids = {selected_cellids_str}')
        print(m,len(selected_cellids))
      except:
        print('NO CELLS IN THIS CLUSTER FOR',m)
        continue
    else:
        # Handle the case where the mouse ID is not found in cellids_percluster_permouse
        print(f'Mouse {m} not found in cellids_percluster_permouse')



#%%


#%%
# Create empty arrays for the data
C_hm_expert = []
W_hm_expert = []
L_hm_expert = []
R_hm_expert = []

# Iterate through each cluster
for cluster, row_numbers in cluster_expert_dict.items():
    for row_number in row_numbers:
        if 0 <= row_number < len(concatenated_data_C_expert):
            c_data = concatenated_data_C_expert[row_number]
            C_hm_expert.append(c_data)

        if 0 <= row_number < len(concatenated_data_W_expert):
            w_data = concatenated_data_W_expert[row_number]
            W_hm_expert.append(w_data)

        if 0 <= row_number < len(concatenated_data_L_expert):
            l_data = concatenated_data_L_expert[row_number]
            L_hm_expert.append(l_data)

        if 0 <= row_number < len(concatenated_data_R_expert):
            r_data = concatenated_data_R_expert[row_number]
            R_hm_expert.append(r_data)

# Convert the lists to numpy arrays if needed
C_hm_expert = np.array(C_hm_expert)
W_hm_expert = np.array(W_hm_expert)
L_hm_expert = np.array(L_hm_expert)
R_hm_expert = np.array(R_hm_expert)

# Create empty arrays for the data
C_hm_naive = []
W_hm_naive = []
L_hm_naive = []
R_hm_naive = []

# Iterate through each cluster
for cluster, row_numbers in cluster_naive_dict.items():
    for row_number in row_numbers:
        if 0 <= row_number < len(concatenated_data_C_naive):
            c_data = concatenated_data_C_naive[row_number]
            C_hm_naive.append(c_data)

        if 0 <= row_number < len(concatenated_data_W_naive):
            w_data = concatenated_data_W_naive[row_number]
            W_hm_naive.append(w_data)

        if 0 <= row_number < len(concatenated_data_L_naive):
            l_data = concatenated_data_L_naive[row_number]
            L_hm_naive.append(l_data)

        if 0 <= row_number < len(concatenated_data_R_naive):
            r_data = concatenated_data_R_naive[row_number]
            R_hm_naive.append(r_data)

# Convert the lists to numpy arrays if needed
C_hm_naive = np.array(C_hm_naive)
W_hm_naive = np.array(W_hm_naive)
L_hm_naive = np.array(L_hm_naive)
R_hm_naive = np.array(R_hm_naive)

# Normalize the 4 arrays

# Convert the input arrays to numeric data type (assuming they are originally strings)
C_hm_expert = C_hm_expert.astype(float)
W_hm_expert = W_hm_expert.astype(float)
L_hm_expert = L_hm_expert.astype(float)
R_hm_expert = R_hm_expert.astype(float)

# Create empty arrays for the normalized data
C_hm_expert_normalized = []
W_hm_expert_normalized = []
L_hm_expert_normalized = []
R_hm_expert_normalized = []

# Iterate through each row in the arrays
for row in C_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row  # Avoid division by zero if max_val == min_val
    C_hm_expert_normalized.append(normalized_row)

for row in W_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    W_hm_expert_normalized.append(normalized_row)

for row in L_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    L_hm_expert_normalized.append(normalized_row)

for row in R_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    R_hm_expert_normalized.append(normalized_row)

# Convert the lists to numpy arrays
C_hm_expert_normalized = np.array(C_hm_expert_normalized)
W_hm_expert_normalized = np.array(W_hm_expert_normalized)
L_hm_expert_normalized = np.array(L_hm_expert_normalized)
R_hm_expert_normalized = np.array(R_hm_expert_normalized)

# Convert the input arrays to numeric data type (assuming they are originally strings)
C_hm_naive = C_hm_naive.astype(float)
W_hm_naive = W_hm_naive.astype(float)
L_hm_naive = L_hm_naive.astype(float)
R_hm_naive = R_hm_naive.astype(float)

# Create empty arrays for the normalized data
C_hm_naive_normalized = []
W_hm_naive_normalized = []
L_hm_naive_normalized = []
R_hm_naive_normalized = []

# Iterate through each row in the arrays
for row in C_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row  # Avoid division by zero if max_val == min_val
    C_hm_naive_normalized.append(normalized_row)

for row in W_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    W_hm_naive_normalized.append(normalized_row)

for row in L_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    L_hm_naive_normalized.append(normalized_row)

for row in R_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    R_hm_naive_normalized.append(normalized_row)

# Convert the lists to numpy arrays
C_hm_naive_normalized = np.array(C_hm_naive_normalized)
W_hm_naive_normalized = np.array(W_hm_naive_normalized)
L_hm_naive_normalized = np.array(L_hm_naive_normalized)
R_hm_naive_normalized = np.array(R_hm_naive_normalized)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(C_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
C_hm_expert_normalized = np.nan_to_num(C_hm_expert_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(W_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
W_hm_expert_normalized = np.nan_to_num(W_hm_expert_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(L_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
L_hm_expert_normalized = np.nan_to_num(L_hm_expert_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(R_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
R_hm_expert_normalized = np.nan_to_num(R_hm_expert_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(C_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
C_hm_naive_normalized = np.nan_to_num(C_hm_naive_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(W_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
W_hm_naive_normalized = np.nan_to_num(W_hm_naive_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(L_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
L_hm_naive_normalized = np.nan_to_num(L_hm_naive_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(R_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
R_hm_naive_normalized = np.nan_to_num(R_hm_naive_normalized, nan=0.0)

#%%
len(C_hm_naive_normalized)

#%%
len(C_hm_expert_normalized)

#%%
# correct v incorrect

#### naive
import numpy as np
import matplotlib.pyplot as plt

naive_correct_traces = []
naive_incorrect_traces = []

naive_correct_mean_traces = []
naive_incorrect_mean_traces = []

naive_correct_sem_traces = []
naive_incorrect_sem_traces = []

naive_correct_ub_traces = []
naive_correct_lb_traces = []
naive_incorrect_ub_traces = []
naive_incorrect_lb_traces = []


correct = C_hm_naive_normalized
incorrect = W_hm_naive_normalized

# for non-normalized:
#correct = C_hm_naive
#incorrect = W_hm_naive

# List of arrays to process with their names
array_names = ["correct", "incorrect"]
arrays_to_process = [correct, incorrect]

# Create a dictionary to store the cluster-specific arrays for each cluster
all_cluster_arrays = {cluster: [] for cluster in cluster_naive_dict}

# Iterate through each array
for array, array_name in zip(arrays_to_process, array_names):
    # Create a dictionary to store the cluster-specific arrays
    cluster_arrays = {}

    # Iterate through cluster_naive_dict
    for cluster, cell_ids in cluster_naive_dict.items():
        # Get the length of the current cluster
        cluster_length = len(cell_ids)

        # Slice the current array
        cluster_array = array[:cluster_length]

        # Store the cluster-specific array in the dictionary
        cluster_arrays[cluster] = cluster_array

        # Update the current array to exclude the sliced rows
        array = array[cluster_length:]

    for cluster, cluster_array in cluster_arrays.items():
      if array_name == "correct":
          naive_correct_traces.append(cluster_array)
      elif array_name == "incorrect":
          naive_incorrect_traces.append(cluster_array)

    # Append the cluster-specific arrays to the corresponding entry in all_cluster_arrays
    for cluster, cluster_array in cluster_arrays.items():
        all_cluster_arrays[cluster].append((array_name, cluster_array))

for i in range(len(naive_correct_traces)):
  mean = np.mean(naive_correct_traces[i], axis=0)
  sem = np.std(naive_correct_traces[i], axis=0) / np.sqrt(len(naive_correct_traces[i]))
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  naive_correct_mean_traces.append(mean)
  naive_correct_sem_traces.append(sem)
  naive_correct_ub_traces.append(upper_bound)
  naive_correct_lb_traces.append(lower_bound)

for i in range(len(naive_incorrect_traces)):
  mean = np.mean(naive_incorrect_traces[i], axis=0)
  sem = np.std(naive_incorrect_traces[i], axis=0) / np.sqrt(len(naive_incorrect_traces[i]))
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  naive_incorrect_mean_traces.append(mean)
  naive_incorrect_sem_traces.append(sem)
  naive_incorrect_ub_traces.append(upper_bound)
  naive_incorrect_lb_traces.append(lower_bound)

#### expert
expert_correct_traces = []
expert_incorrect_traces = []

expert_correct_mean_traces = []
expert_incorrect_mean_traces = []

expert_correct_sem_traces = []
expert_incorrect_sem_traces = []

expert_correct_ub_traces = []
expert_correct_lb_traces = []
expert_incorrect_ub_traces = []
expert_incorrect_lb_traces = []


correct = C_hm_expert_normalized
incorrect = W_hm_expert_normalized

# for non-normalized:
#correct = C_hm_expert
#incorrect = W_hm_expert

# List of arrays to process with their names
array_names = ["correct", "incorrect"]
arrays_to_process = [correct, incorrect]

fig_width, fig_height = 4, 3  # You can adjust these values as needed

# Create a dictionary to store the cluster-specific arrays for each cluster
all_cluster_arrays = {cluster: [] for cluster in cluster_expert_dict}

# Iterate through each array
for array, array_name in zip(arrays_to_process, array_names):
    # Create a dictionary to store the cluster-specific arrays
    cluster_arrays = {}

    # Iterate through cluster_expert_dict
    for cluster, cell_ids in cluster_expert_dict.items():
        # Get the length of the current cluster
        cluster_length = len(cell_ids)

        # Slice the current array
        cluster_array = array[:cluster_length]

        # Store the cluster-specific array in the dictionary
        cluster_arrays[cluster] = cluster_array

        # Update the current array to exclude the sliced rows
        array = array[cluster_length:]

    for cluster, cluster_array in cluster_arrays.items():
      if array_name == "correct":
          expert_correct_traces.append(cluster_array)
      elif array_name == "incorrect":
          expert_incorrect_traces.append(cluster_array)

    # Append the cluster-specific arrays to the corresponding entry in all_cluster_arrays
    for cluster, cluster_array in cluster_arrays.items():
        all_cluster_arrays[cluster].append((array_name, cluster_array))

for i in range(len(expert_correct_traces)):
  mean = np.mean(expert_correct_traces[i],axis=0)
  sem = np.std(expert_correct_traces[i], axis=0) / np.sqrt(len(expert_correct_traces[i]))
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  expert_correct_mean_traces.append(mean)
  expert_correct_sem_traces.append(sem)
  expert_correct_ub_traces.append(upper_bound)
  expert_correct_lb_traces.append(lower_bound)

for i in range(len(expert_incorrect_traces)):
  mean = np.mean(expert_incorrect_traces[i],axis=0)
  sem = np.std(expert_incorrect_traces[i], axis=0) / np.sqrt(len(expert_incorrect_traces[i]))
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  expert_incorrect_mean_traces.append(mean)
  expert_incorrect_sem_traces.append(sem)
  expert_incorrect_ub_traces.append(upper_bound)
  expert_incorrect_lb_traces.append(lower_bound)

#### plot all
# Define colors with different alpha values for different shades of blue
dark_blue = (0, 0, 0.5)  # Dark blue
light_blue = (0, 0, 1, 0.5)  # Light blue with alpha (transparency)
dark_red = (1.0, 0.0, 0.0)  # Dark red (RGB format)
light_red = (1.0, 0.2, 0.4, 0.7)  # Light red with alpha (RGBA format), use blue and red instead

for i in range(len(naive_correct_mean_traces)):
  plt.figure(figsize=(4, 3))
  plt.plot(naive_correct_mean_traces[i], color=dark_red, label='naive_correct')
  plt.fill_between(range(len(naive_correct_mean_traces[i])), naive_correct_lb_traces[i], naive_correct_ub_traces[i], alpha=0.2, color=dark_red)

  plt.plot(naive_incorrect_mean_traces[i], color=light_red, label='naive_incorrect')
  plt.fill_between(range(len(naive_incorrect_mean_traces[i])), naive_incorrect_lb_traces[i], naive_incorrect_ub_traces[i], alpha=0.2, color=light_red)

  plt.plot(expert_correct_mean_traces[i], color=dark_blue, label='expert_correct')
  plt.fill_between(range(len(expert_correct_mean_traces[i])), expert_correct_lb_traces[i], expert_correct_ub_traces[i], alpha=0.2, color=dark_blue)

  plt.plot(expert_incorrect_mean_traces[i], color=light_blue, label='expert_incorrect')
  plt.fill_between(range(len(expert_incorrect_mean_traces[i])), expert_incorrect_lb_traces[i], expert_incorrect_ub_traces[i], alpha=0.2, color=light_blue)


  # Add labels and a legend for the current cluster
  plt.xlabel('Frames')
  plt.ylabel('Average Fluorescence Intensity')
  plt.legend()

  plt.axvline(x=16, color='black', linestyle='--', label='turnframe')
  plt.legend(fontsize='small')

  # Set the title to indicate the cluster
  plt.title('cluster'+f'{i+1}')

  try:
    ### TO SAVE OUTPUTS
    ## Generate the file name dynamically based on the cluster index
    figure_path = r"C:\Users\aozgu\Downloads\cluster" + f'{i+1}' + ".svg"
    #figure_path = r"C:\Users\LurLab\Downloads\cluster" + f'{i+1}' + ".svg"
    ## Save the figure
    plt.savefig(figure_path, format='svg', bbox_inches='tight')
  except:
    pass

  plt.show()


#%%

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# Lists to store the results
expert_scores = []
naive_scores = []

expert_scores_preturn = []
expert_scores_postturn =[]

naive_scores_preturn = []
naive_scores_postturn =[]

correct_scores = []
incorrect_scores = []

# Calculate and store scores for correct traces, all frames
for trace1, trace2 in zip(expert_correct_mean_traces, naive_correct_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    correct_scores.append(result)

# Calculate and store scores for incorrect traces, all frames
for trace1, trace2 in zip(expert_incorrect_mean_traces, naive_incorrect_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    incorrect_scores.append(result)

# Calculate 95% confidence intervals
ci_correct = 1.96 * np.std(correct_scores) / np.sqrt(len(correct_scores))
ci_incorrect = 1.96 * np.std(incorrect_scores) / np.sqrt(len(incorrect_scores))

# Plotting as a bar graph with error bars
groups = ['Correct', 'Incorrect']
average_scores = [np.mean(correct_scores), np.mean(incorrect_scores)]
print(average_scores)
confidence_intervals = [ci_correct, ci_incorrect]
print(confidence_intervals)

plt.bar(groups, average_scores, yerr=confidence_intervals, capsize=10, color=['blue', 'orange'])
plt.xlabel('Group')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences (Across Traces)')

# Set the y-axis to range from 0 to 0.25
plt.ylim(0, 0.3)

# Perform t-test
t_stat, p_value = ttest_ind(correct_scores, incorrect_scores)
print(f'T-test p-value: {p_value}')
#
if p_value < 0.05:
    print('The difference between the groups is statistically significant.')
else:
    print('The difference between the groups is not statistically significant.')

### TO SAVE OUTPUTS
## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\groupedcluster_differences" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\groupedcluster_differences" + ".eps"
## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()


#%%
# left v right

#### naive
import numpy as np
import matplotlib.pyplot as plt

naive_left_traces = []
naive_right_traces = []

naive_left_mean_traces = []
naive_right_mean_traces = []

naive_left_sem_traces = []
naive_right_sem_traces = []

naive_left_ub_traces = []
naive_left_lb_traces = []
naive_right_ub_traces = []
naive_right_lb_traces = []

left = L_hm_naive_normalized
right = R_hm_naive_normalized

# List of arrays to process with their names
array_names = ["left", "right"]
arrays_to_process = [left, right]

# Create a dictionary to store the cluster-specific arrays for each cluster
all_cluster_arrays = {cluster: [] for cluster in cluster_naive_dict}

# Iterate through each array
for array, array_name in zip(arrays_to_process, array_names):
    # Create a dictionary to store the cluster-specific arrays
    cluster_arrays = {}

    # Iterate through cluster_naive_dict
    for cluster, cell_ids in cluster_naive_dict.items():
        # Get the length of the current cluster
        cluster_length = len(cell_ids)

        # Slice the current array
        cluster_array = array[:cluster_length]

        # Store the cluster-specific array in the dictionary
        cluster_arrays[cluster] = cluster_array

        # Update the current array to exclude the sliced rows
        array = array[cluster_length:]

    for cluster, cluster_array in cluster_arrays.items():
        if array_name == "left":
            naive_left_traces.append(cluster_array)
        elif array_name == "right":
            naive_right_traces.append(cluster_array)

    # Append the cluster-specific arrays to the corresponding entry in all_cluster_arrays
    for cluster, cluster_array in cluster_arrays.items():
        all_cluster_arrays[cluster].append((array_name, cluster_array))

for i in range(len(naive_left_traces)):
    mean = np.mean(naive_left_traces[i], axis=0)
    sem = np.std(naive_left_traces[i], axis=0) / np.sqrt(len(naive_left_traces[i]))
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    naive_left_mean_traces.append(mean)
    naive_left_sem_traces.append(sem)
    naive_left_ub_traces.append(upper_bound)
    naive_left_lb_traces.append(lower_bound)

for i in range(len(naive_right_traces)):
    mean = np.mean(naive_right_traces[i], axis=0)
    sem = np.std(naive_right_traces[i], axis=0) / np.sqrt(len(naive_right_traces[i]))
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    naive_right_mean_traces.append(mean)
    naive_right_sem_traces.append(sem)
    naive_right_ub_traces.append(upper_bound)
    naive_right_lb_traces.append(lower_bound)

#### expert
expert_left_traces = []
expert_right_traces = []

expert_left_mean_traces = []
expert_right_mean_traces = []

expert_left_sem_traces = []
expert_right_sem_traces = []

expert_left_ub_traces = []
expert_left_lb_traces = []
expert_right_ub_traces = []
expert_right_lb_traces = []

left = L_hm_expert_normalized
right = R_hm_expert_normalized

# List of arrays to process with their names
array_names = ["left", "right"]
arrays_to_process = [left, right]

fig_width, fig_height = 4, 3  # You can adjust these values as needed

# Create a dictionary to store the cluster-specific arrays for each cluster
all_cluster_arrays = {cluster: [] for cluster in cluster_expert_dict}

# Iterate through each array
for array, array_name in zip(arrays_to_process, array_names):
    # Create a dictionary to store the cluster-specific arrays
    cluster_arrays = {}

    # Iterate through cluster_expert_dict
    for cluster, cell_ids in cluster_expert_dict.items():
        # Get the length of the current cluster
        cluster_length = len(cell_ids)

        # Slice the current array
        cluster_array = array[:cluster_length]

        # Store the cluster-specific array in the dictionary
        cluster_arrays[cluster] = cluster_array

        # Update the current array to exclude the sliced rows
        array = array[cluster_length:]

    for cluster, cluster_array in cluster_arrays.items():
        if array_name == "left":
            expert_left_traces.append(cluster_array)
        elif array_name == "right":
            expert_right_traces.append(cluster_array)

    # Append the cluster-specific arrays to the corresponding entry in all_cluster_arrays
    for cluster, cluster_array in cluster_arrays.items():
        all_cluster_arrays[cluster].append((array_name, cluster_array))

for i in range(len(expert_left_traces)):
    mean = np.mean(expert_left_traces[i], axis=0)
    sem = np.std(expert_left_traces[i], axis=0) / np.sqrt(len(expert_left_traces[i]))
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    expert_left_mean_traces.append(mean)
    expert_left_sem_traces.append(sem)
    expert_left_ub_traces.append(upper_bound)
    expert_left_lb_traces.append(lower_bound)

for i in range(len(expert_right_traces)):
    mean = np.mean(expert_right_traces[i], axis=0)
    sem = np.std(expert_right_traces[i], axis=0) / np.sqrt(len(expert_right_traces[i]))
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    expert_right_mean_traces.append(mean)
    expert_right_sem_traces.append(sem)
    expert_right_ub_traces.append(upper_bound)
    expert_right_lb_traces.append(lower_bound)

#### plot all
# Define colors with different alpha values for different shades of blue
dark_blue = (0, 0, 0.5)  # Dark blue
light_blue = (0, 0, 1, 0.5)  # Light blue with alpha (transparency)
dark_red = (1.0, 0.0, 0.0)  # Dark red (RGB format)
light_red = (1.0, 0.2, 0.4, 0.7)  # Light red with alpha (RGBA format), use blue and red instead


for i in range(len(naive_left_mean_traces)):
    plt.figure(figsize=(4, 3))
    plt.plot(naive_left_mean_traces[i], color=dark_red, label='naive_left')
    plt.fill_between(range(len(naive_left_mean_traces[i])), naive_left_lb_traces[i], naive_left_ub_traces[i], alpha=0.2, color=dark_red)

    plt.plot(naive_right_mean_traces[i], color=light_red, label='naive_right')
    plt.fill_between(range(len(naive_right_mean_traces[i])), naive_right_lb_traces[i], naive_right_ub_traces[i], alpha=0.2, color=light_red)

    plt.plot(expert_left_mean_traces[i], color=dark_blue, label='expert_left')
    plt.fill_between(range(len(expert_left_mean_traces[i])), expert_left_lb_traces[i], expert_left_ub_traces[i], alpha=0.2, color=dark_blue)

    plt.plot(expert_right_mean_traces[i], color=light_blue, label='expert_right')
    plt.fill_between(range(len(expert_right_mean_traces[i])), expert_right_lb_traces[i], expert_right_ub_traces[i], alpha=0.2, color=light_blue)


    # Add labels and a legend for the current cluster
    plt.xlabel('Frames')
    plt.ylabel('Average Fluorescence Intensity')
    plt.legend()

    plt.axvline(x=16, color='black', linestyle='--', label='turnframe')
    plt.legend(fontsize='small')

    # Set the title to indicate the cluster
    plt.title('cluster'+f'{i+1}')

    try:
      ### TO SAVE OUTPUTS
      ## Generate the file name dynamically based on the cluster index
      figure_path = r"C:\Users\aozgu\Downloads\cluster" + f'{i+1}' + ".svg"
      #figure_path = r"C:\Users\LurLab\Downloads\cluster" + f'{i+1}' + ".svg"
      ## Save the figure
      plt.savefig(figure_path, format='svg', bbox_inches='tight')
    except:
      pass

    plt.show()


#%%

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# Lists to store the results
expert_scores = []
naive_scores = []

expert_scores_preturn = []
expert_scores_postturn =[]

naive_scores_preturn = []
naive_scores_postturn =[]

left_scores = []
right_scores = []

# Calculate and store scores for correct traces, all frames
for trace1, trace2 in zip(expert_left_mean_traces, naive_left_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    left_scores.append(result)

# Calculate and store scores for incorrect traces, all frames
for trace1, trace2 in zip(expert_right_mean_traces, naive_right_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    right_scores.append(result)

# Calculate 95% confidence intervals
ci_left = 1.96 * np.std(left_scores) / np.sqrt(len(left_scores))
ci_right = 1.96 * np.std(right_scores) / np.sqrt(len(right_scores))

# Plotting as a bar graph with error bars
groups = ['Left turn', 'Right turn']
average_scores = [np.mean(left_scores), np.mean(right_scores)]
confidence_intervals = [ci_left, ci_right]

plt.bar(groups, average_scores, yerr=confidence_intervals, capsize=10, color=['blue', 'orange'])
plt.xlabel('Group')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences (Across Traces)')

# Set the y-axis to range from 0 to 0.25
plt.ylim(0, 0.3)

# Perform t-test
t_stat, p_value = ttest_ind(left_scores, right_scores)
print(f'T-test p-value: {p_value}')
#
if p_value < 0.05:
    print('The difference between the groups is statistically significant.')
else:
    print('The difference between the groups is not statistically significant.')

### TO SAVE OUTPUTS
## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\groupedcluster_differences" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\groupedcluster_differences" + ".eps"
## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()


#%%
import matplotlib.pyplot as plt

# Create a figure with subplots
fig, axs = plt.subplots(2, 4, figsize=(16, 8))  # 2 rows, 4 columns

# Plot the expert heatmaps with aspect ratio 'auto' and dashed lines
heatmaps_expert = [C_hm_expert_normalized, W_hm_expert_normalized, L_hm_expert_normalized, R_hm_expert_normalized]
titles_expert = ['correct_expert', 'incorrect_expert', 'left_expert', 'right_expert']

for i, (heatmap, title) in enumerate(zip(heatmaps_expert, titles_expert)):
    row = 0  # Plot in the first row
    cax = axs[row, i].matshow(heatmap, cmap='viridis', aspect='auto')

    # Add colorbar for the current heatmap
    #fig.colorbar(cax, ax=axs[row, i])

    # Set title for the current heatmap
    axs[row, i].set_title(title)

    # Iterate through cluster lengths to draw dashed lines for expert data
    rows_drawn = 0
    for cluster_length in cluster_expert_dict.values():
        cluster_length = len(cluster_length)
        if rows_drawn + cluster_length < heatmap.shape[0]:
            axs[row, i].axhline(rows_drawn + cluster_length - 0.5, color='white', linestyle='dashed', linewidth=1)
        rows_drawn += cluster_length

# Plot the naive heatmaps with aspect ratio 'auto' and dashed lines
heatmaps_naive = [C_hm_naive_normalized, W_hm_naive_normalized, L_hm_naive_normalized, R_hm_naive_normalized]
titles_naive = ['correct_naive', 'incorrect_naive', 'left_naive', 'right_naive']

for i, (heatmap, title) in enumerate(zip(heatmaps_naive, titles_naive)):
    row = 1  # Plot in the second row
    cax = axs[row, i].matshow(heatmap, cmap='viridis', aspect='auto')

    # Add colorbar for the current heatmap
    #fig.colorbar(cax, ax=axs[row, i])

    # Set title for the current heatmap
    axs[row, i].set_title(title)

    # Iterate through cluster lengths to draw dashed lines for naive data
    rows_drawn = 0
    for cluster_length in cluster_naive_dict.values():
        cluster_length = len(cluster_length)
        if rows_drawn + cluster_length < heatmap.shape[0]:
            axs[row, i].axhline(rows_drawn + cluster_length - 0.5, color='white', linestyle='dashed', linewidth=1)
        rows_drawn += cluster_length

# Adjust spacing between subplots
plt.tight_layout()

# Show the plot
plt.show()


#%%
import numpy as np
import matplotlib.pyplot as plt

# Assuming heatmaps_expert and heatmaps_naive are provided as complete lists of numpy arrays

# Define the function to sort the rows of a heatmap and return the order
def sort_rows_by_peak_activity(heatmap):
    max_indices = np.argmax(heatmap, axis=1)
    sorted_order = np.argsort(max_indices)
    return sorted_order

# Define the function to rearrange the rows of a naive heatmap based on the expert heatmap's order
def rearrange_naive_to_expert_order(naive_heatmap, expert_order):
    rearranged_naive_heatmap = naive_heatmap[expert_order]
    return rearranged_naive_heatmap

# Get the sorted orders for the expert heatmaps
sorted_orders_expert = [sort_rows_by_peak_activity(hm) for hm in heatmaps_expert]
sorted_heatmaps_expert = [heatmap[order] for heatmap, order in zip(heatmaps_expert, sorted_orders_expert)]

# Apply the rearrangement for each naive heatmap based on its corresponding expert heatmap's sorted order
rearranged_heatmaps_naive = []
for naive_heatmap, expert_order in zip(heatmaps_naive, sorted_orders_expert):
    rearranged_naive_heatmap = rearrange_naive_to_expert_order(naive_heatmap, expert_order)
    rearranged_heatmaps_naive.append(rearranged_naive_heatmap)

# Create figure with subplots for plotting
fig, axs = plt.subplots(2, 4, figsize=(20, 10))  # 2 rows, 4 columns for expert and naive heatmaps

# Plot sorted expert heatmaps
for i, heatmap in enumerate(sorted_heatmaps_expert):
    axs[0, i].imshow(heatmap, aspect='auto', cmap='viridis')
    axs[0, i].set_title(f'expert_{i+1}')  # Or use your specific titles

# Plot rearranged naive heatmaps
for i, heatmap in enumerate(rearranged_heatmaps_naive):
    axs[1, i].imshow(heatmap, aspect='auto', cmap='viridis')
    axs[1, i].set_title(f'naive_{i+1}')  # Or use your specific titles

# Tight layout for better spacing
plt.tight_layout()

### TO SAVE OUTPUTS
## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\heatmaps" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\groupedcluster_differences" + ".eps"
## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()


#%%
sorted_orders_expert

#%%
# Print out the rearranged order for naive
for i, expert_order in enumerate(sorted_orders_expert):
    print(f'Rearranged order for naive_{i+1}: {expert_order}')


#%%
sorted_heatmaps_naive

#%%
len(L_hm_expert_normalized)

#%%
len(L_hm_naive_normalized)

#%%
print("Keys of sorted_cluster_expert_dict:", sorted_cluster_expert_dict.keys())


#%%
print("Keys of sorted_cluster_expert_dict:", sorted_cluster_expert_dict.keys())
print("Keys of sorted_cluster_naive_dict:", sorted_cluster_naive_dict.keys())


#%%


#%%


#%%


#%%
#clustering
### FUNCTIONS
def get_parameters(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    #if want only correct trials
    #clean_df = clean_df[clean_df['correctness']=='correct'] ######## remove

    #if want only left visual trials
    #clean_df = clean_df[clean_df['direction']=='L']

    #cluster by left turn
    #clean_df = clean_df[(clean_df['correctness'] == 'correct') & (clean_df['direction'] == 'L') |
    #                   (clean_df['correctness'] == 'wrong') & (clean_df['direction'] == 'R')]

    ###del when not using all trials
    #clean_df = params
    #a = 0

    #only mixstage
    if init_only == 1:
        clean_df = clean_df[clean_df['mode']=='INIT']

    clean_trial_start_frames = []
    clean_trial_stim_start_frames = []
    clean_trial_end_frames = []
    clean_trial_turn_frames = []
    clean_trial_correctness = []
    clean_trial_direction = []
    count_lc = []
    count_lw = []
    count_rc = []
    count_rw = []

    for i in range(len(clean_df)):
      if a == 1:
        clean_trial_stim_start_frames.append(int(clean_df['stim_start'].iloc[i])) ###      2
        clean_trial_turn_frames.append(int(clean_df['turn_frame'].iloc[i]))          #######       3
      clean_trial_end_frames.append(int(clean_df['trial_end'].iloc[i]))
      clean_trial_start_frames.append(int(clean_df['trial_start'].iloc[i]))

      if clean_df['correctness'].iloc[i]=='correct':
        clean_trial_correctness.append(1)
      else:
        clean_trial_correctness.append(0)

    ## use for turn decoding
      #if clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='correct':
      #  count_lc.append(1)
      #  clean_trial_direction.append(0)
      #elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='wrong':
      #  count_rw.append(1)
      #  clean_trial_direction.append(0)
      #elif clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='wrong':
      #  count_lw.append(1)
      #  clean_trial_direction.append(1)
      #elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='correct':
      #  count_rc.append(1)
      #  clean_trial_direction.append(1)

      ## use only for visual stimuli analysis
      if clean_df['direction'].iloc[i]=='L':
        clean_trial_direction.append(0)
      if clean_df['direction'].iloc[i]=='R':
        clean_trial_direction.append(1)


    ##### MAIN BELOW                                    5
    #ntrials = len(clean_trial_turn_frames)
    ntrials = len(clean_trial_end_frames)

    y_all = []
    if trial_type == 1:
        for i in range(ntrials):
            y_all.append(clean_trial_correctness[i])
    else:
        for i in range(ntrials):
            y_all.append(clean_trial_direction[i])

    num_incorrect_trials = y_all.count(0)

    if a == 1:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_stim_start_frames': clean_trial_stim_start_frames, ###6
      'clean_trial_end_frames': clean_trial_end_frames,
      'clean_trial_turn_frames': clean_trial_turn_frames,       ####7
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }
    else:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_end_frames': clean_trial_end_frames,
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }

    return base_params

def get_trace_average_per_cell(cal,base_params):
    #from scipy.stats import zscore
    t_main = []
    for i in range(base_params['ntrials']):

                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-14:base_params['clean_trial_turn_frames'][i]])
                #tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]-30:base_params['clean_trial_end_frames'][i]-3])
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_turn_frames'][i]])
                tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-15:base_params['clean_trial_turn_frames'][i]+25]) #MAIN
                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-20:base_params['clean_trial_turn_frames'][i]+5])
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_stim_start_frames'][i]+20])

                tframes = np.asarray(tframes)

                t_main.append(tframes)

    t_avg = np.mean(t_main,axis=0)
    #t_avg = np.nanmean(t_main, axis=0)
    return t_avg

### FUNCTIONS
def get_parameters_C(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    #if want only correct trials
    clean_df = clean_df[clean_df['correctness']=='correct']

    #if want to specify direction
    #clean_df = clean_df[clean_df['direction']=='L']


    #only mixstage
    if init_only == 1:
        clean_df = clean_df[clean_df['mode']=='INIT']

    clean_trial_start_frames = []
    clean_trial_stim_start_frames = []
    clean_trial_end_frames = []
    clean_trial_turn_frames = []
    clean_trial_correctness = []
    clean_trial_direction = []
    count_lc = []
    count_lw = []
    count_rc = []
    count_rw = []

    for i in range(len(clean_df)):
      if a == 1:
        clean_trial_stim_start_frames.append(int(clean_df['stim_start'].iloc[i])) ###      2
        clean_trial_turn_frames.append(int(clean_df['turn_frame'].iloc[i]))          #######       3
      clean_trial_end_frames.append(int(clean_df['trial_end'].iloc[i]))
      clean_trial_start_frames.append(int(clean_df['trial_start'].iloc[i]))

      if clean_df['correctness'].iloc[i]=='correct':
        clean_trial_correctness.append(1)
      else:
        clean_trial_correctness.append(0)

    ## use for turn decoding
      if clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='correct':
        count_lc.append(1)
        clean_trial_direction.append(0)
      elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='wrong':
        count_rw.append(1)
        clean_trial_direction.append(0)
      elif clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='wrong':
        count_lw.append(1)
        clean_trial_direction.append(1)
      elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='correct':
        count_rc.append(1)
        clean_trial_direction.append(1)

      ## use only for visual stimuli analysis
      #if clean_df['direction'].iloc[i]=='L':
      # clean_trial_direction.append(0)
      #if clean_df['direction'].iloc[i]=='R':
      #  clean_trial_direction.append(1)

        #del later     4
        #ntrials = len(clean_trial_start_frames)


    ##### MAIN BELOW                                    5
    #ntrials = len(clean_trial_turn_frames)
    ntrials = len(clean_trial_end_frames)


    y_all = []
    if trial_type == 1:
        for i in range(ntrials):
            y_all.append(clean_trial_correctness[i])
    else:
        for i in range(ntrials):
            y_all.append(clean_trial_direction[i])

    num_incorrect_trials = y_all.count(0)

    if a == 1:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_stim_start_frames': clean_trial_stim_start_frames, ###6
      'clean_trial_end_frames': clean_trial_end_frames,
      'clean_trial_turn_frames': clean_trial_turn_frames,       ####7
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }
    else:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_end_frames': clean_trial_end_frames,
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }

    return base_params

### FUNCTIONS
def get_parameters_W(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    #if want only correct trials
    clean_df = clean_df[clean_df['correctness']=='wrong']

    #if want to specify direction
    #clean_df = clean_df[clean_df['direction']=='L']


    #only mixstage
    if init_only == 1:
        clean_df = clean_df[clean_df['mode']=='INIT']

    clean_trial_start_frames = []
    clean_trial_stim_start_frames = []
    clean_trial_end_frames = []
    clean_trial_turn_frames = []
    clean_trial_correctness = []
    clean_trial_direction = []
    count_lc = []
    count_lw = []
    count_rc = []
    count_rw = []

    for i in range(len(clean_df)):
      if a == 1:
        clean_trial_stim_start_frames.append(int(clean_df['stim_start'].iloc[i])) ###      2
        clean_trial_turn_frames.append(int(clean_df['turn_frame'].iloc[i]))          #######       3
      clean_trial_end_frames.append(int(clean_df['trial_end'].iloc[i]))
      clean_trial_start_frames.append(int(clean_df['trial_start'].iloc[i]))

      if clean_df['correctness'].iloc[i]=='correct':
        clean_trial_correctness.append(1)
      else:
        clean_trial_correctness.append(0)

    ## use for turn decoding
      if clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='correct':
        count_lc.append(1)
        clean_trial_direction.append(0)
      elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='wrong':
        count_rw.append(1)
        clean_trial_direction.append(0)
      elif clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='wrong':
        count_lw.append(1)
        clean_trial_direction.append(1)
      elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='correct':
        count_rc.append(1)
        clean_trial_direction.append(1)

      ## use only for visual stimuli analysis
      #if clean_df['direction'].iloc[i]=='L':
      # clean_trial_direction.append(0)
      #if clean_df['direction'].iloc[i]=='R':
      #  clean_trial_direction.append(1)

        #del later     4
        #ntrials = len(clean_trial_start_frames)


    ##### MAIN BELOW                                    5
    #ntrials = len(clean_trial_turn_frames)
    ntrials = len(clean_trial_end_frames)


    y_all = []
    if trial_type == 1:
        for i in range(ntrials):
            y_all.append(clean_trial_correctness[i])
    else:
        for i in range(ntrials):
            y_all.append(clean_trial_direction[i])

    num_incorrect_trials = y_all.count(0)

    if a == 1:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_stim_start_frames': clean_trial_stim_start_frames, ###6
      'clean_trial_end_frames': clean_trial_end_frames,
      'clean_trial_turn_frames': clean_trial_turn_frames,       ####7
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }
    else:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_end_frames': clean_trial_end_frames,
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }

    return base_params

### FUNCTIONS
def get_parameters_L(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    #left turn
    clean_df = clean_df[(clean_df['correctness'] == 'correct') & (clean_df['direction'] == 'L') |
                       (clean_df['correctness'] == 'wrong') & (clean_df['direction'] == 'R')]

    #left vis stim
    #clean_df = clean_df[clean_df['direction'] == 'L']

    #only mixstage
    if init_only == 1:
        clean_df = clean_df[clean_df['mode']=='INIT']

    clean_trial_start_frames = []
    clean_trial_stim_start_frames = []
    clean_trial_end_frames = []
    clean_trial_turn_frames = []
    clean_trial_correctness = []
    clean_trial_direction = []
    count_lc = []
    count_lw = []
    count_rc = []
    count_rw = []

    for i in range(len(clean_df)):
      if a == 1:
        clean_trial_stim_start_frames.append(int(clean_df['stim_start'].iloc[i])) ###      2
        clean_trial_turn_frames.append(int(clean_df['turn_frame'].iloc[i]))          #######       3
      clean_trial_end_frames.append(int(clean_df['trial_end'].iloc[i]))
      clean_trial_start_frames.append(int(clean_df['trial_start'].iloc[i]))

      if clean_df['correctness'].iloc[i]=='correct':
        clean_trial_correctness.append(1)
      else:
        clean_trial_correctness.append(0)

    ## use for turn decoding
      #if clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='correct':
      #  count_lc.append(1)
      #  clean_trial_direction.append(0)
      #elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='wrong':
      #  count_rw.append(1)
      #  clean_trial_direction.append(0)
      #elif clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='wrong':
      #  count_lw.append(1)
      #  clean_trial_direction.append(1)
      #elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='correct':
      #  count_rc.append(1)
      #  clean_trial_direction.append(1)

      ## use only for visual stimuli analysis
      if clean_df['direction'].iloc[i]=='L':
       clean_trial_direction.append(0)
      if clean_df['direction'].iloc[i]=='R':
        clean_trial_direction.append(1)

        #del later     4
        #ntrials = len(clean_trial_start_frames)


    ##### MAIN BELOW                                    5
    #ntrials = len(clean_trial_turn_frames)
    ntrials = len(clean_trial_end_frames)


    y_all = []
    if trial_type == 1:
        for i in range(ntrials):
            y_all.append(clean_trial_correctness[i])
    else:
        for i in range(ntrials):
            y_all.append(clean_trial_direction[i])

    num_incorrect_trials = y_all.count(0)

    if a == 1:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_stim_start_frames': clean_trial_stim_start_frames, ###6
      'clean_trial_end_frames': clean_trial_end_frames,
      'clean_trial_turn_frames': clean_trial_turn_frames,       ####7
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }
    else:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_end_frames': clean_trial_end_frames,
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }

    return base_params

### FUNCTIONS
def get_parameters_R(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    # right turn
    clean_df = clean_df[(clean_df['correctness'] == 'correct') & (clean_df['direction'] == 'R') |
                       (clean_df['correctness'] == 'wrong') & (clean_df['direction'] == 'L')]

    #right vis stim
    clean_df = clean_df[clean_df['direction'] == 'R']

    #only mixstage
    if init_only == 1:
        clean_df = clean_df[clean_df['mode']=='INIT']

    clean_trial_start_frames = []
    clean_trial_stim_start_frames = []
    clean_trial_end_frames = []
    clean_trial_turn_frames = []
    clean_trial_correctness = []
    clean_trial_direction = []
    count_lc = []
    count_lw = []
    count_rc = []
    count_rw = []

    for i in range(len(clean_df)):
      if a == 1:
        clean_trial_stim_start_frames.append(int(clean_df['stim_start'].iloc[i])) ###      2
        clean_trial_turn_frames.append(int(clean_df['turn_frame'].iloc[i]))          #######       3
      clean_trial_end_frames.append(int(clean_df['trial_end'].iloc[i]))
      clean_trial_start_frames.append(int(clean_df['trial_start'].iloc[i]))

      if clean_df['correctness'].iloc[i]=='correct':
        clean_trial_correctness.append(1)
      else:
        clean_trial_correctness.append(0)

    ## use for turn decoding
      #if clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='correct':
      #  count_lc.append(1)
      #  clean_trial_direction.append(0)
      #elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='wrong':
      #  count_rw.append(1)
      #  clean_trial_direction.append(0)
      #elif clean_df['direction'].iloc[i]=='L' and clean_df['correctness'].iloc[i]=='wrong':
      #  count_lw.append(1)
      #  clean_trial_direction.append(1)
      #elif clean_df['direction'].iloc[i]=='R' and clean_df['correctness'].iloc[i]=='correct':
      #  count_rc.append(1)
      #  clean_trial_direction.append(1)

      ## use only for visual stimuli analysis
      if clean_df['direction'].iloc[i]=='L':
       clean_trial_direction.append(0)
      if clean_df['direction'].iloc[i]=='R':
        clean_trial_direction.append(1)

        #del later     4
        #ntrials = len(clean_trial_start_frames)


    ##### MAIN BELOW                                    5
    #ntrials = len(clean_trial_turn_frames)
    ntrials = len(clean_trial_end_frames)


    y_all = []
    if trial_type == 1:
        for i in range(ntrials):
            y_all.append(clean_trial_correctness[i])
    else:
        for i in range(ntrials):
            y_all.append(clean_trial_direction[i])

    num_incorrect_trials = y_all.count(0)

    if a == 1:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_stim_start_frames': clean_trial_stim_start_frames, ###6
      'clean_trial_end_frames': clean_trial_end_frames,
      'clean_trial_turn_frames': clean_trial_turn_frames,       ####7
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }
    else:
      base_params = {
      'clean_trial_start_frames': clean_trial_start_frames,
      'clean_trial_end_frames': clean_trial_end_frames,
      'ntrials':ntrials,
      'y_all':y_all,
      'num_incorrect_trials':num_incorrect_trials,
      'clean_trial_direction': clean_trial_direction,
      'count_lc':len(count_lc),
      'count_lw':len(count_lw),
      'count_rc':len(count_rc),
      'count_rw':len(count_rw)
      }

    return base_params

