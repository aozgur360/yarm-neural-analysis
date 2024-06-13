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
def average_absolute_differences(trace1, trace2):
    # Check if the traces have the same length
    if len(trace1) != len(trace2):
        raise ValueError("Traces must have the same length")

    # Calculate absolute differences for each index
    differences = [abs(trace1[i] - trace2[i]) for i in range(len(trace1))]

    # Calculate the average of the absolute differences
    average_difference = sum(differences) / len(differences)

    return differences

# Example usage:
#trace1 = [1, 2, 3, 4, 5]
#trace2 = [5, 4, 3, 2, 1]
#result = average_absolute_difference(trace1, trace2)
#print("Average Absolute Difference:", result)

#%%


#%%
#CURRENTLY USING TF_PARAMS_PATH
stage = 's2'

#0 for no, 1 for yes
if stage == 's2':
  init_only = 1
else:
  init_only = 0
#0 for turn decoding or visual stimuli decoding, 1 for correctness decoding
# IGNORE TRIAL TYPE
trial_type = 1


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

#TF_params_path = r'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params/1to1_FTP_TF_ne.csv'
TF_params_path = r'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params/1to1_INIT_TF_ne.csv'

#TF_params_path = r'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params/mice_FTP_TF_params_406065_231130.csv'
#TF_params_path = r'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params/mice_INIT_TF_params_406065_231128.csv'
#TF_params_path = r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params\short_test_mice_INIT_TF_params_406065_231128.csv"

#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/init_1n1e_permouse.csv' #5 init mice, 1n/e (no bias, >20 clean trials)
#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/mice_INIT_TF_params.csv' #all init (no bias, >20 clean trials)
#TF_params_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/mice_params/mice_FTP_TF_params.csv' #all ftp (no bias, >20 clean trials)
TF_params = pd.read_csv(TF_params_path)
# Initialize a list to store sessionID values
selected_sessions = []
selected_sessions_expert = []
selected_sessions_naive = []

#crossreg_path = r'/content/drive/MyDrive/Calcium_imaging_yarm_LurLab/mice_datasets/crossreg/m32l/crossreg_230731_230814/mappings.pkl'
#mappings = pd.read_pickle(crossreg_path)
#pairings = mappings.iloc[:, [1, 2]].values.tolist()
#pairings = [pair for pair in pairings if not (np.isnan(pair[0]) or np.isnan(pair[1]))]
#pairings = [[round(pair[0]), round(pair[1])] for pair in pairings]


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

behavior = []
folders_analyzed = []
temp.sort()
param_paths = []
cal_paths_raw = []
cal_paths_spike = []
rejected_cells_paths = []

cell_count_per_sess = []

for file in temp:
    param_paths.append(glob.glob(file + '/trial_parameters' + '*.csv'))
    #cal_paths_spike.append(glob.glob(file + '/spikerate' + '*.csv'))
    cal_paths_raw.append(glob.glob(file + '/calcium' + '*.csv'))
    rejected_cells_paths.append(glob.glob(file + '/rejected' + '*.csv'))

data_transposed_C_arrays_expert = []
data_transposed_W_arrays_expert = []
data_transposed_Lturn_arrays_expert = []
data_transposed_Rturn_arrays_expert = []
data_transposed_Lvis_arrays_expert = []
data_transposed_Rvis_arrays_expert = []

data_transposed_C_arrays_naive = []
data_transposed_W_arrays_naive = []
data_transposed_Lturn_arrays_naive = []
data_transposed_Rturn_arrays_naive = []
data_transposed_Lvis_arrays_naive = []
data_transposed_Rvis_arrays_naive = []

transposed_data_arrays = []
session_info_list = []
for i in range(len(param_paths)):
#for i in range(0,2):
#for i in range(8,11):
    try:
      current_folder = temp[i]
      cal = pd.read_csv(cal_paths_raw[i][0], engine='python')
      #cal_spike = pd.read_csv(cal_paths_spike[i][0])
      params = pd.read_csv(param_paths[i][0], engine='python')
      #reject = pd.read_csv(rejected_cells_paths[i][0])
      #reject = [item for sublist in reject.values.tolist() for item in sublist]


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

      ### save cell counts per session
      cell_count_per_sess.append(len(cal.columns))


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

      window_size = 6 #default 6
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


      parameters_Lturn = get_parameters_Lturn(cal,params)
      if current_folder in selected_sessions_expert:
        Lturn_cal_expert = get_trace_average_per_cell(cal,parameters_Lturn)
        Lturn_cal_expert = zscore(Lturn_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_Lturn_cal_expert = np.zeros_like(Lturn_cal_expert)
        # Apply the moving average to each column
        for i in range(Lturn_cal_expert.shape[1]):
            column = Lturn_cal_expert[:, i]
            for t in range(Lturn_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Lturn_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_Lturn_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        Lturn_cal_expert = smoothed_Lturn_cal_expert
      elif current_folder in selected_sessions_naive:
        Lturn_cal_naive = get_trace_average_per_cell(cal,parameters_Lturn)
        Lturn_cal_naive  = zscore(Lturn_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_Lturn_cal_naive = np.zeros_like(Lturn_cal_naive)
        # Apply the moving average to each column
        for i in range(Lturn_cal_naive.shape[1]):
            column = Lturn_cal_naive[:, i]
            for t in range(Lturn_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Lturn_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_Lturn_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        Lturn_cal_naive = smoothed_Lturn_cal_naive

      parameters_Rturn = get_parameters_Rturn(cal,params)
      if current_folder in selected_sessions_expert:
        Rturn_cal_expert = get_trace_average_per_cell(cal,parameters_Rturn)
        Rturn_cal_expert = zscore(Rturn_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_Rturn_cal_expert = np.zeros_like(Rturn_cal_expert)
        # Apply the moving average to each column
        for i in range(Rturn_cal_expert.shape[1]):
            column = Rturn_cal_expert[:, i]
            for t in range(Rturn_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Rturn_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_Rturn_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        Rturn_cal_expert = smoothed_Rturn_cal_expert
      elif current_folder in selected_sessions_naive:
        Rturn_cal_naive = get_trace_average_per_cell(cal,parameters_Rturn)
        Rturn_cal_naive  = zscore(Rturn_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_Rturn_cal_naive = np.zeros_like(Rturn_cal_naive)
        # Apply the moving average to each column
        for i in range(Rturn_cal_naive.shape[1]):
            column = Rturn_cal_naive[:, i]
            for t in range(Rturn_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Rturn_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_Rturn_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        Rturn_cal_naive = smoothed_Rturn_cal_naive

      parameters_Lvis = get_parameters_Lvis(cal,params)
      if current_folder in selected_sessions_expert:
        Lvis_cal_expert = get_trace_average_per_cell(cal,parameters_Lvis)
        Lvis_cal_expert = zscore(Lvis_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_Lvis_cal_expert = np.zeros_like(Lvis_cal_expert)
        # Apply the moving average to each column
        for i in range(Lvis_cal_expert.shape[1]):
            column = Lvis_cal_expert[:, i]
            for t in range(Lvis_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Lvis_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_Lvis_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        Lvis_cal_expert = smoothed_Lvis_cal_expert
      elif current_folder in selected_sessions_naive:
        Lvis_cal_naive = get_trace_average_per_cell(cal,parameters_Lvis)
        Lvis_cal_naive  = zscore(Lvis_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_Lvis_cal_naive = np.zeros_like(Lvis_cal_naive)
        # Apply the moving average to each column
        for i in range(Lvis_cal_naive.shape[1]):
            column = Lvis_cal_naive[:, i]
            for t in range(Lvis_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Lvis_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_Lvis_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        Lvis_cal_naive = smoothed_Lvis_cal_naive

      parameters_Rvis = get_parameters_Rvis(cal,params)
      if current_folder in selected_sessions_expert:
        Rvis_cal_expert = get_trace_average_per_cell(cal,parameters_Rvis)
        Rvis_cal_expert = zscore(Rvis_cal_expert, axis=0)
        # Create an array to store the result
        smoothed_Rvis_cal_expert = np.zeros_like(Rvis_cal_expert)
        # Apply the moving average to each column
        for i in range(Rvis_cal_expert.shape[1]):
            column = Rvis_cal_expert[:, i]
            for t in range(Rvis_cal_expert.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Rvis_cal_expert.shape[0], t + window_size // 2 + 1)
                smoothed_Rvis_cal_expert[t, i] = np.mean(column[start_idx:end_idx])
        Rvis_cal_expert = smoothed_Rvis_cal_expert
      elif current_folder in selected_sessions_naive:
        Rvis_cal_naive = get_trace_average_per_cell(cal,parameters_Rvis)
        Rvis_cal_naive  = zscore(Rvis_cal_naive , axis=0)
        # Create an array to store the result
        smoothed_Rvis_cal_naive = np.zeros_like(Rvis_cal_naive)
        # Apply the moving average to each column
        for i in range(Rvis_cal_naive.shape[1]):
            column = Rvis_cal_naive[:, i]
            for t in range(Rvis_cal_naive.shape[0]):
                start_idx = max(0, t - window_size // 2)
                end_idx = min(Rvis_cal_naive.shape[0], t + window_size // 2 + 1)
                smoothed_Rvis_cal_naive[t, i] = np.mean(column[start_idx:end_idx])
        Rvis_cal_naive = smoothed_Rvis_cal_naive

      #######################################################################

      # Transpose the data array to have cell IDs as rows and frames as columns
      data_transposed = data.T
      transposed_data_arrays.append(data_transposed)


      if current_folder in selected_sessions_expert:
          data_transposed_C_expert = C_cal_expert.T
          data_transposed_W_expert = W_cal_expert.T
          data_transposed_Lturn_expert = Lturn_cal_expert.T
          data_transposed_Rturn_expert = Rturn_cal_expert.T
          data_transposed_Lvis_expert = Lvis_cal_expert.T
          data_transposed_Rvis_expert = Rvis_cal_expert.T
          data_transposed_C_arrays_expert.append(data_transposed_C_expert)
          data_transposed_W_arrays_expert.append(data_transposed_W_expert)
          data_transposed_Lturn_arrays_expert.append(data_transposed_Lturn_expert)
          data_transposed_Rturn_arrays_expert.append(data_transposed_Rturn_expert)
          data_transposed_Lvis_arrays_expert.append(data_transposed_Lvis_expert)
          data_transposed_Rvis_arrays_expert.append(data_transposed_Rvis_expert)

          copy_shape = data_transposed_C_expert.shape
          filler = "filler"
          filler_array = np.full(copy_shape, filler)
          data_transposed_C_arrays_naive.append(filler_array)
          data_transposed_W_arrays_naive.append(filler_array)
          data_transposed_Lturn_arrays_naive.append(filler_array)
          data_transposed_Rturn_arrays_naive.append(filler_array)
          data_transposed_Lvis_arrays_naive.append(filler_array)
          data_transposed_Rvis_arrays_naive.append(filler_array)

      elif current_folder in selected_sessions_naive:
          data_transposed_C_naive = C_cal_naive.T
          data_transposed_W_naive = W_cal_naive.T
          data_transposed_Lturn_naive = Lturn_cal_naive.T
          data_transposed_Rturn_naive = Rturn_cal_naive.T
          data_transposed_Lvis_naive = Lvis_cal_naive.T
          data_transposed_Rvis_naive = Rvis_cal_naive.T
          data_transposed_C_arrays_naive.append(data_transposed_C_naive)
          data_transposed_W_arrays_naive.append(data_transposed_W_naive)
          data_transposed_Lturn_arrays_naive.append(data_transposed_Lturn_naive)
          data_transposed_Rturn_arrays_naive.append(data_transposed_Rturn_naive)
          data_transposed_Lvis_arrays_naive.append(data_transposed_Lvis_naive)
          data_transposed_Rvis_arrays_naive.append(data_transposed_Rvis_naive)

          copy_shape = data_transposed_C_naive.shape
          filler = "filler"
          filler_array = np.full(copy_shape, filler)
          data_transposed_C_arrays_expert.append(filler_array)
          data_transposed_W_arrays_expert.append(filler_array)
          data_transposed_Lturn_arrays_expert.append(filler_array)
          data_transposed_Rturn_arrays_expert.append(filler_array)
          data_transposed_Lvis_arrays_expert.append(filler_array)
          data_transposed_Rvis_arrays_expert.append(filler_array)

    except Exception as e:  # Catching all exceptions
        print(f"An error occurred: {e}")
        print('skipped: ', current_folder)
        pass


concatenated_data = np.vstack(transposed_data_arrays)

concatenated_data_C_expert = np.vstack(data_transposed_C_arrays_expert)
concatenated_data_W_expert = np.vstack(data_transposed_W_arrays_expert)
concatenated_data_Lturn_expert = np.vstack(data_transposed_Lturn_arrays_expert)
concatenated_data_Rturn_expert = np.vstack(data_transposed_Rturn_arrays_expert)
concatenated_data_Lvis_expert = np.vstack(data_transposed_Lvis_arrays_expert)
concatenated_data_Rvis_expert = np.vstack(data_transposed_Rvis_arrays_expert)

concatenated_data_C_naive = np.vstack(data_transposed_C_arrays_naive)
concatenated_data_W_naive = np.vstack(data_transposed_W_arrays_naive)
concatenated_data_Lturn_naive = np.vstack(data_transposed_Lturn_arrays_naive)
concatenated_data_Rturn_naive = np.vstack(data_transposed_Rturn_arrays_naive)
concatenated_data_Lvis_naive = np.vstack(data_transposed_Lvis_arrays_naive)
concatenated_data_Rvis_naive = np.vstack(data_transposed_Rvis_arrays_naive)

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
        print("removed the following cell index: ", index)

# Remove rows with NaNs
concatenated_data = np.delete(concatenated_data, nan_rows, axis=0)
concatenated_data_C_expert = np.delete(concatenated_data_C_expert, nan_rows, axis=0)
concatenated_data_W_expert = np.delete(concatenated_data_W_expert, nan_rows, axis=0)
concatenated_data_Lturn_expert = np.delete(concatenated_data_Lturn_expert, nan_rows, axis=0)
concatenated_data_Rturn_expert = np.delete(concatenated_data_Rturn_expert, nan_rows, axis=0)
concatenated_data_Lvis_expert = np.delete(concatenated_data_Lvis_expert, nan_rows, axis=0)
concatenated_data_Rvis_expert = np.delete(concatenated_data_Rvis_expert, nan_rows, axis=0)

concatenated_data_C_naive = np.delete(concatenated_data_C_naive, nan_rows, axis=0)
concatenated_data_W_naive = np.delete(concatenated_data_W_naive, nan_rows, axis=0)
concatenated_data_Lturn_naive = np.delete(concatenated_data_Lturn_naive, nan_rows, axis=0)
concatenated_data_Rturn_naive = np.delete(concatenated_data_Rturn_naive, nan_rows, axis=0)
concatenated_data_Lvis_naive = np.delete(concatenated_data_Lvis_naive, nan_rows, axis=0)
concatenated_data_Rvis_naive = np.delete(concatenated_data_Rvis_naive, nan_rows, axis=0)

expert_cellids = []
naive_cellids = []

for index, row in enumerate(concatenated_data_C_expert):
    if "filler" not in row:
        expert_cellids.append(index)
for index, row in enumerate(concatenated_data_C_naive):
    if "filler" not in row:
        naive_cellids.append(index)

#%%
# TESTING

print(cell_count_per_sess)
print(sum(cell_count_per_sess))

#%%
len(concatenated_data)

#%%
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# cluster using all data to get cellIDs
###### cluster
# Calculate the linkage matrix using 'ward' method
#Z = linkage(data_transposed, method='ward')
Z = linkage(concatenated_data, method='ward')

# Set the percentage threshold for determining the number of clusters
percentage_threshold = 0.2  # Adjust this value as needed, 0.2, lower = more clusters 0.05

# Calculate the height threshold based on the percentage
max_d = percentage_threshold * np.max(Z[:, 2])

# Calculate the height threshold for percentage_threshold = 0.05
#threshold_005 = 0.05 * np.max(Z[:, 2])

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
plt.figure()
dendrogram(Z, color_threshold=max_d, above_threshold_color='gray')

# Assign unique integer labels to leaf nodes
leaf_labels = dendrogram(Z, no_plot=True)['ivl']
label_map = {label: idx for idx, label in enumerate(leaf_labels)}
leaf_labels_int = [label_map[label] for label in leaf_labels]

## Apply colors to the leaf nodes (cell IDs) based on the cluster they belong to
leaf_colors = [color_dict[clusters[label]] for label in leaf_labels_int]
ax = plt.gca()
ax.set_xticklabels(ax.get_xticks(), rotation=90)  # Rotate x-axis labels for better readability
plt.bar(range(len(leaf_labels_int)), np.zeros(len(leaf_labels_int)), color=leaf_colors, align='center', width=0.5)
# Plot a horizontal line at the height threshold
plt.axhline(y=max_d, c='red', linestyle='dashed')
# Plot a green dashed line at the height threshold for percentage_threshold = 0.05
#plt.axhline(y=threshold_005, c='green', linestyle='dashed', label='Percentage Threshold = 0.05')



plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Cell IDs')
plt.ylabel('Distance')

try:
  #figure_path = r"C:\Users\aozgu\Downloads\hierarchical_clustering_dendrogram.eps"
  #figure_path = r"C:\Users\aozgu\Downloads\hierarchical_clustering_dendrogram.eps"
  figure_path = r"C:\Users\LurLab\Downloads\hierarchical_clustering_dendrogram" + ".eps"
  plt.savefig(figure_path, format='eps', bbox_inches='tight')
except:
  pass

# Dictionary to store cluster indices
cluster_indices_dict = {}
ncells_per_cluster = []

# Print the cell IDs for each cluster in list format and store indices
for cluster_id in range(1, num_clusters + 1):
    cluster_indices = np.where(clusters == cluster_id)[0]
    cell_ids = [int(leaf_labels_int[idx]) for idx in cluster_indices]
    cluster_name = f"cluster{cluster_id}"
    cluster_indices_dict[cluster_name] = cell_ids
    print(f"{cluster_name} = {cell_ids}")
    ncells_per_cluster.append(len(cell_ids))

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
cluster_expert_dict = {}
cluster_naive_dict = {}

# Iterate through each cluster and check if each value is in expert_cellids
for cluster, values in cluster_indices_dict.items():
    cluster_expert_list = []
    cluster_naive_list = []
    for value in values:
        if value in expert_cellids:
            cluster_expert_list.append(value)
        if value in naive_cellids:
            cluster_naive_list.append(value)
    #cluster_expert_dict[f'cluster{cluster[-1]}_expert'] = cluster_expert_list
    #cluster_naive_dict[f'cluster{cluster[-1]}_naive'] = cluster_naive_list
    cluster_expert_dict[f'{cluster}_expert'] = cluster_expert_list
    cluster_naive_dict[f'{cluster}_naive'] = cluster_naive_list

# Print the results
print("Expert Cellids:")
for key, value in cluster_expert_dict.items():
    print(f'{key}: {value}')

print("\nNaive Cellids:")
for key, value in cluster_naive_dict.items():
    print(f'{key}: {value}')

#%%
# without plotting (faster)
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# cluster using all data to get cellIDs
###### cluster
# Calculate the linkage matrix using 'ward' method
#Z = linkage(data_transposed, method='ward')
Z = linkage(concatenated_data, method='ward')

# Set the percentage threshold for determining the number of clusters
percentage_threshold = 0.2  # Adjust this value as needed, 0.2, lower = more clusters 0.05

# Calculate the height threshold based on the percentage
max_d = percentage_threshold * np.max(Z[:, 2])

# Calculate the height threshold for percentage_threshold = 0.05
#threshold_005 = 0.05 * np.max(Z[:, 2])

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
#plt.figure()
#dendrogram(Z, color_threshold=max_d, above_threshold_color='gray')

# Assign unique integer labels to leaf nodes
leaf_labels = dendrogram(Z, no_plot=True)['ivl']
label_map = {label: idx for idx, label in enumerate(leaf_labels)}
leaf_labels_int = [label_map[label] for label in leaf_labels]

## Apply colors to the leaf nodes (cell IDs) based on the cluster they belong to
leaf_colors = [color_dict[clusters[label]] for label in leaf_labels_int]
#ax = plt.gca()
#ax.set_xticklabels(ax.get_xticks(), rotation=90)  # Rotate x-axis labels for better readability
#plt.bar(range(len(leaf_labels_int)), np.zeros(len(leaf_labels_int)), color=leaf_colors, align='center', width=0.5)
## Plot a horizontal line at the height threshold
#plt.axhline(y=max_d, c='red', linestyle='dashed')
## Plot a green dashed line at the height threshold for percentage_threshold = 0.05
#plt.axhline(y=threshold_005, c='green', linestyle='dashed', label='Percentage Threshold = 0.05')



#plt.title('Hierarchical Clustering Dendrogram')
#plt.xlabel('Cell IDs')
#plt.ylabel('Distance')

#try:
#  figure_path = r"C:\Users\aozgu\Downloads\init65_02_thresh_hierarchical_clustering_dendrogram.eps"
#  plt.savefig(figure_path, format='eps', bbox_inches='tight')
#except:
#  pass

# Dictionary to store cluster indices
cluster_indices_dict = {}
ncells_per_cluster = []

# Print the cell IDs for each cluster in list format and store indices
for cluster_id in range(1, num_clusters + 1):
    cluster_indices = np.where(clusters == cluster_id)[0]
    cell_ids = [int(leaf_labels_int[idx]) for idx in cluster_indices]
    cluster_name = f"cluster{cluster_id}"
    cluster_indices_dict[cluster_name] = cell_ids
    print(f"{cluster_name} = {cell_ids}")
    ncells_per_cluster.append(len(cell_ids))

# Create dictionaries to store the cluster-specific expert and naive lists
cluster_expert_dict = {}
cluster_naive_dict = {}

# Iterate through each cluster and check if each value is in expert_cellids
for cluster, values in cluster_indices_dict.items():
    cluster_expert_list = []
    cluster_naive_list = []
    for value in values:
        if value in expert_cellids:
            cluster_expert_list.append(value)
        if value in naive_cellids:
            cluster_naive_list.append(value)
    cluster_expert_dict[f'{cluster}_expert'] = cluster_expert_list
    cluster_naive_dict[f'{cluster}_naive'] = cluster_naive_list

# Print the results
print("Expert Cellids:")
for key, value in cluster_expert_dict.items():
    print(f'{key}: {value}')

print("\nNaive Cellids:")
for key, value in cluster_naive_dict.items():
    print(f'{key}: {value}')

#%%
# Initialize a dictionary to store the count of clusters for each session
session_cluster_count = {}

# Initialize the start index
start_index = 0

# Iterate over each session's cell count
for session_count in cell_count_per_sess:
    # Initialize dictionary for current session
    session_clusters = {}

    # Calculate end index
    end_index = start_index + session_count

    # Iterate over each cell ID within the current session's range
    for cell_id in range(start_index, end_index):
        # Find the cluster for the current cell ID
        for cluster_id in range(1, num_clusters + 1):
            cluster_indices = np.where(clusters == cluster_id)[0]
            if cell_id in cluster_indices:
                # Increment the count of the cluster for the current session
                session_clusters[cluster_id] = session_clusters.get(cluster_id, 0) + 1

    # Add the cluster counts for the current session to the session_cluster_count dictionary
    session_cluster_count[start_index] = session_clusters

    # Update start_index for the next session
    start_index = end_index

# Print the cluster counts for each session
for session_index, cluster_counts in session_cluster_count.items():
    print(f"Session {session_index}:")
    # Sort the clusters by cluster ID
    sorted_cluster_counts = dict(sorted(cluster_counts.items()))
    for cluster_id, count in sorted_cluster_counts.items():
        print(f"Cluster {cluster_id}: {count} cells")


#%%
# get percentages
# Initialize a dictionary to store the count of clusters for each session
session_cluster_count = {}

# Initialize the start index
start_index = 0

# Initialize a dictionary to store the percentage of clusters for each session
session_cluster_percentage = {}

# Iterate over each session's cell count
for session_index, session_count in enumerate(cell_count_per_sess):
    # Initialize dictionary for current session
    session_clusters = {}

    # Calculate end index
    end_index = start_index + session_count

    # Count total cells in current session
    total_cells_session = end_index - start_index

    # Iterate over each cluster ID within the current session
    for cluster_id in range(1, num_clusters + 1):
        # Initialize count for current cluster
        cluster_count = 0

        # Iterate over each cell ID within the current session's range
        for cell_id in range(start_index, end_index):
            # Find the cluster for the current cell ID
            cluster_indices = np.where(clusters == cluster_id)[0]
            if cell_id in cluster_indices:
                # Increment the count of the cluster for the current session
                cluster_count += 1

        # Calculate percentage for current cluster
        cluster_percentage = (cluster_count / total_cells_session) * 100 if total_cells_session != 0 else 0

        # Add the cluster count for the current session to the session_clusters dictionary
        session_clusters[cluster_id] = cluster_percentage

    # Add the cluster percentages for the current session to the session_cluster_percentage dictionary
    session_cluster_percentage[start_index] = session_clusters

    # Update start_index for the next session
    start_index = end_index

# Print the cluster percentages for each session
for session_index, cluster_percentages in session_cluster_percentage.items():
    print(f"Session {session_index}:")
    # Sort the clusters by cluster ID
    sorted_cluster_percentages = dict(sorted(cluster_percentages.items()))
    for cluster_id, percentage in sorted_cluster_percentages.items():
        print(f"Cluster {cluster_id}: {percentage:.2f}% of cells")


#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Assuming you already have your data and configurations up to this point

# Extract every other session starting from 0 and calculate the sum of percentages for cluster 1
regions = []
sum_lengths = []
should_include = True
for session_index, cluster_percentages in session_cluster_percentage.items():
    if should_include:
        regions.append(session_index)
        sum_length = cluster_percentages.get(6, 0)  # Get the percentage for cluster 1
        sum_lengths.append(sum_length)
    should_include = not should_include  # Toggle the flag

# Create the DataFrame
summary_all = pd.DataFrame({'region': regions, 'sum_length': sum_lengths})

#####

# Some layout stuff ----------------------------------------------

# Initialize layout in polar coordinates
fig, ax = plt.subplots(figsize=(9, 12.6), subplot_kw={"projection": "polar"})

# Set background color to white, both axis and figure.
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Calculate the positions of bars
num_categories = len(summary_all)
category_spacing = 2 * np.pi / num_categories

# Calculate the angles for evenly spaced bars
angles_adjusted = np.linspace(0, 2 * np.pi, len(summary_all), endpoint=False)

# Define bar width
bar_width = category_spacing * 0.8  # Adjust the factor to control the width of bars

ax.set_ylim(-10, 50)

# Add geometries to the plot -------------------------------------

#COLORS = ['#6C5B7B','#C06C84','#F67280','#F8B195']
COLORS = ['#6C5B7B', '#C06C84', '#F67280', '#F8B195', '#355C7D']



# Add bars to represent the cumulative track lengths
ax.bar(angles_adjusted, summary_all['sum_length'], color=COLORS, alpha=0.9, width=bar_width, zorder=10)

# Add labels for the regions -------------------------------------

# Set the labels
ax.set_xticks(angles_adjusted)
ax.set_xticklabels(summary_all['region'], size=12)

# Remove unnecessary guides ---------------------------------------

# Remove lines for polar axis (x)
ax.xaxis.grid(False)

# Remove spines
ax.spines["start"].set_color("none")
ax.spines["polar"].set_color("none")

# Adjust padding of the x axis labels ----------------------------
# This is going to add extra space around the labels for the
# ticks of the x axis.
XTICKS = ax.xaxis.get_major_ticks()
for tick in XTICKS:
    tick.set_pad(10)

# Remove tick labels
ax.set_xticklabels([])

# Remove ticks
ax.set_xticks([])

try:
  ### TO SAVE OUTPUTS
  ## Generate the file name dynamically based on the cluster index
  #figure_path = r"C:\Users\aozgu\Downloads\radial_0" + ".eps"
  figure_path = r"C:\Users\LurLab\Downloads\radial_0" + ".eps"
  #figure_path = r"C:\Users\LurLab\Downloads\cluster" + f'{i+1}' + ".svg"
  ## Save the figure
  plt.savefig(figure_path, format='eps', bbox_inches='tight')
except:
  pass


plt.show()


#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Assuming you already have your data and configurations up to this point

# Extract every other session starting from 1 and calculate the sum of percentages for cluster 1
regions = []
sum_lengths = []
should_include = False
for session_index, cluster_percentages in session_cluster_percentage.items():
    if should_include:
        regions.append(session_index)
        sum_length = cluster_percentages.get(6, 0)  # Get the percentage for cluster 1
        sum_lengths.append(sum_length)
    should_include = not should_include  # Toggle the flag

# Create the DataFrame
summary_all = pd.DataFrame({'region': regions, 'sum_length': sum_lengths})

#####

# Some layout stuff ----------------------------------------------

# Initialize layout in polar coordinates
fig, ax = plt.subplots(figsize=(9, 12.6), subplot_kw={"projection": "polar"})

# Set background color to white, both axis and figure.
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Calculate the positions of bars
num_categories = len(summary_all)
category_spacing = 2 * np.pi / num_categories

# Calculate the angles for evenly spaced bars
angles_adjusted = np.linspace(0, 2 * np.pi, len(summary_all), endpoint=False)

# Define bar width
bar_width = category_spacing * 0.8  # Adjust the factor to control the width of bars

ax.set_ylim(-10, 50)

# Add geometries to the plot -------------------------------------

COLORS = ['#6C5B7B', '#C06C84', '#F67280', '#F8B195', '#355C7D']
# Add bars to represent the cumulative track lengths
ax.bar(angles_adjusted, summary_all['sum_length'], color=COLORS, alpha=0.9, width=bar_width, zorder=10)

# Add labels for the regions -------------------------------------

# Set the labels
ax.set_xticks(angles_adjusted)
ax.set_xticklabels(summary_all['region'], size=12)

# Remove unnecessary guides ---------------------------------------

# Remove lines for polar axis (x)
ax.xaxis.grid(False)

# Remove spines
ax.spines["start"].set_color("none")
ax.spines["polar"].set_color("none")

# Adjust padding of the x axis labels ----------------------------
# This is going to add extra space around the labels for the
# ticks of the x axis.
XTICKS = ax.xaxis.get_major_ticks()
for tick in XTICKS:
    tick.set_pad(10)

# Remove tick labels
ax.set_xticklabels([])

# Remove ticks
ax.set_xticks([])

try:
  ### TO SAVE OUTPUTS
  ## Generate the file name dynamically based on the cluster index
  #figure_path = r"C:\Users\aozgu\Downloads\radial_1" + ".eps"
  figure_path = r"C:\Users\LurLab\Downloads\radial_1" + ".eps"
  #figure_path = r"C:\Users\LurLab\Downloads\cluster" + f'{i+1}' + ".svg"
  ## Save the figure
  plt.savefig(figure_path, format='eps', bbox_inches='tight')
except:
  pass

plt.show()


#%%
summary_all

#%%
# for cell counts per cluster per session

import numpy as np
import matplotlib.pyplot as plt

# Initialize a dictionary to store the count of clusters for each session
session_cluster_count = {}

# Initialize the start index
start_index = 0

# Iterate over each session's cell count
for session_index, session_count in enumerate(cell_count_per_sess):
    # Initialize dictionary for current session
    session_clusters = {}

    # Calculate end index
    end_index = start_index + session_count

    # Iterate over each cell ID within the current session's range
    for cell_id in range(start_index, end_index):
        # Find the cluster for the current cell ID
        for cluster_id in range(1, num_clusters + 1):
            cluster_indices = np.where(clusters == cluster_id)[0]
            if cell_id in cluster_indices:
                # Increment the count of the cluster for the current session
                session_clusters[cluster_id] = session_clusters.get(cluster_id, 0) + 1

    # Add the cluster counts for the current session to the session_cluster_count dictionary
    session_cluster_count[session_index] = session_clusters

    # Update start_index for the next session
    start_index = end_index

# Plot histogram for total cell count for each session
for session_index, cluster_counts in session_cluster_count.items():
    plt.figure()
    plt.bar(range(1, num_clusters + 1), [cluster_counts.get(cluster_id, 0) for cluster_id in range(1, num_clusters + 1)], color='blue', edgecolor='black')
    plt.title(f'Total Cell Count for Each Cluster - Session {session_index + 1}')
    plt.xlabel('Cluster')
    plt.ylabel('Cell Count')
    plt.xticks(range(1, num_clusters + 1))
    #plt.grid(True)
    try:
      ### TO SAVE OUTPUTS
      ## Generate the file name dynamically based on the session index
      figure_path = r"C:\Users\aozgu\Downloads\session" + f'{session_index + 1}' + ".eps"
      ## Save the figure
      plt.savefig(figure_path, format='eps', bbox_inches='tight')
    except:
      pass

    plt.show()

# Plot histogram for each cluster across all sessions
cluster_total_count = {cluster_id: sum(cluster_counts.get(cluster_id, 0) for cluster_counts in session_cluster_count.values()) for cluster_id in range(1, num_clusters + 1)}
plt.figure()
plt.bar(range(1, num_clusters + 1), cluster_total_count.values(), color='blue', edgecolor='black')
plt.title('Total Cell Count for Each Cluster Across Sessions')
plt.xlabel('Cluster')
plt.ylabel('Cell Count')
plt.xticks(range(1, num_clusters + 1))
#plt.grid(True)

try:
  ### TO SAVE OUTPUTS
  ## Generate the file name dynamically based on the session index
  figure_path = r"C:\Users\aozgu\Downloads\ALLsessions" + ".eps"
  ## Save the figure
  plt.savefig(figure_path, format='eps', bbox_inches='tight')
except:
  pass

plt.show()



#%%
import numpy as np
import matplotlib.pyplot as plt

# Determine the number of rows needed for subplot pairs
num_pairs = len(cell_count_per_sess) // 2
if len(cell_count_per_sess) % 2 != 0:
    num_pairs += 1

# Create subplots
fig, axes = plt.subplots(num_pairs, 2, figsize=(12, 6*num_pairs))

# Initialize a dictionary to store the count of clusters for each session
session_cluster_count = {}

# Initialize the start index
start_index = 0

# Iterate over each session's cell count
for session_index, session_count in enumerate(cell_count_per_sess):
    # Initialize dictionary for current session
    session_clusters = {}

    # Calculate end index
    end_index = start_index + session_count

    # Iterate over each cell ID within the current session's range
    for cell_id in range(start_index, end_index):
        # Find the cluster for the current cell ID
        for cluster_id in range(1, num_clusters + 1):
            cluster_indices = np.where(clusters == cluster_id)[0]
            if cell_id in cluster_indices:
                # Increment the count of the cluster for the current session
                session_clusters[cluster_id] = session_clusters.get(cluster_id, 0) + 1

    # Add the cluster counts for the current session to the session_cluster_count dictionary
    session_cluster_count[session_index] = session_clusters

    # Update start_index for the next session
    start_index = end_index

    # Plot histogram for total cell count for each session
    ax = axes[session_index // 2, session_index % 2]
    ax.bar(range(1, num_clusters + 1), [session_clusters.get(cluster_id, 0) for cluster_id in range(1, num_clusters + 1)], color='blue', edgecolor='black')
    ax.set_title(f'Total Cell Count for Each Cluster - Session {session_index + 1}')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Cell Count')
    ax.set_xticks(range(1, num_clusters + 1))

# Hide empty subplots
for i in range(len(cell_count_per_sess), num_pairs * 2):
    axes[i // 2, i % 2].axis('off')

# Adjust layout
plt.tight_layout()

# Save the figure
figure_path = r"C:\Users\aozgu\Downloads\session_pairs.png"
plt.savefig(figure_path, format='png', bbox_inches='tight')

plt.show()


#%%
import numpy as np
import matplotlib.pyplot as plt

# Initialize a dictionary to store the count of clusters for each session
session_cluster_count = {}

# Initialize the start index
start_index = 0

# Iterate over each session's cell count
for session_index, session_count in enumerate(cell_count_per_sess):
    # Initialize dictionary for current session
    session_clusters = {}

    # Calculate end index
    end_index = start_index + session_count

    # Iterate over each cell ID within the current session's range
    for cell_id in range(start_index, end_index):
        # Find the cluster for the current cell ID
        for cluster_id in range(1, num_clusters + 1):
            cluster_indices = np.where(clusters == cluster_id)[0]
            if cell_id in cluster_indices:
                # Increment the count of the cluster for the current session
                session_clusters[cluster_id] = session_clusters.get(cluster_id, 0) + 1

    # Add the cluster counts for the current session to the session_cluster_count dictionary
    session_cluster_count[session_index] = session_clusters

    # Update start_index for the next session
    start_index = end_index

# Calculate the total number of cells for each cluster across all sessions
cluster_total_count = {cluster_id: sum(cluster_counts.get(cluster_id, 0) for cluster_counts in session_cluster_count.values()) for cluster_id in range(1, num_clusters + 1)}

# Calculate the percentage of cells in each cluster for each session
session_cluster_percentage = {}
for session_index, cluster_counts in session_cluster_count.items():
    session_total_cells = sum(cluster_counts.values())
    session_cluster_percentage[session_index] = {cluster_id: (count / session_total_cells) * 100 if session_total_cells != 0 else 0 for cluster_id, count in cluster_counts.items()}

# Calculate the difference in percentage of cells in each cluster between paired sessions
paired_session_diff = {}
for i in range(0, len(session_cluster_percentage), 2):
    session1_index = i
    session2_index = i + 1
    paired_diff = {}
    for cluster_id in range(1, num_clusters + 1):
        percent1 = session_cluster_percentage.get(session1_index, {}).get(cluster_id, 0)
        percent2 = session_cluster_percentage.get(session2_index, {}).get(cluster_id, 0)
        paired_diff[cluster_id] = percent2 - percent1  # Switching the subtraction order
    paired_session_diff[(session1_index, session2_index)] = paired_diff

# Print out the values for every bar
for i, (session1_index, session2_index) in enumerate(sorted(paired_session_diff.keys())):
    paired_diff_data = paired_session_diff[(session1_index, session2_index)]
    print(f'Sessions {session2_index + 1} - {session1_index + 1}:')
    for cluster_id, diff_value in paired_diff_data.items():
        #print(f'Cluster {cluster_id}: {diff_value:.2f}%')
        print(f'{diff_value:.2f}')
    print()


#%%
# for odd (naive) sessions and even (expert) sessions

import numpy as np
import matplotlib.pyplot as plt

# Initialize dictionaries to store the count of clusters for odd and even sessions
odd_session_cluster_count = {}
even_session_cluster_count = {}

# Iterate over each session's cell count
for session_index, session_count in enumerate(cell_count_per_sess):
    # Initialize dictionary for current session
    session_clusters = {}

    # Calculate end index
    end_index = sum(cell_count_per_sess[:session_index+1])
    start_index = end_index - session_count

    # Check if session is odd or even
    if session_index % 2 == 0:
        session_dict = odd_session_cluster_count
    else:
        session_dict = even_session_cluster_count

    # Iterate over each cell ID within the current session's range
    for cell_id in range(start_index, end_index):
        # Find the cluster for the current cell ID
        for cluster_id in range(1, num_clusters + 1):
            cluster_indices = np.where(clusters == cluster_id)[0]
            if cell_id in cluster_indices:
                # Increment the count of the cluster for the current session
                session_dict[cluster_id] = session_dict.get(cluster_id, 0) + 1

# Plot histogram for total cell count for odd sessions
plt.figure()
plt.bar(range(1, num_clusters + 1), [odd_session_cluster_count.get(cluster_id, 0) for cluster_id in range(1, num_clusters + 1)], color='blue', edgecolor='black')
plt.title('Total Cell Count for Each Cluster - Odd Sessions')
plt.xlabel('Cluster')
plt.ylabel('Cell Count')
plt.xticks(range(1, num_clusters + 1))
plt.grid(True)
plt.show()

# Plot histogram for total cell count for even sessions
plt.figure()
plt.bar(range(1, num_clusters + 1), [even_session_cluster_count.get(cluster_id, 0) for cluster_id in range(1, num_clusters + 1)], color='blue', edgecolor='black')
plt.title('Total Cell Count for Each Cluster - Even Sessions')
plt.xlabel('Cluster')
plt.ylabel('Cell Count')
plt.xticks(range(1, num_clusters + 1))
plt.grid(True)
plt.show()


#%%
# Iterate over clusters and print the length of cell IDs for each cluster
for cluster_name, cell_ids in cluster_indices_dict.items():
    cluster_length = len(cell_ids)
    print(f"Length of cell IDs in {cluster_name}: {cluster_length}")


#%%
# Create empty arrays for the data
C_hm_expert = []
W_hm_expert = []
Lturn_hm_expert = []
Rturn_hm_expert = []
Lvis_hm_expert = []
Rvis_hm_expert = []

# Iterate through each cluster
for cluster, row_numbers in cluster_expert_dict.items():
    for row_number in row_numbers:
        if 0 <= row_number < len(concatenated_data_C_expert):
            c_data = concatenated_data_C_expert[row_number]
            C_hm_expert.append(c_data)

        if 0 <= row_number < len(concatenated_data_W_expert):
            w_data = concatenated_data_W_expert[row_number]
            W_hm_expert.append(w_data)

        if 0 <= row_number < len(concatenated_data_Lturn_expert):
            lturn_data = concatenated_data_Lturn_expert[row_number]
            Lturn_hm_expert.append(lturn_data)

        if 0 <= row_number < len(concatenated_data_Rturn_expert):
            rturn_data = concatenated_data_Rturn_expert[row_number]
            Rturn_hm_expert.append(rturn_data)

        if 0 <= row_number < len(concatenated_data_Lvis_expert):
            lvis_data = concatenated_data_Lvis_expert[row_number]
            Lvis_hm_expert.append(lvis_data)

        if 0 <= row_number < len(concatenated_data_Rvis_expert):
            rvis_data = concatenated_data_Rvis_expert[row_number]
            Rvis_hm_expert.append(rvis_data)


# Convert the lists to numpy arrays if needed
C_hm_expert = np.array(C_hm_expert)
W_hm_expert = np.array(W_hm_expert)
Lturn_hm_expert = np.array(Lturn_hm_expert)
Rturn_hm_expert = np.array(Rturn_hm_expert)
Lvis_hm_expert = np.array(Lvis_hm_expert)
Rvis_hm_expert = np.array(Rvis_hm_expert)

# Create empty arrays for the data
C_hm_naive = []
W_hm_naive = []
Lturn_hm_naive = []
Rturn_hm_naive = []
Lvis_hm_naive = []
Rvis_hm_naive = []

# Iterate through each cluster
for cluster, row_numbers in cluster_naive_dict.items():
    for row_number in row_numbers:
        if 0 <= row_number < len(concatenated_data_C_naive):
            c_data = concatenated_data_C_naive[row_number]
            C_hm_naive.append(c_data)

        if 0 <= row_number < len(concatenated_data_W_naive):
            w_data = concatenated_data_W_naive[row_number]
            W_hm_naive.append(w_data)

        if 0 <= row_number < len(concatenated_data_Lturn_naive):
            lturn_data = concatenated_data_Lturn_naive[row_number]
            Lturn_hm_naive.append(lturn_data)

        if 0 <= row_number < len(concatenated_data_Rturn_naive):
            rturn_data = concatenated_data_Rturn_naive[row_number]
            Rturn_hm_naive.append(rturn_data)

        if 0 <= row_number < len(concatenated_data_Lvis_naive):
            lvis_data = concatenated_data_Lvis_naive[row_number]
            Lvis_hm_naive.append(lvis_data)

        if 0 <= row_number < len(concatenated_data_Rvis_naive):
            rvis_data = concatenated_data_Rvis_naive[row_number]
            Rvis_hm_naive.append(rvis_data)

# Convert the lists to numpy arrays if needed
C_hm_naive = np.array(C_hm_naive)
W_hm_naive = np.array(W_hm_naive)
Lturn_hm_naive = np.array(Lturn_hm_naive)
Rturn_hm_naive = np.array(Rturn_hm_naive)
Lvis_hm_naive = np.array(Lvis_hm_naive)
Rvis_hm_naive = np.array(Rvis_hm_naive)

#%%
# ### IF DONT WANT NORMALIZATION...

# C_hm_expert = C_hm_expert.astype(float)
# W_hm_expert = W_hm_expert.astype(float)
# Lturn_hm_expert = Lturn_hm_expert.astype(float)
# Rturn_hm_expert = Rturn_hm_expert.astype(float)
# Lvis_hm_expert = Lvis_hm_expert.astype(float)
# Rvis_hm_expert = Rvis_hm_expert.astype(float)

# C_hm_naive = C_hm_naive.astype(float)
# W_hm_naive = W_hm_naive.astype(float)
# Lturn_hm_naive = Lturn_hm_naive.astype(float)
# Rturn_hm_naive = Rturn_hm_naive.astype(float)
# Lvis_hm_naive = Lvis_hm_naive.astype(float)
# Rvis_hm_naive = Rvis_hm_naive.astype(float)


# C_hm_expert_normalized = C_hm_expert
# W_hm_expert_normalized = W_hm_expert
# Lturn_hm_expert_normalized = Lturn_hm_expert
# Rturn_hm_expert_normalized = Rturn_hm_expert
# Lvis_hm_expert_normalized = Lvis_hm_expert
# Rvis_hm_expert_normalized = Rvis_hm_expert

# C_hm_naive_normalized = C_hm_naive
# W_hm_naive_normalized = W_hm_naive
# Lturn_hm_naive_normalized = Lturn_hm_naive
# Rturn_hm_naive_normalized = Rturn_hm_naive
# Lvis_hm_naive_normalized = Lvis_hm_naive
# Rvis_hm_naive_normalized = Rvis_hm_naive

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(C_hm_expert_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# C_hm_expert_normalized = np.nan_to_num(C_hm_expert_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(W_hm_expert_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# W_hm_expert_normalized = np.nan_to_num(W_hm_expert_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Lturn_hm_expert_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Lturn_hm_expert_normalized = np.nan_to_num(Lturn_hm_expert_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Rturn_hm_expert_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Rturn_hm_expert_normalized = np.nan_to_num(Rturn_hm_expert_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Lvis_hm_expert_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Lvis_hm_expert_normalized = np.nan_to_num(Lvis_hm_expert_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Rvis_hm_expert_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Rvis_hm_expert_normalized = np.nan_to_num(Rvis_hm_expert_normalized, nan=0.0)

# ###########
# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(C_hm_naive_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# C_hm_naive_normalized = np.nan_to_num(C_hm_naive_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(W_hm_naive_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# W_hm_naive_normalized = np.nan_to_num(W_hm_naive_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Lturn_hm_naive_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Lturn_hm_naive_normalized = np.nan_to_num(Lturn_hm_naive_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Rturn_hm_naive_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Rturn_hm_naive_normalized = np.nan_to_num(Rturn_hm_naive_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Lvis_hm_naive_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Lvis_hm_naive_normalized = np.nan_to_num(Lvis_hm_naive_normalized, nan=0.0)

# # Print the row numbers with NaN values (now replaced with 0s)
# nan_indices = np.isnan(Rvis_hm_naive_normalized).any(axis=1)
# for row_index in np.where(nan_indices)[0]:
#     print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
# Rvis_hm_naive_normalized = np.nan_to_num(Rvis_hm_naive_normalized, nan=0.0)


#%%
### IF WANT NORMALIZATION...
# Normalize the 4 arrays

# Convert the input arrays to numeric data type (assuming they are originally strings)
C_hm_expert = C_hm_expert.astype(float)
W_hm_expert = W_hm_expert.astype(float)
Lturn_hm_expert = Lturn_hm_expert.astype(float)
Rturn_hm_expert = Rturn_hm_expert.astype(float)
Lvis_hm_expert = Lvis_hm_expert.astype(float)
Rvis_hm_expert = Rvis_hm_expert.astype(float)

# Create empty arrays for the normalized data
C_hm_expert_normalized = []
W_hm_expert_normalized = []
Lturn_hm_expert_normalized = []
Rturn_hm_expert_normalized = []
Lvis_hm_expert_normalized = []
Rvis_hm_expert_normalized = []

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

for row in Lturn_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Lturn_hm_expert_normalized.append(normalized_row)

for row in Rturn_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Rturn_hm_expert_normalized.append(normalized_row)

for row in Lvis_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Lvis_hm_expert_normalized.append(normalized_row)

for row in Rvis_hm_expert:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Rvis_hm_expert_normalized.append(normalized_row)

# Convert the lists to numpy arrays
C_hm_expert_normalized = np.array(C_hm_expert_normalized)
W_hm_expert_normalized = np.array(W_hm_expert_normalized)
Lturn_hm_expert_normalized = np.array(Lturn_hm_expert_normalized)
Rturn_hm_expert_normalized = np.array(Rturn_hm_expert_normalized)
Lvis_hm_expert_normalized = np.array(Lvis_hm_expert_normalized)
Rvis_hm_expert_normalized = np.array(Rvis_hm_expert_normalized)

# Convert the input arrays to numeric data type (assuming they are originally strings)
C_hm_naive = C_hm_naive.astype(float)
W_hm_naive = W_hm_naive.astype(float)
Lturn_hm_naive = Lturn_hm_naive.astype(float)
Rturn_hm_naive = Rturn_hm_naive.astype(float)
Lvis_hm_naive = Lvis_hm_naive.astype(float)
Rvis_hm_naive = Rvis_hm_naive.astype(float)

# Create empty arrays for the normalized data
C_hm_naive_normalized = []
W_hm_naive_normalized = []
Lturn_hm_naive_normalized = []
Rturn_hm_naive_normalized = []
Lvis_hm_naive_normalized = []
Rvis_hm_naive_normalized = []

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

for row in Lturn_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Lturn_hm_naive_normalized.append(normalized_row)

for row in Rturn_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Rturn_hm_naive_normalized.append(normalized_row)

for row in Lvis_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Lvis_hm_naive_normalized.append(normalized_row)

for row in Rvis_hm_naive:
    min_val = row.min()
    max_val = row.max()
    if max_val != min_val:
        normalized_row = (row - min_val) / (max_val - min_val)
    else:
        normalized_row = row
    Rvis_hm_naive_normalized.append(normalized_row)

# Convert the lists to numpy arrays
C_hm_naive_normalized = np.array(C_hm_naive_normalized)
W_hm_naive_normalized = np.array(W_hm_naive_normalized)
Lturn_hm_naive_normalized = np.array(Lturn_hm_naive_normalized)
Rturn_hm_naive_normalized = np.array(Rturn_hm_naive_normalized)
Lvis_hm_naive_normalized = np.array(Lvis_hm_naive_normalized)
Rvis_hm_naive_normalized = np.array(Rvis_hm_naive_normalized)

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
nan_indices = np.isnan(Lturn_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Lturn_hm_expert_normalized = np.nan_to_num(Lturn_hm_expert_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(Rturn_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Rturn_hm_expert_normalized = np.nan_to_num(Rturn_hm_expert_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(Lvis_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Lvis_hm_expert_normalized = np.nan_to_num(Lvis_hm_expert_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(Rvis_hm_expert_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Rvis_hm_expert_normalized = np.nan_to_num(Rvis_hm_expert_normalized, nan=0.0)

######

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
nan_indices = np.isnan(Lturn_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Lturn_hm_naive_normalized = np.nan_to_num(Lturn_hm_naive_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(Rturn_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Rturn_hm_naive_normalized = np.nan_to_num(Rturn_hm_naive_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(Lvis_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Lvis_hm_naive_normalized = np.nan_to_num(Lvis_hm_naive_normalized, nan=0.0)

# Print the row numbers with NaN values (now replaced with 0s)
nan_indices = np.isnan(Rvis_hm_naive_normalized).any(axis=1)
for row_index in np.where(nan_indices)[0]:
    print(f"Row {row_index} originally contained NaN values and has been replaced with 0s.")
Rvis_hm_naive_normalized = np.nan_to_num(Rvis_hm_naive_normalized, nan=0.0)

#%%


#%%
# new heatmap test:

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Min-Max scaling function
def min_max_scaling(row):
    min_val = np.min(row)
    max_val = np.max(row)
    scaled_row = (row - min_val) / (max_val - min_val)
    return scaled_row

# Apply Min-Max scaling to each row
normalized_data = np.apply_along_axis(min_max_scaling, axis=1, arr=concatenated_data)

# Assume concatenated_data is your heatmap data
heatmap_data = normalized_data

# Create a figure with a single subplot
fig, axs = plt.subplots(1, 1, figsize=(8, 4))

# # Customize Seaborn colormap
# cmap = sns.color_palette("RdYlBu_r", as_cmap=True)

# # Customize Seaborn colormap (use "Blues" colormap for white-to-blue scale)
# cmap = sns.color_palette("Blues", as_cmap=True)

# Customize Seaborn colormap (use "RdBu_r" colormap for red-to-blue scale)
cmap = sns.color_palette("RdBu_r", as_cmap=True)

# # Plot the heatmap with Seaborn
# sns.heatmap(heatmap_data, cmap=cmap, ax=axs, cbar=True, cbar_kws={'label': 'Color Scale'})

# Plot the heatmap with Seaborn (without color scale)
sns.heatmap(heatmap_data, cmap=cmap, ax=axs, cbar=False)

# Set title for the heatmap
#axs.set_title('Concatenated Heatmap')

# Customize cluster lines
rows_drawn = 0
for cluster_length in cluster_indices_dict.values():
    cluster_length = len(cluster_length)
    if rows_drawn + cluster_length < heatmap_data.shape[0]:
        axs.axhline(rows_drawn + cluster_length - 0.5, color='white', linestyle='dashed', linewidth=1)
    rows_drawn += cluster_length

# # Manually set evenly spaced x-axis ticks
# x_ticks = np.linspace(0, heatmap_data.shape[1], num=heatmap_data.shape[1], endpoint=False)
# axs.set_xticks(x_ticks)

# # Manually set x-axis tick labels
# x_tick_labels = [str(int(x)) for x in x_ticks]
# axs.set_xticklabels(x_tick_labels)

# Hide y-axis ticks
axs.set_yticks([])

# Hide x-axis ticks
axs.set_xticks([])

# Add a black vertical line at x = 16
axs.axvline(x=16, color='black', linestyle='-', linewidth=1)

# Adjust spacing between subplots
plt.tight_layout()

try:
  ### TO SAVE OUTPUTS
  ## Generate the file name dynamically based on the cluster index
  figure_path = r"C:\Users\aozgu\Downloads\overall_clustering_heatmap" + ".png"
  #figure_path = r"C:\Users\LurLab\Downloads\init65_0.2thresh_overall_clustering_heatmap" + ".png"
  ## Save the figure
  plt.savefig(figure_path, format='png', bbox_inches='tight')
except:
  pass

# Show the plot
plt.show()


#%%
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Min-Max scaling function
def min_max_scaling(row):
    min_val = np.min(row)
    max_val = np.max(row)
    scaled_row = (row - min_val) / (max_val - min_val)
    return scaled_row

# Apply Min-Max scaling to each row
normalized_data = np.apply_along_axis(min_max_scaling, axis=1, arr=concatenated_data)

# Assume concatenated_data is your heatmap data
heatmap_data = normalized_data

# Create a figure with a single subplot
fig, axs = plt.subplots(1, 1, figsize=(8, 4))

# Customize Seaborn colormap
cmap = sns.color_palette("RdYlBu_r", as_cmap=True)

# Plot the heatmap with Seaborn
sns.heatmap(heatmap_data, cmap=cmap, ax=axs, cbar=True, cbar_kws={'label': 'Color Scale'})

# Set title for the heatmap
axs.set_title('Concatenated Heatmap')

# Customize cluster lines
rows_drawn = 0
for cluster_length in cluster_indices_dict.values():
    cluster_length = len(cluster_length)
    if rows_drawn + cluster_length < heatmap_data.shape[0]:
        axs.axhline(rows_drawn + cluster_length - 0.5, color='white', linestyle='dashed', linewidth=1)
    rows_drawn += cluster_length

# Manually set evenly spaced x-axis ticks
x_ticks = np.linspace(0, heatmap_data.shape[1], num=heatmap_data.shape[1], endpoint=False)
axs.set_xticks(x_ticks)

# Manually set x-axis tick labels
x_tick_labels = [str(int(x)) for x in x_ticks]
axs.set_xticklabels(x_tick_labels)

# Hide y-axis ticks
axs.set_yticks([])

# Add a black vertical line at x = 16
axs.axvline(x=16, color='black', linestyle='-', linewidth=1)

# Adjust spacing between subplots
plt.tight_layout()

# try:
#   ### TO SAVE OUTPUTS
#   ## Generate the file name dynamically based on the cluster index
#   #figure_path = r"C:\Users\aozgu\Downloads\init65_0.2thresh_overall_clustering_heatmap" + ".eps"
#   figure_path = r"C:\Users\LurLab\Downloads\init65_0.2thresh_overall_clustering_heatmap" + ".png"
#   ## Save the figure
#   plt.savefig(figure_path, format='png', bbox_inches='tight')
# except:
#   pass

# Show the plot
plt.show()


#%%


#%%
# Initialize a list to store cluster lengths
cluster_lengths = []

# Iterate over clusters and append the length of cell IDs for each cluster to the list
for cluster_name, cell_ids in cluster_indices_dict.items():
    cluster_length = len(cell_ids)

    # Append the cluster length to the list
    cluster_lengths.append(cluster_length)

# Now, cluster_lengths is a list containing the lengths of cell IDs for each cluster
# You can access the lengths using indexing (e.g., cluster_lengths[0] for the first cluster)
for i, cluster_length in enumerate(cluster_lengths):
    cluster_name = list(cluster_indices_dict.keys())[i]
    print(f"Length of cell IDs in {cluster_name}: {cluster_length}")

# You can also use cluster_lengths list in your code as needed


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

naive_correct_std_traces = []
naive_incorrect_std_traces = []

naive_correct_ub_traces = []
naive_correct_lb_traces = []
naive_incorrect_ub_traces = []
naive_incorrect_lb_traces = []


correct = C_hm_naive_normalized
incorrect = W_hm_naive_normalized

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
  std = np.std(naive_correct_traces[i], axis=0)
  sem = np.std(naive_correct_traces[i], axis=0) / np.sqrt(len(naive_correct_traces[i]))
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  naive_correct_mean_traces.append(mean)
  naive_correct_sem_traces.append(sem)
  naive_correct_std_traces.append(std)
  naive_correct_ub_traces.append(upper_bound)
  naive_correct_lb_traces.append(lower_bound)

for i in range(len(naive_incorrect_traces)):
  mean = np.mean(naive_incorrect_traces[i], axis=0)
  std = np.std(naive_incorrect_traces[i], axis=0)
  sem = np.std(naive_incorrect_traces[i], axis=0) / np.sqrt(len(naive_incorrect_traces[i]))
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  naive_incorrect_mean_traces.append(mean)
  naive_incorrect_sem_traces.append(sem)
  naive_incorrect_std_traces.append(std)
  naive_incorrect_ub_traces.append(upper_bound)
  naive_incorrect_lb_traces.append(lower_bound)

#### expert
expert_correct_traces = []
expert_incorrect_traces = []

expert_correct_mean_traces = []
expert_incorrect_mean_traces = []

expert_correct_sem_traces = []
expert_incorrect_sem_traces = []

expert_correct_std_traces = []
expert_incorrect_std_traces = []

expert_correct_ub_traces = []
expert_correct_lb_traces = []
expert_incorrect_ub_traces = []
expert_incorrect_lb_traces = []


correct = C_hm_expert_normalized
incorrect = W_hm_expert_normalized

# List of arrays to process with their names
array_names = ["correct", "incorrect"]
arrays_to_process = [correct, incorrect]

fig_width, fig_height = 4, 2  # You can adjust these values as needed

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
  std = np.std(expert_correct_traces[i], axis=0)
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  expert_correct_mean_traces.append(mean)
  expert_correct_sem_traces.append(sem)
  expert_correct_std_traces.append(std)
  expert_correct_ub_traces.append(upper_bound)
  expert_correct_lb_traces.append(lower_bound)

for i in range(len(expert_incorrect_traces)):
  mean = np.mean(expert_incorrect_traces[i],axis=0)
  sem = np.std(expert_incorrect_traces[i], axis=0) / np.sqrt(len(expert_incorrect_traces[i]))
  std = np.std(expert_incorrect_traces[i], axis=0)
  upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
  lower_bound = mean - 1.96 * sem
  expert_incorrect_mean_traces.append(mean)
  expert_incorrect_sem_traces.append(sem)
  expert_incorrect_std_traces.append(std)
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

  # Set the y-axis limits
  plt.ylim(0, 1)

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

# Calculate and store scores for expert traces, all frames
for trace1, trace2 in zip(expert_correct_mean_traces, expert_incorrect_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    expert_scores.append(result)

# # Calculate and store scores for expert traces, preturn frames
# for trace1, trace2 in zip(expert_correct_mean_traces, expert_incorrect_mean_traces):
#     result = average_absolute_difference(trace1[0:15], trace2[0:15])
#     expert_scores_preturn.append(result)

# # Calculate and store scores for expert traces, postturn frames
# for trace1, trace2 in zip(expert_correct_mean_traces, expert_incorrect_mean_traces):
#     result = average_absolute_difference(trace1[15:40], trace2[15:40])
#     expert_scores_postturn.append(result)



# Calculate and store scores for expert traces, all frames
for trace1, trace2 in zip(naive_correct_mean_traces, naive_incorrect_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    naive_scores.append(result)

# # Calculate and store scores for expert traces, preturn frames
# for trace1, trace2 in zip(naive_correct_mean_traces, naive_incorrect_mean_traces):
#     result = average_absolute_difference(trace1[0:15], trace2[0:15])
#     naive_scores_preturn.append(result)

# # Calculate and store scores for expert traces, postturn frames
# for trace1, trace2 in zip(naive_correct_mean_traces, naive_incorrect_mean_traces):
#     result = average_absolute_difference(trace1[15:40], trace2[15:40])
#     naive_scores_postturn.append(result)



# Calculate 95% confidence intervals
ci_expert = 1.96 * np.std(expert_scores) / np.sqrt(len(expert_scores))
ci_naive = 1.96 * np.std(naive_scores) / np.sqrt(len(naive_scores))

# ci_expert_preturn = 1.96 * np.std(expert_scores_preturn) / np.sqrt(len(expert_scores_preturn))
# ci_naive_preturn = 1.96 * np.std(naive_scores_preturn) / np.sqrt(len(naive_scores_preturn))

# ci_expert_postturn = 1.96 * np.std(expert_scores_postturn) / np.sqrt(len(expert_scores_postturn))
# ci_naive_postturn = 1.96 * np.std(naive_scores_postturn) / np.sqrt(len(naive_scores_postturn))

# Plotting as a bar graph with error bars
groups = ['Expert', 'Naive']
average_scores = [np.mean(expert_scores), np.mean(naive_scores)]
confidence_intervals = [ci_expert, ci_naive]

plt.bar(groups, average_scores, yerr=confidence_intervals, capsize=10, color=['blue', 'orange'])
plt.xlabel('Group')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences (Across Traces)')

# Perform t-test
t_stat, p_value = ttest_ind(expert_scores, naive_scores)
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

################# PIE CHART
from scipy.stats import ttest_ind_from_stats
import numpy as np
import matplotlib.pyplot as plt

expert_CI_sems = []
naive_CI_sems = []

# Initialize counts for different categories
significant_expert_greater_count = 0
significant_naive_greater_count = 0
non_significant_count = 0

# Calculate the combined standard errors
for i in range(len(expert_correct_sem_traces)):
    expert_CI_sems.append(np.sqrt(np.mean(expert_correct_sem_traces[i])**2 + np.mean(expert_incorrect_sem_traces[i])**2))
    naive_CI_sems.append(np.sqrt(np.mean(naive_correct_sem_traces[i])**2 + np.mean(naive_incorrect_sem_traces[i])**2))

# For each cluster
for i in range(len(expert_CI_sems)):
    mean_expert = expert_scores[i]
    sem_expert = expert_CI_sems[i]
    #n_expert = 22 #init
    #n_expert = 25 #ftp

    n_expert = 5

    mean_naive = naive_scores[i]
    sem_naive = naive_CI_sems[i]
    #n_naive = 15 #init
    #n_naive = 23 #ftp

    n_naive = 5

    # Perform a two-sample t-test
    t_statistic, p_value = ttest_ind_from_stats(mean_expert, sem_expert, n_expert, mean_naive, sem_naive, n_naive)

    # Output the results
    print("t-statistic:", t_statistic)
    print("p-value:", p_value)

    # Check the significance level
    alpha = 0.05
    if p_value < alpha:
        if t_statistic > 0:
            print("Significant (Expert > Naive)")
            significant_expert_greater_count += 1
        elif t_statistic < 0:
            print("Significant (Expert < Naive)")
            significant_naive_greater_count += 1
    else:
        print("Non-Significant")
        non_significant_count += 1

# Pie chart
labels = ['Significant (Expert > Naive)', 'Significant (Expert < Naive)', 'Non-Significant']
sizes = [significant_expert_greater_count, significant_naive_greater_count, non_significant_count]
colors = ['lightcoral', 'lightskyblue', 'lightgreen']
explode = (0.1, 0.1, 0)  # explode 1st and 2nd slices

# Calculate percentages
total_cases = sum(sizes)
percentages = [(count / total_cases) * 100 for count in sizes]

## Display numbers and percentages
for i, label in enumerate(labels):
    print(f"{label}: {sizes[i]} cases ({percentages[i]:.2f}%)")

plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=140)
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title('Significance and Directionality of Results')

### TO SAVE OUTPUTS
## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\piechart_indiv_cluster_differences" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\piechart_indiv_cluster_differences" + ".eps"
## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()
####################################### INDIV CLUSTERS
# Plotting as a bar graph with individual error bars
indices = np.arange(len(expert_scores))
bar_width = 0.35

# Calculate 95% confidence intervals for each bar
ci_expert = [1.96 * sem / np.sqrt(n) for mean, sem, n in zip(expert_scores, expert_CI_sems, [n_expert]*len(expert_scores))]
ci_naive = [1.96 * sem / np.sqrt(n) for mean, sem, n in zip(naive_scores, naive_CI_sems, [n_naive]*len(naive_scores))]

plt.bar(indices, expert_scores, bar_width, yerr=ci_expert, capsize=10, label='Expert Traces', color='blue', alpha=0.7)
plt.bar(indices + bar_width, naive_scores, bar_width, yerr=ci_naive, capsize=10, label='Naive Traces', color='orange', alpha=0.7)

plt.xlabel('Index of Traces')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences')
plt.legend()
plt.xticks(indices + bar_width / 2, [f'Trace {i+1}' for i in range(len(expert_scores))])

### TO SAVE OUTPUTS
## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\indiv_cluster_differences" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\indiv_cluster_differences" + ".eps"


## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()

# # Plotting as a bar graph
# indices = np.arange(len(expert_scores))
# bar_width = 0.35
# plt.bar(indices, expert_scores, bar_width, label='Expert Traces')
# plt.bar(indices + bar_width, naive_scores, bar_width, label='Naive Traces')
# plt.xlabel('Index of Traces')
# plt.ylabel('Average Absolute Difference')
# plt.title('Comparison of Average Absolute Differences')
# plt.legend()
# plt.xticks(indices + bar_width / 2, [f'Trace {i+1}' for i in range(len(expert_scores))])

# ### TO SAVE OUTPUTS
# ## Generate the file name dynamically based on the cluster index
# figure_path = r"C:\Users\aozgu\Downloads\init65_0.05thresh_CI_clustering_indivcluster_differences" + ".eps"
# ## Save the figure
# plt.savefig(figure_path, format='eps', bbox_inches='tight')

# plt.show()



#%%
# FOR MANUAL ANOVA PRISM
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# Calculate and store scores for naive traces, all frames
print('naive:')
for trace1, trace2 in zip(naive_correct_mean_traces, naive_incorrect_mean_traces):
    result = average_absolute_differences(trace1, trace2)
    result_str = ', '.join(map(str, result))
    print(result_str)

# Calculate and store scores for expert traces, all frames
print('expert:')
for trace1, trace2 in zip(expert_correct_mean_traces, expert_incorrect_mean_traces):
    result = average_absolute_differences(trace1, trace2)
    result_str = ', '.join(map(str, result))
    print(result_str)

#%%
for array in expert_incorrect_mean_traces:
    print(', '.join(map(str, array)))

#%%
naive_correct_mean_traces

#%%
naive_incorrect_mean_traces

#%%
expert_correct_mean_traces

#%%
expert_incorrect_mean_traces

#%%
# left v right TURN

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

left = Lturn_hm_naive_normalized
right = Rturn_hm_naive_normalized

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

left_turn = Lturn_hm_expert_normalized
right_turn = Rturn_hm_expert_normalized

# List of arrays to process with their names
array_names = ["left_turn", "right_turn"]
arrays_to_process = [left_turn, right_turn]

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
        if array_name == "left_turn":
            expert_left_traces.append(cluster_array)
        elif array_name == "right_turn":
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
    plt.plot(naive_left_mean_traces[i], color=dark_red, label='naive_left_turn')
    plt.fill_between(range(len(naive_left_mean_traces[i])), naive_left_lb_traces[i], naive_left_ub_traces[i], alpha=0.2, color=dark_red)

    plt.plot(naive_right_mean_traces[i], color=light_red, label='naive_right_turn')
    plt.fill_between(range(len(naive_right_mean_traces[i])), naive_right_lb_traces[i], naive_right_ub_traces[i], alpha=0.2, color=light_red)

    plt.plot(expert_left_mean_traces[i], color=dark_blue, label='expert_left_turn')
    plt.fill_between(range(len(expert_left_mean_traces[i])), expert_left_lb_traces[i], expert_left_ub_traces[i], alpha=0.2, color=dark_blue)

    plt.plot(expert_right_mean_traces[i], color=light_blue, label='expert_right_turn')
    plt.fill_between(range(len(expert_right_mean_traces[i])), expert_right_lb_traces[i], expert_right_ub_traces[i], alpha=0.2, color=light_blue)


    # Add labels and a legend for the current cluster
    plt.xlabel('Frames')
    plt.ylabel('Average Fluorescence Intensity')
    plt.legend()

    plt.axvline(x=16, color='black', linestyle='--', label='turnframe')
    plt.legend(fontsize='small')

    # Set the y-axis limits
    plt.ylim(0, 1)

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
expert_scores_postturn = []

naive_scores_preturn = []
naive_scores_postturn = []

# Calculate and store scores for expert traces, all frames
for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    expert_scores.append(result)

# # Calculate and store scores for expert traces, preturn frames
# for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
#     result = average_absolute_difference(trace1[0:15], trace2[0:15])
#     expert_scores_preturn.append(result)

# # Calculate and store scores for expert traces, postturn frames
# for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
#     result = average_absolute_difference(trace1[15:40], trace2[15:40])
#     expert_scores_postturn.append(result)

# Calculate and store scores for naive traces, all frames
for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    naive_scores.append(result)

# # Calculate and store scores for naive traces, preturn frames
# for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
#     result = average_absolute_difference(trace1[0:15], trace2[0:15])
#     naive_scores_preturn.append(result)

# # Calculate and store scores for naive traces, postturn frames
# for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
#     result = average_absolute_difference(trace1[15:40], trace2[15:40])
#     naive_scores_postturn.append(result)

# Calculate 95% confidence intervals
ci_expert = 1.96 * np.std(expert_scores) / np.sqrt(len(expert_scores))
ci_naive = 1.96 * np.std(naive_scores) / np.sqrt(len(naive_scores))

# ci_expert_preturn = 1.96 * np.std(expert_scores_preturn) / np.sqrt(len(expert_scores_preturn))
# ci_naive_preturn = 1.96 * np.std(naive_scores_preturn) / np.sqrt(len(naive_scores_preturn))

# ci_expert_postturn = 1.96 * np.std(expert_scores_postturn) / np.sqrt(len(expert_scores_postturn))
# ci_naive_postturn = 1.96 * np.std(naive_scores_postturn) / np.sqrt(len(naive_scores_postturn))

# Plotting as a bar graph with error bars
groups = ['Expert', 'Naive']
average_scores = [np.mean(expert_scores), np.mean(naive_scores)]
confidence_intervals = [ci_expert, ci_naive]

plt.bar(groups, average_scores, yerr=confidence_intervals, capsize=10, color=['blue', 'orange'])
plt.xlabel('Group')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences (Across Traces)')

# Perform t-test
t_stat, p_value = ttest_ind(expert_scores, naive_scores)
print(f'T-test p-value: {p_value}')

if p_value < 0.05:
    print('The difference between the groups is statistically significant.')
else:
    print('The difference between the groups is not statistically significant.')

# ### TO SAVE OUTPUTS
# ## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\groupedcluster_differences" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\groupedcluster_differences" + ".eps"
# ## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()

######################## PIE CHART
from scipy.stats import ttest_ind_from_stats
import numpy as np
import matplotlib.pyplot as plt

expert_turn_sems = []
naive_turn_sems = []

# Initialize counts for different categories
significant_expert_greater_count = 0
significant_naive_greater_count = 0
non_significant_count = 0

# Calculate the combined standard errors
for i in range(len(expert_left_sem_traces)):
    expert_turn_sems.append(np.sqrt(np.mean(expert_left_sem_traces[i])**2 + np.mean(expert_right_sem_traces[i])**2))
    naive_turn_sems.append(np.sqrt(np.mean(naive_left_sem_traces[i])**2 + np.mean(naive_right_sem_traces[i])**2))

# For each cluster
for i in range(len(expert_turn_sems)):
    mean_expert = expert_scores[i]
    sem_expert = expert_turn_sems[i]
    #n_expert = 22
    #n_expert = 25

    n_expert = 5

    mean_naive = naive_scores[i]
    sem_naive = naive_turn_sems[i]
    #n_naive = 15
    #n_naive = 23

    n_naive = 5

    # Perform a two-sample t-test
    t_statistic, p_value = ttest_ind_from_stats(mean_expert, sem_expert, n_expert, mean_naive, sem_naive, n_naive)

    # Output the results
    print("t-statistic:", t_statistic)
    print("p-value:", p_value)

    # Check the significance level
    alpha = 0.05
    if p_value < alpha:
        if t_statistic > 0:
            print("Significant (Expert > Naive)")
            significant_expert_greater_count += 1
        elif t_statistic < 0:
            print("Significant (Expert < Naive)")
            significant_naive_greater_count += 1
    else:
        print("Non-Significant")
        non_significant_count += 1

# Pie chart
labels = ['Significant (Expert > Naive)', 'Significant (Expert < Naive)', 'Non-Significant']
sizes = [significant_expert_greater_count, significant_naive_greater_count, non_significant_count]
colors = ['lightcoral', 'lightskyblue', 'lightgreen']
explode = (0.1, 0.1, 0)  # explode 1st and 2nd slices

# Calculate percentages
total_cases = sum(sizes)
percentages = [(count / total_cases) * 100 for count in sizes]

# Display numbers and percentages
for i, label in enumerate(labels):
    print(f"{label}: {sizes[i]} cases ({percentages[i]:.2f}%)")

plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=140)
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title('Significance and Directionality of Results')

# ### TO SAVE OUTPUTS
# ## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\piechart_indiv_cluster_differences" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\indiv_cluster_differences" + ".eps"
# ## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()
####################################### INDIV CLUSTERS
# Plotting as a bar graph with individual error bars
indices = np.arange(len(expert_scores))
bar_width = 0.35

# Calculate 95% confidence intervals for each bar
ci_expert = [1.96 * sem / np.sqrt(n) for mean, sem, n in zip(expert_scores, expert_turn_sems, [n_expert]*len(expert_scores))]
ci_naive = [1.96 * sem / np.sqrt(n) for mean, sem, n in zip(naive_scores, naive_turn_sems, [n_naive]*len(naive_scores))]

plt.bar(indices, expert_scores, bar_width, yerr=ci_expert, capsize=10, label='Expert Traces', color='blue', alpha=0.7)
plt.bar(indices + bar_width, naive_scores, bar_width, yerr=ci_naive, capsize=10, label='Naive Traces', color='orange', alpha=0.7)

plt.xlabel('Index of Traces')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences')
plt.legend()
plt.xticks(indices + bar_width / 2, [f'Trace {i+1}' for i in range(len(expert_scores))])

# ### TO SAVE OUTPUTS
# ## Generate the file name dynamically based on the cluster index
figure_path = r"C:\Users\aozgu\Downloads\indiv_cluster_differences" + ".eps"
#figure_path = r"C:\Users\LurLab\Downloads\piechart_indiv_cluster_differences" + ".eps"
# ## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()


#%%
import os
import csv

# Function to save traces to CSV
def save_traces_to_csv(traces, trace_type, folder_path):
    csv_prefix = f"{trace_type}_cluster"
    csv_suffix = ".csv"

    # Iterate over clusters
    for i, cluster_traces in enumerate(traces):
        # Build the full file path
        csv_file = f"{csv_prefix}{i+1}{csv_suffix}"
        file_path = os.path.join(folder_path, csv_file)

        # Writing to CSV file
        with open(file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(cluster_traces)

        print(f"CSV file '{file_path}' has been created and saved to the specified folder.")


# Specify the folder
folder_path = r'C:\Users\LurLab\Downloads'

# Create the folder if it doesn't exist
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

# Save traces for different types
save_traces_to_csv(naive_right_traces, "naive_right_turn_traces", folder_path)
save_traces_to_csv(naive_left_traces, "naive_left_turn_traces", folder_path)
save_traces_to_csv(expert_right_traces, "expert_right_turn_traces", folder_path)
save_traces_to_csv(expert_left_traces, "expert_left_turn_traces", folder_path)


#%%
# FOR MANUAL ANOVA PRISM
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# Calculate and store scores for naive traces, all frames
print('naive:')
for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
    result = average_absolute_differences(trace1, trace2)
    result_str = ', '.join(map(str, result))
    print(result_str)

# Calculate and store scores for expert traces, all frames
print('expert:')
for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
    result = average_absolute_differences(trace1, trace2)
    result_str = ', '.join(map(str, result))
    print(result_str)

#%%
for array in expert_right_mean_traces:
    print(', '.join(map(str, array)))

#%%
naive_left_mean_traces

#%%
naive_right_mean_traces

#%%
expert_left_mean_traces

#%%
expert_right_mean_traces

#%%
#FOR VISUAL

#### naive
import numpy as np
import matplotlib.pyplot as plt

naive_left_traces = []
naive_right_traces = []

naive_left_mean_traces = []
naive_right_mean_traces = []

naive_left_sem_traces = []
naive_right_sem_traces = []

naive_left_std_traces = []
naive_right_std_traces = []

naive_left_ub_traces = []
naive_left_lb_traces = []
naive_right_ub_traces = []
naive_right_lb_traces = []

left = Lvis_hm_naive_normalized
right = Rvis_hm_naive_normalized

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

naive_left_traces_n = len(naive_left_traces)
naive_right_traces_n = len(naive_right_traces)

for i in range(len(naive_left_traces)):
    mean = np.mean(naive_left_traces[i], axis=0)
    sem = np.std(naive_left_traces[i], axis=0) / np.sqrt(len(naive_left_traces[i]))
    std = np.std(naive_left_traces[i], axis=0)
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    naive_left_mean_traces.append(mean)
    naive_left_sem_traces.append(sem)
    naive_left_std_traces.append(std)
    naive_left_ub_traces.append(upper_bound)
    naive_left_lb_traces.append(lower_bound)

for i in range(len(naive_right_traces)):
    mean = np.mean(naive_right_traces[i], axis=0)
    sem = np.std(naive_right_traces[i], axis=0) / np.sqrt(len(naive_right_traces[i]))
    std = np.std(naive_right_traces[i], axis=0)
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    naive_right_mean_traces.append(mean)
    naive_right_sem_traces.append(sem)
    naive_right_std_traces.append(std)
    naive_right_ub_traces.append(upper_bound)
    naive_right_lb_traces.append(lower_bound)

#### expert
expert_left_traces = []
expert_right_traces = []

expert_left_mean_traces = []
expert_right_mean_traces = []

expert_left_sem_traces = []
expert_right_sem_traces = []

expert_left_std_traces = []
expert_right_std_traces = []

expert_left_ub_traces = []
expert_left_lb_traces = []
expert_right_ub_traces = []
expert_right_lb_traces = []

left = Lvis_hm_expert_normalized
right = Rvis_hm_expert_normalized

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

expert_left_traces_n = len(expert_left_traces)
expert_right_traces_n = len(expert_right_traces)

for i in range(len(expert_left_traces)):
    mean = np.mean(expert_left_traces[i], axis=0)
    sem = np.std(expert_left_traces[i], axis=0) / np.sqrt(len(expert_left_traces[i]))
    std = np.std(expert_left_traces[i], axis=0)
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    expert_left_mean_traces.append(mean)
    expert_left_sem_traces.append(sem)
    expert_left_std_traces.append(std)
    expert_left_ub_traces.append(upper_bound)
    expert_left_lb_traces.append(lower_bound)

for i in range(len(expert_right_traces)):
    mean = np.mean(expert_right_traces[i], axis=0)
    sem = np.std(expert_right_traces[i], axis=0) / np.sqrt(len(expert_right_traces[i]))
    std = np.std(expert_right_traces[i], axis=0)
    upper_bound = mean + 1.96 * sem  # 1.96 corresponds to the Z-score for a 95% CI
    lower_bound = mean - 1.96 * sem
    expert_right_mean_traces.append(mean)
    expert_right_sem_traces.append(sem)
    expert_right_std_traces.append(std)
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
    plt.plot(naive_left_mean_traces[i], color=dark_red, label='naive_left_visual')
    plt.fill_between(range(len(naive_left_mean_traces[i])), naive_left_lb_traces[i], naive_left_ub_traces[i], alpha=0.2, color=dark_red)

    plt.plot(naive_right_mean_traces[i], color=light_red, label='naive_right_visual')
    plt.fill_between(range(len(naive_right_mean_traces[i])), naive_right_lb_traces[i], naive_right_ub_traces[i], alpha=0.2, color=light_red)

    plt.plot(expert_left_mean_traces[i], color=dark_blue, label='expert_left_visual')
    plt.fill_between(range(len(expert_left_mean_traces[i])), expert_left_lb_traces[i], expert_left_ub_traces[i], alpha=0.2, color=dark_blue)

    plt.plot(expert_right_mean_traces[i], color=light_blue, label='expert_right_visual')
    plt.fill_between(range(len(expert_right_mean_traces[i])), expert_right_lb_traces[i], expert_right_ub_traces[i], alpha=0.2, color=light_blue)


    # Add labels and a legend for the current cluster
    plt.xlabel('Frames')
    plt.ylabel('Average Fluorescence Intensity')
    plt.legend()

    plt.axvline(x=16, color='black', linestyle='--', label='turnframe')
    plt.legend(fontsize='small')

    # Set the y-axis limits
    plt.ylim(0, 1)

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
from scipy import stats

# Lists to store the results
expert_scores = []
naive_scores = []

expert_scores_preturn = []
expert_scores_postturn = []

naive_scores_preturn = []
naive_scores_postturn = []

# Calculate and store scores for expert traces, all frames
for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    expert_scores.append(result)

# # Calculate and store scores for expert traces, preturn frames
# for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
#     result = average_absolute_difference(trace1[0:15], trace2[0:15])
#     expert_scores_preturn.append(result)

# # Calculate and store scores for expert traces, postturn frames
# for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
#     result = average_absolute_difference(trace1[15:40], trace2[15:40])
#     expert_scores_postturn.append(result)

# Calculate and store scores for naive traces, all frames
for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
    result = average_absolute_difference(trace1, trace2)
    naive_scores.append(result)

# # Calculate and store scores for naive traces, preturn frames
# for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
#     result = average_absolute_difference(trace1[0:15], trace2[0:15])
#     naive_scores_preturn.append(result)

# # Calculate and store scores for naive traces, postturn frames
# for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
#     result = average_absolute_difference(trace1[15:40], trace2[15:40])
#     naive_scores_postturn.append(result)

# Calculate 95% confidence intervals
ci_expert = 1.96 * np.std(expert_scores) / np.sqrt(len(expert_scores))
ci_naive = 1.96 * np.std(naive_scores) / np.sqrt(len(naive_scores))

# ci_expert_preturn = 1.96 * np.std(expert_scores_preturn) / np.sqrt(len(expert_scores_preturn))
# ci_naive_preturn = 1.96 * np.std(naive_scores_preturn) / np.sqrt(len(naive_scores_preturn))

# ci_expert_postturn = 1.96 * np.std(expert_scores_postturn) / np.sqrt(len(expert_scores_postturn))
# ci_naive_postturn = 1.96 * np.std(naive_scores_postturn) / np.sqrt(len(naive_scores_postturn))

# Plotting as a bar graph with error bars
groups = ['Expert', 'Naive']
average_scores = [np.mean(expert_scores), np.mean(naive_scores)]
confidence_intervals = [ci_expert, ci_naive]

plt.bar(groups, average_scores, yerr=confidence_intervals, capsize=10, color=['blue', 'orange'])
plt.xlabel('Group')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences (Across Traces)')

# Perform t-test
t_stat, p_value = ttest_ind(expert_scores, naive_scores)
print(f'T-test p-value: {p_value}')

if p_value < 0.05:
    print('The difference between the groups is statistically significant.')
else:
    print('The difference between the groups is not statistically significant.')

try:
  ### TO SAVE OUTPUTS
  ## Generate the file name dynamically based on the cluster index
  figure_path = r"C:\Users\aozgu\Downloads\groupedcluster_differences" + ".eps"
  #figure_path = r"C:\Users\LurLab\Downloads\groupedcluster_differences" + ".eps"
  ## Save the figure
  plt.savefig(figure_path, format='eps', bbox_inches='tight')
except:
  pass

plt.show()

######################## PIE CHART
from scipy.stats import ttest_ind_from_stats
import numpy as np
import matplotlib.pyplot as plt

expert_vis_sems = []
naive_vis_sems = []

# Initialize counts for different categories
significant_expert_greater_count = 0
significant_naive_greater_count = 0
non_significant_count = 0

# Calculate the combined standard errors
for i in range(len(expert_left_sem_traces)):
    expert_vis_sems.append(np.sqrt(np.mean(expert_left_sem_traces[i])**2 + np.mean(expert_right_sem_traces[i])**2))
    naive_vis_sems.append(np.sqrt(np.mean(naive_left_sem_traces[i])**2 + np.mean(naive_right_sem_traces[i])**2))

# For each cluster
for i in range(len(expert_vis_sems)):
    mean_expert = expert_scores[i]
    sem_expert = expert_vis_sems[i]
    #n_expert = 22
    #n_expert = 25

    n_expert = 5

    mean_naive = naive_scores[i]
    sem_naive = naive_vis_sems[i]
    #n_naive = 15
    #n_naive = 23

    n_naive = 5

    # Perform a two-sample t-test
    t_statistic, p_value = ttest_ind_from_stats(mean_expert, sem_expert, n_expert, mean_naive, sem_naive, n_naive)
    # Output the results
    print("t-statistic:", t_statistic)
    print("p-value:", p_value)


    # Check the significance level
    alpha = 0.05
    if p_value < alpha:
        if t_statistic > 0:
            print("Significant (Expert > Naive)")
            significant_expert_greater_count += 1
        elif t_statistic < 0:
            print("Significant (Expert < Naive)")
            significant_naive_greater_count += 1
    else:
        print("Non-Significant")
        non_significant_count += 1

# Pie chart
labels = ['Significant (Expert > Naive)', 'Significant (Expert < Naive)', 'Non-Significant']
sizes = [significant_expert_greater_count, significant_naive_greater_count, non_significant_count]
colors = ['lightcoral', 'lightskyblue', 'lightgreen']
explode = (0.1, 0.1, 0)  # explode 1st and 2nd slices

# Calculate percentages
total_cases = sum(sizes)
percentages = [(count / total_cases) * 100 for count in sizes]

# Display numbers and percentages
for i, label in enumerate(labels):
    print(f"{label}: {sizes[i]} cases ({percentages[i]:.2f}%)")

plt.pie(sizes, explode=explode, labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=140)
plt.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title('Significance and Directionality of Results')

try:
  ### TO SAVE OUTPUTS
  ## Generate the file name dynamically based on the cluster index
  figure_path = r"C:\Users\aozgu\Downloads\piechart_indiv_cluster_differences" + ".eps"
  #figure_path = r"C:\Users\LurLab\Downloads\indiv_cluster_differences" + ".eps"
  ## Save the figure
  plt.savefig(figure_path, format='eps', bbox_inches='tight')
except:
  pass

plt.show()
####################################### INDIV CLUSTERS
# Plotting as a bar graph with individual error bars
indices = np.arange(len(expert_scores))
bar_width = 0.35

# Calculate 95% confidence intervals for each bar
ci_expert = [1.96 * sem / np.sqrt(n) for mean, sem, n in zip(expert_scores, expert_vis_sems, [n_expert]*len(expert_scores))]
ci_naive = [1.96 * sem / np.sqrt(n) for mean, sem, n in zip(naive_scores, naive_vis_sems, [n_naive]*len(naive_scores))]

plt.bar(indices, expert_scores, bar_width, yerr=ci_expert, capsize=10, label='Expert Traces', color='blue', alpha=0.7)
plt.bar(indices + bar_width, naive_scores, bar_width, yerr=ci_naive, capsize=10, label='Naive Traces', color='orange', alpha=0.7)

plt.xlabel('Index of Traces')
plt.ylabel('Average Absolute Difference')
plt.title('Comparison of Average Absolute Differences')
plt.legend()
plt.xticks(indices + bar_width / 2, [f'Trace {i+1}' for i in range(len(expert_scores))])

try:
  ### TO SAVE OUTPUTS
  ## Generate the file name dynamically based on the cluster index
  figure_path = r"C:\Users\aozgu\Downloads\indiv_cluster_differences" + ".eps"
  #figure_path = r"C:\Users\LurLab\Downloads\piechart_indiv_cluster_differences" + ".eps"
  ## Save the figure
  plt.savefig(figure_path, format='eps', bbox_inches='tight')
except:
  pass

plt.show()


#%%
# FOR MANUAL ANOVA PRISM
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# Calculate and store scores for naive traces, all frames
print('naive:')
for trace1, trace2 in zip(naive_left_mean_traces, naive_right_mean_traces):
    result = average_absolute_differences(trace1, trace2)
    result_str = ', '.join(map(str, result))
    print(result_str)

# Calculate and store scores for expert traces, all frames
print('expert:')
for trace1, trace2 in zip(expert_left_mean_traces, expert_right_mean_traces):
    result = average_absolute_differences(trace1, trace2)
    result_str = ', '.join(map(str, result))
    print(result_str)

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
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_turn_frames'][i]]) ###MAIN
                tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-15:base_params['clean_trial_turn_frames'][i]+25])
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_stim_start_frames'][i]+20])

                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-20:base_params['clean_trial_turn_frames'][i]])

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
def get_parameters_Lturn(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    # L turn
    clean_df = clean_df[(clean_df['correctness'] == 'correct') & (clean_df['direction'] == 'L') |
                       (clean_df['correctness'] == 'wrong') & (clean_df['direction'] == 'R')]

    # L vis
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

    # use for turn decoding
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

      # ## use only for visual stimuli analysis
      # if clean_df['direction'].iloc[i]=='L':
      #  clean_trial_direction.append(0)
      # if clean_df['direction'].iloc[i]=='R':
      #   clean_trial_direction.append(1)

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
def get_parameters_Rturn(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    # R turn
    clean_df = clean_df[(clean_df['correctness'] == 'correct') & (clean_df['direction'] == 'R') |
                       (clean_df['correctness'] == 'wrong') & (clean_df['direction'] == 'L')]

    # R vis
    #clean_df = clean_df[clean_df['direction']=='R']


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

    # use for turn decoding
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

      # ## use only for visual stimuli analysis
      # if clean_df['direction'].iloc[i]=='L':
      #  clean_trial_direction.append(0)
      # if clean_df['direction'].iloc[i]=='R':
      #   clean_trial_direction.append(1)

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
def get_parameters_Lvis(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    # L turn
    #clean_df = clean_df[(clean_df['correctness'] == 'correct') & (clean_df['direction'] == 'L') |
    #                  (clean_df['correctness'] == 'wrong') & (clean_df['direction'] == 'R')]

    # L vis
    clean_df = clean_df[clean_df['direction']=='L']


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
def get_parameters_Rvis(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    ###del when not using all trials
    #clean_df = params
    #a = 0

    # R turn
    #clean_df = clean_df[(clean_df['correctness'] == 'correct') & (clean_df['direction'] == 'R') |
    #                   (clean_df['correctness'] == 'wrong') & (clean_df['direction'] == 'L')]

    # R vis
    clean_df = clean_df[clean_df['direction']=='R']


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

