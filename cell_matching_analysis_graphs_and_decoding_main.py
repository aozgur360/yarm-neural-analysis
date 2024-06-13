#%%
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import sem

#mice_ids = ['m21r']
mice_ids = ['m2r','m9r','m11n','m12lr','m13l','m14r','m15n','m16l','m17r','m19n','m21r','m22l','m23r','m25lr','m26r','m29n', 'm30lr', 'm31r','m32l','m33r','m34n','m36rr','m37l','m38n','m39rr','m40l', 'm41r','min4r','sn','snr'] #main mice except #m27l
results_dict = {}  # Create an empty dictionary to store the results for each mouse_id

for mouse_id in mice_ids:
    try:
        file_path = rf'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\all_cell_matching_240130\init_e1_e2\{mouse_id}\mappings.pkl'
        #file_path = rf'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\crossreg\init_e1_e2_extras_240216\{mouse_id}\mappings.pkl'
        mappings = pd.read_pickle(file_path)

        print(mouse_id)
        mappings['non_nan_count_c1_c2'] = mappings.apply(lambda row: (not pd.isna(row[mappings.columns[1]]) and not pd.isna(row[mappings.columns[2]])), axis=1)
        count_c1_c2_matched = mappings['non_nan_count_c1_c2'].sum()

        # Store the results in the dictionary
        results_dict[mouse_id] = {
            'count_c1_c2_matched': count_c1_c2_matched,
            'count_c1_c2_unmatched': mappings.apply(lambda row: not pd.isna(row[mappings.columns[2]]) and pd.isna(row[mappings.columns[1]]), axis=1).sum(),
            'count_c2_all': mappings[mappings.columns[2]].count(),
            'percent_active_old_cells_c2': count_c1_c2_matched / mappings[mappings.columns[2]].count(),
            'percent_active_new_cells_c2': mappings.apply(lambda row: not pd.isna(row[mappings.columns[2]]) and pd.isna(row[mappings.columns[1]]), axis=1).sum() / mappings[mappings.columns[2]].count(),
            'c2_matched_cells': mappings.loc[~pd.isna(mappings[mappings.columns[1]]) & ~pd.isna(mappings[mappings.columns[2]]), mappings.columns[2]].astype(int).sort_values().tolist(),
            'c1_matched_cells': mappings.loc[~pd.isna(mappings[mappings.columns[1]]) & ~pd.isna(mappings[mappings.columns[2]]), mappings.columns[1]].astype(int).sort_values().tolist(),
            'c2_unmatched_cells': mappings.loc[pd.isna(mappings[mappings.columns[1]]) & ~pd.isna(mappings[mappings.columns[2]]), mappings.columns[2]].astype(int).sort_values().tolist(),
        }

        # Print or use the results as needed
        print(results_dict[mouse_id])
        print("..........")
    except:
        pass

#stats
# anova to test if there is a sig dif between perccent active old cells c2 and active new cells c2
from scipy.stats import f_oneway
# Extracting data for ANOVA test
percent_active_old_cells = []
percent_active_new_cells = []

for mouse_id in mice_ids:
    if mouse_id in results_dict:
        percent_active_old_cells.append(results_dict[mouse_id]['percent_active_old_cells_c2'])
        percent_active_new_cells.append(results_dict[mouse_id]['percent_active_new_cells_c2'])
# Perform ANOVA test if there are at least two groups
if len(percent_active_old_cells) >= 2 and len(percent_active_new_cells) >= 2:
    statistic, p_value = f_oneway(percent_active_old_cells, percent_active_new_cells)

    # Print the results
    print(f"ANOVA Test: F-statistic = {statistic}, p-value = {p_value}")
else:
    print("Insufficient data for ANOVA test.")


### graphs
#graph for individual mice cell change percentages into c2
# Extracting data for bar graph
mouse_ids = list(results_dict.keys())
percent_active_old_cells = [results_dict[mouse_id]['percent_active_old_cells_c2'] for mouse_id in mouse_ids]
percent_active_new_cells = [results_dict[mouse_id]['percent_active_new_cells_c2'] for mouse_id in mouse_ids]

# Plotting bar graph
bar_width = 0.35
index = np.arange(len(mouse_ids))

fig, ax = plt.subplots()
bar1 = ax.bar(index, percent_active_old_cells, bar_width, label='Percent Active Old Cells')
bar2 = ax.bar(index + bar_width, percent_active_new_cells, bar_width, label='Percent Active New Cells')

ax.set_xlabel('Mouse IDs')
ax.set_ylabel('Percentage')
ax.set_title('Percentage of Active Old and New Cells for Each Mouse')
ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(mouse_ids)
ax.legend()

# graph for average percent active new cells across all mice for pie chart
# Calculate average percent active new cells across all mice for pie chart
average_percent_active_new_cells = np.mean(percent_active_new_cells)

### plotting mean
# Calculate mean and confidence interval for percent active new cells
mean_percent_active_new_cells = np.mean(percent_active_new_cells)
ci_percent_active_new_cells = sem(percent_active_new_cells) * 1.96  # 95% confidence interval

# Calculate mean and confidence interval for percent active old cells
mean_percent_active_old_cells = np.mean(percent_active_old_cells)
ci_percent_active_old_cells = sem(percent_active_old_cells) * 1.96  # 95% confidence interval

# Plotting bar graph with T-shaped error bars
fig, ax = plt.subplots()

# Define bar width and positions
bar_width = 0.35
bar_positions = np.arange(2)

# Plot bars
bar2 = ax.bar(bar_positions[0], mean_percent_active_new_cells, bar_width, label='Active New Cells')
bar1 = ax.bar(bar_positions[1], mean_percent_active_old_cells, bar_width, label='Active Old Cells')

# Plot error bars
ax.errorbar(bar_positions[0], mean_percent_active_new_cells, yerr=ci_percent_active_new_cells, fmt='|', color='k', capsize=5, capthick=2)
ax.errorbar(bar_positions[1], mean_percent_active_old_cells, yerr=ci_percent_active_old_cells, fmt='|', color='k', capsize=5, capthick=2)

# Set labels, title, and legend
ax.set_xticks(bar_positions)
ax.set_xticklabels(['Active New Cells', 'Active Old Cells'])
ax.set_ylabel('Mean Percentage with 95% CI')
ax.set_title('Mean Percentage of Active New and Old Cells with 95% CI')
ax.legend()

### TO SAVE OUTPUTS
## Generate the file name dynamically based on the cluster index
#figure_path = r"C:\Users\aozgu\Downloads\cluster" + f'{i+1}' + ".svg"
figure_path = r"C:\Users\LurLab\Downloads\newvsoldcells" + ".eps"
## Save the figure
plt.savefig(figure_path, format='eps', bbox_inches='tight')

plt.show()

# Plotting pie chart
fig, ax = plt.subplots()
labels = ['Active New Cells', 'Active Same Cells']
sizes = [average_percent_active_new_cells, 1 - average_percent_active_new_cells]
colors = ['lightcoral', 'lightskyblue']

ax.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
ax.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.show()
########################
percent_active_old_cells_str = ', '.join([str(results_dict[mouse_id]['percent_active_old_cells_c2']) for mouse_id in mouse_ids])

print('percent active old cells per mouse')
print(percent_active_old_cells_str)


#%%
for mouse_id in mice_ids:
    if mouse_id in results_dict:
        print(f"Mouse ID: {mouse_id}")
        print(f"Percentage of Active New Cells: {results_dict[mouse_id]['percent_active_new_cells_c2'] * 100:.2f}%")
        print(f"Percentage of Active Old Cells: {results_dict[mouse_id]['percent_active_old_cells_c2'] * 100:.2f}%")
        print("----------")


#%%
percent_active_old_cells_str = ', '.join([str(results_dict[mouse_id]['percent_active_old_cells_c2']) for mouse_id in mouse_ids])
print(percent_active_old_cells_str)


#%%


#%%
c2_matched_cells_per_mouse = []
c2_unmatched_cells_per_mouse = []

for mouse_id in mice_ids:
    if mouse_id in results_dict:
        c2_matched_cells_per_mouse.append(results_dict[mouse_id]['c2_matched_cells'])
        c2_unmatched_cells_per_mouse.append(results_dict[mouse_id]['c2_unmatched_cells'])

# Now, c2_matched_cells_per_mouse and c2_unmatched_cells_per_mouse are lists of lists where each inner list corresponds to the 'c2_matched_cells' and 'c2_unmatched_cells' for a specific mouse, respectively.
# You can access them using indexes, for example: c2_matched_cells_per_mouse[0] for the first mouse.

# Print or use the lists of lists as needed
print("c2_matched_cells_per_mouse:")
print(c2_matched_cells_per_mouse)

print("\nc2_unmatched_cells_per_mouse:")
print(c2_unmatched_cells_per_mouse)


#%%
c1_matched_cells_per_mouse = []

for mouse_id in mice_ids:
    print(mouse_id)
    if mouse_id in results_dict:
        c1_matched_cells_per_mouse.append(results_dict[mouse_id]['c1_matched_cells'])

# Now, c2_matched_cells_per_mouse and c2_unmatched_cells_per_mouse are lists of lists where each inner list corresponds to the 'c2_matched_cells' and 'c2_unmatched_cells' for a specific mouse, respectively.
# You can access them using indexes, for example: c2_matched_cells_per_mouse[0] for the first mouse.

# Print or use the lists of lists as needed
print("c1_matched_cells_per_mouse:")
print(c1_matched_cells_per_mouse)

#%%
print(c1_matched_cells_per_mouse[3])

#%%
c1_matched_cells_per_mouse[1]

#%%
results_dict['m12lr']['c1_matched_cells']

#%%
results_dict['m32l']['c2_unmatched_cells']

#%%
# FOR DECODING ON THE 2ND OF THE PAIR

#CURRENTLY USING TF_PARAMS_PATH
stage = 's2' #s2 for init
#0 for turn decoding, 2 for visual stimuli decoding, 1 for correctness decoding
trial_type = 1
# 0 for ALL CELLS, 1 for matched cells, 2 for unmatched cells
#select_cells = 0

#0 for no, 1 for yes
  #0 for no, 1 for yes
if stage in ['s2', 's5']:
    init_only = 1
else:
    init_only = 0


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

TF_params_path = r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params\sessions_for_cell_matching_init_e1_e2_240130.csv"
#TF_params_path = r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params\sessions_for_cell_matching_init_e1_e2_extras.csv"

TF_params = pd.read_csv(TF_params_path)
# Initialize a list to store sessionID values
selected_sessions = []
selected_sessions_expert = []
selected_sessions_naive = []
mice_accuracies_expert = []
mice_accuracies_naive = []

# for all sessions
# if index % 2 == 0 , decoding for e2
# if index % 2 != 0, decoding for e1

for index, row in TF_params.iterrows():
  # FOR ALL SESSIONS WITH EVEN INDICES ###############
    if index % 2 == 0:
    # Append the value in the 'sessionID' column to the list
      selected_sessions.append(row[2])

# Print the list of selected sessionIDs
#print(selected_sessions)
#########################################
temp = selected_sessions
# Filter out 'NaN' strings and NaN floats from temp
temp = [item for item in temp if (not isinstance(item, float) and item.lower() != 'nan') or (isinstance(item, float) and not math.isnan(item))]

#temp.sort()
# Print the sorted list
print(temp)

### IF USING DESKTOP THEN UNCOMMENT BELOW:
#temp = [path.replace('/content/drive/MyDrive', 'E:/google_drive') for path in temp]

#select_cells_values = [0, 1, 2]
select_cells_values = [0]

# Loop over select_cells values
for select_cells in select_cells_values:


  behavior = []
  folders_analyzed = []
  #temp.sort()
  param_paths = []
  cal_paths = []
  cal_paths_raw = []
  rejected_cells_paths = []
  all_results = []
  score_means = []

  for file in temp:
      #print(file)
      param_paths.append(glob.glob(file + '/trial_parameters' + '*.csv'))
      cal_paths.append(glob.glob(file + '/spikerate' + '*.csv'))
      #cal_paths_raw.append(glob.glob(file + '/calcium' + '*.csv'))
      rejected_cells_paths.append(glob.glob(file + '/rejected' + '*.csv'))

  for i in range(len(cal_paths)):
  #for i in range(4):
  #for i in range(3,4):
    try:

      current_folder = temp[i]
      cal = pd.read_csv(cal_paths[i][0])
      #cal_raw = pd.read_csv(cal_paths_raw[i][0])
      params = pd.read_csv(param_paths[i][0])

      #if len(cal.columns) != len(cal_raw.columns):
      #  print('raw calcium and spike cell count mismatch, skipping')
      #  continue

      reject = pd.read_csv(rejected_cells_paths[i][0])
      reject = [item for sublist in reject.values.tolist() for item in sublist]


      del cal[cal.columns[0]]
      # #blood vessels remove
      #cal.drop(cal.columns[reject], axis=1, inplace=True)
      #reset column indices / cell ids
      #cal = cal.transpose()
      #cal = cal.reset_index(drop=True)
      #cal = cal.transpose()

      if select_cells == 0:
        cal.drop(cal.columns[reject], axis=1, inplace=True)
        #reset column indices / cell ids
        cal = cal.transpose()
        cal = cal.reset_index(drop=True)
        cal = cal.transpose()
      if select_cells == 1:
      # decide what cells to use in analysis
        c2_matched_cells_no_rejects = [num for num in c2_matched_cells_per_mouse[i] if num not in reject]
        cal = cal.iloc[:, c2_matched_cells_no_rejects]
      if select_cells == 2:
        c2_unmatched_cells_no_rejects = [num for num in c2_unmatched_cells_per_mouse[i] if num not in reject]
        cal = cal.iloc[:, c2_unmatched_cells_no_rejects]

      parameters = get_parameters(cal,params)
      spike_sums = get_spikerate_sums_per_cell_per_trial(cal,parameters)
      method = spike_sums

      num_incorrect_use = parameters['num_incorrect_trials'] #use max
      num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)
      if num_trials_tested < 10:
        num_incorrect_use = int(0.2 * parameters['ntrials']) #0.2
        num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)

      if num_trials_tested >= 10:
        try:
          svm_results = get_results(method)
          final_results = get_final_stats(svm_results)
        except:
        #  continue
          num_incorrect_use = int(0.1 * parameters['ntrials']) #0.2
          num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)
          svm_results = get_results(method)
          final_results = get_final_stats(svm_results)
      else:
        continue

      if init_only == 0:
        current_behav = (params['correctness'] == 'correct').sum() / len(params)
        behavior.append((params['correctness'] == 'correct').sum() / len(params))
      else:
        init = params[params['mode']=='INIT']
        current_behav = (init['correctness'] == 'correct').sum() / len(init)
        behavior.append((init['correctness'] == 'correct').sum() / len(init))

      all_results.append(final_results)
      folders_analyzed.append(current_folder)

      #per_lc = parameters['count_lc'] / (parameters['count_lc'] + parameters['count_lw']) ###
      #per_rc = parameters['count_rc'] / (parameters['count_rc'] + parameters['count_rw']) ###

      print("")
      print(current_folder)
      print(final_results)
      print("num_trials_tested",num_trials_tested)
      print("current_behav",current_behav)

      score_mean = final_results['score_mean']  # Assuming 'score_mean' is the column name in final_results
      score_means.append(score_mean)
      #print("percent left correct",per_lc) ###
      #print("percent right correct",per_rc) ###

    except:
     print(current_folder)
     print('unable to run analysis')
     pass

  # Print the score means without nested brackets
  print("############################################################################")
  print("SCORES FOR",select_cells)
  print( ", ".join(map(str, [mean[0] for mean in score_means])))




#%%
# FOR DECODING ON THE 1st OF THE PAIR

#CURRENTLY USING TF_PARAMS_PATH
stage = 's2' #s2 for init
#0 for turn decoding, 2 for visual stimuli decoding, 1 for correctness decoding
trial_type = 2
# 0 for ALL CELLS, 1 for matched cells
#select_cells = 0

#0 for no, 1 for yes
  #0 for no, 1 for yes
if stage in ['s2', 's5']:
    init_only = 1
else:
    init_only = 0


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

TF_params_path = r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params\sessions_for_cell_matching_init_e1_e2_240130.csv"
#TF_params_path = r"E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets\mice_params\sessions_for_cell_matching_init_e1_e2_extras.csv"

TF_params = pd.read_csv(TF_params_path, header=None) #HEADER NONE IS WHAT ALLOWS THIS TO GET ODD ROWS UNDER THE SAME FUNCTION
# Initialize a list to store sessionID values
selected_sessions = []
selected_sessions_expert = []
selected_sessions_naive = []
mice_accuracies_expert = []
mice_accuracies_naive = []

# for all sessions
# if index % 2 == 0 , decoding for e2
# if index % 2 != 0, decoding for e1

for index, row in TF_params.iterrows():
  # FOR ALL SESSIONS WITH EVEN INDICES ###############
    if index % 2 == 0:
    # Append the value in the 'sessionID' column to the list
      selected_sessions.append(row[2])




# Print the list of selected sessionIDs
#print(selected_sessions)
#########################################
temp = selected_sessions
# Filter out 'NaN' strings and NaN floats from temp
temp = [item for item in temp if (not isinstance(item, float) and item.lower() != 'nan') or (isinstance(item, float) and not math.isnan(item))]

#temp.sort()
# Print the sorted list
print(temp)

### IF USING DESKTOP THEN UNCOMMENT BELOW:
#temp = [path.replace('/content/drive/MyDrive', 'E:/google_drive') for path in temp]

select_cells_values = [0,1]

# Loop over select_cells values
for select_cells in select_cells_values:


  behavior = []
  folders_analyzed = []
  #temp.sort()
  param_paths = []
  cal_paths = []
  cal_paths_raw = []
  rejected_cells_paths = []
  all_results = []
  score_means = []

  for file in temp:
      #print(file)
      param_paths.append(glob.glob(file + '/trial_parameters' + '*.csv'))
      cal_paths.append(glob.glob(file + '/spikerate' + '*.csv'))
      #cal_paths_raw.append(glob.glob(file + '/calcium' + '*.csv'))
      rejected_cells_paths.append(glob.glob(file + '/rejected' + '*.csv'))

  for i in range(len(cal_paths)):
  #for i in range(4):
  #for i in range(3,4):
    #try:

      current_folder = temp[i]
      cal = pd.read_csv(cal_paths[i][0])
      #cal_raw = pd.read_csv(cal_paths_raw[i][0])
      params = pd.read_csv(param_paths[i][0])

      #if len(cal.columns) != len(cal_raw.columns):
      #  print('raw calcium and spike cell count mismatch, skipping')
      #  continue

      reject = pd.read_csv(rejected_cells_paths[i][0])
      reject = [item for sublist in reject.values.tolist() for item in sublist]


      del cal[cal.columns[0]]
      # #blood vessels remove
      #cal.drop(cal.columns[reject], axis=1, inplace=True)
      #reset column indices / cell ids
      #cal = cal.transpose()
      #cal = cal.reset_index(drop=True)
      #cal = cal.transpose()

      if select_cells == 0:
        cal.drop(cal.columns[reject], axis=1, inplace=True)
        #reset column indices / cell ids
        cal = cal.transpose()
        cal = cal.reset_index(drop=True)
        cal = cal.transpose()
      if select_cells == 1:
      # decide what cells to use in analysis
        c1_matched_cells_no_rejects = [num for num in c1_matched_cells_per_mouse[i] if num not in reject]
        cal = cal.iloc[:, c1_matched_cells_no_rejects]
      if select_cells == 2:
        c1_unmatched_cells_no_rejects = [num for num in c1_unmatched_cells_per_mouse[i] if num not in reject]
        cal = cal.iloc[:, c1_unmatched_cells_no_rejects]

      parameters = get_parameters(cal,params)
      spike_sums = get_spikerate_sums_per_cell_per_trial(cal,parameters)
      method = spike_sums

      num_incorrect_use = parameters['num_incorrect_trials'] #use max
      num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)
      if num_trials_tested < 10:
        num_incorrect_use = int(0.2 * parameters['ntrials']) #0.2
        num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)

      if num_trials_tested >= 10:
        try:
          svm_results = get_results(method)
          final_results = get_final_stats(svm_results)
        except:
        #  continue
          num_incorrect_use = int(0.1 * parameters['ntrials']) #0.2
          num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)
          svm_results = get_results(method)
          final_results = get_final_stats(svm_results)
      else:
        continue

      if init_only == 0:
        current_behav = (params['correctness'] == 'correct').sum() / len(params)
        behavior.append((params['correctness'] == 'correct').sum() / len(params))
      else:
        init = params[params['mode']=='INIT']
        current_behav = (init['correctness'] == 'correct').sum() / len(init)
        behavior.append((init['correctness'] == 'correct').sum() / len(init))

      all_results.append(final_results)
      folders_analyzed.append(current_folder)

      #per_lc = parameters['count_lc'] / (parameters['count_lc'] + parameters['count_lw']) ###
      #per_rc = parameters['count_rc'] / (parameters['count_rc'] + parameters['count_rw']) ###

      print("")
      print(current_folder)
      print(final_results)
      print("num_trials_tested",num_trials_tested)
      print("current_behav",current_behav)

      score_mean = final_results['score_mean']  # Assuming 'score_mean' is the column name in final_results
      score_means.append(score_mean)
      #print("percent left correct",per_lc) ###
      #print("percent right correct",per_rc) ###

    #except:
     #print(current_folder)
     #print('unable to run analysis')
     #pass

  # Print the score means without nested brackets
  print("############################################################################")
  print("SCORES FOR",select_cells)
  print( ", ".join(map(str, [mean[0] for mean in score_means])))




#%%
TF_params

#%%
c1_matched_cells_no_rejects

#%%
c1_matched_cells_per_mouse[1]

#%%
reject

#%%
c2_matched_cells_no_rejects

#%%
cal

#%%
matched_cells

#%%
for i in cal.columns:
  print(i)

#%%
reject

#%%
results_dict['m21r']['c2_matched_cells']

#%%
### FUNCTIONS
def get_parameters(cal,params):

    clean_df = params[params['clarity']=='clean']    ##### default
    a = 1   ##### default

    #if want only correct trials
    #clean_df = clean_df[clean_df['correctness']=='correct']

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
      if trial_type == 0:
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
      if trial_type == 2:
        if clean_df['direction'].iloc[i]=='L':
          clean_trial_direction.append(0)
        if clean_df['direction'].iloc[i]=='R':
          clean_trial_direction.append(1)

        ##### MAIN BELOW                                    5
    #ntrials = len(clean_trial_turn_frames)
    ntrials = len(clean_trial_end_frames)

    # use only for reward decoding, where you want to take start:stimstart but reward based on "previous" trial so it aligns
    #import math
    #def shift_values_up_with_nan(lst):
    #  lst.insert(0, math.nan)
    #  lst.pop()
    #  return lst
    #shift_values_up_with_nan(y_all)

    #def shift_values_up(lst):
    #  lst.insert(0, lst.pop())
    #  return lst
    #shift_values_up(y_all)

    #clean_trial_start_frames = clean_trial_start_frames.pop(0)
    #clean_trial_stim_start_frames = clean_trial_stim_start_frames.pop(0)
    #clean_trial_end_frames = clean_trial_end_frames.pop(0)
    #clean_trial_turn_frames = clean_trial_turn_frames.pop(0)
    #clean_trial_correctness = clean_trial_correctness.pop(0)
    #clean_trial_direction = clean_trial_direction.pop(0)

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

def get_event_window_lengths(cal,base_params):
  stimstart_turn = []
  turn_end = []
  end_stimstart = []

  for i in range(base_params['ntrials']):
    stimstart_turn.append(len(cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_turn_frames'][i]]))
    turn_end.append(len(cal.iloc[base_params['clean_trial_turn_frames'][i]:base_params['clean_trial_end_frames'][i]]))
    end_stimstart.append(len(cal.iloc[base_params['clean_trial_end_frames'][i]:base_params['clean_trial_stim_start_frames'][i]]))

  event_window_lengths = {
      'stimstart_turn':stimstart_turn,
      'turn_end':turn_end,
      'end_stimstart':end_stimstart
  }
  return event_window_lengths

def get_distances(base_params,df_dist):
  dist = []
  for t in range(parameters['ntrials']):
    a_frames = []
    for i in range(parameters['clean_trial_start_frames'][t] , parameters['clean_trial_turn_frames'][t]):
      if df_dist['label_dist_from_center_along_arm'][i] <= abs(150):
        a_frames.append(i)


    dist.append(a_frames)
  return dist

def get_spikerate_sums_per_cell_per_trial(cal,base_params):
    t_sums = []
    for i in range(base_params['ntrials']):
                #### DEL LATER , MESSY AND CLEAN TRIAL TESTING          #8
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_end_frames'][i]])
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_stim_start_frames'][i]+25])
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_turn_frames'][i]])

                tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-15:base_params['clean_trial_turn_frames'][i]+25])
                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-14:base_params['clean_trial_turn_frames'][i]]) #preturn hard defined frames
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_turn_frames'][i]]) #preturn
                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]:base_params['clean_trial_end_frames'][i]-3])  # MAIN

                #tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]-30:base_params['clean_trial_end_frames'][i]-3])


                #tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]:base_params['clean_trial_end_frames'][i]+60]) #reward

                #for end:stimstart
                #try:
                #  tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]:base_params['clean_trial_stim_start_frames'][i+1]])
                #  #print('working')
                #except:
                #  tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]:base_params['clean_trial_end_frames'][i]+30])
                #  #print('exception')

                #tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]:base_params['clean_trial_end_frames'][i]+60])
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_turn_frames'][i]]) ###MAIN

                ################## RUN PCA FIRST
                #import numpy as np
                #from sklearn.decomposition import PCA
                ## Assuming your 2D array is named 'data'
                ##tframes = np.transpose(tframes)
                #data = tframes
                ## Create an instance of PCA with the desired number of components
                #pca = PCA(n_components=2)
                ## Fit the PCA model to the data and transform it to the new coordinate system
                #tframes = pca.fit_transform(data)
                ##################

                tframes = np.asarray(tframes)
                #nframes = np.shape(tframes)[0]
                t_sum = np.sum(tframes,axis=0)
                #t_sum_avg = t_sum / nframes
                t_sums.append(t_sum)
    return t_sums

def get_train_test_split_matching_num_incorrect_hardset(x_all,parameters,num_incorrect_use):
    y_all = parameters['y_all']
    wrong_idx = []
    correct_idx = []
    for i in range(len(y_all)):
        if y_all[i] == 0:
            wrong_idx.append(i)
        else:
            correct_idx.append(i)

    x_train_idx_wrong = np.random.choice(len(wrong_idx), num_incorrect_use, replace=False)
    x_train_wrong = []
    for i in x_train_idx_wrong:
        x_train_wrong.append(wrong_idx[i])
    x_test_wrong = list(set(wrong_idx).symmetric_difference(set(x_train_wrong)))

    x_train_idx_correct = np.random.choice(len(correct_idx), num_incorrect_use, replace=False)
    x_train_correct = []
    for i in x_train_idx_correct:
        x_train_correct.append(correct_idx[i])
    x_test_correct = list(set(correct_idx).symmetric_difference(set(x_train_correct)))

    # make y_train
    y_train = []
    for i in range(num_incorrect_use):
        y_train.append(0)
    for i in range(num_incorrect_use):
        y_train.append(1)

    # make y_test
    y_test = []
    for i in range(len(x_test_wrong)):
        y_test.append(0)
    for i in range(len(x_test_correct)):
        y_test.append(1)

    # make x_train
    x_train = []
    x_train_idx = x_train_wrong + x_train_correct
    for i in x_train_idx:
        x_train.append(x_all[i])

    # make x_test
    x_test = []
    x_test_idx = x_test_wrong + x_test_correct
    for i in x_test_idx:
      try:
        x_test.append(x_all[i])
      except:
        pass

    train_test_groups = {
        'x_train':x_train,
        'y_train':y_train,
        'x_test':x_test,
        'y_test':y_test
        }
    return train_test_groups

def run_svm():
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import confusion_matrix, accuracy_score
    from sklearn.svm import SVC
    import random

    from sklearn.model_selection import cross_val_score

    train_test_groups = get_train_test_split_matching_num_incorrect_hardset(x_all,parameters,num_incorrect_use)

    x_train = train_test_groups['x_train']
    x_test = train_test_groups['x_test']
    y_train = train_test_groups['y_train']
    y_test = train_test_groups['y_test']

    ############### NEW USING PCA FIRST
    #from sklearn.decomposition import PCA
    #pca = PCA(n_components=3)
    #x_train_pca = pca.fit_transform(x_train)
    #x_test_pca = pca.transform(x_test)
    #classifier = SVC()
    #classifier.fit(x_train_pca, y_train)
    #y_pred = classifier.predict(x_test_pca)
    #score = accuracy_score(y_test, y_pred)
    ###############

    ##### KEEP DEFAULT
    classifier = SVC()
    classifier.fit(x_train, y_train)
    y_pred = classifier.predict(x_test)
    score = accuracy_score(y_test,y_pred)
    return score

def run_svm_shuffled():
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import confusion_matrix, accuracy_score
    from sklearn.svm import SVC
    import random

    from sklearn.model_selection import cross_val_score

    train_test_groups = get_train_test_split_matching_num_incorrect_hardset(x_all,parameters,num_incorrect_use)

    x_train = train_test_groups['x_train']
    x_test = train_test_groups['x_test']
    y_train = train_test_groups['y_train']
    y_test = train_test_groups['y_test']

    y_s = y_train
    y_s = random.sample(y_s, len(y_s))
    y_s = np.asarray(y_s)

    ##### KEEP
    classifier = SVC()
    classifier.fit(x_train, y_s)
    y_pred = classifier.predict(x_test)
    score = accuracy_score(y_test,y_pred)
    return score

def get_results(t_sums):
    window = 1
    for i in range(window):
        global x_all
        x_all = t_sums[i::window]

        svm_results_experimental = []
        svm_results_shuffled = []
        for _ in range(1000):
            svm_results_experimental.append(run_svm())
            svm_results_shuffled.append(run_svm_shuffled())

        svm_results = {
            'svm_results_experimental':svm_results_experimental,
            'svm_results_shuffled':svm_results_shuffled,
            }

    return svm_results

def get_final_stats(svm_results):
    import statistics
    results = {"score_mean":[],"score_std":[],"shuffled_score_mean":[],"shuffled_score_std":[]}

    results['score_mean'].append(statistics.mean(svm_results['svm_results_experimental']))
    results['score_std'].append(statistics.stdev(svm_results['svm_results_experimental']))

    results['shuffled_score_mean'].append(statistics.mean(svm_results['svm_results_shuffled']))
    results['shuffled_score_std'].append(statistics.stdev(svm_results['svm_results_shuffled']))

    return results

