#%%
#SPIKE RATE SVM FOR MICE MAIN PRISM
basepath = r'E:\main_data\Calcium_imaging_yarm_LurLab\mice_datasets' #hd
#basepath = r'\\192.168.1.9\dendrite\Ali\mice_datasets' #dendrdite

mice = ['m2r','m9r','m11n','m12lr','m13l','m14r','m15n','m16l','m17r','m19n','m21r','m22l','m23r','m25lr','m26r','m29n', 'm30lr', 'm31r','m32l','m33r','m34n','m36rr','m37l','m38n','m39rr','m40l', 'm41r','min4r','sn','snr'] #main mice except #m27l

stage = 's1'

mice_accuracies_expert = []
mice_accuracies_naive = []

mice_shuffled_accuracies_expert = []
mice_shuffled_accuracies_naive = []

#mice_std_expert = []
#mice_std_naive = []

for m in mice:

  #0 for no, 1 for yes
  if stage in ['s2', 's5']:
      init_only = 1
  else:
      init_only = 0


  #DEL LATER
  #init_only = 0

  #0 for turn decoding, 2 for visual stimuli decoding, 1 for correctness decoding
  trial_type = 2

  import os
  import glob
  import pandas as pd
  import math
  import numpy as np
  import re

  def getsubfolders(folder):
      ''' returns list of subfolders '''
      return [os.path.join(folder,p) for p in os.listdir(folder) if os.path.isdir(os.path.join(folder,p))]

  basefolders = []
  temp = []
  subfolders=getsubfolders(basepath)

  for subfolder in subfolders:
    subfolders=getsubfolders(subfolder)
    #if m in subfolder:
    pattern = r'\b{}\b'.format(re.escape(m))
    if re.search(pattern, subfolder):
      for subfolder in subfolders:
        subfolders=getsubfolders(subfolder)
        pattern2 = r'\b{}\b'.format(re.escape(stage))
        if re.search(pattern2, subfolder):
          for subfolder in subfolders:
            subfolders=getsubfolders(subfolder)
            temp.append(subfolder)

  behavior = []
  all_results = []
  folders_analyzed = []

  master_n_clean_trials = []
  master_n_incorrect_trials = []
  all_folders = []

  temp.sort()
  param_paths = []
  cal_paths = []
  rejected_cells_paths = []
  for file in temp:
      param_paths.append(glob.glob(file + '/trial_parameters' + '*.csv'))
      cal_paths.append(glob.glob(file + '/spikerate' + '*.csv'))
      #cal_paths.append(glob.glob(file + '/calcium' + '*.csv'))
      rejected_cells_paths.append(glob.glob(file + '/rejected' + '*.csv'))
      #param_paths.append(glob.glob(file + '/with_real_direction' + '*.csv')) ### FOR HABIT ONLY

  for i in range(len(cal_paths)):
  #for i in range(4):
  #for i in range(2,3):
    try:

      current_folder = temp[i]
      cal = pd.read_csv(cal_paths[i][0])
      params = pd.read_csv(param_paths[i][0])
      reject = pd.read_csv(rejected_cells_paths[i][0])
      reject = [item for sublist in reject.values.tolist() for item in sublist]


      del cal[cal.columns[0]]
      #blood vessels remove
      cal.drop(cal.columns[reject], axis=1, inplace=True)
      #reset column indices / cell ids
      cal = cal.transpose()
      cal = cal.reset_index(drop=True)
      cal = cal.transpose()


      parameters = get_parameters(cal,params)
      spike_sums = get_spikerate_sums_per_cell_per_trial(cal,parameters)
      method = spike_sums

      num_incorrect_use = parameters['num_incorrect_trials'] #use max
      num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)
      if num_trials_tested < 10:
        num_incorrect_use = int(0.2 * parameters['ntrials']) #0.2
        num_trials_tested = parameters['ntrials'] - (num_incorrect_use*2)

      if num_trials_tested >= 10:
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

      #print("")
      print(current_folder)
      print(final_results)
      print("num_trials_tested",num_trials_tested)
      print("current_behav",current_behav)
      #print("percent left correct",per_lc) ###
      #print("percent right correct",per_rc) ###

    except:
      pass

  mouse_accuracies = [d['score_mean'] for d in all_results]
  mouse_accuracies = [item[0] for item in mouse_accuracies]

  mouse_std = [d['score_std'] for d in all_results]
  mouse_std = [item[0] for item in mouse_std]

  mouse_accuracies_shuffled = [d['shuffled_score_mean'] for d in all_results]
  mouse_accuracies_shuffled = [item[0] for item in mouse_accuracies_shuffled]

  naive_accuracies = []
  expert_accuracies = []
  naive_accuracies_shuffled = []
  expert_accuracies_shuffled = []

  #naive_std = []
  #expert_std = []
  for i in range(len(mouse_accuracies)):
    #if behavior[i]*100 < 75:
    if 40 < behavior[i]*100 < 60:
    #if behavior[i]*100 <100:         #TEST DEL
      naive_accuracies.append(mouse_accuracies[i])
      naive_accuracies_shuffled.append(mouse_accuracies_shuffled[i])
      #naive_std.append(mouse_std[i])
    #else:
    if behavior[i]*100 > 65:
    #if behavior[i]*100 > 75:
    #if behavior[i]*100 > 100: #TEST DEL
      expert_accuracies.append(mouse_accuracies[i])
      expert_accuracies_shuffled.append(mouse_accuracies_shuffled[i])

      #expert_std.append(mouse_std[i])

  #naive = [item[0] for item in naive]
  #expert = [item[0] for item in expert]
  mouse_accuracy_avg_expert = np.mean(expert_accuracies)
  mouse_accuracy_avg_naive = np.mean(naive_accuracies)

  mouse_shuffled_accuracy_avg_expert = np.mean(expert_accuracies_shuffled)
  mouse_shuffled_accuracy_avg_naive = np.mean(naive_accuracies_shuffled)

  #mouse_std_avg_expert = np.mean(expert_std)
  #mouse_std_avg_naive = np.mean(naive_std)

  #print("expert average accuracy",mouse_accuracy_avg_expert)
  #print("naive average accuracy",mouse_accuracy_avg_naive)

  mice_accuracies_expert.append(mouse_accuracy_avg_expert)
  mice_accuracies_naive.append(mouse_accuracy_avg_naive)

  mice_shuffled_accuracies_expert.append(mouse_shuffled_accuracy_avg_expert)
  mice_shuffled_accuracies_naive.append(mouse_shuffled_accuracy_avg_naive)

  #mice_std_expert.append(mouse_std_avg_expert)
  #mice_std_naive.append(mouse_std_avg_naive)

for i, m in enumerate(mice):
  print(f"Mouse {m}:")
  print(f"  Naive Accuracy: {mice_accuracies_naive[i]}")
  print(f"  Expert Accuracy: {mice_accuracies_expert[i]}")
  print("")

# Print lists of all average naive and expert accuracies
print("All Average Naive Accuracies:")
print(", ".join(map(str, mice_accuracies_naive)))
print("")

# Print All Average Naive Accuracies Shuffled without brackets
print("All Average Naive Accuracies Shuffled:")
print(", ".join(map(str, mice_shuffled_accuracies_naive)))
print("")

# Print All Average Expert Accuracies without brackets
print("All Average Expert Accuracies:")
print(", ".join(map(str, mice_accuracies_expert)))
print("")

# Print All Average Expert Accuracies Shuffled without brackets
print("All Average Expert Accuracies Shuffled:")
print(", ".join(map(str, mice_shuffled_accuracies_expert)))
print("")

print("Final average naive:",np.nanmean(mice_accuracies_naive))
print("Final average expert:",np.nanmean(mice_accuracies_expert))

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

    ## use for turn decoding HAIBT only
      #if clean_df['real_direction'].iloc[i]=='L':
      #  count_lc.append(1)
      #  clean_trial_direction.append(0)
      #elif clean_df['real_direction'].iloc[i]=='R':
      #  count_rc.append(1)
      #  clean_trial_direction.append(1)


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

        #del later     4
        #ntrials = len(clean_trial_start_frames)


    ##### MAIN BELOW                                    5
    #ntrials = len(clean_trial_turn_frames)
    ntrials = len(clean_trial_end_frames)

    #counter_removed_trials = 0
    #for i in range(len(clean_trial_turn_frames)):
    #    if clean_trial_turn_frames[i] > len(cal):
    #        counter_removed_trials = +1
    #        del clean_trial_turn_frames[i]
    #        del clean_trial_stim_start_frames[i]
    #ntrials = ntrials - counter_removed_trials

    y_all = []
    if trial_type == 1:
        for i in range(ntrials):
            y_all.append(clean_trial_correctness[i])
    else:
        for i in range(ntrials):
            y_all.append(clean_trial_direction[i])

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

                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-21:base_params['clean_trial_turn_frames'][i]-15])
                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-14:base_params['clean_trial_turn_frames'][i]]) #preturn hard defined frames
                #tframes = (cal.iloc[base_params['clean_trial_stim_start_frames'][i]:base_params['clean_trial_turn_frames'][i]]) #preturn
                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]:base_params['clean_trial_end_frames'][i]-3])

                #tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]-30:base_params['clean_trial_end_frames'][i]-3])
                tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]-15:base_params['clean_trial_turn_frames'][i]+25]) #MAIN

                #tframes = (cal.iloc[base_params['clean_trial_turn_frames'][i]:base_params['clean_trial_turn_frames'][i]+6])


                #tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]:base_params['clean_trial_end_frames'][i]+60]) #reward

                #tframes = (cal.iloc[base_params['clean_trial_end_frames'][i]:base_params['clean_trial_end_frames'][i]+40]) #reward2

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

