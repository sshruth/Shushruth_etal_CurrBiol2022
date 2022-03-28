function E = eventCodes

E.monkeyID = 1; 
% Charlie 2097 Annie 7208
% Mazurek data: 0 = natasha; 1 = bridget 

E.unitID = 2; 
% DATE.DEPTH i.e. yyyymmdd.xxxxx micrometers
% DATE.TIME for array recordings. Time is the start recording time of file

E.fileID = 3; 
% Trial ID: (B)lock#  and (T)rial # as BBTTT

E.object_in_RF = 4;
% Mazurek data: 0=none,1=T1,2=T2,3=dots
% Annie abstract dots: Expo target ID

E.trial_type = 5;
% Info in plxCodes 

E.dot_diam = 6;
E.coherence = 7;
E.dot_dir = 8;
E.target1_x = 9;
E.target1_y = 10;
E.target2_x = 11;
E.target2_y = 12;
E.dot_duration = 13;
% Also the computer chosen go-time on saccade tasks

E.fixation_x = 14;
E.fixation_y = 15;
E.dots_x = 16;
E.dots_y = 17;
E.dot_speed = 18;

E.pct_novar = 19;
E.seedvar = 20;
E.seed = 21;
E.flicker = 22;

E.time_FP_on = 23;
E.time_fix_acq = 24;
E.time_target_on = 25;
E.time_target_off = 26;
E.time_dots_on = 27;
E.time_dots_off = 28;
E.time_FP_off = 29;
E.time_saccade = 30;
E.time_targ_acq = 31;
E.time_reward = 32;
E.time_end = 33;

E.target_choice = 34;
E.correct_target = 35;
% In saccade tasks, ID of the target in the trial

E.isCorrect = 36;
E.react_time = 37;
E.eye_trl_class = 38;
E.eye_t_sac = 39;
E.eye_rt = 40;
E.RF_x = 41;
E.RF_y = 42;

% Target colors
E.target1_r = 43;
E.target1_g = 44;
E.target1_b = 45;
E.target2_r = 46;
E.target2_g = 47;
E.target2_b = 48;
