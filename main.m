%% Main File Mimimum Work Example
%==== For using optimization-related modules, please use the following two lines to set up the directories of yalmip and Gurobi.
% Yalmip: addpath('./YALMIP-master'); 
% Gurobi PATH for mac: cd('/Library/gurobi801/mac64/matlab'); gurobi_setup; savepath

%==== Clear All Data
function main(app,PATH_dir,figoption)
% close all
clc
clearvars -except app PATH_dir figoption

diary myDiaryFile
%% ==== STEP1: Read in the hazard scenario
% hazard: type 1-earthquake, 2-hurricane 
% Load Input File
run('Data_Supplement_Real.m');
tu = linspace(1,time_horizon,time_horizon);

%--- Set up random samples following uniform distributions with the Latin Hybercube Sampling (LHS) method
[SamplesUniform, SamplesUnifDmg, SamplesUnifDur] = Library.SampleInititalUniformSample(Nsamples, NRun, Total, Seed, active_truerands);
msg = '----STEP1: Complete Reading Input Data----';
disp(msg);

%---- Terminate the program when clicking the "Stop" button
stopVariable = get(app.StopButton, 'Enable');
Library.AbortComputation(stopVariable);

%---- Add the progress bar with the abort button
Library.gui_active(1); % will add an abort button
increment = 0.1;
string = 'Please wait';
titlestr = 'Progress';
myhandle = Library.progressbar( [],increment,string,titlestr );


%% ==== STEP2: assess initial damage(s) under the scenario and determine restoration tasks
[GraphSystem,TotalDamage,RestorationTable,Nsystask] = Library.SampleDamage(EventType, Nsamples, NRun, Hostname, IM, IMx, IMy, ActiveSystem, powerGraph, commGraph, transGraph, Prob_Magic_Battery, SamplesUnifDmg, Seperate_Scheduling, InterdependencePrec);

%---- Save samples of damage scenario and task data
%DamageScenario = TotalDamage;
%save(strcat( deblank(Hostname), '/mat/DamageScenario.mat'), 'DamageScenario', 'GraphSystem');
save(strcat( deblank(Hostname), '/mat/TaskTable.mat'), 'RestorationTable');


%---- Save damage scenarios and the corresponding restoration tasks
msg = '----STEP2: Complete Sampling Damage Scenarios----';
disp(msg);

%---- Terminate the program when clicking the "Stop" button
stopVariable = get(app.StopButton, 'Enable');
Library.AbortComputation(stopVariable);

%---- Add the progress bar
increment = 0.1;
string = '20%';
titlestr = 'Progress';
myhandle = Library.progressbar( myhandle,increment,string,titlestr );


%% ==== STEP3: Planing a restoration plan as a restoration sequence schedule
Nsample = Nsamples*NRun;
progress = 0.2; % 20% = 10% (1st step) + 10% (2nd step);
for isample = 1:Nsample
    SystemSetDamage = TotalDamage{isample,1};
    TableTask = RestorationTable{isample,1};
    TablePrecedence = RestorationTable{isample,2};
    [Schedule, Date] = Library.RestorationPlan(Scheduler_Num, time_horizon, RestorationResource, ...
                active_power, active_comm, active_trans, SystemSetDamage, TableTask, TablePrecedence,... 
                priority_power_num, priority_transportation_num, priority_communication_num, OptimizationChoice);
    TotalSchedule{isample,1} = Schedule;   
    TotalDate{isample,1} = Date; 
    %--- save schedule data in the mat files 
    fnoutput = strcat(deblank(Hostname), '/schedule/Data_Schedule_', num2str(isample)); 
    save(fnoutput, 'Schedule', 'Date'); 
    
    %--- Terminate the program when clicking the "Stop" button
    stopVariable = get(app.StopButton, 'Enable');
    Library.AbortComputation(stopVariable);

    %---- Add the progress bar
    increment = inv(Nsample)*0.3;
    progress = progress + increment; 
    string = strcat(num2str(100*progress), '%');
    titlestr = 'Progress';
    myhandle = Library.progressbar( myhandle,increment,string,titlestr );
    
end
%---- Save the restoration schedule for every damage scenario sample as MAT and TXT files
save(strcat(deblank(Hostname), '/mat/Schedule.mat'), 'TotalSchedule', 'TotalDate');
Library.SaveLogSchedule(Nsample, TotalSchedule, TotalDate, Hostname);

msg = '----STEP3: Complete Developing Restoration Plans----';
disp(msg);




%% ==== STEP4: Simulation the Recovery Process 
%%%% 1. Update Object.Status
%%%% 2. Compute System Functionality
%---- sample the duration sample for every restoration task by mapping
%from uniform samples to new duration samples following the
%corresponding distribution features of task duration 
%---- progress = 0.5; % 50% = 10%+10%+30%
for isample = 1:Nsample

    [TotalDamage] = Library.SampleActualDuration(TotalDamage, SamplesUnifDur, isample);
    
    %--- Terminate the program when clicking the "Stop" button
    stopVariable = get(app.StopButton, 'Enable');
    Library.AbortComputation(stopVariable);

    %---- Add the progress bar
    increment = inv(Nsample)*0.1;
    progress = progress + increment; 
    string = strcat(num2str(100*progress), '%');
    titlestr = 'Progress';
    myhandle = Library.progressbar( myhandle,increment,string,titlestr );

end

%---- Simulate the functionality recovery at system-level by implementing
%different task duration sample in every functionality recovery sample: 
% isample: index of random sample
% kdmg: index of random sample for damage scenario
% jdur: index of random sample for task duration

for isample = 1:Nsample 
    
    kdmg = ceil(isample/NRun);
    jdur = isample - (kdmg-1)*NRun;

    % Set up tables of task and precedence for every task duration sample
    % Set up the repair schedule for every damage scenario sample
    TableTask = RestorationTable{isample,1};
    TablePrecedence = RestorationTable{isample,2}; 
    RepairSchedule = TotalSchedule{isample,1};
    
    % Set up the system set and initial network graph for every system in
    % every task duration sample 
    Total = TotalDamage{isample};
    Graph = GraphSystem{isample};
    
    % At every jumping time step (when a task is completed), update Object.Status and simulate the functionality recovery of every system
    [QSys{isample}, QNbr, trackt, trackCW, trackOS, lookupTable] = Library.Recovery(ActiveSystem, time_horizon, Interdependence_Num, Qtrans0, InterdependenceFunc, ReSchedule_Num, ...
        RestorationResource, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, ...
        RepairSchedule, Total, Graph, System_Dependent_Factor, ...
        Seperate_Scheduling, LinkDirection,Per_New_reschedule,Num_stage,...
        TableTask, TablePrecedence, OptimizationChoice, Diff_unit, Cust_unit);
    
    %---- Document the following intermediate data
    %1. The content of Current Working List at different time steps in the recovery process 
    %2. Object Status update at different time steps in the recovery process 
    %3. Jumping time steps in the recovery process 
    zCurrentWorking{isample} = trackCW;
    zObjectStatus{isample} = trackOS;
    zt{isample} = trackt(:,2); 
    zTaskLookupTable{isample} = lookupTable;
    
    for ii = 1:size(QNbr,1)
        QuNbr{isample}(ii,:) = interp1(zt{isample}, QNbr(ii,:), tu);
    end
    
    %--- Terminate the program when clicking the "Stop" button
    stopVariable = get(app.StopButton, 'Enable');
    Library.AbortComputation(stopVariable);

    %---- Add the progress bar
    increment = inv(Nsample)*0.2;
    progress = progress + increment; 
    string = strcat(num2str(100*progress), '%');
    titlestr = 'Progress';
    myhandle = Library.progressbar( myhandle,increment,string,titlestr );
    
end
    
%---- Interpolate the functionality recovery data on NONUNIFORM time steps
% to compute the functionality recovery data on UNIFORM time steps
% tu = the time vector uniformly starting from t = 1 day to t = time_horizon
% idx1 = the index of element as ones (100%) in the QSys
% n1 = the number of element equal to 1
% isample: index of random sample
% kdmg: index of random sample for damage scenario
% jdur: index of random sample for task duration

for isample = 1:Nsample    
    
    kdmg = ceil(isample/NRun);
    jdur = isample - (kdmg-1)*NRun;

    for isys = 1:length(QSys{isample})
        tmpQ = QSys{isample}{isys};
        nq = 10; % total number of possible system functionality metrics
        for imtr = 1:nq
            tmpq = tmpQ{imtr};
            idx1 = find(tmpq == 1); 
            n1 = length(idx1);
            if gt(n1, 2) % Qsys recovers to 100% at the end in multiple time steps (including Qsys(t=1)=100%).
                QuSys{isys}{imtr}(isample,:) = interp1(zt{isample}, tmpq, tu);
            elseif n1 == 1 % Qsys does not recovers to 100% at the end (t = time_horizon).
                QuSys{isys}{imtr}(isample,:) = interp1(zt{isample}, tmpq, tu);
            elseif n1 == 2 % Qsys recovers to 100% at the last time step (including Qsys(t=1)=100%, and Qsys(t=time_horizon)=100%).
%                 if Diff_unit
%                 azt = zt{isample};
%                 aztp = azt(1:end);
%                 axtpmax = round(max(aztp));
%                 atu = linspace(1, axtpmax, axtpmax);
%                 tmpQuSys = interp1(aztp, tmpq(1:end), atu);
%                 else
                azt = zt{isample};
                aztp = azt(1:end-1);
                axtpmax = round(max(aztp));
                atu = linspace(1, axtpmax, axtpmax);
                tmpQuSys = interp1(aztp, tmpq(1:end-1), atu);
%                 end
                check = find(isnan(tmpQuSys));
                if ~isempty(check)
                    tmpQuSys(end) = tmpq(end);
                end
                
                %a = [isample, isys, imtr]
                
                check2 = time_horizon - axtpmax;
                if gt(check2, 0)
                    QuSys{isys}{imtr}(isample,:) = [tmpQuSys, ones(1, time_horizon - axtpmax)];
                else
                    QuSys{isys}{imtr}(isample,:) = tmpQuSys(1:time_horizon);
                end
                       
                
            end
        end
            
        
    end

    msg = strcat('----Finished Sample No. ', num2str(isample), ' [meaning DamageSample #', num2str(kdmg), ' and DurationSample # ', num2str(jdur), ']-----');
    disp(msg);
     
    %--- Terminate the program when clicking the "Stop" button
    stopVariable = get(app.StopButton, 'Enable');
    Library.AbortComputation(stopVariable);

    %---- Add the progress bar
    increment = inv(Nsample)*0.1;
    progress = progress + increment; 
    string = strcat(num2str(100*progress), '%');
    titlestr = 'Progress';
    myhandle = Library.progressbar( myhandle,increment,string,titlestr );
end   

msg = '----STEP4: Complete the Recovery Process----';
disp(msg);



%% STEP 5 Resilience Assessment and Data Post-Processing

[Resilience] = Library.ComputeSystemResilience(ActiveSystem, nQmetric, ResilienceMetricChoice, time_horizon, QuSys);
[StatisticsFunctionality, StatisticsResilience] = Library.ComputeStatistics(ActiveSystem, nQmetric, QuSys, Resilience);
msg = '----STEP5: Complete the Resilience Assessment----';
disp(msg);

diary off

%% Output data and figure
%--- Export Data Files into MAT files
save(strcat(deblank(Hostname), '/mat/Functionality.mat'), 'QuSys', 'QuNbr','StatisticsFunctionality');
save(strcat(deblank(Hostname), '/mat/Resilience.mat'), 'Resilience', 'StatisticsResilience');
save(strcat(deblank(Hostname), '/mat/OtherIntermediateData.mat'), 'zt', 'zCurrentWorking', 'zObjectStatus', 'zTaskLookupTable', 'QSys');

%--- Save functionality into TXT files
idQmtr{1} = [1]; % functionality metric index for the power system
idQmtr{2} = [1]; % functionality metric index for the communication system
idQmtr{3} = [1]; % functionality metric index for the transportation system
Library.SaveLogFunctionality(QuSys, Hostname, Nsample, NRun, idQmtr);

%--- Plot Figures (this section is not applicable on HPC.)
if ~exist('figoption','var')
    figoption = 'on';
end
Library.PlotFigureFunctionality(Nsamples, NRun, time_horizon, QuSys, Hostname, figoption);
 
%--- Profile for Parallel Computing 
if Profile_Num == 1
    spmd
        mpiprofile('viewer');
        mpiprofile('off');
    end
    profile off
    profsave(profile('info'),strcat(deblank(Hostname), '/profile_log'));
end


%---- Add the progress bar
increment = 0.1;
string = '100%';
titlestr = 'Progress';
myhandle = Library.progressbar( myhandle,increment,string,titlestr );
pause(0.2);
Library.progressbar( myhandle,-1 ); % close the progress bar


%% Clean Data
clear ans index i n Comm_Fun data_index Days m name finish filename Pow_Fun Trans_Fun Total TotalSchedule Transportation_Set Power_Set Communication_Set used used_new ava ava_new;

end

