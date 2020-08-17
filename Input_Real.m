%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input_Real.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==== Set up initial state of flag variables to indicate system/neighborhood being active (true) or not (false)
active_power = false;
active_comm = false;
active_trans = false;
active_neighbor = false;
LinkDirectionChoice = 0;

%==== Set up the system set as empty cells
Power = {{}, {}, {}, {},{}};
Comm = {{}, {}, {}, {}};
Trans = {{}, {}, {}, {}, {}};

%==== Set up the initial resource constraint vector for every
%system/neighborhood
Power_Resource = zeros(1,4);
Comm_Resource = zeros(1,4);
Trans_Resource = zeros(1,4);
Neighborhood = {};

%==== Set up the initial default value for rescheduling (turn OFF)
ReSchedule_Num = 0;

%% set up the Map/Dictionary 
% The following parameters will be called in the function of
% Library.readInput
Dictionary = containers.Map('KeyType','char','ValueType','any');

pow_check = containers.Map('KeyType','char','ValueType','int32'); 
centraloffice_check = containers.Map('KeyType','char','ValueType','int32');
communicationtower_check = containers.Map('KeyType','char','ValueType','int32');
bridge_check = containers.Map('KeyType','char','ValueType','int32');
trafficlight_check = containers.Map('KeyType','char','ValueType','int32');
road_check = containers.Map('KeyType','char','ValueType','int32');

% comm_check = {[], centraloffice_check, communicationtower_check};  
comm_check = {centraloffice_check, communicationtower_check};   
trans_check = {road_check, bridge_check, trafficlight_check};


%% INPUT Data
% Set up values of different variables from the input file
fid = fopen('Input.txt', 'r');
line = fgetl(fid);

while ischar(line)
    try
        tmp = strsplit(line);
        if ~strcmp(tmp(1), '#') 
            
            %=== set up the hazard scenario 
            % (1)Type: 1-earthquake, 2-hurricane wind, 3-flood
            % (2)MapData: IM-intensity measurement, IMx-latitude, IMy-longitude
            if strcmp(tmp(1), 'HazardMap')
                if length(tmp) < 3 
                    fprintf('Input File Error: Map - Missing Input of Hazard Information.\n');
                    return
                else %length(tmp)>=3
                    if strcmp(tmp(2), '1') 
                        EventType = 1; % Earthquake
                        msg = 'Input Setting: Set up an Earthquake hazard scenario.\n';
                        fprintf(msg);
                    elseif strcmp(tmp(2), '2')     
                        EventType = 2; % Hurricane wind
                        msg = 'Input Setting: Set up a Hurricane hazard scenario.\n';
                        fprintf(msg);
                    elseif strcmp(tmp(2), '3')     
                        EventType = 3; % Flood   
                        msg = 'Input Setting: Set up a Flood hazard scenario.\n';
                        fprintf(msg);
                    else
                        fprintf('Input File Error: Map - This type of hazard is not considered yet.\n');
                    end
                    % load the hazard data 
                    fpath = strcat('InputHazard/',char(tmp(3))); 
                    load(fpath); 
                end
            end
            
            %=== set up the Power system
            % (1) turn on(1)/off(0) the infrastructure system
            % (2) the remaining key words are the file names of input data
            % for this infrastructure system.    
            if strcmp(tmp(1), 'Power')
                if length(tmp) == 1
                    fprintf('Input File Error: Power - Missing Active System Input\n');
                    return
                elseif length(tmp)>1
                    system = {};
                    if strcmp(tmp(2), '1')
                        active_power = true; 
                        for ii = 3:length(tmp)
                            [Power, Trans, Comm, Neighborhood] = Library.readInput(char(tmp(ii)), pow_check, comm_check, trans_check, Power, Trans, Comm, Dictionary, Neighborhood);
                        end
                    elseif strcmp(tmp(2), '0')
                        active_power = false;
                        Branch_Set = {};
                        Bus_Set = {};
                        Generator_Set = {};
                        TransmissionTower_Set = {};
                        %Power = {Branch_Set, Bus_Set, Generator_Set, TransmissionTower_Set, {}}; % the last empty cell represents neighborhood
                        
                    else
                        fprintf('Input File Error: Power - Input not valid (0/1)\n');
                        return
                    end           
                   
                end
            end       
             
            %=== set up Communication system
            % (1) turn on(1)/off(1) the infrastructure system
            % (2) the remaining key words are the file names of input data
            % for this infrastructure system.    
            if strcmp(tmp(1), 'Communication')
                if length(tmp) == 1
                    fprintf('Input File Error: Communication - Missing Active System Input\n');
                    return
                else
                    if strcmp(tmp(2), '1')
                        active_comm = true; 
                        for ii = 3:length(tmp)
                            [Power, Trans, Comm, Neighborhood] = Library.readInput(char(tmp(ii)), pow_check, comm_check, trans_check, Power, Trans, Comm, Dictionary, Neighborhood);
                        end
                    elseif strcmp(tmp(2), '0')
                        active_comm = false;
                        Celline_Set = {};
                        Centraloffice_Set = {};
                        CommunicationTower_Set = {};
                        %Comm = {Centraloffice_Set, CommunicationTower_Set, Celline_Set, {}};  % the last empty cell represents neighborhood
                    else
                        fprintf('Input File Error: Communication - Input not valid (0/1)\n');
                        return                   
                    end    
                end
            end        
            
            %=== set up Transporation system
            % (1) turn on(1)/off(1) the infrastructure system
            % (2) the remaining key words are the file names of input data
            % for this infrastructure system.            
            if strcmp(tmp(1), 'Transportation')
                if length(tmp) == 1
                    fprintf('Input File Error: Transportation - Missing Active System Input\n');
                    return
                else
                    if strcmp(tmp(2), '1')
                        active_trans = true;
                        for ii = 3:length(tmp)-1
                            [Power, Trans, Comm, Neighborhood] = Library.readInput(char(tmp(ii)), pow_check, comm_check, trans_check, Power, Trans, Comm, Dictionary, Neighborhood,active_traLights);
                        end

%                         if exist(['./Input/',tmp{end},'.csv'])
%                             active_traLights = true;                
%                         else
%                             active_traLights = false; 
%                         end
                        [Power, Trans, Comm, Neighborhood] = Library.readInput(char(tmp(end)), pow_check, comm_check, trans_check, Power, Trans, Comm, Dictionary, Neighborhood,active_traLights);
                    elseif strcmp(tmp(2), '0')
                        active_trans = false;
                        Road_Set = {};
                        Bridge_Set = {};
                        TrafficLight_Set = {};
                        %Trans = {Road_Set, Bridge_Set, TrafficLight_Set, {}, {}}; % the last empty cell represents neighborhood    
                        
                    else
                        fprintf('Input File Error: Transportation - Input not valid (0/1)\n');
                        return                   
                    end    
                end
            end 
            
            %=== Set up Sub-Objects
            if strcmp(tmp(1), 'SubObject')
                if length(tmp) == 1
                    fprintf('Input File Error: SubComponent - Missing Input to Turn ON(1)/OFF(1) the SubObject Analysis\n');
                    return
                else
                    if strcmp(tmp(2), '1')
                        active_subobj = true;
                        fprintf('Input Setting: Sub-Object Anlysis is Turned ON.\n');
%                         if active_trans && strcmp(tmp(3), 'SubBridge')
%                             for ii = 4:length(tmp)
%                                 iobj = tmp{ii}; 
%                                 tmp2 = strcat('SubBridge_', iobj);
%                                 [Power, Trans, Comm, Neighborhood] = Library.readInput(tmp2, pow_check, comm_check, trans_check, Power, Trans, Comm, Dictionary, Neighborhood);
%                             end
                        if active_trans 

                                [Power, Trans, Comm, Neighborhood] = Library.readInput('SubBridge_1', pow_check, comm_check, trans_check, Power, Trans, Comm, Dictionary, Neighborhood);

                        elseif  ~active_trans   
                            fprintf('Inconsistent Input Setting! Transportation System is turned OFF; therefore, SubBridge Analysis can NOT be turned ON.');
                            return
                        else
                            fprintf('SubObjects for this type of object has not been defined yet. This command is disregarded in the following simulation.');
                            return
                        end 
                        
                    elseif strcmp(tmp(2), '0')
                        active_subobj = false;
                        disp('Input Setting: Sub-Object Anlysis is Turned OFF.');
                    else
                        fprintf('Input File Error: SubObject - Input not valid (0/1)\n');
                        return
                    end

                end
            end
            
            %=== set up the direction of transportation links
            if active_trans && strcmp(tmp(1), 'TransportationLinkDirection')
                if length(tmp) == 1
                    fprintf('Input File Error: TransportationLinkDirection - Missing Input\n');
                    return
                else
                    if strcmp(tmp(2), '1')     % directional in the network graph 
                        LinkDirectionChoice = 1;
                    elseif strcmp(tmp(2), '0') % not directional in the network graph 
                        LinkDirectionChoice = 0; 
                    else
                        fprintf('Input File Error: TransportationLinkDirection - Incompatible Input\n');
                        return
                    end                 
                   
                end
            end 
            
            %=== set up the time horizon
            if strcmp(tmp(1), 'TimeHorizon')
                if length(tmp) == 1
                    fprintf('Input File Error: TimeHorizon - Missing Input\n');
                    return
                else
                    time_horizon = str2num(char(tmp(2)));
                    if time_horizon <10
                        time_horizon=10;
                    end
                end
            end
            
            %=== set up the file name of neighborhood 
            % a neighborhood means a local community as a resident/commercial region.
            % The object of "neighborhoog" has socioeconomic properties, such as population.
            if strcmp(tmp(1), 'DataNeighbor')
                try
                    if length(tmp) == 1
                        fprintf('Input File Error: Neighborhood - Missing Input of the File Name\n');
                        return
                    else
                        for ii = 2:length(tmp)
                            [Power, Trans, Comm, Neighborhood] = Library.readInput(char(tmp(ii)), pow_check, comm_check, trans_check, Power, Trans, Comm, Dictionary, Neighborhood);
                        end
                        active_neighbor = true;
                    end
                catch exception
                    msg = getReport(exception, 'basic');
                    disp(msg);
                    return;
                end
            end

            %=== Number of component damage sampling
            if strcmp(tmp(1), 'Nsamples')
                if length(tmp) == 1
                    fprintf('Input File Error: Nsamples - Missing Number Of Damage Scenario Sample input\n');
                    return
                else
                    Nsamples = str2num(char(tmp(2)));
                end
            end
            
            %=== Number of restoration sampling
            if strcmp(tmp(1), 'NRun')
                if length(tmp) == 1
                    fprintf('Input File Error: NRun - Missing Number Of Task Duration Sample input\n');
                    return
                else
                    NRun = str2num(char(tmp(2)));
                end
            end

            
            %=== set up the profile number (parallel computing?)
            if strcmp(tmp(1), 'Profile_Num')
                if length(tmp) == 1
                    fprintf('Input File Error: Profile_Num - Missing Input\n');
                    return
                else
                    Profile_Num = str2num(char(tmp(2)));
                end
            end
            
            %=== select the scheme of restoration planning
            % (1) Scheme 3A: Scheduler_Num = 1, Policy scheme
            % (2) Scheme 3B: Scheduler_Num = 2, Optimization scheme
            if strcmp(tmp(1), 'Scheduler')
                if length(tmp) == 1
                    fprintf('Input File Error: Scheduler - Missing Scheduler Model Input\n');
                    return
                else
                    Scheduler_Num = str2num(char(tmp(2)));
                end
            end
            
            %=== select the optimization algorithm
            if strcmp(tmp(1), 'OptimizationChoice')
                if length(tmp) == 1
                    fprintf('Input File Error: OptimizationChoice - Missing Input of Optimization Algorithm\n');
                    return
                else
                    OptimizationChoice = str2num(char(tmp(2)));
                end
            end

            %=== Diff_unit: setting up time unit of task
            % (1) uniform time unit: Diff_unit = 0
            %    Power and Communication systems: hour
            %    Transportation system: day
            % (2) non-uniform time unit: Diff_unit = 1
            %    Phase 1: t=0~3 days, TimeUnit = hour
            %    Phase 2: t=3~30 days, TimeUnit = day
            %    Phase 3: t=30~180 days, TimeUnit = week
            %    Phase 4: t=180~timehorizon days, TimeUnit = month
            if strcmp(tmp(1), 'DifferentTimeUnit')
                if length(tmp) == 1  
                    Diff_unit = 0; 
                else                 
                    Diff_unit = str2num(char(tmp(2)));
                end
            end
            
            %=== Cust_unit: setting up disaster manament phases (in the unit of days)
            if strcmp(tmp(1), 'RestorationPhaseDay')
                if length(tmp) == 1 && Diff_unit == 1
                    Cust_unit(1)=3; Cust_unit(2)=28; Cust_unit(3)=168;
                else
                     for ii = 2:length(tmp)
                         Cust_unit(ii-1) = str2num(char(tmp(ii)));
                     end
                end
                if ~issorted(Cust_unit)
                    fprintf('Input File Error: Cust_unit - Input not Valid\n');
                    return
                end
                if Diff_unit
                    if (max(Cust_unit)>time_horizon || min(Cust_unit)<0)
                        fprintf('Input File Error: Cust_unit - Input not Valid\n');
                        return
                    end
                end
            end
                
            %=== set the functionality metric for the power system
            if active_power && strcmp(tmp(1), 'Functionality_Power')
                if length(tmp) == 1
                    fprintf('Input File Error: Functionality_Power - Missing Functionality Metric Input\n');
                    return
                else
                    for ii = 2:length(tmp)
                        Power_Func_Num(ii-1) = str2num(char(tmp(ii)));
                    end
                end
            elseif ~active_power
                Power_Func_Num = 0;
            end
            
            %=== set the functionality metric for the communication system
            if active_comm && strcmp(tmp(1), 'Functionality_Comm')
                if length(tmp) == 1
                    fprintf('Input File Error: Functionality_Comm - Missing Functionality Metric Input\n');
                    return
                else
                    for ii = 2:length(tmp)
                        Comm_Func_Num(ii-1) = str2num(char(tmp(ii)));
                    end
                end
            elseif ~active_comm
                Comm_Func_Num = 0;
            end
            
            %=== set the functionality metric for the transportation system
            if active_trans && strcmp(tmp(1), 'Functionality_Trans')
                if length(tmp) == 1
                    fprintf('Input File Error: Functionality_Transportation - Missing Functionality Metric Input\n');
                    return
                else
                    for ii = 2:length(tmp)
                        Trans_Func_Num(ii-1) = str2num(char(tmp(ii)));
                    end
                end
            elseif ~active_trans
                Trans_Func_Num = 0;
            end
            
            %=== turn on(1)/off(0) rescheduling
            if strcmp(tmp(1), 'ReSchedule')
                if length(tmp) == 1
                    fprintf('Input File Error: ReSchedule - Missing Input\n');
                    return
                else
                    ReSchedule_Num = str2num(char(tmp(2)));
                end
            end
            
            %=== turn on(1)/off(0) the interdependency of the functionality
            % delay effect 
            % (1) Interdependence=1, there is the delaying effect to other
            % systems due to the functionality disruption of the transportation
            % system.
            % (2) Interdependence=0, there is no delaying effect to other
            % systems due to the functionality disruption of the transportation
            % system.
            if strcmp(tmp(1), 'InterdependenceTransDelay')
                if length(tmp) == 1
                    fprintf('Input File Error: InterdependenceTransDelay - Missing Input\n');
                    return
                else
                    if strcmp(tmp(2), '1') 
                        Interdependence_Num = 1;
                        System_Dependent_Factor = tmp{3};
                        Qtrans0 = str2num(tmp{4}); % threshold of transportation functionality to consider restoration delays 
                        msg = strcat('Input Setting: The Delay Effect is turned ON, with Qtrans0=', num2str(Qtrans0),'.');
                        disp(msg);
                    elseif strcmp(tmp(2), '0') 
                        Interdependence_Num = 0;
                        System_Dependent_Factor = 1;
                        Qtrans0 = 0;
                        msg = 'Input Setting: The Delay Effect is turned OFF.';
                        disp(msg);
                    else 
                        fprintf('Input File Error: InterdependeceTransDelay - Invalid Input\n');
                        return
                    end
                end
            end
            
            %=== turn on(1)/off(0) the inter-system functionality dependency 
            if strcmp(tmp(1), 'InterdependenceFunctionality')
                if length(tmp) == 1
                    fprintf('Input File Error: InterdependenceFunctionality - Missing Input\n');
                    return
                else
                    % Inter-system Functionality Dependency is turned on.
                    if strcmp(tmp(2), '1') 
                        InterdependenceFunc = 1;
                        msg = 'Input Setting: Inter-system Functionality Dependency is turned ON.';
                        disp(msg);
                    % Inter-system Functionality Dependency is turned off.    
                    elseif strcmp(tmp(2), '0')
                        InterdependenceFunc = 0;
                        msg = 'Input Setting: Inter-system Functionality Dependency is turned OFF.';
                        disp(msg);
                    else 
                        fprintf('Input File Error: InterdependenceFunctionality - Invalid Input\n');
                        return    
                    end
                end
            end
            
            %=== turn on(1)/off(0) the precedence dependency 
            if strcmp(tmp(1), 'InterdependencePrecedence')
                if length(tmp) == 1
                    fprintf('Input File Error: InterdependencePrecedence - Missing Input\n');
                    return
                else
                    if strcmp(tmp(2), '1') 
                        InterdependencePrec = 1;
                        msg = 'Input Setting: Inter-system Precedence Dependency is turned ON.';
                        disp(msg);
                    elseif strcmp(tmp(2), '0')
                        InterdependencePrec = 0;
                        msg = 'Input Setting: Inter-system Precedence Dependency is turned OFF.';
                        disp(msg);
                    else 
                        fprintf('Input File Error: InterdependencePrecedence - Invalid Input\n');
                        return    
                    end
                end
            end
            
            %=== the avaiability of four types of restoration resources for 
            %the power system
            if active_power && strcmp(tmp(1), 'Power_Resource')
                if length(tmp) == 1
                    fprintf('Input File Error: Power_Resource - Missing Input\n');
                    return
                else
                    for ii = 1 + 1:length(tmp)
                        Power_Resource(ii - 1) = str2num(char(tmp(ii)));
                    end
                end
            end
            
            %=== the avaiability of four types of restoration resources for 
            %the communication system
            if active_comm && strcmp(tmp(1), 'Comm_Resource')
                if length(tmp) == 1
                    fprintf('Input File Error: Comm_Resource - Missing Input\n');
                    return
                else
                    for ii = 1 + 1:length(tmp)
                        Comm_Resource(ii - 1) = str2num(char(tmp(ii)));
                    end
                end
            end
            
            %=== the avaiability of four types of restoration resources for 
            %the transportation system
            if active_trans && strcmp(tmp(1), 'Trans_Resource')
                if length(tmp) == 1
                    fprintf('Input File Error: Trans_Resource - Missing Input\n');
                    return
                else
                    for ii = 1 + 1:length(tmp)
                        Trans_Resource(ii - 1) = str2num(char(tmp(ii)));
                    end
                end
            end

            %===  Set the resilience metric for the three systems 
            if active_power && strcmp(tmp(1), 'Resilience_Power')
                if length(tmp) == 1
                    fprintf('Input File Error: Resilience_Power - Missing Metric Input\n');
                    return
                else
                    for ii = 2:length(tmp)
                        ResiliencePower(ii-1) = str2num(char(tmp(ii)));
                    end
                end
            elseif ~active_power 
                ResiliencePower = 0;
            end
            if active_comm && strcmp(tmp(1), 'Resilience_Comm')
                if length(tmp) == 1
                    fprintf('Input File Error: Resilience_Comm - Missing Metric Input\n');
                    return
                else
                    for ii = 2:length(tmp)
                        ResilienceComm(ii-1) = str2num(char(tmp(ii)));
                    end
                end
            elseif ~active_comm 
                ResilienceComm = 0;    
            end
            if active_trans && strcmp(tmp(1), 'Resilience_Trans')
                if length(tmp) == 1
                    fprintf('Input File Error: Resilience_Trans - Missing Metric Input\n');
                    return
                else
                    for ii = 2:length(tmp)
                        ResilienceTrans(ii-1) = str2num(char(tmp(ii)));
                    end
                end
            elseif ~active_trans 
                ResilienceTrans = 0;    
            end
            
            %=== Probability of the presence of magic battery
            % to relieve the original functionality dependency of components 
            % that are used to use the electricy of the power system by 
            % using the electricity from the magic battery.   
            % By default, Prob_Magic_Battery = 0.9
            if strcmp(tmp(1), 'Prob_Magic_Battery')
                if length(tmp) == 1
                    Prob_Magic_Battery = 0.9;
                else
                    Prob_Magic_Battery = str2num(char(tmp(2)));
                end
            end
            
            if strcmp(tmp(1), 'Save_history')
                if length(tmp) == 1
                    Save_history = 1;
                else
                    Save_history = str2num(char(tmp(2)));
                end
            end            
            
            %=== Do all systems develop restoration plans at Step 3
            %seperately? 
            % (1) Seperate_Scheduling = 1(default), yes, every system develops 
            % the restoration plan seperately.
            % (2) Seperate_Scheduling = 0, no, all restoration tasks from 
            % different systems put together to develop one restoration plan.
            if strcmp(tmp(1), 'Seperate_Scheduling')
                if length(tmp) == 1
                    Seperate_Scheduling = 1;
                else
                    Seperate_Scheduling = str2num(char(tmp(2)));
                end
            end
            
            %=== Scheme 3A (priority): select the priority criterion for the power system 
            if active_power && strcmp(tmp(1), 'Priority_power_num')
                if length(tmp) == 1
                    priority_power_num = 1;
                else
                    priority_power_num = str2num(char(tmp(2)));
                    % So far there are only 2 pority rules implemented,
                    % therefore, any priority number >2 is not accpetable
                    % at the moment.
                    if priority_power_num>2
                        msg = 'Input File Error: Input priority_power_num is not recognized.';
                        disp(msg);
                    end
                end
            elseif ~active_power
                priority_power_num = 0;
            end

            
            %=== Scheme 3A (priority): select the priority criterion for the communication system 
            if active_comm & strcmp(tmp(1), 'Priority_communication_num')
                if length(tmp) == 1
                    priority_communication_num = 1;
                else
                    priority_communication_num = str2num(char(tmp(2)));
                    % So far there are only 2 pority rules implemented,
                    % therefore, any priority number >2 is not accpetable
                    % at the moment.
                    if priority_communication_num>2
                        msg = 'Input File Error: Input priority_communication_num is not recognized.';
                        disp(msg);
                    end
                end
            elseif ~active_comm
                priority_communication_num = 0;
            end
            
            %=== Scheme 3A (priority): select the priority criterion for the transportation system 
            if active_trans && strcmp(tmp(1), 'Priority_transportation_num')
                if length(tmp) == 1
                    priority_transportation_num = 1;
                else
                    priority_transportation_num = str2num(char(tmp(2)));
                    % So far there are only 2 pority rules implemented,
                    % therefore, any priority number >2 is not accpetable
                    % at the moment.
                    if priority_transportation_num>2
                        msg = 'Input File Error: Input priority_transportation_num is not recognized.';
                        disp(msg);
                    end
                end
            elseif ~active_trans
                priority_transportation_num = 0;
            end
            
            %=== Percentage of completed tasks to initiate the rescheduling
            %Percentage = No. of completed tasks at t=t / No. of all tasks at t=0 
            if strcmp(tmp(1), 'RescheduleAtPercentage')
                if size(str2num(char(tmp(end))),1)==0
                    Per_New_reschedule =0;
                else
                    if length(tmp)>=2
                        for ii = 1 + 1:length(tmp)
                        Per_New_reschedule(ii-1) = str2num(char(tmp(ii)));
                        end
                    else
                        Per_New_reschedule = str2num(char(tmp(2)));
                    end
                end
                if issorted(Per_New_reschedule)
                else
                    fprintf('Input File Error: Restoration at Percentage of Finished Input not Valid\n');
                    return
                end
                if (max(Per_New_reschedule)>1 || min(Per_New_reschedule)<0)
                    fprintf('Input File Error: Restoration at Percentage of Finished Input not Valid\n');
                    return
                end
   
            end            
 
            %=== Time to reschedule as different restoration phases (in the time unit of days)
            if strcmp(tmp(1), 'RescheduleAtTime')
                if size(str2num(char(tmp(2))),1)==0
                    Num_stage =0;
                else
                    if length(tmp)>=2
                        for ii = 1 + 1:length(tmp)
                            Num_stage(ii-1) = str2num(char(tmp(ii)));
                        end
                    else
                            Num_stage = str2num(char(tmp(2)));
                    end
                end
                if ~issorted(Num_stage)
                    fprintf('Input File Error: Restoration at Different Time Phase Input not Valid\n');
                    return
                end
                if ReSchedule_Num
                    if (max(Num_stage)>time_horizon || min(Num_stage)<0)
                        fprintf('Input File Error: Restoration at Different Time Phase Input not Valid\n');
                        return
                    end
                end
            end 
            
            %=== Turn on(1)/off(0) random sampling (active_rands) and
            % Set up the random seed (Seed) for random sampling.
            if strcmp(tmp(1), 'RandomChoice')
                if length(tmp) == 1
                    active_truerands = 0;
                else
                    active_truerands = str2num(char(tmp(2)));
                    
                    if active_truerands == 1
                        Seed = [];
                        fprintf('Input Setting: True Randomness is turned ON. \n');
                    elseif active_truerands == 0
                        if ~isempty(tmp(3))
                            Seed = str2num(char(tmp(3)));
                        else
                            Seed = 110;
                        end
                        msg = strcat('Input Setting: Pseudo-randomness is turned ON. The random seed is set as "', num2str(Seed),'".');
                        disp(msg);
                    else
                        Seed = [];
                        msg = 'Input File Error: The input of RandomChoice is not valid, please choose either 0 (OFF) or 1 (ON).'; 
                        disp(msg);
                        return
                    end
                    
                end
            end
            

            
       end
        line = fgetl(fid);
    catch exception
        msg = getReport(exception, 'basic');
        disp(msg);
        break;
    end
    
end

%% Save history results or not
if exist('PATH_dir','var')
    Hostname = 'Output/Results0'; 
    if Save_history==0
        if exist(fullfile(PATH_dir,Hostname), 'dir')
            rmdir('Output','s')
        else
            mkdir(fullfile(PATH_dir,Hostname));
        end
    else
        if ~exist(fullfile(PATH_dir,Hostname), 'dir')
            mkdir(fullfile(PATH_dir,Hostname)); 
        else
            files=dir(fullfile(PATH_dir,Hostname(1:6)));
            all_re = {files([files(:).isdir]).name}; all_re(ismember(all_re,{'.','..'})) = [];
            Hostname=['Output/Results',num2str(size(all_re,2))];
            mkdir(fullfile(PATH_dir,Hostname)); 
        end
    end
    Hostname=fullfile(PATH_dir,Hostname);
else
    Hostname = 'Output/Results0'; 
    if Save_history==0
        % Library.cmd_rmdir('Output')
        % rmdir('Output','s');
        if  ~exist(Hostname, 'dir')
            mkdir(Hostname);
        else
            rmdir Output s;
            mkdir(Hostname);
        end
    else
        if ~exist(Hostname, 'dir')
            mkdir(Hostname); 
        else
            files=dir(Hostname(1:6));
            all_re = {files([files(:).isdir]).name}; all_re(ismember(all_re,{'.','..'})) = [];
            Hostname=['Output/Results',num2str(size(all_re,2))];
            mkdir(Hostname); 
        end
    end    
end    
Library.CleanOldData(Hostname);
Library.CreateFolder(Hostname);

%% Restoration Resource Constraints for Every System
RestorationResource = [Power_Resource; Comm_Resource; Trans_Resource]; 

% if (Per_New_reschedule(1)==0 && Num_stage(1)==0)
%     ReSchedule_Num=0;
% end

%% Flag variable indicating every system is turn ON(1)/OFF(0)
ActiveSystem = [active_power, active_comm, active_trans];


%% Build new objects as needed 
if active_power
    Power = Library.CreateTransmissionTower(Power, Dictionary);
end
if active_comm
    Comm = Library.assignCellLine(Comm, Dictionary);
end
if active_trans
    Trans = Library.assignRoadToRoadNode(Trans,Dictionary);
end
%% Add the information of functionality dependency from communciation and transportation to power
if active_power
    [Trans, Comm] = Library.assignPowerToTransComm(Power, Trans,Comm, Dictionary);
end

%% Add virtual links from every neighborhood to the nearest object in every system
if active_neighbor  
    [Trans, Comm, Power, Neighborhood] = Library.CreateNeighborhoodConnection(Power, Trans, Comm, Dictionary, Neighborhood, active_power, active_trans, active_comm);
end

%% Assign parameters (i.e., mu and sigma of every lognormal curve) of fragility curves for every object
Library.assignFragility(EventType, IM, IMx, IMy, IMmeta, Power, Trans, Comm);


%% Build the initial graph for every system network based on Graph Theory
LinkDirectionChoicePower = 0; % false, meaning not directed
LinkDirectionChoiceComm = 0; % false, meaning not directed
LinkDirectionChoiceTrans = LinkDirectionChoice;
LinkDirection = [LinkDirectionChoicePower, LinkDirectionChoiceComm, LinkDirectionChoiceTrans]; 

if active_power
    powerGraph = Library.BuildGraphPower(Power, Dictionary, LinkDirectionChoicePower);
else
    powerGraph = {};
end
if active_comm
    commGraph = Library.BuildGraphComm(Comm, Dictionary, LinkDirectionChoiceComm);
else
    commGraph = {};
end 
if active_trans
    transGraph = Library.BuildGraphTrans(Trans, Dictionary, LinkDirectionChoiceTrans);
else
    transGraph = {};
end

GraphSystem0{1,1} = powerGraph;
GraphSystem0{1,2} = commGraph;
GraphSystem0{1,3} = transGraph;

ResilienceMetricChoice{1} = ResiliencePower; 
ResilienceMetricChoice{2} = ResilienceComm; 
ResilienceMetricChoice{3} = ResilienceTrans; 

Branch_Set = Power{1};
Bus_Set = Power{2};
Generator_Set = Power{3};
TransmissionTower_Set = Power{4};
Neighborhood_Power_Set = Power{5};

Centraloffice_Set = Comm{1};
CommunicationTower_Set = Comm{2};
Cellline_Set = Comm{3};
Neighborhood_Comm_Set = Comm{4};

Road_Set = Trans{1};
Bridge_Set = Trans{2};
TrafficLight_Set = Trans{3};
RoadNode_Set = Trans{4};
Neighborhood_Trans_Set = Trans{5};

%=== Number of functionality metrics
nQmetric = [length(Power_Func_Num), length(Comm_Func_Num), length(Trans_Func_Num)];

%% prepare input data of the transportation system 
%=== Road Node
mydir = fullfile(myfolder,fn{4});
num = readtable([mydir,'.csv']); 
latitude = num{:,2}';
longtitude = num{:,3}';

%=== Set up the set of road link
for ii = 1:length(Road_Set)
    Road_Set{ii}.Start_Location = [latitude(Road_Set{ii}.Start_Node), longtitude(Road_Set{ii}.Start_Node)];
    Road_Set{ii}.End_Location = [latitude(Road_Set{ii}.End_Node), longtitude(Road_Set{ii}.End_Node)];
end

%% Clean variables 
clear Power Comm Trans fid index system i line tmp antenna_check comm_check pow_check trans_check centraloffice_check bridge_check road_check trafficlight_check communicationtower_check;


