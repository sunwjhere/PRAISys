% Data_Supplement_Real.m
% =====================================================================
% read hazard data and infrastructure data
% & prepare input variables and directories
% ==========set up the directory to call yalmip on Sol (HPC)=================
addpath(genpath('/share/Apps/yalmip/R20181012'))

%% === Set up the directory of output files
%[~, Hostname] = system('hostname');

% clear num txt raw i;
run('read_links.m');
run('Input_Real.m');


priority = 1;
change = 1;
prev = inf;

while change ~= 0 && active_trans
    index = [];
    max = -1;
    change = 0;
    for ii = 1:length(Bridge_Set)
        if Bridge_Set{ii}.Traffic >= max && Bridge_Set{ii}.Traffic < prev
            if Bridge_Set{ii}.Traffic > max
                index = ii;
            else
                index = [index, ii];
            end
            max = Bridge_Set{ii}.Traffic;
            change = 1;
        end
    end
    
    for ii = 1:length(index)
        Bridge_Set{index(ii)}.Priority = priority;
    end
    priority = priority + 1;
    prev = max;
end

Trans_Priority = priority - 1;

if active_trans
    for i = 1:length(Road_Set)
        if isempty(Road_Set{i}.Bridge_Carr)
            Road_Set{i}.Priority = priority - 1;
        else
            tmp = -1;
            for j = 1:length(Road_Set{i}.Bridge_Carr)
                if Bridge_Set{Road_Set{i}.Bridge_Carr(j)}.Priority >= max
                    tmp = Bridge_Set{Road_Set{i}.Bridge_Carr(j)}.Priority;
                end
            end
            Road_Set{i}.Priority = tmp;
        end

        for j = 1:length(Road_Set{i}.TrafficLight)
            TrafficLight_Set{Road_Set{i}.TrafficLight(j)}.Priority = Road_Set{i}.Priority;
        end
    end
end

priority = 1;
change = 1;
prev = inf;

while change ~= 0 && active_power
    index = [];
    max = -1;
    change = 0;
    for ii = 1:length(Bus_Set)
        if Bus_Set{ii}.Capacity >= max && Bus_Set{ii}.Capacity < prev
            if Bus_Set{ii}.Capacity > max
                index = ii;
            else
                index = [index, ii];
            end
            max = Bus_Set{ii}.Capacity;
            change = 1;
        end
    end
    
    for ii = 1:length(index)
        Bus_Set{index(ii)}.Priority = priority;
        if ~isempty(Bus_Set{index(ii)}.Generator)
            temp = extractAfter(Bus_Set{index(ii)}.Generator, 9);
            temp = str2num(temp);
            Generator_Set{temp}.Priority = priority;
        end
        for j = 1:length(Bus_Set{index(ii)}.Branch)
            Branch_Set{Bus_Set{index(ii)}.Branch(j)}.Priority = priority;
        end
    end
    priority = priority + 1;
    prev = max;
end

Power_Priority = priority - 1;
if active_comm
    for ii = 1:length(CommunicationTower_Set)
        CommunicationTower_Set{ii}.Priority = 1;
    end

    for ii = 1:length(Centraloffice_Set)
        Centraloffice_Set{ii}.Priority = 2;
    end
end
Comm_Priority = 3;


% ==== Profile
% Call for Parallel Computing Toolbox
% Execute code in parallel on workers of parallel pool
% spmd: single program, multiple data
if Profile_Num == 1
    profile on
    spmd  
        mpiprofile('on');
    end
end

% ==== Set up infrastructure system(s) 
%--- Save Original Data (Need to use load command to view the data)
Power_Set = {Branch_Set, Bus_Set, Generator_Set, TransmissionTower_Set, Neighborhood_Power_Set};
Communication_Set = {Centraloffice_Set, CommunicationTower_Set, Cellline_Set, Neighborhood_Comm_Set};
Transportation_Set = {Road_Set, Bridge_Set, TrafficLight_Set,RoadNode_Set,Neighborhood_Trans_Set};

Total = {Power_Set, Communication_Set, Transportation_Set, Dictionary, Neighborhood};

%--- save initial infrastructure data
fninfrastructure = strcat( deblank(Hostname), '/Original_Data.mat');
save(fninfrastructure, 'Total');

clear i j index max prev priority temp change tmp latitude longtitude;

%--- For Debugging 
% Nobj = length(Power_Set{1})
% Nobj = Nobj + length(Power_Set{2})
% Nobj = Nobj + length(Power_Set{3})
% Nobj = Nobj + length(Communication_Set{1})
% Nobj = Nobj + length(Communication_Set{2})
% Nobj = Nobj + length(Communication_Set{3})
% Nobj = Nobj + length(Transportation_Set{1})
% Nobj = Nobj + length(Transportation_Set{2})
% Nobj = Nobj + length(Transportation_Set{3})
