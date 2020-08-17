classdef Library
    methods(Static)

        %% Count the total number of variables needed for the random samples following the uniform distribution ([0,1]) with Latin Hybercube Sampling (LHS) (STEP 1)
        function [SamplesUniform, SamplesUnifDmg, SamplesUnifDur] = SampleInititalUniformSample(Ndamagesample, Ntaskduration, Total, Seed, active_truerands)
            %==========================================================================
            % function SampleInititalUniformSample
            %---Input---
            % Ndamagesample = Number of samples for every component damage
            % state
            % Ntaskduration = Number of task druation samples for every damage
            % scenario (a possible combination of component damage state)
            % Total = the cell of Total={Power_Set, Communication_Set, Transportation_Set, Dictionary, Neighborhood};
            %---Output---
            % SampleUniform = samples following the uniform distributions
            %---Computation procedure 
            % X = lhsdesign(n,p)
            % n = No. of all random samples = Ndmgsamples (Nsamplescomponent damage) * Ntsk (Nrun: task duration samples per damage state)
            % p = No. of variables = Nobj + Nobj(>=Ndmgsample)*Ntsk + Nbat
            % = the total number of damaged component + the total number of
            % tasks + the total number of object magic battery
            % Nobj = No. of objects in infrastructures
            % Ndmgsamples: No. of component damage samples  
            % Ntsk = No. of restoration task to repair a damaged component, assuming 
            % that every damaged component has 10 tasks at a damage state at most. 
            % Nbat = No. of magic battery (<=Nobj)
            %==========================================================================
            % SystemPower = Total{1};
            % SystemCommu = Total{2};
            % SystemTrans = Total{3}; 
            %--- count the number of samples NSample to build samples that
            %follow uniform distributions ([0,1]), to be mapped to other distributions
            % NSample = No. of damage samples for components + No. of task duration for components +
            % No.damage scenario samples for sub-components + No. of task duration for sub-components
            % samples
            % Ndmgsample = No. of damage samples 

            Ndmgsample = Ndamagesample;
            Ndur = Ntaskduration; 
            NSample = Ndmgsample*Ndur;     
            
            % assuming that every damaged object at a certain damage state requires at most 10 tasks
            Ntsk = 10; 
            
            %--- count the number of objects Nobj
            Nobj = 0;
            for isys = 1:3
                sys = Total{isys};
                % length(sys)-1: 
                %(1)assuming that neighbors (in the last cell) do not fail
                for iob = 1:length(sys)-1 
                    ObjectSet = sys{iob};
                    Nobj = Nobj + length(ObjectSet);
                end
            end
            
            %--- count the number of sub-objects NobjSub
            % for BridgeSet
            NobjSub = 0;
            isys = 3; iobjtype = 2; 
            ObjSet = Total{isys}{iobjtype};
            nobj = length(ObjSet);
            for iob = 1:nobj
                SubObjSet = ObjSet{iob}.ColumnSet;
                NobjSub = NobjSub + length(SubObjSet);
                SubObjSet = ObjSet{iob}.ColumnFoundSet;
                NobjSub = NobjSub + length(SubObjSet);
                SubObjSet = ObjSet{iob}.AbutmentSet;
                NobjSub = NobjSub + length(SubObjSet);
                SubObjSet = ObjSet{iob}.AbutmentFoundSet;
                NobjSub = NobjSub + length(SubObjSet);
                SubObjSet = ObjSet{iob}.GirderSet;
                NobjSub = NobjSub + length(SubObjSet);
                SubObjSet = ObjSet{iob}.BearingSet;
                NobjSub = NobjSub + length(SubObjSet);
                SubObjSet = ObjSet{iob}.DeckSet;
                NobjSub = NobjSub + length(SubObjSet);
                SubObjSet = ObjSet{iob}.SlabSet;
                NobjSub = NobjSub + length(SubObjSet);
            end
            
            
            %--- count the number of magic battery Nbat
            % determining whether a component has a magic battery or not.
            % Components w/ or w/o a magic battery:
            % In the communciation system
            % 2.1 Central Office
            % 2.2 Communication tower
            % 2.3 Cell line 
            % In the transportation system
            % 3.1 Traffic light
            
            Nbat = 0;
            
            %- 2. Communciation system
            isys = 2; 
            sys = Total{isys};
            for iob = 1:length(sys)-1
                ObjectSet = sys{iob};
                Nbat = Nbat + length(ObjectSet);
            end
            %- 3. Transportation system
            isys = 3; idx = 3;
            TLSet = Total{isys}{idx}; 
            Nbat = Nbat + length(TLSet);
            
            %--- count the total number of variables NVariable
            NVariable = Nobj + NobjSub + Nbat + (Nobj+NobjSub)*Ntsk;
            
                        
            %--- set random seed for generating the same random number
            if ~active_truerands  
                if isempty(Seed)
                    Seed = 110;
                end                    
                rng(Seed);
            end

            %--- generate uniform samples [0,1] with Latin hypercube sampling
            % the default samples follow uniform distributions between 0
            % and 1 
            SamplesUniform = lhsdesign(NSample,NVariable); 
            % Damage scenario samples
            is = 1;
            it = Nobj + NobjSub + Nbat;
            SamplesUnifDmg = repmat(SamplesUniform(1:Ndmgsample, is:it),Ndur,1);
            % Task duration samples 
            is = it+1;
            it = NVariable;
            SamplesUnifDur = SamplesUniform(:, is:it);

        end

        %% Extract uniform samples from Latin Hypercube Sampling before mapping for individual task durations (STEP 4)
        function [SamplesUnifDur] = GetDurationSample(Total, SamplesUniform)
            %=== Field
            % assuming that every damaged object at a certain damage state requires at most 10 tasks
            Ntsk = 10; 
            
            %--- count the number of objects Nobj
            Nobj = 0;
            for isys = 1:3
                sys = Total{isys};
                %-- length(sys)-1: 
                %(1)assuming that neighbors (in the last cell) do not fail
                for iob = 1:length(sys)-1 
                    ObjectSet = sys{iob};
                    Nobj = Nobj + length(ObjectSet);
                end
            end
            
            %--- The original size
            SamplesUnifDur = SamplesUniform(:, end-Nobj*Ntsk+1:end); 
        
        end


        %% revised function for component damage assessment  (STEP 2)
        function [GraphSystem, TotalDamage, RestorationTable, NSysTask] = SampleDamage(EventType, Nsamples, NRun, Hostname, ...
                IM, IMx, IMy, ActiveSystem, ...
                powerGraph, commGraph, transGraph, Prob_Magic_Battery, ...
                SamplesUnifDmg, Seperate_Scheduling, InterdependencePrec)
            %==========================================================================
            % function SampleDamage
            % 1. fragility analysis at component-level with fragility curves (damage status evaluation)
            % 2. cascading failure analysis (stop status evaluation)
            % 3. determining the corresponding restoration tasks for repairing the damaged components
            %--- Input
            % EventType = 1-earthquake, 2-hurricane winds, 3-flood 
            % Nsamples = No. of samples for representing component damage state 
            % NRun = No. of samaples for representing task duration
            % Hostname = the folder name of output files 
            % IM, IMx, IMy = intensity measurement variables and the corresponding locations in the directions of x (longitude) and y (latitude) 
            % active_power, active_comm, active_trans = flags (0-off/1-on) of active systems
            % transGraph, powerGraph, commGraph = network graph of every system at the initial state 
            % Prob_Magic_Battery = the probability of magic battery from the input 
            % SamplesUnifDmg = random damage samples following uniform distributions from Latin Hybercube Sampling 
            %--- Output
            % GraphSystem = {powerGraph, commGraph, transGraph} = network graph of every system after damage assessment
            % TotalDamage = a cell collection of damage scenarios for three systems, the Dictionary, the neighborhoods 
            % RestorationTable = a cell collection of restoration tasks for
            % the damage scenarios of three systems            
            % Nsystask = Total number of restoration task for every damage
            % scenario sample of three systems
            %==========================================================================
            %=== Field
            TotalDamage = [];
            GraphSystem = [];
            RestorationTable = [];
            NSysTask = [0 0 0];
            
            %=== flag variable to indicate whether a system is turn ON(1)/OFF(0)
            active_power = ActiveSystem(1);
            active_comm = ActiveSystem(2);
            active_trans = ActiveSystem(3); 
            
            %=== Set up a hash table (dictionary) 
            [sumTaskHash, sumDamageTaskHash] = Library.SetUpHashTables();
            
            %=== Set up the initial systems (original undamaged infrastructure data) before damage assessment 
            fninput = strcat( deblank(Hostname), '/Original_Data.mat');
            load(fninput);
            
            %--- NDamageVariables = Nobj+Nbat
            samples = SamplesUnifDmg; 
            NSampleDmg = Nsamples*NRun;
            
            %=== Create and save samples of component damage state
            % and creat tables of tasks and precedence relations
            for isample = 1:NSampleDmg
                taskIndex = 1;
                %--- Create Infrastructure Systems at the initial functionality 
                % state, which are saved in "Total" in "Original_Data.mat"
                clearvars Total; a = load(fninput); Total = a.Total;
                [SystemSetPower, SystemSetCommu, SystemSetTrans, Dictionary, Neighborhood] = Library.ResetData(Total);

                %--- Evaluation of Damage and Stop (cascading failure) states
                [powerGraph, commGraph, transGraph, NSysTask(isample,:), taskIndex] = Library.DamageAndStop(EventType, IM, IMx, IMy, active_power, active_comm, active_trans, SystemSetPower, SystemSetCommu, SystemSetTrans, Dictionary, Neighborhood,transGraph, powerGraph, commGraph,sumTaskHash, sumDamageTaskHash,taskIndex,Prob_Magic_Battery, samples, isample);  
                
                if ~any(NSysTask(isample,:))
                    msg = strcat('Function SampleDamage Warning: No damage, therefore no restoration task, in sample #', num2str(isample), '. Stop here!');
                    disp(msg);
                    continue
                end
                
                GraphSystem{isample,1}{1} = powerGraph;
                GraphSystem{isample,1}{2} = commGraph;
                GraphSystem{isample,1}{3} = transGraph;
                
                % Damage samples: update the initial damage state of every object in "Total" 
                Total = {SystemSetPower, SystemSetCommu, SystemSetTrans, Dictionary, Neighborhood};        
                TotalDamage{isample,1} = Total;
                  
                %--- Set up tables of task and precedence
                %-Strong precedence: Spre
                if Seperate_Scheduling == 0
                    [taskTable,index] = Library.createTaskTable(Dictionary);
                    precedenceTable = Library.createPreTable(Dictionary,index, active_power, active_comm, active_trans);
                else
                    taskTable = Interface1.createTaskTableSep(Dictionary, active_power, active_comm, active_trans);
                    precedenceTable = Interface1.createPreTableSep(Dictionary, active_power, active_comm, active_trans);
                end
               
                %-Weak precedence: Wpre, remove all precedence
                % reset the precedence table as zeros to relax precedence dependency
                if InterdependencePrec == 0
                    for isys = 1:length(precedenceTable)
                        n = size(precedenceTable{isys},1); 
                        pre{isys} = num2cell(zeros(n, n)); 
                        for itask = 1:n
                            pre{isys}{1,itask} = precedenceTable{isys}{1,itask};
                            pre{isys}{itask,1} = precedenceTable{isys}{itask,1};
                        end
                    end
                    precedenceTable = pre;
                end
                
                %--- create dictionarys of task and precedence for every
                % system
                [Dic_p_task,Dic_c_task,Dic_t_task,Dic_p_prece,Dic_c_prece,Dic_t_prece] = Library.createDictionaryTaskAndPrecedence(taskTable,precedenceTable);
                
                RestorationTable{isample,1} = taskTable;
                RestorationTable{isample,2} = precedenceTable;
                RestorationTable{isample,1}{2,1} = Dic_p_task;
                RestorationTable{isample,1}{2,2} = Dic_c_task;
                RestorationTable{isample,1}{2,3} = Dic_t_task;
                RestorationTable{isample,2}{2,1} = Dic_p_prece;
                RestorationTable{isample,2}{2,2} = Dic_c_prece;
                RestorationTable{isample,2}{2,3} = Dic_t_prece;

                %--- save data files
                % save initial damage state in the txt file and mat files    
                Library.SaveLogDamage(isample, SystemSetPower, SystemSetCommu, SystemSetTrans, Hostname);      
                fnoutput = strcat( deblank(Hostname), '/mat/Data_DamageScenarioSample_', num2str(isample)); 
                save(fnoutput, 'Total');
                
                % save task-related data in the mat files 
                fnoutput = strcat( deblank(Hostname), '/task/Data_RestorationTable_', num2str(isample)); 
                save(fnoutput, 'RestorationTable');
                
                %--- display the message of finishign damage sampling 
                msg = strcat('----Function SampleDamage: Finish Sampling DamageSample #', num2str(isample), '----');
                disp(msg); 

            end %isample: the index of system damage scenario sample

            
%             DamagedComponent = TotalDamage;
%             save(strcat( deblank(Hostname), '/mat/DamageScenario.mat'), 'DamagedComponent', 'GraphSystem');
%             save(strcat( deblank(Hostname), '/mat/TaskTable.mat'), 'RestorationTable');

        end % function SampleDamage 
        
        
        %% revised function of DamageAndStop (STEP 2)
        function [powerGraph, commGraph, transGraph, NSysTask, taskIndex] = DamageAndStop(EventType, IM, IMx, IMy, active_power, active_comm, active_trans, Power_Set, Communication_Set, Transportation_Set, Dictionary, Neighborhood, transGraph, powerGraph, commGraph, sumTaskHash, sumDamageTaskHash, taskIndex, Prob_Magic_Battery, samples, isample)
            % =====================================================================
            % function DamageAndStop
            % Initial Damage Evaluation, Evaluate Stoped or Cascading Damaged
            % Compoment and Calculate the initial Functionnality for each Road
            % =====================================================================
     
            %=== Field
            % Starting Index (is) and Ending Index (it) of different damage sample in the initial damage samples following uniform distributions from Latin Hypercube Sampling 
            is = 0;
            it = 0;
            
            %=== Initial Damage Evaluation at the Object level
            %---- Power_Set = {Branch_Set, Bus_Set, Generator_Set, TransmissionTower_Set, Neighborhood_Power_Set};
            % damagable components: 1-ranch, 2-bus, 4-transmission tower
            if active_power
                %taskIndex = 1;
                SystemSet = Power_Set;
                idx = [1,2,4];
                for ii = idx
                    if ~isempty(SystemSet{ii})
                        nobj = length(SystemSet{ii});
                        is = it+1; % starting column index = starting variable index 
                        it = is+nobj-1;  % ending column index = ending variable index 
                        sample = samples(isample,is:it); 
                        taskIndex = Library.DamageEval(SystemSet{ii},EventType,IM,IMx,IMy,Dictionary,sumTaskHash, sumDamageTaskHash,taskIndex,Prob_Magic_Battery, sample);
                        
                    end
                end
            end
            % No. of tasks to repair damaged objects in the Power system
            Ntaskp = taskIndex; 
            
            %---- Communication_Set = {Centraloffice_Set, CommunicationTower_Set, Cellline_Set, Neighborhood_Comm_Set};
            % damagable components: 1-central office, 2-communication tower,  3-cell line,
            % therefore idx = [1,2,3];
            if active_comm
                %taskIndex = 1;
                SystemSet = Communication_Set; 
                idx = [1,2,3];
                for ii = idx
                    if ~isempty(SystemSet{ii})
                        nobj = length(SystemSet{ii});
                        nbat = nobj; 
                        is = it+1;
                        it = is+nobj+nbat-1; 
                        sample = samples(isample,is:it); 
                        taskIndex = Library.DamageEval(SystemSet{ii},EventType,IM,IMx,IMy,Dictionary,sumTaskHash, sumDamageTaskHash,taskIndex,Prob_Magic_Battery, sample);
                    end
                end
            end
            % No. of tasks to repair damaged objects in the Communication system
            Ntaskc = taskIndex-Ntaskp; 
            
            %---- Transportation_Set = {Road_Set, Bridge_Set, TrafficLight_Set,RoadNode_Set,Neighborhood_Trans_Set};
            % damagable components: 1-road, 2-bridge, 3-traffic light
            % Index of object set in the system set that can be damaged: 1, 2, 3
            if active_trans
                SystemSet = Transportation_Set; 

                %--- Road Set
                idx = 1;
                for ii = idx 
                    if ~isempty(SystemSet{ii})
                        nobj = length(SystemSet{ii});
                        is = it+1;
                        it = is+nobj-1; 
                        sample = samples(isample,is:it);
                        taskIndex = Library.DamageEval(SystemSet{ii},EventType,IM,IMx,IMy, Dictionary,sumTaskHash, sumDamageTaskHash,taskIndex,Prob_Magic_Battery, sample);
                    end
                end
                
                %--- Bridge Set
                % Considering sub-component analysis
                idx = 2;
                for ii = idx 
                    if ~isempty(SystemSet{ii})
                        nobj = length(SystemSet{ii});
                        nobjsub = 0;
                        for k = 1:nobj
                            Object = SystemSet{ii}{k};
                            if Object.HasSub % Object.HasSub=1, meaning that sub-component analysis is turned ON.  
                                SubObjSet = Object.ColumnSet;
                                nobjsub = nobjsub + length(SubObjSet);
                                SubObjSet = Object.ColumnFoundSet;
                                nobjsub = nobjsub + length(SubObjSet);
                                SubObjSet = Object.AbutmentSet;
                                nobjsub = nobjsub + length(SubObjSet);
                                SubObjSet = Object.AbutmentFoundSet;
                                nobjsub = nobjsub + length(SubObjSet);
                                SubObjSet = Object.GirderSet;
                                nobjsub = nobjsub + length(SubObjSet);
                                SubObjSet = Object.BearingSet;
                                nobjsub = nobjsub + length(SubObjSet);
                                SubObjSet = Object.DeckSet;
                                nobjsub = nobjsub + length(SubObjSet);
                                SubObjSet = Object.SlabSet;
                                nobjsub = nobjsub + length(SubObjSet);
                            end
                        end
                        is = it+1;
                        it = is+nobj+nobjsub-1; 
                        sample = samples(isample,is:it);
                        taskIndex = Library.DamageEval(SystemSet{ii},EventType,IM,IMx,IMy, Dictionary,sumTaskHash, sumDamageTaskHash,taskIndex,Prob_Magic_Battery, sample);
                    end
                end
                
                %--- TrafficLight
                % Index of object set in the system set that can be
                % damaged: 3-traffic light, with the possibility of a magic battery
                idx = 3; 
                for ii = idx
                    if ~isempty(SystemSet{ii})
                        nobj = length(SystemSet{ii});
                        nbat = nobj;
                        is = it+1;
                        it = is+nobj+nbat-1; 
                        sample = samples(isample,is:it);
                        taskIndex = Library.DamageEval(SystemSet{ii},EventType,IM,IMx,IMy, Dictionary,sumTaskHash, sumDamageTaskHash,taskIndex,Prob_Magic_Battery, sample);
                    end
                end
                
            end
            % No. of tasks to repair damaged objects in the Transportation system
            Ntaskt = taskIndex-Ntaskp-Ntaskc;
            
            %=== Nsystask: a vector of the number of tasks to repair damaged objects in the three system
            NSysTask = [Ntaskp-1, Ntaskc, Ntaskt];
            
            %=== Evaluate Stoped or Cascading Damaged Object (or Compoment)
            Library.StopedEval(active_power, active_comm, active_trans, Power_Set, Communication_Set, Transportation_Set, Dictionary, Neighborhood,sumTaskHash, sumDamageTaskHash,taskIndex);
            
            
            %=== Remove non-functional objects (nodes/links) from the
            %network graph for every system
            if active_power
                [powerGraph, NrmNodePower, NrmLinkPower] = Library.RemovePowerGraph(Power_Set, powerGraph);
            end
            if active_comm
                [commGraph, NrmNodeComm, NrmLinkComm] = Library.RemoveCommGraph(Communication_Set, commGraph);
            end
            if active_trans
                [transGraph, NrmLinkTrans] = Library.RemoveTransGraph(Transportation_Set, transGraph);
            end
            
        end % function DamageAndStop 
        
        %% Assess structural component damage based on fragility analysis (STEP 2)
        function taskIndex = DamageEval(cell,EventType,IM,IMx,IMy,Dictionary,sumTaskHash,sumDamageTaskHash,taskIndex,Prob_Magic_Battery,sample)
            % ===============================================================================
            % function DamageEval
            % This function uses fragility curves to determine 
            % the probabilty of exceedance (Prob_Failure) for every structural component 
            % at different damage states in different damage sceanrio samples. It then
            % computes possible damage states as Object.DamageLevel (=0,1,2,3,4)
            % and assigns restoration tasks based on the damage level. 
            % --- damage assessment --- 
            % Object.DamageLevel represents discrete damage states, based on fragility curves.
            % DamageLevel = 0: no damage
            % DamageLevel = 1: slight damage
            % DamageLevel = 2: moderate damage
            % DamageLevel = 3: extensive damage
            % DamageLevel = 4: complete damage
            % --- Output variables
            % (1) taskIndex
            % (2) DamageSamples
            % DamageSamples{1} = component damage samples in the power system 
            % DamageSamples{2} = component damage samples in the communication system
            % DamageSamples{3} = component damage samples in the transportation system
            % DamageSamples{1}{1} = component damage samples of Power_Set{1}
            % DamageSamples{1}{2} = component damage samples of Power_Set{2}
            % DamageSamples{1}{3} = component damage samples of Power_Set{3}
            % DamageSamples{1}{4} = component damage samples of Power_Set{4}
            % DamageSamples{2}{1} = component damage samples of Comm_Set{1}
            % DamageSamples{2}{2} = component damage samples of Comm_Set{2}
            % DamageSamples{2}{3} = component damage samples of Comm_Set{3}
            % DamageSamples{3}{1} = component damage samples of Trans_Set{1}
            % DamageSamples{3}{2} = component damage samples of Trans_Set{2}
            % DamageSamples{3}{3} = component damage samples of Trans_Set{3}
            % DamageSamples{1}{1}{1} = component damage level from the fragility analysis for Power_Set{1}
            % DamageSamples{1}{1}{2} = cascading failure (stop) due to dependencies for Power_Set{1}
            % DamageSamples{1}{1}{3} = component functionality value after
            % the analysis of cascading faulures for Power_Set{1} 
            % ===============================================================================  
            %=== Field
            %--- idv: Index of random variable as the column number in the matrix of "sample", 
            % representing the variable index of magic battery & the variable index of component damage states
            idv = 0;
            
            %=== Set Up Initial Status for Generators and the possible existence of Magic Battery for different objects
            %1. Assuming Generator never gets damage; therefore, there is no task needed for repairing Generator.
            %2. Magic Battery: Unexpected Functionality Support
            % (meaning some unexpected resource providing functionality, which is not physically modeled in this platform) 
            % 'Centraloffice','CommunicationTower','Cellline','TrafficLight': Contigency Plan models that unexpected functionality may happen.
            % Its probability check (comaring rand to Prob_Magic_Batter) determines whether there is any unexpected resource avaiable to relax the original functionality dependency on electricity.
            switch cell{1}.Class
                case {'Generator'}
                    for iobj = 1:size(cell,2)% index of object in the ObjectSet of "cell"
                        Prob_Failure = 0;
                        cell{iobj}.DamageLevel = 0; 
                        cell{iobj}.Functionality = 1; 
                        cell{iobj}.Status='Open';
                    end
                case {'Centraloffice','CommunicationTower','Cellline','TrafficLight'}
                    nobj = size(cell,2);
                    is = idv+1;
                    idv = is+nobj-1;
                    check = double(le(sample(1,is:idv), Prob_Magic_Battery)); 
                    for iobj = 1:nobj % index of object in the ObjectSet of "cell": in these cases, iobj = ibattery
                        cell{iobj}.Battery = check(1,iobj);
                    end
            end

            %=== Objects without Potential Unexpected Functionality 
            % The following objects do not consider any unexpected functionality. 
            %(meaning functionality that is not physically modeled in this platform
            switch cell{1}.Class
                %--- Objects with a point location
                case {'Bus', 'TransmissionTower', 'Centraloffice', 'CommunicationTower','TrafficLight'}
                %case {'TransmissionTower'} 
                    for ii = 1:size(cell,2) % index of object in the ObjectSet of "cell"
                        %--- determine the IM at the location of an object interpolate the event intensity at the location of an object
                        % IMx = longitude;
                        % IMy = latitude; 
                        im = IM{cell{ii}.FragilityIM{1}};
                        Intensity = interp2(IMx,IMy,im,cell{ii}.Location(2),cell{ii}.Location(1));
                        cell{ii}.FragilityIM{3} = Intensity;
                        
                        %--- damage assessment with fragilitu curves
                        % assign the probabilty of exceedance for every structural component at different damage state 
                        Prob_Failure = Library.Prob_Failure(Intensity,cell{ii});
                        
                        %--- get the random sample of the probabality 
                        idv = idv+1;
                        aa = sample(1,idv);
                        
                        %--- compare the Prob_Failure with the random sample to determine the value of Object.DamageLevel
                        % index: the value of Object.DamageLevel
                        bb = aa - [1, Prob_Failure, 0]; 
                        index = find(diff(sign(bb)));
                        dl = min(index)-1; 
                        dl(dl<0) = 0;
                        dl(dl>4) = 4;
                        cell{ii}.DamageLevel = dl; 
 
                        %--- assign tasks based on the Object.DamageLevel
                        % Assuming generators never damage, therefore, there is no task for repairing generator. 
                        % For other types of object, assign tasks based on their damage levels
                        if cell{ii}.DamageLevel > 0
                            cell{ii}.Status='Damaged';
                            cell{ii}.Functionality = 0.0;
                            [cell,taskIndex] = Library.assignTask(cell, ii, sumTaskHash, sumDamageTaskHash, Dictionary, taskIndex);
                        end
                        
                    end
   
                %--- Objects with two point(end) locations
                % These objects are in the form of links, defined by two end points. 
                % How to compute their fragility is a problem. [which intensity value shall we use?]
                % The following lines use the event intensity at the
                % start/end of a link to compute fragility, which may need to be updated.  
                case {'Branch', 'Cellline', 'Road'}
                    for i = 1:size(cell,2)
                        %=== determine the IM at the location of an object interpolate the event intensity at the location of an object
                        % IMx = longitude;
                        % IMy = latitude; 
                        im = IM{cell{i}.FragilityIM{1}};
                        Intensity=interp2(IMx,IMy,im,cell{i}.Start_Location(2),cell{i}.Start_Location(1));
                        % Intensity=0.5*( Intensity + interp2(IMx,IMy,IM,cell{i}.End_Location(1),cell{i}.End_Location(2)) );
                        cell{i}.FragilityIM{3} = Intensity;          
                             
                        %=== damage assessment
                        % assign the probabilty of exceedance for every structural component at different damage state
                        Prob_Failure = Library.Prob_Failure(Intensity,cell{i});
                        %Y = [1, Prob_Failure, 0];
                        
                        %--- get the random sample of the probabality 
                        idv = idv+1;
                        aa = sample(1,idv);
                        
                        %--- compare the Prob_Failure with the randoms
                        % sample to determine the value of Object.DamageLevel
                        % index: the value of Object.DamageLevel
                        bb = aa - [1, Prob_Failure, 0]; 
                        index = find(diff(sign(bb)));
                        dl = min(index)-1; 
                        dl(dl<0) = 0;
                        dl(dl>4) = 4;
                        cell{i}.DamageLevel = dl;
                        
                        %--- assign tasks based on the Object.DamageLevel
                        if cell{i}.DamageLevel > 0
                            cell{i}.Status = 'Damaged';
                            cell{i}.Functionality = 0;
                            [cell,taskIndex] = Library.assignTask(cell,i,sumTaskHash, sumDamageTaskHash, Dictionary,taskIndex);
                        end
                        
                        if isempty(cell{i}.DamageLevel)
                            cell{i}.DamageLevel = 0;
                            disp(strcat('Function DamageEval Warning: empty DamageLevel for ', cell{i}.Class, num2str(i),', now DamageLevel is assigned as 0.'));
                        end
                        
                    end
                    
                %%--- Bridge Objects with a point location
                % if the object is very important, we want to use component-level anslyses of fragility and restoration. 
                % For example, a bridge object with subcomponents.
                case {'Bridge'}
                    for i = 1:size(cell,2)
                        %---- for Bridge Object with subcomponent analysis
                        if cell{i}.HasSub == 1 % sub-component analysis is turned ON.  
                            newcells = [cell{i}.ColumnSet,cell{i}.ColumnFoundSet,cell{i}.AbutmentSet,cell{i}.AbutmentFoundSet,cell{i}.GirderSet,cell{i}.BearingSet,cell{i}.SlabSet,cell{i}.DeckSet]; 
                            newcells = num2cell(newcells); 
                            % interpolate the event intensity at the location of the bridge object
                            %IMx = longtitude, IMy = latitude;
                            %Object.Location = [latitude, longitude]
                            % IMidx: the index of PGA in IMmeta in the earthquake hazard scenario map
                            IMidx = 8; 
                            IMstr = 'PGA(g)';
                            im = IM{IMidx};

                            for j = 1:length(newcells)
                                Intensity = interp2(IMx,IMy,im,newcells{j}.Location(2),newcells{j}.Location(1));
                                newcells{j}.FragilityIM{1} = IMidx; 
                                newcells{j}.FragilityIM{2} = IMstr; 
                                newcells{j}.FragilityIM{3} = Intensity; 
                                
                                % compute the probabilty of exceedance based on fragility curves
                                Prob_Failure = Library.Prob_Failure(Intensity,newcells{j});
                                
                                %--- get the random sample of the probability 
                                idv = idv+1;
                                aa = sample(1,idv);
                                bb = aa - [1, Prob_Failure, 0]; 
                                index = find(diff(sign(bb)));
                                dl = min(index)-1; 
                                dl(dl<0) = 0;
                                dl(dl>4) = 4;

                                newcells{j}.DamageLevel = dl;
                                if isempty(newcells{j}.DamageLevel)
                                    newcells{j}.DamageLevel = 0;
                                end
                                dlcheck(j) = newcells{j}.DamageLevel; 
                                % Assign restoration tasks to the damaged sub-component
                                if newcells{j}.DamageLevel > 0
                                    newcells{j}.Status = 'Damaged';
                                    [newcells,taskIndex] = Library.assignTask(newcells,j,sumTaskHash, sumDamageTaskHash, Dictionary,taskIndex);
                                end   
                                
                            end % j~1:length(newcells)
                            if any(dlcheck) 
                                cell{i}.Status = 'Damaged';
                            end
                        
                        %----%cell{i}.HasSub == 0: for Bridge without the subcomponent analysis
                        % IMidx: the index of PGA in IMmeta in the earthquake hazard scenario map:
                        % IMidx = 2: Sa(0.3second)
                        % IMidx = 3: Sa(1.0second) Spectrum Acceleration at 1.0 second(g)
                        % IMidx = 10: Peak Ground Displacement (PGD) (inch)
                        else  %cell{i}.HasSub == 0
                            idv = idv+1;
                            aa = sample(1,idv);
                            
                            if EventType == 1 % earthquake
                                IMidx = [2 3 10]; 
                                for k = 1:length(IMidx)
                                    kk = IMidx(k);
                                    im = IM{kk}; 
                                    IntensityMeasureVector(k) = interp2(IMx,IMy,im,cell{i}.Location(2),cell{i}.Location(1));
                                end

                                Object = cell{i};
                                SoilClass = [];
                                Prob_Failure = Library.ProbabilityFailureBridge(IntensityMeasureVector,Object,SoilClass);
								
                            elseif EventType == 2 % hurricane
                                im = IM{cell{i}.FragilityIM{1}};
                                Intensity=interp2(IMx,IMy,im,cell{i}.Location(2),cell{i}.Location(1));
                                % Intensity=0.5*( Intensity + interp2(IMx,IMy,IM,cell{i}.End_Location(1),cell{i}.End_Location(2)) );
                                cell{i}.FragilityIM{3} = Intensity;          

                                %=== damage assessment
                                % assign the probabilty of exceedance for every structural component at different damage state
                                Prob_Failure = Library.Prob_Failure(Intensity,cell{i});
                            end
                            
                            
                            
                            for k = 1:size(Prob_Failure,1)
                                p = Prob_Failure(k,:);
                                bb = aa - [1, p, 0]; 
                                
                                index = find(diff(sign(bb)));
                                dl = min(index)-1; 
                                dl(dl<0) = 0;
                                dl(dl>4) = 4;
                                cell{i}.DamageLevel = dl;
                            end
                            
                                
                            if isempty(cell{i}.DamageLevel)
                                cell{i}.DamageLevel = 0;
                                disp(strcat('Function DamageEval Warning: empty DamageLevel for Bridge',num2str(i),', now is assigned as 0.'));
                            end
                            
                            if cell{i}.DamageLevel > 0
                                cell{i}.Status='Damaged';
                                cell{i}.Functionality = 0;
                                [cell,taskIndex] = Library.assignTask(cell,i,sumTaskHash, sumDamageTaskHash, Dictionary,taskIndex);
                                
                            end
                        end
                        
                    end
    
            end

            
        end
        
        %% Damage Assessment for every Object with fragility curves (STEP 2)
        function Prob_Exceedance = Prob_Failure(Intensity,Object)
            % ================================================================
            % function Prob_Failure
            % An rough estimate of the exceeding probability at different damage levels
            % Prob_Exceedance = Exceedance Probability of the damage state [slight, moderate, extensive, complete]
            %--- Description
            % 1. Assuming the cumulative probability density function (CDF) as a lognormal distribution. 
            % 2. The lognormal distribution is defined by mu and sigma.
            % 3. mu and sigma is related to the mean and std 
            % (standdard deviation) of the lognormal distribution.
            % ================================================================
            
            %=== Initial probablity values as zeros
            Prob_Exceedance = zeros(1,4); 
            
            %=== Performing fragility analysis with every fragility curve
            % iIM: the row index of matrix in the matrix of
            % Object.Fragility.  Every row represents 8 values of mu and sigma 
            % for a fragility curve. 
            
            %--- if there are multiple fragility curves corresponding to 
            % different IMs, then Object.Fragility is a cell
            % else, Object.Fragility is a 4x2 matrix. 
            if ~iscell(Object.Fragility)   % only one fragility curve 
                for ids = 1:4 % 4 damage state (ds)
                    clearvars mu sigma;
                    mean = Object.Fragility(ids,1);
                    std = Object.Fragility(ids,2);
                    mu = log(mean) - 0.5*log(1+ (std^2)/(mean^2));
                    sigma = sqrt(log(1+ (std^2)/(mean^2)));
                    
                    %Exceedance Probability
                    Pex(ids) =  cdf('lognormal',Intensity,mu,sigma);
                end
                Prob_Exceedance = Pex; 
                
            else  % multiple fragility curves using different IMs
                for iIM = 1:size(Object.Fragility,1)
                    for ids = 1:4 % 4 damage state (ds)
                        clearvars mu sigma;
                        mean = Object.Fragility{iIM}(ids,1);
                        std = Object.Fragility{iIM}(ids,2);
                        mu = log(mean) - 0.5*log(1+ (std^2)/(mean^2));
                        sigma = sqrt(log(1+ (std^2)/(mean^2)));
                        
                        %Exceedance Probability
                        Pex{iIM}(ids) =  cdf('lognormal',Intensity,mu,sigma);
                    end
                end
                
                %--- Compute the final value of Probability of exeedanace 
                % as the first cell of P_ex determined from different 
                % fragility curves for the same structural component. 
                iIM = 1;
                Prob_Exceedance = Pex{iIM};
            end
            
            %=== postprocess 
            % to ensure computation correct in the damage assessment
            
            %--- (1)remove values outside of [0,1] in the 4-element vector of Prob_Exceedance
            Prob_Exceedance(Prob_Exceedance<0) = 0;
            Prob_Exceedance(Prob_Exceedance>1) = 1;
            
            %--- (2)sort the 4-element vector of Prob_Exceedance in the
            %descending order, which means that the probabality of exceedance
            %decreases as from the slight damage state to complete damage
            %state (from DamageLevel=1 to DamageLevel=4)
            Prob_Exceedance = sort(Prob_Exceedance,'descend');
            
        end
        
		%% Function ProbabilityFailureBridge: compute bridge damage probability based on fragility curves 
        function Prob_Damage = ProbabilityFailureBridge(IntensityMeasureVector,Object, SoilClass)
            %==== function ProbabilityFailureBridge(IntensityMeasureVector,Object)
            % refer to HAZUS-EQ technical manual (hazard type: 1-earthquake)
            % Table 7.7 on page 7-12 
            % By following steps 1-9 on page 7-12~7-16 
			% https://www.fema.gov/media-library-data/20130726-1820-25045-6286/hzmh2_1_eq_tm.pdf
            %---------input--------
            % IntensityMeasureVector = [Sa(0.3s)  Sa(1s) PGD];
            % Object = individual bridge object 
            % PGA: peak ground acceleration(g)
            % PGD: permanent ground deformation (inch)
            % PGV: peak ground velocity (inch/second)
            % Sa03: spectral acceleration [0.3 s] (g)
            % Sa10: spectral acceleration [1.0 s] (g)
            %----------------------
            % Object = bridge
			% hazard type: 1-earthquake
            % without considering soil amplification factor 
            % (table 4.10 on Page 4-16 in HAZUS-EQ technical manual)
			% https://www.fema.gov/media-library-data/20130726-1820-25045-6286/hzmh2_1_eq_tm.pdf
			%----------------------
            % SiteClass: 
            % The National Earthquake Hazards Reduction Program (NEHRP) has 
            % defined 5 soil types (A-E) based on their shear-wave velocity (Vs).
            % 'A': Hard rock, Vs(shear-wave velocity) > 1500 m/sec
            % 'B': Rock, Vs(shear-wave velocity) 750~1500 m/sec
            % 'C': Soft rock and very dense soil, Vs(shear-wave velocity) 350~750 m/sec
            % 'D': Stiff soil, Vs(shear-wave velocity) 200~350 m/sec
            % 'E': Soft clay soil, Vs(shear-wave velocity) <200 m/sec
            % 'F': Site response analysis
            %https://earthquake.usgs.gov/hazards/urban/sfbay/soiltype/
            
            if isempty(SoilClass)
                SoilClass = 'D'; % as default
            end
            
            [FA, FV] = Library.getSoilAmplificationFactor(IntensityMeasureVector, SoilClass);
            
            Sa03 = FA*IntensityMeasureVector(1);
            Sa10 = FV*IntensityMeasureVector(2);
            PGD = IntensityMeasureVector(3);

			% Inferring PGV from Sa(1.0s) based on
            % equation (4-5) in technical manual of HAZUS-earthquake, pp 4-9.  
            % PGV = (386.4*Sa10)*inv(2*pi*1.65);  
            % https://www.fema.gov/media-library-data/20130726-1820-25045-6286/hzmh2_1_eq_tm.pdf
			
            idx = 1; 
            fragility_Sa10 = Object.Fragility{2};
            
            % bridge information
            alpha = Object.SkewAngle*pi*inv(180);
       
            % modification factors fo Sa10 and PGD
            % The Kshape factor is the modifier that converts cases for
            % short periods to an equivalent spectral amplitude at T=1.0 second
            Kskew = sqrt(sin(0.5*pi-alpha));
            Kshape = 2.5*Sa10*inv(Sa03);
            
            [K3D,Ishape,f1,f2] = Library.getBridgeDamageModificationFactor(Object); 
            

            % Ishape: a Boolean indicator
            %When Ishape = 0, the Kshape factor does not apply. When Ishape
            %= 1, the Kshape factor applies.
            if Ishape == 0
                Factor_slight = 1;
            elseif Ishape == 1
                Factor_slight = min(1,Kshape);
            end
            
            Factor_other = Kskew*K3D;
            
            %% Sa: damage probability due to Ground Shaking

            Sa10_NewMedian = [Factor_slight;repmat(Factor_other,3,1)].*fragility_Sa10(:,1); 
            
            % Postprocess Sa10_NewMedian to ensure that all median values
            % are sorted. If not, make the one that violate the sorting
            % rank as a mean value of its two adjacent values
            if ~issorted(Sa10_NewMedian)
                a = [0; Sa10_NewMedian; max(Sa10_NewMedian)];
                b = find(a(1:end-1)>a(2:end));
                Sa10_NewMedian(b) = mean([a(b),a(b+2)]);
            end
                    
            median_Sa10 = Sa10_NewMedian;
            sigma_Sa10 = fragility_Sa10(:,2); 
            mu_Sa10 = log(median_Sa10);     
            p = logncdf(Sa10,mu_Sa10,sigma_Sa10)';
            P_Sa10 = [1, p] - [p,0]; 
            
            %% PGD: damage probability due to Ground Failure 
            fragility_PGD = Object.Fragility{1};
            %Beta_PGD = 0.2; 
            median_PGD = [repmat(f1,3,1); f2].*fragility_PGD(:,1); 
            
            % Postprocess median_PGD to ensure that all median values
            % are sorted. If not, make the one that violate the sorting
            % rank as a mean value of its two adjacent values
            if ~issorted(median_PGD)
                a = [0; median_PGD; max(median_PGD)];
                b = find(a(1:end-1)>a(2:end));
                median_PGD(b) = mean([a(b),a(b+2)]);
            end
            
            sigma_PGD = fragility_PGD(:,2);
            mu_PGD = log(median_PGD);     
            p = logncdf(PGD,mu_PGD,sigma_PGD)';
            P_PGD = [1, p] - [p,0];
            
            Prob_Damage(1,:) = P_Sa10;
            Prob_Damage(2,:) = P_PGD;

        end
        
		%% Function getSoilAmplificationFactor
        function [FA, FV] = getSoilAmplificationFactor(IntensityMeasureVector, SiteClass)
           % function  getSoilAmplificationFactor(IntensityMeasureVector, SiteClass)
           % refer to Table 4.9 on Page 4-16 in HAZUS-EQ technical manual
           Sa03 = IntensityMeasureVector(1);
           Sa10 = IntensityMeasureVector(2);
           
           switch SiteClass
               case 'A'
                   FA = 0.8;
                   FV = 0.8; 
               case 'B'
                   FA = 1.0;
                   FV = 1.0; 
               case 'C'   
                   if le(Sa03,0.5)
                       FA = 1.2;   
                   elseif gt(Sa03,1)
                       FA = 1; 
                   else 
                       FA = 1.1; 
                   end
                   
                   if le(Sa10,0.1)
                       FV = 1.7; 
                   elseif gt(Sa10,0.1)&&le(Sa10,0.2)
                       FV = 1.6;
                   elseif gt(Sa10,0.2)&&le(Sa10,0.3)
                       FV = 1.5;    
                   elseif gt(Sa10,0.3)&&le(Sa10,0.4)
                       FV = 1.4;
                   else 
                       FV = 3; 
                   end     
                   
               case 'D'   
                   if le(Sa03,0.25)
                       FA = 1.6; 
                   elseif gt(Sa03,0.25)&&le(Sa03,0.5)
                       FA = 1.4;
                   elseif gt(Sa03,0.5)&&le(Sa03,0.75)
                       FA = 1.2;    
                   elseif gt(Sa03,0.75)&&le(Sa03,1)
                       FA = 1.1;
                   else 
                       FA = 1; 
                   end   
                   
                   if le(Sa10,0.1)
                       FV = 2.4; 
                   elseif gt(Sa10,0.1)&&le(Sa10,0.2)
                       FV = 2.0;
                   elseif gt(Sa10,0.2)&&le(Sa10,0.3)
                       FV = 1.8;    
                   elseif gt(Sa10,0.3)&&le(Sa10,0.4)
                       FV = 1.6;
                   else 
                       FV = 1.5; 
                   end  
                   
               case 'E'   
                   if le(Sa03,0.25)
                       FA = 2.5; 
                   elseif gt(Sa03,0.25)&&le(Sa03,0.5)
                       FA = 1.7;
                   elseif gt(Sa03,0.5)&&le(Sa03,0.75)
                       FA = 1.2;    
                   else 
                       FA = 0.9; 
                   end   
                   
                   if le(Sa10,0.1)
                       FV = 3.5; 
                   elseif gt(Sa10,0.1)&&le(Sa10,0.2)
                       FV = 3.2;
                   elseif gt(Sa10,0.2)&&le(Sa10,0.3)
                       FV = 2.8;    
                   else 
                       FV = 2.6; 
                   end  
               
           end
           
        end
        
        %% Function getBridgeDamageModificationFactor
        function [K3D,Ishape,f1,f2] = getBridgeDamageModificationFactor(Object)
            % function getBridgeDamageModificationFactor(Object)
            % refer to HAZUS-EQ technical manual
            % Table 7.2 and Table 7.3 on page 7-6 ~7-8 
			% https://www.fema.gov/media-library-data/20130726-1820-25045-6286/hzmh2_1_eq_tm.pdf
			
            %% bridge information
            N = Object.MainSpans; %No of Span
            alpha = Object.SkewAngle*pi*inv(180);
            W = Object.Width;
            L = Object.Length;
            Lmax = Object.MaxSpanLength;
            c = Object.Type; 
                        
            %% compute K3D = 1+A/(N-B)
			% When c = 'HWB1', 'HWB2' 'HWB3', 'HWB4', or 'HWB28'
            A = 0.25; B = 1;
            Ishape = 0;
            f1 = 1; 
            f2 = f1;
			
            if strcmp(c, 'HWB5')
                A = 0.25; B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1;    
            elseif strcmp(c, 'HWB6')
                A = 0.25; B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1;     
            elseif strcmp(c, 'HWB7')
                A = 0.25; B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB8')
                A = 0.33; B = 0;
                Ishape = 0;
                f1 = 1;
                f2 = sin(alpha);  
            elseif strcmp(c, 'HWB9')
                A = 0.33;
                B = 1;
                Ishape = 0;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB10')
                A = 0.33;
                B = 0;    
                Ishape = 1;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB11')
                A = 0.33;
                B = 1;   
                Ishape = 1;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB12')
                A = 0.09;
                B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB13')
                A = 0.09;
                B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB14')
                A = 0.25;
                B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB15')
                A = 0.05;
                B = 0;    
                Ishape = 1;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB16')
                A = 0.33;
                B = 1;   
                Ishape = 1;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB17')
                A = 0.25;
                B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB18')
                A = 0.25;
                B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB19')
                A = 0.25;
                B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB20')
                A = 0.33;
                B = 0;  
                Ishape = 0;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB21')
                A = 0.33;
                B = 1;  
                Ishape = 0;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB22')
                A = 0.33;
                B = 0;
                Ishape = 1;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB23')
                A = 0.33;
                B = 1;
                Ishape = 1;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB24')
                A = 0.20;
                B = 1;
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB25')
                A = 0.20;
                B = 1;    
                Ishape = 0;
                if alpha == 0
                    f1 = 1;
                else
                    f1 = 0.5*L*inv(N*W*sin(alpha)); 
                end
                f2 = f1; 
            elseif strcmp(c, 'HWB26')
                A = 0.20;
                B = 1;   
                Ishape = 1;
                f1 = 1;
                f2 = sin(alpha);
            elseif strcmp(c, 'HWB27')
                A = 0.10;
                B = 0;
                Ishape = 0;
                f1 = 1;
                f2 = sin(alpha);
            end 
            
            %K3D: a factor that modifies the piers’ 2-dimensional capacity
            % to allow for the 3-dimensional arch action in the deck    
            K3D = 1+A*inv(N-B);    
            
            
        end
        

        %% Evaluate cascading failures (Object.Status = Stoped) after initial damage assessment (STEP 2)
        function StopedEval(active_power, active_comm, active_trans, Pow, Comm, Trans, Dictionary,Neighborhood,sumTaskHash, sumDamageTaskHash,taskIndex)
            % ================================================================
            % function StopedEval
            % This function assess the "Stoped" status for every object
            % after the initial damage assessment with fragility curves. 
            % ================================================================
            %==== Field
            Branch= Pow{1};
            Bus= Pow{2};
            Generator= Pow{3};
            TransTower = Pow{4};
            
            Centraloffice = Comm{1};
            CommTower = Comm{2};
            Cellline = Comm{3};

            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};
            
            Branch_list = [];
            Bus_list = [];
            
            TrafficLight_list = [];
            Bridge_list = [];
            Road_list = [];

            Cellline_list = [];
            Centraloffice_list = [];
            CommTower_list = [];
            Neighborhood_list = [];
            
            
            %% Power System
            
            if active_power
                %===== Bus
                % If a Bus is damaged, the connected Branches "stoped" due
                % to functionality dependency.
                for i = 1:length(Bus)
                    if strcmp(Bus{i}.Status, 'Damaged') 
                        for j = 1:length(Bus{i}.Branch)
                            if ~strcmp(Branch{Bus{i}.Branch(j)}.Status, 'Damaged')
                                Branch{Bus{i}.Branch(j)}.Status = 'Stoped';
                            end
                        end
                    end
                    
                end
                
                %==== Transmission Tower
                % If a TransTower is completely damaged (damage level = 4),
                % the connected Branches mechanically break due to mechanical functionality dependency. 
                % If a TransTower is not completely damaged (damage level < 4),
                % the connected Branches have the damage state of "stoped"
                % due to functionality dependency.
                for i = 1:length(TransTower)
                    if strcmp(TransTower{i}.Status, 'Damaged') 
                        if TransTower{i}.DamageLevel == 4
                            for j = 1:length(TransTower{i}.Branch)
                                if ~strcmp(Branch{TransTower{i}.Branch(j)}.Status, 'Damaged')
                                    Branch{TransTower{i}.Branch(j)}.Status = 'Damaged';
                                    Branch{TransTower{i}.Branch(j)}.Functionality = 0;
                                    Branch{TransTower{i}.Branch(j)}.DamageLevel = 4;
                                    [Branch,taskIndex] = Library.assignTask(Branch,TransTower{i}.Branch(j),sumTaskHash, sumDamageTaskHash, Dictionary,taskIndex);
                                end
                                
                            end
                        else
                            for j = 1:length(TransTower{i}.Branch)
                                if ~strcmp(Branch{TransTower{i}.Branch(j)}.Status, 'Damaged')
                                    Branch{TransTower{i}.Branch(j)}.Status = 'Stoped';
                                    Branch{TransTower{i}.Branch(j)}.Functionality = 0;
                                end
                            end
                        end
                    end
                    
                end
                
                %==== Generator connected to a Bus (Generator implemented as a property of Bus)
                % If a Generator is damaged AND other connected objects (Bus and Branch) are structurally safe,
                % the Bus connected to this Generator has the damagae state of "stoped", 
                % and its bus.branch has the damage state of "stoped"  due to functionality dependency.
                for i = 1:length(Bus)
                    if ~isempty(Bus{i}.Generator) && strcmp(Bus{i}.Status, 'Open') 
                        bus = extractAfter(Bus{i}.Generator, 9);
                        bus = str2num(bus);
                        if strcmp(Generator{bus}.Status, 'Damaged')
                            if ~strcmp(Bus{i}.Status, 'Damaged')
                                Bus{i}.Status = 'Stoped';
                                for j = 1:length(Bus{i}.Branch)
                                    if ~strcmp(Branch{Bus{i}.Branch(j)}.Status, 'Damaged')
                                        Branch{Bus{i}.Branch(j)}.Status = 'Stoped';
                                        Branch{Bus{i}.Branch(j)}.Functionality = 0;
                                    end
                                end
                            end
                        end
                    end
                end
                
                %==== Branch
                % If a Branch is not damaged AND any of its connected Bus
                % is damaged, this Branch has the damage state of "stoped"
                % due to functionality dependency.
                for i = 1:length(Branch)
                    obj1 = Dictionary(Branch{i}.connectedObj1);
                    obj2 = Dictionary(Branch{i}.connectedObj2);
                    if ~strcmp(obj1{1}.Status, 'Open') || ~strcmp(obj2{1}.Status, 'Open')
                        if ~strcmp(Branch{i}.Status, 'Damaged')
                            Branch{i}.Status = 'Stoped';
                            Branch{i}.Functionality = 0;
                        end
                    end
                end
                for i = 1:length(Branch)
                    if ~strcmp(Branch{i}.Status, 'Open')
                        Branch_list = [Branch_list, i];
                    end
                end
                
                %% Neighborhood-PowerFunctionality, connected to Bus (by default the nearest Bus UNLESS other specified)
                % if the Bus that connects to the Neighborhood is not "Open", 
                % this Neighborhood object has the status of "stoped" due to 
                % functionality dependency, as its power functionality as 0. 
                for i = 1:length(Bus)
                    if ~strcmp(Bus{i}.Status, 'Open')% "Bus_List": any Bus that is not "Open".  
                        for j = 1:length(Bus{i}.Neighborhood)
                            temp = Dictionary(Bus{i}.Neighborhood{j});
                            temp = temp{1};
                            temp.PowerStatus = 0;
                            temp = Dictionary(Bus{i}.Neighborhood_Power_Link{j});
                            temp = temp{1};
                            %temp.Status = 'Stoped';
                        end
                    end
                end
            end
            
            
            %% Transportation System
            %==== TrafficLight
            % If a TrafficLight is not damaged AND its connected Bus is
            % "stoped", this TrafficLight is "stoped" due to functionality dependency.
            if active_trans && active_power
                for ii = 1:length(TrafficLight)
                    if strcmp(TrafficLight{ii}.Status, 'Damaged')
                        TrafficLight_list = [TrafficLight_list, ii];
                    else
                        bus = extractAfter(TrafficLight{ii}.Bus, 3);
                        bus = str2num(bus);
                        if ismember(bus, Bus_list)
                            TrafficLight{ii}.Status = 'Stoped';
                            TrafficLight{ii}.Functionality = 0;
                            TrafficLight_list = [TrafficLight_list, ii];
                        end
                    end
                end

            end
                
            %==== Road Segment
            % If a Road is not damaged AND the bridge that carries OR
            % crosses this Road segment is damaged, this Raad is "stoped" 
            % due to functionality dependency.
            if active_trans
                for i = 1:length(Road)
                    if strcmp(Road{i}.Status, 'Damaged')
                        Road_list = [Road_list, i];
                    else
                        flag = 0;
                        
                        for j = 1:length(Road{i}.Bridge_Carr)
                            if strcmp(Bridge{Road{i}.Bridge_Carr(j)}.Status, 'Damaged')
                                flag = flag + 1;
                            end
                        end
                        
                        for j = 1:length(Road{i}.Bridge_Cross)
                            if strcmp(Bridge{Road{i}.Bridge_Cross(j)}.Status, 'Damaged')
                                flag = flag + 1;
                            end
                        end
                        
                        if flag ~= 0
                            Road{i}.Status = 'Stoped';
                            Road{i}.Functionality = 0;
                            Road_list = [Road_list, i];
                        end
                    end
                end
                
                %% Neighborhood-TransportationFunctionality (connected to RoadNode by default) 
                % A Neighborhood should have three functionality
                % properties: PowerFunc, CommFunc, TransFunc. The current
                % implementation seems to have some problem to reflect this
                % idea. Therefore, the neighborhood-related codes need updates!!!
                % 
                % If a Road that connected to a Neighborhood virtual object 
                % is not open, this Neighborhood has its TransStatus as "stoped" 
                for i = 1:length(Road)
                    if ~strcmp(Road{i}.Status, 'Open')
                        temp = Dictionary(strcat('RoadNode',num2str(Road{i}.Start_Node)));
                        temp = temp{1};
                        for j = 1:length(temp.Neighborhood)
                            t1 = Dictionary(temp.Neighborhood{j});
                            t1 = t1{1};
                            t1.TransStatus = 0;
                            t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                            t1 = t1{1};
                            t1.Status = 'Stoped';
                        end
                        temp = Dictionary(strcat('RoadNode',num2str(Road{i}.End_Node)));
                        temp = temp{1};
                        for j = 1:length(temp.Neighborhood)
                            t1 = Dictionary(temp.Neighborhood{j});
                            t1 = t1{1};
                            t1.TransStatus = 0;
                            t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                            t1 = t1{1};
                            t1.Status = 'Stoped';
                        end
                    end
                end
            end
            
            %% Communication System
            if active_comm && active_power        
                %==== Communciation Line
                % A Cellline is not damaged AND its connected Bus is not functional (therefore this Bus is in the Bus_list),
                % this Cellline is "stoped" due to functionality dependency. 
                for i = 1:length(Cellline)
                    if strcmp(Cellline{i}.Status, 'Damaged')
                        Cellline_list = [Cellline_list, i];
                    else
                        bus = extractAfter(Cellline{i}.Bus, 3);
                        bus = str2num(bus);
                        if ismember(bus, Bus_list)
                            Cellline{i}.Status = 'Stoped';
                            Cellline{i}.Functionality = 0;
                            Cellline_list = [Cellline_list, i];
                        end
                    end
                end
                
                %==== Communciation Tower
                % If a CommTower is not damaged AND its connected Bus is not functional (therefore this Bus is in the Bus_list),
                % this CommTower is "stoped" due to functionality dependency. 
                for ii = 1:length(CommTower)
                    if strcmp(CommTower{ii}.Status, 'Damaged')
                        CommTower_list = [CommTower_list, ii];
                    else
                        bus = extractAfter(CommTower{ii}.Bus, 3);
                        bus = str2num(bus);
                        if ismember(bus, Bus_list)
                            CommTower{ii}.Status = 'Stoped';
                            CommTower{ii}.Functionality = 0;
                        end
                    end
                end
                
                %==== Central Office
                % If a Centraloffice is damaged, its (connected) Router is "damaged".
                % Otherwise,
                % If a Centraloffice is not damaged AND its connected Bus is not functional (therefore this Bus is in the Bus_list),
                % this Centraloffice is "stoped" due to functionality dependency. 
                for i = 1:length(Centraloffice)
                    if strcmp(Centraloffice{i}.Status, 'Damaged')
                        Centraloffice_list = [Centraloffice_list, i];
                    else
                        bus = extractAfter(Centraloffice{i}.Bus, 3);
                        bus = str2num(bus);
                        if ismember(bus, Bus_list)
                            Centraloffice{i}.Status = 'Stoped';
                            Centraloffice{i}.Functionality = 0;
                            Centraloffice_list = [Centraloffice_list, i];
                        end
                    end
                end               
            end
                
            %==== Neighborhood-CommunicationFunctionality, connected to Centraloffice (by default the nearest Centraloffice UNLESS other specified)
            for i = 1:length(Centraloffice)
                  if ~strcmp(Centraloffice{i}.Status, 'Open')
                       for j = 1:length(Centraloffice{i}.Neighborhood)
                           temp = Dictionary(Centraloffice{i}.Neighborhood{j});
                           temp = temp{1};
                           temp.CommStatus = 0;
                           temp = Dictionary(Centraloffice{i}.Neighborhood_Comm_Link{j});
                           temp = temp{1};
                           %temp.Status = 'Stoped';
                       end
                  end
            end
      
        end
        
        %% Transportation: Remove damaged nodes and edges after the damage assesment (STEP 2)
        function [G, Nrmlink] = RemoveTransGraph(Trans_Set, G)
            % ================================================================
            % function RemoveTransGraph
            % Remove Non-Open edges and the nodes of these Non-Open edges from 
            % the Transpotation network after the initial damage assesment
            % -----------------------------------------------------------------
            % node = RoadNode
            % edge = Road
            % s = start node of a Road segment
            % t = end node of a Road segment
            %  ================================================================

            hash = containers.Map('KeyType','double','ValueType','char');
            road_Set = Trans_Set{1};
            roadnode_Set = Trans_Set{4};
            for ii = 1:length(roadnode_Set)
                hash(roadnode_Set{ii}.nodeID) = roadnode_Set{ii}.uniqueID;
            end
            s = [];
            t = [];
            index = 1;
            for ii = 1:length(road_Set)
                if ~strcmp(road_Set{ii}.Status, 'Open')
                    s{index} = hash(road_Set{ii}.Start_Node);
                    t{index} =  hash(road_Set{ii}.End_Node);
                    index = index + 1;
                end
            end
            G = rmedge(G,s,t);
            Nrmlink = index;
        end

        %% Power: Remove damaged nodes and edges after the damage assesment (STEP 2)
        % Note: assuming that Generator (Power_Set{3}) do not get damaged!
        function [G, Nrmnode, Nrmlink] = RemovePowerGraph(Power_Set, G)
            % ================================================================
            % function RemovePowerGraph
            % Remove Non-Open edges and the nodes of these Non-Open edges from 
            % the Power network after the initial damage assesment
            % -----------------------------------------------------------------
            % node = Bus/TransTower
            % edge = Branch
            % s = start node of a Branch
            % t = end node of a Branch
            % s2 = unique Branch start nodes in the orginal power graph
            % t2 = unique Branch end nodes in the orginal power graph
            % ================================================================
            Branch = Power_Set{1};
            Bus = Power_Set{2};
            TransTower = Power_Set{4};

            s = [];
            t = [];
            st = [];
            tt = [];
            Lst = [];
            s2 = [];
            t2 = [];
            indexNode = 0;
            indexLink = 1;
            
            for i = 1:length(Bus)
                if ~strcmp(Bus{i}.Status, 'Open')
                    G = rmnode(G,Bus{i}.uniqueID);
                    indexNode = indexNode + 1;
                end
            end
            
            for i = 1:length(TransTower)
                if ~strcmp(TransTower{i}.Status, 'Open')
                    G = rmnode(G,TransTower{i}.uniqueID);
                    indexNode = indexNode + 1;
                end
            end
            
            for i = 1:length(Branch)
                if ~strcmp(Branch{i}.Status, 'Open')
                    s{indexLink} = Branch{i}.connectedObj1;
                    t{indexLink} =  Branch{i}.connectedObj2;
                    indexLink = indexLink + 1;
                end
            end
            
            % check Nodes in cells of s and t are in the table of G.Nodes or not
            if ~isempty(s) && ~isempty(t)
                st = cell2table(s'); tt = cell2table(t');
                st.Properties.VariableNames = {'Name'}; tt.Properties.VariableNames = {'Name'};
                Lst = [ismember(st,G.Nodes), ismember(tt,G.Nodes)];

                % find nodes in s and t that are also in G.Nodes 
                id = find(sum(Lst,2) == 2); 
                if length(id)>0
                    for ii = 1:length(id)
                        s2{ii} = s{id(ii)};
                        t2{ii} = t{id(ii)};
                    end
                    G = rmedge(G,s2,t2);
                end
            end
            
            Nrmnode = indexNode;
            Nrmlink = indexLink;

        end
        
        %% Communication: Remove damaged nodes and edges after the damage assesment (STEP 2)
        function [G, Nrmnode, Nrmlink] = RemoveCommGraph(Comm, G)
            % ================================================================
            % function RemoveCommGraph
            % Remove Non-Open edges and the nodes of these Non-Open edges from 
            % the Communication network after the initial damage assesment
            % -----------------------------------------------------------------
            % node = Centraloffice/CommunicationTower
            % edge = Cellline
            % s = start node of a Cellline
            % t = end node of a Cellline
            % s2 = unique link start nodes in the orginal communication graph
            % t2 = unique link end nodes in the orginal communication graph
            % ================================================================

            Centraloffice_Set = Comm{1};
            CommunicationTower_Set = Comm{2};
            Cellline_Set = Comm{3};
            
            s = [];
            t = [];
            st = [];
            tt = [];
            Lst = [];
            s2 = [];
            t2 = [];
            
            indexNode = 0;
            indexLink = 1;
            
            for ii = 1:length(Centraloffice_Set)
                if ~strcmp(Centraloffice_Set{ii}.Status, 'Open')
                    G = rmnode(G,Centraloffice_Set{ii}.uniqueID);
                    indexNode = indexNode+1;
                end
            end
            
            for ii = 1:length(CommunicationTower_Set)
                if ~strcmp(CommunicationTower_Set{ii}.Status, 'Open')
                    G = rmnode(G,CommunicationTower_Set{ii}.uniqueID);
                    indexNode = indexNode+1;
                end
            end
            
            for ii = 1:length(Cellline_Set)
                if ~strcmp(Cellline_Set{ii}.Status, 'Open')
                    s{indexLink} = Cellline_Set{ii}.connectedObj1;
                    t{indexLink} =  Cellline_Set{ii}.connectedObj2;
                    indexLink = indexLink + 1;
                end
            end

            % check Nodes in cellls of s and t in the table of G.Nodes or not
            if ~isempty(s) && ~isempty(t)
                st = cell2table(s'); tt = cell2table(t');
                st.Properties.VariableNames = {'Name'}; tt.Properties.VariableNames = {'Name'};
                Lst = [ismember(st,G.Nodes), ismember(tt,G.Nodes)];

                % find nodes in s and t that are also in G.Nodes 
                id = find(sum(Lst,2) == 2); 
                if length(id)>0
                    for ii = 1:length(id)
                        s2{ii} = s{id(ii)};
                        t2{ii} = t{id(ii)};
                    end
                    G = rmedge(G,s2,t2);
                end  
            end
            
            Nrmnode = indexNode;
            Nrmlink = indexLink;

        end

        %% Select Scheduler Model to develop a restoration plan (STEP 3)
        % 1 - Restoration ranking based on the
        % priority/criticality/importance of individual components
        % 2 - Restoration ranking based on the optimal schedule from a
        % optimization  formulation
        function [Schedule, Date] = RestorationPlan(Scheduler_Num, time_horizon, RestorationResource, ...
                active_power, active_comm, active_trans, SystemSetDamage, TableTask, TablePrecedence,... 
                priority_power_num, priority_transportation_num, priority_communication_num, OptimizationChoice)
            
            %=== prepare input data of task and precedence into tables for every system
            %- Power_Set, Communication_Set, Transportation_Set = system sets;
            %- taskP, taskC, taskT = the table lists of restroation tasks
            % for power, communication, and transportation, respectively;
            %- preP, taskC, taskT = the table lists of restroation tasks
            % for power, communication, and transportation, respectively.
            
            Power_Set = SystemSetDamage{1};
            Communication_Set = SystemSetDamage{2};
            Transportation_Set = SystemSetDamage{3};
            
            taskTable{1} = TableTask{1,1}; %taskP 
            taskTable{2} = TableTask{1,2}; %taskC 
            taskTable{3} = TableTask{1,3}; %taskT
            precedenceTable{1} = TablePrecedence{1,1}; %preP 
            precedenceTable{2} = TablePrecedence{1,2}; %preC 
            precedenceTable{3} = TablePrecedence{1,3}; %preT 
            
            %===== Develop a restoration plan based on Scheduler_Num =====
            % 1-Scheme 3A: priority
            % 2-Scheme 3B: optimization
            switch Scheduler_Num
                %=== Scheme 3A: develop a restoration plan based on the
                %priority ranking of individual object (component)
                case 1 
                    %--- priority_power_num
                    if active_power && ~isempty(taskTable{1})
                        [result_pow, pow_date] = Library.PowerSchedulePriority(priority_power_num,Power_Set);
                    else
                        result_pow = [];
                        pow_date = [];
                    end
                    
                    %--- priority_communication_num
                    if active_comm && ~isempty(taskTable{2})
                        [result_comm, comm_date] = Library.CommSchedulePriority(priority_communication_num,Communication_Set);
                    else
                        result_comm = [];
                        comm_date = [];
                    end      
                    
                    %--- priority_transportation_num
                    % 1: ranking creteria: bridge traffic
                    % 2: ranking creteria: bridge length
                    if active_trans && ~isempty(taskTable{3})
                        [result_trans, trans_date] = Library.TransSchedulePriority(priority_transportation_num,Transportation_Set);
                    else
                        result_trans = [];
                        trans_date = [];
                    end

                    Schedule = {result_pow, result_comm, result_trans};
                    Date = {pow_date, comm_date, trans_date};
                
                %=== Scheme 3B: develop a restoration plan based on the
                %optimal schedule solution of all tasks based on a chosen 
                % optimization formulation (=OptimizationChoice)     
                case 2
                    %resource = [Max_Power; Max_Trans; Max_Comm];
                    resource = RestorationResource;
                    taskP = TableTask{1,1};
                    taskC = TableTask{1,2};
                    taskT = TableTask{1,3};                    
                    
                    if active_power && ~isempty(taskP)
                        [result_pow, pow_date] = Library.RepairScheduleOptimization(OptimizationChoice, 'Power', taskTable, precedenceTable, resource, time_horizon);
                    else
                        result_pow = [];
                        pow_date = [];
                    end
                                        
                    if active_comm && ~isempty(taskC)
                        [result_comm, comm_date] = Library.RepairScheduleOptimization(OptimizationChoice, 'Communication', taskTable, precedenceTable, resource, time_horizon);
                    else
                        result_comm = [];
                        comm_date = [];
                    end
                    
                    if active_trans && ~isempty(taskT)
                        [result_trans, trans_date] = Library.RepairScheduleOptimization(OptimizationChoice, 'Transportation', taskTable, precedenceTable, resource, time_horizon);
                    else
                        result_trans = [];
                        trans_date = [];
                    end

                    Schedule = {result_pow, result_comm, result_trans};
                    Date = {pow_date, comm_date, trans_date};
            end
    
        end
        

        %% Scheme 3A: Priority Scheduler for Transportation System (STEP 3)
        function [Schedule, Return_Date] = TransSchedulePriority(num,Set)
            %==============================================================
            % function TransSchedulePriority
            % This function develops a restoration plan for the transportation system based the object priority. 
            %=== priority creteria
            %- I for Bridge, the alternative ranking creterion are:
            % 0. Randomly
            % 1. AADT (traffic)
            % 2. length
            %- II for Road
            % 0. Randomly
            % 1. AADT (traffic)
            % 2. length
            %- III for TrafficLight
            % 0. Randomly
            % 1. the ranking is based on the traffic (AADT) on the road that
            % this traffic light serves (controlling traffic). 
            %==============================================================
            Return_Date = [];
            Road = Set{1};
            Bridge = Set{2};
            TrafficLight = Set{3};
            Rank = []; RankRoad = []; RankTrafficLight = []; 
            index = 1;
            %===============Bridge
            if num == 0 % develop a restoration plan that repairs the damaged component in a random order
                seed = 310;
                rng(seed);
                N = length(Bridge);
                RandomNumber = rand(1,N);
                
                for ii = 1:N
                    if strcmp(Bridge{ii}.Status, 'Damaged')
                        Rank(end+1,:) = [ii, RandomNumber(ii)];
                    end
                end
            
            elseif num == 1 % bridge traffic
                for ii = 1:length(Bridge)
                    if strcmp(Bridge{ii}.Status, 'Damaged')
                        Rank(end+1,:) = [ii, Bridge{ii}.Traffic];
                    end
                end
            elseif num == 2 % bridge length
                for ii = 1:length(Bridge)
                    if strcmp(Bridge{ii}.Status, 'Damaged')
                        Rank(end+1,:) = [ii, Bridge{ii}.Length];
                    end
                end
            end
            
            if ~isempty(Rank)
                sortedBridge = sortrows(Rank,2);
                for ii = 1:size(sortedBridge,1)
                    bridgeCurrent = sortedBridge(size(sortedBridge,1) - ii + 1,1);
                    tasks = Bridge{bridgeCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('Bridge/',num2str(bridgeCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} = strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end
 

 
            %===============Road
            if num == 0 % develop a restoration plan that repairs the damaged component in a random order
                seed = 320;
                rng(seed);
                N = length(Road);
                RandomNumber = rand(1,N);
                
                for ii = 1:N
                    if strcmp(Road{ii}.Status, 'Damaged')
                        RankRoad(end+1,:) = [ii, RandomNumber(ii)];
                    end
                end
            elseif num == 1 % road traffic
                for ii = 1:length(Road)
                    if strcmp(Road{ii}.Status, 'Damaged')
                        RankRoad(end+1,:) = [ii, Road{ii}.AADT];
                    end
                end
            elseif num == 2 % road length
                for ii = 1:length(Road)
                    if strcmp(Road{ii}.Status, 'Damaged')
                        RankRoad(end+1,:) = [ii, Road{ii}.Length];
                    end
                end
            end
            
            if ~isempty(RankRoad)
                sortedRoad = sortrows(RankRoad,2);
                for ii = 1:size(sortedRoad,1)
                    RoadCurrent = sortedRoad(size(sortedRoad,1) - ii + 1,1);
                    tasks = Road{RoadCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('Road/',num2str(RoadCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end
            %===============TrafficLight
            if num == 0 % develop a restoration plan that repairs the damaged component in a random order
                seed = 320;
                rng(seed);
                N = length(TrafficLight);
                RandomNumber = rand(1,N);
                
                for ii = 1:N
                    if strcmp(TrafficLight{ii}.Status, 'Damaged')
                        RankTrafficLight(end+1,:) = [ii, RandomNumber(ii)];
                    end
                end
            elseif num == 1 % traffic (AADT) of the road that the this traffic light controls
                for ii = 1:length(TrafficLight)
                    if strcmp(TrafficLight{ii}.Status, 'Damaged')
                        iroad = TrafficLight{ii}.RoadLinkID;
                        RankTrafficLight(end+1,:) = [ii, Road{iroad}.AADT];
                    end
                end
            
            end
            
            if ~isempty(RankTrafficLight)
                sortedTrafficLight = sortrows(RankTrafficLight,2);
                for ii = 1:size(sortedTrafficLight,1)
                    TrafficLightCurrent = sortedTrafficLight(size(sortedTrafficLight,1) - ii + 1,1);
                    tasks = TrafficLight{TrafficLightCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('TrafficLight/',num2str(TrafficLightCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end
%             for ii = 1:length(TrafficLight)
%                 if strcmp(TrafficLight{ii}.Status, 'Damaged')
%                     tasks = TrafficLight{ii}.taskUniqueIds;
%                     for j = 1:length(tasks)
%                         Schedule{index} = strcat('TrafficLight/',num2str(ii));
%                         Schedule{index} = strcat(Schedule{index}, '/');
%                         Schedule{index} =  strcat(Schedule{index}, tasks{j});
%                         index = index + 1;
%                     end
%                 end
%             end															 																	  
            
        end % end of function TransSchedulePriority


        %% Priority Scheduler for the Power System (STEP 3)
        % num: a index representing the analyst's choice on the priority policy from the input file to develop the restoration plan: 
        % num = 0 representing that developing a restoration plan by
        % randomly repairing damaged component
        % num = 1 representing the (maximum) voltage that this object serves,
        % num = 2 representing the population that this object serves
        % num = 3 representing the number of household that this object serves
        function [Schedule, Return_Date] = PowerSchedulePriority(num,Set)
            %==============================================================
            % function PowerSchedulePriority
            % develop a restoration plan for the power system based the object priority 
            %=== priority creteria
            %- I. for Bus (substation), the alternative ranking creterion are:
            % 0. randomly
            % 1. capacity (voltage)
            % 2. population served
            % 3. number of household served
            %- II. for Branch (power line)
            % 0. randomly
            % 1&2&3. the ranking is based on the voltage of the branch
            % (power line) carries.
            %- III. for TransmissionTower
            % 0. randomly
            % 1&2&3. the ranking is based on the voltage of the transmission tower
            % carries.
            %==============================================================
            Return_Date = [];
            Branch = Set{1};
            Bus= Set{2};
            Generator = Set{3};
            TransmissionTower = Set{4};
            Rank = []; RankBranch = []; RankTT = [];
            index = 1;
            
            %============ Bus (Substation)
            if num == 0 % develop a restoration plan that repairs the damaged component in a random order
                seed = 110;
                rng(seed);
                N = length(Bus);
                RandomNumber = rand(1,N);               
                
                for ii = 1:N
                    if strcmp(Bus{ii}.Status, 'Damaged')
                        Rank(end+1,:) = [ii, RandomNumber(ii)];
                    end
                end
            elseif num == 1 % voltage
                for ii = 1:length(Bus)
                    if strcmp(Bus{ii}.Status, 'Damaged')
                        Rank(end+1,:) = [ii, Bus{ii}.Capacity];
                    end
                end
            elseif num == 2 % Population served 
                for ii = 1:length(Bus)
                    if strcmp(Bus{ii}.Status, 'Damaged')
                        Rank(end+1,:) = [ii, Bus{ii}.PopulationServed];
                    end                    
                end
            elseif num == 3 % number of household served 
                for ii = 1:length(Bus)
                    if strcmp(Bus{ii}.Status, 'Damaged')
                        Rank(end+1,:) = [ii, Bus{ii}.HouseholdServed];
                    end                    
                end    
            end
            
            if ~isempty(Rank)
                sortedBus = sortrows(Rank,2);
                for ii = 1:size(sortedBus,1)
                    busCurrent = sortedBus(size(sortedBus,1) - ii + 1,1);
                    tasks = Bus{busCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('Bus/',num2str(busCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end
            
            
            %============ Branch (Power Line)
            if num == 0 % develop a restoration plan that repairs the damaged component in a random order
                seed = 120;
                rng(seed);
                N = length(Branch);
                RandomNumber = rand(1,N);               
                
                for ii = 1:N
                    if strcmp(Branch{ii}.Status, 'Damaged')
                        RankBranch(end+1,:) = [ii, RandomNumber(ii)];
                    end
                end
            else % voltage of this branch (power line) carries
                for ii = 1:length(Branch)
                    if strcmp(Branch{ii}.Status, 'Damaged')
                        RankBranch(end+1,:) = [ii, Branch{ii}.Capacity];
                    end
                end
            end
            
            if ~isempty(RankBranch)
                sortedBranch = sortrows(RankBranch, 2);
                for ii = 1:size(sortedBranch,1)
                    BranchCurrent = sortedBranch(size(sortedBranch,1) - ii + 1, 1); 
                    tasks = Branch{BranchCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('Branch/',num2str(BranchCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end % if ~ RankBranch
            
            %============ Transmission Tower
            if num == 0 % develop a restoration plan that repairs the damaged component in a random order
                seed = 120;
                rng(seed);
                N = length(TransmissionTower);
                RandomNumber = rand(1,N);               
                
                for ii = 1:N
                    if strcmp(TransmissionTower{ii}.Status, 'Damaged')
                        RankTT(end+1,:) = [ii, RandomNumber(ii)];
                    end
													  
                end
            else % voltage of the power line that this transmission tower carries 
                for ii = 1:length(TransmissionTower)
                    if strcmp(TransmissionTower{ii}.Status, 'Damaged')
                        TTbranch = TransmissionTower{ii}.Branch;  
                        if length(TTbranch) == 2
                            TTcapacity = max(Branch{TTbranch(1)}.Capacity, Branch{TTbranch(2)}.Capacity);
                        elseif length(TTbranch) == 1
                            TTcapacity = Branch{TTbranch(1)}.Capacity;
                        end
                        RankTT(end+1,:) = [ii, TTcapacity];
                    end
                end %i ~ length(TransmissionTower)
            end
            
            if ~isempty(RankTT)
                sortedTT = sortrows(RankTT, 2);
                for ii = 1:size(sortedTT,1)
                    TTCurrent = sortedTT(size(sortedTT,1) - ii + 1, 1); 
                    tasks = TransmissionTower{TTCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('TransmissionTower/',num2str(TTCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end % RankTT

        end % end of function PowerSchedulePriority
        
        %% Scheme 3A: Priority, Communication System (STEP 3)
        % num: a index representing the analyst's choice on the priority policy from the input file to develop the restoration plan: 
        % num = 0 representing that developing a restoration plan by
        % randomly repairing damaged component
        % num = 1 representing the population that this object serves
        % num = 2 representing the number of household that this object serves
        function [Schedule, Return_Date] = CommSchedulePriority(num,Set)
            %==============================================================
            % function CommSchedulePriority
            % develop a restoration plan for the communication system based the object priority 
            %=== priority creteria
            %- I. for Central Office, the alternative ranking creterion are:
            % 0. randomly
            % 1. population served
            % 2. number of household served
            %- II. for CommunicationTower
            % 0. randomly
            % 1. population served
            % 2. number of household served
            %- III. for CommunicationTower
            % 0. randomly
            % 1. the maximum population served th connected central office/
            % communication tower
            % 2. the maximum number of household served th connected central office/
            % communication tower
            %==============================================================
            
            Return_Date = [];
            
            Centraloffice = Set{1};
            CommunicationTower = Set{2};
            Cellline = Set{3};
  
            Rank = []; RankCT = []; RankCL = [];
            index = 1;
            
            %================= Central Office
            if num == 0 % randomly
                seed = 120;
                rng(seed);
                N = length(Centraloffice);
                RandomNumber = rand(1,N);               

                for iobj = 1:N
                    if strcmp(Centraloffice{iobj}.Status, 'Damaged')
                        Rank(end+1,:) = [iobj, RandomNumber(iobj)];
                    end
                end
            elseif num == 1 % population served
                for iobj = 1:length(Centraloffice)
                    if strcmp(Centraloffice{iobj}.Status, 'Damaged')
                        Rank(end+1,:) = [iobj, Centraloffice{iobj}.PopulationServed];
                    end
                end
            elseif num == 2 % number of household served
                for iobj = 1:length(Centraloffice)
                    if strcmp(Centraloffice{iobj}.Status, 'Damaged')
                        Rank(end+1,:) = [iobj, Centraloffice{iobj}.HouseholdServed];
                    end
                end
            end
            
            if ~isempty(Rank)
                sortedCentral = sortrows(Rank,2);
                for iobj = 1:size(sortedCentral,1)
                    centralCurrent = sortedCentral(size(sortedCentral,1) - iobj + 1,1);
                    tasks = Centraloffice{centralCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('CentralOffice/',num2str(centralCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end
            
            
            %================= Communication tower
            if num == 0 % randomly
                seed = 120;
                rng(seed);
                N = length(CommunicationTower);
                RandomNumber = rand(1,N);               

                for iobj = 1:N
                    if strcmp(CommunicationTower{iobj}.Status, 'Damaged')
                        RankCT(end+1,:) = [iobj, RandomNumber(iobj)];
                    end
                end
            elseif num == 1 % population served
                for iobj = 1:length(CommunicationTower)
                    if strcmp(CommunicationTower{iobj}.Status, 'Damaged')
                        RankCT(end+1,:) = [iobj, CommunicationTower{iobj}.PopulationServed];
                    end
                end
            elseif num == 2 % number of household served
                for iobj = 1:length(CommunicationTower)
                    if strcmp(CommunicationTower{iobj}.Status, 'Damaged')
                        RankCT(end+1,:) = [iobj, CommunicationTower{iobj}.HouseholdServed];
                    end
                end
            end
            
            if ~isempty(RankCT)
                sortedCT = sortrows(RankCT,2);
                for iobj = 1:size(sortedCT,1)
                    CTCurrent = sortedCT(size(sortedCT,1) - iobj + 1,1);
                    tasks = CommunicationTower{CTCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('CommunicationTower/',num2str(CTCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end
            
%             for iobj = 1:length(CommunicationTower)
%                 if strcmp(CommunicationTower{iobj}.Status, 'Damaged')
%                     tasks = CommunicationTower{iobj}.taskUniqueIds;
%                     for j = 1:length(tasks)
%                         Schedule{index} = strcat('CommunicationTower/',num2str(iobj));
%                         Schedule{index} = strcat(Schedule{index}, '/');
%                         Schedule{index} =  strcat(Schedule{index}, tasks{j});
%                         index = index + 1;
%                     end
%                 end
%             end
            
            %================= Celline (Communication line)
            if num == 0 % randomly
                seed = 120;
                rng(seed);
                N = length(Cellline);
                RandomNumber = rand(1,N);               

                for iobj = 1:N
                    if strcmp(Cellline{iobj}.Status, 'Damaged')
                        RankCL(end+1,:) = [iobj, RandomNumber(iobj)];
                    end
                end
            elseif num == 1 % population served
                for iobj = 1:length(Cellline)
                    if strcmp(Cellline{iobj}.Status, 'Damaged')
                        obj1 = Cellline{iobj}.connectedObj1; 
                        obj2 = Cellline{iobj}.connectedObj2;
                        if contains(obj1, 'Centraloffice')
                            pop1 = Centraloffice{iobj}.PopulationServed;
                        elseif contains(obj1, 'CommunicationTower')
                            pop1 = CommunicationTower{iobj}.PopulationServed;
                        end
                        if contains(obj2, 'Centraloffice')
                            pop2 = Centraloffice{iobj}.PopulationServed;
                        elseif contains(obj2, 'CommunicationTower')
                            pop2 = CommunicationTower{iobj}.PopulationServed;
                        end
                        
                        pop = max([pop1, pop2]); 
                        RankCL(end+1,:) = [iobj, pop];
                    end
                end
            elseif num == 2 % number of household served
                for iobj = 1:length(Cellline)
                    if strcmp(Cellline{iobj}.Status, 'Damaged')
                        obj1 = Cellline{iobj}.connectedObj1; 
                        obj2 = Cellline{iobj}.connectedObj2;
                        if contains(obj1, 'Centraloffice')
                            ii = sscanf(obj1,'Centraloffice%d');
                            hh1 = Centraloffice{ii}.HouseholdServed;
                        elseif contains(obj1, 'CommunicationTower')
                            ii = sscanf(obj1,'CommunicationTower%d');
                            hh1 = CommunicationTower{ii}.HouseholdServed;
                        end
                        if contains(obj2, 'Centraloffice')
                            ii = sscanf(obj1,'Centraloffice%d');
                            hh2 = Centraloffice{ii}.HouseholdServed;
                        elseif contains(obj2, 'CommunicationTower')
                            ii = sscanf(obj1,'CommunicationTower%d');
                            hh2 = CommunicationTower{ii}.HouseholdServed;
                        end
                        
                        hh = max([hh1, hh2]); 
                        RankCL(end+1,:) = [iobj, hh];
                    end
                end
            end
            
            if ~isempty(RankCL)
                sortedCL = sortrows(RankCL,2);
                for iobj = 1:size(sortedCL,1)
                    CLCurrent = sortedCL(size(sortedCL,1) - iobj + 1,1);
                    tasks = CommunicationTower{CLCurrent}.taskUniqueIds;
                    for j = 1:length(tasks)
                        Schedule{index} = strcat('Cellline/',num2str(CLCurrent));
                        Schedule{index} = strcat(Schedule{index}, '/');
                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                        index = index + 1;
                    end
                end
            end


%             for iobj = 1:length(Cellline)
%                 if strcmp(Cellline{iobj}.Status, 'Damaged') 
%                     tasks = Cellline{iobj}.taskUniqueIds;
%                     for j = 1:length(tasks)
%                         Schedule{index} = strcat('Cellline/',num2str(iobj));
%                         Schedule{index} = strcat(Schedule{index}, '/');
%                         Schedule{index} =  strcat(Schedule{index}, tasks{j});
%                         index = index + 1;
%                     end
%                 end
%             end
 
        end


        %% Scheme 3A: Create repair schedule according to the importatance of component (STEP 3)
        function [Return_Schedule, Return_Date] = RepairSchedulePriority(System, Set, Max)
            Schedule = {};
            if strcmp(System, 'Power')
                Branch = Set{1};
                Bus = Set{2};
                Generator = Set{3};
                TransmissionTower = Set{4};
                
                current_priority = 1;
                index = 1;
                
                while(current_priority <= Max)
                    change = 0;
                    
                    for i = 1:length(Bus)
                        if strcmp(Bus{i}.Status, 'Damaged') && Bus{i}.Priority == current_priority
                            tasks = Bus{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('Bus/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end
                    for i = 1:length(Generator)
                        if strcmp(Generator{i}.Status, 'Damaged') && Generator{i}.Priority == current_priority
                            Schedule{index} = strcat('Generator/',num2str(i));
                            index = index + 1;
                        end
                    end
                    
                    for i = 1:length(Branch)
                        if strcmp(Branch{i}.Status, 'Damaged') && Branch{i}.Priority == current_priority
                            tasks = Branch{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('Branch/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end

                    for i = 1:length(TransmissionTower)
                        if strcmp(TransmissionTower{i}.Status, 'Damaged') && TransmissionTower{i}.Priority == current_priority
                            tasks = TransmissionTower{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('TransmissionTower/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end
                    
                    current_priority = current_priority + 1;
                end
                
            elseif strcmp(System, 'Communication')
                
                Centraloffice = Set{1};
                CommunicationTower = Set{2};
                Cellline = Set{3};
                
                current_priority = 1;
                index = 1;
                
                while(current_priority <= Max)
                    change = 0;
                    
                    for i = 1:length(Centraloffice)
                        if strcmp(Centraloffice{i}.Status, 'Damaged') && Centraloffice{i}.Priority == current_priority
                            tasks = Centraloffice{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('CentralOffice/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end
                    
                    for i = 1:length(Cellline)
                        if strcmp(Cellline{i}.Status, 'Damaged') && Cellline{i}.Priority == current_priority
                            tasks = Cellline{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('Cellline/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end
                    
                    for i = 1:length(CommunicationTower)
                        if strcmp(CommunicationTower{i}.Status, 'Damaged') && CommunicationTower{i}.Priority == current_priority
                            tasks = CommunicationTower{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('CommunicationTower/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end
                    
                    current_priority = current_priority + 1;
                end
                
            elseif strcmp(System, 'Transportation')
                Road = Set{1};
                Bridge = Set{2};
                TrafficLight = Set{3};
                
                current_priority = 1;
                index = 1;
                
                while(current_priority <= Max)
                    change = 0;
                    for i = 1:length(Road)
                        if strcmp(Road{i}.Status, 'Damaged') && Road{i}.Priority == current_priority
                            tasks = Road{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('Road/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end
                    
                    for i = 1:length(Bridge)
                        if strcmp(Bridge{i}.Status, 'Damaged') && Bridge{i}.Priority == current_priority
                            
                            
%                                 tasks = Bridge{i}.taskUniqueIds;
%                                 for j = 1:length(tasks)
%                                     Schedule{index} = strcat('Bridge/',num2str(i));
%                                     Schedule{index} = strcat(Schedule{index}, '/');
%                                     Schedule{index} =  strcat(Schedule{index}, tasks{j});
%                                     index = index + 1;
%                                 end
                            
                            if Bridge{i}.HasSub == 1
                                for sub_index = 1:length(Bridge{i}.ColumnSet)
                                    tasks = Bridge{i}.ColumnSet(sub_index).taskUniqueIds;
                                    for j = 1:length(tasks)
                                        Schedule{index} = strcat('Column/',num2str(i));
                                        Schedule{index} = strcat(Schedule{index}, '/');
                                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                        index = index + 1;
                                    end
                                end
                                for sub_index = 1:length(Bridge{i}.ColumnFoundSet)
                                    tasks = Bridge{i}.ColumnFoundSet(sub_index).taskUniqueIds;
                                    for j = 1:length(tasks)
                                        Schedule{index} = strcat('ColumnFoundation/',num2str(i));
                                        Schedule{index} = strcat(Schedule{index}, '/');
                                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                        index = index + 1;
                                    end
                                end
                                for sub_index = 1:length(Bridge{i}.AbutmentSet)
                                    tasks = Bridge{i}.AbutmentSet(sub_index).taskUniqueIds;
                                    for j = 1:length(tasks)
                                        Schedule{index} = strcat('Abutment/',num2str(i));
                                        Schedule{index} = strcat(Schedule{index}, '/');
                                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                        index = index + 1;
                                    end
                                end
                                for sub_index = 1:length(Bridge{i}.AbutmentFoundSet)
                                    tasks = Bridge{i}.AbutmentFoundSet(sub_index).taskUniqueIds;
                                    for j = 1:length(tasks)
                                        Schedule{index} = strcat('AbutmentFoundation/',num2str(i));
                                        Schedule{index} = strcat(Schedule{index}, '/');
                                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                        index = index + 1;
                                    end
                                end
                                for sub_index = 1:length(Bridge{i}.BearingSet)
                                    tasks = Bridge{i}.BearingSet(sub_index).taskUniqueIds;
                                    for j = 1:length(tasks)
                                        Schedule{index} = strcat('Bearing/',num2str(i));
                                        Schedule{index} = strcat(Schedule{index}, '/');
                                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                        index = index + 1;
                                    end
                                end
                                for sub_index = 1:length(Bridge{i}.SlabSet)
                                    tasks = Bridge{i}.SlabSet(sub_index).taskUniqueIds;
                                    for j = 1:length(tasks)
                                        Schedule{index} = strcat('ApproachSlab/',num2str(i));
                                        Schedule{index} = strcat(Schedule{index}, '/');
                                        Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                        index = index + 1;
                                    end
                                end
                            else
                                tasks = Bridge{i}.taskUniqueIds;
                                for j = 1:length(tasks)
                                    Schedule{index} = strcat('Bridge/',num2str(i));
                                    Schedule{index} = strcat(Schedule{index}, '/');
                                    Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                    index = index + 1;
                                end
                            end
                        end
                    end
                    
                    for i = 1:length(TrafficLight)
                        if strcmp(TrafficLight{i}.Status, 'Damaged') && TrafficLight{i}.Priority == current_priority
                            tasks = TrafficLight{i}.taskUniqueIds;
                            for j = 1:length(tasks)
                                Schedule{index} = strcat('TrafficLight/',num2str(i));
                                Schedule{index} = strcat(Schedule{index}, '/');
                                Schedule{index} =  strcat(Schedule{index}, tasks{j});
                                index = index + 1;
                            end
                        end
                    end
                    
                    current_priority = current_priority + 1;
                    
                end
                
            end
            Return_Schedule = Schedule;
            Return_Date = [];
        end

        %% Scheme 3B: Create repair schedule according to the optimization  (STEP 3)
        function [Return_Schedule, Return_Date] = RepairScheduleOptimization(OptimizationChoice,System, taskTable, precedenceTable, resource, time_horizon)
            Return_Schedule = {};
            Return_Date = [];
            
            if strcmp(System, 'Power')
                ii = 1;
                task = taskTable{ii};
                resourceConstraint = resource(ii,:); 
                precedence = cell2mat(precedenceTable{ii}(2:end,2:end));
                
                switch OptimizationChoice
                    case 1 % RCPSP_FT
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 2 % RCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                    case 3 % Flexible
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationFlexible(task, precedence, resourceConstraint, time_horizon);
                    case 4 % SuboptimizationRCPSP_FT    
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 5 % SuboptimizationRCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                
                end
                
            elseif strcmp(System, 'Communication')
                ii = 2;
                task = taskTable{ii};
                resourceConstraint = resource(ii,:); 
                precedence = cell2mat(precedenceTable{ii}(2:end,2:end));
                
                switch OptimizationChoice
                    case 1 % RCPSP_FT
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 2 % RCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                    case 3 % Flexible
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationFlexible(task, precedence, resourceConstraint, time_horizon);
                    case 4 % SuboptimizationRCPSP_FT    
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 5 % SuboptimizationRCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                
                end
                
            elseif strcmp(System, 'Transportation')
                ii = 3;              
                task = taskTable{ii};
                resourceConstraint = resource(ii,:); 
                precedence = cell2mat(precedenceTable{ii}(2:end,2:end));
                switch OptimizationChoice
                    case 1 % RCPSP_FT
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 2 % RCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                    case 3 % Flexible
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationFlexible(task, precedence, resourceConstraint, time_horizon);
                    case 4 % SuboptimizationRCPSP_FT    
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 5 % SuboptimizationRCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                
                end
                
            % schedule tasks of all three systems together
            % convert the restoration time following the same unit of DAY 
            elseif strcmp(System, 'All')
                is = 1; newTaskTable{is} = taskTable{is}; 
                for id = 1:size(newTaskTable{is},1)
                    newTaskTable{is}{id,6} = inv(24)*taskTable{is}{id,6};
                    newTaskTable{is}{id,7} = inv(24)*taskTable{is}{id,7};
                    newTaskTable{is}{id,8} = inv(24)*taskTable{is}{id,8};
                end
                
                is = 2; newTaskTable{is} = taskTable{is};
                for id = 1:size(newTaskTable{is},1)
                    newTaskTable{is}{id,6} = inv(24)*taskTable{is}{id,6};
                    newTaskTable{is}{id,7} = inv(24)*taskTable{is}{id,7};
                    newTaskTable{is}{id,8} = inv(24)*taskTable{is}{id,8};
                end   
                    
                task = cat(1, newTaskTable{1}, newTaskTable{2},taskTable{3});

                resourceConstraint = sum(resource); 
                precedence = zeros(size(task,1),size(task,1));
                precedence(1:size(taskTable{1},1)) = cell2mat(precedenceTable{1}(2:end,2:end));
                precedence(1+size(taskTable{1},1):size(taskTable{2},1)+size(taskTable{1},1)) = cell2mat(precedenceTable{2}(2:end,2:end));
                precedence(1+size(taskTable{1},1)+size(taskTable{2},1):end) = cell2mat(precedenceTable{3}(2:end,2:end));
                
                switch OptimizationChoice
                    case 1 % RCPSP_FT
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 2 % RCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                    case 3 % Flexible
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.optimizationFlexible(task, precedence, resourceConstraint, time_horizon);
                    case 4 % SuboptimizationRCPSP_FT    
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon);
                    case 5 % SuboptimizationRCPSP_ST
                        [X,ct,sft,ru,sequenceTask,sequenceObject] = Library.SuboptimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
                    
                end
                
            end
            
            Return_Schedule = sequenceTask; 
            Return_Date = sft(:,1:2); % [st, ft] 
            
        end
        
        %% Scheme 3B: Optimization, RCPSP_FT  (STEP 3)
        function [X,ct,sft,ru,sequenceTask,sequenceObject] = optimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon)
            % ================================================================
            % function optimizationRCPSP_FT
            % updated repair schedule using the formulation of resource-constrained project
            % scheduling problem (RCPSP)
            %Robert Klein (2000), Scheduling of resource-constrained projects, Kluwer
            %Academic Publishers, page 79-80. 
            % binary decision variable x_{it} = 1, when task i finishes at (time) the end of period t;
            % x_{it} = 0 otherwise.
            % the objective to is to minimize the final completion time of
            % all tasks.
            % ================================================================
            %==== Input parameters
            % (1) task = {1-'TaskID in the system',...
            % 2-'Object Type',...
            % 3-'Object ID',...
            % 4-'TaskID in the task library of this Object type', ...
            % 5-'Decription of this task',...
            % 6-'Duration min', 7-'Duration max', 8-'Duration mode', 9-'Duration distribution type',...
            % 10-'Duration sample',...
            % 11-'Resource demand for resource type 1, 2, 3, 4, respectively'}
            % All aforementioned number labels mean the column number in the
            % cell of "task"
            % (2) precedence, binary, = 1, if the task labeled as column
            % number is the precedessory of the task labeled as row number 
            % (3) resourceConstraint = ar = [AR1 AR2, AR3 AR4], avaiability
            % for 4 types of resource, constant over time
            % (4) time_horizon, the entire time horizon of interest
            %==== Output parameters
            % (1) X = opx matrix, size=I*th, 
            % every row represent the decision variable value over t=[1,2, .., th]
            % every column presente the decision variable value at a time snapshot
            % (2) ct = completion time of all tasks
            % (3) sft = [st,ft,taskID]
            % (4) sequenceTask = a ranking list of tasks in terms of restoration time
            % from early to late
            % (5) sequenceObject = a ranking list of damaged objects in terms of restoration time
            % from early to late
            % 
            %%=========================================================
            % d: restoration duration as the mode value in column 8 in task
            % H: time horizon
            % I: number of resource type
            % P: precedence matrix 
            % R: renewable resource constraint
            % r: resource demand as the values in column 11~14 in task
            % z: an intermediate variable as taskID in this function only 
            % from 1 to I
            
            %======================= input parameters===========================
            H = time_horizon;               % time horizon
            K = length(resourceConstraint); % the number of resource type
            I = size(task,1);               % the number of tasks
            R = resourceConstraint;         % resource constraints
            P = precedence;                 % precedence

            %t = linspace(0,H,H+1);
            t = linspace(1,H,H);
            z = linspace(1,I,I); 

            % assign duration d and resource requirement r for every task 
            % during schedule, task duration d = duration mode
                for ii = 1:I
                    d(ii) = max(1, round(task{ii, 8}));
                    if isempty(d(ii)) || isnan(d(ii))
                        d(ii) = max(1, round(0.5*(task{ii,6}+task{ii,7})));
                    end
                    r(ii,:) = task{ii, 11};
                end

            % ====================== decision variables===========================
            %x = binvar(I,H+1); % which time step task i(=1,2,...,I,I+1) is finished within [0, H]
            x = binvar(I,H); % at which time period task i(=1,2,...,I) is finished within [1, H]
            
            % =========================constraints===============================
            constraints = [];

            % c1: every task only executes once. 
            for ii = 1:I
                %constraints = [constraints, sum(x(ii,d(ii)+1:H+1)) - 1 == 0]; 
                constraints = [constraints, sum(x(ii,d(ii):H)) - 1 == 0];
            end

            % c2: resource usage should be less than resource constraint per time step
            for ir = 1:K % K-1
                for it = 1:H %H+1
%                     for ii = 1:I
%                         %ub = min(it+d(ii)-1,H+1);
%                         ub = min(it+d(ii)-1,H);
%                         lb = max(it,d(ii));
%                         if ub-lb>=0
%                             tmp21(ii) = sum(x(ii,lb:ub));
%                             tmp22(ii) = tmp21(ii) * r(ii,ir);  
%                         end
%   
%                     end
%                     if size(tmp22,1)>0
%                         constraints = [constraints, R(ir)-sum(tmp22) >= 0];
%                     end

                    tmp2 = 0;    
                    for ii = 1:I
                        ub = min(it+d(ii)-1,H);
                        lb = max(it,d(ii));
                        tmp2 = tmp2 + sum(x(ii,lb:ub)) * r(ii,ir);                          
                    end
                    
                    constraints = [constraints, sum(tmp2) <= R(ir)];
                    
                end
            end

            % c3: precedence
            for ii = 1:I
                J = P(ii,:).*z; %predecesor task of task i
                J1 = nonzeros(J);
                if isempty(J1)~=0
                    tmp31 = t.*x(ii,:);
                    %fti = sum(tmp31(d(ii)+1:H+1)); % finish time of task i
                    fti = sum(tmp31(d(ii):H)); % finish time of task i
                    for k = 1:length(J1)
                        j = J1(k);
                        tmp32 = t.*x(j,:);
                        %ftj = sum(tmp32(d(j)+1:H+1));  % finish time of task j (predecesor of task i)
                        ftj = sum(tmp32(d(j):H));  % finish time of task j (predecesor of task i)
                        constraints=[constraints, fti-d(ii)-ftj >= 0];
                    end

                end
            end

           % =========================objective=============================
            % minimize the finishing time of the dummy end task
            objective = 0;
            for ii = 1:I
                tmpo = t.*x(ii,:);
                %ft(ii) = sum(tmpo(d(ii)+1:H+1));
                ft(ii) = sum(tmpo(d(ii):H));
            end
            objective = objective + max(ft); 

           % ==========================solve=================================
            ops = sdpsettings('solver','gurobi','verbose',0,'showprogress',1);
            optimize(constraints,objective,ops)
            opx = value(x); % value of variable x 
            obj = value(objective); % value of objective function

            % ==========================post-process=========================
            X = opx;     X(isnan(X))=0;
            XR = zeros(I,H);

            for ii = 1:I                
                id = find(X(ii,:)); 
                if isempty(id)
                    st(ii) = inf;
                    ft(ii) = inf;
                else
                    ft(ii) = min(id);
                    st(ii) = ft(ii) - d(ii)+1;

                    lb = st(ii);
                    ub = ft(ii);
                    XR(ii,lb:ub) = ones(1,d(ii));
                end
            end
%                
%             for ii = 1:I
%                 lb = st(ii);
%                 ub = ft(ii);
%                 XR(ii,lb:ub) = ones(1,d(ii));
%             end
% 
%             for it = 1:H %H+1
%                 ru(it,:) = XR(:,it)'*r; 
%             end
            if any(XR)
                for it = 1:H %H+1
                    ru(it,:) = XR(:,it)'*r; 
                end
            else 
                ru = zeros(H,K)+inf;
            end

            sft = [st',ft', z'];    % starting time (1st column) and finishing time (2nd column)
            ct = max(ft);       % completion time of all tasks
            
            tsort = sortrows(sft, 1); % sort rows based on the value in the first column
            sortTask = tsort(:,3);
                
            % =============== find the object ID and task ID ==============
            % objString = {'Generator','Bus','Branch','Powerline',...
            %                 'Bridge','Road','TrafficLight'...
            %                 'Antenna','Centraloffice','Cellline','CommunicationTower'}; 
            ss1 = task(:,3); ss2 = task(:,1);
            for ii = 1:length(ss1)
                index{ii} = find(isletter(ss1{ii}));
                ss12{ii} = ss1{ii}(1:max(index{ii}));
                ss13{ii} = ss1{ii}(1+max(index{ii}):end); 
            end
            
            for ii = 1:length(ss1)
                allTask{ii,1} = strcat(ss12{ii},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss13{ii});
                allTask{ii,1} = strcat(allTask{ii,1},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss2{ii});
            end
            
            for ii = 1:I % I-1
                id = sortTask(ii);
                sequenceTask{ii} = allTask{id};
                objectList{ii} = task{id,3}; 
            end
            sequenceObject = unique(objectList,'stable');

        end
        
        %% Scheme 3B: Optimization, RCPSP_ST (STEP 3)
        function [X,ct,sft,ru,sequenceTask,sequenceObject] = optimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon)
            % ================================================================
            % function optimizationRCPSP_ST
            % updated optimization algorithm 
            % originally developed by Yinan Liu
            % binary decision variable x_{it} = 1, when task i starts at time t;
            % x_{it} = 0 otherwise. 
            %  ================================================================
            %==== Input parameters
            % (1) task = {1-'TaskID in the system',...
            % 2-'Object Type',...
            % 3-'Object ID',...
            % 4-'TaskID in the task library of this Object type', ...
            % 5-'Decription of this task',...
            % 6-'Duration min', 7-'Duration max', 8-'Duration mode', 9-'Duration distribution type',...
            % 10-'Duration sample',...
            % 11-14'Resource demand for resource type 1, 2, 3, 4, respectively'}
            % All aforementioned number labels mean the column number in the
            % cell of "task"
            % (2) precedence, binary, = 1, if the task labeled as column
            % number is the precedessory of the task labeled as row number 
            % (3) resourceConstraint = ar = [AR1 AR2, AR3 AR4], avaiability
            % for 4 types of resource, constant over time
            % (4) time_horizon, the entire time horizon of interest
            %==== Output parameters
            % (1) X = opx matrix, size=I*th, 
            % every row represent the decision variable value over t=[1,2, .., th]
            % every column presente the decision variable value at a time snapshot
            % (2) ct = completion time of all tasks
            % (3) sft = [st,ft,taskID]
            % (4) sequenceTask = a ranking list of tasks in terms of restoration time
            % from early to late
            % (5) sequenceObject = a ranking list of damaged objects in terms of restoration time
            % from early to late
            % 
            %%=========================================================
            % d: restoration duration as the mode value in column 8 in task
            % H: time horizon
            % I: number of resource type
            % P: precedence matrix 
            % R: renewable resource constraint
            % r: resource demand as the values in column 11~14 in task
            % z: an intermediate variable as taskID in this function only 
            % from 1 to I
            
            %======================= parameters===========================
            H = time_horizon;            
            R = resourceConstraint;      
            K = size(R,1);               
            I = size(task,1);  
            P = precedence;   
            z = linspace(1,I,I); 
%             for ii = 1:I
%                 if isempty(d(ii))
%                     d(ii) = round(0.5*(task{ii,6}+task{ii,7}));
%                 else
%                     d(ii) = round(task{ii, 8});
%                 end
%                 for ir = 1:K                   
%                     r(ii,ir) = task{ii, 10+ir};
%                 end
%             end

                for ii = 1:I
%                     if exist('d(ii)','var')
%                         d(ii) = round(0.5*(task{ii,6}+task{ii,7}));
%                     else
%                         d(ii) = round(task{ii, 8});
%                     end
                    d(ii) = max(1,round(task{ii, 8}));
                    if isempty(d(ii)) || isnan(d(ii))
                        d(ii) = max(1,round(0.5*(task{ii,6}+task{ii,7})));
                    end
                    r(ii,:) = task{ii, 11};
                end
            % ======================== decision variable==================
            x=binvar(I,H); % x=1 if task i is started at time h
            
            % =========================constraints========================
            % c1: each task must be finished within the time horizon with no resource
            % waste
            constraints=[];
            for k=1:K
                for h=1:H
                    %temp7=transpose(x(:,h))*r{kk}(:,k);
                    %constraints=[constraints,temp7<=R{kk}(k)];
                    temp7=0;
                    for i=1:I
                        s=max(h-d(i)+1,1);
                        temp7=temp7+sum(x(i,s:h))*r(i,k);
                    end     
                    constraints=[constraints,temp7<=R(k)];
                end
            end
            
            % c2: precedence constraints
            for i=1:I
                J=P(i,:).*z;
                J1=nonzeros(J);
                if isempty(J1)~=1
                    for j=J1
                        temp=x(i,:)*transpose(1:H);
                        temp3=x(j,:)*transpose((d(j)+1):(d(j)+H));
                        constraints=[constraints,temp>=temp3];
                    end
                end       
            end
            
            % c3: there is only one starting time for each task during the horizon
            % needs
            for i=1:I
                temp4=sum(x,2);
                constraints=[constraints,temp4(i)==1];
            end

            % =========================objective=============================
            objective=0;
            for i=1:I
                temp6=x(i,1:H)*transpose(1:H);
                objective=objective+temp6;
            end

            % ==========================solve=================================
            ops = sdpsettings('solver','gurobi','verbose',0);
            optimize(constraints,objective,ops)
            opx=value(x); % value of variable x
            obj=value(objective); % value of objective function
            
            % ==========================post-process=========================
            
            X = opx;     X(isnan(X))=0;

            for ii = 1:I
                id = find(X(ii,:)); 
                if isempty(id)
                    st(ii) = inf;
                    ft(ii) = inf;
                else
                    st(ii) = min(id);
                    % ft(ii) = max(id);
                    ft(ii) = st(ii) + d(ii)-1;
                end
           end

           for it = 1:H
               ru(it,:) = X(:,it)'*r; 
           end

           sft = [st',ft', z'];
           if isinf(sft)
              ct = inf;
           else
              ct = max(ft);
           end
           tsort = sortrows(sft);
           sortTask = tsort(:,3);
            
           % =============== find the object ID and task ID ==============
           %             objString = {'Generator','Bus','Branch','Powerline',...
           %                 'Bridge','Road','TrafficLight'...
           %                 'Antenna','Centraloffice','Cellline','CommunicationTower'}; 
            ss1 = task(:,3); ss2 = task(:,1);
            for ii = 1:length(ss1)
                index{ii} = find(isletter(ss1{ii}));
                ss12{ii} = ss1{ii}(1:max(index{ii}));
                ss13{ii} = ss1{ii}(1+max(index{ii}):end); 
            end
            
            for ii = 1:length(ss1)
                allTask{ii,1} = strcat(ss12{ii},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss13{ii});
                allTask{ii,1} = strcat(allTask{ii,1},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss2{ii});
            end
            
            for ii = 1:I
                id = sortTask(ii);
                sequenceTask{ii} = allTask{id};
                objectList{ii} = task{id,3}; 
            end
            sequenceObject = unique(objectList,'stable');
            
        end

        %% Scheme 3B: Optimization, Flexible (STEP 3)
        function [X,ct,sft,ru,sequenceTask,sequenceObject] = optimizationFlexible(task, precedence, resourceConstraint, time_horizon)
            % ================================================================
            % function optimizationFlexible
            % updated optimization algorithm allowing flexible duration and
            % varying resource amount in terms of "the set of resource"
            % originally developed by Yinan Liu
            % integer decision variable x_{it} = the set of resources assigned
            % to task i at time t. (>=0)
            % binary variable a_{it}, applied to transfer if-then constraint 
            %  ================================================================
            %==== Input parameters
            % (1) task = {1-'TaskID in the system',...
            % 2-'Object Type',...
            % 3-'Object ID',...
            % 4-'TaskID in the task library of this Object type', ...
            % 5-'Decription of this task',...
            % 6-'Duration min', 7-'Duration max', 8-'Duration mode', 9-'Duration distribution type',...
            % 10-'Duration sample',...
            % 11-14'Resource demand for resource type 1, 2, 3, 4, respectively'}
            % All aforementioned number labels mean the column number in the
            % cell of "task"
            % (2) precedence, binary, = 1, if the task labeled as column
            % number is the precedessory of the task labeled as row number 
            % (3) resourceConstraint = ar = [AR1 AR2, AR3 AR4], avaiability
            % for 4 types of resource, constant over time
            % (4) time_horizon, the entire time horizon of interest
            %==== Output parameters
            % (1) X = opx matrix, size=I*th, 
            % every row represent the decision variable value over t=[1,2, .., th]
            % every column presente the decision variable value at a time snapshot
            % (2) ct = completion time of all tasks
            % (3) sft = [st,ft,taskID]
            % (4) sequenceTask = a ranking list of tasks in terms of restoration time
            % from early to late
            % (5) sequenceObject = a ranking list of damaged objects in terms of restoration time
            % from early to late
            % 
            %%=========================================================
            % d: restoration duration as the mode value in column 8 in task
            % H: time horizon
            % I: number of resource type
            % P: precedence matrix 
            % R: renewable resource constraint
            % r: resource demand as the values in column 11~14 in task
            % z: an intermediate variable as taskID in this function only 
            % from 1 to I
            
            %======================= input parameters===========================
                H = time_horizon;               % time horizon
                K = length(resourceConstraint); % the number of resource type
                I = size(task,1);               % the number of tasks
                R = resourceConstraint;         % resource constraints
                P = precedence;   % precedence

                M = 1e10; % A big enough number for if-then constraint


                t = linspace(1,H,H);

                % assign duration d and resource requirement r for every task 
                % during schedule, task duration d = duration mode
                for ii = 1:I
%                     if exist('d(ii)','var')
%                         d(ii) = round(0.5*(task{ii,6}+task{ii,7}));
%                     else
%                         d(ii) = round(task{ii, 8});
%                     end
                    d(ii) = max(1, round(task{ii, 8}));
                    if isempty(d(ii)) || isnan(d(ii))
                        d(ii) = max(1, round(0.5*(task{ii,6}+task{ii,7})));
                    end
                    r(ii,:) = task{ii, 11};
                end

             %======================= add a dummy end task ======================   
            I = I+1;
            P = [P; ones(1,I-1)]; P = [P,zeros(I,1)]; 
            d = [d,0];
            r = [r; zeros(1,K)];

            z = linspace(1,I,I); 

            % ==========================variable=================================
            x=intvar(I,H); % the number of type k resources that should be sent to task i in time horizon h
            a=binvar(I,H); % applied to transfer if-then constraint

            % =========================constraints================================
            % c1: at any period, the resource used by all tasks cannot exceed the total
            % number of resource available
            constraints=[];
            for k=1:K
                for h=1:H
                    temp7=transpose(x(:,h))*r(:,k);
                    constraints=[constraints,temp7<=R(k)];
                end
            end

            % c2: precedence constraints
            for i=1:I
                J=P(i,:).*z;
                J1=nonzeros(J);
                if isempty(J1)~=0
                    for h=2:H
                        temp=M*a(i,h);
                        temp2=sum(x(:,1:h-1),2);
                        for j=J1
                            temp=temp-t(j)-temp2(j); % check -temp2(j) or +temp2(j)
                        end
                        constraints=[constraints,temp>=0];
                    end
                    temp3=M*a(i,1);
                    for j=J1
                        temp3=temp3-t(j);
                    end
                    constraints=[constraints,temp3>=0];
                    for i=1:I
                        for h=1:H
                            constraints=[constraints,M*(1-a(i,h))-x(i,h)>=0];
                        end
                    end
                end
            end

            % c3: we don't want to assign extra resource to a task than it actually
            % needs
            for i=1:I
                temp4=sum(x,2);
                constraints=[constraints,d(i)-temp4(i)>=0];
            end

            % c4: lowerbound for x
            for i=1:I
                for h=1:H
                    constraints=[constraints,x(i,h)>=0];
                end
            end
            % =========================objective=============================
            objective=0;
            for i=1:I
                temp6=sum(x,2);
                objective=objective+d(i)-temp6(i);
            end

            % ==========================solve=================================
            ops = sdpsettings('solver','gurobi','gurobi.IterationLimit',100,'verbose',0,'showprogress',1);
            optimize(constraints,objective,ops)
            opx=value(x); % value of variable x
            opa=value(a);
            obj=value(objective); % value of objective function

            % ==========================post-process=========================
                X = opx;     X(isnan(X))=0;

                for ii = 1:I
                    id = find(X(ii,:)); 
                    if isempty(id)
                        st(ii) = inf;
                        ft(ii) = inf;
                    else
                        st(ii) = min(id);
                        ft(ii) = max(id);
                    end

                end

                for it = 1:H
                    ru(it,:) = X(:,it)'*r; 
                end

                sft = [st',ft', z(1:end)'];
				sft = sft(1:end-1, :);  % get rid of the dummy end task
                if isinf(sft)
                    ct = inf;
                else
                    ct = max(ft);
                end
                tsort = sortrows(sft);
                sortTask = tsort(:,3);
                
            % =============== find the object ID and task ID ==============
            %             objString = {'Generator','Bus','Branch','Powerline',...
            %                 'Bridge','Road','TrafficLight'...
            %                 'Antenna','Centraloffice','Cellline','CommunicationTower'}; 
            ss1 = task(:,3); ss2 = task(:,1);
            for ii = 1:length(ss1)
                index{ii} = find(isletter(ss1{ii}));
                ss12{ii} = ss1{ii}(1:max(index{ii}));
                ss13{ii} = ss1{ii}(1+max(index{ii}):end); 
            end
            
            for ii = 1:length(ss1)
                allTask{ii,1} = strcat(ss12{ii},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss13{ii});
                allTask{ii,1} = strcat(allTask{ii,1},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss2{ii});
            end
            
            for ii = 1:I-1 % get rid of the dummy end task 
                id = sortTask(ii);
                sequenceTask{ii} = allTask{id};
                objectList{ii} = task{id,3}; 
            end
            sequenceObject = unique(objectList,'stable');
             
        end

 
        %% Scheme 3B: SubOptimization, RCPSP_FT  (STEP 3)       
        function [X,ct,sft,ru,sequenceTask,sequenceObject] = SuboptimizationRCPSP_FT(task, precedence, resourceConstraint, time_horizon)
            %% Input parameters
            % (1) task = {1-'TaskID in the system',...
            % 2-'Object Type',...
            % 3-'Object ID',...
            % 4-'TaskID in the task library of this Object type', ...
            % 5-'Decription of this task',...
            % 6-'Duration min', 7-'Duration max', 8-'Duration mode', 9-'Duration distribution type',...
            % 10-'Duration sample',...
            % 11-14'Resource demand for resource type 1, 2, 3, 4, respectively'}
            % All aforementioned number labels mean the column number in the
            % cell of "task"
            % (2) precedence, binary, = 1, if the task labeled as column
            % number is the precedessory of the task labeled as row number 
            % (3) resourceConstraint = ar = [AR1 AR2, AR3 AR4], avaiability
            % for 4 types of resource, constant over time
            % (4) time_horizon, the entire time horizon of interest
            %% Output parameters
            % (1) X = opx matrix, size=I*th, 
            % every row represent the decision variable value over t=[1,2, .., th]
            % every column presente the decision variable value at a time snapshot
            % (2) ct = completion time of all tasks
            % (3) sft = [st,ft,taskID]
            % (4) sequenceTask = a ranking list of tasks in terms of restoration time
            % from early to late
            % (5) sequenceObject = a ranking list of damaged objects in terms of restoration time
            % from early to late
            % 
            %%=========================================================
            % d: restoration duration as the mode value in column 8 in task
            % H: time horizon
            % I: number of resource type
            % P: precedence matrix 
            % R: renewable resource constraint
            % r: resource demand as the values in column 11~14 in task
            % z: an intermediate variable as taskID in this function only 
            % from 1 to I

            sequenceTask = [];
            sequenceObject = [];
            
            %======================= input parameters===========================
            H = time_horizon;               % time horizon
            K = length(resourceConstraint); % the number of resource type
            I0 = size(task,1);               % the number of tasks
            R = resourceConstraint;         % resource constraints
            P = precedence;                 % precedence

            %t = linspace(0,H,H+1);
            t = linspace(1,H,H);
            z0 = linspace(1,I0,I0); 

            % assign duration d and resource requirement r for every task 
            % during schedule, task duration d = duration mode
            for ii = 1:I0
                d(ii) = max(1, round(task{ii, 8}));  
                if isempty(d(ii)) || isnan(d(ii))
                    d(ii) = max(1, round(0.5*(task{ii,6}+task{ii,7})));
                end
                r(ii,:) = task{ii,11};                 
            end
            %X = zeros(I0,H+1);     
			X = zeros(I0,H);     
            track = [];
        while numel(track)~=I0
            % ========== revise 07/25 =========
            % generate the tasks calculated within the loop
            subtrack=[];
            for ii=1:I0
                if all(ismember(ii,track))==0
                    J=P(ii,:).*z0;
                    J1=nonzeros(J); % predecessor tasks of task ii
                    if all(ismember(J1,track))==1
                        subtrack=[subtrack,ii];
                    end
                end
            end
            % subtrack
            track = [track,subtrack];
            I = numel(subtrack);
            z = linspace(1,I,I); 
            % =================================
            
            % ====================== decision variables===========================
            %x = binvar(I,H+1); % which time step task i(=1,2,...,I,I+1) is finished within [0, H]
			x = binvar(I,H); % which time step task i(=1,2,...,I,I+1) is finished within [1, H]


            % =========================constraints===============================
            constraints = [];

            % c1: every task only executes once. 
            for ii = 1:I
                %constraints = [constraints, sum(x(ii,d(subtrack(ii))+1:H+1)) <=1];
				constraints = [constraints, sum(x(ii,d(subtrack(ii)):H)) <=1];				
            end


            % c2: resource usage should be less than resource constraint per time step         
            for ir = 1:K
                for it = 1:H % H+1
%                     tmp21=zeros(1,I);
%                     tmp22=zeros(1,I);
%                     for ii = 1:I
%                         % ub = min(it+d(subtrack(ii))-1,H+1)
% 						ub = min(it+d(subtrack(ii))-1,H)
%                         lb = max(it,d(subtrack(ii)))
%                         if ub-lb>=0
%                             tmp21(ii) = sum(x(ii,lb:ub));
%                             tmp22(ii) = tmp21(ii) * r(subtrack(ii),ir);  
%                         end
%   
%                     end
%                     if size(tmp22,1)>0
%                       constraints = [constraints, R(ir)-sum(tmp22) >= 0];
%                     end
                   
                    tmp2 = 0;    
                    for ii = 1:I
                        ub = min(it+d(ii)-1,H);
                        lb = max(it,d(ii));
                        tmp2 = tmp2 + sum(x(ii,lb:ub)) * r(ii,ir);                          
                    end
                    
                    constraints = [constraints, sum(tmp2) <= R(ir)];     
                end
            end


           % =========================objective=============================
            % minimize the finishing time of the dummy end task
            objective = 0;
            for ii = 1:I
                tmpo = t.*x(ii,:);
                %ft(ii) = sum(tmpo(d(subtrack(ii))+1:H+1));
				ft(ii) = sum(tmpo(d(subtrack(ii)):H));
            end
            objective = objective + max(ft); 

           % ==========================solve=================================
            ops = sdpsettings('solver','gurobi');
            optimize(constraints,objective,ops)
            opx = value(x); % value of variable x 
            obj = value(objective); % value of objective function
            % ======= revise 07/25 ============
            % post-process
            for ii=1:numel(subtrack)
                X(subtrack(ii),:)=opx(ii,:);
            end
            % =================================
        end
            % ==========================post-process=========================
            % X = opx;     
            %X(isnan(X))=0;
% 
%             for ii = 1:I0
%                 ft(ii) = t*X(ii,:)'; % finishing time of every task
%                 %st(ii) = ft(ii)-d(ii);% starting time of every task
% 				st(ii) = ft(ii)-d(ii)+1;% starting time of every task
%             end
% 
%             %XR = zeros(I0,H+1);
% 			XR = zeros(I0,H);
%             for ii = 1:I0
%                 %lb = ft(ii)-d(ii)+2;
%                 %ub = ft(ii)+1;
%                 lb = st(ii);
%                 ub = ft(ii);
%                 XR(ii,lb:ub) = ones(1,d(ii));
%             end
% 
%             for it = 1:H %H+1
%                 ru(it,:) = XR(:,it)'*r; 
%             end

            X(isnan(X))=0; XR = zeros(I,H);

            for ii = 1:I0                
                id = find(X(ii,:)); 
                if isempty(id)
                    st(ii) = inf;
                    ft(ii) = inf;
                else
                    ft(ii) = min(id);
                    st(ii) = ft(ii) - d(ii)+1;

                    lb = st(ii);
                    ub = ft(ii);
                    XR(ii,lb:ub) = ones(1,d(ii));
                end
            end
               
            if any(XR)
                for it = 1:H %H+1
                    ru(it,:) = XR(:,it)'*r; 
                end
            else 
                ru = zeros(H,K)+inf;
            end


            sft = [st',ft', z0'];    % starting time (1st column) and finishing time (2nd column)
            ct = max(ft);       % completion time of all tasks
            
            tsort = sortrows(sft);
            sortTask = tsort(:,3);
                
            % =============== find the object ID and task ID ==============
            %   objString = {'Generator','Bus','Branch','Powerline',...
            %                 'Bridge','Road','TrafficLight'...
            %                 'Antenna','Centraloffice','Cellline','CommunicationTower'}; 
            ss1 = task(:,3); ss2 = task(:,1);
            for ii = 1:length(ss1)
                index{ii} = find(isletter(ss1{ii}));
                ss12{ii} = ss1{ii}(1:max(index{ii}));
                ss13{ii} = ss1{ii}(1+max(index{ii}):end); 
            end
            
            for ii = 1:length(ss1)
                allTask{ii,1} = strcat(ss12{ii},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss13{ii});
                allTask{ii,1} = strcat(allTask{ii,1},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss2{ii});
            end
            
            for ii = 1:I0 %I-1
                id = sortTask(ii);
                sequenceTask{ii} = allTask{id};
                objectList{ii} = task{id,3}; 
            end
            sequenceObject = unique(objectList,'stable');
            
        end
 
        %% Scheme 3B: SubOptimization, RCPSP_ST  (STEP 3)   
        function [X,ct,sft,ru,sequenceTask,sequenceObject] = SuboptimizationRCPSP_ST(task, precedence, resourceConstraint, time_horizon);
            %% Input parameters
            % (1) task = {1-'TaskID in the system',...
            % 2-'Object Type',...
            % 3-'Object ID',...
            % 4-'TaskID in the task library of this Object type', ...
            % 5-'Decription of this task',...
            % 6-'Duration min', 7-'Duration max', 8-'Duration mode', 9-'Duration distribution type',...
            % 10-'Duration sample',...
            % 11-14'Resource demand for resource type 1, 2, 3, 4, respectively'}
            % All aforementioned number labels mean the column number in the
            % cell of "task"
            % (2) precedence, binary, = 1, if the task labeled as column
            % number is the precedessory of the task labeled as row number 
            % (3) resourceConstraint = ar = [AR1 AR2, AR3 AR4], avaiability
            % for 4 types of resource, constant over time
            % (4) time_horizon, the entire time horizon of interest
            %% Output parameters
            % (1) X = opx matrix, size=I*th, 
            % every row represent the decision variable value over t=[1,2, .., th]
            % every column presente the decision variable value at a time snapshot
            % (2) ct = completion time of all tasks
            % (3) sft = [st,ft,taskID]
            % (4) sequenceTask = a ranking list of tasks in terms of restoration time
            % from early to late
            % (5) sequenceObject = a ranking list of damaged objects in terms of restoration time
            % from early to late
            % 
            %%=========================================================
            % d: restoration duration as the mode value in column 8 in task
            % H: time horizon
            % I: number of resource type
            % P: precedence matrix 
            % R: renewable resource constraint
            % r: resource demand as the values in column 11~14 in task
            % z: an intermediate variable as taskID in this function only 
            % from 1 to I

            sequenceTask = [];
            sequenceObject = [];

            %======================= parameters===========================
            H = time_horizon;               % time horizon
            K = length(resourceConstraint); % the number of resource type
            I0 = size(task,1);               % the number of tasks
            R = resourceConstraint;         % resource constraints
            P = precedence;                 % precedence

            %t = linspace(0,H,H+1);
			t = linspace(1,H,H);
            z0 = linspace(1,I0,I0); 
            for ii = 1:I0
                d(ii) = max(1, round(task{ii, 8})); 
                if isempty(d(ii)) || isnan(d(ii))
                    d(ii) = max(1, round(0.5*(task{ii,6}+task{ii,7})));
                end
                r(ii,:) = task{ii,11};
            end
            X=zeros(I0,H);        
            track=[];
            while numel(track)~=I0
                % generate the tasks calculated within the loop
                subtrack=[];
                for ii=1:I0
                    if all(ismember(ii,track))==0
                        J=P(ii,:).*z0;
                        J1=nonzeros(J); % predecessor tasks of task ii
                        if all(ismember(J1,track))==1
                            subtrack=[subtrack,ii];
                        end
                    end
                end
                track=[track,subtrack];
                I=numel(subtrack);
                z = linspace(1,I,I);
                % ======================== decision variable==================
                x=binvar(I,H); % x=1 if task i is started at time h

                % =========================constraints========================
                % c1: each task must be finished within the time horizon with no resource
                % waste
                constraints=[];
                for k=1:K
                    for h=1:H
                        %temp7=transpose(x(:,h))*r{kk}(:,k);
                        %constraints=[constraints,temp7<=R{kk}(k)];
                        temp7=0;
                        for i=1:I
                            s=max(h-d(i)+1,1);
                            temp7=temp7+sum(x(i,s:h))*r(i,k);
                        end     
                        constraints=[constraints,temp7<=R(k)];
                    end
                end

                % c3: there is only one starting time for each task during the horizon
                % needs
                for i=1:I
                    temp4=sum(x,2);
                    constraints=[constraints,temp4(i)==1];
                end

                % =========================objective=============================
                objective=0;
                for i=1:I
                    temp6=x(i,1:H)*transpose(1:H);
                    objective=objective+temp6;
                end

                % ==========================solve=================================
                ops = sdpsettings('solver','gurobi');
                optimize(constraints,objective,ops)
                opx=value(x); % value of variable x
                obj=value(objective); % value of objective function
                % ======= revise 07/25 ============
                % post-process
                for ii=1:numel(subtrack)
                    X(subtrack(ii),:)=opx(ii,:);
                end
                % =================================
            end
            % ==========================post-process=========================

            X(isnan(X))=0;

            for ii = 1:I0
                id = find(X(ii,:)); 
                if isempty(id)
                    st(ii) = inf;
                    ft(ii) = inf;
                else
                    st(ii) = min(id);
                    ft(ii) = st(ii) + d(ii)-1;
                    %ft(ii) = max(id);
                end
           end

           for it = 1:H
               ru(it,:) = X(:,it)'*r; 
           end

           sft = [st',ft', z0'];
           if isinf(sft)
              ct = inf;
           else
              ct = max(ft);
           end
            
            tsort = sortrows(sft, 1); % sort the sft vector based on the 1st column
            sortTask = tsort(:,3);
                
            % =============== find the object ID and task ID ==============
            % objString = {'Generator','Bus','Branch','Powerline',...
            %                 'Bridge','Road','TrafficLight'...
            %                 'Antenna','Centraloffice','Cellline','CommunicationTower'}; 
            ss1 = task(:,3); ss2 = task(:,1);
            for ii = 1:length(ss1)
                index{ii} = find(isletter(ss1{ii}));
                ss12{ii} = ss1{ii}(1:max(index{ii}));
                ss13{ii} = ss1{ii}(1+max(index{ii}):end); 
            end
            
            for ii = 1:length(ss1)
                allTask{ii,1} = strcat(ss12{ii},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss13{ii});
                allTask{ii,1} = strcat(allTask{ii,1},'/');
                allTask{ii,1} = strcat(allTask{ii,1},ss2{ii});
            end
            
            for ii = 1:I0 %I-1
                id = sortTask(ii);
                sequenceTask{ii} = allTask{id};
                objectList{ii} = task{id,3}; 
            end
            sequenceObject = unique(objectList,'stable');

        end
 
       
        
        
        
        %% Count and get unfixed item from schedule
        function [set, pow_count, comm_count, trans_count] = countSchedule(orig_Schedule)
            % ================================================================
            % function countSchedule
            % Count and get unfixed item from schedule
            % ================================================================
            Power = orig_Schedule{1};
            Comm = orig_Schedule{2};
            Trans = orig_Schedule{3};
            
            power_count = 0;
            comm_count = 0;
            trans_count = 0;
            set = {};
            
            for i = 1:length(Power)
                tem = strsplit(Power{i},'/');
                if length(tem) == 2
                    power_count = power_count + 1;
                    set = [set, Power{i}];
                end
            end
            
            for i = 1:length(Comm)
                tem = strsplit(Comm{i},'/');
                if length(tem) == 2
                    comm_count = comm_count + 1;
                    set = [set, Comm{i}];
                end
            end
            
            for i = 1:length(Trans)
                tem = strsplit(Trans{i},'/');
                if length(tem) == 2 || length(tem) == 3
                    trans_count = trans_count + 1;
                    set = [set, Trans{i}];
                end
            end
        end
        
        %% Simulate the actual recovery process (STEP 4)
        % function Recovery
        % This function carries out the following items.
        % 1. sampling actual task duration by mapping from original samples following uniform distributions
        % 2. based on the resource avaiability and task resource demand, add
        % individual tasks to the current working list, considering
        % precedence constraints
        % 3. remove the working tasks from the remaining task list and subtract
        % resource demand from the available resource at the moment
        % 4. after completing any task, update object status and update the
        % graph network for every system
        % 5. compute the functionality of every system
        function [QSys, QNbr, trackt, trackCW, trackOS, lookupTable] = Recovery(ActiveSystem, time_horizon, Interdependence_Num, Qtrans0, InterdependenceFunc, ReSchedule_Num,...
                RestorationResource, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, ...
                RepairSchedule, Total, Graph, System_Dependent_Factor,  ...
                Seperate_Scheduling, LinkDirection, Per_New_reschedule, Num_stage, ...
                TableTask,TablePrecedence, OptimizationChoice, Diff_unit, Cust_unit)
            %================================================================
            % function Recovery
            % This function simulates the Recovery Process of all systems.
            % This function calls the following major functions: 
            % 1. in the file of Library: AddCurrentWorking, FindMinDays, 
            % UpdateCurrentDuration, SysFuncInterdependence, UpdateStatus
            % 2. in the file of Interface1: Functionality, RepairReSchedule
            %----% This function carries out the following computations.
            % 1. sampling actual task duration by mapping from original samples following uniform distributions
            % 2. based on the resource avaiability and task resource demand, add
            % individual tasks to the current working list, considering
            % precedence constraints
            % 3. remove the working tasks from the remaining task list and subtract
            % resource demand from the available resource at the moment
            % 4. after completing any task, update object status and update the
            % graph network for every system
            % 5. compute the functionality evolution of every system.
            % 
            %================================================================
            %==== Field
            %--- lookupTable: a table that documents all data about all tasks, 
            %such as task description, duration distribution, resource demand, 
            %actually sampled task druation, task uniqueID, etc.
            lookupTable = {};
            Power_Set = Total{1}; 
            Communication_Set = Total{2}; 
            Transportation_Set = Total{3};
            Dictionary = Total{4};
            Neighborhood = Total{5};
            powerGraph = Graph{1};
            commGraph = Graph{2};
            transGraph = Graph{3};
            active_power = ActiveSystem(1);          
            active_comm = ActiveSystem(2);  
            active_trans = ActiveSystem(3);  
            
            
            %--- Set up the initial functionality of every system based on
            %the percentage of population with a certain service at the
            %neighborhood level. 
            [totalPopulation, QNbrPower, QNbrComm, QNbrTrans] = Library.neighbourFunc(Dictionary);
            
            %--- Set up initial variables about the system functionality
            QNbr = [QNbrPower; QNbrComm; QNbrTrans];
            % Prior the hazard scenario, all system are fully functional. 
            QPower = {1,1,1,1,1,1,1,1,1,1};
            QComm = {1,1,1,1,1,1,1,1,1,1};
            QTrans = {1,1,1,1,1,1,1,1,1,1};
            
            %--- Set up a dictionary named funcTable that documents the
            % functionality of every object and every neighborhood if
            % the object/neighbood is the status of "Damaged" or "Stoped". 
            funcTable = containers.Map('KeyType','char','ValueType','any');
            funcTable = Library.CreateDictionaryFunctionality(Power_Set, Communication_Set, Transportation_Set, funcTable, time_horizon, Neighborhood);
            
          
            %--- Set up the working list for every system
            %This list evolves at different time step. 
            % A task can be added to the list if this task satisfies the
            % constraint of resource and precedence at the moment. 
            CurrentWorking_Power = {};
            CurrentWorking_Comm = {};
            CurrentWorking_Trans = {};
                     
            %--- Set up initial values of supporting variables
            %- finish: a flag variable representing whether the restoration
            % is finished or not. (0-No/1-Yes)
            %- Start_Day: time index representing the starting time of a task
            %- End_Day: time index representing the ending time of a task
            %- flagDelay: the flag variable for tracking the transprotatino delaying effect.
            % initially flagDelay = 0, if the delaying effect is turned on,
            % finish = 0 and Qtrans<Qtrans0, flagDelay = 1; then finish = 0
            % and Qtrans>=Qtrans0, flagDelay = 2; else flagDelay = 0. 
            %- need_reschedule: flag variables representing whether a system
            % needs to be rescheduled. 
            % need_reschedule = [need_reschedulePower, need_rescheduleComm, need_rescheduleTrans]
            finish = 0;
            Start_Day = 1;
            End_Day = 1;
            flagDelay = 0;

            need_reschedule = zeros(1,3); 
            
            %--- for the purpose of debugging, track the variable
            %evoluations while the time index jumps 
            indtrack = 1;
            %--- trackt: for the purpose of debugging, a vector of time index 
            % to track which time steps has the time index jumps over in 
            % the entire recovery process.
            trackt(indtrack,:) = [Start_Day, End_Day]; 
            %--- trackCW: for the purpose of debugging, track the variable
            %of CurrentWorking while the time index jumps  
            trackCW{indtrack,1} = CurrentWorking_Power;
            trackCW{indtrack,2} = CurrentWorking_Comm;
            trackCW{indtrack,3} = CurrentWorking_Trans;
            
     
            %--- set up resource constraints for every system
            Max_Power = RestorationResource(1,:);
            Max_Comm = RestorationResource(2,:); 
            Max_Trans = RestorationResource(3,:);
            
            %--- set up the restoration plan for every system (determined from 
            % Step 3, either scheme 3A-priority or scheme 3B-optimization)
            Schedule_Power = RepairSchedule{1};
            Schedule_Comm = RepairSchedule{2};
            Schedule_Trans = RepairSchedule{3};
            
            %--- set up dictionary of task and dictionary of precedence for
            %every system
            Dic_p_task = TableTask{2,1}; 
            Dic_c_task = TableTask{2,2};
            Dic_t_task = TableTask{2,3};
            Dic_p_prece = TablePrecedence{2,1};
            Dic_c_prece = TablePrecedence{2,2};
            Dic_t_prece = TablePrecedence{2,3};

            %--- Other variables
            % total_damaged: a variable couting the total number of tasks to repair damaged components of all three systems in a sample of damage scenario.
            % total_fixed: a variable couting how many tasks have been fixed. 
            total_damaged = length(Schedule_Power) + length(Schedule_Comm) + length(Schedule_Trans);
            total_fixed = 0;
            
            trackOS{1}{indtrack,1} = {}; %Power_Set;
            trackOS{2}{indtrack,1} = {}; %Communication_Set;
            trackOS{3}{indtrack,1} = {}; %Transportation_Set; 
            
            while finish == 0 % meaning there is (are) still some task(s) left.
                
                %==== Original schedule
                ori_schedule.Schedule_Power=Schedule_Power;
                ori_schedule.Schedule_Comm=Schedule_Comm;
                ori_schedule.Schedule_Trans=Schedule_Trans;
                %==== Check if the variable of Remain_schedule exists
                if ~exist('Remain_schedule','var')
                    Remain_schedule=ori_schedule;
                end
  
                %==== Add tasks for reparing any Damaged Component to the 
                % Current Working List Based on the Resource Constraints: 
                % Max_Power, Max_Comm, Max_Trans
                
                %--- time index for debuging
                func_index = Start_Day;
                %ts = Start_Day 
                %te = End_Day
                indtrack = indtrack;
                CurrentTime = Start_Day;
                RestorationPhase = Cust_unit; % a horizontal vector of different time instants in the time unit in "day".
                
                %==== Check whether different time unit is turned on/off 
                % n = the index representing which time unit is used at the current time step.
                % n = either one of (1, 2, 3, 4), 1-hour, 2-day, 3-week, 4-month
                
                if Diff_unit %== 1 
                    RestorationTimePhase = sort([CurrentTime, RestorationPhase]);
                    [~,n] = find(RestorationTimePhase == CurrentTime); n=4;

                    if gt(n,length(RestorationPhase)+1) % n>4 if length(RestorationPhase) = 3 [e.g., RestorationPhase = [3,28,168]];
                        msg = strcat('Function Recovery Error: Time Index of CurrentTime (n) = ', num2str(n), ' > length(RestorationPhase)+1 for ', tasktemp);
                        disp(msg); 
                        return               
                    end
                        
                % Diff_unit is not 1, meaning the choice of using different time unit in recovery is turned off.  
                % All tasks using the same time unit of "day"
                else
                    n = 2; 
                end
%               End_Day
                %==== Add restoration tasks to the current working list for every system  
                %==== Debugging-related to WorkingDays  
                if active_power
                    [CurrentWorking_Power, Schedule_Power, Max_Power, lookupTable] = Library.AddCurrentWorking(Max_Power, CurrentWorking_Power, Schedule_Power, Dictionary, lookupTable, Start_Day, n, indtrack);
                    [DurationP{indtrack}, RdemandP{indtrack}] = Library.GetTaskDuration(Dictionary, CurrentWorking_Power);
                end
                if active_comm
                    [CurrentWorking_Comm, Schedule_Comm, Max_Comm, lookupTable] = Library.AddCurrentWorking(Max_Comm, CurrentWorking_Comm, Schedule_Comm, Dictionary, lookupTable, Start_Day, n, indtrack);
                    [DurationC{indtrack}, RdemandC{indtrack}] = Library.GetTaskDuration(Dictionary, CurrentWorking_Comm);
                end
                if active_trans
                    [CurrentWorking_Trans, Schedule_Trans, Max_Trans, lookupTable] = Library.AddCurrentWorking(Max_Trans, CurrentWorking_Trans, Schedule_Trans, Dictionary, lookupTable, Start_Day, n, indtrack);
                    [DurationT{indtrack}, RdemandT{indtrack}] = Library.GetTaskDuration(Dictionary, CurrentWorking_Trans);
                end
                       
                %==== Find the shortest duration that a task takes among all tasks in the list of CurrentWorking
                % Find Day (the shorest time of a task takes to
                % complete among all tasks in the list of CurrentWorking)
                % so that the time index can jump from the current time
                % step to the new time step, when this task with the shortest 
                % duration in the CurrentWorking list finishes. 
                % Days in the unit of day
                Days = Library.FindMinDays(CurrentWorking_Power, CurrentWorking_Comm, CurrentWorking_Trans, Power_Set, Communication_Set, Transportation_Set,Dictionary);

                if Days > 0 % revised expression 
                    Start_Day = End_Day;    
                    End_Day = End_Day + Days;
                    
                    indtrack = indtrack+1;
                    trackt(indtrack,:) = [Start_Day, End_Day];
                    trackCW{indtrack,1} = CurrentWorking_Power;
                    trackCW{indtrack,2} = CurrentWorking_Comm;
                    trackCW{indtrack,3} = CurrentWorking_Trans;
                else
                    End_Day = time_horizon;
                    finish = 1;
                    
                    indtrack = indtrack+1;
                    trackt(indtrack,:) = [Start_Day, End_Day];
                    trackCW{indtrack,1} = CurrentWorking_Power;
                    trackCW{indtrack,2} = CurrentWorking_Comm;
                    trackCW{indtrack,3} = CurrentWorking_Trans;
                end
                
                %==== Updating the WorkingDays in the list of CurrentWorking 
                if active_power
					Library.UpdateCurrentDuration(CurrentWorking_Power, Dictionary, Days);
                end
                if active_comm
					Library.UpdateCurrentDuration(CurrentWorking_Comm, Dictionary, Days);
                end
                if active_trans
					Library.UpdateCurrentDuration(CurrentWorking_Trans, Dictionary, Days);
                end
              
                %==== Update the list of CurrentWorking for every system by
                %cleaning the task(s) that have been executed, and by
                %adding the amount of used resource due to executing task 
                %in the current list to Max_Power
                if active_power
					[CurrentWorking_Power,Max_Power] = Library.UpdateCurrentResource(CurrentWorking_Power, Dictionary, Max_Power);
                end
                if active_comm
					[CurrentWorking_Comm,Max_Comm] = Library.UpdateCurrentResource(CurrentWorking_Comm, Dictionary, Max_Comm);
                end
                if active_trans
					[CurrentWorking_Trans,Max_Trans] = Library.UpdateCurrentResource(CurrentWorking_Trans, Dictionary, Max_Trans);
                end
    
                %==== Update Object.Status and Functionality for repaired
                %Component/Object in every system
                [Remain_schedule,need_reschedule, total_fixed, transGraph, powerGraph, commGraph,...
                         Power_Set, Communication_Set, Transportation_Set, funcTable, FunctionalityTotal,...
                         totalPopulation, trackOS] = Library.UpdateStatus(ActiveSystem, Power_Set, Communication_Set, Transportation_Set,...
                         total_damaged, total_fixed, need_reschedule,Dictionary,...
                         transGraph, powerGraph, commGraph,funcTable,Start_Day,End_Day,...
                         Per_New_reschedule,Num_stage,Remain_schedule,ori_schedule,...
                         Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirection, InterdependenceFunc, indtrack, trackOS);
                
                %==== Update System functionality and Neighborhood functionality
                FunctionalityPower = FunctionalityTotal{1};
                FunctionalityComm = FunctionalityTotal{2};
                FunctionalityTrans = FunctionalityTotal{3};
                FuncPowerNbr = FunctionalityTotal{4};
                FuncCommNbr = FunctionalityTotal{5};
                FuncTransNbr = FunctionalityTotal{6};

                %==== Revised Expressions
                
                for ii = 1:length(FunctionalityPower) % the index of power functionality metric
                    if active_power
                        QPower{ii}(indtrack) = FunctionalityPower{ii};
                    else
                        QPower{ii}(indtrack) = 1;
                    end
                end
                for ii = 1:length(FunctionalityComm) % the index of communication functionality metric
                    if active_comm
                        QComm{ii}(indtrack) = FunctionalityComm{ii};
                    else
                        QComm{ii}(indtrack) = 1;
                    end
                end
                for ii = 1:length(FunctionalityTrans) % the index of transportation functionality metric
                    if active_trans
                        QTrans{ii}(indtrack) = FunctionalityTrans{ii};
                    else
                        QTrans{ii}(indtrack) = 1;
                    end
                end
                QSys{1} = QPower;
                QSys{2} = QComm;
                QSys{3} = QTrans;
                QNbr(1,indtrack) = FuncPowerNbr;
                QNbr(2,indtrack) = FuncCommNbr;
                QNbr(3,indtrack) = FuncTransNbr;
                    
                %==== Check transport delay effect due to the functionality
                % disruptions of the transportation system to decide
                % whether the duration of remainign tasks in the power
                % andcommunication systems need to be adjusted by
                % multipling a delaying factor > 1 (i.e., variable System_Dependent_Factor). 
                % (1) If Interdependence_Num == 1, the delaying effet is turned on.
                % Simply multiplying the delaying factor to the duration sample of
                % remaining tasks for the power and communication systems,
                % when Qtrans(t)<Qtrans0 (transportation functionality threshold). 
                % (2) If Interdependence_Num == 0, the delaying effet is turned off. 
                % Simply leave all the duration samples as it is for every remaining task for the power and communication systems. 
                if Interdependence_Num == 1 % delaying effect is turned on.
                    %FunctionalityTrans = FunctionalityTotal(3);
                    if flagDelay == 0 && lt(FunctionalityTrans{1}, Qtrans0)
                        Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, System_Dependent_Factor, Dictionary);
                        flagDelay = 1;
                    elseif flagDelay == 1 && ge(FunctionalityTrans{1}, Qtrans0)
                        Restore_Factor = 1 / System_Dependent_Factor;
                        Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, Restore_Factor, Dictionary);
                        flagDelay = 2;
                    else
                        Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, 1, Dictionary);
                        
                    end
                else   % Interdependence_Num == 0, delaying effect is turned off. 
                    Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, 1, Dictionary);
                end
                
                if ~(ReSchedule_Num==0)
                    %==== Check Reschedule effect              
                    % need_reschedule = [need_reschedulePower, need_rescheduleComm, need_rescheduleTrans]
                    need_reschedulePower = need_reschedule(1);
                    need_rescheduleComm = need_reschedule(2);
                    need_rescheduleTrans = need_reschedule(3);
                    if ~exist('Per_New_reschedule_Power','var')
                        Per_New_reschedule_Power=Per_New_reschedule;
                        Per_New_reschedule_Comm=Per_New_reschedule;
                        Per_New_reschedule_Trans=Per_New_reschedule;
                        Num_stage_Power=Num_stage;
                        Num_stage_Comm=Num_stage;
                        Num_stage_Trans=Num_stage;
                    end
                    [need_reschedulePower,Per_New_reschedule_Power,Num_stage_Power] = Library.checkNeedRescheduleNew(need_reschedulePower, FunctionalityPower, Start_Day, End_Day,Per_New_reschedule_Power,Num_stage_Power);
                    [need_rescheduleComm,Per_New_reschedule_Comm,Num_stage_Comm] = Library.checkNeedRescheduleNew(need_rescheduleComm, FunctionalityComm, Start_Day, End_Day,Per_New_reschedule_Comm,Num_stage_Comm);            
                    [need_rescheduleTrans,Per_New_reschedule_Trans,Num_stage_Trans] = Library.checkNeedRescheduleNew(need_rescheduleTrans, FunctionalityTrans, Start_Day, End_Day, Per_New_reschedule_Trans,Num_stage_Trans);
                    need_reschedule = [need_reschedulePower, need_rescheduleComm, need_rescheduleTrans];


                    checkReschedule = find(need_reschedule == 1); % find which system has need_reschedule=1
                    if ~isempty(checkReschedule)

                        % Get cells of remaining tasks and the precedence of remaining tasks
                        [remain_task_p,remain_task_c,remain_task_t,remain_prece_p,remain_prece_c,remain_prece_t]=Library.RemainTaskPrecedence(Dic_p_task,Dic_c_task,Dic_t_task,Dic_p_prece,Dic_c_prece,Dic_t_prece,Remain_schedule);
                        Remain_task={remain_task_p,remain_task_c,remain_task_t};
                        Remain_precedence={remain_prece_p,remain_prece_c,remain_prece_t};

                        resource = [Max_Power; Max_Comm; Max_Trans];

                        if need_reschedule(1) == 1 && active_power && ~isempty(remain_task_p)
                            disp('Reschedule Power')
                            System = 'Power'; 
                            [Power_ReSchedule, Power_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, Remain_task, Remain_precedence, resource, time_horizon);
                            Schedule_Power = Power_ReSchedule;
                        end
                        if need_reschedule(2) == 1 && active_comm && ~isempty(remain_task_c)
                            disp('Reschedule Communication')
                            System = 'Communication';
                            [Comm_ReSchedule, Comm_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, Remain_task, Remain_precedence, resource, time_horizon);
                            Schedule_Comm = Comm_ReSchedule;
                        end
                        if need_reschedule(3) == 1 && active_trans && ~isempty(remain_task_t)
                            disp('Reschedule Transportation')
                            System = 'Transportation';
                            [Trans_ReSchedule, Trans_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, Remain_task, Remain_precedence, resource, time_horizon);
                            Schedule_Trans = Trans_ReSchedule;
                        end               
                    end
                end
                need_reschedule=zeros(1,3);

                Start_Day = End_Day;

               
            end % while finish == 0
            
%             % Prior the hazard scenario, all system are fully functional. 
%             if ~active_power
%                 QPower = {1,1,1,1,1,1,1,1,1,1};
%             end
%             if ~active_comm
%                 QComm = {1,1,1,1,1,1,1,1,1,1};
%             end
%             if ~active_trans
%                 QTrans = {1,1,1,1,1,1,1,1,1,1};
%             end
            
            if QPower{1}(indtrack) < 1
                msg = strcat('Function Recovery Warning: QPower < 100% at indtrack =', num2str(indtrack),'[i.e., t = ', num2str(trackt(indtrack,2)), '].')
                disp(msg);
            end
           
            if QComm{1}(indtrack) < 1
                msg = strcat('Function Recovery Warning: QComm < 100% at indtrack =', num2str(indtrack),'[i.e., t = ', num2str(trackt(indtrack,2)), '].')
                disp(msg);
            end
            
            if QTrans{1}(indtrack) < 1
                msg = strcat('Function Recovery Warning: QTrans < 100% at indtrack =', num2str(indtrack),'[i.e., t = ', num2str(trackt(indtrack,2)), '].')
                disp(msg);
            end
            
        end
        
        %% Initiallize functionality table (func vs time for all objects) (STEP 4) 
        function funcTable = CreateDictionaryFunctionality(PowerSet, CommSet, TransSet, funcTable, time_horizon, Neighborhood)
            %===============================================================
            % function CreateDictionaryFunctionality
            % This function sets up the initial functionality table (a dictionary)
            % for three systems. 
            %===== Input
            % - PowerSet: the set of the Power system, which is a
            % combination of all object sets in the system. 
            % - CommSet: the set of the Communication system, which is a
            % combination of all object sets in the system. 
            % - TransSet: the set of the Transportation system, which is a
            % combination of all object sets in the system. 
            % - funcTable: the initial Dictionary of the object functionality 
            % - time_horizon
            % - Neighborhood: the set of Neighborhood
            %===== Output
            % - funcTable: Dictionary of the object functionality for a
            % system/neighborhood set.
            %===============================================================
            
            %=== set values about functionality Q(t) for every object in the dictionary named funcTable 
            funcTable = Library.initialFuncTableHelper(PowerSet,time_horizon, funcTable);
            funcTable = Library.initialFuncTableHelper(CommSet,time_horizon, funcTable);
            funcTable = Library.initialFuncTableHelper(TransSet,time_horizon, funcTable);
            
            %=== set values about functionality Q(t) for every neighborhood in the dictionary named funcTable 
            for ii = 1:length(Neighborhood)
                iQNbrPower = Neighborhood{ii}.PowerStatus;
                iQNbrComm = Neighborhood{ii}.CommStatus;
                iQNbrTrans = Neighborhood{ii}.TransStatus;
                funcTable(Neighborhood{ii}.uniqueID) = [iQNbrPower; iQNbrComm; iQNbrTrans];
            end
        end
                
        %% get Task.WorkingDays at every jumping time step in the Recovery (STEP 4) 
        function [Duration, Demand] = GetTaskDuration(Dictionary, Current)
            %=====================================================================
            % function GetTaskDuration
            % This function extracts the Task.WorkingDays at every jumping
            % time step in the Recovery process(SETP 4). It can be used in
            % the function of Library.Recovery, in combination with the functions 
            % related to task duration evolutions in case of time step jumping,
            % such as Library.AddCurrentWorking and Library.FindMinDays. 
            %=====================================================================
            ntask = size(Current,2);
            tmpDuration = zeros(ntask, 6+1); 
            tmpDemand = zeros(ntask, 4+1);
            
            Key = 'Task';
            
            for ii = 1:ntask
                if isempty(Current{1,ii})
                    continue;
                end
                tem = strsplit(Current{1,ii} ,'/');
                taskUID = tem{3}; % task unique ID in the Dictionary
                Index = strfind(taskUID, Key);
                taskID = sscanf(taskUID(Index(1) + length(Key):end), '%g', 1);
                
                task = Dictionary(taskUID);
                  
                tmpDuration(ii,:) = [task.WorkingDays, taskID];     
                if iscell(task.Resources)
                    resourceDemand = cell2mat(task.Resources);
                else
                    resourceDemand = task.Resources;
                end
                tmpDemand(ii,:) = [resourceDemand, taskID];   
            end
            
            Duration = tmpDuration;
            Demand = tmpDemand;
        
        end
        
        %% Update the Object.Status when a task is completed in the recovery process (STEP 4)  
        function [Remain_schedule,need_reschedule, total_fixed, transGraph, powerGraph, commGraph, Power, Comm, Trans, funcTable, FunctionalityTotal, totalPopulation, trackOS] = UpdateStatus(ActiveSystem, Power, Comm, Trans, total_damaged, total_fixed, need_reschedule, Dictionary, transGraph, powerGraph, commGraph, funcTable, Start_Day, End_Day, Per_New_reschedule, Num_stage, Remain_schedule, ori_schedule, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirection, InterdependenceFunc, indtrack, trackOS)
            % ============================================================================================
            % function UpdateStatus
            % called by the function Library.Recovery
            % This function updates the Object.Status after finishing a task as needed.
            % In particular, the following steps are performed. 
            % (1) The object status changes from "damaged" to "stoped"/"open",
            % or the object status changes from "stoped" to "open". 
            % (2) The corresponding Object.Functionality changes from 0 to 1. 
            % (3) The system functionality is updated.
            % (4) the graph of every system is updated by adding a node or an edge, 
            % depending on the Object type, when all tasks for repairing this damaged object are completed. 
            % =============================================================================================
            
            %==== Field
            %--- Power System
            Branch = Power{1};
            Bus = Power{2};
            Generator = Power{3};
            TransTower = Power{4};
            
            %--- Communication System
            Centraloffice = Comm{1};
            CommunicationTower = Comm{2};
            Cellline = Comm{3}; 
            
            %--- Transportation System
            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};
            
            %--- Flag variable to indicate whether 
            fixed = false;         
            
            %--- the active flag of every system: turn ON(1) or OFF (0)
            active_power = ActiveSystem(1);
            active_comm = ActiveSystem(2);
            active_trans = ActiveSystem(3);
     
            %--- flag variable of the rescheduling decision of every system: turn ON(1) or OFF (0) 
            need_reschedulePower = need_reschedule(1);
            need_rescheduleComm = need_reschedule(2);
            need_rescheduleTrans = need_reschedule(3);
            
            %--- trackOS: track Object.Status while the time index jumps 
            trackOS{1}{indtrack,1} = {}; %Power_Set;
            trackOS{2}{indtrack,1} = {}; %Communication_Set;
            trackOS{3}{indtrack,1} = {}; %Transportation_Set;
            
            %--- icol: the index that tracks the object status in a system
            icol = 1;
            
            %--- tolerance: the tolerance of task duration
            % If the task duration <=tolerance, considering this task duration = 0
            % and change the Object.Status if all tasks associated with
            % this Object have been executed and all functionality
            % dependencies are satisfied. 
            tolerance = 1e-6;
            
            %--- index_task: the index that tracks the number of tasks
            % have been executed at the moment. 
            index_task = 0;
            
            %==== Power System
            %--- Generator
            clearvars ObjSet; ObjSet = Generator;
            for ii = 1:length(ObjSet)     
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;
                if strcmp(iObj.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))          
                    iObj.Status = 'Open';
                    temp = funcTable(iObj.uniqueID);
                    temp(indtrack) = 1;
                    funcTable(iObj.uniqueID) = temp;
                    iObj.Functionality = 1;
                    total_fixed = total_fixed + 1;
                    index_task = index_task+1;
        
                    %--- track which Object changes status 
                    icol = icol + 1; 
                    trackOS{1}{indtrack,icol} = strcat('Generator', num2str(ii), ', Status changes from Damaged to Open');
                end
            end
            Generator = ObjSet;
            
            %--- Bus (Substation)
            clearvars ObjSet; ObjSet = Bus;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;
                
                %-- For substations (bus) that are connected to a power plant (generator)
                if ~isempty(Bus{ii}.Generator)
                    %-- If a substation is damaged and its task(s) finishes
                    % the restoration (workingdays <=0), change the status
                    % from "Damaged" to "Stoped". 
                    if strcmp(Bus{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
                        Bus{ii}.Status = 'Stoped';
                        total_fixed = total_fixed + 1;
                            
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{1}{indtrack,icol} = strcat('Bus',num2str(ii), ', Status changes from Damaged to Stoped');
                    end
                    
                    %-- If a substation is stoped and its connected generator
                    % (power plant) is open now, change the status of this
                    % substation from "Stoped" to "Open", and update its Functionality. 
                    if strcmp(Bus{ii}.Status, 'Stoped')
                        temp = extractAfter(Bus{ii}.Generator, 9);
                        temp = str2num(temp);
                        if strcmp(Generator{temp}.Status, 'Open')
                            iObj.Status = 'Open';
                            fixed = true;
                            temp = funcTable(Bus{ii}.uniqueID);
                            temp(indtrack) = 1;
                            funcTable(Bus{ii}.uniqueID) = temp;
                            iObj.Functionality = 1;
                            powerGraph = addnode(powerGraph, Bus{ii}.uniqueID);
                            
                            %--- track which Object changes status 
                            icol = icol + 1; 
                            trackOS{1}{indtrack,icol} = strcat('Bus',num2str(ii), ', Status changes from Stoped to Open');
                        end
                    end
                    
                %-- For substations (bus) that are not connected to a power plant (generator)    
                % If a substation is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Open". 
                else % isempty(Bus{ii}.Generator) = 1, meaning that there is no generator connecting to the bus
                    if strcmp(Bus{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
                        Bus{ii}.Status = 'Open';
                        Bus{ii}.Functionality = 1;
                            
                        temp = funcTable(Bus{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(Bus{ii}.uniqueID) = temp;
                            
                        powerGraph = addnode(powerGraph,Bus{ii}.uniqueID);
                        total_fixed = total_fixed + 1;
                            
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{1}{indtrack,icol} = strcat('Bus',num2str(ii), ', Status changes from Damaged to Open');

                    end
                end
                
                %-- If a substation is fixed (fixed=true, "Open"), update the 
                % power functionality property of the neighborhood (Neighborhood{j}.PowerStatus) 
                % that connects to this restored/fixed substation.  
                if fixed % gt(Bus{ii}.Functionality, 0) % meaning fixed = true, meaning strcmp(Bus{ii}.Status, 'Open')=1
                    for j = 1:length(Bus{ii}.Neighborhood)
                        temp = Dictionary(Bus{ii}.Neighborhood{j});
                        temp = temp{1};
                        %disp(temp.uniqueID);
                        temp.PowerStatus = 1;
                        temp2 = Dictionary(Bus{ii}.Neighborhood_Power_Link{j});
                        temp2 = temp2{1};
                        temp2.Status = 'Open'; % Power Status Open for the link (Branch) connecting to this Neighborhood
                        
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{1}{indtrack, icol} = strcat(temp.uniqueID, ': Nbr.Bus and Nbr.PowerLink, Status changes from Stoped to Open');
                        
                    end
                    % prepare the flag variable of fixed for the computation of the next Bus object by setting this flag variable to the initial value of "false". 
                    fixed = false;
                end
            end

            %--- Transmission Tower
            % If a transmission tower is damaged and its task(s) finishes
            % the restoration (workingdays <=0), change the status
            % from "Damaged" to "Open", and set the Object.Functionality= 1.
            clearvars ObjSet; ObjSet = TransTower;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii}; 
                iDuration = iObj.WorkingDays;
                if length(iDuration) == 1
                    iDuration = [24*iDuration, iDuration, inv(7)*iDuration, inv(28)*iDuration]; 
                end
                if strcmp(TransTower{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
                %if strcmp(TransTower{ii}.Status, 'Damaged') && any(lt(abs(iDuration), tolerance))
                   TransTower{ii}.Status = 'Open';
                   TransTower{ii}.Functionality = 1;
                   temp = funcTable(TransTower{ii}.uniqueID);
                   temp(indtrack) = 1;
                   funcTable(TransTower{ii}.uniqueID) = temp;       
                        
                   powerGraph = addnode(powerGraph, TransTower{ii}.uniqueID);
                   total_fixed = total_fixed + 1;   
                   %--- track which Object changes status 
                   icol = icol + 1; 
                   trackOS{1}{indtrack, icol} = strcat('TransTower',num2str(ii), ', from Damaged to Open');
                        
                end
            end
            
            %--- Branch (Power Line) 
            %-- If a power line is damaged and its task(s) finishes
            % the restoration (workingdays <=0), change the status
            % from "Damaged" to "Stoped". 
            clearvars ObjSet; ObjSet = Branch;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;
                if strcmp(Branch{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
                        Branch{ii}.Status = 'Stoped';
                        total_fixed = total_fixed + 1;         
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{1}{indtrack, icol} = strcat('Branch',num2str(ii), ', from Damaged to Stoped');
       
                        % the following lines about rescheduling is
                        % confusing??
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        index_task=index_task+1;
                        if index_task+1> size(Remain_schedule.Schedule_Power,2)
                            New_Schedule_Power=[];
                        else
                        [New_Schedule_Power{1:size(Remain_schedule.Schedule_Power,2)-index_task}] = Remain_schedule.Schedule_Power{index_task+1:end};
                        end
                        Remain_schedule.Schedule_Power=New_Schedule_Power;New_Schedule_Power=[];
   
                end
                
                %-- If a power line is fixed but in the "Stoped" status, 
                % as long as the two end objects are "Open", change the status
                % from "Stoped" to "Open", and set the Object.Functionality
                % = 1. 
                if strcmp(Branch{ii}.Status, 'Stoped')
                    obj1 = Dictionary(Branch{ii}.connectedObj1);
                    obj2 = Dictionary(Branch{ii}.connectedObj2);
                    if strcmp(obj1{1}.Status, 'Open') && strcmp(obj2{1}.Status, 'Open')
                        Branch{ii}.Status = 'Open';
                        Branch{ii}.Functionality = 1;
                        temp = funcTable(Branch{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(Branch{ii}.uniqueID) = temp;
                        powerGraph = Library.addPowerGraph(Power, powerGraph, ii, Dictionary);     
                        
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{1}{indtrack, icol} = strcat('Branch',num2str(ii), ', Status changes from Stoped to Open');
                        
                    end
                end
            end
            
            %==== Communication
            icol = 1;
            
            %--- Central Office
            ObjSet = []; ObjSet = Centraloffice;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;    
                %-- If a central office is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Stoped".
                if strcmp(Centraloffice{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
%                     if active_power % power system is turned on
                        Centraloffice{ii}.Status = 'Stoped';
%                     else % power system is turned off
%                         Centraloffice{ii}.Status = 'Open';
%                         Centraloffice{ii}.Functionality = 1;
%                         temp = funcTable(Centraloffice{ii}.uniqueID);
%                         temp(indtrack) = Centraloffice{ii}.Functionality;
%                         funcTable(Centraloffice{ii}.uniqueID) = temp;
%                         commGraph = addnode(commGraph,Centraloffice{ii}.uniqueID);
%                         fixed = true;
%                         total_fixed = total_fixed + 1;
%                         %--- track which Object changes status 
%                         icol = icol + 1; 
%                         trackOS{2}{indtrack, icol} = strcat('Centraloffice',num2str(ii), ', from Damaged to Open');
%                     end
                end
                
                % InterdependenceFunc == 1: inter-system functionality dependency is turned ON
                if strcmp(Centraloffice{ii}.Status, 'Stoped') && InterdependenceFunc == 1 
                    
                    bus = Dictionary(Centraloffice{ii}.Bus);
                    if strcmp(bus{1}.Status, 'Open')
                        Centraloffice{ii}.Status = 'Open';
                        Centraloffice{ii}.Functionality = 1;
                        temp = funcTable(Centraloffice{ii}.uniqueID);
                        temp(indtrack) = Centraloffice{ii}.Functionality;
                        funcTable(Centraloffice{ii}.uniqueID) = temp;
                        commGraph = addnode(commGraph,Centraloffice{ii}.uniqueID);
                        fixed = true;
                        total_fixed = total_fixed + 1;
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{2}{indtrack, icol} = strcat('Centraloffice',num2str(ii), ', Status changes from Stoped to Open');
                    
                    elseif Centraloffice{ii}.Battery == 1
                        Centraloffice{ii}.Status = 'Open';
                        Centraloffice{ii}.Functionality = 1;
                        temp = funcTable(Centraloffice{ii}.uniqueID);
                        temp(indtrack) = Centraloffice{ii}.Functionality;
                        funcTable(Centraloffice{ii}.uniqueID) = temp;
                        commGraph = addnode(commGraph,Centraloffice{ii}.uniqueID);
                        fixed = true;
                        total_fixed = total_fixed + 1;
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{2}{indtrack, icol} = strcat('Centraloffice',num2str(ii), ', from Stoped to Open');    
                    end
                end
                
                if strcmp(Centraloffice{ii}.Status, 'Stoped') && InterdependenceFunc == 0
                    Centraloffice{ii}.Status = 'Open';
                    Centraloffice{ii}.Functionality = 1;
                    temp = funcTable(Centraloffice{ii}.uniqueID);
                    temp(indtrack) = Centraloffice{ii}.Functionality;
                    funcTable(Centraloffice{ii}.uniqueID) = temp;
                    commGraph = addnode(commGraph,Centraloffice{ii}.uniqueID);
                    fixed = true;
                    total_fixed = total_fixed + 1;
                    %--- track which Object changes status 
                    icol = icol + 1; 
                    trackOS{2}{indtrack, icol} = strcat('Centraloffice',num2str(ii), ', from Damaged to Open');
                end    
                                
                %-- If a central office is fixed (fixed=true, "Open"), update the 
                % communication functionality property of the neighborhood (Neighborhood{j}.CommStatus) 
                % that connects to this restored/fixed central office.  
                if Centraloffice{ii}.Functionality
                    for j = 1:length(Centraloffice{ii}.Neighborhood)
                        temp = Dictionary(Centraloffice{ii}.Neighborhood{j});
                        temp = temp{1};
                        % disp(temp.uniqueID);
                        temp.CommStatus = 1;
                        temp2 = Dictionary(Centraloffice{ii}.Neighborhood_Comm_Link{j});
                        temp2 = temp2{1};
                        temp2.Status = 'Open';
                        
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{2}{indtrack, icol} = strcat(temp.uniqueID,': Nbr.CO, Status changes from Stoped to Open');
                    end
                end
            end
            
            %--- CommunicationTower
            ObjSet = []; ObjSet = CommunicationTower;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;  
     
                %-- If a communication tower is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Stoped". If 
                if strcmp(CommunicationTower{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
%                     if active_power 
                        CommunicationTower{ii}.Status = 'Stoped';
%                     end
                end

                if strcmp(CommunicationTower{ii}.Status, 'Stoped') && InterdependenceFunc == 1
                    if active_power
                    bus = Dictionary(CommunicationTower{ii}.Bus);
                    if strcmp(bus{1}.Status, 'Open')
                        CommunicationTower{ii}.Status = 'Open';
                        temp = funcTable(CommunicationTower{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(CommunicationTower{ii}.uniqueID) = temp;
                        CommunicationTower{ii}.Functionality = 1;
                        commGraph = addnode(commGraph,CommunicationTower{ii}.uniqueID);              
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{2}{indtrack, icol} = strcat('CT', num2str(ii), ', Status changes from Stoped to Open'); 
                    elseif CommunicationTower{ii}.Battery == 1
                        CommunicationTower{ii}.Status = 'Open';
                        temp = funcTable(CommunicationTower{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(CommunicationTower{ii}.uniqueID) = temp;
                        CommunicationTower{ii}.Functionality = 1;
                        commGraph = addnode(commGraph,CommunicationTower{ii}.uniqueID);
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{2}{indtrack, icol} = strcat('CT',num2str(ii), ', Status changes from Stoped to Open');
                    end
                    end

                end
                
                if strcmp(CommunicationTower{ii}.Status, 'Stoped') && InterdependenceFunc == 0
                    CommunicationTower{ii}.Status = 'Open';
                    CommunicationTower{ii}.Functionality = 1;
                    temp = funcTable(CommunicationTower{ii}.uniqueID);
                    temp(indtrack) = CommunicationTower{ii}.Functionality;
                    funcTable(CommunicationTower{ii}.uniqueID) = temp;
                    commGraph = addnode(commGraph,CommunicationTower{ii}.uniqueID);
                    fixed = true;
                    total_fixed = total_fixed + 1;
                    %--- track which Object changes status 
                    icol = icol + 1; 
                    trackOS{2}{indtrack, icol} = strcat('CT',num2str(ii), ', Status changes from Damaged to Open');
                end                    
                    
            end
            
            
            %--- Cell Line
            ObjSet = []; ObjSet = Cellline;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;  
                %-- If a cell line is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Stoped". 
                if strcmp(Cellline{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
%                     if active_power
                        Cellline{ii}.Status = 'Stoped';
%                     end     
                end
                if strcmp(Cellline{ii}.Status, 'Stoped') && InterdependenceFunc == 1 
                    bus = Dictionary(Cellline{ii}.Bus);
                    obj1 = Dictionary(Cellline{ii}.connectedObj1);
                    obj2 = Dictionary(Cellline{ii}.connectedObj2);
                    if strcmp(bus{1}.Status, 'Open') && strcmp(obj1{1}.Status, 'Open') && strcmp(obj2{1}.Status, 'Open')
                        Cellline{ii}.Status = 'Open';
                        temp = funcTable(Cellline{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(Cellline{ii}.uniqueID) = temp;
                        Cellline{ii}.Functionality = 1;
                        commGraph = Library.addCommGraph(Comm, commGraph, ii, Dictionary);
                        total_fixed = total_fixed + 1; 
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{2}{indtrack,icol} = strcat('CL',num2str(ii), ', Status changes from Stoped to Open'); 
                    end
                end
                if strcmp(Cellline{ii}.Status, 'Stoped') && InterdependenceFunc == 0 
                    if strcmp(obj1{1}.Status, 'Open') && strcmp(obj2{1}.Status, 'Open')
                        Cellline{ii}.Status = 'Open';
                        temp = funcTable(Cellline{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(Cellline{ii}.uniqueID) = temp;
                        Cellline{ii}.Functionality = 1;
                        commGraph = Library.addCommGraph(Comm, commGraph, ii,Dictionary);
                        total_fixed = total_fixed + 1; 
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{2}{indtrack,icol} = strcat('CL',num2str(ii), ', Status changes from Stoped to Open'); 
                    end
                end
            end

            %==== Transportation
            icol = 1;
            
            %--- TrafficLight
            ObjSet = []; ObjSet = TrafficLight;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;
                if strcmp(TrafficLight{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
                    if active_power
                        TrafficLight{ii}.Status = 'Stoped';
                    end
                end

                if strcmp(TrafficLight{ii}.Status, 'Stoped') && InterdependenceFunc == 1
                    if TrafficLight{ii}.Battery == 1
                        TrafficLight{ii}.Status = 'Open';
                        temp = funcTable(TrafficLight{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(TrafficLight{ii}.uniqueID) = temp;          
                        TrafficLight{ii}.Functionality = 1;
                        transGraph = addnode(transGraph,TrafficLight{ii}.uniqueID);         
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{3}{indtrack, icol} = strcat('TL',num2str(ii), ', Status changes from Stoped to Open');
                    end
                    temp = extractAfter(TrafficLight{ii}.Bus, 3);
                    temp = str2num(temp);
                    if ~isempty(temp) && strcmp(Bus{temp}.Status, 'Open')
                        TrafficLight{ii}.Status = 'Open';
                        temp = funcTable(TrafficLight{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(TrafficLight{ii}.uniqueID) = temp;
                        TrafficLight{ii}.Functionality = 1;
                        idword{1} = TrafficLight{ii}.uniqueID;
                        table = table2cell(transGraph.Nodes);
                        if ~ismember(table,idword)
                            transGraph = addnode(transGraph,TrafficLight{ii}.uniqueID); 
                        
                            %--- track which Object changes status 
                            icol = icol + 1; 
                            trackOS{3}{indtrack, icol} = strcat('TL',num2str(ii), ', Status changes from Stoped to Open');
                    
                        end
                    end
                end
                
                if strcmp(TrafficLight{ii}.Status, 'Stoped') && InterdependenceFunc == 0
                    TrafficLight{ii}.Status = 'Open';
                    temp = funcTable(TrafficLight{ii}.uniqueID);
                    temp(indtrack) = 1;
                    funcTable(TrafficLight{ii}.uniqueID) = temp;
                    TrafficLight{ii}.Functionality = 1;
                    transGraph = addnode(transGraph,TrafficLight{ii}.uniqueID); 
                    %--- track which Object changes status 
                    icol = icol + 1; 
                    trackOS{3}{indtrack, icol} = strcat('TL',num2str(ii), ', Status changes from Stoped to Open');
                end
                
            end
            
            %--- Bridge
            % Bridges without requiring the sub-component analysis
            ObjSet = []; ObjSet = Bridge;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;
                flag = 1;
                if strcmp(Bridge{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1)), tolerance)) && Bridge{ii}.HasSub == 0 
                    tasks = Bridge{ii}.taskUniqueIds;
                    sumWorkDay = zeros(1,6);
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        iduration = temp.WorkingDays;
                        if any(gt(abs(iduration(1:4)), tolerance))
                            sumWorkDay = sumWorkDay + iduration(1:4);
                            flag = 0;
                        end
                    end
                    Bridge{ii}.WorkingDays = sumWorkDay;
                    if Bridge{ii}.WorkingDays(1) <= 0 || flag
                        Bridge{ii}.Status = 'Open';
                        temp = funcTable(Bridge{ii}.uniqueID);
                        temp(indtrack) = 1;
                        funcTable(Bridge{ii}.uniqueID) = temp;
                        Bridge{ii}.Functionality = 1;
                        total_fixed = total_fixed + 1;                     
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{3}{indtrack, icol} = strcat('Bridge',num2str(ii), ', Status changes from Damaged to Open');
                    end
                end
                % for subcomponent analysis
                if strcmp(Bridge{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1)), tolerance)) && Bridge{ii}.HasSub == 1
                    for sub_index = 1:length(Bridge{ii}.ColumnSet)
                        tasks = Bridge{ii}.ColumnSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            iduration = temp.WorkingDays;
                            if any(gt(abs(iduration(1:4)), tolerance))
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.ColumnFoundSet)
                        tasks = Bridge{ii}.ColumnFoundSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            iduration = temp.WorkingDays;
                            if any(gt(abs(iduration(1)), tolerance))
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.AbutmentSet)
                        tasks = Bridge{ii}.AbutmentSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            iduration = temp.WorkingDays;
                            if any(gt(abs(iduration(1)), tolerance))
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.AbutmentFoundSet)
                        tasks = Bridge{ii}.AbutmentFoundSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            iduration = temp.WorkingDays;
                            if any(gt(abs(iduration(1)), tolerance))
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.BearingSet)
                        tasks = Bridge{ii}.BearingSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            iduration = temp.WorkingDays;
                            if any(gt(abs(iduration(1)), tolerance))
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.SlabSet)
                        tasks = Bridge{ii}.SlabSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            iduration = temp.WorkingDays;
                            if any(gt(abs(iduration(1)), tolerance))
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                   
                    if Bridge{ii}.WorkingDays(1) <= 0 || flag
                        Bridge{ii}.Status = 'Open';
                        temp = funcTable(Bridge{ii}.uniqueID);
                        temp(Start_Day:end) = 1;
                        funcTable(Bridge{ii}.uniqueID) = temp;
                        Bridge{ii}.Functionality = 1;
                        total_fixed = total_fixed + 1;
                        need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end                   
                end

            end
            
            % Road
            ObjSet = []; ObjSet = Road;
            for ii = 1:length(ObjSet)
                iObj = ObjSet{ii};
                iDuration = iObj.WorkingDays;
                if strcmp(Road{ii}.Status, 'Damaged') && any(lt(abs(iDuration(1:4)), tolerance))
                    Road{ii}.Status = 'Stoped';
                end
                if strcmp(Road{ii}.Status, 'Stoped')
                    flag = 0;
                    for j = 1:length(Road{ii}.Bridge_Carr)
                        if strcmp(Bridge{Road{ii}.Bridge_Carr(j)}.Status, 'Damaged')
                            flag = 1;
                            break;
                        end
                    end
                    for j = 1:length(Road{ii}.Bridge_Cross)
                        if strcmp(Bridge{Road{ii}.Bridge_Cross(j)}.Status, 'Damaged')
                            flag = 1;
                            break;
                        end
                    end
                    
                    if flag == 0
                        if active_power % if the power system is ON. 
                            tlID = Road{ii}.TrafficLight;
                            if ~isempty(tlID) % if there is a traffic light associated with this road link.
                                for j = 1:length(tlID)
                                    jj = tlID(j); 
                                    J = strcat('TrafficLight', num2str(jj));
                                    tl = Dictionary(J);
                                    tlStatus{j} = tl{1}.Status;
                                end
                                check = strcmp(tlStatus, 'Open');
                                 
                                if any(check == 0) % if any traffic light on this road link is not functional at the moment.
                                    temp = funcTable(Road{ii}.uniqueID);
                                    temp(indtrack) = 0.5;
                                    funcTable(Road{ii}.uniqueID) = temp;
                                    Road{ii}.Functionality = 0.5;
                                else
                                    Road{ii}.Status = 'Open'; % traffic light is functional.
                                    temp = funcTable(Road{ii}.uniqueID);
                                    temp(indtrack) = 1;
                                    funcTable(Road{ii}.uniqueID) = temp;
                                    Road{ii}.Functionality = 1;
                                    transGraph = Library.addTransGraph(Trans, transGraph, ii);
                                end
                            else
                                Road{ii}.Status = 'Open';
                                temp = funcTable(Road{ii}.uniqueID);
                                temp(indtrack) = 1;
                                funcTable(Road{ii}.uniqueID) = temp;
                                Road{ii}.Functionality = 1;
                                transGraph = Library.addTransGraph(Trans, transGraph, ii);
                            end
                        else % if the power system is OFF.
                            Road{ii}.Status = 'Open';
                            temp = funcTable(Road{ii}.uniqueID);
                            temp(indtrack) = 1;
                            funcTable(Road{ii}.uniqueID) = temp;
                            Road{ii}.Functionality = 1;
                            transGraph = Library.addTransGraph(Trans, transGraph, ii);
                        
                            %--- track which Object changes status 
                            icol = icol + 1; 
                            trackOS{3}{indtrack, icol} = strcat('Road',num2str(ii), ', from Stoped to Open');
                        end
                        
                        temp = Dictionary(strcat('RoadNode', num2str(Road{ii}.Start_Node)));
                        temp = temp{1};
                        for j = 1:length(temp.Neighborhood)
                            t1 = Dictionary(temp.Neighborhood{j});
                            t1 = t1{1};
                            t1.TransStatus = 1;
                            t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                            t1 = t1{1};
                            t1.Status = 'Open';
                            
                            %--- track which Object changes status 
                            icol = icol + 1; 
                            trackOS{3}{indtrack, icol} = strcat(t1.uniqueID, ', Status changes from Stoped to Open');
                        end
                    end

                    temp = Dictionary(strcat('RoadNode', num2str(Road{ii}.End_Node)));
                    temp = temp{1};
                    for j = 1:length(temp.Neighborhood)
                        t1 = Dictionary(temp.Neighborhood{j});
                        t1 = t1{1};
                        t1.TransStatus = 1;
                        t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                        t1 = t1{1};
                        t1.Status = 'Open';
                            
                        %--- track which Object changes status 
                        icol = icol + 1; 
                        trackOS{3}{indtrack, icol} = strcat(t1.uniqueID, ': Nbr.Road, Status changes from Stoped to Open');
                    end
                end
            end
            
            %==== Update Object.Status to every SystemSet
            Power{1} = Branch;
            Power{2} = Bus;
            Power{3} = Generator;
            Power{4} = TransTower;
            
            Comm{1} = Centraloffice;
            Comm{2} = CommunicationTower;
            Comm{3} = Cellline;

            Trans{1} = Road;
            Trans{2} = Bridge;
            Trans{3} = TrafficLight;   
            
            %==== Update system functionality
            [FunctionalityPower, FunctionalityComm, FunctionalityTrans] = Library.ComputeSystemFunctionality(ActiveSystem, LinkDirection, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, Power, Comm, Trans, powerGraph, commGraph, transGraph, Dictionary);
            [totalPopulation, FuncPowerNbr, FuncCommNbr, FuncTransNbr] = Library.neighbourFunc(Dictionary); % Percentage of population that has power/comm/trans service at the Neighborhood level
            %FunctionalityTotal = [FunctionalityPower, FunctionalityComm, FunctionalityTrans, FuncPowerNbr, FuncCommNbr, FuncTransNbr];
            FunctionalityTotal{1} = FunctionalityPower;
            FunctionalityTotal{2} = FunctionalityComm;
            FunctionalityTotal{3} = FunctionalityTrans;
            FunctionalityTotal{4} = FuncPowerNbr;
            FunctionalityTotal{5} = FuncCommNbr;
            FunctionalityTotal{6} = FuncTransNbr;
            
            %==== Calculate System Functionality after Update Object.Status
            % Do we need reschedule? If yes, need_reschedule = 1.             
%             need_reschedulePower = Library.checkNeedRescheduleNew(need_reschedulePower, FunctionalityPower, Start_Day, End_Day, Per_New_reschedule,Num_stage);
%             need_rescheduleComm = Library.checkNeedRescheduleNew(need_rescheduleComm, FunctionalityPower, Start_Day, End_Day, Per_New_reschedule,Num_stage);            
%             need_rescheduleTrans = Library.checkNeedRescheduleNew(need_rescheduleTrans, FunctionalityPower, Start_Day, End_Day, Per_New_reschedule,Num_stage);
%             need_reschedule = [need_reschedulePower,need_rescheduleComm, need_rescheduleTrans];
%             
        end        
        
        %% Add tasks for restoring damaged components to the current working list for every system (STEP 4)   
        % as a way to simulate the execution process of retoration tasks 
        function [Current, Schedule, Max_Resource, lookupTable] = AddCurrentWorking(Max_Resource, Current, Schedule, Dictionary, lookupTable, CurrentTime, n, indtrack)
            %=================================================================
            % function AddCurrentWorking
            % This function evaluates the available resource constraint at the moment
            % and considers precedence constraints of different tasks to
            % add a task/tasks to the list of CurrentWorking. 
            %=================================================================
            Ntask = length(Schedule);  % Number of tasks in the Schedule
            Nr = length(Max_Resource); % Number of resource type
            
            %=== add every task to the current working list named "Current"
            for ii = 1:Ntask

                %--- A flag variable representing whether a variable is a task or
                %not (0-no/1-yes)
                isTask = 0;
                
                %--- Check if the task in the 'tem' after splitting by '/' 
                % satisfies the formats with the length 2 and 3! 
                %- length(tem)==2, meaning tem is some task like
                % 'Bridge/243' (i.e. there is only 1 '/' in the task string).
                %- length(tem)==3, meaning tem is some task like
                % 'Bridge/243/Task1134' (i.e. there are 2 '/' in the task string).
                %- The following "contunue" can Skip tem cells to Next Loop
                % Iteration unless the tem cell contains a string witht he
                % length of 2 or 3 (some task like 'Bridge/243' or 'Bridge/243/Task1134'). 
                %- A side note: 
                % If a task is been executing and not finished yet, this
                % task string of 'tem' would be some string like
                % 'Bridge/243/Task1134/Working/dummy'; therefore,
                % length(tem)=5. This task will be jumped over because of
                % using 'continue', as this task is being executing right
                % now or has been executed. 
                tem = strsplit(Schedule{ii},'/');
                if length(tem) ~= 2 && length(tem) ~= 3
                    continue;
                end
                
                %--- flagResource: A flag variable representing resource
                %availability [before] adding a task in the list of Current
                % (0-no/1-yes)
                flagResource = 0;
                  
                %--- allEmpty: A flag variable representing whether there is 
                % any resource aviable to execute a task immediately [after] 
                % adding a task in the list of Current.  
                % allEmpty=1, there is no resource available at the moment.
                % allEmpty=0, there is resource available at the moment.
                allEmpty = 1;
                
                %--- Set up the initial vector for the resource demand of every task
                resourceNeed = zeros(1,Nr);
 
%                 %=== Revised form about getDependency if there is a 
%                 %problem that the Dictionary cannot find the Key
%                 if Library.getDependency(Schedule{i},Dictionary) == 1
%                     objTemp = Dictionary(Library.getUniqueId(Schedule{i},0));
%                 end
%                 objTemp = Dictionary(Library.getUniqueId(Schedule{ii},0));
%                 %=== Or Revised form
%                 % Library.getDependency(Schedule{ii},Dictionary)==1 meaning either of the following two cases:
%                 % 1. find dependent task(s) for the given task of Schedule{ii}
%                 % 2. there is no dependent task(s) for the given task of Schedule{ii}
%                 % Given the above condition, assign objTemp = the Task cell
%                 % defined in the DictionaryDamage. It should be in the form
%                 % of the following example task.             
%                 if(Library.getDependency(Schedule{ii},Dictionary)==1)
%                      continue;
%                 end
%                 
% %                 Library.getDependency(Schedule{ii},Dictionary)
%                 objTemp = Dictionary(Library.getUniqueId(Schedule{ii},0));
%         
 
                %=== Original form about getDependency if there is no 
                %problem that the Dictionary cannot find the Key
                if(Library.getDependency(Schedule{ii}, Dictionary))
                    continue;
                end
                objTemp = Dictionary(Library.getUniqueId(Schedule{ii},0));
                
                %=== Check if there're enough resources to complete the task
                for j = 1:Nr % the index representing the resource type 
                    % Extracting the resource demand for the j reosurce, 
                    % given a task of Schedule{ii}.
                    resourceNeed(j) = Library.getResource(Schedule{ii},Dictionary,j, 0);
                    %--- flagResource = 1: the avaiable amount of resource at
                    % the moment is less then the resource demand.
                    if (Max_Resource(j) - resourceNeed(j)) < 0 
                        flagResource = 1;
                        break;
                    end
                end
                if flagResource == 1
                    continue;
                end
                %=== Check if the object is a task in the form of
                %tem='Bus/9/Task12' (meaning length(tem) == 3 while spliting by '/'). 
                if length(tem) == 3
                    isTask = 1;
                    tasktemp = Dictionary(tem{3}); % tasktemp: a string of the unique task number, in the form of 'Task1158'.
                    tasktemp.WorkingDays(5) = tasktemp.WorkingDays(2);
                    tasktemp.WorkingDays(6) = n;

 %%%%%%%%%%%% Revised codes as follows about identifying the correct duration value in the right time unit  
                        % for all three systems
                        % time unit = hour in the first 3 days
                        % time unit = day in the days of 3~28
                        % time unit = week in the days of 28~168
                        % n = 1, CurrentTime in the range of [0 3) days
                        % n = 2, CurrentTime in the range of [3 28) days
                        % n = 3, CurrentTime in the range of [28 168) days
                        % n = 4, CurrentTime in the range of [168 time_horizon) days
                        
%                     % Diff_unit==1 meaning the choice of using different time unit in recovery is turned on.        
%                     if Diff_unit == 1 
%                         [~,n]=find(sort([CurrentTime, RestorationPhase]) == CurrentTime); % n = either one of (1,2 3, 4)
% 
%                         if gt(n,length(RestorationPhase)+1) % n>4 if length(RestorationPhase) = 3 [e.g., RestorationPhase = [3,28,168]];
%                             msg = strcat('Function AddCurrentWorking Error: Time Index of CurrentTime (n) = ', num2str(n), ' > length(RestorationPhase)+1 for ', tasktemp);
%                             disp(msg); 
%                             return
%                         else % n<=length(RestorationPhase)+1
%                             tasktemp.WorkingDays(5) = tasktemp.WorkingDays(n);% the time duration based on the selected time unit 
%                             tasktemp.WorkingDays(6) = n; % the index of time unit in the unit vector of [hour,day,week,month]                  
%                         end
%                         
%                     % Diff_unit is not 1, meaning the choice of using different time unit in recovery is turned off.  
%                     % All tasks using the same time unit of "day"
%                     else
%                         n = 2;
%                         tasktemp.WorkingDays(5) = tasktemp.WorkingDays(n);
%                         tasktemp.WorkingDays(6) = n;
% 
%                     end
                        
%%%%%%%%%%%% Original codes as follows about identifying the correct duration value in the right time unit 
%                         %--- find which time range does t=End_Day falls into
%                         % in the time vector of [0 3 28 168 timehorizon](task unit: day for trans and hour for power and communication system)
%                         % considering 1 day = 24 hour, 1 week = 7 days, 1 month = 28 days
%                         % for power and communciation systems
%                         % time unit = hour in the first 3 days
%                         % time unit = day in the days of 3~28
%                         % time unit = week in the days of 28~168
%                         % time unit = month in the days of 168~time horizon
%                         % for transportation system
%                         % time unit = day in the first 3 days
%                         % time unit = day in the days of 3~28
%                         % time unit = week in the days of 28~168
%                         % time unit = month in the days of 168~time horizon
%                         % n = 1, CurrentTime in the range of [0 3) days
%                         % n = 2, CurrentTime in the range of [3 28) days
%                         % n = 3, CurrentTime in the range of [28 168) days
%                         % n = 4, CurrentTime in the range of [168 time_horizon) days
% %                      if Diff_unit==1 % Diff_unit==1 meaning the choice of using different time unit in recovery is turned on.
% %                         if strcmp(system, 'Trans')
% %                            [~,n]=find(sort([CurrentTime,Cust_unit])==End_Day);
% %                         else
% %                            [~,n]=find(sort([CurrentTimeCust_unit*24])==End_Day);
% %                         end
% %                         if strcmp(system, 'Trans')    %, system = 'Trans', transportation system                   
% % %                             if n==2 || n==3
% % %                                tasktemp.WorkingDays(1)=tasktemp.WorkingDays(2);
% %                             if le(n,2) % n<=2, meaning t<=28 days, TimeUnit = day for transportation 
% %                                tasktemp.WorkingDays(5)=tasktemp.WorkingDays(1);
% %                             elseif n==3 % n==3, meaning 28<t<168 days, TimeUnit = week for transportation
% %                                tasktemp.WorkingDays(5)=tasktemp.WorkingDays(2);
% %                             elseif n==4 % n==4, meaning t>168 days, TimeUnit = month for transportation
% %                                tasktemp.WorkingDays(5)=tasktemp.WorkingDays(3);
% %                             end
% %                         else %=== power and communication systems, system = 'Power' or 'Comm'
% %                            % The original time unit in the task library is hour
% %                            % t = [0~3~28~168~time_horizon]
% %                            %[~,n]=find(sort([CurrentTime,Cust_unit])==End_Day)
% %                            % (1)CurrentTime = 2(day), leads to n = 1, time unit is hour
% %                            % (2)CurrentTime = 4(day), leads to n = 2, time unit is day
% %                            % (3)CurrentTime = 30(day), leads to n = 3, time unit is week
% %                            % (4)CurrentTime = 180(day), leads to n = 4, time unit is month
% %                            % WorkingDays = task duration in [hour, day, week, month]
% %                            % here, day = 24 hours, week = 7 days, month = 28 days
% %                            if gt(n,1)
% %                                tasktemp.WorkingDays(5)=tasktemp.WorkingDays(n);
% %                            end
% %                         end
% %                     end
                    
                    %--- A flag variable to represent whether a task has 
                    % any predecessor task or not (0-no/1-yes)
                    flagPredecessor = 0;
                    % Check if dependent tasks are already finished
                    for j = 1:length(tasktemp.predecessorTask)
                        %dpre: the duration vector of the precedessor
                        %representing the same duration in different time
                        %unit.
                        %---- A general explnation about the duration vector in different time unit 
                        %for Task.WorkingDays is as follows.
                        %-[hour,day,week,month, the selected time unit based on the Start_Day falling which time range] for power and communication systems 
                        %-[day,week,month, 0, the selected time unit based on the Start_Day falling which time range] for trnasportation system. 
                        %-The last element in the duration vector (icol = 5)
                        % means the duration in a selected time unit given
                        % the current time step, and this duration value
                        % will be used adn updated in different unit as
                        % the current time step jumps from one recovery phase 
                        % to another in the recovery process at different steps. 
                        %---- reading more details in the above and in the function of SampleActualDuration
                        %
                        dpre = Dictionary(tasktemp.predecessorTask{j}).WorkingDays; 
                        if any(dpre(1:4) > 0) % meaning if any element in the duration vector, i.e. the duration in any time unit, is still great than 0.
                            flagPredecessor = 1;
                            break;
                        end
                    end
                    if flagPredecessor == 1
                        continue;
                    end
                    if ~any(tasktemp.WorkingDays(1:4))
                        disp('---- Function AddCurrentWorking Error: Working Days = 0');
                        disp(tasktemp);
                        continue;
                    end
                end
                
                %=== Find an empty block to fill any task into the list of "Current". 
                %If there is any empty cell (even more than 1 empty cells) in Current,
                %indexCurrent = the index of this empty cell. 
                %If there is no empty cell in the Current list, assign
                %indexCurrent = length(Current) + 1, since there are available resources to execute this task 
                %when [(Max_Resource(j)- resourceNeed(j))>=0] from previous "for loop" for
                %checking the resource availability at the moment. 
                indexCurrent = 0;
                for j = 1:size(Current,2)
                    if isempty(Current{1,j})
                        indexCurrent = j;
                    end
                end

                if indexCurrent == 0
                    indexCurrent = size(Current,2) + 1;
                end
                
                %=== Mark the string of task/object as working by appending
                %'/Working/dummy' at the end of this string.
                % Add the task to the working list of "Current" and
                % Subtract the resource demand of this task.
                Schedule{ii} = strcat(Schedule{ii},'/Working/dummy');
                Current{1,indexCurrent} = Schedule{ii}; % Schedule{ii}: task ID that is in the CurrentWorking list at the current step, and this task ID ends with '/Working/dummy'
                Current{2,indexCurrent} = indtrack;     % track in which round (indtrack) this task of Current{1,indexCurrent} as been added to the CurrentWorking list.
                
                %=== Create a task table named "lookupTable" for all tasks
                %that have been executing (will be executed). 
                lookupTable = Library.addToLookupTable(lookupTable, objTemp, CurrentTime, isTask);
                
                %=== After adding tasks to the Current list, check the resource avaiability at the momenet again.
                % If there is any available resource, allEmpty = 0.
                % If there is no available resource, allEmpty = 1. Exit this function.
          
                for j = 1:Nr
                    Max_Resource(j) = Max_Resource(j) - resourceNeed(j);
                    if Max_Resource(j) ~= 0
                        allEmpty = 0;
                    end
                end
                                
                if allEmpty == 1
                    return;
                end
            end % for loop - Schedule{ii}
        end % function AddCurrentWorking

        %% function UpdateStatusSep (STEP 4) 
        function [Remain_schedule, need_reschedule, total_fixed, transGraph, powerGraph, commGraph, Pow, Comm, Trans, funcTable, FunctionalityTotal, totalPopulation] = UpdateStatusSep(Pow, Comm, Trans, total_damaged, total_fixed, need_reschedule,Dictionary, transGraph, powerGraph, commGraph,funcTable,Start_Day,End_Day,Per_New_reschedule,Num_stage, Remain_schedule,ori_schedule, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirectionChoice)
            %==== Field
            Branch= Pow{1};
            Bus= Pow{2};
            Generator= Pow{3};
            TransTower = Pow{4};
            
            Centraloffice = Comm{1};
            CommunicationTower = Comm{2};
            Cellline = Comm{3};            
            
            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};
            fixed = false;
            index_task=0;
            
            need_reschedulePower = need_reschedule(1);
            need_rescheduleComm = need_reschedule(2);
            need_rescheduleTrans = need_reschedule(3);
            
            %==== Power System
            %---- Generator (Power Plant)
            for i = 1:length(Generator)
                if strcmp(Generator{i}.Status, 'Damaged') && ~isempty(Generator{i}.WorkingDays)
                    if Generator{i}.WorkingDays <= 0
                        Generator{i}.Status = 'Open';
                        
                        temp = funcTable(Generator{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Generator{i}.uniqueID) = temp;
                      
                        Generator{i}.Functionality = 1;
                        total_fixed = total_fixed + 1;
                        index_task = index_task+1;
 
                        [New_Schedule_Power{1:size(ori_schedule.Schedule_Power,2)-i}] = ori_schedule.Schedule_Power{i+1:end};
                        Remain_schedule.Schedule_Power = New_Schedule_Power;
                        New_Schedule_Power = [];
                    end
                end
            end
            
            %---- Bus (Substation)
            for i = 1:length(Bus)
                if ~isempty(Bus{i}.Generator)
                    if strcmp(Bus{i}.Status, 'Damaged') && ~isempty(Bus{i}.WorkingDays)
                        Bus = Library.getWorkDays(Bus, i,Dictionary);
                        if  Bus{i}.WorkingDays <= 0
                            Bus{i}.Status = 'Stoped';
                            total_fixed = total_fixed + 1;
%                             need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        end
                    end
                    if strcmp(Bus{i}.Status, 'Stoped')
                        temp = extractAfter(Bus{i}.Generator, 9);
                        temp = str2num(temp);
                        if strcmp(Generator{temp}.Status, 'Open')
                            Bus{i}.Status = 'Open';
                            fixed = true;
                            temp = funcTable(Bus{i}.uniqueID);
                            temp(End_Day:end) = 1;
                            funcTable(Bus{i}.uniqueID) = temp;
                            
                            Bus{i}.Functionality = 1;
                            powerGraph = addnode(powerGraph,Bus{i}.uniqueID);
                            
                        end
                    end
                else
                    if strcmp(Bus{i}.Status, 'Damaged') && ~isempty(Bus{i}.WorkingDays)
                        if  Bus{i}.WorkingDays <= 0
                            Bus{i}.Status = 'Open';
                            fixed = true;
                            Bus{i}.Functionality = 1;
                            
                            temp = funcTable(Bus{i}.uniqueID);
                            temp(End_Day:end) = 1;
                            funcTable(Bus{i}.uniqueID) = temp;
                            
                            powerGraph = addnode(powerGraph,Bus{i}.uniqueID);
                            total_fixed = total_fixed + 1;
%                             need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        end
                    end
                end
                if fixed
                    for j = 1:length(Bus{i}.Neighborhood)
                        temp = Dictionary(Bus{i}.Neighborhood{j});
                        temp = temp{1};
                        disp( temp.uniqueID);
                        temp.PowerStatus = 1;
                        temp = Dictionary(Bus{i}.Neighborhood_Power_Link{j});
                        temp = temp{1};
                        temp.Status = 'Open';
                    end
                    
                    fixed = false;
                end
            end
            
            %---- Branch (Power Line)
            for i = 1:length(Branch)
                if strcmp(Branch{i}.Status, 'Damaged') && ~isempty(Branch{i}.WorkingDays)      
                    Branch = Library.getWorkDays(Branch, i,Dictionary);          
                    if Branch{i}.WorkingDays <= 0
                        Branch{i}.Status = 'Stoped';
                        total_fixed = total_fixed + 1;
                        index_task=index_task+1;
                        if index_task+1> size(Remain_schedule.Schedule_Power,2)
                            New_Schedule_Power=[];
                        else
                        [New_Schedule_Power{1:size(Remain_schedule.Schedule_Power,2)-index_task}] = Remain_schedule.Schedule_Power{index_task+1:end};
                        end
                        Remain_schedule.Schedule_Power=New_Schedule_Power;New_Schedule_Power=[];
                    end
                end
                if strcmp(Branch{i}.Status, 'Stoped')
                    obj1 = Dictionary(Branch{i}.connectedObj1);
                    obj2 = Dictionary(Branch{i}.connectedObj2);
                    if strcmp(obj1{1}.Status, 'Open')&&strcmp(obj2{1}.Status, 'Open')
                        Branch{i}.Status = 'Open';
                        Branch{i}.Functionality = 1;
                        temp = funcTable(Branch{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Branch{i}.uniqueID) = temp;
                        powerGraph = Library.addPowerGraph(Pow, powerGraph, i, Dictionary);          
                    end
                end
            end
            
            %---- Transmission Tower
            for i = 1:length(TransTower)
                if strcmp(TransTower{i}.Status, 'Damaged') && ~isempty(TransTower{i}.WorkingDays)
                    if  TransTower{i}.WorkingDays <= 0
                        TransTower{i}.Status = 'Open';
                        TransTower{i}.Functionality = 1;
                        temp = funcTable(TransTower{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(TransTower{i}.uniqueID) = temp;
                        
                        powerGraph = addnode(powerGraph,TransTower{i}.uniqueID);
                        total_fixed = total_fixed + 1;   
                    end
                end
            end
            
            %==== Communication
            %---- Central Office
            for i = 1:length(Centraloffice)
                if strcmp(Centraloffice{i}.Status, 'Damaged') && ~isempty(Centraloffice{i}.WorkingDays)
                    Centraloffice = Library.getWorkDays(Centraloffice, i,Dictionary);
                    
                    if Centraloffice{i}.WorkingDays <= 0
                        Centraloffice{i}.Status = 'Open';
                        total_fixed = total_fixed + 1;
                        fixed = true;
                        temp = funcTable(Centraloffice{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Centraloffice{i}.uniqueID) = temp;
                        Centraloffice{i}.Functionality = 1;
                        commGraph = addnode(commGraph,Centraloffice{i}.uniqueID);
                    end
                end
                
                if fixed
                    for j = 1:length(Centraloffice{i}.Neighborhood)
                        temp = Dictionary(Centraloffice{i}.Neighborhood{j});
                        temp = temp{1};
                        disp(temp.uniqueID);
                        temp.CommStatus = 1;
                        temp = Dictionary(Centraloffice{i}.Neighborhood_Comm_Link{j});
                        temp = temp{1};
                        temp.Status = 'Open';
                    end
                    
                    fixed = false;
                end
            end
            
            %---- CommunicationTower
            for i = 1:length(CommunicationTower)
                if strcmp(CommunicationTower{i}.Status, 'Damaged') && ~isempty(CommunicationTower{i}.WorkingDays)
                    CommunicationTower = Library.getWorkDays(CommunicationTower, i,Dictionary);
                    if CommunicationTower{i}.WorkingDays <= 0
                        CommunicationTower{i}.Status = 'Open';
                        total_fixed = total_fixed + 1;
                        
                        temp = funcTable(CommunicationTower{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(CommunicationTower{i}.uniqueID) = temp;
                        CommunicationTower{i}.Functionality = 1;
                        commGraph = addnode(commGraph,CommunicationTower{i}.uniqueID);
                        
                        if i+1> size(ori_schedule.Schedule_Comm,2)
                            New_Schedule_Comm=[];
                        else
                        [New_Schedule_Comm{1:size(ori_schedule.Schedule_Comm,2)-i}] = ori_schedule.Schedule_Comm{i+1:end};
                        end
                        Remain_schedule.Schedule_Comm=New_Schedule_Comm;New_Schedule_Comm=[];
                    end
                end
                
            end
            
            %---- Cell Line
            for i = 1:length(Cellline)
                if strcmp(Cellline{i}.Status, 'Damaged') && ~isempty(Cellline{i}.WorkingDays)
                    Cellline = Library.getWorkDays(Cellline, i,Dictionary);
                    if Cellline{i}.WorkingDays <= 0
                        Cellline{i}.Status = 'Open';
                        total_fixed = total_fixed + 1;
                        temp = funcTable(Cellline{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Cellline{i}.uniqueID) = temp;
                        Cellline{i}.Functionality = 1;
                        commGraph = Library.addCommGraph(Comm, commGraph, i,Dictionary);
                        
                    end
                end

            end

            %==== Transportation
            %---- TrafficLight
            for i = 1:length(TrafficLight)
                if strcmp(TrafficLight{i}.Status, 'Damaged') && ~isempty(TrafficLight{i}.WorkingDays)
                    
                    TrafficLight = Library.getWorkDays(TrafficLight, i,Dictionary);
                    if TrafficLight{i}.WorkingDays <= 0
                        TrafficLight{i}.Status = 'Open';
                        total_fixed = total_fixed + 1;
                        temp = funcTable(TrafficLight{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(TrafficLight{i}.uniqueID) = temp;
                        TrafficLight{i}.Functionality = 1;
                        transGraph = addnode(transGraph,TrafficLight{i}.uniqueID);
                    end
                end
                
            end
            
            %---- Bridge
            for i = 1:length(Bridge)
                flag = 1;
                if strcmp(Bridge{i}.Status, 'Damaged') && ~isempty(Bridge{i}.WorkingDays) && Bridge{i}.HasSub == 0 
                    tasks = Bridge{i}.taskUniqueIds;
                    sumWorkDay = 0;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        if temp.WorkingDays > 0
                            sumWorkDay = sumWorkDay + temp.WorkingDays;
                            flag = 0;
                        end
                    end
                    Bridge{i}.WorkingDays = sumWorkDay;
                    if Bridge{i}.WorkingDays(1) <= 0 || flag
                        Bridge{i}.Status = 'Open';
%                         temp = funcTable(Bridge{i}.uniqueID);
%                         temp(Start_Day:end) = 1;
%                         funcTable(Bridge{i}.uniqueID) = temp;
                        Bridge{i}.Functionality = 1;
                        total_fixed = total_fixed + 1;
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end
                end
                % for subcomponents
                if strcmp(Bridge{i}.Status, 'Damaged') && ~isempty(Bridge{i}.WorkingDays) && Bridge{i}.HasSub == 1
                    for sub_index = 1:length(Bridge{i}.ColumnSet)
                        tasks = Bridge{i}.ColumnSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{i}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{i}.ColumnFoundSet)
                        tasks = Bridge{i}.ColumnFoundSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{i}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{i}.AbutmentSet)
                        tasks = Bridge{i}.AbutmentSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{i}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{i}.AbutmentFoundSet)
                        tasks = Bridge{i}.AbutmentFoundSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{i}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{i}.BearingSet)
                        tasks = Bridge{i}.BearingSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{i}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{i}.SlabSet)
                        tasks = Bridge{i}.SlabSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{i}.WorkingDays = sumWorkDay;
                    end
                   
                    if Bridge{i}.WorkingDays(1) <= 0 || flag
                        Bridge{i}.Status = 'Open';
                        Bridge{i}.Functionality = 1;
                        total_fixed = total_fixed + 1;
                        need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end                   
                end

            end
            
            %---- Road
            for i = 1:length(Road)
                if strcmp(Road{i}.Status, 'Damaged') &&  ~isempty(Road{i}.WorkingDays)
                    Road = Library.getWorkDays(Road, i,Dictionary);
                    if Road{i}.WorkingDays <= 0
                        total_fixed = total_fixed + 1;
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        Road{i}.Status = 'Stoped';
                    end
                end
                if strcmp(Road{i}.Status, 'Stoped')
                    flag = 0;
                    for j = 1:length(Road{i}.Bridge_Carr)
                        if strcmp(Bridge{Road{i}.Bridge_Carr(j)}.Status, 'Damaged')
                            flag = 1;
                            break;
                        end
                    end
                    for j = 1:length(Road{i}.Bridge_Cross)
                        if strcmp(Bridge{Road{i}.Bridge_Cross(j)}.Status, 'Damaged')
                            flag = 1;
                            break;
                        end
                    end
                    
                    if flag == 0
                        Road{i}.Status = 'Open';
                        temp = funcTable(Road{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Road{i}.uniqueID) = temp;
                        Road{i}.Functionality = 1;
                        transGraph = Library.addTransGraph(Trans, transGraph, i);
                        temp = Dictionary(strcat('RoadNode',num2str(Road{i}.Start_Node)));
                        temp = temp{1};
                        for j = 1:length(temp.Neighborhood)
                            t1 = Dictionary(temp.Neighborhood{j});
                            t1 = t1{1};
                            t1.TransStatus = 1;
                            t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                            t1 = t1{1};
                            t1.Status = 'Open';
                        end
                        temp = Dictionary(strcat('RoadNode',num2str(Road{i}.End_Node)));
                        temp = temp{1};
                        for j = 1:length(temp.Neighborhood)
                            t1 = Dictionary(temp.Neighborhood{j});
                            t1 = t1{1};
                            t1.TransStatus = 1;
                            t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                            t1 = t1{1};
                            t1.Status = 'Open';
                        end
                    end
                end
            end
            
            %==== Update Object.Status to SystemSet 
            Pow{1} = Branch;
            Pow{2} = Bus;
            Pow{3} = Generator;
            Pow{4} = TransTower;
            
            Comm{1} = Centraloffice;
            Comm{2} = CommunicationTower;
            Comm{3} = Cellline;            
            
            Trans{1} = Road;
            Trans{2} = Bridge;
            Trans{3} = TrafficLight;
 
            %==== Calculated System Functionality
            % Do we need reschedule? If yes, need_reschedule=1. 
            [FunctionalityPower, FunctionalityComm, FunctionalityTrans] = Interface1.Functionality(Power_Func_Num, Trans_Func_Num, Comm_Func_Num, Pow, Comm, Trans, powerGraph, commGraph, transGraph, Dictionary,LinkDirectionChoice);
            [totalPopulation,FuncPowerNbr,FuncCommNbr,FuncTransNbr] = Library.neighbourFunc(Dictionary);% Propotion of population that has power/comm/trans
            FunctionalityTotal = [FunctionalityPower, FunctionalityComm, FunctionalityTrans, FuncPowerNbr,FuncCommNbr,FuncTransNbr];
         
        end

        %% function checkNeedReschedule: Calculated the need for reschedule (STEP 4)  
        function need_reschedule = checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule)
            % =================================================================
            % function checkNeedReschedule
            % Calculate the need for reschedule
            % =================================================================       
            if need_reschedule == 0
                if isempty(Per_New_reschedule)      
                else
%                     if ~isempty(Per_New_reschedule)
                     if (total_fixed/total_damaged) >= Per_New_reschedule(1)
                           need_reschedule = 1;
%                            Per_New_reschedule(1)=[];
                     end
%                     end
                end
            end
        end
        
        %% function checkNeedRescheduleNew (STEP 4) 
        function [need_reschedule,Per_New_reschedule,Num_stage] = checkNeedRescheduleNew(need_reschedule, FunctionalitySystem, Start_Day, End_Day, Per_New_reschedule,Num_stage)
            if need_reschedule == 0 
                if isempty(Per_New_reschedule) && isempty(Per_New_reschedule)      
                    need_reschedule = -1; 
                else 

                    % this may need to be fixed as Functionality(this time step)<Per_New_reschedule(1)<=Functionality(next time step)
                   if Per_New_reschedule~=0
                       if End_Day> 10
                           End_Day=10;
                       end
                       if  ge(FunctionalitySystem{ceil(End_Day)}, Per_New_reschedule(1))
                           need_reschedule = 1;
                           Per_New_reschedule(1) = [];
                       end
                   end
                    
                    % Start_Day (this time step) < Num_stage(1) <= End_Day (next time step)
                   if Num_stage~=0
                        if  ge(End_Day, Num_stage(1)) 
                            need_reschedule = 1;
                            Num_stage(1) = [];
                        end
                   end
               end
            % if need_reschedule ~=0, do not reschedule, exit this function and continue other computations    
            end
        end
  
        %% create dictionary of task and precedence 
        function [Dic_p_task,Dic_c_task,Dic_t_task,Dic_p_prece,Dic_c_prece,Dic_t_prece] = createDictionaryTaskAndPrecedence(taskTable,precedenceTable)
            
            
            
            power_task=taskTable{1};
            comm_task=taskTable{2};
            trans_task=taskTable{3};
            
            power_prece=precedenceTable{1};
            comm_prece=precedenceTable{2};
            trans_prece=precedenceTable{3};
            
            %==== Set up the initial dictionary of tasks and precedence for three systems
            Dic_p_task=[];
            Dic_c_task=[];
            Dic_t_task=[];
            
            Dic_p_prece=[];
            Dic_c_prece=[];
            Dic_t_prece=[];
            
            % Power
            if ~isempty(power_task)
                [key_p{1:size(power_task,1),1}]=power_task{:,1};
                for i=1:size(power_task,1)
                    val_p{i,:}={power_task{i,2:size(power_task,2)}};
                end
                Dic_p_task = containers.Map(key_p,val_p);
            end

            % Communication
            if ~isempty(comm_task)
                [key_c{1:size(comm_task,1),1}]=comm_task{:,1};
                for i=1:size(comm_task,1)
                    val_c{i,:}={comm_task{i,2:size(comm_task,2)}};
                end
                Dic_c_task = containers.Map(key_c,val_c);
            end

            %Transportation
            if ~isempty(trans_task)
                [key_t{1:size(trans_task,1),1}]=trans_task{:,1};
                for i=1:size(trans_task,1)
                    val_t{i,:}={trans_task{i,2:size(trans_task,2)}};
                end
                Dic_t_task = containers.Map(key_t,val_t);
            end
            clear key_p val_p key_c val_c key_t val_t
            % dictionary of precedence of three sys
            
            % Power
            if ~isempty(power_prece)
                [key_p{1:size(power_prece,1)-1,1}]=power_prece{2:end,1};
                for i=1:size(power_prece,1)-1
                    val_p{i,:}={power_prece{i+1,2:size(power_prece,2)}};
                end
                Dic_p_prece = containers.Map(key_p,val_p);
            end

            % Communication
            if ~isempty(comm_prece)
                [key_c{1:size(comm_prece,1)-1,1}] = comm_prece{2:end,1};
                for i=1:size(comm_prece,1)-1
                    val_c{i,:}={comm_prece{i+1,2:size(comm_prece,2)}};
                end
                Dic_c_prece = containers.Map(key_c,val_c);
            end

            %Transportation
            %if ~isempty(trans_prece)
            if ~isempty(trans_task)
                [key_t{1:size(trans_prece,1)-1,1}] = trans_prece{2:end,1};
                for i=1:size(trans_prece,1)-1
                    val_t{i,:}={trans_prece{i+1,2:size(trans_prece,2)}};
                end
                Dic_t_prece = containers.Map(key_t,val_t);
            end
            clear key_p val_p key_c val_c key_t val_t
        end  
       
        %% Set up the dictionary of remaining tasks and the dictionary of remaining precedence for every system (SETP 4) 
        function[remain_task_p,remain_task_c,remain_task_t,remain_prece_p,remain_prece_c,remain_prece_t]=RemainTaskPrecedence(Dic_p_task,Dic_c_task,Dic_t_task,Dic_p_prece,Dic_c_prece,Dic_t_prece,Remain_schedule)
            
            % The remaining tasks of three systems
            %=== Power
             if ~isempty(Remain_schedule.Schedule_Power)
                  if size(Remain_schedule.Schedule_Power,2)<=size(Dic_p_task.keys,2)
                      for i=1:size(Remain_schedule.Schedule_Power,2)
                          tem_p=strsplit(Remain_schedule.Schedule_Power{i},'/');
                          key_p=tem_p(3);
                          val_p=Dic_p_task(char(key_p));val_p=[key_p,val_p];
                          [remain_task_p{i,1:size(val_p,2)}]=val_p{1,1:size(val_p,2)};
                      end
                  end
             else
                 remain_task_p=[];
             end
            %=== Comm
             if ~isempty(Remain_schedule.Schedule_Comm)
                  if size(Remain_schedule.Schedule_Comm,2)<=size(Dic_c_task.keys,2)
                      for i=1:size(Remain_schedule.Schedule_Comm,2)
                          tem_c=strsplit(Remain_schedule.Schedule_Comm{i},'/');
                          key_c=tem_c(3);
                          val_c=Dic_c_task(char(key_c));val_c=[key_c,val_c];
                          [remain_task_c{i,1:size(val_c,2)}]=val_c{1,1:size(val_c,2)};
                      end
                  end
             else
                 remain_task_c=[];
             end
            %=== Trans  
            if ~isempty(Remain_schedule.Schedule_Trans)
                  if size(Remain_schedule.Schedule_Trans,2)<=size(Dic_t_task.keys,2)
                      for i=1:size(Remain_schedule.Schedule_Trans,2)     
                          tem_t=strsplit(Remain_schedule.Schedule_Trans{i},'/');
                          key_t=tem_t(3);
                          val_t=Dic_t_task(char(key_t));
                          val_t=[key_t,val_t];
                          [remain_task_t{i,1:size(val_t,2)}]=val_t{1,1:size(val_t,2)};
                      end
                  end
            else
                 remain_task_t=[];
            end

          % The remaining precedences of three sys
             % Power
             if ~isempty(Remain_schedule.Schedule_Power) 
              if size(Remain_schedule.Schedule_Power,2)<=size(Dic_p_task.keys,2)
              for i=1:size(Remain_schedule.Schedule_Power,2)
              tem_p=strsplit(Remain_schedule.Schedule_Power{i},'/');
              key_p=tem_p(3);
              all_key_p{i}=char(key_p);
              val_p=Dic_p_prece(char(key_p));val_p=[key_p,val_p];
              [remain_prece_p{i,1:size(val_p,2)}]=val_p{1,1:size(val_p,2)};
              end
              mm=[{1},Dic_p_prece.keys]; mm{1}=[]; remain_prece_p=[mm;remain_prece_p];
              diff_p=setdiff(Dic_p_prece.keys,all_key_p);
              if ~isempty(diff_p)
              for j=1:size(diff_p,2)
              [~,n(j)]=find(strcmp(Dic_p_prece.keys,diff_p{j}));
              end
              remain_prece_p(:,n+1)=[];
              end
              task_num=remain_prece_p(2:end,1);
              for k=1:size(task_num,1)
                 [~,nn(k)]=find(strcmp(remain_prece_p(1,:),task_num{k})); 
              end
              remain_prece_p=remain_prece_p(:,[1,nn]);
              clear mm n nn task_num
              end
             else
                 remain_prece_p=[];
             end

              % Comm
             if ~isempty(Remain_schedule.Schedule_Comm)  
              if size(Remain_schedule.Schedule_Comm,2)<=size(Dic_c_task.keys,2)
              for i=1:size(Remain_schedule.Schedule_Comm,2)
              tem_c=strsplit(Remain_schedule.Schedule_Comm{i},'/');
              key_c=tem_c(3);
              all_key_c{i}=char(key_c);
              val_c=Dic_c_prece(char(key_c));val_c=[key_c,val_c];
              [remain_prece_c{i,1:size(val_c,2)}]=val_c{1,1:size(val_c,2)};
              end
              mm=[{1},Dic_c_prece.keys]; mm{1}=[]; remain_prece_c=[mm;remain_prece_c];
              % delete extra columns
              diff_c=setdiff(Dic_c_prece.keys,all_key_c);
              if ~isempty(diff_c)
              for j=1:size(diff_c,2)
              [~,n(j)]=find(strcmp(Dic_c_prece.keys,diff_c{j}));
              end
              remain_prece_c(:,n+1)=[];
              end
              % sort
              task_num=remain_prece_c(2:end,1); 
              for k=1:size(task_num,1)
                 [~,nn(k)]=find(strcmp(remain_prece_c(1,:),task_num{k})); 
              end
              remain_prece_c=remain_prece_c(:,[1,nn]);
              clear mm n nn task_num
              end
             else
                 remain_prece_c=[];
             end    
              
                  % Trans
             if ~isempty(Remain_schedule.Schedule_Trans) 
               if size(Remain_schedule.Schedule_Trans,2)<=size(Dic_t_task.keys,2)
              for i=1:size(Remain_schedule.Schedule_Trans,2)
              tem_t=strsplit(Remain_schedule.Schedule_Trans{i},'/');
              key_t=tem_t(3);
              all_key_t{i}=char(key_t);
              val_t=Dic_t_prece(char(key_t));val_t=[key_t,val_t];
              [remain_prece_t{i,1:size(val_t,2)}]=val_t{1,1:size(val_t,2)};
              end
              mm=[{1},Dic_t_prece.keys]; mm{1}=[]; remain_prece_t=[mm;remain_prece_t];
              % delete extra columns
              diff_t=setdiff(Dic_t_prece.keys,all_key_t);
              if ~isempty(diff_t)
              for j=1:size(diff_t,2)
              [~,n(j)]=find(strcmp(Dic_t_prece.keys,diff_t{j}));
              end
              remain_prece_t(:,n+1)=[];
              end
              % sort
              task_num=remain_prece_t(2:end,1);
              for k=1:size(task_num,1)
                 [~,nn(k)]=find(strcmp(remain_prece_t(1,:),task_num{k})); 
              end
              remain_prece_t=remain_prece_t(:,[1,nn]);
              clear mm n nn task_num
              end
             else
                 remain_prece_t=[];
             end  

        end

        %% === SampleActualDuration (STEP 4) 
        % sample the duration sample for every restoration task  
        function [TotalDamage] = SampleActualDuration(TotalDamage, SamplesUnifDur, isample)
            % ================================================================
            % function SampleActualDuration
            % Calculate actual restoration duration from its duration
            % distribution using the Latin Hypercube sampling (LHS)
            % ================================================================
            %--- Set Initial Values to Variables

            %--- assuming that every damage object requires at most 10 tasks (Ntsk=10),based on two functions of sampling: 
            % Library.CountTotalSampleNum and Library.GetDurationSample
            %--- icol: the index of column in the uniform sample matrix: SamplesUnifDur
            % icol also means the index of variable from lhsdesign.
            %icol = (iobj-1)*Ntsk + jtasknumber;
            Ntsk = 10;
            nobj = 0;
        
            %==== Field
          
            Power = TotalDamage{isample}{1}; 
            Commu = TotalDamage{isample}{2}; 
            Trans = TotalDamage{isample}{3};
            Dictionary = TotalDamage{isample}{4};
            Neighborhood = TotalDamage{isample}{5};
            
            Branch = Power{1};
            Bus = Power{2};
            Generator = Power{3};
            TransTower = Power{4};
            
            Centraloffice = Commu{1};
            CommunicationTower = Commu{2};
            Cellline = Commu{3};            
            
            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};       
            
            for ii = 1:length(Branch)
                if strcmp(Branch{ii}.Status, 'Damaged')
                    iobj = nobj+ii;
                    sum = 0;
                    tasks = Branch{ii}.taskUniqueIds;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]
                        icol = (iobj-1)*Ntsk + j; 
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode,temp.durationMax]));
                        
                        %%%%%%%%% Revised expression                       
                        if duration == 0
                            duration = 1;
                        end
                        temp.WorkingDays = [1, inv(24), inv(24*7), inv(24*7*4)] * duration;
                        sum = sum + temp.WorkingDays;
                    end
                    Branch{ii}.WorkingDays = sum;                  
                end
            end
            nobj = length(Branch);
            
            for ii = 1:length(Bus)
                if strcmp(Bus{ii}.Status, 'Damaged')
                    sum = 0;
                    tasks = Bus{ii}.taskUniqueIds;
                    %disp('Damaged Bus')
                    iobj = nobj + ii;
                    
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]
                        icol = (iobj-1)*Ntsk + j;
                        samples = SamplesUnifDur(isample,icol);
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin, temp.durationMode, temp.durationMax]));
                        
                        %%%%%%%%% Revised expression                       
                        if duration == 0
                            duration = 1;
                        end
                        temp.WorkingDays = [1, inv(24), inv(24*7), inv(24*7*4)] * duration;
                        
                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration)),max(1,ceil(duration/24)),max(1,ceil(duration/24/7)),max(1,ceil(duration/24/7/4))];
                        sum = sum + temp.WorkingDays;
                    end
                    Bus{ii}.WorkingDays = sum;
                end
            end
            nobj = nobj+length(Bus);
            
            for ii = 1:length(TransTower)
                if strcmp(TransTower{ii}.Status, 'Damaged')
                    sum = 0;
                    tasks = TransTower{ii}.taskUniqueIds;
                    % disp('TransTower')
                    iobj = nobj + ii;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]

                        icol = (iobj-1)*Ntsk + j; 
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode ,temp.durationMax]));
                        %%%%%%%%% Revised expression                       
                        if duration ==0
                            duration = 1;
                        end
                        temp.WorkingDays = [1, inv(24), inv(24*7), inv(24*7*4)] * duration;
                        
                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration)),max(1,ceil(duration/24)),max(1,ceil(duration/24/7)),max(1,ceil(duration/24/7/4))];
                        temp.WorkingDays = duration;
                        sum = sum + temp.WorkingDays;
                    end
                    TransTower{ii}.WorkingDays = sum;
                end
            end
            nobj = nobj+length(TransTower);
            
            % Generator ===================================================
            for ii = 1:length(Generator)
                if strcmp(Generator{ii}.Status, 'Damaged')
                    Generator{ii}.WorkingDays = 0; 
                end
            end
            nobj = nobj+length(Generator);
            
            for ii = 1:length(Centraloffice)
                if strcmp(Centraloffice{ii}.Status, 'Damaged')
                    sum = 0;
                    tasks = Centraloffice{ii}.taskUniqueIds; 
                    iobj = nobj + ii;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]

                        icol = (iobj-1)*Ntsk + j; 
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin, temp.durationMode, temp.durationMax]));
                        %%%%%%%%% Revised expression                       
                        if duration ==0
                            duration = 1;
                        end
                        temp.WorkingDays = [1, inv(24), inv(24*7), inv(24*7*4)] * duration;
                        
                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration)),max(1,ceil(duration/24)),max(1,ceil(duration/24/7)),max(1,ceil(duration/24/7/4))];
                        sum = sum + temp.WorkingDays;
                    end
                    Centraloffice{ii}.WorkingDays = sum;
                end
            end
            nobj = nobj+length(Centraloffice);
            
            
            for ii = 1:length(Cellline)
                if strcmp(Cellline{ii}.Status, 'Damaged')
                    sum = 0;
                    tasks = Cellline{ii}.taskUniqueIds;
                    iobj = nobj + ii;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]

                        icol = (iobj-1)*Ntsk + j; 
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode,temp.durationMax]));
                        %%%%%%%%% Revised expression                       
                        if duration ==0
                            duration = 1;
                        end
                        temp.WorkingDays = [1, inv(24), inv(24*7), inv(24*7*4)] * duration;
                        
                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration)),max(1,ceil(duration/24)),max(1,ceil(duration/24/7)),max(1,ceil(duration/24/7/4))];
                        sum = sum + temp.WorkingDays;
                    end
                    Cellline{ii}.WorkingDays = sum; 
                end
            end
            nobj = nobj+length(Cellline);
            
            for ii = 1:length(CommunicationTower)
                if strcmp(CommunicationTower{ii}.Status, 'Damaged')
                    sum = 0;
                    tasks = CommunicationTower{ii}.taskUniqueIds;
                    iobj = nobj + ii;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]
                        icol = (iobj-1)*Ntsk + j ;
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode,temp.durationMax]));
                        %%%%%%%%% Revised expression                       
                        if duration ==0
                            duration = 1;
                        end
                        temp.WorkingDays = [1, inv(24), inv(24*7), inv(24*7*4)] * duration;
                        
                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration)),max(1,ceil(duration/24)),max(1,ceil(duration/24/7)),max(1,ceil(duration/24/7/4))];
                        sum = sum + temp.WorkingDays;
                    end
                    CommunicationTower{ii}.WorkingDays = sum;
                end
            end
            nobj = nobj+length(CommunicationTower);
            
            for ii = 1:length(Road)
                if strcmp(Road{ii}.Status, 'Damaged')
                    sum = 0;
                    tasks = Road{ii}.taskUniqueIds;
                    iobj = nobj + ii;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]
                        icol = (iobj-1)*Ntsk + j;  
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode,temp.durationMax]));
                        %%%%%%%%% Revised expression                       
                        if duration ==0
                            duration = 1;
                        end
                        temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;
                        
                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))]; 
                        sum = sum + temp.WorkingDays;
                    end
                    Road{ii}.WorkingDays = sum;
                end
            end
            nobj = nobj+length(Road);
            
            
            for ii = 1:length(Bridge)     
                % Bridge that requires sub-component analysis
                if strcmp(Bridge{ii}.Status, 'Damaged') && Bridge{ii}.HasSub == 1
                        all_sum = 0;
                        for sub_index = 1:length(Bridge{ii}.ColumnSet)
                            tasks = Bridge{ii}.ColumnSet(sub_index).taskUniqueIds;
                            sum = 0;
                            iobj = nobj + sub_index;
                            for j = 1:length(tasks)
                                temp = Dictionary(tasks{j});
                                if iscell(temp)
                                    temp = temp{1};
                                end
                                %--- sample the duration of every restoration task, and ensure duration >0
                                % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                                % icol: the column index in SamplesUnifDur 
                                % WorkingDays = the duration value in different units of [hour, day, week, month]
                                
                                icol = (iobj-1)*Ntsk + j;  
                                samples = SamplesUnifDur(isample,icol); 
                                duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode,temp.durationMax]));
                                %%%%%%%%% Revised expression                       
                                if duration == 0
                                    duration = 1;
                                end
                                temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                                %%%%%%%%% Original expression
                                %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))]; 
                                sum = sum + temp.WorkingDays;                          
                            end
                            Bridge{ii}.ColumnSet(sub_index).WorkingDays = sum;
                            all_sum = all_sum + sum;
                        end
                        nobj = nobj+length(Bridge{ii}.ColumnSet);
            
                        
                        for sub_index = 1:length(Bridge{ii}.ColumnFoundSet)
                            tasks = Bridge{ii}.ColumnFoundSet(sub_index).taskUniqueIds;
                            sum = 0;
                            iobj = nobj + sub_index;
                            for j = 1:length(tasks)
                                temp = Dictionary(tasks{j});
                                if iscell(temp)
                                    temp = temp{1};
                                end
                                %--- sample the duration of every restoration task, and ensure duration >0
                                % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                                % icol: the column index in SamplesUnifDur 
                                % WorkingDays = the duration value in different units of [hour, day, week, month]

                                icol = (iobj-1)*Ntsk + j;  
                                samples = SamplesUnifDur(isample,icol); 
                                duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode,temp.durationMax]));
                                %%%%%%%%% Revised expression                       
                                if duration ==0
                                    duration = 1;
                                end
                                temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                                %%%%%%%%% Original expression
                                %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))];  
                                sum = sum + temp.WorkingDays;                              
                            end
                            Bridge{ii}.ColumnFoundSet(sub_index).WorkingDays = sum;
                            all_sum = all_sum + sum;
                        end
                        nobj = nobj+length(Bridge{ii}.ColumnFoundSet);
                        
                        for sub_index = 1:length(Bridge{ii}.AbutmentSet)
                            tasks = Bridge{ii}.AbutmentSet(sub_index).taskUniqueIds;
                            sum = 0;
                            iobj = nobj + sub_index;
                            for j = 1:length(tasks)
                                temp = Dictionary(tasks{j});
                                if iscell(temp)
                                    temp = temp{1};
                                end
                                %--- sample the duration of every restoration task, and ensure duration >0
                                % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                                % icol: the column index in SamplesUnifDur 
                                % WorkingDays = the duration value in different units of [hour, day, week, month]
         
                                icol = (iobj-1)*Ntsk + j;  
                                samples = SamplesUnifDur(isample,icol); 
                                duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode,temp.durationMax]));
                                %%%%%%%%% Revised expression                       
                                if duration ==0
                                    duration = 1;
                                end
                                temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                                %%%%%%%%% Original expression
                                %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))];  
                                sum = sum + temp.WorkingDays;
                            end
                            Bridge{ii}.AbutmentSet(sub_index).WorkingDays = sum;
                            all_sum = all_sum + sum;
                        end
                        nobj = nobj+length(Bridge{ii}.AbutmentSet);
                        
                        for sub_index = 1:length(Bridge{ii}.AbutmentFoundSet)
                            tasks = Bridge{ii}.AbutmentFoundSet(sub_index).taskUniqueIds;
                            sum = 0;
                            iobj = nobj + sub_index;
                            for j = 1:length(tasks)
                                temp = Dictionary(tasks{j});
                                if iscell(temp)
                                    temp = temp{1};
                                end
                                %--- sample the duration of every restoration task, and ensure duration >0
                                % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                                % icol: the column index in SamplesUnifDur 
                                % WorkingDays = the duration value in different units of [hour, day, week, month]
                                icol = (iobj-1)*Ntsk + j;  
                                samples = SamplesUnifDur(isample,icol); 
                                duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode ,temp.durationMax]));
                                %%%%%%%%% Revised expression                       
                                if duration ==0
                                    duration = 1;
                                end
                                temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                                %%%%%%%%% Original expression
                                %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))];  
                                sum = sum + temp.WorkingDays;
                                
                            end
                            Bridge{ii}.AbutmentFoundSet(sub_index).WorkingDays = sum;
                            all_sum = all_sum + sum; 
                        end
                        nobj = nobj+length(Bridge{ii}.AbutmentFoundSet);
                        
                        for sub_index = 1:length(Bridge{ii}.BearingSet)
                            tasks = Bridge{ii}.BearingSet(sub_index).taskUniqueIds;
                            sum = 0;
                            iobj = nobj + sub_index;
                            for j = 1:length(tasks)
                                temp = Dictionary(tasks{j});
                                if iscell(temp)
                                    temp = temp{1};
                                end
                                %--- sample the duration of every restoration task, and ensure duration >0
                                % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                                % icol: the column index in SamplesUnifDur 
                                % WorkingDays = the duration value in different units of [hour, day, week, month]
                                icol = (iobj-1)*Ntsk + j;  
                                samples = SamplesUnifDur(isample,icol); 
                                duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode ,temp.durationMax]));
                                %%%%%%%%% Revised expression                       
                                if duration ==0
                                    duration = 1;
                                end
                                temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                                %%%%%%%%% Original expression
                                %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))];  
                                sum = sum + temp.WorkingDays;  
                            end
                            Bridge{ii}.BearingSet(sub_index).WorkingDays = sum;
                            all_sum = all_sum + sum; 
                        end
                        nobj = nobj+length(Bridge{ii}.BearingSet);
                        
                        for sub_index = 1:length(Bridge{ii}.SlabSet)
                            tasks = Bridge{ii}.SlabSet(sub_index).taskUniqueIds;
                            sum = 0;
                            iobj = nobj + sub_index;
                            for j = 1:length(tasks)
                                temp = Dictionary(tasks{j});
                                if iscell(temp)
                                    temp = temp{1};
                                end
                                %--- sample the duration of every restoration task, and ensure duration >0
                                % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                                % icol: the column index in SamplesUnifDur 
                                % WorkingDays = the duration value in different units of [hour, day, week, month]
                                icol = (iobj-1)*Ntsk + j;  
                                samples = SamplesUnifDur(isample,icol); 
                                duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode ,temp.durationMax]));
                                %%%%%%%%% Revised expression                       
                                if duration ==0
                                    duration = 1;
                                end
                                temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                                %%%%%%%%% Original expression
                                %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))];  
                                sum = sum + temp.WorkingDays;        
                            end
                            Bridge{ii}.SlabSet(sub_index).WorkingDays = sum;
                            all_sum = all_sum + sum;
                        end
                        nobj = nobj+length(Bridge{ii}.SlabSet);
                        
                        Bridge{ii}.WorkingDays = all_sum;
           
                        
                % Bridge requiring object-level analysis, without any sub-component.        
                elseif strcmp(Bridge{ii}.Status, 'Damaged')
%                 if strcmp(Bridge{i}.Status, 'Damaged')

                    iobj = nobj + ii;
                    sum = 0;
                    tasks = Bridge{ii}.taskUniqueIds;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]
                        icol = (iobj-1)*Ntsk + j;  
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin,temp.durationMode ,temp.durationMax]));
                        %%%%%%%%% Revised expression                       
                        if duration ==0
                            duration = 1;
                        end
                        temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))];  
                        sum = sum + temp.WorkingDays;
                    end
                    Bridge{ii}.WorkingDays = sum;
                end
            end
            nobj = nobj+length(Bridge);
            
            for ii = 1:length(TrafficLight)
                if strcmp(TrafficLight{ii}.Status, 'Damaged')
                    iobj = nobj + ii;
                    sum = 0;
                    tasks = TrafficLight{ii}.taskUniqueIds;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        %--- sample the duration of every restoration task, and ensure duration >0
                        % samples size: (No. of row: Ndmgsample*NRun) * (No. of column: Nobj*Ntsk)
                        % icol: the column index in SamplesUnifDur 
                        % WorkingDays = the duration value in different units of [hour, day, week, month]
                        icol = (iobj-1)*Ntsk + j;  
                        samples = SamplesUnifDur(isample,icol); 
                        duration = round(Library.simulatervLHS(samples,temp.durationType, [temp.durationMin, temp.durationMode, temp.durationMax]));
                        %%%%%%%%% Revised expression                       
                        if duration ==0
                            duration = 1;
                        end
                        temp.WorkingDays = [24, 1, inv(7), inv(7*4)] * duration;

                        %%%%%%%%% Original expression
                        %temp.WorkingDays = [max(1,ceil(duration*24)),max(1,ceil(duration)),max(1,ceil(duration/7)),max(1,ceil(duration/7/4))]; 
                        sum = sum + temp.WorkingDays;
                    end
                    TrafficLight{ii}.WorkingDays = sum;
                end
            end
            nobj = nobj+length(TrafficLight); % for debug
            msg = strcat('----Function SampleActualDuration: Finish Sampling DurationSample No.', num2str(isample),'----');
            disp(msg);
            
        end
        
        %% find the minimal task duration (STEP 4)  
        function [dmin, TimeUnitUnq] = FindMinDuration(DurationPower, DurationComm, DurationTrans)
            
            tfactor = [inv(24), 1, 7, 28];
            dmin = zeros(1,4);
            dp = []; 
            dc = []; 
            dt = [];
            TimeUnit = [];
            
            if ~isempty(DurationPower)
                dp = DurationPower(:,5);
                ip = find(dp == min(dp));
                TimeUnit = [TimeUnit;DurationPower(:,6)];
            end
            if ~isempty(DurationComm)
                dc = DurationComm(:,5);
                ic = find(dc == min(dc));
                TimeUnit = [TimeUnit;DurationPower(:,6)];
            end
            if ~isempty(DurationComm)
                dt = DurationTrans(:,5);
                it = find(dt == min(dt));
                TimeUnit = [TimeUnit;DurationPower(:,6)];
            end
            if isempty(TimeUnit)
                TimeUnitUnq = [];
            else
                TimeUnitUnq = unique(TimeUnit);
            end
            
            if length(TimeUnitUnq)>1
                msg = 'Multiple time units exit. Stop.';
                disp(msg);
                dmin = [];
                return;
            else
                if isempty(dp)
                    dminp = 999999;
                else
                    dminp = dp(ip(1));
                end
                if isempty(dc)
                    dminc = 999999;
                else
                    dminc = dc(ic(1));
                end
                if isempty(dt)
                    dmint = 999999;
                else
                    dmint = dt(it(1));
                end
                d = [dminp, dminc, dmint];
                dmin = min(d) * tfactor(TimeUnitUnq); 
            end
            
            
        end
        
        %% Update the remain working days for tasks in the current working list (STEP 4)  
        function UpdateCurrentDuration(Current, Dictionary, MinTaskDuration)
            %===============================================================
            % function UpdateCurrentDuration
            % previously named "UpdateCurrentDuration"
            % This function reads the Current working list and subtracts
            % the Object.WorkingDays by Days (the shortest duration of all 
            % tasks in the Current list.)
            % MinTaskDuration is the task duratino in the unit of "day".
            %===============================================================
            tolerance = 1e-10;
            for ii = 1:size(Current,2)
                if isempty(Current{1,ii})
                    continue;
                end
                uniqueID = Library.getUniqueId(Current{1,ii}, 1);
                temp = Dictionary(uniqueID);
                if iscell(temp)
                    temp = temp{1};
                end
                %temp.WorkingDays(5) = temp.WorkingDays(5) - MinTaskDuration;
                %ind = temp.WorkingDays(6); % the index of selected time unit in the vector of [hour, day, week, month]
        
                %--- a matrix of duration conversion factors 
%                 tfactor = [1, inv(24), inv(24*7), inv(24*4*7);
%                     24, 1, inv(7), inv(7*4);
%                     24*7, 7, 1, inv(4);
%                     24*7*4, 7*4, 4, 1];
%                 ind = 2; % MinTaskDuration is in the unit of "day".

                tinvfactor = [24, 1, inv(7), inv(7*4)]; 
                indt = temp.WorkingDays(6);
                
                temp.WorkingDays(1:4) = temp.WorkingDays(1:4) - tinvfactor * MinTaskDuration;
                temp.WorkingDays(5) = temp.WorkingDays(5) - tinvfactor(indt) * MinTaskDuration;
                
                %---- Reset the task duration as zeros if the duration is 
                % very small (<tolerance). 
                check = le(temp.WorkingDays(1:5), tolerance); 
                
                if any(check)
                    temp.WorkingDays(1:5) = zeros(1,5);
                end
                
            end
        end

        %% Update the CurrentWork list related to the resource demand of every task (STEP 4)  
        % 1. Clean tasks that have completed from the "Current" working list;
        % 2. Relax the resource avaiability by adding resource demand that were taken by the tasks that have completed.
        function [return_Current, Max] = UpdateCurrentResource(Current, Dictionary, Max)
            %==============================================================
            % function UpdateCurrentResource
            % previously named "Clean"
            % This function updates the CurrentWorking list by updating the 
            % amount of resource used for excueting the task at the moment 
            % and clean out this task that has just been executed. 
            %==== Input
            % 1. Current = the CurrentWorking list
            % 2. Dictionary = DictionaryDamage output from Step 2
            % 3. Max = a vector that tracks the amount of resource used
            % before executing the task. This vector has the size of 1*(the
            % total numer of resource types)
            %==== Output
            % 1. return_Current = the updated CurrentWorking list
            % 2. Max = a vector that tracks the amount of resource used
            % after executing the task, by adding the original occupied resource, 
            % which is taken when executing tasks, back to the resource constraint 
            % at the new time step, when the task(s) is/are completed.
            %==============================================================
            %return_Current = {};
            Nr = length(Max); % number of resource type
            
            for itask = 1:size(Current,2)
                if ~isempty(Current{1,itask})
                    temp = Dictionary(Library.getUniqueId(Current{1,itask}, 1));
                    if iscell(temp)
                        temp = temp{1};
                    end
                    
                    if any(le(temp.WorkingDays(1:4), 0)) % task duration <=0
                        %--- for every type of resource:
                        % j: the index of resource type
                        for j = 1:Nr     
                            resourceNeed = Library.getResource(Current{1,itask}, Dictionary, j, 1);
                            %---The following line ensures that the previous occupied
                            %resource amount (in function AddCurrentWorking: there is a line of " Max_Resource(j) = Max_Resource(j) - resourceNeed(j);") 
                            %is now available in the resource constraint of "Max" after finishing this task,
                            % simplying by adding the occupied amount back to the resource constraint.
                            % "Max(j) = Max(j) - resourceNeed;" 
                            Max(j) = Max(j) + resourceNeed; 
                            
                        end
                        Current{1,itask} = [];
                        Current{2,itask} = [];
                    end
                end
            end
            return_Current = Current;
            
        end

        %% Calculate the Functionality of every system based on the Selected Functionality Metric(s) from the Input File (Input.txt)
        function [FunctionalityPower, FunctionalityComm,FunctionalityTrans] = ComputeSystemFunctionality(ActiveSystem, LinkDirection, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, Power_Set, Communication_Set, Transportation_Set, powerGraph, commGraph, transGraph, Dictionary)
            %=================================================================
            % function ComputeSystemFunctionality
            % This function computes the system fucntionality over time, 
            % based on the selected functionality metrics in the input file.
            %=================================================================    
            nQmetric = [length(Power_Func_Num), length(Comm_Func_Num), length(Trans_Func_Num)];
            
            %==== Set up initial values for system funtionality as -1
            FunctionalityPower = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
            FunctionalityComm = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
            FunctionalityTrans = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
            
            %==== flag variables to indicate which system is turned ON(1)/OFF(0)
            active_power = ActiveSystem(1);
            active_comm = ActiveSystem(2);
            active_trans = ActiveSystem(3);
            
            %==== Link direction indicator
            LinkDirectionChoicePower = LinkDirection(1);
            LinkDirectionChoiceComm = LinkDirection(2);
            LinkDirectionChoiceTrans = LinkDirection(3);
            
            flag = 1;
            
            while flag
                change_trans = 0;
                change_pow = 0;
                change_comm = 0;
                
                %%==== Power System   
                if active_power
                    isys = 1;
                    nq = nQmetric(isys);
                    for ii = 1:nq
                        index = Power_Func_Num(ii);
                        switch index
                           case 0
                               FunctionalityPower = {[],[],[],[]};

                           %==== Basic: percentage of functional bus (substation)    
                           case 1 
                               tmp = Library.Functionality_PowerBasic(Power_Set);%Power_Set{2} = Bus(substation)
                               if tmp ~= FunctionalityPower{1}
                                   FunctionalityPower{1} = tmp;
                                   change_pow = change_pow + 1;
                               end

                           %==== Weighted network (w = substation capacity in kV)
                           case 2 
                               SysNum = 1;
                               SystemSet = Power_Set; 
                               WeightNum = 1;
                               tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);
                               if tmp ~= FunctionalityPower{2}
                                   FunctionalityPower{2} = tmp;
                                   change_pow = change_pow + 1;
                               end

                           %==== Percentage of population with power service from substation    
                           case 3
                               WeightNum = 2;
                               tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);
                               if tmp ~= FunctionalityPower{3}
                                   FunctionalityPower{3} = tmp;
                                   change_pow = change_pow + 1;
                               end
                           
                           %==== Percentage of households with power service from substation    
                           case 4
                               WeightNum = 3;
                               tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);
                               if tmp ~= FunctionalityPower{4}
                                   FunctionalityPower{4} = tmp;
                                   change_pow = change_pow + 1;
                               end    
                               
                           %==== Percentage of neighborhoods with power service     
                           case 5
                               [totalPopulation, neighbourPowFunc, neighbourCommFunc, neighbourTransFunc] = Library.neighbourFunc(Dictionary);
                               if neighbourPowFunc ~= FunctionalityPower{5}
                                   FunctionalityPower{5} = neighbourPowFunc;
                                   change_pow = change_pow + 1;
                               end

                            % ======= Graph theory related functionality parameter =======
                            % Gaveragedegree (scalar) - average value of all node degree in the network (0,+inf)
                            % L (scalar) - characteristic path length (0,+inf)
                            % EGlob (scalar) - global efficiency [0,1]

                           case  6 || 7 || 8 || 9 || 10 % graph theory
                               result = Library.Functionality_GraphBasic(powerGraph,LinkDirectionChoicePower); 
                               FunctionalityPower{6} = result(1);
                               FunctionalityPower{7} = result(2);
                               FunctionalityPower{8} = result(3);
                               FunctionalityPower{9} = result(4);
                               FunctionalityPower{10} = result(5);
%                                FunctionalityPower{9} = result(6);
%                                FunctionalityPower{10} = result(7);
                               change_pow = change_pow + 1; % what is this for?

                        end % switch index
                    end % for ii
                end % if actice_power

                %%==== Communication System 
                if active_comm         
                    isys = 2;
                    nq = nQmetric(isys);
                    for ii = 1:nq
                        index = Comm_Func_Num(ii);
                        switch index
                            case 0
                               FunctionalityComm = {[],[],[],[]};
                            %==== Percentage of open central office & communication tower    
                            case 1
                                tmp = Library.Functionality_CommunicationBasic(Communication_Set);
                                if tmp ~= FunctionalityComm{1}
                                    FunctionalityComm{1} = tmp;
                                    change_comm = change_comm + 1;
                                end

                            %==== Weighted network    (equal weight for every comm tower)
                            case 2 
                                SysNum = 2;
                                SystemSet = Communication_Set;
                                WeightNum = 1;
                                tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);                      
                                if tmp ~= FunctionalityComm{2}
                                    FunctionalityComm{2} = tmp;
                                    change_comm = change_comm + 1;
                                end 

                            %==== Weighted network    (% of population served with comm service from comm tower)
                            case 3 
                                SysNum = 2;
                                SystemSet = Communication_Set;
                                WeightNum = 2;
                                tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);                      
                                if tmp ~= FunctionalityComm{3}
                                    FunctionalityComm{3} = tmp;
                                    change_comm = change_comm + 1;
                                end 
                                
                            %==== Weighted network    (% of Household Served with comm service from comm tower)
                            case 4 
                                SysNum = 2;
                                SystemSet = Communication_Set;
                                WeightNum = 3;
                                tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);                      
                                if tmp ~= FunctionalityComm{4}
                                    FunctionalityComm{4} = tmp;
                                    change_comm = change_comm + 1;
                                end 
                                
                            %==== Percentage of neightborhoods with the serice from central office (such as internet)     
                            case 5    
                                [totalPopulation,neighbourPowFunc,neighbourCommFunc,neighbourTransFunc] = Library.neighbourFunc(Dictionary);
                                if neighbourCommFunc ~= FunctionalityComm{5}
                                    FunctionalityComm{5} = neighbourCommFunc;
                                    change_comm = change_comm + 1;
                                end

                            %======= Graph theory related functionality parameter
                            %   Gaveragedegree (scalar) - average value of all node degree in the network (0,+inf)
                            %   L (scalar) - characteristic path length (0,+inf)
                            %   EGlob (scalar) - global efficiency [0,1]
                            case  6 || 7 || 8 || 9 || 10 % graph theory
                                result = Library.Functionality_GraphBasic(commGraph,LinkDirectionChoiceComm); 
                                FunctionalityComm{6} = result(1);
                                FunctionalityComm{7} = result(2);
                                FunctionalityComm{8} = result(3);
                                FunctionalityComm{9} = result(4);
                                FunctionalityComm{10} = result(5);
%                                 FunctionalityComm{9} = result(6);
%                                 FunctionalityComm{10} = result(7);
                                change_comm = change_comm + 1; 
                        end % switch index
                    end % for ii
                end % if actice_comm
              
                %%==== Transportation System                  
                if active_trans
                    isys = 3;
                    nq = nQmetric(isys);
                    for ii = 1:nq
                        index = Trans_Func_Num(ii);
                        switch index
                            case 0
                                FunctionalityTrans = {[],[],[],[]};                  
                            %==== Basic: percentage of open bus (substation)    
                            case 1
                                tmp = Library.Functionality_TransportationBasic(Transportation_Set);
                                if tmp ~= FunctionalityTrans{1}
                                    FunctionalityTrans{1} = tmp;
                                    change_trans = change_trans + 1;
                                end

                            %==== Weighted network   (Length) 
                            case 2 
                                SysNum = 3;
                                SystemSet = Transportation_Set; 
                                WeightNum = 1;
                                tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);                     
                                if tmp ~= FunctionalityTrans{2}
                                    FunctionalityTrans{2} = tmp;
                                    change_trans = change_trans + 1;
                                end
                            %==== Weighted network   (AADT) 
                            case 3 
                                SysNum = 3;
                                SystemSet = Transportation_Set; 
                                WeightNum = 2;
                                tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);                     
                                if tmp ~= FunctionalityTrans{3}
                                    FunctionalityTrans{3} = tmp;
                                    change_trans = change_trans + 1;
                                end
                            %==== Weighted network   (AADT*Length) 
                            case 4 
                                SysNum = 3;
                                SystemSet = Transportation_Set; 
                                WeightNum = 3;
                                tmp = Library.Functionality_WeightNetwork(SysNum,SystemSet,WeightNum);                     
                                if tmp ~= FunctionalityTrans{4}
                                    FunctionalityTrans{4} = tmp;
                                    change_trans = change_trans + 1;
                                end

                            %==== Percentage of neighborhoods with the road access      
                            case 5    
                                [totalPopulation,neighbourPowFunc,neighbourCommFunc,neighbourTransFunc] = Library.neighbourFunc(Dictionary);
                                if neighbourTransFunc ~= FunctionalityTrans{5}
                                    FunctionalityTrans{5} = neighbourTransFunc;
                                    change_trans = change_trans + 1;
                                end

                            % ======= Graph theory related functionality parameter
                            %   Gaveragedegree (scalar) - average value of all node degree in the network (0,+inf)
                            %   L (scalar) - characteristic path length (0,+inf)
                            %   EGlob (scalar) - global efficiency [0,1]
                            case  6 || 7 || 8 || 9 || 10 % graph theory
                                result = Library.Functionality_GraphBasic(transGraph,LinkDirectionChoiceTrans); 
                                FunctionalityTrans{6} = result(1);
                                FunctionalityTrans{7} = result(2);
                                FunctionalityTrans{8} = result(3);
                                FunctionalityTrans{9} = result(4);
                                FunctionalityTrans{10} = result(5);
%                                 FunctionalityTrans{9} = result(6);
%                                 FunctionalityTrans{10} = result(7);
                                change_trans = change_trans + 1; 
                        end
                    end
                end

                if change_pow == 0 && change_trans == 0 && change_comm == 0
                    flag = 0;
                end
              
            end
        end        
        
        %% Calculate the repair time from lognormal distribution
        % this is an older function, not used anymore. 
        function Time_Component = RepairTime(Component)
            m=Component.RecoveryMatrix(1);
            v=Component.RecoveryMatrix(2);
            mu = log((m^2)/sqrt(v+m^2));
            sigma = sqrt(log(v/(m^2)+1));
            Time_Component=round(lognrnd(mu,sigma));
            if(Time_Component <= 0)
                Time_Component = 1;
            end
        end
        
        %% Find the minimum task duration for all tasks in the current working list at the current time step (STEP 4) 
        function MinTaskDuration = FindMinDays(Current_Power, Current_Comm, Current_Trans, Pow, Comm, Trans, Dictionary)
            %===============================================================
            % function FindMinDays
            % This function finds out the shortest task duration among all tasks
            % in the Current Working list. 
            %===============================================================
            %=== Field
            Branch= Pow{1};
            Bus= Pow{2};
            Generator= Pow{3};
            TransTower = Pow{4};         
            
            Centraloffice = Comm{1};
            CommunicationTower = Comm{2};
            Cellline = Comm{3};
            
            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};
            
            flag = 0;
            day = 99999999;
            MinTaskDuration = day;
            
            %--- Index of time unit in the vector of [hour, day, week, month]
            indt = 0;
            %--- A matrix of task duration conversion factors at different
            % time units [hour, day, week, month (considering 1 month=28 days)]
            tfactor = [1, inv(24), inv(24*7), inv(24*4*7);
                    24, 1, inv(7), inv(7*4);
                    24*7, 7, 1, inv(4);
                    24*7*4, 7*4, 4, 1];
            
            
            %=== Power System
            % for every task in the Current list, evlaute the Object
            for ii = 1:size(Current_Power,2) 
                if ~isempty(Current_Power{1,ii})
                    flag = 1;
                    tem = strsplit(Current_Power{1,ii},'/');
                    
                    %--- Branch
                    if strcmp(tem(1), 'Branch') && strcmp(Branch{str2double(tem(2))}.Status, 'Damaged')                  
                        % length(tem) == 4
                        % Meaning tem is something like 'Branch/2/Working/dummy'
                        if length(tem) == 4
                            task = Branch{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];
                            
                        % length(tem) == 5
                        % Meaning the task of "tem" is something like 'Branch/2/Task101/Working/dummy'    
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day, temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0 % Meaning the task of "tem" has completed the restoration
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays: ',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                        
                    %--- Bus     
                    elseif strcmp(tem(1), 'Bus') && strcmp(Bus{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Bus{str2double(tem(2))}; 
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)]; 
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day, temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays:',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                        
                    %--- TransmissionTower    
                    elseif strcmp(tem(1), 'TransmissionTower') && strcmp(TransTower{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = TransTower{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day, temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays:',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                        
                    %--- Generator    
                    elseif strcmp(tem(1), 'Generator') && strcmp(Generator{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Generator{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                            % Asumming Generator does not get damaged,
                            % requiring Generator{indexGenerator}.workingday = 0
                            % If Generator{indexGenerator}.workingday < 0,
                            % display this indexGenerator
                            if lt(Generator{str2double(tem(2))}.WorkingDays(5), 0)
                                s1 = tem{1};
                                s2 = tem{2};
                                ss = strcat('Function FindMinDays:',s1,s2, '.WorkingDays<=0');
                                disp(ss);
                            end
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays:',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                    else
                        msg = strcat('Function FindMinDays Error: This Task/Object is not Recognized: ', tem,' in the Power System.');
                        disp(msg); 
                    end
                end
                
            end
            
            %=== Communication System
            for ii = 1:size(Current_Comm,2)
                if ~isempty(Current_Comm{1,ii})
                    flag = 1;
                    tem = strsplit(Current_Comm{1,ii},'/');
                    
                    %--- Central Office
                    if strcmp(tem(1), 'Centraloffice') && strcmp(Centraloffice{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Centraloffice{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays: ',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end

                    %--- Communication Tower    
                    elseif strcmp(tem(1), 'CommunicationTower')&& strcmp(CommunicationTower{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = CommunicationTower{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays: ',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                        
                    %--- Cell Line    
                    elseif strcmp(tem(1), 'Cellline')&& strcmp(Cellline{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Cellline{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];
                            
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays: ',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                    
                    else
                        msg = strcat('Function FindMinDays Error: This Task/Object is not Recognized: ', tem,' in the Communication System.');
                        disp(msg);                         
                    end
                    
                    
                end
            end
            
            %=== Transportation System
            for ii = 1:size(Current_Trans,2)
                if ~isempty(Current_Trans{1,ii})
                    flag = 1;
                    tem = strsplit(Current_Trans{1,ii},'/');
                    
                    %--- Road
                    if strcmp(tem(1), 'Road')&& strcmp(Road{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Road{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays: ',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                        
                    %--- Traffic Light    
                    elseif strcmp(tem(1), 'TrafficLight')&& strcmp(TrafficLight{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            day = min(day, TrafficLight{str2double(tem(2))}.WorkingDays(5));
                            indt = [indt, TrafficLight{str2double(tem(2))}.WorkingDays(6)];
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                        
                    %--- Bridge    
                    %%% Bridge Object without the sub-component analysis 
                    elseif strcmp(tem(1), 'Bridge')&& strcmp(Bridge{str2double(tem(2))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                s1 = tem{1};
                                s2 = tem{2};
                                s3 = tem{3};
                                ss = strcat('Function FindMinDays: ',s1,'/',s2,'/',s3, '.WorkingDays<=0'); 
                                disp(ss);
                            end
                        end
                        
                    %%% Bridge Object with the sub-component analysis    [These lines need to be revised and improved!]
                    %- Bridge Abutment
                    elseif strcmp(tem(1), 'Bridge') && strcmp(tem(3), 'Abutment')&& strcmp(Bridge{str2double(tem(2))}.AbutmentSet{str2double(tem(4))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.AbutmentSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                        
                    %- Bridge Abutment Foundation   
                    elseif  strcmp(tem(1), 'Bridge') && strcmp(tem(3), 'AbutmentFoundation')&& strcmp(Bridge{str2double(tem(2))}.AbutmentFoundSet{str2double(tem(4))}.Status, 'Damaged')        
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.AbutmentFoundSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                        
                    %- Bridge Approach Slab     
                    elseif  strcmp(tem(1), 'Bridge') && strcmp(tem(3),'ApproachSlab') && strcmp(Bridge{str2double(tem(2))}.SlabSet{str2double(tem(4))}.Status, 'Damaged')  
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.SlabSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                        
                    %- Bridge Bearing   
                    elseif  strcmp(tem(1), 'Bridge') && strcmp(tem(3), 'Bearing') && strcmp(Bridge{str2double(tem(2))}.BearingSet{str2double(tem(4))}.Status, 'Damaged')     
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.BearingSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                        
                    %- Bridge Column     
                    elseif  strcmp(tem(1), 'Bridge') && strcmp(tem(3), 'Column') && strcmp(Bridge{str2double(tem(2))}.ColumnSet{str2double(tem(4))}.Status, 'Damaged')
                        
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.ColumnSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];

                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays:',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                      
                    %- Bridge Column Foundation    
                    elseif strcmp(tem(1), 'Bridge') && strcmp(tem(3), 'ColumnFoundation') && strcmp(Bridge{str2double(tem(2))}.ColumnFoundSet{str2double(tem(4))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.ColumnFoundSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];
                            
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                    
                    %- Bridge Deck    
                    elseif strcmp(tem(1), 'Bridge') && strcmp(tem(3), 'Deck')&& strcmp(Bridge{str2double(tem(2))}.DeckSet{str2double(tem(4))}.Status, 'Damaged')  
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.DeckSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];
                            
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end
                        
                    %- Bridge Girder    
                    elseif strcmp(tem(1), 'Bridge') && strcmp(tem(3), 'Girder')&& strcmp(Bridge{str2double(tem(2))}.GirderSet{str2double(tem(4))}.Status, 'Damaged')
                        if length(tem) == 4
                            task = Bridge{str2double(tem(2))}.GirderSet{str2double(tem(4))};
                            day = min(day, task.WorkingDays(5));
                            indt = [indt, task.WorkingDays(6)];
                            
                        elseif length(tem) == 5
                            temp = Dictionary(tem{3});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            day = min(day,temp.WorkingDays(5));
                            indt = [indt, temp.WorkingDays(6)];
                            if temp.WorkingDays(5) <= 0
                                msg = strcat('Function FindMinDays: ',temp, '.WorkingDays<=0 [',tem{3},']');
                                disp(msg);
                            end
                        end                        
                        

                    else
                        msg = strcat('Function FindMinDays Error: This Task/Object is not Recognized: ', tem,' in the Transportation System.');
                        disp(msg); 
                    end
                    
                end
            end
            
            if day == 99999999
                day = -1;
            end
            
            if (isempty(day)|| day < 0) && flag == 1
                disp('Function FindMinDays Error: MinimalWorkingDays is either empty or <0.');
            else %flag == 0
                % convert the minimal duration variable of day in different time
                % unit to the duration in the time unit of "day".
                ind = unique(indt(2:end));
                
                % flag =0 & day = -1 & ind = [], meaning the Current list for
                % every system is empty, all restoration tasks are finished. 
                if day == -1 && isempty(ind)
                    MinTaskDuration = 0;
                elseif length(ind)==1
                    MinTaskDuration = day*tfactor(ind,2);
                  
                else
                    disp('Function FindMinDays Error: inconsistent time unit is used. Check Functions of "AddCurrentWorking" and "FindMinDays"');
                end
                
            end
            
        end
        
        %% Add back node or edge back to the graph if fixed (STEP 4) 
        function G = addTransGraph(Trans_Set, G, index)
            %===============================================================
            % function addTransGraph
            % This function generates a new graph by adding node/edge to the 
            % old graph, when a damaged Object is fixed/restored.  [Transportation]
            %===============================================================
            hash = containers.Map('KeyType','double','ValueType','char');
            road_Set = Trans_Set{1};
            roadnode_Set = Trans_Set{4};
            for ii = 1:length(roadnode_Set)
                hash(roadnode_Set{ii}.nodeID) = roadnode_Set{ii}.uniqueID;
            end
            s = hash(road_Set{index}.Start_Node);
            t = hash(road_Set{index}.End_Node);
            G = addedge(G, s, t, road_Set{index}.Length);
        end
        
        %% Add back node or edge back to the graph if fixed [I have questions about this function]
        function commGraph = addCommGraph(Comm, commGraph, index, Dictionary)
            %===============================================================
            % function addCommGraph
            % This function generates a new graph by adding node/edge to the 
            % old graph, when a damaged Object is fixed/restored. [Communication]
            %===============================================================
            Cellline = Comm{3};
            s = Cellline{index}.connectedObj1;
            t = Cellline{index}.connectedObj2;
            temp1 = Dictionary(Cellline{index}.connectedObj1);
            temp2 = Dictionary(Cellline{index}.connectedObj2);
            temp1 = temp1{1};
            temp2 = temp2{1};
            weight = distance(temp1.Location(1), temp1.Location(2), temp2.Location(1),temp2.Location(2),referenceSphere('earth','km'));
            commGraph = addedge(commGraph, s, t, weight);
        end
        
        %% Add back node or edge back to the graph if fixed [I have questions about this function]
        function powerGraph = addPowerGraph(Power, powerGraph, index,Dictionary)
            %===============================================================
            % function addPowerGraph
            % This function generates a new graph by adding node/edge to the 
            % old graph, when a damaged Object is fixed/restored.  [Power]
            %===============================================================
            Branch = Power{1};
            s = Branch{index}.connectedObj1;
            t = Branch{index}.connectedObj2;
            temp1 = Dictionary(Branch{index}.connectedObj1);
            temp2 = Dictionary(Branch{index}.connectedObj2);
            temp1 = temp1{1};
            temp2 = temp2{1};
            weight = distance(temp1.Location(1), temp1.Location(2), temp2.Location(1),temp2.Location(2),referenceSphere('earth','km'));
            powerGraph = addedge(powerGraph, s, t, weight);
        end
        

        %% Compute Functionality based on graph status
        function [result] = Functionality_GraphBasic(Graph, DirectionVariable)
            [Graph,Gaveragedegree,L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = Library.ComputeGraphNetwork(Graph, DirectionVariable);
            
            result(1) = Gaveragedegree;
            result(2) = L;
            result(3) = EGlob;
            result(4) = CClosed;
            result(5) = ELocClosed;
            result(6) = COpen;
            result(7) = ELocOpen;
        end
        
        %% compute the graph properties of a network
        function [G,Gaveragedegree,L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = ComputeGraphNetwork(G,DirectionVariable)
            %===================================================================== 
            % function computeGraphNetwork: 
            % computes graph functionality: (1) graph properties
            % related to nodes, such as degree, centrality (betweenness and
            % closeness, etc.) and (2) graph properties related to the network
            % (suing function of graphProperties ).
            % G: graph
            % DirectionVariable: 1=graph(directional), 0=digraph(not directional)
            % Adjacency: adjacency matrix
            %=====================================================================

            EdgeLength = G.Edges.Weight;
            Adjacency = adjacency(G);
 
            switch DirectionVariable
                 case 0 % not directional in the network graph 
                    G.Nodes.degree = degree(G);
                    G.Nodes.closeness = centrality(G,'closeness','Cost',EdgeLength);
                    G.Nodes.betweenness = centrality(G,'betweenness','Cost',EdgeLength);
                    G.Nodes.ShortestDistance = distances(G,'Method','positive');      
                    % system topology functionality 
                    Gaveragedegree = mean(G.Nodes.degree);
                    [L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = Library.ComputeGraphProperty(Adjacency);   
                case 1 % directional in the network graph 
                    G.Nodes.indegree = centrality(G,'indegree');
                    G.Nodes.outdegree = centrality(G,'outdegree');
                    G.Nodes.incloseness = centrality(G,'incloseness','Cost',EdgeLength);
                    G.Nodes.outcloseness = centrality(G,'outcloseness','Cost',EdgeLength);
                    G.Nodes.betweenness = centrality(G,'betweenness','Cost',EdgeLength);
                    G.Nodes.ShortestDistance = distances(G,'Method','positive');                            
                    % system topology functionality 
                    Gaveragedegree = mean(G.Nodes.indegree);
                    [L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = Library.ComputeGraphProperty(Adjacency);
            end         

        end        
        
        %% compute graph properties
        function [ L, EGlob, CClosed, ELocClosed, COpen, ELocOpen ] = ComputeGraphProperty(A)
            %===================================================================== 
            % function ComputeGraphProperty
            % This function computes the following graph properties at the network level.
            % It computes properties of a graph from its adjacency matrix
            % usage: [L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = ComputeGraphProperty(A);
            %-----Input-----
            %   A (nxn) - adjacency matrix of a graph G
            %-----Output-----
            %   L (scalar) - characteristic path length
            %   EGlob (scalar) - global efficiency
            %   CClosed (scalar) - clustering coefficient (closed neighborhood)
            %   ELocClosed (scalar) - local efficiency (closed neighborhood)
            %   COpen (scalar) - clustering coefficient (open neighborhood)
            %   ELocOpen (scalar) - local efficiency (open neighborhood)
            %-----Source-----
            % Sighlt revised, based on the function of graphProperties from the following
            % author: Nathan D. Cahill
            % email: nathan.cahill@rit.edu
            % date: 10 April 2014
            % https://www.mathworks.com/matlabcentral/fileexchange/46084-characteristic-path-length-global-and-local-efficiency-and-clustering-coefficient-of-a-graph?focused=5289221&tab=function
            %===================================================================== 
            
            %==== get number of vertices (nodes)
            n = size(A,1);

            %==== shortest path distances between each node
            D = graphallshortestpaths(A);

            %==== characteristic path length (average path length)
            %L = sum(D(:))/(n*(n-1));
            sumD = sum(D(~isinf(D(:))));
            L = sumD/(n*(n-1));

            %==== global efficiency (average efficiency of a network)
            %EGlob = (sum(sum(1./(D+eye(n)))) - n)/(n*(n-1));
            D2 = 1./(D+eye(n));
            sumD2 = sum(D2(:));
            EGlob = (sumD2 - n)/(n*(n-1));

            %==== subgraphs of G induced by the neighbors of each vertex
            [MClosed,kClosed,MOpen,kOpen] = Library.subgraphs(A);

            %==== local clustering coefficient in each subgraph
            [CLocClosed,CLocOpen] = deal(zeros(n,1));
                for ii = 1:n
                    CLocClosed(ii) = sum(MClosed{ii}(:))/...
                        (numel(kClosed{ii})*(numel(kClosed{ii})-1));
                    CLocOpen(ii) = sum(MOpen{ii}(:))/...
                        (numel(kOpen{ii})*(numel(kOpen{ii})-1));
                end

            %==== clustering coefficients
            %CClosed = mean(CLocClosed);
            %COpen = mean(CLocOpen);
            CClosed = nanmean(CLocClosed);
            COpen = nanmean(CLocOpen);

            %==== local efficiency of each subgraph
            [ELocSGClosed,ELocSGOpen] = deal(zeros(n,1));
                for ii = 1:n
                    % distance matrix and number of vertices for current subgraph
                    DSGClosed = graphallshortestpaths(MClosed{ii});
                    DSGOpen = graphallshortestpaths(MOpen{ii});
                    nSGClosed = numel(kClosed{ii});
                    nSGOpen = numel(kOpen{ii});
                    % efficiency of current subgraph
                    ELocSGClosed(ii) = (sum(sum(1./(DSGClosed+eye(nSGClosed)))) - nSGClosed)/...
                        (nSGClosed*(nSGClosed-1));
                    ELocSGOpen(ii) = (sum(sum(1./(DSGOpen+eye(nSGOpen)))) - nSGOpen)/...
                        (nSGOpen*(nSGOpen-1));
                end

            %==== local efficiency of graph
            %ELocClosed = mean(ELocSGClosed);
            %ELocOpen = mean(ELocSGOpen);
            ELocClosed = nanmean(ELocSGClosed);
            ELocOpen = nanmean(ELocSGOpen);

        end
       
        %% compute adjacency matrices for each vertex in a graph
        function [MClosed,kClosed,MOpen,kOpen] = subgraphs(A)
            %=================================================================
            % function subgraphs: 
            % This function computes adjacency matrices for each vertex in a graph
            % usage: [MClosed,kClosed,MOpen,kOpen] = subgraphs(A);
            %---------Input---------
            %   A - (nxn) adjacency matrix of a graph G
            %---------Output---------
            %   MClosed, MOpen - (nx1) cell arrays containing adjacency matrices of the 
            %       subgraphs corresponding to neighbors of each vertex. For example, 
            %       MClosed{j} is the adjacency matrix of the subgraph of G 
            %       corresponding to the closed neighborhood of the jth vertex of G, 
            %       and kClosed{j} is the list of vertices of G that are in the 
            %       subgraph (and represent the corresponding rows/columns of 
            %       MClosed{j})     
            % author: Nathan D. Cahill
            % email: nathan.cahill@rit.edu
            % date: 10 Apr 2014
            %=================================================================

        % number of vertices
        n = size(A,1);

        % initialize indices of neighboring vertices, and adjacency matrices
        [kClosed,kOpen] = deal(cell(n,1));
        [MClosed,MOpen] = deal(cell(n,1));

        % loop through each vertex, determining its neighbors
            for i=1:n

                % find indices of neighbors of vertex v_i
                k1 = find(A(i,:)); % vertices with edge beginning at v_i
                k2 = find(A(:,i)); % vertices with edge ending at v_i
                kOpen{i} = union(k1,k2); % vertices in neighborhood of v_i
                kClosed{i} = union(kOpen{i},i);

                % extract submatrix describing adjacency of neighbors of v_i
                MOpen{i} = A(kOpen{i},kOpen{i});
                MClosed{i} = A(kClosed{i},kClosed{i});

            end

        end 
   
    
        %% Compute Functionality based on the weighted network 
        function SystemFunctionality = Functionality_WeightNetwork(SysNum, SystemSet, WeightNum) 
            %================================================================================
            % function Functionality_WeightNetwork
            % This function computes the functionality of every system as a
            % weighted network. The weights vary from system to system. 
            % (1) power system: weight =1 the maximum voltage of every bus (substation);
            % weight=2 the population served by every bus (substation)
            % weight=3 the number of household served by every bus
            % (2) comm system: weight = 1 equal weight for every communciation tower
            % weight = 2 the population served by every communication tower
            % weight = 3 the number of household served by every communication tower
            % (3) trans system: weight = 1-length/2-AADT/3-length*AADT for every road segment
            %-----Input-----
            % SysNum = A integer variable represeting which system is
            % selected to compute functionality.
            % SystemSet = the SystemSet including ObjectSets in this system
            %-----Output-----
            % SystemFunctionality = a numerical value of system functionality,
            % with the value between 0 and 1.
            %================================================================================
            count = 0;
            switch SysNum
                
                %==== Power System
                case 1
                    Power_Set = SystemSet;
                    Bus_Set = Power_Set{2}; 
                    n = length(Bus_Set);
                    w = zeros(1,n); 
                    for ii = 1:n
                        if WeightNum == 1
                            w(ii) = Bus_Set{ii}.Capacity; % consider the weight as bus capacity (maximum voltage)
                        elseif WeightNum ==2
                            w(ii) = Bus_Set{ii}.PopulationServed;
                        elseif WeightNum ==3
                            w(ii) = Bus_Set{ii}.HouseholdServed;                                
                        end
                        if strcmp(Bus_Set{ii}.Status,'Open')
                            count = count + w(ii);
                        else 
                            count = count;
                        end
                    end
                
                %==== Communication    
                case 2
                    Communication_Set = SystemSet;
                    CommunicationTower_Set = Communication_Set{2}; 
                    n = length(CommunicationTower_Set);
                    w = zeros(1,n); 
                    for ii = 1:n
                        if WeightNum == 1
                            w(ii) = 1; % consider the weight as the equal capacity of every communication tower
                        elseif WeightNum ==2
                            w(ii) = CommunicationTower_Set{ii}.PopulationServed;
                        elseif WeightNum ==3
                            w(ii) = CommunicationTower_Set{ii}.HouseholdServed;
                        end
                        
                        if strcmp(CommunicationTower_Set{ii}.Status,'Open')
                            count = count + w(ii);
                        else 
                            count = count;
                        end
                    end

                %==== Transportation
                case 3
                    Transportation_Set = SystemSet; 
                    Road_Set = Transportation_Set{1};
                    n = length(Road_Set);
                    w = zeros(1,n); 
                    for ii = 1:n
                        if WeightNum == 1
                            w(ii) = Road_Set{ii}.Length; % (1)consider the weight as as road segment length
                        elseif WeightNum == 2
                            w(ii) = Road_Set{ii}.AADT; %(2) consider the weight as as traffic (AADT) 
                        elseif WeightNum == 3
                            w(ii) = Road_Set{ii}.Length * Road_Set{ii}.AADT; % (3) consider the weight as as traffic*length                    
                        end
                        
                        if strcmp(Road_Set{ii}.Status, 'Open')
                            count = count + w(ii);
                        else 
                            count = count;
                        end
                    end       
            end
            
            SystemFunctionality = count/sum(w);
                    
        end
            

        
        %% Functionality of the Power system = Percentage of Buses that have electricity
        % Calculate the functionality of power
        function Function = Functionality_PowerBasic(Power)
            %================================================================================
            % function Functionality_PowerBasic            
            % This function computes the functionality of the Power System
            % as the percentage of Open Bus(substation) in the power
            % network.
            %================================================================================
            
            Bus = Power{1};  
            n =length(Bus);
            Num_Hous_Open=0;
            for ii = 1:n 
                if strcmp(Bus{ii}.Status,'Open')
                    Num_Hous_Open=Num_Hous_Open+1;
                else
                end
            end
            Function=Num_Hous_Open/n;
        end
        
        %% Functionality of ths Transportation system = Connectivity of each Points on the Map
        % Calculate the functionality of transportation
        function Function = Functionality_TransportationBasic(Trans)
            %================================================================================
            % function Functionality_TransportationBasic            
            % This function computes the functionality of the Transportation System
            % as the percentage of Open Road (road segment) in the
            % transportation network.
            %================================================================================
            score = 0;
            
            %=== functionality = the percentage of road in the "open"/"stoped" status
            Road = Trans{1}; 
            n = length(Road); 
            for ii = 1:length(Road)
                if strcmp(Road{ii}.Status, 'Open')
                     score = score + 1;
                elseif strcmp(Road{ii}.Status, 'Stoped')
                     score = score + 0.5;
                end
            end
            
%             %=== Alternatively, 
%             % functionality = the percentage of bridge in the "open" status
%             Bridge = Trans{2};
%             n = length(Bridge); 
%             for i = 1:length(Bridge)
%                 if strcmp(Bridge{i}.Status, 'Open')
%                     score = score + 1;
%                 else
%                 end
%             end
            
            Function = score/n;
        end
        
        %% Functionality of the Communication system 
        % Calculate the functionality of communication
        function Function = Functionality_CommunicationBasic(Comm)
            %================================================================================
            % function Functionality_TransportationBasic            
            % This function computes the functionality of the Communicationn System
            % as the percentage of communication tower in the
            % communication network.
            %================================================================================
            open = 0;
            
            %=== functionality = the percentage of central office and communication tower in the "open" status
            CO = Comm{1};
            CommTower = Comm{2};
            n = length(CO)+length(CommTower);
            for ii = 1:length(CO)
                if strcmp(CO{ii}.Status, 'Open')
                    open = open + 1;
                end
            end
            for ii = 1:length(CommTower)
                if strcmp(CommTower{ii}.Status, 'Open')
                    open = open + 1;
                end
            end
            Function = open/n;
        end
        
        
        %% computer system resilience based on functionality recovery data
        function [Resilience] = ComputeSystemResilience(ActiveSystem, nQmetric, ResilienceMetricChoice, time_horizon, FunctionalitySystem)
            %=======================================================================================
            % function ComputeSystemResilience
            % This function computes the system resilince based on the system recovery data.  
            % ResilienceMetric_Power = ResilienceMetricChoice{1}; 
            % ResilienceMetric_Comm = ResilienceMetricChoice{2}; 
            % ResilienceMetric_Trans = ResilienceMetricChoice{3}; 
            %=======================================================================================
            
            %=== Flag variables to indicate whether a system is turned ON(1)/OFF(0)
            active_power = ActiveSystem(1);
            active_comm = ActiveSystem(2);
            active_trans = ActiveSystem(3);
            
            isys = 1; Functionality_Power = FunctionalitySystem{isys};
            isys = 2; Functionality_Comm = FunctionalitySystem{isys};
            isys = 3; Functionality_Trans = FunctionalitySystem{isys};

            %=== Power System
            isys = 1;
            if active_power
                ResilienceMetricVector = ResilienceMetricChoice{isys}; 
                for ii = 1:nQmetric(isys)
                    Q = Functionality_Power{ii};
                    Resilience_System{isys}{ii} = Library.ComputeResilience(ResilienceMetricVector, Q, time_horizon); 
                end
            else
                Resilience_System{isys} = {[],[],[],[],[],[],[],[],[],[]};
            end
            
            %=== Communication System
            isys = 2;
            if active_comm
                ResilienceMetricVector = ResilienceMetricChoice{isys}; 
                for ii = 1:nQmetric(isys)
                    Q = Functionality_Power{ii};
                    Resilience_System{isys}{ii} = Library.ComputeResilience(ResilienceMetricVector, Q, time_horizon); 
                end
            else
                Resilience_System{isys} = {[],[],[],[],[],[],[],[],[],[]};
            end
 
            %=== Transportation System
            isys = 3;
            if active_trans
                ResilienceMetricVector = ResilienceMetricChoice{isys}; 
                for ii = 1:nQmetric(isys)
                    Q = Functionality_Power{ii};
                    Resilience_System{isys}{ii} = Library.ComputeResilience(ResilienceMetricVector, Q, time_horizon); 
                end 
            else
                Resilience_System{isys} = {[],[],[],[],[],[],[],[],[],[]};
            end
            
            %=== Output 
            Resilience = Resilience_System;
            
        end
        
        
        %% compute the values of resilience metrics based on the functionality recovery data for a system / an object
        function [Resilience] = ComputeResilience(ResilienceMetricVector, Q, time_horizon)  
            %=============================================================================
            % function ComputeResilience
            %=== INPUT Parameters
            % 1. ResilienceMetricVector = a vector of RessilienceMetricChoice(s)
            % (1) RessilienceMetricChoice = 1 as resilience index (Reed et al. 2009)
            % (2) RessilienceMetricChoice = 2 as resilience loss (Sun et al. 2018)
            % (3) RessilienceMetricChoice = 3 as rapidity (Sun et al. 2018)
            % (4) RessilienceMetricChoice = 4 as restoration time, from the time
            % when event strikes until the time step when the system functionality
            % recovers to the greatest functionality at the end of the time horizon.  
            %=== References about the resilience metrics
            %--- Reed et al. (2009) "Methodology for Assessing the Resilience of
            % Networked Infrastructure", IEEE Ssystems Journal, 3(2), 174-180. 
            % https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4912342
            %--- Sun et al. (2018) "Resilience metrics and measurement methods
            % for transportation infrastructure: the state of the art". 
            % Sustainable and Resilient Infrastructure, DOI: 10.1080/23789689.2018.1448663.
            %=============================================================================
             
            % time horizon 
            th = time_horizon; 
            
            % Nmetric: The total number of resilience metrics selected in the input file
            Nmetric = length(ResilienceMetricVector);

            % total number of functionality recovery samples= number of
            % rows in the functionalit recovery data (Q)
            N = size(Q,1); 
            
            for ii = 1:Nmetric
                index = ResilienceMetricVector(ii); % index of the resilience metric in the vector of ResilienceMetrics
                switch index
                    case 1  % resilience index               
                        Resilience(ii,:) = inv(th)*sum(Q,2)';
                    case 2 % resilience loss
                        for j = 1:N
                            m1 = max(Q(j,2:end),[],2);
                            % tr: time when the functionality recovers to 100% 
                            if m1 == 1
                                tr = min(find(Q(j,2:end) == m1));    
                                Resilience(ii,j) = tr - sum(Q(j,1:tr));
                            else % functionality is not fully recovered, tr is unkown, Resilience loss = -1;  
                                Resilience(ii,j) = -1;
                            end
                        end
                    case 3 % rapidity
                        for j = 1:N
                            m1 = max(Q(j,2:end),[],2);
                            m2 = min(Q(j,2:end),[],2);
                            it1 = min(find(Q(j,2:end) == m1));
                            it2 = min(find(Q(j,2:end) == m2));
                            Resilience(ii,j) = inv(it1-it2)*(m1-m2);
                        end
                    case 4 % restoration time
                        for j = 1:N
                            m1 = max(Q(j,2:end),[],2);
                            m2 = min(Q(j,2:end),[],2);
                            it1 = min(find(Q(j,2:end) == m1));
                            it2 = min(find(Q(j,2:end) == m2));
                            Resilience(ii,j) = it1-it2; 
                        end
                        
                end
            end
            
            
        end
 
        %% compute statistics of system functionality and resilience 
        function [StatisticsFunctionality, StatisticsResilience] = ComputeStatistics(ActiveSystem, nQmetric, QuSys, Resilience)
            %==========================================================================
            % function computeStatistics
            % This function computes the statistics of functionality and resilience.   
            %==========================================================================
            
            %=== Flag variables to indicate whether a system is turned ON(1)/OFF(0)
            active_power = ActiveSystem(1);
            active_comm = ActiveSystem(2);
            active_trans = ActiveSystem(3);
            
            %=== Power System
            isys = 1;
            if active_power
                for imtr = 1:nQmetric(isys)
                    Q = QuSys{isys}{imtr};
                    R = Resilience{isys}{imtr};
                    StatisticsFunctionality{isys}{imtr} = [mean(Q); std(Q); min(Q); prctile(Q,25); prctile(Q,50); prctile(Q,75); max(Q)];
                    StatisticsResilience{isys}{imtr} = [mean(R,2), std(R,[],2), min(R,[],2), prctile(R,25,2),  prctile(R,50,2),  prctile(R,75,2), max(R,[],2)];
                end
            end
            
            %=== Communication System
            isys = 2;
            if active_comm
                for imtr = 1:nQmetric(isys)
                    Q = QuSys{isys}{imtr};
                    R = Resilience{isys}{imtr};
                    StatisticsFunctionality{isys}{imtr} = [mean(Q); std(Q); min(Q); prctile(Q,25); prctile(Q,50); prctile(Q,75); max(Q)];
                    StatisticsResilience{isys}{imtr} = [mean(R,2), std(R,[],2), min(R,[],2), prctile(R,25,2),  prctile(R,50,2),  prctile(R,75,2), max(R,[],2)];
                end
            end            
                    
            %=== Transportation System
            isys = 3;
            if active_trans
                for imtr = 1:nQmetric(isys)
                    Q = QuSys{isys}{imtr};
                    R = Resilience{isys}{imtr};
                    StatisticsFunctionality{isys}{imtr} = [mean(Q); std(Q); min(Q); prctile(Q,25); prctile(Q,50); prctile(Q,75); max(Q)];
                    StatisticsResilience{isys}{imtr} = [mean(R,2), std(R,[],2), min(R,[],2), prctile(R,25,2),  prctile(R,50,2),  prctile(R,75,2), max(R,[],2)];
                end
            end   
            
        end        
        
%         %% compute statistics of system functionality and resilience 
%         function [StatisticsFunctionality, StatisticsResilience] = ComputeStatistics(FunctionalityPower, FunctionalityCommunication, FunctionalityTransportation, Resilience)
%             %==========================================================================
%             % function computeStatistics
%             % This function computes the statistics of functionality and resilience.   
%             %==========================================================================
%         
%             Q{1} = FunctionalityPower;              
%             Q{2} = FunctionalityCommunication;              
%             Q{3} = FunctionalityTransportation;              
%             for ii = 1:length(Q)
%                 StatisticsFunctionality{ii} = [mean(Q{ii},2); std(Q{ii},[],2); min(Q{ii},[],2); prctile(Q{ii},25,2); prctile(Q{ii},50,2); prctile(Q{ii},75,2); max(Q{ii},[],2)];
%             end
% 
%             for ii = 1:length(Q)
%                 StatisticsResilience{ii} = [mean(Resilience{ii},2), std(Resilience{ii},[],2), min(Resilience{ii},[],2), prctile(Resilience{ii},25,2),  prctile(Resilience{ii},50,2),  prctile(Resilience{ii},75,2), max(Resilience{ii},[],2)];
%             end
%             
%         end

        %% Plot and Data Managment
        % Save the Data log to text file
        function SaveLogDamage(num, PowerSet, CommSet, TransSet, name)
            %===================================================================
            % function SaveLogDamage
            % This function saves the initial damage scenario for every
            % system at every damage scenario sample.
            %===================================================================
            % Field
            Branch = PowerSet{1};
            Bus = PowerSet{2};
            Generator = PowerSet{3};
            TransmissionTower = PowerSet{4};
            
            Centraloffice = CommSet{1};
            CommunicationTower = CommSet{2};
            Cellline = CommSet{3};            
 
            Road = TransSet{1};
            Bridge = TransSet{2};
            TrafficLight = TransSet{3};
            
            fileID = fopen(strcat( deblank(name), '/txt/Data_DamageScenarioSample.txt'),'a');
            
            if num > 1
                fprintf(fileID,'\n');
            end
            
            fprintf(fileID,'Sample %3d:\n',num);
            
            %=== Power System
            fprintf(fileID,'Power System:\n');
            
            %--- power plant (generator)
            fprintf(fileID,'Generator:\n');
            fprintf(fileID,'Number  Location  Type  Status Functionality \n');
                        
            for i = 1:length(Generator)
                fprintf(fileID,'%-6d %-8s %-20s %-7s %-10.4f \n', ...
                    Generator{i}.Number, strcat('[', num2str(Generator{i}.Location), ']'),...
                    Generator{i}.Type, Generator{i}.Status, Generator{i}.Functionality);
            end
            
            %--- electric substation (bus)
            fprintf(fileID,'Bus:\n');
            fprintf(fileID,'Number Location Generator BusType Status Functionality \n');
            
            for i = 1:length(Bus)
                if isempty(Bus{i}.Generator)
                    gen = ' ';
                else
                    gen = num2str(Bus{i}.Generator);
                end
                fprintf(fileID,'%-6d %-8s %-9s %-20s %-7s %-10.4f \n', Bus{i}.Number, ...
                    strcat('[', num2str(Bus{i}.Location), ']'), gen, ...
                    Bus{i}.Type, Bus{i}.Status, Bus{i}.Functionality);
            end
            
            %--- transmission tower
            fprintf(fileID,'TransmissionTower:\n');
            fprintf(fileID,'Number Location TowerType Status Functionality \n');
            
            for i = 1:length(TransmissionTower)
                fprintf(fileID,'%-6d %-8s %-20s %-7s %-10.4f\n', ...
                    TransmissionTower{i}.Number, strcat('[', num2str(TransmissionTower{i}.Location), ']'), ...
                    TransmissionTower{i}.Type, TransmissionTower{i}.Status, TransmissionTower{i}.Functionality);
            end
            
            %--- power line (branch)
            fprintf(fileID,'Branch:\n');
            fprintf(fileID,'Number Start_Location End_Location ConnectedObj1  ConnectedObj2 BranchType Status Functionality \n');
            
            for i = 1:length(Branch)
                fprintf(fileID,'%-6d %-14s %-12s %-30s %-30s %-20s %-7s %-10.4f \n', ...
                    Branch{i}.Number, strcat('[', num2str(Branch{i}.Start_Location), ']'), ...
                    strcat('[', num2str(Branch{i}.End_Location), ']'), ...
                    Branch{i}.connectedObj1, Branch{i}.connectedObj2, ...
                    Branch{i}.Type, Branch{i}.Status, Branch{i}.Functionality);
            end
            
            
            %=== Communication System
            fprintf(fileID,'Communication System:\n');
            
            %--- central office
            fprintf(fileID,'Centraloffice:\n');
            fprintf(fileID,'Number Location  BldgType Status Functionality \n');
            
            for i = 1:length(Centraloffice)
                fprintf(fileID,'%-6d %-8s %-20s %-7s %-10.4f \n', Centraloffice{i}.Number,...
                    strcat('[', num2str(Centraloffice{i}.Location), ']'), ...
                    Centraloffice{i}.Type, Centraloffice{i}.Status, Centraloffice{i}.Functionality);
            end
            
            %--- communication tower
            fprintf(fileID,'CommunicationTower:\n');
            fprintf(fileID,'Number Location  TowerType Status Functionality \n');
            
            for i = 1:length(CommunicationTower)
                fprintf(fileID,'%-6d %-8s %-20s %-7s %-10.4d\n', CommunicationTower{i}.Number,...
                    strcat('[', num2str(CommunicationTower{i}.Location), ']'), ...
                    CommunicationTower{i}.Type, CommunicationTower{i}.Status, CommunicationTower{i}.Functionality);
            end

            %--- communication line (Cellline)
            fprintf(fileID,'Cellline:\n');
            fprintf(fileID,'Number Start_Location End_Location ConnectedObj1  ConnectedObj2 Type Status Functionality \n');
            
            for i = 1:length(Cellline)
                fprintf(fileID,'%-6d %-14s %-12s %-30s %-30s %-20s %-7s %-10.4f \n', Cellline{i}.Number, ...
                    strcat('[', num2str(Cellline{i}.Start_Location), ']'), ...
                    strcat('[', num2str(Cellline{i}.End_Location), ']'), ...
                    Cellline{i}.connectedObj1, Cellline{i}.connectedObj2, ...
                    Cellline{i}.Type, Cellline{i}.Status, Cellline{i}.Functionality);
            end
            

            %=== Transportation System
            fprintf(fileID,'Transportation System:\n');
            
            %--- road
            fprintf(fileID,'Road:\n');
            fprintf(fileID,'Number Start_Location End_Location Start_node End_node Type Status Functionality\n');
            
            for i = 1:length(Road)
                fprintf(fileID,'%-6d %-14s %-12s %-5d %-5d %-20s %-15s %-10.4f \n', Road{i}.Number, ...
                    strcat('[', num2str(Road{i}.Start_Location), ']'), ...
                    strcat('[', num2str(Road{i}.End_Location), ']'), ...
                    Road{i}.Start_Node, Road{i}.End_Node, ...
                    Road{i}.Type, Road{i}.Status, Road{i}.Functionality);
            end

            %--- bridge
            fprintf(fileID,'Bridge:\n');
            fprintf(fileID,'Number Location  Type Status Functionality\n');
            
            for i = 1:length(Bridge)
                fprintf(fileID,'%-6d %-8s %-20s %-7s %-10.4f \n', Bridge{i}.Number,...
                    strcat('[', num2str(Bridge{i}.Location), ']'), ...
                    Bridge{i}.Type, Bridge{i}.Status, Bridge{i}.Functionality);
            end
            
            %--- traffic light
            fprintf(fileID,'TrafficLight:\n');
            fprintf(fileID,'Number Location  Type Status Functionality\n');
            if size(TrafficLight,1)~=0
                if size(TrafficLight{1},1)>1
                    for i = 1:length(TrafficLight)
                        fprintf(fileID,'%-6d %-8s %-20s %-7s %-10.4f \n', TrafficLight{i}.Number,...
                            strcat('[', num2str(TrafficLight{i}.Location), ']'), ...
                            TrafficLight{i}.Type, TrafficLight{i}.Status, TrafficLight{i}.Functionality);
                    end
                end
            end
            
            fclose(fileID);
            
        end
        
        %% Save the Schedule log to text file
        function SaveLogSchedule(Ndmgsample, total_schedule, total_date, fname)
            %=====================================================================
            % function SaveScheduleLog
            % save the schedule scample into a txt file named "Data_Schedule.txt"
            %==== Input Variables 
            % 1.Ndmgsample: No. of damage samples
            % 2.total_schedule: the variable of Total_Schedule output from
            % function Library.RestorationPlan at Step 3. This variable
            % documents the restoration plan developed from either the priority
            % ranking of invidial system (Scheme 3A), or based on the optimal 
            % solution of an optimization formation (Scheme 3B). 
            % 3.total_date: the variable of Total_Date output from
            % function Library.RestorationPlan at Step 3. This variable
            % shows the restoration completion date based on the schedule.
            % While using Schem 3A, total_date =[]; while using scheme 3B,
            % total_date is the restoration completion time from the
            % optimization schedule. 
            % 4.fname: the folder name for storing output files.   
            %=====================================================================

            fileID = fopen(strcat( deblank(fname), '/txt/Data_Schedule.txt'),'a');
            for ii = 1:Ndmgsample

                    schedule = total_schedule{ii};
                    
                    pow = schedule{1};
                    comm = schedule{2};
                    trans = schedule{3};
                    
                    if ii > 1 
                        fprintf(fileID,'\n');
                    end
                    
                    fprintf(fileID,'==== Damage Sample %4d:\n', ii);
                    
                    %--- Power system 
                    fprintf(fileID,'---- Power System:\n');
                    %fprintf(fileID,'%-20s %-5s %-5s\n','Compoment', 'Number', 'Date');
                    fprintf(fileID,'%-20s %-10s %-10s %-10s %-10s\n','Compoment', 'Number', 'TaskNumber', 'StartTime', 'EndTime');
                    
                    for k = 1:length(pow)
                        tem = strsplit(pow{k},'/');
                        date = 0;
                        if ~isempty(total_date{ii}{1}) % Power: isys = 1
                            %date = total_date{ii}{1}(k);
                            date = total_date{ii}{1}(k,:);
                        end
                        % fprintf(fileID,'%-20s %-5s %-5s\n',char(tem(1)), char(tem(2)), num2str(date));
                        fprintf(fileID,'%-20s %-10s %-10s %-20s \n',char(tem(1)), char(tem(2)), char(tem(3)), num2str(date));
                    end
                    
                    %--- Communication system
                    fprintf(fileID,'---- Communication System:\n');  
                    %fprintf(fileID,'%-20s %-5s %-5s\n','Compoment', 'Number', 'Date');
                    fprintf(fileID,'%-20s %-10s %-10s %-10s %-10s\n','Compoment', 'Number', 'TaskNumber', 'StartTime', 'EndTime');
                    
                    for k = 1:length(comm)
                        tem = strsplit(comm{k},'/');
                        date = 0;
                        if ~isempty(total_date{ii}{2}) % Comm: isys = 2
                            %date = total_date{ii}{2}(k);
                            date = total_date{ii}{2}(k,:);
                        end
                        %fprintf(fileID,'%-20s %-5s %-5s\n',char(tem(1)), char(tem(2)), num2str(date));
                        fprintf(fileID,'%-20s %-10s %-10s %-20s \n',char(tem(1)), char(tem(2)), char(tem(3)), num2str(date));
                    end
                    
                    %--- Transportation system
                    fprintf(fileID,'---- Transportation System:\n');
                    %fprintf(fileID,'%-20s %-5s %-5s\n','Compoment', 'Number', 'Date');
                    fprintf(fileID,'%-20s %-10s %-10s %-10s %-10s\n','Compoment', 'Number', 'TaskNumber', 'StartTime', 'EndTime');
                    
                    for k = 1:length(trans)
                        tem = strsplit(trans{k},'/');
                        date = 0;
                        if ~isempty(total_date{ii}{3}) % Trans: isys = 3
                            %date = total_date{ii}{3}(k);
                            date = total_date{ii}{3}(k,:);
                        end
                        %fprintf(fileID,'%-20s %-5s %-5s\n',char(tem(1)), char(tem(2)), num2str(date));
                        fprintf(fileID,'%-20s %-10s %-10s %-20s \n',char(tem(1)), char(tem(2)), char(tem(3)), num2str(date));
                    end
            end
            
            fclose(fileID);
        end
        
        %% Save the Functionality log to text file
        function SaveLogFunctionality(FunctionalitySystem, name, Nsample, NRun, qmetric) 
            %%%%%%%%
            %=== Input Parameters
            % qmetric = [functionality metric index of Qpower, functionality metric index of Qcomm, functionality metric index of Qtrans];
            % Nsample: number of random samples
            % th: time_horizon
            
            %--- imtrp, imtrc, imtrt: index of functionality metric(s) for power, communication, and transportation systems, respectively  
            imtrp = qmetric{1};
            imtrc = qmetric{2};
            imtrt = qmetric{3};
            isys = 1; Pow = FunctionalitySystem{isys};
            isys = 2; Comm = FunctionalitySystem{isys};
            isys = 3; Trans = FunctionalitySystem{isys};
            fileID = fopen(strcat( deblank(name), '/txt/Data_Functionality.txt'),'a');           
            
            [~, th] = size(Pow{1});
            
            for isample = 1:Nsample
                
                %ksample = ceil(jsample/NRun);
                kdmg = ceil(isample/NRun);
                jdur = isample - (kdmg-1)*NRun;
                
                fprintf(fileID, '\n');
                fprintf(fileID, 'Sample %4d (meaning DamageSample%4d & DurationSample %4d) :\n', isample, kdmg, jdur);
                fprintf(fileID, '\n');
                    
                fprintf(fileID,'%-17s', 'Time (day): ');
                for k = 1:th
                    fprintf(fileID,'%-7s ', num2str(k));
                end
                fprintf(fileID,'\n');
                    
                for j = 1:length(imtrp)
                    jj = imtrp(j);
                    str = strcat('Power-FunctionalityMetric#', num2str(jj),': ');
                    fprintf(fileID,'%-40s', str);
                    for k = 1:th
                        fprintf(fileID,'%-7s ', num2str(Pow{jj}(isample, k)));
                    end
                    fprintf(fileID,'\n');
                end               
                
                for j = 1:length(imtrc)
                    jj = imtrc(j);
                    str = strcat('Communication-FunctionalityMetric#', num2str(jj),': ');
                    fprintf(fileID,'%-40s', str);
                    for k = 1:th
                        fprintf(fileID,'%-7s ', num2str(Comm{jj}(isample, k)));
                    end
                    fprintf(fileID,'\n');
                end

                for j = 1:length(imtrt)
                    jj = imtrt(j);
                    str = strcat('Transportation-FunctionalityMetric#', num2str(jj),': ');
                    fprintf(fileID,'%-40s', str);
                    for k = 1:th
                        fprintf(fileID,'%-7s ', num2str(Trans{jj}(isample, k)));
                    end
                    fprintf(fileID,'\n');
                end
                
            end
            fclose(fileID);
        end
        
        %% Plot figures of system functionality samples and means
        function PlotFigureFunctionality(Nsamples, NRun, time_horizon, FunctionalitySystem, fname, figoption)
            %===============================================================
            % function PlotFigureFunctionality
            % plot figures of system functionality recovery samples and
            % means over time. 
            %===============================================================
            % Line Width
            lw = 2; 
            %==== Field: initial set up
            isys = 1; Functionality_Power = FunctionalitySystem{isys};  % active_power = ActiveSystem(isys);
            isys = 2; Functionality_Communication = FunctionalitySystem{isys};  % active_comm = ActiveSystem(isys);
            isys = 3; Functionality_Transportation = FunctionalitySystem{isys}; % active_trans = ActiveSystem(isys);

            %==== select the first functionality metric (Functionality_Basic)
            %imtr = 1;            
            imtr = 2;
            
            %==== plot figure 
            fig = figure('visible',figoption, 'Renderer', 'painters', 'Position', [100 100 1200 600]);
            ncol = 3;
            nrow = 1;
            
            Nsample = Nsamples*NRun;
            if Nsample>=2
                nrow = 2;
            end
            
            %==== Power

            subplot(nrow,ncol,1),
            for ii=1:Nsample
                hold on, stairs(Functionality_Power{imtr}(ii,:)); 
            end
            ylim([0,1]);
            box on; grid on; 
            xlabel('Time (day)');
            ylabel('Q_{power}');
            title('Power Functionality Sample');

            %==== Communication
            subplot(nrow,ncol,2),
            for ii=1:Nsample
                hold on, stairs(Functionality_Communication{imtr}(ii,:));
            end
            ylim([0,1]);
            box on; grid on; 
            xlabel('Time (day)');
            ylabel('Q_{comm}');
            title('Communication Functionality Sample');

            %==== Transportation 
            subplot(nrow,ncol,3),
            for ii=1:Nsample
                hold on, stairs(Functionality_Transportation{imtr}(ii,:));
            end
            ylim([0,1]);
            box on; grid on; 
            xlabel('Time (day)');
            ylabel('Q_{trans}');
            title('Transportation Functionality Sample');
            
            
            % Plot Qmean if there are multiple samples 
            if Nsample>=2
                subplot(nrow,ncol,4);
                mf=mean(Functionality_Power{imtr});
                stairs(mf(1,:), 'LineWidth', lw);
                ylim([0,1]);
                box on; grid on; 
                xlabel('Time (day)');
                ylabel('Functionality');
                title('Power Functionality Mean');
                
                subplot(nrow,ncol,5);
                mf=mean(Functionality_Communication{imtr});
                stairs(mf(1,:), 'LineWidth', lw);
                ylim([0,1]);
                box on; grid on; 
                xlabel('Time (day)');
                ylabel('Functionality');
                title('Communication Functionality Mean');
                
                subplot(nrow,ncol,6);
                mf=mean(Functionality_Transportation{imtr});
                stairs(mf(1,:), 'LineWidth', lw);
                ylim([0,1]);
                box on; grid on; 
                xlabel('Time (day)');
                ylabel('Functionality');
                title('Transportation Functionality Mean');       
                
            end
            
            saveas(fig,strcat( deblank(fname), '/plot/Qsystem.jpg'));

        end

        %% Reset data for three systems and the Dictionary, as well as Neighborhood.
        function [Power, Commu, Trans, Dictionary, Neighborhood] = ResetData(Total)
            % ================================================================
            % function ResetData
            % Reset initial infrastructure data to cells of Power, Comm, Trans, Dictionary,Neighborhood
            % ================================================================
            Power = Total{1};
            Commu = Total{2};
            Trans = Total{3};
            Dictionary = Total{4};
            Neighborhood = Total{5};
        end
        
        %% Clean old .txt .jpg .mat files
        function CleanOldData(name)
            
            ans = exist(deblank(name),'dir');
            if ans == 7
                rmdir(deblank(name), 's');
            end
        end
        
        %% Create Folder
        function CreateFolder(name)
            
            x = exist(deblank(name),'dir');
            if x ~= 7
                mkdir(deblank(name));
            end
            
            filename = strcat( deblank(name),'/mat');
            x = exist(filename,'dir');
            if x ~= 7
                mkdir(filename);
            end
            
            filename = strcat( deblank(name),'/task');
            x = exist(filename,'dir');
            if x ~= 7
                mkdir(filename);
            end
            
            filename = strcat( deblank(name),'/schedule');
            x = exist(filename,'dir');
            if x ~= 7
                mkdir(filename);
            end
            
            filename = strcat( deblank(name),'/txt');
            x = exist(filename,'dir');
            if x ~= 7
                mkdir(filename);
            end
            
            filename = strcat( deblank(name),'/plot');
            x = exist(filename,'dir');
            if x ~= 7
                mkdir(filename);
            end
        end
        
        %% Read Input Files in the csv format from the sub-folder named "Input"
        
        function [power, trans, comm, Neighborhood_Set] = readInput(filename, pow_check, comm_check, trans_check, Pow, Trans, Comm, Dictionary, Neighborhood_Set,active_traLights)
            %=================================================================
            % function readInput
            % 1. Read Input Files in csv format 
            % 2. Prepare intital variables and load input data into them 
            %=================================================================
            tmp = strsplit(filename,'.');
            name = tmp{1};
            fpath = 'Input';

            myfolder = strcat(pwd,'/',fpath);
            mydir = fullfile(myfolder,filename);
            if ~strcmp('TrafficLights', name)
                table = readtable([mydir,'.csv']); warning off;
            elseif active_traLights
                table = readtable([mydir,'.csv']); warning off;
            end
            
            %==== Field
            % Assign initial values = [] for all variables of ObjectSet
            Branch_Set = Pow{1};
            Bus_Set = Pow{2};
            Generator_Set = Pow{3};

            Centraloffice_Set = Comm{1};
            CommunicationTower_Set = Comm{2};
            Cellline_Set = Comm{3};
             
            Road_Set = Trans{1};
            Bridge_Set = Trans{2};
            TrafficLight_Set = Trans{3};
            RoadNode_Set = Trans{4};
            
            Task_Set = {};
            
            centraloffice_check = comm_check{1};
            communicationtower_check = comm_check{2};
            
            road_check = trans_check{1};
            bridge_check = trans_check{2};
            trafficlight_check = trans_check{3};

            index_Branch = length(Branch_Set) + 1;
            index_Generator = length(Generator_Set) + 1;
            index_Bus = length(Bus_Set) + 1;
            index_Centraloffice = length(Centraloffice_Set) + 1;
            index_Communicationtower = length(CommunicationTower_Set) + 1;
            index_Road = length(Road_Set) + 1;
            index_Bridge = length(Bridge_Set) + 1;
            index_TrafficLight = length(TrafficLight_Set) + 1;
            index_RoadNode = length(RoadNode_Set) + 1;
            index_Neighborhood = length(Neighborhood_Set) + 1;
            index_task = 1;
            
            %=== PowerPlant
            if strcmp('PowerPlants', name) || strcmp('PowerPlant', name)
                for i = 1:size(table,1)
                    try
                        if ~isKey(pow_check, table.Name{i}) && ~isempty(table.Name{i})
                            keySet = table.Name{i};
                            valueSet = 1;
                            newMap = containers.Map(keySet,valueSet);
                            pow_check = [pow_check; newMap];
                            
							%=== assign the generation capacity of the power plant								 
                            try
                                if ~isnan(table.GenerationCapacity_MW_(i))
                                    cap = table.GenerationCapacity_MW_(i);
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 2');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
							%=== assign the fuel type of the power plant										 
                            try
                                if ~isempty(table.FuelType{i})
                                    type = table.FuelType{i};
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 3');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
 
							%=== assign longitude										  
                            try
                                if ~isnan(table.Longitude(i))
                                    endlocation = table.Longitude(i);
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
							%=== assign latitude											  
                            try
                                if ~isnan(table.Latitude(i))
                                    startlocation = table.Latitude(i);
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 3');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign the total household (TotalHH) served by this power plant
                            try
                                if ~isnan(table.TotalHH(i))
                                    HouseholdServed = table.TotalHH(i);   
                                else
                                    HouseholdServed = 0;
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign total population (TotalPop) served by this power plant
                            try
                                if ~isnan(table.TotalPop(i))
                                    PopulationServed = table.TotalPop(i);  
                                else
                                    PopulationServed = 0;
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                                                        
                            
                            % disp(index_Generator);
                            %Location = [startlocation, endlocation];
                            % SocioEconomicData = [PopulationServed, HouseholdServed];
                            Generator_Set{index_Generator} = Generator(index_Generator, [startlocation, endlocation]);
                            Generator_Set{index_Generator}.PopulationServed = PopulationServed;
                            Generator_Set{index_Generator}.HouseholdServed = HouseholdServed;                          
                            Bus_Set{index_Bus} = Bus(index_Bus, cap, type, keySet, [startlocation, endlocation], PopulationServed, HouseholdServed, Generator_Set{index_Generator}.uniqueID);
                            Generator_Set{index_Generator}.Bus = Bus_Set{index_Bus}.uniqueID;
                            
                            Dictionary(Bus_Set{index_Bus}.uniqueID) = Bus_Set(index_Bus);
                            Dictionary(Generator_Set{index_Generator}.uniqueID) = Generator_Set(index_Generator);
                            
                            
                            index_Bus = index_Bus + 1;
                            index_Generator = index_Generator + 1;
                        else
                            message = strcat('ERROR: Empty Name or Name already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
            end
            
            %=== Bus (Substation)
            if strcmp('Substations', name) || strcmp('Substation', name)
                for i = 1:size(table,1)
                    try
                        if ~isKey(pow_check, table.Name{i}) && ~isempty(table.Name{i})
                            keySet = table.Name{i};
                            valueSet = 1;
                            newMap = containers.Map(keySet,valueSet);
                            pow_check = [pow_check; newMap];
                            
                            %=== assign maximum voltage
                            try
                                if ~isnan(table.MaxVoltage_kV_(i))
                                    %=== assign voltage and type based on the substation type in HAZUS
                                    % cap: kilo-voltage (kV) 
                                    % type: Structural Type
                                    % ESS1 Low voltage (115 kV) with anchored components
                                    % ESS2 Low voltage (115 kV) with unanchored components 
                                    % ESS3 Medium voltage (115~500 kV) with anchored components 
                                    % ESS4 Medium voltage (115~500 kV) with unanchored components 
                                    % ESS5 High voltage (500 kV) with anchored components
                                    % ESS6 High voltage (500 kV) with unanchored components

                                    cap = table.MaxVoltage_kV_(i);
                                    if ge(cap,500)
                                        type = 'ESS5';
                                    elseif ge(cap,115) && lt(cap,500)
                                        type = 'ESS3';
                                    else
                                        type = 'ESS1';
                                    end
                                    
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 2');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign total household (TotalHH) served by this substation
                            try
                                if ~isnan(table.TotalHH(i))
                                    HouseholdServed = table.TotalHH(i);   
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign total population (TotalPop) served by this substation
                            try
                                if ~isnan(table.TotalPop(i))
                                    PopulationServed = table.TotalPop(i);  
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign latitude
                            try
                                if ~isnan(table.Latitude(i))
                                    startlocation = table.Latitude(i);
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign longitude
                            try
                                if ~isnan(table.Longitude(i))
                                    endlocation = table.Longitude(i);
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 5');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end

                            Bus_Set{index_Bus} = Bus(index_Bus, cap, type, keySet, [startlocation, endlocation], PopulationServed, HouseholdServed);
                            
                            Dictionary( Bus_Set{index_Bus}.uniqueID) =  Bus_Set(index_Bus);
                            index_Bus = index_Bus + 1;
                        else
                            message = strcat('ERROR: Empty Name or Name already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
            end
            
            %=== Branch
            if strcmp('Power_Line', name) || strcmp('PowerLine', name)
                for i = 1:size(table,1)
                    try
                        if ~isempty(table.BranchNode1{i}) && ~isempty(table.BranchNode2{i}) && ~isnan(table.Voltage_kV_(i)) && ~isempty(table.TowerType{i})
                            name1 = table.BranchNode1{i};
                            name2 = table.BranchNode2{i};
                            cap = table.Voltage_kV_(i);
                            type = table.TowerType{i};
                            index1 = 0;
                            index2 = 0;
                            
                            for j = 1:length(Bus_Set)
                                if strcmp(Bus_Set{j}.Name, name1)
                                    index1 = j;
                                end
                                
                                if strcmp(Bus_Set{j}.Name, name2)
                                    index2 = j;
                                end
                            end
                            
                            try
                                start_location = Bus_Set{index1}.Location;
                                end_location = Bus_Set{index2}.Location;
                                
                                if strcmp(Bus_Set{index1}.Type, 'Gas') || strcmp(Bus_Set{index1}.Type, 'Nuclear') || strcmp(Bus_Set{index1}.Type, 'Coal') || strcmp(Bus_Set{index2}.Type, 'Gas') || strcmp(Bus_Set{index2}.Type, 'Nuclear') || strcmp(Bus_Set{index2}.Type, 'Coal')
                                    proprity = 1;
                                else
                                    proprity = 2;
                                end
                                
                                Branch_Set{index_Branch} = Branch(index_Branch, start_location, end_location, cap, type, proprity, Bus_Set{index1}.uniqueID,  Bus_Set{index2}.uniqueID);
                                %Bus_Set{index1}.Branch = [Bus_Set{index1}.Branch, index_Branch];
                                %Bus_Set{index2}.Branch = [Bus_Set{index2}.Branch, index_Branch];                                
                                %Dictionary(Branch_Set{index_Branch}.uniqueID) = Branch_Set(index_Branch);
                                index_Branch = index_Branch + 1;
                            catch exception
                                C = {'Missing data!!', name1, 'or', name2, 'does not exist!'};
                                str = strjoin(C);
                                error(str);
                            end
                        else
                            message = strcat('ERROR: Empty Name or Name already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
            end
            
            
            %=== Central Office
            if strcmp('CentralOffices', name)
                for i = 1:height(table)
                    try
                        if  ~isempty(table(i,2)) && ~isempty(table(i,3)) && ~isempty(table(i,4)) && ~isempty(table(i,5)) && ~isKey(centraloffice_check, cell2mat(table2cell(table(i,4))))
                            keySet = table.Label{i};
                            valueSet = 1;
                            newMap = containers.Map(keySet,valueSet);
                            centraloffice_check = [centraloffice_check; newMap];
                            %=== assign total household (TotalHH) served by this central offfice 
                            try
                                if ~isnan(table.TotalHH(i))
                                    HouseholdServed = table.TotalHH(i);   
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign total population (TotalPop) served by this central offfice 
                            try
                                if ~isnan(table.TotalPop(i))
                                    PopulationServed = table.TotalPop(i);  
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            Centraloffice_Set{index_Centraloffice} = Centraloffice(index_Centraloffice, [cell2mat(table2cell(table(i,2))), cell2mat(table2cell(table(i,3)))]);
                            Centraloffice_Set{index_Centraloffice}.Company = cell2mat(table2cell(table(i,5)));
                            Centraloffice_Set{index_Centraloffice}.Label = cell2mat(table2cell(table(i,4)));
                            Centraloffice_Set{index_Centraloffice}.PopulationServed = PopulationServed;
                            Centraloffice_Set{index_Centraloffice}.HouseholdServed = HouseholdServed; 			  
                            
                            Dictionary(Centraloffice_Set{index_Centraloffice}.uniqueID) = Centraloffice_Set(index_Centraloffice);
                            index_Centraloffice = index_Centraloffice + 1;
                            
                        else
                            message = strcat('ERROR: Empty Value or Value already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
            end
            
            %=== Communication Tower
            if strcmp('CommunicationTowers', name)
                for i = 1:height(table)
                    try
                        if  ~isempty(table(i,1)) && ~isempty(table(i,2)) %&& ~isempty(table(i,3)) && ~isempty(table(i,4)) && ~isempty(table(i,6)) && ~isempty(table(i,9)) && ~isempty(table(i,12)) && ~isempty(table(i,15)) && ~isempty(table(i,16)) &&~isempty(table(i,17)) && ~isempty(table(i,18)) && ~isKey(communicationtower_check, cell2mat(table2cell(table(i,1))))
                            keySet = cell2mat(table2cell(table(i,1)));
                            valueSet = 1;
                            newMap = containers.Map(keySet,valueSet);
                            communicationtower_check = [communicationtower_check; newMap];
                            
                            %=== assign total household (TotalHH) served by this communication tower 
                            try
                                if ~isnan(table.TotalHH(i))
                                    HouseholdServed = table.TotalHH(i);   
                                else
                                    HouseholdServed = 0;
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign total population (TotalPop) served by the communication tower 
                            try
                                if ~isnan(table.TotalPop(i))
                                    PopulationServed = table.TotalPop(i); 
                                else
                                    PopulationServed = 0;    
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            %=== assign height by the communication tower 
                            try
                                if ~isnan(table.Height_m_(i))
                                    Height = table.Height_m_(i); 
                                else
                                    Height = 0;    
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end    
                            CommunicationTower_Set{index_Communicationtower} = CommunicationTower(index_Communicationtower, [cell2mat(table2cell(table(i,16))), cell2mat(table2cell(table(i,17)))]);
                            CommunicationTower_Set{index_Communicationtower}.CommunicationTowerID = cell2mat(table2cell(table(i,1)));
                            %CommunicationTower_Set{index_Communicationtower}.Tract = sprintf('%d',cell2mat(table2cell(table(1,3))));
                            %CommunicationTower_Set{index_Communicationtower}.Comment = cell2mat(table2cell(table(i,18)));
                            %CommunicationTower_Set{index_Communicationtower}.BackupPower = cell2mat(table2cell(table(i,15)));
                            %CommunicationTower_Set{index_Communicationtower}.Usage = cell2mat(table2cell(table(i,12)));
                            %CommunicationTower_Set{index_Communicationtower}.Owner = cell2mat(table2cell(table(i,9)));
                            %CommunicationTower_Set{index_Communicationtower}.City = cell2mat(table2cell(table(i,6)));
                            %CommunicationTower_Set{index_Communicationtower}.Name = cell2mat(table2cell(table(i,4)));
                            type = cell2mat(table2cell(table(i,2)));
                            if contains(type,'Lattice')
                                type = 'Lattice';
                            end
                            CommunicationTower_Set{index_Communicationtower}.StructuralType = type;
                            CommunicationTower_Set{index_Communicationtower}.Height = Height;
                            CommunicationTower_Set{index_Communicationtower}.PopulationServed = PopulationServed;
                            CommunicationTower_Set{index_Communicationtower}.HouseholdServed = HouseholdServed;
                            
                            Dictionary(CommunicationTower_Set{index_Communicationtower}.uniqueID) = CommunicationTower_Set(index_Communicationtower);
                            index_Communicationtower = index_Communicationtower + 1;
                        else
                            message = strcat('ERROR: Empty Value or Value already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
            end
            
            %=== Bridge
            if strcmp('Bridges', name)
                for i = 1:size(table,1)
                    try
                        if  ~isempty(table(i,2)) && ~isempty(table(i,3)) && ~isempty(table(i,4)) && ~isempty(table(i,11)) && ~isempty(table(i,10)) && ~isempty(table(i,13)) && ~isempty(table(i,12))
                            keySet = char(table{i,18});
                            valueSet = 1;
                            newMap = containers.Map(keySet,valueSet);
                            bridge_check = [bridge_check; newMap];
                            
                            Bridge_Set{index_Bridge} = Bridge(index_Bridge, [table{i,2}, table{i,3}], char(table{i,18}));

                            if isnan(table{i,66})
                                  table{i,66}=0;
                            end
                            Bridge_Set{index_Bridge}.Width = table{i,66};
                            if isnan(table{i,33})
                                error('MainSpans Input Error');
                            end
                            Bridge_Set{index_Bridge}.MainSpans = table{i,33};
                            Bridge_Set{index_Bridge}.Length=table{i,35};
                            if isnan(table{i,34})
                                error('AppSpans Input Error');
                            end
                            Bridge_Set{index_Bridge}.AppSpans = table{i,34};
                            if isnan(table{i,111})
                                table{i,111}=0;
                            end
							Bridge_Set{index_Bridge}.MaxSpanLength= table{i,111};
							if isnan(table{i,127})
                                error('SkewAngle Input Error');
                            end
                            
                            Bridge_Set{index_Bridge}.Name = char(table.Name{i});
                            Bridge_Set{index_Bridge}.Owner = char(table.MAINT_RESP_DESC{i});
                            Bridge_Set{index_Bridge}.NoOfCarryLinkID=table.NoOfCarryLinkID(i);
                            Bridge_Set{index_Bridge}.CarryLinkID=[table{i,9:11}];
                            Bridge_Set{index_Bridge}.NoOfCrossLinkID=table.NoOfCrossLinkID(i);
                            Bridge_Set{index_Bridge}.CrossLinkID=[table{i,13:15}];                           
                            Bridge_Set{index_Bridge}.CTY_CODE=table.CTY_CODE(i);
                              
                            Bridge_Set{index_Bridge}.SkewAngle = table{i,127};
                            Bridge_Set{index_Bridge}.Year = table{i,29};
                            Bridge_Set{index_Bridge}.Traffic = table{i,135};
                            Bridge_Set{index_Bridge}.Cost = table{i,119};
                            Dictionary(Bridge_Set{index_Bridge}.uniqueID) = Bridge_Set(index_Bridge);
                            
                            index_Bridge = index_Bridge + 1;
                        else
                            message = strcat('ERROR: Empty Value or Value already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg)
                        return;
                    end
                end
            end

            %=== Revised content about Sub-Component analysis for Bridges
            if startsWith(name,'SubBridge')
                sub_id = strsplit(name,'_');
                sub_id = sub_id{2};
                bridgeID = strcat('Bridge', sub_id);
                bridge = Dictionary(bridgeID);
                bridge = bridge{1};
                bridge.HasSub = 1;
                
                abutment_index = 1;
                abutmentFoundation_index = 1;
                bearing_index = 1;
                deck_index = 1;
                columnFoundation_index = 1;
                column_index = 1;
                slab_index = 1;
                girder_index = 1;
                
                %Left abutment foundation
                nlaf = table{1,2}; % # of left abutment foundation
                for ii = 1:nlaf
                    bridge.AbutmentFoundSet = [bridge.AbutmentFoundSet,AbutmentFoundation(ii, bridge.Location, 'left')];
                    Dictionary(bridge.AbutmentFoundSet(abutmentFoundation_index).uniqueID) = bridge.AbutmentFoundSet(abutmentFoundation_index);
                    abutmentFoundation_index = abutmentFoundation_index + 1;
                end
              
                %right abutment foundation
                nraf = table{2,2}; % # of right abutment foundation
                for ii = 1:nraf
                    bridge.AbutmentFoundSet = [bridge.AbutmentFoundSet,AbutmentFoundation(ii, bridge.Location, 'right')];
                    Dictionary(bridge.AbutmentFoundSet(abutmentFoundation_index).uniqueID) = bridge.AbutmentFoundSet(abutmentFoundation_index);
                    abutmentFoundation_index = abutmentFoundation_index + 1;                        
                end      
                
                %Left abutment
                nla =  table{3,2}; % # of left abutment
                for ii = 1:nla
                    bridge.AbutmentSet = [bridge.AbutmentSet,Abutment(ii, bridge.Location, 'left')];
                    Dictionary(bridge.AbutmentSet(abutment_index).uniqueID) = bridge.AbutmentSet(abutment_index);
                    abutment_index = abutment_index + 1;  
                end               
                %Right abutment
                nra = table{4,2}; % # of right abutment 
                for ii = 1:nra
                    bridge.AbutmentSet = [bridge.AbutmentSet,Abutment(ii, bridge.Location, 'right')];
                    Dictionary(bridge.AbutmentSet(abutment_index).uniqueID) = bridge.AbutmentSet(abutment_index);
                    abutment_index = abutment_index + 1;  
                end  
                
                %Rocker bearing 
                nrb = table{5,2}; % # of rocker bearing
                for ii = 1:nrb
                    bridge.BearingSet = [bridge.BearingSet,Bearing(ii,bridge.Location,'Rocker')];
                    Dictionary(bridge.BearingSet(bearing_index).uniqueID) = bridge.BearingSet(bearing_index);
                    bearing_index = bearing_index + 1;
                end
                %Fixed bearing
                nfb = table{6,2}; % # of fixed bearing
                for ii = 1:nfb
                    bridge.BearingSet = [bridge.BearingSet,Bearing(ii,bridge.Location,'Fixed')];
                    Dictionary(bridge.BearingSet(bearing_index).uniqueID) = bridge.BearingSet(bearing_index);
                    bearing_index = bearing_index + 1;
                end     
                
                %Deck
                nd = table{7,2}; % # of decks
                for ii = 1:nd
                    bridge.DeckSet = [bridge.DeckSet,Deck(ii,bridge.Location)];
                    Dictionary(bridge.BearingSet(deck_index).uniqueID) = bridge.DeckSet(deck_index);
                    deck_index = deck_index + 1;                       
                end
                
                %approach slab
                nas = table{8,2} > 0 % # of approach slabs
                for ii = 1:nas
                    bridge.SlabSet = [bridge.SlabSet,ApproachSlab(ii, bridge.Location,'NA')];
                    Dictionary(bridge.SlabSet(slab_index).uniqueID) = bridge.SlabSet(slab_index);
                    slab_index = slab_index + 1;
                end
                
                %left column foundation
                nlcf = table{9,2}; % # of left column foundation
                for ii = 1:nlcf
                    bridge.ColumnFoundSet = [bridge.ColumnFoundSet,ColumnFoundation(ii, bridge.Location,'left')];
                    Dictionary(bridge.ColumnFoundSet(columnFoundation_index).uniqueID) = bridge.ColumnFoundSet(columnFoundation_index);
                    columnFoundation_index = columnFoundation_index + 1;
                end
                
                %right column foundation
                nrcf = table{10,2};
                if nrcf > 0 % # of right column foundation
                    for ii = 1:nrcf
                        bridge.ColumnFoundSet = [bridge.ColumnFoundSet,ColumnFoundation(ii, bridge.Location,'right')];
                        Dictionary(bridge.ColumnFoundSet(columnFoundation_index).uniqueID) = bridge.ColumnFoundSet(columnFoundation_index);
                        columnFoundation_index = columnFoundation_index + 1;
                    end
                end
                
                %left column 
                nlc = table{11,2};
                if nlc > 0 % # of left column
                    for ii = 1:nlc
                        bridge.ColumnSet = [bridge.ColumnSet, Column(ii,bridge.Location,'left')];
                        Dictionary(bridge.ColumnSet(column_index).uniqueID) = bridge.ColumnSet(column_index);
                        column_index = column_index + 1;
                    end
                end 
                
                %right column 
                nrc = table{12,2};
                if nrc > 0 % # of right column
                    for ii = 1:nrc
                        bridge.ColumnSet = [bridge.ColumnSet, Column(ii,bridge.Location,'right')];
                        Dictionary(bridge.ColumnSet(column_index).uniqueID) = bridge.ColumnSet(column_index);
                        column_index = column_index + 1;
                    end
                end 

                %girder 
                ng = table{13,2}; % # of girder
                for ii = 1:ng
                    bridge.GirderSet = [bridge.GirderSet, Girder(ii,bridge.Location,'NA')];
                    Dictionary(bridge.GirderSet(girder_index).uniqueID) = bridge.GirderSet(girder_index);
                    girder_index = girder_index + 1;
                end 
            end 
            
            %=== TrafficLights
            if strcmp('TrafficLights', name)
                if active_traLights
                for ii = 1:size(table,1)
                    try
                        if  ~isempty(table{ii,2}) && ~isempty(table{ii,3}) && ~isempty(table{ii,4}) && ~isempty(table.Link(ii))
                            keySet = num2str(table{ii,9});
                            valueSet = 1;
                            newMap = containers.Map(keySet,valueSet);
                            trafficlight_check = [trafficlight_check; newMap];
                            
                            TrafficLight_Set{index_TrafficLight} = TrafficLight(index_TrafficLight, [table{ii,7}, table{ii,8}]);
                            TrafficLight_Set{index_TrafficLight}.MapKey = num2str(table{ii,9});
                            TrafficLight_Set{index_TrafficLight}.MajorStreet = char(table{ii,14});
                            TrafficLight_Set{index_TrafficLight}.MinorStreet = char(table{ii,16});
							TrafficLight_Set{index_TrafficLight}.RoadLinkID = table.Link(ii);
                            Dictionary(TrafficLight_Set{index_TrafficLight}.uniqueID) = TrafficLight_Set(index_TrafficLight);
                            
                            index_TrafficLight = index_TrafficLight + 1;
                        else
                            message = strcat('ERROR: Empty Value or Value already exist in file: ', filename, ' found at line: ', num2str(ii), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
                end
            end
            

            %=== Road Node
            if strcmp('RoadNode', name)
                for i = 1:size(table,1)
                    %tt = i
                    try
                        if  ~isempty(table{i,1})
                            try
                                if ~isempty(table{i,1})
                                    nodeID = table{i,1};
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 1');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            try
                                if ~isnan(table{i,2})
                                    latitude = table{i,2};
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 2');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            try
                                if ~isnan(table{i,3})
                                    longtitude = table{i,3};
                                else
                                    message = strcat('ERROR: Empty value in file: ', filename, ' found at line: ', num2str(i), ' column: 3');
                                    error(message);
                                end
                            catch exception
                                msg = getReport(exception, 'basic');
                                disp(msg);
                                return;
                            end
                            
                            RoadNode_Set{index_RoadNode} = RoadNode(nodeID, [latitude, longtitude], nodeID);
                            Dictionary(RoadNode_Set{index_RoadNode}.uniqueID) = RoadNode_Set(index_RoadNode);
                            
                            index_RoadNode = index_RoadNode + 1;
                        else
                            message = strcat('ERROR: Empty Name or Name already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
            end
            
            %=== Road 
            if strcmp('RoadLink', name)
                for i = 1:size(table,1)
                    try
                        if  ~isempty(table{i,1})                     
                            keySet = num2str(table{i,1});
                            valueSet = 1;
                            newMap = containers.Map(keySet,valueSet);
                            road_check = [road_check; newMap];
                            
                            Road_Set{index_Road} = Road(index_Road,table{i,2}, table{i,4}, table{i,3}, table{i,5});                      
                            Road_Set{index_Road}.AADT = table{i,6};                
                            Road_Set{index_Road}.numLanes = table{i,7};
                            Road_Set{index_Road}.Length = table{i,8};
                            Road_Set{index_Road}.Speedlimit = table{i,9};
                            %Road_Set{index_Road}.Type = table{i,10};
                            Dictionary(Road_Set{index_Road}.uniqueID) = Road_Set(index_Road);
                            
                            load newroadlinks.mat
                            Road_Set{index_Road}.Bridge_Carr  = bridgeIDcarry{i}';
                            Road_Set{index_Road}.Bridge_Cross = bridgeIDcross{i}';
                            if active_traLights
                            Road_Set{index_Road}.TrafficLight = trafficlightID{i}';
                            end

                            index_Road = index_Road + 1;
                        else
                            message = strcat('ERROR: Empty Value or Value already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                            error(message);
                        end
                    catch exception
                        msg = getReport(exception, 'basic');
                        disp(msg);
                        return;
                    end
                end
            end
       
            %=== Neighborhood
            if strcmp('Neighborhood', name)
                    for i = 1:height(table)
                        try
                            if  ~isempty(table(i,1)) && ~isempty(table(i,2)) && ~isempty(table(i,3)) && ~isempty(table(i,4))
                                Neighborhood_Set{index_Neighborhood} = Neighborhood(index_Neighborhood, [cell2mat(table2cell(table(i,4))), cell2mat(table2cell(table(i,5)))], cell2mat(table2cell(table(i,2))), cell2mat(table2cell(table(i,3))), cell2mat(table2cell(table(i,6))), cell2mat(table2cell(table(i,8))), cell2mat(table2cell(table(i,14))), cell2mat(table2cell(table(i,15))));
                                Dictionary(Neighborhood_Set{index_Neighborhood}.uniqueID) = Neighborhood_Set(index_Neighborhood);
                                index_Neighborhood = index_Neighborhood + 1;
                            else
                                message = strcat('ERROR: Empty Value or Value already exist in file: ', filename, ' found at line: ', num2str(i), ' column: 4');
                                error(message);
                            end
                        catch exception
                            msg = getReport(exception, 'basic');
                            disp(msg);
                            return;
                        end
                    end
    
            end
            
            power = {Branch_Set, Bus_Set, Generator_Set,{},{}};
            comm = {Centraloffice_Set, CommunicationTower_Set, Cellline_Set};
            trans = {Road_Set, Bridge_Set, TrafficLight_Set, RoadNode_Set,{}};
            
        end
        
        %% Creat Transmission Towers for carrying power lines (branch) among substations/power plants (generator), if there is no related input data available.
        function Power = CreateTransmissionTower(Power, Dictionary)
            % =================================================================
            % function CreateTransmissionTower
            % Create and generate transmissionTower based on substation and
            % branch data, based on a certain distance (d), which is
            % assumed as a constant here. 
            % =================================================================
            Branch_Set= Power{1};
            Bus_Set = Power{2};
            TransmissionTower_Set = Power{4};
            NewBranch_Set = {};
            
            TransmissionTower_index = 1;
            NewBranch_index = 1;
            % d: Assuming the distance of 2 adjacent transmission towers as 1 km
            %d = 1;            
            d = 5;
            
            for ii = 1:length(Branch_Set)
                temp = Branch_Set{ii};
                cap = temp.Capacity;
                type = temp.Type;
                proprity = temp.Priority;
                obj1 = temp.connectedObj1;
                obj2 = temp.connectedObj2;
                end_location = temp.End_Location;
                start_location = temp.Start_Location;
                [dis,az] = distance(temp.Start_Location(1), temp.Start_Location(2), temp.End_Location(1), temp.End_Location(2), referenceSphere('earth','km'));
                num = dis/d;
                % number of new branch due to adding transmission towers
                nbranch = floor(num); 
                %disp(num);
                if nbranch > 1
                    for j = 1:nbranch-1
                        %endTemp = [76, 40];
                        for k = 1:2  
                            endTemp(k) = start_location(k)+inv(nbranch)*j*(end_location(k)-start_location(k));
                        end

                        TransmissionTower_Set{TransmissionTower_index} = TransmissionTower(TransmissionTower_index,endTemp);
                        Dictionary(TransmissionTower_Set{TransmissionTower_index}.uniqueID) = TransmissionTower_Set(TransmissionTower_index);
                        if j == 1
                            NewBranch_Set{NewBranch_index} = Branch(NewBranch_index, start_location, endTemp, cap, type, proprity, obj1,  TransmissionTower_Set{TransmissionTower_index}.uniqueID);
                            busTemp = Dictionary(obj1);
                            busTemp{1}.Branch = [busTemp{1}.Branch, NewBranch_index];
                        else %if gt(j,1)*lt(j,nbranch)
                            NewBranch_Set{NewBranch_index} = Branch(NewBranch_index, TransmissionTower_Set{TransmissionTower_index - 1}.Location, endTemp, cap, type, proprity, TransmissionTower_Set{TransmissionTower_index - 1}.uniqueID,  TransmissionTower_Set{TransmissionTower_index}.uniqueID);
                            TransmissionTower_Set{TransmissionTower_index - 1}.Branch = [TransmissionTower_Set{TransmissionTower_index - 1}.Branch,NewBranch_index];
                        end
                        Dictionary(NewBranch_Set{NewBranch_index}.uniqueID) = NewBranch_Set(NewBranch_index);
                        TransmissionTower_Set{TransmissionTower_index}.Branch = [TransmissionTower_Set{TransmissionTower_index}.Branch,NewBranch_index];
                        NewBranch_index = NewBranch_index + 1;
                        TransmissionTower_index = TransmissionTower_index + 1;
                    end
                    % j = nbranch %floor(num), for debug purpose
                    NewBranch_Set{NewBranch_index} = Branch(NewBranch_index, TransmissionTower_Set{TransmissionTower_index - 1}.Location, end_location, cap, type, proprity, TransmissionTower_Set{TransmissionTower_index - 1}.uniqueID,  obj2);
                    Dictionary(NewBranch_Set{NewBranch_index}.uniqueID) = NewBranch_Set(NewBranch_index);
                            %NewBranch_index = NewBranch_index + 1;
                    busTemp = Dictionary(obj2);
                    busTemp{1}.Branch = [busTemp{1}.Branch, NewBranch_index];
                    NewBranch_index = NewBranch_index + 1;
                    
                    %end
                else %nbrach<=1
                    NewBranch_Set{NewBranch_index} = Branch_Set{ii};
                    busTemp = Dictionary(obj1);
                    busTemp{1}.Branch = [busTemp{1}.Branch, NewBranch_index];
                    busTemp = Dictionary(obj2);
                    busTemp{1}.Branch = [busTemp{1}.Branch, NewBranch_index];
                    NewBranch_index = NewBranch_index + 1;
                end
            end
            Power{1} = NewBranch_Set;
            Power{4} = TransmissionTower_Set;
        end
        
        %% count the number of damaged objects
        function number = countDamaged(Set)
            number = 0;
            for iobj = 1:length(Set)
                if strcmp(Set{iobj}.Status,'Damaged')
                    number = number + 1;
                end
            end
        end
        
        %% [Old] Get repairation time at the Object level
        function time = getRepairationTime(num, Object, Power_Set, Communication_Set, Transportation_Set)
            Branch= Power_Set{1};
            Bus= Power_Set{2};
            Generator= Power_Set{3};
            
            Centraloffice = Communication_Set{1};
            CommunicationTower = Communication_Set{2};
            Cellline = Communication_Set{3};
            
            Road = Transportation_Set{1};
            Bridge = Transportation_Set{2};
            TrafficLight = Transportation_Set{3};
            
            tem = strsplit(Object,'/');
            name = tem{1};
            number = str2num(tem{2});
            
            
            if strcmp(name, 'Branch')
                if num == 1
                    time = Branch{number}.Recovery(1);
                elseif num == 2
                    time = Branch{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'Bus')
                if num == 1
                    time = Bus{number}.Recovery(1);
                elseif num == 2
                    time = Bus{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'Generator')
                if num == 1
                    time = Generator{number}.Recovery(1);
                elseif num == 2
                    time = Generator{number}.WorkingDays;
                end
            end
            
            
            if strcmp(name, 'Centraloffice')
                if num == 1
                    time = Centraloffice{number}.Recovery(1);
                elseif num == 2
                    time = Centraloffice{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'Router')
                if num == 1
                    time = Router{number}.Recovery(1);
                elseif num == 2
                    time = Router{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'Cellline')
                if num == 1
                    time = Cellline{number}.Recovery(1);
                elseif num == 2
                    time = Cellline{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'CommunicationTower')
                if num == 1
                    time = CommunicationTower{number}.Recovery(1);
                elseif num == 2
                    time = CommunicationTower{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'Road')
                if num == 1
                    time = Road{number}.Recovery(1);
                elseif num == 2
                    time = Road{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'Bridge')
                if num == 1
                    time = Bridge{number}.Recovery(1);
                elseif num == 2
                    time = Bridge{number}.WorkingDays;
                end
            end
            
            if strcmp(name, 'TrafficLight')
                if num == 1
                    time = TrafficLight{number}.Recovery(1);
                elseif num == 2
                    time = TrafficLight{number}.WorkingDays;
                end
            end
        end
        
        %% Find out the dependent (precedent) task(s), given a task. 
        function result = getDependency(Task, Dictionary)
            %=============================================================
            % function getDependency
            % This function finds the dependent (precedent) task(s) for the
            % given task named "Object" from Dictionary. 
            % 
            %=============================================================
            result = 0;
            
            %=== Library.getUniqueId
            % uniqueId = -99, meaning that this given task string is not recognized.
            % Otherwise, meaning that this given task string is recognized.
            if Library.getUniqueId(Task,0) ~= -99
                taskIDstring = Library.getUniqueId(Task);
                
                TF = isKey(Dictionary, taskIDstring);
                %--- TF = 1, the "Key" is found in "Dictionary". Find out the
                % predecessor component (task) of the given task by looking
                % up "Task.predecessorComponent" (which may contraint one 
                % task or multiple tasks). After that, check whether there 
                % is any precendent task that is not finished yet 
                % (i.e., "Task.predecessorComponent{ii}.WorkingDays>0"). 
                %--- If there is any precendent task that is not finished yet, 
                % mark the flag variable of result = 1.
                if TF  
                    task = Dictionary(Library.getUniqueId(Task));
                    if iscell(task)
                        taskpre = task{1}.predecessorComponent;
                    else
                        taskpre = task.predecessorComponent;
                    end
                    for ii = 1:length(taskpre)
                        taskprei = Dictionary(taskpre(ii));
                        taskpreii = taskprei{1};
                        if any(taskpreii.WorkingDays(1:4))
                            result = 1;
                            return;
                        end
                    end
                    
                %--- TF = 0, cannot find the "Key" from "Dictionary"     
                else
                    %--- For debug 
                    %Library.getUniqueId(Object)  % debug
                    msg = strcat('----Function getDependency Error: Cannot find the Key of', ' ', Task,' in Dictionary.----');
                    disp(msg);
                    %return
                    %--- Setting result = 1 to continue the computation,
                    %disregard the above error (message) for now
                    result = -1;
                
                end
            %=== The form of this task named "Object" in this function is not recognized.    
            % Library.getUniqueId(Object,0) = -99
            else
                result = -2;
            end
        end
        
        %% Get the resource demand of a task by name and index
        function result = getResource(Task, Dictionary, resourceIndex, flag)
            %================================================================
            % function getResource
            % This function extracts the resource demand, given a task. 
            %================================================================
            % flag = 1, meaning this task string is in the form of
            % "Bus/2/Task23/Working/dummy", and temp = 'Task23'. 
            % flag = 0, meaning this task string is in the form of
            % "Bus/2/Task23", and temp = 'Task23'. 
            if flag 
                temp = Dictionary(Library.getUniqueId(Task, 1));
            else
                temp = Dictionary(Library.getUniqueId(Task, 0));
            end
            
            if iscell(temp)
                result = temp{1}.Resources(resourceIndex);
            else
                result = temp.Resources(resourceIndex);
                if iscell(result)
                    result = result{1};
                end
            end
            
        end
        
        %% Get UniqueId from schedule
        function uniqueId = getUniqueId(Task, flag)
            %================================================================
            % function getUniqueId
            % This function extracts the uniqueID, given a task. 
            %================================================================            
            % flag = 1, meaning this task string is in the form of
            % "Bus/2/Task23/Working/dummy", and temp = 'Task23'. 
            % flag = 0, meaning this task string is in the form of
            % "Bus/2/Task23", and temp = 'Task23'. 
            
            tem = strsplit(Task,'/');
            %--- the unique ID of a task in the form of 'Bus2' for Task='Bus/2'
            if length(tem) == 2 
                uniqueId = strcat(tem{1}, tem{2});
            %--- the unique ID of a task in the form of 'Task11' for Task='Bus/2/Task11'
            elseif length(tem) == 3
                uniqueId = tem{3};
            %--- the following two elseifs mean that the unique ID of a task 
            % in the form of either 'Bus2' for Object='Bus/2/Working/dummy' (length(tem) == 4)   
            % or Object='Bus/2/Task11/Working/dummy' (length(tem) == 5). 
            elseif length(tem) == 4 && flag
                uniqueId = strcat(tem{1}, tem{2});
            elseif length(tem) == 5 && flag
                uniqueId = tem{3};
            %--- The uniqueID string of this task is not recognized.     
            else
                uniqueId = -99;
            end
            return;
        end
        
        %% Groubi optimize function
        function [x, fval, exitflag] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub)
            % ===========================================================================
            % function intlinprog
            %INTLINPROG A mixed integer linear programming example using the
            %   Gurobi MATLAB interface
            %
            %   This example is based on the intlinprog interface defined in the
            %   MATLAB Optimization Toolbox. The Optimization Toolbox
            %   is a registered trademark of The MathWorks, Inc.
            %
            %   x = INTLINPROG(f,intcon,A,b) solves the problem:
            %
            %   minimize     f'*x
            %   subject to   A*x <= b
            %                x(j) integer, when j is in the vector
            %                intcon of integer constraints
            %
            %   x = INTLINPROG(f,intcon,A,b,Aeq,beq) solves the problem:
            %
            %   minimize     f'*x
            %   subject to     A*x <= b,
            %                Aeq*x == beq
            %                x(j) integer, where j is in the vector
            %                intcon of integer constraints
            %
            %   x = INTLINPROG(f,intcon,A,b,Aeq,beq,lb,ub) solves the problem:
            %
            %   minimize     f'*x
            %   subject to     A*x <= b,
            %                Aeq*x == beq,
            %          lb <=     x <= ub.
            %                x(j) integer, where j is in the vector
            %                intcon of integer constraints
            %
            %   You can set lb(j) = -inf, if x(j) has no lower bound,
            %   and ub(j) = inf, if x(j) has no upper bound.
            %
            %   [x, fval] = INTLINPROG(f, intcon, A, b) returns the objective value
            %   at the solution. That is, fval = f'*x.
            %
            %   [x, fval, exitflag] = INTLINPROG(f, intcon, A, b) returns an exitflag
            %   containing the status of the optimization. The values for
            %   exitflag and corresponding status codes are:
            %    2 - Solver stopped prematurely. Integer feasible point found.
            %    1 - Optimal solution found.
            %    0 - Solver stopped prematurely. No integer feasible point found.
            %   -2 - No feasible point found.
            %   -3 - Problem is unbounded.
            % ===========================================================================
            if nargin < 4
                error('intlinprog(f, intcon, A, b)')
            end
            
            if nargin > 8
                error('intlinprog(f, intcon, A, b, Aeq, beq, lb, ub)');
            end
            
            if ~isempty(A)
                n = size(A, 2);
            elseif nargin > 5 && ~isempty(Aeq)
                n = size(Aeq, 2);
            else
                error('No linear constraints specified')
            end
            
            if ~issparse(A)
                A = sparse(A);
            end
            
            if nargin > 4 && ~issparse(Aeq)
                Aeq = sparse(Aeq);
            end
            
            model.obj = f;
            model.vtype = repmat('C', n, 1);
            model.vtype(intcon) = 'I';
            
            if nargin < 5
                model.A = A;
                model.rhs = b;
                model.sense = '<';
            else
                model.A = [A; Aeq];
                model.rhs = [b; beq];
                model.sense = [repmat('<', size(A,1), 1); repmat('=', size(Aeq,1), 1)];
            end
            
            if nargin < 7
                model.lb = -inf(n,1);
            else
                model.lb = lb;
            end
            
            if nargin == 8
                model.ub = ub;
            end
            
            params.outputflag = 0;
            params.timelimit = 100;
            
            result = gurobi(model, params);

            if strcmp(result.status, 'OPTIMAL')
                exitflag = 1;
            elseif strcmp(result.status, 'INTERRUPTED')
                if isfield(result, 'x')
                    exitflag = 2;
                else
                    exitflag = 0;
                end
            elseif strcmp(result.status, 'INF_OR_UNBD')
                params.dualreductions = 0;
                result = gurobi(model, params);
                if strcmp(result.status, 'INFEASIBLE')
                    exitflag = -2;
                elseif strcmp(result.status, 'UNBOUNDED')
                    exitflag = -3;
                else
                    exitflag = nan;
                end
            else
                exitflag = nan;
            end

            if isfield(result, 'x')
                x = result.x;
            else
                x = nan(n,1);
            end
            
            if isfield(result, 'objval')
                fval = result.objval;
            else
                fval = nan;
            end
        end
        
        %% Create Hash map for task types
        function M = CreateHashTaskLibrary(matFile)
            % ===================================================================================
            % function CreateHashTaskLibrary
            % This function creates a hash table (or so called dictionary)
            % for all tasks in terms of a table table library of every type of object. 
            % ===================================================================================
            fpath = 'LibraryTask';
            myfolder = strcat(pwd,'/',fpath);
            mydir = fullfile(myfolder,matFile);
            tmp = load(mydir);
            num = tmp.num;
            keySet = [];
            valueSet = [];
            ntask = height(num);
            for ii = 1:ntask
                keySet(ii) = num{ii,1};
                temp = {num{ii,2},num{ii,3},num{ii,4},num{ii,5},num{ii,6},num{ii,7},num{ii,8},num{ii,9},num{ii,10},num{ii,11}};
                valueSet{ii} = temp;
            end
            M = containers.Map(keySet,valueSet);
        end
        
        %% Create Hash map for task types
        function M = CreateHashTaskPerDamageLevel(matFile)
            % ===================================================================================
            % function CreateHashTaskPerDamageLevel
            % This function creates a hash table (or so called dictionary)
            % for all tasks with precedence at every damage level for every type of object. 
            % ===================================================================================            
            fpath = 'LibraryTaskPerDamageLevel';
            myfolder = strcat(pwd,'/',fpath);
            mydir = fullfile(myfolder,matFile);
            tmp = load(mydir);
            a = tmp.a;
            keySet = [];
            valueSet = [];
            for ii = 1:length(a)
                keySet(ii) = ii;
                temp = {};
                for j = 1:length(a{ii})
                    for k = 1:length(a{ii}{j})
                        temp{j}{k} = a{ii}{j}{k};
                    end
                end
                valueSet{ii} = temp;
            end
            M = containers.Map(keySet,valueSet);
        end
        

        %% function Library.assignFragilityParameter
        function assignFragility(EventType,IM,IMx, IMy, IMmeta, Power_Set,  Transportation_Set, Communication_Set)
            %%=================================================================
            % function assignFragility
            % 1. Read MAT files of fragility data for every type of Object and
            % from the subdirectory "LibraryFragility".
            % 2. Assign paramters of fragility curves for every Object 
            % (median and beta in the lognormal cumulative density function)
            % 3. EventType: 1-Earthquake; 2-Hurricane Wind; 3-Flood(not implemented)  
            % =================================================================

            %==== Read Object data from system cell
            Branch = Power_Set{1};
            Bus= Power_Set{2};
            Generator= Power_Set{3};
            TransmissionTower = Power_Set{4};
            
            Centraloffice = Communication_Set{1};
            CommunicationTower = Communication_Set{2};
            Cellline = Communication_Set{3};     
            
            Road = Transportation_Set{1};
            Bridge = Transportation_Set{2};
            TrafficLight = Transportation_Set{3};
            
            %==== Define the subdirectory and keywords for finding fragility-related files 
            strPath = 'LibraryFragility';
            objName = {'PowerPlant','Substation','PowerLine','TransmissionTower',...
                'CommunicationTower','Cellline',...
                'Road','TrafficLight'};
            SubObjBridge = {'BridgeAbutment', ...
                                'BridgeAbutmentFound', ...
                                'BridgeApproachSlab', ...
                                'BridgeBearing', ...
                                'BridgeColumn', ...
                                'BridgeColumnFound', ...
                                'BridgeDeck', ...
                                'BridgeGirder'};
            
            switch EventType
                case 1 %Earthquake
                    % fragility data of Object 

                    str0 = 'fragility_earthquake_';
                    
                    for io = 1:length(objName)
                        obj = objName{io};
                        strName = strcat(str0, obj, '.mat'); 
                        strFull = fullfile(strPath,strName); % the fragility data file number
                        T = load(strFull); % load the fragility data file
                        im = T.IM; % intensity measure metric(s)
                        st = T.StructuralType; % list of string variable of structural type
                        meta = T.metadata;  % metadata of the fragility curve
                        fragility = T.fragility; % median and beta of the fragility curve(s) at slight, medium, extensive, complete damage state
                        
                        if strcmp(obj,'PowerPlant')
                            set = Generator;
                        elseif strcmp(obj,'Substation')
                            set = Bus;
                        elseif strcmp(obj,'PowerLine')
                            set = Branch;
                        elseif strcmp(obj,'TransmissionTower')
                            set = TransmissionTower;
                        elseif strcmp(obj,'CommunicationTower')
                            set = CommunicationTower; 
                        elseif strcmp(obj,'Cellline')
                            set = Cellline;
                        elseif strcmp(obj,'Road')
                            set = Road;   
                        elseif strcmp(obj,'TrafficLight')
                            set = TrafficLight;
                        end
                      
                        for jj = 1:length(set)
                            
                            structuraltype = set{jj}.Type;
                            idx_st = find(contains(st,structuraltype)); % index of the structure type (st) in the list
                            idx_im = find(contains(IMmeta,im)); % index of the intensity measure (im) metric in the list
                            % if the defined intensity measure (im) in the fragility file 
                            % cannot be found from the list of im metrics, then set
                            % idx_im = 1;
                            % if the defined structure type in the fragility file cannot be found from the list of structure types, then set
                            % idx_st = 1;
                            if isempty(idx_im)
                                idx_im = 1;
                                msg = strcat('Library.assignFragility Error: IM=',im, ' is not included in the fragility file (LibraryFragility) for ', obj, num2str(jj), ...
                                    '. The first IM in the IM list in hazard meta is chosen to continue the computation.');
                                disp(msg);
                            end
                            if isempty(idx_st)
                                idx_st = 1;
                                msg = strcat('Library.assignFragility Error: The structure type cannot be found in the fragility file (LibraryFragility) for ', obj, num2str(jj), ...
                                    '. The first one in the structure type list in the fragility file is used to continue the computation.');
                                disp(msg);
                            end
                            
                            % fragilityVector: frgility vector 
                            fragilityVector = fragility{idx_im(1)}(idx_st(1),:); 
                            set{jj}.FragilityIM{1} = idx_im(1);
                            set{jj}.FragilityIM{2} = IMmeta{idx_im(1)};
                            Library.singleFrag(set{jj}, fragilityVector);
                            
                        end
                    end

                    for jj = 1:length(Centraloffice)
                        set = Centraloffice;
                        obj = 'Centraloffice';
                        strName = strcat('fragility_earthquake_',obj,'.mat'); 
                        strFull = fullfile(strPath,strName);
                        T = load(strFull);
                        im = T.IM;
                        st = T.StructuralType;
                        meta = T.metadata; 
                        fragility = T.fragility;
                        code = T.code; 
                                
                        structuraltype = set{jj}.Type;
                        co = set{jj}.Code; 
                        idx_st = find(contains(st,structuraltype));
                        idx_im = find(contains(IMmeta,im));
                        idx_co = find(contains(code,co));

                        % fragilityVector: frgility vector 
                        fragilityVector = fragility{idx_im(1)}{idx_co(1)}(idx_st(1),:);    
                        set{jj}.FragilityIM{1} = idx_im(1);
                        set{jj}.FragilityIM{2} = IMmeta{idx_im(1)};
                        Library.singleFrag(set{jj}, fragilityVector);
                    end
                        
                    for ii = 1:length(Bridge)
                        if Bridge{ii}.HasSub == 0
                            %Library.singleFrag(set{ii}, fragility_Bridge);
                                    
                            set = Bridge;
                            obj = 'Bridge';
                            strName = strcat('fragility_earthquake_',obj,'.mat'); 
                            strFull = fullfile(strPath,strName);
                            T = load(strFull);
                            im = T.IM;
                            st = T.StructuralType;
                            meta = T.metadata; 
                            fragility = T.fragility;
                                    
                            for jj = 1:length(set)
                                structuraltype = set{jj}.Type;
                                idx_st = find(contains(st,structuraltype));
                                idx_im = find(contains(IMmeta,im)); 
                                if isempty(idx_im)
                                    idx_im = 1;
                                end
                                if isempty(idx_st)
                                    idx_st = 1;
                                end

                                % fragilityVector: fragility vector 
                                for k = 1:length(idx_im)
                                    fragilityVectorCell{k} = fragility{idx_im(k)}(idx_st(1),:);   
                                end
                                Library.multipleFrag(set{jj}, fragilityVectorCell);
                                set{jj}.FragilityIM{1} = idx_im;
                                for k = 1:length(idx_im)
                                    set{jj}.FragilityIM{2}{k} = IMmeta{idx_im(k)};
                                end

                             end                                    
                         end

                        if Bridge{ii}.HasSub == 1
                            % fragility data of the sub-object for the Bridge component
                            for isub = 1:length(SubObjBridge)
                                subobj = SubObjBridge{isub};
                                str0 = 'fragility_earthquake_';
                                strName = strcat(str0, subobj,'.mat'); 
                                strFull = fullfile(strPath,strName);
                                T = load(strFull);
                                im = T.IM;
                                st = T.StructuralType;
                                meta = T.metadata; 
                                fragility = T.fragility;
                                
                                if strcmp(subobj,'BridgeAbutment')
                                    set = Bridge{ii}.AbutmentSet;
                                elseif strcmp(subobj,'BridgeAbutmentFound')
                                    set = Bridge{ii}.AbutmentFoundSet;
                                elseif strcmp(subobj,'BridgeApproachSlab')
                                    set = Bridge{ii}.SlabSet;
                                elseif strcmp(subobj,'BridgeBearing')
                                    set = Bridge{ii}.BearingSet;
                                elseif strcmp(subobj,'BridgeColumn')
                                    set = Bridge{ii}.ColumnSet; 
                                elseif strcmp(subobj,'BridgeColumnFound')
                                    set = Bridge{ii}.ColumnFoundSet;
                                elseif strcmp(subobj,'BridgeDeck')
                                    set = Bridge{ii}.DeckSet;   
                                elseif strcmp(subobj,'BridgeGirder')
                                    set = Bridge{ii}.GirderSet;
                                end
                                
                                for jj = 1:length(set)
                                    structualtype = set(jj).Class;
                                    idx_st = find(contains(st,structualtype));
                                    idx_im = find(contains(IMmeta,im)); 
                                    if isempty(idx_im)
                                        idx_im = 1;
                                    end
                                    if isempty(idx_st)
                                        idx_st = 1;
                                    end

                                    % fragilityVector: frgility vector 
                                    fragilityVector = fragility{idx_im(1)}(idx_st(1),:); 
                                    set(jj).FragilityIM{1} = idx_im(1);
                                    set(jj).FragilityIM{2} = IMmeta{idx_im(1)};
                                    Library.singleFrag(set(jj), fragilityVector);

                                end
                            end

                                    
%                             for j = 1:length(Bridge{ii}.ColumnSet)
%                                 Library.singleFrag(Bridge{ii}.ColumnSet(j), fragility_BridgeColumn);
%                             end
%                             for j = 1:length(Bridge{ii}.ColumnFoundSet)
%                                 Library.singleFrag(Bridge{ii}.ColumnFoundSet(j), fragility_BridgeColumnFoundation);
%                             end
%                             for j = 1:length(Bridge{ii}.AbutmentSet)
%                                 Library.singleFrag(Bridge{ii}.AbutmentSet(j), fragility_BridgeAbutment);
%                             end
%                             for j = 1:length(Bridge{ii}.AbutmentFoundSet)
%                                 Library.singleFrag(Bridge{ii}.AbutmentFoundSet(j), fragility_BridgeAbutmentFoundation);
%                             end
%                             for j = 1:length(Bridge{ii}.GirderSet)
%                                  Library.singleFrag(Bridge{ii}.GirderSet(j), fragility_BridgeGirder);
%                             end
%                             for j = 1:length(Bridge{ii}.BearingSet)
%                                 Library.singleFrag(Bridge{ii}.BearingSet(j), fragility_BridgeBearing);
%                             end
%                             for j = 1:length(Bridge{ii}.DeckSet)
%                                 Library.singleFrag(Bridge{ii}.DeckSet(j), fragility_BridgeDeck);
%                             end
%                             for j = 1:length(Bridge{ii}.SlabSet)
%                                 Library.singleFrag(Bridge{ii}.SlabSet(j), fragility_BridgeApproachSlab);
%                             end
                        end  
                    end
                
                case 2 %Hurricane wind
                    % fragility data of Object 
                    
                    %%--------ground terrain type
                    % Terrain = {'Open','Light Suburban','Suburban','Light Urban','Urban'}; 
                    iter = 11; % the index of the IM variable (11th: ground terrain type) in the IM paramter of a wind map
                    if isempty(IM{iter})
                         idx_ter = 1;
                    elseif strcmp(IM{iter},'Open')
                         idx_ter = 1;
                    elseif strcmp(IM{iter},'Light Suburban')
                         idx_ter = 2;
                    elseif strcmp(IM{iter},'Suburban')
                         idx_ter = 3;
                    elseif strcmp(IM{iter},'Light Urban')
                         idx_ter = 4;     
                    elseif strcmp(IM{iter},'Urban')
                         idx_ter = 5;  
                    else 
                        errmsg = 'Input Error in the HazardMap: undefined ground terrain. Please revise this in the 11th cell of the IM paramter of the wind map and choose one of the five strings: Open, Light Suburban, Suburban, Light Urban, Urban';
                        disp(errmsg);
                    end
                    
                    objName = {'PowerPlant','Substation','PowerLine','TransmissionTower','CentralOffice','CommunicationTower','Cellline','Road','TrafficLight','Bridge'};
                    %objName = {'CentralOffice','CommunicationTower','Cellline','Road','TrafficLight','Bridge'};
                    str0 = 'fragility_wind_';
                    
                    for io = 1:length(objName)
                        obj = objName{io};
                        strName = strcat(str0,obj,'.mat'); 
                        strFull = fullfile(strPath,strName);
                        T = load(strFull);
                        im = T.IM;
                        st = T.StructuralType;
                        meta = T.metadata; 
                        fragility = T.fragility;
                        
                        if strcmp(obj,'PowerPlant')
                            set = Generator;
                        elseif strcmp(obj,'Substation')
                            set = Bus;
                        elseif strcmp(obj,'PowerLine')
                            set = Branch;
                        elseif strcmp(obj,'TransmissionTower')
                            set = TransmissionTower;
                        elseif strcmp(obj,'CentralOffice')
                            set = Centraloffice; 
                        elseif strcmp(obj,'CommunicationTower')
                            set = CommunicationTower;     
                        elseif strcmp(obj,'Cellline')
                            set = Cellline;
                        elseif strcmp(obj,'Road')
                            set = Road;   
                        elseif strcmp(obj,'TrafficLight')
                            set = TrafficLight;
                        elseif strcmp(obj,'Bridge')
                            set = Bridge; 
                        end
                      
                        for jj = 1:length(set)
                            idx_im = find(contains(IMmeta, im));
                            %%--------IM type
                            if isempty(idx_im)
                                errmsg = strcat('Library.assignFragility Error: IM=',im, ' is not included in the fragility file (LibraryFragility) for ', obj, num2str(jj),  ...
                                    '. The first IM in the IM list in hazard meta is chosen to continue the computation.');
                                disp(errmsg);
%                                 return
                                idx_im = 1;
                            end
                            
                            %%--------fragilityVectorCell: frgility vector in Cells
                            fragilityVectorCell = fragility{idx_im(1)}{idx_ter};
                            Library.assignObjectFragilityWind(set{jj}, fragilityVectorCell, IM, IMx, IMy, IMmeta, idx_im, idx_ter, st);
                 
                        end
                    end
                    
                    

             end
        end

        %% Assign paramters of fragility curves to every single object 
        function singleFrag(Object, fragVector)
            %%=================================================================
            % function singleFrag
            % Assing fragility parameters to every single object of all types
            % Object in the fragility curve using a single event intensity parameter. 
            % For instance, the fragility curve gives the probability of exceedance
            % of a structure with related to the event intensity, such as 
            % Peak Ground Acceleration (PGA) in a seismic scenario, 
            % or Sustained Wind Speed in a hurricane scneraio.
            % Fragility(rol,col) is in the form of a 4x2 matrix.
            % In this Fragility matrix, every row has two parameters for defining a
            % fragility curve, and four rows represent four fragility curves
            % for four possible damage states (1=slight, 2=moderate, 3=extensive, and 4= complete).
            % =================================================================
            for j = 1:8
                rol = idivide(j+1, int32(2), 'floor');
                col = rem(j,2);
                if col == 0
                    col = 2;
                end
                Object.Fragility(rol,col) = fragVector(j);
            end
        end
        
        %% Assign paramters of fragility curves when multiple paramters are used. 
        function multipleFrag(Object, fragVectorCell)
            Object.Fragility = [];
            for k = 1:length(fragVectorCell)
                for j = 1:8
                    rol = idivide(j+1, int32(2), 'floor');
                    col = rem(j,2);
                    if col == 0
                        col = 2;
                    end
                    
                    Object.Fragility{k}(rol,col) = fragVectorCell{k}(j);
                end
            end
        end
        
        
        %% Specific assignments by type
        function branchFrag(Object, fragVector)
            for j = 1:8
                rol = idivide(j+1, int32(2), 'floor');
                col = rem(j,2);
                if col == 0
                    col = 2;
                end
                switch Object.type
                    case 'null'
                        Object.Fragility(rol,col) = fragVector(j);
                end
            end
        end
        
        %% Assign fragility data to individual objects based on fragility curves under hurricane winds
        function assignObjectFragilityWind(Object, fragilityVectorCell, IM, IMx, IMy, IMmeta, idx_im, idx_ter, st)
            % ==============================================================================================  
            % function assignObjectFragilityWind(Object, fragilityVectorCell, IM, IMx, IMy, IMmeta, idx_im, idx_ter, st)
            % assign fragility data to every object under wind loads
            %%--------wind direction (the 10th cell in the "IM" parameter in the HazardWindMap)
            % ibeta: index of IM in IMmeta, representing the vector of wind
            % direction in the fragility analysis (ibeta = 10)
            % idx_dir: index of wind direction in the vector of wind direction (Beta) in the fragility analysis according to the object location
            % ==============================================================================================      
            
            ibeta = 10; % a numerical matrix of wind direction: the 10th index in the IMmeta parameter in the wind hazard map
            iter = 11; % a string cell of ground terrain: the 11th index in the IMmeta parameter in the wind hazard map
            if ~ismember(ibeta,idx_im)
                idx_dir = 1;
                Beta = 0;
            else
                tf = isprop(Object,'Location'); 
                if tf
                    location = Object.Location;
                else 
                    location = 0.5*(Object.Start_Location+Object.End_Location);
                end
                
                im = IM{ibeta}; 
                Intensity_beta = interp2(IMx,IMy,im,location(2),location(1));

                Beta = [0 15 30 45 60 75 90]; % direction thresholds of wind in the wind fragility curves 
                Beta2 = mean([Beta(1:end-1); Beta(2:end)]); % mean value of the adjacent wind direction angles in the Beta vector
                [~,c] = find(Beta2-Intensity_beta > 0);
                idx_dir = min(c); % the index of wind direction in the vector of wind direction 
            end

            %%--------structural type
            structualtype = Object.Type;
            idx_st = find(contains(st,structualtype));          
            if isempty(idx_st)
                errmsg = strcat('Library.assignObjectFragilityWind Warning: StructuralType=', structualtype, ' is not included in the fragility library.');
                disp(errmsg);
%                 return
                idx_st = 1;
            end
            tmp = fragilityVectorCell{idx_dir}(idx_st,:);
            
            Object.Fragility = reshape(tmp,2,4)';
            Object.FragilityIM{1,1} = idx_im(1);
            Object.FragilityIM{1,2} = IMmeta{idx_im(1)};
            Object.FragilityIM{2,1} = idx_ter(1);
            Object.FragilityIM{2,2} = IM{iter};
            Object.FragilityIM{3,1} = idx_dir;
            Object.FragilityIM{3,2} = Beta(idx_dir);

        end
        
        %% Assign Recovery Matrix for each object (Doesn't affect objects
        %other than generator after the implementation of task system)
        function assignRecovery (Power_Set,  Transportation_Set, Communication_Set)
            Branch = Power_Set{1};
            Bus = Power_Set{2};
            Generator = Power_Set{3};
            
            Centraloffice = Communication_Set{1};
            CommunicationTower = Communication_Set{2};
            Cellline = Communication_Set{3};
            
            Road = Transportation_Set{1};
            Bridge = Transportation_Set{2};
            TrafficLight = Transportation_Set{3};
            
            ros = load('restoration_Branch.mat');
            restoration_Branch = ros.restoration_Branch;
            ros = load('restoration_Bridge.mat');
            restoration_Bridge = ros.restoration_Bridge;
            ros = load('restoration_Bus.mat');
            restoration_Bus = ros.restoration_Bus;
            ros = load('restoration_CellLine.mat');
            restoration_CellLine = ros.restoration_CellLine;
            ros = load('restoration_CentralOffice.mat');
            restoration_CentralOffice = ros.restoration_CentralOffice;
            ros = load('restoration_CommunicationTower.mat');
            restoration_CommunicationTower = ros.restoration_CommunicationTower;
            ros = load('restoration_Generator.mat');
            restoration_Generator = ros.restoration_Generator;
            ros = load('restoration_PowerPlant.mat');
            restoration_PowerPlant = ros.restoration_PowerPlant;
            ros = load('restoration_Road.mat');
            restoration_Road = ros.restoration_Road;
            ros = load('restoration_Router.mat');
            restoration_Router = ros.restoration_Router;
            ros = load('restoration_Substation.mat');
            restoration_Substation = ros.restoration_Substation;
            ros = load('restoration_TrafficLight.mat');
            restoration_TrafficLight = ros.restoration_TrafficLight;
            ros = load('restoration_TransmissionTower.mat');
            restoration_TransmissionTower = ros.restoration_TransmissionTower;
            
            for i = 1:length(TrafficLight)
                Library.singleRec(TrafficLight{i}, restoration_TrafficLight);
            end

            for i = 1:length(Generator)
                Library.singleRec(Generator{i}, restoration_Generator);
            end
            for i = 1:length(CommunicationTower)
                Library.singleRec(CommunicationTower{i}, restoration_CommunicationTower);
            end
            for i = 1:length(Cellline)
                Library.singleRec(Cellline{i}, restoration_CellLine);
            end
            for i = 1:length(Bus)
                Library.singleRec(Bus{i}, restoration_Bus);
            end
            for i = 1:length(Branch)
                %                 Library.branchFrag(Branch{i}, fragility_Branch);
                Library.singleRec(Branch{i}, restoration_Branch);
            end
            for i = 1:length(Centraloffice)
                %                 Library.centralFrag(Centraloffice{i}, fragility_Centraloffice);
                Library.singleRec(Centraloffice{i}, restoration_CentralOffice);
            end
            for i = 1:length(Road)
                %                 Library.roadFrag(Road{i}, fragility_Road);
                Library.singleRec(Road{i}, restoration_Road);
            end
            for i = 1:length(Bridge)
                %                 Library.bridgeFrag(Bridge{i}, fragility_Bridge);
                Library.singleRec(Bridge{i}, restoration_Bridge);
                
                if Bridge{i}.HasSub == 1
                    for j = 1:length(Bridge{i}.ColumnSet)
                        Library.singleRec(Bridge{i}.ColumnSet(j), restoration_Bridge);
                    end
                    for j = 1:length(Bridge{i}.ColumnFoundSet)
                        Library.singleRec(Bridge{i}.ColumnFoundSet(j), restoration_Bridge);
                    end
                    for j = 1:length(Bridge{i}.AbutmentSet)
                        Library.singleRec(Bridge{i}.AbutmentSet(j), restoration_Bridge);
                    end
                    for j = 1:length(Bridge{i}.AbutmentFoundSet)
                        Library.singleRec(Bridge{i}.AbutmentFoundSet(j), restoration_Bridge);
                    end
                    for j = 1:length(Bridge{i}.GirderSet)
                        Library.singleRec(Bridge{i}.GirderSet(j), restoration_Bridge);
                    end
                    for j = 1:length(Bridge{i}.BearingSet)
                        Library.singleRec(Bridge{i}.BearingSet(j), restoration_Bridge);
                    end
                    for j = 1:length(Bridge{i}.DeckSet)
                        Library.singleRec(Bridge{i}.DeckSet(j), restoration_Bridge);
                    end
                    for j = 1:length(Bridge{i}.SlabSet)
                        Library.singleRec(Bridge{i}.SlabSet(j), restoration_Bridge);
                    end
                end         
                
            end
        end
        
        %% Helper function for recovery assignment, when there's only one
        % type for that object
        function singleRec(Object, recVector)
            for j = 1:8
                rol = idivide(j+1, int32(2), 'floor');
                col = rem(j,2);
                if col == 0
                    col = 2;
                end
                Object.RecoveryMatrix(rol,col) = recVector(j);
            end
        end
        
        %% Assign power functional dependecy for each objects in communication system
        function [Transportation_Set, Communication_Set] = assignPowerToTransComm(Power_Set, Transportation_Set, Communication_Set, Dictionary)
            
            Bus = Power_Set{2};
            
            TrafficLight = Transportation_Set{3};
            
            Centraloffice = Communication_Set{1};
            CommunicationTower = Communication_Set{2};
            Cellline = Communication_Set{3};            
            
            Bus_Location = [];
            
            for i = 1: length(Bus)
                Bus_Location = [Bus_Location; Bus{i}.Location];
            end
            
            for i = 1:length(TrafficLight)
                Idx = knnsearch(Bus_Location, TrafficLight{i}.Location,'K',1);
                TrafficLight{i}.Bus = Bus{Idx(1)}.uniqueID;
            end
            
            for i = 1:length(Centraloffice)
                Idx = knnsearch(Bus_Location, Centraloffice{i}.Location,'K',1);
                Centraloffice{i}.Bus = Bus{Idx(1)}.uniqueID;
            end
            
            for i = 1:length(Cellline)
                Idx = knnsearch(Bus_Location, Cellline{i}.Start_Location,'K',1);
                Cellline{i}.Bus = Bus{Idx(1)}.uniqueID;
            end
            
            for i = 1:length(CommunicationTower)
                Idx = knnsearch(Bus_Location, CommunicationTower{i}.Location,'K',1);
                CommunicationTower{i}.Bus = Bus{Idx(1)}.uniqueID;
            end
            
            Transportation_Set{3} = TrafficLight;
            
            Communication_Set{1} = Centraloffice;
            Communication_Set{2} = CommunicationTower;

        end
        
        %% Construct and link cellines between objects in communication system
        function Communication_Set = assignCellLine(Communication_Set, Dictionary)
            
            Centraloffice = Communication_Set{1};
            CommunicationTower = Communication_Set{2};
            Cellline_Set = Communication_Set{3};

            index_Cellline = 1;
            
            Central_location = [];
            
            for i = 1: length(Centraloffice)
                Central_location = [Central_location; Centraloffice{i}.Location];
            end

            %Connect Central office to Central Office
            
            if length(Centraloffice)==1
                Cellline_Set = {};
            elseif length(Centraloffice)==2
                i = 1; Idx = 2; j = 1;
                Cellline_Set{index_Cellline} = Cellline(index_Cellline, CentralOffice{i}.Location, CentralOffice{Idx(j)}.Location, CentralOffice{i}.uniqueID, CentralOffice{Idx(j)}.uniqueID);
                Dictionary(Cellline_Set{index_Cellline}.uniqueID) = Cellline_Set(index_Cellline);

                %record objects that got connected with
                CentralOffice{i}.Cellline = [CentralOffice{i}.Cellline, Cellline_Set{index_Cellline}.uniqueID];
                CentralOffice{Idx(j)}.Cellline = [CentralOffice{Idx(j)}.Cellline, Cellline_Set{index_Cellline}.uniqueID];
                CentralOffice{i}.CentralOffice = [CentralOffice{i}.CentralOffice, CentralOffice{Idx(j)}.uniqueID];
                CentralOffice{Idx(j)}.CentralOffice = [CentralOffice{Idx(j)}.CentralOffice, CentralOffice{i}.uniqueID];
   
            else % length(Centraloffice)>=3
            
                for i = 1:length(Centraloffice)

                    % nnn: number of nearest neighbor
                    nnn = 3;
                    if gt(nnn,length(Centraloffice))
                        nnn = length(Centraloffice)-1;
                    end

                    Idx = knnsearch(Central_location, Centraloffice{i}.Location,'K',nnn);
                    for j = 2:nnn
                        if length(Centraloffice{i}.CentralOffice) > 2
                            break;
                        end
                        if length(Centraloffice{Idx(j)}.CentralOffice) > 2
                            continue;
                        end
                        Cellline_Set{index_Cellline} = Cellline(index_Cellline, Centraloffice{i}.Location, Centraloffice{Idx(j)}.Location, Centraloffice{i}.uniqueID, Centraloffice{Idx(j)}.uniqueID);
                        Dictionary(Cellline_Set{index_Cellline}.uniqueID) = Cellline_Set(index_Cellline);

                        %record objects that got connected with
                        Centraloffice{i}.Cellline = [Centraloffice{i}.Cellline, Cellline_Set{index_Cellline}.uniqueID];
                        Centraloffice{Idx(j)}.Cellline = [Centraloffice{Idx(j)}.Cellline, Cellline_Set{index_Cellline}.uniqueID];
                        Centraloffice{i}.CentralOffice = [Centraloffice{i}.CentralOffice, Centraloffice{Idx(j)}.uniqueID];
                        Centraloffice{Idx(j)}.CentralOffice = [Centraloffice{Idx(j)}.CentralOffice, Centraloffice{i}.uniqueID];

                        index_Cellline = index_Cellline + 1;
                    end
                end
            end
            
            for i = 1:length(CommunicationTower)
                Idx = knnsearch(Central_location, CommunicationTower{i}.Location,'K',2);
                for j = 1:2
                    Cellline_Set{index_Cellline} = Cellline(index_Cellline, CommunicationTower{i}.Location, Centraloffice{Idx(j)}.Location, CommunicationTower{i}.uniqueID, Centraloffice{Idx(j)}.uniqueID);
                    Dictionary(Cellline_Set{index_Cellline}.uniqueID) = Cellline_Set(index_Cellline);
                    
                    Centraloffice{Idx(j)}.Cellline = [Centraloffice{Idx(j)}.Cellline, Cellline_Set{index_Cellline}.uniqueID];
                    CommunicationTower{i}.Cellline = [CommunicationTower{i}.Cellline, Cellline_Set{index_Cellline}.uniqueID];
                    Centraloffice{Idx(j)}.CommTower = [Centraloffice{Idx(j)}.CommTower, CommunicationTower{i}.uniqueID];
                    CommunicationTower{i}.Centraloffice = [CommunicationTower{i}.Centraloffice, Centraloffice{Idx(j)}.uniqueID];
                    index_Cellline = index_Cellline + 1;
                end
            end
            
            Communication_Set{3} = Cellline_Set;
            
        end
        
        %% Help function to compute the delayed task duration due to the significant disrution of transportation system
        % representing the fact that the severe damage in the transportation network causes restoration delays in the other two systems. 
        function SysFuncInterdependenceHelper(Set, Dictionary, DelayFactor)
            %==================================================================
            % function SysFuncInterdependenceHelper
            % This function adjusts the task duration based on the value
            % of DelayFactor from the input file, and the functionality of
            % transportation system.
            %=== Assumptions
            % If the inter-system functionality dependency is turned on,
            % apply the delaying effect. This is to say: 
            % (1) If Qtrans <= Qtrans0, TaskDuration = TaskDuration * DelayFactor
            % (2) Else (Qtrans > Qtrans0), TaskDuration = TaskDuration
            %==================================================================
            
            for ii = 1:length(Set) - 1 % the last cell in every SystemSet is Neighborhood
                
                for j = 1:length(Set{ii})
                    if strcmp(Set{ii}{1}.Class,'RoadNode')
                        break;
                    end
                    if ~strcmp(Set{ii}{j}.Status, 'Damaged')
                        continue;
                    end
                    tasks = Set{ii}{j}.taskUniqueIds;
                    sum = zeros(1,4);
                    for k = 1:length(tasks)
                        temp = Dictionary(tasks{k});
                        if length(temp.WorkingDays) == 1
                            a = temp.WorkingDays;
                            temp.WorkingDays(1:4) = [24*a, a, inv(7)*a, inv(28)*a];
                        end
                        
                        temp.WorkingDays(1:4) = temp.WorkingDays(1:4) * DelayFactor;
                        sum(1:4) = sum(1:4) + temp.WorkingDays(1:4);
                    end
                    Set{ii}{j}.WorkingDays = sum;
                end
            end
        end
        
        %% Modify the task duration if the InterdependenceTransDelay is turned on based on the transportation system functionality
        function SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, DelayFactor, Dictionary)
            Library.SysFuncInterdependenceHelper(Power_Set, Dictionary, DelayFactor);
            Library.SysFuncInterdependenceHelper(Communication_Set, Dictionary, DelayFactor);
            Library.SysFuncInterdependenceHelper(Transportation_Set, Dictionary, DelayFactor);
        end
        
        %% Create graph for transportation system based on road and roadnodes
        function G = BuildGraphTrans(Trans_Set, Dictionary, LinkDirection)
            
            hash = containers.Map('KeyType','double','ValueType','char');
            
            road_Set = Trans_Set{1};
            roadnode_Set = Trans_Set{4};
            Neighborhood = Trans_Set{5};
            for i = 1:length(roadnode_Set)
                hash(roadnode_Set{i}.nodeID) = roadnode_Set{i}.uniqueID;
            end
            s = [];
            t = [];
            weights = [];
            for i = 1:length(road_Set)
                s{i} = hash(road_Set{i}.Start_Node);
                t{i} = hash(road_Set{i}.End_Node);
                weights{i} = road_Set{i}.Length;
            end
            for i = 1:length(Neighborhood)
                s{i + length(road_Set)} = Neighborhood{i}.RoadNode;
                t{i + length(road_Set)} = Neighborhood{i}.Neighborhood;
                temp1 = Dictionary(Neighborhood{i}.Neighborhood);
                temp2 = Dictionary(Neighborhood{i}.RoadNode);
                temp1 = temp1{1};
                temp2 = temp2{1};
                weights{i + length(road_Set)} = distance(temp1.Location(1), temp1.Location(2), temp2.Location(1),temp2.Location(2),referenceSphere('earth','km'));
            end
            
            %==== Build the graph based on the choice of LinkDirection 
            % LinkDirection: 1-directional links, 0-not birectional (bidirectional) links
            if LinkDirection == 1 
                G = digraph(s,t,cell2mat(weights));
            elseif LinkDirection == 0 
                G = graph(s,t,cell2mat(weights));
            end
            
        end
        
        %% Create graph for communication system based on cellines and connected objects based on Graph Theory
        function G = BuildGraphComm(Comm_Set, Dictionary, LinkDirection)
            Cellline = Comm_Set{3};
            Neighborhood = Comm_Set{4};
            s = [];
            t = [];
            weights = [];
            for ii = 1:length(Cellline)
                s{ii} = Cellline{ii}.connectedObj1;
                t{ii} = Cellline{ii}.connectedObj2;
                temp1 = Dictionary(Cellline{ii}.connectedObj1);
                temp2 = Dictionary(Cellline{ii}.connectedObj2);
                temp1 = temp1{1};
                temp2 = temp2{1};
                weights{ii} = distance(temp1.Location(1), temp1.Location(2), temp2.Location(1),temp2.Location(2),referenceSphere('earth','km'));
            end
            
            for ii = 1:length(Neighborhood)
                s{ii + length(Cellline)} = Neighborhood{ii}.Neighborhood;
                t{ii + length(Cellline)} = Neighborhood{ii}.CentralTower;
                temp1 = Dictionary(Neighborhood{ii}.Neighborhood);
                temp2 = Dictionary(Neighborhood{ii}.CentralTower);
                temp1 = temp1{1};
                temp2 = temp2{1};
                weights{ii + length(Cellline)} = distance(temp1.Location(1), temp1.Location(2), temp2.Location(1),temp2.Location(2),referenceSphere('earth','km'));
            end
            
            %==== Build the graph based on the choice of LinkDirection 
            % LinkDirection: 1-directional links, 0-not birectional (bidirectional) links
            if LinkDirection == 1 
                G = digraph(s,t,cell2mat(weights));
            elseif LinkDirection == 0 
                G = graph(s,t,cell2mat(weights));
            end
            
        end
        
        %% Create graph for power system based on branches substations and transTowers
        function G = BuildGraphPower(Power_Set, Dictionary, LinkDirection)
            Branch = Power_Set{1};
            Neighborhood = Power_Set{5};
            s = [];
            t = [];
            weights = [];
            for ii = 1:length(Branch)
                s{ii} = Branch{ii}.connectedObj1;
                t{ii} = Branch{ii}.connectedObj2;
                temp1 = Dictionary(Branch{ii}.connectedObj1);
                temp2 = Dictionary(Branch{ii}.connectedObj2);
                temp1 = temp1{1};
                temp2 = temp2{1};
                weights{ii} = distance(temp1.Location(1), temp1.Location(2), temp2.Location(1),temp2.Location(2),referenceSphere('earth','km'));
            end
            
            for ii = 1:length(Neighborhood)
                s{ii + length(Branch)} = Neighborhood{ii}.Neighborhood;
                t{ii + length(Branch)} = Neighborhood{ii}.Bus;
                temp1 = Dictionary(Neighborhood{ii}.Neighborhood);
                temp2 = Dictionary(Neighborhood{ii}.Bus);
                temp1 = temp1{1};
                temp2 = temp2{1};
                weights{ii + length(Branch)} = distance(temp1.Location(1), temp1.Location(2), temp2.Location(1),temp2.Location(2),referenceSphere('earth','km'));
            end
            %==== Build the graph based on the choice of LinkDirection 
            % LinkDirection: 1-directional links, 0-not birectional (bidirectional) links
            if LinkDirection == 1 
                G = digraph(s,t,cell2mat(weights));
            elseif LinkDirection == 0 
                G = graph(s,t,cell2mat(weights));
            end

        end
        
        %% Link neighborhood with the nearest nodes for each system
        function [Transportation_Set, Communication_Set, Power_Set, Neighborhood_Set] = CreateNeighborhoodConnection(Power_Set, Transportation_Set,Communication_Set, Dictionary, Neighborhood_Set,active_power, active_trans, active_comm)
            
            Centraloffice = Communication_Set{1};

            Bus = Power_Set{2};
            roadNode = Transportation_Set{4};
            
            PowerLink_Set = [];
            TransLink_Set = [];
            CommLink_Set = [];
            
            index_PowerLink = 1;
            index_TransLink = 1;
            index_CommLink = 1;
            
            Central_location = [];
            Bus_location = [];
            roadNode_location = [];
            
            if active_power % Turn ON (1) the Power System
                for ii = 1:length(Bus)
                    Bus_location = [Bus_location; Bus{ii}.Location];
                end
                for ii = 1:length(Neighborhood_Set)
                    Idx = knnsearch(Bus_location, Neighborhood_Set{ii}.Location,'K',1);
                    PowerLink_Set{index_PowerLink} = Neighborhood_Power_Link(index_PowerLink, Neighborhood_Set{ii}.Location, Bus{Idx(1)}.Location, Neighborhood_Set{ii}.uniqueID, Bus{Idx(1)}.uniqueID);
                    Dictionary(PowerLink_Set{index_PowerLink}.uniqueID) = PowerLink_Set(index_PowerLink);

                    Neighborhood_Set{ii}.Bus = Bus{Idx(1)}.uniqueID;
                    Neighborhood_Set{ii}.Neighborhood_Power_Link = PowerLink_Set{index_PowerLink}.uniqueID;
                    Bus{Idx(1)}.Neighborhood{end+1}  = Neighborhood_Set{ii}.uniqueID;
                    Bus{Idx(1)}.PopulationServed = Bus{Idx(1)}.PopulationServed + Neighborhood_Set{ii}.Population;
                    Bus{Idx(1)}.Neighborhood_Power_Link{end+1}  = PowerLink_Set{index_PowerLink}.uniqueID;

                    index_PowerLink = index_PowerLink + 1;
                end          
                Power_Set{5} = PowerLink_Set;
            else
                Power_Set{5} = {};
            end 
           
            if active_comm % Turn ON (1) the Communication System
                for ii = 1: length(Centraloffice)
                    Central_location = [Central_location; Centraloffice{ii}.Location];
                end
                for ii = 1:length(Neighborhood_Set)
                    Idx = knnsearch(Central_location, Neighborhood_Set{ii}.Location,'K',1);
                    CommLink_Set{index_CommLink} = Neighborhood_Comm_Link(index_CommLink, Neighborhood_Set{ii}.Location, Centraloffice{Idx(1)}.Location, Neighborhood_Set{ii}.uniqueID, Centraloffice{Idx(1)}.uniqueID);
                    Dictionary(CommLink_Set{index_CommLink}.uniqueID) = CommLink_Set(index_CommLink);

                    Neighborhood_Set{ii}.Centraloffice = Centraloffice{Idx(1)}.uniqueID;
                    Neighborhood_Set{ii}.Neighborhood_Comm_Link = CommLink_Set{index_CommLink}.uniqueID;
                    Centraloffice{Idx(1)}.PopulationServed = Centraloffice{Idx(1)}.PopulationServed + Neighborhood_Set{ii}.Population;
                    Centraloffice{Idx(1)}.Neighborhood{end+1} = Neighborhood_Set{ii}.uniqueID;
                    Centraloffice{Idx(1)}.Neighborhood_Comm_Link{end+1} = CommLink_Set{index_CommLink}.uniqueID;

                    index_CommLink = index_CommLink + 1;
                end
                Communication_Set{4} = CommLink_Set;
            else
                Communication_Set{4} = {};
            end
           
            if active_trans % Turn ON (1) the Transportation System
                for ii = 1: length(roadNode)
                    roadNode_location = [roadNode_location; roadNode{ii}.Location];
                end
                for ii = 1:length(Neighborhood_Set)
                    Idx = knnsearch(roadNode_location, Neighborhood_Set{ii}.Location,'K',1);
                    TransLink_Set{index_TransLink} = Neighborhood_Trans_Link(index_TransLink, Neighborhood_Set{ii}.Location, roadNode{Idx(1)}.Location, Neighborhood_Set{ii}.uniqueID, roadNode{Idx(1)}.uniqueID);
                    Dictionary(TransLink_Set{index_TransLink}.uniqueID) = TransLink_Set(index_TransLink);

                    Neighborhood_Set{ii}.RoadNode = roadNode{Idx(1)}.uniqueID;
                    Neighborhood_Set{ii}.Neighborhood_Trans_Link = TransLink_Set{index_TransLink}.uniqueID;
                    roadNode{Idx(1)}.Neighborhood{end+1} = Neighborhood_Set{ii}.uniqueID;
                    roadNode{Idx(1)}.Neighborhood_Trans_Link{end+1} = TransLink_Set{index_TransLink}.uniqueID;

                    index_TransLink = index_TransLink + 1;
                end           
                Transportation_Set{5} = TransLink_Set;
            else
                Transportation_Set{5} = {};
            end   
            
        end
        
        %% Assign task(s) to the Object from the Object task library base on the Object.DamageLevel
        function [Set,taskIndex] = assignTask(Set, IndexObject, sumTaskHash, sumDamageTaskHash, Dictionary, taskIndex)
            % ================================================================
            % function assignTask
            % Assign restoration tasks for the specific damaged object based on its
            % task and Object.DamageLevel
            % ================================================================
            task_Set = {};
            predecessorTask = [];
            parentType = Set{IndexObject}.Class;
            parentUniqueId = Set{IndexObject}.uniqueID;
            DamageTempTaskHash = sumDamageTaskHash('Bus');
            tempTaskHash = sumTaskHash('Bus');
                    
            switch Set{IndexObject}.Class
                case 'Bus'
                    DamageTempTaskHash = sumDamageTaskHash('Bus');
                    tempTaskHash = sumTaskHash('Bus');
                case 'Branch'
                    DamageTempTaskHash = sumDamageTaskHash('Branch');
                    tempTaskHash = sumTaskHash('Branch');
                case 'CentralOffice'
                    DamageTempTaskHash = sumDamageTaskHash('CentralOffice');
                    tempTaskHash = sumTaskHash('CentralOffice');
                case 'CommunicationTower'
                    DamageTempTaskHash = sumDamageTaskHash('CommunicationTower');
                    tempTaskHash = sumTaskHash('CommunicationTower');
                case 'CommLine'
                    DamageTempTaskHash = sumDamageTaskHash('CommLine');
                    tempTaskHash = sumTaskHash('CommLine');
                case 'TransmissionTower'
                    DamageTempTaskHash = sumDamageTaskHash('TransmissionTower');
                    tempTaskHash = sumTaskHash('TransmissionTower');
                case 'Bridge'
                    DamageTempTaskHash = sumDamageTaskHash('Bridge');
                    tempTaskHash = sumTaskHash('Bridge');
                case 'Road'
                    DamageTempTaskHash = sumDamageTaskHash('Road');
                    tempTaskHash = sumTaskHash('Road');
                case 'TrafficLight'
                    DamageTempTaskHash = sumDamageTaskHash('TrafficLight');
                    tempTaskHash = sumTaskHash('TrafficLight');
                otherwise
                    msg = strcat('Function assignTask Error: The Object Class is not recognized:', Set{IndexObject}.Class, '.\n')
                    fprintf(msg);
                    return
            end
            
            %=== tasks{L1}{L2}{L3} has 3 layers of cell
            % L1: in the 1st layer, the length of L1 represents possible
            % ways to repair the same damage state.
            % L2: in the 2nd layer, the length of L2 represents the
            % precedence relations between tasks in L3. The 2nd cell of L2
            % is a succesor of the 1st cell of L2, and so forth. 
            % L3: in the 3rd layer, numerical values in L3 represent the
            % task number in the task library of the corresponding object type.
            tasks = DamageTempTaskHash(Set{IndexObject}.DamageLevel);

            %=== find out how many ways avaiable to repair the same damage
            %with one out of two alternative ways:
            %--- (1)choose one mode randomly from multiple posssible modes
%             if length(tasks) == 1 % only 1 mode
%                 tasks = tasks{1};
%             else  % multiple modes
%                 x = max(round(rand * length(tasks)), 1);
%                 tasks = tasks{x};
%             end
            %--- alternatively, (2)always choose the 1st way (i.e., the first mode)
            tasks = tasks{1};
            
            for j = 1:length(tasks) % No. of precedence relations
                for k = 1:length(tasks{j}) % No. of task in one precedence relation
                    idNeed = tasks{j}{k};
                    tempTask = tempTaskHash(idNeed);
                    taskDurationMode = tempTask{5};
                    if iscell(taskDurationMode)
                        taskDurationMode = taskDurationMode{1};
                    end
                        
                    task = Task(taskIndex, idNeed,[tempTask{6},tempTask{7},tempTask{8},tempTask{9}],tempTask{2},tempTask{3},tempTask{4},taskDurationMode,tempTask{1},tempTask{10},parentType,parentUniqueId);
                    tempPredecessorTask{k} = task.uniqueID;
                    task.predecessorTask = predecessorTask;
                    Set{IndexObject}.taskUniqueIds{end+1} = task.uniqueID;
                    Dictionary(task.uniqueID) = task;
                    task_Set{end + 1} = task;
                    taskIndex = taskIndex + 1;
                end
                predecessorTask = horzcat(predecessorTask,tempPredecessorTask);
            end
        end
        
        %% Set up a hash table of tasks for evety type of object
        function  [sumTaskHash, sumDamageTaskHash] = SetUpHashTables()
            % ================================================================
            % function SetUpHashTables
            % This function reads AMT files in libraries for task and hash them together
            % ================================================================

            %==== Power system
            busTaskHash = Library.CreateHashTaskLibrary('busTasks.mat');
            branchTaskHash = Library.CreateHashTaskLibrary('branchTasks.mat');
            transTowerTaskHash = Library.CreateHashTaskLibrary('transTowerTasks.mat');
            %==== Communciation system
            centralOfficeTaskHash = Library.CreateHashTaskLibrary('centralOfficeTasks.mat');
            commTowerTaskHash = Library.CreateHashTaskLibrary('commTowerTasks.mat');
            commLineTaskHash = Library.CreateHashTaskLibrary('commLineTasks.mat');
            %==== Transportation system
            bridgeTaskHash = Library.CreateHashTaskLibrary('bridgeTasks.mat');
            roadTaskHash = Library.CreateHashTaskLibrary('roadTasks.mat');
            trafficLightTaskHash = Library.CreateHashTaskLibrary('trafficLightTasks.mat');
            %==== Sub-Object analysis
            subbridgeTaskHash =Library.CreateHashTaskLibrary('subBridgeTasks.mat');
            
            %==== create a task dictionary for every type of object
            keySet = {'Bridge', 'Road', 'TrafficLight', 'Bus', 'Branch','TransmissionTower','CentralOffice', 'CommunicationTower', 'CommLine'};
            newKeySet = [keySet,'subBridge'];
            valueSet={bridgeTaskHash,roadTaskHash,trafficLightTaskHash,busTaskHash,branchTaskHash,transTowerTaskHash, centralOfficeTaskHash,commTowerTaskHash,commLineTaskHash, subbridgeTaskHash};
            sumTaskHash = containers.Map(newKeySet,valueSet);
            
            %==== create a damage dictionary for every type of object
            % At the 1st level, there are 4 cells. Every cell presents task(s) at a damage state. 
            % At the 2nd level, the cell(s) present the order of all tasks following the precedence relations from left to right. 
            DamageBridgeTaskHash = Library.CreateHashTaskPerDamageLevel('DamageBridgeTasks');
            DamageRoadTaskHash = Library.CreateHashTaskPerDamageLevel('DamageRoadTasks');
            DamageTrafficLightTaskHash = Library.CreateHashTaskPerDamageLevel('DamageTrafficLightTasks');
            DamageBusTaskHash = Library.CreateHashTaskPerDamageLevel('DamageBusTasks');
            DamageBranchTaskHash = Library.CreateHashTaskPerDamageLevel('DamageBranchTasks');
            DamageTransTowerTaskHash = Library.CreateHashTaskPerDamageLevel('DamageTransTowerTasks');
            DamageCentralOfficeTaskHash = Library.CreateHashTaskPerDamageLevel('DamageCentralOfficeTasks');
            DamageCommTowerTaskHash = Library.CreateHashTaskPerDamageLevel('DamageCommTowerTasks');
            DamageCommLineTaskHash = Library.CreateHashTaskPerDamageLevel('DamageCommLineTasks');
            DamageTransTowerTaskHash = Library.CreateHashTaskPerDamageLevel('DamageTransTowerTasks');
            
            %====create a damage dictionary for every subcomponent of bridge
            DamageAbutmentTaskHash = Library.CreateHashTaskPerDamageLevel('DamageAbutment.mat');
            DamageAbutmentFoundationTaskHash = Library.CreateHashTaskPerDamageLevel('DamageAbutmentFoundation.mat');
            DamageApproachSlabTaskHash = Library.CreateHashTaskPerDamageLevel('DamageApproachSlab.mat');
            DamageBearingTaskHash = Library.CreateHashTaskPerDamageLevel('DamageBearing.mat');
            DamageColumnTaskHash = Library.CreateHashTaskPerDamageLevel('DamageColumn.mat');
            DamageColumnFoundationTaskHash = Library.CreateHashTaskPerDamageLevel('DamageColumnFoundation.mat');
            %DamageGirderTaskHash = Library.CreateHashTaskPerDamageLevel('DamageGirder.mat');

            newKeySet = [keySet,'Abutment','AbutmentFoundation','ApproachSlab','Bearing','Column','ColumnFoundation'];
            valueSet = {DamageBridgeTaskHash,DamageRoadTaskHash,DamageTrafficLightTaskHash,DamageBusTaskHash,DamageBranchTaskHash,DamageCentralOfficeTaskHash,DamageCommTowerTaskHash,DamageCommLineTaskHash,DamageTransTowerTaskHash,DamageAbutmentTaskHash,DamageAbutmentFoundationTaskHash,DamageApproachSlabTaskHash,DamageBearingTaskHash,DamageColumnTaskHash,DamageColumnFoundationTaskHash};
            sumDamageTaskHash = containers.Map(newKeySet,valueSet);
            
        end
        
        %% Find out the duration of a task from the Dictionary
        function Set = getWorkDays(Set, ii, Dictionary)
            %===============================================================
            % function getWorkDays
            % This function extracts the task duration, given a task, based
            % on the task uniqueID in the Dictionary. 
            %===============================================================
            tasks = Set{ii}.taskUniqueIds;
            sumWorkDay = 0;
            flag = 0;
            for j = 1:length(tasks)
                temp = Dictionary(tasks{j});
                if iscell(temp)
                    temp = temp{1};
                end
                if temp.WorkingDays > 0
                    if flag == 0
                        Set{ii}.currentWorking = temp.uniqueID;
                        Set{ii}.Functionality = temp.taskFunctionality;
                        flag = 1;
                    end
                    sumWorkDay = sumWorkDay + temp.WorkingDays;
                end
            end
            Set{ii}.WorkingDays = sumWorkDay;
        end
        
        %% Sample a recovery duration for one single task
        function [samples] = simulatervLHS(samples,type,par)
            %===============================================================
            % function simulatervLHS
            % This function maps the samples following unique distributions to 
            % new samples following other distributions based on "type" (distribution type) 
            % and "par" (distribution parameters).  
            %===============================================================            
            switch char(type)
                case 'dete'
                    samples=ones(length(samples),1)*par(1);
                case 'norm'
                    samples = norminv(samples,par(1),par(2));
                case 'logn'
                    sigmaN=sqrt(log(1+(par(2)/par(1))^2));
                    muN=log(par(1))-0.5*sigmaN^2;
                    samples = logninv(samples,muN,sigmaN);
                case 'Uniform/2'
                    samples=samples*(par(2)-par(1))+par(1);
                case 'triangular/1'
                    caso=samples<(par(3)-par(1))/(par(2)-par(1));
                    samples=(par(1)+sqrt(samples*(par(2)-par(1))*(par(3)-par(1))))    .*caso +...
                        (par(2)-sqrt((1-samples)*(par(2)-par(1))*(par(2)-par(3)))).*(~caso);
                case 'Triangular/1'
                    caso=samples<(par(3)-par(1))/(par(2)-par(1));
                    samples=(par(1)+sqrt(samples*(par(2)-par(1))*(par(3)-par(1))))    .*caso +...
                        (par(2)-sqrt((1-samples)*(par(2)-par(1))*(par(2)-par(3)))).*(~caso);
                case 'empi'
                    [ecdf_y,ecdf_x] = ecdf(par);
                    samples = interp1(ecdf_y,ecdf_x,samples,'linear');
                    
                otherwise
                    errordlg('This distribution has not been implemented yet','Code terminated')
            end
        end
        
        %% Creat a table of tasks based on Dictionary
        % Older function?
        function [taskTable,index] = createTaskTable(Dictionary)
            %===============================================================
            % function createTaskTable
            % This function creates a table of all tasks for restoring all 
            % damaged components/objects, based on Dictionary.   
            %===============================================================     
            key = keys(Dictionary);
            taskTable = {};
            for ii = 1:length(key)
                temp = Dictionary(key{ii});
                if iscell(temp)
                    temp = temp{1};
                end
                if strcmp(temp.Class, 'Task')
%                      if ~exist('objTemp','var')
%                         objTemp = temp;
%                      end
%                     taskTable(end + 1,:) = {objTemp.uniqueID,objTemp.parentType,objTemp.parentUniqueID,objTemp.taskID,objTemp.taskDescription,objTemp.durationMin,objTemp.durationMax,objTemp.durationMode,objTemp.durationType,objTemp.WorkingDays,objTemp.Resources{1},objTemp.Resources{2},objTemp.Resources{3},objTemp.Resources{4}};
                    if iscell(temp.Resources)
                        resourceDemand = cell2mat(temp.Resources);
                    else
                        resourceDemand = temp.Resources;
                    end
                    taskTable(end + 1,:) = {temp.uniqueID,temp.parentType,temp.parentUniqueID,temp.taskID,temp.taskDescription,temp.durationMin,temp.durationMax,temp.durationMode,temp.durationType,temp.WorkingDays,resourceDemand};                
                end
            end
            task_power=[];power_indx=[];task_Communication=[];Communication_index=[];task_Transportation=[];Transportation_indx=[];
            for i = 1:size (taskTable,1)
                
                if strcmp(taskTable{i,2},'Branch') | strcmp(taskTable{i,2},'Bus') | strcmp(taskTable{i,2},'TransmissionTower')
                    task_power=[task_power;{taskTable{i,:}}];power_indx=[power_indx,i];
                end
                
                if strcmp(taskTable{i,2},'CentralOffice') | strcmp(taskTable{i,2},'Cellline') | strcmp(taskTable{i,2},'CommunicationTower')
                    task_Communication=[task_Communication;{taskTable{i,:}}];Communication_index=[Communication_index,i];
                end                
            
                 if strcmp(taskTable{i,2},'Road') | strcmp(taskTable{i,2},'Bridge') | strcmp(taskTable{i,2},'TrafficLight')
                    task_Transportation=[task_Transportation;{taskTable{i,:}}];Transportation_indx=[Transportation_indx,i];
                 end  
            end
            taskTable={task_power,task_Communication,task_Transportation}; index={power_indx,Communication_index,Transportation_indx};
        end
        
        %% Creat a table of tasks for an individual type of Object based on Dictionary 
        function taskTable = createTaskTableIndividual(Dictionary,sys)
            taskTable = {};
            if strcmp(sys, 'Power')
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'Branch');
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'Bus');
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'TransmissionTower');
            elseif strcmp(sys, 'Communication')
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'CentralOffice');
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'Cellline');
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'CommunicationTower');    
            elseif strcmp(sys, 'Transportation')
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'Road');
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'Bridge');
                taskTable = Library.createTaskTableIndividualHelper(taskTable, Dictionary , 'TrafficLight');

            end
        end

        
        %% Creat a table of tasks for an individual type of Object based on Dictionary to support the function of "createTaskTableIndividual"
        function taskTable = createTaskTableIndividualHelper(taskTable, Dictionary, objectClass)
            % =================================================================
            % function createTaskTableIndividualHelper
            % This function creates a table of all tasks for restoring all 
            % damaged components/objects, based on Dictionary.  
            % =================================================================
            key = keys(Dictionary);
            for ii = 1:length(key)
                temp = Dictionary(key{ii});
                if iscell(temp)
                    temp = temp{1};
                end
                if strcmp(temp.Class, 'Task') && strcmp(temp.parentType, objectClass)
                    %resourceDemand = temp.Resources;
                    if iscell(temp.Resources)
                        resourceDemand = cell2mat(temp.Resources);
                    else
                        resourceDemand = temp.Resources;
                    end
                    taskTable(end + 1,:) = {temp.uniqueID,temp.parentType,temp.parentUniqueID,temp.taskID,temp.taskDescription,temp.durationMin,temp.durationMax,temp.durationMode,temp.durationType,temp.WorkingDays,resourceDemand};
                end
            end
        end
    
        
        %% create the task list (table) named "lookupTable" in the process of adding tasks to the working list [through the function of AddCrruentWorking]
        function lookupTable = addToLookupTable(lookupTable, task, startDate, isTask)
            %=================================================================
            % function addToLookupTable
            % add the task of "objTemp" to the task list (table) named "lookupTable"
            % The "lookupTable" is a table shows the task information. In
            % this table, the column attributes from left to right are:
            % 1. task uniqueID
            % 2. parentType
            % 3. parentUniqueID
            % 4. taskID
            % 5. taskDescription
            % 6. taskDurationMin
            % 7. taskDurationMax
            % 8. taskDurationMode
            % 9. taskDurationType
            % 10.startDate
            % 11.objTemp.WorkingDays
            % 12.taskResourceDemand
            %=================================================================
            if isTask ~= 1
                task = task{1};
            end
            if strcmp(task.Class,'Task')      
                if iscell(task.Resources)
                    resourceDemand = cell2mat(task.Resources);
                else
                    resourceDemand = task.Resources;
                end
                lookupTable(end + 1,:) = {task.uniqueID, task.parentType, task.parentUniqueID, task.taskID, task.taskDescription, task.durationMin, task.durationMax, task.durationMode, task.durationType, startDate, task.WorkingDays, resourceDemand};
            else
                if iscell(task.Resources)
                    resourceDemand = cell2mat(task.Resources);
                else
                    resourceDemand = task.Resources;
                end
                lookupTable(end + 1,:) = {task.uniqueID, -1, '', '', '', -1, -1, -1, '',startDate,task.WorkingDays,resourceDemand};
            end
        end

        
        %% New function: Set up functionality hash table 
        function funcTable = initialFuncTableHelper(Set, time_horizon, funcTable)
            %===============================================================
            % function initialFuncTableHelper
            % This function sets up the initial functionality table (a dictionary)
            % for three systems based on the object status. 
            % If Object.Status = 'Damaged'/'Stoped',
            % funcTable(Object.uniqueID) = zeros(1,time_horizon). 
            %===============================================================
            for ii = 1:length(Set) - 1 % ii: index of object type in the system set
                for j = 1:length(Set{ii}) % j: index of object in the object set
                    if strcmp(Set{ii}{1}.Class,'RoadNode')
                        break;
                    end
                    if strcmp(Set{ii}{j}.Status, 'Damaged') || strcmp(Set{ii}{j}.Status, 'Stoped')
                        funcTable(Set{ii}{j}.uniqueID) = 0;
                    end
                end
            end
            
        end

        
        %%
        function printedFuncTable = printFuncTable(funcTable, time_horizon)
            key = keys(funcTable);
            printedFuncTable = cell(length(key),time_horizon + 1);
            for i = 1:length(key)
                printedFuncTable{i,1} = key{i};
                temp = funcTable(key{i});
                for j = 1:length(temp)
                    printedFuncTable{i,j + 1} = temp(j);
                end
            end
        end
        
        %%
        function [taskSet,indexTable] = countTask(Dictionary)
            indexTable = containers.Map('KeyType','char','ValueType','double');
            taskSet = {};
            key = keys(Dictionary);
            for i = 1:length(key)
                temp = Dictionary(key{i});
                if iscell(temp)
                    temp = temp{1};
                end
                if strcmp(temp.Class, 'Task')
                    taskSet{end + 1} = temp.uniqueID;
                    indexTable(temp.uniqueID) = length(taskSet);
                end
            end
        end
        
        %%
        function [taskSet,indexTable] = countTaskIndividual(Dictionary, class,taskSet,indexTable)
            key = keys(Dictionary);
            for ii = 1:length(key)
                temp = Dictionary(key{ii});
                if iscell(temp)
                    temp = temp{1};
                end
                if strcmp(temp.Class, 'Task') && strcmp(temp.parentType, class)
                    taskSet{end + 1} = temp.uniqueID;
                    indexTable(temp.uniqueID) = length(taskSet);
                end
            end
        end
        
        %% Creat a table of task precedence
        function precedenceTable = createPreTable(Dictionary,index, active_power, active_comm, active_trans)
            
            [taskSet,indexTable] = Library.countTask(Dictionary);
            taskCount = length(taskSet);
            precedenceTable = cell(taskCount + 1, taskCount + 1);
            for ii = 1:taskCount
                temp = taskSet{ii};
                precedenceTable{1, ii+1} = temp;
                precedenceTable{ii+1, 1} = temp;
            end
            for ii = 1:taskCount
                temp = Dictionary(taskSet{ii});
                if iscell(temp)
                    temp = temp{1};
                end
                precedenceTable(ii+1,2:end) = {0};
                if ~isempty(temp.predecessorTask)
                    tasks = temp.predecessorTask;
                    for j = 1:length(tasks)
                        index2 = indexTable(tasks{j});
                        precedenceTable{ii+1, index2 + 1} = 1;
                    end
                end
            end
            tmp={precedenceTable(2:end,2:end)};
            all=cell2mat(tmp{1});
            power_num=all(index{1},index{1});
            precedenceTable_power=[{precedenceTable{[1,index{1}+1],1}}',[{precedenceTable{1,index{1}+1}};num2cell(power_num)]];
            comm_num=all(index{2},index{2});
            precedenceTable_comm=[{precedenceTable{[1,index{2}+1],1}}',[{precedenceTable{1,index{2}+1}};num2cell(comm_num)]];
            
            trans_num=all(index{3},index{3});
            precedenceTable_trans=[{precedenceTable{[1,index{3}+1],1}}',[{precedenceTable{1,index{3}+1}};num2cell(trans_num)]];     
            if ~active_power
                precedenceTable_power = [];
            end
            
            if ~active_comm
                precedenceTable_comm = [];
            end
            
            if ~active_trans
                precedenceTable_trans = [];
            end
            
            precedenceTable={precedenceTable_power,precedenceTable_comm,precedenceTable_trans};
        end
        
        %% Creat a table of task precedence
        function precedenceTable = createPreTableIndividual(Dictionary, sys)
            indexTable = containers.Map('KeyType','char','ValueType','double');
            taskSet = {};
            if strcmp(sys, 'Power')
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'Branch',taskSet,indexTable);
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'Bus',taskSet,indexTable);
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'TransmissionTower',taskSet,indexTable);
            elseif strcmp(sys, 'Communication')
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'CentralOffice',taskSet,indexTable);
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'Router',taskSet,indexTable);
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'Cellline',taskSet,indexTable);
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'CommunicationTower',taskSet,indexTable);
            elseif strcmp(sys, 'Transportation')
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'Road',taskSet,indexTable);
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'Bridge',taskSet,indexTable);
                [taskSet,indexTable] = Library.countTaskIndividual(Dictionary, 'TrafficLight',taskSet,indexTable);
            end
            taskCount = length(taskSet);
            precedenceTable = cell(taskCount + 1, taskCount + 1);
            for ii = 1:taskCount
                temp = taskSet{ii};
                precedenceTable{1, ii+1} = temp;
                precedenceTable{ii+1, 1} = temp;
            end
            for ii = 1:taskCount
                temp = Dictionary(taskSet{ii});
                if iscell(temp)
                    temp = temp{1};
                end
                precedenceTable(ii+1,2:end) = {0};
                if ~isempty(temp.predecessorTask)
                    tasks = temp.predecessorTask;
                    for j = 1:length(tasks)
                        index = indexTable(tasks{j});
                        precedenceTable{ii+1, index + 1} = 1;
                    end
                end
            end
        end
        
        %% Count the total population at neighborhood level 
        function totalPopulation = countPopulation(Dictionary)
            key = keys(Dictionary);
            totalPopulation = 0;
            for i = 1:length(key)
                temp = Dictionary(key{i});
                if iscell(temp)
                    temp = temp{1};
                end
                if strcmp(temp.Class, 'Neighborhood')
                    totalPopulation = totalPopulation + temp.Population;
                end
            end
        end
        
        %% compute the functionality of neighborhood
        function [totalPopulation, QNbrPower, QNbrComm, QNbrTrans] = neighbourFunc(Dictionary)
            % ===================================================================
            % function neighbourFunc
            % This function computes the functionality of different individual systems based on
            % the percentage of neighborhood population with that service (mean the corresponding status=1).
            % QNbrPower = Number of population at neighborhoods with the power service / total population at all neighborhoods 
            % QNbrCommu = Number of population at neighborhoods with the communication service / total population at all neighborhoods 
            % QNbrTrans = Number of population at neighborhoods with the transportation service / total population at all neighborhoods 
            % ===================================================================
            %==== Field
            % keys in the Dirctionary
            key = keys(Dictionary); 
            
            % initial variables
            totalPopulation = 0; % total population of all neighborhoods
            totalHasPower = 0; % total population of neighborhoods with the power service
            totalHasComm = 0; % total population of neighborhoods with the communication service
            totalHasTrans = 0; % total population of neighborhoods with the transportation service
            
            for ii = 1:length(key)
                temp = Dictionary(key{ii});
                if iscell(temp)
                    temp = temp{1};
                end
                if strcmp(temp.Class, 'Neighborhood')
                    if temp.PowerStatus == 1
                        totalHasPower = totalHasPower + temp.Population;
                    end
                    if temp.CommStatus == 1
                        totalHasComm = totalHasComm + temp.Population;
                    end
                    if temp.TransStatus == 1
                        totalHasTrans = totalHasTrans + temp.Population;
                    end
                    totalPopulation = totalPopulation + temp.Population;
                end
            end
            QNbrPower = totalHasPower / totalPopulation;
            QNbrComm = totalHasComm / totalPopulation;
            QNbrTrans = totalHasTrans / totalPopulation;
        end
        
        %% Creat road links from two ending road nodes
        function Trans = assignRoadToRoadNode(Trans,Dictionary)
            Road_Set = Trans{1};
            for i = 1:length(Road_Set)
                roadTemp = Road_Set{i};
                nodeStart = Dictionary(strcat('RoadNode',num2str(roadTemp.Start_Node)));
                nodeEnd = Dictionary(strcat('RoadNode',num2str(roadTemp.End_Node)));
                nodeStart = nodeStart{1};
                nodeEnd = nodeEnd{1};
                nodeStart.Roads{end+1} = roadTemp.uniqueID;
                nodeEnd.Roads{end+1} = roadTemp.uniqueID;
            end
        end

  


  
        
        %% Old Functions
        function return_Schedule = CleanOld(Current, Pow, Comm, Trans, Max_Power, Max_Comm, Max_Trans, system)
            Branch= Pow{1};
            Bus= Pow{2};
            Generator= Pow{3};
            
                        
            Centraloffice = Comm{1};
            CommunicationTower = Comm{2};
            Cellline = Comm{3};

            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};
            
            tmp = Current;
            
            
            for i = 1:length(tmp)
                if ~isempty(tmp{i})
                    if system == 1
                        for j = 1:length(Max_Power)
                            resourceNeed = Library.getResourceOld(tmp{i}, Pow, Comm, Trans);
                            Max_Power(j) = Max_Power(j) + resourceNeed;
                        end
                    elseif system == 2
                        for j = 1:length(Max_Comm)
                            resourceNeed = Library.getResourceOld(tmp{i}, Pow, Comm, Trans);
                            Max_Comm(j) = Max_Comm(j) + resourceNeed;
                        end
                    elseif system == 3
                        for j = 1:length(Max_Trans)
                            resourceNeed = Library.getResourceOld(tmp{i}, Pow, Comm, Trans);
                            Max_Trans(j) = Max_Trans(j) + resourceNeed;
                        end
                    end
                    tem = strsplit(tmp{i},'/');
                    if strcmp(tem{1}, 'Branch') && Branch{str2double(tem{2})}.WorkingDays <= 0
                        Branch{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Bus') && Bus{str2double(tem{2})}.WorkingDays <= 0
                        Bus{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Generator') && Generator{str2double(tem{2})}.WorkingDays <= 0
                        Generator{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Antenna') && Antenna{str2double(tem{2})}.WorkingDays <= 0
                        Antenna{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Centraloffice') && Centraloffice{str2double(tem{2})}.WorkingDays <= 0
                        Centraloffice{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Router') && Router{str2double(tem{2})}.WorkingDays <= 0
                        Router{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Cellline') && Cellline{str2double(tem{2})}.WorkingDays <= 0
                        Cellline{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'CommunicationTower') && CommunicationTower{str2double(tem{2})}.WorkingDays <= 0
                        CommunicationTower{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Road') && Road{str2double(tem{2})}.WorkingDays <= 0
                        Road{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'Bridge') && Bridge{str2double(tem{2})}.WorkingDays <= 0
                        Bridge{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                    elseif strcmp(tem{1}, 'TrafficLight') && TrafficLight{str2double(tem{2})}.WorkingDays <= 0
                        TrafficLight{str2double(tem{2})}.WorkingDays = 0;
                        Current{i} = [];
                        
                        
                    end
                end
            end
            return_Schedule = Current;
        end
        
        function result = getResourceOld(Object, Power_Set, Communication_Set, Transportation_Set)
            Branch= Power_Set{1};
            Bus= Power_Set{2};
            Generator= Power_Set{3};
            
                       
            Centraloffice = Communication_Set{1};
            CommunicationTower = Communication_Set{2};
            Cellline = Communication_Set{3};
            
            
%             Antenna = Communication_Set{1};
%             Centraloffice = Communication_Set{2};
%             Router = Communication_Set{3};
%             Cellline = Communication_Set{4};
%             CommunicationTower = Communication_Set{5};
            
            Road = Transportation_Set{1};
            Bridge = Transportation_Set{2};
            TrafficLight = Transportation_Set{3};
            
            tem = strsplit(Object,'/');
            name = tem{1};
            number = str2num(tem{2});
            result =1;
        end
        function WorkingProcessOld(Current, Pow, Comm, Trans, Days)
            Branch= Pow{1};
            Bus= Pow{2};
            Generator= Pow{3};
            
                        
            Centraloffice = Comm{1};
            CommunicationTower = Comm{2};
            Cellline = Comm{3};
            
            
%             Antenna = Comm{1};
%             Centraloffice = Comm{2};
%             Router = Comm{3};
%             Cellline = Comm{4};
%             CommunicationTower = Comm{5};
            
            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};
            
            for i = 1:length(Current)
                if ~isempty(Current{i})
                    tem = strsplit(Current{i},'/');
                    if strcmp(tem{1}, 'Branch')
                        Branch{str2double(tem{2})}.WorkingDays = Branch{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Branch %d : %d\n', str2double(tem{2}), Branch{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'Bus')
                        if Bus{str2double(tem(2))}.WorkingDays < Days
                        end
                        Bus{str2double(tem{2})}.WorkingDays = Bus{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Bus %d : %d\n', str2double(tem{2}), Bus{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'Generator')
                        Generator{str2double(tem{2})}.WorkingDays = Generator{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Generator %d : %d\n', str2double(tem{2}), Generator{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'Centraloffice')
                        Centraloffice{str2double(tem{2})}.WorkingDays = Centraloffice{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Centraloffice %d : %d\n', str2double(tem{2}), Centraloffice{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'Router')
                        Router{str2double(tem{2})}.WorkingDays = Router{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Router %d : %d\n', str2double(tem{2}), Router{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'Cellline')
                        Cellline{str2double(tem{2})}.WorkingDays = Cellline{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Cellline %d : %d\n', str2double(tem{2}), Cellline{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'CommunicationTower')
                        CommunicationTower{str2double(tem{2})}.WorkingDays = CommunicationTower{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Cellline %d : %d\n', str2double(tem{2}), Cellline{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'Road')
                        Road{str2double(tem{2})}.WorkingDays = Road{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Road %d : %d\n', str2double(tem{2}), Road{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'Bridge')
                        %Bridge{str2double(tem{2})}.WorkingDays = Bridge{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('Bridge %d : %d\n', str2double(tem{2}), Bridge{str2double(tem{2})}.WorkingDays);
                    elseif strcmp(tem{1}, 'TrafficLight')
                        TrafficLight{str2double(tem{2})}.WorkingDays = TrafficLight{str2double(tem{2})}.WorkingDays - Days;
                        %                         fprintf('TrafficLight %d : %d\n', str2double(tem{2}), TrafficLight{str2double(tem{2})}.WorkingDays);
                    end
                end
            end
            
            
            %                         disp('-----------------------------------');
        end
        function graphTheory(Trans_Set)
            hash = containers.Map('KeyType','double','ValueType','double');
            
            input_node = [];
            input_road = [];
            roadnode_Set = Trans_Set{4};
            road_Set = Trans_Set{1};
            for i = 1:length(roadnode_Set)
                input_node{i, 1} = roadnode_Set{i}.Location(1);
                input_node{i, 2} = roadnode_Set{i}.Location(2);
                hash(roadnode_Set{i}.nodeID) = i;
            end
            for i = 1:length(road_Set)
                index1 = hash(road_Set{i}.Start_Location(1));
                index2 = hash(road_Set{i}.End_Location(1));
                input_road{i, 1} = index1;
                input_road{i, 2} = index2;
            end
            
            %grPlot(cell2mat(input_node),cell2mat(input_road),'g');
        end
        function [returnCurrent, returnSchedule] = AddCurrentWorkingOld(Max, Current, Schedule, Dictionary)
            index = 1;
            for i = 1:Max
                if i <= length(Current)
                    if ~isempty(Current{i})
                        continue;
                    end
                end
                found = 0;
                while found == 0 && index <= length(Schedule)
                    tem = strsplit(Schedule{index},'/');
                    if length(tem) == 3
                        tasktemp = Dictionary(tem{3});
                        flag = 0;
                        for j = 1:length(tasktemp.predecessorTask)
                            if(Dictionary(tasktemp.predecessorTask{j}).WorkingDays > 0)
                                flag = 1;
                                break;
                            end
                        end
                        
                        
                        if flag ~= 0
                            index = index + 1;
                            continue;
                        end
                        if tasktemp.WorkingDays == 0
                            disp('Add Task Temp');
                            disp(tasktemp);
                        end
                    end
                    if length(tem) == 2 || length(tem) == 3
                        Schedule{index} = strcat(Schedule{index},'/Working/dummy');
                        Current{i} = Schedule{index};
                        found = 1;
                    end
                    index = index + 1;
                end
            end
            
            returnCurrent = Current;
            returnSchedule = Schedule;
        end
        
                %% Old function of Library.Recovery
        function [Power, Comm, Trans, trackt, trackCW, trackOS] = RecoveryOld(time_horizon, Interdependence_Num, Qtrans0, InterdependenceFunc, ReSchedule_Num,...
                RestorationResource, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, ...
                RepairSchedule, Power_Set, Communication_Set, Transportation_Set, ...
                Dictionary, System_Dependent_Factor, transGraph, powerGraph, commGraph,Neighborhood, ...
                Seperate_Scheduling, LinkDirectionChoice,Per_New_reschedule,Num_stage,...
                TableTask,TablePrecedence, active_power, active_comm, active_trans, OptimizationChoice,Diff_unit,Cust_unit)
            %================================================================
            % function Recovery
            % This function simulates the Recovery Process of all systems.
            % This function calls the following major functions: 
            % 1. in the file of Library: AddCurrentWorking, FindMinDays, 
            % UpdateCurrentDuration, SysFuncInterdependence, UpdateStatus
            % 2. in the file of Interface1: Functionality, RepairReSchedule
            % 
            %================================================================
            %==== Field
            %--- lookupTable: a table that documents all data about all tasks, 
            %such as task description, duration distribution, resource demand, 
            %actually sampled task druation, task uniqueID, etc.
            lookupTable = {};
            
            %--- Set up the initial functionality of every system based on
            %the percentage of population with a certain service at the
            %neighborhood level. 
            [totalPopulation, FunctionalityPower, FunctionalityComm, FunctionalityTrans] = Library.neighbourFunc(Dictionary);
            
            %--- Set up a dictionary named funcTable that documents the
            % functionality of every object and every neighborhood if
            % the object/neighbood is the status of "Damaged" or "Stoped". 
            funcTable = containers.Map('KeyType','char','ValueType','any');
            %funcTable = LibraryDebug.initialFuncTable(Power_Set, Communication_Set, Transportation_Set, funcTable, time_horizon, Neighborhood);
            
            funcTable = Library.initialFuncTable(Power_Set, Communication_Set, Transportation_Set, funcTable, time_horizon, Neighborhood);
%             
            %--- Set up initial variables about the system functionality
            neighbourPowFunc = zeros(1,time_horizon);
            neighbourCommFunc = zeros(1,time_horizon);
            neighbourTransFunc = zeros(1,time_horizon);
            Power = zeros(1,time_horizon);
            Comm = zeros(1,time_horizon);
            Trans = zeros(1,time_horizon);
            
            %TransTest = zeros(4,time_horizon);
            
            %--- Set up the working list for every system
            %This list evolves at different time step. 
            % A task can be added to the list if this task satisfies the
            % constraint of resource and precedence at the moment. 
            CurrentWorking_Power = {};
            CurrentWorking_Comm = {};
            CurrentWorking_Trans = {};
                     
            %--- Set up initial values of supporting variables
            %- finish: a flag variable representing whether the restoration
            % is finished or not. (0-No/1-Yes)
            %- Start_Day: time index representing the starting time of a task
            %- End_Day: time index representing the ending time of a task
            %- flagDelay: the flag variable for tracking the transprotatino delaying effect.
            % initially flagDelay = 0, if the delaying effect is turned on,
            % finish = 0 and Qtrans<Qtrans0, flagDelay = 1; then finish = 0
            % and Qtrans>=Qtrans0, flagDelay = 2; else flagDelay = 0. 
            %- need_reschedule: flag variables representing whether a system
            % needs to be rescheduled. 
            % need_reschedule = [need_reschedulePower, need_rescheduleComm, need_rescheduleTrans]
            finish = 0;
            Start_Day = 1;
            End_Day = 1;
            flagDelay = 0;
            %ploti= 0;
            need_reschedule = zeros(1,3); 
            
            %--- for the purpose of debugging, track the variable
            %evoluations while the time index jumps 
            indtrack = 1;
            %--- trackt: for the purpose of debugging, a vector of time index 
            % to track which time steps has the time index jumps over in 
            % the entire recovery process.
            trackt(indtrack,:) = [Start_Day, End_Day]; 
            %--- trackCW: for the purpose of debugging, track the variable
            %of CurrentWorking while the time index jumps  
            trackCW{indtrack,1} = CurrentWorking_Power;
            trackCW{indtrack,2} = CurrentWorking_Comm;
            trackCW{indtrack,3} = CurrentWorking_Trans;
            %--- trackOS: track Object.Status while the time index jumps 
            trackOS{1}{indtrack,1} = {}; %Power_Set;
            trackOS{2}{indtrack,1} = {}; %Communication_Set;
            trackOS{3}{indtrack,1} = {}; %Transportation_Set; 
     
            %--- set up resource constraints for every system
            Max_Power = RestorationResource(1,:);
            Max_Comm = RestorationResource(2,:); 
            Max_Trans = RestorationResource(3,:);
            
            %--- set up the restoration plan for every system (determined from 
            % Step 3, either scheme 3A-priority or scheme 3B-optimization)
            Schedule_Power = RepairSchedule{1};
            Schedule_Comm = RepairSchedule{2};
            Schedule_Trans = RepairSchedule{3};
            
            %--- set up dictionary of task and dictionary of precedence for
            %every system
            Dic_p_task = TableTask{2,1}; 
            Dic_c_task = TableTask{2,2};
            Dic_t_task = TableTask{2,3};
            Dic_p_prece = TablePrecedence{2,1};
            Dic_c_prece = TablePrecedence{2,2};
            Dic_t_prece = TablePrecedence{2,3};

            %--- Other variables
            % total_damaged: an variable couting the total number of tasks to repair damaged components of all three systems in a sample of damage scenario.
            % total_fixed: an variable couting how many tasks have been fixed. 
            total_damaged = length(Schedule_Power) + length(Schedule_Comm) + length(Schedule_Trans);
            total_fixed = 0;
            
            while finish == 0 % meaning there is (are) still some task(s) left.
                
                %==== Original schedule
                ori_schedule.Schedule_Power=Schedule_Power;
                ori_schedule.Schedule_Comm=Schedule_Comm;
                ori_schedule.Schedule_Trans=Schedule_Trans;
                %==== Check if the variable of Remain_schedule exists
                if ~exist('Remain_schedule','var')
                    Remain_schedule=ori_schedule;
                end
  
                %==== Add tasks for reparing any Damaged Component to the 
                % Current Working List Based on the Resource Constraints: 
                % Max_Power, Max_Comm, Max_Trans
                
                %--- time index for debuging
                func_index = Start_Day;
                %ts = Start_Day 
                %te = End_Day
                indtrack = indtrack
                CurrentTime = Start_Day
                RestorationPhase = Cust_unit; % a horizontal vector of different time instants in the time unit in "day".
                
                if Diff_unit %== 1 
                    RestorationTimePhase = sort([CurrentTime, RestorationPhase]);
                        [~,n]=find(RestorationTimePhase == CurrentTime); n=4;% n = either one of (1, 2, 3, 4)

                        if gt(n,length(RestorationPhase)+1) % n>4 if length(RestorationPhase) = 3 [e.g., RestorationPhase = [3,28,168]];
                            msg = strcat('Function Recovery Error: Time Index of CurrentTime (n) = ', num2str(n), ' > length(RestorationPhase)+1 for ', tasktemp);
                            disp(msg); 
                            return               
                        end
                        
                    % Diff_unit is not 1, meaning the choice of using different time unit in recovery is turned off.  
                    % All tasks using the same time unit of "day"
                    else
                        n = 2;
                    end
                
                %Rbefore1 = [Max_Power; Max_Comm; Max_Trans];
                [CurrentWorking_Power, Schedule_Power, Max_Power, lookupTable] = Library.AddCurrentWorking(Max_Power, CurrentWorking_Power, Schedule_Power, Dictionary, lookupTable, Start_Day, n, indtrack);
                [CurrentWorking_Comm, Schedule_Comm, Max_Comm, lookupTable] = Library.AddCurrentWorking(Max_Comm, CurrentWorking_Comm, Schedule_Comm, Dictionary, lookupTable, Start_Day, n, indtrack);
                [CurrentWorking_Trans, Schedule_Trans, Max_Trans, lookupTable] = Library.AddCurrentWorking(Max_Trans, CurrentWorking_Trans, Schedule_Trans, Dictionary, lookupTable, Start_Day, n, indtrack);
  
                %Rafter1 = [Max_Power; Max_Comm; Max_Trans];
                %dR1 = Rafter1-Rbefore1
                
%                 % Debugging-related to WorkingDays 
                
%                 [DurationP{indtrack}, RdemandP{indtrack}] = LibraryDebug.GetTaskDuration(Dictionary, CurrentWorking_Power);
%                 [DurationC{indtrack}, RdemandC{indtrack}] = LibraryDebug.GetTaskDuration(Dictionary, CurrentWorking_Comm);
%                 [DurationT{indtrack}, RdemandT{indtrack}] = LibraryDebug.GetTaskDuration(Dictionary, CurrentWorking_Trans);
                
                
                %==== Find the shortest duration that a task takes among all tasks in the list of CurrentWorking
                % Find Day (the shorest time of a task takes to
                % complete among all tasks in the list of CurrentWorking)
                % so that the time index can jump from the current time
                % step to the new time step, when this task with the shortest 
                % duration in the CurrentWorking list finishes. 
                % Days in the unit of day
               Days = Library.FindMinDays(CurrentWorking_Power, CurrentWorking_Comm, CurrentWorking_Trans, Power_Set, Communication_Set, Transportation_Set,Dictionary)
                
                
                %%%%%%%% Revised Expression
                %[Days, TimeUnitUnq] = LibraryDebug.FindMinDuration2(DurationP{indtrack}, DurationC{indtrack}, DurationT{indtrack}); 

                
%                 if Days ~= 0 % original expression
                if Days > 0 % revised expression 
                    Start_Day = End_Day;    
%                   %--- End_Day: revised expression
                    End_Day = End_Day + Days;
                    
                    indtrack = indtrack+1;
                    trackt(indtrack,:) = [Start_Day, End_Day];
                    trackCW{indtrack,1} = CurrentWorking_Power;
                    trackCW{indtrack,2} = CurrentWorking_Comm;
                    trackCW{indtrack,3} = CurrentWorking_Trans;
                    
%                     %--- End_Day: original expression
%                     if Days == 1
%                         End_Day = End_Day + 1;
%                     else
%                         End_Day = End_Day + ceil(Days);
%                     end
  
                else
                    End_Day = time_horizon;
                    finish = 1;
                    
                    indtrack = indtrack+1;
                    trackt(indtrack,:) = [Start_Day, End_Day];
                    trackCW{indtrack,1} = CurrentWorking_Power;
                    trackCW{indtrack,2} = CurrentWorking_Comm;
                    trackCW{indtrack,3} = CurrentWorking_Trans;
                end
                
                %==== Updating the WorkingDays in the list of CurrentWorking
                % by subtracting the original Object.WorkingDays by "Day" from 
                % the above line.                 
%                 Dbefore = [Dictionary('Task1').WorkingDays;
%                 Dictionary('Task35').WorkingDays;
%                 Dictionary('Task97').WorkingDays;
%                 Dictionary('Task106').WorkingDays;
%                 Dictionary('Task1134').WorkingDays;
%                 Dictionary('Task698').WorkingDays];
%                 Rbefore2 = [Max_Power; Max_Comm; Max_Trans]
%%%%%%%%%%%%Revised Expresions that combining both duration and resource                  
%                 [CurrentWorking_Power, Max_Power, TaskComplete_Power] = LibraryDebug.UpdateCurrent(Dictionary, CurrentWorking_Power, Max_Power, Days);
%                 [CurrentWorking_Comm, Max_Comm, TaskComplete_Comm] = LibraryDebug.UpdateCurrent(Dictionary, CurrentWorking_Comm, Max_Comm, Days);
%                 [CurrentWorking_Trans, Max_Trans, TaskComplete_Trans] = LibraryDebug.UpdateCurrent(Dictionary, CurrentWorking_Trans, Max_Trans, Days);
%                 
%%%%%%%%%%%%Original Expresions                
                Library.UpdateCurrentDuration(CurrentWorking_Power, Dictionary, Days);
                Library.UpdateCurrentDuration(CurrentWorking_Comm, Dictionary, Days);
                Library.UpdateCurrentDuration(CurrentWorking_Trans, Dictionary, Days);
                
%             
%%%%%%%%%%%%Original Expresions                
                %==== Update the list of CurrentWorking for every system by
                %cleaning the task(s) that have been executed, and by
                %adding the amount of used resource due to executing task 
                %in the current list to Max_Power
                [CurrentWorking_Power,Max_Power] = Library.UpdateCurrentResource(CurrentWorking_Power, Dictionary, Max_Power);
                [CurrentWorking_Comm,Max_Comm] = Library.UpdateCurrentResource(CurrentWorking_Comm, Dictionary, Max_Comm);
                [CurrentWorking_Trans,Max_Trans] = Library.UpdateCurrentResource(CurrentWorking_Trans, Dictionary, Max_Trans);

                     
                %==== Update Object.Status and Functionality for repaired
                %Component/Object in every system
                
                trackOS{1}{indtrack,1} = {}; %Power_Set;
                trackOS{2}{indtrack,1} = {}; %Communication_Set;
                trackOS{3}{indtrack,1} = {}; %Transportation_Set; 

 
% %%%%%%%%%%%%%%%%%%%%%%%Revised Expression                 
%                 [Remain_schedule,need_reschedule, total_fixed,transGraph, powerGraph, commGraph,...
%                         Power_Set, Communication_Set, Transportation_Set, funcTable,FunctionalityTotal,...
%                          totalPopulation, trackOS] = LibraryDebug.UpdateStatusNew(Power_Set, Communication_Set, Transportation_Set,...
%                          total_damaged, total_fixed, need_reschedule,Dictionary,...
%                          transGraph, powerGraph, commGraph,funcTable,Start_Day,End_Day,...
%                          Per_New_reschedule,Num_stage,Remain_schedule,ori_schedule,...
%                          Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirectionChoice, indtrack, trackOS);
                     

%%%%%%%%%%%%%%%%%%%%%%Original Expression                
                [Remain_schedule,need_reschedule, total_fixed,transGraph, powerGraph, commGraph,...
                        Power_Set, Communication_Set, Transportation_Set, funcTable,FunctionalityTotal,...
                         totalPopulation, trackOS] = Library.UpdateStatusNew(Power_Set, Communication_Set, Transportation_Set,...
                         total_damaged, total_fixed, need_reschedule,Dictionary,...
                         transGraph, powerGraph, commGraph,funcTable,Start_Day,End_Day,...
                         Per_New_reschedule,Num_stage,Remain_schedule,ori_schedule,...
                         Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirectionChoice, indtrack, trackOS);
%                      

                
                % InterdependenceFunc: turn on(1)/off(0) inter-system functionality dependency 
%                 if InterdependenceFunc == 1 
%                     [Remain_schedule,need_reschedule, total_fixed,transGraph, powerGraph, commGraph,...
%                         Power_Set, Communication_Set, Transportation_Set, funcTable,FunctionalityTotal,...
%                          totalPopulation] = Library.UpdateStatusNew(Power_Set, Communication_Set, Transportation_Set,...
%                          total_damaged, total_fixed, need_reschedule,Dictionary,...
%                          transGraph, powerGraph, commGraph,funcTable,Start_Day,End_Day,...
%                          Per_New_reschedule,Num_stage,Remain_schedule,ori_schedule,...
%                          Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirectionChoice);
%                  
%                 elseif  InterdependenceFunc == 0
%                     [Remain_schedule,need_reschedule, total_fixed,transGraph, powerGraph, commGraph,...
%                         Power_Set, Communication_Set, Transportation_Set, funcTable,FunctionalityTotal,...
%                          totalPopulation] = Library.UpdateStatusSep(Power_Set, Communication_Set, Transportation_Set,...
%                          total_damaged, total_fixed, need_reschedule,Dictionary,...
%                          transGraph, powerGraph, commGraph,funcTable,Start_Day,End_Day,...
%                          Per_New_reschedule,Num_stage,Remain_schedule,ori_schedule,...
%                          Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirectionChoice);
%                 end
                
                %==== Update System functionality and Neighborhood functionality
                FunctionalityPower = FunctionalityTotal(1);
                FunctionalityComm = FunctionalityTotal(2);
                FunctionalityTrans = FunctionalityTotal(3);
                FuncPowerNbr = FunctionalityTotal(4);
                FuncCommNbr = FunctionalityTotal(5);
                FuncTransNbr = FunctionalityTotal(6);

%%%%%%%%%%%%%%%%%%%% Revised Expressions
                Power(indtrack)=FunctionalityPower;
                Comm(indtrack) = FunctionalityComm;
                Trans(indtrack) = FunctionalityTrans;
                Qneighbour{1}(indtrack) = FuncPowerNbr;
                Qneighbour{2}(indtrack) = FuncCommNbr;
                Qneighbour{3}(indtrack) = FuncTransNbr;
                
                %ts = [ts,Start_Day]
                %te = [te,End_Day]
%%%%%%%%%%%%%%%%%%%%Original Expressions
%                 Power(Start_Day:End_Day) = FunctionalityPower;
%                 Comm(Start_Day:End_Day) = FunctionalityComm;
%                 Trans(Start_Day:End_Day) = FunctionalityTrans;
%                 
%                 neighbourPowFunc(Start_Day:End_Day) = FuncPowerNbr;
%                 neighbourCommFunc(Start_Day:End_Day) = FuncCommNbr;
%                 neighbourTransFunc(Start_Day:End_Day) = FuncTransNbr;

                %==== Check transport delay effect due to the functionality
                % disruptions of the transportation system to decide
                % whether the duration of remainign tasks in the power
                % andcommunication systems need to be adjusted by
                % multipling a delaying factor (i.e., variable System_Dependent_Factor). 
                % If Interdependence_Num == 1, the delaying effet is turned
                % on. Simply multiplying the delaying factor to the duration sample of
                % remaining tasks for the power and communication systems,
                % when Qtrans(t)<Qtrans0 (transportation functionality threshold). 
                % If Interdependence_Num == 0, the delaying effet is turned off. 
                % Simply leave all the duration samples as it is for every 
                % remaining task for the power and communication systems. 
                if Interdependence_Num == 1 % delaying effect is turned on.
                    %FunctionalityTrans = FunctionalityTotal(3);
                    if flagDelay == 0 && lt(FunctionalityTrans, Qtrans0)
                        Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, System_Dependent_Factor, Dictionary);
                        flagDelay = 1;
                    elseif flagDelay == 1 && ge(FunctionalityTrans, Qtrans0)
                        Restore_Factor = 1 / System_Dependent_Factor;
                        Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, Restore_Factor, Dictionary);
                        flagDelay = 2;
                    else
                        Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, 1, Dictionary);
                        
                    end
                else   % Interdependence_Num == 0, delaying effect is turned off. 
                    Library.SysFuncInterdependence(Power_Set, Communication_Set, Transportation_Set, 1, Dictionary);
                end
                
                if ~(ReSchedule_Num==0)
                    %==== Check Reschedule effect              
                    % need_reschedule = [need_reschedulePower, need_rescheduleComm, need_rescheduleTrans]
                    need_reschedulePower = need_reschedule(1);
                    need_rescheduleComm = need_reschedule(2);
                    need_rescheduleTrans = need_reschedule(3);
                    if ~exist('Per_New_reschedule_Power','var')
                        Per_New_reschedule_Power=Per_New_reschedule;
                        Per_New_reschedule_Comm=Per_New_reschedule;
                        Per_New_reschedule_Trans=Per_New_reschedule;
                        Num_stage_Power=Num_stage;
                        Num_stage_Comm=Num_stage;
                        Num_stage_Trans=Num_stage;
                    end
                    [need_reschedulePower,Per_New_reschedule_Power,Num_stage_Power] = Library.checkNeedRescheduleNew(need_reschedulePower, Power, Start_Day, End_Day,Per_New_reschedule_Power,Num_stage_Power);
                    [need_rescheduleComm,Per_New_reschedule_Comm,Num_stage_Comm] = Library.checkNeedRescheduleNew(need_rescheduleComm, Comm, Start_Day, End_Day,Per_New_reschedule_Comm,Num_stage_Comm);            
                    [need_rescheduleTrans,Per_New_reschedule_Trans,Num_stage_Trans] = Library.checkNeedRescheduleNew(need_rescheduleTrans, Trans, Start_Day, End_Day, Per_New_reschedule_Trans,Num_stage_Trans);
                    need_reschedule = [need_reschedulePower, need_rescheduleComm, need_rescheduleTrans];


                    checkReschedule = find(need_reschedule == 1); % find which system has need_reschedule=1
                    if ~isempty(checkReschedule)

                        % Get cells of remaining tasks and the precedence of
                        % remaining tasks
                        [remain_task_p,remain_task_c,remain_task_t,remain_prece_p,remain_prece_c,remain_prece_t]=Library.RemainTaskPrecedence(Dic_p_task,Dic_c_task,Dic_t_task,Dic_p_prece,Dic_c_prece,Dic_t_prece,Remain_schedule);
                        Remain_task={remain_task_p,remain_task_c,remain_task_t};
                        Remain_precedence={remain_prece_p,remain_prece_c,remain_prece_t};

                        resource = [Max_Power; Max_Comm; Max_Trans];

                        if need_reschedule(1) == 1 && active_power && ~isempty(remain_task_p)
                            disp('Reschedule Power')
                            System = 'Power'; 
                            [Power_ReSchedule, Power_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, Remain_task, Remain_precedence, resource, time_horizon);
                            Schedule_Power = Power_ReSchedule;
                        end
                        if need_reschedule(2) == 1 && active_comm && ~isempty(remain_task_c)
                            disp('Reschedule Communication')
                            System = 'Communication';
                            [Comm_ReSchedule, Comm_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, Remain_task, Remain_precedence, resource, time_horizon);
                            Schedule_Comm = Comm_ReSchedule;
                        end
                        if need_reschedule(3) == 1 && active_trans && ~isempty(remain_task_t)
                            disp('Reschedule Transportation')
                            System = 'Transportation';
                            [Trans_ReSchedule, Trans_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, Remain_task, Remain_precedence, resource, time_horizon);
                            Schedule_Trans = Trans_ReSchedule;
                        end               
                    end
                end
                need_reschedule=zeros(1,3);
%                 %==== record current schedule as original schedule in ori_schedule
%                 ori_schedule.Schedule_Power=Schedule_Power;
%                 ori_schedule.Schedule_Comm=Schedule_Comm;
%                 ori_schedule.Schedule_Trans=Schedule_Trans;
                
                %p = plot(commGraph,'Layout','force');
                %saveas(p,strcat('./test/Comm', num2str(ploti),'.jpg'));
                %p = plot(powerGraph,'Layout','force');
                %saveas(p,strcat('./test/Power', num2str(ploti),'.jpg'));
                %TransTest = Library.ComputeFunc(transGraph, Start_Day, End_Day,TransTest);
                
%                 if active_trans~=0
%                     TransTest = Library.Functionality_GraphBasic(transGraph, LinkDirectionChoice, Start_Day, End_Day, TransTest);
%                 end 
               %ploti = ploti + 1;
                Start_Day = End_Day;
               %Start_Day = End_Day + 1;
            end % while finish == 0
            
            %==== Error Check
%             filename = strcat('lookupTable.mat');
%             save(filename, 'lookupTable');
%             printedFuncTable = Library.printFuncTable(funcTable, time_horizon);
%             filename = strcat('printedFuncTable.mat');
%             save(filename, 'printedFuncTable');
            
            
            if Power(time_horizon) < 1
                disp('Function Recovery Warning: Functionality Power is less than 100% at t = time_horizon.');
            end
           
            if Comm(time_horizon) < 1
                disp('Function Recovery Warning: Functionality Communication is less than 100% at t = time_horizon.');
            end
            
            if Trans(time_horizon) < 1
                disp('Function Recovery Warning: Functionality Transportation is less than 100% at t = time_horizon.');
            end
            
        end
               
        
        
        %% older version of function UpdateStatusNew
        function [Remain_schedule,need_reschedule, total_fixed,  transGraph, powerGraph, commGraph, Pow, Comm, Trans, funcTable, FunctionalityTotal, totalPopulation, trackOS] = UpdateStatusNewOld(Pow, Comm, Trans, total_damaged, total_fixed, need_reschedule, Dictionary, transGraph, powerGraph, commGraph, funcTable, Start_Day, End_Day, Per_New_reschedule, Num_stage, Remain_schedule,ori_schedule, Power_Func_Num, Trans_Func_Num, Comm_Func_Num, LinkDirectionChoice, indtrack, trackOS)
            % ================================================================
            % function UpdateStatusNew
            % called by function Library.Recovery
            % This function updates the Object.Status after finishing a
            % task as needed. The object status changes from "Samaged" to "Stoped"/"Open",
            % or the object status changes from "Stoped" to "Open".
            % =================================================================
            
            %==== Field
            Branch= Pow{1};
            Bus= Pow{2};
            Generator= Pow{3};
            TransTower = Pow{4};
            
            Centraloffice = Comm{1};
            CommunicationTower = Comm{2};
            Cellline = Comm{3};            
            
            Road = Trans{1};
            Bridge = Trans{2};
            TrafficLight = Trans{3};
            
            fixed = false;
            index_task=0;
            
            need_reschedulePower = need_reschedule(1);
            need_rescheduleComm = need_reschedule(2);
            need_rescheduleTrans = need_reschedule(3);
            
            %==== Power System
            %--- Generator
            for ii = 1:length(Generator)              
                if strcmp(Generator{ii}.Status, 'Damaged') && ~isempty(Generator{ii}.WorkingDays)           
                    if Generator{ii}.WorkingDays <= 0
                        Generator{ii}.Status = 'Open';
                        
                        temp = funcTable(Generator{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Generator{ii}.uniqueID) = temp;
                      
                        Generator{ii}.Functionality = 1;
                        total_fixed = total_fixed + 1;
                        index_task = index_task+1;
                        
                        %--- track which Object changes status 
                        trackOS{1}{indtrack,end+1} = strcat('Generator',num2str(ii), ', from Damaged to Open');
                        
                        
                        %--- Rescheduling: 
                        % I have quesitons about these lines. Are the following lines wrong?!!
                        [New_Schedule_Power{1:size(ori_schedule.Schedule_Power,2)-ii}] = ori_schedule.Schedule_Power{ii+1:end};
                        Remain_schedule.Schedule_Power = New_Schedule_Power;
                        New_Schedule_Power=[];
                        
                    end
                end
            end
            
            %--- Bus (Substation)
            for ii = 1:length(Bus)
                
                %-- For substations (bus) that are connected to a power plant (generator)
                if ~isempty(Bus{ii}.Generator)
                    %-- If a substation is damaged and its task(s) finishes
                    % the restoration (workingdays <=0), change the status
                    % from "Damaged" to "Stoped". 
                    if strcmp(Bus{ii}.Status, 'Damaged') && ~isempty(Bus{ii}.WorkingDays)
                        Bus = Library.getWorkDays(Bus, ii,Dictionary);
                        if  Bus{ii}.WorkingDays <= 0
                            Bus{ii}.Status = 'Stoped';
                            total_fixed = total_fixed + 1;
                            
                            %--- track which Object changes status 
                            trackOS{1}{indtrack,end+1} = strcat('Bus',num2str(ii), ', from Damaged to Stoped');
                        
                        end
                    end
                    
                    %-- If a substation is stoped and its connected generator
                    % (power plant) is open now, change the status of this
                    % substation from "Stoped" to "Open", and update its Functionality. 
                    if strcmp(Bus{ii}.Status, 'Stoped')
                        temp = extractAfter(Bus{ii}.Generator, 9);
                        temp = str2num(temp);
                        if strcmp(Generator{temp}.Status, 'Open')
                            Bus{ii}.Status = 'Open';
                            fixed = true;
                            temp = funcTable(Bus{ii}.uniqueID);
                            temp(End_Day:end) = 1;
                            funcTable(Bus{ii}.uniqueID) = temp;
                            
                            Bus{ii}.Functionality = 1;
                            powerGraph = addnode(powerGraph,Bus{ii}.uniqueID);
                            
                            %--- track which Object changes status 
                            trackOS{1}{indtrack,end+1} = strcat('Bus',num2str(ii), ', from Stoped to Open');
                        
                            
                        end
                    end
                    
                %-- For substations (bus) that are not connected to a power plant (generator)    
                % If a substation is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Open". 
                else
                    if strcmp(Bus{ii}.Status, 'Damaged') && ~isempty(Bus{ii}.WorkingDays)
                        if  Bus{ii}.WorkingDays <= 0
                            Bus{ii}.Status = 'Open';
                            fixed = true;
                            Bus{ii}.Functionality = 1;
                            
                            temp = funcTable(Bus{ii}.uniqueID);
                            temp(End_Day:end) = 1;
                            funcTable(Bus{ii}.uniqueID) = temp;
                            
                            powerGraph = addnode(powerGraph,Bus{ii}.uniqueID);
                            total_fixed = total_fixed + 1;
                            
                            %--- track which Object changes status 
                            trackOS{1}{indtrack,end+1} = strcat('Bus',num2str(ii), ', from Damaged to Open');
                        
 
%                             need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        end
                    end
                end
                
                %-- If a substation is fixed (fixed=true, "Open"), update the 
                % power functionality property of the neighborhood (Neighborhood{j}.PowerStatus) 
                % that connects to this restored/fixed substation.  
                if fixed
                    for j = 1:length(Bus{ii}.Neighborhood)
                        temp = Dictionary(Bus{ii}.Neighborhood{j});
                        temp = temp{1};
                        disp( temp.uniqueID);
                        temp.PowerStatus = 1;
                        temp = Dictionary(Bus{ii}.Neighborhood_Power_Link{j});
                        temp = temp{1};
                        temp.Status = 'Open';
                        
                        %--- track which Object changes status 
                        trackOS{1}{indtrack,end+1} = strcat(temp.uniqueID, ': Nbr.Bus, from Stoped to Open');
                        
                    end
                    % prepare the flag variable of fixed for the computation of the next Bus object by setting this flag variable to the initial value of "false". 
                    fixed = false;
                end
            end
            
            %--- Transmission Tower
            % If a transmission tower is damaged and its task(s) finishes
            % the restoration (workingdays <=0), change the status
            % from "Damaged" to "Open", and set the Object.Functionality= 1.
            for ii = 1:length(TransTower)
                if strcmp(TransTower{ii}.Status, 'Damaged') && ~isempty(TransTower{ii}.WorkingDays)
                    if  TransTower{ii}.WorkingDays <= 0
                        TransTower{ii}.Status = 'Open';
                        TransTower{ii}.Functionality = 1;
                        temp = funcTable(TransTower{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(TransTower{ii}.uniqueID) = temp;                      
                        powerGraph = addnode(powerGraph,TransTower{ii}.uniqueID);
                        total_fixed = total_fixed + 1;   
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('TransTower',num2str(ii), ', from Damaged to Open');
                        
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end
                end
            end
            
            %--- Branch (Power Line) 
            for ii = 1:length(Branch)
                %-- If a power line is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Stoped". 
                if strcmp(Branch{ii}.Status, 'Damaged') && ~isempty(Branch{ii}.WorkingDays)
                    Branch = Library.getWorkDays(Branch, ii,Dictionary);             
                    if Branch{ii}.WorkingDays <= 0
                        Branch{ii}.Status = 'Stoped';
                        total_fixed = total_fixed + 1;
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('Branch',num2str(ii), ', from Damaged to Stoped');
                        
                        
                        % the following lines about rescheduling is
                        % confusing??
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        index_task=index_task+1;
                        if index_task+1> size(Remain_schedule.Schedule_Power,2)
                            New_Schedule_Power=[];
                        else
                        [New_Schedule_Power{1:size(Remain_schedule.Schedule_Power,2)-index_task}] = Remain_schedule.Schedule_Power{index_task+1:end};
                        end
                        Remain_schedule.Schedule_Power=New_Schedule_Power;New_Schedule_Power=[];
                    end
                end
                
                %-- If a power line is fixed but in the "Stoped" status, 
                % as long as the two end objects are "Open", change the status
                % from "Stoped" to "Open", and set the Object.Functionality
                % = 1. 
                if strcmp(Branch{ii}.Status, 'Stoped')
                    obj1 = Dictionary(Branch{ii}.connectedObj1);
                    obj2 = Dictionary(Branch{ii}.connectedObj2);
                    if strcmp(obj1{1}.Status, 'Open')&&strcmp(obj2{1}.Status, 'Open')
                        Branch{ii}.Status = 'Open';
                        Branch{ii}.Functionality = 1;
                        temp = funcTable(Branch{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Branch{ii}.uniqueID) = temp;
                        powerGraph = Library.addPowerGraph(Pow, powerGraph, ii, Dictionary);     
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('Branch',num2str(ii), ', from Stoped to Open');
                        
                    end
                end
            end
            

            %==== Communication
            %--- Central Office
            for ii = 1:length(Centraloffice)
                
                %-- If a central office is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Stoped".
                if strcmp(Centraloffice{ii}.Status, 'Damaged') && ~isempty(Centraloffice{ii}.WorkingDays)
                    Centraloffice = Library.getWorkDays(Centraloffice, ii, Dictionary);             
                    if Centraloffice{ii}.WorkingDays <= 0
                        Centraloffice{ii}.Status = 'Stoped';
                        total_fixed = total_fixed + 1;
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('Central Office',num2str(ii), ', from Damaged to Stoped');
                        
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end
                end
                %-- If a central office is in the "Stoped" status, and there is a magic battery for this central office, change the status
                % from "Stoped" to "Open", and set Object.Functionality = 1.           
                if strcmp(Centraloffice{ii}.Status, 'Stoped')
                    if Centraloffice{ii}.Battery == 1
                        Centraloffice{ii}.Status = 'Open';
                        fixed = true;
                        temp = funcTable(Centraloffice{i}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Centraloffice{ii}.uniqueID) = temp;
                        Centraloffice{ii}.Functionality = 1;
                        commGraph = addnode(commGraph,Centraloffice{ii}.uniqueID);
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('Bus',num2str(ii), ', from Stoped to Open');
                        
                        continue;
                    end
                    bus = Dictionary(Centraloffice{ii}.Bus);
                    if strcmp(bus{1}.Status, 'Open')
                        Centraloffice{ii}.Status = 'Open';
                        fixed = true;
                        temp = funcTable(Centraloffice{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Centraloffice{ii}.uniqueID) = temp;
                        Centraloffice{ii}.Functionality = 1;
                        commGraph = addnode(commGraph,Centraloffice{ii}.uniqueID);
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('Bus',num2str(ii), ', from Stoped to Open');
                    end
                end
                
                %-- If a central office is fixed (fixed=true, "Open"), update the 
                % communication functionality property of the neighborhood (Neighborhood{j}.CommStatus) 
                % that connects to this restored/fixed central office.  
                if fixed
                    for j = 1:length(Centraloffice{ii}.Neighborhood)
                        temp = Dictionary(Centraloffice{ii}.Neighborhood{j});
                        temp = temp{1};
                        disp(temp.uniqueID);
                        temp.CommStatus = 1;
                        temp = Dictionary(Centraloffice{ii}.Neighborhood_Comm_Link{j});
                        temp = temp{1};
                        temp.Status = 'Open';
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat(temp.uniqueID,': Nbr.CO, from Stoped to Open');
                    end
                    
                    fixed = false;
                end
            end
            
            %--- CommunicationTower
            for ii = 1:length(CommunicationTower)
                %-- If a communication tower is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Stoped". If 
                if strcmp(CommunicationTower{ii}.Status, 'Damaged') && ~isempty(CommunicationTower{ii}.WorkingDays)
                    CommunicationTower = Library.getWorkDays(CommunicationTower, ii, Dictionary);
                    if CommunicationTower{ii}.WorkingDays <= 0
%                         % original statement, seems wrong
                        if size(Bus,1)~=0 
                            CommunicationTower{ii}.Status = 'Stoped';
                            
                            %--- track which Object changes status 
                            trackOS{2}{indtrack,end+1} = strcat('CT',num2str(ii), ', from Damaged to Stoped');
                        else
                            CommunicationTower{ii}.Status = 'Open';
                            
                            %--- track which Object changes status 
                            trackOS{2}{indtrack,end+1} = strcat('CT',num2str(ii), ', from Damaged to Open');
                        end
                        
%                         % revised statement to consider the magic battery  
%                         bus = Dictionary(CommunicationTower{ii}.Bus);
%                         if strcmp(bus{1}.Status, 'Open')
%                             CommunicationTower{ii}.Status = 'Open';
%                         elseif CommunicationTower{ii}.Battery == 1
%                             CommunicationTower{ii}.Status = 'Open';
%                         else
%                             CommunicationTower{ii}.Status = 'Stoped';
%                         end
                        total_fixed = total_fixed + 1;
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        
                        if ii+1> size(ori_schedule.Schedule_Comm,2)
                            New_Schedule_Comm=[];
                        else
                        [New_Schedule_Comm{1:size(ori_schedule.Schedule_Comm,2)-ii}] = ori_schedule.Schedule_Comm{ii+1:end};
                        end
                        Remain_schedule.Schedule_Comm=New_Schedule_Comm;New_Schedule_Comm=[];
                    end
                end
                
                if size(Bus,1)~=0
                    if strcmp(CommunicationTower{ii}.Status, 'Stoped')
                        if CommunicationTower{ii}.Battery == 1
                            CommunicationTower{ii}.Status = 'Open';
                            temp = funcTable(CommunicationTower{ii}.uniqueID);
                            temp(End_Day:end) = 1;
                            funcTable(CommunicationTower{ii}.uniqueID) = temp;
                            CommunicationTower{ii}.Functionality = 1;
                            commGraph = addnode(commGraph,CommunicationTower{ii}.uniqueID);
                            
                            %--- track which Object changes status 
                            trackOS{2}{indtrack,end+1} = strcat('CT',num2str(ii), ', from Stoped to Open');
                            
                            continue;
                        end

                        bus = Dictionary(CommunicationTower{ii}.Bus);
                        if strcmp(bus{1}.Status, 'Open')
                            CommunicationTower{ii}.Status = 'Open';
                            temp = funcTable(CommunicationTower{ii}.uniqueID);
                            temp(End_Day:end) = 1;
                            funcTable(CommunicationTower{ii}.uniqueID) = temp;
                            CommunicationTower{ii}.Functionality = 1;
                            commGraph = addnode(commGraph,CommunicationTower{ii}.uniqueID);
                            
                            %--- track which Object changes status 
                            trackOS{2}{indtrack,end+1} = strcat('CT', num2str(ii), ', from Stoped to Open');
                        end
                    end
                end
            end
            
            %--- Cell Line
            for ii = 1:length(Cellline)
                %-- If a cell line is damaged and its task(s) finishes
                % the restoration (workingdays <=0), change the status
                % from "Damaged" to "Stoped". 
                if strcmp(Cellline{ii}.Status, 'Damaged') && ~isempty(Cellline{ii}.WorkingDays)
                    Cellline = Library.getWorkDays(Cellline, ii,Dictionary);
                    if Cellline{ii}.WorkingDays <= 0
                        Cellline{ii}.Status = 'Stoped';
                        total_fixed = total_fixed + 1;
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('CL',num2str(ii), ', from Damaged to Stoped');
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end
                end
                if strcmp(Cellline{ii}.Status, 'Stoped')        
                    bus = Dictionary(Cellline{ii}.Bus);
                    obj1 = Dictionary(Cellline{ii}.connectedObj1);
                    obj2 = Dictionary(Cellline{ii}.connectedObj2);
                    
                    if strcmp(bus{1}.Status, 'Open')&& strcmp(obj1{1}.Status, 'Open')&&strcmp(obj2{1}.Status, 'Open')
                        Cellline{ii}.Status = 'Open';
                        temp = funcTable(Cellline{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Cellline{ii}.uniqueID) = temp;
                        Cellline{ii}.Functionality = 1;
                        commGraph = Library.addCommGraph(Comm, commGraph, ii,Dictionary);
                        
                        %--- track which Object changes status 
                        trackOS{2}{indtrack,end+1} = strcat('CL',num2str(ii), ', from Stoped to Open');
                        
                    end
                end
            end
              
            %==== Transportation
            %--- TrafficLight
            for ii = 1:length(TrafficLight)
                if strcmp(TrafficLight{ii}.Status, 'Damaged') && ~isempty(TrafficLight{ii}.WorkingDays)
                    TrafficLight = Library.getWorkDays(TrafficLight, ii,Dictionary);
                    if TrafficLight{ii}.WorkingDays <= 0
                        TrafficLight{ii}.Status = 'Stoped';
                        total_fixed = total_fixed + 1;
                        
                        %--- track which Object changes status 
                        trackOS{3}{indtrack,end+1} = strcat('TL',num2str(ii), ', from Damaged to Stoped');
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end
                end
                
                if strcmp(TrafficLight{ii}.Status, 'Stoped')
                    if TrafficLight{ii}.Battery == 1
                        TrafficLight{ii}.Status = 'Open';
                        temp = funcTable(TrafficLight{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(TrafficLight{ii}.uniqueID) = temp;
                        
                        TrafficLight{ii}.Functionality = 1;
                        transGraph = addnode(transGraph,TrafficLight{ii}.uniqueID);
                        
                        %--- track which Object changes status 
                        trackOS{3}{indtrack,end+1} = strcat('TL',num2str(ii), ', from Stoped to Open');
                        
                        continue;
                    end
                    temp = extractAfter(TrafficLight{ii}.Bus, 3);
                    temp = str2num(temp);
                    if ~isempty(temp) && strcmp(Bus{temp}.Status, 'Open')
                        TrafficLight{ii}.Status = 'Open';
                        temp = funcTable(TrafficLight{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(TrafficLight{ii}.uniqueID) = temp;
                        TrafficLight{ii}.Functionality = 1;
                        
                        %--- track which Object changes status 
                        trackOS{3}{indtrack,end+1} = strcat('TL',num2str(ii), ', from Stoped to Open');
                    end
                end
            end
            
            %--- Bridge
            for ii = 1:length(Bridge)
                flag = 1;
                if strcmp(Bridge{ii}.Status, 'Damaged') && ~isempty(Bridge{ii}.WorkingDays) && Bridge{ii}.HasSub == 0 
                    tasks = Bridge{ii}.taskUniqueIds;
                    sumWorkDay = 0;
                    for j = 1:length(tasks)
                        temp = Dictionary(tasks{j});
                        if iscell(temp)
                            temp = temp{1};
                        end
                        if temp.WorkingDays > 0
                            sumWorkDay = sumWorkDay + temp.WorkingDays;
                            flag = 0;
                        end
                    end
                    Bridge{ii}.WorkingDays = sumWorkDay;
                    if Bridge{ii}.WorkingDays(1) <= 0 || flag
                        Bridge{ii}.Status = 'Open';
                        temp = funcTable(Bridge{ii}.uniqueID);
                        temp(Start_Day:end) = 1;
                        funcTable(Bridge{ii}.uniqueID) = temp;
                        Bridge{ii}.Functionality = 1;
                        total_fixed = total_fixed + 1;
                        
                        %--- track which Object changes status 
                        trackOS{3}{indtrack,end+1} = strcat('Bridge',num2str(ii), ', from Damaged to Open');
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end
                end
                % for subcomponent
                if strcmp(Bridge{ii}.Status, 'Damaged') && ~isempty(Bridge{ii}.WorkingDays) && Bridge{ii}.HasSub == 1
                    for sub_index = 1:length(Bridge{ii}.ColumnSet)
                        tasks = Bridge{ii}.ColumnSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.ColumnFoundSet)
                        tasks = Bridge{ii}.ColumnFoundSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.AbutmentSet)
                        tasks = Bridge{ii}.AbutmentSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.AbutmentFoundSet)
                        tasks = Bridge{ii}.AbutmentFoundSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.BearingSet)
                        tasks = Bridge{ii}.BearingSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                    
                    for sub_index = 1:length(Bridge{ii}.SlabSet)
                        tasks = Bridge{ii}.SlabSet(sub_index).taskUniqueIds;
                        sumWorkDay = 0;
                        for j = 1:length(tasks)
                            temp = Dictionary(tasks{j});
                            if iscell(temp)
                                temp = temp{1};
                            end
                            if temp.WorkingDays > 0
                                sumWorkDay = sumWorkDay + temp.WorkingDays;
                                flag = 0;
                            end
                        end
                         Bridge{ii}.WorkingDays = sumWorkDay;
                    end
                   
                    if Bridge{ii}.WorkingDays(1) <= 0 || flag
                        Bridge{ii}.Status = 'Open';
                        temp = funcTable(Bridge{ii}.uniqueID);
                        temp(Start_Day:end) = 1;
                        funcTable(Bridge{ii}.uniqueID) = temp;
                        Bridge{ii}.Functionality = 1;
                        total_fixed = total_fixed + 1;
                        need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                    end                   
                end

            end
            
            %--- Road
            for ii = 1:length(Road)
                if strcmp(Road{ii}.Status, 'Damaged') &&  ~isempty(Road{ii}.WorkingDays)
                    Road = Library.getWorkDays(Road, ii,Dictionary);
                    if Road{ii}.WorkingDays <= 0
                        total_fixed = total_fixed + 1;
%                         need_reschedule = Library.checkNeedReschedule(total_damaged, total_fixed, need_reschedule,Per_New_reschedule);
                        Road{ii}.Status = 'Stoped';
                    end
                end
                if strcmp(Road{ii}.Status, 'Stoped')
                    flag = 0;
                    for j = 1:length(Road{ii}.Bridge_Carr)
                        if strcmp(Bridge{Road{ii}.Bridge_Carr(j)}.Status, 'Damaged')
                            flag = 1;
                            break;
                        end
                    end
                    for j = 1:length(Road{ii}.Bridge_Cross)
                        if strcmp(Bridge{Road{ii}.Bridge_Cross(j)}.Status, 'Damaged')
                            flag = 1;
                            break;
                        end
                    end
                    
                    if flag == 0
                        Road{ii}.Status = 'Open';
                        temp = funcTable(Road{ii}.uniqueID);
                        temp(End_Day:end) = 1;
                        funcTable(Road{ii}.uniqueID) = temp;
                        Road{ii}.Functionality = 1;
                        transGraph = Library.addTransGraph(Trans, transGraph, ii);
                        
                        %--- track which Object changes status 
                        trackOS{3}{indtrack,end+1} = strcat('Road',num2str(ii), ', from Stoped to Open');
                        
                        temp = Dictionary(strcat('RoadNode',num2str(Road{ii}.Start_Node)));
                        temp = temp{1};
                        for j = 1:length(temp.Neighborhood)
                            t1 = Dictionary(temp.Neighborhood{j});
                            t1 = t1{1};
                            t1.TransStatus = 1;
                            t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                            t1 = t1{1};
                            t1.Status = 'Open';
                            
                            %--- track which Object changes status 
                            trackOS{3}{indtrack,end+1} = strcat(t1.uniqueID, ', from Stoped to Open');
                        end
                        temp = Dictionary(strcat('RoadNode',num2str(Road{ii}.End_Node)));
                        temp = temp{1};
                        for j = 1:length(temp.Neighborhood)
                            t1 = Dictionary(temp.Neighborhood{j});
                            t1 = t1{1};
                            t1.TransStatus = 1;
                            t1 = Dictionary(temp.Neighborhood_Trans_Link{j});
                            t1 = t1{1};
                            t1.Status = 'Open';
                            
                            %--- track which Object changes status 
                            trackOS{3}{indtrack,end+1} = strcat(t1.uniqueID, ': Nbr.Road, from Stoped to Open');
                        end
                    end
                end
            end
            
            %==== Update Object status to System Set
            Pow{1} = Branch;
            Pow{2} = Bus;
            Pow{3} = Generator;
            Pow{4} = TransTower;
            
            Comm{1} = Centraloffice;
            Comm{2} = CommunicationTower;
            Comm{3} = Cellline;
            
            Trans{1} = Road;
            Trans{2} = Bridge;
            Trans{3} = TrafficLight;
            
            %==== Calculated System Functionality
            % Do we need reschedule? If yes, need_reschedule=1. 
            [FunctionalityPower, FunctionalityComm, FunctionalityTrans] = Interface1.Functionality(Power_Func_Num, Trans_Func_Num, Comm_Func_Num, Pow, Comm, Trans, powerGraph, commGraph, transGraph, Dictionary,LinkDirectionChoice);
            [totalPopulation,FuncPowerNbr,FuncCommNbr,FuncTransNbr] = Library.neighbourFunc(Dictionary);% Propotion of population that has power/comm/trans
            FunctionalityTotal = [FunctionalityPower, FunctionalityComm, FunctionalityTrans, FuncPowerNbr,FuncCommNbr,FuncTransNbr];
            
%             need_reschedulePower = Library.checkNeedRescheduleNew(need_reschedulePower, FunctionalityPower, Start_Day, End_Day, Per_New_reschedule,Num_stage);
%             need_rescheduleComm = Library.checkNeedRescheduleNew(need_rescheduleComm, FunctionalityPower, Start_Day, End_Day, Per_New_reschedule,Num_stage);            
%             need_rescheduleTrans = Library.checkNeedRescheduleNew(need_rescheduleTrans, FunctionalityPower, Start_Day, End_Day, Per_New_reschedule,Num_stage);
%             need_reschedule = [need_reschedulePower,need_rescheduleComm, need_rescheduleTrans];
%             
        end
        

        %% Test Functions
        function Count(Set)
            count = 0;
            for i = 1:length(Set)
                for j = 1:length(Set{i})
                    if ~strcmp(Set{i}{j}.Status, 'Open')
                        disp(Set{i}{j});
                        count = count + 1;
                    end
                end
            end
            disp(count);
        end
        function CountTotal(Set1,Set2,Set3)
            Set = {Set1,Set2,Set3};
            for k = 1:length(Set)
                count = 0;
                for i = 1:length(Set{k})
                    for j = 1:length(Set{k}{i})
                        if ~strcmp(Set{k}{i}{j}.Status, 'Open')
                            count = count + 1;
                        end
                    end
                end
                disp(count);
            end
        end
        function CountZeros(Current, Dictionary)
            for i = 1:length(Current)
                if ~isempty(Current{i})
                    temp = Dictionary(Library.getUniqueId(Current{i}, 1));
                    if iscell(temp)
                        temp = temp{1};
                    end
                    if temp.WorkingDays <= 0
                        disp('WARNING....');
                        disp(Current{i});
                        disp(temp);
                    end
                end
            end
        end
        function WrapperCountZeros(a,b,c,Dictionary,num)
            disp(num);
            Library.CountZeros(a,Dictionary);
            Library.CountZeros(b,Dictionary);
            Library.CountZeros(c,Dictionary);
        end
        
        
        function    [ st, msg ] = cmd_rmdir( folderspec )       
            %   cmd_rmdir removes a directory and its contents 
            %   
            %   Removes all directories and files in the specified directory in
            %   addition to the directory itself.  Used to remove a directory tree.
            %   See also: xtests\generic_utilies_test.m

                narginchk( 1, 1 )

                dos_cmd = sprintf( 'rmdir /S /Q "%s"', folderspec );

                [ st, msg ] = system( dos_cmd );
        end

        
                %% Abort the computation when the Stop button is clicked on the GUI
        % Once the Stop botton is clicked, the global variable "stopVariable" changes its value from "false" to "true". 
        function AbortComputation(stopVariable) 
            
            if strcmp(stopVariable, 'off') 
                %%%finish;
                % The following lines are the new command lines that I have
                % wrote based on "finishsav.m"
                %==========================================================================
                disp(getString(message('MATLAB:finishsav:SavingWorkspaceData')));
                % savefilename = 'WorkspaceVariables.mat';
                save WorkspaceVariables.mat;
                exit;
            end        
        end
        
        %% Function for progressbar
        % https://www.mathworks.com/matlabcentral/fileexchange/3607-progressbar
        %gui_active.m
        %progressbar.m
        function handle = progressbar( handle,increment,string,titlestr )
            %
            % progressbar - shows a progress bar dialog based on the function "waitbar"
            %
            % Format: handle = progressbar( handle,increment [,string,titlestr] )
            %
            % Input:    handle      - handle to current progress bar, [] for a new one
            %           increment   - a fraction (0..1) to increment by.
            %                         (-1) signals the function to remove the handle from
            %                         the persistant list and close the progressbar
            %           string      - a string to be replaced in the progress bar (optional)
            %           titlestr    - a title for the dialog box (optional)
            %
            % Output:   handle      - a graphic handle of the dialog box
            %
            %
            % NOTE:     this function uses a persistant list of handles and thier values.
            %           therefore,  to delete a progressbar, please call the function with: 
            %               progressbar( handle,-1 );
            %
            %           an "abort" button is formed inside the progressbar, if the calling process
            %           uses the persistent function "gui_active". when the "abort" button is pressed,
            %           the Callback function "gui_active" changes it's value to zero, which enables 
            %           the application to interactively stop execution
            %
            % Example:  gui_active(1);      % will add an abort button
            %           h           = progressbar( [],0,'my test' );
            %           max_count   = 1e+3;
            %           for idx = 1:max_count
            %               fprintf( '%d\n',idx )';
            %               h = progressbar( h,1/max_count );
            %               if ~gui_active
            %                   break;
            %               end
            %           end
            %           progressbar( h,-1 );
            %


            persistent handle_list counter_list last_handle last_idx;

            % initialize
            % =============
            call_flag = min( nargin,4 );

            % analyze input and decide what to do
            % ====================================
            if isempty( handle )            % create a new dialog
                counter_list(end+1) = 0;
                last_idx            = length( counter_list );
                switch call_flag
                case 2, last_handle = waitbar( increment,'Please Wait...' );
                case 3, last_handle = waitbar( increment,string );
                case 4, last_handle = waitbar( increment,string,'Name',titlestr );
                end        
                handle_list(end+1)  = last_handle;
                handle              = last_handle;
                %check_position( handle_list );      % so that the figures don't hide each other
                if (gui_active)
                    Library.add_button( last_handle );      % add the abort button if the state of the gui_active is set
                end

            elseif ( increment == -1 )      % delete correct handle from the list
                last_handle             = handle;
                last_idx                = find( handle_list == handle );
                handle_list( last_idx ) = [];
                counter_list( last_idx )= [];
                if ishandle( last_handle )      % close the figure, if it's open
                    close( last_handle );       % since user can close it by him self
                end
                last_handle             = [];

            elseif (handle == last_handle)  % update last dialog
                counter_list(last_idx)  = counter_list(last_idx) + increment;
                if ishandle( handle )       % nobody killed my figure
                    switch call_flag
                    case 2, waitbar( counter_list(last_idx),handle );
                    case 3, waitbar( counter_list(last_idx),handle,string );
                    case 4, waitbar( counter_list(last_idx),handle,string,'Name',titlestr );
                    end
                else                        % somebody killed my figure -> so I create it again
                    switch call_flag
                    case 2, handle = waitbar( counter_list(last_idx),'Please Wait...' );
                    case 3, handle = waitbar( counter_list(last_idx),string );
                    case 4, handle = waitbar( counter_list(last_idx),string,'Name',titlestr );
                    end
                    handle_list(last_idx)   = handle;
                    last_handle             = handle;
                    Library.check_position( handle_list );      % so that the figures don't hide each other
                    if (gui_active)
                        add_button( last_handle );      % add the abort button if the state of the gui_active is set
                    end
                end    
            else                            % find the handle inside the list
                last_handle = handle;
                last_idx    = find( handle_list == handle );
                if ~isempty( last_idx )
                    counter_list(last_idx)  = counter_list(last_idx) + increment;
                    switch call_flag
                    case 2, waitbar( counter_list(last_idx),last_handle );
                    case 3, waitbar( counter_list(last_idx),last_handle,string );
                    case 4, waitbar( counter_list(last_idx),last_handle,string,'Name',titlestr );
                    end        
                end
            end

            % update display after all
            % ==========================
            drawnow;

        end
% =======================================================================================
%                               Inner Function Implementation
% =======================================================================================

        function add_button( fig_handle )
            %
            % adds the abort button to the waitbar
            %

            % collect handles and set control units to pixels
            axes_handle     = get( fig_handle,'currentaxes' );
            last_fig_units  = get( fig_handle,'units' );
            last_axes_units = get( axes_handle,'units' );
            set( fig_handle,'units','pixels' );
            set( axes_handle,'units','pixels' );

            % read controls position
            fig_position    = get( fig_handle,'position' );
            axes_position   = get( axes_handle,'position' );
            fig_width       = fig_position(3);
            fig_height      = fig_position(4);
            axes_xcoord     = axes_position(1);
            axes_ycoord     = axes_position(2);
            axes_width      = axes_position(3);
            axes_height     = axes_position(4);

            % load the button icon and create the button
            load( 'gauge_abort_icon' );
            button_width    = ButtonSize16x16(1)+2;
            button_height   = ButtonSize16x16(2)+2;
            button_margin   = 10;
            button_xcoord   = (fig_width + axes_width + axes_xcoord - button_width)/2 - button_margin;
            button_ycoord   = (axes_height - button_height)/2 + axes_ycoord;
            button_handle   = uicontrol( 'Parent',fig_handle,'units','pixels',...
                'Position',[ button_xcoord,button_ycoord,button_width,button_height ],...
                'Callback','gui_active(0);progressbar(get(gcbo,''parent''),-1);close(get(gcbo,''parent''));',...
                'CData',Icon16x16 );

            % resize axis to accommodate the button, and restore axes and figure units back to previous
            axes_position(3) = axes_width - button_width - button_margin;
            set( axes_handle,'position',axes_position );
            set( fig_handle,'units',last_fig_units );
            set( axes_handle,'units',last_axes_units );
        end

% ---------------------------------------------------------------------------------------

         function check_position( handle_list )
            %
            % makes sure that the progressbar does not hide it's nested progressbars
            %

            % process only if there is more than one progress bar on the screen
            if (length(handle_list)>1)
                y_increment = 70;                   % pixels
                x_increment = 30;                   % pixels

                % change units to pixels
                screen_units    = get( 0,'units' );
                last_fig_units  = get( handle_list(end-1),'units' );
                cur_fig_units   = get( handle_list(end),'units' );
                set( 0,'units','pixels' );
                set( handle_list(end-1),'units','pixels' );
                set( handle_list(end),'units','pixels' );

                % get positions, and calc new position for progress bar
                screen_size     = get( 0,'screensize' );
                last_position   = get( handle_list(end-1),'position' );
                cur_position    = get( handle_list(end),'position' );
                cur_position(1) = last_position(1) + x_increment;
                cur_position(2) = last_position(2) - y_increment;

                % check that we don't go outside the screen
                if (cur_position(1)+cur_position(3)>screen_size(3))
                    cur_position(1) = x_increment;
                end
                if (cur_position(2)<screen_size(1))
                    cur_position(2) = screen_size(4) - y_increment - cur_position(4);
                end       

                % store new position and restore units
                set( handle_list(end),'position',cur_position );
                set( 0,'units',screen_units );
                set( handle_list(end-1),'units',last_fig_units );
                set( handle_list(end),'units',cur_fig_units );
            end
         end
        
        function ret = gui_active(optional_input)
            %
            % gui_active - used to implement an abort function in a GUI
            %
            % initially : gui_active(1)
            %
            % in the aplication, occaisionally 
            %  EITHER  1) if (~gui_active) then ... do abort action
            %  OR      2) polling gui_active(-1) will cause an error with the value 'abort_requested'
            %             that can be caught with a "try-catch" block
            %
            % to initiate the abort (for both cases) call:  gui_active(0)
            %

            % the call to drawnow enable other controls to "patch-in" before the
            % rest of this code execute and possibly change the "is_active" state

            persistent    is_active;

            if (nargin > 0)
                if (optional_input == -1) 
                    drawnow;
                    if (~is_active)
                        error('abort_requested');
                    end
                else
                    prev_is_active  = is_active;
                    is_active       = optional_input;
                    if ((prev_is_active > 0) & (is_active == 0))
                        disp('operation aborted');
                    end
                end
            else
                drawnow;
            end

            ret = is_active;
        end
    end
end
