classdef Interface1
    methods(Static)
        % Select Scheduler Model
%         function [Schedule, Date] = RepairSchedule(num, time_horizon, Max_Power, Max_Trans, Max_Comm, Power_Priority, Comm_Priority, Trans_Priority,...
%                 active_power, active_comm, active_trans, Power_Set, Communication_Set, Transportation_Set,...
%                 priority_power_num, priority_transportation_num, priority_communication_num, ...
%                 OptimizationChoice, taskTable, precedenceTable)
%             
%             taskP = taskTable{1};
%             taskC = taskTable{2};
%             taskT = taskTable{3};
%             switch num
%                 case 1
%                     if active_power && ~isempty(taskP)
%                         switch priority_power_num
%                             case 1
%                                 [result_pow, pow_date] = Library.PowerSchedulePriority(1,Power_Set);
%                             case 2
%                                 [result_pow, pow_date] = Library.PowerSchedulePriority(2,Power_Set);
%                         end
%                     else
%                         result_pow = [];
%                         pow_date = [];
%                     end
%                     
%                     if active_trans && ~isempty(taskT)
%                         switch priority_transportation_num
%                             case 1 % ranking creteria: bridge traffic
%                                 [result_trans, trans_date] = Library.TransSchedulePriority(1,Transportation_Set);
%                             case 2 % ranking creteria: bridge length
%                                 [result_trans, trans_date] = Library.TransSchedulePriority(2,Transportation_Set);
%                         end
%                     else
%                         result_trans = [];
%                         trans_date = [];
%                     end
%                     
%                     if active_comm && ~isempty(taskC)
%                         switch priority_communication_num
%                             case 1
%                                 [result_comm, comm_date] = Library.CommSchedulePriority(1,Communication_Set);
%                             case 2
%                                 [result_comm, comm_date] = Library.CommSchedulePriority(2,Communication_Set);
%                         end
%                     else
%                         result_comm = [];
%                         comm_date = [];
%                     end
% 
%                     Schedule = {result_pow, result_comm, result_trans};
%                     Date = {pow_date, comm_date, trans_date};
%                     
%                 case 2
%                     resource = [Max_Power; Max_Trans; Max_Comm];
%                     if active_power
%                         [result_pow, pow_date] = Library.RepairScheduleOptimization(OptimizationChoice, 'Power', taskTable, precedenceTable, resource, time_horizon);
% 
%                     else
%                         result_pow = [];
%                         pow_date = [];
%                     end
%                     
%                     if active_trans
%                         [result_trans, trans_date] = Library.RepairScheduleOptimization(OptimizationChoice, 'Transportation', taskTable, precedenceTable, resource, time_horizon);
% 
%                     else
%                         result_trans = [];
%                         trans_date = [];
%                     end
%                     
%                     if active_comm
%                         [result_comm, comm_date] = Library.RepairScheduleOptimization(OptimizationChoice, 'Communication', taskTable, precedenceTable, resource, time_horizon);
% 
%                     else
%                         result_comm = [];
%                         comm_date = [];
%                     end
%                     Schedule = {result_pow, result_comm, result_trans};
%                     Date = {pow_date, comm_date, trans_date};
%             end
%         end
        
        
        %% ================================================================
        % function RepairReSchedule
        % Select Re-Scheduler Model
        %  ================================================================
        function Return_Schedule = RepairReSchedule(RescheduleChoice, OptimizationChoice, remainTask, remainPrecedenceTable, Max_Power, Max_Trans, Max_Comm, TimeHorizon, active_power, active_comm, active_trans)
            resource = [Max_Power; Max_Comm; Max_Trans];
            switch RescheduleChoice
                case 1
                    if active_power
                        System = 'Power'; 
                        [Power_ReSchedule, Power_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, remainTask, remainPrecedenceTable, resource, TimeHorizon)
                    else
                        Power_ReSchedule = [];
                        Power_Date = [];
                    end
                    if active_comm
                        System = 'Communication';
                        [Comm_ReSchedule, Comm_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, remainTask, remainPrecedenceTable, resource, TimeHorizon)
                    else
                        Comm_ReSchedule = [];
                        Comm_Date = [];
                    end
                    
                    if active_trans
                        System = 'Transportation';
                        [Trans_ReSchedule, Trans_Date] = Library.RepairScheduleOptimization(OptimizationChoice, System, remainTask, remainPrecedenceTable, resource, TimeHorizon)
                    else
                        Trans_ReSchedule = [];
                        Trans_Date = [];
                    end
                 
                
            end
            
            Return_Schedule = {Power_ReSchedule, Comm_ReSchedule, Trans_ReSchedule};
            
        end
            
            
        
        function Schedule = RepairReScheduleOld(num, orig_Schedule, time, Max_Power, Max_Trans, Max_Comm, Power_Set, Communication_Set, Transportation_Set)
            switch num
                case 1
                    Schedule = Library.RepairScheduleReScheduleMean(orig_Schedule, time, Max_Power, Max_Trans, Max_Comm, Power_Set, Communication_Set, Transportation_Set);   
                case 2
                    Schedule = Library.RepairScheduleReScheduleActual(orig_Schedule, time, Max_Power, Max_Trans, Max_Comm, Power_Set, Communication_Set, Transportation_Set);      
            end
        end
        
        % Interdependence Selector
        function InterdependenceFactor(Interdependence_Num, CurrentWorking_Power, CurrentWorking_Comm, CurrentWorking_Trans, Power_Set, Communication_Set, Transportation_Set)
            switch Interdependence_Num
                case 1
                    Library.InterdenpendenceFactorCalculateBasic(CurrentWorking_Power, CurrentWorking_Comm, CurrentWorking_Trans, Power_Set, Communication_Set, Transportation_Set);
            end
        end
        
        % Create the 1/0 precedenceTable as input for scheduler
        function PreTableSep = createPreTableSep(Dictionary, active_power, active_comm, active_trans)
            if active_power
                result_pow = Library.createPreTableIndividual(Dictionary, 'Power');
            else
                result_pow = [];
            end
            
            if active_trans
                result_trans = Library.createPreTableIndividual(Dictionary,'Transportation');
            else
                result_trans = [];
            end
            
            if active_comm
                result_comm = Library.createPreTableIndividual(Dictionary,'Communication');
            else
                result_comm = [];
            end
            PreTableSep = {result_pow, result_comm, result_trans};
        end
        
        % Create the input table that contains all task information for
        % optimial scheduler
        function taskTable = createTaskTableSep(Dictionary, active_power, active_comm, active_trans)
            if active_power
                result_pow = Library.createTaskTableIndividual(Dictionary, 'Power');
            else
                result_pow = [];
            end
            
            if active_trans
                result_trans = Library.createTaskTableIndividual(Dictionary,'Transportation');
            else
                result_trans = [];
            end
            
            if active_comm
                result_comm = Library.createTaskTableIndividual(Dictionary,'Communication');
            else
                result_comm = [];
            end
            taskTable = {result_pow, result_comm, result_trans};
        end
        

      

        
        
    end
end