classdef CommunicationTower<SystemGeneral & handle & matlab.mixin.Copyable
    %COMMUNICATION_TOWER  
    %   Detailed explanation goes here
    
     properties (SetAccess=public)
        Location=[];
        Bus;
        Road;
        Line;
        Centraloffice = {};
        Cellline = {};
        StructuralType;
        Height;
        %Name;
        City;
        Owner;
        %Usage;
        %BackupPower;
        %Comment;
        CommunicationTowerID;
        %Tract;
        taskUniqueIds = [];
        Battery = 0;
        PopulationServed = 0;
        HouseholdServed = 0;
    end
    
    methods
        function getStatement(CT)
            disp('----------------------------')
            disp(['CommunicationTower:', num2str(CT.Number)])
            disp(['Status:', CT.Status])
            disp(['Type:', CT.Type])
            disp(['Location:', num2str(CT.Location)])
            disp(['WorkingDays:', num2str(CT.WorkingDays)])
            disp('----------------------------')
        end
        
%         function Cell=CommunicationTower(Number,Location, Branch,Priority)
%             Cell.Number=Number;
%             Cell.Location=Location;
%             Cell.Branch=Branch;
%             Cell.Priority=Priority;
%         end
        
        function obj=CommunicationTower(Number,Location)
            if nargin == 2
                args{1} = Number;
                args{2} = 'CTP';%'CommunicationTower_SteelLattice';
                args{3} = 'CommunicationTower';
            else
            end
            obj = obj@SystemGeneral(args{:});
            obj.Location=Location;
        end
        
        function addCentraloffice(CommunicationTower,Centraloffice)
            CommunicationTower.Centraloffice=[CommunicationTower.Centraloffice, Centraloffice];
          end
    end
    
end

