classdef Centraloffice<SystemGeneral & handle & matlab.mixin.Copyable
    %COMMUNICATION_TOWER  
    %   Detailed explanation goes here
    
     properties (SetAccess=public)
        Location=[];
        Bus;
        Router;
        Road;
        Cellline = {};
        Label;
        Code;
        Company;
        CommTower = {};
        CentralOffice = {};
        Neighborhood_Comm_Link = [];
        Neighborhood = [];
        taskUniqueIds = [];
        Battery = 0;
        PopulationServed = 0;
        HouseholdServed = 0;
    end
    
    methods
%           function Cell = Centraloffice(Number, Location, Antenna,Priority)
%             Cell.Number=Number;
%             Cell.Location=Location;
%             Cell.Antenna=Antenna;
%             Cell.Priority=Priority;
%           end


          
          function obj = Centraloffice(Number, Location)
              if nargin == 2
                  args{1} = Number;
                  args{2} = 'C1L';
                  args{3} = 'CentralOffice';
              else
              end
              obj = obj@SystemGeneral(args{:});
              obj.Location=Location;
              obj.Code='High-code'; %'High-code', 'Moderate-code', 'Low-code', 'Pre-code'
          end
          
          function getStatement(CO)
              disp('----------------------------')
              disp(['Central Office:', num2str(CO.Number)])
              disp(['Status:', CO.Status])
              disp(['Type:', CO.Type])
              disp(['Location:', num2str(CO.Location)]);
              %disp(['Antenna:', num2str(CO.Antenna)]);
              %disp(['Router:', num2str(CO.Router)]);
              disp('----------------------------')
          end
          
          function addRouter(CO,Router)
              CO.Router=[CO.Router,Router];
          end
    end
    
end