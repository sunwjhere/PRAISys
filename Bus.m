classdef Bus<SystemGeneral & handle & matlab.mixin.Copyable
    properties (SetAccess=public)
        Location=[];
        Name;
        Generator;
        Centraloffice;
        Capacity;
        Line_In;
        Branch = [];
        Road;
        Neighborhood = [];
        Neighborhood_Power_Link = [];
        taskUniqueIds = [];
        PopulationServed = 0;
        HouseholdServed = 0;
    end
    
    
    methods
        %Constructor initializes proerty values with input arguments
%         function Bran=Branch(Number, Location, Priority)
%             Bran.Number=Number;
%             Bran.Location=Location;
%             Bran.Priority=Priority;
%         end
        
        function obj = Bus(varargin)
            if nargin == 8 || nargin == 7
                args{1} = varargin{1};
                args{2} = varargin{3};
                args{3} = 'Bus';
            else
            end
            obj = obj@SystemGeneral(args{:});
            obj.Capacity = varargin{2};
            obj.Type = varargin{3};% Substation classification based on HAZUS: 'ESS1';'ESS3';'ESS5'; etc.
            obj.Name = varargin{4};
            obj.Location = varargin{5};
                       
            obj.PopulationServed = varargin{6};
            obj.HouseholdServed = varargin{7};
            if nargin == 8
                obj.Generator = varargin{8};
            end                        
            
            if nargin == 6
                obj.Generator = varargin{6};
            end
        end
        
        %Display Selected information about the account
        function getStatement(Bus)
            disp('----------------------------')
            disp(['Branch:', num2str(Bus.Number)])
            disp(['CurrentStatus:', Bus.Status])
            disp(['Recovery:', num2str(Bus.Recovery)])
            disp(['DamageLevel:', num2str(Bus.DamageLevel)])
            disp('----------------------------')
        end
        
        function addGenerator(Bus, Generator)
            Bus.Generator=Generator;
        end
        
        function addBus(bus, Bus)
            bus.Bus=[bus.Bus, Bus];
        end
        
    end
end