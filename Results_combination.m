Save_history=1;
if exist('PATH_dir','var')
    Hostname = 'Output/Results0'; 
    if Save_history~=0
            files=dir(fullfile(PATH_dir,Hostname(1:6)));
            all_re = {files([files(:).isdir]).name}; all_re(ismember(all_re,{'.','..'})) = [];
    end
    Hostname=fullfile(PATH_dir,Hostname);
else
    Hostname = 'Output/Results0'; 
    if Save_history~=0
            files=dir(Hostname(1:6));
            all_re = {files([files(:).isdir]).name}; all_re(ismember(all_re,{'.','..'})) = [];
    end    
end    

Hostname=['Output/',all_re{1}];
load(['./',Hostname,'/mat/Functionality.mat'],'QuSys');
tem_QuSys=QuSys;

for i=2:size(all_re,2)
    Hostname=['Output/',all_re{i}];
    load(['./',Hostname,'/mat/Functionality.mat'],'QuSys');
    for j=1:size(QuSys,2)
        for kk=1:size(QuSys{1,j},2)
            tem_QuSys{1,j}{1,kk}=[tem_QuSys{1,j}{1,kk};QuSys{1,j}{1,kk}];
        end
    end
end
fname='Output';
figoption = 'on'; 

PlotFigureFunctionality(tem_QuSys, fname, figoption)
%% Plot figures of system functionality samples and means
        function PlotFigureFunctionality(FunctionalitySystem, fname, figoption)
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
            
            Nsample = size(Functionality_Power{1,1},1);
            if Nsample>=2
                nrow = 2;
            end
            
            %==== Power

            subplot(nrow,ncol,1),
            for ii=1:Nsample
                hold on, stairs(Functionality_Power{imtr}(ii,:)); 
            end
            ylim([0,1]);
            mf=mean(Functionality_Power{imtr});
            stairs(mf(1,:), 'LineWidth', lw);            
            box on; grid on; 
            xlabel('Time (day)');
            ylabel('Q_{power}');
            title('Power Functionality Sample');

            %==== Communication
            subplot(nrow,ncol,2),
            for ii=1:Nsample
                hold on, stairs(Functionality_Communication{imtr}(ii,:));
                io(ii,:)=Functionality_Communication{imtr}(ii,:);
            end
            ylim([0,1]);
            mf=mean(Functionality_Communication{imtr});
            stairs(mf(1,:), 'LineWidth', lw);
            box on; grid on; 
            xlabel('Time (day)');
            ylabel('Q_{comm}');
            title('Communication Functionality Sample');

            
            N = 12;                                      % Number of ?Experiments? In Data Set

yMean = mf;                                    % Mean Of All Experiments At Each Value Of ?x?
ySEM = std(io)/sqrt(N);                              % Compute ?Standard Error Of The Mean? Of All Experiments At Each Value Of ?x?
CI95 = tinv([0.025 0.975], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution
yCI95 = bsxfun(@times, ySEM, CI95(:));              % Calculate 95% Confidence Intervals Of All Experiments At Each Value Of ?x?                                      % Plot Mean Of All Experiments
hold on
kk= yCI95+yMean; 
stairs(kk(1,:), 'LineWidth', lw); % Plot 95% Confidence Intervals Of All Experiments
stairs(kk(2,:), 'LineWidth', lw); % Plot 95% Confidence Intervals Of All Experiments
hold off
grid


            %==== Transportation 
            subplot(nrow,ncol,3),
            for ii=1:Nsample
                hold on, stairs(Functionality_Transportation{imtr}(ii,:));
            end
            ylim([0,1]);
            mf=mean(Functionality_Transportation{imtr});
            stairs(mf(1,:), 'LineWidth', lw);            
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
            
            saveas(fig,strcat( deblank(fname), '/Qsystem_all.jpg'));

        end