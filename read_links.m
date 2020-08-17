fpath = 'Input';
fn = {'RoadLink','TrafficLights','Bridges','RoadNode'};
% myfolder = strcat(pwd,'\',fpath);
myfolder = strcat(pwd,'/',fpath);
mydir = fullfile(myfolder,fn{1});
table1 = readtable([mydir,'.csv']);



mydir = fullfile(myfolder,fn{3});
table3 = readtable([mydir,'.csv']); 

% table1 = readtable('RoadLink.csv'); 
% table2 = readtable('TrafficLights.csv');
% table3 = readtable('Bridges.csv');


%=== Bridge
all_Bridges_links=table3{:,8:15};

for k=2:4
 for i=1:size(table1,1)
[m,n]=find(table1{i,1}==all_Bridges_links(:,k));
  if k==2 
      bridgeIDcarry{i}=m;
  else
      bridgeIDcarry{i}=[bridgeIDcarry{i};m];
  end 
 end
end

for k=6:8
 for i=1:size(table1,1)
[m,n]=find(table1{i,1}==all_Bridges_links(:,k));
  if k==2 
      bridgeIDcross{i}=m;
  else
      bridgeIDcross{i}=[bridgeIDcarry{i};m];
  end 
 end
end


%=== Traffic Light

mydir = fullfile(myfolder,fn{2});
if exist([mydir,'.csv'])
    active_traLights = true;
    table2 = readtable([mydir,'.csv']);
    mm=[table2{:,1},table2{:,6}];
    for i=1:size(table1,1)
        [m,n]=find(table1{i,1}==mm(:,2));
        trafficlightID{i}=m; 
    end

else
    trafficlightID=[];
    active_traLights = false;
end





save newroadlinks.mat bridgeIDcarry bridgeIDcross trafficlightID




