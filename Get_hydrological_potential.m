function get_scheme_locations = get_locations(Area,P_min,P_max,folder_location);

%% Read file containing the points along the streams

location = strcat(folder_location,'\Inputs_perc.csv');

Data=readtable(location);
DataArray=Data{:,:};

DataArray=sortrows(DataArray,[2,1],{'ascend' 'ascend'});

%% Constants

ro=1000; %ro is the water density in kg/m^3
eta=0.85; %eta is the turbine efficiency (no units)
g=9.8; %g is the acceleration due to gravity in m/s^2

n_flow=14; %the column containing the initial Q40 flow

Flow_perc{1,1}='Q1';


for i=5:5:95
    Flow_perc{i,1}=strcat('Q',num2str(i));
end

Flow_perc=Flow_perc(~any(cellfun('isempty',Flow_perc),2),:);

%% Define river segments

%Get unique values which correspond to river segments
counts=unique(DataArray(:,2));
River=cell(length(counts),1);
TempRiver=[];

%Put each river segment in a separate cell

j=1;
lastValue=DataArray(1,2);

for i=1:(height(Data))
    if i<height(Data)
    newValue=DataArray(i+1,2);
    else
        if i==height(Data)
        newValue=DataArray(i,2);
        end
    end
        if newValue==lastValue
           TempRiver(i,:)=DataArray(i,:);
        else
               TempRiver(i,:)=DataArray(i,:);
               TempRiver(all(TempRiver==0,2),:)=[];
               River{j,1}(:,:,:)=TempRiver(:,:);
               TempRiver=[]; 
         j=j+1;
         lastValue=newValue;
        end
end

TempRiver(all(TempRiver==0,2),:)=[];
River{j,1}(:,:,:)=TempRiver(:,:);

%Check if points along the river have the same elevation, if yes replace
%them with "NaN" so the algorithm can calculate power

for i=1:length(River)
   curCell=River{i};
   for j=1:size(curCell,1)-1
       if curCell(j,4)==curCell(j+1,4)
           River{i}(j,:)=NaN;
       end
   end
end

%% Get power at each combination of points along the river

TempP=cell(length(counts),1);
PointID=cell(length(counts),1);
Pmax=cell(length(counts),1);
Pmax2=cell(length(counts),2);

%Start with Q40 and substract Q95

for i=1:length(River)   
    for j=1:size(River{i},1)
        for k=j:size(River{i},1)            
            TempP{i}(k,j)= (g*ro*eta*(River{i}(j,n_flow)-River{i}(j,3))*abs(River{i}(j,25)-River{i}(k,25)))/1000; %Power in kW
            PointID{i}(k,1)=River{i}(k,1);
        end           
    end 
end

P=[PointID TempP];

%Calculate the maximum available power and add it to a database

Intake=zeros(length(Pmax2),6);
PowerHouse=zeros(length(Pmax2),6);

for i=1:length(Pmax2)
    Pmax(i,1)=P(i,1);
    for j=1:size(River{i},1)
        Pmax{i,2}(:,j)=max(P{i,2}(:,j)); 
    end
    M=max(Pmax{i,2}(:,:));
        [r,c]=find(P{i,2}==M);
        if M>0
        Pmax2{i,2}(1,1)=max(Pmax{i,2}(1,:));
            if isscalar(r)==true
            Pmax2{i,1}(1,1)=P{i,1}(c,1);
            Pmax2{i,1}(2,1)= P{i,1}(r,1);
            Intake(i,2)=Pmax2{i,1}(1,1); %Intake ID
            PowerHouse(i,2)=Pmax2{i,1}(2,1); %Power House ID
            Intake(i,3)=Pmax2{i,2}; %Available power
            PowerHouse(i,3)=Pmax2{i,2}; %Available power
            Intake(i,4)=River{i,1}(c,23); %Intake X coord.
            Intake(i,5)=River{i,1}(c,24); %Intake Y coord.
            PowerHouse(i,4)=River{i,1}(r,23); %Power House X coord.
            PowerHouse(i,5)=River{i,1}(r,24); %Power House Y coord.
            else
            r=max(r);
            c=max(c);
            end
        end
    Q_perc{i,1}=Flow_perc(23-n_flow,1);
end

%% Remove combinations that would generate less than P_min or more than P_max or are empty

%Remove combinations that are empty

emptyCells=cellfun(@isempty,Pmax2);

for i=1:length(River)
    if emptyCells(i,1)==1
        River{i,1}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);

[r,~]=find(emptyCells==1);

for i=1:length(Intake)
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

emptyCells=[];

%Remove combinations that would generate less than P_min

for i=1:length(River)
    emptyCells(i,:)=Pmax2{i,2}(:,:)<=P_min;
end

for i=1:length(River)
    if emptyCells(i,1)==1
        River{i,1}=[];
        Pmax2{i,2}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);

[r,~]=find(emptyCells==1);

for i=1:length(Intake)
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

emptyCells=[];

%Remove combinations that would generate more than P_max

for i=1:length(River)
    emptyCells(i,:)=Pmax2{i,2}(:,:)>=P_max;
end

for i=1:length(River)
    if emptyCells(i,1)==1
        River{i,1}=[];
        Pmax2{i,2}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);

[r,~]=find(emptyCells==1);

for i=1:length(Intake)
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

%Create matrices that contain the head and the flow and add them to the
%final Intake and Powerhouse matrices

Head=zeros(length(Pmax2),1);
Flow=zeros(length(Pmax2),1);

[Limitr,~]=size(Pmax2);

for i=1:Limitr
    ID_Intake=Pmax2{i,1}(1,1);
    ID_PowerHouse=Pmax2{i,1}(2,1); 
    [rI,cI]=find(River{i,1}(:,1)==ID_Intake);
    [rPH,cPH]=find(River{i,1}(:,1)==ID_PowerHouse);
    Head(i,1)=abs(River{i,1}(rI,25)-River{i,1}(rPH,25));
    Flow(i,1)=River{i,1}(rI,n_flow);
end

for i=1:Limitr
   Intake(i,6)=Flow(i,1);
   PowerHouse(i,6)=Flow(i,1);
   Intake(i,7)=Head(i,1);
   PowerHouse(i,7)=Head(i,1);
end

counts=Limitr;
Intake(:,1)=1:counts;
PowerHouse(:,1)=1:counts;

%In case the matrices are empty put NaN in them so the write to excel
%doesn't give an error

bool_Intake=isempty(Intake);
bool_PowerHouse=isempty(PowerHouse);

if bool_Intake==1
    Intake='NaN';
    PowerHouse='NaN';
end

%% Save the data to a xls file

%Check if file exists, if yes delete it

name1=strcat(folder_location,'\Intake');
name2=strcat(folder_location,'\PowerHouse');
extension='.xls';

IntakeName=[name1 Area extension];
PowerHouseName=[name2 Area extension];

if exist(IntakeName, 'file')==2
  delete(IntakeName);
end

if exist(PowerHouseName, 'file')==2
  delete(PowerHouseName);
end

col_Header_Intake={'ID','Intake_ID','Power[kW]','Intake_X','Intake_Y','Flow[m^3/s]','Head[m]'};
col_Header_PowerHouse={'ID','PowerHouse_ID','Power[kW]','PowerHouse_X','PowerHouse_Y','Flow[m^3/s]','Head[m]'};

xlswrite(IntakeName,Intake,'Sheet1','A2');
xlswrite(IntakeName,col_Header_Intake,'Sheet1','A1');

xlswrite(PowerHouseName,PowerHouse,'Sheet1','A2');
xlswrite(PowerHouseName,col_Header_PowerHouse,'Sheet1','A1');

end