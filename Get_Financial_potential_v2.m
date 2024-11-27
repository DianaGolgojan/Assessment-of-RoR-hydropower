function get_scheme_locations = get_locations(Area,P_min,P_max,folder_location);

%% Chose if you want to plot figures

print_fig=0; %default value (1=print, 0=not print)
print_NPV_fig=1;

%% Read file containing the points along the streams

location = strcat(folder_location,'\Inputs_perc.csv');

Data=readtable(location);
DataArray=Data{:,:};

DataArray=sortrows(DataArray,[2,1],{'ascend' 'ascend'});

%%Read the file containing penstock information

penstock_file='H:\01.PhD\003.Chapter2\ROR_Search\Penstocks\Penstock.xlsx';
[~,sheet_name]=xlsfinfo(penstock_file);
for k=1:numel(sheet_name)
   penstock_info{k}=xlsread(penstock_file,sheet_name{k});   
end

penstock_price_file='H:\01.PhD\003.Chapter2\ROR_Search\Penstocks\Prices_Matlab.xlsx';
[~,sheet_name2]=xlsfinfo(penstock_price_file);
for k=1:numel(sheet_name2)
   penstock_price{k}=xlsread(penstock_price_file,sheet_name2{k});   
end

Weight_SN5000(:,1)=penstock_info{1,5}(:,1); %diameter in mm
Weight_SN5000(:,2)=penstock_info{1,5}(:,8); %penstock weight in kg

Weight_SN10000(:,1)=penstock_info{1,6}(:,1); %diameter in mm
Weight_SN10000(:,2)=penstock_info{1,6}(:,8); %penstock weight in kg

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

Intake=zeros(length(Pmax2),11);
PowerHouse=zeros(length(Pmax2),11);

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
            Q_perc{i,1}=Flow_perc(23-n_flow,1);
            else
            r=max(r);
            c=max(c);
            end
        end
end

%% Remove combinations that would generate less than P_min or more than P_max or are empty

%Remove combinations that are empty

emptyCells=[];
emptyCells=cellfun(@isempty,Pmax2);

for i=1:length(River)
    if emptyCells(i,1)==1
        River{i,1}=[];
        Q_perc{i,1}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

r=[];
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
        Q_perc{i,1}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

r=[];
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
        Q_perc{i,1}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

r=[];
[r,~]=find(emptyCells==1);

for i=1:length(Intake)
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

%Create matrices that contain the head and the flow

Head=[];
Flow=[];
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

%% Choose the penstock diameter

L=zeros(length(PowerHouse),1); %penstock length in m
diam=zeros(length(PowerHouse),1); %penstock diameter in m
diam_min=zeros(length(PowerHouse),1); %min diameter for choosing nominal diameter in m
diam_max=zeros(length(PowerHouse),1); %max diameter for choosing nominal diameter in m
f=zeros(length(PowerHouse),1); %Darcy friction factor
ftemp=zeros(1000,1)+0.01; %help to calculate Darcy friction factor
diff=zeros(1000,1); %help to calculate Darcy friction factor
LHS=zeros(length(PowerHouse),1); %help to calculate Darcy friction factor
RHS=zeros(length(PowerHouse),1); %help to calculate Darcy friction factor
V=zeros(length(PowerHouse),1); %flow velocity in m/s
Re=zeros(length(PowerHouse),1); %Reynolds number
hfn=zeros(length(PowerHouse),1); %Head loss due to friction along the penstock in m
NetHead=zeros(length(PowerHouse),1); %Net head considering head losses in m
alfa=zeros(length(PowerHouse),1); %penstock wave celerity rating in m/s
hs=zeros(length(PowerHouse),1); %surge head in m
TotalHead=zeros(length(PowerHouse),1); %total head considering surge head in m
Ps=zeros(length(PowerHouse),1); %pressure in pipe in bars
PN=[6 10 16 20 25 32]; %pressure class
SN=[5000 10000]; %rigidity class
n_pipe=0.009; %Manning number for GRP pipe
niu=1.31*10^-6; %kinematic viscosity of water in m^2/s
epsilon=0.029; %roughness coefficient for GRP pipe

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
   ID_Intake=Pmax2{i,1}(1,1);
   ID_PowerHouse=Pmax2{i,1}(2,1); 
   [rI,cI]=find(River{i,1}(:,1)==ID_Intake);
   [rPH,cPH]=find(River{i,1}(:,1)==ID_PowerHouse);
   L(i,1)= sqrt((River{i,1}(rI,23)-River{i,1}(rPH,23)).^2+(River{i,1}(rI,24)-River{i,1}(rPH,24)).^2+(River{i,1}(rI,25)-River{i,1}(rPH,25)).^2); %distance between the Intake and the PowerHouse
   slope(i,1)=abs((River{i,1}(rI,25)-River{i,1}(rPH,25))/(sqrt((River{i,1}(rI,23)-River{i,1}(rPH,23)).^2+(River{i,1}(rI,24)-River{i,1}(rPH,24)).^2)));
   %diam(i,1)=2.83*((L(i,1).^2*Flow(i,1).^2)/Head(i,1)).^0.1876; %pipe diameter
   %diam2(i,1)=nthroot((257.5*(L(i,1)*epsilon^2*Flow(i,1).^2/Head(i,1))).^3,16);
   diam(i,1)=2.69*((epsilon^2*Flow(i,1).^2*L(i,1))/Head(i,1)).^0.1875; %pipe diameter
   
% The maximum diameter for GRP pipe is 4 m so Q needs to be modified to
% reach 4 m diameter
    n_flow=23;
   while diam(i,1)>4&n_flow>3
      n_flow=n_flow-1;
      ID_Intake=Pmax2{i,1}(1,1);
      [rI,cI]=find(River{i,1}(:,1)==ID_Intake);
      Flow(i,1)= River{i,1}(rI,n_flow);
      Q_perc{i,1}=Flow_perc(23-n_flow,1);
      diam(i,1)=2.69*((epsilon^2*Flow(i,1).^2*L(i,1))/Head(i,1)).^0.1875; %pipe diameter
   end
 
   diam_min(i,1)=diam(i,1)-0.1*diam(i,1);
   diam_max(i,1)=diam(i,1)+0.1*diam(i,1);
   
   for k=1:length(penstock_info{1,6}(:,1))
        if (penstock_info{1,6}(k,1)/1000>=diam_min(i,1) & penstock_info{1,6}(k,1)/1000<=diam_max(i,1))
            diam(i,1)=penstock_info{1,6}(k,1)/1000;
        end
   end
   
   if diam(i,1)>=4
      diam(i,1)=5;
   end
   
   %The minimum diameter is 100 mm so the diameter is increased to 100 mm
   %if it's smaller
   
    if diam(i,1)<=0.1
      diam(i,1)=0.1;
   end
   
   V(i,1)=4.*Flow(i,1)/(3.14*diam(i,1).^2); %average water velocity
   Re(i,1)=diam(i,1)*V(i,1)/niu; %Reynolds number
   j=1;
    while j<=1000 %calculate the Darcy friction factor
        LHS(j,1)=1/sqrt(ftemp(j,1));
        RHS(j,1)=-2*log10((epsilon/diam(i,1))./3.7+2.51./Re(i,1)*sqrt(ftemp(j,1)));
        diff(j,1)=abs(LHS(j,1)-RHS(j,1));
        j=j+1;
        ftemp(j,1)=ftemp(j-1,1)+0.0001;
        [rf,cf]=find(diff==min(diff));
    end
    f(i,1)=ftemp(rf,1);
    hfn(i,1)=f(i,1)*((L(i,1)./diam(i,1))*(V(i,1).^2/(2*g))); %headloss across flow range
    NetHead(i,1)=Head(i,1)-hfn(i,1);
    %Chose the wave celerity factor alfa
    Ps_ini(i,1)=Head(i,1)/10;
    for t=1:5
        if diam(i,1)==penstock_info{1,3}(t,2)/1000
            alfa(i,1)=penstock_info{1,3}(t,3)/1000;
        end
    end   
    for y=1:length(penstock_info{1,2}(:,1))
       if (diam(i,1)>= penstock_info{1,2}(y,1)/1000&diam(i,1)<=penstock_info{1,2}(y,2)/1000)
           if Ps_ini(i,1)<=6
          alfa(i,1)= penstock_info{1,2}(y,3);
           else
               if (Ps_ini(i,1)<=10&Ps_ini(i,1)>6)
               alfa(i,1)= penstock_info{1,2}(y,4);   
               else
                   if (Ps_ini(i,1)<=16&Ps_ini(i,1)>10)
                       alfa(i,1)= penstock_info{1,2}(y,5);
                   else
                       if (Ps_ini(i,1)<=25&Ps_ini(i,1)>16)
                       alfa(i,1)= penstock_info{1,2}(y,5);
                       else
                           alfa(i,1)= penstock_info{1,2}(y,3);
                       end
                   end
               end
           end
       end
    end             
    hs(i,1)=alfa(i,1)*V(i,1)/g; %surge head in meters
    TotalHead(i,1)=Head(i,1)+hs(i,1); %total head in meters
    Ps(i,1)=TotalHead(i,1)*0.098; %pressure in bars
    bool=1;
    for z=1:length(PN)
        if bool==1
            if Ps(i,1)<=PN(1,z)
               Ps(i,1)=PN(1,z); 
               bool=0;
            end
        end
    end
end

%Remove locations with a diameter higher than 4 or with the pressure in
%pipe higher than 32

rDiam=[];
rPs=[];
[rDiam,~]=find(diam(:,1)==5);
[rPs,~]=find(Ps(:,1)>32);

for i=1:length(Intake)
   Ps(rDiam,:)=0;
   Ps(rPs,:)=0;
   Intake(rDiam,:)=0;
   Intake(rPs,:)=0;
   PowerHouse(rDiam,:)=0;
   PowerHouse(rPs,:)=0;
   diam(rPs,:)=5;
   Flow(rDiam)=0;
   Flow(rPs)=0;
   NetHead(rDiam)=0;
   NetHead(rPs)=0;
   Head(rDiam)=0;
   Head(rPs)=0;
   slope(rDiam)=9999;
   slope(rPs)=9999;
end

for i=1:length(rDiam)
      River{rDiam(i,1),1}=[];
      Pmax2{rDiam(i,1),2}=[];
      Q_perc{rDiam(i,1),1}=[];
end

for i=1:length(rPs)
    River{rPs(i,1),1}=[];
    Pmax2{rPs(i,1),2}=[];
    Q_perc{rPs(i,1),1}=[];
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

diam(all(diam==5,2),:)=[];
Ps(all(Ps==0,2),:)=[];

Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];
Head(all(Head==0,2),:)=[];
slope(all(slope==9999,2),:)=[];

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

%% Choose the turbine

n_flow=14;

for i=1:length(NetHead)
     if (Flow(i,1)<=2&NetHead(i,1)>100)|(Flow(i,1)<0.5&(NetHead(i,1)>50&NetHead(i,1)<=100))
         Intake(i,6)=1; %1=Pelton
         PowerHouse(i,6)=1; 
     else if ((2<Flow(i,1)&Flow(i,1)<=20)&(100<NetHead(i,1)&NetHead(i,1)<=300))|((0.5<Flow(i,1)&Flow(i,1)<=20)&(20<NetHead(i,1)&NetHead(i,1)<=100))
             Intake(i,6)=2; %2=Francis
             PowerHouse(i,6)=2;            
     else if ((Flow(i,1)>20&Flow(i,1)<=50)&(NetHead(i,1)<=50&NetHead(i,1)>3))|((Flow(i,1)>0.5&Flow(i,1)<=20&(NetHead(i,1)>3&NetHead<=20))|(Flow<=0.5&(NetHead(i,1)>3&NetHead(i,1)<=50)))
              Intake(i,6)=3; %3=Kaplan
              PowerHouse(i,6)=3;
     else if (1.5<=NetHead(i,1)&NetHead(i,1)<=5)&(1<=Flow(i,1)&Flow<=20)&L(i,1)<=25
             Intake(i,6)=4; %4=Archimedean Screw
             PowerHouse(i,6)=4;
     else if (2<=NetHead(i,1)&NetHead(i,1)<=40)&(0.1<=Flow(i,1)&Flow(i,1)<5)
              Intake(i,6)=5; %5=Crossflow
              PowerHouse(i,6)=5;
     else
         Intake(i,6)=-9999; %-9999=Unknown
         PowerHouse(i,6)=-9999;
         n_flow=2;
     while Intake(i,6)==-9999&n_flow<21
           n_flow=n_flow+1;
           [rI,cI]=find(River{i,1}(:,1)==Intake(i,2));
           Flow(i,1)=River{i,1}(rI,n_flow);
           Q_perc{i,1}=Flow_perc(23-n_flow,1);
      if (Flow(i,1)<=2&NetHead(i,1)>100)|(Flow(i,1)<0.5&(NetHead(i,1)>50&NetHead(i,1)<=100))
         Intake(i,6)=1; %1=Pelton
         PowerHouse(i,6)=1; 
      else if ((2<Flow(i,1)&Flow(i,1)<=20)&(100<NetHead(i,1)&NetHead(i,1)<=300))|((0.5<Flow(i,1)&Flow(i,1)<=20)&(20<NetHead(i,1)&NetHead(i,1)<=100))
            Intake(i,6)=2; %2=Francis
            PowerHouse(i,6)=2;            
      else if ((Flow(i,1)>20&Flow(i,1)<=50)&(NetHead(i,1)<=50&NetHead(i,1)>3))|((Flow(i,1)>0.5&Flow(i,1)<=20&(NetHead(i,1)>3&NetHead<=20))|(Flow<=0.5&(NetHead(i,1)>3&NetHead(i,1)<=50)))
                Intake(i,6)=3; %3=Kaplan
                PowerHouse(i,6)=3;
      else if(1.5<=NetHead(i,1)&NetHead(i,1)<=5)&(1<=Flow(i,1)&Flow<=20)&L(i,1)<=25
                     Intake(i,6)=4; %4=Archimedean Screw
                     PowerHouse(i,6)=4;
      else if (2<=NetHead(i,1)&NetHead(i,1)<=40)&(0.1<=Flow(i,1)&Flow(i,1)<5)
                        Intake(i,6)=5; %5=Crossflow
                        PowerHouse(i,6)=5;
      end
      end
      end
      end
      end
      end
      end
      end
      end
      end
      end
end

%Reaction turbines peak efficiency calculation (Francis=2 and Kaplan=3)

it_no=12;
ep=zeros(length(PowerHouse),1); %peak efficiency
Qp=zeros(length(PowerHouse),1); %peak efficiency flow
nq=zeros(length(PowerHouse),1); %specific speed
enq_=zeros(length(PowerHouse),1); %adjustment to turbine peak efficiency
ed_=zeros(length(PowerHouse),1); %adjustment to turbine peak efficiency
d=zeros(length(PowerHouse),1); %runner size
Q=zeros(it_no,1); %flows for calculating efficiency curves
eq=zeros(it_no,1); %efficiency at flows below and/or above peak flow
eg=0.98; %generator efficiency
l_trans=0.02; %percentage losses expected from the transformer
l_para=0.02; %percentage losses expected from the transmiossion lines
l_dt=0.02; %downtime experienced during 1 year of operation

[Limitr,~]=size(PowerHouse);

% Francis turbine
 
for i=1:Limitr
    if PowerHouse(i,6)==2
        if Flow(i,1)<1.8
            k=0.46; %empirical constant
        else
        k=0.41;
        end
    d(i,1)=k*Flow(i,1).^0.473;
    kq=600; %constant determined by turbine type
    nq(i,1)=kq*(NetHead(i,1).^(-0.5));
    enq_(i,1)=((nq(i,1)-56)/256).^2;
    ed_(i,1)=(0.081+enq_(i,1))*(1-0.789*d(i,1).^(-0.2));
    Rm=4.5; %turbine manufacturer/design coefficient, default value
    ep(i,1)=(0.919-enq_(i,1)+ed_(i,1))-0.0305+0.005*Rm;
    if ep(i,1)>=1
           ep(i,1)=0.9; 
    end
    Qp(i,1)=0.65*Flow(i,1)*nq(i,1).^0.05; %flow at peak efficiency
    ep_drop(i,1)=0.0072*nq(i,1)^0.4; %drop in effieicny at full load
    er(i,1)=(1-ep_drop(i,1))*ep(i,1); %efficiency at full load
    for j=1:it_no
        Q(j,1)=(0+j/10)*Flow(i,1);
        x_axis(j,1)=Q(j,1)/Flow(i,1);
        if Q(j,1)<Qp(i,1)
            eq(j,1)=(1-(1.25*((Qp(i,1)-Q(j,1))/Qp(i,1)).^(3.94-0.0195*nq(i,1))))*ep(i,1); %efficiency at flow below peak efficiency flow
        else
            eq(j,1)=ep(i,1)-((((Q(j,1)-Qp(i,1))/(Flow(i,1)-Qp(i,1))).^2)*(ep(i,1)-er(i,1))); %efficiency at flow above peak efficiency flow
        end
    end
    
    P_des(i,1)=(ro*g*Flow(i,1)*NetHead(i,1)*ep(i,1)*eg*(1-l_trans)*(1-l_para))/1000; %power production in kW
    [r_ID,~]=find(River{i}(:,1)==Intake(i,2));
    Pn=zeros(20,1);
    E_avail=zeros(20,1);
    for k=1:20
        if River{i,1}(r_ID,k+2)<Qp(i,1)
            eq(k,1)=(1-(1.25*((Qp(i,1)-River{i,1}(r_ID,k+2))/Qp(i,1)).^(3.94-0.0195*nq(i,1))))*ep(i,1); %efficiency at flow below peak efficiency flow
        else
            eq(k,1)=ep(i,1)-((((River{i,1}(r_ID,k+2)-Qp(i,1))/(Flow(i,1)-Qp(i,1))).^2)*(ep(i,1)-er(i,1))); %efficiency at flow above peak efficiency flow
        end
       Pn(k,1)=(ro*g*River{i,1}(r_ID,k+2)*NetHead(i,1)*eq(k,1)*eg*(1-l_trans)*(1-l_para))/1000; %power production curve in kW
       if k==1
        E_avail(1,1)=Pn(1,1)/2*5/100*8760*(1-l_dt);
       else
           if k>1
                E_avail(k,1)=(Pn(k-1,1)+Pn(k,1))/2*5/100*8760*(1-l_dt);
           end
       end
    end
    %% Figure
    if print_fig==1
        clf %clears figure information
        numb=num2str(PowerHouse(i,2));
        save_name=strcat(folder_location,'\Francis',numb,'.png');
        grid on;
        grid minor;
        flow=num2str(round(Flow(i,1),2));
        nethead=num2str(round(NetHead(i,1)));
        title_name=strcat("Francis Efficiency Curve with Q=",flow," cum/s and Head=",nethead," m");  
        title(title_name);
        xlabel("Q/Qd");
        ylabel("Efficiency");
        hold on
        xlim([0 x_axis(it_no,1)])
        ylim([0 1])
        hold on
        %plot(x_axis(:,1),eq(:,1),'Color','[0.3 0 1]','LineWidth',2,'LineStyle','-','Marker','none');
        line(x_axis(:,1),eq(:,1));
        hold on
        if exist(save_name, 'file')==2
            delete(save_name);
        end
            saveas(gcf,save_name);
    end 
    end
end

% Kaplan turbine

Q=zeros(it_no,1); %flows for calculating efficiency curves
eq=zeros(it_no,1); %efficiency at flows below and/or above peak flow

for i=1:Limitr
    if PowerHouse(i,6)==3
        if Flow(i,1)<1.8
            k=0.46; %empirical constant
        else
        k=0.41;
        end
    d(i,1)=k*Flow(i,1).^0.473;
    kq=800; %constant determined by turbine type
    nq(i,1)=kq*(NetHead(i,1).^(-0.5));
    enq_(i,1)=((nq(i,1)-170)/700).^2;
    ed_(i,1)=(0.095+enq_(i,1))*(1-0.789*d(i,1).^(-0.2));
    Rm=4.5; %turbine manufacturer/design coefficient, default value
    ep(i,1)=(0.905-enq_(i,1)+ed_(i,1))-0.0305+0.005*Rm;
    if ep(i,1)>=1
           ep(i,1)=0.9; 
    end
    Qp(i,1)=0.75*Flow(i,1);
    for j=1:it_no
        Q(j,1)=(0+j/10)*Flow(i,1);
        x_axis(j,1)=Q(j,1)/Flow(i,1);
        eq(j,1)=(1-3.5*((Qp(i,1)-Q(j,1))/Qp(i,1)).^6)*ep(i,1); %efficiency at flows below and above peak efficiency flows   
    end
    %% Figure
    if print_fig==1
    clf %clears figure information
    numb=num2str(PowerHouse(i,2));
    save_name=strcat(folder_location,'\Kaplan',numb,'.png');
    grid on;
    grid minor;
    flow=num2str(round(Flow(i,1),2));
    nethead=num2str(round(NetHead(i,1)));
    title_name=strcat("Kaplan Efficiency Curve with Q=",flow," cum/s and Head=",nethead," m");  
    title(title_name);
    xlabel("Q/Qd");
    ylabel("Efficiency");
    hold on
    xlim([0 x_axis(it_no,1)])
    ylim([0 1])
    hold on
    %plot(x_axis(:,1),eq(:,1),'Color','[0.3 0 1]','LineWidth',2,'LineStyle','-','Marker','none');
    line(x_axis(:,1),eq(:,1));
    hold on
    if exist(save_name, 'file')==2
        delete(save_name);
    end
    saveas(gcf,save_name);
    end
    end
end

% Impulse turbines peak efficiency calculation (Pelton=1)

jets=1; %number of jets for Pelton
n_rot=zeros(length(PowerHouse),1); %rotational speed
Q=zeros(it_no,1); %flows for calculating efficiency curves
eq=zeros(it_no,1); %efficiency at flows below and/or above peak flow

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    if PowerHouse(i,6)==1
        n_rot(i,1)=31*((NetHead(i,1)*Flow(i,1))/jets).^0.5;
        d(i,1)=(49.4*NetHead(i,1).^0.5*jets^0.02)/n_rot(i,1);
        ep(i,1)=0.864*d(i,1).^0.04;
        Qp(i,1)=(0.662+0.001*jets)*Flow(i,1);
        while ep(i,1)<=0.7&jets<6&ep(i,1)<1
            jets=jets+1;
            n_rot(i,1)=((NetHead(i,1)*Flow(i,1))/jets).^0.5;
            d(i,1)=(49.4*NetHead(i,1).^0.5*jets^0.02)/n_rot(i,1);
            ep_temp=0.864*d(i,1).^0.04;
            if ep_temp>ep(i,1)&ep_temp<1
                ep(i,1)=ep_temp;
            else
                ep(i,1)=0.864*d(i,1).^0.04;
            end
            Qp(i,1)=(0.662+0.001*jets)*Flow(i,1);
        end
        for j=1:it_no
            Q(j,1)=(0+j/10)*Flow(i,1);
            x_axis(j,1)=Q(j,1)/Flow(i,1);
            eq(j,1)=(1-((1.31+0.025*jets)*abs((Qp(i,1)-Q(j,1))/Qp(i,1)).^(5.6+0.4*jets)))*ep(i,1); %efficiency at flows above and below peak efficiency flows
        end
        %% Figure
        if print_fig==1
        clf %clears figure information
        numb=num2str(PowerHouse(i,2));
        save_name=strcat(folder_location,'\Pelton',numb,'.png');
        grid on;
        grid minor;
        flow=num2str(round(Flow(i,1),2));
        nethead=num2str(round(NetHead(i,1)));
        title_name=strcat("Pelton Efficiency Curve with Q=",flow," cum/s and Head=",nethead," m");  
        title(title_name);
        xlabel("Q/Qd");
        ylabel("Efficiency");
        hold on
        xlim([0 x_axis(it_no,1)])
        ylim([0 1])
        hold on
        %plot(x_axis(:,1),eq(:,1),'Color','[0.3 0 1]','LineWidth',2,'LineStyle','-','Marker','none');
        line(x_axis(:,1),eq(:,1));
        hold on
        if exist(save_name, 'file')==2
            delete(save_name);
        end
        saveas(gcf,save_name);
        end   
    end
end
jets=1; %number of jets for Pelton needs to be reset for future calculations

% Impulse turbines peak efficiency calculation (Cross-flow=5)

[Limitr,~]=size(PowerHouse);
Q=zeros(it_no,1); %flows for calculating efficiency curves
eq=zeros(it_no,1); %efficiency at flows below and/or above peak flow

for i=1:Limitr
    if PowerHouse(i,6)==5
        ep(i,1)=0.816;
        Qp(i,1)=0.7*Flow(i,1);
        for j=1:it_no
            Q(j,1)=(0+j/10)*Flow(i,1);
            x_axis(j,1)=Q(j,1)/Flow(i,1);
            x=Q(j,1)/Flow(i,1);
            if Q(j,1)<Qp(i,1)*0.103
                eq(j,1)=-925.33841*x^2+188.14434*x-8.95889; %efficiency at flow below peak efficiency flow
            else
                if Q(j,1)>Qp(i,1)*0.103 && Q(j,1)<Qp(i,1)*0.145
                    eq(j,1)=2.38095*x+0.35476; %efficiency at flow below peak efficiency flow
                else
                    if Q(j,1)>Qp(i,1)*0.145
                        eq(j,1)=-0.38935*x^2+0.54132*x+0.63050; %efficiency at flow above peak efficiency flow 
                    end
                end
            end
            if j<=10
            eq_rets(j,1)=0.79-0.15*((Flow(i,1)-Q(j,1))/Flow(i,1))-1.37*((Flow(i,1)-Q(j,1))/Flow(i,1)).^14; %efficiency at flow above and below peak flow fro RETS
            end
        end
     %% Figure
     if print_fig==1
        clf %clears figure information
        numb=num2str(PowerHouse(i,2));
        save_name=strcat(folder_location,'\Cross-flow',numb,'.png');
        grid on;
        grid minor;
        flow=num2str(round(Flow(i,1),2));
        nethead=num2str(round(NetHead(i,1)));
        title_name=strcat("Cross-flow Efficiency Curve with Q=",flow," m^{3}/s and Head=",nethead," m");  
        title(title_name);
        xlabel("Q/Qd");
        ylabel("Efficiency");
        hold on
        xlim([0 x_axis(it_no,1)])
        ylim([0 1])
        hold on
        %plot(x_axis(:,1),eq(:,1),'Color','[0.3 0 1]','LineWidth',2,'LineStyle','-','Marker','none');
        line(x_axis(:,1),eq(:,1));
        line(x_axis((1:10),1),eq_rets(:,1),'Color','red');
        hold off
        legend('Experimental results','Formula based results');
        if exist(save_name, 'file')==2
            delete(save_name);
        end
        saveas(gcf,save_name);
     end
    end
end

% Archimedean Screw=4 turbine peak efficiency calculation

[Limitr,~]=size(PowerHouse);
Q=zeros(it_no,1); %flows for calculating efficiency curves
eq=zeros(it_no,1); %efficiency at flows below and/or above peak flow

for i=1:Limitr
    if PowerHouse(i,6)==4
        % calculate equivalent runner diameter to be able to calculate the
        % cost
        if Flow(i,1)<1.8
            k=0.46; %empirical constant
        else
        k=0.41;
        end
    d(i,1)=k*Flow(i,1).^0.473;
        ep(i,1)=0.828;
        Qp(i,1)=0.78*Flow(i,1);
        for j=4:it_no
            Q(j,1)=(0+j/10)*Flow(i,1);
            x_axis(j,1)=Q(j,1)/Flow(i,1);
            x=Q(j,1)/Flow(i,1);
            if Q(j,1)<Qp(i,1)
                eq(j,1)=-1.31763*x^2+2.0858*x+0.0346; %efficiency at flow below peak efficiency flow
            else
                eq(j,1)=-565.6932290792*x^6 + 3345.5341490459*x^5 - 8204.1750190284*x^4 + 10676.6071388550*x^3 - 7775.4512628411*x^2+3004.2384250625*x-480.2389292428; %efficiency at flow above peak efficiency flow 
            end
        end
         %% Figure
         if print_fig==1
        clf %clears figure information
        numb=num2str(PowerHouse(i,2));
        save_name=strcat(folder_location,'\Archimedean screw',numb,'.png');
        grid on;
        grid minor;
        flow=num2str(round(Flow(i,1),2));
        nethead=num2str(round(NetHead(i,1)));
        title_name=strcat("Archimedean screw Efficiency Curve with Q=",flow," cum/s and Head=",nethead," m");  
        title(title_name);
        xlabel("Q/Qd");
        ylabel("Efficiency");
        hold on
        xlim([0 x_axis(it_no,1)])
        ylim([0 1])
        hold on
        %plot(x_axis(:,1),eq(:,1),'Color','[0.3 0 1]','LineWidth',2,'LineStyle','-','Marker','none');
        line(x_axis(:,1),eq(:,1));
        hold on
        if exist(save_name, 'file')==2
            delete(save_name);
        end
        saveas(gcf,save_name);
         end
    end
end

%% Remove schemes with less than 0.7 peak efficiency or don't have a
%% suitable turbine

%emptyCells=zeros(length(PowerHouse),1);
emptyCells=[];
Limitr=[];

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    emptyCells(i,:)=ep(i,1)>0.7;
end

for i=1:length(emptyCells)
    if emptyCells(i,1)==0
        Q_perc{i,1}=[];
    end
end

r=[];
[r,~]=find(emptyCells==0);

for i=1:length(emptyCells)
    if emptyCells(i,1)==0
        River{i,1}=[];
        Pmax2{i,2}=[];
        Q_perc{i,1}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

r=[];
[r,~]=find(emptyCells==0);

Intake(r,:)=0;
PowerHouse(r,:)=0;
Flow(r,:)=0;
NetHead(r,:)=0;
Head(r,:)=0;
ep(r,:)=0;
diam(r,:)=0;
Ps(r,:)=0;
d(r,:)=9999;
slope(r,:)=9999;

Intake(all(~Intake,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];
Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];
Head(all(Head==0,2),:)=[];
ep(all(ep==0,2),:)=[];
Qp(all(Qp==0,2),:)=[];
diam(all(diam==0,2),:)=[];
Ps(all(Ps==0,2),:)=[];
d(all(d==9999,2),:)=[];
slope(all(slope==9999,2),:)=[];

Intake(:,7)=Flow(:,1);
Intake(:,8)=NetHead(:,1);
Intake(:,9)=ep(:,1);
Intake(:,10)=diam(:,1);

PowerHouse(:,7)=Flow(:,1);
PowerHouse(:,8)=NetHead(:,1);
PowerHouse(:,9)=ep(:,1);
PowerHouse(:,10)=diam(:,1);

%% Recalculate power considering the peak efficiency and flow

Limitr=[];
[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    if ep(i,1)~=0
       Intake(i,3)=ro*g*Qp(i,1)*NetHead(i,1)*ep(i,1)/1000; %power in kW 
       PowerHouse(i,3)=round(abs(ro*g*Qp(i,1)*NetHead(i,1)*ep(i,1))/1000); %power in kW
    end
end

%Eliminate the combinations that would generate less than P_min or more
%than P_max

%Remove combinations that would generate less than P_min

emptyCells=[];
[Limitr,~]=size(Intake);

for i=1:Limitr
    emptyCells(i,:)=Intake(i,3)<=P_min;
end

for i=1:Limitr
    if emptyCells(i,1)==1
        River{i,1}=[];
        Pmax2{i,2}=[];
        Q_perc{i,1}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

r=[];
[r,~]=find(emptyCells==1);

for i=1:Limitr
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
   Qp(r,:)=0;
   Flow(r,:)=0;
   NetHead(r,:)=0;
   Head(r,:)=0;
   ep(r,:)=0;
   diam(r,:)=0;
   Ps(r,:)=0;
   d(r,:)=9999;
   slope(r,:)=9999;   
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];
Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];
Head(all(Head==0,2),:)=[];
ep(all(ep==0,2),:)=[];
Qp(all(Qp==0,2),:)=[];
diam(all(diam==0,2),:)=[];
Ps(all(Ps==0,2),:)=[];
d(all(d==9999,2),:)=[];
slope(all(slope==9999,2),:)=[];

emptyCells=[];

%Remove combinations that would generate more than P_max

Limitr=[];
[Limitr,~]=size(Intake);

for i=1:Limitr
    emptyCells(i,:)=Intake(i,3)>=P_max;
end

for i=1:Limitr
    if emptyCells(i,1)==1
        River{i,1}=[];
        Pmax2{i,2}=[];
        Q_perc{i,1}=[];
    end
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

r=[];
[r,~]=find(emptyCells==1);

for i=1:Limitr
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
   Flow(r,:)=0;
   NetHead(r,:)=0;
   Head(r,:)=0;   
   ep(r,:)=0;
   diam(r,:)=0;
   Qp(r,:)=0;
   Ps(r,:)=0;
   d(r,:)=9999;
   slope(r,:)=9999;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];
Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];
Head(all(Head==0,2),:)=[];
ep(all(ep==0,2),:)=[];
Qp(all(Qp==0,2),:)=[];
diam(all(diam==0,2),:)=[];
Ps(all(Ps==0,2),:)=[];
d(all(d==9999,2),:)=[];
slope(all(slope==9999,2),:)=[];

%% Recalculate the diameter in case the flow has changed to make the scheme
%% fit a turbine

Limitr=[];
[Limitr,~]=size(PowerHouse);

diam_min=[];
diam_max=[];
diam=[];
for i=1:Limitr
   ID_Intake=Pmax2{i,1}(1,1);
   ID_PowerHouse=Pmax2{i,1}(2,1); 
   [rI,cI]=find(River{i,1}(:,1)==ID_Intake);
   [rPH,cPH]=find(River{i,1}(:,1)==ID_PowerHouse);
   L(i,1)= sqrt((River{i,1}(rI,23)-River{i,1}(rPH,23)).^2+(River{i,1}(rI,24)-River{i,1}(rPH,24)).^2+(River{i,1}(rI,25)-River{i,1}(rPH,25)).^2); %distance between the Intake and the PowerHouse
   %diam(i,1)=2.83*((L(i,1).^2*Flow(i,1).^2)/Head(i,1)).^0.1876; %pipe diameter
   %diam2(i,1)=nthroot((257.5*(L(i,1)*epsilon^2*Flow(i,1).^2/Head(i,1))).^3,16);
   diam(i,1)=2.69*((epsilon^2*Flow(i,1).^2*L(i,1))/Head(i,1)).^0.1875; %pipe diameter
    
   diam_min(i,1)=diam(i,1)-0.20*diam(i,1);
   diam_max(i,1)=diam(i,1)+0.20*diam(i,1);
   
   for k=1:length(penstock_info{1,6}(:,1))
        if (penstock_info{1,6}(k,1)/1000>=diam_min(i,1) & penstock_info{1,6}(k,1)/1000<=diam_max(i,1))
            diam(i,1)=penstock_info{1,6}(k,1)/1000;
        end
   end
   
   if diam(i,1)>=4
      diam(i,1)=5;
   end
   
   %The minimum diameter is 100 mm so the diameter is increased to 100 mm
   %if it's smaller
   
    if diam(i,1)<=0.1
      diam(i,1)=0.1;
   end
   
   V(i,1)=4.*Flow(i,1)/(3.14*diam(i,1)); %average water velocity
   Re(i,1)=diam(i,1)*V(i,1)/niu; %Reynolds number
   j=1;
    while j<=1000 %calculate the Darcy friction factor
        LHS(j,1)=1/sqrt(ftemp(j,1));
        RHS(j,1)=-2*log10((epsilon/diam(i,1))./3.7+2.51./Re(i,1)*sqrt(ftemp(j,1)));
        diff(j,1)=abs(LHS(j,1)-RHS(j,1));
        j=j+1;
        ftemp(j,1)=ftemp(j-1,1)+0.0001;
        [rf,cf]=find(diff==min(diff));
    end
    f(i,1)=ftemp(rf,1);
    hfn(i,1)=f(i,1)*((L(i,1)./diam(i,1))*(V(i,1).^2/(2*g))); %headloss across flow range
    NetHead(i,1)=Head(i,1)-hfn(i,1);
    %Chose the wave celerity factor alfa
    Ps_ini(i,1)=Head(i,1)/10;
    for t=1:5
        if diam(i,1)==penstock_info{1,3}(t,2)/1000
            alfa(i,1)=penstock_info{1,3}(t,3)/1000;
        end
    end   
    for y=1:length(penstock_info{1,2}(:,1))
       if (diam(i,1)>= penstock_info{1,2}(y,1)/1000&diam(i,1)<=penstock_info{1,2}(y,2)/1000)
           if Ps_ini(i,1)<=6
          alfa(i,1)= penstock_info{1,2}(y,3);
           else
               if (Ps_ini(i,1)<=10&Ps_ini(i,1)>6)
               alfa(i,1)= penstock_info{1,2}(y,4);   
               else
                   if (Ps_ini(i,1)<=16&Ps_ini(i,1)>10)
                       alfa(i,1)= penstock_info{1,2}(y,5);
                   else
                       if (Ps_ini(i,1)<=25&Ps_ini(i,1)>16)
                       alfa(i,1)= penstock_info{1,2}(y,5);
                       else
                           alfa(i,1)= penstock_info{1,2}(y,3);
                       end
                   end
               end
           end
       end
    end             
    hs(i,1)=alfa(i,1)*V(i,1)/g; %surge head in meters
    TotalHead(i,1)=Head(i,1)+hs(i,1); %total head in meters
    Ps(i,1)=TotalHead(i,1)*0.098; %pressure in bars
    bool=1;
    for z=1:length(PN)
        if bool==1
            if Ps(i,1)<=PN(1,z)
               Ps(i,1)=PN(1,z); 
               bool=0;
            end
        end
    end
end

%Remove locations with a diameter higher than 4 or with the pressure in
%pipe higher than 32

rDiam=[];
rPs=[];
rPs2=[];
[rDiam,~]=find(diam==5);
[rPs,~]=find(Ps(:,1)>32);
[rPs2,~]=find(isnan(Ps(:,1)));

for i=1:length(Intake)
   Ps(rDiam,:)=0;
   Ps(rPs,:)=0;
   Ps(rPs2,:)=0;
   Intake(rDiam,:)=0;
   Intake(rPs,:)=0;
   Intake(rPs2,:)=0;
   PowerHouse(rDiam,:)=0;
   PowerHouse(rPs,:)=0;
   PowerHouse(rPs2,:)=0;
   diam(rPs,:)=5;
   diam(rPs2,:)=5;
   Flow(rDiam)=0;
   Flow(rPs)=0;
   Flow(rPs2)=0;
   NetHead(rDiam)=0;
   NetHead(rPs)=0;
   NetHead(rPs2)=0;
   Head(rDiam)=0;
   Head(rPs)=0;
   Head(rPs2)=0;   
   ep(rDiam)=0;
   ep(rPs)=0;
   ep(rPs2)=0;
   Qp(rDiam)=0;
   Qp(rPs)=0;
   Qp(rPs2)=0;
   d(rDiam)=9999;
   d(rPs)=9999;
   d(rPs2)=9999;
   slope(rDiam)=9999;
   slope(rPs)=9999;
   slope(rPs2)=9999;   
end

for i=1:length(rDiam)
      River{rDiam(i,1),1}=[];
      Pmax2{rDiam(i,1),2}=[];
      Q_perc{rDiam(i,1),1}=[];
end

for i=1:length(rPs)
    River{rPs(i,1),1}=[];
    Pmax2{rPs(i,1),2}=[];
    Q_perc{rPs(i,1),1}=[];
end

for i=1:length(rPs2)
    River{rPs2(i,1),1}=[];
    Pmax2{rPs2(i,1),2}=[];
    Q_perc{rPs2(i,1),1}=[];
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

diam(all(diam==5,2),:)=[];
Ps(all(Ps==0,2),:)=[];
d(all(d==9999,2),:)=[];
slope(all(slope==9999,2),:)=[];

Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];
Head(all(Head==0,2),:)=[];

ep(all(ep==0,2),:)=[];
Qp(all(Qp==0,2),:)=[];

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

[Limitr,~]=size(diam);

counts=Limitr;

%Get the penstock weight

SN_pipe=zeros(counts,1);
W=zeros(counts,1);

for i=1:counts
   if diam(i,1)>=300/1000 & Ps(i,1)<=16 & Ps(i,1)>=6
       SN_pipe(i,1)=SN(1,1);
       [rp,~]=find(Weight_SN5000(:,1)==diam(i,1)*1000);
       W(i,1)=Weight_SN5000(rp,2);
   else
       if diam(i,1)>300/1000 & diam(i,1)<=1600/1000 & Ps(i,1)>=20 & Ps(i,1)<=25
            SN_pipe(i,1)=SN(1,1);
            [rp,~]=find(Weight_SN5000((1:16),1)==diam(i,1)*1000);
            W(i,1)=Weight_SN5000(rp,2);
       else
           if diam(i,1)<300/1000 && Ps(i,1)>=6 && Ps(i,1)<=20
            SN_pipe(i,1)=SN(1,2);
            [rp,~]=find(Weight_SN10000((1:4),1)==diam(i,1)*1000);
            W(i,1)=Weight_SN10000(rp,2);
           else
               if Ps(i,1)>=32 & diam(i,1)>=300/1000 & diam(i,1)<=1600/1000
                    SN_pipe(i,1)=SN(1,2);
                    [rp,~]=find(Weight_SN10000((5:20),1)==diam(i,1)*1000);
                    W(i,1)=Weight_SN10000(rp,2);
               else
                   SN_pipe(i,1)=NaN;
                   W(i,1)=NaN;
               end
           end
       end
   end
end

rPs=[];
[rPs,~]=find(isnan(W(:,1)));

for i=1:length(Intake)
   Ps(rPs,:)=0;
   Intake(rPs,:)=0;
   PowerHouse(rPs,:)=0;
   diam(rPs,:)=5;
   Flow(rPs)=0;
   NetHead(rPs)=0;
   Head(rPs)=0;
   ep(rPs)=0;
   Qp(rPs)=0;
   d(rPs)=9999;
   slope(rPs)=9999;
end

for i=1:length(rPs)
    River{rPs(i,1),1}=[];
    Pmax2{rPs(i,1),2}=[];
    Q_perc{rPs(i,1),1}=[];
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

diam(all(diam==5,2),:)=[];
Ps(all(Ps==0,2),:)=[];
d(all(d==9999,2),:)=[];
slope(all(slope==9999,2),:)=[];

Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];
Head(all(Head==0,2),:)=[];

ep(all(ep==0,2),:)=[];
Qp(all(Qp==0,2),:)=[];

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];
W(isnan(W))=[];

%% Costing

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    i_inf=2.5/100; %inflation rate
    n_inf=0; %number of years to be inflated
    C_ex=0.58; %convertion rate from canadian dollars to british pounds
    C_inf=((1+ i_inf)^n_inf)*C_ex; %conversion factor taking into account the inflation

% Turbine, generator and control
    if P_min==5 || P_min==0.1
        MWu(i,1)=7.53*Flow(i,1)*Head(i,1)/1000; %MW/unit capacity per unit
    else if P_min==100
            MWu(i,1)=7.79*Flow(i,1)*Head(i,1)/1000; %MW/unit capacity per unit
        else if P_min==1000
                MWu(i,1)=8.22*Flow(i,1)*Head(i,1)/1000; %MW/unit capacity per unit
            end
        end
    end
    n_turbine=1; %number of turbines
    Pd(i,1)=Intake(i,3)/1000; %design power in MW
    C1(i,1)=0.82*n_turbine^0.96*(Pd(i,1)/Head(i,1)^0.28)*10^6*C_inf; %generator costs

    if Intake(i,6)==3 % Kaplan
        C2(i,1)=0.27*n_turbine^0.96*d(i,1)^1.47*(1.17*Head(i,1).^0.12+2)*10^6*C_inf; %cost of a Kaplan turbine and governor
    else
        if Intake(i,6)==2 %Francis
            C2(i,1)=0.17*n_turbine^0.96*d(i,1)^1.47*((13+0.01*Head(i,1)).^0.3+3)*10^6*C_inf; %cost of a Francis turbine and governor
        else
            if Intake(i,6)==1 %Pelton
                C2(i,1)=3.47*n_turbine*0.96*(MWu(i,1)/Head(i,1).^0.5).^0.44*10^6*0.5*C_inf; %cost of a Pelton turbine and governor
            else
                if Intake(i,6)==5 %Cross-flow
                    if MWu(i,1)/Head(i,1).^0.5>0.4
                C2(i,1)=0.5*3.47*n_turbine*0.96*(MWu(i,1)/Head(i,1).^0.5).^0.44*10^6*C_inf; %cost of a Cross-flow turbine and governor
                    else
                    C2(i,1)=0.5*5.34*n_turbine*0.96*(MWu(i,1)/Head(i,1).^0.5).^0.91*10^6*C_inf;
                    end
                end
                if Intake(i,6)==4 %Archimedean Screw
                   C2(i,1)=0.9*(0.27*n_turbine^0.96*d(i,1)^1.47*(1.17*Head(i,1).^0.12+2)*10^6*C_inf); %cost of a Kaplan turbine and governor
                end
            end
        end
    end
    C3(i,1)=0.15*(C1(i,1)+C2(i,1)); %installation costs

% Penstock
    if Ps(i,1)==6
        k=1;
    else
        if Ps(i,1)==10
            k=2;
        else
            if Ps(i,1)==16
                k=3;
            else
                if Ps(i,1)==20
                    k=4;
                else
                    if Ps(i,1)==25
                        k=5;
                    else
                        if Ps(i,1)==32
                            k=6;
                        end
                    end
                end
            end
        end
    end
    [diam_price,~]=find(penstock_price{1,k}(:,:)==diam(i,1)*1000);
    Penstock_length(i)=25*abs(Intake(i,2)-PowerHouse(i,2)); %length alonside the riverbed where the penstock would be in m
    Penstock_length(i)=sqrt(Penstock_length(i).^2+(1+slope(i,1).^2));
    C4(i,1)=penstock_price{1,k}(diam_price,2)*Penstock_length(i); %penstock cost

    %Archimedean Screw turbine doesn't require a penstock
    if Intake(i,6)==4 %Archimedean Screw
        C4(i,1)=0;
    end
    
% Weight
    C5(i,1)=5*W(i,1).^0.88*C_inf; %installation costs based on penstock weight

    %Archimedean Screw turbine doesn't require a penstock
    if Intake(i,6)==4 %Archimedean Screw
        C5(i,1)=0;
    end

% Civil Structures
    lb=0.5; %distance to borrow pits in km
    ld=10; %length of crest of the dam or weir in m (assumed)
    C6(i,1)=1.97*n_turbine^-0.04*(Pd(i,1)/Head(i,1).^0.3)*(1+0.01*lb)*(1+0.005*ld/Head(i,1))*10^6*C_inf; %civil structures 
%include the dam/weir, intake works, powerhouse and tailrace

% Substation and transformer
    V=33; %transmission voltage 33 kV
    n_gen=1;
    C7(i,1)=0.0025*n_gen^0.95+0.002*(n_gen+1)*(Pd(i,1)/0.95).^0.9*V^0.3*10^6*C_inf; %substation and transformer costs

    C8(i,1)=0.15*C7(i,1); %installation costs

% Transmission lines
    Ld_line=1; %distance in m to the nearest 33 kV transmission line %%%%%%%%%(need to get this)
    Ld_sub=1; %distance in m to the nearest primary substation %%%%%%%%%(need to get this)
    
    Connection_bool(i,1)=1;
    C9(i,1)=1*73.25*1.25+150000; %considering 1 km as the distance

% Access costs
    T=0.25; %cost reduction factor set at 0.25 for unpaved roads
    A=2; %access difficulty set between 1 and 6
    road_dist=[];
    [road_dist,~]=find(River{i,1}(:,1)==Intake(i,2));
    if isempty(road_dist)
        la(i,1)=1;
    else
        la(i,1)=River{i,1}(road_dist,26)/1000; %access road length in km
    end

    C10(i,1)=0.25*T*A.^2*la(i,1).^0.9*10.^6*C_inf; % access road costs

% Other costs
    E=1; %engineering cost factor 0.67 if existing dam and 1 if no dam
    interest=10/100; %interest rate of finance

    C11(i,1)=0.37*n_turbine^0.1*E*(Pd(i,1)/Head(i,1).^0.3)*10^6; %engineering design cost
    C12(i,1)=0.04*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)); %development costs
    C13(i,1)=0.032*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)+C12(i,1)); %feasibility study cost
    C14(i,1)=0.25*interest*Flow(i,1)^0.35*1.1*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)+C12(i,1))+0.1*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)); %miscellaneous costs
    C_tot(i,1)=C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)+C12(i,1)+C13(i,1)+C14(i,1);
    Cost_per_kW(i,1)=C_tot(i,1)/Intake(i,3); %cost per kW intalled
    
% Get the percetage of costs
    C1_perc(i,1)=C1(i,1)/C_tot(i,1)*100;
    C2_perc(i,1)=C2(i,1)/C_tot(i,1)*100;
    C3_perc(i,1)=C3(i,1)/C_tot(i,1)*100;
    C4_perc(i,1)=C4(i,1)/C_tot(i,1)*100;
    C5_perc(i,1)=C5(i,1)/C_tot(i,1)*100;
    C6_perc(i,1)=C6(i,1)/C_tot(i,1)*100;
    C7_perc(i,1)=C7(i,1)/C_tot(i,1)*100;
    C8_perc(i,1)=C8(i,1)/C_tot(i,1)*100;
    C9_perc(i,1)=C9(i,1)/C_tot(i,1)*100;
    C10_perc(i,1)=C10(i,1)/C_tot(i,1)*100;
    C11_perc(i,1)=C11(i,1)/C_tot(i,1)*100;
    C12_perc(i,1)=C12(i,1)/C_tot(i,1)*100;
    C13_perc(i,1)=C13(i,1)/C_tot(i,1)*100;
    C14_perc(i,1)=C14(i,1)/C_tot(i,1)*100;
        
% Electricity generated on an average year
    Discount_factor=[];
    Discount_costs=[];
    Discount_income_year=[];
    Cumulative_NPV=[];
    NPV=[];
    Lifespan=25+3;
    Capacity_factor=40/100;
    Wholesale_price=45; %price in pounds per MWh
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Discount_rate=4/100; %discount rate - chosen
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CAPEX(i,1)=C_tot(i,1); %capex cost in year -1 (considered 1 in the for loop)
    OPEX(i,1)=3/100*C_tot(i,1); %operational and maitenance costs assumed to be 3% of the total installed costs
    Income_year(i,1)=(Intake(i,3)/1000*24*365)*Capacity_factor*Wholesale_price;%yield per year considering 
    % capacity factor of 40% and the wholesale electricity price 45 pounds/MWh
    for year=1:(Lifespan)
        Discount_factor(year,1)=1/(1+Discount_rate)^(year-3);
        if year==1 %actualy year -2
            Discount_costs(year,1)=0;
            Discount_income_year(year,1)=0;
        else
            if year==2 %actaully year -1
                Discount_costs(year,1)=CAPEX(i,1)*Discount_factor(year,1);
                Discount_income_year(year,1)=0;
            else %actually starting with year 0
                Discount_costs(year,1)=OPEX(i,1)*Discount_factor(year,1);
                Discount_income_year(year,1)=Discount_factor(year,1)*Income_year(i,1); 
            end
        end       
        NPV(year,1)=Discount_income_year(year,1)-Discount_costs(year,1);        
        if year>1 % actually starting year -1          
            Cumulative_NPV(year,1)=Cumulative_NPV(year-1)+NPV(year);
        else % actaully year -2
            Cumulative_NPV(1,1)=0;
        end
    end
    NPV_final(i,1)=Cumulative_NPV(end,1);
    Profitability_index(i,1)=sum(Discount_income_year)/sum(Discount_costs);
    years=[-2:(Lifespan-3)];
    
    % Find the year when NPV becomes above 0
    [r_NPV,~]=find(Cumulative_NPV>0);
    if r_NPV>0
        year_NPV_positive(i,1)=r_NPV(1,1)-3;
    else
        year_NPV_positive(i,1)=9999;
    end
    
% If the NPV after 25 years is less than 0, then add 10 more years 
% of functioning to see if the NPV increaases   
    if NPV_final(i,1)<0
        Discount_factor=[];
        Discount_costs=[];
        Discount_income_year=[];
        Cumulative_NPV=[];
        NPV=[];
        Lifespan=35+3;
        Capacity_factor=40/100;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Wholesale_price=70; %price in pounds per MWh
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Discount_rate=4/100; %discount rate - chosen
        CAPEX(i,1)=C_tot(i,1); %capex cost in year -1 (considered 1 in the for loop)
        OPEX(i,1)=3/100*C_tot(i,1); %operational and maitenance costs assumed to be 3% of the total installed costs
        Income_year(i,1)=(Intake(i,3)/1000*24*365)*Capacity_factor*Wholesale_price;%yield per year considering 
        % capacity factor of 40% and the wholesale electricity price 45 pounds/MWh
        for year=1:(Lifespan)
            Discount_factor(year,1)=1/(1+Discount_rate)^(year-3);
            if year==1 %actualy year -2
                Discount_costs(year,1)=0;
                Discount_income_year(year,1)=0;
            else
                if year==2 %actaully year -1
                    Discount_costs(year,1)=CAPEX(i,1)*Discount_factor(year,1);
                    Discount_income_year(year,1)=0;
                else %actually starting with year 0
                    Discount_costs(year,1)=OPEX(i,1)*Discount_factor(year,1);
                    Discount_income_year(year,1)=Discount_factor(year,1)*Income_year(i,1); 
                end
            end       
            NPV(year,1)=Discount_income_year(year,1)-Discount_costs(year,1);        
            if year>1 % actually starting year -1          
                Cumulative_NPV(year,1)=Cumulative_NPV(year-1)+NPV(year);
            else % actually year -2
                Cumulative_NPV(1,1)=0;
            end
        end
        NPV_final(i,1)=Cumulative_NPV(end,1);
        Profitability_index(i,1)=sum(Discount_income_year)/sum(Discount_costs);
        years=[-2:(Lifespan-3)];
        
        % Find the year when NPV becomes above 0
        [r_NPV,~]=find(Cumulative_NPV>0);
        if r_NPV>0
            year_NPV_positive(i,1)=r_NPV(1,1)-3;
        else
            year_NPV_positive(i,1)=9999;
        end
    end
     
% If the NPV after 35 years is less than 0 for pico or micro hydro,
% consider the initial costs without the power grid connection

    if NPV_final(i,1)<0
        if P_min>=0.1 && P_max<=100
        Connection_bool(i,1)=0;
        %Transmission lines
        C9(i,1)=0; %considering 1 km as the distance
        C12(i,1)=0.04*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)); %development costs
        C13(i,1)=0.032*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)+C12(i,1)); %feasibility study cost
        C14(i,1)=0.25*interest*Flow(i,1)^0.35*1.1*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)+C12(i,1))+0.1*(C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)); %miscellaneous costs
        C_tot(i,1)=C1(i,1)+C2(i,1)+C3(i,1)+C4(i,1)+C5(i,1)+C6(i,1)+C7(i,1)+C8(i,1)+C9(i,1)+C10(i,1)+C11(i,1)+C12(i,1)+C13(i,1)+C14(i,1);
        Cost_per_kW(i,1)=C_tot(i,1)/Intake(i,3); %cost per kW intalled
    
        %Get the percetage of costs
        C1_perc(i,1)=C1(i,1)/C_tot(i,1)*100;
        C2_perc(i,1)=C2(i,1)/C_tot(i,1)*100;
        C3_perc(i,1)=C3(i,1)/C_tot(i,1)*100;
        C4_perc(i,1)=C4(i,1)/C_tot(i,1)*100;
        C5_perc(i,1)=C5(i,1)/C_tot(i,1)*100;
        C6_perc(i,1)=C6(i,1)/C_tot(i,1)*100;
        C7_perc(i,1)=C7(i,1)/C_tot(i,1)*100;
        C8_perc(i,1)=C8(i,1)/C_tot(i,1)*100;
        C9_perc(i,1)=C9(i,1)/C_tot(i,1)*100;
        C10_perc(i,1)=C10(i,1)/C_tot(i,1)*100;
        C11_perc(i,1)=C11(i,1)/C_tot(i,1)*100;
        C12_perc(i,1)=C12(i,1)/C_tot(i,1)*100;
        C13_perc(i,1)=C13(i,1)/C_tot(i,1)*100;
        C14_perc(i,1)=C14(i,1)/C_tot(i,1)*100;
        
        %Electricity generated on an average year
        Discount_factor=[];
        Discount_costs=[];
        Discount_income_year=[];
        Cumulative_NPV=[];
        NPV=[];
        Lifespan=35+3;
        Capacity_factor=40/100;
        Wholesale_price=45; %price in pounds per MWh
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Modify as required
        Discount_rate=4/100; %discount rate - chosen
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CAPEX(i,1)=C_tot(i,1); %capex cost in year -1 (considered 1 in the for loop)
        OPEX(i,1)=3/100*C_tot(i,1); %operational and maitenance costs assumed to be 3% of the total installed costs
        Income_year(i,1)=(Intake(i,3)/1000*24*365)*Capacity_factor*Wholesale_price;%yield per year considering 
        % capacity factor of 40% and the wholesale electricity price 45 pounds/MWh
        for year=1:(Lifespan)
            Discount_factor(year,1)=1/(1+Discount_rate)^(year-3);
            if year==1 %actualy year -2
                Discount_costs(year,1)=0;
                Discount_income_year(year,1)=0;
            else
                if year==2 %actaully year -1
                    Discount_costs(year,1)=CAPEX(i,1)*Discount_factor(year,1);
                    Discount_income_year(year,1)=0;
                else %actually starting with year 0
                    Discount_costs(year,1)=OPEX(i,1)*Discount_factor(year,1);
                    Discount_income_year(year,1)=Discount_factor(year,1)*Income_year(i,1); 
                end
            end       
        NPV(year,1)=Discount_income_year(year,1)-Discount_costs(year,1);        
        if year>1 % actually starting year -1          
            Cumulative_NPV(year,1)=Cumulative_NPV(year-1)+NPV(year);
        else % actaully year -2
            Cumulative_NPV(1,1)=0;
        end
        end
        NPV_final(i,1)=Cumulative_NPV(end,1);
        Profitability_index(i,1)=sum(Discount_income_year)/sum(Discount_costs);
        years=[-2:(Lifespan-3)];
    
        % Find the year when NPV becomes above 0
        [r_NPV,~]=find(Cumulative_NPV>0);
        if r_NPV>0
            year_NPV_positive(i,1)=r_NPV(1,1)-3;
        else
            year_NPV_positive(i,1)=9999;
        end
        end
    end
    
%% Figure NPV
if isnan(NPV_final(i,1))==0
    if print_NPV_fig==1
        clf %clears figure information
        set(gca,'FontSize',20);
        fontsize(gca,14,'pixels');
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        numb=num2str(PowerHouse(i,2));
        save_name=strcat(folder_location,'\CumulativeNPV_',numb,'.png');
        title_name_NPV=strcat("Cumulative Net Present Value for RoR scheme no. ",numb," in Area no. ",Area(5:8));
        title(title_name_NPV, 'FontSize', 20);
        grid on;
        grid minor;        
        xlabel("Years", 'FontSize', 20);
        ylabel("Cumulative NPV [pounds]", 'FontSize', 20);
        hold on
        %plot(x_axis(:,1),eq(:,1),'Color','[0.3 0 1]','LineWidth',2,'LineStyle','-','Marker','none');
        plot(years,Cumulative_NPV,'-x','LineWidth',2);
        hold on
        if exist(save_name, 'file')==2
            delete(save_name);
        end
            saveas(gcf,save_name);
    end 
end

end

%% Remove locations which cannot be realised

[Limitr,~]=size(PowerHouse);

if Limitr>0
r_cost=[];
[r_cost,~]=find(isnan(C_tot));

for i=1:length(Intake)
   Ps(r_cost,:)=0;
   Intake(r_cost,:)=0;
   PowerHouse(r_cost,:)=0;
   diam(r_cost,:)=5;
   Flow(r_cost)=0;
   NetHead(r_cost)=0;
   C_tot(r_cost,1)=0;
   NPV_final(r_cost,1)=0;
   Cost_per_kW(r_cost,1)=0;
   
   C1(r_cost,1)=9999;
   C1_perc(r_cost,1)=9999;
   C2(r_cost,1)=9999;
   C2_perc(r_cost,1)=9999;
   C3(r_cost,1)=9999;
   C3_perc(r_cost,1)=9999;
   C4(r_cost,1)=9999;
   C4_perc(r_cost,1)=9999;
   C5(r_cost,1)=9999;
   C5_perc(r_cost,1)=9999;
   C6(r_cost,1)=9999;
   C6_perc(r_cost,1)=9999;
   C7(r_cost,1)=9999;
   C7_perc(r_cost,1)=9999;
   C8(r_cost,1)=9999;
   C8_perc(r_cost,1)=9999;
   C9(r_cost,1)=9999;
   C9_perc(r_cost,1)=9999;
   C10(r_cost,1)=9999;
   C10_perc(r_cost,1)=9999;
   C11(r_cost,1)=9999;
   C11_perc(r_cost,1)=9999;
   C12(r_cost,1)=9999;
   C12_perc(r_cost,1)=9999;
   C13(r_cost,1)=9999;
   C13_perc(r_cost,1)=9999;
   C14(r_cost,1)=9999;
   C14_perc(r_cost,1)=9999;   
end

for i=1:length(r_cost)
      River{r_cost(i,1),1}=[];
      Pmax2{r_cost(i,1),2}=[];
      Q_perc{r_cost(i,1),1}=[];
end

River=River(~any(cellfun('isempty',River),2),:);
Pmax2=Pmax2(~any(cellfun('isempty',Pmax2),2),:);
Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

diam(all(diam==5,2),:)=[];
Ps(all(Ps==0,2),:)=[];

Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

C_tot(all(C_tot==0,2),:)=[];
NPV_final(all(NPV_final==0,2),:)=[];
Cost_per_kW(all(Cost_per_kW==0,2),:)=[];

C1(all(C1==9999,2),:)=[];
C1_perc(all(C1_perc==9999,2),:)=[];
C2(all(C2==9999,2),:)=[];
C2_perc(all(C2_perc==9999,2),:)=[];
C3(all(C3==9999,2),:)=[];
C3_perc(all(C3_perc==9999,2),:)=[];
C4(all(C4==9999,2),:)=[];
C4_perc(all(C4_perc==9999,2),:)=[];
C5(all(C5==9999,2),:)=[];
C5_perc(all(C5_perc==9999,2),:)=[];
C6(all(C6==9999,2),:)=[];
C6_perc(all(C6_perc==9999,2),:)=[];
C7(all(C7==9999,2),:)=[];
C7_perc(all(C7_perc==9999,2),:)=[];
C8(all(C8==9999,2),:)=[];
C8_perc(all(C8_perc==9999,2),:)=[];
C9(all(C9==9999,2),:)=[];
C9_perc(all(C9_perc==9999,2),:)=[];
C10(all(C10==9999,2),:)=[];
C10_perc(all(C10_perc==9999,2),:)=[];
C11(all(C11==9999,2),:)=[];
C11_perc(all(C11_perc==9999,2),:)=[];
C12(all(C12==9999,2),:)=[];
C12_perc(all(C12_perc==9999,2),:)=[];
C13(all(C13==9999,2),:)=[];
C13_perc(all(C13_perc==9999,2),:)=[];
C14(all(C14==9999,2),:)=[];
C14_perc(all(C14_perc==9999,2),:)=[];
end

%% Add the cost and NPV to the database

%Write the turbine type as string

Intake_cell=cell(length(Intake),1);
PowerHouse_cell=cell(length(PowerHouse),1);

[Limitr,~]=size(PowerHouse);

counts=Limitr;
Intake(:,1)=1:counts;
PowerHouse(:,1)=1:counts;

Q_str=string(Q_perc);

for i=1:Limitr
    Intake_cell{i,1}=Intake(i,1);
    Intake_cell{i,2}=Intake(i,2);
    Intake_cell{i,3}=Intake(i,3);
    Intake_cell{i,4}=Intake(i,4);
    Intake_cell{i,5}=Intake(i,5);
    Intake_cell{i,7}=Intake(i,7);
    Intake_cell{i,8}=Intake(i,8);
    Intake_cell{i,9}=Intake(i,9);
    Intake_cell{i,10}=Intake(i,10);
    Intake_cell{i,11}=Q_str(i,1);
   if Intake(i,6)==1
       Intake_cell{i,6}="Pelton";
   end
   if Intake(i,6)==2
      Intake_cell{i,6}="Francis";
   end
   if Intake(i,6)==3
      Intake_cell{i,6}="Kaplan";
   end
   if Intake(i,6)==4
      Intake_cell{i,6}="Archimedean Screw";
   end   
   if Intake(i,6)==5
      Intake_cell{i,6}="Cross-flow";
   end
   if Intake(i,6)==-9999
      Intake_cell{i,6}="UNK";
   end
   if Connection_bool==1
       Intake_cell{i,16}='Yes';
   else
       Intake_cell{i,16}='No';
   end
end

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    PowerHouse_cell{i,1}=PowerHouse(i,1);
    PowerHouse_cell{i,2}=PowerHouse(i,2);
    PowerHouse_cell{i,3}=PowerHouse(i,3);
    PowerHouse_cell{i,4}=PowerHouse(i,4);
    PowerHouse_cell{i,5}=PowerHouse(i,5);
    PowerHouse_cell{i,7}=PowerHouse(i,7);
    PowerHouse_cell{i,8}=PowerHouse(i,8);
    PowerHouse_cell{i,9}=PowerHouse(i,9);
    PowerHouse_cell{i,10}=PowerHouse(i,10);
    PowerHouse_cell{i,11}=Q_str(i,1);
   if PowerHouse(i,6)==1
       PowerHouse_cell{i,6}="Pelton";
   end
   if PowerHouse(i,6)==2
      PowerHouse_cell{i,6}="Francis";
   end
   if PowerHouse(i,6)==3
      PowerHouse_cell{i,6}="Kaplan";
   end
   if PowerHouse(i,6)==4
      PowerHouse_cell{i,6}="Archimedean Screw";
   end   
   if PowerHouse(i,6)==5
      PowerHouse_cell{i,6}="Cross-flow";
   end
   if PowerHouse(i,6)==-9999
      PowerHouse_cell{i,6}="UNK";
   end
   if Connection_bool==1
       PowerHouse_cell{i,16}='Yes';
   else
       PowerHouse_cell{i,16}='No';
   end
end

%In case the matrices are empty put NaN in them so the write to excel
%doesn't give an error

bool_Intake=isempty(Intake_cell);
bool_PowerHouse=isempty(PowerHouse_cell);

if bool_Intake==1
    Intake_cell='NaN';
    PowerHouse_cell='NaN';
end

Limitr=[];
[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    Intake_cell{i,12}=C_tot(i,1)/10^6;
    PowerHouse_cell{i,12}=C_tot(i,1)/10^6;
    Intake_cell{i,13}=NPV_final(i,1)/10^6;
    PowerHouse_cell{i,13}=NPV_final(i,1)/10^6;
    Intake_cell{i,14}=Cost_per_kW(i,1);
    PowerHouse_cell{i,14}=Cost_per_kW(i,1);
    Intake_cell{i,15}=year_NPV_positive(i,1);
    PowerHouse_cell{i,15}=year_NPV_positive(i,1);
end

% % Write a string matrix with the cost names
% 
% Costs_cell={'Generator cost','Turbine cost','Turbine and generator installation cost','Penstock cost','Penstock installation cost','Civil structures cost','Substation and transformer cost','Substation and transformer installation cost','Transmission line connection cost','Access road cost','Engineering design cost','Development cost','Feasibility study cost','Miscellaneous cost','Total cost'};
% Costs_cell=Costs_cell';
% 
% % Write a matrix containing all cost names and values
% 
% Costs_final_string(:,1)=Costs_cell;
% Costs_final=[];
% 
% Limitr=[];
% [Limitr,~]=size(PowerHouse);
% 
% for i=1:Limitr
% Costs_final(1,1)=C1(i,1);
% Costs_final(2,1)=C2(i,1);
% Costs_final(3,1)=C3(i,1);
% Costs_final(4,1)=C4(i,1);
% Costs_final(5,1)=C5(i,1);
% Costs_final(6,1)=C6(i,1);
% Costs_final(7,1)=C7(i,1);
% Costs_final(8,1)=C8(i,1);
% Costs_final(9,1)=C9(i,1);
% Costs_final(10,1)=C10(i,1);
% Costs_final(11,1)=C11(i,1);
% Costs_final(12,1)=C12(i,1);
% Costs_final(13,1)=C13(i,1);
% Costs_final(14,1)=C14(i,1);
% Costs_final(15,1)=C_tot(i,1);
% 
% Costs_final(1,2)=C1_perc(i,1);
% Costs_final(2,2)=C2_perc(i,1);
% Costs_final(3,2)=C3_perc(i,1);
% Costs_final(4,2)=C4_perc(i,1);
% Costs_final(5,2)=C5_perc(i,1);
% Costs_final(6,2)=C6_perc(i,1);
% Costs_final(7,2)=C7_perc(i,1);
% Costs_final(8,2)=C8_perc(i,1);
% Costs_final(9,2)=C9_perc(i,1);
% Costs_final(10,2)=C10_perc(i,1);
% Costs_final(11,2)=C11_perc(i,1);
% Costs_final(12,2)=C12_perc(i,1);
% Costs_final(13,2)=C13_perc(i,1);
% Costs_final(14,2)=C14_perc(i,1);
% Costs_final(15,2)=C1_perc(i,1)+C2_perc(i,1)+C3_perc(i,1)+C4_perc(i,1)+C5_perc(i,1)+C6_perc(i,1)+C7_perc(i,1)+C8_perc(i,1)+C9_perc(i,1)+C10_perc(i,1)+C11_perc(i,1)+C12_perc(i,1)+C13_perc(i,1)+C14_perc(i,1);
% 
% no=num2str(i);
% name3=strcat(folder_location,'\Costs_',no,'_');
% 
% extension='.xls';
% 
% CostsName=[name3 Area extension];
% 
% if exist(CostsName, 'file')==2
%   delete(CostsName);
% end
% 
% col_Header_Costs={'Cost name','Cost (pounds)','Percentage of total cost'};
% 
% xlswrite(CostsName,Costs_final_string,'Sheet1','A2');
% xlswrite(CostsName,Costs_final,'Sheet1','B2');
% xlswrite(CostsName,col_Header_Costs,'Sheet1','A1');
%end

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

col_Header_Intake={'ID','Intake_ID','Power[kW]','Intake_X','Intake_Y','Turbine_type','Flow[m^3/s]','Net head[m]','Peak efficiency','Diameter [m]','Flow type','Initial total cost [m. pounds]','NPV [m. pounds]','Installed cost per kW [pounds/kW]','Year NPV positive','Connection to the grid considered'};
col_Header_PowerHouse={'ID','PowerHouse_ID','Power[kW]','PowerHouse_X','PowerHouse_Y','Turbine_type','Flow[m^3/s]','Net head[m]','Peak efficiency','Diameter [m]','Flow type','Initial total cost [m. pounds]','NPV [m. pounds]','Installed cost per kW [pounds/kW]','Year NPV positive','Connection to the grid considered'};

xlswrite(IntakeName,Intake_cell,'Sheet1','A2');
xlswrite(IntakeName,col_Header_Intake,'Sheet1','A1');

xlswrite(PowerHouseName,PowerHouse_cell,'Sheet1','A2');
xlswrite(PowerHouseName,col_Header_PowerHouse,'Sheet1','A1');

end