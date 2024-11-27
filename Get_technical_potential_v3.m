function get_scheme_locations = get_locations(Area,P_min,P_max,folder_location);

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

[r,~]=find(emptyCells==1);

for i=1:length(Intake)
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

%Create matrices that contain the head and the flow

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
   V(i,1)=4.*Flow(i,1)/(3.14*diam(i,1));
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

ep=zeros(length(PowerHouse),1); %peak efficiency
Qp=zeros(length(PowerHouse),1); %peak efficiency flow
nq=zeros(length(PowerHouse),1); %specific speed
enq_=zeros(length(PowerHouse),1); %adjustment to turbine peak efficiency
ed_=zeros(length(PowerHouse),1); %adjustment to turbine peak efficiency
d=zeros(length(PowerHouse),1); %runner size

[Limitr,~]=size(PowerHouse);

%Francis turbine
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
    Qp(i,1)=0.65*Flow(i,1)*nq(i,1).^0.05;
    end
end

%Kaplan turbine
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
    end
end

%Impulse turbines peak efficiency calculation (Pelton=1)

jets=1; %number of jets for Pelton
n_rot=zeros(length(PowerHouse),1); %rotational speed

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
    end   
end
jets=1; %number of jets for Pelton needs to be reset for future calculations

%Impulse turbines peak efficiency calculation (Cross-flow=5)

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    if PowerHouse(i,6)==5
        ep(i,1)=0.816;
        Qp(i,1)=0.7*Flow(i,1);
    end   
end

%Archimedean Screw=4 turbine peak efficiency calculation

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    if PowerHouse(i,6)==4
        ep(i,1)=0.828;
        Qp(i,1)=0.78*Flow(i,1);
    end   
end

%% Remove schemes with less than 0.7 peak efficiency or don't have a
%%suitable turbine

emptyCells=zeros(length(PowerHouse),1);

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
    emptyCells(i,:)=ep(i,1)>0.7;
end

for i=1:length(emptyCells)
    if emptyCells(i,1)==0
        Q_perc{i,1}=[];
    end
end

Q_perc=Q_perc(~any(cellfun('isempty',Q_perc),2),:);

[r,~]=find(emptyCells==0);

Intake(r,:)=0;
PowerHouse(r,:)=0;
Flow(r,:)=0;
NetHead(r,:)=0;
ep(r,:)=0;
diam(r,:)=0;

Intake(all(~Intake,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];
Flow(all(Flow==0,2),:)=[];
NetHead(all(NetHead==0,2),:)=[];
ep(all(ep==0,2),:)=[];
Qp(all(Qp==0,2),:)=[];
diam(all(diam==0,2),:)=[];

Intake(:,7)=Flow(:,1);
Intake(:,8)=NetHead(:,1);
Intake(:,9)=ep(:,1);
Intake(:,10)=diam(:,1);

PowerHouse(:,7)=Flow(:,1);
PowerHouse(:,8)=NetHead(:,1);
PowerHouse(:,9)=ep(:,1);
PowerHouse(:,10)=diam(:,1);

%% Recalculate power considering the peak efficiency and flow

[Limitr,Limitc]=size(PowerHouse);

for i=1:Limitr
    if ep(i,1)~=0
       Intake(i,3)=ro*g*Qp(i,1)*NetHead(i,1)*ep(i,1)/1000; %power in kW 
       PowerHouse(i,3)=round(abs(ro*g*Qp(i,1)*NetHead(i,1)*ep(i,1))/1000); %power in kW
    end
end

%Eliminate the combinations that would generate less than P_min or more
%than P_max

%Remove combinations that would generate less than P_min

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

[r,~]=find(emptyCells==1);

for i=1:Limitr
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

emptyCells=[];

%Remove combinations that would generate more than P_max

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

[r,~]=find(emptyCells==1);

for i=1:Limitr
   Intake(r,:)=0;
   PowerHouse(r,:)=0;
end

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

%% Recalculate the diameter in case the flow has changed to make the scheme
%%fit a turbine

[Limitr,~]=size(PowerHouse);

for i=1:Limitr
   ID_Intake=Pmax2{i,1}(1,1);
   ID_PowerHouse=Pmax2{i,1}(2,1); 
   [rI,cI]=find(River{i,1}(:,1)==ID_Intake);
   [rPH,cPH]=find(River{i,1}(:,1)==ID_PowerHouse);
   L(i,1)= sqrt((River{i,1}(rI,23)-River{i,1}(rPH,23)).^2+(River{i,1}(rI,24)-River{i,1}(rPH,24)).^2+(River{i,1}(rI,25)-River{i,1}(rPH,25)).^2); %distance between the Intake and the PowerHouse
   %diam(i,1)=2.83*((L(i,1).^2*Flow(i,1).^2)/Head(i,1)).^0.1876; %pipe diameter
   %diam2(i,1)=nthroot((257.5*(L(i,1)*epsilon^2*Flow(i,1).^2/Head(i,1))).^3,16);
   diam(i,1)=2.69*((epsilon^2*Flow(i,1).^2*L(i,1))/Head(i,1)).^0.1875; %pipe diameter
 
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

Intake(all(Intake==0,2),:)=[];
PowerHouse(all(PowerHouse==0,2),:)=[];

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
end

[Limitr,Limitc]=size(PowerHouse);

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
end

%In case the matrices are empty put NaN in them so the write to excel
%doesn't give an error

bool_Intake=isempty(Intake_cell);
bool_PowerHouse=isempty(PowerHouse_cell);

if bool_Intake==1
    Intake_cell='NaN';
    PowerHouse_cell='NaN';
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

col_Header_Intake={'ID','Intake_ID','Power[kW]','Intake_X','Intake_Y','Turbine_type','Flow[m^3/s]','Net head[m]','Peak efficiency','Diameter [m]','Flow type'};
col_Header_PowerHouse={'ID','PowerHouse_ID','Power[kW]','PowerHouse_X','PowerHouse_Y','Turbine_type','Flow[m^3/s]','Net head[m]','Peak efficiency','Diameter [m]','Flow type'};

xlswrite(IntakeName,Intake_cell,'Sheet1','A2');
xlswrite(IntakeName,col_Header_Intake,'Sheet1','A1');

xlswrite(PowerHouseName,PowerHouse_cell,'Sheet1','A2');
xlswrite(PowerHouseName,col_Header_PowerHouse,'Sheet1','A1');

end