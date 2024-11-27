%% Search for maximum technical potential Version 2

for i=1:97
    i
   no=num2str(i); 
   folder_name=strcat('Area00',no);
   folder_location=strcat('H:\01.PhD\003.Chapter2\ROR_Search\Technical_potential_v2\',folder_name);
   file_number=strcat('Area00',no);
   
   Get_technical_potential_v2(strcat(file_number,'_pico'),0.1,5,folder_location);
   Get_technical_potential_v2(strcat(file_number,'_micro'),5,100,folder_location);
   Get_technical_potential_v2(strcat(file_number,'_mini'),100,1000,folder_location);
   Get_technical_potential_v2(strcat(file_number,'_small'),1000,10000,folder_location);
   
end

%% Search for maximum technical potential Version 3

for i=1:10
    i
   no=num2str(i); 
   folder_name=strcat('Area00',no);
   folder_location=strcat('H:\01.PhD\003.Chapter2\ROR_Search\Technical_potential_v3\',folder_name);
   file_number=strcat('Area00',no);
   
   Get_technical_potential_v3(strcat(file_number,'_pico'),0.1,5,folder_location);
   Get_technical_potential_v3(strcat(file_number,'_micro'),5,100,folder_location);
   Get_technical_potential_v3(strcat(file_number,'_mini'),100,1000,folder_location);
   Get_technical_potential_v3(strcat(file_number,'_small'),1000,10000,folder_location);
   
end

%% Search for maximum hydrological potential

for i=33:97
    i
   no=num2str(i); 
   folder_name=strcat('Area00',no);
   folder_location=strcat('H:\01.PhD\003.Chapter2\ROR_Search\Hydrological_potential\',folder_name);
   file_number=strcat('Area00',no);
   
   Get_hydrological_potential(strcat(file_number,'_pico_hp'),0.1,5,folder_location);
   Get_hydrological_potential(strcat(file_number,'_micro_hp'),5,100,folder_location);
   Get_hydrological_potential(strcat(file_number,'_mini_hp'),100,1000,folder_location);
   Get_hydrological_potential(strcat(file_number,'_small_hp'),1000,10000,folder_location);
   
end

%% Search for turbine 1kW

for i=1
    i
   no=num2str(i); 
   folder_name=strcat('Area00',no);
   folder_location=strcat('H:\01.PhD\003.Chapter2\ROR_Search\Technical_potential\',folder_name);
   file_number=strcat('Area00',no);

   Get_locations_v4(strcat(file_number,'_1kW'),1,'yeet',0.85,folder_location);
   
end

%% Search for maximum economical potential Version 1

for i=27
   i
   no=num2str(i); 
   folder_name=strcat('Area00',no);
   folder_location=strcat('H:\01.PhD\003.Chapter2\ROR_Search\Financial_potential_v1\',folder_name);
   file_number=strcat('Area00',no);
   
   Get_Financial_potential_v1(strcat(file_number,'_pico'),0.1,5,folder_location);
   Get_Financial_potential_v1(strcat(file_number,'_micro'),5,100,folder_location);
   Get_Financial_potential_v1(strcat(file_number,'_mini'),100,1000,folder_location);
   Get_Financial_potential_v1(strcat(file_number,'_small'),1000,10000,folder_location);
   
end

%% Search for maximum economical potential Version 2 - green, amber, red system

for i=1:97
   i
   no=num2str(i); 
   folder_name=strcat('Area00',no);
   folder_location=strcat('H:\01.PhD\003.Chapter2\ROR_Search\Financial_potential_v2\',folder_name);
   file_number=strcat('Area00',no);
   
   Get_Financial_potential_v2(strcat(file_number,'_pico'),0.1,5,folder_location);
   Get_Financial_potential_v2(strcat(file_number,'_micro'),5,100,folder_location);
   Get_Financial_potential_v2(strcat(file_number,'_mini'),100,1000,folder_location);
   Get_Financial_potential_v2(strcat(file_number,'_small'),1000,10000,folder_location);
   
end


