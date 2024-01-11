clc;clear;close all
cd('K:\Battery management systems\DATA THESIS\CALACE')
load('calace_rc_parameter.mat')
%% User selection
% cell selection
% ---------------
% test selection    -> test name
% ---------------
% 1                 -> DST at 50% at 25°C
% 2                 -> DST at 80% at 25°C
% 3                 -> FUDS at 80% at 25°C
% 4                 -> FUDS at 50% at 25°C
% 5                 -> US06 at 50% at 25°C
% 6                 -> US06 at 80% at 25°C
% 7                 -> BJDST at 80% 25°C

test_selection =3;
folder='B:\data_excel';
cd(folder);list = dir('*.xls');names={list.name};

switch test_selection
    case 1
        table = readtable(names{test_selection},'sheet',2);
        time=table.Test_Time_s_;voltage=table.Voltage_V_;current=table.Current_A_;
        range=1:length(time);range=range(2810:end);
        time=time(range);voltage=voltage(range);current=current(range);
        time=time-time(1);
        init_soc=0.5; init_soc_cc=0.5; eta=10^-1;
        CapacityAh=2; dt=1; x_initial=[init_soc; 0];L_gain_initial=[0.08495; 0.1385];
        dx_dLgain_initial=zeros(2,2);
        voltage=timeseries(voltage,time);current=timeseries(current,time);
        sim('ALU_model_1','StartTime','time(1)','StopTime','time(end)','FixedStep','dt');
    
    case 2
        table = readtable(names{test_selection},'sheet',2);
        time=table.Test_Time_s_;voltage=table.Voltage_V_;current=table.Current_A_;
        range=1:length(time);range=range(1914:end);
        time=time(range);voltage=voltage(range);current=current(range);
        time=time-time(1);
        init_soc=0.8; init_soc_cc=0.8; eta=10^-1;
        CapacityAh=2; dt=1; x_initial=[init_soc; 0];L_gain_initial=[0.08495; 0.1385];
        dx_dLgain_initial=zeros(2,2);
        voltage=timeseries(voltage,time);current=timeseries(current,time);
        sim('ALU_model_1','StartTime','time(1)','StopTime','time(end)','FixedStep','dt');
    
    case 3
        table = readtable(names{test_selection},'sheet',2);
        time=table.Test_Time_s_;voltage=table.Voltage_V_;current=table.Current_A_;
        range=1:length(time);range=range(2590:end);
        time=time(range);voltage=voltage(range);current=current(range);
        time=time-time(1);
        init_soc=0.5; init_soc_cc=0.8; eta=10^-1;
        CapacityAh=2; dt=1; x_initial=[init_soc; 0];L_gain_initial=[0.08495; 0.1385];
        dx_dLgain_initial=zeros(2,2);
        voltage=timeseries(voltage,time);current=timeseries(current,time);
        sim('ALU_model_1','StartTime','time(1)','StopTime','time(end)','FixedStep','dt');
            
    case 4
        table = readtable(names{test_selection},'sheet',2);
        time=table.Test_Time_s_;voltage=table.Voltage_V_;current=table.Current_A_;
        range=1:length(time);range=range(2324:end);
        time=time(range);voltage=voltage(range);current=current(range);
        time=time-time(1);
        init_soc=0.5; init_soc_cc=0.5; eta=10^-1;
        CapacityAh=2; dt=1; x_initial=[init_soc; 0];L_gain_initial=[0.08495; 0.1385];
        dx_dLgain_initial=zeros(2,2);
        voltage=timeseries(voltage,time);current=timeseries(current,time);
        sim('ALU_model_1','StartTime','time(1)','StopTime','time(end)','FixedStep','dt');
            
    case 5
        table = readtable(names{test_selection},'sheet',2);
        time=table.Test_Time_s_;voltage=table.Voltage_V_;current=table.Current_A_;
        range=1:length(time);range=range(1433:end);
        time=time(range);voltage=voltage(range);current=current(range);
        time=time-time(1);
        init_soc=0.53; init_soc_cc=0.53; eta=10^-1;
        CapacityAh=2; dt=1; x_initial=[init_soc; 0];L_gain_initial=[0.08495; 0.1385];
        dx_dLgain_initial=zeros(2,2);
        voltage=timeseries(voltage,time);current=timeseries(current,time);
        sim('ALU_model_1','StartTime','time(1)','StopTime','time(end)','FixedStep','dt');

    case 6
        table = readtable(names{test_selection},'sheet',2);
        time=table.Test_Time_s_;voltage=table.Voltage_V_;current=table.Current_A_;
        range=1:length(time);range=range(2324:end);
        time=time(range);voltage=voltage(range);current=current(range);
        time=time-time(1);
        init_soc=0.9; init_soc_cc=0.9; eta=10^-1;
        CapacityAh=5; dt=1; x_initial=[init_soc; 0];L_gain_initial=[10^-9; 10^-9];
        dx_dLgain_initial=zeros(2,2);
        voltage=timeseries(voltage,time);current=timeseries(current,time);
        sim('ALU_model_1','StartTime','time(1)','StopTime','time(end)','FixedStep','dt');

    case 7
        table = readtable(names{test_selection},'sheet',2);
        time=table.Test_Time_s_;voltage=table.Voltage_V_;current=table.Current_A_;
        range=1:length(time);range=range(2324:end);
        time=time(range);voltage=voltage(range);current=current(range);
        time=time-time(1);
        init_soc=0.745; init_soc_cc=0.745; eta=10^-1;
        CapacityAh=2; dt=1; x_initial=[init_soc; 0];L_gain_initial=[0.08495; 0.1385];
        dx_dLgain_initial=zeros(2,2);
        voltage=timeseries(voltage,time);current=timeseries(current,time);
        sim('ALU_model_1','StartTime','time(1)','StopTime','time(end)','FixedStep','dt');
end
