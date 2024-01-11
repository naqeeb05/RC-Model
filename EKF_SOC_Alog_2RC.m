clear
close all
% clc

n=3;

%% load data
load('rc_param')

switch n
    case 1
        load('TerraE_INR21700_50E_intermediate_2')    % data from LYES
        meas_data = TestData.Cell_1.T2.WLTP;
    case 3
        profile = load('_Road_Course_Profile_Battery_Level.txt');
    case 4
        profile = load('_Stacking_Cycle_Profile_Battery_Level.txt');
end

if n==1
    rawtime = meas_data.t;
    rawvoltage = meas_data.U_cell;
    rawcurrent = -1.*meas_data.I_cell;
else
    rawtime = profile(1:end,1);
    rawvoltage = profile(1:end,2)./28;
    rawcurrent = -profile(1:end,3)./60;
end

rawtime = rawtime-rawtime(1);

% up sample to 1s
time = 0:1:round(rawtime(end));
time = time';
Voltage = interp1(rawtime, rawvoltage, time);
Current = interp1(rawtime, rawcurrent, time);

% create Temp & SoC grid
[TempGrid,SoCGrid] = meshgrid(Temp_vec, SoC_vec);

% set initial conditions
dt = 1;
% Temp = data.T2.T_start;
Temp = meas_data.T_chamber(1);
% start_soc = interp2(TempGrid, SoCGrid, Em_table, Temp, Voltage(1), 'pchip', 'extrap');
start_soc = interp1(Em_vec(:,2), SoC_vec, Voltage(1), 'pchip', 'extrap');
x_hat = [0.7; 0; 0];    % zk, i_R1, i_R2


load('tuning_param_2RC')
load('RMSE_2RC')
m = find(RMSE==min(RMSE));

SigmaW = tuning_param(m,1);            % process noise covariance
SigmaV = tuning_param(m,2);            % sensor noise covariance
SigmaX = diag([tuning_param(m,3), tuning_param(m,4), tuning_param(m,5)]);





% plotting
% SoC CC vs EKF
figure()
set(gca, 'fontsize', 12)
hold on
grid on
box on
plot(Time/60, SoC_CC*100, Time/60, SoC_EKF*100, 'linewidth', 1)
plot([Time/60 Time/60], [100*(SoC_EKF+zkbnd) 100*(SoC_EKF-zkbnd)], '--', 'color', 'k')
xlabel('Time / min')
ylabel('SoC / %')
ylim([0 100])
title('SoC Estimation using EKF')
legend('CC', 'EKF', 'Bounds')
fprintf('RMSE=%g%%\n', sqrt(mean((100*(SoC_CC-SoC_EKF)).^2)));

% SoC error
figure()
set(gca, 'fontsize', 12)
hold on
grid on
box on
plot(Time/60, 100*(SoC_CC-SoC_EKF), 'linewidth', 1)
plot([Time/60 Time/60], [100*zkbnd -100*zkbnd], 'linestyle', '--', 'color', 'r')
% ylim([0 100])
% plot(Time/60, error_vec)
xlabel('Time / min')
ylabel('SoC Error / %')
title('SoC Estimation error using EKF')
legend('Estimation error', 'Bounds')


% SoC error
figure()
set(gca, 'fontsize', 12)
hold on
grid on
box on
plot(Time/60, Voltage)
plot(Time/60, V_EKF)
ylim([2.5 4.2])
% plot(TimeInterp/60, error_vec)
ylabel('Voltage / V')
xlabel('Time / min')
legend('V_{true}', 'V_{EKF}')











