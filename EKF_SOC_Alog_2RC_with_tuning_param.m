clear
close all
clc

%% load data
tic
n=1;

%load data
load('rc_param')
load('batemo_data')
load('TerraE_INR21700_50E_intermediate_2')

switch n
    
    % data from LYES
    case 1
        meas_data = TestData.Cell_1.T2.WLTP;
                                                    
        % Read out time, voltage and current vectors
        RawTime = meas_data.t;
        RawCurrent = -1.*meas_data.I_cell;
        RawVoltage = meas_data.U_cell;
        
    case 2
        meas_data = data;
        
        % Read out time, voltage and current vectors
        RawTime = meas_data.T2.Time;
        RawCurrent = -1.*meas_data.T2.Current;
        RawVoltage = meas_data.T2.Voltage;
end

% Read out RC parameters of cell model
Temp_vec = rc_param.Temp(1,:);
SoC_vec = flipud(rc_param.SoC);
Em_vec = flipud(rc_param.Em);
R0_vec = flipud(rc_param.R0);
R1_vec = flipud(rc_param.R1);
R2_vec = flipud(rc_param.R2);
C1_vec = flipud(rc_param.C1);
C2_vec = flipud(rc_param.C2);
dOCVdz_vec = flipud(rc_param.dOCVdz);
Qnom_vec = (rc_param.Qnom).*3600;

% Read out time, voltage and current vectors
% rawTime = meas_data.Time;
% rawCurrent = -1.*meas_data.Current;
% rawVoltage = meas_data.Voltage;

% rawTime = meas_data.t;
% rawCurrent = -1.*meas_data.I_cell;
% rawVoltage = meas_data.U_cell;

% up sample to 1s
Time = linspace(0, RawTime(end-1), RawTime(end-1)+1)';
Current = interp1(RawTime, RawCurrent, Time);
Voltage = interp1(RawTime, RawVoltage, Time);

idx = find(Voltage>3.25);
Time = Time(idx);
Current = Current(idx);
Voltage = Voltage(idx);
t_end = Time(end);

% create Temp & SoC grid
[TempGrid,SoCGrid] = meshgrid(Temp_vec, SoC_vec);

% set initial conditions
dt = 1;
% Temp = data.T2.T_start;
Temp = meas_data.T_chamber(1);
% start_soc = interp2(TempGrid, SoCGrid, Em_table, Temp, Voltage(1), 'pchip', 'extrap');
start_soc = interp1(Em_vec(:,2), SoC_vec, Voltage(1), 'pchip', 'extrap');
x_hat = [start_soc; 0; 0];    % zk, i_R1, i_R2

% vals=[0.9, 0.8, 0.75, 0.6,  0.5, 0.4, 0.25, 0.1,...
%     0.09, 0.08, 0.075, 0.06,  0.05, 0.04, 0.025, 0.01,...
%     0.009, 0.008, 0.0075, 0.006,  0.005, 0.004, 0.0025, 0.001,...
%     0.0009, 0.0008, 0.00075, 0.0006,  0.0005, 0.0004, 0.00025, 0.0001,...
%     0.00009, 0.00008, 0.000075, 0.00006,  0.00005, 0.00004, 0.000025, 0.00001,...
%     0.000009, 0.000008, 0.0000075, 0.000006,  0.000005, 0.000004, 0.0000025, 0.000001,...
%     0.0000009, 0.0000008, 0.00000075, 0.0000006,  0.0000005, 0.0000004, 0.00000025, 0.0000001];

vals = 1e-2:0.01:0.5;

m = 10;
RMSE = zeros(m,1);
MAE = zeros(m,1);
MaxAE = zeros(m,1);
tuning_param = zeros(m,5);

% initialize variables
SoC_CC = zeros(length(Time),1);
SoC_EKF = zeros(length(Time),1);
V_EKF = zeros(length(Time),1);
error_vec = zeros(length(Time),1);
zkbnd = zeros(length(Time),1);


for k=1:m
    
    SigmaX = diag([vals(randperm(50,1)),vals(randperm(50,1)),vals(randperm(50,1))]);
    SigmaW = vals(randperm(50,1));
    SigmaV = vals(randperm(50,1));
    
    tuning_param(k,1) = SigmaW;
    tuning_param(k,2) = SigmaV;
    tuning_param(k,3) = SigmaX(1,1);
    tuning_param(k,4) = SigmaX(2,2);
    tuning_param(k,5) = SigmaX(3,3);
    
    for i=1:length(Time)
        
        % read out Em, R0, R1, C1, R2, C2 specific to SoC & Temp
        Em = interp2(TempGrid, SoCGrid, Em_vec, Temp, x_hat(1));
        R0 = interp2(TempGrid, SoCGrid, R0_vec, Temp, x_hat(1));
        R1 = interp2(TempGrid, SoCGrid, R1_vec, Temp, x_hat(1));
        R2 = interp2(TempGrid, SoCGrid, R2_vec, Temp, x_hat(1));
        C1 = interp2(TempGrid, SoCGrid, C1_vec, Temp, x_hat(1));
        C2 = interp2(TempGrid, SoCGrid, C2_vec, Temp, x_hat(1));
        Qnom = interp1(Temp_vec, Qnom_vec, Temp, 'linear');
        
        % setup A, B, C, D matrices
        A_hat = [1 0 0; 0 exp(-dt/(R1*C1)) 0; 0 0 exp(-dt/(R2*C2))];
        B_hat = [-dt/Qnom; 1-exp(-dt/(R1*C1)); 1-exp(-dt/(R2*C2))];
        B = [B_hat 0*B_hat];
        dOCVdsoc = interp2(TempGrid, SoCGrid, dOCVdz_vec, Temp, x_hat(1));
        C_hat = [dOCVdsoc -R1 -R2];
        D_hat = 1;
        
        u = Current(i);
        ytrue = Voltage(i);
        
        %%% Kalman filter
        % step 1a:
        x_hat = A_hat*x_hat + B_hat*u;
        
        % step 1b:
        SigmaX = A_hat*SigmaX*A_hat' + B_hat*SigmaW*B_hat';
        
        % step 1c:
        OCV = interp2(TempGrid, SoCGrid, Em_vec, Temp, x_hat(1));
        y_hat = OCV - R0*u - R1*x_hat(2) - R2*x_hat(3);
        
        % step 2a:
        SigmaY = C_hat*SigmaX*C_hat' + D_hat*SigmaV*D_hat';
        L_k = SigmaX*C_hat'/SigmaY;
        
        % step 2b:
        error = ytrue - y_hat;
        %     if error^2 > 100*SigmaY
        %         L(:)=0.0;
        %     end
        x_hat = x_hat + L_k*error;
        %     x_hat(1) = min(1.05, max(-0.05, x_hat(1)));
        
        % step 2c:
        SigmaX = SigmaX - L_k*SigmaY*L_k';
        
        %     % P-bump code
        %
        %     if error^2>4*sigma_y
        %         disp('Bumping P\n');
        %         P(1,1) = P(1,1)*1;
        %     end
        %     if ~isnan(P)
        %     [~,S,V] = svd(P);
        %     HH = V*S*V';
        %     P = (P+P'+HH+HH')/4;
        %     end
        
        % Coloumb Counting
        if i==1
            SoC_temp = start_soc - dt*u/Qnom;
        else
            SoC_temp = SoC_temp - dt*u/Qnom;
        end
        
        % store variables
        SoC_CC(i) = SoC_temp;
        SoC_EKF(i) = x_hat(1);
        V_EKF(i) = y_hat;
        error_vec(i) = error;
        zkbnd(i) = 3*sqrt(SigmaX(1,1));
    end
    
    RMSE(k) = sqrt(mean((SoC_CC-SoC_EKF).^2))*100;
    MAE(k) = mean(abs((SoC_CC-SoC_EKF)))*100;
    MaxAE(k) = max(abs((SoC_CC-SoC_EKF)))*100;
    k
end

save('RMSE_2RC.mat', 'RMSE')
save('tuning_param_2RC.mat', 'tuning_param')

toc

% SoC_CC = SoC_CC';
% SoC_EKF = SoC_EKF';
% V_EKF = V_EKF';
% RMSE = sqrt(mean(error_vec.^2));

% close all
% plot(Time, SoC_CC)
% hold on
% plot(Time, SoC_EKF)

% plot(TimeInterp, VoltageInterp)
% hold on
% plot(TimeInterp, V_EKF)
% ylim([2.5 4.2])
% plot(TimeInterp, error_vec)












