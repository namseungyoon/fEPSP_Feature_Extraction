clear all
close all

% load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E17.mat");
% load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E45.mat");
% load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E55.mat");
load("fEPSP_E17.mat")
load("fEPSP_E45.mat")
load("fEPSP_E55.mat")
t_E17 = fEPSP_E17(:,1);
t_E45 = fEPSP_E45(:,1);
t_E55 = fEPSP_E55(:,1);

fEPSP_E17 = fEPSP_E17(:,2:end);
fEPSP_E45 = fEPSP_E45(:,2:end);
fEPSP_E55 = fEPSP_E55(:,2:end);

n_electrodes = ["E17", "E45", "E55"];
ts = [t_E17, t_E45, t_E55];     % ms
fEPSP_all = cat(length(n_electrodes), fEPSP_E17, fEPSP_E45, fEPSP_E55);
max_fEPSP_all = max(max(max(fEPSP_all)))*1.1/1000;
min_fEPSP_all = min(min(min(fEPSP_all)))*1.1/1000;
n=0;

fs = 1/(t_E17(2)-t_E17(1))*1000;    % hz

lp_cutoff = 300;  % 컷오프 주파수
filter_order = 3;  % 필터 차수
lp_fir = fir1(filter_order, lp_cutoff/(fs/2));  % FIR 필터 설계 (Hamming window)

first_fEPSPs = fEPSP_all(:,:,1);
second_fEPSPs = fEPSP_all(:,:,2);
third_fEPSPs = fEPSP_all(:,:,3);
fEPSP_t = t_E17;
[len_fEPSPs, n_fEPSPs] = size(first_fEPSPs);

first_fEPSP_Amps = [];
second_fEPSP_Amps = [];
third_fEPSP_Amps = [];

for n_fEPSP = 1:n_fEPSPs
    first_fEPSP = first_fEPSPs(:,n_fEPSP)/1000;
    second_fEPSP = second_fEPSPs(:,n_fEPSP)/1000;
    third_fEPSP = third_fEPSPs(:,n_fEPSP)/1000;

    %% detection stimulation peak
    [first_fEPSP_stim_peak, first_fEPSP_stim_peak_i] = max(first_fEPSP);
    [second_fEPSP_stim_peak, second_fEPSP_stim_peak_i] = max(second_fEPSP);
    [third_fEPSP_stim_peak, third_fEPSP_stim_peak_i] = max(third_fEPSP);

    %% RoI fEPSP
    % 피크이후 1mS 신호는 artifact noise로 간주
    first_RoI_fEPSP_ts = fEPSP_t(fEPSP_t >= fEPSP_t(first_fEPSP_stim_peak_i) + 1);
    first_RoI_fEPSP = first_fEPSP(fEPSP_t >= fEPSP_t(first_fEPSP_stim_peak_i) + 1);
    second_RoI_fEPSP_ts = fEPSP_t(fEPSP_t >= fEPSP_t(second_fEPSP_stim_peak_i) + 1);
    second_RoI_fEPSP = second_fEPSP(fEPSP_t >= fEPSP_t(second_fEPSP_stim_peak_i) + 1);
    third_RoI_fEPSP_ts = fEPSP_t(fEPSP_t >= fEPSP_t(third_fEPSP_stim_peak_i) + 1);
    third_RoI_fEPSP = third_fEPSP(fEPSP_t >= fEPSP_t(third_fEPSP_stim_peak_i) + 1);

    % LPF
    first_fEPSP_RoI_filtered = filtfilt(lp_fir, 1, first_RoI_fEPSP);
    second_fEPSP_RoI_filtered = filtfilt(lp_fir, 1, second_RoI_fEPSP);
    third_fEPSP_RoI_filtered = filtfilt(lp_fir, 1, third_RoI_fEPSP);

    %% 

    figure(n_fEPSPs);
    subplot(3,2,1);
    plot(fEPSP_t, first_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(first_RoI_fEPSP_ts, first_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(first_fEPSP_stim_peak_i), first_fEPSP_stim_peak, '*');
    hold off;
    subplot(3,2,2);
    plot(first_RoI_fEPSP_ts, first_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold off;

    subplot(3,2,3);
    plot(fEPSP_t, second_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(second_RoI_fEPSP_ts, second_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(second_fEPSP_stim_peak_i), second_fEPSP_stim_peak, '*');
    hold off;
    subplot(3,2,4);
    plot(second_RoI_fEPSP_ts, second_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold off;

    subplot(3,2,5);
    plot(fEPSP_t, third_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(third_RoI_fEPSP_ts, third_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(third_fEPSP_stim_peak_i), third_fEPSP_stim_peak, '*');
    hold off;
    subplot(3,2,6);
    plot(third_RoI_fEPSP_ts, third_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold off;
end

figure();
subplot(3,1,1);
plot(first_fEPSP_Amps);
xlabel("Stimulation intensity(V)");
ylabel("Amplitude(mV)");
title("Electrode 17");
subplot(3,1,2);
plot(second_fEPSP_Amps);
xlabel("Stimulation intensity(V)");
ylabel("Amplitude(mV)");
title("Electrode 45");
subplot(3,1,3);
plot(third_fEPSP_Amps);
xlabel("Stimulation intensity(V)");
ylabel("Amplitude(mV)");
title("Electrode 45");
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
