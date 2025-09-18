clear all
close all

load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E17.mat");
load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E45.mat");
load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E55.mat");

t_E17 = fEPSP_E17(:,1);
t_E45 = fEPSP_E45(:,1);
t_E55 = fEPSP_E55(:,1);

fEPSP_E17 = fEPSP_E17(:,2:end);
fEPSP_E45 = fEPSP_E45(:,2:end);
fEPSP_E55 = fEPSP_E55(:,2:end);

n_electrodes = ["E17", "E45", "E55"];
ts = [t_E17, t_E45, t_E55];
fEPSP_all = cat(length(n_electrodes), fEPSP_E17, fEPSP_E45, fEPSP_E55);
max_fEPSP_all = max(max(max(fEPSP_all)))*1.1;
min_fEPSP_all = min(min(min(fEPSP_all)))*1.1;
n=0;

fs = 1/(t_E17(2)-t_E17(1))*1000;

lp_cutoff = 300;  % 컷오프 주파수
filter_order = 3;  % 필터 차수
lp_fir = fir1(filter_order, lp_cutoff/(fs/2));  % FIR 필터 설계 (Hamming window)

first_fEPSPs = fEPSP_all(:,:,1);
second_fEPSPs = fEPSP_all(:,:,2);
third_fEPSPs = fEPSP_all(:,:,3);
fEPSP_t = t_E17;
[len_fEPSPs, n_fEPSPs] = size(first_fEPSPs);

for n_fEPSP = 1:n_fEPSPs
    first_fEPSP = first_fEPSPs(:,n_fEPSP);
    second_fEPSP = second_fEPSPs(:,n_fEPSP);
    third_fEPSP = third_fEPSPs(:,n_fEPSP);

    %% detection stimulation peak
    [first_fEPSP_stim_peak, first_fEPSP_stim_peak_i] = max(first_fEPSP);
    [second_fEPSP_stim_peak, second_fEPSP_stim_peak_i] = max(second_fEPSP);
    [third_fEPSP_stim_peak, third_fEPSP_stim_peak_i] = max(third_fEPSP);

    %% 1order fEPSP
    fEPSP_1order_t = fEPSP_t(2:end)-0.05;
    first_fEPSP_1order = diff(first_fEPSP);
    second_fEPSP_1order = diff(second_fEPSP);
    third_fEPSP_1order = diff(third_fEPSP);
    % 기울기가 가장 가파르게 떨어지는 지점
    [min_first_fEPSP_1order, min_first_fEPSP_1order_i] = min(first_fEPSP_1order);
    [min_second_fEPSP_1order, min_second_fEPSP_1order_i] = min(second_fEPSP_1order);
    [min_third_fEPSP_1order, min_third_fEPSP_1order_i] = min(third_fEPSP_1order);

    % RoI fEPSP
    first_fEPSP_RoI_t = fEPSP_t(fEPSP_1order_t(min_first_fEPSP_1order_i) < fEPSP_t);
    first_fEPSP_RoI = first_fEPSP(fEPSP_1order_t(min_first_fEPSP_1order_i) < fEPSP_t);
    second_fEPSP_RoI_t = fEPSP_t(fEPSP_1order_t(min_second_fEPSP_1order_i) < fEPSP_t);
    second_fEPSP_RoI = second_fEPSP(fEPSP_1order_t(min_second_fEPSP_1order_i) < fEPSP_t);
    third_fEPSP_RoI_t = fEPSP_t(fEPSP_1order_t(min_third_fEPSP_1order_i) < fEPSP_t);
    third_fEPSP_RoI = third_fEPSP(fEPSP_1order_t(min_third_fEPSP_1order_i) < fEPSP_t);

    % LPF
    first_fEPSP_RoI_filtered = filtfilt(lp_fir, 1, first_fEPSP_RoI);
    second_fEPSP_RoI_filtered = filtfilt(lp_fir, 1, second_fEPSP_RoI);
    third_fEPSP_RoI_filtered = filtfilt(lp_fir, 1, third_fEPSP_RoI);

    
    [first_fEPSP_RoI_filtered_trough, first_fEPSP_RoI_filtered_trough_i] = findpeaks(-first_fEPSP_RoI_filtered, fs, 'MinPeakProminence', 0.5);
    first_fEPSP_RoI_filtered_trough = -first_fEPSP_RoI_filtered_trough;
    first_fEPSP_RoI_filtered_trough_i = first_fEPSP_RoI_filtered_trough_i * 1000 + first_fEPSP_RoI_t(1);
    [second_fEPSP_RoI_filtered_trough, second_fEPSP_RoI_filtered_trough_i] = findpeaks(-second_fEPSP_RoI_filtered, fs, 'MinPeakProminence', 0.5);
    second_fEPSP_RoI_filtered_trough = -second_fEPSP_RoI_filtered_trough;
    second_fEPSP_RoI_filtered_trough_i = second_fEPSP_RoI_filtered_trough_i * 1000 + second_fEPSP_RoI_t(1);
    [third_fEPSP_RoI_filtered_trough, third_fEPSP_RoI_filtered_trough_i] = findpeaks(-third_fEPSP_RoI_filtered, fs, 'MinPeakProminence', 0.5);
    third_fEPSP_RoI_filtered_trough = -third_fEPSP_RoI_filtered_trough;
    third_fEPSP_RoI_filtered_trough_i = third_fEPSP_RoI_filtered_trough_i * 1000 + third_fEPSP_RoI_t(1);

    % peak detection 
    [first_fEPSP_RoI_filtered_peak, first_fEPSP_RoI_filtered_peak_i] = findpeaks(first_fEPSP_RoI_filtered, fs, 'MinPeakProminence', 0.5, 'MinPeakDistance', 0.003);  % 필요시 Threshold 조절
    first_fEPSP_RoI_filtered_peak_i = first_fEPSP_RoI_filtered_peak_i * 1000 + first_fEPSP_RoI_t(1);
    [second_fEPSP_RoI_filtered_peak, second_fEPSP_RoI_filtered_peak_i] = findpeaks(second_fEPSP_RoI_filtered, fs, 'MinPeakProminence', 0.5, 'MinPeakDistance', 0.003);  % 필요시 Threshold 조절
    second_fEPSP_RoI_filtered_peak_i = second_fEPSP_RoI_filtered_peak_i * 1000 + second_fEPSP_RoI_t(1);
    [third_fEPSP_RoI_filtered_peak, third_fEPSP_RoI_filtered_peak_i] = findpeaks(third_fEPSP_RoI_filtered, fs, 'MinPeakProminence', 0.5, 'MinPeakDistance', 0.003);  % 필요시 Threshold 조절
    third_fEPSP_RoI_filtered_peak_i = third_fEPSP_RoI_filtered_peak_i * 1000 + third_fEPSP_RoI_t(1);

    % fEPSP Amplitute
    first_fEPSP_RoI_filtered_amp = first_fEPSP_RoI_filtered_peak(1);
    first_fEPSP_RoI_filtered_amp_t = first_fEPSP_RoI_filtered_peak_i(1);
    second_fEPSP_RoI_filtered_amp = second_fEPSP_RoI_filtered_peak(1);
    second_fEPSP_RoI_filtered_amp_t = second_fEPSP_RoI_filtered_peak_i(1);
    third_fEPSP_RoI_filtered_amp = third_fEPSP_RoI_filtered_peak(1);
    third_fEPSP_RoI_filtered_amp_t = third_fEPSP_RoI_filtered_peak_i(1);
    
    % start fEPSP detection 
    first_fEPSP_start_t = first_fEPSP_RoI_t(first_fEPSP_RoI_t < first_fEPSP_RoI_filtered_amp_t);
    first_fEPSP_start = first_fEPSP_RoI_filtered(first_fEPSP_RoI_t < first_fEPSP_RoI_filtered_amp_t);
    [first_fEPSP_start, first_fEPSP_start_i] = min(first_fEPSP_start);
    first_fEPSP_start_t = first_fEPSP_start_t(first_fEPSP_start_i);
    second_fEPSP_start_t = second_fEPSP_RoI_t(second_fEPSP_RoI_t < second_fEPSP_RoI_filtered_amp_t);
    second_fEPSP_start = second_fEPSP_RoI_filtered(second_fEPSP_RoI_t < second_fEPSP_RoI_filtered_amp_t);
    [second_fEPSP_start, second_fEPSP_start_i] = min(second_fEPSP_start);
    second_fEPSP_start_t = second_fEPSP_start_t(second_fEPSP_start_i);
    third_fEPSP_start_t = third_fEPSP_RoI_t(third_fEPSP_RoI_t < third_fEPSP_RoI_filtered_amp_t);
    third_fEPSP_start = third_fEPSP_RoI_filtered(third_fEPSP_RoI_t < third_fEPSP_RoI_filtered_amp_t);
    [third_fEPSP_start, third_fEPSP_start_i] = min(third_fEPSP_start);
    third_fEPSP_start_t = third_fEPSP_start_t(third_fEPSP_start_i);

     % plot fEPSP
    figure(n_fEPSP);
    subplot(3,1,1);
    plot(fEPSP_t, first_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(first_fEPSP_RoI_t, first_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(first_fEPSP_stim_peak_i), first_fEPSP_stim_peak, '*');
    hold on;
    scatter(first_fEPSP_RoI_filtered_trough_i, first_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    hold on;
    scatter(first_fEPSP_RoI_filtered_peak_i, first_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    hold on;
    scatter(first_fEPSP_RoI_filtered_amp_t, first_fEPSP_RoI_filtered_amp, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [12/255 0 200/255]);
    hold on;
    scatter(first_fEPSP_start_t, first_fEPSP_start, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [1 100/255 0]);
    hold off;
    title(strcat("E17 - ", num2str(n_fEPSP), "V"));
    legend(["raw fEPSP", "LPF fEPSP(RoI)", "Stimulation", "Trough(RoI)", "Peak(RoI)", strcat("Amplitute = ", num2str(first_fEPSP_RoI_filtered_amp)), "Start Point"]);
    ylim([min_fEPSP_all, max_fEPSP_all]);

    subplot(3,1,2);
    plot(fEPSP_t, second_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(second_fEPSP_RoI_t, second_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(second_fEPSP_stim_peak_i), second_fEPSP_stim_peak, '*');
    hold on;
    scatter(second_fEPSP_RoI_filtered_trough_i, second_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    hold on;
    scatter(second_fEPSP_RoI_filtered_peak_i, second_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    hold on;
    scatter(second_fEPSP_RoI_filtered_amp_t, second_fEPSP_RoI_filtered_amp, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [12/255 0 200/255]);
    hold on;
    scatter(second_fEPSP_start_t, second_fEPSP_start, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [1 100/255 0]);
    hold off;
    title(strcat("E45 - ", num2str(n_fEPSP), "V"));
    legend(["raw fEPSP", "LPF fEPSP(RoI)", "Stimulation", "Trough(RoI)", "Peak(RoI)", strcat("Amplitute = ", num2str(first_fEPSP_RoI_filtered_amp)), "Start Point"]);
    ylim([min_fEPSP_all, max_fEPSP_all]);

    subplot(3,1,3);
    plot(fEPSP_t, third_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(third_fEPSP_RoI_t, third_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(third_fEPSP_stim_peak_i), third_fEPSP_stim_peak, '*');
    hold on;
    scatter(third_fEPSP_RoI_filtered_trough_i, third_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    hold on;
    scatter(third_fEPSP_RoI_filtered_peak_i, third_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    hold on;
    scatter(third_fEPSP_RoI_filtered_amp_t, third_fEPSP_RoI_filtered_amp, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [12/255 0 200/255]);
    hold on;
    scatter(third_fEPSP_start_t, third_fEPSP_start, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [1 100/255 0]);
    hold off;
    title(strcat("E55 - ", num2str(n_fEPSP), "V"));
    legend(["raw fEPSP", "LPF fEPSP(RoI)", "Stimulation", "Trough(RoI)", "Peak(RoI)", strcat("Amplitute = ", num2str(first_fEPSP_RoI_filtered_amp)), "Start Point"]);
    ylim([min_fEPSP_all, max_fEPSP_all]);
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);


    figure(11);
    subplot(3,10,n_fEPSP);
    plot(fEPSP_t, first_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(first_fEPSP_RoI_t, first_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(first_fEPSP_stim_peak_i), first_fEPSP_stim_peak, '*');
    hold on;
    scatter(first_fEPSP_RoI_filtered_trough_i, first_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    hold on;
    scatter(first_fEPSP_RoI_filtered_peak_i, first_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    hold on;
    scatter(first_fEPSP_RoI_filtered_amp_t, first_fEPSP_RoI_filtered_amp, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [12/255 0 200/255]);
    hold on;
    scatter(first_fEPSP_start_t, first_fEPSP_start, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [1 100/255 0]);
    hold off;
    title(strcat("E17 - ", num2str(n_fEPSP), "V"));
    ylim([min_fEPSP_all, max_fEPSP_all]);

    subplot(3,10,n_fEPSP+10);
    plot(fEPSP_t, second_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(second_fEPSP_RoI_t, second_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(second_fEPSP_stim_peak_i), second_fEPSP_stim_peak, '*');
    hold on;
    scatter(second_fEPSP_RoI_filtered_trough_i, second_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    hold on;
    scatter(second_fEPSP_RoI_filtered_peak_i, second_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    hold on;
    scatter(second_fEPSP_RoI_filtered_amp_t, second_fEPSP_RoI_filtered_amp, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [12/255 0 200/255]);
    hold on;
    scatter(second_fEPSP_start_t, second_fEPSP_start, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [1 100/255 0]);
    hold off;
    title(strcat("E45 - ", num2str(n_fEPSP), "V"));
    ylim([min_fEPSP_all, max_fEPSP_all]);

    subplot(3,10,n_fEPSP+20);
    plot(fEPSP_t, third_fEPSP, 'Color', [0,0,1], 'LineWidth',0.8);
    hold on;
    plot(third_fEPSP_RoI_t, third_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    hold on;
    scatter(fEPSP_t(third_fEPSP_stim_peak_i), third_fEPSP_stim_peak, '*');
    hold on;
    scatter(third_fEPSP_RoI_filtered_trough_i, third_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    hold on;
    scatter(third_fEPSP_RoI_filtered_peak_i, third_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    hold on;
    scatter(third_fEPSP_RoI_filtered_amp_t, third_fEPSP_RoI_filtered_amp, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [12/255 0 200/255]);
    hold on;
    scatter(third_fEPSP_start_t, third_fEPSP_start, 100, 'magenta', 'o', 'LineWidth',1.8, 'MarkerEdgeColor', [1 100/255 0]);
    hold off;
    title(strcat("E55 - ", num2str(n_fEPSP), "V"));
    ylim([min_fEPSP_all, max_fEPSP_all]);


    % % plot RoI fEPSP
    % subplot(3,2,2);
    % plot(first_fEPSP_RoI_t, first_fEPSP_RoI, 'Color', [0,0,1], 'LineWidth',0.8);
    % hold on;
    % plot(first_fEPSP_RoI_t, first_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    % hold on;
    % scatter(first_fEPSP_RoI_filtered_trough_i, first_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    % hold on;
    % scatter(first_fEPSP_RoI_filtered_peak_i, first_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    % hold off;
    % title("First Electrode RoI");
    % subplot(3,2,4);
    % plot(second_fEPSP_RoI_t, second_fEPSP_RoI, 'Color', [0,0,1], 'LineWidth',0.8);
    % hold on;
    % plot(second_fEPSP_RoI_t, second_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    % hold on;
    % scatter(second_fEPSP_RoI_filtered_trough_i, second_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    % hold on;
    % scatter(second_fEPSP_RoI_filtered_peak_i, second_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    % hold off;
    % title("Second Electrode RoI");
    % subplot(3,2,6);
    % plot(third_fEPSP_RoI_t, third_fEPSP_RoI, 'Color', [0,0,1], 'LineWidth',0.8);
    % hold on;
    % plot(third_fEPSP_RoI_t, third_fEPSP_RoI_filtered, '--', 'Color', [1,0,0], 'LineWidth',1.8);
    % hold on;
    % scatter(third_fEPSP_RoI_filtered_trough_i, third_fEPSP_RoI_filtered_trough, 'v', 'LineWidth',1.8, 'MarkerEdgeColor', [0.9, 0.7, 0]);
    % hold on;
    % scatter(third_fEPSP_RoI_filtered_peak_i, third_fEPSP_RoI_filtered_peak, '^', 'LineWidth',1.8, 'MarkerEdgeColor', [0, 0.7, 0.9]);
    % hold off;
    % title("Third Electrode RoI");

    % sgtitle(strcat(num2str(n_fEPSP), "V stimulation"));

    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);




end


