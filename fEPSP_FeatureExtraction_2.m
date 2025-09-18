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
n=0;

fs = 1/(t_E17(2)-t_E17(1))*1000;

lp_cutoff = 300;  % 컷오프 주파수
filter_order = 3;  % 필터 차수
lp_fir = fir1(filter_order, lp_cutoff/(fs/2));  % FIR 필터 설계 (Hamming window)

for n_electrode = 1:length(n_electrodes)
    fEPSPs = fEPSP_all(:,:,n_electrode);
    [len_fEPSPs, n_fEPSPs] = size(fEPSPs);
    for n_fEPSP = 1:n_fEPSPs
        fEPSP_t = ts(:,n_electrode);
        fEPSP = fEPSPs(:,n_fEPSP);
        [fEPSP_stim_peak, fEPSP_stim_peak_i] = max(fEPSP);

        fEPSP_1order_t = fEPSP_t(2:end)-0.05;
        fEPSP_1order = diff(fEPSP);
        [min_fEPSP_1order, min_fEPSP_1order_i] = min(fEPSP_1order);

        fEPSP_2order_t = fEPSP_t(2:end-1);
        fEPSP_2order = diff(fEPSP_1order);

        [min_fEPSP_2order, min_fEPSP_2order_i] = min(fEPSP_2order);

        % roi_logic = fEPSP_t > fEPSP_1order_t(min_fEPSP_1order_i);
        roi_shift_i = 2;
        roi_t = fEPSP_t(min_fEPSP_1order_i+roi_shift_i:end);
        roi_fEPSP = fEPSP(min_fEPSP_1order_i+roi_shift_i:end);
        [min_roi_fEPSP, min_roi_fEPSP_i] = min(roi_fEPSP);
        [max_roi_fEPSP, max_roi_fEPSP_i] = max(roi_fEPSP);
        
        filtered_roi_fEPSP = filtfilt(lp_fir, 1, roi_fEPSP);
        [peak_values, peak_locs] = findpeaks(filtered_roi_fEPSP, fs, 'MinPeakProminence', 0.5, 'MinPeakDistance', 0.003);  % 필요시 Threshold 조절
        peak_locs = peak_locs * 1000 + roi_t(1);



        [trough_values, trough_locs] = findpeaks(-filtered_roi_fEPSP, fs, 'MinPeakProminence', 0.5);
        trough_values = -trough_values;  % 다시 부호 복원
        trough_locs = trough_locs * 1000 + roi_t(1);

        if (min(filtered_roi_fEPSP) < min(trough_values) )
            [fEPSP_start, fEPSP_start_i] = min(filtered_roi_fEPSP);
            fEPSP_start_t = roi_t(fEPSP_start_i);
        else
            fEPSP_start = trough_values(1);
            fEPSP_start_t = trough_locs(1);
        end


        fEPSP_amplitude = peak_values(1);
        fEPSP_amplitude_t = peak_locs(1);

        fEPSP_slope_start = fEPSP_start + (fEPSP_amplitude - fEPSP_start) * 0.2;
        fEPSP_slope_end = fEPSP_start + (fEPSP_amplitude - fEPSP_start) * 0.8;

        fEPSP_slope = filtered_roi_fEPSP(fEPSP_start_t <= roi_t & roi_t <= fEPSP_amplitude_t);
        fEPSP_slope_t = roi_t(fEPSP_start_t <= roi_t & roi_t <= fEPSP_amplitude_t);
        
        fEPSP_slope_min = min(fEPSP_slope);
        fEPSP_slope_max = max(fEPSP_slope);

        fEPSP_slope_range = fEPSP_slope_max - fEPSP_slope_min;
        fEPSP_slope_10 = fEPSP_slope_min + 0.1 * fEPSP_slope_range;
        fEPSP_slope_90 = fEPSP_slope_min + 0.9 * fEPSP_slope_range;
        in_range_idx = find(fEPSP_slope >= fEPSP_slope_10 & fEPSP_slope <= fEPSP_slope_90); 
        fEPSP_slope_t_10 = fEPSP_slope_t(in_range_idx(1));
        fEPSP_slope_t_90 = fEPSP_slope_t(in_range_idx(end));
        fEPSP_slope_10 = fEPSP_slope(in_range_idx(1));
        fEPSP_slope_90 = fEPSP_slope(in_range_idx(end));
        fEPSP_slope_value = (fEPSP_slope_90 - fEPSP_slope_10) / (fEPSP_slope_t_90 - fEPSP_slope_t_10);
        n = n + 1;
        
        % figure(1);
        % subplot(4,1,1)
        % plot(fEPSP_t, fEPSP);
        % hold on;
        % scatter(fEPSP_t(fEPSP_stim_peak_i), fEPSP_stim_peak, 'red');
        % hold on;
        % scatter(fEPSP_t(min_fEPSP_1order_i+roi_shift_i), fEPSP(min_fEPSP_1order_i+roi_shift_i))
        % hold on;
        % scatter(roi_t(max_roi_fEPSP_i), roi_fEPSP(max_roi_fEPSP_i), 'magenta', '*')
        % hold off;
        % title(["stim peak time = ", num2str(fEPSP_t(fEPSP_stim_peak_i))])
        % 
        % subplot(4,1,2)
        % plot(fEPSP_1order_t, fEPSP_1order);
        % hold off;
        % title(["first order min time = ", num2str(fEPSP_1order_t(min_fEPSP_1order_i))])
        % 
        % subplot(4,1,3)
        % plot(roi_t, roi_fEPSP);
        % hold on;
        % scatter(roi_t(max_roi_fEPSP_i), roi_fEPSP(max_roi_fEPSP_i), 'magenta', '*')
        % hold off
        % 
        % subplot(4,1,4)
        % plot(roi_t, filtered_roi_fEPSP);
        % hold on;
        % scatter(peak_locs, peak_values, '*');
        % hold on;
        % scatter(trough_locs, trough_values, 'x');
        % hold on;
        % scatter(fEPSP_start_t, fEPSP_start, 'o');
        % hold on;
        % scatter(fEPSP_slope_t, fEPSP_slope, 'blue', 'filled', 'o');
        % hold on;
        % plot([fEPSP_slope_t_10, fEPSP_slope_t_90], [fEPSP_slope_10, fEPSP_slope_90], "LineWidth",2, "Color", [1, 0.5, 0.15]);
        % hold off;

        figure(2);
        subplot(3,10,n);
        plot(roi_t, filtered_roi_fEPSP);
        hold on;
        scatter(peak_locs, peak_values, '*');
        hold on;
        scatter(trough_locs, trough_values, 'x');
        hold on;
        scatter(fEPSP_start_t, fEPSP_start, 'o');
        hold on;
        scatter(fEPSP_slope_t, fEPSP_slope, 'blue', 'filled', 'o');
        hold on;
        plot([fEPSP_slope_t_10, fEPSP_slope_t_90], [fEPSP_slope_10, fEPSP_slope_90], "LineWidth",2, "Color", [1, 0.5, 0.15]);
        hold off;
        title(strcat("Amp = ", num2str(fEPSP_amplitude), "Slope = ", num2str(fEPSP_slope_value)));



        



        % figure(1);
        % subplot(length(n_electrodes), n_fEPSPs, n);
        % plot(t, fEPSP,"Marker",".");
        % title(n_electrodes(n_electrode)+" ("+num2str(n_fEPSP)+"V)")
        
        % figure(2);
        % subplot(length(n_electrodes), n_fEPSPs, n);
        % % plot(roi_t(1:50), roi_fEPSP(1:50),"Marker",".");
        % plot(roi_t, roi_fEPSP,"Marker",".");
        % % plot(t, fEPSP,"Marker",".");
        % title(n_electrodes(n_electrode)+" ("+num2str(n_fEPSP)+"V)")

    end
end

