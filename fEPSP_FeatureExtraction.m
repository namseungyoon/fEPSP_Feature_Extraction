clear all
close all

load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E17.mat");
load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E45.mat");
load("D:\Project_2025_2026_HIPPO\Workspace\HippocampalSignalProcessing\DATASET\ETRI\fEPSP_E55.mat");

t = fEPSP_E45(:,1)*10;%ms
dt=mean(diff(t));
fs = 1/dt*10000;
fEPSPs = [fEPSP_E17(:,2:end), fEPSP_E45(:,2:end), fEPSP_E55(:,2:end)];
[~, n_fEPSP] = size(fEPSPs);
figure('Position', [-2000 100 1800 1200])

for i = 1:n_fEPSP
    subplot(3,10,i);
    plot(t, fEPSPs(:,i));
end


for i = 1:n_fEPSP

    fEPSP = fEPSPs(:,i);
    diff_fEPSP = diff(fEPSP);
    diff_fEPSP_t = t(2:end)-dt/2;
    [max_diff_fEPSP, max_idx] = max(diff_fEPSP);
    [min_diff_fEPSP, min_idx] = min(diff_fEPSP);
    
    start_idx = min_idx+1;
    new_fEPSP_t = diff_fEPSP_t(start_idx-1:end);
    new_fEPSP = fEPSP(start_idx:end);
    new_diff_fEPSP = diff(new_fEPSP);
    new_diff_fEPSP_t = new_fEPSP_t(2:end)-dt/2;

    %LPF
    fc = 900;    % 차단 주파수 (Hz), fEPSP는 보통 300Hz 이하 관찰
    [b,a] = butter(4, fc/(fs/2), 'low'); 
    new_ffEPSP = filtfilt(b,a,new_fEPSP);
    new_diff_ffEPSP = diff(new_ffEPSP);

    plot(new_ffEPSP)

    % peak detection
    % ---- 1. 양의 피크 찾기 ----
    [peakVals, peakLocs] = findpeaks(new_ffEPSP, ...
    'MinPeakHeight', 50, ...      % 최소 50 μV 이상
    'MinPeakProminence', 30, ...  % 주변보다 최소 30 μV 돌출
    'MinPeakDistance', 20);        % 5 ms 이상 간격


    % ---- 2. 음의 피크 찾기 ----
    [troughVals, troughLocs] = findpeaks(-new_ffEPSP, ...
    'MinPeakHeight', 50, ...      % 최소 -50 μV 이하
    'MinPeakProminence', 30, ...
    'MinPeakDistance', 20);
    troughVals = -troughVals;   % 부호를 원래대로 되돌림
    
    if ( (length(troughVals) & length(peakVals)) ~= 0 )
        first_peak_t = new_fEPSP_t(peakLocs(1));
        first_peak = peakVals(1);
        first_trough_t = new_fEPSP_t(troughLocs(1));
        first_trough = troughVals(1);
    
        roi_new_ffEPSP = new_ffEPSP(troughLocs(1):peakLocs(1));
        roi_new_ffEPSP_t = new_fEPSP_t(troughLocs(1):peakLocs(1));
    
        roi_new_ffEPSP = roi_new_ffEPSP((first_trough < roi_new_ffEPSP) & (roi_new_ffEPSP < first_peak));
        roi_new_ffEPSP_t = roi_new_ffEPSP_t((first_trough < roi_new_ffEPSP) & (roi_new_ffEPSP < first_peak));
        min_roi_new_ffEPSP  = min(roi_new_ffEPSP);
        start_slope = roi_new_ffEPSP_t(1);
        max_roi_new_ffEPSP = max(roi_new_ffEPSP);
        end_slope = roi_new_ffEPSP_t(end);
    end

    subplot(2,2,1);
    plot(t, fEPSP);
    hold on;
    scatter(max_idx,fEPSP(max_idx));
    hold on;
    scatter(min_idx,fEPSP(min_idx));
    hold on;
    scatter(min_idx+1,fEPSP(min_idx+1),"*");
    hold off;

    subplot(2,2,2);
    plot(new_fEPSP_t, new_fEPSP);
    hold on;
    plot(new_fEPSP_t, new_ffEPSP, "Color",'r');
    hold on;
    if ( (length(troughVals) & length(peakVals)) ~= 0 )
        plot(new_fEPSP_t(peakLocs), peakVals, 'ro','MarkerFaceColor','r','DisplayName','Peak');
        hold on;
        plot(new_fEPSP_t(troughLocs), troughVals, 'bo','MarkerFaceColor','b','DisplayName','Trough');
        hold on;
        plot([start_slope;end_slope],[min_roi_new_ffEPSP; max_roi_new_ffEPSP]);
    end
    hold off;
    xlabel('Time (ms)');              % X축 레이블
    ylabel('Amplitude (\muV)');       % Y축 레이블
    title('fEPSP: Raw, LPF');

    subplot(2,2,3);
    plot(diff_fEPSP_t, diff_fEPSP);
    hold on;
    scatter(max_idx,diff_fEPSP(max_idx));
    hold on;
    scatter(min_idx,diff_fEPSP(min_idx));
    hold on;
    scatter(min_idx+1,diff_fEPSP(min_idx+1),"*");
    hold off;


    subplot(2,2,4);
    plot(new_diff_fEPSP_t, new_diff_fEPSP);
    hold on;
    plot(new_diff_fEPSP_t, new_diff_ffEPSP);
    hold off;
    sgtitle(sprintf('Signal #%d', i));   % 각 subplot의 소제목

end
