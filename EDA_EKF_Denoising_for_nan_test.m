function[state, cond1_array]=EDA_EKF_Denoising_for_nan_test(z, Sudomotor_Response, fs, cond1_low, cond1_high, cond2_low, cond2_high)
% save("noise_signal.mat","z")
n=5;
H=[1 0 1 1 0];

Q=zeros(5,5);
Q(4,4)=1e-3;
Q(4,5)=0.01;
Q(5,4)=0.01;
Q(5,5)=0.01;

kR=100;
R0=5.0;
% kR=10;
% R0=0.5;

x=zeros(n,1);
P=1*eye(n);

dmin=0.01;
dmax=0.5;
U(Sudomotor_Response<dmin | Sudomotor_Response>dmax)=0;
U(Sudomotor_Response>=dmin & Sudomotor_Response<=dmax)=1;

state=[];
cond1_array = [];
for i=1:length(z)
    if i > 1550
        i

        A=Ajacob(x,U(i),fs)
        
        xp=A*x
        Pp=A*P*A'+Q
        
        noise_est=std(z(max([i-fs*10+1,1]):i))
        R=kR*noise_est*+R0
    else
        A=Ajacob(x,U(i),fs);
        
        xp=A*x;
        Pp=A*P*A'+Q;
        
        noise_est=std(z(max([i-fs*10+1,1]):i));
        R=kR*noise_est*+R0;
    end
        

    cond1=abs(z(i)-z(max([i-1,1])));
    cond1_array = [cond1_array;cond1];
    cond2=Sudomotor_Response(i);
    if cond1 > cond1_low && cond1 < cond1_high && cond2 > cond2_low && cond2 < cond2_high %%
%         K=Pp*H'*inv(H*Pp*H'+R);
        K=Pp*H'/(H*Pp*H'+R);

    else
        K=zeros(n,1);
    end

    zz=H*xp;
    x=xp+K*(z(i)-zz);
    P=Pp-K*H*Pp;
    state=[state,x];
%     if isnan(x)
%         subplot(3,1,1);
%         plot(z(1:i));
%         hold on;
%         plot(H*state);
%         ylim([min(z) max(z)])
%         subplot(3,1,2);
%         plot(cond1_array);
%         yline(cond1_low);
%         yline(1);
%         subplot(3,1,3);
%         plot(Sudomotor_Response(1:i));
%         yline(cond2_high);
%         yline(cond2_low);
%         fprintf("kalman cond1 : %d, cond2 : %d\n",cond1_low,cond2_high);
%     end

end
kalman_EDA=(H*state)';
nn = sum(isnan(kalman_EDA));
if nn > 0
    figure(nn);
    subplot(3,1,1);
    plot(z)
    hold on;
    plot(kalman_EDA);
    ylim([min(z), max(z)]);
    subplot(3,1,2);
    plot(cond1_array);
    yline(cond1_low);
    yline(cond1_high);
    hold on;
    legend('cond1 array', num2str(cond1_high), num2str(cond1_low));
    subplot(3,1,3);
    plot(Sudomotor_Response);
    yline(cond2_low);
    yline(cond2_high);
    hold off;
    legend('Sudomotor Response', num2str(cond2_high), num2str(cond2_low));
    fprintf("pause\n")
end
%--------------------------------------------------------------------------
function A=Ajacob(x,U,fs)
%   x=[SCH, kdiff, SC0, SCR, S]';

SCH=x(1);
kdiff=x(2);
SC0=x(3);
SCR=x(4);
S=x(5);

dT=1/fs;
alpha=31;
beta=0.016;
krec=0.071;
A=[1-dT*kdiff -dT*SCH 0 dT*beta 0;...
    0 1 0 0 0;...
    0 0 1 0 0;...
    0 0 0 1-dT*krec dT*alpha;...
    0 0 0 0 U-2*S];






