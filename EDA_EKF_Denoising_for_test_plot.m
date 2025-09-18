function[state, cond1_array]=EDA_EKF_Denoising_for_test_plot(z, Sudomotor_Response, fs, cond1_low, cond1_high, cond2_low, cond2_high, file_name, itr)
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
Rs =[];
noise_ests = [];
kdiffs = [];
SCHs = [];
Ss = [];

idx = 1:length(z);
figure(itr);
set(gcf, 'Position', [2561 40 2561 1320]);
subplot(3,1,1);
plot(idx,z);
hold on;

subplot(3,1,3);
plot(idx,Sudomotor_Response);
hold on;
yline(cond2_high, 'LineWidth',1,'Color','r','LineStyle','--');
yline(cond2_low, 'LineWidth',1,'Color','r','LineStyle','--');
hold off;
str_alpha = strcat('cond2 high = ',num2str(cond2_high));
legend('cond2', str_alpha);


for i=1:length(z)

    A=Ajacob(x,U(i),fs);
    kdiffs = [kdiffs,(-A(1,1)+1)/(1/fs)];
    SCHs = [SCHs,A(1,2)/-(1/fs)];
    Ss = [Ss,(-A(5,5)+U(i))/2];
    
    xp=A*x;
    Pp=A*P*A'+Q;

    noise_est=std(z(max([i-fs*10+1,1]):i));%문제없음
    noise_ests = [noise_ests;noise_est];
    R=kR*noise_est*+R0;%문제없음
    Rs =[Rs;R];

    cond1=abs(z(i)-z(max([i-1,1])));
    cond1_array = [cond1_array;cond1];

    cond2=Sudomotor_Response(i);
    if cond1 > cond1_low && cond1 < cond1_high && cond2 > cond2_low && cond2 < cond2_high
        K=Pp*H'/(H*Pp*H'+R);% K=Pp*H'*inv(H*Pp*H'+R);
    else
        K=zeros(n,1);
    end

    zz=H*xp;
    x=xp+K*(z(i)-zz);
    %constraint 범위를 벗어나지 않게 

    if (x(5) < 0) || isnan(x(5))
        x(5) = 0;
    end

    P=Pp-K*H*Pp;
    state=[state,x];

end
subplot(3,1,2);
plot(idx,cond1_array);
hold on;
yline(cond1_high, 'LineWidth',1,'Color','r','LineStyle','--');
yline(cond1_low, 'LineWidth',1,'Color','g','LineStyle','--');
hold off;
ylim([-0.01,cond1_high+0.1]);
str_high = strcat('cond1 high = ',num2str(cond1_high));
str_low = strcat('cond1 low = ',num2str(cond1_low));
legend('cond1', str_high, str_low);

EDA_kalman=([1 0 1 1 0]*state)';
subplot(3,1,1);
plot(idx, EDA_kalman);
hold off

function A=Ajacob(x,U,fs)
%   x=[SCH, kdiff, SC0, SCR, S]';

SCH=x(1);
kdiff=x(2);
SC0=x(3);
SCR=x(4);
S=x(5);
% if S < 0
%     S = 0;
% end
% S=0;


dT=1/fs;
alpha=31;
beta=0.016;
krec=0.071;
A=[1-dT*kdiff -dT*SCH 0 dT*beta 0;...
    0 1 0 0 0;...
    0 0 1 0 0;...
    0 0 0 1-dT*krec dT*alpha;...
    0 0 0 0 U-2*S];






