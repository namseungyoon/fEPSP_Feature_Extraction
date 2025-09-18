function[state]=EDA_EKF_Denoising_for_test(z, Sudomotor_Response, fs)

n=5;
H=[1 0 1 1 0];

Q=zeros(5,5);
Q(4,4)=1e-3;
Q(4,5)=0.01;
Q(5,4)=0.01;
Q(5,5)=0.01;

kR=100; %WGN 빼면 줄일필요있음
R0=5.0; %WGN 빼면 줄일필요있음

x=zeros(n,1);
P=1*eye(n);

dmin=0.01;
dmax=0.5;
U(Sudomotor_Response<dmin | Sudomotor_Response>dmax)=0;
U(Sudomotor_Response>=dmin & Sudomotor_Response<=dmax)=1;

state=[];
cond1_array = [];
for i=1:length(z)
    cond1=abs(z(i)-z(max([i-1,1])));
    cond1_array = [cond1_array;cond1];
end

subplot(4,1,1);
plot(z)
hold on;
subplot(4,1,2);
plot(Sudomotor_Response);
yline(-0.15);
yline(0.15);
subplot(4,1,3);
plot(cond1_array);
yline(0.0001);
yline(1);

for i=1:length(z)

    A=Ajacob(x,U(i),fs);

    xp=A*x;
    Pp=A*P*A'+Q;

    noise_est=std(z(max([i-fs*10+1,1]):i));
    R=kR*noise_est*+R0;

    cond1=abs(z(i)-z(max([i-1,1])));    % fs와 비례할 수 있다
    cond2=Sudomotor_Response(i);        % 미분의 개념(fs과 관련있음)
    if cond1 > 0.001 && cond1 < 1 && cond2 > -0.2 && cond2 < 0.2
        %노이즈가 있을때
%         subplot(3,1,1);
%         scatter(i,z(i));
%         hold on;
%         cond1
%         cond2

        K=Pp*H'*inv(H*Pp*H'+R);
        
    else
        %노이즈가 없을때
        K=zeros(n,1);
    end

    zz=H*xp;
    x=xp+K*(z(i)-zz);
    P=Pp-K*H*Pp;
    state=[state,x];
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






