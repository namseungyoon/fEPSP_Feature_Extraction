function [pts, idx_alpha, s, L] = arcPointsUniform(y, dt, fracs, t0)
% arcPointsUniform
% 등간격 샘플 신호 y와 샘플 간격 dt가 주어졌을 때,
% 아크 길이 기준 fracs(예: [0.2 0.8])에 해당하는 (t*, y*)를 구합니다.
%
% 입력:
%   y      : 열/행 벡터 (신호값)
%   dt     : 등간격 샘플 간격 (기본 1.0)
%   fracs  : 원하는 비율들 (기본 [0.2 0.8])
%   t0     : 시작 시간 (기본 0.0)
%
% 출력:
%   pts       : [numel(fracs) x 2] 행렬, 각 행이 [t*, y*]
%   idx_alpha : [numel(fracs) x 2], 각 행이 [j, alpha]
%   s         : 누적 아크 길이 벡터(길이 N)
%   L         : 전체 아크 길이(스칼라)

    if nargin < 2 || isempty(dt),   dt = 1.0; end
    if nargin < 3 || isempty(fracs), fracs = [0.2, 0.8]; end
    if nargin < 4 || isempty(t0),   t0 = 0.0; end

    y = y(:);                  % 열벡터화
    N = numel(y);
    if N < 2
        error('y는 최소 2개 이상의 샘플이 필요합니다.');
    end

    % 1) 인접 차분과 구간 아크 길이
    dy = diff(y);
    ds = hypot(dt, dy);        % sqrt(dt^2 + dy^2)보다 수치적으로 안정적

    % 2) 누적 길이 및 전체 길이
    s = [0; cumsum(ds)];
    L = s(end);

    % 모든 세그먼트 길이가 0인 경우 (완전 평탄/중복)
    if L == 0
        t_star = t0;
        y_star = y(1);
        pts = repmat([t_star, y_star], numel(fracs), 1);
        idx_alpha = [ones(numel(fracs),1), zeros(numel(fracs),1)];
        return;
    end

    % 3) 각 비율에 대해 역추정(보간)
    pts = zeros(numel(fracs), 2);
    idx_alpha = zeros(numel(fracs), 2);

    for k = 1:numel(fracs)
        target = fracs(k) * L;

        % s_j <= target <= s_{j+1}인 j 찾기
        % find(...,'last')는 monotonic s에서 안전합니다.
        j = find(s <= target, 1, 'last');
        if j == N           % target이 정확히 s(end)면 마지막 세그먼트로 조정
            j = N - 1;
        end

        seg = s(j+1) - s(j);
        if seg == 0
            alpha = 0.0;    % 길이 0 세그먼트(중복점) 보호
        else
            alpha = (target - s(j)) / seg;
        end

        % 시간·값 보간
        t_star = (t0 + (j-1)*dt) + alpha*dt;
        y_star = y(j) + alpha*(y(j+1) - y(j));

        pts(k,:) = [t_star, y_star];
        idx_alpha(k,:) = [j, alpha];
    end
end
