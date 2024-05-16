clc;
%% 参数初始化
r = size(G, 1);
tau = ones(r * r, r * r);   % 初始信息素矩阵
tau = 8 .* tau;
K = 200;            % 迭代次数
M = 80;             % 蚂蚁数
START = 1;          % 最短路径的起始点
END = r * r;        % 目标点
alpha = 1;          % 信息素重要程度参数
beta = 8;           % 启发式因子重要程度的参数
rho = 0.4;          % 信息素挥发系数
Q = 1;              % 信息素增加强度系数
min_path = inf;
min_path_k = 0;
min_path_m = 0;
D = G2D(G);        % 邻接矩阵
N = size(D, 1);    % N表示问题的规模（像素个数）
a = 1;             % 小方格像素的边长
[ENDx, ENDy] = Getxy(END, r, a);   % 终止点坐标
eta = zeros(N, 1);               % 启发式信息，取为到目标点的直线距离的倒数

%% 初始化启发式信息矩阵
for i = 1: N
    [ix, iy] = Getxy(i, r, a);
    if i ~= END
        eta(i) = 1 / ((ix - ENDx)^2 + (iy - ENDy)^2)^0.5;
    else
        eta(i) = 100;
    end
end

%% 蚂蚁运动过程
routes = cell(K, M);        % 用细胞结构体存储每一代的每一只蚂蚁的爬行路线
path_len = zeros(K, M);     % 用矩阵存储每一代的每一只蚂蚁爬行路线长度
for k = 1:K
    for m = 1:M
        cur = START;        % 当前点为起始点
        path = [START];     % 路线
        len = 0;            % 路径长度
        tabu_list = ones(N);    % 禁忌表初始化
        tabu_list(START) = 0;   % 排除起始点
        adj_matrix = D;         % 邻接矩阵初始化
        [next_index, next_len] = GetNextNode(cur, adj_matrix, tabu_list);   % 下一步可以前往的节点
        % 蚂蚁没遇到事物或者陷入死胡同或觅食停止
        while cur ~= END && next_len >= 1
            % 轮盘赌法选择下一步怎么走
            prob = zeros(next_len);
            for i = 1:next_len
                prob(i) = (tau(cur, next_index(i))^alpha) * ((eta(next_index(i)))^beta);
            end
            sum_prob = sum(prob);
            prob = prob / sum_prob;     % 重新概率分布
            % 计算累积概率
            cum_prob(1) = prob(1);
            for i = 2:next_len
                cum_prob(i) = cum_prob(i - 1) + prob(i);
            end
            select = find(cum_prob >= rand);
            to_visit = next_index(select(1));
            % 状态更新和记录
            path = [path, to_visit];    % 路径
            len = len + adj_matrix(cur, to_visit);    % 路径长度
            cur = to_visit;             % 蚂蚁移到下一个节点
            for node = 1:N
                if tabu_list(node) == 0
                    adj_matrix(cur, node) = 0;
                    adj_matrix(node, cur) = 0;
                end
            end
            tabu_list(cur) = 0;         % 已访问的节点从禁忌表中删除
            [next_index, next_len] = GetNextNode(cur, adj_matrix, tabu_list);
        end
        % 记录每一代每一只蚂蚁的觅食路线和路线长度
        routes{k, m} = path;
        if path(end) == END
            path_len(k, m) = len;
            % 更新最短路径
            if len < min_path
                min_path_k = k;
                min_path_m = m;
                min_path = len;
            end
        else
            path_len(k, m) = 0;
        end
    end
    % 更新信息素 
    delta_tau = zeros(N);
    for m = 1:M
        if path_len(k, m)
            route = routes{k, m};    % 第m个蚂蚁的路线
            ts = length(route) - 1;  % 跳数
            len_tmp = path_len(k,  m);
            for s = 1:ts
                x = route(s);
                y = route(s + 1);
                delta_tau(x, y) = delta_tau(x, y) + Q / len_tmp;
                delta_tau(y, x) = delta_tau(y, x) + Q / len_tmp;
            end
        end
    end
    tau = (1 - rho) .* tau + delta_tau;
end

%% 绘图
plot_if = 1;        % 是否绘图的控制参数
if plot_if == 1
    min_path_len = zeros(K);
    for i = 1:K
        lens = path_len(i, :);
        not_zero = find(lens);
        lens = lens(not_zero);
        min_path_len(i) = min(lens);
    end
    figure(1);
    plot(min_path_len);
    hold on;
    grid on;
    title('收敛曲线变化趋势');
    xlabel('迭代次数'); ylabel('最小路径长度');
    
    plot_map(r, G, 2);
    title('蚂蚁运动轨迹');
    xlabel('坐标x'); ylabel('坐标y');
    route = routes{min_path_k, min_path_m};
    rout_len = length(route);
    rx = route;
    ry = route;
    for ii = 1:rout_len
        [rx(ii), ry(ii)] = Getxy(route(ii), r, a);
    end
    plot(rx, ry);
end

%% 绘制每次迭代蚂蚁最短路径
plot_if1 = 0;
if plot_if1 == 1
    plot_map(r, G, 3);
    for k = 1:K
        lens = path_len(k, :);
        min_len_index = min(lens);
        pos = find(path_len == min_len_index);
        m = pos(1);
        route = routes{k, m};
        rout_len = length(route);
        rx = route;
        ry = route;
        for ii = 1:rout_len
            [rx(ii), ry(ii)] = Getxy(route(ii), r, a);
        end
        plot(rx, ry);
        hold on;
    end  
end


%% 计算每个点到其余空闲点的距离
function D = G2D(G)
    L = size(G, 1);
    D = zeros(L * L, L * L);
    dx = [0 -1 -1 -1 0 1 1 1];
    dy = [1 1 0 -1 -1 -1 0 1];
    for i = 1:L
        for j = 1:L
            if G(i, j) == 0
                for m = 1:8
                    xx = i + dx(m); yy = j + dy(m);
                    if xx >= 1 && xx <= L && yy >= 1 && yy <= L && G(xx, yy) == 0
                        delta_x = abs(i - xx);
                        delta_y = abs(j - yy);
                        D((i - 1) * L + j, (xx - 1) * L + yy) = (delta_x + delta_y)^0.5;
                    end
                end
            end
        end
    end
end
   
%% 求每个点的x坐标和y坐标
function [x, y] = Getxy(num, size, a)
    x = a * (mod(num, size) - 0.5);   % 终止点横坐标
    if x == -0.5
        x = size - 0.5;
    end
    y = a * (size + 0.5 - ceil(num / size));   % 终止点纵坐标 
end

%% 求下一个可以到达的点，考虑邻接矩阵和信息素
function [next_index, next_len] = GetNextNode(cur, adj_matrix, tabu_list)
    next = adj_matrix(cur, :);
    next_temp_index = find(next);    % 返回非0元素的索引未考虑信息素
    for j = 1:length(next_temp_index)
        if tabu_list(next_temp_index(j)) == 0
            next(next_temp_index(j)) = 0;
        end
    end
    next_index = find(next);        % 考虑信息素后返回的下一个可行点的索引
    next_len = length(next_index);
end

%% 绘制地图
function plot_map(r, G, num)
    figure(num);
    axis([0 r 0 r]);
    for i = 1:r
        for j = 1:r
            x1 = j - 1; y1 = r - i;
            x2 = j; y2 = r - i;
            x3 = j; y3 = r - i + 1;
            x4 = j - 1; y4 = r - i + 1;
            if G(i, j) == 1
                fill([x1, x2, x3, x4], [y1, y2, y3, y4], [0.3 0.3 0.3]);
            else
                fill([x1, x2, x3, x4], [y1, y2, y3, y4], [1 1 1]);
            end
            hold on;
        end
    end
end