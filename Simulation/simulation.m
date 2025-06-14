%% setting 4.1

% please use the the data_generation1 to generate the data and 

n = [200 400 600];
rho = [0 0.25 0.5 0.75];
pz = [30 60 90];
px = 5000;
nrep = 1;
tau = [0.25 0.5];
all_time = [];
rep = 1;

for i = 1:3
    for j = 1:4
        for k = 1:3
            timing_results = zeros(1, 3);
            rho_str = num2str(rho(j), '%.15g');
            filename = sprintf("../simu_data/n_%d_rho_%s_pz_%d_rep_%d.csv", n(i), rho_str, pz(k), rep);
            disp(filename)
            dat = readmatrix(filename);
            dat = dat(2:end, 2:end);
            y = dat(:, 1);
            z = dat(:, 2:(pz(k) + 1));
            x = dat(:, (pz(k) + 2):end);

            tic
            [~, ~] = qr_add(z,x,y,0.25,1,'test','kernel'); 
            timing_results(1, 3) = toc;
            all_time = [all_time; timing_results rep i j k];
        end
    end
end


%% setting 4.2
% similar to before but use data_generation2


%% setting 4.3

% similar to before but use data_generation4


for i = 1:nweight
    timing_results = zeros(1, 3);
    filename = sprintf("../simu_data/data%d.csv", i);
    disp(filename)
    dat = readmatrix(filename);
    dat = dat(2:end, 2:end);
    y = dat(:, 1);
    z = dat(:, 2:(pz + 1));
    x = dat(:, (pz + 2):end);


    tic
    for i_inter=1:4000 
        idx = ((i_inter-1) * 5 + 1):(i_inter * 5);
        [beta, p]=qr_standard([z x(:,idx)], y, tau,'test','wald','method','interior');
    end
    timing_results(1, 1) = toc;


    tic
    for i_standard=1:4000 
        idx = ((i_standard-1) * 5 + 1):(i_standard * 5);
        [beta, p] = qr_standard([z x(:,idx)], y, tau,'test','kernel');
    end
    timing_results(1, 2) = toc;


    tic
    [est, pvalue] = qr_add(z, x, y, tau, 5, 'test', 'kernel'); 
    timing_results(1, 3) = toc;

    all_time = [all_time; timing_results];

end


%% setting 4.4

% similar to before but use data_generation3



%% across 

% please see qr_add_across
