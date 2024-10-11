% This script is used for the production of Tables 2 and the real data results in the manuscript
% Please modify the value of OBJ on Line 18 to produce a certain figure or table
% Re-running the simulation program may take a considerable amount of time
% If you want to directly load the previously generated data, please set load_if = 1

% INPUT ARGUMENT
% Please designate the figure or table to be reproduced
% e.g., 'T2' (Table 2), 'SC~', 'PS~' for the manuscript

% Please designate Re-running simulation program (load_if = 0)
% or directly load the previously generated data (load_if = 1)

% Please modify the address of the Transfer_fusion_regression folder
% on the local machine in the NewAddress at lines 20-21.

OBJ='PS~';
load_if = 0;
NewAddress = 'C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\Result';
NewAddress2 = 'C:\Users\张丰川\Documents\MATLAB\Transfer_fusion_regression\real_data\data';




% ------------------------------------Produce  Table 2----------------------------------------
if strcmp(OBJ,'T2') == 1
    if load_if ~= 1
        for which_index = 0:2 %%%%%%%%%%%%%%%%%%%%%%% 可修改 0:辅助域的个数，0是目标域本身
            %-----------------------------辅助域参数估计-------------------------------------%
            [simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_real_data(which_index);
            Aux = Auxiliary_heterogeneity_regression();
            if which_index == 0
                num_partitions_list = 3: 4; %2:5 3:4
                high_if = 1;
            else
                num_partitions_list = 5: 8; %51:100
                high_if = 1;
            end
            random_number_list = 51:100;
            [Results,results,Results_opt,results_opt,...
                initial,Results_list_opt,Results_single_opt,Class_summary, num_partitions_best,combine_principle_best,split_principle_best,random_number_best,partitions] =...
                Aux.Auxiliary_heterogeneous_regression(simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,...
                num_partitions_list,[],[],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],[],[],[],[],[]);
            % Save the estimation results for fusion method
            results_dir = fullfile(NewAddress, ['auxiliary']);
            save(fullfile(results_dir, ['results_' num2str(which_index) '.mat']), 'Results', 'results',...
                'Results_opt', 'results_opt', 'initial', 'Results_list_opt', 'Results_single_opt', 'Class_summary','num_partitions_best',...
                'combine_principle_best','split_principle_best','random_number_best','partitions');
        end
        %----------------------------------------------------------------------------------------------------------------%
        multi_if = 0; %%%%%%%%%%%%%%%%%%%%%%%可修改
        which_index = 2;%%%%%%%%%%%%%%%%%%%%%%%辅助域的个数
        num_partitions_list = 4; %5 4 3:4
        random_number_list = 51:100; %51:100
        high_if = 1;
        [simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_real_data(0);
        Trans = Transfer_heterogeneity_regression_pro();
        Auxiliary_dateset_number = which_index;
        Auxiliary = struct;
        Auxiliary(Auxiliary_dateset_number).position_provide = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        for i = 1:which_index   
            % 本地地址
            Auxiliary(i).position_provide = results_dir;
        end
        c_lambda2_start = 0.1;
        c_lambda2_end = 50;
        [Results,results,Results_opt,results_opt,...
            initial,Results_list_opt,Results_single_opt,Class_summary,num_partitions_best,combine_principle_best,split_principle_best,random_number_best,partitions] =...
            Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
            Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,num_partitions_list,[],...
            [],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],c_lambda2_start,c_lambda2_end,[],[],[],[]);
        % Save the estimation results for transfer method
        results_dir = fullfile(NewAddress, ['transfer']);
        save(fullfile(results_dir, ['results_' num2str(0) '.mat']), 'Results', 'results',...
                'Results_opt', 'results_opt', 'initial', 'Results_list_opt', 'Results_single_opt', 'Class_summary','num_partitions_best',...
                'combine_principle_best','split_principle_best','random_number_best','partitions');
        %----------------------------------------------------------------------------------------------------------------%
        multi_if = 1; %%%%%%%%%%%%%%%%%%%%%%%可修改
        which_index = 2;%%%%%%%%%%%%%%%%%%%%%%%辅助域的个数
        num_partitions_list = 4; %5 4 3:4
        random_number_list = 51:100; %51:100
        high_if = 1;
        [simulation_size, X_target, Y_target, sample_size, p, beta_target_real] = Data_generation_real_data(0);
        Trans = Transfer_heterogeneity_regression_pro();
        Auxiliary_dateset_number = which_index;
        Auxiliary = struct;
        Auxiliary(Auxiliary_dateset_number).position_provide = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        for i = 1:which_index
            % 本地地址
            Auxiliary(i).position_provide = results_dir;
        end
        c_lambda2_start = 0.1;
        c_lambda2_end = 50;
        [Results,results,Results_opt,results_opt,...
            initial,Results_list_opt,Results_single_opt,Class_summary,num_partitions_best,combine_principle_best,split_principle_best,random_number_best,partitions] =...
            Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
            Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,num_partitions_list,[],...
            [],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],c_lambda2_start,c_lambda2_end,[],[],[],[]);
        % Save the estimation results for transfer method
        results_dir = fullfile(NewAddress, ['transfer']);
        save(fullfile(results_dir, ['results_' num2str(1) '.mat']), 'Results', 'results',...
            'Results_opt', 'results_opt', 'initial', 'Results_list_opt', 'Results_single_opt', 'Class_summary','num_partitions_best',...
            'combine_principle_best','split_principle_best','random_number_best','partitions');
        %----------------------------------------------------------------------------------------------------------------%
        Table2 = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']));
        for i = 1: Class_summary.subgroup_number  
            Table2 = [Table2, Class_summary.coef{1, i}];
        end
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']));
        for i = 1: Class_summary.subgroup_number  
            Table2 = [Table2, Class_summary.coef{1, i}];
        end
        disp(Table2);
    else
        Table2 = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']));
        for i = 1: Class_summary.subgroup_number  
            Table2 = [Table2, Class_summary.coef{1, i}];
        end
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']));
        for i = 1: Class_summary.subgroup_number  
            Table2 = [Table2, Class_summary.coef{1, i}];
        end
        disp(Table2);
    end
end
























% ------------------------------------Produce  SC~----------------------------------------
if strcmp(OBJ,'SC~') == 1
    if load_if ~= 1
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']))
        beta_s1 = Class_summary.coef{1, 1};
        beta_s2 = Class_summary.coef{1, 2};
        beta_s3 = Class_summary.coef{1, 3};
        X1 = Class_summary.X{1, 1};
        X2 = Class_summary.X{1, 2};
        X3 = Class_summary.X{1, 3};
        Y1 = Class_summary.y{1, 1};
        Y2 = Class_summary.y{1, 2};
        Y3 = Class_summary.y{1, 3};
        group1_size = Class_summary.subgroup_sample_size{1};
        group2_size = Class_summary.subgroup_sample_size{2};
        group3_size = Class_summary.subgroup_sample_size{3};
        sample_size_or = group1_size + group2_size + group3_size;
        repeat_size = 157;
        matrix_st_nontransfer = zeros(sample_size_or, sample_size_or);
        matrix_st_random = zeros(sample_size_or, sample_size_or, repeat_size);
        matrix_st_panduan = zeros(1, repeat_size);
        index_1_or = Class_summary.index{1};
        index_2_or = Class_summary.index{2};
        index_3_or = Class_summary.index{3};
        for i = 1:(sample_size_or-1)
            for j = (i+1):sample_size_or
                if ismember(i, index_1_or)
                    which_1 = 1;
                elseif ismember(i, index_2_or)
                    which_1 = 2;
                else
                    which_1 = 3;
                end
                if ismember(j, index_1_or)
                    which_2 = 1;
                elseif ismember(j, index_2_or)
                    which_2 = 2;
                else
                    which_2 = 3;
                end
                if which_1 == which_2
                    matrix_st_nontransfer(i, j) = 1;
                else
                    matrix_st_nontransfer(i, j) = 0;
                end
            end
        end
        index_true = [index_1_or, index_2_or, index_3_or];
        Y = readtable(fullfile(NewAddress2, ['Data_target.csv'])).Y;
        X = readtable(fullfile(NewAddress2, ['Data_target.csv']));
        X_target_or = table2array(X(:,{'Data_surg', 'Data_radiation', 'Data_chemotherapy', 'Data_surg_rad_seq',...
            'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
            'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
        Y_target_or = Y;
        for t = 1 : repeat_size
            index_true(:,((t-1)*1 + 1):(1*t)) = [];
            map = zeros(length(index_true), 1);
            for i = 1:length(index_true)
                map(i, 1) = index_true(1, i);
            end
            X_target = X_target_or(index_true, :);
            Y_target = Y_target_or(index_true, :);
            [sample_size, p] = size(X_target);
            simulation_size = 1; % the number of repeat
            mu = 1;
            beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
            beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
            which_index = 0;
            try
                num_partitions_list = 3:4;
                high_if = 1;
                Aux = Auxiliary_heterogeneity_regression();
                random_number_list = 51:100;
                [Results, results, Results_opt, results_opt,...
                    initial, Results_list_opt, Results_single_opt, Class_summary] =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,...
                    num_partitions_list,[],[],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],[],[],[],[],[]);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                end
            catch
                num_partitions_list = 3:4;
                high_if = 0;
                Aux = Auxiliary_heterogeneity_regression();
                random_number_list = 51:100;
                [Results, results, Results_opt, results_opt,...
                    initial, Results_list_opt, Results_single_opt, Class_summary] =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,...
                    num_partitions_list,[],[],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],[],[],[],[],[]);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                end
            end
            for i = 1 : (sample_size-1)
                for j = (i+1):sample_size
                    if ismember(i, index_1)
                        which_1 = 1;
                    elseif ismember(i, index_2)
                        which_1 = 2;
                    else
                        which_1 = 3;
                    end
                    if ismember(j, index_1)
                        which_2 = 1;
                    elseif ismember(j, index_2)
                        which_2 = 2;
                    else
                        which_2 = 3;
                    end
                    if which_1 == which_2
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 1;
                    else
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 0;
                    end
                end
            end
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    sort_result = sort([map(i), map(j)]);
                    if matrix_st_random(sort_result(1), sort_result(2), t) == matrix_st_nontransfer(sort_result(1), sort_result(2))
                        matrix_st_panduan(1, t) = matrix_st_panduan(1, t) + 1;
                    end
                end
            end
            matrix_st_panduan(1, t) = matrix_st_panduan(1, t)/(0.5*sample_size*(sample_size-1));
            index_size;
            index_true = [index_1_or, index_2_or, index_3_or];
        end
        sc = matrix_st_panduan;
        sc_mean = mean(matrix_st_panduan);
        results_dir = fullfile(NewAddress, ['auxiliary']);
        save(fullfile(results_dir, ['sc.mat']), 'sc', 'sc_mean');
        %------------------------------------------------------------------------------------%
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']))
        beta_s1 = Class_summary.coef{1, 1};
        beta_s2 = Class_summary.coef{1, 2};
        beta_s3 = Class_summary.coef{1, 3};
        X1 = Class_summary.X{1, 1};
        X2 = Class_summary.X{1, 2};
        X3 = Class_summary.X{1, 3};
        Y1 = Class_summary.y{1, 1};
        Y2 = Class_summary.y{1, 2};
        Y3 = Class_summary.y{1, 3};
        group1_size = Class_summary.subgroup_sample_size{1};
        group2_size = Class_summary.subgroup_sample_size{2};
        group3_size = Class_summary.subgroup_sample_size{3};
        sample_size_or = group1_size + group2_size + group3_size;
        repeat_size = 157;
        matrix_st_transfer = zeros(sample_size_or, sample_size_or);
        matrix_st_random = zeros(sample_size_or, sample_size_or, repeat_size);
        matrix_st_panduan = zeros(1, repeat_size);
        index_1_or = Class_summary.index{1};
        index_2_or = Class_summary.index{2};
        index_3_or = Class_summary.index{3};
        for i = 1 : (sample_size_or-1)
            for j = (i+1) : sample_size_or
                if ismember(i, index_1_or)
                    which_1 = 1;
                elseif ismember(i, index_2_or)
                    which_1 = 2;
                else
                    which_1 = 3;
                end
                if ismember(j, index_1_or)
                    which_2 = 1;
                elseif ismember(j, index_2_or)
                    which_2 = 2;
                else
                    which_2 = 3;
                end
                if which_1 == which_2

                    matrix_st_transfer(i, j) = 1;
                else
                    matrix_st_transfer(i, j) = 0;
                end
            end
        end
        index_true = [index_1_or, index_2_or, index_3_or];
        Y = readtable(fullfile(NewAddress2, ['Data_target.csv'])).Y;
        X = readtable(fullfile(NewAddress2, ['Data_target.csv']));
        X_target_or = table2array(X(:,{'Data_surg', 'Data_radiation', 'Data_chemotherapy', 'Data_surg_rad_seq',...
            'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
            'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
        Y_target_or = Y;
        for t = 1:repeat_size
            index_true(:,((t-1)*1+1):(1*t)) = [];
            map = zeros(length(index_true), 1);
            for i = 1:length(index_true)
                map(i, 1) = index_true(1, i);
            end
            X_target = X_target_or(index_true, :);
            Y_target = Y_target_or(index_true, :);
            [sample_size, p] = size(X_target);
            simulation_size = 1; % the number of repeat
            mu = 1;
            beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
            beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
            try
                multi_if = 0; %%%%%%%%%%%%%%%%%%%%%%%可修改
                which_index = 2;
                num_partitions_list = 4;
                high_if = 1; %%%%1
                random_number_list = 51:100;
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                for i = 1:which_index
                    % 本地地址
                    results_dir = fullfile(NewAddress, ['auxiliary']);
                    Auxiliary(i).position_provide = results_dir;
                end
                c_lambda2_start = 0.1; %%%%固定最优lambda_2
                c_lambda2_end = 50; %%%%固定最优lambda_2
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary] = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,num_partitions_list,[],...
                    [],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],c_lambda2_start,c_lambda2_end,[],[],[],[]);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                end
            catch
                multi_if = 0; %%%%%%%%%%%%%%%%%%%%%%%可修改
                which_index = 2;
                num_partitions_list = 4;
                high_if = 0;
                random_number_list = 51:100;
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                for i = 1:which_index
                    % 本地地址
                    results_dir = fullfile(NewAddress, ['auxiliary']);
                    Auxiliary(i).position_provide = results_dir;
                end
                c_lambda2_start = 0.1; %%%%固定最优lambda_2
                c_lambda2_end = 50; %%%%固定最优lambda_2
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary] = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,num_partitions_list,[],...
                    [],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],c_lambda2_start,c_lambda2_end,[],[],[],[]);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                end
            end
            for i = 1 : (sample_size-1)
                for j = (i+1):sample_size
                    if ismember(i, index_1)
                        which_1 = 1;
                    elseif ismember(i, index_2)
                        which_1 = 2;
                    else
                        which_1 = 3;
                    end
                    if ismember(j, index_1)
                        which_2 = 1;
                    elseif ismember(j, index_2)
                        which_2 = 2;
                    else
                        which_2 = 3;
                    end
                    if which_1 == which_2
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 1;
                    else
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 0;
                    end
                end
            end
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    sort_result = sort([map(i), map(j)]);
                    if matrix_st_random(sort_result(1), sort_result(2), t) == matrix_st_transfer(sort_result(1), sort_result(2))
                        matrix_st_panduan(1, t) = matrix_st_panduan(1, t) + 1;
                    end
                end
            end
            matrix_st_panduan(1, t) = matrix_st_panduan(1, t)/(0.5*sample_size*(sample_size-1));
            index_true = [index_1_or, index_2_or, index_3_or];
        end
        sc = matrix_st_panduan;
        sc_mean = mean(matrix_st_panduan);
        results_dir = fullfile(NewAddress, ['transfer']);
        save(fullfile(results_dir, ['sc.mat']), 'sc', 'sc_mean');
        SC = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['sc.mat']));
        SC = [SC, sc_mean];
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['sc.mat']));
        SC = [SC, sc_mean];
        disp(SC)
    else
        SC = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['sc.mat']));
        SC = [SC, sc_mean];
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['sc.mat']));
        SC = [SC, sc_mean];
        disp(SC);
    end
end
































% ------------------------------------Produce  PS~----------------------------------------
if strcmp(OBJ,'PS~') == 1
    if load_if ~= 1
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']));
        beta_s1 = Class_summary.coef{1, 1};
        beta_s2 = Class_summary.coef{1, 2};
        beta_s3 = Class_summary.coef{1, 3};
        X1 = Class_summary.X{1, 1};
        X2 = Class_summary.X{1, 2};
        X3 = Class_summary.X{1, 3};
        Y1 = Class_summary.y{1, 1};
        Y2 = Class_summary.y{1, 2};
        Y3 = Class_summary.y{1, 3};
        group1_size = Class_summary.subgroup_sample_size{1};
        group2_size = Class_summary.subgroup_sample_size{2};
        group3_size = Class_summary.subgroup_sample_size{3};
        sample_size_or = group1_size + group2_size + group3_size;
        repeat_size = 157; 
        delate_number = 1;
        real_evaluate = [];
        matrix_st_nontransfer = zeros(sample_size_or, sample_size_or);
        matrix_st_random = zeros(sample_size_or, sample_size_or, repeat_size);
        matrix_st_panduan = zeros(1, repeat_size);
        index_1_or = Class_summary.index{1};
        index_2_or = Class_summary.index{2};
        index_3_or = Class_summary.index{3};
        for i = 1:(sample_size_or-1)
            for j = (i+1):sample_size_or
                if ismember(i, index_1_or)
                    which_1 = 1;
                elseif ismember(i, index_2_or)
                    which_1 = 2;
                else
                    which_1 = 3;
                end
                if ismember(j, index_1_or)
                    which_2 = 1;
                elseif ismember(j, index_2_or)
                    which_2 = 2;
                else
                    which_2 = 3;
                end
                if which_1 == which_2
                    matrix_st_nontransfer(i, j) = 1;
                else
                    matrix_st_nontransfer(i, j) = 0;
                end
            end
        end
        map_or = zeros(1, sample_size_or);
        for i = 1:sample_size_or
            if ismember(i, index_1_or)
                map_or(1, i) = 1;
            elseif ismember(i, index_2_or)
                map_or(1, i) = 2;
            else
                map_or(1, i) = 3;
            end
        end
        matrix_st_panduan2 = zeros(1, repeat_size);
        index_true = [index_1_or, index_2_or, index_3_or];
        index_true_or = [index_1_or, index_2_or, index_3_or];
        Y = readtable(fullfile(NewAddress2, ['Data_target.csv'])).Y;
        X = readtable(fullfile(NewAddress2, ['Data_target.csv']));
        X_target_or = table2array(X(:,{'Data_surg', 'Data_radiation', 'Data_chemotherapy', 'Data_surg_rad_seq',...
            'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
            'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
        Y_target_or = Y;
        for t = 1 : repeat_size
            index_true(:,((t-1)*delate_number + 1):(delate_number*t)) = [];
            map = zeros(length(index_true), 1);
            for i = 1:length(index_true)
                map(i, 1) = index_true(1, i);
            end
            X_target = X_target_or(index_true, :);
            Y_target = Y_target_or(index_true, :);
            [sample_size, p] = size(X_target);
            simulation_size = 1; % the number of repeat
            mu = 1;
            beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
            beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
            which_index = 0;
            try
                num_partitions_list = 3:4;
                high_if = 1;
                Aux = Auxiliary_heterogeneity_regression();
                random_number_list = 51:100; %51:100 87
                [Results, results, Results_opt, results_opt,...
                    initial, Results_list_opt, Results_single_opt, Class_summary] =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,...
                    num_partitions_list,[],[],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],[],[],[],real_evaluate,t);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                    eval([ ['coef_' strrep(num2str(i), '.', '')] ' =  Class_summary.coef{1, i};'])
                end
            catch
                num_partitions_list = 3:4;
                high_if = 0;
                Aux = Auxiliary_heterogeneity_regression();
                random_number_list = 51:100;  %51:100 87
                [Results, results, Results_opt, results_opt,...
                    initial, Results_list_opt, Results_single_opt, Class_summary] =  Aux.Auxiliary_heterogeneous_regression(simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,...
                    num_partitions_list,[],[],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],[],[],[],real_evaluate,t);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                    eval([ ['coef_' strrep(num2str(i), '.', '')] ' =  Class_summary.coef{1, i};'])
                end
            end
            for i = 1 : (sample_size-1)
                for j = (i+1):sample_size
                    if ismember(i, index_1)
                        which_1 = 1;
                    elseif ismember(i, index_2)
                        which_1 = 2;
                    else
                        which_1 = 3;
                    end
                    if ismember(j, index_1)
                        which_2 = 1;
                    elseif ismember(j, index_2)
                        which_2 = 2;
                    else
                        which_2 = 3;
                    end
                    if which_1 == which_2
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 1;
                    else
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 0;
                    end
                end
            end
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    sort_result = sort([map(i), map(j)]);
                    if matrix_st_random(sort_result(1), sort_result(2), t) == matrix_st_nontransfer(sort_result(1), sort_result(2))
                        matrix_st_panduan(1, t) = matrix_st_panduan(1, t) + 1;
                    end
                end
            end
            matrix_st_panduan(1, t) = matrix_st_panduan(1, t)/(0.5*sample_size*(sample_size-1));
            map_match = zeros(1, index_size);
            map_match_2 = zeros(1, 3);
            which_model = zeros(index_size, 3); 
            for i = 1: index_size
                for j = 1: 3 
                    which_model(i, j) = (length(eval(['index_', num2str(j), '_or']))/sample_size_or)*numel(intersect(eval(['index_true(1, index_' strrep(num2str(i), '.', '') ')']),...
                        eval(['index_', num2str(j), '_or'])))/length(eval(['index_', num2str(j), '_or']));%
                end
            end
            for j = 1:3
                [~, w] = max(which_model(:, j));
                map_match_2(1, j) = w;
            end
            map_p = zeros(delate_number, index_size);
            map_estimate = zeros(delate_number, 1);
            for j = 1 : delate_number
                for i = 1 : index_size
                    X_target_part = X_target(eval(['index_' strrep(num2str(i), '.', '')]), :);
                    Y_target_part = Y_target(eval(['index_' strrep(num2str(i), '.', '')]), :);
                    [sample_size, p] = size(X_target);
                    epsilon = Y_target_part - X_target_part*eval(['coef_' strrep(num2str(i), '.', '')]);
                    index_delate = ((t-1)*delate_number + 1):(delate_number*t);
                    epsilon_t = Y_target_or(index_true_or(1,index_delate(j)),:) -...
                        X_target_or(index_true_or(1,index_delate(j)),:)*eval(['coef_' strrep(num2str(i), '.', '')]);
                    map_p(j, i) = (length(eval(['index_' strrep(num2str(i), '.', '')]))/sample_size)...
                        *normpdf(epsilon_t, mean(epsilon), std(epsilon));
                end
                [~,map_estimate(j, 1)] = max(map_p);
            end
            index_delate = ((t-1)*delate_number + 1):(delate_number*t);
            for i = 1:delate_number
                if map_match_2(map_or(1, index_true_or(index_delate(i)))) == map_estimate(i,1)
                    matrix_st_panduan2(1, t) = 1;
                end
            end
            index_size;
            t;
            index_true = [index_1_or, index_2_or, index_3_or];
        end
        mean(matrix_st_panduan2);
        ps = matrix_st_panduan2;
        ps_mean = mean(matrix_st_panduan2);
        results_dir = fullfile(NewAddress, ['auxiliary']);
        save(fullfile(results_dir, ['ps.mat']), 'ps', 'ps_mean');
        %------------------------------------------------------------------------------%
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['results_' num2str(0) '.mat']));
        beta_s1 = Class_summary.coef{1, 1};
        beta_s2 = Class_summary.coef{1, 2};
        beta_s3 = Class_summary.coef{1, 3};
        X1 = Class_summary.X{1, 1};
        X2 = Class_summary.X{1, 2};
        X3 = Class_summary.X{1, 3};
        Y1 = Class_summary.y{1, 1};
        Y2 = Class_summary.y{1, 2};
        Y3 = Class_summary.y{1, 3};
        group1_size = Class_summary.subgroup_sample_size{1};
        group2_size = Class_summary.subgroup_sample_size{2};
        group3_size = Class_summary.subgroup_sample_size{3};
        sample_size_or = group1_size + group2_size + group3_size;
        repeat_size = 157;
        delate_number = 1;
        real_evaluate = [];
        matrix_st_transfer = zeros(sample_size_or, sample_size_or);
        matrix_st_random = zeros(sample_size_or, sample_size_or, repeat_size);
        matrix_st_panduan = zeros(1, repeat_size);
        coef_1_or = Class_summary.coef{1, 1};
        coef_2_or = Class_summary.coef{1, 2};
        coef_3_or = Class_summary.coef{1, 3};
        index_1_or = Class_summary.index{1};
        index_2_or = Class_summary.index{2};
        index_3_or = Class_summary.index{3};
        for i = 1 : (sample_size_or-1)
            for j = (i+1) : sample_size_or
                if ismember(i, index_1_or)
                    which_1 = 1;
                elseif ismember(i, index_2_or)
                    which_1 = 2;
                else
                    which_1 = 3;
                end
                if ismember(j, index_1_or)
                    which_2 = 1;
                elseif ismember(j, index_2_or)
                    which_2 = 2;
                else
                    which_2 = 3;
                end
                if which_1 == which_2
                    matrix_st_transfer(i, j) = 1;
                else
                    matrix_st_transfer(i, j) = 0;
                end
            end
        end
        map_or = zeros(1, sample_size_or);
        for i = 1:sample_size_or
            if ismember(i, index_1_or)
                map_or(1, i) = 1;
            elseif ismember(i, index_2_or)
                map_or(1, i) = 2;
            else
                map_or(1, i) = 3;
            end
        end
        matrix_st_panduan2 = zeros(1, repeat_size);
        index_true = [index_1_or, index_2_or, index_3_or];
        index_true_or = [index_1_or, index_2_or, index_3_or];
        Y = readtable(fullfile(NewAddress2, ['Data_target.csv'])).Y;
        X = readtable(fullfile(NewAddress2, ['Data_target.csv']));
        X_target_or = table2array(X(:,{'Data_surg', 'Data_radiation', 'Data_chemotherapy', 'Data_surg_rad_seq',...
            'Data_age', 'Data_sex', 'Data_site_detail', 'Data_summary_stage', 'Data_time_diag_treat',...
            'Data_first_if', 'Data_total_number', 'Data_marital_status', 'Data_median_income'}));
        Y_target_or = Y;
        for t = 1:repeat_size
            index_true(:,((t-1)*delate_number+1):(delate_number*t)) = [];
            map = zeros(length(index_true), 1);
            for i = 1:length(index_true)
                map(i, 1) = index_true(1, i);
            end
            X_target = X_target_or(index_true, :);
            Y_target = Y_target_or(index_true, :);
            [sample_size, p] = size(X_target);
            simulation_size = 1; % the number of repeat
            mu = 1;
            beta_target_real1 = repmat([mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real2 = repmat([-0.5*mu*ones(floor(0.5*p),1);0.5*mu*ones((p-floor(0.5*p)),1)], (sample_size-2*floor((1/3)*sample_size)), 1);
            beta_target_real3 = repmat([-mu*ones(p,1)], floor((1/3)*sample_size), 1);
            beta_target_real = [beta_target_real1; beta_target_real2; beta_target_real3];
            try
                multi_if = 0; 
                which_index = 2;
                num_partitions_list =  4;
                high_if = 1; 
                random_number_list = 51:100;
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                for i = 1:which_index
                    % 本地地址
                    results_dir = fullfile(NewAddress, ['auxiliary']);
                    Auxiliary(i).position_provide = results_dir;
                end
                c_lambda2_start = 0.1;
                c_lambda2_end = 5;
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary] = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,num_partitions_list,[],...
                    [],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],c_lambda2_start,c_lambda2_end,[],[],real_evaluate,t);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                    eval([ ['coef_' strrep(num2str(i), '.', '')] ' =  Class_summary.coef{1, i};'])
                end
            catch
                multi_if = 0; 
                which_index = 2;
                num_partitions_list = 4;
                high_if = 0;
                random_number_list = 51:100;
                Trans = Transfer_heterogeneity_regression_pro();
                Auxiliary_dateset_number = which_index;
                Auxiliary = struct;
                Auxiliary(Auxiliary_dateset_number).position_provide = [];
                for i = 1:which_index
                    % 本地地址
                    results_dir = fullfile(NewAddress, ['auxiliary']);
                    Auxiliary(i).position_provide = results_dir;
                end
                c_lambda2_start = 0.1;
                c_lambda2_end = 5;
                [Results,results,Results_opt,results_opt,...
                    initial,Results_list_opt,Results_single_opt,Class_summary] = Trans.Transfer_heterogeneous_regression_pro(Auxiliary_dateset_number,...
                    Auxiliary,multi_if,simulation_size,[],X_target,Y_target,sample_size,p,beta_target_real,num_partitions_list,[],...
                    [],[],[],random_number_list,[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],high_if,[],[],[],[],[],[],[],c_lambda2_start,c_lambda2_end,[],[],real_evaluate,t);
                index_size = length(Class_summary.index);
                for i = 1:index_size
                    eval([ ['index_' strrep(num2str(i), '.', '')] ' =  Class_summary.index{i};'])
                    eval([ ['coef_' strrep(num2str(i), '.', '')] ' =  Class_summary.coef{1, i};'])
                end
            end
            for i = 1 : (sample_size - 1)
                for j = (i+1):sample_size
                    if ismember(i, index_1)
                        which_1 = 1;
                    elseif ismember(i, index_2)
                        which_1 = 2;
                    elseif ismember(i, index_3)
                        which_1 = 3;
                    else
                        which_1 = 4;
                    end
                    if ismember(j, index_1)
                        which_2 = 1;
                    elseif ismember(j, index_2)
                        which_2 = 2;
                    elseif ismember(i, index_3)
                        which_2 = 3;
                    else
                        which_2 = 4;
                    end
                    if which_1 == which_2
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 1;
                    else
                        sort_result = sort([map(i),map(j)]);
                        matrix_st_random(sort_result(1), sort_result(2), t) = 0;
                    end
                end
            end
            for i = 1:(sample_size-1)
                for j = (i+1):sample_size
                    sort_result = sort([map(i), map(j)]);
                    if matrix_st_random(sort_result(1), sort_result(2), t) == matrix_st_transfer(sort_result(1), sort_result(2))
                        matrix_st_panduan(1, t) = matrix_st_panduan(1, t) + 1;
                    end
                end
            end
            matrix_st_panduan(1, t) = matrix_st_panduan(1, t)/(0.5*sample_size*(sample_size-1));
            map_match = zeros(1, index_size);
            map_match_2 = zeros(1, 3);
            which_model = zeros(index_size, 3); 
            for i = 1: index_size
                for j = 1: 3 
                    which_model(i, j) = (length(eval(['index_', num2str(j), '_or']))/sample_size_or)*numel(intersect(eval(['index_true(1, index_' strrep(num2str(i), '.', '') ')']),...
                        eval(['index_', num2str(j), '_or'])))/length(eval(['index_', num2str(j), '_or']));
                end
            end
            for j = 1:3
                [~, w] = max(which_model(:, j));
                map_match_2(1, j) = w;
            end
            map_p = zeros(delate_number, index_size);
            map_estimate = zeros(delate_number, 1);
            for j = 1 : delate_number
                for i = 1 : index_size
                    X_target_part = X_target(eval(['index_' strrep(num2str(i), '.', '')]), :);
                    Y_target_part = Y_target(eval(['index_' strrep(num2str(i), '.', '')]), :);
                    [sample_size, p] = size(X_target);
                    epsilon = Y_target_part - X_target_part*eval(['coef_' strrep(num2str(i), '.', '')]);
                    index_delate = ((t-1)*delate_number + 1):(delate_number*t);
                    epsilon_t = Y_target_or(index_true_or(1,index_delate(j)),:) -...
                        X_target_or(index_true_or(1,index_delate(j)),:)*eval(['coef_' strrep(num2str(i), '.', '')]);
                    map_p(j, i) = (length(eval(['index_' strrep(num2str(i), '.', '')]))/sample_size)...
                        *normpdf(epsilon_t, mean(epsilon), std(epsilon));
                end
                [~,map_estimate(j, 1)] = max(map_p);
            end
            index_delate = ((t-1)*delate_number + 1):(delate_number*t);
            for i = 1:delate_number
                if map_match_2(map_or(1, index_true_or(index_delate(i)))) == map_estimate(i,1)
                    matrix_st_panduan2(1, t) = 1;
                end
            end
            index_size;
            t;
            index_true = [index_1_or, index_2_or, index_3_or];
        end
        mean(matrix_st_panduan2);
        ps = matrix_st_panduan2;
        ps_mean = mean(matrix_st_panduan2);
        results_dir = fullfile(NewAddress, ['transfer']);
        save(fullfile(results_dir, ['ps.mat']), 'ps', 'ps_mean');
        PS = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['ps.mat']));
        PS = [PS, ps_mean];
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['ps.mat']));
        PS = [PS, ps_mean];
        disp(PS)
    else
        PS = [];
        results_dir = fullfile(NewAddress, ['auxiliary']);
        load(fullfile(results_dir, ['ps.mat']));
        PS = [PS, ps_mean];
        results_dir = fullfile(NewAddress, ['transfer']);
        load(fullfile(results_dir, ['ps.mat']));
        PS = [PS, ps_mean];
        disp(PS)
    end
end