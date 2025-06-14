%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% important!!!!!!!
%%%% For the use of data please see the instructions: https://www.cancerrxgene.org/.
%%%% All the users should have the permission before usage.
%%%% Also see the paper: Iorio, Francesco, et al. "A landscape of pharmacogenomic interactions in cancer." Cell 166.3 (2016): 740-754.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load drug response data (Y)
drug = readtable("v17.3_fitted_dose_response.csv");

% Check dimensions
disp(size(drug));

% Subset drug data for "Gefitinib"
drug_sub = drug(strcmp(drug.DRUG_NAME, "Gefitinib"), :);
disp(size(drug_sub));

% Count drugs targeting "EGFR"
egfr_drugs = drug(strcmp(drug.PUTATIVE_TARGET, "EGFR"), "DRUG_NAME");
disp(tabulate(egfr_drugs.DRUG_NAME));

% Load RNA expression data (X)
RNA = readtable("sanger1018_brainarray_ensemblgene_rma.txt", 'ReadVariableNames', true);
% Check dimensions
disp(size(RNA));

% Load tissue type data (Z)
cov = readtable("Cell_Lines_Details_CSV.csv");

% Check dimensions
disp(size(cov));

% Count tissue types
disp(tabulate(cov{:,8}));

% Remove empty tissue entries
cov = cov(~cellfun('isempty', cov{:,8}), :);
disp(size(cov));

%% for the loop 
warning('off', 'all')
unique_drug_ids = unique(drug.DRUG_ID);
len = length(unique_drug_ids);
tau = 0.25;

estimates = cell(len, 1);
pvalues = cell(len, 1);

for i = 1:1
    try
        % Example: Analysis for one drug
        drugIDindex = unique_drug_ids(i); % Select drug by ID
        % Subset drug data for selected drug ID
        drugsub = drug(drug.DRUG_ID == drugIDindex, :);
        
        % Extract COSMIC IDs from each dataset
        y_COSMIC_ID = drugsub.COSMIC_ID;
        x_COSMIC_ID = RNA.Properties.VariableNames(2:end); % Exclude first column
        x_COSMIC_ID = regexprep(x_COSMIC_ID, '^x', ''); % Remove leading 'x'
        x_COSMIC_ID = str2double(x_COSMIC_ID);
        z_COSMIC_ID = cov{:,2}; % COSMIC IDs from cov
        
        % Find common IDs across Y, X, and Z
        xyz_ID = intersect(intersect(y_COSMIC_ID, x_COSMIC_ID), z_COSMIC_ID);
        % Extract Y (drug response)
        y_idx = ismember(drugsub.COSMIC_ID, xyz_ID);
        y = drugsub.LN_IC50(y_idx);
        
        % Extract X (gene expression)
        %x_idx = ismember(x_COSMIC_ID, xyz_ID);
        %x = RNA(:, [false x_idx]); % Keep only matched columns

        temRNA = RNA;
        temRNA(:, 1) = [];

        [tf, x_order] = ismember(xyz_ID, x_COSMIC_ID);
        x_order_valid = x_order(tf); % Remove non-matching indices
        x = temRNA(:, x_order_valid);

        % Convert table from n x p to p x n while keeping only numeric data
        tempTable = rows2vars(x);  % Convert rows to variables
        
        % Remove the first column which contains original variable names (cell array)
        tempTable.OriginalVariableNames = []; 
        
        % Convert to numeric matrix and transpose
        T = table2array(tempTable);
        
        % Extract Z (tissue types)
        [tf, idx] = ismember(string(xyz_ID), string(cov{:,2})); % Match IDs
        z_categorical = categorical(cov{idx(tf), 8}); % Convert matched Z values to categorical
        Z_dummy = dummyvar(z_categorical);

        z_levels = categories(z_categorical);
        z_dropped = z_levels{1};

        varNames = [z_dropped string(x.Properties.VariableNames)];

        Z_dummy(:, 1) = [];

        % Run the qr_add function (or any function that might cause an error)
        tic
        [est, pvalue] = qr_add(Z_dummy, T, y, tau, 1, 'test', 'kernel');  
        toc

    catch ME
        % Display error message and skip iteration
        disp(['Error in iteration ', num2str(i), ': ', ME.message])
        continue;
    end
end
