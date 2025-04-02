%% Data Loading and Initial Processing Script
% Description: This script loads cell culture process data from two Excel files,
%              selects necessary columns, and groups the data by Batch ID.
% Author: Dong-Hyuk Choi
% Date: 2025-04-02
%
% Usage:
% 1. Create a subfolder named 'Data' in the same directory as this script file.
% 2. Place the required Excel files inside the 'Data' folder:
%    - 'ELN_Excel data sheet_Bio141 to Bio148.xlsx'
%    - 'processing_ELN_Excel data sheet_Bio141 to Bio148.xlsx'
% 3. Run the script.file_path = '';
%%
file_path = '';
file_name = 'ELN_Excel data sheet_Bio141 to Bio148.xlsx';
sheet_name = 'All data';

% Combine the file path and file name
full_path = fullfile(file_path, file_name);

% Read the Excel file from the specified sheet, starting from the 18th row
df_tmp = readtable(full_path, 'Sheet', sheet_name, 'Range', '18:129', 'ReadVariableNames', true);


% Select the specific columns from the imported data∑
selected_columns = {'BatchID', 'SampleID', 'ViableCells', 'Age_h_', 'TotalVolume', 'BaseVolume', 'EffFeedBVolume', 'GlucFeedVolume'};
df = df_tmp(:, selected_columns);
% Define the row range to select
start_row = 1;
end_row = height(df); % Select all rows from the imported range
% Select the rows from the table
viable_cells = df(start_row:end_row, :);
% Display the imported data
disp(viable_cells);



file_path_meta = "";
file_name_meta = "processing_ELN_Excel data sheet_Bio141 to Bio148.xlsx";
sheet_name_meta = "Amount of uptake and secretion";
sheet_name_gr = "Growth_Rate";
full_path_meta = fullfile(file_path_meta, file_name_meta);

% Read the Excel file from the specified sheet, starting from the 18th row
df_meta = readtable(full_path_meta, 'Sheet', sheet_name_meta, 'Range', '3:129', 'ReadVariableNames', true);
df_gr = readtable(full_path_meta, 'Sheet', sheet_name_gr, 'Range', '1:129', 'ReadVariableNames', true);

% Remove trailing spaces from column names
df_meta.Properties.VariableNames = strtrim(df_meta.Properties.VariableNames);

% Define the row range to select
start_row = 1;
end_row = 123;

% Select the rows from the table
metabolite = df_meta(start_row:end_row, :);

% Group the data by 'Batch ID' for df_meta and 'Batch_ID' for df_gr
batch_dfs_meta = containers.Map;
batch_dfs_raw = containers.Map;

unique_batch_ids_meta = unique(metabolite.BatchID);
unique_batch_ids_gr = unique(df_gr.Batch_ID);

for i = 1:length(unique_batch_ids_meta)
    batch_id = unique_batch_ids_meta{i};
    batch_dfs_meta(batch_id) = metabolite(strcmp(metabolite.BatchID, batch_id), :);
end

for i = 1:length(unique_batch_ids_gr)
    batch_id = unique_batch_ids_gr{i};
    batch_dfs_raw(batch_id) = df_gr(strcmp(df_gr.Batch_ID, batch_id), :);
end

% Display the grouped data for verification
disp('Grouped data for batch_dfs_meta:');
disp(batch_dfs_meta);

disp('Grouped data for batch_dfs_raw:');
disp(batch_dfs_raw);
