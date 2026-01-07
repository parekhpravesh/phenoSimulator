function clusterinfo = get_clusterinfo(fid, iid, varargin)
%% Function to return clusters in the data, given family and individual IDs
% Adapted from FEMA_parse_family function from the FEMA package:
% https://github.com/cmig-research-group/cmig_tools
%% Mandatory inputs:
% fid:          [n x 1]     cell type defining the family IDs per observation
%                           (top level; clusters are defined at this level)
% 
% iid:          [n x 1]     cell type defining the subject IDs per observation
%                           (individuals are nested within families)
% 
%% Optional inputs: should be provided as name-pair values
% age:          [n x 1]     age at measurement; if provided, will be used
%                           to sort data internally within the cluster
% 
% GRM:          [q x q]     matrix of genetic relatedness; this matrix
%                           should be ordered as unique(iid, 'stable')
% 
%% Output:
% clusterinfo:  [1 x c]     cell type containing c clusters; every cluster
%                           in the data is an entry in the cell array,
%                           containing structures inside corresponding to 
%                           random effects design matrix:
%                               * V_F
%                               * V_S
%                               * V_A (if GRM is provided)
%                               * V_D (V_D = V_A ^ 2, if GRM is provided)

%% Check inputs
if ~exist('fid', 'var') || isempty(fid)
    error('Family ID / top level clustering info needs to be provided');
end

if ~exist('iid', 'var') || isempty(iid)
    error('Subject ID needs to be provided');
end

%% Optional inputs
p = inputParser;
addParameter(p, 'age', [], @isnumeric);
addParameter(p, 'GRM', [], @isnumeric);
parse(p, varargin{:});

age = p.Results.age;
GRM = p.Results.GRM;

if isempty(age)
    age = zeros(length(iid),1);
end

%% Main module
[fid_list, ~, IC_fam]  = unique(fid, 'stable'); nfam  = length(fid_list);
[iid_list, ~, IC_subj] = unique(iid, 'stable'); nsubj = length(iid_list);

subj_famtypevec = NaN(length(iid), 1);

% Find all timepoints for each subject
jvecs_subj = cell(1,nsubj);
for subji = 1:length(iid_list)
    jvec_subj = find(IC_subj==subji);
    jvec_subj = reshape(jvec_subj, [1 numel(jvec_subj)]);
    jvecs_subj{subji} = jvec_subj;
end

% Should classify each family as one of multiple types: 
% # of subjects, # of observations each: sort subjects & observations in
% consistent manner
clusterinfo     = cell(1, nfam);
famtypelist     = {}; 
famtypevec      = NaN(1, nfam);
nfmemvec        = NaN(1, nfam);
count           = 1;
str_famtypelist = {};

for fi = 1:nfam
    % Identify all observations (rows) for a given family
    jvec_fam = find(IC_fam == fi);
    jvec_fam = reshape(jvec_fam, [1 numel(jvec_fam)]);

    % Subject number for each observation
    subj_fam = IC_subj(jvec_fam);

    % List of unique subjects
    subj_unique = unique(IC_subj(jvec_fam));

    % Frequency
    freq_unique = NaN(size(subj_unique));
    for j = 1:length(subj_unique)
        freq_unique(j) = sum(IC_subj(jvec_fam) == subj_unique(j));
    end
    [~, si]         = sort(freq_unique);
    freq_unique     = freq_unique(si);
    str_freq_unique = num2str(freq_unique');
    subj_unique     = subj_unique(si); % Re-order subjects to canonical form
    subji_jvec      = NaN(size(subj_fam));
    for subjii = 1:length(subj_unique)
        ivec_tmp = IC_subj(jvec_fam)==subj_unique(subjii);
        subji_jvec(ivec_tmp) = subjii;
    end
    
    [~, si]  = sortrows([subji_jvec age(jvec_fam)]);
    jvec_fam = jvec_fam(si); % Re-order to cannonical form
    ivec     = find(strcmpi(str_freq_unique, str_famtypelist));
    if isempty(ivec)
        famtypelist             = cat(2,famtypelist,{freq_unique});
        str_famtypelist{count}  = str_freq_unique; %#ok<AGROW>
        count                   = count + 1;
        famtypevec(fi)          = length(famtypelist);
    else
        if isscalar(ivec)
            famtypevec(fi) = ivec(1);
        end
    end
    subj_famtypevec(jvec_fam) = famtypevec(fi);
    nfmemvec(fi)              = length(jvec_fam);

    % Initialize structure entries
    V_E = eye(length(jvec_fam),   length(jvec_fam));
    V_F = true(length(jvec_fam),  length(jvec_fam));
    V_S = false(length(jvec_fam), length(jvec_fam));

    % Reorder GRM, if required
    if ~isempty(GRM)
        V_A = GRM(subj_fam(si), subj_fam(si)); % Reorder to canonical form
        V_D = V_A.^2;
    else
        V_A = [];
        V_D = [];
    end

    % Create V_S
    for ji = 1:length(jvec_fam)
        jvec_tmp = jvecs_subj{IC_subj(jvec_fam(ji))};
        ivec_tmp = ismember(jvec_fam,jvec_tmp);
        V_S(ji,ivec_tmp) = true;
    end

    % Create entry in clusterinfo
    clusterinfo{fi}         = struct('V_E', V_E, 'V_S', V_S, 'V_F', V_F, ...
                                     'V_A', V_A, 'V_D', V_D, 'jvec_fam', jvec_fam);
    clusterinfo{fi}.famtype = famtypevec(fi);
end