% GBM Feature Study

%%%%%%%% Flair Image - small tumoral tissue

filename = '/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX FLAIR.nii'; 
img = load_nii(filename);
view_nii(img); % Viewing flair with edema

%%%%%% Post T1 Image - most accurate to see enhancement of tumor

filename2 = '/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX T1 POST.nii';
img2 = load_nii(filename2);
view_nii(img2); % Viewing T1 post

%%%%%% T1 Image - important for hemorrage

filename3 = '/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX T1.nii';
img3 = load_nii(filename3);
view_nii(img3); % Viewing T1 

% Creating a transformation matrix, which will let the reference image turn 30 degree counter-clockwise on XY plane
T = [cos(pi/6) -sin(pi/6) 0; sin(pi/6) cos(pi/6) 0; 0 0 1];

% Get old_xyz information from reference image:
rl = load_untouch_nii('/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX FLAIR.nii');
old_xyz = [rl.hdr.hist.srow_x(1:3);rl.hdr.hist.srow_y(1:3);rl.hdr.hist.srow_z(1:3)];

% Apply transformation matrix, and save new_xyz into a new image:
new_xyz = 1 * old_xyz;
rl.hdr.hist.srow_x(1:3) = new_xyz(1,:);
rl.hdr.hist.srow_y(1:3) = new_xyz(2,:);
rl.hdr.hist.srow_z(1:3) = new_xyz(3,:);
save_untouch_nii(rl, '/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX FLAIR_rl30.nii');

% In order to view this image, you need to reslice it:
reslice_nii('/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX FLAIR_rl30.nii', 'rl30b.nii');
% Now, you can load and view the rotated image:
rl30b = load_nii('rl30b.nii');
view_nii(rl30b);

% Find all images of T1 post and label in our TCGA/TCIA dataset
cd '/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/'
dirData_T1post = dir('**/* POST.nii');
dirData_label = dir('**/*-label.nii');
total_glcm_stats = {};
sample_id = {};

% running for all T1 post and label combos 
for k = 1:size(dirData_T1post)
    T1_post = load_nii(strcat(dirData_T1post(k).folder, '/',dirData_T1post(k).name)); % extract T1_post image
    EEG = load_nii(strcat(dirData_label(k).folder,'/', dirData_label(k).name)); % extract label image
    for i=1:1:size(T1_post.img,3) 
        T1_post_temp = T1_post.img(:,:,i); % obtain slice info for T1 post
        tumor_temp = EEG.img(:,:,i); % obtain slice info for label
        sample_id = T1_post.fileprefix; % find sample ID for current image 
        sample_id = regexp(sample_id,'RPZ.\w*.\w*','match'); % regular expression to parse out name
        masked_tumor_Image = bsxfun(@times, T1_post_temp, cast(tumor_temp, 'like', T1_post_temp)); % create the tumor-image hybrid
        glcm = graycomatrix(masked_tumor_Image,'NumLevels',9,'GrayLimits',[]); % computing graycomatrix values
        stats = graycoprops(glcm,'all'); % find graycoprops - used downstream for analysis
        glcm_stats(i,:)= struct2array(stats); % store glcm values in a structure
    end
    total_glcm_stats(k,:) = {glcm_stats}; % build structure of all glcm characters of all samples
    sample_id_names(k) = sample_id; % keeps tab on which samples are being refernced 
end

% Since the output label only marked the tumor, there were regions from the
% T1 post that did not produce useful GLCM features. We removed the NaN and
% 0 values from our dataset through this loop
for i=1:24
    %idx = find(copy_total_glcm_stats{i}(:,1));
        a = total_glcm_stats(i);
        b = cellfun(@isnan, a, 'UniformOutput',false); 
        idx = find(b{:,1});
        %for j = 1:size(total_glcm_stats{i}, 2)
            %for k = 1:size(a{1}size(,2)
            %idx = union(idx, (find(b{:,1})));
            %end
        %end
        c = a{1};
        size_b = idx-size(b{1},1);
        c(size_b, :) = [];
        total_glcm_stats_tmp{i,:} = {c};
end

% Temporary file creation of GLCM features 
%xlswrite('Updated_GLCM_features.xls',[total_glcm_stats_tmp(1), total_glcm_stats_tmp(2), total_glcm_stats_tmp(3), total_glcm_stats_tmp(4), total_glcm_stats_tmp(5), total_glcm_stats_tmp(6), total_glcm_stats_tmp(7), total_glcm_stats_tmp(8), total_glcm_stats_tmp(9), total_glcm_stats_tmp(10), total_glcm_stats_tmp(11), total_glcm_stats_tmp(12), total_glcm_stats_tmp(13),total_glcm_stats_tmp(14), total_glcm_stats_tmp(15),total_glcm_stats_tmp(16),total_glcm_stats_tmp(17),total_glcm_stats_tmp(18),total_glcm_stats_tmp(19),total_glcm_stats_tmp(20),total_glcm_stats_tmp(21),total_glcm_stats_tmp(22),total_glcm_stats_tmp(23),total_glcm_stats_tmp(24)])

%read in survival data in specific format
formatSpec = '%C%C%f%f%C%C%C%f%f%f';
opts = detectImportOptions('TCGA Clinical Data_ID_Reduced.csv');
opts = setvartype(opts,3,'double');
opts = setvaropts(opts,3,'DecimalSeparator',',');
opts.ExtraColumnsRule = 'ignore'; % There's an extra tab at the end of some of the rows, this prevents that from being imported.
clinical_info = readtable('TCGA Clinical Data_ID_Reduced.csv',opts);

% build table of GLCM features
col_names = {'Contrast','Correlation','Energy','Homogeneity'};
for i =1:24
    total_glcm_stats_tmp_table = cell2mat(total_glcm_stats_tmp{i});
    total_glcm_stats_tmp_table = array2table(total_glcm_stats_tmp_table , 'VariableNames', col_names);
    cont = total_glcm_stats_tmp_table.Contrast;
    corr = total_glcm_stats_tmp_table.Correlation;
    ener = total_glcm_stats_tmp_table.Energy;
    homo = total_glcm_stats_tmp_table.Homogeneity;

    cont_med(i) = median(cont);
    corr_med(i) = median(corr);
    ener_med(i) = median(ener);
    homo_med(i) = median(homo); 
end

% converting datatypes for input into table
header_table = clinical_info.ImagingCode;
clinical_info.Properties.RowNames = clinical_info.bcr_patient_uuid;
clinical_info.bcr_patient_uuid = [];
%clinical_info('2cb9975f-2d1e-4c74-867f-8b151c22b246',:) = [];
clinical_info.Properties.RowNames = clinical_info.ImagingCode;
clinical_info.ImagingCode = [];
%clinical_info('RPZ-GBMTCGA-0000002073',:) = [];

clinical_info_arr = table2cell(clinical_info);
clinical_info_table = cell2table(clinical_info_arr.');
clinical_info_table.Properties.RowNames = clinical_info.Properties.VariableNames;
header_table = strrep(header_table, '-','_');
sample_id_names = strrep(sample_id_names, '-','_');
clinical_info_table.Properties.VariableNames = header_table;
clinical_info_tables_match = clinical_info_table(:, sample_id_names);
clinical_info_tables_match.contrast = cont_med;

% reading in clinical features
clinical_info_tables_match_tmp = transposeTable(clinical_info_tables_match);
categorical_headers = {'gender','birth_days_to','race',	'ethnicity','vital_status',	'death_days_to'	,'karnofsky_score'	,'age_at_initial_pathologic_diagnosis'};
clinical_info_tables_match_tmp.Properties.VariableNames = categorical_headers;
clinical_info_tables_match_tmp.birth_days_to = -clinical_info_tables_match_tmp.birth_days_to
%cont_med = table(cont_med, 'VariableNames',{'Contrast'});

% converting numeric values to cell arrays
clinical_info_tables_match_tmp.Contrast = cell(height(clinical_info_tables_match_tmp),1);
cont_med =  num2cell(cont_med);
clinical_info_tables_match_tmp.Correlation = cell(height(clinical_info_tables_match_tmp),1);
corr_med = num2cell(corr_med);
clinical_info_tables_match_tmp.Energy = cell(height(clinical_info_tables_match_tmp),1);
ener_med = num2cell(ener_med);
clinical_info_tables_match_tmp.Homogeneity = cell(height(clinical_info_tables_match_tmp),1);
homo_med = num2cell(homo_med);
for k=1:24
    clinical_info_tables_match_tmp.Contrast(k) = cont_med(k);
    clinical_info_tables_match_tmp.Correlation(k) = corr_med(k);
    clinical_info_tables_match_tmp.Energy(k) = ener_med(k);
    clinical_info_tables_match_tmp.Homogeneity(k) = homo_med(k);
end

% assigning male/female status
for k = 1:24
    if clinical_info_tables_match_tmp.gender{k} == 0
        clinical_info_tables_match_tmp.gender{k} = 'Male';
    else
        clinical_info_tables_match_tmp.gender{k} = 'Female';
    end
end
    
writetable(clinical_info_tables_match_tmp,'ClinicalInfo_GLCM.txt');
G = findgroups(clinical_info_tables_match_tmp.race);
worst_survival = splitapply(@max,clinical_info_tables_match_tmp.death_days_to,G)

%converting table data to matrix
clinical_info_tables_match_tmp.death_days_to = cell2mat(clinical_info_tables_match_tmp.death_days_to);
clinical_info_tables_match_tmp.birth_days_to = cell2mat(clinical_info_tables_match_tmp.birth_days_to);
clinical_info_tables_match_tmp.karnofsky_score = num2cell(clinical_info_tables_match_tmp.karnofsky_score);
clinical_info_tables_match_tmp.age_at_initial_pathologic_diagnosis = num2cell(clinical_info_tables_match_tmp.age_at_initial_pathologic_diagnosis);
clinical_info_tables_match_tmp.Correlation = cell2mat(clinical_info_tables_match_tmp.Correlation);
clinical_info_tables_match_tmp.Contrast = cell2mat(clinical_info_tables_match_tmp.Contrast);
clinical_info_tables_match_tmp.Energy = cell2mat(clinical_info_tables_match_tmp.Energy);
clinical_info_tables_match_tmp.Homogeneity = cell2mat(clinical_info_tables_match_tmp.Homogeneity);
clinical_info_tables_match_tmp.vital_status = categorical(clinical_info_tables_match_tmp.vital_status);
clinical_info_tables_match_tmp.gender = categorical(clinical_info_tables_match_tmp.gender);
clinical_info_tables_match_tmp.race = categorical(clinical_info_tables_match_tmp.race);
clinical_info_tables_match_tmp.ethnicity = categorical(clinical_info_tables_match_tmp.ethnicity);

% Correlation matrix not used
% CorrelationMatrix(table2array(clinical_info_tables_match_tmp));

clinical_info_tables_match_tmp1 = varfun(@cell2mat,clinical_info_tables_match_tmp);
clinical_info_tables_match_tmp1 = clinical_info_tables_match_tmp;
for i = 1:width(clinical_info_tables_match_tmp1)
    if iscell(clinical_info_tables_match_tmp1.(i))
        clinical_info_tables_match_tmp1.(i) = cell2(clinical_info_tables_match_tmp1.(i));
    end
end

% Computing age as a function of GLCM features
lm = fitlm(clinical_info_tables_match_tmp, 'interactions','ResponseVar', 'Contrast', 'PredictorVars',{'vital_status','Correlation'},'CategoricalVar','vital_status'); 
lm1 = fitlm(clinical_info_tables_match_tmp, 'interactions','ResponseVar', 'birth_days_to','PredictorVars',{'Contrast','Homogeneity','Correlation','Energy','vital_status'},'CategoricalVar','vital_status');
an = anova(lm1);

% Computing days survived/OS as a function of GLCM features
clinical_info_noNan = clinical_info_tables_match_tmp;
clinical_info_noNan=clinical_info_noNan(~any(ismissing(clinical_info_noNan),2),:);
lm_new = fitlm(clinical_info_noNan, 'interactions','ResponseVar', 'Homogeneity','PredictorVars',{'Contrast','Energy','Correlation', 'death_days_to'});
am_new = anova(lm_new);


w1 = linspace(min(clinical_info_tables_match_tmp.Contrast),max(clinical_info_tables_match_tmp.Contrast));
w2 = linspace(min(clinical_info_tables_match_tmp.Energy),max(clinical_info_tables_match_tmp.Energy));
w3 = linspace(min(clinical_info_tables_match_tmp.Correlation),max(clinical_info_tables_match_tmp.Correlation));
w4 = linspace(min(clinical_info_tables_match_tmp.death_days_to),max(clinical_info_tables_match_tmp.death_days_to));

%figure()
%gscatter(clinical_info_tables_match_tmp.Contrast, clinical_info_tables_match_tmp.Correlation , clinical_info_tables_match_tmp.Energy, clinical_info_tables_match_tmp.death_days_to, clinical_info_tables_match_tmp.Homogeneity, 'rbg')
%line(w1,feval(fit,w1,),'Color','r','LineWidth',2)
%title('Fitted Regression Lines by GLCM Contrast Features')

w = linspace(min(clinical_info_tables_match_tmp.Correlation),max(clinical_info_tables_match_tmp.Correlation));

figure()
gscatter(clinical_info_tables_match_tmp.Contrast, clinical_info_tables_match_tmp.Correlation , clinical_info_tables_match_tmp.vital_status,'rg','x.o')
line(w,feval(w,lm,'Dead'),'Color','r','LineWidth',2)
line(w,feval(w,lm,'Alive'),'Color','g','LineWidth',2)
title('Fitted Regression Lines by GLCM Contrast Features')


varfun(@median,clinical_info_tables_match_tmp,'InputVariables','birth_days_to','GroupingVariables','race')

tf = find(cellfun('length',regexp(clinical_info_tables_match_tmp.race,'ASIAN')) == 1);
h1 = histogram(clinical_info_tables_match_tmp.birth_days_to(tf));
hold on
tf = find(cellfun('length',regexp(clinical_info_tables_match_tmp.race,'BLACK OR AFRICAN AMERICAN')) == 1);
tf = [tf-1; tf];
h1 = histogram(clinical_info_tables_match_tmp.birth_days_to(tf));
hold on
tf = find(cellfun('length',regexp(clinical_info_tables_match_tmp.race,'WHITE')) == 1);
h1 = histogram(clinical_info_tables_match_tmp.birth_days_to(tf));
xlabel('Survival in Number of Days');
ylabel('Number of Patients');
legend('ASIAN','BLACK OR AFRICAN AMERICAN','WHITE');
title('Survival in Number of days based on Race');
hold off





% for i=1:24
%     idx = find(copy_total_glcm_stats{i}(:,1));
%         for j = 1:size(total_glcm_stats{j}, 2)
%             %idx = union(idx, find(copy_total_glcm_stats{j}(:,i)));
%             total_glcm_stats_tmp = total_glcm_stats(i);
%             total_glcm_stats_tmp(any(cellfun(@(x) any(isnan(x)),total_glcm_stats_tmp(i), 'UniformOutput',false),1),:) = [];
%         end
% end


% total_glcm_stats(idx, :) = [];
% 
% for i=1:size(copy_total_glcm_stats)
%     copy_total_glcm_stats{i}(any(cellfun(@(x) any(isnan(x)),copy_total_glcm_stats{i}),1),:) = [];
% end


xlswrite('GLCM_features.xls',[total_glcm_stats(1), total_glcm_stats(2), total_glcm_stats(3), total_glcm_stats(4), total_glcm_stats(5), total_glcm_stats(6), total_glcm_stats(7), total_glcm_stats(8), total_glcm_stats(9), total_glcm_stats(10), total_glcm_stats(11), total_glcm_stats(12), total_glcm_stats(13),total_glcm_stats(14), total_glcm_stats(15),total_glcm_stats(16),total_glcm_stats(17),total_glcm_stats(18),total_glcm_stats(19),total_glcm_stats(20),total_glcm_stats(21),total_glcm_stats(22),total_glcm_stats(23),total_glcm_stats(24)])


% All code below was used for testing and troubleshooting; not required for the final result


% T1_post = load_nii('/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX T1 POST.nii');
% EEG = load_nii('/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/Output Image Volume-label.nii');
% 
% opt.setvalue.idx = find(T1_post.img);
% opt.setvalue.val = T1_post.img(opt.setvalue.idx);
% opt.command = 'init';
% opt.useinterp = 1;
% option.setviewpoint = [62 48 46];
% view_nii(EEG, opt); % contrast for tumor
% 
% map1 = T1_post;
% idx1 = find(map1.img < -2 | map1.img > +2);
% 
% for i=1:1:71
%     T1_post_temp = T1_post.img(:,:,i);
%     tumor_temp = EEG.img(:,:,i);
%     masked_tumor_Image = bsxfun(@times, T1_post_temp, cast(tumor_temp, 'like', T1_post_temp));
%     glcm = graycomatrix(masked_tumor_Image,'NumLevels',9,'GrayLimits',[]);
%     stats = graycoprops(glcm,'all');
%     glcm_stats(i,:)= struct2array(stats);
% end
% 
% figure, plot([stats.Correlation]);
% title('Texture Correlation as a function of offset');
% xlabel('Horizontal Offset')
% ylabel('Correlation')
% 
% %map2 = T1_post;
% %idx2 = find(map2.img < -10 | map2.img > +10);
% 
% bg = EEG;
% bg.img(idx1) = map1.img(idx1);
% view_nii(bg);
% save_nii(bg.img, '/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/T1POST_label.tif');
% demo_glcms = graycomatrix(bg.img)
% 
% 
% 
% 
% T1_post = load_nii('/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX T1 POST.nii');
% Flair = load_nii('/Users/amitra2/Documents/CompBioRice17/NIfTI_20140122/RPZ-GBMTCGA-0000002071/New Folder/New Folder/AX FLAIR.nii');
% 
% opt.setvalue.idx = find(T1_post.img);
% opt.setvalue.val = T1_post.img(opt.setvalue.idx);
% opt.command = 'init';
% opt.useinterp = 1;
% option.setviewpoint = [62 48 46];
% view_nii(Flair, opt); % showcases edema

%%%%% Features in GLCM 
% Entropy = degree of disorganization; unequality in area of tumor - high entropy would be highly heterogenous - different types of tumor tissue (necrosis, low/high grade)
% more aggressive tumor, more entropy
% Dissimilarity = simialr to entropy - higher order 
% Homeogeneity = opposite of entropy/dissimilarity 
% Variance - similiar to entropy and dissimilarity
% Cluster of shape = similar to variance/heterogeneity
% Correlation = homogeneity
