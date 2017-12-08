function [] = ExtractFeatures()

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
subdirinfo = cell(length(dirinfo));
for K = 1 : length(dirinfo)
  thisdir = dirinfo(K).name;
  subdirinfo{K} = dir(fullfile(thisdir, '*TI post.nii'));
  subdirinfo2{K} = dir(fullfile(thisdir, '*-label.nii'));
end