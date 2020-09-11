%folder_name = 'D:\VMPO Data\VMPO photometry\analysis\averages\';

start_folder = pwd;
selpath = uigetdir;
cd(selpath);

pre = '_wake_avgs';
during = '_nrem_avgs';
post = '_rem_avgs';

[pre_array, pre_means] = arrayOfMeans(pre);
[during_array, during_means] = arrayOfMeans(during);
[post_array, post_means] = arrayOfMeans(post);

wholeDataSet = [pre_means,during_means,post_means];

cd(start_folder)