clc
clear all;
%导入投影信息
[a,R]=geotiffread('D:\Grassland\Data\PDSI\201812.tif');
info=geotiffinfo('D:\Grassland\Data\PDSI\201812.tif');
[m,n]=size(a);

begin_year=1999;
end_year=2018;
cd=end_year-begin_year+1;
start_mon=4;
end_mon=9;
ms=end_mon-start_mon+1;

%pdsisum=zeros(m,n,cd*ms);
k=1;
for year=begin_year:end_year
    for mon=start_mon:end_mon
        filename=['D:\Grassland\Data\PDSI\',int2str(year),int2str(mon),'.tif'];
        data=importdata(filename);
        data_1=data(:);
        pdsisum(:,k)=data_1;
        k=k+1;
    end
end

%ndvisum=zeros(m,n,cd*ms);
k=1;
for year=begin_year:end_year
    for mon=start_mon:end_mon
        filename=['D:\Grassland\Data\VI\LAI\LAI',int2str(year),'0',int2str(mon),'.tif'];
        data=importdata(filename);
        %data=im2double(data);
        data=double(data);
        data_2=data(:);
        ndvisum(:,k)=data_2; %储存所有的频次
        k=k+1;
    end
end
ndvisum=ndvisum*0.004-0.08;

%建立空矩阵
%c_pdsi_ndvi=zeros(m,n)+nan;
%c_pdsi_ndvi_p=zeros(m,n)+nan;
num = size(ndvisum,1);
for i = 1:num
    pdsidata=pdsisum(i,:);
    ndvidata=ndvisum(i,:);
    %if min(pdsidata)>-10 & min(ndvisum)>=0 %有效性判断
    %相关
    [r1,p1]=corr(pdsidata',ndvidata'); %降水和ndvi相关
    c_pdsi_ndvi(i,1)=r1;
    c_pdsi_ndvi_p(i,1)=p1;
    % else
%end
    if rem(i,10000)==0
        disp(num2str(i/num))
    else
    end
end


c_pdsi_ndvi_r = reshape(c_pdsi_ndvi,[m,n]);
c_pdsi_ndvi_p_r = reshape(c_pdsi_ndvi_p,[m,n]);

%相关系数输出
filename1=['D:\Grassland\Results\correlations\pdsi_LAI_corre.tif'];
geotiffwrite(filename1,c_pdsi_ndvi_r,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

filename2=['D:\Grassland\Results\correlations\pdsi_LAI_corre_p.tif'];
geotiffwrite(filename2,c_pdsi_ndvi_p_r,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

% function [t1] = progress(num,all_num,t)
% if num/all_num > t
%     disp(t);
%     t1 = t + 0.01;
% else
% end
% end

