%研究区平均值，多年时间序列数据准备
clc
clear all;
%% 
%导入投影信息
[a,R]=geotiffread('D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI199901.tif');
info=geotiffinfo('D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI199901.tif');
[m,n]=size(a);

%读入研究年份、时间
%%%修改
begin_year=1999;
end_year=2018;
cd=end_year-begin_year+1;
start_mon=1; 
end_mon=12;
ms=end_mon-start_mon+1;

%% 求NDVI年均值
%逐月读入NDVI
%文件命名规则：年份+月份，对于1-9月，月份前头加个0；eg:1999901
k=1;
for year=begin_year:end_year
    for mon=start_mon:end_mon
    if mon<10
        filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),'0',int2str(mon),'.tif'];
    else
        filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),int2str(mon),'.tif'];
    end
%       filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),int2str(mon),'.tif'];
        data=importdata(filename);
        data_1=reshape(data,m*n,1);
        ndvisum(:,k)=data_1; %
        k=k+1;
    end
end
ndvisum(find(ndvisum<-0.8))=nan; %有效值范围设定

%读入sig-mask
sig_mask=importdata('D:\Grassland\Data_new\bfast\NDVI_prj\0-flag-sigmask.tif');
sig_mask=reshape(sig_mask,m*n,1);
sig_mask=single(sig_mask);%将逻辑值转换为single
sig_mask(find(sig_mask==0))=nan;%对所有0的值赋值为nan

%sig_mask所有有效值,逐月图层数据只有有效值
ndvisum_1=sig_mask.*ndvisum;  %点乘.*

%对非nan列求均值
for col=1:240
    ndvi_col=ndvisum_1(:,col);
    ndvi_mean(1,col)=mean(ndvi_col,'omitnan');
    ndvi_var(1,col)=var(ndvi_col,'omitnan');
    ndvi_std(1,col)=std(ndvi_col,'omitnan');
end
%得到结果ndvi_mean
% a= rmmissing(b);
%% 求PDSI年平均，一张图一个值
%D:\Grassland\Data_new\2_PDSI\PDSI_prj_albers
k=1;
for year=begin_year:end_year
    for mon=start_mon:end_mon
    if mon<10
        filename=['D:\Grassland\Data_new\2_PDSI\PDSI_prj_albers\',int2str(year),'0',int2str(mon),'.tif'];
    else
        filename=['D:\Grassland\Data_new\2_PDSI\PDSI_prj_albers\',int2str(year),int2str(mon),'.tif'];
    end
%       filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),int2str(mon),'.tif'];
        data=importdata(filename);
        data_1=reshape(data,m*n,1);
        pdsisum(:,k)=data_1; %储存所有的频次
        k=k+1;
    end
end
pdsisum(find(pdsisum<-10))=nan;

pdsisum_1=sig_mask.*pdsisum;

for col=1:240
    pdsi_col=pdsisum_1(:,col);
    pdsi_mean(1,col)=mean(pdsi_col,'omitnan');
%     ndvi_var(1,col)=var(ndvi_col,'omitnan');
%     ndvi_std(1,col)=std(ndvi_col,'omitnan');
end
%得到结果pdsi_mean

%% 将NDVI和PDSI年均值写在一个矩阵中
% t1 = datetime(1999,1,1);
% t=dateshift(t1,'start','month',0:239);
% t=double(t);
t=linspace(1999+1/365,2018+(365-30)/365,240);
s_ndvi_pdsi=[t',ndvi_mean',pdsi_mean'];
csvwrite('D:\Grassland\Data_new\Result\s_ndvi_pdsi.csv',s_ndvi_pdsi);
%% 年平均值合成
%导入投影信息
[a,R]=geotiffread('D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI199901.tif');
info=geotiffinfo('D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI199901.tif');
[m,n]=size(a);

%读入研究年份、时间
%%%修改
begin_year=1999;
end_year=2018;
cd=end_year-begin_year+1;
start_mon=5; 
end_mon=9;
ms=end_mon-start_mon+1;
%% NDVI年平均值合成，生长季月份4-10月，求平均值
k=1;
for year=begin_year:end_year
    for mon=start_mon:end_mon
    if mon<10
        filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),'0',int2str(mon),'.tif'];
    else
        filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),int2str(mon),'.tif'];
    end
%       filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),int2str(mon),'.tif'];
        data=importdata(filename);
        data_1=reshape(data,m*n,1);
        ndvisum(:,k)=data_1; %
        k=k+1;
    end
    y=year-1998;
    ndvisum(find(ndvisum<-0.8))=nan; %有效值范围设定
    temp=mean(ndvisum','omitnan');
    ndvi_year(:,y)=temp';
    k=1;
end

%读入sig-mask
sig_mask=importdata('D:\Grassland\Data_new\bfast\NDVI_prj\0-flag-sigmask.tif');
sig_mask=reshape(sig_mask,m*n,1);
sig_mask=single(sig_mask);%将逻辑值转换为single
sig_mask(find(sig_mask==0))=nan;%对所有0的值赋值为nan

%sig_mask所有有效值,
%ndvi_year
ndvi_year=sig_mask.*ndvi_year;  %点乘.* 
ndvi_year_ave=mean(ndvi_year,'omitnan');

for y=1:20
    year=y+1998;
    filename1=['D:\Grassland\Data_new\VI\NDVI_year_ave_5_9\NDVI',int2str(year),'.tif'];
    ndvi_y=ndvi_year(:,y);
    ndvi_y_2=reshape(ndvi_y,[m,n]);
    geotiffwrite(filename1,ndvi_y_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
end
%ndvi_year,研究区20年所有数据
%ndvi_year_ave，是20年的每年平均

%% PDSI平均值
%逐月读入生长季，并计算年均值
k=1;
for year=begin_year:end_year
    for mon=start_mon:end_mon
    if mon<10
        filename=['D:\Grassland\Data_new\2_PDSI\PDSI_prj_albers\',int2str(year),'0',int2str(mon),'.tif'];
    else
        filename=['D:\Grassland\Data_new\2_PDSI\PDSI_prj_albers\',int2str(year),int2str(mon),'.tif'];
    end
%       filename=['D:\Grassland\Data_new\VI\NDVI_rectangle_prj\NDVI',int2str(year),int2str(mon),'.tif'];
        data=importdata(filename);
        data_1=reshape(data,m*n,1);
        pdsisum(:,k)=data_1; %
        k=k+1;
    end
    y=year-1998;
    pdsisum(find(pdsisum<-10))=nan; %有效值范围设定
    temp=mean(pdsisum','omitnan');
    pdsi_year(:,y)=temp';
    k=1;
end
pdsi_year=sig_mask.*pdsi_year;  %点乘.* 
% pdsi_year=double(pdsi_year);
pdsi_year_ave=mean(pdsi_year,'omitnan');
pdsi_year_ave=double(pdsi_year_ave);

for y=1:20
    year=y+1998;
    filename1=['D:\Grassland\Data_new\2_PDSI\Year_ave_5_9\PDSI',int2str(year),'.tif'];
    pdsi_y=pdsi_year(:,y);
    pdsi_y_2=reshape(pdsi_y,[m,n]);
    geotiffwrite(filename1,pdsi_y_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
end
%pdsi_year,研究区20年所有数据
%pdsi_year_ave，是20年的每年平均,1*20
%% 写出数据 
s_ndvi_pdsi=[ndvi_year_ave',pdsi_year_ave'];
csvwrite('D:\Grassland\Data_new\Result\s_ndvi_pdsi_year_5_9.csv',s_ndvi_pdsi);
%% 求相关系数
num = size(ndvi_year,1);
for i = 1:num
    pdsidata=pdsi_year(i,:);
    ndvidata=ndvi_year(i,:);
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
filename1=['D:\Grassland\Data_new\Result\Correlation\pdsi_NDVI_corre-2.tif'];
geotiffwrite(filename1,c_pdsi_ndvi_r,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
filename2=['D:\Grassland\Data_new\Result\Correlation\pdsi_NDVI_corre_p-2.tif'];
geotiffwrite(filename2,c_pdsi_ndvi_p_r,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
%% 读取7种变化类型
corr=importdata('D:\Grassland\Data_new\Result\Correlation\pdsi_NDVI_corre_sig.tif');

for i=0:8
    if i==3 | i==4
        continue
    else
        filename=['D:\Grassland\Data_new\bfast\NDVI_flag_new\n_NDVI_Type',num2str(i),'.tif']
        data=importdata(filename);
        Type(:,i+1)=reshape(data,m*n,1);
    end
end

Type=single(Type);
Type(find(Type==15))=nan;
Type(find(Type==3))=nan;