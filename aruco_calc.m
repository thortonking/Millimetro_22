function [GT_distance,GT_angle,GT_x,GT_y,GT_time]=aruco_calc(foldername)
dirGT=dir([foldername '/*.csv']);
GT=[];GT_distance=[];GT_angle=[];GT_x=[];GT_y=[];GT_time=[];
for fileInd=1:size(dirGT,1)
    GT = importCSV([foldername '/' dirGT(fileInd).name]);
    %if contains(GT{1,10})==1
        GT_timet=datetime(str2double(GT{:,10}), 'convertfrom', 'posixtime', 'Format', 'MM/dd/yy HH:mm:ss.SSS','TimeZone','America/New_York');
    %else
    %    GT_timet=datetime(GT{:,10},'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSSSSS','Format','yyyy-MM-dd HH:mm:ss.SSS');
    %end
    

    GT_distancet=GT{:,1};
    GT_anglet=GT{:,2};
    GT_xt=GT{:,4};
    GT_yt=GT{:,6};
    GT_distance=[GT_distance;GT_distancet];
    GT_angle=[GT_angle;GT_anglet];
    GT_x=[GT_x;GT_xt];
    GT_y=[GT_y;GT_yt];
    GT_time=[GT_time;GT_timet];
end