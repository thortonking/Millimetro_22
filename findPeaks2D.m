function [indices,pVal]=findPeaks2D(Pmusic,L)

 % Get AoAs that have peaks
 %{
    binary_peaks_pmusic = imregionalmax(Pmusic);
   
    aoa_indices = linspace(1, size(binary_peaks_pmusic, 1), size(binary_peaks_pmusic, 1));
    ToF_indices = linspace(1, size(binary_peaks_pmusic, 2), size(binary_peaks_pmusic, 2));
    
    aoa_peaks_binary_vector = any(binary_peaks_pmusic, 2);
    aoa_peak_indices = aoa_indices(aoa_peaks_binary_vector);
    aoa_peak_indices=aoa_peak_indices(aoa_peak_indices>10);
    aoa_peak_indices=aoa_peak_indices(aoa_peak_indices<size(theta,2)-10);
    
    tof_peak_indices=zeros(0,length(aoa_peak_indices));
    estimated_aoas = theta((aoa_peak_indices));
    % Get ToFs that have peaks
    %time_peak_indices = zeros(length(aoa_peak_indices), size(Pmusic,2));
    % AoA loop (only looping over peaks in AoA found above)
    ii=1;
    indiceSize=length(aoa_peak_indices);
    while ii <=indiceSize
        aoa_index = aoa_peak_indices(ii);
        binary_tof_peaks_vector = binary_peaks_pmusic(aoa_index, :);
        %matching_tofs = tau(binary_tof_peaks_vector);
        ind=ToF_indices(binary_tof_peaks_vector);
        if ind(1)>10 && ind(1)<size(Pmusic,2)-10
            tof_peak_indices(ii)=ind(1);
             ii=ii+1;
        else
            estimated_aoas(ii)=[];
            aoa_peak_indices(ii)=[];
            indiceSize=length(aoa_peak_indices);
        end
        % Pad ToF rows with -1s to have non-jagged matrix
        %negative_ones_for_padding = -1 * ones(1, length(tau) - length(matching_tofs));
        %time_peak_indices(ii, :) = horzcat(matching_tofs, negative_ones_for_padding);
       
    end
%}


BW = imregionalmax(Pmusic);
[angle_idx,delay_idx, delta_idx] = ind2sub(size(BW), find(BW));
indices=[angle_idx,delay_idx];

%indices=indices(5<indices(:,2),:);
%indices=indices(indices(:,2)<size(Pmusic,2)-5,:);
removeIndices=find(5>indices(:,2) | indices(:,2)>size(Pmusic,2)-5 | 5>indices(:,1) | indices(:,1)>size(Pmusic,1)-5);
[sortedX,sortingIndices] = sort(Pmusic(BW),'descend');
sortedX=sortedX(~ismember(sortingIndices,removeIndices));
sortingIndices=sortingIndices(~ismember(sortingIndices,removeIndices));

if L<length(sortedX)
    maxValues = sortedX(1:L);
    maxValueIndices = sortingIndices(1:L);
else
    [sortedX,sortingIndices] = sort(Pmusic(BW),'descend');
    if L<length(sortedX)
        maxValues = sortedX(1:L);
        maxValueIndices = sortingIndices(1:L);
    else
    maxValues = sortedX;
    maxValueIndices = sortingIndices;
    end
end
indices=indices(maxValueIndices,:);
pVal=sortedX(1:size(indices,1),:);
%indices=indices(10<indices(:,1),:);
%indices=indices(indices(:,1)<size(Pmusic,1)-10,:);


%sortrows(indices,2);
%for i=1:length(angle_idx)
%    if angle_idx(i)

