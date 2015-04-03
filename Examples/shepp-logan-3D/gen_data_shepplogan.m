%Script to create test X-ray tomo data
clc; clear; close all;

proj_cols = 256;%2560
proj_num = 256;%2048;
proj_rows = 16;
noise_sigma = 0.5;
dose = 5000;
time_samples = 1;
write_file = 1;

angles = 0:180/proj_num:180*time_samples-180/proj_num;
P=phantom(proj_cols);
P=P./1000; %Values are between 0 - .001 inverse distance

Obj = repmat(P,[1,1,proj_rows]);
data = zeros(proj_cols,proj_rows,proj_num*time_samples);
for i = 1:proj_rows
    R = radon(squeeze(Obj(:,:,i)),angles);
    R = R';
    [m n] = size(R);
    Ax = R(:,int16(n/2)-proj_cols/2:int16(n/2)+proj_cols/2-1);
    data(:,i,:) = (dose.*exp(-Ax))';
end

noisy_data = data + noise_sigma.*sqrt(data).*randn(size(data));
projection = log(dose./noisy_data);

if(write_file)
    proj_time = 1:proj_num*time_samples;
    recon_time =1:proj_num:proj_num*time_samples+1;
    recon_time(end)=proj_num*time_samples;
    
    fid = fopen('projections.bin', 'wb');
    fwrite(fid, projection, 'float');
    fclose(fid);
    
    fid = fopen('weights.bin', 'wb');
    fwrite(fid, noisy_data, 'float');
    fclose(fid);
    
    fid = fopen('proj_angles.bin', 'wb');
    fwrite(fid, angles*pi/180, 'float');
    fclose(fid);
    
    fid = fopen('proj_times.bin', 'wb');
    fwrite(fid, proj_time, 'float');
    fclose(fid);
    
    fid = fopen('recon_times.bin', 'wb');
    fwrite(fid, recon_time, 'float');
    fclose(fid);
end
