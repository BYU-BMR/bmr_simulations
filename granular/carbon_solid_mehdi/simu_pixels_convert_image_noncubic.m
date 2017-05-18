function simu_pixels_convert_image_noncubic()

% ***this function convert simulation output to a serious 2d images, avaliable for uniform and non-uniform length of boxes.
clc;
clear;

%NX=66;
%NY=33;
%NZ=26;

if ~exist('2dslices', 'dir')
  mkdir('2dslices');
end

data=importdata('pic_property_calculation_large_length_vol_change.txt');

%Array_size=size(data);
%row_size=Array_size(1);
%length=ceil(power(row_size,1/3));

data2=importdata('noncubic_box_length.txt');

NX=data2(1,1);
NY=data2(2,1);
NZ=data2(3,1);

site=0;
for k=1:NX
   
    for j=1:NY
        
        
       for i=1:NZ
           site=site+1;
       %   data2(site)=data((k-1)+1+NZ*(i-1)+NY*NZ*(j-1));
           data2(site)=data(site);
       end
       
    end
    
end




site=0;
for i=1:NX
x1=sprintf('%2.3i',i); 
    
 FILENAME3 = ['image_convert_',[x1,'.txt']]; 
  fid_61=fopen(strcat('2dslices/',FILENAME3),'wt');
  
     for j = 1: NY
      for k = 1: NZ
         % site = (i-1)*NY*NZ + (j-1)*NZ + k; % set up sites from 1 at first slice, 2 at second slice ....
          site=site+1;                
  
            fprintf(fid_61,'%g\t',data2(site));            % input domain identities at each node
            
      end
      
       fprintf(fid_61,'\n');
     end 
    
     
end






for z=1:NX
    
x1=sprintf('%2.3i', z); 
 FILENAME1 = ['image_convert_',[x1,'.txt']]; 
 
 FILENAME2 = ['pixels_convert_image_simulation',[x1,'.tif']]; 

new_image= dlmread(strcat('2dslices/',FILENAME1));%add the text file that has the numbers

         thresh = multithresh(new_image,2);
         seg_I = imquantize(new_image, thresh);
         RGB = label2rgb(seg_I);       
         imwrite(RGB,strcat('2dslices/',FILENAME2));
 
end
 
         
         
fprintf(1,'%s \n', 'converting pixels to images successfully ');
         
end