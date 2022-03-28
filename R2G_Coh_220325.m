%% RGB_to_gray 2021.12.18  2022.03.12

clc;
clear;
close all; 

% Image_raw=imread('F:\matlab_project\2021.12.16_Coherence_PVSK_PhC\ForMatlab\sunyu-10um-50x.tif');
fr_aera_coherence=[]; area_ave_test =[];  fr_aera_P11 =[]; P11=[];
um=0.05613;% 20um = 1.31cm  5440��3648pix 20��13.41cm  % ���ֶ���ȡ���ݶ�Ӧ�����ص�
CheckPosition = 0; % ��һ������ȷ����������λ�� 1ִ��ȷ��λ�� 0 ȷ��������
Show_Pic = 0; % �Ƿ���ʾͼ 1��ʾ 0 ����ʾ
Save_Pic = 0; % �Ƿ��ͼ 1�� 0 ����
Save_Data = 1; % �Ƿ������ 1�� 0 ����

           Boundary = [121; 556; 506; 146];  % Whole  [183; 762; 749; 175];
        %   Boundary = [121; 556; 173; 146];  % Part  [183; 762; 228; 180];
    boundary_x1=Boundary(1,:); boundary_y1=Boundary(2,:);  % Ŀ������pix����
    boundary_x2=Boundary(3,:); boundary_y2=Boundary(4,:); 
    
file=dir('F:\matlab_project\Year2022-Q1\2022.03.25_Coherence_PVSK_PhC\S28-Process\*.tif');
n_max = length(file);
  
for n =1:n_max
  strposition = strfind(file(n).name, 'S28-');
  if ~isempty (strposition)
      Image_raw = imread(file(n).name);  % �������ļ�����
      Fr_Gray=rgb2gray(Image_raw);

%% Fig1_ Raw_Pic
                if Save_Pic == 1         
                   figure('Name','Fig1_Raw_Pic','NumberTitle','off');
                end
          if (1)
          set(gcf,'Colormap',bone);
          set(gcf,'Position',[0 600 400 300]);  % Figure_1 [������Ļ��ߵ����� ������Ļ�±ߵ����� �������� ��������]
          set(gca,'Position',[.18 .17 .6 .67]);
          colorbar('position',[0.82 0.17 0.03 0.73]); 
          end
            
          pcolor(Fr_Gray); shading interp%α��ɫͼ
              xlabel('x-axis (pixel)','FontName','Times newman','FontSize',12);
              ylabel('y-axis (pixel)','FontName','Times newman','FontSize',12);%����xy���ǩ���ݺ�����
          %imshow(Gray);
          title(file(n).name,'�Ҷ�ͼ��')
          hold on
%% Fig1_2_Selected_Aera
    boundary_x=[boundary_x1,boundary_x1,boundary_x2,boundary_x2];  % ��ע�ĺ��λ��
    boundary_y=[boundary_y1,boundary_y2,boundary_y2,boundary_y1];

    plot(boundary_x, boundary_y,'+r','linewidth', 2, 'markersize', 4);

   if Save_Pic == 1;
    print([num2str(file(n).name),'-Fig1_Raw_Pic','.tif'] ,'-dtiffn','-r150');   %,'-2'Ϊ��׺
   end
    end  
   
%% Fig2_Selected_Aera
    if Show_Pic == 1; 
    figure('Name','Fig2_Selected_Aera','NumberTitle','off');
    
    %Gray_target = Gray(boundary_y1:boundary_y2,boundary_x1:boundary_x2,:);
    Gray_target = Fr_Gray(boundary_y2:boundary_y1,boundary_x1:boundary_x2);
    pcolor(Gray_target); shading interp%α��ɫͼ
    %imshow(Gray_target);
              xlabel('x-axis (pixel)','FontName','Times newman','FontSize',12);
              ylabel('y-axis (pixel)','FontName','Times newman','FontSize',12);%����xy���ǩ���ݺ�����
       title(file(n).name,'Ŀ������')
    
          if(1)
            set(gcf,'Colormap',bone);
            set(gcf,'Position',[0 100 400 300]); % Figure_2
            set(gca,'Position',[.18 .17 .6 .67]);
          end
          
    if Save_Pic == 1;  
    print([num2str(file(n).name),'-Fig2_Selected_Aera','.tif'] ,'-dtiffn','-r150');   %,'-2'Ϊ��׺
    end      
    end
    
    
if CheckPosition == 0;    
        range_x=abs(boundary_x1-boundary_x2)+1;
        range_y=abs(boundary_y1-boundary_y2)+1;
        area_cau=zeros(range_y,range_x);
        i_x=min(boundary_x):max(boundary_x);
        i_y=min(boundary_y):max(boundary_y);
        j_x=1:range_x;
        j_y=1:range_y;
        j_x_position = j_x*um;
        j_y_position = j_y*um;
        area_cau(j_y,j_x)=Fr_Gray(i_y,i_x);
        area_cau_y=area_cau';
        % area_sum_x=zeros(1,range_x);
        area_ave_x=sum(area_cau)./range_y;    % x����ƽ��ǿ�ȷֲ�
        area_ave_y=sum(area_cau_y)./range_x;  % y����ƽ��ǿ�ȷֲ�
    
%% Fig3_Slice_Selected
    if Show_Pic == 1; 
    figure('Name','Fig3_Slice_Selected','NumberTitle','off');
    end
    
    % figure 4-1    
    subplot(2,1,1)
        plot(j_x_position,area_ave_x);
              xlabel('x-axis (pixel)','FontName','Times newman','FontSize',12);
              ylabel('Intensity (a.u.)','FontName','Times newman','FontSize',12);%����xy���ǩ���ݺ�����
    % figure 4-2
    subplot(2,1,2)
        plot(j_y_position,area_ave_y);
        % area_ave_test = area_ave_x;
        % j_test = j_x;% x �����ź�
              xlabel('x-axis (pixel)','FontName','Times newman','FontSize',12);
              ylabel('Intensity (a.u.)','FontName','Times newman','FontSize',12);%����xy���ǩ���ݺ�����
        set(gcf,'Position',[450 100 600 800]);
        
    if Save_Pic == 1;     
    print([num2str(file(n).name),'-Fig3_Slice_Selected','.tif'] ,'-dtiffn','-r150');   %,'-2'Ϊ��׺
    end              

        area_ave_test = area_ave_y;
        j_test = j_y;% y �����ź�
        %%
        Fs = round(1/um);
        L=length(area_ave_test);%�źų���
        t= (0:L-1)*um;

    
%% Fig4_Mode_Fitting_Selected
        if Show_Pic == 1; 
           figure('Name','Fig4_Mode_Fitting','NumberTitle','off');
        end
    % figure 4-1
    subplot(2,1,1) 
        ax1=plot(j_test*um,area_ave_test);
        set(gca,'LineWidth',1,'FontSize',12);
              xlabel('x-axis (pixel)','FontName','Times newman','FontSize',12);
              ylabel('Intensity (a.u.)','FontName','Times newman','FontSize',12);%����xy���ǩ���ݺ�����
              
        ax1.LineWidth=2;
        Y = fft(area_ave_test);       % ����Ҷ�任
           % Y = fft(area_ave_test-mean(area_ave_test));
        P22 = abs(Y/L);
           % P11 = P22(1:L/2+1);
        P11 = P22(1:round(L/2));
        P11(2:end-1) = 2*P11(2:end-1);
        f = Fs*(0:(L/2))/L;     % Ƶ��

        
    % figure 4-2
    if Show_Pic == 1; 
    subplot(2,1,2);
        ax2=plot(f,P11);
        ax2.LineWidth=2;
       
        title(file(n).name,'Single-Sided Amplitude Spectrum of Test(t)')
        xlabel('f (um^-^1)','FontName','Times newman','FontSize',12)
        ylabel('|P1(f)|','FontName','Times newman','FontSize',12)
        set(gca,'LineWidth',1,'FontSize',12);
        set(gcf,'Position',[1100 100 600 800]);
    end
    
    if Save_Pic == 1;    
    print([num2str(file(n).name),'-Fig4_Mode_Fitting_Selected','.tif'] ,'-dtiffn','-r150');   %,'-2'Ϊ��׺
    end
  
    %% Fig5_Mode_Fitting_Filter
   if (0)
    figure('Name','Fig5_Mode_Fitting_Filter','NumberTitle','off');
    % figure 5-1   S26---24 45
%             P11(:,1:9) = 0;
%             P11(:,19:end) = 0;
            P33 = abs(ifft(P11));
            ax3 = plot(f,P33);
   end
            
            
  fr_aera_coherence=[(j_test*um)', fr_aera_coherence,area_ave_test'];
%   fr_aera_coherence_data = 
%   save('S4-Coherence-data.txt','fr_aera_coherence','-ascii');
  fr_aera_P11      =[f',           fr_aera_P11      ,P11'          ];
%   save('S4-Coherence-data-fft.txt','fr_aera_P11','-ascii');
    end
end


   if Save_Data == 1;
      %���ԭʼ����
      fr_aera_coherence_data = fr_aera_coherence(:,n:end);  
      fr_aera_coherence_data = [(0:1:n);fr_aera_coherence_data];
         save('S28-5um-W_Coh_data.txt','fr_aera_coherence_data','-ascii');
%       fr_aera_coherence_axis = [0;t'];   
%          save('S93-6-2um_Coh_data_axis1.txt','fr_aera_coherence_axis','-ascii');
      %���FFT����   
      fr_aera_P11_data       = fr_aera_P11(:,n:end);
      fr_aera_P11_data       = [(0:1:n);fr_aera_P11_data];
         save('S28-5um-W_Coh_data_fft.txt','fr_aera_P11_data','-ascii');
%       fr_aera_P11_axis       = [0;f'];   
%          save('S93-6-2um_Coh_data_fft_axis1.txt','fr_aera_P11_axis','-ascii');
   end
