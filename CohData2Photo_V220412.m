%                       _oo0oo_
%                      o8888888o
%                      88" . "88
%                      (| -_- |)
%                      0\  =  /0
%                    ___/`---'\___
%                  .' \\|     |// '.
%                 / \\|||  :  |||// \
%                / _||||| -:- |||||- \
%               |   | \\\  -  /// |   |
%               | \_|  ''\---/''  |_/ |
%               \  .-\__  '-'  ___/-. /
%             ___'. .'  /--.--\  `. .'___
%          ."" '<  `.___\_<|>_/___.' >' "".
%         | | :  `- \`.;`\ _ /`;.`/ - ` : | |
%         \  \ `_.   \_ __\ /__ _/   .-` /  /
%     =====`-.____`.___ \_____/___.-`___.-'=====
%                       `=---='
%
%
%     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%               ���汣��         ����BUG

% 20220412��
% ��CohData������ݽ��л�ͼ
% ������Զ���ȡ���ֵ�������ݼ�ֵ������ʱ��ĺ���

clc;
clear;
close all; 
file=dir('F:\matlab_project\Year2022-Q2\2022.04.04_Coh_PVSK_PhC_PowDep\Result_FFT\*.txt');  %�����ļ��ĵ�ַ
n_max = length(file);
Freq_FFT_Zero = 0.5;   %  �㼶Ƶ�ʵķ�Χ������0~0.5�� (1ɾ��/0����)
    ZeroExcavate = 1;
    CurveFitting = 1;  %  �˴��ж��Ƿ�ɾ���㼶Ƶ�� (1ɾ��/0����)
    
   %% Step0_ Add Time ���ʱ��
        fr_time = importdata('TimeScale-1-10um-31��4.txt');                  %ת������No��λ��Position
        % fr_time_min = min(min(fr_time));
        fr_time = (fr_time - (ones(size(fr_time)).* min(min(fr_time)))).*3.33;  %DelayLine��λ�ú�ʱ���ת��
        fr_time = [0;reshape(fr_time,[],1)];                                    %�ع�DelayTime

   %% Step1_��ʼ��������
for n =1:n_max
     strposition = strfind(file(n).name, 'S09-09003-P_Coh_FFTdata');  % �ؼ���ѡȡ
  if ~isempty (strposition)
     fr_0_Ex = importdata(file(n).name);  % �������ļ�����
    
          [size_x,size_y]=size(fr_0_Ex);
          fr_time_1 = ( fr_time(1:size_y,:))';
          fr_whole_rescaled = [fr_time_1;fr_0_Ex(2:end,:)];   % ����ʱ������ 
          fr_Ex_sorted = sortrows(fr_whole_rescaled',1);      % ����ʱ���������� (sortrows)
          fr_0_Ex = fr_Ex_sorted;
          
     fr_0_x = fr_0_Ex(1,2:end);
     fr_0_y = fr_0_Ex(2:end,1);
     fr_0_Ex = fr_0_Ex(2:end,2:end);
     
         Size_fr_0.x = length(fr_0_x);
         Size_fr_0.y = length(fr_0_y);
         xx_res =Size_fr_0.x;             %x������ֱ���
         yy_res =Size_fr_0.y;             %y������ֱ���
         n_max_res = xx_res * yy_res; %����ֵ���� 
         
    %%ת��
        fr_0_Ex_one = ones(size(fr_0_Ex));   %������Ex��ͬ�ṹ��ȫ1����
        fr_1_x = fr_0_Ex_one' .* fr_0_x';          %fr_0_Ex_one' .* fr_0_x';
        fr_1_y = fr_0_y' .* (fr_0_Ex_one)';
        fr_1_Ex = fr_0_Ex';
        
    %%ɾ���㼶����    
      if ZeroExcavate==1;  
        % fr_0_x ��С��N������
            [~, Position_Freq_FFT_Zero]=min(abs(fr_0_x - Freq_FFT_Zero));     % Ѱ�����Freq_FFT_Zeroλ�õ�λ��
            %(1:6)��ɾ�� / ɾ��0.5Ƶ��ǰ������
            fr_1_x = fr_1_x((Position_Freq_FFT_Zero+1):end,:);
            fr_1_y = fr_1_y((Position_Freq_FFT_Zero+1):end,:);
            fr_1_Ex = fr_1_Ex((Position_Freq_FFT_Zero+1):end,:);
      end
      
      %%Ѱ��Matrix�����ֵ��λ��
        fr_1_Ex_Max = max(max(fr_1_Ex));
        [fr_1_Ex_x fr_1_Ex_y] = find (fr_1_Ex_Max == fr_1_Ex);
        
      %%����Matrix�����ֵȷ����Ӧ��Ƶ��λ�ã���ȡ��Ӧ��Ƶ�����ߣ��������
      if CurveFitting==1;    
            fr_Curve.Ex = fr_1_Ex((fr_1_Ex_x-2):(fr_1_Ex_x+2),:); % curve����
            fr_Curve.x = fr_1_x((fr_1_Ex_x-2):(fr_1_Ex_x+2),1);   % curve��Ӧ��Ƶ����Ϣ
            fr_Curve.y = fr_1_y(1,:);                             % curve��Ӧ��ʱ����Ϣ
            fr_Curve_total = [fr_Curve.y;fr_Curve.Ex];

%             FileName_a = file(n).name;
%             FileName_b = '-Fitting.txt';                          % �ļ�������Ӻ�׺�����������������⸲��
%             FileName_c = strcat(FileName_a,FileName_b);
%              save(FileName_c,'fr_Curve_total','-ascii');             % save(file(n).name,'fr_Ex_sorted','-ascii');[Old Version]
            
           if(1) % ���
                   [Line Row] = size (fr_Curve_total);     % n_max = Line;  % ��ȡ��������
                   Autoguess = 1; % ʹ���Զ���ֵ
                   % A= 1; Omega=2; Theta=0; C=0;
                   Para_x_c= 150; Para_w=500; Para_y_0 =0; Para_Aera=10000;  %
                   Parameter_Sum = [];
                   Fitting_Curve_Sum = [];

                if (1)
                    syms t; %�趨�Ա���
                    n_fitting_max = Line;
                    for n_fitting =2:n_fitting_max;        % ��һ����ʱ������
                         x = fr_Curve_total (1,:) ;    % ��ȡʱ������ 
                         y = fr_Curve_total (n_fitting,:);   % ��ȡ��Ҫ��ϵ�����

                   %% �������� ����������Ϊ�Զ���ĺ������Ա���t����������    % ������ߣ�A*sin(omega.*t + theta) + C
                         % Gauss
                          ft=fittype('y_0 + A/(w.*(sqrt(pi./2)))* exp(-2*((t-x_c)^2)./((w)^2))','independent','t','coefficients',{'y_0','w','x_c','A'});  
                         % Lorentz
                         % ft=fittype('y_0 + 2A/pi*(w/(4*(t-x_c)^2 + w^2))','independent','t','coefficients',{'y_0','w','x_c','A'});  % Lorentz

                        if(1) % ������߳�ֵ�趨
                                  if Autoguess == 1;    
                                      Para.y_0 = 0 ;  
                                      Para.w =  Para_w;
                                      Para.x_c =  fr_1_y(1,fr_1_Ex_y); % Para_x_c;   
                                      Para.A = Para_Aera;                
                                  else
                                      Para.A = A;             %��Сֵ�����ֵ����һ����Ϊ��ֵ
                                      Para.Omega = Omega;
                                      Para.Theta = Theta;
                                      Para.C = C;  
                                  end
                            opts = fitoptions(ft); 
                            opts.StartPoint=[Para.y_0 Para.w Para.x_c Para.A]; % ��ʼ��������
                        end    

                        % ��ʼ���
                        fittedfun=fit(x',y',ft,opts);%�����ϵõ��ĺ���  
                            xi = [0:(max(x)./1000):max(x)];   % 1000����ϵ�
                            yi = fittedfun(xi');              % ��Ͻ��

                         Parameter_Sum_i = [fittedfun.x_c, fittedfun.w, fittedfun.y_0, fittedfun.A];  % ��ϲ���
                         Fitting_Sum_i = [yi']; %�������
                         % �ֲ��洢�������
                         Parameter_Sum = [Parameter_Sum; Parameter_Sum_i];
                         Fitting_Curve_Sum = [Fitting_Curve_Sum; Fitting_Sum_i];
                    end
            %%�����������
                    Parameter_Sum = [fr_Curve.x,Parameter_Sum];
                    Fitting_Curve_Sum = [xi;Fitting_Curve_Sum]';
                    fr_Curve_total=fr_Curve_total';
                       % ��ͬ�ߴ������������ϵ�ͬһ��������  (��� fr_Curve_total��
                           % [Size_Para_Sum.x Size_Para_Sum.y] = size(Parameter_Sum);        
                       [Size_Fit_Cure_Sum.x Size_Fit_Cure_Sum.y] = size(Fitting_Curve_Sum);  
                       [Size_fr_Cure_Tot.x Size_fr_Cure_Tot.y] = size(fr_Curve_total);
                       fr_Cure_Output =zeros([max(Size_Fit_Cure_Sum.x,Size_fr_Cure_Tot.x) (Size_Fit_Cure_Sum.y+Size_fr_Cure_Tot.y)]);
                       fr_Cure_Output(1:Size_fr_Cure_Tot.x,1:Size_fr_Cure_Tot.y) = fr_Curve_total ;  
                       fr_Cure_Output(1:Size_Fit_Cure_Sum.x,(Size_fr_Cure_Tot.y+1):end) = Fitting_Curve_Sum ; 
            %%��������   
                FileName_Cure_Output.a = file(n).name;
                FileName_Cure_Output.b = '-Cure_Output.txt';                         
                FileName_Cure_Output.c = strcat(FileName_Cure_Output.a,FileName_Cure_Output.b);
                save(FileName_Cure_Output.c,'fr_Cure_Output','-ascii');   
                
                FileName_Para_Sum.a = file(n).name;
                FileName_Para_Sum.b = '-Fitting_Parameter_Sum.txt';                         
                FileName_Para_Sum.c = strcat(FileName_Para_Sum.a,FileName_Para_Sum.b);                
                save(FileName_Para_Sum.c,'Parameter_Sum','-ascii'); 
                    
                end
           end
      end    
                  
      %%������߹�һ��           
      [fr_1_Ex_x, fr_1_Ex_y]= size(fr_1_Ex);
%          n_max_y = length(fr_0_y);
%          xx_res =n_max_x;             %x������ֱ���
%          yy_res =n_max_y;             %y������ֱ���
         n_max_res_2 = fr_1_Ex_x * fr_1_Ex_y;
            
            fr_2_x =reshape(fr_1_x,n_max_res_2,1);
            fr_2_y =reshape(fr_1_y,n_max_res_2,1);
            fr_2_Ex =reshape(fr_1_Ex,n_max_res_2,1);
    %%�ϲ�
        fr_total = [[fr_2_x';fr_2_y';fr_2_Ex']'];
        
   %%��ͼ        
      figure
      %yy_region = fr_x_0 * fr_Ex_one fr_x_0;
      yy_region = fr_1_x;
      xx_region = fr_1_y;
      zz_region = fr_1_Ex;
      pcolor(yy_region,xx_region,zz_region);shading interp%α��ɫͼ
      set(gcf,'Colormap',hot)
       % xlim([440 570])  ylim([-10 1000]) 

      xlabel('x (��m)','FontName','Arial','FontSize',11);
      ylabel('Delay time (fs)','FontName','Arial','FontSize',11);%����xy���ǩ���ݺ�����
      set(gca, 'Fontname', 'Arial', 'Fontsize', 11);%����xy����������ͺ��ֺŴ�С��
      set(gca,'Position',[.18 .17 .6 .75]);%���������xy����ͼƬ��ռ�ı�����������Ҫ�Լ�΢����
                 % zz_region_1 = reshape(zz_region,1,n_max_res_2);
                 % Zmax= max(zz_region_1);
                 % set(gca,'CLim',[0 Zmax]); %����color bar ��ʾ��Χ
      set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);%����ǽ��߿��Ϊ2
      %set(gca,'Ytick',[-30:30:60]);%����xy��Ŀ̶���ʾ���� %[Grating600_Center530nm]
      %set(gca,'Xtick',[-60:30:60]);
%               xlim([0 2])
%               set(gca,'Xtick',[0:0.5:2]);
      set(gcf,'Position',[100 100 400 320]);%��������û�ͼ�Ĵ�С������Ҫ��word���ٵ�����С���Ҹ��Ĳ�����ͼ�Ĵ�С��7cm
      colorbar('position',[0.82 0.17 0.035 0.75]); %���ò�ɫ����λ�úʹ�С
      %caxis([0,10]); %����ɫ�ȷ�Χ����ʾ��ʵ����ɫͼ�ж�Ӧ��ֵ�ķ�Χ
      %set( h,'ticks',(1:5:10),'fontsize',8,'ticklabels',{'<0',(1:5:10),'>15'});%����ɫ�������ϵĿ̶�ֵ
      %title(file(n).name);
        print([num2str(file(n).name),'-2','.tif'] ,'-dtiffn','-r150');
            FileName_New.a = file(n).name;
            FileName_New.b = '-Resorted.txt';                          % �ļ�������Ӻ�׺�����������������⸲��
            FileName_New.c = strcat(FileName_New.a,FileName_New.b);
             save(FileName_New.c,'fr_0_Ex','-ascii');             % save(file(n).name,'fr_Ex_sorted','-ascii');[Old Version]
      
  end
end


if(0)
    %% Scale
        fr_time = importdata('TimeScale.txt');
        fr_time_min = min(min(fr_time));
        fr_time = (fr_time - ones(size(fr_time)).* fr_time_min).*3.33;
        fr_time = [0;reshape(fr_time,[],1)];
    %% Scale
        fr_whole = importdata('S15-5um-P_Coh_data_fft.txt');
        fr_time_1 = ( fr_time(1:428,:))';
        fr_whole_rescaled = [fr_time_1;fr_whole(2:end,:)];
        % sortrows(data,1)
        fr_whole_sorted = sortrows(fr_whole_rescaled',1);
        
    save('S15-5um-P_Coh_data_fft2.txt','fr_whole_sorted','-ascii');

end  