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
%               佛祖保佑         永无BUG

% 20220412版
% 对CohData输出数据进行画图
% 添加了自动提取最大值，及根据极值拟合相干时间的函数

clc;
clear;
close all; 
file=dir('F:\matlab_project\Year2022-Q2\2022.04.04_Coh_PVSK_PhC_PowDep\Result_FFT\*.txt');  %处理文件的地址
n_max = length(file);
Freq_FFT_Zero = 0.5;   %  零级频率的范围（比如0~0.5） (1删除/0保留)
    ZeroExcavate = 1;
    CurveFitting = 1;  %  此处判定是否删除零级频率 (1删除/0保留)
    
   %% Step0_ Add Time 添加时间
        fr_time = importdata('TimeScale-1-10um-31×4.txt');                  %转换次序No与位置Position
        % fr_time_min = min(min(fr_time));
        fr_time = (fr_time - (ones(size(fr_time)).* min(min(fr_time)))).*3.33;  %DelayLine上位置和时间的转换
        fr_time = [0;reshape(fr_time,[],1)];                                    %重构DelayTime

   %% Step1_开始批量处理
for n =1:n_max
     strposition = strfind(file(n).name, 'S09-09003-P_Coh_FFTdata');  % 关键字选取
  if ~isempty (strposition)
     fr_0_Ex = importdata(file(n).name);  % 导入主文件数据
    
          [size_x,size_y]=size(fr_0_Ex);
          fr_time_1 = ( fr_time(1:size_y,:))';
          fr_whole_rescaled = [fr_time_1;fr_0_Ex(2:end,:)];   % 插入时间坐标 
          fr_Ex_sorted = sortrows(fr_whole_rescaled',1);      % 根据时间数据重排 (sortrows)
          fr_0_Ex = fr_Ex_sorted;
          
     fr_0_x = fr_0_Ex(1,2:end);
     fr_0_y = fr_0_Ex(2:end,1);
     fr_0_Ex = fr_0_Ex(2:end,2:end);
     
         Size_fr_0.x = length(fr_0_x);
         Size_fr_0.y = length(fr_0_y);
         xx_res =Size_fr_0.x;             %x轴坐标分辨率
         yy_res =Size_fr_0.y;             %y轴坐标分辨率
         n_max_res = xx_res * yy_res; %总数值个数 
         
    %%转换
        fr_0_Ex_one = ones(size(fr_0_Ex));   %构建和Ex相同结构的全1矩阵
        fr_1_x = fr_0_Ex_one' .* fr_0_x';          %fr_0_Ex_one' .* fr_0_x';
        fr_1_y = fr_0_y' .* (fr_0_Ex_one)';
        fr_1_Ex = fr_0_Ex';
        
    %%删除零级数据    
      if ZeroExcavate==1;  
        % fr_0_x 中小于N的数据
            [~, Position_Freq_FFT_Zero]=min(abs(fr_0_x - Freq_FFT_Zero));     % 寻找最靠近Freq_FFT_Zero位置的位置
            %(1:6)行删除 / 删除0.5频率前的数据
            fr_1_x = fr_1_x((Position_Freq_FFT_Zero+1):end,:);
            fr_1_y = fr_1_y((Position_Freq_FFT_Zero+1):end,:);
            fr_1_Ex = fr_1_Ex((Position_Freq_FFT_Zero+1):end,:);
      end
      
      %%寻找Matrix中最大值的位置
        fr_1_Ex_Max = max(max(fr_1_Ex));
        [fr_1_Ex_x fr_1_Ex_y] = find (fr_1_Ex_Max == fr_1_Ex);
        
      %%根据Matrix中最大值确定对应的频率位置，提取对应的频率曲线，进行拟合
      if CurveFitting==1;    
            fr_Curve.Ex = fr_1_Ex((fr_1_Ex_x-2):(fr_1_Ex_x+2),:); % curve数据
            fr_Curve.x = fr_1_x((fr_1_Ex_x-2):(fr_1_Ex_x+2),1);   % curve对应的频率信息
            fr_Curve.y = fr_1_y(1,:);                             % curve对应的时间信息
            fr_Curve_total = [fr_Curve.y;fr_Curve.Ex];

%             FileName_a = file(n).name;
%             FileName_b = '-Fitting.txt';                          % 文件名后添加后缀方便批量处理，并避免覆盖
%             FileName_c = strcat(FileName_a,FileName_b);
%              save(FileName_c,'fr_Curve_total','-ascii');             % save(file(n).name,'fr_Ex_sorted','-ascii');[Old Version]
            
           if(1) % 拟合
                   [Line Row] = size (fr_Curve_total);     % n_max = Line;  % 读取数据行数
                   Autoguess = 1; % 使用自动估值
                   % A= 1; Omega=2; Theta=0; C=0;
                   Para_x_c= 150; Para_w=500; Para_y_0 =0; Para_Aera=10000;  %
                   Parameter_Sum = [];
                   Fitting_Curve_Sum = [];

                if (1)
                    syms t; %设定自变量
                    n_fitting_max = Line;
                    for n_fitting =2:n_fitting_max;        % 第一行是时间坐标
                         x = fr_Curve_total (1,:) ;    % 提取时间坐标 
                         y = fr_Curve_total (n_fitting,:);   % 读取需要拟合的数据

                   %% 构建函数 括号内依次为自定义的函数、自变量t、待定参数    % 拟合曲线：A*sin(omega.*t + theta) + C
                         % Gauss
                          ft=fittype('y_0 + A/(w.*(sqrt(pi./2)))* exp(-2*((t-x_c)^2)./((w)^2))','independent','t','coefficients',{'y_0','w','x_c','A'});  
                         % Lorentz
                         % ft=fittype('y_0 + 2A/pi*(w/(4*(t-x_c)^2 + w^2))','independent','t','coefficients',{'y_0','w','x_c','A'});  % Lorentz

                        if(1) % 拟合曲线初值设定
                                  if Autoguess == 1;    
                                      Para.y_0 = 0 ;  
                                      Para.w =  Para_w;
                                      Para.x_c =  fr_1_y(1,fr_1_Ex_y); % Para_x_c;   
                                      Para.A = Para_Aera;                
                                  else
                                      Para.A = A;             %最小值和最大值相距的一般作为幅值
                                      Para.Omega = Omega;
                                      Para.Theta = Theta;
                                      Para.C = C;  
                                  end
                            opts = fitoptions(ft); 
                            opts.StartPoint=[Para.y_0 Para.w Para.x_c Para.A]; % 初始参数代入
                        end    

                        % 开始拟合
                        fittedfun=fit(x',y',ft,opts);%求出拟合得到的函数  
                            xi = [0:(max(x)./1000):max(x)];   % 1000个拟合点
                            yi = fittedfun(xi');              % 拟合结果

                         Parameter_Sum_i = [fittedfun.x_c, fittedfun.w, fittedfun.y_0, fittedfun.A];  % 拟合参数
                         Fitting_Sum_i = [yi']; %拟合曲线
                         % 分步存储拟合数据
                         Parameter_Sum = [Parameter_Sum; Parameter_Sum_i];
                         Fitting_Curve_Sum = [Fitting_Curve_Sum; Fitting_Sum_i];
                    end
            %%整理拟合数据
                    Parameter_Sum = [fr_Curve.x,Parameter_Sum];
                    Fitting_Curve_Sum = [xi;Fitting_Curve_Sum]';
                    fr_Curve_total=fr_Curve_total';
                       % 不同尺寸的数据如何整合到同一个矩阵中  (添加 fr_Curve_total）
                           % [Size_Para_Sum.x Size_Para_Sum.y] = size(Parameter_Sum);        
                       [Size_Fit_Cure_Sum.x Size_Fit_Cure_Sum.y] = size(Fitting_Curve_Sum);  
                       [Size_fr_Cure_Tot.x Size_fr_Cure_Tot.y] = size(fr_Curve_total);
                       fr_Cure_Output =zeros([max(Size_Fit_Cure_Sum.x,Size_fr_Cure_Tot.x) (Size_Fit_Cure_Sum.y+Size_fr_Cure_Tot.y)]);
                       fr_Cure_Output(1:Size_fr_Cure_Tot.x,1:Size_fr_Cure_Tot.y) = fr_Curve_total ;  
                       fr_Cure_Output(1:Size_Fit_Cure_Sum.x,(Size_fr_Cure_Tot.y+1):end) = Fitting_Curve_Sum ; 
            %%保存数据   
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
                  
      %%拟合曲线归一化           
      [fr_1_Ex_x, fr_1_Ex_y]= size(fr_1_Ex);
%          n_max_y = length(fr_0_y);
%          xx_res =n_max_x;             %x轴坐标分辨率
%          yy_res =n_max_y;             %y轴坐标分辨率
         n_max_res_2 = fr_1_Ex_x * fr_1_Ex_y;
            
            fr_2_x =reshape(fr_1_x,n_max_res_2,1);
            fr_2_y =reshape(fr_1_y,n_max_res_2,1);
            fr_2_Ex =reshape(fr_1_Ex,n_max_res_2,1);
    %%合并
        fr_total = [[fr_2_x';fr_2_y';fr_2_Ex']'];
        
   %%画图        
      figure
      %yy_region = fr_x_0 * fr_Ex_one fr_x_0;
      yy_region = fr_1_x;
      xx_region = fr_1_y;
      zz_region = fr_1_Ex;
      pcolor(yy_region,xx_region,zz_region);shading interp%伪彩色图
      set(gcf,'Colormap',hot)
       % xlim([440 570])  ylim([-10 1000]) 

      xlabel('x (μm)','FontName','Arial','FontSize',11);
      ylabel('Delay time (fs)','FontName','Arial','FontSize',11);%设置xy轴标签内容和字体
      set(gca, 'Fontname', 'Arial', 'Fontsize', 11);%设置xy轴的字体类型和字号大小的
      set(gca,'Position',[.18 .17 .6 .75]);%这句是设置xy轴在图片中占的比例，可能需要自己微调。
                 % zz_region_1 = reshape(zz_region,1,n_max_res_2);
                 % Zmax= max(zz_region_1);
                 % set(gca,'CLim',[0 Zmax]); %设置color bar 显示范围
      set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);%这句是将线宽改为2
      %set(gca,'Ytick',[-30:30:60]);%设置xy轴的刻度显示步长 %[Grating600_Center530nm]
      %set(gca,'Xtick',[-60:30:60]);
%               xlim([0 2])
%               set(gca,'Xtick',[0:0.5:2]);
      set(gcf,'Position',[100 100 400 320]);%这句是设置绘图的大小，不需要到word里再调整大小。我给的参数，图的大小是7cm
      colorbar('position',[0.82 0.17 0.035 0.75]); %设置彩色条的位置和大小
      %caxis([0,10]); %设置色度范围，表示真实的着色图中对应的值的范围
      %set( h,'ticks',(1:5:10),'fontsize',8,'ticklabels',{'<0',(1:5:10),'>15'});%设置色度条边上的刻度值
      %title(file(n).name);
        print([num2str(file(n).name),'-2','.tif'] ,'-dtiffn','-r150');
            FileName_New.a = file(n).name;
            FileName_New.b = '-Resorted.txt';                          % 文件名后添加后缀方便批量处理，并避免覆盖
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