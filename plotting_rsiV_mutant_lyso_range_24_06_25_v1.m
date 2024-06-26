function plotting_rsiV_mutant_lyso_range_24_06_25_v1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%What to plotswitch_frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all;
plot_mean_individual_strain=1; %plots mean data with each strain in a subplot, different conditions are plotted in different colours, repeats are plotted as individual lines with different line specs
plot_mean_all=1; %as plot_mean_individual_strain but plots mean of all repeats with each strain in a subolot, in different colours, only the mean for a condition and strain is plotted.

plot_mean_individual=0;
plot_single_trace=0;
plot_single_trace2=0;
plot_mean_mean=0;
plot_mean_mean_zoom=0;
plot_10xV=0;
plot_all_data=0;
plot_all_single=1;%%%%%%%%%%%%%%%
plot_all_data_part=0;
plot_hetero_cum=1;
plot_hetero_hist=0;
plot_mean_fraction=0;
plot_single_trace_15x=0;
plot_single_trace_15x_single=0;

export_data=0;

%%%%%%%%%%%%%%%%%%%%%%
% Swtich input
%%%%%%%%%%%%%%%%%%%%%%%%%
switch_frame=[101,101,102,101,101,...
    107,107,107,107,107,107,107];%frames when switching happend
adjust_switch=1; %1 when for different switching is adjusted 0 otherwise
adjust_to=101; %Frame to which to adjust
last_frame=220;%289 %last frame of movie/ to analyize

%%%%%%%%%%%%%%%%%%%%%%%%%
% Kill input
%%%%%%%%%%%%%%%%%%%%%%%%%
kill=1; %1= Removing repeats with problems 0 otherwise
max_cells_rep=40; %min number of cells in last frame
max_cells_short_gr=10;%Use to calculate how many slow growing cells are acceptable in one repeat


%Getting path with all data
p = mfilename('fullpath');
f=strfind(p,'\');
data_path_main=[p(1:f(end-2))];

%%%%%%%%%%%%%%%%%%%%%%%%%
% Data loading characteristics
%%%%%%%%%%%%%%%%%%%%%%%%

%Name of days with repeats
data_day={'2024-03-11','2024-03-19','2024-03-22','2024-04-01','2024-04-02','2024-04-23','2024-05-21','2024-05-28','2024-05-30','2024-06-03','2024-06-10','2024-06-12'};
%data_day={'2024-03-22','2024-04-01','2024-04-02','2024-04-23','2024-05-21','2024-05-28','2024-05-30','2024-06-03','2024-06-10','2024-06-12'};
%data_day={'2024-04-23','2024-05-21','2024-05-28','2024-05-30','2024-06-03','2024-06-10','2024-06-12'};%data 0.1 0.5
%
%data_day={'2024-05-21'};
%data_day={'2024-04-23'};
%data_day={'2024-05-28'};
%data_day={'2024-06-03'};
%data_day={'2024-06-10'};
%data_day={'2024-06-12'};
max_repeats=21; %maximum number of repeats
repeat_name={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21'}; %name of repeats
condition={'01','05'}; %name of conditions 01 stands for 0.1 ug/ml lyzosyme and 05 for 0.5 ug/ml lysozyme
condition={'01','05','1'};
%condition={'1'};
strains={'JLB130','JLB293','JLB327','JLB82'}; %Name of strains in strain list
t_names= {'WT','5xrsiV amyE','10xrsiV +sigV','15xrsiV +sigV'};%More useful name of strains


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_what='MY'; %Which property to plot
plot_what='MR';%Which property to plot
plot_what='elong_rate';%Which property to plot
cmap=distinguishable_colors(10); %colormap for plott
repeat_line={'-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.','-','--',':','-.','-'};%Defining Line prop
axis_y_max=[1500,3000,4000,500,1500,2500,1000,3000,3500,600,1000,2000]; %Max y for individual plots
title_names_simple=1; %1 for simple title; 0 for more info into title

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Defining variables
all_data=cell(max_repeats,length(strains),length(condition));
all_data_names=cell(max_repeats,length(strains),length(condition));
ind_kill=zeros(max_repeats,length(strains),length(condition));
num_loaded_cond=zeros(length(strains),length(condition));

%For loops over day, strain and condition
for day_now=1:length(data_day)
    for strain_now=1:length(strains);
        for cond_now=1:length(condition);

            %Check if condition exists
            D=dir([data_path_main,data_day{day_now},'\Data\',strains{strain_now},'_','*','_',condition{cond_now},'*']);
            if isempty(D);
                %If condition does not exist skip the rest
                %disp([data_path_main,data_day{day_now},'\Data\',strains{strain_now}])
                continue; 
            else

                %Going through all repeats of the day
                for repeat_do=1:length(D)
                    %displaying what is loaded
                    disp([data_path_main,data_day{day_now},'\Data\',D(repeat_do).name]); 

                    %Setting Repeat index
                    num_loaded_cond(strain_now,cond_now)=num_loaded_cond(strain_now,cond_now)+1;
                    rep_now=num_loaded_cond(strain_now,cond_now);

                    %Saving_name in two different formats depending on
                    %all_data_names
                    if title_names_simple==0
                        all_data_names{rep_now,strain_now,cond_now}={[data_day{day_now},': ',strrep(D(repeat_do).name(1:end-9),'_',' ')];[' S:',strains{strain_now}(4:end),' R:',num2str(rep_now),' c:',condition{cond_now}]};
                    else
                        all_data_names{rep_now,strain_now,cond_now}={[' S:',strains{strain_now}(4:end),' R:',num2str(rep_now),' c:',condition{cond_now}]};
                    end
                    
                    %Adjusting loaded data for different switching times
                    if adjust_switch==0
                        %Loading data for no switching
                        all_data{rep_now,day_now,strain_now}=load([data_path_main,data_day{day_now},'\Data\',D(repeat_do).name]);
                    else
                        %Actual loading of data
                        data_temp=load([data_path_main,data_day{day_now},'\Data\',D(repeat_do).name]);

                        %aligning all data
                        fnames=fieldnames(data_temp); %Getting all fieldnames
                        adjust=switch_frame(day_now)-adjust_to; %Determining needed adjustment for this repeat

                        for fn=1:length(fnames)
                            %Getting all field names
                            temp_field_data=data_temp.(fnames{fn});
%                             MY_temp=data_temp.MY;
                            %calculating adjustment
                            temp_matrix=nan(size(temp_field_data));%defining variable
                            temp_matrix(1:(289-adjust),:)=temp_field_data((1+adjust):289,:);
                            data_temp.(fnames{fn})=temp_matrix;
                            clear temp_matrix;
                        end
    
                            %saving corrected data
                            all_data{rep_now,strain_now,cond_now}=data_temp;
                        
                        %making a matrix with conditions to kill:
                        data_now=all_data{rep_now,strain_now,cond_now}.(plot_what)(1:last_frame,:);
                        data_gr_now=all_data{rep_now,strain_now,cond_now}.elong_rate(1:last_frame,:);
                        m_gr_now=nanmean(data_gr_now,2);
                        goodones=~isnan(data_now(last_frame,:));
                        %two conditions 1. I need a min number of cells
                        % 2. a min ratio of normal growing cells
                        if sum(goodones)>max_cells_rep&&sum(m_gr_now<0.3)/sum(goodones)<max_cells_short_gr/max_cells_rep
                           ind_kill(rep_now,strain_now,cond_now)=1;
                        else
                           ind_kill(rep_now,strain_now,cond_now)=0;
                        end
                    end
                end
            end
        end
    end
end
%Killing bad repeats
if kill==1;
    all_data_temp=cell(max_repeats,length(strains),length(condition));
    all_data_names_temp=cell(max_repeats,length(strains),length(condition));
    ind_kill=logical(ind_kill);
    for strain_now=1:length(strains)
        for cond_now=1:length(condition)
            ind=1;
            for rep_now=1:size(all_data,1);
                if ind_kill(rep_now,strain_now,cond_now)==1
                    all_data_temp{ind,strain_now,cond_now}=all_data{rep_now,strain_now,cond_now};
                    all_data_names_temp{ind,strain_now,cond_now}=all_data_names{rep_now,strain_now,cond_now};
                    ind=ind+1;
                end
            end
        end
    end
    %saving corrected data to old variable
    all_data=all_data_temp;
    all_data_names=all_data_names_temp;
end

% disp('test');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting mean by condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_mean_individual_strain==1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Descrption
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plots mean data with each strain in a subplot, different conditions 
    % are plotted in different colours, repeats are plotted as individual 
    % lines with different line specs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Setting Figure
    figure('units','normalized','outerposition',[0 0 1 1]);

   %looping over conditions
    for strain_now=1:length(strains)
        for rep_now=1:length(repeat_name)
            for cond_now=1:length(condition)

                subplot(2,2,strain_now);
                %Setting subplot figure
                if rep_now==1
                    title(t_names{strain_now});
                    xlabel('Frames');
                    ylabel(plot_what);
                    box on;
                    %Different axis values depending on what is plotted
                    if strcmp(plot_what,'MY');
                        axis([0,250,0,2500]);%MY;
                    elseif strcmp(plot_what,'MR');
                        axis([0,250,0,1000]);
                    elseif strcmp(plot_what,'elong_rate')
                        axis([0,250,0.2,1]);%Elong
                    end   
                end
                hold on;
                %checking if there is somehting to plot
                if isempty(all_data{rep_now,strain_now,cond_now})
                    continue;
                end
                %Actual Plotting
                plot(nanmean(all_data{rep_now,strain_now,cond_now}.(plot_what)(1:last_frame,:),2),'color',cmap(cond_now,:),'LineStyle', repeat_line{rep_now});
            end
        end
    %Setting Figure Title
    sgtitle(strrep(plot_what,'_',' '));
    if export_data==1
        exportgraphics(gcf,[data_path_main,'\Figures\mean_individual.png'],'Resolution',1200);
    end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting Mean of all data per condition and mutant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_mean_all==1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Description
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % As plot_mean_individual_strain but plots mean of all repeats with 
    % each strain in a subolot, in different colours, only the mean for a 
    % condition and strain is plotted.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Setting figure
    figure('units','normalized','outerposition',[0 0 1 1]);
    
    %loopting over conditions
    for strain_now=1:length(strains)
        for cond_now=1:length(condition)
            %Defining temp data and new index var
            temp_data=nan(1100,50);
            ind_now=1;
            for rep_now=1:length(repeat_name)
                subplot(2,2,strain_now);
                %Setting up figure
                if cond_now==1
                    title(t_names{strain_now});
                    xlabel('Frames');
                    ylabel(plot_what);
                    box on;
                    if strcmp(plot_what,'MY');
                        axis([0,250,0,2500]);%MY;
                    elseif strcmp(plot_what,'MR');
                        axis([0,250,0,1000]);
                    elseif strcmp(plot_what,'elong_rate')
                        axis([0,250,0.2,1]);%Elong
                    end
                end
                %only continue if there is data
                if isempty(all_data{rep_now,strain_now,cond_now})
                    continue;
                end
                
                %Calculating mean trace for this Repeat, strain and condition
                temp_data(1:last_frame,ind_now)=nanmean(all_data{rep_now,strain_now,cond_now}.(plot_what)(1:last_frame,:),2);
                ind_now=ind_now+1;
    
            end
            
            % Calculating mean of condition and plotting this it with errorbars
            hold on;
            errorbar(mean(temp_data,2,'omitnan'),std(temp_data,0,2,'omitnan')/sqrt(sum(~isnan(temp_data(1,:)))),'color',cmap(cond_now,:))
        end
    sgtitle(strrep(plot_what,'_',' '));

    if export_data==1
        exportgraphics(gcf,[data_path_main,'\Figures\mean_individual.png'],'Resolution',1200);
    end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting individal repeats in subplots with single cell traces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_all_single==1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Description
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting individual repats in their own subplot. Each plot has single
    % cell traces.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %parameters
    good_traces_only=1;%plots only traces longer than cutoff
    %cut_off=40;
    cut_off=max_cells_rep;
    thresh_gr=0.3;
    num_loaded_good_now=zeros(length(strains),length(condition));
    data_do_now=cell(10,length(strains),length(condition));
   
    if good_traces_only==1
        for strain_now=1:length(strains)
            for cond_now=1:length(condition)
                for rep_now=1:num_loaded_cond(strain_now,cond_now);
                    data_in=all_data{rep_now,strain_now,cond_now};
                    %data_now=data_in.MY;
                    data_now=data_in.(plot_what);

                    %Checking that growth rate is good
                    data_gr=data_in.elong_rate;
                    mean_gr=nanmean(data_gr,2);
                    fg=(mean_gr<thresh_gr);

                    goodones=~isnan(data_now(last_frame,:));
                    if sum(goodones)>cut_off&&sum(fg)/sum(goodones)<max_cells_short_gr/50;
                        num_loaded_good_now(strain_now,cond_now)=num_loaded_good_now(strain_now,cond_now)+1;
                        data_do_now{num_loaded_good_now(strain_now,cond_now),strain_now,cond_now}=data_in;
                        all_data_names_now{num_loaded_good_now(strain_now,cond_now),strain_now,cond_now}=all_data_names{rep_now,strain_now,cond_now};
                    else
                        %disp('shit');
                    end
                end
            end
        end
    end

    %First checking max number of repeats
    max_done_repeats=max(num_loaded_cond(:));
    max_cond=length(strains)*length(condition);
    ind_cond=0;
    %keep_rep=0;
    goodones=0;

    if good_traces_only==1
        num_loaded_cond_now=num_loaded_good_now;
        all_data_now=data_do_now;
        max_done_repeats_here=max(num_loaded_cond_now(:));
    else
        all_data_now=all_data;
        num_loaded_cond_now=num_loaded_cond;
        max_done_repeats_here=max_done_repeats;
    end
    
    figure;
   
    for strain_now=1:length(strains)
        for cond_now=1:length(condition)
            ind_cond=ind_cond+1;
            %ind=ind_cond;
            keep_rep=0;
            for rep_now=1:num_loaded_cond_now(strain_now,cond_now);
%                 %Setting index
%                 if keep_rep~=0
%                     ind=keep_rep;
%                 else
%                     ind=ind_cond+(rep_now-1)*max_cond;
%                 end
               ind=ind_cond+(rep_now-1)*max_cond;

                %Getting data
                data_in=all_data_now{rep_now,strain_now,cond_now};
                %data_now=data_in.MY;
                data_now=data_in.(plot_what);

                %Cleaning data
                goodones=~isnan(data_now(last_frame,:));
                data_now=data_now(:,goodones);

                %Plotting
                
%                 if sum(goodones)<30
%                     keep_rep=ind;
%                     continue;
%                 else
%                     keep_rep=0;
%                 end
                subplot(max_done_repeats_here,max_cond,ind);
                plot(data_now);
                if strcmp(plot_what,'MY')
                    axis([0,270,0,axis_y_max(ind_cond)]);
                elseif strcmp(plot_what,'MR')
                    axis([0,270,0,1500]);
                elseif strcmp(plot_what,'elong_rate')
                    axis([0,270,0,1]);
                end
                %Setting title
                %title(all_data_names{rep_now,strain_now,cond_now});
                title(all_data_names_now{rep_now,strain_now,cond_now});

                text(0.1,0.8,['n: ',num2str(sum(goodones))],'Unit','normalize');
                
                %Setting labels
                if mod(ind,max_cond)==1
                    ylabel('MY');
                end
                if ind>max_cond*max_done_repeats_here-max_cond
                    xlabel('frames');
                end
                xline(107,'r');


            end
        end
    end
%     sgtitle(t_name3{i});
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting mean by condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_mean_individual==1;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    % do_now=[1,7:10]
    %do_now=[1,7,8];
    for day_now=1:length(data_day);
        for strain_now=1:length(strains);
            subplot(5,2,strain_now);
            if day_now==1
                title(t_names{strain_now});
                xlabel('Frames');
                ylabel(plot_what);
                box on;
                if strcmp(plot_what,'MY');
                    axis([0,250,0,2500]);%MY;
                elseif strcmp(plot_what,'MR');
                    axis([0,250,0,1000]);
                elseif strcmp(plot_what,'elong_rate')
                    axis([0,250,0.2,1]);%Elong
                end
                
            end
            hold on;
            if isempty(all_data{day_now,strain_now})
                continue;
            end
            %plot(nanmean(all_data{day_now,strain_now}.MY,2),'color',cmap(strain_now,:),'LineStyle', repeat_line{day_now});
            %plot(nanmean(all_data{day_now,strain_now}.elong_rate,2),'color',cmap(strain_now,:),'LineStyle', repeat_line{day_now});
            plot(nanmean(all_data{day_now,strain_now}.(plot_what),2),'color',cmap(strain_now,:),'LineStyle', repeat_line{day_now});
            %if day_now==length(data_day)

        end
    end
    sgtitle(strrep(plot_what,'_',' '));
    if export_data==1
        exportgraphics(gcf,[data_path_main,'\Figures\mean_individual.png'],'Resolution',1200);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting mean by day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_mean_mean==1;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    % do_now=[1,7:10]
    %do_now=[1,7,8];
    m=nan(length(data_day),last_frame,length(strains));
    for day_now=1:length(data_day);
        for strain_now=1:length(strains);
            %subplot(5,2,strain_now)
            hold on;
            m(day_now,:,strain_now)=nanmean(all_data{day_now,strain_now}.MY(1:last_frame,:),2);
        end
    end

    %Plotting
    for strain_now=1:length(strains);
        hold on;
        shadedErrorBar(1:last_frame,mean(m(:,:,strain_now),1,'omitnan'),std(m(:,:,strain_now),1,'omitnan'),{'markerfacecolor',cmap(strain_now,:)},1);
        %shadedErrorBar(1:last_frame,mean(m(:,:,strain_now),1,'omitnan'),std(m(:,:,strain_now),1,'omitnan'));
    end
    box on;
    set(gca,'Linewidth',2,'FontWeight','bold')
%     legend(t_name3(do_now));
%     legend(t_names);
    xlabel('frames');
    ylabel('MY');
    title('Mean Fluorescence');
    h=get(gca,'children');
    legend(flip([h(1:4:end)]),t_names,'location','northwest');
    sgtitle('Mean by day');
    if export_data==1
        exportgraphics(gcf,[data_path_main,'\Figures\mean_mean.png'],'Resolution',1200);
    end
end

%%%%%%%%%%%%%%%%5
% Mean by day zoom
%%%%%%%%%%%%%%%%%%%
if plot_mean_mean_zoom==1;
       
    figure('units','normalized','outerposition',[0 0 1 1]);
    % do_now=[1,7:10]
    %do_now=[1,7,8];
    for day_now=1:length(data_day);
        for strain_now=1:length(strains);
            subplot(5,2,strain_now)
            hold on;
            plot(nanmean(all_data{day_now,strain_now}.MY-200,2),'color',cmap(strain_now,:),'LineStyle', repeat_line{day_now});
            
            if day_now==length(data_day)
                title(t_names{strain_now});
                xlabel('frames');
                ylabel('MY');
                box on;
                axis([90 110 0 200]);
            end
        end
    end
    if export_data==1
        exportgraphics(gcf,[data_path_main,'\Figures\mean_individual_zoom.png'],'Resolution',1200);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Checking 10x rsiV psigV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_10xV==1;
       
    figure('units','normalized','outerposition',[0 0 1 1]);
        for day_now=1:length(data_day);
            for strain_now=8
                MY_now=all_data{day_now,strain_now}.MY;

                goodones=~isnan(MY_now(last_frame,:));
                MY_now=MY_now(1:last_frame,goodones);

                subplot(3,1,day_now)
                plot(MY_now);
                axis([0,250,0,5000])

                ylabel('MY');
                xlabel('Frames');
                if day_now==1
                    sgtitle(t_names{strain_now});
                end
                title(data_day{day_now});
   
    
            end
        end
        if export_data==1
            exportgraphics(gcf,[data_path_main,'\Figures\10xV.png'],'Resolution',1200);
        end

        %Plotting histogram in last frame
        figure('units','normalized','outerposition',[0 0 1 1]);
        for day_now=1:length(data_day);
            for strain_now=8
                MY_now=all_data{day_now,strain_now}.MY;

                goodones=~isnan(MY_now(last_frame,:));
                MY_now=MY_now(last_frame,goodones);

                subplot(3,1,day_now)
                histogram(MY_now,[0:100:3000],'Normalization','probability');

                axis([0,3000,0,0.6])

                ylabel('MY');
                xlabel('Frames');
                if day_now==1
                    sgtitle(t_names{strain_now});
                end
                title(data_day{day_now});
    
            end
        end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting single trace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_single_trace==1;
    zoom_in=0;
    
    for i=1:2
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        ind=1;
        for day_now=1:length(data_day);
            for strain_now=1:length(strains);
                if isempty(all_data{day_now,strain_now})
                    continue;
                end
                MY_now=all_data{day_now,strain_now}.MY;
                if i>1
                    goodones=~isnan(MY_now(last_frame,:));
                    MY_now=MY_now(1:last_frame,goodones);
                end
                subplot(6,5,ind)
                plot(MY_now);
                if zoom_in==0
                    axis([0,290,0,5000]);
                elseif zoom_in==1;
                    axis([90,110,0,1000]);
                end
                n_s=sum(~isnan(MY_now(1,:)));
                n_l=sum(~isnan(MY_now(last_frame,:)));
                
                a=axis;
                text(a(2)*0.05,a(4)*0.8,['#start: ',num2str(n_s)],'FontSize',9);
                text(a(2)*0.05,a(4)*0.6,['#last: ',num2str(n_l)],'FontSize',9);
                if i==1
                    n_e=sum(~isnan(MY_now(289,:)));
                    text(a(2)*0.05,a(4)*0.4,['#end: ',num2str(n_e)],'FontSize',9);
                end
                
                if mod(ind,5)==1
                    ylabel([data_day{day_now},' MY']);
                end
                if ind>25
                    xlabel('Frames');
                end
                if ind<6
                    title(t_names{strain_now});
                elseif ind>15&&ind<21
                    title(t_names{strain_now});
                end
    
                if strain_now==5
                    ind=ind+11;
                elseif strain_now==10
                    ind=ind-14;
                else
                    ind=ind+1;
                end
    
    
    %             if day_now==length(data_day)
    %                 title(t_names{strain_now});
    %                 xlabel('frames');
    %                 ylabel('MY');
    %                 box on;
    %                 axis([0,250,0,2500])
    %             end
            end
        end
        if export_data==1&&i==1
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_all.png'],'Resolution',1200);
        elseif export_data==1&&i==2
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_goodones.png'],'Resolution',1200);
        end
        if i==1
            sgtitle('All Data');
        else
            sgtitle(['Data survive ', num2str(last_frame)]);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting single trace more data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_single_trace2==1;
    zoom_in=0;
    
%     for i=1:2
    for i=2
        figure('units','normalized','outerposition',[0 0 1 1]);
        
%         ind=-4;
        ind=1;
        strain_done=1;
        start_pos=1;
        for strain_now=1:length(strains);
            for day_now=1:length(data_day);

%                 if strain_now==6&&day_now==1
%                     figure; 
%                     ind=1;
%                     start_pos=1;
%                     strain_done=strain_done+1;
%                 elseif strain_now>strain_done
%                     ind=start_pos+1;
%                     strain_done=strain_done+1;
%                     start_pos=start_pos+1;
%                 else
%                     ind=ind+5;
%                 end

                if isempty(all_data{day_now,strain_now})
                    continue;
                end
                MY_now=all_data{day_now,strain_now}.MY;
                if i>1
                    goodones=~isnan(MY_now(last_frame,:));
                    MY_now=MY_now(1:last_frame,goodones)-200;
                end
%                 subplot(5,5,ind)
                subplot(3,2,ind)
                ind=ind+1; %added for this
                plot(MY_now);
                 xline(99,'-','10^{-4} ug');%t0
                 xline(140,'-','10^{-3} ug');%1h
                 xline(177,'-','10^{-2} ug');%2h
                 xline(268,'-','10^{-1} ug');%3h
                 xline(412,'-','2*10^{-2} ug');%4h
                 xline(539,'-','5*10^{-2} ug');%5h


                if zoom_in==0
%                     axis([0,290,0,1000]);
                      axis([0,570,0,800]);
                elseif zoom_in==1;
                    axis([90,110,0,1000]);
                end
                n_s=sum(~isnan(MY_now(1,:)));
                n_l=sum(~isnan(MY_now(last_frame,:)));
                
                a=axis;
                text(a(2)*0.05,a(4)*0.8,['#start: ',num2str(n_s)],'FontSize',9);
                text(a(2)*0.05,a(4)*0.6,['#last: ',num2str(n_l)],'FontSize',9);
                text(a(2)*0.05,a(4)*0.4,['day: ',num2str(day_now)],'FontSize',9);
                text(a(2)*0.05,a(4)*0.2,['strain: ',num2str(strain_now)],'FontSize',9);
                if i==1
                    n_e=sum(~isnan(MY_now(289,:)));
                    text(a(2)*0.05,a(4)*0.4,['#end: ',num2str(n_e)],'FontSize',9);
                end
                
%                 if mod(ind,5)==1;
%                     ylabel([data_day{day_now},' MY']);
%                 end
%                 if ind>25
%                     xlabel('Frames');
%                 end
%                 if ind<6
%                     title(t_names{strain_now});
%                 elseif ind>15&&ind<21
%                     title(t_names{strain_now});
%                 end
                
                title(t_names{strain_now});
                xlabel('Frames');
                ylabel([data_day{day_now},' MY']);

    
    
    %             if day_now==length(data_day)
    %                 title(t_names{strain_now});
    %                 xlabel('frames');
    %                 ylabel('MY');
    %                 box on;
    %                 axis([0,250,0,2500])
    %             end
            end
        end
        if export_data==1&&i==1!
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_all.png'],'Resolution',1200);
        elseif export_data==1&&i==2
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_goodones.png'],'Resolution',1200);
        end
        if i==1
            sgtitle('All Data');
        else
            sgtitle(['Data survive ', num2str(last_frame)]);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_all_data==1;

    figure; 
    for i=1:length(D); 
        a=load([data_path,D(i).name]);
        hold on; 
        subplot(5,2,i)
        data_all=a.MY-200;
        goodones=~isnan(data_all(289,:));
        plot(data_all(:,goodones));
    
        xlabel('frames');
        ylabel('MY');
        title(t_name3{i});
        axis([0,300,0,5000]);
        xline(103);
        %plotting off
        data_all=a.MY-200;
        data_part=data_all(1:switch_frame,:);
        m=nanmean(data_part(:));
        s=nanmean(data_part(:));
        thresh=m+3*s;
        yline(thresh);
    
    end
end



        




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All data WT and 8x only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_all_data_part==1;
    
    figure; 
    ind=1;
    for i=[1,10]
        a=load([data_path,D(i).name]);
        hold on; 
        subplot(1,2,ind)
        ind=ind+1;
        data_all=a.MY-200;
        goodones=~isnan(data_all(289,:));
        plot(data_all(:,goodones));
    
        xlabel('frames');
        ylabel('MY');
        title(t_name3{i});
        axis([0,300,0,1000]);
        xline(103);
        %plotting off
        data_all=a.MY-200;
        data_part=data_all(1:switch_frame,:);
        m=nanmean(data_part(:));
        s=nanmean(data_part(:));
        thresh=m+3*s;
        %yline(thresh);
    
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turn on dynamics cummulative plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_hetero_cum==1;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    ii=0;
    for strain_now=1:length(strains)
        for cond_now=1:length(condition)
            ii=ii+1;
            for rep_now=1:length(repeat_name)
                %subplot(2,2,strain_now);
                subplot(4,3,ii);

                if rep_now==1
                    %title(t_names{strain_now});
                    title(['S: ',strains{strain_now},' C:',condition{cond_now}]);
                    xlabel('frames');
                    ylabel('Fraction');
                    box on;
                    axis([0,150,0,1])
                end
                if isempty(all_data{rep_now,strain_now,cond_now});
                    continue;
                end
                goodones=~isnan(data_now(last_frame,:));
                data_gr_now=all_data{rep_now,strain_now,cond_now}.elong_rate(1:last_frame,:);
                m_gr_now=nanmean(data_gr_now,2);
                
                if sum(goodones)<=cut_off&&sum(fg)/sum(goodones)>=max_cells_short_gr/50
                    continue;
                end
                data_all=all_data{rep_now,strain_now,cond_now}.MY-200;
                data_part=data_all(1:switch_frame(day_now),:);
                if strain_now==1
                    m=nanmean(data_part(:));
                    s=nanmean(data_part(:));
                    thresh=m+3*s;
                end
                
                data_part_on=data_all(switch_frame(day_now):last_frame,:);
                goodones=~isnan(data_part_on(end,:));
                data_part_on2=data_part_on(:,goodones);
                on_temp=sum(data_part_on2>thresh,2)/sum(goodones);
    
    
                hold on;
                plot(on_temp,'color',cmap(cond_now,:),'LineStyle', repeat_line{rep_now});
    %             if day_now==1
    %                 title(t_names{strain_now});
    %                 xlabel('frames');
    %                 ylabel('Fraction');
    %                 box on;
    %                 axis([0,150,0,1])
    %             end
                 sgtitle('Cumulative Activation');
            end
        end
    end
    if export_data==1
        exportgraphics(gcf,[data_path_main,'\Figures\mean_individual.png'],'Resolution',1200);
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%turn histogram  activation times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_hetero_hist==1;
    
    figure; 
    %do_now=[1,7:10];
    for i=1:length(D); 
    %for i=do_now
    
        a=load([data_path,D(i).name]);
        %plotting off
        data_all=a.MY-200;
        data_part=data_all(1:switch_frame,:);
        m=nanmean(data_part(:));
        s=nanmean(data_part(:));
        thresh=m+3*s;
    
        data_part_on=data_all(switch_frame:last_frame,:);
        goodones=~isnan(data_part_on(end,:));
        data_part_on2=data_part_on(:,goodones);
        data_temp=data_part_on2>thresh;
        
        on=ones(size(data_temp,2),1)*last_frame;
        for j=1:size(data_temp,2)
            data_min=min(find(data_temp(:,j)));
            if ~isempty(data_min)
                on(j)=min(find(data_temp(:,j)));
            end
        end
        %on_temp=sum(data_part_on2>thresh,2)/sum(goodones);
        subplot(5,2,i);
        histogram(on,[0:10:300],'normalization','probability');
    %     plot(on_temp,c{i},'Linewidth',2);
        set(gca,'Linewidth',2,'FontWeight','bold');
        xlabel('Frames after Switch');
        ylabel('Fraction on');
        a=axis;
        axis([0,a(2),0,0.8])
        title(t_name3{i});
        box on;
    
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting mean activity and fraction of activation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_mean_fraction==1;
    ind=1;
    %Calculating threshhold
    for day_now=1:length(data_day);
        for strain_now=1
            if isempty(all_data{day_now,strain_now});
                continue;
            end

            data_all=all_data{day_now,strain_now}.MY-200;
            data_part=data_all(1:switch_frame(day_now),:);

            m=nanmean(data_part(:));
            s=nanmean(data_part(:));
            thresh_temp(ind)=m+3*s;
            ind=ind+1;
        end
    end
    m_data=nan(length(data_day),length(strains));
    n1=nan(length(data_day),length(strains));
    %s_data=nan(length(data_day),length(strains));
    thresh=mean(thresh);
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    for day_now=1:length(data_day);
        for strain_now=1:length(strains);
            %subplot(5,2,strain_now)
            if isempty(all_data{day_now,strain_now});
                continue;
            end

            data_all=all_data{day_now,strain_now}.MY-200;
            data_part=data_all(1:switch_frame(day_now),:);
%             if strain_now==1
%                 m=nanmean(data_part(:));
%                 s=nanmean(data_part(:));
%                 thresh=m+3*s;
%             end
            
            data_part_on=data_all(switch_frame(day_now):last_frame,:);
            goodones=~isnan(data_part_on(end,:));
            data_part_on2=data_part_on(:,goodones);

            data_temp=data_part_on2>thresh;
            
            on=ones(size(data_temp,2),1)*last_frame;
            for j=1:size(data_temp,2)
                data_min=min(find(data_temp(:,j)));
                if ~isempty(data_min)
                    on(j)=min(find(data_temp(:,j)));
                end
            end
            n=histcounts(on,[0:10:300],'normalization','probability');
            n1(day_now,strain_now)=n(1);

            data_temp=data_all(last_frame-50:last_frame,:);
            m_data(day_now,strain_now)=nanmean(data_temp(:));


        end
    end

    plot_now=[1:7];
    plot_now=[1:10];
    %plot_now=[1,7,8,10];
    %Plotting steady state
    subplot(2,1,1);
    data_x=[0,1,3,5,8,10,15,16,17,18];

    %correcting for outliers
    %10x second repeat
    %m_data(2,6)=nan;


    m_data_plot=mean(m_data,1,'omitnan');
    s_data_plot=std(m_data,'omitnan');
    n_plot=sum(~isnan(m_data));
    errorbar([data_x(plot_now)],m_data_plot(plot_now)/m_data_plot(1),s_data_plot(plot_now)/[m_data_plot(1)],'-x','Linewidth',2);
    a=axis;
    %axis([0,15,0,3000]);
    if length(plot_now)==7
        axis([0,15,0,1.2]);
        set(gca,'XTick',[0:1:15]);
    elseif length(plot_now)==10
        axis([0,18,0,1.2]);
        set(gca,'XTick',[0:1:18]);
        labelx=get(gca,'XTicklabel');
        labelx(17:19)={'10+V','15+V1','15+V2'};
        set(gca,'XTickLabel',labelx)
    end


    grid on;
    title('Steady State Expression');
    xlabel('Copies rsiV');
    ylabel('MY(au)');
    
    %Plotting Fraction
    subplot(2,1,2);

    %correcting for outliers
    %10x second repeat
    %n1(2,6)=nan;

    n1_m_plot=mean(n1,1,'omitnan');
    n1_s_plot=std(n1,1,'omitnan');

    errorbar([data_x(plot_now)],n1_m_plot(plot_now),n1_s_plot(plot_now),'-x','Linewidth',2);
    a=axis;
   if length(plot_now)==7
        axis([0,15,0,1.2]);
        set(gca,'XTick',[0:1:15]);
    elseif length(plot_now)==10
        axis([0,18,0,1.2]);
        set(gca,'XTick',[0:1:18]);
        labelx=get(gca,'XTicklabel');
        labelx(17:19)={'10+V','15+V1','15+V2'};
        set(gca,'XTickLabel',labelx)
   end
    %set(gca, 'Linewidth',2)
    grid on;
    title('Fraction active after 100 min');
    xlabel('Copies rsiV');
    ylabel('Fraction (au)');

    if export_data==1
        exportgraphics(gcf,[data_path_main,'\Figures\mean_individual.png'],'Resolution',1200);
    end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting 15x rsiV sigV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_single_trace_15x==1;
    zoom_in=1;
    
    %CV=nan(3,2);

    for i=1
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        ind=1;
        for day_now=1:length(data_day);
            for strain_now=9:10
                MY_now=all_data{day_now,strain_now}.MY;
                if i>1
                    goodones=~isnan(MY_now(last_frame,:));
                    MY_now=MY_now(1:last_frame,goodones);
                end
                subplot(3,2,ind)
                plot(MY_now-200);

                if zoom_in==0
                    axis([0,290,0,5000]);
                elseif zoom_in==1;
                    axis([90,150,0,100]);
                end
                n_s=sum(~isnan(MY_now(1,:)));
                n_l=sum(~isnan(MY_now(last_frame,:)));
                

                    xline(101);

                
%                 a=axis;
%                 text(a(2)*0.05,a(4)*0.8,['#start: ',num2str(n_s)],'FontSize',9);
%                 text(a(2)*0.05,a(4)*0.6,['#last: ',num2str(n_l)],'FontSize',9);
%                 if i==1
%                     n_e=sum(~isnan(MY_now(289,:)));
%                     text(a(2)*0.05,a(4)*0.4,['#end: ',num2str(n_e)],'FontSize',9);
%                 end
                

                    ylabel([data_day{day_now},' MY']);


                    xlabel('Frames');

                if ind<3
                    title(t_names{strain_now});
                end
                ind=ind+1;
    
    
    
    %             if day_now==length(data_day)
    %                 title(t_names{strain_now});
    %                 xlabel('frames');
    %                 ylabel('MY');
    %                 box on;
    %                 axis([0,250,0,2500])
    %             end
            end
        end
        if export_data==1&&i==1
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_all.png'],'Resolution',1200);
        elseif export_data==1&&i==2
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_goodones.png'],'Resolution',1200);
        end

%        figure; 
%         if i==1
%             sgtitle('All Data');
%         else
%             sgtitle(['Data survive ', num2str(last_frame)]);
%         end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting 15x rsiV sigV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plot_single_trace_15x_single==1;
    
    
    %CV=nan(3,2);

    for i=1
        ind=1;
        for day_now=1:length(data_day);
            for strain_now=[1,9:10];
                figure('units','normalized','outerposition',[0 0 1 1]);
                ind=1;
                MY_now=all_data{day_now,strain_now}.MY;

                    goodones=~isnan(MY_now(last_frame,:));
                    MY_now=MY_now(1:289,goodones);
                end
                subplot(8,8,ind)
                for cell_num=1:64
                    subplot(8,8,cell_num)
                    plot(MY_now(101:289,cell_num)-200);
                    if mod(ind,8)==1;
                        ylabel(['MY']);
                    end
                    if ind>56
                        xlabel('Frames');
                    end
                end

                sgtitle([data_day{day_now},' ', strains{strain_now}]);

                
%                 a=axis;
%                 text(a(2)*0.05,a(4)*0.8,['#start: ',num2str(n_s)],'FontSize',9);
%                 text(a(2)*0.05,a(4)*0.6,['#last: ',num2str(n_l)],'FontSize',9);
%                 if i==1
%                     n_e=sum(~isnan(MY_now(289,:)));
%                     text(a(2)*0.05,a(4)*0.4,['#end: ',num2str(n_e)],'FontSize',9);
%                 end
                

                    ylabel([data_day{day_now},' MY']);


                    xlabel('Frames');

                if ind<3
                    title(t_names{strain_now});
                end
                ind=ind+1;
    
    
    
    %             if day_now==length(data_day)
    %                 title(t_names{strain_now});
    %                 xlabel('frames');
    %                 ylabel('MY');
    %                 box on;
    %                 axis([0,250,0,2500])
    %             end
            end
        end
        if export_data==1&&i==1
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_all.png'],'Resolution',1200);
        elseif export_data==1&&i==2
            exportgraphics(gcf,[data_path_main,'\Figures\all_single_cell_goodones.png'],'Resolution',1200);
        end

%         figure; 
%         if i==1
%             sgtitle('All Data');
%         else
%             sgtitle(['Data survive ', num2str(last_frame)]);
%         end
    end
end
