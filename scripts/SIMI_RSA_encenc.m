%% As a response to reviewer#2. RSA between encoding stimuli

savepath = '';
low_pass = 0.5; % preprocessed data with two sets of parameters 0.1-30Hz and 0.5-30Hz

for loaddata=1

savename = strcat(savepath,'ERP_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs

savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs

savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save rejection index


clearvars savepath savename loaddata
end 



for s = 1:30 %
    disp(s);
    data_img_cl = ERP(s).img;

    
    data_img_cl([2,31,34],:,:)=[];

    N= 50;
    for ichan=1:59
        for itrial=1:size(data_img_cl,3)
            data_img_cl_smo(ichan,:,itrial)=smoothdata(data_img_cl (ichan,:,itrial),'Gaussian',N);
        end
        

    end
    
    for i=1:size(data_img_cl_smo,3)
        %%%%%%%% enc %%%%%%%%%%%
        for j=1:59
            data_img_cl_smo_down(j,:,i)=downsample(squeeze(data_img_cl_smo(j,:,i)),5); % downsample and removal of the baseline
        end
    end
    clearvars i j

    
    %% Get all
    % only high rem trial
    Data_S11 = data_img_cl_smo_down(:,1:205,1:3:216);
    Data_S22 = data_img_cl_smo_down(:,1:205,2:3:216);
    Data_S33 = data_img_cl_smo_down(:,1:205,3:3:216);

    
    n_trial = size(Data_S11,3); %all trial, separate later on

    for i = 1:n_trial % all trials
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S11 with S22 %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S11 = Data_S11(:,:,i); S22 = Data_S22(:,:,i);
        R = fast_corr(repmat(flip(S11,2),1,205),repelem(S22,1,205));
        if R==0
            R = zeros(1,42025);
        end
        SIMI_enc(s).S11S22(:,:,i) = reshape(R,[205,205]); 
        clearvars RSA_shuff enc trl_off R
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S22 with S33 %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S33 = Data_S33(:,:,i); S22 = Data_S22(:,:,i);
        R = fast_corr(repmat(flip(S22,2),1,205),repelem(S33,1,205));
        if R==0
            R = zeros(1,42025);
        end
        SIMI_enc(s).S22S33(:,:,i) = reshape(R,[205,205]); 
        clearvars RSA_shuff enc trl_off R
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S33 with S22 %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S33 = Data_S33(:,:,i); S22 = Data_S22(:,:,i);
        R = fast_corr(repmat(flip(S33,2),1,205),repelem(S22,1,205));
        if R==0
            R = zeros(1,42025);
        end
        SIMI_enc(s).S33S22(:,:,i) = reshape(R,[205,205]); 
        clearvars RSA_shuff enc trl_off R
    end
    
end

savename = strcat(savepath,'SIMI_Enc_(Gaussian_5).mat');
save (savename,'SIMI_enc','-v7.3');%save SIMI


%% Check the data 

for loaddata=1
savename = strcat(savepath,'SIMI_Enc_(Gaussian_5).mat');
load (savename);%load ERPs

savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs

savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save rejection index


clearvars savepath savename loaddata
end 

for s=1:30
    rej_indS11 = reject_index_img(s,1:3:216);rej_indS22 = reject_index_img(s,2:3:216);rej_indS33 = reject_index_img(s,3:3:216);
    
    S11S22_high(:,:,s) = mean(SIMI_enc(s).S11S22(:,:,num_recall(s,:)==2 & rej_indS11==1 & rej_indS22==1),3);
    S11S22_low(:,:,s) =  mean(SIMI_enc(s).S11S22(:,:,(num_recall(s,:)==0 | num_recall(s,:)==1) & rej_indS11==1 & rej_indS22==1),3);
    S11S22_all(:,:,s) =  mean(SIMI_enc(s).S11S22(:,:,rej_indS11==1 & rej_indS22==1),3);
    
    S22S33_high(:,:,s) = mean(SIMI_enc(s).S22S33(:,:,num_recall(s,:)==2 & rej_indS33==1 & rej_indS22==1),3);
    S22S33_low(:,:,s) =  mean(SIMI_enc(s).S22S33(:,:,(num_recall(s,:)==0 | num_recall(s,:)==1) & rej_indS33==1 & rej_indS22==1),3);
    S22S33_all(:,:,s) =  mean(SIMI_enc(s).S22S33(:,:,rej_indS33==1 & rej_indS22==1),3);
    
    S33S22_high(:,:,s) = mean(SIMI_enc(s).S33S22(:,:,num_recall(s,:)==2 & rej_indS33==1 & rej_indS22==1),3);
    S33S22_low(:,:,s) =  mean(SIMI_enc(s).S33S22(:,:,(num_recall(s,:)==0 | num_recall(s,:)==1) & rej_indS33==1 & rej_indS22==1),3);
    S33S22_all(:,:,s) =  mean(SIMI_enc(s).S33S22(:,:,rej_indS33==1 & rej_indS22==1),3);
end 



[clusters, p_values, t_sums, permutation_distribution ] = permutest( S22S33_all, S11S22_all, true,0.05,1000,true);


Cluster=zeros(205,205);
Cluster(clusters{1})=1;
%Cluster(clusters{2})=2;
%Cluster(clusters{3})=3;





for fig=2
    

[~,~,~,t] = ttest(reshape(S22S33_all,[size(S22S33_all,1)*size(S22S33_all,2),size(S22S33_all,3)]),reshape(S11S22_all,[size(S11S22_all,1)*size(S11S22_all,2),size(S11S22_all,3)]),'Dim',2);


    figure(fig);
    subplot(1,3,1);
    imagesc(mean(S11S22_all,3));
    xticks([0 0.5 1 1.5 2 ]*102.4)
    xticklabels({'0','500','1000','1500','2000'})
    yticks([0 0.5 1 1.5 2]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('S2 [ms]')
    ylabel('S1 [ms]')
    title('S1 with S2')
    caxis([-0.45,0.45])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

    
    subplot(1,3,2);
    imagesc(mean(S22S33_all,3));
    xticks([0 0.5 1 1.5 2 ]*102.4)
    xticklabels({'0','500','1000','1500','2000'})
    yticks([0 0.5 1 1.5 2]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('S3 [ms]')
    ylabel('S2 [ms]')
    title('S2 with S3')
    caxis([-0.45,0.45])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    
    %[c_row1,c_col1] = find(Cluster==1);
    %[c_row2,c_col2] = find(Cluster==2);
    %[c_row3,c_col3] = find(Cluster==3);
    
    

    %b1 = boundary(c_row1,c_col1,1);
    %b2 = boundary(c_row2,c_col2,1);
    %b3 = boundary(c_row3,c_col3);
    
    [B,L,N,A] = bwboundaries(Cluster); 
    
    subplot(1,3,3);
    imagesc(reshape(t.tstat,[size(S22S33_all,1),size(S11S22_all,2)]));hold on
    
    for k=1:length(B)
        boundary = B{k}; 
        plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
    end
    %plot(c_col1(b1),c_row1(b1),'LineWidth',3,'color','k');
    %plot(c_col2(b2),c_row2(b2),'LineWidth',3,'color','k');
    %plot(c_col3(b3),c_row3(b3),'LineWidth',3,'color','k');
    hold off   
    xticks([0 0.5 1 1.5 2 ]*102.4)
    xticklabels({'0','500','1000','1500','2000'})
    yticks([0 0.5 1 1.5 2]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    ylabel('S11/S22 [ms]')
    xlabel('S22/S33 [ms]')
    title('t')
    caxis([-3,3])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
end 


% focusing on the cluster and ANOVA

for s=1:30
    data = S11S22_high(:,:,s);
    avecluster(s,1) = mean(data(Cluster==1)); clearvars data
    data = S22S33_high(:,:,s);
    avecluster(s,2) = mean(data(Cluster==1)); clearvars data
    data = S11S22_low(:,:,s);
    avecluster(s,3) = mean(data(Cluster==1)); clearvars data
    data = S22S33_low(:,:,s);
    avecluster(s,4) = mean(data(Cluster==1)); clearvars data
    
end 

for fig=1 % mean cluster value
data = [reshape(avecluster,[120,1])];
lable_stim = repmat([repmat(1,30,1);repmat(2,30,1)],2,1);
remfor_stim = [repmat(1,30*2,1);repmat(2,30*2,1)];

figure(10);clf;
C = repmat([1,0,0;0,0,1],2,1);
boxplot(data,{lable_stim,remfor_stim},'FactorGap',[30,5],'Widths',1.2,'OutlierSize',4,'Symbol','');hold on
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),C(j,:),'FaceAlpha',.6);
end
set(findobj(gca,'type','line'),'linew',1.5)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% add indi point rem
for stim=1:2
    indi_data=data(30*(stim-1)+1:30*stim);
    scatter(ones(size(indi_data)).*(stim*3.05-2.05+(rand(size(indi_data))-0.5)/5),indi_data,'b','filled');
end
% add indi point for
for stim=1:2
    indi_data=data(30*(stim-1)+61:30*stim+60);
    scatter(ones(size(indi_data)).*(stim*3.05-0.9+(rand(size(indi_data))-0.5)/5),indi_data,'r','filled');
end
hold off
set(lines, 'Color', 'k');
xticks([1.775 1.775+5.85*1  ])
xticklabels({'S1','S2'})
ylim([-0.1,0.5])
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
end 