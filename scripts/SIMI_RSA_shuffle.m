%% Shuffle of RSA results
% shuffle offset for trial specificity exploration, with correct cluster
% and fast corr
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

%%
savename = strcat(savepath,'RSA_cluster_aftersmo0.mat');
load(savename)

Cluster=zeros(205,359);
Cluster(clusters{1})=1;

%%
nperm = 100;
for s = 1:30 % 
    disp(s);
    data_img_cl = ERP(s).img;
    data_off_cl = ERP(s).off;
    
    data_img_cl([2,31,34],:,:)=[];
    data_off_cl([2,31,34],:,:)=[];
    
    N= 50;
    for ichan=1:59
        for itrial=1:size(data_img_cl,3)
            data_img_cl_smo(ichan,:,itrial)=smoothdata(data_img_cl (ichan,:,itrial),'Gaussian',N);
        end
        
        for itrial=1:size(data_off_cl,3)
            data_off_cl_smo(ichan,:,itrial)=smoothdata(data_off_cl(ichan,:,itrial),'Gaussian',N);
        end
    end
    
    for i=1:size(data_img_cl_smo,3)
        %%%%%%%% enc %%%%%%%%%%%
        for j=1:59
            data_img_cl_smo_down(j,:,i)=downsample(squeeze(data_img_cl_smo(j,:,i)),5); % downsample and removal of the baseline
        end
    end
    clearvars i j
    
    for i=1:size(data_off_cl_smo,3)
        %%%%%%%% off %%%%%%%%%%%
        for j=1:59
            data_off_cl_smo_down(j,:,i)=downsample(squeeze(data_off_cl_smo(j,:,i)),5); % downsample and removal of the baseline
        end
    end
    clearvars i j
    
    
    %% shuffle
    % only high rem trial
    Data_S11 = data_img_cl_smo_down(:,1:205,1:3:216);Data_S11=Data_S11(:,:,num_recall(s,:)==1 |num_recall(s,:)==0);
    Data_S22 = data_img_cl_smo_down(:,1:205,2:3:216);Data_S22=Data_S22(:,:,num_recall(s,:)==1 |num_recall(s,:)==0 );
    Data_S33 = data_img_cl_smo_down(:,1:205,3:3:216);Data_S33=Data_S33(:,:,num_recall(s,:)==1 |num_recall(s,:)==0 );
    Data_off = data_off_cl_smo_down(:,:,num_recall(s,:)==1 |num_recall(s,:)==0 );
    
    n_trial = size(Data_off,3); %all trial, separate later on
    
    % 2s only, shuffle
    for i_perm = 1:nperm
        fprintf('Subject %d Perm %d',s,i_perm);
        off_perm = randperm(n_trial);
        perm_register{s}(i_perm,:)=off_perm;
        
        Data_offperm = Data_off(:,:,off_perm);
        for i = 1:n_trial % all trials
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%% S11 with offset %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            enc = Data_S11(:,:,i); trl_off = Data_offperm(:,:,i);
            R = fast_corr(repmat(flip(enc,2),1,359),repelem(trl_off,1,205));
            if (sum(isnan(enc(:))) + sum(isnan(trl_off(:))))>0
                RSA_shuff = zeros(205,359);
            else
                RSA_shuff = reshape(R,[205,359]); 
            end
            Cluster_zone{s}(1,i,i_perm) = mean(RSA_shuff(Cluster==1));
            clearvars RSA_shuff enc trl_off R
 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%% S22 with offset %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            enc = Data_S22(:,:,i); trl_off = Data_offperm(:,:,i);
            R = fast_corr(repmat(flip(enc,2),1,359),repelem(trl_off,1,205));
            if (sum(isnan(enc(:))) + sum(isnan(trl_off(:))))>0
                RSA_shuff = zeros(205,359);
            else
                RSA_shuff = reshape(R,[205,359]); 
            end
            Cluster_zone{s}(2,i,i_perm) = mean(RSA_shuff(Cluster==1));
            clearvars RSA_shuff enc trl_off R

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%% S33 with offset %%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            enc = Data_S33(:,:,i); trl_off = Data_offperm(:,:,i);
            R = fast_corr(repmat(flip(enc,2),1,359),repelem(trl_off,1,205));
            if (sum(isnan(enc(:))) + sum(isnan(trl_off(:))))>0
                RSA_shuff = zeros(205,359);
            else
                RSA_shuff = reshape(R,[205,359]); 
            end  
            Cluster_zone{s}(3,i,i_perm) = mean(RSA_shuff(Cluster==1));
            clearvars RSA_shuff enc trl_off R

        end
     
        
    end
    
end

% data saving
save ([savepath 'RSA_shuffle_perm3_forg.mat'],'Cluster_zone','perm_register')


%%
% Real value
for realvalue =1
low_pass = 0.5; % preprocessed data with two sets of parameters 0.1-30Hz and 0.5-30Hz
for loaddata=1

savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs
savename = strcat(savepath,'SIMI_all_(Gaussian)_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save ERPs

savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save rejection index

savename = strcat(savepath,'RSA_cluster_aftersmo0.mat');
load(savename)

Cluster=zeros(205,359);
Cluster(clusters{1})=1;

clearvars savepath savename loaddata
end 

% recall num for each stim
for s=1:size(SIMI_all,2)
    recall_indperstim(s,:) = rep_num(num_recall(s,:),3,1); 
    reject_index_off_perstim(s,:) = rep_num(reject_index_off(s,:),3,1); 
end
% stim index
stim_ind = rep_num([1,2,3],1,72); % first second and third img always;

t_select=52:256;
S11_off_n0 = zeros(205,359,size(SIMI_all,2));S11_off_n1 = zeros(205,359,size(SIMI_all,2));S11_off_n2 = zeros(205,359,size(SIMI_all,2));
S22_off_n0 = zeros(205,359,size(SIMI_all,2));S22_off_n1 = zeros(205,359,size(SIMI_all,2));S22_off_n2 = zeros(205,359,size(SIMI_all,2));
S33_off_n0 = zeros(205,359,size(SIMI_all,2));S33_off_n1 = zeros(205,359,size(SIMI_all,2));S33_off_n2 = zeros(205,359,size(SIMI_all,2));
S11_off_low = zeros(205,359,size(SIMI_all,2));S22_off_low = zeros(205,359,size(SIMI_all,2));S33_off_low = zeros(205,359,size(SIMI_all,2));

for s=1:size(SIMI_all,2)
    S11_off_n0(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 0 & stim_ind == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S11_off_n1(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 1 & stim_ind == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S11_off_n2(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 2 & stim_ind == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    S22_off_n0(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 0 & stim_ind == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S22_off_n1(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 1 & stim_ind == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S22_off_n2(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 2 & stim_ind == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    S33_off_n0(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 0 & stim_ind == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S33_off_n1(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 1 & stim_ind == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S33_off_n2(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) == 2 & stim_ind == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    S11_off_low(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) ~= 2 & stim_ind == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S22_off_low(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) ~= 2 & stim_ind == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    S33_off_low(:,:,s) = mean(SIMI_all{s}(t_select,:,recall_indperstim(s,:) ~= 2 & stim_ind == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
end 

off_n0 =  (S11_off_n0 + S22_off_n0 + S33_off_n0)./3;
off_n1 =  (S11_off_n1 + S22_off_n1 + S33_off_n1)./3;
off_n2 =  (S11_off_n2 + S22_off_n2 + S33_off_n2)./3;
off_low = (S11_off_low + S22_off_low + S33_off_low)./3;

for s=1:30
    S1 = S11_off_n2(:,:,s);
    real_C_high(1,s) = mean(S1(Cluster==1));
    S2 = S22_off_n2(:,:,s);
    real_C_high(2,s) = mean(S2(Cluster==1));
    S3 = S33_off_n2(:,:,s);
    real_C_high(3,s) = mean(S3(Cluster==1));
    clearvars S1 S2 S3
    
    S1 = off_low(:,:,s);
    real_C_low(1,s) = mean(S1(Cluster==1));
    S2 = off_low(:,:,s);
    real_C_low(2,s) = mean(S2(Cluster==1));
    S3 = off_low(:,:,s);
    real_C_low(3,s) = mean(S3(Cluster==1));
    clearvars S1 S2 S3
    
  
end 

end 


%% Shuffle
for loadhide=1
    
    %% Rem
    load ([savepath 'RSA_shuffle_perm1_high.mat'])
    C1P1_high = Cluster_zone; O1_high = perm_register;
    load ([savepath 'RSA_shuffle_perm3_high.mat'])
    C1P2_high = Cluster_zone; O2_high = perm_register;
    %% Forg
    load ([savepath 'RSA_shuffle_perm1_forg.mat'])
    C1P1_low = Cluster_zone; O1_low = perm_register;
    load ([savepath 'RSA_shuffle_perm3_forg.mat'])
    C1P2_low = Cluster_zone; O2_low = perm_register;
    
end 


% High
for s = 1:30
    disp(s)
    perm_Cluster_high =cat(3,C1P1_high{s},C1P2_high{s});
    perm_r = cat(1,O1_high{s},O2_high{s});
    
    %% rejection index
    rej_S11 = reject_index_img(s,1:3:216);
    rej_S22 = reject_index_img(s,2:3:216);
    rej_S33 = reject_index_img(s,3:3:216);
    
    rej_S11 = rej_S11(num_recall(s,:)==2)==1;
    rej_S22 = rej_S22(num_recall(s,:)==2)==1;
    rej_S33 = rej_S33(num_recall(s,:)==2)==1;
    
    rej_off = reject_index_off(s,num_recall(s,:)==2);
    
    for i_perm = 1: size(perm_r,1)
        rej_off = rej_off(perm_r(i_perm,:))==1;
        
        Perm_C_high(1,s,i_perm)=mean(perm_Cluster_high(1,rej_S11 & rej_off,i_perm),2);
        Perm_C_high(2,s,i_perm)=mean(perm_Cluster_high(2,rej_S22 & rej_off,i_perm),2);
        Perm_C_high(3,s,i_perm)=mean(perm_Cluster_high(3,rej_S33 & rej_off,i_perm),2);
        
    end 
end 

% Low
for s = 1:30
    disp(s)
    perm_Cluster_low =cat(3,C1P1_low{s},C1P2_low{s});
    perm_r = cat(1,O1_low{s},O2_low{s});
    
    %% rejection index
    rej_S11 = reject_index_img(s,1:3:216);
    rej_S22 = reject_index_img(s,2:3:216);
    rej_S33 = reject_index_img(s,3:3:216);
    
    rej_S11 = rej_S11(num_recall(s,:)==1|num_recall(s,:)==0)==1;
    rej_S22 = rej_S22(num_recall(s,:)==1|num_recall(s,:)==0)==1;
    rej_S33 = rej_S33(num_recall(s,:)==1|num_recall(s,:)==0)==1;
    
    rej_off = reject_index_off(s,num_recall(s,:)==1|num_recall(s,:)==0);
    
    for i_perm = 1: size(perm_r,1)
        rej_off = rej_off(perm_r(i_perm,:))==1;
        
        Perm_C_low(1,s,i_perm)=mean(perm_Cluster_low(1,rej_S11 & rej_off,i_perm),2);
        Perm_C_low(2,s,i_perm)=mean(perm_Cluster_low(2,rej_S22 & rej_off,i_perm),2);
        Perm_C_low(3,s,i_perm)=mean(perm_Cluster_low(3,rej_S33 & rej_off,i_perm),2);
        
    end 
end 

Perm_ave_high =mean(mean(Perm_C_high,3),1);
Perm_ave_low =mean(mean(Perm_C_low,3),1);


%% 
figure(34);clf;
plot_data2 = [mean(real_C_high,1);Perm_ave_high;mean(real_C_low,1);Perm_ave_low]';
for fig=1
err = std(plot_data2);
b=bar(mean(plot_data2,1),0.8);hold on
er = errorbar(mean(plot_data2,1),err); 
scatter(ones(size(plot_data2(:,1))).*(1+(rand(size(plot_data2(:,1)))-0.5)/5),plot_data2(:,1),'k','filled');
scatter(ones(size(plot_data2(:,2))).*(2+(rand(size(plot_data2(:,2)))-0.5)/5),plot_data2(:,2),'k','filled');
scatter(ones(size(plot_data2(:,3))).*(3+(rand(size(plot_data2(:,3)))-0.5)/5),plot_data2(:,3),'k','filled');
scatter(ones(size(plot_data2(:,4))).*(4+(rand(size(plot_data2(:,4)))-0.5)/5),plot_data2(:,4),'k','filled');hold off
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.LineWidth = 3;
% b.FaceColor = 'flat';
% b.FaceAlpha = 0.6;
%     b.CData(1,:) = [0,0,1];
%     b.CData(2,:) = [1,0,0];
xticks([1 2 3 4])
xticklabels({'Real High','Perm High','Real Low','Perm Low'})
ylim([-0.2,0.2])
ylabel('Pearson´s R')
title('Cluster ')
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);
end 