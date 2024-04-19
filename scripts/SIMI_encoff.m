%% SIMI analysis enc off
% Xiongbo 2021-04-17
savepath = '';
low_pass = 0.5; % preprocessed data with two sets of parameters 0.1-30Hz and 0.5-30Hz

for loaddata=1

savename = strcat(savepath,'ERP_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs

savename = strcat(savepath,'SIMI_all_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save SIMI

clearvars savepath savename loaddata
end 

% alternative (faster)
for s = 1:30 % update newest participants only
    disp(['Subject_',int2str(s)])
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

    
    for i = 1:72 % all trials
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S11 with offset %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S11 = (i-1)*3+1;
        enc = data_img_cl_smo_down(:,1:205,ind_S11); trl_off = data_off_cl_smo_down(:,:,i);
        R = fast_corr(repmat(flip(enc,2),1,359),repelem(trl_off,1,205));
        if (sum(isnan(enc(:))) + sum(isnan(trl_off(:))))>0
            SIMI(:,:,ind_S11) = zeros(205,359);
        else
            SIMI(:,:,ind_S11) = reshape(R,[205,359]);
        end
        clearvars ind_S11 A1 A2 R enc trl_off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S22 with offset %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S22 = (i-1)*3+2;
        enc = data_img_cl_smo_down(:,1:205,ind_S22); trl_off = data_off_cl_smo_down(:,:,i);
        R = fast_corr(repmat(flip(enc,2),1,359),repelem(trl_off,1,205));
        if (sum(isnan(enc(:))) + sum(isnan(trl_off(:))))>0
            SIMI(:,:,ind_S22) = zeros(205,359);
        else
            SIMI(:,:,ind_S22) = reshape(R,[205,359]);
        end
        clearvars ind_S11 A1 A2 R enc trl_off
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S33 with offset %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S33 = (i-1)*3+3;
        enc = data_img_cl_smo_down(:,1:205,ind_S33); trl_off = data_off_cl_smo_down(:,:,i);
        R = fast_corr(repmat(flip(enc,2),1,359),repelem(trl_off,1,205));
        if (sum(isnan(enc(:))) + sum(isnan(trl_off(:))))>0
            SIMI(:,:,ind_S33) = zeros(205,359);
        else
            SIMI(:,:,ind_S33) = reshape(R,[205,359]);
        end
        clearvars ind_S11 A1 A2 R enc trl_off
        
    end
    
    SIMI_all{s}=SIMI;clearvars SIMI
end

% alternative (fastest) concat trials also
for s = 1:30 % update newest participants only
    disp(['Subject_',int2str(s)])
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

    enc_S11 = [];enc_S22 = [];enc_S33 = []; trl_off = [];
    for i = 1:72 % concat all trials
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S11  %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S11 = (i-1)*3+1;
        enc = data_img_cl_smo_down(:,1:205,ind_S11);  
        if sum(isnan(enc(:))) >0
            enc_S11 = horzcat(enc_S11,repmat(flip(zeros(size(enc)),2),1,359));
        else
            enc_S11 = horzcat(enc_S11,repmat(flip(enc,2),1,359));
        end
        clearvars ind_S11  enc 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S22  %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S22 = (i-1)*3+2;
        enc = data_img_cl_smo_down(:,1:205,ind_S22);  
        if sum(isnan(enc(:))) >0
            enc_S22 = horzcat(enc_S22,repmat(flip(zeros(size(enc)),2),1,359));
        else
            enc_S22 = horzcat(enc_S22,repmat(flip(enc,2),1,359));
        end
        clearvars ind_S22  enc 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S33  %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S33 = (i-1)*3+3;
        enc = data_img_cl_smo_down(:,1:205,ind_S33);  
        if sum(isnan(enc(:))) >0
            enc_S33 = horzcat(enc_S33,repmat(flip(zeros(size(enc)),2),1,359));
        else
            enc_S33 = horzcat(enc_S33,repmat(flip(enc,2),1,359));
        end
        clearvars ind_S33  enc 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% offset  %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        off = data_off_cl_smo_down(:,:,i);
        if sum(isnan(trl_off(:)))>0
            trl_off = horzcat(trl_off,repelem(zeros(size(off)),1,205));
        else
            trl_off = horzcat(trl_off,repelem(off,1,205));
        end
        clearvars  off
    end
    SIMI_all{s} = zeros(205,359,72);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S11 with offset %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = reshape(fast_corr(enc_S11,trl_off),[205,359,72]);
    SIMI_all{s}(:,:,1:3:216)=R;
    clearvars R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S22 with offset %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = reshape(fast_corr(enc_S22,trl_off),[205,359,72]);
    SIMI_all{s}(:,:,2:3:216)=R;
    clearvars R
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% S33 with offset %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R = reshape(fast_corr(enc_S33,trl_off),[205,359,72]);
    SIMI_all{s}(:,:,3:3:216)=R;
    clearvars R
    
    
end


% traditional method
for s = 1:30 % update newest participants only
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

    
    for i = 1:72 % all trials
        disp(['Subject_',int2str(s),' trial ',int2str(i)])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S11 with offset %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S11 = (i-1)*3+1;
        % first step(1:250 for both)
        A1=corr_matrix(data_off_cl_smo_down(:,1:256,i),data_img_cl_smo_down(:,:,ind_S11));
        % second step (1:250 enc with 51:300 for off)
        A2=corr_matrix(data_off_cl_smo_down(:,104:359,i),data_img_cl_smo_down(:,:,ind_S11));
        % third step: combine the matrix
        SIMI(:,:,ind_S11)=[A1 A2(:,154:256)]; clearvars ind_S11 A1 A2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S22 with offset %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S22 = (i-1)*3+2;
        % first step(1:250 for both)
        A1=corr_matrix(data_off_cl_smo_down(:,1:256,i),data_img_cl_smo_down(:,:,ind_S22));
        % second step (1:250 enc with 51:300 for off)
        A2=corr_matrix(data_off_cl_smo_down(:,104:359,i),data_img_cl_smo_down(:,:,ind_S22));
        % third step: combine the matrix
        SIMI(:,:,ind_S22)=[A1 A2(:,154:256)]; clearvars ind_S11 A1 A2
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% S33 with offset %%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ind_S33 = (i-1)*3+3;
        % first step(1:250 for both)
        A1=corr_matrix(data_off_cl_smo_down(:,1:256,i),data_img_cl_smo_down(:,:,ind_S33));
        % second step (1:250 enc with 51:300 for off)
        A2=corr_matrix(data_off_cl_smo_down(:,104:359,i),data_img_cl_smo_down(:,:,ind_S33));
        % third step: combine the matrix
        SIMI(:,:,ind_S33)=[A1 A2(:,154:256)]; clearvars ind_S11 A1 A2
        
    end
    
    SIMI_all{s}=SIMI;clearvars SIMI
end


% data saving
savepath = 'HERE PLACE YOUR PATH\';
savename = strcat(savepath,'SIMI_all_(Gaussian)_',num2str(low_pass),'Hzlowpass_group.mat');
save (savename,'SIMI_all','-v7.3');%save SIMI

%% analysis
low_pass = 0.5; % preprocessed data with two sets of parameters 0.1-30Hz and 0.5-30Hz
for loaddata=1

savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs
savename = strcat(savepath,'SIMI_all_(Gaussian)_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save ERPs

savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save rejection index

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



off_high_smo = zeros(0);off_low_smo = zeros(0);
N=0;
for s = 1:30
    switch N
        case 0
            off_high_smo(:,:,s) = off_n2(:,:,s);
            off_low_smo(:,:,s) = off_low(:,:,s);
        otherwise
            off_high_smo(:,:,s) = smooth2a(off_n2(:,:,s),N);
            off_low_smo(:,:,s) = smooth2a(off_low(:,:,s),N);
    end
end


%[Cluster,T_Cluster,T_per,p_Cluster]=xb_cluster2D_perm_ori(off_high_smo,off_low_smo,1000,'within','two',0.05);

%[Cluster,T_Cluster,T_per,p_Cluster]=xb_cluster2D_perm(off_high_smo,off_low_smo,1000,'within','two',0.1,1);

[clusters, p_values, t_sums, permutation_distribution ] = permutest( off_high_smo, off_low_smo, true,0.05,1000,true);
[Cluster2,T_Cluster2,T_per2,p_Cluster2]=clusterperm2D(off_high_smo,off_low_smo,1000,'within','two',0.05);

savename = strcat(savepath,'RSA_cluster_aftersmo',int2str(N),'.mat');
save (savename,'clusters','p_values','t_sums','permutation_distribution','-v7.3');

Cluster=zeros(205,359);
Cluster(clusters{1})=1;
%Cluster(clusters{2})=2;

% SIMI
for fig=1
    figure(11);
    imagesc(mean(off_high_smo,3));
    xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
    xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
    yticks([0 0.5 1 1.5 2 ]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('High Mem')
    caxis([-0.08,0.08])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

%     figure(22);
%     imagesc(mean(off_n1,3));
%     xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
%     xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
%     yticks([0 0.5 1 1.5 2 ]*102.4)
%     yticklabels({'2000','1500','1000','500','0'})
%     xlabel('Offset [ms]')
%     ylabel('Encoding [ms]')
%     title('Recall 1')
%     caxis([-0.08,0.08])
%     set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(33);
    imagesc(mean(off_low_smo,3));
    xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
    xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
    yticks([0 0.5 1 1.5 2 ]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Low Mem')
    caxis([-0.08,0.08])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    clearvars fig
    
    
    
%     figure(35);
%     imagesc(mean(off_low,3));
%     xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
%     xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
%     yticks([0 0.5 1 1.5 2 ]*102.4)
%     yticklabels({'2000','1500','1000','500','0'})
%     xlabel('Offset [ms]')
%     ylabel('Encoding [ms]')
%     title('Recall low')
%     caxis([-0.08,0.08])
%     set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
%     clearvars fig
    
    [~,~,~,t] = ttest(reshape(off_high_smo,[size(off_high_smo,1)*size(off_high_smo,2),size(off_high_smo,3)]),reshape(off_low_smo,[size(off_low_smo,1)*size(off_low_smo,2),size(off_low_smo,3)]),'Dim',2);
    
% 
%     [c_row1,c_col1] = find(Cluster==1);
%     [c_row2,c_col2] = find(Cluster==2);
%     %[c_row3,c_col3] = find(Cluster==3);
%     
%     b1 = boundary(c_row1,c_col1,1);
%     b2 = boundary(c_row2,c_col2,1);
    %b3 = boundary(c_row3,c_col3);
     [B,L,N,A] = bwboundaries(Cluster); 

    figure(37);
    imagesc(reshape(t.tstat,[size(off_low_smo,1),size(off_low_smo,2)]));hold on
     for k=1:length(B)
        boundary = B{k}; 
        plot(boundary(:,2), boundary(:,1), 'LineWidth',3,'color','k');
    end
    %plot(c_col2(b2),c_row2(b2),'LineWidth',3,'color','k');
    %plot(c_col3(b3),c_row3(b3),'LineWidth',3,'color','k');
    hold off
    xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
    xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
    yticks([0 0.5 1 1.5 2]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Rem vs For (t)')
    caxis([-2.5,2.5])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    clearvars fig

end 

t_map = reshape(t.tstat,[size(off_low_smo,1),size(off_low_smo,2)]);
tmean_cluster1 = mean(t_map(Cluster==1));
tmean_cluster2 = mean(t_map(Cluster==2));
tmax_cluster1 = max(t_map(Cluster==1));
tmax_cluster2 = max(t_map(Cluster==2));

High_SIMI = mean(off_high_smo,3);Low_SIMI = mean(off_low_smo,3);
savename = strcat(savepath,'SIMI_forPyhthonplot.mat');
save (savename,'High_SIMI','Low_SIMI','t_map','-v7.3');



% Mean cluster value plot
for fig=1


for i=1:30
    rem_S11=S11_off_n2(:,:,i);for_S11=S11_off_low(:,:,i);
    cluster_rem(i,1)=mean(rem_S11(clusters{1}));cluster_for(i,1)=mean(for_S11(clusters{1}));
    clearvars rem_S11 for_S11
    
    rem_S22=S22_off_n2(:,:,i);for_S22=S22_off_low(:,:,i);
    cluster_rem(i,2)=mean(rem_S22(clusters{1}));cluster_for(i,2)=mean(for_S22(clusters{1}));
    clearvars rem_S22 for_S22
    
    rem_S33=S33_off_n2(:,:,i);for_S33=S33_off_low(:,:,i);
    cluster_rem(i,3)=mean(rem_S33(clusters{1}));cluster_for(i,3)=mean(for_S33(clusters{1}));
    clearvars rem_S33 for_S33
    
end

data = [reshape(cluster_rem,[30*3,1]);reshape(cluster_for,[30*3,1])];
lable_stim = repmat([repmat(1,30,1);repmat(2,30,1);repmat(3,30,1)],2,1);
remfor_stim = [repmat(1,30*3,1);repmat(2,30*3,1)];

figure(10);clf;
C = repmat([1,0,0;0,0,1],3,1);
boxplot(data,{lable_stim,remfor_stim},'FactorGap',[30,5],'Widths',1.2,'OutlierSize',4,'Symbol','');hold on
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),C(j,:),'FaceAlpha',.6);
end
set(findobj(gca,'type','line'),'linew',1.5)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
% add indi point rem
for stim=1:3
    indi_data=data(30*(stim-1)+1:30*stim);
    scatter(ones(size(indi_data)).*(stim*3.75-2.75+(rand(size(indi_data))-0.5)/5),indi_data,'b','filled');
end
% add indi point rem
for stim=1:3
    indi_data=data(30*(stim-1)+91:30*stim+90);
    scatter(ones(size(indi_data)).*(stim*3.75-1.5+(rand(size(indi_data))-0.5)/5),indi_data,'r','filled');
end
hold off
set(lines, 'Color', 'k');
xticks([1.775 1.775+5.85*1 1.775+5.85*2 ])
xticklabels({'S1','S2','S3'})
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
end 




%% Category
F_off_cue = zeros(256,359,size(SIMI_all,2));F_off_rem = zeros(256,359,size(SIMI_all,2));F_off_for = zeros(256,359,size(SIMI_all,2));
P_off_cue = zeros(256,359,size(SIMI_all,2));P_off_rem = zeros(256,359,size(SIMI_all,2));P_off_for = zeros(256,359,size(SIMI_all,2));
O_off_cue = zeros(256,359,size(SIMI_all,2));O_off_rem = zeros(256,359,size(SIMI_all,2));O_off_for = zeros(256,359,size(SIMI_all,2));

for s=1:size(SIMI_all,2)
    F_off_cue(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 3 & enc_cat_index(s,:) == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    F_off_rem(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 1 & enc_cat_index(s,:) == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    F_off_for(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 0 & enc_cat_index(s,:) == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    P_off_cue(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 3 & enc_cat_index(s,:) == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    P_off_rem(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 1 & enc_cat_index(s,:) == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    P_off_for(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 0 & enc_cat_index(s,:) == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    O_off_cue(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 3 & enc_cat_index(s,:) == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    O_off_rem(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 1 & enc_cat_index(s,:) == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    O_off_for(:,:,s) = mean(SIMI_all{s}(:,:,item_recall(s,:) == 0 & enc_cat_index(s,:) == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    
end 


% Stim order
stim1st_off_rem = zeros(256,359,size(SIMI_all,2));stim1st_off_for = zeros(256,359,size(SIMI_all,2));
stim2nd_off_rem = zeros(256,359,size(SIMI_all,2));stim2nd_off_for = zeros(256,359,size(SIMI_all,2));
stim3rd_off_rem = zeros(256,359,size(SIMI_all,2));stim3rd_off_for = zeros(256,359,size(SIMI_all,2));

for s=1:size(SIMI_all,2)
    stim1st_off_rem(:,:,s) = mean(SIMI_all{s}(:,:,recall_indperstim(s,:) == 2 & stim_ind == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    stim2nd_off_rem(:,:,s) = mean(SIMI_all{s}(:,:,recall_indperstim(s,:) == 2 & stim_ind == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    stim3rd_off_rem(:,:,s) = mean(SIMI_all{s}(:,:,recall_indperstim(s,:) == 2 & stim_ind == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    stim1st_off_for(:,:,s) = mean(SIMI_all{s}(:,:,recall_indperstim(s,:) ~= 2 & stim_ind == 1 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    stim2nd_off_for(:,:,s) = mean(SIMI_all{s}(:,:,recall_indperstim(s,:) ~= 2 & stim_ind == 2 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    stim3rd_off_for(:,:,s) = mean(SIMI_all{s}(:,:,recall_indperstim(s,:) ~= 2 & stim_ind == 3 & reject_index_img(s,:) == 1 & reject_index_off_perstim(s,:) == 1),3);
    
    
end 

%% plot
% Cluster
for fig=1
    figure(1);
    imagesc(Cluster==1);
    
    figure(2);
    imagesc(Cluster==2);
    
    figure(3);
    imagesc(Cluster==3);
end 
% SIMI
for fig=1
    figure(11);
    imagesc(mean(off_high_smo,3));
    xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
    xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
    yticks([0 0.5 1 1.5 2 ]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('High Mem')
    caxis([-0.08,0.08])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

%     figure(22);
%     imagesc(mean(off_n1,3));
%     xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
%     xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
%     yticks([0 0.5 1 1.5 2 ]*102.4)
%     yticklabels({'2000','1500','1000','500','0'})
%     xlabel('Offset [ms]')
%     ylabel('Encoding [ms]')
%     title('Recall 1')
%     caxis([-0.08,0.08])
%     set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    figure(33);
    imagesc(mean(off_low_smo,3));
    xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
    xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
    yticks([0 0.5 1 1.5 2 ]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Low Mem')
    caxis([-0.08,0.08])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    clearvars fig
    
    
    
%     figure(35);
%     imagesc(mean(off_low,3));
%     xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
%     xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
%     yticks([0 0.5 1 1.5 2 ]*102.4)
%     yticklabels({'2000','1500','1000','500','0'})
%     xlabel('Offset [ms]')
%     ylabel('Encoding [ms]')
%     title('Recall low')
%     caxis([-0.08,0.08])
%     set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
%     clearvars fig
    
    [~,~,~,t] = ttest(reshape(off_high_smo,[size(off_high_smo,1)*size(off_high_smo,2),size(off_high_smo,3)]),reshape(off_low_smo,[size(off_low_smo,1)*size(off_low_smo,2),size(off_low_smo,3)]),'Dim',2);
    

    [c_row1,c_col1] = find(Cluster==1);
    [c_row2,c_col2] = find(Cluster==2);
    [c_row3,c_col3] = find(Cluster==3);

    b1 = boundary(c_row1,c_col1);
    b2 = boundary(c_row2,c_col2);
    b3 = boundary(c_row3,c_col3);

    figure(37);
    imagesc(reshape(t.tstat,[size(off_low_smo,1),size(off_low_smo,2)]));hold on
    plot(c_col1(b1),c_row1(b1),'LineWidth',3,'color','k');
    plot(c_col2(b2),c_row2(b2),'LineWidth',3,'color','k');
    plot(c_col3(b3),c_row3(b3),'LineWidth',3,'color','k');hold off
    xticks([0 0.5 1 1.5 2 2.5 3 3.5]*102.4)
    xticklabels({'0','500','1000','1500','2000','2500','3000','3500'})
    yticks([0 0.5 1 1.5 2]*102.4)
    yticklabels({'2000','1500','1000','500','0'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Rem vs For (t)')
    caxis([-3,3])
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    clearvars fig

end 

% individual level
for i=1:28
    figure(41);
    subplot(5,6,i);
    imagesc(off_n0(:,:,i))
    caxis([-0.1,0.1])
    
    figure(42);
    subplot(5,6,i);
    imagesc(off_n2(:,:,i))
    caxis([-0.1,0.1])
    
    figure(43);
    subplot(5,6,i);
    imagesc(off_n2(:,:,i)-off_n0(:,:,i))
    caxis([-0.1,0.1])
end 


% categroy
for fig=2
Cax = [-0.08,0.08];
    figure(51);
    subplot(3,3,1)
    imagesc(mean(F_off_cue,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Face cue')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,2)
    imagesc(mean(F_off_rem,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Face rem')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,3)
    imagesc(mean(F_off_for,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Face for')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,4)
    imagesc(mean(P_off_cue,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Place cue')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,5)
    imagesc(mean(P_off_rem,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Place rem')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,6)
    imagesc(mean(P_off_for,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Place for')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,7)
    imagesc(mean(O_off_cue,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Object cue')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,8)
    imagesc(mean(O_off_rem,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Object rem')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,3,9)
    imagesc(mean(O_off_for,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('Object for')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);

    
end 


% stim seq
% categroy
for fig=2
Cax = [-0.08,0.08];
    figure(61);
    subplot(3,2,1)
    imagesc(mean(stim1st_off_rem,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('1st-off (Epi Rem)')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,2,2)
    imagesc(mean(stim1st_off_for,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('1st-off (Epi For)')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,2,3)
    imagesc(mean(stim2nd_off_rem,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('2nd-off (Epi Rem)')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,2,4)
    imagesc(mean(stim2nd_off_for,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('2nd-off (Epi For)')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    
    subplot(3,2,5)
    imagesc(mean(stim3rd_off_rem,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('3rd-off (Epi Rem)')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
    
    subplot(3,2,6)
    imagesc(mean(stim3rd_off_for,3));
    xticks([ 0.5  1.5 2.5  3.5]*102.4)
    xticklabels({'500','1500','2500','3500'})
    yticks([0  1  2 ]*102.4)
    yticklabels({'2500','1500','500'})
    xlabel('Offset [ms]')
    ylabel('Encoding [ms]')
    title('3rd-off (Epi For)')
    caxis(Cax)
    set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','normal', 'LineWidth', 2);
end 

%% LMM on the result
savepath = '';
low_pass = 0.5; % preprocessed data with two sets of parameters 0.1-30Hz and 0.5-30Hz
for loaddata=1


savename = strcat(savepath,'Beh&index_retordered_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%load ERPs

savename = strcat(savepath,'SIMI_all_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save ERPs

savename = strcat(savepath,'reject_ind_enc_',num2str(low_pass),'Hzlowpass_group.mat');
load (savename);%save rejection index


clearvars savepath savename loaddata
end 

% recall num for each stim
for s=1:size(SIMI_all,2)
    recall_indperstim(s,:) = rep_num(num_recall(s,:),3,1); 
    reject_index_off_perstim(s,:) = rep_num(reject_index_off(s,:),3,1); 
end
% stim index
stim_ind = rep_num([1,2,3],1,72); % first second and third img always;


for prepdata=1
    
for s = 1:30% same subjects used to identify the cluster and have coh data
    ind_trial = reject_index_img(s,:)==1 & reject_index_off_perstim(s,:)==1;
    data_2D = squeeze(SIMI_all{s}(:,:,ind_trial));
    
    data_all{s} = reshape(data_2D,[size(data_2D,1)*size(data_2D,2),size(data_2D,3)]);
    
end
Ydata_alltp = cat(2,data_all{:});

       
% other values

subject_all = zeros(0); img_seq = zeros(0);remfor=zeros(0);dif_rate=zeros(0);

for s = 1:30% same subjects used to identify the cluster and have coh data

    n_trial = sum(reject_index_img(s,:)==1 & reject_index_off_perstim(s,:)==1);
    
    % subject
    subject_all(length(subject_all)+1:length(subject_all)+n_trial) = s;
    
    % img_seq
    seq = rep_num(1:3,1,72);
    img_seq(length(img_seq)+1:length(img_seq)+n_trial) = seq(reject_index_img(s,:)==1 & reject_index_off_perstim(s,:)==1);
    
    % rem_for
    rf = recall_indperstim(s,:);
    remfor(length(remfor)+1:length(remfor)+n_trial) = rf(reject_index_img(s,:)==1 & reject_index_off_perstim(s,:)==1);
    
    % diff rating
    df = rep_num(coh_rating(s,:),3,1);
    dif_rate(length(dif_rate)+1:length(dif_rate)+n_trial) = df(reject_index_img(s,:)==1 & reject_index_off_perstim(s,:)==1);
end

% take out diff_rating nan trials
Ydata_alltp=Ydata_alltp(:,~isnan(dif_rate));
subject_all=subject_all(:,~isnan(dif_rate));
img_seq=img_seq(:,~isnan(dif_rate));
remfor=remfor(:,~isnan(dif_rate));
dif_rate=dif_rate(:,~isnan(dif_rate));

savename = strcat(savepath,'regression_SIMIencoff_p2p_notnorm_data.mat');
save (savename,'Ydata_alltp','subject_all','img_seq','remfor','dif_rate','-v7.3');
    
end 


%% LMM

%% True regression 
savename = strcat(savepath,'regression_SIMIencoff_p2p_notnorm_data.mat');
load (savename)

r_data = zeros(0);
for i_trial = 1:size(Ydata_alltp,2)
    disp([int2str(i_trial) '/5951'])
    data_ori = reshape(Ydata_alltp(:,i_trial),[256,359]);
    data = smooth2a(reshape(Ydata_alltp(:,i_trial),[256,359]),10);
    r_data(:,:,i_trial) = imresize(data,0.4,'bilinear');
    
end
  
t_map_real = zeros(4,size(r_data,3)); %1,intercept; 2, Seq   3, Remfor 4, diff rating
beta_map_real = zeros(4,size(r_data,3)); %1,intercept; 2, Seq   3, Remfor 4, diff rating
p_map_real = zeros(4,size(r_data,3)); %1,intercept; 2, Seq   3, Remfor 4, diff rating
Ydata_alltp_rs = reshape(r_data,[size(r_data,1)*size(r_data,2),size(Ydata_alltp,2)]);
% vif_result=zeros(2,2704);
tic
for i_time =1:size(Ydata_alltp_rs,1)
    if floor(i_time/size(r_data,1))==ceil(i_time/size(r_data,1))
        disp([int2str(i_time/size(Ydata_alltp_rs,1)*100) '% complete'])
    end
    tbl_time = table(Ydata_alltp_rs(i_time,:)',img_seq(:),subject_all(:),remfor(:),dif_rate(:),'VariableNames',{'Dis','Img_seq','Subject_ID','RemFor','DF_rate'});
    lme_time = fitlme(tbl_time,'Dis~Img_seq+RemFor+DF_rate+(1|Subject_ID)');
    
%     matrix_result = table2array(tbl_time);
%     vif_result(:,i_time) = vif(matrix_result(:,2:5));
    
    t_map_real(1,i_time) = lme_time.Coefficients{1,4};
    t_map_real(2,i_time) = lme_time.Coefficients{2,4};
    t_map_real(3,i_time) = lme_time.Coefficients{3,4};
    t_map_real(4,i_time) = lme_time.Coefficients{4,4};


    beta_map_real(1,i_time) = lme_time.Coefficients{1,2};
    beta_map_real(2,i_time) = lme_time.Coefficients{2,2};
    beta_map_real(3,i_time) = lme_time.Coefficients{3,2};
    beta_map_real(4,i_time) = lme_time.Coefficients{4,2};
    
    
    p_map_real(1,i_time) = lme_time.Coefficients{1,6};
    p_map_real(2,i_time) = lme_time.Coefficients{2,6};
    p_map_real(3,i_time) = lme_time.Coefficients{3,6};
    p_map_real(4,i_time) = lme_time.Coefficients{4,6};
    
    clearvars tbl_time lme_time
end
toc

savename = strcat(savepath,'regression_SIMIencoff_notnorm_0.4res_10ds.mat');

save (savename,'t_map_real','beta_map_real','p_map_real','-v7.3');



p_seq = reshape(p_map_real(2,:),[103,144]);
p_mem = reshape(p_map_real(3,:),[103,144]);
p_diff =reshape(p_map_real(4,:),[103,144]);


t_seq = reshape(t_map_real(2,:),[103,144]);
t_mem = reshape(t_map_real(3,:),[103,144]);
t_diff =reshape(t_map_real(4,:),[103,144]);


%% Regression with binary memory
savename = strcat(savepath,'regression_SIMIencoff_p2p_notnorm_data.mat');
load (savename)

r_data = zeros(0);
for i_trial = 1:size(Ydata_alltp,2)
    disp([int2str(i_trial) '/5951'])
    data_ori = reshape(Ydata_alltp(:,i_trial),[256,359]);
    data = smooth2a(reshape(Ydata_alltp(:,i_trial),[256,359]),10);
    r_data(:,:,i_trial) = imresize(data,0.4,'bilinear');
    
end
  
t_map_real = zeros(3,size(r_data,3)); %1,intercept; 2, Seq   3, Remfor 4, diff rating
beta_map_real = zeros(3,size(r_data,3)); %1,intercept; 2, Seq   3, Remfor 4, diff rating
p_map_real = zeros(3,size(r_data,3)); %1,intercept; 2, Seq   3, Remfor 4, diff rating
Ydata_alltp_rs = reshape(r_data,[size(r_data,1)*size(r_data,2),size(Ydata_alltp,2)]);
% vif_result=zeros(2,2704);
tic
for i_time =1:size(Ydata_alltp_rs,1)
    if floor(i_time/size(r_data,1))==ceil(i_time/size(r_data,1))
        disp([int2str(i_time/size(Ydata_alltp_rs,1)*100) '% complete'])
    end
    tbl_time = table(Ydata_alltp_rs(i_time,:)',img_seq(:),subject_all(:),remfor(:)==2,'VariableNames',{'Dis','Img_seq','Subject_ID','RemFor'});
    lme_time = fitlme(tbl_time,'Dis~Img_seq+RemFor+(1|Subject_ID)');
    
%     matrix_result = table2array(tbl_time);
%     vif_result(:,i_time) = vif(matrix_result(:,2:5));
    
    t_map_real(1,i_time) = lme_time.Coefficients{1,4};
    t_map_real(2,i_time) = lme_time.Coefficients{2,4};
    t_map_real(3,i_time) = lme_time.Coefficients{3,4};


    beta_map_real(1,i_time) = lme_time.Coefficients{1,2};
    beta_map_real(2,i_time) = lme_time.Coefficients{2,2};
    beta_map_real(3,i_time) = lme_time.Coefficients{3,2};

    
    
    p_map_real(1,i_time) = lme_time.Coefficients{1,6};
    p_map_real(2,i_time) = lme_time.Coefficients{2,6};
    p_map_real(3,i_time) = lme_time.Coefficients{3,6};

    
    clearvars tbl_time lme_time
end
toc


savename = strcat(savepath,'regression_SIMIencoff_3factor_bMem_notnorm_0.4res_10ds.mat'); % three factor and binary memory

save (savename,'t_map_real','beta_map_real','p_map_real','-v7.3');



p_seq = reshape(p_map_real(2,:),[103,144]);
p_mem = reshape(p_map_real(3,:),[103,144]);



t_seq = reshape(t_map_real(2,:),[103,144]);
t_mem = reshape(t_map_real(3,:),[103,144]);
