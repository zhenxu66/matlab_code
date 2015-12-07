clear all;
close all;
load('UTA all data.mat')
%Compare acceleration period and constant running period, fast training and
%with smaller error find with logic
%WATSON_CHRISTA Subject 1 F
%DAWOOD_WAJID Subject 2 rest 31:69 running  31:467
%AGUIRRE_SERGIO Subject 3  rest 31:51 running  31:738
%HULLA_RYAN Subject 4
%JENKINS_ERIC Subject 5 time 1549

%%
sig1_t=cell2mat(sig1(1));
sig2_t=cell2mat(sig2(1));
sig3_t=cell2mat(sig3(1));
sig4_t=cell2mat(sig4(1));
sig5_t=cell2mat(sig5(1));


sig1HR_F=cell2mat(sig1(2))*60;
sig2HR_F=cell2mat(sig2(2))*60;
sig3HR_F=cell2mat(sig3(2))*60;
sig4HR_F=cell2mat(sig4(2))*60;
sig5HR_F=cell2mat(sig5(2))*60;

% sig1HR_F=cell2mat(sig1(2));
% sig2HR_F=cell2mat(sig2(2));
% sig3HR_F=cell2mat(sig3(2));
% sig4HR_F=cell2mat(sig4(2));
% sig5HR_F=cell2mat(sig5(2));

sig1EE_F=cell2mat(sig1(3));
sig2EE_F=cell2mat(sig2(3));
sig3EE_F=cell2mat(sig3(3));
sig4EE_F=cell2mat(sig4(3));
sig5EE_F=cell2mat(sig5(3));

% sig1HR=cell2mat(sig1(4));
% sig2HR=cell2mat(sig2(4));
% sig3HR=cell2mat(sig3(4));
% sig4HR=cell2mat(sig4(4));
% sig5HR=cell2mat(sig5(4));

sig1HR=cell2mat(sig1(4))*60;
sig2HR=cell2mat(sig2(4))*60;
sig3HR=cell2mat(sig3(4))*60;
sig4HR=cell2mat(sig4(4))*60;
sig5HR=cell2mat(sig5(4))*60;

sig1EE=cell2mat(sig1(5));
sig2EE=cell2mat(sig2(5));
sig3EE=cell2mat(sig3(5));
sig4EE=cell2mat(sig4(5));
sig5EE=cell2mat(sig5(5));

%All transition
HR_total_p_2=[sig2HR(1:7);sig2HR(122:127);sig2HR(260:280);sig2HR(417:end)];
HR_total_p_3=[sig3HR(162:167);sig3HR(250:270);sig3HR(671:end)];
HR_total_p_4=[sig4HR(1:7);sig4HR(144:149);sig4HR(325:345);sig4HR(550:end)];
HR_total_p_5=[sig5HR(1:7);sig5HR(120:125);sig5HR(281:301);sig5HR(492:end)];
HR_total_p=[HR_total_p_2;HR_total_p_3;HR_total_p_4;HR_total_p_5];

EE_total_p_2=[sig2EE(1:7);sig2EE(122:127);sig2EE(260:280);sig2EE(417:end)];
EE_total_p_3=[sig3EE(162:167);sig3EE(250:270);sig3EE(671:end)];
EE_total_p_4=[sig4EE(1:7);sig4EE(144:149);sig4EE(325:345);sig4EE(550:end)];
EE_total_p_5=[sig5EE(1:7);sig5EE(120:125);sig5EE(281:301);sig5EE(492:end)];
EE_total_p=[EE_total_p_2;EE_total_p_3;EE_total_p_4;EE_total_p_5];
%%
% for S=1:4
%     if S==1
%         range=[1,2,3,119];
%         speed='0-3mph';
%     elseif S==2 
%         range=[120,121,122,284];
%         speed='3-4mph';
%     elseif S==3  
%         range=[284,286,287,494];
%         speed='4-5mph';
%     else 
%         range=[495,length(sig5HR),287,491];
% %        range=[495,501,287,491];
%         speed='5mph-rest';
%     end
%         
%         
%     HR_5_p2=sig5HR(range(3):range(4));
%     HR_5_p1=sig5HR(range(1):range(2));
%     
%     EE_5_p2=sig5EE(range(3):range(4));
%     EE_5_p1=sig5EE(range(1):range(2));
%     
%     Window=10;
%     
%     for j=1:floor(length(HR_5_p2)/Window),
%         k=7+Window*(j-1);
%         reg5_2=LinearModel.fit(sig5HR(k:k+Window),sig5EE((k:k+Window)));
%         RMSE5_2(j)=reg5_2.RMSE;
%         Slope5_2(j)=double(reg5_2.Coefficients(1,1));
%         Intercept5_2(j)=double(reg5_2.Coefficients(2,1));
%     end;
%     
%     reg5_1=LinearModel.fit(HR_5_p1,EE_5_p1);
%     RMSE5_1=reg5_2.RMSE;
%     Slope5_1=double(reg5_1.Coefficients(1,1));
%     Intercept5_1=double(reg5_1.Coefficients(2,1));
%     
%     reg5=LinearModel.fit((sig5HR),sig5EE);
%     RMSE5=reg5.RMSE;
%     Slope5=double(reg5.Coefficients(1,1));
%     Intercept5=double(reg5.Coefficients(2,1));
%     
%     Array1=ones(1,floor(length(HR_5_p2)/Window));
%     
%     figure;
%     subplot(311);
%     plot(1:floor(length(HR_5_p2)/Window),Slope5_2,'-.s',1:floor(length(HR_5_p2)/Window),Array1*Slope5,'--mo',1:floor(length(HR_5_p2)/Window),Array1*Slope5_1,'-ch');
%     title('5th People Analysis (sliced into same size time windows)','FontSize',20);set(gca,'FontSize',16);
%     xlabel('Time slot','FontSize',16); ylabel('Slope','FontSize',16);
%     hleg3 = legend(['Constant Running',speed],'Whole range 3-5mph',['Speed Changing',speed]);
%     set(hleg3,'Location','NorthWest');
%     
%     subplot(312);
%     plot(1:floor(length(HR_5_p2)/Window),Intercept5_2,'-.s',1:floor(length(HR_5_p2)/Window),Array1*Intercept5,'--mo',1:floor(length(HR_5_p2)/Window),Array1*Intercept5_1,'-ch');
%     xlabel('Time slot','FontSize',16); ylabel('Intercept','FontSize',16);
%     hleg3 = legend(['Constant Running',speed],'Whole range 3-5mph',['Speed Changing',speed]);
%     set(hleg3,'Location','NorthWest');
%     
%     subplot(313);
%     plot(1:floor(length(HR_5_p2)/Window),RMSE5_2,'-.s',1:floor(length(HR_5_p2)/Window),Array1*RMSE5,'--mo',1:floor(length(HR_5_p2)/Window),Array1*RMSE5_1,'-ch');
%     xlabel('Time slot','FontSize',16); ylabel('RMSE','FontSize',16);
%     hleg3 = legend(['Constant Running',speed],'Whole range 3-5mph',['Speed Changing',speed]);
%     set(hleg3,'Location','NorthWest');
% end;
%%
%Smooth the data
sig1HR_F_S = sgolayfilt(sig1HR_F,2,13);
sig2HR_F_S = sgolayfilt(sig2HR_F,2,13);
sig3HR_F_S = sgolayfilt(sig3HR_F,2,13);
sig4HR_F_S = sgolayfilt(sig4HR_F,2,13);
sig5HR_F_S = sgolayfilt(sig5HR_F,2,13);

sig1HR_S = sgolayfilt(sig1HR,2,13);
sig2HR_S = sgolayfilt(sig2HR,2,13);
sig3HR_S = sgolayfilt(sig3HR,2,13);
sig4HR_S = sgolayfilt(sig4HR,2,13);
sig5HR_S = sgolayfilt(sig5HR,2,13);

% figure;
% plot(sig1_t,sig1HR_F,'-.s',sig1_t,sig1HR_F_S,'--mo');
% figure;
% plot(sig2_t,sig2HR_F,'-.s',sig2_t,sig2HR_F_S,'--mo');
% figure;
% plot(sig3_t,sig3HR_F,'-.s',sig3_t,sig3HR_F_S,'--mo');
% figure;
% plot(sig4_t,sig4HR_F,'-.s',sig4_t,sig4HR_F_S,'--mo');
% figure;
% plot(sig5_t,sig5HR_F,'-.s',sig5_t,sig5HR_F_S,'--mo');

%%  5

for j=1:size(sig5HR_F_S)-3,
    l1=((sig5HR_F_S(j+1)-sig5HR_F_S(j))/(sig5_t(j+2)-sig5_t(j+1))>1);
    l2=((sig5HR_F_S(j+2)-sig5HR_F_S(j+1))/(sig5_t(j+3)-sig5_t(j+2))>1);
    if l1&&l2
        Selected_5(j)=j;
    end;
end;

Selected_5(Selected_5==0) = [];

for j=1:3
    Selected_5_2(3*j-2)=Selected_5(j);
    Selected_5_2(3*j-1)=Selected_5(j)+1;
    Selected_5_2(3*j)=Selected_5(j)+2;
end;

Selected_5_2=unique(Selected_5_2);
Selected_5_2=Selected_5_2(1:3);% Point 4 of HR is defect keep same HR not take

% reg5_s=LinearModel.fit(sig5HR(Selected_5_2),sig5EE(Selected_5_2));
reg5_s=LinearModel.fit(sig5HR_F_S(Selected_5_2+1),sig5EE(Selected_5_2)); 
% must smooth to get clean Intercept and slope
%RMSE5_s=reg5_s.RMSE;
Intercept5_s=double(reg5_s.Coefficients(1,1));
Slope5_s=double(reg5_s.Coefficients(2,1));
EstiEE5=Intercept5_s+Slope5_s*sig5HR;
RMSE5_s = sqrt( sum( (sig5EE(:)-EstiEE5(:)).^2) / length(sig5EE) );

sig5EE_T_esti(1)=0;
for j=1:size(EstiEE5),
    sig5EE_T_esti(j+1)= sig5EE_T_esti(j)+(sig5_t(j+1)-sig5_t(j))*EstiEE5(j);
end

reg5=LinearModel.fit(sig5HR,sig5EE);
RMSE5=reg5.RMSE;
Intercept5=double(reg5.Coefficients(1,1));
Slope5=double(reg5.Coefficients(2,1));

%% 4

for j=1:size(sig4HR)-3,
    l1=((sig4HR(j+1)-sig4HR(j))/(sig4_t(j+2)-sig4_t(j+1))>1);
    l2=((sig4HR(j+2)-sig4HR(j+1))/(sig4_t(j+3)-sig4_t(j+2))>1);
    if l1&&l2
        Selected_4(j)=j;
    end;
end;

Selected_4(Selected_4==0) = [];
for j=1:2
    Selected_4_2(3*j-2)=Selected_4(j);
    Selected_4_2(3*j-1)=Selected_4(j)+1;
    Selected_4_2(3*j)=Selected_4(j)+2;
end;

Selected_4_2=unique(Selected_4_2);
Selected_4_2=Selected_4_2(1:6);% Add more Points as HR is too smooth 
% Selected_4_2=[1;2;3;143;144;145];

reg4_s=LinearModel.fit(sig4HR(Selected_4_2),sig4EE(Selected_4_2));
% reg4_s=LinearModel.fit(sig4HR_F_S(Selected_4_2+1),sig5EE(Selected_4_2)); 
%Smooth of 4 is not good, remove lots of info.
% RMSE4_s=reg4_s.RMSE;
Intercept4_s=double(reg4_s.Coefficients(1,1));
Slope4_s=double(reg4_s.Coefficients(2,1));
EstiEE4=Intercept4_s+Slope4_s*sig4HR;
RMSE4_s = sqrt( sum( (sig4EE(:)-EstiEE4(:)).^2) / length(sig4EE) );

sig4EE_T_esti(1)=0;
for j=1:size(EstiEE4),
    sig4EE_T_esti(j+1)= sig4EE_T_esti(j)+(sig4_t(j+1)-sig4_t(j))*EstiEE4(j);
end

reg4=LinearModel.fit(sig4HR,sig4EE);
RMSE4=reg4.RMSE;
Intercept4=double(reg4.Coefficients(1,1));
Slope4=double(reg4.Coefficients(2,1));

%% 3

for j=1:size(sig3HR)-3,
    l1=((sig3HR(j+1)-sig3HR(j))/(sig3_t(j+2)-sig3_t(j+1))>0.6);
    l2=((sig3HR(j+2)-sig3HR(j+1))/(sig3_t(j+3)-sig3_t(j+2))>0.6);
    if l1&&l2
        Selected_3(j)=j;
    end;
end;

Selected_3(Selected_3==0) = [];
for j=1:2
    Selected_3_2(3*j-2)=Selected_3(j);
    Selected_3_2(3*j-1)=Selected_3(j)+1;
    Selected_3_2(3*j)=Selected_3(j)+2;
end;

Selected_3_2=unique(Selected_3_2);
Selected_3_2=Selected_3_2(1:6);% Add more Points as HR is too smooth 
% Selected_3_2=(163:169);

reg3_s=LinearModel.fit(sig3HR(Selected_3_2),sig3EE(Selected_3_2));
% RMSE3_s=reg3_s.RMSE;
Intercept3_s=double(reg3_s.Coefficients(1,1));
Slope3_s=double(reg3_s.Coefficients(2,1));
EstiEE3=Intercept3_s+Slope3_s*sig3HR;
RMSE3_s = sqrt( sum( (sig3EE(:)-EstiEE3(:)).^2) / length(sig3EE) );



sig3EE_T_esti(1)=0;
for j=1:size(EstiEE3),
    sig3EE_T_esti(j+1)= sig3EE_T_esti(j)+(sig3_t(j+1)-sig3_t(j))*EstiEE3(j);
end

reg3=LinearModel.fit(sig3HR,sig3EE);
RMSE3=reg3.RMSE;
Intercept3=double(reg3.Coefficients(1,1));
Slope3=double(reg3.Coefficients(2,1));

%% 2

for j=1:size(sig2HR)-2,
    l1=((sig2HR(j+1)-sig2HR(j))/(sig2_t(j+2)-sig2_t(j+1))>1);
    l2=((sig2HR(j+2)-sig2HR(j+1))/(sig2_t(j+2)-sig2_t(j+2))>1);
    if l1&&l2
        Selected_2(j)=j;
    end;
end;

Selected_2(Selected_2==0) = [];
for j=1:2
    Selected_2_2(3*j-2)=Selected_2(j);
    Selected_2_2(3*j-1)=Selected_2(j)+1;
    Selected_2_2(3*j)=Selected_2(j)+2;
end;

Selected_2_2=unique(Selected_2_2);
Selected_2_2=Selected_2_2(1:6);% Add more Points as HR is too much artifact
% Selected_2_2=[29;30;31];
% Selected_2_2=[121;122;123];

reg2_s=LinearModel.fit(sig2HR(Selected_2_2),sig2EE(Selected_2_2));
% RMSE2_s=reg2_s.RMSE;
Intercept2_s=double(reg2_s.Coefficients(1,1));
Slope2_s=double(reg2_s.Coefficients(2,1));

EstiEE2=Intercept2_s+Slope2_s*sig2HR;
RMSE2_s = sqrt( sum( (sig2EE(:)-EstiEE2(:)).^2) / length(sig2EE) );

sig2EE_T_esti(1)=0;
for j=1:size(EstiEE2),
    sig2EE_T_esti(j+1)= sig2EE_T_esti(j)+(sig2_t(j+1)-sig2_t(j))*EstiEE2(j);
end

reg2=LinearModel.fit(sig2HR,sig2EE);
RMSE2=reg2.RMSE;
Intercept2=double(reg2.Coefficients(1,1));
Slope2=double(reg2.Coefficients(2,1));

%% 1

for j=1:size(sig1HR)-2,
    l1=((sig1HR(j+1)-sig1HR(j))/(sig1_t(j+2)-sig1_t(j+1))>1);
    l2=((sig1HR(j+2)-sig1HR(j+1))/(sig1_t(j+2)-sig1_t(j+2))>1);
    if l1&&l2
        Selected_1(j)=j;
    end;
end;

Selected_1(Selected_1==0) = [];
for j=1:1
    Selected_1_2(3*j-2)=Selected_1(j);
    Selected_1_2(3*j-1)=Selected_1(j)+1;
    Selected_1_2(3*j)=Selected_1(j)+2;
end;

% Selected_1_2=[5,6,7];
reg1_s=LinearModel.fit(sig1HR_F_S(Selected_1_2),sig1EE(Selected_1_2),'RobustOpts','on');
% RMSE1_s=reg1_s.RMSE;
Intercept1_s=double(reg1_s.Coefficients(1,1));
Slope1_s=double(reg1_s.Coefficients(2,1));

% EstiEE1=Intercept2_s+Slope2_s*sig1HR_F_S(2:end);
EstiEE1=Intercept2_s+Slope2_s*sig1HR_F_S(2:end);
RMSE1_s = sqrt( sum( (sig1EE(:)-EstiEE1(:)).^2) / length(sig1EE) );

sig1EE_T_esti(1)=0;
for j=1:size(EstiEE1),
    sig1EE_T_esti(j+1)= sig1EE_T_esti(j)+(sig1_t(j+1)-sig1_t(j))*EstiEE1(j);
end

reg1=LinearModel.fit(sig1HR,sig1EE);
RMSE1=reg1.RMSE;
Intercept1=double(reg1.Coefficients(1,1));
Slope1=double(reg1.Coefficients(2,1));





%% 2345
sig2345HR=[sig2HR;sig3HR;sig4HR;sig5HR];
sig2345EE=[sig2EE;sig3EE;sig4EE;sig5EE];
sig2345_t=[sig2_t;sig3_t;sig4_t;sig5_t];

for j=1:size(sig2345HR)-3,
    l1=((sig2345HR(j+1)-sig2345HR(j))/(sig2345_t(j+2)-sig2345_t(j+1))>0.04);
    l2=((sig2345HR(j+2)-sig2345HR(j+1))/(sig2345_t(j+3)-sig2345_t(j+2))>0.04);
    if l1&&l2
        Selected_2345(j)=j;
    end;
end;

Selected_2345(Selected_2345==0) = [];
for j=1:length(Selected_2345)
    Selected_2345_2(3*j-2)=Selected_2345(j);
    Selected_2345_2(3*j-1)=Selected_2345(j)+1;
    Selected_2345_2(3*j)=Selected_2345(j)+2;
end;

reg2345_s=LinearModel.fit(sig2345HR(Selected_2345_2),sig2345EE(Selected_2345_2));
RMSE2345_s=reg2345_s.RMSE;
Intercept2345_s=double(reg2345_s.Coefficients(1,1));
Slope2345_s=double(reg2345_s.Coefficients(2,1));

reg2345=LinearModel.fit(sig2345HR,sig2345EE);
RMSE2345=reg2345.RMSE;
Intercept2345=double(reg2345.Coefficients(1,1));
Slope2345=double(reg2345.Coefficients(2,1));

sig1EE_T_err=sig1EE_T_esti'-sig1EE_F;
sig2EE_T_err=sig2EE_T_esti'-sig2EE_F;
sig3EE_T_err=sig3EE_T_esti'-sig3EE_F;
sig4EE_T_err=sig4EE_T_esti'-sig4EE_F;
sig5EE_T_err=sig5EE_T_esti'-sig5EE_F;



% figure;
% subplot(211);
% plot(sig1_t,sig1EE_F,'--go',sig1_t,sig1EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE ','FontSize',20);set(gca,'FontSize',16);
% hleg1 = legend('Test data','Regression from selected data');
% set(hleg1,'Location','NorthWest')
% set(hleg1,'Interpreter','none')
% 
% subplot(212);
% plot(sig1_t(2:end),sig1EE,'--go',sig1_t(2:end),EstiEE1,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE Kcal/Min','FontSize',20);set(gca,'FontSize',16);
% hleg2 = legend(['Test','*RMSE=',num2str(RMSE1)],['Regression from selected points','*RMSE=',num2str(RMSE1_s)]);
% set(hleg2,'Location','NorthWest')
% set(hleg2,'Interpreter','none')
% 
% figure;
% subplot(211);
% plot(sig2_t,sig2EE_F,'--go',sig2_t,sig2EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE ','FontSize',20);set(gca,'FontSize',16);
% hleg1 = legend('Test data','Regression from selected data');
% set(hleg1,'Location','NorthWest');
% set(hleg1,'Interpreter','none');
% 
% subplot(212);
% plot(sig2_t(2:end),sig2EE,'--go',sig2_t(2:end),EstiEE2,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
% title('EE Kcal/Min','FontSize',20);set(gca,'FontSize',16);
% hleg2 = legend(['Test','*RMSE=',num2str(RMSE2)],['Regression from selected points','*RMSE=',num2str(RMSE2_s)]);
% set(hleg2,'Location','NorthWest');
% set(hleg2,'Interpreter','none');
% 
% figure;
% subplot(211);
% plot(sig3_t,sig3EE_F,'--go',sig3_t,sig3EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE ','FontSize',20);set(gca,'FontSize',16);
% hleg1 = legend('Test data','Regression from selected data');
% set(hleg1,'Location','NorthWest');
% set(hleg1,'Interpreter','none');
% 
% subplot(212);
% plot(sig3_t(2:end),sig3EE,'--go',sig3_t(2:end),EstiEE3,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
% title('EE Kcal/Min','FontSize',20);set(gca,'FontSize',16);
% hleg2 = legend(['Test','*RMSE=',num2str(RMSE3)],['Regression from selected points','*RMSE=',num2str(RMSE3_s)]);
% set(hleg2,'Location','NorthWest');
% set(hleg2,'Interpreter','none');
% 
% figure;
% subplot(211);
% plot(sig4_t,sig4EE_F,'--go',sig4_t,sig4EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE ','FontSize',20);set(gca,'FontSize',16);
% hleg1 = legend('Test data','Regression from selected data');
% set(hleg1,'Location','NorthWest');
% set(hleg1,'Interpreter','none');
% 
% subplot(212);
% plot(sig4_t(2:end),sig4EE,'--go',sig4_t(2:end),EstiEE4,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
% title('EE Kcal/Min','FontSize',20);set(gca,'FontSize',16);
% hleg2 = legend(['Test','*RMSE=',num2str(RMSE4)],['Regression from selected points','*RMSE=',num2str(RMSE4_s)]);
% set(hleg2,'Location','NorthWest');
% set(hleg2,'Interpreter','none');
% 
% figure;
% subplot(211);
% plot(sig5_t,sig5EE_F,'--go',sig5_t,sig5EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE ','FontSize',20);set(gca,'FontSize',16);
% hleg1 = legend('Test data','Regression from selected data');
% set(hleg1,'Location','NorthWest');
% set(hleg1,'Interpreter','none');
% 
% subplot(212);
% plot(sig5_t(2:end),sig5EE,'--go',sig5_t(2:end),EstiEE5,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
% title('EE Kcal/Min','FontSize',20);set(gca,'FontSize',16);
% hleg2 = legend(['Test','*RMSE=',num2str(RMSE5)],['Regression from selected points','*RMSE=',num2str(RMSE5_s)]);
% set(hleg2,'Location','NorthWest');
% set(hleg2,'Interpreter','none');


% figure;
% subplot(211);
% plot(sig1_t,sig1EE_F,'--go',sig1_t,sig1EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE ','FontSize',20);set(gca,'FontSize',16);
% hleg1 = legend('Test data','Regression from selected data');
% set(hleg1,'Location','NorthWest')
% set(hleg1,'Interpreter','none')
% 
% subplot(212);
% plot(sig1_t(2:end),sig1EE,'--go',sig1_t(2:end),EstiEE1,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title('EE Kcal/Min','FontSize',20);set(gca,'FontSize',16);
% hleg2 = legend(['Test','*RMSE=',num2str(RMSE1)],['Regression from selected points','*RMSE=',num2str(RMSE1_s)]);
% set(hleg2,'Location','NorthWest')
% set(hleg2,'Interpreter','none')
%% 
%2 1,5  3 2,3,4
figure;
set(gcf,'color','w');
% subplot(221);
% plot(sig1_t,sig1EE_F,'--go',sig1_t,sig1EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title(['No.1 EE from',mat2str(Selected_1_2)],'FontSize',20);set(gca,'FontSize',10);
% hleg1 = legend('1st Test data','Reg selected');
% set(hleg1,'Location','NorthWest');
% set(hleg1,'Interpreter','none');
% 
% subplot(222);
% plot(sig1_t(2:end),sig1EE,'--go',sig1_t(2:end),EstiEE1,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
% hleg2 = legend(['*RMSE(T)=',num2str(RMSE1)],['*RMSE(S)=',num2str(RMSE1_s)]);
% set(hleg1,'Location','NorthWest');
% set(hleg1,'Interpreter','none');
% 
% subplot(223);
% plot(sig5_t,sig5EE_F,'--go',sig5_t,sig5EE_T_esti,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
% title(['No.5 EE from',mat2str(Selected_5_2)],'FontSize',20);set(gca,'FontSize',10);
% hleg1 = legend('Test data','Reg selected');
% set(hleg1,'Location','NorthWest');
% set(hleg1,'Interpreter','none');
% 
% subplot(224);
% plot(sig5_t(2:end),sig5EE,'--go',sig5_t(2:end),EstiEE5,':r*');
% xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
% hleg2 = legend(['Test','*RMSE(T)=',num2str(RMSE5)],['RMSE(S)=',num2str(RMSE5_s)]);
% set(hleg2,'Location','NorthWest');
% set(hleg2,'Interpreter','none');

subplot(231);
plot(sig1_t,sig1EE_F,'--go',sig1_t,sig1EE_T_esti,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
title(['No.1 EE from',mat2str(Selected_1_2)],'FontSize',20);set(gca,'FontSize',10);
hleg1 = legend('1st Test data','Reg selected');
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');
set(gca,'FontName','Times New Roman','FontSize',19)

subplot(232);
plot(sig1_t(2:end),sig1EE,'--go',sig1_t(2:end),EstiEE1,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
hleg2 = legend(['*RMSE(T)=',num2str(RMSE1)],['*RMSE(S)=',num2str(RMSE1_s)]);
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

subplot(233);
plot(sig1_t(2:end),sig1HR,'--go',sig1_t(Selected_1_2+1),sig1HR(Selected_1_2),'r*');
xlabel('Time (second)','FontSize',16); ylabel('HR (BPM)','FontSize',16);
hleg2 = legend(['HR'],['HR',mat2str(Selected_1_2)]);
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

subplot(234);
plot(sig5_t,sig5EE_F,'--go',sig5_t,sig5EE_T_esti,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
title(['No.5 EE from',mat2str(Selected_5_2)],'FontSize',20);set(gca,'FontSize',10);
hleg1 = legend('Test data','Reg selected');
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

subplot(235);
plot(sig5_t(2:end),sig5EE,'--go',sig5_t(2:end),EstiEE5,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
hleg2 = legend(['*RMSE(T)=',num2str(RMSE5)],['RMSE(S)=',num2str(RMSE5_s)]);
set(hleg2,'Location','NorthWest');
set(hleg2,'Interpreter','none');

subplot(236);
plot(sig5_t(2:end),sig5HR,'--go',sig5_t(Selected_5_2+1),sig5HR(Selected_5_2),'r*');
xlabel('Time (second)','FontSize',16); ylabel('HR (BPM)','FontSize',16);
hleg2 = legend(['HR'],['HR',mat2str(Selected_5_2)]);
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

figure;
set(gcf,'color','w');
subplot(331);
plot(sig2_t,sig2EE_F,'--go',sig2_t,sig2EE_T_esti,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
title(['No.2 EE from',mat2str(Selected_2_2)],'FontSize',20);set(gca,'FontSize',10);
hleg1 = legend('Test data','Reg selected');
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

subplot(332);
plot(sig2_t(2:end),sig2EE,'--go',sig2_t(2:end),EstiEE2,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
hleg2 = legend(['RMSE(T)=',num2str(RMSE2)],['RMSE(S)=',num2str(RMSE2_s)]);
set(hleg2,'Location','NorthWest');
set(hleg2,'Interpreter','none');

subplot(333);
plot(sig2_t(2:end),sig2HR,'--go',sig2_t(Selected_2_2+1),sig2HR(Selected_2_2),'r*');
xlabel('Time (second)','FontSize',16); ylabel('HR (BPM)','FontSize',16);
hleg2 = legend(['HR'],['HR',mat2str(Selected_2_2)]);
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');


subplot(334);
plot(sig3_t,sig3EE_F,'--go',sig3_t,sig3EE_T_esti,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
title(['No.3 EE from',mat2str(Selected_3_2)],'FontSize',20);set(gca,'FontSize',10);
hleg1 = legend('Test data','Reg selected');
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

subplot(335);
plot(sig3_t(2:end),sig3EE,'--go',sig3_t(2:end),EstiEE3,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
hleg2 = legend(['RMSE(T)=',num2str(RMSE3)],['RMSE(S)=',num2str(RMSE3_s)]);
set(hleg2,'Location','NorthWest');
set(hleg2,'Interpreter','none');

subplot(336);
plot(sig3_t(2:end),sig3HR,'--go',sig3_t(Selected_3_2+1),sig3HR(Selected_3_2),'r*');
xlabel('Time (second)','FontSize',16); ylabel('HR (BPM)','FontSize',16);
hleg2 = legend(['HR'],['HR',mat2str(Selected_3_2)]);
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');


subplot(337);
plot(sig4_t,sig4EE_F,'--go',sig4_t,sig4EE_T_esti,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal)','FontSize',16);
title(['No.4 EE from',mat2str(Selected_4_2)],'FontSize',20);set(gca,'FontSize',10);
hleg1 = legend('Test data','Reg selected');
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

subplot(338);
plot(sig4_t(2:end),sig4EE,'--go',sig4_t(2:end),EstiEE4,':r*');
xlabel('Time (second)','FontSize',16); ylabel('EE (Kcal/Min)','FontSize',16);
% title('EE Kcal/Min','FontSize',20);set(gca,'FontSize',10);
hleg2 = legend(['RMSE(T)=',num2str(RMSE4)],['RMSE(S)=',num2str(RMSE4_s)]);
set(hleg2,'Location','NorthWest');
set(hleg2,'Interpreter','none');

subplot(339);
plot(sig4_t(2:end),sig4HR,'--go',sig4_t(Selected_4_2+1),sig4HR(Selected_4_2),'r*');
xlabel('Time (second)','FontSize',16); ylabel('HR (BPM)','FontSize',16);
hleg2 = legend(['HR'],['HR',mat2str(Selected_4_2)]);
set(hleg1,'Location','NorthWest');
set(hleg1,'Interpreter','none');

% 
% 
% 
% 
% 
