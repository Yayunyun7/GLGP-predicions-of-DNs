
% data_ppg: ppg signal
% tppg: time span of ppg signal (seconds)

% PPG R peaks detection and DNs detection 
[R, Q]= PPG_peakdetection(data_ppg, 200);
[aloc, bloc, eloc, nloc, p2loc, NonDN, DN, DNloc] = PPG_abcde_detection2(data_ppg, R, Q, 200); 

% Synchorsqueezing transform of PPG signal
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(data_ppg',0, 0.1, 0.0002, 10, 2001, 1, 6, 1, 1, 0);

% phase extraction and reconstruct fundamental (F) and harmonic components (M1, M2, M3)
idx = find(tfrsqtic*200 < 3) ;     
idx = idx(2:end);
[c1] = CurveExt_M(abs(tfrsq(idx,:)'), 0.2) ;
[h,Dh,tt] = hermf(2001,1,6);
[f1] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 20, c1, 0.06, h(200*5+1)) ;
[f2] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 20, c1*2, 0.06, h(200*5+1)) ;
[f3] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 20, c1*3, 0.06, h(200*5+1)) ;
[f4] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 20, c1*4, 0.06, h(200*5+1)) ;

t0 = [1/200:1/200*10:tppg]; 
t1 = linspace(0, tppg ,200*tppg);
F = interp1(t0,f1,t1,"spline");
M1 = interp1(t0,f2,t1,"spline");
M2 = interp1(t0,f3,t1,"spline");
M3 = interp1(t0,f4,t1,"spline");

%%
% Figure 1: plot DNs on PPG
%   subplot1: PPG with DNs predicted by APPG; 
%   subplot2: PPG with both DNs predicted by APPG and true DNs
figure(1);
subplot(2,1,1);
h(1) = plot(t1, data_ppg, "k"); hold on; grid on;
h(2) = scatter(t1(eloc),data_ppg(eloc),30, "filled", "r" );
xlim([6645 6660]);
xlabel("Time(seconds)", 'FontSize',18);
ylabel("PPG(au)", 'FontSize',18);
legend(h([2]),"DN by APPG", 'FontSize',16);

subplot(2,1,2);
plot(t1, data_ppg, "k"); hold on; grid on;
q(1)= scatter(t1(eloc),data_ppg(eloc),30, "filled", "r" );
q(2)= scatter(t1(DNloc),data_ppg(DNloc),30,[0.1000 0.6000 0.3000], '^',"filled");
xlim([6645 6660]);
xlabel("Time(seconds)", 'FontSize',18);
ylabel("PPG(au)", 'FontSize',18);
legend(q([1,2]),"DN by APPG", "True DN", 'FontSize',16);


%% Apply GLGP model to predict DNs for beats without DNs
%  All beats are sperated into two groups
%   G1: beats with true DN 
%   G2: beats without DNs

% Denote the location of R peak, DNs for beats G1 and G2
r1 = R(DN)';
d1 = DNloc;
r2 = R(NonDN)';

% Collect features and labels in G1 as training dataset (W1)
%   Dk: instant phase change of fundamental and harmonic components at R peaks in G1(radians)
%   AK: instant amplitude of fundamental and harmonic components at R peaks  in G1
%   Y1: time difference between R peaks and true DNs (seconds)

Dk = zeros(4, length(r1)); 
for i = 1:length(r1)
    Dk(1,i) = R_phase(r1(i), r1(i)-1, F);
    Dk(2,i) = R_phase(r1(i), r1(i)-1, M1);
    Dk(3,i) = R_phase(r1(i), r1(i)-1, M2);
    Dk(4,i) = R_phase(r1(i), r1(i)-1, M3);
end  

Ak = zeros(4, length(r1));
for i = 1:length (r1)
    Ak(1,i)= abs(F(r1(i))-F(r1(i)-1));
    Ak(2,i)= abs(M1(r1(i))-M1(r1(i)-1));
    Ak(3,i)= abs(M2(r1(i))-M2(r1(i)-1));
    Ak(4,i)= abs(M3(r1(i))-M3(r1(i)-1));   
end

W1 = zeros(7,length(r1));
W1(1,:)= Dk(2,:)-Dk(1,:)*2;
W1(2,:)= Dk(3,:)-Dk(1,:)*3;
W1(3,:)= Dk(4,:)-Dk(1,:)*4;
W1(4:7,:)= Ak;

Y1 = t1(d1)-t1(r1);

%%
% Collect features in G2 for GLGP model as testing dataset(W2)

Dk2 = zeros(4, length(r2)); 

for i = 1:length(r2)
    Dk2(1,i) = R_phase(r2(i), r2(i)-1, F);
    Dk2(2,i) = R_phase(r2(i), r2(i)-1, M1);
    Dk2(3,i) = R_phase(r2(i), r2(i)-1, M2);
    Dk2(4,i) = R_phase(r2(i), r2(i)-1, M3);
end 

Ak2 = zeros(4, length(r2));
for i = 1:length (r2)
    Ak2(1,i)= abs(F(r2(i))-F(r2(i)-1));
    Ak2(2,i)= abs(M1(r2(i))-M1(r2(i)-1));
    Ak2(3,i)= abs(M2(r2(i))-M2(r2(i)-1));
    Ak2(4,i)= abs(M3(r2(i))-M3(r2(i)-1));   
end

W2 = zeros(7,length(r2));
W2(1,:)= Dk2(2,:)-Dk2(1,:)*2;
W2(2,:)= Dk2(3,:)-Dk2(1,:)*3;
W2(3,:)= Dk2(4,:)-Dk2(1,:)*4;
W2(4:7,:)= Ak2;

% Principal component analysis for all beats and G2, respectively
W = [W1 W2];
[coeffw,score,latent,~,explained] = pca(W');
[coeffw2,score2,latent2,~,explained2] = pca(W2');
%%
% Predict DNs for G2 with GLGP model
pair_w = W1';
dist_w = pdist(pair_w);
eu = prctile(dist_w, 25);
[Y2] = CovMatrix(W1',W2',Y1',100,eu,1,0.1);
t2 = t1(r2) + Y2';
n2 = interp1(t1, data_ppg, t2, "spline");

%%
% Figure2: PCA plot
%   subplot-1: PCA plot for labeled G1 and unlabled G2
%   subplot-2: PCA plot for G2 colored with prediction values
n1 = length(Y1);
n2 = length(Y2);
r = [-300 300];

figure(2);
subplot(1,2,1);
p1(1) = scatter3(score(1:n1,1),score(1:n1,2),score(1:n1,3),5,Y1',"filled"); colorbar;
hold on
p1(2)= scatter3(score(n1+1:end,1),score(n1+1:end,2),score(n1+1:end,3),5,"k", "filled");
axis equal
legend("Group1 (labeled)", "Group2 (unlabeled)",  'FontSize',18);
xlabel(sprintf('PC1 (%.1f%%)',explained(1)), 'FontSize',20)
ylabel(sprintf('PC2 (%.1f%%) ',explained(2)), 'FontSize',20)
zlabel(sprintf('PC3 (%.1f%%)',explained(3)), 'FontSize',20)
title(sprintf("PCA of all beats"), 'FontSize',20) ;
title(colorbar,"y(sec)",'FontSize',14);
caxis([0 0.5]);
xlim(r);
ylim(r);
zlim(r);

subplot(1,2,2);
scatter3(score2(:,1),score2(:,2),score2(:,3),5,Y2, "filled"); colorbar;
axis equal
legend("Group2 (predicted)", 'FontSize',18);
xlabel(sprintf('PC1 (%.1f%%)',explained2(1)), 'FontSize',20)
ylabel(sprintf('PC2 (%.1f%%) ',explained2(2)), 'FontSize',20)
zlabel(sprintf('PC3 (%.1f%%)',explained2(3)), 'FontSize',20)
title(sprintf("PCA of Group 2"), 'FontSize',20) ;
title(colorbar," Predicted y^*(sec)",'FontSize',14);
caxis([0 0.5]);
xlim(r);
ylim(r);
zlim(r);


%%
% Figure.3: Mark both GLGP predicted DNs and APPG predicted DNs for beats in G2 on PPG signal
figure(3);
plot(t1, data_ppg, "k"); hold on; grid on;
x(1)= scatter(t1(DNloc),data_ppg(DNloc),40,[0.1000 0.6000 0.3000], '^',"filled");
x(2)= scatter(t1(eloc(NonDN)),data_ppg(eloc(NonDN)),40, "filled", "r" );
x(3)= scatter(t2, n2, 40, "blue", "*");
xlim([6645 6660]);
xlabel("Time(seconds)", 'FontSize',16);
ylabel("PPG(au)", 'FontSize',16);
legend(x([1,2,3]),"G1: True DN","G2: APPG predicted DN","G2: GL-GP predicted DN", 'FontSize',14);

%% 
% Test GLGP model on G1
% Randomly divide G1 into two groups G1-1 and G1-2 as training and testing datasets
rd = randsample(n1,round(0.5*n1));
W12 = zeros(7, length(rd));
for i = 1:length(rd)
    W12(:,i) = W1(:,rd(i));
end
d12 = d1(rd);
W11 = W1;
Y11 = Y1;
d11 = d1;
W11(:,rd) =[];
Y11(:,rd) = [];
d11(:,rd) = [];

% Predict DNs for G1-2 by GLGP model with G1-1
pair1 = W11';
dist1 = pdist(pair1);
e1 = prctile(dist1, 25);
[Y12] = CovMatrix(W11',W12',Y11',40,e1,1,0.1);
t12 = t1(r1(rd))+ + Y12';
n12 = interp1(t1, data_ppg, t12, "spline");
%%
% Figure. 4 Mark predicted DNs on beats in G1-2 
figure(4);
plot(t1,data_ppg, "k"); hold on; grid on;
g(1)=scatter(t1(d11(1,:)), data_ppg(d11(1,:)),30, "r", "filled");
g(2)=scatter(t1(d12(1,:)), data_ppg(d12(1,:)), 30,[0.1000 0.6000 0.3000],"filled");
g(3)=scatter(t12, n12, 40,"b", "*");
xlabel("Time (seconds)","FontSize", 16);
ylabel("PPG (au)","FontSize", 16);
legend(g([1:3]),"G1-1: True DN", "G1-2: True DN", "G1-2: GLGP-predicted DN", "FontSize", 14)
xlim([6613 6628]);



function [theta]= R_phase(a,b, f)
v1 = f(b);
v2 = f(a);
zz = v1*conj(v2)/abs(v1)/abs(v2); 
zz = zz * exp(1i * pi);
theta = angle(zz);

end