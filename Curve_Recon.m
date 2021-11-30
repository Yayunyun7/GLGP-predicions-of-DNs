%% 10-15 mins
ts = 10;
[tfr, tfrtic, tfrsq, ConceFT, tfrsqtic] = ConceFT_sqSTFT_C(data_ppg(1,1+ts*60*200:(ts+10)*60*200)',0, 0.1, 0.0002, 1, 2001, 1, 6, 1, 1, 0);

%%
idx = find(tfrsqtic*fs2 < 3);    

%% hop =1
[h,Dh,tt] = hermf(2001,1,6);
[recon] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 200, c1, 0.06, h(200*5+1)) ;
PPGfund = recon;
[recon2] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 200, c1*2, 0.06, h(200*5+1)) ;
PPGmulti1 = recon2;
[recon3] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 200, c1*3, 0.06, h(200*5+1)) ;
PPGmulti2 = recon3;

%%
figure(2)
plot(timesq, PPGfund,"r", "LineWidth", 2.0);
hold on;
plot(timesq, PPGmulti1,"b", "LineWidth", 2.0);
hold on;
plot(timesq, PPGmulti2,"k", "LineWidth", 2.0);

yu = -3*10^4;
yl = 3*10^4;



xlim([3600 3610]);
ylim([yu yl]);
xlabel("Time(sec)");
ylabel("PPG (au)");
xticks([3600:0.5:3610]);


%%
[recon3] = Recon_sqSTFT_v2(tfrsq, tfrsqtic, 10, c1*3, 0.03, h(200*6+1)) ;
PPGmulti3 = recon3;

hold on;
plot(timesq, PPGmulti2, "LineWidth", 2.0);
hold on;
plot(timesq, PPGmulti3, "LineWidth", 2.0);

%%%% 
figure(1);
imageSQ(timesq,tfrsqtic*200,abs(tfrsq),.995),colormap(1-gray); hold on;
plot(timesq,tfrsqtic(c1')*200,"LineWidth",2.0);
xlim([251 255]);
ylim([0 5]);
yticks([0:0.1:5]);
