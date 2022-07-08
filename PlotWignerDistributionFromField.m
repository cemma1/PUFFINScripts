%% Plot Wigner Distribution for Harmonic Field
load xf
load taxis_field
 
dt = abs(taxis_power(2)-taxis_power(1));
xfcut = xf(1000:1350);% Cut edges of E field signal
t = [0:1:length(xfcut)-1]*dt;
 
plotWignerDist(xfcut,t,1,0);
 
%% Load the field at different locations and plot the wigner Dist
%load xfield
%load yfield
load xfieldNw10to40
load yfieldNw10to40
% Set the slippage manually 
%slippage(9:19) = -20*[1:11];% one wavelength per period for the first 19 periods
%slippage(20:39) = -20*11-round(20/5)*[1:20];% one fifth harmonic wavlength per period for the next 20 periods
slippage(10:20) = -20*[1:11];% one wavelength per period for the first 19 periods
slippage(21:40) = -20*11-round(20/5)*[1:20];% one fifth harmonic wavlength per period for the next 20 periods
close all
 
for ij =10:1:size(xfield,1)
    idxcut = 450:1350;
    xfcut = xfield(ij,idxcut);
    yfcut = yfield(ij,idxcut);
    totalField = xfcut + yfcut;
    totalFieldShifted = fliplr(circshift(totalField,slippage(ij)));
    t = [0:1:length(xfcut)-1]*dt;
    plotWignerDist(totalFieldShifted,t,ij-9,1);
end
%% Plot the field at different point along the undulator to figure out the slippage
% for ij = 9:1:10
%      idxcut = 750:1500;
%    figure(1)
%    plot(circshift(xfield(ij,idxcut),(ij-8)*-1.0*20));hold on
% end
%%
function plotWignerDist(field,t,ij,makeGif)
[wignerDist,fw,tw] = wvd(field,'smoothedPseudo');
% Normalize the wigner Dist to its peak value
wignerDist = wignerDist/max(max(wignerDist));
% Re-scale the frequency axis so the fundamental is at fw = 1
fw = fw*11;
% Re-scale the time axis so it goes from 0 to max(t)
tw = tw * max(t)/max(tw);
% Caculate the field in the fundamental and fifth harmonics
idx = fw<2;
ExFundamental = sum(wignerDist(idx,:),1);
idx = fw>4;
ExFifthHarm= sum(wignerDist(idx,:),1);
 
figure(1);
set(gcf,'Position',[100 100 700 900])
subplot(3,1,1);plot(t,field);xlim([0,t(end)]);xlabel('t [fs]');ylabel('|E (t)| [a.u.]');ylim([-1.1,1.1]*max(field));title(sprintf(['Nw = ',num2str(ij+8)]))
subplot(3,1,2);imagesc(tw,fw,wignerDist);xlabel('t [fs]');ylabel('\omega/\omega_0');ylim([0,6]);colormap(parula(7));
title('Wigner distribution')%colorbar
 
subplot(3,1,3);
yyaxis left
    plot(tw,ExFifthHarm);ylim([0,1.1*max(max(ExFifthHarm),max(ExFundamental))])
    ylabel('|E (t)|^2 (5\omega_0) [a.u.]');hold off;
    yyaxis right
    plot(tw,ExFundamental);ylim([0,1.1*max(ExFundamental)])
    xlim([0,tw(end)]);xlabel('t [fs]')
    ylabel('|E (t)|^2 (\omega_0) [a.u.]');hold off;
 
    drawnow   
    if makeGif
% make a Gif
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ij ==1
        imwrite(imind,cm,'WignerDistGif','gif','Loopcount',inf);
    else
        imwrite(imind,cm,'WignerDistGif','gif','WriteMode','append');
    end
    end
    
end
%% Plot the spectrogram for the field
% [spectroGram,fSpec,tSpec] = spectrogram(xfcut,'power');
% 
% figure
% imagesc(fSpec,tSpec,abs(spectroGram));
% colorbar
