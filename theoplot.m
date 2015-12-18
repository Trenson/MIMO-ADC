%plot for theory and simulation
hf = figure;
set( hf, 'color', 'white');
SNR=0:2.5:15;
semilogy( SNR, BER(1,:), '^r','LineWidth',1.5);
hold on;
semilogy( SNR, BER(2,:), 'ob','LineWidth',1.5);
hold on;
semilogy( SNR, BER(3,:), '^m','LineWidth',1.5);
hold on;
semilogy( SNR, BER(4,:), 'og','LineWidth',1.5);
% hold on;
% semilogy( SNR, BER(5,:), '--r','LineWidth',1.5);
% hold on;
% semilogy( SNR, BER(6,:), '--b','LineWidth',1.5);
% hold on;
% semilogy( SNR, BER(7,:), '--m','LineWidth',1.5);
% hold on;
% semilogy( SNR, BER(8,:), '--g','LineWidth',1.5);
set(gca,'XTick',[0 2.5 5 7.5 10 12.5 15]);
% set(gca,'XTick',0:2.5:15);  
% set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})  
set(hf,'position',[0,0,760,468])
legend('1-bit LMMSE simulation','1-bit AMP simulation','3-bit LMMSE simulation','3-bit AMP simulation',...
'1-bit LMMSE theory','1-bit AMP theory','3-bit LMMSE theory','3-bit AMP theory');
xlabel('\fontsize{14}SNR(dB)');
ylabel('\fontsize{14}BER');