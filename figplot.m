% clear all;
hf = figure;
set( hf, 'color', 'white');
SNR=0:2.5:15;
% load BER;
%BER = zeros(4,7);
% 1-bit
% BER(1,:) = [0.181	0.1958	0.2088	0.2145	0.2162	0.2165	0.2165];
% BER(2,:) = [0.2122	0.2032	0.1984	0.1932	0.1885	0.1855	0.1842];
% BER(3,:) = [0.2092	0.1902	0.1752	0.1609	0.1482	0.1394	0.1374];
% BER(4,:) = [0.1845	0.1519	0.1199	0.0832	0.0325	0.0033	0.0005];
% 2-bit
% BER(1,:) = [0.0836	0.043	0.0168	0.0056	0.0019	0.0008	0.0004];
% BER(2,:) = [0.0813	0.0395	0.014	0.004	0.0011	0.0004	0.0002];
% BER(3,:) = [0.0787	0.0361	1.17E-02	0.0029	0.0006	0.0002	0.0001];
% BER(4,:) = [0.0727	0.0299	8.00E-03	0.0015	0.0002	0	0];
%save BER BER;
semilogy( SNR, BER(1,:), '-r','LineWidth',1.5);
hold on;
semilogy( SNR, BER(2,:), '--r','LineWidth',1.5);
hold on;
semilogy( SNR, BER(3,:), '-b','LineWidth',1.5);
hold on;
semilogy( SNR, BER(4,:), '--b','LineWidth',1.5);
% hold on;
% semilogy( SNR, BER(3,:), '-.k','LineWidth',1.5);
% hold on;
% semilogy( SNR, BER(4,:), ':k','LineWidth',1.5);

set(hf,'position',[0,0,760,468])
%grid on;
% legend('100% Mixed Adc','95% Mixed Adc','90% Mixed Adc','80% Mixed Adc')
legend('80% Mixed 1-bit','80% Mixed 1-bit Opt','70% Mixed 1-bit','70% Mixed 1-bit Opt')
xlabel('\fontsize{14}SNR(dB)');
ylabel('\fontsize{14}BER');