


legend(...
    '$eig(\hat{\mathcal{L}})$',...
    '$eig(\overline{\mathbf{F}}_{M,0})$',...
    '$eig(\overline{\mathbf{F}}_{M,2})$',...
    '$eig(\overline{\mathbf{F}}_{M,4})$',...
    '$eig(\overline{\mathbf{F}}_{M,6})$',...
    '$eig(\overline{\mathbf{F}}_{M,8})$',...
    'Interpreter', 'latex', 'FontWeight', 'bold','FontSize',20);




legend(...
    '$eig(\hat{\mathcal{L}})$',...
    '$eig(\overline{\mathbf{F}}_{M,1})$',...
    '$eig(\overline{\mathbf{F}}_{M,3})$',...
    '$eig(\overline{\mathbf{F}}_{M,5})$',...
    '$eig(\overline{\mathbf{F}}_{M,7})$',...
    'Interpreter', 'latex', 'FontWeight', 'bold','FontSize',20);


set(gca, 'FontName', 'Arial');
xlabel('Real','Interpreter', 'latex','FontWeight', 'bold');
ylabel('Imag',...
    'Interpreter', 'latex','FontWeight', 'bold');

% ylabel('$ln{\langle}e^{\frac{1}{L}\sum_{j=1}^{L}a_{j}^{\dagger}a_{j}}{\rangle}_{vac}$',...
%     'Interpreter', 'latex','FontWeight', 'bold');


box on;
title('L=4, J=1, \gamma_l=\gamma_g=0.1, \gamma_t=1');
set(gca,'FontSize',35,'LineWidth',3); 
xlim([-3,0.2])

set(gcf, 'PaperPositionMode', 'auto'); % 自动设置
set(gcf, 'PaperSize', [8 6]); % 设置纸张大小
set(gcf, 'PaperPosition', [0 0 8 6]); % 设置位置为零，去除边距
exportgraphics(gcf,'Spectrum_odd_L_4_gam1_0_5_gamL_0_5.eps','ContentType','vector','BackgroundColor','none');