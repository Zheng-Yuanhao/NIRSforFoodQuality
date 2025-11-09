%%IntensityFigure
function IntensityFigure()

    xlabel('Wavelength (nm)', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial')
    ylabel('Intensity', 'FontSize', 18, 'FontWeight', 'bold', 'FontName', 'Arial')
    
    box on
    set(gca, 'LineWidth', 3); % 设置坐标轴边框宽度
    % 调整X轴标注的字体大小和加粗
    set(gca, 'FontSize', 14, 'FontWeight', 'bold');
    
    grid on; % 启用网格线
    set(gca, 'GridLineStyle', '--'); % 设置网格线为虚线

end