function [] = plot_test_bias(stim,resp,varargin)
    
    title = "Testing Bias";
    
    figure;
    hold on
    bar([resp.match,resp.center]);
    
    legend('match','center');
    
end