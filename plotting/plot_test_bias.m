function [] = plot_test_bias(match, center, W)
    
    title = "Testing Bias";
    
    figure;
    hold on
    bar(1:size(W,3), match);
    bar(1:size(W,3), center);
    
end