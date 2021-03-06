function [Ch_diff,gradient_Ch_diff] = Lsquare_Ch_difference (design_variable,enforceChdiff,epsilon,initial_value)

global post_info
ncomponents = 6;
Ch_diff = 0;
gradient_Ch_diff = zeros(length(design_variable),1);
for i = 1:ncomponents
    [costval,gradientval] = enforceChdiff(design_variable,0,i);
    Ch_diff = Ch_diff + costval;
    gradient_Ch_diff = gradient_Ch_diff + gradientval;
end
ct = 1;
Ch_diff = ct*Ch_diff/initial_value - epsilon;
gradient_Ch_diff = ct*gradient_Ch_diff/initial_value;

end