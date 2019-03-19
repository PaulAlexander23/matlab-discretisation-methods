function DY = diff_fd2(x, Y, D)
    shape = size(Y);
    
    y = reshape(Y,numel(Y),1);
    
    dy = cell2mat(cellfun(@(D) D*y, D, 'UniformOutput', false));
    
    DY = reshape(dy,[shape,length(D)]);
end