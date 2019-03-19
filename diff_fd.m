function dy = diff_fd(x, y, D)
    
    
    dy = cell2mat(arrayfun(@(D) D*y, D,'UniformOutput',false));
    
    
end
