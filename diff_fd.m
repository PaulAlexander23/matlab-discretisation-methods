function dy = diff_fd(x, y, degree)
    global D
    
    %init_fd(x, 2);
    
    dy = cell2mat(arrayfun(@(d) D{d}*y,degree,'UniformOutput',false));
end
