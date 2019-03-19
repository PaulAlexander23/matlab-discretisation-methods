function dy = diff_fd_2d(x, y, degree)
    global D
    
    %init_fd(x, 2);
    dy = zeros([size(y),length(degree)]);
    shape = size(y);
    y = reshape(y,[numel(y),1]);
    for di = 1:length(degree)
        dy(:,:,di) = reshape(...
            D{degree(1,di)+1,degree(2,di)+1}*y,...
            shape);
    end
end
