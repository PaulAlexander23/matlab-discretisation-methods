function dy = diff_fd(x, Y, deg, D, problemDeg)
    shape = cellfun(@length,x)';
    diffnum = size(deg,2);
    
    y = reshape(Y,prod(shape),[]);
    
    dy = cell(diffnum,1);
    for j = 1:diffnum
        [~,index] = ismember(deg(:,j)',problemDeg','rows');
        
        dy{j} = reshape(D{index}*y,[shape,size(y,2)]);
    end
end