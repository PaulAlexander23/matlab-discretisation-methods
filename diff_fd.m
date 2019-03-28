function dy = diff_fd(x, Y, deg, D, problemDeg)
    shape = size(Y);
    diffnum = size(deg,2);
    
    y = reshape(Y,numel(Y),1);
    
    dy = cell(diffnum,1);
    for j = 1:diffnum
        [~,index] = ismember(deg(:,j)',problemDeg','rows');
        
        dy{j} = reshape(D{index}*y,shape);
    end
end