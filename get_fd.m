function out = get_fd(deg, D, problemDeg)
    diffnum = size(deg,2);
    
    out = cell(diffnum,1);
    for j = 1:diffnum
        [~,index] = ismember(deg(:,j)',problemDeg','rows');
        out{j} = D{index};
    end
end