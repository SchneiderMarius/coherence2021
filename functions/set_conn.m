function [W] = set_conn_log(conpar)
    
    d = ones(length(conpar.cell),1); % The diagonal values
    t = triu(bsxfun(@min,d,d.').*rand(sum(conpar.cellnum)),1); 


    % The upper trianglar random values
    M = diag(d)+t+t.'; % Put them together in a symmetric matrix
    W = zeros(length(conpar.cell),length(conpar.cell));
    for ct1 = 1 : length(conpar.celltype)
        for ct2 = 1 : length(conpar.celltype)
            X = M(find(conpar.cell==ct1),find(conpar.cell==ct2))<conpar.prob(ct1,ct2);
            W(find(conpar.cell==ct1),find(conpar.cell==ct2)) = X;  
            W(find(conpar.cell==ct1),find(conpar.cell==ct2)) = W(find(conpar.cell==ct1),find(conpar.cell==ct2))*conpar.syn(ct1,ct2);  
        end
    end

    W = W - diag(diag(W));

end