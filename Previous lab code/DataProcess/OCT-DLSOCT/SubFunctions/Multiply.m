function ret = Multiply(A, B)

	szA = size(A);
	szB = size(B);

    szret = [szA(1:end-1) szB(2:end)];
    ret = zeros([numel(A)/szA(end) numel(B)/szB(1)]);

    sqA = reshape(A, [numel(A)/szA(end) szA(end)]);
    sqB = reshape(B, [szB(1) numel(B)/szB(1)]);

    for (j=1:szA(end))
        ret = ret + sqA(:,j) * sqB(j,:);
    end

    ret = reshape(ret,szret);

