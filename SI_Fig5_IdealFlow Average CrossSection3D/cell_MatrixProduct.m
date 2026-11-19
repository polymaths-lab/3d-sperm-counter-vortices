% matrix product, but using  cell arrays

%the two matrix must be in cell format

function [P] = cell_MatrixProduct(Rot,V)

[m n]   = size(Rot); % arclenth index
[k l]   = size(V);
[np nf] = size(V{1});

if (n~=k)
    P= ('Error: cell dimension must agree')
    return
end

P = cell(m,l);
for i=1:m
    for j=1:l
        P{i,j}=zeros(np,nf);
    end
end

for i=1:m
    for j=1:l
        for p=1:n
        P{i,j} = P{i,j} + (Rot{i,p}.*V{p,j});
        end
    end
end

