
function sample=gigrndHandle(p, a, b_vec)
n=numel(b_vec);
sample = zeros(1,n);
for i=1:n
    sample(i) = gigrnd(p, a, b_vec(i), 1);
end
end


