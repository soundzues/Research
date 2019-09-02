function z=prbs(init,g)
% z=prbs(init,g)
% 2^n-1-bit PRBS based on initial string 'init'
% and polynomial represented by vector g (e.g., g=[7 1] => x^7+x+1).
z=init;
n=length(init);
for i=(n+1):(2^n-1)
    q=z(i-g(1));
    for j=2:length(g)
        q=xor(q,z(i-g(j)));
    end
z=[z q];
end