function eta = learnrate(iter, ns, eta_min, eta_max)
j = fix(iter/(2*ns));
if  (2*j*ns <= iter) && (iter <= (2*j+1)*ns)
        eta = eta_max-(iter-(2*j)*ns)*((eta_max-eta_min)/ns);
        %eta = eta_min+(iter-2*j*ns)*((eta_max-eta_min)/ns);
else
        eta = eta_max-(iter-(2*j+1)*ns)*((eta_max-eta_min)/ns);
end
end