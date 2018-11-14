
for ke1 = 1:nel
    for je1 = 1:12
        ke2 = ke1; %need to fix for cases with more than one element
        PlasList(je1,ke1,step+1) = ElemP(je1,ke1); %#ok<SAGROW>
    end
end