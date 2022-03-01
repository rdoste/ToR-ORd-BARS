% Emulate the ternary operator x = (cond) ? iftrue : iffalse
function y = ifthenelse(cond, iftrue, iffalse)
    if (cond)
        y = iftrue;
    else
        y = iffalse;
    end
end
