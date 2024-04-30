function [wc, Kp, taui, taud, ok] = findpid(G1, gm, Ni, al, w)

    phii = -90 + atan(Ni) * 180/pi;
    phiM = asin((1-al)/(1+al))*180/pi;
    PG = -180 + gm - phiM - phii;
    if nargin < 5
       [M,P0,W] = bode(G1);
    else
        [M,P0,W] = bode(G1,w);
    end

    N = 0;
    while P0(1) - N * 360 > 180
        N = N + 1;
    end

    P = P0 - N * 360;
    
    [m1,p1] = bode(G1,W(1));
    vs = p1 - PG - N * 360;
    n = 0;

    for i = 2:size(W,1)
        v = P(i) - PG;
        if v*vs < 0
            n = n+1;
            wcs(n) = W(i) - (v/(v-vs)*(W(i)-W(i-1)));
        end
        vs = v;
    end
   
    for j = 1:n
        td(j) = 1/(sqrt(al)*wcs(j));
        ti(j) = Ni/wcs(j);
        Cd = tf([td(j) 1], [td(j)*al 1]);
        Ci = tf([ti(j) 1], [ti(j) 0]);
        [mj,pj] = bode(G1*Cd*Ci,wcs(j));
        kp(j) = 1/mj;
        Gol(j) = kp(j)*Cd*G1*Ci;
    end

    nOK = 0;
    ok = 0;
    for k = 1:n
        Gcl(k) = feedback(Gol(k),1);

        try
            if (isstable(Gcl(k)))
                nOK = k;
                ok = ok + 1;
            end
        catch ME
            nOK = k;
            ok = ok + 1;
        end
    end
    if ok > 0
        Kp = kp(nOK);
        taud = td(nOK);
        taui = ti(nOK);
        wc = wcs(nOK);
    else
        if n > 0
            Kp = kp(n);
            taud = td(n);
            taui = ti(n);
            wc = wcs(n);
        else
            Kp = 1;
            taud = 0;
            taui = 1;
            wc = 1;
        end
    end
    sprintf('Found %d stable soultion(s) out of %d phase crossing(s)' ...
        , ok, n)
end

