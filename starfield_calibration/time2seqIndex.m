function seqIndex = time2seqIndex(time)

    time_th = 30e-4; % 30 sec 
    A= pdist2(time,time) < time_th;
    G = graph(A, 'OmitSelfLoops'); 
    seqIndex  = conncomp(G);

end