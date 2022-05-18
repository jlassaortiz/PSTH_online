function dy = takens(t, y , env_song)
% myOde() 1D linear dynamical system for computing song envelope

dy = (-1/0.001)*y + abs(env_song(uint32(t),1));

end

