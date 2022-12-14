% código para obtener las efemérides de cada planeta

function [position, velocity] = Initialconditions_solarsystem(Save)
    if nargin < 1
        Save = "False";
    end
    Today = juliandate(2022,12,1);
    Planets = ["Mercury" "Venus" "Earth" "Moon" "Mars" "Jupiter" "Saturn" "Uranus" "Neptune" "Pluto"];
    
    position = zeros(length(Planets)+1,3);
    velocity = zeros(length(Planets)+1,3);
    
    % Position and velocity of the Sun
    position(1,:) = [0;0;0];
    velocity(1,:) = [0;0;0];
    
    for i = 1:length(Planets)
        
        [position(i+1,:),velocity(i+1,:)] = planetEphemeris(Today,"Sun",Planets(i));
    end
    
    if Save == "True"
        save("Initialposition_solarsystem.txt","position",'-ascii')
        save("Initialvelocity_solarsystem.txt","velocity",'-ascii')
    end
end