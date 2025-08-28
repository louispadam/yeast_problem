function return_data = influence(mass,alpha)
%INFLUENCE calculates the impact of cells in the signalling region on cells
%in the responsive region
arguments
    mass        % fraction of all cells in signalling region
    alpha       % strength of linear interaction
end

    return_data = -alpha*mass;

end