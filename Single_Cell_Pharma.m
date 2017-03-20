%% Plot individual drug concentrations for each tracked region

% Log drug concentrations in a new variable
tot_tracks = length(tracksFinal);

for i = 1:tot_tracks
    for j = 1:length(tracksFinal(i,1).tracksCoordAmpCG)/8
        drug_props(i,j) = tracksFinal(i).tracksCoordAmpCG(1,(j-1)*8+4);

% Clear 0s from graphs
    drug_props(drug_props==0)=NaN; 
    end
end

% Plot all drug concentrations
figure;
hold all
for i=1:tot_tracks
  set(0,'DefaultAxesColorOrder',jet(5));
    plot(drug_props(i,1:nframes));
end
hold off


%% Create/ Plot all Average Values

for i=1:nframes
averages(i) = nanmean(drug_props(i,1:10));
standard_deviation(i) = mean(nanstd(drug_props(i,1:10)));
end
figure; 
plot(1:nframes, averages, 1:nframes, standard_deviation);

%% Plot min-max curves

for i = 1:nframes
    min_value(i) = min(drug_props(i,:));
    max_value(i) = max(drug_props(i,:));
end

figure;
plot(1:nframes, min_value, 1:nframes, max_value);
    
    
%% Save Output from object recognition phase in mfile
save('movie_data_zoomed_out.mat');