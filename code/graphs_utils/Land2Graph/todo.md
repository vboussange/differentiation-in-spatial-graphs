# TODO

## Script
- Read the 1km dataset
- slice it in windows of 100x100km
- perform a loop over those windows, from which graphs are extracted
    - let's first consider only one type of habitats (say `lvl2_frac_1km_ver004/iucn_habitatclassification_fraction_lvl2__104_Forest – Temperate__ver004.tif`)
    - [x] extract the graphs in an array, by considering that nodes are connected if neighbours (consider 8 neihbours, i.e. direct neighbours and diagonal)
    > for this, use LightGraphs.jl. See example
```julia
using LightGraphs
A = [0 1; 1 0] # adjacency matrix
g = SimpleGraph(A)
```
    - compute the mean betweenness_centrality for each of the graph, and then average over all graphs for the window
```julia
using Statistics
mean(betweenness_centrality(g))
```
- recreate a raster with average of betweenness_centrality for each 100x100km window



## questions
- What is the difference between `iucn_habitatclassification_composite_lvl1_ver004` and `lvl1_frac_1km_ver004`?
> For me, composite shows all data and frac_1km shows data for each habitat aggregated at 1 km

- Whate is `lvl2_changemasks_ver004.zip` and `lvl2_frac_change_1km_ver004.zip`?


- what are the values inside the rasters? I do not understand

# Benji meeting
> level2_frac1km

## Hengudan mountains

coordinates : 
Top left :  latitude 31°50'4.05"N; longitude :  94° 7'55.83"E (decimal : 31.834458, 94.132175)
bottom left :  22° 8'3.44"N;  94° 9'52.28"E
bottom right:  22° 5'52.20"N; 101°52'25.13"E
top right:  31°49'48.83"N; 101°52'4.15"E
Description ([NASA](https://earthobservatory.nasa.gov/images/8098/hengduan-mountains-china))
In southwestern China, just northeast of the borders of India and Myanmar, lies a mountain range known as the Hengduan Shan (sometimes Heng-tuan Shan). The mountains run generally north-south with a slight tilt toward the west at their northern end. Like all mountain ranges, the Hengduan consist of peaks and valleys, but these mountains provide an exceptional range of climates, due to changes of both elevation and latitude. While some valleys are frost-free throughout the year, some mountain peaks sport glaciers.

The Moderate Resolution Imaging Spectroradiometer (MODIS) on NASA’s Terra satellite captured this image of the mountain range and the surrounding area on December 25, 2000. Taken near the start of the northern hemisphere winter, the picture shows plenty of snow-capped mountains, the icy whiteness forming dendritic patterns atop ridges. Just below the snowcaps are forested slopes. The valley floors appear primarily brown, contrasting with the lush greenness of neighboring lowlands in Myanmar. Yet the seeming brownness in the Hengduan Mountains can be misleading. These mountains actually support some of the world’s greatest diversity of flowering plants.

The vegetation of the Hengduan Mountains includes many species that would be surprisingly familiar to gardeners from across the world. These mountains are where many of the flowers so popular in gardens first evolved—geraniums, lilies, and lady slipper orchids, to name just a few. In fact, the mountains of southwestern China are home to an estimated 12,000 species, and roughly 3,500 of them are endemic, found naturally nowhere else. Although prized by botanists, the Hengduan Mountains’s biodiversity is under threat from overgrazing, firewood collection, and development.