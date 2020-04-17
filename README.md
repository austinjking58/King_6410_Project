# King_6410_Project

In this project, I apply an angular two point correlation function to stellar data
with the goal of finding which color of stars (red or blue) are more strongly correlated.
This was motivated by a similar process applies to galaxies. It was found that red 
galaxies are more strongly clustered than blue. I believe that we will see this same
result for stars, although blue stars will likely have a fairly strong correlation as 
well. 

I used data from both SDSS and Gaia to carry out this project. I took red stars as being 
those with an effective temperature below 3500K and blue stars with effective temperature
above 7500K. I then calcualted the angular correlation using astroML's two_point_angular
function. I repeated this process for specific regions as well. In galactic coordinates, I 
took stars that have galactic latitude between -40 degrees and 40 degrees and called these
my Disk Stars. I then took the remaining stars, those below -40 degrees and above 40 degrees
and called these my Halo Stars. 

Overall, my results were as expected. All plots showed a stronger correlation of the postitions
of red stars than blue stars except for the halo stars from Gaia data. I expect this is simply a 
sample error and is due to the fact that we did not have enough stars in this region from the 
Gaia data to make a more accurate estimate. I applied error bars to this plot to show that the 
very high coorelation for these blue stars has a significant error associated with it. I also 
found that stars in the disk have stronger correlation than stars in the halo, which is expected
due to the incresed stellar density in the disk.
