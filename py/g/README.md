Using publically available data we can plot the monthly worldwide temperature anomaly from 1880-2011 from a variety of data sources. Three sources (GISS, HAD, and NOAA) provide data for the entire range of dates. Two other sources (RSS, UAH) provide data from 1977 onward. The temperature anomaly is the difference between the monthly average temperature and the average temperature for a reference period or baseline. Different datasets have different baselines. GISS: 1951-1980, HAD: 1961-1990, NOAA: 1961-1990, RSS: 1979-1998, UAH: 1979-1998. 

I have not yet converted these to a common baseline, but will do so in the future to improve comparisons.

This plot takes advantage of 
1. mouseover events in d3, 
2. connecting text mouseover to change path color, and 
3. reordering lines so the selected line is always on top.