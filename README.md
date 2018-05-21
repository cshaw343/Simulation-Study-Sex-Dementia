# Simulation-Study-Sex-Dementia

**Variable Definitions**

Cij:   Cognitive function for person i at time j<br />
U:     Unmeasured/underlying variable<br />
Sij:   Survival for person i at time j<br />

**Model for Cognitive Function**

![Cij Equation](https://latex.codecogs.com/gif.latex?Cij%20%3D%20%5Cbegin%7Bcases%7D%20b00%20&plus;%20z0i%20&plus;%20b01sex%20&plus;%20b02age0c50%20&plus;%20b03U%20&plus;%20%5C%5C%28b10a%20&plus;%20z1i%20&plus;%20b11sex%20&plus;%20b12age0c50%20&plus;%20b13U%29j%20&plus;%20%5Cvarepsilon_%7Bij%7D%20%26%200%20%5Cleq%20j%20%3C%2020%5C%5C%20%26%5C%5C%20b00%20&plus;%20z0i%20&plus;%20b01sex%20&plus;%20b02age0c50%20&plus;%20b03U%20&plus;%2020%28b10a%20-%20b10b%29%20%5C%5C%28b10b%20&plus;%20z1i%20&plus;%20b11sex%20&plus;%20b12age0c50%20&plus;%20b13U%29j%20&plus;%20%5Cvarepsilon_%7Bij%7D%20%26%2020%20%5Cleq%20j%20%3C%2035%5C%5C%20%26%5C%5C%20b00%20&plus;%20z0i%20&plus;%20b01sex%20&plus;%20b02age0c50%20&plus;%20b03U%20&plus;%2020%28b10a%20-%20b10b%29%20&plus;%2035%28b10b%20-%20b10c%29%20%5C%5C%28b10c%20&plus;%20z1i%20&plus;%20b11sex%20&plus;%20b12age0c50%20&plus;%20b13U%29j%20&plus;%20%5Cvarepsilon_%7Bij%7D%26%2035%20%5Cleq%20j%20%5Cend%7Bcases%7D)

**Definitions of regression parameters used to generate models for Cij**

b00:   Group mean cognitive intercept (baseline level of cognitive function) for females<br />
b01:   Effect of sex on level of cognitive function at baseline<br />
b02:   Effect of a 1-year change in baseline age (in years, centered at 50) on level of cognitive function at baseline<br />
b03:   Effect of a 1-unit change in U on level of cognitive function at baseline<br />
b10a:  Group mean cognitive slope (annual rate of cognitive change) for females aged 50-70<br />
b10b:  Group mean cognitive slope (annual rate of cognitive change) for females aged 70-85<br />
b10c:  Group mean cognitive slope (annual rate of cognitive change) for females aged 85+<br />
b11:   Effect of sex on annual rate of cognitive change<br />
b12:   Effect of a 1-year change in baseline age (in years, centered at 50) on annual rate of cognitive change<br />
b13:   Effect of a 1-unit change in U on annual rate of cognitive change<br />

**Generation of survival times for individuals**  

Each person's underlying time to death is generated for each time interval between cognitive assessments under an exponential survival distribution conditional on the past, provided the person has not died in a previous interval.  If the person's generated survival time exceeds the length of the interval between cognitive assessments j and j + 1, they are considered alive at cognitive assessment j + 1. A new survival time is generated for the next interval conditional on history up to the start of the interval. This process is repeated until the person's survival time falls within a given interval or the end of the study, whichever comes first. <br />

A person's survival time for a given time interval at risk is generated using the inverse cumulative hazard function transformation formula described by Bender et al. (Stat Med 2011)


