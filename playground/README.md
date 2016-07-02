# Playground

This is an experiment around lines with angular representations, see `angular_line`. 

## Sampling

It uses `angularrnd` to generate points on a line segment. It uses the so-called normal form of a line (angle, distance to origin) to do so. 

    % LP, the closest point on the line w.r.t. the origin 
    LP = [T.p*cos(T.theta) T.p*sin(T.theta)];

It will transform a line segment from [d1-d2] on the x-axis according to these line coordinates. And it will add noise perpendicular to the line.

## Probability density function

It uses `logangularpdf` to calculate the log-likelihood for an individual combination of angle and distance for N observations.

## Visualization

To visualize the pdf, every combination of angle and distance is calculated with a picked resolution around these known values. This is only to show that the probability density function is correct:

![Visualization line segment](figures/line_from_5_to_8_with_theta_p_random.png)

## Notes

One thing that has be kept in mind is that the current visualization visualizes the likelihood, not the log-likelihood. Due to the extra exponential it can happen that the calculations return `Inf` for the likelihood calculations. It is recommended to use only log-likelihood in MCMC sampling.
