# Dirichlet Process Mixture

The Dirichlet Process is a (prior) distribution over distributions. You throw in a probability distribution (a so-called base distribution) and out you get a series of probability distributions. In the form of an equation: `G ~ DP(H, alpha)`. The scalar alpha defines how often exactly the same distribution is sampled from H.


```
                _______________
               |               |
               |               |
               |               |
       H  -->  |   Dirichlet   | --> G
    alpha -->  |    Process    |
               |               |
               |               |
               |_______________|

```

So, suppose we sample multiple G0, G1, G2, ..., then Gi and Gj can be exactly the same distribution even if their parameters come from a continuous space. For example Gi can be the normal distribution N(mu=1.120391948, sigma=0.23874921111), and Gj can be exactly the same distribution! 


# The cool stuff

Check [this presentation](http://perso.telecom-paristech.fr/~gfort/Slides/Barcelone14.pdf) by Fort for some ways to cope with very multimodal distributions. Their example uses 20 multinormal distributions.

```
pi = sum_{i=1}^20 N(mu_i, Sigma_i)
```

Our example is similar, but for 100 of these distributions, with Sigma homogenous, and most importantly `mu` for each distribution derived from a deterministic function depending on `i`:

```
pi = sum_{i=1}^100 N(f_i(mu), Sigma)
```

The approaches to cope with multimodality:

1. Wang-Laundau sampler. A few directions of meta-stability are defined beforehand. It is biasing the potential.

2. Equi-Energy sampler. These are tempering methods.


# Conventions

Due to idiosyncrasies of octave/matlab the code gravitated towards the following conventions:

* Data is stored column-wise: `model(1).mu = [0; 1]`. This allows extraction as one matrix by: `[model.mu]`.
* Data in nd-arrays is ordered by last index: `R(:,:,1) = rand(2,2)`, not as `S(1,:,:)` because that's harder to extract. 
* Data is displayed by explicitly using `disp`. 

## Dependencies

Install matlab or octave, in the latter case:

	sudo aptitude install octave octave-statistics

Make sure it is loaded as well, for example by adding to your `~/.octaverc` file the following command:

	pkg load statistics



