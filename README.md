# PhD research

## Simulate data

### Simulate identifiers (name, family name, country, gender, birth year)

- Get names, family names from 32 European countries using [these datasets](https://github.com/philipperemy/name-dataset) (Facebook Data Leak, 2019)
- Match countries 2 letters code based on [this dict](https://gist.github.com/mlisovyi/e8df5c907a8250e14cc1e5933ed53ffd)
- Generate age distributions based on: [the World Bank](https://data.worldbank.org/) (helped with Chatgpt), [the World Health Organization](https://www.who.int/countries/) (helped with Chatgpt), [the World Factbook](https://www.cia.gov/the-world-factbook/countries/)
- Use realistic population sizes from: [United Nations Department of Economic and Social Affairs](https://www.un.org/development/desa/pd/data-landing-page), [Eurostat](https://ec.europa.eu/eurostat/web/main/data/database) and [the World Bank](https://data.worldbank.org/) (helped with Chatgpt)

Only in python, see:
```
.
└── simulate data
    └── python
        └── generate identifiers.ipynb
```

### Generate covariates, treatment and outcome

- Generate covariates, (age + 4 gaussian variables with different parameters)
- Generate treatment $T \sim Bernoulli(p=0.4)$
- Generate outcome (either fixed or non-fixed treatment effect over individuals)
    - $Y = -10 + a T + b \exp(X_{4}) + c X_{3} X_{1} + d X_{5}$
    - $Y = -10 + a T X_{2} + b \exp(X_{4}) + c X_{3} X_{1} + d X_{5}$

In python and R, see:
```
.
└── simulate data
    ├── python
    |   └── generate association.ipynb
    └── R
        └── generate association.R
```

## Replicate **Estimate-Tethered Stopping Rule algorithm**

- A **Minimum Estimated Variance algorithm** developed in [Simultaneous record linkage and causal inference with propensity score subclassification](https://onlinelibrary.wiley.com/doi/10.1002/sim.7911)

In python and R, see:
```
.
└── replicate MEV
    ├── python
    |   └── launch MEV algo.ipynb
    └── R
        └── launch MEV algo.R
```

## Results 

The context studied in the [paper](https://onlinelibrary.wiley.com/doi/10.1002/sim.7911) is the one of an additive treatment effect (treatment effect is the same for all individuals).

On ```python``` we can see the evolution of the estimated variance of the treatment effect trough the linked records sets (in decreasing links confidence order) and the evolution of the estimated average treatment effect.

Solid blue lines represent the average evolution through the sets for 10 rounds (simulate data + apply MEV algorithm), surounding shaded areas represent the 95% confidence intervals (over the 10 rounds) through the sets of linked records. Solid orange line represents the tretment effect value we are trying to recover.

See images:
```
.
└── images
    ├── fixed treatment effect
    |   ├── ate
    |   └── variance
    └── non-fixed treatment effect
        ├── ate
        └── variance
```

We first observe that the algorithm does not work in a non-fixed treatment effect setting (when treatment effect is different for each individual).

We then notice an elbow on the variance plots that points out the best average treatment effect estimation, that we get for an early set (relying almost only on true links).