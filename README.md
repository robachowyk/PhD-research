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
├── simulate data
    ├── python
        └── generate identifiers.ipynb
```

### Generate covariates, treatment and outcome

- Generate covariates, (age + 4 gaussian variables with different parameters)
- Generate treatment $T \sim Bernoulli(p=0.4)$
- Generate outcome $Y = -10 + a T X_{1} + b \log(X_{4}) + c X_{2} X_{3} + d X_{5}$

In python and R, see:
```
.
├── simulate data
    ├── python
    |   └── generate association.ipynb
    └── R
        └── generate association.R
```

## Replicate **Estimate-Tethered Stopping Rule algorithm**

- A modified **Minimum Estimated Variance algorithm** developed in [Simultaneous record linkage and causal inference with propensity score subclassification](https://onlinelibrary.wiley.com/doi/10.1002/sim.7911)
- 
```
In python and R, see:
.
├── replicate ETSR
    ├── python
    |   └── ETSR.ipynb
    └── R
        └── ETSR.R
```

blablabla
