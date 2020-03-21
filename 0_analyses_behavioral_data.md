## Figure 1 & Supplementary Figure 1 - ANALYSES  ##
Alexander E. Hausmann, March 2020

## Initial Setup

Set working directory (change to respective device)
```{r}
setwd("C:/Users/Hausmann/Desktop/Rossi_et_al_2020/")
```

Call `brms` library (for categorical/multinomial Bayesian models)
```{r}
suppressMessages(suppressWarnings(library(brms)))
```

Read in data
```{r}
tern_pref<-read.csv("ternary_final.csv",header=T,
                    stringsAsFactors = F)
```

## Data Handling

Calculate number of trials with response
```{r}
tern_pref$total_trials_with_response<-rowSums(tern_pref[,c("total_trials_melp_court_only",
                                                           "total_trials_cydno_court_only",
                                                           "total_trials_court_both")])
```

Discard males with 0 responses
```{r}
tern_pref<-tern_pref[tern_pref$total_trials_with_response>0,]
```

Calculate proportions for each male
```{r}
tern_pref$prop_trials_melp_court_only<-tern_pref$total_trials_melp_court_only/tern_pref$total_trials_with_response
tern_pref$prop_trials_cydno_court_only<-tern_pref$total_trials_cydno_court_only/tern_pref$total_trials_with_response
tern_pref$prop_trials_court_both<-tern_pref$total_trials_court_both/tern_pref$total_trials_with_response
```

Bloat up to observation level data
```{r}
# Create empty data frame
tern_pref_model<-data.frame(insectary_id=character(),Type=character(),
                            geno1=character(),geno17=character(),geno18=character(),
                            decision=character())
# Go row by row through table
for(expander in 1:nrow(tern_pref)){
  # Categorical decision
  decision<-c(rep("melp_only",tern_pref$total_trials_melp_court_only[expander]),
              rep("cyd_only",tern_pref$total_trials_cydno_court_only[expander]),
              rep("both",tern_pref$total_trials_court_both[expander]))
  # Bind together (repeat insectary ID, Type and genotypes accordingly)
  tern_pref_model<-rbind(tern_pref_model,
                         cbind(insectary_id=rep(tern_pref$insectary_id[expander],length(decision)),
                               Type=rep(tern_pref$Type[expander],length(decision)),
                               geno1=rep(tern_pref$geno1[expander],length(decision)),
                               geno17=rep(tern_pref$geno17[expander],length(decision)),
                               geno18=rep(tern_pref$geno18[expander],length(decision)),
                               decision=decision))
}
```

Create new dataset for only *cydno* backcross individuals. We have to drop unused levels. This is only for modelling purposes.
```{r}
tern_pref_model_cydBC<-tern_pref_model[tern_pref_model$Type=="CPx(CPxMP)",]
tern_pref_model_cydBC<-droplevels(tern_pref_model_cydBC)
```

Same subset for plotting (subset the initial non-bloated table)
```{r}
tern_pref_cydBC<-tern_pref[tern_pref$Type=="CPx(CPxMP)",]
tern_pref_cydBC<-droplevels(tern_pref_cydBC)
```


## Statistical modelling

We fit 2 models

* Model 1: *cydno*, *melpomene*, F1s, backcross to *melpomene* & backcross to *cydno* males (irrespective of their genotype in backcross to *cydno*) (dataset `tern_pref_model`)
* Model 2: backcross to *cydno* males with different genotypes on the chromosome 1 peak, different genotypes on the chromosome 17 peak and different genotypes on the chromosome 18 peak (dataset `tern_pref_model_cydBC`)

Recode contrasts (necessary to retrieve predictors from categorical brm model). First, save old contrasts
```{r}
old_contr<-options("contrasts")
options(contrasts = c("contr.sum", "contr.poly"))
```

Run the two models. Use `insectary_id` as random factor.

### Model 1

As for the later models as well, note that we always set a seed within the `brm` function.
```{r}
mod1 <- suppressMessages(brm(decision ~ Type + (1|insectary_id), data = tern_pref_model, 
                             family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=42))
```

### Model 2
```{r}
mod2 <- suppressMessages(brm(decision ~ geno1 + geno17 + geno18 + (1|insectary_id), data = tern_pref_model_cydBC, 
                             family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=43))
```


## Approximate leave-one-out cross-validation (LOO) for fixed effects of models

### Model 1 (LOO)

Check importance of `Type` as fixed effect in `mod1`. Define reduced model (drop `Type`).
```{r}
mod1_r <- suppressMessages(brm(decision ~ 1 + (1|insectary_id), data = tern_pref_model, 
                               family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=44))
```

Perform approximate leave-one-out cross-validation (LOO) with the full and the reduced variation of model 1.
```{r}
loo_1_full<-loo(mod1)
loo_1_r<-loo(mod1_r)
```

Compare loos between full model and reduced model
```{r}
loo_compare(list(loo_1_full,loo_1_r))
```

As a rule of thumb, dividing the difference in ELPD by its standard error (SE) can give an idea whether a fixed effect is important. A typical cutoff for calling a fixed effect relevant is if ELPD difference >= 1.96 units of its SE (95% confidence interval, as based on the normal distribution; with small sample sizes where normal approximation is inadequate, some people argue to go up to 4SE; since we have high sample sizes, this shouldn't be a problem).

What we can see from this output is:
Dropping `Type` leads to a big difference in units of standard error (ELPD difference = 7.72SE units). `Type` seems therefore a highly important factor!

### Model 2 (LOO)

Check importance of the different fixed effects in `mod2`. Define reduced models (drop always one of the fixed effects)
```{r}
mod2_r_geno1 <- suppressMessages(brm(decision ~ 1 + geno17 + geno18 + (1|insectary_id), data = tern_pref_model_cydBC, 
                                     family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=45))
mod2_r_geno17 <- suppressMessages(brm(decision ~ geno1 + 1 + geno18 + (1|insectary_id), data = tern_pref_model_cydBC, 
                                      family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=46))
mod2_r_geno18 <- suppressMessages(brm(decision ~ geno1 + geno17 + 1 + (1|insectary_id), data = tern_pref_model_cydBC, 
                                      family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=47))
```

Perform approximate leave-one-out cross-validation (LOO) with the full and the reduced variations of model 2.
```{r}
loo_2_full<-loo(mod2)
loo_2_r_geno1<-loo(mod2_r_geno1)
loo_2_r_geno17<-loo(mod2_r_geno17)
loo_2_r_geno18<-loo(mod2_r_geno18)
```

Compare loos between full model and reduced models
```{r}
loo_compare(list(loo_2_full,loo_2_r_geno1,loo_2_r_geno17,loo_2_r_geno18))
```

What we can see from this output is:

* A) Dropping `geno17` does not seem to affect predicitve capability of the model very much (ELPD difference = 0.70 SE units). This is clearly the least relevant factor. It will not be considered for producing the graphs in the graph markdown.
* B) Dropping `geno18` makes the model worse (ELPD difference = 2.14 SE units). `geno18` is an important factor!
* C) Dropping `geno1` makes the model worse (ELPD difference = 2.34 SE units). `geno1` is an important factor!

We drop `geno17` and perform again approximate leave-one-out cross-validation on this model, testing the effect `geno1` and `geno18` in the reduced model (`mod2_r_geno17`).

Define models reduced by either `geno1` or `geno18`.

```{r}
mod2_r_geno17_r_geno1 <- suppressMessages(brm(decision ~ 1 + 1 + geno18 + (1|insectary_id), data = tern_pref_model_cydBC, 
                                              family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=48))

mod2_r_geno17_r_geno18 <- suppressMessages(brm(decision ~ geno1 + 1 + 1 + (1|insectary_id), data = tern_pref_model_cydBC, 
                                               family="categorical",chains=5, iter=6000, warmup=3000, refresh=0,silent = TRUE,seed=49))
```

Perform approximate leave-one-out cross-validation (LOO) with the even more reduced variations of the model (where `geno17` and either `geno1` or `geno18` were dropped).
```{r}
loo_2_r_geno17_r_geno1<-loo(mod2_r_geno17_r_geno1)
loo_2_r_geno17_r_geno18<-loo(mod2_r_geno17_r_geno18)
```

Compare LOOs between model reduced by just `geno17` and models reduced by additionally `geno1` or `geno18`.
```{r}
loo_compare(list(loo_2_r_geno17,loo_2_r_geno17_r_geno1,loo_2_r_geno17_r_geno18))
```

What we can see from this output is:

* A) Dropping `geno18` leads to the same ELPD difference and SE as before (ELPD diff = 2.14 SE units). `geno18` is an important factor!
* B) Dropping `geno1` leads to a slightly bigger ELPD and SE as before (ELPD diff = 2.38 SE units). `geno1` is an important factor!


## Retrieve conditional effects

In all calls, we use `re_formula=NA` in order not to condition of the group-level effects.

### Model 1 (conditional effects)

Retrieve conditional effects from model 1.
```{r}
cond_eff_1<-conditional_effects(mod1,categorical=T)
```

### Model 2 (conditional effects)

Retrieve conditional effects from model 2. We used the model reduced by `geno17` (`mod2_r_geno17`), since `geno17` turned out to be of little importance in the model.

* We retrieve the conditional effects for genotype 18 (A vs. B) while keeping `geno1` at the grand mean, by setting it to `NA` (therefore only extracting the effect of genotype at the chromosome 18 peak alone). These predictors will be used for figure 1. 
* Second, we extract the predictors for genotype at QTL on chromosome 18 when genotype at QTL on chromosome 1 is homozygous. These predictors will be used for supplementary figure 1.
* Last, we extract the predictors for genotype at QTL on chromosome 18 when genotype at QTL on chromosome 1 is heterozygous. These predictors will be used for supplementary figure 1.
```{r}
cond_eff_2_geno18<-conditional_effects(mod2_r_geno17,  effects = "geno18", re_formula = NA,categorical=T,
                                       conditions = data.frame(geno1 = NA))

cond_eff_2_geno18_with_geno1_homo<-conditional_effects(mod2_r_geno17,  effects = "geno18", re_formula = NA,categorical=T,
                                                       conditions = data.frame(geno1 = "A"))

cond_eff_2_geno18_with_geno1_hetero<-conditional_effects(mod2_r_geno17,  effects = "geno18", re_formula = NA,categorical=T,
                                                         conditions = data.frame(geno1 = "B"))
```

Set contrasts back to old setting
```{r}
options(contrasts = old_contr$contrasts)
```


## Produce predictor tables for later plots

### Model 1 (predictor tables)

For model 1, we want a list with 5 table entries. Table 1 is for *cydno*, table 2 for *melpomene*, table 3 for F1s, table 4 for *cydno* backcross (actually not used for later graphs), table 5 for *melpomene* backcross.

Each table has three rows and three columns. The columns show predictor, lower 2.5% threshold of posterior, upper 97.5% threshold of posterior.
```{r}
estimator_table_1<-list()
temp_t<-cond_eff_1$`Type:cats__`[cond_eff_1$`Type:cats__`$effect1__=="CP",c("effect2__","estimate__","lower__","upper__")]
estimator_table_1[[1]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                 (1:3)[temp_t$effect2__=="both"]),2:4]
temp_t<-cond_eff_1$`Type:cats__`[cond_eff_1$`Type:cats__`$effect1__=="MP",c("effect2__","estimate__","lower__","upper__")]
estimator_table_1[[2]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                 (1:3)[temp_t$effect2__=="both"]),2:4]
temp_t<-cond_eff_1$`Type:cats__`[cond_eff_1$`Type:cats__`$effect1__=="CPxMP",c("effect2__","estimate__","lower__","upper__")]
estimator_table_1[[3]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                 (1:3)[temp_t$effect2__=="both"]),2:4]
temp_t<-cond_eff_1$`Type:cats__`[cond_eff_1$`Type:cats__`$effect1__=="CPx(CPxMP)",c("effect2__","estimate__","lower__","upper__")]
estimator_table_1[[4]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                 (1:3)[temp_t$effect2__=="both"]),2:4]
temp_t<-cond_eff_1$`Type:cats__`[cond_eff_1$`Type:cats__`$effect1__=="MPx(CPxMP)",c("effect2__","estimate__","lower__","upper__")]
estimator_table_1[[5]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                 (1:3)[temp_t$effect2__=="both"]),2:4]
```

### Model 2 (predictor tables)

For model 2, we want three lists, each with 2 table entries. The first list is for `geno18`, the second list is for `geno18` when `geno1` is homozygous, the third list is for `geno18` when `geno1` is heterozygous.

In each list, the first table is for backcross to *cydno* individuals homozygous at the chromosome 18 peak ("A"), the second table for individuals heterozygous at the chromosome 18 peak ("B"). Columns and rows are the same as before.
```{r}
estimator_table_2_geno18<-list()
temp_t<-cond_eff_2_geno18$`geno18:cats__`[cond_eff_2_geno18$`geno18:cats__`$effect1__=="A",c("effect2__","estimate__","lower__","upper__")]
estimator_table_2_geno18[[1]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                        (1:3)[temp_t$effect2__=="both"]),2:4]
temp_t<-cond_eff_2_geno18$`geno18:cats__`[cond_eff_2_geno18$`geno18:cats__`$effect1__=="B",c("effect2__","estimate__","lower__","upper__")]
estimator_table_2_geno18[[2]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                        (1:3)[temp_t$effect2__=="both"]),2:4]

estimator_table_2_geno18_with_geno1_homo<-list()
temp_t<-cond_eff_2_geno18_with_geno1_homo$`geno18:cats__`[cond_eff_2_geno18_with_geno1_homo$`geno18:cats__`$effect1__=="A",c("effect2__","estimate__","lower__","upper__")]
estimator_table_2_geno18_with_geno1_homo[[1]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                                        (1:3)[temp_t$effect2__=="both"]),2:4]
temp_t<-cond_eff_2_geno18_with_geno1_homo$`geno18:cats__`[cond_eff_2_geno18_with_geno1_homo$`geno18:cats__`$effect1__=="B",c("effect2__","estimate__","lower__","upper__")]
estimator_table_2_geno18_with_geno1_homo[[2]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                                        (1:3)[temp_t$effect2__=="both"]),2:4]

estimator_table_2_geno18_with_geno1_hetero<-list()
temp_t<-cond_eff_2_geno18_with_geno1_hetero$`geno18:cats__`[cond_eff_2_geno18_with_geno1_hetero$`geno18:cats__`$effect1__=="A",c("effect2__","estimate__","lower__","upper__")]
estimator_table_2_geno18_with_geno1_hetero[[1]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                                          (1:3)[temp_t$effect2__=="both"]),2:4]
temp_t<-cond_eff_2_geno18_with_geno1_hetero$`geno18:cats__`[cond_eff_2_geno18_with_geno1_hetero$`geno18:cats__`$effect1__=="B",c("effect2__","estimate__","lower__","upper__")]
estimator_table_2_geno18_with_geno1_hetero[[2]]<-temp_t[c((1:3)[temp_t$effect2__=="cyd_only"],(1:3)[temp_t$effect2__=="melp_only"],
                                                          (1:3)[temp_t$effect2__=="both"]),2:4]
```

## Save produced tables/lists

Save relevant tables for plotting Markdown
```{r}
save(list=c("tern_pref",
            "estimator_table_1",
            "tern_pref_cydBC",
            "estimator_table_2_geno18",
            "estimator_table_2_geno18_with_geno1_homo",
            "estimator_table_2_geno18_with_geno1_hetero"), file="analyses_fig1_suppl_fig_1.RData")
```
