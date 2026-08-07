#'Suggested analysis pipeline for your paper

#'Given the mechanistic nature of your simulator, I would present the results in the following order:

#' Linear model (ANOVA + effect sizes) to provide interpretable baseline effects.
#' Generalized Additive Model (GAM) to identify nonlinear responses and thresholds.
#' Gaussian Process Regression as a high-fidelity surrogate model for the simulation.
#' Random Forest to assess nonlinear variable importance.
#' SHAP values to explain how each predictor influences each response.
#' Partial dependence (and optionally ICE) plots for the most influential variables.
#' Sobol global sensitivity indices, preferably computed from the Gaussian Process surrogate, to quantify both first-order and total effects.
#' Cross-validation metrics (RMSE, MAE, R2) to compare predictive performance across all models.
#' This workflow combines interpretability, predictive accuracy, and rigorous global sensitivity analysis,
#' making it well suited for publication in computational biology, ecological modeling, or evolutionary rescue studies.

# resposes to loop over
responses <- c(
  "final_g_bar",
  "final_z_bar",
  "final_var_g",
  "final_N",
  "persistence_time"
)
# 1. Packages
packages <- c(
  "tidyverse",
  "caret",
  "kernlab",
  "randomForest",
  "mgcv",
  "pdp",
  "iml",
  "fastshap",
  "DALEX",
  "sensitivity",
  "lhs",
  "DiceKriging"
)

invisible(lapply(packages, function(p){
    if(!require(p, character.only=TRUE))
        install.packages(p)
    library(p, character.only=TRUE)
}))
2. Prepare data
dsd$tipo_plantilla <- factor(dsd$tipo_plantilla)
dsd$modo_temporal  <- factor(dsd$modo_temporal)

responses <- c(
    "final_g_bar",
    "final_z_bar",
    "final_var_g",
    "final_N",
    "persistence_time"
)

predictors <- c(

"A_max",
"omega2",
"theta1",
"K_theta",

"g_inicial",
"N0",
"p_max",

"paso_introduccion",
"mutacion_sd",
"sigma_e",
"D",
"k",
"periodo_temporal",

"tipo_plantilla",
"modo_temporal"
)
3. Linear model
response <- "final_g_bar"

formula <- as.formula(

paste(
response,
"~",
paste(predictors,
collapse=" + ")
)
)

lm_fit <- lm(formula,data=dsd)

summary(lm_fit)

anova(lm_fit)

For publication I recommend

car::Anova(lm_fit,type=2)

instead of the default ANOVA.

4. Gaussian Process Regression

I recommend DiceKriging, which is much better than kernlab::gausspr for response surfaces.

gp_fit <- km(

formula = formula,

design = dsd[,predictors],

response = dsd[[response]],

covtype="gauss"

)

summary(gp_fit)

Prediction

pred <- predict(gp_fit,newdata=dsd)

head(pred$mean)
5. Random Forest
rf_fit <- randomForest(

formula,

data=dsd,

importance=TRUE,

ntree=1000

)

importance(rf_fit)

varImpPlot(rf_fit)

For publication

use

caret::varImp(rf_fit)
6. Partial dependence plots
partial(

rf_fit,

pred.var="A_max"

) %>%

autoplot()

For two variables

partial(

rf_fit,

pred.var=c("A_max","omega2")

)
7. SHAP values

I recommend fastshap

X <- dsd[,predictors]

shap <- fastshap::explain(

rf_fit,

X=X,

pred_wrapper=function(object,newdata){

predict(object,newdata)

},

nsim=100

)

Mean SHAP importance

colMeans(abs(shap))

Visualization

library(shapviz)

sv <- shapviz(shap,X)

plot(sv)
8. GAM

This is one of the strongest models for publication.

gam_formula <-

as.formula(

paste(

response,

"~",

paste(

paste0("s(",predictors[predictors!="tipo_plantilla" &
predictors!="modo_temporal"],")"),

collapse=" + "

),

"+",

"tipo_plantilla",

"+",

"modo_temporal"

)

)

gam_fit <-

gam(

gam_formula,

data=dsd,

method="REML"

)

summary(gam_fit)

plot(gam_fit)
9. Sobol sensitivity indices

Since your simulator already exists, define

model <- function(X){

predict(

rf_fit,

newdata=as.data.frame(X)

)

}

Generate samples

X1 <- randomLHS(1000,length(predictors)-2)

X2 <- randomLHS(1000,length(predictors)-2)

Compute

sobol <- sobolSalt(

model=NULL,

X1=X1,

X2=X2,

nboot=100

)

tell(

sobol,

model(rbind(X1,X2))

)

print(sobol)

For a Gaussian Process emulator, replace rf_fit with gp_fit in the prediction function to estimate Sobol indices from the surrogate model, which is a common approach in uncertainty quantification.

10. Model comparison
pred_lm <- predict(lm_fit)

pred_rf <- predict(rf_fit)

pred_gp <- predict(gp_fit)$mean

pred_gam <- predict(gam_fit)

RMSE <- function(obs,pred){

sqrt(mean((obs-pred)^2))

}

data.frame(

Model=c("LM","RF","GP","GAM"),

RMSE=c(

RMSE(dsd[[response]],pred_lm),

RMSE(dsd[[response]],pred_rf),

RMSE(dsd[[response]],pred_gp),

RMSE(dsd[[response]],pred_gam)

)

)
