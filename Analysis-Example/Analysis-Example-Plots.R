#-------------------------------------------------------------------#
# A Factored Regression Approach to Modeling Latent Variable
#  Interactions and Nonlinear Effects
#
# Analysis Example in `rblimp`, including plots
#
# Analysis-Example-rblimp.R
#
# Copyright Brian T. Keller 2025 all rights reserved
#-------------------------------------------------------------------#


#### Setup ####

# Load packages
library(rblimp)
library(ggplot2)

## Set ggplot2's global theme
theme_set(
    theme_bw(base_size = 12, base_family = 'serif')
    + theme(
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.key = element_blank(),
        legend.title= element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(
            color = 'black', linewidth = 0.5,
            lineend = 'square'
        )
    )
)

#### Read data ####
mydata <- read.table(
    'data_cat.txt',
    col.names = c(
        paste0('consci', 1:10),
        paste0('emosta', 1:10),
        paste0('orgcon', 1:10),
        paste0('cwb',    1:23),
        'bin_orgcon', 'cat_orgcon'
    )
)

#### Fit the No Interaction Model ####

## Specify model
m1 <- list(
    cwb_model = c(
        'cwb ~ consci orgcon',
        'cwb ~~ cwb@1'
    ),
    predictor_model = c(
        'consci ~~ orgcon',
        'orgcon ~~ orgcon@1',
        'consci ~~ consci@1'
    ),
    measurement_consci = 'consci -> consci1@l_con consci2:consci10',
    measurement_orgcon = 'orgcon -> orgcon1@l_org orgcon2:orgcon10',
    measurement_cwb    = 'cwb -> cwb1@l_cwb cwb2:cwb23'
)

## Fit model
f1 <- rblimp(
    m1, mydata,
    ordinal = 'consci1:consci10 orgcon1:orgcon10 cwb1:cwb23',
    latent = ~ consci + orgcon + cwb,
    parameters = c(
        'l_con ~ truncate(0, Inf)',
        'l_org ~ truncate(0, Inf)',
        'l_cwb ~ truncate(0, Inf)'
    ),
    burn = 20000,
    iter = 50000,
    seed = 1239871,
    chains = 10,
    nimps = 20
)

### Create Plot of loess line with +/- 2 sd around it in red

## Residual plot
resid1 <- (
    # Use `rblimp::residual_plot` function to create plot
    residual_plot(f1, 'cwb', nsigma = 2)
    # Add custom labels
    + labs(
        title = 'Residuals v. Predicted for Model without Interaction',
        subtitle = 'Averaged over 20 imputations',
        x     = 'Predicted Scores',
        y     = 'Residual Scores'
    )
)


#### Fit the Two-way Interaction Model ####

## Specify model
m2 <- list(
    cwb_model = c(
        'cwb ~ consci orgcon consci*orgcon',
        'cwb ~~ cwb@1'
    ),
    predictor_model = c(
        'consci ~~ orgcon',
        'orgcon ~~ orgcon@1',
        'consci ~~ consci@1'
    ),
    measurement_consci = 'consci -> consci1@l_con consci2:consci10',
    measurement_orgcon = 'orgcon -> orgcon1@l_org orgcon2:orgcon10',
    measurement_cwb    = 'cwb -> cwb1@l_cwb cwb2:cwb23'
)

## Fit model
f2 <- rblimp(
    m2, mydata,
    ordinal = 'consci1:consci10 orgcon1:orgcon10 cwb1:cwb23',
    latent = ~ consci + orgcon + cwb,
    parameters = c(
        'l_con ~ truncate(0, Inf)',
        'l_org ~ truncate(0, Inf)',
        'l_cwb ~ truncate(0, Inf)'
    ),
    simple = c(
        'consci | orgcon @ -1.0',
        'consci | orgcon @  0.0',
        'consci | orgcon @  1.0'
    ),
    burn = 20000,
    iter = 50000,
    seed = 298722,
    chains = 10,
    nimps = 20
)

### Create Plot of loess line with +/- 2 sd around it in red

## Residual plot
resid2 <- (
    # Use `rblimp::residual_plot` function to create plot
    residual_plot(f2, 'cwb', nsigma = 2)
    # Add custom labels
    + labs(
        title = 'Residuals v. Predicted for Model with Interaction',
        subtitle = 'Averaged over 20 imputations',
        x     = 'Predicted Scores',
        y     = 'Residual Scores'
    )
)


#### Fit the Three-way Interaction Model ####

## Specify model
m3 <- list(
    cwb_model = c(
        'cwb ~ consci@b1 orgcon@b2 emosta@b3 
              consci*orgcon@b4 consci*emosta@b5 orgcon*emosta@b6
              consci*orgcon*emosta@b7',
        'cwb ~~ cwb@1'
    ),
    predictor_model = c(
        'consci orgcon emosta ~~ consci orgcon emosta',
        'orgcon ~~ orgcon@1',
        'consci ~~ consci@1',
        'emosta ~~ emosta@1'
    ),
    measurement_consci = 'consci -> consci1@l_con consci2:consci10',
    measurement_orgcon = 'orgcon -> orgcon1@l_org orgcon2:orgcon10',
    measurement_emosta = 'emosta -> emosta1@l_emo emosta2:emosta10',
    measurement_cwb    = 'cwb -> cwb1@l_cwb cwb2:cwb23'
)

## Fit model
f3 <- rblimp(
    m3, mydata,
    ordinal = 'consci1:consci10 orgcon1:orgcon10 cwb1:cwb23 emosta1:emosta10',
    latent = ~ consci + orgcon + cwb + emosta,
    parameters = c(
        'l_con ~ truncate(0, Inf)',
        'l_org ~ truncate(0, Inf)',
        'l_emo ~ truncate(0, Inf)',
        'l_cwb ~ truncate(0, Inf)',
        # Obtain two-way conditional effects
        'inter.emo_m1 = b4 + b7 * (-1)',
        'inter.emo_0  = b4 + b7 * ( 0)',
        'inter.emo_p1 = b4 + b7 * (+1)'
    ),
    burn = 20000,
    iter = 50000,
    seed = 548973,
    chains = 10,
    nimps = 20
)

### Create Plot of loess line with +/- 2 sd around it in red

## Residual plot
resid3 <- (
    # Use `rblimp::residual_plot` function to create plot
    residual_plot(f3, 'cwb', nsigma = 2)
    # Add custom labels
    + labs(
        title = 'Residuals v. Predicted with Threeway Latent Interactions',
        subtitle = 'Averaged over 20 imputations',
        x     = 'Predicted Scores',
        y     = 'Residual Scores'
    )
)

#### Display Plots ####

# Set custom limits
limits <- c(
    scale_x_continuous(
        breaks = c(-2, -1, 0, 1, 2, 3),
        limits = c(-2.25, 3.5)
    ),
    scale_y_continuous(
        breaks = c(-2, -1, 0, 1, 2, 3),
        limits = c(-2.5, 3.25)
    )
)

## Display Residual Plots
resid1 + limits
resid2 + limits
resid3 + limits


#### Generate Simple Effects Plots for Latent Interaction #### 

## Make Conditional Effects Plot
(
    # Create Plot with `rblimp::simple_plot` function
    simple_plot(cwb ~ consci | orgcon, f2, xvals = c(-3, 3)) 
    # Add custom labeling
    + scale_x_continuous(
        'Latent `consci` Scores',
        breaks = seq(-4, 4, by = 1)
    )
    + scale_y_continuous(
        'Latent `cwb` Scores',
        breaks = seq(-1.5, 3, by = .5)
    )
    + labs(
        title = 'Plot of Conditional Regressions',
        subtitle = '`orgcon` as moderator',
        x = 'Latent `consci` Scores',
        y = 'Latent `cwb` Scores'
    )
    # Change display colors
    + scale_color_brewer(palette = 'Set2')
    + scale_fill_brewer(palette = 'Set2')
)

#### Generate Johnson-Neyman Plot for Latent Interaction #### 

## Make Johnson-Neyman Plot

(
    # Create Plot with `rblimp::jn_plot` function
    jn_plot(cwb ~ consci | orgcon, f2)
    + scale_x_continuous(
        'Latent `orgcon` Scores',
        breaks = seq(-3, 3, by = 1),
        limits = c(-3, 3)
    )
    # Add custom labeling and limits
    + scale_y_continuous(
        'Latent `cwb` Scores | Latent `consci` Scores',
        breaks = seq(-1.5, .75, by = .25),
        limits = c(-1.25, 0.75)
    )
)

#### Create Simple Slope Plots For Three-way Interaction ####
## NOTE: These are not currently available in `rblimp`

# Function to cmopute conditional effects based on two moderators
cond_eff3 <- function(x, m1 , m2) {
    cbind(
        b0 = x[,2] * m1 + x[,3] * m2 + m1 * m2 * x[,6],
        b1 = x[,1] + x[,4] * m1 + x[,5] * m2 + x[,7] * m1 * m2,
        m1 = m1,
        m2 = m2
    )
}

## Obtain Matrix of coefficients from Focal Model
mat_p3 <- as.matrix(f3@iterations)[,2:8]

## Obtain average conditional effect
mydata3 <- as.data.frame(rbind(
    cond_eff3(mat_p3, m1 = -1, m2 = -1),
    cond_eff3(mat_p3, m1 =  0, m2 = -1),
    cond_eff3(mat_p3, m1 =  1, m2 = -1),
    cond_eff3(mat_p3, m1 = -1, m2 =  0),
    cond_eff3(mat_p3, m1 =  0, m2 =  0),
    cond_eff3(mat_p3, m1 =  1, m2 =  0),
    cond_eff3(mat_p3, m1 = -1, m2 =  1),
    cond_eff3(mat_p3, m1 =  0, m2 =  1),
    cond_eff3(mat_p3, m1 =  1, m2 =  1)
))

mydata3$m1f <- as.factor(mydata3$m1)
mydata3$m2f <- as.factor(mydata3$m2)

## Generate predicted scores

# Create xvals
xvals <- seq(-3, 3, length.out = 1000)

pred_score <- \(d) d[1] + d[2] * xvals

pred <- lapply(split(mydata3[,1:2], ~ mydata3$m1f + mydata3$m2f), \(x) apply(x, 1, pred_score))
quan <- lapply(pred, \(x) apply(x, 1, quantile, p = c(0.025, 0.5, 0.975)))
rib_data3  <- do.call(rbind, lapply(names(quan), \(x) {
    data.frame( l = quan[[x]][1,], fit = quan[[x]][2,], h = quan[[x]][3,], x = xvals, m = x)
}))

# Create overall factor
rib_data3$mf <- factor(
    rib_data3$m,
    levels = c(
        '1.1',
        '1.0',
        '1.-1',
        '0.1',
        '0.0',
        '0.-1',
        '-1.1',
        '-1.0',
        '-1.-1'
    )
)

# Split into two factors
rib_data3$m1f <- factor(
    gsub('\\..+', '', rib_data3$m),
    levels = c('-1', '0', '1'),
    labels = c(
        '-1 SD of `orgcon`',
        'Mean of `orgcon`',
        '+1 SD of `orgcon`'
    )
)
rib_data3$m2f <- factor(
    gsub('.+\\.', '', rib_data3$m),
    levels = c('-1', '0', '1'),
    labels = c(
        '-1 SD of `emosta`',
        'Mean of `emosta`',
        '+1 SD of `emosta`'
    )
)

## Make simple slope plot
(
    ggplot(rib_data3, aes(x, color = m1f, fill = m1f))
    + geom_ribbon(aes(ymin = l, ymax = h), alpha = 0.25)
    + geom_line(aes(y = fit), linewidth = 1.1)
    + facet_grid(. ~ m2f, axes = 'all')
    + scale_x_continuous(
        'Latent `consci` Scores',
        breaks = seq(-3, 3, by = 1),
        limits = c(-3, 3)
    )
    + scale_y_continuous(
        'Latent `cwb` Scores',
        breaks = seq(-2, 3.5, by = .5),
        limits = c(-2.1, 3.6)
    )
    + scale_color_brewer(palette = 'Set2', direction = -1)
    + scale_fill_brewer(palette = 'Set2', direction = -1)
    + ggtitle(
        'Plot of Conditional Regressions',
        '`orgcon` and `emosta` as moderators'
    )
)

#### Generate Johnson-Neyman Plot for the Two-way Interaction Effect #### 

# Get parameter draws from three-way model: Remove out \beta_4 and \beta_7
params_f3 <- as.matrix(f3)[, c(5, 8)]

# Create Johnson-Neyman Plot 
(
    # Create Plot with `rblimp::jn_plot_func` function
    jn_plot_func(
        # Use `rblimp::compute_condeff` function with \beta_4 and \beta_7
        compute_condeff(params_f3[,1], params_f3[,2]),
        xrange = c(-3, 3)
    )
    # Custom labels and title
    + scale_x_continuous(
        'Latent `emosta` Scores',
        breaks = seq(-3, 3, by = 1),
        limits = c(-3, 3)
    )
    + scale_y_continuous(
        '`cwb` | `consci` * `orgcon`',
        breaks = seq(-1, 0.5, by = .25),
        limits = c(-1, 0.5)
    )
    + ggtitle(
        'Johnson-Neyman Plot for `consci` * `orgcon` Moderated by `emotsa`',
        subtitle = 'Red area represents 0 within 95% interval'
    )

)
