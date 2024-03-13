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


## Plot based on fitted and residuals
residual_plot <- function(
        model, variable, nsigma = 1,
        point_col = 'black', horz_line = 'black',
        col1 = '#0571b0', col2 = '#ca0020',
        ...
) {
    o <- mapply(
        \(x, y) {
            # Specific for this model
            coord_data <- xy.coords(x, y)
            x <- coord_data$x; y <- coord_data$y; x0 <- sort(x);
            mod <- loess(y ~ x, ...)
            pred <- predict(mod, data.frame(x = x0), se = T)
            yfit <- pred$fit
            var <- pred$se.fit^2 
            list(x = x, y = y, yfit = yfit, var = var)
        },
        x = lapply(predict(model), \(.) .[, paste0(variable, '.predicted')]),
        y = lapply(residuals(model), \(.) .[, paste0(variable, '.residual')]),
        SIMPLIFY = F
    )
    
    # Pool
    vars <- sapply(o, \(.) .$var)
    varW <- rowMeans(vars)
    varB <- (1 / (NCOL(vars) - 1))*apply(vars, 1, var)
    sd <- sqrt(varW + varB + varB/NCOL(vars))
    
    yfit <- rowMeans(sapply(o, \(.) .$yfit))
    x <- rowMeans(sapply(o, \(.) .$x))
    y <- rowMeans(sapply(o, \(.) .$y))
    x0 <- sort(x)
    
    # Return Plot information
    return(list(
        geom_point(aes(x, y), alpha = .25, color = point_col),
        geom_hline(yintercept = 0, color = horz_line),
        geom_line(aes(x0, yfit), color = col1, linewidth = 1.1),
        geom_line(aes(x0, yfit + nsigma * sd), color = col2, linewidth = 1.1),
        geom_line(aes(x0, yfit - nsigma * sd), color = col2, linewidth = 1.1)
    ))
}

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
    ggplot()
    + residual_plot(f1, 'cwb', nsigma = 2) # Add custom create loess
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
    ggplot()
    + residual_plot(f2, 'cwb', nsigma = 2) # Add custom create loess
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
    measurement_consci = '{ consci1:consci10 } ~ consci',
    measurement_orgcon = '{ orgcon1:orgcon10 } ~ orgcon',
    measurement_consci = 'consci -> consci1@l_con consci2:consci10',
    measurement_emosta = 'emosta -> emosta1@l_emo emosta2:emosta10',
    measurement_orgcon = 'orgcon -> orgcon1@l_org orgcon2:orgcon10',
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
        'l_cwb ~ truncate(0, Inf)'
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
    ggplot()
    + residual_plot(f3, 'cwb') # Add custom create loess
    + labs(
        title = 'Residuals v. Predicted with Threeway Latent Interactions',
        subtitle = 'Averaged over 20 imputations',
        x     = 'Predicted Scores',
        y     = 'Residual Scores'
    )
)

#### Display Plots ####

# Set limits
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


#### Generate Simple Effects Plots #### 


## Extract Parameter draws for regression coefficients
mat_p <- as.matrix(f2@iterations[,c(2, 3, 4)])

# Function to Compute conditional Effects
cond_eff <- function(x, m) {
    cbind(b0 = x[,2] * m, b1 = x[,1] + x[,3] * m, m = m )
}

# Compute all conditional effects into data.frame
simple_data <- as.data.frame(rbind(
    cond_eff(mat_p, m = -1),
    cond_eff(mat_p, m =  0),
    cond_eff(mat_p, m =  1)
))

# Create factor based on levels of m
simple_data$mf <- as.factor(simple_data$m)

## Generate predicted scores

# Create range of predictor scores
xvals <- seq(-3, 3, by = 0.1)
pred_score <- \(d) d[1] + d[2] * xvals

# Compute predicted scores
pred <- lapply(split(simple_data[,1:2], simple_data$mf), \(x) apply(x, 1, pred_score))
# Compute quantiles (2.5%, 50%, 97.5%)
quan <- lapply(pred, \(x) apply(x, 1, quantile, p = c(0.025, 0.5, 0.975)))

# Combine into data.frame
rib_data  <- do.call('rbind', lapply(names(quan), \(x) {
    data.frame( l = quan[[x]][1,], fit = quan[[x]][2,], h = quan[[x]][3,], x = xvals, m = x)
}))

# Create factor with labels
rib_data$mf <- factor(
    rib_data$m,
    levels = c(1, 0, -1),
    labels = c(
        '+1 SD',
        'Mean',
        '-1 SD'
    )
)

## Make Conditional Effects Plot
cond_plot <- (
    ggplot(rib_data, aes(x, color = mf, fill = mf))
    + geom_ribbon(aes(ymin = l, ymax = h), alpha = 0.25)
    + geom_line(aes(y = fit), linewidth = 1.1)
)
    
## Print Plot with labels
(
    cond_plot 
    + scale_x_continuous(
        'Latent `consci` Scores',
        breaks = seq(-3, 3, by = 1),
        limits = c(-3, 3)
    )
    + scale_y_continuous(
        'Latent `cwb` Scores',
        breaks = seq(-1.5, 3, by = .5),
        limits = c(-1.5, 2.6)
    )
    + scale_color_brewer(palette = 'Set2')
    + scale_fill_brewer(palette = 'Set2')
    + ggtitle(
        'Plot of Conditional Regressions',
        '`orgcon` as moderator'
    )
)

#### Generate Johnson-Neyman Plot #### 

## Extract Parameter draws for regression coefficients
mat_p <- as.matrix(f2@iterations[,c(2, 3, 4)])

# Function to compute conditional effects
cond_eff <- function(m, b, func = mean, ...) {
    o <- sapply(m, \(x) b[,1] + b[,3] * x, simplify = T)
    apply(o, 2, func, ...)
}

# Function to Filter if significant or Not
set_group <- function(x){ with(rle(x), {
    unlist(lapply(seq_along(lengths), \(i) rep(i, lengths[i])))
})}

## Create plot with grouping ribbon
jn_plot <- (
    ggplot()
    # Set 0 value line
    + geom_hline(yintercept = 0, linewidth = 1.5)
    # Create Ribbon
    + stat_function(
        fun = \(m) {
            # Check if 0 is within the interval (product will be negative)
            apply(cond_eff(m, mat_p, quantile, p = c(0.025, 0.975)), 2, prod) < 0
        },
        aes(
            # Draw ribbon along lower
            ymin = cond_eff(after_stat(x), mat_p, quantile, p = 0.025),
            # Draw ribbon along upper
            ymax = cond_eff(after_stat(x), mat_p, quantile, p = 0.975),
            # Set color based on 0 being in the interval
            fill = after_stat(y), group = set_group(after_stat(y))
        ),
        # Draws a ribbon transparency and 
        geom = 'ribbon', alpha = 0.25, n = 1000
    )
    # Line for 2.5%
    + geom_function(
        fun = \(m) cond_eff(m, mat_p, quantile, p = 0.025),
        color = 'black', linetype = 'dashed', linewidth = 1.5
    )
    # Line for 97.5%
    + geom_function(
        fun = \(m) cond_eff(m, mat_p, quantile, p = 0.975),
        color = 'black', linetype = 'dashed', linewidth = 1.5
    )
    # Line for Mean
    + geom_function(
        fun = \(m) cond_eff(m, mat_p, mean),
        color = 'black', linewidth = 1.5
    )
)

## Print Plot with labels
(
    jn_plot
    + scale_x_continuous(
        'Latent `orgcon` Scores',
        breaks = seq(-3, 3, by = 1),
        limits = c(-3, 3)
    )
    + scale_y_continuous(
        'Latent `cwb` Scores | Latent `consci` Scores',
        breaks = seq(-1.5, .75, by = .25),
        limits = c(-1.25, 0.75)
    )
    + scale_fill_manual(guide = 'none', values = c('#0571b0', '#ca0020'))
    + ggtitle(
        'Johnson-Neyman Plot of Conditional Slope',
        'Red area represents 0 within 95% interval'
    )
)

#### Generate Johnson-Neyman Plot #### 

## Obtain Matrix of coefficients
mat_p3 <- as.matrix(f3@iterations[,2:8])

#### Create Simple Slope Plots For Three-way Interaction ####

# Compute conditional effects
cond_eff3 <- function(x, m1 , m2) {
    cbind(
        b0 = x[,2] * m1 + x[,3] * m2 + m1 * m2 * x[,6],
        b1 = x[,1] + x[,4] * m1 + x[,5] * m2 + x[,7] * m1 * m2,
        m1 = m1,
        m2 = m2
    )
}

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
threeway_simple_plot <- (
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

#### Generate Johnson-Neyman Plot for Two-way Interaction#### 

# Compute conditional slopes
cond_slope3 <- function(m, b, func = mean, ...) {
    o <- sapply(m, \(x) b[,4] + b[,7] * x, simplify = T)
    apply(o, 2, func, ...)
}

## Plot JN
jn_plot3 <- (
    ggplot()
    # Set 0 value line
    + geom_hline(yintercept = 0, linewidth = 1.5)
    # Create Ribbon
    + stat_function(
        fun = \(m) {
            # Check if 0 is within the interval (product will be negative)
            apply(cond_slope3(m, mat_p3, quantile, p = c(0.025, 0.975)), 2, prod) < 0
        },
        aes(
            # Draw ribbon along lower
            ymin = cond_slope3(after_stat(x), mat_p3, quantile, p = 0.025),
            # Draw ribbon along upper
            ymax = cond_slope3(after_stat(x), mat_p3, quantile, p = 0.975),
            # Set color based on 0 being in the interval
            fill = after_stat(y), group = set_group(after_stat(y))
        ),
        # Draws a ribbon transparency and 
        geom = 'ribbon', alpha = 0.25, n = 1000
    )
    # Line for 2.5%
    + geom_function(
        fun = \(m) cond_slope3(m, mat_p3, quantile, p = 0.025),
        color = 'black', linetype = 'dashed', linewidth = 1.5
    )
    # Line for 97.5%
    + geom_function(
        fun = \(m) cond_slope3(m, mat_p3, quantile, p = 0.975),
        color = 'black', linetype = 'dashed', linewidth = 1.5
    )
    # Line for Mean
    + geom_function(
        fun = \(m) cond_slope3(m, mat_p3, mean),
        color = 'black', linewidth = 1.5
    )
)

## Make pretty
(
    jn_plot3
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
    + scale_fill_manual(guide = 'none', values = c(col1, col2))
    + ggtitle(
        'Johnson-Neyman Plot for `consci` * `orgcon` Moderated by `emotsa`',
        subtitle = 'Red area represents 0 within 95% interval'
    )
)

