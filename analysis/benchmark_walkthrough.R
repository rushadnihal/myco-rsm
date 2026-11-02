# ==============================================================
# benchmark_walkthrough.R -- step-by-step analysis behind the app
# Run top to bottom from the repository root (sections build on each other).
# Sections: 1 setup | 2 Block-1 lack of fit | 3 blocked 2nd-order model |
#   4 canonical analysis | 5 steepest ascent | 6 desirability |
#   7 CCD generator | 8 rotatability + BBD | 9 generic pipeline |
#   10 recovery test (simulated) | 11 setpoint card | 12 surface plots |
#   13 upload validator
# app.R contains the same functions, consolidated, plus the Shiny UI.
# ==============================================================


# ==== 1. SETUP ================================================

pkgs <- c("rsm", "plotly", "DT", "bslib", "shiny", "dplyr")
missing_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) install.packages(missing_pkgs)

suppressPackageStartupMessages({
  library(rsm)
  library(dplyr)
})

# ChemReact coding: Time (min), Temp (deg C); center = (85, 175), step = 5
CR_CENTER <- c(Time = 85, Temp = 175)
CR_STEP   <- c(Time = 5,  Temp = 5)

to_coded   <- function(X, center, step) (X - center) / step
to_natural <- function(x, center, step) center + step * x

# Label each run as Factorial / Axial / Center
label_point_type <- function(coded_df, tol = 1e-6) {
  xs <- as.matrix(coded_df)
  is_center <- apply(abs(xs) < tol, 1, all)
  is_cube   <- apply(abs(abs(xs) - 1) < tol, 1, all)
  ifelse(is_center, "Center", ifelse(is_cube, "Factorial", "Axial"))
}

load_chemreact <- function() {
  d <- as.data.frame(rsm::ChemReact)   # lazy-loaded inside rsm
  stopifnot(all(c("Time", "Temp", "Block", "Yield") %in% names(d)))
  d$x1 <- to_coded(d$Time, CR_CENTER["Time"], CR_STEP["Time"])
  d$x2 <- to_coded(d$Temp, CR_CENTER["Temp"], CR_STEP["Temp"])
  d$Point_Type <- label_point_type(d[, c("x1", "x2")])
  d
}

cr <- load_chemreact()
print(table(cr$Block, cr$Point_Type))


# ==== 2. BLOCK 1 ONLY: FO + TWI FIT AND LACK OF FIT ===========

# Formal lack-of-fit F-test: model vs one-mean-per-design-point model.
# Pure error comes from the replicated center points.
lack_of_fit <- function(formula, data, coords = c("x1", "x2")) {
  data$.cell <- factor(do.call(paste, c(data[coords], sep = "|")))
  reduced <- lm(formula, data = data)
  full    <- lm(reformulate(".cell", response = all.vars(formula)[1]),
                data = data)
  a <- anova(reduced, full)
  list(table = a, F = a$F[2], df1 = a$Df[2],
       df2 = a$Res.Df[2], p = a$`Pr(>F)`[2])
}

b1   <- subset(cr, Block == "B1")
fit1 <- rsm(Yield ~ FO(x1, x2) + TWI(x1, x2), data = b1)
print(summary(fit1))

lof1 <- lack_of_fit(Yield ~ x1 * x2, data = b1)
print(lof1$table)
cat(sprintf("\nF(%d, %d) = %.2f, p = %.4f -> %s\n",
            lof1$df1, lof1$df2, lof1$F, lof1$p,
            ifelse(lof1$p < 0.05,
                   "SIGNIFICANT curvature: add axial points (Block 2)",
                   "no significant lack of fit")))

m_fact <- mean(b1$Yield[b1$Point_Type == "Factorial"])
m_cent <- mean(b1$Yield[b1$Point_Type == "Center"])
cat(sprintf("Mean factorial = %.2f | Mean center = %.2f | curvature = %.2f\n",
            m_fact, m_cent, m_cent - m_fact))


# ==== 3. FULL BLOCKED SECOND-ORDER MODEL ======================

# Block enters as a fixed nuisance effect; SO() = FO + TWI + PQ.
fit_full <- rsm(Yield ~ Block + SO(x1, x2), data = cr)
print(summary(fit_full))   # coefficients, FO/TWI/PQ anova, lack of fit

# Block shift (B2 relative to B1), with 95% CI
blk <- summary(fit_full)$coefficients["BlockB2", ]
blk_ci <- confint(fit_full)["BlockB2", ]
cat(sprintf("\nBlock shift (B2 - B1) = %.2f  [95%% CI %.2f, %.2f], p = %.4g\n",
            blk["Estimate"], blk_ci[1], blk_ci[2], blk["Pr(>|t|)"]))

# Lack of fit of the full model. Cells are Block x design point,
# so each block's 3 centers give 2 + 2 = 4 df of pure error.
lof_full <- lack_of_fit(
  Yield ~ Block + x1 * x2 + I(x1^2) + I(x2^2),
  data   = cr,
  coords = c("Block", "x1", "x2")
)
print(lof_full$table)
cat(sprintf("\nF(%d, %d) = %.2f, p = %.4f -> %s\n",
            lof_full$df1, lof_full$df2, lof_full$F, lof_full$p,
            ifelse(lof_full$p > 0.05,
                   "second-order model is adequate",
                   "lack of fit: consider a higher-order model")))


# ==== 4. CANONICAL ANALYSIS ===================================

# Classify the surface from the eigenvalues of the quadratic matrix B.
# tol: eigenvalues with |lambda| < tol * max|lambda| count as "near zero".
classify_surface <- function(eig, tol = 0.05) {
  near0 <- abs(eig) < tol * max(abs(eig))
  if (any(near0))  return("Ridge (near-zero eigenvalue): stationary point poorly determined")
  if (all(eig < 0)) return("True maximum")
  if (all(eig > 0)) return("True minimum")
  "Saddle point"
}

# Stationary point, eigenvalues, classification, location check
canonical_report <- function(fit, coded_names, center, step, alpha) {
  can <- canonical(fit)
  xs  <- can$xs[coded_names]
  eig <- can$eigen$values
  list(
    xs_coded   = xs,
    xs_natural = to_natural(xs, center, step),
    eigenvalues = eig,
    topology   = classify_surface(eig),
    radius     = sqrt(sum(xs^2)),
    inside     = sqrt(sum(xs^2)) <= alpha + 1e-9
  )
}

alpha2 <- 2^(2/4)   # rotatable CCD, k = 2 -> sqrt(2) = 1.414
rep_cr <- canonical_report(fit_full, c("x1", "x2"),
                           center = CR_CENTER, step = CR_STEP,
                           alpha = alpha2)

cat("\nStationary point (coded):   ", round(rep_cr$xs_coded, 3), "\n")
cat("Stationary point (natural): ", round(rep_cr$xs_natural, 2),
    " (Time min, Temp C)\n")
cat("Eigenvalues:                ", round(rep_cr$eigenvalues, 2), "\n")
cat("Surface topology:           ", rep_cr$topology, "\n")
cat(sprintf("Radius from center = %.3f (alpha = %.3f) -> %s\n",
            rep_cr$radius, alpha2,
            ifelse(rep_cr$inside, "inside the design region",
                   "OUTSIDE: use steepest ascent / ridge analysis")))

# Predicted yield at the stationary point (reference block = B1)
nd <- data.frame(Block = factor("B1", levels = levels(cr$Block)),
                 x1 = rep_cr$xs_coded[1], x2 = rep_cr$xs_coded[2])
cat(sprintf("Predicted yield at optimum (Block B1 basis): %.2f\n",
            predict(fit_full, newdata = nd)))

# Same model via lm() -> predictable coefficient names
fit_chk <- lm(Yield ~ Block + x1 + x2 + I(x1^2) + I(x2^2) + x1:x2, data = cr)
cf <- coef(fit_chk)
print(round(cf, 3))

B <- matrix(c(cf["I(x1^2)"],  cf["x1:x2"] / 2,
              cf["x1:x2"] / 2, cf["I(x2^2)"]), 2, 2)
b <- cf[c("x1", "x2")]

print(B)
cat("Eigenvalues (hand):", round(eigen(B)$values, 2), "\n")
cat("Eigenvalues (rsm): ", round(rep_cr$eigenvalues, 2), "\n")
cat("Stationary point (hand):", round(-0.5 * solve(B, b), 3), "\n")
cat("Gradient at rsm xs (should be ~0):",
    round(b + 2 * B %*% rep_cr$xs_coded, 6), "\n")

# ==== 5. PATH OF STEEPEST ASCENT ==============================


# Move along the gradient of a first-order fit, from the design center.
steepest_ascent <- function(fit, coded_names, center, step,
                            dist = seq(0, 3, by = 0.5)) {
  cf <- coef(fit)
  b  <- sapply(coded_names, function(n) unname(cf[endsWith(names(cf), n)][1]))
  stopifnot(!anyNA(b))
  dir <- b / sqrt(sum(b^2))                      # unit direction (coded)
  xs  <- outer(dist, dir)                        # one row per distance
  colnames(xs) <- coded_names
  nat <- sweep(sweep(xs, 2, step, `*`), 2, center, `+`)
  colnames(nat) <- names(center)
  pred <- predict(fit, newdata = as.data.frame(xs))
  data.frame(dist = dist, round(xs, 3), round(nat, 2),
             pred_yield = round(pred, 2))
}

# First-order model on Block 1 (factorial + centers)
fit_fo <- rsm(Yield ~ FO(x1, x2), data = b1)
print(summary(fit_fo))
print(names(coef(fit_fo)))   # how rsm names its coefficients

path <- steepest_ascent(fit_fo, c("x1", "x2"),
                        center = CR_CENTER, step = CR_STEP)
print(path)

# ==== 6. DERRINGER-SUICH DESIRABILITY =========================

# Individual desirabilities (0 = unacceptable, 1 = ideal); s = shape exponent
d_max <- function(y, L, U, s = 1)           # larger is better
  ifelse(y <= L, 0, ifelse(y >= U, 1, ((y - L) / (U - L))^s))

d_min <- function(y, L, U, s = 1)           # smaller is better
  ifelse(y <= L, 1, ifelse(y >= U, 0, ((U - y) / (U - L))^s))

d_target <- function(y, L, T, U, s = 1, t = 1)   # hit a target T
  ifelse(y < L | y > U, 0,
         ifelse(y <= T, ((y - L) / (T - L))^s, ((U - y) / (U - T))^t))

# Composite D = geometric mean, weighted: D = (prod d_i^w_i)^(1 / sum w_i)
composite_D <- function(D, w = rep(1, ncol(D)))
  exp(rowSums(sweep(log(D), 2, w, `*`)) / sum(w))

# --- DEMO on ChemReact: maximize Yield, minimize an ILLUSTRATIVE cost ---
g <- seq(-alpha2, alpha2, length.out = 81)
grid <- expand.grid(x1 = g, x2 = g)
grid <- grid[sqrt(grid$x1^2 + grid$x2^2) <= alpha2, ]   # stay inside the design

grid$Time <- unname(to_natural(grid$x1, CR_CENTER["Time"], CR_STEP["Time"]))
grid$Temp <- unname(to_natural(grid$x2, CR_CENTER["Temp"], CR_STEP["Temp"]))
grid$Yield <- predict(fit_full,
                      newdata = transform(grid, Block = factor("B1", levels = levels(cr$Block))))
grid$Cost  <- 0.05 * (grid$Time - 80) + 0.03 * (grid$Temp - 170)  # made-up

D <- cbind(
  d_Yield = d_max(grid$Yield, L = 78, U = max(grid$Yield)),
  d_Cost  = d_min(grid$Cost,  L = min(grid$Cost), U = max(grid$Cost))
)
grid$D <- composite_D(D, w = c(2, 1))     # yield weighted twice as heavily

best <- grid[which.max(grid$D), ]
print(round(best, 3))

# How does the answer move with the weights?
for (w in list(c(1, 1), c(2, 1), c(1, 3))) {
  Dw <- composite_D(D, w)
  b  <- grid[which.max(Dw), ]
  cat(sprintf("w = (%g, %g): Time %.1f, Temp %.1f, Yield %.2f, Cost %.2f, D = %.3f\n",
              w[1], w[2], b$Time, b$Temp, b$Yield, b$Cost, max(Dw)))
}

# ==== 7. DESIGN GENERATOR =====================================

# Coding for the 4-factor bioprocess mode (center, step from the project doc)
coding4 <- list(
  x1 ~ (Carbon   - 30)  / 10,    # g/L   : 10, 20, 30, 40, 50
  x2 ~ (Nitrogen - 6)   / 2,     # g/L   : 2, 4, 6, 8, 10
  x3 ~ (pH       - 5.5) / 0.5,   # pH    : 4.5 ... 6.5
  x4 ~ (RPM      - 250) / 50     # RPM   : 150 ... 350
)

set.seed(42)   # reproducible randomization

# Rotatable CCD, blocked: cube block + star block, 3 centers in each
d_ccd <- ccd(basis = 4, n0 = c(3, 3), alpha = "rotatable",
             randomize = TRUE, coding = coding4)

print(nrow(d_ccd))                       # expect 30
print(table(d_ccd$Block))                # expect 19 (16 + 3) and 11 (8 + 3)
print(unique(abs(as.data.frame(d_ccd)$x1)))   # levels: 0, 1, 2 (alpha = 2)

# Run-sheet in natural units, in randomized execution order
runsheet <- decode.data(d_ccd)
print(head(runsheet, 10))

# Save for the wet lab
write.csv(as.data.frame(runsheet), "data/example_ccd4_runsheet.csv", row.names = FALSE)

# ==== 8. ROTATABILITY CHECK + BBD COMPARISON ==================

# Scaled prediction variance: N * x' (X'X)^-1 x for a 2nd-order model
spv_fun <- function(design_coded, xs) {
  f  <- ~ (x1 + x2 + x3 + x4)^2 + I(x1^2) + I(x2^2) + I(x3^2) + I(x4^2)
  X  <- model.matrix(f, as.data.frame(design_coded))
  Mi <- solve(crossprod(X))
  Z  <- model.matrix(f, as.data.frame(xs))
  nrow(X) * rowSums((Z %*% Mi) * Z)
}

# Random directions on the sphere of radius r
set.seed(1)
rand_pts <- function(r, n = 300) {
  z <- matrix(rnorm(n * 4), n, 4)
  z <- z / sqrt(rowSums(z^2)) * r
  colnames(z) <- paste0("x", 1:4)
  z
}

spv_range <- function(design, radii = c(1, 1.5)) {
  for (r in radii) {
    v <- spv_fun(design, rand_pts(r))
    cat(sprintf("  radius %.1f: SPV min %.2f, max %.2f, spread %.3f\n",
                r, min(v), max(v), max(v) - min(v)))
  }
}

cat("Rotatable CCD (30 runs):\n")
spv_range(d_ccd)

# Box-Behnken for comparison (4 factors, 3 centers)
set.seed(42)
d_bbd <- bbd(k = 4, n0 = 3, coding = coding4)
cat("\nBBD:", nrow(d_bbd), "runs\n")
print(table(d_bbd$Block))
cat("Box-Behnken:\n")
spv_range(d_bbd)

# Cube corner is outside the BBD region: compare at a fixed radius only

# ==== 9. GENERIC ANALYSIS PIPELINE ============================

analyze_rsm <- function(data, response, coded, block = NULL,
                        center, step, alpha = 2^(length(coded) / 4)) {
  data <- as.data.frame(data)
  so   <- paste0("SO(", paste(coded, collapse = ", "), ")")
  rhs  <- if (is.null(block)) so else paste(block, "+", so)
  fit  <- rsm(as.formula(paste(response, "~", rhs)), data = data)
  
  # canonical analysis
  can <- canonical(fit)
  xs  <- can$xs[coded]
  eig <- can$eigen$values
  
  # lack of fit vs pure error from replicated points
  quad <- paste0("I(", coded, "^2)", collapse = " + ")
  lin  <- paste0("(", paste(coded, collapse = " + "), ")^2")
  f_lof <- as.formula(paste(response, "~",
                            paste(c(block, lin, quad), collapse = " + ")))
  lof <- lack_of_fit(f_lof, data, coords = c(block, coded))
  
  radius <- sqrt(sum(xs^2))
  list(fit = fit,
       lof = lof,
       xs_coded    = xs,
       xs_natural  = to_natural(xs, center, step),
       eigenvalues = eig,
       topology    = classify_surface(eig),
       radius      = radius,
       alpha       = alpha,
       inside      = radius <= alpha + 1e-9)
}

print_report <- function(r) {
  cat(sprintf("Lack of fit: F(%d, %d) = %.2f, p = %.4f\n",
              r$lof$df1, r$lof$df2, r$lof$F, r$lof$p))
  cat("Stationary point (coded):  ", round(r$xs_coded, 3), "\n")
  cat("Stationary point (natural):", round(r$xs_natural, 2), "\n")
  cat("Eigenvalues:               ", round(r$eigenvalues, 2), "\n")
  cat("Topology:                  ", r$topology, "\n")
  cat(sprintf("Radius %.3f vs alpha %.3f -> %s\n", r$radius, r$alpha,
              ifelse(r$inside, "inside region", "OUTSIDE region")))
}

# Test: must reproduce the ChemReact results from sections 3 and 4
res_cr <- analyze_rsm(cr, "Yield", c("x1", "x2"), block = "Block",
                      center = CR_CENTER, step = CR_STEP)
print_report(res_cr)

# ==== 10. RECOVERY TEST ON SIMULATED 4-FACTOR DATA ============

coded4  <- paste0("x", 1:4)
center4 <- c(Carbon = 30, Nitrogen = 6, pH = 5.5, RPM = 250)
step4   <- c(Carbon = 10, Nitrogen = 2, pH = 0.5, RPM = 50)

# Known truth: Y = ymax + (x - xs)' B (x - xs)  (+ block shift + noise)
xs_true <- c(0.4, 0.3, -0.2, 0.1)
B_true  <- diag(c(-3, -2, -1.5, -2.5))
B_true[1, 2] <- B_true[2, 1] <- 0.8        # carbon x nitrogen interaction
ymax    <- 40
blk_eff <- c(`1` = 0, `2` = -1.5)
sigma   <- 0.3

simulate_ccd <- function(design, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  d  <- as.data.frame(design)
  d$Block <- factor(d$Block)               # must be a factor, not numeric
  X  <- as.matrix(d[, coded4])
  dx <- sweep(X, 2, xs_true)
  d$Y <- ymax + rowSums((dx %*% B_true) * dx) +
    blk_eff[as.character(d$Block)] + rnorm(nrow(d), 0, sigma)
  d
}

# --- one run, full report ---
sim1 <- simulate_ccd(d_ccd, seed = 1)
res1 <- analyze_rsm(sim1, "Y", coded4, block = "Block",
                    center = center4, step = step4)
print_report(res1)
cat("\nTRUE stationary point:", xs_true, "\n")
cat("TRUE eigenvalues:     ", round(eigen(B_true)$values, 2), "\n")

# --- 200 repeated experiments: bias and false-alarm rate ---
n_sim <- 200
est   <- matrix(NA, n_sim, 4)
lof_p <- numeric(n_sim)
for (i in seq_len(n_sim)) {
  r <- analyze_rsm(simulate_ccd(d_ccd, seed = 1000 + i), "Y", coded4,
                   block = "Block", center = center4, step = step4)
  est[i, ]  <- r$xs_coded
  lof_p[i]  <- r$lof$p
}
cat("\nMean estimated xs:", round(colMeans(est), 3), "\n")
cat("SD of estimated xs:", round(apply(est, 2, sd), 3), "\n")
cat(sprintf("Lack-of-fit p > 0.05 in %.1f%% of runs (expect ~95%%)\n",
            100 * mean(lof_p > 0.05)))

# ==== 11. MULTI-OBJECTIVE OPTIMIZER + SETPOINT CARD ===========
# (d_max, d_min, composite_D come from section 6; run that section first)

# --- simulate 2 responses on d_ccd (illustrative truth) ---
sim_surface <- function(design, xs, B, ymax, blk, sigma, seed) {
  set.seed(seed)
  d  <- as.data.frame(design)
  d$Block <- factor(d$Block)
  dx <- sweep(as.matrix(d[, coded4]), 2, xs)
  ymax + rowSums((dx %*% B) * dx) + blk[as.character(d$Block)] +
    rnorm(nrow(d), 0, sigma)
}

sim_data <- as.data.frame(d_ccd)
sim_data$Block <- factor(sim_data$Block)
sim_data$Y1 <- sim_surface(d_ccd, xs_true, B_true, 40,            # biomass g/L
                           c(`1` = 0, `2` = -1.5), 0.3, seed = 11)
sim_data$Y2 <- sim_surface(d_ccd, c(-0.2, 0.5, 0.1, -0.3),        # protein %
                           diag(c(-2, -3, -1, -1.5)), 52,
                           c(`1` = 0, `2` = 0.5), 0.8, seed = 12)

fits <- list(
  Y1 = analyze_rsm(sim_data, "Y1", coded4, "Block", center4, step4)$fit,
  Y2 = analyze_rsm(sim_data, "Y2", coded4, "Block", center4, step4)$fit
)

# --- illustrative deterministic cost (EUR/L) from natural units ---
media_cost <- function(nat) 0.00035 * nat[1] + 0.0015 * nat[2]

# --- composite desirability at a coded point ---
alpha4 <- 2
w      <- c(2, 1, 1)           # biomass, protein, cost
newd   <- function(x) {
  nd <- as.data.frame(t(setNames(x, coded4)))
  nd$Block <- factor("1", levels = levels(sim_data$Block))   # reference block
  nd
}
D_at <- function(x) {
  nd  <- newd(x)
  y1  <- predict(fits$Y1, nd)
  y2  <- predict(fits$Y2, nd)
  y3  <- media_cost(to_natural(x, center4, step4))
  d   <- c(d_max(y1, 30, 41), d_max(y2, 35, 50), d_min(y3, 0.0065, 0.0325))
  prod(d^w)^(1 / sum(w))
}

# --- multi-start Nelder-Mead, penalized to stay inside radius alpha ---
obj <- function(x) -D_at(x) + 10 * max(0, sqrt(sum(x^2)) - alpha4)
set.seed(7)
starts <- matrix(runif(40 * 4, -1, 1), 40, 4)
runs <- lapply(seq_len(nrow(starts)),
               function(i) optim(starts[i, ], obj, method = "Nelder-Mead"))
best <- runs[[which.min(sapply(runs, `[[`, "value"))]]
x_opt <- best$par

# --- Optimal Setpoint Card ---
nat <- to_natural(x_opt, center4, step4)
nd  <- newd(x_opt)
p1  <- predict(fits$Y1, nd, interval = "prediction")
p2  <- predict(fits$Y2, nd, interval = "prediction")

cat("\n===== OPTIMAL BIOREACTOR SETPOINT CARD (simulated) =====\n")
cat(sprintf("Carbon   %.1f g/L\nNitrogen %.2f g/L\npH       %.2f\nAgitation %.0f RPM\n",
            nat[1], nat[2], nat[3], nat[4]))
cat(sprintf("Biomass  %.2f g/L  [95%% PI %.2f, %.2f]\n", p1[1], p1[2], p1[3]))
cat(sprintf("Protein  %.1f %%    [95%% PI %.1f, %.1f]\n", p2[1], p2[2], p2[3]))
cat(sprintf("Cost     %.4f EUR/L\n", media_cost(nat)))
cat(sprintf("Composite D = %.3f | radius %.2f (alpha %.1f)\n",
            D_at(x_opt), sqrt(sum(x_opt^2)), alpha4))

# --- does the answer depend on the start? ---
vals <- sapply(runs, function(r) -r$value)
cat(sprintf("Best D across 40 starts: %.4f | starts within 0.001 of best: %d\n",
            max(vals), sum(vals > max(vals) - 0.001)))

# ==== 12. RESPONSE SURFACE PLOTS ==============================
library(plotly)

# Predicted surface over two coded factors; others held at `hold` (coded)
surface_grid <- function(fit, coded, vars, hold, center, step,
                         block = NULL, block_levels = NULL,
                         alpha = 2^(length(coded) / 4), n = 60) {
  g    <- seq(-alpha, alpha, length.out = n)
  grid <- expand.grid(g, g)
  names(grid) <- vars
  for (o in setdiff(coded, vars)) grid[[o]] <- hold[[o]]
  if (!is.null(block)) grid[[block]] <- factor(block_levels[1], levels = block_levels)
  z <- predict(fit, grid)
  rad <- sqrt(rowSums(as.matrix(grid[, coded])^2))
  z[rad > alpha + 1e-9] <- NA                       # mask extrapolation
  idx <- match(vars, coded)
  list(x = center[idx[1]] + step[idx[1]] * g,       # natural units
       y = center[idx[2]] + step[idx[2]] * g,
       z = t(matrix(z, n, n)),                      # rows = y, cols = x
       labs = names(center)[idx])
}

plot_surface <- function(sg, type = c("contour", "surface"), mark = NULL,
                         zlab = "Response") {
  type <- match.arg(type)
  p <- plot_ly(x = sg$x, y = sg$y, z = sg$z, type = type) |>
    layout(xaxis = list(title = sg$labs[1]), yaxis = list(title = sg$labs[2]),
           title = zlab)
  if (type == "surface")
    p <- p |> layout(scene = list(xaxis = list(title = sg$labs[1]),
                                  yaxis = list(title = sg$labs[2]),
                                  zaxis = list(title = zlab)))
  if (type == "contour" && !is.null(mark))
    p <- p |> add_markers(x = mark[1], y = mark[2], inherit = FALSE,
                          marker = list(color = "red", size = 10, symbol = "x"),
                          name = "optimum")
  p
}

# --- 2-factor benchmark: yield surface, Block B1 basis, optimum marked ---
sg_cr <- surface_grid(res_cr$fit, c("x1", "x2"), c("x1", "x2"),
                      hold = list(), center = CR_CENTER, step = CR_STEP,
                      block = "Block", block_levels = levels(cr$Block))
plot_surface(sg_cr, "contour", mark = res_cr$xs_natural, zlab = "Yield (%)")
plot_surface(sg_cr, "surface", zlab = "Yield (%)")

# --- 4-factor: biomass over Carbon x Nitrogen, pH and RPM at the optimum ---
hold <- as.list(setNames(x_opt, coded4))
sg4 <- surface_grid(fits$Y1, coded4, c("x1", "x2"), hold = hold,
                    center = center4, step = step4,
                    block = "Block", block_levels = levels(sim_data$Block))
plot_surface(sg4, "contour", mark = nat[1:2], zlab = "Biomass (g/L)")

# ==== 13. UPLOAD VALIDATOR ====================================

# Infer center & step from a factor's levels (CCD: 5 levels, BBD: 3 levels)
infer_coding <- function(x) {
  lv <- sort(unique(round(x, 8)))
  n  <- length(lv)
  if (!n %in% c(3, 5))
    stop(sprintf("expected 3 (BBD) or 5 (CCD) distinct levels, found %d", n))
  mid  <- lv[(n + 1) / 2]
  step <- lv[(n + 1) / 2 + 1] - mid
  sym  <- abs((mid - lv[1]) - (lv[n] - mid)) < 1e-6 * max(1, abs(mid))
  list(center = mid, step = step, alpha = (lv[n] - mid) / step, symmetric = sym)
}

validate_design <- function(df, factors, responses, block = NULL) {
  msgs <- character()
  err  <- function(m) msgs <<- c(msgs, paste("ERROR:", m))
  wrn  <- function(m) msgs <<- c(msgs, paste("WARNING:", m))
  
  need <- c(factors, responses, block)
  miss <- setdiff(need, names(df))
  if (length(miss)) { err(paste("missing columns:", paste(miss, collapse = ", ")))
    return(list(ok = FALSE, messages = msgs)) }
  
  for (v in c(factors, responses))
    if (!is.numeric(df[[v]])) err(paste(v, "is not numeric"))
  if (anyNA(df[, need])) err("missing values (NA) in the design or responses")
  if (any(msgs != "")) if (any(grepl("^ERROR", msgs)))
    return(list(ok = FALSE, messages = msgs))
  
  cod <- lapply(df[factors], infer_coding)
  center <- sapply(cod, `[[`, "center"); step <- sapply(cod, `[[`, "step")
  alpha  <- sapply(cod, `[[`, "alpha")
  if (!all(sapply(cod, `[[`, "symmetric"))) wrn("some factor levels are not symmetric about the center")
  
  coded <- as.data.frame(mapply(function(x, c0, s) (x - c0) / s,
                                df[factors], center, step, SIMPLIFY = FALSE))
  names(coded) <- paste0("x", seq_along(factors))
  is_ctr <- rowSums(abs(as.matrix(coded)) < 1e-6) == length(factors)
  
  k <- length(factors)
  p <- 1 + k + choose(k, 2) + k + if (is.null(block)) 0 else
    (length(unique(df[[block]])) - 1)
  if (nrow(df) <= p) err(sprintf("%d runs but the model has %d parameters", nrow(df), p))
  if (sum(is_ctr) < 2) err("fewer than 2 replicated center points: no pure error, no lack-of-fit test")
  if (!is.null(block)) {
    ctr_by_blk <- table(factor(df[[block]])[is_ctr])
    if (any(ctr_by_blk < 2)) wrn("a block has fewer than 2 center points")
  }
  if (diff(range(alpha)) > 1e-6) wrn("axial distances differ between factors")
  
  ok <- !any(grepl("^ERROR", msgs))
  list(ok = ok, messages = msgs, center = center, step = step, alpha = alpha,
       n_runs = nrow(df), n_center = sum(is_ctr),
       data = if (ok) cbind(df, coded) else NULL)
}

# --- test 1: the 30-run design (should pass and recover center4/step4) ---
rs <- as.data.frame(decode.data(d_ccd))
rs$Y1 <- sim_data$Y1
v <- validate_design(rs, c("Carbon", "Nitrogen", "pH", "RPM"), "Y1", "Block")
print(v[c("ok", "messages", "center", "step", "alpha", "n_runs", "n_center")])

# --- test 2: drop the center points (should fail) ---
v2 <- validate_design(subset(rs, !(Carbon == 30 & Nitrogen == 6 & pH == 5.5 & RPM == 250)),
                      c("Carbon", "Nitrogen", "pH", "RPM"), "Y1", "Block")
print(v2$messages)

