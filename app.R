# ==============================================================
# myco-rsm: Response Surface Methodology & multi-objective media
# optimizer for submerged mycelial fermentation  (single-file Shiny app)
#
# Sections: 1. helpers and statistical functions
#           2. app-only helpers
#           3. UI
#           4. server
# Run locally:  shiny::runApp()
# ==============================================================
library(shiny)
library(bslib)
library(rsm)
library(plotly)
library(DT)

# ---- 1. helpers and statistical functions ---------------------------
to_coded   <- function(X, center, step) (X - center) / step
to_natural <- function(x, center, step) center + step * x

CR_CENTER <- c(Time = 85, Temp = 175)
CR_STEP   <- c(Time = 5,  Temp = 5)

label_point_type <- function(coded_df, tol = 1e-6) {
  xs <- as.matrix(coded_df)
  is_center <- apply(abs(xs) < tol, 1, all)
  is_cube   <- apply(abs(abs(xs) - 1) < tol, 1, all)
  ifelse(is_center, "Center", ifelse(is_cube, "Factorial", "Axial"))
}

load_chemreact <- function() {
  d <- as.data.frame(rsm::ChemReact)
  d$x1 <- to_coded(d$Time, CR_CENTER["Time"], CR_STEP["Time"])
  d$x2 <- to_coded(d$Temp, CR_CENTER["Temp"], CR_STEP["Temp"])
  d$Point_Type <- label_point_type(d[, c("x1", "x2")])
  d
}

# ---- lack of fit + surface classification ----------------------
lack_of_fit <- function(formula, data, coords) {
  data$.cell <- factor(do.call(paste, c(data[coords], sep = "|")))
  reduced <- lm(formula, data = data)
  full    <- lm(reformulate(".cell", response = all.vars(formula)[1]),
                data = data)
  a <- anova(reduced, full)
  list(table = a, F = a$F[2], df1 = a$Df[2],
       df2 = a$Res.Df[2], p = a$`Pr(>F)`[2])
}

classify_surface <- function(eig, tol = 0.05) {
  near0 <- abs(eig) < tol * max(abs(eig))
  if (any(near0))   return("Ridge (near-zero eigenvalue): stationary point poorly determined")
  if (all(eig < 0)) return("True maximum")
  if (all(eig > 0)) return("True minimum")
  "Saddle point"
}

# ---- generic blocked second-order analysis ---------------------
analyze_rsm <- function(data, response, coded, block = NULL,
                        center, step, alpha = 2^(length(coded) / 4)) {
  data <- as.data.frame(data)
  so   <- paste0("SO(", paste(coded, collapse = ", "), ")")
  rhs  <- if (is.null(block)) so else paste(block, "+", so)
  fit  <- rsm(as.formula(paste(response, "~", rhs)), data = data)

  can <- canonical(fit)
  xs  <- can$xs[coded]
  eig <- can$eigen$values

  quad  <- paste0("I(", coded, "^2)", collapse = " + ")
  lin   <- paste0("(", paste(coded, collapse = " + "), ")^2")
  f_lof <- as.formula(paste(response, "~",
                            paste(c(block, lin, quad), collapse = " + ")))
  lof <- lack_of_fit(f_lof, data, coords = c(block, coded))

  radius <- sqrt(sum(xs^2))
  list(fit = fit, lof = lof, xs_coded = xs,
       xs_natural = to_natural(xs, center, step),
       eigenvalues = eig, topology = classify_surface(eig),
       radius = radius, alpha = alpha, inside = radius <= alpha + 1e-9)
}

# ---- desirability ----------------------------------------------
d_max <- function(y, L, U, s = 1)
  ifelse(y <= L, 0, ifelse(y >= U, 1, ((y - L) / (U - L))^s))
d_min <- function(y, L, U, s = 1)
  ifelse(y <= L, 1, ifelse(y >= U, 0, ((U - y) / (U - L))^s))
d_target <- function(y, L, T, U, s = 1, t = 1)
  ifelse(y < L | y > U, 0,
         ifelse(y <= T, ((y - L) / (T - L))^s, ((U - y) / (U - T))^t))
composite_D <- function(D, w = rep(1, ncol(D)))
  exp(rowSums(sweep(log(D), 2, w, `*`)) / sum(w))

# ---- surface plots ---------------------------------------------
surface_grid <- function(fit, coded, vars, hold, center, step,
                         block = NULL, block_levels = NULL,
                         alpha = 2^(length(coded) / 4), n = 60) {
  g    <- seq(-alpha, alpha, length.out = n)
  grid <- expand.grid(g, g)
  names(grid) <- vars
  for (o in setdiff(coded, vars)) grid[[o]] <- hold[[o]]
  if (!is.null(block)) grid[[block]] <- factor(block_levels[1], levels = block_levels)
  z   <- predict(fit, grid)
  rad <- sqrt(rowSums(as.matrix(grid[, coded])^2))
  z[rad > alpha + 1e-9] <- NA
  idx <- match(vars, coded)
  list(x = center[idx[1]] + step[idx[1]] * g,
       y = center[idx[2]] + step[idx[2]] * g,
       z = t(matrix(z, n, n)),
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

# ---- upload validator ------------------------------------------
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
  if (length(miss)) {
    err(paste("missing columns:", paste(miss, collapse = ", ")))
    return(list(ok = FALSE, messages = msgs))
  }
  for (v in c(factors, responses))
    if (!is.numeric(df[[v]])) err(paste(v, "is not numeric"))
  if (anyNA(df[, need])) err("missing values (NA) in the design or responses")
  if (any(grepl("^ERROR", msgs))) return(list(ok = FALSE, messages = msgs))

  cod    <- lapply(df[factors], function(x) tryCatch(infer_coding(x),
                                                     error = function(e) e))
  bad    <- sapply(cod, inherits, "error")
  if (any(bad)) {
    err(paste(names(cod)[bad], ":", sapply(cod[bad], conditionMessage),
              collapse = "; "))
    return(list(ok = FALSE, messages = msgs))
  }
  center <- sapply(cod, `[[`, "center"); step <- sapply(cod, `[[`, "step")
  alpha  <- sapply(cod, `[[`, "alpha")
  if (!all(sapply(cod, `[[`, "symmetric")))
    wrn("some factor levels are not symmetric about the center")

  coded <- as.data.frame(mapply(function(x, c0, s) (x - c0) / s,
                                df[factors], center, step, SIMPLIFY = FALSE))
  names(coded) <- paste0("x", seq_along(factors))
  is_ctr <- rowSums(abs(as.matrix(coded)) < 1e-6) == length(factors)

  k <- length(factors)
  p <- 1 + k + choose(k, 2) + k +
       if (is.null(block)) 0 else (length(unique(df[[block]])) - 1)
  if (nrow(df) <= p) err(sprintf("%d runs but the model has %d parameters", nrow(df), p))
  if (sum(is_ctr) < 2) err("fewer than 2 replicated center points: no pure error, no lack-of-fit test")
  if (!is.null(block) && any(table(factor(df[[block]])[is_ctr]) < 2))
    wrn("a block has fewer than 2 center points")
  if (diff(range(alpha)) > 1e-6) wrn("axial distances differ between factors")

  ok <- !any(grepl("^ERROR", msgs))
  list(ok = ok, messages = msgs, center = center, step = step, alpha = alpha,
       n_runs = nrow(df), n_center = sum(is_ctr),
       data = if (ok) cbind(df, coded) else NULL)
}

# ---- simulated 4-factor "upload" for the demo dataset ----------
make_demo4 <- function() {
  coding4 <- list(x1 ~ (Carbon - 30) / 10, x2 ~ (Nitrogen - 6) / 2,
                  x3 ~ (pH - 5.5) / 0.5,   x4 ~ (RPM - 250) / 50)
  set.seed(42)
  d  <- ccd(basis = 4, n0 = c(3, 3), alpha = "rotatable",
            randomize = TRUE, coding = coding4)
  co <- as.data.frame(d)
  X  <- as.matrix(co[, paste0("x", 1:4)])

  surf <- function(xs, B, ymax, blk, sigma, seed) {
    set.seed(seed)
    dx <- sweep(X, 2, xs)
    ymax + rowSums((dx %*% B) * dx) + blk[as.character(co$Block)] +
      rnorm(nrow(X), 0, sigma)
  }
  B1 <- diag(c(-3, -2, -1.5, -2.5)); B1[1, 2] <- B1[2, 1] <- 0.8
  nat <- as.data.frame(decode.data(d))
  out <- nat[, c("Carbon", "Nitrogen", "pH", "RPM")]
  out$Block   <- factor(co$Block)
  out$Biomass <- surf(c(0.4, 0.3, -0.2, 0.1), B1, 40, c(`1` = 0, `2` = -1.5), 0.3, 11)
  out$Protein <- surf(c(-0.2, 0.5, 0.1, -0.3), diag(c(-2, -3, -1, -1.5)),
                      52, c(`1` = 0, `2` = 0.5), 0.8, 12)
  out
}
# ---- 2. app-only helpers ---------------------------------------------

# Steepest-ascent path from a first-order fit on the cube + center runs
steepest_path <- function(data, response, coded, block, center, step,
                          dist = seq(0, 3, by = 0.5)) {
  X    <- as.matrix(data[, coded])
  cube <- data[rowSums(abs(X) > 1 + 1e-6) == 0, ]
  use_blk <- !is.null(block) && length(unique(cube[[block]])) > 1
  rhs <- paste0("FO(", paste(coded, collapse = ", "), ")")
  if (use_blk) rhs <- paste(block, "+", rhs)
  fit <- rsm(as.formula(paste(response, "~", rhs)), data = cube)

  cf <- coef(fit)
  b  <- sapply(coded, function(n) unname(cf[endsWith(names(cf), n)][1]))
  if (anyNA(b)) stop("first-order coefficients could not be estimated")
  dir <- b / sqrt(sum(b^2))
  xs  <- outer(dist, dir); colnames(xs) <- coded
  nat <- sweep(sweep(xs, 2, step, `*`), 2, center, `+`)
  colnames(nat) <- names(center)

  nd <- as.data.frame(xs)
  if (use_blk) nd[[block]] <- factor(levels(cube[[block]])[1],
                                     levels = levels(cube[[block]]))
  data.frame(dist = dist, round(xs, 3), round(nat, 2),
             predicted = round(predict(fit, nd), 2), check.names = FALSE)
}

block_levels <- function(d) if (is.null(d$block)) NULL else levels(d$data[[d$block]])

# One-row new-data frame at coded point x, in the reference (first) block
ref_newd <- function(d, x) {
  nd <- as.data.frame(t(setNames(as.numeric(x), d$coded)))
  if (!is.null(d$block))
    nd[[d$block]] <- factor(block_levels(d)[1], levels = block_levels(d))
  nd
}

# ==============================================================
# UI
# ==============================================================
ui <- page_navbar(
  title = "myco-rsm",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(
    width = 320,
    radioButtons("src", "Dataset",
                 choices = c("ChemReact benchmark (2 factors)" = "chem",
                             "Simulated 4-factor demo"         = "demo",
                             "Upload CSV"                      = "upload")),
    conditionalPanel(
      "input.src == 'upload'",
      fileInput("file", "CSV file", accept = ".csv"),
      uiOutput("map_ui")
    ),
    uiOutput("resp_ui"),
    hr(),
    helpText("Models are second-order (SO) with a fixed Block effect.",
             "Predictions use the first block as reference.")
  ),

  nav_panel("Data",
    card(card_header("Design check"),
         verbatimTextOutput("design_info"), verbatimTextOutput("msgs")),
    card(card_header("Data"), DTOutput("tbl"))
  ),

  nav_panel("Model",
    card(card_header("Lack-of-fit test (vs pure error from replicated centers)"),
         verbatimTextOutput("lof_txt")),
    card(card_header("Second-order fit (rsm summary)"),
         verbatimTextOutput("fit_txt"))
  ),

  nav_panel("Optimum",
    card(card_header("Stationary point"), tableOutput("opt_tbl")),
    card(card_header("Surface type"), verbatimTextOutput("eig_txt")),
    card(card_header("Path of steepest ascent"),
         checkboxInput("show_sa", "Show path (first-order fit on cube + center runs)", FALSE),
         tableOutput("sa_tbl"))
  ),

  nav_panel("Surfaces",
    card(uiOutput("surf_ui"), plotlyOutput("surf", height = "560px"))
  ),

  nav_panel("Optimizer",
    card(card_header("Desirability settings"), uiOutput("opt_ui")),
    card(card_header("Optimal setpoint"),
         tableOutput("opt_set"), tableOutput("opt_pred"),
         verbatimTextOutput("opt_txt"))
  ),

  nav_panel("Design",
    layout_columns(
      col_widths = c(4, 8),
      card(card_header("Run-sheet generator"),
           radioButtons("d_type", "Design", c("Rotatable CCD" = "CCD",
                                              "Box-Behnken"   = "BBD"), inline = TRUE),
           numericInput("d_k", "Number of factors (CCD 2-5, BBD 3-5)", 4, 2, 5),
           numericInput("d_n0", "Center points (per block for CCD)", 3, 1, 10),
           textInput("d_names",  "Factor names (comma-separated, no spaces)",
                     "Carbon,Nitrogen,pH,RPM"),
           textInput("d_center", "Center values", "30,6,5.5,250"),
           textInput("d_step",   "Step sizes",    "10,2,0.5,50"),
           numericInput("d_seed", "Random seed", 42),
           actionButton("gen", "Generate", class = "btn-primary"),
           downloadButton("dl", "Download CSV")),
      card(card_header("Run-sheet (natural units, randomized)"),
           verbatimTextOutput("d_info"), DTOutput("d_tbl"))
    )
  )
)

# ==============================================================
# SERVER
# ==============================================================
server <- function(input, output, session) {

  # ---- data sources ---------------------------------------------
  demo_data <- reactive(make_demo4())

  uploaded <- reactive({
    req(input$file)
    read.csv(input$file$datapath)
  })

  output$map_ui <- renderUI({
    df  <- uploaded()
    num <- names(df)[sapply(df, is.numeric)]
    tagList(
      selectizeInput("up_factors", "Factors (2 or more)", choices = num, multiple = TRUE),
      selectizeInput("up_resp", "Responses", choices = num, multiple = TRUE),
      selectInput("up_block", "Block column", choices = c("(none)", names(df)))
    )
  })

  ds <- reactive({
    if (input$src == "chem") {
      d <- load_chemreact()
      coded <- c("x1", "x2")
      return(list(ok = TRUE, data = d, factors = c("Time", "Temp"), coded = coded,
                  resp = "Yield", block = "Block",
                  center = CR_CENTER, step = CR_STEP,
                  rmax = max(sqrt(rowSums(as.matrix(d[, coded])^2))),
                  msgs = character()))
    }
    if (input$src == "demo") {
      raw     <- demo_data()
      factors <- c("Carbon", "Nitrogen", "pH", "RPM")
      resp    <- c("Biomass", "Protein")
      block   <- "Block"
    } else {
      raw     <- uploaded()
      factors <- input$up_factors
      resp    <- input$up_resp
      block   <- if (is.null(input$up_block) || input$up_block == "(none)") NULL
                 else input$up_block
      req(length(factors) >= 2, length(resp) >= 1)
    }
    v <- validate_design(raw, factors, resp, block)
    if (!v$ok) return(list(ok = FALSE, msgs = v$messages))
    d <- v$data
    coded <- paste0("x", seq_along(factors))
    if (!is.null(block)) d[[block]] <- factor(d[[block]])
    list(ok = TRUE, data = d, factors = factors, coded = coded,
         resp = resp, block = block, center = v$center, step = v$step,
         rmax = max(sqrt(rowSums(as.matrix(d[, coded])^2))), msgs = v$messages)
  })

  models <- reactive({
    d <- ds(); req(d$ok)
    out <- tryCatch(
      setNames(lapply(d$resp, function(r)
        analyze_rsm(d$data, r, d$coded, d$block, d$center, d$step, d$rmax)),
        d$resp),
      error = function(e) {
        showNotification(paste("Model error:", conditionMessage(e)),
                         type = "error", duration = NULL)
        NULL
      })
    req(out)
    out
  })

  cur <- reactive({ req(input$resp); models()[[input$resp]] })

  output$resp_ui <- renderUI({
    d <- ds(); req(d$ok)
    selectInput("resp", "Response to display", d$resp)
  })

  # ---- Data tab ---------------------------------------------------
  output$design_info <- renderText({
    d <- ds()
    if (!d$ok) return("Design failed validation.")
    X <- as.matrix(d$data[, d$coded])
    sprintf("%d runs | %d factors | %d center runs | design region radius %.3f (coded)",
            nrow(X), ncol(X), sum(rowSums(abs(X) < 1e-6) == ncol(X)), d$rmax)
  })
  output$msgs <- renderText({
    d <- ds()
    if (length(d$msgs)) paste(d$msgs, collapse = "\n") else "Design checks passed."
  })
  output$tbl <- renderDT({
    d <- ds(); req(d$ok)
    datatable(d$data, options = list(pageLength = 15, scrollX = TRUE))
  })

  # ---- Model tab --------------------------------------------------
  output$lof_txt <- renderText({
    m <- cur()
    sprintf("F(%d, %d) = %.2f, p = %.4f -> %s",
            m$lof$df1, m$lof$df2, m$lof$F, m$lof$p,
            ifelse(m$lof$p > 0.05, "second-order model is adequate",
                   "lack of fit: consider a higher-order model or more runs"))
  })
  output$fit_txt <- renderPrint({ summary(cur()$fit) })

  # ---- Optimum tab ------------------------------------------------
  output$opt_tbl <- renderTable({
    d <- ds(); m <- cur()
    data.frame(Factor = d$factors,
               Coded = round(as.numeric(m$xs_coded), 3),
               Natural = round(as.numeric(m$xs_natural), 3))
  })

  output$eig_txt <- renderText({
    d <- ds(); m <- cur()
    pred <- predict(m$fit, ref_newd(d, m$xs_coded))
    txt <- c(
      sprintf("Eigenvalues: %s", paste(round(m$eigenvalues, 2), collapse = ", ")),
      sprintf("Topology:    %s", m$topology),
      sprintf("Radius %.3f vs region %.3f -> %s", m$radius, d$rmax,
              ifelse(m$inside, "inside the experimental region",
                     "OUTSIDE the region (extrapolation): use the steepest-ascent path")),
      sprintf("Predicted %s at stationary point (reference block): %.2f",
              input$resp, pred)
    )
    if (input$src == "chem" && input$resp == "Yield") {
      ok <- all(abs(as.numeric(m$xs_coded) - c(0.372, 0.334)) < 0.002) &&
            all(abs(round(m$eigenvalues, 2) - c(-0.92, -1.32)) < 0.011)
      txt <- c(txt, "",
               sprintf("Benchmark check (0.372, 0.334; eigenvalues -0.92, -1.32): %s",
                       ifelse(ok, "PASS", "FAIL")))
    }
    paste(txt, collapse = "\n")
  })

  output$sa_tbl <- renderTable({
    req(input$show_sa)
    d <- ds()
    tryCatch(steepest_path(d$data, input$resp, d$coded, d$block, d$center, d$step),
             error = function(e) data.frame(Note = conditionMessage(e)))
  })

  # ---- Surfaces tab -----------------------------------------------
  output$surf_ui <- renderUI({
    d <- ds(); req(d$ok)
    fluidRow(
      column(4, selectInput("sv1", "X axis", d$factors, selected = d$factors[1])),
      column(4, selectInput("sv2", "Y axis", d$factors, selected = d$factors[2])),
      column(4, radioButtons("stype", "Plot", c("Contour" = "contour",
                                                "3D surface" = "surface"),
                             inline = TRUE))
    )
  })

  output$surf <- renderPlotly({
    d <- ds(); m <- cur()
    req(input$sv1, input$sv2, input$stype)
    validate(need(input$sv1 != input$sv2, "Pick two different factors."))
    sel  <- c(input$sv1, input$sv2)
    vars <- d$coded[match(sel, d$factors)]
    hold <- if (m$inside) as.list(m$xs_coded)
            else as.list(setNames(rep(0, length(d$coded)), d$coded))
    sg <- surface_grid(m$fit, d$coded, vars, hold, d$center, d$step,
                       block = d$block, block_levels = block_levels(d),
                       alpha = d$rmax)
    mark <- if (m$inside) as.numeric(m$xs_natural)[match(sel, d$factors)] else NULL
    plot_surface(sg, input$stype, mark = mark, zlab = input$resp)
  })

  # ---- Optimizer tab ----------------------------------------------
  output$opt_ui <- renderUI({
    d <- ds(); req(d$ok)
    tagList(
      lapply(seq_along(d$resp), function(i) {
        rg <- range(d$data[[d$resp[i]]])
        fluidRow(
          column(12, strong(d$resp[i])),
          column(3, selectInput(paste0("goal_", i), "Goal",
                                c("Maximize", "Minimize", "Target"))),
          column(2, numericInput(paste0("L_", i), "Low",    signif(rg[1], 4))),
          column(2, numericInput(paste0("T_", i), "Target", signif(mean(rg), 4))),
          column(2, numericInput(paste0("U_", i), "High",   signif(rg[2], 4))),
          column(3, numericInput(paste0("w_", i), "Weight", 1, min = 0.1, step = 0.5))
        )
      }),
      hr(),
      checkboxInput("use_cost",
                    "Include media cost (assumes factors 1 and 2 are concentrations in g/L)",
                    FALSE),
      conditionalPanel("input.use_cost",
        fluidRow(
          column(4, numericInput("p1", "Price, factor 1 (EUR/kg)", 0.35)),
          column(4, numericInput("p2", "Price, factor 2 (EUR/kg)", 1.5)),
          column(4, numericInput("wc", "Cost weight", 1, min = 0.1, step = 0.5))
        )),
      actionButton("run_opt", "Optimize", class = "btn-primary")
    )
  })

  opt <- reactiveVal(NULL)

  observeEvent(input$run_opt, {
    d <- ds(); mods <- models()
    nR <- length(d$resp); k <- length(d$coded)
    spec <- lapply(seq_len(nR), function(i) list(
      goal = input[[paste0("goal_", i)]], L = input[[paste0("L_", i)]],
      T = input[[paste0("T_", i)]], U = input[[paste0("U_", i)]],
      w = input[[paste0("w_", i)]]))
    req(!is.null(spec[[1]]$goal))

    use_cost <- isTRUE(input$use_cost)
    cost_at  <- function(x) {
      nat <- to_natural(as.numeric(x), as.numeric(d$center), as.numeric(d$step))
      (input$p1 * nat[1] + input$p2 * nat[2]) / 1000
    }
    lo <- pmax(d$center - d$rmax * d$step, 0)
    hi <- d$center + d$rmax * d$step
    Lc <- (input$p1 * lo[1] + input$p2 * lo[2]) / 1000
    Uc <- (input$p1 * hi[1] + input$p2 * hi[2]) / 1000

    d_i <- function(i, y) {
      s <- spec[[i]]
      switch(s$goal,
             Maximize = d_max(y, s$L, s$U),
             Minimize = d_min(y, s$L, s$U),
             Target   = d_target(y, s$L, s$T, s$U))
    }
    D_at <- function(x) {
      nd <- ref_newd(d, x)
      dv <- vapply(seq_len(nR), function(i)
        as.numeric(d_i(i, predict(mods[[i]]$fit, nd))), numeric(1))
      w <- vapply(spec, function(s) s$w, numeric(1))
      if (use_cost) { dv <- c(dv, d_min(cost_at(x), Lc, Uc)); w <- c(w, input$wc) }
      prod(dv^w)^(1 / sum(w))
    }
    obj <- function(x) {
      v <- D_at(x)
      if (!is.finite(v)) v <- 0
      -v + 10 * max(0, sqrt(sum(x^2)) - d$rmax)
    }

    set.seed(7)
    n_start <- 30
    starts  <- matrix(runif(n_start * k, -1, 1), n_start, k)
    runs <- vector("list", n_start)
    withProgress(message = "Optimizing...", value = 0, {
      for (i in seq_len(n_start)) {
        runs[[i]] <- optim(starts[i, ], obj, method = "Nelder-Mead",
                           control = list(maxit = 500))
        incProgress(1 / n_start)
      }
    })
    vals <- sapply(runs, function(r) -r$value)
    best <- runs[[which.max(vals)]]
    x    <- best$par
    nd   <- ref_newd(d, x)

    preds <- do.call(rbind, lapply(seq_len(nR), function(i) {
      p <- predict(mods[[i]]$fit, nd, interval = "prediction")
      data.frame(Response = d$resp[i], Predicted = round(p[1, 1], 3),
                 PI_low = round(p[1, 2], 3), PI_high = round(p[1, 3], 3))
    }))

    opt(list(
      setpoint = data.frame(Factor = d$factors,
                            Value = round(to_natural(x, as.numeric(d$center),
                                                     as.numeric(d$step)), 3)),
      preds = preds,
      cost  = if (use_cost) cost_at(x) else NULL,
      D = D_at(x), radius = sqrt(sum(x^2)), rmax = d$rmax,
      n_agree = sum(vals > max(vals) - 1e-3), n_start = n_start
    ))
  })

  output$opt_set  <- renderTable({ req(opt()); opt()$setpoint })
  output$opt_pred <- renderTable({ req(opt()); opt()$preds })
  output$opt_txt  <- renderText({
    o <- req(opt())
    c1 <- if (!is.null(o$cost)) sprintf("Media cost: %.4f EUR/L\n", o$cost) else ""
    warn <- if (o$D < 1e-6)
      "\nNo point has every desirability above 0: widen the Low/High limits.\n" else ""
    sprintf(paste0("%sComposite D = %.3f | radius %.2f (region %.2f)\n",
                   "%d of %d random starts reached the same optimum.\n",
                   "Intervals are 95%% prediction intervals for one new run in the ",
                   "reference block.%s"),
            c1, o$D, o$radius, o$rmax, o$n_agree, o$n_start, warn)
  })

  # ---- Design tab -------------------------------------------------
  gen <- eventReactive(input$gen, {
    k  <- input$d_k
    nm <- trimws(strsplit(input$d_names, ",")[[1]])
    ce <- as.numeric(strsplit(input$d_center, ",")[[1]])
    st <- as.numeric(strsplit(input$d_step, ",")[[1]])
    validate(need(length(nm) == k && length(ce) == k && length(st) == k,
                  "Names, centers and steps must each have exactly k entries."))
    validate(need(!anyNA(ce) && !anyNA(st), "Centers and steps must be numbers."))
    validate(need(input$d_type == "CCD" || k >= 3, "Box-Behnken needs at least 3 factors."))
    coding <- lapply(seq_len(k), function(i)
      as.formula(sprintf("x%d ~ (%s - (%s)) / %s", i, nm[i], ce[i], st[i])))
    set.seed(input$d_seed)
    des <- tryCatch(
      if (input$d_type == "CCD")
        ccd(basis = k, n0 = c(input$d_n0, input$d_n0), alpha = "rotatable",
            randomize = TRUE, coding = coding)
      else
        bbd(k = k, n0 = input$d_n0, randomize = TRUE, coding = coding),
      error = function(e) {
        showNotification(paste("Design error:", conditionMessage(e)), type = "error")
        NULL
      })
    req(des)
    list(nat = as.data.frame(decode.data(des)), coded = as.data.frame(des))
  })

  output$d_info <- renderText({
    g <- gen()
    sprintf("%d runs | blocks: %s", nrow(g$nat),
            paste(names(table(g$nat$Block)), table(g$nat$Block), sep = "=",
                  collapse = ", "))
  })
  output$d_tbl <- renderDT({
    datatable(gen()$nat, options = list(pageLength = 15, scrollX = TRUE))
  })
  output$dl <- downloadHandler(
    filename = function() paste0("runsheet_", tolower(input$d_type), ".csv"),
    content  = function(file) write.csv(gen()$nat, file, row.names = FALSE)
  )
}

shinyApp(ui, server)
