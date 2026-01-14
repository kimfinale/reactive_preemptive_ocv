# R/dashboard_app.R
library(shiny)
library(dplyr)
library(ggplot2)
library(tidyr)

# Source all utility functions in the same directory
# specific exclusion of this file to prevent recursion structure
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
if (length(current_dir) == 0 || current_dir == ".") current_dir <- "R"
if (!dir.exists(current_dir)) current_dir <- "."

r_files <- list.files(current_dir, pattern = "\\.R$", full.names = TRUE)
r_files <- r_files[!grepl("dashboard_app.R", r_files)]
r_files <- r_files[!grepl("check_setup.R", r_files)]

for (f in r_files) {
    tryCatch(
        {
            source(f, local = TRUE)
        },
        error = function(e) {
            message("Error sourcing ", f, ": ", e$message)
        }
    )
}

# Define UI
ui <- navbarPage(
    "Vaccination Strategy Dashboard",

    # Tab 1: Analytical Model
    tabPanel(
        "Analytical Model",
        sidebarLayout(
            sidebarPanel(
                h4("Parameters"),
                sliderInput("R_param", "Reproductive Number (R)",
                    min = 0.5, max = 5.0, value = 2.0, step = 0.1
                ),
                helpText("Normalized cost of outbreak"),
                sliderInput("r_param", "Reactive Effectiveness (r)",
                    min = 0.1, max = 1.0, value = 0.6, step = 0.05
                ),
                sliderInput("nu_param", "Pre-emptive Effectiveness (nu)",
                    min = 0.1, max = 1.0, value = 0.8, step = 0.05
                )
            ),
            mainPanel(
                tabsetPanel(
                    tabPanel(
                        "Single Population",
                        br(),
                        h4("Cost Comparison"),
                        plotOutput("costPlot_one"),
                        br(),
                        h4("Threshold Probability"),
                        tableOutput("threshTable_one")
                    ),
                    tabPanel(
                        "Multi-Population",
                        br(),
                        p("Explore mixed strategies with high-risk (f) populations."),
                        sidebarLayout(
                            sidebarPanel(
                                sliderInput("f_param", "Fraction High Risk (f)",
                                    min = 0.1, max = 0.9, value = 0.2, step = 0.05
                                ),
                                sliderInput("alpha_param", "Vaccination Allocation (alpha)",
                                    min = 0, max = 1, value = 0.5, step = 0.05
                                )
                            ),
                            mainPanel(
                                plotOutput("costPlot_multi")
                            )
                        )
                    )
                )
            )
        )
    ),

    # Tab 2: Simulation
    tabPanel(
        "Simulation",
        sidebarLayout(
            sidebarPanel(
                actionButton("runSim", "Run Simulation", class = "btn-primary"),
                hr(),
                numericInput("n_regions", "Number of Regions", value = 200, min = 50, max = 2000),
                sliderInput("outbreak_prob_mean", "Mean Outbreak Prob", min = 0, max = 1, value = 0.3),
                hr(),
                h4("Strategy Parameters"),
                sliderInput("cov_param", "Coverage", min = 0, max = 1, value = 0.7),
                h5("Reactive"),
                numericInput("max_targets", "Max RV Targets", value = 10),
                h5("Pre-emptive"),
                numericInput("B_pre_pct", "Budget (% of Pop)", value = 10, min = 0, max = 100)
            ),
            mainPanel(
                h4("Impact Comparison"),
                p("Comparison of strategies: Reactive Only, Pre-emptive Only, and Combined (Pre-emptive + Reactive on remainder)."),
                plotOutput("simPlot"),
                h4("Summary Statistics"),
                tableOutput("simSummary")
            )
        )
    )
)

# Define Server
server <- function(input, output, session) {
    # --- Single Population Analytical Logic ---

    output$costPlot_one <- renderPlot({
        p_seq <- seq(0, 1, length.out = 100)

        c_pre <- cost_pre_one(p = p_seq, R = input$R_param, r = input$r_param, nu = input$nu_param)
        c_react <- cost_react_one(p = p_seq, R = input$R_param, r = input$r_param)

        df <- data.frame(p = p_seq, Preemptive = c_pre, Reactive = c_react) %>%
            pivot_longer(cols = c("Preemptive", "Reactive"), names_to = "Strategy", values_to = "Cost")

        # Calculate intersection
        p_star <- p_star_one(R = input$R_param, r = input$r_param, nu = input$nu_param)
        p_star <- ifelse(p_star < 0 | p_star > 1, NA, p_star)

        g <- ggplot(df, aes(x = p, y = Cost, color = Strategy)) +
            geom_line(linewidth = 1.2) +
            theme_minimal(base_size = 14) +
            scale_color_brewer(palette = "Set1") +
            labs(x = "Probability of Outbreak (p)", y = "Normalized Cost")

        if (!is.na(p_star)) {
            c_star <- cost_pre_one(p_star, input$R_param, input$r_param, input$nu_param)
            g <- g +
                geom_vline(xintercept = p_star, linetype = "dashed") +
                annotate("text", x = p_star, y = max(df$Cost), label = paste0("p*=", round(p_star, 2)), vjust = -1)
        }
        g
    })

    output$threshTable_one <- renderTable({
        p_star <- p_star_one(R = input$R_param, r = input$r_param, nu = input$nu_param)
        interp <- if (p_star < 0) "Pre-emptive always worse" else if (p_star > 1) "Pre-emptive always better" else "Switch point exists"

        data.frame(
            Metric = c("Threshold p*", "Interpretation"),
            Value = c(round(p_star, 4), interp)
        )
    })

    # --- Multi Population Analytical Logic ---

    output$costPlot_multi <- renderPlot({
        p_seq <- seq(0, 1, length.out = 100)

        c_pre_m <- cost_pre_multi(p_seq, input$R_param, input$r_param, input$f_param, input$nu_param)
        c_react_m <- cost_react_multi(p_seq, input$R_param, input$r_param, input$f_param)

        df <- data.frame(p = p_seq, Preemptive = c_pre_m, Reactive = c_react_m) %>%
            pivot_longer(cols = c("Preemptive", "Reactive"), names_to = "Strategy", values_to = "Cost")

        ggplot(df, aes(x = p, y = Cost, color = Strategy)) +
            geom_line(linewidth = 1.2) +
            theme_minimal(base_size = 14) +
            labs(title = paste0("Multi-Pop Cost (f=", input$f_param, ")"))
    })

    # --- Simulation Logic ---

    sim_data <- eventReactive(input$runSim, {
        # 1. Setup Regions
        n <- input$n_regions
        # Outbreak prob with some spread
        prob <- pmin(1, pmax(0, rnorm(n, input$outbreak_prob_mean, 0.15)))

        regions <- tibble(
            id = 1:n,
            pop_size = sample(1000:10000, n, replace = TRUE),
            outbreak_prob = prob,
            Y = runif(n) < prob,
            # For pre-emptive score, assume it correlates with true prob
            S_hat = prob
        )

        # 2. Trajectories (truth)
        outbreaks <- assign_outbreak_trajectory(regions)

        # 3. Pre-emptive Allocation
        total_pop <- sum(regions$pop_size)
        B_pre <- round(total_pop * (input$B_pre_pct / 100))
        patch_df <- regions %>% select(id, m = pop_size, S_hat)
        alloc <- allocate_preemptive_patches(patch_df, B_pre = B_pre, patch_id = regions$id)

        # Add pre-vax status to outbreaks
        # Ensure order is preserved or use join
        outbreaks <- outbreaks %>%
            left_join(tibble(id = regions$id, pre_vax = alloc$pre_vax), by = "id")

        # 4. Impact Calculations

        # A) Pre-emptive Only
        # Reduction = nu * coverage IF pre_vax
        outbreaks <- outbreaks %>%
            mutate(
                cases_no_vax = attack_rate_per_capita * pop_size,

                # Pre-emptive Logic
                doses_pre = ifelse(pre_vax, pop_size * input$cov_param, 0),
                red_factor_pre = ifelse(pre_vax, input$cov_param * input$nu_param, 0),
                cases_averted_pre = cases_no_vax * red_factor_pre
            )

        # B) Reactive Only
        rv_only <- assign_reactive_vaccination(
            outbreaks,
            t_now = 50, max_targets = input$max_targets, coverage = input$cov_param
        )
        rv_only <- compute_frac_reduction_ci(rv_only, nu = input$nu_param)

        # Map back to main df
        # Note: assign_reactive_vaccination keeps all rows
        outbreaks$cases_averted_rv <- outbreaks$cases_no_vax * rv_only$frac_reduction_ci
        outbreaks$doses_rv <- ifelse(rv_only$reactive_vax, outbreaks$pop_size * input$cov_param, 0)

        # C) Combined (Pre-emptive + Reactive)
        # Reactive targets are chosen from those NOT already pre-vaccinated (or we assume boost? let's assume ignore pre-vaxed)
        # We need to filter potential reactive targets to !pre_vax

        # We can hack this by setting eligible=FALSE for pre-vaxed in assign_reactive_vaccination
        # But assign_reactive_vaccination calculates eligibility internally based on t_now.
        # We can modify 'ongoing' input? No.
        # We can run it, then just zero out reactive_vax if pre_vax is TRUE?
        # But that messes up "max_targets" allocation (we might waste slots on pre-vaxed).

        # Better: Filter input to assign_reactive?
        # assign_reactive checks 'ongoing'. If we dummy out 'ongoing' for pre-vax?

        outbreaks_for_combined <- outbreaks
        # If pre-vax, it's not "ongoing" in the sense of needing RV (simplified assumption)
        outbreaks_for_combined$t0[outbreaks_for_combined$pre_vax] <- NA # Remove from outbreak list effectively?
        # assign_reactive checks !is.na(t0) for has_outbreak.

        rv_comb <- assign_reactive_vaccination(
            outbreaks_for_combined,
            t_now = 50, max_targets = input$max_targets, coverage = input$cov_param
        )
        rv_comb_red <- compute_frac_reduction_ci(rv_comb, nu = input$nu_param)

        # If pre-vax, we get pre-emptive reduction.
        # If reactive (and not pre-vax), we get reactive reduction.
        # (rv_comb assumes 0 reduction if not reactive_vax)

        outbreaks$doses_comb_rv <- ifelse(rv_comb$reactive_vax, outbreaks$pop_size * input$cov_param, 0)
        outbreaks$cases_averted_comb_rv <- outbreaks$cases_no_vax * rv_comb_red$frac_reduction_ci

        outbreaks <- outbreaks %>%
            mutate(
                doses_comb = doses_pre + doses_comb_rv,
                # Total averted = Pre-emptive Averted + Reactive Averted (mutually exclusive in this logic)
                cases_averted_comb = cases_averted_pre + cases_averted_comb_rv
            )

        outbreaks
    })

    output$simSummary <- renderTable({
        req(sim_data())
        df <- sim_data()

        stats <- tibble(
            Strategy = c("Reactive Only", "Pre-emptive Only", "Combined"),
            Doses_Used = c(
                sum(df$doses_rv),
                sum(df$doses_pre),
                sum(df$doses_comb)
            ),
            Cases_Averted = c(
                sum(df$cases_averted_rv, na.rm = TRUE),
                sum(df$cases_averted_pre, na.rm = TRUE),
                sum(df$cases_averted_comb, na.rm = TRUE)
            )
        ) %>%
            mutate(
                Cases_per_1000_Doses = (Cases_Averted / Doses_Used) * 1000
            )

        stats
    })

    output$simPlot <- renderPlot({
        req(sim_data())
        df <- sim_data()

        plot_df <- df %>%
            select(id, cases_no_vax, cases_averted_rv, cases_averted_pre, cases_averted_comb) %>%
            pivot_longer(
                cols = starts_with("cases_averted"),
                names_to = "Strategy", values_to = "CasesAverted"
            ) %>%
            mutate(Strategy = recode(Strategy,
                "cases_averted_rv" = "Reactive",
                "cases_averted_pre" = "Pre-emptive",
                "cases_averted_comb" = "Combined"
            ))

        ggplot(plot_df, aes(x = Strategy, y = CasesAverted, fill = Strategy)) +
            geom_boxplot() +
            theme_minimal() +
            labs(title = "Cases Averted per Region by Strategy", y = "Cases Averted") +
            scale_fill_brewer(palette = "Pastel1")
    })
}

# Run the application
shinyApp(ui = ui, server = server)
