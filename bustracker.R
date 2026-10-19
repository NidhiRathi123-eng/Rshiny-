

library(shiny)

# ---------------- USER INTERFACE ----------------

ui <- fluidPage(
  titlePanel("City Bus Arrival Delay Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      h4("1. Bus Delay Data"),
      helpText("Enter delays in minutes, separated by commas."),
      textAreaInput(
        "delay", "All Bus Delays",
        value = "2,3,4,5,6,3,4,2,7,4"
      ),
      
      numericInput(
        "sigma", "Known Population SD (minutes)",
        value = 2, min = 0.01
      ),
      
      numericInput(
        "reference", "Reference Delay (minutes)",
        value = 5, min = 0
      ),
      
      h4("2. Late Buses Per Hour"),
      helpText("Enter late-bus counts for each hour."),
      textAreaInput(
        "late", "Hourly Late-Bus Counts",
        value = "2,3,1,4,2,3,2,1,3,4"
      ),
      
      numericInput(
        "k", "Number of late buses (k)",
        value = 3, min = 0, step = 1
      ),
      
      h4("3. Peak and Off-Peak Delays"),
      textAreaInput(
        "peak", "Peak-Hour Delays",
        value = "6,8,7,5,9,6,8"
      ),
      
      textAreaInput(
        "offpeak", "Off-Peak Delays",
        value = "2,4,3,5,2,3,4"
      ),
      
      actionButton("analyze", "Analyze Data")
    ),
    
    mainPanel(
      h3("Analysis Results"),
      
      h4("P4: One-Sample Z-Test"),
      verbatimTextOutput("zresult"),
      
      h4("P1: Poisson Distribution"),
      verbatimTextOutput("poissonresult"),
      
      h4("P3: 95% Confidence Interval"),
      verbatimTextOutput("ciresult"),
      
      h4("P8: Sign Test"),
      verbatimTextOutput("signresult"),
      
      h4("P9: Mann-Whitney U Test"),
      verbatimTextOutput("mwresult")
    )
  )
)

# ---------------- SERVER ----------------

server <- function(input, output, session) {
  
  # Convert comma-separated input into numeric data
  parse_data <- function(x) {
    x <- gsub("[;\n\r\t ]+", ",", trimws(x))
    values <- unlist(strsplit(x, ",", fixed = TRUE))
    values <- values[nzchar(values)]
    
    if (length(values) == 0) return(NULL)
    
    numbers <- suppressWarnings(as.numeric(values))
    
    if (anyNA(numbers) || any(!is.finite(numbers))) {
      return(NULL)
    }
    
    numbers
  }
  
  # Validate all data when Analyze is clicked
  data <- eventReactive(input$analyze, {
    
    d <- parse_data(input$delay)
    l <- parse_data(input$late)
    p <- parse_data(input$peak)
    o <- parse_data(input$offpeak)
    
    validate(
      need(!is.null(d) && length(d) >= 2 &&
             all(d >= 0),
           "Enter at least two valid non-negative bus delays."),
      need(!is.null(l) && all(l >= 0) &&
             all(l == floor(l)),
           "Hourly late-bus counts must be non-negative integers."),
      need(!is.null(p) && length(p) >= 1 &&
             all(p >= 0),
           "Enter valid peak-hour delays."),
      need(!is.null(o) && length(o) >= 1 &&
             all(o >= 0),
           "Enter valid off-peak delays."),
      need(is.finite(input$sigma) && input$sigma > 0,
           "Population SD must be greater than zero."),
      need(is.finite(input$reference) &&
             input$reference >= 0,
           "Reference delay must be non-negative."),
      need(is.finite(input$k) && input$k >= 0 &&
             input$k == floor(input$k),
           "k must be a non-negative integer.")
    )
    
    list(delay = d, late = l, peak = p, offpeak = o)
  })
  
  # P4: One-sample left-tailed Z-test
  output$zresult <- renderPrint({
    x <- data()$delay
    n <- length(x)
    mu <- mean(x)
    z <- (mu - 5) / (input$sigma / sqrt(n))
    pv <- pnorm(z)
    
    cat("H0: Mean delay >= 5 minutes\n")
    cat("H1: Mean delay < 5 minutes\n")
    cat("Sample size:", n, "\n")
    cat("Sample mean:", round(mu, 3), "minutes\n")
    cat("Z statistic:", round(z, 4), "\n")
    cat("P-value:", round(pv, 4), "\n")
    
    if (pv < 0.05) {
      cat("Conclusion: Reject H0. Evidence supports a mean delay below 5 minutes.")
    } else {
      cat("Conclusion: Fail to reject H0. Insufficient evidence that mean delay is below 5 minutes.")
    }
  })
  
  # P1: Poisson distribution
  output$poissonresult <- renderPrint({
    x <- data()$late
    lambda <- mean(x)
    k <- input$k
    
    cat("Estimated lambda:", round(lambda, 4), "\n")
    cat("Average late buses per hour:", round(lambda, 4), "\n")
    cat("P(X = 0):", round(dpois(0, lambda), 4), "\n")
    cat("P(X =", k, "):",
        round(dpois(k, lambda), 4), "\n")
    cat("P(X <=", k, "):",
        round(ppois(k, lambda), 4), "\n")
  })
  
  # P3: 95% confidence interval for population mean
  output$ciresult <- renderPrint({
    x <- data()$delay
    
    validate(
      need(length(x) >= 2 && sd(x) > 0,
           "At least two observations with some variation are needed.")
    )
    
    result <- t.test(x, conf.level = 0.95)
    
    cat("Sample mean:", round(mean(x), 4), "minutes\n")
    cat("95% Confidence Interval:\n")
    cat("Lower limit:",
        round(result$conf.int[1], 4), "minutes\n")
    cat("Upper limit:",
        round(result$conf.int[2], 4), "minutes\n")
  })
  
  # P8: Exact sign test for median delay
  output$signresult <- renderPrint({
    x <- data()$delay
    ref <- input$reference
    x <- x[x != ref]
    
    cat("H0: Population median =", ref, "minutes\n")
    cat("H1: Population median differs from", ref, "minutes\n")
    cat("Ties excluded:", sum(data()$delay == ref), "\n")
    
    if (length(x) == 0) {
      cat("Test cannot be performed: all observations are tied.")
    } else {
      positive <- sum(x > ref)
      result <- binom.test(
        positive, length(x), p = 0.5,
        alternative = "two.sided"
      )
      
      cat("Positive signs:", positive, "\n")
      cat("Negative signs:", sum(x < ref), "\n")
      cat("P-value:", round(result$p.value, 4), "\n")
      
      if (result$p.value < 0.05) {
        cat("Conclusion: Median delay differs significantly from the reference.")
      } else {
        cat("Conclusion: No significant evidence that the median differs from the reference.")
      }
    }
  })
  
  # P9: Mann-Whitney U test
  output$mwresult <- renderPrint({
    p <- data()$peak
    o <- data()$offpeak
    
    result <- wilcox.test(
      p, o,
      alternative = "two.sided",
      paired = FALSE,
      exact = FALSE
    )
    
    cat("H0: Peak and off-peak delay distributions are equal.\n")
    cat("H1: Peak and off-peak delay distributions differ.\n")
    cat("Peak mean delay:", round(mean(p), 3), "\n")
    cat("Off-peak mean delay:", round(mean(o), 3), "\n")
    cat("W statistic:", unname(result$statistic), "\n")
    cat("P-value:", round(result$p.value, 4), "\n")
    
    if (result$p.value < 0.05) {
      cat("Conclusion: Significant difference between peak and off-peak delays.")
    } else {
      cat("Conclusion: No statistically significant difference detected.")
    }
  })
}

# ---------------- RUN APPLICATION ----------------

shinyApp(ui = ui, server = server)
