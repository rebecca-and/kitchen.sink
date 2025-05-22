inla_summary <- function(inlamod_summary, model_description, hyper) {
  df <- data.frame(inlamod_summary$fixed) %>%
    select(mean = mean, lower = X0.025quant, upper = X0.975quant)

  if (hyper) {
    df %>%
    bind_rows(data.frame(inlamod_summary$hyperpar) %>%
                select(mean = mean, lower = X0.025quant, upper = X0.975quant)) %>%
    mutate(model = model_description) %>%
    rownames_to_column()
  } else {
    df %>%
      mutate(model = model_description) %>%
      rownames_to_column()
  }
}
