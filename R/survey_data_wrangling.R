single_year_to_five_year <- function (df, fifteen_to_49 = TRUE)  {
  df <- df %>% dplyr::mutate(age_group_label = cut(age, c(0,
                                                          seq(5, 85, 5) - 1), c(paste0(seq(0, 79, 5), "-", seq(5,
                                                                                                               80, 5) - 1), "80+"), include.lowest = TRUE)) %>% dplyr::left_join(naomi::get_age_groups() %>%
                                                                                                                                                                                   select(age_group, age_group_label)) %>% dplyr::select(-age_group_label)
  if (fifteen_to_49) {
    df %>% dplyr::filter(age %in% 15:49)
  }
  else {
    df
  }
}
