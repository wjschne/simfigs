# Primary parameters
list(
  tar_target(myfont, "Roboto Condensed"),
  tar_target(threshold_s, 85),
  tar_target(threshold_a, 85),
  tar_target(threshold_g, 90),
  tar_target(buffer, 5),
  tar_target(meaningful_difference, 10),
  tar_target(n, 100000),
  tar_target(v_true, c("g", "s", "a")),
  tar_target(v_observed, c("G", "S", "A")),
  tar_target(v_names, c(v_true, v_observed)),
  tar_target(
    decision_labels,
    c("SLD-Unlikely", "SLD-Possible", "SLD-Likely")
  ),
  tar_target(m, make_model()),
  tar_target(m_cov, lavaan::fitted(lavaan::sem(m))$cov[v_names, v_names]),
  tar_target(m_cor, cov2cor(m_cov)),
  tar_target(
    m_cor_table,
    knitr::kable(m_cor, digits = 2, caption = "Model-Implied Correlations")
  ),
  tar_target(
    d,
    make_data(
      m = m,
      n = n,
      threshold_s = threshold_s,
      threshold_a = threshold_a,
      threshold_g = threshold_g,
      buffer = buffer,
      meaningful_difference = meaningful_difference,
      decision_labels = decision_labels
    )
  ),
  tar_target(n_counts, deframe(count(d, outcome))),
  tar_target(accuracy_stats, diagnostic_accuracy(n_counts)),
  tar_target(maxdistplot, make_max_dist_plot(d, family = myfont)),
  tar_target(
    criteriaplot,
    make_criteria_plot(
      d,
      threshold_g = threshold_g,
      threshold_s = threshold_s,
      threshold_a = threshold_a,
      buffer = buffer,
      meaningful_difference = meaningful_difference,
      m_cor,
      family = myfont,
      decision_labels = decision_labels
    )
  ),
  tar_target(ppvnpv, make_ppvnpv_plot()),
  tar_target(
    certain_case,
    conditional_ppv(
      General = 100,
      Specific = 74,
      Academic = 77,
      threshold_g = 90,
      threshold_s = 85,
      threshold_a = 85,
      r_gg = 0.97,
      r_ss = 0.92,
      r_aa = 0.92,
      buffer = 5,
      meaningful_difference = 10,
      b_a.g = .5,
      b_a.s = .3,
      b_s.g = .65,
      filename = "certain_case",
      d = d
    )
  ),
  tar_target(
    uncertain_case,
    conditional_ppv(
      General = 97,
      Specific = 83,
      Academic = 84,
      threshold_g = 90,
      threshold_s = 85,
      threshold_a = 85,
      r_gg = 0.97,
      r_ss = 0.92,
      r_aa = 0.92,
      buffer = 5,
      meaningful_difference = 10,
      b_a.g = .5,
      b_a.s = .3,
      b_s.g = .65,
      filename = "uncertain_case",
      d = d
    )
  ),
  tar_target(
    uncertain_case_follow_up,
    conditional_ppv(
      General = 97,
      Specific = 83,
      Academic = 84,
      threshold_g = 90,
      threshold_s = 85,
      threshold_a = 85,
      r_gg = 0.97,
      r_ss = 0.96,
      r_aa = 0.96,
      buffer = 5,
      meaningful_difference = 10,
      b_a.g = .5,
      b_a.s = .3,
      b_s.g = .65,
      filename = "uncertain_case_follow_up",
      d = d
    )
  ),
  tar_target(
    mild_case,
    conditional_ppv(
      General = 99,
      Specific = 87,
      Academic = 86,
      threshold_g = 90,
      threshold_s = 85,
      threshold_a = 85,
      r_gg = 0.97,
      r_ss = 0.92,
      r_aa = 0.92,
      buffer = 5,
      meaningful_difference = 10,
      b_a.g = .5,
      b_a.s = .3,
      b_s.g = .65,
      filename = "mild_case",
      d = d
    )
  ),
  tar_target(
    mild_case_follow_up,
    conditional_ppv(
      General = 99,
      Specific = 87,
      Academic = 86,
      threshold_g = 90,
      threshold_s = 85,
      threshold_a = 85,
      r_gg = 0.97,
      r_ss = 0.96,
      r_aa = 0.96,
      buffer = 5,
      meaningful_difference = 10,
      b_a.g = .5,
      b_a.s = .3,
      b_s.g = .65,
      filename = "mild_case_follow_up",
      d = d
    )
  )
)
